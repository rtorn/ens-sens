import os, glob
import sys, getopt
import importlib
import json
import time
import shutil
import tarfile
import numpy as np
import xarray as xr
import configparser
import datetime as dt
import atcf_tools as atools
import trop_cyclone as tc
import fcst_metrics_tc as fmtc
import compute_tc_fields as tcf
import nhc_sens as sens

#  Routine to read configuration file
def read_config(datea, storm, filename):

    confin = configparser.ConfigParser()
    confin.read(filename)

    config = {}
    config['vitals_plot'] = confin['vitals_plot']
    config['metric']      = confin['metric']
    config['fields']      = confin['fields']
    config['sens']        = confin['sens']
#    config['display']     = confin['display']
    config.update(confin['model'])
    config.update(confin['locations'])

    config['work_dir']   = config['work_dir']   + "/" + datea + "." + storm
    config['output_dir'] = config['output_dir'] + "/" + datea + "." + storm
    config['storm']      = storm

    if not os.path.isdir(config['work_dir']):
      try:
        os.makedirs(config['work_dir'])
      except OSError as e:
        raise e

    if not os.path.isdir(config['output_dir']):
      try:
        os.makedirs(config['output_dir'])
      except OSError as e:
        raise e

    the_dict = {}
    for section in confin.sections():
        the_dict[section] = {}
        for key, val in confin.items(section):
            the_dict[section][key] = val

    return(config)


def main():

    paramfile = 'example.parm'

    #  Read the initialization time and storm from command line
    try:
      opts, args = getopt.getopt(sys.argv[1:],"hu:d:s:p:",["infile="])
    except getopt.GetoptError:
      print('run_NHC_sens.py -d <yyyymmddhh> -s <storm_name>')
      sys.exit(2)
    for opt, arg in opts:
      if opt == '-h':
        print('run_NHC_sens.py -d <yyyymmddhh> -s <storm_name>')
        sys.exit()
      elif opt in ("-d", "--init_date"):
        datea = arg
      elif opt in ("-s", "--storm_name"):
        storm = arg
      elif opt in ("-p", "--param_file"):
        paramfile = arg

    config = read_config(datea, storm, paramfile)

    dpp = importlib.import_module(config['io_module'])

    os.chdir(config['work_dir'])

    #  Set the domain parameters based on basin
    if storm[-1] == "l":
      bbl = "al"
      config['sens']['min_lat'] = config['sens'].get('min_lat','8.0')
      config['sens']['max_lat'] = config['sens'].get('max_lat','65.0')
      config['sens']['min_lon'] = config['sens'].get('min_lon','-140.0')
      config['sens']['max_lon'] = config['sens'].get('max_lon','-20.0')
    elif storm[-1] == "e":
      bbl = "ep"
      config['sens']['min_lat'] = config['sens'].get('min_lat','8.0')
      config['sens']['max_lat'] = config['sens'].get('max_lat','65.0')
      config['sens']['min_lon'] = config['sens'].get('min_lon','-180.0')
      config['sens']['max_lon'] = config['sens'].get('max_lon','-80.0')
    elif storm[-1] == "w":
      bbl = "wp"

    bbnnyyyy = "{0}{1}{2}".format(bbl, storm[-3:-1], datea[0:4])

    if ( config.get('webpage','True') == 'True' ):
      if not os.path.isdir(config['html_dir'] + '/' + datea + "." + storm):
         os.makedirs(config['html_dir'] + '/' + datea + "." + storm)
      if not os.path.islink(config['work_dir'] + '/' + datea + "." + storm):
         os.symlink(config['html_dir'] + '/' + datea + '.' + storm, config['work_dir'] + '/' + datea + "." + storm)

    print("STARTING SENSITIVITIES for {0} on {1}".format(bbnnyyyy, str(datea)))

    #  Copy grib and ATCF data to the work directory
    dpp.stage_grib_files(datea, config)
    dpp.stage_atcf_files(datea, bbnnyyyy, config)

    #  Read ATCF data into dictionary
    atcf = atools.ReadATCFData(config['work_dir']+"/atcf_*.dat")
    btk = atcf.get_best_data(bbnnyyyy)

    #  Plot the ensemble forecast
    config['vitals_plot']['track_output_dir'] = config['vitals_plot'].get('track_output_dir', config['work_dir'] + "/" + str(datea) + '.' + storm)
    config['vitals_plot']['int_output_dir'] = config['vitals_plot'].get('int_output_dir', config['work_dir'] + "/" + str(datea) + '.' + storm)
    tc.atcf_ens_tc_vitals(atcf, btk, datea, storm, config['model_src'], 50, config['vitals_plot'], config['work_dir'])

    #  Compute forecast metrics
    fmtc.ComputeForecastMetrics(datea, atcf, config)

    #  Compute forecast fields to use in sensitivity calculation
    for fhr in range(0,int(config['fcst_hour_max'])+int(config['fcst_hour_int']),int(config['fcst_hour_int'])):

      tcf.ComputeTCFields(datea, fhr, atcf, config)

#    tcf.ComputeTCFields(datea, 48, atcf, config) 

    #  Compute sensitivity of each metric to forecast fields at earlier times
    metlist = [e.strip() for e in config['sens']['metrics'].split(',')]
    for i in range(len(metlist)):

      a = metlist[i].split('_')
      fhrstr = a[0]
      fhrmax = int(np.min([float(fhrstr[1:4]),float(config['fcst_hour_max'])]))

      for fhr in range(0,fhrmax+int(config['fcst_hour_int']),int(config['fcst_hour_int'])):

        sens.ComputeSensitivity(datea, fhr, metlist[i], atcf, config)
#        print(metlist[i],fhr)

#    sens.ComputeSensitivity(datea, 48, 'f120_intmajtrack', atcf, config)

    #  Save some of the files, if needed
    if ( config.get('archive_metric','False') == 'True' ):
      print("Add capability")

    if ( config.get('archive_fields','False') == 'True' ):
      os.rename(config['work_dir'] + '/\*_ens.nc', config['output_dir'] + '/.')

    #  Create a tar file of gridded sensitivity files, if needed
    os.chdir(config['work_dir'])

    tarout = config['outgrid_dir'] + '/' + datea + '.tar'
    if ( os.path.isfile(tarout) and tarfile.is_tarfile(tarout) ):
      tar = tarfile.open(tarout) 
      tar.extractall()
      tar.close()

    for f in ['usteer', 'vsteer', 'masteer', 'misteer']:
      for g in glob.glob(config['work_dir'] + "/" + datea + '.' + storm + '/*_intmajtrack/sens/' + f + '/*.nc'): 
        shutil.copy(g, datea + "/" + bbl + storm[-3:-1] + "/.")

    tar = tarfile.open(tarout, 'w')
    for f in glob.glob(datea + '/*/*.nc'):
      tar.add(f)
    tar.close()

    #  Clean up work directory
#    os.chdir(config['work_dir'] + "../")
#    if ( config.get('save_work_dir','False') == 'False' ):
#      shutil.rmtree(config['work_dir'])

if __name__ == '__main__':

   main()
