import os, glob
import sys
import argparse
import importlib
import json
import shutil
import tarfile
import numpy as np
import configparser
import logging
import atcf_tools as atools
import trop_cyclone as tc
import fcst_metrics_tc as fmtc
import compute_tc_fields as tcf
import nhc_sens as sens

#  Routine to read configuration file
def read_config(datea, storm, filename):
    '''
    This function reads a configuration file, and puts all of the appropriate variables into 
    a nested dictionary that can be passed around the appropriate scripts.  The result is the 
    configuration dictionary.

    Attributes:
        datea  (string):  The initialization time of the forecast (yyyymmddhh)
        storm  (string):  TC name, where XXXXXXXNNB, where XXXXXXXX is the name, NN is the number, B is the basin
        filename (dict):  The configuration file with all of the parameters
    '''

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

    #  Modify work and output directory for specific case/time
    config['work_dir']   = '{0}/{1}.{2}'.format(config['work_dir'],datea,storm)
    config['output_dir'] = '{0}/{1}.{2}'.format(config['output_dir'],datea,storm)
    config['storm']      = storm

    #  Create appropriate directories
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

    return(config)


def main():
    '''
    This is the main routine that calls all of the steps needed to compute ensemble-based
    sensitivity for TC forecasts.  The script can be called from the command line, where the
    user inputs the forecast initialization date, and storm name.  The user can also add the
    path to the parameter file.

    Important:  within the parameter file, the user needs to set the variable io_module, which
    contains information for how to read and use grib and ATCF data from a specific source and
    model.  The module specified in this variable will be used to get all input data.

    From command line:

    python run_NHC_sens.py -init yyyymmddhh --storm XXXXXXNNB --param paramfile

      where:

        -init is the initialization date in yyyymmddhh format
        -storm is the TC name (XXXXXX is the storm name, NN is the number, B is the basin)
        -param is the parameter file path (optional, otherwise goes to default values in default.parm)
    '''

    #  Read the initialization time and storm from the command line
    exp_parser = argparse.ArgumentParser()
    exp_parser.add_argument('--init',  action='store', type=str, required=True)
    exp_parser.add_argument('--storm', action='store', type=str, required=True)
    exp_parser.add_argument('--param', action='store', type=str)

    args = exp_parser.parse_args()

    datea = args.init
    storm = args.storm

    if args.param:
       paramfile = args.param
    else:
       paramfile = 'example.parm'

    #  Read the configuration file and set up for usage later
    config = read_config(datea, storm, paramfile)

    #  Import the module that contains routines to read ATCF and Grib data specific to the model
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
    logging.basicConfig(filename='{0}/{1}_{2}.log'.format(config.get('log_dir','.'),str(datea),storm), \
                               filemode='w', format='%(message)s')
    logging.warning("STARTING SENSITIVITIES for {0} on {1}".format(bbnnyyyy, str(datea)))


    #  Create directories to put output graphics 
    if ( config.get('webpage','True') == 'True' ):
       if not os.path.isdir('{0}/{1}.{2}'.format(config['html_dir'],datea,storm)):
          os.makedirs('{0}/{1}.{2}'.format(config['html_dir'],datea,storm))
       if not os.path.islink('{0}/{1}.{2}'.format(config['work_dir'],datea,storm)):
          os.symlink('{0}/{1}.{2}'.format(config['html_dir'],datea,storm), \
                     '{0}/{1}.{2}'.format(config['work_dir'],datea,storm))
    
    #  Copy grib and ATCF data to the work directory
    dpp.stage_grib_files(datea, config)
    dpp.stage_atcf_files(datea, bbnnyyyy, config)


    #  Read ATCF data into dictionary
    atcf = atools.ReadATCFData('{0}/atcf_*.dat'.format(config['work_dir']))
    btk = atcf.get_best_data(bbnnyyyy)


    #  Plot the ensemble forecast
    config['vitals_plot']['track_output_dir'] = config['vitals_plot'].get('track_output_dir', '{0}/{1}.{2}'.format(config['work_dir'],str(datea),storm))
    config['vitals_plot']['int_output_dir'] = config['vitals_plot'].get('int_output_dir', '{0}/{1}.{2}'.format(config['work_dir'],str(datea),storm))
    tc.atcf_ens_tc_vitals(atcf, btk, datea, storm, config['model_src'], 50, config['vitals_plot'], config['work_dir'])


    #  Compute TC-related forecast metrics
    fmtc.ComputeForecastMetrics(datea, storm, atcf, config)


    #  Compute forecast fields at each desired time to use in sensitivity calculation
    for fhr in range(0,int(config['fcst_hour_max'])+int(config['fcst_hour_int']),int(config['fcst_hour_int'])):

       tcf.ComputeTCFields(datea, fhr, atcf, config)


    #  Compute sensitivity of each metric to forecast fields at earlier times, as specified by the user
    metlist = [e.strip() for e in config['sens']['metrics'].split(',')]
    for i in range(len(metlist)):

       #  Limit loop over time to forecast metric lead time (i.e., for a 72 h forecast, do not compute 
       #  the sensitivity to fields beyond 72 h
       a = metlist[i].split('_')
       fhrstr = a[0]
       fhrmax = int(np.min([float(fhrstr[1:4]),float(config['fcst_hour_max'])]))

       for fhr in range(0,fhrmax+int(config['fcst_hour_int']),int(config['fcst_hour_int'])):

          sens.ComputeSensitivity(datea, fhr, metlist[i], atcf, config)


    #  Save some of the files, if needed
    if ( config.get('archive_metric','False') == 'True' ):
       print("Add capability")

    if ( config.get('archive_fields','False') == 'True' ):
       os.rename('{0}/\*_ens.nc'.format(config['work_dir']), '{0}/.'.format(config['output_dir']))


    #  Create a tar file of gridded sensitivity files, if needed
    os.chdir(config['work_dir'])

    tarout = '{0}/{1}.tar'.format(config['outgrid_dir'],datea) 
    if ( os.path.isfile(tarout) and tarfile.is_tarfile(tarout) ):
       tar = tarfile.open(tarout) 
       tar.extractall()
       tar.close()

    for f in ['usteer', 'vsteer', 'masteer', 'misteer']:
       for g in glob.glob('{0}/{1}.{2}/*_intmajtrack/sens/{3}/*.nc'.format(config['work_dir'],datea,storm,f)): 
          shutil.copy(g, '{0}/{1}{2}/.'.format(datea,bbl,storm[-3:-1]))

    tar = tarfile.open(tarout, 'w')
    for f in glob.glob('{0}/*/*.nc'.format(datea)):
       tar.add(f)
    tar.close()


    #  Clean up work directory, if desired
#    os.chdir(config['work_dir'] + "../")
#    if ( config.get('save_work_dir','False') == 'False' ):
#      shutil.rmtree(config['work_dir'])

if __name__ == '__main__':

   main()
