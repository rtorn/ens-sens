import os, glob
import sys
import argparse
import importlib
import json
import shutil
import tarfile
import numpy as np
import datetime as dt
import configparser
import logging

import matplotlib
from IPython.core.pylabtools import figsize, getfigs
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs

import atcf_tools as atools
import trop_cyclone as tc
import fcst_metrics_tc as fmtc
import compute_tc_fields as tcf
import nhc_sens as sens
from SensPlotRoutines import background_map

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
    config['figure_dir'] = '{0}/{1}.{2}'.format(config['figure_dir'],datea,storm)
    config['storm']      = storm

    #  Create appropriate directories
    if not os.path.isdir(config['work_dir']):
      try:
        os.makedirs(config['work_dir'])
      except OSError as e:
        raise e

    if (eval(config.get('archive_metric','False')) or eval(config.get('archive_metric','False')) ) and \
               (not os.path.isdir(config['output_dir'])):
      try:
        os.makedirs(config['output_dir'])
      except OSError as e:
        raise e

    if not os.path.isdir(config['figure_dir']):
      try:
        os.makedirs(config['figure_dir'])
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

    for handler in logging.root.handlers[:]:
       logging.root.removeHandler(handler)
    logging.basicConfig(filename='{0}/{1}_{2}.log'.format(config.get('log_dir','.'),str(datea),storm), \
                               filemode='w', format='%(asctime)s;%(message)s')
    logging.warning("STARTING SENSITIVITIES for {0} on {1}".format(bbnnyyyy, str(datea)))


    #  Copy grib and ATCF data to the work directory
    logging.info("Staging Grib Files")
    dpp.stage_grib_files(datea, config)
    logging.info("Staging ATCF Files")
    dpp.stage_atcf_files(datea, bbnnyyyy, config)
#    dpp.stage_best_file(bbnnyyyy, config)


    #  Read ATCF data into dictionary
    logging.info("Reading ATCF Files")
    atcf = atools.ReadATCFData('{0}/atcf_*.dat'.format(config['work_dir']))
    atcf.read_best_data('{0}/b{1}.dat'.format(config['work_dir'],bbnnyyyy))


    #  Plot the ensemble forecast
    config['vitals_plot']['track_output_dir'] = config['vitals_plot'].get('track_output_dir', config['figure_dir'])
    config['vitals_plot']['int_output_dir'] = config['vitals_plot'].get('int_output_dir', config['figure_dir'])
    tc.plot_ens_tc_track(atcf, storm, datea, config) 
    tc.plot_ens_tc_intensity(atcf, storm, datea, config)

    #  Plot the precipitation forecast
#    fhr1 = json.loads(config['vitals_plot'].get('precip_hour_1'))
#    fhr2 = json.loads(config['vitals_plot'].get('precip_hour_2'))
    fhr1 = [48]
    fhr2 = [72]

    for h in range(len(fhr1)):
       precipitation_ens_maps(datea, int(fhr1[h]), int(fhr2[h]), config)


    #  Compute TC-related forecast metrics
    logging.info("Computing forecast Metrics")
    fmtc.ComputeForecastMetrics(datea, storm, atcf, config)


    #  Compute forecast fields at each desired time to use in sensitivity calculation
    for fhr in range(0,int(config['fcst_hour_max'])+int(config['fcst_hour_int']),int(config['fcst_hour_int'])):

       logging.debug(f"Computing Fields {fhr}")
       tcf.ComputeTCFields(datea, fhr, atcf, config)


    #  Compute sensitivity of each metric to forecast fields at earlier times, as specified by the user
    logging.info("Computing Sensitivity")
    metlist = [e.strip() for e in config['sens']['metrics'].split(',')]
    for i in range(len(metlist)):

       #  Limit loop over time to forecast metric lead time (i.e., for a 72 h forecast, do not compute 
       #  the sensitivity to fields beyond 72 h
       a = metlist[i].split('_')
       fhrstr = a[0]
       fhrmax = int(np.min([float(fhrstr[1:4]),float(config['fcst_hour_max'])]))

       for fhr in range(0,fhrmax+int(config['fcst_hour_int']),int(config['fcst_hour_int'])):

          sens.ComputeSensitivity(datea, fhr, metlist[i], atcf, config)


    with open('{0}/metric_list'.format(config['work_dir']), 'w') as f:
       for item in metlist:
          f.write("%s\n" % item)
    f.close()


    #  Save some of the files, if needed
    if ( config.get('archive_metric','False') == 'True' ):
       for met in metlist:
          os.rename('{0}/{1}_{2}.nc'.format(config['work_dir'],datea,met), '{0}/.'.format(config['output_dir']))

    if ( config.get('archive_fields','False') == 'True' ):
       os.rename('{0}/\*_ens.nc'.format(config['work_dir']), '{0}/.'.format(config['output_dir']))


    #  Create a tar file of gridded sensitivity files, if needed
    os.chdir(config['work_dir'])

    tarout = '{0}/{1}.tar'.format(config['outgrid_dir'],datea) 
    if ( os.path.isfile(tarout) and tarfile.is_tarfile(tarout) ):
       os.system('tar --skip-old-files -xf {0}'.format(tarout))

    tar = tarfile.open(tarout, 'w')
    for f in glob.glob('{0}/*/*.nc'.format(datea)):
       tar.add(f)
    tar.close()


    #  Clean up work directory, if desired
    os.chdir('{0}/..'.format(config['work_dir']))
    if not eval(config.get('save_work_dir','False')):
       shutil.rmtree(config['work_dir'])


def precipitation_ens_maps(datea, fhr1, fhr2, config):
    '''
    Function that plots the ensemble precipitation forecast between two forecast hours.

    Attributes:
        datea (string):  initialization date of the forecast (yyyymmddhh format)
        fhr1     (int):  starting forecast hour of the window
        fhr2     (int):  ending forecast hour of the window
        config (dict.):  dictionary that contains configuration options (read from file)
    '''

    dpp = importlib.import_module(config['io_module'])

    lat1 = float(config['vitals_plot'].get('min_lat_precip','22.'))
    lat2 = float(config['vitals_plot'].get('max_lat_precip','50.'))
    lon1 = float(config['vitals_plot'].get('min_lon_precip','-100.'))
    lon2 = float(config['vitals_plot'].get('max_lon_precip','-65.'))

    fff1 = '%0.3i' % fhr1
    fff2 = '%0.3i' % fhr2
    datea_1   = dt.datetime.strptime(datea, '%Y%m%d%H') + dt.timedelta(hours=fhr1)
    date1_str = datea_1.strftime("%Y%m%d%H")
    datea_2   = dt.datetime.strptime(datea, '%Y%m%d%H') + dt.timedelta(hours=fhr2)
    date2_str = datea_2.strftime("%Y%m%d%H")

    #  Read the total precipitation for the beginning of the window
    g1 = dpp.ReadGribFiles(datea, fhr1, config)

    vDict = {'latitude': (lat1, lat2), 'longitude': (lon1, lon2),
             'description': 'precipitation', 'units': 'mm', '_FillValue': -9999.}
    vDict = g1.set_var_bounds('precipitation', vDict)

    g2 = dpp.ReadGribFiles(datea, fhr2, config)

    ensmat = g2.create_ens_array('precipitation', g2.nens, vDict)

    for n in range(g2.nens):
       ens1 = np.squeeze(g1.read_grib_field('precipitation', n, vDict))
       ens2 = np.squeeze(g2.read_grib_field('precipitation', n, vDict))
       ensmat[n,:,:] = ens2[:,:] - ens1[:,:]

    if hasattr(ens2, 'units'):
       if ens2.units == "m":
          vscale = 1000.
       else:
          vscale = 1.
    else:
       vscale = 1.

    #  Scale all of the rainfall to mm and to a 24 h precipitation
    ensmat[:,:,:] = ensmat[:,:,:] * vscale * 24. / float(fhr2-fhr1)

    e_mean = np.mean(ensmat, axis=0)
    e_std  = np.std(ensmat, axis=0)

    #  Create basic figure, including political boundaries and grid lines
    fig = plt.figure(figsize=(11,6.5), constrained_layout=True)

    colorlist = ("#FFFFFF", "#00ECEC", "#01A0F6", "#0000F6", "#00FF00", "#00C800", "#009000", "#FFFF00", \
                 "#E7C000", "#FF9000", "#FF0000", "#D60000", "#C00000", "#FF00FF", "#9955C9")

    plotBase = config.copy()
    plotBase['subplot']       = 'True'
    plotBase['subrows']       = 1
    plotBase['subcols']       = 2
    plotBase['subnumber']     = 1
    plotBase['grid_interval'] = config['vitals_plot'].get('grid_interval', 5)
    plotBase['left_labels'] = 'True'
    plotBase['right_labels'] = 'None'
    ax1 = background_map(config.get('projection', 'PlateCarree'), lon1, lon2, lat1, lat2, plotBase)

    #  Plot the mean precipitation map
    mpcp = [0.0, 0.25, 0.50, 1., 1.5, 2., 4., 6., 8., 12., 16., 24., 32., 64., 96., 97.]
    norm = matplotlib.colors.BoundaryNorm(mpcp,len(mpcp))
    pltf1 = plt.contourf(ensmat.longitude.values,ensmat.latitude.values,e_mean,mpcp,norm=norm,extend='max', \
                         cmap=matplotlib.colors.ListedColormap(colorlist), transform=ccrs.PlateCarree())

    cbar = plt.colorbar(pltf1, fraction=0.15, aspect=45., pad=0.04, orientation='horizontal', ticks=mpcp)
    cbar.set_ticks(mpcp[1:(len(mpcp)-1):2])

    plt.title('Mean')

    plotBase['subnumber']     = 2
    plotBase['left_labels'] = 'None'
    plotBase['right_labels'] = 'None'
    ax2 = background_map(config.get('projection', 'PlateCarree'), lon1, lon2, lat1, lat2, plotBase)

    #  Plot the standard deviation of the ensemble precipitation
    spcp = [0., 3., 6., 9., 12., 15., 18., 21., 24., 27., 30., 33., 36., 39., 42., 43.]
    norm = matplotlib.colors.BoundaryNorm(spcp,len(spcp))
    pltf2 = plt.contourf(ensmat.longitude.values,ensmat.latitude.values,e_std,spcp,norm=norm,extend='max', \
                         cmap=matplotlib.colors.ListedColormap(colorlist), transform=ccrs.PlateCarree())

    cbar = plt.colorbar(pltf2, fraction=0.15, aspect=45., pad=0.04, orientation='horizontal', ticks=spcp)
    cbar.set_ticks(spcp[1:(len(spcp)-1)])

    plt.title('Standard Deviation')

    fig.suptitle('F{0}-F{1} Precipitation ({2}-{3})'.format(fff1, fff2, date1_str, date2_str), fontsize=16)

    outdir = '{0}/std/pcp'.format(config['figure_dir'])
    if not os.path.isdir(outdir):
       try:
          os.makedirs(outdir)
       except OSError as e:
          raise e

    plt.savefig('{0}/{1}_f{2}_pcp24h_std.png'.format(outdir,datea,fff2),format='png',dpi=120,bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':

   main()
