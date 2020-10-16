import matplotlib
from IPython.core.pylabtools import figsize, getfigs
import os

import sys
import netCDF4 as nc
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import numpy as np
from matplotlib import colors
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from math import radians, degrees, sin, cos, asin, acos, sqrt
from SensPlotRoutines import plotVecSens, plotScalarSens, computeSens, writeSensFile

####   Class for computing TC-related sensitivity fields for a single metric and forecast lead time
class ComputeSensitivity:

    def __init__(self, datea, fhr, metname, atcf, config):

      self.df_files = atcf.atcf_files
      self.atcf_array = atcf.atcf_array
      self.nens = int(len(atcf.atcf_files))

      plotDict = {}
      for key in config['sens']:
        plotDict[key] = config['sens'][key]

      fhrt = '%0.3i' % fhr

      print('Sensitivity of ' + metname + ' to F' + fhrt)

      if config['storm'][-1] == "l":
        bbnn = "al" + config['storm'][-3:-1]
      elif config['storm'][-1] == "e":
        bbnn = "ep" + config['storm'][-3:-1]
      elif storm[-1] == "w":
        bbnn = "wp" + config['storm'][-3:-1]

      elat, elon = atcf.ens_lat_lon_time(fhr)

      m_lat = 0.0
      m_lon = 0.0
      e_cnt = 0.0
      for n in range(self.nens):
        if elat[n] != 0.0 and elon[n] != 0.0:
          e_cnt = e_cnt + 1
          m_lat = m_lat + elat[n]
          m_lon = m_lon + elon[n]

      plotDict['tcLat']     = m_lat / e_cnt
      plotDict['tcLon']     = m_lon / e_cnt
      plotDict['plotTitle'] = datea + " F" + fhrt
      plotDict['fileTitle'] = "TEST JHT-Torn ECMWF Sensitivity"
      plotDict['initDate']  = datea[0:4] + "-" + datea[4:6] + "-" + datea[6:8] + " " + datea[8:10] + ":00:00"

      #  Obtain the metric information (here read from file)
      mfile = nc.Dataset(config['output_dir'] + '/' + datea + '_' + metname + '.nc')

      metric = mfile.variables['fore_met_init'][:]
      nens   = len(metric)
      metric = metric[:] - np.mean(metric, axis=0)

      plotDict['sensmax'] = np.std(metric) * 0.9

      #  Read major axis direction if appropriate  
      if hasattr(mfile.variables['fore_met_init'],'units'):
        plotDict['metricUnits'] = mfile.variables['fore_met_init'].units

      if hasattr(mfile,'X_DIRECTION_VECTOR') and hasattr(mfile,'Y_DIRECTION_VECTOR'):
        ivec = mfile.X_DIRECTION_VECTOR
        jvec = mfile.Y_DIRECTION_VECTOR
      else:
        ivec = 1.0
        jvec = 0.0

      stceDict = plotDict.copy()
      stceDict['output_sens']=False
      stceDict['range_rings']='True'
      stceDict['min_lat']=float(plotDict['tcLat'])-float(plotDict.get('storm_center_radius', 10.))
      stceDict['max_lat']=float(plotDict['tcLat'])+float(plotDict.get('storm_center_radius', 10.))
      stceDict['min_lon']=float(plotDict['tcLon'])-float(plotDict.get('storm_center_radius', 10.))
      stceDict['max_lon']=float(plotDict['tcLon'])+float(plotDict.get('storm_center_radius', 10.))
      stceDict['grid_interval']=3.
      stceDict['barb_interval']=3
      stceDict["figsize"]=(8.5,11)

      #  Read the ensemble zonal and meridional steering wind, compute ensemble mean
      ufile = nc.Dataset(config['work_dir'] + "/" + datea + "_f" + fhrt + "_usteer_ens.nc")
      uVar = ufile.variables

      lat  = uVar['latitude'][:]
      lon  = uVar['longitude'][:]
      uens = np.squeeze(uVar['ensemble_data'][:])
      umea = np.mean(uens, axis=0)
      uvar = np.var(uens, axis=0)

      vfile = nc.Dataset(config['work_dir'] + "/" + datea + "_f" + fhrt + "_vsteer_ens.nc")
      vVar = vfile.variables

      vens = np.squeeze(vVar['ensemble_data'][:])
      vmea = np.mean(vens, axis=0)
      vvar = np.var(vens, axis=0)


      #  Compute sensitivity with respect to zonal steering wind
      sens, sigv = computeSens(uens, uvar, metric)
      sens[:,:] = sens[:,:] * np.sqrt(uvar[:,:])

      outdir = config['work_dir'] + "/" + datea + '.' + config['storm'] + '/' + metname + '/sens/usteer'
      if not os.path.isdir(outdir):
        os.makedirs(outdir)

      if plotDict.get('output_sens', False):
        if not os.path.isdir(datea + '/' + bbnn):
          os.makedirs(datea + '/' + bbnn)
        writeSensFile(lat, lon, fhr, umea, sens, sigv, datea + '/' + bbnn + '/' + datea + "_f" + fhrt + "_usteer_sens.nc", plotDict)

      plotVecSens(lat, lon, sens, umea, vmea, sigv, outdir + '/' + datea + "_f" + fhrt + "_usteer_sens.png", plotDict)


      #  Compute sensitivity with respect to meridional steering wind
      sens, sigv = computeSens(vens, vvar, metric)
      sens[:,:] = sens[:,:] * np.sqrt(vvar[:,:])

      outdir = config['work_dir'] + "/" + datea + '.' + config['storm'] + '/' + metname + '/sens/vsteer'
      if not os.path.isdir(outdir):
        os.makedirs(outdir)

      if plotDict.get('output_sens', False):
        writeSensFile(lat, lon, fhr, umea, sens, sigv, datea + '/' + bbnn + '/' + datea + "_f" + fhrt + "_vsteer_sens.nc", plotDict)

      plotVecSens(lat, lon, sens, umea, vmea, sigv, outdir + '/' + datea + "_f" + fhrt + "_vsteer_sens.png", plotDict)


      #  Rotate wind into major axis direction, compute sensitivity to steering wind in that direction
      wens = ivec * uens[:,:,:] + jvec * vens[:,:,:]
      emea = np.mean(wens, axis=0)
      evar = np.var(wens, axis=0)
      sens, sigv = computeSens(wens, evar, metric)
      sens[:,:] = sens[:,:] * np.sqrt(evar[:,:])

      outdir = config['work_dir'] + "/" + datea + '.' + config['storm'] + '/' + metname + '/sens/masteer'
      if not os.path.isdir(outdir):
        os.makedirs(outdir)

      if plotDict.get('output_sens', False):
        writeSensFile(lat, lon, fhr, emea, sens, sigv, datea + '/' + bbnn + '/' + datea + "_f" + fhrt + "_masteer_sens.nc", plotDict)

      plotVecSens(lat, lon, sens, umea, vmea, sigv, outdir + '/' + datea + "_f" + fhrt + "_masteer_sens.png", plotDict)

      outdir = config['work_dir'] + "/" + datea + '.' + config['storm'] + '/' + metname + '/sens_sc/masteer'
      if not os.path.isdir(outdir):
        os.makedirs(outdir)

      plotVecSens(lat, lon, sens, umea, vmea, sigv, outdir + '/' + datea + "_f" + fhrt + "_masteer_sens.png", stceDict)


      #  Read 250 hPa PV, compute sensitivity to that field
      if os.path.isfile(config['work_dir'] + "/" + datea + "_f" + fhrt + "_pv250_ens.nc"): 
      
         efile = nc.Dataset(config['work_dir'] + "/" + datea + "_f" + fhrt + "_pv250_ens.nc")
         lat   = efile.variables['latitude'][:]
         lon   = efile.variables['longitude'][:]
         ens   = np.squeeze(efile.variables['ensemble_data'][:])
         emea  = np.mean(ens, axis=0)
         emea.units = efile.variables['ensemble_data'].units
         evar = np.var(ens, axis=0)

         sens, sigv = computeSens(ens, evar, metric)
         sens[:,:] = sens[:,:] * np.sqrt(evar[:,:])

         outdir = config['work_dir'] + "/" + datea + '.' + config['storm'] + '/' + metname + '/sens/pv250'
         if not os.path.isdir(outdir):
            os.makedirs(outdir)

         if plotDict.get('output_sens', False):
            writeSensFile(lat, lon, fhr, emea, sens, sigv, datea + '/' + bbnn + '/' + datea + "_f" + fhrt + "_pv250hPa_sens.nc", plotDict)

         plotDict['meanCntrs'] = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
         plotScalarSens(lat, lon, sens, emea, sigv, outdir + '/' + datea + "_f" + fhrt + "_pv250hPa_sens.png", plotDict)


      #  Read 500 hPa height, compute sensitivity to that field
      efile = nc.Dataset(config['work_dir'] + "/" + datea + "_f" + fhrt + "_h500_ens.nc")
      eVar = efile.variables

      lat  = eVar['latitude'][:]
      lon  = eVar['longitude'][:]
      ens  = np.squeeze(eVar['ensemble_data'][:])
      emea = np.mean(ens, axis=0)
      emea.units = eVar['ensemble_data'].units
      evar = np.var(ens, axis=0)

      sens, sigv = computeSens(ens, evar, metric)
      sens[:,:] = sens[:,:] * np.sqrt(evar[:,:])

      outdir = config['work_dir'] + "/" + datea + '.' + config['storm'] + '/' + metname + '/sens/h500'
      if not os.path.isdir(outdir):
        os.makedirs(outdir)

      if plotDict.get('output_sens', False):
        writeSensFile(lat, lon, fhr, emea, sens, sigv, datea + '/' + bbnn + '/' + datea + "_f" + fhrt + "_h500hPa_sens.nc", plotDict)

      plotDict['meanCntrs'] = np.array([5400, 5460, 5520, 5580, 5640, 5700, 5760, 5820, 5880, 5940])
      plotScalarSens(lat, lon, sens, emea, sigv, outdir + '/' + datea + "_f" + fhrt + "_h500hPa_sens.png", plotDict)

#      outdir = config['work_dir'] + "/" + datea + '.' + config['storm'] + '/' + metname + '/sens_sc/h500'
#      if not os.path.isdir(outdir):
#        os.makedirs(outdir)

