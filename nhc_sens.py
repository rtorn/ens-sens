import os
import logging
import sys
import netCDF4 as nc
import numpy as np
from SensPlotRoutines import plotVecSens, plotScalarSens, computeSens, writeSensFile

class ComputeSensitivity:
    '''
    This class is the workhorse of the code because it computes the sensitivity of a given forecast
    metric, which is computed for each ensemble member earlier in the code, to a set of forecast 
    fields at a given forecast hour.  These forecast fields were also computed and placed in separate
    netCDF files in an earlier routine.  The result of this routine is a set of sensitivity graphics
    that can be used for various purposes, and if desired netCDF files that can be ingested into
    AWIPS, or into the traveling salesman software for flight planning.  Note that this routine is 
    called several times within the main code; therefore, it could be parallelized.

    Attributes:
        datea   (string):  initialization date of the forecast (yyyymmddhh format)
        fhr        (int):  forecast hour
        metname (string):  name of forecast metric to compute sensitivity for
        atcf     (class):  ATCF class object that includes ensemble information
        config   (dict.):  dictionary that contains configuration options (read from file)
    '''

    def __init__(self, datea, fhr, metname, atcf, config):

      self.df_files = atcf.atcf_files
      self.atcf_array = atcf.atcf_array
      self.nens = int(len(atcf.atcf_files))

      plotDict = {}
      for key in config['sens']:
        plotDict[key] = config['sens'][key]

      fhrt = '%0.3i' % fhr

      logging.warning('Sensitivity of {0} to F{1}'.format(metname,fhrt))

      if config['storm'][-1] == "l":
        bbnn = 'al{0}'.format(config['storm'][-3:-1])
      elif config['storm'][-1] == "e":
        bbnn = 'ep{0}'.format(config['storm'][-3:-1])
      elif storm[-1] == "w":
        bbnn = 'wp{0}'.format(config['storm'][-3:-1])

      #  Compute the ensemble-mean lat/lon for plotting
      elat, elon = atcf.ens_lat_lon_time(fhr)

      m_lat = 0.0
      m_lon = 0.0
      e_cnt = 0.0
      for n in range(self.nens):
        if elat[n] != atcf.missing and elon[n] != atcf.missing:
          e_cnt = e_cnt + 1
          m_lat = m_lat + elat[n]
          m_lon = m_lon + elon[n]

      plotDict['tcLat']     = m_lat / e_cnt
      plotDict['tcLon']     = m_lon / e_cnt
      plotDict['plotTitle'] = '{0} F{1}'.format(datea,fhrt)
      plotDict['fileTitle'] = 'TEST JHT-Torn ECMWF Sensitivity'
      plotDict['initDate']  = '{0}-{1}-{2} {3}:00:00'.format(datea[0:4],datea[4:6],datea[6:8],datea[8:10])

      #  Obtain the metric information (here read from file)
      try:
         mfile = nc.Dataset('{0}/{1}_{2}.nc'.format(config['work_dir'],datea,metname))
      except IOError:
         logging.error('{0}/{1}_{2}.nc does not exist'.format(config['work_dir'],datea,metname))
         return

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

      #  Create dictionary for storm-centered figures
      stceDict = plotDict.copy()
      stceDict['output_sens']=False
      stceDict['range_rings']='True'
      stceDict['ring_center_lat']=float(plotDict['tcLat'])
      stceDict['ring_center_lon']=float(plotDict['tcLon'])
      stceDict['min_lat']=float(plotDict['tcLat'])-float(plotDict.get('storm_center_radius', 10.))
      stceDict['max_lat']=float(plotDict['tcLat'])+float(plotDict.get('storm_center_radius', 10.))
      stceDict['min_lon']=float(plotDict['tcLon'])-float(plotDict.get('storm_center_radius', 10.))
      stceDict['max_lon']=float(plotDict['tcLon'])+float(plotDict.get('storm_center_radius', 10.))
      stceDict['grid_interval']=3.
      stceDict['barb_interval']=3
      stceDict["figsize"]=(8.5,11)

      #  Read the ensemble zonal and meridional steering wind, compute ensemble mean
      ufile = nc.Dataset('{0}/{1}_f{2}_usteer_ens.nc'.format(config['work_dir'],datea,fhrt))
      uVar = ufile.variables

      lat  = uVar['latitude'][:]
      lon  = uVar['longitude'][:]
      uens = np.squeeze(uVar['ensemble_data'][:])
      umea = np.mean(uens, axis=0)
      uvar = np.var(uens, axis=0)

      vfile = nc.Dataset('{0}/{1}_f{2}_vsteer_ens.nc'.format(config['work_dir'],datea,fhrt))
      vVar = vfile.variables

      vens = np.squeeze(vVar['ensemble_data'][:])
      vmea = np.mean(vens, axis=0)
      vvar = np.var(vens, axis=0)


      #  Compute sensitivity with respect to zonal steering wind
      sens, sigv = computeSens(uens, uvar, metric)
      sens[:,:] = sens[:,:] * np.sqrt(uvar[:,:])

      outdir = '{0}/{1}.{2}/{3}/sens/usteer'.format(config['work_dir'],datea,config['storm'],metname)
      if not os.path.isdir(outdir):
        os.makedirs(outdir)

      if plotDict.get('output_sens', False):
        if not os.path.isdir('{0}/{1}'.format(datea,bbnn)):
          os.makedirs('{0}/{1}'.format(datea,bbnn))
        writeSensFile(lat, lon, fhr, umea, sens, sigv, '{0}/{1}/{0}_f{2}_usteer_sens.nc'.format(datea,bbnn,fhrt), plotDict)

      plotVecSens(lat, lon, sens, umea, vmea, sigv, '{0}/{1}_f{2}_usteer_sens.png'.format(outdir,datea,fhrt), plotDict)


      #  Compute sensitivity with respect to meridional steering wind
      sens, sigv = computeSens(vens, vvar, metric)
      sens[:,:] = sens[:,:] * np.sqrt(vvar[:,:])

      outdir = '{0}/{1}.{2}/{3}/sens/vsteer'.format(config['work_dir'],datea,config['storm'],metname)
      if not os.path.isdir(outdir):
        os.makedirs(outdir)

      if plotDict.get('output_sens', False):
        writeSensFile(lat, lon, fhr, umea, sens, sigv, '{0}/{1}/{0}_f{2}_vsteer_sens.nc'.format(datea,bbnn,fhrt), plotDict)

      plotVecSens(lat, lon, sens, umea, vmea, sigv, '{0}/{1}_f{2}_vsteer_sens.png'.format(outdir,datea,fhrt), plotDict)


      #  Rotate wind into major axis direction, compute sensitivity to steering wind in that direction
      wens = ivec * uens[:,:,:] + jvec * vens[:,:,:]
      emea = np.mean(wens, axis=0)
      evar = np.var(wens, axis=0)
      sens, sigv = computeSens(wens, evar, metric)
      sens[:,:] = sens[:,:] * np.sqrt(evar[:,:])

      outdir = '{0}/{1}.{2}/{3}/sens/masteer'.format(config['work_dir'],datea,config['storm'],metname)
      if not os.path.isdir(outdir):
        os.makedirs(outdir)

      if plotDict.get('output_sens', False):
        writeSensFile(lat, lon, fhr, emea, sens, sigv, '{0}/{1}/{0}_f{2}_masteer_sens.nc'.format(datea,bbnn,fhrt), plotDict)

      plotVecSens(lat, lon, sens, umea, vmea, sigv, '{0}/{1}_f{2}_masteer_sens.png'.format(outdir,datea,fhrt), plotDict)

      outdir = '{0}/{1}.{2}/{3}/sens_sc/masteer'.format(config['work_dir'],datea,config['storm'],metname)
      if not os.path.isdir(outdir):
        os.makedirs(outdir)

      plotVecSens(lat, lon, sens, umea, vmea, sigv, '{0}/{1}_f{2}_masteer_sens.png'.format(outdir,datea,fhrt), stceDict)


      #  Read 250 hPa PV, compute sensitivity to that field, if the file exists
      ensfile = '{0}/{1}_f{2}_pv250_ens.nc'.format(config['work_dir'],datea,fhrt)
      if os.path.isfile(ensfile): 
      
         efile = nc.Dataset(ensfile)
         lat   = efile.variables['latitude'][:]
         lon   = efile.variables['longitude'][:]
         ens   = np.squeeze(efile.variables['ensemble_data'][:])
         emea  = np.mean(ens, axis=0)
         emea.units = efile.variables['ensemble_data'].units
         evar = np.var(ens, axis=0)

         sens, sigv = computeSens(ens, evar, metric)
         sens[:,:] = sens[:,:] * np.sqrt(evar[:,:])

         outdir = '{0}/{1}.{2}/{3}/sens/pv250'.format(config['work_dir'],datea,config['storm'],metname)
         if not os.path.isdir(outdir):
            os.makedirs(outdir)

         if plotDict.get('output_sens', False):
            writeSensFile(lat, lon, fhr, emea, sens, sigv, '{0}/{1}/{0}_f{2}_pv250hPa_sens.nc'.format(datea,bbnn,fhrt), plotDict)

         plotDict['meanCntrs'] = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0])
         plotScalarSens(lat, lon, sens, emea, sigv, '{0}/{1}_f{2}_pv250hPa_sens.png'.format(outdir,datea,fhrt), plotDict)


      #  Read 500 hPa height, compute sensitivity to that field
      efile = nc.Dataset('{0}/{1}_f{2}_h500_ens.nc'.format(config['work_dir'],datea,fhrt))
      eVar = efile.variables

      lat  = eVar['latitude'][:]
      lon  = eVar['longitude'][:]
      ens  = np.squeeze(eVar['ensemble_data'][:])
      emea = np.mean(ens, axis=0)
      emea.units = eVar['ensemble_data'].units
      evar = np.var(ens, axis=0)

      sens, sigv = computeSens(ens, evar, metric)
      sens[:,:] = sens[:,:] * np.sqrt(evar[:,:])

      outdir = '{0}/{1}.{2}/{3}/sens/h500'.format(config['work_dir'],datea,config['storm'],metname)
      if not os.path.isdir(outdir):
        os.makedirs(outdir)

      if plotDict.get('output_sens', False):
        writeSensFile(lat, lon, fhr, emea, sens, sigv, '{0}/{1}/{0}_f{2}_h500hPa_sens.nc'.format(datea,bbnn,fhrt), plotDict)

      plotDict['meanCntrs'] = np.array([5400, 5460, 5520, 5580, 5640, 5700, 5760, 5820, 5880, 5940])
      plotScalarSens(lat, lon, sens, emea, sigv, '{0}/{1}_f{2}_h500hPa_sens.png'.format(outdir,datea,fhrt), plotDict)

#      outdir = config['work_dir'] + "/" + datea + '.' + config['storm'] + '/' + metname + '/sens_sc/h500'
#      if not os.path.isdir(outdir):
#        os.makedirs(outdir)

