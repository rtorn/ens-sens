import matplotlib
from IPython.core.pylabtools import figsize, getfigs
import json
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
import numpy as np
from matplotlib import colors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.path as mpath
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from math import radians, degrees, sin, cos, asin, acos, sqrt


def addDrop(dropfile, plt, plotDict):
  '''
  This function adds markers that represent dropsonde locations based on the *_drop.txt files
  produced by NHC.  Users can specify the style of the markers using configuration file 
  parameters.

  Attributes:
      dropfile (string):  name of dropsonde file to plot
      plt      (object):  matplotlib object to add the dropsonde markers to
      plotDict (dict.):   dictionary that contains configuration options
  '''

  msize  = plotDict.get('drop_mark_size', 6)
  mcolor = plotDict.get('drop_mark_color', 'black') 
  mtype  = plotDict.get('drop_mark_type', '+')

  try:
    with open(dropfile, 'r') as fdrop:
      intext = fdrop.readlines()

      ndrops = len(intext)-14
      droplat = np.zeros(ndrops)
      droplon = np.zeros(ndrops)

      #  Loop over dropsonde locations, convert text to lat/lon, plot
      for i in range(ndrops):
        str1 = intext[13+i]
        droplat[i] = float(str1[7:9]) + float(str1[10:12])/60.
        droplon[i] = float(str1[16:20]) + np.sign(float(str1[16:20]))*float(str1[21:23])/60.

      plt.plot(droplon, droplat, mtype, color=mcolor, markersize=msize)
      fdrop.close()

  except FileNotFoundError:
    pass


def addRangeRings(cen_lat, cen_lon, lat, lon, plt, plotDict):
  '''
  Function that adds range rings from a center lat/lon point to a sensitivity plot.
  Users can specify the interval of the rings and their characteristics.

  Attributes:
      cen_lat  (float):  center latitude value
      cen_lon  (float):  center longitude value
      lat      (float):  vector of model latitude grid values
      lon      (float):  vector of model longitude grid values
      plt     (object):  matplotlib object to add the rawinsonde markers to
      plotDict (dict.):  dictionary that contains configuration options
  '''

  lonarr, latarr = np.meshgrid(lon, lat)
  tcdist = great_circle(cen_lon, cen_lat, lonarr, latarr)

  if 'ring_values' in plotDict:
    rrings = json.loads(plotDict.get('ring_values'))
  else:
    rrings = [500., 1000., 1500.]

  pltrr = plt.contour(lon[:],lat[:],tcdist,rrings,linewidths=1.0, colors='gray')
  lab = plt.clabel(pltrr, [])


def addRawin(rawinfile, plt, plotDict):
  '''
  This function adds markers that represent rawinsonde locations based on a file that is 
  taken from AWIPS.  Users can specify the style of the markers using configuration file 
  parameters.

  Attributes:
      rawinfile (string):  name of rawinsonde file to plot
      plt       (object):  matplotlib object to add the rawinsonde markers to
      plotDict   (dict.):  dictionary that contains configuration options
  '''

  msize  = plotDict.get('rawin_mark_size', 6)
  mcolor = plotDict.get('rawin_mark_color', 'gray')

  try:
    with open(rawinfile,"r") as frawin:
      intext = frawin.readlines()
      rawlat = np.zeros(len(intext))
      rawlon = np.zeros(len(intext))

      for i in range(len(intext)):
        str1 = intext[i]
        rawlat[i] = float(str1[55:60])*0.01
        rawlon[i] = float(str1[61:67])*0.01

      plt.plot(rawlon, rawlat, 'o', color=mcolor, markersize=msize, zorder=5);
      frawin.close()

  except FileNotFoundError:
    pass


def background_map(proj, lon1, lon2, lat1, lat2, DomDict):
  '''
  Function that creates a background map for plotting.

  Attributes:
      proj (string):  String that contains name of projection for plot
      lon1 (float):   minimum longitude of the plot
      lon2 (float):   maximum longitude of the plot
      lat1 (float):   minimum latitude of the plot
      lat2 (float):   maximum latitude of the plot
      DomDict(dict.):  Dictionary that contains plot preferences
  '''

  if proj == 'NorthPolarStereo':
     projinfo = ccrs.NorthPolarStereo(central_longitude=0.0)
  elif (lon1 < 180. and lon2 > 180.):
     projinfo = ccrs.PlateCarree(central_longitude=180.)
  else:
     projinfo = ccrs.PlateCarree()

  if eval(DomDict.get('subplot','False')):

     ax = plt.subplot(DomDict['subrows'], DomDict['subcols'], \
                       DomDict['subnumber'], projection=projinfo)

  else:

     ax = plt.axes(projection=projinfo)
 
  states = NaturalEarthFeature(category="cultural", scale="50m",
                               facecolor="none", name="admin_1_states_provinces")
  ax.add_feature(states, linewidth=0.5, edgecolor="black")
  ax.coastlines('50m', linewidth=1.0)
  ax.add_feature(cartopy.feature.LAKES, facecolor='None', linewidth=1.0, edgecolor='black')
  ax.add_feature(cartopy.feature.BORDERS, facecolor='None', linewidth=1.0, edgecolor='black')

  if proj == 'NorthPolarStereo':

     ax.set_extent([lon1, lon2, lat2, lat1], ccrs.PlateCarree())

     # specifying xlocs/ylocs yields number of meridian/parallel lines
     dmeridian = float(DomDict.get('grid_lon', 30))  # spacing for lines of meridian
     dparallel = float(DomDict.get('grid_lat', 10))  # spacing for lines of parallel 
     num_merid = 360/dmeridian + 1
     num_parra = 30/dparallel + 1
     gl = ax.gridlines(crs=ccrs.PlateCarree(), xlocs=np.linspace(-180, 180, int(num_merid)),
            ylocs=np.linspace(60, 90, int(num_parra)), linestyle="-", linewidth=1, color='gray', alpha=0.5)

     theta = np.linspace(0, 2*np.pi, 120)
     verts = np.vstack([np.sin(theta), np.cos(theta)]).T
     center, radius = [0.5, 0.5], 0.5
     circle = mpath.Path(verts * radius + center)

     ax.set_boundary(circle, transform=ax.transAxes)  #without this; get rect bound 

  else:

     if lon1 > 180. or lon2 > 180.:
       lon1 = lon1 - 360.
       lon2 = lon2 - 360.

     ax.set_extent([lon1, lon2, lat1, lat2], ccrs.PlateCarree())

     gridInt = float(DomDict.get('grid_interval', 10.))

     gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                     linewidth=1, color='gray', alpha=0.5, linestyle='-')
     gl.top_labels = None
     gl.left_labels = eval(DomDict.get('left_labels','None'))
     gl.right_labels = eval(DomDict.get('right_labels','True'))
     gl.xlocator = mticker.FixedLocator(np.arange(-180.,180.,gridInt))
     gl.xformatter = LONGITUDE_FORMATTER
     gl.xlabel_style = {'size': 12, 'color': 'gray'}
     gl.ylocator = mticker.FixedLocator(np.arange(-90.+gridInt,90.,gridInt))
     gl.yformatter = LATITUDE_FORMATTER
     gl.ylabel_style = {'size': 12, 'color': 'gray'}

  return ax


def computeSens(ens, emea, evar, metric):
  '''
  Function that computes the sensitivity of the forecast metric to the provided
  field using the ensemble-based sensitivity technique and the confidence bounds 
  on the regression coefficient at each grid point.

  Attributes:
      ens    (float):  array that contains the ensemble estimate of a field
      evar   (float):  array that is the ensemble variance of forecast field
      metric (float):  vector of ensemble estimates of forecast metric
  '''

  sens = np.zeros(evar.shape)
  sigv = np.zeros(evar.shape)
  nens = len(metric)

  for n in range(nens):
    ens[n,:,:] = ens[n,:,:] - emea[:,:]
    sens[:,:] = sens[:,:] + ens[n,:,:]*metric[n] 
  sens[:,:] = sens[:,:] / (float(nens-1) * evar[:,:])

  for n in range(nens):
    sigv[:,:] = sigv[:,:] + (metric[n]-sens[:,:]*ens[n,:,:])**2
  sigv[:,:] = np.sqrt(sigv[:,:] / float(nens-2)) / np.sqrt(evar[:,:]*float(nens-1))
  sigv[:,:] = abs(sens[:,:]) / sigv[:,:]

  return sens, sigv


def great_circle(lon1, lat1, lon2, lat2):
    '''
    Function that computes the distance between two lat/lon pairs.  The result of this function 
    is the distance in kilometers.

    Attributes
        lon1 (float): longitude of first point
        lat1 (float): latitude of first point
        lon2 (float): longitude of second point.  Can be an array
        lat2 (float): latitude of second point.  Can be an array
    '''

    dist = np.empty(lon2.shape)

    lon1 = np.radians(lon1)
    lat1 = np.radians(lat1)
    lon2[:] = np.radians(lon2[:])
    lat2[:] = np.radians(lat2[:])

    dist[:] = np.sin(lat1) * np.sin(lat2[:]) + np.cos(lat1) * np.cos(lat2[:]) * np.cos(lon1 - lon2[:])

    return 6371. * np.arccos(np.minimum(dist,1.0))


def writeSensFile(lat, lon, fhr, emea, sens, sigv, sensfile, plotDict):
  '''
  Routine that writes a netCDF file that includes the ensemble-mean forecast field, 
  sensitivity, and measure of statistical significance.  This file can be ingested into
  other programs for plotting, such as AWIPS, or could be used within the traveling 
  salesman flight planning software that NHC uses.

  Attributes:
      lat       (float):  Vector of latitude values
      lon       (float):  Vector of longitude values
      fhr         (int):  Forecast hour
      emea      (float):  2D array of the ensemble-mean field
      sens      (float):  2D array of sensitivity field
      sigv      (float):  2D array of statistical significance
      sensfile (string):  Name of netCDF file
      plotDict  (dict.):  Dictionary that contains configuration options
  '''

  #  Create file and dimensions
  ncfile = nc.Dataset(sensfile, mode='w')
  lat_dim = ncfile.createDimension('lat', len(lat))
  lon_dim = ncfile.createDimension('lon', len(lon))
  tim_dim = ncfile.createDimension('time', 1)

  if 'fileTitle' in plotDict:
     ncfile.title       = plotDict.get('fileTitle')
  if 'initDate' in plotDict:
     ncfile.time_origin = plotDict.get('initDate')

  #  Create coordinate variables
  lat_out = ncfile.createVariable('lat', np.float32, ('lat',))
  lat_out.units = 'degrees_north'
  lat_out.long_name = 'latitude'
  lon_out = ncfile.createVariable('lon', np.float32, ('lon',))
  lon_out.units = 'degrees_east'
  lon_out.long_name = 'longitude'
  fhr_out = ncfile.createVariable('forecast_hour', np.float32, ('time',))
  fhr_out.units = 'h'
  fhr_out.long_name = 'forecast_hour'

  #  Create other variables
  emea_out = ncfile.createVariable('ensemble_mean',np.float32,('lat','lon')) #,fill_value=fVar['ensemble_data']._FillValue)
  emea_out.description = 'ensemble mean'
  if hasattr(emea, 'units'):
    emea_out.units = emea.units
  sens_out = ncfile.createVariable('sensitivity',np.float32,('lat','lon'))
  sens_out.description = 'regression coefficient'
  if 'metricUnits' in plotDict:
    sens_out.units = plotDict.get('metricUnits')
  sigv_out = ncfile.createVariable('z_score',np.float32,('lat','lon'))
  sigv_out.description = 'regression coefficient z score'
  sigv_out.units       = ''

  if plotDict.get('nhc_sens', False):
    asen_out = ncfile.createVariable('sensitivity_track_software',np.float32,('lat','lon'))
    asen_out.description = 'abs. value of regression coefficient'
    if 'metricUnits' in plotDict:
      asen_out.units = plotDict.get('metricUnits')

  #  Write variables to a file
  lat_out[:]    = lat
  lon_out[:]    = lon
  fhr_out[:]    = fhr

  emea_out[:,:] = emea
  sens_out[:,:] = sens
  sigv_out[:,:] = sigv
  if plotDict.get('nhc_sens', False):
    asen_out[:,:] = abs(sigv)

  ncfile.close()


def plotScalarSens(lat, lon, sens, emea, sigv, fileout, plotDict):
  '''
  Function that plots the sensitivity of a forecast metric to a scalar field, along
  with the ensemble mean field in contours, and the statistical significance in 
  stippling.  The user has the option to add customized elements to the plot, including
  range rings, locations of rawinsondes/dropsondes, titles, etc.  These are all turned
  on or off using the configuration file.

  Attributes:
      lat      (float):  Vector of latitude values
      lon      (float):  Vector of longitude values
      sens     (float):  2D array of sensitivity field
      emea     (float):  2D array of the ensemble-mean field
      sigv     (float):  2D array of statistical significance
      fileout (string):  Name of output figure in .png format
      plotDict (dict.):  Dictionary that contains configuration options
  '''

#  if np.max(lon) > 180.:
#     lon[:] = (lon[:] + 180.) % 360. - 180.

  minLat = float(plotDict.get('min_lat', np.amin(lat)))
  maxLat = float(plotDict.get('max_lat', np.amax(lat)))
  minLon = float(plotDict.get('min_lon', np.amin(lon)))
  maxLon = float(plotDict.get('max_lon', np.amax(lon)))

  tcLat      = plotDict.get('tcLat', -9999.)
  tcLon      = plotDict.get('tcLon', -9999.)

  gridInt    = float(plotDict.get('grid_interval', 10.))

  colorlist = ("#9A32CD","#00008B","#3A5FCD","#00BFFF","#B0E2FF","#FFFFFF","#FFEC8B","#FFA500","#FF4500","#B22222","#FF82AB")
  cmap = matplotlib.colors.ListedColormap(colorlist)
  compd_range = np.array([ -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2 ]) * plotDict.get('sensmax', 1.)

  #  Create basic figure, including political boundaries and grid lines
  fig = plt.figure(figsize=plotDict.get('figsize',(11,8.5)))

  ax = background_map(plotDict.get('projection', 'PlateCarree'), minLon, maxLon, minLat, maxLat, plotDict)

  addRawin(plotDict.get("rawinsonde_file","null"), plt, plotDict)

  #  create ensemble-mean contours, sensitivity field, and stat. sig
  if plotDict.get('zero_non_sig_sens','False') == 'True':
     sens[sigv < 2.007] = 0.

  pltf = plt.contourf(lon[:],lat[:],sens,compd_range,cmap=cmap,extend='both',transform=ccrs.PlateCarree())
  pltm = plt.contour(lon[:],lat[:],emea,plotDict.get('meanCntrs'),linewidths=1.5, colors='k', zorder=10,transform=ccrs.PlateCarree())
  plts = plt.contour(lon[:],lat[:],sigv,[-2.007, 2.007], linewidths=1.0, colors='k',transform=ccrs.PlateCarree())
  plth = plt.contourf(lon[:],lat[:],sigv,[-2.007, 2.007],hatches=['..', None, '..'], colors='none', extend='both',transform=ccrs.PlateCarree())

  if 'plotTitle' in plotDict:
    plt.title(plotDict['plotTitle'])

  if tcLat != -9999. and tcLon != -9999.:
    plt.plot(tcLon, tcLat, 'o', color='black', markersize=10);

  #  Add range rings to the file if desired
  if plotDict.get('range_rings', 'False')=='True' and \
     'ring_center_lat' in plotDict and 'ring_center_lon' in plotDict:
    addRangeRings(plotDict['ring_center_lat'], plotDict['ring_center_lon'], lat, lon, plt, plotDict)

  addDrop(plotDict.get("dropsonde_file","null"), plt, plotDict)

  #  Add colorbar to the plot
  cbar = plt.colorbar(pltf, fraction=0.15, aspect=45., pad=0.04, orientation='horizontal')
  cbar.set_ticks(compd_range[1:11])
  cb = plt.clabel(pltm, inline_spacing=0.0, fontsize=12, fmt=plotDict.get('clabel_fmt', "%1.0f"))

  plt.savefig(fileout,format='png',dpi=120,bbox_inches='tight')
  plt.close(fig)


def plotVecSens(lat, lon, sens, umea, vmea, sigv, fileout, plotDict):
  '''
  Function that plots the sensitivity of a forecast metric to the component of a vector, 
  along with the ensemble mean wind barbs, and the statistical significance in
  stippling.  The user has the option to add customized elements to the plot, including
  range rings, locations of rawinsondes/dropsondes, titles, etc.  These are all turned
  on or off using the configuration file.

  Attributes:
      lat      (float):  Vector of latitude values
      lon      (float):  Vector of longitude values
      sens     (float):  2D array of sensitivity field
      umea     (float):  2D array of the ensemble-mean zonal wind
      vmea     (float):  2D array of the ensemble-mean meridional wind
      sigv     (float):  2D array of statistical significance
      fileout (string):  Name of output figure in .png format
      plotDict (dict.):  Dictionary that contains configuration options
  '''

#  if np.max(lon) > 180.:
#     lon[:] = (lon[:] + 180.) % 360. - 180.

  minLat = float(plotDict.get('min_lat', np.amin(lat)))
  maxLat = float(plotDict.get('max_lat', np.amax(lat)))
  minLon = float(plotDict.get('min_lon', np.amin(lon)))
  maxLon = float(plotDict.get('max_lon', np.amax(lon)))

  tcLat      = plotDict.get('tcLat', -9999.)
  tcLon      = plotDict.get('tcLon', -9999.)

  gridInt    = float(plotDict.get('grid_interval', 10.))
  barbInt    = int(plotDict.get('barb_interval', 6))

  colorlist = ("#9A32CD","#00008B","#3A5FCD","#00BFFF","#B0E2FF","#FFFFFF","#FFEC8B","#FFA500","#FF4500","#B22222","#FF82AB")
  cmap = matplotlib.colors.ListedColormap(colorlist)
  compd_range = np.array([ -1.2, -1.0, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2 ]) * plotDict.get('sensmax', 1.)

  #  Create basic figure, including political boundaries and grid lines
  fig = plt.figure(figsize=plotDict.get('figsize',(11,8.5)))

  ax = background_map(plotDict.get('projection', 'PlateCarree'), minLon, maxLon, minLat, maxLat, plotDict)

  addRawin(plotDict.get("rawinsonde_file","null"), plt, plotDict)

  #  create ensemble-mean vectors, sensitivity field, and stat. sig
  if plotDict.get('zero_non_sig_sens','False') == 'True':
     sens[sigv < 2.007] = 0.

  pltf = plt.contourf(lon[:],lat[:],sens,compd_range,cmap=cmap,extend='both',transform=ccrs.PlateCarree())
  pltm = plt.barbs(lon[::barbInt], lat[::barbInt], umea[::barbInt,::barbInt]*1.94, vmea[::barbInt,::barbInt]*1.94, \
                    pivot='middle', length=6, linewidths=0.2, zorder=10, transform=ccrs.PlateCarree())
  plts = plt.contour(lon[:],lat[:],sigv,[-2.007, 2.007], linewidths=0.5, colors='k', transform=ccrs.PlateCarree())
  plth = plt.contourf(lon[:],lat[:],sigv,[-2.007, 2.007],hatches=['..', None, '..'], colors='none', extend='both', transform=ccrs.PlateCarree())

  if 'plotTitle' in plotDict:
    plt.title(plotDict['plotTitle'])

  if tcLat != -9999. and tcLon != -9999.:
    plt.plot(tcLon, tcLat, 'o', color='black', markersize=10)

  #  Add range rings to the file if desired
  if plotDict.get('range_rings', 'False')=='True' and \
     'ring_center_lat' in plotDict and 'ring_center_lon' in plotDict:
    addRangeRings(plotDict['ring_center_lat'], plotDict['ring_center_lon'], lat, lon, plt, plotDict)    

  addDrop(plotDict.get("dropsonde_file","null"), plt, plotDict)

  #  Add colorbar to the plot
  cbar = plt.colorbar(pltf, fraction=0.10, aspect=45., pad=0.04, shrink=1.0, orientation='horizontal')
  cbar.set_ticks(compd_range[1:11])

  plt.savefig(fileout,format='png',dpi=120,bbox_inches='tight')
  plt.close(fig)
