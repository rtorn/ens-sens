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
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from math import radians, degrees, sin, cos, asin, acos, sqrt

#####   Function that adds dropsonde markers to a plot
def addDrop(dropfile, plt, plotDict):

  msize  = plotDict.get('drop_mark_size', 6)
  mcolor = plotDict.get('drop_mark_color', 'black') 
  mtype  = plotDict.get('drop_mark_type', '+')

  try:
    with open(dropfile, 'r') as fdrop:
      intext = fdrop.readlines()

      ndrops = len(intext)-14
      droplat = np.zeros(ndrops)
      droplon = np.zeros(ndrops)

      for i in range(ndrops):
        str1 = intext[13+i]
        droplat[i] = float(str1[7:9]) + float(str1[10:12])/60.
        droplon[i] = float(str1[16:20]) + np.sign(float(str1[16:20]))*float(str1[21:23])/60.

      plt.plot(droplon, droplat, mtype, color=mcolor, markersize=msize)
      fdrop.close()

  except FileNotFoundError:
    pass

#####   Function that adds rawinsonde markers to a plot
def addRawin(rawinfile, plt, plotDict):

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

#####  Function that computes ensemble-based sensitivity between metric and field 
def computeSens(ens, evar, metric):

  sens = np.zeros(evar.shape)
  sigv = np.zeros(evar.shape)
  nens = len(metric)

  for i in range(len(evar[0,:])):
    for j in range(len(evar[:,0])):

      #  Remove mean, compute regression, and t value
      ens[:,j,i] = ens[:,j,i] - np.mean(ens[:,j,i],axis=0)
      sens[j,i]  = np.sum(ens[:,j,i]*metric[:]) / (float(nens-1) * evar[j,i])
      sy         = sum((metric[:]-sens[j,i]*ens[:,j,i])**2)
      sy         = np.sqrt(sy / float(nens-2)) / np.sqrt(evar[j,i]*float(nens-1))
      sigv[j,i]  = abs(sens[j,i]) / sy

  return sens, sigv

def exists(var):
     var_exists = var in locals() or var in globals()
     return var_exists

def great_circle(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    return 6371. * (acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2)))

#####  Function to write sensitivity fields to a netcdf file
def writeSensFile(lat, lon, fhr, emea, sens, sigv, sensfile, plotDict):

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

#####  Function to plot sensitivity of a scalar (contoured ensemble mean)
def plotScalarSens(lat, lon, sens, emea, sigv, fileout, plotDict):

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

  ax = plt.axes(projection=ccrs.PlateCarree())
  states = NaturalEarthFeature(category="cultural", scale="50m",
                               facecolor="none",
                               name="admin_1_states_provinces_shp")
  ax.add_feature(states, linewidth=0.5, edgecolor="black")
  ax.coastlines('50m', linewidth=1.0)
  ax.add_feature(cartopy.feature.LAKES, facecolor='None', linewidth=1.0, edgecolor='black')
  ax.add_feature(cartopy.feature.BORDERS, facecolor='None', linewidth=1.0, edgecolor='black')

  gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=1, color='gray', alpha=0.5, linestyle='-')
  gl.top_labels = None
  gl.left_labels = None
  gl.xlocator = mticker.FixedLocator(np.arange(10.*np.floor(0.1*minLon),10.*np.ceil(0.1*maxLon)+1.,gridInt))
  gl.xformatter = LONGITUDE_FORMATTER
  gl.xlabel_style = {'size': 12, 'color': 'gray'}
  gl.ylocator = mticker.FixedLocator(np.arange(10.*np.floor(0.1*minLat),10.*np.ceil(0.1*maxLat)+1.,gridInt))
  gl.yformatter = LATITUDE_FORMATTER
  gl.ylabel_style = {'size': 12, 'color': 'gray'}

  ax.set_extent([minLon, maxLon, minLat, maxLat], ccrs.PlateCarree())

  addRawin(plotDict.get("rawinsonde_file","null"), plt, plotDict)

  #  create ensemble-mean contours, sensitivity field, and stat. sig
  pltf = plt.contourf(lon[:],lat[:],sens,compd_range,cmap=cmap,extend='both')
  pltm = plt.contour(lon[:],lat[:],emea,plotDict.get('meanCntrs'),linewidths=1.5, colors='k', zorder=10)
  plts = plt.contour(lon[:],lat[:],sigv,[-2.007, 2.007], linewidths=1.0, colors='k')
  plth = plt.contourf(lon[:],lat[:],sigv,[-2.007, 2.007],hatches=['..', None, '..'], colors='none', extend='both')

  if 'plotTitle' in plotDict:
    plt.title(plotDict['plotTitle'])

  if tcLat != -9999. and tcLon != -9999.:
    plt.plot(tcLon, tcLat, 'o', color='black', markersize=10);

  #  Add range rings to the file if desired
  if plotDict.get('range_rings', 'False')=='True' and tcLat != -9999. and tcLon != -9999:
    tcdist = np.zeros(sens.shape)
    tccen = (tcLat, tcLon)
    for i in range(len(lon)):
      for j in range(len(lat)):
        tcdist[j,i] = great_circle(tcLon, tcLat, lon[i], lat[j])

    if 'ring_values' in plotDict:
      rrings = json.loads(plotDict.get('ring_values'))
    else:
      rrings = [500., 1000., 1500.]

    pltrr = plt.contour(lon[:],lat[:],tcdist,rrings,linewidths=1.0, colors='gray')
    lab = plt.clabel(pltrr, [])

  addDrop(plotDict.get("dropsonde_file","null"), plt, plotDict)

  #  Add colorbar to the plot
  cbar = plt.colorbar(pltf, fraction=0.15, aspect=45., pad=0.04, orientation='horizontal')
  cbar.set_ticks(compd_range[1:11])
  cb = plt.clabel(pltm, inline_spacing=0.0, fontsize=12, fmt="%1.0f")

  plt.savefig(fileout,format='png',dpi=120,bbox_inches='tight')
  plt.close(fig)


#####  Function to plot sensitivity of a vector component (vector ensemble mean)
def plotVecSens(lat, lon, sens, umea, vmea, sigv, fileout, plotDict):

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

  ax = plt.axes(projection=ccrs.PlateCarree())
  states = NaturalEarthFeature(category="cultural", scale="50m",
                               facecolor="none",
                               name="admin_1_states_provinces_shp")
  ax.add_feature(states, linewidth=0.5, edgecolor="black")
  ax.coastlines('50m', linewidth=1.0)
  ax.add_feature(cartopy.feature.LAKES, facecolor='None', linewidth=1.0, edgecolor='black')
  ax.add_feature(cartopy.feature.BORDERS, facecolor='None', linewidth=1.0, edgecolor='black')

  gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                    linewidth=1, color='gray', alpha=0.5, linestyle='-')
  gl.top_labels = None
  gl.left_labels = None
  gl.xlocator = mticker.FixedLocator(np.arange(10.*np.floor(0.1*minLon),10.*np.ceil(0.1*maxLon)+1.,gridInt))
  gl.xformatter = LONGITUDE_FORMATTER
  gl.xlabel_style = {'size': 12, 'color': 'gray'}
  gl.ylocator = mticker.FixedLocator(np.arange(10.*np.floor(0.1*minLat),10.*np.ceil(0.1*maxLat)+1.,gridInt))
  gl.yformatter = LATITUDE_FORMATTER
  gl.ylabel_style = {'size': 12, 'color': 'gray'}

  ax.set_extent([minLon, maxLon, minLat, maxLat], ccrs.PlateCarree())

  addRawin(plotDict.get("rawinsonde_file","null"), plt, plotDict)

  #  create ensemble-mean vectors, sensitivity field, and stat. sig
  pltf = plt.contourf(lon[:],lat[:],sens,compd_range,cmap=cmap,extend='both')
  pltm = plt.barbs(lon[::barbInt], lat[::barbInt], umea[::barbInt,::barbInt]*1.94, vmea[::barbInt,::barbInt]*1.94, pivot='middle', length=6, linewidths=0.2, zorder=10)
  plts = plt.contour(lon[:],lat[:],sigv,[-2.007, 2.007], linewidths=0.5, colors='k')
  plth = plt.contourf(lon[:],lat[:],sigv,[-2.007, 2.007],hatches=['..', None, '..'], colors='none', extend='both')

  if 'plotTitle' in plotDict:
    plt.title(plotDict['plotTitle'])

  if tcLat != -9999. and tcLon != -9999.:
    plt.plot(tcLon, tcLat, 'o', color='black', markersize=10);

  #  Add range rings to the file if desired
  if plotDict.get('range_rings', 'False')=='True' and tcLat != -9999. and tcLon != -9999:
    tcdist = np.zeros(sens.shape)
    tccen = (tcLat, tcLon)
    for i in range(len(lon)):
      for j in range(len(lat)):
        tcdist[j,i] = great_circle(tcLon, tcLat, lon[i], lat[j])

    if 'ring_values' in plotDict:
      rrings = json.loads(plotDict.get('ring_values'))
    else:
      rrings = [500., 1000., 1500.]

    pltrr = plt.contour(lon[:],lat[:],tcdist,rrings,linewidths=1.0, colors='gray')
    lab = plt.clabel(pltrr, [])

  addDrop(plotDict.get("dropsonde_file","null"), plt, plotDict)

  #  Add colorbar to the plot
  cbar = plt.colorbar(pltf, fraction=0.15, aspect=45., pad=0.04, shrink=1.0, orientation='horizontal')
  cbar.set_ticks(compd_range[1:11])

  plt.savefig(fileout,format='png',dpi=120,bbox_inches='tight')
  plt.close(fig)

