import numpy as np
import os
import datetime as dt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import math
from SensPlotRoutines import background_map

def plot_ens_tc_track(atcf, storm, datea, config):
    '''
    Routine that creates a plot of the ensemble TC tracks over an entire forecast.  The plot can be
    customized using the items in the vitals_plot section of the configuration file.  The result of this
    routine is a plot of the ensemble tracks that will be placed in the graphics output directory.

    Attributes:
        atcf   (class):  ATCF class object that includes ensemble information
        storm (string):  TC name that will go into plot
        datea (string):  Initialization date (yyyymmddhh format)
        config (dict.):  dictionary that contains configuration options (read from file)
    '''

    earth_radius = 6378.388
    deg2rad      = np.radians(180.)/180.

    output_dir = config['vitals_plot'].get('track_output_dir', '.')
    fhrint  = float(config['vitals_plot'].get('fhrint',6))
    fhrmax  = float(config['vitals_plot'].get('forecast_hour_max',120))

    ntimes = int(fhrmax / fhrint) + 1
    nens   = len(atcf.atcf_files)  # total ensembles

    subcolors = ["Blue", "DarkOrange"]
    ellcol = ["#551A8B", "#00FFFF", "#00EE00", "#FF0000", "#FF00FF", "#551A8B", "#00FFFF", "#00EE00", "#FF0000"]
    plot_best = eval(config['vitals_plot'].get('plot_best', 'True'))

    minLat =  90.
    maxLat = -90.
    minLon = 360.
    maxLon = -180.

    all_lat  = np.ones([nens, ntimes]) * atcf.missing
    all_lon  = np.ones([nens, ntimes]) * atcf.missing
    fhrvec = np.empty(ntimes) 

    for t in range(ntimes):
       fhr = fhrint * t
       all_lat[:,t], all_lon[:,t] = atcf.ens_lat_lon_time(fhr)
       fhrvec[t] = fhr

    #  Figure out the min/max latitude and longitude values
    for n in range(nens):
      for t in range(ntimes):
        if all_lat[n,t] != atcf.missing and all_lon[n,t] != atcf.missing:
          minLat = min([minLat, all_lat[n,t]])
          maxLat = max([maxLat, all_lat[n,t]])
          minLon = min([minLon, all_lon[n,t]])
          maxLon = max([maxLon, all_lon[n,t]]) 

    minLat = minLat - 2.5
    maxLat = maxLat + 2.5
    minLon = minLon - 2.5
    maxLon = maxLon + 2.5

    trackDict = {}
    trackDict['grid_interval']=config['vitals_plot'].get('grid_interval', 5)

    #  Create basic map figure
    fig = plt.figure(figsize=(11,8.5))
    ax = background_map(config['vitals_plot'].get('projection', 'PlateCarree'), minLon, maxLon, minLat, maxLat, trackDict)

    #  Plot each of the ensemble members
    for n in range(nens):
        x = []
        y = []
        for t in range(ntimes):
            if all_lat[n,t] != atcf.missing and all_lon[n,t] != atcf.missing:
                y.append(all_lat[n,t])
                x.append(all_lon[n,t])
        ax.plot(x, y, color='gray', zorder=1, transform=ccrs.PlateCarree())

    # Plot best track position
    fhrbest = []
    bestlat = []
    bestlon = []
    if atcf.has_best and plot_best:

       init = dt.datetime.strptime(datea, '%Y%m%d%H')

       for t in range(ntimes):

          fhr = fhrint * t
          datef = init + dt.timedelta(hours=fhr)

          lat, lon = atcf.best_lat_lon_time(datef.strftime("%Y%m%d%H"))
          if lat != atcf.missing:

             fhrbest.append(fhr)
             bestlat.append(lat)
             bestlon.append(lon)

       ax.plot(bestlon, bestlat, linestyle='-', color='black', zorder=1, transform=ccrs.PlateCarree()) 
        

#    while t <= bcnt and plot_best:
#        if bvital.get_value([t, 6]) != btype or t == bcnt:
#            if btype > 0.5:
#                linestyle = '-'
#            else:
#                linestyle = '-.'
#            x = []
#            y = []
#            for ib in range(bstart, t):
#                if float(bvital.get_value([ib, 0])) != 0.0 and float(bvital.get_value([ib, 1])) != 0.0:
#                    y.append(bvital.get_value([ib, 0]))
#                    x.append(bvital.get_value([ib, 1]))
#
#            ax.plot(x, y, linestyle=linestyle, color='black', zorder=1, transform=ccrs.Geodetic())
#            bstart = t
#            btype = bvital.get_value([t, 6])
#        t = t + 1

    #  Plot individual member positions and ellipse if desired
    if eval(config['vitals_plot'].get('plot_ellipse', 'True')):

       color_index = 0
       x_ell = np.empty(361)
       y_ell = np.empty(361)
       e_lat = np.empty(nens)
       e_lon = np.empty(nens)
       ell_freq = float(config['vitals_plot'].get('ellipse_frequency', 24))
       pb = np.empty([2,2])

       for t in range(ntimes):
          if (fhrvec[t] % ell_freq) == 0 and fhrvec[t] > 0:

             #  Compute the ensemble-mean lat/lon for the members that have position, plot members
             m_lat = 0.0
             m_lon = 0.0
             pcnt  = 0
             for n in range(nens):
                 if all_lat[n,t] != atcf.missing and all_lon[n,t] != atcf.missing:
                     e_lat[pcnt] = all_lat[n,t]
                     e_lon[pcnt] = all_lon[n,t]
                     m_lat = m_lat + all_lat[n,t]
                     m_lon = m_lon + all_lon[n,t]
                     pcnt = pcnt + 1
             if pcnt > 0:
                ax.scatter(e_lon[0:(pcnt-1)], e_lat[0:(pcnt-1)], s=4, marker='o', color=ellcol[color_index], zorder=12)
             if pcnt <= 2:
                break

             #  Compute the deviations from the ensemble mean lat/lon, including covariance
             m_lat = m_lat / pcnt
             m_lon = m_lon / pcnt
             pb[:,:] = 0.0
             for n in range(pcnt):
                fx      = np.radians(e_lon[n]-m_lon) * earth_radius * np.cos(np.radians(0.5*(e_lat[n] + m_lat)))
                fy      = np.radians(e_lat[n]-m_lat) * earth_radius
                pb[0,0] = pb[0,0] + fx**2
                pb[1,1] = pb[1,1] + fy**2
                pb[1,0] = pb[1,0] + fx*fy
             pb[0,1] = pb[1,0]

             pb[:,:] = pb[:,:] / float(pcnt-1)
             rho = pb[1,0] / (math.sqrt(pb[0,0]) * math.sqrt(pb[1,1]))
             sigma_x = math.sqrt(pb[0,0])
             sigma_y = math.sqrt(pb[1,1])
             fac = 1. / (2. * (1. - rho * rho))

             #  Loop over each radian, find the radius that is consistent with the 90% contour
             rdex = 0
             for rad in range(int(math.degrees(2 * math.pi))+1):
                x_start = math.cos(math.radians(rad))
                y_start = math.sin(math.radians(rad))
                for r_distance in range(4000):
                   x_loc = x_start * r_distance
                   y_loc = y_start * r_distance
                   prob = math.exp(-1.0 * fac * ((x_loc / sigma_x) ** 2 + (y_loc / sigma_y) ** 2 -
                                               2.0 * rho * (x_loc / sigma_x) * (y_loc / sigma_y)))
                   if prob < 0.256:
                      x_ell[rdex] = m_lon + x_loc / (deg2rad*earth_radius*np.cos(np.radians(m_lat)))
                      y_ell[rdex] = m_lat + y_loc / (deg2rad*earth_radius)
                      rdex = rdex + 1
                      break
             ax.plot(x_ell, y_ell, color=ellcol[color_index], zorder=12, transform=ccrs.PlateCarree())
             color_index = color_index + 1
   
    plt.title("{0} {1} forecast of {2}".format(str(datea), config.get('model_src',''), storm))
    
    try:   # Create target Directory
        os.makedirs(output_dir)
    except FileExistsError:
        pass

    #  Create the output plot, which is the result of this script
    plt.savefig('{0}/{1}'.format(output_dir,config['vitals_plot'].get('trackfile','{0}_{1}_track.png'.format(str(datea),str(storm)))), \
                        format='png',dpi=150,bbox_inches='tight') 
    plt.close()


def plot_ens_tc_intensity(atcf, storm, datea, config):
    '''
    Routine that creates a plot of the ensemble TC minimum SLP and maximum wind speed as a function of
    forecast lead time.  The plot can be customized using the items in the vitals_plot section of the 
    configuration file.  The result of this routine is a plot of the ensemble min. SLP and maximum wind
    that will be placed in the graphics output directory.

    Attributes:
        atcf   (class):  ATCF class object that includes ensemble information
        storm (string):  TC name that will go into plot
        datea (string):  Initialization date (yyyymmddhh format)
        config (dict.):  dictionary that contains configuration options (read from file)
    '''

    output_dir = config['vitals_plot'].get('int_output_dir', '.')
    fhrint     = float(config['vitals_plot'].get('fhrint',6))
    fhrmax     = float(config['vitals_plot'].get('forecast_hour_max',120))
    plot_best  = eval(config['vitals_plot'].get('plot_best', 'True'))

    ntimes = int(fhrmax / fhrint) + 1
    nens   = len(atcf.atcf_files)  # total ensembles

    all_slp  = np.ones([nens, ntimes]) * atcf.missing
    all_wnd  = np.ones([nens, ntimes]) * atcf.missing
    fhrvec = np.empty(ntimes)

    for t in range(ntimes):
       fhr = fhrint * t
       all_slp[:,t], all_wnd[:,t] = atcf.ens_intensity_time(fhr)
       fhrvec[t] = fhr

    fhrbest = []
    bestslp = []
    bestwnd = []
    if atcf.has_best and plot_best:

       init = dt.datetime.strptime(datea, '%Y%m%d%H')

       for t in range(ntimes):
       
          fhr = fhrint * t
          datef = init + dt.timedelta(hours=fhr)

          wnd, slp = atcf.best_intensity_time(datef.strftime("%Y%m%d%H"))
          if wnd != atcf.missing:
        
             fhrbest.append(fhr)
             bestslp.append(slp)
             bestwnd.append(wnd)


    fig = plt.figure(figsize=(6, 10))
    grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)

    ax0 = fig.add_subplot(grid[0, 0:])

    #  Plot each ensemble member's minimum SLP trace
    minval = 10000000000
    maxval = -10000000000.
    for n in range(nens):
       sens_x = []
       sens_y = []
       for t in range(ntimes):
          if all_slp[n,t] != atcf.missing:
             sens_x.append(fhrvec[t])
             sens_y.append(all_slp[n,t])
             minval = min([minval, all_slp[n,t]])
             maxval = max([maxval, all_slp[n,t]])
       ax0.plot(sens_x, sens_y, color='gray')

    #  Plot the best track minimum sea-level pressure, if it exists
    if atcf.has_best and plot_best:
       ax0.plot(fhrbest, bestslp, color='black')

   #  Add plot labels and proper tick marks 
    ax0.set_xlabel("Forecast Hour")
    ax0.set_ylabel("Minimum Pressure (hPa)")
    plt.title("{0} {1} forecast of {2}".format(str(datea), config.get('model_src',''), storm))
    plt.xticks(range(0,240,24))
    plt.xlim(0, fhrmax)

    ax1 = fig.add_subplot(grid[1, 0:])

    #  Plot each ensemble member's maximum wind speed
    for n in range(nens):
       sens_x1 = []
       sens_y1 = []
       for t in range(ntimes):
          if all_wnd[n,t] != atcf.missing:
             sens_x1.append(fhrvec[t])
             sens_y1.append(all_wnd[n,t])
             minval = min([minval, all_wnd[n,t]])
             maxval = max([maxval, all_wnd[n,t]])
       ax1.plot(sens_x1, sens_y1, color='gray')

    #  Plot the best track maximum wind speed, if it exists
    if atcf.has_best and plot_best:
       ax1.plot(fhrbest, bestwnd, color='black')

    #  Add plot labels and proper tick marks
    ax1.set_xlabel("Forecast Hour")
    ax1.set_ylabel("Maximum Wind Speed (knots)")
    plt.xticks(range(0,240,24))
    plt.xlim(0, fhrmax)    

    try:
        os.makedirs(output_dir)
    except FileExistsError:
        pass

    #  save the figure to a .png file in output graphics directory
    plt.savefig('{0}/{1}'.format(output_dir,config['vitals_plot'].get('intfile','{0}_{1}_intensity.png'.format(str(datea),str(storm)))), \
                  format='png',dpi=150,bbox_inches='tight')
    plt.close()
