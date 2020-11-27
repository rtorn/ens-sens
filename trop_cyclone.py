import numpy as np
import os
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy
import cartopy.feature as cf
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import copy
import math
import atcf_tools as atools
import track as tr

def plot_ens_tc_track(atcf, bvital, storm, datea, config):
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
    
    output_dir = config.get('track_output_dir', '.')
    fhrint  = float(config.get('fhrint',6))
    fhrmax  = float(config.get('forecast_hour_max',120))

    ntimes = int(fhrmax / fhrint) + 1
    nens   = len(atcf.atcf_files)  # total ensembles

    plot_ellipse = config.get('plot_ellipse', 'True')
    subcolors = ["Blue", "DarkOrange"]
    ellcol = ["#551A8B", "#00FFFF", "#00EE00", "#FF0000", "#FF00FF", "#551A8B", "#00FFFF", "#00EE00", "#FF0000"]
    plot_best = True

    #  Create basic map figure
    fig = plt.figure(figsize=(11,8.5))

    ax = plt.axes(projection=ccrs.PlateCarree())
    states = NaturalEarthFeature(category="cultural", scale="50m",
                                 facecolor="none",
                                 name="admin_1_states_provinces_shp")
    ax.add_feature(states, linewidth=0.5, edgecolor="black")
    ax.coastlines('50m', linewidth=1.0)
    ax.add_feature(cartopy.feature.LAKES, facecolor='None', linewidth=1.0, edgecolor='black')
    ax.add_feature(cartopy.feature.BORDERS, facecolor='None', linewidth=1.0, edgecolor='black')

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

    #  Add lat/lon grid lines
    grid_int = float(config.get('grid_interval', 5))
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.top_labels = None
    gl.right_labels = None
    gl.xlocator = mticker.FixedLocator(np.arange(10.*np.floor(0.1*minLon),10.*np.ceil(0.1*maxLon)+1.,grid_int))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.xlabel_style = {'size': 12, 'color': 'gray'}
    gl.ylocator = mticker.FixedLocator(np.arange(10.*np.floor(0.1*minLat),10.*np.ceil(0.1*maxLat)+1.,grid_int))
    gl.yformatter = LATITUDE_FORMATTER
    gl.ylabel_style = {'size': 12, 'color': 'gray'}

    ax.set_extent([minLon, maxLon, minLat, maxLat])

    #  Plot each of the ensemble members
    for n in range(nens):
        x = []
        y = []
        for t in range(ntimes):
            if all_lat[n,t] != atcf.missing and all_lon[n,t] != atcf.missing:
                y.append(all_lat[n,t])
                x.append(all_lon[n,t])
        ax.plot(x, y, color='gray', zorder=1, transform=ccrs.Geodetic())
    t = 1
    bstart = 0
    d = bvital.dim[0]
    bcnt = 0

    # Plot best track position if desired (needs updating)
    for i in range(d):
        if bvital.get_value([i, 6]) == math.inf:
            bcnt = bcnt + 1
    btype = bvital.get_value([0, 6])

    while t <= bcnt and plot_best:
        if bvital.get_value([t, 6]) != btype or t == bcnt:
            if btype > 0.5:
                linestyle = '-'
            else:
                linestyle = '-.'
            x = []
            y = []
            for ib in range(bstart, t):
                if float(bvital.get_value([ib, 0])) != 0.0 and float(bvital.get_value([ib, 1])) != 0.0:
                    y.append(bvital.get_value([ib, 0]))
                    x.append(bvital.get_value([ib, 1]))

            ax.plot(x, y, linestyle=linestyle, color='black', zorder=1, transform=ccrs.Geodetic())
            bstart = t
            btype = bvital.get_value([t, 6])
        t = t + 1

    #  Plot individual member positions and ellipse if desired
    if plot_ellipse == 'True':

       color_index = 0
       x_ell = np.empty(360)
       y_ell = np.empty(360)
       e_lat = np.empty(nens)
       e_lon = np.empty(nens)
       ell_freq = float(config.get('ellipse_frequency', 24))
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
                pb[0,0] = pb[0,0] + (e_lon[n] - m_lon)**2
                pb[1,1] = pb[1,1] + (e_lat[n] - m_lat)**2 
                pb[1,0] = pb[1,0] + (e_lon[n] - m_lon) * (e_lat[n] - m_lat)
             pb[0,1] = pb[1,0]

             pb[:,:] = pb[:,:] / float(pcnt-1)
             rho = pb[1,0] / (math.sqrt(pb[0,0]) * math.sqrt(pb[1,1]))
             sigma_x = math.sqrt(pb[0,0])
             sigma_y = math.sqrt(pb[1,1])
             fac = 1. / (2. * (1. - rho * rho))

             #  Loop over each radian, find the radius that is consistent with the 90% contour
             rdex = 0
             for rad in range(int(math.degrees(2 * math.pi))):
                x_start = math.cos(math.radians(rad))
                y_start = math.sin(math.radians(rad))
                for r_distance in range(2400):
                   x_loc = x_start * r_distance / 80.0
                   y_loc = y_start * r_distance / 80.0
                   prob = math.exp(-1.0 * fac * ((x_loc / sigma_x) ** 2 + (y_loc / sigma_y) ** 2 -
                                               2.0 * rho * (x_loc / sigma_x) * (y_loc / sigma_y)))
                   if prob < 0.256:
                      x_ell[rdex] = x_loc + m_lon
                      y_ell[rdex] = y_loc + m_lat
                      rdex = rdex + 1
                      break
             ax.plot(x_ell, y_ell, color=ellcol[color_index], zorder=12, transform=ccrs.Geodetic())
             color_index = color_index + 1
    
    plt.title("{0} ECMWF forecast of {1}".format(str(datea), storm))
    
    try:   # Create target Directory
        os.makedirs(output_dir)
    except FileExistsError:
        pass

    #  Create the output plot, which is the result of this script
    plt.savefig('{0}/{1}'.format(output_dir,config.get('trackfile','{0}_{1}_track.png'.format(str(datea),str(storm)))), \
                        format='png',dpi=150,bbox_inches='tight') 
    plt.close()


def plot_ens_tc_intensity(atcf, b_vital, storm, datea, config):
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

    output_dir = config.get('int_output_dir', '.')
    fhrint  = float(config.get('fhrint',6))
    fhrmax  = float(config.get('forecast_hour_max',120))

    ntimes = int(fhrmax / fhrint) + 1
    nens   = len(atcf.atcf_files)  # total ensembles

    dim_bvital = b_vital.dim
#    nsub = subens[0]

#    tmax = ntime - 1
#    for t in range(ntime):
#        for n in range(nens):
#            if (fvital.get_value([n, t, 3]) is None) and (b_vital.get_value([t, 3]) is None):
#                tmax = t - 1
#                break
#    subfhr = tr.Track([nens, ntime])
#    subfhr.create_dimension(copy.copy([nens, ntime]))
#    subslp = tr.Track([nsub, nens, ntime])
#    subslp.create_dimension(copy.copy([nsub, nens, ntime]))
#    subwnd = tr.Track([nsub, nens, ntime])
#    subwnd.create_dimension(copy.copy([nsub, nens, ntime]))

    all_slp  = np.ones([nens, ntimes]) * atcf.missing
    all_wnd  = np.ones([nens, ntimes]) * atcf.missing
    fhrvec = np.empty(ntimes)

    for t in range(ntimes):
       fhr = fhrint * t
       all_slp[:,t], all_wnd[:,t] = atcf.ens_intensity_time(fhr)
       fhrvec[t] = fhr

    fig = plt.figure(figsize=(6, 10))
    grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)

    # print(subens)
#    for t in range(ntime):
#        for s in range(nens):
#            subfhr.set_value([s, t], fvital.get_value([0, t, 5]))
#    for s in range(nsub):
#        cnt = 0
#        for n in range(nens):
#            for k in range(subens[1]):
#                for t in range(ntime):
#                    subslp.set_value([s, cnt, t], fvital.get_value([n, t, 2]))
#                    subwnd.set_value([s, cnt, t], fvital.get_value([n, t, 3]))
#                    cnt = cnt + 1
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

#    for s in range(nsub):
#        mean_slp = []
#        mean_y = []
#        for t in range(ntime):
#            mean_arr = []
#            for n in range(nens):
#                mean_arr.append(subslp.get_value([s, n, t]))
#            mean_slp.append(np.mean(mean_arr))
#            mean_y.append(fvital.get_value([0, t, 5]))
#        ax0.plot(mean_slp, mean_y)

    #  Plot the best track minimum sea-level pressure, if it exists
    b_x_pres = []
    b_y_pres = []
    for i in range(dim_bvital[0]):
        b_x_pres.append(b_vital.get_value([i, 4]))
        b_y_pres.append(b_vital.get_value([i, 2]))
    ax0.plot(b_x_pres, b_y_pres, color='black')

   #  Add plot labels and proper tick marks 
    ax0.set_xlabel("Forecast Hour")
    ax0.set_ylabel("Minimum Pressure (hPa)")
    plt.title("{0} ECMWF forecast of {1}".format(str(datea), storm))
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
    b_x_wind = []
    b_y_wind = []
    for i in range(dim_bvital[0]):
        b_x_wind.append(b_vital.get_value([i, 4]))
        b_y_wind.append(b_vital.get_value([i, 3]))
    ax1.plot(b_x_wind, b_y_wind, color='black')

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
    plt.savefig('{0}/{1}'.format(output_dir,config.get('intfile','{0}_{1}_intensity.png'.format(str(datea),str(storm)))), \
                  format='png',dpi=150,bbox_inches='tight')
    plt.close()


def atcf_ens_tc_vitals(atcf, best_file, datea, storm, fstid, bestmax, config, output_dir):
    '''
    Generic routine that can be called to create generic plots of the ensemble TC track and intensity 
    information.  This routine calls the two indivdual routines that actually create the plots.    

    Attributes:
        atcf   (class):  ATCF class object that includes ensemble information
        datea (string):  Initialization date (yyyymmddhh format)
        storm (string):  TC name that will go into plot
        fstid (string):  Name of the model forecast that is being plotted
        config (dict.):  dictionary that contains configuration options (read from file)
    '''
    
    maxbest = config.get('maxbest',120)
    fhrint  = float(config.get('fhrint',6))
    fhrmax  = float(config.get('forecast_hour_max',120))
    
    maxfhr = int(fhrmax / fhrint) + 1
    nens = len(atcf.atcf_files)  # total ensembles
    bestfile = best_file
    b34wrad = {}
    bvital = tr.Track([maxbest, 7])
    bvital.create_dimension(copy.copy([maxbest, 7]))
    for i in range(maxbest):
        b34wrad.update({i: []})
        for j in range(5):
            if j < 5:
                b34wrad.get(i).append(0.0)
    gdate_p = -1.0
    max_wind = -1.0
    gdatea = datetime.strptime(str(datea), '%Y%m%d%H')
    gdate_bm = gdatea.day + bestmax

    #  This code is old and will likely be removed in future versions once the proper best track reading and query code is in atcf_tools.py
    if bestfile is not None:
        # with open(bestfile, 'r') as fo:
        best_data = bestfile
        total_rows = best_data.shape[0]
        bcnt = -1
        for i in range(total_rows):
            date_f = str(best_data[i][8:18])
            gdate_f = datetime.strptime(date_f, '%Y%m%d%H')
            fhr = (gdate_f - gdatea).days * 24.0 + (gdate_f.hour - gdatea.hour)

            tctype = str(best_data[i][59:61])
            if (gdate_f >= gdatea) and ((fhr % fhrint) == 0) and (fhr <= fhrmax):
                if gdate_f != gdate_p:
                    bcnt = bcnt + 1
                    for n in range(4):
                        b34wrad.get(bcnt)[n] == 0.0
                bvital.set_value([bcnt, 0], float(best_data[i][35:38]) * 0.1)
                if str(best_data[i][45]) == "W":
                    bvital.set_value([bcnt, 1], -float(best_data[i][41:45]) * 0.1)
                else:
                    bvital.set_value([bcnt, 1], float(best_data[i][41:45]) * 0.1)
                bvital.set_value([bcnt, 2], float(best_data[i][53:57]))
                bvital.set_value([bcnt, 3], float(best_data[i][48:51]))
                bvital.set_value([bcnt, 4], fhr)
                bvital.set_value([bcnt, 5], float(best_data[i][16:18]))
                if (tctype == "TD") or (tctype == "SD") or (tctype == "TS") or (
                        tctype == "SS") or (tctype == "HU") or (tctype == "TY") or (tctype == "ST"):
                    bvital.set_value([bcnt, 6], 1.0)
                else:
                    bvital.set_value([bcnt, 6], 0.0)
                if float(best_data[i][63:66]) == 34.0:
                    b34wrad.get(bcnt)[1] = float(best_data[i][73:77])
                    b34wrad.get(bcnt)[2] = float(best_data[i][79:83])
                    b34wrad.get(bcnt)[3] = float(best_data[i][85:89])
                    b34wrad.get(bcnt)[4] = float(best_data[i][91:95])
                gdate_p = gdate_f
                gdate_bm = gdate_f.day + bestmax
    else:
        RuntimeError("File Not Found")

#    f34wrad = {}
#    for i in range(nens):
#        f34wrad.update({i: {}})
#        for j in range(int(maxfhr)):
#            f34wrad.get(i).update({j: []})
#            for k in range(5):
#                f34wrad.get(i).get(j).append(0.0)
#    fhr = 0
#    for t in range(int(maxfhr)):
#        for key in f34wrad.keys():
#            f34wrad.get(key).get(t)[0] = fhr
#        fhr = fhr + float(fhrint)

    plot_ens_tc_track(atcf, bvital, storm, datea, config)
    plot_ens_tc_intensity(atcf, bvital, storm, datea, config)

