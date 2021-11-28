import os
import numpy as np
import xarray as xr
import json
import numpy as np
import datetime as dt
import logging

import matplotlib
from IPython.core.pylabtools import figsize, getfigs
import importlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from eofs.standard import Eof
from eofs.xarray import Eof as Eof_xarray

from SensPlotRoutines import background_map

#####   Function to compute the great circle distance between two points
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


class ComputeForecastMetrics:
    '''
    Function that computes ensemble-based estimates of TC forecast metrics based on the information
    within the configuration file.  Each of these metrics is stored in a seperate netCDF file that
    is used to compute the sensitivity.

    Attributes:
        datea (string): initialization date of the forecast (yyyymmddhh format)
        atcf   (class):  ATCF class object that includes ensemble information
        config (dict.):  dictionary that contains configuration options (read from file)
    '''

    def __init__(self, datea, storm, atcf, config):

        #  Define class-specific variables
        self.fhr = None
        self.deg2rad = 0.01745
        self.earth_radius = 6378388.
        self.missing = -9999.
        self.deg2km = self.earth_radius * np.radians(1)

        self.nens = int(len(atcf.atcf_files))
        fhr_list = json.loads(config['metric']['metric_hours'])
        self.datea_str = datea
        self.datea = dt.datetime.strptime(datea, '%Y%m%d%H')
        self.datea_s = self.datea.strftime("%m%d%H%M")
        self.outdir = config['work_dir']
        self.storm  = storm

        self.dpp = importlib.import_module(config['io_module'])

        self.config = config
        self.atcf   = atcf
        for self.fhr in fhr_list:

            self.fff = str(self.fhr + 1000)[1:]
            logging.warning('  Computing Forecast Metrics for F{0}'.format(self.fff))

            #  Obtain the TC latitude/longitude during, before and after time
            self.ens_lat,  self.ens_lon  = self.atcf.ens_lat_lon_time(self.fhr)
            self.ens_lat1, self.ens_lon1 = self.atcf.ens_lat_lon_time(self.fhr - 6)
            self.ens_lat2, self.ens_lon2 = self.atcf.ens_lat_lon_time(self.fhr + 6)

            e_cnt = 0
            for n in range(self.nens):
               if self.ens_lat[n] != self.atcf.missing and self.ens_lon[n] != self.atcf.missing:
                   e_cnt = e_cnt + 1
 
            #  Go to the next potential forecast hour if there are not enough members
            if e_cnt <= 1:
               continue

            #  Calculate the distance along major axis and along/across track distance
            self.forecast_maj_track, self.forecast_min_track = self.__f_metric_tc_el_track()
            if self.fhr == 0.0:
                self.forecast_maj_track_0 = self.forecast_maj_track
            self.forecast_al_track, self.forecast_ax_track = self.__f_metric_tc_ax_track()

            #  Determine the along/across track direction relative to major variability axis.
            #  Reorient, so that positive major axis is along and to right of track
            alx = self.forecast_al_track['attrs']['X_DIRECTION_VECTOR']
            aly = self.forecast_al_track['attrs']['Y_DIRECTION_VECTOR']

            axx = self.forecast_ax_track['attrs']['X_DIRECTION_VECTOR']
            axy = self.forecast_ax_track['attrs']['Y_DIRECTION_VECTOR']

            elx = self.forecast_maj_track['attrs']['X_DIRECTION_VECTOR']
            ely = self.forecast_maj_track['attrs']['Y_DIRECTION_VECTOR']

            almag = alx * elx + aly * ely
            axmag = alx * elx + aly * ely

            if abs(almag) > abs(axmag):
               if almag < 0:
                  fmet = self.forecast_maj_track['data_vars']["fore_met_init"]['data']
                  self.forecast_maj_track['data_vars']["fore_met_init"]['data'] = list(-1.0 * np.array(fmet))
                  self.forecast_maj_track['attrs']['X_DIRECTION_VECTOR'] = -elx
                  self.forecast_maj_track['attrs']['Y_DIRECTION_VECTOR'] = -ely
                  fmet = self.forecast_min_track['data_vars']["fore_met_init"]['data']
                  self.forecast_min_track['data_vars']["fore_met_init"]['data'] = list(-1.0 * np.array(fmet))
                  self.forecast_min_track['attrs']['X_DIRECTION_VECTOR'] = -elx
                  self.forecast_min_track['attrs']['Y_DIRECTION_VECTOR'] = -ely
            else:
               if axmag < 0:
                  fmet = self.forecast_maj_track['data_vars']["fore_met_init"]['data']
                  self.forecast_maj_track['data_vars']["fore_met_init"]['data'] = list(-1.0 * np.array(fmet))
                  self.forecast_maj_track['attrs']['X_DIRECTION_VECTOR'] = -elx
                  self.forecast_maj_track['attrs']['Y_DIRECTION_VECTOR'] = -ely
                  fmet = self.forecast_min_track['data_vars']["fore_met_init"]['data']
                  self.forecast_min_track['data_vars']["fore_met_init"]['data'] = list(-1.0 * np.array(fmet))
                  self.forecast_min_track['attrs']['X_DIRECTION_VECTOR'] = -elx
                  self.forecast_min_track['attrs']['Y_DIRECTION_VECTOR'] = -ely

            #  Compute metric that is distance from the ensemble-mean position
            fore_met_min = self.forecast_min_track['data_vars']["fore_met_init"]['data']
            vmin = np.var(fore_met_min)
            fore_met_maj = self.forecast_maj_track['data_vars']["fore_met_init"]['data']
            vmaj = np.var(fore_met_maj)

            #out_stat = list(np.zeros(4))
            #out_stat[0] = np.sqrt(vmaj) / (np.sqrt(vmin) + 0.0001)
            f_max_miss = 0.6667
            f_missing = len(np.where(fore_met_maj == 0.0)[0]) / len(fore_met_maj)
#            if f_missing < f_max_miss or self.fhr <= 0.0:
#               if self.fhr > 0.0:
#                  f_met0 = self.forecast_maj_track_0['data_vars']["fore_met_init"]['data']
#               mdist = np.mean(f_met0)
#               vf_met0 = np.var(f_met0)
#               out_stat[1] = max(np.log(np.sqrt(vmaj) / (max(np.sqrt(vf_met0), 10.0))), 1.0)

            f_met = list(np.zeros(len(fore_met_min)))
            f_met[:] = [np.sqrt((i ** 2) + (j ** 2)) for i, j in zip(fore_met_maj, fore_met_min)]
            self.forecast_m_dist = {'coords': {},
                          'attrs': {'FORECAST_METRIC_LEVEL': '',
                                    'FORECAST_METRIC_NAME': 'ensemble distance',
                                    'FORECAST_METRIC_SHORT_NAME': 'emdist',
                                    'FORECAST_VALID_DATE': str(self.datea)},
                          'dims': {'num_ens': self.nens},
                          'data_vars': {'fore_met_init': {'dims': ('num_ens',),
                                                      'attrs': {'units': 'km',
                                                             'description': 'distance from ensemble-mean position'},
                                                   'data': f_met}}}

            #  Write track-related metrics to netcdf files for future use.
            xr.Dataset.from_dict(self.forecast_maj_track).to_netcdf(
                self.outdir + "/{1}_f{0}_majtrack.nc".format(self.fff, str(self.datea_str)), encoding={'fore_met_init': {'dtype': 'float32'}})
            xr.Dataset.from_dict(self.forecast_min_track).to_netcdf(
                self.outdir + "/{1}_f{0}_mintrack.nc".format(self.fff, str(self.datea_str)), encoding={'fore_met_init': {'dtype': 'float32'}})
            xr.Dataset.from_dict(self.forecast_al_track).to_netcdf(
                self.outdir + "/{1}_f{0}_altrack.nc".format(self.fff, str(self.datea_str)), encoding={'fore_met_init': {'dtype': 'float32'}})
            xr.Dataset.from_dict(self.forecast_ax_track).to_netcdf(
                self.outdir + "/{1}_f{0}_axtrack.nc".format(self.fff, str(self.datea_str)), encoding={'fore_met_init': {'dtype': 'float32'}})
            xr.Dataset.from_dict(self.forecast_m_dist).to_netcdf(
                self.outdir + "/{1}_f{0}_dist.nc".format(self.fff, str(self.datea_str)), encoding={'fore_met_init': {'dtype': 'float32'}})

            #  Calculate various intensity-related metrics
            self.__intensity_metrics() 

        #  Compute integrated position EOF metric
        self.__position_eof()

        #  Compute integrated intensity EOF metric
        self.__intensity_eof()

        #  Compute precipitation EOF metric
        if self.config['metric'].get('wind_speed_eof_metric', 'False') == 'True':
           self.__wind_speed_eof()

        #  Compute mean precipitation metric
        if self.config['metric'].get('precipitation_metric', 'False') == 'True':
           self.__precipitation_mean()

        #  Compute precipitation EOF metric
        if self.config['metric'].get('precipitation_eof_metric', 'False') == 'True':
           self.__precipitation_eof()


    def __f_metric_tc_el_track(self):
        '''
        Function that computes the ensemble member's displacement from the ensemble-mean position in the direction
        of the largest ensemble position variability and in the direction that is normal to it.  The result of
        this function is two xarray objects with the ensemble estimates of these two forecast metrics (these are
        saved into netCDF files in the main routine). 
        '''

        e_cnt = 0.0
        x_mean = 0.0
        y_mean = 0.0
        m_lat = 0.0
        m_lon = 0.0
        x_var = 0.0
        y_var = 0.0
        xy_cov = 0.0
        fx_dir = list(np.zeros(self.nens))
        fy_dir = list(np.zeros(self.nens))
        f_met_maj = list(np.zeros(self.nens))
        f_met_min = list(np.zeros(self.nens))

        #  Compute the ensemble-mean if lat/lon pair is not missing
        for n in range(self.nens):
            if self.ens_lat[n] != self.atcf.missing and self.ens_lon[n] != self.atcf.missing:
                e_cnt = e_cnt + 1
                m_lat = m_lat + self.ens_lat[n]
                m_lon = m_lon + self.ens_lon[n]
        m_lon = m_lon / e_cnt
        m_lat = m_lat / e_cnt

        #  Compute the distance in zonal and meridonal direction from mean
        for n in range(self.nens):
            if self.ens_lat[n] != self.atcf.missing and self.ens_lon[n] != self.atcf.missing:
                fx_dir[n] = (self.ens_lon[n] - m_lon) * \
                            self.deg2km * np.cos(np.radians(0.5 * (self.ens_lat[n] + m_lat)))
                fy_dir[n] = (self.ens_lat[n] - m_lat) * self.deg2km
                x_mean = x_mean + fx_dir[n]
                y_mean = y_mean + fy_dir[n]
            else:
                fx_dir[n] = 0.0
                fy_dir[n] = 0.0

        #  Compute variance and covariance in zonal/meridional distance
        x_mean = x_mean / e_cnt
        y_mean = y_mean / e_cnt
        for n in range(self.nens):
            if self.ens_lat[n] != self.atcf.missing and self.ens_lon[n] != self.atcf.missing:
                x_var = x_var + (fx_dir[n] - x_mean) ** 2
                y_var = y_var + (fy_dir[n] - y_mean) ** 2
                xy_cov = xy_cov + (fx_dir[n] - x_mean) * (fy_dir[n] - y_mean)
        x_var = x_var / (e_cnt - 1.0)
        y_var = y_var / (e_cnt - 1.0)
        xy_cov = xy_cov / (e_cnt - 1.0)

        #  Compute major and minor axis based on variance/covariance
        m_trace = x_var + y_var
        m_det = (x_var * y_var) - (xy_cov * xy_cov)
        eval1 = max(0.5 * m_trace + np.sqrt(m_trace * m_trace / 4.0 - m_det),
                    0.5 * m_trace - np.sqrt(m_trace * m_trace / 4.0 - m_det))
        maj_ax = [0.0, 0.0]
        min_ax = [0.0, 0.0]
        if abs(xy_cov) > 0:
            maj_ax[0] = eval1 - y_var
            maj_ax[1] = xy_cov
        else:
            maj_ax[0] = 1.0
            maj_ax[1] = 0.0
        vec_len = np.sqrt((maj_ax[0] * maj_ax[0] + maj_ax[1] * maj_ax[1]))
        maj_ax[0] = maj_ax[0] / vec_len
        maj_ax[1] = maj_ax[1] / vec_len
        min_ax[0] = maj_ax[1]
        min_ax[1] = -maj_ax[0]

        rand1 = np.random.normal(0.0, 0.1, len(fx_dir))

        #  Compute distance between major/minor axis direction for each member
        for n in range(self.nens):
           f_met_maj[n] = maj_ax[0] * (fx_dir[n] - x_mean) + maj_ax[1] * (fy_dir[n] - y_mean) + rand1[n]
           f_met_min[n] = min_ax[0] * (fx_dir[n] - x_mean) + min_ax[1] * (fy_dir[n] - y_mean) + rand1[n]

        forecast_maj_track = {'coords': {},
                              'attrs': {'FORECAST_METRIC_SHORT_NAME': 'tc_maj_pos',
                                        'FORECAST_METRIC_NAME': 'major track error',
                                        'FORECAST_METRIC_LEVEL': '',
                                        'VERIFICATION': 0.0,
                                        'X_DIRECTION_VECTOR': maj_ax[0],
                                        'Y_DIRECTION_VECTOR': maj_ax[1]},
                              'dims': {'num_ens': self.nens},
                              'data_vars': {'fore_met_init': {'dims': ('num_ens',),
                                                              'attrs': {'description': 'TC major axis track error',
                                                                        'units': 'km', '_FillValue': self.missing},
                                                              'data': f_met_maj}}}
        forecast_min_track = {'coords': {},
                              'attrs': {'FORECAST_METRIC_SHORT_NAME': 'tc_min_pos',
                                        'FORECAST_METRIC_NAME': 'minor track error',
                                        'FORECAST_METRIC_LEVEL': '',
                                        'VERIFICATION': 0.0,
                                        'X_DIRECTION_VECTOR': min_ax[0],
                                        'Y_DIRECTION_VECTOR': min_ax[1]},
                              'dims': {'num_ens': self.nens},
                              'data_vars': {'fore_met_init': {'dims': ('num_ens',),
                                                              'attrs': {'description': 'TC minor axis track error',
                                                                        'units': 'km', '_FillValue': self.missing},
                                                              'data': f_met_min}}}

        return forecast_maj_track, forecast_min_track


    def __f_metric_tc_ax_track(self):
        '''
        Function that computes the ensemble member's displacement from the ensemble-mean position in the along
        and across direction.  The result of this function is two xarray objects with the ensemble estimates of 
        these two forecast metrics (these are saved into netCDF files in the main routine). 
        '''

        m_lat = 0.0
        m_lon = 0.0
        e_cnt = 0.0
        f_met_ax = list(np.zeros(self.nens))
        f_met_al = list(np.zeros(self.nens))

        #  Compute the ensemble-mean position at center time
        for n in range(self.nens):
            if self.ens_lat[n] != self.atcf.missing and self.ens_lon[n] != self.atcf.missing:
                e_cnt = e_cnt + 1
                m_lat = m_lat + self.ens_lat[n]
                m_lon = m_lon + self.ens_lon[n]

        m_lon = m_lon / e_cnt
        m_lat = m_lat / e_cnt

        m_lat1 = 0.0
        m_lon1 = 0.0
        e_cnt = 0.0

        #  Compute mean lat/lon at time before
        for n in range(self.nens):
            if self.ens_lat1[n] != self.atcf.missing and self.ens_lon1[n] != self.atcf.missing:
                e_cnt = + 1
                m_lat1 = + self.ens_lat1[n]
                m_lon1 = + self.ens_lon1[n]
        if e_cnt > 0:
            m_lat1 = m_lat1 / e_cnt
            m_lon1 = m_lon1 / e_cnt
        else:
            m_lat1 = m_lat
            m_lon1 = m_lon

        m_lat2 = 0.0
        m_lon2 = 0.0
        e_cnt = 0.0

        #  Compute mean lat/lon at time after
        for n in range(self.nens):
            if self.ens_lat2[n] != self.atcf.missing and self.ens_lon2[n] != self.atcf.missing:
                e_cnt = e_cnt + 1
                m_lat2 = m_lat2 + self.ens_lat2[n]
                m_lon2 = m_lon2 + self.ens_lon2[n]
        if e_cnt > 0:
            m_lat2 = m_lat2 / e_cnt
            m_lon2 = m_lon2 / e_cnt
        else:
            m_lat2 = m_lat
            m_lon2 = m_lon

        #  Compute the along - track direction
        x_dir = (m_lon2 - m_lon1) * self.deg2km * np.cos(0.5 * np.radians(m_lat1 + m_lat2))
        y_dir = (m_lat2 - m_lat1) * self.deg2km
        v_len = max(np.sqrt(x_dir * x_dir + y_dir * y_dir), 0.00001)

        #  Compute unit vectors in along/across directions
        alx_dir = x_dir / v_len
        aly_dir = y_dir / v_len
        axx_dir = aly_dir
        axy_dir = -alx_dir

        rand1 = np.random.normal(0.0, 0.1, len(self.ens_lat2))

        #  Compute the distance in the along/across directions
        for n in range(self.nens):
            if self.ens_lat[n] != self.atcf.missing and self.ens_lon[n] != self.atcf.missing:
                x_dir = (self.ens_lon[n] - m_lon) \
                        * self.deg2km * np.cos(0.5 * np.radians(self.ens_lat[n] + m_lat))
                y_dir = (self.ens_lat[n] - m_lat) * self.deg2km

                f_met_al[n] = x_dir * alx_dir + y_dir * aly_dir + rand1[n]
                f_met_ax[n] = x_dir * axx_dir + y_dir * axy_dir + rand1[n]
            else:
                f_met_al[n] = 0.0
                f_met_ax[n] = 0.0


        forecast_al_track = {'coords': {},
                             'attrs': {'FORECAST_METRIC_SHORT_NAME': 'tc_al_track',
                                       'FORECAST_METRIC_NAME': 'along track error',
                                       'FORECAST_METRIC_LEVEL': '',
                                       'VERIFICATION': 0.0,
                                       'X_DIRECTION_VECTOR': alx_dir,
                                       'Y_DIRECTION_VECTOR': aly_dir,
                                       'TC_LATITUDE': 0.0,
                                       'TC_LONGITUDE': 0.0},
                             'dims': {'num_ens': self.nens},
                             'data_vars': {'fore_met_init': {'dims': ('num_ens',),
                                                             'attrs': {'description': 'TC along track error',
                                                                       'units': 'km', '_FillValue': self.missing},
                                                             'data': f_met_al}}}
        forecast_ax_track = {'coords': {},
                             'attrs': {'FORECAST_METRIC_SHORT_NAME': 'tc_ax_track',
                                       'FORECAST_METRIC_NAME': 'across track error',
                                       'FORECAST_METRIC_LEVEL': '',
                                       'VERIFICATION': 0.0,
                                       'X_DIRECTION_VECTOR': axx_dir,
                                       'Y_DIRECTION_VECTOR': axy_dir,
                                       'TC_LATITUDE': 0.0,
                                       'TC_LONGITUDE': 0.0},
                             'dims': {'num_ens': self.nens},
                             'data_vars': {'fore_met_init': {'dims': ('num_ens',),
                                                             'attrs': {'description': 'TC across track error',
                                                                       'units': 'km', '_FillValue': self.missing},
                                                             'data': f_met_ax}}}
        
        return forecast_al_track, forecast_ax_track


    def __intensity_metrics(self):
        """
        Routine that computes TC intensity-based forecast metrics from each ensemble member.  Currently, the 
        code computes the minimum SLP (always) and kinetic energy within a certain distance of the TC center
        on a certain pressure level (optional).  The result of this function is the ensemble forecast metrics,
        which are saved to netCDF files.
        """

        #  Read the ATCF information for this lead time
        g1 = self.dpp.ReadGribFiles(self.datea_str, self.fhr, self.config)
        lat_vec, lon_vec = self.atcf.ens_lat_lon_time(self.fhr)

        #  Compute the mean latitude and longitude, if missing, replace with ensemble mean
        e_cnt = 0
        m_lat = 0.0
        m_lon = 0.0
        for n in range(self.nens):
           if lat_vec[n] != self.atcf.missing and lon_vec[n] != self.atcf.missing:
              e_cnt = e_cnt + 1
              m_lat = m_lat + lat_vec[n]
              m_lon = m_lon + lon_vec[n]

        m_lon = m_lon / e_cnt
        m_lat = m_lat / e_cnt

        #  Replace missing lat/lon with the ensemble mean.
        for n in range(self.nens):
           if lat_vec[n] == self.atcf.missing or lon_vec[n] == self.atcf.missing:
              lat_vec[n] = m_lat
              lon_vec[n] = m_lon

        mslp_dll = 2.0
        f_met_slp = list(np.zeros(self.nens))
        for n in range(self.nens):

           #  Read SLP field, compute the minimum SLP within a specified distance of the center
           vDict = {'latitude': (lat_vec[n]-mslp_dll, lat_vec[n]+mslp_dll), 'longitude': (lon_vec[n]-mslp_dll,lon_vec[n]+mslp_dll)}
           vDict = g1.set_var_bounds('sea_level_pressure', vDict)
           f_met_slp[n] = np.min(g1.read_grib_field('sea_level_pressure', n, vDict))*0.01

#        if self.fhr > 0.0:
#            f_met0_slp = xr.open_dataset(self.outdir + '/' + str(self.datea_str) + '_f000_minslp.nc').to_dict()['data_vars']["fore_met_init"]['data']
#            mdist_slp = np.mean(f_met0_slp)
#            vf_met0_slp = np.var(f_met0_slp)

#            v_fmet_slp = np.var(f_met_slp)
#            out_stat[2] = max((np.sqrt(v_fmet_slp) /
#                               (max(np.sqrt(vf_met0_slp), 2.0))
#                               ), 1.0)
        
        f_met_slp_nc = {'coords': {},
                        'attrs': {'FORECAST_METRIC_LEVEL': '',
                                  'FORECAST_METRIC_NAME': 'minimum SLP',
                                  'FORECAST_METRIC_SHORT_NAME': 'minslp',
                                  'FORECAST_VALID_DATE': str(self.datea)},
                        'dims': {'num_ens': self.nens},
                        'data_vars': {'fore_met_init': {'dims': ('num_ens',),
                                                        'attrs': {'units': 'hPa',
                                                                  'description': 'minimum sea-level pressure'},
                                                        'data': f_met_slp}}}

        xr.Dataset.from_dict(f_met_slp_nc).to_netcdf(
                self.outdir + "/{1}_f{0}_minslp.nc".format(self.fff, str(self.datea_str)), encoding={'fore_met_init': {'dtype': 'float32'}})

        if self.config['metric'].get('kinetic_energy_metric', 'True') == 'True':

           ke_dll = 4.0
           ke_radius = self.config['metric'].get('kinetic_energy_radius',200.)
           ke_level  = self.config['metric'].get('kinetic_energy_level',1000.)

           logging.warning('    Computing {0} hPa Kinetic Energy'.format(str(ke_level)))

           fmet_kmetric = np.zeros(self.nens)

           for n in range(self.nens):

              vDict = {'latitude': (lat_vec[n]-ke_dll, lat_vec[n]+ke_dll), 
                       'longitude': (lon_vec[n]-ke_dll, lon_vec[n]+ke_dll), 'isobaricInhPa': (ke_level, ke_level)}
              vDict = g1.set_var_bounds('zonal_wind', vDict)              

              #  Read zonal and meridonal wind within certain distance of TC center
              ul = g1.read_grib_field('zonal_wind', n, vDict).squeeze()
              vl = g1.read_grib_field('meridional_wind', n, vDict).squeeze()

              nlat = len(ul.latitude.values)
              nlon = len(ul.longitude.values)

              lonarr, latarr = np.meshgrid(ul.longitude.values, ul.latitude.values)

              dist = great_circle(lon_vec[n], lat_vec[n], lonarr, latarr)

              #  Compute lat/lon weights, replace with zeros where greater than radius
              awght = np.zeros(dist.shape)
              for j in range(nlat):
                awght[j,:] = np.cos(np.radians(ul.latitude.values[j]))

              awght = np.where(dist <= ke_radius, awght, 0.) 
 
              #  Compute the kinetic energy
              fmet_kmetric[n] = 0.5 * np.sum(awght[:,:] * (ul[:,:]**2 + vl[:,:]**2)) / np.sum(awght)

#           if self.fhr > 0.0:
#
#             f_met0_kmetric = xr.open_dataset(self.outdir + '/' + str(self.datea_str) + '_f000' +
#                                               '_ke_10m.nc').to_dict()['data_vars']["fore_met_init"]['data']
#             vf_met0_kmetric = np.var(f_met0_kmetric)
#
#             v_fmet_kmetric = np.var(fmet_kmetric)
#             out_stat[3] = max((np.sqrt(v_fmet_kmetric) /
#                                (max(np.sqrt(vf_met0_kmetric), 2.0))), 1.0)

           f_met_kmetric_nc = {'coords': {},
                               'attrs': {'FORECAST_METRIC_LEVEL': '',
                                         'FORECAST_METRIC_NAME': 'Kinetic Energy',
                                         'FORECAST_METRIC_SHORT_NAME': 'ke',
                                         'FORECAST_VALID_DATE': str(self.datea)},
                               'dims': {'num_ens': self.nens},
                               'data_vars': {'fore_met_init': {'dims': ('num_ens',),
                                             'attrs': {'units': 'm2/s2',
                                                       'description': '10 m Kinetic Energy averaged '
                                                                      'within 200 km of TC center'},
                                                       'data': fmet_kmetric}}}

           xr.Dataset.from_dict(f_met_kmetric_nc).to_netcdf(
                   self.outdir + "/{1}_f{0}_ke_10m.nc".format(self.fff, str(self.datea_str)), encoding={'fore_met_init': {'dtype': 'float32'}})


    def __position_eof(self):
        '''
        Function that computes time-integrated track metric, which is calculated by taking the EOF of 
        the ensemble latitude and longitude for the lead times specified.  The resulting forecast metric is the 
        principal component of the EOF.  The function also plots a figure showing the TC tracks and the 
        track perturbation that is consistent with the first EOF. 
        '''

        logging.warning('  Computing time-integrated track metric')

        ens_min = int(float(self.config['metric'].get('track_eof_member_frac', 0.5))*float(self.nens))

        ellfreq = 24.0

        esign = self.config['metric'].get('track_eof_esign', 1.0)

        fhr1 = self.config['metric'].get('track_eof_hour_init', 24)
        fint = self.config['metric'].get('track_eof_hour_int', 6)
        fhr2 = self.config['metric'].get('track_eof_hour_final', 120)

        ntimes = int((fhr2-fhr1) / fint) + 1

        p1     = -2
        ensvec = np.zeros((self.nens, 2*ntimes))       

        for t in range(ntimes):

           fhr=fhr1+t*fint
           lat, lon=self.atcf.ens_lat_lon_time(fhr)

           #  Compute the ensemble mean for members that have lat/lon values at this time
           e_cnt   = 0
           m_lat_t = 0.0
           m_lon_t = 0.0
           for n in range(self.nens):
              if lat[n] != self.atcf.missing and lon[n] != self.atcf.missing:
                e_cnt = e_cnt + 1
                m_lat_t = m_lat_t + lat[n]
                m_lon_t = m_lon_t + lon[n]

           #  Only consider this time if a critical number of members are present
           if e_cnt >= ens_min:

              m_lon_t = m_lon_t / e_cnt
              m_lat_t = m_lat_t / e_cnt

              #  Compute distance in x/y directions if member is not missing
              p1 = p1 + 2
              p2 = p1 + 1
              for n in range(self.nens):
                 if lat[n] != self.atcf.missing and lon[n] != self.atcf.missing:
                    ensvec[n,p1] = np.radians(lat[n]-m_lat_t)*self.earth_radius
                    ensvec[n,p2] = np.radians(lon[n]-m_lon_t)*self.earth_radius*np.cos(np.radians(m_lat_t))
                 else:
                    ensvec[n,p1] = 0.0
                    ensvec[n,p2] = 0.0

        #  Compute EOF/PCs of the track perturbations
        solver = Eof(ensvec[:,0:(p2+1)])
        pc1    = np.squeeze(solver.pcs(npcs=1, pcscaling=1))

        pc1[:] = pc1[:] / np.std(pc1)

        f1 = 0
        f2 = 120
        ntimes = int((f2-f1) / 6.) + 1

        m_lat   = np.zeros(ntimes)
        m_lon   = np.zeros(ntimes)
        dx      = np.zeros(ntimes)
        dy      = np.zeros(ntimes)
        ens_lat = np.zeros((self.nens, ntimes))
        ens_lon = np.zeros((self.nens, ntimes))

        #  Loop over all times, determine the perturbation distance in x/y for a 1.0 unit PC
        for t in range(ntimes):

           fhr=f1+t*6
           ens_lat[:,t], ens_lon[:,t]=self.atcf.ens_lat_lon_time(fhr)

           e_cnt = 0
           for n in range(self.nens):
              if ens_lat[n,t] != self.atcf.missing and ens_lon[n,t] != self.atcf.missing:
                e_cnt = e_cnt + 1
                m_lat[t] = m_lat[t] + ens_lat[n,t]
                m_lon[t] = m_lon[t] + ens_lon[n,t]

           if e_cnt > 2:

              m_lon[t] = m_lon[t] / e_cnt
              m_lat[t] = m_lat[t] / e_cnt

              for n in range(self.nens):
                 if ens_lat[n,t] != self.atcf.missing and ens_lon[n,t] != self.atcf.missing:
                    dy[t] = dy[t] + np.radians(ens_lat[n,t]-m_lat[t])*self.earth_radius * pc1[n]
                    dx[t] = dx[t] + np.radians(ens_lon[n,t]-m_lon[t])*self.earth_radius*np.cos(np.radians(m_lat[t])) * pc1[n]

              dy[t] = dy[t] / e_cnt
              dx[t] = dx[t] / e_cnt


        imsum = 0.
        jmsum = 0.
        alsum = 0.
        axsum = 0.

        #  Determine the extent to which the PC is aligned with the along/right of track direction
        for t in range(ntimes):

           fhr = f1+t*6

           if fhr >= fhr1 and fhr <= fhr2:

             t1 = max((t-1,0))
             t2 = min((t+1,ntimes-1))

             aloi = np.radians(m_lon[t2]-m_lon[t1])*self.earth_radius*np.cos(np.radians(0.5*(m_lat[t1]+m_lat[t2])))
             aloj = np.radians(m_lat[t2]-m_lat[t1])*self.earth_radius

             veclen = np.sqrt(aloi*aloi + aloj*aloj)
             aloi   = aloi / veclen
             aloj   = aloj / veclen
             acri   = aloj
             acrj   = -aloi

             veclen = np.sqrt(dx[t]**2 + dy[t]**2)
             peri   = dx[t] / np.max([veclen,0.000001])
             perj   = dy[t] / np.max([veclen,0.000001])

             adist  = aloi*peri + aloj*perj
             xdist  = acri*peri + acrj*perj

             alsum = alsum + adist
             axsum = axsum + xdist

             if abs(adist) > abs(xdist):
               if adist < 0:
                 imsum = imsum - peri
                 jmsum = jmsum - perj
               else:
                 imsum = imsum + peri
                 jmsum = jmsum + perj
             else:
               if xdist < 0:
                 imsum = imsum - peri
                 jmsum = jmsum - perj
               else:
                 imsum = imsum + peri
                 jmsum = jmsum + perj

        #  Flip the sign of the EOF, so positive values are along and to right of track
        veclen = np.sqrt(imsum*imsum + jmsum*jmsum)
        imsum  = imsum / veclen
        jmsum  = jmsum / veclen

        if abs(alsum) >= abs(axsum):
          if alsum < 0.0:
            esign = -esign
        else:
          if axsum < 0.0:
            esign = -esign

        pc1[:] = esign * pc1[:]
        dx[:]  = esign * dx[:]
        dy[:]  = esign * dy[:]

        #  Compute perturbed lat/lon for plotting track EOF
        p_lat   = np.zeros(ntimes)
        p_lon   = np.zeros(ntimes)
        for t in range(ntimes):
          p_lat[t] = m_lat[t] + dy[t] / (self.deg2rad*self.earth_radius)
          p_lon[t] = m_lon[t] + dx[t] / (self.deg2rad*self.earth_radius*np.cos(np.radians(m_lat[t])))


        plot_ellipse = self.config['vitals_plot'].get('plot_ellipse',True)
        ell_freq = self.config['vitals_plot'].get('ellipse_frequency', 24)
        ellcol = ["#551A8B", "#00FFFF", "#00EE00", "#FF0000", "#FF00FF", "#551A8B", "#00FFFF", "#00EE00", "#FF0000"]

        minLat =  90.
        maxLat = -90.
        minLon = 360.
        maxLon = -180.

        #  Determine range of figure
        for n in range(self.nens):
          for t in range(ntimes):
            if ens_lat[n,t] != self.atcf.missing and ens_lon[n,t] != self.atcf.missing:
              minLat = min([minLat, ens_lat[n,t]])
              maxLat = max([maxLat, ens_lat[n,t]])
              minLon = min([minLon, ens_lon[n,t]])
              maxLon = max([maxLon, ens_lon[n,t]])

        minLat = minLat - 2.5
        maxLat = maxLat + 2.5
        minLon = minLon - 2.5
        maxLon = maxLon + 2.5

        trackDict = {}
        trackDict['grid_interval'] = self.config['vitals_plot'].get('grid_interval', 5)
        trackDict['left_labels'] = 'True'
        trackDict['right_labels'] = 'None'

        #  Create basic figure plotting options
        fig = plt.figure(figsize=(11,8.5))
        ax = background_map(self.config['vitals_plot'].get('projection', 'PlateCarree'), minLon, maxLon, minLat, maxLat, trackDict)

        x_ell = np.zeros(360)
        y_ell = np.zeros(360)
        pb    = np.zeros((2, 2))

        #  Plot the individual ensemble members
        for n in range(self.nens):
          x = []
          y = []
          for t in range(ntimes):
            if ens_lat[n,t] != self.atcf.missing and ens_lon[n,t] != self.atcf.missing:
              y.append(ens_lat[n,t])
              x.append(ens_lon[n,t])
          ax.plot(x, y, color='lightgray', zorder=1, transform=ccrs.PlateCarree())

        #  Plot the ensemble mean and track perturbation
        ax.plot(m_lon, m_lat, color='black', linewidth=3, zorder=15, transform=ccrs.PlateCarree())        
        ax.plot(p_lon, p_lat, '--', color='black', linewidth=3, zorder=15, transform=ccrs.PlateCarree())

        #  Plot the ellipses and points 
        color_index = 0
        for t in range(ntimes):
          fhr   = f1+t*6

          if (fhr % ell_freq) == 0 and fhr > 0:
            x_ens = []
            y_ens = []
            e_cnt = 0
            for n in range(self.nens):
              if ens_lat[n,t] != self.atcf.missing or ens_lon[n,t] != self.atcf.missing:
                e_cnt = e_cnt + 1
                y_ens.append(ens_lat[n,t])
                x_ens.append(ens_lon[n,t])

            if e_cnt > 2:
              ax.scatter(x_ens, y_ens, s=2, color=ellcol[color_index], zorder=20)
              ax.scatter(m_lon[t], m_lat[t], s=14, color=ellcol[color_index], zorder=20)
              ax.scatter(p_lon[t], p_lat[t], s=14, color=ellcol[color_index], zorder=20)
            else:
              break

            pb[:,:] = 0.0
            for n in range(len(x_ens)):
              pb[0,0] = pb[0,0] + (x_ens[n]-m_lon[t])**2
              pb[1,1] = pb[1,1] + (y_ens[n]-m_lat[t])**2
              pb[1,0] = pb[1,0] + (x_ens[n]-m_lon[t])*(y_ens[n]-m_lat[t])

            pb[0,1] = pb[1,0]
            pb[:,:] = pb[:,:] / float(e_cnt-1)
            rho = pb[1,0] / (np.sqrt(pb[0,0]) * np.sqrt(pb[1,1]))
            sigma_x = np.sqrt(pb[0,0])
            sigma_y = np.sqrt(pb[1,1])
            fac = 1. / (2. * (1. - rho * rho))

            rdex = 0
            for rad in range(int(np.degrees(2*np.pi))):
              x_start = np.cos(np.radians(rad))
              y_start = np.sin(np.radians(rad))
              for r_distance in range(2400):
                x_loc = x_start * r_distance / 80.0
                y_loc = y_start * r_distance / 80.0
                prob = np.exp(-1.0 * fac * ((x_loc / sigma_x) ** 2 + (y_loc / sigma_y) ** 2 -
                                  2.0 * rho * (x_loc / sigma_x) * (y_loc / sigma_y)))
                if prob < 0.256:
                  x_ell[rdex] = x_loc + m_lon[t]
                  y_ell[rdex] = y_loc + m_lat[t]
                  rdex = rdex + 1
                  break

            ax.plot(x_ell, y_ell, color=ellcol[color_index], zorder=20, transform=ccrs.PlateCarree())

            color_index += 1            

        plt.title("{0} {1} forecast of {2}".format(self.datea_str, self.config.get('model_src',''), self.storm))

        outdir = '{0}/{1}.{2}/f{3}_intmajtrack'.format(self.config['work_dir'],self.datea_str,self.storm,'%0.3i' % fhr2)
        if not os.path.isdir(outdir):
           try:
              os.makedirs(outdir)
           except OSError as e:
              raise e

        plt.savefig('{0}/metric.png'.format(outdir), format='png', dpi=150, bbox_inches='tight')
        plt.close(fig)

        #  Create xarray object of forecast metric, write to file.
        f_met_trackeof_nc = {'coords': {},
                             'attrs': {'FORECAST_METRIC_LEVEL': '',
                                       'FORECAST_METRIC_NAME': 'integrated track PC',
                                       'FORECAST_METRIC_SHORT_NAME': 'majtrack',
                                       'X_DIRECTION_VECTOR': imsum,
                                       'Y_DIRECTION_VECTOR': jmsum},
                             'dims': {'num_ens': self.nens},
                             'data_vars': {'fore_met_init': {'dims': ('num_ens',),
                                                            'attrs': {'units': '',
                                                                      'description': 'integrated track PC'},
                                                            'data': np.squeeze(pc1)}}}

        xr.Dataset.from_dict(f_met_trackeof_nc).to_netcdf(
            self.outdir + "/{0}_f{1}_intmajtrack.nc".format(str(self.datea_str), '%0.3i' % fhr2), encoding={'fore_met_init': {'dtype': 'float32'}})


    def __intensity_eof(self):
        '''
        Function that computes time-integrated minimum SLP metric, which is calculated by taking the EOF of 
        the ensemble minimum SLP forecast.  The resulting forecast metric is the principal component of the
        EOF.  The function also plots a figure showing the TC minimum SLP and maximum wind, along with the 
        min. SLP and max. wind perturbation that is consistent with the first EOF. 
        '''

        logging.warning('  Computing time-integrated intensity metric')

        ens_min = int(float(self.config['metric'].get('intensity_eof_member_frac', 0.5))*float(self.nens))

        esign=1.0

        fhr1 = int(self.config['metric'].get('intensity_eof_hour_init', 24))
        fint = int(self.config['metric'].get('intensity_eof_hour_int', 6))
        fhr2 = int(self.config['metric'].get('intensity_eof_hour_final', 96))

        ntimes = int((fhr2-fhr1) / fint) + 1

        ensvec = np.zeros((self.nens, ntimes))
        tt = -1

        #  Loop over all times, calculate ensemble-mean SLP
        for t in range(ntimes):

           fhr=fhr1+t*fint
           slp, wnd=self.atcf.ens_intensity_time(fhr)

           e_cnt   = 0
           m_slp_t = 0.0
           for n in range(self.nens):
              if slp[n] != self.atcf.missing:
                e_cnt = e_cnt + 1
                m_slp_t = m_slp_t + slp[n]

           #   Only consider times where at least half of members have storm for EOF
           if e_cnt >= ens_min:

              m_slp_t = m_slp_t / e_cnt
              tt      = tt + 1

              for n in range(self.nens):
                 if slp[n] != self.atcf.missing:
                    ensvec[n,tt] = slp[n]-m_slp_t
                 else:
                    ensvec[n,tt] = 0.0

        #  Compute the EOF of the MSLP time series
        solver = Eof(ensvec)
        pc1    = np.squeeze(solver.pcs(npcs=1, pcscaling=1))

        pc1[:] = pc1[:] / np.std(pc1)        

        f1 = 0
        f2 = 120
        ntimes = int((f2-f1) / 6.) + 1

        m_fhr   = np.zeros(ntimes)
        m_slp   = np.zeros(ntimes)
        m_wnd   = np.zeros(ntimes)
        dslp    = np.zeros(ntimes)
        dwnd    = np.zeros(ntimes)
        ens_slp = np.zeros((self.nens, ntimes))
        ens_wnd = np.zeros((self.nens, ntimes))
        sumslp  = 0.

        #  Loop over all times, get MSLP, compute mean and EOF perturbation
        for t in range(ntimes):

           fhr=f1+t*6
           ens_slp[:,t], ens_wnd[:,t]=self.atcf.ens_intensity_time(fhr)

           e_cnt = 0
           for n in range(self.nens):
              if ens_slp[n,t] != self.atcf.missing:
                e_cnt = e_cnt + 1
                m_slp[t] = m_slp[t] + ens_slp[n,t]
                m_wnd[t] = m_wnd[t] + ens_wnd[n,t]

           m_fhr[t] = fhr
           if e_cnt > 1:
             m_slp[t] = m_slp[t] / e_cnt
             m_wnd[t] = m_wnd[t] / e_cnt
           else:
             m_slp[t] = None
             m_wnd[t] = None

           #  Compute the MSLP trace associated with a 1 PC perturbation
           for n in range(self.nens):
              if ens_slp[n,t] != self.atcf.missing:
                 dslp[t] = dslp[t] + (ens_slp[n,t]-m_slp[t]) * pc1[n]
                 dwnd[t] = dwnd[t] + (ens_wnd[n,t]-m_wnd[t]) * pc1[n]

           dslp[t] = dslp[t] / e_cnt
           dwnd[t] = dwnd[t] / e_cnt
           sumslp  = sumslp + dslp[t]

        #  Make sure that positive PC is always associated with intensification
        if sumslp > 0.:
          esign = -esign

        pc1[:]  = esign * pc1[:]
        dslp[:] = esign * dslp[:]
        dwnd[:] = esign * dwnd[:]


        #  Create plots of MSLP and maximum wind for each member, mean and EOF perturbation
        fig = plt.figure(figsize=(6, 10))
        grid = plt.GridSpec(2, 2, hspace=0.2, wspace=0.2)

        ax0 = fig.add_subplot(grid[0, 0:])

        minval = 10000000000
        maxval = -10000000000.
        for n in range(self.nens):
          sens_x = []
          sens_y = []
          for t in range(ntimes):
            if ens_slp[n,t] != self.atcf.missing:
              sens_x.append(f1+t*6)
              sens_y.append(ens_slp[n,t])
              minval = min([minval, ens_slp[n,t]])
              maxval = max([maxval, ens_slp[n,t]])
          ax0.plot(sens_x, sens_y, color='lightgray')

        ax0.plot(m_fhr, m_slp, color='black', linewidth=3)
        ax0.plot(m_fhr, m_slp[:]+dslp[:], '--', color='black', linewidth=3)

        ax0.set_xlabel("Forecast Hour")
        ax0.set_ylabel("Minimum Pressure (hPa)")
        plt.title("{0} ECMWF forecast of {1}".format(self.datea_str, self.config['storm']))
        plt.xticks(range(0,240,24))
        plt.xlim(0, 120)

        minval = 10000000000
        maxval = -10000000000.
        ax1 = fig.add_subplot(grid[1, 0:])
        for n in range(self.nens):
          sens_x = []
          sens_y = []
          for t in range(ntimes):
            if ens_wnd[n,t] != self.atcf.missing:
              sens_x.append(f1+t*6)
              sens_y.append(ens_wnd[n,t])
              minval = min([minval, ens_wnd[n,t]])
              maxval = max([maxval, ens_wnd[n,t]])
          ax1.plot(sens_x, sens_y, color='lightgray')

        ax1.plot(m_fhr, m_wnd, color='black', linewidth=3)
        ax1.plot(m_fhr, m_wnd[:]+dwnd[:], '--', color='black', linewidth=3)
 
        ax1.set_xlabel("Forecast Hour")
        ax1.set_ylabel("Maximum Wind Speed (knots)")
        plt.xticks(range(0,240,24))
        plt.xlim(0, 120)

        outdir = '{0}/{1}.{2}/f{3}_intmslp'.format(self.config['work_dir'],self.datea_str,self.storm,'%0.3i' % fhr2)
        if not os.path.isdir(outdir):
           try:
              os.makedirs(outdir)
           except OSError as e:
              raise e

        plt.savefig('{0}/metric.png'.format(outdir), format='png', dpi=150, bbox_inches='tight')
        plt.close(fig)

        f_met_inteneof_nc = {'coords': {},
                             'attrs': {'FORECAST_METRIC_LEVEL': '',
                                       'FORECAST_METRIC_NAME': 'integrated min. SLP PC',
                                       'FORECAST_METRIC_SHORT_NAME': 'intslp'},
                             'dims': {'num_ens': self.nens},
                             'data_vars': {'fore_met_init': {'dims': ('num_ens',),
                                                            'attrs': {'units': '',
                                                                      'description': 'integrated min. SLP PC'},
                                                            'data': pc1}}}

        xr.Dataset.from_dict(f_met_inteneof_nc).to_netcdf(
            self.outdir + "/{0}_f{1}_intmslp.nc".format(str(self.datea_str), '%0.3i' % fhr2), encoding={'fore_met_init': {'dtype': 'float32'}})


    def __precipitation_mean(self):
        '''
        Function that computes precipitation EOF metric, which is calculated by taking the EOF of 
        the ensemble precipitation forecast over a domain defined by the user in a text file.  
        The resulting forecast metric is the principal component of the
        EOF.  The function also plots a figure showing the ensemble-mean precipitation pattern 
        along with the precipitation perturbation that is consistent with the first EOF. 
        '''

        search_max = 150.

        infile = self.config['metric'].get('precip_metric_file').format(self.datea_str,self.storm)
        try:
           f = open(infile, 'r')
        except IOError:
           logging.warning('{0} does not exist.  Cannot compute precip EOF'.format(infile))
           return None

        #  Read the text file that contains information on the precipitation metric
        fhr1 = int(f.readline())
        fhr2 = int(f.readline())
        lat1 = float(f.readline())
        lon1 = float(f.readline())
        lat2 = float(f.readline())
        lon2 = float(f.readline())
        try:
           latc = float(f.readline())
           lonc = float(f.readline())
           maxf = float(f.readline())
        except IOError:
           logging.warning('Mean metric location does not exist.  Cannot compute mean precipitation metric')
           return None

        f.close()

        fff1 = '%0.3i' % fhr1
        fff2 = '%0.3i' % fhr2
        datea_1   = self.datea + dt.timedelta(hours=fhr1)
        date1_str = datea_1.strftime("%Y%m%d%H")
        datea_2   = self.datea + dt.timedelta(hours=fhr2)
        date2_str = datea_2.strftime("%Y%m%d%H")

        #  Read the total precipitation for the beginning of the window
        g1 = self.dpp.ReadGribFiles(self.datea_str, fhr2, self.config)

        vDict = {'latitude': (lat1, lat2), 'longitude': (lon1, lon2),
                 'description': 'precipitation', 'units': 'mm', '_FillValue': -9999.}
        vDict = g1.set_var_bounds('precipitation', vDict)

        ensmat = g1.create_ens_array('precipitation', self.nens, vDict)

        for n in range(self.nens):
           eout = g1.read_grib_field('precipitation', n, vDict).squeeze()
           ensmat[n,:,:] = eout[:,:]

        g1 = self.dpp.ReadGribFiles(self.datea_str, fhr1, self.config)
        if fhr1 > 0.:
           ensmati = g1.create_ens_array('precipitation', self.nens, vDict)
           for n in range(self.nens):
              ensmati[n,:,:] = np.squeeze(g1.read_grib_field('precipitation', n, vDict))
        else:
           ensmati = np.zeros(ensmat.shape)
           ensmati[:,:,:] = 0.

        if eout.units == "m":
           vscale = 1000.
        else:
           vscale = 1.

        #  Scale all of the rainfall to mm and to a 24 h precipitation
        ensmat[:,:,:] = (ensmat[:,:,:] - ensmati[:,:,:]) * vscale * 24. / float(fhr2-fhr1)

        lonarr, latarr = np.meshgrid(eout.longitude.values, eout.latitude.values)
        cdist = great_circle(lonc, latc, lonarr, latarr)
        nlon  = len(eout.longitude.values)
        nlat  = len(eout.latitude.values)

        e_mean = np.mean(ensmat, axis=0)
        e_std  = np.std(ensmat, axis=0)

        #  Blank out all locations outside of search radius
        cdist = np.where(cdist <= search_max, 1.0, 0.0)
        estd_mask = e_std[:,:] * cdist[:,:]

        stdmax = estd_mask.max()
        maxloc = np.where(estd_mask == stdmax)
        icen   = int(maxloc[1])
        jcen   = int(maxloc[0])

        stdmax = stdmax * maxf

        fmgrid = np.zeros(e_mean.shape)
        fmgrid[jcen,icen] = 1.0

        #  start from location of max. SD, search for contiguous points above threshold
        for r in range(1,max(nlon,nlat)):

           i1 = max(icen-r,0)
           i2 = min(icen+r,nlon-1)
           j1 = max(jcen-r,0)
           j2 = min(jcen+r,nlat-1)

           nring = 0
           for i in range(i1+1, i2):
              if ( e_std[j1,i] >= stdmax and fmgrid[j1+1,i] == 1. ):
                 nring = nring + 1
                 fmgrid[j1,i] = 1.
              if ( e_std[j2,i] >= stdmax and fmgrid[j2-1,i] == 1. ):
                 nring = nring + 1
                 fmgrid[j2,i] = 1.

           for j in range(j1+1, j2):
              if ( e_std[j,i1] >= stdmax and fmgrid[j,i1+1] == 1. ):
                 nring = nring + 1
                 fmgrid[j,i1] = 1.
              if ( e_std[j,i2] >= stdmax and fmgrid[j,i2-1] == 1. ):
                 nring = nring + 1
                 fmgrid[j,i2] = 1.

           if ( e_std[j1,i1] >= stdmax and (fmgrid[j1,i1+1] == 1. or fmgrid[j1+1,i1] == 1.)):
              nring = nring + 1
              fmgrid[j1,i1] = 1.

           if ( e_std[j1,i2] >= stdmax and (fmgrid[j1,i2-1] == 1. or fmgrid[j1+1,i2] == 1.)):
              nring = nring + 1
              fmgrid[j1,i2] = 1.

           if ( e_std[j2,i1] >= stdmax and (fmgrid[j2,i1+1] == 1. or fmgrid[j2-1,i1] == 1.)):
              nring = nring + 1
              fmgrid[j2,i1] = 1.

           if ( e_std[j2,i2] >= stdmax and (fmgrid[j2,i2-1] == 1. or fmgrid[j2-1,i2] == 1.)):
              nring = nring + 1
              fmgrid[j2,i2] = 1.

        fmout = np.zeros(g1.nens)
        npts  = np.sum(fmgrid)

        #  Average precipitation
        for n in range(g1.nens):

           fmout[n] = np.sum(fmgrid[:,:]*ensmat[n,:,:]) / npts

        #  Create basic figure, including political boundaries and grid lines
        fig = plt.figure(figsize=(11,6.5), constrained_layout=True)

        colorlist = ("#FFFFFF", "#00ECEC", "#01A0F6", "#0000F6", "#00FF00", "#00C800", "#009000", "#FFFF00", \
                     "#E7C000", "#FF9000", "#FF0000", "#D60000", "#C00000", "#FF00FF", "#9955C9")

        ax1 = plt.subplot(1, 2, 1, projection=ccrs.PlateCarree())
        ax1 = background_map(self.config['metric'].get('projection', 'PlateCarree'), lon1, lon2, lat1, lat2, self.config['metric'])

        #  Plot the ensemble-mean precipitation on the left panel
        mpcp = [0.0, 0.25, 0.50, 1., 1.5, 2., 4., 6., 8., 12., 16., 24., 32., 64., 96., 97.]
        norm = matplotlib.colors.BoundaryNorm(mpcp,len(mpcp))
        pltf1 = plt.contourf(ensmat.longitude.values,ensmat.latitude.values,e_mean,mpcp, \
                              cmap=matplotlib.colors.ListedColormap(colorlist), norm=norm, extend='max')

        cbar = plt.colorbar(pltf1, fraction=0.15, aspect=45., pad=0.04, orientation='horizontal', ticks=mpcp)
        cbar.set_ticks(mpcp[1:(len(mpcp)-1):2])

        plt.title('Mean')

        ax2 = plt.subplot(1, 2, 2, projection=ccrs.PlateCarree())
        ax2 = background_map(self.config['metric'].get('projection', 'PlateCarree'), lon1, lon2, lat1, lat2, self.config['metric'])

        #  Plot the ensemble standard deviation precipitation on the right panel
        spcp = [0., 3., 6., 9., 12., 15., 18., 21., 24., 27., 30., 33., 36., 39., 42., 43.]
        norm = matplotlib.colors.BoundaryNorm(spcp,len(spcp))
        pltf2 = plt.contourf(ensmat.longitude.values,ensmat.latitude.values,e_std,spcp, \
                              cmap=matplotlib.colors.ListedColormap(colorlist), norm=norm, extend='max')

        pltm = plt.contour(ensmat.longitude.values,ensmat.latitude.values,fmgrid,[0.49, 0.51],linewidths=2.5, colors='w')

        cbar = plt.colorbar(pltf2, fraction=0.15, aspect=45., pad=0.04, orientation='horizontal', ticks=spcp)
        cbar.set_ticks(spcp[1:(len(spcp)-1)])

        plt.title('Standard Deviation')

        fig.suptitle('F{0}-F{1} Precipitation ({2}-{3})'.format(fff1, fff2, date1_str, date2_str), fontsize=16)

        outdir = '{0}/{1}.{2}/f{3}_precip'.format(self.config['work_dir'],self.datea_str,self.storm,'%0.3i' % fhr2)
        if not os.path.isdir(outdir):
           try:
              os.makedirs(outdir)
           except OSError as e:
              raise e

        plt.savefig('{0}/metric.png'.format(outdir), format='png', dpi=120, bbox_inches='tight')
        plt.close(fig)

        f_metric = {'coords': {},
                    'attrs': {'FORECAST_METRIC_LEVEL': '',
                              'FORECAST_METRIC_NAME': 'mean precipitation',
                              'FORECAST_METRIC_SHORT_NAME': 'pcp'},
                    'dims': {'num_ens': g1.nens},
                             'data_vars': {'fore_met_init': {'dims': ('num_ens',),
                                                             'attrs': {'units': 'mm',
                                                                       'description': 'precipitation'},
                                                             'data': fmout}}}

        xr.Dataset.from_dict(f_metric).to_netcdf(
            "{0}/{1}_f{2}_precip.nc".format(self.outdir,str(self.datea_str),'%0.3i' % fhr2), encoding={'fore_met_init': {'dtype': 'float32'}})


    def __precipitation_eof(self):
        '''
        Function that computes precipitation EOF metric, which is calculated by taking the EOF of 
        the ensemble precipitation forecast over a domain defined by the user in a text file.  
        The resulting forecast metric is the principal component of the
        EOF.  The function also plots a figure showing the ensemble-mean precipitation pattern 
        along with the precipitation perturbation that is consistent with the first EOF. 
        '''

        infile = self.config['metric'].get('precip_metric_file').format(self.datea_str,self.storm)
        try:
           f = open(infile, 'r')
        except IOError:
           logging.warning('{0} does not exist.  Cannot compute precip EOF'.format(infile))
           return None

        #  Read the text file that contains information on the precipitation metric
        fhr1 = int(f.readline())
        fhr2 = int(f.readline())
        lat1 = float(f.readline())
        lon1 = float(f.readline())
        lat2 = float(f.readline())
        lon2 = float(f.readline())

        f.close()

        #  Read the total precipitation for the beginning of the window
        g1 = self.dpp.ReadGribFiles(self.datea_str, fhr2, self.config)

        vDict = {'latitude': (lat1, lat2), 'longitude': (lon1, lon2),
                 'description': 'precipitation', 'units': 'mm', '_FillValue': -9999.}
        vDict = g1.set_var_bounds('precipitation', vDict)

        ensmat = g1.create_ens_array('precipitation', self.nens, vDict)

        for n in range(self.nens):
           eout = g1.read_grib_field('precipitation', n, vDict).squeeze()
           ensmat[n,:,:] = eout[:,:]

        g1 = self.dpp.ReadGribFiles(self.datea_str, fhr1, self.config)
        if fhr1 > 0.:
           ensmati = g1.create_ens_array('precipitation', self.nens, vDict)
           for n in range(self.nens):
              ensmati[n,:,:] = np.squeeze(g1.read_grib_field('precipitation', n, vDict))
        else:
           ensmati = np.zeros(ensmat.shape)
           ensmati[:,:,:] = 0.

        if eout.units == "m":
           vscale = 1000.
        else:
           vscale = 1.

        #  Scale all of the rainfall to mm and to a 24 h precipitation
        ensmat[:,:,:] = (ensmat[:,:,:] - ensmati[:,:,:]) * vscale * 24. / float(fhr2-fhr1)

        e_mean = np.mean(ensmat, axis=0)
        ensmat = ensmat - e_mean

        #  Compute the EOF of the precipitation pattern and then the PCs
        coslat = np.cos(np.deg2rad(ensmat.latitude.values)).clip(0., 1.)
        wgts = np.sqrt(coslat)[..., np.newaxis]

        solver = Eof_xarray(ensmat.rename({'ensemble': 'time'}), weights=wgts)
        pcout  = solver.pcs(npcs=3, pcscaling=1)
        pc1 = np.squeeze(pcout[:,0])
        pc1[:] = pc1[:] / np.std(pc1)

        #  Compute the precipitation pattern associated with a 1 PC perturbation
        dpcp = np.zeros(e_mean.shape)

        for n in range(self.nens):
           dpcp[:,:] = dpcp[:,:] + ensmat[n,:,:] * pc1[n]

        dpcp[:,:] = dpcp[:,:] / float(self.nens)
  
        if np.sum(dpcp) < 0.:
           pc1[:]    = -pc1[:]
           dpcp[:,:] = -dpcp[:,:]

        gridInt = 5

        #  Create basic figure, including political boundaries and grid lines
        fig = plt.figure(figsize=(11,8.5))

        colorlist = ("#FFFFFF", "#00ECEC", "#01A0F6", "#00BFFF", "#00FF00", "#00C800", "#009000", "#FFFF00", \
                     "#E7C000", "#FF9000", "#FF0000", "#D60000", "#C00000", "#FF00FF", "#9955C9")

        ax = background_map(self.config['metric'].get('projection', 'PlateCarree'), lon1, lon2, lat1, lat2, self.config['metric'])

        mpcp = [0.0, 0.25, 0.50, 1., 1.5, 2., 4., 6., 8., 12., 16., 24., 32., 64., 96., 97.]
        norm = matplotlib.colors.BoundaryNorm(mpcp,len(mpcp))
        pltf = plt.contourf(ensmat.longitude.values,ensmat.latitude.values,e_mean,mpcp, \
                             cmap=matplotlib.colors.ListedColormap(colorlist), norm=norm, alpha=0.5, antialiased=True, extend='max')
        
        pcpfac = np.ceil(np.max(dpcp) / 5.0)
        cntrs = np.array([-5., -4., -3., -2., -1., 1., 2., 3., 4., 5]) * pcpfac
        pltm = plt.contour(ensmat.longitude.values,ensmat.latitude.values,dpcp,cntrs,linewidths=1.5, colors='k', zorder=10)

        #  Add colorbar to the plot
        cbar = plt.colorbar(pltf, fraction=0.15, aspect=45., pad=0.04, orientation='horizontal', ticks=mpcp)
        cbar.set_ticks(mpcp[1:(len(mpcp)-1)])
        cb = plt.clabel(pltm, inline_spacing=0.0, fontsize=12, fmt="%1.0f")

        fracvar = '%4.3f' % solver.varianceFraction(neigs=1)
        plt.title("{0} {1}-{2} hour Precipitation, {3} of variance".format(str(self.datea_str),fhr1,fhr2,fracvar))

        outdir = '{0}/{1}.{2}/f{3}_pcpeof'.format(self.config['work_dir'],self.datea_str,self.storm,'%0.3i' % fhr2)
        if not os.path.isdir(outdir):
           try:
              os.makedirs(outdir)
           except OSError as e:
              raise e

        plt.savefig('{0}/metric.png'.format(outdir), format='png', dpi=120, bbox_inches='tight')
        plt.close(fig)

        f_met_pcpeof_nc = {'coords': {},
                           'attrs': {'FORECAST_METRIC_LEVEL': '',
                                     'FORECAST_METRIC_NAME': 'precipitation PC',
                                     'FORECAST_METRIC_SHORT_NAME': 'pcpeof'},
                             'dims': {'num_ens': self.nens},
                             'data_vars': {'fore_met_init': {'dims': ('num_ens',),
                                                            'attrs': {'units': '',
                                                                      'description': 'precipitation PC'},
                                                            'data': pc1.data}}}

        xr.Dataset.from_dict(f_met_pcpeof_nc).to_netcdf(
            "{0}/{1}_f{2}_pcpeof.nc".format(self.outdir,str(self.datea_str),'%0.3i' % fhr2), encoding={'fore_met_init': {'dtype': 'float32'}})


    def __wind_speed_eof(self):
        '''
        Function that computes wind speed EOF metric, which is calculated by taking the EOF of 
        the ensemble maximum wind speed forecast over a domain defined by the user in a text file.  
        The resulting forecast metric is the principal component of the
        EOF.  The function also plots a figure showing the ensemble-mean max. wind speed pattern 
        along with the wind speed perturbation that is consistent with the first EOF. 
        '''

        infile = self.config['metric'].get('wind_metric_file').format(self.datea_str,self.storm)
        try:
           f = open(infile, 'r')
        except IOError:
           logging.warning('{0} does not exist.  Cannot compute wind EOF'.format(infile))
           return None

        #  Read the text file that contains information on the precipitation metric
        fhr1 = int(f.readline())
        fhr2 = int(f.readline())
        lat1 = float(f.readline())
        lon1 = float(f.readline())
        lat2 = float(f.readline())
        lon2 = float(f.readline())

        f.close()

        #  Read the total precipitation for the beginning of the window
        g1 = self.dpp.ReadGribFiles(self.datea_str, fhr1, self.config)

        vDict = {'latitude': (lat1, lat2), 'longitude': (lon1, lon2),
                 'description': 'wind speed', 'units': 'm/s', '_FillValue': -9999.}
        vDict = g1.set_var_bounds('zonal_wind_10m', vDict)

        ensmat = g1.create_ens_array('zonal_wind_10m', self.nens, vDict)

        for fhr in range(fhr1, fhr2+int(self.config['fcst_hour_int']), int(self.config['fcst_hour_int'])):

           g1 = self.dpp.ReadGribFiles(self.datea_str, fhr, self.config)

           for n in range(self.nens):
              uwnd = g1.read_grib_field('zonal_wind_10m', n, vDict).squeeze()
              vwnd = g1.read_grib_field('meridional_wind_10m', n, vDict).squeeze()
              ensmat[n,:,:] = np.maximum(ensmat[n,:,:],np.sqrt(uwnd[:,:]**2 + vwnd[:,:]**2))

        e_mean = np.mean(ensmat, axis=0)
        ensmat = ensmat - e_mean

        #  Compute the EOF of the precipitation pattern and then the PCs
        coslat = np.cos(np.deg2rad(ensmat.latitude.values)).clip(0., 1.)
        wgts = np.sqrt(coslat)[..., np.newaxis]

        solver = Eof_xarray(ensmat.rename({'ensemble': 'time'}), weights=wgts)
        pcout  = solver.pcs(npcs=3, pcscaling=1)
        pc1 = np.squeeze(pcout[:,0])
        pc1[:] = pc1[:] / np.std(pc1)

        #  Compute the precipitation pattern associated with a 1 PC perturbation
        dwnd = np.zeros(e_mean.shape)

        for n in range(self.nens):
           dwnd[:,:] = dwnd[:,:] + ensmat[n,:,:] * pc1[n]

#        dwnd[:,:] = (dwnd[:,:] / float(self.nens) * units('m/s')) * to('knots')
#        e_mean = (e_mean * units('m/s')) * to('knots') 
        dwnd[:,:] = dwnd[:,:] / float(self.nens) * 1.94
#        e_mean = (e_mean * units('m/s')) * to('knots') 
        e_mean[:,:] = e_mean[:,:] * 1.94

        if np.sum(dwnd) < 0.:
           pc1[:]    = -pc1[:]
           dwnd[:,:] = -dwnd[:,:]

        gridInt = 5

        #  Create basic figure, including political boundaries and grid lines
        fig = plt.figure(figsize=(11,8.5))

        colorlist = ("#FFFFFF", "#00ECEC", "#01A0F6", "#00BFFF", "#00FF00", "#00C800", "#009000", "#FFFF00", \
                     "#E7C000", "#FF9000", "#FF0000", "#D60000", "#C00000", "#FF00FF", "#9955C9")

        ax = background_map(self.config['metric'].get('projection', 'PlateCarree'), lon1, lon2, lat1, lat2, self.config['metric'])

        mwnd = [0.0, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 100, 110]
        norm = matplotlib.colors.BoundaryNorm(mwnd,len(mwnd))
        pltf = plt.contourf(ensmat.longitude.values,ensmat.latitude.values,e_mean,mwnd, \
                             cmap=matplotlib.colors.ListedColormap(colorlist), alpha=0.5, antialiased=True, norm=norm, extend='max')

        wndfac = np.ceil(np.max(dwnd) / 5.0)
        cntrs = np.array([-5., -4., -3., -2., -1., 1., 2., 3., 4., 5]) * wndfac
        pltm = plt.contour(ensmat.longitude.values,ensmat.latitude.values,dwnd,cntrs,linewidths=1.5, colors='k', zorder=10)

        #  Add colorbar to the plot
        cbar = plt.colorbar(pltf, fraction=0.12, aspect=45., pad=0.04, orientation='horizontal', ticks=mwnd)
        cbar.set_ticks(mwnd[1:(len(mwnd)-1)])
        cb = plt.clabel(pltm, inline_spacing=0.0, fontsize=12, fmt="%1.0f")

        fracvar = '%4.3f' % solver.varianceFraction(neigs=1)
        plt.title("{0} {1}-{2} hour Max. Wind Speed, {3} of variance".format(str(self.datea_str),fhr1,fhr2,fracvar))

        outdir = '{0}/{1}.{2}/f{3}_wndeof'.format(self.config['work_dir'],self.datea_str,self.storm,'%0.3i' % fhr2)
        if not os.path.isdir(outdir):
           try:
              os.makedirs(outdir)
           except OSError as e:
              raise e

        plt.savefig('{0}/metric.png'.format(outdir), format='png', dpi=120, bbox_inches='tight')
        plt.close(fig)

        f_met_wndeof_nc = {'coords': {},
                           'attrs': {'FORECAST_METRIC_LEVEL': '',
                                     'FORECAST_METRIC_NAME': 'wind speed PC',
                                     'FORECAST_METRIC_SHORT_NAME': 'wndeof'},
                             'dims': {'num_ens': self.nens},
                             'data_vars': {'fore_met_init': {'dims': ('num_ens',),
                                                            'attrs': {'units': '',
                                                                      'description': 'wind speed PC'},
                                                            'data': pc1.data}}}

        xr.Dataset.from_dict(f_met_wndeof_nc).to_netcdf(
            "{0}/{1}_f{2}_wndeof.nc".format(self.outdir,str(self.datea_str),'%0.3i' % fhr2), encoding={'fore_met_init': {'dtype': 'float32'}})


if __name__ == "__main__":
    src1 = "/Users/parthpatwari/RA_Atmospheric_Science/Old_Code/atcf_data"
    grib_src = "/Users/parthpatwari/RA_Atmospheric_Science/GRIB_files"
    atcf = dpp.Readatcfdata(src1)
    atcf_data = atcf.atcf_files
    no_files = atcf.no_atcf_files
    # g1 = dpp.ReadGribFiles(grib_src, '2019082900', 180)
    ct = ComputeForecastMetrics("ECMWF", '2019082900', atcf.atcf_files, atcf.atcf_array, grib_src)
