import os
import importlib
import math
import netCDF4 as nc
import cfgrib
import sys
import numpy as np
import datetime as dt
#from windspharm.standard import VectorWind
#from windspharm.tools import prep_data
import metpy.constants as mpcon
import metpy.calc as mpcalc
from metpy.units import units

####   Class for computing forecast fields for TC sensitivity calculations
class ComputeTCFields:
    def __init__(self, datea, fhr, atcf, config):
        # Forecast fields to compute

        wnd_lev_1 = [250, 500]
        wnd_lev_2 = [350, 500]
        n_wnd_lev = len(wnd_lev_1)

        # Read steering flow parameters, or use defaults
        # steer info
        steerp1 = 300
        steerp2 = 850.0
        tcradius = 333.0

        # lat_lon info
        lat1 = float(config['fields'].get('min_lat','0.'))
        lat2 = float(config['fields'].get('max_lat','65.'))
        lon1 = float(config['fields'].get('min_lon','-180.'))
        lon2 = float(config['fields'].get('max_lon','-10.'))
        self.fhr = fhr
        self.deg2rad = 0.01745
        self.earth_radius = 6378
        self.deg2km = self.earth_radius * math.radians(1)
        self.atcf_files = atcf.atcf_files
        self.config     = config       
        self.nens = int(len(self.atcf_files))
        df_files = {}
        self.datea_str = datea
        self.datea = dt.datetime.strptime(datea, '%Y%m%d%H')
        self.datea_s = self.datea.strftime("%m%d%H%M")
        self.fff = str(self.fhr + 1000)[1:]
        datea_1 = self.datea + dt.timedelta(hours=self.fhr)
        datea_1 = datea_1.strftime("%m%d%H%M")

        self.dpp = importlib.import_module(config['io_module'])

        print("Computing hour " + self.fff + " ensemble fields")

        #  Obtain the ensemble lat/lon information, replace missing values with mean
        self.ens_lat, self.ens_lon = atcf.ens_lat_lon_time(self.fhr)

        e_cnt = 0
        m_lat = 0.0
        m_lon = 0.0
        for n in range(self.nens):
            if self.ens_lat[n] != atcf.missing and self.ens_lon[n] != atcf.missing:
                e_cnt = e_cnt + 1
                m_lat = m_lat + self.ens_lat[n]
                m_lon = m_lon + self.ens_lon[n]
        m_lon = m_lon / e_cnt
        m_lat = m_lat / e_cnt

        for n in range(self.nens):
            if self.ens_lat[n] == atcf.missing or self.ens_lon[n] == atcf.missing:
                self.ens_lat[n] = m_lat
                self.ens_lon[n] = m_lon

        g1 = self.dpp.ReadGribFiles(self.datea_str, self.fhr, self.config)

        dencode = {'ensemble_data': {'dtype': 'float32'}, 'latitude': {'dtype': 'float32'},
                   'longitude': {'dtype': 'float32'}, 'ensemble': {'dtype': 'int32'}}

        uoutfile=config['work_dir'] + "/" + str(self.datea_str) + '_f' + self.fff + '_usteer_ens.nc'
        voutfile=config['work_dir'] + "/" + str(self.datea_str) + '_f' + self.fff + '_vsteer_ens.nc'
        if (not os.path.isfile(uoutfile) or not os.path.isfile(voutfile)) and config['fields'].get('calc_uvsteer','True') == 'True':

          print("  Computing steering wind information")

          inpDict = {'isobaricInhPa': (steerp1, steerp2)}
          inpDict = g1.set_var_bounds('zonal_wind', inpDict)

          #  Create output arrays
          outDict = {'latitude': (lat1, lat2), 'longitude': (lon1, lon2), 
                     'description': 'zonal steering wind', 'units': 'm/s', '_FillValue': -9999.}
          outDict = g1.set_var_bounds('zonal_wind', outDict)
          uensmat = g1.create_ens_array('zonal_wind', self.nens, outDict)          

          outDict = {'latitude': (lat1, lat2), 'longitude': (lon1, lon2),
                     'description': 'meridional steering wind', 'units': 'm/s', '_FillValue': -9999.}
          outDict = g1.set_var_bounds('meridional_wind', outDict)
          vensmat = g1.create_ens_array('meridional_wind', self.nens, outDict)       

          outDict = {'latitude': (lat1, lat2), 'longitude': (lon1, lon2),
                     'description': 'steering wind vorticity', 'units': '1/s', '_FillValue': -9999.}
          outDict = g1.set_var_bounds('zonal_wind', outDict)
          vortmat = g1.create_ens_array('zonal_wind', self.nens, outDict)

          wencode = {'latitude': {'dtype': 'float32'}, 'longitude': {'dtype': 'float32'}}

          for n in range(self.nens):

             #  Read global zonal and meridional wind, write to file
             uwnd = g1.read_grib_field('zonal_wind', n, inpDict).rename('u')
             vwnd = g1.read_grib_field('meridional_wind', n, inpDict).rename('v')

             uwnd.to_netcdf('wind_info.nc', mode='w', encoding=wencode, format='NETCDF3_CLASSIC')
             vwnd.to_netcdf('wind_info.nc', mode='a', encoding=wencode, format='NETCDF3_CLASSIC')

             #  Call NCL to remove TC winds, read result from file
             os.system('ncl -Q ' + config['script_dir'] + '/tc_steer.ncl tclat=' + str(self.ens_lat[n]) + 
                       ' tclon=' + str(self.ens_lon[n]) + ' tcradius=' + str(tcradius))

             wfile     = nc.Dataset('wind_info.nc')
             uwnd[:,:] = wfile.variables['u'][:,:]
             vwnd[:,:] = wfile.variables['v'][:,:]

             os.remove('wind_info.nc')

             #  Integrate the winds over the layer to obtain the steering wind
#             pres      = np.array(uwnd.isobaricInhPa.data)
             pres,lat,lon = uwnd.indexes.values()
             nlev      = len(pres)

             uint      = uwnd[0,:,:]
             uint[:,:] = 0.0
             vint      = vwnd[0,:,:]
             vint[:,:] = 0.0

             for k in range(nlev-1):

               uint[:,:] = uint[:,:] + 0.5 * (uwnd[k,:,:]+uwnd[k+1,:,:]) * abs(pres[k+1]-pres[k])
               vint[:,:] = vint[:,:] + 0.5 * (vwnd[k,:,:]+vwnd[k+1,:,:]) * abs(pres[k+1]-pres[k])

             if lat[0] > lat[-1]:
               slat1 = lat2
               slat2 = lat1
             else:
               slat1 = lat1
               slat2 = lat2

             uensmat[n,:,:] = np.squeeze(uint.sel(latitude=slice(slat1, slat2), longitude=slice(lon1, lon2))) / abs(pres[nlev-1]-pres[0])
             vensmat[n,:,:] = np.squeeze(vint.sel(latitude=slice(slat1, slat2), longitude=slice(lon1, lon2))) / abs(pres[nlev-1]-pres[0])

             #  Compute the vorticity associated with the steering wind

#             circ = VectorWind(unew, vnew).vorticity() * 1.0e5

#             vortmat[n,:,:] = np.squeeze(circ.sel(latitude=slice(lat2, lat1), longitude=slice(lon1, lon2)))

          uensmat.to_netcdf(uoutfile, encoding=dencode)
          vensmat.to_netcdf(voutfile, encoding=dencode) 
#          vortmat.to_netcdf(vortfile, encoding=dencode)

        else:

          print("  Obtaining steering wind information from file")

        #  Read 500 hPa geopotential height from file, if ensemble file is not present
        outfile=config['work_dir'] + "/" + str(self.datea_str) + '_f' + self.fff + '_h500_ens.nc'
        if (not os.path.isfile(outfile) and config['fields'].get('calc_h500hPa','True') == 'True'):

          print("  Computing 500 hPa height")

          vDict = {'latitude': (lat1, lat2), 'longitude': (lon1, lon2), 'isobaricInhPa': (500, 500), 
                   'description': '500 hPa height', 'units': 'm', '_FillValue': -9999.}
          vDict = g1.set_var_bounds('geopotential_height', vDict)
          ensmat = g1.create_ens_array('geopotential_height', self.nens, vDict)

          for n in range(self.nens):
             ensmat[n,:,:] = np.squeeze(g1.read_grib_field('geopotential_height', n, vDict))

          ensmat.to_netcdf(outfile, encoding=dencode)

        else:

          print("  Obtaining 500 hPa height data from " + outfile)

        #  Compute 250 hPa PV if the file does not exist
        outfile=config['work_dir'] + "/" + str(self.datea_str) + '_f' + self.fff + '_pv250_ens.nc'
        if (not os.path.isfile(outfile) and config['fields'].get('calc_pv250hPa','True') == 'True'):

          print("  Computing 250 hPa PV")

          vDict = {'latitude': (lat1, lat2), 'longitude': (lon1, lon2), 'isobaricInhPa': (200, 300),
                   'description': '250 hPa Potential Vorticity', 'units': 'PVU', '_FillValue': -9999.}
          vDict = g1.set_var_bounds('zonal_wind', vDict)

          ensmat = g1.create_ens_array('zonal_wind', self.nens, vDict)

          for n in range(self.nens):

            #  Read all the necessary files from file, smooth fields, so sensitivities are useful
            tmpk = g1.read_grib_field('temperature', n, vDict) * units('K')

            lats = tmpk.latitude.values
            lons = tmpk.longitude.values
            pres = tmpk.isobaricInhPa.values * units('hPa')

            tmpk = mpcalc.smooth_n_point(tmpk, 9, 4)

            thta = mpcalc.potential_temperature(pres[:, None, None], tmpk)

            uwnd = mpcalc.smooth_n_point(g1.read_grib_field('zonal_wind', n, vDict) * units('m/s'), 9, 4)
            vwnd = mpcalc.smooth_n_point(g1.read_grib_field('meridional_wind', n, vDict) * units('m/s'), 9, 4)

            dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)

            #  Compute PV and place in ensemble array
            pvout = mpcalc.potential_vorticity_baroclinic(thta, pres[:, None, None], uwnd, vwnd,
                                           dx[None, :, :], dy[None, :, :], lats[None, :, None] * units('degrees'))

            ensmat[n,:,:] = np.squeeze(pvout[np.where(pres == 250 * units('hPa'))[0],:,:]) * 1.0e6
 
          ensmat.to_netcdf(outfile, encoding=dencode)

        elif os.path.isfile(outfile):

          print("  Obtaining 250 hPa PV data from " + outfile)

        outfile=config['work_dir'] + "/" + str(self.datea_str) + '_f' + self.fff + '_e700_ens.nc'
        if (not os.path.isfile(outfile) and config['fields'].get('calc_the700hPa','False') == 'True'):

          print("  Computing 700 hPa Theta-E")

          vDict = {'latitude': (lat1, lat2), 'longitude': (lon1, lon2), 'isobaricInhPa': (700, 700),
                   'description': '700 hPa Equivalent Potential Temperature', 'units': 'K', '_FillValue': -9999.}
          vDict = g1.set_var_bounds('temperature', vDict)

          ensmat = g1.create_ens_array('temperature', len(self.atcf_files), vDict)

          for n in range(self.nens):

            tmpk = np.squeeze(g1.read_grib_field('temperature', n, vDict))
            relh = np.squeeze(g1.read_grib_field('relative_humidity', n, vDict))
            relh[:,:] = np.minimum(np.maximum(0.01 * relh[:,:], 0.0001), 1.0)
            relh.attrs['units'] = 'dimensionless'

            tdew = mpcalc.dewpoint_from_relative_humidity(tmpk, relh).to(units.K)

            pres = tmpk.isobaricInhPa.values * units('hPa')

            ensmat[n,:,:] = np.squeeze(mpcalc.equivalent_potential_temperature(pres[None, None], tmpk, tdew))

          ensmat.to_netcdf(outfile, encoding=dencode)

        else:

          print("  Obtaining 700 hPa Theta-e data from " + outfile)

        #  Compute the 500-850 hPa water vapor mixing ratio, or read from file
        outfile=config['work_dir'] + "/" + str(self.datea_str) + '_f' + self.fff + '_q500-850_ens.nc'
        if (not os.path.isfile(outfile) and config['fields'].get('calc_q500-850hPa','False') == 'True'):

          print("  Computing 500-850 hPa Water Vapor")

          vDict = {'latitude': (lat1, lat2), 'longitude': (lon1, lon2),
                   'description': '500-850 hPa Integrated Water Vapor', 'units': 'hPa', '_FillValue': -9999.}
          vDict = g1.set_var_bounds('temperature', vDict)

          ensmat = g1.create_ens_array('temperature', len(self.atcf_files), vDict)

          vDict = {'latitude': (lat1, lat2), 'longitude': (lon1, lon2), 'isobaricInhPa': (500, 850),
                   'description': '500-850 hPa Integrated Water Vapor', 'units': 'hPa', '_FillValue': -9999.}
          vDict = g1.set_var_bounds('temperature', vDict)

          for n in range(self.nens):

            tmpk = np.squeeze(g1.read_grib_field('temperature', n, vDict))
            relh = np.squeeze(g1.read_grib_field('relative_humidity', n, vDict))
            relh[:,:,:] = np.minimum(np.maximum(0.01 * relh[:,:,:], 0.0001), 1.0)
            relh.attrs['units'] = 'dimensionless'

            pres = tmpk.isobaricInhPa.values

            qvap = mpcalc.mixing_ratio_from_relative_humidity(relh, tmpk, pres[:,None,None] * units('hPa'))

            ensmat[n,:,:] = 0.0

            for k in range(len(pres)-1):
              ensmat[n,:,:] = ensmat[n,:,:] + 0.5 * (qvap[k,:,:]+qvap[k+1,:,:]) * abs(pres[k]-pres[k+1])

            ensmat[n,:,:] = ensmat[n,:,:] * (100.0 / mpcon.earth_gravity)

          ensmat.to_netcdf(outfile, encoding=dencode)

        else:

          print("  Obtaining 500-850 hPa water vapor data from " + outfile)


if __name__ == "__main__":
    ComputeTrackFields('2019082900', 12)
