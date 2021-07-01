import types

import numpy as np
import xarray as xr

from windspharm.xarray import VectorWind, _find_latitude_coordinate, _find_longitude_coordinate, _reverse
from windspharm._common import get_apiorder, inspect_gridtype, to3d


def getuv(self, vorticity, divergence):
    """Compute vector winds from vorticity and divergence fields.
        Method to attach to windspharm.standard.vectorWind
        This is kinda hacky but windspharm development seems dead.

    **Argument:**

    *vorticity*
        A scalar field of vorticity. Its shape must be either (nlat, nlon) or
        (nlat, nlon, nfields) where nlat and nlon are the same
        as those for the vector wind components that initialized the
        `VectorWind` instance.

    *divergence*
        A scalar field of divergence. Its shape must be either (nlat, nlon) or
        (nlat, nlon, nfields) where nlat and nlon are the same
        as those for the vector wind components that initialized the
        `VectorWind` instance.

    **Optional argument:**

    **Returns:**

    *u*, *v*
        Zonal and meridional wind components respectively. Their types will 
        match input types to passed to `VectorWind` instance. 
    """
    vrspec = self.s.grdtospec(vorticity)
    dvspec = self.s.grdtospec(divergence)
    ugrd, vgrid = self.s.getuv(vrspec, dvspec)
    return ugrd, vgrid



class ExtendedVectorWind(VectorWind):
    def __init__(self, u, v, rsphere=6.3712e6, legfunc='stored'):
        super().__init__(u, v, rsphere, legfunc)
        self._api.getuv = types.MethodType(getuv, self._api)

    def getuv(self, vorticity, divergence):
         """Compute vector winds from vorticity and divergence fields.

         **Argument:**

         *vorticity*
             A scalar field of vorticity. Its shape must be either (nlat, nlon) or
             (nlat, nlon, nfields) where nlat and nlon are the same
             as those for the vector wind components that initialized the
             `VectorWind` instance.

         *divergence*
             A scalar field of divergence. Its shape must be either (nlat, nlon) or
             (nlat, nlon, nfields) where nlat and nlon are the same
             as those for the vector wind components that initialized the
             `VectorWind` instance.

         **Returns:**

         *u*, *v*
             Zonal and meridional wind components respectively. Their types will 
             match input types to passed to `VectorWind` instance. 
         """
         def clean_array(field):
             if not isinstance(field, xr.DataArray):
                 raise TypeError('scalar field must be an xarray.DataArray')
             name = field.name
             lat, lat_dim = _find_latitude_coordinate(field)
             lon, lon_dim = _find_longitude_coordinate(field)
             if (lat.values[0] < lat.values[1]):
                 # need to reverse latitude dimension
                 field = _reverse(field, lat_dim)
                 lat, lat_dim = _find_latitude_coordinate(field)
             apiorder, _ = get_apiorder(field.ndim, lat_dim, lon_dim)
             apiorder = [field.dims[i] for i in apiorder]
             reorder = field.dims
             field = field.copy().transpose(*apiorder)
             ishape = field.shape
             coords = [field.coords[n] for n in field.dims]
             field = to3d(field.values)
             return field

         vortic = clean_array(vorticity)
         diverg = clean_array(divergence)

         ugrd, vgrd = self._api.getuv(vortic, diverg)
         ugrd = self._metadata(ugrd, 'u',
                     units='m s**-1',
                     standard_name='eastward_wind',
                     long_name='eastward_component_of_wind')
         vgrd = self._metadata(vgrd, 'v',
                     units='m s**-1',
                     standard_name='northward_wind',
                     long_name='northward_component_of_wind')
         return ugrd, vgrd



def haversine_distance_angle(lat1, lon1, lat2, lon2):
    '''Simple implementation of calculating great circle distance and
    angle between points using haversine, hopefully numpy allows use of
    series rather than singular values'''
    R = 6371.  # Earths Circumferences  [km]
    rlat1 = np.radians(lat1)
    rlon1 = np.radians(lon1)
    rlat2 = np.radians(lat2)
    rlon2 = np.radians(lon2)
    distance = np.arccos(
        np.sin(rlat1) * np.sin(rlat2) +
        np.cos(rlat1) * np.cos(rlat2) * np.cos(rlon2 - rlon1)) * R
    return distance

def remove_TC_circulation(u_wind, v_wind, tc_location, radius_km):
    ws_wnd = ExtendedVectorWind(u_wind, v_wind)
    
    vr = ws_wnd.vorticity()
    dv = ws_wnd.divergence()

    tc_lat, tc_lon = tc_location
    distances = haversine_distance_angle(tc_lat, tc_lon, vr.coords['latitude'], vr.coords['longitude'])
    tc_mask =  distances <= radius_km

    dv_masked = dv.where(~tc_mask, 0)
    dv_masked.attrs['long_name'] = "environment_horizontal_divergence"
    vr_masked = vr.where(~tc_mask, 0)
    vr_masked.attrs['long_name'] = "environment_relative_vorticity"
    
    ugrd, vgrd = ws_wnd.getuv(vr_masked, dv_masked)
    return ugrd, vgrd