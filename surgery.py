import numpy as np
from windspharm.xarray import VectorWind

def haversine_distance_angle(lat1, lon1, lat2, lon2):
  '''
  Function that computes the distance between two lat/lon pairs.  The result of this function 
  is the distance in kilometers.

  Attributes
      lon1 (float): longitude of first point
      lat1 (float): latitude of first point
      lon2 (float): longitude of second point.  Can be an array
      lat2 (float): latitude of second point.  Can be an array
  '''

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
  '''
  Function to remove the TC from the 3D wind field based on the method of 
  Galarneau and Davis (2013).  The function computes the vorticity and divergence
  within a certain distance of the TC, inverts that wind field (assumed to be the 
  TC-related wind), and subtracts this from the total wind.

  Attributes:
      u_wind (xarray DataArray):  zonal component of the wind (pressure, lat, lon order)
      v_wind (xarray DataArray):  meridional component of the wind (pressure, lat, lon order) 
      tc_location (float):        latitude/longitude of the TC in the model
      radius_km (float):          radius over which to remove TC vorticity and divergence
  '''

  #  Compute the vorticity and divergence fields using spherical harmonics
  ws_wnd = VectorWind(u_wind, v_wind)
  vr = ws_wnd.vorticity()
  dv = ws_wnd.divergence()

  #  Compute the distance from each grid point to the TC center.  Mask all points greater than the radius
  tc_lat, tc_lon = tc_location
  distances = haversine_distance_angle(tc_lat, tc_lon, vr.coords['latitude'], vr.coords['longitude'])
  tc_mask =  distances <= radius_km

  #  Mask vorticity and divergence outside of radius
  dv_masked = dv.where(tc_mask, 0)
  dv_masked.attrs['long_name'] = "tc_horizontal_divergence"
  vr_masked = vr.where(tc_mask, 0)
  vr_masked.attrs['long_name'] = "tc_relative_vorticity"

  #  Invert the vorticity and divergence within the TC radius, 
  # remove from the total wind to get environmental wind
  uout, vout = ws_wnd.getuv(vr_masked, dv_masked)
  u_wind[:,:,:] = u_wind[:,:,:] - uout[:,:,:]
  v_wind[:,:,:] = v_wind[:,:,:] - vout[:,:,:]

  return u_wind, v_wind
