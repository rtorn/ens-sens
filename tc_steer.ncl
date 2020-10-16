;  tc_steer.ncl - NCAR Command Language function to remove the TC from the
;                 3D wind field based on spherical harmonic method of 
;                 Galarneau and Davis (2013).  This capability is not in 
;                 python at this time.
;
;  Inputs:
;
;    usage:  ncl tc_steer.ncl tclat=lat tclon=lon tcradius=radius
;
;    tclat - TC center latitude
;    tclon - TC center longitude
;    tcradius - radius of the TC circulation to remove
;

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

filename = "wind_info.nc"
uwndname = "u"
vwndname = "v"

ncid = addfile(filename, "w")

dNames = getfilevardims(ncid,uwndname)
latname  = dNames(1)
lonname  = dNames(2)

;  Make sure lat is increasing with index (NCL requirement)
fliplat = ncid->$latname$(0) .gt. ncid->$latname$(1)

if ( fliplat ) then
  uwnd = ncid->$uwndname$(:,::-1,:)
  vwnd = ncid->$vwndname$(:,::-1,:)
  lat  = ncid->$latname$(::-1)
  lon  = ncid->$lonname$(:)
else
  uwnd = ncid->$uwndname$(:,:,:)
  vwnd = ncid->$vwndname$(:,:,:)
  lat  = ncid->$latname$(:)
  lon  = ncid->$lonname$(:)
end if

nlat = dimsizes(lat)
nlon = dimsizes(lon)

latcen = uwnd
loncen = uwnd
latarr = uwnd
lonarr = uwnd

do i = 0, nlon-1
do j = 0, nlat-1
  lonarr(:,j,i) = (/ lon(i) /)
  latarr(:,j,i) = (/ lat(j) /)
end do
end do

;  Compute the divergence and vorticity
div = uv2dvF_Wrap(uwnd,vwnd)
vrt = uv2vrF_Wrap(uwnd,vwnd)

if ( tclon .lt. 0.0 ) then
  tclon = tclon + 360.
end if
latcen(:,:,:) = tclat
loncen(:,:,:) = tclon
radius = (/ gc_latlon(latarr,lonarr,latcen,loncen,2,4) /)

;  zero out vorticity and divergence for radii greater than value
newdiv=where(radius.gt.tcradius,0.0,div)
newvrt=where(radius.gt.tcradius,0.0,vrt)

uout = uwnd
vout = vwnd

;  Compute divergence associated with TC, remove from total wind
dv2uvf(newdiv,uout,vout)
uwnd(:,:,:) = uwnd(:,:,:) - uout(:,:,:)
vwnd(:,:,:) = vwnd(:,:,:) - vout(:,:,:)

;  Compute vorticity associated with TC, remove from total wind
vr2uvf(newvrt,uout,vout)
uwnd(:,:,:) = uwnd(:,:,:) - uout(:,:,:)
vwnd(:,:,:) = vwnd(:,:,:) - vout(:,:,:)

;  Write zonal/meridonal wind to file
if ( fliplat ) then
  ncid->$uwndname$(:,::-1,:) = (/ uwnd /)
  ncid->$vwndname$(:,::-1,:) = (/ vwnd /)
else
  ncid->$uwndname$(:,:,:) = (/ uwnd /)
  ncid->$vwndname$(:,:,:) = (/ vwnd /)
end if

end
