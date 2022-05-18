subroutine calc_circ_llgrid(vort, circrad, lat, lon, global, nx, ny, circ)

   implicit none

   real, parameter     :: pid    = atan(1.0) / 45.0
   real, parameter     :: Rearth = 6378.388

   integer, intent(in) :: nx, ny
   logical, intent(in) :: global
   real, intent(in)    :: vort(ny,nx), lat(ny), lon(nx), circrad
   real, intent(out)   :: circ(ny,nx)
   real, allocatable   :: dgrid(:,:)

   integer :: xint, yint, i, ii, iii, j, jj, jjj, k, i1, i2, j1, j2, ngy
   real    :: csum, asum, abox(ny,nx), dlon, dlat, gc_dist

   dlon = abs(lon(2)-lon(1))
   dlat = abs(lat(2)-lat(1))
   yint = ceiling(circrad / (Rearth*pid*dlat))

   !  Compute the area of a latitutde/longitude box
   do i = 1, nx
     do j = 1, ny
       abox(j,i) = abs(dlon*pid*cos(pid*lat(j))*dlat*pid)
     end do
   end do 

   do jj = 1, ny

     xint = min(ceiling(abs(circrad / (dlon*pid*Rearth*cos(pid*lat(jj))))),nx/2)

     j1  = max(jj-yint,1)
     j2  = min(jj+yint,ny)
     ngy = j2-j1+1 

     !  Figure out all of the points that are within radius circrad
     allocate(dgrid(ngy,2*xint+1))
     do i = 1, 2*xint+1
       do j = 1, ngy
         if ( gc_dist(lat(j+j1-1),lon(i),lat(jj),lon(xint+1)) <= circrad ) then
           dgrid(j,i) = 1.0
         else
           dgrid(j,i) = 0.0
         end if
       end do
     end do

     do ii = 1, nx

       csum = 0.0
       asum = 0.0

       if ( global ) then
         i1 = ii - xint
         i2 = ii + xint
       else
         i1 = max(ii-xint, 1)
         i2 = min(ii+xint, nx)
       end if

       !  Loop over all points near a specific center point
       do i = 1, i2-i1+1

         iii = i + i1 - 1
         if ( iii < 1   .and. global )  iii = iii + nx
         if ( iii > nx  .and. global )  iii = iii - nx

         !  Compute the area-average quantity and total area
         do j = 1, ngy
           jjj = j + j1 - 1
           csum = csum + vort(jjj,iii)*abox(jjj,iii)*dgrid(j,i)
           asum = asum + abox(jjj,iii)*dgrid(j,i)
         end do

       end do

       !  Divide by area to get area-average vorticity
       if ( asum > 0.0 ) then
         circ(jj,ii) = csum / asum
       end if

     end do

     deallocate(dgrid)

   end do

   return
end subroutine calc_circ_llgrid


function gc_dist(lat1,lon1,lat2,lon2)

implicit none

real, parameter     :: pid    = atan(1.0) / 45.0
real, parameter     :: Rearth = 6378.388

real, intent(in) :: lat1, lon1, lat2, lon2

real :: gc_dist

gc_dist = Rearth * acos( min(sin(lat1*pid) * sin(lat2*pid) + &
                             cos(lat1*pid) * cos(lat2*pid) * &
                             cos((lon2-lon1)*pid), 1.0) )

return
end function gc_dist
