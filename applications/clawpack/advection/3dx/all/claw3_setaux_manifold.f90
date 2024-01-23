subroutine claw3_setaux_manifold(mbc,mx,my,mz, & 
        xlower,ylower,zlower,dx,dy,dz,maux, & 
        aux,blockno, xd,yd,zd,xp, yp, zp, volume, & 
        faceareas)
   implicit none

   integer mbc, mx, my, mz, maux
   integer blockno
   double precision dx,dy, dz, xlower, ylower, zlower
   double precision  aux(1-mbc:mx+mbc,1-mbc:my+mbc, 1-mbc:mz+mbc,maux)

   double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2, -mbc:mx+mbc+2)
   double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2, -mbc:my+mbc+2)
   double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2, -mbc:mz+mbc+2)

   double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
   double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
   double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)

   double precision volume(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc)
   double precision faceareas(-mbc:mx+mbc+1,-mbc:my+mbc+1, -mbc:mz+mbc+2,3)

   integer i,j, k
   double precision dxdydz, t

   dxdydz = dx*dy*dz

   do i = 1-mbc,mx+mbc
      do j = 1-mbc,my+mbc
         do k = 1-mbc,mz+mbc
            aux(i,j,k,1) = volume(i,j,k)/dxdydz
         end do
      end do
   end do

   t = 0
   call claw3_set_velocity_manifold(mx,my,mz,mbc, & 
            dx,dy,dz,t,blockno,xd,yd,zd,xp, yp, zp, aux,maux, & 
            faceareas)

   return
end subroutine claw3_setaux_manifold


!! This routine relies on a streamfunction and so 
!! assumes that the velocity does not vary in z
subroutine claw3_set_velocity_manifold(mx,my,mz,mbc, & 
           dx,dy,dz,t,blockno,xd,yd,zd,xp, yp, zp, aux,maux, & 
           faceareas)
   implicit none

   integer mx,my,mz,mbc,maux, blockno
   double precision dx,dy,dz,t

   double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2, -mbc:mz+mbc+2)
   double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2, -mbc:mz+mbc+2)
   double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2, -mbc:mz+mbc+2)

   double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
   double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)
   double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+1)

   double precision aux(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,maux)
   double precision faceareas(-mbc:mx+mbc+1,-mbc:my+mbc+1, -mbc:mz+mbc+2,3)

   double precision,dimension(0:1,0:1) :: xd1,yd1,zd1

   integer i,j, k,ii, jj,kk
   double precision vn, g
   logical ispillowsphere

   do k = 1-mbc,mz+mbc
      do i = 1-mbc,mx+mbc
         do j = 1-mbc,my+mbc
            !! # Velocity at center of x-face
            !! # although z values are at bottom of cell    
            do jj = 0,1
               do kk = 0,1 
                  xd1(jj,kk) = xd(i,j+jj,k+kk)
                  yd1(jj,kk) = yd(i,j+jj,k+kk)
                  zd1(jj,kk) = zd(i,j+jj,k+kk)  
               end do
            end do

            call get_psi_vel(xd1,yd1,zd1,vn,t)

            !! Scale by gamma : facearea/(dy*dz)
            g = faceareas(i,j,k,1)/(dy*dz)
            !! g = 1.0/dz
            vn = g*vn

            if (ispillowsphere()) then
               if (blockno == 1) then
                  vn = -vn
               endif
            endif
            aux(i,j,k,2) = vn
         end do
      end do
   end do

   do k = 1-mbc,mz+mbc
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            !!# y-faces
            do ii = 0,1
               do kk = 0,1
                  xd1(ii,kk) = xd(i+ii,j,k+kk)
                  yd1(ii,kk) = yd(i+ii,j,k+kk)
                  zd1(ii,kk) = zd(i+ii,j,k+kk)
               end do
            end do

            call get_psi_vel(xd1,yd1,zd1,vn,t)

            !! Scale by gamma : facearea/(dx*dz)
            g = faceareas(i,j,k,2)/(dx*dz)
            !! g = 1.d0/dz

            vn = g*vn

            if (ispillowsphere()) then
               if (blockno == 1) then
                  vn = -vn
               endif
            endif

            aux(i,j,k,3) = -vn
         end do
      end do
   end do


   !! Assume zero velocity in z direction
   do j = 1-mbc,my+mbc
      do i = 1-mbc,mx+mbc
         do k = 1-mbc,mz+mbc
            aux(i,j,k,4) = 0
         end do
      end do
   end do

end subroutine claw3_set_velocity_manifold


!! This takes a 2d psi and averages it in z to get 
!! velocity field for the 3d extruded mesh. 
subroutine get_psi_vel(xd1,yd1,zd1,vn,t)
   implicit none

   double precision, dimension(0:1,0:1) :: xd1,yd1,zd1
   double precision vn,t

   double precision, dimension(0:1,0:1) :: xd1_com,yd1_com,zd1_com
   double precision t_com
   common /com_psi/ xd1_com, yd1_com, zd1_com, t_com

   double precision psi_integrate, p1, p2, p3

   xd1_com = xd1
   yd1_com = yd1
   zd1_com = zd1
   t_com = t

   !! Quadrature points : ri = (0, 0.5, 1.0) in [0,1]
   p1 = psi_integrate(0.d0)
   p2 = psi_integrate(0.5d0)
   p3 = psi_integrate(1.d0)

   !! Use Simpson's rule to integrate in r over [0,1]
   !! Weights are dr/6, 4*dr/6, dr/6, with dr=1
   vn = (p1 + 4*p2 + p3)/6.0

end subroutine get_psi_vel


double precision function psi_integrate(r)

   implicit none
   double precision r

   double precision, dimension(0:1,0:1) :: xd1,yd1,zd1
   double precision t
   common /com_psi/ xd1, yd1, zd1,t

   double precision x0_low, x0_hi, x1_low, x1_hi
   double precision y0_low, y0_hi, y1_low, y1_hi
   double precision z0_low, z0_hi, z1_low, z1_hi
   double precision x0, x1, y0, y1, z0, z1
   double precision psi, u, d

   x0_low = xd1(0,0)
   x0_hi  = xd1(0,1)

   x1_low = xd1(1,0)
   x1_hi  = xd1(1,1)

   y0_low = yd1(0,0)
   y0_hi  = yd1(0,1)

   y1_low = yd1(1,0)
   y1_hi  = yd1(1,1)

   z0_low = zd1(0,0)
   z0_hi  = zd1(0,1)

   z1_low = zd1(1,0)
   z1_hi  = zd1(1,1)

   !! Integrating with respect to r \in [0,1]
   x0 = x0_low + r*(x0_hi - x0_low)
   x1 = x1_low + r*(x1_hi - x1_low)
   y0 = y0_low + r*(y0_hi - y0_low)
   y1 = y1_low + r*(y1_hi - y1_low)
   z0 = z0_low + r*(z0_hi - z0_low)
   z1 = z1_low + r*(z1_hi - z1_low)

   d = sqrt((x1-x0)**2 + (y1-y0)**2 + (z1-z0)**2)

   !! Average velocity along 'ds' direction.
   !! Note : Divide by physical distance between (x0,y0,z0)
   !! (x1,y1,z1) rather than use 2d trick where we only 
   !! need to divide by dx or dy.  
   !! The final integral is scaled by facearea/(ds*dz) which
   !! cancels d here.
   
   u = (psi(x1,y1,z1,t) - psi(x0,y0,z0,t))/d

   psi_integrate = u

end function psi_integrate
