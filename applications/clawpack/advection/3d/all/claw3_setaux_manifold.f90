subroutine claw3_setaux_manifold(mbc,mx,my,mz, & 
        xlower,ylower,zlower,dx,dy,dz,maux, & 
        aux,blockno, xd,yd,zd,xp, yp, zp, volume)
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
            dx,dy,dz,t,blockno,xd,yd,zd,xp, yp, zp, aux,maux)

   return
end subroutine claw3_setaux_manifold


!! This routine relies on a streamfunction and so 
!! assumes that the velocity does not vary in z
subroutine claw3_set_velocity_manifold(mx,my,mz,mbc, & 
           dx,dy,dz,t,blockno,xd,yd,zd,xp, yp, zp, aux,maux)
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

   double precision xd1(3),xd2(3)

   integer i,j, k
   double precision vn
   logical ispillowsphere

   do i = 1-mbc,mx+mbc
      do j = 1-mbc,my+mbc
         do k = 1-mbc,mz+mbc
            !! # Velocity at center of x-face
            !! # although z values are at bottom of cell
            xd1(1) = xd(i,j+1,k)
            xd1(2) = yd(i,j+1,k)
            xd1(3) = zp(i,j,k)  

            xd2(1) = xd(i,j,k)
            xd2(2) = yd(i,j,k)
            xd2(3) = zp(i,j,k)    

            call get_psi_vel(xd1,xd2,dy,vn,t)
            if (ispillowsphere()) then
               if (blockno == 1) then
                  vn = -vn
               endif
            endif
            aux(i,j,k,2) = vn
         end do
      end do
   end do

   do j = 1-mbc,my+mbc
      do i = 1-mbc,mx+mbc
         do k = 1-mbc,mz+mbc
            !!# y-faces
            xd1(1) = xd(i+1,j,k)
            xd1(2) = yd(i+1,j,k)
            xd1(3) = zp(i,j,k)

            xd2(1) = xd(i,j,k)
            xd2(2) = yd(i,j,k)
            xd2(3) = zp(i,j,k)

            call get_psi_vel(xd1,xd2,dx,vn,t)
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
