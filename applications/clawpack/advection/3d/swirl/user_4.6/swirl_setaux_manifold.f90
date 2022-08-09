subroutine swirl_setaux_manifold(mbc,mx,my,mz, & 
               xlower,ylower,zlower,dx,dy,dz,maux, & 
               aux,blockno, xrot, yrot, zrot, volume, & 
               faceareas)
   implicit none

   integer mbc, mx, my, mz, maux
   integer blockno
   double precision dx,dy, dz, xlower, ylower, zlower
   double precision  aux(1-mbc:mx+mbc,1-mbc:my+mbc, 1-mbc:mz+mbc,maux)

   double precision xrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)
   double precision yrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)
   double precision zrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)

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
   call swirl_set_velocity_manifold(mx,my,mz,mbc, & 
            dx,dy,dz,t,blockno,xlower,ylower,zlower, & 
            xrot, yrot, zrot, faceareas,aux,maux)

   return
end subroutine swirl_setaux_manifold


!!     ==================================================================
subroutine swirl_set_velocity_manifold(mx,my,mz,mbc, & 
           dx,dy,dz,t,blockno,xlower,ylower,zlower, & 
           xrot, yrot, zrot, faceareas,aux, maux)
!!     ==================================================================

    !!
    !!   # set auxiliary arrays
    !!
    !!   # advection
    !!   #    aux(i,j,k,1) is u velocity on left face of cell
    !!   #    aux(i,j,k,2) is v velocity on bottom face of cell
    !!   #    aux(i,j,k,3) is w velocity on back face of cell
    !!

   implicit none

   integer  :: mx, my, mz, mbc, maux, blockno
   double precision :: xlower, ylower, zlower, dx, dy, dz, t
   double precision :: aux(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,maux)
   double precision faceareas(-mbc:mx+mbc+1,-mbc:my+mbc+1, -mbc:mz+mbc+2,3)

   double precision xrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)
   double precision yrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)
   double precision zrot(-mbc:mx+mbc+1,-mbc:my+mbc+1,-mbc:mz+mbc+2,3,3)

   double precision :: xc(1-mbc:mx+mbc)
   double precision :: yc(1-mbc:my+mbc)
   double precision :: zc(1-mbc:mz+mbc)

   integer*8 cont, fclaw_map_get_context

   integer :: i,j,k, mcapa
   double precision :: dx2, dy2, dz2, compute_u, compute_v, compute_w
   double precision :: dxdy, dxdz, dydz, u,v,w,g, xp,yp,zp
   double precision :: nv(3), vn

   integer jj

   cont = fclaw_map_get_context()

   do i = 1-mbc,mx+mbc
      xc(i) = xlower + (i-0.5d0)*dx
   enddo

   do j = 1-mbc,my+mbc
      yc(j) = ylower + (j-0.5d0)*dy
   end do

   do k = 1-mbc,mz+mbc
      zc(k) = zlower + (k-0.5d0)*dz
   end do

   dx2 = 0.5d0*dx
   dy2 = 0.5d0*dy
   dz2 = 0.5d0*dz

   dxdy = dx*dy
   dxdz = dx*dz
   dydz = dy*dz

   mcapa = 1
   do  k = 1-mbc,mz+mbc
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc

            !! Get velocity normal to x-face
            call fclaw3d_map_c2m(cont,blockno,xc(i)-dx2,yc(j),zc(k),xp,yp,zp)
            g = faceareas(i,j,k,1)/dydz
            u = compute_u(xp,yp,zp)
            v = compute_v(xp,yp,zp)
            w = compute_w(xp,yp,zp)

            do jj = 1,3
               nv(jj) = xrot(i,j,k,1,jj)
            end do
            vn = nv(1)*u + nv(2)*v + nv(3)*w
            aux(i,j,k,mcapa + 1) = g*vn

            !! Get velocity normal to y-face
            call fclaw3d_map_c2m(cont,blockno,xc(i),yc(j)-dy2,zc(k),xp,yp,zp)
            g = faceareas(i,j,k,2)/dxdz
            u = compute_u(xp,yp,zp)
            v = compute_v(xp,yp,zp)
            w = compute_w(xp,yp,zp)

            do jj = 1,3
               nv(jj) = yrot(i,j,k,1,jj)
            end do
            vn = nv(1)*u + nv(2)*v + nv(3)*w
            aux(i,j,k,mcapa + 2) = g*vn

            !! Get velocity normal to z-face
            call fclaw3d_map_c2m(cont,blockno,xc(i),yc(j),zc(k)-dz2,xp,yp,zp)
            g = faceareas(i,j,k,3)/dxdy
            u = compute_u(xp,yp,zp)
            v = compute_v(xp,yp,zp)
            w = compute_w(xp,yp,zp)

            do jj = 1,3
               nv(jj) = zrot(i,j,k,1,jj)
            end do
            vn = nv(1)*u + nv(2)*v + nv(3)*w
            aux(i,j,k,mcapa + 3) = g*vn
         enddo
      enddo
   enddo

   return
end subroutine swirl_set_velocity_manifold

