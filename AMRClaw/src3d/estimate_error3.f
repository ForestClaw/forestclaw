      subroutine estimate_error3(mx,my,mz,mbc,meqn,
     &      xlower, ylower, zlower, dx,dy,dz,t, level,
     &      isBoundary, tol, q, error_measure)
      implicit none

      integer mx,my,mz, meqn, level,mbc
      double precision dx,dy,dz, xlower, ylower, zlower
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc, 1-mbc:mz+mbc,meqn)
      double precision error_measure(1:mx,1:my, 1:mz)

      double precision gx,gy,gz, tol, t, divx, divy, divz, err
      integer isBoundary(6), ip, im, jp, jm, kp, km
      integer imin, imax, jmin, jmax, kmin, kmax
      integer i,j,k,m

      imin = isBoundary(1)
      imax = mx + 1-isBoundary(2)

      jmin = isBoundary(3)
      jmax = my + 1 - isBoundary(4)

      kmin = isBoundary(5)
      kmax = mz + 1 - isBoundary(6)

      do i = 1,mx
         ip = min(i+1,imax)
         im = max(i-1,imin)
         divx = real(ip - im)
         do j = 1,my
            jp = min(j+1,imax)
            jm = max(j-1,imin)
            divy = real(jp - jm)
            do k = 1,mz
               kp = min(k+1,kmax)
               km = max(k-1,kmin)
               divz = real(kp - km)

c              # This is what AMRClaw does
               err = 0.d0
               do m = 1,meqn
                  gx = abs(q(ip,j,k,m) - q(im,j,k,m))
                  gy = abs(q(i,jp,k,m) - q(i,jm,k,m))
                  gz = abs(q(i,j,kp,m) - q(i,j,km,m))
                  err = dmax1(gx,gy,gz,err)
               enddo
               error_measure(i,j,k) = err


c               gx = (q(ip,j,k,1) - q(im,j,k,1))/divx
c               gy = (q(i,jp,k,1) - q(i,jm,k,1))/divy
c               gz = (q(i,j,kp,1) - q(i,j,km,1))/divz
c               error_measure(i,j,k) = sqrt(gx*gx + gy*gy + gz*gz)
            enddo
         enddo
      enddo
      end
