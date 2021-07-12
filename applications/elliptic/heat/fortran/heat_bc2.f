      subroutine heat_fort_bc2(meqn,mbc,
     &      mx,my,xlower,ylower,dx,dy,q,t,dt,intersects_bc)
      implicit none

      integer meqn, mbc, mx, my, intersects_bc(0:3)
      double precision xlower, ylower, dx, dy, t, dt

      double precision q(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn)

      integer m, i, j, ibc, jbc

c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
c     # zero-order extrapolation:
      if (intersects_bc(0) .ne. 0) then
          do m=1,meqn
              do ibc=1,mbc
                  do j = 1-mbc,my+mbc
                      q(1-ibc,j,m) = q(1,j,m)
                  end do
              end do
          end do
      endif

c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      if (intersects_bc(1) .ne. 0) then
c         # zero-order extrapolation:
          do m=1,meqn
             do ibc=1,mbc
                do j = 1-mbc,my+mbc
                   q(mx+ibc,j,m) = q(mx,j,m)
                end do
            end do
         end do
      endif


c-------------------------------------------------------
c     # bottom boundary:
c-------------------------------------------------------
      if (intersects_bc(2) .ne. 0) then
c         # zero-order extrapolation:
          do m=1,meqn
             do jbc=1,mbc
                do i = 1-mbc,mx+mbc
                   q(i,1-jbc,m) = q(i,1,m)
                end do 
             end do
          end do
      endif


c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------

      if (intersects_bc(3) .ne. 0) then
c         # zero-order extrapolation:
          do m=1,meqn
             do jbc=1,mbc
                do i = 1-mbc,mx+mbc
                   q(i,my+jbc,m) = q(i,my,m)
                end do
            end do
        end do
      end if


      return
      end
