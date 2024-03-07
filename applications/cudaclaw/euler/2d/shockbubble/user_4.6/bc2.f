      subroutine clawpack46_bc2(maxmx,maxmy,meqn,mbc,
     &      mx,my,xlower,ylower,
     &      dx,dy,q,maux,aux,t,dt,mthbc)
c     =====================================================
c
c     # Standard boundary condition choices for claw2
c
c     # Modified for shock-bubble problem to impose inflow at left edge
c     # Values of inflow density,velocity,energy must be set in setprob.f
c     # and passed in the common block   /cominf/ rinf,vinf,einf
c
c     # At each boundary  k = 1 (left),  2 (right),  3 (top), 4 (bottom):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
c     #                  component of q.
c     ------------------------------------------------
c
c     # Extend the data from the interior cells (1:mx, 1:my)
c     # to the ghost cells outside the region:
c     #   (i, 1-jbc)   for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (i, my+jbc)  for jbc = 1,mbc,  i = 1-mbc, mx+mbc
c     #   (1-ibc, j)   for ibc = 1,mbc,  j = 1-mbc, my+mbc
c     #   (mx+ibc, j)  for ibc = 1,mbc,  j = 1-mbc, my+mbc
c
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx,my,maux
      double precision xlower, ylower, dx, dy, t, dt
      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
      integer mthbc(4)

      double precision rinf, vinf, einf
      common /cominf/ rinf,vinf,einf

      integer ibc, jbc, i, j, m

c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1

      goto 199
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
      do ibc=1,mbc
         do j=1-mbc,my+mbc
            q(1-ibc,j,1) = rinf
            q(1-ibc,j,2) = rinf*vinf
            q(1-ibc,j,3) = 0.d0
            q(1-ibc,j,4) = einf
         end do
      end do
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do m=1,meqn
         do ibc=1,mbc
            do j = 1-mbc, my+mbc
               q(1-ibc,j,m) = q(1,j,m)
            end do
         end do
      end do
      go to 199

  120 continue
c     # periodic:
      do m=1,meqn
         do ibc=1,mbc
            do j = 1-mbc, my+mbc
               q(1-ibc,j,m) = q(mx+1-ibc,j,m)
            end do
         end do
      end do
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do m=1,meqn
         do ibc=1,mbc
            do j = 1-mbc, my+mbc
               q(1-ibc,j,m) = q(ibc,j,m)
            end do
         end do
      end do
c     # negate the normal velocity:
      do ibc=1,mbc
         do j = 1-mbc, my+mbc
            q(1-ibc,j,2) = -q(ibc,j,2)
         end do
      end do
      go to 199

  199 continue
c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1
      goto 299
c
  200 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc2'
      stop
      go to 299

  210 continue
c     # zero-order extrapolation:
      do m=1,meqn
         do ibc=1,mbc
            do j = 1-mbc, my+mbc
               q(mx+ibc,j,m) = q(mx,j,m)
            end do
         end do
      end do
      go to 299

  220 continue
c     # periodic:
      do m=1,meqn
         do ibc=1,mbc
            do j = 1-mbc, my+mbc
               q(mx+ibc,j,m) = q(ibc,j,m)
            end do
         end do
      end do
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do m=1,meqn
         do ibc=1,mbc
            do j = 1-mbc, my+mbc
               q(mx+ibc,j,m) = q(mx+1-ibc,j,m)
            end do
         end do
      end do
c     # negate the normal velocity:
      do ibc=1,mbc
         do j = 1-mbc, my+mbc
            q(mx+ibc,j,2) = -q(mx+1-ibc,j,2)
         end do
      end do
      go to 299

  299 continue
c
c-------------------------------------------------------
c     # bottom boundary:
c-------------------------------------------------------
      go to (300,310,320,330) mthbc(3)+1
      goto 399
c
  300 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc2'
      stop
      go to 399
c
  310 continue
c     # zero-order extrapolation:
      do m=1,meqn
         do jbc=1,mbc
            do i = 1-mbc, mx+mbc
               q(i,1-jbc,m) = q(i,1,m)
            end do
         end do
      end do
      go to 399

  320 continue
c     # periodic:
      do m=1,meqn
         do jbc=1,mbc
            do i = 1-mbc, mx+mbc
               q(i,1-jbc,m) = q(i,my+1-jbc,m)
            end do
         end do
      end do
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do m=1,meqn
         do jbc=1,mbc
            do i = 1-mbc, mx+mbc
               q(i,1-jbc,m) = q(i,jbc,m)
            end do
         end do
      end do
c     # negate the normal velocity:
      do jbc=1,mbc
         do i = 1-mbc, mx+mbc
            q(i,1-jbc,3) = -q(i,jbc,3)
         end do
      end do
      go to 399

  399 continue
c
c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------
      goto (400,410,420,430) mthbc(4)+1
      goto 499
c
  400 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2'
      stop
      go to 499

  410 continue
c     # zero-order extrapolation:
      do m=1,meqn
         do jbc=1,mbc
            do i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,my,m)
            end do
         end do
      end do
      go to 499

  420 continue
c     # periodic:
      do m=1,meqn
         do jbc=1,mbc
            do i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,jbc,m)
            end do
         end do
      end do
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do m=1,meqn
         do jbc=1,mbc
            do i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,my+1-jbc,m)
            end do
         end do
      end do
c     # negate the normal velocity:
      do jbc=1,mbc
         do i = 1-mbc, mx+mbc
            q(i,my+jbc,3) = -q(i,my+1-jbc,3)
         end do
      end do
      go to 499

  499 continue

      return
      end
