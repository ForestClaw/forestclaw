      subroutine clawpack46_bc3_default(maxmx,maxmy, maxmz, meqn,mbc,
     &      mx,my,mz, xlower,ylower,zlower, 
     &      dx,dy,dz, q,maux,aux,t,dt,mthbc)
      implicit none

      integer maxmx, maxmy, maxmz, meqn, mbc, mx, my, mz 
      integer maux, mthbc(4)
      double precision xlower, ylower, zlower, dx, dy, dz, t, dt

      double precision q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

      integer m, i, j, ibc, jbc

c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1
c     this is how we skip over this side... if (mthbc(1)+1
c     is not 1,2,3 or 4, then the goto above falls through to here...
      goto 199

  100 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(1)=0 and no BCs specified in bc2'
      stop
      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do m=1,meqn
         do ibc=1,mbc
            do j = 1-mbc,my+mbc
               q(1-ibc,j,m) = q(1,j,m)
            end do
          end do
      end do
      go to 199

  120 continue
c     # periodic:
      write(6,*) 'bc3 : Set periodic BCs using period_x'
      stop

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do  m=1,meqn
          do  ibc=1,mbc
              do  j = 1-mbc, my+mbc
                  q(1-ibc,j,m) = q(ibc,j,m)
              end do
          end do
      end do
c     # negate the normal velocity:
      do  ibc=1,mbc
         do  j = 1-mbc, my+mbc
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
      do  m=1,meqn
         do  ibc=1,mbc
            do  j = 1-mbc,my+mbc
               q(mx+ibc,j,m) = q(mx,j,m)
            end do
         end do
      end do
      go to 299

  220 continue
c     # periodic:
c     # periodic:
      write(6,*) 'bc3 : Set periodic BCs using period_x'
      stop
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do  m=1,meqn
         do  ibc=1,mbc
            do  j = 1-mbc, my+mbc
               q(mx+ibc,j,m) = q(mx+1-ibc,j,m)
             end do
           end do
         end do

c     # negate the normal velocity:
      do  ibc=1,mbc
         do  j = 1-mbc, my+mbc
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
      do  m=1,meqn
         do  jbc=1,mbc
            do  i = 1-mbc,mx+mbc
               q(i,1-jbc,m) = q(i,1,m)
             end do
           end do
         end do
      go to 399

  320 continue
c     # periodic:
c     # periodic:
      write(6,*) 'bc3 : Set periodic BCs using period_y'
      stop
  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do  m=1,meqn
         do  jbc=1,mbc
            do  i = 1-mbc,mx+mbc
               q(i,1-jbc,m) = q(i,jbc,m)
             end do
           end do
         end do
c     # negate the normal velocity:
      do  jbc=1,mbc
         do  i = 1-mbc, mx+mbc
            q(i,1-jbc,3) = -q(i,jbc,3)
          end do
        end do
      go to 399

  399 continue
c
c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------
      go to (400,410,420,430) mthbc(4)+1
      goto 499
c
  400 continue
c     # user-specified boundary conditions go here in place of error output
      write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc2'
      stop
      go to 499

  410 continue
c     # zero-order extrapolation:
      do  m=1,meqn
         do  jbc=1,mbc
            do  i = 1-mbc,mx+mbc
               q(i,my+jbc,m) = q(i,my,m)
             end do
           end do
         end do
      go to 499

  420 continue
c     # periodic:
c     # periodic:
      write(6,*) 'bc3 : Set periodic BCs using period_y'
      stop

      do  m=1,meqn
         do  jbc=1,mbc
            do  i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,jbc,m)
             end do
           end do
         end do
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do  m=1,meqn
         do  jbc=1,mbc
            do  i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,my+1-jbc,m)
             end do
          end do
      end do
c     # negate the normal velocity:
      do  jbc=1,mbc
         do  i = 1-mbc, mx+mbc
            q(i,my+jbc,3) = -q(i,my+1-jbc,3)
          end do
      end do
      
      go to 499

  499 continue

      return
      end
