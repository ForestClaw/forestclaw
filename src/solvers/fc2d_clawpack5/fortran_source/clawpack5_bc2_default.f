      subroutine clawpack5_bc2_default(meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux,t,dt,mthbc)
      implicit none

      integer meqn, mbc, mx, my, maux, mthbc(4)
      double precision xlower, ylower, dx, dy, t, dt

      double precision q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
      double precision aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)

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
               q(m,1-ibc,j) = q(m,1,j)
            end do
         end do
      end do
      go to 199

  120 continue
c     # periodic:
      write(6,*) 'clawpack5_bc2_default : ', 
     &            'Periodic BCs should not be set in bc2'
      stop
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do m=1,meqn
         do ibc=1,mbc
            do j = 1-mbc, my+mbc
               q(m,1-ibc,j) = q(m,ibc,j)
            end do
         end do
      end do
c     # negate the normal velocity:
      do ibc=1,mbc
         do j = 1-mbc, my+mbc
            q(2,1-ibc,j) = -q(2,ibc,j)
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
            do j = 1-mbc,my+mbc
               q(m,mx+ibc,j) = q(m,mx,j)
            end do
         end do
      end do
      go to 299

  220 continue
c     # periodic:
      write(6,*) 'clawpack5_bc2_default : Periodic ',
     &              'BCs should not be set in bc2'
      stop
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do m=1,meqn
         do ibc=1,mbc
            do j = 1-mbc, my+mbc
               q(m,mx+ibc,j) = q(m,mx+1-ibc,j)
            end do
         end do
      end do
c     # negate the normal velocity:
      do ibc=1,mbc
         do j = 1-mbc, my+mbc
            q(2,mx+ibc,j) = -q(2,mx+1-ibc,j)
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
            do i = 1-mbc,mx+mbc
               q(m,i,1-jbc) = q(m,i,1)
            end do
         end do
      end do
      go to 399

  320 continue
c     # periodic:
      write(6,*) 'clawpack5_bc2_default : ', 
     &           'Periodic BCs should not be set in bc2'
      stop
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do m=1,meqn
         do jbc=1,mbc
            do i = 1-mbc,mx+mbc
               q(m,i,1-jbc) = q(m,i,jbc)
            end do
         end do
      end do
c     # negate the normal velocity:
      do jbc=1,mbc
         do i = 1-mbc, mx+mbc
            q(3,i,1-jbc) = -q(3,i,jbc)
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
      do m=1,meqn
         do jbc=1,mbc
            do i = 1-mbc,mx+mbc
               q(m,i,my+jbc) = q(m,i,my)
            end do
         end do
      end do
      go to 499

  420 continue
c     # periodic:
      write(6,*) 'clawpack5_bc2_default : ', 
     &           'Periodic BCs should not be set in bc2'
      stop
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do m=1,meqn
         do jbc=1,mbc
            do i = 1-mbc, mx+mbc
               q(m,i,my+jbc) = q(m,i,my+1-jbc)
            end do
         end do
      end do
  435       continue
c     # negate the normal velocity:
      do jbc=1,mbc
         do i = 1-mbc, mx+mbc
            q(3,i,my+jbc) = -q(3,i,my+1-jbc)
         end do
      end do
      go to 499

  499 continue

      return
      end
