      subroutine clawpack_bc3_default(maxmx,maxmy,meqn,mbc,
     &      mx,my,xlower,ylower,dx,dy,q,maux,aux,t,dt,mthbc)
      implicit none

      integer maxmx, maxmy, meqn, mbc, mx, my, maux, mthbc(4)
      double precision xlower, ylower, dx, dy, t, dt

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
      do 115 m=1,meqn
         do 115 ibc=1,mbc
            do 115 j = 1-mbc,my+mbc
               q(1-ibc,j,m) = q(1,j,m)
  115       continue
      go to 199

  120 continue
c     # periodic:
      do 125 m=1,meqn
         do 125 ibc=1,mbc
            do 125 j = 1-mbc, my+mbc
               q(1-ibc,j,m) = q(mx+1-ibc,j,m)
  125       continue
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 135 m=1,meqn
         do 135 ibc=1,mbc
            do 135 j = 1-mbc, my+mbc
               q(1-ibc,j,m) = q(ibc,j,m)
  135       continue
c     # negate the normal velocity:
      do 136 ibc=1,mbc
         do 136 j = 1-mbc, my+mbc
            q(1-ibc,j,2) = -q(ibc,j,2)
  136    continue
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
      do 215 m=1,meqn
         do 215 ibc=1,mbc
            do 215 j = 1-mbc,my+mbc
               q(mx+ibc,j,m) = q(mx,j,m)
  215       continue
      go to 299

  220 continue
c     # periodic:
      do 225 m=1,meqn
         do 225 ibc=1,mbc
            do 225 j = 1-mbc, my+mbc
               q(mx+ibc,j,m) = q(ibc,j,m)
  225       continue
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do 235 m=1,meqn
         do 235 ibc=1,mbc
            do 235 j = 1-mbc, my+mbc
               q(mx+ibc,j,m) = q(mx+1-ibc,j,m)
  235       continue
c     # negate the normal velocity:
      do 236 ibc=1,mbc
         do 236 j = 1-mbc, my+mbc
            q(mx+ibc,j,2) = -q(mx+1-ibc,j,2)
  236    continue
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
      do 315 m=1,meqn
         do 315 jbc=1,mbc
            do 315 i = 1-mbc,mx+mbc
               q(i,1-jbc,m) = q(i,1,m)
  315       continue
      go to 399

  320 continue
c     # periodic:
      do 325 m=1,meqn
         do 325 jbc=1,mbc
            do 325 i = 1-mbc, mx+mbc
               q(i,1-jbc,m) = q(i,my+1-jbc,m)
  325       continue
      go to 399

  330 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 335 m=1,meqn
         do 335 jbc=1,mbc
            do 335 i = 1-mbc,mx+mbc
               q(i,1-jbc,m) = q(i,jbc,m)
  335       continue
c     # negate the normal velocity:
      do 336 jbc=1,mbc
         do 336 i = 1-mbc, mx+mbc
            q(i,1-jbc,3) = -q(i,jbc,3)
  336    continue
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
      do 415 m=1,meqn
         do 415 jbc=1,mbc
            do 415 i = 1-mbc,mx+mbc
               q(i,my+jbc,m) = q(i,my,m)
  415       continue
      go to 499

  420 continue
c     # periodic:
      do 425 m=1,meqn
         do 425 jbc=1,mbc
            do 425 i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,jbc,m)
  425       continue
      go to 499

  430 continue
c     # solid wall (assumes 3'rd component is velocity or momentum in y):
      do 435 m=1,meqn
         do 435 jbc=1,mbc
            do 435 i = 1-mbc, mx+mbc
               q(i,my+jbc,m) = q(i,my+1-jbc,m)
  435       continue
c     # negate the normal velocity:
      do 436 jbc=1,mbc
         do 436 i = 1-mbc, mx+mbc
            q(i,my+jbc,3) = -q(i,my+1-jbc,3)
  436    continue
      go to 499

  499 continue

      return
      end
