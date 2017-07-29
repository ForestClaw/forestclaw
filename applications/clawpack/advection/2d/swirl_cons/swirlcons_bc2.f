      subroutine swirlcons_bc2(maxmx,maxmy,meqn,mbc,
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
      go to (110) mthbc(1)+1
c     this is how we skip over this side... if (mthbc(1)+1
c     is not 1,2,3 or 4, then the goto above falls through to here...
      goto 199

  110 continue
c     # zero-order extrapolation:
     
      do 115 j = 1-mbc, my+mbc
         do 115 ibc=1,mbc
c            aux(1-ibc,j,1) = aux(1,j,1)
            do 115 m=1,meqn
               q(1-ibc,j,m) = q(1,j,m)
  115       continue
      go to 199


  199 continue
c
c-------------------------------------------------------
c     # right boundary:
c-------------------------------------------------------
      go to (210) mthbc(2)+1
      goto 299
c
  210 continue
c     # zero-order extrapolation:
      do 215 j = 1-mbc, my+mbc
         do 215 ibc=1,mbc
c            aux(mx+ibc,j,1) = aux(mx,j,1)
            do 215 m=1,meqn
               q(mx+ibc,j,m) = q(mx,j,m)
  215       continue
      go to 299

  299 continue
c
c-------------------------------------------------------
c     # bottom boundary:
c-------------------------------------------------------
      go to (310) mthbc(3)+1
      goto 399
c
  310 continue
c     # zero-order extrapolation:
      do 315 jbc=1,mbc
         do 315 i = 1-mbc, mx+mbc
c            aux(i,1-jbc,2) = aux(i,1,2)
            do 315 m=1,meqn
               q(i,1-jbc,m) = q(i,1,m)
  315       continue
      go to 399


  399 continue
c
c-------------------------------------------------------
c     # top boundary:
c-------------------------------------------------------
      go to (410) mthbc(4)+1
      goto 499
c
  410 continue
c     # zero-order extrapolation:
      do 415 jbc=1,mbc
         do 415 i = 1-mbc, mx+mbc
c            aux(i,my+jbc,2) = aux(i,my,2)
            do 415 m=1,meqn
               q(i,my+jbc,m) = q(i,my,m)
  415       continue
      go to 499

  499 continue

      return
      end
