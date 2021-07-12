      subroutine bc1(meqn,mbc,mx,xlower,dx,q,maux,aux,t,dt,mthbc)
c     =================================================================
c
c     # Standard boundary condition choices for claw2
c
c     # At each boundary  k = 1 (left),  2 (right):
c     #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
c     #            =  1  for zero-order extrapolation
c     #            =  2  for periodic boundary coniditions
c     #            =  3  for solid walls, assuming this can be implemented
c     #                  by reflecting the data about the boundary and then
c     #                  negating the 2'nd component of q.
c     ------------------------------------------------
c
c     # Extend the data from the computational region
c     #      i = 1, 2, ..., mx2
c     # to the virtual cells outside the region, with
c     #      i = 1-ibc  and   i = mx+ibc   for ibc=1,...,mbc
c
    
      implicit none

      integer meqn, mbc,mx,maux,mthbc(4)

      double precision xlower, dx, t, dt
      double precision q(meqn,1-mbc:mx+mbc)
      double precision aux(maux,1-mbc:mx+mbc)

      integer ibc, m
      double precision eta, b, h, u

c
c
c-------------------------------------------------------
c     # left boundary:
c-------------------------------------------------------
      go to (100,110,120,130) mthbc(1)+1
c
  100 continue
c     # user-specified boundary conditions go here in place of error output
      do ibc = 1,mbc
        eta = q(1,ibc) + aux(1,ibc)
        b = aux(1,1-ibc)
        h = eta - b
        u = q(2,ibc)/q(1,ibc)
        q(1,1-ibc) = h
        q(2,1-ibc) = u*h
      end do

      go to 199
c
  110 continue
c     # zero-order extrapolation:
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,1)
         end do
      end do
      go to 199

  120 continue
c     # periodic:  
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,mx+1-ibc)
         end do
      end do
      go to 199

  130 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do ibc=1,mbc
         do m=1,meqn
            q(m,1-ibc) = q(m,ibc)
         end do
c        # negate the normal velocity:
         q(2,1-ibc) = -q(2,ibc)
      end do
      go to 199

  199 continue

c 
c -------------------------------------------------------
c      # right boundary:
c -------------------------------------------------------
      go to (200,210,220,230) mthbc(2)+1

200   continue
c     # user-specified boundary conditions go here in place of error output
      do ibc = 1,mbc
          eta = q(1,mx-ibc+1) + aux(1,mx-ibc+1)
          b = aux(1,mx+ibc)  !! bathymetry in ghost cell
          h = eta - b
          u = q(2,mx-ibc+1)/q(1,mx-ibc+1)
          q(1,mx+ibc) = h
          q(2,mx+ibc) = u*h
      end do

      go to 299

  210 continue
c     # zero-order extrapolation:
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,mx)
         end do
      end do
      go to 299

  220 continue
c     # periodic:  
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,ibc)
         end do
      end do
      go to 299

  230 continue
c     # solid wall (assumes 2'nd component is velocity or momentum in x):
      do ibc=1,mbc
         do m=1,meqn
            q(m,mx+ibc) = q(m,mx+1-ibc)
         end do
         q(2,mx+ibc) = -q(2,mx+1-ibc)
      end do
      go to 299

  299 continue
c
      return
      end

