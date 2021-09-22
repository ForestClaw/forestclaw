      subroutine clawpack46_rpt2cons_manifold(ixy,maxm,meqn,
     &                       mwaves,maux,mbc,mx,
     &                       ql,qr,aux1,aux2,aux3,imp,asdq,
     &                       bmasdq,bpasdq)
      
      implicit none

      integer ixy, maxm, meqn,mwaves,mbc,mx,imp, maux

      double precision     ql(1-mbc:maxm+mbc, meqn)
      double precision     qr(1-mbc:maxm+mbc, meqn)
      double precision   asdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision   aux1(1-mbc:maxm+mbc, maux)
      double precision   aux2(1-mbc:maxm+mbc, maux)
      double precision   aux3(1-mbc:maxm+mbc, maux)


      integer i, i1, idir
      double precision vrrot, vlrot, g, vhat

c     # Transverse direction
      idir = 2 - ixy

      do i = 2-mbc, mx+mbc
          i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq

c         # -----------------------------------------
c         # Lower faces - cell centered velocities
c         # -----------------------------------------
           
c         # 6-7    Edge lengths (x-face, y-face)
          g = aux2(i1,6+idir)

c         # left-right : 2,3
c         # top-bottom : 4,5         
          vrrot = g*aux2(i,  2 + 2*idir)   !! Left edge of right cell
          vlrot = g*aux2(i-1,3 + 2*idir)   !! Right edge of left cell

          vhat = (vrrot + vlrot)/2.0

          bmasdq(i,1) = min(vhat,0.d0)*asdq(i,1)

c         # -----------------------------------------
c         # Upper faces - cell centered velocities
c         # -----------------------------------------

          g = aux3(i1,6+idir)

c         # left-right : 2,3
c         # top-bottom : 4,5         
          vrrot = g*aux3(i,  2 + 2*idir)   !! Left edge of right cell
          vlrot = g*aux3(i-1,3 + 2*idir)   !! Right edge of left cell

          vhat = (vrrot + vlrot)/2.0

          bpasdq(i,1) = max(vhat,0.d0)*asdq(i,1)

      enddo


      return
      end
