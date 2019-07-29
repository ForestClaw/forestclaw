      subroutine rpt2cons_manifold(ixy,imp,maxm,meqn,mwaves,maux,
     &               mbc,mx,ql,qr,
     &               aux1,aux2,aux3,asdq,bmasdq,bpasdq)
      
      implicit none

      integer ixy, maxm, meqn,mwaves,mbc,mx,maux,imp

      double precision     ql(meqn,1-mbc:maxm+mbc)
      double precision     qr(meqn,1-mbc:maxm+mbc)
      double precision   asdq(meqn,1-mbc:maxm+mbc)
      double precision bmasdq(meqn,1-mbc:maxm+mbc)
      double precision bpasdq(meqn,1-mbc:maxm+mbc)
      double precision   aux1(maux,1-mbc:maxm+mbc)
      double precision   aux2(maux,1-mbc:maxm+mbc)
      double precision   aux3(maux,1-mbc:maxm+mbc)


      integer i, i1, k, idir, iface
      double precision vrrot, vlrot, g, vhat

c     # ixy = 1 --> idir = 1
c     # ixy = 2 --> idir = 0
      idir = 2-ixy

c     # iface = 2*(2-ixy)
      if (ixy .eq. 1) then
         iface = 2  !! Bottom edge
      else
         iface = 0  !! left edge
      endif

      do i = 2-mbc, mx+mbc
c         # imp = 1 : i-1   (for apdq)
c         # imp = 2 : i     (for apdq)
          i1 = i-2+imp   

c         # -----------------------------------------
c         # Lower faces - cell centered velocities
c         # -----------------------------------------
           
c         # 6-7 (left-right)
c         # 8-9 (bottom-top)
          g = aux2(6 + 2*idir,i1)

c         # left-right : 2,3
c         # bottom-top : 4,5         
          vrrot = g*aux2(2 + 2*idir,i1)   !! Left edge of right cell
          vlrot = g*aux1(3 + 2*idir,i1)   !! Right edge of left cell

          vhat = (vrrot + vlrot)/2.0

          bmasdq(1,i) = min(vhat,0.d0)*asdq(1,i)

c         # -----------------------------------------
c         # Upper faces - cell centered velocities
c         # -----------------------------------------

          g = aux3(6+2*idir,i1)

c         # left-right : 2,3
c         # bottom-top : 4,5         
          vrrot = g*aux3(2 + 2*idir,i1)   !! Left edge of right cell
          vlrot = g*aux2(3 + 2*idir,i1)   !! Right edge of left cell

          vhat = (vrrot + vlrot)/2.0

          bpasdq(1,i) = max(vhat,0.d0)*asdq(1,i)

      enddo


      return
      end
