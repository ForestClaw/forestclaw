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
      double precision nv(2), ul(2), ur(2)

c     # idir = direction of normal solve
c     # ixy = 1 --> idir = 0
c     # ixy = 2 --> idir = 1
      idir = ixy-1

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
           
c         # left/right : 12,13
c         # bottom/top : 14,15
          g = aux2(14-2*idir,i1)

          nv(1) = aux2(8-4*idir,i1)
          nv(2) = aux2(9-4*idir,i1)

          ur(1) = aux2(2,i1)
          ur(2) = aux2(3,i1)

          ul(1) = aux1(2,i1)
          ul(2) = aux1(3,i1)

          vrrot = g*(nv(1)*ur(1) + nv(2)*ur(2))
          vlrot = g*(nv(1)*ul(1) + nv(2)*ul(2))            

          vhat = (vrrot + vlrot)/2.0

          bmasdq(1,i) = min(vhat,0.d0)*asdq(1,i)

c         # -----------------------------------------
c         # Upper faces - cell centered velocities
c         # -----------------------------------------

          g = aux3(14-2*idir,i1)
          nv(1) = aux3(8-4*idir,i1)
          nv(2) = aux3(9-4*idir,i1)
          ur(1) = aux2(2,i1)
          ur(2) = aux2(3,i1)

          ul(1) = aux1(2,i1)
          ul(2) = aux1(3,i1)

          vrrot = g*(nv(1)*ur(1) + nv(2)*ur(2))
          vlrot = g*(nv(1)*ul(1) + nv(2)*ul(2))            

          vhat = (vrrot + vlrot)/2.0

          bpasdq(1,i) = max(vhat,0.d0)*asdq(1,i)

      enddo


      return
      end
