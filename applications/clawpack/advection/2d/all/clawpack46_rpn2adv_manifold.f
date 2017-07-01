      subroutine clawpack46_rpn2adv_manifold(ixy,maxm,meqn,
     &      mwaves,mbc,mx,ql,qr,
     &      auxl,auxr,wave,s,amdq,apdq)
c
      implicit none

      integer ixy, mx,maxm, meqn,mwaves,mbc

      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision apdq(1-mbc:maxm+mbc, meqn)
      double precision amdq(1-mbc:maxm+mbc, meqn)
      double precision auxl(1-mbc:maxm+mbc, *)
      double precision auxr(1-mbc:maxm+mbc, *)

      integer meqn1
      parameter(meqn1 = 10)
      double precision delta(meqn1)

      integer i,m,mw, iface, m1, m2, get_vflag, vflag

      if (meqn1 .lt. meqn) then
         write(6,*) 'rpn2noncons : meqn1 .lt. meqn'
         stop
      endif

      iface = ixy
      do i = 2-mbc, mx+mbc
         do m1 = 1,meqn
            do m2 = 1,meqn
               wave(i,m1,m2) = 0
            enddo
         enddo

         do m = 1,meqn
            delta(m) = ql(i,m) - qr(i-1,m)
            wave(i,m,m) = delta(m)
            s(i,m) = auxl(i,1 + iface)
         enddo

         do m = 1,meqn
            amdq(i,m) = 0
            apdq(i,m) = 0
            do mw = 1,mwaves
               amdq(i,m) = amdq(i,m) + min(s(i,mw), 0.d0) * wave(i,m,mw)
               apdq(i,m) = apdq(i,m) + max(s(i,mw), 0.d0) * wave(i,m,mw)
            enddo
         enddo
      enddo

      return
      end
