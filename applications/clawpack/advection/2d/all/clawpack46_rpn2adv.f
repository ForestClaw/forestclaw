      subroutine clawpack46_rpn2adv(ixy,maxm,meqn,mwaves,
     &      mbc,mx,ql,qr, auxl,auxr,wave,s,amdq,apdq)
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

      integer i,mq,mw, iface

      if (meqn1 .lt. meqn) then
         write(6,*) 'rpn2noncons : meqn1 .lt. meqn'
         stop
      endif

      iface = ixy
      do i = 2-mbc, mx+mbc
         do mq = 1,meqn
            do mw = 1,mwaves
               wave(i,mq,mw) = 0
            enddo
         enddo

         do mq = 1,meqn
            delta(mq) = ql(i,mq) - qr(i-1,mq)
            do mw = 1,mwaves
               wave(i,mq,mw) = delta(mq)
            end do
         end do

         do mw = 1,mwaves
c           # This assumes that meqn == mwaves            
            s(i,mw) = auxl(i,iface)
         enddo

         do mq = 1,meqn
            amdq(i,mq) = 0
            apdq(i,mq) = 0
            do mw = 1,mwaves
               amdq(i,mq) = amdq(i,mq) + min(s(i,mw), 0.d0) * 
     &                wave(i,mq,mw)
               apdq(i,mq) = apdq(i,mq) + max(s(i,mw), 0.d0) * 
     &                wave(i,mq,mw)
            enddo
         enddo
      enddo

      return
      end
