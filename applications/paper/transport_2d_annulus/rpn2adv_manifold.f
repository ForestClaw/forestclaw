c      subroutine rpn2(ixy,maxm,meqn,
c     &      mwaves,mbc,mx,ql,qr,
c     &      auxl,auxr,wave,s,amdq,apdq)
      subroutine rpn2adv_manifold(ixy,maxm,meqn,mwaves,maux,
     &                mbc,mx,ql,qr,
     &                auxl,auxr,wave,s,amdq,apdq)
c
      implicit none

      integer ixy, mx,maxm, meqn,mwaves,mbc, maux

      double precision wave(meqn,mwaves,1-mbc:maxm+mbc)
      double precision    s(mwaves,1-mbc:maxm+mbc)
      double precision   ql(meqn,1-mbc:maxm+mbc)
      double precision   qr(meqn,1-mbc:maxm+mbc)
      double precision apdq(meqn,1-mbc:maxm+mbc)
      double precision amdq(meqn,1-mbc:maxm+mbc)
      double precision auxl(maux,1-mbc:maxm+mbc)
      double precision auxr(maux,1-mbc:maxm+mbc)

      integer meqn1
      parameter(meqn1 = 10)
      double precision delta(meqn1)

      integer i,m

c     # Must use edge velocities
      integer color_equation
      common /eqn_comm/ color_equation

c     # used only for edge velocities
      integer use_stream
      common /velocity_comm/ use_stream

      if (color_equation .ne. 1) then
         write(6,*) 'rpn2adv : Riemann solver used only ',
     &       'for color equation'
         stop
      endif

      if (use_stream .ne. 1) then
         write(6,*) 'rpn2adv : Velocities must be defined using a ',
     &       'stream function'
         stop
      endif

      if (meqn1 .lt. meqn) then
         write(6,*) 'rpn2noncons : meqn1 .lt. meqn'
         stop
      endif

      do i = 2-mbc, mx+mbc
         do m = 1,meqn
            wave(m,1,i) = ql(m,i) - qr(m,i-1)
         enddo

         s(1,i) = auxl(1 + ixy,i)

         do m = 1,meqn
            amdq(m,i) = min(s(1,i), 0.d0) * wave(m,1,i)
            apdq(m,i) = max(s(1,i), 0.d0) * wave(m,1,i)
         enddo
c         write(6,*) s(i,1), wave(i,1,1)         
      enddo
c      write(6,*) ' '

      return
      end
