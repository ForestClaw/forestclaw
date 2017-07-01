      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &			auxl,auxr,wave,s,amdq,apdq)
      implicit none

      integer maxm, mbc,mwaves,meqn,maux, mx
      integer ixy

      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision amdq(1-mbc:maxm+mbc, meqn)
      double precision apdq(1-mbc:maxm+mbc, meqn)
      double precision auxl(1-mbc:maxm+mbc,*)
      double precision auxr(1-mbc:maxm+mbc,*)


      integer i, iface, m, mq, mw
      double precision ustar, fi,fim
      double precision nv(3), ul(3), ur(3), vl,vr,v
      double precision fstar


      double precision dtcom, dxcom, dycom, tcom
      integer icom, jcom
      integer vflag, get_vflag

      common/comxyt/dtcom,dxcom,dycom,tcom,icom,jcom

      vflag = get_vflag()
      if (vflag .ne. 3) then
         write(6,*) 'rpn3fwave : You must supply cell ',
     &         'centered velocities for this Riemann ',
     &         'solver;  set vflag = 3 in setprob.data'
         write(6,100) 'vflag', vflag
         stop
  100    format(A10,I3)
      endif



      iface = ixy
      do i = 2-mbc, mx+mbc

         do mq = 1,meqn
            do mw = 1,mwaves
               wave(i,mq,mw) = 0
            enddo
         enddo

         do mq = 1,meqn
c           # Get normal (times area) at this edge
            vl = 0
            vr = 0
            do m = 1,3
               nv(m) = auxl(i,1+(iface-1)*3+m)
               ul(m) = auxr(i-1,10+(mq-1)*3 + m)
               ur(m) = auxl(i,10 + (mq-1)*3 + m)

c              # Use the same normal for both left and right
c              # cell centered velocities
               vl = vl + nv(m)*ul(m)
               vr = vr + nv(m)*ur(m)
            enddo

            fi = vr*ql(i,mq)
            fim = vl*qr(i-1,mq)

            ustar = (vl + vr)/2.d0

            wave(i,mq,mq) = fi - fim
            s(i,mq) = ustar
            if (ustar .lt. 0) then
               amdq(i,mq) = wave(i,mq,mq)
               apdq(i,mq) = 0.d0
            else
               amdq(i,mq) = 0.d0
               apdq(i,mq) = wave(i,mq,mq)
            endif
         enddo
      enddo


      return
      end
