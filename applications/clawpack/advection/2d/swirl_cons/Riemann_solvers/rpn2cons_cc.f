      subroutine rpn2cons_cc(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &                       auxl,auxr,wave,s,amdq,apdq)

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
      double precision qstar, qi, qim
      double precision nv(3), ul(3), ur(3), vl,vr,v
      double precision fstar


      double precision dtcom, dxcom, dycom, tcom
      integer icom, jcom, vflag, get_vflag

      common/comxyt/dtcom,dxcom,dycom,tcom,icom,jcom

C       vflag = get_vflag()
C       if (vflag .ne. 3) then
C          write(6,*) 'rpn3cons_cc : You must supply cell ',
C      &         'centered velocities for this Riemann ',
C      &         'solver;  set vflag = 3; vflag = ', vflag
C          stop
C       endif


      iface = ixy
      do i = 2-mbc, mx+mbc

         do mq = 1,meqn
            do mw = 1,meqn
               wave(i,mq,mw) = 0
            enddo
         enddo

         do mq = 1,meqn
c           # Get normal (times area) at this edge
c             vl = 0
c             vr = 0
c             do m = 1,3
c c              # Components of normal to edge
c                nv(m) = auxl(i,1 + (iface-1)*3 + m)
c
c c              # Components of left and right velocities
c                ul(m) = auxr(i-1,10 + (mq-1)*3 + m)
c                ur(m) = auxl(i,  10 + (mq-1)*3 + m)
c
c c              # Use the same normal for both left and right
c c              # cell centered velocities
c                vl = vl + nv(m)*ul(m)
c                vr = vr + nv(m)*ur(m)
c             enddo


c           # In each cell, we have :

c            -----------------
c            |      vt       |
c            |               |
c            |vl   (i,j)  vr |
c            |               |
c            |      vb       |
c            -----------------
c
c           Each velocity is a project of the center velocity onto each of the
c           four faces
c           In the aux array (vl,vr,vb,vt) are stored in postions (2,3,4,5)


c           # At a cell edge, we need vr from the left cell and vl from the
c           right cell.

            vr = auxl(i,iface)
            vl = auxr(i-1,iface)

            qi = ql(i,mq)
            qim = qr(i-1,mq)

            if (vr .ge. 0.d0) then
               if(vl .ge. 0.d0) then
                  if (vr .eq. 0) then
                     qstar = qim
                  else
                     qstar = vl*qim/vr
                  endif
                  wave(i,mq,mq) = qi - qstar
                  s(i,mq) = vr
                  amdq(i,mq) = 0.d0
                  apdq(i,mq) = vr*qi - vl*qim
               else
C                 # Note: Not sure what is wave, s for this case 
                  qstar = 0
                  wave(i,mq,mq) = qi - qstar
                  s(i,mq) = vr
                  amdq(i,mq) = -vl*qim
                  apdq(i,mq) = vr*qi                  
               endif
c               fstar = vl*qim
            else
               if(vl .le. 0.d0) then
                  if (vl .eq. 0) then
                     qstar = qi
                  else
                     qstar = vr*qi/vl
                  endif
                  wave(i,mq,mq) = qstar - qim
                  s(i,mq) = vl
                  amdq(i,mq) = vr*qi - vl*qim
                  apdq(i,mq) = 0.d0
               else
                  write(*,*) "vl > 0, vr < 0"
                  qstar = vr*qi/vl
                  wave(i,mq,mq) = qstar - qim
                  s(i,mq) = vl
C                   amdq(i,mq) = 0.5*(vr*qi - vl*qim)
C                   apdq(i,mq) = 0.5*(vr*qi - vl*qim)
                  amdq(i,mq) = vr*qi - vl*qim
                  apdq(i,mq) = 0.d0
               endif
c               fstar = vr*qi
            endif

c           # If we define apdq/amdq this way, we have essentially
c           # the edge-centered scheme, but using upwinded velocites
c           # from cell centers.
c           amdq(i,1) = fstar
c           apdq(i,1) = -fstar
         enddo
      enddo

      return
      end
