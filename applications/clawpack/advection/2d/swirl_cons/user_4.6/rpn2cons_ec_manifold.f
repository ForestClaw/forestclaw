      subroutine rpn2cons_ec_manifold(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &      auxl,auxr,wave,s,amdq,apdq)
      implicit none

      integer maxm, mbc,mwaves,meqn,mx
      integer ixy

      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision amdq(1-mbc:maxm+mbc, meqn)
      double precision apdq(1-mbc:maxm+mbc, meqn)
      double precision auxl(1-mbc:maxm+mbc,*)
      double precision auxr(1-mbc:maxm+mbc,*)


      integer i, iface
      double precision fstar, uhat, qll,qrr, ul,ur
      double precision urc,ulc,vrc,vlc,sigma

      iface = ixy
      do i = 2-mbc, mx+mbc

c        # Cell-centered velocities         
         urc = auxl(i,1)
         ulc = auxr(i-1,1)    

         vrc = auxl(i,2)
         vlc = auxr(i-1,2)    


c        # Assume auxl, auxr are the same (which they most certainly are)
         if (ixy .eq. 1) then
c           # Project velocity onto normal at the x-face
            ur = auxr(i,6)*urc + auxr(i,7)*vrc            
            ul = auxr(i,6)*ulc + auxr(i,7)*vlc
            sigma = auxr(i,4)
         else
c            # Project velocity onto normal at the y-face              
            ur = auxr(i,8)*urc + auxr(i,9)*vrc            
            ul = auxr(i,8)*ulc + auxr(i,9)*vlc
            sigma = auxr(i,5)
         endif

         qrr = ql(i,1)
         qll = qr(i-1,1)


         uhat = (ur + ul)/2.d0

         if (uhat .ge. 0.d0) then
            fstar = uhat*qll
         else
            fstar = uhat*qrr
         endif
         amdq(i,1) = sigma*(fstar - ul*qll)
         apdq(i,1) = sigma*(ur*qrr - fstar)

         wave(i,1,1) = qrr - qll
         s(i,1) = uhat         
      enddo


      return
      end
