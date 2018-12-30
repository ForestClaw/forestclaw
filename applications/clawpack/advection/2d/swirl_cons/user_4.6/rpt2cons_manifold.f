      subroutine rpt2cons_manifold(ixy,maxm,meqn,mwaves,mbc,mx,
     &                       ql,qr,aux1,aux2,aux3,imp,asdq,
     &                       bmasdq,bpasdq)
      
      implicit none

      integer ixy, icoor, maxm, meqn,mwaves,mbc,mx,maux,imp

      double precision     ql(1-mbc:maxm+mbc, meqn)
      double precision     qr(1-mbc:maxm+mbc, meqn)
      double precision   asdq(1-mbc:maxm+mbc, meqn)
      double precision bmasdq(1-mbc:maxm+mbc, meqn)
      double precision bpasdq(1-mbc:maxm+mbc, meqn)
      double precision   aux1(1-mbc:maxm+mbc, *)
      double precision   aux2(1-mbc:maxm+mbc, *)
      double precision   aux3(1-mbc:maxm+mbc, *)


      integer iuv, iface,i,j,i1
      double precision vl,vr,vhat,ulc,urc,vrc,vlc,sigma

      iface = 3-ixy
      do i = 2-mbc, mx+mbc
         i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq


c        # Lower faces - cell centered velocities
         ulc = aux1(i1,1)    
         vlc = aux1(i1,2)    

         urc = aux2(i1,1)
         vrc = aux2(i1,2)


         if (ixy .eq. 1) then
c           # Project velocity onto normal at the lower y-face
            vr = aux2(i1,8)*urc + aux2(i1,9)*vrc            
            vl = aux2(i1,8)*ulc + aux2(i1,9)*vlc
            sigma = aux2(i1,5)
         else
c            # Project velocity onto normal at the lower x-face              
            vr = aux2(i1,6)*urc + aux2(i1,7)*vrc            
            vl = aux2(i1,6)*ulc + aux2(i1,7)*vlc
            sigma = aux2(i1,4)
         endif
         vhat = sigma*(vl + vr)/2.d0
         bmasdq(i,1) = min(vhat,0.d0)*asdq(i,1)


c        # Upper faces - cell centered velocities
         ulc = aux2(i1-1,1)    
         urc = aux3(i1,1)

         vlc = aux2(i1-1,2)    
         vrc = aux3(i1,2)

         if (ixy .eq. 1) then
c           # Project velocity onto normal at the upper y-face
            vr = aux3(i1,8)*urc + aux3(i1,9)*vrc            
            vl = aux3(i1,8)*ulc + aux3(i1,9)*vlc
            sigma = aux3(i1,5)
         else
c            # Project velocity onto normal at the upper x-face              
            vr = aux3(i1,6)*urc + aux3(i1,7)*vrc            
            vl = aux3(i1,6)*ulc + aux3(i1,7)*vlc
            sigma = aux3(i1,4)
         endif
         vhat = sigma*(vl + vr)/2.d0
         bpasdq(i,1) = max(vhat,0.d0)*asdq(i,1)
      enddo


      return
      end
