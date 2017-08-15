      subroutine  rpn2_cons_update_manifold(meqn,maux, idir,
     &                ql,qr,auxl,auxr,sl,sr,fluxdiff)

      implicit none

      integer meqn,maux,idir, ixy
      double precision ql(meqn), qr(meqn), auxl(maux), auxr(maux)
      double precision fluxdiff(meqn)

      double precision ur,ul,urc,ulc,vrc,vlc, sl,sr

c     # Flux function : ur*qr - ul*ql = apdq + amdq

      urc = auxr(1)
      vrc = auxr(2)

      ulc = auxl(1)    
      vlc = auxl(2)    


      ixy = idir + 1

c     # Assume auxl, auxr are the same (which they most certainly are)
c     # Normals are taken from right cell only (at cell interface)
      if (ixy .eq. 1) then
c        # Project velocity onto normal at the x-face
         ur = auxr(6)*urc + auxr(7)*vrc            
         ul = auxr(6)*ulc + auxr(7)*vlc
      else
c        # Project velocity onto normal at the y-face              
         ur = auxr(8)*urc + auxr(9)*vrc            
         ul = auxr(8)*ulc + auxr(9)*vlc
      endif

      
      fluxdiff(1) = (sr*ur*qr(1) - sl*ul*ql(1))

      end
