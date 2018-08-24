      subroutine  rpn2_cons_update_manifold(meqn,maux, idir, q,aux,flux)

      implicit none

      integer meqn,maux,idir
      double precision q(meqn), aux(maux), flux(meqn)
      double precision u,v,urot
      integer iface

c     # Cell-centered velocities         
      u = aux(1)
      v = aux(2)

c     # x-face normal : (6,7)
c     # y-face normal : (8,9)       
      iface = idir+1  
      urot = aux(5+iface)*u + aux(7+iface)*v

c     # xnormals      
      flux(1) = urot*q(1)

      end
