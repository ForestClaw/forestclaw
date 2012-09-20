      subroutine qinit_mapped(mx,my,meqn,mbc,xp, yp, zp, q_claw,
     &      maux,aux)
      implicit none

      integer mx,my,meqn,mbc, maux
      double precision q_claw(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision    aux(1-mbc:mx+mbc,1-mbc:my+mbc,maux)

      double precision xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)


      end
