      subroutine mapc2m(blockno,xc,yc,xp,yp,zp,alpha)
     &      bind(c,name="mapc2m")
      implicit none
      integer blockno
      double precision xc,yc,xp,yp,zp,alpha

      call mapc2m_fivepatch(blockno,xc,yc,xp,yp,zp,alpha)

      end
