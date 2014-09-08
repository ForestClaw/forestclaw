      subroutine timeinterp_interior(mx,my,mbc,meqn,qcurr,qlast,
     &      qinterp,alpha)
      implicit none

      integer mx,my,mbc,meqn
      double precision alpha
      double precision   qcurr(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision   qlast(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qinterp(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, m

      do m = 1,meqn
         do i = 1,mx
            do j = 1,my
               qinterp(i,j,m) = qlast(i,j,m) +
     &               alpha*(qcurr(i,j,m)-qlast(i,j,m))
            enddo
         enddo
      enddo

      end
