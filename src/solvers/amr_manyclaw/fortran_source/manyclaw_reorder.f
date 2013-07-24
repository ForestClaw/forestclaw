      subroutine manyclaw_reorder2new(mx,my,mbc,meqn,qin,qout)
      implicit none

      integer mx, my, mbc,meqn
      double precision qin(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision qout(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j,m

      do j = 1-mbc,mx+mbc
         do i = 1-mbc,my+mbc
            do m = 1,meqn
               qout(m,i,j) = qin(i,j,m)
            enddo
         enddo
      enddo

      end


      subroutine manyclaw_reorder2old(mx,my,mbc,meqn,qin,qout)
      implicit none

      integer mx, my, mbc,meqn
      double precision qin(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision qout(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,m

      do j = 1-mbc,mx+mbc
         do i = 1-mbc,my+mbc
            do m = 1,meqn
               qout(i,j,m) = qin(m,i,j)
            enddo
         enddo
      enddo

      end
