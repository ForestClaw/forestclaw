      subroutine print_array(mx,my,mbc,meqn,q)
      implicit none

      integer mx,my,meqn,mbc
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer i,j,m

      do m = 1,meqn
         do i = 1-mbc,mx+mbc
            do j=1-mbc,my+mbc
               write(6,100) i,j,q(1,i,j)
            enddo
            write(6,*) ' '
         enddo
      enddo
  100 format(2I5,E30.16)

      end
