c     # Write out brick vertices for Matlab
      subroutine write_brick_data(n,mi,mj,xv,yv)
      implicit none
      integer n,mi,mj
      double precision xv(n), yv(n)
      integer i

      open(10,file='brick.dat')
      write(10,100) mi,mj
      do i = 1,n
         write(10,110) xv(i),yv(i)
      enddo
  100 format(2I8)
  110 format(2F8.0)

      end
