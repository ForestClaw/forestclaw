      subroutine assign_array(m,a,val)
      implicit none 

      integer m
      double precision val
      double precision a(m)

      integer i,j

      do i = 1,m
        a(i) = val
      enddo

      end