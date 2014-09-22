      subroutine assign(s,n,x) bind(c,name="assign")
      use iso_c_binding, only : c_ptr
      implicit none

c     # Headers for c routine 'mult'
      include 'c_routines.i'

      type(c_ptr) s
      double precision x(n)
      integer n, i

      do i = 1,n
         x(i) = i
      enddo

      call c_mult(s,n,x)

      end
