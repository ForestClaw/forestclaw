      module c_routines
      implicit none
      interface
         subroutine mult(s,n,x) bind(c,name="mult")
         use iso_c_binding, only : c_ptr, c_double
         implicit none

         type(c_ptr) s
         integer n
         real(c_double) x(n)   !! double precision works too
         end subroutine
      end interface
      end module


      subroutine assign(s,n,x) bind(c,name="assign")
      use iso_c_binding, only : c_ptr, c_double
      use c_routines, only : mult
      implicit none

      type(c_ptr) s
      real(c_double) x(n)
      integer n, i

      do i = 1,n
           x(i) = i
      enddo

      call mult(s,n,x)

      end
