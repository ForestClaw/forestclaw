      module c_routines
      implicit none
      interface
         subroutine mult(s,n,x) bind(c,name="mult")
         use iso_c_binding, only : c_ptr
         implicit none

         type(c_ptr) s
         integer n
         double precision x(n)
         end subroutine
      end interface
      end module


      subroutine assign(s,n,x) bind(c,name="assign")
      use iso_c_binding, only : c_ptr
      use c_routines
      implicit none

      type(c_ptr) s
      double precision x(n)
      integer n, i

      do i = 1,n
           x(i) = i
      enddo

      call mult(s,n,x)

      end
