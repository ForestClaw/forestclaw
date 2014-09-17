      module c_routines
      implicit none
      interface
         subroutine mult(n,x) bind(c,name="mult")
         implicit none
         integer n
         double precision x(n)
         end subroutine
      end interface
      end module


      subroutine assign(n,x) bind(c,name="assign")
      use iso_c_binding
      use c_routines
      implicit none
      double precision x(n)
      integer n, i

      do i = 1,n
           x(i) = i
      enddo

      call mult(n,x)

      end
