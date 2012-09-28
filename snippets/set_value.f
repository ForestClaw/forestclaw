


      subroutine set_to_2(a)
      implicit none

      integer a

      a = 2
      end

      subroutine set_to_5(a)
      implicit none

      integer a

      a = 5
      end

      subroutine fortran_cannot_square(f,a)
      implicit none

      integer a

c     # fortran doesn't know how to multiply!
      call f(a)

      end
