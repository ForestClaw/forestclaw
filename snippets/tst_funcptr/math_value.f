      subroutine math_value(f,a,b)
      implicit none

      external f
      integer a,b

      call f(a,b)

      end
