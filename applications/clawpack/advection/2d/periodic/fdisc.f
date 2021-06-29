      double precision function fdisc(blockno, x, y)
      implicit none

      double precision x,y
      integer blockno

      fdisc = sqrt(x*x + y*y) - 0.25

      end

