      subroutine clawpack46_set_corners(block_corner_count)
      implicit none

      integer block_corner_count(0:3)
      integer i
      integer fclaw2d_corner_com(0:3)

      common /comcorners/ fclaw2d_corner_com

      do i = 0,3
         fclaw2d_corner_com(i) = block_corner_count(i)
      enddo

      end

      subroutine clawpack46_get_corners(block_corner_count)
      implicit none

      integer block_corner_count(0:3)

      integer i
      integer fclaw2d_corner_com(0:3)

      common /comcorners/ fclaw2d_corner_com

      do i = 0,3
         block_corner_count(i) = fclaw2d_corner_com(i)
      enddo

      end
