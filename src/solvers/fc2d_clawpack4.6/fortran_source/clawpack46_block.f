      subroutine clawpack46_set_block(blockno)
      implicit none

      integer blockno

      integer blockno_com
      common /comblock/ blockno_com

      blockno_com = blockno
      end

      integer function fc2d_clawpack46_get_block()
      implicit none

      integer blockno_com
      common /comblock/ blockno_com

      fc2d_clawpack46_get_block = blockno_com
      return
      end

      subroutine clawpack46_unset_block()
      implicit none

      integer blockno_com
      common /comblock/ blockno_com

      blockno_com = -1
      end
