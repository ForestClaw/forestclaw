      subroutine clawpack46_set_block(blockno)
      implicit none

      integer blockno, blockno_com
      common /comblock/ blockno_com

      blockno_com = blockno
      end

      integer function clawpack46_get_block()
      implicit none

      integer blockno_com
      common /comblock/ blockno_com

      clawpack46_get_block = blockno_com
      return
      end

      subroutine clawpack46_unset_block()
      implicit none

      integer blockno, blockno_com
      common /comblock/ blockno_com

      blockno_com = -1
      end
