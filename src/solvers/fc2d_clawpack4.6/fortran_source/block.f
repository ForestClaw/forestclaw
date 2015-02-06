      subroutine clawpack_set_block(blockno)
      implicit none

      integer blockno, blockno_com
      common /comblock/ blockno_com

      blockno_com = blockno
      end

      integer function clawpack_get_block()
      implicit none

      integer blockno_com
      common /comblock/ blockno_com

      clawpack_get_block = blockno_com
      return
      end
