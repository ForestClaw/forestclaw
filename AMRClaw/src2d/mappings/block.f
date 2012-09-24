      subroutine set_block(blockno)
      implicit none

      integer blockno, blockno_com
      common /comblock/ blockno_com

      blockno_com = blockno
      end

      integer function get_block()
      implicit none

      integer blockno_com
      common /comblock/ blockno_com

      get_block = blockno_com
      return
      end
