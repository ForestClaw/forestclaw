      subroutine cudaclaw_set_block(blockno)
      implicit none

      integer blockno, blockno_com
      common /comblock/ blockno_com

      blockno_com = blockno
      end

      integer function fc2d_cudaclaw_get_block()
      implicit none

      integer blockno_com
      common /comblock/ blockno_com

      fc2d_cudaclaw_get_block = blockno_com
      return
      end

      subroutine cudaclaw_unset_block()
      implicit none

      integer blockno_com
      common /comblock/ blockno_com

      blockno_com = -1
      end
