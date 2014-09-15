c      subroutine set_block(blockno)
c      implicit none
c
c      integer blockno, blockno_com
c      common /comblock/ blockno_com
c
c      blockno_com = blockno
c      end
c
c      integer function get_block()
c      implicit none
c
c      integer blockno_com
c      common /comblock/ blockno_com
c
c      get_block = blockno_com
c      return
c      end
