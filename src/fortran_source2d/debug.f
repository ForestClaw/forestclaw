c     # Debug routines that can be called from either C code
c     # or Fortran.  To call from C/C++, include "clawpack_fort.H"
c     #

      subroutine set_debug_flag(iflag)
      implicit none
      include "debug.i"
      integer iflag
      debug_flag_com = iflag
      call set_debug_on()
      end

      integer function get_debug_flag()
      implicit none
      include "debug.i"
      get_debug_flag = debug_flag_com
      end

      subroutine set_debug_on()
      implicit none
      include "debug.i"
      debug_on_com = .true.
      debug_off_com = .not. debug_on_com
      end

      subroutine set_debug_off()
      implicit none
      include "debug.i"
      debug_off_com = .true.
      debug_on_com = .not. debug_off_com
      end

      logical function debug_is_on()
      implicit none
      include "debug.i"
      debug_is_on = debug_on_com
      end

      logical function debug_is_off()
      implicit none
      include "debug.i"
      debug_is_off = debug_off_com
      end


      subroutine set_debug_info(block_idx, patch_idx, level)
      implicit none
      include "debug.i"
      block_idx_com = block_idx
      patch_idx_com = patch_idx
      level_com = level
      end

      subroutine reset_debug_info()
      implicit none
      include "debug.i"
      block_idx_com = -1
      patch_idx_com = -1
      level_com = -1
      end


      integer function get_block_idx()
      implicit none
      include "debug.i"
      get_block_idx = block_idx_com
      end

      integer function get_patch_idx()
      implicit none
      include "debug.i"
      get_patch_idx = patch_idx_com
      end

      integer function get_level()
      implicit none
      include "debug.i"
      get_level = level_com
      end


      subroutine dump_patch(mx,my,mbc,meqn,mq,q)
      implicit none
      integer mx,my,mbc,meqn,mq
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j

      do j = my+mbc,1-mbc,-1
         write(6,100) j, (q(i,j,mq),i = 1-mbc,mx+mbc)
      enddo
      write(6,*) ' '
  100 format(I5,50F12.4)

      end
