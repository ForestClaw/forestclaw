c     -*- Fortran -*-

      integer debug_flag
      integer debug_flag_com
      common /comdebug1/ debug_flag_com

      logical debug_on, debug_off
      logical debug_on_com, debug_off_com
      common /comdebug2/ debug_on_com, debug_off_com
c      logical debug_is_on, debug_is_off

      integer block_idx, patch_idx, level
      integer block_idx_com, patch_idx_com, level_com
      common /comdebug3/ block_idx_com, patch_idx_com, level_com
c      integer get_block_idx, get_patch_idx, get_level
