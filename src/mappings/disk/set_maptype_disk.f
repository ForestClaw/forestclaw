c     # This file gets called from driver routines.
      subroutine set_maptype()

      call set_disk()

      end

      logical function isflat()
      implicit none

      isflat = .true.
      end
