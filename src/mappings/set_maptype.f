c     # Set up new mappings here.   The user should supply
c     #   1) a file called 'set_maptype_<newmap>.f', which calls a routine
c     #      below.
c     #   2) a routine set_<newmap>() and is<newmap>()
c     #
c     # The file 'set_maptype_<newmap>.f' will be selectively compiled.  See Makefiles.
c     # See file 'set_maptype_disk.f' and routines set_disk() and isdisk(), below.
c     #
c     # ----------------------
c     # Map numbering :
c     # ----------------------
c     # cart = 1
c     # diamond = 2
c     # disk = 3
c     # hemisphere = 4
c     # sphere = 5
c     # rotated square = 6
c     # bilinear quad = 7
c     # ----------------------

c     # ----------------------------
c     # CART = 1
c     # ----------------------------
      subroutine set_cart()
      implicit none
      call set_map_value(1)
      call set_map_defaults()
      end

      logical function iscart()
      implicit none
      logical get_map_value
      iscart = get_map_value(1)
      end

      logical function isflat_cart()
      implicit none

      isflat_cart = .true.
      end


c     # ----------------------------
c     # DIAMOND = 2
c     # ----------------------------
      subroutine set_diamond()
      implicit none
      call set_map_value(2)
      call set_map_defaults()
      end

      logical function isdiamond()
      implicit none
      logical get_map_value
      isdiamond = get_map_value(2)
      end

      logical function isflat_diamond()
      implicit none

      isflat_diamond = .true.
      end


c     # ----------------------------
c     # DISK = 3
c     # ----------------------------
      subroutine set_disk()
      implicit none
      call set_map_value(3)
      call set_map_defaults()
      end

      logical function isdisk()
      implicit none
      logical get_map_value
      isdisk = get_map_value(3)
      end

      logical function isflat_disk()
      implicit none

      isflat_disk = .true.
      end


c     # ----------------------------
c     # HEMISPHERE = 4
c     # ----------------------------
      subroutine set_hemisphere()
      implicit none
      call set_map_value(4)
      call set_map_defaults()
      end

      logical function ishemisphere()
      implicit none
      logical get_map_value
      ishemisphere = get_map_value(4)
      end

      logical function isflat_hemisphere()
      implicit none

      isflat_hemisphere = .false.
      end


c     # ----------------------------
c     # SPHERE = 5
c     # ----------------------------
      subroutine set_sphere()
      implicit none
      call set_map_value(5)
      call set_map_defaults()
      end

      logical function issphere()
      implicit none
      logical get_map_value
      issphere = get_map_value(5)
      end

      logical function isflat_sphere()
      implicit none

      isflat_sphere = .false.
      end


c     # ----------------------------
c     # ROTSQ = 6
c     # ----------------------------
      subroutine set_rotsq()
      implicit none
      call set_map_value(6)
      call set_map_defaults()
      end

      logical function isrotsq()
      implicit none
      logical get_map_value
      isrotsq = get_map_value(6)
      end

      logical function isflat_rotsq()
      implicit none

      isflat_rotsq = .true.
      end


c     # ----------------------------
c     # BIQUAD = 7
c     # ----------------------------
      subroutine set_biquad()
      implicit none
      call set_map_value(7)
      call set_map_defaults()
      end

      logical function isbiquad()
      implicit none
      logical get_map_value
      isbiquad = get_map_value(7)
      end

      logical function isflat_biquad()
      implicit none

      isflat_biquad = .true.
      end


c     # ------------------------------------
c     # Cartesian BC
c     # ------------------------------------
c
c     # If we have cartesian BCs then we can skip
c     # the tridiagonal solve for boundary conditions.
      subroutine set_cart_bc(cart_bc)
      implicit none

      logical cart_bc, cart_bc_com
      common /comcartbc/ cart_bc_com

      cart_bc_com = cart_bc
      end

      logical function has_cart_bc()
      implicit none

      logical cart_bc_com
      common /comcartbc/ cart_bc_com

      has_cart_bc = cart_bc_com

      end


c     # ------------------------------------------
c     # Mapping routines
c     # ------------------------------------------
      subroutine set_all_maps_false()
      implicit none

      integer max_maps
      parameter(max_maps = 20)
      logical ismap_com(max_maps)
      integer n_maps_com, k

      common /comsurf1/ n_maps_com
      common /comsurf2/ ismap_com

      n_maps_com = 7
      if (n_maps_com .gt. max_maps) then
         write(6,*) 'set_all_maps_false (setup_mesh.f) : ',
     &         'Increase size of max_maps to ', n_maps_com
         stop
      endif
      do k = 1,n_maps_com
         ismap_com(k) = .false.
      enddo

      end

      subroutine set_map_value(imap)
      implicit none

      integer imap

      integer max_maps
      parameter(max_maps = 20)
      logical ismap_com(max_maps)
      integer n_maps_com, k

      common /comsurf1/ n_maps_com
      common /comsurf2/ ismap_com

      call set_all_maps_false()

      if (imap .gt. n_maps_com) then
         write(6,*) 'set_map_value : Too many maps; ',imap
         stop
      endif

      ismap_com(imap) = .true.

      end

      logical function get_map_value(imap)
      implicit none

      integer imap

      integer max_maps
      parameter(max_maps = 20)
      logical ismap_com(max_maps)
      integer n_maps_com, k

      common /comsurf1/ n_maps_com
      common /comsurf2/ ismap_com

      if (imap .gt. n_maps_com) then
         write(6,*) 'get_map_value : Too many maps; ', imap
         stop
      endif

      get_map_value = ismap_com(imap)

      end
