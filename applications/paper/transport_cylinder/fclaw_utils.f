      integer*8 function get_context()
      implicit none

      get_context = 0

      return

      end


      subroutine fclaw2d_map_c2m_basis(cont, 
     &                 x,y, t, tinv, tderivs, flag)
      implicit none

      integer*8 cont
      double precision x,y,t(3,2), tinv(3,2), tderivs(4)
      integer flag

      call cylinder_basis_complete(x,y,t,tinv,tderivs,flag)

      end
