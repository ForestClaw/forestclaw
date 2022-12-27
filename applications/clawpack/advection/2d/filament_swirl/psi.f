      double precision function psi(x,y)
      implicit none

      double precision x,y

      double precision filament_psi
      double precision swirl_psi

      integer solver
      common /solver/ solver

      if (solver == 1) then
            psi = filament_psi(x,y)
      else if (solver == 2) then
            psi = swirl_psi(x,y)
      else
            error stop "psi.f: invalid solver"
      end if

      return
      end

      subroutine get_psi_vel(xd1,xd2,ds,vn,t)
      implicit none

      double precision xd1(3),xd2(3), ds, vn, t

      integer solver
      common /solver/ solver

      if (solver == 1) then
            call filament_get_psi_vel(xd1,xd2,ds,vn,t)
      else
            error stop "psi.f: invalid solver"
      end if

      end
