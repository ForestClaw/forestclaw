!> Called before each call to step2.
!! Use to set time-dependent aux arrays or perform other tasks.
!!
!! This default version does nothing. 

subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux)

 
    implicit none
    integer, intent(in) :: mbc,mx,my,meqn,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy,t,dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer example
    common /example_comm/ example

    integer i,j
    double precision xc,yc, vel(3)

    if (.not. (example .ge. 2 .and. example .le. 4)) then
        return
    endif

    !! # Cell-centered velocities : entries (4,5,6)
    DO i = 1-mbc,mx+mbc
        DO j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            CALL annulus_center_velocity(xc,yc,t, vel)

            aux(2,i,j) = vel(1)
            aux(3,i,j) = vel(2)
        end do
    end do


end subroutine b4step2
