
subroutine clawpack46_b4step3(mbc,mx,my,mz,meqn,q,xlower,ylower,zlower, &
    dx,dy,dz,t,dt,maux,aux)

    ! Called before each call to step3.
    ! Use to set time-dependent aux arrays or perform other tasks.
    !
    ! This default version does nothing. 
 
    implicit none
    integer  :: mbc,mx,my,mz,meqn,maux
    double precision  :: xlower,ylower,zlower,dx,dy,dz,t,dt
    double precision  :: q(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision  :: aux(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,maux)

end subroutine clawpack46_b4step3
