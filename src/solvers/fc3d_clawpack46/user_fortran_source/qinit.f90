subroutine clawpack46_qinit(meqn,mbc,mx,my,mz,xlower,ylower,zlower,& 
    dx,dy,dz,q,maux,aux)

    !! Set initial conditions for the q array.
    !! This default version prints an error message since it should
    !! not be used directly.  Copy this to an application directory and
    !! loop over all grid cells to set values of q(1:meqn, 1:mx, 1:my, 1:mz).

    implicit none
    
    integer :: mbc,mx,my,mz,maux
    double precision :: xlower,ylower,zlower,dx,dy,dz
    double precision :: aux(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,maux)
    double precision :: q(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    write(6,*) '*** Error -- you must set initial conditions'
    stop

end subroutine clawpack46_qinit
