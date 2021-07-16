subroutine clawpack46)setaux(mbc,mx,my,mz,xlower,ylower,zlower, & 
    dx,dy,dz,maux,aux)

    !! Called at start of computation before calling qinit, and
    !! when AMR is used, also called every time a new grid patch is created.
    !! Use to set auxiliary arrays 
    !!   aux(1:maux, 1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc).
    !! Note that ghost cell values may need to be set if the aux arrays
    !! are used by the Riemann solver(s).
    !!
    !! This default version does nothing. 
 
    implicit none

    integer :: mbc,mx,my,mz,maux
    double precision :: xlower,ylower,zlower,dx,dy,dz
    double precision :: aux(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,maux)

end subroutine clawpack46_setaux
