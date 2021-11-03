SUBROUTINE clawpack46_src3(meqn,mbc,mx,my,mz,xlower,ylower,zlower, & 
       dx,dy,dz,q,maux,aux,t,dt)

    !! Called to update q by solving source term equation 
    !! $q_t = \psi(q)$ over time dt starting at time t.
    !!
    !! This default version does nothing. 
 
    IMPLICIT NONE

    INTEGER :: MEQN,mbc,mx,my,mz,maux
    DOUBLE PRECISION :: xlower,ylower,zlower,dx,dy,dz,t,dt
    DOUBLE PRECISION :: aux(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,maux)
    DOUBLE PRECISION :: q(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

eND SUBROUTINE clawpack46_src3
