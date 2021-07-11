
SUBROUTINE torus5_setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
    implicit none

    integer mbc, mx,my, maux
    double precision dx,dy, xlower, ylower
    double precision aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    INTEGER blockno, fc2d_clawpack5_get_block

    integer color_equation
    common /eqn_comm/ color_equation      

!!     # ----------------------------------------------------------------
!!     # Color equation (edge velocities)
!!     # 1      capacity
!!     # 2-3    Edge velocities
!!     #
!!     # Conservative form (cell-centered velocities)
!!     # 2-5    Cell-centered velocities projected onto four edge normals
!!     # 6-7    Edge lengths (x-face, y-face)
!!     # ----------------------------------------------------------------


    blockno = fc2d_clawpack5_get_block()


    write(6,*) 'torus5_setaux : Transport equation not implemented in version 5'
    stop

    return
end subroutine torus5_setaux

