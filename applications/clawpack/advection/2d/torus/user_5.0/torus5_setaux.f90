
SUBROUTINE torus5_setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
    implicit none

    integer mbc, mx,my, meqn, maux
    double precision dx,dy, xlower, ylower
    double precision aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer i,j, k
    double precision dxdy

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

    if (color_equation .eq. 1) then
!!      # Edge velocities : entries (2-3)      
        call torus5_set_edge_velocities(mx,my,mbc,dx,dy,  &
                                        blockno,xlower,ylower,aux,maux)
    else
        write(6,*) 'torus5_setaux : Transport equation not implemented in version 5'
        stop
    endif

    return
end subroutine torus5_setaux

subroutine torus5_set_edge_velocities(mx,my,mbc,  &
    dx,dy,blockno,xlower,ylower,aux,maux)
    implicit none

    integer mx,my,mbc,maux,blockno
    double precision dx,dy, xlower,ylower

    double precision xc,yc
    double precision xc1, yc1, zc1, xc2, yc2, zc2
    double precision aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    integer*8 cont, get_context

    integer i,j
    double precision vn

    cont = get_context()

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
!!          # x-face - upper vertex
            xc = xlower + (i-1)*dx
            yc = ylower + j*dy

!!          # Map the brick to a unit square      
            call fclaw2d_map_brick2c(cont, blockno,xc,yc,xc1,yc1,zc1)

!!          # x-face - lower vertex
            xc = xlower + (i-1)*dx
            yc = ylower + (j-1)*dy

!!          # Map the brick to a unit square      
            call fclaw2d_map_brick2c(cont, blockno,xc,yc,xc2,yc2,zc2)

            call torus_edge_velocity(xc1,yc1,xc2,yc2,dy,vn)
            aux(2,i,j) = vn
        end do
    end do

    do j = 1-mbc,my+mbc
       do i = 1-mbc,mx+mbc

!!         # Map (xc,yc,blockno) to [0,1]x[0,1]

!!         # y-face - right vertex
           xc = xlower + i*dx
           yc = ylower + (j-1)*dy
           call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)

!!         # y-face - left vertex
           xc = xlower + (i-1)*dx
           yc = ylower + (j-1)*dy
           call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc2,yc2,zc2)

           call torus_edge_velocity(xc1,yc1,xc2,yc2,dx,vn)
           aux(3,i,j) = -vn
       end do
    end do

end subroutine torus5_set_edge_velocities

