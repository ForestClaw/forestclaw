SUBROUTINE setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
    IMPLICIT NONE

    INTEGER mbc, mx,my, meqn, maux
    DOUBLE PRECISION dx,dy, xlower, ylower
    DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    INTEGER :: blockno, ghost_only, quadsize
    INTEGER :: i,j, k, m
    DOUBLE PRECISION :: dxdy, xc, yc, nv(3), t, u, v, w
    DOUBLE PRECISION, allocatable, DIMENSION(:) :: quadstore
    DOUBLE PRECISION :: sum, dx0
    INTEGER :: level
    logical :: area_check

    INTEGER maxlevel, rfactor, grid_mx, mi, mj
    COMMON /amr_comm/ maxlevel, rfactor, grid_mx, mi, mj

    INCLUDE "metric_terms.i"

    area_check = .false.

    !! Figure out level of current grid so we can compute areas
    !! in a conservative manner
    dx0 = 1.d0/(mi*grid_mx)
    level = log((dx0+1d-8)/dx)/log(real(rfactor)) + 1

    !!     # ----------------------------------------------------------------
    !!     # Color equation (edge velocities)
    !!     # 1      capacity
    !!     # 2-3    Edge velocities
    !!     #
    !!     # Conservative form (cell-centered velocities)
    !!     # 2-5    Cell-centered velocities projected onto four edge normals
    !!     # 6-7    Edge lengths (x-face, y-face)
    !!     # ----------------------------------------------------------------

    blockno = 0
    ghost_only = 0

    CALL fclaw2d_fort_compute_mesh(mx,my,mbc,&
                                   xlower,ylower,dx,dy,blockno,xp,yp,zp,xd,yd,zd)

    !! from ForestClaw
    m = rfactor**(maxlevel-level)
    ALLOCATE(quadstore(3*(m+1)*(m+1)))
    quadsize = m
    CALL fclaw2d_fort_compute_area(mx,my,mbc,dx,dy, xlower, ylower, blockno, & 
                                   area, quadsize, quadstore, ghost_only)
    DEALLOCATE(quadstore)

    CALL fclaw2d_fort_compute_normals(mx,my,mbc, xp,yp,zp,xd,yd,zd, &
                                      xnormals,ynormals)

    CALL fclaw2d_fort_compute_tangents(mx,my,mbc, xd,yd,zd, xtangents,ytangents, &
                                       edgelengths)

    CALL fclaw2d_fort_compute_surf_normals(mx,my,mbc, xnormals,ynormals, &
                                           edgelengths,curvature, surfnormals,area)

    t = 0
    dxdy = dx*dy

    !!     # Capacity : entry (1)
    DO i = 1-mbc,mx+mbc
        DO j = 1-mbc,my+mbc
            aux(1,i,j) = area(i,j)/dxdy
        END DO
    END DO

    if (area_check) then
        sum = 0
        do i = 1,mx
            do j = 1,my
                sum = sum + area(i,j)
            end do
        end do
        write(6,101) level,sum
        101 format(I5,F24.16) 
    endif

    !! # Center velocities : entries (2-5)
    t = 0
    CALL annulus46_set_center_velocities(mx,my,mbc,dx,dy, blockno,xlower,ylower, &
          edgelengths,xnormals,ynormals,surfnormals, t, aux, maux)

    !! # Needed to scale speeds in Riemann solver when using
    !! # cell-centered velocities
    DO i = 1-mbc,mx+mbc
        DO j = 1-mbc,my+mbc
           !! Need to store all four edge lengths in each cell, since they are all
           !! needed in QAD fix.            
           aux(12,i,j) = edgelengths(i,j,1)/dy
           aux(13,i,j) = edgelengths(i+1,j,1)/dy
           aux(14,i,j) = edgelengths(i,j,2)/dx
           aux(15,i,j) = edgelengths(i,j+1,2)/dx
        ENDDO
    ENDDO

    open(10,file='metric.dat')
    if (mx .eq. 8) then        
        !! coarse grid
        write(10,*) ' '
        write(10,*) '% coarse grid'
        do i = 1,mx
            write(10,105) i,aux(2,i,5),aux(3,i,5)
        end do
    endif
105 format(I5,3F24.16)            
    close(10)

    open(10,file='metric.dat',position='append')
    if (mx .eq. 16) then
        !! fine grid
        write(10,*) ' '
        write(10,*) '% fine grid'
        do i = 1,mx
            write(10,105) i,aux(2,i,0),aux(3,i,0)
        end do
    endif
    close(10)

END SUBROUTINE setaux


SUBROUTINE annulus46_set_center_velocities(mx,my,mbc,dx,dy,blockno, &
     xlower,ylower, edgelengths,xnormals,ynormals,surfnormals, t, aux, maux)
    IMPLICIT NONE

    INTEGER mx,my,mbc,maux,blockno
    DOUBLE PRECISION dx,dy, xlower,ylower
    DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    DOUBLE PRECISION xc, yc, x,y, t
    DOUBLE PRECISION vel(3)

    DOUBLE PRECISION nl(3), nb(3), nt(3), nr(3)

    INTEGER i,j, k

    INCLUDE "metric_terms.i"


    !! # Cell-centered velocities : entries (4,5,6)
    DO i = 1-mbc,mx+mbc
        DO j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            CALL annulus_center_velocity(xc,yc,t, vel)

            aux(2,i,j) = vel(1)
            aux(3,i,j) = vel(2)

            DO k = 1,3
                nl(k) = xnormals(i,j,k)
                nr(k) = xnormals(i+1,j,k)
                nb(k) = ynormals(i,j,k)
                nt(k) = ynormals(i,j+1,k)
            ENDDO

            aux(4,i,j) = nl(1)
            aux(5,i,j) = nl(2)
            aux(6,i,j) = nr(1)
            aux(7,i,j) = nr(2)
            aux(8,i,j) = nb(1)
            aux(9,i,j) = nb(2)
            aux(10,i,j) = nt(1)
            aux(11,i,j) = nt(2)
        END DO
    END DO
END SUBROUTINE annulus46_set_center_velocities
