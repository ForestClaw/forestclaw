subroutine heat_init(blockno, mbc,mx,my,meqn, & 
                           xlower,ylower,dx,dy,q)
    IMPLICIT NONE

    INTEGER mbc,mx,my, meqn
    DOUBLE PRECISION xlower,ylower,dx,dy
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

    INTEGER i,j, m
    DOUBLE PRECISION xc,yc, heat_qexact
    INTEGER blockno

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            q(i,j,1) = heat_qexact(xc,yc)
        end do
    end do

end subroutine heat_init


