subroutine mgtest_fort_rhs(blockno, mbc,mx,my,xlower,ylower,dx,dy,q)
    IMPLICIT NONE

    INTEGER mbc,mx,my
    DOUBLE PRECISION xlower,ylower,dx,dy
    DOUBLE PRECISION q(1-mbc:mx+mbc,1-mbc:my+mbc)

    INTEGER i,j
    DOUBLE PRECISION xc,yc, xc1, yc1, zc1, mgtest_qexact_rhs
    INTEGER blockno

    do i = 1-mbc,mx+mbc
        do j = 1-mbc,my+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            !! This should only be required for mapped domains
            !! call fclaw2d_map_brick2c(cont,blockno,xc,yc,xc1,yc1,zc1)

            q(i,j) =  mgtest_qexact_rhs(xc,yc)
        end do
    end do

end subroutine mgtest_fort_rhs


