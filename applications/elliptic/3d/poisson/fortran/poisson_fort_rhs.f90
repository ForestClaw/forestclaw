subroutine poisson_fort_rhs(blockno, mbc, mx, my, mz, mfields, & 
                           xlower, ylower, zlower, dx, dy, dz, rhs)
    IMPLICIT NONE

    INTEGER mbc, mx, my, mz, mfields
    DOUBLE PRECISION xlower, ylower, zlower, dx, dy, dz
    DOUBLE PRECISION rhs(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc, mfields)
    double precision q

    INTEGER i, j, k, m
    DOUBLE PRECISION xc, yc, zc, poisson_qexact_rhs
    INTEGER blockno

    do i = 1-mbc, mx+mbc
        do j = 1-mbc, my+mbc
            do k = 1-mbc, mz+mbc
                xc = xlower + (i-0.5)*dx
                yc = ylower + (j-0.5)*dy
                zc = zlower + (k-0.5)*dz
                q = poisson_qexact_rhs(xc, yc, zc)
                do m = 1, mfields
                    rhs(i, j, k, m) =  q
                end do
            end do
        end do
    end do

end subroutine poisson_fort_rhs
