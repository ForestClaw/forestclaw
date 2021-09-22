subroutine clawpack46_setaux(mbc,mx,my,mz, & 
          xlower,ylower,zlower, dx,dy,dz, maux,aux)
    implicit none

    integer :: mbc, mx, my, mz, maux
    double precision :: xlower, ylower, zlower, dx, dy, dz
    double precision  :: aux(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,maux)

    integer :: i, j, k
    double precision :: xc(1-mbc:mx+mbc)
    double precision :: yc(1-mbc:my+mbc)
    double precision :: xll, yll, psi

    do i = 1-mbc,mx+mbc
        xc(i) = xlower + (i-0.5d0)*dx
    enddo

    do j = 1-mbc,my+mbc
        yc(j) = ylower + (j-0.5d0)*dy
    end do


    do k = 1-mbc,mz+mbc
        do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
                xll = xc(i)
                yll = yc(j)
                !! # difference stream function psi to get normal velocities:
                aux(i,j,k,1) = (psi(xll,     yll+dy) - psi(xll,yll)) / dy
                aux(i,j,k,2) = -(psi(xll+dx, yll)    - psi(xll,yll)) / dx
                aux(i,j,k,3) = 0
            enddo
        enddo
    end do

    return
end subroutine clawpack46_setaux
