!!     ==================================================================
subroutine clawpack46_setaux(mbc,mx,my,mz,xlower,ylower,zlower, & 
                  dx,dy,dz,maux,aux)
!!     ==================================================================

    !!
    !!   # set auxiliary arrays
    !!
    !!   # advection
    !!   #    aux(i,j,k,1) is u velocity on left face of cell
    !!   #    aux(i,j,k,2) is v velocity on bottom face of cell
    !!   #    aux(i,j,k,3) is w velocity on back face of cell
    !!

    implicit none

    integer  :: mx, my, mz, mbc, maux
    double precision :: xlower, ylower, zlower, dx, dy, dz
    double precision :: aux(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,maux)

    double precision :: xc(1-mbc:mx+mbc)
    double precision :: yc(1-mbc:my+mbc)
    double precision :: zc(1-mbc:mz+mbc)

    integer :: i,j,k
    double precision :: dx2, dy2, dz2, compute_u, compute_v, compute_w

    do i = 1-mbc,mx+mbc
        xc(i) = xlower + (i-0.5d0)*dx
    enddo

    do j = 1-mbc,my+mbc
        yc(j) = ylower + (j-0.5d0)*dy
    end do

    do k = 1-mbc,mz+mbc
        zc(k) = zlower + (k-0.5d0)*dz
    end do

    dx2 = 0.5d0*dx
    dy2 = 0.5d0*dy
    dz2 = 0.5d0*dz

    do  k = 1-mbc,mz+mbc
        do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
                aux(i,j,k,1) = compute_u(xc(i)-dx2, yc(j),zc(k))
                aux(i,j,k,2) = compute_v(xc(i),yc(j)-dy2, zc(k))
                aux(i,j,k,3) = compute_w(xc(i),yc(j), zc(k)-dz2)
            enddo
        enddo
    enddo

    return
end subroutine clawpack46_setaux

double precision function compute_u(x,y,z)
    implicit none
    double precision :: x,y,z,pi, pi2

    common /compi/ pi, pi2

    compute_u = 2.d0 * dsin(pi*x)**2 * dsin(pi2*y) * dsin(pi2*z)
end function compute_u

double precision function compute_v(x,y,z)
    implicit none
    double precision :: x,y,z,pi, pi2

    common /compi/ pi, pi2

    compute_v = -dsin(pi2*x) * dsin(pi*y)**2 * dsin(pi2*z)
end function compute_v

double precision function compute_w(x,y,z)
    implicit none
    double precision :: x,y,z,pi, pi2

    common /compi/ pi, pi2

    compute_w = -dsin(pi2*x) * dsin(pi2*y) * dsin(pi*z)**2
end function compute_w
