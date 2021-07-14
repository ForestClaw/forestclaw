subroutine clawpack46_qinit(meqn,mbc, mx,my,mz, & 
            xlower,ylower,zlower, dx,dy,dz, q,maux,aux)
    implicit none

    integer :: meqn, mbc, mx, my, mz, maux
    double precision :: xlower, ylower, zlower, dx, dy, dz
    double precision ::   q(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc,meqn)
    double precision :: aux(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc,maux)

    integer :: i, j, k, mq
    double precision :: xc

    do mq = 1,meqn
        do i = 1-mbc,mx+mbc
            do j = 1-mbc,my+mbc
                do k = 1-mbc, mz + mbc
                    xc = xlower + (i-0.5d0)*dx
                    if (xc .lt. 0.5d0) then
                        q(i,j,k,mq) = 1.d0
                    else
                        q(i,j,k,mq) = 0
                    endif
                end do
            end do
        end do
    end do

    return
end subroutine clawpack46_qinit
