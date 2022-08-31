subroutine clawpack46_qinit(meqn,mbc, mx,my,mz, & 
            xlower,ylower,zlower, dx,dy,dz, q,maux,aux)
    implicit none

    integer :: meqn, mbc, mx, my, mz, maux
    double precision :: xlower, ylower, zlower, dx, dy, dz
    double precision ::   q(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc,meqn)
    double precision :: aux(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc,maux)

    integer example
    common /com_swirl/ example

    integer :: i, j, k, mq
    double precision :: xc, xlow, ylow, zlow,w

    integer :: blockno, fc3d_clawpack46_get_block

    blockno = fc3d_clawpack46_get_block()

    do mq = 1,meqn
        do i = 1-mbc,mx+mbc
            do j = 1-mbc,my+mbc
                do k = 1-mbc, mz + mbc
                    xc = xlower + (i-0.5d0)*dx
                    if (example .eq. 0) then
                        if (xc .lt. 0.5d0) then
                            q(i,j,k,mq) = 1.d0
                        else
                            q(i,j,k,mq) = 0
                        endif
                    else
                        xlow = xlower + (i-1)*dx
                        ylow = ylower + (j-1)*dy
                        zlow = zlower + (k-1)*dz
                        CALL cellave3(blockno,xlow,ylow,zlow,dx,dy,dz,w)
                        q(i,j,k,1) = w
                    endif
                end do
            end do
        end do
    end do

    return
end subroutine clawpack46_qinit
