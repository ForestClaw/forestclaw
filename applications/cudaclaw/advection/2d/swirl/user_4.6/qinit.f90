subroutine cudaclaw_qinit(maxmx,maxmy, meqn,mbc, & 
    mx,my,xlower,ylower,dx,dy,q,maux,aux)
    implicit none

    integer :: meqn, mbc, mx, my, maux, maxmx, maxmy
    double precision :: xlower, ylower, dx, dy
    double precision :: q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
    double precision :: aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

    integer :: i, j, mq
    double precision :: xc,yc

    do mq = 1,meqn
        do i = 1-mbc,mx+mbc
            xc = xlower + (i-0.5d0)*dx
            do j = 1-mbc,my+mbc
                yc = ylower + (j-0.5d0)*dy
                if (xc .lt. 0.5d0) then
                    q(i,j,mq) = 1
                else
                    q(i,j,mq) = 0
                endif
            enddo
        enddo
    enddo

    return
end
