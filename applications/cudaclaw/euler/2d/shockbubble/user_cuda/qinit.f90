subroutine cudaclaw_qinit(maxmx,maxmy,meqn,mbc, & 
    mx,my,xlower,ylower,dx,dy,q,maux,aux)
    !! =====================================================
    !!
    !! # Set initial conditions for q.

    IMPLICIT NONE

    integer :: maxmx, maxmy, meqn, mbc, mx, my, maux
    double precision :: xlower, ylower, dx,dy

    double precision :: q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)    
    double precision :: aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, maux)

    double precision :: qin(5), qout(5)
    common /comic/ qin,qout

    double precision :: gamma, gamma1
    common /cparam/  gamma,gamma1

    integer :: idisc
    double precision :: x0, y0, alf, beta, r0
    common/cdisc/ x0,y0,alf,beta,r0,idisc

    double precision :: rinf, vinf, einf
    common /cominf/ rinf,vinf,einf

    integer :: i,j, m, idisc_example
    double precision :: xlow, ylow, win
    ! double precision :: rhoin, rhoout, pout, pin, pinf

    integer :: blockno, fc2d_cudaclaw_get_block

    blockno = fc2d_cudaclaw_get_block()

    if (idisc .ne. 2) then
        write(6,*) 'qinit.f : idisc must be 2 for this example.'
        stop
    endif

    if (meqn .ne. 4) then
        write(6,*) 'qinit.f : meqn must be 4 for this example.'
        stop
    endif

    idisc_example = idisc
    do  i = 1-mbc,mx+mbc
        xlow = xlower + (i-1)*dx
        do j = 1-mbc,my+mbc
            ylow = ylower + (j-1)*dy
            !! # set (xlow,ylow) to lower left corner of grid cell:

            call cellave2(blockno,xlow,ylow,dx,dy,win)
            !! # win is now the fraction of the cell that lies inside the
            !! # circle
            do m = 1,meqn
                q(i,j,m) = win*qin(m) + (1.d0-win)*qout(m)
            enddo
        enddo

        !! # behind shock:
        do  j=1-mbc,my+mbc
            if (xlow .lt. 0.2d0) then
                q(i,j,1) = rinf
                q(i,j,2) = rinf*vinf
                q(i,j,3) = 0.d0
                q(i,j,4) = einf
            endif
        enddo
    enddo

    return
end subroutine cudaclaw_qinit
