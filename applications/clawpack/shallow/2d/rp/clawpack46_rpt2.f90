subroutine clawpack46_rpt2(ixy,maxm,meqn,mwaves,mbc,mx, & 
           ql,qr,aux1,aux2,aux3,ilr,asdq,bmasdq,bpasdq)
    !! =====================================================
    implicit none

    !! # Riemann solver in the transverse direction for the shallow water
    !! equations .
    !! # Split asdq (= A^* \Delta q, where * = + or -)
    !! # into down-going flux difference bmasdq (= B^- A^* \Delta q)
    !! #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
    !!
    !! # Uses Roe averages and other quantities which were
    !! # computed in rpn2sh and stored in the common block comroe.
    !!

    integer ixy,maxm,meqn,mwaves,mbc,mx,ilr
    double precision, dimension(1-mbc:maxm+mbc, meqn) :: qr, ql, asdq, bmasdq, bpasdq
    double precision, dimension(1-mbc:maxm+mbc,*) :: aux1, aux2, aux3

    double precision grav
    common /cparam/  grav    !# gravitational parameter

    double precision :: waveb(3,3),sb(3)
    double precision, dimension(2-mbc:maxm+mbc) :: u,v,a,h

    integer i,j,mw, m, mu,mv
    double precision a1, a2, a3, hsqrtl, hsqrtr, hsq2


    if (ixy.eq.1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    do i = 2-mbc, mx+mbc
        h(i) = (qr(i-1,1)+ql(i,1))/2.d0
        hsqrtl = dsqrt(qr(i-1,1))
        hsqrtr = dsqrt(ql(i,1))
        hsq2 = hsqrtl + hsqrtr
        u(i) = (qr(i-1,mu)/hsqrtl + ql(i,mu)/hsqrtr) / hsq2
        v(i) = (qr(i-1,mv)/hsqrtl + ql(i,mv)/hsqrtr) / hsq2
        a(i) = dsqrt(grav*h(i))
    enddo

!!      if (-2.gt.1-mbc .or. maxm2 .lt. maxm+mbc) then
!!          write(6,*) 'need to increase maxm2 in rpB'
!!          stop
!!      endif

    do i = 2-mbc, mx+mbc
        a1 = ((v(i)+a(i))*asdq(i,1)-asdq(i,mv))/(2*a(i))
        a2 = asdq(i,mu) - u(i)*asdq(i,1)
        a3 = (-(v(i)-a(i))*asdq(i,1)+asdq(i,mv))/(2*a(i))

        waveb(1,1) = a1
        waveb(mu,1) = a1*u(i)
        waveb(mv,1) = a1*(v(i)-a(i))
        sb(1) = v(i) - a(i)

        waveb(1,2) = 0.0d0
        waveb(mu,2) = a2
        waveb(mv,2) = 0.0d0
        sb(2) = v(i)

        waveb(1,3) = a3
        waveb(mu,3) = a3*u(i)
        waveb(mv,3) = a3*(v(i)+a(i))
        sb(3) = v(i) + a(i)

        !! # compute the flux differences bmasdq and bpasdq

        do m=1,meqn
            bmasdq(i,m) = 0.d0
            bpasdq(i,m) = 0.d0
            do mw=1,mwaves
                bmasdq(i,m) = bmasdq(i,m) + & 
                    dmin1(sb(mw), 0.d0) * waveb(m,mw)
                bpasdq(i,m) = bpasdq(i,m) + & 
                    dmax1(sb(mw), 0.d0) * waveb(m,mw)
            end do
        end do
    end do 

    return
end subroutine clawpack46_rpt2
