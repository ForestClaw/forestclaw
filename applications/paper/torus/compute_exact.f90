subroutine compute_exact()
    implicit none

    double precision pi, pi2
    common /compi/ pi, pi2

    double precision xc, yc, xp, yp, zp, tfinal
    double precision x0,y0,x1,y1
    integer N, i

    open(10,file='ring.dat')
    read(10,*) N
    read(10,*) tfinal

    open(20,file='ring.out')
    do i = 1,N
        read(10,*) xp, yp, zp

        call map2comp(xp,yp,zp,x0,y0)

        call evolve_point(x0,y0,tfinal,x1,y1)
        write(20,*) x1,y1
    end do
    close(10)
    close(20)

end subroutine compute_exact


subroutine evolve_point(x0,y0,tfinal,x1,y1)
    implicit none

    external qexact_rhs
    external solout

    double precision x0,y0,tfinal,x1,y1

    double precision sigma(3), rtol, atol
    integer itol, iout

    integer Nmax, lwork,nrdens, liwork
    parameter(Nmax=4, nrdens=0)
    parameter(lwork=8*Nmax+5*nrdens+21)
    parameter(liwork=nrdens+21)

    double precision work(lwork), rpar(1)
    integer iwork(liwork), ipar(2), idid

    double precision t0
    integer i


!!  # ------------------------------------------
!!  # Numerical parameters
!!  # ------------------------------------------

    itol = 0
    rtol = 1.d-14
    atol = 1.d-14
    iout = 0

    do i = 1,20
        work(i) = 0
        iwork(i) = 0
    end do

!!  # Evolve from t=t0 to t=tfinal
    t0 = 0

!!  # Initial conditions for ODE

    sigma(1) = x0
    sigma(2) = y0
    sigma(3) = 0

!!  # Needed for backward trace      
    rpar(1) = tfinal

!!  # Trace forwards    
    ipar(1) = 1

!!  # This traces the velocity field back to the origin.

    call dopri5(3,qexact_rhs,t0,sigma,tfinal, &
                 rtol,atol,itol, &
                 solout,iout, work,lwork,iwork,liwork, &
                 rpar,ipar,idid)

    if (idid .ne. 1) then
        write(6,*) 'DOPRI5 : idid .ne. 1'
        stop
    endif

!!  # Initial position in [0,1]x[0,1]

    x1 = sigma(1)
    y1 = sigma(2)

end subroutine evolve_point    






