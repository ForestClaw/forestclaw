program phase_iterative
    implicit none

!!  # Grid parameters
    integer maxmx, maxmy
    parameter(maxmx = 2**10, maxmy = 2**10)


!!   # Input parameters
    double precision ax_in, ay_in, dx_in, dy_in, t_in, tstart
    double precision domain_length
    integer mx_in, my_in, nstp, izero, nchar
    character*100 fname
    logical restart

    double precision phi(maxmx, maxmy)
    double precision u(maxmx, maxmy)
    double precision theta, r, rp, amp

!!  # More grid parameters
    integer mx, my, i, j, n, nout, nstep
    double precision ax, ay, bx, by, dx, dy, dt, tfinal, t
    double precision x(maxmx), y(maxmy)

!!  # Misc
    integer iframe, iframe_start
    double precision pi, dsgn

    integer mbdcnd, nbdcnd, ierror, idimf
    double precision lambda, lambda_phi, pertrb
    double precision bda(maxmy), bdb(maxmy)
    double precision bdc(maxmx), bdd(maxmx)
    double precision f(maxmx, maxmy)

!!  # Temporary variables
    double precision u_n, phi_n, g0, g

    double precision phik(maxmx, maxmy)
    double precision uk(maxmx, maxmy)
    double precision phikp1(maxmx, maxmy)
    double precision ukp1(maxmx, maxmy)

    double precision f1(maxmx, maxmy)
    double precision f2(maxmx, maxmy)
    double precision S1(maxmx, maxmy)
    double precision S2(maxmx, maxmy)
    double precision S3(maxmx, maxmy)
    double precision err(2), errmax(2), tol

!!  # Model parameters
    double precision S, alpha, mparm, xi, gamma

    double precision r0, x0, y0, kanio, pdel
    common /parm_comm/ r0, x0, y0, kanio, pdel

    !! Temporary variables
    double precision beta, Tinv

    integer wsize, wsize_t
!!         13M + 4N + M*INT(LOG2(N))
!!             log2(N) <= log2(maxmy) ~ 10
    parameter(wsize = 4*maxmy + (13 + 12)*maxmx)
    double precision work(wsize)

    double precision xlow, ylow, w

    integer k, kmax, ktotal, avg_iterations, method
    logical prt_iterations
    integer jacobi, gauss_seidel
    parameter(jacobi=0, gauss_seidel=1)

    wsize_t = 4*(maxmy) + & 
            (13 + int(log(real(my))/log(2.d0))*mx)
    if (wsize < wsize_t) then
        write(6,*) 'Increase wsize from ', wsize, ' to ',wsize_t
        stop
    endif
    pi = 4.d0*datan(1.d0)


!!  # --------------------------------
!!  # Set up grid
!!  # --------------------------------

    open(10,file='phase_iterative.dat')
    read(10,*) mx, my
!!      read(10,*) ax, bx, ay, by
    read(10,*) domain_length
    read(10,*) nstep, nout
    read(10,*) dt

    if (mx > maxmx .or. my > maxmy) then
        write(6,*) 'mx or my too big'
        stop
    endif

    ax = -domain_length/2
    bx = domain_length/2
    ay = -domain_length/2
    by = domain_length/2


    dx = (bx - ax)/mx
    dy = (by - ay)/my

    do i = 1,mx
        x(i) = ax + (i-0.5)*dx
    enddo

    do j = 1,my
        y(j) = ay + (j-0.5)*dy
    enddo


!!  # ------------------------------------
!!  # Numerical parameters - Jacobi method
!!  # ------------------------------------
    read(10,*) kmax
    read(10,*) tol
    read(10,*) method
    read(10,*) prt_iterations

!!  # ----------------------------
!!  # model parameters
!!  # ----------------------------
    read(10,*) S
    read(10,*) alpha
    read(10,*) mparm
    read(10,*) xi
    read(10,*) kanio
    read(10,*) r0
    read(10,*) pdel
    read(10,*) gamma

!!  # ----------------------------
!!  # Get restart information
!!  # ----------------------------
    read(10,*) restart
    if (restart) then
       read(10,*) fname
    endif
    close(10)

    write(6,*) 'Parameters used for this run : '
    write(6,*) '---------------------------'
    write(6,*) 'Non-dimensional parameters'
    write(6,*) '---------------------------'
    write(6,200) 'S       = ', S
    write(6,200) 'alpha   = ', alpha
    write(6,200) 'm       = ', mparm
    write(6,200) 'xi      = ', xi
    write(6,200) 'k       = ', kanio
    write(6,200) 'r0      = ', r0
    write(6,200) 'pdel    = ', pdel
    write(6,200) 'gamma   = ', gamma
    write(6,*) ' '
    write(6,*) '---------------------------'
    write(6,*) 'Grid parameters'
    write(6,*) '---------------------------'
    write(6,220) 'mx   = ', mx
    write(6,220) 'my   = ', my
    write(6,200) 'ax   = ', ax
    write(6,200) 'bx   = ', bx
    write(6,200) 'ay   = ', ay
    write(6,200) 'by   = ', by
    write(6,200) 'dx   = ', dx
    write(6,200) 'dy   = ', dy
    write(6,220) 'kmax = ',kmax
    write(6,225) 'tol  = ',tol

    if (restart) then
        write(6,*) 'Restart from file ', fname
    else
        write(6,*) 'Initializing run from scratch...'
    endif

 200  format(a,F10.5)
 210  format(a,E12.6)
 225  format(a,E12.1)
 220  format(a,I5)

!!   # --------------------------------
!!   # Initialize phi
!!   # --------------------------------
    write(6,*) 'Initializing phi and u...'

    if (restart) then
        open(20,file=fname)
        read(20,*) mx_in
        read(20,*) my_in
        read(20,*) ax_in
        read(20,*) ay_in
        read(20,*) dx_in
        read(20,*) dy_in
        read(20,*) t_in

        if (mx_in .ne. mx .or. my_in .ne. my) then
            write(6,*) 'Problem with restart: mx or my not right.'
            stop
        elseif (ax_in .ne. ax .or. ay_in .ne. ay) then
            write(6,*) 'Problem with restart: ax or ay not right'
            stop
        endif
        izero = iachar('0')
        tstart = t_in
        iframe_start = 0
        nstp = 1
        do i = 10,7,-1
            nchar = iachar(fname(i:i)) - izero
            iframe_start = iframe_start + nchar*nstp
            nstp = nstp*10
        enddo

        do j = 1,my
            do i = 1,mx
                read(20,*) phi(i,j), u(i,j)
            enddo
        enddo
        close(20)
    else

        x0 = (ax + bx)/2.d0
        y0 = (ay + by)/2.d0
        tstart = 0.0
        iframe = 0
        do j = 1,my
            do i = 1,mx
                xlow = ax + (i-1)*dx
                ylow = ay + (j-1)*dy
!!              rp = sqrt((xlow)**2 + (ylow)**2)
 
                call cellave(xlow,ylow,dx,dy,w)
                u(i,j) = w-1.d0
                phi(i,j) = 1-w

!!              r = sqrt((amp*(x(i) - x0))**2 + (y(j) - y0)**2)
!!              theta = atan2(y(j) - y0,x(i) - x0)
!!              rp = r0*(1.d0 + pdel*cos(kanio*theta))
!!              d = tanh((r-rp)/xi)
!!              u(i,j) = -0.5d0*(d + 1)  !! In [-1,0]
!!              phi(i,j) = 0 0.5d0*(d + 1)  !! In [0,1]
            enddo
        enddo
        write(6,*) 'Done with initialization'

!!      # Output frame only if we have initialized data from scratch.
        write(6,600) 0,tstart,0
        call out2(maxmx,maxmy,mx,my,ax,ay,dx,dy,phi,u,tstart, & 
              iframe_start)
    endif


    !!  % --------------------------------------
    !!  % Time step through results
    !!  % --------------------------------------

    do j = 1,my
        bda(j) = 0.d0
        bdb(j) = 0.d0
    enddo

    do i = 1,mx
        bdc(i) = 0.d0
        bdd(i) = 0.d0
    enddo

    lambda = -1.d0/dt
    idimf = maxmx
    mbdcnd = 3
    nbdcnd = 3

    beta = xi**2/mparm
    Tinv = 1.d0/xi**2
    lambda_phi = lambda*beta*Tinv

    avg_iterations = 0
    iframe = iframe_start
    do n = 1,nout
        t = tstart + n*dt

        !!  write(6,100) n,t
100     format('Step',I5,' at time t = ',1PE12.5)

        !! From right hand side for block linear system
        do j = 1,my
            do i = 1,mx
                phi_n = phi(i,j) !! old phi(i,j)
                u_n = u(i,j)

                g0 = phi_n*(1-phi_n)
                g = g0*g0

                S1(i,j) = 30*g/S
                S2(i,j) = 30*g*xi*alpha*S
                S3(i,j) = g0*(phi_n - 0.5d0)

                f1(i,j) = lambda*(u_n + S1(i,j)*phi_n);
                f2(i,j) = lambda*beta*phi_n - S3(i,j)
            end do
        end do

        !!  # Start Jacobi iterations
        do j = 1,my
            do i = 1,mx
                uk(i,j) = u(i,j)
                phik(i,j) = phi(i,j)
            end do
        end do
        do k = 1,kmax

            !! ----------------------------------------------
            !!  # Set up right hand side for u
            !! ----------------------------------------------
            do j = 1,my
                do i = 1,mx
                    !! this could be  computed outside of k-loop
!!                    phi_n = phi(i,j)
!!                    g0 = phi_n*(1-phi_n)
!!                    g = g0*g0
!!                    S1 = 30*g/S

                    f(i,j) = -lambda*S1(i,j)*phik(i,j) + f1(i,j)
                end do
            end do

            !!  # Solve equation for u; solution will be in f            
            call hstcrt(ax,bx,mx,mbdcnd,bda,bdb,ay,by,my,nbdcnd, & 
                        bdc,bdd,lambda,f,idimf,pertrb,ierror,work)

            if (work(1) .gt. wsize) then
                write(6,*) 'hstcrt : Allocate more work for wsize; ', wsize, work(1)
            endif

            if (ierror .ne. 0) then
                write(6,*) 'hstcrt : ierror /= 0', ierror
                stop
            endif

            if (abs(pertrb) .gt. 0) then
                write(6,*) 'pertrb > 0;  ', pertrb
                stop
            endif

            do j = 1,my
                do i = 1,mx
                    ukp1(i,j) = f(i,j)
                end do
            end do


            !! ----------------------------------------------
            !! Set-up right hand side for phi
            !! ----------------------------------------------
            do j = 1,my
                do i = 1,mx
!!                    u_n = u(i,j)
!!                    phi_n = phi(i,j)
!!                    g0 = phi_n*(1-phi_n)
!!                    g = g0*g0
!!                    S2 = 30*g*xi*alpha*S
                    if (method .eq. jacobi) then
                        f(i,j) = Tinv*(-S2(i,j)*uk(i,j) + f2(i,j))
                    else
                        f(i,j) = Tinv*(-S2(i,j)*ukp1(i,j) + f2(i,j))
                    endif
                end do
            end do

            !!  # Solve equation for phi; solution will be in f
            call hstcrt(ax,bx,mx,mbdcnd,bda,bdb,ay,by,my,nbdcnd, & 
                        bdc,bdd,lambda_phi,f,idimf,pertrb,ierror,work)

            if (ierror .ne. 0) then
                write(6,*) 'hstcrt : ierror /= 0; ', ierror
                stop
            endif

            do j = 1,my
                do i = 1,mx
                    phikp1(i,j) = f(i,j)
                end do
            end do

            errmax(1) = 0
            errmax(2) = 0
            do j = 1,my
                do i = 1,mx
                    err(1) = abs(ukp1(i,j) - uk(i,j))
                    err(2) = abs(phikp1(i,j) - phik(i,j))
                    errmax(1) = max(errmax(1),err(1))
                    errmax(2) = max(errmax(2),err(2))
                end do
            end do
            if (prt_iterations) then
!!                write(6,601) 'Step ', n, 'Iteration ',k,'Residual(u)', errmax(1), & 
!!                                        'Residual(phi)', errmax(2)
            endif 
            ktotal = k
            if (max(errmax(1), errmax(2)) < tol) then
                exit  !! Exit Jacobi iteration
            endif

            !! Swap k and kp1
            do j = 1,my
                do i = 1,mx
                    uk(i,j) = ukp1(i,j)
                    phik(i,j) = phikp1(i,j)
                end do
            end do


        !! end of Jacobi iteration
        end do
        if (prt_iterations) then
            write(6,601) 'Step ', n, 'Iterations/step ',ktotal,'Residual(u)', errmax(1), & 
                     'Residual(phi)', errmax(2)
!!            write(6,*) ' '
        endif
        avg_iterations = avg_iterations + ktotal 

        !!  # ---------------------------------
        !!  # Output results
        !!  # ---------------------------------


        if (nstep*(n/nstep) == n) then
            iframe = iframe + 1
            write(6,600) iframe, t, avg_iterations/nstep
            call out2(maxmx,maxmy,mx,my,ax,ay,dx,dy,phi,u,t,iframe)
            avg_iterations = 0
        endif

600     format('Writing frame ',I5,' at time t = ',1PE12.5,' (avg. iterations = ',I3,')')
601     format(A,I5,A20,I4,A15,E12.4,A15,E12.4)                         

        !! Update time level n solutions
        do j = 1,my
            do i = 1,mx
                u(i,j) = ukp1(i,j)
                phi(i,j) = phikp1(i,j)
            end do
        end do


    !! End of time  loop
    end do 

end program phase_iterative

double precision function fdisc(x,y)
    implicit none

    double precision x,y, r, rp, theta

    double precision r0, x0, y0, kanio, pdel
    common /parm_comm/ r0, x0, y0, kanio, pdel

    r = sqrt((x-x0)**2 + (y-y0)**2)
    theta = atan2(y - y0,x - x0)
    rp = r0*(1.d0 + pdel*cos(kanio*theta))
    fdisc = r - rp

end function fdisc
