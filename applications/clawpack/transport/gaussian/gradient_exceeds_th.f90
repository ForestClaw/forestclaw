!! # check to see if value exceeds threshold

logical function gradient_exceeds_th(blockno,& 
                                     qval,qmin,qmax,quad, & 
                                     dx,dy,xc,yc,threshhold)
    implicit none
    
    double precision :: qval,qmin,qmax,threshhold
    double precision :: quad(-1:1,-1:1)
    double precision :: dx,dy, xc, yc
    integer :: blockno
    logical :: refine

    double precision :: dqx, dqy
    double precision :: xp, yp, zp, xpp, ypp, zpp, xpm, ypm, zpm
    double precision, dimension(3) :: t1, t2, t1inv, t2inv
    double precision, dimension(2,2) :: gmat, gmatinv
    double precision :: grad(3), dx2, dy2, d, ds

    integer*8 :: cont, get_context
    logical :: fclaw2d_map_is_used

    double precision :: gradient_dot

    integer :: m

    cont = get_context()

    dx2 = 2*dx
    dy2 = 2*dy

    dqx = (quad(1,0) - quad(-1,0))/dx2
    dqy = (quad(0,1) - quad(0,-1))/dy2

    if (fclaw2d_map_is_used(cont)) THEN
        CALL fclaw2d_map_c2m(cont,blockno,xc,yc,xp,yp,zp)
        CALL fclaw2d_map_c2m(cont,blockno,xc+dx,yc,xpp,ypp,zpp)
        CALL fclaw2d_map_c2m(cont,blockno,xc-dx,yc,xpm,ypm,zpm)
        t1(1) = (xpp - xpm)/dx2
        t1(2) = (ypp - ypm)/dx2
        t1(3) = (zpp - zpm)/dx2

        CALL fclaw2d_map_c2m(cont,blockno,xc,yc+dy,xpp,ypp,zpp)
        CALL fclaw2d_map_c2m(cont,blockno,xc,yc-dy,xpm,ypm,zpm)
        t2(1) = (xpp - xpm)/dy2
        t2(2) = (ypp - ypm)/dy2
        t2(3) = (zpp - zpm)/dy2

        gmat(1,1) = gradient_dot(t1,t1)
        gmat(1,2) = gradient_dot(t1,t2)
        gmat(2,1) = gmat(1,2)
        gmat(2,2) = gradient_dot(t2,t2)

        d = gmat(1,1)*gmat(2,2) - gmat(1,2)*gmat(2,1)

        gmatinv(1,1) = gmat(2,2)/d
        gmatinv(2,2) = gmat(1,1)/d
        gmatinv(1,2) = -gmat(1,2)/d
        gmatinv(2,1) = -gmat(2,1)/d

        do m = 1,3
            t1inv(m) = gmatinv(1,1)*t1(m) + gmatinv(1,2)*t2(m)
            t2inv(m) = gmatinv(2,1)*t1(m) + gmatinv(2,2)*t2(m)
            grad(m) = dqx*t1inv(m) + dqy*t2inv(m)
        end do
    else
        grad(1) = dqx
        grad(2) = dqy
        grad(3) = 0
    endif
    refine = .false.

    ds = sqrt(gradient_dot(grad,grad))

    if (ds .gt. threshhold) then
        refine = .true.
    endif

    gradient_exceeds_th = refine

end function gradient_exceeds_th

double precision function gradient_dot(u,v)
    implicit none

    double precision :: u(3), v(3)

    gradient_dot = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)

end function gradient_dot

