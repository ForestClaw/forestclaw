!> @file
!! check to see if value exceeds threshold

!  --------------------------------------------------------------
!> @brief Check if the gradient exceeds the threshold
!!
!! @param[in] blockno the block number
!! @param[in] qval the 
!! @param[in] qmin the minimum q value
!! @param[in] qmax the maximum q value
!! @param[in] quad the value and adjacent values of q
!! @param[in] dx, dy the spacing in the x and y directions
!! @param[in] xc, yc the coordinate of the cell
!! @param[in] threshold the threshold
!! @param[in] init_flag true if in init stage
!! @param[in] is_ghost true if cell is a ghost cell
!! @return 1 if exceeds threshold, 0 if not, -1 if inconclusive.
!  --------------------------------------------------------------
integer function fclaw2d_clawpatch_gradient_exceeds_th(blockno, meqn, & 
                                     qval,qmin,qmax,quad, & 
                                     dx,dy,xc,yc,ivar_threshold, &
                                     threshold, init_flag, is_ghost)
    implicit none
    
    integer :: meqn, ivar_threshold
    double precision :: qval(meqn),qmin(meqn),qmax(meqn),threshold
    double precision :: quad(-1:1,-1:1,meqn)
    double precision :: dx,dy, xc, yc
    integer :: blockno, init_flag
    logical(kind=4) :: is_ghost

    double precision :: dqx, dqy
    double precision :: xp, yp, zp, xpp, ypp, zpp, xpm, ypm, zpm
    double precision, dimension(3) :: t1, t2, t1inv, t2inv
    double precision, dimension(2,2) :: gmat, gmatinv
    double precision :: grad(3), dx2, dy2, d, ds


    integer*8 :: cont, fclaw_map_get_context
    integer :: fclaw_map_is_used

    double precision :: clawpatch_gradient_dot
    integer :: refine
    integer :: m, mq

    mq = ivar_threshold

    if (is_ghost) then
!!      # quad may have uninitialized values.  Test is inconclusive
        fclaw2d_clawpatch_gradient_exceeds_th = -1
        return
    endif

    cont = fclaw_map_get_context()

    dx2 = 2*dx
    dy2 = 2*dy

    dqx = (quad(1,0,mq) - quad(-1,0,mq))/dx2
    dqy = (quad(0,1,mq) - quad(0,-1,mq))/dy2


    refine = 0
    if (fclaw_map_is_used(cont) .ne. 0) THEN
        CALL fclaw_map_2d_c2m(cont,blockno,xc,yc,xp,yp,zp)
        CALL fclaw_map_2d_c2m(cont,blockno,xc+dx,yc,xpp,ypp,zpp)
        CALL fclaw_map_2d_c2m(cont,blockno,xc-dx,yc,xpm,ypm,zpm)
        t1(1) = (xpp - xpm)/dx2
        t1(2) = (ypp - ypm)/dx2
        t1(3) = (zpp - zpm)/dx2

        CALL fclaw_map_2d_c2m(cont,blockno,xc,yc+dy,xpp,ypp,zpp)
        CALL fclaw_map_2d_c2m(cont,blockno,xc,yc-dy,xpm,ypm,zpm)
        t2(1) = (xpp - xpm)/dy2
        t2(2) = (ypp - ypm)/dy2
        t2(3) = (zpp - zpm)/dy2

        gmat(1,1) = clawpatch_gradient_dot(t1,t1)
        gmat(1,2) = clawpatch_gradient_dot(t1,t2)
        gmat(2,1) = gmat(1,2)
        gmat(2,2) = clawpatch_gradient_dot(t2,t2)

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

    ds = sqrt(clawpatch_gradient_dot(grad,grad))

    if (ds .gt. threshold) then
        refine = 1
    endif

    fclaw2d_clawpatch_gradient_exceeds_th = refine

end function fclaw2d_clawpatch_gradient_exceeds_th

!  --------------------------------------------------------------
!> @brief Dot product of two vector
!!
!! @param[in] u, v the vectors
!! @return the dot product
!  --------------------------------------------------------------
double precision function clawpatch_gradient_dot(u,v)
    implicit none

    double precision :: u(3), v(3)

    clawpatch_gradient_dot = u(1)*v(1) + u(2)*v(2) + u(3)*v(3)

end function clawpatch_gradient_dot

