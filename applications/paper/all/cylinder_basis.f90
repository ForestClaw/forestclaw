SUBROUTINE cylinder_basis_complete(x,y, t, tinv,tderivs, flag)
    IMPLICIT NONE
      
    DOUBLE PRECISION x,y 
    INTEGER flag
    DOUBLE PRECISION t(3,2), tinv(3,2), tderivs(3,2,2)

    DOUBLE PRECISION pi, pi2
    COMMON /compi/ pi, pi2

    DOUBLE PRECISION r_cyl, h_cyl
    COMMON /cylinder_comm/ r_cyl, h_cyl

    DOUBLE PRECISION r1, r1x, r1xx, r1y, r1yy, r1xy
    DOUBLE PRECISION R,  Rx,  Rxx,  Ry,  Ryy,  Rxy
    DOUBLE PRECISION f(3), fx(3), fy(3), fxx(3), fyy(3), fxy(3)
    DOUBLE PRECISION g(3), gx(3), gy(3), gxx(3), gyy(3), gxy(3)
    DOUBLE PRECISION pi4
    DOUBLE PRECISION map_dot

    DOUBLE PRECISION theta, z
    DOUBLE PRECISION thetax, thetay, zx, zy

    INTEGER k, kk, i
    LOGICAL compute_covariant, compute_contravariant
    LOGICAL compute_derivatives, b(32)

    IF (flag > 7) THEN
        WRITE(6,*) 'cylinder_basis_complete : flag > 7'
        STOP
    ENDIF

!!  # flag = 0      NA
!!  # flag = 1      Covariant basis only
!!  # flag = 2      NA
!!  # flag = 3      Covariant + contravariant basis
!!  # flag = 4      Derivatives only
!!  # flag = 5      Derivatives + covariant basis
!!  # flag = 6      NA
!!  # flag = 7      Covariant + contravariant + derivatives


    DO i = 1,bit_size(flag)
        b(i) = btest(flag,i-1)          
    END DO

    compute_covariant     = b(1) .or. b(2)
    compute_contravariant = b(2)
    compute_derivatives   = b(3)

    pi4 = pi2*pi2

!!    # Cylinder mapping
!!    #
!!    #    xp =  r_cyl*cos(pi2*x)
!!    #    yp =  r_cyl*sin(pi2*x)
!!    #    zp =  h_cyl*y
!!    #
!!    #  Express X,Y,Z as 
!!    # 
!!    #     X(x,y) = g(x,y)*f(x,y) == g(1)*f(1)
!!    #     Y(x,y) = g(x,y)*f(x,y) == g(2)*f(2)
!!    #     Z(x,y) = g(x,y)*f(x,y) == g(3)*f(3)
!!    # 
!!    # Compute derivatives and use product rule to get higher order
!!    # derivatives, i.e. 
!!    #
!!    #     dX/dx = g(1)*fx(1) + gx(1)*f(1)
!!    #     dY/dx = g(2)*fx(2) + gx(2)*f(2)
!!    #     dZ/dx = g(3)*fx(3) + gx(3)*f(3)
!!    #

  
!!  # Map (x,y) in [0,1]x[0,1] to (theta, z). 
!!  # Compute derivatives dtheta/dx, dtheta/dy, etc.
    CALL map_comp2cylinder_derivs(x,y,theta,z, & 
                                  thetax, thetay, zx, zy)

    R     = r_cyl
    Rx    = 0  
    Rxx   = 0 
    Ry    = 0  
    Ryy   = 0 
    Rxy   = 0

!!  X = g(1)*f(1) = R*cos(theta)
    g(1)  = R
    gx(1) = Rx
    gy(1) = Ry

    f(1)  = cos(theta)
    fx(1) = -sin(theta)*thetax
    fy(1) = 0;

!!  Y = g(2)*f(2) = R*sin(theta)
    g(2)  = R
    gx(2) = Rx
    gy(2) = Ry

    f(2)  = sin(theta)
    fx(2) = cos(theta)*thetax
    fy(2) = 0;

!!  Z = g(3)*f(3) = h_cyl*y
    g(3)  = h_cyl
    gx(3) = 0
    gy(3) = 0

    f(3)  = y
    fx(3) = 0
    fy(3) = 1

    IF (compute_covariant) THEN
        DO k = 1,3
            t(k,1) = gx(k)*f(k) + g(k)*fx(k);
            t(k,2) = gy(k)*f(k) + g(k)*fy(k);
        END DO
    ENDIF

    IF (compute_contravariant) then
        CALL map_contravariant(t,tinv)
    ENDIF

    IF (compute_derivatives) then
        fxx(1) = -cos(theta)*thetax**2
        fyy(1) = 0
        fxy(1) = 0

        fxx(2) = -sin(theta)*thetax**2
        fyy(2) = 0
        fxy(2) = 0

        fxx(3) = 0
        fxy(3) = 0
        fyy(3) = 0

        gxx(1) = Rxx
        gxx(2) = Rxx
        gxx(3) = 0

        gyy(1) = Ryy
        gyy(2) = Ryy
        gyy(3) = 0

        gxy(1) = Rxy
        gxy(2) = Rxy
        gxy(3) = 0

        DO k = 1,3
!!          # d(t1)/dx = d(g*fx + gx*f)/dx
            tderivs(k,1,1) = g(k)*fxx(k) + 2*fx(k)*gx(k) + gxx(k)*f(k)

!!          # d(t1)/dy = d(g*fx + gx*f)/dy       
            tderivs(k,1,2) = g(k)*fxy(k) + gy(k)*fx(k) + & 
                            gx(k)*fy(k) + gxy(k)*f(k)

!!          # d(t2)/dx = d(g*fy + gy*f)/dx       
            tderivs(k,2,1) = g(k)*fxy(k) + gx(k)*fy(k) + & 
                         gy(k)*fx(k) + gxy(k)*f(k)

!!          # d(t2)/dy = d(g*fy + gy*f)/dy         
            tderivs(k,2,2) = g(k)*fyy(k) + 2*fy(k)*gy(k) + gyy(k)*f(k)
          END DO
    ENDIF

END SUBROUTINE cylinder_basis_complete

SUBROUTINE map_comp2cylinder(xc,yc,theta,z)
    IMPLICIT NONE

    DOUBLE PRECISION xc,yc,theta, z

    DOUBLE PRECISION thetax, thetay, zx, zy

!!    # Map xc in [0,1] to theta in [-pi,pi]
!!    # Map yc in [0,1] to phi in [-pi/2,pi/2]      

    call map_comp2cylinder_derivs(xc,yc,theta,z, & 
                    thetax, thetay, zx, zy)

END SUBROUTINE map_comp2cylinder

SUBROUTINE map_comp2cylinder_derivs(xc,yc,theta,z, & 
              thetax, thetay, zx, zy)
    IMPLICIT NONE

    DOUBLE PRECISION xc,yc,theta, z
    DOUBLE PRECISION  thetax, thetay, zx, zy

    DOUBLE PRECISION r_cyl, h_cyl
    COMMON /cylinder_comm/ r_cyl, h_cyl

    DOUBLE PRECISION pi, pi2
    COMMON /compi/ pi, pi2

    DOUBLE PRECISION tr1, tr2, zr1, zr2


!!  # Map xc in [0,1] to theta in [0,2*pi]
!!  # Map yc in [0,1] to phi in [-pi/2,pi/2]      


    tr1 = 0
    tr2 = pi2

    theta = tr1 + (tr2-tr1)*xc
    thetax = tr2-tr1
    thetay = 0

    zr1 = 0
    zr2 = h_cyl

    z = zr1 + (zr2-zr1)*yc
    zx = 0
    zy = zr2-zr1

END SUBROUTINE map_comp2cylinder_derivs

SUBROUTINE map_cylinder2comp(theta,z,xc,yc)
    IMPLICIT NONE

    DOUBLE PRECISION xc,yc,theta,z

    DOUBLE PRECISION pi, pi2
    COMMON /compi/ pi, pi2

    DOUBLE PRECISION r_cyl, h_cyl
    COMMON /cylinder_comm/ r_cyl, h_cyl

    DOUBLE PRECISION tr1, tr2, zr1, zr2

!!  # Map xc in [0,1] to theta in [0,2*pi]
!!  # Map yc in [0,1] to phi in [-pi/2,pi/2]      

    tr1 = 0
    tr2 = pi2

    zr1 = 0
    zr2 = h_cyl

    xc = (theta-tr1)/(tr2-tr1)
    yc = (z-zr1)/(zr2-zr1)

END SUBROUTINE map_cylinder2comp

SUBROUTINE map2cylinder(xp,yp,zp,theta,z)
    IMPLICIT NONE

    DOUBLE PRECISION xp,yp,zp,theta,z

    DOUBLE PRECISION pi, pi2
    COMMON /compi/ pi, pi2

    DOUBLE PRECISION r_cyl, h_cyl
    common /cylinder_comm/ r_cyl, h_cyl

    DOUBLE PRECISION r1, r

!!  xp = R*cos(pi2*xc)
!!  yp = R*sin(pi2*xc)
!!  zp = z

  
    theta = atan2(yp,xp)    !! returns value in [-pi, pi]
    if (theta < 0) then
        theta = theta + pi2
    endif

    z = zp

END SUBROUTINE map2cylinder

SUBROUTINE map2comp(xp,yp,zp,xc,yc)
    IMPLICIT NONE

    DOUBLE PRECISION xp,yp,zp, xc,yc

    DOUBLE PRECISION pi, pi2
    COMMON /compi/ pi, pi2

    DOUBLE PRECISION theta, z

    CALL map2cylinder(xp,yp,zp,theta,z)

    CALL map_cylinder2comp(theta,z,xc,yc)

END SUBROUTINE map2comp

