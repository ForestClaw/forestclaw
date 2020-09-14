SUBROUTINE mapc2m_cylinder(xc,yc,xp,yp,zp)
    IMPLICIT NONE

    DOUBLE PRECISION xc,yc,xp,yp,zp

    DOUBLE PRECISION pi, pi2
    COMMON /compi/ pi, pi2

    DOUBLE PRECISION R_cyl, h_cyl
    COMMON /cylinder_comm/ r_cyl, h_cyl

    DOUBLE PRECISION r1, R, theta, z

    CALL map_comp2cylinder(xc,yc,theta,z)

    
    xp = R_cyl*cos(theta)
    yp = R_cyl*sin(theta)
    zp = h_cyl*z

END SUBROUTINE MAPC2M_CYLINDER

