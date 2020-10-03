SUBROUTINE mapc2m_cylinder(xc,yc,xp,yp,zp)
    IMPLICIT NONE

    DOUBLE PRECISION xc,yc,xp,yp,zp

    DOUBLE PRECISION pi, pi2
    COMMON /compi/ pi, pi2

    DOUBLE PRECISION R_cyl, h_cyl
    COMMON /cylinder_comm/ r_cyl, h_cyl

    DOUBLE PRECISION r1, R, theta, z, yc1


    !! finite cylinder : 
    if (yc .lt. 0 .or. yc .gt. 1) then
!!        write(6,*) 'yc not in [0,1];  yc = ', yc
    endif

    !!yc1 = mod(yc,1.0)
    
    CALL map_comp2cylinder(xc,yc,theta,z)
    
    xp = R_cyl*cos(theta)
    yp = R_cyl*sin(theta)
    zp = z

END SUBROUTINE MAPC2M_CYLINDER

