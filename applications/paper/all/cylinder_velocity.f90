!!   # ------------------------------------------------------------
!!   # Prescribes velocity fields for the unit sphere.
!!   # 
!!   # ------------------------------------------------------------
   
SUBROUTINE cylinder_velocity_derivs(x,y,t, u,vcart,derivs, cart_flag)
    IMPLICIT NONE

    DOUBLE PRECISION x, y, u(2), vcart(3), t
    DOUBLE PRECISION derivs(4)
    INTEGER cart_flag

    DOUBLE PRECISION pi, pi2
    COMMON /compi/ pi, pi2

    INTEGER example
    COMMON /example_comm/ example  

    DOUBLE PRECISION revs_per_s, v_speed
    COMMON /stream_comm/ revs_per_s, v_speed

    DOUBLE PRECISION r_cyl, h_cyl
    COMMON /cylinder_comm/ r_cyl, h_cyl

    DOUBLE PRECISION s, a, xc1, yc1
    DOUBLE PRECISION theta, thetax, thetay
    DOUBLE PRECISION z, zx, zy

!!    CALL map_comp2cylinder_derivs(x,y,theta,z,thetax,thetay, zx, zy)

    IF (example .eq. 0) THEN
        cart_flag = 0
        u(1) = revs_per_s
        u(2) = v_speed
        derivs(1) = 0
        derivs(2) = 0
        derivs(3) = 0
        derivs(4) = 0
    ENDIF

END SUBROUTINE cylinder_velocity_derivs



!!    # ----------------------------------------------------------------
!!    #                       Public interface
!!    # ----------------------------------------------------------------


!!    # ---------------------------------------------
!!    # Called from setaux
!!    # 
!!    #    -- used to compute velocity at faces
!!    # ---------------------------------------------

SUBROUTINE velocity_components_cart(x,y,t,vcart)
    IMPLICIT NONE

    DOUBLE PRECISION x,y,t, u(2), vcart(3)

    INTEGER mapping 
    COMMON /mapping_comm/ mapping

    DOUBLE PRECISION t1(3), t2(3), uderivs(4)
    INTEGER cart_flag, k

    if (mapping .eq. 0) then
        CALL cylinder_velocity_derivs(x,y,t, u,vcart,uderivs,cart_flag)
    else
        CALL latlong_velocity_derivs(x,y,t, u,vcart,uderivs,cart_flag)
    endif        

    IF (cart_flag .eq. 0) THEN
!!      # Velocity components are given in spherical components
!!      # and must be converted to Cartesian
        CALL map_covariant_basis(x, y, t1,t2)

        DO k = 1,3
            vcart(k) = u(1)*t1(k) + u(2)*t2(k)
        END DO
    ENDIF

END SUBROUTINE velocity_components_cart


!!     # ------------------------------------------------------------
!!     # Called from map_divergence
!!     # 
!!     #    -- Needed to define ODE system to get exact solution
!!     # ------------------------------------------------------------

SUBROUTINE velocity_derivs(x,y,t, u, vcart, derivs, cart_flag)
    IMPLICIT NONE

    DOUBLE PRECISION x,y,t, u(2), vcart(3), derivs(4)

    INTEGER mapping 
    COMMON /mapping_comm/ mapping

    DOUBLE PRECISION t1(3), t2(3), t1n2, t2n2, map_dot
    DOUBLE PRECISION t1inv(3), t2inv(3)
    INTEGER cart_flag

    if (mapping .eq. 0) then
        CALL cylinder_velocity_derivs(x,y,t, u,vcart,derivs,cart_flag)
    else
        CALL latlong_velocity_derivs(x,y,t, u,vcart,derivs,cart_flag)
    endif

    IF (cart_flag .eq. 1) THEN
!!      # Velocity components are given in Cartesian components
        CALL map_contravariant_basis(x, y, t1inv,t2inv)
        u(1) = map_dot(vcart,t1inv)
        u(2) = map_dot(vcart,t2inv)

!!      # Need to convert derivatives to derivatives with respect
!!      # to computational coordinates.        
    
    ENDIF
END SUBROUTINE velocity_derivs

!!    # ------------------------------------------------------------
!!    # Called from qexact
!!    # 
!!    #  -- components relative to basis are needed.
!!    # ------------------------------------------------------------

SUBROUTINE velocity_components(x,y,t,u)
    IMPLICIT NONE

    DOUBLE PRECISION x,y,t, u(2), vcart(3), derivs(4)
    INTEGER cart_flag

    CALL velocity_derivs(x,y,t, u,vcart,derivs,cart_flag)

END SUBROUTINE velocity_components



