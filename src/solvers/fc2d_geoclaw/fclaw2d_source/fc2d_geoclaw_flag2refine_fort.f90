SUBROUTINE fc2d_geoclaw_flag2refine(blockno, meqn,maux,qvec, auxvec, dx,dy, & 
    xc,yc,t,level, maxlevel, init_flag, is_coarsening)

    USE geoclaw_module, ONLY : dry_tolerance, sea_level
    USE refinement_module, only : wave_tolerance, speed_tolerance
    
    IMPLICIT NONE

    INTEGER :: blockno, init_flag,level, meqn, maux, maxlevel
    DOUBLE PRECISION :: qvec(meqn),auxvec(maux), xc,yc,dx,dy,t
    logical :: is_coarsening

    INTEGER :: flag_patch, m, max_num_speeds
    DOUBLE PRECISION :: eta, th_factor, speed
    LOGICAL :: allowflag

    !! Don't coarsen on initial refinement
    if (init_flag .ne. 0 .and. is_coarsening) then
        flag_patch = 0
        return
    endif


    if (.not. allowflag(xc,yc,t,level)) then
        !! It isn't clear what this means;  do we not refine the entire patch
        !! if a single cell cannot be refined? 
        flag_patch = -1
        return
    endif

    max_num_speeds = min(size(speed_tolerance),maxlevel)

    th_factor = 1
    if (is_coarsening) then
        !! Coarsening factor should be 0.5 refinement factor.
        th_factor = 0.5
    endif

    flag_patch = 0
    if (qvec(1) > dry_tolerance) then

        !! Check wave height criteria
        eta = qvec(1) + auxvec(1)
        if (abs(eta - sea_level) > th_factor*wave_tolerance) then
            flag_patch = 1
            return
        endif

        !! Check speed criteria
        speed = sqrt(qvec(2)**2 + qvec(3)**2) / qvec(1)
        DO m = 1,max_num_speeds
            IF (speed > th_factor*speed_tolerance(m) .AND. level <= m) THEN
                flag_patch = 1
                RETURN
            ENDIF
        ENDDO
    ENDIF

end subroutine fc2d_geoclaw_flag2refine