!! # check to see if value exceeds threshold

integer function euler3d_pressure_exceeds_th(blockno,meqn,& 
        qval,qmin,qmax,quad, & 
        dx,dy,dz,xc,yc,zc,ivar_threshold, threshold, &
        init_flag, is_ghost)
    use setprob_mod, only : gamma1
    implicit none
    
    integer :: meqn, ivar_threshold
    double precision :: qval(meqn),qmin(meqn),qmax(meqn),threshold
    double precision :: quad(-1:1,-1:1,-1:1,meqn)
    double precision :: dx,dy, dz,xc, yc,zc
    integer :: blockno, init_flag, refine
    logical(kind=4) :: is_ghost

    double precision rho, energy, kinetic, pressure    


    rho = qval(1)
    energy = qval(5)
    kinetic = 0.5*(qval(2)**2 + qval(3)**2 + qval(4)**2)/rho
    pressure = gamma1*(energy - kinetic)

    refine = 0
    if (pressure .gt. threshold) then
        refine = 1
    endif

    euler3d_pressure_exceeds_th = refine

end function euler3d_pressure_exceeds_th
