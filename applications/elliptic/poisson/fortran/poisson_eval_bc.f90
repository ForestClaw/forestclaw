double precision function poisson_fort_eval_bc(iface,t,x,y)
    implicit none

    integer iface
    double precision x,y,t

    integer bctype(0:3)
    common /comm_bc/ bctype

    integer normals(0:3,2)
    double precision q,grad(2), qn, a,b

    data normals /-1,1,0,0,0,0,-1,1/

    call poisson_qexact_gradient(x,y,q,grad)

    !! q_n = grad q \cdot n
    qn = normals(iface,1)*grad(1) + normals(iface,2)*grad(2)

    !! bc_type is set in options .ini file as [multigrid] boundary_conditions
    if (bctype(iface) .eq. 1) then
        a = 1
        b = 0
    elseif (bctype(iface) .eq. 2) then
        a = 0
        b = 1
    endif

    poisson_fort_eval_bc = a*q + b*qn

    return
    
end function poisson_fort_eval_bc