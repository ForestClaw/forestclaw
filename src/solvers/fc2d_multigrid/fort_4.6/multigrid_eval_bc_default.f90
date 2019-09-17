
double precision function multigrid_fort_eval_bc_default(iface,t,x,y)
    implicit none

    integer iface
    double precision t, x,y

    multigrid_fort_eval_bc_default = 0

    return
    
end function multigrid_fort_eval_bc_default