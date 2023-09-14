subroutine clawpack46_bc3_default(meqn,mbc,mx,my,mz, & 
        xlower, ylower,zlower, dx,dy,dz,q,maux,aux,t,dt,mthbc)
    !!  =====================================================
    !!
    !!  # Standard boundary condition choices for claw3
    !!
    !!  # At each boundary  k = 1 (xlower),  2 (xupper), 
    !!  #                       3 (ylower),  4 (yupper),
    !!  #                       5 (zlower),  6 (zupper):
    !!  #   mthbc(k) =  0  for user-supplied BC's (must be inserted!)
    !!  #            =  1  for zero-order extrapolation
    !!  #            =  2  for periodic boundary coniditions
    !!  #            =  3  for solid walls, assuming this can be implemented
    !!  #                  by reflecting the data about the boundary and then
    !!  #                  negating the 2'nd (for k=1,2) or 3'rd (for k=3,4)
    !!  #                  or 4'th (for k=5,6) component of q.
    !!  ------------------------------------------------
    !!
    !!  # Extend the data from the interior cells (1:mx, 1:my, 1:mz)
    !!  # to a layer of mbc ghost cells outside the region.

    implicit none

    integer :: meqn,mbc,mx,my,mz,maux
    double precision :: xlower, ylower, zlower, dx,dy,dz, t, dt

    double precision :: q(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc, meqn)
    double precision :: aux(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc, maux)
    integer :: mthbc(6)

    integer :: i,j,k,m, ibc, jbc, kbc

!! -------------------------------------------------------
!!      # left boundary (xlower):
!! -------------------------------------------------------

    go to (100,110,120,130) mthbc(1)+1
    !! fall through if mthbc(1) = -1. In this case, we are at an internal
    !! patch boundary and shouldn't apply any physical BCs.
    goto 199

100  continue
    !! # user-specified boundary conditions go here in place of error output
    write(6,*) '*** ERROR *** mthbc(1)=0 and no BCs specified in bc3'
    stop
    go to 199

110   continue
    !! # zero-order extrapolation:
    do m = 1,meqn
        do ibc = 1,mbc
            do j = 1-mbc, my+mbc
                do k = 1-mbc, mz+mbc
                    q(1-ibc,j,k,m) = q(1,j,k,m)
                end do
            end do
        end do
    end do
    go to 199

120   continue
    !! Periodic
    write(6,*) 'bc3 : Horizontal periodic BCs handled by p4est.'
    go to 199

130   continue
    !! # solid wall (assumes 2'nd component is velocity or momentum in x):
    do m = 1,meqn
        do ibc = 1,mbc
            do j = 1-mbc, my+mbc
                do k = 1-mbc, mz+mbc
                    q(1-ibc,j,k,m) = q(ibc,j,k,m)
                end do
            end do
        end do
    end do

    !! # negate the normal velocity:
    do ibc = 1,mbc
        do j = 1-mbc, my+mbc
            do k = 1-mbc, mz+mbc
                q(1-ibc,j,k,2) = -q(ibc,j,k,2)
            end do
        end do
    end do
    go to 199

199   continue
    
!! -------------------------------------------------------
!! # right boundary (xupper):
!! -------------------------------------------------------
    go to (200,210,220,230) mthbc(2)+1
    !!  Fall through if mthbc(2) = -1.   Boundary is an internal patch
    !! boundary.
    goto 299

200 continue
    !! # user-specified boundary conditions go here in place of error output
    write(6,*) '*** ERROR *** mthbc(2)=0 and no BCs specified in bc3'
    stop
    go to 299

210 continue
    !! # zero-order extrapolation:
    do m = 1,meqn
        do ibc = 1,mbc
            do j = 1-mbc, my+mbc
                do k = 1-mbc, mz+mbc
                    q(mx+ibc,j,k,m) = q(mx,j,k,m)
                end do
            end do
        end do
    end do
    go to 299

220 continue
    !! # periodic:  
    write(6,*) 'bc3 : Horizontal periodic BCs handled by p4est.'
    goto 299

230 continue
    !! # solid wall (assumes 2'nd component is velocity or momentum in x):
    do m = 1,meqn
        do ibc = 1,mbc
            do j = 1-mbc, my+mbc
                do k = 1-mbc, mz+mbc
                    q(mx+ibc,j,k,m) = q(mx+1-ibc,j,k,m)
                end do
            end do
        end do
    end do

    !! # negate the normal velocity:
    do ibc = 1,mbc
        do j = 1-mbc, my+mbc
            do k = 1-mbc, mz+mbc
                q(mx+ibc,j,k,2) = -q(mx+1-ibc,j,k,2)
            end do
        end do
    end do
    go to 299

299 continue

!! -------------------------------------------------------
!!      # front boundary (ylower):
!! -------------------------------------------------------
    go to (300,310,320,330) mthbc(3)+1
    !!  Fall through if mthbc(3) = -1.   Boundary is an internal patch
    !! boundary.
    goto 399

300 continue
    !! # user-specified boundary conditions go here in place of error output
    write(6,*) '*** ERROR *** mthbc(3)=0 and no BCs specified in bc3'
    stop
    go to 399

310 continue
    !! # zero-order extrapolation:
    do m = 1,meqn
        do jbc = 1,mbc
            do i = 1-mbc, mx+mbc
                do k = 1-mbc, mz+mbc
                    q(i,1-jbc,k,m) = q(i,1,k,m)
                end do
            end do
        end do
    end do
    go to 399

320 continue
    !! # periodic:  
    write(6,*) 'bc3 : Horizontal periodic BCs handled by p4est.'
    go to 399

330 continue
    !! # solid wall (assumes 3'rd component is velocity or momentum in y):
    do m = 1,meqn
        do jbc = 1,mbc
            do i = 1-mbc, mx+mbc
                do k = 1-mbc, mz+mbc
                  q(i,1-jbc,k,m) = q(i,jbc,k,m)
                end do
            end do
        end do
    end do
    !! # negate the normal velocity:
    do jbc = 1,mbc
        do i = 1-mbc, mx+mbc
            do k = 1-mbc, mz+mbc
               q(i,1-jbc,k,3) = -q(i,jbc,k,3)
            end do
        end do
    end do
    go to 399

  399 continue

!! -------------------------------------------------------
!!      # back boundary (yupper):
!! -------------------------------------------------------
    go to (400,410,420,430) mthbc(4)+1
    !!  Fall through if mthbc(3) = -1.   Boundary is an internal patch
    !! boundary.
    goto 499

400 continue
    !! # user-specified boundary conditions go here in place of error output
    write(6,*) '*** ERROR *** mthbc(4)=0 and no BCs specified in bc3'
    stop
    go to 499

  410 continue
    !! # zero-order extrapolation:
    do m = 1,meqn
        do jbc = 1,mbc
            do i = 1-mbc, mx+mbc
                do k = 1-mbc, mz+mbc
                    q(i,my+jbc,k,m) = q(i,my,k,m)
                end do
            end do
        end do
    end do
    go to 499

420 continue
    !! # periodic:  
    write(6,*) 'bc3 : Horizontal periodic BCs handled by p4est.'
    goto 499

430 continue
    !! # solid wall (assumes 3'rd component is velocity or momentum in y):
    do m = 1,meqn
        do jbc = 1,mbc
            do i = 1-mbc, mx+mbc
                do k = 1-mbc, mz+mbc
                    q(i,my+jbc,k,m) = q(i,my+1-jbc,k,m)
                end do
            end do
        end do
    end do

    !! # negate the normal velocity:
    do jbc = 1,mbc
        do i = 1-mbc, mx+mbc
            do k = 1-mbc, mz+mbc
                q(i,my+jbc,k,3) = -q(i,my+1-jbc,k,3)
            end do
        end do
    end do
    go to 499

499 continue


!! -------------------------------------------------------
!!      # bottom boundary (zlower):
!! -------------------------------------------------------
    goto (500,510,520,530) mthbc(5)+1
    !!  Fall through if mthbc(3) = -1.   Boundary is an internal patch
    !! boundary.
    goto 599

500 continue
    !! # user-specified boundary conditions go here in place of error output
    write(6,*) '*** ERROR *** mthbc(5)=0 and no BCs specified in bc3'
    stop
    go to 599

510 continue
    !! # zero-order extrapolation:
    do m = 1,meqn
        do kbc = 1,mbc
            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    q(i,j,1-kbc,m) = q(i,j,1,m)
                end do
            end do
        end do
    end do
    go to 599

520 continue
    !! # periodic:  
    do m = 1,meqn
        do kbc = 1,mbc
            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    q(i,j,1-kbc,m) = q(i,j,mz+1-kbc,m)
                end do
            end do
        end do
    end do
    go to 599

530 continue
    !! # solid wall (assumes 4'rd component is velocity or momentum in y):
    do m = 1,meqn
        do kbc = 1,mbc
            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    q(i,j,1-kbc,m) = q(i,j,kbc,m)
                end do
            end do
        end do
    end do
    !! # negate the normal velocity:
    do kbc = 1,mbc
        do i = 1-mbc, mx+mbc
            do j = 1-mbc, my+mbc
                q(i,j,1-kbc,4) = -q(i,j,kbc,4)
            end do
        end do
    end do
    go to 599

  599 continue

!! -------------------------------------------------------
!!      # boundary (zupper):
!! -------------------------------------------------------

    goto (600,610,620,630) mthbc(6)+1
    !!  Fall through if mthbc(3) = -1.   Boundary is an internal patch
    !! boundary.
    goto 699

600 continue
    !! # user-specified boundary conditions go here in place of error output
    write(6,*) '*** ERROR *** mthbc(6)=0 and no BCs specified in bc3'
    stop
    goto 699

610 continue
    !! # zero-order extrapolation:
    do m = 1,meqn
        do kbc = 1,mbc
            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    q(i,j,mz+kbc,m) = q(i,j,mz,m)
                end do
            end do
        end do
    end do
    go to 699

620 continue
    !! # periodic:  
    do m = 1,meqn
        do kbc = 1,mbc
            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    q(i,j,mz+kbc,m) = q(i,j,kbc,m)
                end do
            end do
        end do
    end do
    go to 699

630 continue
    !! # solid wall (assumes 3'rd component is velocity or momentum in y):
    do m = 1,meqn
        do kbc = 1,mbc
            do i = 1-mbc, mx+mbc
                do j = 1-mbc, my+mbc
                    q(i,j,mz+kbc,m) = q(i,j,mz+1-kbc,m)
                end do
            end do
        end do
    end do
    !! # negate the normal velocity:
    do kbc = 1,mbc
       do i = 1-mbc, mx+mbc
            do j = 1-mbc, my+mbc
                q(i,j,mz+kbc,4) = -q(i,j,mz+1-kbc,4)
            end do
        end do
    end do
    go to 699

699 continue

    return
end subroutine clawpack46_bc3_default

