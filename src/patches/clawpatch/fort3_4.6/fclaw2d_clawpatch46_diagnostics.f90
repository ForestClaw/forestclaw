!! # ----------------------------------------------------------------------------------
!! # Output and diagnostics
!! # ----------------------------------------------------------------------------------
subroutine fclaw2d_clawpatch46_fort3_conservation_check(mx,my,mz, & 
    mbc,mfields,dx,dy,dz,area,q,sum,c_kahan)
    implicit none

    integer :: mx,my,mz,mbc,mfields
    double precision :: dx, dy,dz
    double precision :: sum(mfields), c_kahan(mfields)
    double precision :: q(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,mfields)
    double precision :: area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

    integer :: i,j,k,m
    double precision :: t, y, area_ij, vol_ij
    integer*8 :: cont, get_context
    logical :: fclaw2d_map_is_used

    logical :: use_kahan

    use_kahan = .true.

    cont = get_context()

    area_ij = dx*dy  !! Area in each mesh cell is constant
    do m = 1,mfields
        do k = 1,mz
            do j = 1,my
                do i = 1,mx
                    if (fclaw2d_map_is_used(cont)) then
                        area_ij = area(i,j)  !! Area varies
                    endif
                    vol_ij = area_ij*dz
                    if (use_kahan) then
                        y = q(i,j,k,m)*vol_ij - c_kahan(m)
                        t = sum(m) + y
                        c_kahan(m) = (t-sum(m)) - y
                        sum(m) = t
                    else
                        sum(m) = sum(m) + q(i,j,k,m)*vol_ij
                    endif
                enddo 
            enddo
        end do
    end do
end subroutine fclaw2d_clawpatch46_fort3_conservation_check

!! # Compute area of a patch
double precision function fclaw2d_clawpatch46_fort3_compute_patch_area( & 
         mx,my, mz, mbc,dx,dy,dz, area)
    implicit none

    integer :: mx,my, mz, mbc
    double precision :: dx, dy, dz
    double precision :: area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

    integer :: i,j, k
    integer*8 :: cont, get_context
    logical :: fclaw2d_map_is_used
    double precision :: sum

    cont = get_context()

    if (fclaw2d_map_is_used(cont)) then
        sum = 0       
        do j = 1,my
            do i = 1,mx
                sum = sum + area(i,j)
            end do
        end do
    else
        sum = dx*dy*mx*my
    endif

    fclaw2d_clawpatch46_fort3_compute_patch_area = sum

end function fclaw2d_clawpatch46_fort3_compute_patch_area


subroutine fclaw2d_clawpatch46_fort3_compute_error_norm( & 
     blockno, mx,my,mz,mbc,mfields,dx,dy,dz,area,error,error_norm)
    implicit none

    integer :: mx,my,mz, mbc,mfields, blockno
    double precision :: dx, dy, dz
    double precision :: error_norm(mfields,3)
    double precision :: error(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,mfields)
    double precision :: area(-mbc:mx+mbc+1,-mbc:my+mbc+1)

    integer :: i,j,k,m
    double precision :: dxdydz, eij, vol_ij

    integer*8 :: cont, get_context
    logical :: fclaw2d_map_is_used

    cont = get_context()

    !! # error_norm(:) comes in with values;  do not initialize it here!
    dxdydz = dx*dy*dz
    do m = 1,mfields
       if (fclaw2d_map_is_used(cont)) then
            do k = 1,mz
                do j = 1,my
                    do i = 1,mx
                        vol_ij = area(i,j)*dz
                        eij = abs(error(i,j,k,m))
                        error_norm(m,1) = error_norm(m,1) + & 
                            eij*vol_ij
                        error_norm(m,2) = error_norm(m,2) + & 
                            eij**2*vol_ij
                        error_norm(m,3) = max(eij,error_norm(m,3))
                    end do
                end do
            end do
        else
            do k = 1,mz
                do j = 1,my
                    do i = 1,mx
                        eij = abs(error(i,j,k,m))
                        error_norm(m,1) = error_norm(m,1) + & 
                            eij*dxdydz
                        error_norm(m,2) = error_norm(m,2) + & 
                            eij**2*dxdydz
                        error_norm(m,3) = max(eij,error_norm(m,3))
                    end do
                end do
            end do
        endif
    end do

end subroutine fclaw2d_clawpatch46_fort3_compute_error_norm
