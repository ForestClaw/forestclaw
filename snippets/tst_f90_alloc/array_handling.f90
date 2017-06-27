SUBROUTINE store_ptrs(fc_x_comp_sp_ptr,fc_concen_ptr)
    use iso_c_binding
    use test_mod
    IMPLICIT NONE

    TYPE(c_ptr) :: fc_x_comp_sp_ptr, fc_concen_ptr

    fc_x_comp_sp_ptr = c_loc(x_comp_sp)
    fc_concen_ptr = c_loc(concen)
end

SUBROUTINE copy_ptrs2mod(fc_x_comp_sp_ptr,fc_concen_ptr)
    use iso_c_binding
    use test_mod
    IMPLICIT NONE

    TYPE(c_ptr) :: fc_x_comp_sp_ptr, fc_concen_ptr
    double precision, pointer :: xp(:)   !! Needed only for 2d and higher dim arrays
    INTEGER :: i,j

    CALL c_f_pointer(fc_x_comp_sp_ptr,x_comp_sp,[1])

    CALL c_f_pointer(fc_concen_ptr,xp,[1])
    concen(1:nx,1:nx) => xp

    do i = 1,nx
        write(6,'(F12.0)') x_comp_sp(i)
    end do

end

subroutine deallocate_arrays(fc_x_comp_sp_ptr,fc_concen_ptr)
    use iso_c_binding
    use test_mod
    implicit none

    type(c_ptr) :: fc_x_comp_sp_ptr, fc_concen_ptr
    double precision, pointer :: xp(:)

    CALL c_f_pointer(fc_concen_ptr,xp,[1])
    concen(1:nx,1:nx) => xp
    DEALLOCATE(concen)  !! May also be okay to just de-allocate xp

    CALL c_f_pointer(fc_x_comp_sp_ptr,xp,[1])
    x_comp_sp(1:nx) => xp
    DEALLOCATE(x_comp_sp)  

end

