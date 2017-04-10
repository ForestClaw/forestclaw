MODULE test_mod
    IMPLICIT NONE

    DOUBLE PRECISION, POINTER :: x_comp_sp(:) 
    INTEGER :: nx_submet = 5

END MODULE test_mod

MODULE fc_arrays
    USE iso_c_binding        
    IMPLICIT NONE

    DOUBLE PRECISION, POINTER :: fc_x_comp_sp(:) 

contains

END MODULE fc_arrays

SUBROUTINE assign_ptrs(fc_ptr)
    use iso_c_binding
    use fc_arrays
    IMPLICIT NONE

    TYPE(c_ptr) :: fc_ptr

    call c_f_pointer(fc_ptr,fc_x_comp_sp,[1]);

end


SUBROUTINE use_array(value)
    use iso_c_binding
    USE test_mod
    USE fc_arrays
    IMPLICIT NONE

    double precision value
    integer i

    ALLOCATE(fc_x_comp_sp(nx_submet))
    x_comp_sp(1:nx_submet) => fc_x_comp_sp

    do i = 1,nx_submet
        x_comp_sp(i) = value
    enddo 


end





