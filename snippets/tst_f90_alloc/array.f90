MODULE test_mod
    IMPLICIT NONE

    DOUBLE PRECISION, POINTER :: x_comp_sp(:) 
    INTEGER :: nx_submet = 5

END MODULE test_mod

SUBROUTINE store_ptrs(fc_ptr)
    use iso_c_binding
    use test_mod
    IMPLICIT NONE

    TYPE(c_ptr) :: fc_ptr

    fc_ptr = c_loc(x_comp_sp)
end

SUBROUTINE copy_ptrs2mod(fc_ptr)
    use iso_c_binding
    use test_mod
    IMPLICIT NONE

    TYPE(c_ptr) :: fc_ptr
    INTEGER :: i

    CALL c_f_pointer(fc_ptr,x_comp_sp,[1])

    DO i = 1,nx_submet
        WRITE(6,*) x_comp_sp(i)
    ENDDO 

end

subroutine deallocate_arrays(fc_ptr)
    use iso_c_binding
    use test_mod
    implicit none

    type(c_ptr) :: fc_ptr

    CALL c_f_pointer(fc_ptr,x_comp_sp,[1])

    DEALLOCATE(x_comp_sp)

end


SUBROUTINE use_array(value)
    use iso_c_binding
    USE test_mod
    IMPLICIT NONE

    DOUBLE PRECISION value
    INTEGER i

    ALLOCATE(x_comp_sp(nx_submet))

    do i = 1,nx_submet
        x_comp_sp(i) = value
    enddo 

end





