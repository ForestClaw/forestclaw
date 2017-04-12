MODULE test_mod
    IMPLICIT NONE

    DOUBLE PRECISION, POINTER :: x_comp_sp(:)     
    DOUBLE PRECISION, POINTER :: concen(:,:) 

    INTEGER :: nx = 5

END MODULE test_mod


SUBROUTINE use_array(value1,value2)
    USE test_mod
    IMPLICIT NONE

    DOUBLE PRECISION value1, value2
    INTEGER i,j

    ALLOCATE(x_comp_sp(nx), concen(nx,nx))

    do i = 1,nx
        x_comp_sp(i) = value1
        do j = 1,nx
            concen(i,j) = value2
        end do
    end do 



end





