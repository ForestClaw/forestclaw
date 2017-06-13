MODULE test_mod
    implicit none
    DOUBLE PRECISION, POINTER :: concen(:,:)
    INTEGER :: nsize

END MODULE test_mod

SUBROUTINE set_ptr(fc_concen,n)
    USE iso_c_binding
    USE test_mod

    DOUBLE PRECISION, POINTER :: fc_ptr(:)
    INTEGER :: n
    TYPE(c_ptr)  :: fc_concen

    nsize = n

    CALL c_f_POINTER(fc_concen,fc_ptr,[1])
    concen(0:n-1,0:n-1)  => fc_ptr

END SUBROUTINE

SUBROUTINE assign_array()    
    USE test_mod, ONLY: nsize, concen
    IMPLICIT NONE

    INTEGER :: i,j

    WRITE(6,*) 'Shape : ', SHAPE(concen)

    DO j = 0,nsize-1
        DO i = 0,nsize-1
            concen(i,j) = nsize*i + j
        end do
    ENDDO

END SUBROUTINE

SUBROUTINE print_array()
    USE test_mod,  ONLY  : nsize, concen
    IMPLICIT none

    INTEGER :: i,j

    DO i = 0,nsize-1
        DO j = 0,nsize-1
            WRITE(6,'(2F24.3)') concen(i,j)
        END DO
    END DO

END SUBROUTINE print_array
