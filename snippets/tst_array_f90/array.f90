MODULE test_mod

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: array
  integer :: nsize

END MODULE test_mod

SUBROUTINE assign_array(x,n)
  USE test_mod
  IMPLICIT NONE

  INTEGER :: n
  DOUBLE PRECISION :: x(n)

  CALL test_allocated()
  array = x
  CALL test_allocated()

  nsize = n

END SUBROUTINE assign_array


SUBROUTINE print_array()
  USE test_mod, ONLY: nsize, array
  IMPLICIT NONE

  INTEGER :: i

  DO i = 1,nsize
     WRITE(6,'(F24.16)') array(i)
  END DO

END SUBROUTINE print_array

SUBROUTINE test_allocated()
  USE test_mod
  IMPLICIT NONE

  IF (ALLOCATED(array)) THEN
     WRITE(6,*) 'Array is allocated'
     WRITE(6,*) 'size is ', SIZE(array)
  ELSE
     WRITE(6,*) 'Array is NOT allocated'
  END IF

END SUBROUTINE test_allocated
