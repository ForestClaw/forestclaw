! ==============================================================================
!  Regions Module
!   Module containing data structures and setup routines for region refinement.
!
! ==============================================================================
module regions_module

    implicit none
    save

    ! Region type definition
    type region_type
        integer :: min_level,max_level
        real(kind=8) :: x_low,y_low,x_hi,y_hi,t_low,t_hi
    end type region_type

    logical, private :: module_setup

    integer :: num_regions
    type(region_type), allocatable :: regions(:)
      
contains

    subroutine set_regions(fname)

        use amr_module
      
        implicit none
      
        ! Function Arguments
        character(len=*), optional, intent(in) :: fname
      
        ! Locals
        integer, parameter :: unit = 7
        integer :: i

        if (.not. module_setup) then

            write(parmunit,*) ' '
            write(parmunit,*) '--------------------------------------------'
            write(parmunit,*) 'REGIONS:'
            write(parmunit,*) '-----------'

            if (present(fname)) then
                call opendatafile(unit,fname)
            else
                call opendatafile(unit,'regions.data')
            endif

            read(unit,"(i2)") num_regions
            if (num_regions == 0) then
                write(parmunit,*) '  No regions specified for refinement'
                
            else
                ! Refinement region data
                allocate(regions(num_regions))
                do i=1,num_regions
                    read(unit,*) regions(i)%min_level, regions(i)%max_level, &
                                 regions(i)%t_low, regions(i)%t_hi, &
                                 regions(i)%x_low, regions(i)%x_hi, &
                                 regions(i)%y_low, regions(i)%y_hi
                enddo
            endif
            close(unit)

            module_setup = .true.
        end if

    end subroutine set_regions

end module regions_module
