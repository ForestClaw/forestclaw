integer function geoclaw_gauges_getnum(fname)

    implicit none

    ! Input
    character(len=20), intent(in), optional :: fname
    integer :: num_gauges

    ! Locals
    integer, parameter :: iunit = 7

    write(*,*) "writing fname"
    write(*,*) fname
    ! Open file
    if (present(fname)) then
      call opendatafile(iunit,fname)
    else
      call opendatafile(iunit,'gauges.data')
    endif

    read(iunit,*) num_gauges
    geoclaw_gauges_getnum = num_gauges
    close(iunit)
end function

SUBROUTINE geoclaw_gauges_init(restart, meqn, num_gauges, gauges, fname)

    ! use amr_module
    implicit none
    type gauge_type
      integer:: blockno
      double precision :: xc, yc, t1, t2
      integer num;
      double precision, pointer :: buffer
    end type gauge_type

    ! Input
    character(len=20), intent(in), optional :: fname
    logical, intent(in)  :: restart
    integer, intent(in) :: meqn, num_gauges
    integer :: num_gauges_not_use

    type(gauge_type), dimension(num_gauges) :: gauges

    ! Locals
    integer :: i, ipos, idigit, OUTGAUGEUNIT
    integer, parameter :: iunit = 7
    character*14 ::  fileName

    ! Open file
    if (present(fname)) then
      call opendatafile(iunit,fname)
    else
      call opendatafile(iunit,'gauges.data')
    endif

    read(iunit,*) num_gauges_not_use

    ! allocate(xgauge(num_gauges), ygauge(num_gauges))
    ! allocate(t1gauge(num_gauges), t2gauge(num_gauges))
    ! allocate(mbestsrc(num_gauges), mbestorder(num_gauges))
    ! allocate(igauge(num_gauges))
    ! allocate(mbestg1(maxgr), mbestg2(maxgr))

    ! allocate(nextLoc(num_gauges))
    ! allocate(gaugeArray(nvar + 2,MAXDATA,num_gauges))
    ! allocate(levelArray(MAXDATA,num_gauges))

    do i=1,num_gauges
      read(iunit,*) gauges(i)%num,gauges(i)%xc,gauges(i)%yc,gauges(i)%t1,gauges(i)%t2
    enddo

    close(iunit)


    do i = 1, num_gauges
     fileName = 'gaugexxxxx.txt'    ! NB different name convention too
     ! inum = igauge(i)
     ! do ipos = 10,6,-1              ! do this to replace the xxxxx in the name
      ! idigit = mod(inum,10)
      ! fileName(ipos:ipos) = char(ichar('0') + idigit)
      ! inum = inum / 10
    ! end do

    !          status unknown since might be a restart run. maybe need to test and rewind?
    if (restart) then
      open(unit=OUTGAUGEUNIT, file=fileName, status='old',        &
       position='append', form='formatted')
    else
      open(unit=OUTGAUGEUNIT, file=fileName, status='unknown',        &
       position='append', form='formatted')
      rewind OUTGAUGEUNIT
      ! write(OUTGAUGEUNIT,100) igauge(i), xgauge(i), ygauge(i), 4
      ! 100              format("# gauge_id= ",i5," location=( ",1e15.7," ",1e15.7," ) num_eqn= ",i2)
      write(OUTGAUGEUNIT,101)
      101              format("# Columns: level time h    hu    hv    eta")
    endif

    close(OUTGAUGEUNIT)

end do

END SUBROUTINE geoclaw_gauges_init
