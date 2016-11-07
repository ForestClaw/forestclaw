integer function geoclaw_gauges_getnum(fname)

    implicit none

    !! Input
    !! NOTE : Length specified must match exact length of file name.
    CHARACTER(len=12), INTENT(in), OPTIONAL :: fname
    integer :: num_gauges

    ! Locals
    integer, parameter :: iunit = 7

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
      integer:: patchno
      integer:: location_in_results
      double precision :: xc, yc, t1, t2
      integer num;
      !! DOUBLE PRECISION, POINTER :: buffer
    end type gauge_type

    ! Input
    character(len=12), intent(in), optional :: fname
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
       READ(iunit,*) gauges(i)%num,gauges(i)%xc,gauges(i)%yc,gauges(i)%t1,gauges(i)%t2
       gauges(i)%location_in_results = -1
       gauges(i)%blockno = -1
       gauges(i)%patchno = -1
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

SUBROUTINE geoclaw_update_gauge (mx,my,mbc,meqn,xlower,ylower,dx,dy,q,maux,aux,&
                                 xc,yc,var,eta)

      use geoclaw_module, only: dry_tolerance

      implicit none

      integer :: mx, my, mbc, meqn, maux
      real(kind=8) :: xlower, ylower, dx, dy, xc, yc
      real(kind=8) :: q(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
      real(kind=8) :: aux(maux, 1-mbc:mx+mbc, 1-mbc:my+mbc)
      real(kind=8) :: var(meqn), eta


      ! local variables:
      real(kind=8) :: xcent,ycent,xoff,yoff
      integer :: ioff,joff,iindex,jindex,mq
      real(kind=8) :: h(4),drytol2,topo
      ! integer :: index


      iindex =  int((xc-xlower)/dx) + 1
      jindex =  int((yc-ylower)/dy) + 1

      xcent  = xlower + (iindex-.5d0)*dx
      ycent  = ylower + (jindex-.5d0)*dy
      xoff   = (xc-xcent)/dx
      yoff   = (yc-ycent)/dy

      drytol2 = 0.1d0 * dry_tolerance

      h(1) = q(1,iindex,jindex)
      h(2) = q(1,iindex+1,jindex)
      h(3) = q(1,iindex,jindex+1)
      h(4) = q(1,iindex+1,jindex+1)


      if ((h(1) < drytol2) .or.  &
          (h(2) < drytol2) .or.  &
          (h(3) < drytol2) .or.  &
          (h(4) < drytol2)) then
          !! One of the cells is dry, so just use value from grid cell
          !! that contains gauge rather than interpolating

          ! icell = int(1.d0 + (xc - xlower) / hx)
          ! jcell = int(1.d0 + (yc - ylower) / hy)
          do mq=1,3
              var(mq) = q(mq,iindex,jindex)
          enddo
          !! This is the bottom layer and we should figure out the
          !! topography
          topo = aux(1,iindex,jindex)
      else
          !! Linear interpolation between four cells
          do mq=1,3
              var(mq) = (1.d0 - xoff) * (1.d0 - yoff) &
                        * q(mq,iindex,jindex)  &
                        + xoff*(1.d0 - yoff) * q(mq,iindex+1,jindex)  &
                        + (1.d0 - xoff) * yoff * q(mq,iindex,jindex+1)  &
                        + xoff * yoff * q(mq,iindex+1,jindex+1)
          enddo
          topo = (1.d0 - xoff) * (1.d0 - yoff)  &
                  * aux(1,iindex,jindex)  &
                  + xoff * (1.d0 - yoff) * aux(1,iindex+1,jindex)  &
                  + (1.d0 - xoff) * yoff * aux(1,iindex,jindex+1)  &
                  + xoff * yoff * aux(1,iindex+1,jindex+1)
        endif

        ! Extract surfaces
        eta = var(1) + topo

        ! Zero out tiny values to prevent later problems reading data,
        ! as done in valout.f
        do mq = 1,3
           if (abs(var(mq)) < 1d-90) var(mq) = 0.d0
           end do
        if (abs(eta) < 1d-90) eta = 0.d0

        ! ! save info for this time
        ! index = nextLoc(ii)

        ! levelArray(index,ii) = level
        ! gaugeArray(1,index,ii) = tgrid
        ! gaugeArray(2,index,ii) = var(1)
        ! gaugeArray(3,index,ii) = var(2)
        ! gaugeArray(4,index,ii) = var(3)
        ! gaugeArray(5,index,ii) = eta

        ! nextLoc(ii) = nextLoc(ii) + 1
        ! if (nextLoc(ii) .gt. MAXDATA) then
        !   call print_gauges_and_reset_nextLoc(ii, nvar)
        ! endif

END SUBROUTINE geoclaw_update_gauge
