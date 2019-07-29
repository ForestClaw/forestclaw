SUBROUTINE cudaclaw5_step2(maxm,meqn,maux,mbc,mx,my,qold,aux,dx,dy,dt, &
     cflgrid,fm,fp,gm,gp,rpn2,rpt2,block_corner_count)
!
!     clawpack routine ...  modified for AMRCLAW
!
!     Take one time step, updating q.
!     On entry, qold gives
!        initial data for this step
!        and is unchanged in this version.
!
!     fm, fp are fluxes to left and right of single cell edge
!     See the flux2 documentation for more information.
!
!     Converted to f90 2012-1-04 (KTM)
!

    use cudaclaw5_amr_module

    implicit none

    external rpn2, rpt2

    ! Arguments
    integer, intent(in) :: maxm,meqn,maux,mbc,mx,my
    real(kind=8), intent(in) :: dx,dy,dt
    real(kind=8), intent(inout) :: cflgrid
    real(kind=8), intent(inout) :: qold(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fm(meqn, 1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: fp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gm(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
    real(kind=8), intent(inout) :: gp(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)

    ! Local storage for flux accumulation
    real(kind=8) :: faddm(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: faddp(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: gaddm(meqn,1-mbc:maxm+mbc,2)
    real(kind=8) :: gaddp(meqn,1-mbc:maxm+mbc,2)

    ! Scratch storage for Sweeps and Riemann problems
    real(kind=8) ::  q1d(meqn,1-mbc:maxm+mbc)
    real(kind=8) :: aux1(maux,1-mbc:maxm+mbc)
    real(kind=8) :: aux2(maux,1-mbc:maxm+mbc)
    real(kind=8) :: aux3(maux,1-mbc:maxm+mbc)
    real(kind=8) :: dtdx1d(1-mbc:maxm+mbc)
    real(kind=8) :: dtdy1d(1-mbc:maxm+mbc)

    real(kind=8) ::  wave(meqn, mwaves, 1-mbc:maxm+mbc)
    real(kind=8) ::     s(mwaves, 1-mbc:maxm + mbc)
    real(kind=8) ::  amdq(meqn,1-mbc:maxm + mbc)
    real(kind=8) ::  apdq(meqn,1-mbc:maxm + mbc)
    real(kind=8) ::  cqxx(meqn,1-mbc:maxm + mbc)
    real(kind=8) :: bmadq(meqn,1-mbc:maxm + mbc)
    REAL(kind=8) :: bpadq(meqn,1-mbc:maxm + mbc)

    INTEGER block_corner_count(0:3), ibc,jbc,m


    ! Looping scalar storage
    integer :: i,j,thread_num
    real(kind=8) :: dtdx,dtdy,cfl1d

    ! Common block storage
    integer :: icom,jcom
    real(kind=8) :: dtcom,dxcom,dycom,tcom
    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

    ! Store mesh parameters in common block
    dxcom = dx
    dycom = dy
    dtcom = dt

    cflgrid = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy

    fm = 0.d0
    fp = 0.d0
    gm = 0.d0
    gp = 0.d0

    ! ============================================================================
    ! Perform X-Sweeps
    DO j = 0,my+1
       ! Copy old q into 1d slice
       q1d(:,1-mbc:mx+mbc) = qold(:,1-mbc:mx+mbc,j)

       IF (j .EQ. 0) THEN
          DO m = 1,meqn
             IF (block_corner_count(0) .EQ. 3) THEN
                DO ibc = 1,mbc
                   q1d(m,1-ibc) = qold(m,ibc,0)
                ENDDO
             ENDIF
             IF (block_corner_count(1) .EQ. 3) THEN
                DO ibc = 1,mbc
                   q1d(m,mx+ibc) = qold(m,mx-ibc+1,0)
                ENDDO
             ENDIF
          ENDDO
       ELSE IF (j .EQ. my+1) THEN
          DO m = 1,meqn
             IF (block_corner_count(2) .EQ. 3) THEN
                DO ibc = 1,mbc
                   q1d(m,1-ibc) = qold(m,ibc,my+1)
                ENDDO
             ENDIF
             IF (block_corner_count(3) .EQ. 3) THEN
                DO ibc = 1,mbc
                   q1d(m,mx+ibc) = qold(m,mx-ibc+1,my+1)
                ENDDO
             ENDIF
          ENDDO
       ENDIF



        ! Set dtdx slice if a capacity array exists
        if (mcapa > 0)  then
            dtdx1d(1-mbc:mx+mbc) = dtdx / aux(mcapa,1-mbc:mx+mbc,j)
        else
            dtdx1d = dtdx
        endif

        ! Copy aux array into slices
        if (maux > 0) then
            aux1(:,1-mbc:mx+mbc) = aux(:,1-mbc:mx+mbc,j-1)
            aux2(:,1-mbc:mx+mbc) = aux(:,1-mbc:mx+mbc,j  )
            aux3(:,1-mbc:mx+mbc) = aux(:,1-mbc:mx+mbc,j+1)
        endif

        ! Store value of j along the slice into common block
        jcom = j

        ! Compute modifications fadd and gadd to fluxes along this slice:
        !!! Change from flux2 to cudaclaw5_flux2
        call cudaclaw5_flux2(1,maxm,meqn,maux,mbc,mx,q1d,dtdx1d,aux1,aux2,aux3, &
                   faddm,faddp,gaddm,gaddp,cfl1d,wave,s, &
                   amdq,apdq,cqxx,bmadq,bpadq,rpn2,rpt2)
        cflgrid = max(cflgrid,cfl1d)

        ! Update fluxes
        fm(:,1:mx+1,j) = fm(:,1:mx+1,j) + faddm(:,1:mx+1)
        fp(:,1:mx+1,j) = fp(:,1:mx+1,j) + faddp(:,1:mx+1)
        gm(:,1:mx+1,j) = gm(:,1:mx+1,j) + gaddm(:,1:mx+1,1)
        gp(:,1:mx+1,j) = gp(:,1:mx+1,j) + gaddp(:,1:mx+1,1)
        gm(:,1:mx+1,j+1) = gm(:,1:mx+1,j+1) + gaddm(:,1:mx+1,2)
        gp(:,1:mx+1,j+1) = gp(:,1:mx+1,j+1) + gaddp(:,1:mx+1,2)

    enddo

    ! ============================================================================
    !  y-sweeps
    !
    do i = 0,mx+1

       ! Copy data along a slice into 1d arrays:
       q1d(:,1-mbc:my+mbc) = qold(:,i,1-mbc:my+mbc)

       IF (i .EQ. 0) THEN
          DO m = 1,meqn
             IF (block_corner_count(0) .EQ. 3) THEN
                DO jbc = 1,mbc
                   q1d(m,1-jbc) = qold(m,0,jbc)
                ENDDO
             ENDIF
             IF (block_corner_count(2) .EQ. 3) THEN
                DO jbc = 1,mbc
                   q1d(m,my+jbc) = qold(m,0,my-jbc+1)
                ENDDO
             ENDIF
          ENDDO
       ELSE IF (i .EQ. mx+1) THEN
          DO m = 1,meqn
             IF (block_corner_count(1) .EQ. 3) THEN
                DO jbc = 1,mbc
                   q1d(m,1-jbc) = qold(m,mx+1,jbc)
                ENDDO
             ENDIF
             IF (block_corner_count(3) .EQ. 3) THEN
                DO jbc= 1,mbc
                   q1d(m,my+jbc) = qold(m,mx+1,my-jbc+1)
                ENDDO
             ENDIF
          ENDDO
       ENDIF



        ! Set dt/dy ratio in slice
        if (mcapa > 0) then
            dtdy1d(1-mbc:my+mbc) = dtdy / aux(mcapa,i,1-mbc:my+mbc)
        else
            dtdy1d = dtdy
        endif

        ! Copy aux slices
        if (maux .gt. 0)  then
            aux1(:,1-mbc:my+mbc) = aux(:,i-1,1-mbc:my+mbc)
            aux2(:,1-mbc:my+mbc) = aux(:,i,1-mbc:my+mbc)
            aux3(:,1-mbc:my+mbc) = aux(:,i+1,1-mbc:my+mbc)
        endif

        ! Store the value of i along this slice in the common block
        icom = i

        ! Compute modifications fadd and gadd to fluxes along this slice
        !!! Change from flux2 to cudaclaw5_flux2
        call cudaclaw5_flux2(2,maxm,meqn,maux,mbc,my,q1d,dtdy1d,aux1,aux2,aux3, &
                   faddm,faddp,gaddm,gaddp,cfl1d,wave,s, &
                   amdq,apdq,cqxx,bmadq,bpadq,rpn2,rpt2)
        cflgrid = max(cflgrid,cfl1d)

        ! Update fluxes
        gm(:,i,1:my+1) = gm(:,i,1:my+1) + faddm(:,1:my+1)
        gp(:,i,1:my+1) = gp(:,i,1:my+1) + faddp(:,1:my+1)
        fm(:,i,1:my+1) = fm(:,i,1:my+1) + gaddm(:,1:my+1,1)
        fp(:,i,1:my+1) = fp(:,i,1:my+1) + gaddp(:,1:my+1,1)
        fm(:,i+1,1:my+1) = fm(:,i+1,1:my+1) + gaddm(:,1:my+1,2)
        fp(:,i+1,1:my+1) = fp(:,i+1,1:my+1) + gaddp(:,1:my+1,2)

    enddo


end subroutine cudaclaw5_step2
