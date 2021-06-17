! =====================================================
SUBROUTINE clawpack46_rpn2_manifold(ixy,maxm,meqn,mwaves,mbc, &
     mx,ql,qr,auxl,auxr,wave,s,amdq,apdq)
! =====================================================

! Roe-solver for the 2D shallow water equations
!  on the sphere, using 3d Cartesian representation of velocities

! waves: 3
! equations: 4
! aux fields: 16

! Conserved quantities:
!       1 depth
!       2 x_momentum
!       3 y_momentum
!       4 z_momentum

! Auxiliary variables:
!         1  kappa
!         2  enx
!         3  eny
!         4  enz
!         5  etx
!         6  ety
!         7  etz
!         8  enx
!         9  eny
!        10  enz
!        11  etx
!        12  ety
!        13  etz
!        14  erx
!        15  ery
!        16  erz

! solve Riemann problems along one slice of data.

! On input, ql contains the state vector at the left edge of each cell
!           qr contains the state vector at the right edge of each cell

! This data is along a slice in the x-direction if ixy=1
!                            or the y-direction if ixy=2.
! On output, wave contains the waves, s the speeds,
! and amdq, apdq the decomposition of the flux difference
!   f(qr(i-1)) - f(ql(i))
! into leftgoing and rightgoing parts respectively.
! With the Roe solver we have
!    amdq  =  A^- \Delta q    and    apdq  =  A^+ \Delta q
! where A is the Roe matrix.  An entropy fix can also be incorporated
! into the flux differences.

! Note that the i'th Riemann problem has left state qr(i-1,:)
!                                    and right state ql(i,:)
! From the basic clawpack routines, this routine is called with ql = qr


    implicit none

    integer ixy, maxm,meqn,mwaves,mbc,mx
    double precision  wave(1-mbc:maxm+mbc,meqn,mwaves)
    double precision     s(1-mbc:maxm+mbc,mwaves)
    double precision    ql(1-mbc:maxm+mbc,meqn)
    double precision    qr(1-mbc:maxm+mbc,meqn)
    double precision  apdq(1-mbc:maxm+mbc,meqn)
    double precision  amdq(1-mbc:maxm+mbc,meqn)
    double precision  auxl(1-mbc:maxm+mbc,*)
    double precision  auxr(1-mbc:maxm+mbc,*)

!     local arrays
!     ------------
    double precision delta(3)
    logical :: efix

    integer maxm2
    parameter (maxm2 = 1800)
    double precision u(1-mbc:maxm+mbc),v(1-mbc:maxm+mbc),a(1-mbc:maxm+mbc)
    double precision h(1-mbc:maxm+mbc)

    double precision dtcom, dxcom, dycom, tcom
    integer icom, jcom
    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

    DOUBLE PRECISION grav
    COMMON /cparam/  grav

    integer i, m, mw, ioff
    double precision enx, eny, enz, etx,ety,etz
    double precision hunl, hunr, hutl, hutr, hl,hr,hsqr,hsql,hsq
    double precision a1,a2,a3, gamma, amn, apn, df, dy
    double precision erx, ery, erz, h1, h3, hi, him1, hu1, hu3
    double precision s0, s03, s1, s3, sfract


    data efix /.true./    !# use entropy fix for transonic rarefactions

    if(ixy == 1) then
        dy = dycom
    else
        dy = dxcom
    endif

!     The aux array has the following elements:
!         1  kappa = ratio of cell area to dxc*dyc
!         2  enx = x-component of normal vector to left edge in tangent plane
!         3  eny = y-component of normal vector to left edge in tangent plane
!         4  enz = z-component of normal vector to left edge in tangent plane
!         5  etx = x-component of tangent vector to left edge in tangent plane
!         6  ety = y-component of tangent vector to left edge in tangent plane
!         7  etz = z-component of tangent vector to left edge in tangent plane
!         8  enx = x-component of normal vector to bottom edge in tangent plane
!         9  eny = y-component of normal vector to bottom edge in tangent plane
!        10  enz = z-component of normal vector to bottom edge in tangent plane
!        11  etx = x-component of tangent vector to bottom edge in tangent plane
!        12  ety = y-component of tangent vector to bottom edge in tangent plane
!        13  etz = z-component of tangent vector to bottom edge in tangent plane
!        14  erx = x-component of unit vector in radial direction at cell ctr
!        15  ery = y-component of unit vector in radial direction at cell ctr
!        16  erz = z-component of unit vector in radial direction at cell ctr

!     # offset to index into aux array for enx, eny, etx, ety, gamma
!     #    depends on whether ixy=1 (left edge) or ixy=2 (bottom edge).
    ioff = 6*(ixy-1) + 1


!     # find a1 thru a3, the coefficients of the 3 eigenvectors:

    do i = 2-mbc, mx+mbc

        enx =   auxl(i,ioff+1)
        eny =   auxl(i,ioff+2)
        enz =   auxl(i,ioff+3)
        etx =   auxl(i,ioff+4)
        ety =   auxl(i,ioff+5)
        etz =   auxl(i,ioff+6)
        gamma = dsqrt(etx**2 + ety**2 + etz**2)
        etx =   etx / gamma
        ety =   ety / gamma
        etz =   etz / gamma


        !!  # compute normal and tangential momentum at cell edge:
        hunl = enx*ql(i,2)   + eny*ql(i,3)   + enz*ql(i,4)
        hunr = enx*qr(i-1,2) + eny*qr(i-1,3) + enz*qr(i-1,4)

        hutl = etx*ql(i,2)   + ety*ql(i,3)   + etz*ql(i,4)
        hutr = etx*qr(i-1,2) + ety*qr(i-1,3) + etz*qr(i-1,4)

        !! # compute the Roe-averaged variables needed in the Roe solver.
        !! # These are stored in the common block comroe since they are
        !! # later used in routine rpt2 to do the transverse wave splitting.

        hl = ql(i,1)
        hr = qr(i-1,1)
        h(i) = (hl+hr)*0.50d0
        hsqr = dsqrt(hr)
        hsql = dsqrt(hl)
        hsq = hsqr + hsql
        u(i) = (hunr/hsqr + hunl/hsql) / hsq
        v(i) = (hutr/hsqr + hutl/hsql) / hsq
        a(i) = dsqrt(grav*h(i))

        !! # Split the jump in q at each interface into waves
        delta(1) = hl - hr
        delta(2) = hunl - hunr
        delta(3) = hutl - hutr

        a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
        a2 = -v(i)*delta(1) + delta(3)
        a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))


        !! # Compute the waves.

        wave(i,1,1) = a1
        wave(i,2,1) = a1*(u(i)-a(i))*enx + a1*v(i)*etx
        wave(i,3,1) = a1*(u(i)-a(i))*eny + a1*v(i)*ety
        wave(i,4,1) = a1*(u(i)-a(i))*enz + a1*v(i)*etz
        s(i,1) = (u(i)-a(i)) * gamma/dy

        wave(i,1,2) = 0.0d0
        wave(i,2,2) = a2*etx
        wave(i,3,2) = a2*ety
        wave(i,4,2) = a2*etz
        s(i,2) = u(i) * gamma/dy

        wave(i,1,3) = a3
        wave(i,2,3) = a3*(u(i)+a(i))*enx + a3*v(i)*etx
        wave(i,3,3) = a3*(u(i)+a(i))*eny + a3*v(i)*ety
        wave(i,4,3) = a3*(u(i)+a(i))*enz + a3*v(i)*etz
        s(i,3) = (u(i)+a(i)) * gamma/dy
!!        281 format(2i4,5d12.4)
!!        283 format(8x,5d12.4)
    END DO


    !! # compute flux differences amdq and apdq.
    !! ---------------------------------------

    if (efix) go to 110

    !! # no entropy fix
    !! ----------------

    !! # amdq = SUM s*wave   over left-going waves
    !! # apdq = SUM s*wave   over right-going waves

    DO m=1,meqn
        do i=2-mbc, mx+mbc
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do mw=1,mwaves
                if (s(i,mw) < 0.d0) then
                    amdq(i,m) = amdq(i,m) + s(i,mw)*wave(i,m,mw)
                else
                    apdq(i,m) = apdq(i,m) + s(i,mw)*wave(i,m,mw)
                endif
            END DO
        end do
    END DO

    !! # project momentum components of amdq and apdq onto tangent plane:

    do i=2-mbc,mx+mbc
        erx = auxr(14,i-1)
        ery = auxr(15,i-1)
        erz = auxr(16,i-1)
        amn = erx*amdq(i,2) + ery*amdq(i,3) + erz*amdq(i,4)
        amdq(i,2) = amdq(i,2) - amn*erx
        amdq(i,3) = amdq(i,3) - amn*ery
        amdq(i,4) = amdq(i,4) - amn*erz

        erx = auxl(i,14)
        ery = auxl(i,15)
        erz = auxl(i,16)
        apn = erx*apdq(i,2) + ery*apdq(i,3) + erz*apdq(i,4)
        apdq(i,2) = apdq(i,2) - apn*erx
        apdq(i,3) = apdq(i,3) - apn*ery
        apdq(i,4) = apdq(i,4) - apn*erz
    enddo
    go to 900

!-----------------------------------------------------

    110 continue

!     # With entropy fix
!     ------------------

!    # compute flux differences amdq and apdq.
!    # First compute amdq as sum of s*wave for left going waves.
!    # Incorporate entropy fix by adding a modified fraction of wave
!    # if s should change sign.

    do 200 i=2-mbc,mx+mbc
        do m=1, meqn
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
        enddo
        enx =   auxl(i,ioff+1)
        eny =   auxl(i,ioff+2)
        enz =   auxl(i,ioff+3)
        etx =   auxl(i,ioff+4)
        ety =   auxl(i,ioff+5)
        etz =   auxl(i,ioff+6)
        gamma = dsqrt(etx**2 + ety**2 + etz**2)
        etx =   etx / gamma
        ety =   ety / gamma
        etz =   etz / gamma
    !           # compute normal and tangential momentum at cell edge:
        hunl = enx*ql(i,2)   + eny*ql(i,3) + enz*ql(i,4)
        hunr = enx*qr(i-1,2) + eny*qr(i-1,3) + enz*qr(i-1,4)

    !           check 1-wave
        him1 = qr(i-1,1)
        s0 =  (hunr/him1 - dsqrt(grav*him1)) * gamma / dy
    !           check for fully supersonic case :
        if (s0 > 0.0d0 .AND. s(i,1) > 0.0d0) then
            do 60 m=1,4
                amdq(i,m)=0.0d0
            60 END DO
            goto 200
        endif

        h1 = qr(i-1,1)+wave(i,1,1)
        hu1= hunr + enx*wave(i,2,1) + eny*wave(i,3,1) &
                       + enz*wave(i,4,1)
        s1 = (hu1/h1 - dsqrt(grav*h1))*gamma/dy  !speed just to right of 1-wave
        if (s0 < 0.0d0 .AND. s1 > 0.0d0) then
        !              transonic rarefaction in 1-wave
            sfract = s0*((s1-s(i,1))/(s1-s0))
        else if (s(i,1) < 0.0d0) then
        !              1-wave is leftgoing
            sfract = s(i,1)
        else
        !              1-wave is rightgoing
            sfract = 0.0d0
        endif
        do 120 m=1,4
            amdq(i,m) = sfract*wave(i,m,1)
        120 END DO
    !           check 2-wave
        if (s(i,2) > 0.0d0) then
        !	       #2 and 3 waves are right-going
            go to 200
        endif

        do 140 m=1,4
            amdq(i,m) = amdq(i,m) + s(i,2)*wave(i,m,2)
        140 END DO

    !           check 3-wave

        hi = ql(i,1)
        s03 = (hunl/hi + dsqrt(grav*hi)) * gamma/dy
        h3=ql(i,1)-wave(i,1,3)
        hu3=hunl-(enx*wave(i,2,3)+eny*wave(i,3,3)+enz*wave(i,4,3))
        s3=(hu3/h3 + dsqrt(grav*h3)) * gamma/dy
        if (s3 < 0.0d0 .AND. s03 > 0.0d0) then
        !              transonic rarefaction in 3-wave
            sfract = s3*((s03-s(i,3))/(s03-s3))
        else if (s(i,3) < 0.0d0) then
        !              3-wave is leftgoing
            sfract = s(i,3)
        else
        !              3-wave is rightgoing
            goto 200
        endif
        do 160 m=1,4
            amdq(i,m) = amdq(i,m) + sfract*wave(i,m,3)
        160 END DO
    200 END DO

!           compute rightgoing flux differences :

    do 220 i = 2-mbc,mx+mbc
        do 222 m=1,4
            df = 0.0d0
            do 210 mw=1,mwaves
                df = df + s(i,mw)*wave(i,m,mw)
            210 END DO
            apdq(i,m)=df - amdq(i,m)
        222 END DO

    !                 project momentum components onto tangent plane

        erx = auxr(i-1,14)
        ery = auxr(i-1,15)
        erz = auxr(i-1,16)
        amn = erx*amdq(i,2)+ery*amdq(i,3)+erz*amdq(i,4)
        amdq(i,2) = amdq(i,2) - amn*erx
        amdq(i,3) = amdq(i,3) - amn*ery
        amdq(i,4) = amdq(i,4) - amn*erz

        erx = auxl(i,14)
        ery = auxl(i,15)
        erz = auxl(i,16)
        apn = erx*apdq(i,2)+ery*apdq(i,3)+erz*apdq(i,4)
        apdq(i,2) = apdq(i,2) - apn*erx
        apdq(i,3) = apdq(i,3) - apn*ery
        apdq(i,4) = apdq(i,4) - apn*erz

    220 END DO


    900 continue
    return
  END SUBROUTINE clawpack46_rpn2_manifold
