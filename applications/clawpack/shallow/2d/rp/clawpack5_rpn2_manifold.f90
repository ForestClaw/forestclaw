! =====================================================
SUBROUTINE clawpack5_rpn2_manifold(ixy,maxm,meqn,mwaves,maux,mbc, &
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


    implicit double precision (a-h,o-z)

    dimension wave(meqn, mwaves,1-mbc:maxm+mbc)
    dimension    s(mwaves, 1-mbc:maxm+mbc)
    dimension   ql(meqn, 1-mbc:maxm+mbc)
    dimension   qr(meqn, 1-mbc:maxm+mbc)
    dimension  apdq(meqn, 1-mbc:maxm+mbc)
    dimension  amdq(meqn, 1-mbc:maxm+mbc)
    dimension auxl(maux, 1-mbc:maxm+mbc)
    dimension auxr(maux, 1-mbc:maxm+mbc)

!     local arrays
!     ------------
    dimension delta(3)
    logical :: efix

    parameter (maxm2 = 1800)
    dimension u(1-mbc:maxm+mbc),v(1-mbc:maxm+mbc),a(1-mbc:maxm+mbc), &
    h(1-mbc:maxm+mbc)
    common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

    DOUBLE PRECISION grav
    COMMON /cparam/  grav


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

        enx =   auxl(ioff+1,i)
        eny =   auxl(ioff+2,i)
        enz =   auxl(ioff+3,i)
        etx =   auxl(ioff+4,i)
        ety =   auxl(ioff+5,i)
        etz =   auxl(ioff+6,i)
        gamma = dsqrt(etx**2 + ety**2 + etz**2)
        etx =   etx / gamma
        ety =   ety / gamma
        etz =   etz / gamma


        g = grav

    !        # compute normal and tangential momentum at cell edge:
        hunl = enx*ql(2,i) + eny*ql(3,i) + enz*ql(4,i)
        hunr = enx*qr(2,i-1) + eny*qr(3,i-1) + enz*qr(4,i-1)

        hutl = etx*ql(2,i) + ety*ql(3,i) + etz*ql(4,i)
        hutr = etx*qr(2,i-1) + ety*qr(3,i-1) + etz*qr(4,i-1)

    !        # compute the Roe-averaged variables needed in the Roe solver.
    !        # These are stored in the common block comroe since they are
    !        # later used in routine rpt2 to do the transverse wave splitting.

        hl = ql(1,i)
        hr = qr(1,i-1)
        h(i) = (hl+hr)*0.50d0
        hsqr = dsqrt(hr)
        hsql = dsqrt(hl)
        hsq = hsqr + hsql
        u(i) = (hunr/hsqr + hunl/hsql) / hsq
        v(i) = (hutr/hsqr + hutl/hsql) / hsq
        a(i) = dsqrt(grav*h(i))

    !        # Split the jump in q at each interface into waves
        delta(1) = hl - hr
        delta(2) = hunl - hunr
        delta(3) = hutl - hutr

        a1 = ((u(i)+a(i))*delta(1) - delta(2))*(0.50d0/a(i))
        a2 = -v(i)*delta(1) + delta(3)
        a3 = (-(u(i)-a(i))*delta(1) + delta(2))*(0.50d0/a(i))


    !        # Compute the waves.

        wave(1,1,i) = a1
        wave(2,1,i) = a1*(u(i)-a(i))*enx + a1*v(i)*etx
        wave(3,1,i) = a1*(u(i)-a(i))*eny + a1*v(i)*ety
        wave(4,1,i) = a1*(u(i)-a(i))*enz + a1*v(i)*etz
        s(1,i) = (u(i)-a(i)) * gamma/dy

        wave(1,2,i) = 0.0d0
        wave(2,2,i) = a2*etx
        wave(3,2,i) = a2*ety
        wave(4,2,i) = a2*etz
        s(2,i) = u(i) * gamma/dy

        wave(1,3,i) = a3
        wave(2,3,i) = a3*(u(i)+a(i))*enx + a3*v(i)*etx
        wave(3,3,i) = a3*(u(i)+a(i))*eny + a3*v(i)*ety
        wave(4,3,i) = a3*(u(i)+a(i))*enz + a3*v(i)*etz
        s(3,i) = (u(i)+a(i)) * gamma/dy
!!        281 format(2i4,5d12.4)
!!        283 format(8x,5d12.4)
    END DO


!    # compute flux differences amdq and apdq.
!    ---------------------------------------

    if (efix) go to 110

!     # no entropy fix
!     ----------------

!     # amdq = SUM s*wave   over left-going waves
!     # apdq = SUM s*wave   over right-going waves

    DO m=1,meqn
        do i=2-mbc, mx+mbc
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
            do mw=1,mwaves
                if (s(mw,i) < 0.d0) then
                    amdq(m,i) = amdq(m,i) + s(mw,i)*wave(m,mw,i)
                else
                    apdq(m,i) = apdq(m,i) + s(mw,i)*wave(m,mw,i)
                endif
            END DO
        end do
    END DO

!     # project momentum components of amdq and apdq onto tangent plane:

    do i=2-mbc,mx+mbc
        erx = auxr(14,i-1)
        ery = auxr(15,i-1)
        erz = auxr(16,i-1)
        amn = erx*amdq(2,i) + ery*amdq(3,i) + erz*amdq(4,i)
        amdq(2,i) = amdq(2,i) - amn*erx
        amdq(3,i) = amdq(3,i) - amn*ery
        amdq(4,i) = amdq(4,i) - amn*erz

        erx = auxl(14,i)
        ery = auxl(15,i)
        erz = auxl(16,i)
        apn = erx*apdq(2,i) + ery*apdq(3,i) + erz*apdq(4,i)
        apdq(2,i) = apdq(2,i) - apn*erx
        apdq(3,i) = apdq(3,i) - apn*ery
        apdq(4,i) = apdq(4,i) - apn*erz
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

    do i=2-mbc,mx+mbc
        do m=1, meqn
            amdq(m,i) = 0.d0
            apdq(m,i) = 0.d0
        enddo
        enx =   auxl(ioff+1,i)
        eny =   auxl(ioff+2,i)
        enz =   auxl(ioff+3,i)
        etx =   auxl(ioff+4,i)
        ety =   auxl(ioff+5,i)
        etz =   auxl(ioff+6,i)
        gamma = dsqrt(etx**2 + ety**2 + etz**2)
        etx =   etx / gamma
        ety =   ety / gamma
        etz =   etz / gamma
    !           # compute normal and tangential momentum at cell edge:
        hunl = enx*ql(2,i) + eny*ql(3,i) + enz*ql(4,i)
        hunr = enx*qr(2,i-1) + eny*qr(3,i-1) + enz*qr(4,i-1)

    !           check 1-wave
        him1 = qr(1,i-1)
        s0 =  (hunr/him1 - dsqrt(g*him1)) * gamma / dy
    !           check for fully supersonic case :
        if (s0 > 0.0d0 .AND. s(1,i) > 0.0d0) then
            do m=1,4
                amdq(m,i)=0.0d0
            END DO
            goto 200
        endif

        h1 = qr(1,i-1)+wave(1,1,i)
        hu1= hunr + enx*wave(2,1,i) + eny*wave(3,1,i) &
        + enz*wave(4,1,i)
        s1 = (hu1/h1 - dsqrt(g*h1))*gamma/dy  !speed just to right of 1-wave
        if (s0 < 0.0d0 .AND. s1 > 0.0d0) then
        !              transonic rarefaction in 1-wave
            sfract = s0*((s1-s(1,i))/(s1-s0))
        else if (s(1,i) < 0.0d0) then
        !              1-wave is leftgoing
            sfract = s(1,i)
        else
        !              1-wave is rightgoing
            sfract = 0.0d0
        endif
        do m=1,4
            amdq(m,i) = sfract*wave(m,1,i)
        END DO
    !           check 2-wave
        if (s(2,i) > 0.0d0) then
        !	       #2 and 3 waves are right-going
            go to 200
        endif

        do m=1,4
            amdq(m,i) = amdq(m,i) + s(2,i)*wave(m,2,i)
        END DO

    !           check 3-wave

        hi = ql(1,i)
        s03 = (hunl/hi + dsqrt(g*hi)) * gamma/dy
        h3=ql(1,i)-wave(1,3,i)
        hu3=hunl-(enx*wave(2,3,i)+eny*wave(3,3,i)+enz*wave(4,3,i))
        s3=(hu3/h3 + dsqrt(g*h3)) * gamma/dy
        if (s3 < 0.0d0 .AND. s03 > 0.0d0) then
        !              transonic rarefaction in 3-wave
            sfract = s3*((s03-s(3,i))/(s03-s3))
        else if (s(3,i) < 0.0d0) then
        !              3-wave is leftgoing
            sfract = s(3,i)
        else
        !              3-wave is rightgoing
            goto 200
        endif
        do m=1,4
            amdq(m,i) = amdq(m,i) + sfract*wave(m,3,i)
        END DO
    200 END DO

!           compute rightgoing flux differences :

    do i = 2-mbc,mx+mbc
        do m=1,4
            df = 0.0d0
            do mw=1,mwaves
                df = df + s(mw,i)*wave(m,mw,i)
            END DO
            apdq(m,i)=df - amdq(m,i)
        END DO

    !                 project momentum components onto tangent plane

        erx = auxr(14,i-1)
        ery = auxr(15,i-1)
        erz = auxr(16,i-1)
        amn = erx*amdq(2,i)+ery*amdq(3,i)+erz*amdq(4,i)
        amdq(2,i) = amdq(2,i) - amn*erx
        amdq(3,i) = amdq(3,i) - amn*ery
        amdq(4,i) = amdq(4,i) - amn*erz

        erx = auxl(14,i)
        ery = auxl(15,i)
        erz = auxl(16,i)
        apn = erx*apdq(2,i)+ery*apdq(3,i)+erz*apdq(4,i)
        apdq(2,i) = apdq(2,i) - apn*erx
        apdq(3,i) = apdq(3,i) - apn*ery
        apdq(4,i) = apdq(4,i) - apn*erz

    END DO


    900 continue
    return
  END SUBROUTINE clawpack5_rpn2_manifold
