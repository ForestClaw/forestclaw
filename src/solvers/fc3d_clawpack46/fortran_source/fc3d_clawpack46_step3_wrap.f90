subroutine clawpack46_step3_wrap(maxm, meqn, maux, mbc, &
      method, mthlim, mcapa, mwaves, mx, my, mz, qold, aux, &
      dx, dy, dz, dt,cfl, work, mwork,xlower,ylower,zlower, &
      level, t, fp,fm, gp, gm, hm, hp, rpn3, rpt3, rptt3, &
      use_fwaves, ierror)

    implicit none

    external :: rpn3,rpt3, rptt3, flux3

    INTEGER :: maxm,meqn,maux,mbc,mcapa,mwaves,mx,my, mz, mwork
    INTEGER :: level, ierror, use_fwaves
    INTEGER :: method(7), mthlim(mwaves)
!!      integer block_corner_count(0:3)

    DOUBLE PRECISION :: dx,dy,dz, dt,cfl, xlower, ylower, zlower, t
    DOUBLE PRECISION :: work(mwork)

    DOUBLE PRECISION :: qold(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION ::  aux(1-mbc:mx+mbc, 1-mbc:my+mbc, 1-mbc:mz+mbc,maux)

    DOUBLE PRECISION :: fp(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION :: fm(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION :: gp(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION :: gm(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION :: hp(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    DOUBLE PRECISION :: hm(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)


    !!  # Local variables
    integer :: i0faddm, i0faddp, i0gadd, i0hadd
    integer :: i0q1d, i0dtdx1, i0dtdy1, i0dtdz1
    integer :: i0aux1, i0aux2, i0aux3, i0next, mused, mwork1
    integer :: i0wave, i0s, i0amdq, i0apdq, i0ql, i0qr
    integer :: i0auxl, i0auxr

    double precision :: dtdx, dtdy, dtdz
    integer :: i,j,k, m

    !! Needed by Riemann solvers.  This should be fixed later by a 'context'
    !! for a Riemann solver.
    double precision :: dtcom, dxcom,dycom,dzcom, tcom
    integer :: icom, jcom, kcom
    common /comxyt/ dtcom,dxcom,dycom,dzcom, tcom,icom,jcom, kcom

    ierror = 0

    !! This should be set to actual time, in case the user wants it
    !! it for some reason in the Riemann solver.
    tcom = t
    dxcom = dx
    dycom = dy
    dzcom = dz

    !!  # Set up work arrays

    !! # From stepgrid.f
    !! mwork = msize*(46*meqn + (meqn+1)*mwaves + 9*maux + 3)

    i0faddm = 1
    i0faddp = i0faddm   + (maxm+2*mbc)*meqn
    i0gadd  = i0faddp   + (maxm+2*mbc)*meqn
    i0hadd  = i0gadd  + 6*(maxm+2*mbc)*meqn
    i0q1d   = i0hadd  + 6*(maxm+2*mbc)*meqn
    i0dtdx1 = i0q1d     + (maxm+2*mbc)*meqn
    i0dtdy1 = i0dtdx1   + (maxm+2*mbc)
    i0dtdz1 = i0dtdy1   + (maxm+2*mbc)
    i0aux1 = i0dtdz1    + (maxm+2*mbc)
    i0aux2 = i0aux1     + (maxm+2*mbc)*maux*3
    i0aux3 = i0aux2     + (maxm+2*mbc)*maux*3

    i0next = i0aux3     + (maxm+2*mbc)*maux*3    !# next free space
    mused  = i0next - 1                    !# space already used
    mwork1 = mwork - mused           !# remaining space (passed to step2)

    if (mused .gt. mwork) then
       ierror = 1
       return
    endif


    !! # Include four additional arguments to avoid need for
    !! # global array

    call clawpack46_step3(maxm,meqn,maux, mbc, &
         mx,my, mz, qold,aux,dx,dy,dz,dt, &
         cfl,fm,fp,gm,gp, hm, hp, &
         work(i0faddm),work(i0faddp), &
         work(i0gadd),work(i0hadd), &
         work(i0q1d),work(i0dtdx1),work(i0dtdy1), work(i0dtdz1), &
         work(i0aux1),work(i0aux2),work(i0aux3), &
         work(i0next),mwork1,rpn3,rpt3, rptt3, &
         mwaves,mcapa,method,mthlim,use_fwaves, ierror)

!!  # update q
    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz
    do m = 1,meqn
        do i = 1,mx
            do j = 1,my
                do k = 1,mz
                    if (mcapa .eq. 0) then
                        !! # no capa array.  Standard flux differencing:
                        qold(i,j,k, m) = qold(i,j,k,m) &
                            - dtdx * (fm(i+1,j,k,m) - fp(i,j,k,m))  &
                            - dtdy * (gm(i,j+1,k,m) - gp(i,j,k,m))  & 
                            - dtdz * (hm(i,j+1,k,m) - hp(i,j,k,m))
                    else
                        !! # with capa array.
                        qold(i,j,k, m) = qold(i,j,k, m) & 
                             -(dtdx*(fm(i+1,j,k,m) - fp(i,j,k,m)) &
                             + dtdy*(gm(i,j+1,k,m) - gp(i,j,k,m)) &
                             + dtdz*(hm(i,j+1,k,m) - hp(i,j,k,m)))/aux(i,j,k,mcapa)
                    endif
                enddo
            enddo
        enddo
    end do

end subroutine clawpack46_step3_wrap
