subroutine clawpack46_step3(maxm,meqn,maux,mbc,mx,my,mz, & 
    qold,aux,dx,dy,dz,dt,cflgrid, & 
    fm,fp,gm,gp,hm,hp, & 
    faddm,faddp,gadd,hadd, & 
    q1d,dtdx1d,dtdy1d,dtdz1d, & 
    aux1,aux2,aux3,work,mwork,rpn3,rpt3,rptt3, &
    mwaves,mcapa,method,mthlim,use_fwaves,block_corner_count, ierror)

    !!  ==================================================================

    !!  # clawpack routine ...  modified for AMRCLAW

    !!
    !!  # Take one time step, updating q.
    !!  # On entry, qold gives
    !!  #    initial data for this step
    !!  #    and is unchanged in this version.
    !!
    !!  # fm, fp are fluxes to left and right of single cell edge
    !!
    !!  # See the flux3 documentation for more information.
    !!
    !!      
    implicit none
    external rpn3,rpt3,rptt3

    integer :: maxm, meqn, maux, mbc, use_fwaves
    integer :: mx,my,mz, mwaves, mcapa, ierror
    integer :: mthlim(mwaves), method(7)
    double precision :: dx,dy,dz,dt,cflgrid

    !! Claw 4 indexing
    double precision :: qold(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)
    double precision :: aux(1-mbc:mx+mbc, 1-mbc:my+mbc,1-mbc:mz+mbc, maux)

    !! Claw 5 indexing
    double precision :: q1d(meqn,1-mbc:maxm+mbc)
    double precision ::  fm(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    double precision ::  fp(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    double precision ::  gm(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    double precision ::  gp(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    double precision ::  hm(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)
    double precision ::  hp(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)

    double precision :: faddm(meqn,1-mbc:maxm+mbc)
    double precision :: faddp(meqn,1-mbc:maxm+mbc)
    double precision ::  gadd(meqn,1-mbc:maxm+mbc, 2, -1:1)
    double precision ::  hadd(meqn,1-mbc:maxm+mbc, 2, -1:1)

    double precision :: aux1(maux,1-mbc:maxm+mbc, 3)
    double precision :: aux2(maux,1-mbc:maxm+mbc, 3)
    double precision :: aux3(maux,1-mbc:maxm+mbc, 3)
    double precision :: dtdx1d(1-mbc:maxm+mbc)
    double precision :: dtdy1d(1-mbc:maxm+mbc)
    double precision :: dtdz1d(1-mbc:maxm+mbc)
    double precision :: work(mwork)

    double precision :: dtcom, dxcom, dycom, dzcom, tcom
    integer :: icom,jcom,kcom
    common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

    integer :: i0wave, i0s, i0amdq, i0apdq, i0cqxx, i0bmamdq, i0bmapdq
    integer :: i0bpamdq, i0bpapdq, i0cmamdq, i0cmapdq, i0cpamdq
    integer :: i0cpapdq, i0cmamdq2, i0cmapdq2, i0cpamdq2, i0cpapdq2
    integer :: i0bmcqxxm, i0bmcqxxp, i0bpcqxxm, i0bpcqxxp, i0cmcqxxm
    integer :: i0cmcqxxp, i0cpcqxxm, i0cpcqxxp

    integer :: i0bmcmamdq, i0bmcmapdq, i0bpcmamdq, i0bpcmapdq
    integer :: i0bmcpamdq, i0bmcpapdq, i0bpcpamdq, i0bpcpapdq, iused
    integer :: mwork

    double precision :: dtdx, dtdy, dtdz, cfl1d
    integer :: i,j,k,m, ma

    integer :: block_corner_count(0:3), sweep_dir

    if (maxm .lt. max(mx,my,mz)) then
        write(6,*) 'step3 : maxm < max(mx,my,mz)'
        write(6,*) 'mx,my,mz : ', mx,my,mz
        write(6,*) 'maxm : ', maxm
        stop
    endif

    !!  # store mesh parameters that may be needed in Riemann solver but not
    !!  # passed in...
    dxcom = dx
    dycom = dy
    dtcom = dt

    !!  # partition work array into pieces needed for local storage in
    !!  # flux3 routine.  Find starting index of each piece:

    i0wave     = 1
    i0s        = i0wave     + (maxm+2*mbc)*meqn*mwaves
    i0amdq     = i0s        + (maxm+2*mbc)*mwaves
    i0apdq     = i0amdq     + (maxm+2*mbc)*meqn
    i0cqxx     = i0apdq     + (maxm+2*mbc)*meqn
    i0bmamdq   = i0cqxx     + (maxm+2*mbc)*meqn
    i0bmapdq   = i0bmamdq   + (maxm+2*mbc)*meqn
    i0bpamdq   = i0bmapdq   + (maxm+2*mbc)*meqn
    i0bpapdq   = i0bpamdq   + (maxm+2*mbc)*meqn
    i0cmamdq   = i0bpapdq   + (maxm+2*mbc)*meqn
    i0cmapdq   = i0cmamdq   + (maxm+2*mbc)*meqn
    i0cpamdq   = i0cmapdq   + (maxm+2*mbc)*meqn
    i0cpapdq   = i0cpamdq   + (maxm+2*mbc)*meqn
    i0cmamdq2  = i0cpapdq   + (maxm+2*mbc)*meqn
    i0cmapdq2  = i0cmamdq2  + (maxm+2*mbc)*meqn
    i0cpamdq2  = i0cmapdq2  + (maxm+2*mbc)*meqn
    i0cpapdq2  = i0cpamdq2  + (maxm+2*mbc)*meqn

    i0bmcqxxm   = i0cpapdq2  + (maxm+2*mbc)*meqn
    i0bmcqxxp   = i0bmcqxxm  + (maxm+2*mbc)*meqn
    i0bpcqxxm   = i0bmcqxxp  + (maxm+2*mbc)*meqn
    i0bpcqxxp   = i0bpcqxxm  + (maxm+2*mbc)*meqn
    i0cmcqxxm   = i0bpcqxxp  + (maxm+2*mbc)*meqn
    i0cmcqxxp   = i0cmcqxxm  + (maxm+2*mbc)*meqn
    i0cpcqxxm   = i0cmcqxxp  + (maxm+2*mbc)*meqn
    i0cpcqxxp   = i0cpcqxxm  + (maxm+2*mbc)*meqn

    work(i0bmcqxxm) = i0bmcqxxm
    work(i0bmcqxxp) = i0bmcqxxp
    work(i0bpcqxxm) = i0bpcqxxm
    work(i0bpcqxxp) = i0bpcqxxp
    work(i0cmcqxxm) = i0cmcqxxm
    work(i0cmcqxxp) = i0cmcqxxp
    work(i0cpcqxxm) = i0cpcqxxm
    work(i0cpcqxxp) = i0cpcqxxp

    i0bmcmamdq = i0cpcqxxp  + (maxm+2*mbc)*meqn
    i0bmcmapdq = i0bmcmamdq + (maxm+2*mbc)*meqn
    i0bpcmamdq = i0bmcmapdq + (maxm+2*mbc)*meqn
    i0bpcmapdq = i0bpcmamdq + (maxm+2*mbc)*meqn
    i0bmcpamdq = i0bpcmapdq + (maxm+2*mbc)*meqn
    i0bmcpapdq = i0bmcpamdq + (maxm+2*mbc)*meqn
    i0bpcpamdq = i0bmcpapdq + (maxm+2*mbc)*meqn
    i0bpcpapdq = i0bpcpamdq + (maxm+2*mbc)*meqn
    iused      = i0bpcpapdq + (maxm+2*mbc)*meqn - 1

    if (iused .gt. mwork) then
        !! # This shouldn't happen if mwork is set properly in stepgrid3
        write(6,*) '*** not enough work space in step3'
        write(6,*) '*** check parameter mwork in fc3d_clawpack46_step3'
        write(6,*) '*** iused = ', iused, '   mwork =',mwork
        stop
    endif

    cflgrid = 0.d0
    dtdx = dt/dx
    dtdy = dt/dy
    dtdz = dt/dz

    do m = 1,meqn
        do i = 1-mbc,mx+mbc
          do j = 1-mbc,my+mbc
            do k = 1-mbc,mz+mbc
                fm(m,i,j,k) = 0.d0
                fp(m,i,j,k) = 0.d0
                gm(m,i,j,k) = 0.d0
                gp(m,i,j,k) = 0.d0
                hm(m,i,j,k) = 0.d0
                hp(m,i,j,k) = 0.d0
             end do
          end do
       end do
    end do

    if (mcapa.eq.0) then
        !! # no capa array:
        do i = 1-mbc,maxm+mbc
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
            dtdz1d(i) = dtdz
        end do
    end if

    !! # Cubed sphere : Set corners for an x-sweep
    !! # This does nothing for non-cubed-sphere grids. 
    sweep_dir = 0
    call fc3d_clawpack46_fix_corners(mx,my,mz, mbc,meqn,qold,sweep_dir, & 
          block_corner_count)


    !!  ==================
    !!  # perform x-sweeps
    !!  ==================

    kx_loop : do k = 0,mz+1
        jx_loop : do j = 0,my+1

            !! # Store the value of j and k along this slice in the common block
            !! # comxyt in case it is needed in the Riemann solver (for
            !!  # variable coefficient problems)

            jcom = j
            kcom = k

            do m = 1,meqn
                do i = 1-mbc, mx+mbc
                    !! # copy data along a slice into 1d array:
                    q1d(m,i) = qold(i,j,k,m)
                end do
            end do

            if (mcapa .gt. 0)  then
                do i = 1-mbc, mx+mbc
                    dtdx1d(i) = dtdx / aux(i,j,k,mcapa)
                    dtdy1d(i) = dtdy / aux(i,j,k,mcapa)
                    dtdz1d(i) = dtdz / aux(i,j,k,mcapa)
                end do
            endif

            if (maux .gt. 0)  then
                do ma = 1,maux
                    do i = 1-mbc, mx+mbc
                        aux1(ma,i,1) = aux(i,j-1,k-1,ma)
                        aux1(ma,i,2) = aux(i,j-1,k,ma)
                        aux1(ma,i,3) = aux(i,j-1,k+1,ma)
                        aux2(ma,i,1) = aux(i,j,k-1,ma)
                        aux2(ma,i,2) = aux(i,j,k,ma)
                        aux2(ma,i,3) = aux(i,j,k+1,ma)
                        aux3(ma,i,1) = aux(i,j+1,k-1,ma)
                        aux3(ma,i,2) = aux(i,j+1,k,ma)
                        aux3(ma,i,3) = aux(i,j+1,k+1,ma)
                    end do
                end do
            endif

            !! # compute modifications fadd, gadd and hadd to fluxes along
            !! # this slice:

            call flux3(1,maxm,meqn,maux,mbc,mx,  & 
                    q1d,dtdx1d,dtdy1d,dtdz1d,aux1,aux2,aux3,  & 
                    faddm,faddp,gadd,hadd,cfl1d,  & 
                    work(i0wave),work(i0s),work(i0amdq),  & 
                    work(i0apdq),work(i0cqxx),  & 
                    work(i0bmamdq),work(i0bmapdq),   & 
                    work(i0bpamdq),work(i0bpapdq),   & 
                    work(i0cmamdq),work(i0cmapdq),   & 
                    work(i0cpamdq),work(i0cpapdq),   & 
                    work(i0cmamdq2),work(i0cmapdq2), & 
                    work(i0cpamdq2),work(i0cpapdq2), & 

                    work(i0bmcqxxm),work(i0bmcqxxp), & 
                    work(i0bpcqxxm),work(i0bpcqxxp), &
                    work(i0cmcqxxm),work(i0cmcqxxp), &
                    work(i0cpcqxxm),work(i0cpcqxxp), &

                    work(i0bmcmamdq),work(i0bmcmapdq), &
                    work(i0bpcmamdq),work(i0bpcmapdq), &
                    work(i0bmcpamdq),work(i0bmcpapdq), &
                    work(i0bpcpamdq),work(i0bpcpapdq), &
                    rpn3,rpt3,rptt3,mwaves,mcapa,method,mthlim,use_fwaves)

            cflgrid = dmax1(cflgrid,cfl1d)

            !! # update fluxes for use in AMR:

            do m = 1,meqn
                do i = 2-mbc,mx+mbc
                    fm(m,i,j,k) = fm(m,i,j,k) + faddm(m,i)
                    fp(m,i,j,k) = fp(m,i,j,k) + faddp(m,i)

                    gm(m,i,j  ,k-1) = gm(m,i,j  ,k-1) + gadd(m,i,1,-1)
                    gp(m,i,j  ,k-1) = gp(m,i,j  ,k-1) + gadd(m,i,1,-1)
                    gm(m,i,j  ,k  ) = gm(m,i,j  ,k  ) + gadd(m,i,1, 0)
                    gp(m,i,j  ,k  ) = gp(m,i,j  ,k  ) + gadd(m,i,1, 0)
                    gm(m,i,j  ,k+1) = gm(m,i,j  ,k+1) + gadd(m,i,1, 1)
                    gp(m,i,j  ,k+1) = gp(m,i,j  ,k+1) + gadd(m,i,1, 1)

                    gm(m,i,j+1,k-1) = gm(m,i,j+1,k-1) + gadd(m,i,2,-1)
                    gp(m,i,j+1,k-1) = gp(m,i,j+1,k-1) + gadd(m,i,2,-1)
                    gm(m,i,j+1,k  ) = gm(m,i,j+1,k  ) + gadd(m,i,2, 0)
                    gp(m,i,j+1,k  ) = gp(m,i,j+1,k  ) + gadd(m,i,2, 0)
                    gm(m,i,j+1,k+1) = gm(m,i,j+1,k+1) + gadd(m,i,2, 1)
                    gp(m,i,j+1,k+1) = gp(m,i,j+1,k+1) + gadd(m,i,2, 1)

                    hm(m,i,j-1,k) = hm(m,i,j-1,k) + hadd(m,i,1,-1)
                    hp(m,i,j-1,k) = hp(m,i,j-1,k) + hadd(m,i,1,-1)
                    hm(m,i,j  ,k) = hm(m,i,j,  k) + hadd(m,i,1, 0)
                    hp(m,i,j  ,k) = hp(m,i,j,  k) + hadd(m,i,1, 0)
                    hm(m,i,j+1,k) = hm(m,i,j+1,k) + hadd(m,i,1, 1)
                    hp(m,i,j+1,k) = hp(m,i,j+1,k) + hadd(m,i,1, 1)

                    hm(m,i,j-1,k+1) = hm(m,i,j-1,k+1) + hadd(m,i,2,-1)
                    hp(m,i,j-1,k+1) = hp(m,i,j-1,k+1) + hadd(m,i,2,-1)
                    hm(m,i,j  ,k+1) = hm(m,i,j,  k+1) + hadd(m,i,2, 0)
                    hp(m,i,j  ,k+1) = hp(m,i,j,  k+1) + hadd(m,i,2, 0)
                    hm(m,i,j+1,k+1) = hm(m,i,j+1,k+1) + hadd(m,i,2, 1)
                    hp(m,i,j+1,k+1) = hp(m,i,j+1,k+1) + hadd(m,i,2, 1)
                end do
            end do
        end do jx_loop
    end do kx_loop


    !!  ==================
    !!  # perform y sweeps
    !!  ==================


    !! # Cubed sphere : Set corners for an y-sweep
    !! # This does nothing for non-cubed-sphere grids. 
    sweep_dir = 1
    call fc3d_clawpack46_fix_corners(mx,my,mz, mbc,meqn,qold,sweep_dir, & 
                                block_corner_count)

    ky_loop : do k = 0, mz+1
        iy_loop : do i = 0, mx+1

            !! # Store the value of i and k along this slice in the common block
            !! # comxyzt in case it is needed in the Riemann solver (for
            !! # variable coefficient problems)

            icom = i
            kcom = k

            do m = 1,meqn
                do j = 1-mbc, my+mbc
                    !! # copy data along a slice into 1d array:
                    q1d(m,j) = qold(i,j,k,m)
                end do
            end do

            if (mcapa .gt. 0)  then
                do j = 1-mbc, my+mbc
                    dtdx1d(j) = dtdx / aux(i,j,k,mcapa)
                    dtdy1d(j) = dtdy / aux(i,j,k,mcapa)
                    dtdz1d(j) = dtdz / aux(i,j,k,mcapa)
                end do
            endif

            if (maux .gt. 0)  then
                do ma = 1,maux
                    do j = 1-mbc, my+mbc
                        aux1(ma,j,1) = aux(i-1,j,k-1,ma)
                        aux1(ma,j,2) = aux(i,  j,k-1,ma)
                        aux1(ma,j,3) = aux(i+1,j,k-1,ma)
                        aux2(ma,j,1) = aux(i-1,j,k,  ma)
                        aux2(ma,j,2) = aux(i,  j,k,  ma)
                        aux2(ma,j,3) = aux(i+1,j,k,  ma)
                        aux3(ma,j,1) = aux(i-1,j,k+1,ma)
                        aux3(ma,j,2) = aux(i,  j,k+1,ma)
                        aux3(ma,j,3) = aux(i+1,j,k+1,ma)
                    end do
                end do
            endif

            !! # compute modifications fadd, gadd and hadd to fluxes along this
            !! # slice:

            call flux3(2,maxm,meqn,maux,mbc,my, & 
                    q1d,dtdy1d,dtdz1d,dtdx1d,aux1,aux2,aux3, & 
                    faddm,faddp,gadd,hadd,cfl1d, & 
                    work(i0wave),work(i0s),work(i0amdq), & 
                    work(i0apdq),work(i0cqxx), & 
                    work(i0bmamdq),work(i0bmapdq), & 
                    work(i0bpamdq),work(i0bpapdq), & 
                    work(i0cmamdq),work(i0cmapdq), & 
                    work(i0cpamdq),work(i0cpapdq), & 
                    work(i0cmamdq2),work(i0cmapdq2), & 
                    work(i0cpamdq2),work(i0cpapdq2), & 

                    work(i0bmcqxxm),work(i0bmcqxxp),  &
                    work(i0bpcqxxm),work(i0bpcqxxp),  &
                    work(i0cmcqxxm),work(i0cmcqxxp),  &
                    work(i0cpcqxxm),work(i0cpcqxxp),  &

                    work(i0bmcmamdq),work(i0bmcmapdq), &
                    work(i0bpcmamdq),work(i0bpcmapdq), &
                    work(i0bmcpamdq),work(i0bmcpapdq), &
                    work(i0bpcpamdq),work(i0bpcpapdq), &
                    rpn3,rpt3,rptt3,mwaves,mcapa,method,mthlim, & 
                    use_fwaves)

            cflgrid = dmax1(cflgrid,cfl1d)

            !! # update fluxes for use in AMR:

            !! # Note that the roles of the flux updates are changed.
            !! # fadd - modifies the g-fluxes
            !! # gadd - modifies the h-fluxes
            !! # hadd - modifies the f-fluxes

            do  m = 1,meqn
                do  j = 1,my+1
                    gm(m,i,j,k) = gm(m,i,j,k) + faddm(m,j)
                    gp(m,i,j,k) = gp(m,i,j,k) + faddp(m,j)

                    hm(m,i-1,j,k) = hm(m,i-1,j,k) + gadd(m,j,1,-1)
                    hp(m,i-1,j,k) = hp(m,i-1,j,k) + gadd(m,j,1,-1)
                    hm(m,i  ,j,k) = hm(m,i  ,j,k) + gadd(m,j,1, 0)
                    hp(m,i  ,j,k) = hp(m,i  ,j,k) + gadd(m,j,1, 0)
                    hm(m,i+1,j,k) = hm(m,i+1,j,k) + gadd(m,j,1, 1)
                    hp(m,i+1,j,k) = hp(m,i+1,j,k) + gadd(m,j,1, 1)

                    hm(m,i-1,j,k+1) = hm(m,i-1,j,k+1) + gadd(m,j,2,-1)
                    hp(m,i-1,j,k+1) = hp(m,i-1,j,k+1) + gadd(m,j,2,-1)
                    hm(m,i  ,j,k+1) = hm(m,i  ,j,k+1) + gadd(m,j,2, 0)
                    hp(m,i  ,j,k+1) = hp(m,i  ,j,k+1) + gadd(m,j,2, 0)
                    hm(m,i+1,j,k+1) = hm(m,i+1,j,k+1) + gadd(m,j,2, 1)
                    hp(m,i+1,j,k+1) = hp(m,i+1,j,k+1) + gadd(m,j,2, 1)

                    fm(m,i  ,j,k-1) = fm(m,i  ,j,k-1) + hadd(m,j,1,-1)
                    fp(m,i  ,j,k-1) = fp(m,i  ,j,k-1) + hadd(m,j,1,-1)
                    fm(m,i  ,j,k  ) = fm(m,i  ,j,k  ) + hadd(m,j,1, 0)
                    fp(m,i  ,j,k  ) = fp(m,i  ,j,k  ) + hadd(m,j,1, 0)
                    fm(m,i  ,j,k+1) = fm(m,i  ,j,k+1) + hadd(m,j,1, 1)
                    fp(m,i  ,j,k+1) = fp(m,i  ,j,k+1) + hadd(m,j,1, 1)

                    fm(m,i+1,j,k-1) = fm(m,i+1,j,k-1) + hadd(m,j,2,-1)
                    fp(m,i+1,j,k-1) = fp(m,i+1,j,k-1) + hadd(m,j,2,-1)
                    fm(m,i+1,j,k  ) = fm(m,i+1,j,k  ) + hadd(m,j,2, 0)
                    fp(m,i+1,j,k  ) = fp(m,i+1,j,k  ) + hadd(m,j,2, 0)
                    fm(m,i+1,j,k+1) = fm(m,i+1,j,k+1) + hadd(m,j,2, 1)
                    fp(m,i+1,j,k+1) = fp(m,i+1,j,k+1) + hadd(m,j,2, 1)
                end do
            end do
        end do iy_loop
    end do ky_loop


    !! ==================
    !! # perform z sweeps
    !! ==================

    !! # Cubed sphere : Set corners for an y-sweep
    !! # This does nothing for non-cubed-sphere grids. 
!!    sweep_dir = 2
!!    call fc3d_clawpack46_fix_corners(mx,my,mz,mbc,meqn,qold,sweep_dir, & 
!!                                block_corner_count)


    jz_loop : do j = 0, my+1
        iz_loop : do i = 0, mx+1

            !! # Store the value of i and j along this slice in the common block
            !! # comxyzt in case it is needed in the Riemann solver (for
            !! # variable coefficient problems)

            icom = i
            jcom = j

            do m = 1,meqn
                do k = 1-mbc, mz+mbc
                    !! # copy data along a slice into 1d array:
                    q1d(m,k) = qold(i,j,k,m)
                end do
            end do

            if (mcapa .gt. 0)  then
                do k = 1-mbc, mz+mbc
                    dtdx1d(k) = dtdx / aux(i,j,k,mcapa)
                    dtdy1d(k) = dtdy / aux(i,j,k,mcapa)
                    dtdz1d(k) = dtdz / aux(i,j,k,mcapa)
                end do
            endif

            if (maux .gt. 0)  then
                do  ma = 1,maux
                    do  k = 1-mbc, mz+mbc
                        aux1(ma,k,1) = aux(i-1,j-1,k,ma)
                        aux1(ma,k,2) = aux(i-1,j,  k,ma)
                        aux1(ma,k,3) = aux(i-1,j+1,k,ma)
                        aux2(ma,k,1) = aux(i,j-1,  k,ma)
                        aux2(ma,k,2) = aux(i,j,    k,ma)
                        aux2(ma,k,3) = aux(i,j+1,  k,ma)
                        aux3(ma,k,1) = aux(i+1,j-1,k,ma)
                        aux3(ma,k,2) = aux(i+1,j,  k,ma)
                        aux3(ma,k,3) = aux(i+1,j+1,k,ma)
                    end do
                end do
            endif

            !! # compute modifications fadd, gadd and hadd to fluxes along this
            !! # slice:

            call flux3(3,maxm,meqn,maux,mbc,mz,  & 
                    q1d,dtdz1d,dtdx1d,dtdy1d,aux1,aux2,aux3, & 
                    faddm,faddp,gadd,hadd,cfl1d, & 
                    work(i0wave),work(i0s),work(i0amdq), & 
                    work(i0apdq),work(i0cqxx), &
                    work(i0bmamdq),work(i0bmapdq), & 
                    work(i0bpamdq),work(i0bpapdq), &
                    work(i0cmamdq),work(i0cmapdq), & 
                    work(i0cpamdq),work(i0cpapdq), & 
                    work(i0cmamdq2),work(i0cmapdq2), & 
                    work(i0cpamdq2),work(i0cpapdq2), &

                    work(i0bmcqxxm),work(i0bmcqxxp), &
                    work(i0bpcqxxm),work(i0bpcqxxp), & 
                    work(i0cmcqxxm),work(i0cmcqxxp), & 
                    work(i0cpcqxxm),work(i0cpcqxxp), &

                    work(i0bmcmamdq),work(i0bmcmapdq), &
                    work(i0bpcmamdq),work(i0bpcmapdq), &
                    work(i0bmcpamdq),work(i0bmcpapdq), &
                    work(i0bpcpamdq),work(i0bpcpapdq), &
                    rpn3,rpt3,rptt3,mwaves,mcapa,method,mthlim, & 
                    use_fwaves)

            cflgrid = dmax1(cflgrid,cfl1d)


            !! # update fluxes for use in AMR:
            !! # Note that the roles of the flux updates are changed.
            !! # fadd - modifies the h-fluxes
            !! # gadd - modifies the f-fluxes
            !! # hadd - modifies the g-fluxes

             do m = 1,meqn
                do k = 1,mz+1
                    hm(m,i,j,k) = hm(m,i,j,k) + faddm(m,k)
                    hp(m,i,j,k) = hp(m,i,j,k) + faddp(m,k)

                    fm(m,i  ,j-1,k) = fm(m,i  ,j-1,k) + gadd(m,k,1,-1)
                    fp(m,i  ,j-1,k) = fp(m,i  ,j-1,k) + gadd(m,k,1,-1)
                    fm(m,i  ,j  ,k) = fm(m,i  ,j  ,k) + gadd(m,k,1, 0)
                    fp(m,i  ,j  ,k) = fp(m,i  ,j  ,k) + gadd(m,k,1, 0)
                    fm(m,i  ,j+1,k) = fm(m,i  ,j+1,k) + gadd(m,k,1, 1)
                    fp(m,i  ,j+1,k) = fp(m,i  ,j+1,k) + gadd(m,k,1, 1)

                    fm(m,i+1,j-1,k) = fm(m,i+1,j-1,k) + gadd(m,k,2,-1)
                    fp(m,i+1,j-1,k) = fp(m,i+1,j-1,k) + gadd(m,k,2,-1)
                    fm(m,i+1,j  ,k) = fm(m,i+1,j  ,k) + gadd(m,k,2, 0)
                    fp(m,i+1,j  ,k) = fp(m,i+1,j  ,k) + gadd(m,k,2, 0)
                    fm(m,i+1,j+1,k) = fm(m,i+1,j+1,k) + gadd(m,k,2, 1)
                    fp(m,i+1,j+1,k) = fp(m,i+1,j+1,k) + gadd(m,k,2, 1)

                    gm(m,i-1,j  ,k) = gm(m,i-1,j  ,k) + hadd(m,k,1,-1)
                    gp(m,i-1,j  ,k) = gp(m,i-1,j  ,k) + hadd(m,k,1,-1)
                    gm(m,i  ,j  ,k) = gm(m,i  ,j  ,k) + hadd(m,k,1, 0)
                    gp(m,i  ,j  ,k) = gp(m,i  ,j  ,k) + hadd(m,k,1, 0)
                    gm(m,i+1,j  ,k) = gm(m,i+1,j  ,k) + hadd(m,k,1, 1)
                    gp(m,i+1,j  ,k) = gp(m,i+1,j  ,k) + hadd(m,k,1, 1)

                    gm(m,i-1,j+1,k) = gm(m,i-1,j+1,k) + hadd(m,k,2,-1)
                    gp(m,i-1,j+1,k) = gp(m,i-1,j+1,k) + hadd(m,k,2,-1)
                    gm(m,i  ,j+1,k) = gm(m,i  ,j+1,k) + hadd(m,k,2, 0)
                    gp(m,i  ,j+1,k) = gp(m,i  ,j+1,k) + hadd(m,k,2, 0)
                    gm(m,i+1,j+1,k) = gm(m,i+1,j+1,k) + hadd(m,k,2, 1)
                    gp(m,i+1,j+1,k) = gp(m,i+1,j+1,k) + hadd(m,k,2, 1)                    
                end do 
            end do
        end do iz_loop
    end do jz_loop

    return
end subroutine clawpack46_step3


!! #  See 'cubed_sphere_corners.ipynb'
!! 
!! NOTE : 'set_block_corner_count' is only getting set when 
!! filling ghost cells.  This means that when this routine 
!! is called for the first time, block_corner_count is set to all
!! all zeros.
subroutine fc3d_clawpack46_fix_corners(mx,my,mz, mbc,meqn,q,sweep_dir, & 
            block_corner_count)
    implicit none

    integer :: mx,my,mz, mbc,meqn,sweep_dir
    integer :: block_corner_count(0:3)
    double precision :: q(1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc,meqn)

    integer :: k,m,idata,jdata
    double precision :: ihat(0:3),jhat(0:3)
    integer :: i1, j1, ibc, jbc, kz
    logical :: use_b

    !! # Lower left corner
    ihat(0) = 0.5
    jhat(0) = 0.5

    !! # Lower right corner
    ihat(1) = mx+0.5
    jhat(1) = 0.5

    !! # Upper left corner
    ihat(2) = 0.5
    jhat(2) = my+0.5

    !! # Upper right corner
    ihat(3) = mx+0.5
    jhat(3) = my+0.5

    do kz = 1-mbc,mz+mbc
        do k = 0,3
            if (block_corner_count(k) .ne. 3) then
                cycle
            endif
            use_b = sweep_dir .eq. 0 .and. (k .eq. 0 .or. k .eq. 3) & 
                .or.  sweep_dir .eq. 1 .and. (k .eq. 1 .or. k .eq. 2) 
            do ibc = 1,mbc
                do jbc = 1,mbc
                    !! # Average fine grid corners onto coarse grid ghost corners
                    if (k .eq. 0) then
                        i1 = 1-ibc
                        j1 = 1-jbc
                    elseif (k .eq. 1) then
                        i1 = mx+ibc
                        j1 = 1-jbc
                    elseif (k .eq. 2) then
                        i1 = 1-ibc
                        j1 = my+jbc
                    else if (k .eq. 3) then
                        i1 = mx+ibc
                        j1 = my+jbc
                    else
                        write(6,*) 'fix_corners (3d) : i1,j2 not defined'
                        stop
                    endif

                    if (use_b) then
                        !! # Transform involves B                
                        idata =  j1 + int(ihat(k) - jhat(k))
                        jdata = -i1 + int(ihat(k) + jhat(k))
                    else
                        !! # Transform involves B.transpose()             
                        idata = -j1 + int(ihat(k) + jhat(k))
                        jdata =  i1 - int(ihat(k) - jhat(k))
                    endif 
                    do m = 1,meqn
                        q(i1,j1,kz,m) = q(idata,jdata,kz,m)
                    end do   !!meqn
                end do       !! jbc
            end do           !! ibc           
        end do               !! k corner
    end do                   !! kz loop

end
