subroutine clawpack46_rpn3_mapped(ixyz,maxm,meqn,mwaves,maux,mbc,mx,& 
    ql_cart,qr_cart, auxl,auxr,wave_cart,s,amdq_cart,apdq_cart)
!!
!!     # Roe-solver for the Euler equations
!!      -----------------------------------------------------------
!!
!!     # solve Riemann problems along one slice of data.
!!     # This data is along a slice in the x-direction if ixyz=1
!!     #                               the y-direction if ixyz=2.
!!     #                               the z-direction if ixyz=3.
!!
!!     # On input, ql contains the state vector at the left edge of each cell
!!     #           qr contains the state vector at the right edge of each cell
!!
!!     # On output, wave contains the waves, s the speeds,
!!     # and amdq, apdq the left-going and right-going flux differences,
!!     # respectively.
!!
!!     # Note that the i'th Riemann problem has left state qr(i-1,:)
!!     #                                    and right state ql(i,:)
!!     # From the basic clawpack routines, this routine is called with ql = qr
!!

    use setprob_mod, only : gamma, gamma1, mcapa
    implicit none

    !! Riemann solvers all use (meqn,i,j,k) indexing (5.x indexing)
    integer ixyz, maxm, meqn, mwaves, mbc, mx, maux
    double precision wave_cart(meqn,mwaves,1-mbc:maxm+mbc)
    double precision         s(mwaves,1-mbc:maxm+mbc)
    double precision   ql_cart(meqn,1-mbc:maxm+mbc)
    double precision   qr_cart(meqn,1-mbc:maxm+mbc)
    double precision amdq_cart(meqn,1-mbc:maxm+mbc)
    double precision apdq_cart(meqn,1-mbc:maxm+mbc)
    double precision      auxl(maux,1-mbc:maxm+mbc)
    double precision      auxr(maux,1-mbc:maxm+mbc)

    double precision dtcom, dxcom, dycom, dzcom, tcom
    integer icom, jcom, kcom
    common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

    double precision ql(5),qr(5), s_rot(3), wave(5,3)
    double precision amdq(5), apdq(5), rot(9), uvw(3)


    !!  local arrays -- common block comroe is passed to rpt3eu

    double precision delta(5)
    logical efix
    integer i, j, ii, m, info, mws, locrot, locarea
    double precision pr, pl, enth, rhsq2, rhsqrtl, rhsqrtr
    double precision rho_im1, pim1, cim1, s0
    double precision rho1, rhou1, rhov1, rhow1, en1, p1, c1
    double precision s1, sfract, rhoi, pi, ci, s3
    double precision rho2, rhou2, rhov2, rhow2, en2, p2, c2
    double precision s2, df, area, uvw2

    logical debug

    data efix /.false./    !# use entropy fix for transonic rarefactions

    debug = .false.
    if (ixyz .eq. 1 .and. jcom .eq. 4 .and. kcom .eq. 4) then
        debug = .true.
    endif
    debug = .false.

    if (mwaves .ne. 3) then
        write(6,*) '*** Should have mwaves=3 for this Riemann solver'
        stop
    endif

    !! Locations in aux arrays
    call get_aux_locations_n(ixyz,mcapa,locrot,locarea)

    do i = 2-mbc, mx+mbc

        do m = 1,meqn
            ql(m) = ql_cart(m,i)
            qr(m) = qr_cart(m,i-1)
        enddo

        !!# Use Roe averaged values for normal solves
        if (qr(1) .lt. 0) then
            write(6,*) 'qr(1) < 0; ', qr(1)
            stop
        endif
         if (ql(1) .lt. 0) then
            write(6,*) 'ql(1) < 0; ', ql(1)
            stop
        endif

        !! We can use unrotated values here, since we only use 
        !! norm of velocity.
        rhsqrtl = sqrt(qr(1))
        rhsqrtr = sqrt(ql(1))
        rhsq2 = rhsqrtl + rhsqrtr

        uvw2 = qr(2)**2 + qr(3)**2 + qr(4)**2
        pl = gamma1*(qr(5) - 0.5d0*uvw2/qr(1))

        uvw2 = ql(2)**2 + ql(3)**2 + ql(4)**2
        pr = gamma1*(ql(5) - 0.5d0*(uvw2)/ql(1))

        enth = (((qr(5)+pl)/rhsqrtl + (ql(5)+pr)/rhsqrtr)) / rhsq2


        !! --------------- Use rotated values in Riemann solve ------------
        do j = 1,9
            rot(j) = auxl(locrot+j-1,i)
        enddo

        do j = 1,3
            uvw(j) = (qr(j+1)/rhsqrtl + ql(j+1)/rhsqrtr) / rhsq2
        enddo
        call rotate3(rot,uvw)

        do j = 1,meqn
            delta(j) = ql(j) - qr(j)
        enddo
        call rotate3(rot,delta(2))

        !! # Solve normal Riemann problem
        call solve_riemann(uvw, enth, delta, wave,s_rot,info)

        if (info > 0) then
            write(6,'(A)') 'In rpn3'
            write(6,1001) 'rhol = ',qr_cart(1,i-1)
            write(6,1001) 'rhor = ',ql_cart(1,i)
            write(6,1001) 'ul   = ',qr_cart(2,i-1)/qr_cart(1,i-1)
            write(6,1001) 'ur   = ',ql_cart(2,i)/ql_cart(1,i)
            write(6,1001) 'pl   = ',pl
            write(6,1001) 'pr   = ',pr
            do ii = 1,5
               write(6,1001) 'ql = ', ql(ii)
            enddo
            write(6,*) ' '
            do ii = 1,5
               write(6,1001) 'qr = ', qr(ii)
            enddo
            write(6,*) ' '
            do ii = 1,9
               write(6,1001) 'rot = ', rot(ii)
            enddo
            write(6,*) 'locarea = ', locarea
            write(6,*) 'locrot = ', locrot
            write(6,*) 'ixyz = ', ixyz
            write(6,'(A,I5)') 'i = ', i
            stop
        endif
1001    format(A,2E16.8)

        if (efix) then

            !! # check 1-wave:
            !! ---------------

            rho_im1 = qr(1)
            uvw2 = qr(2)**2 + qr(3)**2 + qr(4)**2
            pim1 = gamma1*(qr(5) - 0.5d0*uvw2/rho_im1)
            cim1 = sqrt(gamma*pim1/rho_im1)

            !! # u-c in left state (cell i-1)
            s0 = qr(2)/rho_im1 - cim1 


            !! # check for fully supersonic case:
            !! right_going = .false.
            if (s0 .ge. 0.d0 .and. s_rot(1) .gt. 0.d0) then
                !! # everything is right-going
                amdq = 0.d0
                !! right_going = .true.
            else
                rho1 = qr(1) + wave(1,1)
                rhou1 = qr(2) + wave(2,1)
                rhov1 = qr(3) + wave(3,1)
                rhow1 = qr(4) + wave(4,1)
                en1 = qr(5) + wave(5,1)
                uvw2 = rhou1**2 + rhov1**2 + rhow1**2
                p1 = gamma1*(en1 - 0.5d0*uvw2/rho1)
                c1 = sqrt(gamma*p1/rho1)

                !! # u-c to right of 1-wave_rot
                s1 = rhou1/rho1 - c1     

                !!left_going = .false.
                if (s0 .lt. 0.d0 .and. s1 .gt. 0.d0) then
                    !! # transonic rarefaction in the 1-wave
                    sfract = s0 * (s1-s_rot(1)) / (s1-s0)
                else if (s_rot(1) .lt. 0.d0) then
                    !! # 1-wave is leftgoing
                    sfract = s_rot(1)
                else
                    !! # 1-wave is right-going
                    !! # this shouldn't happen since s0 < 0
                    sfract = 0.d0         
                endif
                do m = 1,meqn
                    amdq(m) = sfract*wave(m,1)
                end do
            endif

            !! # check 2-wave:
            !! ---------------

            if (s_rot(2) .ge. 0.d0) then 
                !! # 2-,3- and 4- waves are rightgoing
                go to 200 
            endif

            do m=1,meqn
                amdq(m) = amdq(m) + s_rot(2)*wave(m,2)
            enddo

           !! # check 3-wave:
           !! ---------------

            rhoi = ql(1)
            uvw2 = ql(2)**2 + ql(3)**2 + ql(4)**2
            pi = gamma1*(ql(5) - 0.5d0*(uvw2) / rhoi)
            ci = dsqrt(gamma*pi/rhoi)

            !! # u+c in right state  (cell i)
            s3 = ql(2)/rhoi + ci     

            rho2 = ql(1) - wave(1,3)
            rhou2 = ql(2) - wave(2,3)
            rhov2 = ql(3) - wave(3,3)
            rhow2 = ql(4) - wave(4,3)
            en2 = ql(5) - wave(5,3)

            uvw2 = rhou2**2 + rhov2**2 + rhow2**2
            p2 = gamma1*(en2 - 0.5d0*(uvw2)/rho2)
            c2 = dsqrt(gamma*p2/rho2)
            s2 = rhou2/rho2 + c2        !# u+c to left of 3-wave
            if (s2 .lt. 0.d0 .and. s3 .gt. 0.d0 ) then
                !! # transonic rarefaction in the 3-wave
                sfract = s2 * (s3-s_rot(3)) / (s3-s2)
            else if (s_rot(3) .lt. 0.d0) then
                !! # 3-wave is leftgoing
                sfract = s_rot(3)
            else
                !! # 3-wave is rightgoing
                go to 200
            endif

            do m = 1,5
                amdq(m) = amdq(m) + sfract*wave(m,3)
            enddo

  200       continue

           !! # compute the rightgoing flux differences:
           !! # df = SUM s*wave   is the total flux difference and apdq = df -
           !! # amdq

            do m = 1,meqn
               df = 0.d0
               do mws = 1,mwaves
                    df = df + s_rot(mws)*wave(m,mws)
               enddo
               apdq(m) = df - amdq(m)
            enddo
        else
            !! # No entropy fix
            do mws = 1,mwaves
                call rotate3_tr(rot,wave(2,mws))
                do m = 1,meqn
                    wave_cart(m,mws,i) = wave(m,mws)
                enddo 
            enddo

            area = auxl(locarea,i)
            do mws = 1,mwaves
                s(mws,i) = area*s_rot(mws)
            enddo

            do m=1,meqn
                amdq_cart(m,i) = 0.d0
                apdq_cart(m,i) = 0.d0
                do  mws = 1,mwaves
                    amdq_cart(m,i) = amdq_cart(m,i) + min(s(mws,i),0.d0)*wave(m,mws)
                    apdq_cart(m,i) = apdq_cart(m,i) + max(s(mws,i),0.d0)*wave(m,mws)
                enddo
            enddo
        endif  !! end of efix

        if (debug) then
            write(6,211) 1, i, (ql_cart(j,i),j=1,5)
            !!write(6,211) i, (qr_cart(j,i),j=1,5)
            !!write(6,211) i, (delta(j),j=1,5)
            !!write(6,211) i, (amdq_cart(j,i)/area,j=1,5)
            !!write(6,211) i, (apdq_cart(j,i)/area,j=1,5)
            !!write(6,*) ' '
        endif
211 format(2I5,5E16.8)

    enddo  !! end of i loop over 1d sweep array

    return
end subroutine clawpack46_rpn3_mapped


subroutine solve_riemann(uvw,enth,delta,wave,s,info)
    use setprob_mod, only : gamma1
    implicit none

    double precision enth, uvw(3), delta(5)
    double precision wave(5,3), s(3)
    double precision u2v2w2, a2, a, g1a2, euv
    double precision a1, a3, a4, a5,u,v,w
    integer info

    info = 0

    u = uvw(1)
    v = uvw(2)
    w = uvw(3)

    u2v2w2 = u**2 + v**2 + w**2
    a2 = gamma1*(enth - 0.5d0*u2v2w2)
    if (a2 .lt. 0.d0) then
        write(6,*) 'solve_riemann : '
        write(6,*) 'a2 .lt. 0; ', a2
        write(6,*) enth, u2v2w2
        info = 1
        return
    endif
    a = sqrt(a2)
    g1a2 = gamma1 / a2
    euv = enth - u2v2w2

    a4 = g1a2 * (euv*delta(1) + u*delta(2) + v*delta(3) + w*delta(4) & 
                 - delta(5))
    a2 = delta(3) - v*delta(1)
    a3 = delta(4) - w*delta(1)
    a5 = (delta(2) + (a-u)*delta(1) - a*a4) / (2.d0*a)
    a1 = delta(1) - a4 - a5

    !! # Compute the waves.
    !! # Note that the 2-wave, 3-wave and 4-wave travel at the same speed
    !! # and are lumped together in wave(.,.,2).  The 5-wave is then stored
    !! # in wave(.,.,3).

    wave(1,1) = a1
    wave(2,1) = a1*(u-a)
    wave(3,1) = a1*v
    wave(4,1) = a1*w
    wave(5,1) = a1*(enth - u*a)
    s(1) = u - a

    wave(1,2) = a4
    wave(2,2) = a4*u
    wave(3,2) = a4*v + a2
    wave(4,2) = a4*w + a3
    wave(5,2) = a4*0.5d0*u2v2w2  + a2*v + a3*w
    s(2) = u

    wave(1,3) = a5
    wave(2,3) = a5*(u+a)
    wave(3,3) = a5*v
    wave(4,3) = a5*w
    wave(5,3) = a5*(enth+u*a)
    s(3) = u + a

end subroutine solve_riemann
