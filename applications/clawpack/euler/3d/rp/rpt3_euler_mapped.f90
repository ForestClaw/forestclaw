subroutine clawpack46_rpt3_mapped(ixyz,icoor,ilr,maxm,meqn,mwaves,maux,mbc, & 
    mx,ql_cart,qr_cart,aux1,aux2,aux3,asdq_cart,bmasdq_cart,bpasdq_cart)

    !!==================================================================
    !! # Riemann solver in the transverse direction for the
    !! # Euler equations.
    !! #
    !! # Uses Roe averages and other quantities which were
    !! # computed in rpn3eu and stored in the common block comroe.
    !! #
    !! #
    !! # On input,
    !! #
    !! #    ql,qr is the data along some one-dimensional slice, as in rpn3
    !! #         This slice is
    !! #             in the x-direction if ixyz=1,
    !! #             in the y-direction if ixyz=2, or
    !! #             in the z-direction if ixyz=3.
    !! #    asdq is an array of flux differences (A^*\Dq).
    !! #         asdq(i,:) is the flux difference propagating away from
    !! #         the interface between cells i-1 and i.
    !! #    Note that asdq represents B^*\Dq if ixyz=2 or C^*\Dq if ixyz=3.
    !! #
    !! #    ixyz indicates the direction of the original Riemann solve,
    !! #         called the x-like direction in the table below:
    !! #
    !! #                  normal dir         icoor=2            icoor=3
    !! #               x-like direction   y-like direction   z-like direction
    !! #      ixyz=1:        x                  y                  z
    !! #      ixyz=2:        y                  z                  x
    !! #      ixyz=3:        z                  x                  y
    !! #
    !! #    icoor indicates direction in which the transverse solve should
    !! #         be performed.
    !! #      icoor=2: split in the y-like direction.
    !! #      icoor=3: split in the z-like direction.
    !! #
    !! #    For example,
    !! #      ixyz=1, icoor=2 means asdq=A^*\Dq, and should be split in y into
    !! #                        bmasdq = B^-A^*\Dq,
    !! #                        bpasdq = B^+A^*\Dq.
    !! #
    !! #      ixyz=2, icoor=2 means asdq=B^*\Dq, and should be split in z into
    !! #                        bmasdq = C^-B^*\Dq,
    !! #                        bpasdq = C^+B^*\Dq.
    !! #
    !! #    The parameter imp is generally needed only if aux
    !! #    arrays are being used, in order to access the appropriate
    !! #    variable coefficients.


    use setprob_mod, only : gamma1, mcapa

    implicit none

    integer ixyz, icoor, ilr, maxm, meqn, mwaves,maux,mbc,mx
    double precision     ql_cart(meqn,1-mbc:maxm+mbc)
    double precision     qr_cart(meqn,1-mbc:maxm+mbc)
    double precision   asdq_cart(meqn,1-mbc:maxm+mbc)
    double precision bmasdq_cart(meqn,1-mbc:maxm+mbc)
    double precision bpasdq_cart(meqn,1-mbc:maxm+mbc)
    double precision        aux1(maux,1-mbc:maxm+mbc,3)
    double precision        aux2(maux,1-mbc:maxm+mbc,3)
    double precision        aux3(maux,1-mbc:maxm+mbc,3)

    double precision dtcom, dxcom, dycom, dzcom, tcom
    integer icom, jcom, kcom
    common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

    double precision wave(5,3),s_rot(3), asdq(5), uvw(3)
    double precision uvw_cart(3), rot(9)

    integer i, j, m, mws, i1, locrot, locarea
    double precision uvw2, pres, enth, area
    integer info, ii

    integer mu, mv, mw

    logical debug, debug_check_rpt

    debug = debug_check_rpt(ixyz,icoor,ilr,icom,jcom,kcom)


    IF (ixyz == 1) THEN
       mu = 2
       mv = 3
       mw = 4
    ELSE IF (ixyz == 2) THEN
       mu = 3
       mv = 4
       mw = 2
    ELSE
       mu = 4
       mv = 2
       mw = 3
    ENDIF


    !! # This just tells us where to find the particular rotation vectors
    !! # not in which array they are (aux1, aux2, aux3).  That we do below.
    call get_aux_locations_t(ixyz,icoor,mcapa,locrot,locarea)

    !! # Should we have maxm == mx ? 
    do i = 2-mbc, mx+mbc
        i1 = i + ilr - 2

        !! # compute values needed for Jacobian
        uvw_cart(1) = ql_cart(2,i)/ql_cart(1,i)
        uvw_cart(2) = ql_cart(3,i)/ql_cart(1,i)
        uvw_cart(3) = ql_cart(4,i)/ql_cart(1,i)

        uvw2 = uvw_cart(1)**2 + uvw_cart(2)**2 + uvw_cart(3)**2
        pres = gamma1*(ql_cart(5,i)  - 0.5d0*uvw2*ql_cart(1,i))
        enth = (ql_cart(5,i) + pres)/ql_cart(1,i)

        !! # -------------------------------------------------------
        !! # Compute bmasdq
        !! # -------------------------------------------------------

        !! # Using aux2 works for both y-like and z-like directions, 
        !! # since we are solving at the front/bottom faces.
        do j = 1,9
            rot(j) = aux2(locrot+j-1,i1,2)
        enddo        
        area = aux2(locarea,i1,2)

        do m = 1,meqn
            asdq(m) = asdq_cart(m,i)
        enddo
        call rotate3(rot,asdq(2))

        do j = 1,3
            uvw(j) = uvw_cart(j)
        enddo
        call rotate3(rot,uvw)
 
        call solve_riemann(uvw,enth,asdq, wave,s_rot,info)

        if (info > 0) then
            write(6,*) 'Calling from transverse solve; B^-'
            do j = 1,5
                write(6,1005) j,asdq(j)
            enddo
            write(6,*) ' '
            write(6,*) 'uvw'
            write(6,1001) (uvw(j),j=1,3)
            write(6,*) ' '
            write(6,*) 'uvw_cart'
            write(6,1001) (uvw_cart(j),j=1,3)
            write(6,*) ' '
            write(6,1001) enth
            write(6,*) ' '
            do ii = 1,9
               write(6,1001) 'rot = ', rot(ii)
            enddo
            write(6,*) ' '
            write(6,1002) 'ixyz = ', ixyz
            write(6,1002) 'icoor = ', icoor
            write(6,1002) 'locrot = ', locrot
            write(6,*) 'i = ', i
            write(6,1003) icom, jcom, kcom
            write(6,*) ' '
            stop
        endif
        do mws = 1,mwaves
            call rotate3_tr(rot,wave(2,mws))
        enddo

        do m=1,meqn
            bmasdq_cart(m,i) = 0.d0
            do mws = 1,mwaves
                bmasdq_cart(m,i) = bmasdq_cart(m,i) & 
                    + area*min(s_rot(mws), 0.d0) * wave(m,mws)
            enddo
        enddo

        block
            integer m
            if (debug) then
                write(6,108) 'Minus (mapped rpt) : ', ixyz,icoor,ilr, i 
                write(6,109) (area*s_rot(m),m=1,3)
                write(6,109) (bmasdq_cart(m,i),m=1,meqn)
            endif 
        end block


        !! # -------------------------------------------------------
        !! # Compute bpasdq
        !! # -------------------------------------------------------

        area = aux3(locarea,i1,2)
        if (icoor .eq. 2) then
            !! # y-like direction (back of the cell, so we need aux3)
            do j = 1,9
                rot(j) = aux3(locrot+j-1,i1,2)
            enddo
            area = aux3(locarea,i1,2)
        else
            !! # z-like direction  (top of the cell)
            do j = 1,9
                rot(j) = aux2(locrot+j-1,i1,3)
            enddo
            area = aux2(locarea,i1,3)
        endif

        do m = 1,meqn
            asdq(m) = asdq_cart(m,i)
        enddo
        call rotate3(rot,asdq(2))


        do j = 1,3
            uvw(j) = uvw_cart(j)
        enddo
        call rotate3(rot,uvw)

        call solve_riemann(uvw,enth,asdq, wave,s_rot,info)

        if (info > 0) then
            write(6,*) 'Calling from transverse solve; B^+'
            do j = 1,5
                write(6,1005) j, asdq(j)
            enddo
            write(6,*) ' '
            write(6,*) 'uvw'
            write(6,1001) (uvw(j),j=1,3)
            write(6,*) ' '
            write(6,1001) 'uvw_cart', (uvw_cart(j),j=1,3)
            write(6,*) ' '
            write(6,1001) enth
            write(6,*) ' '
            do ii = 1,9
               write(6,1001) 'rot = ', rot(ii)
            enddo
            write(6,*) ' '
            write(6,1002) 'ixyz = ', ixyz
            write(6,1002) 'icoor = ', icoor
            write(6,1002) 'locrot = ', locrot
            write(6,*) 'i = ', i
            write(6,1003) icom, jcom, kcom
            write(6,*) ' '
            stop
        endif


        do mws = 1,mwaves
            !! rotate waves back to Cartesian
            call rotate3_tr(rot,wave(2,mws))
        enddo

        do m = 1,meqn
            bpasdq_cart(m,i) = 0.d0
            do mws=1,mwaves
               bpasdq_cart(m,i) = bpasdq_cart(m,i) & 
               + area*max(s_rot(mws),0.d0) * wave(m,mws)
            enddo
        enddo

!!        block
!!            integer m
!!            double precision area
!!            if (debug) then
!!                write(6,108) 'Plus : ', ixyz,icoor,i 
!!                write(6,109) (s_rot(m)/area,m=1,3)
!!                do mws=1,mwaves
!!                    write(6,109) (bpasdq_cart(m,mws),m=1,meqn)
!!                end do
!!                write(6,*) ' '
!!            endif 
!!        end block


    enddo  !! end of i loop

108     format(A,'ixyz=',I2,'; icoor=',I2,'; ilr = ',I2,';  i=',I2)
109     format(5E24.16)

1001    format(A,3E24.16)
1002    format(A,I5)
1003    format('icom = ',I5,'; jcom = ',I5,'; kcom = ',I5)
1005    format('asdq(',I2,') = ',F24.16)




    return
end subroutine clawpack46_rpt3_mapped
