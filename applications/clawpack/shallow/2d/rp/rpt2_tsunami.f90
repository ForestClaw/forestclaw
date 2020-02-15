!! =====================================================

SUBROUTINE clawpack46_rpt2(ixy,maxm,meqn,mwaves,maux,mbc,mx, & 
                   ql,qr,aux1,aux2,aux3,imp,asdq,bmasdq,bpasdq)
    
    IMPLICIT NONE
!!
!!     # Riemann solver in the transverse direction using an einfeldt
!!     Jacobian.

!!-----------------------last modified 1/10/05----------------------

    INTEGER ixy,maxm,meqn,maux,mwaves,mbc,mx,imp

    DOUBLE PRECISION      ql(1-mbc:maxm+mbc,meqn)
    DOUBLE PRECISION      qr(1-mbc:maxm+mbc,meqn)
    DOUBLE PRECISION    asdq(1-mbc:maxm+mbc,meqn)
    DOUBLE PRECISION  bmasdq(1-mbc:maxm+mbc,meqn)
    DOUBLE PRECISION  bpasdq(1-mbc:maxm+mbc,meqn)
    DOUBLE PRECISION    aux1(1-mbc:maxm+mbc,*)
    DOUBLE PRECISION    aux2(1-mbc:maxm+mbc,*)
    DOUBLE PRECISION    aux3(1-mbc:maxm+mbc,*)

    DOUBLE PRECISION  s(3)
    DOUBLE PRECISION  r(3,3)
    DOUBLE PRECISION  beta(3)
    DOUBLE PRECISION  abs_tol
    DOUBLE PRECISION  hl,hr,hul,hur,hvl,hvr,vl,vr,ul,ur,bl,br
    DOUBLE PRECISION  uhat,vhat,hhat,roe1,roe3,s1,s2,s3,s1l,s3r
    DOUBLE PRECISION  delf1,delf2,delf3,dxdcd,dxdcu
    DOUBLE PRECISION  dxdcm,dxdcp,topo1,topo3,eta
    double precision g

    DOUBLE PRECISION :: grav, dry_tolerance, sea_level
    COMMON /common_swe/ grav, dry_tolerance, sea_level

    double precision tol
    double precision z, sz, szm(3), szp(3)


    integer i,m,mw,mu,mv, mq

    tol = dry_tolerance
    abs_tol = tol

    g = grav

    if (ixy .eq. 1) then
        mu = 2
        mv = 3
    else
        mu = 3
        mv = 2
    endif

    do i = 2-mbc,mx+mbc

        hl  = qr(i-1,1) 
        hr  = ql(i,1) 
        hul = qr(i-1,mu) 
        hur = ql(i,mu) 
        hvl = qr(i-1,mv) 
        hvr = ql(i,mv)

!!===========determine velocity from momentum===========================

        if (hl .lt. abs_tol) then
            hl = 0.d0
            ul = 0.d0
            vl = 0.d0
        else
            ul = hul/hl
            vl = hvl/hl
        endif

        if (hr .lt. abs_tol) then
            hr = 0.d0
            ur = 0.d0
            vr = 0.d0
        else
            ur = hur/hr
            vr = hvr/hr
        endif

        do mw=1,mwaves
            s(mw) = 0.d0
            beta(mw) = 0.d0
            do m=1,meqn
                r(m,mw) = 0.d0
            enddo
        enddo
        dxdcp = 1.d0
        dxdcm = 1.d0

        if (hl <= tol .and. hr <= tol) go to 90

!!      !check and see if cell that transverse waves are going in is high and dry

        if (imp .eq. 1) then
            eta = qr(i-1,1)  + aux2(i-1,1)
            topo1 = aux1(i-1,1)
            topo3 = aux3(i-1,1)
!!            s1 = vl-sqrt(g*hl)
!!            s3 = vl+sqrt(g*hl)
!!            s2 = 0.5d0*(s1+s3)
        else
            eta = ql(i,1) + aux2(i,1)
            topo1 = aux1(i,1)
            topo3 = aux3(i,1)
!!            s1 = vr-sqrt(g*hr)
!!            s3 = vr+sqrt(g*hr)
!!            s2 = 0.5d0*(s1+s3)
        endif
        if (eta .lt. max(topo1,topo3)) go to 90

!!      if (coordinate_system.eq.2) then
!!         if (ixy.eq.2) then
!!             dxdcp=(earth_radius*deg2rad)
!!            dxdcm = dxdcp
!!         else
!!            if (imp.eq.1) then
!!               dxdcp = earth_radius*cos(aux3(3,i-1))*deg2rad
!!               dxdcm = earth_radius*cos(aux1(3,i-1))*deg2rad
!!            else
!!               dxdcp = earth_radius*cos(aux3(3,i))*deg2rad
!!               dxdcm = earth_radius*cos(aux1(3,i))*deg2rad
!!            endif
!!         endif
!!      endif

!!=====Determine some speeds necessary for the Jacobian=================

        vhat = (vr*dsqrt(hr))/(dsqrt(hr) + dsqrt(hl)) + & 
               (vl*dsqrt(hl))/(dsqrt(hr) + dsqrt(hl))

        uhat = (ur*dsqrt(hr))/(dsqrt(hr) + dsqrt(hl)) + & 
               (ul*dsqrt(hl))/(dsqrt(hr) + dsqrt(hl))

        hhat = (hr + hl)/2.d0

        roe1 = vhat - dsqrt(g*hhat)
        roe3 = vhat + dsqrt(g*hhat)

        s1l = vl - dsqrt(g*hl)
        s3r = vr + dsqrt(g*hr)

        s1 = dmin1(roe1,s1l)
        s3 = dmax1(roe3,s3r)

        !!s2 = 0.5d0*(s1+s3)
        s2 = uhat  !! This is what the shallow water eqns do

        s(1) = s1
        s(2) = s2
        s(3) = s3
!!=======================Determine asdq decomposition (beta)============
        delf1 = asdq(i,1)
        delf2 = asdq(i,mu)
        delf3 = asdq(i,mv)

        beta(1) = (s3*delf1/(s3-s1))-(delf3/(s3-s1))
        beta(2) = -s2*delf1 + delf2
        beta(3) = (delf3/(s3-s1))-(s1*delf1/(s3-s1))
!!======================End =================================================

!!=====================Set-up eigenvectors===================================
        r(1,1) = 1.d0
        r(mu,1) = s2
        r(mv,1) = s1

        r(1,2) = 0.d0
        r(mu,2) = 1.d0
        r(mv,2) = 0.d0

        r(1,3) = 1.d0
        r(mu,3) = s2
        r(mv,3) = s3
!!============================================================================
90      continue
!!============= compute fluctuations==========================================


        do mw = 1,mwaves
            z = sign(1.d0,s(mw))
            szm(mw) = (1-z)/2*s(mw)*beta(mw)
            szp(mw) = (1+z)/2*s(mw)*beta(mw)
        end do

        do mq = 1,meqn
            bmasdq(i,mq) =                szm(1)*r(mq,1)
            bmasdq(i,mq) = bmasdq(i,mq) + szm(2)*r(mq,2)
            bmasdq(i,mq) = bmasdq(i,mq) + szm(3)*r(mq,3)

            bpasdq(i,mq) =                szp(1)*r(mq,1)
            bpasdq(i,mq) = bpasdq(i,mq) + szp(2)*r(mq,2)
            bpasdq(i,mq) = bpasdq(i,mq) + szp(3)*r(mq,3)
        enddo 

!!========================================================================
    enddo

      return
      end
