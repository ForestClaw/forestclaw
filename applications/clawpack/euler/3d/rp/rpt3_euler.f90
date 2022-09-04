!-----------------------------------------------------------------------------------
! Name: rpt3(ixyz,icoor,imp,maxm,meqn,mwaves,maux,mbc,mx,ql,qr,aux1,aux2,aux3,asdq,
!            bmasdq,bpasdq)
!
! Description: Roe solver in the transverse direction for the 3D Euler Equations
!
! Inputs: 
!         ixyz <INTEGER>   : direction to take slice x-direction if ixyz=1
!                                                    y-direction if ixyz=2.
!                                                    z-direction if ixyz=3.
!         icoor <INTEGER>  : direction in which the transverse solve should be
!                            performed.  icoor=2 (split in y-like direction)
!                                        icoor=3 (split in z-like direction)
!                            ixyz=1, icoor=2 means asdq=A^*\Dq, split in y into:
!                                            bmasdq = B^-A^*\Dq,
!                                            bpasdq = B^+A^*\Dq.
!                            ixyz=2, icoor=2 means asdq=B^*\Dq, split in z into: 
!                                            bmasdq = C^-B^*\Dq,
!                                            bpasdq = C^+B^*\Dq.
!         imp <INTEGER>    : index for aux arrays
!         maxm <INTEGER>   : max number of grid cells (less the ghost cells)
!         meqn <INTEGER>   : number of equations in the system
!         mwaves <INTEGER> : nuber of waves in the system
!         maux <INTEGER>   : number of auxilary equations
!         mbc <INTEGER>    : number of ghost cells on either end
!         mx <INTEGER>     : number of elements
!         ql <REAL>        : state vector at left edge of each cell
!                            Note the i'th Riemann problem has left state qr(:,i-1)
!         qr <REAL>        : state vector at right edge of each cell
!                            Note the i'th Riemann problem has right state ql(:,i)
!         aux1 <REAL>      : 
!         aux2 <REAL>      :
!         aux3 <REAL>      :
!         asdq <REAL>      : array of flux differences (A^*\Dq).  asdq(:,i) is the
!                            flux difference propagating away from the interface
!                            between cells i-1 and i.  Note: asdq represents
!                            B^*\Dq if ixyz=2 or C^*\Dq if ixyz=3.
!         
! Outputs: 
!         bmasdq <REAL>    : left-going flux differences
!         bmasdq <REAL>    : right-going flux differences
!
! Adapted from rpt3_euler.f90 in $CLAWHOME/riemann/src
!-----------------------------------------------------------------------------------
SUBROUTINE clawpack46_rpt3(ixyz,icoor,imp,maxm,meqn,mwaves,maux,mbc,& 
      mx,ql,qr,aux1,aux2,aux3,asdq,bmasdq,bpasdq)

    USE setprob_mod, only : gamma1
    IMPLICIT NONE

    ! Input
    INTEGER, INTENT(IN) :: ixyz,icoor,imp,maxm,meqn,mwaves,maux,mbc,mx
    REAL(kind=8), DIMENSION(meqn,1-mbc:maxm+mbc), INTENT(IN) :: ql,qr,asdq
    REAL(kind=8), DIMENSION(maux,1-mbc:maxm+mbc,3), INTENT(IN) :: aux1,aux2,aux3

    ! Output
    REAL(kind=8), DIMENSION(meqn,1-mbc:maxm+mbc), INTENT(INOUT) :: bmasdq,bpasdq
  
    ! Local Storage
    INTEGER, PARAMETER :: maxmrp = 1002
    REAL(kind=8), DIMENSION(5) :: alpha
    REAL(kind=8), DIMENSION(-1:maxmrp) :: u2v2w2,u,v,w,enth,a,g1a2,euv
    REAL(kind=8) :: asqrd !! ,pl,pr,rhsq2,rhsqrtl,rhsqrtr
    REAL(kind=8), DIMENSION(5,3) :: waveb
    REAL(kind=8), DIMENSION(3) :: sb
    INTEGER :: i,j,mu,mv,mw,mws
    
    double precision dtcom, dxcom, dycom, dzcom, tcom
    integer icom, jcom, kcom
    common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

    double precision pres

    IF (maxmrp < maxm+mbc)THEN
       WRITE(6,*) 'need to increase maxmrp in rpt3_euler.f90'
       WRITE(6,*) 'maxmrp: ',maxmrp,' maxm: ',maxm,' mbc: ',mbc
       WRITE(6,*) 'maxm+mbc=',maxm+mbc
       STOP
    ENDIF

    IF (mwaves /= 3) THEN
       WRITE(6,*) '*** Should have mwaves=3 for this Riemann solver'
       STOP
    ENDIF
    
    ! Set mu to point to  the component of the system that corresponds to 
    ! momentum in the direction of this slice, mv and mw to the orthogonal 
    ! momentums.
    !
    ! ixyz indicates the direction of the original Riemann solve,
    ! called the x-like direction in the table below:
    ! 
    !         x-like direction   y-like direction   z-like direction
    ! ixyz=1:        x                  y                  z         
    ! ixyz=2:        y                  z                  x         
    ! ixyz=3:        z                  x                  y      

    IF(ixyz == 1)THEN
       mu = 2
       mv = 3
       mw = 4
    ELSE IF(ixyz == 2)THEN
       mu = 3
       mv = 4
       mw = 2
    ELSE
       mu = 4
       mv = 2
       mw = 3
    ENDIF

    ! Note that notation for u,v, and w reflects assumption that the
    !   Riemann problems are in the x-direction with u in the normal
    !   direction and v and w in the orthogonal directions, but with the
    !   above definitions of mu, mv, and mw the routine also works with
    !   ixyz=2 and ixyz = 3
    !   and returns, for example, f0 as the Godunov flux g0 for the
    !   Riemann problems u_t + g(u)_y = 0 in the y-direction.

    ! Compute the Roe-averaged variables needed in the Roe solver.

    ! Loop over grid cell interfaces
    DO i = 2-mbc, mx+mbc
       IF (qr(1,i-1) <= 0.d0 .OR. ql(1,i) <= 0.d0) THEN
          WRITE(*,*) i, mu, mv, mw
          WRITE(*,'(5e12.4)') (qr(j,i-1),j=1,5)
          WRITE(*,'(5e12.4)') (ql(j,i),j=1,5)
          IF (ixyz == 1) WRITE(6,*) '*** rho <= 0 in x-sweep at ',i
          IF (ixyz == 2) WRITE(6,*) '*** rho <= 0 in y-sweep at ',i
          IF (ixyz == 3) WRITE(6,*) '*** rho <= 0 in z-sweep at ',i
          WRITE(6,*) 'stopped with rho <= 0...'
          STOP
       ENDIF
!!       # Do we really need Roe averaged values here?        
!!       rhsqrtl = SQRT(qr(1,i-1))
!!       rhsqrtr = SQRT(ql(1,i))
!!       pl = gamma1*(qr(5,i-1) - 0.5d0*(qr(mu,i-1)**2 + &
!!            qr(mv,i-1)**2 + qr(mw,i-1)**2)/qr(1,i-1))
!!       pr = gamma1*(ql(5,i) - 0.5d0*(ql(mu,i)**2 + &
!!            ql(mv,i)**2 + ql(mw,i)**2)/ql(1,i))
!!       rhsq2 = rhsqrtl + rhsqrtr
!!       u(i) = (qr(mu,i-1)/rhsqrtl + ql(mu,i)/rhsqrtr) / rhsq2
!!       v(i) = (qr(mv,i-1)/rhsqrtl + ql(mv,i)/rhsqrtr) / rhsq2
!!       w(i) = (qr(mw,i-1)/rhsqrtl + ql(mw,i)/rhsqrtr) / rhsq2
!!       enth(i) = (((qr(5,i-1)+pl)/rhsqrtl &
!!            + (ql(5,i)+pr)/rhsqrtr)) / rhsq2
!!       u2v2w2(i) = u(i)**2 + v(i)**2 + w(i)**2
!!       asqrd = gamma1*(enth(i) - .5d0*u2v2w2(i))

        u(i) = ql(mu,i)/ql(1,i)
        v(i) = ql(mv,i)/ql(1,i)
        w(i) = ql(mw,i)/ql(1,i)
        u2v2w2(i) = u(i)**2 + v(i)**2 + w(i)**2
        pres = gamma1*(ql(5,i)  - 0.5d0*u2v2w2(i)*ql(1,i))
        enth(i) = (ql(5,i)+pres) / ql(1,i)
        asqrd = gamma1*(enth(i) - 0.5d0*u2v2w2(i))

       IF (i>=0 .AND. i<=mx .AND. asqrd <= 0.d0) THEN
          IF (ixyz == 1) WRITE(6,*) '*** a**2 <= 0 in x-sweep at ',i
          IF (ixyz == 2) WRITE(6,*) '*** a**2 <= 0 in y-sweep at ',i
          IF (ixyz == 3) WRITE(6,*) '*** a**2 <= 0 in z-sweep at ',i
          WRITE(6,*) 'stopped with a**2 < 0...'
          STOP
       ENDIF
       a(i) = SQRT(asqrd)
       g1a2(i) = gamma1 / asqrd
       euv(i) = enth(i) - u2v2w2(i)
    END DO

    ! Solve Riemann problem in the second coordinate direction
    IF(icoor == 2)THEN
       DO i = 2-mbc, mx+mbc
          alpha(4) = g1a2(i) &
               * (euv(i)*asdq(1,i) + u(i)*asdq(mu,i) + v(i)*asdq(mv,i) &
               + w(i)*asdq(mw,i) - asdq(5,i))
          alpha(2) = asdq(mu,i) - u(i)*asdq(1,i)
          alpha(3) = asdq(mw,i) - w(i)*asdq(1,i)
          alpha(5) = (asdq(mv,i) &
               + (a(i)-v(i))*asdq(1,i) - a(i)*alpha(4))/(2.d0*a(i))
          alpha(1) = asdq(1,i) - alpha(4) - alpha(5)
        
          waveb(1,1)  = alpha(1)
          waveb(mu,1) = alpha(1)*u(i)
          waveb(mv,1) = alpha(1)*(v(i)-a(i))
          waveb(mw,1) = alpha(1)*w(i)
          waveb(5,1)  = alpha(1)*(enth(i) - v(i)*a(i))
          sb(1) = v(i) - a(i)
          
          waveb(1,2)  = alpha(4)
          waveb(mv,2) = alpha(4)*v(i)
          waveb(mw,2) = alpha(4)*w(i) + alpha(3)
          waveb(mu,2) = alpha(4)*u(i) + alpha(2)
          waveb(5,2)  = alpha(4)*0.5d0*u2v2w2(i) + alpha(2)*u(i) + alpha(3)*w(i)
          sb(2) = v(i)
          
          waveb(1,3)  = alpha(5)
          waveb(mu,3) = alpha(5)*u(i)
          waveb(mv,3) = alpha(5)*(v(i) + a(i))
          waveb(mw,3) = alpha(5)*w(i)
          waveb(5,3)  = alpha(5)*(enth(i) + v(i)*a(i))
          sb(3) = v(i) + a(i)
          
          bmasdq(:,i) = 0.d0
          bpasdq(:,i) = 0.d0
          DO mws = 1,mwaves
             bmasdq(:,i) = bmasdq(:,i) + MIN(sb(mws), 0.d0)*waveb(:,mws)
             bpasdq(:,i) = bpasdq(:,i) + MAX(sb(mws), 0.d0)*waveb(:,mws)
          END DO
          if (ixyz .eq. 1 .and. jcom .eq. 4 .and. kcom .eq. 4) then
              !!write(6,201) i, jcom, enth(i), (bpasdq(j+1,i),j=1,3)
           endif
       END DO

    ! Solve Riemann problem in the third coordinate direction
    ELSEIF(icoor == 3)THEN
       DO i = 2-mbc, mx+mbc
        
          alpha(4) = g1a2(i) &
               * (euv(i)*asdq(1,i) + u(i)*asdq(mu,i) + v(i)*asdq(mv,i) &
               + w(i)*asdq(mw,i) - asdq(5,i))
          alpha(2) = asdq(mu,i) - u(i)*asdq(1,i)
          alpha(3) = asdq(mv,i) - v(i)*asdq(1,i)
          alpha(5) = (asdq(mw,i) &
               + (a(i)-w(i))*asdq(1,i) - a(i)*alpha(4))/(2.d0*a(i))
          alpha(1) = asdq(1,i) - alpha(4) - alpha(5)
          
          waveb(1,1)  = alpha(1)
          waveb(mu,1) = alpha(1)*u(i)
          waveb(mv,1) = alpha(1)*v(i)
          waveb(mw,1) = alpha(1)*(w(i) - a(i))
          waveb(5,1)  = alpha(1)*(enth(i) - w(i)*a(i))
          sb(1) = w(i) - a(i)
          
          waveb(1,2)  = alpha(4)
          waveb(mu,2) = alpha(2) + alpha(4)*u(i)
          waveb(mv,2) = alpha(3) + alpha(4)*v(i)
          waveb(mw,2) = alpha(4)*w(i)
          waveb(5,2)  = alpha(4)*0.5d0*u2v2w2(i) + alpha(2)*u(i) + alpha(3)*v(i)
          sb(2) = w(i)
          
          waveb(1,3)  = alpha(5)
          waveb(mu,3) = alpha(5)*u(i)
          waveb(mv,3) = alpha(5)*v(i)
          waveb(mw,3) = alpha(5)*(w(i) + a(i))
          waveb(5,3)  = alpha(5)*(enth(i) + w(i)*a(i))
          sb(3) = w(i) + a(i)
          
          bmasdq(:,i) = 0.d0
          bpasdq(:,i) = 0.d0
          DO mws = 1,mwaves
             bmasdq(:,i) = bmasdq(:,i) + MIN(sb(mws), 0.d0)*waveb(:,mws)
             bpasdq(:,i) = bpasdq(:,i) + MAX(sb(mws), 0.d0)*waveb(:,mws)
          END DO
          if (ixyz .eq. 3 .and. icom .eq. 4 .and. jcom .eq. 4) then
             !!write(6,201) i, kcom, enth(i), (bmasdq(j+1,i),j=1,3)
          endif
201    format(2I5,5E24.16)        

       END DO
    END IF

END SUBROUTINE clawpack46_rpt3
