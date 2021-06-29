!! For each coarse-fine interface, a Riemann problem between an inner
!! ghost cell value on the fine grid and cell value in the adjacent coarse
!! cell must be solved and added to corresponding location in
!! **node(ffluxptr, mptr)** for conservative fix later
!!
!! -------------------------------------------------------------
!!

SUBROUTINE qad_new(valbig,mitot,mjtot,nvar, &
   svdflx,qc1d,lenbc,lratiox,lratioy,hx,hy, &
   maux,auxbig,auxc1d,delt,mptr)

  USE amr_module, only : timemult, nghost, max1d, maxaux, mwaves,nestlevel,rnode, &
       auxtype, node, method, max1d
  IMPLICIT NONE

  INTEGER mitot, mjtot, nvar, lenbc, lratiox, lratioy, maux
  INTEGER mptr
  DOUBLE PRECISION hx,hy,delt

  DOUBLE PRECISION, TARGET :: valbig(nvar,mitot,mjtot)
  DOUBLE PRECISION, TARGET :: auxbig(maux,mitot,mjtot)
  DOUBLE PRECISION qc1d(nvar,lenbc)
  DOUBLE PRECISION svdflx(nvar,lenbc)
  DOUBLE PRECISION auxc1d(maux,lenbc)

  DOUBLE PRECISION, POINTER :: q(:,:,:)
  DOUBLE PRECISION, POINTER :: aux(:,:,:)

  !!
  !! ::::::::::::::::::::::::::: QAD ::::::::::::::::::::::::::::::::::
  !!  are added in to coarse grid value, as a conservation fixup.
  !!  Done each fine grid time step. If source terms are present, the
  !!  coarse grid value is advanced by source terms each fine time step too.

  !!  No change needed in this sub. for spherical mapping: correctly
  !!  mapped vals already in bcs on this fine grid and coarse saved
  !!  vals also properly prepared
  !!
  !! Side 1 is the left side of the fine grid patch.  Then go around clockwise.
  !! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  !!
  !!      # local storage

  !! max1d is largest possible grid size (defined in amr_module.f90)
  INTEGER :: max1dp1
  DOUBLE PRECISION wave(nvar,mwaves,max1d+1)
  DOUBLE PRECISION s(mwaves,max1d+1)
  DOUBLE PRECISION amdq(nvar,max1d+1)  
  DOUBLE PRECISION apdq(nvar,max1d+1)

  DOUBLE PRECISION, TARGET :: auxlbig(maxaux*(max1d+1))
  DOUBLE PRECISION, TARGET :: auxrbig(maxaux*(max1d+1))

  DOUBLE PRECISION, POINTER :: auxf(:,:)
  DOUBLE PRECISION, POINTER :: auxc(:,:)

  !! Temporary variables (used instead of ql,qr)
  DOUBLE PRECISION qf(nvar,max1d+1)
  DOUBLE PRECISION qc(nvar,max1d+1)

  DOUBLE PRECISION tgrid
  INTEGER nc, nr, level, index, l
  INTEGER i,j,ma,lind, ncrse, ic, jc, ifine, jfine
  INTEGER iaux

  INTEGER mx,my,mbc,meqn, mxc, myc, mq
  DOUBLE PRECISION dt, dx, dy, delta_fix
  INTEGER iface, idir
  LOGICAL prt

  max1dp1 = max1d + 1

  tgrid = rnode(timemult, mptr)
  nr = mitot-2*nghost
  nc = mjtot-2*nghost
  level = node(nestlevel, mptr)
  index = 0

  !! Rename variables to use Clawpack convention
  mbc = nghost
  mx = nr
  my = nc
  meqn = nvar
  dt = delt
  dx = hx
  dy = hy

  mxc = mx/lratiox
  myc = my/lratioy

  !! Redimension arrays to use indexing that starts at 1-mbc, etc
  q(1:meqn,1-mbc:mx+mbc,1-mbc:my+mbc) => valbig

  if (maux .gt. 0) then
     aux(1:maux,1-mbc:mx+mbc,1-mbc:my+mbc) => auxbig
     auxf(1:maux,1:max1dp1) => auxlbig
     auxc(1:maux,1:max1dp1) => auxrbig
  endif


  !! Counter for indexing into 1d arrays of coarse grid values 
  index = 0

  !! ----------------------------------------------
  !!  Side 1 
  !! ----------------------------------------------
  !!  
  !! iface = 0 (left edge)
  !! 
  !! Fine grid is on the right;  coarse grid is on the left
  !!
  !! auxf, qf : data for right (fine) cell
  !! auxc, qc : data for left (coarse) cell
  !!

  DO j = 1,my
     IF (maux .GT. 0) THEN
        DO ma = 1,maux
           IF (auxtype(ma) .EQ. "xleft") THEN
              !! # Assuming velocity at left-face, this fix
              !! # preserves conservation in incompressible flow:
              auxf(ma,j) = aux(ma,1,j)
           ELSE
              !! # Normal case -- we set the aux arrays
              !! # from the cell corresponding  to q
              auxf(ma,j) = aux(ma,0,j)
           ENDIF
        ENDDO
     ENDIF
     DO mq = 1,meqn
        qf(mq,j) = q(mq,0,j)
     ENDDO
  ENDDO

  !! Side 1
  DO jc = 1,myc
     DO l = 1,lratioy
        jfine = (jc-1)*lratioy + l
        IF (maux .GT. 0) THEN
           DO ma = 1,maux
              auxc(ma,jfine) = auxc1d(ma,index + jc)
           ENDDO
        ENDIF
        DO mq = 1,meqn
           qc(mq,jfine) = qc1d(mq,index + jc)
        ENDDO
     ENDDO
  ENDDO

  idir = 0
  iface = 0    !! Face of Cartesian grid
  CALL rpn2qad_new(my,meqn,maux,mbc, idir, iface, &
                    qf,qc,auxf,auxc,amdq,apdq)

  !! Side 1
  DO jc = 1,myc
     DO l = 1,lratioy
        jfine = (jc-1)*lratioy + l
        DO  mq = 1,meqn
           delta_fix = amdq(mq,jfine) + apdq(mq,jfine)
           svdflx(mq,index+jc) = svdflx(mq,index+jc) + dy*dt*delta_fix
        ENDDO
     ENDDO
  ENDDO
  index  = myc

  !! ----------------------------------------------
  !!  Side 2
  !! ----------------------------------------------
  !!
  !! iface = 4 (top edge)
  !!
  !! Fine grid is on the left (bottom);  coarse grid is on the right (top)
  !!
  !! auxf, qf : data for bottom (fine) cell
  !! auxc, qc : data for top (coarse) cell
  !!

  IF (my .EQ. 1) THEN
     !! # a single row of interior cells only happens when using the
     !! # 2d amrclaw code to do a 1d problem with refinement.
     !! # (feature added in Version 4.3)
     !! # skip over sides 2 and 4 in this case
     go to 299
  ENDIF


  !! Side 2
  DO i = 1,mx
     IF (maux .GT. 0) THEN
        DO ma = 1,maux
           auxf(ma,i) = aux(ma,i,my+1)
        ENDDO
     ENDIF
     DO mq = 1,meqn
        qf(mq,i) = q(mq,i,my+1)
     ENDDO
  ENDDO

  !! Side 2
  DO ic = 1, mxc
     DO l = 1,lratiox
        ifine = (ic-1)*lratiox + l
        IF (maux .GT. 0) THEN
           DO  ma = 1,maux
              IF (auxtype(ma) .EQ. "yleft") THEN
                 !! This only makes sense if coarse and fine grids align at their 
                 !! shared boundary.
                 auxc(ma,ifine) = aux(ma,ifine,my+1)
              ELSE
                 auxc(ma,ifine) = auxc1d(ma,index + ic)
              ENDIF
           ENDDO
        ENDIF
        DO  mq = 1,meqn
           qc(mq,ifine) = qc1d(mq,index + ic)
        ENDDO
     ENDDO
  ENDDO

  idir = 1
  iface = 3
  CALL rpn2qad_new(mx,meqn,maux,mbc, idir, iface, &
                    qf,qc,auxf,auxc,amdq,apdq)

  !! Side 2
  DO ic = 1,mxc
     DO l = 1,lratiox
        ifine = (ic-1)*lratiox + l
        DO mq = 1,meqn
           delta_fix = amdq(mq,ifine) + apdq(mq,ifine)
           svdflx(mq,index+ic) = svdflx(mq,index+ic) + dx*dt*delta_fix
        ENDDO
     ENDDO
  ENDDO
  index  = index + mxc

299 continue


  !! ----------------------------------------------
  !!  Side 3
  !! ----------------------------------------------
  !!
  !! iface = 1 (right edge)
  !!
  !! Fine grid is on the left;  coarse grid is on the right
  !!
  !! auxf, qf : data for left (fine) cell
  !! auxc, qc : data for right (coarse) cell
  !!

  DO j = 1, my
     IF (maux .GT. 0) THEN
        DO ma = 1,maux
           auxf(ma,j) = aux(ma,mx+1,j)
        ENDDO
     ENDIF
     DO mq = 1, meqn
        qf(mq,j) = q(mq,mx+1,j)
     ENDDO
  ENDDO

  !! Side 3
  DO jc = 1,myc
     DO l = 1,lratioy
        jfine = (jc-1)*lratioy + l
        IF (maux .GT. 0) THEN
           DO ma = 1,maux
              IF (auxtype(ma).EQ."xleft") THEN
                 !! This only makes sense if coarse and fine grids align at their 
                 !! shared boundary.
                 auxc(ma,jfine) = aux(ma,mx+1,jfine)
              ELSE
                 auxc(ma,jfine) = auxc1d(ma,index + jc)
              ENDIF
           ENDDO
        ENDIF
        DO mq = 1, meqn
           qc(mq,jfine) = qc1d(mq,index + jc)
        ENDDO
     ENDDO
  ENDDO


  idir = 0
  iface = 1
  CALL rpn2qad_new(my,meqn,maux,mbc, idir, iface, &
                    qf,qc,auxf,auxc,amdq,apdq)

  !! Side 3
  DO jc = 1, myc
     DO l = 1, lratioy
        jfine = (jc-1)*lratioy + l
        DO mq = 1, meqn
           delta_fix = amdq(mq,jfine) + apdq(mq,jfine)
           svdflx(mq,index + jc) = svdflx(mq,index + jc) + dy*dt*delta_fix
        ENDDO
     ENDDO
  ENDDO
  index  = index + myc



  !! ----------------------------------------------
  !!  Side 4
  !! ----------------------------------------------
  !! 
  !! iface = 2 (bottom edge)
  !!
  !! Fine grid is on the top (right);  coarse grid is on the bottom (left)
  !!
  !! auxf, qf : data for  right (fine) cell
  !! auxc, qc : data for left (coarse) cell
  !!

  IF (my .EQ. 1) THEN
     !! # a single row of interior cells only happens when using the
     !! # 2d amrclaw code to do a 1d problem with refinement.
     !! # (feature added in Version 4.3)
     !! # skip over sides 2 and 4 in this case
     go to 499
  ENDIF

  !! Side 4
  DO  i = 1, mx
     IF (maux .GT. 0) THEN
        !! Is this conditional needed?  Loop won't do anything if maux == 0
        DO ma = 1,maux
           IF (auxtype(ma) .EQ. "yleft") THEN
              auxf(ma,i) = aux(ma,i,1)
           ELSE
              auxf(ma,i) = aux(ma,i,0)
           ENDIF
        ENDDO
     ENDIF
     DO mq = 1, meqn
        qf(mq,i) = q(mq,i,0)
     ENDDO
  ENDDO

  !! Side 4
  DO ic = 1,mxc
     DO l = 1,lratiox
        ifine = (ic-1)*lratiox + l
        IF (maux .GT. 0) THEN
           DO ma = 1,maux
              auxc(ma,ifine) = auxc1d(ma,index + ic)
           ENDDO
        ENDIF
        DO  mq = 1, meqn
           qc(mq,ifine) = qc1d(mq,index + ic)
        ENDDO
     ENDDO
  ENDDO

  idir = 1
  iface = 2
  CALL rpn2qad_new(mx,meqn,maux,mbc, idir, iface, &                
                    qf,qc,auxf,auxc,amdq,apdq)

  !! Side 4
  DO ic = 1,mxc
     DO l = 1,lratiox
        ifine = (ic-1)*lratiox + l
        DO mq = 1,meqn
           delta_fix = amdq(mq,ifine) + apdq(mq,ifine)
           svdflx(mq,index + ic) = svdflx(mq,index + ic) + dx*dt*delta_fix
        ENDDO
     ENDDO
  ENDDO

499 continue

  !!      # for source terms:
  IF (method(5) .NE. 0) THEN   ! should I test here if index=0 and all skipped?
     !!   call src1d(meqn,nghost,lenbc,qc1d,maux,auxc1d,tgrid,delt)
     !! # how can this be right - where is the integrated src term used?
  ENDIF

  RETURN
END SUBROUTINE qad_new


SUBROUTINE rpn2qad_new(mx,meqn,maux,mbc, idir, iface, &
                       qf,qc,auxf,auxc,amdq,apdq)

  IMPLICIT NONE
  INTEGER mx, meqn, maux, mbc, idir, iface
  DOUBLE PRECISION qf(meqn,mx), qc(meqn,mx)
  DOUBLE PRECISION auxf(maux,mx), auxc(maux,mx)

  DOUBLE PRECISION amdq(meqn,mx)  
  DOUBLE PRECISION apdq(meqn,mx)

  !! automatically allocated
  DOUBLE PRECISION qvc(meqn), qvf(meqn)
  DOUBLE PRECISION auxvc(maux), auxvf(maux)
  DOUBLE PRECISION fluxf(meqn), fluxc(meqn)

  DOUBLE PRECISION fd
  INTEGER m,i, iface_cell, sgn

  !! idir refers to direction of the Riemann solver
  !! iface refers to the face of the Cartesian grid (not face of cell)

  DO i = 1,mx
    do m = 1,maux
        auxvf(m) = auxf(m,i)
        auxvc(m) = auxc(m,i)
    end do

    do m = 1,meqn
        qvf(m) = qf(m,i)
        qvc(m) = qc(m,i)
    end do


    !! Get face relative to ghost cell
    !! Left face of Cartesian grid   --> right edge of ghost cell
    !! Right face of Cartesian grid  --> left edge of ghost cell
    !! Bottom face of Cartesian grid --> top edge of ghost cell
    !! Top face of Cartesian grid    --> bottom edge of ghost cell
    if (idir .eq. 0) then
        iface_cell = 1-iface    !! Swap left and right edges
    else
        iface_cell = 5-iface    !! Swap bottom and top edges
    endif

    !! Call user defined function to compute fluxes.  The resulting
    !! flux should be those projected onto face 'iface_cell' and 
    !! scaled by edgelength/dx or edgelength/dy.
    call rpn2qad_flux(meqn,maux,idir,iface_cell,qvf,auxvf,fluxf)
    call rpn2qad_flux(meqn,maux,idir,iface_cell,qvc,auxvc,fluxc)

    do m = 1,meqn
        fd = fluxf(m) - fluxc(m)
        apdq(m,i) = 0.5*fd  
        amdq(m,i) = 0.5*fd
    end do

  ENDDO

END SUBROUTINE rpn2qad_new

