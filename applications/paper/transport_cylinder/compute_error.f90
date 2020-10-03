SUBROUTINE amrclaw_compute_error(level,naux,nvar,time)

    use amr_module, only: alloc, node, rnode, lstart, hxposs, hyposs, store1, storeaux, & 
               nestlevel, ndihi, ndjhi,ndilo,ndjlo,nestlevel, & 
               timemult, nghost, cornxlo,cornylo,mcapa, levelptr
    IMPLICIT NONE

    INTEGER naux, nvar, level
    DOUBLE PRECISION time

!!      double precision error(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
!!      double precision soln(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

    INTEGER mx,my,mbc,meqn, blockno, maux
    DOUBLE PRECISION dx, dy, xlower, ylower, t, err(3), total_area

    INTEGER mitot, mjtot, loc, locaux
    integer iaddq, iaddaux

    integer i,j,m, mptr

    if (level .ne. 1) then
        write(6,*) 'amrclaw_compute_error : level .ne. 1'
        stop
    endif

    !! convert from AMRClaw parameters to Clawpack parameters

    meqn  = nvar
    maux  = naux
    mbc   = nghost

    dx     = hxposs(level)
    dy     = hyposs(level)
    t      = time

    do m = 1,3
        err(m) = 0
    end do
    total_area = 0

!!  Loop over all grids at level 'level'

    mptr = lstart(level)
20  if (mptr .eq. 0) go to 85

        mx = node(ndihi,mptr) - node(ndilo,mptr) + 1
        my = node(ndjhi,mptr) - node(ndjlo,mptr) + 1

        xlower = rnode(cornxlo,mptr)
        ylower = rnode(cornylo,mptr)

        !! Compute location of q and aux in allocation array 'alloc'
        !! Use hints from 'conck.f'.

        loc    = node(store1,mptr)
        locaux = node(storeaux,mptr)

        !! Accumulate error
        CALL cylinder_compute_error(blockno,mx,my,mbc,meqn, maux, mcapa, & 
                                    dx,dy,xlower,ylower,t, alloc(locaux), & 
                                    alloc(loc),err,total_area)

        mptr = node(levelptr,mptr) !! Linked list? 
        GO TO 20

85  continue
    err(1) = err(1)/total_area
    err(2) = sqrt(err(2)/total_area)    

    do m = 1,meqn
        write(6,100) m, err(1), err(2), err(3)
    end do

100 format('error[',I1,']',3E12.4)

END SUBROUTINE amrclaw_compute_error


SUBROUTINE cylinder_compute_error(blockno,mx,my,mbc,meqn, maux, mcapa, & 
           dx,dy,xlower,ylower,t,aux,q,err,total_area)

    IMPLICIT NONE

    INTEGER mx,my,mbc,meqn, maux, blockno, mcapa
    DOUBLE PRECISION dx, dy, xlower, ylower, t
    DOUBLE PRECISION q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    DOUBLE PRECISION err(3), total_area

    INTEGER example
    COMMON /example_comm/ example    

    !!integer*8 cont, get_context

    INTEGER i,j,m
    DOUBLE PRECISION xc,yc, qexact
    DOUBLE PRECISION xc1, yc1, zc1, x,y

    DOUBLE PRECISION cylinder_divergence
    INTEGER flow_flag

    DOUBLE PRECISION soln, errval, dxdy, area

    IF (example .eq. 0) THEN
!!      # non-divergent flow         
        flow_flag = 0
    ELSE
!!       # Divergent flow         
         flow_flag = 1
    ENDIF

!!  # Assume a single field variable only
    dxdy = dx*dy
    DO j = 1,my
        DO i = 1,mx
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy

            IF (t .eq. 0) THEN
               soln = q(1,i,j)
            ELSE
               !! Pass computational coordinates to qexact
               soln = qexact(xc,yc,t,flow_flag)
            ENDIF

            errval = abs(q(1,i,j) - soln)
            IF (mcapa .eq. 0) THEN
               area = dxdy
            ELSE
               area = aux(mcapa,i,j)*dxdy
            ENDIF

            err(1) = err(1) + errval*area
            err(2) = err(2) + errval**2*area
            err(3) = max(err(3),errval)
            total_area = total_area + area
        END DO
    END DO


END SUBROUTINE cylinder_compute_error
