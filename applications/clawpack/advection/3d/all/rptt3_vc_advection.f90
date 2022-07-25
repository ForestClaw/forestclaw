subroutine clawpack46_rptt3(ixyz,icoor,imp,impt,maxm,meqn,mwaves,& 
    maux,mbc,mx,ql,qr,aux1,aux2,aux3,bsasdq,cmbsasdq,cpbsasdq)
    ! Riemann solver in the transverse direction for the
    ! advection equations.
    !
    ! On input,
    !
    !    ql,qr is the data along some one-dimensional slice, as in rpn3
    !         This slice is
    !             in the x-direction if ixyz=1,
    !             in the y-direction if ixyz=2, or
    !             in the z-direction if ixyz=3.
    !
    !    bsasdq is an array of flux differences that result from a
    !         transverse splitting (a previous call to rpt3).  
    !         This stands for B^* A^* \Dq but could represent any of 
    !         6 possibilities, e.g.  C^* B^* \Dq, as specified by ixyz
    !         and icoor (see below).
    !         Moreover, each * represents either + or -, as specified by
    !         imp and impt.
    !
    !    ixyz indicates the direction of the original Riemann solve,
    !         called the x-like direction in the table below:
    !
    !               x-like direction   y-like direction   z-like direction
    !      ixyz=1:        x                  y                  z
    !      ixyz=2:        y                  z                  x
    !      ixyz=3:        z                  x                  y
    !
    !    icoor indicates direction in which the transverse solve should
    !         be performed.
    !      icoor=2: split in the y-like direction.
    !      icoor=3: split in the z-like direction.
    !
    !    For example,
    !        ixyz=1, icoor=3 means bsasdq=B^*A^*\Dq, and should be
    !                        split in z into
    !                           cmbsasdq = C^-B^*A^*\Dq,
    !                           cpbsasdq = C^+B^*A^*\Dq.
    !
    !        ixyz=2, icoor=3 means bsasdq=C^*B^*\Dq, and should be
    !                        split in x into
    !                           cmbsasdq = A^-C^*B^*\Dq,
    !                           cpbsasdq = A^+C^*B^*\Dq.
    !
    !    The parameters imp and impt are generally needed only if aux
    !    arrays are being used, in order to access the appropriate
    !    variable coefficients:
    !
    !    imp =  1 if bsasdq = B^*A^- \Dq, a left-going flux difference
    !           2 if bsasdq = B^*A^+ \Dq, a right-going flux difference
    !    impt = 1 if bsasdq = B^-A^* \Dq, a down-going flux difference
    !           2 if bsasdq = B^+A^* \Dq, an up-going flux difference
    !
    !    aux2(:,:,2) is a 1d slice of the aux array along the row
    !                 where the data ql, qr lie.   
    !    aux1(:,:,2) and aux3(:,:,2) are neighboring rows in the 
    !                 y-like direction
    !    aux2(:,:,1) and aux2(:,:,3) are neighboring rows in the 
    !                z-like direction

    implicit none

    integer, intent(in) :: ixyz, icoor, maxm, meqn, mwaves, mbc, mx, maux, imp, impt
    double precision, intent(in) :: ql, qr, bsasdq, aux1, aux2, aux3
    double precision, intent(out) :: cmbsasdq, cpbsasdq

    dimension       ql(meqn,1-mbc:maxm+mbc)
    dimension       qr(meqn,1-mbc:maxm+mbc)
    dimension   bsasdq(meqn,1-mbc:maxm+mbc)
    dimension cmbsasdq(meqn,1-mbc:maxm+mbc)
    dimension cpbsasdq(meqn,1-mbc:maxm+mbc)
    dimension     aux1(maux,1-mbc:maxm+mbc,3)
    dimension     aux2(maux,1-mbc:maxm+mbc,3)
    dimension     aux3(maux,1-mbc:maxm+mbc,3)

    integer manifold
    common /com_manifold/ manifold

    integer :: i, i1, iuvw, mcapa

    ! set iuvw = 1 for u, 2 for v, 3 for w component of velocity
    ! depending on transverse direction:
    mcapa = manifold
    iuvw = ixyz + icoor - 1
    if (iuvw.gt.3) iuvw = iuvw-3

    do i=2-mbc,mx+mbc
        i1 = i-2+imp    !#  =  i-1 for amdq,  i for apdq
        if (icoor == 2) then
            !! double-transverse dir. is the y-like direction
            if (impt == 1) then
                cmbsasdq(1,i) = dmin1(aux2(iuvw + mcapa,i1,1),0.d0)*bsasdq(1,i)
                cpbsasdq(1,i) = dmax1(aux3(iuvw + mcapa,i1,1),0.d0)*bsasdq(1,i)
            elseif (impt == 2) then
                cmbsasdq(1,i) = dmin1(aux2(iuvw + mcapa,i1,3),0.d0)*bsasdq(1,i)
                cpbsasdq(1,i) = dmax1(aux3(iuvw + mcapa,i1,3),0.d0)*bsasdq(1,i)
            endif
        else
            !! double-transverse dir. is the z-like direction
            if (impt == 1) then
                !! bmasdq is split into cmbmasdq and cpbmasdq
                cmbsasdq(1,i) = dmin1(aux1(iuvw + mcapa,i1,2),0.d0)*bsasdq(1,i)
                cpbsasdq(1,i) = dmax1(aux1(iuvw + mcapa,i1,3),0.d0)*bsasdq(1,i)
            elseif (impt == 2) then
                !! bpasdq is split into cmbpasdq and cpbpasdq
                cmbsasdq(1,i) = dmin1(aux3(iuvw + mcapa,i1,2),0.d0)*bsasdq(1,i)
                cpbsasdq(1,i) = dmax1(aux3(iuvw + mcapa,i1,3),0.d0)*bsasdq(1,i)
            endif
        endif
    end do

    return
end subroutine clawpack46_rptt3
