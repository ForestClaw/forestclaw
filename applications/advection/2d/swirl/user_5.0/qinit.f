      subroutine clawpack5_qinit(meqn,mbc,mx,my,
     &      xlower,ylower,
     &      dx,dy,q,maux,aux)

       implicit double precision (a-h,o-z)
       dimension q(meqn,1-mbc:mx+mbc, 1-mbc:my+mbc)
       dimension aux(maux,1-mbc:mx+mbc, 1-mbc:my+mbc)
c
       do 20 i=1-mbc,mx+mbc
          xi = xlower + (i-0.5d0)*dx
          do 20 j=1-mbc,my+mbc
             yj = ylower + (j-0.5d0)*dy
             if (xi.lt.0.5d0) then
c            if (xi.lt.0.5d0 .and. xi.gt.0.1d0 .and. yj.gt.0.1d0 .and.
c    &           yj.lt.0.3d0) then
                     q(1,i,j) = 1.d0
                   else
                     q(1,i,j) = 0.d0
                   endif
  20         continue
       return
       end
