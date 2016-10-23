      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
       implicit double precision (a-h,o-z)
       dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
c
       do 20 i = 1-mbc,mx+mbc
          xi = xlower + (i-0.5d0)*dx
          do 20 j=1-mbc,my+mbc
             yj = ylower + (j-0.5d0)*dy
             if (xi.gt.0.1d0 .and. xi.lt.0.6d0 .and.
     &           yj.gt.0.1d0 .and. yj.lt.0.6d0) then
                     q(i,j,1) = 1.d0
                   else
                     q(i,j,1) = 0.1d0
                   endif
  20         continue
       return
       end
