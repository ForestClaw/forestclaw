      subroutine clawpack46_qinit(maxmx,maxmy,meqn,mbc,mx,my,
     &      xlower,ylower,dx,dy,q,maux,aux)
       implicit double precision (a-h,o-z)
       dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)
c
       do 20 i=1-mbc,mx+mbc
          xi = xlower + (i-0.5d0)*dx
          do 20 j=1-mbc,my+mbc
             yj = ylower + (j-0.5d0)*dy
             r = dsqrt((xi-0.25d0)**2 + (yj-0.4d0)**2)
             q(i,j,1) = p0(r)
             q(i,j,2) = 0.d0
             q(i,j,3) = 0.d0
  20         continue
       return
       end
