      subroutine clawpack46_setaux(maxmx,maxmy,mbc,mx,my,
     &      xlower,ylower,dx,dy,maux,aux)
      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, 2)
      common /comaux/ rhol,cl,rhor,cr

      do 30 j=1-mbc,my+mbc
       do 20 i=1-mbc,mx+mbc
          xcell = xlower + (i-0.5d0)*dx
          if (xcell .lt. 0.5d0) then
              aux(i,j,1) = rhol
              aux(i,j,2) = cl
            else
              aux(i,j,1) = rhor
              aux(i,j,2) = cr
            endif
   20     continue
   30    continue

       return
       end
