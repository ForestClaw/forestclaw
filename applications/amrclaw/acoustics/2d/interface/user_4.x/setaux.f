c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays 
c     # variable coefficient acoustics
c     #  aux(i,j,1) = density rho in (i,j) cell
c     #  aux(i,j,2) = sound speed c in (i,j) cell
c
c     # Piecewise constant medium with single interface at x=0.5
c     # Density and sound speed to left and right are set in setprob.f

c
c     
      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, 2)
      common /comaux/ rhol,cl,rhor,cr
c

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
