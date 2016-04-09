c     ============================================
      subroutine setaux(mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays 

c     #   aux(1,i,j) is edge velocity at "left" boundary of grid point (i,j)
c     #   aux(2,i,j) is edge velocity at "bottom" boundary of grid point (i,j)
c     #   aux(3,i,j) is kappa if a mapped grid is used.

c
c     
      implicit double precision (a-h,o-z)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
c
c     # constant velocities which are used if tperiod=0 is specified
c     # in setprob.data

      do 20 i=1-mbc,mx+mbc
         do 20 j=1-mbc,my+mbc

c           # coordinates of lower left corner of grid cell:
            xll = xlower + (i-1)*dx
            yll = ylower + (j-1)*dy

c           # difference stream function psi to get normal velocities:
            aux(1,i,j) = -(psi(xll, yll+dy) - psi(xll,yll)) / dy
            aux(2,i,j) =  (psi(xll+dx, yll) - psi(xll,yll)) / dx
   20       continue

c
       return
       end
