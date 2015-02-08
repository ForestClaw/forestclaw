c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
c
c     # set auxiliary arrays 
c     # aux(i,j,1) = y coordinate of cell center for cylindrical source terms
c
c     
      implicit double precision (a-h,o-z)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)
c
      do 20 i=1-mbc,mx+mbc
         do 10 j=1-mbc,my+mbc
            yj = ylower + (j-0.5d0)*dy
            aux(i,j,1) = yj
   10       continue
   20    continue
c
       return
       end
