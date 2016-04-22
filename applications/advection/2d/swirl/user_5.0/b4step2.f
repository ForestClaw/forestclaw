c     ============================================
      subroutine b4step2(mbc,mx,my,meqn,q,
     &            xlower,ylower,dx,dy,time,dt,maux,aux)
c     ============================================
c
c     # called before each call to step
c     # use to set time-dependent aux arrays or perform other tasks.
c
c     # make velocity time dependent, reversing flow.
c
c     
      implicit double precision (a-h,o-z)
      dimension q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      dimension aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      common /comvt/ tperiod,pi2
c
      if (tperiod .eq. 0.d0) then
c        # special case --- indication that velocities specified in 
c        # setaux should be used for all time.
         return
         endif

       vt = dcos(pi2*(time+dt/2.d0)/tperiod)
c
       do 200 i=1-mbc,mx+mbc
          do 100 j=1-mbc,my+mbc
c           # coordinates of lower left corner of grid cell:
            xll = xlower + (i-1)*dx
            yll = ylower + (j-1)*dy

c           # difference stream function psi to get normal velocities:
            aux(1,i,j) = -(psi(xll, yll+dy) - psi(xll,yll)) / dy
            aux(2,i,j) =  (psi(xll+dx, yll) - psi(xll,yll)) / dx
c
c           # multiply by time-factor:
            aux(1,i,j) = vt * aux(1,i,j)
            aux(2,i,j) = vt * aux(2,i,j)
  100        continue
  200     continue
       return
       end
