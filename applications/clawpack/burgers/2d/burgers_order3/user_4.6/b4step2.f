c     ============================================
      subroutine clawpack46_b4step2(maxmx,maxmy,mbc,
     &            mx,my,meqn,q,
     &            xlower,ylower,dx,dy,t,dt,maux,aux)
c     ============================================
c
c     # called from claw2 before each call to step2.
c     # use to set time-dependent aux arrays or perform other tasks
c     # which must be done every time step.

c
c     # dummy routine 

       
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)
c
c     computes point values of the conserved quantities    
      do i=2-mbc, mx+mbc-1
         do j=2-mbc, my+mbc-1
            aux(i,j,1) = q(i,j,1)
     &             -(q(i-1,j,1)-2.d0*q(i,j,1)+q(i+1,j,1))/24d0
     &            - (q(i,j-1,1)-2.d0*q(i,j,1)+q(i,j+1,1))/24d0
         enddo
      enddo

c     # Use cell-averaged values instead of pointwise values at 
c     # last layers of ghost cells (where we cannot compute the
c     # Laplacian.)
      do j = 1-mbc,my+mbc
          aux(1-mbc,j,1) = q(1-mbc,j,1)
          aux(mx+mbc,j,1) = q(mx+mbc,j,1)
      enddo 

      do i = 1-mbc,mx+mbc
          aux(i,1-mbc,1) = q(i,1-mbc,1)
          aux(i,my+mbc,1) = q(i,my+mbc,1)
      enddo 
     
      return
      end
