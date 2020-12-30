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
      double precision fp(-1:1,-1:1), qavg(-1:1,-1:1)
      double precision qp, qlap, flap

      integer i,j,ii,jj
c
c     computes point values of the conserved quantities    
      do i=2-mbc, mx+mbc-1
         do j=2-mbc, my+mbc-1
            do ii = -1,1
                do jj = -1,1
                  qavg(ii,jj) = q(i+ii,j+jj,1)
                end do
            end do
            qlap = (qavg(-1,0) + qavg(1,0) + qavg(0,-1) + 
     &                    qavg(0,1) - 4*qavg(0,0))

c           # Compute pointwise from average values          
            aux(i,j,1) = qavg(0,0) - qlap/24.
        end do
      end do

c     # Use cell-averaged values instead of pointwise values at 
c     # last layers of ghost cells (where we cannot compute the
c     # Laplacian.)
c     #
c     # Instead, use : 
c     # (2, -5, 4, -1) for second order one-sided derivative 
c     # approximation to second derivative 
      do j = 1-mbc,my+mbc
          aux(1-mbc,j,1) = q(1-mbc,j,1)
          aux(mx+mbc,j,1) = q(mx+mbc,j,1)
      enddo 

      do i = 1-mbc,mx+mbc
          aux(i,1-mbc,1) = q(i,1-mbc,1)
          aux(i,my+mbc,1) = q(i,my+mbc,1)
      enddo 
     

      do i=2-mbc, mx+mbc-1
         do j=2-mbc, my+mbc-1
            do ii = -1,1
                do jj = -1,1
                  qp = aux(i+ii,j+jj,1)
                  fp(ii,jj) = 0.5*qp**2
                end do
            end do
c           # Compute the average flux from pointwise values            
            flap = (fp(-1,0) + fp(1,0) + fp(0,-1) + 
     &                  fp(0,1) - 4*fp(0,0))
            aux(i,j,2) = fp(0,0) + flap/24.
        end do
      end do
            

      return
      end
