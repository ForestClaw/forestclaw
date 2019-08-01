      subroutine periodic_compute_error(blockno, mx,my,mbc,meqn,
     &      dx,dy,xlower,ylower,t,q,error, soln)
      implicit none

      external f

      integer mx,my,mbc,meqn, blockno
      double precision dx, dy, xlower, ylower, t
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision soln(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

c      integer*8 cont, get_context

      integer mapping
      common /mapping_comm/ mapping

      integer example
      common /example_comm/ example  

      double precision ubar, vbar
      common /comrp/ ubar,vbar

      integer i,j,m, k, l
      double precision xi, yj, xk, yk, a(3,3), q0


      a(1,1) = 1d0
      a(1,2) = 4d0
      a(1,3) = 1d0
      a(2,1) = 4d0
      a(2,2) = 16d0
      a(2,3) = 4d0
      a(3,1) = 1d0
      a(3,2) = 4d0
      a(3,3) = 1d0

c      cont = get_context()


      do j = 1-mbc,my+mbc
         do  i = 1-mbc,mx+mbc
            soln(i,j,1) = 0d0
         end do
      end do

c     # Assume a single field variable only
      do j = 1,my
         do i = 1,mx
c            xc = xlower + (i-0.5)*dx
c            yc = ylower + (j-0.5)*dy

            if (t .eq. 0) then               
               soln(i,j,1) = q(i,j,1)
            else
               if (example .eq. 0) then
                  xi = xlower + (i-1d0)*dx
                  yj = ylower + (j-1d0)*dy
                  do k = 0,2
                     do l = 0,2     
                        xk = xi + k*dx/2 - ubar*t
                        yk = yj + l*dy/2 - vbar*t             
                        soln(i,j,1) = soln(i,j,1) + 
     &                      (1d0/36d0)*a(k+1,l+1)*q0(xk,yk)
                      end do
                   end do
                endif
             endif
             error(i,j,1) = q(i,j,1) - soln(i,j,1)
         enddo
      enddo


      end



