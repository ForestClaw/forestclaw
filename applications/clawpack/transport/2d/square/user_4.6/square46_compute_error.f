      subroutine square46_compute_error(blockno, mx,my,mbc,meqn,
     &      dx,dy,xlower,ylower,t,q,error, soln)
      implicit none

      integer mx,my,mbc,meqn, blockno
      double precision dx, dy, xlower, ylower, t
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision error(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision soln(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer*8 cont, fclaw_map_get_context

      integer example
      common /example_comm/ example

      integer i,j
      double precision xc,yc, qexact
      double precision xc1, yc1, zc1
      integer flow_flag

      cont = fclaw_map_get_context()

      if (example .eq. 0) then
c          # Divergent flow         
           flow_flag = 0
      else
c          # non-divergent flow         
           flow_flag = 1
      endif

c     # Assume a single field variable only
      do j = 1,my
         yc = ylower + (j-0.5)*dy
         do i = 1,mx
            xc = xlower + (i-0.5)*dx

            if (t .eq. 0) then
               soln(i,j,1) = q(i,j,1)
            else
c              # Map computational coordinates to physical coordinates
c              # In this case, physical and computational are the same.
               call fclaw_map_2d_c2m(cont,blockno,xc,yc,
     &                              xc1,yc1,zc1)
               soln(i,j,1) = qexact(xc1,yc1,t,flow_flag)
            endif
            error(i,j,1) = q(i,j,1) - soln(i,j,1)
         enddo
      enddo


      end
