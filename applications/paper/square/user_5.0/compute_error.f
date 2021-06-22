      subroutine square5_compute_error(blockno, mx,my,mbc,meqn,
     &      dx,dy,xlower,ylower,t,q,error, soln)
      implicit none

      integer mx,my,mbc,meqn, blockno
      double precision dx, dy, xlower, ylower, t
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision error(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision soln(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)

      integer*8 cont, get_context

      integer example
      common /example_comm/ example

      integer i,j
      double precision xc,yc, qexact
      double precision xc1, yc1, zc1
      integer flow_flag

      cont = get_context()

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
               soln(1,i,j) = q(1,i,j)
            else
c              # Map computational coordinates to physical coordinates
c              # In this case, physical and computational are the same.
               call fclaw2d_map_c2m(cont,blockno,xc,yc,
     &                              xc1,yc1,zc1)
               soln(1,i,j) = qexact(xc1,yc1,t,flow_flag)
            endif
            error(1,i,j) = q(1,i,j) - soln(1,i,j)
         enddo
      enddo


      end
