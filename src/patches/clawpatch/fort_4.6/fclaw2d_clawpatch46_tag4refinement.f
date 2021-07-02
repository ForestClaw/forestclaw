      subroutine fclaw2d_clawpatch46_fort_tag4refinement(mx,my,mbc,
     &      meqn, xlower,ylower,dx,dy,blockno,
     &      q, tag_threshold, init_flag, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch, init_flag
      integer blockno
      double precision xlower, ylower, dx, dy
      double precision tag_threshold
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax

      logical exceeds_th, fclaw2d_clawpatch_minmax_exceeds_th
      integer ii,jj
      double precision xc,yc,quad(-1:1,-1:1)

c     # Assume that we won't refine      
      tag_patch = 0

c     # Default : Refinement based only on first variable in system.  
c     # Users can modify this by creating a local copy of this routine
c     # and the corresponding tag4coarsening routine.
      mq = 1

      qmin = q(1,1,mq)
      qmax = q(1,1,mq)
      do j = 1-mbc,my+mbc
         do i = 1-mbc,mx+mbc
            xc = xlower + (i-0.5)*dx
            yc = ylower + (j-0.5)*dy
            qmin = min(qmin,q(i,j,mq))
            qmax = max(qmax,q(i,j,mq))
            do ii = -1,1
               do jj = -1,1
                  quad(ii,jj) = q(i+ii,j+jj,mq)
               end do
            end do
            exceeds_th = fclaw2d_clawpatch_minmax_exceeds_th(
     &             blockno, q(i,j,mq),qmin,qmax,quad, dx,dy,xc,yc,
     &             tag_threshold)
            if (exceeds_th) then
c              # Refine this patch               
               tag_patch = 1
               return
            endif
         enddo
      enddo

      end