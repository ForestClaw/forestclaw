      subroutine fc3d_clawpack5_fort_tag4refinement(mx,my,mz,mbc,
     &      meqn,xlower,ylower,zlower,dx,dy,dz,blockno,
     &      q, tag_threshold, init_flag,tag_patch)
      implicit none

      integer mx,my,mz,mbc,meqn,tag_patch,init_flag
      integer blockno
      double precision xlower,ylower,zlower,dx,dy,dz
      double precision tag_threshold
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc,1-mbc:mz+mbc)

      integer i,j,k,mq
      double precision qmin, qmax

      tag_patch = 0

c     # Refine based only on first variable in system.
      mq = 1
      qmin = q(mq,1,1,1)
      qmax = q(mq,1,1,1)
      do k = 1-mbc,mz+mbc
         do j = 1-mbc,my+mbc
            do i = 1-mbc,mx+mbc
               qmin = min(q(mq,i,j,k),qmin)
               qmax = max(q(mq,i,j,k),qmax)
               if (qmax - qmin .gt. tag_threshold) then
                  tag_patch = 1
                  return
               endif
            enddo
         enddo
      enddo
      end
