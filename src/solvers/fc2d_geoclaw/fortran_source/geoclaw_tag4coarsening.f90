c     # Template function for setting coarsening criteria.  The
c     # user can copy this file to their directory.  To
c     # indicate that this file should be used, set :
c     #
c     #      fclaw2d_vtable_t vt;
c     #      /* .... */
c     #      vt.fort_tag4coarsening = &tag4coarsening;
c     #      fclaw2d_set_vtable(domain,&vt);
c     #
c     # in virtual tables (typically set in <application>_user.cpp, in a
c     # a routine link '<application>_link_solvers(domain)'
c     #
c     # See also 'tag4refinement.f'

      subroutine tag4coarsening(mx,my,mbc,meqn,
     &      xlower,ylower,dx,dy, blockno, q0, q1, q2, q3,
     &      coarsen_threshold, tag_patch)
      implicit none

      integer mx,my, mbc, meqn, tag_patch
      integer blockno
      double precision xlower(0:3), ylower(0:3), dx, dy
      double precision coarsen_threshold
      double precision q0(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q1(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q2(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision q3(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j, mq
      double precision qmin, qmax

      tag_patch = 0
      qmin = q(1,1,1)
      qmax = q(1,1,1)
      call get_minmax(mx,my,mbc,meqn,q0,qmin,qmax)
      call get_minmax(mx,my,mbc,meqn,q1,qmin,qmax)
      call get_minmax(mx,my,mbc,meqn,q2,qmin,qmax)
      call get_minmax(mx,my,mbc,meqn,q3,qmin,qmax)
      if (qmax - qmin .lt. coarsen_threshold) then
         tag_patch = 1
         return
      endif

      end

      subroutine get_minmax(mx,my,mbc,meqn,q,qmin,qmax)

      implicit none
      integer mx,my,mbc,meqn
      double precision qmin,qmax
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      integer i,j,mq

      mq = 1
      do i = 1,mx
         do j = 1,my
            qmin = min(q(i,j,mq),qmin)
            qmax = max(q(i,j,mq),qmax)
         enddo
      enddo

      end
