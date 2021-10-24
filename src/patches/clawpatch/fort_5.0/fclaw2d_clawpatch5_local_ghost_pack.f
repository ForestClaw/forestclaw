c> @file
c> ghost pack routine for clawpack 5

c--------------------------------------------------------------------
c> @brief @copybrief ::clawpatch_fort_local_ghost_pack_t
c>
c> Implementation for clawpack 5.
c>
c> @details @copydetails ::clawpatch_fort_local_ghost_pack_t
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch5_fort_local_ghost_pack(mx,my,mbc,
     &      meqn, mint,qdata,area,qpack,psize,packmode,ierror)

      implicit none
      integer mx,my,mbc,meqn,psize, mint
      integer packmode, ierror
      double precision qdata(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision qpack(psize)

      integer packq, unpackq, packarea, unpackarea
      parameter(packq = 0, unpackq = 1, packarea = 2,
     &      unpackarea = 3)

      integer i,j,mq,k, ibc,jbc, kfinal
      integer nghost
      logical packdata

      ierror = 0
      if (packmode .ne. packq .and.
     &      packmode .ne. unpackq .and.
     &      packmode .ne. packarea .and.
     &      packmode .ne. unpackarea) then
         ierror = 1
         return
      endif
      packdata = packmode .eq. packq .or. packmode .eq. packarea

      nghost = 2
      k = 1
c     # Face 0
      do mq = 1,meqn
         do j = 1-nghost,my-mint
            do ibc = 1-nghost,mint
               if (packdata) then
                  qpack(k) = qdata(mq,ibc,j)
               else
                  qdata(mq,ibc,j) = qpack(k)
               endif
               k = k + 1
            enddo
         enddo

c        # Face 2
         do jbc = 1-nghost,mint
            do i = mint+1,mx+nghost
               if (packdata) then
                  qpack(k) = qdata(mq,i,jbc)
               else
                  qdata(mq,i,jbc) = qpack(k)
               endif
               k = k + 1
            enddo
         enddo

c        # Face 1
         do j = mint+1,my+nghost
            do ibc = mx-mint+1,mx+nghost
               if (packdata) then
                  qpack(k) = qdata(mq,ibc,j)
               else
                  qdata(mq,ibc,j) = qpack(k)
               endif
               k = k + 1
            enddo
         enddo

c        # Face 3
         do jbc = my-mint+1,my+nghost
            do i = 1-nghost,mx-mint
               if (packdata) then
                  qpack(k) = qdata(mq,i,jbc)
               else
                  qdata(mq,i,jbc) = qpack(k)
               endif
               k = k + 1
            enddo
         enddo
      enddo

      kfinal = k - 1
      if (packmode .ne. packarea .and. packmode .ne. unpackarea) then
         if (kfinal .ne. psize) then
            write(6,*) 'kfinal = ',kfinal
            write(6,*) 'psize = ',psize
            ierror = 2
         endif
         return
      endif


c     # Face 0
      do j = 1-nghost,my-mint
         do ibc = 1-nghost,mint
            if (packdata) then
               qpack(k) = area(ibc,j)
            else
               area(ibc,j) = qpack(k)
            endif
            k = k + 1
         enddo
      enddo

c     # Face 2
      do jbc = 1-nghost,mint
         do i = mint+1,mx+nghost
            if (packdata) then
               qpack(k) = area(i,jbc)
            else
               area(i,jbc) = qpack(k)
            endif
            k = k + 1
         enddo
      enddo

c     # Face 1
      do j = mint+1,my+nghost
         do ibc = mx-mint+1,mx+nghost
            if (packdata) then
               qpack(k) = area(ibc,j)
            else
               area(ibc,j) = qpack(k)
            endif
            k = k + 1
         enddo
      enddo

c     # Face 3
      do jbc = my-mint+1,my+nghost
         do i = 1-nghost,mx-mint
            if (packdata) then
               qpack(k) = area(i,jbc)
            else
               area(i,jbc) = qpack(k)
            endif
            k = k + 1
         enddo
      enddo

      kfinal = k-1
      if (kfinal .ne. psize) then
         ierror = 2
      endif

      end


c--------------------------------------------------------------------
c> @brief Set the boundaries (faces) to a value
c>
c> @param[in] mx, my the number of cells in the x and y directions
c> @param[in] mbc the number of ghost cells
c> @param[in] meqn the number of equations
c> @param[in,out] q the solution
c> @param[in] val the value to set the boundaries to
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch5_set_boundary_to_value(mx,my,mbc,
     &      meqn,q,val)
      implicit none

      integer mx,my,mbc,meqn
      double precision q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision val

      integer i,j,ibc,jbc,mq

c     # set boundary in strips, as with the packing routines

c     # Face 0
      do mq = 1,meqn
         do j = 1-mbc,my
            do ibc = 1,mbc
               q(mq,1-ibc,j) = val
            enddo
         enddo

c        # Face 2
         do jbc = 1,mbc
            do i = 1,mx+mbc
               q(mq,i,1-jbc) = val
            enddo
         enddo

c        # Face 1
         do j = 1,my+mbc
            do ibc = 1,mbc
               q(mq,mx+ibc,j) = val
            enddo
         enddo

c        # Face 3
         do jbc = 1,mbc
            do i = 1-mbc,mx
               q(mq,i,my+jbc) = val
            enddo
         enddo
      enddo

      end

c--------------------------------------------------------------------
c> @brief Set the corners to a value
c>
c> @param[in] mx, my the number of cells in the x and y directions
c> @param[in] mbc the number of ghost cells
c> @param[in] meqn the number of equations
c> @param[in,out] q the solution
c> @param[in] val the value to set the boundaries to
c--------------------------------------------------------------------
      subroutine fclaw2d_clawpatch5_set_corners_to_value(mx,my,mbc,
     &      meqn, q,value)
      implicit none

      integer mx,my,mbc,meqn
      double precision value
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer mq,ibc,jbc

      do mq = 1,meqn
         do ibc = mbc,mbc
            do jbc = mbc,mbc
               q(mq,1-ibc,1-jbc) = value
               q(mq,mx+ibc,1-jbc) = value
               q(mq,mx+ibc,my+jbc) = value
               q(mq,1-ibc,my+jbc) = value
            enddo
         enddo
      enddo


      end
