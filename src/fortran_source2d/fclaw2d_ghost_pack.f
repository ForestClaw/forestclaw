      subroutine fclaw2d_ghost_pack(mx,my,mbc,meqn,
     &      mint,qdata,area,qpack,psize,packmode,
     &      ierror)

      implicit none
      integer mx,my,mbc,meqn,psize, mint
      integer pack_area, packmode, ierror
      double precision qdata(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision area(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision qpack(psize)

      integer packq, unpackq, packarea, unpackarea
      parameter(packq = 0, unpackq = 1, packarea = 2,
     &      unpackarea = 3)

      integer i,j,mq,k, ibc,jbc, kfinal

      ierror = 0
      if (packmode .ne. packq .and.
     &      packmode .ne. unpackq .and.
     &      packmode .ne. packarea .and.
     &      packmode .ne. unpackarea) then
         ierror = 1
         return
      endif


      k = 1
c     # Face 0
      do mq = 1,meqn
         do j = 1-mbc,my-mint
            do ibc = 1-mbc,mint
               if (packmode .eq. packq) then
                  qpack(k) = qdata(ibc,j,mq)
               elseif (packmode .eq. unpackq) then
                  qdata(ibc,j,mq) = qpack(k)
               endif
               k = k + 1
            enddo
         enddo

c        # Face 2
         do jbc = 1-mbc,mint
            do i = mint+1,mx+mbc
               if (packmode .eq. packq) then
                  qpack(k) = qdata(i,jbc,mq)
               elseif (packmode .eq. unpackq) then
                  qdata(i,jbc,mq) = qpack(k)
               endif
               k = k + 1
            enddo
         enddo

c        # Face 1
         do j = mint+1,my+mbc
            do ibc = mx-mint+1,mx+mbc
               if (packmode .eq. packq) then
                  qpack(k) = qdata(ibc,j,mq)
               elseif (packmode .eq. unpackq) then
                  qdata(ibc,j,mq) = qpack(k)
               endif
               k = k + 1
            enddo
         enddo

c        # Face 3
         do jbc = my-mint+1,my+mbc
            do i = 1-mbc,mx-mint
               if (packmode .eq. packq) then
                  qpack(k) = qdata(i,jbc,mq)
               elseif (packmode .eq. unpackq) then
                  qdata(i,jbc,mq) = qpack(k)
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
      do j = 1-mbc,my-mint
         do ibc = 1-mbc,mint
            if (packmode .eq. packarea) then
               qpack(k) = area(ibc,j)
            elseif (packmode .eq. unpackarea) then
               area(ibc,j) = qpack(k)
            endif
            k = k + 1
         enddo
      enddo

c     # Face 2
      do jbc = 1-mbc,mint
         do i = mint+1,mx+mbc
            if (packmode .eq. packarea) then
               qpack(k) = area(i,jbc)
            elseif (packmode .eq. unpackarea) then
               area(i,jbc) = qpack(k)
            endif
            k = k + 1
         enddo
      enddo

c     # Face 1
      do j = mint+1,my+mbc
         do ibc = mx-mint+1,mx+mbc
            if (packmode .eq. packarea) then
               qpack(k) = area(ibc,j)
            elseif (packmode .eq. unpackarea) then
               area(ibc,j) = qpack(k)
            endif
            k = k + 1
         enddo
      enddo

c     # Face 3
      do jbc = my-mint+1,my+mbc
         do i = 1-mbc,mx-mint
            if (packmode .eq. packarea) then
               qpack(k) = area(i,jbc)
            elseif (packmode .eq. unpackarea) then
               area(i,jbc) = qpack(k)
            endif
            k = k + 1
         enddo
      enddo

      kfinal = k-1;
      if (kfinal .ne. psize) then
         ierror = 2
      endif

      end
