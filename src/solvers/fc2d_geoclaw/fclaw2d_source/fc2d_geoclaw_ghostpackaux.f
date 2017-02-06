      subroutine fc2d_geoclaw_fort_ghostpackaux(mx,my,mbc,maux,
     &           mint,auxdata,auxpack,auxsize,packmode,ierror)

      implicit none
      integer mx,my,mbc,maux,auxsize,mint
      integer packmode,ierror
      double precision auxdata(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
      double precision auxpack(auxsize)

      integer packq, unpackq, packarea, unpackarea
      parameter(packq = 0, unpackq = 1, packarea = 2,
     &          unpackarea = 3)

      integer i,j,mq,k, ibc,jbc, kfinal
      logical packdata

      ierror = 0
      if (packmode .ne. packq .and.
     &    packmode .ne. unpackq .and.
     &    packmode .ne. packarea .and.
     &    packmode .ne. unpackarea) then
         ierror = 1
         return
      endif

      packdata = packmode .eq. packq .or. packmode .eq. packarea

      k = 1

      do mq = 1,maux
c        # Face 0
         do j = 0,my-mint
            do ibc = 0,mint
               if (packdata) then
                  auxpack(k) = auxdata(mq,ibc,j)
               else
                  auxdata(mq,ibc,j) = auxpack(k)
               endif
               k = k + 1
            enddo
         enddo

c        # Face 2
         do jbc = 0,mint
            do i = mint+1,mx+1
               if (packdata) then
                  auxpack(k) = auxdata(mq,i,jbc)
               else
                  auxdata(mq,i,jbc) = auxpack(k)
               endif
               k = k + 1
            enddo
         enddo

c        # Face 1
         do j = mint+1,my+1
            do ibc = mx-mint+1,mx+1
               if (packdata) then
                  auxpack(k) = auxdata(mq,ibc,j)
               else
                  auxdata(mq,ibc,j) = auxpack(k)
               endif
               k = k + 1
            enddo
         enddo

c        # Face 3
         do jbc = my-mint+1,my+1
            do i = 0,mx-mint
               if (packdata) then
                  auxpack(k) = auxdata(mq,i,jbc)
               else
                  auxdata(mq,i,jbc) = auxpack(k)
               endif
               k = k + 1
            enddo
         enddo
      enddo

      kfinal = k - 1
      if (kfinal .ne. auxsize) then
         write(6,*) 'kfinal = ',kfinal
         write(6,*) 'auxsize = ',auxsize
         ierror = 2
      endif

      end