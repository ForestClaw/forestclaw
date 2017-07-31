      subroutine fc2d_geoclaw_local_ghost_pack_aux(mx,my,mbc,maux,
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
      integer nghost
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

      nghost = 2
      k = 1

      do mq = 1,maux
c        # Face 0
         do j = 1-nghost,my-mint
            do ibc = 1-nghost,mint
               if (packdata) then
                  auxpack(k) = auxdata(mq,ibc,j)
               else
                  auxdata(mq,ibc,j) = auxpack(k)
               endif
               k = k + 1
            enddo
         enddo

c        # Face 2
         do jbc = 1-nghost,mint
            do i = mint+1,mx+nghost
               if (packdata) then
                  auxpack(k) = auxdata(mq,i,jbc)
               else
                  auxdata(mq,i,jbc) = auxpack(k)
               endif
               k = k + 1
            enddo
         enddo

c        # Face 1
         do j = mint+1,my+nghost
            do ibc = mx-mint+1,mx+nghost
               if (packdata) then
                  auxpack(k) = auxdata(mq,ibc,j)
               else
                  auxdata(mq,ibc,j) = auxpack(k)
               endif
               k = k + 1
            enddo
         enddo

c        # Face 3
         do jbc = my-mint+1,my+nghost
            do i = 1-nghost,mx-mint
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