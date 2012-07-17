      subroutine fixcapaq2(mxc,myc,mxf,myf,mbc,meqn,qcoarse,
     &      kappacoarse,qfine,kappafine,refratio)
      implicit none

      integer mxc,myc,mxf,myf,mbc,meqn, refratio

      double precision qcoarse(0:mxc+1,0:myc+1,meqn)
      double precision kappacoarse(mxc,myc)

      double precision qfine(1-mbc:mxf+mbc,1-mbc:myf+mbc,meqn)
      double precision kappafine(mxf,myf)

      integer i,j,ir, jr, ifine, jfine, m
      double precision kappaf, kappac, cons_diff, sum, ref_scale,
     &      maxdiff, sum1

      ref_scale = 1.d0/real(refratio*refratio)
      maxdiff = 0.d0
      do i = 1,mxc
         do j = 1,myc
            kappac = kappacoarse(i,j)
            do m = 1,meqn
               sum = 0.d0
               sum1 = 0.d0
               do ir = 1,refratio
                  do jr = 1,refratio
                     ifine = (i-1)*refratio + ir
                     jfine = (j-1)*refratio + jr
                     kappaf = kappafine(ifine,jfine)
                     sum = sum + qfine(ifine,jfine,m)*kappaf
                     sum1 = sum1 + kappaf
                  enddo
               enddo

               cons_diff = qcoarse(i,j,m)*kappac - sum*ref_scale
               maxdiff = max(maxdiff,abs(sum1*ref_scale - kappac))

               do ir = 1,refratio
                  do jr = 1,refratio
                     ifine  = (i-1)*refratio + ir
                     jfine  = (j-1)*refratio + jr
                     kappaf = kappafine(ifine,jfine)
                     qfine(ifine,jfine,m) = qfine(ifine,jfine,m) +
     &                     cons_diff/kappaf
                  enddo
               enddo
            enddo  !! end of meqn
         enddo
      enddo

c      write(6,'(A,E24.16)') 'fixcapaq : maxdiff = ',maxdiff
      end
