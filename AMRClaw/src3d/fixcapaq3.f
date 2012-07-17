      subroutine fixcapaq3(mxc,myc,mzc,mxf,myf,mzf,mbc,meqn,
     &      qcoarse,kappacoarse,qfine,kappafine,refratio)

      implicit none

      integer mxc,myc,mzc,mxf,myf,mzf,mbc,meqn,refratio

      double precision qcoarse(0:mxc+1,0:myc+1,0:mzc+1,meqn)
      double precision kappacoarse(mxc,myc,mzc)

      double precision qfine(1-mbc:mxf+mbc,1-mbc:myf+mbc,
     &      1-mbc:mzf + mbc,meqn)
      double precision kappafine(mxf,myf,mzf)

      integer i,j, k, ir, jr, kr, ifine, jfine, kfine, m
      double precision kappaf, kappac, cons_diff, ref_scale, sum
      double precision sum1,maxdiff

      ref_scale = 1.d0/(refratio**3)
      maxdiff = 0.d0
      do i = 1,mxc
         do j = 1,myc
            do k = 1,mzc
               kappac = kappacoarse(i,j,k)
               do m = 1,meqn
                  sum = 0.d0
                  sum1 = 0.d0
                  do ir = 1,refratio
                     do jr = 1,refratio
                        do kr = 1,refratio
                           ifine = (i-1)*refratio + ir
                           jfine = (j-1)*refratio + jr
                           kfine = (k-1)*refratio + kr
                           kappaf = kappafine(ifine,jfine,kfine)
                           sum = sum + qfine(ifine,jfine,kfine,m)*kappaf
                           sum1 = sum1 + kappaf
                        enddo
                     enddo
                  enddo

                  cons_diff = qcoarse(i,j,k,m)*kappac - sum*ref_scale
                  maxdiff = max(maxdiff,abs(sum1*ref_scale - kappac))

                  do ir = 1,refratio
                     do jr = 1,refratio
                        do kr = 1,refratio
                           ifine  = (i-1)*refratio + ir
                           jfine  = (j-1)*refratio + jr
                           kfine  = (k-1)*refratio + kr
                           kappaf = kappafine(ifine,jfine,kfine)
                           qfine(ifine,jfine,kfine,m) =
     &                           qfine(ifine,jfine,kfine,m) +
     &                           cons_diff/kappaf
                        enddo
                     enddo
                  enddo
               enddo                    !! end of meqn
            enddo !! k loop
         enddo  !! j loop
      enddo !! i loop

c      write(6,'(A,E24.16)') 'fixcapaq : maxdiff = ',maxdiff

      end
