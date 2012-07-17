      subroutine average2(mxc,myc,mxf,myf,mbc,meqn,qcoarse,qfine,
     &      kappacoarse,kappafine,refratio)
      implicit none

      integer mxc,myc,mxf,myf,mbc,meqn, refratio

      double precision qcoarse(mxc,myc,meqn)
      double precision   qfine(1-mbc:mxf + mbc,1-mbc:myf+mbc,meqn)

      double precision kappacoarse(mxc,myc)
      double precision kappafine(mxf,myf)

      integer i,j,m, ir, jr, ifine,jfine
      double precision ref_scale, sum, kappaf, kappac,maxdiff,
     &      sum1

      ref_scale = 1.d0/(refratio*refratio)
      maxdiff = 0.d0
      do i = 1,mxc
         do j = 1,myc
            do m = 1,meqn
               sum = 0.d0
               sum1 = 0.d0
               do ir = 1,refratio
                  do jr = 1,refratio
                     ifine = (i-1)*refratio + ir
                     jfine  = (j-1)*refratio + jr
                     kappaf = kappafine(ifine,jfine)
                     sum = sum + qfine(ifine,jfine,m)*kappaf
                     sum1 = sum1 + kappaf
                  enddo
               enddo
               kappac = kappacoarse(i,j)
               maxdiff = dmax1(maxdiff,abs(sum1*ref_scale - kappac))
               qcoarse(i,j,m) = sum*ref_scale/kappac
            enddo !! end m
         enddo  !! end i
      enddo  !! end j

c      write(6,'(A,E24.16)') 'maxdiff = ',maxdiff

      end
