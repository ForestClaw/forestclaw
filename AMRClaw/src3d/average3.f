      subroutine average3(mxc,myc,mzc,mxf,myf,mzf,mbc,meqn,qcoarse,
     &      qfine,kappacoarse,kappafine,refratio)

      integer mxc,myc,mzc,mxf,myf,mzf,mbc,meqn

      double precision qcoarse(mxc,myc,mzc,meqn)
      double precision qfine(1-mbc:mxf+mbc,1-mbc:myf+mbc,
     &      1-mbc:mzf+mbc,meqn)

      double precision kappacoarse(mxc,myc,mzc)
      double precision kappafine(mxf,myf,mzf)

      integer i,j,k,m, refratio, ir, jr, kr, ifine,jfine,kfine
      double precision ref_scale, sum, kappaf, kappac

      ref_scale = 1.d0/real(refratio**3)
      do i = 1,mxc
         do j = 1,myc
            do k = 1,mzc
               do m = 1,meqn
                  sum = 0.d0
                  do ir = 1,refratio
                     do jr = 1,refratio
                        do kr = 1,refratio
                           ifine  = (i-1)*refratio + ir
                           jfine  = (j-1)*refratio + jr
                           kfine  = (k-1)*refratio + kr
                           kappaf = kappafine(ifine,jfine,kfine)
                           sum = sum + qfine(ifine,jfine,kfine,m)*kappaf
                        enddo
                     enddo
                  enddo
                  kappac = kappacoarse(i,j,k)
                  qcoarse(i,j,k,m) = sum*ref_scale/kappac
               enddo
            enddo
         enddo
      enddo

      end
