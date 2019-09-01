c
c
c =========================================================
      subroutine limiter3(maxmx,meqn,mwaves,mbc,mx,wave,s,
     &                   dtdx,phi,mthlim)
c =========================================================
c
c     # this routine limits the high-order correction term phi
c
c
      implicit double precision (a-h,o-z)
      dimension   wave(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension    dtdx(1-mbc:maxmx+mbc)
      dimension      s(1-mbc:maxmx+mbc, mwaves)
      dimension    phi(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension mthlim(mwaves)
      
      dimension      d(1-mbc:maxmx+mbc, meqn, mwaves)

c
c
c     ===============================================================
c     # TVD limiter
c     ===============================================================
c
      do 10 i=1,mx
        do 20 m=1,meqn
            do 30 mw=1,mwaves
              if(mthlim(mw) .eq. 1) then
               if(s(i,mw) .ge. 0d0) then
                phi(i,m,mw) = dmax1(0.d0,
     &          dmin1(2d0/(1d0-dtdx(i)*dabs(s(i,mw))),
     &          phi(i,m,mw),2.*wave(i-1,m,mw)
     &          /(dtdx(i)*dabs(s(i,mw))*wave(i,m,mw))))
               else
                phi(i,m,mw) = dmax1(0.d0,
     &          dmin1(2d0/(1d0-dtdx(i)*dabs(s(i,mw))),
     &          phi(i,m,mw),2.*wave(i+1,m,mw)
     &          /(dtdx(i)*dabs(s(i,mw))*wave(i,m,mw))))
               endif
              endif
   30          continue
   20       continue
   10    continue


c
      return
      end
