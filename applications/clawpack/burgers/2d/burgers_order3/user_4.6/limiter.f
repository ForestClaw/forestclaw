c
c
c =========================================================
      subroutine limiter(maxmx,meqn,mwaves,mbc,mx,wave,s,
     &                   dtdx,phi,mthlim)
c =========================================================
c
c     # this routine limits the high-order correction term phi
c
c
      implicit double precision (a-h,o-z)
c     # extern variables
      dimension   wave(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension      s(1-mbc:maxmx+mbc, mwaves)
      dimension    phi(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension   dtdx(1-mbc:maxmx+mbc)
      dimension mthlim(mwaves)
c     # intern variables
      dimension      d(1-mbc:maxmx+mbc, meqn, mwaves)
      dimension    snu(1-mbc:maxmx+mbc, mwaves)
c
c     ===============================================================
c     # wave limiter
c     ===============================================================
c
      do 50 mw=1,mwaves
         if (mthlim(mw) .ge. 1) go to 50
         dotr = 0.d0
         do 40 i = 2-mbc, mx+mbc-1
            wnorm2 = 0.d0
            dotl = dotr
            dotr = 0.d0
            do 20 m=1,meqn
               wnorm2 = wnorm2 + wave(i,m,mw)**2
               dotr = dotr + wave(i,m,mw)*wave(i+1,m,mw)
   20          continue
            if (i.eq.0) go to 40
            if (wnorm2.eq.0.d0) go to 40
c
            if (s(i,mw) .gt. 0.d0) then
                wlimitr = philim(wnorm2, dotl, mthlim(mw))
              else
                wlimitr = philim(wnorm2, dotr, mthlim(mw))
              endif
c
            do 30 m=1,meqn
               phi(i,m,mw) = wlimitr * phi(i,m,mw)
   30          continue
   40       continue
   50    continue
c
c     # computing courant-number
      do 5 i=3-mbc,mx+mbc-2
         do 6 mw=1,mwaves
            if(s(i,mw) .ge. 0d0) then
                dtdxave=dtdx(i)
            else
                dtdxave=dtdx(i-1)
            endif
            snu(i,mw)=dtdxave*dabs(s(i,mw))
    6       continue
    5    continue
c
c     ===============================================================
c     # TVD limiter
c     ===============================================================
c
      do 12 i=1,mx
         do 22 m=1,meqn
            do 32 mw=1,mwaves
              if(mthlim(mw) .eq. 1) then
               if(s(i,mw) .ge. 0d0) then
               phi(i,m,mw)=dmax1(0d0,dmin1(2d0/(1d0-snu(i,mw)),
     &            phi(i,m,mw),
     &                 (2d0 * wave(i-1,m,mw))
     &                   /(snu(i,mw) * wave(i,m,mw)+1.d-15)
     &            * (1d0-snu(i-1,mw))/(1d0-snu(i,mw))
     &            ))
               else
               phi(i,m,mw)=dmax1(0d0,dmin1(2d0/(1d0-snu(i,mw)),
     &            phi(i,m,mw),
     &                 (2d0 * wave(i+1,m,mw))
     &                   /(snu(i,mw) * wave(i,m,mw)+1.d-15)
     &            * (1d0-snu(i+1,mw))/(1d0-snu(i,mw))
     &            ))
               endif
              endif
   32          continue
   22       continue
   12    continue
c
c
c     ===============================================================
c     # MP limiter
c     ===============================================================
c
c     # computing curvature
      do 13 i=1,mx
         do 23 m=1,meqn
            do 33 mw=1,mwaves
               d1=wave(i+1,m,mw)-wave(i,m,mw)
               d2=wave(i,m,mw)-wave(i-1,m,mw)
               d(i,m,mw)=dmax1(0d0, dmin1(d1,d2))
   33          continue
   23       continue
   13    continue
c
c     # computing MP-limiter
      do 14 i=1,mx
         do 24 m=1,meqn
            do 34 mw=1,mwaves
              if(mthlim(mw) .eq. 2) then
               if(s(i,mw) .ge. 0d0) then
                  r=wave(i-1,m,mw)/(wave(i,m,mw)+1.d-15)
               else
                  r=wave(i+1,m,mw)/(wave(i,m,mw)+1.d-15)
               endif
               phimd=0.5d0*(wave(i,m,mw)-d(i,m,mw))
               philc=1d0/snu(i,mw) * (r
     &           +  d(i-1,m,mw)/(wave(i,m,mw)+1.d-15))
               phimin=dmax1(dmin1(0d0,phimd),
     &                dmin1(0d0,2d0*r/snu(i,mw),
     &                      philc))
               phimax=dmin1(dmax1(2d0/(1d0-snu(i,mw)),phimd),
     &                dmax1(0d0,2d0*r/snu(i,mw),
     &                      philc))
               phi(i,m,mw)=dmax1(phimin,dmin1(phi(i,m,mw),phimax))
              endif
   34          continue
   24       continue
   14    continue
c
      return
      end
