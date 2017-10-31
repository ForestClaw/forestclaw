C -*- FORTRAN -*-

c     # Include this file, or cut and paste dimensioning to get
c     # correct dimensioning in Fortran

c     # Cell centered quantities
      double precision          xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision          yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision          zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision        area(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   curvature(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision surfnormals(-mbc:mx+mbc+1,-mbc:my+mbc+1,3)

c     # Node based quantities
      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

c     # Edge based quantities
      double precision     xnormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision     ynormals(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision    xtangents(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision    ytangents(-mbc:mx+mbc+2,-mbc:my+mbc+2,3)
      double precision  edgelengths(-mbc:mx+mbc+2,-mbc:my+mbc+2,2)
