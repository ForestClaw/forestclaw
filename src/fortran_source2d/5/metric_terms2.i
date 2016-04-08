C -*- FORTRAN -*-

c     # Include this file, or cut and paste dimensioning to get
c     # correct dimensioning in Fortran

c     # Cell centered quantities
      double precision          xp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision          yp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision          zp(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision        area(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision   curvature(-mbc:mx+mbc+1,-mbc:my+mbc+1)
      double precision surfnormals(3,-mbc:mx+mbc+1,-mbc:my+mbc+1)

c     # Node based quantities
      double precision xd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision yd(-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision zd(-mbc:mx+mbc+2,-mbc:my+mbc+2)

c     # Edge based quantities
      double precision     xnormals(3,-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision     ynormals(3,-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision    xtangents(3,-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision    ytangents(3,-mbc:mx+mbc+2,-mbc:my+mbc+2)
      double precision edge_lengths(2,-mbc:mx+mbc+2,-mbc:my+mbc+2)
