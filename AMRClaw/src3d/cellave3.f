      subroutine cellave3(xlow,ylow,zlow,dx,dy,dz,w)
      implicit none

      external fs
      double precision  xlow,ylow,zlow,dx,dy,dz, w

c     # cellave3.f
c
c     # This computes the volume fraction of a 3d mesh cell that lies
c     # on the negative side of a surface
c
c     # This is the main routine which the user calls
c     # It is assumed that the user has defined a level set function
c     # 'fdisc', which implicitly describes the surface which cuts the
c     # mesh.  'fdisc' should have the following form :
c
c         double precision function fdisc(x,y,z)
c
c     # fdisc >= 0 is the positive side of the surface; fdisc < 0
c     # is the negative side.
c
c     # This is the only routine in which 'fs' (and therefore 'fdisc')
c     # is called.   If the user has direct access to 'cube' (values of
c     # the level set at 8 corners of mesh cube) and 'edges', the
c     # connectivity matrix describing surface/edge intersections (see below),
C     # then then the user can call 'cellave3_main' directly.  Doing so
c     # should significantly reduce computational time for situations
c     # in which cell averages are computed more than once.
c
c     # As a pre-processing step (for initialization purposes, for example)
c     # calling 'cellave3' for each mesh cell should be okay.

      integer icube(0:7)
      double precision edges(0:7,0:7)

      integer ipt, ipt_cnt, iedge, ipt1, ipt2, i,j
      double precision xl,xh,yl,yh,zl,zh, tol, s
      double precision corners(0:7,3)

      logical all_in, all_out

c     # functions
      double precision zeroin, fs

c     # For function fs
      common /comfs/ xl,xh,yl,yh,zl,zh

c     # Evaluate function at all 8 corners of cube.
      all_in = .true.
      all_out = .true.
      do ipt = 0,7
c        # xl,yl,zl,xh,yh,zh are corners that are passed in common
c        # block to function fs.
         xl = xlow + ibits(ipt,0,1)*dx
         yl = ylow + ibits(ipt,1,1)*dy
         zl = zlow + ibits(ipt,2,1)*dz
         xh = xl
         yh = yl
         zh = zl
         corners(ipt,1) = xl
         corners(ipt,2) = yl
         corners(ipt,3) = zl
         icube(ipt) = sign(1.0d0,fs(0.d0))
         all_in = all_in .and. icube(ipt) .lt. 0
         all_out = all_out .and. icube(ipt) .ge. 0
c         icube(ipt) = sign(1,fs(0.d0))   !! sign(1,0) = 1
c         all_in = all_in .and. icube(ipt) .gt. 0
c         all_out = all_out .and. icube(ipt) .lt. 0
      enddo

c     # Immediately eliminate those cubes which are not cut by the surface
      if (all_in) then
         w = 1.d0
         return
      elseif (all_out) then
         w = 0.d0
         return
      endif

c     # Determine intersections at cube edges.
      do i = 0,7
         do j = 0,7
            edges(i,j) = -1.d0
         enddo
      enddo

      tol = 1.0e-8
      ipt_cnt = 0
      do iedge = 0,11
         call get_endpoints(iedge,ipt1,ipt2)
         if (icube(ipt1)*icube(ipt2) .lt. 0) then
            xl = corners(ipt1,1)
            yl = corners(ipt1,2)
            zl = corners(ipt1,3)
            xh = corners(ipt2,1)
            yh = corners(ipt2,2)
            zh = corners(ipt2,3)
            s = zeroin(0.d0,1.d0,fs,tol)
            edges(ipt1,ipt2) = s
            edges(ipt2,ipt1) = 1.d0-s
         endif
      enddo

      call cellave3_main(icube, edges, w)

      end

c     # -----------------------------------------------------------
c     # cellave_main.f
c     #
c     # Main routine which computes the volume fraction, using only
c     # data in arrays 'cube' and 'edges'.
c     # Functions 'fs' (and 'fdisc') are not called from this routine
      subroutine cellave3_main(icube,edges,w)
      implicit none

      double precision  w

c     # Main info stored here
      integer icube(0:7)
      integer icorners(0:7,3)
      integer faces(0:5,0:5)
      double precision edges(0:7,0:7)

c     # intersection point list and adjacent faces list
c     # isect_cnt is length of list
c     # isect_list and face_list are not in any particular order
c     # But surf_list is an order list of points describing the
c     # surface
      integer isect_cnt, surf_cnt, poly_cnt
      integer surf_list(10)
      integer edge_list(10)
      integer face_list(10,2)

c     # Computing area of polygon at x = 1
      double precision isect_list(10,3)
      double precision ypoly(10), zpoly(10), area
      integer iv_list(5)

      integer ipt, iedge, ipt1, ipt2, i, j, k
      double precision s, v(3)

c     # For face matrix
      integer if1, if2
      integer if1_start, if1_new,if2_start,if2_new
      integer d1, d2
      logical done

c     # For triangulation of surface
      double precision p0(3), p1(3), p2(3)
      double precision dw, tri_normal(3), v1(3), v2(3)
      integer sgn, idir

c     # Get corner locations
      do ipt = 0,7
         icorners(ipt,1) = ibits(ipt,0,1)
         icorners(ipt,2) = ibits(ipt,1,1)
         icorners(ipt,3) = ibits(ipt,2,1)
      enddo

c     # Construct list of intersection points and adjacent faces
      isect_cnt = 0
      do iedge = 0,11
         call get_endpoints(iedge,ipt1,ipt2)
         s = edges(ipt1,ipt2)
         if (s .ge. 0.d0) then
            isect_cnt = isect_cnt + 1  !! This number shouldn't exceed 10

            do j = 1,3
               v(j) = icorners(ipt2,j) - icorners(ipt1,j)
               isect_list(isect_cnt,j) = icorners(ipt1,j) + s*v(j)
            enddo
            call get_faces(iedge,if1,if2)
            edge_list(isect_cnt) = iedge
            face_list(isect_cnt,1) = if1
            face_list(isect_cnt,2) = if2
         endif
      enddo

c     # Find area of intersection of surface with x = 1 face.
c     # plist returns four corners of face at x = 1
      call get_face_vertices(1,iv_list)
      iv_list(5) = iv_list(1)

c     # Traverse points on face
      ipt1 = iv_list(1)
      if (icube(ipt1) .lt. 0) then
c        # Add first corner to list of polygon points
         ypoly(1) = icorners(ipt1,2)
         zpoly(1) = icorners(ipt1,3)
         poly_cnt = 1
      else
         poly_cnt = 0
      endif

      do i = 2,5
         ipt2 = iv_list(i)

c        # First check edge between ipt1 and ipt2
c        # An intersection point will always be added to the list
         s = edges(ipt1,ipt2)
         if (s .ge. 0.d0) then
            poly_cnt = poly_cnt + 1
            v(2) = icorners(ipt2,2) - icorners(ipt1,2)
            v(3) = icorners(ipt2,3) - icorners(ipt1,3)
            ypoly(poly_cnt) = icorners(ipt1,2) + s*v(2)
            zpoly(poly_cnt) = icorners(ipt1,3) + s*v(3)
         endif

c        # Now check if corner point should be added to list
         if (icube(ipt2) .lt. 0) then
            poly_cnt = poly_cnt + 1
            ypoly(poly_cnt) = icorners(ipt2,2)
            zpoly(poly_cnt) = icorners(ipt2,3)
         endif
         ipt1 = ipt2
      enddo                             !! end loop over corners of face

      area = 0.d0
      if (poly_cnt .gt. 2) then
c        # Compute area of polygon on x = xlow+ival*dx face.
         ypoly(poly_cnt+1) = ypoly(1)
         zpoly(poly_cnt+1) = zpoly(1)
         do  i = 1,poly_cnt
            area = area +
     &            0.5d0*(zpoly(i)+zpoly(i+1))*(ypoly(i+1)-ypoly(i))
         enddo
      endif

c     # First contribution to volume fraction
      w = abs(area)

c     # ------------------------------------------------------
c     # Get the integral (x*nx) for each triangle that
c     # makes up irregular cut surface

c     # Reorder points on isect_list so that we traverse the surface
c     # in the proper order

c     # Connectivity matrix for faces.  Two faces are connected
c     # if there is a surface intersection point between them
      do i = 0,5
         do j = 0,5
            faces(i,j) = 0
         enddo
      enddo

      do k = 1,isect_cnt
        if1 = face_list(k,1)
        if2 = face_list(k,2)
        faces(if1,if2) = k
        faces(if2,if1) = k
      enddo

c     # Faces adjacent to starting edge value
      if1_start = face_list(1,1)
      if2_start = face_list(1,2)
      if1_new = if1_start
      if2_new = if2_start

c     # Traverse edges by successively checking rows and columns
c     # of the connectivity matrix for positive entries
      d1 = 1  !! search directions.
      d2 = 0
      surf_cnt = 0
      done = .false.
      do while (.not. done)
c        # loop over row (if1; d1 == 1) or column (if2; d2 = 1)
         do i = 1,5
            if1 = mod(if1_new+i*mod(d1,2),6)
            if2 = mod(if2_new+i*mod(d2,2),6)
c           # k will be index into isect_list, the list of
c           # intersection points
            k = faces(if1,if2)
            if (k .gt. 0) then
               surf_cnt = surf_cnt + 1
               surf_list(surf_cnt) = k
               if1_new = if1
               if2_new = if2
c              # Portland compiler doesn't like 'exit'
c               exit
               goto 10
            endif
         enddo
   10    continue
c        # if dN is even, it becomes odd, and vice. versa.
         d1 = d1 + 1
         d2 = d2 + 1
         done = ((if1 .eq. if1_start .and. if2 .eq. if2_start) .or.
     &           (if2 .eq. if1_start .and. if1 .eq. if2_start))
      enddo

c     #  Compute integral (x*nx) over each triangle

c     # p0 is root vertex for surface triangulation
      do j = 1,3
         p0(j) = isect_list(surf_list(1),j)
      enddo

c     # Form triangles successively : (p0,p1,p2), where p1, p2 are
c     # incremented at each cycle through loop
      do i = 2,surf_cnt - 1
         do j = 1,3
            p1(j) = isect_list(surf_list(i),j)
            p2(j) = isect_list(surf_list(i+1),j)
            v1(j) = p1(j) - p0(j)
            v2(j) = p2(j) - p0(j)
         enddo
         tri_normal(1) =   v1(2)*v2(3) - v1(3)*v2(2)
         tri_normal(2) = -(v1(1)*v2(3) - v1(3)*v2(1))
         tri_normal(3) =   v1(1)*v2(2) - v1(2)*v2(1)
         dw = (p0(1) + p1(1) + p2(1))/6.d0

c        # 'get_updir' returns a vector dir = sgn*(1,0,0), sgn*(0,1,0) or
c        # sgn*(0,0,1).  'idir' indicates where the '1' is.
c        # Dot product of dir with tri_normal should give us
c        # the orientation of the surface.  Furthermore, because
c        # current triangle (p0,p1,p2) has one vertex at edge 'iedge'
c        # tri_normal(idir+1) is non-zero. If 'iedge' were to lie exactly
c        # on the surface,  then both endpoints of edge 'iedge' are on
c        # same side of surface, (the '+' side), and there would
c        # be no intersection at that edge.

         call get_updir(icube,edge_list(i),idir,sgn)
         w = w + sgn*sign(1.0d0,tri_normal(idir+1))
     &         *tri_normal(1)*dw
      enddo

      end

      integer function get_edge_dir(iedge)
      implicit none

      integer iedge
      get_edge_dir = iedge/4
      end

      subroutine get_endpoints(iedge,ipt1,ipt2)
      implicit none

      integer get_edge_dir
      integer iedge, ipt1, ipt2, idir, n
      integer kvec(0:1,0:2), k, i

      data kvec /1, 2, 0, 2, 0, 1/

      idir = get_edge_dir(iedge)
      n = iedge - 4*idir

      ipt1 = 0
      ipt2 = 0
      ipt2 = ibset(ipt2,idir)
      do i = 0,1
         if (ibits(n,i,1) .eq. 1) then
            k = kvec(i,idir)
            ipt1 = ibset(ipt1,k)
            ipt2 = ibset(ipt2,k)
         endif
      enddo
      end


      subroutine get_faces(iedge, if1,if2)
      implicit none

      integer iedge, if1, if2, ipt1,ipt2, icnt
      integer b1, b2,ifaces(2), idir

      call get_endpoints(iedge,ipt1,ipt2)

      icnt = 0
      do idir = 0,2
         b1 = ibits(ipt1,idir,1)
         b2 = ibits(ipt2,idir,1)
         if (b1 .eq. b2) then
            icnt = icnt + 1
            ifaces(icnt) = 2*idir + b1
         endif
      enddo
      if1 = ifaces(1)
      if2 = ifaces(2)

      end

      subroutine get_face_vertices(iface,plist)
      implicit none

      integer idir,ival, plist(5), iface
      integer ip, ipt, b, ptmp

      idir = iface/2
      ival = iface - 2*idir

      ip = 0
      do ipt = 0,7
         b = ibits(ipt,idir,1)
         if (b .eq. ival) then
            ip = ip + 1
            plist(ip) = ipt
         endif
      enddo

c     # Swap last two elements so that points traverse face, rather
c     # than criss-cross it.
      ptmp = plist(3)
      plist(3) = plist(4)
      plist(4) = ptmp

      end

      subroutine get_updir(icube,iedge,idir,sgn)
      implicit none

      integer get_edge_dir
      integer ipt1, ipt2, iedge,j,idir, diff
      integer icube(0:7), sgn

      call get_endpoints(iedge,ipt1,ipt2)

      sgn = sign(1,icube(ipt2) - icube(ipt1))
      idir = get_edge_dir(iedge)

      end

      double precision function fs(s)
      implicit none

      double precision x,y,z, fdisc
      double precision xl,xh,yl,yh,zl,zh,s

      common /comfs/ xl,xh,yl,yh,zl,zh

      x = xl + s*(xh - xl)
      y = yl + s*(yh - yl)
      z = zl + s*(zh - zl)

      fs = fdisc(x,y,z)

      end



c     =================================================
      double precision function zeroin(ax,bx,f,tol)
c     =================================================
      implicit none
      external f
c
c      a zero of the function  f(x)  is computed in the interval ax,bx .
c      (Standard routine from netlib)
c
c  input..
c
c  ax     left endpoint of initial interval
c  bx     right endpoint of initial interval
c  f      function subprogram which evaluates f(x) for any x in
c         the interval  ax,bx
c  tol    desired length of the interval of uncertainty of the
c         final result ( .ge. 0.0)
c
c
c  output..
c
c  zeroin abcissa approximating a zero of  f  in the interval ax,bx
c
c
c      it is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
c  without  a  check.  zeroin  returns a zero  x  in the given interval
c  ax,bx  to within a tolerance  4*macheps*dabs(x) + tol, where macheps
c  is the relative machine precision.
c      this function subprogram is a slightly  modified  translation  of
c  the algol 60 procedure  zero  given in  richard brent, algorithms for
c  minimization without derivatives, prentice - hall, inc. (1973).
c
c
c
c  compute eps, the relative machine precision
c

      double precision ax,bx,tol, f, eps, tol1, a,b,c,d,e
      double precision fa,fb,fc, xm, p, q, r, s

      eps = 1.0
   10 eps = eps/2.0
      tol1 = 1.0 + eps
      if (tol1 .gt. 1.0) go to 10
c
c initialization
c
      a = ax
      b = bx
      fa = f(a)
      fb = f(b)
c
c begin step
c
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 if (dabs(fc) .ge. dabs(fb)) go to 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
c
c convergence test
c
   40 tol1 = 2.0*eps*dabs(b) + 0.5*tol
      xm = .5*(c - b)
      if (dabs(xm) .le. tol1) go to 90
      if (fb .eq. 0.0) go to 90
c
c is bisection necessary
c
      if (dabs(e) .lt. tol1) go to 70
      if (dabs(fa) .le. dabs(fb)) go to 70
c
c is quadratic interpolation possible
c
      if (a .ne. c) go to 50
c
c linear interpolation
c
      s = fb/fa
      p = 2.0*xm*s
      q = 1.0 - s
      go to 60
c
c inverse quadratic interpolation
c
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0*xm*q*(q - r) - (b - a)*(r - 1.0))
      q = (q - 1.0)*(r - 1.0)*(s - 1.0)
c
c adjust signs
c
   60 if (p .gt. 0.0) q = -q
      p = dabs(p)
c
c is interpolation acceptable
c
      if ((2.0*p) .ge. (3.0*xm*q - dabs(tol1*q))) go to 70
      if (p .ge. dabs(0.5*e*q)) go to 70
      e = d
      d = p/q
      go to 80
c
c bisection
c
   70 d = xm
      e = d
c
c complete step
c
   80 a = b
      fa = fb
      if (dabs(d) .gt. tol1) b = b + d
      if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
      fb = f(b)
      if ((fb*(fc/dabs(fc))) .gt. 0.0) go to 20
      go to 30
c
c done
c
   90 zeroin = b
      return
      end
