      double precision function cosine_bell_sum(x,y,z)
      implicit none

      double precision x,y,z

      double precision wcloc(3,2)
      common /location_parms/ wcloc

      double precision r, hmax, b, c
      common /cosinebell_parms/ r, hmax, b, c

      integer m,k
      double precision w(3), q,qv, cosbell

      q = 0
      do m = 1,2
         do k = 1,3
            w(k) = wcloc(k,m)
         end do
         qv = cosbell(x,y,z,w,r,hmax)
         q  = q  + qv
      enddo
      cosine_bell_sum = b + c*q

      end


      double precision function cosbell(x,y,z,w,r,hmax)
      implicit none

      double precision x,y,z,w(3),r,hmax
      double precision get_gc_distance, h1, ri

      double precision pi, pi2
      common /compi/ pi, pi2


      ri = get_gc_distance(x,y,z,w)

      if (ri .le. r) then
         h1 = hmax*(1 + cos(pi*ri/r))
      else
         h1 = 0.d0
      endif
      cosbell = h1

      end

