      double precision function compute_slopes(sl,sr)
      implicit none

      double precision sl,sr, sc
      integer mth

c     # Use AMRClaw slopes  (use minimum in absolute value;  sign is
c     # chosen from centered (sc) slope
      sc = (sl + sr)/2.d0
      compute_slopes = min(abs(sl),abs(sr),abs(sc))*
     &      max(0.d0,sign(1.d0,sl*sr))*sign(1.d0,sc)

c     # Do this to guarantee that ghost cells are used; this is a check
c     # on the ghost-fill procedures.  Could raise an exception if face
c     # a patch communicates with more two or more procs.  If this
c     # is uncommented, also uncomment warning in fclaw2d_ghost_fill.cpp
c     compute_slopes = sc

      end
