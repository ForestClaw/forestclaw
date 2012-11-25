c     # -----------------------------------------------------
c     # parabolic_step_rkc
c     #
c     # Advance the parabolic step from time 't' to time 't + dt'.
c     #
c     #      nv(mq) = .true.  ==> diffuse species mq
c     #      nv(mq) = .false. ==> do not diffuse species mq
c     # mq = 1,..,meqn
c     #
c     # Normally, you have 'nv(mq) == plist(mq)'.  The only case
c     # in which this is not a good choice is the situation in which
c     # the sources terms containing the diffusion terms are
c     # completely decoupled.  In this case, it is much faster to
c     # solve each diffusion term separately.
c     #
c     # After each call, Check ierror to make sure that the RKC
c     # solver returned a valid solution.
c     #
c     #  ierror == 0   Normal exit
c     #
c     #  ierror  > 0   Problem with RKC solver.   Usually this
c     #                means that your time step is too large
c     #                for the requested accuracy.  See the
c     #                description of 'idid' in the RKC documentation.
c     #
c     # This routine is called either from the 'src2' routine in
c     # clawpack, or directly from 'parabolic_main' (if the clawpack
c     # driver is not being used.
c     # -----------------------------------------------------


       subroutine rkc_stats(tip1,dt)
       implicit none

       double precision tip1,dt

       integer nfe,nsteps,naccpt,nrejct,nfesig,maxm,stats(6)

       integer ts_counter
       integer get_time_step_counter

       integer lc_total, lc_inc
       integer get_laplacian_calls_total
       integer get_laplacian_calls_elapsed

c     # Common block for RKC solver - this is set by solver.
       common /rkcdid/   nfe,nsteps,naccpt,nrejct,nfesig,maxm

c       ts_counter = get_time_step_counter()
c       lc_total = get_laplacian_calls_total()
c       lc_inc = get_laplacian_calls_elapsed()

c       call write_line(6,47)
c       write(6,100) 'Step ', ts_counter
       write(6,101) 'Time ', tip1
       write(6,102) 'dt ',dt
c       write(6,100) 'Number of Lap. calls (inc)', lc_inc
c       write(6,100) 'Number of Lap. calls (total)', lc_total
       write(6,100) 'Number of f evaluations', nfe
       write(6,100) 'Number of steps', nsteps
       write(6,100) 'Number of accepted steps', naccpt
       write(6,100) 'Number of rejected steps', nrejct
       write(6,100) 'Number of evals. for Jacobian', nfesig
       write(6,100) 'Maximum number of stages used', maxm
       write(6,100) 'Number of fevals (curr. step)', nfesig + nfe
c       write(6,100) 'Total number of fevals', lc_total

c       call write_line(6,47)
       write(6,*) ' '


  100  format(A35,I12)
  101  format(A35,F12.6)
  102  format(A35,E12.4)

       end

      integer function get_stages_rkc()
      implicit none

      integer nfe,nsteps,naccpt,nrejct,nfesig,maxm,stats(6)

c     # Common block for RKC solver - this is set by solver.
      common /rkcdid/   nfe,nsteps,naccpt,nrejct,nfesig,maxm

      get_stages_rkc = maxm

      end

      integer function get_fevals_rkc()
      implicit none

      integer nfe,nsteps,naccpt,nrejct,nfesig,maxm,stats(6)

c     # Common block for RKC solver - these are set by solver.
      common /rkcdid/   nfe,nsteps,naccpt,nrejct,nfesig,maxm

      get_fevals_rkc = nfesig + nfe

      end

      double precision function spcrad(neqn,t,y)
      implicit none
      integer          neqn
      double precision t,y(neqn)

      spcrad = 0

      return
      end
