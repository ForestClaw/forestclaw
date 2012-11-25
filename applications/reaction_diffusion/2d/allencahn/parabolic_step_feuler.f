c     # -----------------------------------------------------
c     # parabolic_step_feuler
c     #
c     # Advance the parabolic step from time 't' to time 't + dt'
c     # using Forward Euler
c     #
c     #      nv(mq) = .true.  ==> diffuse species mq
c     #      nv(mq) = .false. ==> do not diffuse species mq
c     #                           mq = 1,..,meqn
c     #
c     # Normally, you have 'nv(mq) == plist(mq)'.  The only case
c     # in which this is not a good choice is the situation in which
c     # the sources terms containing the diffusion terms are
c     # completely decoupled.  In this case, it is much faster to
c     # solve each diffusion term separately.
c     #
c     # The variable 'ierror' will always return 0 (i.e. successful return)
c     #
c     # This routine is called either from the 'src2' routine in
c     # clawpack, or directly from 'parabolic_main' (if the clawpack
c     # driver is not being used.
c     # -----------------------------------------------------


      subroutine feuler(f_rhs,neqn,q_amr,t,dt,rhs_work)
      implicit none

      external f_rhs
      integer neqn
      double precision t,dt, dx
      double precision q_amr(neqn)
      double precision rhs_work(neqn)

      double precision t_inner, frac_dt,dt_exp, dt_stab
      double precision dt_inner, tip1
      integer m, dim

      t_inner = t
      tip1 = t + dt

      call f_rhs(neqn,t_inner,q_amr,rhs_work)

c     # This updates the ghost cells as well, but oh well.
      do m = 1,neqn
         q_amr(m) = q_amr(m) + dt*rhs_work(m)
      enddo

      end

c       subroutine feuler_stats(tip1,dt)
c       implicit none
c
c       double precision tip1, dt
c
c       integer ts_counter
c       integer get_time_step_counter
c
c       integer lc_total, lc_inc
c       integer get_laplacian_calls_elapsed
c       integer get_laplacian_calls_total
c
c       integer verbose
c       integer get_verbosity
c
c       lc_inc = get_laplacian_calls_elapsed()
c       lc_total = get_laplacian_calls_total()
c       verbose  = get_verbosity()
c       ts_counter = get_time_step_counter()
c
c       write(6,100) 'Step ', ts_counter
c       write(6,101) 'Time ', tip1
c       write(6,102) 'dt ',dt
c       write(6,100) 'Laplacian calls (this step) ', lc_inc
c       write(6,100) 'Laplacian calls (total) ', lc_total
c       write(6,100) 'Number of internal time steps', lc_inc
c       write(6,100)
c
c   100 format(A35,I12)
c   101 format(A35,F12.6)
c   102 format(A35,E12.4)
c
c       end
