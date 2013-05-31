       subroutine rkc_stats(tip1,dt)
       implicit none

       double precision tip1,dt

       integer nfe,nsteps,naccpt,nrejct,nfesig,maxm,stats(6)

c     # Common block for RKC solver - this is set by solver.
       common /rkcdid/   nfe,nsteps,naccpt,nrejct,nfesig,maxm

       write(6,101) 'Time ', tip1
       write(6,102) 'dt ',dt
       write(6,100) 'Number of f evaluations', nfe
       write(6,100) 'Number of steps', nsteps
       write(6,100) 'Number of accepted steps', naccpt
       write(6,100) 'Number of rejected steps', nrejct
       write(6,100) 'Number of evals. for Jacobian', nfesig
       write(6,100) 'Maximum number of stages used', maxm
       write(6,100) 'Number of fevals (curr. step)', nfesig + nfe
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
