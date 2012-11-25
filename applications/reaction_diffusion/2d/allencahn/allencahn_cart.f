c     # ------------------------------------------------------
c     # Details of the AllenCahn flow.
c     #
c     # Complete definition of source term.  This is only needed
c     # if diffusion coefficients were not included in the
c     # computation of the laplacian ('apply_lb') and/or you want
c     # to include reaction terms
c     # ------------------------------------------------------
      subroutine src_diffusion(mx,my,mbc,meqn,t,lap)
      implicit none

      integer mx,my,meqn,mbc
      double precision t
      double precision lap(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

c     # Do nothing, since our diffusion coefficient is one

      return
      end

      subroutine src_reaction(mx,my,mbc,meqn,t,q,rhs)
      implicit none

      integer mx,my,meqn,mbc
      double precision t
      double precision q(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision rhs(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j
      double precision u,lap_u

      double precision D
      common /comsrc/ D

      do j = 1,my
         do i = 1,mx
            u = q(i,j,1)

c           # Treat only diffusion term explicitly
            rhs(i,j,1) = (u-u**3)/D**2
         enddo
      enddo

      return
      end

      subroutine compute_rhs_exp(mx,my,mbc,meqn,diff,reac,rhs)
      implicit none

      integer mx,my,meqn, mbc
      double precision diff(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision reac(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision rhs(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,mq,m

      do mq = 1,meqn
         m = mq
         do j = 1,my
            do i = 1,mx
               rhs(i,j,m) = diff(i,j,m) + reac(i,j,m)
            enddo
         enddo
      enddo

      end

      subroutine set_zeros(mx,my,mbc,meqn,rhs)
      implicit none

      integer mx,my,meqn, mbc
      double precision rhs(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)

      integer i,j,m,mq, ibc, jbc


c     # --------------------------------------------
c     # Primary cells
c     # --------------------------------------------
c     # Left and right ghost cells
      do mq = 1,meqn
         m = mq
         do j = 1,my
            do ibc = 1,mbc
               rhs(1-ibc,j,m)  = 0
               rhs(mx+ibc,j,m) = 0
            enddo
         enddo

c        # Bottom and top ghost cells
         do i = 1,mx
            do jbc = 1,mbc
               rhs(i,1-jbc,m)  = 0
               rhs(i,my+jbc,m) = 0
            enddo
         enddo

c        # corners
         do ibc = 1,mbc
            do jbc = 1,mbc
               rhs(1-ibc,1-jbc,m) = 0
               rhs(1-ibc,my+jbc,m) = 0
               rhs(mx+ibc,1-jbc,m) = 0
               rhs(mx+ibc,my+jbc,m) = 0
            enddo
         enddo
      enddo

      end


      subroutine src_reaction_irkc(i,j,mpara,t,q_pt,rhs_pt,
     &      compute_jac,jac)
      implicit none

      integer mx,my,meqn,mpara
      double precision t,q_pt(mpara),rhs_pt(mpara)
      integer i,j
      double precision jac(mpara,mpara)
      logical compute_jac

      double precision u

      double precision D
      common /comsrc/ D

      u = q_pt(1)
      rhs_pt(1) = (u-u**3)/D**2
      if (compute_jac) then
         jac(1,1) = (1.d0-3*u**2)/D**2
      endif

      end

c     # ------------------------------------------------------
c     # ASSIGN_BC
c     #
c     # Set physical boundary conditions for primal cells.  Assume
c     # that boundary conditions are of the form
c     #
c     #            a(x,y) u + b(x,y) u_n = c(x,y,t)
c     #
c     # c(x,y,t) can be updated in 'update_bc", called from RHS
c     # routine to RKC solver.  See below
c     # ------------------------------------------------------
      subroutine assign_bc(mx,my,meqn,xlower,ylower,dx,dy,
     &      a_bc, b_bc, c_bc,ldbc)
      implicit none

      integer mx,my,ldbc,meqn
      double precision xlower,ylower,dx,dy

      double precision a_bc(ldbc,4,meqn)
      double precision b_bc(ldbc,4,meqn)
      double precision c_bc(ldbc,4,meqn)

      double precision a,b,c
      integer i,j,m

c     # Boundary conditions are assigned as follows :
c     #
c     # V_bc(:,1,m) : left edge of equation m
c     # V_bc(:,2,m) : bottom edge of equation m
c     # V_bc(:,3,m) : right edge of equation m
c     # V_bc(:,4,m) : top edge of equation m
c     #
c     # Where 'V' stands for 'a', 'b', or 'c'
c     #


c     # Example : We set no-flux boundary conditions.
      a = 0
      b = 1
      c = 0

      do m = 1,meqn
c        # left and right edges
         do j = 1,my
            a_bc(j,1,m) = a
            a_bc(j,3,m) = a

            b_bc(j,1,m) = b
            b_bc(j,3,m) = b

            c_bc(j,1,m) = c
            c_bc(j,3,m) = c
         enddo

c        # Bottom and top edges
         do i = 1,mx
            a_bc(i,2,m) = a
            a_bc(i,4,m) = a

            b_bc(i,2,m) = b
            b_bc(i,4,m) = b

            c_bc(i,2,m) = c
            c_bc(i,4,m) = c
         enddo
      enddo

      end


c     # -------------------------------------------------------
c     # UPDATE_BC
c     #
c     # Set time  dependent terms in inhomogenous boundary data,
c     #
c     # i.e. the "c" in A*q + B*q_n = C(t)
c     #
c     # If C does not depend on time, this routine can be written
c     # to return without doing anything.
c     # -------------------------------------------------------
      subroutine update_bc(mx,my,meqn,xlower,ylower,dx,dy,t,
     &      c_bc,ldbc)
      implicit none

      integer mx,my, meqn, ldbc
      double precision xlower,ylower,dx,dy,t
      double precision c_bc(ldbc,4,meqn)

      integer i,j,m

      return

      end


c     # -------------------------------------------------------
c     # APPLY_LB
c     #
c     # Apply the diffusion operator. Boundary data has aleady
c     # been supplied to qp and qd.
c     # -------------------------------------------------------
      subroutine apply_lb(mx,my,mbc,meqn,dx,dy,qp,Ap)
      implicit none

      integer mx,my, meqn, mbc

      double precision qp(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision ap(1-mbc:mx+mbc,1-mbc:my+mbc,meqn)
      double precision dx,dy

      integer i,j,m,mq
      double precision fimh, gjmh
      double precision c1, c2

      c1 = 1.d0
      c2 = 1.d0

c     # Loop over each equation
      do mq = 1,meqn
         m = mq
c        # Loop over left edges
         do j = 1,my
            do i = 1,mx+1
c               write(6,*) qp(i,j,m)
               fimh = c1*(qp(i,j,m)   - qp(i-1,j,m))
               Ap(i,j,m)   = -fimh
               Ap(i-1,j,m) = Ap(i-1,j,m) + fimh
            enddo
         enddo

c        % Loop over bottom edges
         do j = 1,my+1
            do i = 1,mx
               gjmh = c2*(qp(i,j,m)   - qp(i,j-1,m))
               Ap(i,j,m)   = Ap(i,j,m) - gjmh
               Ap(i,j-1,m) = Ap(i,j-1,m) + gjmh
            enddo
         enddo

c        # Divide through by area to get an approximation
c        # to the pointwise Laplacian
         do j = 1,my
            do i = 1,mx
               Ap(i,j,m) = Ap(i,j,m)/(dx*dy)
            enddo
         enddo
      enddo

      end

c     # -------------------------------------------------------
c     # COMPUTE_AVG
c     #
c     # Apply the averaging  operator to get dual values
c     # -------------------------------------------------------
      subroutine compute_dual_phys(mx,my,meqn,qp,qd)
      implicit none

      integer mx,my,meqn

      double precision qp(0:mx+1,0:my+1,*)
      double precision qd(mx+1,my+1,*)

      integer i,j,m

c     # Dual interior cells.
      do m = 1,meqn
         do j = 2,my
            do i = 2,mx
               if (i .eq. j) then
                  qd(i,j,m) = 0.5*(qp(i,j,m) + qp(i-1,j-1,m))
               elseif (abs(i+j) .eq. mx+2) then
                  qd(i,j,m) = 0.5d0*(qp(i-1,j,m) + qp(i,j-1,m))
               else
                  qd(i,j,m) = 0.25d0*(qp(i-1,j,m) +
     &                  qp(i-1,j-1,m) + qp(i,j-1,m) + qp(i,j,m))
               endif
            enddo
         enddo
      enddo

      end
