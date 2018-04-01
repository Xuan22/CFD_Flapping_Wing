!   
! fortran_dialect=f90
!
!               Accelerating Reference Frame (ARF) Module
!
! DESC:
!   Provides procedures for solving Navier-Stokes in an accelerating
!   reference frame. In short, the requisitic modifications to a NS
!   solver are the addition of 'source' terms on the moment and energy
!   equations 

!   RECALL: NS momentum equations is nothing but Newton's
!   second law applied to a continuum + the appropriate constitutive
!   relation to close the relationship between the motion of the
!   body and internal stress. Hence, kinematic quantities such
!   as velocity and acceleration must be properly transformed from
!   an inertial reference to the working Accelerating Reference Frame.
!   In addition, the kinetic energy of a fluid particle must also be
!   transformed from an inertial reference to the ARF for consistency
!   between observers in inertial and accelerating frames. This is the
!   basis for the foregoing momentum and energy 'source' terms. Of course,
!   if the motion of the reference frame consists of a steady translation
!   without rotation, then the source terms vanish to zero.
!
! DEVELOPER(S):
!   Benjamin Silbaugh
!
! LAST MODIFIED:
!   March 12, 2008
!
!==========================================================================

module arf_mod

  use params_global

  implicit none

  !Solver flag
  integer :: arf_opt = 0

  !Frame translational velocity and acceleration:
  real, dimension(3) :: frame_vel, frame_acc

  !Frame angular velocity and acceleration:
  real, dimension(3) :: frame_avel

contains

subroutine init_arf_mod()
!DESC: Module initialization procedure
  frame_vel = 0.0
  frame_acc = 0.0
  frame_avel = 0.0
  print*, "ARF module initialized"
end subroutine init_arf_mod

subroutine add_arf_source(q, s)
!DESC: Computes source terms due to frame acceleration and adds
!      them to the right-hand-side vector 's'.
!
! NOTE: Assumes frame_vel, frame_acc, etc have been properly
!       non-dimensionalized in a manner consistent with conservative
!       flow quantities in q. Also, assuming that q has already been scaled
!       by 1/J.
  real, dimension(jmax,kmax,lmax,nd), intent(in) :: q
  real, dimension(jmax,kmax,lmax,nv), intent(inout) :: s
  !Local var:
  integer :: j, k, l
  real :: rut, rvt, rwt, rudot, rvdot, rwdot, a1, a2, a3

  !Components of total acceleration:
   a1 = frame_acc(1) &
      + frame_avel(2)*frame_vel(3) - frame_avel(3)*frame_vel(2)
   a2 = frame_acc(2) &
      + frame_avel(3)*frame_vel(1) - frame_avel(1)*frame_vel(3)
   a3 = frame_acc(3) &
      + frame_avel(1)*frame_vel(2) - frame_avel(2)*frame_vel(1)

  do l = 2, lmax - 1
    do k = 2, kmax - 1
      do j = 2, jmax - 1

        ! density times frame velocity + flow velocity components:
        rut = q(j,k,l,1)*frame_vel(1) + q(j,k,l,2)
        rvt = q(j,k,l,1)*frame_vel(2) + q(j,k,l,3)
        rwt = q(j,k,l,1)*frame_vel(3) + q(j,k,l,4)

        ! density times frame partial-acceleration components:
        rudot = q(j,k,l,1)*frame_acc(1)
        rvdot = q(j,k,l,1)*frame_acc(2)
        rwdot = q(j,k,l,1)*frame_acc(3)

        ! Add momentum source:
        s(j,k,l,2) = s(j,k,l,2) &
                   - (rudot + frame_avel(2)*rwt - frame_avel(3)*rvt)
        s(j,k,l,3) = s(j,k,l,3) &
                   - (rvdot - frame_avel(1)*rwt + frame_avel(3)*rut)
        s(j,k,l,4) = s(j,k,l,4) &
                   - (rwdot + frame_avel(1)*rvt - frame_avel(2)*rut)

        ! Add energy source:
        s(j,k,l,5) = s(j,k,l,5) &
                   - (q(j,k,l,2)*a1 + q(j,k,l,3)*a2 + q(j,k,l,4)*a3)

      end do
    end do
  end do
end subroutine add_arf_source

subroutine arf_momentum_source(qjkl, sx, sy, sz)
!DESC: Computes momentum source cartesian components; used in
!      wall bc routine. Assumes qijk to be the vector of conserved
!      variables for the cell (j,k,l).
!
  real, dimension(nd), intent(in) :: qjkl
  real, intent(out) :: sx, sy, sz
  !local var
  real :: rut, rvt, rwt, rudot, rvdot, rwdot

  ! density times frame velocity + flow velocity components:
  rut = qjkl(1)*frame_vel(1) + qjkl(2)
  rvt = qjkl(1)*frame_vel(2) + qjkl(3)
  rwt = qjkl(1)*frame_vel(3) + qjkl(4)

  ! density times frame partial-acceleration components:
  rudot = qjkl(1)*frame_acc(1)
  rvdot = qjkl(1)*frame_acc(2)
  rwdot = qjkl(1)*frame_acc(3)

  ! momentum source:
  sx = -(rudot + frame_avel(2)*rwt - frame_avel(3)*rvt)
  sy = -(rvdot - frame_avel(1)*rwt + frame_avel(3)*rut)
  sz = -(rwdot + frame_avel(1)*rvt - frame_avel(2)*rut)

end subroutine arf_momentum_source

subroutine add_arf_gridvel(x, y, z, ug, vg, wg)
!DESC: Adds the additional grid velocity terms due to unsteady frame
!      motion. NOTE: This "adds" the additional velocity terms - be
!      sure to set the "usual" grid velocities before calling this
!      routine.
!
!      *** THIS IS AN ALTERNATIVE TO THE SOURCE TERMS ABOVE ***
!                        DO NOT USE BOTH
!
!
  real, dimension(jmax,kmax,lmax), intent(in) :: x, y, z
  real, dimension(jmax,kmax,lmax), intent(inout) :: ug, vg, wg
  !local var
  integer :: j, k, l

  do l = 1, lmax
    do k = 1, kmax
      do j = 1, jmax
        ug(j,k,l) = ug(j,k,l) + frame_vel(1) &
                  + frame_avel(2)*z(j,k,l) - frame_avel(3)*y(j,k,l)
        vg(j,k,l) = vg(j,k,l) + frame_vel(2) &
                  - frame_avel(1)*z(j,k,l) + frame_avel(3)*x(j,k,l)
        wg(j,k,l) = wg(j,k,l) + frame_vel(3) &
                  + frame_avel(1)*y(j,k,l) - frame_avel(2)*x(j,k,l)
      end do
    end do
  end do

end subroutine add_arf_gridvel

end module arf_mod
