#include "../src/cpp_macros.h"
! Test the different implementations of interpolation
program test_interpolation
  use m_af_all

  implicit none

  type(af_t)       :: tree
  integer :: nc, counter, dummy
  real(dp)  :: p(NDIM), tril(NDIM), ord1(NDIM)
  logical :: success
  real(dp) :: t0, t1

  print *, "== Unit vectors =="
  print *, af_unit_vector(:, 1)
  print *, af_unit_vector(:, 2)
#if NDIM == 3
  print *, af_unit_vector(:, 3)
#endif

  call af_add_cc_variable(tree, "phi")
  call af_add_fc_variable(tree, "fc_val")
  call af_set_cc_methods(tree, 1, af_bc_neumann_zero)
  call af_init(tree, 8, [DTIMES(8.0_dp)], [DTIMES(8)])
  call af_loop_box(tree, init_grad)

  print *, "== box_stats =="
  print *, "dr = ", tree%boxes(1)%dr
  print *, "r_min = ", tree%boxes(1)%r_min
  print *, "n_cell = ", tree%boxes(1)%n_cell

  print *, "== corners low =="
  print *, af_get_fc_corner_low(tree%boxes(1), 1)
  print *, af_get_fc_corner_low(tree%boxes(1), 2)
#if NDIM == 3
  print *, af_get_fc_corner_low(tree%boxes(1), 3)
#endif

  print *, "== corners high =="
  print *, af_get_fc_corner_hi(tree%boxes(1), 1)
  print *, af_get_fc_corner_hi(tree%boxes(1), 2)
#if NDIM == 3
  print *, af_get_fc_corner_hi(tree%boxes(1), 3)
#endif

  nc = tree%boxes(1)%n_cell

  ! print *, "== fc values =="
  ! print *, tree%boxes(1)%fc(DTIMES(0:nc+2), 2, 1)

  print *, "== interpolation errors =="
  !Now check the relative difference between interpolations
  do counter = 2, 100
    p = 8.0_dp / counter
    call CPU_TIME(t0)
    do dummy = 1, 10000
      tril = af_interp_trilinear_fc(tree, p, 1, success)
    end do
    call CPU_TIME(t1)
    print *, "CPU time tril: ", t1-t0

    call CPU_TIME(t0)
    do dummy = 1, 10000
      ord1 = af_interp1_fc(tree, p, 1, success)
    end do
    call CPU_TIME(t1)
    print *, "CPU time ord1: ", t1-t0
    ! print *, tril
    ! print *, ord1
    print *, (tril-ord1)**2
    print *, "=="
  end do

contains

  subroutine init_constant(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc

    nc                    = box%n_cell
    box%cc(DTIMES(1:nc), 1) = 0
    box%fc(DTIMES(0:nc+2), :, 1) = 1.5_dp
  end subroutine init_constant

  subroutine init_grad(box)
    type(box_t), intent(inout) :: box
    integer                    :: nc, ii, jj, kk

    nc                    = box%n_cell
    box%cc(DTIMES(1:nc), 1) = 0

#if NDIM == 2
    do jj = 0, nc+2
      do ii = 0, nc+2
        box%fc(ii, jj, :, 1) = ii + jj
      end do
    end do
#elif NDIM == 3
    do kk = 0, nc+2
      do jj = 0, nc+2
        do ii = 0, nc+2
          box%fc(ii, jj, kk, :, 1) = ii+ jj + kk
        end do
      end do
    end do
#endif
  end subroutine init_grad

end program
