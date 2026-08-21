#include "../src/cpp_macros.h"
!> \example benchmark_silo.f90
!>
!> This example shows how to create an AMR tree, perform random refinement,
!> and benchmark the writing of Silo output files.
program benchmark_silo
  use m_af_all

  implicit none

  type(af_t)         :: tree
  integer            :: grid_size(NDIM)
  real(dp)           :: domain_size(NDIM)
  logical            :: periodic(NDIM) = .true.
  integer            :: iter, boxes_used
  integer, parameter :: coord_type     = af_xyz ! af_xyz or af_cyl

  integer, parameter :: box_size = 8
  integer, parameter :: i_phi = 1
  type(ref_info_t)   :: ref_info
  character(len=100) :: fname, arg
  integer            :: count_rate, t_start, t_end
  integer            :: max_ref_lvl
  real(dp)           :: t_total, bytes_per_write, bandwidth
  integer            :: n_writes


  write(*,'(A,I0,A)') 'program benchmark_silo_', NDIM, "d"

  print *, "Number of threads", af_get_max_threads()

  if (command_argument_count() /= 1) then
     error stop "Usage: ./benchmark_silo max_refinement_level"
  end if
  call get_command_argument(1, arg)
  read(arg, *) max_ref_lvl

  call af_add_cc_variable(tree, "phi")

  call af_set_cc_methods(tree, 1, af_bc_dirichlet_zero, &
       prolong=af_prolong_linear)

  ! Initialize tree
  grid_size(1) = 2 * box_size
  grid_size(2:) = box_size
  domain_size(1) = 2 * acos(-1.0_dp)
  domain_size(2:) = acos(-1.0_dp)

  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       domain_size, &
       grid_size, &
       periodic=periodic, &
       coord=coord_type)

  call af_print_info(tree)

  ! Set variables on base by using the helper functions af_loop_box(tree, sub)
  call af_loop_box(tree, set_init_cond)

  ! Fill ghost cells for phi.
  call af_gc_tree(tree, [i_phi])

  ! Refine the tree up to the maximum refinement level
  boxes_used = 1
  do
     call af_adjust_refinement(tree, ref_routine, ref_info, ref_buffer=0)
     boxes_used = boxes_used + ref_info%n_add - ref_info%n_rm
     write(*,'(3(3x,A,1x,i6))') "# new boxes", ref_info%n_add, &
          "# boxes used", boxes_used, " highest level", tree%highest_lvl
     if (tree%highest_lvl >= max_ref_lvl) exit
  end do

  call af_print_info(tree)

  ! Benchmark writing Silo files
  call system_clock(t_start, count_rate)
  do iter = 1, 10
     write(fname, "(A,I0)") "output/benchmark_silo_" // DIMNAME // "_", iter
     call af_write_silo(tree, trim(fname), n_cycle=iter)
  end do
  call system_clock(t_end, count_rate)


  n_writes        = iter - 1
  t_total         = (t_end - t_start) / real(count_rate, dp)
  bytes_per_write = af_num_leaves_used(tree) * &
       real((box_size+1)**NDIM * tree%n_var_cell, dp) * 8.0_dp
  bandwidth       = n_writes * bytes_per_write / t_total

  write(*, '(A,i3,1x,A,f8.2,1x,A,/)') &
       ' Wall-clock time after ', n_writes, &
       ' Silo writes: ', t_total, ' seconds'
  write(*, '(A,f10.2,A)') &
       ' Estimated size per file: ', bytes_per_write / (1024.0_dp**2), ' MB'
  write(*, '(A,f10.2,A)') &
       ' Estimated bandwidth: ', bandwidth / (1024.0_dp**2), ' MB/s'

  call af_destroy(tree)

contains

  ! Return the refinement flag for boxes(id)
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box ! A list of all boxes in the tree
    integer, intent(out)     :: cell_flags(DTIMES(box%n_cell))
    real(dp)                 :: rr

    ! Draw a [0, 1) random number
    call random_number(rr)

    if (rr < 0.5_dp**NDIM .and. box%lvl < max_ref_lvl) then
       cell_flags = af_do_ref   ! Add refinement
    else
       cell_flags = af_keep_ref ! Keep refinement
    end if
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box_t), intent(inout) :: box
    integer                     :: IJK, nc
    real(dp)                    :: rr(NDIM)

    nc = box%n_cell
    do KJI_DO(1,nc)
       ! Get the coordinate of the cell center at i,j
       rr = af_r_cc(box, [IJK])

       ! Set the values at each cell according to some function
#if NDIM > 1
       box%cc(IJK, i_phi) = sin(0.5_dp * rr(1))**2 * cos(rr(2))**2
#else
       box%cc(IJK, i_phi) = sin(0.5_dp * rr(1))**2
#endif
    end do; CLOSE_DO
  end subroutine set_init_cond

end program benchmark_silo
