!> \example random_refinement_2d.f90

! This example shows how to create an AMR tree, perform random refinement, write
! output files, and how to fill ghost cells.
program random_refinement_2d
  use m_a2_types
  use m_a2_core
  use m_a2_utils
  use m_a2_ghostcell
  use m_a2_output

  implicit none

  type(a2_t)           :: tree
  integer              :: id, iter, boxes_used
  integer, parameter   :: n_boxes_base = 5
  integer              :: ix_list(2, n_boxes_base)
  integer              :: nb_list(4, n_boxes_base)
  integer, parameter   :: box_size     = 8
  integer, parameter   :: i_phi        = 1
  type(ref_info_t)     :: ref_info
  real(dp)             :: dr
  character(len=40)    :: fname
  integer              :: count_rate,t_start,t_end

  write(*,'(A)') 'program random_refinement_2d'

  print *, "Number of threads", af_get_max_threads()

  ! The cell spacing at the coarsest grid level
  dr = 2 * acos(-1.0_dp) / box_size ! 2 * pi / box_size

  ! Initialize tree
  call a2_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       1, &            ! Number of cell-centered variables
       0, &            ! Number of face-centered variables
       dr, &           ! Distance between cells on base level
       cc_names = ["phi"])      ! Optional: names of cell-centered variables

  ! Set up geometry. 
  ! Neighbors for the boxes are stored in nb_list
  nb_list(:, :) = af_no_box     ! Default value

  ! Spatial indices are used to define the coordinates of a box. By default the
  ! box at [1,1] touches the origin (x,y) = (0,0)
  do id = 1, n_boxes_base
     ix_list(:, id)               = [id, 1] ! Boxes at (1,1), (2,1), (3,1) etc.
     nb_list(a2_neighb_highy, id) = id      ! Periodic in y-direction
  end do

  ! Periodic boundaries only have to be specified from one side. The lower-x
  ! neighbor of box 1 is the last box
  nb_list(a2_neighb_lowx, 1)  = n_boxes_base

  ! Create the base mesh, using the box indices and their neighbor information
  call a2_set_base(tree, ix_list, nb_list)
  call a2_print_info(tree)

  ! Set variables on base by using the helper functions a2_loop_box(tree, sub)
  ! and a2_loop_boxes(tree, sub). These functions call the subroutine sub for
  ! each box in the tree, with a slightly different syntax.
  call a2_loop_box(tree, set_init_cond)

  ! Fill ghost cells for phi. The third argument is a subroutine that fills
  ! ghost cells near refinement boundaries, and the fourth argument fill ghost
  ! cells near physical boundaries.
  call a2_gc_tree(tree, i_phi, a2_gc_interp, a2_bc_dirichlet_zero)

  call system_clock(t_start, count_rate)
  boxes_used = 1
  do iter = 1, 50
     ! This writes a VTK output file containing the cell-centered values of the
     ! leaves of the tree (the boxes not covered by refinement).
     ! Variables are the names given as the third argument.
     write(fname, "(A,I0)") "random_refinement_2d_", iter
     call a2_write_vtk(tree, trim(fname), dir="output")

     ! This updates the refinement of the tree, by at most one level per call.
     ! The second argument is a subroutine that is called for each box that can
     ! be refined or derefined, and it should set refinement flags. Information
     ! about the changes in refinement are returned in the third argument.
     !
     ! Within each (de)refinement step the subroutine a2_tidy_up is called, with
     !    max_hole_frac = 0.5_dp
     !    n_clean_min    = 1000
     ! This means that every now and then whether too much holes appear in the
     ! tree box list. Therefore reordering of tree is not necessary at the end
     ! of the iteration process
     call a2_adjust_refinement(tree, ref_routine, ref_info)
     boxes_used = boxes_used + ref_info%n_add - ref_info%n_rm
     write(*,'(4(3x,A,1x,i6))') "# new     boxes", ref_info%n_add, &
                                "# removed boxes", ref_info%n_rm,  &
                                "# boxes used   ", boxes_used,     &
                                " highest level ", tree%highest_lvl

     ! Newly added boxes should be initialized, which is done in this routine.
     call prolong_to_new_children(tree, ref_info)
  end do
  call system_clock(t_end, count_rate)

  write(*, '(A,i3,1x,A,f8.2,1x,A,/)') &
           ' Wall-clock time after ',iter, &
           ' iterations: ', (t_end-t_start) / real(count_rate, dp), &
           ' seconds'

  call a2_print_info(tree)

  ! This call is not really necessary here, but cleaning up the data in a tree
  ! is important if your program continues with other tasks.
  call a2_destroy(tree)

contains

  ! Return the refinement flag for boxes(id)
  subroutine ref_routine(box, cell_flags)
    type(box2_t), intent(in) :: box ! A list of all boxes in the tree
    integer, intent(out)     :: cell_flags(box%n_cell, box%n_cell)
    real(dp)                 :: rr

    ! Draw a [0, 1) random number
    call random_number(rr)

    if (rr < 0.25_dp .and. box%lvl < 7) then
       cell_flags = af_do_ref   ! Add refinement
    else
       cell_flags = af_rm_ref   ! Ask to remove this box, which will not always
                              ! happen (see documentation)
    end if
  end subroutine ref_routine

  ! This routine sets the initial conditions for each box
  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 1, nc
       do i = 1, nc
          ! Get the coordinate of the cell center at i,j
          xy = a2_r_cc(box, [i,j])

          ! Set the values at each cell according to some function
          box%cc(i, j, i_phi) = sin(0.5_dp * xy(1)) * cos(xy(2))
       end do
    end do
  end subroutine set_init_cond

  ! Set values on newly added boxes
  subroutine prolong_to_new_children(tree, ref_info)
    use m_a2_prolong
    type(a2_t), intent(inout)    :: tree
    type(ref_info_t), intent(in) :: ref_info
    integer                      :: lvl, i, id, p_id

    do lvl = 1, tree%highest_lvl
       ! For each newly added box ...
       do i = 1, size(ref_info%lvls(lvl)%add)
          ! Use prolongation to set the values of variable i_phi
          id = ref_info%lvls(lvl)%add(i)
          p_id = tree%boxes(id)%parent

          call a2_prolong_quadratic(tree%boxes(p_id), tree%boxes(id), i_phi)
       end do

       do i = 1, size(ref_info%lvls(lvl)%add)
          id = ref_info%lvls(lvl)%add(i)
          ! After values have been set on this level, fill ghost cells
          call a2_gc_box(tree%boxes, id, i_phi, &
               a2_gc_interp, a2_bc_dirichlet_zero)
       end do
    end do
  end subroutine prolong_to_new_children

end program random_refinement_2d
