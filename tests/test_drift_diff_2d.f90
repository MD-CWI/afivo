program test_drift_diff
  use m_afivo_2d

  implicit none

  integer, parameter :: dp = kind(0.0d0)
  integer, parameter :: box_size    = 8
  integer, parameter :: i_phi       = 1
  integer, parameter :: i_phi_old   = 2
  integer, parameter :: i_flux      = 1

  type(a2_t)         :: tree
  integer            :: i, n, n_steps, id
  integer            :: ix_list(2, 1)
  integer            :: nb_list(4, 1)
  integer            :: n_boxes_init = 10*1000
  integer            :: n_lvls_max  = 20
  integer            :: n_changes, output_cnt
  real(dp)           :: dr, dt, time, end_time
  real(dp)           :: dt_adapt, dt_output
  real(dp)           :: diff_coeff, vel_x, vel_y, dr_min(2)
  character(len=40)  :: fname

  logical            :: write_out
  integer            :: time_step_method = 2

  dr = 4.0_dp / box_size

  print *, "Initialize tree"
  call a2_init(tree, n_lvls_max, n_boxes_init, box_size, n_var_cell=2, &
       n_var_face=2, dr = dr, r_min = [0.0_dp, 0.0_dp], coarsen_to=-1)

  ! Set up geometry
  id             = 1
  ix_list(:, id) = [1,1] ! Set index of box
  nb_list(:, id) = id    ! Box is periodic, so its own neighbor

  print *, "Create the base mesh"
  call a2_set_base(tree, ix_list, nb_list)

  print *, "Set up the initial conditions"
  do i = 1, 20
     ! We should only set the finest level, but this also works
     call a2_loop_box(tree, set_init_cond)
     call a2_gc_sides(tree, i_phi, a2_sides_interp, have_no_bc)
     call a2_adjust_refinement(tree, ref_func, n_changes)
     if (n_changes == 0) exit
  end do

  print *, "Restrict the initial conditions"
  call a2_restrict_tree(tree, i_phi)
  call a2_gc_sides(tree, i_phi, a2_sides_interp, have_no_bc)

  i          = 0
  output_cnt = 0
  time       = 0
  dt_adapt   = 0.01_dp
  dt_output  = 0.05_dp
  end_time   = 2.0_dp
  diff_coeff = 0.0_dp
  vel_x      = 2.0_dp
  vel_y      = 0.0_dp

  !$omp parallel private(n)
  do
     !$omp single
     i       = i + 1
     dr_min  = a2_min_dr(tree)
     dt      = 0.5_dp / (2 * diff_coeff * sum(1/dr_min**2) + &
          sum( abs([vel_x, vel_y]) / dr_min ) + epsilon(1.0_dp))
     n_steps = ceiling(dt_adapt/dt)
     dt      = dt_adapt / n_steps
     time    = time + dt_adapt

     if (output_cnt * dt_output < time) then
        write_out = .true.
        output_cnt = output_cnt + 1
        write(fname, "(A,I0,A)") "test_drift_diff_2d_", output_cnt, ".vtu"
     else
        write_out = .false.
     end if
     !$omp end single

     if (write_out) call a2_write_tree(tree, trim(fname), &
          (/"phi", "tmp"/), output_cnt, time)

     if (time > end_time) exit

     select case (time_step_method)
     case (1)
        ! Forward Euler
        do n = 1, n_steps
           call a2_loop_boxes_arg(tree, fluxes_koren, [diff_coeff, vel_x, vel_y])
           call a2_consistent_fluxes(tree, [1])
           call a2_loop_box_arg(tree, update_solution, [dt])
           call a2_restrict_tree(tree, i_phi)
           call a2_gc_sides(tree, i_phi, a2_sides_interp, have_no_bc)
        end do
     case (2)
        do n = 1, n_steps
           ! Copy previous solution
           call a2_tree_copy_cc(tree, i_phi, i_phi_old)

           ! Two forward Euler steps over dt
           do i = 1, 2
              call a2_loop_boxes_arg(tree, fluxes_koren, [diff_coeff, vel_x, vel_y])
              call a2_consistent_fluxes(tree, [1])
              call a2_loop_box_arg(tree, update_solution, [dt])
              call a2_restrict_tree(tree, i_phi)
              call a2_gc_sides(tree, i_phi, a2_sides_interp, have_no_bc)
           end do

           ! Take average of phi_old and phi
           call a2_loop_box(tree, average_phi)
        end do
     end select

     call a2_restrict_tree(tree, i_phi)
     call a2_gc_sides(tree, i_phi, a2_sides_interp, have_no_bc)

     call a2_adjust_refinement(tree, ref_func)
     call a2_loop_boxes(tree, prolong_to_new_children)
     call a2_gc_sides(tree, i_phi, a2_sides_interp, have_no_bc)
     call a2_tidy_up(tree, 0.8_dp, 0.5_dp, 10000, .false.)
  end do
  !$omp end parallel

  call a2_destroy(tree)

contains

  integer function ref_func(box)
    type(box2_t), intent(in) :: box
    real(dp)                 :: diff
    integer                  :: i, j, nc

    nc   = box%n_cell
    diff = max( &
         maxval(abs(box%cc(1:nc+1, 1:nc, i_phi) - box%cc(0:nc, 1:nc, i_phi))), &
         maxval(abs(box%cc(1:nc, 1:nc+1, i_phi) - box%cc(1:nc, 0:nc, i_phi))))

    if (box%lvl < 4 .or. diff > 0.15_dp) then
       ref_func = a5_do_ref
    else
       ref_func = a5_rm_ref
    end if
  end function ref_func

  subroutine set_init_cond(box)
    type(box2_t), intent(inout) :: box
    integer                     :: i, j, nc
    real(dp)                    :: xy(2)

    nc = box%n_cell
    do j = 0, nc+1
       do i = 0, nc+1
          xy = a2_r_cc(box, [i,j])
          if (norm2(xy - 2) < 1) then
             box%cc(i, j, i_phi) = 1
          else if (norm2(xy - 2) < 1.1_dp) then
             box%cc(i, j, i_phi) = (1.1_dp - norm2(xy - 2)) * 10
          else
             box%cc(i, j, i_phi) = 0
          end if
       end do
    end do
  end subroutine set_init_cond

  subroutine fluxes_upwind1(boxes, id, flux_args)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp), intent(in)        :: flux_args(:)
    real(dp)                    :: inv_dr
    integer                     :: nc

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    ! Diffusion
    boxes(id)%fx(:,:,i_phi) = boxes(id)%cc(0:nc, 1:nc, i_phi) - boxes(id)%cc(1:nc+1, 1:nc, i_phi)
    boxes(id)%fx(:,:,i_phi) = boxes(id)%fx(:,:,i_phi) * flux_args(1) * inv_dr
    boxes(id)%fy(:,:,i_phi) = boxes(id)%cc(1:nc, 0:nc, 1) - boxes(id)%cc(1:nc, 1:nc+1, 1)
    boxes(id)%fy(:,:,i_phi) = boxes(id)%fy(:,:,i_phi) * flux_args(1) * inv_dr

    ! Drift (1st order upwind, which is very diffusive!)
    if (flux_args(2) > 0) then
       boxes(id)%fx(:,:,i_phi) = boxes(id)%fx(:,:,i_phi) + &
            flux_args(2) * boxes(id)%cc(0:nc, 1:nc, 1)
    else
       boxes(id)%fx(:,:,i_phi) = boxes(id)%fx(:,:,i_phi) + &
            flux_args(2) * boxes(id)%cc(1:nc+1, 1:nc, 1)
    end if
    if (flux_args(3) > 0) then
       boxes(id)%fy(:,:,i_phi) = boxes(id)%fy(:,:,i_phi) + &
            flux_args(3) * boxes(id)%cc(1:nc, 0:nc, 1)
    else
       boxes(id)%fy(:,:,i_phi) = boxes(id)%fy(:,:,i_phi) + &
            flux_args(3) * boxes(id)%cc(1:nc, 1:nc+1, 1)
    end if
  end subroutine fluxes_upwind1

  elemental function limiter_koren(theta)
    real(dp), intent(in) :: theta
    real(dp)             :: limiter_koren
    real(dp), parameter  :: one_sixth = 1.0_dp / 6.0_dp
    limiter_koren = max(0.0d0, min(1.0_dp, theta, (1.0_dp + 2.0_dp * theta) * one_sixth))
  end function limiter_koren

  elemental function ratio(numerator, denominator)
    real(dp), intent(in) :: numerator, denominator
    real(dp)             :: ratio
    if (denominator /= 0.0_dp) then
       ratio = numerator / denominator
    else
       ratio = numerator * huge(1.0_dp)
    end if
  end function ratio

  subroutine fluxes_koren(boxes, id, flux_args)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id
    real(dp), intent(in)        :: flux_args(:)
    real(dp)                    :: inv_dr, theta
    real(dp)                    :: gradp, gradc, gradn
    integer                     :: i, j, nc, nb_id

    nc     = boxes(id)%n_cell
    inv_dr = 1/boxes(id)%dr

    ! x-fluxes interior, advective part with flux limiter
    do j = 1, nc
       do i = 1, nc+1
          gradc = boxes(id)%cc(i, j, i_phi) - boxes(id)%cc(i-1, j, i_phi)
          if (flux_args(2) < 0.0_dp) then

             if (i == nc+1) then
                nb_id = boxes(id)%neighbors(a2_nb_hx)
                if (nb_id > a5_no_box) then
                   gradn = boxes(nb_id)%cc(2, j, i_phi) - boxes(id)%cc(i, j, i_phi)
                else
                   gradn = 0
                end if
             else
                gradn = boxes(id)%cc(i+1, j, i_phi) - boxes(id)%cc(i, j, i_phi)
             end if

             theta = ratio(gradc, gradn)
             boxes(id)%fx(i, j, i_phi) = flux_args(2) * &
                  (boxes(id)%cc(i, j, i_phi) - limiter_koren(theta) * gradn)
          else                  ! flux_args(2) > 0

             if (i == 1) then
                nb_id = boxes(id)%neighbors(a2_nb_lx)
                if (nb_id > a5_no_box) then
                   gradp = boxes(id)%cc(i-1, j, i_phi) - boxes(nb_id)%cc(nc-1, j, i_phi)
                else
                   gradp = 0
                end if
             else
                gradp = boxes(id)%cc(i-1, j, i_phi) - boxes(id)%cc(i-2, j, i_phi)
             end if

             theta = ratio(gradc, gradp)
             boxes(id)%fx(i, j, i_phi) = flux_args(2) * &
                  (boxes(id)%cc(i-1, j, i_phi) + limiter_koren(theta) * gradp)
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
          boxes(id)%fx(i, j, i_phi) = boxes(id)%fx(i, j, i_phi) - &
               flux_args(1) * gradc * inv_dr
       end do
    end do

    ! y-fluxes interior, advective part with flux limiter
    do j = 1, nc+1
       do i = 1, nc
          gradc = boxes(id)%cc(i, j, i_phi) - boxes(id)%cc(i, j-1, i_phi)
          if (flux_args(3) < 0.0_dp) then

             if (j == nc+1) then
                nb_id = boxes(id)%neighbors(a2_nb_hy)
                if (nb_id > a5_no_box) then
                   gradn = boxes(nb_id)%cc(i, 2, i_phi) - boxes(id)%cc(i, j, i_phi)
                else
                   gradn = 0
                end if
             else
                gradn = boxes(id)%cc(i, j+1, i_phi) - boxes(id)%cc(i, j, i_phi)
             end if

             theta = ratio(gradc, gradn)
             boxes(id)%fy(i, j, i_phi) = flux_args(3) * &
                  (boxes(id)%cc(i, j, i_phi) - limiter_koren(theta) * gradn)
          else                  ! flux_args(3) > 0

             if (j == 1) then
                nb_id = boxes(id)%neighbors(a2_nb_ly)
                if (nb_id > a5_no_box) then
                   gradp = boxes(id)%cc(i, j-1, i_phi) - boxes(nb_id)%cc(i, nc-1, i_phi)
                else
                   gradp = 0
                end if
             else
                gradp = boxes(id)%cc(i, j-1, i_phi) - boxes(id)%cc(i, j-2, i_phi)
             end if

             theta = ratio(gradc, gradp)
             boxes(id)%fy(i, j, i_phi) = flux_args(3) * &
                  (boxes(id)%cc(i, j-1, i_phi) + limiter_koren(theta) * gradp)
          end if

          ! Diffusive part with 2-nd order explicit method. dif_f has to be scaled by 1/dx
          boxes(id)%fy(i, j, i_phi) = boxes(id)%fy(i, j, i_phi) - &
               flux_args(1) * gradc * inv_dr
       end do
    end do
  end subroutine fluxes_koren

  subroutine update_solution(box, dt)
    type(box2_t), intent(inout) :: box
    real(dp), intent(in)        :: dt(:)
    real(dp)                    :: inv_dr
    integer                     :: nc, i, j
    real(dp), parameter :: eps = epsilon(1.0_dp)

    nc                    = box%n_cell
    inv_dr                = 1/box%dr
    box%cc(1:nc, 1:nc, i_phi) = box%cc(1:nc, 1:nc, i_phi) + dt(1) * ( &
         (box%fx(1:nc, :, i_phi) - box%fx(2:nc+1, :, i_phi)) * inv_dr + &
         (box%fy(:, 1:nc, i_phi) - box%fy(:, 2:nc+1, i_phi)) * inv_dr)
  end subroutine update_solution

  subroutine average_phi(box)
    type(box2_t), intent(inout) :: box
    box%cc(:, :, i_phi) = 0.5_dp * (box%cc(:, :, i_phi) + box%cc(:, :, i_phi_old))
  end subroutine average_phi

  subroutine prolong_to_new_children(boxes, id)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in) :: id

    if (btest(boxes(id)%tag, a5_bit_new_children)) then
       call a2_prolong1_from(boxes, id, i_phi, .true.)
       boxes(id)%tag = ibclr(boxes(id)%tag, a5_bit_new_children)
    end if
  end subroutine prolong_to_new_children

  subroutine have_no_bc(boxes, id, i, iv)
    type(box2_t), intent(inout) :: boxes(:)
    integer, intent(in)         :: id, i, iv
    stop "We have no boundary conditions in this example"
  end subroutine have_no_bc

end program test_drift_diff
