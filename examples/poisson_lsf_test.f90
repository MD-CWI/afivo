#include "../src/cpp_macros.h"
!> \example poisson_lsf_test.f90
!>
!> Test of the level-set functionality of the Poisson solver
program poisson_lsf_test
  use m_af_all

  implicit none

  integer            :: box_size         = 8
  integer, parameter :: n_iterations     = 10
  integer            :: max_refine_level = 2
  integer            :: i_phi
  integer            :: i_rhs
  integer            :: i_tmp
  integer            :: i_lsf
  integer            :: i_error
  integer            :: i_field
  integer            :: i_field_norm
  integer            :: i_normgradlsf
  integer            :: i_lsfmask

  ! Which shape to use, 1 = circle, 2 = heart, 3 = rhombus, 4-5 = triangle
  integer             :: shape             = 1

  ! Sharpness_triangle, can choose a larger number for more acute angle
  integer             :: sharpness_t       = 10
  real(dp), parameter :: boundary_value    = 1.0_dp
  real(dp), parameter :: solution_coeff    = 1.0_dp
  real(dp), parameter :: solution_radius   = 0.05_dp
  real(dp)            :: solution_r0(NDIM) = 0.5_dp ! .51625_dp is a problematic one

  type(af_t)         :: tree
  type(ref_info_t)   :: ref_info
  integer            :: n, mg_iter, coord, n_args
  real(dp)           :: residu, max_error, max_field
  character(len=100) :: fname, argv
  type(mg_t)         :: mg

  call af_add_cc_variable(tree, "phi", ix=i_phi)
  call af_add_cc_variable(tree, "rhs", ix=i_rhs)
  call af_add_cc_variable(tree, "tmp", ix=i_tmp)
  call af_add_cc_variable(tree, "lsf", ix=i_lsf)
  call af_add_cc_variable(tree, "error", ix=i_error)
  call af_add_cc_variable(tree, "field_norm", ix=i_field_norm)
  call af_add_fc_variable(tree, "field", ix=i_field)
  call af_add_cc_variable(tree, "normgradlsf", ix=i_normgradlsf)
  call af_add_cc_variable(tree, "lsfmask", ix=i_lsfmask)

  call af_set_cc_methods(tree, i_lsf, funcval=set_lsf)
  call af_set_cc_methods(tree, i_field_norm, af_bc_neumann_zero)

  ! If an argument is given, switch to cylindrical coordinates in 2D
  n_args = command_argument_count()
  if (n_args > 6) then
     stop "Usage: ./poisson_lsf_test [cyl] [max_refinement_level] [shape_type] "&
          "[size of box] [sharpness_triangle]"
  end if

  coord = af_xyz
  do n = 1, n_args
     call get_command_argument(n, argv)
     if (argv == 'cyl') then
        coord = af_cyl
        ! Place solution on axis
        solution_r0(1) = 0.0_dp
     else if (n == 2) then
        read(argv, *) max_refine_level
     else if (n == 3) then
        read(argv, *) shape
     else if (n == 4) then
        read(argv, *) box_size
     else if (n == 5) then
        read(argv, *) sharpness_t
     else if (n == 6) then
        read(argv, *) solution_r0(NDIM)
     end if
  end do

  mg%i_phi    = i_phi
  mg%i_rhs    = i_rhs
  mg%i_tmp    = i_tmp
  mg%i_lsf    = i_lsf
  mg%sides_bc => bc_solution
  mg%lsf_boundary_value = boundary_value
  mg%lsf_dist => golden_lsf
  mg%lsf => get_lsf

  ! Initialize tree
  call af_init(tree, & ! Tree to initialize
       box_size, &     ! A box contains box_size**DIM cells
       [DTIMES(1.0_dp)], &
       4 * [DTIMES(box_size)], &
       coord=coord)

  do n = 1, 100
     call af_adjust_refinement(tree, ref_routine, ref_info)
     if (ref_info%n_add == 0) exit
  end do

  call af_print_info(tree)

  call mg_init(tree, mg)

  do mg_iter = 1, n_iterations
     call mg_fas_fmg(tree, mg, .true., mg_iter>1)
     call mg_compute_phi_gradient(tree, mg, i_field, 1.0_dp, i_field_norm)
     call af_gc_tree(tree, [i_field_norm])
     call af_loop_box(tree, set_error)

     ! Determine the minimum and maximum residual and error
     call af_tree_maxabs_cc(tree, i_tmp, residu)
     call af_tree_maxabs_cc(tree, i_error, max_error)
     call af_tree_maxabs_cc(tree, i_field_norm, max_field)
     write(*, "(I8,3E14.5)") mg_iter, residu, max_error, max_field
     write(fname, "(A,I0)") "output/poisson_lsf_test_" // DIMNAME // "_", mg_iter
     call af_write_silo(tree, trim(fname))
  end do

  call af_stencil_print_info(tree)

contains

  ! Set refinement flags for box
  subroutine ref_routine(box, cell_flags)
    type(box_t), intent(in) :: box
    integer, intent(out)    :: cell_flags(DTIMES(box%n_cell))
    integer                 :: nc
    integer, parameter      :: refinement_type = 1

    nc = box%n_cell

    select case (refinement_type)
    case (1)
       ! Uniform refinement
       if (box%lvl < max_refine_level) then
          cell_flags(DTIMES(:)) = af_do_ref
       else
          cell_flags(DTIMES(:)) = af_keep_ref
       end if

    case (2)
       ! Only refine at boundary
       if (box%lvl < max_refine_level .and. &
            minval(box%cc(DTIMES(:), i_lsf)) * &
            maxval(box%cc(DTIMES(:), i_lsf)) <= 0) then
          cell_flags(DTIMES(:)) = af_do_ref
       else
          cell_flags(DTIMES(:)) = af_keep_ref
       end if

    case (3)
       ! 'Bad' refinement to test the method
       if (norm2(box%r_min - solution_r0) < solution_radius .and. &
            box%lvl < max_refine_level) then
          cell_flags(DTIMES(:)) = af_do_ref
       else
          cell_flags(DTIMES(:)) = af_keep_ref
       end if
    end select
  end subroutine ref_routine

  real(dp) function solution(r)
    real(dp), intent(in) :: r(NDIM)
    real(dp) :: distance, lsf

    select case (shape)
    case (1)
       ! Relative distance
       distance = norm2(r-solution_r0) / solution_radius

       ! Let values increase for distance -> infinity
       if (distance < 1.0_dp) then
          solution = boundary_value
       else if (NDIM == 1) then
          solution = boundary_value + solution_coeff * (distance - 1.0_dp)
       else if (NDIM == 2 .and. coord == af_xyz) then
          solution = boundary_value + solution_coeff * log(distance)
       else
          solution = boundary_value + solution_coeff * (1 - 1/distance)
       end if
    case (4, 5)
       ! Triangle
       lsf = get_lsf(r)
       if (lsf <= 0.0_dp) then
          solution = boundary_value
       else
          solution = boundary_value * exp(-lsf)
       end if
    case default
       solution = 0.0_dp
    end select
  end function solution

  ! This routine sets the level set function
  subroutine set_lsf(box, iv)
    type(box_t), intent(inout) :: box
    integer, intent(in)        :: iv
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM), norm_dr

    nc = box%n_cell
    norm_dr = norm2(box%dr)

    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, iv) = get_lsf(rr)
       box%cc(IJK, i_normgradlsf) = numerical_gradient_amplitude(get_lsf, rr)

       if (abs(box%cc(IJK, iv)) < norm_dr * box%cc(IJK, i_normgradlsf) * &
            mg%lsf_gradient_safety_factor) then
          box%cc(IJK, i_lsfmask) = 1.0_dp
       else
          box%cc(IJK, i_lsfmask) = 0.0_dp
       end if
    end do; CLOSE_DO
  end subroutine set_lsf

  real(dp) function get_lsf(rr)
    real(dp), intent(in) :: rr(NDIM)
    real(dp)             :: distance
#if NDIM > 1
    real(dp)             :: qq(NDIM), dist1, dist2
#endif

    select case (shape)
    case (1)
       distance = norm2(rr-solution_r0) / solution_radius
       get_lsf = distance - 1.0_dp
#if NDIM > 1
    case (2)
       ! Heart centered on r0
       qq = (rr - solution_r0) * 4.0_dp
       get_lsf = (qq(1)**2 + qq(2)**2 - 1)**3 - &
            qq(1)**2 * qq(2)**3
    case (3)
       ! Rhombus or astroid
       qq = (rr-solution_r0)*4.0_dp
       get_lsf = ((qq(1)**2)**(1.0_dp/3)/0.8) + &
            ((qq(2)**2)**(1.0_dp/3)/1.5) - 0.8_dp
    case (4)
       ! sharpness_t -> for sharpness of the triangle top angle,
       ! larger sharpness_t equals more acute angle
       qq = rr-solution_r0
       get_lsf = sharpness_t * abs(qq(1)) + qq(2)
    case (5)
       ! Triangle v2, uses signed distance from the triangle
       dist1 = GM_dist_line(rr, [solution_r0(1) - solution_r0(2)/sharpness_t, 0.0_dp], &
            solution_r0, 2)
       dist2 = GM_dist_line(rr, [solution_r0(1) + solution_r0(2)/sharpness_t, 0.0_dp], &
            solution_r0, 2)

       ! Determine sign of lsf function
       qq = rr - solution_r0
       get_lsf = sharpness_t * abs(qq(1)) + qq(2)

       ! Use sign in front of minimum distance
       get_lsf = sign(min(dist1, dist2), get_lsf)
#endif
    case default
       error stop "Invalid case"
    end select

  end function get_lsf

  subroutine set_error(box)
    type(box_t), intent(inout) :: box
    integer                    :: IJK, nc
    real(dp)                   :: rr(NDIM)

    nc = box%n_cell
    do KJI_DO(0,nc+1)
       rr = af_r_cc(box, [IJK])
       box%cc(IJK, i_error) = box%cc(IJK, i_phi) - solution(rr)
    end do; CLOSE_DO
  end subroutine set_error

  ! This routine sets boundary conditions for a box
  subroutine bc_solution(box, nb, iv, coords, bc_val, bc_type)
    type(box_t), intent(in) :: box
    integer, intent(in)     :: nb
    integer, intent(in)     :: iv
    real(dp), intent(in)    :: coords(NDIM, box%n_cell**(NDIM-1))
    real(dp), intent(out)   :: bc_val(box%n_cell**(NDIM-1))
    integer, intent(out)    :: bc_type
    integer                 :: n

    bc_type = af_bc_dirichlet

    do n = 1, box%n_cell**(NDIM-1)
       bc_val(n) = solution(coords(:, n))
    end do
  end subroutine bc_solution

  !> Compute distance vector between point and its projection onto a line
  !> between r0 and r1
  subroutine GM_dist_vec_line(r, r0, r1, n_dim, dist_vec, frac)
    integer, intent(in)   :: n_dim
    real(dp), intent(in)  :: r(n_dim), r0(n_dim), r1(n_dim)
    real(dp), intent(out) :: dist_vec(n_dim)
    real(dp), intent(out) :: frac !< Fraction [0,1] along line
    real(dp)              :: line_len2

    line_len2 = sum((r1 - r0)**2)
    frac = sum((r - r0) * (r1 - r0))

    if (frac <= 0.0_dp) then
       frac = 0.0_dp
       dist_vec = r - r0
    else if (frac >= line_len2) then
       frac = 1.0_dp
       dist_vec = r - r1
    else
       dist_vec = r - (r0 + frac/line_len2 * (r1 - r0))
       frac = sqrt(frac / line_len2)
    end if
  end subroutine GM_dist_vec_line

  function GM_dist_line(r, r0, r1, n_dim) result(dist)
    integer, intent(in)  :: n_dim
    real(dp), intent(in) :: r(n_dim), r0(n_dim), r1(n_dim)
    real(dp)             :: dist, dist_vec(n_dim), frac
    call GM_dist_vec_line(r, r0, r1, n_dim, dist_vec, frac)
    dist = norm2(dist_vec)
  end function GM_dist_line

  function numerical_gradient_amplitude(f, r) result(normgrad)
    procedure(mg_func_lsf) :: f
    real(dp), intent(in)   :: r(NDIM)
    real(dp), parameter    :: sqrteps      = sqrt(epsilon(1.0_dp))
    real(dp), parameter    :: min_stepsize = epsilon(1.0_dp)
    real(dp)               :: r_eval(NDIM), gradient(NDIM)
    real(dp)               :: stepsize(NDIM), flo, fhi, normgrad
    integer                :: idim

    stepsize = max(min_stepsize, sqrteps * abs(r))
    r_eval = r

    do idim = 1, NDIM
       ! Sample function at (r - step_idim) and (r + step_idim)
       r_eval(idim) = r(idim) - stepsize(idim)
       flo = f(r_eval)

       r_eval(idim) = r(idim) + stepsize(idim)
       fhi = f(r_eval)

       ! Use central difference scheme
       gradient(idim) = (fhi - flo)/(2 * stepsize(idim))

       ! Reset to original coordinate
       r_eval(idim) = r(idim)
    end do

    normgrad = norm2(gradient)
  end function numerical_gradient_amplitude

  !> Find root of f in the interval [a, b]. If f(a) and f(b) have different
  !> signs, apply bisection directly. Else, first find the (assumed to be)
  !> unique local minimum/maximum to determine a bracket. Return relative
  !> location of root, or 1 if there is no root.
  real(dp) function golden_lsf(box_a, IJK_(a), box_b, IJK_(b), mg)
    type(box_t), intent(in) :: box_a   !< Box a (start point)
    integer, intent(in)     :: IJK_(a) !< Cell-centered index in box a
    type(box_t), intent(in) :: box_b   ! Box b (end point)
    integer, intent(in)     :: IJK_(b) !< Cell-centered index in box b
    type(mg_t), intent(in)  :: mg
    real(dp)                :: a(NDIM), b(NDIM), bracket(NDIM, 2)
    real(dp)                :: r_root(NDIM)
    real(dp), parameter     :: tol      = 1e-8_dp
    integer, parameter      :: max_iter = 100

    a = af_r_cc(box_a, [IJK_(a)])
    b = af_r_cc(box_b, [IJK_(b)])

    golden_lsf = 1.0_dp

    if (get_lsf(a) * get_lsf(b) < 0) then
       r_root = bisection(get_lsf, a, b, tol, max_iter)
    else
       ! Determine bracket by finding local minimum/maximum
       bracket = gss(get_lsf, a, b, &
            minimization=(get_lsf(a) >= 0), tol=tol)

       if (get_lsf(bracket(:, 1)) * get_lsf(a) > 0) then
          return                ! No root
       else
          r_root = bisection(get_lsf, a, bracket(:, 1), tol, max_iter)
       end if
    end if

    golden_lsf = norm2(r_root - a)/norm2(b-a)
  end function golden_lsf

  !> Simple bisection
  function bisection(f, in_a, in_b, tol, max_iter) result(c)
    procedure(mg_func_lsf) :: f
    real(dp), intent(in)   :: in_a(NDIM), in_b(NDIM), tol
    integer, intent(in)    :: max_iter
    real(dp)               :: a(NDIM), b(NDIM), c(NDIM), fc
    integer                :: n

    a = in_a
    b = in_b

    do n = 1, max_iter
       c = 0.5 * (a + b)
       fc = f(c)
       if (0.5 * norm2(b-a) < tol .or. abs(fc) <= 0) exit

       if (fc * f(a) >= 0) then
          a = c
       else
          b = c
       end if
    end do
  end function bisection

  !> Golden-section search. Given a function f with a single local minimum in
  !> the interval [a,b], gss returns a subset interval [c,d] that contains the
  !> minimum with d-c <= tol. Copied from
  !> https://en.wikipedia.org/wiki/Golden-section_search
  function gss(f, in_a, in_b, minimization, tol) result(bracket)
    procedure(mg_func_lsf) :: f
    real(dp), intent(in)  :: in_a(NDIM), in_b(NDIM), tol
    logical, intent(in)   :: minimization
    real(dp)              :: bracket(NDIM, 2)
    real(dp)              :: a(NDIM), b(NDIM), c(NDIM), d(NDIM)
    real(dp)              :: h(NDIM), yc, yd
    real(dp)              :: invphi, invphi2
    integer               :: n, k

    invphi  = (sqrt(5.0_dp) - 1) / 2 ! 1 / phi
    invphi2 = (3 - sqrt(5.0_dp)) / 2 ! 1 / phi^2

    a = in_a
    b = in_b
    h = b - a

    if (norm2(h) <= tol) then
       bracket(:, 1) = a
       bracket(:, 2) = b
       return
    end if

    ! Required steps to achieve tolerance
    n = int(ceiling(log(tol / norm2(h)) / log(invphi)))

    c = a + invphi2 * h
    d = a + invphi * h
    yc = f(c)
    yd = f(d)

    do k = 1, n-1
       if ((yc < yd) .eqv. minimization) then
          b = d
          d = c
          yd = yc
          h = invphi * h
          c = a + invphi2 * h
          yc = f(c)
       else
          a = c
          c = d
          yc = yd
          h = invphi * h
          d = a + invphi * h
          yd = f(d)
       end if
    end do

    if ((yc < yd) .eqv. minimization) then
       bracket(:, 1) = a
       bracket(:, 2) = d
    else
       bracket(:, 1) = c
       bracket(:, 2) = b
    end if
  end function gss

end program poisson_lsf_test
