!-------------------------------------------------------------------------------
!  1D Sod shock tube solver for the compressible Euler equations
!
!  Numerical methods:
!    - godunov : first-order finite-volume method using Rusanov flux
!    - muscl   : second-order MUSCL reconstruction with Van Albada limiter
!    - hybrid  : MUSCL in smooth regions, first-order near strong gradients
!
!  Author: Mohammad Shahid
!-------------------------------------------------------------------------------
program shock_tube_solver
    implicit none

    integer, parameter :: rk = selected_real_kind(15, 300)
    integer, parameter :: neq = 3
    integer :: nx, i, nsteps, argc
    real(rk) :: gamma, x_left, x_right, x0, final_time, cfl, dx, time, dt
    character(len=32) :: scheme
    character(len=256) :: output_file, arg
    real(rk), allocatable :: x(:), u(:,:), unew(:,:), flux(:,:)

    gamma = 1.4_rk
    x_left = 0.0_rk
    x_right = 1.0_rk
    x0 = 0.5_rk
    nx = 201
    final_time = 0.25_rk
    cfl = 0.5_rk
    scheme = 'godunov'
    output_file = 'results/godunov_n201.csv'

    argc = command_argument_count()
    if (argc >= 1) then
        call get_command_argument(1, scheme)
        scheme = adjustl(scheme)
    end if
    if (argc >= 2) then
        call get_command_argument(2, arg)
        read(arg, *) nx
    end if
    if (argc >= 3) then
        call get_command_argument(3, arg)
        read(arg, *) final_time
    end if
    if (argc >= 4) then
        call get_command_argument(4, output_file)
        output_file = adjustl(output_file)
    end if

    if (nx < 20) then
        write(*,*) 'ERROR: nx must be at least 20.'
        stop 1
    end if

    allocate(x(nx), u(neq, nx+4), unew(neq, nx+4), flux(neq, nx+3))

    dx = (x_right - x_left) / real(nx - 1, rk)
    do i = 1, nx
        x(i) = x_left + real(i - 1, rk) * dx
    end do

    call initialize_sod(u, nx, x, x0, gamma)

    time = 0.0_rk
    nsteps = 0
    do while (time < final_time)
        call apply_zero_gradient_bc(u, nx)
        dt = compute_dt(u, nx, dx, cfl, gamma)
        if (time + dt > final_time) dt = final_time - time

        select case (trim(scheme))
        case ('godunov')
            call compute_fluxes(u, flux, nx, gamma, .false., .false.)
        case ('muscl')
            call compute_fluxes(u, flux, nx, gamma, .true., .false.)
        case ('hybrid', 'godmuscl')
            call compute_fluxes(u, flux, nx, gamma, .true., .true.)
        case default
            write(*,*) 'ERROR: unknown scheme: ', trim(scheme)
            write(*,*) 'Use: godunov, muscl, or hybrid'
            stop 1
        end select

        unew = u
        do i = 3, nx+2
            unew(:, i) = u(:, i) - (dt / dx) * (flux(:, i) - flux(:, i-1))
            call enforce_physical_state(unew(:, i), gamma)
        end do

        u = unew
        time = time + dt
        nsteps = nsteps + 1
    end do

    call apply_zero_gradient_bc(u, nx)
    call write_solution(output_file, x, u, nx, gamma, trim(scheme), nsteps, time)
    write(*,'(A,A,A,I0,A,F10.6,A,I0)') 'Finished ', trim(scheme), ' with nx=', nx, ', time=', time, ', steps=', nsteps

contains

    subroutine initialize_sod(u, nx, x, x0, gamma)
        integer, intent(in) :: nx
        real(rk), intent(out) :: u(neq, nx+4)
        real(rk), intent(in) :: x(nx), x0, gamma
        integer :: i, j
        real(rk) :: rho, vel, p

        u = 0.0_rk
        do i = 1, nx
            if (x(i) < x0) then
                rho = 1.0_rk
                vel = 0.0_rk
                p = 1.0_rk
            else
                rho = 0.125_rk
                vel = 0.0_rk
                p = 0.1_rk
            end if
            j = i + 2
            call primitive_to_conserved(rho, vel, p, gamma, u(:, j))
        end do
        call apply_zero_gradient_bc(u, nx)
    end subroutine initialize_sod

    subroutine primitive_to_conserved(rho, vel, p, gamma, q)
        real(rk), intent(in) :: rho, vel, p, gamma
        real(rk), intent(out) :: q(neq)
        q(1) = rho
        q(2) = rho * vel
        q(3) = p / (gamma - 1.0_rk) + 0.5_rk * rho * vel * vel
    end subroutine primitive_to_conserved

    subroutine conserved_to_primitive(q, gamma, rho, vel, p)
        real(rk), intent(in) :: q(neq), gamma
        real(rk), intent(out) :: rho, vel, p
        rho = max(q(1), 1.0e-12_rk)
        vel = q(2) / rho
        p = (gamma - 1.0_rk) * (q(3) - 0.5_rk * rho * vel * vel)
        p = max(p, 1.0e-12_rk)
    end subroutine conserved_to_primitive

    subroutine flux_euler(q, gamma, f)
        real(rk), intent(in) :: q(neq), gamma
        real(rk), intent(out) :: f(neq)
        real(rk) :: rho, vel, p
        call conserved_to_primitive(q, gamma, rho, vel, p)
        f(1) = rho * vel
        f(2) = rho * vel * vel + p
        f(3) = vel * (q(3) + p)
    end subroutine flux_euler

    function max_wave_speed(q, gamma) result(speed)
        real(rk), intent(in) :: q(neq), gamma
        real(rk) :: speed, rho, vel, p, a
        call conserved_to_primitive(q, gamma, rho, vel, p)
        a = sqrt(gamma * p / rho)
        speed = abs(vel) + a
    end function max_wave_speed

    subroutine apply_zero_gradient_bc(u, nx)
        integer, intent(in) :: nx
        real(rk), intent(inout) :: u(neq, nx+4)
        u(:, 1) = u(:, 3)
        u(:, 2) = u(:, 3)
        u(:, nx+3) = u(:, nx+2)
        u(:, nx+4) = u(:, nx+2)
    end subroutine apply_zero_gradient_bc

    function compute_dt(u, nx, dx, cfl, gamma) result(dt)
        integer, intent(in) :: nx
        real(rk), intent(in) :: u(neq, nx+4), dx, cfl, gamma
        real(rk) :: dt, smax
        integer :: i
        smax = 0.0_rk
        do i = 3, nx+2
            smax = max(smax, max_wave_speed(u(:, i), gamma))
        end do
        dt = cfl * dx / max(smax, 1.0e-12_rk)
    end function compute_dt

    function van_albada(r) result(phi)
        real(rk), intent(in) :: r
        real(rk) :: phi
        if (r <= 0.0_rk) then
            phi = 0.0_rk
        else
            phi = (r*r + r) / (r*r + 1.0_rk)
        end if
    end function van_albada

    subroutine limited_slopes(u, slope, nx)
        integer, intent(in) :: nx
        real(rk), intent(in) :: u(neq, nx+4)
        real(rk), intent(out) :: slope(neq, nx+4)
        integer :: i, k
        real(rk) :: dl, dr, r, eps

        eps = 1.0e-12_rk
        slope = 0.0_rk
        do i = 2, nx+3
            do k = 1, neq
                dl = u(k, i) - u(k, i-1)
                dr = u(k, i+1) - u(k, i)
                if (abs(dr) < eps) then
                    r = 0.0_rk
                else
                    r = dl / dr
                end if
                slope(k, i) = van_albada(r) * dr
            end do
        end do
    end subroutine limited_slopes

    subroutine compute_fluxes(u, flux, nx, gamma, use_muscl, use_hybrid)
        integer, intent(in) :: nx
        real(rk), intent(in) :: u(neq, nx+4), gamma
        real(rk), intent(out) :: flux(neq, nx+3)
        logical, intent(in) :: use_muscl, use_hybrid
        real(rk) :: slope(neq, nx+4), ql(neq), qr(neq), fl(neq), fr(neq), alpha
        real(rk) :: sensor, p_left, p_right, rho_tmp, vel_tmp
        integer :: i

        slope = 0.0_rk
        if (use_muscl) call limited_slopes(u, slope, nx)

        do i = 2, nx+2
            if (use_muscl) then
                ql = u(:, i) + 0.5_rk * slope(:, i)
                qr = u(:, i+1) - 0.5_rk * slope(:, i+1)

                if (.not. is_physical(ql, gamma)) ql = u(:, i)
                if (.not. is_physical(qr, gamma)) qr = u(:, i+1)

                if (use_hybrid) then
                    call conserved_to_primitive(u(:, i), gamma, rho_tmp, vel_tmp, p_left)
                    call conserved_to_primitive(u(:, i+1), gamma, rho_tmp, vel_tmp, p_right)
                    sensor = abs(p_right - p_left) / max(abs(p_right + p_left), 1.0e-12_rk)
                    if (sensor > 0.12_rk) then
                        ql = u(:, i)
                        qr = u(:, i+1)
                    end if
                end if
            else
                ql = u(:, i)
                qr = u(:, i+1)
            end if

            call flux_euler(ql, gamma, fl)
            call flux_euler(qr, gamma, fr)
            alpha = max(max_wave_speed(ql, gamma), max_wave_speed(qr, gamma))
            flux(:, i) = 0.5_rk * (fl + fr) - 0.5_rk * alpha * (qr - ql)
        end do
    end subroutine compute_fluxes

    function is_physical(q, gamma) result(ok)
        real(rk), intent(in) :: q(neq), gamma
        logical :: ok
        real(rk) :: rho, vel, p
        call conserved_to_primitive(q, gamma, rho, vel, p)
        ok = (q(1) > 1.0e-10_rk .and. p > 1.0e-10_rk)
    end function is_physical

    subroutine enforce_physical_state(q, gamma)
        real(rk), intent(inout) :: q(neq)
        real(rk), intent(in) :: gamma
        real(rk) :: rho, vel, p
        call conserved_to_primitive(q, gamma, rho, vel, p)
        rho = max(rho, 1.0e-10_rk)
        p = max(p, 1.0e-10_rk)
        call primitive_to_conserved(rho, vel, p, gamma, q)
    end subroutine enforce_physical_state

    subroutine write_solution(filename, x, u, nx, gamma, scheme, nsteps, time)
        integer, intent(in) :: nx, nsteps
        character(len=*), intent(in) :: filename, scheme
        real(rk), intent(in) :: x(nx), u(neq, nx+4), gamma, time
        integer :: i, unit_no
        real(rk) :: rho, vel, p

        unit_no = 20
        open(unit=unit_no, file=trim(filename), status='replace', action='write')
        write(unit_no,'(A)') '# 1D Sod shock tube solution'
        write(unit_no,'(A,A)') '# scheme = ', trim(scheme)
        write(unit_no,'(A,I0)') '# nsteps = ', nsteps
        write(unit_no,'(A,F16.10)') '# time = ', time
        write(unit_no,'(A)') 'x,rho,u,p,E'
        do i = 1, nx
            call conserved_to_primitive(u(:, i+2), gamma, rho, vel, p)
            write(unit_no,'(F14.8,4(",",ES18.10))') x(i), rho, vel, p, u(3, i+2)
        end do
        close(unit_no)
    end subroutine write_solution

end program shock_tube_solver
