program main
    use sbe_gs
    use sbe_solver
    use input_parameter
    use test
    use em_field
    implicit none

    type(s_sbe_bloch_solver) :: sbe
    type(s_sbe_gs) :: gs
    real(8) :: t,  Ac(3), jmat(3)
    real(8), allocatable :: Ac_ext_t(:, :)
    integer :: it
    real(8) ::energy0, energy

    call read_input()

    ! Read ground state electronic system:
    call init_sbe_gs(gs, sysname, gs_directory, &
        & nkgrid, nstate, nelec, &
        & al_vec1, al_vec2, al_vec3, &
        & (read_sbe_gs_bin .eq. 'y'))

    ! Calculate dielectric spectra and save as SYSNAME_dielec.data:
    if (trim(theory) == 'perturb_dielec') then
        call calc_dielec(sysname, base_directory, gs, nenergy, de, gamma)
        stop
    end if

    if (nstate_sbe < 1) nstate_sbe = nstate

    ! Initialization of SBE solver and density matrix:
    call init_sbe(sbe, gs, nstate_sbe)

    ! Prepare external pulse
    allocate(Ac_ext_t(1:3, -1:nt+1))
    call calc_Ac_ext_t(0.0d0, dt, 0, nt, Ac_ext_t)

    write(*, '("#",99(1x,a))') "1:Step", "2:Time[au]", "3:Ac_x", "4:Ac_y", "5:Ac_z", &
        & "6:E_x", "7:E_y", "8:E_z", "9:Jmat_x", "10:Jmat_y", "11:Jmat_z", "12:n_v", "13:n_all" 

    ! Realtime calculation

    open(unit=999, file=trim(base_directory)//trim(sysname)//"_rt_energy.data")
    energy0 = calc_energy(sbe, gs, Ac_ext_t(:, 0))
    write(999, '(f12.6,2(es24.15e3))') 0.0d0, energy0, 0.0d0

    do it = 1, nt
        t = dt * it
        call dt_evolve_bloch(sbe, gs, Ac_ext_t(:, it), dt)
        if (mod(it, 10) == 0) then
            call calc_current_bloch(sbe, gs, Ac, Jmat)
            write(*, "(10(es24.15e3))") t, Ac_ext_t(:, it), &
            & - (Ac_ext_t(:, it+1)-Ac_ext_t(:, it-1)) / (2*dt), &
            & jmat! calc_trace(sbe, gs, sbe%nb) , calc_trace(sbe, gs, sbe%nb/2)
        end if

        if (mod(it, 10) == 0) then
            energy = calc_energy(sbe, gs, Ac_ext_t(:, it))
            write(999, '(f12.6,2(es24.15e3))') t, energy, energy-energy0
        end if



    !     if (mod(it, out_rt_step) == 0) then
    !         call calc_cos2_pulse(t, pulse_tw1, &
    !             & rlaser_int_wcm2_1, omega1, phi_cep1, epdir_re1, epdir_im1, &
    !             & Ac, E)
    !         call calc_current(sbe, gs, Ac, jmat)
    !         write(*, '(i6,1x,f9.3,99(1x,es23.15e3))') &
    !         & it, t, Ac, E, jmat, calc_trace(sbe, sbe%nb/2), calc_trace(sbe, sbe%nb)
    !     end if
    end do

    close(999)

    stop
end program 
