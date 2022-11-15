program main
    use mpi
    use sbe_gs
    use sbe_solver
    use input_parameter
    use test
    use em_field
    implicit none

    type(s_sbe_bloch_solver) :: sbe
    type(s_sbe_gs) :: gs
    real(8) :: t,  E(3), jmat(3)
    real(8), allocatable :: Ac_ext_t(:, :)
    integer :: it
    real(8) :: energy0, energy
    real(8) :: tr_all, tr_vb
    integer :: icomm, nproc, irank, ierr

    icomm = MPI_COMM_WORLD
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(icomm, nproc, ierr)
    call MPI_COMM_RANK(icomm, irank, ierr)
    if (irank == 0) write(*, "(a,i6)") "# MPI initialized: nproc=", nproc

    call read_input(MPI_COMM_WORLD)

    if (0.0d0 < al(1)) al_vec1(1:3) = (/ al(1), 0.0d0, 0.0d0 /)
    if (0.0d0 < al(2)) al_vec2(1:3) = (/ 0.0d0, al(2), 0.0d0 /)
    if (0.0d0 < al(3)) al_vec3(1:3) = (/ 0.0d0, 0.0d0, al(3) /)
    if (nstate_sbe < 1) nstate_sbe = nstate

    ! Read ground state electronic system:
    call init_sbe_gs(gs, sysname, gs_directory, &
        & nkgrid, nstate, nelec, &
        & al_vec1, al_vec2, al_vec3, &
        & .false., MPI_COMM_WORLD)

        
    

    ! Calculate dielectric spectra and save as SYSNAME_dielec.data:
    if (trim(theory) == 'lr_dielec') then
        if (irank == 0) call calc_dielec(sysname, base_directory, gs, nenergy, de, gamma)
        stop
    end if


    ! Initialization of SBE solver and density matrix:
    call init_sbe(sbe, gs, nstate_sbe, MPI_COMM_WORLD)

    ! Prepare external pulse
    allocate(Ac_ext_t(1:3, -1:nt+1))
    call calc_Ac_ext_t(0.0d0, dt, 0, nt, Ac_ext_t)

    ! Realtime calculation
    open(unit=100, file=trim(base_directory)//trim(sysname)//"_sbe_rt.data")
    write(100, '(4a)') "# 1:Time[a.u.] 2:Ac_ext_x[a.u.] 3:Ac_ext_y[a.u.] 4:Ac_ext_z[a.u.] ", &
        & "5:E_ext_x[a.u.] 6:E_ext_y[a.u.] 7:E_ext_z[a.u.] 8:Ac_tot_x[a.u.] ", &
        & "9:Ac_tot_y[a.u.] 10:Ac_tot_z[a.u.] 11:E_tot_x[a.u.] 12:E_tot_y[a.u.] ", &
        & "13:E_tot_z[a.u.]  14:Jm_x[a.u.] 15:Jm_y[a.u.] 16:Jm_z[a.u.]"

    open(unit=101, file=trim(base_directory)//trim(sysname)//"_sbe_rt_energy.data")
    write(101, '(a)') "# 1:Time[a.u.] 2:Eall[a.u.] 3:Eall-Eall0[a.u.]"

    open(unit=102, file=trim(base_directory)//trim(sysname)//"_sbe_nex.data")
    write(102, '(a)') "# 1:Time[a.u.] 2:nelec[a.u.] 3:nhole[a.u.]"

    energy0 = calc_energy(sbe, gs, Ac_ext_t(:, 0), MPI_COMM_WORLD)
    write(101, '(f12.6,2(es24.15e3))') 0.0d0, energy0, 0.0d0

    do it = 1, nt
        t = dt * it
        call dt_evolve_bloch(sbe, gs, Ac_ext_t(:, it), dt)
        
        if (mod(it, 10) == 0) then
            E(:) = (Ac_ext_t(:, it + 1) - Ac_ext_t(:, it - 1)) / (2 * dt)
            call calc_current_bloch(sbe, gs, Ac_ext_t(:, it), Jmat, MPI_COMM_WORLD)
            write(100, '(f12.6,15(es24.15e3))') t, Ac_ext_t(:, it), E(:), Ac_ext_t(:, it), E(:), Jmat(:)

            energy = calc_energy(sbe, gs, Ac_ext_t(:, it), MPI_COMM_WORLD)
            write(101, '(f12.6,2(es24.15e3))') t, energy, energy-energy0

            tr_all = calc_trace(sbe, gs, nstate_sbe, MPI_COMM_WORLD)
            tr_vb = calc_trace(sbe, gs, nelec / 2, MPI_COMM_WORLD)
            write(102, '(f12.6,2(es24.15e3))') t, tr_all - tr_vb, nelec - tr_vb 
            
            write(*, '(f12.6,es24.15e3)') t, tr_all
        end if

    end do

    close(100)
    close(101)
    close(102)


    call MPI_FINALIZE(ierr)
    stop "Bye!"

    
end program 
