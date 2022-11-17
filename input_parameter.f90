! This file is automatically created by input_parameter.py
module input_parameter
    use salmon_file, only: open_filehandle
    implicit none

    character(256) :: theory
    character(256) :: sysname
    character(256) :: base_directory
    character(256) :: gs_directory
    character(256) :: read_sbe_gs_bin
    character(256) :: write_sbe_gs_bin
    real(8), dimension(3) :: al
    real(8), dimension(3) :: al_vec1
    real(8), dimension(3) :: al_vec2
    real(8), dimension(3) :: al_vec3
    integer :: nstate
    integer :: nelec
    integer :: nstate_sbe
    integer, dimension(3) :: nkgrid
    integer :: nt
    real(8) :: dt
    real(8) :: e_impulse
    character(256) :: ae_shape1
    character(256) :: ae_shape2
    real(8), dimension(3) :: epdir_re1
    real(8), dimension(3) :: epdir_re2
    real(8), dimension(3) :: epdir_im1
    real(8), dimension(3) :: epdir_im2
    real(8) :: phi_cep1
    real(8) :: phi_cep2
    real(8) :: E_amplitude1
    real(8) :: E_amplitude2
    real(8) :: I_wcm2_1
    real(8) :: I_wcm2_2
    real(8) :: tw1
    real(8) :: tw2
    real(8) :: omega1
    real(8) :: omega2
    real(8) :: t1_t2
    real(8) :: t1_start
    integer :: nenergy
    real(8) :: de
    real(8) :: gamma


contains

    subroutine read_input(icomm)
        use mpi
        implicit none
        integer, intent(in) :: icomm
        integer :: ret, fh
        character(256) :: tmp
        integer :: irank, ierr


        namelist/calculation/ &
        & theory

        namelist/control/ &
        & sysname, &
        & base_directory, &
        & gs_directory, &
        & read_sbe_gs_bin, &
        & write_sbe_gs_bin

        namelist/system/ &
        & al, &
        & al_vec1, &
        & al_vec2, &
        & al_vec3, &
        & nstate, &
        & nelec, &
        & nstate_sbe

        namelist/kgrid/ &
        & nkgrid

        namelist/tgrid/ &
        & nt, &
        & dt

        namelist/emfield/ &
        & e_impulse, &
        & ae_shape1, &
        & ae_shape2, &
        & epdir_re1, &
        & epdir_re2, &
        & epdir_im1, &
        & epdir_im2, &
        & phi_cep1, &
        & phi_cep2, &
        & E_amplitude1, &
        & E_amplitude2, &
        & I_wcm2_1, &
        & I_wcm2_2, &
        & tw1, &
        & tw2, &
        & omega1, &
        & omega2, &
        & t1_t2, &
        & t1_start

        namelist/analysis/ &
        & nenergy, &
        & de, &
        & gamma



        theory = 'perturb_dielec'
        sysname = 'test'
        base_directory = './'
        gs_directory = './'
        read_sbe_gs_bin = 'n'
        write_sbe_gs_bin = 'y'
        al = (/0.0, 0.0, 0.0/)
        al_vec1 = (/0.0, 0.0, 0.0/)
        al_vec2 = (/0.0, 0.0, 0.0/)
        al_vec3 = (/0.0, 0.0, 0.0/)
        nstate = 0
        nelec = 0
        nstate_sbe = 0
        nkgrid = (/0, 0, 0/)
        nt = 1000
        dt = 1.0d-2
        e_impulse = 0.0d0
        ae_shape1 = 'none'
        ae_shape2 = 'none'
        epdir_re1 = (/0.0, 0.0, 0.0/)
        epdir_re2 = (/0.0, 0.0, 0.0/)
        epdir_im1 = (/0.0, 0.0, 0.0/)
        epdir_im2 = (/0.0, 0.0, 0.0/)
        phi_cep1 = 0.0d0
        phi_cep2 = 0.0d0
        E_amplitude1 = 0.0d0
        E_amplitude2 = 0.0d0
        I_wcm2_1 = 0.0d0
        I_wcm2_2 = 0.0d0
        tw1 = 0.0d0
        tw2 = 0.0d0
        omega1 = 0.0d0
        omega2 = 0.0d0
        t1_t2 = 0.0d0
        t1_start = 0.0d0
        nenergy = 1000
        de = 1.0d-3
        gamma = 5.0d-3

        call MPI_COMM_RANK(icomm, irank, ierr)

        if (irank == 0) then
            fh = open_filehandle('.namelist.tmp')
            do while (.true.)
                read(*, '(a)', iostat=ret) tmp
                if (ret < 0) exit ! End of file
                tmp = adjustl(tmp)
                if (tmp(1:1) .ne. '!') write(fh, '(a)') trim(tmp)
            end do
            rewind(fh); read(fh, nml=calculation, iostat=ret)
            rewind(fh); read(fh, nml=control, iostat=ret)
            rewind(fh); read(fh, nml=system, iostat=ret)
            rewind(fh); read(fh, nml=kgrid, iostat=ret)
            rewind(fh); read(fh, nml=tgrid, iostat=ret)
            rewind(fh); read(fh, nml=emfield, iostat=ret)
            rewind(fh); read(fh, nml=analysis, iostat=ret)
            close(fh)

            write(*, '(a)') '# calculation'
            write(*, '(a, a)') '# theory=', trim(theory)
            write(*, '(a)') '# control'
            write(*, '(a, a)') '# sysname=', trim(sysname)
            write(*, '(a, a)') '# base_directory=', trim(base_directory)
            write(*, '(a, a)') '# gs_directory=', trim(gs_directory)
            write(*, '(a, a)') '# read_sbe_gs_bin=', trim(read_sbe_gs_bin)
            write(*, '(a, a)') '# write_sbe_gs_bin=', trim(write_sbe_gs_bin)
            write(*, '(a)') '# system'
            write(*, '(a, 9(es23.15e3))') '# al=', al
            write(*, '(a, 9(es23.15e3))') '# al_vec1=', al_vec1
            write(*, '(a, 9(es23.15e3))') '# al_vec2=', al_vec2
            write(*, '(a, 9(es23.15e3))') '# al_vec3=', al_vec3
            write(*, '(a, 9(i9))') '# nstate=', nstate
            write(*, '(a, 9(i9))') '# nelec=', nelec
            write(*, '(a, 9(i9))') '# nstate_sbe=', nstate_sbe
            write(*, '(a)') '# kgrid'
            write(*, '(a, 9(i9))') '# nkgrid=', nkgrid
            write(*, '(a)') '# tgrid'
            write(*, '(a, 9(i9))') '# nt=', nt
            write(*, '(a, 9(es23.15e3))') '# dt=', dt
            write(*, '(a)') '# emfield'
            write(*, '(a, 9(es23.15e3))') '# e_impulse=', e_impulse
            write(*, '(a, a)') '# ae_shape1=', trim(ae_shape1)
            write(*, '(a, a)') '# ae_shape2=', trim(ae_shape2)
            write(*, '(a, 9(es23.15e3))') '# epdir_re1=', epdir_re1
            write(*, '(a, 9(es23.15e3))') '# epdir_re2=', epdir_re2
            write(*, '(a, 9(es23.15e3))') '# epdir_im1=', epdir_im1
            write(*, '(a, 9(es23.15e3))') '# epdir_im2=', epdir_im2
            write(*, '(a, 9(es23.15e3))') '# phi_cep1=', phi_cep1
            write(*, '(a, 9(es23.15e3))') '# phi_cep2=', phi_cep2
            write(*, '(a, 9(es23.15e3))') '# E_amplitude1=', E_amplitude1
            write(*, '(a, 9(es23.15e3))') '# E_amplitude2=', E_amplitude2
            write(*, '(a, 9(es23.15e3))') '# I_wcm2_1=', I_wcm2_1
            write(*, '(a, 9(es23.15e3))') '# I_wcm2_2=', I_wcm2_2
            write(*, '(a, 9(es23.15e3))') '# tw1=', tw1
            write(*, '(a, 9(es23.15e3))') '# tw2=', tw2
            write(*, '(a, 9(es23.15e3))') '# omega1=', omega1
            write(*, '(a, 9(es23.15e3))') '# omega2=', omega2
            write(*, '(a, 9(es23.15e3))') '# t1_t2=', t1_t2
            write(*, '(a, 9(es23.15e3))') '# t1_start=', t1_start
            write(*, '(a)') '# analysis'
            write(*, '(a, 9(i9))') '# nenergy=', nenergy
            write(*, '(a, 9(es23.15e3))') '# de=', de
            write(*, '(a, 9(es23.15e3))') '# gamma=', gamma
        end if

        call MPI_BCAST(theory, 256, MPI_CHARACTER, icomm, 0, ierr)
        call MPI_BCAST(sysname, 256, MPI_CHARACTER, icomm, 0, ierr)
        call MPI_BCAST(base_directory, 256, MPI_CHARACTER, icomm, 0, ierr)
        call MPI_BCAST(gs_directory, 256, MPI_CHARACTER, icomm, 0, ierr)
        call MPI_BCAST(read_sbe_gs_bin, 256, MPI_CHARACTER, icomm, 0, ierr)
        call MPI_BCAST(write_sbe_gs_bin, 256, MPI_CHARACTER, icomm, 0, ierr)
        call MPI_BCAST(al, 3, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(al_vec1, 3, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(al_vec2, 3, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(al_vec3, 3, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(nstate, 1, MPI_INTEGER, icomm, 0, ierr)
        call MPI_BCAST(nelec, 1, MPI_INTEGER, icomm, 0, ierr)
        call MPI_BCAST(nstate_sbe, 1, MPI_INTEGER, icomm, 0, ierr)
        call MPI_BCAST(nkgrid, 3, MPI_INTEGER, icomm, 0, ierr)
        call MPI_BCAST(nt, 1, MPI_INTEGER, icomm, 0, ierr)
        call MPI_BCAST(dt, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(e_impulse, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(ae_shape1, 256, MPI_CHARACTER, icomm, 0, ierr)
        call MPI_BCAST(ae_shape2, 256, MPI_CHARACTER, icomm, 0, ierr)
        call MPI_BCAST(epdir_re1, 3, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(epdir_re2, 3, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(epdir_im1, 3, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(epdir_im2, 3, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(phi_cep1, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(phi_cep2, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(E_amplitude1, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(E_amplitude2, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(I_wcm2_1, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(I_wcm2_2, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(tw1, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(tw2, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(omega1, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(omega2, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(t1_t2, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(t1_start, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(nenergy, 1, MPI_INTEGER, icomm, 0, ierr)
        call MPI_BCAST(de, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
        call MPI_BCAST(gamma, 1, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
    end subroutine read_input
end module input_parameter
