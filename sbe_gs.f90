! Ground State Date Storage Module:

module sbe_gs
    use salmon_math, only: pi
    implicit none

    type s_sbe_gs
        !Lattice information
        real(8) :: a_matrix(1:3, 1:3)
        real(8) :: b_matrix(1:3, 1:3)
        real(8) :: volume

        !Ground state (GS) electronic system information
        integer :: nk, nb, ne
        real(8), allocatable :: kpoint(:, :), kweight(:)
        real(8), allocatable :: eigen(:, :)
        real(8), allocatable :: occup(:, :)
        real(8), allocatable :: delta_omega(:, :, :)
        complex(8), allocatable :: p_matrix(:, :, :, :)
        complex(8), allocatable :: d_matrix(:, :, :, :)
        complex(8), allocatable :: tm_p_matrix(:, :, :, :)
        complex(8), allocatable :: rvnl_matrix(:, :, :, :)
        ! complex(8), allocatable :: prod_dk(:, :, :, :, :, :)

        !k-space grid and geometry information
        !NOTE: prepred for uniformally distributed k-grid....
        integer :: nkgrid(1:3)
    end type


contains


subroutine init_sbe_gs(gs, sysname, gs_directory, nkgrid, nb, ne, a1, a2, a3, read_bin, icomm)
    use mpi
    use salmon_file, only: open_filehandle, get_filehandle
    implicit none
    type(s_sbe_gs), intent(inout) :: gs
    character(*), intent(in) :: sysname
    character(*), intent(in) :: gs_directory
    integer, intent(in) :: nkgrid(1:3)
    integer, intent(in) :: nb
    integer, intent(in) :: ne
    real(8), intent(in) :: a1(1:3), a2(1:3), a3(1:3)
    logical, intent(in) :: read_bin
    integer, intent(in) :: icomm 
    integer :: nk, irank, ierr

    call MPI_COMM_RANK(icomm, irank, ierr)

    nk = nkgrid(1) * nkgrid(2) * nkgrid(3)

    gs%nk = nk
    gs%nb = nb
    gs%ne = ne
    gs%nkgrid(1:3) = nkgrid(1:3)

    !Calculate b_matrix, volume_cell and volume_bz from a1..a3 vector.
    call calc_lattice_info()

    allocate(gs%kpoint(1:3, 1:nk))
    allocate(gs%kweight(1:nk))
    allocate(gs%eigen(1:nb, 1:nk))
    allocate(gs%occup(1:nb, 1:nk))
    allocate(gs%delta_omega(1:nb, 1:nb, 1:nk))
    allocate(gs%p_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%d_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%tm_p_matrix(1:nb, 1:nb, 1:3, 1:nk))
    allocate(gs%rvnl_matrix(1:nb, 1:nb, 1:3, 1:nk))

    if (irank == 0) then
        if (read_bin) then
            !Retrieve all data from binray
            write(*,*) "# read_sbe_gs_bin"
            call read_sbe_gs_bin()
        else
            !Retrieve eigenenergies from 'SYSNAME_eigen.data':
            write(*,*) "# read_eigen_data"
            call read_eigen_data()
            !Retrieve k-points from 'SYSNAME_k.data':
            write(*,*) "# read_k_data"
            call read_k_data()
            !Retrieve transition matrix from 'SYSNAME_tm.data':
            write(*,*) "# read_tm_data"
            call read_tm_data()
            !Export all data from binray
            write(*,*) "# save_sbe_gs_bin"
            call save_sbe_gs_bin()
        end if
    end if
    call MPI_BCAST(gs%kpoint, 3*nk, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
    call MPI_BCAST(gs%kweight, nk, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
    call MPI_BCAST(gs%eigen, nb*nk, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
    call MPI_BCAST(gs%occup, nb*nk, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
    call MPI_BCAST(gs%delta_omega, nb*nb*nk, MPI_DOUBLE_PRECISION, icomm, 0, ierr)
    call MPI_BCAST(gs%p_matrix, nb*nb*3*nk, MPI_DOUBLE_COMPLEX, icomm, 0, ierr)
    call MPI_BCAST(gs%d_matrix, nb*nb*3*nk, MPI_DOUBLE_COMPLEX, icomm, 0, ierr)
    call MPI_BCAST(gs%tm_p_matrix, nb*nb*3*nk, MPI_DOUBLE_COMPLEX, icomm, 0, ierr)
    call MPI_BCAST(gs%rvnl_matrix, nb*nb*3*nk, MPI_DOUBLE_COMPLEX, icomm, 0, ierr)
    
    !Calculate omega and d_matrix (neglecting diagonal part):
    if (irank == 0) write(*,*) "# prepare_matrix"
    call prepare_matrix()

    !Initial Occupation Number
    gs%occup(:,:) = 0d0 !!Experimental!!
    gs%occup(1:(ne/2),:) = 2d0 !!Experimental!!

contains

    ! Calculate lattice and reciprocal vectors
    subroutine calc_lattice_info()
        implicit none
        real(8) :: a12(1:3), a23(1:3), a31(1:3), volume
        real(8) :: b1(1:3), b2(1:3), b3(1:3)

        a12(1) = a1(2) * a2(3) - a1(3) * a2(2)
        a12(2) = a1(3) * a2(1) - a1(1) * a2(3)
        a12(3) = a1(1) * a2(2) - a1(2) * a2(1)
        a23(1) = a2(2) * a3(3) - a2(3) * a3(2)
        a23(2) = a2(3) * a3(1) - a2(1) * a3(3)
        a23(3) = a2(1) * a3(2) - a2(2) * a3(1)
        a31(1) = a3(2) * a1(3) - a3(3) * a1(2)
        a31(2) = a3(3) * a1(1) - a3(1) * a1(3)
        a31(3) = a3(1) * a1(2) - a3(2) * a1(1)
        volume = dot_product(a12, a3)
        b1(1:3) = (2d0 * pi / volume) * a23(1:3)
        b2(1:3) = (2d0 * pi / volume) * a31(1:3)
        b3(1:3) = (2d0 * pi / volume) * a12(1:3)
        
        gs%a_matrix(1:3, 1) = a1(1:3)
        gs%a_matrix(1:3, 2) = a2(1:3)
        gs%a_matrix(1:3, 3) = a3(1:3)    
        gs%b_matrix(1, 1:3) = b1(1:3)
        gs%b_matrix(2, 1:3) = b2(1:3)
        gs%b_matrix(3, 1:3) = b3(1:3)
        gs%volume = volume
    end subroutine calc_lattice_info


    ! Read k-point coordinates from SALMON's output file
    subroutine read_k_data()
        implicit none
        character(256) :: dummy
        integer :: fh, ik, iik
        real(8) :: tmp(4)
        fh = open_filehandle(trim(gs_directory) // trim(sysname) // '_k.data', 'old')
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)    
        do ik = 1, nk
            read(fh, *) iik, tmp(1:4)
            if (ik .ne. iik) stop "ik mismatch"
            gs%kpoint(1:3, ik) = tmp(1:3)
            gs%kweight(ik) = tmp(4)
        end do
        close(fh)
    end subroutine read_k_data


    ! Read eigenvalue data from SALMON's output file
    subroutine read_eigen_data()
        implicit none
        character(256) :: dummy
        integer :: fh, i, ik, iik, iib, ib
        real(8) :: tmp(2)

        fh = open_filehandle(trim(gs_directory) // trim(sysname) // '_eigen.data', 'old')
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        do ik = 1, nk
            read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
            do ib = 1, nb
                read(fh, *) iib, tmp(1:2)
                if (ib .ne. iib) stop "ib mismatch"
                gs%eigen(ib, ik) = tmp(1)
                ! gs%occup(ib, ik) = ctmp(2)
            end do
        end do
        close(fh)
    end subroutine read_eigen_data




    ! Read transition dipole moment from SALMON's output file
    subroutine read_tm_data()
        implicit none
        character(256) :: dummy
        integer :: fh, i, ik, ib, jb, iik, iib, jjb
        real(8) :: tmp(1:6)

        fh = open_filehandle(trim(gs_directory) // trim(sysname) // '_tm.data', 'old')
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        do ik = 1, nk
            do ib = 1, nb
                do jb = 1, nb
                    read(fh, *) iik, iib, jjb, tmp(1:6)
                    if (ik .ne. iik) stop "ik mismatch"
                    if (ib .ne. iib) stop "ib mismatch"
                    if (jb .ne. jjb) stop "jb mismatch"
                    gs%tm_p_matrix(ib, jb, 1, ik) = dcmplx(tmp(1), tmp(2))
                    gs%tm_p_matrix(ib, jb, 2, ik) = dcmplx(tmp(3), tmp(4))
                    gs%tm_p_matrix(ib, jb, 3, ik) = dcmplx(tmp(5), tmp(6))
                end do
            end do
        end do
        read(fh, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        do ik = 1, nk
            do ib = 1, nb
                do jb = 1, nb
                    read(fh, *) iik, iib, jjb, tmp(1:6)
                    if (ik .ne. iik) stop "ik mismatch"
                    if (ib .ne. iib) stop "ib mismatch"
                    if (jb .ne. jjb) stop "jb mismatch"
                    gs%rvnl_matrix(ib, jb, 1, ik) = dcmplx(tmp(1), tmp(2))
                    gs%rvnl_matrix(ib, jb, 2, ik) = dcmplx(tmp(3), tmp(4))
                    gs%rvnl_matrix(ib, jb, 3, ik) = dcmplx(tmp(5), tmp(6))
                end do
            end do
        end do
        gs%p_matrix = gs%tm_p_matrix! + gs%rvnl_matrix 
        close(fh)
    end subroutine read_tm_data


    subroutine read_sbe_gs_bin()
        implicit none
        integer :: fh
        ! fh = get_filehandle()
        ! open(fh, file=trim(gs_directory) // trim(sysname) // '_sbe_gs.bin', form='unformatted', status='old')
        ! read(fh) gs%kpoint
        ! read(fh) gs%kweight
        ! read(fh) gs%eigen
        ! read(fh) gs%p_matrix
        ! read(fh) gs%rvnl_matrix
        ! ! read(fh) gs%prod_dk
        ! close(fh)
        ! return
    end subroutine read_sbe_gs_bin


    subroutine save_sbe_gs_bin()
        implicit none
        integer :: fh
        ! fh = get_filehandle()
        ! open(fh, file=trim(gs_directory) // trim(sysname) // '_sbe_gs.bin', form='unformatted', status='replace')
        ! write(fh) gs%kpoint
        ! write(fh) gs%kweight
        ! write(fh) gs%eigen
        ! write(fh) gs%p_matrix
        ! write(fh) gs%rvnl_matrix
        ! ! write(fh) gs%prod_dk
        ! close(fh)
        ! return
    end subroutine save_sbe_gs_bin




    subroutine prepare_matrix()
        implicit none
        integer :: ik, ib, jb
        real(8), parameter :: omega_eps = 1d-9
        complex(8), parameter :: zi = dcmplx(0d0, 1d0)
        do ik=1, nk
            do ib=1, nb
                do jb=1, nb
                    gs%delta_omega(ib, jb, ik) = gs%eigen(ib, ik) - gs%eigen(jb, ik)
                    if (omega_eps < abs(gs%delta_omega(ib, jb, ik))) then
                        ! gs%d_matrix(ib, jb, 1:3, ik) = &
                        !     & (zi * gs%p_matrix(ib, jb, 1:3, ik) - gs%rvnl_matrix(ib, jb, 1:3, ik)) &
                        !     & / gs%delta_omega(ib, jb, ik)  
                        gs%d_matrix(ib, jb, 1:3, ik) = &
                            & zi * (gs%p_matrix(ib, jb, 1:3, ik) +  gs%rvnl_matrix(ib, jb, 1:3, ik)) &
                            & / gs%delta_omega(ib, jb, ik)  
                        ! gs%d_matrix(ib, jb, 1:3, ik) = &
                        !     & zi * (gs%p_matrix(ib, jb, 1:3, ik)) &
                        !     & / gs%delta_omega(ib, jb, ik)  
                    else
                        gs%d_matrix(ib, jb, 1:3, ik) = 0d0
                    end if
                end do
            end do
        end do
    end subroutine prepare_matrix

    
end subroutine init_sbe_gs

end module sbe_gs

