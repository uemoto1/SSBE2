module ground_state
    implicit none
    
    type :: gs_data
        ! Number of k-points
        integer :: nk
        ! Number of states
        integer :: nstate
        ! Reduced coordinate of kpoints
        real(8), allocatable :: kpoint(:, :)
        ! Integration weight of kpoints
        real(8), allocatable :: kweight(:)
        ! Eigenvalue
        real(8), allocatable :: eigen(:, :)
        ! Occupation
        real(8), allocatable :: occup(:, :)
        ! Momentum matrix
        complex(8), allocatable :: pmatrix(:, :, :, :)
        ! [r, Vnl] matrix element
        complex(8), allocatable :: rvnl(:, :, :, :)
        logical :: use_nonlocal_potential
        ! Volume
        real(8) :: volume
    end type

contains

subroutine load_salmon_data(gs, nk, nstate, sysname, base_directory)
    implicit none
    type(gs_data), intent(inout) :: gs
    integer, intent(in) :: nk
    integer, intent(in) :: nstate
    character(*), intent(in) :: sysname
    character(*), intent(in) :: base_directory
    character(99) :: file_k_data
    character(99) :: file_eigen_data
    character(99) :: file_tm_data
    character(99) :: dummy
    integer :: ik, iik, ib, iib, jb, jjb
    integer :: itmp
    real(8) :: ctmp(6)

    gs%nk = nk
    gs%nstate = nstate
    allocate(gs%kpoint(3, nk))
    allocate(gs%kweight(nk))
    allocate(gs%eigen(nstate, nk))
    allocate(gs%occup(nstate, nk))
    allocate(gs%pmatrix(nstate, nstate, 3, nk))
    allocate(gs%rvnl(nstate, nstate, 3, nk))

    file_k_data = trim(base_directory) // "/" // trim(sysname) // "_k.data"
    write(*, "(a,a)") "# Open:", trim(file_k_data); flush(0)

    open(10, file=trim(file_k_data), action="read")
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    do ik = 1, nk
        read(10, *) iik, ctmp(1:4)
        if (ik .ne. iik) stop "ik mismatch"
        gs%kpoint(1:3, ik) = ctmp(1:3)
        gs%kweight(ik) = ctmp(4)
    end do
    close(10)

    file_eigen_data = trim(base_directory) // "/" // trim(sysname) // "_eigen.data"
    write(*, "(a,a)") "# Open:", trim(file_eigen_data); flush(0)

    open(10, file=trim(file_eigen_data), action="read")
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    do ik = 1, nk
        read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        do ib = 1, nstate
            read(10, *) iib, ctmp(1:2)
            if (ib .ne. iib) stop "ib mismatch"
            gs%eigen(ib, ik) = ctmp(1)
            gs%occup(ib, ik) = ctmp(2)
        end do
    end do
    close(10)

    file_tm_data = trim(base_directory) // "/" // trim(sysname) // "_tm.data"
    write(*, "(a,a)") "# Open:", trim(file_tm_data); flush(0)

    open(10, file=trim(file_tm_data), action="read")
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    do ik = 1, nk
        do ib = 1, nstate
            do jb = 1, nstate
                read(10, *) iik, iib, jjb, ctmp(1:6)
                if (ik .ne. iik) stop "ik mismatch"
                if (ib .ne. iib) stop "ib mismatch"
                if (jb .ne. jjb) stop "jb mismatch"
                gs%pmatrix(ib, jb, 1, ik) = dcmplx(ctmp(1), ctmp(2))
                gs%pmatrix(ib, jb, 2, ik) = dcmplx(ctmp(3), ctmp(4))
                gs%pmatrix(ib, jb, 3, ik) = dcmplx(ctmp(5), ctmp(6))
            end do
        end do
    end do
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    gs%use_nonlocal_potential = .true.
    do ik = 1, nk
        do ib = 1, nstate
            do jb = 1, nstate
                read(10, *) iik, iib, jjb, ctmp(1:6)
                if (ik .ne. iik) stop "ik mismatch"
                if (ib .ne. iib) stop "ib mismatch"
                if (jb .ne. jjb) stop "jb mismatch"
                gs%rvnl(ib, jb, 1, ik) = dcmplx(ctmp(1), ctmp(2))
                gs%rvnl(ib, jb, 2, ik) = dcmplx(ctmp(3), ctmp(4))
                gs%rvnl(ib, jb, 3, ik) = dcmplx(ctmp(5), ctmp(6))
            end do
        end do
    end do
    close(10)

    return
end subroutine load_salmon_data



subroutine load_elk_data(gs, nk, nstate, volume, base_directory)
    implicit none
    type(gs_data), intent(inout) :: gs
    integer, intent(in) :: nk
    integer, intent(in) :: nstate
    real(8), intent(in) :: volume
    character(*), intent(in) :: base_directory
    character(99) :: file_kpoints
    character(99) :: file_eigval
    character(99) :: file_pmat
    character(99) :: dummy
    integer :: ik, iik, ib, iib, jb, jjb
    integer :: itmp, recl
    real(8) :: ctmp(6)
    complex(8), allocatable :: ptmp(:, :, :)

    gs%nk = nk
    gs%nstate = nstate
    gs%volume = volume
    allocate(gs%kpoint(3, nk))
    allocate(gs%kweight(nk))
    allocate(gs%eigen(nstate, nk))
    allocate(gs%occup(nstate, nk))
    allocate(gs%pmatrix(nstate, nstate, 3, nk))
    allocate(gs%rvnl(nstate, nstate, 3, nk))

    file_kpoints = trim(base_directory) // "/" // "KPOINTS.OUT"
    write(*, "(a,a)") "# Open:", trim(file_kpoints); flush(0)

    open(10, file=trim(file_kpoints), action="read")
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    do ik = 1, nk
        read(10, *) iik, ctmp(1:4), itmp
        if (ik .ne. iik) stop "ik mismatch!"
        gs%kpoint(:, ik) = (/ ctmp(1), ctmp(2), ctmp(3) /)
        gs%kweight(ik) = ctmp(4)
    end do
    close(10)

    file_eigval = trim(base_directory) // "/" // "EIGVAL.OUT"
    write(*, "(a,a)") "# Open:", trim(file_eigval); flush(0)

    open(10, file=trim(file_eigval), action="read")
    read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
    do ik = 1, nk
        read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        read(10, "(a)") dummy; write(*, "('#>',4x,a)") trim(dummy)
        do ib = 1, nstate
            read(10, *) iib, ctmp(1:2)
            if (ib .ne. iib) stop "ib mismatch!"
            gs%eigen(ib, ik) = ctmp(1)
            gs%occup(ib, ik) = ctmp(2)
        end do
    end do

    file_pmat = trim(base_directory) // "/" // "PMAT.OUT"
    write(*, "(a,a)") "# Open:", trim(file_pmat); flush(0)

    allocate(ptmp(nstate, nstate, 3))
    inquire(iolength=recl) ctmp(1:3), itmp, ptmp
    
    open(10, file=trim(file_pmat), action="read", form="unformatted", access="direct", recl=recl)
    do ik = 1, nk
        read(10, rec=ik) ctmp(1:3), itmp, ptmp
        if (itmp .ne. nstate) stop "nstate mismatch!"
        if (abs(ctmp(1) - gs%kpoint(1, ik)) > 1d-9) stop "k1 mismatch"
        if (abs(ctmp(2) - gs%kpoint(2, ik)) > 1d-9) stop "k2 mismatch"
        if (abs(ctmp(3) - gs%kpoint(3, ik)) > 1d-9) stop "k3 mismatch"
        gs%pmatrix(:, :, :, ik) = ptmp(:, :, :)
    end do
    gs%rvnl(:, :, :, :) = 0.0d0
    gs%use_nonlocal_potential = .false.
    deallocate(ptmp)
    close(10)
end subroutine
end module ground_state




