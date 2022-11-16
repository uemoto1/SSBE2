module sbe_solver
    use salmon_math, only: pi
    use sbe_gs
    implicit none



    type s_sbe_bloch_solver
        !k-points for real-time SBE calculation
        integer :: nk, nb
        integer :: ik_max, ik_min
        complex(8), allocatable :: rho(:, :, :)
    end type



contains



subroutine init_sbe(sbe, gs, nb_sbe, icomm)
    use mpi
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    integer, intent(in) :: nb_sbe
    integer, intent(in) :: icomm
    integer :: ik, ib, nk_proc, irank, nproc, ierr

    call MPI_COMM_SIZE(icomm, nproc, ierr)
    call MPI_COMM_RANK(icomm, irank, ierr)

    sbe%nk = gs%nk
    sbe%nb = nb_sbe

    nk_proc = (sbe%nk - 1) / nproc + 1
    sbe%ik_min = irank * nk_proc + 1
    sbe%ik_max = min(sbe%ik_min + nk_proc - 1, sbe%nk)

    allocate(sbe%rho(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max))
    
    sbe%rho(:, :, :) = 0d0
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, sbe%nb
            sbe%rho(ib, ib, ik) = gs%occup(ib, ik)
        end do
    end do
end subroutine


subroutine calc_current_bloch(sbe, gs, Ac, jmat, icomm)
    use mpi
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(out) :: jmat(1:3)
    integer, intent(in) :: icomm
    integer :: ik, idir, ib, jb, ierr
    complex(8) :: jtmp(1:3)
    complex(8), parameter :: zI = dcmplx(0.0d0, 1.0d0)

    jtmp(1:3) = 0d0

    !$omp parallel do default(shared) private(ik,ib,jb,idir) reduction(+:jtmp)
    do ik = sbe%ik_min, sbe%ik_max
        do idir = 1, 3
            do ib = 1, sbe%nb
                do jb = 1, sbe%nb
                    jtmp(idir) = jtmp(idir) + gs%kweight(ik) * sbe%rho(jb, ib, ik) * ( &
                        & gs%p_matrix(ib, jb, idir, ik) & !+ zI * gs%rvnl_matrix(ib, jb, idir, ik) &
                        & )
                end do
            end do
        end do
    end do
    !$omp end parallel do
    call MPI_ALLREDUCE(MPI_IN_PLACE, jtmp, 3, MPI_DOUBLE_COMPLEX, MPI_SUM, icomm, ierr)
    
    jtmp(1:3) = jtmp(1:3) / sum(gs%kweight(:))

    jmat(:) = (real(jtmp(:)) + Ac * calc_trace(sbe, gs, sbe%nb, MPI_COMM_WORLD)) / gs%volume    
    !jmat(1:3) = (real(jtmp(1:3))) / gs%volume    
    return
end subroutine calc_current_bloch










subroutine dt_evolve_bloch(sbe, gs, Ac, dt)
    implicit none
    type(s_sbe_bloch_solver), intent(inout) :: sbe
    type(s_sbe_gs), intent(inout) :: gs
    real(8), intent(in) :: Ac(1:3)
    real(8), intent(in) :: dt
    complex(8), parameter :: zi = dcmplx(0d0, 1d0)
    integer :: nb, nk

    complex(8) :: hrho1(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max)
    complex(8) :: hrho2(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max)
    complex(8) :: hrho3(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max)
    complex(8) :: hrho4(1:sbe%nb, 1:sbe%nb, sbe%ik_min:sbe%ik_max)

    nb = sbe%nb
    nk = sbe%nk

    call calc_hrho_bloch(sbe%rho, hrho1)
    call calc_hrho_bloch(hrho1, hrho2)
    call calc_hrho_bloch(hrho2, hrho3)
    call calc_hrho_bloch(hrho3, hrho4)

    sbe%rho = sbe%rho + hrho1 * (- zi * dt)
    sbe%rho = sbe%rho + hrho2 * (- zi * dt) ** 2 * (1d0 / 2d0)
    sbe%rho = sbe%rho + hrho3 * (- zi * dt) ** 3 * (1d0 / 6d0)
    sbe%rho = sbe%rho + hrho4 * (- zi * dt) ** 4 * (1d0 / 24d0)
    return

contains


    !Calculate [H, rho] commutation:
    subroutine calc_hrho_bloch(rho, hrho)
        implicit none
        complex(8), intent(in)     :: rho(nb, nb, sbe%ik_min:sbe%ik_max)
        complex(8), intent(out)    :: hrho(nb, nb, sbe%ik_min:sbe%ik_max)
        integer :: ik, idir
        !$omp parallel do default(shared) private(ik,idir)
        do ik = sbe%ik_min, sbe%ik_max
            hrho(1:nb, 1:nb, ik) = gs%delta_omega(1:nb, 1:nb, ik) * rho(1:nb, 1:nb, ik)
            !hrho = hrho + Ac(t) * (p * rho - rho * p)
            do idir=1, 3 !1:x, 2:y, 3:z
                ! hrho(1:nb, 1:nb, ik) = hrho(1:nb, 1:nb, ik) + Ac(idir) * (&
                ! & + matmul(gs%p_matrix(1:nb, 1:nb, idir, ik), rho(1:nb, 1:nb, ik)) &
                ! & - matmul(rho(1:nb, 1:nb, ik), gs%p_matrix(1:nb, 1:nb, idir, ik)) &
                ! & )


                call ZGEMM("N","N", sbe%nb, sbe%nb, sbe%nb, &
                    dcmplx(+Ac(idir), 0d0), &
                    gs%p_matrix(:, :, idir, ik),sbe%nb, &
                    rho(:, :, ik), sbe%nb, &
                    dcmplx(1d0, 0d0),hrho(:, :, ik),sbe%nb)

                call ZGEMM("N","N", sbe%nb, sbe%nb, sbe%nb, &
                    dcmplx(-Ac(idir), 0d0), &
                    rho(:, :, ik), sbe%nb, &
                    gs%p_matrix(:, :, idir, ik),sbe%nb, &
                    dcmplx(1d0, 0d0), hrho(:, :, ik), sbe%nb)
                    
            end do !idir
        end do !ik
        !$omp end parallel do
        return
    end subroutine calc_hrho_bloch
end subroutine

function calc_trace(sbe, gs, nb_max, icomm) result(tr)
    use mpi
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    integer, intent(in) :: icomm
    integer, intent(in) :: nb_max
    complex(8), parameter :: zi = dcmplx(0d0, 1d0)
    integer :: ik, ib, ierr
    real(8) :: tr, tr_tmp
    tr_tmp = 0d0
    !$omp parallel do default(shared) private(ik, ib) reduction(+: tr_tmp) collapse(2) 
    do ik = sbe%ik_min, sbe%ik_max
        do ib = 1, nb_max
            tr_tmp = tr_tmp + real(sbe%rho(ib, ib, ik)) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    call MPI_ALLREDUCE(MPI_IN_PLACE, tr_tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr)
    tr = tr_tmp / sum(gs%kweight)
    return
end function calc_trace


function calc_energy(sbe, gs, Ac, icomm) result(energy)
    use mpi
    implicit none
    type(s_sbe_bloch_solver), intent(in) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    integer, intent(in) :: icomm
    real(8), intent(in) :: Ac(1:3)
    integer :: ik, ib, jb, idir, ierr
    real(8) :: energy
    ! real(8) :: kvec(1:3)
    energy = 0d0
    !$omp parallel do default(shared) private(ik, ib, jb, idir) reduction(+: energy)
    do ik = sbe%ik_min, sbe%ik_max
        ! kvec(1:3) = gs%kpoint(1, ik) * gs%b_matrix(1, 1:3) &
        !     & + gs%kpoint(2, ik) * gs%b_matrix(2, 1:3) &
        !     & + gs%kpoint(3, ik) * gs%b_matrix(3, 1:3)
        do ib = 1, sbe%nb
            do idir = 1, 3
                do jb = 1, sbe%nb
                    energy = energy &
                        & + Ac(idir) * real(sbe%rho(ib, jb, ik) * gs%p_matrix(jb, ib, idir, ik)) * gs%kweight(ik)
                end do
            end do
            energy = energy &
                & + real(sbe%rho(ib, ib, ik)) * ( &
                & + gs%eigen(ib, ik) &
                !& + dot_product(kvec(:), Ac(:))
                & + 0.5 * dot_product(Ac, Ac) &
                & ) * gs%kweight(ik)
        end do
    end do
    !$omp end parallel do
    call MPI_ALLREDUCE(MPI_IN_PLACE, energy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, icomm, ierr)
    energy = energy / sum(gs%kweight)

    return 
end function calc_energy


end module



