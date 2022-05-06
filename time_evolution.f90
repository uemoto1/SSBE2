module time_evolution
    use ground_state, only: gs_data
    implicit none

    type :: rt_data
        integer :: nk
        integer :: nstate
        real(8) :: dt
        complex(8), allocatable :: rho(:, :, :)
    end type

contains

subroutine init_bloch(rt, gs)
    implicit none
    type(rt_data), intent(out) :: rt
    type(gs_data), intent(in) :: gs
    integer :: ik, ib
    rt%nk = gs%nk
    rt%nstate = gs%nstate
    allocate(rt%rho(rt%nstate, rt%nstate, rt%nk))
    do ik = 1, rt%nk
        do ib = 1, rt%nstate
            rt%rho(ib, ib, ik) = gs%occup(ib, ik)
        end do
    end do
    return
end subroutine init_bloch

subroutine dt_evolve_bloch(rt, gs, dt, Ac0, Ac1)
    implicit none
    type(rt_data), intent(inout) :: rt
    type(gs_data), intent(in) :: gs
    real(8), intent(in) :: dt, Ac0(3), Ac1(3)
    complex(8) :: drho_k1(rt%nstate, rt%nstate)
    complex(8) :: rho_k2(rt%nstate, rt%nstate)
    complex(8) :: drho_k2(rt%nstate, rt%nstate)
    ! allocate(drho1(rt%nstate, rt%nstate, rt%nk))
    integer :: ik

    do ik = 1, rt%nk
        call calc_drho_k(drho_k1, rt%rho(:, :, ik), &
            & gs%pmatrix(:, :, :, ik), gs%rvnl(:, :, :, ik), Ac0)
        rho_k2 = rt%rho(:, :, ik) + dt * drho_k1
        call calc_drho_k(drho_k2, rho_k2, &
            & gs%pmatrix(:, :, :, ik), gs%rvnl(:, :, :, ik), Ac1)
        rt%rho(:, :, ik) = rt%rho(:, :, ik) + 0.5 * dt * (drho_k1 + drho_k2)
    end do


!     write(*,*) -2; flush(0)
!     rho2 = rt%rho + dt * drho1
!     call calc_drho(drho2, rho2, Ac0)
!     rt%rho = rt%rho + 0.5d0 * dt * (drho1 + drho2)
contains

subroutine calc_drho_k(drho_k, rho_k, p_k, rvnl_k, Ac)
    implicit none
    complex(8), intent(out) :: drho_k(rt%nstate, rt%nstate)
    complex(8), intent(in) :: rho_k(rt%nstate, rt%nstate)
    complex(8), intent(in) :: p_k(rt%nstate, rt%nstate, 3)
    complex(8), intent(in) :: rvnl_k(rt%nstate, rt%nstate, 3)
    real(8), intent(in) :: Ac(3)
    integer :: i, j, l, n
    real(8) :: t

    !!!$omp parallel do collapse(2) default(share) private(i, j)
    do j = 1, rt%nstate
        do i = 1, rt%nstate
            drho_k(i, j) = 0.0d0
            do n = 1, 3
                do l = 1, rt%nstate
                    drho_k(i, j) = drho_k(i, j) + dcmplx(0.0, -1.0) * Ac(n) &
                        & * ((p_k(i, l, n) + rvnl_k(i, l, n)) * rho_k(l, j) &
                        & -  rho_k(i, l) * (p_k(l, j, n) + rvnl_k(l, j, n)))
                end do
            end do
        end do
    end do
    !!!$omp end parallel do
    
end subroutine calc_drho_k
end subroutine dt_evolve_bloch


real(8) function calc_total(rt, gs)
    implicit none
    type(rt_data), intent(in) :: rt
    type(gs_data), intent(in) :: gs
    integer :: ik, i

    calc_total = 0.0d0
    do ik = 1, rt%nk
        do i = 1, rt%nstate
            calc_total = calc_total + gs%kweight(ik) * real(rt%rho(i, i, ik))
        end do
    end do
    return
end function
                




end module time_evolution