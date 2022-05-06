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
    complex(8) :: drho1(rt%nstate, rt%nstate, rt%nk)
    complex(8) :: rho2(rt%nstate, rt%nstate, rt%nk)
    complex(8) :: drho2(rt%nstate, rt%nstate, rt%nk)

    ! Modified Euler
    call calc_drho(drho1, rt%rho, Ac0)
    rho2 = rt%rho + dt * drho1
    call calc_drho(drho2, rho2, Ac0)
    rt%rho = rt%rho + 0.5d0 * dt * (drho1 + drho2)
contains

subroutine calc_drho(drho, rho, Ac)
    implicit none
    complex(8), intent(out) :: drho(rt%nstate, rt%nstate, rt%nk)
    complex(8), intent(in) :: rho(rt%nstate, rt%nstate, rt%nk)
    real(8), intent(in) :: Ac(3)
    integer :: i, j, l, n, ik
    
    drho(:, :, :) = 0.0d0
    !!$omp parallel do default(shared) private(ik,n,i,j,l) collapse(3)
    do ik = 1, rt%nk
        do i = 1, rt%nstate
            do j = 1, rt%nstate
                do n = 1, 3
                    do l = 1, rt%nstate
                        drho(i, j, ik) = drho(i, j, ik) + dcmplx(0.0, -1.0) * Ac(n) &
                            & * ((gs%pmatrix(i, l, n, ik) + gs%rvnl(i, l, n, ik)) * rho(l, j, ik) &
                            & -  rho(i, l, ik) * (gs%pmatrix(l, j, n, ik) + gs%rvnl(l, j, n, ik)))
                    end do
                end do
            end do
        end do
    end do
    !!$omp end parallel do
    return
end subroutine calc_drho
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