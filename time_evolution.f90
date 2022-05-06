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

    !$omp parallel do default(shared) private(ik, drho_k1, rho_k2, drho_k2)
    do ik = 1, rt%nk
        call calc_drho_k(drho_k1, rt%rho(:, :, ik), &
            & gs%eigen(:, ik), gs%pmatrix(:, :, :, ik), gs%rvnl(:, :, :, ik), Ac0)
        rho_k2 = rt%rho(:, :, ik) + dt * drho_k1
        call calc_drho_k(drho_k2, rho_k2, &
            & gs%eigen(:, ik), gs%pmatrix(:, :, :, ik), gs%rvnl(:, :, :, ik), Ac1)
        rt%rho(:, :, ik) = rt%rho(:, :, ik) + 0.5 * dt * (drho_k1 + drho_k2)
    end do
    !$omp end parallel do

contains

subroutine calc_drho_k(drho_k, rho_k, e_k, p_k, rvnl_k, Ac)
    implicit none
    complex(8), intent(out) :: drho_k(rt%nstate, rt%nstate)
    complex(8), intent(in) :: rho_k(rt%nstate, rt%nstate)
    real(8), intent(in) :: e_k(rt%nstate)
    complex(8), intent(in) :: p_k(rt%nstate, rt%nstate, 3)
    complex(8), intent(in) :: rvnl_k(rt%nstate, rt%nstate, 3)
    real(8), intent(in) :: Ac(3)
    integer :: i, j, l, n

    !$omp parallel do default(shared) private(i,j) collapse(2)
    do j = 1, rt%nstate
        do i = 1, rt%nstate
            drho_k(i, j) = dcmplx(0.0, -1.0) * (e_k(i) - e_k(j)) * rho_k(i, j)
            do n = 1, 3
                do l = 1, rt%nstate
                    drho_k(i, j) = drho_k(i, j) + dcmplx(0.0, -1.0) * Ac(n) * ( &
                        & (p_k(i, l, n) + rvnl_k(i, l, n)) * rho_k(l, j) &
                        & -  rho_k(i, l) * (p_k(l, j, n) + rvnl_k(l, j, n)))
                end do
            end do
        end do
    end do
    !$omp end parallel do
    
end subroutine calc_drho_k
end subroutine dt_evolve_bloch

subroutine current(jcur, rt, gs, Ac)
    implicit none
    real(8), intent(out) :: jcur(3)
    type(rt_data), intent(in) :: rt
    type(gs_data), intent(in) :: gs
    real(8), intent(in) :: Ac(3)
    integer :: ik, ib, jb

    jcur(:) = 0.0d0
    !$omp parallel do default(shared) private(ik,ib,jb) reduction(+:jcur) collapse(2)
    do ik = 1, rt%nk
        do ib = 1, rt%nstate
            do jb = 1, rt%nstate
                jcur(:) = jcur(:) + (gs%kweight(ik) / gs%volume) * real( &
                    & gs%pmatrix(ib, jb, :, ik) * rt%rho(jb, ib, ik) &
                    & + gs%rvnl(ib, jb, :, ik) * rt%rho(jb, ib, ik))
            end do
            jcur(:) = jcur(:) + Ac(:) * (gs%kweight(ik) / gs%volume) * real(rt%rho(ib, ib, ik))
        end do
    end do
    !$omp end parallel do
    return
end subroutine
    

real(8) function calc_total(rt, gs)
    implicit none
    type(rt_data), intent(in) :: rt
    type(gs_data), intent(in) :: gs
    integer :: ik, ib
    real(8) :: tmp
    tmp = 0.0d0
    !$omp parallel do default(shared) private(ik,ib) reduction(+:tmp)
    do ik = 1, rt%nk
        do ib = 1, rt%nstate
            tmp = tmp + gs%kweight(ik) * real(rt%rho(ib, ib, ik))
        end do
    end do
    !$omp end parallel do
    calc_total = tmp
    return
end function

end module time_evolution