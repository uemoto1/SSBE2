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
    complex(8) :: rho_k_tmp(rt%nstate, rt%nstate)
    complex(8) :: drho1_k(rt%nstate, rt%nstate)
    complex(8) :: drho2_k(rt%nstate, rt%nstate)
    complex(8) :: drho3_k(rt%nstate, rt%nstate)
    complex(8) :: drho4_k(rt%nstate, rt%nstate)
    integer :: ik

    !$omp parallel do default(shared) private(ik,rho_k_tmp,drho1_k,drho2_k,drho3_k,drho4_k)
    do ik = 1, rt%nk
        ! 4th-order Runge-Kutta method
        call calc_drho_k(rt%nstate, drho1_k, rt%rho(:, :, ik), &
            & gs%omega(:, :, ik), gs%pmatrix(:, :, :, ik), Ac0)
        rho_k_tmp = rt%rho(:, :, ik) + 0.5 * dt * drho1_k
        call calc_drho_k(rt%nstate, drho2_k, rho_k_tmp, &
            & gs%omega(:, :, ik), gs%pmatrix(:, :, :, ik), 0.5*(Ac0+Ac1))
        rho_k_tmp = rt%rho(:, :, ik) + 0.5 * dt * drho2_k
        call calc_drho_k(rt%nstate, drho3_k, rho_k_tmp, &
            & gs%omega(:, :, ik), gs%pmatrix(:, :, :, ik), 0.5*(Ac0+Ac1))
        rho_k_tmp = rt%rho(:, :, ik) + dt * drho3_k
        call calc_drho_k(rt%nstate, drho4_k, rho_k_tmp, &
            & gs%omega(:, :, ik), gs%pmatrix(:, :, :, ik), Ac1)
        rt%rho(:, :, ik) = rt%rho(:, :, ik) &
            & + (dt/6) * (drho1_k(:, :) + 2*drho2_k(:, :) + 2*drho3_k(:, :) + drho4_k(:, :))
    end do
    !$omp end parallel do

contains

subroutine calc_drho_k(nstate, drho_k, rho_k, omega_k, p_k, Ac)
    implicit none
    integer, intent(in) :: nstate
    complex(8), intent(out) :: drho_k(nstate, nstate)
    complex(8), intent(in) :: rho_k(nstate, nstate)
    real(8), intent(in) :: omega_k(nstate, nstate)
    complex(8), intent(in) :: p_k(nstate, nstate, 3)
    real(8), intent(in) :: Ac(3)
    integer :: i

    drho_k(:, :) = dcmplx(0.0, -1.0) * omega_k(:, :) * rho_k(:, :)
    do i = 1, 3
        drho_k(:, :) = drho_k(:, :) + dcmplx(0.0, -Ac(i)) * ( &
            & matmul(p_k(:, :, i), rho_k(:, :)) - matmul(rho_k(:, :), p_k(:, :, i)))
    end do
end subroutine calc_drho_k
end subroutine dt_evolve_bloch

subroutine current(jcur, qtot, rt, gs, Ac)
    implicit none
    real(8), intent(out) :: jcur(3)
    real(8), intent(out) :: qtot
    type(rt_data), intent(in) :: rt
    type(gs_data), intent(in) :: gs
    real(8), intent(in) :: Ac(3)
    integer :: ik, ib, jb

    jcur(:) = 0.0d0
    qtot = 0.0d0
    !$omp parallel do default(shared) private(ik,ib,jb) reduction(+:qtot,jcur) collapse(2)
    do ik = 1, rt%nk
        do ib = 1, rt%nstate
            qtot = qtot + gs%kweight(ik) * real(rt%rho(ib, ib, ik))
            do jb = 1, rt%nstate
                jcur(:) = jcur(:) + (gs%kweight(ik) / gs%volume) * real( &
                    & gs%pmatrix(ib, jb, :, ik) * rt%rho(jb, ib, ik))
            end do
        end do
    end do
    !$omp end parallel do
    jcur(:) = jcur(:) + Ac(:) * qtot / gs%volume
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
            tmp = tmp + gs%kweight(ik) * gs%eigen(ib, ik) * real(rt%rho(ib, ib, ik))
        end do
    end do
    !$omp end parallel do
    calc_total = tmp
    return
end function

end module time_evolution