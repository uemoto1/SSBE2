module epsilon
    implicit none

contains

    subroutine calc_epsilon(gs)
        use ground_state
        implicit none
        type(gs_data), intent(in) :: gs
        integer :: ik, ib, jb
        real(8) :: de
        complex(8), allocatable :: mu_tdm(:, :, :, :)
        character(99) :: file_epsilon
    
        real(8), parameter :: ax = 10.26
        real(8), parameter :: ay = 10.26
        real(8), parameter :: az = 10.26
        integer, parameter :: nomega = 1000
        real(8), parameter :: domega = 1d-3
        real(8), parameter :: gamma = 2.5d-3
        real(8), parameter :: vol = ax*ay*az
        integer :: i, j, iomega
        real(8) :: omega
        complex(8) ::eps
    
    
        allocate(mu_tdm(3, gs%nstate, gs%nstate, gs%nk))
        do ik = 1, gs%nk
            do ib = 1, gs%nstate
                do jb = 1, gs%nstate
                    de = gs%eigen(jb, ik) - gs%eigen(ib, ik)
                    if (abs(de) > 1d-6) then
                        if (gs%use_nonlocal_potential) then
                            mu_tdm(:, ib, jb, ik) = &
                                & dcmplx(0.0, 1.0) / de &
                                & * (gs%pmatrix(:, ib, jb, ik) &
                                & + gs%rvnl(:, ib, jb, ik))
                        else
                            mu_tdm(:, ib, jb, ik) = &
                                & dcmplx(0.0, 1.0) / de &
                                & * (gs%pmatrix(:, ib, jb, ik))                    
                        end if
                    else
                        mu_tdm(:, ib, jb, ik) = 0.0d0
                    end if
                end do
            end do
        end do
    
        do i = 1, 3
            do j = 1, 3
                write(file_epsilon, "(a,i1,i1,a)") "epsilon", i, j, ".txt"
                write(*, "(a,a)") "# Open:", trim(file_epsilon) 
                open(10, file=trim(file_epsilon), action="write")
                do iomega = 1, nomega
                    omega = domega * iomega
                    if (i == j) then
                        eps = 1.0d0
                    else
                        eps = 0.0d0
                    end if
                    do ik = 1, gs%nk
                        do ib = 1, gs%nstate
                            do jb = ib+1, gs%nstate
                                eps = eps + 4.0*3.14159/vol*gs%kweight(ik) &
                                    & * (gs%occup(ib, ik) - gs%occup(jb, ik)) &
                                    & * mu_tdm(i, ib, jb, ik) * mu_tdm(j, jb, ib, ik) &
                                    & / (gs%eigen(jb, ik) - gs%eigen(ib, ik) - omega - dcmplx(0.0, gamma))
                            end do
                        end do
                    end do
                    write(10, "(f12.6,2es25.15e4)") omega, real(eps), aimag(eps)
                end do        
                close(10)
            end do
        end do
        close(10)
    end subroutine

end module epsilon