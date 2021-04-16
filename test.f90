module test
    implicit none
    contains

subroutine calc_dielec(sysname, base_directory, gs, nenergy, de, gamma)
    use sbe_solver
    use salmon_math, only: pi
    use salmon_file, only: open_filehandle
    implicit none
    character(*), intent(in) :: sysname
    character(*), intent(in) :: base_directory
    type(s_sbe_gs), intent(in) :: gs
    integer, intent(in) :: nenergy
    real(8), intent(in) :: de
    real(8), intent(in) :: gamma
    real(8) :: e
    complex(8) :: eps

    integer :: fh, ie, ik, ib, jb, i, j
    character(256) :: file_dielec_data

    do i = 1, 3
        do j = 1, 3
            write(file_dielec_data, '(a, a, a, i1, i1, a)') &
                trim(base_directory), trim(sysname), "_epsilon", i, j, ".data"
            
            fh = open_filehandle(file_dielec_data, 'replace')

            do ie = 1, nenergy
                e = de * ie
                if (i == j) then
                    eps = 1.0d0
                else
                    eps = 0.0d0
                end if
                do ik = 1, gs%nk
                    do ib = 1, gs%nb
                        do jb = 1, gs%nb
                            if (gs%occup(ib, ik) < 1.0 .and. 1.0 < gs%occup(jb, ik)) then
                                eps = eps + 2 * (4.0 * pi) / (gs%volume * gs%nk) &
                                & * gs%d_matrix(ib, jb, i, ik) &
                                & * gs%d_matrix(jb, ib, j, ik) &
                                & / (gs%omega(ib, jb, ik) - e - dcmplx(0d0, gamma)) 
                            end if
                        end do
                    end do
                end do
                write(fh, '(f9.4, 2(1x,e23.15e3))') e, real(eps), aimag(eps)
            end do
        end do
    end do
    close(fh)
end subroutine calc_dielec


end module test