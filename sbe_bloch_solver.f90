
module sbe_solver
    use salmon_math, only: pi
    use sbe_gs
    implicit none



    type s_sbe_solver
        !k-points for real-time SBE calculation
        integer :: nk, nb
        complex(8), allocatable :: rho(:, :, :)
        complex(8), allocatable :: pmatrix(:, :, :)
        real(8), allocatable :: domega(:, :, :)
        ! integer :: ik_s, ik_e, icomm_k
    end type



contains



subroutine init_sbe(sbe, gs, nkgrid)
    implicit none
    type(s_sbe_solver), intent(inout) :: sbe
    type(s_sbe_gs), intent(in) :: gs
    integer, intent(in) :: nkgrid(1:3)
    integer :: ik, ib
    integer :: nk, nb

    nk = nkgrid(1) * nkgrid(2) * nkgrid(3)
    nb = gs%nb

    sbe%nk = nk
    sbe%nb = nb

    allocate(sbe%rho(1:nb, 1:nb, 1:nk))
    allocate(sbe%pmatrix(1:nb, 1:nb, 1:nk))
    allocate(sbe%domega(1:nb, 1:nb, 1:nk))
    
    sbe%rho = 0d0

    do ik = 1, sbe%nk
        do ib = 1, sbe%nb
            sbe%rho(ib, ib, ik) = gs%occup(ib, ik)
        end do
    end do   


end subroutine








! subroutine calc_current(sbe, gs, Ac, jmat)
!     implicit none
!     type(s_sbe_solver), intent(in) :: sbe
!     type(s_sbe_gs), intent(in) :: gs
!     real(8), intent(in) :: Ac(1:3)
!     real(8), intent(out) :: jmat(1:3)
!     complex(8) :: pk(sbe%nb, sbe%nb, 3)
!     integer :: ik, idir, ib, jb
!     real(8) :: jtot(1:3)

!     jtot(1:3) = 0d0

!     !$omp parallel do default(shared) private(ik,ib,jb,idir,pk) reduction(+:jtot)
!     do ik=1, sbe%nk
!         call interp_gs(gs, sbe%kvec(1:3, ik) + Ac(1:3), p_k=pk)
!         do idir=1, 3
!             do jb=1, sbe%nb
!                 do ib=1, sbe%nb
!                     jtot(idir) = jtot(idir) + sbe%kweight(ik) &
!                         & * real(pk(ib, jb, idir) * sbe%rho(jb, ib, ik)) 
!                 end do
!             end do
!         end do
!     end do
!     !$omp end parallel do
!     jmat(:) =  (jtot(:) ) / gs%volume

!     return
! end subroutine calc_current








! subroutine calc_current_bloch(sbe, gs, Ac, jmat)
!     implicit none
!     type(s_sbe_solver), intent(in) :: sbe
!     type(s_sbe_gs), intent(in) :: gs
!     real(8), intent(in) :: Ac(1:3)
!     real(8), intent(out) :: jmat(1:3)
!     complex(8) :: pk(sbe%nb, sbe%nb, 3)
!     integer :: ik, idir, ib, jb
!     real(8) :: jtot(1:3)

!     jtot(1:3) = 0d0

!     !$omp parallel do default(shared) private(ik,ib,jb,idir,pk) reduction(+:jtot)
!     do ik=1, sbe%nk
!         call interp_gs(gs, sbe%kvec(1:3, ik) + Ac(1:3), p_k=pk)
!         do idir=1, 3
!             do jb=1, sbe%nb
!                 do ib=1, sbe%nb
!                     jtot(idir) = jtot(idir) + sbe%kweight(ik) &
!                         & * real(pk(ib, jb, idir) * sbe%rho(jb, ib, ik)) 
!                 end do
!             end do
!         end do
!     end do
!     !$omp end parallel do
!     jmat(:) =  (jtot(:) + gs%ne * ac(:)) / gs%volume
    

!     return
! end subroutine calc_current_bloch










! subroutine dt_evolve(sbe, gs, E, Ac, dt)
!     implicit none
!     type(s_sbe_solver), intent(inout) :: sbe
!     type(s_sbe_gs), intent(in) :: gs
!     real(8), intent(in) :: E(1:3)
!     real(8), intent(in) :: Ac(1:3)
!     real(8), intent(in) :: dt
!     complex(8), parameter :: zi = dcmplx(0d0, 1d0)

!     complex(8) :: hrho1(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)
!     complex(8) :: hrho2(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)
!     complex(8) :: hrho3(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)
!     complex(8) :: hrho4(1:sbe%nb, 1:sbe%nb, 1:sbe%nk)

!     call calc_hrho(sbe%rho, hrho1)
!     call calc_hrho(hrho1, hrho2)
!     call calc_hrho(hrho2, hrho3)
!     call calc_hrho(hrho3, hrho4)

!     sbe%rho = sbe%rho + hrho1 * (- zi * dt)
!     sbe%rho = sbe%rho + hrho2 * (- zi * dt) ** 2 * (1d0 / 2d0)
!     sbe%rho = sbe%rho + hrho3 * (- zi * dt) ** 3 * (1d0 / 6d0)
!     sbe%rho = sbe%rho + hrho4 * (- zi * dt) ** 4 * (1d0 / 24d0)
!     return

! contains

    ! !Calculate [H, rho] commutation:
    ! subroutine calc_hrho(rho, hrho)
    !     implicit none
    !     complex(8), intent(in)     :: rho(sbe%nb, sbe%nb, sbe%nk)
    !     complex(8), intent(out)    :: hrho(sbe%nb, sbe%nb, sbe%nk)
    !     real(8) :: ek(sbe%nb)
    !     complex(8) :: dk(sbe%nb, sbe%nb, 3)
    !     integer :: ik, idir, ib, jb
    !     !$omp parallel do default(shared) private(ik,ib,jb,idir,ek,dk)
    !     do ik=1, sbe%nk
    !         call interp_gs(gs, sbe%kvec(1:3, ik) + Ac(1:3),  e_k=ek, d_k=dk)
    !         !hrho(k) = omega(k + A/c) * rho(k) 
    !         do ib = 1, sbe%nb
    !             do jb = 1, sbe%nb
    !                 hrho(ib, jb, ik) = (ek(ib) - ek(jb)) * rho(ib, jb, ik)
    !             end do
    !         end do
    !         !hrho = hrho - E(t) * (d * rho - rho * d)
    !         do idir=1, 3 !1:x, 2:y, 3:z
    !             call ZHEMM('L', 'U', sbe%nb, sbe%nb, &
    !                 & dcmplx(-E(idir), 0d0), &
    !                 & dk(:, :, idir), sbe%nb, &
    !                 & rho(:, :, ik), sbe%nb, &
    !                 & dcmplx(1d0, 0d0), hrho(:, :, ik), sbe%nb)
    !             call ZHEMM('L', 'U', sbe%nb, sbe%nb, &
    !                 & dcmplx(+E(idir), 0d0), &
    !                 & rho(:, :, ik), sbe%nb, &
    !                 & dk(:, :, idir), sbe%nb, &
    !                 & dcmplx(1d0, 0d0), hrho(:, :, ik), sbe%nb)
    !         end do !idir
    !     end do !ik
    !     !$omp end parallel do
    !     return
    ! end subroutine calc_hrho


    !Calculate [H, rho] commutation:
!     subroutine calc_hrho_bloch(rho, hrho)
!         implicit none
!         complex(8), intent(in)     :: rho(sbe%nb, sbe%nb, sbe%nk)
!         complex(8), intent(out)    :: hrho(sbe%nb, sbe%nb, sbe%nk)
!         integer :: ik, ib, jb, idir
!         real(8) :: ek(sbe%nb)
!         complex(8) :: pk(sbe%nb, sbe%nb, 3)
!         !$omp parallel do default(shared) private(ik,ib,jb,idir,ek,pk)
!         do ik=1, sbe%nk
!             call interp_gs(gs, sbe%kvec(1:3, ik),  e_k=ek, p_k=pk)
!             !hrho(k) = omega(k + A/c) * rho(k) 
!             do ib = 1, sbe%nb
!                 do jb = 1, sbe%nb
!                     hrho(ib, jb, ik) = (ek(ib) - ek(jb)) * rho(ib, jb, ik)
!                 end do
!             end do
!             !hrho = hrho + Ac(t) * (p * rho - rho * p)
!             do idir=1, 3 !1:x, 2:y, 3:z
!                 call ZHEMM('L', 'U', sbe%nb, sbe%nb, &
!                     & dcmplx(+Ac(idir), 0d0), &
!                     & pk, sbe%nb, &
!                     & rho(:, :, ik), sbe%nb, &
!                     & dcmplx(1d0, 0d0), hrho(:, :, ik), sbe%nb)
!                 call ZHEMM('L', 'U', sbe%nb, sbe%nb, &
!                     & dcmplx(-Ac(idir), 0d0), &
!                     & rho(:, :, ik), sbe%nb, &
!                     & pk, sbe%nb, &
!                     & dcmplx(1d0, 0d0), hrho(:, :, ik), sbe%nb)
!             end do !idir
!         end do !ik
!         !$omp end parallel do
!         return
!     end subroutine calc_hrho_bloch
! end subroutine

! function calc_trace(sbe, nb_max) result(tr)
!     implicit none
!     type(s_sbe_solver), intent(in) :: sbe
!     integer, intent(in) :: nb_max
!     complex(8), parameter :: zi = dcmplx(0d0, 1d0)
!     integer :: ik, ib
!     real(8) :: tr
!     tr = 0d0
!     !$omp parallel do default(shared) private(ik, ib) reduction(+: tr) collapse(2) 
!     do ik = 1, sbe%nk
!         do ib = 1, nb_max
!             tr = tr + real(sbe%rho(ib, ib, ik)) * sbe%kweight(ik)
!         end do
!     end do
!     !$omp end parallel do
!     tr = tr / sum(sbe%kweight)
!     return 
! end function calc_trace


end module



