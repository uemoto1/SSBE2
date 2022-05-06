module em_field
  implicit none

contains

!===================================================================================================================================

Subroutine calc_Ac_ext(t,Ac_ext)
  real(8),intent(in) :: t
  real(8)            :: Ac_ext(1:3)
  real(8)            :: Ac_ext_t(1:3, 0:0)

  call calc_Ac_ext_t(t, 0d0, 0, 0, Ac_ext_t)
  Ac_ext(1:3) = Ac_ext_t(1:3, 0)
  return
end subroutine calc_Ac_ext

!===================================================================================================================================

Subroutine calc_Ac_ext_t(t0, delta_t, is, ie, Ac_ext_t)
  use math_constants,only: zi,pi
  ! use salmon_global, only: I_wcm2_1,I_wcm2_2,E_amplitude1,E_amplitude2,ae_shape1,ae_shape2, &
  !                        & epdir_re1,epdir_re2,epdir_im1,epdir_im2,tw1,tw2,t1_start,omega1,omega2, &
  !                        & phi_CEP1,phi_CEP2,T1_T2,e_impulse,file_input1
  use input_parameter, only: I_wcm2_1,I_wcm2_2,E_amplitude1,E_amplitude2,ae_shape1,ae_shape2, &
  & epdir_re1,epdir_re2,epdir_im1,epdir_im2,tw1,tw2,t1_start,omega1,omega2, &
  & phi_CEP1,phi_CEP2,T1_T2,e_impulse!,file_input1

  implicit none
  real(8),intent(in) :: t0
  real(8),intent(in) :: delta_t
  integer,intent(in) :: is
  integer,intent(in) :: ie
  real(8),intent(out) :: Ac_ext_t(1:3, is:ie)
  !
  integer :: i,npower
  real(8) :: t,f0_1,f0_2,tt,T1_T2_tmp
  
  Ac_ext_t(:,:) = 0d0
!  if(t0 < 0d0) return

  if(I_wcm2_1 < 0d0)then
    f0_1 = E_amplitude1
  else
    f0_1=5.338d-9*sqrt(I_wcm2_1)      ! electric field in a.u.
  end if
  if(I_wcm2_2 < 0d0)then
    f0_2 = E_amplitude2
  else
    f0_2=5.338d-9*sqrt(I_wcm2_2)      ! electric field in a.u.
  end if
  
  T1_T2_tmp = T1_T2

  select case(AE_shape1)
  case('impulse')
    do i=is, ie
      t=t0+i*delta_t
      if (0d0 <= t) then
        Ac_ext_t(1,i) = epdir_re1(1)*e_impulse
        Ac_ext_t(2,i) = epdir_re1(2)*e_impulse
        Ac_ext_t(3,i) = epdir_re1(3)*e_impulse
      end if
    end do
    
  case('Acos2','Acos3','Acos4','Acos6','Acos8')
  
    select case(ae_shape1)
    case('Acos2'); npower = 2
    case('Acos3'); npower = 3
    case('Acos4'); npower = 4
    case('Acos6'); npower = 6
    case('Acos8'); npower = 8
    case default
      stop 'Error in init_rt.f90'
    end select

    do i=is, ie
      t=t0+i*delta_t
      tt = t - 0.5d0*tw1 - t1_start
      if (abs(tt)<0.5d0*tw1) then
        Ac_ext_t(:,i) = -f0_1/omega1*(cos(pi*tt/tw1))**npower &
          *aimag( (epdir_re1(:) + zI*epdir_im1(:)) &
          *exp(zI*(omega1*tt+phi_CEP1*2d0*pi))  &
          )
      end if
    end do
    T1_T2_tmp = T1_T2 + t1_start

  case('Ecos2')
  
    if(phi_CEP1 /= 0.75d0)then
      stop "Error: phi_cep1 should be 0.75 when ae_shape1 is 'Ecos2'."
    end if
    if(sum(abs(epdir_im1(:)))>1.0d-8)then
      stop "Error: ae_shape1 should be 'Acos2' when epdir_im1 is used."
    end if
    do i=is, ie
      t=t0+i*delta_t
      tt = t - 0.5d0*tw1 - t1_start
      if (abs(tt)<0.5d0*tw1) then
        Ac_ext_t(:,i) = -epdir_re1(:)*f0_1/(8d0*pi**2*omega1 - 2d0*tw1**2*omega1**3) &
          *( &
          (-4d0*pi**2+tw1**2*omega1**2 + tw1**2*omega1**2*cos(2d0*pi*tt/tw1))*cos(omega1*tt) &
          +2d0*pi*(2d0*pi*cos(tw1*omega1/2d0) &
          +tw1*omega1*sin(2d0*pi*tt/tw1)*sin(omega1*tt)))
      end if
    end do
    T1_T2_tmp = T1_T2 + t1_start

  case('Esin2sin')
  
    stop "Esin2sin is not implemented"
    
  case('Asin2cos')
  
    ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
    ! pump laser
    do i=is, ie
      t=t0+i*delta_t
      tt = t
      if (tt<tw1) then
        Ac_ext_t(:,i) = -epdir_re1(:)*f0_1/omega1*(sin(pi*tt/tw1))**2*cos(omega1*tt+phi_CEP1*2d0*pi)
      end if
    end do
    
  ! case('input')
     
  !   call add_Ac_from_file(trim(file_input1), t0, delta_t, is, ie, Ac_ext_t)
  !   !  stop "ae_shape1='input' is not implemented"
    
  case('Asin2_cw')
  
    ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
    ! pump laser
    do i=is, ie
      t=t0+i*delta_t
      tt = t
      if (tt<tw1*0.5d0) then
        Ac_ext_t(:,i) = -Epdir_re1(:)*f0_1/omega1*(sin(pi*tt/tw1))**2*cos(omega1*tt+phi_CEP1*2d0*pi)
      else
        Ac_ext_t(:,i) = -Epdir_re1(:)*f0_1/omega1*cos(omega1*tt+phi_CEP1*2d0*pi)
      end if
    end do
    
  case('none')
  
    Ac_ext_t(:,:) = 0d0
    
  case default
  
    stop "Invalid pulse_shape_1 parameter!"
    
  end select

! Probe
  select case(ae_shape2)
  case('impulse')

    do i=is, ie
      t=t0+i*delta_t
      tt = t - 0.5d0*tw1 - T1_T2_tmp
      if(tt > 0d0)then
        Ac_ext_t(1,i) = Ac_ext_t(1,i) + epdir_re2(1)*e_impulse
        Ac_ext_t(2,i) = Ac_ext_t(2,i) + epdir_re2(2)*e_impulse
        Ac_ext_t(3,i) = Ac_ext_t(3,i) + epdir_re2(3)*e_impulse
      end if
    end do
    
  case('Acos2','Acos3','Acos4','Acos6','Acos8')
  
    select case(ae_shape2)
    case('Acos2'); npower = 2
    case('Acos3'); npower = 3
    case('Acos4'); npower = 4
    case('Acos6'); npower = 6
    case('Acos8'); npower = 8
    case default
      stop 'Error in init_Ac.f90'
    end select

    do i=is, ie
      t=t0+i*delta_t
      tt = t - 0.5d0*tw1 - T1_T2_tmp
      if (abs(tt)<0.5d0*tw2) then
        Ac_ext_t(:,i)=Ac_ext_t(:,i) &
          -f0_2/omega2*(cos(pi*tt/tw2))**npower &
          *aimag( (epdir_re2(:) + zI*epdir_im2(:)) &
          *exp(zI*(omega2*tt+phi_CEP2*2d0*pi))  &
          )
      end if
    end do

  case('Ecos2')
  
    if(phi_CEP2 /= 0.75d0)then
      stop "Error: phi_cep2 should be 0.75 when ae_shape2 is 'Ecos2'."
    end if
    if(sum(abs(epdir_im2(:)))>1.0d-8)then
      stop "Error: ae_shape2 should be 'Acos2' when epdir_im2 is used."
    end if

    do i=is, ie
      t=t0+i*delta_t
      tt = t - 0.5d0*tw1 - T1_T2_tmp
      if (abs(tt)<0.5d0*tw2) then
        Ac_ext_t(:,i)=Ac_ext_t(:,i) &
          -epdir_re2(:)*f0_2/(8d0*pi**2*omega2 - 2d0*tw2**2*omega2**3) &
          *( &
          (-4d0*pi**2+tw2**2*omega2**2 + tw2**2*omega2**2*cos(2d0*pi*tt/tw2))*cos(omega2*tt) &
          +2d0*pi*(2d0*pi*cos(tw2*omega2/2d0) &
          +tw2*omega2*sin(2d0*pi*tt/tw2)*sin(omega2*tt)))
      end if
    end do

  case('Esin2sin')
      
    stop "Esin2sin is not implemented"
    
  case('Asin2cos')
  
      ! pulse shape : A(t)=f0/omega*sin(Pi t/T)**2 *cos (omega t+phi_CEP*2d0*pi) 
    ! probe laser

    do i=is, ie
      t=t0+i*delta_t
      tt = t
      if ( (tt-T1_T2_tmp>0d0) .and. (tt-T1_T2_tmp<tw2) ) then
        Ac_ext_t(:,i) = Ac_ext_t(:,i) &
          &-Epdir_re2(:)*f0_2/omega2*(sin(pi*(tt-T1_T2_tmp)/tw2))**2*cos(omega2*(tt-T1_T2_tmp)+phi_CEP2*2d0*pi)
      endif
    end do

  case('input')
  case('Asin2_cw')
  case('none')
  case default
  
    stop "Invalid pulse_shape_2 parameter!"
    
  end select

  return
End Subroutine calc_Ac_ext_t

end module em_field
!--------10--------20--------30--------40--------50--------60--------70--------80--------90--------100-------110-------120-------130
