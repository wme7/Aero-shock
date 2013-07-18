module RK_Var
  implicit none
  !parameter for RK carpenter version
  real(kind=8), save             :: CFL
  real(kind=8), save             :: Time_Final
  real(kind=8), save, allocatable:: RK_para_a(:),RK_para_b(:),RK_para_c(:)
  real(kind=8), allocatable:: const_a_E(:,:),const_a_I(:,:)
  real(kind=8), allocatable:: const_b(:),const_b_hat(:),const_c(:)
  real(kind=8), save, allocatable:: alpha(:)
  real(kind=8), save             :: dt
  real(kind=8), save             :: time, time_loc 
  integer, save                  :: ks
  real(kind=8)                   :: dt_min,dt_local
  integer, save                  :: First_Comp, Last_Comp
  integer:: RK_Stage
  INTEGER, DIMENSION (1:2) :: ORDER2 = (/ 2, 1 /)
!  real(kind=8), allocatable:: A(:,:)
!  real(kind=8), allocatable:: c(:) ,b(:)

contains 

  !------------------------------------------------------------------------------

  subroutine Input_RK_parameter
    implicit none
    
    integer:: lid
    
    lid=20
    open(lid,file='RK_Para.in',form='formatted', status='old')
    read(lid,*) !==============================================!
    read(lid,*) CFL
    read(lid,*) Time_Final
    read(lid,*) !==============================================!
    
    close(lid) 
    
    return 
    
  end subroutine Input_RK_parameter
  
  !------------------------------------------------------------------------------
  !Calculation of the coefficients alpha for SSPRK
  subroutine Init_SSPRK(rk)
    implicit none
    integer:: rk,i,j,k,ierr
    real(kind=8)::temp

    allocate(alpha(1:rk), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate rk_alpha'
       stop
    endif

    alpha(1) = 1d0
    do i=1,rk
       do j=(i-1),1,-1
          temp = dble(j)
          temp = dble(1d0/temp)
          alpha(j+1) = temp*alpha(j)
       enddo
       temp = 1.d0
       do k=1,i
          temp = temp*dble(k)
       enddo
       alpha(i) = 1d0/dble(temp)
       alpha(1) = 1d0-sum(alpha(2:i))
    enddo
  end subroutine Init_SSPRK
    

!=====================================================================

  subroutine Init_RK4S5()
    implicit none
    
    integer:: RK_Stage
    integer i,ierr
    
    ! This is RK4S5 of 2N storage version
    RK_Stage=5
    allocate(Rk_para_a(0:RK_Stage), &
             Rk_para_b(0:RK_Stage), &
             Rk_para_c(0:RK_Stage), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate RK_para'
       write(*,*)'Abort!'
       stop
    endif
    !
    RK_para_a(1) =              0.0d0
    RK_para_a(2) =  -567301805773.0d0/ 1357537059087.0d0
    RK_para_a(3) = -2404267990393.0d0/ 2016746695238.0d0
    RK_para_a(4) = -3550918686646.0d0/ 2091501179385.0d0
    RK_para_a(5) = -1275806237668.0d0/  842570457699.0d0
    !
    RK_para_b(1) =  1432997174477.0d0/ 9575080441755.0d0
    RK_para_b(2) =  5161836677717.0d0/13612068292357.0d0
    RK_para_b(3) =  1720146321549.0d0/ 2090206949498.0d0
    RK_para_b(4) =  3134564353537.0d0/ 4481467310338.0d0
    RK_para_b(5) =  2277821191437.0d0/14882151754819.0d0
    !
    RK_para_c(1) =              0.0d0
    RK_para_c(2) =  1432997174477.0d0/ 9575080441755.0d0
    RK_para_c(3) =  2526269341429.0d0/ 6820363962896.0d0
    RK_para_c(4) =  2006345519317.0d0/ 3224310063776.0d0
    RK_para_c(5) =  2802321613138.0d0/ 2924317926251.0d0

  end subroutine Init_RK4S5
  !------------------------------------------------------------------------------

  !Calculation of the coefficients  for ARK
  subroutine Init_ARK()
    implicit none
    real(kind=8) :: ARK_tmp(6*6)
!    integer:: RK_Stage
    integer i,ierr

    ! This is RK4S5 of 2N storage version
    RK_Stage=6
    allocate(const_a_E(RK_Stage,RK_Stage), &
             const_a_I(RK_Stage,RK_Stage), &
             const_b(RK_Stage),const_b_hat(RK_Stage), &
             const_c(RK_Stage), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate coefficients for ARK'
       write(*,*)'Abort!'
       stop
    endif
!============ coeffice of ARK4_6 ============
    ARK_tmp = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
      1d0/2d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
      13861d0/62500d0, 6889d0/62500d0, 0d0, 0d0, 0d0, 0d0, &
      -116923316275d0/2393684061468d0, -2731218467317d0/15368042101831d0, &
        9408046702089d0/11113171139209d0, 0d0, 0d0, 0d0, &
      -451086348788d0/2902428689909d0,  -2682348792572d0/7519795681897d0, &
        12662868775082d0/11960479115383d0, 3355817975965d0/11060851509271d0,0d0, 0d0, &
      647845179188d0/3216320057751d0,     73281519250d0/8382639484533d0,  &
        552539513391d0/3454668386233d0,  3354512671639d0/8306763924573d0, 4040d0/17871d0, 0d0 /)
    const_a_E = RESHAPE( ARK_tmp, (/ 6, 6 /), ORDER = ORDER2)

    ARK_tmp = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
      1d0/4d0, 1d0/4d0, 0d0, 0d0, 0d0, 0d0, &
      8611d0/62500d0, -1743d0/31250d0, 1d0/4d0, 0d0, 0d0, 0d0, &
      5012029d0/34652500d0, -654441d0/2922500d0, 174375d0/388108d0, 1d0/4d0, 0d0, 0d0, &
      15267082809d0/155376265600d0, -71443401d0/120774400d0, &
        730878875d0/902184768d0, 2285395d0/8070912d0, 1d0/4d0,   0d0, &
      82889d0/524892d0, 0d0, 15625d0/83664d0, 69875d0/102672d0, -2260d0/8211d0, 1d0/4d0 /)
    const_a_I = RESHAPE( ARK_tmp, (/ 6, 6 /), ORDER = ORDER2)

    const_b = (/ 82889d0/524892d0, 0d0, 15625d0/83664d0, 69875d0/102672d0, -2260d0/8211d0, &
      1d0/4d0  /)

    const_b_hat = (/ 4586570599d0/29645900160d0, 0d0, 178811875d0/945068544d0, &
      814220225d0/1159782912d0, -3700637d0/11593932d0, 61727d0/225920d0 /)

    const_c = (/ 0d0, 1d0/2d0, 83d0/250d0, 31d0/50d0, 17d0/20d0, 1d0 /)

  end subroutine Init_ARK

  
  
end module RK_VAR
