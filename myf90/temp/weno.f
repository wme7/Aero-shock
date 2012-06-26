!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!  Fifth-order weno reconstruction 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!    subroutine wenoReconst( tu, tuL,tuR ) 
!    use mainvar 
!    implicit none 
!    real :: tu(0:NN+1), tuL(0:NN), tuR(0:NN) 
!    real :: PC(3,0:NN), LCoef(3), d(3), beta1,beta2,beta3,alp1,alp2,alp3,salp,we1,we2,we3 
!    integer :: i 
!     
!    d=(/0.1,0.6,0.3/) 
!    !! initial polynomial 
!    do i=1, NN   !! p= PCoef(1)+PCoef(2)*(x-xc)/dx*2+PCoef(3)*(()**2-1/3) 
!      PC(1,i)= tu(i) 
!      PC(2,i)= (tu(i+1)-tu(i-1))/4. 
!      PC(3,i)= (tu(i+1)+tu(i-1)-2*tu(i))/8. 
!    enddo 
!    !! smooth 
!    do i=2,NN-1 
!      beta1= 13./12*( (tu(i-2)-2*tu(i-1)+tu(i))**2 )+ 1./4*( (tu(i-2)-4*tu(i-1)+3*tu(i))**2 ) 
!      beta2= 13./12*( (tu(i-1)-2*tu(i)+tu(i+1))**2 )+ 1./4*( (tu(i-1)-tu(i+1))**2 ) 
!      beta3= 13./12*( (tu(i)-2*tu(i+1)+tu(i+2))**2 )+ 1./4*( (3*tu(i)-4*tu(i+1)+tu(i+2))**2 ) 
!      alp1= d(1)/(EPSILON+beta1)**2   !!????? 
!      alp2= d(2)/(EPSILON+beta2)**2 
!      alp3= d(3)/(EPSILON+beta3)**2 
!      salp= alp1+alp2+alp3 
!      we1 = alp1/salp 
!      we2 = alp2/salp 
!      we3 = alp3/salp 
!      LCoef(1)= we1*(PC(1,i-1)+2*PC(2,i-1)+4*PC(3,i-1))+we2*PC(1,i) +we3*(PC(1,i+1)-2*PC(2,i+1)+4*PC(3,i+1)) 
!      LCoef(2)= we1*(PC(2,i-1)+4*PC(3,i-1))            +we2*PC(2,i) +we3*(PC(2,i+1)-4*PC(3,i+1)) 
!      LCoef(3)= we1* PC(3,i-1)                         +we2*PC(3,i) +we3* PC(3,i+1) 
!      tuL(i)  = LCoef(1)+ LCoef(2)+ 2./3*LCoef(3) 
! 
! 
!      beta1= 13./12*( (tu(i-2)-2*tu(i-1)+tu(i))**2 )+ 1./4*( (tu(i-2)-4*tu(i-1)+3*tu(i))**2 ) 
!      beta2= 13./12*( (tu(i-1)-2*tu(i)+tu(i+1))**2 )+ 1./4*( (tu(i-1)-tu(i+1))**2 ) 
!      beta3= 13./12*( (tu(i)-2*tu(i+1)+tu(i+2))**2 )+ 1./4*( (3*tu(i)-4*tu(i+1)+tu(i+2))**2 ) 
!      alp1= d(3)/(EPSILON+beta1)**2     !! ???? 
!      alp2= d(2)/(EPSILON+beta2)**2 
!      alp3= d(1)/(EPSILON+beta3)**2 
!      salp= alp1+alp2+alp3 
!      we1 = alp1/salp 
!      we2 = alp2/salp 
!      we3 = alp3/salp 
!      LCoef(1)= we1*(PC(1,i-1)+2*PC(2,i-1)+4*PC(3,i-1))+we2*PC(1,i) +we3*(PC(1,i+1)-2*PC(2,i+1)+4*PC(3,i+1)) 
!      LCoef(2)= we1*(PC(2,i-1)+4*PC(3,i-1))            +we2*PC(2,i) +we3*(PC(2,i+1)-4*PC(3,i+1)) 
!      LCoef(3)= we1* PC(3,i-1)                         +we2*PC(3,i) +we3* PC(3,i+1) 
!      tuR(i-1)= LCoef(1)- LCoef(2)+ 2./3*LCoef(3) 
!    enddo 
! 
!	tuL(1)= tu(1) 
!	tuR(0)= tu(1) 
!	tuL(NN)  = tu(NN) 
!	tuR(NN-1)= tu(NN) 
! 
!    end subroutine wenoReconst 
! 
! 
!!!!!!  use left and right state to demonstrate !!!!!! 
!    subroutine wenoReconst2( tv, tvL,tvR ) 
!    use mainvar 
!    implicit none 
!    real,parameter :: eps=1.0e-6 
!    real :: vr0,vr1,vr2 
!    real :: tv(0:NN+1), tvL(0:NN), tvR(0:NN) 
!    real :: beta1,beta2,beta0,omega0,omega1,omega2,omega 
!    integer :: i 
! 
!    !! use eigen varialbes to reconstr 
! 
! 
!    do i=2,NN-1 
!      !! 单元内左 
!      Beta0 = 13.0*(tv(i) - 2.0*tv(i+1) + tv(i+2) )**2.0 /12.0+& 
!		        0.25*(3.0*tv(i) - 4.0*tv(i+1) + tv(i+2) )**2.0 
!	  Beta1 = 13.0*(tv(i-1) -2.0*tv(i) + tv(i+1) )**2.0 /12.0 +& 
!		        0.25*(tv(i+1) -    tv(i-1)           )**2.0 
!	  Beta2 = 13.0*(tv(i-2) -2.0*tv(i-1) + tv(i) )**2.0 /12.0 +& 
!		        0.25*(tv(i-2) -4.0*tv(i-1) +3.0*tv(i) )**2.0 
!      omega0 = 0.1 /(eps+Beta0)**2.0 
!      omega1 = 0.6 /(eps+Beta1)**2.0 
!      omega2 = 0.3 /(eps+Beta2)**2.0 
!      omega  = omega0 + omega1 + omega2 
! 
!      omega0 = omega0/omega 
!      omega1 = omega1/omega 
!      omega2 = omega2/omega 
! 
!      vR0 = 11.0*tv( i)/6.0 - 7.0*tv(i+1) / 6.0 + tv(i+2) / 3.0 
!	  vR1 =      tv( i-1)/3.0 + 5.0*tv(i) / 6.0 - tv(i+1) / 6.0 
!	  vR2 =     -tv(i-2)/6.0 + 5.0*tv(i-1) / 6.0 + tv(i) / 3.0 
! 
!       
!	  tvR(i-1)  = omega0 * vR0 + omega1 * vR1 + omega2 * vR2 
!     enddo 
! 
!	 do i=2,NN-1 
! 
!      !! 单元内右 
!      Beta0 = 13.0*(tv(i) -2.0*tv(i+1) + tv(i+2) )**2.0 /12.0 +& 
!		        0.25*(3.0*tv(i) - 4.0*tv(i+1) + tv(i+2) )**2.0 
!      Beta1 = 13.0*(tv(i-1) -2.0*tv(i) + tv(i+1) )**2.0 /12.0 +& 
!		        0.25*(tv(i-1) - tv(i+1) )**2.0 
!      Beta2 = 13.0*(tv(i-2) -2.0*tv(i-1) + tv(i) )**2.0 /12.0 +& 
!		        0.25*(tv(i-2) - 4.0*tv(i-1) + 3.0*tv(i) )**2.0 
!      omega0 = 0.3 /(eps+Beta0)**2.0 
!      omega1 = 0.6 /(eps+Beta1)**2.0 
!      omega2 = 0.1 /(eps+Beta2)**2.0 
!      omega  = omega0 + omega1 + omega2 
! 
!      omega0 = omega0/omega 
!      omega1 = omega1/omega 
!      omega2 = omega2/omega 
! 
!      vR0 = tv(i )/3.0 + 5.0*tv( i+1) / 6.0 -     tv(i+2) / 6.0 
!      vR1 =-tv(i-1)/6.0 + 5.0*tv( i) / 6.0 +     tv(i+1) / 3.0 
!      vR2 = tv(i-2)/3.0 - 7.0*tv(i-1) / 6.0 +11.0*tv(i) / 6.0 
!      tvL(i)  = omega0 * vR0 + omega1 * vR1 + omega2 * vR2 
!    enddo 
! 
!    tvL(1)= tv(1) 
!	tvR(0)= tv(1) 
!	tvL(NN)  = tv(NN) 
!	tvR(NN-1)= tv(NN) 
!    end subroutine wenoReconst2 
 
 
!!!!!  use left and right state to demonstrate !!!!!! 
    subroutine wenoReconst( tv, tvL,tvR ) 
    use mainvar 
    implicit none 
    real,parameter :: eps=1.0e-6 
    real :: vr0,vr1,vr2 
    real :: tv(6), tvL, tvR 
    real :: beta1,beta2,beta0,omega0,omega1,omega2,omega 
    integer :: i 
 
    !! (i+1/2)+ 
      Beta0 = 13.0*(tv(4) - 2.0*tv(5) + tv(6) )**2.0 /12.0+& 
		        0.25*(3.0*tv(4) - 4.0*tv(5) + tv(6) )**2.0 
	  Beta1 = 13.0*(tv(3) -2.0*tv(4) + tv(5) )**2.0 /12.0 +& 
		        0.25*(tv(5) -    tv(3)           )**2.0 
	  Beta2 = 13.0*(tv(2) -2.0*tv(3) + tv(4) )**2.0 /12.0 +& 
		        0.25*(tv(2) -4.0*tv(3) +3.0*tv(4) )**2.0 
      omega0 = 0.1 /(eps+Beta0)**2.0 
      omega1 = 0.6 /(eps+Beta1)**2.0 
      omega2 = 0.3 /(eps+Beta2)**2.0 
      omega  = omega0 + omega1 + omega2 
 
      omega0 = omega0/omega 
      omega1 = omega1/omega 
      omega2 = omega2/omega 
 
      vR0 = 11.0*tv( 4)/6.0 - 7.0*tv(5) / 6.0 + tv(6) / 3.0 
	  vR1 =      tv( 3)/3.0 + 5.0*tv(4) / 6.0 - tv(5) / 6.0 
	  vR2 =     -tv(2)/6.0 + 5.0*tv(3) / 6.0 + tv(4) / 3.0 
 
	  tvR  = omega0 * vR0 + omega1 * vR1 + omega2 * vR2 
 
 
    !! (i+1/2)- 
      Beta0 = 13.0*(tv(3) -2.0*tv(4) + tv(5) )**2.0 /12.0 +& 
		        0.25*(3.0*tv(3) - 4.0*tv(4) + tv(5) )**2.0 
      Beta1 = 13.0*(tv(2) -2.0*tv(3) + tv(4) )**2.0 /12.0 +& 
		        0.25*(tv(2) - tv(4) )**2.0 
      Beta2 = 13.0*(tv(1) -2.0*tv(2) + tv(3) )**2.0 /12.0 +& 
		        0.25*(tv(1) - 4.0*tv(2) + 3.0*tv(3) )**2.0 
      omega0 = 0.3 /(eps+Beta0)**2.0 
      omega1 = 0.6 /(eps+Beta1)**2.0 
      omega2 = 0.1 /(eps+Beta2)**2.0 
      omega  = omega0 + omega1 + omega2 
 
      omega0 = omega0/omega 
      omega1 = omega1/omega 
      omega2 = omega2/omega 
 
      vR0 = tv(3 )/3.0 + 5.0*tv( 4) / 6.0 -     tv(5) / 6.0 
      vR1 =-tv(2)/6.0 + 5.0*tv( 3) / 6.0 +     tv(4) / 3.0 
      vR2 = tv(1)/3.0 - 7.0*tv(2) / 6.0 +11.0*tv(3) / 6.0 
      tvL  = omega0 * vR0 + omega1 * vR1 + omega2 * vR2 
 
    end subroutine wenoReconst