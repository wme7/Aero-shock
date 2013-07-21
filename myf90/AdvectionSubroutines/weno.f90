!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  Fifth-order weno reconstruction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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