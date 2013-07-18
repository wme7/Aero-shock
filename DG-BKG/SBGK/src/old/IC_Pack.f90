subroutine Initial_Field(demo_case)
use State_Var
use MD2D_Grid
implicit none
integer::demo_case

    select case(demo_case)
    case(0) 
       ! do nothing
       
    case(1) 
       call demo_001

    end select


return 
end subroutine Initial_Field

subroutine demo_001
  use State_Var
  use MD2D_Grid
  use Material_Var
  implicit none
  integer::i,j
  integer(kind=8):: omega
  real(kind=8)::x_coor, y_coor, t_coor
  real(kind=8)::EL_demo001
  
  
  t_coor=0.d0
  do DDK=1,TotNum_DM 
     ND1=PolyDegN_DM(1,DDK)
     ND2=PolyDegN_DM(2,DDK)
     
     do j=0,ND2
        do i=0,ND1 
           x_coor=x1(i,j,DDK)
           y_coor=x2(i,j,DDK)
           
           
           v1(i,j,DDK) = &
                EL_demo001(rho(DDK),Lame_lambda(DDK),Lame_mu(DDK), &
                x_coor,y_coor,t_coor,1)
           
           v2(i,j,DDK) = &
                EL_demo001(rho(DDK),Lame_lambda(DDK),Lame_mu(DDK), &
                x_coor,y_coor,t_coor,2)
           
           T11(i,j,DDK) = &
                EL_demo001(rho(DDK),Lame_lambda(DDK),Lame_mu(DDK), &
                x_coor,y_coor,t_coor,3)
           
           T12(i,j,DDK) = &
                EL_demo001(rho(DDK),Lame_lambda(DDK),Lame_mu(DDK), &
                x_coor,y_coor,t_coor,4)
           
           T22(i,j,DDK) = & 
                EL_demo001(rho(DDK),Lame_lambda(DDK),Lame_mu(DDK), &
                x_coor,y_coor,t_coor,5)

!             fs1(i,j,DDK)=EM_demo001(x_coor,y_coor,z_coor,0.d0,6)
!             fs2(i,j,DDK)=EM_demo001(x_coor,y_coor,z_coor,0.d0,7)
           
        enddo
     enddo
  enddo



end subroutine demo_001

function EL_demo001(rho,lambda,mu,x,y,t,EL001_id)

  !=============================================================
  ! 
  ! 
  !=============================================================

  implicit none
  real(kind=8):: EL_demo001
  real(kind=8):: rho, lambda, mu, omega, x, y, t
  integer:: EL001_ID
  
  real(kind=8), parameter:: pi=3.14159265358979d0
  real(kind=8):: omega


  ! subroutine begins
  !
  if ( (EL001_ID .lt. 1) .or. (EL001_ID .gt. 7))  then 
     write(*,*)'IC_Pack.f90'
     write(*,*)'em_id is out of range'
     stop
  endif

  omega=2.d0*pi

  select case(EL001_ID)
  case(1) ! v1 field

     EL_demo001 = dsin(pi*x)*dsin(pi*y)*dcos(omega*t)

  case(2) ! v2 field

     EL_demo001 = dsin(pi*x)*dsin(pi*y)*dcos(omega*t)

  case(3) ! T11 field

     EL_demo001 = pi/omega*dsin(omega*t)* &
          (lambda*dsin(pi*(x+y))+2.d0*mu*dcos(pi*x)*dsin(pi*y))

  case(4) ! T12 Field
     
     EL_demo001 = mu*pi/omega*dsin(omega*t)* dsin(pi*(x+y))
          

  case(5)! T22 field

     EL_demo001 = pi/omega*dsin(omega*t)* &
          (lambda*dsin(pi*(x+y))+2.d0*mu*dsin(pi*x)*dcos(pi*y))

  case(6) !f1 field

     EL_demo001 = dsin(omega*t) * (&
          (2.d0*pi*pi*mu/omega-rho*omega)*dsin(pi*x)*dsin(pi*y) -&
           pi*pi/omega*(lambda+mu)*dcos(pi*(x+y)) )

  case(7) !f2 field

     EL_demo001 = dsin(omega*t) * (&
          (2.d0*pi*pi*mu/omega-rho*omega)*dsin(pi*x)*dsin(pi*y) -&
           pi*pi/omega*(lambda+mu)*dcos(pi*(x+y)) )

  end select

 
end function EL_demo001
