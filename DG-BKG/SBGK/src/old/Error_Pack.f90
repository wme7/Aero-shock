subroutine Error_001
use MD2D_Grid
use State_Var
use RK_Var
use Material_Var
implicit none 
integer i,j,k 
real(kind=8):: x_coor, y_coor, z_coor
real(kind=8):: err_loc(1:5), err_max(1:5)
real(kind=8):: EL_demo001
real(kind=8):: error
integer :: lid


     lid=85
     if (First_Comp .eq. 1) then
        open(lid,file='Error001.dat')
        write(lid,*)'VARIABLES = "Time", "v1", "v2", "T11", "T12", "T22"'
     endif


     err_max(1:5)=0.d0
     do DDK=1,TotNum_DM
        ND1=PolyDegN_DM(1,DDK)
        ND2=PolyDegN_DM(2,DDK)

        
        do j=0,ND2 
           do i=0,ND1
              
              x_coor=x1(i,j,DDK)
              y_coor=x2(i,j,DDK)
              
              err_loc(1)= abs( v1(i,j,DDK) - &
                EL_demo001(rho(DDK),Lame_lambda(DDK),Lame_mu(DDK),&
                           x_coor,y_coor,time,1))

              err_loc(2)= abs( v2(i,j,DDK) - &
                EL_demo001(rho(DDK),Lame_lambda(DDK),Lame_mu(DDK),&
                           x_coor,y_coor,time,2))

              err_loc(3)= abs(T11(i,j,DDK) - &
                EL_demo001(rho(DDK),Lame_lambda(DDK),Lame_mu(DDK),&
                           x_coor,y_coor,time,3))

              err_loc(4)= abs(T12(i,j,DDK) - &
                EL_demo001(rho(DDK),Lame_lambda(DDK),Lame_mu(DDK),&
                           x_coor,y_coor,time,4))

              err_loc(5)= abs(T22(i,j,DDK) - &
                EL_demo001(rho(DDK),Lame_lambda(DDK),Lame_mu(DDK),&
                           x_coor,y_coor,time,5))

                 
              err_max(1)=max(err_max(1),err_loc(1))
              err_max(2)=max(err_max(2),err_loc(2))
              err_max(3)=max(err_max(3),err_loc(3))
              err_max(4)=max(err_max(4),err_loc(4))
              err_max(5)=max(err_max(5),err_loc(5))


           enddo
        enddo
     enddo


     
     write(lid,1000)time, err_max(1:5), dsqrt(sum(err_max(1:5)*err_max(1:5)))
     write(*,1000)time,err_max(1:5),dsqrt(sum(err_max(1:5)*err_max(1:5)))

     if (Last_Comp .eq. 1) then 
        close(lid)
     endif





1000 format(1e14.6,1x,6e12.4)


end subroutine
