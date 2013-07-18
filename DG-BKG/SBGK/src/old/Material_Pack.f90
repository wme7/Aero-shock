subroutine Init_Material_Parameters
  use MD2D_Grid
  use Material_Var
  implicit none
  ! declare local arguments
  integer:: lid, count
  
  ! allocate memory for storing epsilon and mu
  call  alloc_mem_material(maxval(PolyDegN_Max(1:2)),TotNum_DM)

  lid=80
!  open(lid,file='Material.in',form='formatted', status='unknown')
  open(lid,file='Material.in',form='formatted', status='unknown')
  read(lid,*) !'============================================'
  do DDK=1,TotNum_DM
     read(lid,1000) count, rho(DDK), lame_mu(DDK), lame_lambda(DDK)
     write(*,*) count, rho(DDK), lame_mu(DDK), lame_lambda(DDK)
  enddo
  read(lid,*) !'============================================'
  close(lid)
  write(*,*)'complete Initializing Material Parameters'



1000 format(i6,3f10.4)
1001 format(A4)

end subroutine Init_Material_Parameters

subroutine Init_Hinv_Symmetrizer
  use MD2D_Grid
  use Material_Var
  implicit none
  ! declare local arguments
  integer:: i,j
  real(kind=8):: density, mu, lambda
  real(kind=8):: hinv_11, hinv_22, hinv_33, hinv_44, hinv_55
  real(kind=8):: hinv_35, hinv_53

  do DDK=1, TotNum_DM
     
     ND1 = PolyDegN_DM(1,DDK)
     ND2 = PolyDegN_DM(2,DDK)

     density=rho(DDK)
     mu=lame_mu(DDK)
     lambda=lame_lambda(DDK)

     hinv_11=1.d0/density;  
     hinv_22=1.d0/density
     hinv_33=lambda+2.d0*mu
     hinv_35=lambda
     hinv_44=mu
     hinv_53=lambda
     hinv_55=lambda+2.d0*mu

     do Edge_Num=1,4
        
        ND=ND1; if(mod(Edge_Num,2) .eq. 0) ND=ND2

        do i=0,ND
           H_inv(1:5,1:5,i,Edge_Num,DDK) = reshape( & 
                (/hinv_11 ,   0.d0  ,   0.d0  ,   0.d0  ,   0.d0  , &
                    0.d0  , hinv_22 ,   0.d0  ,   0.d0  ,   0.d0  , &
                    0.d0  ,   0.d0  , hinv_33 ,   0.d0  , hinv_53 , &
                    0.d0  ,   0.d0  ,   0.d0  , hinv_44 ,   0.d0  , &
                    0.d0  ,   0.d0  , hinv_35 ,   0.d0  , hinv_55   &
                /) , (/5,5/) )
        enddo

     enddo ! Edge_Num
     
  enddo ! DDK

  return 


end subroutine Init_Hinv_Symmetrizer
