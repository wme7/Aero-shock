module MD2D_Grid
!--------------------------------------------------------------------------------
!
! This module contains functions and subroutines to initialize required 
! geometrical parameters for 3 dimensional multidomain pseudospectral 
! scheme
!
!--------------------------------------------------------------------------------

  implicit none
  !Declare Domain Variable

  integer, save :: TotNum_DM,TotNumglo_edg,TotNumin_edg,TotNumout_edg,Numin_edg1,Numin_edg2
  integer, save :: TotNum_EG,TotNum_in_EG,in1_EG,in2_EG,TotNum_out_EG,out1_EG,out2_EG,out3_EG,out4_EG
  integer, save :: PolyNodes(2),PND1,PND2,PDeg1,PDeg2

  integer, save, allocatable :: DM_Connect(:,:,:)!(2,4,TotNum_DM)
  integer, save, allocatable :: gloedg_Connect(:,:)!(Totedg,5)
  integer, save, allocatable :: inedg_Connect(:,:)!(Totinedg,5)
  integer, save, allocatable :: outedg_Connect(:,:)!(Totoutedg,5)
  integer, save, allocatable :: DMBnd_Shape(:,:)      !(1:4,TotNumDM)
  integer, save, allocatable :: EG_glo_Connect(:,:)!(1:5,1:TotNum_EG)
  integer, save, allocatable :: EG_in_Connect(:,:)!(1:5,1:TotNum_in_EG)
  integer, save, allocatable :: EG_out_Connect(:,:)!(1:5,1:TotNum_out_EG)
  integer, save, allocatable :: DGPolyDeg(:,:) !(2,TotNum_DM)

  character(len=5) :: NumDM

  real(kind=8), save :: a(1:5),lx(1:5),ly(1:5),T1(1:5),lambda,mu,rho,c_s,c_p
  real(kind=8), save :: error_L2_tau11,error_L2_tau22,error_L2_tau12,error_L2_v1,error_L2_v2
  real(kind=8), save :: time_start,time_end,sqTDM,ttau_scale,dx,CCFL,cfl_in

  real(kind=8), save, allocatable :: DM_Vertex(:,:,:) !(4,2,TotNum_DM)
  real(kind=8), save, allocatable :: x1(:,:,:)      !(0:PolyDegN_Max(1:2),TotNum_DM)
  real(kind=8), save, allocatable :: x2(:,:,:)      !(0:PolyDegN_Max(1:2),TotNum_DM)
  real(kind=8), save, allocatable :: DMBnd_Shape_Par(:,:,:)

  integer :: DDK, ND, ND1, ND2, ND3, Surf_Num, Edge_Num  ! [ND: Degree Index, 
                                                         !DDK: Domain Index]
  integer :: DDK_Connect, Edge_Connect, Patch_Type, edg

  integer :: DM_Range
  integer :: side_num
! Description of variables
!================================================================================
! [Name] :: PolyDegN_DM(Index1,Index2)
!  ^^^^
!    [Size]    :: (1:2,1:TotNum_DM)
!     ^^^^
!    [Purpose] :: Store Polynomial Degree used in each domain
!     ^^^^^^^
!    [Detail]  :: Index2 is used to represent domain number
!     ^^^^^^
!                 Index1 is used to distinguish the number of Polynomial
!                 Degree used in x1 and x2 directions
!
!     Example :: PolyDeg_DM(1,DDK) stores the following:
!                    The degree of the polynomial used in the x1 direction 
!                    of Domain DDK.
!                  
!                PolyDeg_DM(2,DDK) stores the following:
!                    The degree of the polynomial used in the x2 direction
!                    of Domain DDK. 
! 
! [Name] :: DM_Vertex(Index1,Index2)
! ^^^^
!     Size    :: (1:5,1:2,1:TotNum_DM) 
!     ^^^^
!     Purpose :: To store the 4 physical coordinates of each domain  
!     ^^^^^^^
!     Detail  :: The first index indicate the four vertex the fifth one 
!     ^^^^^^     is a copy of the first one for coding convience
!                The second index stands for x and y (1 for x, and 2 for y)
!                The third index is the number of the Domain.
!
! [Name] :: DM_Connect(index1,index2,index3) for 2D problem
! ^^^^
!     Size    ::  (3,4,1:TotNum_DM)
!     ^^^^
!     Purpose :: This variable is used to store the required info for  
!     ^^^^^^^    patching bc between two connecting domains
!              
!     Detail  :: Index3 is used to denoted the domain number 
!     ^^^^^^     
!                Index2 is used to denote the 4 edges: 
!                    1: xi_2= -1
!                    2: xi_1=  1
!                    3: xi_2=  1
!                    4: xi_1= -1
!                    
!
!                Index1 is used to denote the store info:
!                    1 for the number of the connected domain
!                    2 for the connecting edge
!                    3 for patching direction
!
!     Example :: DM_Connect(1,1,DDK) stores the following:
!                    The number of the connected domain on the 
!                    1st edge of domain DDK
!                  
!                DM_Connect(2,1,DDK) stores the following:
!                    The edge number of the connected domain on the 
!                    1st edge of domain DDK
!
!                DM_Connect(3,1,DDK) stores the following:
!                    The patching direction of the connected domain on the 
!                    1st edge of domain DDK. 1 for the same and 2 for opposite
!================================================================================

contains
  

  !==============================================================================
  subroutine Input_Domain_Grid_Parameters()
  use State_Var
  use Kinetic_Var
    implicit none
    integer:: lid, ierr, my_choice,PND,PDeg,ls

!    character(len=5) :: NumDM
    character(len=30) :: DomainInFile
  
    write(*,*) "Number of Domains="
    read(*,"(5a)") NumDM
  

     
    lid=20
    IT=0
  open(lid,file= './Parameter/my.in',form='formatted', status='old')
 read(lid, * ) !'====================== Input Parameters ==================='
 read(lid, * )  my_choice
 if (my_choice .gt. 0) then
 read(lid, * )   PND
 read(lid, * )   PDeg
 read(lid, * )   IT
 read(lid, * )   I_Case
 read(lid, * )   cfl_in
 read(lid, * )   DM_Range
endif
close(lid)

    select case(my_choice)
    case(1) !Old Case
    DomainInFile = './Parameter/'//NumDM//'dom.in'
    I_prob=1
    case(2) !Periodic
    DomainInFile = './Parameter/'//NumDM//'Pdom.in'
    I_prob=2
    case(3) !Neumann BC
    DomainInFile = './Parameter/'//NumDM//'Ndom.in'
    I_prob=3
    end select

    
    write(*,*) DomainInFile

    open(lid,file=DomainInFile,form='formatted', status='old', IOSTAT=ierr)
    if (ierr .gt. 0) then
        write(*,*)'Message from MD2D_Grid.f90'
        write(*,*)'Open file '//DomainInFile//' failure!'
        write(*,*)'The file may not exist.'
        write(*,*)'Stop!'
        close(lid)
        stop
    endif
  


    !read in domain limit parameter
    read(lid, * ) !'====================== Input Parameters ==================='
    read(lid, * )  TotNum_DM
    read(lid, * )  PolyNodes(1)
    read(lid, * )  PolyNodes(2)
    ! 
    !---allocate memory for Geometric Parameter
 if (my_choice .gt. 0)  PolyNodes(:) = PND

    !PolyNodes(:) = 16
    PND1 = PolyNodes(1)
    PND2 = PolyNodes(2)

    call alloc_mem_domain_grid_variables
    !                                        
    !---Read in Domain Geometric Parameter
    do DDK=1,TotNum_DM
       read(lid, * ) !'-----------------------------------------------------------'
       read(lid, * ) !'Parameters of Domain ### ----------------------------------'
       read(lid, * ) !'-----------------------------------------------------------'
       read(lid, * )  DGPolyDeg(1,DDK), DGPolyDeg(2,DDK)
       read(lid, * ) !'-----------------------------------------------------------' 
       read(lid, * )  DM_Vertex(1,1,DDK), DM_Vertex(1,2,DDK)
       read(lid, * )  DM_Vertex(2,1,DDK), DM_Vertex(2,2,DDK)
       read(lid, * )  DM_Vertex(3,1,DDK), DM_Vertex(3,2,DDK)
       read(lid, * )  DM_Vertex(4,1,DDK), DM_Vertex(4,2,DDK)
       read(lid, * ) !'-----------------------------------------------------------' 
       read(lid, * )  DM_Connect(1,1,DDK),DM_Connect(2,1,DDK) 
       read(lid, * )  DM_Connect(1,2,DDK),DM_Connect(2,2,DDK)
       read(lid, * )  DM_Connect(1,3,DDK),DM_Connect(2,3,DDK)
       read(lid, * )  DM_Connect(1,4,DDK),DM_Connect(2,4,DDK)
    enddo
    read(lid, * ) !'-----------------------------------------------------------'
    close(lid)
    ! finish domain parameters computation
 if (DM_Range .gt. 0) then
    do DDK=1,TotNum_DM
     do ls=1,4 
      DM_Vertex(ls,1,DDK)=(DM_Vertex(ls,1,DDK)+1d0)/2d0
      DM_Vertex(ls,2,DDK)=(DM_Vertex(ls,2,DDK)+1d0)/2d0
     enddo
    enddo
endif


 if (my_choice .gt. 0) DGPolyDeg(:,:) = PDeg

    !DGPolyDeg(:,:) = 15
    PDeg1 = DGPolyDeg(1,1)
    PDeg2 = DGPolyDeg(2,1)

    select case(my_choice)
    case(1) !Old Case
    DomainInFile = './Parameter/'//NumDM//'edg.in'
    case(2) !Periodic
    DomainInFile = './Parameter/'//NumDM//'Pedg.in'
    case(3) !Neumann BC
    DomainInFile = './Parameter/'//NumDM//'Nedg.in'
    end select

    
    write(*,*) DomainInFile
     
    lid=20
    open(lid,file=DomainInFile,form='formatted', status='old', IOSTAT=ierr)
    if (ierr .gt. 0) then
        write(*,*)'Message from MD2D_Grid.f90'
        write(*,*)'Open file '//DomainInFile//' failure!'
        write(*,*)'The file may not exist.'
        write(*,*)'Stop!'
        close(lid)
        stop
    endif

    read(lid, * ) !'====================== Input Parameters ==================='
    read(lid, * )  TotNum_EG
    read(lid, * ) !'-----------------------------------------------------------' 
    call alloc_mem_edge_variables1                                        
    !---Read in Domain Geometric Parameter
    do DDK=1,TotNum_EG
       read(lid, * )  EG_glo_Connect(1,DDK),EG_glo_Connect(2,DDK),EG_glo_Connect(3,DDK),EG_glo_Connect(4,DDK),EG_glo_Connect(5,DDK)
    enddo
    read(lid, * ) !'-----------------------------------------------------------'
    read(lid, * )  TotNum_in_EG
    read(lid, * )  in1_EG
    read(lid, * )  in2_EG
    read(lid, * ) !'-----------------------------------------------------------'

    call alloc_mem_edge_variables2 
    do DDK=1,TotNum_in_EG
       read(lid, * )  EG_in_Connect(1,DDK),EG_in_Connect(2,DDK),EG_in_Connect(3,DDK),EG_in_Connect(4,DDK),EG_in_Connect(5,DDK)
    enddo
    read(lid, * ) !'-----------------------------------------------------------'
    read(lid, * )  TotNum_out_EG
    read(lid, * )  out1_EG
    read(lid, * )  out3_EG
    read(lid, * )  out2_EG
    read(lid, * )  out4_EG
    read(lid, * ) !'-----------------------------------------------------------'

    call alloc_mem_edge_variables3 
    do DDK=1,TotNum_out_EG
       read(lid, * )  EG_out_Connect(1,DDK),EG_out_Connect(2,DDK),EG_out_Connect(3,DDK),EG_out_Connect(4,DDK),EG_out_Connect(5,DDK)
    enddo
    read(lid, * ) !'-----------------------------------------------------------'
    close(lid)

    return
200 format(1x,i5,1x,i5)
202 format(3x,f13.8,4x,f13.8)
204 format(1x,i5,1x,i5,1x,i5,1x,i5,1x,i5)

  end subroutine Input_Domain_Grid_Parameters

  !==============================================================================
  Subroutine Geometric_Parameters_On_Screen( )
    implicit none
    
    ! Print Computation control parameters on screen
    write(*,* ) '====================== Input Parameters ==================='
    write(*,204)'ToTal Number of Domain : ', TotNum_DM
    write(*,205)'Maxium Degree of Polynomial :' , PolyNodes(1),&
                                                  PolyNodes(2)
    
    ! Set all approximation polynomial to same degree
    !   do DDK=2,TotNumDomain
    !      PolyDegreeN_Domain(DDK)=PolyDegreeN_Domain(1)
    !   enddo
    
    do DDK=1,TotNum_DM
       write(*, *)  '-----------------------------------------------------------'
       write(*, 203)'---------- Parameters of Domain', DDK,'--------------------'
       write(*, *)  '-----------------------------------------------------------'
       write(*, 207)  DGPolyDeg(1,DDK), DGPolyDeg(2,DDK)
       write(*, *)  '-----------------------------------------------------------'
       write(*, 201)  DM_Vertex(1,1,DDK), DM_Vertex(1,2,DDK)
       write(*, 201)  DM_Vertex(2,1,DDK), DM_Vertex(2,2,DDK)
       write(*, 201)  DM_Vertex(3,1,DDK), DM_Vertex(3,2,DDK)
       write(*, 201)  DM_Vertex(4,1,DDK), DM_Vertex(4,2,DDK)
       write(*, *)  '-----------------------------------------------------------'  
       write(*, 207)  DM_Connect(1,1,DDK),DM_Connect(2,1,DDK)
       write(*, 207)  DM_Connect(1,2,DDK),DM_Connect(2,2,DDK)
       write(*, 207)  DM_Connect(1,3,DDK),DM_Connect(2,3,DDK)
       write(*, 207)  DM_Connect(1,4,DDK),DM_Connect(2,4,DDK)
    enddo
       write(*, *)  '-----------------------------------------------------------'


    ! Check consistency of PolyDegN_DM and  PolyDegN_max
    do DDK=1,TotNum_DM
       if ( (DGPolyDeg(1,DDK) .gt. PolyNodes(1)) .or. &
            (DGPolyDeg(2,DDK) .gt. PolyNodes(2)) ) then 
   
            write(*,*)'MD3D_Grid.f90:'
            write(*,*)'PolyDegN_DM > PolyDegN_Max for DDK=',DDK
            write(*,*)'Abort!'
            stop

       endif
    enddo

    return

!200 format(1x,i5)
201 format(3x,f13.8,4x,f13.8)
203 format(A32,1x,i3,2x,A22)
204 format(A32,1x,i4)
205 format(A32,1x,i4,1x,i4,1x,i4)
207 format(1x,i5,1x,i5)

  end subroutine Geometric_parameters_On_Screen

  !==============================================================================

  subroutine Init_patch_direction()
    implicit none
     
    DM_Connect(3,1:6,1:TotNum_DM)=1
    
    do DDK=1,TotNum_DM
     
       do Edge_Num=1,2
       
          select case (DM_Connect(2,Edge_Num,DDK)) ! side_number 1,2 of DDK
          case(1,2) ! connect to side_num=1,2: reverse
             DM_Connect(3,Edge_Num,DDK)=-1
             
          end select
          
       enddo
       
       do Edge_Num=3,4
          
          select case (DM_Connect(2,Edge_Num,DDK)) ! side_number 3,4 of DDK
          case(3,4) ! connect to side_num=3,4: reverse
             DM_Connect(3,Edge_Num,DDK)=-1
             
          end select
          
       enddo
       
    enddo
    
    return 
    
  end subroutine Init_patch_direction
 
  !==============================================================================

  subroutine alloc_mem_domain_grid_variables
    !--Declare subroutine argument
    implicit none
    !   integer PolyDegreeN_max,TotNumDomain
    !
    !--Declare local argument
    integer ierr
    !
    !--Start subroutine
    !allocate required memory for basic domain parameters
    allocate(DGPolyDeg(2,TotNum_DM), stat=ierr)
    if ( ierr .ne. 0) then
       write(*,*)'Can not allocatable PolyDegreeN_Domain variables'
    endif
    DGPolyDeg=0
    
    allocate(DM_Vertex(4,2,TotNum_DM), stat=ierr)
    if ( ierr .ne. 0) then
       write(*,*)'Can not allocatable DM_Vertex variables'
    endif
    DM_Vertex=0.d0
     
    allocate( x1(0:PND1,0:PND2,1:TotNum_DM), &
              x2(0:PND1,0:PND2,1:TotNum_DM), stat=ierr)

    if (ierr .ne. 0) then
       write(*,*)'Can not allocate variable x and y '
    endif
    x1=0.d0; x2=0.d0

    allocate(DM_Connect(1:3,1:4,1:TotNum_DM), stat=ierr)
    if (ierr .ne. 0) then 
       write(*,*)'Can not allocate veriable DM_Connect'
    endif
    DM_Connect=0
    
    
    allocate(DMBnd_Shape(4,TotNum_DM), stat=ierr)
    if ( ierr .ne. 0) then 
      write(*,*)'Can not allocate DMBnd_Shape variable'
    endif

    allocate(DMBnd_Shape_Par(1:7,4,TotNum_DM), stat=ierr)
    if (ierr .ne. 0) then 
       write(*,*)'Can not allocate DMB_Shape_Par variable'
    endif
    
    
    return
    
  end subroutine alloc_mem_domain_grid_variables

  !==============================================================================
  subroutine alloc_mem_edge_variables1

  integer ierr

  allocate(EG_glo_Connect(1:5,1:TotNum_EG), stat=ierr)
  if (ierr .ne. 0) then
     write(*,*)'Can not allocate veriable EG_glo_Connect'
  endif

  end subroutine alloc_mem_edge_variables1
  !==============================================================================
  subroutine alloc_mem_edge_variables2

  integer ierr

  allocate(EG_in_Connect(1:5,1:TotNum_in_EG), stat=ierr)
  if (ierr .ne. 0) then
     write(*,*)'Can not allocate veriable EG_glo_Connect'
  endif

  end subroutine alloc_mem_edge_variables2
  !==============================================================================
  subroutine alloc_mem_edge_variables3

  integer ierr

  allocate(EG_out_Connect(1:5,1:TotNum_out_EG), stat=ierr)
  if (ierr .ne. 0) then
     write(*,*)'Can not allocate veriable EG_glo_Connect'
  endif

  end subroutine alloc_mem_edge_variables3
  !==============================================================================
  subroutine alloc_mem_edge_grid_variables

    integer ierr

    allocate(gloedg_Connect(1:TotNumglo_edg,1:5), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate veriable DM_Connect'
    endif
    gloedg_Connect=0
  
    if (TotNumin_edg .GT. 0) then
        allocate(inedg_Connect(1:TotNumin_edg,1:5), stat=ierr)
        if (ierr .ne. 0) then
           write(*,*)'Can not allocate veriable DM_Connect'
        endif
    endif
    inedg_Connect=0

    if (TotNumout_edg .GT. 0) then
        allocate(outedg_Connect(1:TotNumout_edg,1:5), stat=ierr)
        if (ierr .ne. 0) then
            write(*,*)'Can not allocate veriable DM_Connect'
        endif
    endif
    outedg_Connect=0

  end subroutine alloc_mem_edge_grid_variables

  !==============================================================================
  
end module MD2D_Grid
