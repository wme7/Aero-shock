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
  integer, save :: TotNum_DM
  integer, save :: PolyDegN_max(2)
  integer, save, allocatable      :: PolyDegN_DM(:,:) !(2,TotNum_DM)
  real(kind=8), save, allocatable :: DM_Vertex(:,:,:) !(4,2,TotNum_DM)
  real(kind=8), save, allocatable :: x1(:,:,:)      !(0:PolyDegN_Max(1:2),TotNum_DM)
  real(kind=8), save, allocatable :: x2(:,:,:)      !(0:PolyDegN_Max(1:2),TotNum_DM)
  integer,      save, allocatable :: DM_Connect(:,:,:)!(2,4,TotNum_DM)

  integer :: DDK, ND, ND1, ND2, ND3, Surf_Num, Edge_Num  ! [ND: Degree Index, 
                                                         !DDK: Domain Index]

  integer :: DDK_Connect, Surf_Connect, Patch_Type, Edge_Connect

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
    implicit none
    !  Declare local variable
    !  integer :: DDK
    
    !read in domain limit parameter
    open(20,file='domain.in',form='formatted',status='old')
    read(20, * ) !'====================== Input Parameters ==================='
    read(20, * )  TotNum_DM
    read(20, * )  PolyDegN_max(1)
    read(20, * )  PolyDegN_max(2)
    ! 
    !---allocate memory for Geometric Parameter
    call alloc_mem_domain_grid_variables()
    !                                        
    !---Read in Domain Geometric Parameter
    do DDK=1,TotNum_DM
       read(20, * ) !'-----------------------------------------------------------'
       read(20, * ) !'Parameters of Domain ### ----------------------------------'
       read(20, * ) !'-----------------------------------------------------------'
       read(20,200)  PolyDegN_DM(1,DDK), PolyDegN_DM(2,DDK)
       read(20, * ) !'-----------------------------------------------------------' 
       read(20,202)  DM_Vertex(1,1,DDK), DM_Vertex(1,2,DDK)
       read(20,202)  DM_Vertex(2,1,DDK), DM_Vertex(2,2,DDK)
       read(20,202)  DM_Vertex(3,1,DDK), DM_Vertex(3,2,DDK)
       read(20,202)  DM_Vertex(4,1,DDK), DM_Vertex(4,2,DDK)
       read(20, * ) !'-----------------------------------------------------------' 
       read(20,204)  DM_Connect(1,1,DDK),DM_Connect(2,1,DDK) 
       read(20,204)  DM_Connect(1,2,DDK),DM_Connect(2,2,DDK)
       read(20,204)  DM_Connect(1,3,DDK),DM_Connect(2,3,DDK)
       read(20,204)  DM_Connect(1,4,DDK),DM_Connect(2,4,DDK)
    enddo
       read(20, * ) !'-----------------------------------------------------------'
    close(20)
    ! finish domain parameters computation

    return
200 format(1x,i5,1x,i5)
202 format(3x,f9.4,4x,f9.4)
204 format(1x,i5,1x,i5)

  end subroutine Input_Domain_Grid_Parameters

  !==============================================================================
  Subroutine Geometric_Parameters_On_Screen( )
    implicit none
    
    ! Print Computation control parameters on screen
    write(*,* ) '====================== Input Parameters ==================='
    write(*,204)'ToTal Number of Domain : ', TotNum_DM
    write(*,205)'Maxium Degree of Polynomial :' , PolyDegN_max(1),&
                                                  PolyDegN_max(2)
    
    ! Set all approximation polynomial to same degree
    !   do DDK=2,TotNumDomain
    !      PolyDegreeN_Domain(DDK)=PolyDegreeN_Domain(1)
    !   enddo
    
    do DDK=1,TotNum_DM
       write(*, *)  '-----------------------------------------------------------'
       write(*, 203)'---------- Parameters of Domain', DDK,'--------------------'
       write(*, *)  '-----------------------------------------------------------'
       write(*, 207)  PolyDegN_DM(1,DDK), PolyDegN_DM(2,DDK)
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
       if ( (PolyDegN_DM(1,DDK) .gt. PolyDegN_Max(1)) .or. &
            (PolyDegN_DM(2,DDK) .gt. PolyDegN_Max(2)) ) then 
   
            write(*,*)'MD3D_Grid.f90:'
            write(*,*)'PolyDegN_DM > PolyDegN_Max for DDK=',DDK
            write(*,*)'Abort!'
            stop

       endif
    enddo

    return

!200 format(1x,i5)
201 format(3x,f8.4,4x,f8.4)
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

  subroutine alloc_mem_domain_grid_variables()
    !--Declare subroutine argument
    implicit none
    !   integer PolyDegreeN_max,TotNumDomain
    !
    !--Declare local argument
    integer ierr
    !
    !--Start subroutine
    !allocate required memory for basic domain parameters
    allocate(PolyDegN_DM(2,TotNum_DM), stat=ierr)
    if ( ierr .ne. 0) then
       write(*,*)'Can not allocatable PolyDegreeN_Domain variables'
    endif
    PolyDegN_DM=0
    
    allocate(DM_Vertex(4,2,TotNum_DM), stat=ierr)
    if ( ierr .ne. 0) then
       write(*,*)'Can not allocatable DM_Vertex variables'
    endif
    DM_Vertex=0.d0
     
    allocate( x1(0:PolyDegN_Max(1),0:PolyDegN_Max(2),1:TotNum_DM), &
              x2(0:PolyDegN_Max(1),0:PolyDegN_Max(2),1:TotNum_DM), stat=ierr)

    if (ierr .ne. 0) then
       write(*,*)'Can not allocate variable x and y '
    endif
    x1=0.d0; x2=0.d0

    allocate(DM_Connect(1:3,1:4,1:TotNum_DM), stat=ierr)
    if (ierr .ne. 0) then 
       write(*,*)'Can not allocate veriable DM_Connect'
    endif
    DM_Connect=0
    
    return
    
  end subroutine alloc_mem_domain_grid_variables


  !==============================================================================
  


end module MD2D_Grid
