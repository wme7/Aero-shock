module Legendre
  implicit none
  real(kind=8), save, allocatable:: LG_grids(:),LG_weights(:)
  real(kind=8), save, allocatable:: ValuesOfPolyNatGrids(:)
  real(kind=8), save, allocatable:: LG_DiffMat(:,:)

  real(kind=8), save, allocatable:: F_alt(:,:,:,:,:),F_t(:,:,:,:,:),FS(:,:,:,:,:)

  real(kind=8), save, allocatable:: F_0(:,:,:,:,:)
  real(kind=8), save, allocatable:: F_loc(:,:,:,:,:)
  real(kind=8), save, allocatable:: F_new(:,:,:,:,:),exa_solij(:,:,:,:,:)
  real(kind=8), save, allocatable:: flx_F(:,:,:,:,:)
  real(kind=8), save, allocatable:: FS_tmp(:,:)


  real(kind=8), save, allocatable:: GH(:),GHW(:),Vx(:),Vy(:)
  integer:: IGH=20


!  real(kind=8), save, allocatable:: u_alt(:,:,:),u_t(:,:,:),u_new(:,:,:),exa_sol(:,:,:)
!  real(kind=8), save, allocatable:: flx_u(:,:,:)

  real(kind=8), save, allocatable:: B_tal1x(:,:),B_tal1y(:,:),B_tal2x(:,:),A1(:,:,:),Bq(:,:,:)
!  real(kind=8), save, allocatable:: B(:,:)
  real(kind=8), save, allocatable:: B_gen(:,:,:)

  
  real(kind=8), save, allocatable:: D1_DiffMat(:,:,:) !(0:N_max,0:N_max,Degree_max)
  real(kind=8), save, allocatable:: D1_Vertex(:,:,:)  !(0:N_max,2,Degree_max)
  real(kind=8), save, allocatable:: D2_DiffMat(:,:,:)
  real(kind=8), save, allocatable:: LGLCoord(:,:)     !(0:N_max,Degree_Max)
  real(kind=8), save, allocatable:: LGLWeights(:,:)   !(0:N_max,Degree_max)
  real(kind=8), save, allocatable:: leg_tb(:,:),Dleg_tb(:,:)

  real(kind=8), save, allocatable:: Leg_Grid_xi1(:,:,:),Leg_Grid_xi2(:,:,:)
  real(kind=8), save, allocatable:: DLeg_Grid_xi1(:,:,:),DLeg_Grid_xi2(:,:,:)
  real(kind=8), save, allocatable:: LGLWeights_Grid_xi1(:,:),LGLWeights_Grid_xi2(:,:)

  real(kind=8), save, allocatable:: Diff_xi1(:,:,:)
  real(kind=8), save, allocatable:: Diff_xi2(:,:,:)

contains

  !==============================================================================
  subroutine Init_LGL_Diff_Matrix_23D(N_max,N_DM)
    implicit none 
    integer:: N_max, N_DM
    
    integer:: N13_max
    integer:: N
    integer:: ierr

    allocate(  LGLCoord(0:N_max,N_max), &
             LGLWeights(0:N_max,N_max), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate memory for LGLCoord and LGLWeights'
       stop
    endif
    LGLCoord=0.d0
    LGLWeights=0.d0

    ! allocate memory for differential matrix in xi1 direction
    allocate(Diff_xi1(0:N_Max,0:N_Max,N_Max), &
             Diff_xi2(0:N_Max,0:N_Max,N_Max), stat=ierr)

    if (ierr .ne. 0) then
       write(*,*)'Legendre.f90:' 
       write(*,*)'Fail to allocate memory for Diff_xi1 and Diff_xi2'
       write(*,*)'Abort!'
       stop 
    endif

    Diff_xi1=0.d0
    Diff_xi2=0.d0

    !allocate variables for LGL_grid and initialize it

    call Alloc_Mem_LG_grid(N_max) 

    call Alloc_Mem_Leg_Diff_Mat(N_max)
    
    ! set up LGL grid and differential matrix computational space
    
    do N=1,N_max
       
       call init_LGL_grids(N)
       
       LGLCoord(0:N,N)   = LG_grids(0:N)
       LGLWeights(0:N,N) = LG_weights(0:N)

       call init_diff_matrix(N,N_max)  ! out put LG_DiffMat
 
       Diff_xi1(0:N,0:N,N)=LG_DiffMat(0:N,0:N)

       Diff_xi2(0:N,0:N,N)=transpose(LG_DiffMat(0:N,0:N))

    enddo           
    
    call MLPN(N_max,N_DM)
  
    return

  end subroutine Init_LGL_Diff_Matrix_23D
  !=============================================================================
  subroutine MLPN(N_max,DOF)
    implicit none
    integer:: i,N_max,DOF,ierr
    real(kind=8):: x,PD(0:DOF),PN(0:DOF)   

    allocate(  leg_tb(0:N_max,0:DOF), &
               Dleg_tb(0:N_max,0:DOF), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate memory for leg_tb'
       stop
    endif

    do i=0,N_max
       x=LGLCoord(i,N_max)
       call LPN(DOF,x,PN,PD)
       leg_tb(i,0:DOF) = PN(0:DOF) 
       Dleg_tb(i,0:DOF) = PD(0:DOF)
    enddo

  end subroutine MLPN
  !==============================================================================
  subroutine LPN(N,X,PN,PD)
    implicit none
    integer::N,k
    real(kind=8)::PN(0:N),PD(0:N),PF,P0,P1,X
    
    PN(0)=1
    PN(1)=X
    PD(0)=0
    PD(1)=1
    
    do k=2,N
       PN(k) = ((2*k-1)*X*PN(k-1)-(k-1)*PN(k-2))/k
       PD(k) = (2*k-1)*PN(k-1)+PD(k-2)
    enddo

  end subroutine LPN
  !==============================================================================
  subroutine Init_LGL_Diff_Matrix_1D(Degree_Max)
    implicit none
    integer:: Degree_Max
    
    integer:: N       ! Degree of the Polynomial
    integer:: ierr
    
    ! Differential Matrix of all Degree
    allocate(D1_DiffMat(0:Degree_Max,0:Degree_Max,Degree_Max),stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate memory for D1_DiffMat'
       stop
    endif
    D1_DiffMat=0.d0
    
    allocate(D2_DiffMat(0:Degree_Max,0:Degree_Max,Degree_Max),stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate memory for D2_DiffMat'
       stop
    endif
    D2_DiffMat=0.d0
    
    allocate(D1_Vertex(0:Degree_Max,1:2,Degree_Max),stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate memory for D1_DiffMat_Vertex'
       stop
    endif
    D1_Vertex=0.d0

    allocate(  LGLCoord(0:Degree_Max,Degree_Max), &
             LGLWeights(0:Degree_Max,Degree_Max), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate memory for LGLCoord and LGLWeights'
       stop
    endif
    LGLCoord=0.d0
    LGLWeights=0.d0
    
    !allocate variables for LGL_grid and initialize it

    call Alloc_Mem_LG_grid(Degree_Max)

    call Alloc_Mem_Leg_Diff_Mat(Degree_Max)
    
    ! set up LGL grid and differential matrix computational space
    
    do N=1,Degree_Max
       
       call init_LGL_grids(N)
       
         LGLCoord(0:N,N) =   LG_grids(0:N);

       LGLWeights(0:N,N) = LG_weights(0:N)
       
       call init_diff_matrix(N,Degree_Max)  ! out put LG_DiffMat
       
       ! set up D1 matrix and  boundary operator
      
       D1_DiffMat( 0:N , 0:N , N )=LG_DiffMat(0:N,0:N)
        D1_Vertex( 0:N ,  1  , N )=LG_DiffMat(0,0:N)
        D1_Vertex( 0:N ,  2  , N )=LG_DiffMat(N,0:N)
      
            
      ! set up D2 matrix

       D2_DiffMat( 0:N , 0:N , N ) = Matmul(LG_DiffMat(0:N,0:N),&
                                            LG_DiffMat(0:N,0:N))

      
   
    enddo

    return
    
  end subroutine Init_LGL_Diff_Matrix_1D

  !==============================================================================
  subroutine init_diff_matrix(N,NM)
    implicit none
    integer N,NM !NM is the Maximum Degree
    
    !     subroutine begins
    call DMLEGL(N,NM,LG_grids(0),ValuesOfPolyNatGrids(0),LG_DiffMat(0,0))
    
    return
    
  end subroutine init_diff_matrix

  !==============================================================================
  subroutine init_LGL_grids(N)
    implicit none
    integer N
    
    !     subroutine begins
    call ZELEGL(N, LG_grids, ValuesOfPolyNatGrids)
    call WELEGL(N, LG_grids, ValuesOfPolyNatGrids, LG_weights)
    
    return
    
  end subroutine init_LGL_grids

  !==============================================================================
  subroutine Alloc_Mem_LG_grid(N)
    implicit none
    integer N ! degree of polynomial N and total number of grids is N+1
    
    integer ierr
    
    !     subroutine begin
    allocate(LG_grids(0:N), LG_weights(0:N), ValuesOfPolyNatGrids(0:N), &
         stat=ierr)
    if ( ierr .ne. 0 ) then
       write(*,*)'Can not allocate Legendre grids'
       write(*,*)'Abort! '
       stop
    endif

    return
    
  end subroutine Alloc_Mem_LG_grid
  !==============================================================================
  subroutine Alloc_Mem_Leg_Diff_Mat(N)
    implicit none
    integer N
    
    integer ierr
    
    !     subroutine begins
    allocate(LG_DiffMat(0:N,0:N), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)' Can not allocatable Memory of variables LG_DiffMat '
       write(*,*)' Abort! '
       stop
    endif
    
  end subroutine Alloc_Mem_Leg_Diff_Mat
  !==============================================================================
  subroutine alloc_Mem_Leg_Grid(ND1,ND2,PD1)
  implicit none 
  integer::ierr,ND1,ND2,PD1

  allocate( Leg_Grid_xi1(0:ND1,0:ND1,0:PD1),Leg_Grid_xi2(0:ND1,0:ND2,0:PD1), &
            DLeg_Grid_xi1(0:ND1,0:ND1,0:PD1),DLeg_Grid_xi2(0:ND1,0:ND2,0:PD1), &    
            LGLWeights_Grid_xi1(0:ND1,0:ND2),LGLWeights_Grid_xi2(0:ND1,0:ND2),stat=ierr)
  
  if (ierr .ne. 0) then
     write(*,*) "Can't allocate Leg_Grid"
     stop
  endif
  
  end subroutine alloc_Mem_Leg_Grid 

  !==============================================================================
  subroutine alloc_mem_u(PD1,PD2,PN1,PN2,TDM)
  implicit none
  integer::ierr,PD1,PD2,TDM,PN1,PN2

 allocate(  flx_F(1:IGH,1:IGH,0:PD1,0:PD2,1:TDM) )
  allocate( F_alt(1:IGH,1:IGH,0:PD1,0:PD2,1:TDM),&
            F_new(1:IGH,1:IGH,0:PD1,0:PD2,1:TDM), &
              F_t(1:IGH,1:IGH,0:PD1,0:PD2,1:TDM),&
               FS(1:IGH,1:IGH,0:PD1,0:PD2,1:TDM),&
              F_0(1:IGH,1:IGH,0:PN1,0:PN2,1:TDM),&
            F_loc(1:IGH,1:IGH,0:PN1,0:PN2,1:TDM),&
            FS_tmp(0:PN1,0:PN2),&
            exa_solij(1:IGH,1:IGH,0:PN1,0:PN2,1:TDM), stat=ierr)

  if (ierr .ne. 0) then
     write(*,*) "Can't allocate F"
     stop
  endif

    allocate(GH(1:IGH),GHW(1:IGH),Vx(1:IGH),Vy(1:IGH),B_gen(2,1:IGH,1:IGH), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)' Message from Legendre.f90'
       write(*,*)' Cannot allocate memory for variable GH'
       write(*,*)' Abort!'
       stop
    endif

!  allocate( flx_u(0:PD1,0:PD2,1:TDM) )
!  allocate( u_alt(0:PD1,0:PD2,1:TDM),u_new(0:PD1,0:PD2,1:TDM), &
!            u_t(0:PD1,0:PD2,1:TDM),exa_sol(0:PN1,0:PN2,1:TDM)

  allocate( B_tal1x(0:PD1,0:PD2), &
            B_tal1y(0:PD1,0:PD2),B_tal2x(0:PD1,0:PD2), &
            A1(0:PD1,0:PD2,1:TDM),&
            Bq(1:2,1:2,1:TDM),stat=ierr)

  if (ierr .ne. 0) then
     write(*,*) "Can't allocate u"
     stop
  endif

!allocate( B(1:2,1:TDM))
  



  end subroutine alloc_mem_u
  !==============================================================================

end module Legendre
