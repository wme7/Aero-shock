module filter
implicit none
real(kind=8), save :: cutoff
integer, save :: filter_order
real(kind=8), save :: filter_alpha
real(kind=8), save, allocatable :: filter_sigma(:)
real(kind=8), save, allocatable :: ValPolyAtGrid(:,:)
real(kind=8), save, allocatable :: smoothing_oper_xi1(:,:,:)
real(kind=8), save, allocatable :: smoothing_oper_xi2(:,:,:)
real(kind=8), save, allocatable :: Q(:,:), Qt(:,:)
!
integer :: DegCutOff

contains

subroutine Read_In_Filter_Profile(Filter_OnOff)
implicit none
integer Filter_OnOff

integer :: lid

   lid = 21
   open(lid,file='filter.in', form='formatted', status='old')
   read(lid,*)
   read(lid,*)filter_order
   read(lid,*)cutoff		! number between [0 1 ]
   read(lid,*)
   close(lid)

   return

end subroutine Read_In_Filter_Profile

subroutine alloc_mem_filter_variables(PolyDegMax)
implicit none
integer :: PolyDegMax
!
integer :: ierr
!

   allocate (filter_sigma(0:PolyDegMax), &
	     ValPolyAtGrid(0:PolyDegMax,0:PolyDegMax), &
	     smoothing_oper_xi1(0:PolyDegMax,0:PolyDegMax,0:PolyDegMax), &
	     smoothing_oper_xi2(0:PolyDegMax,0:PolyDegMax,0:PolyDegMax), stat=ierr )

   if ( ierr .ne. 0 ) then
      write(*,*)'Can not allocate memory for filter variables'
      write(*,*)'Abort!'
      stop
   endif

   filter_sigma=0.d0
   ValPolyAtGrid=0.d0
   smoothing_oper_xi1=0.d0
   smoothing_oper_xi2=0.d0

   allocate (  Q(0:PolyDegMax,0:PolyDegMax), &
	      Qt(0:PolyDegMax,0:PolyDegMax), stat=ierr )
   if (ierr .ne. 0) then
      write(*,*) 'Can not allocate memory for Qt'
      write(*,*) 'abort!'
      stop
   endif
   Q=0.d0; Qt=0.d0

   return

end subroutine alloc_mem_filter_variables

subroutine filter_profile(PolyDeg )
implicit none
integer :: PolyDeg
!
!declare local arguments
integer :: i,j
real(kind=8) :: machine_zero, tmp

   machine_zero=1.0E-14
   filter_alpha=-log(machine_zero)
   DegCutOff=int(CutOff*dfloat(PolyDeg))
!   write(*,*)filter_order,CutOff,PolyDeg; pause
   do i=0,PolyDeg

      if ( i .ge. DegCutOff ) then
	 tmp=(dfloat(i-DegCutOff)/dfloat(PolyDeg-DegCutOff))**filter_order
	 filter_sigma(i)=dexp(-filter_alpha*tmp)
      else
	 filter_sigma(i)=1.d0
      endif

!	  filter_sigma(i)=1.d0
   enddo

   do i=0, PolyDeg
      if ( filter_sigma(i) .le. 1.0E-14 )  filter_sigma(i)=0.d0
   enddo

!   do i=0,PolyDeg
!      write(*,*)i,filter_sigma(i)
!   enddo

!   pause
   return

end subroutine filter_profile



subroutine Init_Filter_Oper(ND_max)
use Legendre
implicit none
integer :: ND_max
real(kind=8), allocatable :: PnAtgrid(:),dPnAtgrid(:),d2PnAtgrid(:)
real(kind=8), allocatable :: gamma_N(:)
!real(kind=8), save, allocatable :: Q(:,:), Qt(:,:)

integer :: ierr, i, j ,l ,k, ND1

   !allocate variables for LGL_grid and initialize it
!   ND1=PolyDegN_Max(1); ND2=PolyDegN_Max(2); ND_Max=max(ND1,ND2)
   ! Take the max of PolyDegreeN_max as leading dimention to set up
   ! the legendre gauss labatto grid and its differential matrix for
   ! each degree
!   call Alloc_Mem_LG_grid(ND_Max)
!   call Alloc_Mem_Leg_Diff_Mat(ND_Max)

   ! set up LGL grid and differential matrix computational space
   allocate(  PnAtgrid(0:ND_max), &
	     dPnAtgrid(0:ND_max), &
	    d2PnAtgrid(0:ND_max), &
	    gamma_N(0:ND_max), stat=ierr)
   if (ierr .ne. 0) then
      write(*,*) 'Can not allocate memory for ValPolyDegNAtgrid'
      write(*,*) 'abort!'
      stop
   endif
   PnAtgrid=0.d0; dPnAtgrid=0.d0; d2PnAtgrid=0.d0
   gamma_N(0:ND_max)=0.d0;

   do ND1=4,ND_Max  ! compute filter for each degree

      call filter_profile(ND1)

      call init_LGL_grids(ND1)	! find the LGL grid for each degree
				! LG_weights(0:ND1) is also computed

      do l=0,ND1

	 do i=0,ND1
	    gamma_N(i)=2.d0/(2.d0*dfloat(i)+1.d0)
	 enddo
	    gamma_N(ND1)=2.d0/dfloat(ND1)

	 do k=0,ND1
	    call VALEPO(l,LG_grids(k),PnAtgrid(k),dPnAtgrid(k),d2PnAtgrid(k))
	 enddo

	 do i=0,ND1
	    ValPolyAtGrid(i,l)=PnAtgrid(i)
!	     write(*,*)i,PnAtgrid(i),LG_weights(i)
	 enddo
!	  pause
      enddo

      Q(0:ND1,0:ND1)=ValPolyAtGrid(0:ND1,0:ND1)
      Qt(0:ND1,0:ND1)=transpose(ValPolyAtGrid(0:ND1,0:ND1))

      do j=0,ND1
	 do i=0,ND1

	     Q(i,j) = filter_sigma(j)*Q(i,j)/gamma_N(j)
	    Qt(i,j) = LG_weights(j)*Qt(i,j)
!	     write(*,*)i,j,Q(i,j),Qt(i,j)

	 enddo
      enddo

      smoothing_oper_xi1(0:ND1,0:ND1,ND1)=matmul(Q(0:ND1,0:ND1),Qt(0:ND1,0:ND1))
      smoothing_oper_xi2(0:ND1,0:ND1,ND1)= &
		     transpose(smoothing_oper_xi1(0:ND1,0:ND1,ND1))

!      do j=0,ND1
!	  do i=0,ND1
!	     write(*,*)i,j,smoothing_oper_xi1(i,j,ND1)
!	  enddo
!      enddo

   enddo


!   deallocate()

   return

end subroutine Init_Filter_Oper

end module
