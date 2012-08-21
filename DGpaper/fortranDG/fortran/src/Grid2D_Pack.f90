!=======================================================
subroutine Init_Physical_Grid_2D  !.. assign coordination value
  use MD2D_Grid
  use Legendre
  use universal_const
  implicit none

  ! Declare local arguments
  ! Boundary grid for setting physical space grid by transfinite blending
  ! Function : Need to deallocate before leaving this subroutine
  
  real(kind=8), allocatable :: Edge_x(:,:), Edge_y(:,:)
  integer ierr
  integer N, i, j, k
  
  !Linear mapping parameters
  real(kind=8) :: x_start, x_end
  real(kind=8) :: y_start, y_end
  !..Circular arc mapping parameters
  real(kind=8) :: x_cen, y_cen, radius
  real(kind=8) :: v1x, v1y, Lv1, v2x, v2y, Lv2
  real(kind=8) :: angle_ref, cross_ref, cross, angle, t
  real(kind=8) :: focus1_x, focus1_y, focus2_x, focus2_y
  real(kind=8) :: maj_axis, min_axis, dis_focus, rot_angle
  
  integer:: ND_Max
  real(kind=8) :: xi1, xi2 
  
  !  Start subroutine
  !  Find the max of degree in all directions

  ND1=PND1
  ND2=PND2
  ND_Max=max(ND1,ND2)

  allocate( Edge_x(0:ND_Max,1:12), &
       	    Edge_y(0:ND_Max,1:12), stat= ierr)
  
  if ( ierr .ne. 0 ) then
     write(*,*) 'Grid2D_Pack.f90:'
     write(*,*) 'Can not allocate Edge_x Edge_y and in Grid2D_Pack'
     write(*,*) 'Abort!'
     stop
  endif
  
  !!call Input_Domain_Shape_Parameters

  DMBnd_Shape=1  ! Assign straight edge directly
  
! Construct Computational Grids by Transfinite Blending Mapping
! Step 1: Construct Points on Egdes (Totally 4 Edges)
! Step 2: Construct Points on Surfaces (Total 6 Surfaces) based on the Edges

  do DDK=1,TotNum_DM
     
     
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
     !
     ! Set up grid points on Edges 
     !
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
     
     ! Edge 1: v1 (-1,-1) to v2 ( 1,-1)
     ! ^^^^^^
     call AssignPts(1,ND1,DDK,Edge_x(:,1),Edge_y(:,1))
     
     !-------------------------------------------------------------!
     
     ! Edge 2: v2 ( 1,-1) to v3 ( 1, 1)
     ! ^^^^^^
     call AssignPts(2,ND2,DDK,Edge_x(:,2),Edge_y(:,2))
     
     !-------------------------------------------------------------!
     
     ! Edge 3: v4 (-1, 1) to v3 ( 1, 1) 
     ! ^^^^^^ 
     call AssignPts(3,ND1,DDK,Edge_x(:,3),Edge_y(:,3))
     
     !-------------------------------------------------------------!
     
     ! Edge 4: v1 (-1,-1) to v4 (-1, 1)
     ! ^^^^^^
     call AssignPts(4,ND2,DDK,Edge_x(:,4),Edge_y(:,4))
               
     !-------------------------------------------------------------!
     
     !-------------------------------------------------------------!
     ! This completes the points on Edges
     !-------------------------------------------------------------!
     
     !=============================================================!
     
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
     !
     ! Set up grid points on Surfaces 
     !
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
     
     do j=0,ND2
        
        xi2 = LGLCoord(j,ND2)
        
        do i=0,ND1
           
           xi1 = LGLCoord(i,ND1)
           
           ! Original Two D version
           x1(i,j,DDK) = &
                 (1.d0-xi2)/2.d0 * Edge_x(i,1) + &
                 (1.d0+xi2)/2.d0 * Edge_x(i,3) + &
                 (1.d0-xi1)/2.d0 * Edge_x(j,4) + &
                 (1.d0+xi1)/2.d0 * Edge_x(j,2) - &
                 (1.d0-xi1)*(1.d0-xi2)/4.d0 * DM_Vertex(1,1,DDK)- &
                 (1.d0+xi1)*(1.d0-xi2)/4.d0 * DM_Vertex(2,1,DDK)- &
                 (1.d0+xi1)*(1.d0+xi2)/4.d0 * DM_Vertex(3,1,DDK)- &
                 (1.d0-xi1)*(1.d0+xi2)/4.d0 * DM_Vertex(4,1,DDK)

           x2(i,j,DDK) = &
                 (1.d0-xi2)/2.d0 * Edge_y(i,1) + &
                 (1.d0+xi2)/2.d0 * Edge_y(i,3) + &
                 (1.d0-xi1)/2.d0 * Edge_y(j,4) + &
                 (1.d0+xi1)/2.d0 * Edge_y(j,2) - &
                 (1.d0-xi1)*(1.d0-xi2)/4.d0 * DM_Vertex(1,2,DDK)- &
                 (1.d0+xi1)*(1.d0-xi2)/4.d0 * DM_Vertex(2,2,DDK)- &
                 (1.d0+xi1)*(1.d0+xi2)/4.d0 * DM_Vertex(3,2,DDK)- &
                 (1.d0-xi1)*(1.d0+xi2)/4.d0 * DM_Vertex(4,2,DDK)

        enddo ! enddo i
        
     enddo ! enddo j
     
     
     !-------------------------------------------------------------!
     ! This completes the points on Quadalateron
     !-------------------------------------------------------------!
     
     !=============================================================!
     
  enddo ! enddo DDK

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!  Deallocate variable Edge_x and Edge_y
  if ( allocated( Edge_x ) ) Deallocate( Edge_x, stat=ierr )
  if ( ierr .ne. 0 ) then
     write(*,*) 'Can not deallocate Bnd_x and in Grid'
     write(*,*) 'Abort'
     stop
  endif
  
  if ( allocated( Edge_y ) ) Deallocate( Edge_y, stat=ierr )
  if ( ierr .ne. 0 ) then
     write(*,*) 'Can not deallocate Bnd_y and in Grid'
     write(*,*) 'Abort'
     stop
  endif
  
  return
  
end subroutine Init_Physical_Grid_2D

!======================================================================

subroutine AssignPts(SN,NumDeg,DK,Eg_x,Eg_y)
  use MD2D_Grid
  use Legendre
  use universal_const
  implicit none
  
  integer :: SN
  integer :: v1_num, v2_num
  integer :: NumDeg, DK
  real(kind=8) :: Eg_x(0:NumDeg), Eg_y(0:NumDeg)
  
  integer N, i, j, k
  
  !Linear mapping parameters
  real(kind=8) :: x_start, x_end
  real(kind=8) :: y_start, y_end
  !..Circular arc mapping parameters
  real(kind=8) :: x_cen, y_cen, radius
  real(kind=8) :: v1x, v1y, Lv1, v2x, v2y, Lv2
  real(kind=8) :: angle_ref, cross_ref, cross, angle, t
  real(kind=8) :: focus1_x, focus1_y, focus2_x, focus2_y
  real(kind=8) :: maj_axis, min_axis, dis_focus, rot_angle
  
  select case(SN)
    case(1) !..Side 1: v1-->v2
      v1_num= 1
      v2_num= 2      
    case(2) !..Side 2: v2-->v3
      v1_num= 2
      v2_num= 3      
    case(3) !..Side 3: v4-->v3
      v1_num= 4
      v2_num= 3      
    case(4) !..Side 4: v1-->v4
      v1_num= 1
      v2_num= 4
  end select !SN

     select case(DMBnd_Shape(SN,DK))
     case(1) !..straight side linear transform
     	x_start = DM_Vertex(v1_num,1,DK)
     	y_start = DM_vertex(v1_num,2,DK) 
     
     	x_end   = DM_Vertex(v2_num,1,DK)
     	y_end   = DM_vertex(v2_num,2,DK)
     	
     	     
     	Eg_x(0:NumDeg) = x_start + (LGLCoord(0:NumDeg,NumDeg)+1.d0)/2.d0 * &
        			(x_end-x_start)
     
    	Eg_y(0:NumDeg) = y_start + (LGLCoord(0:NumDeg,NumDeg)+1.d0)/2.d0 * &
        	  		(y_end-y_start)
        	  	
     case(2) !..circular arc
        x_cen  = DMBnd_Shape_Par(1,SN,DK)
        y_cen  = DMBnd_Shape_Par(2,SN,DK) 
        radius = DMBnd_Shape_Par(3,SN,DK)
         
        v1x = DM_Vertex(v1_num,1,DK) - x_cen
        v1y = DM_Vertex(v1_num,2,DK) - y_cen
        Lv1 = dsqrt(v1x**2+v1y**2)
        v2x = DM_Vertex(v2_num,1,DK) - x_cen
        v2y = DM_Vertex(v2_num,2,DK) - y_cen
        Lv2 = dsqrt(v2x**2+v2y**2)
        angle_ref = dacos(v1x/Lv1)
        cross_ref = v1y
        if (cross_ref.lt.0.0d0) angle_ref = 2.0d0*pi - angle_ref
        cross = v1x*v2y - v1y*v2x
        angle = dacos((v1x*v2x+v1y*v2y)/(Lv1*Lv2))
        if (cross.lt.0.0d0) angle=-angle
        if (dabs(angle).lt.1.0d-6) angle=2.0d0*pi

        do i=0,NumDeg
           t  = (LGLCoord(i,NumDeg)+1.d0)/2.d0*angle + angle_ref
           Eg_x(i)= x_cen + radius * dcos(t)
           Eg_y(i)= y_cen + radius * dsin(t)
        end do
        
     case(3) !..elliptic arc
        focus1_x = DMBnd_Shape_Par(1,SN,DK)
        focus1_y = DMBnd_Shape_Par(2,SN,DK)
        focus2_x = DMBnd_Shape_Par(3,SN,DK)
        focus2_y = DMBnd_Shape_Par(4,SN,DK)
        maj_axis = DMBnd_Shape_Par(5,SN,DK)
        min_axis = DMBnd_Shape_Par(6,SN,DK)
        dis_focus = dsqrt( (focus2_x-focus1_x)**2 + (focus2_y-focus1_y)**2 )
        
        rot_angle = dacos( (focus2_x-focus1_x)/dis_focus )
        
        if ( (focus2_y-focus1_y) .lt. 0.0d0) rot_angle = 2.0d0*pi - rot_angle
         
        x_cen  = (focus1_x+focus2_x)/2.d0
        y_cen  = (focus1_y+focus2_y)/2.d0
         
        v1x = DM_Vertex(v1_num,1,DK) - x_cen
        v1y = DM_Vertex(v1_num,2,DK) - y_cen
        Lv1 = dsqrt(v1x**2+v1y**2)
        v2x = DM_Vertex(v2_num,1,DK) - x_cen
        v2y = DM_Vertex(v2_num,2,DK) - y_cen
        Lv2 = dsqrt(v2x**2+v2y**2)
        angle_ref = dacos(v1x/Lv1)
        cross_ref = v1y
        if (cross_ref.lt.0.0d0) angle_ref = 2.0d0*pi - angle_ref
        cross = v1x*v2y - v1y*v2x
        angle = dacos((v1x*v2x+v1y*v2y)/(Lv1*Lv2))
        if (cross.lt.0.0d0) angle=-angle
        if (dabs(angle).lt.1.0d-6) angle=2.0d0*pi

        do i=0,NumDeg
           t = (LGLCoord(i,NumDeg)+1.d0)/2.d0*angle + (angle_ref-rot_angle)
           Eg_x(i)= x_cen + maj_axis*dcos(t)*dcos(rot_angle) + min_axis*dsin(t)*(-dsin(rot_angle))
           Eg_y(i)= y_cen + maj_axis*dcos(t)*dsin(rot_angle) + min_axis*dsin(t)*(dcos(rot_angle))
        end do

     end select !...Edge assignment completed
     
     return
     
end subroutine AssignPts

!===============================================================================
subroutine Init_intept_Grid_2D
  use Legendre
  use MD2D_Grid
  implicit none
  integer::i,j,ierr 

  call alloc_Mem_Leg_Grid(PND1,PND2,PDeg1)

  do j = 0,PDeg1
     do i = 0,ND1 
        Leg_Grid_xi1(i,0:PND1,j) = leg_tb(i,j)
        DLeg_Grid_xi1(i,0:PND1,j) = Dleg_tb(i,j)        
     enddo
     Leg_Grid_xi2(0:PND1,0:PND1,j) = transpose(Leg_Grid_xi1(0:PND1,0:PND1,j))
     DLeg_Grid_xi2(0:PND1,0:PND1,j) = transpose(DLeg_Grid_xi1(0:PND1,0:PND1,j))
  enddo

  do i=0,ND1
     LGLWeights_Grid_xi1(i,0:PND1) = LGLWeights(i,PND1)
  enddo
  LGLWeights_Grid_xi2 = transpose(LGLWeights_Grid_xi1)

end subroutine Init_intept_Grid_2D

