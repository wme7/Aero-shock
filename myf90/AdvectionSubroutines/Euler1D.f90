    MODULE MAINVAR
    IMPLICIT NONE
    real,parameter :: small=1.0e-10, gama=1.4, EPSILON=1.0e-6
    INTEGER :: NN,step, i_step, n_output, &   !! grid number
               rk   !! runge-kutta step
    REAL,ALLOCATABLE :: RO(:), U(:), P(:),     &   !! natural var, conservative var better???
            ROU(:),ROE(:),   &
            ro1(:),rou1(:),roe1(:) , & !! for temporary use in runge-kutta
            F1(:), F2(:),F3(:),    &   !! FLUX
            ROL(:),UL(:),PL(:), ROR(:),UR(:),PR(:)
    real :: dx,dt,deltaT_X

    END MODULE MAINVAR
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! MAIN PROGRAM 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    program main
    use mainvar
    
    call initiate
    do i_step=1,step
      write(*,*) "computing ",i_step,"-th step"
      do rk=1,3
      call ghostcell
      call reconstruction
      call boundarycondition
      call flux_HLL
      call EVOLUTION
      call UpdateNatureVar
      enddo
      if(mod(i_step,n_output)==0)then
        write(*,*) "backup..."
        call output
      endif
    enddo
    call output
    end program

!!  inittiate
    subroutine initiate
    use mainvar
    implicit none
    integer :: i,j
    real :: a,b,c, xx
	!! read init data
    open(unit=101,file="init.dat", status="old")
    read(101,*) NN, dx,dt ,step, n_output
    deltaT_X=dt/dx

    !! allocate 
    allocate(RO(0:NN+1))
    allocate(U(0:NN+1))
    allocate(P(0:NN+1))
    allocate(ROU(0:NN+1))
    allocate(ROE(0:NN+1))
    allocate(RO1(0:NN+1))
    allocate(ROU1(0:NN+1))
    allocate(ROE1(0:NN+1))

    allocate(F1(0:NN))
    allocate(F2(0:NN))
    allocate(F3(0:NN))
    allocate(ROL(0:NN))
    allocate(UL(0:NN))
    allocate(PL(0:NN))
    allocate(ROR(0:NN))
    allocate(UR(0:NN))
    allocate(PR(0:NN))

    !! shock tube
     read(101,*) a,b,c
     do i=1,NN/2 !! left state
       RO(i)=a
       U(i) =b
       P(I) =c
       call NatureToConserv(ro(i),u(i),p(i), rou(i),roe(i))
     enddo
     read(101,*) a,b,c
     do i=NN/2+1, NN !! right state
       RO(i)=a
       U(i) =b
       P(I) =c
       call NatureToConserv(ro(i),u(i),p(i), rou(i),roe(i))
     enddo

	!! shu-osher 
!	NN=200
!	dx=0.05
!	dt=0.001
!	step=1800
!    do i=1,NN
!      xx=(i-0.5)*dx
!      if(xx<1.)then
!        RO(i)= 3.857
!        U(i) = 2.629
!        P(i) = 10.333
!      else
!        RO(i)= 1+0.2*sin(5.*(xx-5.))
!        U(i) = 0.0
!        P(i) = 1.0
!      endif
!      call NatureToConserv(ro(i),u(i),p(i), rou(i),roe(i))
!    enddo

    end subroutine

    subroutine NatureToConserv(tro,tu,tp, trou,troe)
    use mainvar
    real :: tro,tu,tp, trou,troe
    trou=tro*tu
    troe= tro*(0.5*tu*tu+tp/(gama-1.)/tro)
    end subroutine

    subroutine ConservToNature(tro,trou,troe, tu,tp)
    use mainvar
    real :: tro,trou,troe, tu,tp
    tu= trou/tro
    tp= (troe-0.5*tro*tu*tu)*(gama-1.)
    end subroutine

!!  ghost cell, just for the MUSCL reconst,
    subroutine ghostcell
    use mainvar
    implicit none
    RO(0)= RO(1) 
    U(0) = -U(1)
    p(0) = P(1)
    RO(NN+1)= RO(NN)
    U(NN+1) = -U(NN)
    P(NN+1) = P(NN)
    end subroutine ghostcell

    subroutine boundarycondition
    use mainvar
    implicit none
    !! wall bc ???
    ROL(0) = roR(0)
    UL(0)  = uR(0)
    PL(0)  = pR(0)

    ROR(NN)= roL(NN)
    UR(NN) = uL(NN)
    PR(NN) = pL(NN)
    end subroutine boundarycondition

    real function dot(x,y)
    real :: x(3),y(3)
    dot=x(1)*y(1)+x(2)*y(2)+x(3)*y(3)
    return
    end
    subroutine reconstruction
    use mainvar
	integer :: i
	real :: Lv1(3),Lv2(3),Lv3(3), Rv1(3),Rv2(3),Rv3(3), vv(6,3), vl(3),vr(3),uu(3)

    do i=2, NN-2   !! boundary reocnst ????    (i+1/2)+, (i+1/2)-
      tro= 0.5*(ro(i)+ro(i+1))
      tu = 0.5*(u(i) +u(i+1) )
      tp = 0.5*(p(i) +p(i+1) )
      tH=gama/(gama-1.)*tp/tro+ 0.5*tu*tu
      ta = sqrt(gama*tp/tro)
      Lv1=(/tu*tu/2.+tu*ta/(gama-1.), -tu-ta/(gama-1.), 1./) *(gama-1.)/2./ta/ta
      Lv2=(/2*ta*ta/(gama-1.)-tu*tu, 2*tu, -2. /)*(gama-1.)/2./ta/ta
      Lv3=(/tu*tu/2-tu*ta/(gama-1.), -tu+ta/(gama-1.), 1./) *(gama-1.)/2./ta/ta
      Rv1=(/1.,1.,1./)
      Rv2=(/tu-ta,tu,tu+ta/)
      Rv3=(/tH-tu*ta,0.5*tu*tu,tH+tu*ta/)

      !!! 将邻近的6点在左特征向量上投影,??? 应用使用的是守恒变量,,,,,,,,,,,,
      do k=1,6
        uu=(/ ro(i-3+k), rou(i-3+k), roe(i-3+k) /)
        vv(k,1:3)= (/dot(Lv1,uu), dot(Lv2,uu), dot(Lv3,uu) /)
      enddo

      call WenoReconst(vv(1:6,1),v1L,v1R)
      call WenoReconst(vv(1:6,2),v2L,v2R)
      call WenoReconst(vv(1:6,3),v3L,v3R)

      vL=(/v1L,v2L,v3L/)
      roL(i)= dot(Rv1,vL)
      uL(i) = dot(Rv2,vL)/roL(i)
      pL(i) = (dot(Rv3,vL)-0.5*roL(i)*uL(i)*uL(i))*(gama-1.)

      vR=(/v1R,v2R,v3R/)
      roR(i)= dot(Rv1,vR)
      uR(i) = dot(Rv2,vR)/roR(i)
      pR(i) = (dot(Rv3,vR)-0.5*roR(i)*uR(i)*uR(i))*(gama-1.)

      !! 得出守恒变量的重构后，再变回原始变量
    enddo

    roR(0)= ro(1)
    uR(0) = u(1)
    pR(0) = p(1)
    roR(1)= ro(2)
    uR(1) = u(2)
    pR(1) = p(2)
    roL(1)= ro(1)
    uL(1) = u(1)
    pL(1) = p(1)

    roR(NN-1)= ro(NN)
    uR(NN-1) = u(NN)
    pR(NN-1) = p(NN)
    roL(NN-1)= ro(NN-1)
    uL(NN-1) = u(NN-1)
    pL(NN-1) = p(NN-1)
    roL(NN)= ro(NN)
    uL(NN) = u(NN)
    pL(NN) = p(NN)

!      do i=1,NN !! use 1-th order
!      roL(i)= ro(i)
!      uL(i) = u(i)
!      pL(i) = p(i)
! 
!      roR(i-1)= ro(i)
!      uR(i-1) = u(i)
!      pR(i-1) = p(i)
!    enddo
    end subroutine

!!!!!!!!!!!!!!!! FLUX , HLL Flux !!!!!!!!!!!!!!
    subroutine FLUX_HLL
    use mainvar
    implicit none
    integer :: i
    real :: ro_L,u_L,p_L, ro_R,u_R,p_R, al,ar,DL,DR,F1L,F2L,f3l,f1r,f2r,f3r
    do i=0,NN
      !! get left and right state
      ro_L= max(roL(i),small)
      u_L = uL(i)
      p_L = max(pL(i),small)
      ro_R= max(roR(i),small)
      u_R = uR(i)
      p_R = max(pR(i),small)

      al=sqrt(gama*p_L/ro_L)
      ar=sqrt(gama*p_R/ro_R)
      DL=min(u_L-al, u_R-ar)
      DR=max(u_L+al, u_R+ar)

      F1L= ro_L*u_L
      F2L= ro_L*u_L*u_L + p_L
      F3L= u_L*p_L*gama/(gama-1.) + 0.5*ro_L*u_L*u_L*u_L

      F1R= ro_R*u_R
      F2R= ro_R*u_R*u_R + p_R
      F3R= u_R*p_R*gama/(gama-1.) + 0.5*ro_R*u_R*u_R*u_R

      !! compute the flux
	  if(DL>0.) then
	    F1(i)=F1L
		F2(i)=F2L
		F3(i)=F3L
      else if(DR<0) then
	    F1(i)=F1R
		F2(i)=F2R
		F3(i)=F3R
	  else
        F1(i) = ( DR*F1L-DL*F1R+DR*DL*(ro_R-ro_L) )/(DR-DL)
        F2(i) = ( DR*F2L-DL*F2R+DR*DL*(u_R - u_L) )/(DR-DL)
        F3(i) = ( DR*F3L-DL*F3R+DR*DL*(p_R - p_L) )/(DR-DL)
	  endif
    enddo
    end subroutine FLUX_HLL

    subroutine flux_LAX
    use mainvar
    implicit none
	integer :: i
    real :: ro_L,u_L,p_L, ro_R,u_R,p_R, al,ar,DL,DR,F1L,F2L,f3l,f1r,f2r,f3r, rom(3), &
            rou_L,roe_L, rou_R,roe_R

    do i=0,NN
      ro_L= max(roL(i),small)
      u_L = uL(i)
      p_L = max(pL(i),small)
      rou_L=ro_L*u_L
      roe_L=p_L/(gama-1.)+0.5*ro_L*u_L*u_L
      ro_R= max(roR(i),small)
      u_R = uR(i)
      p_R = max(pR(i),small)
      rou_R=ro_R*u_R
      roe_R=p_R/(gama-1.)+0.5*ro_R*u_R*u_R

      al=sqrt(gama*p_L/ro_L)
      ar=sqrt(gama*p_R/ro_R)
      rom(1)=max(abs(u_L-al), abs(u_R-ar))
      rom(2)=max(abs(u_L), abs(u_R))
      rom(3)=max(abs(u_L+al), abs(u_R+ar))

      F1L= ro_L*u_L
      F2L= ro_L*u_L*u_L + p_L
      F3L= u_L*p_L*gama/(gama-1.) + 0.5*ro_L*u_L*u_L*u_L

      F1R= ro_R*u_R
      F2R= ro_R*u_R*u_R + p_R
      F3R= u_R*p_R*gama/(gama-1.) + 0.5*ro_R*u_R*u_R*u_R

      F1(i)=0.5*(F1L+F1R - rom(1)*(ro_R-ro_L))
      F2(i)=0.5*(F2L+F2R - rom(2)*(rou_R-rou_L))
      F3(i)=0.5*(F3L+F3R - rom(3)*(roe_R-roe_L))
    enddo
    end subroutine flux_LAX

!!!!!!!!!!!!! Evolution !!!!!!!!!!!!
    subroutine evolution
    use mainvar
    !! use third-order runge-kutta 
    if(rk==1)then
     do i=1,NN
      ro1(i)  = ro(i)
      rou1(i) = rou(i)
      roe1(i) = roe(i)

      ro(i) = ro(i) + deltaT_X*(F1(i-1)-F1(i))
      rou(i)= rou(i)+ deltaT_X*(F2(i-1)-F2(i))
      roe(i)= roe(i)+ deltaT_X*(F3(i-1)-F3(i))
     enddo
    else if(rk==2)then
     do i=1,NN
      ro(i) = 0.75*ro1(i) + 0.25*(ro(i) + deltaT_X*(F1(i-1)-F1(i)))
      rou(i)= 0.75*rou1(i)+ 0.25*(rou(i)+ deltaT_X*(F2(i-1)-F2(i)))
      roe(i)= 0.75*roe1(i)+ 0.25*(roe(i)+ deltaT_X*(F3(i-1)-F3(i)))
     enddo
    else  !! rk==3
     do i=1,NN
      ro(i) = 0.333333333333*ro1(i) + 0.6666666666667*( ro(i) + deltaT_X*(F1(i-1)-F1(i)) )
      rou(i)= 0.333333333333*rou1(i)+ 0.6666666666667*( rou(i)+ deltaT_X*(F2(i-1)-F2(i)) )
      roe(i)= 0.333333333333*roe1(i)+ 0.6666666666667*( roe(i)+ deltaT_x*(F3(i-1)-F3(i)) )
     enddo
    endif
    end subroutine

    subroutine UpdateNatureVar
    use mainvar
    do i=1,NN
      !! ro(i)
      u(i)= rou(i)/ro(i)
      p(i)=(roe(i)-0.5*ro(i)*u(i)*u(i))*(gama-1.)
    enddo
    end subroutine UpdateNatureVar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! output !!!!!!!!!!!!!!!!!!!!
    subroutine output
    use mainvar
    character*20 outfile
    write(outfile,*) i_step/n_output,".dat"
    open(unit=103,file=outfile)
    write(103,*) "title= shocktube"
    write(103,*) 'variables="x","ro","u","p"'
    write(103,*) 'zone i=',NN," f= POINT"
    do i=1,NN
      xx=dx*(i-0.5)
      write(103,*) xx, ro(i),u(i),p(i)
    enddo
    close(103)
    end subroutine output
 
