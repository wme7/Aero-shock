PROGRAM INIT
IMPLICIT REAL(8) (A-H,O-Z)
REAL, DIMENSION (100,100,200,200)   :: F
REAL, DIMENSION (200,200)           :: UX, UY, Z, T, R, ET, P, FN, tha,thb,phia,phib
REAL, DIMENSION (200,200)           :: fll,flh,flcd,fl, fhl,fhh,fhcd,fh,gll,glh,glcd,gl,ghl,ghh,ghcd,gh2
REAL, DIMENSION (200)               :: Y, X, C, V 
REAL, DIMENSION (20)                :: GH, W
OPEN (UNIT = 10, FILE = 'ABSCISSAS.TEC', STATUS = 'UNKNOWN')
OPEN (UNIT = 20, FILE = 'XY.TEC', STATUS = 'UNKNOWN')
OPEN (UNIT = 30, FILE = 'INIT.TEC', STATUS = 'UNKNOWN')
OPEN (UNIT = 40, FILE = 'F.TEC', STATUS = 'UNKNOWN')
OPEN (UNIT = 50, FILE = 'RECIPROCALZT.TEC', STATUS = 'UNKNOWN')
OPEN (UNIT = 60, FILE = 'U2D.TEC', STATUS = 'UNKNOWN')
OPEN (UNIT = 70, FILE = 'U2D1dx.TEC', STATUS = 'UNKNOWN')
OPEN (UNIT = 80, FILE = 'U2D1dy.TEC', STATUS = 'UNKNOWN')

NX          = 198
NY          = 199
CFL         = 0.4
OUTTIME     = 0.26
NV = 20
PI = ATAN2(1.,1.)*4.

!
! GAUSS-HERMITE QUADRATURE
!
GH(:) =(/-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738,-2.78880605843,-2.25497400209,-1.73853771212,-1.2340762154,-0.737473728545,-0.245340708301,0.245340708301,0.737473728545,1.2340762154,1.73853771212,2.25497400209,2.78880605843,3.34785456738,3.94476404012,4.60368244955,5.38748089001/)
W(:)  =(/0.898591961453,0.704332961176,0.62227869619,0.575262442852,0.544851742366,0.524080350949,0.509679027117,0.499920871336,0.493843385272,0.490921500667,0.490921500667,0.493843385272,0.499920871336,0.509679027117,0.524080350949,0.544851742366,0.575262442852,0.62227869619,0.704332961176,0.898591961453/)
    DO I = 1, NV
    C(I)    = W(I)
    V(I)    = GH(I)
    WRITE (10,*) C(I), V(I)    
    END DO
!
! INITIAL INPUT
!
    Z1   = 0.078
    UX1  = 0.3578
    UY1  = -0.6197
    T1   = 1.2586
    
    Z2   = 0.1
    UX2  = 0.
    UY2  = 0.
    T2   = 0.5
    
    Z3   = 0.078656
    UX3  = 0.7156
    UY3  = 0.7156
    T3   = 1.2586
    
    Z4   = 0.1
    UX4  = 0.
    UY4  = 0.
    T4   = 0.5
    
    
    DX = 1./DFLOAT(NX-1)
    DY = 1./DFLOAT(NY-1)
    X(1)  = 0
    Y(1)  = 0
    DO I = 2,NX
     X(I) = X(I-1) + DX
    END DO
    DO J = 2, NY
     Y(J) = Y(J-1) + DY   
    END DO
    DO I = 1, NX
    WRITE (20,*) X(I), Y(I)
    END DO

DO I = 1,NX
    DO J = 2,NY
    
    if (x(I).le.(0.1+(y(j)*0.57735))) then
            Z(I,J) = Z1
            UX(I,J)= UX1
            UY(I,J)= UY1
            T(I,J) = T1 
            else
            Z(I,J) = Z2
            UX(I,J)= UX2
            UY(I,J)= UY2
            T(I,J) = T2
    end if 
           
    
    DO K = 1, NV
    DO L = 1, NV
    PP = ( (V(K)-UX(I,J))**2 + (V(L)-UY(I,J))**2 ) / T(I,J)
    F(k,l,i,j) = 1/((EXP(PP)/Z(I,J)) + 1)
    !PPBa = ( (V(K)+UX(I,ny-1))**2 + (V(L)+UY(I,ny-1))**2 ) / T(I,ny-1)
    PPBb = ( (V(K)+UX(I,2))**2 + (V(L)+UY(I,2))**2 ) / T(I,2)
    !F(k,l,i,ny) = 1/((EXP(PPBa)/Z(I,ny-1)) + 1)
    F(k,l,i,1) = 1/((EXP(PPBb)/Z(I,2)) + 1)
    END DO
    END DO
    
    END DO
END DO

DO I = 1, NX
DO J = 2, NY
    SR  = 0
    SUX = 0
    SUY = 0
    SE  = 0
    DO K = 1, NV
    DO L = 1, NV
    SR  = SR + C(K)*C(L) * F(k,l,i,j)
    SUX = SUX + C(K)*C(L) * F(k,l,i,j) * V(K)
    SUY = SUY + C(K)*C(L) * F(k,l,i,j) * V(L)
    SE  = SE + C(K)*C(L) * F(k,l,i,j) * (0.5 * (V(K)*V(K) + V(L)*V(L)))
    END DO
    END DO
    R(I,J)    = SR
    UX(I,J)   = SUX/SR 
    UY(I,J)   = SUY/SR
    ET(I,J)   = SE  
    WRITE (30,*) X(I),Y(J),r(I,J),ux(I,J),uy(I,J),et(I,J)     
END DO
END DO

ITER  = 1
TIME  = 0
ISTOP = 0

    1000 CONTINUE
    
    DT = min(DX,DY) * CFL/V(NV)
    TIME = TIME + DT
    DTDX = DT/DX
    DTDY = DT/DY
IF (TIME .GT. OUTTIME) THEN
    DTCFL = OUTTIME - (TIME - DTCFL)  
    TIME = OUTTIME
    DT = DTCFL/CFL
    DTDX = DTCFL / DX
    DTDY = DTCFL / DY
    ISTOP = 1
END IF

!
! PROJECTION
!
DO I = 1, NX
DO J = 2, NY
    DO K = 1, NV
    DO L = 1, NV
    PP = ( (V(K)-UX(I,J))**2 + (V(L)-UY(I,J))**2 ) / T(I,J)
    F(k,l,i,j) = 1/((EXP(PP)/Z(I,J)) + 1)
   ! PPBa = ( (V(K)+UX(I,ny-1))**2 + (V(L)+UY(I,ny-1))**2 ) / T(I,ny-1)
    PPBb = ( (V(K)+UX(I,2))**2 + (V(L)+UY(I,2))**2 ) / T(I,2)
   ! F(k,l,i,ny) = 1/((EXP(PPBa)/Z(I,ny-1)) + 1)
    F(k,l,i,1) = 1/((EXP(PPBb)/Z(I,2)) + 1)
    END DO
    END DO  
END DO
END DO


DO K = 1, NV
DO L = 1, NV
    vxp = max(v(k),0.)
    vxm = min(v(k),0.) 
    vyp = max(v(l),0.)
    vym = min(v(l),0.) 
    !---------------theta
    do i = 2, nx-1
    do j = 2, ny-1
        if (f(k,l,i+1,j) .eq. f(k,l,i,j)) then
        tha(i,j) =  0
        else
        tha(i,j) = (f(k,l,i-sign(1.,v(k))+1,j) - f(k,l,i-sign(1.,v(k)),j))/(f(k,l,i+1,j)-f(k,l,i,j))
        end if   
        
        if (f(k,l,i,j+1) .eq. f(k,l,i,j)) then
        thb(i,j) =  0
        else
        thb(i,j) = (f(k,l,i,j-sign(1.,v(l))+1) - f(k,l,i,j-sign(1.,v(l))))/(f(k,l,i,j+1)-f(k,l,i,j))
        end if 
        
    end do
    end do
  !  write (*,*) tha(i,j)
  !  pause
    
    do i = 1, nx
    thb(i,1)=1
    thb(i,ny)=1
    end do
    do j = 1, ny
    tha(1,j)=1
    tha(nx,j)=1    
    end do
    !--------------phi
    do i = 1,nx
    do j = 1,ny
     if (tha(i,j) .le. 0) then
     phia(i,j)=0
     else
     phia(i,j)= (abs(tha(i,j))+tha(i,j))/(1+abs(tha(i,j)))
     end if
     if (thb(i,j) .le. 0) then
     phib(i,j)=0
     else
     phib(i,j)= (abs(thb(i,j))+thb(i,j))/(1+abs(thb(i,j)))
     end if    
    end do
    end do
    !-----------------------TVD SCHEME
    
    DO I = 2, NX-1
    DO J = 2, NY-1 
    fll(i,j) = vxp*f(k,l,i-1,j)+vxm*f(k,l,i,j)
    flh(i,j) = 0.5d0*v(k)*(f(k,l,i-1,j)+f(k,l,i,j))-0.5d0*dtdx*v(k)**2*(f(k,l,i,j)-f(k,l,i-1,j)) 
    !flcd(i,j)= 0.5d0*dtdy*(vxm*vym*(f(k,l,i,j+1)-f(k,l,i,j))+vxp*vym*(f(k,l,i-1,j+1)-f(k,l,i-1,j))+vxm*vyp*(f(k,l,i,j)-f(k,l,i,j-1))+vxp*vyp*(f(k,l,i-1,j)-f(k,l,i-1,j-1)))
    fl(i,j) = fll(i,j)+ (phia(i-1,j)*(flh(i,j)-fll(i,j)))
    !+ (phia(i-1,j)*(flh(i,j)-fll(i,j)))- flcd(i,j)
    
    fhl(i,j) = vxp*f(k,l,i,j)+vxm*f(k,l,i+1,j)
    fhh(i,j) = 0.5d0*v(k)*(f(k,l,i,j)+f(k,l,i+1,j))-0.5d0*dtdx*v(k)**2*(f(k,l,i+1,j)-f(k,l,i,j)) 
    !fhcd(i,j)= 0.5d0*dtdy*(vxm*vym*(f(k,l,i+1,j+1)-f(k,l,i+1,j))+vxp*vym*(f(k,l,i,j+1)-f(k,l,i,j))+vxm*vyp*(f(k,l,i+1,j)-f(k,l,i+1,j-1))+vxp*vyp*(f(k,l,i,j)-f(k,l,i,j-1)))
    fh(i,j)  = fhl(i,j)+ (phia(i,j)*(fhh(i,j)-fhl(i,j)))
    !+ (phia(i,j)*(fhh(i,j)-fhl(i,j)))- fhcd(i,j)
    
    gll(i,j) = vyp*f(k,l,i,j-1)+vym*f(k,l,i,j)
    glh(i,j) = 0.5d0*v(l)*(f(k,l,i,j-1)+f(k,l,i,j))-0.5d0*dtdy*v(l)**2*(f(k,l,i,j)-f(k,l,i,j-1))
    !glcd(i,j)= 0.5d0*dtdx*(vxm*vym*(f(k,l,i+1,j)-f(k,l,i,j))+vxm*vyp*(f(k,l,i+1,j-1)-f(k,l,i,j-1))+vxp*vym*(f(k,l,i,j)-f(k,l,i-1,j))+vxp*vyp*(f(k,l,i,j-1)-f(k,l,i-1,j-1)))
    gl(i,j)  = gll(i,j)+ (phib(i,j-1)*(glh(i,j)-gll(i,j)))
    !+ (phib(i,j-1)*(glh(i,j)-gll(i,j)))- glcd(i,j)
    
    ghl(i,j) = vyp*f(k,l,i,j)+vym*f(k,l,i,j+1)
    ghh(i,j) = 0.5d0*v(l)*(f(k,l,i,j)+f(k,l,i,j+1))-0.5d0*dtdy*v(l)**2*(f(k,l,i,j+1)-f(k,l,i,j))
    !ghcd(i,j)= 0.5d0*dtdx*(vxm*vym*(f(k,l,i+1,j+1)-f(k,l,i,j+1))+vxm*vyp*(f(k,l,i+1,j)-f(k,l,i,j))+vxp*vym*(f(k,l,i,j+1)-f(k,l,i-1,j+1))+vxp*vyp*(f(k,l,i,j)-f(k,l,i-1,j)))
    gh2(i,j)  = ghl(i,j)+ (phib(i,j)*(ghh(i,j)-ghl(i,j)))
    !+ (phib(i,j)*(ghh(i,j)-ghl(i,j)))- ghcd(i,j)
    
    fn(i,j) = f(k,l,i,j) - dtdx*(fh(i,j)-fl(i,j)) - dtdy*(gh2(i,j)-gl(i,j))
    f(k,l,i,j) = fn(i,j)
    ! boundary
    !f(k,l,i,1) = fn(i,3)
    !f(k,l,i,ny)= fn(i,ny-2)
      if (x(i).le. 0.1+0.57735+(2.3094*(TIME-10*dt))) then
      pbound = ( (V(K)-UX1)**2 + (V(L)-UY1)**2 ) / T1 
      f(k,l,i,ny)= 1/((EXP(pbound)/Z1) + 1)
      else
      pbound = ( (V(K)-UX2)**2 + (V(L)-UY2)**2 ) / T2 
      f(k,l,i,ny)= 1/((EXP(pbound)/Z2) + 1)
      end if
      
    f(k,l,1,j) =fn(3,j)
    f(k,l,nx,j)=fn(nx-2,j)   
    END DO
    END DO 
 
END DO
END DO

DO I = 1, NX
DO J = 2, NY
    SR  = 0
    SUX = 0
    SUY = 0
    SE  = 0
    DO K = 1, NV
    DO L = 1, NV
    SR  = SR + C(K)*C(L) * F(k,l,i,j)
    SUX = SUX + C(K)*C(L) * F(k,l,i,j) * V(K)
    SUY = SUY + C(K)*C(L) * F(k,l,i,j) * V(L)
    SE  = SE + C(K)*C(L) * F(k,l,i,j) * (0.5 * (V(K)*V(K) + V(L)*V(L)))
    END DO
    END DO
    R(I,J)    = SR
    UX(I,J)   = SUX/SR 
    UY(I,J)   = SUY/SR
    ET(I,J)   = SE  
     
END DO
END DO

DO I = 1, NX
DO J = 2, NY
ZA = 0.01
ZB = 0.99
DO WHILE (ABS(ZA-ZB) .GT. 0.00001)
    GA1 = 0
    GB1 = 0
    GA2 = 0
    GB2 = 0
    DO L = 1, 50
    GA1 = GA1 + (ZA**L) * (-1)**(L-1)/L
    GB1 = GB1 + (ZB**L) * (-1)**(L-1)/L
    GA2 = GA2 + (ZA**L) * (-1)**(L-1)/(L**2)
    GB2 = GB2 + (ZB**L) * (-1)**(L-1)/(L**2)
    END DO
    PSIA = 2*ET(I,J) - (GA2*(R(I,J)/GA1)**2)/PI - R(I,J)*(UX(I,J)*UX(I,J)+UY(I,J)*UY(I,J))
    PSIB = 2*ET(I,J) - (GB2*(R(I,J)/GB1)**2)/PI - R(I,J)*(UX(I,J)*UX(I,J)+UY(I,J)*UY(I,J))
    ZC = (ZA+ZB)/2
    GC1 = 0
    GC2 = 0
    DO L = 1, 50
    GC1 = GC1 + (ZC**L) * (-1)**(L-1)/L
    GC2 = GC2 + (ZC**L) * (-1)**(L-1)/(L**2)
    END DO
    PSIC = 2*ET(I,J) - (GC2*(R(I,J)/GC1)**2)/PI - R(I,J)*(UX(I,J)*UX(I,J)+UY(I,J)*UY(I,J))
    IF ((PSIA*PSIC) .LT. 0) THEN
    ZB = ZC
    ELSE
    ZA = ZC
    END IF
END DO
Z(I,J) = ZC
T(I,J) = R(I,J)/(PI*GC1)
P(I,j) = ET(i,j) - 0.5*r(i,j)*(ux(i,j)**2+uy(i,j)**2)
END DO
END DO    

WRITE (*,*) ITER, R(NX/2,NY/2), F(NV/2,NV/2,NX/2,NY/2)
IF (ISTOP .EQ. 1) GOTO 2000
ITER = ITER + 1
GOTO 1000
PAUSE
2000 CONTINUE
DO I = 1, NX
DO J = 2, NY
WRITE (60,*) X(I), Y(J), R(I,J), p(I,J), ux(I,J)
END DO 
END DO
skw=sqrt(2.)
DO I = 1, NX
WRITE (70,*) X(I), R(I,(0.5*ny))
END DO
DO J = 1, NY-1
WRITE (80,*) Y(J), R((0.95*nx),J), R((0.5*nx),J)
END DO

STOP
END PROGRAM    