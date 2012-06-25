PROGRAM MAIN
IMPLICIT REAL(8) (A-H,O-Z)
REAL, DIMENSION (1000,1000) :: F,FEQ, FFF
REAL, DIMENSION (1000) :: X, C, V, R, RI,P, U, T, TI, ET, EI, FN, FA, FB, FC, AV
REAL, DIMENSION (20)   :: GH, W
REAL, DIMENSION (1000) :: Z, THETA, PSI, VIS,G1
OPEN (UNIT = 10, FILE = 'POSTITER.TEC', STATUS = 'UNKNOWN')
write(10,*) 'variables = "x","r","p","z","t","tau" '
OPEN (UNIT = 20, FILE = 'INITCHECK.TEC', STATUS = 'UNKNOWN')
NX          = 101
NXP1        = NX + 1
GHNC        = 1
CFL         = 0.5
OUTTIME     = 0.1
IT       = 1
IF (GHNC .EQ. 0) THEN
NV = 20
ELSE
NV = 201
END IF
PI = ATAN2(1.,1.)*4.
GH(:) =(/-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738,-2.78880605843,-2.25497400209,-1.73853771212,-1.2340762154,-0.737473728545,-0.245340708301,0.245340708301,0.737473728545,1.2340762154,1.73853771212,2.25497400209,2.78880605843,3.34785456738,3.94476404012,4.60368244955,5.38748089001/)
W(:)  =(/0.898591961453,0.704332961176,0.62227869619,0.575262442852,0.544851742366,0.524080350949,0.509679027117,0.499920871336,0.493843385272,0.490921500667,0.490921500667,0.493843385272,0.499920871336,0.509679027117,0.524080350949,0.544851742366,0.575262442852,0.62227869619,0.704332961176,0.898591961453/)
IF (GHNC .EQ. 1) GO TO 10
    DO I = 1, NV
    C(I)    = W(I)
    V(I)    = GH(I)    
    END DO
    GO TO 30
10  CONTINUE
    V1  = -10
    V2  = 10
    DV  = (V2-V1)/(NV-1.)
    DO I = 1, NV
    V(I) = V1 + DV * (I-1.)
    END DO
    DO I = 2, (NV-1)
    C(I) = 64./45. * DV
    IF(MOD(I,4).EQ.1) C(I) = 28./45. * DV
    IF(MOD(I,4).EQ.3) C(I) = 24./45. * DV
    END DO
    C(1)    = 14./45. * DV
    C(NV)   = C(1)
30  CONTINUE
    
    !UL  = 0.
    !TL  = 4.38385
    !ZL  = 0.2253353
    !UR  = 0.
    !TR  = 8.972544
    !ZR  =  0.1204582
    
    UL  = 0.
    TL  = 4.38385
    ZL  = 0.2253353
    UR  = 0.
    TR  = 8.972544
    ZR  = 0.1204582
    
    
    DX = 1./DFLOAT(NX-1)
    X(1)  = -0.5 *DX
    DO J = 2,NXP1
    X(J) = X(J-1) + DX
    END DO
    
    DO J = 1, NXP1
     IF (X(J).LE.0.5)THEN
     U(J) = UL
     T(J) = TL
     Z(J) = ZL
     ELSE
     U(J) = UR
     T(J) = TR
     Z(J) = ZR
     END IF       

        DO K = 1, NV
         PP       = (V(K)-U(J))*(V(K)-U(J))/T(J)
         f(k,j)   = 1/((EXP(PP)/Z(J)) + IT) 
        END DO       
    END DO
    
    DO J = 1, NXP1
    SR = 0
    SU = 0
    SE = 0
    SAV= 0
    DO K = 1, NV
    SR = SR + C(K) * F(K,J)
    SU = SU + C(K) * F(K,J) * V(K)
    SE = SE + C(K) * F(K,J) * (0.5 * V(K) * V(K))
    SAV = SAV + C(K) * F(K,J) * abs(V(K))
    END DO
    R(J)    = SR
    U(J)    = SU/SR 
    ET(J)   = SE
    AV(j)   = SAV
    WRITE (20,*) X(J), R(J), Z(J), ET(J)
    END DO
ITER  = 1
TIME  = 0
ISTOP = 0

1000 CONTINUE
    DO J = 1, NXP1
     c_kn = 0.001
     tau = 0.001
     
     vis(j) = c_kn / av(j)
    END DO
    
    DT = DX * CFL/V(NV)
    TIME = TIME + DT
    DTDX = DT/DX
    
    IF (TIME .GT. OUTTIME) THEN
    DTCFL = OUTTIME - (TIME - DTCFL)  
    TIME = OUTTIME
    DT = DTCFL/CFL
    DTDX = DTCFL / DX
    ISTOP = 1
    END IF       
    DO J = 1, NXP1
    DO K = 1, NV
        PP      = (V(K)-U(J)) * (V(K)-U(J)) /T(J)
        Feq(K,J)   = 1/((EXP(PP)/Z(J)) + IT )
    END DO    
    END DO

    
    DO K = 1, NV             
             ANU  = V(K)*DTDX 
             DO j = 2,(NXP1-1)
                !IF (F(K,J+1).EQ.F(K,J)) THEN
                THETA(j) = 0.
                !ELSE
                !THETA(j) = (F(K,J-SIGN(1.,ANU)+1)-F(K,J-SIGN(1.,ANU)))/(F(K,J+1)-F(K,J))
             !END IF        
             END DO
             THETA(1) = 1
             THETA(NXP1) = 1
             DO j = 3,(nxp1-2)
                if (THETA(j) .LE. 0) then
                PSI(j) = 0          
                else
                PSI(J) = (ABS(THETA(J))+THETA(J))/(1+ABS(THETA(J)))
                end if
             END DO
            PSI(1) = 1
            PSI(2) = 1
            PSI(NXP1) = 1
            PSI(NXP1-1) = 1
             DO J=3,(NXP1-2)
             FN(J)  =  F(K,J) - DTDX * ( 0.5D0*V(K)*( ( F(K,J) + F(K,J+1) ) - ( F(K,J-1)+F(K,J) ) ) - 0.5D0 * SIGN(1.,ANU)*V(K) * ( ( F(K,J+1)-F(K,J) ) - ( F(K,J)-F(K,J-1) ) )+ 0.5D0 * ( (PSI(J) * (SIGN(1.,ANU)-ANU) * V(K)*(F(K,J+1)-F(K,J)) ) - (PSI(J-1) * (SIGN(1.,ANU)-ANU)*V(K)*( F(K,J)-F(K,J-1) ) ) ) ) + (DT/VIS(J))*(FEQ(K,J)-F(K,J))
             END DO 
            !FN(1) = F(K,1)
            !FN(2) = F(K,2)
            !FN(NXP1) = F(K,NXP1)
            !FN(NXP1-1)= F(K,NXP1-1) 
            DO j = 1,NXP1
            F(K,J) = FN(J)  
            F(K,1) = fn(4)
            F(K,2) = fn(3) 
            F(K,NXP1) = fn(nxp1-3)      
            F(K,NXP1-1) = fn(nxp1-2)        
            END DO

    END DO    
    DO J = 1, NXP1
    SR = 0
    SU = 0
    SE = 0
    SAV= 0
    DO K = 1, NV
    SR = SR + C(K) * F(K,J)
    SU = SU + C(K) * F(K,J) * V(K)
    SE = SE + C(K) * F(K,J) * (0.5 * V(K) * V(K))
    SAV = SAV + C(K) * F(K,J) * abs(V(K))
    END DO
    R(J)    = SR
    U(J)    = (SU/SR)  
    ET(J)   = SE  
    AV(J)   = SAV   
    END DO   
    if (IT.eq.0) go to 1111
    
    do j = 1, nxp1
    za = 0.0001
    zb = 0.99
    do while (abs(za-zb) .gt. 0.00001)
    ga12 = 0
    gb12 = 0
    ga32 = 0
    gb32 = 0
        do l = 1, 50
        if (IT.eq.1) then 
        ga12 = ga12 + (za**l)*(-1)**(l-1)/(l**0.5)
        gb12 = gb12 + (zb**l)*(-1)**(l-1)/(l**0.5)
        ga32 = ga32 + (za**l)*(-1)**(l-1)/(l**1.5)
        gb32 = gb32 + (zb**l)*(-1)**(l-1)/(l**1.5)
        else
        ga12 = ga12 + (za**l)/(l**0.5)
        gb12 = gb12 + (zb**l)/(l**0.5)
        ga32 = ga32 + (za**l)/(l**1.5)
        gb32 = gb32 + (zb**l)/(l**1.5)
        end if
        end do
   psia = 2*et(j) - ga32*(r(j)/ga12)**3/(2*pi) - r(j)*u(j)**2
   psib = 2*et(j) - gb32*(r(j)/gb12)**3/(2*pi) - r(j)*u(j)**2
        zc = (za + zb)/2
        gc12 = 0
        gc32 = 0
        gc52 = 0
        do l = 1, 50
        if (IT.eq.1) then
        gc12 = gc12 + (zc**l)*(-1)**(l-1)/(l**0.5)
        gc32 = gc32 + (zc**l)*(-1)**(l-1)/(l**1.5)
        gc52 = gc52 + (zc**l)*(-1)**(l-1)/(l**2.5)
        else
        gc12 = gc12 + (zc**l)/(l**0.5)
        gc32 = gc32 + (zc**l)/(l**1.5)
        gc52 = gc52 + (zc**l)/(l**2.5)
        end if
        end do
   psic = 2*et(j) - gc32*(r(j)/gc12)**3/(2*pi) - r(j)*u(j)**2
   
        if ((psia*psic) .lt. 0) then
        zb = zc
        else
        za = zc
        end if
   end do
   z(j) = zc
   t(j) = r(j)**2 / (pi*gc12**2 )  
   p(j) = et(j) - 0.5 * r(j) * u(j)**2 
   
   end do  
    go to 1112
1111 continue
      do j = 1, nxp1
      t(j)    = 4*et(j)/r(j) - 2*u(j)*u(j) 
      z(j)    = r(j) / SQRT(PI * t(j))    
      p(j) = et(j) - 0.5 * r(j) * u(j)**2  
      end do  
1112 continue        
        
        WRITE(*,*) ITER, X(NXP1/2), R(NXP1/2), Z(NXP1/2)
        IF (ISTOP .EQ. 1) GO TO 2000
        ITER = ITER + 1
        GO TO 1000
2000 CONTINUE

      DO K = 1, NV
        WRITE (40,*) V(K), F(K,nxp1/3)
      END DO
    DO J = 1, NXP1
    WRITE (10,*) X(J), R(J), p(j), z(j), t(j), vis(j)
    END DO
STOP    
END PROGRAM




