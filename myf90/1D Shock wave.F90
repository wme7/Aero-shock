PROGRAM MAIN
IMPLICIT REAL*8 (A-H,O-Z)
REAL, DIMENSION (1000,1000) :: F,FEQ
REAL, DIMENSION (1000) :: X, C, V, R,P, U, T,ET,FN,AV,WX
REAL, DIMENSION (20)   :: GH, W
REAL, DIMENSION (1000) :: Z, THETA, PSI, VIS 
INTEGER (KIND=1) :: MMM
OPEN (UNIT = 10, FILE = 'POSTITER10.0001.plt', STATUS = 'UNKNOWN')
WRITE(10,*) 'VARIABLES = "X","R","P","Z","T","U" '

OPEN (UNIT = 20, FILE = 'INITCHECK.TXT', STATUS = 'UNKNOWN')
NX          = 101
NXP1        = NX + 1
GHNC        = 1
CFL         = 0.05
OUTTIME     = 0.1
TAU			= 0.0001 !RELAXATION TIME
WRITE(10,*) "zone I=",Nx
IT       =0 ! 0 for MB, -1 for BE, +1 FD
IF (GHNC .EQ. 0) THEN
NV = 20
ELSE
NV = 201
END IF
PI = ATAN2(1.,1.)*4.
GH(1:5) =(/-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738,-2.78880605843/)
GH(6:10) =(/-2.25497400209,-1.73853771212,-1.2340762154,-0.737473728545,-0.245340708301/)
GH(11:15) =(/0.245340708301,0.737473728545,1.2340762154,1.73853771212,2.25497400209/)
GH(16:20) =(/2.78880605843,3.34785456738,3.94476404012,4.60368244955,5.38748089001/)
W(1:5)  =(/0.898591961453,0.704332961176,0.62227869619,0.575262442852,0.544851742366/)
W(6:10) =(/0.524080350949,0.509679027117,0.499920871336,0.493843385272,0.490921500667/)
W(11:15) =(/0.490921500667,0.493843385272,0.499920871336,0.509679027117,0.524080350949/)
W(16:20) =(/0.544851742366,0.575262442852,0.62227869619,0.704332961176,0.898591961453/)
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

    !Case 1 Sod's shock tube problem
	! Classical IC
	RL=1.0  ! Density
	UL=0.75
	PL=1.0
	! Semi-classical
	TL  = 4
    ZL  = 0.2821
    
    ! Classical IC
    RR=0.1
	UR=0
	PR=0.125
    ! Semi-classical
    TR  = 3.2
    ZR  = 0.0394

    !Case 2
    !UL  = 0.
    !TL  = 4.38385
    !ZL  = 0.2253353
    !UR  = 0.
    !TR  = 8.972544
    !ZR  = 0.1204582	       
	
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
     P(J) = PL
     R(J) = RL
     WX(J) = PL
     ELSE
     U(J) = UR
     T(J) = TR
     Z(J) = ZR
     P(J) = PR
     R(J) = RR
     WX(J) = PR
     END IF       
     
        DO K = 1, NV
        
        PP      =(V(K)-U(J)) * (V(K)-U(J))/ T(J)
        F(K,J)   = 1/((EXP(PP)/Z(J)) + IT )
         !PP       = (V(K)-U(J))*(V(K)-U(J))/(1*T(J))
         !F(K,J)   = 1/((EXP(PP)/Z(J)) + IT) 
        END DO       
    END DO
    
!!!!!!!!!!!D.O.M
    DO J = 1, NXP1
    SR = 0
    SU = 0
    SE = 0
    SAV= 0
    SW_xx = 0
    DO K = 1, NV
    SR = SR + C(K) * F(K,J)
    SU = SU + C(K) * F(K,J) * V(K)
    SE = SE + C(K) * F(K,J) * (0.5 * V(K) * V(K))
    SAV = SAV + C(K) * F(K,J) * ABS(V(K))
    SW_xx = SW_xx + C(K)*F(K,J)*(V(K)-U(J))*(V(K)-U(J))
    END DO
    R(J)    = SR
    U(J)    = SU/SR 
    ET(J)   = SE
    AV(J)   = SAV
    WX (J)  = SW_xx
    WRITE (20,*) X(J), R(J), Z(J), ET(J)
    END DO
ITER  = 1
TIME  = 0
ISTOP = 0

1000 CONTINUE
    DO J = 1, NXP1
  
     
     VIS(J) = TAU
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
          
        PP      = (V(K)-U(J)) * (V(K)-U(J))/ T(J)
        FEQ(K,J)   = 1/((EXP(PP)/Z(J)) + IT )
        
        !PP      = (V(K)-U(J)) * (V(K)-U(J)) /(1*T(J))
        !FEQ(K,J)   = 1/((EXP(PP)/Z(J)) + IT )
    END DO    
    END DO
    
    
!!!!!!!!!!!!!!T.V.D.     
    DO K = 1, NV             
       ANU  = V(K)*DTDX 
       DO J = 2,(NXP1-1)
          IF (F(K,J+1).EQ.F(K,J)) THEN
          THETA(J) = 0.
          ELSE
          MMM = SIGN(1d0,ANU)
          THETA(J) = (F(K,J-MMM+1)-F(K,J-MMM))/(F(K,J+1)-F(K,J))
       END IF        
       END DO
       THETA(1) = 1
       THETA(NXP1) = 1
       DO J = 3,(NXP1-2)
          IF (THETA(J) .LE. 0) THEN
          PSI(J) = 0          
          ELSE
          PSI(J) = (ABS(THETA(J))+THETA(J))/(1+ABS(THETA(J)))
          END IF
       END DO
      PSI(1) = 1
      PSI(2) = 1
      PSI(NXP1) = 1
      PSI(NXP1-1) = 1
       DO J=3,(NXP1-2)
       FN(J)  =  F(K,J) - DTDX * ( 0.5D0*V(K)*( ( F(K,J) + F(K,J+1) ) - ( F(K,J-1)+F(K,J) ) ) &
                - 0.5D0 * SIGN(1.,ANU)*V(K) * ( ( F(K,J+1)-F(K,J) ) - ( F(K,J)-F(K,J-1) ) ) & 
                + 0.5D0 * ( (PSI(J) * (SIGN(1.,ANU)-ANU) * V(K)*(F(K,J+1)-F(K,J)) ) & 
                - (PSI(J-1) * (SIGN(1.,ANU)-ANU)*V(K)*( F(K,J)-F(K,J-1) ) ) ) ) & 
                + (DT/VIS(J))*(FEQ(K,J)-F(K,J))
       END DO 
      DO J = 1,NXP1
      F(K,J) = FN(J)  
      F(K,1) = FN(4)
      F(K,2) = FN(3) 
      F(K,NXP1) = FN(NXP1-3)      
      F(K,NXP1-1) = FN(NXP1-2)        
      END DO

    END DO    
    DO J = 1, NXP1
    SR = 0
    SU = 0
    SE = 0
    SAV= 0
    SW_xx= 0
    DO K = 1, NV
    SR = SR + C(K) * F(K,J)
    SU = SU + C(K) * F(K,J) * V(K)
    SE = SE + C(K) * F(K,J) * (0.5 * V(K) * V(K))
    SAV = SAV + C(K) * F(K,J) * ABS(V(K))
    SW_xx = SW_xx + C(K)*F(K,J)*(V(K)-U(J))*(V(K)-U(J))
    END DO
    R(J)    = SR
    U(J)    = (SU/SR)  
    ET(J)   = SE  
    AV(J)   = SAV   
    WX(J)    = SW_xx 
    END DO   
    IF (IT.EQ.0) GO TO 1111
    
    ! TEST MYSELF
       
    
    !DO J=1,NXP1
     !ZC = 0.15
     !GA012 = 0
     !GA12 = 0
     !GA32 = 0
     !DO L=1,15
      !IF(IT.EQ.1) THEN
      !GA012 = GA012 + Z(J)**L*(-1)**(L-1)*L**0.5
      !GA12 = GA12 + Z(J)**L*(-1)**(L-1)/L**0.5
      !GA32 = GA31 + Z(J)**L*(-1)**(L-1)/L**1.5
      !ELSE
      !GA012 = GA012 + Z(J)**L*L**0.5
      !GA12 = GA12 + Z(J)**L/L**0.5
      !GA32 = GA32 + Z(J)**L/L**1.5
      !END IF
     !END DO
     !PSIA = 2*PI*P(J)/R(J)**3 - GA32/GA12**3
     !PSIB = 3*GA32*GA012/(Z(J)*GA12**4) - 1/Z(J)*GA12**2
     !PSIA = 2*ET(J) - GA32*(R(J)/GA12)**3/(2*PI) - R(J)*U(J)**2
     !PSIB = (1/(Z(J)*GA12**2) - 3*GA32*GA012/(Z(J)*GA12**4))*(R(J)**3/(2*PI))
     !HA = PSIA/PSIB
     !Z(J) = Z(J) - HA            
     !T(J) = R(J)**2 / (PI*GA12**2 )  
     !P(J) = (ET(J) - 0.5 * R(J) * U(J)**2)
     
     
!!!!!!!!!!!!!!!!!!!!!!!!!! orignal code    
!!!!!! BISECTION METHOD 
     DO J = 1, NXP1
     ZA = 0.0001
     ZB = 0.99
     DO WHILE (ABS(ZA-ZB) .GT. 0.00001)
     GA12 = 0
     GB12 = 0
     GA32 = 0
     GB32 = 0
       DO L = 1, 9
       IF (IT.EQ.1) THEN 
        GA12 = GA12 + (ZA**L)*(-1)**(L-1)/(L**0.5)
        GB12 = GB12 + (ZB**L)*(-1)**(L-1)/(L**0.5)
        GA32 = GA32 + (ZA**L)*(-1)**(L-1)/(L**1.5)
        GB32 = GB32 + (ZB**L)*(-1)**(L-1)/(L**1.5)
        ELSE
        GA12 = GA12 + (ZA**L)/(L**0.5)
        GB12 = GB12 + (ZB**L)/(L**0.5)
        GA32 = GA32 + (ZA**L)/(L**1.5)
        GB32 = GB32 + (ZB**L)/(L**1.5)
        END IF
        END DO
  !PSIA = 2*ET(J) - GA32*(R(J)/GA12)**3/(2*PI) - R(J)*U(J)**2
  !PSIB = 2*ET(J) - GB32*(R(J)/GB12)**3/(2*PI) - R(J)*U(J)**2
   PSIA = 2*PI*WX(J)/R(J)**3 - GA32/(GA12)**3
   PSIB = 2*PI*WX(J)/R(J)**3 - GB32/(GB12)**3
   
        ZC = (ZA + ZB)/2
        GC12 = 0
        GC32 = 0
        GC52 = 0
        DO L = 1, 50
        IF (IT.EQ.1) THEN
        GC12 = GC12 + (ZC**L)*(-1)**(L-1)/(L**0.5)
        GC32 = GC32 + (ZC**L)*(-1)**(L-1)/(L**1.5)
        GC52 = GC52 + (ZC**L)*(-1)**(L-1)/(L**2.5)
        ELSE
        GC12 = GC12 + (ZC**L)/(L**0.5)
        GC32 = GC32 + (ZC**L)/(L**1.5)
        GC52 = GC52 + (ZC**L)/(L**2.5)
        END IF
        END DO
    !PSIC = 2*ET(J) - GC32*(R(J)/GC12)**3/(2*PI) - R(J)*U(J)**2
    PSIC = 2*PI*WX(J)/R(J)**3 - GC32/(GC12)**3
    
        IF ((PSIA*PSIC) .LT. 0) THEN
        ZB = ZC
        ELSE
        ZA = ZC
        END IF
   END DO
   Z(J) = ZC
  !T(J) = (P(J)/R(J))*(GC12/GC32)
   T(J) = R(J)**2 / (PI*GC12**2 )  
   P(J) = (ET(J) - 0.5 * R(J) * U(J)**2)
   
    END DO
    GO TO 1112
1111 CONTINUE
      DO J = 1, NXP1
      T(J)    = 4*ET(J)/R(J) - 2*U(J)*U(J) 
      Z(J)    = R(J) / SQRT(PI * T(J))    
      P(J) = (ET(J) - 0.5 * R(J) * U(J)**2)
      END DO  
1112 CONTINUE        
        
WRITE(*,777) TIME, R(NXP1/2)
777	FORMAT (1X,'ELAPSED TIME:',F7.4,4X, 'DENSITY AT X=4.0,Y=5.:',F7.4)


        IF (ISTOP .EQ. 1) GO TO 2000
        ITER = ITER + 1
        GO TO 1000
2000 CONTINUE
    DO J = 1, NXP1
    WRITE (10,*) X(J), R(J), P(J), Z(J), T(J), U(J) !VIS(J)
    END DO
STOP    
END PROGRAM





