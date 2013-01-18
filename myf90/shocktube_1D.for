C====================================================================C
C                                                                    C
C     FINITE DIFFERENCE METHODS FOR 1-D SHOCK TUBE PROBLEMS          C
C                                                                    C
C            HSU CHIANG-AN       3, AUG.  1989                       C
C                                7, JAN.  1991                       C
C                               26, JUNE  1991                       C
C                                                                    C
C     COMPARATIVE STUDY OF THE STATE-OF-THE-ART HIGH RESOLUTION      C
C     SHOCK-CAPTURE SCHEMES                                          C
C                                                                    C
C            HSU CHIANG-AN      19, JUNE  1992                       C
C                               30, JUNE  1992                       C
C                               20, SEP.  1997                       C
C                                                                    C
C====================================================================C
      PARAMETER( MAXI=407,NV=3 )
      COMMON /A1/ RQ(MAXI),UQ(MAXI),PQ(MAXI),EQ(MAXI),CQ(MAXI),HQ(MAXI)
      DIMENSION UA(MAXI),RA(MAXI),CA(MAXI),HA(MAXI)
      DIMENSION Q(MAXI,NV),F(MAXI,NV),FN(MAXI,NV),X(MAXI)
      DIMENSION EIG(MAXI,NV),ALFA(MAXI,NV),BET(MAXI,NV)
      DIMENSION SGM(MAXI,NV),SIGMB(MAXI,NV),DA(MAXI,NV)
      DIMENSION GA(MAXI,NV),GB(MAXI,NV),DB(MAXI,NV),PSI(MAXI,NV)
      DIMENSION TK(MAXI,NV,NV),TKI(MAXI,NV,NV),SM(MAXI,NV)
      DIMENSION THETA(MAXI,NV),WL(NV),SGMT(MAXI,NV),SGMH(MAXI,NV)
      COMMON /C3/ XEXA(103),PEXA(103),REXA(103),UEXA(103),EEXA(103)
      COMMON /D1/ GAMMA,GM1,EPSI,LEVEL,DTX,IPROB
      CHARACTER*1 IANS
C
C --- INITIAL CONDITION
C
      DATA GAMMA,GM1,X0,WL /1.4, 0.4, 0.5, 0., 0., 0./
C
C     FOR SOD's PROBLEM
C
C     DATA DENSL,UVEL,PREL /1.,0.,1./
C     DATA DENSR,UVER,PRER,TOTIME /0.125,0.,0.1,0.24/
C
C     FOR LAX's PROBLEM
C
C     DATA DENSL,UVEL,PREL /0.445,0.698,3.528/
C     DATA DENSR,UVER,PRER,TOTIME /0.5,0.,0.571,0.12/
C
C     FOR WOODWARD & COLELLA PROBLEM
C
C     DATA DENSL,UVEL,PREL/1.,0.,1000./
C     DATA DENSM,UVEM,PREM/1.,0.,0.01/
C     DATA DENSR,UVER,PRER,TOTIME/1.,0.,100.,0.038/
C
      WRITE(6,*)
      WRITE(6,*) '*************** METHOD OPTION: ***************'
      WRITE(6,*)
      WRITE(6,*) '  (0) 2nd ORDER SYMMETRIC TVD (A CLASS)'
      WRITE(6,*) '  (1) 2nd ORDER UPWIND TVD (A CLASS)'
      WRITE(6,*) '  (2) 2nd ORDER ENO (RD=2 , A CLASS)'
      WRITE(6,*) '  (3) 2nd ORDER UNIFORMLY ACCURACY TVD'
      WRITE(6,*) '  (4) 2nd ORDER TVB'
      WRITE(6,*) '  (5) 3rd ORDER ENO (RP=3 , A CLASS)'
      WRITE(6,*) '  (6) 3rd ORDER TVD (MAKE QUICKEST SCHEME TVD)'
      WRITE(6,*) ' (11) MUSCL-ROE TYPE CHAKRAVARTHY-OSHER TVD (A CLASS)'
      WRITE(6,*) ' (12) MUSCL-ROE TYPE FNO3 (MODIFIED FROM (7) By Pan)'
      WRITE(6,*) ' (13) MUSCL-ROE TYPE ENO2 '
      WRITE(6,*) ' (14) MUSCL-ROE TYPE TVD2 (A CLASS OF SLOPE LIMITER)'
      WRITE(6,*) ' (15) MUSCL-ROE TYPE ENO3 '
      WRITE(6,*) ' (21) MUSCL-LLF TYPE CHAKRAVARTHY-OSHER TVD (A CLASS)'
      WRITE(6,*) ' (22) MUSCL-LLF TYPE FNO3 (MODIFIED FROM (7) By Pan)'
      WRITE(6,*) ' (23) MUSCL-LLF TYPE ENO2 '
      WRITE(6,*) ' (24) MUSCL-LLF TYPE TVD2 (A CLASS OF SLOPE LIMITER)'
      WRITE(6,*) ' (25) MUSCL-LLF TYPE ENO3 '
      WRITE(6,*) ' (31) HANCOCK TYPE CHAKRAVARTHY-OSHER TVD (A CLASS)'
      WRITE(6,*) ' (32) HANCOCK TYPE FNO3 (MODIFIED FROM (7) By Pan)'
      WRITE(6,*) ' (33) HANCOCK TYPE ENO2 '
      WRITE(6,*) ' (34) HANCOCK TYPE TVD2 (A CLASS OF SLOPE LIMITER)'
      WRITE(6,*) ' (35) HANCOCK TYPE ENO3 '
      WRITE(6,*) ' (41) PREDICTOR-CORRECTOR MacCormack (original)'
      WRITE(6,*) ' (42) PREDICTOR-CORRECTOR MacCormack (C. S. Song)'
      WRITE(6,*) ' (43) PREDICTOR-CORRECTOR MacCormack TVD '
      WRITE(6,*) ' (44) DISSIPATIVE TWO-FOUR (Gottlieb & Turkel)'
      WRITE(6,*)
      WRITE(6,'(10H CHOOSE ? ,$)')
      READ (5,*) IMETH
C
      IF( IMETH .LE. 2 .OR. IMETH .EQ. 14 .OR. IMETH .EQ. 43 .OR.
     &    IMETH .EQ. 5 .OR. IMETH .EQ. 24 .OR. IMETH .EQ. 34 ) THEN
          WRITE(6,*) 'INPUT LIMITER PARAMETER : IMOD (1-9) '
          READ (5,*)  IMOD
      END IF
      IF( IMETH .EQ. 4 ) THEN
          WRITE(6,*) 'INPUT COEFFICIENT OF TVB, B (1-3) & M (50-200)'
          READ (5,*)  COEB,COEM
      END IF
      IF( IMETH .EQ. 11 .OR. IMETH .EQ. 21 .OR. IMETH .EQ. 31 ) THEN
          WRITE(6,*) 'INPUT: MUSCL TYPE PARAMETERS : CK '
          WRITE(6,*) '    CK = -1  ==> FULLY UPWIND SCHEME'
          WRITE(6,*) '    CK =  0  ==> Fromm SCHEME'
          WRITE(6,*) '    CK = 1/3 ==> 3rd ORDER UPWIND-BIASED SCHEME'
          WRITE(6,*) '    CK =  1  ==> THREE POINT CENTRAL SCHEME'
          READ (5,*)      CK
          OMIGA = 2.0
      END IF
      IF( IMETH .EQ. 12 .OR. IMETH .EQ. 22 .OR. IMETH .EQ. 32 ) THEN
          CK = 1./3.
          OMIGA= 2.0
      END IF
      IF( IMETH .GT. 10 .AND. IMETH .LE. 30 ) THEN
          WRITE(6,*) 'INPUT: TVD RUNGE-KUTTA INTEGRATION LEVEL '
          WRITE(6,*) '       LEVEL = 1  ==> 1st ORDER TIME ACCURACY'
          WRITE(6,*) '       LEVEL = 2  ==> 2nd ORDER TIME ACCURACY'
          WRITE(6,*) '       LEVEL = 3  ==> 3rd ORDER TIME ACCURACY'
          WRITE(6,*) '       LEVEL = 4  ==> 4rd ORDER TIME ACCURACY'
          READ (5,*)         LEVEL
      END IF
      IF( IMETH .GT. 40 ) LEVEL = 2
      IF( IMETH .GT. 10 .AND. IMETH .LE. 30 ) THEN
          WRITE(6,*) 'INPUT: OPTION FOR LIMITING VARIABLES '
          WRITE(6,*) '       IOPT = 1  ==> CHARACTERISTIC VARIABLES'
          WRITE(6,*) '       IOPT = 2  ==> PRIMITIVE VARIABLES'
          WRITE(6,*) '       IOPT = 3  ==> CONSERVATIVE VARIABLES'
          READ (5,*)         IOPT
      END IF
      IF( IMETH .GE. 1 .AND. IMETH .LE. 4 ) THEN
          WRITE(6,*) 'INPUT A.C.M. COEFFICIENT: W1, W2, W3 (0 - 2)'
          READ (5,*)  WL(1),WL(2),WL(3)
      END IF
      WRITE(6,*) 'INPUT : GRID POINTS (<= 401),  CFL NUMBER (<=1) '
      READ (5,*)  NP,CFL
      WRITE(6,*) 'INPUT : ENTROPY COEFFICIENT (0, 0.05, 0.1)'
      READ (5,*)  EPSI
c     EPSI = 0.01
      WRITE(6,*) 'INPUT  LEFT/RIGHT STATE PROPERTIES'
      WRITE(6,*) '   (1) SOD PROBLEM '
      WRITE(6,*) '   (2) LAX PROBLEM '
      WRITE(6,*) '   (3) USER SPECIFY '
      WRITE(6,*) '   (4) TWO BLAST WAVE INTERACTION PROBLEM '
      WRITE(6,'(16H CHOOSE (1-4) ? ,$)')
      READ (5,*) IPROB
      IF( IPROB .EQ. 3 ) THEN
          WRITE(6,*) 'INPUT : PRER,UVER,DENSR,PREL,UVEL,DENSL,TOTIME '
          READ (5,*)          PRER,UVER,DENSR,PREL,UVEL,DENSL,TOTIME
      ELSEIF( IPROB .EQ. 2 ) THEN
          PREL  = 3.528
          UVEL  = 0.698
          DENSL = 0.445
          PRER  = 0.571
          UVER  = 0.0
          DENSR = 0.5
          TOTIME= 0.12
      ELSEIF( IPROB .EQ. 1 ) THEN
          PREL  = 1.
          UVEL  = 0.
          DENSL = 1.
          PRER  = 0.1
          UVER  = 0.
          DENSR = 0.125
          TOTIME= 0.24
      ELSE
          PREL  = 1000.
          UVEL  = 0.
          DENSL = 1.
          PRER  = 100.
          UVER  = 0.
          DENSR = 1.
          PREM  = 0.01
          UVEM  = 0.
          DENSM = 1.
          TOTIME= 0.038
      END IF
C
      DX  = 1./(NP-1)
      NP  = NP + 3
      NP1 = NP + 1
      NP2 = NP + 2
      NP3 = NP + 3
      W123= WL(1)+WL(2)+WL(3)
C
C     TIME1 = CPUTIME(0.0)
      CALL CLOCK@(TIME1)
      DO 10 I = 1,NP3
      X(I) = DX*(I-4)
   10 CONTINUE
      IF( IPROB .LE. 3 ) THEN
          DO 15 I = 1,NP3
          IF( X(I) .LE. X0 ) THEN
              RQ(I)= DENSL
              UQ(I)= UVEL
              PQ(I)= PREL
          ELSE
              RQ(I)= DENSR
              UQ(I)= UVER
              PQ(I)= PRER
          END IF
   15     CONTINUE
      ELSE
          DO 20 I = 1,NP3
          IF( X(I) .LE. 0.1 ) THEN
              RQ(I)= DENSL
              UQ(I)= UVEL
              PQ(I)= PREL
          ELSEIF( X(I) .LT. 0.9 .AND. X(I) .GT. 0.1 ) THEN
              RQ(I)= DENSM
              UQ(I)= UVEM
              PQ(I)= PREM
          ELSE
              RQ(I)= DENSR
              UQ(I)= UVER
              PQ(I)= PRER
          END IF
   20     CONTINUE
      END IF
C
      DO I = 1,3
      RQ(4-I) = RQ(4+I)
      PQ(4-I) = PQ(4+I)
      RQ(NP+I)= RQ(NP-I)
      PQ(NP+I)= PQ(NP-I)
      IF( IPROB .GT. 3 ) THEN    ! REFLECTION BC
          UQ(4-I) =-UQ(4+I)
          UQ(NP+I)=-UQ(NP-I)
      ELSE                       ! NO BC
          UQ(4-I) = UQ(4+I)
          UQ(NP+I)= UQ(NP-I)
      END IF
      END DO
C
C     CALCULATE E AND H AND C
C
      DO 30 I = 1,NP3
      EQ(I)   = PQ(I)/GM1/RQ(I)
      CQ(I)   = SQRT( GAMMA*PQ(I)/RQ(I) )
      HQ(I)   = GAMMA*PQ(I)/GM1/RQ(I) + 0.5*UQ(I)**2
   30 CONTINUE
C
C     FILL Q COLUMN VECTOR
C
      DO 40 I = 1,NP3
      Q(I,1)  = RQ(I)
      Q(I,2)  = RQ(I)*UQ(I)
      Q(I,3)  = RQ(I)*EQ(I) + 0.5*RQ(I)*UQ(I)*UQ(I)
   40 CONTINUE
C
      DO 45 K = 1,NV
      DO 45 I = 1,NP3
      PSI(I,K)= 1.
   45 CONTINUE
C
      TIME  = 0.
      ITER  = 0
      FAC_LF= 1.1
C
C --- BEGIN ITERATION ---
C
 1000 CONTINUE
C
C     CALCULATE DT
C
      SIGMAX = 0.
      DO 50 I = 1,NP3
      SIGMAX = MAX( SIGMAX , ABS(UQ(I)) + CQ(I) )
   50 CONTINUE
      DT  = CFL*DX/SIGMAX
      ITER= ITER + 1
      TIME= TIME + DT
      IF( TIME .GT. TOTIME ) THEN
          DT  = TOTIME - ( TIME - DT )
          TIME= TOTIME
      END IF
      DTX = DT/DX
C
      IF( IMETH .GT. 10 .AND. IMETH .LE. 30 .OR. IMETH .GT. 40 ) THEN
          CALL RUNGE( 1,NP3,CK,OMIGA,Q,IMETH,IMOD,IOPT )
          GOTO 240
      END IF
      IF( IMETH .GT. 30 .AND. IMETH .LE. 40 ) THEN
          CALL HANCOCK( 1,NP3,CK,OMIGA,Q,IMETH,IMOD )
          GOTO 210
      END IF
C
C     CALCULATE rho, u, p,... at i+1/2
C
C --- by ROE'S AVERAGE ---
C
      DO 55 I = 1,NP2
      DD   = SQRT( RQ(I+1)/RQ(I) )
      RA(I)= DD*RQ(I)
      UA(I)= (DD*UQ(I+1)+UQ(I))/(1.+DD)
      HA(I)= (DD*HQ(I+1)+HQ(I))/(1.+DD)
      CA(I)= SQRT( GM1*(HA(I)-0.5*UA(I)**2) )
   55 CONTINUE
C
C     FILL FLUX VECTOR
C
      DO 60 I= 1,NP3
      F(I,1) = RQ(I)*UQ(I)
      F(I,2) = RQ(I)*UQ(I)*UQ(I) + PQ(I)
      F(I,3) = UQ(I)*( Q(I,3) + PQ(I) )
   60 CONTINUE
C
C     CALCULATE EIGENVALUE AT I+1/2
C
      DO 65 I = 1,NP2
      EIG(I,1)= UA(I)
      EIG(I,2)= UA(I) + CA(I)
      EIG(I,3)= UA(I) - CA(I)
   65 CONTINUE
C
C     CALCULATE AFAR  AT I+1/2
C
      CALL RLVECT( MAXI,NP2,UA,CA,HA,RA,TK,TKI )
C
      DO 70 I = 1,NP2
      Q1 = Q(I+1,1) - Q(I,1)
      Q2 = Q(I+1,2) - Q(I,2)
      Q3 = Q(I+1,3) - Q(I,3)
      ALFA(I,1) = TKI(I,1,1)*Q1 + TKI(I,1,2)*Q2 + TKI(I,1,3)*Q3
      ALFA(I,2) = TKI(I,2,1)*Q1 + TKI(I,2,2)*Q2 + TKI(I,2,3)*Q3
      ALFA(I,3) = TKI(I,3,1)*Q1 + TKI(I,3,2)*Q2 + TKI(I,3,3)*Q3
   70 CONTINUE
C
      IF( IMETH .EQ. 0 ) THEN
C
C     2nd ORDER SYMMETRIC TVD
C
          CALL QFUNCT(IMOD,1,NP2,ALFA,GB)
C
          DO 75 K = 1,NV
          DO 75 I = 1,NP2
          OO      = PHI(EIG(I,K),EPSI)
          O2      = DTX*EIG(I,K)*EIG(I,K)
          BET(I,K)= -( (O2-OO)*GB(I,K) + OO*ALFA(I,K) )
   75     CONTINUE
      END IF
C
      IF( IMETH .LE. 5 .AND. IMETH .GE. 1 ) THEN
C
C     CALCULATE SGM
C
          DO 78 K = 1,NV
          DO 78 I = 1,NP2
          Z       = EIG(I,K)
          SGM(I,K)= 0.5*(ABS(Z)-DTX*Z*Z)
   78     CONTINUE
      END IF
C
      IF( IMETH .EQ. 1 .OR. IMETH .EQ. 2 ) THEN
C
C     2nd ORDER UPWIND TVD/ENO
C
          CALL GFUNCT(IMOD,1,NP2,ALFA,GB,IMETH)
C
      END IF
C
      IF( IMETH .EQ. 3 ) THEN
C
C     2nd ORDER UNIFORMLY UPWIND TVD
C
          DO 81 K=1,NV
          DO 80 I=2,NP2
          Z1 = ALFA( I ,K)
          Z2 = ALFA(I-1,K)
          IF( Z1*Z2 .GE. 0. ) THEN
              GB(I,K)= Z1
              IF( ABS(Z2) .LE. ABS(Z1) ) GB(I,K)= Z2
          ELSE
              GB(I,K)= -Z1
          END IF
   80     CONTINUE
          GB( 1 ,K)= GB(  2,K)
          GB(NP3,K)= GB(NP2,K)
   81     CONTINUE
      END IF
C
      IF( IMETH .EQ. 4 ) THEN
C
C     2nd ORDER UPWIND TVB
C
          DO 83 K=1,NV
          DO 82 I=2,NP2
          Z1 = ALFA( I ,K)
          Z2 = ALFA(I-1,K)
          D1 = COEB*Z2 + COEM*DX*DX*SGN(Z1)
          D2 = COEB*Z1 + COEM*DX*DX*SGN(Z2)
          Y1 = 0.5*( SGN(Z1)+SGN(D1) )*MIN( ABS(Z1),ABS(D1) )
          Y2 = 0.5*( SGN(Z2)+SGN(D2) )*MIN( ABS(Z2),ABS(D2) )
          GB(I,K) = 0.5*( Y1 + Y2 )
   82     CONTINUE
          GB( 1 ,K)= GB(  2,K)
          GB(NP3,K)= GB(NP2,K)
   83     CONTINUE
      END IF
C
      IF( IMETH .LE. 4 .AND. IMETH .GE. 1 ) THEN
C
C     ARTIFICIAL COMPRESSION METHOD (OPTION)
C
          IF( W123 .EQ. 0. ) GOTO 88
          DO 86 K=1,NV
          DO 85 I=2,NP2
          T1= ABS( ALFA(I,K) - ALFA(I-1,K) )
          T2= ABS( ALFA(I,K) ) + ABS( ALFA(I-1,K) )
          IF( T2 .EQ. 0. ) THEN
              THETA(I,K)= 0.
          ELSE
              THETA(I,K)= T1/T2
          ENDIF
   85     CONTINUE
C
C     B.C. for THETA (zero order extrapolation)
C
          THETA( 1 ,K)= THETA(  2,K)
          THETA(NP3,K)= THETA(NP2,K)
   86     CONTINUE
C
          DO 87 K = 1,NV
          DO 87 I = 1,NP2
          PSI(I,K)= 1.+WL(K)*MAX( THETA(I,K),THETA(I+1,K) )
   87     CONTINUE
   88     CONTINUE
C
C     CALCULATE GA FOR TVD/TVB/ENO ( SECOND CHARACTERISTIC SPEED)
C
          DO 90 K = 1,NV
          DO 90 I = 1,NP2
          IF( ALFA(I,K) .EQ. 0. ) THEN
              GA(I,K)= 0.
          ELSE
              GA(I,K)= PSI(I,K)*SGM(I,K)*(GB(I+1,K)-GB(I,K))/ALFA(I,K)
          END IF
   90     CONTINUE
C
          DO 92 K = 1,NV
          DO 92 I = 1,NP2
          BET(I,K)= PSI(I,K)*SGM(I,K)*( GB(I,K)+GB(I+1,K) )
     &            - PHI( EIG(I,K)+GA(I,K),EPSI )*ALFA(I,K)
   92     CONTINUE
      END IF
C
      IF( IMETH .EQ. 5 ) THEN
C
C     THIRD ORDER ENO
C
C                         n          _____     n
C     CALCULATE SIGMA( EIG     )  &  SIGMA( EIG     )
C                        J+1/2                j+1/2
C
          DO 105 K = 1,NV
          DO 105 I = 1,NP2
          Z  = EIG(I,K)
          AZ = ABS(Z)
          T1 = ( DTX*DTX*Z*Z - 3.*DTX*AZ+ 2. )/6.
          T2 = ( DTX*DTX*Z*Z - 1. )/6.
          IF( ABS(ALFA(I-1,K)) .LE. ABS(ALFA(I,K)) ) THEN
              SIGMB(I,K)= 0.5*( (Z+AZ)*T1 + (Z-AZ)*T2 )
          ELSE
              SIGMB(I,K)= 0.5*( (Z+AZ)*T2 + (Z-AZ)*T1 )
          END IF
  105     CONTINUE
C
C                  n                 n           n
C     CALCULATE  GB    = MINMOD{ ALFA      , ALFA      }
C                  j                 j-1/2       j+1/2
C
          CALL GFUNCT(IMOD,1,NP2,ALFA,GB,IMETH)
C
          DO 111 K = 1,NV
          DO 110 I = 3,NP1
          P1= ALFA(I+1,K) - ALFA( I ,K)
          P2= ALFA( I ,K) - ALFA(I-1,K)
          P4= ALFA(I-1,K) - ALFA(I-2,K)
C
C...  The Following Statements for Third Order  ...
C
C                  n                    n              n
C     CALCULATE  DB    = MINMOD{ D  ALFA      , D  ALFA      }
C                  j              -     j-1/2    +     j-1/2
C
C        OR
C                  n                    n              n
C                DB    = MINMOD{ D  ALFA      , D  ALFA      }
C                  j              -     j+1/2    +     j+1/2
C
          IF( ABS(ALFA(I-1,K)) .LE. ABS(ALFA(I,K)) ) THEN
C             DB(I,K)= 0.5*(SGN(P4)+SGN(P2))*MIN(ABS(P4),ABS(P2))
              DB(I,K)= (P4*P2*(P2+P4)+1.E-6*(P4+P2))/(P4*P4+P2*P2+2.E-6)
C             DB(I,K)= 0.
C             IF( P4+P2 .NE. 0. ) DB(I,K) = (P4*P2+ABS(P4*P2))/(P4+P2)
          ELSE
C             DB(I,K)= 0.5*(SGN(P2)+SGN(P1))*MIN(ABS(P2),ABS(P1))
              DB(I,K)= (P1*P2*(P2+P1)+1.E-6*(P1+P2))/(P1*P1+P2*P2+2.E-6)
C             DB(I,K) = 0.
C             IF( P1+P2 .NE. 0. ) DB(I,K) = (P1*P2+ABS(P1*P2))/(P1+P2)
          END IF
  110     CONTINUE
C
C     B.C. for  DB (zero order extrapolation)
C
          DB(2,K)  = DB(3,K)
          DB(1,K)  = DB(2,K)
          DB(NP2,K)= DB(NP1,K)
          DB(NP3,K)= DB(NP2,K)
  111     CONTINUE
C
C     CALCULATE GA & DA (2nd & 3rd CHARACTERISTIC SPEEDS) AT J+1/2
C
          DO 120 K = 1,NV
          DO 120 I = 1,NP2
          IF( ALFA(I,K) .EQ. 0. ) THEN
              GA(I,K)= 0.
              DA(I,K)= 0.
          ELSE
              GA(I,K)=   SGM(I,K)*(GB(I+1,K)-GB(I,K))/ALFA(I,K)
              DA(I,K)= SIGMB(I,K)*(DB(I+1,K)-DB(I,K))/ALFA(I,K)
          END IF
  120     CONTINUE
C
          DO 125 K = 1,NV
          DO 125 I = 1,NP2
                 Z = EIG(I,K) + GA(I,K) + DA(I,K)
          BET(I,K) =     SGM(I,K)*( GB(I,K)+GB(I+1,K) )
     &               + SIGMB(I,K)*( DB(I,K)+DB(I+1,K) )
     &               - PHI(Z,EPSI)*ALFA(I,K)
  125     CONTINUE
      END IF
C
      IF( IMETH .EQ. 6 ) THEN
C
C     THIRD ORDER TVD ( MODIFIED FROM QUICKEST SCHEME )
C
C     CALCULATE SGM, tilde \sigma & \hat \sigma ..
C
          DO 150 K = 1,NV
          DO 150 I = 1,NP2
          Z        = EIG(I,K)
          OO       = ABS(Z)
          SGM(I,K) = 0.5*(PHI(Z,EPSI)-DTX*Z*Z)
          SGMT(I,K)= OO*(DTX**2*Z*Z - 3.*DTX*OO + 2.)/6.
          SGMH(I,K)= OO*(1. - DTX**2*Z*Z)/6.
  150     CONTINUE
C
C     VAN LEER LIMITER
C
          DO 155 K = 1,NV
          DO 155 I = 2,NP
          Z1= ABS(ALFA( I ,K))
          Z2= ABS(ALFA(I-1,K))
          ABC = Z1 - Z2
          DEF = Z1 + Z2
C
          IF( DEF .EQ. 0. ) THEN
              SM(I,K) = 0.
          ELSE
              SM(I,K) = ABC/DEF
          END IF
  155     CONTINUE
C
C     CALCULATE GBAR
C
          DO 161 K = 1,NV
          DO 160 I = 2,NP2
          GB(I,K)  = SGMT( I ,K)*(1. - SM(I,K))*ALFA( I ,K)
     &             + SGMH(I-1,K)*(1. + SM(I,K))*ALFA(I-1,K)
C
  160     CONTINUE
C
C     B.C. for GBAR (zero order extrapolation)
C
          GB( 1 ,K)= GB(  2,K)
          GB(NP3,K)= GB(NP2,K)
  161     CONTINUE
C
C     CALCULATE GA ( SECOND CHARACTERISTIC SPEED)
C
          DO 170 K = 1,NV
          DO 170 I = 1,NP2
          IF( ALFA(I,K) .EQ. 0. ) THEN
              GA(I,K)= 0.
          ELSE
              GA(I,K)= ( GB(I+1,K)-GB(I,K) )/ALFA(I,K)
          END IF
  170     CONTINUE
C
          DO 180 K = 1,NV
          DO 180 I = 1,NP2
          BET(I,K) = GB(I,K) + GB(I+1,K)
     &             - PHI( EIG(I,K)+GA(I,K),EPSI )*ALFA(I,K)
  180     CONTINUE
      END IF
C
C     CALCULATE NUMERICAL FLUX
C
      DO 190 K= 1,NV
      DO 190 I= 1,NP2
      CCC     = TK(I,K,1)*BET(I,1) + TK(I,K,2)*BET(I,2)
     &        + TK(I,K,3)*BET(I,3)
      FN(I,K) = 0.5*( F(I+1,K) + F(I,K) + CCC )
  190 CONTINUE
C
C     UPDATE
C
      DO 200 K= 1,NV
      DO 200 I= 2,NP2
      Q(I,K)  = Q(I,K) - DTX*( FN(I,K)-FN(I-1,K) )
  200 CONTINUE
C
  210 CALL PRIMER( 1,NP3,Q )
C
  240 CONTINUE
C
      IF( MOD(ITER,10) .NE. 0 ) THEN
          WRITE(6,'(13HITERATION # =,I4,5X,7HCFL # =,F6.4,5X,6HTIME =,
     &          F9.5,5X,4HDT =,F8.6)') ITER,CFL,TIME,DT
      END IF
      IF( TIME .LT. TOTIME ) GOTO 1000
C
C     TIME2  = CPUTIME(0.)
      CALL CLOCK@(TIME2)
      TIMECPU= TIME2 - TIME1
      WRITE(6,'(//5X,16HTOTAL CPUTIME = ,F8.2)') TIMECPU
C
      IF( IPROB .EQ. 4 ) THEN
          OPEN (4,FILE='BLAST.DAT',STATUS='UNKNOWN')
          WRITE(4,'(19HTITLE = "BLAST.PLT")')
          WRITE(4,'(21HVARIABLES = X R U P E)')
          WRITE(4,'(20HZONE T = "CALCU", I=,I3, 14H, J=1, F=POINT)')NP-3
          DO I = 4,NP
          WRITE(4,'(5E15.6)') X(I),RQ(I),UQ(I),PQ(I),EQ(I)
          ENDDO
          STOP
      END IF
C
      CALL EXACT(PRER,UVER,DENSR,PREL,UVEL,DENSL,TOTIME)
      OPEN(10,FILE='fort.10',STATUS='UNKNOWN')
      OPEN(11,FILE='fort.11',STATUS='UNKNOWN')
      DO 260 I = 4,NP
      WRITE(10,'(5E15.6)') X(I),PQ(I),RQ(I),UQ(I),EQ(I)
  260 CONTINUE
      DO 270 I = 1,101
      WRITE(11,'(5E15.6)') XEXA(I),PEXA(I),REXA(I),UEXA(I),EEXA(I)
  270 CONTINUE
C
C     POSTPROCESSOR
C
      WRITE(6,'(/,27H PLOT THE RESULTS ?  (Y/N) ,$)')
      READ (5,'(A1)') IANS
      IF( IANS .NE. 'N' .and. IANS .NE. 'n') THEN
          CALL GRAPH( X(4),RQ(4),UQ(4),PQ(4),EQ(4),NP1-4)
      END IF
      STOP
      END
C
      SUBROUTINE RLVECT(MM,N,UA,CA,HA,RA,TK,TKI)
      DIMENSION UA(MM),CA(MM),HA(MM),RA(MM),TK(MM,3,3),TKI(MM,3,3)
      GAMMA= 1.4
      GM1  = GAMMA-1.
C
      DO 1000 I = 1,N
      B0        = 1./CA(I)
      B1        = 0.5*RA(I)*B0
      B2        = GM1*B0*B0
      B3        = 0.5*UA(I)*UA(I)
      B5        = 1./RA(I)
      B4        = GM1*B5*B0
      TK(I,1,1) = 1.
      TK(I,1,2) = B1
      TK(I,1,3) =-B1
      TK(I,2,1) = UA(I)
      TK(I,2,2) = ( UA(I) + CA(I) )*B1
      TK(I,2,3) =-( UA(I) - CA(I) )*B1
      TK(I,3,1) = B3
      TK(I,3,2) = ( HA(I) + UA(I)*CA(I) )*B1
      TK(I,3,3) =-( HA(I) - UA(I)*CA(I) )*B1
      TKI(I,1,1)= 1. - B2*B3
      TKI(I,1,2)= B2*UA(I)
      TKI(I,1,3)=-B2
      TKI(I,2,1)= B3*B4 - UA(I)*B5
      TKI(I,2,2)= B5 - B4*UA(I)
      TKI(I,2,3)= B4
      TKI(I,3,1)=-B3*B4 - UA(I)*B5
      TKI(I,3,2)= B5 + B4*UA(I)
      TKI(I,3,3)=-B4
 1000 CONTINUE
      RETURN
      END
C
      FUNCTION PHI(Z,EPSI)
C
C     FOR ENTROPY CONDITION
C
      IF( ABS(Z) .GE. EPSI ) THEN
          PHI= ABS(Z)
      ELSE
          PHI= (Z*Z+EPSI**2)/(2.*EPSI)
      END IF
      RETURN
      END
C
      SUBROUTINE EXACT(PRER,UVER,DENSR,PREL,UVEL,DENSL,TIME)
      PARAMETER (IMAX=101,XB=0.,XE=1.,DX=(XE-XB)/(IMAX-1),X0=0.5)
      COMMON /C1/ PR,UR,DR,AR,PL,UL,DL,AL,G,TOL,IERROR
      COMMON /C3/ X(103),P(103),D(103),U(103),E(103)
      DATA G,TOL,IERROR /1.4, 1.E-6, 0/
C
      PR= PRER
      UR= UVER
      DR= DENSR
      PL= PREL
      UL= UVEL
      DL= DENSL
      PI= 4.*ATAN(1.0)
      CALL RIEMANN(PS,US)
      IF( IERROR .EQ. 1 ) THEN
          WRITE(6,'(/,21H ERROR IN INTERATION.)')
      ELSEIF( IERROR .EQ. 2 ) THEN
          WRITE(6,'(/,17H NOT CONVERGE !! )')
      END IF
      IF( US .EQ. UR .AND. US .EQ. UL ) STOP
C...  RCS CASE
      IF( US .GT. UL .AND. US .GE. UR ) CALL RCS(TIME,PS,US)
C...  SCR CASE
      IF( US .LE. UL .AND. US .LT. UR ) CALL SCR(TIME,PS,US)
C...  SCS CASE
      IF( US .LE. UL .AND. US .GE. UR ) CALL SCS(TIME,PS,US)
C...  RCR CASE
      IF( US .GT. UL .AND. US .LT. UR ) CALL RCR(TIME,PS,US)
      RETURN
      END
C
      SUBROUTINE RIEMANN(PS,US)
      COMMON /C1/ PR,UR,DR,AR,PL,UL,DL,AL,G,TOL,IERROR
      COMMON /C2/ GP1,GM1,CL,CR,T1,T2,T3L,T3R,T4,T5
      IF( PR .EQ. PL .AND. UR .EQ. UL ) THEN
          PS = PR
          US = UR
      ELSE
          AR = SQRT(G*PR/DR)
          AL = SQRT(G*PL/DL)
          GP1= G+1.
          GM1= G-1.
          T1 = 2./GM1
          T2 = 1./G/T1
          ZZ = AR/AL*(PL/PR)**T2
          U1 = UL+T1*AL
          U2 = UR-T1*AR
          UI = (U1*ZZ+U2)/(1.+ZZ)
          CALL ITERATE(UI,PS,US)
      END IF
      RETURN
      END
C
      SUBROUTINE ITERATE(U1,PS,US)
      COMMON /C1/ PR,UR,DR,AR,PL,UL,DL,AL,G,TOL,IERROR
      COMMON /C2/ GP1,GM1,CL,CR,T1,T2,T3L,T3R,T4,T5
C
      T3L= GP1/4./AL
      T3R= GP1/4./AR
      T4 = GM1/2.
      T5 = 1./T2
      CL = G*PL/AL
      CR = G*PR/AR
C
      ITER = 0
C
  10  CONTINUE
      ITER = ITER+1
      IF( ITER .GT. 30 ) THEN
          IERROR = 2
          RETURN
      END IF
C
      CALL PLR(U1,PSR,PSL,PSDR,PSDL)
      FUN = PSL - PSR
      FUND= PSDL - PSDR
      IF( IERROR .EQ. 1 ) GOTO 20
      UNEW= U1-FUN/FUND
      U0  = U1
      U1  = UNEW
      WRITE(6,'(/5X,6HITER= ,I2,3X,6HOLDU= ,F11.6,3X,6HNEWU= ,F11.6)')
     &      ITER,U0,U1
      IF( ABS((PSR-PSL)/PSR) .GT. TOL ) GOTO 10
  20  CONTINUE
      PS = PSR
      US = U1
      RETURN
      END
C
      SUBROUTINE PLR(U1,PSR,PSL,PSDR,PSDL)
      COMMON /C1/ PR,UR,DR,AR,PL,UL,DL,AL,G,TOL,IERROR
      COMMON /C2/ GP1,GM1,CL,CR,T1,T2,T3L,T3R,T4,T5
C
      IF( U1 .LE. UL ) THEN
          TM  = T3L*(U1-UL)
          WL  = TM - SQRT(1+TM*TM)
          PSL = PL + CL*WL*(U1-UL)
          TT  = WL*WL
          PSDL= 2.*CL*TT*WL/(1.+TT)
      ELSE
          ASL = AL - T4*(U1-UL)
          PSL = PL*(ASL/AL)**T5
          PSDL= -G*PSL/ASL
      END IF
C
      IF( U1 .GE. UR ) THEN
          TM  = T3R*(U1-UR)
          WR  = TM + SQRT(1+TM*TM)
          PSR = PR + CR*WR*(U1-UR)
          TT  = WR*WR
          PSDR= 2.*CR*TT*WR/(1.+TT)
      ELSE
          ASR = AR + T4*(U1-UR)
          PSR = PR*(ASR/AR)**T5
          PSDR= G*PSR/ASR
      END IF
      RETURN
      END
C
      SUBROUTINE RCS(T,PS,US)
      PARAMETER (IMAX=101,XB=0.,XE=1.,DX=(XE-XB)/(IMAX-1),X0=0.5)
      COMMON /C1/ PR,UR,DR,AR,PL,UL,DL,AL,G,TOL,IERROR
      COMMON /C2/ GP1,GM1,CL,CR,T1,T2,T3L,T3R,T4,T5
      COMMON /C3/ X(103),P(103),D(103),U(103),E(103)
C
      TM  = T3R*(US-UR)
      WR  = TM + SQRT(1+TM*TM)
      ASR = AR*SQRT( (GP1+GM1*PS/PR)/(GP1+GM1*PR/PS) )
      ASL = AL - T4*(US-UL)
      VSK = UR + AR*WR
      VRH = UL - AL
      VRT = US - ASL
      XSK = X0 + VSK*T
      XCD = X0 + US*T
      XRH = X0 + VRH*T
      XRT = X0 + VRT*T
      DO 10 I=1,IMAX
      X(I)= XB + DX*(I-1)
      IF( X(I) .GT. XSK ) THEN
          P(I) = PR
          U(I) = UR
          D(I) = DR
      ELSEIF( X(I) .LE. XSK .AND. X(I) .GT. XCD ) THEN
          P(I) = PS
          U(I) = US
          D(I) = G*PS/ASR/ASR
      ELSEIF( X(I) .LE. XCD .AND. X(I) .GT. XRT ) THEN
          P(I) = PS
          U(I) = US
          D(I) = G*PS/ASL/ASL
      ELSEIF( X(I) .LE. XRH ) THEN
          P(I) = PL
          U(I) = UL
          D(I) = DL
      ELSE
          UU   = (X(I)-X0)/T
          AS4  = 2./GP1 - GM1/GP1*(UU-UL)/AL
          P(I) = PL*AS4**T5
          U(I) = 2.*AL*(1.-AS4)/GM1 + UL
          D(I) = G*P(I)/(AS4*AL)/(AS4*AL)
      END IF
      E(I) = P(I)/GM1/D(I)
   10 CONTINUE
      DO 20 I=1,IMAX-1
      IF( X(I) .LT. XSK .and. X(I+1) .GT. XSK ) THEN
          X(I)= XSK
          X(I+1)= XSK
      END IF
      IF( X(I) .LT. XCD .and. X(I+1) .GT. XCD ) THEN
          X(I)= XCD
          X(I+1)= XCD
      END IF
   20 CONTINUE
      RETURN
      END
C
      SUBROUTINE SCR(T,PS,US)
      PARAMETER (IMAX=101,XB=0.,XE=1.,DX=(XE-XB)/(IMAX-1),X0=0.5)
      COMMON /C1/ PR,UR,DR,AR,PL,UL,DL,AL,G,TOL,IERROR
      COMMON /C2/ GP1,GM1,CL,CR,T1,T2,T3L,T3R,T4,T5
      COMMON /C3/ X(103),P(103),D(103),U(103),E(103)
C
      TM  = T3L*(US-UL)
      WL  = TM - SQRT(1+TM*TM)
      ASL = AL*SQRT( (GP1+GM1*PS/PL)/(GP1+GM1*PL/PS) )
      ASR = AR + T4*(US-UR)
      VSK = UL + AL*WL
      VRH = UR + AR
      VRT = US + ASR
      XSK = X0 + VSK*T
      XCD = X0 + US*T
      XRH = X0 + VRH*T
      XRT = X0 + VRT*T
      DO 10 I=1,IMAX
      X(I)= XB + DX*(I-1)
      IF( X(I) .LT. XSK ) THEN
          P(I) = PL
          U(I) = UL
          D(I) = DL
      ELSEIF( X(I) .GE. XSK .AND. X(I) .LT. XCD ) THEN
          P(I) = PS
          U(I) = US
          D(I) = G*PS/ASL/ASL
      ELSEIF( X(I) .GE. XCD .AND. X(I) .LT. XRT ) THEN
          P(I) = PS
          U(I) = US
          D(I) = G*PS/ASR/ASR
      ELSEIF( X(I) .GE. XRH ) THEN
          P(I) = PR
          U(I) = UR
          D(I) = DR
      ELSE
          UU   = (X(I)-X0)/T
          AS1  = 2./GP1 + GM1/GP1*(UU-UR)/AR
          P(I) = PR*AS1**T5
          U(I) = 2.*AR*(AS1-1.)/GM1 + UR
          D(I) = G*P(I)/(AS1*AR)/(AS1*AR)
      END IF
      E(I) = P(I)/GM1/D(I)
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE SCS(T,PS,US)
      PARAMETER (IMAX=101,XB=0.,XE=1.,DX=(XE-XB)/(IMAX-1),X0=0.5)
      COMMON /C1/ PR,UR,DR,AR,PL,UL,DL,AL,G,TOL,IERROR
      COMMON /C2/ GP1,GM1,CL,CR,T1,T2,T3L,T3R,T4,T5
      COMMON /C3/ X(103),P(103),D(103),U(103),E(103)
C
      TM  = T3R*(US-UR)
      WR  = TM + SQRT(1+TM*TM)
      ASR = AR*SQRT( (GP1+GM1*PS/PR)/(GP1+GM1*PR/PS) )
      TM  = T3L*(US-UL)
      WL  = TM - SQRT(1+TM*TM)
      ASL = AL*SQRT( (GP1+GM1*PS/PL)/(GP1+GM1*PL/PS) )
      VSKR= UR + AR*WR
      VSKL= UL + AL*WL
      XSKR= X0 + VSKR*T
      XSKL= X0 + VSKL*T
      XCD = X0 + US*T
      DO 10 I=1,IMAX
      X(I)= XB + DX*(I-1)
      IF( X(I) .GT. XSKR ) THEN
          P(I) = PR
          U(I) = UR
          D(I) = DR
      ELSEIF( X(I) .LE. XSKR .AND. X(I) .GT. XCD ) THEN
          P(I) = PS
          U(I) = US
          D(I) = G*PS/ASR/ASR
      ELSEIF( X(I) .LE. XCD .AND. X(I) .GT. XSKL ) THEN
          P(I) = PS
          U(I) = US
          D(I) = G*PS/ASL/ASL
      ELSE
          P(I) = PL
          U(I) = UL
          D(I) = DL
      END IF
      E(I) = P(I)/GM1/D(I)
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE RCR(T,PS,US)
      PARAMETER (IMAX=101,XB=0.,XE=1.,DX=(XE-XB)/(IMAX-1),X0=0.5)
      COMMON /C1/ PR,UR,DR,AR,PL,UL,DL,AL,G,TOL,IERROR
      COMMON /C2/ GP1,GM1,CL,CR,T1,T2,T3L,T3R,T4,T5
      COMMON /C3/ X(103),P(103),D(103),U(103),E(103)
C
      ASR = AR + T4*(US-UR)
      ASL = AL - T4*(US-UL)
      VRHR= UR + AR
      VRTR= US + ASR
      VRHL= UL - AL
      VRTL= US - ASL
      XCD = X0 + US*T
      XRHR= X0 + VRHR*T
      XRTR= X0 + VRTR*T
      XRHL= X0 + VRHL*T
      XRTL= X0 + VRTL*T
      DO 10 I=1,IMAX
      X(I)= XB + DX*(I-1)
      IF( X(I) .GT. XRHR ) THEN
          P(I) = PR
          U(I) = UR
          D(I) = DR
      ELSEIF( X(I) .LE. XRTR .AND. X(I) .GT. XCD ) THEN
          P(I) = PS
          U(I) = US
          D(I) = G*PS/ASR/ASR
      ELSEIF( X(I) .LE. XCD .AND. X(I) .GT. XRTL ) THEN
          P(I) = PS
          U(I) = US
          D(I) = G*PS/ASL/ASL
      ELSEIF( X(I) .LE. XRHL ) THEN
          P(I) = PL
          U(I) = UL
          D(I) = DL
      ELSEIF( X(I) .LE. XRTL .AND. X(I) .GT. XRHL ) THEN
          UU   = (X(I)-X0)/T
          AS4  = 2./GP1 - GM1/GP1*(UU-UL)/AL
          P(I) = PL*AS4**T5
          U(I) = 2.*AL*(1.-AS4)/GM1 + UL
          D(I) = G*P(I)/(AS4*AL)/(AS4*AL)
      ELSE
          UU   = (X(I)-X0)/T
          AS1  = 2./GP1 + GM1/GP1*(UU-UR)/AR
          P(I) = PR*AS1**T5
          U(I) = 2.*AR*(AS1-1.)/GM1 + UR
          D(I) = G*P(I)/(AS1*AR)/(AS1*AR)
      END IF
      E(I) = P(I)/GM1/D(I)
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE GRAPH(X,RQ,UQ,PQ,EQ,MAXI)
      DIMENSION X(1),RQ(1),UQ(1),PQ(1),EQ(1)
      COMMON /C3/ XEXA(103),PEXA(103),REXA(103),UEXA(103),EEXA(103)
      DATA PAGEY,PAGEX,AREAY,AREAX,OY,OX / 25.5,19.,6.,5.,8.5,5.0 /
    1 NEXA = 101
C
      CALL PLOTST
      CALL WINDOW( 0.,PAGEX,0.,PAGEY )
      CALL PLOT( OX,OY,-3 )
      CALL NEWPEN(3)
      YMIN = 1.E20
      YMAX =-1.E20
      DO 10 I = 1,NEXA
      IF( UEXA(I) .LE. YMIN ) YMIN= UEXA(I)
      IF( UEXA(I) .GE. YMAX ) YMAX= UEXA(I)
   10 CONTINUE
      YMIN = YMIN - 0.15*(YMAX-YMIN)
      YMAX = YMAX + 0.15*(YMAX-YMIN)
      YINC = (YMAX-YMIN)/AREAY
      XEXA(NEXA+1) = 0.
      XEXA(NEXA+2) = 0.2
      X(MAXI+1)    = 0.
      X(MAXI+2)    = 0.2
      UEXA(NEXA+1) = YMIN
      UEXA(NEXA+2) = YINC
      UQ(MAXI+1)   = YMIN
      UQ(MAXI+2)   = YINC
      CALL SETCLR(15)
      CALL AXIS(0.,0.,'X',-1,AREAX, 0.,0.0,0.2)
      CALL AXIS(0.,0.,'VELOCITY', 8,AREAY,90.,YMIN,YINC)
      CALL SETCLR( 2)
      CALL LINEM(XEXA,UEXA,NEXA,1,0,0,0.6)
      CALL SETCLR( 5)
      CALL LINEM(X,UQ,MAXI,1,-1,1,0.6)
C
      CALL PLOT(AREAX+2.5,0.,-3)
      YMIN = 1.E20
      YMAX =-1.E20
      DO 20 I = 1,NEXA
      IF( EEXA(I) .LE. YMIN ) YMIN= EEXA(I)
      IF( EEXA(I) .GE. YMAX ) YMAX= EEXA(I)
   20 CONTINUE
      YMIN = YMIN - 0.15*(YMAX-YMIN)
      YMAX = YMAX + 0.15*(YMAX-YMIN)
      YINC = (YMAX-YMIN)/AREAY
      EEXA(NEXA+1) = YMIN
      EEXA(NEXA+2) = YINC
      EQ(MAXI+1)   = YMIN
      EQ(MAXI+2)   = YINC
      CALL SETCLR(15)
      CALL AXIS(0.,0.,'X',-1,AREAX, 0.,0.0,0.2)
      CALL AXIS(0.,0.,'INTERNAL ENERGY',15,AREAY,90.,YMIN,YINC)
      CALL SETCLR( 2)
      CALL LINEM(XEXA,EEXA,NEXA,1,0,0,0.6)
      CALL SETCLR( 5)
      CALL LINEM(X,EQ,MAXI,1,-1,1,0.6)
C
      CALL PLOT(-AREAX-2.5,AREAY+2.5,-3)
      YMIN= 1.E20
      YMAX=-1.E20
      DO 30 I=1,NEXA
      IF( REXA(I) .LE. YMIN ) YMIN= REXA(I)
      IF( REXA(I) .GE. YMAX ) YMAX= REXA(I)
   30 CONTINUE
      YMIN = YMIN - 0.15*(YMAX-YMIN)
      YMAX = YMAX + 0.15*(YMAX-YMIN)
      YINC = (YMAX-YMIN)/AREAY
      REXA(NEXA+1) = YMIN
      REXA(NEXA+2) = YINC
      RQ(MAXI+1)   = YMIN
      RQ(MAXI+2)   = YINC
      CALL SETCLR(15)
      CALL AXIS(0.,0.,'X',-1,AREAX, 0.,0.0,0.2)
      CALL AXIS(0.,0.,'DENSITY',7,AREAY,90.,YMIN,YINC)
      CALL SETCLR( 2)
      CALL LINEM(XEXA,REXA,NEXA,1,0,0,0.6)
      CALL SETCLR( 5)
      CALL LINEM(X,RQ,MAXI,1,-1,1,0.6)
C
      CALL PLOT(AREAX+2.5,0.,-3)
      YMIN = 1.E20
      YMAX =-1.E20
      DO 40 I = 1,NEXA
      IF( PEXA(I) .LE. YMIN ) YMIN= PEXA(I)
      IF( PEXA(I) .GE. YMAX ) YMAX= PEXA(I)
   40 CONTINUE
      YMIN = YMIN - 0.15*(YMAX-YMIN)
      YMAX = YMAX + 0.15*(YMAX-YMIN)
      YINC = (YMAX-YMIN)/AREAY
      PEXA(NEXA+1) = YMIN
      PEXA(NEXA+2) = YINC
      PQ(MAXI+1)   = YMIN
      PQ(MAXI+2)   = YINC
      CALL SETCLR(15)
      CALL AXIS(0.,0.,'X',-1,AREAX, 0.,0.0,0.2)
      CALL AXIS(0.,0.,'PRESSURE',8,AREAY,90.,YMIN,YINC)
      CALL SETCLR( 2)
      CALL LINEM(XEXA,PEXA,NEXA,1,0,0,0.6)
      CALL SETCLR( 5)
      CALL LINEM(X,PQ,MAXI,1,-1,1,0.6)
      CALL PLOTND
      GOTO 1
      RETURN
      END
C
      SUBROUTINE GFUNCT(IMOD,I1,I2,A,B,IMTH)
      PARAMETER( MAXI=407,NV=3 )
      DIMENSION A(MAXI,NV),B(MAXI,NV)
      SM   = 1.E-6
      I1P1 = I1 + 1
      I2M1 = I2 - 1
      I1P2 = I1 + 2
      I2M2 = I2 - 2
C
      IF( IMTH .EQ. 2 ) GOTO 200               ! ENO2
C
      GOTO (10,20,30,40,50,60,70,80,90),IMOD
C
  10  CONTINUE
      DO 12 N= 1,NV
      DO 12 I= I1P1,I2M1
      Z1     = A( I ,N)
      Z2     = A(I-1,N)
      B(I,N) = 0.5*(SGN(Z1)+SGN(Z2))*MIN(ABS(Z1),ABS(Z2))
  12  CONTINUE
      GOTO 100
C
  20  CONTINUE
      DO 22 N= 1,NV
      DO 22 I= I1P1,I2M1
      Z1     = A( I ,N)
      Z2     = A(I-1,N)
      B(I,N) = ( Z1*(Z2*Z2+SM)+Z2*(Z1*Z1+SM) )/(Z1*Z1+Z2*Z2+2.*SM)
  22  CONTINUE
      GOTO 100
C
  30  CONTINUE
      DO 32 N= 1,NV
      DO 32 I= I1P1,I2M1
      Z1     = A( I ,N)*2.
      Z2     = A(I-1,N)*2.
      Z3     = (A(I,N)+A(I-1,N))/2.
      Z4     = (SGN(Z1)+SGN(Z2)+SGN(Z3))/3.
      B(I,N) = 0.
      IF(ABS(Z4) .EQ. 1.) B(I,N)= Z4*MIN(ABS(Z1),ABS(Z2),ABS(Z3))
  32  CONTINUE
      GOTO 100
C
  40  CONTINUE
      DO 42 N= 1,NV
      DO 42 I= I1P1,I2M1
      Z1     = A( I ,N)
      Z2     = A(I-1,N)
      B(I,N) = 0.
      IF( Z1+Z2 .NE. 0. ) B(I,N) = (Z1*Z2+ABS(Z1*Z2))/(Z1+Z2)
  42  CONTINUE
      GOTO 100
C
  50  CONTINUE
      DO 52 N = 1,NV
      DO 52 I = I1P1,I2M1
      Z1      = A( I ,N)
      Z2      = A(I-1,N)
      SS      = SGN(Z1)
      B(I,N)  = SS*MAX(0.,MIN(2.*SS*Z1,SS*Z2),
     &                    MIN(SS*Z1,2.*SS*Z2) )
  52  CONTINUE
      GOTO 100
C
  60  CONTINUE
      DO 62 N = 1,NV
      DO 62 I = I1P1,I2M1
      Z1      = A( I ,N)
      Z2      = A(I-1,N)
      IF( N .LE. NV-2 ) THEN
          SS    = SGN(Z1)
          B(I,N)= SS*MAX(0.,MIN(2.*SS*Z1,SS*Z2),
     &                      MIN(SS*Z1,2.*SS*Z2) )
      ELSE
          B(I,N)= 0.5*(SGN(Z1)+SGN(Z2))*MIN(ABS(Z1),ABS(Z2))
      END IF
  62  CONTINUE
      GOTO 100
C
  70  CONTINUE
      DO 72 N = 1,NV
      DO 72 I = I1P1,I2M1
      Z1      = A( I ,N)
      Z2      = A(I-1,N)
      IF( N .LE. NV-2 ) THEN
          SS    = SGN(Z1)
          B(I,N)= SS*MAX(0.,MIN(2.*SS*Z1,SS*Z2),
     &                      MIN(SS*Z1,2.*SS*Z2) )
      ELSE
          B(I,N)= ( Z1*(Z2*Z2+SM)+Z2*(Z1*Z1+SM) )/(Z1*Z1+Z2*Z2+2.*SM)
      END IF
  72  CONTINUE
      GOTO 100
C
  80  CONTINUE
      DO 82 N = 1,NV
      DO 82 I = I1P1,I2M1
      Z1      = A( I ,N)
      Z2      = A(I-1,N)
      IF( N .LE. NV-2 ) THEN
          SS    = SGN(Z1)
          B(I,N)= SS*MAX(0.,MIN(2.*SS*Z1,SS*Z2),
     &                      MIN(SS*Z1,2.*SS*Z2) )
      ELSE
          Z3    = (A(I,N)+A(I-1,N))/2.
          Z4    = (SGN(Z1)+SGN(Z2)+SGN(Z3))/3.
          B(I,N)= 0.
          IF(ABS(Z4) .EQ. 1.) B(I,N)= Z4*MIN(ABS(Z1),ABS(Z2),ABS(Z3))
      END IF
  82  CONTINUE
      GOTO 100
C
  90  CONTINUE
      DO 92 N = 1,NV
      DO 92 I = I1P1,I2M1
      Z1      = A( I ,N)
      Z2      = A(I-1,N)
      IF( N .LE. NV-2 ) THEN
          SS    = SGN(Z1)
          B(I,N)= SS*MAX(0.,MIN(2.*SS*Z1,SS*Z2),
     &                      MIN(SS*Z1,2.*SS*Z2) )
      ELSE
          B(I,N)= 0.
          IF( Z1+Z2 .NE. 0. ) B(I,N) = (Z1*Z2+ABS(Z1*Z2))/(Z1+Z2)
      END IF
  92  CONTINUE
      GOTO 100
C
 100  CONTINUE
C
C     B.C. for B  (zero order extrapolation)
C
      DO 105 N= 1,NV
      B(I1,N) = B(I1P1,N)
      B(I2,N) = B(I2M1,N)
 105  CONTINUE
      RETURN
C
 200  CONTINUE
      DO 210 N = 1,NV
      DO 205 I = I1P2,I2M2
      P1       = A(I+1,N) - A( I ,N)
      P2       = A( I ,N) - A(I-1,N)
      P4       = A(I-1,N) - A(I-2,N)
      TMP1     = P2
      TMP2     = P2
      IF( ABS(P1) .LT. ABS(P2) ) TMP1= P1
      IF( ABS(P4) .LT. ABS(P2) ) TMP2= P4
C     TMP1= 0.5*(SGN(P1)+SGN(P2))*MIN(ABS(P1),ABS(P2))
C     TMP2= 0.5*(SGN(P4)+SGN(P2))*MIN(ABS(P4),ABS(P2))
      Z1  = A( I ,N) - 0.5*TMP1
      Z2  = A(I-1,N) + 0.5*TMP2
      IF( IMOD .EQ. 1 ) THEN
          B(I,N) = 0.5*(SGN(Z1)+SGN(Z2))*MIN(ABS(Z1),ABS(Z2))
      ELSEIF( IMOD .EQ. 2 ) THEN
          B(I,N) = ( Z1*(Z2*Z2+SM)+Z2*(Z1*Z1+SM) )/(Z1*Z1+Z2*Z2+2.*SM)
      ELSEIF( IMOD .EQ. 3 ) THEN
          Z3     = (A(I,N)+A(I-1,N))/2.
          Z4     = (SGN(Z1)+SGN(Z2)+SGN(Z3))/3.
          B(I,N) = 0.
          IF(ABS(Z4) .EQ. 1.) B(I,N)= Z4*MIN(ABS(Z1),ABS(Z2),ABS(Z3))
      ELSEIF( IMOD .EQ. 4 ) THEN
          B(I,N) = 0.
          IF( Z1+Z2 .NE. 0. ) B(I,N)= (Z1*Z2+ABS(Z1*Z2))/(Z1+Z2)
      ELSEIF( IMOD .EQ. 5 ) THEN
          SS      = SGN(Z1)
          B(I,N)  = SS*MAX(0.,MIN(2.*SS*Z1,SS*Z2),
     &                        MIN(SS*Z1,2.*SS*Z2) )
      ELSE
          IF( N .LE. NV-2 ) THEN
              SS    = SGN(Z1)
              B(I,N)= SS*MAX(0.,MIN(2.*SS*Z1,SS*Z2),
     &                          MIN(SS*Z1,2.*SS*Z2) )
          ELSE
              IF( IMOD .EQ. 6 ) THEN
                  B(I,N)=0.5*(SGN(Z1)+SGN(Z2))*MIN(ABS(Z1),ABS(Z2))
              ELSEIF( IMOD .EQ. 7 ) THEN
                  Z3 = Z1*Z1+SM
                  Z4 = Z2*Z2+SM
                  B(I,N)=( Z1*Z4 + Z2*Z3 )/( Z3+Z4 )
              ELSEIF( IMOD .EQ. 8 ) THEN
                  Z3     = (A(I,N)+A(I-1,N))/2.
                  Z4     = (SGN(Z1)+SGN(Z2)+SGN(Z3))/3.
                  B(I,N) = 0.
                  IF( ABS(Z4) .EQ. 1. )
     &                B(I,N)= Z4*MIN(ABS(Z1),ABS(Z2),ABS(Z3))
              ELSE
                  B(I,N)= 0.
                  IF( Z1+Z2 .NE. 0. ) B(I,N)=(Z1*Z2+ABS(Z1*Z2))/(Z1+Z2)
              END IF
          END IF
      END IF
 205  CONTINUE
      Z1       = A(I2M2,N)
      Z2       = A(I2M1,N)
      B(I2M1,N)= 0.
      IF( Z1+Z2 .NE. 0. ) B(I2M1,N)= (Z1*Z2+ABS(Z1*Z2))/(Z1+Z2)
      Z1       = A(I1P1,N)
      Z2       = A(I1  ,N)
      B(I1P1,N)= 0.
      IF( Z1+Z2 .NE. 0. ) B(I1P1,N)= (Z1*Z2+ABS(Z1*Z2))/(Z1+Z2)
 210  CONTINUE
C
C     B.C. for B  (zero order extrapolation)
C
      DO 215 N = 1,NV
      B(I1,N)  = B(I1P1,N)
      B(I2,N)  = B(I2M1,N)
 215  CONTINUE
      RETURN
      END
C
      SUBROUTINE QFUNCT(IMOD,I1,I2,A,B)
      PARAMETER( MAXI=407,NV=3 )
      DIMENSION A(MAXI,NV),B(MAXI,NV)
      I1P = I1 + 1
      I2M = I2 - 1
      GOTO (10,20,30,40,50,10,20,30,40),IMOD
  10  CONTINUE
      DO 12 N= 1,NV
      DO 12 I= I1P,I2M
      C1     = A(I-1,N)
      C2     = A( I ,N)
      C3     = A(I+1,N)
      C4     = 0.5*( SGN(C1)+SGN(C2) )*MIN( ABS(C1),ABS(C2) )
      C5     = 0.5*( SGN(C3)+SGN(C2) )*MIN( ABS(C3),ABS(C2) )
      B(I,N) = C4 + C5 - C2
  12  CONTINUE
      GOTO 60
  20  CONTINUE
      DO 22 N= 1,NV
      DO 22 I= I1P,I2M
      C1     = A(I-1,N)
      C2     = A( I ,N)
      C3     = A(I+1,N)
      B(I,N) = 0.
      C4     = ( SGN(C1)+SGN(C2)+SGN(C3) )/3.
      IF(ABS(C4) .EQ. 1.) B(I,N) = C4*MIN(ABS(C1),ABS(C2),ABS(C3))
  22  CONTINUE
      GOTO 60
  30  CONTINUE
      DO 32 N= 1,NV
      DO 32 I= I1P,I2M
      C1     = A(I-1,N)*2.
      C2     = A( I ,N)*2.
      C3     = A(I+1,N)*2.
      C4     = (A(I-1,N)+A(I+1,N))/2.
      B(I,N) = 0.
      C5     = (SGN(C1)+SGN(C2)+SGN(C3)+SGN(C4))/4.
      IF(ABS(C5) .EQ. 1.) B(I,N)=
     &                    C5*MIN(ABS(C1),ABS(C2),ABS(C3),ABS(C4))
  32  CONTINUE
      GOTO 60
  40  CONTINUE
      DO 42 N= 1,NV
      DO 42 I= I1P,I2M
      C1     = A(I-1,N)
      C2     = A( I ,N)
      C3     = A(I+1,N)
      C4     = 0.
      C5     = 0.
      IF( C1+C2 .NE. 0. ) C4 = (C1*C2+ABS(C1*C2))/(C1+C2)
      IF( C3+C2 .NE. 0. ) C5 = (C3*C2+ABS(C3*C2))/(C3+C2)
      B(I,N) = C4 + C5 - C2
  42  CONTINUE
      GOTO 60
  50  CONTINUE
      DO 52 N= 1,NV-2
      DO 52 I= I1P,I2M
      C1     = A(I-1,N)
      C2     = A( I ,N)
      C3     = A(I+1,N)
      S      = SGN(C2)
      C2A    = ABS(C2)
      C4     = S*MAX( 0.,MIN(2.*C2A,C1*S),MIN(C2A,2.*C1*S) )
      C5     = S*MAX( 0.,MIN(2.*C2A,C3*S),MIN(C2A,2.*C3*S) )
      B(I,N) = C4 + C5 - C2
  52  CONTINUE
      DO 54 N= NV-1,NV
      DO 54 I= I1P,I2M
      C1     = A(I-1,N)
      C2     = A( I ,N)
      C3     = A(I+1,N)
      C4     = 0.
      C5     = 0.
      IF( C1+C2 .NE. 0. ) C4 = (C1*C2+ABS(C1*C2))/(C1+C2)
      IF( C3+C2 .NE. 0. ) C5 = (C3*C2+ABS(C3*C2))/(C3+C2)
      B(I,N) = C4 + C5 - C2
  54  CONTINUE
  60  CONTINUE
      DO 70 N= 1,NV
      B(I1,N)= B(I1P,N)
      B(I2,N)= B(I2M,N)
  70  CONTINUE
      RETURN
      END
C
      SUBROUTINE RUNGE(IB,IE,CK,OMIGA,Q,IMTH,IMOD,IOPT)
      PARAMETER ( MAXI=407,NV=3 )
      COMMON /A1/ RQ(MAXI),UQ(MAXI),PQ(MAXI),EQ(MAXI),CQ(MAXI),HQ(MAXI)
      COMMON /D1/ GAMMA,GM1,EPSI,LEVEL,DTX,IPROB
      DIMENSION Q(MAXI,NV),S(MAXI,NV,4),S0(MAXI,NV),A(4)
      DATA A/1., 1., 1., 0.5/
C
*****************************************************************
*     Runge-Kutta in time
*     Level = 1 ---> 1st order Runge-Kutta (TVD)
*     Level = 2 ---> 2nd order Runge-Kutta (TVD)
*     Level = 3 ---> 3rd order Runge-Kutta (TVD)
*     Level = 4 ---> 4th order Runge-Kutta (non-TVD)
*****************************************************************
C
      IF( IMTH .LE. 30 ) THEN
          IF( IOPT .EQ. 1 ) THEN
              CALL MUSCL1(IB,IE,CK,OMIGA,Q,S(1,1,1),IMTH,IMOD)
          ELSEIF( IOPT .EQ. 2 ) THEN
              CALL MUSCL2(IB,IE,CK,OMIGA,Q,S(1,1,1),IMTH,IMOD)
          ELSE
              CALL MUSCL3(IB,IE,CK,OMIGA,Q,S(1,1,1),IMTH,IMOD)
          END IF
      ELSEIF( IMTH .EQ. 44 ) THEN
          CALL TFOUR(0,IB,IE,Q,S(1,1,1),IMTH)
      ELSE
          CALL CENTR(0,IB,IE,Q,S(1,1,1),IMTH)
      END IF
      IF( LEVEL .GE. 2 .AND. IMTH .EQ. 43 ) THEN
          CALL DAMPER(IB,IE,Q,S0,IMTH,IMOD)
      END IF
C
C     UPDATE Q WITH CORRECTIONS (STAGE 1)
C
      DO 10 I = IB+1,IE-1
      Q(I,1)  = Q(I,1) + A(LEVEL)*DTX*S(I,1,1)
      Q(I,2)  = Q(I,2) + A(LEVEL)*DTX*S(I,2,1)
      Q(I,3)  = Q(I,3) + A(LEVEL)*DTX*S(I,3,1)
   10 CONTINUE
C
      CALL PRIMER( IB,IE,Q )
      IF( LEVEL .EQ. 1 ) RETURN
C
      IF( IMTH .LE. 30 ) THEN
          IF( IOPT .EQ. 1 ) THEN
              CALL MUSCL1(IB,IE,CK,OMIGA,Q,S(1,1,2),IMTH,IMOD)
          ELSEIF( IOPT .EQ. 2 ) THEN
              CALL MUSCL2(IB,IE,CK,OMIGA,Q,S(1,1,2),IMTH,IMOD)
          ELSE
              CALL MUSCL3(IB,IE,CK,OMIGA,Q,S(1,1,2),IMTH,IMOD)
          END IF
      ELSEIF( IMTH .EQ. 44 ) THEN
          CALL TFOUR(1,IB,IE,Q,S(1,1,2),IMTH)
      ELSE
          CALL CENTR(1,IB,IE,Q,S(1,1,2),IMTH)
      END IF
C
C     UPDATE Q WITH CORRECTIONS (STAGE 2)
C
      IF( LEVEL .EQ. 2 .OR. LEVEL .EQ. 4 ) THEN
          DO 20 I = IB+1,IE-1
          Q(I,1)  = Q(I,1) + 0.5*DTX*( -S(I,1,1) +S(I,1,2) )
          Q(I,2)  = Q(I,2) + 0.5*DTX*( -S(I,2,1) +S(I,2,2) )
          Q(I,3)  = Q(I,3) + 0.5*DTX*( -S(I,3,1) +S(I,3,2) )
   20     CONTINUE
          IF( LEVEL .EQ. 2 .AND. IMTH .EQ. 43 ) THEN
              DO 25 I = IB+1,IE-1
              Q(I,1)  = Q(I,1) + S0(I,1)
              Q(I,2)  = Q(I,2) + S0(I,2)
              Q(I,3)  = Q(I,3) + S0(I,3)
   25         CONTINUE
          END IF
      END IF
      IF( LEVEL .EQ. 3 ) THEN
          DO 30 I = IB+1,IE-1
          Q(I,1)  = Q(I,1) + 0.25*DTX*( -3.*S(I,1,1) +S(I,1,2) )
          Q(I,2)  = Q(I,2) + 0.25*DTX*( -3.*S(I,2,1) +S(I,2,2) )
          Q(I,3)  = Q(I,3) + 0.25*DTX*( -3.*S(I,3,1) +S(I,3,2) )
   30     CONTINUE
      END IF
C
      CALL PRIMER( IB,IE,Q )
      IF( LEVEL .EQ. 2 ) RETURN
C
      IF( IMTH .LE. 30 ) THEN
          IF( IOPT .EQ. 1 ) THEN
              CALL MUSCL1(IB,IE,CK,OMIGA,Q,S(1,1,3),IMTH,IMOD)
          ELSEIF( IOPT .EQ. 2 ) THEN
              CALL MUSCL2(IB,IE,CK,OMIGA,Q,S(1,1,3),IMTH,IMOD)
          ELSE
              CALL MUSCL3(IB,IE,CK,OMIGA,Q,S(1,1,3),IMTH,IMOD)
          END IF
      END IF
C
C     UPDATE Q WITH CORRECTIONS (STAGE 3)
C
      IF( LEVEL .EQ. 3 ) THEN
          DO 40 I = IB+1,IE-1
          Q(I,1)  = Q(I,1) +DTX/12.*( -S(I,1,1) -S(I,1,2) +8.*S(I,1,3) )
          Q(I,2)  = Q(I,2) +DTX/12.*( -S(I,2,1) -S(I,2,2) +8.*S(I,2,3) )
          Q(I,3)  = Q(I,3) +DTX/12.*( -S(I,3,1) -S(I,3,2) +8.*S(I,3,3) )
   40     CONTINUE
      END IF
      IF( LEVEL .EQ. 4 ) THEN
          DO 50 I = IB+1,IE-1
          Q(I,1)  = Q(I,1) +DTX*( -0.5*S(I,1,2) +S(I,1,3) )
          Q(I,2)  = Q(I,2) +DTX*( -0.5*S(I,2,2) +S(I,2,3) )
          Q(I,3)  = Q(I,3) +DTX*( -0.5*S(I,3,2) +S(I,3,3) )
   50     CONTINUE
      END IF
C
      CALL PRIMER( IB,IE,Q )
      IF( LEVEL .EQ. 3 ) RETURN
C
      IF( IMTH .LE. 30 ) THEN
          IF( IOPT .EQ. 1 ) THEN
              CALL MUSCL1(IB,IE,CK,OMIGA,Q,S(1,1,4),IMTH,IMOD)
          ELSEIF( IOPT .EQ. 2 ) THEN
              CALL MUSCL2(IB,IE,CK,OMIGA,Q,S(1,1,4),IMTH,IMOD)
          ELSE
              CALL MUSCL3(IB,IE,CK,OMIGA,Q,S(1,1,4),IMTH,IMOD)
          END IF
      END IF
C
C     UPDATE Q WITH CORRECTIONS (STAGE 4)
C
      IF( LEVEL .EQ. 4 ) THEN
          DO 60 I = IB+1,IE-1
          Q(I,1)  = Q(I,1) +DTX/6.*( S(I,1,1) +2.*S(I,1,2) +
     &                               S(I,1,4) -4.*S(I,1,3) )
          Q(I,2)  = Q(I,2) +DTX/6.*( S(I,2,1) +2.*S(I,2,2) +
     &                               S(I,2,4) -4.*S(I,2,3) )
          Q(I,3)  = Q(I,3) +DTX/6.*( S(I,3,1) +2.*S(I,3,2) +
     &                               S(I,3,4) -4.*S(I,3,3) )
   60     CONTINUE
      END IF
      RETURN
      END
C
      SUBROUTINE MUSCL1( JB,JE,CK,OMIGA,Q,S,IMTH,IMOD )
      PARAMETER ( MAXI=407,NV=3 )
      COMMON /A1/ RQ(MAXI),UQ(MAXI),PQ(MAXI),EQ(MAXI),CQ(MAXI),HQ(MAXI)
      DIMENSION FN(MAXI,NV),BET(MAXI,NV),EIG(MAXI,NV),TT(MAXI,NV,NV)
      DIMENSION UA(MAXI),CA(MAXI),Q(MAXI,NV),S(MAXI,NV),BB(NV),CC(NV)
      DIMENSION QR(MAXI,NV),QL(MAXI,NV),FR(MAXI,NV),FL(MAXI,NV)
      DIMENSION DF(MAXI,NV),PL(MAXI),PR(MAXI),TK(MAXI,NV,NV)
      DIMENSION GB(MAXI,NV),GC(MAXI,NV),DQ(MAXI,NV,-2:1),TKI(MAXI,NV,NV)
      COMMON /D1/ GAMMA,GM1,EPSI,LEVEL,DTX,IPROB
C
C     LIMITER ON CHARACTERISTIC VARIABLES
C
      CKM   = 1. - CK
      CKP   = 1. + CK
      JUP   = JE - 1
      JLOW  = JB + 1
C
      DO 10 N = 1,NV
      QL(JB,N)= Q(JB,N)
      QR(JB,N)= Q(JB,N)
      QL(JE,N)= Q(JE,N)
      QR(JE,N)= Q(JE,N)
      QL(JLOW,N) = Q(JLOW,N)
      QR(JLOW,N) = Q(JLOW+1,N)
      QL(JUP ,N) = Q(JUP,N)
      QR(JUP ,N) = Q(JE,N)
   10 CONTINUE
C
C     CALCULATE EIGENVALUE
C
      DO 15 J = JB,JE
      EIG(J,1)= UQ(J)
      EIG(J,2)= UQ(J) + CQ(J)
      EIG(J,3)= UQ(J) - CQ(J)
   15 CONTINUE
C
C     CALCULATE AFAR
C
      CALL RLVECT( MAXI,JE,UQ,CQ,HQ,RQ,TK,TKI )
C
      DO 20 J = JB,JUP
      Q1         = Q(J+1,1) - Q(J,1)
      Q2         = Q(J+1,2) - Q(J,2)
      Q3         = Q(J+1,3) - Q(J,3)
      DQ(J,1,0)  = TKI(J,1,1)*Q1 + TKI(J,1,2)*Q2 + TKI(J,1,3)*Q3
      DQ(J,2,0)  = TKI(J,2,1)*Q1 + TKI(J,2,2)*Q2 + TKI(J,2,3)*Q3
      DQ(J,3,0)  = TKI(J,3,1)*Q1 + TKI(J,3,2)*Q2 + TKI(J,3,3)*Q3
   20 CONTINUE
      DO 21 J = JB,JUP-1
      Q1         = Q(J+2,1) - Q(J+1,1)
      Q2         = Q(J+2,2) - Q(J+1,2)
      Q3         = Q(J+2,3) - Q(J+1,3)
      DQ(J,1,1)  = TKI(J,1,1)*Q1 + TKI(J,1,2)*Q2 + TKI(J,1,3)*Q3
      DQ(J,2,1)  = TKI(J,2,1)*Q1 + TKI(J,2,2)*Q2 + TKI(J,2,3)*Q3
      DQ(J,3,1)  = TKI(J,3,1)*Q1 + TKI(J,3,2)*Q2 + TKI(J,3,3)*Q3
   21 CONTINUE
      DO 22 J = JLOW,JE
      Q1         = Q(J,1) - Q(J-1,1)
      Q2         = Q(J,2) - Q(J-1,2)
      Q3         = Q(J,3) - Q(J-1,3)
      DQ(J,1,-1) = TKI(J,1,1)*Q1 + TKI(J,1,2)*Q2 + TKI(J,1,3)*Q3
      DQ(J,2,-1) = TKI(J,2,1)*Q1 + TKI(J,2,2)*Q2 + TKI(J,2,3)*Q3
      DQ(J,3,-1) = TKI(J,3,1)*Q1 + TKI(J,3,2)*Q2 + TKI(J,3,3)*Q3
   22 CONTINUE
      DO 23 J = JLOW+1,JE
      Q1         = Q(J-1,1) - Q(J-2,1)
      Q2         = Q(J-1,2) - Q(J-2,2)
      Q3         = Q(J-1,3) - Q(J-2,3)
      DQ(J,1,-2) = TKI(J,1,1)*Q1 + TKI(J,1,2)*Q2 + TKI(J,1,3)*Q3
      DQ(J,2,-2) = TKI(J,2,1)*Q1 + TKI(J,2,2)*Q2 + TKI(J,2,3)*Q3
      DQ(J,3,-2) = TKI(J,3,1)*Q1 + TKI(J,3,2)*Q2 + TKI(J,3,3)*Q3
   23 CONTINUE
C
      IF( IMTH .EQ. 11 .OR. IMTH .EQ. 21 ) THEN  ! CHAKRAVARTHY-OSHER TVD
          DO 31 J = JLOW,JUP                     ! OK 09/21/1997
          DO 30 N = 1,NV
          Z1      = SGN( DQ(J,N,-1) )/2.
          Z2      = SGN( DQ(J,N, 0) )/2.
          Z3      = (Z1+Z2)*MIN( ABS(DQ(J,N,-1)),ABS(OMIGA*DQ(J,N,0)) )
          Z4      = (Z1+Z2)*MIN( ABS(OMIGA*DQ(J,N,-1)),ABS(DQ(J,N,0)) )
          BB(N)   = 0.25*( CKM*Z3 + CKP*Z4 )
          CC(N)   = 0.25*( CKP*Z3 + CKM*Z4 )
   30     CONTINUE
          DO 31 N = 1,NV
          QL(J  ,N) = Q(J,N) + TK(J,N,1)*BB(1) + TK(J,N,2)*BB(2)
     &                       + TK(J,N,3)*BB(3)
          QR(J-1,N) = Q(J,N) - TK(J,N,1)*CC(1) - TK(J,N,2)*CC(2)
     &                       - TK(J,N,3)*CC(3)
   31     CONTINUE
      END IF
C
      IF( IMTH .EQ. 12 .OR. IMTH .EQ. 13 .OR.
     &    IMTH .EQ. 22 .OR. IMTH .EQ. 23 ) THEN
C
C     CALCULATE MINMOD FOR GBAR
C
          DO 36 N= 1,NV
          DO 35 J= JLOW+1,JUP-1
          P1   = DQ(J,N, 1) - DQ(J,N, 0)
          P2   = DQ(J,N, 0) - DQ(J,N,-1)
          P4   = DQ(J,N,-1) - DQ(J,N,-2)
C         TMP1 = P2
C         TMP2 = P4
C         IF( ABS(P1) .LE. ABS(P2) ) TMP1= P1
C         IF( ABS(P2) .LE. ABS(P4) ) TMP2= P2
          TMP1 = 0.5*(SGN(P1)+SGN(P2))*MIN(ABS(P1),ABS(P2))
          TMP2 = 0.5*(SGN(P4)+SGN(P2))*MIN(ABS(P4),ABS(P2))
          Z1   = DQ(J,N, 0) - 0.5*TMP1
          Z2   = DQ(J,N,-1) + 0.5*TMP2
          GB(J,N)= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
   35     CONTINUE
C
C     B.C. for GBAR (zero order extrapolation)
C
          Z1        = DQ(JUP,N, 0)
          Z2        = DQ(JUP,N,-1)
          GB(JUP,N) = 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
          Z1        = DQ(JLOW,N, 0)
          Z2        = DQ(JLOW,N,-1)
          GB(JLOW,N)= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
          GB(JB,N)  = GB(JLOW,N)
          GB(JE,N)  = GB(JUP ,N)
   36     CONTINUE
          IF( IMTH .EQ. 13 .OR. IMTH .EQ. 23 ) THEN  ! ENO2
              DO 40 N = 1,NV
              DO 40 J = JLOW,JUP
              SLOPE_J = TK(J,N,1)*GB(J,1) +
     &                  TK(J,N,2)*GB(J,2) + TK(J,N,3)*GB(J,3)
              QR(J-1,N) = Q(J,N) - 0.5*SLOPE_J
              QL(J  ,N) = Q(J,N) + 0.5*SLOPE_J
   40         CONTINUE
          ELSE                                       ! FNO3
              DO 46 J = JLOW,JUP
              DO 45 N = 1,NV
              Z1      = DQ(J,N, 0)
              Z2      = DQ(J,N,-1)
              IF( Z1*Z2 .GT. 0. ) THEN
                  BP = MAX( 2.,GB(J,N)/Z1 )
                  BM = MAX( 2.,GB(J,N)/Z2 )
                  CP = (2.*BP-CKP)/CKM
                  CM = (2.*BM-CKP)/CKM
                  Z3 = SGN(Z2)*MIN( ABS(Z2),CP*ABS(Z1) )
                  Z4 = SGN(Z1)*MIN( ABS(Z1),CM*ABS(Z2) )
              ELSE
                  Z3 = SGN( ABS(GB(J,N)) )*Z2
                  Z4 = SGN( ABS(GB(J,N)) )*Z1
              END IF
              BB(N)= 0.25*( CKM*Z3 + CKP*Z4 )
              CC(N)= 0.25*( CKP*Z3 + CKM*Z4 )
   45         CONTINUE
              DO 46 N = 1,NV
              QL(J  ,N) = Q(J,N) + TK(J,N,1)*BB(1) + TK(J,N,2)*BB(2)
     &                           + TK(J,N,3)*BB(3)
              QR(J-1,N) = Q(J,N) - TK(J,N,1)*CC(1) - TK(J,N,2)*CC(2)
     &                           - TK(J,N,3)*CC(3)
   46         CONTINUE
          END IF
      ENDIF
C
      IF( IMTH .EQ. 15 .OR. IMTH .EQ. 25 ) THEN      ! ENO3
          DO 56 N = 1,NV
          DO 55 J = JLOW,JUP
          Z1      = DQ(J,N, 0)
          Z2      = DQ(J,N,-1)
          GB(J,N) = 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
   55     CONTINUE
          GB(JB,N)= GB(JLOW,N)
          GB(JE,N)= GB(JUP ,N)
   56     CONTINUE
          DO 61 N = 1,NV
          DO 60 J = JLOW+1,JUP-1
          Z1      = DQ(J,N, 0)
          Z2      = DQ(J,N,-1)
          P2      = Z1 - Z2
          IF( ABS(Z2) .GT. ABS(Z1) ) THEN
              P1    = DQ(J,N,1) - Z1
              TMP1  = P2
              IF( ABS(P1) .LE. ABS(P2) ) TMP1= P1
              GC(J,N)= -TMP1
          ELSE
              P4    = Z2 - DQ(J,N,-2)
              TMP2  = P4
              IF( ABS(P2) .LE. ABS(P4) ) TMP2= P2
              GC(J,N)= 2.*TMP2
          END IF
   60     CONTINUE
          GC(JLOW,N)= GC(JLOW+1,N)
          GC(JB  ,N)= GC(JLOW  ,N)
          GC(JUP ,N)= GC(JUP-1 ,N)
          GC(JE  ,N)= GC(JUP   ,N)
   61     CONTINUE
          DO 65 N = 1,NV
          DO 65 J = JLOW,JUP
          SLOPE_J = TK(J,N,1)*(0.5*GB(J,1) + GC(J,1)/6.) +
     &              TK(J,N,2)*(0.5*GB(J,2) + GC(J,2)/6.) +
     &              TK(J,N,3)*(0.5*GB(J,3) + GC(J,3)/6.)
          QR(J-1,N) = Q(J,N) - SLOPE_J
          QL(J  ,N) = Q(J,N) + SLOPE_J
   65     CONTINUE
      END IF
C
      IF( IMTH .EQ. 14 .OR. IMTH .EQ. 24 ) THEN  ! TVD2
          DO 76 N = 1,NV
          IF( IMOD .EQ. 0 ) THEN
              DO 70 J = JLOW,JUP
              Z1      = DQ(J,N, 0)
              Z2      = DQ(J,N,-1)
              GB(J,N)= Z1
              IF( ABS(Z2) .LE. ABS(Z1) ) GB(J,N)= Z2
   70         CONTINUE
          ELSEIF( IMOD .EQ. 1 ) THEN
              DO 71 J = JLOW,JUP
              Z1      = DQ(J,N, 0)
              Z2      = DQ(J,N,-1)
              GB(J,N)= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
   71         CONTINUE
          ELSEIF( IMOD .EQ. 2 ) THEN
              DO 72 J = JLOW,JUP
              Z1      = DQ(J,N, 0)
              Z2      = DQ(J,N,-1)
              GB(J,N)= (Z1*Z2+0.0005)*(Z1+Z2)/( Z1*Z1+Z2*Z2+0.001 )
   72         CONTINUE
          ELSEIF( IMOD .EQ. 3 ) THEN
              DO 73 J = JLOW,JUP
              Z1      = DQ(J,N, 0)
              Z2      = DQ(J,N,-1)
              Z3= ( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
              Z4= 0.5*( Z1+Z2 )
              GB(J,N)= 0.5*( SGN(Z3)+SGN(Z4) )*MIN( ABS(Z3),ABS(Z4) )
   73         CONTINUE
          ELSEIF( IMOD .EQ. 4 ) THEN
              DO 74 J = JLOW,JUP
              Z1      = DQ(J,N, 0)
              Z2      = DQ(J,N,-1)
              GB(J,N)= 0.
              IF( Z1+Z2 .NE. 0. ) GB(J,N)= (ABS(Z1*Z2)+Z1*Z2)/(Z1+Z2)
   74         CONTINUE
          ELSE
              DO 75 J = JLOW,JUP
              Z1      = DQ(J,N, 0)
              Z2      = DQ(J,N,-1)
              Z3= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( 2.*ABS(Z1),ABS(Z2) )
              Z4= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( 2.*ABS(Z2),ABS(Z1) )
              GB(J,N)= SGN(Z1)*MAX( ABS(Z3),ABS(Z4) )
   75         CONTINUE
          END IF
          GB(JB,N)= GB(JLOW,N)
          GB(JE,N)= GB(JUP ,N)
   76     CONTINUE
          DO 77 N = 1,NV
          DO 77 J = JLOW,JUP
          SLOPE_J = TK(J,N,1)*GB(J,1) +
     &              TK(J,N,2)*GB(J,2) + TK(J,N,3)*GB(J,3)
          QR(J-1,N) = Q(J,N) - 0.5*SLOPE_J
          QL(J  ,N) = Q(J,N) + 0.5*SLOPE_J
   77     CONTINUE
      END IF
C
      DO 80 J = JB,JUP
      UL      = QL(J,2)/QL(J,1)
      PL(J)   =(QL(J,3)-0.5*QL(J,1)*UL*UL)*GM1
      FL(J,1) = QL(J,2)
      FL(J,2) = FL(J,1)*UL + PL(J)
      FL(J,3) = ( QL(J,3) + PL(J) )*UL
C
      UR      = QR(J,2)/QR(J,1)
      PR(J)   =(QR(J,3)-0.5*QR(J,1)*UR*UR)*GM1
      FR(J,1) = QR(J,2)
      FR(J,2) = FR(J,1)*UR + PR(J)
      FR(J,3) = ( QR(J,3) + PR(J) )*UR
   80 CONTINUE
C
      IF( IMTH .LE. 20 ) GOTO 84
C
      DO 81 J = JB,JUP
      Q1      = QR(J,1) - QL(J,1)
      Q2      = QR(J,2) - QL(J,2)
      Q3      = QR(J,3) - QL(J,3)
C
      UL      = QL(J,2)/QL(J,1)
      UR      = QR(J,2)/QR(J,1)
      CL      = SQRT( PR(J)*GAMMA/QR(J,1) )
      CR      = SQRT( PL(J)*GAMMA/QL(J,1) )
      EMAX1   = 1.1*MAX( ABS(UL)    , ABS(UR)    )
      EMAX2   = 1.1*MAX( ABS(UL+CL) , ABS(UR+CR) )
      EMAX3   = 1.1*MAX( ABS(UL-CL) , ABS(UR-CR) )
C
C     EMAX    = MAX(0.7,EPSI)*( ABS(UL+UR) + CL+CR )*0.5
C
      BET(J,1)= Q1*EMAX1
      BET(J,2)= Q2*EMAX2
      BET(J,3)= Q3*EMAX3
   81 CONTINUE
      GOTO 99
C
   84 DO 85 J = JB,JUP
      RBAR    = SQRT( QR(J,1)/QL(J,1) )
      ORBAR   = 1./(1. + RBAR)
      UA(J)   = ( QR(J,2)/QR(J,1)*RBAR + QL(J,2)/QL(J,1) )*ORBAR
      HJR     = (PR(J)+QR(J,3))/QR(J,1)
      HJL     = (PL(J)+QL(J,3))/QL(J,1)
      HB      = ( HJR*RBAR + HJL )*ORBAR
      CSQ     = ( HB - 0.5*UA(J)*UA(J) )*GM1
      CA(J)   = SQRT(CSQ)
   85 CONTINUE
C
      DO 90 J = JB,JUP
      EIG(J,1)= PHI( UA(J),EPSI )
      EIG(J,2)= PHI( UA(J)+CA(J),EPSI )
      EIG(J,3)= PHI( UA(J)-CA(J),EPSI )
   90 CONTINUE
C
      CALL RDIGL(JB,JUP,UA,CA,EIG,TT,MAXI,GAMMA)
C
      DO 95 J = JB,JUP
      Q1      = QR(J,1) - QL(J,1)
      Q2      = QR(J,2) - QL(J,2)
      Q3      = QR(J,3) - QL(J,3)
      BET(J,1)= TT(J,1,1)*Q1+TT(J,1,2)*Q2+TT(J,1,3)*Q3
      BET(J,2)= TT(J,2,1)*Q1+TT(J,2,2)*Q2+TT(J,2,3)*Q3
      BET(J,3)= TT(J,3,1)*Q1+TT(J,3,2)*Q2+TT(J,3,3)*Q3
   95 CONTINUE
C
C     CALCULATE NUMERICAL FLUX
C
   99 DO 100 N = 1,NV
      DO 100 J = JB,JUP
      FN(J,N)  = 0.5*( FL(J,N) + FR(J,N) - BET(J,N) )
  100 CONTINUE
      DO 110 N = 1,NV
      DO 110 J = JLOW,JUP
      S(J,N)   = FN(J-1,N) - FN(J,N)
  110 CONTINUE
      RETURN
      END
C
      SUBROUTINE MUSCL2( JB,JE,CK,OMIGA,Q,S,IMTH,IMOD )
      PARAMETER ( MAXI=407,NV=3 )
      COMMON /A1/ RQ(MAXI),UQ(MAXI),PQ(MAXI),EQ(MAXI),CQ(MAXI),HQ(MAXI)
      DIMENSION FN(MAXI,NV),BET(MAXI,NV),EIG(MAXI,NV),TT(MAXI,NV,NV)
      DIMENSION UA(MAXI),CA(MAXI),Q(MAXI,NV),S(MAXI,NV)
      DIMENSION QR(MAXI,NV),QL(MAXI,NV),FR(MAXI,NV),FL(MAXI,NV)
      DIMENSION QQ(MAXI,NV),EL(MAXI),ER(MAXI)
      DIMENSION DQ(MAXI,NV),GB(MAXI,NV),GC(MAXI,NV)
      COMMON /D1/ GAMMA,GM1,EPSI,LEVEL,DTX,IPROB
C
C     LIMITER ON PRIMITIVE VARIABLES
C
      CKM   = 1. - CK
      CKP   = 1. + CK
      JUP   = JE - 1
      JLOW  = JB + 1
C
      DO 10 J = JB,JE
      QQ(J,1) = RQ(J)
      QQ(J,2) = UQ(J)
      QQ(J,3) = PQ(J)
   10 CONTINUE
      DO 15 N = 1,NV
      QL(JB,N)= QQ(JB,N)
      QR(JB,N)= QQ(JB,N)
      QL(JE,N)= QQ(JE,N)
      QR(JE,N)= QQ(JE,N)
      QL(JLOW,N) = QQ(JLOW,N)
      QR(JLOW,N) = QQ(JLOW+1,N)
      QL(JUP ,N) = QQ(JUP,N)
      QR(JUP ,N) = QQ(JE,N)
   15 CONTINUE
      DO 20 J = JB,JUP
      DQ(J,1) = RQ(J+1) - RQ(J)
      DQ(J,2) = UQ(J+1) - UQ(J)
      DQ(J,3) = PQ(J+1) - PQ(J)
   20 CONTINUE
C
      IF( IMTH .EQ. 11 .OR. IMTH .EQ. 21 ) THEN
          DO 25 N = 1,NV
          DO 25 J = JLOW,JUP
          Z1      = SGN( DQ(J-1,N) )/2.
          Z2      = SGN( DQ( J ,N) )/2.
          Z3      = (Z1+Z2)*MIN( ABS(DQ(J-1,N)),ABS(OMIGA*DQ(J,N)) )
          Z4      = (Z1+Z2)*MIN( ABS(OMIGA*DQ(J-1,N)),ABS(DQ(J,N)) )
          QL(J  ,N) = QQ(J,N) + .25*( CKM*Z3 + CKP*Z4 )
          QR(J-1,N) = QQ(J,N) - .25*( CKP*Z3 + CKM*Z4 )
   25     CONTINUE
      END IF
C
      IF( IMTH .EQ. 12 .OR. IMTH .EQ. 13 .OR.
     &    IMTH .EQ. 22 .OR. IMTH .EQ. 23 ) THEN
C
C     CALCULATE MINMOD FOR GBAR
C
          DO 36 N= 1,NV
          DO 35 J= JLOW+1,JUP-1
          P1   = DQ(J+1,N) - DQ( J ,N)
          P2   = DQ( J ,N) - DQ(J-1,N)
          P4   = DQ(J-1,N) - DQ(J-2,N)
C         TMP1 = P2
C         TMP2 = P4
C         IF( ABS(P1) .LE. ABS(P2) ) TMP1= P1
C         IF( ABS(P2) .LE. ABS(P4) ) TMP2= P2
          TMP1 = 0.5*(SGN(P1)+SGN(P2))*MIN(ABS(P1),ABS(P2))
          TMP2 = 0.5*(SGN(P4)+SGN(P2))*MIN(ABS(P4),ABS(P2))
          Z1   = DQ( J ,N) - 0.5*TMP1
          Z2   = DQ(J-1,N) + 0.5*TMP2
          GB(J,N)= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
   35     CONTINUE
C
C     B.C. for GBAR (zero order extrapolation)
C
          Z1        = DQ(JUP  ,N)
          Z2        = DQ(JUP-1,N)
          GB(JUP,N) = 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
          Z1        = DQ(JLOW,N)
          Z2        = DQ(JB  ,N)
          GB(JLOW,N)= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
          GB(JB,N)  = GB(JLOW,N)
          GB(JE,N)  = GB(JUP ,N)
   36     CONTINUE
          IF( IMTH .EQ. 13 .OR. IMTH .EQ. 23 ) THEN
              DO 40 N = 1,NV
              DO 40 J = JLOW,JUP
              QR(J-1,N) = QQ(J,N) - 0.5*GB(J,N)
              QL(J  ,N) = QQ(J,N) + 0.5*GB(J,N)
   40         CONTINUE
          ELSE
              DO 45 N = 1,NV
              DO 45 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              IF( Z1*Z2 .GT. 0. ) THEN
                  BP = MAX( 2.,GB(J,N)/Z1 )
                  BM = MAX( 2.,GB(J,N)/Z2 )
                  CP = (2.*BP-CKP)/CKM
                  CM = (2.*BM-CKP)/CKM
                  Z3 = SGN(Z2)*MIN( ABS(Z2),CP*ABS(Z1) )
                  Z4 = SGN(Z1)*MIN( ABS(Z1),CM*ABS(Z2) )
              ELSE
                  Z3 = SGN( ABS(GB(J,N)) )*Z2
                  Z4 = SGN( ABS(GB(J,N)) )*Z1
              END IF
              QL(J  ,N)= QQ(J,N) + .25*( CKM*Z3 + CKP*Z4 )
              QR(J-1,N)= QQ(J,N) - .25*( CKP*Z3 + CKM*Z4 )
   45         CONTINUE
          END IF
      ENDIF
C
      IF( IMTH .EQ. 15 .OR. IMTH .EQ. 25 ) THEN
          DO 56 N = 1,NV
          DO 55 J = JLOW,JUP
          Z1      = DQ( J ,N)
          Z2      = DQ(J-1,N)
          GB(J,N) = 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
   55     CONTINUE
          GB(JB,N)= GB(JLOW,N)
          GB(JE,N)= GB(JUP ,N)
   56     CONTINUE
          DO 61 N = 1,NV
          DO 60 J = JLOW+1,JUP-1
          Z1      = DQ( J ,N)
          Z2      = DQ(J-1,N)
          P2      = Z1 - Z2
          IF( ABS(Z2) .GT. ABS(Z1) ) THEN
              P1    = DQ(J+1,N) - Z1
              TMP1  = P2
              IF( ABS(P1) .LE. ABS(P2) ) TMP1= P1
              GC(J,N)= -TMP1
          ELSE
              P4    = Z2 - DQ(J-2,N)
              TMP2  = P4
              IF( ABS(P2) .LE. ABS(P4) ) TMP2= P2
              GC(J,N)= 2.*TMP2
          END IF
   60     CONTINUE
          GC(JLOW,N)= GC(JLOW+1,N)
          GC(JB  ,N)= GC(JLOW  ,N)
          GC(JUP ,N)= GC(JUP-1 ,N)
          GC(JE  ,N)= GC(JUP   ,N)
   61     CONTINUE
          DO 65 N = 1,NV
          DO 65 J = JLOW,JUP
          QR(J-1,N) = QQ(J,N) - 0.5*GB(J,N) - GC(J,N)/6.
          QL(J  ,N) = QQ(J,N) + 0.5*GB(J,N) + GC(J,N)/6.
   65     CONTINUE
      END IF
C
      IF( IMTH .EQ. 14 .OR. IMTH .EQ. 24 ) THEN
          DO 76 N = 1,NV
          IF( IMOD .EQ. 0 ) THEN
              DO 70 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              GB(J,N)= Z1
              IF( ABS(Z2) .LE. ABS(Z1) ) GB(J,N)= Z2
   70         CONTINUE
          ELSEIF( IMOD .EQ. 1 ) THEN
              DO 71 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              GB(J,N)= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
   71         CONTINUE
          ELSEIF( IMOD .EQ. 2 ) THEN
              DO 72 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              GB(J,N)= (Z1*Z2+0.0005)*(Z1+Z2)/( Z1*Z1+Z2*Z2+0.001 )
   72         CONTINUE
          ELSEIF( IMOD .EQ. 3 ) THEN
              DO 73 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              Z3= ( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
              Z4= 0.5*( Z1+Z2 )
              GB(J,N)= 0.5*( SGN(Z3)+SGN(Z4) )*MIN( ABS(Z3),ABS(Z4) )
   73         CONTINUE
          ELSEIF( IMOD .EQ. 4 ) THEN
              DO 74 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              GB(J,N)= 0.
              IF( Z1+Z2 .NE. 0. ) GB(J,N)= (ABS(Z1*Z2)+Z1*Z2)/(Z1+Z2)
   74         CONTINUE
          ELSE
              DO 75 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              Z3= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( 2.*ABS(Z1),ABS(Z2) )
              Z4= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( 2.*ABS(Z2),ABS(Z1) )
              GB(J,N)= SGN(Z1)*MAX( ABS(Z3),ABS(Z4) )
   75         CONTINUE
          END IF
          GB(JB,N)= GB(JLOW,N)
          GB(JE,N)= GB(JUP ,N)
   76     CONTINUE
          DO 77 N = 1,NV
          DO 77 J = JLOW,JUP
          QR(J-1,N) = QQ(J,N) - 0.5*GB(J,N)
          QL(J  ,N) = QQ(J,N) + 0.5*GB(J,N)
   77     CONTINUE
      END IF
C
      DO 80 J = JB,JUP
      EL(J)   = QL(J,3)/GM1 + 0.5*QL(J,1)*QL(J,2)*QL(J,2)
      FL(J,1) = QL(J,1)*QL(J,2)
      FL(J,2) = FL(J,1)*QL(J,2) + QL(J,3)
      FL(J,3) = ( EL(J) + QL(J,3) )*QL(J,2)
C
      ER(J)   = QR(J,3)/GM1 + 0.5*QR(J,1)*QR(J,2)*QR(J,2)
      FR(J,1) = QR(J,1)*QR(J,2)
      FR(J,2) = FR(J,1)*QR(J,2) + QR(J,3)
      FR(J,3) = ( ER(J) + QR(J,3) )*QR(J,2)
   80 CONTINUE
C
      IF( IMTH .LE. 20 ) GOTO 84
C
      DO 81 J = JB,JUP
      Q1      = QR(J,1) - QL(J,1)
      Q2      = QR(J,2)*QR(J,1) - QL(J,2)*QL(J,1)
      Q3      = ER(J) - EL(J)
C
      UL      = QL(J,2)
      UR      = QR(J,2)
      CL      = SQRT( QR(J,3)*GAMMA/QR(J,1) )
      CR      = SQRT( QL(J,3)*GAMMA/QL(J,1) )
      EMAX1   = 1.1*MAX( ABS(UL)    , ABS(UR)    )
      EMAX2   = 1.1*MAX( ABS(UL+CL) , ABS(UR+CR) )
      EMAX3   = 1.1*MAX( ABS(UL-CL) , ABS(UR-CR) )
C
C     EMAX    = MAX(0.7,EPSI)*( ABS(UL+UR) + CL+CR )*0.5
C
      BET(J,1)= Q1*EMAX1
      BET(J,2)= Q2*EMAX2
      BET(J,3)= Q3*EMAX3
   81 CONTINUE
      GOTO 99
C
   84 DO 85 J = JB,JUP
      RBAR    = SQRT( QR(J,1)/QL(J,1) )
      ORBAR   = 1./(1. + RBAR)
      UA(J)   = ( QR(J,2)*RBAR + QL(J,2) )*ORBAR
      HJR     = (ER(J)+QR(J,3))/QR(J,1)
      HJL     = (EL(J)+QL(J,3))/QL(J,1)
      HB      = ( HJR*RBAR + HJL )*ORBAR
      CSQ     = ( HB - 0.5*UA(J)*UA(J) )*GM1
      CA(J)   = SQRT(CSQ)
   85 CONTINUE
C
      DO 90 J = JB,JUP
      EIG(J,1)= PHI( UA(J),EPSI )
      EIG(J,2)= PHI( UA(J)+CA(J),EPSI )
      EIG(J,3)= PHI( UA(J)-CA(J),EPSI )
   90 CONTINUE
C
      CALL RDIGL(JB,JUP,UA,CA,EIG,TT,MAXI,GAMMA)
C
      DO 95 J = JB,JUP
      Q1      = QR(J,1) - QL(J,1)
      Q2      = QR(J,2)*QR(J,1) - QL(J,2)*QL(J,1)
      Q3      = ER(J) - EL(J)
      BET(J,1)= TT(J,1,1)*Q1+TT(J,1,2)*Q2+TT(J,1,3)*Q3
      BET(J,2)= TT(J,2,1)*Q1+TT(J,2,2)*Q2+TT(J,2,3)*Q3
      BET(J,3)= TT(J,3,1)*Q1+TT(J,3,2)*Q2+TT(J,3,3)*Q3
   95 CONTINUE
C
C     CALCULATE NUMERICAL FLUX
C
   99 DO 100 N = 1,NV
      DO 100 J = JB,JUP
      FN(J,N)  = 0.5*( FL(J,N) + FR(J,N) - BET(J,N) )
  100 CONTINUE
      DO 110 N = 1,NV
      DO 110 J = JLOW,JUP
      S(J,N)   = FN(J-1,N) - FN(J,N)
  110 CONTINUE
      RETURN
      END
C
      SUBROUTINE MUSCL3( JB,JE,CK,OMIGA,Q,S,IMTH,IMOD )
      PARAMETER ( MAXI=407,NV=3 )
      COMMON /A1/ RQ(MAXI),UQ(MAXI),PQ(MAXI),EQ(MAXI),CQ(MAXI),HQ(MAXI)
      DIMENSION FN(MAXI,NV),BET(MAXI,NV),EIG(MAXI,NV),TT(MAXI,NV,NV)
      DIMENSION UA(MAXI),CA(MAXI),Q(MAXI,NV),S(MAXI,NV)
      DIMENSION QR(MAXI,NV),QL(MAXI,NV),FR(MAXI,NV),FL(MAXI,NV)
      DIMENSION PL(MAXI),PR(MAXI),DQ(MAXI,NV),GB(MAXI,NV),GC(MAXI,NV)
      COMMON /D1/ GAMMA,GM1,EPSI,LEVEL,DTX,IPROB
C
C     LIMITER ON CONSERVATIVE VARIABLES
C
      CKM   = 1. - CK
      CKP   = 1. + CK
      JUP   = JE - 1
      JLOW  = JB + 1
C
      DO 10 N = 1,NV
      QL(JB,N)= Q(JB,N)
      QR(JB,N)= Q(JB,N)
      QL(JE,N)= Q(JE,N)
      QR(JE,N)= Q(JE,N)
      QL(JLOW,N) = Q(JLOW,N)
      QR(JLOW,N) = Q(JLOW+1,N)
      QL(JUP ,N) = Q(JUP,N)
      QR(JUP ,N) = Q(JE,N)
   10 CONTINUE
      DO 20 J = JB,JUP
      DQ(J,1) = Q(J+1,1) - Q(J,1)
      DQ(J,2) = Q(J+1,2) - Q(J,2)
      DQ(J,3) = Q(J+1,3) - Q(J,3)
   20 CONTINUE
C
      IF( IMTH .EQ. 11 .OR. IMTH .EQ. 21 ) THEN
          DO 25 N = 1,NV
          DO 25 J = JLOW,JUP
          Z1      = SGN( DQ(J-1,N) )/2.
          Z2      = SGN( DQ( J ,N) )/2.
          Z3      = (Z1+Z2)*MIN( ABS(DQ(J-1,N)),ABS(OMIGA*DQ(J,N)) )
          Z4      = (Z1+Z2)*MIN( ABS(OMIGA*DQ(J-1,N)),ABS(DQ(J,N)) )
          QL(J  ,N) = Q(J,N) + .25*( CKM*Z3 + CKP*Z4 )
          QR(J-1,N) = Q(J,N) - .25*( CKP*Z3 + CKM*Z4 )
   25     CONTINUE
      END IF
C
      IF( IMTH .EQ. 12 .OR. IMTH .EQ. 13 .OR.
     &    IMTH .EQ. 22 .OR. IMTH .EQ. 23 ) THEN
C
C     CALCULATE MINMOD FOR GBAR
C
          DO 36 N= 1,NV
          DO 35 J= JLOW+1,JUP-1
          P1   = DQ(J+1,N) - DQ( J ,N)
          P2   = DQ( J ,N) - DQ(J-1,N)
          P4   = DQ(J-1,N) - DQ(J-2,N)
C         TMP1 = P2
C         TMP2 = P4
C         IF( ABS(P1) .LE. ABS(P2) ) TMP1= P1
C         IF( ABS(P2) .LE. ABS(P4) ) TMP2= P2
          TMP1 = 0.5*(SGN(P1)+SGN(P2))*MIN(ABS(P1),ABS(P2))
          TMP2 = 0.5*(SGN(P4)+SGN(P2))*MIN(ABS(P4),ABS(P2))
          Z1   = DQ( J ,N) - 0.5*TMP1
          Z2   = DQ(J-1,N) + 0.5*TMP2
          GB(J,N)= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
   35     CONTINUE
C
C     B.C. for GBAR (zero order extrapolation)
C
          Z1        = DQ(JUP  ,N)
          Z2        = DQ(JUP-1,N)
          GB(JUP,N) = 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
          Z1        = DQ(JLOW,N)
          Z2        = DQ(JB  ,N)
          GB(JLOW,N)= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
          GB(JB,N)  = GB(JLOW,N)
          GB(JE,N)  = GB(JUP ,N)
   36     CONTINUE
          IF( IMTH .EQ. 13 .OR. IMTH .EQ. 23 ) THEN
              DO 40 N = 1,NV
              DO 40 J = JLOW,JUP
              QR(J-1,N) = Q(J,N) - 0.5*GB(J,N)
              QL(J  ,N) = Q(J,N) + 0.5*GB(J,N)
   40         CONTINUE
          ELSE
              DO 45 N = 1,NV
              DO 45 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              IF( Z1*Z2 .GT. 0. ) THEN
                  BP = MAX( 2.,GB(J,N)/Z1 )
                  BM = MAX( 2.,GB(J,N)/Z2 )
                  CP = (2.*BP-CKP)/CKM
                  CM = (2.*BM-CKP)/CKM
                  Z3 = SGN(Z2)*MIN( ABS(Z2),CP*ABS(Z1) )
                  Z4 = SGN(Z1)*MIN( ABS(Z1),CM*ABS(Z2) )
              ELSE
                  Z3 = SGN( ABS(GB(J,N)) )*Z2
                  Z4 = SGN( ABS(GB(J,N)) )*Z1
              END IF
              QL(J  ,N)= Q(J,N) + .25*( CKM*Z3 + CKP*Z4 )
              QR(J-1,N)= Q(J,N) - .25*( CKP*Z3 + CKM*Z4 )
   45         CONTINUE
          END IF
      ENDIF
C
      IF( IMTH .EQ. 15 .OR. IMTH .EQ. 25 ) THEN
          DO 56 N = 1,NV
          DO 55 J = JLOW,JUP
          Z1      = DQ( J ,N)
          Z2      = DQ(J-1,N)
          GB(J,N) = 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
   55     CONTINUE
          GB(JB,N)= GB(JLOW,N)
          GB(JE,N)= GB(JUP ,N)
   56     CONTINUE
          DO 61 N = 1,NV
          DO 60 J = JLOW+1,JUP-1
          Z1      = DQ( J ,N)
          Z2      = DQ(J-1,N)
          P2      = Z1 - Z2
          IF( ABS(Z2) .GT. ABS(Z1) ) THEN
              P1    = DQ(J+1,N) - Z1
              TMP1  = P2
              IF( ABS(P1) .LE. ABS(P2) ) TMP1= P1
              GC(J,N)= -TMP1
          ELSE
              P4    = Z2 - DQ(J-2,N)
              TMP2  = P4
              IF( ABS(P2) .LE. ABS(P4) ) TMP2= P2
              GC(J,N)= 2.*TMP2
          END IF
   60     CONTINUE
          GC(JLOW,N)= GC(JLOW+1,N)
          GC(JB  ,N)= GC(JLOW  ,N)
          GC(JUP ,N)= GC(JUP-1 ,N)
          GC(JE  ,N)= GC(JUP   ,N)
   61     CONTINUE
          DO 65 N = 1,NV
          DO 65 J = JLOW,JUP
          QR(J-1,N) = Q(J,N) - 0.5*GB(J,N) - GC(J,N)/6.
          QL(J  ,N) = Q(J,N) + 0.5*GB(J,N) + GC(J,N)/6.
   65     CONTINUE
      END IF
C
      IF( IMTH .EQ. 14 .OR. IMTH .EQ. 24 ) THEN
          DO 76 N = 1,NV
          IF( IMOD .EQ. 0 ) THEN
              DO 70 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              GB(J,N)= Z1
              IF( ABS(Z2) .LE. ABS(Z1) ) GB(J,N)= Z2
   70         CONTINUE
          ELSEIF( IMOD .EQ. 1 ) THEN
              DO 71 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              GB(J,N)= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
   71         CONTINUE
          ELSEIF( IMOD .EQ. 2 ) THEN
              DO 72 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              GB(J,N)= (Z1*Z2+0.0005)*(Z1+Z2)/( Z1*Z1+Z2*Z2+0.001 )
   72         CONTINUE
          ELSEIF( IMOD .EQ. 3 ) THEN
              DO 73 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              Z3= ( SGN(Z1)+SGN(Z2) )*MIN( ABS(Z1),ABS(Z2) )
              Z4= 0.5*( Z1+Z2 )
              GB(J,N)= 0.5*( SGN(Z3)+SGN(Z4) )*MIN( ABS(Z3),ABS(Z4) )
   73         CONTINUE
          ELSEIF( IMOD .EQ. 4 ) THEN
              DO 74 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              GB(J,N)= 0.
              IF( Z1+Z2 .NE. 0. ) GB(J,N)= (ABS(Z1*Z2)+Z1*Z2)/(Z1+Z2)
   74         CONTINUE
          ELSE
              DO 75 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              Z3= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( 2.*ABS(Z1),ABS(Z2) )
              Z4= 0.5*( SGN(Z1)+SGN(Z2) )*MIN( 2.*ABS(Z2),ABS(Z1) )
              GB(J,N)= SGN(Z1)*MAX( ABS(Z3),ABS(Z4) )
   75         CONTINUE
          END IF
          GB(JB,N)= GB(JLOW,N)
          GB(JE,N)= GB(JUP ,N)
   76     CONTINUE
          DO 77 N = 1,NV
          DO 77 J = JLOW,JUP
          QR(J-1,N) = Q(J,N) - 0.5*GB(J,N)
          QL(J  ,N) = Q(J,N) + 0.5*GB(J,N)
   77     CONTINUE
      END IF
C
      DO 80 J = JB,JUP
      UL      = QL(J,2)/QL(J,1)
      PL(J)   =(QL(J,3)-0.5*QL(J,1)*UL*UL)*GM1
      FL(J,1) = QL(J,2)
      FL(J,2) = FL(J,1)*UL + PL(J)
      FL(J,3) = ( QL(J,3) + PL(J) )*UL
C
      UR      = QR(J,2)/QR(J,1)
      PR(J)   =(QR(J,3)-0.5*QR(J,1)*UR*UR)*GM1
      FR(J,1) = QR(J,2)
      FR(J,2) = FR(J,1)*UR + PR(J)
      FR(J,3) = ( QR(J,3) + PR(J) )*UR
   80 CONTINUE
C
      IF( IMTH .LE. 20 ) GOTO 84
C
      DO 81 J = JB,JUP
      Q1      = QR(J,1) - QL(J,1)
      Q2      = QR(J,2) - QL(J,2)
      Q3      = QR(J,3) - QL(J,3)
C
      UL      = QL(J,2)/QL(J,1)
      UR      = QR(J,2)/QR(J,1)
      CL      = SQRT( PR(J)*GAMMA/QR(J,1) )
      CR      = SQRT( PL(J)*GAMMA/QL(J,1) )
      EMAX1   = 1.1*MAX( ABS(UL)    , ABS(UR)    )
      EMAX2   = 1.1*MAX( ABS(UL+CL) , ABS(UR+CR) )
      EMAX3   = 1.1*MAX( ABS(UL-CL) , ABS(UR-CR) )
C
C     EMAX    = MAX(0.7,EPSI)*( ABS(UL+UR) + CL+CR )*0.5
C
      BET(J,1)= Q1*EMAX1
      BET(J,2)= Q2*EMAX2
      BET(J,3)= Q3*EMAX3
   81 CONTINUE
      GOTO 99
C
   84 DO 85 J = JB,JUP
      RBAR    = SQRT( QR(J,1)/QL(J,1) )
      ORBAR   = 1./(1. + RBAR)
      UA(J)   = ( QR(J,2)/QR(J,1)*RBAR + QL(J,2)/QL(J,1) )*ORBAR
      HJR     = (PR(J)+QR(J,3))/QR(J,1)
      HJL     = (PL(J)+QL(J,3))/QL(J,1)
      HB      = ( HJR*RBAR + HJL )*ORBAR
      CSQ     = ( HB - 0.5*UA(J)*UA(J) )*GM1
      CA(J)   = SQRT(CSQ)
   85 CONTINUE
C
      DO 90 J = JB,JUP
      EIG(J,1)= PHI( UA(J),EPSI )
      EIG(J,2)= PHI( UA(J)+CA(J),EPSI )
      EIG(J,3)= PHI( UA(J)-CA(J),EPSI )
   90 CONTINUE
C
      CALL RDIGL(JB,JUP,UA,CA,EIG,TT,MAXI,GAMMA)
C
      DO 95 J = JB,JUP
      Q1      = QR(J,1) - QL(J,1)
      Q2      = QR(J,2) - QL(J,2)
      Q3      = QR(J,3) - QL(J,3)
      BET(J,1)= TT(J,1,1)*Q1+TT(J,1,2)*Q2+TT(J,1,3)*Q3
      BET(J,2)= TT(J,2,1)*Q1+TT(J,2,2)*Q2+TT(J,2,3)*Q3
      BET(J,3)= TT(J,3,1)*Q1+TT(J,3,2)*Q2+TT(J,3,3)*Q3
   95 CONTINUE
C
C     CALCULATE NUMERICAL FLUX
C
   99 DO 100 N = 1,NV
      DO 100 J = JB,JUP
      FN(J,N)  = 0.5*( FL(J,N) + FR(J,N) - BET(J,N) )
  100 CONTINUE
      DO 110 N = 1,NV
      DO 110 J = JLOW,JUP
      S(J,N)   = FN(J-1,N) - FN(J,N)
  110 CONTINUE
      RETURN
      END
C
      FUNCTION SGN(Z)
      IF( Z .GT. 0. ) THEN
          SGN = 1.
      ELSEIF( Z .LT. 0. ) THEN
          SGN = -1.
      ELSE
          SGN = 0.
      END IF
      RETURN
      END
C
      SUBROUTINE RDIGL(I1,I2,UA,CA,EIG,AT,MAXI,GAMMA)
      DIMENSION AT(MAXI,3,3),UA(MAXI),CA(MAXI),EIG(MAXI,3)
C
      GM1 = GAMMA - 1.
      DO 10 I = I1,I2
      E1 = EIG(I,1)
      E2 = 0.5*( EIG(I,2) + EIG(I,3) )
      E3 = 0.5*( EIG(I,2) - EIG(I,3) )
      E4 = E2 - E1
      UB = UA(I)
      Q22= 0.5*UB*UB
      CB = CA(I)
      CC = CB*CB
C
      A1 = E3/CB
      A2 = GM1*E4/CC
      B1 = UB*A1
      B2 = Q22*A1
      C1 = UB*A2
      C2 = Q22*A2
      D1 = E3*CB/GM1
C
      AT(I,1,1) = E1 - B1 + C2
      AT(I,1,2) = A1 - C1
      AT(I,1,3) = A2
C
      AT(I,2,1) = AT(I,1,1)*UB - UB*E2 + GM1*B2
      AT(I,2,2) = AT(I,1,2)*UB + E2 - GM1*B1
      AT(I,2,3) = AT(I,1,3)*UB + GM1*A1
C
      AT(I,3,1) = AT(I,1,1)*Q22 - E2*Q22 + GM1*B1*Q22 - UB*D1
      AT(I,3,2) = AT(I,1,2)*Q22 + D1 - GM1*B1*UB
      AT(I,3,3) = AT(I,1,3)*Q22 + E2 + GM1*B1
C
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE CENTR(JJ,IS,IE,Q,S,IMTH)
      PARAMETER ( MAXI=407,NV=3 )
      COMMON /A1/ RQ(MAXI),UQ(MAXI),PQ(MAXI),EQ(MAXI),CQ(MAXI),HQ(MAXI)
      DIMENSION F(MAXI,NV),Q(MAXI,NV),S(MAXI,NV)
      I1 = JJ + IS
C
C     FILL FLUX VECTOR
C
      IF( IMTH .NE. 42 ) THEN
          DO 10 I= IS,IE
          F(I,1) = RQ(I)*UQ(I)
          F(I,2) = RQ(I)*UQ(I)*UQ(I) + PQ(I)
          F(I,3) = UQ(I)*( Q(I,3) + PQ(I) )
   10     CONTINUE
      ELSE
          N = 2*IS - I1
          DO 20 I= I1,IE
          N = N + 1
          QA= UQ(I)
          QB= UQ(N)
          IF( (QA-QB)/(I-N) .GT. 0. ) QA = 0.5*(QA+QB)
          F(I,1) = RQ(I)*QA
          F(I,2) = RQ(I)*QA*UQ(I) + PQ(I)
          F(I,3) = QA*( Q(I,3) + PQ(I) )
   20     CONTINUE
      END IF
      M = IS
      DO 90 I= I1+1,IE
      M = M + 1
      S(M,1) = F(I-1,1) - F(I,1)
      S(M,2) = F(I-1,2) - F(I,2)
      S(M,3) = F(I-1,3) - F(I,3)
   90 CONTINUE
      RETURN
      END
C
      SUBROUTINE TFOUR(JJ,IS,IE,Q,S,IMTH)
      PARAMETER ( MAXI=407,NV=3 )
      COMMON /A1/ RQ(MAXI),UQ(MAXI),PQ(MAXI),EQ(MAXI),CQ(MAXI),HQ(MAXI)
      DIMENSION F(MAXI,NV),Q(MAXI,NV),S(MAXI,NV)
C
C     FILL FLUX VECTOR
C
      DO 10 I= IS,IE
      F(I,1) = RQ(I)*UQ(I)
      F(I,2) = RQ(I)*UQ(I)*UQ(I) + PQ(I)
      F(I,3) = UQ(I)*( Q(I,3) + PQ(I) )
   10 CONTINUE
C
      IF( JJ .EQ. 0 ) THEN
          S(IS+1,1) = F(IS,1) - F(IS+1,1)
          S(IS+1,2) = F(IS,2) - F(IS+1,2)
          S(IS+1,3) = F(IS,3) - F(IS+1,3)
          DO I = IS+2,IE-1
          S(I,1) =  - ( 7.*F(I,1) -8.*F(I-1,1) +F(I-2,1) )/6.
          S(I,2) =  - ( 7.*F(I,2) -8.*F(I-1,2) +F(I-2,2) )/6.
          S(I,3) =  - ( 7.*F(I,3) -8.*F(I-1,3) +F(I-2,3) )/6.
          ENDDO
      ELSE
          S(IE-1,1) = F(IE-1,1) - F(IE,1)
          S(IE-1,2) = F(IE-1,2) - F(IE,2)
          S(IE-1,3) = F(IE-1,3) - F(IE,3)
          DO I = IS+1,IE-2
          S(I,1) = ( 7.*F(I,1) -8.*F(I+1,1) +F(I+2,1) )/6.
          S(I,2) = ( 7.*F(I,2) -8.*F(I+1,2) +F(I+2,2) )/6.
          S(I,3) = ( 7.*F(I,3) -8.*F(I+1,3) +F(I+2,3) )/6.
          ENDDO
      END IF
      RETURN
      END
C
      SUBROUTINE DAMPER(IB,IE,Q,S,IMTH,IMOD)
      PARAMETER ( MAXI=407,NV=3 )
      COMMON /A1/ RQ(MAXI),UQ(MAXI),PQ(MAXI),EQ(MAXI),CQ(MAXI),HQ(MAXI)
      DIMENSION ALFA(MAXI,NV),EIG(MAXI,NV),FN(MAXI,NV),BET(MAXI,NV)
      DIMENSION UA(MAXI),CA(MAXI),HA(MAXI),RA(MAXI),Q(MAXI,NV)
      DIMENSION GB(MAXI,NV),TK(MAXI,NV,NV),TKI(MAXI,NV,NV),S(MAXI,NV)
      COMMON /D1/ GAMMA,GM1,EPSI,LEVEL,DTX,IPROB
C
      NP = IE-1
C
C     CALCULATE rho, u, p,... at i+1/2
C
C --- by ROE'S AVERAGE ---
C
      DO 10 I=IB,NP
      DD    = SQRT( RQ(I+1)/RQ(I) )
      RA(I) = DD*RQ(I)
      UA(I) = (DD*UQ(I+1)+UQ(I))/(1.+DD)
      HA(I) = (DD*HQ(I+1)+HQ(I))/(1.+DD)
      CA(I) = SQRT( GM1*(HA(I)-0.5*UA(I)**2) )
   10 CONTINUE
C
C     CALCULATE EIGENVALUE AT I+1/2
C
      DO 15 I = IB,NP
      EIG(I,1)= UA(I)*DTX
      EIG(I,2)= ( UA(I)+CA(I) )*DTX
      EIG(I,3)= ( UA(I)-CA(I) )*DTX
   15 CONTINUE
C
C     CALCULATE AFAR  AT I+1/2
C
      CALL RLVECT(MAXI,NP,UA,CA,HA,RA,TK,TKI)
C
      DO 20 I = IB,NP
      Q1 = Q(I+1,1) - Q(I,1)
      Q2 = Q(I+1,2) - Q(I,2)
      Q3 = Q(I+1,3) - Q(I,3)
      ALFA(I,1)= TKI(I,1,1)*Q1 + TKI(I,1,2)*Q2 + TKI(I,1,3)*Q3
      ALFA(I,2)= TKI(I,2,1)*Q1 + TKI(I,2,2)*Q2 + TKI(I,2,3)*Q3
      ALFA(I,3)= TKI(I,3,1)*Q1 + TKI(I,3,2)*Q2 + TKI(I,3,3)*Q3
   20 CONTINUE
C
C     2nd ORDER SYMMETRIC TVD
C
      CALL QFUNCT(IMOD,IB,IE,ALFA,GB)
C
      DO 25 K = 1 ,NV
      DO 25 I = IB,NP
C     OO      = PHI( EIG(I,K),EPSI )
      OO      = MAX( ABS(EIG(I,K)),EPSI )
      O2      = EIG(I,K)*EIG(I,K)
      BET(I,K)= (OO-O2)*( ALFA(I,K) - GB(I,K) )
   25 CONTINUE
C
C     CALCULATE NUMERICAL FLUX
C
      DO 30 I = IB,NP
      FN(I,1) = TK(I,1,1)*BET(I,1) + TK(I,1,2)*BET(I,2)
     &        + TK(I,1,3)*BET(I,3)
      FN(I,2) = TK(I,2,1)*BET(I,1) + TK(I,2,2)*BET(I,2)
     &        + TK(I,2,3)*BET(I,3)
      FN(I,3) = TK(I,3,1)*BET(I,1) + TK(I,3,2)*BET(I,2)
     &        + TK(I,3,3)*BET(I,3)
   30 CONTINUE
C
C     UPDATE
C
      DO 35 I = IB+1,NP
      S(I,1)  = 0.5*( FN(I,1)-FN(I-1,1) )
      S(I,2)  = 0.5*( FN(I,2)-FN(I-1,2) )
      S(I,3)  = 0.5*( FN(I,3)-FN(I-1,3) )
   35 CONTINUE
      RETURN
      END
C
      SUBROUTINE HANCOCK(JB,JE,CK,OMIGA,Q,IMTH,IMOD)
      PARAMETER ( MAXI=407,NV=3 )
      DIMENSION FN(MAXI,NV),BET(MAXI,NV),EIG(MAXI,NV),TT(MAXI,NV,NV)
      DIMENSION UA(MAXI),CA(MAXI),Q(MAXI,NV)
      DIMENSION QR(MAXI,NV),QL(MAXI,NV),FR(MAXI,NV),FL(MAXI,NV)
      DIMENSION DQ(MAXI,NV),GB(MAXI,NV),GC(MAXI,NV)
      COMMON /D1/ GAMMA,GM1,EPSI,LEVEL,DTX,IPROB
C
      CKM   = 1. - CK
      CKP   = 1. + CK
      JUP   = JE - 1
      JLOW  = JB + 1
C
      DO 10 N = 1,NV
      QL(JB,N)= Q(JB,N)
      QR(JB,N)= Q(JB,N)
      QL(JE,N)= Q(JE,N)
      QR(JE,N)= Q(JE,N)
      QL(JLOW,N) = Q(JLOW,N)
      QR(JLOW,N) = Q(JLOW+1,N)
      QL(JUP ,N) = Q(JUP,N)
      QR(JUP ,N) = Q(JE,N)
   10 CONTINUE
      DO 20 J = JB,JUP
      DQ(J,1) = Q(J+1,1) - Q(J,1)
      DQ(J,2) = Q(J+1,2) - Q(J,2)
      DQ(J,3) = Q(J+1,3) - Q(J,3)
   20 CONTINUE
C
      IF( IMTH .EQ. 31 ) THEN
          DO 25 N = 1,NV
          DO 25 J = JLOW,JUP
          Z1      = SGN( DQ(J-1,N) )/2.
          Z2      = SGN( DQ( J ,N) )/2.
          Z3      = (Z1+Z2)*AMIN1( ABS(DQ(J-1,N)),ABS(OMIGA*DQ(J,N)) )
          Z4      = (Z1+Z2)*AMIN1( ABS(OMIGA*DQ(J-1,N)),ABS(DQ(J,N)) )
          QL(J  ,N) = Q(J,N) + .25*( CKM*Z3 + CKP*Z4 )
          QR(J-1,N) = Q(J,N) - .25*( CKP*Z3 + CKM*Z4 )
   25     CONTINUE
      END IF
C
      IF( IMTH .EQ. 32 .OR. IMTH .EQ. 33 ) THEN
C
C     CALCULATE MINMOD FOR GBAR
C
          DO 36 N= 1,NV
          DO 35 J= JLOW+1,JUP-1
          P1   = DQ(J+1,N) - DQ( J ,N)
          P2   = DQ( J ,N) - DQ(J-1,N)
          P4   = DQ(J-1,N) - DQ(J-2,N)
C         TMP1 = P2
C         TMP2 = P4
C         IF( ABS(P1) .LE. ABS(P2) ) TMP1= P1
C         IF( ABS(P2) .LE. ABS(P4) ) TMP2= P2
          TMP1 = 0.5*(SGN(P1)+SGN(P2))*AMIN1(ABS(P1),ABS(P2))
          TMP2 = 0.5*(SGN(P4)+SGN(P2))*AMIN1(ABS(P4),ABS(P2))
          Z1   = DQ( J ,N) - 0.5*TMP1
          Z2   = DQ(J-1,N) + 0.5*TMP2
          GB(J,N)= 0.5*( SGN(Z1)+SGN(Z2) )*AMIN1( ABS(Z1),ABS(Z2) )
   35     CONTINUE
C
C     B.C. for GBAR (zero order extrapolation)
C
          Z1        = DQ(JUP  ,N)
          Z2        = DQ(JUP-1,N)
          GB(JUP,N) = 0.5*( SGN(Z1)+SGN(Z2) )*AMIN1( ABS(Z1),ABS(Z2) )
          Z1        = DQ(JLOW,N)
          Z2        = DQ(JB  ,N)
          GB(JLOW,N)= 0.5*( SGN(Z1)+SGN(Z2) )*AMIN1( ABS(Z1),ABS(Z2) )
          GB(JB,N)  = GB(JLOW,N)
          GB(JE,N)  = GB(JUP ,N)
   36     CONTINUE
          IF( IMTH .EQ. 33 ) THEN
              DO 40 N = 1,NV
              DO 40 J = JLOW,JUP
              QR(J-1,N) = Q(J,N) - 0.5*GB(J,N)
              QL(J  ,N) = Q(J,N) + 0.5*GB(J,N)
   40         CONTINUE
          ELSE
              DO 45 N = 1,NV
              DO 45 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              IF( Z1*Z2 .GT. 0. ) THEN
                  BP = AMAX1( 2.,GB(J,N)/Z1 )
                  BM = AMAX1( 2.,GB(J,N)/Z2 )
                  CP = (2.*BP-CKP)/CKM
                  CM = (2.*BM-CKP)/CKM
                  Z3 = SGN(Z2)*AMIN1( ABS(Z2),CP*ABS(Z1) )
                  Z4 = SGN(Z1)*AMIN1( ABS(Z1),CM*ABS(Z2) )
              ELSE
                  Z3 = SGN( ABS(GB(J,N)) )*Z2
                  Z4 = SGN( ABS(GB(J,N)) )*Z1
              END IF
              QL(J  ,N)= Q(J,N) + .25*( CKM*Z3 + CKP*Z4 )
              QR(J-1,N)= Q(J,N) - .25*( CKP*Z3 + CKM*Z4 )
   45         CONTINUE
          END IF
      END IF
C
      IF( IMTH .EQ. 35 ) THEN
          DO 56 N = 1,NV
          DO 55 J = JLOW,JUP
          Z1      = DQ( J ,N)
          Z2      = DQ(J-1,N)
          GB(J,N) = 0.5*( SGN(Z1)+SGN(Z2) )*AMIN1( ABS(Z1),ABS(Z2) )
   55     CONTINUE
          GB(JB,N)= GB(JLOW,N)
          GB(JE,N)= GB(JUP ,N)
   56     CONTINUE
          DO 61 N = 1,NV
          DO 60 J = JLOW+1,JUP-1
          Z1      = DQ( J ,N)
          Z2      = DQ(J-1,N)
          P2      = Z1 - Z2
          IF( ABS(Z2) .GT. ABS(Z1) ) THEN
              P1    = DQ(J+1,N) - Z1
              TMP1  = P2
              IF( ABS(P1) .LE. ABS(P2) ) TMP1= P1
              GC(J,N)= -TMP1
          ELSE
              P4    = Z2 - DQ(J-2,N)
              TMP2  = P4
              IF( ABS(P2) .LE. ABS(P4) ) TMP2= P2
              GC(J,N)= 2.*TMP2
          END IF
   60     CONTINUE
          GC(JLOW,N)= GC(JLOW+1,N)
          GC(JB  ,N)= GC(JLOW  ,N)
          GC(JUP ,N)= GC(JUP-1 ,N)
          GC(JE  ,N)= GC(JUP   ,N)
   61     CONTINUE
          DO 65 N = 1,NV
          DO 65 J = JLOW,JUP
          QR(J-1,N) = Q(J,N) - 0.5*GB(J,N) - GC(J,N)/6.
          QL(J  ,N) = Q(J,N) + 0.5*GB(J,N) + GC(J,N)/6.
   65     CONTINUE
      END IF
C
      IF( IMTH .EQ. 34 ) THEN
          DO 76 N = 1,NV
          IF( IMOD .EQ. 0 ) THEN
              DO 70 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              GB(J,N)= Z1
              IF( ABS(Z2) .LE. ABS(Z1) ) GB(J,N)= Z2
   70         CONTINUE
          ELSEIF( IMOD .EQ. 1 ) THEN
              DO 71 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              GB(J,N)= 0.5*( SGN(Z1)+SGN(Z2) )*AMIN1( ABS(Z1),ABS(Z2) )
   71         CONTINUE
          ELSEIF( IMOD .EQ. 2 ) THEN
              DO 72 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              GB(J,N)= (Z1*Z2+0.0005)*(Z1+Z2)/( Z1*Z1+Z2*Z2+0.001 )
   72         CONTINUE
          ELSEIF( IMOD .EQ. 3 ) THEN
              DO 73 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              Z3= ( SGN(Z1)+SGN(Z2) )*AMIN1( ABS(Z1),ABS(Z2) )
              Z4= 0.5*( Z1+Z2 )
              GB(J,N)= 0.5*( SGN(Z3)+SGN(Z4) )*AMIN1( ABS(Z3),ABS(Z4) )
   73         CONTINUE
          ELSEIF( IMOD .EQ. 4 ) THEN
              DO 74 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              GB(J,N)= 0.
              IF( Z1+Z2 .NE. 0. ) GB(J,N)= (ABS(Z1*Z2)+Z1*Z2)/(Z1+Z2)
   74         CONTINUE
          ELSE
              DO 75 J = JLOW,JUP
              Z1      = DQ( J ,N)
              Z2      = DQ(J-1,N)
              Z3= 0.5*( SGN(Z1)+SGN(Z2) )*AMIN1( 2.*ABS(Z1),ABS(Z2) )
              Z4= 0.5*( SGN(Z1)+SGN(Z2) )*AMIN1( 2.*ABS(Z2),ABS(Z1) )
              GB(J,N)= SGN(Z1)*AMAX1( ABS(Z3),ABS(Z4) )
   75         CONTINUE
          END IF
          GB(JB,N)= GB(JLOW,N)
          GB(JE,N)= GB(JUP ,N)
   76     CONTINUE
          DO 77 N = 1,NV
          DO 77 J = JLOW,JUP
          QR(J-1,N) = Q(J,N) - 0.5*GB(J,N)
          QL(J  ,N) = Q(J,N) + 0.5*GB(J,N)
   77     CONTINUE
      END IF
C
      DO 80 J = JB,JUP
      UL = QL(J,2)/QL(J,1)
      PL = GM1*( QL(J,3) - 0.5*QL(J,1)*UL*UL )
      FL(J,1) = QL(J,1)*UL
      FL(J,2) = QL(J,2)*UL + PL
      FL(J,3) = ( QL(J,3) + PL )*UL
C
      UR = QR(J,2)/QR(J,1)
      PR = GM1*( QR(J,3) - 0.5*QR(J,1)*UR*UR )
      FR(J,1) = QR(J,1)*UR
      FR(J,2) = QR(J,2)*UR + PR
      FR(J,3) = ( QR(J,3) + PR )*UR
   80 CONTINUE
C
      DO 82 J = JLOW,JUP     ! OK 09/22/1997
      TMP1 = 0.5*DTX*( FL(J,1) - FR(J-1,1) )
      TMP2 = 0.5*DTX*( FL(J,2) - FR(J-1,2) )
      TMP3 = 0.5*DTX*( FL(J,3) - FR(J-1,3) )
      QL(J  ,1) = QL(J  ,1) - TMP1
      QL(J  ,2) = QL(J  ,2) - TMP2
      QL(J  ,3) = QL(J  ,3) - TMP3
      QR(J-1,1) = QR(J-1,1) - TMP1
      QR(J-1,2) = QR(J-1,2) - TMP2
      QR(J-1,3) = QR(J-1,3) - TMP3
   82 CONTINUE
C
      DO 85 J = JB,JUP
      RBAR    = SQRT( QR(J,1)/QL(J,1) )
      ORBAR   = 1./(1. + RBAR)
      UR      = QR(J,2)/QR(J,1)
      UL      = QL(J,2)/QL(J,1)
      UA(J)   = ( UR*RBAR + UL )*ORBAR
      HJR     = GAMMA*QR(J,3)/QR(J,1) - 0.5*GM1*UR*UR
      HJL     = GAMMA*QL(J,3)/QL(J,1) - 0.5*GM1*UL*UL
      HB      = ( HJR*RBAR + HJL )*ORBAR
      CSQ     = ( HB - 0.5*UA(J)*UA(J) )*GM1
      CA(J)   = SQRT(CSQ)
   85 CONTINUE
C
      DO 90 J = JB,JUP
      EIG(J,1)= PHI( UA(J),EPSI )
      EIG(J,2)= PHI( UA(J)+CA(J),EPSI )
      EIG(J,3)= PHI( UA(J)-CA(J),EPSI )
   90 CONTINUE
C
      CALL RDIGL(JB,JUP,UA,CA,EIG,TT,MAXI,GAMMA)
C
      DO 95 J = JB,JUP
      Q1      = QR(J,1) - QL(J,1)
      Q2      = QR(J,2) - QL(J,2)
      Q3      = QR(J,3) - QL(J,3)
      BET(J,1)= TT(J,1,1)*Q1+TT(J,1,2)*Q2+TT(J,1,3)*Q3
      BET(J,2)= TT(J,2,1)*Q1+TT(J,2,2)*Q2+TT(J,2,3)*Q3
      BET(J,3)= TT(J,3,1)*Q1+TT(J,3,2)*Q2+TT(J,3,3)*Q3
   95 CONTINUE
C
C     CALCULATE NUMERICAL FLUX
C
      DO 100 N = 1,NV
      DO 100 J = JB,JUP
      FN(J,N)  = 0.5*( FL(J,N) + FR(J,N) - BET(J,N) )
  100 CONTINUE
      DO 110 N = 1,NV
      DO 110 J = JLOW,JUP
      Q(J,N)   = Q(J,N) - DTX*( FN(J,N) - FN(J-1,N) )
  110 CONTINUE
      RETURN
      END
C
      SUBROUTINE PRIMER( IB,IE,Q )
      PARAMETER ( MAXI=407,NV=3 )
      COMMON /A1/ RQ(MAXI),UQ(MAXI),PQ(MAXI),EQ(MAXI),CQ(MAXI),HQ(MAXI)
      COMMON /D1/ GAMMA,GM1,EPSI,LEVEL,DTX,IPROB
      DIMENSION Q(MAXI,NV)
C
C     TOTALLY REFLECTIVE WALL BC FOR BLAST WAVE INTERACTION PROBLEM
C
      NP = IE - 3
C
      DO I = 1, 3
      Q( 4-I,1) = Q( 4+I,1)
      Q( 4-I,3) = Q( 4+I,3)
      Q(NP+I,1) = Q(NP-I,1)
      Q(NP+I,3) = Q(NP-I,3)
      IF( IPROB .GT. 3 ) THEN
          Q( 4-I,2) =-Q( 4+I,2)
          Q(NP+I,2) =-Q(NP-I,2)
      ELSE
          Q( 4-I,2) = Q( 4+I,2)
          Q(NP+I,2) = Q(NP-I,2)
      END IF
      END DO
C
C     CALCULATE RHO,U,P,E,H,C
C
      DO 25 I = IB,IE
      RQ(I) = ABS( Q(I,1) )
      UQ(I) = Q(I,2)/RQ(I)
      PQ(I) = GM1*ABS( Q(I,3) - 0.5*RQ(I)*UQ(I)*UQ(I) )
      EQ(I) = PQ(I)/GM1/RQ(I)
      HQ(I) = GAMMA*EQ(I) + 0.5*UQ(I)**2
      CQ(I) = SQRT( GAMMA*PQ(I)/RQ(I) )
   25 CONTINUE
      RETURN
      END
