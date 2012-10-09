      PROGRAM A_RELA2D
C
C     BEAM TYPE SCHEME SOLVE EULER LIMIT BOLTZMANN EQUATION
C
C       DEVELOPED BY  YI-NAN TSAI & CJW
C       VERSION X  START AT 1996  9 20
C                  FOR ARBITRARY BOUNDARY SET
C
C       MAXWK : MAX WORKING SPACE
C       MAXIJP: MAX POINT IN I OR J DIRECTION
C       MAXGRD: MAX GRID POINT
C       MAXLIN: MAX LU SWEEP LINE (DIAG. LINE)
C       MAXBCP: MAX B.C. IDEX SPACE
C       MAXLUI: MAX LU IDEX SPACE
C
      PARAMETER (MAXIP=371,MAXJP=261,MAXIJP=371)
c     PARAMETER (MAXIP=517,MAXJP=517,MAXIJP=517)
      PARAMETER (MAXWK= 66*MAXIP*MAXJP+2)
      PARAMETER (MAXGRD=MAXIP*MAXJP,MAXLIN=MAXIP+MAXJP)
      PARAMETER (MAXBCP=MAXIJP*10*2*4,MAXLUI=MAXGRD*2+MAXLIN*2+1)
C
      COMMON /TAPE / INAME,IGRID,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /NEXT / NEXTLOC
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /PNT  / LX,LY,LXIX,LXIY,LETAX,LETAY,LJAC,
     1               LR,LE,LT,LU,LV,LP,LXM,LPSI,LRNP,LRUNP,LRVNP,LRENP
     1              ,LUS,LVS,LRS,LRUS,LRVS,LRES,LDRS,LDRUS,LDRVS,LDRES
     1               ,LUSP,LUSM,LVSP,LVSM,LRNPO
      COMMON /PNTI / LKDI,LKDJ,LLKD,LKDP,MDA,MDB,NDIG
      DIMENSION S(MAXWK),MLU(MAXLUI),IJP(MAXBCP)
C
      MAXW  = MAXWK
      MAXI  = MAXIJP
      MAXG  = MAXGRD
      MAXL  = MAXLIN
      MAXBP = MAXBCP
      MAXLI = MAXLUI
      NEXTLOC = 1
      IGRID  = 11
      INAME  = 12
      IFORC  = 15
      IRSTA  = 16
      IRESD  = 8
      ISHOW  = 6
C
      READ(IGRID,*)NI,NJ
      CALL SETUP(IJP(1))
      CALL POINT
C
C     CALL LUPNT(MLU(LKDI),MLU(LKDJ),MLU(LLKD),MLU(LKDP))
C
      CALL SOLV(IJP,MLU,S(LX),S(LY),S(LXIX),S(LXIY),S(LETAX),S(LETAY)
     1  ,S(LJAC),S(LR),S(LE),S(LT),S(LU),S(LV),S(LP),S(LXM),S(LPSI)
     1  ,S(LUS),S(LVS),S(LRS),S(LRUS),S(LRVS),S(LRES),S(LDRS),S(LRNP)
     1  ,S(LRUNP),S(LRVNP),S(LRENP),S(LRNPO))
C
      STOP
      END
      SUBROUTINE BC(IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P,X)
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /CNST2/ R0,U0,V0,P0,T0,E0
      COMMON /CNST4/ DR0,UR0,VR0,PR0,TR0,ER0,
     1               DL0,UL0,VL0,PL0,TL0,EL0
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /BC01 / IDBC(30,20),FDBC(10,20)
      DIMENSION IJP(MXBP)
      DIMENSION  XIX(-1:NI2,-1:NJ2),   XIY(-1:NI2,-1:NJ2)
      DIMENSION ETAX(-1:NI2,-1:NJ2),  ETAY(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2),     E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2),     V(-1:NI2,-1:NJ2)
      DIMENSION    T(-1:NI2,-1:NJ2),     P(-1:NI2,-1:NJ2)
      DIMENSION    X(-1:NI2,-1:NJ2)
C..
C     IBTYPE = 0     ZERO GRADIENT
C            = 1     WALL
C            = 2     SYMMETRY I WITH Y=0 ( Cylinder or airfoil )
C            = 3     SYMMETRY J WITH Y=0
C            = 4     WAKE (FOR C GRID WAKE MATCH)
C            = 5     FARFIELD
      DO 100 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IF(IBTYPE .EQ. 0)
     1   CALL BCZERO(IBC,IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P)
      IF(IBTYPE .EQ. 1)
     1   CALL BCWALL(IBC,IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P)
      IF(IBTYPE .EQ. 2 .OR. IBTYPE .EQ. 3)
     1   CALL BCSYMY(IBC,IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P)
      IF(IBTYPE .EQ. 4 .OR. IBTYPE .EQ. 7)
     1   CALL BCWAKE(IBC,IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P)
      IF(IBTYPE .EQ. 5)
     1   CALL BCFARD(IBC,IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P,X)
  100 CONTINUE
      RETURN
      END
      SUBROUTINE BCFARD(IBC,IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P,X)
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /CNST2/ R0,U0,V0,P0,T0,E0
      COMMON /CNST4/ DR0,UR0,VR0,PR0,TR0,ER0,
     1               DL0,UL0,VL0,PL0,TL0,EL0
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /BC01 / IDBC(30,20),FDBC(10,20)
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
C
      DIMENSION IJP(MXBP)
      DIMENSION  XIX(-1:NI2,-1:NJ2),   XIY(-1:NI2,-1:NJ2)
      DIMENSION ETAX(-1:NI2,-1:NJ2),  ETAY(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2),     E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2),     V(-1:NI2,-1:NJ2)
      DIMENSION    T(-1:NI2,-1:NJ2),     P(-1:NI2,-1:NJ2)
      DIMENSION    X(-1:NI2,-1:NJ2)
C
      IDA        = IDBC( 2,IBC)
      ISA        = IDBC( 3,IBC)
      ICA        = IDBC( 4,IBC)
      IAS        = IDBC( 5,IBC)
      IAE        = IDBC( 6,IBC)
      IPA        = IDBC( 7,IBC)
      LENA       = (IAE - IAS)*ICA + 1
      ISGN       = -1
      IF(ISA .EQ. 1) ISGN = 1
      DO 100 II  = 1,LENA
        IA       = IJP(IPA + LENA*0 + II)
        JA       = IJP(IPA + LENA*1 + II)
        IB       = IJP(IPA + LENA*2 + II)
        JB       = IJP(IPA + LENA*3 + II)
        IC       = IJP(IPA + LENA*4 + II)
        JC       = IJP(IPA + LENA*5 + II)
        ID       = IJP(IPA + LENA*6 + II)
        JD       = IJP(IPA + LENA*7 + II)
        IE       = IJP(IPA + LENA*8 + II)
        JE       = IJP(IPA + LENA*9 + II)
      IF(IDA .EQ. 1) THEN
          E00    = SQRT(XIX(IC,JC)*XIX(IC,JC)+XIY(IC,JC)*XIY(IC,JC))
          E1     = XIX(IC,JC)/(E00 + 1.0E-25)
          E2     = XIY(IC,JC)/(E00 + 1.0E-25)
      ELSE
          E00    = SQRT(ETAX(IC,JC)*ETAX(IC,JC)+ETAY(IC,JC)*ETAY(IC,JC))
          E1     = ETAX(IC,JC)/(E00+ 1.0E-25)
          E2     = ETAY(IC,JC)/(E00+ 1.0E-25)
      ENDIF
      RDF        = ABS( R(ID,JD) )
      PDF        = ABS( P(ID,JD) )
      UDF        =      U(ID,JD)
      VDF        =      V(ID,JD)
      UNDF       = E1*UDF + E2*VDF
C
      RCF        = ABS( R(IC,JC) )
      PCF        = ABS( P(IC,JC) )
      UCF        =      U(IC,JC)
      VCF        =      V(IC,JC)
      UNCF       = E1*UCF + E2*VCF
      UTOT       = SQRT( UCF*UCF + VCF*VCF)
C
       if( lstd .eq. 0 )then
        P(IC,JC)  = PDF
        U(IC,JC)  = UDF
        V(IC,JC)  = VDF
        R(IC,JC)  = RDF
        E(IC,JC)  = P(IC,JC)/(GAMMA-1.) + R(IC,JC)
       endif
C
       P(IB,JB)  =  P(IC,JC)
       T(IB,JB)  =  T(IC,JC)
       R(IB,JB)  =  R(IC,JC)
       U(IB,JB)  =  U(IC,JC)
       V(IB,JB)  =  V(IC,JC)
       E(IB,JB)  =  E(IC,JC)
C
       P(IA,JA)  =  P(IC,JC)
       T(IA,JA)  =  T(IC,JC)
       R(IA,JA)  =  R(IC,JC)
       U(IA,JA)  =  U(IC,JC)
       V(IA,JA)  =  V(IC,JC)
       E(IA,JA)  =  E(IC,JC)
C
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE BCSYMY(IBC,IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P)
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /CNST2/ R0,U0,V0,P0,T0,E0
      COMMON /CNST4/ DR0,UR0,VR0,PR0,TR0,ER0,
     1               DL0,UL0,VL0,PL0,TL0,EL0
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /BC01 / IDBC(30,20),FDBC(10,20)
      DIMENSION IJP(MXBP)
      DIMENSION  XIX(-1:NI2,-1:NJ2),   XIY(-1:NI2,-1:NJ2)
      DIMENSION ETAX(-1:NI2,-1:NJ2),  ETAY(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2),     E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2),     V(-1:NI2,-1:NJ2)
      DIMENSION    T(-1:NI2,-1:NJ2),     P(-1:NI2,-1:NJ2)
C
      IDA        = IDBC( 2,IBC)
      ISA        = IDBC( 3,IBC)
      ICA        = IDBC( 4,IBC)
      IAS        = IDBC( 5,IBC)
      IAE        = IDBC( 6,IBC)
      IPA        = IDBC( 7,IBC)
      LENA       = (IAE - IAS)*ICA + 1
      ISGN       = -1
      IF(ISA .EQ. 1) ISGN = 1
      DO 120 II  = 1,LENA
        IA       = IJP(IPA + LENA*0 + II)
        JA       = IJP(IPA + LENA*1 + II)
        IB       = IJP(IPA + LENA*2 + II)
        JB       = IJP(IPA + LENA*3 + II)
        IC       = IJP(IPA + LENA*4 + II)
        JC       = IJP(IPA + LENA*5 + II)
        ID       = IJP(IPA + LENA*6 + II)
        JD       = IJP(IPA + LENA*7 + II)
        IE       = IJP(IPA + LENA*8 + II)
        JE       = IJP(IPA + LENA*9 + II)
       P(IB,JB)  =  P(ID,JD)
       R(IB,JB)  =  R(ID,JD)
       U(IB,JB)  =  U(ID,JD)
       V(IB,JB)  = -V(ID,JD)
       E(IB,JB)  =  E(ID,JD)
       P(IA,JA)  =  P(IE,JE)
       R(IA,JA)  =  R(IE,JE)
       U(IA,JA)  =  U(IE,JE)
       V(IA,JA)  = -V(IE,JE)
       E(IA,JA)  =  E(IE,JE)
       R(IC,JC) =0.5*(R(IB,JB)         + R(ID,JD)         )
       E(IC,JC) =0.5*(E(IB,JB)         + E(ID,JD)         )
       U(IC,JC) =0.5*(U(IB,JB)         + U(ID,JD)         )
       V(IC,JC)  = 0.
       P(IC,JC)  = GM1*( E(IC,JC)-R(IC,JC) )
c
c      RU2      =0.5*(R(IB,JB)*U(IB,JB)+ R(ID,JD)*U(ID,JD))
c      RE2      =0.5*(R(IB,JB)*E(IB,JB)+ R(ID,JD)*E(ID,JD))
c      E(IC,JC)  = RE2/R(IC,JC)
c      U(IC,JC)  = RU2/R(IC,JC)
c    1                    -(U(IC,JC)*U(IC,JC)+V(IC,JC)*V(IC,JC)))
  120 CONTINUE
      RETURN
      END
      SUBROUTINE BCWAKE(IBC,IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P)
C
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /CNST2/ R0,U0,V0,P0,T0,E0
      COMMON /CNST4/ DR0,UR0,VR0,PR0,TR0,ER0,
     1               DL0,UL0,VL0,PL0,TL0,EL0
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /BC01 / IDBC(30,20),FDBC(10,20)
      DIMENSION IJP(MXBP)
      DIMENSION  XIX(-1:NI2,-1:NJ2),   XIY(-1:NI2,-1:NJ2)
      DIMENSION ETAX(-1:NI2,-1:NJ2),  ETAY(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2),     E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2),     V(-1:NI2,-1:NJ2)
      DIMENSION    T(-1:NI2,-1:NJ2),     P(-1:NI2,-1:NJ2)
C
      IDA        = IDBC( 2,IBC)
      ISA        = IDBC( 3,IBC)
      ICA        = IDBC( 4,IBC)
      IAS        = IDBC( 5,IBC)
      IAE        = IDBC( 6,IBC)
      IPA        = IDBC( 7,IBC)
      LENA       = (IAE - IAS)*ICA + 1
      IDB        = IDBC( 12,IBC)
      ISB        = IDBC( 13,IBC)
      ICB        = IDBC( 14,IBC)
      IBS        = IDBC( 15,IBC)
      IBE        = IDBC( 16,IBC)
      IPB        = IDBC(  8,IBC)
      LENB       = (IBE - IBS)*ICB + 1
      ISGN       = -1
      IF(ISA .EQ. 1) ISGN = 1
      DO 100 II  = 1,LENA
        IA1      = IJP(IPA + LENA*0 + II)
        JA1      = IJP(IPA + LENA*1 + II)
        IB1      = IJP(IPA + LENA*2 + II)
        JB1      = IJP(IPA + LENA*3 + II)
        IC1      = IJP(IPA + LENA*4 + II)
        JC1      = IJP(IPA + LENA*5 + II)
        ID1      = IJP(IPA + LENA*6 + II)
        JD1      = IJP(IPA + LENA*7 + II)
        IE1      = IJP(IPA + LENA*8 + II)
        JE1      = IJP(IPA + LENA*9 + II)
C
        IA2      = IJP(IPB + LENA*0 + II)
        JA2      = IJP(IPB + LENA*1 + II)
        IB2      = IJP(IPB + LENA*2 + II)
        JB2      = IJP(IPB + LENA*3 + II)
        IC2      = IJP(IPB + LENA*4 + II)
        JC2      = IJP(IPB + LENA*5 + II)
        ID2      = IJP(IPB + LENA*6 + II)
        JD2      = IJP(IPB + LENA*7 + II)
        IE2      = IJP(IPB + LENA*8 + II)
        JE2      = IJP(IPB + LENA*9 + II)
C
       P(IC1,JC1)  = 0.5*( P(ID1,JD1) + P(ID2,JD2) )
       U(IC1,JC1)  = 0.5*( U(ID1,JD1) + U(ID2,JD2) )
       V(IC1,JC1)  = 0.5*( V(ID1,JD1) + V(ID2,JD2) )
       R(IC1,JC1)  = 0.5*( R(ID1,JD1) + R(ID2,JD2) )
c      T(IC1,JC1)  = 0.5*( T(ID1,JD1) + T(ID2,JD2) )
       E(IC1,JC1)  = 0.5*( E(ID1,JD1) + E(ID2,JD2) )
C
       P(IC2,JC2)  = P(IC1,JC1)
c      T(IC2,JC2)  = T(IC1,JC1)
       R(IC2,JC2)  = R(IC1,JC1)
       U(IC2,JC2)  = U(IC1,JC1)
       V(IC2,JC2)  = V(IC1,JC1)
       E(IC2,JC2)  = E(IC1,JC1)
C
       P(IA1,JA1)  = P(IE2,JE2)
       P(IB1,JB1)  = P(ID2,JD2)
       P(IA2,JA2)  = P(IE1,JE1)
       P(IB2,JB2)  = P(ID1,JD1)
C
c      T(IA1,JA1)  = T(IE2,JE2)
c      T(IB1,JB1)  = T(ID2,JD2)
c      T(IA2,JA2)  = T(IE1,JE1)
c      T(IB2,JB2)  = T(ID1,JD1)
C
       R(IA1,JA1)  = R(IE2,JE2)
       R(IB1,JB1)  = R(ID2,JD2)
       R(IA2,JA2)  = R(IE1,JE1)
       R(IB2,JB2)  = R(ID1,JD1)
C
       U(IA1,JA1)  = U(IE2,JE2)
       U(IB1,JB1)  = U(ID2,JD2)
       U(IA2,JA2)  = U(IE1,JE1)
       U(IB2,JB2)  = U(ID1,JD1)
C
       V(IA1,JA1)  = V(IE2,JE2)
       V(IB1,JB1)  = V(ID2,JD2)
       V(IA2,JA2)  = V(IE1,JE1)
       V(IB2,JB2)  = V(ID1,JD1)
C
       E(IA1,JA1)  = E(IE2,JE2)
       E(IB1,JB1)  = E(ID2,JD2)
       E(IA2,JA2)  = E(IE1,JE1)
       E(IB2,JB2)  = E(ID1,JD1)
100    CONTINUE
C
      RETURN
      END
      SUBROUTINE BCWALL(IBC,IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P)
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /CNST2/ R0,U0,V0,P0,T0,E0
      COMMON /CNST4/ DR0,UR0,VR0,PR0,TR0,ER0,
     1               DL0,UL0,VL0,PL0,TL0,EL0
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /BC01 / IDBC(30,20),FDBC(10,20)
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
C
      DIMENSION IJP(MXBP)
      DIMENSION  XIX(-1:NI2,-1:NJ2),   XIY(-1:NI2,-1:NJ2)
      DIMENSION ETAX(-1:NI2,-1:NJ2),  ETAY(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2),     E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2),     V(-1:NI2,-1:NJ2)
      DIMENSION    T(-1:NI2,-1:NJ2),     P(-1:NI2,-1:NJ2)
C
      IDA        = IDBC( 2,IBC)
      ISA        = IDBC( 3,IBC)
      ICA        = IDBC( 4,IBC)
      IAS        = IDBC( 5,IBC)
      IAE        = IDBC( 6,IBC)
      IPA        = IDBC( 7,IBC)
      LENA       = (IAE - IAS)*ICA + 1
      ISGN       = -1
      IF(ISA .EQ. 1) ISGN = 1
      DO 100 II  = 1,LENA
        IA       = IJP(IPA + LENA*0 + II)
        JA       = IJP(IPA + LENA*1 + II)
        IB       = IJP(IPA + LENA*2 + II)
        JB       = IJP(IPA + LENA*3 + II)
        IC       = IJP(IPA + LENA*4 + II)
        JC       = IJP(IPA + LENA*5 + II)
        ID       = IJP(IPA + LENA*6 + II)
        JD       = IJP(IPA + LENA*7 + II)
        IE       = IJP(IPA + LENA*8 + II)
        JE       = IJP(IPA + LENA*9 + II)
      IF(IDA .EQ. 1) THEN
          E00    = SQRT(XIX(IC,JC)*XIX(IC,JC)+XIY(IC,JC)*XIY(IC,JC))
          E1     = XIX(IC,JC)/(E00 + 1.0E-25)
          E2     = XIY(IC,JC)/(E00 + 1.0E-25)
      ELSE
          E00    = SQRT(ETAX(IC,JC)*ETAX(IC,JC)+ETAY(IC,JC)*ETAY(IC,JC))
          E1     = ETAX(IC,JC)/(E00+ 1.0E-25)
          E2     = ETAY(IC,JC)/(E00+ 1.0E-25)
      ENDIF
       RCF        =  ABS(R(ID,JD))
       PCF        =  ABS(P(ID,JD))
       ECF        =  E(ID,JD)
       UCF        =  U(ID,JD)
       VCF        =  V(ID,JD)
       UNCF       = E1*UCF + E2*VCF
       UTCF       = E2*UCF - E1*VCF
C
C       XC         = PCF/RCF**GAMMA
C
C       GOTO 13
C
       XC         = PCF/RCF**GAMMA
       XK         = 1./GM1
       XM         = 1./(GAMMA*XC)
       XRHOS      = SQRT(XM*RCF**(-0.4)+XK)
       X0         = SQRT(XK)
       X1         = XRHOS - X0
       X2         = XRHOS + X0
       X3         = 1. + UNCF
       X4         = 1. - UNCF
       XPHIS      = -5./2./X0*LOG(ABS(X1/X2))
       XUS        = XPHIS - 0.5*LOG(ABS(X3/X4))
       XUES       = EXP( -2./5.*X0*XUS )
       XRHO       = X0*(1.+XUES)/(1.-XUES)
C
       UCB        = UTCF * E2
       VCB        =-UTCF * E1
C       IF( LSTD .EQ. 1 )THEN
C       ELSE
        RCB       = (XM/(XRHO*XRHO-XK))**2.5
        PCB       = XC*RCB**GAMMA
C       ENDIF
C
       P(IC,JC)  = PCB
       U(IC,JC)  = UCB
       V(IC,JC)  = VCB
       R(IC,JC)  = RCB
       E(IC,JC)  = PCB*XK+RCB
C
C  13   CONTINUE
C       U(IC,JC) =  UTCF * E2
C       V(IC,JC) = -UTCF * E1
C       R(IC,JC) = ABS(2.*R(ID,JD)-R(IE,JE)) 
C       P(IC,JC) = XC*R(IC,JC)**GAMMA
C       E(IC,JC) = P(IC,JC)/(GAMMA-1.)+R(IC,JC)
C
       U(IB,JB)  = U(IC,JC)
       V(IB,JB)  = V(IC,JC)
       T(IB,JB)  = T(IC,JC)
C       R(IB,JB)  = ABS(2.*R(IC,JC)-R(ID,JD))
       R(IB,JB)  = R(IC,JC)
C       P(IB,JB) = XC*R(IB,JB)**GAMMA
       P(IB,JB) = P(IC,JC)
       E(IB,JB) = P(IB,JB)/(GAMMA-1.)+R(IB,JB)
C
       U(IA,JA)  = U(IC,JC)
       V(IA,JA)  = V(IC,JC)
       T(IA,JA)  = T(IC,JC)
C       R(IA,JA)  = ABS(2.*R(IB,JB)-R(IC,JC))
       R(IA,JA)  = R(IC,JC)
C       P(IA,JA) = XC*R(IA,JA)**GAMMA
       P(IA,JA) = P(IC,JC)
       E(IA,JA) = P(IA,JA)/(GAMMA-1.)+R(IA,JA)
C
  100 CONTINUE
      RETURN
      END
      SUBROUTINE BCZERO(IBC,IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P)
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /CNST2/ R0,U0,V0,P0,T0,E0
      COMMON /CNST4/ DR0,UR0,VR0,PR0,TR0,ER0,
     1               DL0,UL0,VL0,PL0,TL0,EL0
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /BC01 / IDBC(30,20),FDBC(10,20)
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
C
      DIMENSION IJP(MXBP)
      DIMENSION  XIX(-1:NI2,-1:NJ2),   XIY(-1:NI2,-1:NJ2)
      DIMENSION ETAX(-1:NI2,-1:NJ2),  ETAY(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2),     E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2),     V(-1:NI2,-1:NJ2)
      DIMENSION    T(-1:NI2,-1:NJ2),     P(-1:NI2,-1:NJ2)
C
      IDA        = IDBC( 2,IBC)
      ISA        = IDBC( 3,IBC)
      ICA        = IDBC( 4,IBC)
      IAS        = IDBC( 5,IBC)
      IAE        = IDBC( 6,IBC)
      IPA        = IDBC( 7,IBC)
      LENA       = (IAE - IAS)*ICA + 1
      ISGN       = -1
      IF(ISA .EQ. 1) ISGN = 1
      DO 100 II  = 1,LENA
        IA       = IJP(IPA + LENA*0 + II)
        JA       = IJP(IPA + LENA*1 + II)
        IB       = IJP(IPA + LENA*2 + II)
        JB       = IJP(IPA + LENA*3 + II)
        IC       = IJP(IPA + LENA*4 + II)
        JC       = IJP(IPA + LENA*5 + II)
        ID       = IJP(IPA + LENA*6 + II)
        JD       = IJP(IPA + LENA*7 + II)
        IE       = IJP(IPA + LENA*8 + II)
        JE       = IJP(IPA + LENA*9 + II)
      IF(IDA .EQ. 1) THEN
          E00    = SQRT(XIX(IC,JC)*XIX(IC,JC)+XIY(IC,JC)*XIY(IC,JC))
          E1     = XIX(IC,JC)/(E00 + 1.0E-25)
          E2     = XIY(IC,JC)/(E00 + 1.0E-25)
      ELSE
          E00    = SQRT(ETAX(IC,JC)*ETAX(IC,JC)+ETAY(IC,JC)*ETAY(IC,JC))
          E1     = ETAX(IC,JC)/(E00+ 1.0E-25)
          E2     = ETAY(IC,JC)/(E00+ 1.0E-25)
      ENDIF
       RCF        = ABS( R(ID,JD) )
       PCF        = ABS( P(ID,JD) )
       UCF        = U(ID,JD)
       VCF        = V(ID,JD)
c
       UNCF       = E1*UCF + E2*VCF
       UTCF       = E2*UCF - E1*VCF
       CCF        = SQRT( GAMMA*PCF/RCF )
       UCB        = UTCF * E2
       VCB        =-UTCF * E1
C
       U(IC,JC)  = UCB
       V(IC,JC)  = VCB
       R(IC,JC)  = R(ID,JD)
       E(IC,JC)  = E(ID,JD)
       P(IC,JC)  = GM1*( E(IC,JC) - R(IC,JC) )
c
       P(IB,JB)  = P(IC,JC)
       U(IB,JB)  = U(IC,JC)
       V(IB,JB)  = V(IC,JC)
       T(IB,JB)  = T(IC,JC)
       R(IB,JB)  = R(IC,JC)
       E(IB,JB)  = E(IC,JC)
C
       P(IA,JA)  = P(IC,JC)
       U(IA,JA)  = U(IC,JC)
       V(IA,JA)  = V(IC,JC)
       T(IA,JA)  = T(IC,JC)
       R(IA,JA)  = R(IC,JC)
       E(IA,JA)  = E(IC,JC)
C
  100 CONTINUE
      RETURN
      END
      SUBROUTINE BEAM5(R,E,RS,RUS,RVS,RES,U,V,P,US,VS,
     1                 XIX,XIY,ETAX,ETAY,QJ)
      COMMON /TAPE / INAME,IGRID,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      DIMENSION  RUS(-1:NI2,-1:NJ2,5), RVS(-1:NI2,-1:NJ2,5)
      DIMENSION   RS(-1:NI2,-1:NJ2,5), RES(-1:NI2,-1:NJ2,5)
      DIMENSION   US(-1:NI2,-1:NJ2,5),  VS(-1:NI2,-1:NJ2,5)
      DIMENSION    U(-1:NI2,-1:NJ2),     V(-1:NI2,-1:NJ2)
      DIMENSION    P(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2),     E(-1:NI2,-1:NJ2)
      DIMENSION  XIX(-1:NI2,-1:NJ2),   XIY(-1:NI2,-1:NJ2)
      DIMENSION ETAX(-1:NI2,-1:NJ2),  ETAY(-1:NI2,-1:NJ2)
      DIMENSION   QJ(-1:NI2,-1:NJ2)
C
      XGAMAD = 1./(3.*GAMMA-4.)
      DO 10 IJ     = -1, MIJ2
      I         = MOD(IJ+1,MI)-1
      J         = (IJ+1)/MI -1
      IF(ABS(U(I,J)) .LE. 1.E-8 .AND. ABS(V(I,J)) .LE. 1.E-8)THEN
C       XGAMAD = (E(IJ,-1)-R(IJ,-1))/(3.*P(IJ,-1)-(E(IJ,-1)-R(IJ,-1)))
C       XGAMAD = 1./(3.*GAMMA-4.)
       XDV    = SQRT(XGAMAD**2.-1.)/XGAMAD
       XB     = ABS(E(IJ,-1)/R(IJ,-1)-1.)/6./XGAMAD/(XGAMAD-1.)
       XA     = 1.-6.*XB*XGAMAD
       X1     = R(IJ,-1)*(XA+2.*XB*XGAMAD)
       X2     = R(IJ,-1)*(XA+2.*XB*XGAMAD*XGAMAD)
       X3     = R(IJ,-1)*XB*XGAMAD
       X4     = R(IJ,-1)*XB*XGAMAD*XGAMAD*XDV
       X5     = R(IJ,-1)*XB*XGAMAD*XGAMAD
C
       US(IJ,-1,1) = 0.
       US(IJ,-1,2) = XDV*XIX(IJ,-1)
       US(IJ,-1,3) =-XDV*XIX(IJ,-1)
       US(IJ,-1,4) = XDV*XIY(IJ,-1)
       US(IJ,-1,5) =-XDV*XIY(IJ,-1)
C
       VS(IJ,-1,1) = 0.
       VS(IJ,-1,2) = XDV*ETAX(IJ,-1)
       VS(IJ,-1,3) =-XDV*ETAX(IJ,-1)
       VS(IJ,-1,4) = XDV*ETAY(IJ,-1)
       VS(IJ,-1,5) =-XDV*ETAY(IJ,-1)
C
       RS(IJ,-1,1) = X1/QJ(IJ,-1)
       RS(IJ,-1,2) = X3/QJ(IJ,-1)
       RS(IJ,-1,3) = X3/QJ(IJ,-1)
       RS(IJ,-1,4) = X3/QJ(IJ,-1)
       RS(IJ,-1,5) = X3/QJ(IJ,-1)
C
       RUS(IJ,-1,1) = 0.
       RUS(IJ,-1,2) = X4/QJ(IJ,-1)
       RUS(IJ,-1,3) =-X4/QJ(IJ,-1)
       RUS(IJ,-1,4) = 0.
       RUS(IJ,-1,5) = 0.
C
       RVS(IJ,-1,1) = 0.
       RVS(IJ,-1,2) = 0.
       RVS(IJ,-1,3) = 0.
       RVS(IJ,-1,4) = X4/QJ(IJ,-1)
       RVS(IJ,-1,5) =-X4/QJ(IJ,-1)
C
       RES(IJ,-1,1) = X2/QJ(IJ,-1)
       RES(IJ,-1,2) = X5/QJ(IJ,-1)
       RES(IJ,-1,3) = X5/QJ(IJ,-1)
       RES(IJ,-1,4) = X5/QJ(IJ,-1)
       RES(IJ,-1,5) = X5/QJ(IJ,-1)
      ELSE
C      XGAMAD = (E(IJ,-1)-R(IJ,-1))/(3.*P(IJ,-1)-(E(IJ,-1)-R(IJ,-1)))
C      XGAMAD = 1./(3.*GAMMA-4.)
      XDV    = SQRT(XGAMAD**2.-1.)/XGAMAD
      XU2    = U(IJ,-1)*U(IJ,-1) + V(IJ,-1)*V(IJ,-1)
      XGAMAU = 1./SQRT(ABS(1.-XU2))
      XB     = ABS(E(IJ,-1)/R(IJ,-1)-1.)/6./XGAMAD/(XGAMAD-1.)
      XA     = 1.-6.*XB*XGAMAD
      X1     = XGAMAU*XGAMAD*U(IJ,-1)
      X2     = (XGAMAU*U(IJ,-1)**2.+V(IJ,-1)**2.)*XGAMAD*XDV/XU2
      X3     = XGAMAU*XGAMAD*(1. + U(IJ,-1)*XDV)
      X4     = XGAMAU*XGAMAD*(1. - U(IJ,-1)*XDV)
      X5     = XGAMAU*XGAMAD*(1. + V(IJ,-1)*XDV)
      X6     = XGAMAU*XGAMAD*(1. - V(IJ,-1)*XDV)
      X7     = (XGAMAU-1.)*U(IJ,-1)*V(IJ,-1)*XGAMAD*XDV/XU2
      X8     = XGAMAU*XGAMAD*V(IJ,-1)
      X9     = (XGAMAU*V(IJ,-1)**2.+U(IJ,-1)**2.)*XGAMAD*XDV/XU2
      X11    = R(IJ,-1)*XGAMAU
      X12    = XA+2.*XB*XGAMAD*XGAMAD
      X13    = XB*R(IJ,-1)*XGAMAU*XGAMAD*(1. + U(IJ,-1)*XDV)
      X14    = XB*R(IJ,-1)*XGAMAU*XGAMAD*(1. - U(IJ,-1)*XDV)
      X15    = XB*R(IJ,-1)*XGAMAU*XGAMAD*(1. + V(IJ,-1)*XDV)
      X16    = XB*R(IJ,-1)*XGAMAU*XGAMAD*(1. - V(IJ,-1)*XDV)
      X17    = XB*R(IJ,-1)*XGAMAU*XGAMAU*XGAMAD*XGAMAD*
     1         (1.+U(IJ,-1)*XDV)*(1.+U(IJ,-1)*XDV)
      X18    = XB*R(IJ,-1)*XGAMAU*XGAMAU*XGAMAD*XGAMAD*
     1         (1.-U(IJ,-1)*XDV)*(1.-U(IJ,-1)*XDV)
      X19    = XB*R(IJ,-1)*XGAMAU*XGAMAU*XGAMAD*XGAMAD*
     1         (1.+V(IJ,-1)*XDV)*(1.+V(IJ,-1)*XDV)
      X20    = XB*R(IJ,-1)*XGAMAU*XGAMAU*XGAMAD*XGAMAD*
     1         (1.-V(IJ,-1)*XDV)*(1.-V(IJ,-1)*XDV)
C
       US(IJ,-1,1) = U(IJ,-1)*XIX(IJ,-1)+V(IJ,-1)*XIY(IJ,-1)
       US(IJ,-1,2) = (X1+X2)/X3*XIX(IJ,-1) + (X8+X7)/X3*XIY(IJ,-1)
       US(IJ,-1,3) = (X1-X2)/X4*XIX(IJ,-1) + (X8-X7)/X4*XIY(IJ,-1)
       US(IJ,-1,4) = (X1+X7)/X5*XIX(IJ,-1) + (X8+X9)/X5*XIY(IJ,-1)
       US(IJ,-1,5) = (X1-X7)/X6*XIX(IJ,-1) + (X8-X9)/X6*XIY(IJ,-1)
C
       VS(IJ,-1,1) = U(IJ,-1)*ETAX(IJ,-1)+V(IJ,-1)*ETAY(IJ,-1)
       VS(IJ,-1,2) = (X1+X2)/X3*ETAX(IJ,-1) + (X8+X7)/X3*ETAY(IJ,-1)
       VS(IJ,-1,3) = (X1-X2)/X4*ETAX(IJ,-1) + (X8-X7)/X4*ETAY(IJ,-1)
       VS(IJ,-1,4) = (X1+X7)/X5*ETAX(IJ,-1) + (X8+X9)/X5*ETAY(IJ,-1)
       VS(IJ,-1,5) = (X1-X7)/X6*ETAX(IJ,-1) + (X8-X9)/X6*ETAY(IJ,-1)
C
       RS(IJ,-1,1) = X11*(XA+2.*XB*XGAMAD)/QJ(IJ,-1)
       RS(IJ,-1,2) = X13/QJ(IJ,-1)
       RS(IJ,-1,3) = X14/QJ(IJ,-1)
       RS(IJ,-1,4) = X15/QJ(IJ,-1)
       RS(IJ,-1,5) = X16/QJ(IJ,-1)
c
      RUS(IJ,-1,1) = X11*XGAMAU*U(IJ,-1)*X12/QJ(IJ,-1)
      RUS(IJ,-1,2) = X13*(X1+X2)/QJ(IJ,-1)
      RUS(IJ,-1,3) = X14*(X1-X2)/QJ(IJ,-1)
      RUS(IJ,-1,4) = X15*(X1+X7)/QJ(IJ,-1)
      RUS(IJ,-1,5) = X16*(X1-X7)/QJ(IJ,-1)
C
      RVS(IJ,-1,1) = X11*XGAMAU*V(IJ,-1)*X12/QJ(IJ,-1)
      RVS(IJ,-1,2) = X13*(X8+X7)/QJ(IJ,-1)
      RVS(IJ,-1,3) = X14*(X8-X7)/QJ(IJ,-1)
      RVS(IJ,-1,4) = X15*(X8+X9)/QJ(IJ,-1)
      RVS(IJ,-1,5) = X16*(X8-X9)/QJ(IJ,-1)
C
      RES(IJ,-1,1) = X11*XGAMAU*X12/QJ(IJ,-1)
      RES(IJ,-1,2) = X17/QJ(IJ,-1)
      RES(IJ,-1,3) = X18/QJ(IJ,-1)
      RES(IJ,-1,4) = X19/QJ(IJ,-1)
      RES(IJ,-1,5) = X20/QJ(IJ,-1)
      ENDIF
   10 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE BMUPW(RS,US,VS,DRS)
C
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /CNST2/ R0,U0,V0,P0,T0,E0
      COMMON /CNST4/ DR0,UR0,VR0,PR0,TR0,ER0,
     1               DL0,UL0,VL0,PL0,TL0,EL0
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      DIMENSION RS(-1:NI2,-1:NJ2,5)
      DIMENSION DRS(-1:NI2,-1:NJ2,5)
      DIMENSION US(-1:NI2,-1:NJ2,5),VS(-1:NI2,-1:NJ2,5)
C
      DO 20 IJ     = -1, MIJ2
      I            = MOD(IJ+1,MI)-1
      J            = (IJ+1)/MI -1
      IF( I.GE.0.AND.I.LE.NI1.AND.J.GE.0.AND.J.LE.NJ1 )THEN
      USP1  = 0.5 *(US(IJ,-1,1)+ABS(US(IJ,-1,1)))*DTCFL
      USP2  = 0.5 *(US(IJ,-1,2)+ABS(US(IJ,-1,2)))*DTCFL
      USP3  = 0.5 *(US(IJ,-1,3)+ABS(US(IJ,-1,3)))*DTCFL
      USP4  = 0.5 *(US(IJ,-1,4)+ABS(US(IJ,-1,4)))*DTCFL
      USP5  = 0.5 *(US(IJ,-1,5)+ABS(US(IJ,-1,5)))*DTCFL
      USM1  = 0.5 *(US(IJ,-1,1)-ABS(US(IJ,-1,1)))*DTCFL
      USM2  = 0.5 *(US(IJ,-1,2)-ABS(US(IJ,-1,2)))*DTCFL
      USM3  = 0.5 *(US(IJ,-1,3)-ABS(US(IJ,-1,3)))*DTCFL
      USM4  = 0.5 *(US(IJ,-1,4)-ABS(US(IJ,-1,4)))*DTCFL
      USM5  = 0.5 *(US(IJ,-1,5)-ABS(US(IJ,-1,5)))*DTCFL
c
      VSP1  = 0.5 *(VS(IJ,-1,1)+ABS(VS(IJ,-1,1)))*DTCFL
      VSP2  = 0.5 *(VS(IJ,-1,2)+ABS(VS(IJ,-1,2)))*DTCFL
      VSP3  = 0.5 *(VS(IJ,-1,3)+ABS(VS(IJ,-1,3)))*DTCFL
      VSP4  = 0.5 *(VS(IJ,-1,4)+ABS(VS(IJ,-1,4)))*DTCFL
      VSP5  = 0.5 *(VS(IJ,-1,5)+ABS(VS(IJ,-1,5)))*DTCFL
      VSM1  = 0.5 *(VS(IJ,-1,1)-ABS(VS(IJ,-1,1)))*DTCFL
      VSM2  = 0.5 *(VS(IJ,-1,2)-ABS(VS(IJ,-1,2)))*DTCFL
      VSM3  = 0.5 *(VS(IJ,-1,3)-ABS(VS(IJ,-1,3)))*DTCFL
      VSM4  = 0.5 *(VS(IJ,-1,4)-ABS(VS(IJ,-1,4)))*DTCFL
      VSM5  = 0.5 *(VS(IJ,-1,5)-ABS(VS(IJ,-1,5)))*DTCFL
C
      DRS(I,J,1)   = DRS(I,J,1)
     1               - (USP1*RS(IJ,-1,1)-USM1*RS(IJ,-1,1)+
     1                  VSP1*RS(IJ,-1,1)-VSM1*RS(IJ,-1,1))
      DRS(I,J,2)   = DRS(I,J,2)
     1               - (USP2*RS(IJ,-1,2)-USM2*RS(IJ,-1,2)+
     1                  VSP2*RS(IJ,-1,2)-VSM2*RS(IJ,-1,2))
      DRS(I,J,3)   = DRS(I,J,3)
     1               - (USP3*RS(IJ,-1,3)-USM3*RS(IJ,-1,3)+
     1                  VSP3*RS(IJ,-1,3)-VSM3*RS(IJ,-1,3))
      DRS(I,J,4)   = DRS(I,J,4)
     1               - (USP4*RS(IJ,-1,4)-USM4*RS(IJ,-1,4)+
     1                  VSP4*RS(IJ,-1,4)-VSM4*RS(IJ,-1,4))
      DRS(I,J,5)   = DRS(I,J,5)
     1               - (USP5*RS(IJ,-1,5)-USM5*RS(IJ,-1,5)+
     1                  VSP5*RS(IJ,-1,5)-VSM5*RS(IJ,-1,5))
C
      DRS(I+1,J ,1) = DRS(I+1,J  ,1) +  USP1*RS(IJ,-1,1)
      DRS(I+1,J ,2) = DRS(I+1,J  ,2) +  USP2*RS(IJ,-1,2)
      DRS(I+1,J ,3) = DRS(I+1,J  ,3) +  USP3*RS(IJ,-1,3)
      DRS(I+1,J ,4) = DRS(I+1,J  ,4) +  USP4*RS(IJ,-1,4)
      DRS(I+1,J ,5) = DRS(I+1,J  ,5) +  USP5*RS(IJ,-1,5)
      DRS(I-1,J ,1) = DRS(I-1,J  ,1) -  USM1*RS(IJ,-1,1)
      DRS(I-1,J ,2) = DRS(I-1,J  ,2) -  USM2*RS(IJ,-1,2)
      DRS(I-1,J ,3) = DRS(I-1,J  ,3) -  USM3*RS(IJ,-1,3)
      DRS(I-1,J ,4) = DRS(I-1,J  ,4) -  USM4*RS(IJ,-1,4)
      DRS(I-1,J ,5) = DRS(I-1,J  ,5) -  USM5*RS(IJ,-1,5)
      DRS(I  ,J+1 ,1) = DRS(I  ,J+1,1) +  VSP1*RS(IJ,-1,1)
      DRS(I  ,J+1 ,2) = DRS(I  ,J+1,2) +  VSP2*RS(IJ,-1,2)
      DRS(I  ,J+1 ,3) = DRS(I  ,J+1,3) +  VSP3*RS(IJ,-1,3)
      DRS(I  ,J+1 ,4) = DRS(I  ,J+1,4) +  VSP4*RS(IJ,-1,4)
      DRS(I  ,J+1 ,5) = DRS(I  ,J+1,5) +  VSP5*RS(IJ,-1,5)
      DRS(I  ,J-1 ,1) = DRS(I  ,J-1,1) -  VSM1*RS(IJ,-1,1)
      DRS(I  ,J-1 ,2) = DRS(I  ,J-1,2) -  VSM2*RS(IJ,-1,2)
      DRS(I  ,J-1 ,3) = DRS(I  ,J-1,3) -  VSM3*RS(IJ,-1,3)
      DRS(I  ,J-1 ,4) = DRS(I  ,J-1,4) -  VSM4*RS(IJ,-1,4)
      DRS(I  ,J-1 ,5) = DRS(I  ,J-1,5) -  VSM5*RS(IJ,-1,5)
      ENDIF
   20 CONTINUE
C
      RETURN
      END
      SUBROUTINE DTIME(US,VS)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
      COMMON /U1/ LTP
      DIMENSION   US(-1:NI2,-1:NJ2,5),  VS(-1:NI2,-1:NJ2,5)
C
      DTI = 9999.
      DTJ = 9999.
      DO 100  K =  1,5
      DO 100 IJ = -1,MIJ2
      I         = MOD(IJ+1,MI)-1
      J         = (IJ+1)/MI -1
      IF(I.GT.0.AND.I.LT.NI1.AND.J.GT.0.AND.J.LT.NJ1)THEN
      EI        = 1./(ABS(US(IJ,-1,K))+1.E-20)
      DTI       = MIN(DTI,EI)
      EJ        = 1./(ABS(VS(IJ,-1,K))+1.E-20)
      DTJ       = MIN(DTJ,EJ)
      ENDIF
  100 CONTINUE
C
      IF((MOD(ITER,IPRINT).EQ.0).OR.(ITER.EQ.NSTEP)
     1    .OR. LTP .EQ. 1) THEN
      WRITE(*,*)'TIME STEP FOR I  DTI =  ',DTI
      WRITE(*,*)'TIME STEP FOR J  DTJ =  ',DTJ
      ENDIF
C
      DT        = MIN(DTI,DTJ)
C
      IF(DTFIX .GT.0.) THEN
        CFL  = DTFIX/DT
        CFL1 = CFL/5.
        CFL2 = CFL
        WRITE(*,*)' USE FIX DT ,  DTFIX = ',DTFIX
        WRITE(*,*)'               CFL   = ',CFL
      ENDIF
      RETURN
      END
c
c
c
      SUBROUTINE FDITVD2(C,QJ,RG,DQG,L)
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /ENTR / EPS,EPSS,EPS2
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2)
      DIMENSION    C(-1:NI2,-1:NJ2)
C
      DIMENSION  FXG(-1:MAXIJP)
      DIMENSION  FMG(-1:MAXIJP)
      DIMENSION  FNG(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  DPG(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP)
      DIMENSION SCXG(-1:MAXIJP), CCPG(-1:MAXIJP),  CCMG(-1:MAXIJP)
      DIMENSION DDAG(-1:MAXIJP)
      DIMENSION  EAG(-1:MAXIJP)
      DIMENSION   EG(-1:MAXIJP)
C
      ETA = 0.
      XI  = 0.
      LSTD= 1
C
      DO 10 J = -1,NJ2
C
      DO 20 I = -1,NI2
      FXG(I)   = C(I,J)*RG(I,J)
   20 CONTINUE
C     
      DO 30 I = -1,NI1
      DPG(I)=2.*(RG(I+1,J)*QJ(I+1,J)-RG(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
      DPFG(I) = FXG(I+1)   - FXG(I)
   30 CONTINUE
      DPG(NI2)  = DPG(NI1)
      DPFG(NI2) = DPFG(NI1)
C
      DO 35 I = -1,NI1
      IF(DPG(I) .NE. 0.) THEN
      CXG(I) = DPFG(I)/DPG(I)
      ELSE
      CXG(I) = 0.5*(C(I+1,J) + C(I,J))
      ENDIF
   35 CONTINUE
      CXG(NI2) = C(NI2,J)
      CXG(-1)  = C(-1,J)
C
      DO 40 I = -1,NI2
      ACXG(I) = ABS(CXG(I))
      IF(CXG(I).GT.EPS)THEN
        SCXG(I) = 1.
        CCPG(I) = 1.
        CCMG(I) = 0.
      ELSEIF(CXG(I) .LT.-EPS) THEN
        SCXG(I) = -1.
        CCPG(I) = 0.
        CCMG(I) = 1.
      ELSE
        SCXG(I) = 0.
        CCPG(I) = 0.5
        CCMG(I) = 0.5
      ENDIF
      DDAG(I)    = 0.
      IF( LSTD .EQ. 0)THEN
      DDAG(I)    = DTCFL*ACXG(I)
      ENDIF
   40 CONTINUE
C
      DO 60 I = -1,NI2
      EAG(I)   = SCXG(I)*(1.      -    DDAG(I)     )*DPFG(I)/2.
   60 CONTINUE
C
      DO 110 I = 0,NI2
      E1G = EAG(I  )
      E2G = EAG(I-1)
      A1G = -1.
      A2G = -1.
      IF(E1G.GE.0.)THEN
      A1G =1.
      ENDIF
      IF(E2G.GE.0.)THEN
      A2G =1.
      ENDIF
      EG(I) = 0.5*(A1G+A2G)*MIN(A1G*E1G,A2G*E2G)
  110 CONTINUE
      EG(-1) = EG(0)
C
      DO 120 I = -1,NI2
      FMG(I) = FXG(I) + EG(I)
  120 CONTINUE
C
      DO 130 I = -1,NI1
      FNG(I) = FMG(I+1) - CCPG(I)*(FMG(I+1) - FMG(I))
  130 CONTINUE
      FNG(NI2) = FMG(NI1)
C
      DO 140 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTCFL*(FNG(I) - FNG(I-1))
  140 CONTINUE
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE FDJTVD2(C,QJ,RG,DQG,L)
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /ENTR / EPS,EPSS,EPS2
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2)
      DIMENSION    C(-1:NI2,-1:NJ2)
C
      DIMENSION  FXG(-1:MAXIJP)
      DIMENSION  FMG(-1:MAXIJP)
      DIMENSION  FNG(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  DPG(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP)
      DIMENSION SCXG(-1:MAXIJP), CCPG(-1:MAXIJP),  CCMG(-1:MAXIJP)
      DIMENSION DDAG(-1:MAXIJP)
      DIMENSION  EAG(-1:MAXIJP)
      DIMENSION   EG(-1:MAXIJP)
C
      ETA = 0.
      XI  = 0.
      LSTD= 1
C
      DO 10 I = -1,NI2
C
      DO 20 J = -1,NJ2
      FXG(J)   = C(I,J)*RG(I,J)
   20 CONTINUE
C     
      DO 30 J = -1,NJ1
      DPG(J)=2.*(RG(I,J+1)*QJ(I,J+1)-RG(I,J)*QJ(I,J))
     1           /(QJ(I,J+1)+QJ(I,J))
      DPFG(J) = FXG(J+1)   - FXG(J)
   30 CONTINUE
      DPG(NJ2)  = DPG(NJ1)
      DPFG(NJ2) = DPFG(NJ1)
C
      DO 35 J = -1,NJ1
      IF(DPG(J) .NE. 0.) THEN
      CXG(J) = DPFG(J)/DPG(J)
      ELSE
      CXG(J) = 0.5*(C(I,J+1) + C(I,J))
      ENDIF
   35 CONTINUE
      CXG(NJ2) = C(I,NJ2)
      CXG(-1)  = C(I,-1)
C
      DO 40 J = -1,NJ2
      ACXG(J) = ABS(CXG(J))
      IF(CXG(J).GT.EPS)THEN
        SCXG(J) = 1.
        CCPG(J) = 1.
        CCMG(J) = 0.
      ELSEIF(CXG(J) .LT.-EPS) THEN
        SCXG(J) = -1.
        CCPG(J) = 0.
        CCMG(J) = 1.
      ELSE
        SCXG(J) = 0.
        CCPG(J) = 0.5
        CCMG(J) = 0.5
      ENDIF
      DDAG(J)    = 0.
      IF( LSTD .EQ. 0)THEN
      DDAG(J)    = DTCFL*ACXG(J)
      ENDIF
   40 CONTINUE
C
      DO 60 J = -1,NJ2
      EAG(J)   = SCXG(J)*(1.      -    DDAG(J)     )*DPFG(J)/2.
   60 CONTINUE
C
      DO 110 J = 0,NJ2
      E1G = EAG(J  )
      E2G = EAG(J-1)
      A1G = -1.
      A2G = -1.
      IF(E1G.GE.0.)THEN
      A1G =1.
      ENDIF
      IF(E2G.GE.0.)THEN
      A2G =1.
      ENDIF
      EG(J) = 0.5*(A1G+A2G)*MIN(A1G*E1G,A2G*E2G)
  110 CONTINUE
      EG(-1) = EG(0)
C
      DO 120 J = -1,NJ2
      FMG(J) = FXG(J) + EG(J)
  120 CONTINUE
C
      DO 130 J = -1,NJ1
      FNG(J) = FMG(J+1) - CCPG(J)*(FMG(J+1) - FMG(J))
  130 CONTINUE
      FNG(NJ2) = FMG(NJ1)
C
      DO 140 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTCFL*(FNG(J) - FNG(J-1))
  140 CONTINUE
   10 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE  FIEEN2(C,QJ,RG,DQG)
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2,5)
      DIMENSION    C(-1:NI2,-1:NJ2,5)
C
      DIMENSION  FP(-1:MAXIJP),FM(-1:MAXIJP)
      DIMENSION  DFP(-1:MAXIJP),DFM(-1:MAXIJP)
      DIMENSION  DDFP(-1:MAXIJP),DDFM(-1:MAXIJP)
      DIMENSION  XMBP(-1:MAXIJP),XMBM(-1:MAXIJP)
      DIMENSION  DIP(-1:MAXIJP),DIM(-1:MAXIJP)
      DIMENSION  FPP(-1:MAXIJP),FMM(-1:MAXIJP)
      DIMENSION  FG(-1:MAXIJP)
C
      DO 10 J = -1,NJ2
C
      AC      = 0.
      DO 15 I = -1,NI2
      ACTEMP  =  ABS(C(I,J,1))
      AC      =  MAX(AC,ACTEMP)
   15 CONTINUE
      DO 20 I = -1,NI2
C     AC       = ABS(C(I,J,1))
      FP(I)    = 0.5*(C(I,J,1)+AC)*RG(I,J,1)
      FM(I)    = 0.5*(C(I,J,1)-AC)*RG(I,J,1)
   20 CONTINUE
      DO 21 L =  2,5
      AC      = 0.
      DO 16 I = -1,NI2
      ACTEMP  =  ABS(C(I,J,L))
      AC      =  MAX(AC,ACTEMP)
   16 CONTINUE
      DO 21 I = -1,NI2
C     AC       = ABS(C(I,J,L))
      FP(I)    = FP(I) + 0.5*(C(I,J,L)+AC)*RG(I,J,L)
      FM(I)    = FM(I) + 0.5*(C(I,J,L)-AC)*RG(I,J,L)
21    CONTINUE
C    
      DO 24 I = -1,NI1
       DFP(I) = FP(I+1) - FP(I)
       DFM(I) = FM(I+1) - FM(I)
  24  CONTINUE
      DFP(NI2) = DFP(NI1)
      DFM(NI2) = DFM(NI1)
C
      DO 25 I = -1,NI1
       DDFP(I) = DFP(I+1) - DFP(I)
       DDFM(I) = DFM(I+1) - DFM(I)
  25  CONTINUE
      DDFP(NI2) = DDFP(NI1)
      DDFM(NI2) = DDFM(NI1)
C
      DO 28 I= 0,NI2
       IF( ABS(DDFP(I-1)) .LE. ABS(DDFP(I)) )THEN
        XMBP(I) = DDFP(I-1)
       ELSE
        XMBP(I) = DDFP(I)
       ENDIF
       IF( ABS(DDFM(I-1)) .LE. ABS(DDFM(I)) )THEN
        XMBM(I) = DDFM(I-1)
       ELSE
        XMBM(I) = DDFM(I)
       ENDIF
  28  CONTINUE
C
      XMBP(-1) = XMBP(0)
      XMBM(-1) = XMBM(0)
C
      DO 30 I = 0,NI2
       IF( ABS(DFP(I-1)) .LE. ABS(DFP(I)) )THEN
        DIP(I) = XMBP(I-1)/3.
       ELSE
        DIP(I) =-XMBP(I)/6.
       ENDIF
       IF( ABS(DFM(I-1)) .LE. ABS(DFM(I)) )THEN
        DIM(I) =-XMBM(I-1)/6.
       ELSE
        DIM(I) = XMBM(I)/3.
       ENDIF
30    CONTINUE
C
      DIP(-1) = DIP(0)
      DIM(-1) = DIM(0)
C
      DO 35 I = 0,NI2
       IF( ABS(DFP(I-1)) .LE. ABS(DFP(I)) )THEN
        XMBP(I) = DFP(I-1)/2.
       ELSE
        XMBP(I) = DFP(I)/2.
       ENDIF
       IF( ABS(DFM(I-1)) .LE. ABS(DFM(I)) )THEN
        XMBM(I) = DFM(I-1)/2.
       ELSE
        XMBM(I) = DFM(I)/2.
       ENDIF
35    CONTINUE
      XMBP(-1) = XMBP(0)
      XMBM(-1) = XMBM(0)
C
      DO 40 I = -1,NI2
       FPP(I) = FP(I) + XMBP(I) + DIP(I)
40    CONTINUE
C
      DO 50 I = -1,NI1
       FMM(I) = FM(I+1) - XMBM(I+1) + DIM(I+1)
50    CONTINUE
      FMM(NI2) = FMM(NI1)
C
      DO 60 I = -1,NI2
       FG(I) = FPP(I) + FMM(I)
60    CONTINUE
C
      DO 110 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTCFL*(FG(I) - FG(I-1))
  110 CONTINUE
   10 CONTINUE
c
      RETURN
      END
C
      SUBROUTINE  FIEEN3(C,QJ,RG,DQG)
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2,5)
      DIMENSION    C(-1:NI2,-1:NJ2,5)
C
      DIMENSION  FP(-1:MAXIJP),FM(-1:MAXIJP)
      DIMENSION  DFP(-1:MAXIJP),DFM(-1:MAXIJP)
      DIMENSION  DDFP(-1:MAXIJP),DDFM(-1:MAXIJP)
      DIMENSION  DDDFP(-1:MAXIJP),DDDFM(-1:MAXIJP)
      DIMENSION  XMBP(-1:MAXIJP),XMBM(-1:MAXIJP)
      DIMENSION  XMMBP(-1:MAXIJP),XMMBM(-1:MAXIJP)
      DIMENSION  DIP(-1:MAXIJP),DIM(-1:MAXIJP)
      DIMENSION  DDIP(-1:MAXIJP),DDIM(-1:MAXIJP)
      DIMENSION  FPP(-1:MAXIJP),FMM(-1:MAXIJP)
      DIMENSION  FG(-1:MAXIJP)
C
      DO 10 J = -1,NJ2
C
      AC      = 0.
      DO 15 I = -1,NI2
      ACTEMP  =  ABS(C(I,J,1))
      AC      =  MAX(AC,ACTEMP)
   15 CONTINUE
      DO 20 I = -1,NI2
C     AC       = ABS(C(I,J,1))
      FP(I)    = 0.5*(C(I,J,1)+AC)*RG(I,J,1)
      FM(I)    = 0.5*(C(I,J,1)-AC)*RG(I,J,1)
   20 CONTINUE
      DO 21 L =  2,5
      AC      = 0.
      DO 16 I = -1,NI2
      ACTEMP  =  ABS(C(I,J,L))
      AC      =  MAX(AC,ACTEMP)
   16 CONTINUE
      DO 21 I = -1,NI2
C     AC       = ABS(C(I,J,L))
      FP(I)    = FP(I) + 0.5*(C(I,J,L)+AC)*RG(I,J,L)
      FM(I)    = FM(I) + 0.5*(C(I,J,L)-AC)*RG(I,J,L)
21    CONTINUE
C    
      DO 24 I = -1,NI1
       DFP(I) = FP(I+1) - FP(I)
       DFM(I) = FM(I+1) - FM(I)
  24  CONTINUE
      DFP(NI2) = DFP(NI1)
      DFM(NI2) = DFM(NI1)
C
      DO 25 I = -1,NI1
       DDFP(I) = DFP(I+1) - DFP(I)
       DDFM(I) = DFM(I+1) - DFM(I)
  25  CONTINUE
      DDFP(NI2) = DDFP(NI1)
      DDFM(NI2) = DDFM(NI1)
C
      DO 26 I = -1,NI1
       DDDFP(I) = DDFP(I+1) - DDFP(I)
       DDDFM(I) = DDFM(I+1) - DDFM(I)
  26  CONTINUE
      DDDFP(NI2) = DDDFP(NI1)
      DDDFM(NI2) = DDDFM(NI1)
C
      DO 28 I= 0,NI2
       IF( ABS(DDFP(I-1)) .LE. ABS(DDFP(I)) )THEN
        XMBP(I) = DDFP(I-1)
       ELSE
        XMBP(I) = DDFP(I)
       ENDIF
       IF( ABS(DDFM(I-1)) .LE. ABS(DDFM(I)) )THEN
        XMBM(I) = DDFM(I-1)
       ELSE
        XMBM(I) = DDFM(I)
       ENDIF
  28  CONTINUE
C
      XMBP(-1) = XMBP(0)
      XMBM(-1) = XMBM(0)
C
      DO 31 I= 0,NI2
       IF( ABS(DDDFP(I-1)) .LE. ABS(DDDFP(I)) )THEN
        XMMBP(I) = DDDFP(I-1)
       ELSE
        XMMBP(I) = DDDFP(I)
       ENDIF
       IF( ABS(DDDFM(I-1)) .LE. ABS(DDDFM(I)) )THEN
        XMMBM(I) = DDDFM(I-1)
       ELSE
        XMMBM(I) = DDDFM(I)
       ENDIF
  31  CONTINUE
C
      XMMBP(-1) = XMMBP(0)
      XMMBM(-1) = XMMBM(0)
C
      DO 32 I = 0,NI2
       IF( ABS(DFP(I-1)) .LE. ABS(DFP(I)) )THEN
        DIP(I) = XMBP(I-1)/3.
       ELSE
        DIP(I) =-XMBP(I)/6.
       ENDIF
       IF( ABS(DFM(I-1)) .LE. ABS(DFM(I)) )THEN
        DIM(I) =-XMBM(I-1)/6.
       ELSE
        DIM(I) = XMBM(I)/3.
       ENDIF
32    CONTINUE
C
      DIP(-1) = DIP(0)
      DIM(-1) = DIM(0)
C
      DO 33 I = 0,NI2
       IF( ABS(DFP(I-1)) .LE. ABS(DFP(I)) )THEN
        XMBP(I) = DFP(I-1)/2.
       ELSE
        XMBP(I) = DFP(I)/2.
       ENDIF
       IF( ABS(DFM(I-1)) .LE. ABS(DFM(I)) )THEN
        XMBM(I) = DFM(I-1)/2.
       ELSE
        XMBM(I) = DFM(I)/2.
       ENDIF
33    CONTINUE
      XMBP(-1) = XMBP(0)
      XMBM(-1) = XMBM(0)
C
      DO 35 I = 0,NI1
      IF( ABS(DFP(I-1)) .LE. ABS(DFP(I)) )THEN
       IF( ABS(DDFP(I-1)) .LE. ABS(DDFP(I)) )THEN
        DDIP(I) = XMMBP(I-1)/4.
       ELSE
        DDIP(I) =-XMMBP(I)/12.
       ENDIF
      ELSE
       IF( ABS(DDFP(I)) .LE. ABS(DDFP(I+1)) )THEN
        DDIP(I) =-XMMBP(I)/12.
       ELSE
        DDIP(I) = XMMBP(I+1)/12.
       ENDIF
      ENDIF
      IF( ABS(DFM(I-1)) .LE. ABS(DFM(I)) )THEN
       IF( ABS(DDFM(I-1)) .LE. ABS(DDFM(I)) )THEN
        DDIM(I) =-XMMBM(I-1)/12.
       ELSE
        DDIM(I) = XMMBM(I)/12.
       ENDIF
      ELSE
       IF( ABS(DDFM(I)) .LE. ABS(DDFM(I+1)) )THEN
        DDIM(I) = XMMBM(I)/12.
       ELSE
        DDIM(I) =-XMMBM(I+1)/4.
       ENDIF
      ENDIF
35    CONTINUE
C
      DDIP(-1)  = DDIP(0)
      DDIP(NI2) = DDIP(NI1)
      DDIM(-1) = DDIM(0)
      DDIM(NI2) = DDIM(NI1)
C
      DO 40 I = -1,NI2
       FPP(I) = FP(I) + XMBP(I) + DIP(I) + DDIP(I)
40    CONTINUE
C
      DO 50 I = -1,NI1
       FMM(I) = FM(I+1) - XMBM(I+1) + DIM(I+1) + DDIM(I+1)
50    CONTINUE
      FMM(NI2) = FMM(NI1)
C
      DO 60 I = -1,NI2
       FG(I) = FPP(I) + FMM(I)
60    CONTINUE
C
      DO 110 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTCFL*(FG(I) - FG(I-1))
  110 CONTINUE
   10 CONTINUE
c
      RETURN
      END
C
      SUBROUTINE FIWEN2(C,QJ,RG,DQG)
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2,5)
      DIMENSION    C(-1:NI2,-1:NJ2,5)
C
      DIMENSION  FP(-1:MAXIJP),FM(-1:MAXIJP)
      DIMENSION  QP1(-1:MAXIJP),QP2(-1:MAXIJP)
      DIMENSION  QM1(-1:MAXIJP),QM2(-1:MAXIJP)
      DIMENSION  AP0(-1:MAXIJP),AP1(-1:MAXIJP)
      DIMENSION  AM0(-1:MAXIJP),AM1(-1:MAXIJP)
      DIMENSION    FG(-1:MAXIJP)
C
      DO 10 J = -1,NJ2
C
      AC = 1.E-6
      DO 15 I = -1,NI2
      ACTEMP   = ABS(C(I,J,1))
      AC       = MAX(AC,ACTEMP)
   15 CONTINUE
      AC       = AC*1.01
      DO 20 I = -1,NI2
C      AC       = ABS(C(I,J,1))
      FP(I)    = 0.5*(C(I,J,1)+AC)*RG(I,J,1)
      FM(I)    = 0.5*(C(I,J,1)-AC)*RG(I,J,1)
   20 CONTINUE
      DO 21 L =  2,5
      AC = 1.E-6
      DO 16 I = -1,NI2
      ACTEMP   = ABS(C(I,J,L))
      AC       = MAX(AC,ACTEMP)
   16 CONTINUE
      AC       = AC*1.01
      DO 21 I = -1,NI2
C      AC       = ABS(C(I,J,L))
      FP(I)    = FP(I) + 0.5*(C(I,J,L)+AC)*RG(I,J,L)
      FM(I)    = FM(I) + 0.5*(C(I,J,L)-AC)*RG(I,J,L)
21    CONTINUE
C    
        DO 24 I=0,NI1
          QP1(I) = (3.*FP(I)-FP(I-1))/2.
          QP2(I) = (FP(I+1)+FP(I))/2.
          QM1(I) = (FM(I+1)+FM(I))/2.
          SP0    = (FP(I)-FP(I-1))*(FP(I)-FP(I-1))
          SP1    = (FP(I+1)-FP(I))*(FP(I+1)-FP(I))
          SM0    = (FM(I+1)-FM(I))*(FM(I+1)-FM(I))
          AP0(I) = 1./3./(EPS+SP0)/(EPS+SP0)
          AP1(I) = 2./3./(EPS+SP1)/(EPS+SP1)
          AM0(I) = 2./3./(EPS+SM0)/(EPS+SM0)
  24    CONTINUE
C        QP1(-1) = (FP(0)+FP(-1))/2.
C        QP2(-1) = (FP(0)+FP(-1))/2.
C        QM1(-1) = (FM(0)+FM(-1))/2.
        QP1(-1) = QP1(0)
        QP2(-1) = QP2(0)
        QM1(-1) = QM1(0)
        AP0(-1)  = AP0(0)
C        AP1(-1)  = AP0(0)
        AP1(-1)  = AP1(0)
C        SM0    = (FM(0)-FM(-1))*(FM(0)-FM(-1))
C        AM0(-1) = 2./3./(EPS+SM0)/(EPS+SM0)
        AM0(-1)  = AM0(0)
        QP1(NI2) = QP1(NI1)
        QP2(NI2) = QP2(NI1)
        QM1(NI2) = QM1(NI1)
        AP0(NI2) = AP0(NI1)
        AP1(NI2) = AP1(NI1)
        AM0(NI2) = AM0(NI1)
C
        DO 25 I=-1,NI
          QM2(I) = (3.*FM(I+1)-FM(I+2))/2.
          SM1    = (FM(I+2)-FM(I+1))*(FM(I+2)-FM(I+1))
          AM1(I) = 1./3./(EPS+SM1)/(EPS+SM1)
  25    CONTINUE
        QM2(NI1) = QM2(NI)
        AM1(NI1) = AM1(NI)
        QM2(NI2) = QM2(NI)
        AM1(NI2) = AM1(NI)
C
        DO 28 I=-1,NI2
         WP    = AP0(I)+AP1(I)
         WM    = AM0(I)+AM1(I)
         WP0   = AP0(I)/WP
         WP1   = AP1(I)/WP
         WM0   = AM0(I)/WM
         WM1   = AM1(I)/WM
         FPP   = WP0*QP1(I)+WP1*QP2(I)
         FMM   = WM0*QM1(I)+WM1*QM2(I)
         FG(I) = FPP + FMM
  28    CONTINUE
C
      DO 110 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTCFL*(FG(I) - FG(I-1))
  110 CONTINUE
   10 CONTINUE
c
      RETURN
      END
C
      SUBROUTINE  FIWEN3(C,QJ,RG,DQG)
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2,5)
      DIMENSION    C(-1:NI2,-1:NJ2,5)
C
      DIMENSION  FP(-1:MAXIJP),FM(-1:MAXIJP)
      DIMENSION  QP1(-1:MAXIJP),QP2(-1:MAXIJP),QP3(-1:MAXIJP)
      DIMENSION  QM1(-1:MAXIJP),QM2(-1:MAXIJP),QM3(-1:MAXIJP)
      DIMENSION  AP0(-1:MAXIJP),AP1(-1:MAXIJP),AP2(-1:MAXIJP)
      DIMENSION  AM0(-1:MAXIJP),AM1(-1:MAXIJP),AM2(-1:MAXIJP)
      DIMENSION    FG(-1:MAXIJP)
C
      DO 10 J = -1,NJ2
C
      AC = 1.E-6
      DO 15 I = -1,NI2
      ACTEMP  = ABS(C(I,J,1))
      AC    = MAX(AC,ACTEMP)
  15  CONTINUE
      AC = 1.01*AC
      DO 20 I = -1,NI2
C      AC       = ABS(C(I,J,1))
      FP(I)    = 0.5*(C(I,J,1)+AC)*RG(I,J,1)
      FM(I)    = 0.5*(C(I,J,1)-AC)*RG(I,J,1)
   20 CONTINUE
      DO 21 L =  2,5
      AC   = 1.E-6
      DO 16 I  = -1,NI2
      ACTEMP = ABS(C(I,J,L))
      AC     = MAX(AC,ACTEMP)
   16 CONTINUE
      AC = 1.01*AC
      DO 21 I = -1,NI2
C      AC       = ABS(C(I,J,L))
      FP(I)    = FP(I) + 0.5*(C(I,J,L)+AC)*RG(I,J,L)
      FM(I)    = FM(I) + 0.5*(C(I,J,L)-AC)*RG(I,J,L)
  21  CONTINUE
C    
        DO 24 I =  1,NI
          QP1(I) = (2.*FP(I-2)-7.*FP(I-1)+11.*FP(I))/6.
          QP2(I) = (-FP(I-1)+5.*FP(I)+2.*FP(I+1))/6.
          QP3(I) = (2.*FP(I)+5.*FP(I+1)-FP(I+2))/6.
          QM1(I) = (-FM(I-1)+5.*FM(I)+2.*FM(I+1))/6.
          QM2(I) = (2.*FM(I)+5.*FM(I+1)-FM(I+2))/6.
          SP0    = 13./12.*(FP(I-2)-2.*FP(I-1)+FP(I))*
     1                     (FP(I-2)-2.*FP(I-1)+FP(I))+
     2                0.25*(FP(I-2)-4.*FP(I-1)+3.*FP(I))*
     3                     (FP(I-2)-4.*FP(I-1)+3.*FP(I))
          SP1    = 13./12.*(FP(I-1)-2.*FP(I)+FP(I+1))*
     1                     (FP(I-1)-2.*FP(I)+FP(I+1))+
     2                0.25*(FP(I-1)-FP(I+1))*(FP(I-1)-FP(I+1))
          SP2    = 13./12.*(FP(I+2)-2.*FP(I+1)+FP(I))*
     1                     (FP(I+2)-2.*FP(I+1)+FP(I))+
     2                0.25*(FP(I+2)-4.*FP(I+1)+3.*FP(I))*
     3                     (FP(I+2)-4.*FP(I+1)+3.*FP(I))
          SM0    = 13./12.*(FM(I-1)-2.*FM(I)+FM(I+1))*
     1                     (FM(I-1)-2.*FM(I)+FM(I+1))+
     2                0.25*(FM(I-1)-4.*FM(I)+3.*FM(I+1))*
     3                     (FM(I-1)-4.*FM(I)+3.*FM(I+1))
          SM1    = 13./12.*(FM(I)-2.*FM(I+1)+FM(I+2))*
     1                     (FM(I)-2.*FM(I+1)+FM(I+2))+
     2                0.25*(FM(I)-FM(I+2))*(FM(I)-FM(I+2))
          AP0(I) = 0.1/(EPS+SP0)/(EPS+SP0)
          AP1(I) = 0.6/(EPS+SP1)/(EPS+SP1)
          AP2(I) = 0.3/(EPS+SP2)/(EPS+SP2)
          AM0(I) = 0.3/(EPS+SM0)/(EPS+SM0)
          AM1(I) = 0.6/(EPS+SM1)/(EPS+SM1)
  24    CONTINUE
        QP1(-1)  = QP1(1)
        QP2(-1)  = QP2(1)
        QP3(-1)  = QP3(1)
        QM1(-1)  = QM1(1)
        QM2(-1)  = QM2(1)
        AP0(-1)  = AP0(1)
        AP1(-1)  = AP1(1)
        AP2(-1)  = AP2(1)
        AM0(-1)  = AM0(1)
        AM1(-1)  = AM1(1)
        QP1( 0)  = QP1(1)
        QP2( 0)  = QP2(1)
        QP3( 0)  = QP3(1)
        QM1( 0)  = QM1(1)
        QM2( 0)  = QM2(1)
        AP0( 0)  = AP0(1)
        AP1( 0)  = AP1(1)
        AP2( 0)  = AP2(1)
        AM0( 0)  = AM0(1)
        AM1( 0)  = AM1(1)
        QP1(NI2) = QP1(NI)
        QP2(NI2) = QP2(NI)
        QP3(NI2) = QP3(NI)
        QM1(NI2) = QM1(NI)
        QM2(NI2) = QM2(NI)
        AP0(NI2) = AP0(NI)
        AP1(NI2) = AP1(NI)
        AP2(NI2) = AP2(NI)
        AM0(NI2) = AM0(NI)
        AM1(NI2) = AM1(NI)
        QP1(NI1) = QP1(NI)
        QP2(NI1) = QP2(NI)
        QP3(NI1) = QP3(NI)
        QM1(NI1) = QM1(NI)
        QM2(NI1) = QM2(NI)
        AP0(NI1) = AP0(NI)
        AP1(NI1) = AP1(NI)
        AP2(NI1) = AP2(NI)
        AM0(NI1) = AM0(NI)
        AM1(NI1) = AM1(NI)
C
        DO 25 I=-1,NI-1
          QM3(I) = (11.*FM(I+1)-7.*FM(I+2)+2.*FM(I+3))/6.
          SM2    = 13./12.*(FM(I+1)-2.*FM(I+2)+FM(I+3))*
     1                     (FM(I+1)-2.*FM(I+2)+FM(I+3))+
     2                0.25*(3.*FM(I+1)-4.*FM(I+2)+FM(I+3))*
     3                     (3.*FM(I+1)-4.*FM(I+2)+FM(I+3))
          AM2(I) = 0.1/(EPS+SM2)/(EPS+SM2)
  25    CONTINUE
        QM3(NI ) = QM3(NI-1)
        AM2(NI ) = AM2(NI-1)
        QM3(NI1) = QM3(NI-1)
        AM2(NI1) = AM2(NI-1)
        QM3(NI2) = QM3(NI-1)
        AM2(NI2) = AM2(NI-1)
C
        DO 28 I=-1,NI2
         WP    = AP0(I)+AP1(I)+AP2(I)
         WM    = AM0(I)+AM1(I)+AM2(I)
         WP0   = AP0(I)/WP
         WP1   = AP1(I)/WP
         WP2   = AP2(I)/WP
         WM0   = AM0(I)/WM
         WM1   = AM1(I)/WM
         WM2   = AM2(I)/WM
         FPP   = WP0*QP1(I)+WP1*QP2(I)+WP2*QP3(I)
         FMM   = WM0*QM1(I)+WM1*QM2(I)+WM2*QM3(I)
         FG(I) = FPP + FMM
  28    CONTINUE
C
      DO 110 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTCFL*(FG(I) - FG(I-1))
  110 CONTINUE
   10 CONTINUE
c
      RETURN
      END
C
      SUBROUTINE  FJEEN2(C,QJ,RG,DQG)
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2,5)
      DIMENSION    C(-1:NI2,-1:NJ2,5)
C
      DIMENSION  FP(-1:MAXIJP),FM(-1:MAXIJP)
      DIMENSION  DFP(-1:MAXIJP),DFM(-1:MAXIJP)
      DIMENSION  DDFP(-1:MAXIJP),DDFM(-1:MAXIJP)
      DIMENSION  XMBP(-1:MAXIJP),XMBM(-1:MAXIJP)
      DIMENSION  DIP(-1:MAXIJP),DIM(-1:MAXIJP)
      DIMENSION  FPP(-1:MAXIJP),FMM(-1:MAXIJP)
      DIMENSION  FG(-1:MAXIJP)
C
      DO 10 I = -1,NI2
C
      AC = 0.
      DO 15 J = -1,NJ2
      ACTEMP = ABS(C(I,J,1))
      AC     = MAX(AC,ACTEMP)
  15  CONTINUE
      DO 20 J = -1,NJ2
C      AC       = ABS(C(I,J,1))
      FP(J)    = 0.5*(C(I,J,1)+AC)*RG(I,J,1)
      FM(J)    = 0.5*(C(I,J,1)-AC)*RG(I,J,1)
   20 CONTINUE
      DO 21 L =  2,5
      AC = 0.
      DO 16 J = -1,NJ2
      ACTEMP   = ABS(C(I,J,L))
      AC      = MAX(AC,ACTEMP)
  16  CONTINUE
      DO 21 J = -1,NJ2
C      AC       = ABS(C(I,J,L))
      FP(J)    = FP(J) + 0.5*(C(I,J,L)+AC)*RG(I,J,L)
      FM(J)    = FM(J) + 0.5*(C(I,J,L)-AC)*RG(I,J,L)
21    CONTINUE
C    
      DO 24 J = -1,NJ1
       DFP(J) = FP(J+1) - FP(J)
       DFM(J) = FM(J+1) - FM(J)
  24  CONTINUE
      DFP(NJ2) = DFP(NJ1)
      DFM(NJ2) = DFM(NJ1)
C
      DO 25 J = -1,NJ1
       DDFP(J) = DFP(J+1) - DFP(J)
       DDFM(J) = DFM(J+1) - DFM(J)
  25  CONTINUE
      DDFP(NJ2) = DDFP(NJ1)
      DDFM(NJ2) = DDFM(NJ1)
C
      DO 28 J= 0,NJ2
       IF( ABS(DDFP(J-1)) .LE. ABS(DDFP(J)) )THEN
        XMBP(J) = DDFP(J-1)
       ELSE
        XMBP(J) = DDFP(J)
       ENDIF
       IF( ABS(DDFM(J-1)) .LE. ABS(DDFM(J)) )THEN
        XMBM(J) = DDFM(J-1)
       ELSE
        XMBM(J) = DDFM(J)
       ENDIF
  28  CONTINUE
C
      XMBP(-1) = XMBP(0)
      XMBM(-1) = XMBM(0)
C
      DO 30 J = 0,NJ2
       IF( ABS(DFP(J-1)) .LE. ABS(DFP(J)) )THEN
        DIP(J) = XMBP(J-1)/3.
       ELSE
        DIP(J) =-XMBP(J)/6.
       ENDIF
       IF( ABS(DFM(J-1)) .LE. ABS(DFM(J)) )THEN
        DIM(J) =-XMBM(J-1)/6.
       ELSE
        DIM(J) = XMBM(J)/3.
       ENDIF
30    CONTINUE
C
      DIP(-1) = DIP(0)
      DIM(-1) = DIM(0)
C
      DO 35 J = 0,NJ2
       IF( ABS(DFP(J-1)) .LE. ABS(DFP(J)) )THEN
        XMBP(J) = DFP(J-1)/2.
       ELSE
        XMBP(J) = DFP(J)/2.
       ENDIF
       IF( ABS(DFM(J-1)) .LE. ABS(DFM(J)) )THEN
        XMBM(J) = DFM(J-1)/2.
       ELSE
        XMBM(J) = DFM(J)/2.
       ENDIF
35    CONTINUE
      XMBP(-1) = XMBP(0)
      XMBM(-1) = XMBM(0)
C
      DO 40 J = -1,NJ2
       FPP(J) = FP(J) + XMBP(J) + DIP(J)
40    CONTINUE
C
      DO 50 J = -1,NJ1
       FMM(J) = FM(J+1) - XMBM(J+1) + DIM(J+1)
50    CONTINUE
      FMM(NJ2) = FMM(NJ1)
C
      DO 60 J = -1,NJ2
       FG(J) = FPP(J) + FMM(J)
60    CONTINUE
C
      DO 110 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTCFL*(FG(J) - FG(J-1))
  110 CONTINUE
   10 CONTINUE
c
      RETURN
      END
C
      SUBROUTINE  FJEEN3(C,QJ,RG,DQG)
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2,5)
      DIMENSION    C(-1:NI2,-1:NJ2,5)
C
      DIMENSION  FP(-1:MAXIJP),FM(-1:MAXIJP)
      DIMENSION  DFP(-1:MAXIJP),DFM(-1:MAXIJP)
      DIMENSION  DDFP(-1:MAXIJP),DDFM(-1:MAXIJP)
      DIMENSION  DDDFP(-1:MAXIJP),DDDFM(-1:MAXIJP)
      DIMENSION  XMBP(-1:MAXIJP),XMBM(-1:MAXIJP)
      DIMENSION  XMMBP(-1:MAXIJP),XMMBM(-1:MAXIJP)
      DIMENSION  DIP(-1:MAXIJP),DIM(-1:MAXIJP)
      DIMENSION  DDIP(-1:MAXIJP),DDIM(-1:MAXIJP)
      DIMENSION  FPP(-1:MAXIJP),FMM(-1:MAXIJP)
      DIMENSION  FG(-1:MAXIJP)
C
      DO 10 I = -1,NI2
C
      AC  = 0.
      DO 15 J = -1,NJ2
      ACTEMP = ABS(C(I,J,1))
      AC     = MAX(AC,ACTEMP)
   15 CONTINUE
      DO 20 J = -1,NJ2
C      AC       = ABS(C(I,J,1))
      FP(J)    = 0.5*(C(I,J,1)+AC)*RG(I,J,1)
      FM(J)    = 0.5*(C(I,J,1)-AC)*RG(I,J,1)
   20 CONTINUE
      DO 21 L =  2,5
      AC = 0.
      DO 16 J = -1,NJ2
      ACTEMP = ABS(C(I,J,L))
      AC     = MAX(AC,ACTEMP)
   16 CONTINUE
      DO 21 J = -1,NJ2
C      AC       = ABS(C(I,J,L))
      FP(J)    = FP(J) + 0.5*(C(I,J,L)+AC)*RG(I,J,L)
      FM(J)    = FM(J) + 0.5*(C(I,J,L)-AC)*RG(I,J,L)
21    CONTINUE
C    
      DO 24 J = -1,NJ1
       DFP(J) = FP(J+1) - FP(J)
       DFM(J) = FM(J+1) - FM(J)
  24  CONTINUE
      DFP(NJ2) = DFP(NJ1)
      DFM(NJ2) = DFM(NJ1)
C
      DO 25 J = -1,NJ1
       DDFP(J) = DFP(J+1) - DFP(J)
       DDFM(J) = DFM(J+1) - DFM(J)
  25  CONTINUE
      DDFP(NJ2) = DDFP(NJ1)
      DDFM(NJ2) = DDFM(NJ1)
C
      DO 26 J = -1,NJ1
       DDDFP(J) = DDFP(J+1) - DDFP(J)
       DDDFM(J) = DDFM(J+1) - DDFM(J)
  26  CONTINUE
      DDDFP(NJ2) = DDDFP(NJ1)
      DDDFM(NJ2) = DDDFM(NJ1)
C
      DO 28 J= 0,NJ2
       IF( ABS(DDFP(J-1)) .LE. ABS(DDFP(J)) )THEN
        XMBP(J) = DDFP(J-1)
       ELSE
        XMBP(J) = DDFP(J)
       ENDIF
       IF( ABS(DDFM(J-1)) .LE. ABS(DDFM(J)) )THEN
        XMBM(J) = DDFM(J-1)
       ELSE
        XMBM(J) = DDFM(J)
       ENDIF
  28  CONTINUE
C
      XMBP(-1) = XMBP(0)
      XMBM(-1) = XMBM(0)
C
      DO 31 J= 0,NJ2
       IF( ABS(DDDFP(J-1)) .LE. ABS(DDDFP(J)) )THEN
        XMMBP(J) = DDDFP(J-1)
       ELSE
        XMMBP(J) = DDDFP(J)
       ENDIF
       IF( ABS(DDDFM(J-1)) .LE. ABS(DDDFM(J)) )THEN
        XMMBM(J) = DDDFM(J-1)
       ELSE
        XMMBM(J) = DDDFM(J)
       ENDIF
  31  CONTINUE
C
      XMMBP(-1) = XMMBP(0)
      XMMBM(-1) = XMMBM(0)
C
      DO 32 J = 0,NJ2
       IF( ABS(DFP(J-1)) .LE. ABS(DFP(J)) )THEN
        DIP(J) = XMBP(J-1)/3.
       ELSE
        DIP(J) =-XMBP(J)/6.
       ENDIF
       IF( ABS(DFM(J-1)) .LE. ABS(DFM(J)) )THEN
        DIM(J) =-XMBM(J-1)/6.
       ELSE
        DIM(J) = XMBM(J)/3.
       ENDIF
32    CONTINUE
C
      DIP(-1) = DIP(0)
      DIM(-1) = DIM(0)
C
      DO 33 J = 0,NJ2
       IF( ABS(DFP(J-1)) .LE. ABS(DFP(J)) )THEN
        XMBP(J) = DFP(J-1)/2.
       ELSE
        XMBP(J) = DFP(J)/2.
       ENDIF
       IF( ABS(DFM(J-1)) .LE. ABS(DFM(J)) )THEN
        XMBM(J) = DFM(J-1)/2.
       ELSE
        XMBM(J) = DFM(J)/2.
       ENDIF
33    CONTINUE
      XMBP(-1) = XMBP(0)
      XMBM(-1) = XMBM(0)
C
      DO 35 J = 0,NJ1
      IF( ABS(DFP(J-1)) .LE. ABS(DFP(J)) )THEN
       IF( ABS(DDFP(J-1)) .LE. ABS(DDFP(J)) )THEN
        DDIP(J) = XMMBP(J-1)/4.
       ELSE
        DDIP(J) =-XMMBP(J)/12.
       ENDIF
      ELSE
       IF( ABS(DDFP(J)) .LE. ABS(DDFP(J+1)) )THEN
        DDIP(J) =-XMMBP(J)/12.
       ELSE
        DDIP(J) = XMMBP(J+1)/12.
       ENDIF
      ENDIF
      IF( ABS(DFM(J-1)) .LE. ABS(DFM(J)) )THEN
       IF( ABS(DDFM(J-1)) .LE. ABS(DDFM(J)) )THEN
        DDIM(J) =-XMMBM(J-1)/12.
       ELSE
        DDIM(J) = XMMBM(J)/12.
       ENDIF
      ELSE
       IF( ABS(DDFM(J)) .LE. ABS(DDFM(J+1)) )THEN
        DDIM(J) = XMMBM(J)/12.
       ELSE
        DDIM(J) =-XMMBM(J+1)/4.
       ENDIF
      ENDIF
35    CONTINUE
C
      DDIP(-1)  = DDIP(0)
      DDIP(NJ2) = DDIP(NJ1)
      DDIM(-1) = DDIM(0)
      DDIM(NJ2) = DDIM(NJ1)
C
      DO 40 J = -1,NJ2
       FPP(J) = FP(J) + XMBP(J) + DIP(J) + DDIP(J)
40    CONTINUE
C
      DO 50 J = -1,NJ1
       FMM(J) = FM(J+1) - XMBM(J+1) + DIM(J+1) + DDIM(J+1)
50    CONTINUE
      FMM(NJ2) = FMM(NJ1)
C
      DO 60 J = -1,NJ2
       FG(J) = FPP(J) + FMM(J)
60    CONTINUE
C
      DO 110 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTCFL*(FG(J) - FG(J-1))
  110 CONTINUE
   10 CONTINUE
c
      RETURN
      END
C
      SUBROUTINE  FJWEN2(C,QJ,RG,DQG)
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2,5)
      DIMENSION    C(-1:NI2,-1:NJ2,5)
C
      DIMENSION  FP(-1:MAXIJP),FM(-1:MAXIJP)
      DIMENSION  QP1(-1:MAXIJP),QP2(-1:MAXIJP)
      DIMENSION  QM1(-1:MAXIJP),QM2(-1:MAXIJP)
      DIMENSION  AP0(-1:MAXIJP),AP1(-1:MAXIJP)
      DIMENSION  AM0(-1:MAXIJP),AM1(-1:MAXIJP)
      DIMENSION    FG(-1:MAXIJP)
C
      DO 10 I = -1,NI2
C
      AC  = 1.E-6
      DO 15 J = -1,NJ2
      ACTEMP = ABS(C(I,J,1))
      AC     = MAX(AC,ACTEMP)
  15  CONTINUE
       AC    = AC*1.01
      DO 20 J = -1,NJ2
C      AC       = ABS(C(I,J,1))
      FP(J)    = 0.5*(C(I,J,1)+AC)*RG(I,J,1)
      FM(J)    = 0.5*(C(I,J,1)-AC)*RG(I,J,1)
   20 CONTINUE
      DO 21 L =  2,5
      AC  = 1.E-6
      DO 16 J = -1,NJ2
      ACTEMP = ABS(C(I,J,L))
      AC     = MAX(AC,ACTEMP)
   16 CONTINUE
       AC    = AC*1.01
      DO 21 J = -1,NJ2
C      AC       = ABS(C(I,J,L))
      FP(J)    = FP(J) + 0.5*(C(I,J,L)+AC)*RG(I,J,L)
      FM(J)    = FM(J) + 0.5*(C(I,J,L)-AC)*RG(I,J,L)
21    CONTINUE
C    
        DO 24 J = 0,NJ1
          QP1(J) = (3.*FP(J)-FP(J-1))/2.
          QP2(J) = (FP(J+1)+FP(J))/2.
          QM1(J) = (FM(J+1)+FM(J))/2.
          SP0    = (FP(J)-FP(J-1))*(FP(J)-FP(J-1))
          SP1    = (FP(J+1)-FP(J))*(FP(J+1)-FP(J))
          SM0    = (FM(J+1)-FM(J))*(FM(J+1)-FM(J))
          AP0(J) = 1./3./(EPS+SP0)/(EPS+SP0)
          AP1(J) = 2./3./(EPS+SP1)/(EPS+SP1)
          AM0(J) = 2./3./(EPS+SM0)/(EPS+SM0)
  24    CONTINUE
C        QP1(-1) = (FP(0)+FP(-1))/2.
C        QP2(-1) = (FP(0)+FP(-1))/2.
C        QM1(-1) = (FM(0)+FM(-1))/2.
        QP1(-1) = QP1(0)
        QP2(-1) = QP2(0)
        QM1(-1) = QM1(0)
        AP0(-1)  = AP0(0)
        AP1(-1)  = AP1(0)
        AM0(-1)  = AM0(0)
C        AP1(-1)  = AP0(0)
C        SM0    = (FM(0)-FM(-1))*(FM(0)-FM(-1))
C        AM0(-1) = 2./3./(EPS+SM0)/(EPS+SM0)
        QP1(NJ2) = QP1(NJ1)
        QP2(NJ2) = QP2(NJ1)
        QM1(NJ2) = QM1(NJ1)
        AP0(NJ2) = AP0(NJ1)
        AP1(NJ2) = AP1(NJ1)
        AM0(NJ2) = AM0(NJ1)
C
C
        DO 25 J = -1,NJ
          QM2(J) = (3.*FM(J+1)-FM(J+2))/2.
          SM1    = (FM(J+2)-FM(J+1))*(FM(J+2)-FM(J+1))
          AM1(J) = 1./3./(EPS+SM1)/(EPS+SM1)
  25    CONTINUE
        QM2(NJ1) = QM2(NJ)
        AM1(NJ1) = AM1(NJ)
        QM2(NJ2) = QM2(NJ)
        AM1(NJ2) = AM1(NJ)
C
        DO 28 J = -1,NJ2
         WP    = AP0(J)+AP1(J)
         WM    = AM0(J)+AM1(J)
         WP0   = AP0(J)/WP
         WP1   = AP1(J)/WP
         WM0   = AM0(J)/WM
         WM1   = AM1(J)/WM
         FPP   = WP0*QP1(J)+WP1*QP2(J)
         FMM   = WM0*QM1(J)+WM1*QM2(J)
         FG(J) = FPP + FMM
  28    CONTINUE
C
      DO 110 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTCFL*(FG(J) - FG(J-1))
  110 CONTINUE
   10 CONTINUE
c
      RETURN
      END
C
      SUBROUTINE  FJWEN3(C,QJ,RG,DQG)
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2,5)
      DIMENSION    C(-1:NI2,-1:NJ2,5)
C
      DIMENSION  FP(-1:MAXIJP),FM(-1:MAXIJP)
      DIMENSION  QP1(-1:MAXIJP),QP2(-1:MAXIJP),QP3(-1:MAXIJP)
      DIMENSION  QM1(-1:MAXIJP),QM2(-1:MAXIJP),QM3(-1:MAXIJP)
      DIMENSION  AP0(-1:MAXIJP),AP1(-1:MAXIJP),AP2(-1:MAXIJP)
      DIMENSION  AM0(-1:MAXIJP),AM1(-1:MAXIJP),AM2(-1:MAXIJP)
      DIMENSION    FG(-1:MAXIJP)
C
      DO 10 I = -1,NI2
C
      AC   = 1.E-6
      DO 15 J = -1,NJ2
      ACTEMP = ABS(C(I,J,1))
      AC     = MAX(AC,ACTEMP)
  15  CONTINUE
      AC  = 1.01*AC
      DO 20 J = -1,NJ2
C      AC       = ABS(C(I,J,1))
      FP(J)    = 0.5*(C(I,J,1)+AC)*RG(I,J,1)
      FM(J)    = 0.5*(C(I,J,1)-AC)*RG(I,J,1)
   20 CONTINUE
      DO 21 L =  2,5
      AC   = 1.E-6
      DO 16 J = -1,NJ2
      ACTEMP  = ABS(C(I,J,L))
      AC      = MAX(AC,ACTEMP)
   16 CONTINUE
      AC  = 1.01*AC
      DO 21 J = -1,NJ2
C      AC       = ABS(C(I,J,L))
      FP(J)    = FP(J) + 0.5*(C(I,J,L)+AC)*RG(I,J,L)
      FM(J)    = FM(J) + 0.5*(C(I,J,L)-AC)*RG(I,J,L)
21    CONTINUE
C    
        DO 24 J =  1,NJ
          QP1(J) = (2.*FP(J-2)-7.*FP(J-1)+11.*FP(J))/6.
          QP2(J) = (-FP(J-1)+5.*FP(J)+2.*FP(J+1))/6.
          QP3(J) = (2.*FP(J)+5.*FP(J+1)-FP(J+2))/6.
          QM1(J) = (-FM(J-1)+5.*FM(J)+2.*FM(J+1))/6.
          QM2(J) = (2.*FM(J)+5.*FM(J+1)-FM(J+2))/6.
          SP0    = 13./12.*(FP(J-2)-2.*FP(J-1)+FP(J))*
     1                     (FP(J-2)-2.*FP(J-1)+FP(J))+
     2                0.25*(FP(J-2)-4.*FP(J-1)+3.*FP(J))*
     3                     (FP(J-2)-4.*FP(J-1)+3.*FP(J))
          SP1    = 13./12.*(FP(J-1)-2.*FP(J)+FP(J+1))*
     1                     (FP(J-1)-2.*FP(J)+FP(J+1))+
     2                0.25*(FP(J-1)-FP(J+1))*(FP(J-1)-FP(J+1))
          SP2    = 13./12.*(FP(J+2)-2.*FP(J+1)+FP(J))*
     1                     (FP(J+2)-2.*FP(J+1)+FP(J))+
     2                0.25*(FP(J+2)-4.*FP(J+1)+3.*FP(J))*
     3                     (FP(J+2)-4.*FP(J+1)+3.*FP(J))
          SM0    = 13./12.*(FM(J-1)-2.*FM(J)+FM(J+1))*
     1                     (FM(J-1)-2.*FM(J)+FM(J+1))+
     2                0.25*(FM(J-1)-4.*FM(J)+3.*FM(J+1))*
     3                     (FM(J-1)-4.*FM(J)+3.*FM(J+1))
          SM1    = 13./12.*(FM(J)-2.*FM(J+1)+FM(J+2))*
     1                     (FM(J)-2.*FM(J+1)+FM(J+2))+
     2                0.25*(FM(J)-FM(J+2))*(FM(J)-FM(J+2))
          AP0(J) = 0.1/(EPS+SP0)/(EPS+SP0)
          AP1(J) = 0.6/(EPS+SP1)/(EPS+SP1)
          AP2(J) = 0.3/(EPS+SP2)/(EPS+SP2)
          AM0(J) = 0.3/(EPS+SM0)/(EPS+SM0)
          AM1(J) = 0.6/(EPS+SM1)/(EPS+SM1)
  24    CONTINUE
        QP1(-1)  = QP1(1)
        QP2(-1)  = QP2(1)
        QP3(-1)  = QP3(1)
        QM1(-1)  = QM1(1)
        QM2(-1)  = QM2(1)
        AP0(-1)  = AP0(1)
        AP1(-1)  = AP1(1)
        AP2(-1)  = AP2(1)
        AM0(-1)  = AM0(1)
        AM1(-1)  = AM1(1)
        QP1( 0)  = QP1(1)
        QP2( 0)  = QP2(1)
        QP3( 0)  = QP3(1)
        QM1( 0)  = QM1(1)
        QM2( 0)  = QM2(1)
        AP0( 0)  = AP0(1)
        AP1( 0)  = AP1(1)
        AP2( 0)  = AP2(1)
        AM0( 0)  = AM0(1)
        AM1( 0)  = AM1(1)
        QP1(NJ2) = QP1(NJ)
        QP2(NJ2) = QP2(NJ)
        QP3(NJ2) = QP3(NJ)
        QM1(NJ2) = QM1(NJ)
        QM2(NJ2) = QM2(NJ)
        AP0(NJ2) = AP0(NJ)
        AP1(NJ2) = AP1(NJ)
        AP2(NJ2) = AP2(NJ)
        AM0(NJ2) = AM0(NJ)
        AM1(NJ2) = AM1(NJ)
        QP1(NJ1) = QP1(NJ)
        QP2(NJ1) = QP2(NJ)
        QP3(NJ1) = QP3(NJ)
        QM1(NJ1) = QM1(NJ)
        QM2(NJ1) = QM2(NJ)
        AP0(NJ1) = AP0(NJ)
        AP1(NJ1) = AP1(NJ)
        AP2(NJ1) = AP2(NJ)
        AM0(NJ1) = AM0(NJ)
        AM1(NJ1) = AM1(NJ)
C
        DO 25 J=-1,NJ-1
          QM3(J) = (11.*FM(J+1)-7.*FM(J+2)+2.*FM(J+3))/6.
          SM2    = 13./12.*(FM(J+1)-2.*FM(J+2)+FM(J+3))*
     1                     (FM(J+1)-2.*FM(J+2)+FM(J+3))+
     2                0.25*(3.*FM(J+1)-4.*FM(J+2)+FM(J+3))*
     3                     (3.*FM(J+1)-4.*FM(J+2)+FM(J+3))
          AM2(J) = 0.1/(EPS+SM2)/(EPS+SM2)
  25    CONTINUE
        QM3(NJ ) = QM3(NJ-1)
        AM2(NJ ) = AM2(NJ-1)
        QM3(NJ1) = QM3(NJ-1)
        AM2(NJ1) = AM2(NJ-1)
        QM3(NJ2) = QM3(NJ-1)
        AM2(NJ2) = AM2(NJ-1)
C
        DO 28 J=-1,NJ2
         WP    = AP0(J)+AP1(J)+AP2(J)
         WM    = AM0(J)+AM1(J)+AM2(J)
         WP0   = AP0(J)/WP
         WP1   = AP1(J)/WP
         WP2   = AP2(J)/WP
         WM0   = AM0(J)/WM
         WM1   = AM1(J)/WM
         WM2   = AM2(J)/WM
         FPP   = WP0*QP1(J)+WP1*QP2(J)+WP2*QP3(J)
         FMM   = WM0*QM1(J)+WM1*QM2(J)+WM2*QM3(J)
         FG(J) = FPP + FMM
  28    CONTINUE
C
      DO 110 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTCFL*(FG(J) - FG(J-1))
  110 CONTINUE
   10 CONTINUE
c
      RETURN
      END
      SUBROUTINE FVIEN2(C,QJ,RG,DQG,L)
      PARAMETER (MAXIJP=1029)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2)
      DIMENSION    C(-1:NI2,-1:NJ2)
C
      DIMENSION    FG(-1:MAXIJP)
      DIMENSION   GI(-1:MAXIJP), DG(-1:MAXIJP)
      DIMENSION   FXG(-1:MAXIJP),CXG(-1:MAXIJP)
      DIMENSION   DPFG(-1:MAXIJP)
      DIMENSION  GAG(-1:MAXIJP),GJP(-1:MAXIJP), GJM(-1:MAXIJP)
      DIMENSION  SIG(-1:MAXIJP),DDG(-1:MAXIJP),XMB1(-1:MAXIJP)
      DIMENSION  ADDG(-1:MAXIJP)
      DIMENSION   AP(-1:MAXIJP), AM(-1:MAXIJP)
      DIMENSION   PSG(-1:MAXIJP)
C
      DO 10 J = -1,NJ2
C
      DO 20 I = -1,NI2
      FXG(I)   = C(I,J)*RG(I,J)
   20 CONTINUE
C     
      DO 30 I = -1,NI1
      DG(I)=2.*(RG(I+1,J)*QJ(I+1,J)-RG(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
      DPFG(I) = FXG(I+1)   - FXG(I)
   30 CONTINUE
      DG(NI2)  = DG(NI1)
      DPFG(NI2) = DPFG(NI1)
C
      DO 31 I = -1,NI1
C      IF(DG(I) .NE. 0.) THEN
C      CXG(I) = DPFG(I)/DG(I)
C      ELSE
C      CXG(I) = C(I,J) 
      CXG(I) = 0.5*(C(I+1,J) + C(I,J))
C      CXG(I) = (QJ(I+1,J)*C(I+1,J) + QJ(I,J)*C(I,J))
C     1           /(QJ(I+1,J)+QJ(I,J))
C      ENDIF
   31 CONTINUE
      CXG(NI2) = CXG(NI1)
C
      DO 33 I  = -1,NI1
      DDG(I)   = ABS(DG(I+1) - DG(I)   )
C      DDG(I)   = DG(I+1) - DG(I)   
   33 CONTINUE
      DDG(NI2) = DDG(NI1)
C
      DO 35 I  = 0,NI2
      IF(DDG(I-1) .LE. DDG(I))THEN
      XMB1(I)  = DDG(I-1)
      ELSE
      XMB1(I)  = DDG(I)
      ENDIF
   35 CONTINUE
      XMB1(-1) = XMB1(0)
C
      DO 37 I = 0,NI2
      GJP(I)  = DG (I)   - 0.5 * XMB1(I)
      GJM(I)  = DG (I-1) + 0.5 * XMB1(I-1)
   37 CONTINUE
C
      GJP(-1)  = GJP (0)
      GJM(-1)  = GJM (0)
C
      DO 40 I  = -1,NI2
      IF(GJP(I).GT.0.)THEN
        AP(I)  = 1.
      ELSE
        AP(I)  = -1.
      ENDIF
      IF(GJM(I).GT.0.)THEN
        AM(I)  = 1.
      ELSE
        AM(I)  = -1.
      ENDIF
      PSG(I) = ABS(CXG(I))
      IF(ABS(CXG(I)) .LT. EPS) PSG(I) = (CXG(I)*CXG(I)+EPSS)/EPS2
   40 CONTINUE
      DO 42 I  = -1,NI2
      GI(I)= 0.5*(AP(I)+AM(I))*MIN(AP(I)*GJP(I),AM(I)*GJM(I))
   42 CONTINUE
C
      DO 50 I = -1,NI2
C      IF(LSTD.EQ.0)THEN
      SIG(I)  = 0.5 * (PSG(I)
     1               -DTCFL*CXG(I)*CXG(I))
C      ELSE
C      SIG(I)  = 0.5 * PSG(I)
C      ENDIF
   50 CONTINUE
C
      DO 60 I = -1,NI1
      IF(DG(I).EQ.0.)THEN
        GAG(I) = 0.
      ELSE
        GAG(I) = SIG(I)*(GI(I+1) - GI(I))/DG(I)
      ENDIF
   60 CONTINUE
      GAG(NI2) = GAG(NI1)
C
      DO 62 I = -1,NI2
      CCG     = CXG(I) + GAG(I)
      PSG(I)  = ABS(CCG)
      IF(ABS(CCG).LT.EPS) PSG(I) = (CCG*CCG + EPSS)/EPS2
62    CONTINUE
c
      DO 65 I = -1,NI1
      FG(I)   = 0.5*(FXG(I)+FXG(I+1)
     &            +SIG(I)*(GI(I)+GI(I+1))-PSG(I)*DG(I))
   65 CONTINUE
C
      DO 70 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTCFL*(FG(I) - FG(I-1))
   70 CONTINUE

   10 CONTINUE
      RETURN
      END
      SUBROUTINE FVIENO(C,QJ,RG,DQ,XI,ETA,ILMT,L)
      PARAMETER (NMAX=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION    DQ(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2)
      DIMENSION    C(-1:NI2,-1:NJ2)
C
      DIMENSION     CXG(-1:NMAX),    FXG(-1:NMAX),      DPFG(-1:NMAX)
      DIMENSION     DPQ(-1:NMAX),    DMQ(-1:NMAX)
      DIMENSION    ADPQ(-1:NMAX),   ADMQ(-1:NMAX)
      DIMENSION   DPDPQ(-1:NMAX),  DMDMQ(-1:NMAX),     DPDMQ(-1:NMAX)
      DIMENSION  ADPDPQ(-1:NMAX), ADMDMQ(-1:NMAX),    ADPDMQ(-1:NMAX)
      DIMENSION    AMB1(-1:NMAX),   AMB2(-1:NMAX)
      DIMENSION      AD(-1:NMAX),     AG(-1:NMAX)
      DIMENSION      FN(-1:NMAX)
C
      IS        = -1
      IE        = NI2
      ISP1      = IS + 1
      IEM1      = IE - 1
C
      DO 10 J   = -1,NJ2
C
      DO 15 I   = IS,IE
      FXG(I)    = C(I,J)*RG(I,J)
   15 CONTINUE
C
      DO 20 I   = IS,IEM1
      DPQ (I)   = 2.*(RG(I+1,J)*QJ(I+1,J)-RG(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
      DPFG(I)   = FXG(I+1)   - FXG(I)
      ADPQ(I)   = ABS(DPQ(I))
   20 CONTINUE
      DPQ (IE) = DPQ (IEM1)
      ADPQ(IE) = ADPQ(IEM1)
      DPFG(IE) = DPFG(IEM1)
C
      DO 22 I  = IS,IEM1
      CXG(I)   = 0.5*(C(I+1,J) + C(I,J))
   22 CONTINUE
      CXG(IE) = CXG(IEM1)
C
      DO 25 I   = ISP1,IE
      DMQ (I)   = DPQ (I-1)
      ADMQ(I)   = ADPQ(I-1)
   25 CONTINUE
      DMQ (IS)  = DMQ (ISP1)
      ADMQ(IS)  = ADMQ(ISP1)
C
C     DO 30 I   = IS,IEM1
C     DPDPQ (I) = DPQ(I+1) - DPQ(I)
C     ADPDPQ(I) = ABS(DPDPQ(I))
C  30 CONTINUE
C     DPDPQ (IE)= DPDPQ (IEM1)
C     ADPDPQ(IE)= ADPDPQ(IEM1)
      DO 30 I   = IS,IEM1
      DPDPQ (I) = 2.*(RG(I+2,J)*QJ(I+2,J)-2.*RG(I+1,J)*QJ(I+1,J)
     1            + RG(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
      ADPDPQ(I) = ABS(DPDPQ(I))
   30 CONTINUE
      DPDPQ (IE)= DPDPQ (IEM1)
      ADPDPQ(IE)= ADPDPQ(IEM1)
C
      DO 35 I   = ISP1,IE
      DPDMQ (I) = DPDPQ(I-1)
      ADPDMQ(I) = ADPDPQ(I-1)
   35 CONTINUE
      DPDMQ (IS) = DPDMQ(ISP1)
      ADPDMQ(IS) = ADPDMQ(ISP1)
C
      DO 36 I   = ISP1,IE
      DMDMQ (I) = DPDMQ (I-1)
      ADMDMQ(I) = ADPDMQ(I-1)
   36 CONTINUE
      DMDMQ (IS) = DMDMQ (ISP1)
      ADMDMQ(IS) = ADMDMQ(ISP1)
C
      DO 40 I = IS,IE
      IF(ADPDMQ(I).GT.ADPDPQ(I))THEN
        AMB2(I) = DPDPQ(I)
      ELSE
        AMB2(I) = DPDMQ(I)
      ENDIF
   40 CONTINUE
C
      DO 50 I = IS,IE
      IF(ADMDMQ(I).GT.ADPDMQ(I))THEN
        AMB1(I) = DPDMQ(I)
      ELSE
        AMB1(I) = DMDMQ(I)
      ENDIF
   50 CONTINUE
C
      DO 60 I = IS,IE
      IF(ADMQ(I).GT.ADPQ(I))THEN
      AD(I) = AMB2(I)
      ELSE
      AD(I) = AMB1(I)
      ENDIF
   60 CONTINUE
C
      IF(ILMT .EQ. 1) THEN
      DO 70 I = IS,IE
      AG1 = DPQ(I) - ETA*AMB2(I)
      AG2 = DMQ(I) + ETA*AMB1(I)
      IF(AG1.GT.0.)THEN
      AG11 = 1.
      ELSE
      AG11 = -1.
      ENDIF
      IF(AG2.GT.0.) THEN
      AG22 = 1.
      ELSE
      AG22 = -1.
      ENDIF
      AG(I) = 0.5*(AG11 + AG22)*MIN(AG11*AG1,AG22*AG2)
   70 CONTINUE
      ELSEIF(ILMT .EQ. 2) THEN
      DEL2  = 1.E-6
      DO 76 I = IS,IE
      AG(I)   = (DPQ(I)*DMQ(I)+DEL2)*(DPQ(I)+DMQ(I))
     1          /(DPQ(I)*DPQ(I)+DMQ(I)*DMQ(I)+2.*DEL2)
   76 CONTINUE
      ELSEIF(ILMT .EQ. 3) THEN
      DO 78 I = IS,IE
      AG1     = 0.5*(DPQ(I)+ DMQ(I))
      IF(DPQ(I) .GT. 0.) THEN
       AG21   = 1.
      ELSE
       AG21   = -1.
      ENDIF
      IF(DMQ(I) .GT. 0.) THEN
       AG22   = 1.
      ELSE
       AG22   = -1.
      ENDIF
      AG2     = (AG21 + AG22) *MIN(AG21*DPQ(I),AG22*DMQ(I))
      IF(AG1 .GT. 0.) THEN
      AG31    = 1.
      ELSE
      AG31    = -1.
      ENDIF
      IF(AG2 .GT. 0.) THEN
      AG32    = 1.
      ELSE
      AG32    = -1.
      ENDIF
      AG(I)   = 0.5*(AG31+AG32)*MIN(AG31*AG1,AG32*AG2)
   78 CONTINUE
      ELSEIF(ILMT .EQ. 4) THEN
      DO 74 I = IS,IE
      AG1     = DPQ(I) + DMQ(I)
      IF(ABS(AG1) .GT. 1.E-6) THEN
      AG2     = DPQ(I)*DMQ(I) + ADPQ(I)*ADMQ(I)
      AG(I)   = AG2/AG1
      ELSE
      AG(I)   = 0.
      ENDIF
   74 CONTINUE
      ELSEIF(ILMT .EQ. 5) THEN
      DO 72 I = IS,IE
      S       = SIGN(ADPQ(I),DPQ(I))
C     S       = SIGN(1.,DPQ(I))
      IF(ABS(S) .GT. 1.E-6)THEN
        S=S/ABS(S)
      ENDIF
      AG1     = MIN(2.*ADPQ(I),S*DMQ(I))
      AG2     = MIN(ADPQ(I),2.*S*DMQ(I))
      AG(I)   = S*MAX(0.,AG1,AG2)
   72 CONTINUE
      ENDIF
C
      DO 80 I = IS,IEM1
      ACX   = ABS(CXG(I))
      IF(ACX.GT.EPS)THEN
        PSI = ACX
      ELSE
        PSI = (ACX*ACX + EPSS)/EPS2
      ENDIF
      IF( LSTD .EQ. 0 )THEN
        SGM   = (PSI - DTCFL*CXG(I)*CXG(I))*0.5
      ELSE
        SGM   = PSI *0.5
      ENDIF
      IF( LSTD .EQ. 0 )THEN
        ACX2  = ACX*ACX
        ACX3  = ACX2*ACX
        DTCFL2 = DTCFL*DTCFL
        SGMB1 = (DTCFL2*ACX3-ACX)/6.
        SGMB2 = (DTCFL2*ACX3-3.*DTCFL*ACX2+2.*ACX)/6.
      ELSE
        SGMB1 = -ACX/6.
        SGMB2 =  ACX/3.
      ENDIF
      IF(ADMQ(I).GT.ADPQ(I))THEN
      SGMB  = SGMB1
      ELSE
      SGMB  = SGMB2
      ENDIF
      GAM   = 0.
      DLA   = 0.
      IF(ADPQ(I).NE. 0.)THEN
      GAM   = SGM *(AG(I+1) - AG(I))/DPQ(I)
      DLA   = SGMB*(AD(I+1) - AD(I))/DPQ(I)
      ENDIF
      CXB   = CXG(I) + GAM + XI*DLA
      ACXB  = ABS(CXB)
      IF(ACXB.GT.EPS)THEN
        PSIB = ACXB
      ELSE
        PSIB = (ACXB*ACXB + EPSS)/EPS2
      ENDIF
C
      PHI    = SGM*(AG(I)+AG(I+1)) + XI*SGMB*(AD(I) + AD(I+1))
     1         - PSIB*DPQ(I)
      FN(I)  = 0.5*(FXG(I+1)+FXG(I) + PHI)
   80 CONTINUE
C
      DO 90 I = ISP1,IEM1
      DQ(I,J) = DQ(I,J) - DTCFL*(FN(I) - FN(I-1))
   90 CONTINUE
C
      I = IS
      IF(CXG(I).LT.0.)THEN
         DQ(I,J)    = DQ(I,J) - DTCFL*(FXG(I+1)-FXG(I))
      ENDIF
C
      I = IE
      IF(CXG(I).GT.0.)THEN
         DQ(I,J)      = DQ(I,J)  - DTCFL*(FXG(I)-FXG(I-1))
      ENDIF
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE FVITVD2(C,QJ,RG,DQG,L)
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2)
      DIMENSION    C(-1:NI2,-1:NJ2)
C
      DIMENSION  GI(-1:MAXIJP),DG(-1:MAXIJP)
      DIMENSION  ADG(-1:MAXIJP)
      DIMENSION  AS(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP)
      DIMENSION  FXG(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  PSG(-1:MAXIJP)
      DIMENSION  SIG(-1:MAXIJP)
      DIMENSION  GAG(-1:MAXIJP)
      DIMENSION  AGAG(-1:MAXIJP)
      DIMENSION    FG(-1:MAXIJP)
C
      DO 10 J = -1,NJ2
C
      DO 20 I = -1,NI2
      FXG(I)   = C(I,J)*RG(I,J)
   20 CONTINUE
C     
      DO 30 I = -1,NI1
      DG(I)=2.*(RG(I+1,J)*QJ(I+1,J)-RG(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
      DPFG(I) = FXG(I+1)   - FXG(I)
   30 CONTINUE
      DG(NI2)  = DG(NI1)
      DPFG(NI2) = DPFG(NI1)
C
      DO 40 I = -1,NI1
c      IF(DG(I) .NE. 0.) THEN
c      CXG(I) = DPFG(I)/DG(I)
c      ELSE
        CXG(I) = 0.5*(C(I+1,J) + C(I,J))
c      ENDIF
   40 CONTINUE
      CXG(NI2) = CXG(NI1)
C
      DO 50 I = -1,NI2
      IF(DG(I).GE.0.)THEN
      AS(I) = 1.
      ELSE
      AS(I) =-1.
      ENDIF
      ADG(I) = ABS(DG(I))
      ACXG(I) = ABS(CXG(I))
      PSG(I)  = ACXG(I)
      IF(ACXG(I) .LT. EPS) PSG(I) = (CXG(I)*CXG(I)+EPSS)/EPS2
   50 CONTINUE
      DO 60 I = 0,NI2
      GI(I) = 0.5*(AS(I-1)+AS(I))*MIN(ADG(I-1),ADG(I))
   60 CONTINUE
      GI( -1) = GI(  0)
C
      DO 70 I = -1,NI2
      IF(LSTD.EQ.0)THEN
       SIG(I) = 0.5*(PSG(I)
     1              - DTCFL*CXG(I)*CXG(I))
       ELSE
       SIG(I) = 0.5*PSG(I)
       ENDIF
   70 CONTINUE
C
      DO 80 I = -1,NI1
      IF(DG(I).NE.0.)THEN
      GAG(I)  = SIG(I)*(GI(I+1)-GI(I))/DG(I)
      ELSE
      GAG(I)  = 0.
      ENDIF
   80 CONTINUE
      GAG(NI2) = GAG(NI1)
C
      DO 90 I = -1,NI2
      AGAG(I)  = ABS(CXG(I) + GAG(I))
      PSG(I)   = AGAG(I)
      IF(AGAG(I) .LT. EPS) PSG(I) = (AGAG(I)*AGAG(I)+EPSS)/EPS2
   90 CONTINUE
C
      DO 100 I = -1,NI1
      FG(I)  = 0.5*(FXG(I+1)+FXG(I)+SIG(I)*(GI(I+1)+GI(I))
     1              -PSG(I)*DG(I))
  100 CONTINUE
      DO 110 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTCFL*(FG(I) - FG(I-1))
  110 CONTINUE
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FVJEN2(C,QJ,RG,DQG,L)
      PARAMETER (MAXIJP=1029)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2)
      DIMENSION    C(-1:NI2,-1:NJ2)
C
      DIMENSION    FG(-1:MAXIJP)
      DIMENSION   GI(-1:MAXIJP), DG(-1:MAXIJP)
      DIMENSION   FXG(-1:MAXIJP),CXG(-1:MAXIJP)
      DIMENSION   DPFG(-1:MAXIJP)
      DIMENSION  GAG(-1:MAXIJP),GJP(-1:MAXIJP), GJM(-1:MAXIJP)
      DIMENSION  SIG(-1:MAXIJP),DDG(-1:MAXIJP),XMB1(-1:MAXIJP)
      DIMENSION  ADDG(-1:MAXIJP)
      DIMENSION   AP(-1:MAXIJP), AM(-1:MAXIJP)
      DIMENSION   PSG(-1:MAXIJP)
C
      DO 10 I = -1,NI2
C
      DO 20 J = -1,NJ2
      FXG(J)   = C(I,J)*RG(I,J)
   20 CONTINUE
C     
      DO 30 J = -1,NJ1
      DG(J)=2.*(RG(I,J+1)*QJ(I,J+1)-RG(I,J)*QJ(I,J))
     1           /(QJ(I,J+1)+QJ(I,J))
      DPFG(J) = FXG(J+1)   - FXG(J)
   30 CONTINUE
      DG(NJ2)  = DG(NJ1)
      DPFG(NJ2) = DPFG(NJ1)
C
      DO 31 J = -1,NJ1
C      IF(DG(J) .NE. 0.) THEN
C      CXG(J) = DPFG(J)/DG(J)
C      ELSE
      CXG(J) = 0.5*(C(I,J+1) + C(I,J))
C      CXG(J) = (QJ(I,J+1)*C(I,J+1) + QJ(I,J)*C(I,J))
C     1           /(QJ(I,J+1)+QJ(I,J))
C      ENDIF
      CXG(J) =  C(I,J)
   31 CONTINUE
      CXG(NJ2) = CXG(NJ1)
C
      DO 33 J  = -1,NJ1
      DDG(J)   = ABS(DG(J+1) - DG(J)   )
C      DDG(J)   = DG(J+1) - DG(J)   
   33 CONTINUE
      DDG(NJ2) = DDG(NJ1)
C
      DO 35 J  = 0,NJ2
      IF(DDG(J-1) .LE. DDG(J))THEN
      XMB1(J)  = DDG(J-1)
      ELSE
      XMB1(J)  = DDG(J)
      ENDIF
   35 CONTINUE
      XMB1(-1) = XMB1(0)
C
      DO 37 J = 0,NJ2
      GJP(J)  = DG (J)   - 0.5 * XMB1(J)
      GJM(J)  = DG (J-1) + 0.5 * XMB1(J-1)
   37 CONTINUE
C
      GJP(-1)  = GJP (0)
      GJM(-1)  = GJM (0)
C
      DO 40 J  = -1,NJ2
      IF(GJP(J).GT.0.)THEN
        AP(J)  = 1.
      ELSE
        AP(J)  = -1.
      ENDIF
      IF(GJM(J).GT.0.)THEN
        AM(J)  = 1.
      ELSE
        AM(J)  = -1.
      ENDIF
      PSG(J) = ABS(CXG(J))
      IF(ABS(CXG(J)) .LT. EPS) PSG(J) = (CXG(J)*CXG(J)+EPSS)/EPS2
   40 CONTINUE
      DO 42 J  = -1,NJ2
      GI(J)= 0.5*(AP(J)+AM(J))*MIN(AP(J)*GJP(J),AM(J)*GJM(J))
   42 CONTINUE
C
      DO 50 J = -1,NJ2
C      IF(LSTD.EQ.0)THEN
      SIG(J)  = 0.5 * (PSG(J)
     1               -DTCFL*CXG(J)*CXG(J))
C      ELSE
C      SIG(J)  = 0.5 * PSG(J)
C      ENDIF
   50 CONTINUE
C
      DO 60 J = -1,NJ1
      IF(DG(J).EQ.0.)THEN
        GAG(J) = 0.
      ELSE
        GAG(J) = SIG(J)*(GI(J+1) - GI(J))/DG(J)
      ENDIF
   60 CONTINUE
      GAG(NJ2) = GAG(NJ1)
C
      DO 62 J = -1,NJ2
      CCG     = CXG(J) + GAG(J)
      PSG(J)  = ABS(CCG)
      IF(ABS(CCG).LT.EPS) PSG(J) = (CCG*CCG + EPSS)/EPS2
62    CONTINUE
c
      DO 65 J = -1,NJ1
      FG(J)   = 0.5*(FXG(J)+FXG(J+1)
     &            +SIG(J)*(GI(J)+GI(J+1))-PSG(J)*DG(J))
   65 CONTINUE
C
      DO 70 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTCFL*(FG(J) - FG(J-1))
   70 CONTINUE

   10 CONTINUE
      RETURN
      END
      SUBROUTINE FVJENO(C,QJ,RG,DQ,XI,ETA,ILMT,L)
      PARAMETER (NMAX=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION    DQ(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2)
      DIMENSION    C(-1:NI2,-1:NJ2)
C
      DIMENSION     CXG(-1:NMAX),    FXG(-1:NMAX),      DPFG(-1:NMAX)
      DIMENSION     DPQ(-1:NMAX),    DMQ(-1:NMAX)
      DIMENSION    ADPQ(-1:NMAX),   ADMQ(-1:NMAX)
      DIMENSION   DPDPQ(-1:NMAX),  DMDMQ(-1:NMAX),     DPDMQ(-1:NMAX)
      DIMENSION  ADPDPQ(-1:NMAX), ADMDMQ(-1:NMAX),    ADPDMQ(-1:NMAX)
      DIMENSION    AMB1(-1:NMAX),   AMB2(-1:NMAX)
      DIMENSION      AD(-1:NMAX),     AG(-1:NMAX)
      DIMENSION      FN(-1:NMAX)
C
      JS        = -1
      JE        = NJ2
      JSP1      = JS + 1
      JEM1      = JE - 1
C
      DO 10 I   = -1,NI2
C
      DO 15 J   = JS,JE
      FXG(J)    = C(I,J)*RG(I,J)
   15 CONTINUE
C
      DO 20 J   = JS,JEM1
      DPQ (J)   = 2.*(RG(I,J+1)*QJ(I,J+1)-RG(I,J)*QJ(I,J))
     1           /(QJ(I,J+1)+QJ(I,J))
      DPFG(J)   = FXG(J+1)   - FXG(J)
      ADPQ(J)   = ABS(DPQ(J))
   20 CONTINUE
      DPQ (JE) = DPQ (JEM1)
      ADPQ(JE) = ADPQ(JEM1)
      DPFG(JE) = DPFG(JEM1)
C
      DO 22 J  = JS,JEM1
      CXG(J)   = 0.5*(C(I,J+1) + C(I,J))
   22 CONTINUE
      CXG(JE) = CXG(JEM1)
C
      DO 25 J   = JSP1,JE
      DMQ (J)   = DPQ (J-1)
      ADMQ(J)   = ADPQ(J-1)
   25 CONTINUE
      DMQ (JS)  = DMQ (JSP1)
      ADMQ(JS)  = ADMQ(JSP1)
C
C     DO 30 J   = JS,JEM1
C     DPDPQ (J) = DPQ(J+1) - DPQ(J)
C     ADPDPQ(J) = ABS(DPDPQ(J))
C  30 CONTINUE
C     DPDPQ (JE)= DPDPQ (JEM1)
C     ADPDPQ(JE)= ADPDPQ(JEM1)
C
      DO 30 J   = JS,JEM1
      DPDPQ(J)  = 2.*(RG(I,J+2)*QJ(I,J+2)-2.*RG(I,J+1)*QJ(I,J+1)
     1            +RG(I,J)*QJ(I,J))
     1           /(QJ(I,J+1)+QJ(I,J))
      ADPDPQ(J) = ABS(DPDPQ(J))
   30 CONTINUE
      DPDPQ (JE)= DPDPQ (JEM1)
      ADPDPQ(JE)= ADPDPQ(JEM1)
C
      DO 35 J   = JSP1,JE
      DPDMQ (J) = DPDPQ(J-1)
      ADPDMQ(J) = ADPDPQ(J-1)
   35 CONTINUE
      DPDMQ (JS) = DPDMQ(JSP1)
      ADPDMQ(JS) = ADPDMQ(JSP1)
C
      DO 36 J   = JSP1,JE
      DMDMQ (J) = DPDMQ (J-1)
      ADMDMQ(J) = ADPDMQ(J-1)
   36 CONTINUE
      DMDMQ (JS) = DMDMQ (JSP1)
      ADMDMQ(JS) = ADMDMQ(JSP1)
C
      DO 40 J = JS,JE
      IF(ADPDMQ(J).GT.ADPDPQ(J))THEN
        AMB2(J) = DPDPQ(J)
      ELSE
        AMB2(J) = DPDMQ(J)
      ENDIF
   40 CONTINUE
C
      DO 50 J = JS,JE
      IF(ADMDMQ(J).GT.ADPDMQ(J))THEN
        AMB1(J) = DPDMQ(J)
      ELSE
        AMB1(J) = DMDMQ(J)
      ENDIF
   50 CONTINUE
C
      DO 60 J = JS,JE
      IF(ADMQ(J).GT.ADPQ(J))THEN
      AD(J) = AMB2(J)
      ELSE
      AD(J) = AMB1(J)
      ENDIF
   60 CONTINUE
C
      IF(ILMT .EQ. 1) THEN
      DO 70 J = JS,JE
      AG1 = DPQ(J) - ETA*AMB2(J)
      AG2 = DMQ(J) + ETA*AMB1(J)
      IF(AG1.GT.0.)THEN
      AG11 = 1.
      ELSE
      AG11 = -1.
      ENDIF
      IF(AG2.GT.0.) THEN
      AG22 = 1.
      ELSE
      AG22 = -1.
      ENDIF
      AG(J) = 0.5*(AG11 + AG22)*MIN(AG11*AG1,AG22*AG2)
   70 CONTINUE
      ELSEIF(ILMT .EQ. 2) THEN
      DEL2  = 1.E-6
      DO 76 J = JS,JE
      AG(J)   = (DPQ(J)*DMQ(J)+DEL2)*(DPQ(J)+DMQ(J))
     1          /(DPQ(J)*DPQ(J)+DMQ(J)*DMQ(J)+2.*DEL2)
   76 CONTINUE
      ELSEIF(ILMT .EQ. 3) THEN
      DO 78 J = JS,JE
      AG1     = 0.5*(DPQ(J)+ DMQ(J))
      IF(DPQ(J) .GT. 0.) THEN
       AG21   = 1.
      ELSE
       AG21   = -1.
      ENDIF
      IF(DMQ(J) .GT. 0.) THEN
       AG22   = 1.
      ELSE
       AG22   = -1.
      ENDIF
      AG2     = (AG21 + AG22) *MIN(AG21*DPQ(J),AG22*DMQ(J))
      IF(AG1 .GT. 0.) THEN
      AG31    = 1.
      ELSE
      AG31    = -1.
      ENDIF
      IF(AG2 .GT. 0.) THEN
      AG32    = 1.
      ELSE
      AG32    = -1.
      ENDIF
      AG(J)   = 0.5*(AG31+AG32)*MIN(AG31*AG1,AG32*AG2)
   78 CONTINUE
      ELSEIF(ILMT .EQ. 4) THEN
      DO 74 J = JS,JE
      AG1     = DPQ(J) + DMQ(J)
      IF(ABS(AG1) .GT. 1.E-6) THEN
      AG2     = DPQ(J)*DMQ(J) + ADPQ(J)*ADMQ(J)
      AG(J)   = AG2/AG1
      ELSE
      AG(J)   = 0.
      ENDIF
   74 CONTINUE
      ELSEIF(ILMT .EQ. 5) THEN
      DO 72 J = JS,JE
      S       = SIGN(ADPQ(J),DPQ(J))
C     S       = SIGN(1.,DPQ(J))
      IF(ABS(S) .GT. 1.E-6)THEN
        S=S/ABS(S)
      ENDIF
      AG1     = MIN(2.*ADPQ(J),S*DMQ(J))
      AG2     = MIN(ADPQ(J),2.*S*DMQ(J))
      AG(J)   = S*MAX(0.,AG1,AG2)
   72 CONTINUE
      ENDIF
C
      DO 80 J = JS,JEM1
      ACX   = ABS(CXG(J))
      IF(ACX.GT.EPS)THEN
        PSI = ACX
      ELSE
        PSI = (ACX*ACX + EPSS)/EPS2
      ENDIF
      IF( LSTD .EQ. 0 )THEN
        SGM   = (PSI - DTCFL*CXG(J)*CXG(J))*0.5
      ELSE
        SGM   = PSI *0.5
      ENDIF
      IF( LSTD .EQ. 0 )THEN
        ACX2  = ACX*ACX
        ACX3  = ACX2*ACX
        DTCFL2 = DTCFL*DTCFL
        SGMB1 = (DTCFL2*ACX3-ACX)/6.
        SGMB2 = (DTCFL2*ACX3-3.*DTCFL*ACX2+2.*ACX)/6.
      ELSE
        SGMB1 = -ACX/6.
        SGMB2 =  ACX/3.
      ENDIF
      IF(ADMQ(J).GT.ADPQ(J))THEN
      SGMB  = SGMB1
      ELSE
      SGMB  = SGMB2
      ENDIF
      GAM   = 0.
      DLA   = 0.
      IF(ADPQ(J).NE.0.)THEN
      GAM   = SGM *(AG(J+1) - AG(J))/DPQ(J)
      DLA   = SGMB*(AD(J+1) - AD(J))/DPQ(J)
      ENDIF
      CXB   = CXG(J) + GAM + XI*DLA
      ACXB  = ABS(CXB)
      IF(ACXB.GT.EPS)THEN
        PSIB = ACXB
      ELSE
        PSIB = (ACXB*ACXB + EPSS)/EPS2
      ENDIF
C
      PHI    = SGM*(AG(J)+AG(J+1)) + XI*SGMB*(AD(J) + AD(J+1))
     1         - PSIB*DPQ(J)
      FN(J)  = 0.5*(FXG(J+1)+FXG(J) + PHI)
   80 CONTINUE
C
      DO 90 J = JSP1,JEM1
      DQ(I,J) = DQ(I,J) - DTCFL*(FN(J) - FN(J-1))
   90 CONTINUE
C
      J = JS
      IF(CXG(J).LT.0.)THEN
         DQ(I,J)    = DQ(I,J) - DTCFL*(FXG(J+1)-FXG(J))
      ENDIF
C
      J = JE
      IF(CXG(J).GT.0.)THEN
         DQ(I,J)      = DQ(I,J)  - DTCFL*(FXG(J)-FXG(J-1))
      ENDIF
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE FVJTVD2(C,QJ,RG,DQG,L)
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
C
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2)
      DIMENSION    C(-1:NI2,-1:NJ2)
C
      DIMENSION  GI(-1:MAXIJP),DG(-1:MAXIJP)
      DIMENSION  ADG(-1:MAXIJP)
      DIMENSION  AS(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP)
      DIMENSION  FXG(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  PSG(-1:MAXIJP)
      DIMENSION  SIG(-1:MAXIJP)
      DIMENSION  GAG(-1:MAXIJP)
      DIMENSION  AGAG(-1:MAXIJP)
      DIMENSION    FG(-1:MAXIJP)
C
      DO 10 I = -1,NI2
C
      DO 20 J = -1,NJ2
      FXG(J)   = C(I,J)*RG(I,J)
   20 CONTINUE
C     
      DO 30 J = -1,NJ1
      DG(J)=2.*(RG(I,J+1)*QJ(I,J+1)-RG(I,J)*QJ(I,J))
     1           /(QJ(I,J+1)+QJ(I,J))
      DPFG(J) = FXG(J+1)   - FXG(J)
   30 CONTINUE
      DG(NJ2)  = DG(NJ1)
      DPFG(NJ2) = DPFG(NJ1)
C
      DO 40 J = -1,NJ1
c      IF(DG(J) .NE. 0.) THEN
c      CXG(J) = DPFG(J)/DG(J)
c      ELSE
        CXG(J) = 0.5*(C(I,J+1) + C(I,J))
c      ENDIF
   40 CONTINUE
      CXG(NJ2) = CXG(NJ1)
C
      DO 50 J = -1,NJ2
      IF(DG(J).GE.0.)THEN
      AS(J) = 1.
      ELSE
      AS(J) =-1.
      ENDIF
      ADG(J) = ABS(DG(J))
      ACXG(J) = ABS(CXG(J))
      PSG(J)  = ACXG(J)
      IF(ACXG(J) .LT. EPS) PSG(J) = (CXG(J)*CXG(J)+EPSS)/EPS2
   50 CONTINUE
      DO 60 J = 0,NJ2
      GI(J) = 0.5*(AS(J-1)+AS(J))*MIN(ADG(J-1),ADG(J))
   60 CONTINUE
      GI( -1) = GI(  0)
C
      DO 70 J = -1,NJ2
      IF(LSTD.EQ.0)THEN
       SIG(J) = 0.5*(PSG(J)
     1              - DTCFL*CXG(J)*CXG(J))
       ELSE
       SIG(J) = 0.5*PSG(J)
       ENDIF
   70 CONTINUE
C
      DO 80 J = -1,NJ1
      IF(DG(J).NE.0.)THEN
      GAG(J)  = SIG(J)*(GI(J+1)-GI(J))/DG(J)
      ELSE
      GAG(J)  = 0.
      ENDIF
   80 CONTINUE
      GAG(NJ2) = GAG(NJ1)
C
      DO 90 J = -1,NJ2
      AGAG(J)  = ABS(CXG(J) + GAG(J))
      PSG(J)   = AGAG(J)
      IF(AGAG(J) .LT. EPS) PSG(J) = (AGAG(J)*AGAG(J)+EPSS)/EPS2
   90 CONTINUE
C
      DO 100 J = -1,NJ1
      FG(J)  = 0.5*(FXG(J+1)+FXG(J)+SIG(J)*(GI(J+1)+GI(J))
     1              -PSG(J)*DG(J))
  100 CONTINUE
      DO 110 J = 2,NJ1
      DQG(I,J) = DQG(I,J) - DTCFL*(FG(J) - FG(J-1))
  110 CONTINUE
   10 CONTINUE
      RETURN
      END
      SUBROUTINE GRIDIN( X , Y ,IJP)
      COMMON /TAPE / INAME,IGRID,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /NEXT / NEXTLOC
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /BC01 / IDBC(30,20),FDBC(10,20)
C
      DIMENSION X(-1:NI2,-1:NJ2),Y(-1:NI2,-1:NJ2)
      DIMENSION IJP(MXBP)
C
      READ(IGRID,*) ((X(I,J),I=1,NI),J=1,NJ),
     1              ((Y(I,J),I=1,NI),J=1,NJ)
C
      DO 100 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IDA    = IDBC( 2,IBC)
      ISA    = IDBC( 3,IBC)
      ICA    = IDBC( 4,IBC)
      I2S    = IDBC( 5,IBC)
      I2E    = IDBC( 6,IBC)
      IPA    = IDBC( 7,IBC)
      LENA   = (I2E - I2S)*ICA + 1
      IF(IBTYPE .NE. 4) THEN
      DO 110 I2 = 1,LENA
      IAI    = IJP(IPA + LENA*0 + I2)
      IAJ    = IJP(IPA + LENA*1 + I2)
      IBI    = IJP(IPA + LENA*2 + I2)
      IBJ    = IJP(IPA + LENA*3 + I2)
      ICI    = IJP(IPA + LENA*4 + I2)
      ICJ    = IJP(IPA + LENA*5 + I2)
      IDI    = IJP(IPA + LENA*6 + I2)
      IDJ    = IJP(IPA + LENA*7 + I2)
      IEI    = IJP(IPA + LENA*8 + I2)
      IEJ    = IJP(IPA + LENA*9 + I2)
      IF(IBTYPE .EQ. 0 .OR. IBTYPE .EQ. 1 .OR. IBTYPE .EQ. 5)THEN
      X(IBI,IBJ) = 2.*X(ICI,ICJ) - X(IDI,IDJ)
      Y(IBI,IBJ) = 2.*Y(ICI,ICJ) - Y(IDI,IDJ)
      X(IAI,IAJ) = 2.*X(IBI,IBJ) - X(ICI,ICJ)
      Y(IAI,IAJ) = 2.*Y(IBI,IBJ) - Y(ICI,ICJ)
      ENDIF
      IF(IBTYPE .EQ. 2 .OR. IBTYPE .EQ. 3)THEN
      X(IAI,IAJ) =   X(IEI,IEJ)
      Y(IAI,IAJ) = - Y(IEI,IEJ)
      X(IBI,IBJ) =   X(IDI,IDJ)
      Y(IBI,IBJ) = - Y(IDI,IDJ)
      Y(ICI,ICJ) = 0.
      ENDIF
  110 CONTINUE
      IF( IBTYPE .EQ. 7 )THEN
      IDB    = IDBC(12,IBC)
      ISB    = IDBC(13,IBC)
      ICB    = IDBC(14,IBC)
      I3S    = IDBC(15,IBC)
      I3E    = IDBC(16,IBC)
      IPB    = IDBC( 8,IBC)
      LENB   = (I3E - I3S)*ICB + 1
      DO 105 I2 = 1,LENA
      IAI    = IJP(IPA + LENA*0 + I2)
      IAJ    = IJP(IPA + LENA*1 + I2)
      IBI    = IJP(IPA + LENA*2 + I2)
      IBJ    = IJP(IPA + LENA*3 + I2)
      ICI    = IJP(IPA + LENA*4 + I2)
      ICJ    = IJP(IPA + LENA*5 + I2)
      IDI    = IJP(IPA + LENA*6 + I2)
      IDJ    = IJP(IPA + LENA*7 + I2)
      IEI    = IJP(IPA + LENA*8 + I2)
      IEJ    = IJP(IPA + LENA*9 + I2)
C
      JAI    = IJP(IPB + LENB*0 + I2)
      JAJ    = IJP(IPB + LENB*1 + I2)
      JBI    = IJP(IPB + LENB*2 + I2)
      JBJ    = IJP(IPB + LENB*3 + I2)
      JCI    = IJP(IPB + LENB*4 + I2)
      JCJ    = IJP(IPB + LENB*5 + I2)
      JDI    = IJP(IPB + LENB*6 + I2)
      JDJ    = IJP(IPB + LENB*7 + I2)
      JEI    = IJP(IPB + LENB*8 + I2)
      JEJ    = IJP(IPB + LENB*9 + I2)
      WRITE(9,72)IAI,IAJ,IBI,IBJ,ICI,ICJ,IDI,IDJ,IEI,IEJ,
     1            JAI,JAJ,JBI,JBJ,JCI,JCJ,JDI,JDJ,JEI,JEJ
   72 FORMAT(20I4)
      X(IBI,IBJ) = 2.*X(ICI,ICJ) - X(IDI,IDJ)
      Y(IBI,IBJ) = 2.*Y(ICI,ICJ) - Y(IDI,IDJ)
      X(IAI,IAJ) = 2.*X(IBI,IBJ) - X(ICI,ICJ)
      Y(IAI,IAJ) = 2.*Y(IBI,IBJ) - Y(ICI,ICJ)
      X(JBI,JBJ) = 2.*X(JCI,JCJ) - X(JDI,JDJ)
      Y(JBI,JBJ) = 2.*Y(JCI,JCJ) - Y(JDI,JDJ)
      X(JAI,JAJ) = 2.*X(JBI,JBJ) - X(JCI,JCJ)
      Y(JAI,JAJ) = 2.*Y(JBI,JBJ) - Y(JCI,JCJ)
  105 CONTINUE
      ENDIF
      ELSE
      IDB    = IDBC(12,IBC)
      ISB    = IDBC(13,IBC)
      ICB    = IDBC(14,IBC)
      I3S    = IDBC(15,IBC)
      I3E    = IDBC(16,IBC)
      IPB    = IDBC( 8,IBC)
      LENB   = (I3E - I3S)*ICB + 1
      DO 120 I2 = 1,LENA
      IAI    = IJP(IPA + LENA*0 + I2)
      IAJ    = IJP(IPA + LENA*1 + I2)
      IBI    = IJP(IPA + LENA*2 + I2)
      IBJ    = IJP(IPA + LENA*3 + I2)
      ICI    = IJP(IPA + LENA*4 + I2)
      ICJ    = IJP(IPA + LENA*5 + I2)
      IDI    = IJP(IPA + LENA*6 + I2)
      IDJ    = IJP(IPA + LENA*7 + I2)
      IEI    = IJP(IPA + LENA*8 + I2)
      IEJ    = IJP(IPA + LENA*9 + I2)
C
      JAI    = IJP(IPB + LENB*0 + I2)
      JAJ    = IJP(IPB + LENB*1 + I2)
      JBI    = IJP(IPB + LENB*2 + I2)
      JBJ    = IJP(IPB + LENB*3 + I2)
      JCI    = IJP(IPB + LENB*4 + I2)
      JCJ    = IJP(IPB + LENB*5 + I2)
      JDI    = IJP(IPB + LENB*6 + I2)
      JDJ    = IJP(IPB + LENB*7 + I2)
      JEI    = IJP(IPB + LENB*8 + I2)
      JEJ    = IJP(IPB + LENB*9 + I2)
      WRITE(9,71)IAI,IAJ,IBI,IBJ,ICI,ICJ,IDI,IDJ,IEI,IEJ,
     1            JAI,JAJ,JBI,JBJ,JCI,JCJ,JDI,JDJ,JEI,JEJ
   71 FORMAT(20I4)
      X(IAI,IAJ) =   X(JEI,JEJ)
      Y(IAI,IAJ) =   Y(JEI,JEJ)
      X(IBI,IBJ) =   X(JDI,JDJ)
      Y(IBI,IBJ) =   Y(JDI,JDJ)
      X(JAI,JAJ) =   X(IEI,IEJ)
      Y(JAI,JAJ) =   Y(IEI,IEJ)
      X(JBI,JBJ) =   X(IDI,IDJ)
      Y(JBI,JBJ) =   Y(IDI,IDJ)
  120 CONTINUE
      ENDIF
  100 CONTINUE
      I = NI1
      J = NJ1
      XM1        = 0.5*(X(I-1,J  ) + X(I  ,J-1))
      YM1        = 0.5*(Y(I-1,J  ) + Y(I  ,J-1))
      X(I  ,J  ) =   2.*XM1        - X(I-1,J-1)
      Y(I  ,J  ) =   2.*YM1        - Y(I-1,J-1)
      XM2        = 0.5*(X(I-1,J+1) + X(I+1,J-1))
      YM2        = 0.5*(Y(I-1,J+1) + Y(I+1,J-1))
      X(I+1,J+1) =   2.*XM2        - X(I-1,J-1)
      Y(I+1,J+1) =   2.*YM2        - Y(I-1,J-1)
      X(I  ,J+1) = 0.5*(X(I+1,J+1) + X(I-1,J+1))
      Y(I  ,J+1) = 0.5*(Y(I+1,J+1) + Y(I-1,J+1))
      X(I+1,J  ) = 0.5*(X(I+1,J-1) + X(I+1,J+1))
      Y(I+1,J  ) = 0.5*(Y(I+1,J-1) + Y(I+1,J+1))
      I = 0
      J = NJ1
      XM1        = 0.5*(X(I+1,J  ) + X(I  ,J-1))
      YM1        = 0.5*(Y(I+1,J  ) + Y(I  ,J-1))
      X(I  ,J  ) =   2.*XM1        - X(I+1,J-1)
      Y(I  ,J  ) =   2.*YM1        - Y(I+1,J-1)
      XM2        = 0.5*(X(I-1,J-1) + X(I+1,J+1))
      YM2        = 0.5*(Y(I-1,J-1) + Y(I+1,J+1))
      X(I-1,J+1) =   2.*XM2        - X(I+1,J-1)
      Y(I-1,J+1) =   2.*YM2        - Y(I+1,J-1)
      X(I  ,J+1) = 0.5*(X(I+1,J+1) + X(I-1,J+1))
      Y(I  ,J+1) = 0.5*(Y(I+1,J+1) + Y(I-1,J+1))
      X(I-1,J  ) = 0.5*(X(I-1,J-1) + X(I-1,J+1))
      Y(I-1,J  ) = 0.5*(Y(I-1,J-1) + Y(I-1,J+1))
      I = NI1
      J = 0
      XM1        = 0.5*(X(I-1,J  ) + X(I  ,J+1))
      YM1        = 0.5*(Y(I-1,J  ) + Y(I  ,J+1))
      X(I  ,J  ) =   2.*XM1        - X(I-1,J+1)
      Y(I  ,J  ) =   2.*YM1        - Y(I-1,J+1)
      XM2        = 0.5*(X(I-1,J-1) + X(I+1,J+1))
      YM2        = 0.5*(Y(I-1,J-1) + Y(I+1,J+1))
      X(I+1,J-1) =   2.*XM2        - X(I-1,J+1)
      Y(I+1,J-1) =   2.*YM2        - Y(I-1,J+1)
      X(I  ,J-1) = 0.5*(X(I+1,J-1) + X(I-1,J-1))
      Y(I  ,J-1) = 0.5*(Y(I+1,J-1) + Y(I-1,J-1))
      X(I+1,J  ) = 0.5*(X(I+1,J-1) + X(I+1,J+1))
      Y(I+1,J  ) = 0.5*(Y(I+1,J-1) + Y(I+1,J+1))
      I = 0
      J = 0
      XM1        = 0.5*(X(I+1,J  ) + X(I  ,J+1))
      YM1        = 0.5*(Y(I+1,J  ) + Y(I  ,J+1))
      X(I  ,J  ) =   2.*XM1        - X(I+1,J+1)
      Y(I  ,J  ) =   2.*YM1        - Y(I+1,J+1)
      XM2        = 0.5*(X(I+1,J-1) + X(I-1,J+1))
      YM2        = 0.5*(Y(I+1,J-1) + Y(I-1,J+1))
      X(I-1,J-1) =   2.*XM2        - X(I+1,J+1)
      Y(I-1,J-1) =   2.*YM2        - Y(I+1,J+1)
      X(I  ,J-1) = 0.5*(X(I+1,J-1) + X(I-1,J-1))
      Y(I  ,J-1) = 0.5*(Y(I+1,J-1) + Y(I-1,J-1))
      X(I-1,J  ) = 0.5*(X(I-1,J-1) + X(I-1,J+1))
      Y(I-1,J  ) = 0.5*(Y(I-1,J-1) + Y(I-1,J+1))
c     WRITE(1,*)MI,MJ
c     WRITE(1,*)((X(I,J),I=-1,MI-2),J=-1,MJ-2),
c    1           ((Y(I,J),I=-1,MI-2),J=-1,MJ-2)
c     CLOSE(1)
      RETURN
      END
      SUBROUTINE INITQ(X,Y,R,T,U,V,P,E,PSI)
      COMMON /TAPE / INAME,IGRID,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /CNST2/ R0,U0,V0,P0,T0,E0
      COMMON /CNST4/ DR0,UR0,VR0,PR0,TR0,ER0,
     1               DL0,UL0,VL0,PL0,TL0,EL0
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /PSII/ ICYL,PSI0,AMPLI,NCURVE,URATIO,YY0
      DIMENSION    U(-1:NI2,-1:NJ2),   V(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2),   T(-1:NI2,-1:NJ2)
      DIMENSION    P(-1:NI2,-1:NJ2),   E(-1:NI2,-1:NJ2)
      DIMENSION    X(-1:NI2,-1:NJ2),   Y(-1:NI2,-1:NJ2)
      DIMENSION    PSI(-1:NI2,-1:NJ2)
C
      IF( ICYL .EQ. 1 .AND. LSTD .EQ. 0 )GOTO 100
      IF( LSTD .EQ. 1 )THEN
      DO 10 IJ  =-1,MIJ2
       R(IJ,-1) = R0
       U(IJ,-1) = U0
       V(IJ,-1) = V0
       E(IJ,-1) = E0
       T(IJ,-1) = T0
       P(IJ,-1) = P0
   10 CONTINUE
      ELSE
      X1  = PSI0
      X2  =-X1
      DO 40 I = -1,NI2
        IF( NCURVE .EQ. 1 )THEN
c        YY = 1.0 + AMPLI*(-0.5+SIN(PI*X(I,(NJ+1)/2)))
         YY = YY0 + AMPLI*COS(PI*X(I,(NJ+1)/2))
        ELSEIF( NCURVE .EQ. 2 )THEN
c        YY = 0.5 + AMPLI*COS(4.*PI*X(I,(NJ+1)/2))
         YY = YY0 + AMPLI*COS(4.*PI*X(I,(NJ+1)/2))
        ELSEIF( NCURVE .EQ. 3 )THEN
         YYT = YY0 - AMPLI*SIN(PI*X(I,(NJ+1)/2))
         YYB =-YY0 + AMPLI*SIN(PI*X(I,(NJ+1)/2))
        ELSEIF( NCURVE .EQ. 4 )THEN
         YY  = YY0 - AMPLI*SIN(PI*X(I,(NJ+1)/2))
c        YY  =-YY0 + AMPLI*SIN(PI*X(I,(NJ+1)/2))
        ENDIF
         YZ = 10000.
         YZB= 10000.
         YZT= 10000.
      DO 45 J = -1,NJ1
            YYY = YY - Y(I,J)
            DY  = (Y(I,J+1) - Y(I,J))
            IF( ABS(YYY) .LE. YZ )THEN
               YZ = ABS(YYY)
               JJ = J
            ENDIF
          IF( NCURVE .EQ. 3 )THEN
            YYYT = YYT - Y(I,J)
            YYYB = YYB - Y(I,J)
            IF( ABS(YYYB) .LE. YZB )THEN
               YZB = ABS(YYYB)
               JJB = J
            ENDIF
            IF( ABS(YYYT) .LE. YZT )THEN
               YZT = ABS(YYYT)
               JJT = J
            ENDIF
          ENDIF
45       CONTINUE
         IF( ABS(Y(I,JJ+1)-YY) .LT. DY )THEN
            PSI(I,JJ) = X2
              U(I,JJ) = U0*URATIO
         ELSE
            PSI(I,JJ) = X1
              U(I,JJ) = U0
         ENDIF
c        IF( NCURVE .EQ. 3 )THEN
c         IF( ABS(Y(I,JJ+1)-YYT) .LT. DY )THEN
c           PSI(I,JJT) = X1
c             U(I,JJT) = U0*URATIO
c         ELSE
c           PSI(I,JJT) = X2
c             U(I,JJT) = U0
c         ENDIF
c         IF( ABS(Y(I,JJ+1)-YYB) .LT. DY )THEN
c           PSI(I,JJB) = X1
c             U(I,JJB) = U0*URATIO
c         ELSE
c           PSI(I,JJB) = X2
c             U(I,JJB) = U0
c         ENDIF
c        ENDIF
         DO 46 J = JJ-1,-1,-1
            PSI(I,J) = X2
              U(I,J) = U0*URATIO
46       CONTINUE
         DO 47 J = JJ+1,NJ2
            PSI(I,J) = X1
              U(I,J) = U0
47       CONTINUE
        IF( NCURVE .EQ. 3 )THEN
         DO 50 J = JJB-1,-1,-1
            PSI(I,J) = X1
              U(I,J) = U0
50       CONTINUE
         DO 52 J = JJB,JJT
            PSI(I,J) = X2
              U(I,J) = U0*URATIO
52       CONTINUE
         DO 55 J = JJT+1,NJ2
            PSI(I,J) = X1
              U(I,J) = U0
55       CONTINUE
        ENDIF
40    CONTINUE
      DO 59 J = -1,NJ2
         PSI(-1,J) = PSI(NI-2,J)
         PSI( 0,J) = PSI(NI-1,J)
         PSI(NI  ,J) = PSI(1,J)
         PSI(NI+1,J) = PSI(2,J)
         PSI(NI+2,J) = PSI(3,J)
           U(-1,J) = U(NI-2,J)
           U( 0,J) = U(NI-1,J)
           U(NI  ,J) = U(1,J)
           U(NI+1,J) = U(2,J)
           U(NI+2,J) = U(3,J)
59    CONTINUE
      DO 60 I = -1,NI2
         PSI(I,-1) = PSI(I,1)
         PSI(I, 0) = PSI(I,1)
         PSI(I,NJ  ) = PSI(I,NJ)
         PSI(I,NJ+1) = PSI(I,NJ)
         PSI(I,NJ+2) = PSI(I,NJ)
           U(I,-1) = U(I,1)
           U(I, 0) = U(I,1)
           U(I,NJ  ) = U(I,NJ)
           U(I,NJ+1) = U(I,NJ)
           U(I,NJ+2) = U(I,NJ)
60    CONTINUE
      DO 30 IJ  =-1,MIJ2
       R(IJ,-1) = R0
       V(IJ,-1) = V0
       T(IJ,-1) = T0
       P(IJ,-1) = P0
       E(IJ,-1) = E0
   30 CONTINUE
      ENDIF
C
100   CONTINUE
      IF( ICYL .EQ. 1 .AND. LSTD .EQ. 0 )THEN
      DO 20 IJ  =-1,MIJ2
       IF(X(IJ,-1) .LE. XST) THEN
       R(IJ,-1) = DL0
       U(IJ,-1) = UL0
       V(IJ,-1) = VL0
       E(IJ,-1) = EL0
       P(IJ,-1) = PL0
      ELSE
       R(IJ,-1) = DR0
       U(IJ,-1) = UR0
       V(IJ,-1) = VR0
       E(IJ,-1) = ER0
       P(IJ,-1) = PR0
      ENDIF
   20 CONTINUE
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE ITFACE(IJP,R,U,V,QJ,PSI,XIX,XIY,ETAX,ETAY,DRS)
C
      PARAMETER (MAXIJP=516)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
      COMMON /BC01 / IDBC(30,20),FDBC(10,20)
      COMMON /JSY/ JSYM,JTHR
      DIMENSION IJP(MXBP)
      DIMENSION    U(-1:NI2,-1:NJ2),V(-1:NI2,-1:NJ2)
      DIMENSION  XIX(-1:NI2,-1:NJ2), XIY(-1:NI2,-1:NJ2)
      DIMENSION ETAX(-1:NI2,-1:NJ2),ETAY(-1:NI2,-1:NJ2)
      DIMENSION PSII(-1:MAXIJP,-1:MAXIJP),QJ(-1:NI2,-1:NJ2)
      DIMENSION PSI(-1:NI2,-1:NJ2),DPSI(-1:MAXIJP)
      DIMENSION PSIII(-1:MAXIJP,-1:MAXIJP)
      DIMENSION DRS(-1:NI2,-1:NJ2,5),R(-1:NI2,-1:NJ2)
C
C*****  Starting of the JSYM=1 for nonconservation form with Osher's method ******
      IF( JSYM .EQ. 1 )THEN
      DO 10 J = -1,NJ2
      DO 5  I =  0,NI1
         PSIF    = ABS( PSI(I+1,J)-PSI(I  ,J) )
         PSIB    = ABS( PSI(I  ,J)-PSI(I-1,J) )
         DPSI(I) = MIN( PSIF,PSIB )
5     CONTINUE
      DPSI(1) = DPSI(2)
      DPSI(NI)= DPSI(NI-1)
      DO 15 I =  0,NI1
         UU     = U(I,J)*XIX(I,J)+V(I,J)*XIY(I,J)
         UP     = 0.5*( UU + ABS(UU) )
         UM     = 0.5*( UU - ABS(UU) )
       PSI1     = PSI(I  ,J)+0.5*DPSI(I)
       PSI2     = PSI(I-1,J)+0.5*DPSI(I-1)
       PSI3     = PSI(I+1,J)-0.5*DPSI(I+1)
       PSI4     = PSI(I  ,J)-0.5*DPSI(I)
      PSII(I,J) =-UP*( PSI1 - PSI2 ) - UM*( PSI3 - PSI4 )
15    CONTINUE
10    CONTINUE
C
      DO 20 I = -1,NI2
      DO 25 J =  0,NJ1
         PSIF    = ABS( PSI(I,J+1)-PSI(I,  J) )
         PSIB    = ABS( PSI(I,  J)-PSI(I,J-1) )
         DPSI(J) = MIN( PSIF,PSIB )
25    CONTINUE
      DPSI(1 ) = DPSI(2)
      DPSI(NJ) = DPSI(NJ-1)
      DO 30 J = 0,NJ1
         VV     = U(I,J)*ETAX(I,J)+V(I,J)*ETAY(I,J)
         VP     = 0.5*( VV + ABS(VV) )
         VM     = 0.5*( VV - ABS(VV) )
         PSI1   = PSI(I,  J)+0.5*DPSI(J)
         PSI2   = PSI(I,J-1)+0.5*DPSI(J-1)
         PSI3   = PSI(I,J+1)-0.5*DPSI(J+1)
         PSI4   = PSI(I,  J)-0.5*DPSI(J)
        PSII(I,J) = PSII(I,J)-VP*( PSI1 - PSI2 )-VM*( PSI3 - PSI4 )
30    CONTINUE
20    CONTINUE
C
      DO 35 J = 0,NJ1
      DO 35 I = 0,NI1
         PSII(I,J) = PSI(I,J) + DTCFL*PSII(I,J)
35    CONTINUE
C
      DO 37 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IDA    = IDBC( 2,IBC)
      ISA    = IDBC( 3,IBC)
      ICA    = IDBC( 4,IBC)
      IAS    = IDBC( 5,IBC)
      IAE    = IDBC( 6,IBC)
      IPA    = IDBC( 7,IBC)
      LENA   = (IAE - IAS)*ICA + 1
      IF(IBTYPE .EQ. 4 .OR. IBTYPE .EQ. 7) THEN
      IDB    = IDBC(12,IBC)
      ISB    = IDBC(13,IBC)
      ICB    = IDBC(14,IBC)
      IBS    = IDBC(15,IBC)
      IBE    = IDBC(16,IBC)
      IPB    = IDBC( 8,IBC)
      LENB   = (IBE - IBS)*ICB + 1
      ENDIF
C
      DO 38 II = 1 , LENA
          IC    =  IJP(IPA+LENA*4+II)
          JC    =  IJP(IPA+LENA*5+II)
          IA    =  IJP(IPA+LENA*0+II)
          JA    =  IJP(IPA+LENA*1+II)
          IB    =  IJP(IPA+LENA*2+II)
          JB    =  IJP(IPA+LENA*3+II)
          ID    =  IJP(IPA+LENA*6+II)
          JD    =  IJP(IPA+LENA*7+II)
          IE    =  IJP(IPA+LENA*8+II)
          JE    =  IJP(IPA+LENA*9+II)
      IF(IBTYPE .EQ. 0)THEN
        PSII(IB,JB) = PSII(IC,JC)
        PSII(IA,JA) = PSII(IC,JC)
      ELSEIF(IBTYPE .EQ. 1 ) THEN
        PSII(IB,JB) = PSII(IC,JC)
        PSII(IA,JA) = PSII(IC,JC)
      ELSEIF(IBTYPE .EQ. 6 ) THEN
        PSII(IB,JB) = PSII(ID,JD)
        PSII(IA,JA) = PSII(IE,JE)
        PSII(IC,JC) = 0.5*( PSII(IB,JB) + PSII(ID,JD) )
      ELSEIF(IBTYPE .EQ. 2 .OR. IBTYPE .EQ. 3) THEN
        PSII(IB,JB) = PSII(ID,JD)
        PSII(IA,JA) = PSII(IE,JE)
        PSII(IC,JC) = 0.5*( PSII(IB,JB) + PSII(ID,JD) )
      ELSEIF(IBTYPE .EQ. 4 .OR. IBTYPE .EQ. 7) THEN
          IC2   =  IJP(IPB+LENA*4+II)
          JC2   =  IJP(IPB+LENA*5+II)
          IA2   =  IJP(IPB+LENA*0+II)
          JA2   =  IJP(IPB+LENA*1+II)
          IB2   =  IJP(IPB+LENA*2+II)
          JB2   =  IJP(IPB+LENA*3+II)
          ID2   =  IJP(IPB+LENA*6+II)
          JD2   =  IJP(IPB+LENA*7+II)
          IE2   =  IJP(IPB+LENA*8+II)
          JE2   =  IJP(IPB+LENA*9+II)
        PSII(IB,JB)   = PSII(ID2,JD2)
        PSII(IA,JA)   = PSII(IE2,JE2)
        PSII(IB2,JB2) = PSII(ID,JD)
        PSII(IA2,JA2) = PSII(IE,JE)
C
C     FOR SIMPLE SET EQUAL
C
        PSII(IC,JC)   = 0.5*( PSII(ID,JD) + PSII(ID2,JD2) )
        PSII(IC2,JC2) = PSII(IC,JC)
C
      ENDIF
  38  CONTINUE
  37  CONTINUE
C
      DO 40 J = -1,NJ2
      DO 45 I =  0,NI1
         PSIF    = ABS( PSII(I+1,J)-PSII(I  ,J) )
         PSIB    = ABS( PSII(I  ,J)-PSII(I-1,J) )
         DPSI(I) = MIN( PSIF,PSIB )
45    CONTINUE
      DPSI(1)  = DPSI(2)
      DPSI(NI) = DPSI(NI-1)
      DO 50 I =  0,NI1
         UU     = U(I,J)*XIX(I,J)+V(I,J)*XIY(I,J)
         UP     = 0.5*( UU + ABS(UU) )
         UM     = 0.5*( UU - ABS(UU) )
       PSI1     = PSII(I  ,J)+0.5*DPSI(I)
       PSI2     = PSII(I-1,J)+0.5*DPSI(I-1)
       PSI3     = PSII(I+1,J)-0.5*DPSI(I+1)
       PSI4     = PSII(I  ,J)-0.5*DPSI(I)
      PSIII(I,J) =-UP*( PSI1 - PSI2 ) - UM*( PSI3 - PSI4 )
50    CONTINUE
40    CONTINUE
C
      DO 60 I = -1,NI2
      DO 65 J =  0,NJ1
         PSIF    = ABS( PSII(I,J+1)-PSII(I,  J) )
         PSIB    = ABS( PSII(I,  J)-PSII(I,J-1) )
         DPSI(J) = MIN( PSIF,PSIB )
65    CONTINUE
      DPSI(1)  = DPSI(2)
      DPSI(NJ) = DPSI(NJ-1)
      DO 70 J = 0,NJ1
         VV     = U(I,J)*ETAX(I,J)+V(I,J)*ETAY(I,J)
         VP     = 0.5*( VV + ABS(VV) )
         VM     = 0.5*( VV - ABS(VV) )
         PSI1   = PSII(I,  J)+0.5*DPSI(J)
         PSI2   = PSII(I,J-1)+0.5*DPSI(J-1)
         PSI3   = PSII(I,J+1)-0.5*DPSI(J+1)
         PSI4   = PSII(I,  J)-0.5*DPSI(J)
        PSIII(I,J) = PSIII(I,J)-VP*( PSI1 - PSI2 )-VM*( PSI3 - PSI4 )
70    CONTINUE
60    CONTINUE
C
      DO 80 J = 0,NJ1
      DO 80 I = 0,NI1
         PSI(I,J) = 0.5*(PSI(I,J)+PSII(I,J))+0.5*DTCFL*PSIII(I,J)
80    CONTINUE
c
      DO 85 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IDA    = IDBC( 2,IBC)
      ISA    = IDBC( 3,IBC)
      ICA    = IDBC( 4,IBC)
      IAS    = IDBC( 5,IBC)
      IAE    = IDBC( 6,IBC)
      IPA    = IDBC( 7,IBC)
      LENA   = (IAE - IAS)*ICA + 1
      IF(IBTYPE .EQ. 4 .OR. IBTYPE .EQ. 7) THEN
      IDB    = IDBC(12,IBC)
      ISB    = IDBC(13,IBC)
      ICB    = IDBC(14,IBC)
      IBS    = IDBC(15,IBC)
      IBE    = IDBC(16,IBC)
      IPB    = IDBC( 8,IBC)
      LENB   = (IBE - IBS)*ICB + 1
      ENDIF
C
      DO 90 II = 1 , LENA
          IC    =  IJP(IPA+LENA*4+II)
          JC    =  IJP(IPA+LENA*5+II)
          IA    =  IJP(IPA+LENA*0+II)
          JA    =  IJP(IPA+LENA*1+II)
          IB    =  IJP(IPA+LENA*2+II)
          JB    =  IJP(IPA+LENA*3+II)
          ID    =  IJP(IPA+LENA*6+II)
          JD    =  IJP(IPA+LENA*7+II)
          IE    =  IJP(IPA+LENA*8+II)
          JE    =  IJP(IPA+LENA*9+II)
      IF(IBTYPE .EQ. 0)THEN
        PSI(IB,JB) = PSI(IC,JC)
        PSI(IA,JA) = PSI(IC,JC)
      ELSEIF(IBTYPE .EQ. 1 ) THEN
        PSI(IB,JB) = PSI(IC,JC)
        PSI(IA,JA) = PSI(IC,JC)
      ELSEIF(IBTYPE .EQ. 6 ) THEN
        PSI(IB,JB) = PSI(ID,JD)
        PSI(IA,JA) = PSI(IE,JE)
        PSI(IC,JC) = 0.5*( PSI(IB,JB) + PSI(ID,JD) )
      ELSEIF(IBTYPE .EQ. 2 .OR. IBTYPE .EQ. 3) THEN
        PSI(IB,JB) = PSI(ID,JD)
        PSI(IA,JA) = PSI(IE,JE)
        PSI(IC,JC) = 0.5*( PSI(IB,JB) + PSI(ID,JD) )
      ELSEIF(IBTYPE .EQ. 4 .OR. IBTYPE .EQ. 7) THEN
          IC2   =  IJP(IPB+LENA*4+II)
          JC2   =  IJP(IPB+LENA*5+II)
          IA2   =  IJP(IPB+LENA*0+II)
          JA2   =  IJP(IPB+LENA*1+II)
          IB2   =  IJP(IPB+LENA*2+II)
          JB2   =  IJP(IPB+LENA*3+II)
          ID2   =  IJP(IPB+LENA*6+II)
          JD2   =  IJP(IPB+LENA*7+II)
          IE2   =  IJP(IPB+LENA*8+II)
          JE2   =  IJP(IPB+LENA*9+II)
        PSI(IB,JB)   = PSI(ID2,JD2)
        PSI(IA,JA)   = PSI(IE2,JE2)
        PSI(IB2,JB2) = PSI(ID,JD)
        PSI(IA2,JA2) = PSI(IE,JE)
C
C     FOR SIMPLE SET EQUAL
C
        PSI(IC,JC)   = 0.5*( PSI(ID,JD) + PSI(ID2,JD2) )
        PSI(IC2,JC2) = PSI(IC,JC)
C
      ENDIF
  90  CONTINUE
  85  CONTINUE
C
      ENDIF
C****************   END OF THE JSYM=1       *************************
C
C****** Starting of the JSYM=2 for nonconservation form with TVD's method *********
      IF( JSYM .EQ. 2 )THEN
        DO 200 I = -1,NI2
        DO 200 J = -1,NJ2
           DRS(I,J,1) = 1.
           DRS(I,J,2) = 0.
200     CONTINUE
        DO 210 I = -1,NI2
        DO 210 J = -1,NJ2
         DRS(I,J,3) = U(I,J)*XIX(I,J) +V(I,J)*XIY(I,J)
         DRS(I,J,4) = U(I,J)*ETAX(I,J)+V(I,J)*ETAY(I,J)
210     CONTINUE
        DO 220 J = -1,NJ2
        DO 225 I =  0,NI1
          DRS(I,J,2) = DRS(I,J,2)+PSI(I,J)*DTCFL*
     1                 (DRS(I+1,J,3)-DRS(I-1,J,3))/2.
225     CONTINUE
          I = -1
          DRS(I,J,2) = DRS(I,J,2)+PSI(I,J)*DTCFL*
     1                (-DRS(I+2,J,3)+4.*DRS(I+1,J,3)-3.*DRS(I,J,3))/2.
          I = NI2
          DRS(I,J,2) = DRS(I,J,2)+PSI(I,J)*DTCFL*
     1                ( DRS(I-2,J,3)-4.*DRS(I-1,J,3)+3.*DRS(I,J,3))/2.
 220    CONTINUE
        DO 230 I = -1,NI2
        DO 235 J =  0,NJ1
          DRS(I,J,2) = DRS(I,J,2)+PSI(I,J)*DTCFL*
     1                 (DRS(I,J+1,4)-DRS(I,J-1,4))/2.
 235    CONTINUE
          J = -1
          DRS(I,J,2) = DRS(I,J,2)+PSI(I,J)*DTCFL*
     1                (-DRS(I,J+2,4)+4.*DRS(I,J+1,4)-3.*DRS(I,J,4))/2.
          J = NJ2
          DRS(I,J,2) = DRS(I,J,2)+PSI(I,J)*DTCFL*
     1                ( DRS(I,J-2,4)-4.*DRS(I,J-1,4)+3.*DRS(I,J,4))/2.
 230    CONTINUE
C
        IF( JTHR .EQ. 2 )THEN
         CALL FVITVD2(DRS(-1,-1,3), DRS(-1,-1,1), PSI, DRS(-1,-1,2), 1)
         CALL FVJTVD2(DRS(-1,-1,4), DRS(-1,-1,1), PSI, DRS(-1,-1,2), 1)
        ENDIF
        IF( JTHR .EQ. 3 )THEN
         CALL FVIEN2(DRS(-1,-1,3), DRS(-1,-1,1), PSI,  DRS(-1,-1,2), 1)
         CALL FVJEN2(DRS(-1,-1,4), DRS(-1,-1,1), PSI,  DRS(-1,-1,2), 1)
        ENDIF
      IF(JTHR .GE. 5 .AND. JTHR .LE. 7) THEN
C
      IF(JTHR .EQ. 5) THEN
       XI     = 0.
       ETA    = 0.
      ELSEIF(JTHR .EQ. 6) THEN
       XI     = 0.
       ETA    = 0.5
      ELSEIF(JTHR .EQ. 7) THEN
       XI     = 1.0
       ETA    = 0.
      ENDIF
        CALL FVIENO(DRS(-1,-1,3), DRS(-1,-1,1), PSI, DRS(-1,-1,2),
     1              XI,ETA,ILMT,L1)
        CALL FVJENO(DRS(-1,-1,4), DRS(-1,-1,1), PSI, DRS(-1,-1,2),
     1              XI,ETA,ILMT,L1)
      ENDIF
      IF( JTHR .EQ. 11 )THEN
        CALL FDITVD2(DRS(-1,-1,3), DRS(-1,-1,1), PSI, DRS(-1,-1,2), 1)
        CALL FDJTVD2(DRS(-1,-1,3), DRS(-1,-1,1), PSI, DRS(-1,-1,2), 1)
      ENDIF
C
        DO 250 I = 0,NI1
        DO 250 J = 0,NJ1
           PSI(I,J) = PSI(I,J) + DRS(I,J,2)
250     CONTINUE
C
      DO 260 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IDA    = IDBC( 2,IBC)
      ISA    = IDBC( 3,IBC)
      ICA    = IDBC( 4,IBC)
      IAS    = IDBC( 5,IBC)
      IAE    = IDBC( 6,IBC)
      IPA    = IDBC( 7,IBC)
      LENA   = (IAE - IAS)*ICA + 1
      IF(IBTYPE .EQ. 4 .OR. IBTYPE .EQ. 7) THEN
      IDB    = IDBC(12,IBC)
      ISB    = IDBC(13,IBC)
      ICB    = IDBC(14,IBC)
      IBS    = IDBC(15,IBC)
      IBE    = IDBC(16,IBC)
      IPB    = IDBC( 8,IBC)
      LENB   = (IBE - IBS)*ICB + 1
      ENDIF
C
      DO 270 II = 1 , LENA
          IC    =  IJP(IPA+LENA*4+II)
          JC    =  IJP(IPA+LENA*5+II)
          IA    =  IJP(IPA+LENA*0+II)
          JA    =  IJP(IPA+LENA*1+II)
          IB    =  IJP(IPA+LENA*2+II)
          JB    =  IJP(IPA+LENA*3+II)
          ID    =  IJP(IPA+LENA*6+II)
          JD    =  IJP(IPA+LENA*7+II)
          IE    =  IJP(IPA+LENA*8+II)
          JE    =  IJP(IPA+LENA*9+II)
      IF(IBTYPE .EQ. 0)THEN
        PSI(IB,JB) = PSI(IC,JC)
        PSI(IA,JA) = PSI(IC,JC)
      ELSEIF(IBTYPE .EQ. 1 ) THEN
        PSI(IB,JB) = PSI(IC,JC)
        PSI(IA,JA) = PSI(IC,JC)
      ELSEIF(IBTYPE .EQ. 6 ) THEN
        PSI(IB,JB) = PSI(ID,JD)
        PSI(IA,JA) = PSI(IE,JE)
        PSI(IC,JC) = 0.5*( PSI(IB,JB) + PSI(ID,JD) )
      ELSEIF(IBTYPE .EQ. 2 .OR. IBTYPE .EQ. 3) THEN
        PSI(IB,JB) = PSI(ID,JD)
        PSI(IA,JA) = PSI(IE,JE)
        PSI(IC,JC) = 0.5*( PSI(IB,JB) + PSI(ID,JD) )
      ELSEIF(IBTYPE .EQ. 4 .OR. IBTYPE .EQ. 7) THEN
          IC2   =  IJP(IPB+LENA*4+II)
          JC2   =  IJP(IPB+LENA*5+II)
          IA2   =  IJP(IPB+LENA*0+II)
          JA2   =  IJP(IPB+LENA*1+II)
          IB2   =  IJP(IPB+LENA*2+II)
          JB2   =  IJP(IPB+LENA*3+II)
          ID2   =  IJP(IPB+LENA*6+II)
          JD2   =  IJP(IPB+LENA*7+II)
          IE2   =  IJP(IPB+LENA*8+II)
          JE2   =  IJP(IPB+LENA*9+II)
        PSI(IB,JB)   = PSI(ID2,JD2)
        PSI(IA,JA)   = PSI(IE2,JE2)
        PSI(IB2,JB2) = PSI(ID,JD)
        PSI(IA2,JA2) = PSI(IE,JE)
        PSI(IC,JC)   = 0.5*( PSI(ID,JD) + PSI(ID2,JD2) )
        PSI(IC2,JC2) = PSI(IC,JC)
      ENDIF
C
  270 CONTINUE
  260 CONTINUE
      ENDIF
C****************   END OF THE JSYM=2       *************************
C
C*****************  Starting of the JSYM=3 for Conservation form *******************
      IF( JSYM .EQ. 3 )THEN
        DO 300 I = -1,NI2
        DO 300 J = -1,NJ2
           DRS(I,J,2) = 0.
           XU2        = U(I,J)*U(I,J) + V(I,J)*V(I,J)
           DRS(I,J,5) = 1./SQRT(1. - XU2)
300     CONTINUE
        DO 310 I = -1,NI2
        DO 310 J = -1,NJ2
         DRS(I,J,3) = U(I,J)*XIX(I,J) +V(I,J)*XIY(I,J)
         DRS(I,J,4) = U(I,J)*ETAX(I,J)+V(I,J)*ETAY(I,J)
         DRS(I,J,1) = PSI(I,J)*DRS(I,J,5)*R(I,J)/QJ(I,J)
310     CONTINUE
C
        IF( JTHR .EQ. 2 )THEN
         CALL FVITVD2(DRS(-1,-1,3), QJ, DRS(-1,-1,1), DRS(-1,-1,2), 1)
         CALL FVJTVD2(DRS(-1,-1,4), QJ, DRS(-1,-1,1), DRS(-1,-1,2), 1)
        ENDIF
        IF( JTHR .EQ. 3 )THEN
         CALL FVIEN2(DRS(-1,-1,3), QJ, DRS(-1,-1,1),  DRS(-1,-1,2), 1)
         CALL FVJEN2(DRS(-1,-1,4), QJ, DRS(-1,-1,1),  DRS(-1,-1,2), 1)
        ENDIF
      IF(JTHR .GE. 5 .AND. JTHR .LE. 7) THEN
C
      IF(JTHR .EQ. 5) THEN
       XI     = 0.
       ETA    = 0.
      ELSEIF(JTHR .EQ. 6) THEN
       XI     = 0.
       ETA    = 0.5
      ELSEIF(JTHR .EQ. 7) THEN
       XI     = 1.0
       ETA    = 0.
      ENDIF
        CALL FVIENO(DRS(-1,-1,3), QJ, DRS(-1,-1,1), DRS(-1,-1,2),
     1              XI,ETA,ILMT,L1)
        CALL FVJENO(DRS(-1,-1,4), QJ, DRS(-1,-1,1), DRS(-1,-1,2),
     1              XI,ETA,ILMT,L1)
      ENDIF
      IF( JTHR .EQ. 11 )THEN
        CALL FDITVD2(DRS(-1,-1,3), QJ, DRS(-1,-1,1), DRS(-1,-1,2), 1)
        CALL FDJTVD2(DRS(-1,-1,3), QJ, DRS(-1,-1,1), DRS(-1,-1,2), 1)
      ENDIF
C
        DO 320 I = 0,NI1
        DO 320 J = 0,NJ1
           DRS(I,J,1) = DRS(I,J,1) + DRS(I,J,2)
           PSI(I,J) = DRS(I,J,1)*QJ(I,J)/R(I,J)/DRS(I,J,5)
320     CONTINUE
C
      DO 330 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IDA    = IDBC( 2,IBC)
      ISA    = IDBC( 3,IBC)
      ICA    = IDBC( 4,IBC)
      IAS    = IDBC( 5,IBC)
      IAE    = IDBC( 6,IBC)
      IPA    = IDBC( 7,IBC)
      LENA   = (IAE - IAS)*ICA + 1
      IF(IBTYPE .EQ. 4 .OR. IBTYPE .EQ. 7) THEN
      IDB    = IDBC(12,IBC)
      ISB    = IDBC(13,IBC)
      ICB    = IDBC(14,IBC)
      IBS    = IDBC(15,IBC)
      IBE    = IDBC(16,IBC)
      IPB    = IDBC( 8,IBC)
      LENB   = (IBE - IBS)*ICB + 1
      ENDIF
C
      DO 340 II = 1 , LENA
          IC    =  IJP(IPA+LENA*4+II)
          JC    =  IJP(IPA+LENA*5+II)
          IA    =  IJP(IPA+LENA*0+II)
          JA    =  IJP(IPA+LENA*1+II)
          IB    =  IJP(IPA+LENA*2+II)
          JB    =  IJP(IPA+LENA*3+II)
          ID    =  IJP(IPA+LENA*6+II)
          JD    =  IJP(IPA+LENA*7+II)
          IE    =  IJP(IPA+LENA*8+II)
          JE    =  IJP(IPA+LENA*9+II)
      IF(IBTYPE .EQ. 0)THEN
        PSI(IB,JB) = PSI(IC,JC)
        PSI(IA,JA) = PSI(IC,JC)
      ELSEIF(IBTYPE .EQ. 1 ) THEN
        PSI(IB,JB) = PSI(IC,JC)
        PSI(IA,JA) = PSI(IC,JC)
      ELSEIF(IBTYPE .EQ. 6 ) THEN
        PSI(IB,JB) = PSI(ID,JD)
        PSI(IA,JA) = PSI(IE,JE)
        PSI(IC,JC) = 0.5*( PSI(IB,JB) + PSI(ID,JD) )
      ELSEIF(IBTYPE .EQ. 2 .OR. IBTYPE .EQ. 3) THEN
        PSI(IB,JB) = PSI(ID,JD)
        PSI(IA,JA) = PSI(IE,JE)
        PSI(IC,JC) = 0.5*( PSI(IB,JB) + PSI(ID,JD) )
      ELSEIF(IBTYPE .EQ. 4 .OR. IBTYPE .EQ. 7) THEN
          IC2   =  IJP(IPB+LENA*4+II)
          JC2   =  IJP(IPB+LENA*5+II)
          IA2   =  IJP(IPB+LENA*0+II)
          JA2   =  IJP(IPB+LENA*1+II)
          IB2   =  IJP(IPB+LENA*2+II)
          JB2   =  IJP(IPB+LENA*3+II)
          ID2   =  IJP(IPB+LENA*6+II)
          JD2   =  IJP(IPB+LENA*7+II)
          IE2   =  IJP(IPB+LENA*8+II)
          JE2   =  IJP(IPB+LENA*9+II)
        PSI(IB,JB)   = PSI(ID2,JD2)
        PSI(IA,JA)   = PSI(IE2,JE2)
        PSI(IB2,JB2) = PSI(ID,JD)
        PSI(IA2,JA2) = PSI(IE,JE)
        PSI(IC,JC)   = 0.5*( PSI(ID,JD) + PSI(ID2,JD2) )
        PSI(IC2,JC2) = PSI(IC,JC)
      ENDIF
C
  340 CONTINUE
  330 CONTINUE
      ENDIF
C****************   END OF THE JSYM=3       *************************
C
      RETURN
      END
      SUBROUTINE METJAC(X,Y,XIX,XIY,ETAX,ETAY,QJ,IJP)
      COMMON /TAPE / INAME,IGRID,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /NEXT / NEXTLOC
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /BC01 / IDBC(30,20),FDBC(10,20)
      DIMENSION    X(-1:NI2,-1:NJ2),   Y(-1:NI2,-1:NJ2)
      DIMENSION  XIX(-1:NI2,-1:NJ2), XIY(-1:NI2,-1:NJ2)
      DIMENSION ETAX(-1:NI2,-1:NJ2),ETAY(-1:NI2,-1:NJ2)
      DIMENSION   QJ(-1:NI2,-1:NJ2)
      DIMENSION IJP(MXBP)
C
      DO 10 I = 0,NI1
      DO 10 J = 0,NJ1
        XI  = (X(I+1,J) - X(I-1,J))*0.5
        YI  = (Y(I+1,J) - Y(I-1,J))*0.5
        XJ  = (X(I,J+1) - X(I,J-1))*0.5
        YJ  = (Y(I,J+1) - Y(I,J-1))*0.5
        VOL = XI*YJ - XJ*YI
          QJ(I,J)  = 1./(VOL + 1.E-25)
         XIX(I,J)  =  YJ*QJ(I,J)
         XIY(I,J)  = -XJ*QJ(I,J)
        ETAX(I,J)  = -YI*QJ(I,J)
        ETAY(I,J)  =  XI*QJ(I,J)
   10 CONTINUE
C
      J       = -1
      DO 20 I = 0 , NI1
         XI   =  (X(I+1,J) - X(I-1,J))*0.5
         YI   =  (Y(I+1,J) - Y(I-1,J))*0.5
         XJ   = -(3.0*X(I,J) -4.0*X(I,J+1) +X(I,J+2))*0.5
         YJ   = -(3.0*Y(I,J) -4.0*Y(I,J+1) +Y(I,J+2))*0.5
         VOL  =   XI*YJ - XJ*YI
          QJ(I,J)  = 1./(VOL + 1.E-25)
         XIX(I,J)  =  YJ*QJ(I,J)
         XIY(I,J)  = -XJ*QJ(I,J)
        ETAX(I,J)  = -YI*QJ(I,J)
        ETAY(I,J)  =  XI*QJ(I,J)
   20 CONTINUE
C
      J       = NJ2
      DO 30 I = 0 , NI1
         XI   =  (X(I+1,J) - X(I-1,J))*0.5
         YI   =  (Y(I+1,J) - Y(I-1,J))*0.5
         XJ   =  (3.0*X(I,J) -4.0*X(I,J-1) +X(I,J-2))*0.5
         YJ   =  (3.0*Y(I,J) -4.0*Y(I,J-1) +Y(I,J-2))*0.5
         VOL  =   XI*YJ - XJ*YI
          QJ(I,J)  = 1./(VOL + 1.E-25)
         XIX(I,J)  =  YJ*QJ(I,J)
         XIY(I,J)  = -XJ*QJ(I,J)
        ETAX(I,J)  = -YI*QJ(I,J)
        ETAY(I,J)  =  XI*QJ(I,J)
   30 CONTINUE
      I       = -1
C
      DO 40 J = 0 , NJ1
         XI   = -(3.0*X(I,J) -4.0*X(I+1,J) +X(I+2,J))*0.5
         YI   = -(3.0*Y(I,J) -4.0*Y(I+1,J) +Y(I+2,J))*0.5
         XJ   =  (X(I,J+1) - X(I,J-1))*0.5
         YJ   =  (Y(I,J+1) - Y(I,J-1))*0.5
         VOL  =   XI*YJ - XJ*YI
          QJ(I,J)  = 1./(VOL + 1.E-25)
         XIX(I,J)  =  YJ*QJ(I,J)
         XIY(I,J)  = -XJ*QJ(I,J)
        ETAX(I,J)  = -YI*QJ(I,J)
        ETAY(I,J)  =  XI*QJ(I,J)
   40 CONTINUE
C
      I       = NI2
      DO 50 J = 0 , NJ1
         XI   =  (3.0*X(I,J) -4.0*X(I-1,J) +X(I-2,J))*0.5
         YI   =  (3.0*Y(I,J) -4.0*Y(I-1,J) +Y(I-2,J))*0.5
         XJ   =  (X(I,J+1) - X(I,J-1))*0.5
         YJ   =  (Y(I,J+1) - Y(I,J-1))*0.5
         VOL  =   XI*YJ - XJ*YI
          QJ(I,J)  = 1./(VOL + 1.E-25)
         XIX(I,J)  =  YJ*QJ(I,J)
         XIY(I,J)  = -XJ*QJ(I,J)
        ETAX(I,J)  = -YI*QJ(I,J)
        ETAY(I,J)  =  XI*QJ(I,J)
   50 CONTINUE
C
      I = -1
      J = -1
         XI   = -(3.0*X(I,J) -4.0*X(I+1,J) +X(I+2,J))*0.5
         YI   = -(3.0*Y(I,J) -4.0*Y(I+1,J) +Y(I+2,J))*0.5
         XJ   = -(3.0*X(I,J) -4.0*X(I,J+1) +X(I,J+2))*0.5
         YJ   = -(3.0*Y(I,J) -4.0*Y(I,J+1) +Y(I,J+2))*0.5
         VOL = XI*YJ - XJ*YI
          QJ(I,J)  = 1./(VOL + 1.E-25)
         XIX(I,J)  =  YJ*QJ(I,J)
         XIY(I,J)  = -XJ*QJ(I,J)
        ETAX(I,J)  = -YI*QJ(I,J)
        ETAY(I,J)  =  XI*QJ(I,J)
C
      I = NI2
      J = -1
         XI   =  (3.0*X(I,J) -4.0*X(I-1,J) +X(I-2,J))*0.5
         YI   =  (3.0*Y(I,J) -4.0*Y(I-1,J) +Y(I-2,J))*0.5
         XJ   = -(3.0*X(I,J) -4.0*X(I,J+1) +X(I,J+2))*0.5
         YJ   = -(3.0*Y(I,J) -4.0*Y(I,J+1) +Y(I,J+2))*0.5
         VOL = XI*YJ - XJ*YI
          QJ(I,J)  = 1./(VOL + 1.E-25)
         XIX(I,J)  =  YJ*QJ(I,J)
         XIY(I,J)  = -XJ*QJ(I,J)
        ETAX(I,J)  = -YI*QJ(I,J)
        ETAY(I,J)  =  XI*QJ(I,J)
C
      I = -1
      J = NJ2
         XI   = -(3.0*X(I,J) -4.0*X(I+1,J) +X(I+2,J))*0.5
         YI   = -(3.0*Y(I,J) -4.0*Y(I+1,J) +Y(I+2,J))*0.5
         XJ   =  (3.0*X(I,J) -4.0*X(I,J-1) +X(I,J-2))*0.5
         YJ   =  (3.0*Y(I,J) -4.0*Y(I,J-1) +Y(I,J-2))*0.5
         VOL = XI*YJ - XJ*YI
          QJ(I,J)  = 1./(VOL + 1.E-25)
         XIX(I,J)  =  YJ*QJ(I,J)
         XIY(I,J)  = -XJ*QJ(I,J)
        ETAX(I,J)  = -YI*QJ(I,J)
        ETAY(I,J)  =  XI*QJ(I,J)
C
      I = NI2
      J = NJ2
         XI   =  (3.0*X(I,J) -4.0*X(I-1,J) +X(I-2,J))*0.5
         YI   =  (3.0*Y(I,J) -4.0*Y(I-1,J) +Y(I-2,J))*0.5
         XJ   =  (3.0*X(I,J) -4.0*X(I,J-1) +X(I,J-2))*0.5
         YJ   =  (3.0*Y(I,J) -4.0*Y(I,J-1) +Y(I,J-2))*0.5
         VOL = XI*YJ - XJ*YI
          QJ(I,J)  = 1./(VOL + 1.E-25)
         XIX(I,J)  =  YJ*QJ(I,J)
         XIY(I,J)  = -XJ*QJ(I,J)
        ETAX(I,J)  = -YI*QJ(I,J)
        ETAY(I,J)  =  XI*QJ(I,J)
C
      DO 100 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IDA    = IDBC( 2,IBC)
      ISA    = IDBC( 3,IBC)
      ICA    = IDBC( 4,IBC)
      I2S    = IDBC( 5,IBC)
      I2E    = IDBC( 6,IBC)
      IPA    = IDBC( 7,IBC)
      LENA   = (I2E - I2S)*ICA + 1
      IF(IBTYPE .EQ. 2 ) THEN
      DO 110 I2 = 1,LENA
      IAI    = IJP(IPA + LENA*0 + I2)
      IAJ    = IJP(IPA + LENA*1 + I2)
      IBI    = IJP(IPA + LENA*2 + I2)
      IBJ    = IJP(IPA + LENA*3 + I2)
      ICI    = IJP(IPA + LENA*4 + I2)
      ICJ    = IJP(IPA + LENA*5 + I2)
      IDI    = IJP(IPA + LENA*6 + I2)
      IDJ    = IJP(IPA + LENA*7 + I2)
      IEI    = IJP(IPA + LENA*8 + I2)
      IEJ    = IJP(IPA + LENA*9 + I2)
          QJ(IAI,IAJ)  =    QJ(IEI,IEJ)
         XIX(IAI,IAJ)  =  -XIX(IEI,IEJ)
         XIY(IAI,IAJ)  =   XIY(IEI,IEJ)
        ETAX(IAI,IAJ)  =  ETAX(IEI,IEJ)
        ETAY(IAI,IAJ)  = -ETAY(IEI,IEJ)
          QJ(IBI,IBJ)  =    QJ(IDI,IDJ)
         XIX(IBI,IBJ)  =  -XIX(IDI,IDJ)
         XIY(IBI,IBJ)  =   XIY(IDI,IDJ)
        ETAX(IBI,IBJ)  =  ETAX(IDI,IDJ)
        ETAY(IBI,IBJ)  = -ETAY(IDI,IDJ)
         XIX(ICI,ICJ)  =  0.
        ETAY(ICI,ICJ)  =  0.
  110 CONTINUE
      ENDIF
      IF(IBTYPE .EQ. 3) THEN
      DO 115 I2 = 1,LENA
      IAI    = IJP(IPA + LENA*0 + I2)
      IAJ    = IJP(IPA + LENA*1 + I2)
      IBI    = IJP(IPA + LENA*2 + I2)
      IBJ    = IJP(IPA + LENA*3 + I2)
      ICI    = IJP(IPA + LENA*4 + I2)
      ICJ    = IJP(IPA + LENA*5 + I2)
      IDI    = IJP(IPA + LENA*6 + I2)
      IDJ    = IJP(IPA + LENA*7 + I2)
      IEI    = IJP(IPA + LENA*8 + I2)
      IEJ    = IJP(IPA + LENA*9 + I2)
          QJ(IAI,IAJ)  =    QJ(IEI,IEJ)
         XIX(IAI,IAJ)  =   XIX(IEI,IEJ)
         XIY(IAI,IAJ)  =  -XIY(IEI,IEJ)
        ETAX(IAI,IAJ)  = -ETAX(IEI,IEJ)
        ETAY(IAI,IAJ)  =  ETAY(IEI,IEJ)
          QJ(IBI,IBJ)  =    QJ(IDI,IDJ)
         XIX(IBI,IBJ)  =   XIX(IDI,IDJ)
         XIY(IBI,IBJ)  =  -XIY(IDI,IDJ)
        ETAX(IBI,IBJ)  = -ETAX(IDI,IDJ)
        ETAY(IBI,IBJ)  =  ETAY(IDI,IDJ)
         XIY(ICI,ICJ)  =  0.
        ETAX(ICI,ICJ)  =  0.
  115 CONTINUE
      ENDIF
      IF(IBTYPE .EQ. 4) THEN
      IDB    = IDBC(12,IBC)
      ISB    = IDBC(13,IBC)
      ICB    = IDBC(14,IBC)
      I3S    = IDBC(15,IBC)
      I3E    = IDBC(16,IBC)
      IPB    = IDBC( 8,IBC)
      LENB   = (I3E - I3S)*ICB + 1
      DO 120 I2 = 1,LENA
      IAI    = IJP(IPA + LENA*0 + I2)
      IAJ    = IJP(IPA + LENA*1 + I2)
      IBI    = IJP(IPA + LENA*2 + I2)
      IBJ    = IJP(IPA + LENA*3 + I2)
      ICI    = IJP(IPA + LENA*4 + I2)
      ICJ    = IJP(IPA + LENA*5 + I2)
      IDI    = IJP(IPA + LENA*6 + I2)
      IDJ    = IJP(IPA + LENA*7 + I2)
      IEI    = IJP(IPA + LENA*8 + I2)
      IEJ    = IJP(IPA + LENA*9 + I2)
      JAI    = IJP(IPB + LENB*0 + I2)
      JAJ    = IJP(IPB + LENB*1 + I2)
      JBI    = IJP(IPB + LENB*2 + I2)
      JBJ    = IJP(IPB + LENB*3 + I2)
      JCI    = IJP(IPB + LENB*4 + I2)
      JCJ    = IJP(IPB + LENB*5 + I2)
      JDI    = IJP(IPB + LENB*6 + I2)
      JDJ    = IJP(IPB + LENB*7 + I2)
      JEI    = IJP(IPB + LENB*8 + I2)
      JEJ    = IJP(IPB + LENB*9 + I2)
          QJ(IAI,IAJ)  =    QJ(JEI,JEJ)
         XIX(IAI,IAJ)  =   XIX(JEI,JEJ)
         XIY(IAI,IAJ)  =   XIY(JEI,JEJ)
        ETAX(IAI,IAJ)  =  ETAX(JEI,JEJ)
        ETAY(IAI,IAJ)  =  ETAY(JEI,JEJ)
          QJ(IBI,IBJ)  =    QJ(JDI,JDJ)
         XIX(IBI,IBJ)  =   XIX(JDI,JDJ)
         XIY(IBI,IBJ)  =   XIY(JDI,JDJ)
        ETAX(IBI,IBJ)  =  ETAX(JDI,JDJ)
        ETAY(IBI,IBJ)  =  ETAY(JDI,JDJ)
          QJ(JAI,JAJ)  =    QJ(IEI,IEJ)
         XIX(JAI,JAJ)  =   XIX(IEI,IEJ)
         XIY(JAI,JAJ)  =   XIY(IEI,IEJ)
        ETAX(JAI,JAJ)  =  ETAX(IEI,IEJ)
        ETAY(JAI,JAJ)  =  ETAY(IEI,IEJ)
          QJ(JBI,JBJ)  =    QJ(IDI,IDJ)
         XIX(JBI,JBJ)  =   XIX(IDI,IDJ)
         XIY(JBI,JBJ)  =   XIY(IDI,IDJ)
        ETAX(JBI,JBJ)  =  ETAX(IDI,IDJ)
        ETAY(JBI,JBJ)  =  ETAY(IDI,IDJ)
  120 CONTINUE
      ENDIF
  100 CONTINUE
C
c     CALL TWRITE(-1,NI2,-1,NJ2,X,'X         ')
c     CALL TWRITE(-1,NI2,-1,NJ2,Y,'Y         ')
c     CALL TWRITE(-1,NI2,-1,NJ2,QJ,'QJ        ')
c     CALL TWRITE(-1,NI2,-1,NJ2,XIX,'XIX       ')
c     CALL TWRITE(-1,NI2,-1,NJ2,XIY,'XIY       ')
c     CALL TWRITE(-1,NI2,-1,NJ2,ETAX,'ETAX      ')
c     CALL TWRITE(-1,NI2,-1,NJ2,ETAY,'ETAY      ')
C
      DO 200 IJ = 1, MI*MJ
      I          = MOD(IJ-1,MI)-1
      J          = (IJ-1)/MI-1
      IF(QJ(I,J) .LE. 0.) THEN
        WRITE(*,*)'QJ .LE. 0.  AT IJ=',IJ,I,J
        STOP
      ENDIF
  200 CONTINUE
      RETURN
      END
      SUBROUTINE POINT
C
      COMMON /TAPE / INAME,IGRID,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /NEXT / NEXTLOC
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /PNT  /  LX,LY,LXIX,LXIY,LETAX,LETAY,LJAC,
     1                LR,LE,LT,LU,LV,LP,LXM,LPSI,LRNP,LRUNP,LRVNP,LRENP
     1               ,LUS,LVS,LRS,LRUS,LRVS,LRES,LDRS,LDRUS,LDRVS,LDRES
     1               ,LUSP,LUSM,LVSP,LVSM,LRNPO
      COMMON /PNTI / LKDI,LKDJ,LLKD,LKDP,MDA,MDB,NDIG
C...
      NI1     = NI + 1
      NI2     = NI + 2
      NJ1     = NJ + 1
      NJ2     = NJ + 2
C
      MI       = NI + 4
      MJ       = NJ + 4
      MIJ2     = MI*MJ-2
      MXI      = MAX(MI,MJ)
      IS       = -1
      IE       = NI2
      JS       = -1
      JE       = NJ2
      IQS      = IS
      IQE      = 20*(IE-IS+1) + IS - 1
      JQS      = JS
      JQE      = 20*(JE-JS+1) + JS - 1
C
      MXL = MI + MJ
      IF(MXL.GT.MAXL)THEN
       WRITE(*,*)' ERROR : MXL .GT. MAXL (POINT) INSUFF. LU DIAG PLAN'
       WRITE(*,*)'         MXL = ',MXL,'   MAXL=',MAXL
       STOP
      ENDIF
      IF(MXI.GT.MAXI)THEN
       WRITE(*,*)' ERROR : MXI .GT. MAXI (POINT) INSUFF. POINT LENGTH'
       WRITE(*,*)'         MXI = ',MXI,'   MAXI=',MAXI
       STOP
      ENDIF
C
      IS2D   = MI * MJ
      MXG    = IS2D
C
      MDA    = IS2D
      MDB    = MI + MJ
      LKDI   = 1
      LKDJ   = LKDI + MDA
      LLKD   = LKDJ + MDA
      LKDP   = LLKD + MDB
      MXLI   = LKDP + MDB
C
      IF(MXG.GT.MAXG)THEN
       WRITE(*,*)' ERROR : MXG .GT. MAXG (POINT) INSUFF. MAXGRD POINT'
       WRITE(*,*)'         MXG = ',MXG,'    MAXG=',MAXG
       STOP
      ENDIF
      IF(MXLI.GT.MAXLI)THEN
       WRITE(*,*)' ERROR : MXLI.GT.MAXLUI (POINT)INSUF. LU IJ INDEX '
       WRITE(*,*)'         MXLI = ',MXLI,'    MAXLUI=',MAXLI
       STOP
      ENDIF
C
      LX      = NEXTLOC
      LY      = LX     + IS2D
      NEXTLOC = LY     + IS2D
C
      LXIX    = NEXTLOC
      LXIY    = LXIX   + IS2D
      LETAX   = LXIY   + IS2D
      LETAY   = LETAX  + IS2D
      NEXTLOC = LETAY  + IS2D
C
      LJAC    = NEXTLOC
      NEXTLOC = LJAC   + IS2D
C
      LR      = NEXTLOC
      LE      = LR     + IS2D
      LT      = LE     + IS2D
      LU      = LT     + IS2D
      LV      = LU     + IS2D
      LP      = LV     + IS2D
      LXM     = LP     + IS2D
      LPSI    = LXM    + IS2D
      NEXTLOC = LPSI   + IS2D
C
      LRNP    = NEXTLOC
      LRUNP   = LRNP   + IS2D * 3
      LRVNP   = LRUNP  + IS2D * 3
      LRENP   = LRVNP  + IS2D * 3
      NEXTLOC = LRENP  + IS2D * 3
C
      LUS     = NEXTLOC
      LVS     = LUS    + IS2D * 5
      LRS     = LVS    + IS2D * 5
      LRUS    = LRS    + IS2D * 5
      LRVS    = LRUS   + IS2D * 5
      LRES    = LRVS   + IS2D * 5
      LDRS    = LRES   + IS2D * 5
      LRNPO   = LDRS   + IS2D * 5
      NEXTLOC = LRNPO  + IS2D * 4
C
      LDRUS   = NEXTLOC+ 0
      LDRVS   = LDRUS  + 0
      LDRES   = LDRVS  + 0
      NEXTLOC = LDRES  + 0
C
      LUSP    = NEXTLOC
      LUSM    = LUSP   + 0
      LVSP    = LUSM   + 0
      LVSM    = LVSP   + 0
      NEXTLOC = LVSM   + 0
C
      WRITE(ISHOW,*) ' *** WORKING SIZE TOTAL *** '
      WRITE(ISHOW,*) '     NEXTLOC = ' , NEXTLOC
C
      RAT = FLOAT(NEXTLOC)/FLOAT(IS2D)
      WRITE(ISHOW,*) '     RATIO OF NEXLOC/GRID PTS = ',RAT
C
      IF( NEXTLOC .GT. MAXW ) THEN
         WRITE(ISHOW,*) ' *** ERROR WORKING SIZE *** '
         WRITE(ISHOW,*) '     NEXTLOC = ',NEXTLOC,'   MAXW =',MAXW
         STOP
      END IF
C
      RETURN
      END
      SUBROUTINE QPLT3D(IFILE,R,U,V,P,E,XM,X,XIX,XIY,ETAX,ETAY,PSI,
     1                  DRS,Y)
      PARAMETER (MAXIJP=516)
      COMMON /TAPE / INAME,IGRID,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /FINP3/ ALPHA,UINF
      COMMON /CNST4/ DR0,UR0,VR0,PR0,TR0,ER0,
     1               DL0,UL0,VL0,PL0,TL0,EL0
      COMMON /CNST2/ R0,U0,V0,P0,T0,E0
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /PSII/ ICYL,PSI0,AMPLI,NCURVE,URATIO,YY0
      CHARACTER*80 FASN1,FASNR
      CHARACTER CU1,CU2,CU3
C
      DIMENSION    R(-1:NI2,-1:NJ2), DRS(-1:NI2,-1:NJ2,5)
      DIMENSION    U(-1:NI2,-1:NJ2),   V(-1:NI2,-1:NJ2)
      DIMENSION    P(-1:NI2,-1:NJ2),   E(-1:NI2,-1:NJ2)
      DIMENSION   XM(-1:NI2,-1:NJ2),   X(-1:NI2,-1:NJ2)
      DIMENSION  XIX(-1:NI2,-1:NJ2), XIY(-1:NI2,-1:NJ2)
      DIMENSION  ETAX(-1:NI2,-1:NJ2), ETAY(-1:NI2,-1:NJ2)
      DIMENSION  PSI(-1:NI2,-1:NJ2),Y(-1:NI2,-1:NJ2)
C
      GOTO 59
      DO 20 J = 1,NJ
      I = 1
      DRS(I,J,1) = -(3.0*V(I,J) -4.0*V(I+1,J) +V(I+2,J))*0.5*XIX(I,J)
      DO 25 I = 2,NI-1
      DRS(I,J,1) = (V(I+1,J)-V(I-1,J))*0.5*XIX(I,J)
25    CONTINUE
      I = NI
      DRS(I,J,1) = (3.0*V(I,J) -4.0*V(I-1,J) +V(I-2,J))*0.5*XIX(I,J)
20    CONTINUE
      DO 30 I = 1,NI
      J = 1
      DRS(I,J,1) = DRS(I,J,1)-(3.*V(I,J)-4.*V(I,J+1)+V(I,J+2))
     1             *0.5*ETAX(I,J)
      DO 35 J = 2,NJ-1
      DRS(I,J,1) = DRS(I,J,1)+(V(I,J+1)-V(I,J-1))*0.5*ETAX(I,J)
35    CONTINUE
      J = NJ
      DRS(I,J,1) = DRS(I,J,1)+(3.*V(I,J)-4.*V(I,J-1)+V(I,J-2))
     1             *0.5*ETAX(I,J)
30    CONTINUE
      DO 40 J = 1,NJ
      I = 1
      DRS(I,J,1) = DRS(I,J,1)+(3.*U(I,J)-4.*U(I+1,J)+U(I+2,J))
     1             *0.5*XIY(I,J)
      DO 45 I = 2,NI-1
      DRS(I,J,1) = DRS(I,J,1)-(U(I+1,J)-U(I-1,J))*0.5*XIY(I,J)
45    CONTINUE
      I = NI
      DRS(I,J,1) = DRS(I,J,1)-(3.*U(I,J)-4.*U(I-1,J)+U(I-2,J))
     1             *0.5*XIY(I,J)
40    CONTINUE
      DO 50 I = 1,NI
      J = 1
      DRS(I,J,1) = DRS(I,J,1)+(3.*U(I,J)-4.*U(I,J+1)+U(I,J+2))
     1              *0.5*ETAY(I,J)
      DO 55 J = 2,NJ-1
      DRS(I,J,1) = DRS(I,J,1)-(U(I,J+1)-U(I,J-1))*0.5*ETAY(I,J)
55    CONTINUE
      J = NJ
      DRS(I,J,1) = DRS(I,J,1)-(3.*U(I,J)-4.*U(I,J-1)+U(I,J-2))
     1             *0.5*ETAY(I,J)
50    CONTINUE
59    CONTINUE
C
C     psi is a streamline variable for shock over a cylinder
C
      DO 60 I = 1,NI
       PSI(I,1) = 0.
      DO 60 J = 2,NJ
       DX = X(I,J) - X(I,J-1)
       DY = Y(I,J) - Y(I,J-1)
       RAVE = 0.5*( R(I,J) + R(I,J-1) )
       UAVE = 0.5*( U(I,J) + U(I,J-1) )
       VAVE = 0.5*( V(I,J) + V(I,J-1) )
       PSI(I,J) = PSI(I,J-1) + RAVE*( UAVE*DY-VAVE*DX )
60    CONTINUE
C
C     drs is a Mach number variable for shock over a cylinder
C
      DO 70 I = 1,NI
      DO 70 J = 1,NJ
C         AUP   = GAMMA*GM1*(E(I,J)-R(I,J))
C         ADOWN = R(I,J)+GAMMA*(E(I,J)-R(I,J))
         AUP   = GAMMA*ABS(P(I,J))
         ADOWN = ABS(R(I,J))+AUP/(GAMMA-1.)
         UV2   = SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))
         DRS(I,J,1) = UV2/SQRT(AUP/ADOWN)
70    CONTINUE
C
      IF( IFILE .LE. 0 .OR. IFILE .GE. 1000 )THEN
        WRITE(6,*)'UNIT NUMBER ERROR !!!'
        STOP
      ENDIF
      IU1 = IFILE/10
      IF( IU1 .EQ. 0 )THEN
        CU1   = CHAR(IFILE+48)
*        FASN1 = 'assign -F f77 -N ieee u:'//CU1
*        FASNR = 'assign -R u:'//CU1
       ELSEIF( IU1 .GE. 10 .AND. IU1 .LE. 99 )THEN
        IU2 = MOD(IFILE,10)
        CU1 = CHAR(IU2+48)
        IU3 = MOD(IU1,10)
        CU2 = CHAR(IU3+48)
        IU4 = IU1/10
        CU3 = CHAR(IU4+48)
*        FASN1 = 'assign -F f77 -N ieee u:'//CU3//CU2//CU1
*        FASNR = 'assign -R u:'//CU3//CU2//CU1
       ELSEIF( IU1 .GE. 1 .AND. IU1 .LE. 9 )THEN
        CU2 = CHAR(IU1+48)
        IU2 = MOD(IFILE,10)
        CU1 = CHAR(IU2+48)
*        FASN1 = 'assign -F f77 -N ieee u:'//CU2//CU1
*        FASNR = 'assign -R u:'//CU2//CU1
       ENDIF
*      CALL ASSIGN( FASN1 )
c
      WRITE(IFILE)NI,NJ
      WRITE(IFILE)UINF,ALPHA,GAMMA,TIME,ITER
      WRITE(IFILE)(( R(I,J),I=1,NI),J=1,NJ),
     1           (( U(I,J),I=1,NI),J=1,NJ),
     1           (( V(I,J),I=1,NI),J=1,NJ),
     1           (( P(I,J),I=1,NI),J=1,NJ),
     1           (( PSI(I,J),I=1,NI),J=1,NJ),
     1           (( DRS(I,J,1),I=1,NI),J=1,NJ)
*      CALL ASSIGN( FASNR )
C
      IRSTA = 16
*      CALL ASSIGN('assign -F f77 -N ieee u:16')
      WRITE(IRSTA)MI,MJ
      WRITE(IRSTA)UINF,ALPHA,GAMMA,TIME,ITER
      WRITE(IRSTA)(( R(I,J),I=-1,MI-2),J=-1,MJ-2),
     1           (( U(I,J),I=-1,MI-2),J=-1,MJ-2),
     1           (( V(I,J),I=-1,MI-2),J=-1,MJ-2),
     1           (( P(I,J),I=-1,MI-2),J=-1,MJ-2),
     1           (( PSI(I,J),I=-1,MI-2),J=-1,MJ-2)
*      CALL ASSIGN('assign -R u:16')
C
       IF( ICYL .NE. 0 .AND. MOD(ITER,10) .EQ. 0)THEN
         WRITE(15,99)TIME,P(1,1),P(61,1),P(121,1),P(181,1),P(241,1),
     1               P(301,1),P(361,1)
       ENDIF
99     FORMAT(1X,8(F9.4,1X))
C
c     CLOSE(15)
      CLOSE (IFILE)
      CLOSE (IRSTA)
      RETURN
      END
      SUBROUTINE RSTQ(R,U,V,P,E,PSI,ITER0,TIME1)
      COMMON /TAPE / INAME,IGRID,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /CNST2/ R0,U0,V0,P0,T0,E0
      COMMON /CNST4/ DR0,UR0,VR0,PR0,TR0,ER0,
     1               DL0,UL0,VL0,PL0,TL0,EL0
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      DIMENSION    R(-1:NI2,-1:NJ2),   E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2),   V(-1:NI2,-1:NJ2)
      DIMENSION    P(-1:NI2,-1:NJ2),   T(-1:NI2,-1:NJ2)
      DIMENSION   XM(-1:NI2,-1:NJ2), PSI(-1:NI2,-1:NJ2)
C
      READ(IRSTA)MI,MJ
C
      IF(MI .NE. NI+4 .OR. MJ .NE. NJ+4) THEN
        WRITE(*,*)'ERROR : RSTQ,  MI OR MJ NOT MATCH'
        STOP
      ENDIF
C
      READ(IRSTA)XMINF,ALPHA,GAMMA,TIME1,ITER0
      READ(IRSTA)(( R(I,J),I=-1,MI-2),J=-1,MJ-2),
     1             (( U(I,J),I=-1,MI-2),J=-1,MJ-2),
     1             (( V(I,J),I=-1,MI-2),J=-1,MJ-2),
     1             (( P(I,J),I=-1,MI-2),J=-1,MJ-2),
     1             (( PSI(I,J),I=-1,MI-2),J=-1,MJ-2)
C
      DO 10 IJ  =-1,MIJ2
      E(IJ,-1)     = P(IJ,-1)/GM1 + R(IJ,-1)
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE RUNGE2(IJP,MLU,X,Y,R,E,U,V,T,P,QJ,XIX,XIY,ETAX,ETAY,
     1 RS,RUS,RVS,RES,US,VS,DRS,RNP,RUNP,RVNP,RENP,RNPO)
C
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /RES01/ RDMX,IDMX,JDMX,RTMX,ITMX,JTMX,RDL2,RTL2
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
C
      DIMENSION IJP(MXBP),MLU(MAXLI)
      DIMENSION    X(-1:NI2,-1:NJ2),     Y(-1:NI2,-1:NJ2)
      DIMENSION   QJ(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2), E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2), V(-1:NI2,-1:NJ2)
      DIMENSION    T(-1:NI2,-1:NJ2), P(-1:NI2,-1:NJ2)
      DIMENSION  DRS(-1:NI2,-1:NJ2,5)
      DIMENSION  RUS(-1:NI2,-1:NJ2,5), RVS(-1:NI2,-1:NJ2,5)
      DIMENSION   RS(-1:NI2,-1:NJ2,5), RES(-1:NI2,-1:NJ2,5)
      DIMENSION   US(-1:NI2,-1:NJ2,5),  VS(-1:NI2,-1:NJ2,5)
      DIMENSION  RNP(-1:NI2,-1:NJ2,3),RUNP(-1:NI2,-1:NJ2,3)
      DIMENSION RVNP(-1:NI2,-1:NJ2,3),RENP(-1:NI2,-1:NJ2,3)
      DIMENSION  XIX(-1:NI2,-1:NJ2),XIY(-1:NI2,-1:NJ2)
      DIMENSION ETAX(-1:NI2,-1:NJ2),ETAY(-1:NI2,-1:NJ2)
      DIMENSION RNPO(-1:NI2,-1:NJ2,4)
C
      DO 5 IJ = -1,MIJ2
       DRS(IJ,-1,2) = RS(IJ,-1,1)+RS(IJ,-1,2)+RS(IJ,-1,3)+
     1                RS(IJ,-1,4)+RS(IJ,-1,5)
       DRS(IJ,-1,3) = RUS(IJ,-1,1)+RUS(IJ,-1,2)+RUS(IJ,-1,3)+
     1                RUS(IJ,-1,4)+RUS(IJ,-1,5)
       DRS(IJ,-1,4) = RVS(IJ,-1,1)+RVS(IJ,-1,2)+RVS(IJ,-1,3)+
     1                RVS(IJ,-1,4)+RVS(IJ,-1,5)
       DRS(IJ,-1,5) = RES(IJ,-1,1)+RES(IJ,-1,2)+RES(IJ,-1,3)+
     1                RES(IJ,-1,4)+RES(IJ,-1,5)
5     CONTINUE
      DO 10 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
10    CONTINUE
C
      IF( MTHR .EQ. 9 )THEN
       CALL FIWEN2( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJWEN2( VS,QJ, RS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 11 )THEN
       CALL FIEEN2( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJEEN2( VS,QJ, RS, DRS(-1,-1,1) )
      ENDIF
      DO 15 IJ = -1,MIJ2
       RNP(IJ,-1,1) = QJ(IJ,-1)*(DRS(IJ,-1,2)+DRS(IJ,-1,1))
15    CONTINUE
      DO 20 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
20    CONTINUE
C
      IF( MTHR .EQ. 9 )THEN
       CALL FIWEN2( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJWEN2( VS,QJ, RUS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 11 )THEN
       CALL FIEEN2( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJEEN2( VS,QJ, RUS, DRS(-1,-1,1) )
      ENDIF
      DO 25 IJ = -1,MIJ2
       RUNP(IJ,-1,1) = QJ(IJ,-1)*(DRS(IJ,-1,3)+DRS(IJ,-1,1))
25    CONTINUE
      DO 30 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
30    CONTINUE
C
      IF( MTHR .EQ. 9 )THEN
       CALL FIWEN2( US,QJ, RVS, DRS(-1,-1,1) )
       CALL FJWEN2( VS,QJ, RVS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 11 )THEN
       CALL FIEEN2( US,QJ, RVS, DRS(-1,-1,1) )
       CALL FJEEN2( VS,QJ, RVS, DRS(-1,-1,1) )
      ENDIF
      DO 35 IJ = -1,MIJ2
       RVNP(IJ,-1,1) = QJ(IJ,-1)*(DRS(IJ,-1,4)+DRS(IJ,-1,1))
35    CONTINUE
      DO 40 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
40    CONTINUE
C
      IF( MTHR .EQ. 9 )THEN
       CALL FIWEN2( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJWEN2( VS,QJ, RES, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 11 )THEN
       CALL FIEEN2( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJEEN2( VS,QJ, RES, DRS(-1,-1,1) )
      ENDIF
      DO 45 IJ = -1,MIJ2
       RENP(IJ,-1,1) = QJ(IJ,-1)*(DRS(IJ,-1,5)+DRS(IJ,-1,1))
45    CONTINUE
C
      CALL SOLUV(R,E,U,V,P,RNP(-1,-1,1),RUNP(-1,-1,1),RVNP(-1,-1,1),
     1           RENP(-1,-1,1))
      CALL BC(IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P,X)
      CALL BEAM5(R,E,RS,RUS,RVS,RES,U,V,P,US,VS,
     1           XIX,XIY,ETAX,ETAY,QJ)
      DO 50 IJ = -1,MIJ2
       RNP(IJ,-1,1)  = RS(IJ,-1,1)+RS(IJ,-1,2)+RS(IJ,-1,3)+
     1                 RS(IJ,-1,4)+RS(IJ,-1,5)
       RUNP(IJ,-1,1) = RUS(IJ,-1,1)+RUS(IJ,-1,2)+RUS(IJ,-1,3)+
     1                 RUS(IJ,-1,4)+RUS(IJ,-1,5)
       RVNP(IJ,-1,1) = RVS(IJ,-1,1)+RVS(IJ,-1,2)+RVS(IJ,-1,3)+
     1                 RVS(IJ,-1,4)+RVS(IJ,-1,5)
       RENP(IJ,-1,1) = RES(IJ,-1,1)+RES(IJ,-1,2)+RES(IJ,-1,3)+
     1                 RES(IJ,-1,4)+RES(IJ,-1,5)
50    CONTINUE
C
      DO 55 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
55    CONTINUE
C
      IF( MTHR .EQ. 9 )THEN
       CALL FIWEN2( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJWEN2( VS,QJ, RS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 11 )THEN
       CALL FIEEN2( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJEEN2( VS,QJ, RS, DRS(-1,-1,1) )
      ENDIF
      DO 60 IJ = -1,MIJ2
       RNP(IJ,-1,2) = 0.25*QJ(IJ,-1)*(3.*DRS(IJ,-1,2)+RNP(IJ,-1,1)+
     1                DRS(IJ,-1,1))
60    CONTINUE
      DO 65 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
65    CONTINUE
C
      IF( MTHR .EQ. 9 )THEN
       CALL FIWEN2( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJWEN2( VS,QJ, RUS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 11 )THEN
       CALL FIEEN2( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJEEN2( VS,QJ, RUS, DRS(-1,-1,1) )
      ENDIF
      DO 70 IJ = -1,MIJ2
       RUNP(IJ,-1,2) = 0.25*QJ(IJ,-1)*(3.*DRS(IJ,-1,3)+RUNP(IJ,-1,1)+
     1                 DRS(IJ,-1,1))
70    CONTINUE
      DO 75 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
75    CONTINUE
C
      IF( MTHR .EQ. 9 )THEN
       CALL FIWEN2( US,QJ, RVS, DRS(-1,-1,1) )
       CALL FJWEN2( VS,QJ, RVS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 11 )THEN
       CALL FIEEN2( US,QJ, RVS, DRS(-1,-1,1) )
       CALL FJEEN2( VS,QJ, RVS, DRS(-1,-1,1) )
      ENDIF
      DO 80 IJ = -1,MIJ2
       RVNP(IJ,-1,2)=0.25*QJ(IJ,-1)*(3.*DRS(IJ,-1,4)+RVNP(IJ,-1,1)+
     1               DRS(IJ,-1,1))
80    CONTINUE
      DO 85 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
85    CONTINUE
C
      IF( MTHR .EQ. 9 )THEN
       CALL FIWEN2( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJWEN2( VS,QJ, RES, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 11 )THEN
       CALL FIEEN2( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJEEN2( VS,QJ, RES, DRS(-1,-1,1) )
      ENDIF
      DO 90 IJ = -1,MIJ2
       RENP(IJ,-1,2)=0.25*QJ(IJ,-1)*(3.*DRS(IJ,-1,5)+RENP(IJ,-1,1)+
     1               DRS(IJ,-1,1))
90    CONTINUE
C
      CALL SOLUV(R,E,U,V,P,RNP(-1,-1,2),RUNP(-1,-1,2),RVNP(-1,-1,2),
     1           RENP(-1,-1,2))
      CALL BC(IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P,X)
      CALL BEAM5(R,E,RS,RUS,RVS,RES,U,V,P,US,VS,
     1           XIX,XIY,ETAX,ETAY,QJ)
C
      DO 95 IJ = -1,MIJ2
       RNP(IJ,-1,2)  = RS(IJ,-1,1)+RS(IJ,-1,2)+RS(IJ,-1,3)+
     1                 RS(IJ,-1,4)+RS(IJ,-1,5)
       RUNP(IJ,-1,2) = RUS(IJ,-1,1)+RUS(IJ,-1,2)+RUS(IJ,-1,3)+
     1                 RUS(IJ,-1,4)+RUS(IJ,-1,5)
       RVNP(IJ,-1,2) = RVS(IJ,-1,1)+RVS(IJ,-1,2)+RVS(IJ,-1,3)+
     1                 RVS(IJ,-1,4)+RVS(IJ,-1,5)
       RENP(IJ,-1,2) = RES(IJ,-1,1)+RES(IJ,-1,2)+RES(IJ,-1,3)+
     1                 RES(IJ,-1,4)+RES(IJ,-1,5)
95    CONTINUE
      DO 100 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
100   CONTINUE
C
      IF( MTHR .EQ. 9 )THEN
       CALL FIWEN2( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJWEN2( VS,QJ, RS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 11 )THEN
       CALL FIEEN2( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJEEN2( VS,QJ, RS, DRS(-1,-1,1) )
      ENDIF
      DO 105 IJ = -1,MIJ2
       RNPO(IJ,-1,1) = 1./3.*QJ(IJ,-1)*(DRS(IJ,-1,2)+2.*RNP(IJ,-1,2)+
     1                 2.*DRS(IJ,-1,1))
105   CONTINUE
      DO 110 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
110   CONTINUE
C
      IF( MTHR .EQ. 9 )THEN
       CALL FIWEN2( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJWEN2( VS,QJ, RUS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 11 )THEN
        CALL FIEEN2( US,QJ, RUS, DRS(-1,-1,1) )
        CALL FJEEN2( VS,QJ, RUS, DRS(-1,-1,1) )
      ENDIF
      DO 115 IJ = -1,MIJ2
       RNPO(IJ,-1,2) = 1./3.*QJ(IJ,-1)*(DRS(IJ,-1,3)+2.*RUNP(IJ,-1,2)+
     1                 2.*DRS(IJ,-1,1))
115   CONTINUE
      DO 120 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
120   CONTINUE
C
      IF( MTHR .EQ. 9 )THEN
       CALL FIWEN2( US,QJ, RVS, DRS(-1,-1,1) )
       CALL FJWEN2( VS,QJ, RVS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 11 )THEN
        CALL FIEEN2( US,QJ, RVS, DRS(-1,-1,1) )
        CALL FJEEN2( VS,QJ, RVS, DRS(-1,-1,1) )
      ENDIF
      DO 125 IJ = -1,MIJ2
       RNPO(IJ,-1,3) = 1./3.*QJ(IJ,-1)*(DRS(IJ,-1,4)+2.*RVNP(IJ,-1,2)+
     1                 2.*DRS(IJ,-1,1))
125     CONTINUE
      DO 130 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
130   CONTINUE
C
      IF( MTHR .EQ. 9 )THEN
       CALL FIWEN2( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJWEN2( VS,QJ, RES, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 11 )THEN
        CALL FIEEN2( US,QJ, RES, DRS(-1,-1,1) )
        CALL FJEEN2( VS,QJ, RES, DRS(-1,-1,1) )
      ENDIF
      DO 135 IJ = -1,MIJ2
       RNPO(IJ,-1,4) = 1./3.*QJ(IJ,-1)*(DRS(IJ,-1,5)+2.*RENP(IJ,-1,2)+
     1                 2.*DRS(IJ,-1,1))
135   CONTINUE
C
      SUMD         = 0.
      SUMT         = 0.
      RDMX         = -1.
      RTMX         = -1.
      DO 200 I     = 2,NI-1
      DO 200 J     = 2,NJ-1
      RNP1  = RNPO(I,J,1)
      RUNP1 = RNPO(I,J,2)
      RVNP1 = RNPO(I,J,3)
      RENP1 = RNPO(I,J,4)
       XU           = 1./SQRT(1.-U(I,J)**2-V(I,J)**2)
       XU2          = XU*XU
       EP           = E(I,J) + P(I,J)
       DSQD         = (RNP1 - R(I,J)*XU)**2
       DSQT         = DSQD + (RUNP1 - XU2*EP*U(I,J))**2
     1                     + (RVNP1 - XU2*EP*V(I,J))**2
     1                     + (RENP1 - XU2*EP+P(I,J))**2
       SUMD         = SUMD + DSQD
       SUMT         = SUMT + DSQT
       IF(SQRT(DSQD) .GT. RDMX)THEN
         RDMX       = SQRT(DSQD)
         IDMX       = I
         JDMX       = J
       ENDIF
       IF(SQRT(DSQT) .GT. RTMX)THEN
         RTMX       = SQRT(DSQT)
         ITMX       = I
         JTMX       = J
       ENDIF
200   CONTINUE
      RDL2         = SQRT(SUMD)
      RTL2         = SQRT(SUMT)
C
      CALL SOLUV(R,E,U,V,P,RNPO(-1,-1,1),RNPO(-1,-1,2),RNPO(-1,-1,3),
     1           RNPO(-1,-1,4))
C
      RETURN
      END
C
      SUBROUTINE RUNGE3(IJP,MLU,X,Y,R,E,U,V,T,P,QJ,XIX,XIY,ETAX,ETAY,
     1 RS,RUS,RVS,RES,US,VS,DRS,RNP,RUNP,RVNP,RENP,RNPO)
C
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /RES01/ RDMX,IDMX,JDMX,RTMX,ITMX,JTMX,RDL2,RTL2
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
C
      DIMENSION IJP(MXBP),MLU(MAXLI)
      DIMENSION    X(-1:NI2,-1:NJ2),     Y(-1:NI2,-1:NJ2)
      DIMENSION   QJ(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2), E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2), V(-1:NI2,-1:NJ2)
      DIMENSION    T(-1:NI2,-1:NJ2), P(-1:NI2,-1:NJ2)
      DIMENSION  DRS(-1:NI2,-1:NJ2,5)
      DIMENSION  RUS(-1:NI2,-1:NJ2,5), RVS(-1:NI2,-1:NJ2,5)
      DIMENSION   RS(-1:NI2,-1:NJ2,5), RES(-1:NI2,-1:NJ2,5)
      DIMENSION   US(-1:NI2,-1:NJ2,5),  VS(-1:NI2,-1:NJ2,5)
      DIMENSION  RNP(-1:NI2,-1:NJ2,3),RUNP(-1:NI2,-1:NJ2,3)
      DIMENSION RVNP(-1:NI2,-1:NJ2,3),RENP(-1:NI2,-1:NJ2,3)
      DIMENSION  XIX(-1:NI2,-1:NJ2),XIY(-1:NI2,-1:NJ2)
      DIMENSION ETAX(-1:NI2,-1:NJ2),ETAY(-1:NI2,-1:NJ2)
      DIMENSION RNPO(-1:NI2,-1:NJ2,4)
C
      DO 5 IJ = -1,MIJ2
       DRS(IJ,-1,2) = RS(IJ,-1,1)+RS(IJ,-1,2)+RS(IJ,-1,3)+
     1                RS(IJ,-1,4)+RS(IJ,-1,5)
       DRS(IJ,-1,3) = RUS(IJ,-1,1)+RUS(IJ,-1,2)+RUS(IJ,-1,3)+
     1                RUS(IJ,-1,4)+RUS(IJ,-1,5)
       DRS(IJ,-1,4) = RVS(IJ,-1,1)+RVS(IJ,-1,2)+RVS(IJ,-1,3)+
     1                RVS(IJ,-1,4)+RVS(IJ,-1,5)
       DRS(IJ,-1,5) = RES(IJ,-1,1)+RES(IJ,-1,2)+RES(IJ,-1,3)+
     1                RES(IJ,-1,4)+RES(IJ,-1,5)
5     CONTINUE
      DO 10 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
10    CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RS, DRS(-1,-1,1) )
      ENDIF
      DO 15 IJ = -1,MIJ2
       RNP(IJ,-1,1) = QJ(IJ,-1)*(DRS(IJ,-1,2)+0.5*DRS(IJ,-1,1))
15    CONTINUE
      DO 20 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
20    CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RUS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RUS, DRS(-1,-1,1) )
      ENDIF
      DO 25 IJ = -1,MIJ2
       RUNP(IJ,-1,1) = QJ(IJ,-1)*(DRS(IJ,-1,3)+0.5*DRS(IJ,-1,1))
25    CONTINUE
      DO 30 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
30    CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RVS, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RVS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
        CALL FIEEN3( US,QJ, RVS, DRS(-1,-1,1) )
        CALL FJEEN3( VS,QJ, RVS, DRS(-1,-1,1) )
      ENDIF
      DO 35 IJ = -1,MIJ2
       RVNP(IJ,-1,1) = QJ(IJ,-1)*(DRS(IJ,-1,4)+0.5*DRS(IJ,-1,1))
35    CONTINUE
      DO 40 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
40    CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RES, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RES, DRS(-1,-1,1) )
      ENDIF
      DO 45 IJ = -1,MIJ2
       RENP(IJ,-1,1) = QJ(IJ,-1)*(DRS(IJ,-1,5)+0.5*DRS(IJ,-1,1))
45    CONTINUE
C
      CALL SOLUV(R,E,U,V,P,RNP(-1,-1,1),RUNP(-1,-1,1),RVNP(-1,-1,1),
     1           RENP(-1,-1,1))
      CALL BC(IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P,X)
      CALL BEAM5(R,E,RS,RUS,RVS,RES,U,V,P,US,VS,
     1           XIX,XIY,ETAX,ETAY,QJ)
      DO 50 IJ = -1,MIJ2
       RNP(IJ,-1,1)  = RS(IJ,-1,1)+RS(IJ,-1,2)+RS(IJ,-1,3)+
     1                 RS(IJ,-1,4)+RS(IJ,-1,5)
       RUNP(IJ,-1,1) = RUS(IJ,-1,1)+RUS(IJ,-1,2)+RUS(IJ,-1,3)+
     1                 RUS(IJ,-1,4)+RUS(IJ,-1,5)
       RVNP(IJ,-1,1) = RVS(IJ,-1,1)+RVS(IJ,-1,2)+RVS(IJ,-1,3)+
     1                 RVS(IJ,-1,4)+RVS(IJ,-1,5)
       RENP(IJ,-1,1) = RES(IJ,-1,1)+RES(IJ,-1,2)+RES(IJ,-1,3)+
     1                 RES(IJ,-1,4)+RES(IJ,-1,5)
50    CONTINUE
C
      DO 55 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
55    CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RS, DRS(-1,-1,1) )
      ENDIF
      DO 60 IJ = -1,MIJ2
       RNP(IJ,-1,2) = QJ(IJ,-1)*(DRS(IJ,-1,2)+0.5*DRS(IJ,-1,1))
60    CONTINUE
      DO 65 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
65    CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RUS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RUS, DRS(-1,-1,1) )
      ENDIF
      DO 70 IJ = -1,MIJ2
       RUNP(IJ,-1,2) = QJ(IJ,-1)*(DRS(IJ,-1,3)+0.5*DRS(IJ,-1,1))
70    CONTINUE
      DO 75 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
75    CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RVS, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RVS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RVS, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RVS, DRS(-1,-1,1) )
      ENDIF
      DO 80 IJ = -1,MIJ2
       RVNP(IJ,-1,2) = QJ(IJ,-1)*(DRS(IJ,-1,4)+0.5*DRS(IJ,-1,1))
80    CONTINUE
      DO 85 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
85    CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RES, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RES, DRS(-1,-1,1) )
      ENDIF
      DO 90 IJ = -1,MIJ2
       RENP(IJ,-1,2) = QJ(IJ,-1)*(DRS(IJ,-1,5)+0.5*DRS(IJ,-1,1))
90    CONTINUE
      DO 95 IJ = -1,MIJ2
       DRS(IJ,-1,1) = 0.
95    CONTINUE
C
      CALL SOLUV(R,E,U,V,P,RNP(-1,-1,2),RUNP(-1,-1,2),RVNP(-1,-1,2),
     1           RENP(-1,-1,2))
      CALL BC(IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P,X)
      CALL BEAM5(R,E,RS,RUS,RVS,RES,U,V,P,US,VS,
     1           XIX,XIY,ETAX,ETAY,QJ)
C
      DO 100 IJ = -1,MIJ2
       RNP(IJ,-1,2)  = RS(IJ,-1,1)+RS(IJ,-1,2)+RS(IJ,-1,3)+
     1                 RS(IJ,-1,4)+RS(IJ,-1,5)
       RUNP(IJ,-1,2) = RUS(IJ,-1,1)+RUS(IJ,-1,2)+RUS(IJ,-1,3)+
     1                 RUS(IJ,-1,4)+RUS(IJ,-1,5)
       RVNP(IJ,-1,2) = RVS(IJ,-1,1)+RVS(IJ,-1,2)+RVS(IJ,-1,3)+
     1                 RVS(IJ,-1,4)+RVS(IJ,-1,5)
       RENP(IJ,-1,2) = RES(IJ,-1,1)+RES(IJ,-1,2)+RES(IJ,-1,3)+
     1                 RES(IJ,-1,4)+RES(IJ,-1,5)
100   CONTINUE
      DO 105 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
105   CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RS, DRS(-1,-1,1) )
      ENDIF
      DO 110 IJ = -1,MIJ2
       RNP(IJ,-1,3) = QJ(IJ,-1)*(DRS(IJ,-1,2)+DRS(IJ,-1,1))
110   CONTINUE
      DO 115 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
115   CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RUS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RUS, DRS(-1,-1,1) )
      ENDIF
      DO 120 IJ = -1,MIJ2
       RUNP(IJ,-1,3) = QJ(IJ,-1)*(DRS(IJ,-1,3)+DRS(IJ,-1,1))
120   CONTINUE
      DO 125 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
125   CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RVS, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RVS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RVS, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RVS, DRS(-1,-1,1) )
      ENDIF
      DO 130 IJ = -1,MIJ2
       RVNP(IJ,-1,3) = QJ(IJ,-1)*(DRS(IJ,-1,4)+DRS(IJ,-1,1))
130   CONTINUE
      DO 135 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
135   CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RES, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RES, DRS(-1,-1,1) )
      ENDIF
      DO 140 IJ = -1,MIJ2
       RENP(IJ,-1,3) = QJ(IJ,-1)*(DRS(IJ,-1,5)+DRS(IJ,-1,1))
140   CONTINUE
C
      CALL SOLUV(R,E,U,V,P,RNP(-1,-1,3),RUNP(-1,-1,3),RVNP(-1,-1,3),
     1           RENP(-1,-1,3))
      CALL BC(IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P,X)
      CALL BEAM5(R,E,RS,RUS,RVS,RES,U,V,P,US,VS,
     1           XIX,XIY,ETAX,ETAY,QJ)
C
      DO 145 IJ = -1,MIJ2
       RNP(IJ,-1,3)  = RS(IJ,-1,1)+RS(IJ,-1,2)+RS(IJ,-1,3)+
     1                 RS(IJ,-1,4)+RS(IJ,-1,5)
       RUNP(IJ,-1,3) = RUS(IJ,-1,1)+RUS(IJ,-1,2)+RUS(IJ,-1,3)+
     1                 RUS(IJ,-1,4)+RUS(IJ,-1,5)
       RVNP(IJ,-1,3) = RVS(IJ,-1,1)+RVS(IJ,-1,2)+RVS(IJ,-1,3)+
     1                 RVS(IJ,-1,4)+RVS(IJ,-1,5)
       RENP(IJ,-1,3) = RES(IJ,-1,1)+RES(IJ,-1,2)+RES(IJ,-1,3)+
     1                 RES(IJ,-1,4)+RES(IJ,-1,5)
145   CONTINUE
      DO 150 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
150   CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RS, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RS, DRS(-1,-1,1) )
      ENDIF
      DO 155 IJ = -1,MIJ2
      RNPO(IJ,-1,1)=1./6.*QJ(IJ,-1)*(-2.*DRS(IJ,-1,2)+2.*RNP(IJ,-1,1)+
     1              4.*RNP(IJ,-1,2)+2.*RNP(IJ,-1,3)+DRS(IJ,-1,1))
155   CONTINUE
      DO 160 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
160   CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RUS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RUS, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RUS, DRS(-1,-1,1) )
      ENDIF
      DO 165 IJ = -1,MIJ2
      RNPO(IJ,-1,2)=1./6.*QJ(IJ,-1)*(-2.*DRS(IJ,-1,3)+2.*RUNP(IJ,-1,1)+
     1              4.*RUNP(IJ,-1,2)+2.*RUNP(IJ,-1,3)+DRS(IJ,-1,1))
165   CONTINUE
      DO 170 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
170   CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RVS, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RVS, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RVS, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RVS, DRS(-1,-1,1) )
      ENDIF
      DO 175 IJ = -1,MIJ2
      RNPO(IJ,-1,3)=1./6.*QJ(IJ,-1)*(-2.*DRS(IJ,-1,4)+2.*RVNP(IJ,-1,1)+
     1              4.*RVNP(IJ,-1,2)+2.*RVNP(IJ,-1,3)+DRS(IJ,-1,1))
175   CONTINUE
      DO 180 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
180   CONTINUE
C
      IF( MTHR .EQ. 10 )THEN
       CALL FIWEN3( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJWEN3( VS,QJ, RES, DRS(-1,-1,1) )
      ELSEIF( MTHR .EQ. 12 )THEN
       CALL FIEEN3( US,QJ, RES, DRS(-1,-1,1) )
       CALL FJEEN3( VS,QJ, RES, DRS(-1,-1,1) )
      ENDIF
      DO 185 IJ = -1,MIJ2
      RNPO(IJ,-1,4)=1./6.*QJ(IJ,-1)*(-2.*DRS(IJ,-1,5)+2.*RENP(IJ,-1,1)+
     1              4.*RENP(IJ,-1,2)+2.*RENP(IJ,-1,3)+DRS(IJ,-1,1))
185   CONTINUE
      DO 190 IJ      = -1,MIJ2
       DRS(IJ,-1,1) = 0.
190   CONTINUE
C
      SUMD         = 0.
      SUMT         = 0.
      RDMX         = -1.
      RTMX         = -1.
      DO 200 I     = 2,NI-1
      DO 200 J     = 2,NJ-1
      RNP1  = RNPO(I,J,1)
      RUNP1 = RNPO(I,J,2)
      RVNP1 = RNPO(I,J,3)
      RENP1 = RNPO(I,J,4)
       XU           = 1./SQRT(ABS(1.-U(I,J)**2-V(I,J)**2))
       XU2          = XU*XU
       EP           = E(I,J) + P(I,J)
       DSQD         = (RNP1 - R(I,J)*XU)**2
       DSQT         = DSQD + (RUNP1 - XU2*EP*U(I,J))**2
     1                     + (RVNP1 - XU2*EP*V(I,J))**2
     1                     + (RENP1 - XU2*EP+P(I,J))**2
       SUMD         = SUMD + DSQD
       SUMT         = SUMT + DSQT
       IF(SQRT(DSQD) .GT. RDMX)THEN
         RDMX       = SQRT(DSQD)
         IDMX       = I
         JDMX       = J
       ENDIF
       IF(SQRT(DSQT) .GT. RTMX)THEN
         RTMX       = SQRT(DSQT)
         ITMX       = I
         JTMX       = J
       ENDIF
200   CONTINUE
      RDL2         = SQRT(SUMD)
      RTL2         = SQRT(SUMT)
C
      CALL SOLUV(R,E,U,V,P,RNPO(-1,-1,1),RNPO(-1,-1,2),RNPO(-1,-1,3),
     1           RNPO(-1,-1,4))
C
      RETURN
      END
      SUBROUTINE SETUP(IJP)
      CHARACTER*80 BCNAME
      COMMON /TAPE / INAME,IGRID,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /NEXT / NEXTLOC
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
      COMMON /FINP3/ ALPHA,UINF
      COMMON /IKTT/ KTT
      COMMON /CNST2/ R0,U0,V0,P0,T0,E0
      COMMON /CNST4/ DR0,UR0,VR0,PR0,TR0,ER0,
     1               DL0,UL0,VL0,PL0,TL0,EL0
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /BC01 / IDBC(30,20),FDBC(10,20)
      COMMON /AXIS / IAXIS
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /PSII/ ICYL,PSI0,AMPLI,NCURVE,URATIO,YY0
      COMMON /JSY/ JSYM,JTHR
      COMMON /NEWTON/ERRN
      COMMON /T6/ AXISYM
      LOGICAL AXISYM
C..
      DIMENSION KAS(2),KAE(2),KBS(2),KBE(2),IDUM(2)
      DIMENSION IJP(MAXBP)
      DIMENSION NDIM(2)
C
      INTEGER CY(2),DD(2,2)
      DATA  CY / 2 , 1 /
      DATA  DD / 1 , 0 , 0 , 1 /
C...
C     NAMELIST/PINP/IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILOCDT,
C    1              CFL1,CFL2,NUP,DTFIX,IAXIS,ENTEPS
C     NAMELIST/FINP/XMINF,ALPHA,GAMMA
C
C     NAMELIST/BCIN/BCNAME,IBTYPE,KAS,KAE,KBS,KBE
C...
C     DATA  IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILOCDT, CFL1, CFL2,
C    1         NUP,DTFIX, IAXIS,ENTEPS
C    1     /     0,   10,    10,   0,   1,     0,  0.9,  0.9,
C    1         100, -0.1,   0, 0.05/
C     DATA  XMINF,ALPHA,GAMMA
C    1     /  1.8,   0.,  1.666667/
C     DATA  BCNAME,IBTYPE,      KAS,      KAE,    KBS,    KBE
C    1     /'NULL',     1,   1,   1,    1,  1,  1,  1,  1,  1/
C...
C
C     READ(INAME,NML=PINP)              ! work in SGI & CONVEX
C     READ(INAME,NML=FINP)              ! work in SGI & CONVEX
C
C     READ(INAME,PINP)                  ! work in CRAY
C     READ(INAME,FINP)                  ! work in CRAY
C...
      READ(12,*)    IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT,
     1              CFL1,CFL2,NUP,DTFIX,ENTEPS,ERRN
      READ(12,*)    UINF,ALPHA,GAMMA,LSTD
      READ(12,*)    ICYL,PSI0,AMPLI,NCURVE,URATIO,YY0,JSYM,JTHR
      IF( LSTD .EQ. 0 )THEN
       READ(12,*)    XST,MTIM,KTT
       READ(12,*)    (TPP(I),I=1,MTIM)
      ENDIF
C
      PI     = 3.1415926535
      SPI    = SQRT(3.1415926535)
c     UTX    = SQRT(2./GAMMA)
      UTX    = 1.
      EPS    = ENTEPS
      EPSS   = EPS*EPS
      EPS2   = 2.*EPS
      GM1    = GAMMA - 1.
      RGM1   = 1./GM1
      RGM2   = 2.*RGM1
C
      IF( LSTD .EQ. 1 )THEN
       R0   = 1.
       T0   = 1.
       P0   = R0*T0
	 U0   = XMINF*COS(ALPHA*PI/180)/UTX
       V0   = XMINF*SIN(ALPHA*PI/180)/UTX
C       U0   = XMINF*COS(ALPHA/57.3)/UTX
C       V0   = XMINF*SIN(ALPHA/57.3)/UTX
       E0   = P0/GM1 + R0
      ELSE
      IF( ICYL .EQ. 0 )THEN
       XMINF= 0.
       XMS  = XMINF
       R0   = 1.
       T0   = 1.
c      P0   = R0*T0
       P0   = 1.
c      U0   = XMINF*COS(ALPHA/57.3)/UTX
       U0   = UINF
c      V0   = XMINF*SIN(ALPHA/57.3)/UTX
       V0   = 0.
c      E0   = 0.5*(T0/(GAMMA-1.)+U0*U0+V0*V0)
       E0   = P0/GM1 + R0
      ELSEIF( ICYL .EQ. 1 )THEN
       UR0    = 0.
       VR0    = 0.
       PR0    = 1.
       DR0    = 1.
       ER0    = PR0/(GAMMA-1.) + DR0
       VL0    = 0.
       READ(12,*)UL0,PL0,DL0
       EL0    = PL0/(GAMMA-1.) + DL0
       READ(12,*) AXISYM
      ENDIF
      ENDIF
C
      WRITE(*,99)
   99 FORMAT(80('-'))
      WRITE(*,*) '    FREE STREAM STATE'
      IF( LSTD .EQ. 0 )THEN
       IF( ICYL .EQ. 1 )THEN
       WRITE(*,*) 'SHOCK MACH NO.              =  ',UINF
       WRITE(*,*) 'ANGLE OF ATTACK             =  ',ALPHA
       WRITE(*,*) 'DENSITY              =  ',DL0,DR0
       WRITE(*,*) 'VELOCITY U           =  ',UL0,UR0
       WRITE(*,*) 'VELOCITY V           =  ',VL0,VR0
       WRITE(*,*) 'PRESSURE             =  ',PL0,PR0
c      WRITE(*,*) 'TEMPERATURE          =  ',TL0,TR0
       WRITE(*,*) 'ENERGY               =  ',EL0,ER0
      ENDIF
      WRITE(*,*) 'MACH NO.             =  ',XMINF
      WRITE(*,*) 'ANGLE OF ATTACK      =  ',ALPHA
      WRITE(*,*) 'DENSITY              =  ',R0
      WRITE(*,*) 'VELOCITY U           =  ',U0
      WRITE(*,*) 'VELOCITY V           =  ',V0
      WRITE(*,*) 'PRESSURE             =  ',P0
      WRITE(*,*) 'TEMPERATURE          =  ',T0
      WRITE(*,*) 'ENERGY               =  ',E0
      ELSE
      WRITE(*,*) 'MACH NO.             =  ',XMINF
      WRITE(*,*) 'ANGLE OF ATTACK      =  ',ALPHA
      WRITE(*,*) 'DENSITY              =  ',R0
      WRITE(*,*) 'VELOCITY U           =  ',U0
      WRITE(*,*) 'VELOCITY V           =  ',V0
      WRITE(*,*) 'PRESSURE             =  ',P0
      WRITE(*,*) 'TEMPERATURE          =  ',T0
      WRITE(*,*) 'ENERGY               =  ',E0
      ENDIF
C
C
      NDIM(1)  = NI
      NDIM(2)  = NJ
C
      NBC      = 0
      NCOUNT   = 0
      NEXTPNT  = 1
C
 1000 CONTINUE
C..
C     IBTYPE = 0     ZERO GRADIENT
C            = 1     WALL
C            = 2     SYMMETRY I WITH Y=0
C            = 3     SYMMETRY J WITH Y=0
C            = 4     WAKE (FOR C GRID WAKE MATCH)
C            = 5     FARFIELD
C
      NCOUNT = NCOUNT + 1
C
C     READ(INAME,NML=BCIN)               ! work in SGI & CONVEX
C
C     READ(INAME,BCIN)                   ! work in CRAY
C
      READ(12,*) IBTYPE,KAS(1),KAE(1),KAS(2),KAE(2),
     1                  KBS(1),KBE(1),KBS(2),KBE(2)
C..
      IF( IBTYPE .LE. -1 .OR. BCNAME .EQ. 'QUIT' ) GO TO 2000
c     WRITE(*,1300) NCOUNT,BCNAME
      WRITE(*,1300) NCOUNT,IBTYPE
      WRITE(*,1400) IBTYPE,KAS,KAE
C     IF( IBTYPE.EQ.1 ) WRITE(6,1500) TWALL
      IF( IBTYPE.EQ.4 ) WRITE(6,1600) KBS,KBE
C
      IF( IBTYPE .GT. 7 ) THEN
         WRITE(6,*) ' ** ERROR ** IBTYPE AT NBC = ' , NBC
         STOP
      END IF
C
      NBC  = NBC + 1
C
      IDA   = 0
      DO 100 ID  = 1 , 2
         IF( KAS(ID) .EQ. KAE(ID) ) THEN
            IDA  = ID
            ISA   = 2
            IF( KAS(ID) .EQ. 1 ) ISA = 1
         END IF
  100 CONTINUE
C
      IF(IDA.EQ.0) THEN
        WRITE(*,*)' BC SET ERROR  AT NBC = ',NBC
        STOP
      ENDIF
C
      ICA = -1
      IF(KAE(CY(IDA)) .GT. KAS(CY(IDA)) ) ICA = 1
C
      IDBC(1,NBC) = IBTYPE
      IDBC(2,NBC) = IDA
      IDBC(3,NBC) = ISA
      IDBC(4,NBC) = ICA
      IDBC(5,NBC) = KAS(CY(IDA))
      IDBC(6,NBC) = KAE(CY(IDA))
C------------------------------------------------------------------------
C     ADDITION WALL BC PARAMETER
C
      IF(IBTYPE .EQ. 1) THEN
c      READ(12,*)RWALL,TWALL,UWALL,VWALL
c       FDBC( 1,NBC) = RWALL
c       FDBC( 2,NBC) = TWALL
c       FDBC( 3,NBC) = UWALL
c       FDBC( 4,NBC) = VWALL
      ENDIF
C
      IF( ISA .EQ. 1 ) THEN
         IIA   = -1
         IIB   =  0
         IIC   =  1
         IID   =  2
         IIE   =  3
      ELSE
         IIA   =  NDIM(IDA) + 2
         IIB   =  NDIM(IDA) + 1
         IIC   =  NDIM(IDA)
         IID   =  NDIM(IDA) - 1
         IIE   =  NDIM(IDA) - 2
      END IF
C
      I2S      = KAS(CY(IDA))
      I2E      = KAE(CY(IDA))
      LENG     = (I2E - I2S)*ICA + 1
      IP       = NEXTPNT - 1
C$DOIT SCALAR
      WRITE(9,*)'NBC = ',NBC
      WRITE(9,*)'IBTYPE = ',IBTYPE
      WRITE(9,*)'IDA = ',IDA
      WRITE(9,*)'ISA = ',ISA
      WRITE(9,*)'ICA = ',ICA
      WRITE(9,*)'I2S = ',I2S
      WRITE(9,*)'I2E = ',I2E
      WRITE(9,*)'LENG = ',LENG
      WRITE(9,*)'IP = ',IP
      DO 200 I2 = I2S , I2E ,ICA
         II    = I2 - I2S + 1
         IAI   = IIA*DD(1,IDA) + I2*DD(1,CY(IDA))
         IAJ   = IIA*DD(2,IDA) + I2*DD(2,CY(IDA))
         IBI   = IIB*DD(1,IDA) + I2*DD(1,CY(IDA))
         IBJ   = IIB*DD(2,IDA) + I2*DD(2,CY(IDA))
         ICI   = IIC*DD(1,IDA) + I2*DD(1,CY(IDA))
         ICJ   = IIC*DD(2,IDA) + I2*DD(2,CY(IDA))
         IDI   = IID*DD(1,IDA) + I2*DD(1,CY(IDA))
         IDJ   = IID*DD(2,IDA) + I2*DD(2,CY(IDA))
         IEI   = IIE*DD(1,IDA) + I2*DD(1,CY(IDA))
         IEJ   = IIE*DD(2,IDA) + I2*DD(2,CY(IDA))
         IJP(IP + LENG*0 + II ) = IAI
         IJP(IP + LENG*1 + II ) = IAJ
         IJP(IP + LENG*2 + II ) = IBI
         IJP(IP + LENG*3 + II ) = IBJ
         IJP(IP + LENG*4 + II ) = ICI
         IJP(IP + LENG*5 + II ) = ICJ
         IJP(IP + LENG*6 + II ) = IDI
         IJP(IP + LENG*7 + II ) = IDJ
         IJP(IP + LENG*8 + II ) = IEI
         IJP(IP + LENG*9 + II ) = IEJ
         WRITE(9,71)IAI,IAJ,IBI,IBJ,ICI,ICJ,IDI,IDJ,IEI,IEJ
   71    FORMAT(1X,10I5)
  200 CONTINUE
C
      IDBC(7,NBC) = IP
      NEXTPNT      = NEXTPNT + LENG*10
C
C     FOR WAKE BC
C
      IF(IBTYPE .EQ. 4 .OR. IBTYPE .EQ. 7) THEN
C
      NBC4  = NBC4 + 1
C
      NPABC = ABS(KAE(CY(IDA)) - KAS(CY(IDA))) + 1
C
      IDB  = 0
      DO 110 ID = 1 , 2
         IF( KBS(ID) .EQ. KBE(ID) ) THEN
            IDB = ID
            ISB = 2
            IF( KBS(ID) .EQ. 1 ) ISB = 1
         END IF
  110 CONTINUE
C
      IF(IDB.EQ.0) THEN
        WRITE(*,*)' BC SET ERROR  AT NBC = ',NBC
        STOP
      ENDIF
C
      ICB = -1
      IF(KBE(CY(IDB)) .GT. KBS(CY(IDB)) ) ICB = 1
C
      IDBC(12,NBC) = IDB
      IDBC(13,NBC) = ISB
      IDBC(14,NBC) = ICB
      IDBC(15,NBC) = KBS(CY(IDB))
      IDBC(16,NBC) = KBE(CY(IDB))
C
      IF( ISB .EQ. 1 ) THEN
         IIA   = -1
         IIB   =  0
         IIC   =  1
         IID   =  2
         IIE   =  3
      ELSE
         IIA   =  NDIM(IDB) + 2
         IIB   =  NDIM(IDB) + 1
         IIC   =  NDIM(IDB)
         IID   =  NDIM(IDB) - 1
         IIE   =  NDIM(IDB) - 2
      END IF
C
      I3S      = KBS(CY(IDB))
      I3E      = KBE(CY(IDB))
      LENB     = (I3E - I3S)*ICB + 1
      IPB      = NEXTPNT - 1
      WRITE(9,*)'NBC = ',NBC
      WRITE(9,*)'IBTYPE = ',IBTYPE
      WRITE(9,*)'IDA = ',IDA,'   IDB = ',IDB
      WRITE(9,*)'ISA = ',ISA,'   ISB = ',ISB
      WRITE(9,*)'ICA = ',ICA,'   ICB = ',ICB
      WRITE(9,*)'I2S = ',I2S,'   I3S = ',I3S
      WRITE(9,*)'I2E = ',I2E,'   I3E = ',I3E
      WRITE(9,*)'LENG = ',LENG,'   LENB = ',LENB
      WRITE(9,*)'IP = ',IP,'   IPB = ',IPB
C$DOIT SCALAR
      DO 210 I2 = I3S , I3E ,ICB
         II    = ICB*(I2 - I3S) + 1
         IAI   = IIA*DD(1,IDB) + I2*DD(1,CY(IDB))
         IAJ   = IIA*DD(2,IDB) + I2*DD(2,CY(IDB))
         IBI   = IIB*DD(1,IDB) + I2*DD(1,CY(IDB))
         IBJ   = IIB*DD(2,IDB) + I2*DD(2,CY(IDB))
         ICI   = IIC*DD(1,IDB) + I2*DD(1,CY(IDB))
         ICJ   = IIC*DD(2,IDB) + I2*DD(2,CY(IDB))
         IDI   = IID*DD(1,IDB) + I2*DD(1,CY(IDB))
         IDJ   = IID*DD(2,IDB) + I2*DD(2,CY(IDB))
         IEI   = IIE*DD(1,IDB) + I2*DD(1,CY(IDB))
         IEJ   = IIE*DD(2,IDB) + I2*DD(2,CY(IDB))
         IJP(IPB+ LENB*0 + II ) = IAI
         IJP(IPB+ LENB*1 + II ) = IAJ
         IJP(IPB+ LENB*2 + II ) = IBI
         IJP(IPB+ LENB*3 + II ) = IBJ
         IJP(IPB+ LENB*4 + II ) = ICI
         IJP(IPB+ LENB*5 + II ) = ICJ
         IJP(IPB+ LENB*6 + II ) = IDI
         IJP(IPB+ LENB*7 + II ) = IDJ
         IJP(IPB+ LENB*8 + II ) = IEI
         IJP(IPB+ LENB*9 + II ) = IEJ
         WRITE(9,71)IAI,IAJ,IBI,IBJ,ICI,ICJ,IDI,IDJ,IEI,IEJ
  210 CONTINUE
C
      IDBC(8,NBC) = IPB
      NEXTPNT = NEXTPNT + LENG*10
C
      NPBBC = ABS(KBE(CY(IDA)) - KBS(CY(IDA))) + 1
C
      IF(NPBBC .NE. NPABC ) THEN
        WRITE(*,*)' BC MATCH ERROR   NBC = ',NBC,NPBBC,NPABC,IDA,IDB
        STOP
      ENDIF
C
      ENDIF
C
      GO TO 1000
C
 2000 CONTINUE
C
      MXB = NBC
      MXBP = NEXTPNT
      MAXB = 20
C
      IF(MXB .GT. MAXB)THEN
        WRITE(*,*)' ERROR : MXB > MAXB (SETUP) INSUFF. BC NO.'
        WRITE(*,*)'         MXB = ',MXB,'     MAXB = ',MAXB
        STOP
      ENDIF
C
      IF(MXBP .GT. MAXBP)THEN
        WRITE(*,*)' ERROR : MXBP > MAXBP (SETUP) INSUFF. BC POINTS'
        WRITE(*,*)'         MXBP = ',MXBP,'     MAXBP = ',MAXBP
        STOP
      ENDIF
C
 1300 FORMAT( /2X,'SURFACE NO. ',I4,2X,I4,2X)
 1400 FORMAT(  5X,'BCTYPE = ',I5,
     1        /5X,'IBS    = ',I5,2X,'    JBS=',I5
     1        /5X,'IBE    = ',I5,2X,'    JBE=',I5)
 1500 FORMAT(  5X,'WALL TEMP. =',F5.2)
 1600 FORMAT(  5X,'IPS    = ',I5,2X,'    JPS=',I5
     1        /5X,'IPE    = ',I5,2X,'    JPE=',I5)
      RETURN
      END
C
      SUBROUTINE SOLU2(R,E,U,V,P,RNP,RUNP,RVNP,RENP,I,J)
C
      COMMON /FINP3/ ALPHA,UINF
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /NEWTON/ERRN
C
C      DIMENSION    R(-1:NI2,-1:NJ2), E(-1:NI2,-1:NJ2)
C      DIMENSION    U(-1:NI2,-1:NJ2), V(-1:NI2,-1:NJ2)
C      DIMENSION    P(-1:NI2,-1:NJ2)
C      DIMENSION  RNP(-1:NI2,-1:NJ2),RUNP(-1:NI2,-1:NJ2)
C      DIMENSION RVNP(-1:NI2,-1:NJ2),RENP(-1:NI2,-1:NJ2)
C
       XR       = RNP
       XM       = RUNP
       XN       = RVNP
       XE       = RENP
C
       IF( ABS(XM/XE) .LT. 1.E-8) THEN
         U = 0.
         IF( ABS(XN/XE) .LT. 1.E-8) THEN
          V = 0.
          R = XR
          E = XE
          XGAMAU = 1.
         ELSE
          VA     = XN/XE
          VB     = 0.
          XUV    = VA*VA
          IF( XUV .GT. 1. ) THEN
C            WRITE(6,*)'THE INITIAL POINT :'
C            WRITE(6,*) 'VA = ',VA
C            WRITE(6,*) 'XR,XM,XN,XE :'
C            WRITE(6,*)  XR,XM,XN,XE
            VA = V
           ENDIF
          RLGA   = 1./SQRT(ABS(1.-VA*VA))
          DA     = XR/RLGA
          EA     = XE-XN*VA
          PA     = GM1*(EA-DA)
          PB     = GM1*(XE-XR)
          FA     = (XE+PA)/(EA+PA)-RLGA*RLGA
          FB     = (XE+PB)/(XE+PB)-1.
          DO 10 K = 1,300
           VNEW = (VA+VB)/2.
           XUV    = VNEW*VNEW
C          IF( XUV .GT. 1. ) THEN
C            WRITE(6,*)'THE DO_LOOP:'
C            WRITE(6,*) 'VA = ',VNEW
C            WRITE(6,*) 'XR,XM,XN,XE :'
C            WRITE(6,*)  XR,XM,XN,XE
C            STOP
C           ENDIF
           RLGAN= 1./SQRT(ABS(1.-VNEW*VNEW))
           DNEW = XR/RLGAN
           ENEW = XE-XN*VNEW
           PNEW = GM1*(ENEW-DNEW)
           FNEW = (XE+PNEW)/(ENEW+PNEW)-RLGAN*RLGAN
           IF ( ABS(FNEW) .LE. ERRN ) GOTO 10
           IF ( FA*FNEW .GT. 0. ) THEN
            VA  = VNEW
            FA  = FNEW
           ELSE
            VB  = VNEW
c           FA  = FNEW
           ENDIF
10        CONTINUE
          V= VNEW
           XUV    = VNEW*VNEW
C          IF( XUV .GT. 1. ) THEN
C            WRITE(6,*)'THE OUTPUT:'
C            WRITE(6,*) 'VA = ',VNEW
C            WRITE(6,*) 'XR,XM,XN,XE :'
C            WRITE(6,*)  XR,XM,XN,XE
C            STOP
C           ENDIF
          R = ABS(XR*SQRT(ABS(1.-V*V)))
C          IF( R .LE. -0.1 )THEN
C           WRITE(6,*)'DENSITY ERROR AT I,J=',I,J,R
C           STOP
C          ENDIF
          E = ABS(XE-XN*V)
         ENDIF
        ELSE
         VA   = XM/XE
         VB   = 0.
         VY   = XN/XM*VA
         UV2  = VA*VA + VY*VY
         IF( UV2 .GE. 1.0 )THEN
C          WRITE(6,*)'U,V =',VA,VY
C          WRITE(6,*)'Error AT FUNCTION VEL2D'
          VA = U
         ENDIF
         VY = XN/XM*VA
         RLGA = 1./SQRT(ABS(1.-VA*VA-VY*VY))
         DA   = XR/RLGA
         EA   = XE-XM*VA-XN*VY
         PA   = GM1*(EA-DA)
         PB   = GM1*(XE-XR)
         FA   = (XE+PA)/(EA+PA)-RLGA*RLGA
         FB   = (XE+PB)/(XE+PB)-1.
         DO 20 K = 1,300
          VNEW = (VA+VB)/2.
          VY   = XN/XM*VNEW
           XUV    = VNEW*VNEW+VY*VY
C          IF( XUV .GT. 1. ) THEN
C            WRITE(6,*)'THE DO_LOOPIJ:I,J=',I,J
C            WRITE(6,*) 'VNEW,VY = ',VNEW,VY
C            WRITE(6,*) 'XR,XM,XN,XE :'
C            WRITE(6,*)  XR,XM,XN,XE
C            STOP
C           ENDIF
          RLGAN= 1./SQRT(ABS(1.-VNEW*VNEW-VY*VY))
          DNEW = XR/RLGAN
          ENEW = XE-XM*VNEW-XN*VY
          PNEW = GM1*(ENEW-DNEW)
          FNEW = (XE+PNEW)/(ENEW+PNEW)-RLGAN*RLGAN
          IF (ABS(FNEW) .LE. ERRN ) GOTO 20
          IF (FA*FNEW .GT. 0.) THEN
           VA  = VNEW
           FA  = FNEW
          ELSE
           VB  = VNEW
c          FA  = FNEW
          ENDIF
20       CONTINUE
         U = VNEW
         V = XN/XM*U
           XUV    = VNEW*VNEW+V*V
C          IF( XUV .GT. 1. ) THEN
C            WRITE(6,*)'THE OUTPUTIJ:I,J=',I,J
C            WRITE(6,*) 'U,V = ',VNEW,V
C            WRITE(6,*) 'XR,XM,XN,XE :'
C            WRITE(6,*)  XR,XM,XN,XE
C            STOP
C           ENDIF
         XGAMAU = SQRT(ABS(1.-U*U-V*V))
         R = ABS(XR*XGAMAU)
C         IF( R .LE. -0.1 )THEN
C          WRITE(6,*)'DENSITY ERROR AT I,J=',I,J,R
C          STOP
C         ENDIF
         E = ABS(XE-XM*U-XN*V)
        ENDIF
      P     = GM1*ABS( E-R )
C
      RETURN
      END
      SUBROUTINE SOLUV(R,E,U,V,P,RNP,RUNP,RVNP,RENP)
C
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /NEWTON/ERRN
C
      DIMENSION    R(-1:NI2,-1:NJ2), E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2), V(-1:NI2,-1:NJ2)
      DIMENSION    P(-1:NI2,-1:NJ2)
      DIMENSION  RNP(-1:NI2,-1:NJ2),RUNP(-1:NI2,-1:NJ2)
      DIMENSION RVNP(-1:NI2,-1:NJ2),RENP(-1:NI2,-1:NJ2)
C
      DO 200 I =2,NI2
      DO 200 J =2,NJ2
       XR       = ABS(RNP(I,J))
       XM       = RUNP(I,J)
       XN       = RVNP(I,J)
       XE       = ABS(RENP(I,J))
       RNP(I,J) = ABS(RNP(I,J))
       RENP(I,J) = ABS(RENP(I,J))
C
       IF( ABS(XM) .LT. 1.E-10) THEN
         U(I,J) = 0.
         IF( ABS(XN) .LT. 1.E-10) THEN
          V(I,J) = 0.
          R(I,J) = XR
          E(I,J) = XE
          XGAMAU = 1.
         ELSE
CCC *** use C.D.MUNZ's formula to obtain the initial V0 ***
C
          DELTA = 1.E-6
          TOL   = 1.E-6
          AXN   = ABS(XN)
          VA    = MIN(1.,AXN/XE+DELTA)
          VD    = SQRT(ABS((GAMMA*XE)**2-4.*GM1*AXN**2))
          VB    = (GAMMA*XE-VD)/(2.*AXN*GM1)
          Z     = 0.
          IF (ABS(VB) .GT. DELTA) Z=ABS(1.-XR/XE)*(VB-VA)/2.
          V0    =(VA+VB)/2.+Z
          IF( XN .LT. 0. ) V0 = -V0
C
C          VA    = MAX(-1.,XN/XE-DELTA)
C          VD    = SQRT(ABS((GAMMA*XE)**2-4.*GM1*XN**2))
C          VB    = ABS(GAMMA*XE-VD)/(2.*XN*GM1)
C          Z     = 0.
C          IF (ABS(VB) .GT. DELTA) Z=ABS(1.-XR/XE)*(VB-VA)/2.
C          V0    = (VA+VB)/2.+Z
C          ENDIF 
C
          ITER  = 0
 10       CONTINUE
          ITER  = ITER + 1
          U0    = 0.
          X1    = 1.-V0**2
          E0    = XE-XN*V0
          D0    = XR*SQRT(ABS(X1))
          P0    = GM1*(E0-D0)
CCC ****** ((((  Newton's Method )))) *******
          FUNH   = GAMMA*V0*E0-XN*X1
          FUNK   = X1*(V0*GM1*XR)**2
          DH     = GAMMA*XE-2.*GM1*XN*V0
          DK     = 2.*V0*(1.-2.*V0**2)*(GM1*XR)**2
          FUNF   = FUNH**2 - FUNK
          DF     = 2.*FUNH*DH-DK
          DELTAV = -FUNF/DF
          VNEW   = V0 + DELTAV
          XERR   = ABS((VNEW-V0)/V0)
          V0     = VNEW
          IF (XERR .GT. TOL .AND. ITER .LT. 100) GOTO 10
          IF( ITER .GT. 50 ) THEN
c            WRITE(6,*) 'VELOCITY AT I,J=',I,J,V0
          ENDIF
          V(I,J) = V0
          X1     = 1.-V(I,J)**2
          IF( X1 .LT. 0.) THEN
c           WRITE(6,*)'VELOCITY ERROR AT I,J=',I,J,V(I,J)
           STOP
          ENDIF
C          AX1    = ABS(X1) 
C          R(I,J) = XR*SQRT(AX1)
          R(I,J) = XR*SQRT(X1)
          IF( R(I,J) .LE. -1.0 )THEN
           WRITE(6,*)'DENSITY ERROR AT I,J=',I,J,R(I,J)
           STOP
          ENDIF
          E(I,J) = XE-XN*V(I,J)
         ENDIF
        ELSE
CCC *** use C.D.MUNZ's formula to obtain the initial U0 ***
C
         AXM   = ABS(XM)
         XNM    = XN/AXM
         XA    = AXM+XNM*XN
         DELTA = 1.E-6
         TOL   = 1.E-6
         UA    = MIN(1.,AXM/XE+DELTA)
         UD    = SQRT(ABS((GAMMA*XE)**2-4.*GM1*XA*AXM))
         UB    = (GAMMA*XE-UD)/(2.*XA*GM1)
         Z     = 0.
         IF (ABS(UB) .GT. DELTA) Z=ABS(1-XR/XE)*(UB-UA)/2.
         U0    =(UA+UB)/2.+Z
         IF( XM .LT. 0.) U0 = -U0
C         ELSE
C         UA    = MAX(-1.,XM/XE-DELTA)
C         UD    = SQRT(ABS((GAMMA*XE)**2-4.*GM1*XA*XM))
C         UB    = ABS(GAMMA*XE-UD)/(2.*XA*GM1)
C         Z     = 0.
C         IF (ABS(UB) .GT. DELTA) Z=ABS(1-XR/XE)*(UB-UA)/2.
C         U0    =(UA+UB)/2.+Z
C         ENDIF
C
         XNM    = XN/XM
C
         ITER  = 0
 15      CONTINUE
         ITER  = ITER + 1
         V0    = XNM*U0
         X1    = 1.-U0**2-V0**2
         E0    = XE-XM*U0-XN*V0
         D0    = XR*SQRT(ABS(X1))
         P0    = GM1*(E0-D0)
CCC ****** ((((  Newton's Method )))) *******
         FUNH   = GAMMA*U0*E0-XM*X1
         FUNK   = X1*(U0*GM1*XR)**2
         X2     = 1.+XNM**2
         DH     = GAMMA*(E0-U0*XA)+2.*XM*U0*X2
         DK     = 2.*U0*(1-2.*X2*U0**2)*(GM1*XR)**2
         FUNF   = FUNH**2 - FUNK
         DF     = 2.*FUNH*DH-DK
         DELTAV = -FUNF/DF
         UNEW   = U0 + DELTAV
         XERR   = ABS((UNEW-U0)/U0)
         U0     = UNEW
         IF (XERR .GT. TOL .AND. ITER .LT. 100) GOTO 15
         U(I,J) = U0
         V(I,J) = XNM*U(I,J)
          IF( ITER .GT. 50 ) THEN
c            WRITE(6,*) 'VELOCITY AT I,J=',I,J,U0,V(I,J)
          ENDIF
         X1     = 1.-U(I,J)*U(I,J)-V(I,J)*V(I,J)
         IF( X1 .LT. 0. ) THEN
c          WRITE(6,*) 'X1 =',X1
c          WRITE(6,*)'VELOCITY ERROR AT I,J=',I,J,U(I,J),V(I,J)
          STOP
         ENDIF
C         AX1    = ABS(1.-U(I,J)*U(I,J)-V(I,J)*V(I,J))
C         XGAMAU = SQRT(AX1)
         X1    = 1.-U(I,J)*U(I,J)-V(I,J)*V(I,J)
         XGAMAU = SQRT(X1)
         R(I,J) = XR*XGAMAU
         IF( R(I,J) .LE. -1.0 )THEN
          WRITE(6,*)'DENSITY ERROR AT I,J=',I,J,R(I,J)
          STOP
         ENDIF
         E(I,J) = XE-XM*U(I,J)-XN*V(I,J)
        ENDIF
      P(I,J)     = GM1*( E(I,J)-R(I,J) )
  200 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE SOLV(IJP,MLU,X,Y,XIX,XIY,ETAX,ETAY,QJ,R,E,T,U,V,P
     1 ,XM,PSI,US,VS,RS,RUS,RVS,RES,DRS,RNP,RUNP,RVNP,RENP,RNPO)
C
      COMMON /TAPE / INAME,IGRID,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXG,MAXL,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /PNTI / LKDI,LKDJ,LLKD,LKDP,MDA,MDB,NDIG
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
      COMMON /FINP3/ ALPHA,UINF
      COMMON /IKTT/ KTT
      COMMON /CNST2/ R0,U0,V0,P0,T0,E0
      COMMON /CNST4/ DR0,UR0,VR0,PR0,TR0,ER0,
     1               DL0,UL0,VL0,PL0,TL0,EL0
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /BC01 / IDBC(30,20),FDBC(10,20)
      COMMON /AXIS / IAXIS
      COMMON /ENTR / EPS,EPSS,EPS2
      COMMON /RES01/ RDMX,IDMX,JDMX,RTMX,ITMX,JTMX,RDL2,RTL2
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
      COMMON /PSII/ ICYL,PSI0,AMPLI,NCURVE,URATIO,YY0
      COMMON /T6/ AXISYM
      COMMON /U1/ LTP
      LOGICAL AXISYM
C
      DIMENSION IJP(MXBP),MLU(MAXLI)
      DIMENSION    X(-1:NI2,-1:NJ2),     Y(-1:NI2,-1:NJ2)
      DIMENSION  XIX(-1:NI2,-1:NJ2),   XIY(-1:NI2,-1:NJ2)
      DIMENSION ETAX(-1:NI2,-1:NJ2),  ETAY(-1:NI2,-1:NJ2)
      DIMENSION   QJ(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2),     E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2),     V(-1:NI2,-1:NJ2)
      DIMENSION    P(-1:NI2,-1:NJ2),     T(-1:NI2,-1:NJ2)
      DIMENSION   XM(-1:NI2,-1:NJ2),   PSI(-1:NI2,-1:NJ2)
c     DIMENSION DRUS(-1:NI2,-1:NJ2,5),DRVS(-1:NI2,-1:NJ2,5)
c     DIMENSION DRES(-1:NI2,-1:NJ2,5)
      DIMENSION  DRS(-1:NI2,-1:NJ2,5)
      DIMENSION  RUS(-1:NI2,-1:NJ2,5), RVS(-1:NI2,-1:NJ2,5)
      DIMENSION   RS(-1:NI2,-1:NJ2,5), RES(-1:NI2,-1:NJ2,5)
      DIMENSION   US(-1:NI2,-1:NJ2,5),  VS(-1:NI2,-1:NJ2,5)
      DIMENSION  RNP(-1:NI2,-1:NJ2,3),RUNP(-1:NI2,-1:NJ2,3)
      DIMENSION RVNP(-1:NI2,-1:NJ2,3),RENP(-1:NI2,-1:NJ2,3)
      DIMENSION RNPO(-1:NI2,-1:NJ2,4)
C
      CALL GRIDIN(X,Y,IJP)
C
      CALL METJAC(X,Y,XIX,XIY,ETAX,ETAY,QJ,IJP)
C
      IF(IRSTRT .NE. 0) THEN
*      CALL ASSIGN('assign -F f77 -N ieee u:16')
      CALL RSTQ(R,U,V,P,E,PSI,ITER0,TIME1)
*      CALL ASSIGN('assign -R u:16')
      ITER0 = 0
      ELSE
      CALL INITQ(X,Y,R,T,U,V,P,E,PSI)
      ITER0 = 0
      TIME1 = 0.
      ENDIF
C
      CALL BEAM5(R,E,RS,RUS,RVS,RES,U,V,P,US,VS,
     1           XIX,XIY,ETAX,ETAY,QJ)
C
      CALL DTIME(US,VS)
C
      ITER     = ITER0
      TIME     = TIME1
       KT      =  1
      IF( LSTD .EQ. 0 .AND. IRSTRT .EQ. 0 )THEN
        KTT = 30
      ENDIF
      IFILE = KTT
C
      CALL QPLT3D(IFILE,R,U,V,P,E,XM,X,XIX,XIY,ETAX,ETAY,PSI,DRS,Y)
C
 1000 CONTINUE
c      T1 = SECOND()
       call cpu_time(t1)
C
      CALL BEAM5(R,E,RS,RUS,RVS,RES,U,V,P,US,VS,
     1           XIX,XIY,ETAX,ETAY,QJ)
C
      CALL DTIME(US,VS)
C
       LTP      = 0
      CFL      = (CFL2 - CFL1)/NUP*ITER + CFL1
      IF(CFL .GE. CFL2) CFL = CFL2
      DTCFL    = CFL * DT
      IF( LSTD .EQ. 0 )THEN
       IF(TIME.GE.TPP(KT))THEN
         DTCFL = TPP(KT) - (TIME - DTCFL)
c        DTT   = DT
c        DT    = DTCFL/CFL
c        DT    = DTT
         LTP   = 1
         TIME  = TPP(KT)
         KT    = KT + 1
       ENDIF
      ENDIF
c      T2 = SECOND()
       call cpu_time(t2)
C      CALL BEAM5(R,E,RS,RUS,RVS,RES,U,V,P,US,VS,
C     1           XIX,XIY,ETAX,ETAY,QJ)
C
      IF( MTHR .EQ. 9 .OR. MTHR .EQ. 11 )THEN
       CALL RUNGE2(IJP,MLU,X,Y,R,E,U,V,T,P,QJ,XIX,XIY,ETAX,ETAY,
     1  RS,RUS,RVS,RES,US,VS,DRS,RNP,RUNP,RVNP,RENP,RNPO)
      ELSEIF( MTHR .EQ. 10 .OR. MTHR .EQ. 12 )THEN
       CALL RUNGE3(IJP,MLU,X,Y,R,E,U,V,T,P,QJ,XIX,XIY,ETAX,ETAY,
     1  RS,RUS,RVS,RES,US,VS,DRS,RNP,RUNP,RVNP,RENP,RNPO)
      ELSE
       IF( MOD(ITER,2) .EQ. 1 )THEN
C  ---
C      SOLVE THE SOURCE TERM  FOR AXISYMMETRIC FLOW
C  ---
      IF( AXISYM ) THEN
      CALL SOURCE( X,Y,R,E,U,V,T,P,QJ,
     1            RS,RUS,RVS,RES,US,VS,DRS,RNPO)
C
      CALL BC(IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P,X)
      CALL BEAM5(R,E,RS,RUS,RVS,RES,U,V,P,US,VS,
     1           XIX,XIY,ETAX,ETAY,QJ)
C
      ENDIF
        CALL STEP(MLU(LKDI),MLU(LKDJ),MLU(LLKD),MLU(LKDP),
     1   IJP,R,E,U,V,T,P,QJ,RS,RUS,RVS,RES,US,VS,DRS,RNPO)
       ELSE
        CALL STEP1(MLU(LKDI),MLU(LKDJ),MLU(LLKD),MLU(LKDP),
     1   IJP,R,E,U,V,T,P,QJ,RS,RUS,RVS,RES,US,VS,DRS,RNPO)
C  ---
C      SOLVE THE SOURCE TERM  FOR AXISYMMETRIC FLOW
C  ---
       IF( AXISYM ) THEN
       CALL BC(IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P,X)
       CALL SOURCE( X,Y,R,E,U,V,T,P,QJ,
     1            RS,RUS,RVS,RES,US,VS,DRS,RNPO)
C
       ENDIF
       ENDIF
      ENDIF
c       T3 = SECOND()
        call cpu_time(t3)
c       write(6,*)'Time of sub. step=',t3-t2
C
      CALL BC(IJP,XIX,XIY,ETAX,ETAY,R,E,U,V,T,P,X)
C
      IF( ICYL .EQ. 0 )THEN
       CALL ITFACE(IJP,R,U,V,QJ,PSI,XIX,XIY,ETAX,ETAY,DRS)
      ENDIF
C
      XMMX      = -1.
      NSUP      = 0
c     DO 110 IJ = -1,MIJ2
c     I         = MOD(IJ+1,MI)-1
c     J         = (IJ+1)/MI -1
c     XM(IJ,-1) = SQRT(U(IJ,-1)*U(IJ,-1)+V(IJ,-1)*V(IJ,-1))*UTX
c     XM(IJ,-1) = SQRT( (U(IJ,-1)*U(IJ,-1)+V(IJ,-1)*V(IJ,-1))
c    +                 /T(IJ,-1) )*UTX
c     IF(I.GT.0.AND.I.LT.NI1.AND.J.GT.0.AND.J.LT.NJ1)THEN
c     IF(XM(IJ,-1) .GT. 1.) NSUP = NSUP + 1
c     IF(XM(IJ,-1) .GT. XMMX) THEN
c        XMMX   = XM(IJ,-1)
c        IMMX   = I
c        JMMX   = J
c     ENDIF
c     ENDIF
c 110 CONTINUE
C
       IF((MOD(ITER,IPRINT).EQ.0).OR.(ITER.EQ.NSTEP)) THEN
        IF( LSTD .EQ. 1 )THEN
         WRITE(*,11)ITER,DTCFL,TIME,LOG10(RDL2),LOG10(RTL2),NSUP,
     1              RDMX,IDMX,JDMX,RTMX,ITMX,JTMX
   11    FORMAT(I4,1X,F6.4,1X,F6.3,2(1X,F6.3),1X,I5,
     1           2(1X,F5.3,1X,I3,1X,I3))
C
        ELSE
         WRITE(*,13)ITER,DTCFL,TIME,NSUP
   13    FORMAT(I5,1X,1PE13.5,1X,E13.5,I7)
C
        ENDIF
       ENDIF
C
      IF((MOD(ITER,IPRINT).EQ.0).OR.(ITER.EQ.NSTEP)
     1    .OR. LTP .EQ. 1) THEN
c     IF((MOD(ITER,IPRINT).EQ.0).OR.(ITER.EQ.NSTEP)) THEN
C
      IF( LSTD .EQ. 0 )THEN
       IF((MOD(ITER,IPRINT).EQ.0) .AND. LTP .EQ. 0)IFILE = KTT+KT
       IF(LTP .EQ. 1) IFILE = KTT + KT - 1
      WRITE(*,99)IFILE
   99 FORMAT(1X,'OUTPUT Q DATA TO FILE fort.',I2)
      ELSE
        IFILE = 4
      ENDIF
      CALL QPLT3D(IFILE,R,U,V,P,E,XM,X,XIX,XIY,ETAX,ETAY,PSI,DRS,Y)
C
c     IF(MTHL .EQ. 1)THEN
c      DO 120 IT = 1,ITER,5
c      WRITE(IRESD,12)IT-1,LOG10(RDL22(IT-1)),LOG10(RTL22(IT-1))
c 120  CONTINUE
c  12  FORMAT(1X,I10,2(2X,F9.5))
c      CLOSE(IRESD)
c     ENDIF
C
      ENDIF
C
      IF( LSTD .EQ. 0 )THEN
        IF(ITER.GE.NSTEP .OR. TIME .GE. TPP(MTIM) ) GO TO 1100
      ELSE
        IF(ITER.GE.NSTEP ) GO TO 1100
      ENDIF
      ITER = ITER + 1
      TIME = TIME + DTCFL
c      T8 = SECOND()
       call cpu_time(t8)
c     write(6,*)'Time of each step=',t8-t1
      GO TO 1000
 1100 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE SOURCE(X,Y,R,E,U,V,T,P,QJ,
     1      RS,RUS,RVS,RES,US,VS,DRS,RNPO)
C
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /PNTI / LKDI,LKDJ,LLKD,LKDP,MDA,MDB,NDIG
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /RES01/ RDMX,IDMX,JDMX,RTMX,ITMX,JTMX,RDL2,RTL2
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
C
      DIMENSION    X(-1:NI2,-1:NJ2),     Y(-1:NI2,-1:NJ2)
      DIMENSION   QJ(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2),     E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2),     V(-1:NI2,-1:NJ2)
      DIMENSION    T(-1:NI2,-1:NJ2),     P(-1:NI2,-1:NJ2)
      DIMENSION  DRS(-1:NI2,-1:NJ2,5)
      DIMENSION  RUS(-1:NI2,-1:NJ2,5), RVS(-1:NI2,-1:NJ2,5)
      DIMENSION   RS(-1:NI2,-1:NJ2,5), RES(-1:NI2,-1:NJ2,5)
      DIMENSION   US(-1:NI2,-1:NJ2,5),  VS(-1:NI2,-1:NJ2,5)
      DIMENSION RNPO(-1:NI2,-1:NJ2,4), S(-1:NI2,-1:NJ2,4)
      DO 10 I=-1,NI2
      DO 10 J=-1,NJ2
       IF( ABS(Y(I,J)) .LE. 1.E-6 ) GOTO 10
       XGAMMAU     = 1./SQRT(ABS(1.-U(I,J)**2-V(I,J)**2))
       YINV        = 1./ABS(Y(I,J))
       XGD         = XGAMMAU**2
       XEP         = E(I,J)+P(I,J)
       RNPO(I,J,1)  =  XGAMMAU*R(I,J)
       RNPO(I,J,2)  =  XGD*XEP*U(I,J)
       RNPO(I,J,3)  =  XGD*XEP*V(I,J)
       RNPO(I,J,4)  =  XGD*XEP-P(I,J)
C
       S(I,J,1)  = RNPO(I,J,1)*V(I,J)*YINV
       S(I,J,2)  = RNPO(I,J,2)*V(I,J)*YINV
       S(I,J,3)  = RNPO(I,J,3)*V(I,J)*YINV
       S(I,J,4)  =(RNPO(I,J,4)+P(I,J))*V(I,J)*YINV
 10   CONTINUE
C ---
C  UPDATE
C ---
      DO 20 I=-1,NI2
      DO 20 J=-1,NJ2
      DO 20 K=1,4
       IF( ABS(Y(I,J)) .LE. 1.E-6 ) GOTO 20
       RNPO(I,J,K) = RNPO(I,J,K) - DTCFL*S(I,J,K)
 20   CONTINUE
      CALL SOLUV(R,E,U,V,P,RNPO(-1,-1,1),RNPO(-1,-1,2),RNPO(-1,-1,3),
     1           RNPO(-1,-1,4))
      RETURN
      END
C
      SUBROUTINE STEP(KDI,KDJ,LKD,KDP,IJP,R,E,U,V,T,P,QJ,
     1      RS,RUS,RVS,RES,US,VS,DRS,RNPO)
C
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /PNTI / LKDI,LKDJ,LLKD,LKDP,MDA,MDB,NDIG
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /RES01/ RDMX,IDMX,JDMX,RTMX,ITMX,JTMX,RDL2,RTL2
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
C
      DIMENSION IJP(MXBP),KDI(MDA),KDJ(MDA),LKD(MDB),KDP(MDB)
      DIMENSION   QJ(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2),     E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2),     V(-1:NI2,-1:NJ2)
      DIMENSION    T(-1:NI2,-1:NJ2),     P(-1:NI2,-1:NJ2)
      DIMENSION  DRS(-1:NI2,-1:NJ2,5)
      DIMENSION  RUS(-1:NI2,-1:NJ2,5), RVS(-1:NI2,-1:NJ2,5)
      DIMENSION   RS(-1:NI2,-1:NJ2,5), RES(-1:NI2,-1:NJ2,5)
      DIMENSION   US(-1:NI2,-1:NJ2,5),  VS(-1:NI2,-1:NJ2,5)
      DIMENSION RNPO(-1:NI2,-1:NJ2,4)
C
       DO 10 IJ      = -1,MIJ2
        DRS(IJ,-1,1) = 0.
        DRS(IJ,-1,2) = 0.
        DRS(IJ,-1,3) = 0.
        DRS(IJ,-1,4) = 0.
        DRS(IJ,-1,5) = 0.
10     CONTINUE
      IF(MTHR .EQ. 1) THEN
       CALL BMUPW(RS,US,VS,DRS)
      ELSEIF(MTHR .EQ. 2)THEN
cmic$ do all shared(qj) autoscope
        DO 12 NE  = 1,5
        CALL FVITVD2(US(-1,-1,NE),QJ,RS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVJTVD2(VS(-1,-1,NE),QJ,RS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
12      CONTINUE
      ELSEIF(MTHR .EQ. 3)THEN
cmic$ do all shared(qj) autoscope
        DO 14 NE  = 1,5
        CALL FVIEN2(US(-1,-1,NE),QJ, RS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVJEN2(VS(-1,-1,NE),QJ, RS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
14      CONTINUE
      ENDIF
       DO 16 IJ    =  -1,MIJ2
        RS(IJ,-1,1)  = RS(IJ,-1,1) + DRS(IJ,-1,1)
        RS(IJ,-1,2)  = RS(IJ,-1,2) + DRS(IJ,-1,2)
        RS(IJ,-1,3)  = RS(IJ,-1,3) + DRS(IJ,-1,3)
        RS(IJ,-1,4)  = RS(IJ,-1,4) + DRS(IJ,-1,4)
        RS(IJ,-1,5)  = RS(IJ,-1,5) + DRS(IJ,-1,5)
16     CONTINUE
       DO 20 IJ      = -1,MIJ2
        DRS(IJ,-1,1) = 0.
        DRS(IJ,-1,2) = 0.
        DRS(IJ,-1,3) = 0.
        DRS(IJ,-1,4) = 0.
        DRS(IJ,-1,5) = 0.
20     CONTINUE
      IF(MTHR .EQ. 1) THEN
       CALL BMUPW(RUS,US,VS,DRS)
      ELSEIF(MTHR .EQ. 2) THEN
cmic$ do all shared(qj) autoscope
        DO 22 NE  = 1,5
        CALL FVITVD2(US(-1,-1,NE),QJ,RUS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVJTVD2(VS(-1,-1,NE),QJ,RUS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
22      CONTINUE
      ELSEIF(MTHR .EQ. 3) THEN
cmic$ do all shared(qj) autoscope
        DO 24 NE  = 1,5
        CALL FVIEN2(US(-1,-1,NE),QJ, RUS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVJEN2(VS(-1,-1,NE),QJ, RUS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
24      CONTINUE
      ENDIF
       DO 26 IJ    =  -1,MIJ2
        RUS(IJ,-1,1)  = RUS(IJ,-1,1) + DRS(IJ,-1,1)
        RUS(IJ,-1,2)  = RUS(IJ,-1,2) + DRS(IJ,-1,2)
        RUS(IJ,-1,3)  = RUS(IJ,-1,3) + DRS(IJ,-1,3)
        RUS(IJ,-1,4)  = RUS(IJ,-1,4) + DRS(IJ,-1,4)
        RUS(IJ,-1,5)  = RUS(IJ,-1,5) + DRS(IJ,-1,5)
26     CONTINUE
       DO 30 IJ      = -1,MIJ2
        DRS(IJ,-1,1) = 0.
        DRS(IJ,-1,2) = 0.
        DRS(IJ,-1,3) = 0.
        DRS(IJ,-1,4) = 0.
        DRS(IJ,-1,5) = 0.
30     CONTINUE
      IF(MTHR .EQ. 1) THEN
       CALL BMUPW(RVS,US,VS,DRS)
      ELSEIF(MTHR .EQ. 2) THEN
cmic$ do all shared(qj) autoscope
        DO 32 NE  = 1,5
        CALL FVITVD2(US(-1,-1,NE),QJ,RVS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVJTVD2(VS(-1,-1,NE),QJ,RVS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
32      CONTINUE
      ELSEIF(MTHR .EQ. 3) THEN
cmic$ do all shared(qj) autoscope
        DO 34 NE  = 1,5
        CALL FVIEN2(US(-1,-1,NE),QJ, RVS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVJEN2(VS(-1,-1,NE),QJ, RVS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
34      CONTINUE
      ENDIF
       DO 36 IJ    =  -1,MIJ2
        RVS(IJ,-1,1)  = RVS(IJ,-1,1) + DRS(IJ,-1,1)
        RVS(IJ,-1,2)  = RVS(IJ,-1,2) + DRS(IJ,-1,2)
        RVS(IJ,-1,3)  = RVS(IJ,-1,3) + DRS(IJ,-1,3)
        RVS(IJ,-1,4)  = RVS(IJ,-1,4) + DRS(IJ,-1,4)
        RVS(IJ,-1,5)  = RVS(IJ,-1,5) + DRS(IJ,-1,5)
36     CONTINUE
       DO 40 IJ      = -1,MIJ2
        DRS(IJ,-1,1) = 0.
        DRS(IJ,-1,2) = 0.
        DRS(IJ,-1,3) = 0.
        DRS(IJ,-1,4) = 0.
        DRS(IJ,-1,5) = 0.
40     CONTINUE
      IF(MTHR .EQ. 1) THEN
       CALL BMUPW(RES,US,VS,DRS)
      ELSEIF(MTHR .EQ. 2) THEN
cmic$ do all shared(qj) autoscope
        DO 42 NE  = 1,5
        CALL FVITVD2(US(-1,-1,NE),QJ,RES(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVJTVD2(VS(-1,-1,NE),QJ,RES(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
42      CONTINUE
      ELSEIF(MTHR .EQ. 3) THEN
cmic$ do all shared(qj) autoscope
        DO 44 NE  = 1,5
        CALL FVIEN2(US(-1,-1,NE),QJ, RES(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVJEN2(VS(-1,-1,NE),QJ, RES(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
44      CONTINUE
      ENDIF
       DO 46 IJ    =  -1,MIJ2
        RES(IJ,-1,1)  = RES(IJ,-1,1) + DRS(IJ,-1,1)
        RES(IJ,-1,2)  = RES(IJ,-1,2) + DRS(IJ,-1,2)
        RES(IJ,-1,3)  = RES(IJ,-1,3) + DRS(IJ,-1,3)
        RES(IJ,-1,4)  = RES(IJ,-1,4) + DRS(IJ,-1,4)
        RES(IJ,-1,5)  = RES(IJ,-1,5) + DRS(IJ,-1,5)
46     CONTINUE
C
      DO 100 IJ = -1,MIJ2
       RNPO(IJ,-1,1) = QJ(IJ,-1)*( RS(IJ,-1,1)+RS(IJ,-1,2)+
     1                RS(IJ,-1,3)+RS(IJ,-1,4)+RS(IJ,-1,5) )
       RNPO(IJ,-1,2) = QJ(IJ,-1)*( RUS(IJ,-1,1)+RUS(IJ,-1,2)+
     1                RUS(IJ,-1,3)+RUS(IJ,-1,4)+RUS(IJ,-1,5) )
       RNPO(IJ,-1,3) = QJ(IJ,-1)*( RVS(IJ,-1,1)+RVS(IJ,-1,2)+
     1                RVS(IJ,-1,3)+RVS(IJ,-1,4)+RVS(IJ,-1,5) )
       RNPO(IJ,-1,4) = QJ(IJ,-1)*( RES(IJ,-1,1)+RES(IJ,-1,2)+
     1                RES(IJ,-1,3)+RES(IJ,-1,4)+RES(IJ,-1,5) )
100   CONTINUE
C
      SUMD         = 0.
      SUMT         = 0.
      RDMX         = -1.
      RTMX         = -1.
      DO 200 I     = 2,NI-1
      DO 200 J     = 2,NJ-1
       RNP1         = RNPO(I,J,1)
       RUNP1        = RNPO(I,J,2)
       RVNP1        = RNPO(I,J,3)
       RENP1        = RNPO(I,J,4)
C
       XU           = 1./SQRT(ABS(1.-U(I,J)**2-V(I,J)**2))
       XU2          = XU*XU
       EP           = E(I,J) + P(I,J)
       DSQD         = (RNP1 - R(I,J)*XU)**2
       DSQT         = DSQD + (RUNP1 - XU2*EP*U(I,J))**2
     1                     + (RVNP1 - XU2*EP*V(I,J))**2
     1                     + (RENP1 - XU2*EP+P(I,J))**2
       SUMD         = SUMD + DSQD
       SUMT         = SUMT + DSQT
       IF(SQRT(DSQD) .GT. RDMX)THEN
         RDMX       = SQRT(DSQD)
         IDMX       = I
         JDMX       = J
       ENDIF
       IF(SQRT(DSQT) .GT. RTMX)THEN
         RTMX       = SQRT(DSQT)
         ITMX       = I
         JTMX       = J
       ENDIF
  200 CONTINUE
C
      RDL2         = SQRT(SUMD)
      RTL2         = SQRT(SUMT)
C
      CALL SOLUV(R,E,U,V,P,RNPO(-1,-1,1),RNPO(-1,-1,2),RNPO(-1,-1,3),
     1           RNPO(-1,-1,4))
C
      RETURN
      END
C
      SUBROUTINE STEP1(KDI,KDJ,LKD,KDP,IJP,R,E,U,V,T,P,QJ,
     1      RS,RUS,RVS,RES,US,VS,DRS,RNPO)
C
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,MI,MJ,MIJ2,IQS,IQE,JQS,JQE
      COMMON /PNTI / LKDI,LKDJ,LLKD,LKDP,MDA,MDB,NDIG
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILMT
      COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      COMMON /TIME0/ ITER,TIME,DT,DTI,DTJ,DTCFL,
     1               CFL,CFL1,CFL2,NUP,DTFIX
      COMMON /RES01/ RDMX,IDMX,JDMX,RTMX,ITMX,JTMX,RDL2,RTL2
      COMMON /UNST1/ LSTD,XMS,XST,MTIM,TPP(70)
C
      DIMENSION IJP(MXBP),KDI(MDA),KDJ(MDA),LKD(MDB),KDP(MDB)
      DIMENSION   QJ(-1:NI2,-1:NJ2)
      DIMENSION    R(-1:NI2,-1:NJ2),     E(-1:NI2,-1:NJ2)
      DIMENSION    U(-1:NI2,-1:NJ2),     V(-1:NI2,-1:NJ2)
      DIMENSION    T(-1:NI2,-1:NJ2),     P(-1:NI2,-1:NJ2)
      DIMENSION  DRS(-1:NI2,-1:NJ2,5)
      DIMENSION  RUS(-1:NI2,-1:NJ2,5), RVS(-1:NI2,-1:NJ2,5)
      DIMENSION   RS(-1:NI2,-1:NJ2,5), RES(-1:NI2,-1:NJ2,5)
      DIMENSION   US(-1:NI2,-1:NJ2,5),  VS(-1:NI2,-1:NJ2,5)
      DIMENSION RNPO(-1:NI2,-1:NJ2,4)
C
       DO 10 IJ      = -1,MIJ2
        DRS(IJ,-1,1) = 0.
        DRS(IJ,-1,2) = 0.
        DRS(IJ,-1,3) = 0.
        DRS(IJ,-1,4) = 0.
        DRS(IJ,-1,5) = 0.
10     CONTINUE
      IF(MTHR .EQ. 1) THEN
       CALL BMUPW(RS,US,VS,DRS)
      ELSEIF(MTHR .EQ. 2)THEN
cmic$ do all shared(qj) autoscope
        DO 12 NE  = 1,5
        CALL FVJTVD2(VS(-1,-1,NE),QJ,RS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVITVD2(US(-1,-1,NE),QJ,RS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
12      CONTINUE
      ELSEIF(MTHR .EQ. 3)THEN
cmic$ do all shared(qj) autoscope
        DO 14 NE  = 1,5
        CALL FVJEN2(VS(-1,-1,NE),QJ, RS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVIEN2(US(-1,-1,NE),QJ, RS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
14      CONTINUE
      ENDIF
       DO 16 IJ    =  -1,MIJ2
        RS(IJ,-1,1)  = RS(IJ,-1,1) + DRS(IJ,-1,1)
        RS(IJ,-1,2)  = RS(IJ,-1,2) + DRS(IJ,-1,2)
        RS(IJ,-1,3)  = RS(IJ,-1,3) + DRS(IJ,-1,3)
        RS(IJ,-1,4)  = RS(IJ,-1,4) + DRS(IJ,-1,4)
        RS(IJ,-1,5)  = RS(IJ,-1,5) + DRS(IJ,-1,5)
16     CONTINUE
       DO 20 IJ      = -1,MIJ2
        DRS(IJ,-1,1) = 0.
        DRS(IJ,-1,2) = 0.
        DRS(IJ,-1,3) = 0.
        DRS(IJ,-1,4) = 0.
        DRS(IJ,-1,5) = 0.
20     CONTINUE
      IF(MTHR .EQ. 1) THEN
       CALL BMUPW(RUS,US,VS,DRS)
      ELSEIF(MTHR .EQ. 2) THEN
cmic$ do all shared(qj) autoscope
        DO 22 NE  = 1,5
        CALL FVJTVD2(VS(-1,-1,NE),QJ,RUS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVITVD2(US(-1,-1,NE),QJ,RUS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
22      CONTINUE
      ELSEIF(MTHR .EQ. 3) THEN
cmic$ do all shared(qj) autoscope
        DO 24 NE  = 1,5
        CALL FVJEN2(VS(-1,-1,NE),QJ, RUS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVIEN2(US(-1,-1,NE),QJ, RUS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
24      CONTINUE
      ENDIF
       DO 26 IJ    =  -1,MIJ2
        RUS(IJ,-1,1)  = RUS(IJ,-1,1) + DRS(IJ,-1,1)
        RUS(IJ,-1,2)  = RUS(IJ,-1,2) + DRS(IJ,-1,2)
        RUS(IJ,-1,3)  = RUS(IJ,-1,3) + DRS(IJ,-1,3)
        RUS(IJ,-1,4)  = RUS(IJ,-1,4) + DRS(IJ,-1,4)
        RUS(IJ,-1,5)  = RUS(IJ,-1,5) + DRS(IJ,-1,5)
26     CONTINUE
       DO 30 IJ      = -1,MIJ2
        DRS(IJ,-1,1) = 0.
        DRS(IJ,-1,2) = 0.
        DRS(IJ,-1,3) = 0.
        DRS(IJ,-1,4) = 0.
        DRS(IJ,-1,5) = 0.
30     CONTINUE
      IF(MTHR .EQ. 1) THEN
       CALL BMUPW(RVS,US,VS,DRS)
      ELSEIF(MTHR .EQ. 2) THEN
cmic$ do all shared(qj) autoscope
        DO 32 NE  = 1,5
        CALL FVJTVD2(VS(-1,-1,NE),QJ,RVS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVITVD2(US(-1,-1,NE),QJ,RVS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
32      CONTINUE
      ELSEIF(MTHR .EQ. 3) THEN
cmic$ do all shared(qj) autoscope
        DO 34 NE  = 1,5
        CALL FVJEN2(VS(-1,-1,NE),QJ, RVS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVIEN2(US(-1,-1,NE),QJ, RVS(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
34      CONTINUE
      ENDIF
       DO 36 IJ    =  -1,MIJ2
        RVS(IJ,-1,1)  = RVS(IJ,-1,1) + DRS(IJ,-1,1)
        RVS(IJ,-1,2)  = RVS(IJ,-1,2) + DRS(IJ,-1,2)
        RVS(IJ,-1,3)  = RVS(IJ,-1,3) + DRS(IJ,-1,3)
        RVS(IJ,-1,4)  = RVS(IJ,-1,4) + DRS(IJ,-1,4)
        RVS(IJ,-1,5)  = RVS(IJ,-1,5) + DRS(IJ,-1,5)
36     CONTINUE
       DO 40 IJ      = -1,MIJ2
        DRS(IJ,-1,1) = 0.
        DRS(IJ,-1,2) = 0.
        DRS(IJ,-1,3) = 0.
        DRS(IJ,-1,4) = 0.
        DRS(IJ,-1,5) = 0.
40     CONTINUE
      IF(MTHR .EQ. 1) THEN
       CALL BMUPW(RES,US,VS,DRS)
      ELSEIF(MTHR .EQ. 2) THEN
cmic$ do all shared(qj) autoscope
        DO 42 NE  = 1,5
        CALL FVJTVD2(VS(-1,-1,NE),QJ,RES(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVITVD2(US(-1,-1,NE),QJ,RES(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
42      CONTINUE
      ELSEIF(MTHR .EQ. 3) THEN
cmic$ do all shared(qj) autoscope
        DO 44 NE  = 1,5
        CALL FVJEN2(VS(-1,-1,NE),QJ, RES(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
        CALL FVIEN2(US(-1,-1,NE),QJ, RES(-1,-1,NE),
     1               DRS(-1,-1,NE),NE)
44      CONTINUE
      ENDIF
       DO 46 IJ    =  -1,MIJ2
        RES(IJ,-1,1)  = RES(IJ,-1,1) + DRS(IJ,-1,1)
        RES(IJ,-1,2)  = RES(IJ,-1,2) + DRS(IJ,-1,2)
        RES(IJ,-1,3)  = RES(IJ,-1,3) + DRS(IJ,-1,3)
        RES(IJ,-1,4)  = RES(IJ,-1,4) + DRS(IJ,-1,4)
        RES(IJ,-1,5)  = RES(IJ,-1,5) + DRS(IJ,-1,5)
46     CONTINUE
C
      DO 100 IJ = -1,MIJ2
       RNPO(IJ,-1,1) = QJ(IJ,-1)*( RS(IJ,-1,1)+RS(IJ,-1,2)+
     1                RS(IJ,-1,3)+RS(IJ,-1,4)+RS(IJ,-1,5) )
       RNPO(IJ,-1,2) = QJ(IJ,-1)*( RUS(IJ,-1,1)+RUS(IJ,-1,2)+
     1                RUS(IJ,-1,3)+RUS(IJ,-1,4)+RUS(IJ,-1,5) )
       RNPO(IJ,-1,3) = QJ(IJ,-1)*( RVS(IJ,-1,1)+RVS(IJ,-1,2)+
     1                RVS(IJ,-1,3)+RVS(IJ,-1,4)+RVS(IJ,-1,5) )
       RNPO(IJ,-1,4) = QJ(IJ,-1)*( RES(IJ,-1,1)+RES(IJ,-1,2)+
     1                RES(IJ,-1,3)+RES(IJ,-1,4)+RES(IJ,-1,5) )
100   CONTINUE
C
      SUMD         = 0.
      SUMT         = 0.
      RDMX         = -1.
      RTMX         = -1.
      DO 200 I     = 2,NI-1
      DO 200 J     = 2,NJ-1
       RNP1         = RNPO(I,J,1)
       RUNP1        = RNPO(I,J,2)
       RVNP1        = RNPO(I,J,3)
       RENP1        = RNPO(I,J,4)
C
       XU           = 1./SQRT(ABS(1.-U(I,J)**2-V(I,J)**2))
       XU2          = XU*XU
       EP           = E(I,J) + P(I,J)
       DSQD         = (RNP1 - R(I,J)*XU)**2
       DSQT         = DSQD + (RUNP1 - XU2*EP*U(I,J))**2
     1                     + (RVNP1 - XU2*EP*V(I,J))**2
     1                     + (RENP1 - XU2*EP+P(I,J))**2
       SUMD         = SUMD + DSQD
       SUMT         = SUMT + DSQT
       IF(SQRT(DSQD) .GT. RDMX)THEN
         RDMX       = SQRT(DSQD)
         IDMX       = I
         JDMX       = J
       ENDIF
       IF(SQRT(DSQT) .GT. RTMX)THEN
         RTMX       = SQRT(DSQT)
         ITMX       = I
         JTMX       = J
       ENDIF
  200 CONTINUE
C
      RDL2         = SQRT(SUMD)
      RTL2         = SQRT(SUMT)
C
      CALL SOLUV(R,E,U,V,P,RNPO(-1,-1,1),RNPO(-1,-1,2),RNPO(-1,-1,3),
     1           RNPO(-1,-1,4))
C
      RETURN
      END
