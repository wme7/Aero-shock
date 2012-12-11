      PROGRAM B2D
C
C     FINITE DIFFERENCE SOLVER FOR BGK TYPE BOLTZMANN EQUATION
C
C       DEVELOPED BY  HJC
C       VERSION X  START AT 1994 03 16
C                  FOR ARBITRARY BOUNDARY SET
C
C       MAXWK : MAX WORKING SPACE
C       MAXIJP: MAX POINT IN I OR J DIRECTION
C       MAXBC : MAX B.C. SETTING
C       MAXGRD: MAX GRID POINT
C       MAXLIN: MAX LU SWEEP LINE (DIAG. LINE)
C       MAXVEL: MAX VELOCITY POINTS IN PHASE PLACE
C       MAXBCP: MAX B.C. IDEX SPACE
C       MAXLUI: MAX LU IDEX SPACE
C
      PARAMETER (MAXWK=32000000,MAXIJP=305,MAXBC=20)
C     PARAMETER (MAXWK=11500000,MAXIJP=305,MAXBC=20)
      PARAMETER (MAXGRD=305*155,MAXLIN=305+155,MAXVEL=8000)
      PARAMETER (MAXBCP=305*10*2*2, MAXLUI=MAXGRD*2+MAXLIN*2+1)
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION S(MAXWK),MLU(MAXLUI),IJP(MAXBCP)
C
      COMMON /TAPE / INAME,IGRID,IWRIT,IREAD,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXB,MAXG,MAXL,MAXV,MAXBP,MAXLI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /NEXT / NEXTLOC
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /FINP0/ IQX,IQY,IOX,IOY,MVX,MVY,DVX,DVY,VXS,VYS
      COMMON /FINP1/ KMOD,MMOD,ZETA,PR,VISS,XLDA
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR
      COMMON /FINP3/ ALPHA,XKNINF,TSTM,XMINF
      COMMON /CNST0/ DM,VISREF,VISTEM,VISFAC,BO,AVO,RGAS,TREF
      COMMON /CNST1/ TINF,XMPINF,XMUINF,XNUINF,DNINF,SSINF,UDINF,
     1               RHOINF,PINF,EINF,CINF,REINF
      COMMON /CNST2/ D0,T0,P0,U0,V0,TAU0,GAMMA,PI,SPI
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST0/ URGHT,VRGHT,DRGHT,TRGHT,PRGHT,
     1               ULEFT,VLEFT,DLEFT,TLEFT,PLEFT
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /IMPL0/ NSOU
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
      COMMON /BC02 / IDA,ISA,ICA,IAS,IAE,IPA,LENA,
     1               IDB,ISB,ICB,IBS,IBE,IPB,LENB
      COMMON /PNT  / LVX,LVY,LWX,LWY,LCXI,LCET,LCXP,LCXM,LCEP,LCEM,
     &               LX,LY,LXIX,LXIY,LETAX,LETAY,LJAC,LDT,
     &               LG,LH,LD,LT,LU,LV,LP,LXM,LQX,LQY,LA1,LA2,LDD,
     &               LFG,LFH,LRG,LRH,LDQG,LDQH,LTXX,LTYY,LTXY,LCP,
     &                LSXX,LSYY,LSZZ,LSXY
      COMMON /PNTI / LKDI,LKDJ,LLKD,LKDP,MDA,MDB,NDIG
      COMMON /RES0/ SUMD,Z1MAX,M1DI,M1DJ,Z2MAX,M2DI,M2DJ,SUMDG,SUMDH
      COMMON /EQUL / ISOU,IEQUL
C
      DIMENSION RES(10000,6)
C
      MAXW  = MAXWK
      MAXI  = MAXIJP
      MAXB  = MAXBC
      MAXG  = MAXGRD
      MAXL  = MAXLIN
      MAXV  = MAXVEL
      MAXBP = MAXBCP
      MAXLI = MAXLUI
C
      KT      = 1
      NEXTLOC = 1
C
      IGRID  = 11
      IVGRD  = 22
      INAME  = 12
      IWRIT  = 13
      IREAD  = 14
      IFORC  = 15
      IRSTA  = 16
      IRESD  = 18
      ISHOW  = 6
	ICONV  = 66
      OPEN(IGRID,FILE='B11.GRD',STATUS='UNKNOWN')
      OPEN(INAME,FILE='B12.NAM',STATUS='UNKNOWN')
      OPEN(IFORC,FILE='B15.FOC',STATUS='UNKNOWN')
      OPEN(IRSTA,FILE='B16.RST',STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(   17,FILE='B17.WAL',STATUS='UNKNOWN')
      OPEN(IRESD,FILE='B18.RES',STATUS='UNKNOWN')
      OPEN(   67,FILE='B67.VND',STATUS='UNKNOWN')
      OPEN(   69,FILE='B69.BND',STATUS='UNKNOWN')
      OPEN(   71,FILE='B70.TEC',STATUS='UNKNOWN')
      OPEN(   66,FILE='B66.TEC',STATUS='UNKNOWN')
C
      READ(IGRID,*)NI,NJ
C
      CALL SETUP(IJP(1))
C
      IF(IQX .EQ. 5 .AND. IQY .EQ. 5)THEN
      OPEN(IVGRD,FILE='B22.VGR',STATUS='UNKNOWN')
      READ(IVGRD,*)MVX,MVY
      NVX   = MVX
      NVY   = MVY
      MXV   = NVX*NVY
      ENDIF
C
      CALL POINT
C
      CALL LUPNT(MLU(LKDI),MLU(LKDJ),MLU(LLKD),MLU(LKDP))
C
      CALL VELDAT(IQX,IOX,NVX,S(LWX),S(LVX),DVX,VXS,2,MVX,IVGRD)
      CALL VELDAT(IQY,IOY,NVY,S(LWY),S(LVY),DVY,VYS,2,MVY,IVGRD)
C
      CALL VELCOE(S(LVX),S(LVY),S(LWX),S(LWY),MVX,MVY,MXV,NVX,NVY)
C
      CALL GRIDIN(S(LX),S(LY),IJP(1))
C
      CALL METJAC(S(LX),S(LY),S(LXIX),S(LXIY),S(LETAX),S(LETAY),
     1            S(LJAC),IJP(1))
C
      CALL DTIME(S(LDT),S(LVX),S(LVY),
     1           S(LXIX),S(LXIY),S(LETAX),S(LETAY),S(LJAC))
C
         ISAV   = 1
      IF(IRSTRT .LT. 0) THEN
         ISAV   = 2
         IRSTRT = ABS(IRSTRT)
      ENDIF
      IF(IEQUL .EQ. 1) ISAV = 2
C       
      IF(IRSTRT.EQ.1)THEN
        READ(IRSTA)KS2D,NI,NJ,NVX,NVY
        CALL XREAD(S(LG),KS2D,IRSTA)
        CALL XREAD(S(LH),KS2D,IRSTA)
        READ(IRSTA)TIME,ITER0,KT
        CALL INITQ(S(LX),S(LD),S(LT),S(LU),S(LV),S(LP),S(LQX),S(LQY))
      ELSEIF(IRSTRT.EQ.2)THEN
        CALL INITQ(S(LX),S(LD),S(LT),S(LU),S(LV),S(LP),S(LQX),S(LQY))
        CALL QRST(S(LD),S(LT),S(LU),S(LV),S(LP),S(LXM),4,ITER0)
        CALL RSTGH(S(LVX),S(LVY),S(LD),S(LU),S(LV),S(LT),S(LP),
     1           S(LQX),S(LQY),S(LG),S(LH))
        CALL BCG(IJP(1),
     1           S(LX),S(LY),S(LXIX),S(LXIY),S(LETAX),S(LETAY),
     1           S(LJAC),S(LD),S(LU),S(LV),
     1           S(LT),S(LG),S(LH),
     1           S(LVX),S(LVY),S(LWX),S(LWY),S(LA1))
      ELSE
        CALL INITQ(S(LX),S(LD),S(LT),S(LU),S(LV),S(LP),S(LQX),S(LQY))
        CALL INITG(S(LX),S(LVX),S(LVY),S(LD),S(LT),S(LU),S(LV),S(LP),
     1             S(LQX),S(LQY),S(LG),S(LH),S(LJAC))
        TIME = 0.
        ITER0 = 0
      ENDIF
      CALL MACRO(S(LDT),S(LVX),S(LVY),S(LWX),S(LWY),
     1             S(LX),S(LY),S(LXIX),S(LXIY),S(LETAX),S(LETAY),
     1             S(LJAC),S(LD),S(LU),S(LV),
     1             S(LT),S(LP),S(LQX),S(LQY),S(LG),S(LH),IJP(1),
     1             S(LSXX),S(LSYY),S(LSZZ),S(LSXY))
C
      CALL TECPLT(S(LX),S(LY),S(LD),S(LU),S(LV),S(LT),S(LP),S(LXM),
     &            S(LQX),S(LQY),71)
C
      ITER = ITER0+1
 1000 CONTINUE
C
      LTP  = 0
      SUMDG = 0.
      SUMDH = 0.
C
      CFL  = (CFL2 - CFL1)/NUP * ITER+ CFL1
      IF(CFL .GE. CFL2) CFL = CFL2
      DTCFL = DT*CFL
      TIME = TIME + DTCFL
      IF(ILOCDT.EQ.0)THEN
        CALL LDTSET(S(LDT))
      ENDIF
      IF(TIME.GE.TPP(KT).AND.LSTD.EQ.0)THEN
        DTCFL = TPP(KT) - (TIME - DTCFL)
        DTT   = DT
        DT    = DTCFL/CFL
        CALL  LDTSET(S(LDT))
        DT    = DTT
        LTP   = 1
        TIME  = TPP(KT)
        KT    = KT + 1
      ENDIF
C
      DO 200 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IF(IBTYPE .EQ. 2 .OR. IBTYPE . EQ. 3) THEN
      IDA    = IDBC( 2,IBC)
      ISA    = IDBC( 3,IBC)
      ICA    = IDBC( 4,IBC)
      IAS    = IDBC( 5,IBC)
      IAE    = IDBC( 6,IBC)
      IPG    = IDBC( 7,IBC)
      IPH    = IDBC( 8,IBC)
      IPA    = IDBC(25,IBC)
      LENA   = (IAE - IAS)*ICA + 1
      CALL PUTBC(IJP(IPA+1),S(LG),S(LH),S(IPG),S(IPH))
      ENDIF
  200 CONTINUE
C
      CALL STEP(S(LDT),S(LVX),S(LVY),S(LWX),S(LWY),S(LX),S(LY),
     1          S(LXIX),S(LXIY),S(LETAX),S(LETAY),S(LJAC),
     1          S(LG),S(LH),S(LRG),S(LRH),S(LDQG),S(LDQH),
     1          S(LD),S(LT),S(LU),S(LV),S(LP),S(LQX),S(LQY),
     1          S(LCXI),S(LCET),S(LCXP),S(LCXM),S(LCEP),S(LCEM),
     1          S(LFG),S(LFH),S(LA1),S(LA2),S(LDD),
     1          MLU(LKDI),MLU(LKDJ),MLU(LLKD),MLU(LKDP),
     1          IJP(1),S(1),S(LSXX),S(LSYY),S(LSZZ),S(LSXY))
C
      CALL MACRO(S(LDT),S(LVX),S(LVY),S(LWX),S(LWY),
     1           S(LX),S(LY),S(LXIX),S(LXIY),S(LETAX),S(LETAY),
     1           S(LJAC),S(LD),S(LU),S(LV),
     1           S(LT),S(LP),S(LQX),S(LQY),S(LG),S(LH),IJP(1),
     1           S(LSXX),S(LSYY),S(LSZZ),S(LSXY))
C
      IF(IEQUL .EQ. 1) THEN
        CALL EQULG(S(LVX),S(LVY),S(LD),S(LT),S(LU),S(LV),S(LP),
     1             S(LQX),S(LQY),S(LG),S(LH),S(LJAC))
        CALL BCG(IJP(1),
     1           S(LX),S(LY),S(LXIX),S(LXIY),S(LETAX),S(LETAY),
     1           S(LJAC),S(LD),S(LU),S(LV),
     1           S(LT),S(LG),S(LH),
     1           S(LVX),S(LVY),S(LWX),S(LWY),S(LA1))
      ENDIF
      RES(ITER,1) = DTCFL
      RES(ITER,2) = SQRT(SUMD)
      RES(ITER,3) = Z1MAX
      RES(ITER,4) = Z2MAX
      RES(ITER,5) = SQRT(SUMDG)
      RES(ITER,6) = SQRT(SUMDH)
C
      WRITE(*,11)ITER,DTCFL,TIME,M1DI,M1DJ,LOG10(Z1MAX),
     1           M2DI,M2DJ,LOG10(Z2MAX),LOG10(RES(ITER,2)),
     1           LOG10(RES(ITER,5)),LOG10(RES(ITER,6))
      WRITE(66,11)ITER,DTCFL,TIME,M1DI,M1DJ,LOG10(Z1MAX),
     1           M2DI,M2DJ,LOG10(Z2MAX),LOG10(RES(ITER,2)),
     1           LOG10(RES(ITER,5)),LOG10(RES(ITER,6))
   11 FORMAT(I5,1X,F8.5,1X,F7.3,2(1X,2I3,1X,F7.3),3(1X,F7.3))
C
      IF((MOD(ITER,IPRINT).EQ.0).OR.(ITER.EQ.NSTEP).OR.
     1    LTP.EQ.1) THEN
C
      IF(LTP.EQ.1) IQFIL1 = (KT-1) + 30
C
      CALL TECPLT(S(LX),S(LY),S(LD),S(LU),S(LV),S(LT),S(LP),S(LXM),
     &            S(LQX),S(LQY),71)
C
      CALL FORCE(S(LX),S(LY),S(LD),S(LU),S(LV),S(LT),S(LP),
     &           S(LG),IJP(1))
C
      REWIND (IRSTA)
      REWIND (IRESD)
C
      IF(ISAV .EQ. 1) THEN
        WRITE(IRSTA)KS2D,NI,NJ,NVX,NVY
        CALL XWRITE(S(LG),KS2D,IRSTA)
        CALL XWRITE(S(LH),KS2D,IRSTA)
        WRITE(IRSTA)TIME,ITER,KT
      ENDIF
C
      DO 100 IT = 1,ITER
      WRITE(IRESD,12) IT,RES(IT,1),RES(IT,2),RES(IT,3),
     1                   RES(IT,4),RES(IT,5),RES(IT,6)
  100 CONTINUE
   12 FORMAT(1X,I10,6(2X,E15.6))
      WRITE(IRESD,*)' TIME = ',TIME
C     CLOSE (IRESD)
C
C  ENDIF OF IPRINT
C
C     CLOSE (IRSTA)
C
      ENDIF
C
      IF(ITER.GE.NSTEP) STOP
      ITER = ITER + 1
      GO TO 1000
      END
      SUBROUTINE BCDQ(LU,CXI,CET,DQG,DQH,IJP)
      PARAMETER (MAXBC=20)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
      COMMON /AXIS / IAXIS
      DIMENSION  IJP(MXBP)
      DIMENSION  DQG(-2:NI3,-2:NJ3), DQH(-2:NI3,-2:NJ3)
      DIMENSION  CXI(-2:NI3,-2:NJ3), CET(-2:NI3,-2:NJ3)
C
      DO 100 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IDA    = IDBC( 2,IBC)
      ISA    = IDBC( 3,IBC)
      ICA    = IDBC( 4,IBC)
      IAS    = IDBC( 5,IBC)
      IAE    = IDBC( 6,IBC)
      IPA    = IDBC(25,IBC)
      LENA   = (IAE - IAS)*ICA + 1
      ISGN   = -1
      IF(ISA .EQ. 1) ISGN = 1
      IF(ISA .NE. LU) GO TO 100
      IF(IBTYPE .EQ. 0)THEN
C
C     SETTING DQG AND DQH AT 0 POINT 
C       1. NOT SPECIFY BOUNDARY
C       2. Equal Upstream           ------------- Now try this.
C
      DO 130 II = 1,LENA
      I2A   =  IJP(IPA+LENA*2+II)
      J2A   =  IJP(IPA+LENA*3+II)
      I3A   =  IJP(IPA+LENA*4+II)
      J3A   =  IJP(IPA+LENA*5+II)
C
C     IF(IDA .EQ. 1) THEN
C       CC0 =  CXI(I2A,J2A) * ISGN
C     ELSE
C       CC0 =  CET(I2A,J2A) * ISGN
C     ENDIF
C     IF(LU .EQ. 2) THEN
C     DQG(I2A,J2A) = DQG(I3A,J3A)
C     DQH(I2A,J2A) = DQH(I3A,J3A)
C     ELSE
C
      DQG(I2A,J2A) = 0.
      DQH(I2A,J2A) = 0.
C     ENDIF
  130 CONTINUE
      ELSEIF(IBTYPE .EQ. 2 .OR. IBTYPE .EQ. 3) THEN
C
C     SETTING DQG AND DQH AT 0 POINT 
C         1.  NOT SPECIFY BOUNDARY
C         2.  SET DQ AT 0 AND NI1 POINTS EQUAL 0  ---- NOW TRY
C
      IF(IAXIS .EQ. 0 ) THEN
      DO 140 II = 1,LENA
          I2A   =  IJP(IPA+LENA*2+II)
          J2A   =  IJP(IPA+LENA*3+II)
      DQG(I2A,J2A) = 0.
      DQH(I2A,J2A) = 0.
  140 CONTINUE
      ELSE
      DO 145 II = 1,LENA
          I2A   =  IJP(IPA+LENA*2+II)
          J2A   =  IJP(IPA+LENA*3+II)
          I3A   =  IJP(IPA+LENA*4+II)
          J3A   =  IJP(IPA+LENA*5+II)
      DQG(I2A,J2A) = 0.
      DQH(I2A,J2A) = 0.
      DQG(I3A,J3A) = 0.
      DQH(I3A,J3A) = 0.
  145 CONTINUE
      ENDIF
      ELSEIF(IBTYPE .EQ. 4)THEN
      IDB    = IDBC(12,IBC)
      ISB    = IDBC(13,IBC)
      ICB    = IDBC(14,IBC)
      IBS    = IDBC(15,IBC)
      IBE    = IDBC(16,IBC)
      IPB    = IDBC(26,IBC) + 1
      LENB   = (IBE - IBS)*ICB + 1
      DO 110 II = 1,LENA
          I2A   =  IJP(IPA+LENA*2+II)
          J2A   =  IJP(IPA+LENA*3+II)
          I3A   =  IJP(IPA+LENA*4+II)
          J3A   =  IJP(IPA+LENA*5+II)
          I4A   =  IJP(IPA+LENA*6+II)
          J4A   =  IJP(IPA+LENA*7+II)
          I2B   =  IJP(IPB+LENA*2+II)
          J2B   =  IJP(IPB+LENA*3+II)
          I3B   =  IJP(IPB+LENA*4+II)
          J3B   =  IJP(IPB+LENA*5+II)
          I4B   =  IJP(IPB+LENA*6+II)
          J4B   =  IJP(IPB+LENA*7+II)
C     DQG(I2A,J2A) = DQG(I4B,J4B)
C     DQG(I2B,J2B) = DQG(I4A,J4A)
C     DQH(I2A,J2A) = DQH(I4B,J4B)
C     DQH(I2B,J2B) = DQH(I4A,J4A)
      DQG(I2A,J2A) = 0.
      DQG(I2B,J2B) = 0.
      DQH(I2A,J2A) = 0.
      DQH(I2B,J2B) = 0.
C
C  Test at 1994 7 21 for airfoil
C
      DQG(I3A,J3A) = 0.
      DQG(I3B,J3B) = 0.
      DQH(I3A,J3A) = 0.
      DQH(I3B,J3B) = 0.
  110 CONTINUE
      ELSEIF(IBTYPE .EQ. 1)THEN
C
C   TRY SETTING DQ BC
C   1. NOT SPECIFY ANYTHING
C   2. DQG(IA2,JA2) = 0 AND DQG(IA3,JA3)=0 FOR INEFFECT --- NOW
C
      DO 120 II = 1,LENA
          I2A   =  IJP(IPA+LENA*2+II)
          J2A   =  IJP(IPA+LENA*3+II)
          I3A   =  IJP(IPA+LENA*4+II)
          J3A   =  IJP(IPA+LENA*5+II)
C
C     IF(IDA .EQ. 1) THEN
C       CC2 =  CXI(I3A,J3A) * ISGN
C     ELSE
C       CC2 =  CET(I3A,J3A) * ISGN
C     ENDIF
C     IF(CC2 .GE. 0.) THEN
C for bc1
C       DQG(I2A,J2A) = 0.
C       DQH(I2A,J2A) = 0.
C       DQG(I3A,J3A) = 0.
C       DQH(I3A,J3A) = 0.
C for bc2
        DQG(I2A,J2A) = 0.
        DQH(I2A,J2A) = 0.
C for bc3
C       DQG(I2A,J2A) = 0.
C       DQH(I2A,J2A) = 0.
C for bc4
C       DQG(I2A,J2A) = 0.
C       DQH(I2A,J2A) = 0.
C       DQG(I3A,J3A) = 0.
C       DQH(I3A,J3A) = 0.
C     ENDIF
  120 CONTINUE
      ELSEIF(IBTYPE .EQ. 5)THEN
C
      DO 150 II = 1,LENA
          I2A   =  IJP(IPA+LENA*2+II)
          J2A   =  IJP(IPA+LENA*3+II)
          I3A   =  IJP(IPA+LENA*4+II)
          J3A   =  IJP(IPA+LENA*5+II)
        DQG(I2A,J2A) = 0.
        DQH(I2A,J2A) = 0.
        DQG(I3A,J3A) = 0.
        DQH(I3A,J3A) = 0.
  150 CONTINUE
      ENDIF
  100 CONTINUE
      RETURN
      END
C
      SUBROUTINE BCG(IJP,X,Y,XIX,XIY,ETAX,ETAY,QJ,D,U,V,T,
     1              G,H,A,B,WA,WB,SA1)
      PARAMETER (MAXBC=20)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
      COMMON /BC02 / IDA,ISA,ICA,IAS,IAE,IPA,LENA,
     1               IDB,ISB,ICB,IBS,IBE,IPB,LENB
      COMMON /CNST2/ D0,T0,P0,U0,V0,TAU0,GAMMA,PI,SPI
      COMMON /AXIS / IAXIS
      DIMENSION  IJP(MXBP)
      DIMENSION    X(-2:NI3,-2:NJ3),   Y(-2:NI3,-2:NJ3)
      DIMENSION  XIX(-2:NI3,-2:NJ3), XIY(-2:NI3,-2:NJ3)
      DIMENSION ETAX(-2:NI3,-2:NJ3),ETAY(-2:NI3,-2:NJ3)
      DIMENSION   QJ(-2:NI3,-2:NJ3)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION    G(NI,NJ,NVX,NVY),   H(NI,NJ,NVX,NVY)
      DIMENSION  SA1(-2:NI3,-2:NJ3)
      DIMENSION    A(NVX),B(NVY),WA(NVX),WB(NVY)
C
      DO 100 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IDA    = IDBC( 2,IBC)
      ISA    = IDBC( 3,IBC)
      ICA    = IDBC( 4,IBC)
      IAS    = IDBC( 5,IBC)
      IAE    = IDBC( 6,IBC)
      IPG    = IDBC( 7,IBC)
      IPH    = IDBC( 8,IBC)
      IPA    = IDBC(25,IBC)
      LENA   = (IAE - IAS)*ICA + 1
      ISGN = -1
      IF(ISA .EQ. 1) ISGN = 1
C
      IF(IBTYPE .EQ. 0 ) THEN
        DO 110 II = 1,LENA
        IC    = IJP(IPA + LENA*4 + II)
        JC    = IJP(IPA + LENA*5 + II)
        ID    = IJP(IPA + LENA*6 + II)
        JD    = IJP(IPA + LENA*7 + II)
        DO 115 KL   = 1,MXV
         K      = MOD(KL-1,NVX)+1
         L      = (KL-1)/NVX + 1
         CX     = A(K)
         CY     = B(L)
        IF(IDA .EQ. 1) THEN
          CC        = ISGN *(CX * XIX(IC,JC) + CY* XIY(IC,JC) )
        ELSE
          CC        = ISGN *(CX *ETAX(IC,JC) + CY*ETAY(IC,JC) )
        ENDIF
        IF(CC .LT. 0.) THEN
          G(IC,JC,K,L)  =  G(ID,JD,K,L)
          H(IC,JC,K,L)  =  H(ID,JD,K,L)
        ELSE
C
C     FARFIELD G AND H SETTING 
C        1.  WITH LOCAL MAXWELLIAN    ----------TRY NOW
C        2.  FREE STREAM MAXWELLIAN   
C
          UC                    = U(IC,JC)
          VC                    = V(IC,JC)
C         IF(IDA .EQ. 1) THEN
C         CCU       = ISGN *(UC * XIX(IC,JC) + VC* XIY(IC,JC) )
C         ELSE
C         CCU       = ISGN *(UC *ETAX(IC,JC) + VC*ETAY(IC,JC) )
C         ENDIF
C         IF(CCU .LT. 0.) THEN
          TC                    = T(IC,JC)
          TUEXP                 = EXP(-((CX-UC)*(CX-UC)
     1                             +(CY-VC)*(CY-VC))/TC)
          G(IC,JC,K,L)          = D(IC,JC)/(PI*TC)*TUEXP
          H(IC,JC,K,L)          = 0.5*TC*G(IC,JC,K,L)
C         ELSE
C         TUEXP                 = EXP(-((CX-U0)*(CX-U0)
C    1                             +(CY-V0)*(CY-V0))/T0)
C         G(IC,JC,K,L) = D0/(PI*T0)*EXP(TUEXP)
C         H(IC,JC,K,L) = 0.5 * T0 * G(IC,JC,K,L)
C         ENDIF
        ENDIF
  115 CONTINUE
  110 CONTINUE
C
      ELSEIF(IBTYPE .EQ. 1) THEN
        DO 120 II = 1,LENA
        IC    = IJP(IPA + LENA*4 + II)
        JC    = IJP(IPA + LENA*5 + II)
        ID    = IJP(IPA + LENA*6 + II)
        JD    = IJP(IPA + LENA*7 + II)
C
        DO 125 KL   = 1,MXV
         K      = MOD(KL-1,NVX)+1
         L      = (KL-1)/NVX + 1
         CX     = A(K)
         CY     = B(L)
        IF(IDA .EQ. 2) THEN
          TX      = ETAX(IC,JC)
          TY      = ETAY(IC,JC)
        ELSE
          TX      =  XIX(IC,JC)
          TY      =  XIY(IC,JC)
        ENDIF
        AREA      = SQRT(TX*TX + TY*TY)
        ENX       = TX / AREA
        ENY       = TY / AREA
        CC        = ISGN *(CX * ENX + CY* ENY )
        IF(CC .LE. 0.) THEN
          G(IC,JC,K,L)  =  G(ID,JD,K,L)
          H(IC,JC,K,L)  =  H(ID,JD,K,L)
        ENDIF
  125 CONTINUE
  120 CONTINUE
      ELSEIF(IBTYPE .EQ. 2 .OR. IBTYPE .EQ. 3) THEN
      IF(IAXIS .EQ. 1) THEN
        DO 130 II = 1,LENA
        IC    = IJP(IPA + LENA*4 + II)
        JC    = IJP(IPA + LENA*5 + II)
        ID    = IJP(IPA + LENA*6 + II)
        JD    = IJP(IPA + LENA*7 + II)
        IE    = IJP(IPA + LENA*8 + II)
        JE    = IJP(IPA + LENA*9 + II)
        DO 135 KL   = 1,MXV
         K      = MOD(KL-1,NVX)+1
         L      = (KL-1)/NVX + 1
         CX     = A(K)
         CY     = B(L)
        IF(L .LT. NVY/2+1) THEN
        FF           = (Y(IE,JE)/Y(ID,JD))**2
        G(IC,JC,K,L) = (FF*G(ID,JD,K,L)-G(IE,JE,K,L))/(FF-1.)
        H(IC,JC,K,L) = (FF*H(ID,JD,K,L)-H(IE,JE,K,L))/(FF-1.)
        ELSE
        LL = NVY - L + 1
        G(IC,JC,K,L) = G(IC,JC,K,LL)
        H(IC,JC,K,L) = H(IC,JC,K,LL)
        ENDIF
  135   CONTINUE
  130   CONTINUE
      ENDIF
      ELSEIF(IBTYPE .EQ. 4) THEN
      IDB    = IDBC( 12,IBC)
      ISB    = IDBC( 13,IBC)
      ICB    = IDBC( 14,IBC)
      IBS    = IDBC( 15,IBC)
      IBE    = IDBC( 16,IBC)
      IPB    = IDBC(26,IBC)
      LENB   = (IBE - IBS)*ICB + 1
        DO 140 II = 1,LENA
        IC1    = IJP(IPA + LENA*4 + II)
        JC1    = IJP(IPA + LENA*5 + II)
        ID1    = IJP(IPA + LENA*6 + II)
        JD1    = IJP(IPA + LENA*7 + II)
        IC2    = IJP(IPB + LENA*4 + II)
        JC2    = IJP(IPB + LENA*5 + II)
        ID2    = IJP(IPB + LENA*6 + II)
        JD2    = IJP(IPB + LENA*7 + II)
        DO 145 KL   = 1,MXV
         K      = MOD(KL-1,NVX)+1
         L      = (KL-1)/NVX + 1
         CX     = A(K)
         CY     = B(L)
C
C   Test at 1994 7 21 for airfoil
C
C       G(IC1,JC1,K,L) = 0.5*(G(IC1,JC1,K,L) + G(IC2,JC2,K,L))
C       H(IC1,JC1,K,L) = 0.5*(H(IC1,JC1,K,L) + H(IC2,JC2,K,L))
        G(IC1,JC1,K,L) = 0.5*(G(ID1,JD1,K,L) + G(ID2,JD2,K,L))
        H(IC1,JC1,K,L) = 0.5*(H(ID1,JD1,K,L) + H(ID2,JD2,K,L))
        G(IC2,JC2,K,L) = G(IC1,JC1,K,L)
        H(IC2,JC2,K,L) = H(IC1,JC1,K,L)
  145 CONTINUE
  140   CONTINUE
      ELSEIF(IBTYPE .EQ. 5) THEN
        DO 150 II = 1,LENA
        IC    = IJP(IPA + LENA*4 + II)
        JC    = IJP(IPA + LENA*5 + II)
        ID    = IJP(IPA + LENA*6 + II)
        JD    = IJP(IPA + LENA*7 + II)
        DO 155 KL   = 1,MXV
         K      = MOD(KL-1,NVX)+1
         L      = (KL-1)/NVX + 1
         CX     = A(K)
         CY     = B(L)
        G(IC,JC,K,L) = G(ID,JD,K,L)
        H(IC,JC,K,L) = H(ID,JD,K,L)
  155 CONTINUE
  150 CONTINUE
      ENDIF
C
  100 CONTINUE
      RETURN
      END
      SUBROUTINE BCQWAL(ID,EX,EY,QJ,D,U,V,T,P,QX,QY,G,H,A,B,WA,WB,
     1                  SLXX,SLYY,SLZZ,SLXY)
      PARAMETER (NV=8000)
      PARAMETER (MAXBC=20)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
      COMMON /BC02 / IDA,ISA,ICA,IAS,IAE,IPA,LENA,
     1               IDB,ISB,ICB,IBS,IBE,IPB,LENB
      COMMON /BC03 / TWAL,TACC,UACC
      COMMON /VELC / W0(NV),W1(NV),W2(NV),W3(NV),
     1               A3(NV),A4(NV),A5(NV),A6(NV)
      COMMON /FINP1/ KMOD,MMOD,ZETA,PR,VISS,XLDA
      DIMENSION  ID(LENA,10)
      DIMENSION  EX(-2:NI3,-2:NJ3),EY(-2:NI3,-2:NJ3)
      DIMENSION  QJ(-2:NI3,-2:NJ3)
      DIMENSION   QX(-2:NI3,-2:NJ3),  QY(-2:NI3,-2:NJ3)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION    P(-2:NI3,-2:NJ3)
      DIMENSION   G(NI,NJ,NVX,NVY),   H(NI,NJ,NVX,NVY)
      DIMENSION   A(NVX),B(NVY),WA(NVX),WB(NVY)
      DIMENSION SLXX(-2:NI3,-2:NJ3),SLYY(-2:NI3,-2:NJ3)
      DIMENSION SLZZ(-2:NI3,-2:NJ3),SLXY(-2:NI3,-2:NJ3)

      ISGN = -1
      IF(ISA .EQ. 1) ISGN = 1
      PI        = 3.1415926535
C
      DO 110 II = 1,LENA
      TX        = EX(ID(II,5),ID(II,6))
      TY        = EY(ID(II,5),ID(II,6))
      AREA      = SQRT(TX*TX + TY*TY)
      ENX       =  TX/AREA
      ENY       =  TY/AREA
      I1        =  ID(II,1)
      J1        =  ID(II,2)
      I2        =  ID(II,3)
      J2        =  ID(II,4)
      I3        =  ID(II,5)
      J3        =  ID(II,6)
      I4        =  ID(II,7)
      J4        =  ID(II,8)
      I5        =  ID(II,9)
      J5        =  ID(II,10)
C
      SD        = 0.
      DO 120 KL = 1,MXV
      K         = MOD(KL-1,NVX)+1
      L         = (KL-1)/NVX + 1
      CN0    = A(K) * ENX + B(L) *ENY
      IF(CN0 .LE. 0.)THEN
         SD     = SD  + CN0*WA(K)*WB(L)*G(I3,J3,K,L)
      ENDIF
  120 CONTINUE
C
      D(I3,J3) = -2.*SD*SQRT(PI/TWAL)
      DO 130 KL = 1,MXV
      K         = MOD(KL-1,NVX)+1
      L         = (KL-1)/NVX + 1
C
      CN0    = A(K) * ENX + B(L) *ENY
      IF(CN0 .GT. 0.)THEN
      AA1  = -(A(K)*A(K)+B(L)*B(L))/TWAL
      G(I3,J3,K,L) = D(I3,J3)/(PI*TWAL)*EXP(AA1)
      H(I3,J3,K,L) = 0.5*TWAL*G(I3,J3,K,L)
      ENDIF
  130 CONTINUE
C
C     CHECK
C
      IF(D(I3,J3).LE.0.)THEN
        WRITE(*,*)'ERROR : DENSITY < 0  (BCQWALL) (I,J)=',
     1              I3,J3,D(I3,J3)
        STOP
      ENDIF
C
      SD     = 0.
      SU     = 0.
      SV     = 0.
      ST     = 0.
      SQX    = 0.
      SQY    = 0.
      SAA    = 0.
      SBB    = 0.
      SAB    = 0.
      DO 230 KL = 1 , MXV
      K      = MOD(KL-1,NVX)+1
      L      = (KL-1)/NVX + 1
      SD     = SD   + W0(KL)*G(I3,J3,K,L)
      SU     = SU   + W1(KL)*G(I3,J3,K,L)
      SV     = SV   + W2(KL)*G(I3,J3,K,L)
      ST     = ST   + W3(KL)*G(I3,J3,K,L) + W0(KL)*H(I3,J3,K,L)
      AA0    = H(I3,J3,K,L)   + A5(KL)*G(I3,J3,K,L)
      SQX    = SQX  + W1(KL)*AA0
      SQY    = SQY  + W2(KL)*AA0
      SAA    = SAA  + A3(KL)*G(I3,J3,K,L)
      SBB    = SBB  + A4(KL)*G(I3,J3,K,L)
      SAB    = SAB  + A6(KL)*G(I3,J3,K,L)
  230 CONTINUE
      U(I3,J3) = SU / SD
      V(I3,J3) = SV / SD
      T(I3,J3) = (ST - SD*(U(I3,J3)*U(I3,J3) +
     1           V(I3,J3)*V(I3,J3)))/SD*2./3.
      P(I3,J3) = D(I3,J3) * T(I3,J3)
      UV2      = U(I3,J3) *U(I3,J3) + V(I3,J3)*V(I3,J3) - 1.5*T(I3,J3)
      QX(I3,J3)= SQX - 2.*(SAA*U(I3,J3) + SAB*V(I3,J3)) + SU*UV2
      QY(I3,J3)= SQY - 2.*(SBB*V(I3,J3) + SAB*U(I3,J3)) + SV*UV2
C
      D(I2,J2) = D(I3,J3)
      U(I2,J2) = U(I3,J3)
      V(I2,J2) = V(I3,J3)
      T(I2,J2) = T(I3,J3)
      P(I2,J2) = P(I3,J3)
      QX(I2,J2) =QX(I3,J3)
      QY(I2,J2) =QY(I3,J3)
      D(I1,J1) = D(I2,J2)
      U(I1,J1) = U(I2,J2)
      V(I1,J1) = V(I2,J2)
      T(I1,J1) = T(I2,J2)
      P(I1,J1) = P(I2,J2)
      QX(I1,J1) =QX(I2,J2)
      QY(I1,J1) =QY(I2,J2)
C
      IF(KMOD.EQ.3)THEN
         SSXX   = 0.
         SSYY   = 0.
         SSZZ   = 0.
         SSXY   = 0.
      DO 140 KL = 1,MXV
         K      = MOD(KL-1,NVX)+1
         L      = (KL-1)/NVX + 1
         U1     = A(K)  - U(I3,J3)
         V1     = B(L)  - V(I3,J3)
         SSXX   = SSXX  + A3(KL)*G(I3,J3,K,L)
         SSYY   = SSYY  + A4(KL)*G(I3,J3,K,L)
         SSZZ   = SSZZ  + W0(KL)*H(I3,J3,K,L)
         SSXY   = SSXY  + A6(KL)*G(I3,J3,K,L)
  140 CONTINUE
         SLXX(I3,J3) = SSXX
         SLYY(I3,J3) = SSYY
         SLZZ(I3,J3) = SSZZ
         SLXY(I3,J3) = SSXY
      SLXX(I1,J1) =SLXX(I3,J3)
      SLYY(I1,J1) =SLYY(I3,J3)
      SLZZ(I1,J1) =SLZZ(I3,J3)
      SLXY(I1,J1) =SLXY(I3,J3)
      SLXX(I2,J2) =SLXX(I3,J3)
      SLYY(I2,J2) =SLYY(I3,J3)
      SLZZ(I2,J2) =SLZZ(I3,J3)
      SLXY(I2,J2) =SLXY(I3,J3)
      ENDIF
  110 CONTINUE
      RETURN
      END
C
      SUBROUTINE BCRG(IJP,X,Y,XIX,XIY,ETAX,ETAY,QJ,D,U,V,T,
     1              G,H,RG,RH,A,B,WA,WB,CX,CY,K,L,SA1,S)
      PARAMETER (MAXBC=20)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
      COMMON /BC02 / IDA,ISA,ICA,IAS,IAE,IPA,LENA,
     1               IDB,ISB,ICB,IBS,IBE,IPB,LENB
      COMMON /AXIS / IAXIS
      DIMENSION  IJP(MXBP)
      DIMENSION    X(-2:NI3,-2:NJ3),   Y(-2:NI3,-2:NJ3)
      DIMENSION  XIX(-2:NI3,-2:NJ3), XIY(-2:NI3,-2:NJ3)
      DIMENSION ETAX(-2:NI3,-2:NJ3),ETAY(-2:NI3,-2:NJ3)
      DIMENSION   QJ(-2:NI3,-2:NJ3)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION    G(NI,NJ,NVX,NVY),   H(NI,NJ,NVX,NVY)
      DIMENSION   RG(-2:NI3,-2:NJ3),  RH(-2:NI3,-2:NJ3)
      DIMENSION  SA1(-2:NI3,-2:NJ3)
      DIMENSION    A(NVX),B(NVY),WA(NVX),WB(NVY)
      DIMENSION    S(1)
C
      DO 100 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IDA    = IDBC( 2,IBC)
      ISA    = IDBC( 3,IBC)
      ICA    = IDBC( 4,IBC)
      IAS    = IDBC( 5,IBC)
      IAE    = IDBC( 6,IBC)
      IPG    = IDBC( 7,IBC)
      IPH    = IDBC( 8,IBC)
      IPA    = IDBC(25,IBC)
      LENA   = (IAE - IAS)*ICA + 1
      IF(IBTYPE .EQ. 0 .OR. IBTYPE .EQ. 1) THEN
      DO 110 II = 1,LENA
      IA    = IJP(IPA + LENA*0 + II)
      JA    = IJP(IPA + LENA*1 + II)
      IB    = IJP(IPA + LENA*2 + II)
      JB    = IJP(IPA + LENA*3 + II)
      IC    = IJP(IPA + LENA*4 + II)
      JC    = IJP(IPA + LENA*5 + II)
      RG(IB,JB) = RG(IC,JC)
      RH(IB,JB) = RH(IC,JC)
      RG(IA,JA) = RG(IC,JC)
      RH(IA,JA) = RH(IC,JC)
  110 CONTINUE
C
      ELSEIF(IBTYPE .EQ. 2 .OR. IBTYPE .EQ. 3) THEN
      IF(IAXIS .EQ. 0) THEN
      CALL GETBC(IJP(IPA+1),RG,RH,QJ,S(IPG),S(IPH),K,L)
      ELSE
C
       DO 120 II = 1,LENA
       IA    = IJP(IPA + LENA*0 + II)
       JA    = IJP(IPA + LENA*1 + II)
       IB    = IJP(IPA + LENA*2 + II)
       JB    = IJP(IPA + LENA*3 + II)
       IC    = IJP(IPA + LENA*4 + II)
       JC    = IJP(IPA + LENA*5 + II)
       ID    = IJP(IPA + LENA*6 + II)
       JD    = IJP(IPA + LENA*7 + II)
       IE    = IJP(IPA + LENA*8 + II)
       JE    = IJP(IPA + LENA*9 + II)
C
C      RG(IA,JA) =  - RG(IE,JE)
C      RG(IB,JB) =  - RG(ID,JD)
       RG(IA,JA) =  0.
       RG(IB,JB) =  0.
       RG(IC,JC) =  0.
C      RH(IA,JA) =  - RH(IE,JE)
C      RH(IB,JB) =  - RH(ID,JD)
       RH(IA,JA) =  0.
       RH(IB,JB) =  0.
       RH(IC,JC) =  0.
  120 CONTINUE
      ENDIF
C
      ELSEIF(IBTYPE .EQ. 4) THEN
      IDB    = IDBC(12,IBC)
      ISB    = IDBC(13,IBC)
      ICB    = IDBC(14,IBC)
      IBS    = IDBC(15,IBC)
      IBE    = IDBC(16,IBC)
      IPB    = IDBC(26,IBC)
      LENB   = (IBE - IBS)*ICB + 1
      LL        = NVY - L + 1
      DO 130 II = 1,LENA
      IA1   = IJP(IPA + LENA*0 + II)
      JA1   = IJP(IPA + LENA*1 + II)
      IB1   = IJP(IPA + LENA*2 + II)
      JB1   = IJP(IPA + LENA*3 + II)
      ID1   = IJP(IPA + LENA*6 + II)
      JD1   = IJP(IPA + LENA*7 + II)
      IE1   = IJP(IPA + LENA*8 + II)
      JE1   = IJP(IPA + LENA*9 + II)
      IA2   = IJP(IPB + LENA*0 + II)
      JA2   = IJP(IPB + LENA*1 + II)
      IB2   = IJP(IPB + LENA*2 + II)
      JB2   = IJP(IPB + LENA*3 + II)
      ID2   = IJP(IPB + LENA*6 + II)
      JD2   = IJP(IPB + LENA*7 + II)
      IE2   = IJP(IPB + LENA*8 + II)
      JE2   = IJP(IPB + LENA*9 + II)
C
      RG(IA1,JA1) =  G(IE2,JE2,K,L)/QJ(IA1,JA1)
      RG(IB1,JB1) =  G(ID2,JD2,K,L)/QJ(IB1,JB1)
      RH(IA1,JA1) =  H(IE2,JE2,K,L)/QJ(IA1,JA1)
      RH(IB1,JB1) =  H(ID2,JD2,K,L)/QJ(IB1,JB1)
      RG(IA2,JA2) =  G(IE1,JE1,K,L)/QJ(IA2,JA2)
      RG(IB2,JB2) =  G(ID1,JD1,K,L)/QJ(IB2,JB2)
      RH(IA2,JA2) =  H(IE1,JE1,K,L)/QJ(IA2,JA2)
      RH(IB2,JB2) =  H(ID1,JD1,K,L)/QJ(IB2,JB2)
  130 CONTINUE
      ELSEIF(IBTYPE .EQ. 5) THEN
      DO 150 II = 1,LENA
      IA    = IJP(IPA + LENA*0 + II)
      JA    = IJP(IPA + LENA*1 + II)
      IB    = IJP(IPA + LENA*2 + II)
      JB    = IJP(IPA + LENA*3 + II)
      IC    = IJP(IPA + LENA*4 + II)
      JC    = IJP(IPA + LENA*5 + II)
      RG(IB,JB) = RG(IC,JC)
      RH(IB,JB) = RH(IC,JC)
      RG(IA,JA) = RG(IC,JC)
      RH(IA,JA) = RH(IC,JC)
  150 CONTINUE
      ENDIF
C
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DTIME(DTLOC,A,B,XIX,XIY,ETAX,ETAY,QJ)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX0 / MAXW,MAXI,MAXB,MAXG,MAXL,MAXV,MAXBP,MAXLUI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      DIMENSION  XIX(-2:NI3,-2:NJ3), XIY(-2:NI3,-2:NJ3)
      DIMENSION ETAX(-2:NI3,-2:NJ3),ETAY(-2:NI3,-2:NJ3)
      DIMENSION   QJ(-2:NI3,-2:NJ3)
      DIMENSION DTLOC(-2:NI3,-2:NJ3)
      DIMENSION    A(NVX),B(NVY)
C
      DT  = 9999.
      DTI = 9999.
      DTJ = 9999.
      DTX = 9999.
      DTY = 9999.
      MIJ = MI*MJ-2
C
      DO 300 J = 1,NJ
      DO 300 I = 1,NI
C
      DO 200 L = 1 , NVY
      DO 200 K = 1 , NVX
         AX    = XIX(I,J)
         AY    = XIY(I,J)
         CXI   = (AX * A(K) + AY * B(L))
         EI    = 1./(ABS(CXI)+1.E-20)
C
         AX    = ETAX(I,J)
         AY    = ETAY(I,J)
         CXJ   = (AX * A(K) + AY * B(L))
         EJ    = 1./(ABS(CXJ)+1.E-20)
C
         DTX = MIN(DTX,EI)
         DTY = MIN(DTY,EJ)
C
  200 CONTINUE
      DTV    = MIN(DTX,DTY)
      DT     = MIN(DTV,DT)
      DTI    = MIN(DTX,DTI)
      DTJ    = MIN(DTY,DTJ)
  300 CONTINUE
C
      WRITE(*,*)'TIME STEP FOR I  DTI =  ',DTI
      WRITE(*,*)'TIME STEP FOR J  DTJ =  ',DTJ
C
      IF(DTFIX .GT.0.) THEN
        CFL  = DTFIX/DT
        CFL1 = CFL/5.
        CFL2 = CFL
        WRITE(*,*)' USE FIX DT ,  DTFIX = ',DTFIX
        WRITE(*,*)'               CFL   = ',CFL
      ENDIF
C
      DO 400 IJ = -2,MIJ   !有問題
      DTLOC(IJ,-2) = DT
  400 CONTINUE
C
      IF(ILOCDT.EQ.1.AND.DTFIX.LE.0.)THEN
       WRITE(*,*)'USE LOCAL TIME STEP MUST SET DTFIX GT 0.'
       STOP
      ENDIF
C
      IF(ILOCDT.EQ.1.AND.DTFIX.GT.0.) THEN
      WRITE(*,*)'USE LOCAL TIME STEP'
      DTMIN = 9999.
      DTMAX = -9999.
      DO 500 IJ = -2,MIJ
      DTLOC(IJ,-2) = DTFIX/(1.+ SQRT(QJ(IJ,-2)))
      DTMIN = MIN(DTMIN,DTLOC(IJ,-2))
      DTMAX = MAX(DTMAX,DTLOC(IJ,-2))
  500 CONTINUE
      CFL = DTMIN/DT
      CFL1= CFL
      CFL2= CFL
      CFLM= DTMAX/DT
        WRITE(*,*)' USE FIX & LOC DT ,  DTFIX = ',DTFIX
        WRITE(*,*)'                     DTMIN = ',DTMIN,'   CFL=',CFL
        WRITE(*,*)'                     DTMAX = ',DTMAX,'   CFL=',CFLM
      ENDIF
      RETURN
      END

      SUBROUTINE EQULG(A,B,D,T,U,V,P,QX,QY,G,H,QJ)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX0 / MAXW,MAXI,MAXB,MAXG,MAXL,MAXV,MAXBP,MAXLUI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /CNST2/ D0,T0,P0,U0,V0,TAU0,GAMMA,PI,SPI
      COMMON /UNST0/ URGHT,VRGHT,DRGHT,TRGHT,PRGHT,
     1               ULEFT,VLEFT,DLEFT,TLEFT,PLEFT
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION   QX(-2:NI3,-2:NJ3),  QY(-2:NI3,-2:NJ3)
      DIMENSION    P(-2:NI3,-2:NJ3)
      DIMENSION   QJ(-2:NI3,-2:NJ3)
      DIMENSION    G(NI,NJ,NVX,NVY),   H(NI,NJ,NVX,NVY)
      DIMENSION    A(NVX),B(NVY)
C
C     NI4        = NI+4
C     NJ4        = NJ+4
C     NIJ        = NI4*NJ4
C     NIJK       = NI4*NJ4*NVX
C     NIJKL      = NI4*NJ4*NVX*NVY
C     DO 10 IJKL = 1,NIJKL
C     I          = MOD(IJKL-1,NI4)-1
C     J          = MOD(IJKL-1,NIJ)/NI4-1
C     K          = MOD((IJKL-1)/NIJ,NVX)+1
C     L          = (IJKL-1)/NIJK + 1
C     PP   = -1./T(I,J)*
C    1       ((A(K)-U(I,J))*(A(K)-U(I,J))+(B(L)-V(I,J))*(B(L)-V(I,J)))
C     G(I,J,K,L) = D(I,J)/(PI*T(I,J))*EXP(PP)
C     H(I,J,K,L) = 0.5*T(I,J)*G(I,J,K,L)
C  10 CONTINUE
      NIJ        = NI*NJ
      NIJK       = NI*NJ*NVX
      NIJKL      = NI*NJ*NVX*NVY
      DO 10 IJKL = 1,NIJKL
      I          = MOD(IJKL-1,NI)+1
      J          = MOD(IJKL-1,NIJ)/NI+1
      K          = MOD((IJKL-1)/NIJ,NVX)+1
      L          = (IJKL-1)/NIJK + 1
      PP   = -1./T(I,J)*
     1       ((A(K)-U(I,J))*(A(K)-U(I,J))+(B(L)-V(I,J))*(B(L)-V(I,J)))
      G(I,J,K,L) = D(I,J)/(PI*T(I,J))*EXP(PP)
      H(I,J,K,L) = 0.5*T(I,J)*G(I,J,K,L)
   10 CONTINUE
C
C
      RETURN
      END
      SUBROUTINE FDIENO2(DTLOC,CX,CY,CXI,XIX,XIY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   XIX(-1:NI2,-1:NJ2), XIY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CXI(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION   FXG(-1:MAXIJP),     FXH(-1:MAXIJP)
      DIMENSION   FMG(-1:MAXIJP),     FMH(-1:MAXIJP)
      DIMENSION   FNG(-1:MAXIJP),     FNH(-1:MAXIJP)
      DIMENSION   CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION   CXH(-1:MAXIJP), ACXH(-1:MAXIJP)
      DIMENSION   DPG(-1:MAXIJP),  DPH(-1:MAXIJP)
      DIMENSION  DPFG(-1:MAXIJP), DPFH(-1:MAXIJP)
      DIMENSION  SCXG(-1:MAXIJP), CCPG(-1:MAXIJP),  CCMG(-1:MAXIJP)
      DIMENSION  SCXH(-1:MAXIJP), CCPH(-1:MAXIJP),  CCMH(-1:MAXIJP)
      DIMENSION   DDAG(-1:MAXIJP),    DDA2G(-1:MAXIJP)
      DIMENSION   DDAH(-1:MAXIJP),    DDA2H(-1:MAXIJP)
      DIMENSION   EAG(-1:MAXIJP),     EAH(-1:MAXIJP)
      DIMENSION    EG(-1:MAXIJP),      EH(-1:MAXIJP)
      DIMENSION DPEAG(-1:MAXIJP),   DPEAH(-1:MAXIJP)
      DIMENSION DMEAG(-1:MAXIJP),   DMEAH(-1:MAXIJP)
      DIMENSION  DEBG(-1:MAXIJP),    DEBH(-1:MAXIJP)
C
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = 2.*EPS
C
      ETA = 0.5
      XI  = 0.
C
      DO 10 J = -1,NJ2
C
      DO 20 I = -1,NI2
      FXG(I)   = CXI(I,J)*RG(I,J)
      FXH(I)   = CXI(I,J)*RH(I,J)
   20 CONTINUE
C     
      DO 30 I = -1,NI1
      DPG(I)  = RG(I+1,J)  - RG(I,J)
      DPH(I)  = RH(I+1,J)  - RH(I,J)
      DPFG(I)  = FXG(I+1)   - FXG(I)
      DPFH(I)  = FXH(I+1)   - FXH(I)
   30 CONTINUE
      DPG(NI2)  = DPG(NI1)
      DPH(NI2)  = DPH(NI1)
      DPFG(NI2) = DPFG(NI1)
      DPFH(NI2) = DPFH(NI1)
C
      DO 35 I = -1,NI1
      IF(DPG(I) .NE. 0.) THEN
      CXG(I) = DPFG(I)/DPG(I)
      ELSE
      CXG(I) = 0.5*(CXI(I+1,J) + CXI(I,J))
      ENDIF
      IF(DPH(I) .NE. 0.) THEN
      CXH(I) = DPFH(I)/DPH(I)
      ELSE
      CXH(I) = 0.5*(CXI(I+1,J) + CXI(I,J))
      ENDIF
   35 CONTINUE
      CXG(NI2) = CXI(NI2,J)
      CXH(NI2) = CXI(NI2,J)
C
      DO 40 I = -1,NI2
      ACXG(I) = ABS(CXG(I))
      ACXH(I) = ABS(CXH(I))
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
      IF(CXH(I).GT.EPS)THEN
        SCXH(I) = 1.
        CCPH(I) = 1.
        CCMH(I) = 0.
      ELSEIF(CXH(I) .LT.-EPS) THEN
        SCXH(I) = -1.
        CCPH(I) = 0.
        CCMH(I) = 1.
      ELSE
        SCXH(I) = 0.
        CCPH(I) = 0.5
        CCMH(I) = 0.5
      ENDIF
      DDAG(I)    = 0.
      DDAH(I)    = 0.
      IF( LSTD .EQ. 0)THEN
      DDAG(I)    = DTLOC(I,J)*ACXG(I)
      DDAH(I)    = DTLOC(I,J)*ACXH(I)
      ENDIF
      DDA2G(I)   = DDAG(I)*DDAG(I)
      DDA2H(I)   = DDAH(I)*DDAH(I)
   40 CONTINUE
C
      DO 60 I = -1,NI2
      EAG(I)   = SCXG(I)*(1.      -    DDAG(I)     )*DPFG(I)/2.
      EAH(I)   = SCXH(I)*(1.      -    DDAH(I)     )*DPFH(I)/2.
   60 CONTINUE
C
      DO 70 I = -1,NI1
      DPEAG(I) = EAG(I+1) - EAG(I)
      DPEAH(I) = EAH(I+1) - EAH(I)
   70 CONTINUE
      DPEAG(NI2) = DPEAG(NI1)
      DPEAH(NI2) = DPEAH(NI1)
C
      DO 80 I = 0,NI2
      DMEAG(I) = DPEAG(I-1)
      DMEAH(I) = DPEAH(I-1)
   80 CONTINUE
      DMEAG(-1) = DMEAG(0)
      DMEAH(-1) = DMEAH(0)
C
      DO 150 I = -1,NI2
      IF(ABS(DPEAG(I)).LE.ABS(DMEAG(I)))THEN
       DEBG(I) = DPEAG(I)
      ELSE
       DEBG(I) = DMEAG(I)
      ENDIF
      IF(ABS(DPEAH(I)).LE.ABS(DMEAH(I)))THEN
       DEBH(I) = DPEAH(I)
      ELSE
       DEBH(I) = DMEAH(I)
      ENDIF
  150 CONTINUE
C
      DO 110 I = 0,NI2
      E1G = EAG(I  ) - ETA *DEBG(I  )
      E2G = EAG(I-1) + ETA *DEBG(I-1)
      E1H = EAH(I  ) - ETA *DEBH(I  )
      E2H = EAH(I-1) + ETA *DEBH(I-1)
      A1G = -1.
      A2G = -1.
      IF(E1G.GE.0.)THEN
      A1G =1.
      ENDIF
      IF(E2G.GE.0.)THEN
      A2G =1.
      ENDIF
      A1H = -1.
      A2H = -1.
      IF(E1H.GE.0.)THEN
      A1H =1.
      ENDIF
      IF(E2H.GE.0.)THEN
      A2H =1.
      ENDIF
      EG(I) = 0.5*(A1G+A2G)*MIN(A1G*E1G,A2G*E2G)
      EH(I) = 0.5*(A1H+A2H)*MIN(A1H*E1H,A2H*E2H)
  110 CONTINUE
      EG(-1) = 0.
      EH(-1) = 0.
C
      DO 120 I = -1,NI2
      FMG(I) = FXG(I) + EG(I)
      FMH(I) = FXH(I) + EH(I)
  120 CONTINUE
C
      DO 130 I = -1,NI1
      IF(CCPG(I).GT.0.6)THEN
       FNG(I) = FMG(I)
      ELSEIF(CCPG(I) .LT. 0.4) THEN
       FNG(I) = FMG(I+1)
      ELSE
       FNG(I) = 0.5*(FMG(I+1) + FMG(I))
      ENDIF
      IF(CCPH(I).GT.0.6)THEN
       FNH(I) = FMH(I)
      ELSEIF(CCPH(I) .LT. 0.4) THEN
       FNH(I) = FMH(I+1)
      ELSE
       FNH(I) = 0.5*(FMH(I+1) + FMH(I))
      ENDIF
  130 CONTINUE
      FNG(NI2) = FMG(NI1)
      FNH(NI2) = FMH(NI1)
C
      DO 140 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FNG(I) - FNG(I-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FNH(I) - FNH(I-1))
  140 CONTINUE
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE FDIENO3(DTLOC,CX,CY,CXI,XIX,XIY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   XIX(-1:NI2,-1:NJ2), XIY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CXI(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION   FXG(-1:MAXIJP),     FXH(-1:MAXIJP)
      DIMENSION   FMG(-1:MAXIJP),     FMH(-1:MAXIJP)
      DIMENSION   FNG(-1:MAXIJP),     FNH(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  CXH(-1:MAXIJP), ACXH(-1:MAXIJP)
      DIMENSION  DPG(-1:MAXIJP), DPH(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP), DPFH(-1:MAXIJP)
      DIMENSION SCXG(-1:MAXIJP), CCPG(-1:MAXIJP),  CCMG(-1:MAXIJP)
      DIMENSION SCXH(-1:MAXIJP), CCPH(-1:MAXIJP),  CCMH(-1:MAXIJP)
      DIMENSION   DDAG(-1:MAXIJP),    DDA2G(-1:MAXIJP)
      DIMENSION   DDAH(-1:MAXIJP),    DDA2H(-1:MAXIJP)
      DIMENSION   EAG(-1:MAXIJP),     EAH(-1:MAXIJP)
      DIMENSION    EG(-1:MAXIJP),      EH(-1:MAXIJP)
      DIMENSION  ADPG(-1:MAXIJP),    ADPH(-1:MAXIJP)
      DIMENSION   DMG(-1:MAXIJP),     DMH(-1:MAXIJP)
      DIMENSION  DMFG(-1:MAXIJP),    DMFH(-1:MAXIJP)
      DIMENSION  ADMG(-1:MAXIJP),    ADMH(-1:MAXIJP)
      DIMENSION   DAG(-1:MAXIJP),     DAH(-1:MAXIJP)
      DIMENSION   DBG(-1:MAXIJP),     DBH(-1:MAXIJP)
      DIMENSION DPDAG(-1:MAXIJP),   DPDAH(-1:MAXIJP)
      DIMENSION DPDBG(-1:MAXIJP),   DPDBH(-1:MAXIJP)
      DIMENSION DMDAG(-1:MAXIJP),   DMDAH(-1:MAXIJP)
      DIMENSION DMDBG(-1:MAXIJP),   DMDBH(-1:MAXIJP)
      DIMENSION    DG(-1:MAXIJP),      DH(-1:MAXIJP)
C
      !EPS=0.01
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = 2.*EPS
C
      ETA = 0.
      XI  = 1.
C
      DO 10 J = -1,NJ2
C
      DO 20 I = -1,NI2
      FXG(I)   = CXI(I,J)*RG(I,J)
      FXH(I)   = CXI(I,J)*RH(I,J)
   20 CONTINUE
C     
      DO 30 I  = -1,NI1
      DPG(I)   = RG(I+1,J) - RG(I,J)
      DPH(I)   = RH(I+1,J) - RH(I,J)
      DPFG(I)  = FXG(I+1)   - FXG(I)
      DPFH(I)  = FXH(I+1)   - FXH(I)
   30 CONTINUE
      DPFG(NI2) = DPFG(NI1)
      DPFH(NI2) = DPFH(NI1)
      DPG(NI2)  = DPG(NI1)
      DPH(NI2)  = DPH(NI1)
C
      DO 35 I = -1,NI1
      IF(DPG(I) .NE. 0.) THEN
      CXG(I) = DPFG(I)/DPG(I)
      ELSE
      CXG(I) = 0.5*(CXI(I+1,J) + CXI(I,J))
      ENDIF
      IF(DPH(I) .NE. 0.) THEN
      CXH(I) = DPFH(I)/DPH(I)
      ELSE
      CXH(I) = 0.5*(CXI(I+1,J) + CXI(I,J))
      ENDIF
   35 CONTINUE
      CXG(NI2) = CXI(NI2,J)
      CXH(NI2) = CXI(NI2,J)
C
      DO 40 I = -1,NI2
      ACXG(I) = ABS(CXG(I))
      ACXH(I) = ABS(CXH(I))
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
      IF(CXH(I).GT.EPS)THEN
        SCXH(I) = 1.
        CCPH(I) = 1.
        CCMH(I) = 0.
      ELSEIF(CXH(I) .LT.-EPS) THEN
        SCXH(I) = -1.
        CCPH(I) = 0.
        CCMH(I) = 1.
      ELSE
        SCXH(I) = 0.
        CCPH(I) = 0.5
        CCMH(I) = 0.5
      ENDIF
      DDAG(I)    = 0.
      DDAH(I)    = 0.
	!LSTD=1.
      IF( LSTD .EQ. 0)THEN
      DDAG(I)    = DTLOC(I,J)*ACXG(I)
      DDAH(I)    = DTLOC(I,J)*ACXH(I)
      ENDIF
      DDA2G(I)   = DDAG(I)*DDAG(I)
      DDA2H(I)   = DDAH(I)*DDAH(I)
      ADPG(I)   = ABS(DPG(I))
      ADPH(I)   = ABS(DPH(I))
   40 CONTINUE
C
      DO 50 I   = 0,NI2
      DMG (I)   = DPG (I-1)
      DMH (I)   = DPH (I-1)
      DMFG(I)   = DPFG(I-1)
      DMFH(I)   = DPFH(I-1)
      ADMG(I)   = ADPG(I-1)
      ADMH(I)   = ADPH(I-1)
   50 CONTINUE
      DMG (-1)   = DMG (0)
      DMFG(-1)   = DMFG(0)
      ADMG(-1)   = ADMG(0)
      DMH (-1)   = DMH (0)
      DMFH(-1)   = DMFH(0)
      ADMH(-1)   = ADMH(0)
C
      !E(3.2-27E,F)
      DO 60 I = -1,NI2
      EAG(I)   = SCXG(I)*(1.      -    DDAG(I)     )*DPFG(I)/2.
      EAH(I)   = SCXH(I)*(1.      -    DDAH(I)     )*DPFH(I)/2.
      DAG(I)   = SCXG(I)*(DDA2G(I) - 3.*DDAG(I) + 2.)*DPFG(I)/6.
      DBG(I)   = SCXG(I)*(DDA2G(I)              - 1.)*DPFG(I)/6.
      DAH(I)   = SCXH(I)*(DDA2H(I) - 3.*DDAH(I) + 2.)*DPFH(I)/6.
      DBH(I)   = SCXH(I)*(DDA2H(I)              - 1.)*DPFH(I)/6.
   60 CONTINUE
C
      DO 70 I = -1,NI1
      DPDAG(I) = DAG(I+1) - DAG(I)
      DPDBG(I) = DBG(I+1) - DBG(I)
      DPDAH(I) = DAH(I+1) - DAH(I)
      DPDBH(I) = DBH(I+1) - DBH(I)
   70 CONTINUE
      DPDAG(NI2) = DPDAG(NI1)
      DPDBG(NI2) = DPDBG(NI1)
      DPDAH(NI2) = DPDAH(NI1)
      DPDBH(NI2) = DPDBH(NI1)
C
      DO 80 I = 0,NI2
      DMDAG(I) = DPDAG(I-1)
      DMDBG(I) = DPDBG(I-1)
      DMDAH(I) = DPDAH(I-1)
      DMDBH(I) = DPDBH(I-1)
   80 CONTINUE
      DMDAG(-1) = DPDAG(0)
      DMDBG(-1) = DPDBG(0)
      DMDAH(-1) = DPDAH(0)
      DMDBH(-1) = DPDBH(0)
C
      DO 90 I = -1,NI2
      IF(ADPG(I) .GE. ADMG(I) ) THEN
          A1G = -1.
        IF(DMDAG(I-1).GE.0.)THEN
          A1G = 1.
        ENDIF
          A2G = -1.
        IF(DPDAG(I-1).GE.0.)THEN
          A2G = 1.
        ENDIF
        DG(I) = 0.5*(A1G+A2G)*MIN(A1G*DMDAG(I-1),A2G*DPDAG(I-1))
      ELSE
          A1G = -1.
        IF(DMDBG(I).GE.0.)THEN
          A1G = 1.
        ENDIF
          A2G = -1.
        IF(DPDBG(I).GE.0.)THEN
          A2G = 1.
        ENDIF
        DG(I) = 0.5*(A1G+A2G)*MIN(A1G*DMDBG(I),A2G*DPDBG(I))
      ENDIF
      IF(ADPH(I) .GE. ADMH(I) ) THEN
          A1H = -1.
        IF(DMDAH(I-1).GE.0.)THEN
          A1H = 1.
        ENDIF
          A2H = -1.
        IF(DPDAH(I-1).GE.0.)THEN
          A2H = 1.
        ENDIF
        DH(I) = 0.5*(A1H+A2H)*MIN(A1H*DMDAH(I-1),A2H*DPDAH(I-1))
      ELSE
          A1H = -1.
        IF(DMDBG(I).GE.0.)THEN
          A1H = 1.
        ENDIF
          A2H = -1.
        IF(DPDBH(I).GE.0.)THEN
          A2H = 1.
        ENDIF
        DH(I) = 0.5*(A1H+A2H)*MIN(A1H*DMDBH(I),A2H*DPDBH(I))
      ENDIF
   90 CONTINUE
C
      DO 110 I = 0,NI2
      E1G = EAG(I  )
      E2G = EAG(I-1)
      E1H = EAH(I  )
      E2H = EAH(I-1)
      A1G = -1.
      A2G = -1.
      IF(E1G.GE.0.)THEN
      A1G =1.
      ENDIF
      IF(E2G.GE.0.)THEN
      A2G =1.
      ENDIF
      A1H = -1.
      A2H = -1.
      IF(E1H.GE.0.)THEN
      A1H =1.
      ENDIF
      IF(E2H.GE.0.)THEN
      A2H =1.
      ENDIF
      EG(I) = 0.5*(A1G+A2G)*MIN(A1G*E1G,A2G*E2G)
      EH(I) = 0.5*(A1H+A2H)*MIN(A1H*E1H,A2H*E2H)
  110 CONTINUE
      EG(-1) = EG(0)
      EH(-1) = EH(0)
C
      DO 120 I = -1,NI2
      FMG(I) = FXG(I) + EG(I) + XI*DG(I)
      FMH(I) = FXH(I) + EH(I) + XI*DH(I)
  120 CONTINUE
C
      DO 130 I = -1,NI1
      IF(CCPG(I).GT.0.6)THEN
       FNG(I) = FMG(I)
      ELSEIF(CCPG(I) .LT. 0.4) THEN
       FNG(I) = FMG(I+1)
      ELSE
       FNG(I) = 0.5*(FMG(I+1) + FMG(I))
      ENDIF
      IF(CCPH(I).GT.0.6)THEN
       FNH(I) = FMH(I)
      ELSEIF(CCPH(I) .LT. 0.4) THEN
       FNH(I) = FMH(I+1)
      ELSE
       FNH(I) = 0.5*(FMH(I+1) + FMH(I))
      ENDIF
  130 CONTINUE
      FNG(NI2) = FMG(NI1)
      FNH(NI2) = FMH(NI1)
C
      DO 140 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FNG(I) - FNG(I-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FNH(I) - FNH(I-1))
  140 CONTINUE
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE FDITVD2(DTLOC,CX,CY,CXI,XIX,XIY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   XIX(-1:NI2,-1:NJ2), XIY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CXI(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION  FXG(-1:MAXIJP),  FXH(-1:MAXIJP)
      DIMENSION  FMG(-1:MAXIJP),  FMH(-1:MAXIJP)
      DIMENSION  FNG(-1:MAXIJP),  FNH(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  CXH(-1:MAXIJP), ACXH(-1:MAXIJP)
      DIMENSION  DPG(-1:MAXIJP), DPH(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP), DPFH(-1:MAXIJP)
      DIMENSION DDAG2(-1:MAXIJP)
      DIMENSION DDAH2(-1:MAXIJP)
      DIMENSION  EAG(-1:MAXIJP),  EAH(-1:MAXIJP)
      DIMENSION   EG(-1:MAXIJP),   EH(-1:MAXIJP)
      DIMENSION  PSG(-1:MAXIJP), PSH(-1:MAXIJP)
C
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = 2.*EPS
C
      ETA = 0.
      XI  = 0.
C
      DO 10 J = -1,NJ2
C
      DO 20 I = -1,NI2
      FXG(I)   = CXI(I,J)*RG(I,J)
      FXH(I)   = CXI(I,J)*RH(I,J)
   20 CONTINUE
C     
      DO 30 I = -1,NI1
      DPG(I)=2.*(RG(I+1,J)*QJ(I+1,J)-RG(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
      DPH(I)=2.*(RH(I+1,J)*QJ(I+1,J)-RH(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
C     DPG(I)  = RG(I+1,J)  - RG(I,J)
C     DPH(I)  = RH(I+1,J)  - RH(I,J)
      DPFG(I) = FXG(I+1)   - FXG(I)
      DPFH(I) = FXH(I+1)   - FXH(I)
   30 CONTINUE
      DPG(NI2)  = DPG(NI1)
      DPH(NI2)  = DPH(NI1)
      DPFG(NI2) = DPFG(NI1)
      DPFH(NI2) = DPFH(NI1)
C
      DO 35 I = -1,NI1
      IF(DPG(I) .NE. 0.) THEN
      CXG(I) = DPFG(I)/DPG(I)
      ELSE
      CXG(I) = 0.5*(CXI(I+1,J) + CXI(I,J))
      ENDIF
      IF(DPH(I) .NE. 0.) THEN
      CXH(I) = DPFH(I)/DPH(I)
      ELSE
      CXH(I) = 0.5*(CXI(I+1,J) + CXI(I,J))
      ENDIF
   35 CONTINUE
      CXG(NI2) = CXI(NI2,J)
      CXH(NI2) = CXI(NI2,J)
C
      DO 40 I = -1,NI2
      ACXG(I) = ABS(CXG(I))
      ACXH(I) = ABS(CXH(I))
C
      PSG(I)  = ACXG(I)
      IF(ACXG(I) .LT. EPS) PSG(I) = (CXG(I)*CXG(I)+EPSS)/EPS2
      PSH(I)  = ACXH(I)
      IF(ACXH(I) .LT. EPS) PSH(I) = (CXH(I)*CXH(I)+EPSS)/EPS2
C
      DDAG2(I)    = 0.
      DDAH2(I)    = 0.
      IF( LSTD .EQ. 0)THEN
      DDAG2(I)    = DTLOC(I,J)*ACXG(I)*ACXG(I)
      DDAH2(I)    = DTLOC(I,J)*ACXH(I)*ACXH(I)
      ENDIF
   40 CONTINUE
C
      DO 60 I = -1,NI2
      EAG(I)   = (PSG(I)      -    DDAG2(I)     )*DPG(I)/2.
      EAH(I)   = (PSH(I)      -    DDAH2(I)     )*DPH(I)/2.
   60 CONTINUE
C
      DO 110 I = 0,NI2
      E1G = EAG(I  )
      E2G = EAG(I-1)
      E1H = EAH(I  )
      E2H = EAH(I-1)
      A1G = -1.
      A2G = -1.
      IF(E1G.GE.0.)THEN
      A1G =1.
      ENDIF
      IF(E2G.GE.0.)THEN
      A2G =1.
      ENDIF
      A1H = -1.
      A2H = -1.
      IF(E1H.GE.0.)THEN
      A1H =1.
      ENDIF
      IF(E2H.GE.0.)THEN
      A2H =1.
      ENDIF
      EG(I) = 0.5*(A1G+A2G)*MIN(A1G*E1G,A2G*E2G)
      EH(I) = 0.5*(A1H+A2H)*MIN(A1H*E1H,A2H*E2H)
  110 CONTINUE
      EG(-1) = EG(0)
      EH(-1) = EH(0)
C
      DO 120 I = -1,NI2
      FMG(I) = FXG(I) + EG(I)
      FMH(I) = FXH(I) + EH(I)
  120 CONTINUE
C
      DO 130 I = -1,NI1
      FNG(I) = 0.5*(FMG(I+1)+FMG(I)
     1        -PSG(I)*DPG(I)-CXG(I)/PSG(I)*(EG(I+1)-EG(I)))
      FNH(I) = 0.5*(FMH(I+1)+FMH(I)
     1        -PSH(I)*DPH(I)-CXH(I)/PSH(I)*(EH(I+1)-EH(I)))
  130 CONTINUE
      FNG(NI2) = FMG(NI1)
      FNH(NI2) = FMH(NI1)
C
      DO 140 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FNG(I) - FNG(I-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FNH(I) - FNH(I-1))
  140 CONTINUE
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE FDJENO2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   ETAX(-1:NI2,-1:NJ2), ETAY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CET(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION  FXG(-1:MAXIJP),  FXH(-1:MAXIJP)
      DIMENSION  FMG(-1:MAXIJP),  FMH(-1:MAXIJP)
      DIMENSION  FNG(-1:MAXIJP),  FNH(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  CXH(-1:MAXIJP), ACXH(-1:MAXIJP)
      DIMENSION  DPG(-1:MAXIJP), DPH(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP), DPFH(-1:MAXIJP)
      DIMENSION SCXG(-1:MAXIJP), CCPG(-1:MAXIJP),  CCMG(-1:MAXIJP)
      DIMENSION SCXH(-1:MAXIJP), CCPH(-1:MAXIJP),  CCMH(-1:MAXIJP)
      DIMENSION  DDAG(-1:MAXIJP), DDA2G(-1:MAXIJP)
      DIMENSION  DDAH(-1:MAXIJP), DDA2H(-1:MAXIJP)
      DIMENSION  EAG(-1:MAXIJP),  EAH(-1:MAXIJP)
      DIMENSION   EG(-1:MAXIJP),   EH(-1:MAXIJP)
      DIMENSION DPEAG(-1:MAXIJP),   DPEAH(-1:MAXIJP)
      DIMENSION DMEAG(-1:MAXIJP),   DMEAH(-1:MAXIJP)
      DIMENSION  DEBG(-1:MAXIJP),    DEBH(-1:MAXIJP)
C
      EPS = ENTEPS
      ETA = 0.5
      XI  = 0.
C
      DO 10 I = -1,NI2
C
      DO 20 J = -1,NJ2
      FXG(J)   = CET(I,J)*RG(I,J)
      FXH(J)   = CET(I,J)*RH(I,J)
   20 CONTINUE
C     
      DO 30 J = -1,NJ1
      DPG(J)  = RG(I,J+1)  - RG(I,J)
      DPH(J)  = RH(I,J+1)  - RH(I,J)
      DPFG(J)  = FXG(J+1)   - FXG(J)
      DPFH(J)  = FXH(J+1)   - FXH(J)
   30 CONTINUE
      DPG(NJ2) = DPG(NJ1)
      DPH(NJ2) = DPH(NJ1)
      DPFG(NJ2) = DPFG(NJ1)
      DPFH(NJ2) = DPFH(NJ1)
C
      DO 35 J = -1,NJ1
      IF(DPG(J) .NE. 0.) THEN
      CXG(J) = DPFG(J)/DPG(J)
      ELSE
      CXG(J) = 0.5*(CET(I,J+1) + CET(I,J))
      ENDIF
      IF(DPH(J) .NE. 0.) THEN
      CXH(J) = DPFH(J)/DPH(J)
      ELSE
      CXH(J) = 0.5*(CET(I,J+1) + CET(I,J))
      ENDIF
   35 CONTINUE
      CXG(NJ2) = CET(I,NJ2)
      CXH(NJ2) = CET(I,NJ2)
C
      DO 40 J = -1,NJ2
      ACXG(J) = ABS(CXG(J))
      ACXH(J) = ABS(CXH(J))
      IF(CXG(J).GT.EPS)THEN
        SCXG(J) = 1.
        CCPG(J) = 1.
        CCMG(J) = 0.
      ELSEIF(CXG(J).LT.-EPS)THEN
        SCXG(J) = -1.
        CCPG(J) = 0.
        CCMG(J) = 1.
      ELSE
        SCXG(J) = 0.
        CCPG(J) = 0.5
        CCMG(J) = 0.5
      ENDIF
      IF(CXH(J).GT.EPS)THEN
        SCXH(J) = 1.
        CCPH(J) = 1.
        CCMH(J) = 0.
      ELSEIF(CXH(J).LT.-EPS)THEN
        SCXH(J) = -1.
        CCPH(J) = 0.
        CCMH(J) = 1.
      ELSE
        SCXH(J) = 0.
        CCPH(J) = 0.5
        CCMH(J) = 0.5
      ENDIF
      DDAG(J)    = 0.
      DDAH(J)    = 0.
      IF(LSTD .EQ. 0) THEN
      DDAG(J)    = DTLOC(I,J)*ACXG(J)
      DDAH(J)    = DTLOC(I,J)*ACXH(J)
      ENDIF
      DDA2G(J)   = DDAG(J)*DDAG(J)
      DDA2H(J)   = DDAH(J)*DDAH(J)
   40 CONTINUE
C
      DO 60 J = -1,NJ2
      EAG(J)   = SCXG(J)*(1.      -    DDAG(J)     )*DPFG(J)/2.
      EAH(J)   = SCXH(J)*(1.      -    DDAH(J)     )*DPFH(J)/2.
   60 CONTINUE
C
      DO 70 J = -1,NJ1
      DPEAG(J) = EAG(J+1) - EAG(J)
      DPEAH(J) = EAH(J+1) - EAH(J)
   70 CONTINUE
      DPEAG(NJ2) = DPEAG(NJ1)
      DPEAH(NJ2) = DPEAH(NJ1)
C
      DO 80 J = 0,NJ2
      DMEAG(J) = DPEAG(J-1)
      DMEAH(J) = DPEAH(J-1)
   80 CONTINUE
      DMEAG(-1) = DMEAG(0)
      DMEAH(-1) = DMEAH(0)
C
      DO 150 J = -1,NJ2
      IF(ABS(DPEAG(J)).LE.ABS(DMEAG(J)))THEN
       DEBG(J) = DPEAG(J)
      ELSE
       DEBG(J) = DMEAG(J)
      ENDIF
      IF(ABS(DPEAH(J)).LE.ABS(DMEAH(J)))THEN
       DEBH(J) = DPEAH(J)
      ELSE
       DEBH(J) = DMEAH(J)
      ENDIF
  150 CONTINUE
C
      DO 110 J = 0,NJ2
      E1G = EAG(J  ) - ETA *DEBG(J  )
      E2G = EAG(J-1) + ETA *DEBG(J-1)
      E1H = EAH(J  ) - ETA *DEBH(J  )
      E2H = EAH(J-1) + ETA *DEBH(J-1)
      A1G = -1.
      A2G = -1.
      IF(E1G.GE.0.)THEN
      A1G =1.
      ENDIF
      IF(E2G.GE.0.)THEN
      A2G =1.
      ENDIF
      A1H = -1.
      A2H = -1.
      IF(E1H.GE.0.)THEN
      A1H =1.
      ENDIF
      IF(E2H.GE.0.)THEN
      A2H =1.
      ENDIF
      EG(J) = 0.5*(A1G+A2G)*MIN(A1G*E1G,A2G*E2G)
      EH(J) = 0.5*(A1H+A2H)*MIN(A1H*E1H,A2H*E2H)
  110 CONTINUE
      EG(-1) = 0.
      EH(-1) = 0.
C
      DO 120 J = -1,NJ2
      FMG(J) = FXG(J) + EG(J)
      FMH(J) = FXH(J) + EH(J)
  120 CONTINUE
C
      DO 130 J = -1,NJ1
      IF(CCPG(J).GT.0.6)THEN
       FNG(J) = FMG(J)
      ELSEIF(CCPG(J) .LT. 0.4) THEN
       FNG(J) = FMG(J+1)
      ELSE
       FNG(J) = 0.5*(FMG(J+1) + FMG(J))
      ENDIF
      IF(CCPH(J).GT.0.6)THEN
       FNH(J) = FMH(J)
      ELSEIF(CCPH(J) .LT. 0.4) THEN
       FNH(J) = FMH(J+1)
      ELSE
       FNH(J) = 0.5*(FMH(J+1) + FMH(J))
      ENDIF
  130 CONTINUE
      FNG(NJ2) = FMG(NJ1)
      FNH(NJ2) = FMH(NJ1)
C
      DO 140 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FNG(J) - FNG(J-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FNH(J) - FNH(J-1))
  140 CONTINUE
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE FDJENO3(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   ETAX(-1:NI2,-1:NJ2), ETAY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CET(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION   FXG(-1:MAXIJP),     FXH(-1:MAXIJP)
      DIMENSION   FMG(-1:MAXIJP),     FMH(-1:MAXIJP)
      DIMENSION   FNG(-1:MAXIJP),     FNH(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  CXH(-1:MAXIJP), ACXH(-1:MAXIJP)
      DIMENSION  DPG(-1:MAXIJP), DPH(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP), DPFH(-1:MAXIJP)
      DIMENSION SCXG(-1:MAXIJP), CCPG(-1:MAXIJP),  CCMG(-1:MAXIJP)
      DIMENSION SCXH(-1:MAXIJP), CCPH(-1:MAXIJP),  CCMH(-1:MAXIJP)
      DIMENSION   DDAG(-1:MAXIJP),    DDA2G(-1:MAXIJP)
      DIMENSION   DDAH(-1:MAXIJP),    DDA2H(-1:MAXIJP)
      DIMENSION   EAG(-1:MAXIJP),     EAH(-1:MAXIJP)
      DIMENSION    EG(-1:MAXIJP),      EH(-1:MAXIJP)
      DIMENSION  ADPG(-1:MAXIJP),    ADPH(-1:MAXIJP)
      DIMENSION   DMG(-1:MAXIJP),     DMH(-1:MAXIJP)
      DIMENSION  DMFG(-1:MAXIJP),    DMFH(-1:MAXIJP)
      DIMENSION  ADMG(-1:MAXIJP),    ADMH(-1:MAXIJP)
      DIMENSION   DAG(-1:MAXIJP),     DAH(-1:MAXIJP)
      DIMENSION   DBG(-1:MAXIJP),     DBH(-1:MAXIJP)
      DIMENSION DPDAG(-1:MAXIJP),   DPDAH(-1:MAXIJP)
      DIMENSION DPDBG(-1:MAXIJP),   DPDBH(-1:MAXIJP)
      DIMENSION DMDAG(-1:MAXIJP),   DMDAH(-1:MAXIJP)
      DIMENSION DMDBG(-1:MAXIJP),   DMDBH(-1:MAXIJP)
      DIMENSION    DG(-1:MAXIJP),      DH(-1:MAXIJP)
C
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = 2.*EPS
C
      ETA = 0.
      XI  = 1.
C
      DO 10 I = -1,NI2
C
      DO 20 J = -1,NJ2
      FXG(J)   = CET(I,J)*RG(I,J)
      FXH(J)   = CET(I,J)*RH(I,J)
   20 CONTINUE
C     
      DO 30 J  = -1,NJ1
      DPG(J)   = RG(I,J+1) - RG(I,J)
      DPH(J)   = RH(I,J+1) - RH(I,J)
      DPFG(J)  = FXG(J+1)   - FXG(J)
      DPFH(J)  = FXH(J+1)   - FXH(J)
   30 CONTINUE
      DPFG(NJ2) = DPFG(NJ1)
      DPFH(NJ2) = DPFH(NJ1)
      DPG(NJ2)  = DPG(NJ1)
      DPH(NJ2)  = DPH(NJ1)
C
      DO 35 J = -1,NJ1
      IF(DPG(J) .NE. 0.) THEN
      CXG(J) = DPFG(J)/DPG(J)
      ELSE
      CXG(J) = 0.5*(CET(I,J+1) + CET(I,J))
      ENDIF
      IF(DPH(J) .NE. 0.) THEN
      CXH(J) = DPFH(J)/DPH(J)
      ELSE
      CXH(J) = 0.5*(CET(I,J+1) + CET(I,J))
      ENDIF
   35 CONTINUE
      CXG(NJ2) = CET(I,NJ2)
      CXH(NJ2) = CET(I,NJ2)
C
      DO 40 J = -1,NJ2
      ACXG(J) = ABS(CXG(J))
      ACXH(J) = ABS(CXH(J))
      IF(CXG(J).GT.EPS)THEN
        SCXG(J) = 1.
        CCPG(J) = 1.
        CCMG(J) = 0.
      ELSEIF(CXG(J).LT.-EPS)THEN
        SCXG(J) = -1.
        CCPG(J) = 0.
        CCMG(J) = 1.
      ELSE
        SCXG(J) = 0.
        CCPG(J) = 0.5
        CCMG(J) = 0.5
      ENDIF
      IF(CXH(J).GT.EPS)THEN
        SCXH(J) = 1.
        CCPH(J) = 1.
        CCMH(J) = 0.
      ELSEIF(CXH(J).LT.-EPS)THEN
        SCXH(J) = -1.
        CCPH(J) = 0.
        CCMH(J) = 1.
      ELSE
        SCXH(J) = 0.
        CCPH(J) = 0.5
        CCMH(J) = 0.5
      ENDIF
      DDAG(J)    = 0.
      DDAH(J)    = 0.
      IF(LSTD .EQ. 0) THEN
      DDAG(J)    = DTLOC(I,J)*ACXG(J)
      DDAH(J)    = DTLOC(I,J)*ACXH(J)
      ENDIF
      DDA2G(J)   = DDAG(J)*DDAG(J)
      DDA2H(J)   = DDAH(J)*DDAH(J)
      ADPG(J)   = ABS(DPG(J))
      ADPH(J)   = ABS(DPH(J))
   40 CONTINUE
C
      DO 50 J   = 0,NJ2
      DMG (J)   = DPG (J-1)
      DMH (J)   = DPH (J-1)
      DMFG(J)   = DPFG(J-1)
      DMFH(J)   = DPFH(J-1)
      ADMG(J)   = ADPG(J-1)
      ADMH(J)   = ADPH(J-1)
   50 CONTINUE
      DMG (-1)   = DMG (0)
      DMFG(-1)   = DMFG(0)
      ADMG(-1)   = ADMG(0)
      DMH (-1)   = DMH (0)
      DMFH(-1)   = DMFH(0)
      ADMH(-1)   = ADMH(0)
C
      DO 60 J = -1,NJ2
      EAG(J)   = SCXG(J)*(1.      -    DDAG(J)     )*DPFG(J)/2.
      EAH(J)   = SCXH(J)*(1.      -    DDAH(J)     )*DPFH(J)/2.
      DAG(J)   = SCXG(J)*(DDA2G(J) - 3.*DDAG(J) + 2.)*DPFG(J)/6.
      DBG(J)   = SCXG(J)*(DDA2G(J)              - 1.)*DPFG(J)/6.
      DAH(J)   = SCXH(J)*(DDA2H(J) - 3.*DDAH(J) + 2.)*DPFH(J)/6.
      DBH(J)   = SCXH(J)*(DDA2H(J)              - 1.)*DPFH(J)/6.
   60 CONTINUE
C
      DO 70 J = -1,NJ1
      DPDAG(J) = DAG(J+1) - DAG(J)
      DPDBG(J) = DBG(J+1) - DBG(J)
      DPDAH(J) = DAH(J+1) - DAH(J)
      DPDBH(J) = DBH(J+1) - DBH(J)
   70 CONTINUE
      DPDAG(NJ2) = DPDAG(NJ1)
      DPDBG(NJ2) = DPDBG(NJ1)
      DPDAH(NJ2) = DPDAH(NJ1)
      DPDBH(NJ2) = DPDBH(NJ1)
C
      DO 80 J = 0,NJ2
      DMDAG(J) = DPDAG(J-1)
      DMDBG(J) = DPDBG(J-1)
      DMDAH(J) = DPDAH(J-1)
      DMDBH(J) = DPDBH(J-1)
   80 CONTINUE
      DMDAG(-1) = DPDAG(0)
      DMDBG(-1) = DPDBG(0)
      DMDAH(-1) = DPDAH(0)
      DMDBH(-1) = DPDBH(0)
C
      DO 90 J = -1,NJ2
      IF(ADPG(J) .GE. ADMG(J) ) THEN
          A1G = -1.
        IF(DMDAG(J-1).GE.0.)THEN
          A1G = 1.
        ENDIF
          A2G = -1.
        IF(DPDAG(J-1).GE.0.)THEN
          A2G = 1.
        ENDIF
        DG(J) = 0.5*(A1G+A2G)*MIN(A1G*DMDAG(J-1),A2G*DPDAG(J-1))
      ELSE
          A1G = -1.
        IF(DMDBG(J).GE.0.)THEN
          A1G = 1.
        ENDIF
          A2G = -1.
        IF(DPDBG(J).GE.0.)THEN
          A2G = 1.
        ENDIF
        DG(J) = 0.5*(A1G+A2G)*MIN(A1G*DMDBG(J),A2G*DPDBG(J))
      ENDIF
      IF(ADPH(J) .GE. ADMH(J) ) THEN
          A1H = -1.
        IF(DMDAH(J-1).GE.0.)THEN
          A1H = 1.
        ENDIF
          A2H = -1.
        IF(DPDAH(J-1).GE.0.)THEN
          A2H = 1.
        ENDIF
        DH(J) = 0.5*(A1H+A2H)*MIN(A1H*DMDAH(J-1),A2H*DPDAH(J-1))
      ELSE
          A1H = -1.
        IF(DMDBG(J).GE.0.)THEN
          A1H = 1.
        ENDIF
          A2H = -1.
        IF(DPDBH(J).GE.0.)THEN
          A2H = 1.
        ENDIF
        DH(J) = 0.5*(A1H+A2H)*MIN(A1H*DMDBH(J),A2H*DPDBH(J))
      ENDIF
   90 CONTINUE
C
      DO 110 J = 0,NJ2
      E1G = EAG(J  )
      E2G = EAG(J-1)
      E1H = EAH(J  )
      E2H = EAH(J-1)
      A1G = -1.
      A2G = -1.
      IF(E1G.GE.0.)THEN
      A1G =1.
      ENDIF
      IF(E2G.GE.0.)THEN
      A2G =1.
      ENDIF
      A1H = -1.
      A2H = -1.
      IF(E1H.GE.0.)THEN
      A1H =1.
      ENDIF
      IF(E2H.GE.0.)THEN
      A2H =1.
      ENDIF
      EG(J) = 0.5*(A1G+A2G)*MIN(A1G*E1G,A2G*E2G)
      EH(J) = 0.5*(A1H+A2H)*MIN(A1H*E1H,A2H*E2H)
  110 CONTINUE
      EG(-1) = EG(0)
      EH(-1) = EH(0)
C
      DO 120 J = -1,NJ2
      FMG(J) = FXG(J) + EG(J) + XI*DG(J)
      FMH(J) = FXH(J) + EH(J) + XI*DH(J)
  120 CONTINUE
C
      DO 130 J = -1,NJ1
      IF(CCPG(J).GT.0.6)THEN
       FNG(J) = FMG(J)
      ELSEIF(CCPG(J) .LT. 0.4) THEN
       FNG(J) = FMG(J+1)
      ELSE
       FNG(J) = 0.5*(FMG(J+1) + FMG(J))
      ENDIF
      IF(CCPH(J).GT.0.6)THEN
       FNH(J) = FMH(J)
      ELSEIF(CCPH(J) .LT. 0.4) THEN
       FNH(J) = FMH(J+1)
      ELSE
       FNH(J) = 0.5*(FMH(J+1) + FMH(J))
      ENDIF
  130 CONTINUE
      FNG(NJ2) = FMG(NJ1)
      FNH(NJ2) = FMH(NJ1)
C
      DO 140 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FNG(J) - FNG(J-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FNH(J) - FNH(J-1))
  140 CONTINUE
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE FDJTVD2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   ETAX(-1:NI2,-1:NJ2), ETAY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CET(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION  FXG(-1:MAXIJP),  FXH(-1:MAXIJP)
      DIMENSION  FMG(-1:MAXIJP),  FMH(-1:MAXIJP)
      DIMENSION  FNG(-1:MAXIJP),  FNH(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  CXH(-1:MAXIJP), ACXH(-1:MAXIJP)
      DIMENSION  DPG(-1:MAXIJP), DPH(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP), DPFH(-1:MAXIJP)
      DIMENSION DDAG2(-1:MAXIJP)
      DIMENSION DDAH2(-1:MAXIJP)
      DIMENSION  EAG(-1:MAXIJP),  EAH(-1:MAXIJP)
      DIMENSION   EG(-1:MAXIJP),   EH(-1:MAXIJP)
      DIMENSION  PSG(-1:MAXIJP), PSH(-1:MAXIJP)
C
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = 2.*EPS
C
      ETA = 0.
      XI  = 0.
C
      DO 10 I = -1,NI2
C
      DO 20 J = -1,NJ2
      FXG(J)   = CET(I,J)*RG(I,J)
      FXH(J)   = CET(I,J)*RH(I,J)
   20 CONTINUE
C     
      DO 30 J = -1,NJ1
      DPG(J) = 2.*(RG(I,J+1)*QJ(I,J+1)-RG(I,J)*QJ(I,J))
     1                    /(QJ(I,J+1)+QJ(I,J))
      DPH(J) = 2.*(RH(I,J+1)*QJ(I,J+1)-RH(I,J)*QJ(I,J))
     1                    /(QJ(I,J+1)+QJ(I,J))
C     DPG(J)  = RG(I,J+1)  - RG(I,J)
C     DPH(J)  = RH(I,J+1)  - RH(I,J)
      DPFG(J) = FXG(J+1)   - FXG(J)
      DPFH(J) = FXH(J+1)   - FXH(J)
   30 CONTINUE
      DPG(NJ2) = DPG(NJ1)
      DPH(NJ2) = DPH(NJ1)
      DPFG(NJ2) = DPFG(NJ1)
      DPFH(NJ2) = DPFH(NJ1)
C
      DO 35 J = -1,NJ1
      IF(DPG(J) .NE. 0.) THEN
      CXG(J) = DPFG(J)/DPG(J)
      ELSE
      CXG(J) = 0.5*(CET(I,J+1) + CET(I,J))
      ENDIF
      IF(DPH(J) .NE. 0.) THEN
      CXH(J) = DPFH(J)/DPH(J)
      ELSE
      CXH(J) = 0.5*(CET(I,J+1) + CET(I,J))
      ENDIF
   35 CONTINUE
      CXG(NJ2) = CET(I,NJ2)
      CXH(NJ2) = CET(I,NJ2)
C
      DO 40 J = -1,NJ2
      ACXG(J) = ABS(CXG(J))
      ACXH(J) = ABS(CXH(J))
C
      PSG(J)  = ACXG(J)
      IF(ACXG(J) .LT. EPS) PSG(J) = (CXG(J)*CXG(J)+EPSS)/EPS2
      PSH(J)  = ACXH(J)
      IF(ACXH(J) .LT. EPS) PSH(J) = (CXH(J)*CXH(J)+EPSS)/EPS2
C
      DDAG2(J)    = 0.
      DDAH2(J)    = 0.
      IF(LSTD .EQ. 0) THEN
      DDAG2(J)    = DTLOC(I,J)*ACXG(J)*ACXG(J)
      DDAH2(J)    = DTLOC(I,J)*ACXH(J)*ACXH(J)
      ENDIF
   40 CONTINUE
C
      DO 60 J = -1,NJ2
      EAG(J)   = (PSG(J)      -    DDAG2(J)     )*DPG(J)/2.
      EAH(J)   = (PSH(J)      -    DDAH2(J)     )*DPH(J)/2.
   60 CONTINUE
C
      DO 110 J = 0,NJ2
      E1G = EAG(J  )
      E2G = EAG(J-1)
      E1H = EAH(J  )
      E2H = EAH(J-1)
      A1G = -1.
      A2G = -1.
      IF(E1G.GE.0.)THEN
      A1G =1.
      ENDIF
      IF(E2G.GE.0.)THEN
      A2G =1.
      ENDIF
      A1H = -1.
      A2H = -1.
      IF(E1H.GE.0.)THEN
      A1H =1.
      ENDIF
      IF(E2H.GE.0.)THEN
      A2H =1.
      ENDIF
      EG(J) = 0.5*(A1G+A2G)*MIN(A1G*E1G,A2G*E2G)
      EH(J) = 0.5*(A1H+A2H)*MIN(A1H*E1H,A2H*E2H)
  110 CONTINUE
      EG(-1) = 0.
      EH(-1) = 0.
C
      DO 120 J = -1,NJ2
      FMG(J) = FXG(J) + EG(J)
      FMH(J) = FXH(J) + EH(J)
  120 CONTINUE
C
      DO 130 J = -1,NJ1
      FNG(J) = 0.5*(FMG(J+1)+FMG(J)
     1        -PSG(J)*DPG(J)-CXG(J)/PSG(J)*(EG(J+1)-EG(J)))
      FNH(J) = 0.5*(FMH(J+1)+FMH(J)
     1        -PSH(J)*DPH(J)-CXH(J)/PSH(J)*(EH(J+1)-EH(J)))
  130 CONTINUE
      FNG(NJ2) = FMG(NJ1)
      FNH(NJ2) = FMH(NJ1)
C
      DO 140 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FNG(J) - FNG(J-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FNH(J) - FNH(J-1))
  140 CONTINUE
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE FEITVD2(DTLOC,CX,CY,CXI,XIX,XIY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
C
C     Original FDS formular
C
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   XIX(-1:NI2,-1:NJ2), XIY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CXI(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION  FXG(-1:MAXIJP),  FXH(-1:MAXIJP)
      DIMENSION  FMG(-1:MAXIJP),  FMH(-1:MAXIJP)
      DIMENSION  FNG(-1:MAXIJP),  FNH(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  CXH(-1:MAXIJP), ACXH(-1:MAXIJP)
      DIMENSION  DPG(-1:MAXIJP), DPH(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP), DPFH(-1:MAXIJP)
      DIMENSION SCXG(-1:MAXIJP), CCPG(-1:MAXIJP),  CCMG(-1:MAXIJP)
      DIMENSION SCXH(-1:MAXIJP), CCPH(-1:MAXIJP),  CCMH(-1:MAXIJP)
      DIMENSION DDAG(-1:MAXIJP)
      DIMENSION DDAH(-1:MAXIJP)
      DIMENSION  EAG(-1:MAXIJP),  EAH(-1:MAXIJP)
      DIMENSION   EG(-1:MAXIJP),   EH(-1:MAXIJP)
C
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = 2.*EPS
C
      ETA = 0.
      XI  = 0.
C
      DO 10 J = -1,NJ2
C
      DO 20 I = -1,NI2
      FXG(I)   = CXI(I,J)*RG(I,J)
      FXH(I)   = CXI(I,J)*RH(I,J)
   20 CONTINUE
C     
      DO 30 I = -1,NI1
      DPG(I)=2.*(RG(I+1,J)*QJ(I+1,J)-RG(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
      DPH(I)=2.*(RH(I+1,J)*QJ(I+1,J)-RH(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
C     DPG(I)  = RG(I+1,J)  - RG(I,J)
C     DPH(I)  = RH(I+1,J)  - RH(I,J)
      DPFG(I) = FXG(I+1)   - FXG(I)
      DPFH(I) = FXH(I+1)   - FXH(I)
   30 CONTINUE
      DPG(NI2)  = DPG(NI1)
      DPH(NI2)  = DPH(NI1)
      DPFG(NI2) = DPFG(NI1)
      DPFH(NI2) = DPFH(NI1)
C
      DO 35 I = -1,NI1
      IF(DPG(I) .NE. 0.) THEN
      CXG(I) = DPFG(I)/DPG(I)
      ELSE
      CXG(I) = 0.5*(CXI(I+1,J) + CXI(I,J))
      ENDIF
      IF(DPH(I) .NE. 0.) THEN
      CXH(I) = DPFH(I)/DPH(I)
      ELSE
      CXH(I) = 0.5*(CXI(I+1,J) + CXI(I,J))
      ENDIF
   35 CONTINUE
      CXG(NI2) = CXI(NI2,J)
      CXH(NI2) = CXI(NI2,J)
C
      DO 40 I = -1,NI2
      ACXG(I) = ABS(CXG(I))
      ACXH(I) = ABS(CXH(I))
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
      IF(CXH(I).GT.EPS)THEN
        SCXH(I) = 1.
        CCPH(I) = 1.
        CCMH(I) = 0.
      ELSEIF(CXH(I) .LT.-EPS) THEN
        SCXH(I) = -1.
        CCPH(I) = 0.
        CCMH(I) = 1.
      ELSE
        SCXH(I) = 0.
        CCPH(I) = 0.5
        CCMH(I) = 0.5
      ENDIF
      DDAG(I)    = 0.
      DDAH(I)    = 0.
      IF( LSTD .EQ. 0)THEN
      DDAG(I)    = DTLOC(I,J)*ACXG(I)
      DDAH(I)    = DTLOC(I,J)*ACXH(I)
      ENDIF
   40 CONTINUE
C
      DO 60 I = -1,NI2
      EAG(I)   = SCXG(I)*(1.      -    DDAG(I)     )*DPFG(I)/2.
      EAH(I)   = SCXH(I)*(1.      -    DDAH(I)     )*DPFH(I)/2.
   60 CONTINUE
C
      DO 110 I = 0,NI2
      E1G = EAG(I  )
      E2G = EAG(I-1)
      E1H = EAH(I  )
      E2H = EAH(I-1)
      A1G = -1.
      A2G = -1.
      IF(E1G.GE.0.)THEN
      A1G =1.
      ENDIF
      IF(E2G.GE.0.)THEN
      A2G =1.
      ENDIF
      A1H = -1.
      A2H = -1.
      IF(E1H.GE.0.)THEN
      A1H =1.
      ENDIF
      IF(E2H.GE.0.)THEN
      A2H =1.
      ENDIF
      EG(I) = 0.5*(A1G+A2G)*MIN(A1G*E1G,A2G*E2G)
      EH(I) = 0.5*(A1H+A2H)*MIN(A1H*E1H,A2H*E2H)
  110 CONTINUE
      EG(-1) = EG(0)
      EH(-1) = EH(0)
C
      DO 120 I = -1,NI2
      FMG(I) = FXG(I) + EG(I)
      FMH(I) = FXH(I) + EH(I)
  120 CONTINUE
C
      DO 130 I = -1,NI1
      FNG(I) = FMG(I+1) - CCPG(I)*(FMG(I+1) - FMG(I))
      FNH(I) = FMH(I+1) - CCPH(I)*(FMH(I+1) - FMH(I))
C     IF(CCPG(I).GT.0.6)THEN
C      FNG(I) = FMG(I)
C     ELSEIF(CCPG(I) .LT. 0.4) THEN
C      FNG(I) = FMG(I+1)
C     ELSE
C      FNG(I) = 0.5*(FMG(I+1) + FMG(I))
C     ENDIF
C     IF(CCPH(I).GT.0.6)THEN
C      FNH(I) = FMH(I)
C     ELSEIF(CCPH(I) .LT. 0.4) THEN
C      FNH(I) = FMH(I+1)
C     ELSE
C      FNH(I) = 0.5*(FMH(I+1) + FMH(I))
C     ENDIF
  130 CONTINUE
      FNG(NI2) = FMG(NI1)
      FNH(NI2) = FMH(NI1)
C
      DO 140 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FNG(I) - FNG(I-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FNH(I) - FNH(I-1))
  140 CONTINUE
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE FEJTVD2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   ETAX(-1:NI2,-1:NJ2), ETAY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CET(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION  FXG(-1:MAXIJP),  FXH(-1:MAXIJP)
      DIMENSION  FMG(-1:MAXIJP),  FMH(-1:MAXIJP)
      DIMENSION  FNG(-1:MAXIJP),  FNH(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  CXH(-1:MAXIJP), ACXH(-1:MAXIJP)
      DIMENSION  DPG(-1:MAXIJP), DPH(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP), DPFH(-1:MAXIJP)
      DIMENSION SCXG(-1:MAXIJP), CCPG(-1:MAXIJP),  CCMG(-1:MAXIJP)
      DIMENSION SCXH(-1:MAXIJP), CCPH(-1:MAXIJP),  CCMH(-1:MAXIJP)
      DIMENSION DDAG(-1:MAXIJP)
      DIMENSION DDAH(-1:MAXIJP)
      DIMENSION  EAG(-1:MAXIJP),  EAH(-1:MAXIJP)
      DIMENSION   EG(-1:MAXIJP),   EH(-1:MAXIJP)
C
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = 2.*EPS
C
      ETA = 0.
      XI  = 0.
C
      DO 10 I = -1,NI2
C
      DO 20 J = -1,NJ2
      FXG(J)   = CET(I,J)*RG(I,J)
      FXH(J)   = CET(I,J)*RH(I,J)
   20 CONTINUE
C     
      DO 30 J = -1,NJ1
      DPG(J) = 2.*(RG(I,J+1)*QJ(I,J+1)-RG(I,J)*QJ(I,J))
     1                    /(QJ(I,J+1)+QJ(I,J))
      DPH(J) = 2.*(RH(I,J+1)*QJ(I,J+1)-RH(I,J)*QJ(I,J))
     1                    /(QJ(I,J+1)+QJ(I,J))
C     DPG(J)  = RG(I,J+1)  - RG(I,J)
C     DPH(J)  = RH(I,J+1)  - RH(I,J)
      DPFG(J) = FXG(J+1)   - FXG(J)
      DPFH(J) = FXH(J+1)   - FXH(J)
   30 CONTINUE
      DPG(NJ2) = DPG(NJ1)
      DPH(NJ2) = DPH(NJ1)
      DPFG(NJ2) = DPFG(NJ1)
      DPFH(NJ2) = DPFH(NJ1)
C
      DO 35 J = -1,NJ1
      IF(DPG(J) .NE. 0.) THEN
      CXG(J) = DPFG(J)/DPG(J)
      ELSE
      CXG(J) = 0.5*(CET(I,J+1) + CET(I,J))
      ENDIF
      IF(DPH(J) .NE. 0.) THEN
      CXH(J) = DPFH(J)/DPH(J)
      ELSE
      CXH(J) = 0.5*(CET(I,J+1) + CET(I,J))
      ENDIF
   35 CONTINUE
      CXG(NJ2) = CET(I,NJ2)
      CXH(NJ2) = CET(I,NJ2)
C
      DO 40 J = -1,NJ2
      ACXG(J) = ABS(CXG(J))
      ACXH(J) = ABS(CXH(J))
      IF(CXG(J).GT.EPS)THEN
        SCXG(J) = 1.
        CCPG(J) = 1.
        CCMG(J) = 0.
      ELSEIF(CXG(J).LT.-EPS)THEN
        SCXG(J) = -1.
        CCPG(J) = 0.
        CCMG(J) = 1.
      ELSE
        SCXG(J) = 0.
        CCPG(J) = 0.5
        CCMG(J) = 0.5
      ENDIF
      IF(CXH(J).GT.EPS)THEN
        SCXH(J) = 1.
        CCPH(J) = 1.
        CCMH(J) = 0.
      ELSEIF(CXH(J).LT.-EPS)THEN
        SCXH(J) = -1.
        CCPH(J) = 0.
        CCMH(J) = 1.
      ELSE
        SCXH(J) = 0.
        CCPH(J) = 0.5
        CCMH(J) = 0.5
      ENDIF
      DDAG(J)    = 0.
      DDAH(J)    = 0.
      IF(LSTD .EQ. 0) THEN
      DDAG(J)    = DTLOC(I,J)*ACXG(J)
      DDAH(J)    = DTLOC(I,J)*ACXH(J)
      ENDIF
   40 CONTINUE
C
      DO 60 J = -1,NJ2
      EAG(J)   = SCXG(J)*(1.      -    DDAG(J)     )*DPFG(J)/2.
      EAH(J)   = SCXH(J)*(1.      -    DDAH(J)     )*DPFH(J)/2.
   60 CONTINUE
C
      DO 110 J = 0,NJ2
      E1G = EAG(J  )
      E2G = EAG(J-1)
      E1H = EAH(J  )
      E2H = EAH(J-1)
      A1G = -1.
      A2G = -1.
      IF(E1G.GE.0.)THEN
      A1G =1.
      ENDIF
      IF(E2G.GE.0.)THEN
      A2G =1.
      ENDIF
      A1H = -1.
      A2H = -1.
      IF(E1H.GE.0.)THEN
      A1H =1.
      ENDIF
      IF(E2H.GE.0.)THEN
      A2H =1.
      ENDIF
      EG(J) = 0.5*(A1G+A2G)*MIN(A1G*E1G,A2G*E2G)
      EH(J) = 0.5*(A1H+A2H)*MIN(A1H*E1H,A2H*E2H)
  110 CONTINUE
      EG(-1) = 0.
      EH(-1) = 0.
C
      DO 120 J = -1,NJ2
      FMG(J) = FXG(J) + EG(J)
      FMH(J) = FXH(J) + EH(J)
  120 CONTINUE
C
      DO 130 J = -1,NJ1
      FNG(J)   = FMG(J+1) - CCPG(J)*(FMG(J+1)-FMG(J))
      FNH(J)   = FMH(J+1) - CCPH(J)*(FMH(J+1)-FMH(J))
C     IF(CCPG(J).GT.0.6)THEN
C      FNG(J) = FMG(J)
C     ELSEIF(CCPG(J) .LT. 0.4) THEN
C      FNG(J) = FMG(J+1)
C     ELSE
C      FNG(J) = 0.5*(FMG(J+1) + FMG(J))
C     ENDIF
C     IF(CCPH(J).GT.0.6)THEN
C      FNH(J) = FMH(J)
C     ELSEIF(CCPH(J) .LT. 0.4) THEN
C      FNH(J) = FMH(J+1)
C     ELSE
C      FNH(J) = 0.5*(FMH(J+1) + FMH(J))
C     ENDIF
  130 CONTINUE
      FNG(NJ2) = FMG(NJ1)
      FNH(NJ2) = FMH(NJ1)
C
      DO 140 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FNG(J) - FNG(J-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FNH(J) - FNH(J-1))
  140 CONTINUE
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE FHITVD2(DTLOC,CX,CY,CXI,XIX,XIY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
C
      DIMENSION   XIX(-1:NI2,-1:NJ2), XIY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CXI(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION  GI(-1:MAXIJP),GJ(-1:MAXIJP),DG(-1:MAXIJP)
      DIMENSION  AS(-1:MAXIJP),BS(-1:MAXIJP)
      DIMENSION  HI(-1:MAXIJP),HJ(-1:MAXIJP),DH(-1:MAXIJP)
      DIMENSION  AH(-1:MAXIJP),HS(-1:MAXIJP)
C
      EPS = 0.20
C
      DO 10 J = -1,NJ2
      DO 20 I = -1,NI1
      DG(I)=2.*(RG(I+1,J)*QJ(I+1,J)-RG(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
      DH(I)=2.*(RH(I+1,J)*QJ(I+1,J)-RH(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
      IF(LSTD.EQ.0)THEN
      GJ(I) = 0.5*(ABS(CXI(I,J))
     1             - DTLOC(I,J)*CXI(I,J)*CXI(I,J))*DG(I)
      HJ(I) = 0.5*(ABS(CXI(I,J))
     1             - DTLOC(I,J)*CXI(I,J)*CXI(I,J))*DH(I)
      ELSE
      GJ(I) = 0.5*ABS(CXI(I,J))*DG(I)
      HJ(I) = 0.5*ABS(CXI(I,J))*DH(I)
      ENDIF
      BS(I) = ABS(GJ(I))
      HS(I) = ABS(HJ(I))
   20 CONTINUE
      DG(NI2) = DG(NI1)
      GJ(NI2) = GJ(NI1)
      BS(NI2) = BS(NI1)
      DH(NI2) = DH(NI1)
      HJ(NI2) = HJ(NI1)
      HS(NI2) = HS(NI1)
C
      DO 30 I = -1,NI2
      IF(GJ(I).GE.0.)THEN
      AS(I) = 1.
      ELSE
      AS(I) =-1.
      ENDIF
      IF(HJ(I).GE.0.)THEN
      AH(I) = 1.
      ELSE
      AH(I) =-1.
      ENDIF
   30 CONTINUE
      DO 32 I = 0,NI2
      GI(I) = 0.5*(AS(I-1)+AS(I))*MIN(BS(I-1),BS(I))
      HI(I) = 0.5*(AH(I-1)+AH(I))*MIN(HS(I-1),HS(I))
   32 CONTINUE
      GI( -1) = GI(  0)
      HI( -1) = HI(  0)
C
      DO 40 I = -1,NI1
      IF(DG(I).EQ.0.)THEN
       GAM  = 0.
      ELSE
       GAM  = (GI(I+1) - GI(I))/DG(I)
      ENDIF
      IF(DH(I).EQ.0.)THEN
       GAH  = 0.
      ELSE
       GAH  = (HI(I+1) - HI(I))/DH(I)
      ENDIF
      CC0   = CXI(I,J) + GAM
      CCH   = CXI(I,J) + GAH
      PSI   = ABS(CC0)
      PSH   = ABS(CCH)
      IF(PSI.LT.EPS)THEN
      PSI   = 0.5*(CC0*CC0 + EPS*EPS)/EPS
      ENDIF
      IF(PSH.LT.EPS)THEN
      PSH   = 0.5*(CCH*CCH + EPS*EPS)/EPS
      ENDIF
      FG(I)  = 0.5*(CXI(I,J)*RG(I,J)
     &                 +CXI(I+1,J)*RG(I+1,J)
     &                 +GI(I)+GI(I+1)-PSI*DG(I))
      FH(I)  = 0.5*(CXI(I,J)*RH(I,J)
     &                 +CXI(I+1,J)*RH(I+1,J)
     &                 +HI(I)+HI(I+1)-PSH*DH(I))
   40 CONTINUE
      DO 50 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FG(I) - FG(I-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FH(I) - FH(I-1))
   50 CONTINUE
   10 CONTINUE
      RETURN
      END

      SUBROUTINE FHJTVD2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,
     1                            FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
C
      DIMENSION  ETAX(-1:NI2,-1:NJ2),ETAY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2),  RH(-1:NI2,-1:NJ2)
      DIMENSION    CET(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION  GI(-1:MAXIJP),GJ(-1:MAXIJP),DG(-1:MAXIJP)
      DIMENSION  AS(-1:MAXIJP),BS(-1:MAXIJP)
      DIMENSION  HI(-1:MAXIJP),HJ(-1:MAXIJP),DH(-1:MAXIJP)
      DIMENSION  AH(-1:MAXIJP),HS(-1:MAXIJP)
C
      EPS = 0.2
C
      DO 10 I = -1,NI2
      DO 20 J = -1,NJ1
      DG(J) = 2.*(RG(I,J+1)*QJ(I,J+1)-RG(I,J)*QJ(I,J))
     1                    /(QJ(I,J+1)+QJ(I,J))
      DH(J) = 2.*(RH(I,J+1)*QJ(I,J+1)-RH(I,J)*QJ(I,J))
     1                    /(QJ(I,J+1)+QJ(I,J))
C     DG(J) = RG(I,J+1) - RG(I,J)
C     DH(J) = RH(I,J+1) - RH(I,J)
      IF(LSTD.EQ.0)THEN
      GJ(J) = 0.5*(ABS(CET(I,J))
     1             - DTLOC(I,J)*CET(I,J)*CET(I,J))*DG(J)
      HJ(J) = 0.5*(ABS(CET(I,J))
     1             - DTLOC(I,J)*CET(I,J)*CET(I,J))*DH(J)
      ELSE
      GJ(J) = 0.5*ABS(CET(I,J)) *DG(J)
      HJ(J) = 0.5*ABS(CET(I,J)) *DH(J)
      ENDIF
      BS(J) = ABS(GJ(J))
      HS(J) = ABS(HJ(J))
   20 CONTINUE
      DG(NJ2) = DG(NJ1)
      GJ(NJ2) = GJ(NJ1)
      BS(NJ2) = BS(NJ1)
      DH(NJ2) = DH(NJ1)
      HJ(NJ2) = HJ(NJ1)
      HS(NJ2) = HS(NJ1)
C
      DO 30 J = -1,NJ2
      IF(GJ(J).GT.0.)THEN
      AS(J) = 1.
      ELSE
      AS(J) = -1.
      ENDIF
      IF(HJ(J).GT.0.)THEN
      AH(J) = 1.
      ELSE
      AH(J) = -1.
      ENDIF
   30 CONTINUE
      DO 32 J = 0,NJ2
      GI(J) = 0.5*(AS(J-1)+AS(J))*MIN(BS(J-1),BS(J))
      HI(J) = 0.5*(AH(J-1)+AH(J))*MIN(HS(J-1),HS(J))
   32 CONTINUE
      GI(-1) = GI(  0)
      HI(-1) = HI(  0)
C
      DO 40 J = -1,NJ1
      IF(DG(J).NE.0.)THEN
      GAM  = (GI(J+1) - GI(J))/DG(J)
      ELSE
      GAM  = 0.
      ENDIF
      IF(DH(J).NE.0.)THEN
      GAH  = (HI(J+1) - HI(J))/DH(J)
      ELSE
      GAH  = 0.
      ENDIF
      CC0   = CET(I,J) + GAM
      CCH   = CET(I,J) + GAH
      PSI   = ABS(CC0)
      PSH   = ABS(CCH)
      IF(PSI.LT.EPS)THEN
      PSI   = 0.5*(CC0*CC0 + EPS*EPS)/EPS
      ENDIF
      IF(PSH.LT.EPS)THEN
      PSH   = 0.5*(CCH*CCH + EPS*EPS)/EPS
      ENDIF
      FG(J)  = 0.5*(CET(I,J)*RG(I,J)
     &         +CET(I,J+1)*RG(I,J+1)
     &         +GI(J)+GI(J+1)-PSI*DG(J))
      FH(J)  = 0.5*(CET(I,J)*RH(I,J)
     &         +CET(I,J+1)*RH(I,J+1)
     &         +HI(J)+HI(J+1)-PSH*DH(J))
   40 CONTINUE
      DO 50 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FG(J) - FG(J-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FH(J) - FH(J-1))
   50 CONTINUE
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FORCE(X,Y,D,U,V,T,P,G,IJP)
      PARAMETER (NV=8000)
      PARAMETER (MAXBC=20)
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /TAPE/ INAME,IGRID,IWRIT,IREAD,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /VELC / W0(NV),W1(NV),W2(NV),W3(NV),
     1               A3(NV),A4(NV),A5(NV),A6(NV)
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
      COMMON /CNST1/ TINF,XMPINF,XMUINF,XNUINF,DNINF,SSINF,UDINF,
     1               RHOINF,PINF,EINF,CINF,REINF
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
C
      DIMENSION    X(-2:NI3,-2:NJ3),   Y(-2:NI3,-2:NJ3)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION    P(-2:NI3,-2:NJ3)
      DIMENSION    G(NI,NJ,NVX,NVY)
      DIMENSION    TXX(501),TYY(501),TXY(501),CP(501)
      DIMENSION    IJP(MXBP)
C
      QINF = 0.5*RHOINF*CINF*CINF
      RINF = 0.5*RHOINF*UDINF*UDINF
      SAER = 1.
      DO 500 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IDA    = IDBC( 2,IBC)
      ISA    = IDBC( 3,IBC)
      ICA    = IDBC( 4,IBC)
      IAS    = IDBC( 5,IBC)
      IAE    = IDBC( 6,IBC)
      IPA    = IDBC(25,IBC)
      LENA   = (IAE - IAS)*ICA + 1
C
      IF(IBTYPE .EQ. 1) THEN
      PPL = 0.
      PPD = 0.
      SSL = 0.
      SSD = 0.
      DO 100 II = 1 , LENA
      I    =  IJP(IPA+LENA*4+II)
      J    =  IJP(IPA+LENA*5+II)
C
      SXX = 0.
      SXY = 0.
      SYY = 0.
C
      DN   = D(I,J)
      UW   = U(I,J)
      VW   = V(I,J)
      TW   = T(I,J)
      PW   = DN*TW
C
      DO 110 KL = 1,MXV
       K      = MOD(KL-1,NVX)+1
       L      = (KL-1)/NVX + 1
      SXX = SXX + G(I,1,K,L)*A3(KL)
      SYY = SYY + G(I,1,K,L)*A4(KL)
      SXY = SXY + G(I,1,K,L)*A6(KL)
  110 CONTINUE
      TXX(II) = 2.*(SXX - DN*UW*UW - 0.5*PW)
      TYY(II) = 2.*(SYY - DN*VW*VW - 0.5*PW)
      TXY(II) = 2.*(SXY - DN*UW*VW)
       CP(II) = (PW - 1.)*QINF/RINF
  100 CONTINUE
      DO 200 II = 1 , LENA-1
      I1   =  IJP(IPA+LENA*4+II)
      J1   =  IJP(IPA+LENA*5+II)
      I2   =  IJP(IPA+LENA*4+II+1)
      J2   =  IJP(IPA+LENA*5+II+1)
      PW     = 0.5*(P(I1,J1) + P(I2,J2))
      STXX   = 0.5*(TXX(II+1) + TXX(II))
      STYY   = 0.5*(TYY(II+1) + TYY(II))
      STXY   = 0.5*(TXY(II+1) + TXY(II))
      EX   =   Y(I2,J2) - Y(I1,J1)
      EY   = -(X(I2,J2) - X(I1,J1))
      PPL  = PPL  + (PW - 1.)*EY
      PPD  = PPD  + (PW - 1.)*EX
      SSL  = SSL  + STYY*EY + STXY*EX
      SSD  = SSD  + STXX*EX + STXY*EY
  200 CONTINUE
C
      FLP  = PPL*QINF
      FDP  = PPD*QINF
      FLS  = SSL*QINF
      FDS  = SSD*QINF
C
      CNP  = FLP/RINF/SAER
      CAP  = FDP/RINF/SAER
      CNS  = FLS/RINF/SAER
      CAS  = FDS/RINF/SAER
C
      CN   = CNP + CNS
      CA   = CAP + CAS
C
      WRITE(IFORC,900)ITER,TIME,CN,CNP,CNS,CA,CAP,CAS
      WRITE(*,920)ITER,TIME,CN,CNP,CNS,CA,CAP,CAS
C
      REWIND (17)
      WRITE(17,*)'NBC=',IBC,'LEN=',LENA,'IDA=',IDA,'ISA=',ISA
      WRITE(17,*)'X,Y,D,U,V,T,P,TXX,TYY,TXY,CP'
      DO 300 II = 1,LENA
      I    =  IJP(IPA+LENA*4+II)
      J    =  IJP(IPA+LENA*5+II)
      WRITE(17,910)X(I,J),Y(I,J),D(I,J),U(I,J),V(I,J),
     1             T(I,J),P(I,J),TXX(II),TYY(II),TXY(II),CP(II)
  300 CONTINUE
      ENDIF
  500 CONTINUE
C     CLOSE(17)
  900 FORMAT(1X,I5,7E13.4)
  910 FORMAT(3(1X,F8.4),8(1X,E12.3))
  920 FORMAT('   ITER   = ',I5,5X,
     1       '   TIME   = ',E12.3,/,
     2       '   CN     = ',E12.3,5X,
     2       '   CNP    = ',E12.3,5X,
     2       '   CNS    = ',E12.3,/,
     2       '   CA     = ',E12.3,5X,
     2       '   CAP    = ',E12.3,5X,
     2       '   CAS    = ',E12.3,/)
C
      RETURN
      END

      SUBROUTINE FVIENO2(DTLOC,CX,CY,CXI,XIX,XIY,QJ,
     1                               FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
C
      DIMENSION   XIX(-1:NI2,-1:NJ2), XIY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CXI(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION   GI(-1:MAXIJP), GJ(-1:MAXIJP),  DGP(-1:MAXIJP)
      DIMENSION   AS(-1:MAXIJP), BS(-1:MAXIJP)
      DIMENSION  GAM(-1:MAXIJP),GJP(-1:MAXIJP), GJM(-1:MAXIJP)
      DIMENSION  GJJ(-1:MAXIJP),DDG(-1:MAXIJP),XMB1(-1:MAXIJP)
      DIMENSION   AP(-1:MAXIJP), AM(-1:MAXIJP)
C
      DIMENSION   HI(-1:MAXIJP), HJ(-1:MAXIJP),  DHP(-1:MAXIJP)
      DIMENSION   AH(-1:MAXIJP), BH(-1:MAXIJP)
      DIMENSION  GAH(-1:MAXIJP),HJP(-1:MAXIJP), HJM(-1:MAXIJP)
      DIMENSION  HJJ(-1:MAXIJP),DDH(-1:MAXIJP),XMH1(-1:MAXIJP)
      DIMENSION   HP(-1:MAXIJP), HM(-1:MAXIJP)
C
      EPS = 0.2
C
      DO 10 J  = -1,NJ2
C
      DO 20 I  = -1,NI1
      DGP(I)   = RG(I+1,J) - RG(I,J)
      DHP(I)   = RH(I+1,J) - RH(I,J)
   20 CONTINUE
      DGP(NI2) = DGP (NI1)
      DHP(NI2) = DHP (NI1)
C
      DO 30 I  = -1,NI1
      DDG(I)   = ABS(DGP(I+1) - DGP(I)   )
      DDH(I)   = ABS(DHP(I+1) - DHP(I)   )
   30 CONTINUE
      DDG(NI2) = DDG(NI1)
      DDH(NI2) = DDH(NI1)
C
      DO 32 I  = 0,NI2
      IF(DDG(I-1) .LE. DDG(I))THEN
      XMB1(I)  = DDG(I-1)
      ELSE
      XMB1(I)  = DDG(I)
      ENDIF
      IF(DDH(I-1) .LE. DDH(I))THEN
      XMH1(I)  = DDH(I-1)
      ELSE
      XMH1(I)  = DDH(I)
      ENDIF
   32 CONTINUE
C
C      I = -1
C
      DDG2  = DDG(-1)
      DDG4  = 0.
      IF(DDG4.LE.DDG2)THEN
        XMB1(-1) = DDG4
      ELSE
        XMB1(-1) = DDG2
      ENDIF
      DDH2  = DDH(-1)
      DDH4  = 0.
      IF(DDH4.LE.DDH2)THEN
        XMH1(-1) = DDH4
      ELSE
        XMH1(-1) = DDH2
      ENDIF
C
      DO 34 I = 0,NI2
      GJP(I)  = DGP (I)   - 0.5 * XMB1(I)
      GJM(I)  = DGP (I-1) + 0.5 * XMB1(I-1)
      HJP(I)  = DHP (I)   - 0.5 * XMH1(I)
      HJM(I)  = DHP (I-1) + 0.5 * XMH1(I-1)
   34 CONTINUE
C
      GJP(-1)  = GJP (0)
      GJM(-1)  = GJM (0)
      HJP(-1)  = HJP (0)
      HJM(-1)  = HJM (0)
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
      IF(HJP(I).GT.0.)THEN
        HP(I)  = 1.
      ELSE
        HP(I)  = -1.
      ENDIF
      IF(HJM(I).GT.0.)THEN
        HM(I)  = 1.
      ELSE
        HM(I)  = -1.
      ENDIF
   40 CONTINUE
      DO 42 I  = -1,NI2
      GJJ(I)= 0.5*(AP(I)+AM(I))*MIN(AP(I)*GJP(I),AM(I)*GJM(I))
      HJJ(I)= 0.5*(HP(I)+HM(I))*MIN(HP(I)*HJP(I),HM(I)*HJM(I))
   42 CONTINUE
C
      IF(LSTD.EQ.0)THEN
      DO 50 I = -1,NI2
      GI(I)   = 0.5 * (ABS(CXI(I,J))
     1               -DTLOC(I,J)*CXI(I,J)*CXI(I,J))*GJJ(I)
      HI(I)   = 0.5 * (ABS(CXI(I,J))
     1               -DTLOC(I,J)*CXI(I,J)*CXI(I,J))*HJJ(I)
   50 CONTINUE
      ELSE
      DO 52 I = -1,NI2
      GI(I)   = 0.5 * ABS(CXI(I,J))*GJJ(I)
      HI(I)   = 0.5 * ABS(CXI(I,J))*HJJ(I)
   52 CONTINUE
      ENDIF
C
      DO 60 I = -1,NI1
      IF(DGP(I).EQ.0.)THEN
        GAM(I) = 0.
      ELSE
        GAM(I) = (GI(I+1) - GI(I))/DGP(I)
      ENDIF
      IF(DHP(I).EQ.0.)THEN
        GAH(I) = 0.
      ELSE
        GAH(I) = (HI(I+1) - HI(I))/DHP(I)
      ENDIF
   60 CONTINUE
      DO 62 I = -1,NI1
      CC0     = CXI(I,J) + GAM(I)
      PSI     = ABS(CC0)
      CCH     = CXI(I,J) + GAH(I)
      PSH     = ABS(CCH)
      IF(PSI.LT.EPS)THEN
      PSI     = 0.5*(CC0*CC0 + EPS*EPS)/EPS
      ENDIF
      IF(PSH.LT.EPS)THEN
      PSH     = 0.5*(CCH*CCH + EPS*EPS)/EPS
      ENDIF
      FG(I)   = 0.5*(CXI(I,J)*RG(I,J)
     &            +CXI(I+1,J)*RG(I+1,J)
     &            +GI(I)+GI(I+1)-PSI*DGP(I))
      FH(I)   = 0.5*(CXI(I,J)*RH(I,J)
     &            +CXI(I+1,J)*RH(I+1,J)
     &            +HI(I)+HI(I+1)-PSI*DHP(I))
   62 CONTINUE
C
      DO 70 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FG(I) - FG(I-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FH(I) - FH(I-1))
   70 CONTINUE
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE FVITVD2(DTLOC,CX,CY,CXI,XIX,XIY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   XIX(-1:NI2,-1:NJ2), XIY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CXI(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION  GI(-1:MAXIJP),DG(-1:MAXIJP)
      DIMENSION  HI(-1:MAXIJP),DH(-1:MAXIJP)
      DIMENSION  ADG(-1:MAXIJP),ADH(-1:MAXIJP)
      DIMENSION  AS(-1:MAXIJP),AH(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP), DPFH(-1:MAXIJP)
      DIMENSION  FXG(-1:MAXIJP),  FXH(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  CXH(-1:MAXIJP), ACXH(-1:MAXIJP)
      DIMENSION  PSG(-1:MAXIJP), PSH(-1:MAXIJP)
      DIMENSION  SIG(-1:MAXIJP), SIH(-1:MAXIJP)
      DIMENSION  GAG(-1:MAXIJP), GAH(-1:MAXIJP)
      DIMENSION  AGAG(-1:MAXIJP), AGAH(-1:MAXIJP)
C
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = EPS*2.
C
      DO 10 J = -1,NJ2
C
      DO 20 I = -1,NI2
      FXG(I)   = CXI(I,J)*RG(I,J)
      FXH(I)   = CXI(I,J)*RH(I,J)
   20 CONTINUE
C     
      DO 30 I = -1,NI1
      DG(I)=2.*(RG(I+1,J)*QJ(I+1,J)-RG(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
      DH(I)=2.*(RH(I+1,J)*QJ(I+1,J)-RH(I,J)*QJ(I,J))
     1           /(QJ(I+1,J)+QJ(I,J))
C     DH(I) = RH(I+1,J) - RH(I,J)
      DPFG(I) = FXG(I+1)   - FXG(I)
      DPFH(I) = FXH(I+1)   - FXH(I)
   30 CONTINUE
      DG(NI2)  = DG(NI1)
      DH(NI2)  = DH(NI1)
      DPFG(NI2) = DPFG(NI1)
      DPFH(NI2) = DPFH(NI1)
C
      DO 40 I = -1,NI1
      IF(DG(I) .NE. 0.) THEN
      CXG(I) = DPFG(I)/DG(I)
      ELSE
      CXG(I) = 0.5*(CXI(I+1,J) + CXI(I,J))
      ENDIF
      IF(DH(I) .NE. 0.) THEN
      CXH(I) = DPFH(I)/DH(I)
      ELSE
      CXH(I) = 0.5*(CXI(I+1,J) + CXI(I,J))
      ENDIF
   40 CONTINUE
      CXG(NI2) = CXG(NI1)
      CXH(NI2) = CXG(NI1)
C
      DO 50 I = -1,NI2
      IF(DG(I).GE.0.)THEN
      AS(I) = 1.
      ELSE
      AS(I) =-1.
      ENDIF
      IF(DH(I).GE.0.)THEN
      AH(I) = 1.
      ELSE
      AH(I) =-1.
      ENDIF
      ADG(I) = ABS(DG(I))
      ADH(I) = ABS(DH(I))
      ACXG(I) = ABS(CXG(I))
      ACXH(I) = ABS(CXH(I))
      PSG(I)  = ACXG(I)
      PSH(I)  = ACXH(I)
      IF(ACXG(I) .LT. EPS) PSG(I) = (CXG(I)*CXG(I)+EPSS)/EPS2
      IF(ACXH(I) .LT. EPS) PSH(I) = (CXH(I)*CXH(I)+EPSS)/EPS2
   50 CONTINUE
      DO 60 I = 0,NI2
      GI(I) = 0.5*(AS(I-1)+AS(I))*MIN(ADG(I-1),ADG(I))
      HI(I) = 0.5*(AH(I-1)+AH(I))*MIN(ADH(I-1),ADH(I))
   60 CONTINUE
      GI( -1) = GI(  0)
      HI( -1) = HI(  0)
C
      DO 70 I = -1,NI2
      IF(LSTD.EQ.0)THEN
      SIG(I) = 0.5*(PSG(I)
     1             - DTLOC(I,J)*CXG(I)*CXG(I))*DG(I)
      SIH(I) = 0.5*(PSH(I)
     1             - DTLOC(I,J)*CXH(I)*CXH(I))*DH(I)
      ELSE
      SIG(I) = 0.5*PSG(I)*DG(I)
      SIH(I) = 0.5*PSH(I)*DH(I)
      ENDIF
   70 CONTINUE
C
      DO 80 I = -1,NI1
      IF(DG(I).NE.0.)THEN
      GAG(I)  = SIG(I)*(GI(I+1)-GI(I))/DG(I)
      ELSE
      GAG(I)  = 0.
      ENDIF
      IF(DH(I).NE.0.)THEN
      GAH(I)  = SIH(I)*(HI(I+1)-HI(I))/DH(I)
      ELSE
      GAH(I)  = 0.
      ENDIF
   80 CONTINUE
      GAG(NI2) = GAG(NI1)
      GAH(NI2) = GAH(NI1)
C
      DO 90 I = -1,NI2
      AGAG(I)  = ABS(CXG(I) + GAG(I))
      AGAH(I)  = ABS(CXH(I) + GAH(I))
      PSG(I)   = AGAG(I)
      PSH(I)   = AGAH(I)
      IF(AGAG(I) .LT. EPS) PSG(I) = (AGAG(I)*AGAG(I)+EPSS)/EPS2
      IF(AGAH(I) .LT. EPS) PSH(I) = (AGAH(I)*AGAH(I)+EPSS)/EPS2
   90 CONTINUE
C
      DO 100 I = -1,NI1
      FG(I)  = 0.5*(FXG(I+1)+FXG(I)+SIG(I)*(GI(I+1)+GI(I))
     1              -PSG(I)*DG(I))
      FH(I)  = 0.5*(FXH(I+1)+FXH(I)+SIH(I)*(HI(I+1)+HI(I))
     1              -PSH(I)*DH(I))
  100 CONTINUE
      DO 110 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FG(I) - FG(I-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FH(I) - FH(I-1))
  110 CONTINUE
   10 CONTINUE
      RETURN
      END

      SUBROUTINE FVIUPW(DTLOC,CX,CY,CXI,XIX,XIY,QJ,
     1                         FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
C
      DIMENSION   XIX(-1:NI2,-1:NJ2), XIY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2),  RH(-1:NI2,-1:NJ2)
      DIMENSION   CXI(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
      DO 20 J = -1,NJ2
      DO 10 I = -1,NI1
      DG1 = (RG(I+1,J) + RG(I,J))
      DG2 = (RG(I+1,J) - RG(I,J))
      DH1 = (RH(I+1,J) + RH(I,J))
      DH2 = (RH(I+1,J) - RH(I,J))
      FG(I)  = 0.5*(CXI(I,J)*DG1 - ABS(CXI(I,J))*DG2)
      FH(I)  = 0.5*(CXI(I,J)*DH1 - ABS(CXI(I,J))*DH2)
   10 CONTINUE
C
C     QJ2 = 2./(QJ(I+1,J)+QJ(I,J))
C     DG2 = (RG(I+1,J)*QJ(I+1,J) - RG(I,J)*QJ(I,J))*QJ2
C     DH2 = (RH(I+1,J)*QJ(I+1,J) - RH(I,J)*QJ(I,J))*QJ2
C     DFG = (RG(I+1,J)*QJ(I+1,J)*CXI(I+1,J)
C    1      -RG(I,J)*QJ(I,J)*CXI(I,J))*QJ2
C     DFH = (RH(I+1,J)*QJ(I+1,J)*CXI(I+1,J)
C    1      -RH(I,J)*QJ(I,J)*CXI(I,J))*QJ2
C     IF(DG2 .NE. 0.)THEN
C     CXG = DFG/DG2
C     ELSE
C     CXG = 0.5*(CXI(I+1,J)+CXI(I,J))
C     ENDIF
C     IF(DH2 .NE. 0.)THEN
C     CXH = DFH/DH2
C     ELSE
C     CXH = 0.5*(CXI(I+1,J)+CXI(I,J))
C     ENDIF
C     IF(CXG .GE. 0.)THEN
C     FG(I)  = CXI(I,J)*RG(I,J)*QJ(I,J)*QJ2
C     ELSE
C     FG(I)  = CXI(I,J)*RG(I+1,J)*QJ(I+1,J)*QJ2
C     ENDIF
C     IF(CXH .GE. 0.)THEN
C     FG(I)  = CXI(I,J)*RH(I,J)*QJ(I,J)*QJ2
C     ELSE
C     FH(I)  = CXI(I,J)*RH(I+1,J)*QJ(I+1,J)*QJ2
C     ENDIF
C  10 CONTINUE
C
C     DO 20 J = -1,NJ2
C     DO 10 I = -1,NI1
C
C     QJ2 = 2./(QJ(I+1,J)+QJ(I,J))
C     DG1 = (RG(I+1,J)*QJ(I+1,J) + RG(I,J)*QJ(I,J))*QJ2
C     DG2 = (RG(I+1,J)*QJ(I+1,J) - RG(I,J)*QJ(I,J))*QJ2
C     DH1 = (RH(I+1,J)*QJ(I+1,J) + RH(I,J)*QJ(I,J))*QJ2
C     DH2 = (RH(I+1,J)*QJ(I+1,J) - RH(I,J)*QJ(I,J))*QJ2
C
C     DFG = (RG(I+1,J)*QJ(I+1,J)*CXI(I+1,J)
C    1      -RG(I,J)*QJ(I,J)*CXI(I,J))*QJ2
C     DFH = (RH(I+1,J)*QJ(I+1,J)*CXI(I+1,J)
C    1      -RH(I,J)*QJ(I,J)*CXI(I,J))*QJ2
C     IF(DG2 .NE. 0.)THEN
C     CXG = DFG/DG2
C     ELSE
C     CXG = 0.5*(CXI(I+1,J)+CXI(I,J))
C     ENDIF
C     IF(DH2 .NE. 0.)THEN
C     CXH = DFH/DH2
C     ELSE
C     CXH = 0.5*(CXI(I+1,J)+CXI(I,J))
C     ENDIF
C
C     FG(I)  = 0.5*(CXI(I,J)*DG1 - ABS(CXI(I,J))*DG2)
C     FH(I)  = 0.5*(CXI(I,J)*DH1 - ABS(CXI(I,J))*DH2)
C     FG(I)  = 0.5*(CXG*DG1 - ABS(CXG)*DG2)
C     FH(I)  = 0.5*(CXH*DH1 - ABS(CXH)*DH2)
C  10 CONTINUE
C
      DO 30 I = 0,NI1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FG(I) - FG(I-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FH(I) - FH(I-1))
   30 CONTINUE
   20 CONTINUE
C
      RETURN
      END
      SUBROUTINE FVJENO2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,
     1                               FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
C
      DIMENSION  ETAX(-1:NI2,-1:NJ2),ETAY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2),  RH(-1:NI2,-1:NJ2)
      DIMENSION    CET(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
C
      DIMENSION   GI(-1:MAXIJP), GJ(-1:MAXIJP),  DGP(-1:MAXIJP)
      DIMENSION   AS(-1:MAXIJP), BS(-1:MAXIJP)
      DIMENSION  GJP(-1:MAXIJP), GJM(-1:MAXIJP)
      DIMENSION  GJJ(-1:MAXIJP),DDG(-1:MAXIJP),XMB1(-1:MAXIJP)
C
      DIMENSION   HI(-1:MAXIJP), HJ(-1:MAXIJP),  DHP(-1:MAXIJP)
      DIMENSION   AH(-1:MAXIJP), BH(-1:MAXIJP)
      DIMENSION  HJP(-1:MAXIJP), HJM(-1:MAXIJP)
      DIMENSION  HJJ(-1:MAXIJP),DDH(-1:MAXIJP),XMH1(-1:MAXIJP)
C
      EPS      = 0.2
C
      DO 10 I  = -1,NI2
C
      DO 20 J  = -1,NJ1
      DGP (J)  = RG(I,J+1) - RG(I,J)
      DHP (J)  = RH(I,J+1) - RH(I,J)
   20 CONTINUE
      DGP(NJ2) = DGP (NJ1)
      DHP(NJ2) = DHP (NJ1)
C
      DO 30 J  = -1,NJ1
      DDG(J)   = ABS(DGP(J+1) - DGP(J)   )
      DDH(J)   = ABS(DHP(J+1) - DHP(J)   )
   30 CONTINUE
      DDG(NJ2) = DDG(NJ1)
      DDH(NJ2) = DDH(NJ1)
C
      DO 32 J  = 0,NJ2
      IF(DDG(J-1).LE.DDG(J))THEN
      XMB1(J)  = DDG(J-1)
      ELSE
      XMB1(J)  = DDG(J)
      ENDIF
      IF(DDH(J-1).LE.DDH(J))THEN
      XMH1(J)  = DDH(J-1)
      ELSE
      XMH1(J)  = DDH(J)
      ENDIF
   32 CONTINUE
C
C      J = -1
C
      DDG2  = DDG(-1)
      DDG4  = 0.
      IF(DDG4.LE.DDG2)THEN
        XMB1(-1) = DDG4
      ELSE
        XMB1(-1) = DDG2
      ENDIF
      DDH2  = DDH(-1)
      DDH4  = 0.
      IF(DDH4.LE.DDH2)THEN
        XMH1(-1) = DDH4
      ELSE
        XMH1(-1) = DDH2
      ENDIF
C
      DO 34 J = 0,NJ2
      GJP(J)  = DGP (J)   - 0.5 * XMB1(J)
      GJM(J)  = DGP (J-1) + 0.5 * XMB1(J-1)
      HJP(J)  = DHP (J)   - 0.5 * XMH1(J)
      HJM(J)  = DHP (J-1) + 0.5 * XMH1(J-1)
   34 CONTINUE
C
      GJP (-1) = GJP (0)
      GJM (-1) = GJM (0)
      HJP (-1) = HJP (0)
      HJM (-1) = HJM (0)
C
      DO 40 J = -1,NJ2
      IF(GJP(J).GT.0.)THEN
        AP    = 1.
      ELSE
        AP    = -1.
      ENDIF
      IF(GJM(J).GT.0.)THEN
        AM    = 1.
      ELSE
        AM    = -1.
      ENDIF
      APP     = AP*GJP(J)
      AMM     = AM*GJM(J)
      GJJ(J)  = 0.5*(AP+AM)*MIN(APP,AMM)
      IF(HJP(J).GT.0.)THEN
        HP    = 1.
      ELSE
        HP    = -1.
      ENDIF
      IF(GJM(J).GT.0.)THEN
        HM    = 1.
      ELSE
        HM    = -1.
      ENDIF
      HPP     = HP*HJP(J)
      HMM     = HM*HJM(J)
      HJJ(J)  = 0.5*(HP+HM)*MIN(HPP,HMM)
   40 CONTINUE
      IF(LSTD.EQ.0)THEN
      DO 50 J = -1,NJ2
      GI(J)   = 0.5 * (ABS(CET(I,J))
     1               -DTLOC(I,J)*CET(I,J)*CET(I,J))*GJJ(J)
      HI(J)   = 0.5 * (ABS(CET(I,J))
     1               -DTLOC(I,J)*CET(I,J)*CET(I,J))*HJJ(J)
   50 CONTINUE
      ELSE
      DO 52 J = -1,NJ2
      GI(J)   = 0.5 * ABS(CET(I,J))*GJJ(J)
      HI(J)   = 0.5 * ABS(CET(I,J))*HJJ(J)
   52 CONTINUE
      ENDIF
C
      DO 60 J = -1,NJ1
      IF(DGP(J).EQ.0.)THEN
       GAM  = 0.
      ELSE
       GAM  = (GI(J+1) - GI(J))/DGP(J)
      ENDIF
      CC0   = CET(I,J) + GAM
      PSI0  = ABS(CC0)
      PSI1  = 0.5*(CC0*CC0 + EPS*EPS)/EPS
      PSI   = CVMGT(PSI0,PSI1,ABS(CC0) .GE. EPS)
      FG(J) = 0.5*(CET(I,J)*RG(I,J)
     &            +CET(I,J+1)*RG(I,J+1)
     &            +GI(J)+GI(J+1)-PSI*DGP(J))
      IF(DHP(J).EQ.0.)THEN
       GAH  = 0.
      ELSE
       GAH  = (HI(J+1) - HI(J))/DHP(J)
      ENDIF
      CCH   = CET(I,J) + GAH
      PSH0  = ABS(CCH)
      PSH1  = 0.5*(CCH*CCH + EPS*EPS)/EPS
      PSH   = CVMGT(PSH0,PSH1,ABS(CCH) .GE. EPS)
      FH(J) = 0.5*(CET(I,J)*RH(I,J)
     &            +CET(I,J+1)*RH(I,J+1)
     &            +HI(J)+HI(J+1)-PSI*DHP(J))
   60 CONTINUE
C
      DO 70 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FG(J) - FG(J-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FH(J) - FH(J-1))
   70 CONTINUE
C
   10 CONTINUE
      RETURN
      END

      SUBROUTINE FVJTVD2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,
     1                            FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR /ENTEPS
C
      DIMENSION  ETAX(-1:NI2,-1:NJ2),ETAY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2),  RH(-1:NI2,-1:NJ2)
      DIMENSION    CET(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION  GI(-1:MAXIJP),DG(-1:MAXIJP)
      DIMENSION  HI(-1:MAXIJP),DH(-1:MAXIJP)
      DIMENSION  ADG(-1:MAXIJP),ADH(-1:MAXIJP)
      DIMENSION  AS(-1:MAXIJP),AH(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP), DPFH(-1:MAXIJP)
      DIMENSION  FXG(-1:MAXIJP),  FXH(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  CXH(-1:MAXIJP), ACXH(-1:MAXIJP)
      DIMENSION  PSG(-1:MAXIJP), PSH(-1:MAXIJP)
      DIMENSION  SIG(-1:MAXIJP), SIH(-1:MAXIJP)
      DIMENSION  GAG(-1:MAXIJP), GAH(-1:MAXIJP)
      DIMENSION  AGAG(-1:MAXIJP), AGAH(-1:MAXIJP)
C
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = EPS*2.
C
      DO 10 I = -1,NI2
C
      DO 20 J = -1,NJ2
      FXG(J)   = CET(I,J)*RG(I,J)
      FXH(J)   = CET(I,J)*RH(I,J)
   20 CONTINUE
C     
      DO 30 J = -1,NJ1
      DG(J) = 2.*(RG(I,J+1)*QJ(I,J+1)-RG(I,J)*QJ(I,J))
     1                    /(QJ(I,J+1)+QJ(I,J))
      DH(J) = 2.*(RH(I,J+1)*QJ(I,J+1)-RH(I,J)*QJ(I,J))
     1                    /(QJ(I,J+1)+QJ(I,J))
C     DH(J) = RH(I,J+1) - RH(I,J)
      DPFG(J) = FXG(J+1)   - FXG(J)
      DPFH(J) = FXH(J+1)   - FXH(J)
   30 CONTINUE
      DG(NJ2)  = DG(NJ1)
      DH(NJ2)  = DH(NJ1)
      DPFG(NJ2) = DPFG(NJ1)
      DPFH(NJ2) = DPFH(NJ1)
C
      DO 40 J = -1,NJ1
      IF(DG(J) .NE. 0.) THEN
      CXG(J) = DPFG(J)/DG(J)
      ELSE
      CXG(J) = 0.5*(CET(I,J+1) + CET(I,J))
      ENDIF
      IF(DH(J) .NE. 0.) THEN
      CXH(J) = DPFH(J)/DH(J)
      ELSE
      CXH(J) = 0.5*(CET(I,J+1) + CET(I,J))
      ENDIF
   40 CONTINUE
      CXG(NJ2) = CXG(NJ1)
      CXH(NJ2) = CXG(NJ1)
C
      DO 50 J = -1,NJ2
      IF(DG(J).GE.0.)THEN
      AS(J) = 1.
      ELSE
      AS(J) =-1.
      ENDIF
      IF(DH(J).GE.0.)THEN
      AH(J) = 1.
      ELSE
      AH(J) =-1.
      ENDIF
      ADG(J) = ABS(DG(J))
      ADH(J) = ABS(DH(J))
      ACXG(J) = ABS(CXG(J))
      ACXH(J) = ABS(CXH(J))
      PSG(J)  = ACXG(J)
      PSH(J)  = ACXH(J)
      IF(ACXG(J) .LT. EPS) PSG(J) = (CXG(J)*CXG(J)+EPSS)/EPS2
      IF(ACXH(J) .LT. EPS) PSH(J) = (CXH(J)*CXH(J)+EPSS)/EPS2
   50 CONTINUE
      DO 60 J = 0,NJ2
      GI(J) = 0.5*(AS(J-1)+AS(J))*MIN(ADG(J-1),ADG(J))
      HI(J) = 0.5*(AH(J-1)+AH(J))*MIN(ADH(J-1),ADH(J))
   60 CONTINUE
      GI( -1) = GI(  0)
      HI( -1) = HI(  0)
C
      DO 70 J = -1,NJ2
      IF(LSTD.EQ.0)THEN
      SIG(J) = 0.5*(PSG(J)
     1             - DTLOC(I,J)*CXG(J)*CXG(J))*DG(J)
      SIH(J) = 0.5*(PSH(J)
     1             - DTLOC(I,J)*CXH(J)*CXH(J))*DH(J)
      ELSE
      SIG(J) = 0.5*PSG(J)*DG(J)
      SIH(J) = 0.5*PSH(J)*DH(J)
      ENDIF
   70 CONTINUE
C
      DO 80 J = -1,NJ1
      IF(DG(J).NE.0.)THEN
      GAG(J)  = SIG(J)*(GI(J+1)-GI(J))/DG(J)
      ELSE
      GAG(J)  = 0.
      ENDIF
      IF(DH(J).NE.0.)THEN
      GAH(J)  = SIH(J)*(HI(J+1)-HI(J))/DH(J)
      ELSE
      GAH(J)  = 0.
      ENDIF
   80 CONTINUE
      GAG(NJ2) = GAG(NJ1)
      GAH(NJ2) = GAH(NJ1)
C
      DO 90 J = -1,NJ2
      AGAG(J)  = ABS(CXG(J) + GAG(J))
      AGAH(J)  = ABS(CXH(J) + GAH(J))
      PSG(J)   = AGAG(J)
      PSH(J)   = AGAH(J)
      IF(AGAG(J) .LT. EPS) PSG(J) = (AGAG(J)*AGAG(J)+EPSS)/EPS2
      IF(AGAH(J) .LT. EPS) PSH(J) = (AGAH(J)*AGAH(J)+EPSS)/EPS2
   90 CONTINUE
C
      DO 100 J = -1,NJ1
      FG(J)  = 0.5*(FXG(J+1)+FXG(J)+SIG(J)*(GI(J+1)+GI(J))
     1              -PSG(J)*DG(J))
      FH(J)  = 0.5*(FXH(J+1)+FXH(J)+SIH(J)*(HI(J+1)+HI(J))
     1              -PSH(J)*DH(J))
  100 CONTINUE
      DO 110 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FG(J) - FG(J-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FH(J) - FH(J-1))
  110 CONTINUE
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FVJUPW(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,
     1                                      FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
C
      DIMENSION  ETAX(-1:NI2,-1:NJ2),ETAY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2),  RH(-1:NI2,-1:NJ2)
      DIMENSION   CET(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DO  40 I = -1,NI2
      DO  50 J = -1,NJ1
C
      DG1 = RG(I,J+1) + RG(I,J)
      DG2 = RG(I,J+1) - RG(I,J)
      DH1 = RH(I,J+1) + RH(I,J)
      DH2 = RH(I,J+1) - RH(I,J)
C
      FG(J)  = 0.5*(CET(I,J)*DG1 - ABS(CET(I,J))*DG2)
      FH(J)  = 0.5*(CET(I,J)*DH1 - ABS(CET(I,J))*DH2)
   50 CONTINUE
C
C     DO  40 I = -1,NI2
C     DO  50 J = -1,NJ1
C
C     QJ2 = 2./(QJ(I,J+1)+QJ(I,J))
C     DG2 = (RG(I,J+1)*QJ(I,J+1) - RG(I,J)*QJ(I,J))*QJ2
C     DH2 = (RH(I,J+1)*QJ(I,J+1) - RH(I,J)*QJ(I,J))*QJ2
C     DFG = (RG(I,J+1)*QJ(I,J+1)*CET(I,J+1) 
C    1      - RG(I,J)*QJ(I,J)*CET(I,J))*QJ2
C     DFH = (RH(I,J+1)*QJ(I,J+1)*CET(I,J+1) 
C    1      - RH(I,J)*QJ(I,J)*CET(I,J))*QJ2
C     IF(DG2 .NE. 0.)THEN
C     CYG = DFG/DG2
C     ELSE
C     CYG = 0.5*(CET(I,J+1)+CET(I,J))
C     ENDIF
C     IF(DH2 .NE. 0.)THEN
C     CYH = DFH/DH2
C     ELSE
C     CYH = 0.5*(CET(I,J+1)+CET(I,J))
C     ENDIF
C
C     IF(CYG .GE. 0.)THEN
C     FG(J)  = CET(I,J)*RG(I,J)*QJ(I,J)*QJ2
C     ELSE
C     FG(J)  = CET(I,J+1)*RG(I,J+1)*QJ(I,J+1)*QJ2
C     ENDIF
C     IF(CYG .GE. 0.)THEN
C     FH(J)  = CET(I,J)*RH(I,J)*QJ(I,J)*QJ2
C     ELSE
C     FH(J)  = CET(I,J+1)*RH(I,J+1)*QJ(I,J+1)*QJ2
C     ENDIF
C  50 CONTINUE
C
C     DO  40 I = -1,NI2
C     DO  50 J = -1,NJ1
C
C     QJ2 = 2./(QJ(I,J+1)+QJ(I,J))
C     DG1 = (RG(I,J+1)*QJ(I,J+1) + RG(I,J)*QJ(I,J))*QJ2
C     DG2 = (RG(I,J+1)*QJ(I,J+1) - RG(I,J)*QJ(I,J))*QJ2
C     DH1 = (RH(I,J+1)*QJ(I,J+1) + RH(I,J)*QJ(I,J))*QJ2
C     DH2 = (RH(I,J+1)*QJ(I,J+1) - RH(I,J)*QJ(I,J))*QJ2
C     DFG = (RG(I,J+1)*QJ(I,J+1)*CET(I,J+1) 
C    1      - RG(I,J)*QJ(I,J)*CET(I,J))*QJ2
C     DFH = (RH(I,J+1)*QJ(I,J+1)*CET(I,J+1) 
C    1      - RH(I,J)*QJ(I,J)*CET(I,J))*QJ2
C     IF(DG2 .NE. 0.)THEN
C     CYG = DFG/DG2
C     ELSE
C     CYG = 0.5*(CET(I,J+1)+CET(I,J))
C     ENDIF
C     IF(DH2 .NE. 0.)THEN
C     CYH = DFH/DH2
C     ELSE
C     CYH = 0.5*(CET(I,J+1)+CET(I,J))
C     ENDIF
C
C     FG(J)  = 0.5*(CYG*DG1 - ABS(CYG)*DG2)
C     FH(J)  = 0.5*(CYH*DH1 - ABS(CYH)*DH2)
C  50 CONTINUE
C
      DO 60 J = 0,NJ1
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FG(J) - FG(J-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FH(J) - FH(J-1))
   60 CONTINUE
C
   40 CONTINUE
      RETURN
      END
      SUBROUTINE GETBC(ID,RG,RH,QJ,GBC,HBC,K,L)
      PARAMETER (MAXBC=20)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
      COMMON /BC02 / IDA,ISA,ICA,IAS,IAE,IPA,LENA,
     1               IDB,ISB,ICB,IBS,IBE,IPB,LENB
      DIMENSION   ID(LENA,10)
      DIMENSION   RG(-2:NI3,-2:NJ3),  RH(-2:NI3,-2:NJ3)
      DIMENSION   QJ(-2:NI3,-2:NJ3)
      DIMENSION  GBC(LENA,2,NVX,NVY),HBC(LENA,2,NVX,NVY)
C
      DO 100 II = 1 , LENA
      I1        = ID(II,1)
      J1        = ID(II,2)
      I2        = ID(II,3)
      J2        = ID(II,4)
      RG(I1,J1) = GBC(II,1,K,L)/QJ(I1,J1)
      RG(I2,J2) = GBC(II,2,K,L)/QJ(I2,J2)
      RH(I1,J1) = HBC(II,1,K,L)/QJ(I1,J1)
      RH(I2,J2) = HBC(II,2,K,L)/QJ(I2,J2)
  100 CONTINUE
      RETURN
      END

      SUBROUTINE GRIDIN( X , Y ,IJP)
      PARAMETER (MAXBC=20)
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /TAPE / INAME,IGRID,IWRIT,IREAD,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
C
      DIMENSION X(-2:NI3,-2:NJ3),Y(-2:NI3,-2:NJ3)
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
      IPA    = IDBC(25,IBC)
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
      ELSE
      IDB    = IDBC(12,IBC)
      ISB    = IDBC(13,IBC)
      ICB    = IDBC(14,IBC)
      I3S    = IDBC(15,IBC)
      I3E    = IDBC(16,IBC)
      IPB    = IDBC(26,IBC)
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
      WRITE(69,71)IAI,IAJ,IBI,IBJ,ICI,ICJ,IDI,IDJ,IEI,IEJ,
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
C     X(I  ,J  ) =   2.*X(I-1,J-1) - X(I-2,J-2)
C     Y(I  ,J  ) =   2.*Y(I-1,J-1) - Y(I-2,J-2)
C     X(I+1,J+1) =   2.*X(I  ,J  ) - X(I-1,J-1)
C     Y(I+1,J+1) =   2.*Y(I  ,J  ) - Y(I-1,J-1)
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
C     X(I  ,J  ) =   2.*X(I+1,J-1) - X(I+2,J-2)
C     Y(I  ,J  ) =   2.*Y(I+1,J-1) - Y(I+2,J-2)
C     X(I-1,J+1) =   2.*X(I  ,J  ) - X(I+1,J-1)
C     Y(I-1,J+1) =   2.*Y(I  ,J  ) - Y(I+1,J-1)
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
C     X(I  ,J  ) =   2.*X(I-1,J+1) - X(I-2,J+2)
C     Y(I  ,J  ) =   2.*Y(I-1,J+1) - Y(I-2,J+2)
C     X(I+1,J-1) =   2.*X(I  ,J  ) - X(I-1,J+1)
C     Y(I+1,J-1) =   2.*Y(I  ,J  ) - Y(I-1,J+1)
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
C     
C     X(I  ,J  ) =   2.*X(I+1,J+1) - X(I+2,J+2)
C     Y(I  ,J  ) =   2.*Y(I+1,J+1) - Y(I+2,J+2)
C     X(I-1,J-1) =   2.*X(I  ,J  ) - X(I+1,J+1)
C     Y(I-1,J-1) =   2.*Y(I  ,J  ) - Y(I+1,J+1)
C
      RETURN
      END

      SUBROUTINE INITG(X,A,B,D,T,U,V,P,QX,QY,G,H,QJ)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX0 / MAXW,MAXI,MAXB,MAXG,MAXL,MAXV,MAXBP,MAXLUI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /CNST2/ D0,T0,P0,U0,V0,TAU0,GAMMA,PI,SPI
      COMMON /UNST0/ URGHT,VRGHT,DRGHT,TRGHT,PRGHT,
     1               ULEFT,VLEFT,DLEFT,TLEFT,PLEFT
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      DIMENSION    X(-2:NI3,-2:NJ3)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION   QX(-2:NI3,-2:NJ3),  QY(-2:NI3,-2:NJ3)
      DIMENSION    P(-2:NI3,-2:NJ3)
      DIMENSION   QJ(-2:NI3,-2:NJ3)
      DIMENSION    G(NI,NJ,NVX,NVY),   H(NI,NJ,NVX,NVY)
      DIMENSION    A(NVX),B(NVY)
C
      IF(LSTD.EQ.0)THEN
      DO 20 I = 1,NI
      DO 20 J = 1,NJ
      DO 20 K = 1,NVX
      DO 20 L = 1,NVY
      IF(X(I,J).GE.XST)THEN
      PP   = -1./TRGHT*((A(K)-URGHT)*(A(K)-URGHT)+
     1                  (B(L)-VRGHT)*(B(L)-VRGHT))
      GINF = DRGHT/(PI*TRGHT)*EXP(PP)
      HINF = 0.5 * TRGHT * GINF
      G(I,J,K,L) = GINF
      H(I,J,K,L) = HINF
      ELSE
      PP   = -1./TLEFT*((A(K)-ULEFT)*(A(K)-ULEFT)
     1                 +(B(L)-VLEFT)*(B(L)-VLEFT))
      GINF = DLEFT/(PI*TLEFT)*EXP(PP)
      HINF = 0.5 * TLEFT * GINF
      G(I,J,K,L) = GINF
      H(I,J,K,L) = HINF
      ENDIF
   20 CONTINUE
      ELSE
      DO 10 K = 1,NVX
      DO 10 L = 1,NVY
      PP   = -1./T0*((A(K)-U0)*(A(K)-U0)+(B(L)-V0)*(B(L)-V0))
      GINF = D0/(PI*T0)*EXP(PP)
      HINF = 0.5 * T0 * GINF
      DO 10 I = 1,NI
      DO 10 J = 1,NJ
      G(I,J,K,L) = GINF
      H(I,J,K,L) = HINF
   10 CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE INITQ(X,D,T,U,V,P,QX,QY)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX0 / MAXW,MAXI,MAXB,MAXG,MAXL,MAXV,MAXBP,MAXLUI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /CNST2/ D0,T0,P0,U0,V0,TAU0,GAMMA,PI,SPI
      COMMON /UNST0/ URGHT,VRGHT,DRGHT,TRGHT,PRGHT,
     1               ULEFT,VLEFT,DLEFT,TLEFT,PLEFT
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      DIMENSION    X(-2:NI3,-2:NJ3)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION   QX(-2:NI3,-2:NJ3),  QY(-2:NI3,-2:NJ3)
      DIMENSION    P(-2:NI3,-2:NJ3)
C
      IF(LSTD.EQ.0)THEN
      DO 20 I=-2,NI3
      DO 20 J=-2,NJ3
      IF(X(I,J).GE.XST)THEN
      D(I,J) = DRGHT
      T(I,J) = TRGHT
      U(I,J) = URGHT
      V(I,J) = VRGHT
      P(I,J) = PRGHT
      QX(I,J) = 0.
      QY(I,J) = 0.
      ELSE
      D(I,J) = DLEFT
      T(I,J) = TLEFT
      U(I,J) = ULEFT
      V(I,J) = VLEFT
      P(I,J) = PLEFT
      QX(I,J) = 0.
      QY(I,J) = 0.
      ENDIF
   20 CONTINUE
      ELSE
      DO 10 I=-2,NI3
      DO 10 J=-2,NJ3
      D(I,J) = D0
      T(I,J) = T0
      U(I,J) = U0
      V(I,J) = V0
      P(I,J) = P0
      QX(I,J) = 0.
      QY(I,J) = 0.
   10 CONTINUE
      ENDIF
C
      RETURN
      END
      SUBROUTINE LDTSET(DTLOC)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      DIMENSION DTLOC(-2:NI3,-2:NJ3)
      MIJ = MI*MJ-2  !有問題
      DO 10 IJ = -2,MIJ
      DTLOC(IJ,-2) = DTCFL
   10 CONTINUE
      RETURN
      END

      SUBROUTINE LUPNT(KDI,KDJ,LKD,KDP)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /PNTI / LKDI,LKDJ,LLKD,LKDP,MDA,MDB,NDIG
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      DIMENSION KDI(MDA),KDJ(MDA),LKD(MDB),KDP(MDB)
C
      NEXT = 1
      NDIG = MI + MJ -1
      MIS   = -2
      MJS   = -2
      MIE   = NI3
      MJE   = NJ3
      DO 20 L = 1 , NDIG
      LSUM    = L + MIS + MJS - 1
      LENG    = 0
      KDP(L) = NEXT
      DO 10 J = MJS,MJE
      DO 10 I = MIS,MIE
      IF( (I+J) .EQ. LSUM ) THEN
          LENG = LENG + 1
          KDI(NEXT+LENG-1) = I
          KDJ(NEXT+LENG-1) = J
      ENDIF
   10 CONTINUE
      LKD(L)  = LENG
      NEXT    = NEXT + LENG
   20 CONTINUE
      RETURN
      END
      SUBROUTINE MACRO(DTLOC,A,B,WA,WB,X,Y,XIX,XIY,ETAX,ETAY,QJ,
     1                D,U,V,T,P,QX,QY,G,H,IJP,SLXX,SLYY,SLZZ,SLXY)
      PARAMETER (NV=8000)
      PARAMETER (MAXBC=20)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /VELC / W0(NV),W1(NV),W2(NV),W3(NV),
     1               A3(NV),A4(NV),A5(NV),A6(NV)
      COMMON /RES0/ SUMD,Z1MAX,M1DI,M1DJ,Z2MAX,M2DI,M2DJ,SUMDG,SUMDH
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
      COMMON /BC02 / IDA,ISA,ICA,IAS,IAE,IPA,LENA,
     1               IDB,ISB,ICB,IBS,IBE,IPB,LENB
      COMMON /BC03 / TWAL,TACC,UACC
      COMMON /FINP1/ KMOD,MMOD,ZETA,PR,VISS,XLDA
      DIMENSION  IJP(MXBP)
      DIMENSION    X(-2:NI3,-2:NJ3),   Y(-2:NI3,-2:NJ3)
      DIMENSION  XIX(-2:NI3,-2:NJ3), XIY(-2:NI3,-2:NJ3)
      DIMENSION ETAX(-2:NI3,-2:NJ3),ETAY(-2:NI3,-2:NJ3)
      DIMENSION   QJ(-2:NI3,-2:NJ3)
      DIMENSION   QX(-2:NI3,-2:NJ3),  QY(-2:NI3,-2:NJ3)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION    P(-2:NI3,-2:NJ3)
      DIMENSION DTLOC(-2:NI3,-2:NJ3)
      DIMENSION    G(NI,NJ,NVX,NVY),   H(NI,NJ,NVX,NVY)
      DIMENSION    A(NVX),B(NVY),WA(NVX),WB(NVY)
      DIMENSION SLXX(-2:NI3,-2:NJ3),SLYY(-2:NI3,-2:NJ3)
      DIMENSION SLZZ(-2:NI3,-2:NJ3),SLXY(-2:NI3,-2:NJ3)
C
      SUMD = 0.
      Z1MAX = -1.E5
      Z2MAX = -1.E5
      NIM1  = NI - 1
      NJM1  = NJ - 1
      DO 100 J = 2 , NJM1
      DO 100 I = 2 , NIM1
         SD     = 0.
         SU     = 0.
         SV     = 0.
         ST     = 0.
         SQX    = 0.
         SQY    = 0.
         SAA    = 0.
         SBB    = 0.
         SAB    = 0.
      DO 110 KL = 1,MXV
         K      = MOD(KL-1,NVX)+1
         L      = (KL-1)/NVX + 1
         SD     = SD   + W0(KL)*G(I,J,K,L)
         SU     = SU   + W1(KL)*G(I,J,K,L)
         SV     = SV   + W2(KL)*G(I,J,K,L)
         ST     = ST   + W3(KL)*G(I,J,K,L) + W0(KL)*H(I,J,K,L)
         AA0    = H(I,J,K,L)   + A5(KL)*G(I,J,K,L)
         SQX    = SQX  + W1(KL)*AA0
         SQY    = SQY  + W2(KL)*AA0
         SAA    = SAA  + A3(KL)*G(I,J,K,L)
         SBB    = SBB  + A4(KL)*G(I,J,K,L)
         SAB    = SAB  + A6(KL)*G(I,J,K,L)
  110 CONTINUE
         DD0     = ABS(SD - D(I,J))
         DD1     = (SD - D(I,J))**2
         DD2     = DD0/SD/DTLOC(I,J)
         IF(DD0.GT.Z1MAX)THEN
           Z1MAX = DD0
           M1DI  = I
           M1DJ  = J
         ENDIF
         IF(DD2.GT.Z2MAX)THEN
           Z2MAX = DD2
           M2DI  = I
           M2DJ  = J
         ENDIF
         SUMD   = SUMD  + DD1
         D(I,J) = SD
         U(I,J) = SU/SD
         V(I,J) = SV/SD
         T(I,J) = (ST - SD*(U(I,J)*U(I,J) + V(I,J)*V(I,J)))/SD*2./3.
         P(I,J) = D(I,J) * T(I,J)
         UV2    = U(I,J) *U(I,J) + V(I,J)*V(I,J) - 1.5*T(I,J)
         QX(I,J)= SQX - 2.*(SAA*U(I,J) + SAB*V(I,J)) + SU*UV2
         QY(I,J)= SQY - 2.*(SBB*V(I,J) + SAB*U(I,J)) + SV*UV2
      IF(KMOD.EQ.3)THEN
         SSXX   = 0.
         SSYY   = 0.
         SSZZ   = 0.
         SSXY   = 0.
      DO 130 KL = 1,MXV
         K      = MOD(KL-1,NVX)+1
         L      = (KL-1)/NVX + 1
         U1     = A(K)  - U(I,J)
         V1     = B(L)  - V(I,J)
         SSXX   = SSXX  + A3(KL)*G(I,J,K,L)
         SSYY   = SSYY  + A4(KL)*G(I,J,K,L)
         SSZZ   = SSZZ  + W0(KL)*H(I,J,K,L)
         SSXY   = SSXY  + A6(KL)*G(I,J,K,L)
  130 CONTINUE
         SLXX(I,J) = SSXX
         SLYY(I,J) = SSYY
         SLZZ(I,J) = SSZZ
         SLXY(I,J) = SSXY
      ENDIF
  100 CONTINUE
C
      DO 200 IBC = 1,MXB
      IBTYPE = IDBC( 1,IBC)
      IDA    = IDBC( 2,IBC)
      ISA    = IDBC( 3,IBC)
      ICA    = IDBC( 4,IBC)
      IAS    = IDBC( 5,IBC)
      IAE    = IDBC( 6,IBC)
      IPA    = IDBC(25,IBC)
      LENA   = (IAE - IAS)*ICA + 1
      IF(IBTYPE .EQ. 4) THEN
      IDB    = IDBC(12,IBC)
      ISB    = IDBC(13,IBC)
      ICB    = IDBC(14,IBC)
      IBS    = IDBC(15,IBC)
      IBE    = IDBC(16,IBC)
      IPB    = IDBC(26,IBC)
      LENB   = (IBE - IBS)*ICB + 1
      ENDIF
      IF(IBTYPE .EQ. 1 ) GO TO 250
      DO 210 II = 1 , LENA
          IC    =  IJP(IPA+LENA*4+II)
          JC    =  IJP(IPA+LENA*5+II)
           I    =  IC
           J    =  JC
         SD     = 0.
         SU     = 0.
         SV     = 0.
         ST     = 0.
         SQX    = 0.
         SQY    = 0.
         SAA    = 0.
         SBB    = 0.
         SAB    = 0.
      DO 220 KL = 1,MXV
         K      = MOD(KL-1,NVX)+1
         L      = (KL-1)/NVX + 1
         SD     = SD   + W0(KL)*G(I,J,K,L)
         SU     = SU   + W1(KL)*G(I,J,K,L)
         SV     = SV   + W2(KL)*G(I,J,K,L)
         ST     = ST   + W3(KL)*G(I,J,K,L) + W0(KL)*H(I,J,K,L)
         AA0    = H(I,J,K,L)   + A5(KL)*G(I,J,K,L)
         SQX    = SQX  + W1(KL)*AA0
         SQY    = SQY  + W2(KL)*AA0
         SAA    = SAA  + A3(KL)*G(I,J,K,L)
         SBB    = SBB  + A4(KL)*G(I,J,K,L)
         SAB    = SAB  + A6(KL)*G(I,J,K,L)
  220 CONTINUE
         D(I,J) = SD
         U(I,J) = SU/SD
         V(I,J) = SV/SD
         T(I,J) = (ST - SD*(U(I,J)*U(I,J) + V(I,J)*V(I,J)))/SD*2./3.
         P(I,J) = D(I,J) * T(I,J)
         UV2    = U(I,J) *U(I,J) + V(I,J)*V(I,J) - 1.5*T(I,J)
         QX(I,J)= SQX - 2.*(SAA*U(I,J) + SAB*V(I,J)) + SU*UV2
         QY(I,J)= SQY - 2.*(SBB*V(I,J) + SAB*U(I,J)) + SV*UV2
  210 CONTINUE
      IF(KMOD.EQ.3)THEN
         SSXX   = 0.
         SSYY   = 0.
         SSZZ   = 0.
         SSXY   = 0.
      DO 120 KL = 1,MXV
         K      = MOD(KL-1,NVX)+1
         L      = (KL-1)/NVX + 1
         U1     = A(K)  - U(I,J)
         V1     = B(L)  - V(I,J)
         SSXX   = SSXX  + A3(KL)*G(I,J,K,L)
         SSYY   = SSYY  + A4(KL)*G(I,J,K,L)
         SSZZ   = SSZZ  + W0(KL)*H(I,J,K,L)
         SSXY   = SSXY  + A6(KL)*G(I,J,K,L)
  120 CONTINUE
         SLXX(I,J) = SSXX
         SLYY(I,J) = SSYY
         SLZZ(I,J) = SSZZ
         SLXY(I,J) = SSXY
      ENDIF
C
      DO 230 II = 1 , LENA
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
      D(IB,JB)  =  D(IC,JC)
      U(IB,JB)  =  U(IC,JC)
      V(IB,JB)  =  V(IC,JC)
      T(IB,JB)  =  T(IC,JC)
      P(IB,JB)  =  P(IC,JC)
      QX(IB,JB) = QX(IC,JC)
      QY(IB,JB) = QY(IC,JC)
      D(IA,JA)  =  D(IC,JC)
      U(IA,JA)  =  U(IC,JC)
      V(IA,JA)  =  V(IC,JC)
      T(IA,JA)  =  T(IC,JC)
      P(IA,JA)  =  P(IC,JC)
      QX(IA,JA) = QX(IC,JC)
      QY(IA,JA) = QY(IC,JC)
      IF(KMOD.EQ.3)THEN
      SLXX(IB,JB) = SLXX(IC,JC)
      SLYY(IB,JB) = SLYY(IC,JC)
      SLZZ(IB,JB) = SLZZ(IC,JC)
      SLXY(IB,JB) = SLXY(IC,JC)
      SLXX(IA,JA) = SLXX(IC,JC)
      SLYY(IA,JA) = SLYY(IC,JC)
      SLZZ(IA,JA) = SLZZ(IC,JC)
      SLXY(IA,JA) = SLXY(IC,JC)
      ENDIF
      ELSEIF(IBTYPE .EQ. 2 .OR. IBTYPE .EQ. 3) THEN
      D(IB,JB)  =  D(ID,JD)
      U(IB,JB)  =  U(ID,JD)
      V(IB,JB)  = -V(ID,JD)
      T(IB,JB)  =  T(ID,JD)
      P(IB,JB)  =  P(ID,JD)
      QX(IB,JB) = QX(ID,JD)
      QY(IB,JB) =-QY(ID,JD)
      D(IA,JA)  =  D(IE,JE)
      U(IA,JA)  =  U(IE,JE)
      V(IA,JA)  = -V(IE,JE)
      T(IA,JA)  =  T(IE,JE)
      P(IA,JA)  =  P(IE,JE)
      QX(IA,JA) = QX(IE,JE)
      QY(IA,JA) =-QY(IE,JE)
      V(IC,JC)  = 0.
      QY(IC,JC) = 0.
      IF(KMOD.EQ.3)THEN
      SLXX(IB,JB) = SLXX(IC,JC)
      SLYY(IB,JB) = SLYY(IC,JC)
      SLZZ(IB,JB) = SLZZ(IC,JC)
      SLXY(IB,JB) = SLXY(IC,JC)
      SLXX(IA,JA) = SLXX(IC,JC)
      SLYY(IA,JA) = SLYY(IC,JC)
      SLZZ(IA,JA) = SLZZ(IC,JC)
      SLXY(IA,JA) = SLXY(IC,JC)
      ENDIF
      ELSEIF(IBTYPE .EQ. 4) THEN
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
      D(IB,JB)  =  D(ID2,JD2)
      U(IB,JB)  =  U(ID2,JD2)
      V(IB,JB)  =  V(ID2,JD2)
      T(IB,JB)  =  T(ID2,JD2)
      P(IB,JB)  =  P(ID2,JD2)
      QX(IB,JB) = QX(ID2,JD2)
      QY(IB,JB) = QY(ID2,JD2)
      D(IA,JA)  =  D(IE2,JE2)
      U(IA,JA)  =  U(IE2,JE2)
      V(IA,JA)  =  V(IE2,JE2)
      T(IA,JA)  =  T(IE2,JE2)
      P(IA,JA)  =  P(IE2,JE2)
      QX(IA,JA) = QX(IE2,JE2)
      QY(IA,JA) = QY(IE2,JE2)
      D(IA2,JA2)  =  D(IE,JE)
      U(IA2,JA2)  =  U(IE,JE)
      V(IA2,JA2)  =  V(IE,JE)
      T(IA2,JA2)  =  T(IE,JE)
      P(IA2,JA2)  =  P(IE,JE)
      QX(IB2,JA2) = QX(IE,JE)
      QY(IB2,JA2) = QY(IE,JE)
      D(IB2,JB2)  =  D(ID,JD)
      U(IB2,JB2)  =  U(ID,JD)
      V(IB2,JB2)  =  V(ID,JD)
      T(IB2,JB2)  =  T(ID,JD)
      P(IB2,JB2)  =  P(ID,JD)
      QX(IB2,JB2) = QX(ID,JD)
      QY(IB2,JB2) = QY(ID,JD)
C
C     FOR SIMPLE SET EQUAL
C
        D(IC,JC)   = 0.5*(  D(ID,JD)+  D(ID2,JD2))
        U(IC,JC)   = 0.5*(  U(ID,JD)+  U(ID2,JD2))
        V(IC,JC)   = 0.5*(  V(ID,JD)+  V(ID2,JD2))
        P(IC,JC)   = 0.5*(  P(ID,JD)+  P(ID2,JD2))
        T(IC,JC)   = 0.5*( T(ID,JD)+ T(ID2,JD2))
       QX(IC,JC)   = 0.5*(QX(ID,JD)+QX(ID2,JD2))
       QY(IC,JC)   = 0.5*(QY(ID,JD)+QY(ID2,JD2))

      D(IC2,JC2)  =  D(IC,JC)
      U(IC2,JC2)  =  U(IC,JC)
      V(IC2,JC2)  =  V(IC,JC)
      T(IC2,JC2)  =  T(IC,JC)
      P(IC2,JC2)  =  P(IC,JC)
      QX(IC2,JC2) = QX(IC,JC)
      QY(IC2,JC2) = QY(IC,JC)
C
      IF(KMOD .EQ. 3) THEN
      SLXX(IB,JB) = SLXX(ID2,JD2)
      SLYY(IB,JB) = SLYY(ID2,JD2)
      SLZZ(IB,JB) = SLZZ(ID2,JD2)
      SLXY(IB,JB) = SLXY(ID2,JD2)
      SLXX(IA,JA) = SLXX(IE2,JE2)
      SLYY(IA,JA) = SLYY(IE2,JE2)
      SLZZ(IA,JA) = SLZZ(IE2,JE2)
      SLXY(IA,JA) = SLXY(IE2,JE2)
      SLXX(IB2,JB2) = SLXX(ID,JD)
      SLYY(IB2,JB2) = SLYY(ID,JD)
      SLZZ(IB2,JB2) = SLZZ(ID,JD)
      SLXY(IB2,JB2) = SLXY(ID,JD)
      SLXX(IA2,JA2) = SLXX(IE,JE)
      SLYY(IA2,JA2) = SLYY(IE,JE)
      SLZZ(IA2,JA2) = SLZZ(IE,JE)
      SLXY(IA2,JA2) = SLXY(IE,JE)
      SLXX(IC,JC)   = 0.5*(SLXX(ID,JD)+SLXX(ID2,JD2))
      SLYY(IC,JC)   = 0.5*(SLYY(ID,JD)+SLYY(ID2,JD2))
      SLZZ(IC,JC)   = 0.5*(SLZZ(ID,JD)+SLZZ(ID2,JD2))
      SLXY(IC,JC)   = 0.5*(SLXY(ID,JD)+SLXY(ID2,JD2))
      SLXX(IC2,JC2) = SLXX(IC,JC)
      SLYY(IC2,JC2) = SLYY(IC,JC)
      SLZZ(IC2,JC2) = SLZZ(IC,JC)
      SLXY(IC2,JC2) = SLXY(IC,JC)
      ENDIF
      ENDIF
  230 CONTINUE
      GO TO 200
  250 CONTINUE
C
      TWAL   = FDBC( 1,IBC)
      TACC   = FDBC( 2,IBC)
      UACC   = FDBC( 3,IBC)
      IDTW   = IDBC(21,IBC)
      IDUW   = IDBC(22,IBC)
C
      IF(IDA .EQ. 1)
     1  CALL BCQWAL(IJP(IPA+1),XIX,XIY,QJ,
     1         D,U,V,T,P,QX,QY,G,H,A,B,WA,WB,SLXX,SLYY,SLZZ,SLXY)
C
      IF(IDA .EQ. 2)
     1  CALL BCQWAL(IJP(IPA+1),ETAX,ETAY,QJ,
     1         D,U,V,T,P,QX,QY,G,H,A,B,WA,WB,SLXX,SLYY,SLZZ,SLXY)
C
  200 CONTINUE
C
      RETURN
      END


      SUBROUTINE METJAC(X,Y,XIX,XIY,ETAX,ETAY,QJ,IJP)
      PARAMETER (MAXBC=20)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
      COMMON /AXIS / IAXIS
      DIMENSION    X(-2:NI3,-2:NJ3),   Y(-2:NI3,-2:NJ3)
      DIMENSION  XIX(-2:NI3,-2:NJ3), XIY(-2:NI3,-2:NJ3)
      DIMENSION ETAX(-2:NI3,-2:NJ3),ETAY(-2:NI3,-2:NJ3)
      DIMENSION   QJ(-2:NI3,-2:NJ3)
      DIMENSION IJP(MXBP)
C
      DO 10 I = -1,NI2
      DO 10 J = -1,NJ2
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
      J       = -2
      DO 20 I = -1 , NI2
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
      J       = NJ3
      DO 30 I = -1 , NI2
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
      I       = -2
C
      DO 40 J = -1 , NJ2
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
      I       = NI3
      DO 50 J = -1 , NJ2
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
      I = -2
      J = -2
         XI   = -(3.0*X(I,J) -4.0*X(I+1,J) +X(I+2,J))*0.5
         YI   = -(3.0*Y(I,J) -4.0*Y(I+1,J) +Y(I+2,J))*0.5
         XJ   = -(3.0*X(I,J) -4.0*X(I,J+1) +X(I,J+2))*0.5
         YJ   = -(3.0*Y(I,J) -4.0*Y(I,J+1) +Y(I,J+2))*0.5
          QJ(I,J)  = 1./(VOL + 1.E-25)
         XIX(I,J)  =  YJ*QJ(I,J)
         XIY(I,J)  = -XJ*QJ(I,J)
        ETAX(I,J)  = -YI*QJ(I,J)
        ETAY(I,J)  =  XI*QJ(I,J)
C
      I = NI3
      J = -2
         XI   =  (3.0*X(I,J) -4.0*X(I-1,J) +X(I-2,J))*0.5
         YI   =  (3.0*Y(I,J) -4.0*Y(I-1,J) +Y(I-2,J))*0.5
         XJ   = -(3.0*X(I,J) -4.0*X(I,J+1) +X(I,J+2))*0.5
         YJ   = -(3.0*Y(I,J) -4.0*Y(I,J+1) +Y(I,J+2))*0.5
          QJ(I,J)  = 1./(VOL + 1.E-25)
         XIX(I,J)  =  YJ*QJ(I,J)
         XIY(I,J)  = -XJ*QJ(I,J)
        ETAX(I,J)  = -YI*QJ(I,J)
        ETAY(I,J)  =  XI*QJ(I,J)
C
      I = -2
      J = NJ3
         XI   = -(3.0*X(I,J) -4.0*X(I+1,J) +X(I+2,J))*0.5
         YI   = -(3.0*Y(I,J) -4.0*Y(I+1,J) +Y(I+2,J))*0.5
         XJ   =  (3.0*X(I,J) -4.0*X(I,J-1) +X(I,J-2))*0.5
         YJ   =  (3.0*Y(I,J) -4.0*Y(I,J-1) +Y(I,J-2))*0.5
          QJ(I,J)  = 1./(VOL + 1.E-25)
         XIX(I,J)  =  YJ*QJ(I,J)
         XIY(I,J)  = -XJ*QJ(I,J)
        ETAX(I,J)  = -YI*QJ(I,J)
        ETAY(I,J)  =  XI*QJ(I,J)
C
      I = NI3
      J = NJ3
         XI   =  (3.0*X(I,J) -4.0*X(I-1,J) +X(I-2,J))*0.5
         YI   =  (3.0*Y(I,J) -4.0*Y(I-1,J) +Y(I-2,J))*0.5
         XJ   =  (3.0*X(I,J) -4.0*X(I,J-1) +X(I,J-2))*0.5
         YJ   =  (3.0*Y(I,J) -4.0*Y(I,J-1) +Y(I,J-2))*0.5
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
      IPA    = IDBC(25,IBC)
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
      IF(IBTYPE .EQ. 1 ) THEN
      DO 130 I2 = 1,LENA
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
C         QJ(IAI,IAJ)  =    QJ(IEI,IEJ)
C        XIX(IAI,IAJ)  =  -XIX(IEI,IEJ)
C        XIY(IAI,IAJ)  =   XIY(IEI,IEJ)
C       ETAX(IAI,IAJ)  =  ETAX(IEI,IEJ)
C       ETAY(IAI,IAJ)  = -ETAY(IEI,IEJ)
C         QJ(IBI,IBJ)  =    QJ(IDI,IDJ)
C        XIX(IBI,IBJ)  =  -XIX(IDI,IDJ)
C        XIY(IBI,IBJ)  =   XIY(IDI,IDJ)
C       ETAX(IBI,IBJ)  =  ETAX(IDI,IDJ)
C       ETAY(IBI,IBJ)  = -ETAY(IDI,IDJ)
C         QJ(IAI,IAJ)  =    QJ(ICI,ICJ)
C        XIX(IAI,IAJ)  =   XIX(ICI,ICJ)
C        XIY(IAI,IAJ)  =   XIY(ICI,ICJ)
C       ETAX(IAI,IAJ)  =  ETAX(ICI,ICJ)
C       ETAY(IAI,IAJ)  =  ETAY(ICI,ICJ)
C         QJ(IBI,IBJ)  =    QJ(ICI,ICJ)
C        XIX(IBI,IBJ)  =   XIX(ICI,ICJ)
C        XIY(IBI,IBJ)  =   XIY(ICI,ICJ)
C       ETAX(IBI,IBJ)  =  ETAX(ICI,ICJ)
C       ETAY(IBI,IBJ)  =  ETAY(ICI,ICJ)
  130 CONTINUE
      ENDIF
      IF(IBTYPE .EQ. 4) THEN
      IDB    = IDBC(12,IBC)
      ISB    = IDBC(13,IBC)
      ICB    = IDBC(14,IBC)
      I3S    = IDBC(15,IBC)
      I3E    = IDBC(16,IBC)
      IPB    = IDBC(26,IBC)
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
      NIM3=NI-3
	NIM2=NI-2
	NIM1=NI-1
	NJM3=NJ-3
	NJM2=NJ-2
	NJM1=NJ-1

      DO J=-2,NJ3
      QJ(NI3,J)=QJ(NIM3,J)
	QJ(NI2,J)=QJ(NIM2,J)
	QJ(NI1,J)=QJ(NIM1,J)
	QJ(-2,J)=QJ(4,J)
	QJ(-1,J)=QJ(3,J)
	QJ( 0,J)=QJ(2,J)
	ENDDO
      DO I=-2,NI3
      QJ(I,NJ3)=QJ(I,NJM3)
	QJ(I,NJ2)=QJ(I,NJM2)
	QJ(I,NJ1)=QJ(I,NJM1)
	QJ(I,-2)=QJ(I,4)
	QJ(I,-1)=QJ(I,3)
	QJ(I, 0)=QJ(I,2)
	ENDDO
C
      IF(IAXIS .EQ. 1) THEN
      DO 300 IJ = -2,MI*MJ-2    !有問題
      IF(Y(IJ,-1) .EQ. 0.) THEN
        QJ(IJ,-1) = QJ(IJ,-1)/1.E-25
      ELSE
        QJ(IJ,-1) = QJ(IJ,-1)/ABS(Y(IJ,-1))
      ENDIF
  300 CONTINUE
      ENDIF
C     CALL TWRITE(-1,NI2,-1,NJ2,X,'X         ')
C     CALL TWRITE(-1,NI2,-1,NJ2,Y,'Y         ')
C     CALL TWRITE(-1,NI2,-1,NJ2,QJ,'QJ        ')
C     CALL TWRITE(-1,NI2,-1,NJ2,XIX,'XIX       ')
C     CALL TWRITE(-1,NI2,-1,NJ2,XIY,'XIY       ')
C     CALL TWRITE(-1,NI2,-1,NJ2,ETAX,'ETAX      ')
C     CALL TWRITE(-1,NI2,-1,NJ2,ETAY,'ETAY      ')
C
C     CALL TWRITE(104,109,-1,3,X,'X         ')
C     CALL TWRITE(104,109,-1,3,Y,'Y         ')
C     CALL TWRITE(104,109,-1,3,QJ,'QJ        ')
C     CALL TWRITE(104,109,-1,3,XIX,'XIX       ')
C     CALL TWRITE(104,109,-1,3,XIX,'XIX       ')
C     CALL TWRITE(104,109,-1,3,ETAX,'ETAX      ')
C     CALL TWRITE(104,109,-1,3,ETAY,'ETAY      ')
C     CALL TWRITE(13,18,-1,3,X,'X         ')
C     CALL TWRITE(13,18,-1,3,Y,'Y         ')
C     CALL TWRITE(13,18,-1,3,QJ,'QJ        ')
C     CALL TWRITE(13,18,-1,3,XIX,'XIX       ')
C     CALL TWRITE(13,18,-1,3,XIX,'XIX       ')
C     CALL TWRITE(13,18,-1,3,ETAX,'ETAX      ')
C     CALL TWRITE(13,18,-1,3,ETAY,'ETAY      ')
      DO 200 IJ = 1, MI*MJ
      I          = MOD(IJ-1,MI)-2
      J          = (IJ-1)/MI-2
      IF(QJ(I,J) .LE. 0.) THEN
        WRITE(*,*)'QJ .LE. 0.  AT IJ=',IJ,I,J
C       STOP
      ENDIF
  200 CONTINUE
      RETURN
      END

      SUBROUTINE POINT
      IMPLICIT REAL*8(A-H,O-Z)
C
      COMMON /TAPE/ INAME,IGRID,IWRIT,IREAD,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXB,MAXG,MAXL,MAXV,MAXBP,MAXLUI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /NEXT / NEXTLOC
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /PNT  /  LVX,LVY,LWX,LWY,LCXI,LCET,LCXP,LCXM,LCEP,LCEM,
     &                LX,LY,LXIX,LXIY,LETAX,LETAY,LJAC,LDT,
     &                LG,LH,LD,LT,LU,LV,LP,LXM,LQX,LQY,LA1,LA2,LDD,
     &                LFG,LFH,LRG,LRH,LDQG,LDQH,LTXX,LTYY,LTXY,LCP,
     &                LSXX,LSYY,LSZZ,LSXY
      COMMON /PNTI / LKDI,LKDJ,LLKD,LKDP,MDA,MDB,NDIG
C...
      write(*,*)'nextloc=',nextloc
      IERR = 0
C
      NI1     = NI + 1
      NI2     = NI + 2
      NJ1     = NJ + 1
      NJ2     = NJ + 2
	NI3     = NI + 3
	NJ3     = NJ + 3
C
      ISTA     = 0
      IEND     = NI1
      JSTA     = 0
      JEND     = NJ1
C
      MI       = NI + 6
      MJ       = NJ + 6
      MXI      = MAX(MI,MJ)
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
      IS2D  = MI * MJ
      KS2D  = NI * NJ * NVX * NVY
      MXG   = IS2D
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
      IF(MXLI.GT.MAXLUI)THEN
       WRITE(*,*)' ERROR : MXLI.GT.MAXLUI (POINT)INSUF. LU IJ INDEX '
       WRITE(*,*)'         MXLI = ',MXLI,'    MAXLUI=',MAXLUI
       STOP
      ENDIF
C
      LVX     = NEXTLOC
      LVY     = LVX + NVX
      LWX     = LVY + NVY
      LWY     = LWX + NVX
      NEXTLOC = LWY + NVY
C
      LX      = NEXTLOC
      LY      = LX   + IS2D
      NEXTLOC = LY   + IS2D
C
      LXIX    = NEXTLOC
      LXIY    = LXIX  + IS2D
      LETAX   = LXIY  + IS2D
      LETAY   = LETAX + IS2D
      NEXTLOC = LETAY + IS2D
C
      LJAC    = NEXTLOC
      NEXTLOC = LJAC   + IS2D
C
      LG      = NEXTLOC
      LH      = LG + KS2D
      NEXTLOC = LH + KS2D
C
      LDT     = NEXTLOC
      LD      = LDT + IS2D
      LT      = LD  + IS2D
      LU      = LT  + IS2D
      LV      = LU  + IS2D
      LP      = LV  + IS2D
      LXM     = LP  + IS2D
      LQX     = LXM + IS2D
      LQY     = LQX + IS2D
      NEXTLOC = LQY + IS2D
C
      LFG     = NEXTLOC
      LFH     = LFG  + MXI
      LRG     = LFH  + MXI
      LRH     = LRG  + IS2D
      LDQG    = LRH  + IS2D
      LDQH    = LDQG + IS2D
      LA1     = LDQH + IS2D
      LA2     = LA1  + IS2D
      LDD     = LA2  + IS2D
      LCXI    = LDD  + IS2D
      LCET    = LCXI + IS2D
      LCXP    = LCET + IS2D
      LCXM    = LCXP + IS2D
      LCEP    = LCXM + IS2D
      LCEM    = LCEP + IS2D
      NEXTLOC = LCEM + IS2D
C
      LSXX    = NEXTLOC
      LSYY    = LSXX + IS2D
      LSZZ    = LSYY + IS2D
      LSXY    = LSZZ + IS2D
      NEXTLOC = LSXY + IS2D
C
      LTXX    = NEXTLOC
      LTYY    = LTXX + MXI
      LTXY    = LTYY + MXI
      LCP     = LTXY + MXI
      NEXTLOC = LCP  + MXI
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
      SUBROUTINE PUTBC(ID,G,H,GBC,HBC)
      PARAMETER (MAXBC=20)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
      COMMON /BC02 / IDA,ISA,ICA,IAS,IAE,IPA,LENA,
     1               IDB,ISB,ICB,IBS,IBE,IPB,LENB
      DIMENSION   ID(LENA,10)
      DIMENSION    G(NI,NJ,NVX,NVY),   H(NI,NJ,NVX,NVY)
      DIMENSION    GBC(LENA,2,NVX,NVY),HBC(LENA,2,NVX,NVY)
C
      DO 100 II = 1 , LENA
      I4        = ID(II, 7)
      J4        = ID(II, 8)
      I5        = ID(II, 9)
      J5        = ID(II,10)
      DO 100 L  = 1 , NVY
      LL        = NVY - L + 1
      DO 100 K  = 1 , NVX
      GBC(II,1,K,LL) = G(I5,J5,K,L)
      GBC(II,2,K,LL) = G(I4,J4,K,L)
      HBC(II,1,K,LL) = H(I5,J5,K,L)
      HBC(II,2,K,LL) = H(I4,J4,K,L)
  100 CONTINUE
      RETURN
      END

      SUBROUTINE TECPLT(X,Y,D,U,V,T,P,XM,QX,QY,IQFIL2)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /FINP3/ ALPHA,XKNINF,TSTM,XMINF
      COMMON /CNST0/ DM,VISREF,VISTEM,VISFAC,BO,AVO,RGAS,TREF
      COMMON /CNST1/ TINF,XMPINF,XMUINF,XNUINF,DNINF,SSINF,UDINF,
     1               RHOINF,PINF,EINF,CINF,REINF
      COMMON /CNST2/ D0,T0,P0,U0,V0,TAU0,GAMMA,PI,SPI
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      DIMENSION    X(-2:NI3,-2:NJ3),   Y(-2:NI3,-2:NJ3)
      DIMENSION   QX(-2:NI3,-2:NJ3),  QY(-2:NI3,-2:NJ3)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION    P(-2:NI3,-2:NJ3),  XM(-2:NI3,-2:NJ3)
C
      DO 10 I = 1, NI
      DO 10 J = 1, NJ
      SS  = SQRT(GAMMA*RGAS*T(I,J)*TREF)
      XM(I,J) = SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))*CINF/SS
   10 CONTINUE
C
      REWIND (IQFIL2)
      WRITE(IQFIL2,999)
  999 FORMAT('VARIABLES="X","Y","D","M","P","T","U","V"')
      WRITE(IQFIL2,998) NI,NJ
  998 FORMAT('ZONE  I=',I4,'    J=',I4,'   F=BLOCK')
      WRITE(IQFIL2,*)(( X(I,J),I=1,NI),J=1,NJ),
     1               (( Y(I,J),I=1,NI),J=1,NJ),
     1               (( D(I,J),I=1,NI),J=1,NJ),
     1               ((XM(I,J),I=1,NI),J=1,NJ),
     1               (( P(I,J),I=1,NI),J=1,NJ),
     1               (( T(I,J),I=1,NI),J=1,NJ),
     1               (( U(I,J),I=1,NI),J=1,NJ),
     1               (( V(I,J),I=1,NI),J=1,NJ)
      RETURN
      END
      SUBROUTINE QPLT3D(D,U,V,T,P,XM,QX,QY,IQFIL1)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /FINP3/ ALPHA,XKNINF,TSTM,XMINF
      COMMON /CNST0/ DM,VISREF,VISTEM,VISFAC,BO,AVO,RGAS,TREF
      COMMON /CNST1/ TINF,XMPINF,XMUINF,XNUINF,DNINF,SSINF,UDINF,
     1               RHOINF,PINF,EINF,CINF,REINF
      COMMON /CNST2/ D0,T0,P0,U0,V0,TAU0,GAMMA,PI,SPI
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      DIMENSION   QX(-2:NI3,-2:NJ3),  QY(-2:NI3,-2:NJ3)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION    P(-2:NI3,-2:NJ3),  XM(-2:NI3,-2:NJ3)
C
      DO 10 I = 1, NI
      DO 10 J = 1, NJ
      SS  = SQRT(GAMMA*RGAS*T(I,J)*TREF)
      XM(I,J) = SQRT(U(I,J)*U(I,J)+V(I,J)*V(I,J))*CINF/SS
   10 CONTINUE
C
      REWIND (IQFIL1)
      WRITE(IQFIL1,*) NI,NJ
      WRITE(IQFIL1,*) XMINF,REINF,XKNINF,TIME,ITER
      WRITE(IQFIL1,*)(( D(I,J),I=1,NI),J=1,NJ),
     1               ((XM(I,J),I=1,NI),J=1,NJ),
     1               (( P(I,J),I=1,NI),J=1,NJ),
     1               (( T(I,J),I=1,NI),J=1,NJ),
     1               (( U(I,J),I=1,NI),J=1,NJ),
     1               (( V(I,J),I=1,NI),J=1,NJ)
C      CLOSE (IQFIL1)
      RETURN
      END
      SUBROUTINE QRST(D,T,U,V,P,XM,IQFIL1,ITER0)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /CNST1/ TINF,XMPINF,XMUINF,XNUINF,DNINF,SSINF,UDINF,
     1               RHOINF,PINF,EINF,CINF,REINF
      COMMON /FINP3/ ALPHA,XKNINF,TSTM,XMINF
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION    P(-2:NI3,-2:NJ3),  XM(-2:NI3,-2:NJ3)
C
      READ(IQFIL1,*)NI0,NJ0
      IF(NI0 .NE. NI .OR. NJ0 .NE. NJ)THEN
        WRITE(*,*)'Restart ERROR : IRSTA=2, NI0.ne.NI or NJ0 .ne. NJ'
        STOP
      ENDIF
      READ(IQFIL1,*)XMINF,REINF,XKNINF,TIME,ITER0
      READ(IQFIL1,*)(( D(I,J),I=1,NI),J=1,NJ),
     1               ((XM(I,J),I=1,NI),J=1,NJ),
     1               (( P(I,J),I=1,NI),J=1,NJ),
     1               (( T(I,J),I=1,NI),J=1,NJ),
     1               (( U(I,J),I=1,NI),J=1,NJ),
     1               (( V(I,J),I=1,NI),J=1,NJ)
      RETURN
      END


      SUBROUTINE RSTGH(A,B,D,U,V,T,P,QX,QY,G,H)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /CNST2/ D0,T0,P0,U0,V0,TAU0,GAMMA,PI,SPI
      DIMENSION   QX(-2:NI3,-2:NJ3),  QY(-2:NI3,-2:NJ3)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION    P(-2:NI3,-2:NJ3)
      DIMENSION    G(NI,NJ,NVX,NVY),   H(NI,NJ,NVX,NVY)
      DIMENSION    A(NVX),B(NVY)
C
      DO 10 K = 1,NVX
      DO 10 L = 1,NVY
      DO 10 I = 1,NI
      DO 10 J = 1,NJ
      PP   = -1./T(I,J)
     1       *((A(K)-U(I,J))*(A(K)-U(I,J))+(B(L)-V(I,J))*(B(L)-V(I,J)))
      G(I,J,K,L) = D(I,J)/(PI*T(I,J))*EXP(PP)
      H(I,J,K,L) = 0.5 * T(I,J) * G(I,J,K,L)
   10 CONTINUE
C
      RETURN
      END

      SUBROUTINE SETUP(IJP)
      PARAMETER (MAXBC=20)
      IMPLICIT REAL*8(A-H,O-Z)
C..
      CHARACTER*80 BCNAME
      COMMON /TAPE / INAME,IGRID,IWRIT,IREAD,IFORC,ISHOW,IRSTA,IRESD
      COMMON /MAX0 / MAXW,MAXI,MAXB,MAXG,MAXL,MAXV,MAXBP,MAXLUI
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /NEXT / NEXTLOC
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /FINP0/ IQX,IQY,IOX,IOY,MVX,MVY,DVX,DVY,VXS,VYS
      COMMON /FINP1/ KMOD,MMOD,ZETA,PR,VISS,XLDA
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR
      COMMON /FINP3/ ALPHA,XKNINF,TSTM,XMINF
      COMMON /CNST0/ DM,VISREF,VISTEM,VISFAC,BO,AVO,RGAS,TREF
      COMMON /CNST1/ TINF,XMPINF,XMUINF,XNUINF,DNINF,SSINF,UDINF,
     1               RHOINF,PINF,EINF,CINF,REINF
      COMMON /CNST2/ D0,T0,P0,U0,V0,TAU0,GAMMA,PI,SPI
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST0/ URGHT,VRGHT,DRGHT,TRGHT,PRGHT,
     1               ULEFT,VLEFT,DLEFT,TLEFT,PLEFT
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /IMPL0/ NSOU
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
      COMMON /AXIS / IAXIS
      COMMON /ENTR / ENTEPS
      COMMON /EQUL / ISOU,IEQUL
C..
      DIMENSION KAS(2),KAE(2),KBS(2),KBE(2),IDUM(2)
      DIMENSION IJP(MAXBP)
C
      INTEGER CY(2),DD(2,2)
      DATA  CY / 2 , 1 /
      DATA  DD / 1 , 0 , 0 , 1 /
C...
      NAMELIST/VINP/IQX,IQY,IOX,IOY,MVX,MVY,DVX,DVY,VXS,VYS
      NAMELIST/MINP/KMOD,MMOD,ZETA,XLDA,GAMMA
      NAMELIST/PINP/IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILOCDT,
     1              CFL1,CFL2,NUP,DTFIX,LSTD,NSOU,IAXIS,ENTEPS,
     1              ISOU,IEQUL
      NAMELIST/FINP/XMINF,XKNINF,ALPHA,TSTM
      NAMELIST/SINP/XMS,XST,MTIM,TPP
C
      NAMELIST/BCIN/BCNAME,IBTYPE,KAS,KAE,KBS,KBE,
     1              UFAC,TFAC,TWALL
C...
      DATA  IQX,IQY,IOX,IOY,MVX,MVY,DVX,DVY,VXS,VYS
     1     /  1,  1, 4 , 4, 20, 20,  0.5,  0.5,0.,0./
      DATA  KMOD,MMOD,ZETA,XLDA,GAMMA
     1     /   1,   2,  5.,  0.,1.667/
      DATA  IRSTRT,NSTEP,IPRINT,MTHL,MTHR,ILOCDT, CFL1, CFL2,
     1         NUP,DTFIX,  LSTD,  NSOU,IAXIS,ENTEPS,ISOU,IEQUL
     1     /     0,   10,    10,   0,   1,     0,  0.9,  0.9,
     1         100, -0.1,     1,      2  , 0  , 0.05,1,0/
      DATA  XMINF,XKNINF,ALPHA,TSTM
     1     /  1.8,   0.1,   0.,  1./
      DATA    XMS,   XKNR, XST, MTIM, TPP
     1     /   1.,    1.,    1,   70,  70*0./
      DATA  BCNAME,IBTYPE,      KAS,      KAE,    KBS,    KBE,
     1      UFAC,TFAC,TWALL
     1     /'NULL',     1,   1,   1,    1,  1,  1,  1,  1,  1,
     1        1.,  1.,   1./
C...
C
C     READ(INAME,NML=VINP)              ! work in SGI & CONVEX
C     READ(INAME,NML=MINP)              ! work in SGI & CONVEX
C     READ(INAME,NML=PINP)              ! work in SGI & CONVEX
C     READ(INAME,NML=FINP)              ! work in SGI & CONVEX
C
      READ(INAME,VINP)                  ! work in CRAY
      READ(INAME,MINP)                  ! work in CRAY
      READ(INAME,PINP)                  ! work in CRAY
      READ(INAME,FINP)                  ! work in CRAY
      READ(INAME,SINP)                  ! work in CRAY
C...
C     ONE MOLE MOLECULAR AT STP CONDITION 2.68699E25
C
C     GAS PROPERTIES AT STP ( 101325Pa 0 DEG C)
C
C     GAS   MOLE(G)   GAM   TREF   MU0(NS/M/M)  W
C
C     H2    2.016    1.41   273  8.4110E-6  0.680
C     HE    4.003    1.67   293  1.9604E-5  0.657
C     N2    28.01    1.40   273  1.6630E-5  0.670
C     AIR   28.97    1.40   273  1.716PE-5  0.666
C     ARGON 39.994   1.67   273  2.1250E-5  0.720
C
C     SAMPLE GAS : ARGON
C

      DMM    = 39.994
      VISREF = 2.1250E-5
      VISTEM = 273.
      VISFAC = 0.720
C
C     SETTING GAS CONSTANTS
C
      BO     = 1.380622E-23
      AVO    = 6.0225E23
      DM     = DMM/AVO*1.E-3
      RGAS   = BO/DM
      PI     = 3.1415926535
      SPI    = SQRT(3.1415926535)
C
C     SETTING COLLISION AND MOLECULAR PARAMETER
C
      IF(KMOD.EQ.1)THEN
      PR   = 1.
      ELSE IF(KMOD.EQ.2)THEN
      PR   = 2./3.
      ELSE IF(KMOD.EQ.3)THEN
      IF(XLDA.EQ.1)THEN
       WRITE(*,*)'IN ELLIPSOIDAL MODEL, THE PARRMETER XLDA MUST NE 1'
       WRITE(*,*)'IN GENERAL -0.5 <= XLDA <= 0.   '
       STOP
      ENDIF
      PR   = 1./(1.-XLDA)
      ELSE
      WRITE(*,*)'KINETIC MODEL BEYOND SETTING, KMOD=',KMOD
      STOP
      ENDIF
C
      IF(MMOD.EQ.1)THEN
        VISS = VISFAC
      ELSEIF(MMOD.EQ.2) THEN
        VISS = 1. - (ZETA -5.)/(2.*(ZETA - 1.))
      ELSE
      WRITE(*,*)'MOLECULAR MODEL BEYOND SETTING, MMOD=',MMOD
      STOP
      ENDIF
C
C     CAL. FREE STREAM PARAMETERS (DIMENSION)
C
      TREF   = 273.
      CL     = 1.
      TINF   = TSTM*TREF
      XMPINF = XKNINF * CL
      XMUINF = VISREF*(TINF/VISTEM)**VISS
      DNINF  = 16./5.*XMUINF/(DM*XMPINF*SQRT(2.*PI*RGAS*TINF))
      XNUINF = DNINF*BO*TINF/XMUINF
      SSINF  = SQRT(GAMMA*RGAS*TINF)
      UDINF  = XMINF * SSINF
      RHOINF = DM*DNINF
      PINF   = RHOINF*RGAS*TINF
      EINF   = 1.5*RGAS*TINF
      CINF   = SQRT(2.*RGAS*TINF)
C
C     CAL. NONDIMENSION PARAMETERS
C
      VXS   = VXS*SSINF/CINF
      VYS   = VYS*SSINF/CINF
C
      REINF  = RHOINF*UDINF*CL/XMUINF
      TAU0   = XMPINF/(2.*CINF/SPI)/(CL/CINF)
C
      D0   = 1.
      T0   = TINF/TREF
      P0   = D0*T0
      U0   = UDINF*COS(ALPHA/57.3)/CINF
      V0   = UDINF*SIN(ALPHA/57.3)/CINF
C
      WRITE(*,99)
      WRITE(*,*) '    FREE STREAM STATE'
      WRITE(*,*) 'MOLECULAR MASS              =  ',DM
      WRITE(*,*) 'DENSITY (KG/M**3)           =  ',RHOINF
      WRITE(*,*) 'NUMBER DENSITY              =  ',DNINF
      WRITE(*,*) 'PRESSURE (N/M**2)           =  ',PINF
      WRITE(*,*) 'TEMPERATURE (K)             =  ',TINF
      WRITE(*,*) 'INTERNAL ENERGY (M/S)**2    =  ',EINF
      WRITE(*,*) 'SONIC SPEED  (M/S)          =  ',SSINF
      WRITE(*,*) 'VELOCITY  C0                =  ',CINF
      WRITE(*,*) 'VELOCITY  UD/C0             =  ',UDINF/CINF
      WRITE(*,*) 'VISCOUS COEFF               =  ',XMUINF
      WRITE(*,*) 'MEAN FREE PATH              =  ',XMPINF
      WRITE(*,*) 'COLLISION FREQUENCE         =  ',XNUINF
      WRITE(*,*) 'COLLISION TIME              =  ',TAU0
      WRITE(*,*) 'KNUDSEN NO             KN   =  ',XKNINF
      WRITE(*,*) 'REYNOLD NO             RE   =  ',REINF
      WRITE(*,*) 'MACH NO.                    =  ',XMINF
C
      IF(LSTD.EQ.0)THEN
C
C     INITIAL PARAMETERS FOR UNSTEADY PROBLEM (SHOCK REFLECTION PROBLEM)
C
        URGHT    = 0.
        VRGHT    = 0.
        DRGHT    = 1.
        TRGHT    = 1.
        PRGHT    = 1.
        XMPRD    = XKNR * CL
        XMURD    = VISREF*(TRGHT*TREF/VISTEM)**VISS
        DRD      = 16./5.*XMURD/(DM*XMPRD*SQRT(2.*PI*RGAS*TRGHT*TREF))
        XNURD    = DRD*BO*TINF/XMUINF
        SSRD     = SQRT(GAMMA*RGAS*TRGHT*TREF)
        CRD      = SQRT(2.*RGAS*TRGHT*TREF)
        RER      = 0.
        TAUR     = XMPRD/(2.*CRD/SPI)/(CL/CRD)
C
        ULEFT    = GAMMA*(XMS*XMS-1.)/((GAMMA+1.)*XMS)
        VLEFT    = 0.
        PLEFT    = PRGHT*(1.+2.*GAMMA/(GAMMA+1.)*(XMS*XMS-1.))
        DLEFT    = DRGHT*(GAMMA+1.)*XMS*XMS/(2.+(GAMMA-1.)*XMS*XMS)
        TLEFT    = PLEFT/DLEFT
        TLD      = TLEFT*TREF
        DLD      = DLEFT*DRD
        XMULD    = VISREF*(TLD/VISTEM)**VISS
        XMPLD    = 16./5.*XMULD/(DM*DLD*SQRT(2.*PI*RGAS*TLD))
        XNULD    = DLD*BO*TLD/XMULD
        SSLD     = SQRT(GAMMA*RGAS*TLEFT*TREF)
        CLD      = SQRT(2.*RGAS*TLEFT*TREF)
        REL      = DM*DLD*ULEFT*CRD*CL/XMULD
        TAUL     = XMPLD/(2.*CRD/SPI)/(CL/CRD)
        XKNL     = XMPLD/CL
C
        CINF     = CRD
        SSINF    = SSRD
        TAU0     = TAUR
        XKNINF   = XKNR
        TINF     = TREF*TRGHT
        XMPINF   = XMPRD
        XMUINF   = XMURD
        XNUINF   = XNURD
        DNINF    = DRD
        UDINF    = 0.
        RHOINF   = DM*DNINF
        PINF     = RHOINF*RGAS*TINF
        EINF     = 1.5*RGAS*TINF
        CINF     = SQRT(2.*RGAS*TINF)
C
      WRITE(*,99)
      WRITE(*,*) 'SHOCK MACH NO.              =  ',XMS
      WRITE(*,*) 'SHOCK POSITION              =  ',XST
      WRITE(*,*) '    LEFT SIDE STATE      '
      WRITE(*,*) 'NUMBER DENSITY              =  ',DLEFT
      WRITE(*,*) 'VELOCITY  UD/C0             =  ',ULEFT
      WRITE(*,*) 'MACH NO.                    =  ',ULEFT*CRD/SSRD
      WRITE(*,*) 'TEMPERATURE                 =  ',TLEFT
      WRITE(*,*) 'PRESSURE                    =  ',PLEFT
      WRITE(*,*) 'VISCOUS COEFF               =  ',XMULD
      WRITE(*,*) 'MEAN FREE PATH              =  ',XMPLD
      WRITE(*,*) 'COLLISION FREQUENCE         =  ',XNULD
      WRITE(*,*) 'COLLISION TIME              =  ',TAUL
      WRITE(*,*) 'KNUDSEN NO             KN   =  ',XKNL
      WRITE(*,*) 'REYNOLD NO             RE   =  ',REL
      WRITE(*,*) '    RIGHT SIDE STATE      '
      WRITE(*,*) 'NUMBER DENSITY              =  ',DRGHT
      WRITE(*,*) 'VELOCITY  UD/C0             =  ',URGHT
      WRITE(*,*) 'MACH NO.                    =  ',URGHT*CRD/SSRD
      WRITE(*,*) 'TEMPERATURE                 =  ',TRGHT
      WRITE(*,*) 'PRESSURE                    =  ',PRGHT
      WRITE(*,*) 'VISCOUS COEFF               =  ',XMURD
      WRITE(*,*) 'MEAN FREE PATH              =  ',XMPRD
      WRITE(*,*) 'COLLISION FREQUENCE         =  ',XNURD
      WRITE(*,*) 'COLLISION TIME              =  ',TAUR
      WRITE(*,*) 'KNUDSEN NO             KN   =  ',XKNR
      WRITE(*,*) 'REYNOLD NO             RE   =  ',RER
      ENDIF
C
   99 FORMAT(80('-'))
C...
      NVX  = MVX
      NVY  = MVY
      MXV  = NVX*NVY
C
      IF(MXV .GT. MAXV)THEN
        WRITE(*,*)' ERROR : MXV > MAXV (SETUP) INSUFF. VEL POINTS'
        WRITE(*,*)'         MXV = ',MXV,'     MAXV = ',MAXV
        STOP
      ENDIF
C
      NDIM(1)  = NI
      NDIM(2)  = NJ
C
      NBC      = 0
      NBC4     = 0
      NCOUNT   = 0
      NEXTPNT  = 1
C
 1000 CONTINUE
C..
C     IBTYPE = 0     FARFIELD
C            = 1     WALL
C            = 2     SYMMETRY I WITH Y=0
C            = 3     SYMMETRY J WITH Y=0
C            = 4     SYMMETRY J WITH Y=0 (FOR C GRID WAKE MATCH)
C            = 5     ZERO GRADIENT
C
      NCOUNT = NCOUNT + 1
C
C     READ(INAME,NML=BCIN)               ! work in SGI & CONVEX
C
      READ(INAME,BCIN)                   ! work in CRAY
C..
      IF( IBTYPE .LE. -1 .OR. BCNAME .EQ. 'QUIT' ) GO TO 2000
      WRITE(*,1300) NCOUNT,BCNAME
      WRITE(*,1400) IBTYPE,KAS,KAE
      IF( IBTYPE.EQ.1 ) WRITE(6,1500) TWALL
      IF( IBTYPE.EQ.4 ) WRITE(6,1600) KBS,KBE
C
      IF( IBTYPE .GT. 5 ) THEN
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
C
C     ADDITION WALL BC PARAMETER
C
      IF(IBTYPE .EQ. 1) THEN
        FDBC( 1,NBC) = TWALL
        FDBC( 2,NBC) = UFAC
        FDBC( 3,NBC) = TFAC
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
      WRITE(69,*)'NBC = ',NBC
      WRITE(69,*)'IBTYPE = ',IBTYPE
      WRITE(69,*)'IDA = ',IDA
      WRITE(69,*)'ISA = ',ISA
      WRITE(69,*)'ICA = ',ICA
      WRITE(69,*)'I2S = ',I2S
      WRITE(69,*)'I2E = ',I2E
      WRITE(69,*)'LENG = ',LENG
      WRITE(69,*)'IP = ',IP
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
         WRITE(69,71)IAI,IAJ,IBI,IBJ,ICI,ICJ,IDI,IDJ,IEI,IEJ
   71    FORMAT(1X,10I5)
  200 CONTINUE
C
      IDBC(25,NBC) = IP
      NEXTPNT      = NEXTPNT + LENG*10
C
      IF(IBTYPE .EQ. 2 .OR. IBTYPE .EQ. 3) THEN
      NPBC         = ABS(I2E-I2S) + 1
      IDBC(7,NBC)  = NEXTLOC
      IDBC(8,NBC)  = IDBC(7,NBC)  + NPBC * 2 * NVX * NVY
      NEXTLOC      = IDBC(8,NBC)  + NPBC * 2 * NVX * NVY
      ENDIF
      WRITE(*,*)'NEXTLOC for SYMM. BC ',NEXTLOC
C
C     FOR WAKE BC
C
      IF(IBTYPE .EQ. 4) THEN
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
      WRITE(69,*)'NBC = ',NBC
      WRITE(69,*)'IBTYPE = ',IBTYPE
      WRITE(69,*)'IDA = ',IDA,'   IDB = ',IDB
      WRITE(69,*)'ISA = ',ISA,'   ISB = ',ISB
      WRITE(69,*)'ICA = ',ICA,'   ICB = ',ICB
      WRITE(69,*)'I2S = ',I2S,'   I3S = ',I3S
      WRITE(69,*)'I2E = ',I2E,'   I3E = ',I3E
      WRITE(69,*)'LENG = ',LENG,'   LENB = ',LENB
      WRITE(69,*)'IP = ',IP,'   IPB = ',IPB
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
         WRITE(69,71)IAI,IAJ,IBI,IBJ,ICI,ICJ,IDI,IDJ,IEI,IEJ
  210 CONTINUE
C
      IDBC(26,NBC) = IPB
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
 1300 FORMAT( /2X,'SURFACE NO. ',I4,2X,':',2X,A60 )
 1400 FORMAT(  5X,'BCTYPE = ',I5,
     1        /5X,'IBS    = ',I5,2X,'    JBS=',I5
     1        /5X,'IBE    = ',I5,2X,'    JBE=',I5)
 1500 FORMAT(  5X,'WALL TEMP. =',F5.2)
 1600 FORMAT(  5X,'IPS    = ',I5,2X,'    JPS=',I5
     1        /5X,'IPE    = ',I5,2X,'    JPE=',I5)
      RETURN
      END

      SUBROUTINE SOURCE(K,L,DTLOC,CX,CY,QJ,D,T,U,V,P,QX,QY,
     1                  RG,RH,DQG,DQH,RR,WW,A1,A2,SLXX,SLYY,SLZZ,SLXY)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /FINP0/ IQX,IQY,IOX,IOY,MVX,MVY,DVX,DVY,VXS,VYS
      COMMON /FINP1/ KMOD,MMOD,ZETA,PR,VISS,XLDA
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR
      COMMON /FINP3/ ALPHA,XKNINF,TSTM,XMINF
      COMMON /CNST0/ DM,VISREF,VISTEM,VISFAC,BO,AVO,RGAS,TREF
      COMMON /CNST1/ TINF,XMPINF,XMUINF,XNUINF,DNINF,SSINF,UDINF,
     1               RHOINF,PINF,EINF,CINF,REINF
      COMMON /CNST2/ D0,T0,P0,U0,V0,TAU0,GAMMA,PI,SPI
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      DIMENSION   QX(-2:NI3,-2:NJ3),  QY(-2:NI3,-2:NJ3)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION    P(-2:NI3,-2:NJ3)
      DIMENSION   RG(-2:NI3,-2:NJ3),  RH(-2:NI3,-2:NJ3)
      DIMENSION  DQG(-2:NI3,-2:NJ3), DQH(-2:NI3,-2:NJ3)
      DIMENSION   QJ(-2:NI3,-2:NJ3),DTLOC(-2:NI3,-2:NJ3)
      DIMENSION   A1(-2:NI3,-2:NJ3),  A2(-2:NI3,-2:NJ3)
      DIMENSION SLXX(-2:NI3,-2:NJ3),SLYY(-2:NI3,-2:NJ3)
      DIMENSION SLZZ(-2:NI3,-2:NJ3),SLXY(-2:NI3,-2:NJ3)
C
      NIJ      = NI *NJ
      MIJ      = MI *MJ
      MI2      = MI - 2
      MJ2      = MJ - 2
      IF(KMOD.EQ.1)THEN
C     DO 10 IJ   = 1,NIJ
C      I         = MOD(IJ-1,NI) + 1
C      J         = (IJ-1)/NI + 1
C
C     AP = -((CX-U(I,J))**2 + (CY-V(I,J))**2)/T(I,J)
C     BG = A1(I,J)*EXP(AP)
C     DQG(I,J) =    DQG(I,J)
C    1            + DTLOC(I,J)*(           BG-RG(I,J))*A2(I,J)
C     DQH(I,J) =    DQH(I,J)
C    1            + DTLOC(I,J)*(0.5*T(I,J)*BG-RH(I,J))*A2(I,J)
C  10 CONTINUE
      DO 10 IJ = -2,MIJ-2
C
      AP = -((CX-U(IJ,-2))**2 + (CY-V(IJ,-2))**2)/T(IJ,-2)
      BG = A1(IJ,-2)*EXP(AP)
      DQG(IJ,-2) =    DQG(IJ,-2)
     1            + DTLOC(IJ,-2)*(             BG-RG(IJ,-2))*A2(IJ,-2)
      DQH(IJ,-2) =    DQH(IJ,-2)
     1            + DTLOC(IJ,-2)*(0.5*T(IJ,-2)*BG-RH(IJ,-2))*A2(IJ,-2)
   10 CONTINUE
C

C
      ENDIF
C
      IF(KMOD.EQ.2)THEN
C
      OPR = 1. - PR
C
C     DO 20 IJ   = 1,MIJ
C      I         = MOD(IJ-3,MI2) + 1 
C      J         = (IJ-3)/MI2 + 1 
C
C     if(t(i,j).le.0.)then
C       write(*,*)'t',i,j,t(i,j)
C     endif
C     if(p(i,j).le.0.)then
C       write(*,*)'p',i,j,p(i,j)
C     endif
C     U1    = CX-U(I,J)
C     V1    = CY-V(I,J)
C     AP    = -(U1**2 + V1**2)/T(I,J)
C     G1    = A1(I,J)*EXP(AP)
C     C5    = 2./(5.*P(I,J))/T(I,J)
C     QC    = U1*QX(I,J)+V1*QY(I,J)
C     G2    = 1. + OPR*QC*(-2.*AP - 4.)*C5
C     H2    = 1. + OPR*QC*(-2.*AP - 2.)*C5
C     BG    = G1*G2
C     BH    = 0.5*T(I,J)*G1*H2
C
C     DQG(I,J) =   DQG(I,J)
C    1            +DTLOC(I,J)*(             BG-RG(I,J))*A2(I,J)
C     DQH(I,J) =   DQH(I,J)
C    1            +DTLOC(I,J)*(0.5*T(I,J)*BG-RH(I,J))*A2(I,J)
C
C  20 CONTINUE
      DO 20 IJ   = -2,MIJ-2
       I         = MOD(IJ+1,MI) - 1
       J         = (IJ+1)/MI - 1
      IF(I.GT.1.AND.I.LT.NI.AND.J.GT.1.AND.J.LT.NJ)THEN
      U1    = CX-U(IJ,-2)
      V1    = CY-V(IJ,-2)
      AP    = -(U1**2 + V1**2)/T(IJ,-2)
      G1    = A1(IJ,-2)*EXP(AP)
      C5    = 2./(5.*P(IJ,-2))/T(IJ,-2)
      QC    = U1*QX(IJ,-2)+V1*QY(IJ,-2)
      G2    = 1. + OPR*QC*(-2.*AP - 4.)*C5
      H2    = 1. + OPR*QC*(-2.*AP - 2.)*C5
      BG    = G1*G2
      BH    = 0.5*T(IJ,-2)*G1*H2
C
      DQG(IJ,-2) =   DQG(IJ,-2)
     1            +DTLOC(IJ,-2)*(             BG-RG(IJ,-2))*A2(IJ,-2)
      DQH(IJ,-2) =   DQH(IJ,-2)
     1            +DTLOC(IJ,-2)*(0.5*T(IJ,-2)*BG-RH(IJ,-2))*A2(IJ,-2)
      ELSE
      AP = -((CX-U(IJ,-2))**2 + (CY-V(IJ,-2))**2)/T(IJ,-2)
      BG = A1(IJ,-2)*EXP(AP)
      DQG(IJ,-2) =    DQG(IJ,-2)
     1            + DTLOC(IJ,-2)*(             BG-RG(IJ,-2))*A2(IJ,-2)
      DQH(IJ,-2) =    DQH(IJ,-2)
     1            + DTLOC(IJ,-2)*(0.5*T(IJ,-2)*BG-RH(IJ,-2))*A2(IJ,-2)
      ENDIF
C
   20 CONTINUE
C
      ENDIF
C
      IF(KMOD.EQ.3)THEN
C
C     DO 30 IJ   = 1,MIJ
C      I         = MOD(IJ-1,MI) - 1
C      J         = (IJ-1)/MI - 1
C
C     XLXX      = 0.5*(1.-XLDA)*T(I,J) -
C    1            XLDA*U(I,J)*U(I,J)+XLDA*SLXX(I,J)/D(I,J)
C     XLYY      = 0.5*(1.-XLDA)*T(I,J) -
C    1            XLDA*V(I,J)*V(I,J)+XLDA*SLYY(I,J)/D(I,J)
C     XLXY      = -XLDA*U(I,J)*V(I,J) +XLDA*SLXY(I,J)/D(I,J)
C
C     XLYY      = 0.5*(1.-XLDA)*T(I,J)+XLDA*SLYY(I,J)/D(I,J)
C     XLXY      = XLDA*SLXY(I,J)/D(I,J)
C
C     XLZZ      = 0.5*(1.-XLDA)*T(I,J)+XLDA*SLZZ(I,J)/D(I,J)
C     XA1       = ABS(XLXX*XLYY-XLXY*XLXY)
C     IF(XA1.LT.0.)THEN
C      write(*,*) i,j,k,l,xa1,xlxx,xlyy,xlxy
C     endif
C     XA2       = A1(I,J)*0.5*T(I,J)/SQRT(XA1)
C     XA3       = -0.5*XLYY/XA1
C     XA4       = -0.5*XLXX/XA1
C     XA5       = XLXY/XA1
C        U1     = CX - U(I,J)
C        V1     = CY - V(I,J)
C        U2     = U1*U1
C        V2     = V1*V1
C        UV     = U1*V1
C        BG     = XA2*EXP(U2*XA3+V2*XA4+XA5*UV)
C        BH     = XLZZ*BG
C
C     DQG(I,J) =   DQG(I,J)
C    1            +DTLOC(I,J)*(BG-RG(I,J))*A2(I,J)
C     DQH(I,J) =   DQH(I,J)
C    1            +DTLOC(I,J)*(BH-RH(I,J))*A2(I,J)
C  30 CONTINUE
C
      DO 30 IJ   = -2,MIJ-2
       I         = MOD(IJ+1,MI) - 1
       J         = (IJ+1)/MI - 1
      IF(I.GT.1.AND.I.LT.NI.AND.J.GT.1.AND.J.LT.NJ)THEN
      XLXX      = 0.5*(1.-XLDA)*T(IJ,-2) -
     1            XLDA*U(IJ,-2)*U(IJ,-2)+XLDA*SLXX(IJ,-2)/D(IJ,-2)
      XLYY      = 0.5*(1.-XLDA)*T(IJ,-2) -
     1            XLDA*V(IJ,-2)*V(IJ,-2)+XLDA*SLYY(IJ,-2)/D(IJ,-2)
      XLXY      = -XLDA*U(IJ,-2)*V(IJ,-2) +XLDA*SLXY(IJ,-2)/D(IJ,-2)
      XLZZ      = 0.5*(1.-XLDA)*T(IJ,-2)+XLDA*SLZZ(IJ,-2)/D(IJ,-2)
      XA1       = ABS(XLXX*XLYY-XLXY*XLXY)
      XA2       = A1(IJ,-2)*0.5*T(IJ,-2)/SQRT(XA1)
      XA3       = -0.5*XLYY/XA1
      XA4       = -0.5*XLXX/XA1
      XA5       = XLXY/XA1
         U1     = CX - U(IJ,-2)
         V1     = CY - V(IJ,-2)
         U2     = U1*U1
         V2     = V1*V1
         UV     = U1*V1
         BG     = XA2*EXP(U2*XA3+V2*XA4+XA5*UV)
         BH     = XLZZ*BG
C
      DQG(IJ,-2) =   DQG(IJ,-2)
     1            +DTLOC(IJ,-2)*(BG-RG(IJ,-2))*A2(IJ,-2)
      DQH(IJ,-2) =   DQH(IJ,-2)
     1            +DTLOC(IJ,-2)*(BH-RH(IJ,-2))*A2(IJ,-2)
      ELSE
      AP = -((CX-U(IJ,-2))**2 + (CY-V(IJ,-2))**2)/T(IJ,-2)
      BG = A1(IJ,-2)*EXP(AP)
      DQG(IJ,-2) =    DQG(IJ,-2)
     1            + DTLOC(IJ,-2)*(             BG-RG(IJ,-2))*A2(IJ,-2)
      DQH(IJ,-2) =    DQH(IJ,-2)
     1            + DTLOC(IJ,-2)*(0.5*T(IJ,-2)*BG-RH(IJ,-2))*A2(IJ,-2)
      ENDIF
   30 CONTINUE
C
      ENDIF
      RETURN
      END
      SUBROUTINE STEP(DTLOC,A,B,WA,WB,X,Y,XIX,XIY,ETAX,ETAY,QJ,
     1                G,H,RG,RH,DQG,DQH,D,T,U,V,P,QX,QY,
     1                CXI,CET,CXIP,CXIM,CETP,CETM,FG,FH,SA1,SA2,DD,
     1                KDI,KDJ,LKD,KDP,IJP,S,SLXX,SLYY,SLZZ,SLXY)
      PARAMETER (MAXBC=20)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /NEXT / NEXTLOC
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /FINP0/ IQX,IQY,IOX,IOY,MVX,MVY,DVX,DVY,VXS,VYS
      COMMON /FINP1/ KMOD,MMOD,ZETA,PR,VISS,XLDA
      COMMON /FINP2/ IRSTRT,NSTEP,IPRINT,MTHL,MTHR
      COMMON /FINP3/ ALPHA,XKNINF,TSTM,XMINF
      COMMON /CNST0/ DM,VISREF,VISTEM,VISFAC,BO,AVO,RGAS,TREF
      COMMON /CNST1/ TINF,XMPINF,XMUINF,XNUINF,DNINF,SSINF,UDINF,
     1               RHOINF,PINF,EINF,CINF,REINF
      COMMON /CNST2/ D0,T0,P0,U0,V0,TAU0,GAMMA,PI,SPI
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST0/ URGHT,VRGHT,DRGHT,TRGHT,PRGHT,
     1               ULEFT,VLEFT,DLEFT,TLEFT,PLEFT
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /IMPL0/ NSOU
      COMMON /BC01 / IDBC(30,MAXBC),FDBC(10,MAXBC)
      COMMON /PNTI / LKDI,LKDJ,LLKD,LKDP,MDA,MDB,NDIG
      COMMON /RES0/ SUMD,Z1MAX,M1DI,M1DJ,Z2MAX,M2DI,M2DJ,SUMDG,SUMDH
      COMMON /EQUL / ISOU,IEQUL
      DIMENSION  IJP(MXBP)
      DIMENSION   QX(-2:NI3,-2:NJ3),  QY(-2:NI3,-2:NJ3)
      DIMENSION    U(-2:NI3,-2:NJ3),   V(-2:NI3,-2:NJ3)
      DIMENSION    D(-2:NI3,-2:NJ3),   T(-2:NI3,-2:NJ3)
      DIMENSION    P(-2:NI3,-2:NJ3)
      DIMENSION   RG(-2:NI3,-2:NJ3),  RH(-2:NI3,-2:NJ3)
      DIMENSION  DQG(-2:NI3,-2:NJ3), DQH(-2:NI3,-2:NJ3)
      DIMENSION    X(-2:NI3,-2:NJ3),   Y(-2:NI3,-2:NJ3)
      DIMENSION  XIX(-2:NI3,-2:NJ3), XIY(-2:NI3,-2:NJ3)
      DIMENSION ETAX(-2:NI3,-2:NJ3),ETAY(-2:NI3,-2:NJ3)
      DIMENSION   QJ(-2:NI3,-2:NJ3)
      DIMENSION  CET(-2:NI3,-2:NJ3),  CXI(-2:NI3,-2:NJ3)
      DIMENSION CETP(-2:NI3,-2:NJ3), CXIP(-2:NI3,-2:NJ3)
      DIMENSION CETM(-2:NI3,-2:NJ3), CXIM(-2:NI3,-2:NJ3)
      DIMENSION  SA1(-2:NI3,-2:NJ3),  SA2(-2:NI3,-2:NJ3)
      DIMENSION   DD(-2:NI3,-2:NJ3)
      DIMENSION DTLOC(-2:NI3,-2:NJ3)
      DIMENSION    FG(-2:MXI-2),FH(-2:MXI-2)
      DIMENSION  G(NI,NJ,NVX,NVY), H(NI,NJ,NVX,NVY)
      DIMENSION A(NVX),B(NVY),WA(NVX),WB(NVY)
      DIMENSION KDI(MDA),KDJ(MDA),LKD(MDB),KDP(MDB)
      DIMENSION S(1)
      DIMENSION SLXX(-2:NI3,-2:NJ3),SLYY(-2:NI3,-2:NJ3)
      DIMENSION SLZZ(-2:NI3,-2:NJ3),SLXY(-2:NI3,-2:NJ3)
	DIMENSION    DG(-2:NI3),      DH(-2:NI3)
	DIMENSION    GG(-2:NI3),      HH(-2:NI3)
	DIMENSION  DFXG(-2:NI3),    DFXH(-2:NI3)
	DIMENSION  GG2(-2:NI3,2),      HH2(-2:NI3,2)
	DIMENSION GG1(-2:NI3,4,2),   HH1(-2:NI3,4,2)
	DIMENSION   FGP(-2:NI3),     FHP(-2:NI3)
	DIMENSION   FGM(-2:NI3),     FHM(-2:NI3)
C
      RR    = 8./5./SPI/XKNINF
      WW    = 1.-VISS
      MIJ   = MI*MJ
      NIJ   = NI*NJ
      RSOU  = 1./NSOU
      EPS   = 0.05
C
      DO 2  IJ = -2,MIJ-2
      SA2(IJ,-2) = PR*D(IJ,-2)*(ABS(T(IJ,-2))**WW)*RR
      SA1(IJ,-2) = D(IJ,-2)/( PI*T(IJ,-2))/QJ(IJ,-2)
    2 CONTINUE
C
      DO 10 K = 1,NVX
      DO 10 L = 1,NVY
C
      CX        = A(K)
      CY        = B(L)
      CC0       = SQRT(CX*CX + CY*CY)
      CC1       = CC0*EPS
C
      DO 20 IJ  = -2,MIJ-2
      DQG(IJ,-2) = 0.
      DQH(IJ,-2) = 0.
      CXI(IJ,-2) = CX * XIX(IJ,-2) + CY * XIY(IJ,-2)
      CET(IJ,-2) = CX *ETAX(IJ,-2) + CY *ETAY(IJ,-2)
   20 CONTINUE
      DO 25 IJ   = 1,NIJ
       I         = MOD(IJ-1,NI) + 1
       J         = (IJ-1)/NI + 1
       RG(I,J)   = G(I,J,K,L)/QJ(I,J)
       RH(I,J)   = H(I,J,K,L)/QJ(I,J)
   25 CONTINUE
C
      CALL BCRG(IJP,X,Y,XIX,XIY,ETAX,ETAY,QJ,D,U,V,T,
     1              G,H,RG,RH,A,B,WA,WB,CX,CY,K,L,SA1,S(1))
C
       IF(ISOU .NE. 0) THEN
       CALL SOURCE(K,L,DTLOC,CX,CY,QJ,D,T,U,V,P,QX,QY,
     1      RG,RH,DQG,DQH,RR,WW,SA1,SA2,SLXX,SLYY,SLZZ,SLXY)
       ENDIF
C
       IF(MTHR .EQ. 1) THEN
       CALL FVIUPW(DTLOC,CX,CY,CXI, XIX, XIY,QJ,FG,FH,RG,RH,DQG,DQH)
       CALL FVJUPW(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,FG,FH,RG,RH,DQG,DQH)
       ENDIF
C
       IF(MTHR .EQ. 2) THEN
       CALL FVITVD2(DTLOC,CX,CY,CXI, XIX, XIY,QJ,FG,FH,RG,RH,DQG,DQH)
       CALL FVJTVD2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,FG,FH,RG,RH,DQG,DQH)
       ENDIF
C
       IF(MTHR .EQ. 3) THEN
       CALL FVIENO2(DTLOC,CX,CY,CXI, XIX, XIY,QJ,FG,FH,RG,RH,DQG,DQH)
       CALL FVJENO2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,FG,FH,RG,RH,DQG,DQH)
       ENDIF
C
       IF(MTHR .EQ. 4) THEN
       CALL FDITVD2(DTLOC,CX,CY,CXI, XIX, XIY,QJ,FG,FH,RG,RH,DQG,DQH)
       CALL FDJTVD2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,FG,FH,RG,RH,DQG,DQH)
       ENDIF
C
       IF(MTHR .EQ. 5) THEN
       CALL FDIENO2(DTLOC,CX,CY,CXI, XIX, XIY,QJ,FG,FH,RG,RH,DQG,DQH)
       CALL FDJENO2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,FG,FH,RG,RH,DQG,DQH)
       ENDIF
C
       IF(MTHR .EQ. 6) THEN
       CALL FDIENO3(DTLOC,CX,CY,CXI, XIX, XIY,QJ,FG,FH,RG,RH,DQG,DQH)
       CALL FDJENO3(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,FG,FH,RG,RH,DQG,DQH)
       ENDIF
C
       IF(MTHR .EQ. 7) THEN
       CALL FEITVD2(DTLOC,CX,CY,CXI, XIX, XIY,QJ,FG,FH,RG,RH,DQG,DQH)
       CALL FEJTVD2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,FG,FH,RG,RH,DQG,DQH)
       ENDIF
C
       IF(MTHR .EQ. 8) THEN
       CALL FHITVD2(DTLOC,CX,CY,CXI, XIX, XIY,QJ,FG,FH,RG,RH,DQG,DQH)
       CALL FHJTVD2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,FG,FH,RG,RH,DQG,DQH)
       ENDIF
C
       IF(MTHR .EQ. 9) THEN
       CALL FIWENO3(DTLOC,CX,CY,CXI, XIX, XIY,QJ,FG,FH,RG,RH,DQG,DQH)
       CALL FJWENO3(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,FG,FH,RG,RH,DQG,DQH)
       ENDIF
C
       IF(MTHR .EQ. 10) THEN
       CALL FIWENO2(DTLOC,CX,CY,CXI, XIX, XIY,QJ,FG,FH,RG,RH,DQG,DQH)
       CALL FJWENO2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,FG,FH,RG,RH,DQG,DQH)
       ENDIF
C
        IF(MTHL .EQ. 1)THEN
        DO 130 IJ  = -2,MIJ-2
        CXIP(IJ,-2)= 0.5*(CXI(IJ,-2) + ABS(CXI(IJ,-2)))
        CXIM(IJ,-2)= 0.5*(CXI(IJ,-2) - ABS(CXI(IJ,-2)))
        CETP(IJ,-2)= 0.5*(CET(IJ,-2) + ABS(CET(IJ,-2)))
        CETM(IJ,-2)= 0.5*(CET(IJ,-2) - ABS(CET(IJ,-2)))
        AAA        = 1./(1. + RSOU*DTLOC(IJ,-2)*SA2(IJ,-2))
        BBB        = AAA**NSOU
        DQG(IJ,-2) = DQG(IJ,-2) *BBB
        DQH(IJ,-2) = DQH(IJ,-2) *BBB
  130   CONTINUE
C
      CALL BCDQ(1,CXI,CET,DQG,DQH,IJP)
C
        DO 140 ND  = 5,NDIG-4
        NPNT = LKD(ND)
        LP   = KDP(ND) - 1
CDIR@ IVDEP
        DO 140 NP  = 1,NPNT
        I = KDI(LP + NP)
        J = KDJ(LP + NP)
        IF(I.GE.1.AND.J.GE.1)THEN
        FG(NP)     = 1./(1. + DTLOC(I,J)*(CXIP(I,J) + CETP(I,J)))
        DQG(I,J)  = ( DQG(I,J)
     1          + DTLOC(I,J)*(CXIP(I-1,J)*DQG(I-1,J)+
     1          CETP(I,J-1)*DQG(I,J-1)))* FG(NP)
        DQH(I,J)  = ( DQH(I,J)
     1          + DTLOC(I,J)*(CXIP(I-1,J)*DQH(I-1,J)+
     1          CETP(I,J-1)*DQH(I,J-1)))* FG(NP)
        ENDIF
  140 CONTINUE
C
      CALL BCDQ(2,CXI,CET,DQG,DQH,IJP)
C
        DO 150 ND  = NDIG-4,5,-1
        NPNT = LKD(ND)
        LP   = KDP(ND) - 1
CDIR@ IVDEP
        DO 150 NP  = 1,NPNT
        I = KDI(LP + NP)
        J = KDJ(LP + NP)
        IF(I.LE.NI.AND.J.LE.NJ)THEN
        FG(NP)     = 1./(1. - DTLOC(I,J)*(CXIM(I,J) + CETM(I,J)))
        DQG(I,J)  = ( DQG(I,J)
     1          - DTLOC(I,J)*(CXIM(I+1,J)*DQG(I+1,J)+
     1          CETM(I,J+1)*DQG(I,J+1)))* FG(NP)
        DQH(I,J)  = ( DQH(I,J)
     1          - DTLOC(I,J)*(CXIM(I+1,J)*DQH(I+1,J)+
     1          CETM(I,J+1)*DQH(I,J+1)))* FG(NP)
        ENDIF
  150   CONTINUE
      ENDIF
C
      DO 90 IJ   = 1,NIJ
       I         = MOD(IJ-1,NI) + 1
       J         = (IJ-1)/NI  + 1
      DGQ        = QJ(I,J)    * DQG(I,J)
      DHQ        = QJ(I,J)    * DQH(I,J)
      IF(I.GT.1.AND.J.GT.1.AND.I.LT.NI.AND.J.LT.NJ)THEN
      SUMDG      = SUMDG      + DGQ*DGQ
      SUMDH      = SUMDH      + DHQ*DHQ
      ENDIF
      G(I,J,K,L) = G(I,J,K,L) + DGQ
      H(I,J,K,L) = H(I,J,K,L) + DHQ
   90 CONTINUE
C
   10 CONTINUE
      CALL BCG(IJP,X,Y,XIX,XIY,ETAX,ETAY,QJ,D,U,V,T,
     1              G,H,A,B,WA,WB,SA1)
C
      RETURN
      END

      SUBROUTINE TWRITE(IIS,IIE,JJS,JJE,A,TEXT)
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER * 10 TEXT
      DIMENSION A(IIS:IIE,JJS:JJE)
      N  = IIE - IIS + 1
      KK = N/10 + 1
      WRITE(66,13)TEXT
   13 FORMAT(///,30('-'),A10,80('-'))
      DO 10 LL = 1,KK
      IS = (LL-1)*10+IIS
      IE = LL*10+IIS
      IF(IE .GT. IIE) IE = IIE
      WRITE(66,12) (I,I=IS,IE)
      DO 10 J = JJS,JJE
      WRITE(66,11) J,(A(I,J),I=IS,IE)
   10 CONTINUE
   11 FORMAT(1X,I5,10(1X,E12.3))
   12 FORMAT(/////,6X,10(I13))
      RETURN
      END
  
      SUBROUTINE VCOTES(NORD,NV,W,A,DA,AS,IS,MV)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(MV),A(MV)
C
C     NEWTON-COTES INTEGRATING
C
C     NEWTON-COTES NORD=2
C
      IF(NORD .EQ. 2) THEN
C
      NN = NV - 2
      IF(MOD(NN,2).NE.1)THEN
         WRITE(*,*)'NV ERROR FOR NEWTON-COTES NORD=2,  NV=',NV
      ENDIF
      DO 110 I = 2,NV-1
      W(I)     = 2./3.*DA
      IF(MOD(I,2).EQ.0) W(I) = 4./3.*DA
  110 CONTINUE
      W(1)     = 1./3.*DA
      W(NV)    = W(1)
C
      NVM      = NV - 1
      A0       = -0.5*NVM*DA + AS
      IF(IS .EQ. 1) A0 = AS
      DO 120 I = 1,NV
      A(I)     = A0 + (I-1.)*DA
  120 CONTINUE
C
C     NEWTON-COTES NORD=3
C
      ELSEIF(NORD .EQ. 3) THEN
C
      NN = NV - 4
      IF(MOD(NN,3).NE.0)THEN
         WRITE(*,*)'NV ERROR FOR NEWTON-COTES NORD=3,  NV=',NV
      ENDIF
C
      DO 210 I = 2,NV-1
      W(I)     = 9./8.*DA
      IF(MOD(I-1,3).EQ.0) W(I) = 3./4.*DA
  210 CONTINUE
      W(1)     = 3./8.*DA
      W(NV)    = W(1)
C
      NVM      = NV - 1
      A0       = -0.5*NVM*DA + AS
      IF(IS .EQ. 1) A0 = AS
      DO 220 I = 1,NV
      A(I)     = A0 + (I-1.)*DA
  220 CONTINUE
C
C     NEWTON-COTES NORD=3
C
      ELSEIF(NORD .EQ. 4) THEN
C
      NN = NV - 5
      IF(MOD(NN,4).NE.0)THEN
         WRITE(*,*)'NV ERROR FOR NEWTON-COTES NORD=4,  NV=',NV
      ENDIF
C
      DO 310 I = 2,NV-1
      W(I)     = 64./45.*DA
      IF(MOD(I-1,4).EQ.0) W(I) = 28./45.*DA
      IF(MOD(I-1,4).EQ.2) W(I) = 24./45.*DA
  310 CONTINUE
      W(1)     = 14./45.*DA
      W(NV)    = W(1)
C
      NVM      = NV - 1
      A0       = -0.5*NVM*DA + AS
      IF(IS .EQ. 1) A0 = AS
      DO 320 I = 1,NV
      A(I)     = A0 + (I-1.)*DA
  320 CONTINUE
C
      ELSE
       WRITE(*,*)'NCOTES ERROR : NORD NOT MATCH   NORD=',NORD
      ENDIF
C
      RETURN
      END
      SUBROUTINE VELCOE(A,B,WA,WB,MVX,MVY,MXV,NVX,NVY)
      PARAMETER (NV = 8000)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /VELC/ W0(NV),W1(NV),W2(NV),W3(NV),
     1              A3(NV),A4(NV),A5(NV),A6(NV)
      DIMENSION A(MVX),B(MVY),WA(MVX),WB(MVY)
C
      DO 100 KL = 1,MXV
      K      = MOD(KL-1,NVX)+1
      L      = (KL-1)/NVX + 1
      W0(KL)   = WA(K)*WB(L)
      W1(KL)   = W0(KL)*A(K)
      W2(KL)   = W0(KL)*B(L)
      W3(KL)   = W1(KL)*A(K) + W2(KL)*B(L)
      A3(KL)   = W1(KL)*A(K)
      A4(KL)   = W2(KL)*B(L)
      A5(KL)   = A(K)*A(K) + B(L)*B(L)
      A6(KL)   = W1(KL)*B(L)
  100 CONTINUE
      RETURN
      END

      SUBROUTINE VELDAT(IQ,NORD,NV,W,A,DA,AS,IS,MV,IVGRD)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(MV),W(MV)
      IF(IQ.EQ.1) THEN
        CALL VHERMI(NV,W,A,AS,IS,MV)
      ELSEIF(IQ.EQ.2) THEN
        CALL VLAGUE(NV,W,A,AS,IS,MV)
      ELSEIF(IQ.EQ.3) THEN
        CALL VCOTES(NORD,NV,W,A,DA,AS,IS,MV)
      ELSEIF(IQ.EQ.4) THEN
        CALL VHUANG(NV,W,A,AS,IS,MV)
      ELSEIF(IQ.EQ.5) THEN
C     READ FROM DATA FILE fort.22
      DO 12 K = 1,NV
      READ(IVGRD,*)KK,A(K),W(K)
   12 CONTINUE
      ENDIF
      WRITE(67,*)'DISCRETE NODE and WEIGHT',NV
      DO 10 I = 1,NV
      WRITE(67,11)I,A(I),W(I)
   10 CONTINUE
   11 FORMAT(I5,2E18.9)
      RETURN
      END
      SUBROUTINE VHERMI(NV,W,A,AS,IS,MV)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(MV),A(MV)
      DIMENSION A7(4),A8(4),A9(5),A10(5),A11(6),A12(6),A13(7),A14(7),
     1          A15(8),A16(8),A17(9),A18(9),A19(10),A20(10)
      DIMENSION W7(4),W8(4),W9(5),W10(5),W11(6),W12(6),W13(7),W14(7),
     1          W15(8),W16(8),W17(9),W18(9),W19(10),W20(10)
C...
      DATA (A7(I),W7(I),I=1,4)/
     1   0.000000000000000,0.8102646175568,
     1   0.816287882858965,0.4256072526101,
     1   1.673551628767471,0.5451558281913E-1,
     1   2.651961356825223,0.9717812450995E-3/
C
      DATA (A8(I),W8(I),I=1,4)/
     1    0.381186990207322,0.6611470125582,
     2    1.157193712446780,0.2078023258149,
     3    1.981656756695843,0.1707798300741E-1,
     4    2.930637420257244,0.1996040722114E-3/
C
      DATA (A9(I),W9(I),I=1,5)/
     1    0.000000000000000,0.7202352156061,
     2    0.723551018752838,0.4326515590026,
     3    1.468553289216668,0.8847452739438E-1,
     4    2.266580584531843,0.4943624275537E-2,
     5    3.190993201781528,0.3960697726326E-4/
C
      DATA (A10(I),W10(I),I=1,5)/
     1    0.342901327223705,0.6108626337353,
     2    1.036610829789514,0.2601386110823,
     3    1.756683649299882,0.3387439445548E-1,
     4    2.532731674232790,0.1343645746781E-2,
     5    3.436159118837738,0.7640432855233E-5/
C
      DATA (A11(I),W11(I),I=1,6)/
     1    0.000000000000000,0.6547592869146,
     1    0.656809566882100,0.4293597523561,
     1    1.326557084494933,0.1172278751677,
     1    2.025948015825755,0.1191139544491E-1,
     1    2.783290099781652,0.3468194663233E-3,
     1    3.668470846559583,0.1439560393714E-5/
C
      DATA (A12(I),W12(I),I=1,6)/
     1    0.314240376254359,0.5701352362625,
     2    0.947788391240164,0.4293597523561,
     3    1.597682635152605,0.5160798561588E-1,
     4    2.279507080501060,0.3905390584629E-2,
     5    3.020637025120890,0.8573687043588E-4,
     6    3.889724897869782,0.2658551684356E-6/
C
      DATA (A13(I),W13(I),I=1,7)/
     1    0.00000000000000,0.6043931879211,
     2    0.605763879171060,0.4216162968985,
     3    1.220055036590748,0.1403233206870,
     4    1.853107651601512,0.2086277529617E-1,
     5    2.519735685678238,0.1207459992719E-2,
     6    3.246608978372410,0.2043036040271E-4,
     7    4.101337596178640,0.4825731850073E-7/
C
      DATA (A14(I),W14(I),I=1,7)/
     1    0.29174551067256,0.5364059097121,
     2    0.87871378732940,0.2731056090642,
     3    1.47668273114114,0.6850553422347E-1,
     4    2.09518325850772,0.7850054726458E-2,
     5    2.74874072498540,0.3550926135519E-3,
     6    3.46265693360227,0.4716484335019E-5,
     7    4.30444857047363,0.8628591168125E-8/
C
      DATA (A15(I),W15(I),I=1,8)/
     1    0.00000000000000,0.5641003087264,
     1    0.56506958325558,0.4120286874989,
     1    1.13611558521092,0.1584889157959,
     1    1.71999257518649,0.3078003387255E-1,
     1    2.32573248617386,0.2778068842913E-2,
     1    2.96716692790560,0.1000044412325E-3,
     1    3.66995037340445,0.1059115547711E-5,
     1    4.49999070730939,0.1522475804254E-8/
C
      DATA (A16(I),W16(I),I=1,8)/
     1    0.27348104613815,0.5079294790166,
     2    0.82295144914466,0.2806474585285,
     3    1.38025853919888,0.8381004139899E-1,
     4    1.95178799091625,0.1288031153551E-1,
     5    2.54620215784748,0.9322840086242E-3,
     6    3.17699916197996,0.2711860092538E-4,
     7    3.86944790486012,0.2320980844865E-6,
     8    4.68873893930582,0.2654807474011E-9/
C
      DATA (A17(I),W17(I),I=1,9)/
     1  0.0000000000000,0.5309179376249,
     2  0.5316330013427,0.4018264694704,
     3  1.0676487257435,0.1726482976701,
     4  1.6129243142212,0.4092003414976E-1,
     5  2.1735028266666,0.5067349957628E-2,
     6  2.7577629157039,0.2986432866978E-3,
     7  3.3789320911415,0.7112289140021E-5,
     8  4.0619466758755,0.4977078981631E-7,
     9  4.8713451936744,0.4580578930799E-10/
C
      DATA (A18(I),W18(I),I=1,9)/
     1  0.2582677505191,0.4834056947255,
     1  0.7766829192674,0.2848072856700,
     1  1.3009208583896,0.9730174764132E-1,
     1  1.8355316042616,0.1864004238754E-1,
     1  2.3862990891667,0.1888522630268E-2,
     1  2.9613775055316,0.9181126867929E-4,
     1  3.5737690684863,0.1810654481093E-5,
     1  4.2481178735681,0.1046720579579E-7,
     1  5.0483640088745,0.7828199772116E-11/
C
      DATA (A19(I),W19(I),I=1,10)/
     1  0.0000000000000,0.5029748882762,
     2  0.5035201634239,0.3916089886130,
     3  1.0103683871343,0.1836327013070,
     4  1.5241706193935,0.5081038690905E-1,
     5  2.0492317098506,0.7988866777723E-2,
     6  2.5911337897945,0.6708775214072E-3,
     7  3.1578488183476,0.2720919776316E-4,
     8  3.7621873519640,0.4488243147223E-6,
     9  4.4283528066038,0.2163051009864E-8,
     1  5.2202716905375,0.1326297094499E-11/
C
       DATA (A20(I),W20(I),I=1,10)/
     1    0.2453407083009,0.4622436696006,
     1    0.7374737285454,0.2866755053628,
     1    1.2340762153953,0.1090172060200,
     1    1.7385377121166,0.2481052088746E-1,
     1    2.2549740020893,0.3243773342238E-2,
     1    2.7888060584281,0.2283386360163E-3,
     1    3.3478545673832,0.7802556478532E-5,
     1    3.9447640401156,0.1086069370769E-6,
     1    4.6036824495507,0.4399340992273E-9,
     1    5.3874808900112,0.2229393645534E-12/
C
C     GAUSS HERMITE QUADRATURE
C         IS = 1 HALF RANGE (MUST AS=0)
C         IS = 2 FULL RANGE
C
      IF(NV .EQ. 7) THEN
      DO 210 I = 2,4
      I1       = 4 + I-1
      I2       = 4 - I+1
      A(I1)    = AS + A7(I)
      A(I2)    = AS - A7(I)
      W(I1)    = W7(I)*EXP(A7(I)*A7(I))
      W(I2)    = W7(I)*EXP(A7(I)*A7(I))
  210 CONTINUE
      A(4)     = AS
      W(4)     = W7(1)*EXP(A7(1)*A7(1))
C
      ELSEIF(NV .EQ. 8) THEN
      DO 212 I = 1,4
      I1       = 4 + I
      I2       = 5 - I
      A(I1)    = AS + A8(I)
      A(I2)    = AS - A8(I)
      W(I1)    = W8(I)*EXP(A8(I)*A8(I))
      W(I2)    = W8(I)*EXP(A8(I)*A8(I))
  212 CONTINUE
C
      ELSEIF(NV .EQ. 9) THEN
      DO 214 I = 2,5
      I1       = 5 + I-1
      I2       = 5 - I+1
      A(I1)    = AS + A9(I)
      A(I2)    = AS - A9(I)
      W(I1)    = W9(I)*EXP(A9(I)*A9(I))
      W(I2)    = W9(I)*EXP(A9(I)*A9(I))
  214 CONTINUE
      A(5)     = AS
      W(5)     = W9(1)*EXP(A9(1)*A9(1))
C
      ELSEIF(NV .EQ. 10) THEN
      DO 216 I = 1,5
      I1       = 5 + I
      I2       = 6 - I
      A(I1)    = AS + A10(I)
      A(I2)    = AS - A10(I)
      W(I1)    = W10(I)*EXP(A10(I)*A10(I))
      W(I2)    = W10(I)*EXP(A10(I)*A10(I))
  216 CONTINUE
C
      ELSEIF(NV .EQ. 11) THEN
      DO 218 I = 2,6
      I1       = 6 + I-1
      I2       = 6 - I+1
      A(I1)    = AS + A11(I)
      A(I2)    = AS - A11(I)
      W(I1)    = W11(I)*EXP(A11(I)*A11(I))
      W(I2)    = W11(I)*EXP(A11(I)*A11(I))
  218 CONTINUE
      A(6)     = AS
      W(6)     = W11(1)*EXP(A11(1)*A11(1))
C
      ELSEIF(NV .EQ. 12) THEN
      DO 220 I = 1,6
      I1       = 6 + I
      I2       = 7 - I
      A(I1)    = AS + A12(I)
      A(I2)    = AS - A12(I)
      W(I1)    = W12(I)*EXP(A12(I)*A12(I))
      W(I2)    = W12(I)*EXP(A12(I)*A12(I))
  220 CONTINUE
C
      ELSEIF(NV .EQ. 13) THEN
      DO 222 I = 2,7
      I1       = 7 + I-1
      I2       = 7 - I+1
      A(I1)    = AS + A13(I)
      A(I2)    = AS - A13(I)
      W(I1)    = W13(I)*EXP(A13(I)*A13(I))
      W(I2)    = W13(I)*EXP(A13(I)*A13(I))
  222 CONTINUE
      A(7)     = AS
      W(7)     = W13(1)*EXP(A13(1)*A13(1))
C
      ELSEIF(NV .EQ. 14) THEN
      DO 224 I = 1,7
      I1       = 7 + I
      I2       = 8 - I
      A(I1)    = AS + A14(I)
      A(I2)    = AS - A14(I)
      W(I1)    = W14(I)*EXP(A14(I)*A14(I))
      W(I2)    = W14(I)*EXP(A14(I)*A14(I))
  224 CONTINUE
C
      ELSEIF(NV .EQ. 15) THEN
      DO 226 I = 2,8
      I1       = 8 + I-1
      I2       = 8 - I+1
      A(I1)    = AS + A15(I)
      A(I2)    = AS - A15(I)
      W(I1)    = W15(I)*EXP(A15(I)*A15(I))
      W(I2)    = W15(I)*EXP(A15(I)*A15(I))
  226 CONTINUE
      A(8)     = AS
      W(8)     = W15(1)*EXP(A15(1)*A15(1))
C
      ELSEIF(NV .EQ. 16) THEN
      DO 228 I = 1,8
      I1       = 8 + I
      I2       = 9 - I
      A(I1)    = AS + A16(I)
      A(I2)    = AS - A16(I)
      W(I1)    = W16(I)*EXP(A16(I)*A16(I))
      W(I2)    = W16(I)*EXP(A16(I)*A16(I))
  228 CONTINUE
C
      ELSEIF(NV .EQ. 17) THEN
      DO 230 I = 2,9
      I1       = 9 + I -1
      I2       = 9 - I +1
      A(I1)    = AS + A17(I)
      A(I2)    = AS - A17(I)
      W(I1)    = W17(I)*EXP(A17(I)*A17(I))
      W(I2)    = W17(I)*EXP(A17(I)*A17(I))
  230 CONTINUE
      A(9)     = AS
      W(9)     = W17(1)*EXP(A17(1)*A17(1))
C
      ELSEIF(NV .EQ. 18) THEN
      DO 232 I = 1,9
      I1       = 9 + I
      I2       = 10- I
      A(I1)    = AS + A18(I)
      A(I2)    = AS - A18(I)
      W(I1)    = W18(I)*EXP(A18(I)*A18(I))
      W(I2)    = W18(I)*EXP(A18(I)*A18(I))
  232 CONTINUE
C
      ELSEIF(NV .EQ. 19) THEN
      DO 234 I = 2,10
      I1       =10 + I -1
      I2       =10 - I +1
      A(I1)    = AS + A19(I)
      A(I2)    = AS - A19(I)
      W(I1)    = W19(I)*EXP(A19(I)*A19(I))
      W(I2)    = W19(I)*EXP(A19(I)*A19(I))
  234 CONTINUE
      A(10)    = AS
      W(10)    = W19(1)*EXP(A19(1)*A19(1))
C
      ELSEIF(NV .EQ. 20) THEN
      DO 236 I = 1,10
      I1       = 10+ I
      I2       = 11- I
      A(I1)    = AS + A20(I)
      A(I2)    = AS - A20(I)
      W(I1)    = W20(I)*EXP(A20(I)*A20(I))
      W(I2)    = W20(I)*EXP(A20(I)*A20(I))
  236 CONTINUE
      ELSE
       WRITE(*,*)'VHERMIT ERROR : NV NOT MATCH   NV=',NV
       STOP
      ENDIF
C
      IF(IS.EQ.1) THEN
       MODNV2 = MOD(NV,2)
       IF(MODNV2 .EQ. 1) THEN
        IV1   = (NV - 1)/2 + 1
       ELSE
        IV1   = NV/2 + 1
       ENDIF
        IV2   =  NV
        NV    = IV2 - IV1 + 1
        DO 310 I = IV1,IV2
        II    = IV1 - I + 1
        A(II) = A(I)
        W(II) = W(I)
  310   CONTINUE
      WRITE(*,*)'HALF RANGE GAUSS HERMIT QUADRATURE, NV=',NV
      ELSE
      WRITE(*,*)'FULL RANGE GAUSS HERMIT QUADRATURE, NV=',NV
      ENDIF
      RETURN
      END
      SUBROUTINE VHUANG(NV,W,A,AS,IS,MV)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(MV),A(MV)
      DIMENSION A4(2),A6(3),A8(4),A10(5),A12(6),A14(7),A16(8),A18(9),
     1          A20(10),A22(11),A24(12),A26(13),A28(14)
      DIMENSION W4(2),W6(3),W8(4),W10(5),W12(6),W14(7),W16(8),W18(9),
     1          W20(10),W22(11),W24(12),W26(13),W28(14)
C...
      DATA (A4(I),W4(I),I=1, 2)/
     1        .3001939310608E+00 , .6405291796844E+00,
     1        .1252421045334E+01 , .2456977457684E+00/
C
      DATA (A6(I),W6(I),I=1, 3)/
     1        .1905541497982E+00 , .4460297704667E+00,
     1        .8482518675446E+00 , .3964682669983E+00,
     1        .1799776578416E+01 , .4372888798776E-01/
C 
      DATA (A8(I),W8(I),I=1, 4)/
     1        .1337764469963E+00 , .3253029997573E+00,
     1        .6243246901877E+00 , .4211071018520E+00,
     1        .1342537825646E+01 , .1334425003573E+00,
     1        .2262664477011E+01 , .6374323486246E-02/
C
      DATA (A14(I),W14(I),I=1, 7)/
     1        .6371648483377E-01 , .1606099656840E+00,
     1        .3181920198031E+00 , .3063198086382E+00,
     1        .7241989908086E+00 , .2755271414712E+00,
     1        .1238035601152E+01 , .1206301926098E+00,
     1        .1838528222377E+01 , .2189228617719E-01,
     1        .2531488153459E+01 , .1236446714932E-02,
     1        .3373456432188E+01 , .1108415743492E-04/
C
      DATA (A10(I),W10(I),I=1, 5)/
     1        .1002421519705E+00 , .2484061520333E+00,
     1        .4828139660535E+00 , .3923310666530E+00,
     1        .1060949821536E+01 , .2114181930720E+00,
     1        .1779729418532E+01 , .3324666034997E-01,
     1        .2669760356100E+01 , .8248533444629E-03/
C
      DATA (A12(I),W12(I),I=1, 6)/
     1        .7860065943877E-01 , .1968496755464E+00,
     1        .3867394103619E+00 , .3491542015517E+00,
     1        .8664294718156E+00 , .2572595205372E+00,
     1        .1465698049811E+01 , .7601313755155E-01,
     1        .2172707797082E+01 , .6851918620771E-02,
     1        .3036820169458E+01 , .9847164512269E-04/
C
      DATA (A16(I),W16(I),I=1, 8)/
     1        .5297864565837E-01 , .1341091926237E+00,
     1        .2673983796319E+00 , .2683307595127E+00,
     1        .6163028979244E+00 , .2759533971747E+00,
     1        .1064246330536E+01 , .1574482776307E+00,
     1        .1588855883434E+01 , .4481410705415E-01,
     1        .2183921175421E+01 , .5367935237021E-02,
     1        .2863133906128E+01 , .2020636231371E-03,
     1        .3686007184393E+01 , .1192596732610E-05/
C
      DATA (A18(I),W18(I),I=1, 9)/
     1        .4493904292488E-01 , .1140889997309E+00,
     1        .2286053587131E+00 , .2359408305979E+00,
     1        .5321959434768E+00 , .2664254770376E+00,
     1        .9272808796501E+00 , .1832516454048E+00,
     1        .1392924011305E+01 , .7134402055945E-01,
     1        .1918843266340E+01 , .1398140961583E-01,
     1        .2506248004313E+01 , .1163851727278E-02,
     1        .3172692301930E+01 , .3056698870312E-04,
     1        .3978899031889E+01 , .1237903482316E-06/
C
      DATA (A20(I),W20(I),I=1,10)/
     1        .3873864884860E-01 , .9855240332800E-01,
     1        .1982338654186E+00 , .2086785207157E+00,
     1        .4652021935390E+00 , .2520518280321E+00,
     1        .8168633858776E+00 , .1986840280664E+00,
     1        .1234543087049E+01 , .9719804076584E-01,
     1        .1706800040403E+01 , .2702424603217E-01,
     1        .2229942013479E+01 , .3804616962480E-02,
     1        .2809105645124E+01 , .2288837837905E-03,
     1        .3463874253042E+01 , .4345289131007E-05,
     1        .4255363544009E+01 , .1247718489704E-07/
C
      DATA (A22(I),W22(I),I=1,11)/
     1        .3384061403093E-01 , .8622390849070E-01,
     1        .1739617277093E+00 , .1857725460643E+00,
     1        .4108858967918E+00 , .2358288434434E+00,
     1        .7262892878472E+00 , .2058479131548E+00,
     1        .1103884754493E+01 , .1195764488002E+00,
     1        .1532319076505E+01 , .4314132896931E-01,
     1        .2005808290066E+01 , .8866762220730E-02,
     1        .2524378003661E+01 , .9270210486154E-03,
     1        .3095377423956E+01 , .4156528700232E-04,
     1        .3739503719383E+01 , .5867468163773E-06,
     1        .4517860062176E+01 , .1226874754296E-08/
C
      DATA (A24(I),W24(I),I=1,12)/
     1        .2991099939464E-01 , .7629916824372E-01,
     1        .1543052290047E+00 , .1665373089828E+00,
     1        .3663499091549E+00 , .2194552234052E+00,
     1        .6511867037134E+00 , .2069930336427E+00,
     1        .9947501731010E+00 , .1371838706119E+00,
     1        .1386327248372E+01 , .6043861774174E-01,
     1        .1819316107995E+01 , .1652649253512E-01,
     1        .2291325454080E+01 , .2580415530996E-02,
     1        .2804582984020E+01 , .2056754688152E-03,
     1        .3367751822205E+01 , .7043545098565E-05,
     1        .4002153024932E+01 , .7562694895932E-07,
     1        .4768667906262E+01 , .1176801140394E-09/
C
      DATA (A26(I),W26(I),I=1,13)/
     1        .2687539893042E-01 , .6860743874552E-01,
     1        .1389664314025E+00 , .1510724295434E+00,
     1        .3310895348420E+00 , .2044087887178E+00,
     1        .5908275292527E+00 , .2040591385184E+00,
     1        .9059199991909E+00 , .1492987609926E+00,
     1        .1266464844884E+01 , .7656911208622E-01,
     1        .1665690249468E+01 , .2601931635476E-01,
     1        .2100047720488E+01 , .5487550755778E-02,
     1        .2569279370181E+01 , .6623725822131E-03,
     1        .3077004838886E+01 , .4092137715937E-04,
     1        .3632719447215E+01 , .1086769114328E-05,
     1        .4258100209939E+01 , .8999216253810E-08,
     1        .5013695982972E+01 , .1056181764226E-10/
C
      DATA (A28(I),W28(I),I=1,14)/
     1        .2574593750171E-01 , .6573088665062E-01,
     1        .1331473976273E+00 , .1450007948865E+00,
     1        .3172834649517E+00 , .1976903952423E+00,
     1        .5662379126244E+00 , .2012086859914E+00,
     1        .8681075880846E+00 , .1528864568113E+00,
     1        .1213086106429E+01 , .8355201801999E-01,
     1        .1594180474269E+01 , .3130804321888E-01,
     1        .2007226518418E+01 , .7620883072174E-02,
     1        .2450765117455E+01 , .1130613659204E-02,
     1        .2926155234545E+01 , .9416408715712E-04,
     1        .3438309154336E+01 , .3916031412192E-05,
     1        .3997895360339E+01 , .6744233894962E-07,
     1        .4628038787602E+01 , .3391774320172E-09,
     1        .5392407922630E+01 , .2070921821819E-12/
C
C     GAUSS HERMITE QUADRATURE
C         IS = 1 HALF RANGE (MUST AS=0)
C         IS = 2 FULL RANGE
C
      IF(NV .EQ. 4) THEN
      DO 210 I = 1,2
      I1       = 2 + I
      I2       = 2 - I + 1
      A(I1)    = AS + A4(I)
      A(I2)    = AS - A4(I)
      W(I1)    = W4(I)*EXP(A4(I)*A4(I))
      W(I2)    = W4(I)*EXP(A4(I)*A4(I))
  210 CONTINUE
C
      ELSEIF(NV .EQ. 6) THEN
      DO 212 I = 1,3
      I1       = 3 + I
      I2       = 3 - I + 1
      A(I1)    = AS + A6(I)
      A(I2)    = AS - A6(I)
      W(I1)    = W6(I)*EXP(A6(I)*A6(I))
      W(I2)    = W6(I)*EXP(A6(I)*A6(I))
  212 CONTINUE
C
      ELSEIF(NV .EQ. 8) THEN
      DO 214 I = 1,4
      I1       = 4 + I
      I2       = 4 - I + 1
      A(I1)    = AS + A8(I)
      A(I2)    = AS - A8(I)
      W(I1)    = W8(I)*EXP(A8(I)*A8(I))
      W(I2)    = W8(I)*EXP(A8(I)*A8(I))
  214 CONTINUE
C
      ELSEIF(NV .EQ.10) THEN
      DO 216 I = 1,5
      I1       = 5 + I
      I2       = 5 - I + 1
      A(I1)    = AS + A10(I)
      A(I2)    = AS - A10(I)
      W(I1)    = W10(I)*EXP(A10(I)*A10(I))
      W(I2)    = W10(I)*EXP(A10(I)*A10(I))
  216 CONTINUE
C
      ELSEIF(NV .EQ.12) THEN
      DO 218 I = 1,6
      I1       = 6 + I
      I2       = 6 - I + 1
      A(I1)    = AS + A12(I)
      A(I2)    = AS - A12(I)
      W(I1)    = W12(I)*EXP(A12(I)*A12(I))
      W(I2)    = W12(I)*EXP(A12(I)*A12(I))
  218 CONTINUE
C
      ELSEIF(NV .EQ.14) THEN
      DO 220 I = 1,7
      I1       = 7 + I
      I2       = 7 - I + 1
      A(I1)    = AS + A14(I)
      A(I2)    = AS - A14(I)
      W(I1)    = W14(I)*EXP(A14(I)*A14(I))
      W(I2)    = W14(I)*EXP(A14(I)*A14(I))
  220 CONTINUE
C
      ELSEIF(NV .EQ.16) THEN
      DO 222 I = 1,8
      I1       = 8 + I
      I2       = 8 - I + 1
      A(I1)    = AS + A16(I)
      A(I2)    = AS - A16(I)
      W(I1)    = W16(I)*EXP(A16(I)*A16(I))
      W(I2)    = W16(I)*EXP(A16(I)*A16(I))
  222 CONTINUE
C
      ELSEIF(NV .EQ.18) THEN
      DO 224 I = 1,9
      I1       = 9 + I
      I2       = 9 - I + 1
      A(I1)    = AS + A18(I)
      A(I2)    = AS - A18(I)
      W(I1)    = W18(I)*EXP(A18(I)*A18(I))
      W(I2)    = W18(I)*EXP(A18(I)*A18(I))
  224 CONTINUE
C
      ELSEIF(NV .EQ.20) THEN
      DO 226 I = 1,10
      I1       = 10 + I
      I2       = 10 - I + 1
      A(I1)    = AS + A20(I)
      A(I2)    = AS - A20(I)
      W(I1)    = W20(I)*EXP(A20(I)*A20(I))
      W(I2)    = W20(I)*EXP(A20(I)*A20(I))
  226 CONTINUE
C
      ELSEIF(NV .EQ.22) THEN
      DO 228 I = 1,11
      I1       = 11 + I
      I2       = 11 - I + 1
      A(I1)    = AS + A22(I)
      A(I2)    = AS - A22(I)
      W(I1)    = W22(I)*EXP(A22(I)*A22(I))
      W(I2)    = W22(I)*EXP(A22(I)*A22(I))
  228 CONTINUE
C
      ELSEIF(NV .EQ.24) THEN
      DO 230 I = 1,12
      I1       = 12 + I
      I2       = 12 - I + 1
      A(I1)    = AS + A24(I)
      A(I2)    = AS - A24(I)
      W(I1)    = W24(I)*EXP(A24(I)*A24(I))
      W(I2)    = W24(I)*EXP(A24(I)*A24(I))
  230 CONTINUE
C
      ELSEIF(NV .EQ.26) THEN
      DO 232 I = 1,13
      I1       = 13 + I
      I2       = 13 - I + 1
      A(I1)    = AS + A26(I)
      A(I2)    = AS - A26(I)
      W(I1)    = W26(I)*EXP(A26(I)*A26(I))
      W(I2)    = W26(I)*EXP(A26(I)*A26(I))
  232 CONTINUE
C
      ELSEIF(NV .EQ.28) THEN
      DO 234 I = 1,14
      I1       = 14 + I
      I2       = 14 - I + 1
      A(I1)    = AS + A28(I)
      A(I2)    = AS - A28(I)
      W(I1)    = W28(I)*EXP(A28(I)*A28(I))
      W(I2)    = W28(I)*EXP(A28(I)*A28(I))
  234 CONTINUE
C
      ELSE
       WRITE(*,*)'VHUANG ERROR : NV NOT MATCH   NV=',NV
       STOP
      ENDIF
C
      IF(IS.EQ.1) THEN
        IV1   = NV/2 + 1
        IV2   =  NV
        NV    = IV2 - IV1 + 1
        DO 310 I = IV1,IV2
        II    = IV1 - I + 1
        A(II) = A(I)
        W(II) = W(I)
  310   CONTINUE
      WRITE(*,*)'HALF RANGE HUANG QUADRATURE, NV=',NV
      ELSE
      WRITE(*,*)'FULL RANGE HUANG QUADRATURE, NV=',NV
      ENDIF
      RETURN

      END


















      SUBROUTINE VLAGUE(NV,W,A,AS,IS,MV)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION W(MV),A(MV)
      DIMENSION A2(2),A3(3),A4(4),A5(5),A6(6),A10(10),A15(15)
      DIMENSION W2(2),W3(3),W4(4),W5(5),W6(6),W10(10),W15(15)
C
      DATA (A2(I),W2(I),I=1,2)/
     1    0.585786437627,0.853553390593,
     1    3.414213562373,0.146446609407/
C
      DATA (A3(I),W3(I),I=1,3)/
     1    0.415774556783,0.711093009929,
     1    2.294280360279,0.278517733569,
     1    6.289945082937,0.103892565016E-1/
C
      DATA (A4(I),W4(I),I=1,4)/
     1    0.322547689619,0.603154104342,
     1    1.745761101158,0.357418692438,
     1    4.536620296921,0.388879085150E-1,
     1    9.395070912301,0.539294705561E-3/
C
      DATA (A5(I),W5(I),I=1,5)/
     1    0.263560319718,0.521755610583,
     1    1.423403059107,0.398666811083,
     1    3.596425771041,0.759424496817E-1,
     1    7.085810005859,0.361175867992E-2,
     1   12.640800844276,0.233699723858E-4/
C
      DATA (A6(I),W6(I),I=1,6)/
     1    0.222846604179,0.458964673950,
     1    1.188932101673,0.417000830772,
     1    2.992736326059,0.113373382074,
     1    5.775143569105,0.103991974531E-1,
     1    9.877467418383,0.261017202815E-3,
     1   15.982873980602,0.898547906430E-6/
C
      DATA (A10(I),W10(I),I=1,10)/
     1    0.137793470540,0.308441115765,
     2    0.729454549503,0.401119929155,
     3    1.808342901740,0.218068287612,
     4    3.401433697855,0.620874560987E-1,
     5    5.552496140064,0.950151697518E-2,
     6    8.330152746764,0.753008388588E-3,
     7   11.843785837900,0.282592334960E-4,
     8   16.279257831378,0.424931398496E-6,
     9   21.996585811981,0.183956482398E-8,
     1   29.920697012274,0.991182721961E-12/
C
      DATA (A15(I),W15(I),I=1,15)/
     1    0.093307812017,0.218234885940,
     2    0.492691740302,0.342210177923,
     3    1.215595412071,0.263027577942,
     4    2.269949526204,0.126425818106,
     5    3.667622721751,0.402068649210E-1,
     6    5.425336627414,0.856387780361E-2,
     7    7.565916226613,0.121243614721E-2,
     8   10.120228568019,0.111674392344E-3,
     9   13.130282482176,0.645992676202E-5,
     1   16.654407708330,0.222631690710E-6,
     1   20.776478899449,0.422743038498E-8,
     2   25.623894226729,0.392189726704E-10,
     3   31.407519169754,0.145651526407E-12,
     4   38.530683306486,0.148302705111E-15,
     5   48.026085572686,0.160059490621E-19/
C
C     GAUSS Laguerre QUADRATURE
C
      IF(NV .EQ. 2) THEN
      DO 15 I = 1,NV
      A(I)    = A2(I) + AS
      W(I)    = W2(I)*EXP(A2(I))
   15 CONTINUE
C
      ELSEIF(NV .EQ. 3) THEN
      DO 25 I = 1,NV
      A(I)    = A3(I) + AS
      W(I)    = W3(I)*EXP(A3(I))
   25 CONTINUE
C
      ELSEIF(NV .EQ. 4) THEN
      DO 35 I = 1,NV
      A(I)    = A4(I) + AS
      W(I)    = W4(I)*EXP(A4(I))
   35 CONTINUE
C
      ELSEIF(NV .EQ. 5) THEN
      DO 45 I = 1,NV
      A(I)    = A5(I) + AS
      W(I)    = W5(I)*EXP(A5(I))
   45 CONTINUE
C
      ELSEIF(NV .EQ. 6) THEN
      DO 55 I = 1,NV
      A(I)    = A6(I) + AS
      W(I)    = W6(I)*EXP(A6(I))
   55 CONTINUE
C
      ELSEIF(NV .EQ. 10) THEN
      DO 65 I = 1,NV
      A(I)    = A10(I) + AS
      W(I)    = W10(I)*EXP(A10(I))
   65 CONTINUE
C
      ELSEIF(NV .EQ. 15) THEN
      DO 75 I = 1,NV
      A(I)    = A15(I) + AS
      W(I)    = W15(I)*EXP(A15(I))
   75 CONTINUE
      ELSE
       WRITE(*,*)'VLAGURE ERROR : NV  NOT MATCH   NV=',NV
       STOP
      ENDIF
C
      IF(IS .EQ. 2) THEN
      NV2     = NV*2
      DO 100 I= 1,NV
      II      = NV + I
      A(II)   = -A(I)
      W(II)   = W(I)
  100 CONTINUE
      NV      = NV2
      WRITE(*,*)' FULL RANGE GAUSS LAGUERRE QUADRTURE POINTS NV=',NV
        IF(NV .GT. MV)THEN
          WRITE(*,*)' VLAGUE ERROR : NV > MV,   NV=',NV,'   MV=',MV
          STOP
        ENDIF
      ELSE
      WRITE(*,*)' HALF RANGE GAUSS LAGUERRE QUADRTURE POINTS NV=',NV
      ENDIF
      RETURN
      END
      SUBROUTINE XREAD( A , N , ITAPE )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N)
C
      READ(ITAPE,IOSTAT=IS) A
C
      RETURN
      END
      SUBROUTINE XWRITE( A , N , ITAPE )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N)
C
      WRITE(ITAPE,IOSTAT=IS) A
C
      RETURN
      END
      FUNCTION CVMGT(A,B,I)
      IMPLICIT REAL*8(A-H,O-Z)
      IF (I) 100, 200, 100
  100 CVMGT = A
      RETURN
  200 CVMGT = B
      RETURN
      END

      SUBROUTINE FIWENO3(DTLOC,CX,CY,CXI,XIX,XIY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   ETAX(-2:NI3,-2:NJ3), ETAY(-2:NI3,-2:NJ3)
      DIMENSION    QJ(-2:NI3,-2:NJ3)
      DIMENSION DTLOC(-2:NI3,-2:NJ3)
      DIMENSION   DQG(-2:NI3,-2:NJ3), DQH(-2:NI3,-2:NJ3)
      DIMENSION    RG(-2:NI3,-2:NJ3), RH(-2:NI3,-2:NJ3)
      DIMENSION   CXI(-2:NI3,-2:NJ3)
      DIMENSION    FG(-2:MXI-2),FH(-2:MXI-2)
C
      DIMENSION   FXG(-2:MAXIJP),     FXH(-2:MAXIJP)
      DIMENSION   FMG(-2:MAXIJP),     FMH(-2:MAXIJP)
      DIMENSION   FNG(-2:MAXIJP),     FNH(-2:MAXIJP)
      DIMENSION  CXG(-2:MAXIJP), ACXG(-2:MAXIJP)
      DIMENSION  CXH(-2:MAXIJP), ACXH(-2:MAXIJP)
      DIMENSION  DPG(-2:MAXIJP), DPH(-2:MAXIJP)
      DIMENSION DPFG(-2:MAXIJP), DPFH(-2:MAXIJP)
      DIMENSION SCXG(-2:MAXIJP), CCPG(-2:MAXIJP),  CCMG(-2:MAXIJP)
      DIMENSION SCXH(-2:MAXIJP), CCPH(-2:MAXIJP),  CCMH(-2:MAXIJP)
      DIMENSION   DDAG(-2:MAXIJP),    DDA2G(-2:MAXIJP)
      DIMENSION   DDAH(-2:MAXIJP),    DDA2H(-2:MAXIJP)
      DIMENSION   EAG(-2:MAXIJP),     EAH(-2:MAXIJP)
      DIMENSION    EG(-2:MAXIJP),      EH(-2:MAXIJP)
      DIMENSION  ADPG(-2:MAXIJP),    ADPH(-2:MAXIJP)
      DIMENSION   DMG(-2:MAXIJP),     DMH(-2:MAXIJP)
      DIMENSION  DMFG(-2:MAXIJP),    DMFH(-2:MAXIJP)
      DIMENSION  ADMG(-2:MAXIJP),    ADMH(-2:MAXIJP)
      DIMENSION   DAG(-2:MAXIJP),     DAH(-2:MAXIJP)
      DIMENSION   DBG(-2:MAXIJP),     DBH(-2:MAXIJP)
      DIMENSION DPDAG(-2:MAXIJP),   DPDAH(-2:MAXIJP)
      DIMENSION DPDBG(-2:MAXIJP),   DPDBH(-2:MAXIJP)
      DIMENSION DMDAG(-2:MAXIJP),   DMDAH(-2:MAXIJP)
      DIMENSION DMDBG(-2:MAXIJP),   DMDBH(-2:MAXIJP)
      DIMENSION    DG(-2:MAXIJP),      DH(-2:MAXIJP)
	DIMENSION    GG(-2:MAXIJP),      HH(-2:MAXIJP)
	DIMENSION  DFXG(-2:MAXIJP),    DFXH(-2:MAXIJP)
	DIMENSION  GG2(-2:MAXIJP,2),      HH2(-2:MAXIJP,2)
	DIMENSION GG1(-2:MAXIJP,4,2),   HH1(-2:MAXIJP,4,2)
	DIMENSION   FGP(-2:MAXIJP),     FHP(-2:MAXIJP)
	DIMENSION   FGM(-2:MAXIJP),         FHM(-2:MAXIJP)
C
      !EPS=0.01
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = 2.*EPS
	epweno = 1.e-6
	NIM3 =NI-3
	NIM2 =NI-2
	NIM1 =NI-1
C
      DO  J = -2,NJ3
C
      DO  I = -2,NI3
      FXG(I)   = CXI(I,J)*RG(I,J)
      FXH(I)   = CXI(I,J)*RH(I,J)
      enddo
	do i = -2,NI2 
		DFXG(i) = FXG(i+1) - FXG(i)
          DFXH(i) = FXH(i+1) - FXH(i)
		DG(i) = RG(I+1,J) - RG(I,J)
          DH(i) = RH(I+1,J) - RH(I,J)
	   GG2(i,1) = 0.5 * (DFXG(i) + abs(CXI(I,J)) * DG(i))
	   GG2(i,2) = 0.5 * (-1.*DFXG(i) + abs(CXI(I,J)) * DG(i))
         HH2(i,1) = 0.5 * (DFXH(i) + abs(CXI(I,J)) * DH(i))
	   HH2(i,2) = 0.5 * (-1.*DFXH(i) + abs(CXI(I,J)) * DH(i))	   
	end do
	   DFXG(NI3)=DFXG(NI2)
         DFXH(NI3)=DFXH(NI2)
	   DG(NI3)=DG(NI2)
         DH(NI3)=DH(NI2)
	   GG2(NI3,1)=GG2(NI2,1)
         HH2(NI3,1)=HH2(NI2,1)
	   GG2(NI3,2)=GG2(NI2,2)
         HH2(NI3,2)=HH2(NI2,2)

	do m1 = 1, 4
		k0 = m1 -  3		! (-2, -1,  0,  1)
		k1 =  3 - m1		! ( 2,  1,  0, -1)
			
	  do  i = 0, NI1
		   GG1(i,m1,1) = GG2(i+k0,1)   
		   GG1(i,m1,2) = GG2(i+k1,2)
             HH1(i,m1,1) = HH2(i+k0,1)   
		   HH1(i,m1,2) = HH2(i+k1,2)
	  end do
	end do
	       GG1(NI2,m1,1) = GG1(NI1,M1,1)   
		   GG1(NI2,m1,2) = GG1(NI1,M1,2)
             HH1(NI2,m1,1) = HH1(NI1,M1,1)   
		   HH1(NI2,m1,2) = HH1(NI1,M1,2)
	       GG1(NI3,m1,1) = GG1(NI,M1,1)   
		   GG1(NI3,m1,2) = GG1(NI,M1,2)
             HH1(NI3,m1,1) = HH1(NI,M1,1)   
		   HH1(NI3,m1,2) = HH1(NI,M1,2)
	       GG1(-1,m1,1) = GG1(0,M1,1)   
		   GG1(-1,m1,2) = GG1(0,M1,2)
             HH1(-1,m1,1) = HH1(0,M1,1)   
		   HH1(-1,m1,2) = HH1(0,M1,2)
	       GG1(-2,m1,1) = GG1(1,M1,1)   
		   GG1(-2,m1,2) = GG1(1,M1,2)
             HH1(-2,m1,1) = HH1(1,M1,1)   
		   HH1(-2,m1,2) = HH1(1,M1,2)
C		
	do i = -2, NI3
		FNG(I) = 0.
	    FNH(I) = 0.
	end do
		
	do m1 = 1, 2
	  do i = 0, NI
		TG1 = GG1(i,1,m1) - GG1(i,2,m1)
		TG2 = GG1(i,2,m1) - GG1(i,3,m1)
		TG3 = GG1(i,3,m1) - GG1(i,4,m1)
		TH1 = HH1(i,1,m1) - HH1(i,2,m1)
		TH2 = HH1(i,2,m1) - HH1(i,3,m1)
		TH3 = HH1(i,3,m1) - HH1(i,4,m1)
				
	TTG1 = 13./12. * TG1*TG1 + 1./4. 
	1      * (    GG1(i,1,m1) - 3.*GG1(i,2,m1) )**2
	TTG2 = 13./12. * TG2*TG2 + 1./4.
	1      * (    GG1(i,2,m1) +    GG1(i,3,m1) )**2
	TTG3 = 13./12. * TG3*TG3 + 1./4. 
	1      * ( 3.*GG1(i,3,m1) -    GG1(i,4,m1) )**2
	TTH1 = 13./12. * TH1*TH1 + 1./4. 
	1      * (    HH1(i,1,m1) - 3.*HH1(i,2,m1) )**2
	TTH2 = 13./12. * TH2*TH2 + 1./4.
	1      * (    HH1(i,2,m1) +    HH1(i,3,m1) )**2
	TTH3 = 13./12. * TH3*TH3 + 1./4. 
	1      * ( 3.*HH1(i,3,m1) -    HH1(i,4,m1) )**2
	!tt1=IS0
	!tt2=IS1
	!tt3=IS2
	if(m1 .eq. 1)then		
	TTTG1 =  1./(10.*( epweno + TTG1 )**2)
	TTTG2 =  6./(10.*( epweno + TTG2 )**2)
	TTTG3 =  3./(10.*( epweno + TTG3 )**2)
	TTTG=    TTTG1+TTTG2+TTTG3
	SG1=     TTTG1/TTTG
	SG2=     TTTG2/TTTG
	SG3=     TTTG3/TTTG

	TTTH1 =  1./(10.*( epweno + TTH1 )**2)
	TTTH2 =  6./(10.*( epweno + TTH2 )**2)
	TTTH3 =  3./(10.*( epweno + TTH3 )**2)
	TTTH=    TTTH1+TTTH2+TTTH3
	SH1=     TTTH1/TTTH
	SH2=     TTTH2/TTTH
	SH3=     TTTH3/TTTH
      
	
	fGpm2=0.5*(FXG(i-2)+abs(CXI(I,J))*RG(i-2,J))
      fGpm1=0.5*(FXG(i-1)+abs(CXI(I,J))*RG(i-1,J))
	fGpm0=0.5*(FXG(i)+abs(CXI(I,J))*RG(i,J))
	fGpp1=0.5*(FXG(i+1)+abs(CXI(I,J))*RG(i+1,J))
	fGpp2=0.5*(FXG(i+2)+abs(CXI(I,J))*RG(i+2,J))
	eGp1=1./6.*(2.*fGpm2-7.*fGpm1+11.*fGpm0)
	eGp2=1./6.*(-1.*fGpm1+5.*fGpm0+2.*fGpp1)
	eGp3=1./6.*(2.*fGpm0+5.*fGpp1-1.*fGpp2)

	FGP(i)=SG1*eGp1+SG2*eGp2+SG3*eGp3
C
	fHpm2=0.5*(FXH(i-2)+abs(CXI(I,J))*RH(i-2,J))
      fHpm1=0.5*(FXH(i-1)+abs(CXI(I,J))*RH(i-1,J))
	fHpm0=0.5*(FXH(i)+abs(CXI(I,J))*RH(i,J))
	fHpp1=0.5*(FXH(i+1)+abs(CXI(I,J))*RH(i+1,J))
	fHpp2=0.5*(FXH(i+2)+abs(CXI(I,J))*RH(i+2,J))
	eHp1=1./6.*(2.*fHpm2-7.*fHpm1+11.*fHpm0)
	eHp2=1./6.*(-1.*fHpm1+5.*fHpm0+2.*fHpp1)
	eHp3=1./6.*(2.*fHpm0+5.*fHpp1-1.*fHpp2)

	FHP(i)=SH1*eHp1+SH2*eHp2+SH3*eHp3
	else
	TTTG1 =  3./(10.*( epweno + TTG3 )**2)
	TTTG2 =  6./(10.*( epweno + TTG2 )**2)
	TTTG3 =  1./(10.*( epweno + TTG1 )**2)
	TTTG=    TTTG1+TTTG2+TTTG3
	SG1=     TTTG1/TTTG
	SG2=     TTTG2/TTTG
	SG3=     TTTG3/TTTG

	TTTH1 =  3./(10.*( epweno + TTH3 )**2)
	TTTH2 =  6./(10.*( epweno + TTH2 )**2)
	TTTH3 =  1./(10.*( epweno + TTH1 )**2)
	TTTH=    TTTH1+TTTH2+TTTH3
	SH1=     TTTH1/TTTH
	SH2=     TTTH2/TTTH
	SH3=     TTTH3/TTTH

	fGmm1=0.5*(FXG(i-1)-abs(CXI(I,J))*RG(i-1,J))
	fGmm0=0.5*(FXG(i)-abs(CXI(I,J))*RG(i,J))
	fGmp1=0.5*(FXG(i+1)-abs(CXI(I,J))*RG(i+1,J))
	fGmp2=0.5*(FXG(i+2)-abs(CXI(I,J))*RG(i+2,J))
	fGmp3=0.5*(FXG(i+3)-abs(CXI(I,J))*RG(i+3,J))
	eGm1=1./6.*(-1.*fGmm1+5.*fGmm0+2.*fGmp1)
	eGm2=1./6.*(2.*fGmm0+5.*fGmp1-1.*fGmp2)
	eGm3=1./6.*(11.*fGmp1-7.*fGmp2+2.*fGmp3)

	FGM(i)=SG1*eGm1+SG2*eGm2+SG3*eGm3
C
	fHmm1=0.5*(FXH(i-1)-abs(CXI(I,J))*RH(i-1,J))
	fHmm0=0.5*(FXH(i)-abs(CXI(I,J))*RH(i,J))
	fHmp1=0.5*(FXH(i+1)-abs(CXI(I,J))*RH(i+1,J))
	fHmp2=0.5*(FXH(i+2)-abs(CXI(I,J))*RH(i+2,J))
	fHmp3=0.5*(FXH(i+3)-abs(CXI(I,J))*RH(i+3,J))
	eHm1=1./6.*(-1.*fHmm1+5.*fHmm0+2.*fHmp1)
	eHm2=1./6.*(2.*fHmm0+5.*fHmp1-1.*fHmp2)
	eHm3=1./6.*(11.*fHmp1-7.*fHmp2+2.*fHmp3)

	FHM(i)=SH1*eHm1+SH2*eHm2+SH3*eHm3
	endif
C
	  end do
	end do
 
      do I=0,NI
	FNG(I)=FGP(I)+FGM(I)
	FNH(I)=FHP(I)+FHM(I)
	end do
C
      FNG(-1)=FNG(0)
	FNH(-1)=FNH(0)
	FNG(-2)=FNG(1)
	FNH(-2)=FNH(1)
      FNG(NI1)=FNG(NI)
	FNH(NI1)=FNH(NI)
	FNG(NI2)=FNG(NIM1)
      FNH(NI2)=FNH(NIM1)
	FNG(NI3)=FNG(NIM2)
	FNH(NI3)=FNH(NIM2)
C
      DO  I = -1,NI3
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FNG(I) - FNG(I-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FNH(I) - FNH(I-1))
      enddo
	DQG(-2,J)=DQG(-1,J)
	DQH(-2,J)=DQH(-1,J)
      enddo
C
      
!      I = 4
!      IF(CX.LT.0.)THEN
!         DGG      = u(I+1,K) - u(I,K)
!         rhs(i)    = rhs(i) - DTDX*CX*DGG
!      ENDIF
c
!      I = NIM3
!      IF(CX.GT.0.)THEN
!         DGG        = u(I,K) - u(I-1,K)
!         rhs(i)      = rhs(i)  - DTDX*CX*DGG
!      ENDIF
!      rhs(1)=rhs(7)
!	rhs(2)=rhs(6)
!     rhs(3)=rhs(5)
!	rhs(MI)=rhs(NIM6)
!	rhs(NIM1)=rhs(NIM5)
!	rhs(NIM2)=rhs(NIM4)
      RETURN
      END
C
C
C
      SUBROUTINE FJWENO3(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   ETAX(-2:NI3,-2:NJ3), ETAY(-2:NI3,-2:NJ3)
      DIMENSION    QJ(-2:NI3,-2:NJ3)
      DIMENSION DTLOC(-2:NI3,-2:NJ3)
      DIMENSION   DQG(-2:NI3,-2:NJ3), DQH(-2:NI3,-2:NJ3)
      DIMENSION    RG(-2:NI3,-2:NJ3), RH(-2:NI3,-2:NJ3)
      DIMENSION   CET(-2:NI3,-2:NJ3)
      DIMENSION    FG(-2:MXI-2),FH(-2:MXI-2)
C
      DIMENSION   FXG(-2:MAXIJP),     FXH(-2:MAXIJP)
      DIMENSION   FMG(-2:MAXIJP),     FMH(-2:MAXIJP)
      DIMENSION   FNG(-2:MAXIJP),     FNH(-2:MAXIJP)
      DIMENSION  CXG(-2:MAXIJP), ACXG(-2:MAXIJP)
      DIMENSION  CXH(-2:MAXIJP), ACXH(-2:MAXIJP)
      DIMENSION  DPG(-2:MAXIJP), DPH(-2:MAXIJP)
      DIMENSION DPFG(-2:MAXIJP), DPFH(-2:MAXIJP)
      DIMENSION SCXG(-2:MAXIJP), CCPG(-2:MAXIJP),  CCMG(-2:MAXIJP)
      DIMENSION SCXH(-2:MAXIJP), CCPH(-2:MAXIJP),  CCMH(-2:MAXIJP)
      DIMENSION   DDAG(-2:MAXIJP),    DDA2G(-2:MAXIJP)
      DIMENSION   DDAH(-2:MAXIJP),    DDA2H(-2:MAXIJP)
      DIMENSION   EAG(-2:MAXIJP),     EAH(-2:MAXIJP)
      DIMENSION    EG(-2:MAXIJP),      EH(-2:MAXIJP)
      DIMENSION  ADPG(-2:MAXIJP),    ADPH(-2:MAXIJP)
      DIMENSION   DMG(-2:MAXIJP),     DMH(-2:MAXIJP)
      DIMENSION  DMFG(-2:MAXIJP),    DMFH(-2:MAXIJP)
      DIMENSION  ADMG(-2:MAXIJP),    ADMH(-2:MAXIJP)
      DIMENSION   DAG(-2:MAXIJP),     DAH(-2:MAXIJP)
      DIMENSION   DBG(-2:MAXIJP),     DBH(-2:MAXIJP)
      DIMENSION DPDAG(-2:MAXIJP),   DPDAH(-2:MAXIJP)
      DIMENSION DPDBG(-2:MAXIJP),   DPDBH(-2:MAXIJP)
      DIMENSION DMDAG(-2:MAXIJP),   DMDAH(-2:MAXIJP)
      DIMENSION DMDBG(-2:MAXIJP),   DMDBH(-2:MAXIJP)
      DIMENSION    DG(-2:MAXIJP),      DH(-2:MAXIJP)
	DIMENSION    GG(-2:MAXIJP),      HH(-2:MAXIJP)
	DIMENSION  DFXG(-2:MAXIJP),    DFXH(-2:MAXIJP)
	DIMENSION  GG2(-2:MAXIJP,2),      HH2(-2:MAXIJP,2)
	DIMENSION GG1(-2:MAXIJP,4,2),   HH1(-2:MAXIJP,4,2)
	DIMENSION   FGP(-2:MAXIJP),     FHP(-2:MAXIJP)
	DIMENSION   FGM(-2:MAXIJP),         FHM(-2:MAXIJP)
C
      !EPS=0.01
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = 2.*EPS
	epweno = 1.e-6
      NJM3=NJ-3
	NJM2=NJ-2
	NJM1=NJ-1
C
      DO  I = -2,NI3
C
      DO  J = -2,NJ3
      FXG(J)   = CET(I,J)*RG(I,J)
      FXH(J)   = CET(I,J)*RH(I,J)
      enddo
	do J = -2,NJ2 
		DFXG(J) = FXG(J+1) - FXG(J)
          DFXH(J) = FXH(J+1) - FXH(J)
		DG(J) = RG(I+1,J) - RG(I,J)
          DH(J) = RH(I+1,J) - RH(I,J)
	   GG2(J,1) = 0.5 * (DFXG(J) + abs(CET(I,J)) * DG(J))
	   GG2(J,2) = 0.5 * (-1.*DFXG(J) + abs(CET(I,J)) * DG(J))
         HH2(J,1) = 0.5 * (DFXH(J) + abs(CET(I,J)) * DH(J))
	   HH2(J,2) = 0.5 * (-1.*DFXH(J) + abs(CET(I,J)) * DH(J))	   
	end do
	   DFXG(NJ3)=DFXG(NJ2)
         DFXH(NJ3)=DFXH(NJ2)
	   DG(NJ3)=DG(NJ2)
         DH(NJ3)=DH(NJ2)
	   GG2(NJ3,1)=GG2(NJ2,1)
         HH2(NJ3,1)=HH2(NJ2,1)
	   GG2(NJ3,2)=GG2(NJ2,2)
         HH2(NJ3,2)=HH2(NJ2,2)

	do m1 = 1, 4
		k0 = m1 -  3		! (-2, -1,  0,  1)
		k1 =  3 - m1		! ( 2,  1,  0, -1)
			
	  do  J = 0, NJ1
		   GG1(J,m1,1) = GG2(J+k0,1)   
		   GG1(J,m1,2) = GG2(J+k1,2)
             HH1(J,m1,1) = HH2(J+k0,1)   
		   HH1(J,m1,2) = HH2(J+k1,2)
	  end do
	end do
	       GG1(NJ2,m1,1) = GG1(NJ1,M1,1)   
		   GG1(NJ2,m1,2) = GG1(NJ1,M1,2)
             HH1(NJ2,m1,1) = HH1(NJ1,M1,1)   
		   HH1(NJ2,m1,2) = HH1(NJ1,M1,2)
	       GG1(NJ3,m1,1) = GG1(NJ,M1,1)   
		   GG1(NJ3,m1,2) = GG1(NJ,M1,2)
             HH1(NJ3,m1,1) = HH1(NJ,M1,1)   
		   HH1(NJ3,m1,2) = HH1(NJ,M1,2)
	       GG1(-1,m1,1) = GG1(0,M1,1)   
		   GG1(-1,m1,2) = GG1(0,M1,2)
             HH1(-1,m1,1) = HH1(0,M1,1)   
		   HH1(-1,m1,2) = HH1(0,M1,2)
	       GG1(-2,m1,1) = GG1(1,M1,1)   
		   GG1(-2,m1,2) = GG1(1,M1,2)
             HH1(-2,m1,1) = HH1(1,M1,1)   
		   HH1(-2,m1,2) = HH1(1,M1,2)
C
	do J = -2, NJ3
		FNG(J) = 0.
	    FNH(J) = 0.
	end do
		
	do m1 = 1, 2
	  do J = 0, NJ
		TG1 = GG1(J,1,m1) - GG1(J,2,m1)
		TG2 = GG1(J,2,m1) - GG1(J,3,m1)
		TG3 = GG1(J,3,m1) - GG1(J,4,m1)
		TH1 = HH1(J,1,m1) - HH1(J,2,m1)
		TH2 = HH1(J,2,m1) - HH1(J,3,m1)
		TH3 = HH1(J,3,m1) - HH1(J,4,m1)
				
	TTG1 = 13./12. * TG1*TG1 + 1./4. 
	1      * (    GG1(J,1,m1) - 3.*GG1(J,2,m1) )**2
	TTG2 = 13./12. * TG2*TG2 + 1./4.
	1      * (    GG1(J,2,m1) +    GG1(J,3,m1) )**2
	TTG3 = 13./12. * TG3*TG3 + 1./4. 
	1      * ( 3.*GG1(J,3,m1) -    GG1(J,4,m1) )**2
	TTH1 = 13./12. * TH1*TH1 + 1./4. 
	1      * (    HH1(J,1,m1) - 3.*HH1(J,2,m1) )**2
	TTH2 = 13./12. * TH2*TH2 + 1./4.
	1      * (    HH1(J,2,m1) +    HH1(J,3,m1) )**2
	TTH3 = 13./12. * TH3*TH3 + 1./4. 
	1      * ( 3.*HH1(J,3,m1) -    HH1(J,4,m1) )**2
	!tt1=IS0
	!tt2=IS1
	!tt3=IS2
	if(m1 .eq. 1)then		
	TTTG1 =  1./(10.*( epweno + TTG1 )**2)
	TTTG2 =  6./(10.*( epweno + TTG2 )**2)
	TTTG3 =  3./(10.*( epweno + TTG3 )**2)
	TTTG=    TTTG1+TTTG2+TTTG3
	SG1=     TTTG1/TTTG
	SG2=     TTTG2/TTTG
	SG3=     TTTG3/TTTG

	TTTH1 =  1./(10.*( epweno + TTH1 )**2)
	TTTH2 =  6./(10.*( epweno + TTH2 )**2)
	TTTH3 =  3./(10.*( epweno + TTH3 )**2)
	TTTH=    TTTH1+TTTH2+TTTH3
	SH1=     TTTH1/TTTH
	SH2=     TTTH2/TTTH
	SH3=     TTTH3/TTTH
      
	
	fGpm2=0.5*(FXG(J-2)+abs(CET(I,J))*RG(i,J-2))
      fGpm1=0.5*(FXG(J-1)+abs(CET(I,J))*RG(i,J-1))
	fGpm0=0.5*(FXG(J)+abs(CET(I,J))*RG(i,J))
	fGpp1=0.5*(FXG(J+1)+abs(CET(I,J))*RG(i,J+1))
	fGpp2=0.5*(FXG(J+2)+abs(CET(I,J))*RG(i,J+2))
	eGp1=1./6.*(2.*fGpm2-7.*fGpm1+11.*fGpm0)
	eGp2=1./6.*(-1.*fGpm1+5.*fGpm0+2.*fGpp1)
	eGp3=1./6.*(2.*fGpm0+5.*fGpp1-1.*fGpp2)

	FGP(J)=SG1*eGp1+SG2*eGp2+SG3*eGp3
C
	fHpm2=0.5*(FXH(J-2)+abs(CET(I,J))*RH(i,J-2))
      fHpm1=0.5*(FXH(J-1)+abs(CET(I,J))*RH(i,J-1))
	fHpm0=0.5*(FXH(J)+abs(CET(I,J))*RH(i,J))
	fHpp1=0.5*(FXH(J+1)+abs(CET(I,J))*RH(i,J+1))
	fHpp2=0.5*(FXH(J+2)+abs(CET(I,J))*RH(i,J+2))
	eHp1=1./6.*(2.*fHpm2-7.*fHpm1+11.*fHpm0)
	eHp2=1./6.*(-1.*fHpm1+5.*fHpm0+2.*fHpp1)
	eHp3=1./6.*(2.*fHpm0+5.*fHpp1-1.*fHpp2)

	FHP(J)=SH1*eHp1+SH2*eHp2+SH3*eHp3
	else
	TTTG1 =  3./(10.*( epweno + TTG3 )**2)
	TTTG2 =  6./(10.*( epweno + TTG2 )**2)
	TTTG3 =  1./(10.*( epweno + TTG1 )**2)
	TTTG=    TTTG1+TTTG2+TTTG3
	SG1=     TTTG1/TTTG
	SG2=     TTTG2/TTTG
	SG3=     TTTG3/TTTG

	TTTH1 =  3./(10.*( epweno + TTH3 )**2)
	TTTH2 =  6./(10.*( epweno + TTH2 )**2)
	TTTH3 =  1./(10.*( epweno + TTH1 )**2)
	TTTH=    TTTH1+TTTH2+TTTH3
	SH1=     TTTH1/TTTH
	SH2=     TTTH2/TTTH
	SH3=     TTTH3/TTTH

	fGmm1=0.5*(FXG(J-1)-abs(CET(I,J))*RG(i,J-1))
	fGmm0=0.5*(FXG(J)-abs(CET(I,J))*RG(i,J))
	fGmp1=0.5*(FXG(J+1)-abs(CET(I,J))*RG(i,J+1))
	fGmp2=0.5*(FXG(J+2)-abs(CET(I,J))*RG(i,J+2))
	fGmp3=0.5*(FXG(J+3)-abs(CET(I,J))*RG(i,J+3))
	eGm1=1./6.*(-1.*fGmm1+5.*fGmm0+2.*fGmp1)
	eGm2=1./6.*(2.*fGmm0+5.*fGmp1-1.*fGmp2)
	eGm3=1./6.*(11.*fGmp1-7.*fGmp2+2.*fGmp3)

	FGM(J)=SG1*eGm1+SG2*eGm2+SG3*eGm3
C
	fHmm1=0.5*(FXH(J-1)-abs(CET(I,J))*RH(i,J-1))
	fHmm0=0.5*(FXH(J)-abs(CET(I,J))*RH(i,J))
	fHmp1=0.5*(FXH(J+1)-abs(CET(I,J))*RH(i,J+1))
	fHmp2=0.5*(FXH(J+2)-abs(CET(I,J))*RH(i,J+2))
	fHmp3=0.5*(FXH(J+3)-abs(CET(I,J))*RH(i,J+3))
	eHm1=1./6.*(-1.*fHmm1+5.*fHmm0+2.*fHmp1)
	eHm2=1./6.*(2.*fHmm0+5.*fHmp1-1.*fHmp2)
	eHm3=1./6.*(11.*fHmp1-7.*fHmp2+2.*fHmp3)

	FHM(J)=SH1*eHm1+SH2*eHm2+SH3*eHm3
	endif
C
	  end do
	end do
 
      do J=0,NJ
	FNG(J)=FGP(J)+FGM(J)
	FNH(J)=FHP(J)+FHM(J)
	end do

      FNG(-1)=FNG(0)
	FNH(-1)=FNH(0)
	FNG(-2)=FNG(1)
	FNH(-2)=FNH(1)
      FNG(NJ1)=FNG(NJ)
	FNH(NJ1)=FNH(NJ)
	FNG(NJ2)=FNG(NJM1)
      FNH(NJ2)=FNH(NJM1)
	FNG(NJ3)=FNG(NJM2)
	FNH(NJ3)=FNH(NJM2)

      DO  J = -1,NJ3
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FNG(J) - FNG(J-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FNH(J) - FNH(J-1))
      enddo
	      DQG(I,-2)=DQG(I,-1)
	      DQH(I,-2)=DQH(I,-1)
      enddo
C

!      I = 4
!      IF(CX.LT.0.)THEN
!         DGG      = u(I+1,K) - u(I,K)
!         rhs(i)    = rhs(i) - DTDX*CX*DGG
!      ENDIF
c
!      I = NIM3
!      IF(CX.GT.0.)THEN
!         DGG        = u(I,K) - u(I-1,K)
!         rhs(i)      = rhs(i)  - DTDX*CX*DGG
!      ENDIF
!      rhs(1)=rhs(7)
!	rhs(2)=rhs(6)
!     rhs(3)=rhs(5)
!	rhs(MI)=rhs(NIM6)
!	rhs(NIM1)=rhs(NIM5)
!	rhs(NIM2)=rhs(NIM4)
      RETURN
      END
C
C
C
      SUBROUTINE FIWENO2(DTLOC,CX,CY,CXI,XIX,XIY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   ETAX(-1:NI2,-1:NJ2), ETAY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CXI(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION   FXG(-1:MAXIJP),     FXH(-1:MAXIJP)
      DIMENSION   FMG(-1:MAXIJP),     FMH(-1:MAXIJP)
      DIMENSION   FNG(-1:MAXIJP),     FNH(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  CXH(-1:MAXIJP), ACXH(-1:MAXIJP)
      DIMENSION  DPG(-1:MAXIJP), DPH(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP), DPFH(-1:MAXIJP)
      DIMENSION SCXG(-1:MAXIJP), CCPG(-1:MAXIJP),  CCMG(-1:MAXIJP)
      DIMENSION SCXH(-1:MAXIJP), CCPH(-1:MAXIJP),  CCMH(-1:MAXIJP)
      DIMENSION   DDAG(-1:MAXIJP),    DDA2G(-1:MAXIJP)
      DIMENSION   DDAH(-1:MAXIJP),    DDA2H(-1:MAXIJP)
      DIMENSION   EAG(-1:MAXIJP),     EAH(-1:MAXIJP)
      DIMENSION    EG(-1:MAXIJP),      EH(-1:MAXIJP)
      DIMENSION  ADPG(-1:MAXIJP),    ADPH(-1:MAXIJP)
      DIMENSION   DMG(-1:MAXIJP),     DMH(-1:MAXIJP)
      DIMENSION  DMFG(-1:MAXIJP),    DMFH(-1:MAXIJP)
      DIMENSION  ADMG(-1:MAXIJP),    ADMH(-1:MAXIJP)
      DIMENSION   DAG(-1:MAXIJP),     DAH(-1:MAXIJP)
      DIMENSION   DBG(-1:MAXIJP),     DBH(-1:MAXIJP)
      DIMENSION DPDAG(-1:MAXIJP),   DPDAH(-1:MAXIJP)
      DIMENSION DPDBG(-1:MAXIJP),   DPDBH(-1:MAXIJP)
      DIMENSION DMDAG(-1:MAXIJP),   DMDAH(-1:MAXIJP)
      DIMENSION DMDBG(-1:MAXIJP),   DMDBH(-1:MAXIJP)
      DIMENSION    DG(-1:MAXIJP),      DH(-1:MAXIJP)
	DIMENSION    GG(-1:MAXIJP),      HH(-1:MAXIJP)
	DIMENSION  DFXG(-1:MAXIJP),    DFXH(-1:MAXIJP)
	DIMENSION  GG2(MAXIJP,2),      HH2(MAXIJP,2)
	DIMENSION GG1(MAXIJP,3,2),   HH1(MAXIJP,3,2)
	DIMENSION   FGP(-1:MAXIJP),     FHP(-1:MAXIJP)
	DIMENSION   FGM(-1:MAXIJP),    FHM(-1:MAXIJP)
C
      !EPS=0.01
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = 2.*EPS
	epweno = 1.e-6
C
      DO  J = -1,NJ2
C
      DO  I = -1,NI2
      FXG(I)   = CXI(I,J)*RG(I,J)
      FXH(I)   = CXI(I,J)*RH(I,J)
      enddo
	do i = -1,NI1 
		DFXG(i) = FXG(i+1) - FXG(i)
          DFXH(i) = FXH(i+1) - FXH(i)
		DG(i) = RG(I+1,J) - RG(I,J)
          DH(i) = RH(I+1,J) - RH(I,J)
	   GG2(i,1) = 0.5 * (DFXG(i) + abs(CXI(I,J)) * DG(i))
	   GG2(i,2) = 0.5 * (-1.*DFXG(i) + abs(CXI(I,J)) * DG(i))
         HH2(i,1) = 0.5 * (DFXH(i) + abs(CXI(I,J)) * DH(i))
	   HH2(i,2) = 0.5 * (-1.*DFXH(i) + abs(CXI(I,J)) * DH(i))	   
	end do
	   DFXG(NI2)=DFXG(NI1)
         DFXH(NI2)=DFXH(NI1)
	   DG(NI2)=DG(NI1)
         DH(NI2)=DH(NI1)
	   GG2(NI2,1)=GG2(NI1,1)
         HH2(NI2,1)=HH2(NI1,1)
	   GG2(NI2,2)=GG2(NI1,2)
         HH2(NI2,2)=HH2(NI1,2)

	do m1 = 1, 3
		k0 = m1 -  2		! ( -1,  0,  1)
		k1 =  2 - m1		! (  1,  0, -1)
			
	  do  i = 0, NI1
		   GG1(i,m1,1) = GG2(i+k0,1)   
		   GG1(i,m1,2) = GG2(i+k1,2)
             HH1(i,m1,1) = HH2(i+k0,1)   
		   HH1(i,m1,2) = HH2(i+k1,2)
	  end do
	       GG1(NI2,M1,1)=GG1(NI1,M1,1)
	       GG1(NI2,M1,2)=GG1(NI1,M1,2)
             HH1(NI2,M1,1)=HH1(NI1,M1,1)
	       HH1(NI2,M1,2)=HH1(NI1,M1,2)
	       GG1(-1,M1,1)=GG1(0,M1,1)
	       GG1(-1,M1,2)=GG1(0,M1,2)
             HH1(-1,M1,1)=HH1(0,M1,1)
	       HH1(-1,M1,2)=HH1(0,M1,2)
	end do
		
	do i = -1, NI2
		FNG(I) = 0.
	    FNH(I) = 0.
	end do
		
	do m1 = 1, 2
	  do i = 0, NI			
	if(m1 .eq. 1)then		
	TTTG1 =  1./(3.*( epweno + (GG1(I,1,1))**2 )**2)
	TTTG2 =  2./(3.*( epweno + (GG1(I,2,1))**2 )**2)
	TTTG=    TTTG1+TTTG2
	SG1=     TTTG1/TTTG
	SG2=     TTTG2/TTTG
	TTTH1 =  1./(3.*( epweno + (HH1(I,1,1))**2 )**2)
	TTTH2 =  2./(3.*( epweno + (HH1(I,2,1))**2 )**2)
	TTTH=    TTTH1+TTTH2
	SH1=     TTTH1/TTTH
	SH2=     TTTH2/TTTH	      

      fGpm1=0.5*(FXG(i-1)+abs(CXI(I,J))*RG(i-1,J))
	fGpm0=0.5*(FXG(i)+abs(CXI(I,J))*RG(i,J))
	fGpp1=0.5*(FXG(i+1)+abs(CXI(I,J))*RG(i+1,J))
	eGp1=1./2.*(3.*fGpm0-1.*fGpm1)
	eGp2=1./2.*(fGpm0+fGpp1)

	FGP(i)=SG1*eGp1+SG2*eGp2
!	WRITE(*,*) FGP(3),FGP(20)
!	PAUSE
C
      fHpm1=0.5*(FXH(i-1)+abs(CXI(I,J))*RH(i-1,J))
	fHpm0=0.5*(FXH(i)+abs(CXI(I,J))*RH(i,J))
	fHpp1=0.5*(FXH(i+1)+abs(CXI(I,J))*RH(i+1,J))
	eHp1=1./2.*(3.*fHpm0-1.*fHpm1)
	eHp2=1./2.*(fHpm0+fHpp1)

	FHP(i)=SH1*eHp1+SH2*eHp2
	else
	TTTG1 =  2./(3.*( epweno + (GG1(I,2,2))**2 )**2)
	TTTG2 =  1./(3.*( epweno + (GG1(I,1,2))**2 )**2)
	TTTG=    TTTG1+TTTG2
	SG1=     TTTG1/TTTG
	SG2=     TTTG2/TTTG
	TTTH1 =  2./(3.*( epweno + (HH1(I,2,2))**2 )**2)
	TTTH2 =  1./(3.*( epweno + (HH1(I,1,2))**2 )**2)
	TTTH=    TTTH1+TTTH2
	SH1=     TTTH1/TTTH
	SH2=     TTTH2/TTTH
C
	fGmm0=0.5*(FXG(i)-abs(CXI(I,J))*RG(i,J))
	fGmp1=0.5*(FXG(i+1)-abs(CXI(I,J))*RG(i+1,J))
	fGmp2=0.5*(FXG(i+2)-abs(CXI(I,J))*RG(i+2,J))
	eGm1=1./2.*(fGmm0+1.*fGmp1)
	eGm2=1./2.*(3.*fGmp1-1.*fGmp2)

	FGM(i)=SG1*eGm1+SG2*eGm2
C
	fHmm0=0.5*(FXH(i)-abs(CXI(I,J))*RH(i,J))
	fHmp1=0.5*(FXH(i+1)-abs(CXI(I,J))*RH(i+1,J))
	fHmp2=0.5*(FXH(i+2)-abs(CXI(I,J))*RH(i+2,J))
	eHm1=1./2.*(fHmm0+1.*fHmp1)
	eHm2=1./2.*(3.*fHmp1-1.*fHmp2)

	FHM(i)=SH1*eHm1+SH2*eHm2
	endif

	  end do
	end do
 
      do I=0,NI
	FNG(I)=FGP(I)+FGM(I)
	FNH(I)=FHP(I)+FHM(I)
	end do

      FNG(-1)=FNG(0)
	FNH(-1)=FNH(0)
      FNG(NI1)=FNG(NI)
	FNH(NI1)=FNH(NI)
	FNG(NI2) = FNG(NI1)
      FNH(NI2) = FNH(NI1)
      DO  I = 0,NI2
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FNG(I) - FNG(I-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FNH(I) - FNH(I-1))
      enddo
      enddo
C
      RETURN
      END
C
C
C
      SUBROUTINE FJWENO2(DTLOC,CX,CY,CET,ETAX,ETAY,QJ,
     1                          FG,FH,RG,RH,DQG,DQH)
      PARAMETER (MAXIJP=305)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON /MAX1 / MXW ,MXI ,MXB ,MXG , MXL,MXV ,MXBP ,MXLI
      COMMON /DIM  / NI,NJ,NI1,NJ1,NI2,NJ2,NI3,NJ3,NDIM(2),
     &               NVX,NVY,IS2D,KS2D,ISTA,IEND,JSTA,JEND,MI,MJ
      COMMON /TIME0I/ ITER,NUP,ILOCDT
      COMMON /TIME0R/ TIME,DT,DTI,DTJ,DTCFL,CFL,CFL1,CFL2,DTFIX
      COMMON /UNST1I/ LSTD,MTIM
      COMMON /UNST1R/ XMS,XST,TPP(70)
      COMMON /ENTR / ENTEPS
C
      DIMENSION   ETAX(-1:NI2,-1:NJ2), ETAY(-1:NI2,-1:NJ2)
      DIMENSION    QJ(-1:NI2,-1:NJ2)
      DIMENSION DTLOC(-1:NI2,-1:NJ2)
      DIMENSION   DQG(-1:NI2,-1:NJ2), DQH(-1:NI2,-1:NJ2)
      DIMENSION    RG(-1:NI2,-1:NJ2), RH(-1:NI2,-1:NJ2)
      DIMENSION   CET(-1:NI2,-1:NJ2)
      DIMENSION    FG(-1:MXI-2),FH(-1:MXI-2)
C
      DIMENSION   FXG(-1:MAXIJP),     FXH(-1:MAXIJP)
      DIMENSION   FMG(-1:MAXIJP),     FMH(-1:MAXIJP)
      DIMENSION   FNG(-1:MAXIJP),     FNH(-1:MAXIJP)
      DIMENSION  CXG(-1:MAXIJP), ACXG(-1:MAXIJP)
      DIMENSION  CXH(-1:MAXIJP), ACXH(-1:MAXIJP)
      DIMENSION  DPG(-1:MAXIJP), DPH(-1:MAXIJP)
      DIMENSION DPFG(-1:MAXIJP), DPFH(-1:MAXIJP)
      DIMENSION SCXG(-1:MAXIJP), CCPG(-1:MAXIJP),  CCMG(-1:MAXIJP)
      DIMENSION SCXH(-1:MAXIJP), CCPH(-1:MAXIJP),  CCMH(-1:MAXIJP)
      DIMENSION   DDAG(-1:MAXIJP),    DDA2G(-1:MAXIJP)
      DIMENSION   DDAH(-1:MAXIJP),    DDA2H(-1:MAXIJP)
      DIMENSION   EAG(-1:MAXIJP),     EAH(-1:MAXIJP)
      DIMENSION    EG(-1:MAXIJP),      EH(-1:MAXIJP)
      DIMENSION  ADPG(-1:MAXIJP),    ADPH(-1:MAXIJP)
      DIMENSION   DMG(-1:MAXIJP),     DMH(-1:MAXIJP)
      DIMENSION  DMFG(-1:MAXIJP),    DMFH(-1:MAXIJP)
      DIMENSION  ADMG(-1:MAXIJP),    ADMH(-1:MAXIJP)
      DIMENSION   DAG(-1:MAXIJP),     DAH(-1:MAXIJP)
      DIMENSION   DBG(-1:MAXIJP),     DBH(-1:MAXIJP)
      DIMENSION DPDAG(-1:MAXIJP),   DPDAH(-1:MAXIJP)
      DIMENSION DPDBG(-1:MAXIJP),   DPDBH(-1:MAXIJP)
      DIMENSION DMDAG(-1:MAXIJP),   DMDAH(-1:MAXIJP)
      DIMENSION DMDBG(-1:MAXIJP),   DMDBH(-1:MAXIJP)
      DIMENSION    DG(-1:MAXIJP),      DH(-1:MAXIJP)
	DIMENSION    GG(-1:MAXIJP),      HH(-1:MAXIJP)
	DIMENSION  DFXG(-1:MAXIJP),    DFXH(-1:MAXIJP)
	DIMENSION  GG2(-1:MAXIJP,2),      HH2(-1:MAXIJP,2)
	DIMENSION GG1(-1:MAXIJP,3,2),   HH1(-1:MAXIJP,3,2)
	DIMENSION   FGP(-1:MAXIJP),     FHP(-1:MAXIJP)
	DIMENSION   FGM(-1:MAXIJP),         FHM(-1:MAXIJP)
C
      !EPS=0.01
      EPS = ENTEPS
      EPSS = EPS*EPS
      EPS2 = 2.*EPS
	epweno = 1.e-6
C
      DO  I = -1,NI2
C
      DO  J = -1,NJ2
      FXG(J)   = CET(I,J)*RG(I,J)
      FXH(J)   = CET(I,J)*RH(I,J)
      enddo
	do J = -1,NJ1 
		DFXG(J) = FXG(J+1) - FXG(J)
          DFXH(J) = FXH(J+1) - FXH(J)
		DG(J) = RG(I,J+1) - RG(I,J)
          DH(J) = RH(I,J+1) - RH(I,J)
	   GG2(J,1) = 0.5 * (DFXG(J) + abs(CET(I,J)) * DG(J))
	   GG2(J,2) = 0.5 * (-1.*DFXG(J) + abs(CET(I,J)) * DG(J))
         HH2(J,1) = 0.5 * (DFXH(J) + abs(CET(I,J)) * DH(J))
	   HH2(J,2) = 0.5 * (-1.*DFXH(J) + abs(CET(I,J)) * DH(J))	   
	end do
	   DFXG(NJ2)=DFXG(NJ1)
         DFXH(NJ2)=DFXH(NJ1)
	   DG(NJ2)=DG(NJ1)
         DH(NJ2)=DH(NJ1)
	   GG2(NJ2,1)=GG2(NJ1,1)
         HH2(NJ2,1)=HH2(NJ1,1)
	   GG2(NJ2,2)=GG2(NJ1,2)
         HH2(NJ2,2)=HH2(NJ1,2)

	do m1 = 1, 3
		k0 = m1 -  2		! (-1,  0,  1)
		k1 =  2 - m1		! ( 1,  0, -1)
			
	  do  J = 0, NJ1
		   GG1(J,m1,1) = GG2(J+k0,1)   
		   GG1(J,m1,2) = GG2(J+k1,2)
             HH1(J,m1,1) = HH2(J+k0,1)   
		   HH1(J,m1,2) = HH2(J+k1,2)
	  end do
	       GG1(NJ2,M1,1)=GG1(NJ1,M1,1)
	       GG1(NJ2,M1,2)=GG1(NJ1,M1,2)
             HH1(NJ2,M1,1)=HH1(NJ1,M1,1)
	       HH1(NJ2,M1,2)=HH1(NJ1,M1,2)
	       GG1(-1,M1,1)=GG1(0,M1,1)
	       GG1(-1,M1,2)=GG1(0,M1,2)
             HH1(-1,M1,1)=HH1(0,M1,1)
	       HH1(-1,M1,2)=HH1(0,M1,2)
	end do
		
	do J = -1, NJ2
		FNG(J) = 0.
	    FNH(J) = 0.
	end do
		
	do m1 = 1, 2
	  do J = 0, NJ
	if(m1 .eq. 1)then		
	TTTG1 =  1./(3.*( epweno + (GG1(J,1,1))**2 )**2)
	TTTG2 =  2./(3.*( epweno + (GG1(J,2,1))**2 )**2)
	TTTG=    TTTG1+TTTG2
	SG1=     TTTG1/TTTG
	SG2=     TTTG2/TTTG
	TTTH1 =  1./(3.*( epweno + (HH1(J,1,1))**2 )**2)
	TTTH2 =  2./(3.*( epweno + (HH1(J,2,1))**2 )**2)
	TTTH=    TTTH1+TTTH2
	SH1=     TTTH1/TTTH
	SH2=     TTTH2/TTTH	      

      fGpm1=0.5*(FXG(J-1)+abs(CET(I,J))*RG(i,J-1))
	fGpm0=0.5*(FXG(J)+abs(CET(I,J))*RG(i,J))
	fGpp1=0.5*(FXG(J+1)+abs(CET(I,J))*RG(i,J+1))
	eGp1=1./2.*(3.*fGpm0-1.*fGpm1)
	eGp2=1./2.*(fGpm0+fGpp1)

	FGP(J)=SG1*eGp1+SG2*eGp2
C
      fHpm1=0.5*(FXH(J-1)+abs(CET(I,J))*RH(i,J-1))
	fHpm0=0.5*(FXH(J)+abs(CET(I,J))*RH(i,J))
	fHpp1=0.5*(FXH(J+1)+abs(CET(I,J))*RH(i,J+1))
	eHp1=1./2.*(3.*fHpm0-1.*fHpm1)
	eHp2=1./2.*(fHpm0+fHpp1)

	FHP(J)=SH1*eHp1+SH2*eHp2
	else
	TTTG1 =  2./(3.*( epweno + (GG1(J,2,2))**2 )**2)
	TTTG2 =  1./(3.*( epweno + (GG1(J,1,2))**2 )**2)
	TTTG=    TTTG1+TTTG2
	SG1=     TTTG1/TTTG
	SG2=     TTTG2/TTTG
	TTTH1 =  2./(3.*( epweno + (HH1(J,2,2))**2 )**2)
	TTTH2 =  1./(3.*( epweno + (HH1(J,1,2))**2 )**2)
	TTTH=    TTTH1+TTTH2
	SH1=     TTTH1/TTTH
	SH2=     TTTH2/TTTH
C
	fGmm0=0.5*(FXG(J)-abs(CET(I,J))*RG(i,J))
	fGmp1=0.5*(FXG(J+1)-abs(CET(I,J))*RG(i,J+1))
	fGmp2=0.5*(FXG(J+2)-abs(CET(I,J))*RG(I,J+2))
	eGm1=1./2.*(fGmm0+1.*fGmp1)
	eGm2=1./2.*(3.*fGmp1-1.*fGmp2)

	FGM(J)=SG1*eGm1+SG2*eGm2
!	WRITE(*,*) FGM(3),FGM(20)
!	PAUSE
C
	fHmm0=0.5*(FXH(J)-abs(CET(I,J))*RH(i,J))
	fHmp1=0.5*(FXH(J+1)-abs(CET(I,J))*RH(i,J+1))
	fHmp2=0.5*(FXH(J+2)-abs(CET(I,J))*RH(i,J+2))
	eHm1=1./2.*(fHmm0+1.*fHmp1)
	eHm2=1./2.*(3.*fHmp1-1.*fHmp2)

	FHM(J)=SH1*eHm1+SH2*eHm2
	endif

	  end do
	end do
 
      do J=0,NJ
	FNG(J)=FGP(J)+FGM(J)
	FNH(J)=FHP(J)+FHM(J)
	end do
	FNG(-1)=FNG(0)
	FNH(-1)=FNH(0)
	FNG(NJ1)=FNG(NJ)
	FNH(NJ1)=FNH(NJ)
	FNG(NJ2)=FNG(NJ1)
	FNH(NJ2)=FNH(NJ1)

      DO  J = 0,NJ2
      DQG(I,J) = DQG(I,J) - DTLOC(I,J)*(FNG(J) - FNG(J-1))
      DQH(I,J) = DQH(I,J) - DTLOC(I,J)*(FNH(J) - FNH(J-1))
      enddo
      enddo
C
      RETURN
      END