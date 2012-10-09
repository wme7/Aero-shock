       SUBROUTINE VEL2D(XR,XM,XN,XE,R,E,U,V,P)
       COMMON /FINP3/ ALPHA,UINF
       COMMON /CNST3/ GAMMA,GM1,RGM1,RGM2,PI,SPI,UTX
      VX   = XM/XE
      XA   = XN/XM
      VY   = XA*VX
      UV2  = VX*VX + VY*VY
      IF( UV2 .GE. 1.0 )THEN
        WRITE(6,*)'U,V =',VX,VY
        WRITE(6,*)'Error AT FUNCTION VEL2D'
        IF( VX .GE. 0. )THEN 
         VX = UINF
        ELSE
         VX =-UINF
        ENDIF
       VY = XA*VX
      ENDIF
      RLGA = 1./SQRT(1.-VX*VX-VY*VY)
      DA   = XR/RLGA
      EA   = XE-XM*VX-XN*VY
      PA   = GM1*(EA-DA)
      XF=DA*DA*(XE+PA)-XR*XR*(EA+PA)
      XG=XA*VX*(XE+PA)-XN
c
      IF( XF .LE. 1.E-20 .AND. XG .LE. 1.E-20 ) THEN
        U = VX
        V = XA*VX
        R = DA
        E = XE-(XM+XA*XN)*VX
        P = GM1*(E-R)
      ELSE
       ITER = 0
  10   CONTINUE
       ITER = ITER+1
       IF( ITER .GT. 100 ) GOTO 20
  5    CONTINUE
       UV2  = VX*VX + VY*VY
       IF( UV2 .GE. 1 ) THEN
         VX=VX/2.
         VY=VY/2.
         GOTO 5
       ENDIF
       RLGA = 1./SQRT(1.-VX*VX-VY*VY)
       DA   = XR/RLGA
       EA   = XE-(XM+XA*XN)*VX
       PA   = GM1*(EA-DA)
       VY   = XA*VX
       XEU = -(XM+XA*XN)
       XFU = (DA*DA*GM1-XR*XR*GAMMA)*XEU
       XFR = 2.*DA*(XE+PA)-(DA*DA-XR*XR)*GM1
       XGU = XA*(XE+PA)+XA*VX*GM1*XEU
       XGR = -XA*VX*GM1
       DELTA = XFU*XGR-XFR*XGU
       XD1 =  XGR/DELTA
       XD2 = -XGU/DELTA
       XD3 = -XFR/DELTA
       XD4 =  XFU/DELTA
       XF=DA*DA*(XE+PA)-XR*XR*(EA+PA)
       XG=XA*VX*(XE+PA)-XN
       IF( XF .LE. 1.E-20 .AND. XG .LE. 1.E-20 ) GOTO 20
       UNEW = VX-(XD1*XF+XD3*XG)
       DNEW = DA-(XD2*XF+XD4*XG)
       UOLD = VX
       DOLD = DA
       VX = UNEW
       DA = DNEW
       DIS = (UOLD-UNEW)**2+(DOLD-DNEW)**2
       ERR = SQRT(ABS(DIS))
       VY   = XA*VX
       UV2  = VX*VX + VY*VY
       IF( UV2 .GE. 1 ) GOTO 5
       IF( ERR .GT. 1.E-10) GOTO 10
 20    CONTINUE
        U = UNEW
        V = XA*UNEW
        R = DNEW
        E = XE-(XM+XA*XN)*U
        P = GM1*(E-R)
        ENDIF
      RETURN
      END 
C
       
