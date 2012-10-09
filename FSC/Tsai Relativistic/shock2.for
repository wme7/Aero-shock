      PROGRAM SHOCK_STATE
      COMMON /VAR/ T1,T2,T3,FD,EAP,GMI,GMT
      COMMON /WAR/ S1,S2,P1,R1,R2,X1
      WRITE(*,*)' INPUT THE LEFT_STATE :'
      WRITE(*,*) 'V1,P1,D1,XM'
      READ(*,*) V1,P1,D1,XM
      WRITE(*,*) 'GAMMA'
      READ(*,*) GAMMA
      ERROR = 1.E-10
      write(*,*) 'input initial  P0'
       read(*,*) PS
C***** LEFT_STATE *****
      i = 0
      GM1 = GAMMA - 1.
      GMI = 1./GM1
      GMT = GAMMA*GMI
      E1 = GMI*P1+D1
      DW = D1*D1
      EAP = E1+P1
      S1 = EAP/DW
      S2 = EAP*EAP/DW
      write(*,*) 'E1,S1,S2'
      write(*,*) e1,s1,s2
C***** INITIAL POINT GAUSS *****
      A1W = (GAMMA*P1)/(D1+GMT*P1)
      A1 = SQRT(A1W)
      US = XM*A1
      write(*,13) A1,US
 13   format('A1 = ',F13.8,1X,'US = ',F13.8)
C
  5   continue
      T1 = US*US
      T2 = T1 - 1.
      T3 = P1 + T1*E1
      FD = T3/T2
C
      R1 = ROE(P0)
      R2 = R1*R1
      X1 = GMT*P0+R1
C
       XF = XFUN(P0)
      IF( ABS(XF) .LT. ERROR ) GOTO 20
      ITER = 0
 10   CONTINUE
      ITER = ITER + 1
      IF( ITER .GT. 500 ) GOTO 20
C
       DXF = XDFUN(P0)
       IF( ABS(DXF) .LT. 1.E-19 ) STOP
       PNEW = P0-(XF/DXF)
       P0 = PNEW
C
C     T1 = US*US
C     T2 = T1 - 1.
C     T3 = P1 + T1*E1
C     FD = T3/T2
C
      R1 = ROE(P0)
      R2 = R1*R1
      X1 = GMT*P0+R1
C
       XF = XFUN(P0)
      IF( ABS(XF) .LT. ERROR ) GOTO 20
C 
      GOTO 10
 20   CONTINUE
      D0 = R1
      if( (P0-P1) .lt. error .or. D0 .LT. 0.) then
c         write(*,*) 'D0=D1, P0=P1'
         i = i+1
         if( i .gt. 10000 ) stop
         P0 = PS + i*0.3
         pinit = p0
         goto 5
      endif
         write(*,2) pinit
 2       format('p_init = ' f13.8)
      write(6,14) iter
 14   format('ITER = ',I5)
         write(*,*) ' result :'
      E = GMI*P0+D0
      U1=SQRT((P0-P1)*(E+P1)/(E-E1)/(E1+P0))
      U2=SQRT((P0-P1)*(E1+P0)/(E-E1)/(E+P1))

      WRITE(*,17) U1,U2
 17   FORMAT('U1 = ',F13.9,2X,'U2 = ',F13.9)
      V = (U2-U1)/(1.-U1*U2)
      IF( ABS(V) .GT. 1. ) THEN
        WRITE(*,*) '  VELOCITY ERROR !'
        STOP
      ENDIF
      WRITE(*,11) V,P0,D0,E
 11   FORMAT('V = ',F13.9,2X,'P = ',F13.9,2X,'D = '
     1      ,F13.9,2X,'E = ',F15.9)
      RGAM = 1./SQRT(1.-U2*U2)
      EP1 = E1+P1
      EP2 = E + P0
      RGAM1 = 1./SQRT(1.-U1*U1)
      ERR1 = D1*RGAM1*U1 - D0*RGAM*U2
      ERR2 = EP1*RGAM1**2*U1**2+P1-
     1      ( EP2*RGAM**2*U2**2+P0 )
      ERR3 = EP1*RGAM1**2*U1-EP2*RGAM**2*U2
      WRITE(*,12) ERR1,ERR2,ERR3
 12   FORMAT('ERR1 = ',E17.9,2X,'ERR2 = ',E17.9,2X
     1       ,'ERR3 = ',E17.9)
C      ENDIF
      STOP
      END
C
      FUNCTION ROE(P)
      COMMON /VAR/ T1,T2,T3,FD,EAP,GMI,GMT
      XR = T1/T2*EAP*EAP
      XDF = -T2*P-T3
      IF( ABS(XDF) .LT. 1.E-18 ) STOP
      ROE = -GMI*P+FD+XR/XDF
      RETURN
      END
C
      FUNCTION XROE(P)
      COMMON /VAR/ T1,T2,T3,FD,EAP,GMI,GMT
      XR = T1/T2*EAP*EAP
      XDF = (-T2*P-T3)**2
      IF( ABS(XDF) .LT. 1.E-18 ) STOP
      XROE = -GMI+T2*XR/XDF
      RETURN
      END
C
      FUNCTION XFUN(P)
      COMMON /WAR/ S1,S2,P1,R1,R2,X1
      XFUN = (P-P1)*(R2*S1+X1)-X1*X1+S2*R2
      RETURN
      END
C
      FUNCTION XDFUN(P)
      REAL P
      COMMON /VAR/ T1,T2,T3,FD,EAP,GMI,GMT
      COMMON /WAR/ S1,S2,P1,R1,R2,X1
      A1 = R2*S1 + X1
      DR1 = XROE(P)
      A2 = (P-P1)*(2.*R1*DR1*S1+GMT+DR1)
      A3 = -2.*X1*(GMT+DR1)+2.*S2*R1*DR1
      XDFUN = A1+A2+A3
      RETURN
      END
