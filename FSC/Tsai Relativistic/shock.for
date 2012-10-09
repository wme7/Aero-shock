      PROGRAM SHOCK_STATE
      WRITE(*,*)' INPUT THE LEFT_STATE :'
      WRITE(*,*) 'V1,P1,D1,XM'
      READ(*,*) V1,P1,D1,XM
      WRITE(*,*) 'GAMMA'
      READ(*,*) GAMMA
      ERROR = 1.E-8
      write(*,*) 'input initial D0, P0'
       read(*,*) D0,P0
C***** LEFT_STATE *****
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
      write(*,*) 'A1 ='
      write(*,*) A1
      US = XM*A1
      write(*,*) 'US = '
      write(*,*) US
c      P0 = P1+1.
c      D0 = D1+0.5
C
      T1 = GMI*P0+D0-E1
      T2 = E1+P0
      T3 = US*US
      U1 = GMT*P0+D0
      U2 = P0-P1
      U3 = GMI*P0+D0+P1
      W1 = S1*D0*D0
      W2 = S2*D0*D0
C
       XF = U2*(W1+U1)+W2-U1*U1
       XG = T3*T1*T2-U2*U3
      IF( ABS(XF) .LT. ERROR .AND. ABS(XG) .LT. ERROR) 
     1 GOTO 20          
      ITER = 0
 10   CONTINUE
      ITER = ITER + 1
      IF( ITER .GT. 100 ) GOTO 20
C 
      XFP = (W1+U1)+U2*GMT-2.*GMT*U1
      XFD = U2*(2.*S1*D0+1.)+2.*S2*D0-2.*U1
      XGP = T3*GMI*T2+T3*T1-U3-U2*GMI
      XGD = T3*T2-U2
C***** DETERMINANT *****
       xdel = XFP*XGD-XFD*XGP
C       write(*,*) 'determinant : Jacobia'
C       write(*,*) xdel
      DELI = 1./xdel
       XD1 =  XGD*DELI
       XD2 = -XGP*DELI
       XD3 = -XFD*DELI
       XD4 =  XFP*DELI
       PNEW = P0-(XD1*XF+XD3*XG) 
       DNEW = D0-(XD2*XF+XD4*XG) 
       POLD = P0
       DOLD = D0
       P0 = PNEW
       D0 = DNEW
C 
      T1 = GMI*P0+D0-E1
      T2 = E1+P0
      T3 = US*US
      U1 = GMT*P0+D0
      U2 = P0-P1
      U3 = GMI*P0+D0+P1
      W1 = S1*D0*D0
      W2 = S2*D0*D0
C
       XF = U2*(W1+U1)+W2-U1*U1
       XG = T3*T1*T2-U2*U3
c       write(*,*) 'xf,xg,p0,d0'
c       write(*,*) xf,xg,p0,d0
      IF( ABS(XF) .LT. ERROR .AND. ABS(XG) .LT. ERROR) 
     1 GOTO 20   
      GOTO 10
 20   CONTINUE
      write(*,*) iter
      E = GMI*P0+D0
      U1=SQRT((P0-P1)*(E+P1)/(E-E1)/(E1+P0))
      U2=SQRT((P0-P1)*(E1+P0)/(E-E1)/(E+P1))
      V = (U2-U1)/(1.-U1*U2)
      IF( ABS(V) .GT. 1. ) THEN
        WRITE(*,*) '  VELOCITY ERROR !'
        STOP
      ENDIF
      RGAM = 1./SQRT(1.-U2*U2)
      EP1 = E1+P1
      EP2 = E + P0
      RGAM1 = 1./SQRT(1.-U1*U1)
      ERR1 = D1*RGAM1*U1 - D0*RGAM*U2
      ERR2 = EP1*RGAM1**2*U1**2+P1-
     1      ( EP2*RGAM**2*U2**2+P0 ) 
      ERR3 = EP1*RGAM1**2*U1-EP2*RGAM**2*U2    
      WRITE(*,11) V,P0,D0,E
 11   FORMAT('V = ',F13.9,2X,'P = ',F13.9,2X,'D = '
     1      ,F13.9,2X,'E = ',F13.9)    
      WRITE(*,12) ERR1,ERR2,ERR3
 12   FORMAT('ERR1 = ',F13.9,2X,'ERR2 = ',F13.9,2X 
     1       ,'ERR3 = ',F13.9)  
C      ENDIF
      STOP
      END  
