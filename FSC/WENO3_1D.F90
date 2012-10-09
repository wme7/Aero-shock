PROGRAM MAIN
IMPLICIT REAL(8) (A-H,O-Z)
REAL, DIMENSION (1000,1000) :: F,FEQ, FFF
REAL, DIMENSION (1000) :: X, C, V, R, RI,P, U, T, TI, ET, EI, FN, FA, FB, FC, AV
REAL, DIMENSION (20)   :: GH, W
REAL, DIMENSION (1000) :: Z, THETA, PSI, VIS,G1
OPEN (UNIT = 10, FILE = 'POSTITER.TEC', STATUS = 'UNKNOWN')
write(10,*) 'variables = "x","r","p","z","t" '
OPEN (UNIT = 20, FILE = 'INITCHECK.TEC', STATUS = 'UNKNOWN')
NX          = 101
NXP1        = NX
GHNC        = 1
CFL         = 0.8
OUTTIME     = 0.1
IT       = 0
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
    V1  = -20
    V2  = 20
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
        !Feq(K,J)   = 1/((EXP(PP)/Z(J)) + IT )
        F(K,J)   = 1/((EXP(PP)/Z(J)) + IT )
    END DO    
    END DO

    eps = 1E-25
    c1 = - 1./6.
    c2 = 2./6.
    c3 = 5./6.
    c4 = - 7./6.
    c5 = 11./6. 
    
    DO K = 1, nv  
             vxp = max(v(k),0.)
             vxm = min(v(k),0.)             
             do j = 4, nxp1-3
             sl0m = vxm*(13./12.)*(f(k,j-2)-2*f(k,j-1)+f(k,j))**2 + vxm*(1./4.)*(f(k,j-2)-4.*f(k,j-1)+3.*f(k,j))**2
             sl1m = vxm*(13./12.)*(f(k,j-1)-2*f(k,j)+f(k,j+1))**2 + vxm*(1./4.)*(f(k,j-1)-f(k,j+1))**2
             sl2m = vxm*(13./12.)*(f(k,j)-2*f(k,j+1)+f(k,j+2))**2 + vxm*(1./4.)*(3.*f(k,j)-4.*f(k,j+1)+f(k,j+2))**2
             sl0p = vxp*(13./12.)*(f(k,j-3)-2*f(k,j-2)+f(k,j-1))**2 + vxp*(1./4.)*(f(k,j-3)-4.*f(k,j-2)+3.*f(k,j-1))**2
             sl1p = vxp*(13./12.)*(f(k,j-2)-2*f(k,j-1)+f(k,j))**2 + vxp*(1./4.)*(f(k,j-2)-f(k,j))**2
             sl2p = vxp*(13./12.)*(f(k,j-1)-2*f(k,j)+f(k,j+1))**2 + vxp*(1./4.)*(3.*f(k,j-1)-4.*f(k,j)+f(k,j+1))**2
             
             sr0m = vxm*(13./12.)*(f(k,j-1)-2*f(k,j)+f(k,j+1))**2 + vxm*(1./4.)*(f(k,j-1)-4.*f(k,j)+3.*f(k,j+1))**2
             sr1m = vxm*(13./12.)*(f(k,j)-2*f(k,j+1)+f(k,j+2))**2 + vxm*(1./4.)*(f(k,j)-f(k,j+2))**2
             sr2m = vxm*(13./12.)*(f(k,j+1)-2*f(k,j+2)+f(k,j+3))**2 + vxm*(1./4.)*(3.*f(k,j+1)-4.*f(k,j+2)+f(k,j+3))**2
             sr0p = vxp*(13./12.)*(f(k,j-2)-2*f(k,j-1)+f(k,j))**2 + vxp*(1./4.)*(f(k,j-2)-4.*f(k,j-1)+3.*f(k,j))**2
             sr1p = vxp*(13./12.)*(f(k,j-1)-2*f(k,j)+f(k,j+1))**2 + vxp*(1./4.)*(f(k,j-1)-f(k,j+1))**2
             sr2p = vxp*(13./12.)*(f(k,j)-2*f(k,j+1)+f(k,j+2))**2 + vxp*(1./4.)*(3.*f(k,j)-4.*f(k,j+1)+f(k,j+2))**2
             
             al0m  = 1. / (10. * (eps + sl0m))**2
             al1m  = 6. / (10. * (eps + sl1m))**2
             al2m  = 3. / (10. * (eps + sl2m))**2
             al0p  = 1. / (10. * (eps + sl0p))**2
             al1p  = 6. / (10. * (eps + sl1p))**2
             al2p  = 3. / (10. * (eps + sl2p))**2
             
             ar0m  = 3. / (10. * (eps + sr0m))**2
             ar1m  = 6. / (10. * (eps + sr1m))**2
             ar2m  = 1. / (10. * (eps + sr2m))**2
             ar0p  = 3. / (10. * (eps + sr0p))**2
             ar1p  = 6. / (10. * (eps + sr1p))**2
             ar2p  = 1. / (10. * (eps + sr2p))**2
             
             wl0m  = al0m / (al0m+al1m+al2m)
             wl1m  = al1m / (al0m+al1m+al2m)
             wl2m  = al2m / (al0m+al1m+al2m)
             wl0p  = al0p / (al0p+al1p+al2p)
             wl1p  = al1p / (al0p+al1p+al2p)
             wl2p  = al2p / (al0p+al1p+al2p)
             
             wr0m  = ar0m / (ar0m+ar1m+ar2m)
             wr1m  = ar1m / (ar0m+ar1m+ar2m)
             wr2m  = ar2m / (ar0m+ar1m+ar2m)
             wr0p  = ar0p / (ar0p+ar1p+ar2p)
             wr1p  = ar1p / (ar0p+ar1p+ar2p)
             wr2p  = ar2p / (ar0p+ar1p+ar2p)
             
             flm=vxm * (wl0m*(c1*f(k,j-2)+c3*f(k,j-1)+c2*f(k,j))+wl1m*(c2*f(k,j-1)+c3*f(k,j)+c1*f(k,j+1))+wl2m*(c5*f(k,j)+c4*f(k,j+1)+c2*f(k,j+2))) 
             flp=vxp * (wl0p*(c2*f(k,j-3)+c4*f(k,j-2)+c5*f(k,j-1))+wl1p*(c1*f(k,j-2)+c3*f(k,j-1)+c2*f(k,j))+wl2p*(c2*f(k,j-1)+c3*f(k,j)+c1*f(k,j+1)))
             frm=vxm * (wr0m*(c1*f(k,j-1)+c3*f(k,j)+c2*f(k,j+1))+wr1m*(c2*f(k,j)+c3*f(k,j+1)+c1*f(k,j+2))+wr2m*(c5*f(k,j+1)+c4*f(k,j+2)+c2*f(k,j+3)))
             frp=vxp * (wr0p*(c2*f(k,j-2)+c4*f(k,j-1)+c5*f(k,j))+wr1p*(c1*f(k,j-1)+c3*f(k,j)+c2*f(k,j+1))+wr2p*(c2*f(k,j)+c3*f(k,j+1)+c1*f(k,j+2)))
             
             fl = flm + flp
             fr = frm + frp
           
             fn(j)  =  f(k,j) - dt/dx * (fr - fl)
             END DO 

            DO j = 1,NXP1   
            F(K,J) = FN(J)     
                  
            F(K,1) = fn(4)
            F(K,2) = fn(4) 
            F(k,3) = fn(4)
            F(K,NXP1) = fn(nxp1-3)      
            F(K,NXP1-1) = fn(nxp1-3) 
            F(K,NXP1-2) = fn(nxp1-3)                  
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
    WRITE (10,*) X(J), R(J), p(j), z(j), t(j)
    END DO
STOP    
END PROGRAM




