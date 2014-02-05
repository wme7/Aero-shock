PROGRAM test_app_2
IMPLICIT REAL(8) (A-H,O-Z)
REAL, DIMENSION (20,20,400,400)     :: f
REAL, DIMENSION (400,400)           :: ux,uy,z,t,r,et,p,fn,tha,thb,phia,phib
REAL, DIMENSION (400,400)           :: fll,flh,fl,fhl,fhh,fh,gll,glh,gl,ghl,ghh,gh2
REAL, DIMENSION (400)               :: x,y
REAL, DIMENSION (20)                :: c,v,gh,w

    NX          = 50
    NY          = 50
    CFL         = 0.4
    OUTTIME     = 0.2
    sv_point1 = (1./3.)*outtime   
    sv_point2 = (2./3.)*outtime 
    theta       = 1
    ! MB = 0, FD = 1, BE = -1
    NV = 20 
    PI = ATAN2(1.,1.)*4.

OPEN (UNIT = 10, FILE = 'abscissas.TEC', STATUS = 'UNKNOWN')
OPEN (UNIT = 20, FILE = 'sv_point(1).TEC', STATUS = 'UNKNOWN')
OPEN (UNIT = 30, FILE = 'sv_point(2).TEC', STATUS = 'UNKNOWN')
    write(20,888)
    write(20,*) 'zone T="sv_point(1)",I=',nx,',J=',ny
    write(30,888)
    write(30,*) 'zone T="sv_point(2)",I=',nx,',J=',ny
OPEN (UNIT = 60, FILE = 'results.TEC', STATUS = 'UNKNOWN')
    write(60,888)
    write(60,*) 'zone T="res_2d",I=',nx,',J=',ny
OPEN (UNIT = 90, FILE = 'f_initial.TEC', STATUS = 'UNKNOWN')
OPEN (UNIT = 100, FILE = 'f_time.TEC', STATUS = 'UNKNOWN')
    write(90,*)'variables = "Vx","Vy","fq1","fq2","fq3","fq4"'       
    write(90,*)'zone T="f_2d",I=',nv,',J=',nv
    write(100,*) 'variables = "Vx","Vy","fq1","fq2","fq3","fq4"'
    write(100,*)'zone T="f_2d",I=',nv,',J=',nv    
888 format ('variables = "x","y","n","p","z"')    

!
! GAUSS-HERMITE QUADRATURE
!
GH(:)=(/-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738, & 
	-2.78880605843,-2.25497400209,-1.73853771212,-1.2340762154, &
	-0.737473728545,-0.245340708301,0.245340708301,0.737473728545, &
	1.2340762154,1.73853771212,2.25497400209,2.78880605843, &
	3.34785456738,3.94476404012,4.60368244955,5.38748089001/)
W(:) =(/0.898591961453,0.704332961176,0.62227869619,0.575262442852, & 
	0.544851742366,0.524080350949,0.509679027117,0.499920871336, &
	0.493843385272,0.490921500667,0.490921500667,0.493843385272, &
	0.499920871336,0.509679027117,0.524080350949,0.544851742366, &
	0.575262442852,0.62227869619,0.704332961176,0.898591961453/)
    DO I = 1, NV
        C(I)    = W(I)
        V(I)    = GH(I)
        WRITE (10,*) C(I),V(I)    
    END DO

    Z1   = 0.142
    UX1  = -0.75
    UY1  = -0.5
    T1   = 2.078
    
    Z2   = 0.4253
    UX2  = -0.75
    UY2  = 0.5
    T2   = 1.1494
    
    Z3   = 0.142
    UX3  = 0.75
    UY3  = 0.5
    T3   = 2.078
    
    Z4   = 0.6635
    UX4  = 0.75
    UY4  = -0.5
    T4   = 0.87685
    
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

DO I = 1,NX
DO J = 1,NY    
           if ((Y(J) .lt. 0.5) .and. (x(i) .gt. 0.5))  then
                Z(I,J) = Z4
                UX(I,J)= UX4
                UY(I,J)= UY4
                T(I,J) = T4   
            else if ((Y(J) .gt. 0.5) .and. (x(i) .lt. 0.5)) then
                Z(I,J) = Z2
                UX(I,J)= UX2
                UY(I,J)= UY2
                T(I,J) = T2
            else if ((Y(J) .ge. 0.5) .and. (x(i) .ge. 0.5)) then
                Z(I,J) = Z1
                UX(I,J)= UX1
                UY(I,J)= UY1
                T(I,J) = T1
            else
                Z(I,J) = Z3
                UX(I,J)= UX3
                UY(I,J)= UY3
                T(I,J) = T3    
            end if                
               
    DO K = 1, NV
    DO L = 1, NV
        PP = ( (V(K)-UX(I,J))**2 + (V(L)-UY(I,J))**2 ) / T(I,J)
        F(k,l,i,j) = 1/((EXP(PP)/Z(I,J)) + theta)
    END DO
    END DO    
END DO
END DO

do k = 1,nv
do l = 1,nv
    write (90,*) v(k),v(l),f(k,l,nx*3/4,ny*3/4),f(k,l,nx/4,ny*3/4),f(k,l,nx/4,ny/4),f(k,l,nx*3/4,ny/4)
end do
end do

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

DO I = 1, NX
DO J = 1, NY
    DO K = 1, NV
    DO L = 1, NV
        PP = ( (V(K)-UX(I,J))**2 + (V(L)-UY(I,J))**2 ) / T(I,J)
        F(k,l,i,j) = 1/((EXP(PP)/Z(I,J))+theta)
    END DO
    END DO  
END DO
END DO

    eps = 1E-10
    c1 = - 1./6.
    c2 = 2./6.
    c3 = 5./6.
    c4 = - 7./6.
    c5 = 11./6. 

DO K = 1, NV
DO L = 1, NV
    vxp = max(v(k),0.)
    vxm = min(v(k),0.) 
    vyp = max(v(l),0.)
    vym = min(v(l),0.) 
    
    DO I = 4, NX-3
    DO J = 4, NY-3     
             sl0mx = vxm*(13./12.)*(f(k,l,i-2,j)-2*f(k,l,i-1,j)+f(k,l,i,j))**2 &
		+ vxm*(1./4.)*(f(k,l,i-2,j)-4.*f(k,l,i-1,j)+3.*f(k,l,i,j))**2
             sl1mx = vxm*(13./12.)*(f(k,l,i-1,j)-2*f(k,l,i,j)+f(k,l,i+1,j))**2 &
		+ vxm*(1./4.)*(f(k,l,i-1,j)-f(k,l,i+1,j))**2
             sl2mx = vxm*(13./12.)*(f(k,l,i,j)-2*f(k,l,i+1,j)+f(k,l,i+2,j))**2 &
		+ vxm*(1./4.)*(3.*f(k,l,i,j)-4.*f(k,l,i+1,j)+f(k,l,i+2,j))**2
             sl0px = vxp*(13./12.)*(f(k,l,i-3,j)-2*f(k,l,i-2,j)+f(k,l,i-1,j))**2 &
		+ vxp*(1./4.)*(f(k,l,i-3,j)-4.*f(k,l,i-2,j)+3.*f(k,l,i-1,j))**2
             sl1px = vxp*(13./12.)*(f(k,l,i-2,j)-2*f(k,l,i-1,j)+f(k,l,i,j))**2 &
		+ vxp*(1./4.)*(f(k,l,i-2,j)-f(k,l,i,j))**2
             sl2px = vxp*(13./12.)*(f(k,l,i-1,j)-2*f(k,l,i,j)+f(k,l,i+1,j))**2 &
		+ vxp*(1./4.)*(3.*f(k,l,i-1,j)-4.*f(k,l,i,j)+f(k,l,i+1,j))**2
             sr0mx = vxm*(13./12.)*(f(k,l,i-1,j)-2*f(k,l,i,j)+f(k,l,i+1,j))**2 &
		+ vxm*(1./4.)*(f(k,l,i-1,j)-4.*f(k,l,i,j)+3.*f(k,l,i+1,j))**2
             sr1mx = vxm*(13./12.)*(f(k,l,i,j)-2*f(k,l,i+1,j)+f(k,l,i+2,j))**2 &
		+ vxm*(1./4.)*(f(k,l,i,j)-f(k,l,i+2,j))**2
             sr2mx = vxm*(13./12.)*(f(k,l,i+1,j)-2*f(k,l,i+2,j)+f(k,l,i+3,j))**2 &
		+ vxm*(1./4.)*(3.*f(k,l,i+1,j)-4.*f(k,l,i+2,j)+f(k,l,i+3,j))**2
             sr0px = vxp*(13./12.)*(f(k,l,i-2,j)-2*f(k,l,i-1,j)+f(k,l,i,j))**2 &
		+ vxp*(1./4.)*(f(k,l,i-2,j)-4.*f(k,l,i-1,j)+3.*f(k,l,i,j))**2
             sr1px = vxp*(13./12.)*(f(k,l,i-1,j)-2*f(k,l,i,j)+f(k,l,i+1,j))**2 &
		+ vxp*(1./4.)*(f(k,l,i-1,j)-f(k,l,i+1,j))**2
             sr2px = vxp*(13./12.)*(f(k,l,i,j)-2*f(k,l,i+1,j)+f(k,l,i+2,j))**2 &
		+ vxp*(1./4.)*(3.*f(k,l,i,j)-4.*f(k,l,i+1,j)+f(k,l,i+2,j))**2
             
             sl0my = vym*(13./12.)*(f(k,l,i,j-2)-2*f(k,l,i,j-1)+f(k,l,i,j))**2 &
		+ vym*(1./4.)*(f(k,l,i,j-2)-4.*f(k,l,i,j-1)+3.*f(k,l,i,j))**2
             sl1my = vym*(13./12.)*(f(k,l,i,j-1)-2*f(k,l,i,j)+f(k,l,i,j+1))**2 &
		+ vym*(1./4.)*(f(k,l,i,j-1)-f(k,l,i,j+1))**2
             sl2my = vym*(13./12.)*(f(k,l,i,j)-2*f(k,l,i,j+1)+f(k,l,i,j+2))**2 &
		+ vym*(1./4.)*(3.*f(k,l,i,j)-4.*f(k,l,i,j+1)+f(k,l,i,j+2))**2
             sl0py = vyp*(13./12.)*(f(k,l,i,j-3)-2*f(k,l,i,j-2)+f(k,l,i,j-1))**2 &
		+ vyp*(1./4.)*(f(k,l,i,j-3)-4.*f(k,l,i,j-2)+3.*f(k,l,i,j-1))**2
             sl1py = vyp*(13./12.)*(f(k,l,i,j-2)-2*f(k,l,i,j-1)+f(k,l,i,j))**2 &
		+ vyp*(1./4.)*(f(k,l,i,j-2)-f(k,l,i,j))**2
             sl2py = vyp*(13./12.)*(f(k,l,i,j-1)-2*f(k,l,i,j)+f(k,l,i,j+1))**2 &
		+ vyp*(1./4.)*(3.*f(k,l,i,j-1)-4.*f(k,l,i,j)+f(k,l,i,j+1))**2
             sr0my = vym*(13./12.)*(f(k,l,i,j-1)-2*f(k,l,i,j)+f(k,l,i,j+1))**2 &
		+ vym*(1./4.)*(f(k,l,i,j-1)-4.*f(k,l,i,j)+3.*f(k,l,i,j+1))**2
             sr1my = vym*(13./12.)*(f(k,l,i,j)-2*f(k,l,i,j+1)+f(k,l,i,j+2))**2 &
		+ vym*(1./4.)*(f(k,l,i,j)-f(k,l,i,j+2))**2
             sr2my = vym*(13./12.)*(f(k,l,i,j+1)-2*f(k,l,i,j+2)+f(k,l,i,j+3))**2 &
		+ vym*(1./4.)*(3.*f(k,l,i,j+1)-4.*f(k,l,i,j+2)+f(k,l,i,j+3))**2
             sr0py = vyp*(13./12.)*(f(k,l,i,j-2)-2*f(k,l,i,j-1)+f(k,l,i,j))**2 &
		+ vyp*(1./4.)*(f(k,l,i,j-2)-4.*f(k,l,i,j-1)+3.*f(k,l,i,j))**2
             sr1py = vyp*(13./12.)*(f(k,l,i,j-1)-2*f(k,l,i,j)+f(k,l,i,j+1))**2 &
		+ vyp*(1./4.)*(f(k,l,i,j-1)-f(k,l,i,j+1))**2
             sr2py = vyp*(13./12.)*(f(k,l,i,j)-2*f(k,l,i,j+1)+f(k,l,i,j+2))**2 &
		+ vyp*(1./4.)*(3.*f(k,l,i,j)-4.*f(k,l,i,j+1)+f(k,l,i,j+2))**2
             
             al0mx  = 1. / (10. * (eps + sl0mx))**2
             al1mx  = 6. / (10. * (eps + sl1mx))**2
             al2mx  = 3. / (10. * (eps + sl2mx))**2
             al0px  = 1. / (10. * (eps + sl0px))**2
             al1px  = 6. / (10. * (eps + sl1px))**2
             al2px  = 3. / (10. * (eps + sl2px))**2
             ar0mx  = 3. / (10. * (eps + sr0mx))**2
             ar1mx  = 6. / (10. * (eps + sr1mx))**2
             ar2mx  = 1. / (10. * (eps + sr2mx))**2
             ar0px  = 3. / (10. * (eps + sr0px))**2
             ar1px  = 6. / (10. * (eps + sr1px))**2
             ar2px  = 1. / (10. * (eps + sr2px))**2
             
             al0my  = 1. / (10. * (eps + sl0my))**2
             al1my  = 6. / (10. * (eps + sl1my))**2
             al2my  = 3. / (10. * (eps + sl2my))**2
             al0py  = 1. / (10. * (eps + sl0py))**2
             al1py  = 6. / (10. * (eps + sl1py))**2
             al2py  = 3. / (10. * (eps + sl2py))**2
             ar0my  = 3. / (10. * (eps + sr0my))**2
             ar1my  = 6. / (10. * (eps + sr1my))**2
             ar2my  = 1. / (10. * (eps + sr2my))**2
             ar0py  = 3. / (10. * (eps + sr0py))**2
             ar1py  = 6. / (10. * (eps + sr1py))**2
             ar2py  = 1. / (10. * (eps + sr2py))**2
             !weightings x             
             wl0mx  = al0mx / (al0mx+al1mx+al2mx)
             wl1mx  = al1mx / (al0mx+al1mx+al2mx)
             wl2mx  = al2mx / (al0mx+al1mx+al2mx)
             wl0px  = al0px / (al0px+al1px+al2px)
             wl1px  = al1px / (al0px+al1px+al2px)
             wl2px  = al2px / (al0px+al1px+al2px)
             wr0mx  = ar0mx / (ar0mx+ar1mx+ar2mx)
             wr1mx  = ar1mx / (ar0mx+ar1mx+ar2mx)
             wr2mx  = ar2mx / (ar0mx+ar1mx+ar2mx)
             wr0px  = ar0px / (ar0px+ar1px+ar2px)
             wr1px  = ar1px / (ar0px+ar1px+ar2px)
             wr2px  = ar2px / (ar0px+ar1px+ar2px)
             !weightings y
             wl0my  = al0my / (al0my+al1my+al2my)
             wl1my  = al1my / (al0my+al1my+al2my)
             wl2my  = al2my / (al0my+al1my+al2my)
             wl0py  = al0py / (al0py+al1py+al2py)
             wl1py  = al1py / (al0py+al1py+al2py)
             wl2py  = al2py / (al0py+al1py+al2py)
             wr0my  = ar0my / (ar0my+ar1my+ar2my)
             wr1my  = ar1my / (ar0my+ar1my+ar2my)
             wr2my  = ar2my / (ar0my+ar1my+ar2my)
             wr0py  = ar0py / (ar0py+ar1py+ar2py)
             wr1py  = ar1py / (ar0py+ar1py+ar2py)
             wr2py  = ar2py / (ar0py+ar1py+ar2py) 
             !negative & positive fluxes x
             flmx=vxm * (wl0mx*(c1*f(k,l,i-2,j)+c3*f(k,l,i-1,j)+c2*f(k,l,i,j)) &
		+wl1mx*(c2*f(k,l,i-1,j)+c3*f(k,l,i,j)+c1*f(k,l,i+1,j)) &
		+wl2mx*(c5*f(k,l,i,j)+c4*f(k,l,i+1,j)+c2*f(k,l,i+2,j))) 
             flpx=vxp * (wl0px*(c2*f(k,l,i-3,j)+c4*f(k,l,i-2,j)+c5*f(k,l,i-1,j)) &
		+wl1px*(c1*f(k,l,i-2,j)+c3*f(k,l,i-1,j)+c2*f(k,l,i,j)) & 
		+wl2px*(c2*f(k,l,i-1,j)+c3*f(k,l,i,j)+c1*f(k,l,i+1,j)))
             frmx=vxm * (wr0mx*(c1*f(k,l,i-1,j)+c3*f(k,l,i,j)+c2*f(k,l,i+1,j)) &
		+wr1mx*(c2*f(k,l,i,j)+c3*f(k,l,i+1,j)+c1*f(k,l,i+2,j)) &
		+wr2mx*(c5*f(k,l,i+1,j)+c4*f(k,l,i+2,j)+c2*f(k,l,i+3,j)))
             frpx=vxp * (wr0px*(c2*f(k,l,i-2,j)+c4*f(k,l,i-1,j)+c5*f(k,l,i,j)) &
		+wr1px*(c1*f(k,l,i-1,j)+c3*f(k,l,i,j)+c2*f(k,l,i+1,j)) &
		+wr2px*(c2*f(k,l,i,j)+c3*f(k,l,i+1,j)+c1*f(k,l,i+2,j)))
             !negative & positive fluxes y
             flmy=vym * (wl0my*(c1*f(k,l,i,j-2)+c3*f(k,l,i,j-1)+c2*f(k,l,i,j)) &
		+wl1my*(c2*f(k,l,i,j-1)+c3*f(k,l,i,j)+c1*f(k,l,i,j+1)) &
		+wl2my*(c5*f(k,l,i,j)+c4*f(k,l,i,j+1)+c2*f(k,l,i,j+2))) 
             flpy=vyp * (wl0py*(c2*f(k,l,i,j-3)+c4*f(k,l,i,j-2)+c5*f(k,l,i,j-1)) &
		+wl1py*(c1*f(k,l,i,j-2)+c3*f(k,l,i,j-1)+c2*f(k,l,i,j)) &
		+wl2py*(c2*f(k,l,i,j-1)+c3*f(k,l,i,j)+c1*f(k,l,i,j+1)))
             frmy=vym * (wr0my*(c1*f(k,l,i,j-1)+c3*f(k,l,i,j)+c2*f(k,l,i,j+1)) &
		+wr1my*(c2*f(k,l,i,j)+c3*f(k,l,i,j+1)+c1*f(k,l,i,j+2)) &
		+wr2my*(c5*f(k,l,i,j+1)+c4*f(k,l,i,j+2)+c2*f(k,l,i,j+3)))
             frpy=vyp * (wr0py*(c2*f(k,l,i,j-2)+c4*f(k,l,i,j-1)+c5*f(k,l,i,j)) &
		+wr1py*(c1*f(k,l,i,j-1)+c3*f(k,l,i,j)+c2*f(k,l,i,j+1)) &
		+wr2py*(c2*f(k,l,i,j)+c3*f(k,l,i,j+1)+c1*f(k,l,i,j+2)))
             !fluxes x
             flx = flmx + flpx
             frx = frmx + frpx
             !fluxes y
             fly = flmy + flpy
             fry = frmy + frpy
             
             fn(i,j)  =  f(k,l,i,j) - dt/dx * (frx - flx) - dt/dy * (fry - fly)
             f(k,l,i,j) = fn(i,j)
        
        f(k,l,i,1) = fn(i,6)
        f(k,l,i,2) = fn(i,5)
        f(k,l,i,3) = fn(i,4)

        f(k,l,i,ny)  = fn(i,ny-6)
        f(k,l,i,ny-1)= fn(i,ny-5) 
        f(k,l,i,ny-2)= fn(i,ny-4)
        
        f(k,l,1,j) = fn(6,j)
        f(k,l,2,j) = fn(5,j)
        f(k,l,3,j) = fn(4,j)
        
        f(k,l,nx,j)  = fn(nx-6,j)
        f(k,l,nx-1,j)= fn(nx-5,j)
        f(k,l,nx-2,j)= fn(nx-4,j) 
    END DO
    END DO  
END DO
END DO

DO I = 1, NX
DO J = 1, NY
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

if (theta .eq. 0.) go to 1100
      
DO I = 1, NX
DO J = 1, NY
  ZA = 0.1
  ZB = 0.9
  DO WHILE (ABS(ZA-ZB) .GT. 0.0001)
    GA1 = 0
    GB1 = 0
    GA2 = 0
    GB2 = 0
        DO L = 1, 25
            if (theta .eq. 1.) then      
            GA1 = GA1 + (ZA**L) * (-1)**(L-1)/L
            GB1 = GB1 + (ZB**L) * (-1)**(L-1)/L
            GA2 = GA2 + (ZA**L) * (-1)**(L-1)/(L**2)
            GB2 = GB2 + (ZB**L) * (-1)**(L-1)/(L**2)
            else
            GA1 = GA1 + (ZA**L) /L
            GB1 = GB1 + (ZB**L) /L
            GA2 = GA2 + (ZA**L) /(L**2)
            GB2 = GB2 + (ZB**L) /(L**2)
            end if    
        END DO
    PSIA = 2*ET(I,J) - (GA2*(R(I,J)/GA1)**2)/PI - R(I,J)*(UX(I,J)*UX(I,J)+UY(I,J)*UY(I,J))
    PSIB = 2*ET(I,J) - (GB2*(R(I,J)/GB1)**2)/PI - R(I,J)*(UX(I,J)*UX(I,J)+UY(I,J)*UY(I,J))
    ZC = (ZA+ZB)/2
    GC1 = 0
    GC2 = 0
        DO L = 1, 25
            if (theta .eq. 1.) then
            GC1 = GC1 + (ZC**L) * (-1)**(L-1)/L
            GC2 = GC2 + (ZC**L) * (-1)**(L-1)/(L**2)
            else 
            GC1 = GC1 + (ZC**L)/L
            GC2 = GC2 + (ZC**L)/(L**2)
            end if
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
go to 1300

1100 continue
DO I = 1, NX
DO J = 1, NY
  T(I,J) = (2* ET(I,J)/R(I,J))-(UX(I,J)**2+Uy(I,J)**2)
  Z(I,J) = R(I,J)/(pi*T(i,j))
  P(I,j) = r(i,j)*T(i,j)/2
END DO
END DO    

1300 continue
DO I = 1, NX
DO J = 1, NY
   if (iter.eq.floor(sv_point1/DT)) then
     WRITE (20,*) X(I), Y(J), r(I,J), p(I,J), z(I,J)
   else if (iter.eq.floor(sv_point2/dt)) then
     WRITE (30,*) X(I), Y(J), r(I,J), p(I,J), z(I,J)
   end if
END DO 
END DO
     
WRITE (*,*) ITER,time, R(NX/2,NY/2), F(NV/2,NV/2,NX/2,NY/2)
IF (ISTOP .EQ. 1) GOTO 2000
ITER = ITER + 1
GOTO 1000
PAUSE

2000 CONTINUE
DO I = 1, NX
DO J = 1, NY
   WRITE (60,*) X(I), Y(J), r(I,J), p(I,J), z(I,J)
END DO 
END DO


STOP
END PROGRAM    
