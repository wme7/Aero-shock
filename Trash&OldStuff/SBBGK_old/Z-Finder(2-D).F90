program z
!this program finds z from et,r,ux,uy
IMPLICIT REAL (A-H,O-Z)
OPEN (UNIT = 10, FILE = 'z.TEC', STATUS = 'UNKNOWN')
PI = ATAN2(1.,1.)*4.
!INPUT
et = 0.22263
r  = 0.521295
ux = -0.42308
uy = 0.
!iteration
ZA = 0.001
ZB = 0.999
DO WHILE (ABS(ZA-ZB) .GT. 0.00001)
    GA1 = 0
    GB1 = 0
    GA2 = 0
    GB2 = 0
    DO L = 1, 50
    GA1 = GA1 + (ZA**L) * (-1)**(L-1)/L
    GB1 = GB1 + (ZB**L) * (-1)**(L-1)/L
    GA2 = GA2 + (ZA**L) * (-1)**(L-1)/(L**2)
    GB2 = GB2 + (ZB**L) * (-1)**(L-1)/(L**2)
    END DO
    PSIA = 2*ET - (GA2*(R/GA1)**2)/PI - R*(UX*UX+UY*UY)
    PSIB = 2*ET - (GB2*(R/GB1)**2)/PI - R*(UX*UX+UY*UY)
    ZC = (ZA+ZB)/2
    GC1 = 0
    GC2 = 0
    DO L = 1, 50
    GC1 = GC1 + (ZC**L) * (-1)**(L-1)/L
    GC2 = GC2 + (ZC**L) * (-1)**(L-1)/(L**2)
    END DO
    PSIC = 2*ET - (GC2*(R/GC1)**2)/PI - R*(UX*UX+UY*UY)    
    IF ((PSIA*PSIC) .LT. 0) THEN
    ZB = ZC
    ELSE
    ZA = ZC
    END IF
END DO
Z1 = ZC 
write (10,*) z1, gc2, gc1
end  