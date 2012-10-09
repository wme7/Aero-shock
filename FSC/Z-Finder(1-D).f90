program z
!this program finds z from et,r,ux,uy
IMPLICIT REAL (A-H,O-Z)
OPEN (UNIT = 10, FILE = 'z.TEC', STATUS = 'UNKNOWN')
PI = ATAN2(1.,1.)*4.
!INPUT
et = 2.4
r  = 1.  
u = 2.

it = 1
!iteration
ZA = 0.001
ZB = 0.999
 do while (abs(za-zb) .gt. 0.000001)
    ga12 = 0
    gb12 = 0
    ga32 = 0
    gb32 = 0
        do l = 1,100
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
   psia = 2*et - ga32*(r/ga12)**3/(2*pi) - r*u**2
   psib = 2*et - gb32*(r/gb12)**3/(2*pi) - r*u**2
        zc = (za + zb)/2
        gc12 = 0
        gc32 = 0
        gc52 = 0
        do l = 1, 100
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
   psic = 2*et - gc32*(r/gc12)**3/(2*pi) - r*u**2
   
        if ((psia*psic) .lt. 0) then
        zb = zc
        else
        za = zc
        end if
   end do
Z1 = ZC 
write (10,*) z1, gc12, gc32
end  