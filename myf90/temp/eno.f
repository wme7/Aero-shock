    subroutine enoReconst( tv, tvL,tvR ) 
    use mainvar 
    implicit none 
    real :: tv(6), tvL, tvR 
    real :: a,b,aa,bb 
    integer :: i 
 
    !! (i+1/2)+ 
    a=abs(tv(5)-tv(4)) 
    b=abs(tv(4)-tv(3)) 
    if( a>b )then 
      aa=abs(tv(5)-2*tv(4)+tv(3)) 
      bb=abs(tv(4)-2*tv(3)+tv(2)) 
      if(aa > bb)then 
        tvR= -tv(2)/6.0 + 5.0*tv(3) / 6.0 + tv(4) / 3.0 
      else 
        tvR= tv( 3)/3.0 + 5.0*tv(4) / 6.0 - tv(5) / 6.0 
      endif 
    else 
      aa=abs(tv(6)-2.*tv(5)+tv(4)) 
      bb=abs(tv(5)-2.*tv(4)+tv(3)) 
      if(aa>bb) then 
        tvR= tv( 3)/3.0 + 5.0*tv(4) / 6.0 - tv(5) / 6.0 
      else 
        tvR= 11.0*tv( 4)/6.0 - 7.0*tv(5) / 6.0 + tv(6) / 3.0 
      endif 
    endif 
 
    !! (i+1/2)- 
    a=abs(tv(4)-tv(3)) 
    b=abs(tv(3)-tv(2)) 
    if( a>b )then 
      aa=abs(tv(4)-2*tv(3)+tv(2)) 
      bb=abs(tv(3)-2*tv(2)+tv(1)) 
      if(aa > bb)then 
        tvL= tv(1)/3.0 - 7.0*tv(2) / 6.0 +11.0*tv(3) / 6.0 
      else 
        tvL=-tv(2)/6.0 + 5.0*tv( 3) / 6.0 +     tv(4) / 3.0 
      endif 
    else 
      aa=abs(tv(5)-2.*tv(4)+tv(3)) 
      bb=abs(tv(4)-2.*tv(3)+tv(2)) 
      if(aa>bb) then 
        tvL=-tv(2)/6.0 + 5.0*tv( 3) / 6.0 +     tv(4) / 3.0 
      else 
        tvL= tv(3 )/3.0 + 5.0*tv( 4) / 6.0 -     tv(5) / 6.0 
      endif 
    endif 
 
    end subroutine enoReconst