      SUBROUTINE coefb(n) 
      use mod_wave 
      INTEGER n 
      INTEGER jj,m 
      allocate(dl(n,n))    
      allocate(cbin(n,n))               
	do m=1,n 
         dl(1,m)=0.d0 
         dl(2,m)=0.d0 
         cbin(1,m)=0.d0 
         cbin(2,m)=0.d0 
      end do 
 
      if (nbc .eq. 1) then 
 
      dl(1,1)=1.d0 
      dl(2,1)=1.d0 
      dl(2,2)=1.d0 
	cbin(1,1)=1.d0 
      cbin(2,1)=1.d0 
      cbin(2,2)=1.d0 
 
      do jj=2,n-1 
      dl(jj+1,1)=1.d0 
	cbin(jj+1,1)=1.d0 
      do m=jj+1,n-1 
      dl(jj+1,m+1)=0.d0 
	cbin(jj+1,m+1)=0.d0 
      end do 
      do m=1,jj 
      dl(jj+1,m+1)=dl(jj-1,m+1)+(2*jj-1)*dl(jj,m) 
	cbin(jj+1,m+1)=cbin(jj,m+1)+cbin(jj,m) 
      end do 
          end do 
      else 
 
 
        dl(1,1)=1.d0 
        do j=2,n 
        dl(j,1)=dl(j-1,1)+j 
        dl(j,j)=dl(j-1,j-1)*j 
		end do 
		do j=2,n 
        do i=j+1,n 
        dl(i,j)=dl(i-1,j)+dl(i-1,j-1)*i 
        end do 
        end do 
 
    	cbin(1,1)=1.d0 
        cbin(2,1)=1.d0 
        cbin(2,2)=1.d0 
 
        do jj=2,n-1 
	     cbin(jj+1,1)=1.d0 
          do m=jj+1,n-1 
           cbin(jj+1,m+1)=0.d0 
          end do 
          do m=1,jj 
           cbin(jj+1,m+1)=cbin(jj,m+1)+cbin(jj,m) 
          end do 
        end do 
      end if 
!	  write(6,*) "nbc=",nbc 
!	  write(6,*) dl 
!      write(6,*) (cbin(4,jj),jj=1,n) 
      return 
      end subroutine coefb
