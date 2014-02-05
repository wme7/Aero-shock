      subroutine calc
!     ***************
      use mod_wave
!
!      --- the upwinding-scheme   
!      --- ******************
      if (nsch.eq.1) then
         call rk
      else
!
      write(06,10) nsch
   10 format(t5,' the scheme number ', i3, ' has not been implemented!')  
      stop
!
      endif 
!          
      end subroutine calc
