      subroutine calc
c     ***************
      include 'com'
c
c      --- the upwinding-scheme   
c      --- ******************
      if (nsch.eq.1) then
         call rk
      else
c
      write(06,10) nsch
   10 format(t5,' the scheme number ', i3, ' has not been implemented!')  
      stop
c
      endif 
c          
      end
