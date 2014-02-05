program factorial
  implicit none

  ! variables and intial values
  integer :: nfact = 1
  integer :: n

  ! compute factorials
  do n = 1, 100
     nfact=nfact*n
     write(*,*) n, nfact
     if (n>10) exit
  end do

end program factorial
