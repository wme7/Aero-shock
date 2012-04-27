program testset

  ! include gloval data
  use dg
  !make available the setval rouitne
  use setval

  real :: x, y

  ! read values from terminal
  read(*,*) a_g, b_g

  ! set x and y and check output
  call setv(x, y)
  write(*,*) x, y

end program testset
