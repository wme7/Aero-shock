! A Discontinuous Galerkin implementation for
! Semiclassical Botlzmann ES-BGK for one-dimensional space.
!
! coded and modified by,
! Manuel Diaz, NTU, 2013.07.13
! f99543083'at'.ntu.edu.tw
!
program hello
    ! Load modules
    use mathmodule      ! linear Algebra functions
    use dispmodule      ! matlab display functions
    use tecplotmodule   ! write to tecplot functions

    ! Define Variables
    implicit none
    integer, parameter :: n=3, m=3
    real, dimension(n,m) :: A,B !Matrix Arrays

    ! Print message
    print *, 'This is the beginning of DG-ESBGK'

    ! Calculations and displays
    A = transpose( reshape( (/1,2,3,4,5,6,7,8,9/),shape(A)))
    B = transpose( reshape( (/2,0,0,0,2,0,0,0,2/),shape(B)))
    call disp('A = ',A)
    call disp('B = ',B)

    ! write to tecplot file


end program

