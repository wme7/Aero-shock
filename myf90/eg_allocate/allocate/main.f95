! A fortran95 example on using allocation of uknown size arrays
! coded by Manuel Diaz, NTU, 2013.07.21
!
program main
    ! load modules
    use dispmodule

    ! define no implicit variables
    implicit none

    ! define variables and arrays
    integer, allocatable :: N(:,:), E(:,:)
    integer, allocatable :: u(:)
    integer :: i,j,c
    integer :: g(18)

    ! now lets compute the final size for N
    print *, 'input matrix size: '
    read (*,*) c

    ! allocate matrix of unknown sise
    allocate(N(c,2*c)) ! first matrix
    allocate(E(c,2*c)) ! second matrix
    allocate(u(size(g)))   ! for initial condition

    ! make sure the matrixes are clean, we don't want surprices from the compiler!
    N = 0
    E = 0

    g(10:18) = 4
    g(1:9) = 1
    where (g > 2)
        u = 1
    elsewhere
        u = 0
    end where

    ! print result
    call disp('u = ',u)

    forall (i = 1:c)
        N(i,i) = 3
        E(i,i+c) = 2
    end forall

    ! print result
    call disp("N = ",N)
    print *, ' '
    call disp("E = ",E)

    N(:,1:2) = 1
    N(:,3:4) = 2
    N(:,5:6) = 0
    E = 0
    where(N <= 1)
        E = 3
    end where

    ! print result
    call disp("N = ",N)
    print *, ' '
    call disp("E = ",E)

end
