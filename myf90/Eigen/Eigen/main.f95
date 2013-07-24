!Computing EigenValues and EigenVectors
!coded by Manuel Diaz, NTU, 2013.07.13
!
program main
  USE dispmodule
  USE eigenmodule

implicit none

    real :: T1,T2
    real, dimension(6,6) :: a,z
    real, dimension(6) :: d,e,f

    a = reshape( &
        (/0.0,0.5774,0.0,0.0,0.0,0.0, &
          0.5774,0.0,0.5164,0.0,0.0,0.0, &
          0.0,0.5164,0.0,0.5071,0.0,0.0, &
          0.0,0.0,0.5071,0.0,0.5040,0.0, &
          0.0,0.0,0.0,0.5040,0.0,0.5025, &
          0.0,0.0,0.0,0.0,0.5025,0.0/), &
          (/6,6/) )


    call disp('a = ',a)

    call cpu_time(T1) ! time 1

    call tred2(a,d,e) ! convert the matrix into tridiagonal matrix
    call tqli(d,e,z)  ! compute eigen data. Based on Numerical Recepies.

    call cpu_time(T2) ! time 2
    print *, 'Time to compute EigenData: ', T1-T2, 'seconds.'

    call disp('z = ',z) ! EigenVectors Matrix

    ! Normalize the Quadrature Weigths, such that, sum(w) = 2.
    f = (2/sum(z(1,:)*z(1,:)))*(z(1,:)*z(1,:))

    ! Sort with respect to array d
    call sort2(d,f)

    ! display in matlab form the result
    call disp('d = ',d)
    call disp('f = ',f)
    ! This gives us exactly the same result as in my Matlab
    ! function for gauss legendre: [d,f] = GaussLegendre(6)
    ! It took a night to built this program even using NR!
    ! Using fortran for this is ridiculus!! :(

end program main
