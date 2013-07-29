module lgl_interp

  implicit none

contains

  subroutine legendre_poly(L0,L0_1,L0_2,p,x)
    integer, intent(in) :: p
    real(kind=4), intent(in) :: x
    real(kind=4), intent(out) :: L0,L0_1,L0_2
    integer :: i
    real(kind=4) :: a, b
    real(kind=4) :: L1, L1_1, L1_2,L2,L2_1,L2_2

    L1 = 0.; L1_1 = 0.; L1_2 = 0.;
    L0 = 1.; L0_1 = 0.; L0_2 = 0.;

    do i= 1,p
       L2 = L1; L2_1 = L1_1; L2_2 = L1_2;
       L1 = L0; L1_1 = L0_1; L1_2 = L0_2;
       a = (2.0*i-1.0)/i
       b = (i-1.)/i
       L0 = a*x*L1 - b*L2
       L0_1 = a*(L1 + x*L1_1) - b*L2_1
       L0_2 = a*(2*L1_1 + x*L1_2) - b*L2_2

    end do

  end subroutine legendre_poly

  subroutine legendre_gauss_lobatto(xgl, wgl, ngl)
    integer, intent(in) :: ngl
    real(kind=4),dimension(ngl), intent(out) :: xgl,wgl
    integer :: p, ph
    integer i, k
    real(kind=4) :: x, dx
    real(kind=4) :: pi
    real(kind=4) :: L0,L0_1,L0_2

    pi = 4.*atan(1.0)
    p = ngl - 1
    ph = floor((p+1.)/2.)
    xgl = 0.0; wgl = 0.0;

    do i = 1,ph
       x = cos((2.0*i-1)*pi/(2.0*p+1))
       do k = 1,20
          call legendre_poly(L0,L0_1,L0_2,p,x)
! get new Newton Iteration
          dx = -(1-x**2)*L0_1/(-2*x*L0_1 + (1-x**2)*L0_2)
          x = x+dx
          if (abs(dx) .lt. 1.0e-20) then
             exit ! end the loop
          endif
       end do
       xgl(p+2-i)=x
       wgl(p+2-i)=2./(p*(p+1)*(L0**2))
    end do

    
    ! check for zero root
    if ((p+1) .ne. 2*ph) then
       x=0.0
       call legendre_poly(L0,L0_1,L0_2,p,x)
       xgl(ph+1)=x
       wgl(ph+1)=2./(p*(p+1)*(L0**2))

    endif

    ! find remainder of roots via symmetry
    
    do i = 1,ph
       xgl(i) = -xgl(p+2-i)
       wgl(i) = wgl(p+2-i)
    end do



  end subroutine legendre_gauss_lobatto

  subroutine lagrange_basis(psi, dpsi, xnq, wnq,ngl,nq,xgl)
    integer, intent(in) :: ngl, nq
    real(kind=4), dimension(ngl),intent(in) :: xgl
    real(kind=4), dimension(nq), intent(inout) :: xnq, wnq
    real(kind=4), dimension(ngl,nq), intent(inout) :: psi, dpsi
    integer :: l, i, j, k
    real(kind=4) :: xl, xi, xj, ddpsi, xk
    
    call legendre_gauss_lobatto(xnq,wnq,nq)
    psi = 0.0; dpsi = 0.0;
    !Perform Quadrature
    do l = 1,nq
       xl = xnq(l)
       do i = 1,ngl
          xi = xgl(i)
          psi(i,l)=1.0
          dpsi(i,l) = 0.0
          do j = 1,ngl
             xj = xgl(j)
             if (j .ne. i) then
                psi(i,l) = psi(i,l)*(xl-xj)/(xi-xj)
             endif !(j .ne. i)
             ddpsi = 1.
             if (j .ne. i) then
                do k = 1,ngl
                   xk = xgl(k)
                   if((k .ne. i) .and. (k .ne. j)) then
                      ddpsi = ddpsi*(xl-xk)/(xi-xk)
                   endif !((k .ne. i) ...
                end do ! k
                dpsi(i,l) = dpsi(i,l)+ddpsi/(xi-xj)
             endif ! (j .ne. i)...
          end do ! j
       end do ! i
    end do ! l

  end subroutine lagrange_basis

end module lgl_interp
