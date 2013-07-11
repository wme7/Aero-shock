! Example program of dispmodule demonstrating use from Matlab mex files
!
! To try this program, set up Fortran mexing (on Windows one may use
! Gnumex; see gnumex.sourceforge.net) and make a version of dispmodule
! (say in dispmodule_mex.f90) by uncommenting the putstrmodule that
! calls mexPrintf, and commenting the dummy putstrmodule. Then issue
! from the Matlab prompt:
!
!     >> mex -c dispmodule_mex.f90
!     >> mex mexdispdemo.f90 dispmodule.obj
!     >> mexdispdemo(hilb(3), 1:4)
!
! If all is well one should obtain the output:
!
! A = 1.00000  0.50000  0.33333   b = 1.00000
!     0.50000  0.33333  0.25000       2.00000
!     0.33333  0.25000  0.20000       3.00000
!                                     4.00000
!
! This version works with the g95 compiler. With the gfortran compiler 
! a different version of the interface-block must be used, namely:
!
!   interface
!     function mxgetpr(pm) bind(c, name = 'MXGETPR')
!       import c_int, c_ptr
!       integer(c_int) :: pm
!       type(c_ptr) mxgetpr
!     end function mxgetpr
!     !
!     function mxgetm(pm) bind(c, name = 'MXGETM')
!       import c_int
!       integer(c_int) :: pm, mxgetm
!     end function mxgetm
!     !
!     function mxgetn(pm) bind(c, name = 'MXGETN')
!       import c_int
!       integer(c_int) pm, mxgetn
!     end function mxgetn
!     !
!   end interface

subroutine mexfunction(nlhs, plhs, nrhs, prhs)
  implicit none
  interface
    function mxgetpr(pm)
      integer :: pm
      integer, pointer :: mxgetpr
    end function mxgetpr
    !
    function mxgetm(pm)
      integer :: pm, mxgetm
    end function mxgetm
    !
    function mxgetn(pm)
      integer pm, mxgetn
    end function mxgetn
    !
  end interface
  integer :: nlhs, nrhs, plhs(nlhs), prhs(nrhs), m, n, nb
  integer, pointer :: Ap, bp
  m = mxgetm(prhs(1))
  n = mxgetn(prhs(1))
  nb = max(mxgetm(prhs(2)),mxgetn(prhs(2)))
  Ap => mxgetpr(prhs(1))
  bp => mxgetpr(prhs(2))
  call mexdispdemo(Ap, m, n, bp, nb)
end subroutine mexfunction

subroutine mexdispdemo(A, m, n, b, nb)
  use dispmodulemex
  implicit none
  integer m, n, nb
  double precision A(m, n), b(nb)
  call disp('A = ', A, advance = 'no')
  call disp('b = ', b)
end subroutine mexdispdemo
