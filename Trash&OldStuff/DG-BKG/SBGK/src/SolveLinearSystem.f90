module LinearSolvers
   implicit none

   ! The default value for the smallest pivot that will be accepted
   ! using the LinearSolvers subroutines.  Pivots smaller than this 
   ! threshold will cause premature termination of the linear equation 
   ! solver and return false as the return value of the function.

   real(kind=8), parameter :: DEFAULT_SMALLEST_PIVOT = 1.0d-13

contains

   ! Use Gaussian elimination to calculate the solution to the linear 
   ! system, A x = b.  No partial pivoting is done.  If the threshold 
   ! argument is present, it is used as the smallest allowable pivot 
   ! encountered in the computation; otherwise, DEFAULT_SMALLEST_PIVOT, 
   ! defined in this module, is used as the default threshold.  The status
   ! of the computation is a logical returned by the function indicating
   ! the existence of a unique solution (.true.), or the nonexistence of
   ! a unique solution or threshold passed (.false.).

   ! Note that this is an inappropriate method for some linear systems.
   ! In particular, the linear system, M x = b, where M = 10e-12 I, will 
   ! cause this routine to fail due to the presence of small pivots.  
   ! However, this system is perfectly conditioned, with solution x = b.

   function gaussianElimination( A, b, x, threshold )
      implicit none
      logical gaussianElimination
      real(kind=8), dimension( :, : ), intent( in ) :: A   ! Assume the shape of A.
      real(kind=8), dimension( : ), intent( in ) ::  b     ! Assume the shape of b.
      real(kind=8), dimension( : ), intent( out ) :: x     ! Assume the shape of x.

      ! The optional attribute specifies that the indicated argument
      ! is not required to be present in a call to the function.  The
      ! presence of optional arguments, such as threshold, may be checked
      ! using the intrinsic logical function, present (see below).

      real(kind=8), optional, intent( in ) :: threshold

      integer i, j   ! Local index variables.
      integer N      ! Order of the linear system.
      real(kind=8) m         ! Multiplier.
      real(kind=8) :: smallestPivot = DEFAULT_SMALLEST_PIVOT

      ! Pointers to the appropriate rows of the matrix during the elmination.
      real(kind=8), dimension( : ), pointer :: pivotRow
      real(kind=8), dimension( : ), pointer :: currentRow

      ! Copies of the input arguments.  These copies are modified during
      ! the computation.
      ! The target attribute is used to indicate that the specified 
      ! variable may be the target of a pointer.  Rows of ACopy are targets
      ! of pivotRow and currentRow, defined above.

      real(kind=8), dimension( size( A, 1 ), size( A, 2) ), target :: ACopy
      real(kind=8), dimension( size( b ) ) :: bCopy

      ! Status of the computation.  The return value of the function.
      logical successful

      ! Change the smallestPivot if the threshold argument was included.
      if ( present( threshold ) ) smallestPivot = abs( threshold )

      ! Setup the order of the system by using the intrinsic function size.
      ! size returns the number of elements in the specified dimension of
      ! an array or the total number of elements if the dimension is not
      ! specified.  Also assume that a unique solution exists initially.

      N = size( b )   
      ACopy = A
      bCopy = b
      successful = .true.

      ! Begin the Gaussian elimination algorithm.
      ! Note the use of array sections in the following loops.  These 
      ! eliminate the need for many do loops that are common in Fortran 
      ! 77 code.
      ! Pointers are also used below and enhance the readability of the
      ! elimination process.

      ! Begin with the first row.
      i = 1

      ! Reduce the system to upper triangular.
      do while ( ( successful ) .and. ( i <= N-1 ) )

         ! The following statement is called pointer assignment and uses
         ! the pointer assignment operator =>'.  This causes pivotRow 
         ! to be an alias for the ith row of ACopy.  Note that this does
         ! not cause any movement of data.

         ! Assign the pivot row.
         pivotRow => ACopy( i, : )

         ! Verify that the current pivot is not smaller than smallestPivot.
         successful = abs( pivotRow( i ) ) >= smallestPivot

         if ( successful ) then

	    ! Eliminate the entries in the pivot column below the pivot row.

	    do j = i+1, N
	       ! Assign the current row.
	       currentRow => ACopy( j, : )

               ! Calculate the multiplier.
               m = currentRow( i ) / pivotRow( i ) 

               ! Perform the elimination step on currentRow and right 
               ! hand side, bCopy.
               currentRow = m * pivotRow - currentRow
               bCopy( j ) = m * bCopy( i ) - bCopy( j )
	    end do

         end if

         ! Move to the next row.
         i = i + 1

      end do

      ! Check the last pivot.
      pivotRow => ACopy( N, : )
      if ( successful ) successful = abs( pivotRow( N ) ) >= smallestPivot

      if ( successful ) then
         do i = N, 2, -1   ! Backward substitution.

            ! Determine the ith unknown, x( i ).
            x( i ) = bCopy( i ) / ACopy( i, i )

            ! Substitute the now known value of x( i ), reducing the order of 
            ! the system by 1.
            bCopy = bCopy - x( i ) * ACopy( :, i )

         end do
      end if

      ! Determine the value of x( 1 ) as a special case.
      if ( successful ) x( 1 ) = bCopy( 1 ) / ACopy( 1, 1 )

      ! Prepare the return value of the function.
      gaussianElimination = successful

   end function gaussianElimination


   ! The LU decomposition of a matrix may be represented in a compact form
   ! existing in a single matrix, M,  if the assignments M=L and M=U are 
   ! done (in that order).  The diagonal entries in L are assumed to be 
   ! unity so that no storage space is necessary.  Instead, the diagonal
   ! of M is used to hold the diagonal entries of U.  This is a common 
   ! method of storing the LU decomposition of a matrix.

   ! The algorithm belows makes an additional assumption concerning the 
   ! pivots or diagonal elements of U.  Computation terminates if one of
   ! these pivots is smaller than the given or default threshold.  In this
   ! case, the LU decomposition is not formed.  Note that this algorithm 
   ! successfully terminates if such an LU can be computed.  In this case
   ! the coefficient matrix, A, is nonsingular.  (No attempt for recovery,
   ! such as permutation of rows, is done.)

   ! Compute the LU decomposition of A, storing the result in LU so that
   ! A is not overwritten.  If the threshold argument is present, it is used 
   ! as the smallest allowable pivot encountered in the computation;
   ! otherwise, DEFAULT_SMALLEST_PIVOT, defined in this module, is used as
   ! the default threshold during the computation.  The status of the
   ! computation is a logical returned by the function indicating the
   ! success (.true.) or failure (.false.) of the factorization
   ! After the computation, LU will contain the multipliers below the main
   ! diagonal (L) and the result after elimination on and above the main
   ! diagonal (U), so that A = L * U.

   function LUFactor ( A, LU, threshold ) 
      implicit none
      logical LUFactor
      real(kind=8), dimension( :, : ), intent( in ) :: A
      real(kind=8), dimension( :, : ), intent( out ) :: LU 
      real(kind=8), optional, intent( in ) :: threshold

      integer k, i
      integer N
      logical successful   ! Status of the computation.
      real(kind=8) :: smallestPivot = DEFAULT_SMALLEST_PIVOT

      ! Reassign the smallestPivot, set the order of the system, and 
      ! copy A into LU as it will be written to during the factorization.

      if ( present( threshold ) ) smallestPivot = abs( threshold )
      N = size( A, 1 )   
      LU = A

      ! Begin the LU factorization algorithm.
      ! The status of the computation is initially successful.
      successful = .true.

      k = 1   ! Begin with the first column.
      do while ( ( successful ) .and. ( k <= N-1 ) )

         ! Verify that the kth pivot is not smaller than smallestPivot.
         successful = abs( LU( k, k ) ) >= smallestPivot

         if ( successful ) then
            ! Calculate the multipliers (L) for the current column.
            LU( k+1:N, k ) = LU( k+1:N, k ) / LU( k, k )

            ! Perform elimination on the upper portion of the matrix (U). 
            do i = k+1, N
               LU( i, k+1:N ) = LU( i, k+1:N ) - LU( i, k ) * LU( k, k+1:N )
            enddo

            k = k + 1   ! Move to the next column.
         end if

      enddo

      ! Prepare the return value of the function.
      LUFactor = successful

   end function LUFactor


   ! Let A = L*U where LU represents the LU decomposition of A stored in the
   ! format produced by LUFactor, A, L, U in R**(NxN).
   ! Solve the linear system, A x = b, using the LU decomposition of A stored
   ! in LU.  Since LU is the LU decomposition of A, A is nonsingular.  
   ! Consequently, the columns of A constitute a basis for R**N.   So, there 
   ! must exist a unique solution to the linear system A x = b.
   ! LUSolve returns the solution to this linear system.

   function LUSolve( LU, b ) result( x )
      implicit none
      real(kind=8), dimension( :, : ), intent( in ) :: LU
      real(kind=8), dimension( : ), intent( in ) :: b
      real(kind=8), dimension( size( b ) ) :: x

      integer k
      integer N
      real(kind=8), dimension( size( b ) ) :: bCopy

      ! Determine the order of the system and store a copy of b in bCopy
      ! as it is written during the computation.
      N = size( b )
      bCopy = b

      ! Assume LU is in the form of LU and solve the system in two steps.
      ! First, using forward elmination to solve L y = b, then using
      ! backward elmination to solve U x = y.  In both cases, the right
      ! hand side is overwritten with the solution as it is computed.

      ! Forward elimination.  Store the solution into the right hand side.
      do k = 1, N-1
         bCopy( k+1:N ) = bCopy( k+1:N ) - bCopy( k ) * LU( k+1:N, k )
      end do

      ! Backward elimination.  Store the solution into the right hand side.
      do k = N, 2, -1
         bCopy( k ) = bcopy( k ) / LU( k, k )
         bCopy( 1:k-1 ) = bCopy( 1:k-1 ) - bCopy( k ) * LU( 1:k-1, k )
      end do

      ! Solve for the 1st unknown as a special case.
      bCopy( 1 ) = bCopy( 1 ) / LU( 1, 1 )

      ! Assign a return value for the function via its result variable, x.
      x = bCopy

   end function LUSolve


   ! Output A in Matlab format, using name in the Matlab assignment statement.
   subroutine printMatrix( A, name )
      implicit none
      real(kind=8), dimension( :, : ) :: A   ! Assume the shape of A.
      character name  ! Name for use in assignment, ie, name = ......

      integer n, m, i, j

      n = size( A, 1 )
      m = size( A, 2 )

      write( *, fmt="(a1,a5)", advance = "no" ) name, ' = [ '

      ! Output the matrix, except for the last row, which needs no ;'.
      do i = 1, n-1

         ! Output current row.
         do j = 1, m-1
            write( *, fmt="(f16.8,a2)", advance = "no" ) A( i, j ), ', '
         end do 

         ! Output last element in row and end current row.
         write( *, fmt="(f16.8,a1)" ) A( i, m ), ';'

      end do 

      ! Output the last row.
      do j = 1, m-1
         write( *, fmt="(f16.8,a2)", advance = "no" ) A( i, j ), ', '
      end do 

      ! Output last element in row and end.
      write( *, fmt="(f10.6,a1)" ) A( i, m ), ']'

   end subroutine printMatrix


   ! Output b in Matlab format, using name in the Matlab assignment statement.
   subroutine printVector( b, name )
      implicit none
      real(kind=8), dimension( : ) :: b   ! Assume the shape of b.
      character name   ! Name for use in assignment, ie, name = ......

      integer n, i

      n = size( b )

      write( *, fmt="(a1,a5)", advance = "no" ) name, ' = [ '

      do i = 1, n-1
         write( *, fmt = "(f10.6,a2)", advance = "no" ) b( i ), ', '
      end do

      write( *, fmt = "(f10.6,a2)" ) b( n ), ']'''

   end subroutine printVector

end module LinearSolvers




! A program to solve linear systems using the LinearSolvers module.

program SolveLinearSystem

   ! Include the module for the various linear solvers.

   use LinearSolvers
   implicit none

   integer, parameter :: N =  5 ! Order of the linear system.
   real(kind=8), parameter    :: TOO_SMALL = 1.0d-12  ! Threshold for pivots.

   ! Declare the necessary arrays and vectors to solve the linear system
   ! A x = b.

   real(kind=8), dimension( N, N ) :: A   ! Coefficient matrix.
   real(kind=8), dimension( N ) :: x, b   ! Vector of unknowns, and right hand side.
   real(kind=8), dimension( N, N ) :: LU  ! Matrix for LU factorization of A.

   logical successful   ! Status of computations.


   ! The intrinsic subroutine, random_number, fills a real(kind=8) array or scalar,
   ! with uniformly distributed random variates in the interval [0,1).

   call random_number( A )   ! Initialize the coefficient matrix.
   call random_number( b )   ! Initialize the right-hand side.

   ! Output the matrix in Matlab format for ease of checking the solution.
   call printMatrix( A, 'A' )
   call printVector( b, 'b')

   ! Use Gaussian elmination to calcuate the solution of the linear system.
   ! The call below uses the default threshold specified in the 
   ! LinearSolvers module by omitting the optional argument.

   successful = gaussianElimination( A, b, x )

   print *, '===================================='
   print *, 'Gaussian Elimination:'
   print *, '------------------------------------'
   if ( successful ) then
      call printVector( x, 'x' )
      print *, 'Infinity Norm of Difference = ',  &
	 maxval( abs ( matmul( A, x ) - b ) )
   else
      print *, 'No unique solution or threshold passed.'
   end if

   ! Compute the LU decomposition of A.

   successful = LUFactor( A, LU ) 

   ! Calculate the solution of the linear system given the LU decomposition.
   ! Output the results.

   print *
   print *, '===================================='
   print *, 'LU Factorization:'
   print *, '------------------------------------'
   if ( successful ) then
      x = LUSolve(LU, b)

      print *, 'LU Decomposition:'
      call printMatrix( LU, 'M' )
      call printVector( x, 'x' )
      print *, 'Infinity Norm of Difference = ', &
	 maxval( abs ( matmul( A, x ) - b ) )
   else
      print *, 'No unique solution or threshold passed.'
   end if

end program SolveLinearSystem
