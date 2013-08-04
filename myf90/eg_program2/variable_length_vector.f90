! Variable-length array output :D 
! By NickFort 23 July 2012, modified by Manuel Diaz
! NTU, 2013.07.31. Original code available at www.tek-tips.com
! use gfortran for compiling 
!
program variable_length_vector
implicit none

interface
    subroutine random_vector(vector)
    implicit none
    real, allocatable, intent(out) :: vector(:)
    end subroutine
end interface

integer :: length_call,k
real, allocatable :: vector_call(:)

call random_vector(vector_call)

length_call = size(vector_call)

print *, "length_call: ", length_call
do k=1,size(vector_call)
    print *, vector_call(k)
end do

deallocate(vector_call)

end program variable_length_vector

!==============================================

subroutine random_vector(vector)
implicit none
real, allocatable, intent(out) :: vector(:)
integer :: length
real :: length_real

if (allocated(vector)) deallocate(vector)
call init_random_seed()
call random_number(length_real)
length = int(length_real*10)+1
allocate (vector(length))
call random_number(vector)

end subroutine random_vector

!==============================================

subroutine init_random_seed()
integer :: i, n, clock
integer, allocatable :: seed(:)
call random_seed(size = n)
allocate(seed(n))
call system_clock(count=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
call random_seed(put = seed)
deallocate(seed)

end subroutine init_random_seed

