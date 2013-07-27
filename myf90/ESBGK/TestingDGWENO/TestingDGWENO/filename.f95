! For Testing the DG&WENO subroutines
! By Manuel Diaz

program main
    implicit none
    ! Define
	integer :: mo,mt,kflux,phase_quadrature_rule,init_value,ierror,n
    integer :: k = 50
	real    :: cflc,tau,tprint,xmmm
    ! Namelist
	namelist /proj_list/ mt,kflux,cflc,tau,phase_quadrature_rule,init_value,ierror,tprint,n,xmmm,mo
    ! character names
    character(len=12) :: filename
    character(len=12) :: format_string
    character(len=12) :: name
    integer :: counter

    do counter=1, 10
        if (counter < 10) then
            format_string = "(A5,I1)"
        else
            format_string = "(A5,I2)"
        endif

        write (filename,format_string) "hello", counter
        print *, trim(filename)
    end do

    write(name,"(A5,I2,A5)") "hello",10,"folks" ! <-- this is the key!!!
    print *,'this is a filename: ', name

    ! Prepare to read file
    open(10, file="proj.in", form="FORMATTED", action="READ")
    read(10, nml=proj_list)

    ! Write namelist output
    write(*, nml=proj_list)

   print *, 'this is mt: ',mt

end
