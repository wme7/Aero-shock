! A Discontinuous Galerkin implementation for
! Semiclassical Botlzmann ES-BGK for one-dimensional space.
!
! coded and modified by,
! Manuel Diaz, NTU, 2013.07.13
! f99543083'at'.ntu.edu.tw
!
program main
    ! Load modules
    use mathmodule          ! linear Algebra functions
    use dispmodule          ! matlab display functions
    use tecplotmodule       ! write to tecplot functions
    use quadraturemodule    ! Gauss Hermite and Newton cotes abscissas and weights

    ! No implicit variables allowed
    implicit none

    ! Define variables
    integer, parameter  :: nx=3, ny=3       ! XY-grid size
    integer, parameter  :: ngh=10, nnc=200  ! GaussHermite points, Newton-Cotes points
    real, dimension(ngh):: wi,xi
    real, dimension(nx,ny) :: f,feq,f_next  ! Matrix Arrays
    real, dimension(nx) :: x,y   !Vector arrays
    integer :: idnumber,np
    character(len = 100) :: output_file = 'myresults.plt'
    character(len = 100) :: output_file2 = 'myresults.plt'
    real T1, T2

    ! Print message
    print *, 'This is the beginning of DGWENO_ESBGK program'

    call gausshermite(ngh,xi,wi)
    call disp('xi = ',xi); call disp('wi = ',wi);

    ! Calculate CPU time
    call cpu_time(T1)
    print *, 'time for calcualations', T1, 'seconds.'

    ! write to tecplot file
    !call tecplot_write_open(output_file,idnumber) ! open output file and identified No.
        print *, ' '
        print *, 'Opening output file with id No.: ',idnumber
    !call tecplot_write_header(idnumber,'Line data','"X","Y"') ! write header.
        print *, ' '
        print *, 'Writing the Tecplot header'
    !call tecplot_write_xy_line(idnumber,np,x,y) ! write data to file
        print *, ' '
        print *, 'Write data "X","Y" to file'
    !call tecplot_write_close(idnumber) ! close file
        print *, ' '
        print *, 'Number of data point writen was: ',np

    ! Calculate CPU Time
    call cpu_time(T2)
    print *, 'time to write file', T2-T1, 'seconds.'

end program main
