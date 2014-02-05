program hello
    ! Load modules
    use tecplotmodule   ! write to tecplot functions

    ! Define Variables
    real, dimension(10)  :: x,f !Vector arrays
    integer :: idnumber,np
    character(len = 100) :: output_file = 'myresults.plt'

    ! Calculations and displays
    x = (/0,1,2,3,4,5,6,7,8,9/)
    f = sin(x)
    np = size(f)

    ! write to tecplot file
    call tecplot_write_open(output_file,idnumber) ! open output file and identified No.
        print *, ' '
        print *, 'Opening output file with id No.: ',idnumber
    call tecplot_write_header(idnumber,'Line data','"X","F"') ! write header.
        print *, ' '
        print *, 'Writing the Tecplot header'
    call tecplot_write_xy_line(idnumber,np,x,f) ! write data to file
        print *, ' '
        print *, 'Write data "X","F" to file'
    call tecplot_write_close(idnumber) ! close file
        print *,
        print *, 'Number of data point writen was: ',np

end program
