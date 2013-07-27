module IDmodule
    !
    ! Define Globals here
    !
contains
    subroutine ID_name(name,theta,nx,P_deg,RK_stages,tau,IC_case,fmodel,f_case,method,IDn,IDf)
    ! Define the inputs and outputs
    character(len=*), intent(in)  :: name
    integer, intent(in)           :: theta,nx,P_deg,RK_stages,IC_case,fmodel,f_case,method
    real, intent(in)              :: tau
    character(len=*), intent(out) :: IDn,IDf
    ! Define working variables
    character(len=20) :: name1,name2,name3,statistic,feq,advec,p_degree,elements,RKs,IC,omega,f
    ! Define format for ID files names
    character(len=80) :: format_string,format_string2

    ! Break Name into useful parts "SBBGK1d"
    name1=name(1:2)
    name2=name(3:5)
    name3=name(6:7)

    ! Define statistic
    select case (theta)
        case (-1)
            statistic = "BE"
        case (0)
            statistic = "MB"
        case (1)
            statistic = "FD"
        case default
            print *, 'Not a valid statistic'
    end select

    ! Define equilibrium distirbution model useq
    select case (fmodel)
        case (1) ! UU
            feq = "-UU"
        case (2) ! ES
            feq = "-ES"
        case default
            print *, 'Model not available'
    end select

    ! Define the Method to use
    select case (method)
        case (1)
            advec = "Upwind"
            P_degree = "1"
        case (2)
            advec = "*TVD**"
            P_degree = "1"
        case (3)
            advec = "WENO_3"
            P_degree = "3"
        case (4)
            advec = "WENO_5"
            P_degree = "5"
        case (5)
            advec = "DGWENO"
            P_degree = char(P_deg+48)
        case default
            print *, 'Advection method not available'
    end select

    ! Define the number of cells to be used,
    write(elements,"(A1,I3)") 'X',nx

    ! Define the number of RK stages
    write(Rks,"(I1)") RK_stages

    ! Define the ID of Initial Condition used
    write(IC,"(I1)") IC_case

    ! Define how to evolve f,
        select case (f_case)
            case (1) ! with Collision term-BGK approx
                write(omega,"(A1,I5)") 'W',ceiling(1.0/tau)
            case (2) ! no-collison-Euler limit
                omega = "EL"
            case default
                print *, 'case not available'
        end select

    ! Define Format for file
    f = ".plt"

    format_string = "(A2,A3,A3,A1,A2,A6,A2,A2,I3,A1,I1,A2,I1,A1,A6,A2,A1,A4)"
    write(IDf,format_string) trim(name1),trim(feq),trim(name2),'-', &
                             trim(statistic),trim(advec),trim(name3), &
                             '-X',nx,'P',P_deg,'RK',RK_stages,'-',omega,'IC',IC,trim(f)
    format_string2 = "(A2,A3,A3,A1,A2,A6,A2,A2,I3,A1,I1,A2,I1,A1,A6,A2,A1)"
    write(IDn,format_string2) trim(name1),trim(feq),trim(name2),' ', &
                             trim(statistic),trim(advec),trim(name3), &
                             ' X',nx,'P',P_deg,'RK',RK_stages,' ',omega,'IC',IC
    end subroutine ID_name

    subroutine output_names
        ! I desided to preserve this idea for the future
        implicit none

        ! character names
        character(len=12) :: filename
        character(len=12) :: format_string
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

    end subroutine output_names

end module IDmodule
