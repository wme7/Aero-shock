module initmodule
implicit none
    !
    ! Globals go here
    !
contains

    subroutine Initial_Condition(id,rho_l,u_l,p_l,rho_r,u_r,p_r)
    implicit none
    ! Choose the initial condition to work by it ID number
    integer, intent(in) :: id
    real, intent(out)   :: rho_l,u_l,p_l,rho_r,u_r,p_r
    ! Table with cases
    select case (id)
        case (1)
           print *, 'Shocktube problem of G.A. Sod, JCP 27:1, 1978 ';
            p_l   = 1.00; p_r   = 0.100;
            u_l   = 0.75; u_r   = 0.000;
            rho_l = 1.00; rho_r = 0.125;
        case (2)
            print *, 'Lax test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997 ';
            p_l   = 3.528; p_r   = 0.571;
            u_l   = 0.698; u_r   = 0.000;
            rho_l = 0.445; rho_r = 0.500;
        case (3)
            print *, 'Mach = 3 test case: M. Arora and P.L. Roe: JCP 132:3-11, 1997 ';
            p_l   = 10.333; p_r   = 1.00;
            u_l   =  0.920; u_r   = 3.55;
            rho_l =  3.857; rho_r = 1.00;
        case (4)
            print *, 'Shocktube problem with supersonic zone ';
            p_l   = 1.0;  p_r   = 0.02;
            u_l   = 0.0;  u_r   = 0.00;
            rho_l = 1.0;  rho_r = 0.02;
        case (5)
            print *, 'Contact discontinuity ';
            p_l   = 0.5; p_r   = 0.5;
            u_l   = 0.0; u_r   = 0.0;
            rho_l = 1.0; rho_r = 0.6;
        case (6)
            print *, 'Stationary shock ';
            p_l   =  1.0; p_r   = 0.1;
            u_l   = -2.0; u_l   = -2.0;
            rho_l =  1.0; rho_r = 0.125;
        case default
            print *, 'IC not available';
    end select
    print *, ' ';
    end subroutine Initial_Condition

    !subroutine Build_Elements()

end module initmodule
