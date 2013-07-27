! For Testing the DG&WENO Subroutines
! By Manuel Diaz, NTU, 2013.06.16

program main
        ! Load modules
        use IDmodule

        ! No implicit definitions allowed
        implicit none

        ! Name of simulation
        character(len=7) :: name = "SBBGK1d"

        ! Define manually inputs
        !integer :: theta = 0
        !integer :: nx = 100
        !integer :: P_deg = 3
        !integer :: RK_stages = 3
        !integer :: IC_case  = 1
        !integer :: fmodel = 1
        !integer :: f_case = 1
        !integer :: method = 1
        !real    :: tau = 1.0/10000

        ! Expected Inputs
        integer :: theta,nx,P_deg,RK_stages,IC_case,fmodel,f_case,method,kflux
        real    :: tau,tEnd,MM

        ! IDs for files
        character(len=100) :: IDn, IDf

        ! Name list (comment if parameters where to be setup manually)
        namelist /parameters_list/ name,theta,nx,P_deg,RK_stages,IC_case,fmodel,f_case,method,kflux,tau,tEnd,MM
        ! this list need not to be in order

        ! Read file with parameters
        open(10, file="configuration.in", form="FORMATTED", action="READ")
        read(10, nml=parameters_list)

        ! Checking parametes read
        write(*, nml=parameters_list)

        ! Create IDs for simulation and result file
        call ID_name(name,theta,nx,P_deg,RK_stages,tau,IC_case,fmodel,f_case,method,IDn,IDf)

        ! Checking IDs
        print *, 'IDn: ',IDn
        print *, 'IDf: ',IDf

end program
