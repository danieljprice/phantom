module einsteintk_wrapper
!
!
! This module is a "wrapper" for the hydro evol + communication with ET
! Subroutines here should be called by ET rather than calling phantom subroutines
! directly 
!
    implicit none 
    contains

    subroutine init_et2phantom(infilestart,dt_et)
        ! Wrapper that intialises phantom
        ! Intended to hide all of the inner works of phantom from ET
        ! Majority of the code from HelloHydro_init has been moved here 

        use io,              only:id,master,nprocs,set_io_unit_numbers,die
        use mpiutils,        only:init_mpi,finalise_mpi
        use initial,         only:initialise,finalise,startrun,endrun
        use evolve,          only:evol_init
        use tmunu2grid
        use einsteintk_utils
        

        implicit none
        character(len=*),  intent(in) :: infilestart
        real,          intent(in) :: dt_et 
        !character(len=500) :: logfile,evfile,dumpfile,path
        integer :: i,j,k,pathstringlength
    
        ! For now we just hardcode the infile, to see if startrun actually works!
        ! I'm not sure what the best way to actually do this is? 
        ! Do we store the phantom.in file in par and have it read from there?
        !infile = "/Users/spencer/phantomET/phantom/test/flrw.in"
        !infile = trim(infile)//'.in'
        !print*, "phantom_path: ", phantom_path   
        !infile = phantom_path // "flrw.in" 
        !infile = trim(path) // "flrw.in"
        !infile = 'flrw.in'
        !infile = trim(infile)
        !print*, "Phantom path is: ", path 
        !print*, "Infile is: ", infile
        ! Use system call to copy phantom files to simulation directory
        ! This is a digusting temporary fix
        !call SYSTEM('cp ~/phantomET/phantom/test/flrw* ./')

        ! The infile from ET
        infilestor = infilestart
        
        ! We should do everything that is done in phantom.f90
        
        ! Setup mpi
        id=0
        call init_mpi(id,nprocs)
        ! setup io 
        call set_io_unit_numbers
        ! routine that starts a phantom run
        print*, "Start run called!"
        ! Do we want to pass dt in here??
        call startrun(infilestor,logfilestor,evfilestor,dumpfilestor)
        print*, "Start run finished!"
        print*, "tmunugrid: ", tmunugrid(1,1,6,6,6)
        !stop 
        ! Intialises values for the evol routine: t, dt, etc..
        call evol_init(infilestor,logfilestor,evfilestor,dumpfilestor,dt_et)
        print*, "Evolve init finished!"
        ! Calculate the stress energy tensor for each particle
        ! Might be better to do this in evolve init 
        !call get_tmunugrid_all

        
    end subroutine init_et2phantom

    subroutine init_et2phantomgrid(nx,ny,nz,originx,originy,originz,dx,dy,dz)
        use einsteintk_utils
        integer,            intent(in) :: nx,ny,nz ! The maximum values of the grid in each dimension
        real(8),            intent(in) :: originx, originy, originz ! The origin of grid 
        real(8),            intent(in) :: dx, dy, dz ! Grid spacing in each dimension
        !integer,            intent(in) :: boundsizex, boundsizey, boundsizez

        ! Setup metric grid 
        call init_etgrid(nx,ny,nz,originx,originy,originz,dx,dy,dz)

    end subroutine init_et2phantomgrid

    subroutine init_phantom2et()
        ! Subroutine 
    end subroutine init_phantom2et

    subroutine et2phantom(rho,nx,ny,nz)
        integer, intent(in) :: nx, ny, nz
        real, intent(in) :: rho(nx,ny,nz)
        
        print*, "Grid limits: ", nx, ny, nz
        print*, "rho 1-10: ", rho(1:10,1,1)
        ! get mpi thread number
        ! send grid limits 
    end subroutine et2phantom

    subroutine step_et2phantom(infile,dt_et)
        use einsteintk_utils
        use evolve,          only:evol_step
        use tmunu2grid
        character(len=*),  intent(in) :: infile
        real,          intent(inout) :: dt_et
        character(len=500) :: logfile,evfile,dumpfile,path
         
       
        ! Print the values of logfile, evfile, dumpfile to check they are sensible
        !print*, "logfile, evfile, dumpfile: ", logfile, evfile, dumpfile
        print*, "stored values of logfile, evfile, dumpfile: ", logfilestor, evfilestor, dumpfilestor
        
        ! Interpolation stuff 
        ! Call et2phantom (construct global grid, metric, metric derivs, determinant)
        ! Run phantom for a step 
        call evol_step(infile,logfilestor,evfilestor,dumpfilestor,dt_et)
        ! Interpolation stuff back to et
        !call get_tmunugrid_all()
        ! call phantom2et (Tmunu_grid)
    
    end subroutine step_et2phantom
    
    subroutine phantom2et()
        ! should take in the cctk_array for tmunu??
        ! Is it better if this routine is just 
        ! Calculate stress energy tensor for each particle 

        ! Perform kernel interpolation from particles to grid positions 
    
    end subroutine phantom2et
end module einsteintk_wrapper 
