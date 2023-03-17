module einsteintk_wrapper
!
!
! This module is a "wrapper" for the hydro evol + communication with ET
! Subroutines here should be called by ET rather than calling phantom subroutines
! directly 
!
    implicit none 
    contains

    subroutine init_et2phantom(infilestart,dt_et,nophantompart,dtout)
        ! Wrapper that intialises phantom
        ! Intended to hide all of the inner works of phantom from ET
        ! Majority of the code from HelloHydro_init has been moved here 

        use io,              only:id,master,nprocs,set_io_unit_numbers,die
        use mpiutils,        only:init_mpi,finalise_mpi
        use initial,         only:initialise,finalise,startrun,endrun
        use evolve,          only:evol_init
        use tmunu2grid
        use einsteintk_utils
        use extern_gr
        use metric 
        use part, only:xyzh,vxyzu,dens,metricderivs, metrics, npart, tmunus
        

        implicit none
        character(len=*),  intent(in) :: infilestart
        real,          intent(in) :: dt_et 
        integer,       intent(inout) :: nophantompart
        real,          intent(out)   :: dtout
        !character(len=500) :: logfile,evfile,dumpfile,path
        integer :: i,j,k,pathstringlength
        integer :: xlower,ylower,zlower,xupper,yupper,zupper
        real :: pos(3), gcovpart(0:3,0:3)
        !real :: dtout
    
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
        !print*, "tmunugrid: ", tmunugrid(1,1,6,6,6)
        !stop 
        ! Intialises values for the evol routine: t, dt, etc..
        call evol_init(infilestor,logfilestor,evfilestor,dumpfilestor,dt_et,nophantompart)
        print*, "Evolve init finished!"
        ! Calculate the stress energy tensor for each particle
        ! Might be better to do this in evolve init 
        !call get_tmunugrid_all
        ! Calculate the stress energy tensor
        call get_metricderivs_all(dtout,dt_et) ! commented out to try and fix prim2cons
        !call get_tmunu_all(npart,xyzh,metrics,vxyzu,metricderivs,dens,tmunus) ! commented out to try and fix prim2cons
        !call get_tmunu_all_exact(npart,xyzh,metrics,vxyzu,metricderivs,dens,tmunus)
        ! Interpolate stress energy tensor from particles back 
        ! to grid
        !call get_tmunugrid_all(npart,xyzh,vxyzu,tmunus,calc_cfac=.true.) ! commented out to try and fix cons2prim

        call get_phantom_dt(dtout)
        
        
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

    subroutine step_et2phantom_MoL(infile,dt_et,dtout)
        use part, only:xyzh,vxyzu,fxyzu,pxyzu,dens,metricderivs, metrics, npart, tmunus,eos_vars
        use cons2prim, only: cons2primall
        use deriv
        use extern_gr
        use tmunu2grid
        use einsteintk_utils, only: get_phantom_dt,gcovgrid
        character(len=*),  intent(in) :: infile
        real,          intent(inout) :: dt_et
        real,          intent(out)   :: dtout
        real                         :: vbefore,vafter

        ! Metric should have already been passed in
        ! and interpolated
        ! Call get_derivs global
        call get_derivs_global
        
        ! Get metric derivs 
        call get_metricderivs_all(dtout,dt_et)
        ! Store our particle quantities somewhere / send them to ET
        ! Cons2prim after moving the particles with the external force
        vbefore = vxyzu(1,1)
        call cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars)
        vafter = vxyzu(1,1)

        ! Does get_derivs_global perform a stress energy calc??
        ! If not do that here

        ! Perform the calculation of the stress energy tensor 
        ! Interpolate the stress energy tensor back to the ET grid!
        ! Calculate the stress energy tensor
        call get_tmunu_all(npart,xyzh,metrics,vxyzu,metricderivs,dens,tmunus)
        ! Interpolate stress energy tensor from particles back 
        ! to grid
        call get_tmunugrid_all(npart,xyzh,vxyzu,tmunus) 
        call get_phantom_dt(dtout)


    end subroutine step_et2phantom_MoL

    subroutine et2phantom_tmunu()
        use part,   only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
        Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,&
        dustfrac,ddustevol,eos_vars,pxyzu,dens,metrics,dustevol,tmunus,metricderivs,&
        massoftype,igas,rhoh,alphaind,dvdx,gradh
        !use part, only:xyzh,vxyzu,fxyzu,pxyzu,dens,metricderivs, metrics, npart, tmunus,eos_vars
        use cons2prim, only: cons2primall
        use deriv
        use extern_gr
        use tmunu2grid
        use einsteintk_utils, only: get_phantom_dt,gcovgrid,rhostargrid,tmunugrid
        use metric_tools, only:init_metric
        use densityforce, only:densityiterate
        use linklist,     only:set_linklist

        real :: stressmax
        real(kind=16) :: cfac

        stressmax = 0. 

        ! Also probably need to pack the metric before I call things
        call init_metric(npart,xyzh,metrics)
         ! Might be better to just do this in get derivs global with a number 2 call?
        ! Rebuild the tree
        call set_linklist(npart,npart,xyzh,vxyzu)
        ! Apparently init metric needs to be called again???
        !call init_metric(npart,xyzh,metrics)
        ! Calculate the cons density
        call densityiterate(1,npart,npart,xyzh,vxyzu,divcurlv,divcurlB,Bevol,&
                        stressmax,fxyzu,fext,alphaind,gradh,rad,radprop,dvdx) 
        ! Get primative variables for tmunu 
        call cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars)

        call get_tmunu_all(npart,xyzh,metrics,vxyzu,metricderivs,dens,tmunus)
        ! Interpolate stress energy tensor from particles back 
        ! to grid
        call get_tmunugrid_all(npart,xyzh,vxyzu,tmunus)
        
        ! Interpolate density to grid
        call phantom2et_rhostar

        ! Density check vs particles
        call check_conserved_dens(rhostargrid,cfac)

        ! Correct Tmunu 
        tmunugrid = cfac*tmunugrid


    end subroutine et2phantom_tmunu

    subroutine phantom2et_consvar()
        use part,   only:npart,xyzh,vxyzu,fxyzu,fext,divcurlv,divcurlB,&
        Bevol,dBevol,rad,drad,radprop,dustprop,ddustprop,&
        dustfrac,ddustevol,eos_vars,pxyzu,dens,metrics,dustevol,tmunus,metricderivs,&
        massoftype,igas,rhoh,alphaind,dvdx,gradh
        use densityforce, only:densityiterate
        use metric_tools, only:init_metric
        use linklist,     only:set_linklist
        use einsteintk_utils, only:rhostargrid,pxgrid
        use tmunu2grid, only:check_conserved_dens

        real :: stressmax
        real(kind=16) :: cfac

        ! Init metric
        call init_metric(npart,xyzh,metrics)
        
        ! Might be better to just do this in get derivs global with a number 2 call?
        ! Rebuild the tree
        call set_linklist(npart,npart,xyzh,vxyzu)
        ! Apparently init metric needs to be called again???
        call init_metric(npart,xyzh,metrics)
        ! Calculate the cons density
        call densityiterate(1,npart,npart,xyzh,vxyzu,divcurlv,divcurlB,Bevol,&
                         stressmax,fxyzu,fext,alphaind,gradh,rad,radprop,dvdx) 
        
        ! Interpolate density to grid
        call phantom2et_rhostar
        
        ! Interpolate momentum to grid
        call phantom2et_momentum


        ! Conserved quantity checks + corrections

        ! Density check vs particles
        call check_conserved_dens(rhostargrid,cfac)

        ! Momentum check vs particles

        ! Correct momentum and Density 
        rhostargrid = cfac*rhostargrid
        pxgrid = cfac*pxgrid
        !entropygrid = cfac*entropygrid


    end subroutine phantom2et_consvar

    subroutine phantom2et_rhostar()
        use part, only:xyzh,vxyzu,fxyzu,pxyzu,dens,metricderivs, metrics, npart, tmunus,eos_vars,&
        igas, massoftype,rhoh
        use cons2prim, only: cons2primall
        use deriv
        use extern_gr
        use tmunu2grid
        use einsteintk_utils, only: get_phantom_dt,rhostargrid
        use metric_tools, only:init_metric
        real :: dat(npart), h, pmass,rho
        integer :: i  


        ! Get new cons density from new particle positions somehow (maybe)? 
        ! Set linklist to update the tree for neighbour finding
        ! Calculate the density for the new particle positions
        ! Call density iterate 

        ! Interpolate from particles to grid
        ! This can all go into its own function as it will essentially
        ! be the same thing for all quantites  
        ! get particle data 
        ! get rho from xyzh and rhoh 
        ! Get the conserved density on the particles 
        dat = 0. 
        do i=1, npart
            ! Get the smoothing length 
            h = xyzh(4,i)
            ! Get pmass
            pmass = massoftype(igas)
            rho = rhoh(h,pmass)
            dat(i) = rho
        enddo 
        rhostargrid = 0. 
        call interpolate_to_grid(rhostargrid,dat)

    end subroutine phantom2et_rhostar

    subroutine phantom2et_momentum()
        use part, only:xyzh,vxyzu,fxyzu,pxyzu,dens,metricderivs, metrics, npart, tmunus,eos_vars,&
        igas,massoftype,alphaind,dvdx,gradh
        use cons2prim, only: cons2primall
        use deriv
        use extern_gr
        use tmunu2grid
        use einsteintk_utils, only: get_phantom_dt,gcovgrid,pxgrid
        use metric_tools, only:init_metric
        real :: dat(3,npart)
        integer :: i 


        ! Pi is directly updated at the end of each MoL add  

        ! Interpolate from particles to grid 
        ! get particle data for the x component of momentum
        dat = 0. 
        do i=1, npart
            dat(1,i) = pxyzu(1,i)
            dat(2,i) = pxyzu(2,i)
            dat(3,i) = pxyzu(3,i)
        enddo 
        
        pxgrid = 0. 
        ! call interpolate 3d 
        ! In this case call it 3 times one for each vector component
        ! px component
        call interpolate_to_grid(pxgrid(1,:,:,:), dat(1,:))
        ! py component
        call interpolate_to_grid(pxgrid(2,:,:,:), dat(2,:))
        ! pz component 
        call interpolate_to_grid(pxgrid(3,:,:,:),dat(3,:))
        
        

    end subroutine phantom2et_momentum



    ! Subroutine for performing a phantom dump from einstein toolkit 
    subroutine et2phantom_dumphydro(time,dt_et)
        use cons2prim, only:cons2primall
        use part, only:npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars
        use einsteintk_utils
        use evwrite,          only:write_evfile,write_evlog
        use readwrite_dumps,  only:write_smalldump,write_fulldump
        use fileutils,        only:getnextfilename
        real, intent(in)  :: time, dt_et
        !character(len=20) :: logfile,evfile,dumpfile
        
        ! Call cons2prim since values are updated with MoL
        !call cons2primall(npart,xyzh,metrics,pxyzu,vxyzu,dens,eos_vars) 
        ! Write EV_file 
        call write_evfile(time,dt_et)

        evfilestor  = getnextfilename(evfilestor)
        logfilestor = getnextfilename(logfilestor)
        dumpfilestor = getnextfilename(dumpfilestor)

        !print*, "Evfile: ", evfilestor
        !print*, "logfile: ", logfilestor
        !print*, "dumpfle: ", dumpfilestor
        ! Write full dump
        call write_fulldump(time,dumpfilestor)

    end subroutine et2phantom_dumphydro

    ! Provides the RHS derivs for a particle at index i 
    subroutine phantom2et_rhs(index, vx,vy,vz,fx,fy,fz,e_rhs)
        use einsteintk_utils
        real, intent(inout) :: vx,vy,vz,fx,fy,fz, e_rhs
        integer, intent(in) :: index

        call get_particle_rhs(index,vx,vy,vz,fx,fy,fz,e_rhs)

    end subroutine phantom2et_rhs

    subroutine phantom2et_initial(index,x,y,z,px,py,pz,e)
        use einsteintk_utils
        real, intent(inout) :: x,y,z,px,py,pz,e
        integer, intent(in) :: index

        call get_particle_val(index,x,y,z,px,py,pz,e)

    end subroutine phantom2et_initial

    subroutine et2phantom_setparticlevars(index,x,y,z,px,py,pz,e)
        use einsteintk_utils
        real, intent(inout) :: x,y,z,px,py,pz,e
        integer, intent(in) :: index

        call set_particle_val(index,x,y,z,px,py,pz,e)

    end subroutine et2phantom_setparticlevars
    
    ! I really HATE this routine being here but it needs to be to fix dependency issues. 
    subroutine get_metricderivs_all(dtextforce_min,dt_et)
      use einsteintk_utils, only: metricderivsgrid
      use part, only:npart, xyzh,vxyzu,fxyzu,metrics,metricderivs,dens,fext
      use timestep, only:bignumber,C_force
      use extern_gr, only:get_grforce
      use metric_tools, only:pack_metricderivs
      real, intent(out) :: dtextforce_min
      real, intent(in)  :: dt_et
      integer :: i
      real :: pri,dtf

      pri = 0.
      dtextforce_min = bignumber
      
      !$omp parallel do default(none) &
      !$omp shared(npart, xyzh,metrics,metricderivs,vxyzu,dens,C_force,fext) &
      !$omp firstprivate(pri) & 
      !$omp private(i,dtf) & 
      !$omp reduction(min:dtextforce_min)
      do i=1, npart 
        call pack_metricderivs(xyzh(1:3,i),metricderivs(:,:,:,i))
        call get_grforce(xyzh(:,i),metrics(:,:,:,i),metricderivs(:,:,:,i), &
             vxyzu(1:3,i),dens(i),vxyzu(4,i),pri,fext(1:3,i),dtf)
        dtextforce_min = min(dtextforce_min,C_force*dtf)
      enddo 
      !$omp end parallel do 
      ! manually add v contribution from gr
    !    do i=1, npart
    !      !fxyzu(:,i) = fxyzu(:,i) + fext(:,i)
    !      vxyzu(1:3,i) = vxyzu(1:3,i) + fext(:,i)*dt_et
    !    enddo 
    end subroutine get_metricderivs_all


end module einsteintk_wrapper 
