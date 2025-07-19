!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2025 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module easter_egg
!
! easter_egg
!
! :References: None
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: io, random
!
 implicit none

 public :: bring_the_egg

 logical, public, save :: egged = .false.

 private
 integer :: iseed  = -1234
 logical :: backed = .false.
 logical :: cleared

contains

subroutine bring_the_egg
 integer, allocatable :: grid_real(:,:),grid_mask(:,:)
 integer, allocatable :: grid_mem(:,:,:)
 integer :: ndim,idir,iter,istore
 logical :: stop_game,mem_stored
 character(len=20) :: command,quit_in

 call init(ndim,grid_real,grid_mask,grid_mem)
 call write_to_shell(ndim,grid_real)
 stop_game  = .false.
 mem_stored = .false.
 iter = 0

 tile_smash: do while (.not.stop_game)
    iter = iter + 1
    idir = 0
    !- Interpret input
    read(*,'(a)') command
    command = trim(adjustl(command))
    if (command == 'q') then
       write(*,'(a)',advance='no') 'Quit game (y/n)? '
       read(*,'(a)') quit_in
       if (trim(adjustl(quit_in)) == 'y') exit
    elseif (command == 'b') then
       if (mem_stored) then
          if (.not.backed) then
             backed = .true.
             !- fetch the second-last stored grid
             istore = 1
             if (mod(iter,2) == 0) istore = 2
             grid_real(1:ndim,1:ndim) = grid_mem(istore,1:ndim,1:ndim)
             call write_to_shell(ndim,grid_real)
             cycle tile_smash
          else
             print*,'You can only go back once @_@'
             cycle tile_smash
          endif
       else
          print*,'Go back to what ?_?'
          cycle tile_smash
       endif
    else
       backed = .false.
       if (command == 'a') then
          idir = 0
       elseif (command == 'w') then
          idir = 1
       elseif (command == 'd') then
          idir = 2
       elseif (command == 's') then
          idir = 3
       else
          print*,'da heck you mean'
          cycle tile_smash
       endif
    endif
    !- Smash tile
    call update_grid(ndim,grid_real,grid_mask,idir)

    !- Output and store to memory
    call write_to_shell(ndim,grid_real)
    istore = 1
    if (mod(iter,2) == 0) istore = 2
    grid_mem(istore,1:ndim,1:ndim) = grid_real(1:ndim,1:ndim)
    mem_stored = .true.

    call check_end(ndim,grid_real,stop_game)
 enddo tile_smash
 print*,'- GET BACK TO WORK YOU SLACKER -'
 call sleep(1)

end subroutine bring_the_egg

 !
 ! Initialise grid
 !
subroutine init(ndim,grid_real,grid_mask,grid_mem)
 integer, allocatable, intent(inout) :: grid_real(:,:)
 integer, allocatable, intent(inout) :: grid_mask(:,:)
 integer, allocatable, intent(inout) :: grid_mem(:,:,:)
 integer, intent(out) :: ndim
 integer :: i,j
 character(len=20) :: ndim_in

 call print_2048
 print*,'Intructions: w - up; s - down; a - left; d - right; b - back; q - exit'
 print*,''
 write(*,'(a)',advance='no') 'Enter dimension: [default 4] '
 read(*,'(a)') ndim_in
 if (len(trim(adjustl(ndim_in))) == 0) then
    ndim = 4
 else
    read(ndim_in,*) ndim
    ndim = max(ndim,3)
    ndim = min(ndim,10)
 endif
 allocate(grid_real(ndim,ndim))
 allocate(grid_mask(ndim,ndim))
 allocate(grid_mem(2,ndim,ndim))
 do i = 1,ndim
    do j = 1,ndim
       grid_real(i,j) = -1
    enddo
 enddo
 call add_extra(ndim,grid_real)

 cleared = .false.

end subroutine init

 !
 ! Checks to end the game
 !
subroutine check_end(ndim,grid,end_game)
 integer, intent(in) :: ndim
 integer, intent(in) :: grid(ndim,ndim)
 logical, intent(inout) :: end_game
 integer :: i,j,num
 character(len=20) :: cont_in

 end_game = .true.

 !- Check 2048 or empty cells
 each_cell: do i = 1,ndim
    do j = 1,ndim
       if (grid(i,j) == 2048) then
          if (.not.cleared) then
             print*,'- YOU WON -'
             write(*,'(a)',advance='no') 'Continue game ([y]/n)? '
             read(*,'(a)') cont_in
             end_game = .false.
             if (trim(adjustl(cont_in)) == 'n') then
                end_game = .true.
                exit each_cell
             else
                cleared = .true.
                print*,'You are the king of procrastination.'
             endif
          endif
       elseif (grid(i,j) == -1) then
          end_game = .false.
          exit each_cell
       endif
    enddo
 enddo each_cell

 !- Check if further moves are possible
 do i = 1,ndim
    do j = 1,ndim
       num = grid(i,j)
       if (i == 1) then
          if (j == 1) then        ! only i+1 and j+1
             if (grid(i+1,j) == num .or. grid(i,j+1) == num) end_game = .false.
          elseif (j == ndim) then ! only i+1 and j-1
             if (grid(i+1,j) == num .or. grid(i,j-1) == num) end_game = .false.
          endif
       elseif (i == ndim) then
          if (j == 1) then        ! only i-1 and j+1
             if (grid(i-1,j) == num .or. grid(i,j+1) == num) end_game = .false.
          elseif (j == ndim) then ! only i-1 and j-1
             if (grid(i-1,j) == num .or. grid(i,j-1) == num) end_game = .false.
          endif
       else
          if (grid(i+1,j) == num .or. grid(i-1,j) == num .or. &
              grid(i,j+1) == num .or. grid(i,j-1) == num) end_game = .false.
       endif
    enddo
 enddo

end subroutine check_end

 !
 ! Routine to add a cell of 2 or 4 at random location
 !
subroutine add_extra(ndim,grid)
 use random, only:ran2
 use io,     only:fatal
 integer, intent(in)    :: ndim
 integer, intent(inout) :: grid(ndim,ndim)
 integer :: num_to_add,niter,i,j
 real    :: rani,ranj,prob
 real    :: probof2 = 0.9
 logical :: got_it

 niter = 0
 got_it = .false.
 do while (.not.got_it)
    rani = ran2(iseed)
    ranj = ran2(iseed)
    i = floor(rani*ndim)+1
    j = floor(ranj*ndim)+1
    i = min(i,ndim)
    j = min(j,ndim)
    if (grid(i,j) == -1) then
       prob = ran2(iseed)
       num_to_add = 4
       if (prob < probof2) num_to_add = 2
       grid(i,j) = num_to_add
       got_it = .true.
    endif
    niter = niter + 1
    if (niter > 100) call fatal('egg.f90','broken egg :( ')
 enddo

end subroutine add_extra

 !
 ! Main routine-
 ! Fix the direction coords within smashing code; Rotate the
 ! grid before going into smash, then revert back to original.
 !
subroutine update_grid(ndim,grid_real,grid_mask,idir)
 integer, intent(in) :: ndim,idir
 integer, intent(inout) :: grid_mask(ndim,ndim)
 integer, intent(inout) :: grid_real(ndim,ndim)
 integer :: i,j
 integer :: grid_prev(ndim,ndim)
 logical :: can_add

 grid_prev = grid_real
 call rotate(ndim,grid_real,grid_mask,idir)
 call smash(ndim,grid_mask)
 call backrotate(ndim,grid_mask,grid_real,idir)

 !- add extra 2 only if grid has changed
 can_add = .false.
 do i = 1,ndim
    do j = 1,ndim
       if (grid_real(i,j) /= grid_prev(i,j)) then
          can_add = .true.
          exit
       endif
    enddo
 enddo
 if (can_add) call add_extra(ndim,grid_real)

end subroutine update_grid

 !
 ! Rotate grid to make it swiping right to left
 !
subroutine rotate(ndim,grid_in,grid_out,idir)
 integer, intent(in) :: ndim,idir
 integer, intent(in) :: grid_in(ndim,ndim)
 integer, intent(out) :: grid_out(ndim,ndim)
 integer :: num,i,j

 if (idir == 1) then !- 90deg anticlockwise
    do i = 1,ndim
       do j = 1,ndim
          num = grid_in(j,(ndim+1)-i)
          grid_out(i,j) = num
       enddo
    enddo
 elseif (idir == 2) then !- 180deg flip
    do i = 1,ndim
       do j = 1,ndim
          num = grid_in(i,(ndim+1)-j)
          grid_out(i,j) = num
       enddo
    enddo
 elseif (idir == 3) then !- 90deg clockwise
    do i = 1,ndim
       do j = 1,ndim
          num = grid_in((ndim+1)-j,i)
          grid_out(i,j) = num
       enddo
    enddo
 else
    grid_out = grid_in
 endif

end subroutine rotate

 !
 ! Rotate grid back to original position
 !
subroutine backrotate(ndim,grid_in,grid_out,idir)
 integer, intent(in) :: ndim,idir
 integer, intent(in) :: grid_in(ndim,ndim)
 integer, intent(out) :: grid_out(ndim,ndim)
 integer :: num,i,j

 if (idir == 1) then !- 90deg clockwise
    do i = 1,ndim
       do j = 1,ndim
          num = grid_in((ndim+1)-j,i)
          grid_out(i,j) = num
       enddo
    enddo
 elseif (idir == 2) then !- 180deg flip
    do i = 1,ndim
       do j = 1,ndim
          num = grid_in(i,(ndim+1)-j)
          grid_out(i,j) = num
       enddo
    enddo
 elseif (idir == 3) then !- 90deg anticlockwise
    do i = 1,ndim
       do j = 1,ndim
          num = grid_in(j,(ndim+1)-i)
          grid_out(i,j) = num
       enddo
    enddo
 else
    grid_out = grid_in
 endif

end subroutine backrotate

 !
 ! Smash a given grid from right to left
 !
subroutine smash(ndim,grid)
 integer, intent(in) :: ndim
 integer, intent(inout) :: grid(ndim,ndim)
 integer :: i,j
 logical :: merged(ndim,ndim)
 logical :: no_more_moves

 !- marks whether a cell had already been merged with another cell
 !- follows each non-empty cell in grid
 do i = 1,ndim
    do j = 1,ndim
       merged(i,j) = .false.
    enddo
 enddo

 no_more_moves = .false.
 do while (.not.no_more_moves)
    no_more_moves = .true.
    each_row: do i = 1,ndim
       each_cell: do j = 2,ndim
          if (grid(i,j) /= -1) then
             if (grid(i,j-1) == -1) then  !-shift to the right
                no_more_moves = .false.
                grid(i,j-1) = grid(i,j)
                merged(i,j-1) = merged(i,j)
                if (j < ndim) then
                   merged(i,j) = merged(i,j+1)
                else
                   merged(i,j) = .false.
                endif
                grid(i,j) = -1
             elseif (grid(i,j) == grid(i,j-1)) then !-merge
                if (merged(i,j).eqv..false.) then
                   no_more_moves = .false.
                   grid(i,j-1) = grid(i,j-1)*2
                   merged(i,j-1) = .true.
                   grid(i,j) = -1
                endif
             endif
          endif
       enddo each_cell
    enddo each_row
 enddo

end subroutine smash

 !
 ! Output current grid_real
 !
subroutine write_to_shell(ndim,grid)
 integer, intent(in) :: ndim
 integer, intent(in) :: grid(ndim,ndim)
 integer :: i,j

 do i = 1,ndim
    write(*,'(3x,a)',advance='no') ''
    do j = 1,ndim
       write(*,'(a)',advance='no') '+-------'
    enddo
    write(*,'(a)',advance='no') '+'
    write(*,'(a)') ''
    write(*,'(3x,a)',advance='no') '|'
    do j = 1,ndim
       write(*,'(6x,a2)',advance='no') ' |'
    enddo
    write(*,'(a)') ''
    write(*,'(3x,a)',advance='no') '|'
    do j = 1,ndim
       if (grid(i,j) /= -1) then
          if (grid(i,j) < 999) then
             write(*,'(i5,1x,a2)',advance='no') grid(i,j), ' |'
          else
             write(*,'(i7,a1)',advance='no') grid(i,j), '|'
          endif
       else
          write(*,'(6x,a2)',advance='no') ' |'
       endif
    enddo
    write(*,'(a)') ''
    write(*,'(3x,a)',advance='no') '|'
    do j = 1,ndim
       write(*,'(6x,a2)',advance='no') ' |'
    enddo
    write(*,'(a)') ''
 enddo
 write(*,'(3x,a)',advance='no') ''
 do j = 1,ndim
    write(*,'(a)',advance='no') '+-------'
 enddo
 write(*,'(a)',advance='no') '+'
 write(*,'(a)') ''

end subroutine write_to_shell

subroutine print_2048
 write(*,'(a)') ''
 write(*,'(a)') ''
 write(*,'(a)') ' .-----.    .-----.     -----.   .------. '
 write(*,'(a)') '/ /`` \ \  / / ``\ \   / / | |  / / `` \ \ '
 write(*,'(a)') '--     | | | |   | |  / /  | |  | |     | | '
 write(*,'(a)') '      / /  | |   | | / /___| |_  \ \ __/ / '
 write(*,'(a)') '    / /    | |   | | |_____   _| / /    \ \ '
 write(*,'(a)') ' ./ /_____ \ \__/  /       | |   | |     | | '
 write(*,'(a)') '|________|  \___ _/        |_|   \ \ ___/ /  '
 write(*,'(a)') '.                                  .-----.  of Phantom'
 write(*,'(a)') ''
end subroutine print_2048

end module easter_egg
