
module random_util

contains

subroutine init_random_seed(seed1)

    ! Seed the random number generator.
    ! If seed1 = 0 then set seed1 randomly using the system clock.
    !              This will give different sequences of random numbers
    !              from different runs, so results are not reproducible.
    ! If seed1 > 0 then results will be reproducible since calling this
    !              twice with the same seed will initialize the random
    !              number generator so the same sequence is generated.
    ! Once seed1 is set, set the other elements of the seed array by adding
    ! multiples of 37 as suggested in the documentation.
    ! The length of the seed array is determined by calling  
    !    
    integer, intent(inout) :: seed1

    integer :: nseed, clock 
    integer, dimension(:), allocatable :: seed

    call random_seed(size = nseed)  ! determine how many numbers needed to seed
    allocate(seed(nseed)) 

    if (seed1 == 0) then
        ! randomize the seed: not repeatable
        call system_clock(count = clock)
        seed1 = clock
      endif

    do i=1,nseed
        seed(i) = seed1 + 37*(i-1)  
        enddo

    print *, "seed1 for random number generator:", seed1

    call random_seed(put = seed)   ! seed the generator
    deallocate(seed)

end subroutine init_random_seed

end module random_util
