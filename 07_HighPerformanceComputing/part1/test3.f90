
program test3

    use mpi

    use quadrature, only: trapezoid
    use functions, only: f, fevals_proc, k

    implicit none
    real(kind=8) :: a,b,int_true, int_approx, int_result, dx_sub
    real(kind=8) , dimension(2) :: ab_sub
    real(kind=8), allocatable, dimension(:) :: intgrals

    integer :: proc_num, num_procs, ierr, n, fevals_total,j, jj, nerr, nsub, numsent, sender, nextsub
    integer :: fevals_proc_total
    integer, dimension(MPI_STATUS_SIZE) :: status
    

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)

    ! All processes set these values so we don't have to broadcast:
    k = 1.d3   ! functions module variable 
    a = 0.d0
    b = 2.d0
    int_true = (b-a) + (b**4 - a**4) / 4.d0 - (1.d0/k) * (cos(k*b) - cos(k*a))
    n = 1000

    ! Each process keeps track of number of fevals:
    fevals_proc = 0


    if (proc_num==0) then
        if (num_procs == 1) then
            print *, "*** Error, this version requires more than 1 process  ***"
            nerr = 1
        endif
    
    endif

    ! if nerr == 1 then all processes must stop:
    call MPI_BCAST(nerr, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (nerr == 1) then
        ! Note that error message already printed by Process 0
        ! All processes must execute the MPI_FINALIZE 
        ! (Could also just have "go to 99" here.)
        call MPI_FINALIZE(ierr)
        stop
        endif

    ! obtain from user the number of sub intervals
    if (proc_num==0) then
        print *, "How many subintervals? "
        read *, nsub
    endif

    call MPI_BCAST(nsub, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    

    ! -----------------------------------------
    ! code for Master (Processor 0):
    ! -----------------------------------------
    
    if (proc_num == 0) then

        print '("Using ",i3," processes")', num_procs
        print '("true integral: ", es22.14)', int_true
        print *, " "  ! blank line
        print '("fevals by Process ",i2,": ",i13)',  proc_num, fevals_proc

        allocate(intgrals(nsub))    ! to hold integral of each sub interval in MPI_RECV

        dx_sub = (b-a) / (nsub)

        numsent = 0 ! keep track of how many columns sent

        ! send the first batch to get all workers working:
        do j=1,min(num_procs-1,nsub)
            ab_sub(1) = a + (j-1)*dx_sub
            ab_sub(2) = a + j*dx_sub
            call MPI_SEND(ab_sub, 2, MPI_DOUBLE_PRECISION, j, j, &
                          MPI_COMM_WORLD, ierr)
            numsent = numsent + 1
        enddo

        do j=1,nsub
            call MPI_RECV(int_result, 1, MPI_DOUBLE_PRECISION, &
                        MPI_ANY_SOURCE, MPI_ANY_TAG, &
                        MPI_COMM_WORLD, status, ierr)

            sender = status(MPI_SOURCE)
            jj = status(MPI_TAG) 
            intgrals(jj) = int_result

            if (numsent < nsub) then
                ! still more work to do, the next sub will be sent and
                ! this index also used as the tag:
                nextsub = numsent + 1 
                

                ab_sub(1) = a + (nextsub-1)*dx_sub
                ab_sub(2) = a + nextsub*dx_sub
                
                call MPI_SEND(ab_sub, 2, MPI_DOUBLE_PRECISION,&
                            sender, nextsub, MPI_COMM_WORLD, ierr)
                numsent = numsent + 1
            else
                ! send an empty message with tag=0 to indicate this worker
                ! is done:
                call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,&
                            sender, 0, MPI_COMM_WORLD, ierr)
            endif

        enddo

    endif

    
    ! -----------------------------------------
    ! code for Workers (Processors 1, 2, ...):
    ! -----------------------------------------
    if (proc_num /= 0) then

        if (proc_num > nsub) go to 99   ! no work expected

        do while (.true.)
            ! repeat until message with tag==0 received...
            call MPI_RECV(ab_sub, 2, MPI_DOUBLE_PRECISION, &
                      0, MPI_ANY_TAG, &
                      MPI_COMM_WORLD, status, ierr)

            j = status(MPI_TAG)   ! this is the column number
                                  ! may not be proc_num in general

            if (j==0) then
                print '("total fevals by Process ",i2,": ",i13)',  proc_num, fevals_proc_total 
                go to 99    ! received "done" message
            endif
        
            int_approx = trapezoid(f, ab_sub(1), ab_sub(2), n)
        
            call MPI_SEND(int_approx, 1, MPI_DOUBLE_PRECISION, &
                    0, j, MPI_COMM_WORLD, ierr)

            fevals_proc_total = fevals_proc_total + fevals_proc
        enddo

    endif

    
    ! -----------------------------------------
    ! Final code for Master (Processor 0):
    ! -----------------------------------------
    
    if (proc_num==0) then 
        print '("Trapezoid approximation with ",i8," total points: ",es22.14)',&
            (nsub)*n, sum(intgrals)
        print '("Total number of fevals: ",i10)', nsub * n
    endif

    99  continue   ! might jump to here if finished early
    call MPI_FINALIZE(ierr)

end program test3
