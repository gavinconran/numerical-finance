
module mc_walk

    use mpi

    implicit none
    real(kind=8), parameter :: ax = 0.0d0
    real(kind=8), parameter :: bx = 1.0d0
    real(kind=8), parameter :: ay = 0.4d0
    real(kind=8), parameter :: by = 1.0d0
    integer, parameter :: nx = 19
    integer, parameter :: ny = 11
    real(kind=8), parameter :: dx = (bx - ax) / (nx+1)
    real(kind=8), parameter :: dy = (by - ay) / (ny+1)
    integer :: nwalks
    integer, allocatable, dimension(:) :: nwalksArray

contains

subroutine initialise(num_procs)
    ! initialses nwalksArray used to count random walks per process
    implicit none
    integer, intent(in) :: num_procs
    allocate(nwalksArray(num_procs))
    nwalksArray = 0

end subroutine initialise

subroutine print_nwalks(num_procs)
    ! prints out the number of nwalks per process
    implicit none
    integer, intent(in) :: num_procs
    integer :: i

    ! terminate all nwalks per processor
    do i=0,num_procs-1
        print '("Walks performed by Process",i2,": ",i10)',  i, nwalksArray(i)        
    enddo
    
end subroutine print_nwalks

function utrue(x, y)

    ! True solution for comparison, if known.

    implicit none
    real(kind=8), intent(in) :: x,y
    real(kind=8) :: utrue

    utrue = x**2 - y**2

end function utrue

function uboundary(x, y)

    ! Return u(x,y) assuming (x,y) is a boundary point.

    implicit none
    real(kind=8), intent(in) :: x,y
    real(kind=8) :: uboundary

    if ((x-ax)*(x-bx)*(y-ay)*(y-by) .ne. 0.d0) then
        print *, "*** Error -- called uboundary at non-boundary point"
        stop
        endif

    uboundary = utrue(x,y)   ! assuming we know this

end function uboundary

subroutine many_walks(i0, j0, max_steps, n_mc, u_mc, n_success, terminate)
     
    implicit none
    integer, intent(in) :: n_mc, max_steps
    integer, intent(out) :: n_success
    real(kind=8), intent(out) :: u_mc
    integer, intent(in) :: i0, j0, terminate
    ! loacal variables
    integer :: i, j, k, iabort, iabort_result, bottom = 1
    integer :: numsent = 0, sender, ierr, jj, proc_num, num_procs, nextsub=0, total_nwalks, mc_count=0
    real(kind=8) :: ub, ub_sum, ub_result
    integer, dimension(MPI_STATUS_SIZE) :: status
    

    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)

     nwalks=0
     
    ! -----------------------------------------
    ! code for Master (Processor 0):
    ! -----------------------------------------
    
    if (proc_num == 0) then

    ub_sum = 0.d0   ! to accumulate boundary values reached from all walks
    n_success = 0    ! to keep track of how many walks reached boundary
    ub_result = 0.d0  !variable for boundary value returned
    numsent = 0 ! to keep track of number of jobs sent to slave processes
    

    ! send the first batch to get all workers working:
    do j=1,min(num_procs-1,n_mc)
        call MPI_SEND(bottom, 1, MPI_INTEGER, j, j, &
                      MPI_COMM_WORLD, ierr)
        numsent = numsent + 1
    enddo

    do j=1,n_mc
        call MPI_RECV(ub_result, 1, MPI_DOUBLE_PRECISION, &
                    MPI_ANY_SOURCE, MPI_ANY_TAG, &
                    MPI_COMM_WORLD, status, ierr)

        sender = status(MPI_SOURCE) 
        
        if (.NOT.isnan(ub_result)) then
            ! use this result unless walk didn't reach boundary
            ub_sum = ub_sum + ub_result
            n_success = n_success + 1
            nwalksArray(sender) = nwalksArray(sender) + 1
        endif

        if (numsent < n_mc) then
                ! still more work to do, the next sub will be sent and
                nextsub = numsent + 1
                call MPI_SEND(bottom, 1, MPI_INTEGER,&
                        sender, nextsub, MPI_COMM_WORLD, ierr)
                numsent = numsent + 1

        else
            ! send an empty message with tag=0 to indicate this worker
            ! is done:
            call MPI_SEND(MPI_BOTTOM, 0, MPI_DOUBLE_PRECISION,&
                        sender, 0, MPI_COMM_WORLD, ierr)     
        endif
    enddo

    u_mc = ub_sum / n_success   ! average over successful walks
    endif

    ! -----------------------------------------
    ! code for Workers (Processors 1, 2, ...):
    ! -----------------------------------------
    if (proc_num /= 0) then
       ub_result = 0.d0
        do while (.true.)
            ! repeat until message with tag==0 received...
            call MPI_RECV(bottom, 1, MPI_INTEGER, &
                      MPI_ANY_SOURCE, MPI_ANY_TAG, &
                      MPI_COMM_WORLD, status, ierr)

            j = status(MPI_TAG)   ! this is the column number
                                  ! may not be proc_num in general

            if (j==0) then
                EXIT    ! received "done" message
            endif

            ! continue on as normal
            call random_walk(i0, j0, max_steps, ub_result, iabort)
        
            call MPI_SEND(ub_result, 1, MPI_DOUBLE_PRECISION, &
                    0, proc_num, MPI_COMM_WORLD, ierr)

        enddo    
    endif

end subroutine many_walks

subroutine random_walk(i0, j0, max_steps, ub, iabort)

    ! Take one random walk starting at (i0,j0) until we reach the boundary or
    ! exceed max_steps steps.
    ! Return the value at the boundary point reached, or nan if we failed.

    implicit none
    integer, intent(in) :: i0, j0, max_steps
    integer, intent(out) :: iabort
    real(kind=8), intent(out) :: ub

    ! local variables
    integer :: i, j, istep
    real(kind=4) :: r(max_steps)
    real(kind=8) :: xb, yb, NaN

    iabort=0 !set to ok
    
    ! starting point:
    i = i0
    j = j0

    ! generate as many random numbers as we could possibly need
    ! for this walk, since this is much faster than generating one at a time:
    call random_number(r)

    do istep=1,max_steps
    
        ! Take the next random step with equal probability in each direction:

        if (r(istep) < 0.25) then
            i = i-1   ! step left
        elseif (r(istep) < 0.5) then
            i = i+1   ! step right
        elseif (r(istep) < 0.75) then
            j = j-1   ! step down
        else   
            j = j+1   ! step up
        endif

        ! check if we hit the boundary:
        if (i*j*(nx+1-i)*(ny+1-j) == 0) then
            xb = ax + i*dx
            yb = ay + j*dy
            ub = uboundary(xb, yb)

            EXIT ! end the walk
        endif


        if (istep==(max_steps)) then
            NaN = 0.
            NaN = NaN / NaN
            iabort = 1
            ub = NaN
        endif
    enddo

    iabort = 0

end subroutine random_walk    

end module mc_walk
