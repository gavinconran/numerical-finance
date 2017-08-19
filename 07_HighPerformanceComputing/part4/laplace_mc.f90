
program laplace_mc

!Random walk approximate solution to Laplace's equation u_{xx} + u{yy} = 0.

!Set demo==True to plot a few random walks interactively.
!With demo==False many more walks are used to estimate the solution.

!The boundary conditions are set for this test problem by evaluating the
!true solution u(x,y) = x^2 - y^2 of Laplace's equation on the
!boundary of the domain.

!Moreover, the exact solution to the discrete equations
!  U_{i-1,j} + U_{i+1,j} + U_{i,j-1} + U_{i,j+1} - 4u_{ij} = 0
!with boundary values obtained in this way is easily computed. It is simply
!given by evaluating the exact solution at the grid points,
!  U_{ij} = x_i^2 - y_j^2
!This is because the centered difference approximation to the second
!derivative is exact when applied to a quadratic function.

!This code implements a random walk on a lattice (rectangular grid) where in
!each step the walk goes in one of 4 directions based on the value of a
!random number that's uniformly distributed in [0,1].

    use mpi
    use mc_walk, only: many_walks, utrue, uboundary, nwalks, initialise, print_nwalks
    use random_util, only: init_random_seed

    implicit none
    integer, parameter :: ndim = 20
    real(kind=8), dimension(ndim) :: a, b
    real(kind=8) :: ax = 0.d0, bx = 1.d0, ay = 0.4, by = 1.d0
    real(kind=8) :: dx, dy, x0, y0, u_mc, error, u_sum_new, u_sum_old, u_true, u_mc_total
    integer :: i, seed1, n_total, nx = 19, ny = 11, max_steps, n_mc, n_success, i0, j0
    integer :: ierr, num_procs, proc_num, nerr, total_nwalks

    call MPI_INIT(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, proc_num, ierr)

    if (proc_num==0) then
        open(unit=25, file='mc_laplace_error.txt', status='unknown')
        if (num_procs == 1) then
            print *, "*** Error, this version requires more than 1 process  ***"
            nerr = 1
        endif
    
    endif

    ! if nerr == 1 then all processes must stop:
    call MPI_BCAST(nerr, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    if (nerr == 1) then
        ! Note that error message already printed by Process 0
        ! All processes must execute the MPI_FINALIZE 
        ! (Could also just have "go to 99" here.)
        call MPI_FINALIZE(ierr)
        stop
    endif

    
    

    ! initialise nwalks
    nwalks = 0

    ! problem description:
    dx = (bx-ax)/float(nx+1)
    dy = (by-ay)/float(ny+1)

    ! Try it out from a specific (x0,y0):
    x0 = 0.9
    y0 = 0.6

    i0 = nint((x0-ax)/dx)
    j0 = nint((y0-ay)/dy)

    ! shift (x0,y0) to a grid point if it wasn't already:
    x0 = ax + i0*dx
    y0 = ay + j0*dy

    u_true = utrue(x0,y0)

    seed1 = 12345
    seed1 = seed1 + 97*proc_num  ! unique for each process
    call init_random_seed(seed1)

    ! initialise nwalksArray that captures the random walks per process
    call initialise(num_procs)
        

    ! maximum number of step in each before giving up:
    max_steps = 100*max(nx, ny) !10

    ! initial number of Monte-Carlo walks to take:
    n_mc = 10

    call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success, 0)

    if (proc_num==0) then
    ! relative error gives estimate of number of correct digits:
        error = abs((u_mc - u_true) / u_true)
        print '(i10, es22.14, es22.14, i5)',n_success, u_mc, error
        write(25,'(i10,e23.15,e15.6)') n_success, u_mc, error
    endif


    ! start accumulating totals:
    u_mc_total = u_mc
    n_total = n_success

    ! loop to add successively double the number of points used:
    do i=1,12
        u_sum_old = u_mc_total * n_total
        call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success, 0)
        u_sum_new = u_mc * n_success
        n_total = n_total + n_success
        u_mc_total = (u_sum_old + u_sum_new) / n_total
        error = abs((u_mc_total - u_true) / u_true)  
        if (proc_num == 0)then 
            print '(i10,es22.14, es22.14, i5)',n_total, u_mc_total, error
            write(25,'(i10,e23.15,e15.6)') n_total, u_mc_total, error
        endif 
        n_mc = 2*n_mc   ! double number of trials for next iteration

        enddo

    ! -----------------------------------------
    ! Final code for Master (Processor 0):
    ! -----------------------------------------
    !print '("i: ", i5, ",: ", i5)', i, proc_num
    if (proc_num == 0)then 
        print '("Final approximation to u(x0,y0): ", es22.14)', u_mc_total
        print '("Total walks performed by all processes: ",i13)',  n_total
          
    endif

    if (proc_num /= 0) goto 99
    ! print nwalks per process and terminate all salve processes  
    call print_nwalks(num_procs)
    !call many_walks(i0, j0, max_steps, 1, u_mc, n_success, 1000000)

    99  continue   ! might jump to here if finished early
!print '("i: ", i5, ", terminating: ", i5)', i, proc_num
    call MPI_FINALIZE(ierr)
end program laplace_mc
