
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


    use mc_walk, only: many_walks, utrue, uboundary, nwalks
    use random_util, only: init_random_seed

    implicit none
    integer, parameter :: ndim = 20
    real(kind=8), dimension(ndim) :: a, b
    real(kind=8) :: ax = 0.d0, bx = 1.d0, ay = 0.4, by = 1.d0
    real(kind=8) :: dx, dy, x0, y0, u_mc, error, u_sum_new, u_sum_old, u_true, u_mc_total
    integer :: i, seed1, n_total, nx = 19, ny = 11, max_steps, n_mc, n_success, i0, j0

    open(unit=25, file='mc_laplace_error.txt', status='unknown')

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

    !print '("True solution of PDE: u(",es22.14,", ", es22.14") = ",es22.14)',  x0, y0, u_true
    !print '("Note: with solution used in demo this is also the solution to the")'
    !print '("      the finite-difference equations on the same grid.")'


    seed1 = 12345   ! or set to 0 for random seed
    call init_random_seed(seed1)

    ! maximum number of step in each before giving up:
    max_steps = 100*max(nx, ny) !10

    ! initial number of Monte-Carlo walks to take:
    n_mc = 10
    
    call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)

    ! relative error gives estimate of number of correct digits:
    error = abs((u_mc - u_true) / u_true)
    print '(i10, es22.14, es22.14)',n_success, u_mc, error
    write(25,'(i10,e23.15,e15.6)') n_success, u_mc, error


    ! start accumulating totals:
    u_mc_total = u_mc
    n_total = n_success

    ! loop to add successively double the number of points used:

    do i=1,12

        u_sum_old = u_mc_total * n_total
        call many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)
        u_sum_new = u_mc * n_success
        n_total = n_total + n_success
        u_mc_total = (u_sum_old + u_sum_new) / n_total
        error = abs((u_mc_total - u_true) / u_true)

        print '(i10,es22.14, es22.14)',n_total, u_mc_total, error
        write(25,'(i10,e23.15,e15.6)') n_total, u_mc_total, error
        n_mc = 2*n_mc   ! double number of trials for next iteration

        enddo
    print '("Final approximation to u(x0,y0): ", es22.14)', u_mc_total
    print '("Total number of random walks:", i5)', nwalks
end program laplace_mc
