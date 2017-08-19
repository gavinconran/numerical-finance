
program test_quad_mc

    use functions, only: g, gevals
    use quadrature_mc, only: quad_mc
    use random_util, only: init_random_seed

    implicit none
    integer, parameter :: ndim = 20
    real(kind=8), dimension(ndim) :: a, b
    real(kind=8) :: volume, int_true, int_mc, int_mc_total, error
    integer :: i, seed1, n_total
    integer :: num_chunks, points_per_chunk, npoints

    gevals = 0
    open(unit=25, file='mc_quad_error.txt', status='unknown')

    do i=1,ndim
        a(i) = 2.d0
        b(i) = 4.d0
        enddo

    ! compute the true integral for special case where
    ! g(x) = sum(x**2) over all dimensions  -- think about why this works!
    volume = product(b-a)  ! =  product of b(i)-a(i) of ndim dimensions
    int_true = volume * sum((b**3 - a**3) / (3.d0*(b-a)))

    print '("Testing Monte Carlo quadrature in ",i2," dimensions")', ndim
    print '("True integral: ", es22.14)', int_true


    seed1 = 12345   ! or set to 0 for random seed
    call init_random_seed(seed1)

    ! Start with Monte Carlo using only a few points:
    npoints = 10
    int_mc = quad_mc(g,a,b,ndim,npoints)

    ! start accumulating totals:
    int_mc_total = int_mc
    n_total = npoints

    ! relative error gives estimate of number of correct digits:
    error = abs(int_mc_total - int_true) / abs(int_true)
    write(25,'(i10,e23.15,e15.6)') n_total, int_mc_total, error

    ! loop to add successively double the number of points used:

    do i=1,17
        ! compute approximation with new chunk of points
        int_mc = quad_mc(g,a,b,ndim,npoints)

        ! update estimate of integral based on doubling number:
        int_mc_total = 0.5*(int_mc_total + int_mc)
        n_total = n_total + npoints

        ! relative error:
        error = abs(int_mc_total - int_true) / abs(int_true)

        write(25,'(i10,e23.15,e15.6)') n_total, int_mc_total, error

        npoints = 2*npoints
        enddo

    print '("Final approximation to integral: ",es22.14)',int_mc_total
    print *, "Total g evaluations: ",gevals

end program test_quad_mc
