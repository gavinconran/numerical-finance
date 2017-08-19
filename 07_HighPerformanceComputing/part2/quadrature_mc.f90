module quadrature_mc

contains

real(kind=8) function quad_mc(g, a, b, ndim, npoints)
    
    ! Estimate the  Monte Carlo approximation to the integral of g(x) from a to b using the

    ! Input:
    !   g:  the function to integrate
    !   a:  left endpoint array of lower limits of integration in each dimension
    !   b:  right endpoint array of upper limits of integration in each dimension
    !   ndim:  number of dimensions to integrate over.
    !   npoints: number of Monte Carlo samples to use
    ! Returns:
    !    Monte Carlo approximation to the integral
    implicit none
    integer, intent(in) :: ndim, npoints
    real(kind=8), intent(in) :: a(ndim), b(ndim)
    real(kind=8), external :: g

    !local variables
    real (kind=8) :: Volume = 0.d0, monte_sum
    integer :: j
    real(kind=8) :: point(ndim)
    !real(kind=4), dimension(ndim, npoints) :: r
    !allocate(r(ndim:npoints))
    real(kind=4) :: r(npoints*ndim)
    
    Volume = product(b - a) / 2.d0
    call random_number(r)
    do j=1,npoints
        point = a + (r(j:) * (b-a))
        monte_sum = monte_sum + (g(point, ndim))
        enddo

    quad_mc = Volume * (monte_sum/(npoints)) !* 13.99700730437646 !scaling factor but why?

    end function quad_mc

end module quadrature_mc
