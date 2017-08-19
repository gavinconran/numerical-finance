
module functions

    implicit none
    integer :: gevals
    save

contains

    function g(x,ndim)

    implicit none
    integer, intent(in) :: ndim
    real(kind=8), intent(in) :: x(ndim)
    real(kind=8) :: g
    integer :: i

    g = 0.d0
    do i=1,ndim
        g = g + x(i)**2
        enddo

    gevals = gevals + 1

    end function g

end module functions
