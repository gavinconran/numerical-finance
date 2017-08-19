
module functions

    implicit none
    integer :: fevals_proc
    real(kind=8) :: k
    save

contains

    real(kind=8) function f(x)
        implicit none
        real(kind=8), intent(in) :: x 
        integer proc_num, ierr

        ! keep track of number of function evaluations by
        ! each process:
        fevals_proc = fevals_proc + 1
        
        f = 1.d0 + x**3 + sin(k*x)
        
    end function f

end module functions
