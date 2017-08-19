
module mc_walk

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

contains

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

subroutine many_walks(i0, j0, max_steps, n_mc, u_mc, n_success)
    
    
    implicit none
    integer, intent(in) :: n_mc, max_steps
    integer, intent(out) :: n_success
    real(kind=8), intent(out) :: u_mc
    integer, intent(in) :: i0, j0
    ! loacal variables
    integer :: i, j, k, iabort
    real(kind=8) :: ub, ub_sum

    ub_sum = 0.d0   ! to accumulate boundary values reached from all walks
    n_success = 0    ! to keep track of how many walks reached boundary

    do k=1,n_mc
        i = i0
        j = j0
        call random_walk(i0, j0, max_steps, ub, iabort)
        if (.NOT.isnan(ub)) then
            ! use this result unless walk didn't reach boundary
            ub_sum = ub_sum + ub
            n_success = n_success + 1
        endif
    enddo

    u_mc = ub_sum / n_success   ! average over successful walks

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

    ! increment nwalk each time random_walk is called
    nwalks = nwalks + 1

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
            ub = NaN
        endif
    enddo

end subroutine random_walk    

end module mc_walk
