module fb_findloc
    implicit none
    private
    public :: findloc

    interface findloc
        module procedure :: findloc_char_1d
        module procedure :: findloc_char_1d_dim
    end interface

contains

    function findloc_char_1d(array, value) result(loc)
        character(*), intent(in) :: array(:)
        character(*), intent(in) :: value
        integer :: loc(1)
        integer :: i

        loc = 0
        do i = 1, size(array)
            if (array(i) == value) then
                loc(1) = i
                return
            end if
        end do
    end function

    integer function findloc_char_1d_dim(array, value, dim) result(loc)
        character(*), intent(in) :: array(:)
        character(*), intent(in) :: value
        integer, intent(in) :: dim
        integer :: i

        loc = 0
        if (dim /= 1) return
        do i = 1, size(array)
            if (array(i) == value) then
                loc = i
                return
            end if
        end do
    end function

end module
