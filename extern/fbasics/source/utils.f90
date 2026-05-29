module fb_utils
    ! A collection of simple utility functions

    use fb_logging
    use fb_parameters, only: wp, empty_int, empty_real, tol
    use iso_c_binding, only: c_char, c_int, c_null_char
    implicit none
    integer, parameter :: IK = kind(0)
    integer, parameter :: SK = kind("a")
    integer, parameter :: LK = kind(.false.)

    interface is_empty
        module procedure :: is_empty_i
        module procedure :: is_empty_r
        module procedure :: is_empty_s
    end interface
    interface to_str
        module procedure :: to_str_i
        module procedure :: to_str_r
        module procedure :: to_str_l
    end interface
    interface alloc_size
        module procedure :: asize_r1, asize_r2, asize_r3, asize_r4
        module procedure :: asize_i1, asize_i2, asize_i3, asize_i4
    end interface
    interface ones
        module procedure :: ones1, ones2, ones3, ones4
    end interface
    interface zeros
        module procedure :: zeros1, zeros2, zeros3, zeros4
    end interface
    interface flat
        module procedure :: flat_r2, flat_r3, flat_r4
        module procedure :: flat_i2, flat_i3, flat_i4
    end interface
    interface
        function cea_bindc_abort_recovery_is_active() bind(c)
            import c_int
            integer(c_int) :: cea_bindc_abort_recovery_is_active
        end function

        subroutine cea_bindc_abort_recover(message, message_len) bind(c)
            import c_char, c_int
            character(kind=c_char), intent(in) :: message(*)
            integer(c_int), value :: message_len
        end subroutine
    end interface

    private :: set_ios, &
               is_empty_i, is_empty_r, is_empty_s,         &
               to_str_i,   to_str_r,   to_str_l,           &
               ones1,      ones2,      ones3,      ones4,  &
               zeros1,     zeros2,     zeros3,     zeros4, &
               flat_r2,    flat_r3,    flat_r4,            &
               flat_i2,    flat_i3,    flat_i4,            &
               asize_r1, asize_r2, asize_r3, asize_r4,     &
               asize_i1, asize_i2, asize_i3, asize_i4

contains

    ! Returns `.true.` if it is a desired character.
    function isDesired(char) result(desired)
        character(1, SK), intent(in) :: char
        logical(LK) :: desired
        desired =   (SK_"0" <= char .and. char <= SK_"9") .or. &
                    (SK_"A" <= char .and. char <= SK_"Z") .or. &
                    (SK_"a" <= char .and. char <= SK_"z")
    end function

    function replace(str, isDesired) result(strrep)
        character(*, SK), intent(in)    :: str
        character(:, SK), allocatable   :: strrep
        procedure(logical(LK))          :: isDesired
        integer(IK)                     :: i, counter
        allocate(character(len(str), SK) :: strrep)
        counter = 0_IK
        do i = 1, len(str, kind = IK)
            if (.not. isDesired(str(i:i))) cycle
            counter = counter + 1_IK
            strrep(counter:counter) = str(i:i)
        end do
        strrep = strrep(1:counter)
    end function

    !-----------------------------------------------------------------------
    ! Program Control
    !-----------------------------------------------------------------------
    subroutine abort(msg)
        ! Exit with error code 1, optionally displaying error message.
        character(*), intent(in), optional :: msg
        character(kind=c_char), allocatable :: cmsg(:)
        integer :: i, nmsg

        if (present(msg)) call log_critical(msg)

        if (cea_bindc_abort_recovery_is_active() /= 0) then
            if (present(msg)) then
                nmsg = len_trim(msg)
                allocate(cmsg(nmsg + 1))
                do i = 1, nmsg
                    cmsg(i) = msg(i:i)
                end do
                cmsg(nmsg + 1) = c_null_char
                call cea_bindc_abort_recover(cmsg, int(nmsg, c_int))
            else
                allocate(cmsg(1))
                cmsg(1) = c_null_char
                call cea_bindc_abort_recover(cmsg, 0_c_int)
            end if
        end if

        stop 1
    end subroutine

    subroutine assert(cond,msg)
        ! Abort if condition is not true.
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        if (.not. cond) call abort(msg)
    end subroutine

    function getenv(var) result(val)
        ! Returns the specified environment variable; aborts on failure.

        character(*), intent(in) :: var
        character(:), allocatable :: val
        character(512) :: buffer
        integer :: length, stat

        call get_environment_variable(var,buffer,length,stat)
        if (stat==-1) call abort('fb-getenv: Environment variable '//trim(var)//' is too long (>512 chars).')
        if (stat== 1) call abort('fb-getenv: Environment variable '//trim(var)//' is not set.')
        if (stat== 2) call abort('fb-getenv: OS does not support environment variables.')
        val = buffer(1:length)

    end function

    !-----------------------------------------------------------------------
    ! Empty Checking
    !-----------------------------------------------------------------------
    function is_empty_i(val) result(tf)
        integer, intent(in) :: val
        logical :: tf
        tf = (val == empty_int)
    end function

    function is_empty_r(val) result(tf)
        real(wp), intent(in) :: val
        logical :: tf
        tf = (val >= empty_real) ! Use >= to avoid compiler warning
    end function

    function is_empty_s(str) result(tf)
        character(*), intent(in) :: str
        logical :: tf
        tf = (len_trim(replace(str, isDesired)) == 0)
    end function

    !-----------------------------------------------------------------------
    ! File Paths
    !-----------------------------------------------------------------------
    function basename(path) result(bname)
        ! Returns the non-path part of the file path, with extension
        character(*), intent(in)  :: path
        character(:), allocatable :: bname
        integer :: last_slash
        last_slash = scan(trim(path), "/", BACK=.true.)
        if (last_slash > 0) then
            bname = trim(path(last_slash+1:))
        else
            bname = path
        end if
    end function

    function dirname(path) result(dir)
        ! Returns the directory part of the file path, with trailing slash
        character(*), intent(in) :: path
        character(:), allocatable :: dir
        integer :: last_slash
        last_slash = scan(trim(path), "/", BACK=.true.)
        if (last_slash > 0) then
            dir = path(:last_slash)
        else
            dir = "./"
        end if
    end function

    logical function exists(path)
        ! Check if a file path exists
        character(*), intent(in) :: path
        inquire(file=trim(path), exist=exists)
    end function

    logical function is_readable(path)
        ! Check if a file path exists and is readable
        character(*), intent(in) :: path
        integer :: io, ios
        open(newunit=io, file=trim(path), status='old', action='read', iostat=ios)
        is_readable = (ios == 0)
        if (is_readable) close(io)
    end function

    logical function is_writeable(path)
        ! Check if a file path exists and is writeable
        character(*), intent(in) :: path
        integer :: io, ios
        open(newunit=io, file=trim(path), status='old', action='write', position='append', iostat=ios)
        is_writeable = (ios == 0)
        if (is_writeable) close(io)
    end function

    function locate(path_in,candidates) result(path_out)
        ! Return path of a file that may exists in one of several candidate directories
        !  1. Check if input path is already valid; if so, return it
        !  2. Check if file exists in candidate directories; return first match
        !  3. If no match, return empty string

        character(*), intent(in)  :: path_in
        character(*), intent(in), optional  :: candidates(:)
        character(:), allocatable :: path_out
        integer :: i

        ! Check if already valid
        if (exists(path_in)) then
            path_out = path_in
            return
        end if

        ! Check if in candidates
        if (present(candidates)) then
            do i = 1,size(candidates)
                path_out = trim(candidates(i)) // '/' // path_in
                if (exists(path_out)) return
            end do
        end if

        ! Return empty on failure
        path_out = ''

    end function

    !-----------------------------------------------------------------------
    ! String Analysis / Manipulation
    !-----------------------------------------------------------------------
    function substring(str,substr,istart) result(is_substring)
        ! True if input contains substring
        character(*), intent(in) :: str, substr
        integer, intent(out), optional :: istart  ! Starting index of substr
        logical :: is_substring
        integer :: is

        is = index(str,substr)
        is_substring = (is /= 0)
        if (present(istart)) istart=is

    end function

    logical function startswith(str,substr)
        ! Check if string starts with given substring
        character(*), intent(in) :: str, substr
        integer :: n,m

        m = len(str)
        n = len(substr)
        startswith = (n <= m)
        if (.not. startswith) return
        startswith = (str(1:n) == substr)

    end function

    logical function endswith(str,substr)
        ! Check if string ends with given substring
        character(*), intent(in) :: str, substr
        integer :: n,m

        m  = len(str)
        n  = len(substr)
        endswith = (n <= m)
        if (.not. endswith) return
        endswith = (str(m-n+1:) == substr)

    end function

    function lower(str) result (lstr)
        ! Return a lower case copy of the input string

        character(*), intent(in) :: str
        character(len(str))      :: lstr
        integer :: ic, i

        character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

        ! Lowercase each letter if it is capitalized
        lstr = str
        do i = 1,len_trim(str)
            ic = index(cap, str(i:i))
            if (ic > 0) lstr(i:i) = low(ic:ic)
        end do

    end function

    function upper(str) result (ustr)
        ! Return an upper case copy of the input string

        character(*), intent(in) :: str
        character(len(str))      :: ustr
        integer :: ic, i

        character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

        ! Capitalize each letter if it is lowecase
        ustr = str
        do i = 1,len_trim(str)
            ic = index(low, str(i:i))
            if (ic > 0) ustr(i:i) = cap(ic:ic)
        end do

    end function

    !-----------------------------------------------------------------------
    ! String Formatting/Parsing
    !-----------------------------------------------------------------------
    function to_str_i(val) result(s)
        ! Return string representation of integer
        integer, intent(in) :: val
        character(:), allocatable :: s
        character(256) :: buffer
        write(buffer,'(i0)') val
        s = trim(buffer)
    end function

    function to_str_r(val,digits) result(s)
        ! Return string representation of integer
        real(wp), intent(in) :: val
        integer, intent(in), optional :: digits
        character(:), allocatable :: fmt, s
        character(256) :: buffer
        fmt = '(1p,g0.6)'
        if (present(digits)) then
            fmt = '(1p,g0.'//to_str(digits)//')'
        end if
        write(buffer,fmt) val
        s = trim(buffer)
    end function

    function to_str_l(val) result(s)
        ! Return string representation of integer
        logical, intent(in) :: val
        character(:), allocatable :: s
        character(256) :: buffer
        write(buffer,'(l4)') val
        s = trim(buffer)
    end function

    function to_int(str,ios) result(val)
        ! Convert a string-formatted integer into integer variable
        character(*), intent(in) :: str
        integer, optional, intent(out) :: ios
        integer :: val
        integer :: ios_
        val = empty_int
        read(str,*,iostat=ios_) val
        call set_ios(ios_, ios, 'Could not parse string as int.')
    end function

    function to_real(str,ios) result(val)
        ! Convert a string-formatted real into real variable
        character(*), intent(in) :: str
        integer, optional, intent(out) :: ios
        real(wp) :: val
        integer :: ios_
        val = empty_real
        read(str,*,iostat=ios_) val
        call set_ios(ios_, ios, 'Could not parse string as real')
    end function

    function to_logical(str,ios) result(val)
        ! Convert a string-formatted logical into logical variable
        character(*), intent(in) :: str
        integer, optional, intent(out) :: ios
        logical :: val
        integer :: ios_
        val = .false.
        read(str,*,iostat=ios_) val
        call set_ios(ios_, ios, 'Could not parse string as logical')
    end function

    !-----------------------------------------------------------------------
    ! Arrays
    !-----------------------------------------------------------------------
    function asize_r1(arr,dim) result(asize)
        ! Check size of a potentially un-allocated array.
        ! Returns -1 is array is unallocated. Based on:
        ! https://fortran-lang.discourse.group/t/size-of-an-unallocated-array/4633/2
        real(wp), intent(in), allocatable :: arr(:)
        integer, intent(in), optional :: dim
        integer :: asize
        asize = -1
        if (allocated(arr)) asize = size(arr,dim)
    end function

    function asize_r2(arr,dim) result(asize)
        real(wp), intent(in), allocatable :: arr(:,:)
        integer, intent(in), optional :: dim
        integer :: asize
        asize = -1
        if (allocated(arr)) asize = size(arr,dim)
    end function

    function asize_r3(arr,dim) result(asize)
        real(wp), intent(in), allocatable :: arr(:,:,:)
        integer, intent(in), optional :: dim
        integer :: asize
        asize = -1
        if (allocated(arr)) asize = size(arr,dim)
    end function

    function asize_r4(arr,dim) result(asize)
        real(wp), intent(in), allocatable :: arr(:,:,:,:)
        integer, intent(in), optional :: dim
        integer :: asize
        asize = -1
        if (allocated(arr)) asize = size(arr,dim)
    end function

    function asize_i1(arr,dim) result(asize)
        integer, intent(in), allocatable :: arr(:)
        integer, intent(in), optional :: dim
        integer :: asize
        asize = -1
        if (allocated(arr)) asize = size(arr,dim)
    end function

    function asize_i2(arr,dim) result(asize)
        integer, intent(in), allocatable :: arr(:,:)
        integer, intent(in), optional :: dim
        integer :: asize
        asize = -1
        if (allocated(arr)) asize = size(arr,dim)
    end function

    function asize_i3(arr,dim) result(asize)
        integer, intent(in), allocatable :: arr(:,:,:)
        integer, intent(in), optional :: dim
        integer :: asize
        asize = -1
        if (allocated(arr)) asize = size(arr,dim)
    end function

    function asize_i4(arr,dim) result(asize)
        integer, intent(in), allocatable :: arr(:,:,:,:)
        integer, intent(in), optional :: dim
        integer :: asize
        asize = -1
        if (allocated(arr)) asize = size(arr,dim)
    end function

    function flat_r2(arr) result(flat)
        real(wp), intent(in)  :: arr(:,:)
        real(wp) :: flat(size(arr))
        flat = pack(arr,.true.)
    end function

    function flat_r3(arr) result(flat)
        real(wp), intent(in)  :: arr(:,:,:)
        real(wp) :: flat(size(arr))
        flat = pack(arr,.true.)
    end function

    function flat_r4(arr) result(flat)
        real(wp), intent(in)  :: arr(:,:,:,:)
        real(wp) :: flat(size(arr))
        flat = pack(arr,.true.)
    end function

    function flat_i2(arr) result(flat)
        integer, intent(in)  :: arr(:,:)
        integer :: flat(size(arr))
        flat = pack(arr,.true.)
    end function

    function flat_i3(arr) result(flat)
        integer, intent(in)  :: arr(:,:,:)
        integer :: flat(size(arr))
        flat = pack(arr,.true.)
    end function

    function flat_i4(arr) result(flat)
        integer, intent(in)  :: arr(:,:,:,:)
        integer :: flat(size(arr))
        flat = pack(arr,.true.)
    end function

    function ones1(n1) result(x)
        ! Returns array filled with 1's
        integer, intent(in) :: n1
        real(wp) :: x(n1)
        x = 1.0d0
    end function

    function ones2(n1,n2) result(x)
        ! Returns array filled with 1's
        integer, intent(in) :: n1,n2
        real(wp) :: x(n1,n2)
        x = 1.0d0
    end function

    function ones3(n1,n2,n3) result(x)
        ! Returns array filled with 1's
        integer, intent(in) :: n1,n2,n3
        real(wp) :: x(n1,n2,n3)
        x = 1.0d0
    end function

    function ones4(n1,n2,n3,n4) result(x)
        ! Returns array filled with 1's
        integer, intent(in) :: n1,n2,n3,n4
        real(wp) :: x(n1,n2,n3,n4)
        x = 1.0d0
    end function

    function zeros1(n1) result(x)
        ! Returns array filled with 0's
        integer, intent(in) :: n1
        real(wp) :: x(n1)
        x = 0.0d0
    end function

    function zeros2(n1,n2) result(x)
        ! Returns array filled with 0's
        integer, intent(in) :: n1,n2
        real(wp) :: x(n1,n2)
        x = 0.0d0
    end function

    function zeros3(n1,n2,n3) result(x)
        ! Returns array filled with 0's
        integer, intent(in) :: n1,n2,n3
        real(wp) :: x(n1,n2,n3)
        x = 0.0d0
    end function

    function zeros4(n1,n2,n3,n4) result(x)
        ! Returns array filled with 0's
        integer, intent(in) :: n1,n2,n3,n4
        real(wp) :: x(n1,n2,n3,n4)
        x = 0.0d0
    end function

    !-----------------------------------------------------------------------
    ! Helper Functions
    !-----------------------------------------------------------------------
    subroutine set_ios(value, ios, message)
        ! Set an optional iostat flag. If no ios and value /= 0, aborts.
        integer, intent(in) :: value
        integer, intent(out), optional :: ios
        character(*), intent(in), optional :: message
        character(:), allocatable :: msg
        if (present(ios)) then
            ios = value
        else if (value /= 0) then
            msg = 'Operation failed. '
            if (present(message)) msg = msg // message
            call abort(msg)
        end if
    end subroutine


end module
