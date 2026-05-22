!-------------------------------------------------------------------------------
! Logging Module
!-------------------------------------------------------------------------------
! This module provides a basic abiltiy to add logging statements throughout the
! program. The module is modeled on Python's logging module and provides five
! logging levels: log_debug, log_info, log_warning, log_error and log_critical.
!-------------------------------------------------------------------------------

module fb_logging
    use fb_parameters, only: stdout
    implicit none
    private

    ! Unit where logging is written
    integer :: lunit = stdout

    ! Current logging level
    integer :: current_level = 30

    ! Defined logging levels
    type :: log_level_enum
        integer :: critical = 50
        integer :: error    = 40
        integer :: warning  = 30
        integer :: info     = 20
        integer :: debug    = 10
        integer :: none     = 0
    end type
    type(log_level_enum), parameter :: log_levels = log_level_enum()

    ! Scratch area for writing formatted log messages
    character(len=512) :: log_buffer

    ! Public interface
    public  :: log_buffer,   &
               log_levels,   &
               log_debug,    &
               log_info,     &
               log_warning,  &
               log_error,    &
               log_critical, &
               set_log_unit, &
               set_log_level

contains

    subroutine set_log_level(level)
        integer, intent(in) :: level
        current_level = level
    end subroutine

    subroutine set_log_unit(unit)
        integer, intent(in) :: unit
        lunit = unit
    end subroutine

    subroutine log_debug(message)
        character(len=*), intent(in) :: message
        call log(log_levels%debug, "DEBUG: ", message)
    end subroutine

    subroutine log_info(message)
        character(len=*), intent(in) :: message
        call log(log_levels%info, "INFO: ", message)
    end subroutine

    subroutine log_warning(message)
        character(len=*), intent(in) :: message
        call log(log_levels%warning, "WARNING: ", message)
    end subroutine

    subroutine log_error(message)
        character(len=*), intent(in) :: message
        call log(log_levels%error, "ERROR: ", message)
    end subroutine

    subroutine log_critical(message)
        character(len=*), intent(in) :: message
        call log(log_levels%critical, "CRITICAL: ", message)
    end subroutine

    subroutine log(level, prefix, message)
        integer, intent(in) :: level
        character(len=*), intent(in) :: prefix, message
        if (current_level == log_levels%none) then
            return
        end if
        if (level >= current_level) then
            write(lunit,'(a)') prefix // message
        end if
    end subroutine

end module
