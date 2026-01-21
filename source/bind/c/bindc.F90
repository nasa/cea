#include "cea_enum.h"
module cea_bindc
    use cea, snl => species_name_len, &
             enl => element_name_len, &
             wp => real_kind
    use cea_param, only: empty_dp, gas_constant, get_data_search_dirs
    use cea_input, only: ReactantInput
    use cea_mixture, only: names_match
    use iso_c_binding
    use fb_logging
    use fb_utils, only: assert, locate, is_empty, to_str
    implicit none


    !-----------------------------------------------------------------
    ! Enumerations (mirror source/bind/c/cea_enum.h)
    !-----------------------------------------------------------------
    enum, bind(c)
        enumerator :: CEA_TP = 0
        enumerator :: CEA_HP = 1
        enumerator :: CEA_SP = 2
        enumerator :: CEA_TV = 3
        enumerator :: CEA_UV = 4
        enumerator :: CEA_SV = 5
    end enum

    enum, bind(c)
        enumerator :: CEA_NUM_REACTANTS = 0
        enumerator :: CEA_NUM_PRODUCTS = 1
        enumerator :: CEA_NUM_GAS = 2
        enumerator :: CEA_NUM_CONDENSED = 3
        enumerator :: CEA_NUM_ELEMENTS = 4
        enumerator :: CEA_MAX_EQUATIONS = 5
    end enum

    enum, bind(c)
        enumerator :: CEA_TEMPERATURE = 0
        enumerator :: CEA_PRESSURE = 1
        enumerator :: CEA_VOLUME = 2
        enumerator :: CEA_DENSITY = 3
        enumerator :: CEA_M = 4
        enumerator :: CEA_MW = 5
        enumerator :: CEA_ENTHALPY = 6
        enumerator :: CEA_ENERGY = 7
        enumerator :: CEA_ENTROPY = 8
        enumerator :: CEA_GIBBS_ENERGY = 9
        enumerator :: CEA_GAMMA_S = 10
        enumerator :: CEA_FROZEN_CP = 11
        enumerator :: CEA_FROZEN_CV = 12
        enumerator :: CEA_EQUILIBRIUM_CP = 13
        enumerator :: CEA_EQUILIBRIUM_CV = 14
        enumerator :: CEA_VISCOSITY = 15
        enumerator :: CEA_FROZEN_CONDUCTIVITY = 16
        enumerator :: CEA_EQUILIBRIUM_CONDUCTIVITY = 17
        enumerator :: CEA_FROZEN_PRANDTL = 18
        enumerator :: CEA_EQUILIBRIUM_PRANDTL = 19
    end enum

    enum, bind(c)
        enumerator :: CEA_ROCKET_TEMPERATURE = 0
        enumerator :: CEA_ROCKET_PRESSURE = 1
        enumerator :: CEA_ROCKET_VOLUME = 2
        enumerator :: CEA_ROCKET_DENSITY = 3
        enumerator :: CEA_ROCKET_M = 4
        enumerator :: CEA_ROCKET_MW = 5
        enumerator :: CEA_ROCKET_ENTHALPY = 6
        enumerator :: CEA_ROCKET_ENERGY = 7
        enumerator :: CEA_ROCKET_ENTROPY = 8
        enumerator :: CEA_ROCKET_GIBBS_ENERGY = 9
        enumerator :: CEA_ROCKET_GAMMA_S = 10
        enumerator :: CEA_ROCKET_FROZEN_CP = 11
        enumerator :: CEA_ROCKET_FROZEN_CV = 12
        enumerator :: CEA_ROCKET_EQUILIBRIUM_CP = 13
        enumerator :: CEA_ROCKET_EQUILIBRIUM_CV = 14
        enumerator :: CEA_MACH = 15
        enumerator :: CEA_SONIC_VELOCITY = 16
        enumerator :: CEA_AE_AT = 17
        enumerator :: CEA_C_STAR = 18
        enumerator :: CEA_COEFFICIENT_OF_THRUST = 19
        enumerator :: CEA_ISP = 20
        enumerator :: CEA_ISP_VACUUM = 21
        enumerator :: CEA_ROCKET_VISCOSITY = 22
        enumerator :: CEA_ROCKET_FROZEN_CONDUCTIVITY = 23
        enumerator :: CEA_ROCKET_EQUILIBRIUM_CONDUCTIVITY = 24
        enumerator :: CEA_ROCKET_FROZEN_PRANDTL = 25
        enumerator :: CEA_ROCKET_EQUILIBRIUM_PRANDTL = 26
    end enum

    enum, bind(c)
        enumerator :: CEA_SHOCK_TEMPERATURE = 0
        enumerator :: CEA_SHOCK_PRESSURE = 1
        enumerator :: CEA_SHOCK_VELOCITY = 2
        enumerator :: CEA_SHOCK_MACH = 3
        enumerator :: CEA_SHOCK_SONIC_VELOCITY = 4
        enumerator :: CEA_SHOCK_RHO12 = 5
        enumerator :: CEA_SHOCK_RHO52 = 6
        enumerator :: CEA_SHOCK_P21 = 7
        enumerator :: CEA_SHOCK_P52 = 8
        enumerator :: CEA_SHOCK_T21 = 9
        enumerator :: CEA_SHOCK_T52 = 10
        enumerator :: CEA_SHOCK_M21 = 11
        enumerator :: CEA_SHOCK_M52 = 12
        enumerator :: CEA_SHOCK_V2 = 13
        enumerator :: CEA_SHOCK_U5_P_V2 = 14
        enumerator :: CEA_SHOCK_VOLUME = 15
        enumerator :: CEA_SHOCK_DENSITY = 16
        enumerator :: CEA_SHOCK_M = 17
        enumerator :: CEA_SHOCK_MW = 18
        enumerator :: CEA_SHOCK_ENTHALPY = 19
        enumerator :: CEA_SHOCK_ENERGY = 20
        enumerator :: CEA_SHOCK_ENTROPY = 21
        enumerator :: CEA_SHOCK_GIBBS_ENERGY = 22
        enumerator :: CEA_SHOCK_GAMMA_S = 23
        enumerator :: CEA_SHOCK_FROZEN_CP = 24
        enumerator :: CEA_SHOCK_FROZEN_CV = 25
        enumerator :: CEA_SHOCK_EQUILIBRIUM_CP = 26
        enumerator :: CEA_SHOCK_EQUILIBRIUM_CV = 27
        enumerator :: CEA_SHOCK_VISCOSITY = 28
        enumerator :: CEA_SHOCK_FROZEN_CONDUCTIVITY = 29
        enumerator :: CEA_SHOCK_EQUILIBRIUM_CONDUCTIVITY = 30
        enumerator :: CEA_SHOCK_FROZEN_PRANDTL = 31
        enumerator :: CEA_SHOCK_EQUILIBRIUM_PRANDTL = 32
    end enum

    enum, bind(c)
        enumerator :: CEA_DETONATION_P1 = 0
        enumerator :: CEA_DETONATION_T1 = 1
        enumerator :: CEA_DETONATION_H1 = 2
        enumerator :: CEA_DETONATION_M1 = 3
        enumerator :: CEA_DETONATION_GAMMA1 = 4
        enumerator :: CEA_DETONATION_V_SONIC1 = 5
        enumerator :: CEA_DETONATION_PRESSURE = 6
        enumerator :: CEA_DETONATION_TEMPERATURE = 7
        enumerator :: CEA_DETONATION_DENSITY = 8
        enumerator :: CEA_DETONATION_ENTHALPY = 9
        enumerator :: CEA_DETONATION_ENERGY = 10
        enumerator :: CEA_DETONATION_GIBBS_ENERGY = 11
        enumerator :: CEA_DETONATION_ENTROPY = 12
        enumerator :: CEA_DETONATION_MACH = 13
        enumerator :: CEA_DETONATION_VELOCITY = 14
        enumerator :: CEA_DETONATION_SONIC_VELOCITY = 15
        enumerator :: CEA_DETONATION_GAMMA = 16
        enumerator :: CEA_DETONATION_P_P1 = 17
        enumerator :: CEA_DETONATION_T_T1 = 18
        enumerator :: CEA_DETONATION_M_M1 = 19
        enumerator :: CEA_DETONATION_RHO_RHO1 = 20
        enumerator :: CEA_DETONATION_FROZEN_CP = 21
        enumerator :: CEA_DETONATION_FROZEN_CV = 22
        enumerator :: CEA_DETONATION_EQUILIBRIUM_CP = 23
        enumerator :: CEA_DETONATION_EQUILIBRIUM_CV = 24
        enumerator :: CEA_DETONATION_M = 25
        enumerator :: CEA_DETONATION_MW = 26
        enumerator :: CEA_DETONATION_VISCOSITY = 27
        enumerator :: CEA_DETONATION_FROZEN_CONDUCTIVITY = 28
        enumerator :: CEA_DETONATION_EQUILIBRIUM_CONDUCTIVITY = 29
        enumerator :: CEA_DETONATION_FROZEN_PRANDTL = 30
        enumerator :: CEA_DETONATION_EQUILIBRIUM_PRANDTL = 31
    end enum

    enum, bind(c)
        enumerator :: CEA_SUCCESS = 0
        enumerator :: CEA_INVALID_FILENAME = 1
        enumerator :: CEA_INVALID_PROPERTY_TYPE = 2
        enumerator :: CEA_INVALID_EQUILIBRIUM_TYPE = 3
        enumerator :: CEA_INVALID_ROCKET_TYPE = 4
        enumerator :: CEA_INVALID_EQUILIBRIUM_SIZE_TYPE = 5
        enumerator :: CEA_INVALID_INDEX = 6
        enumerator :: CEA_INVALID_SIZE = 7
        enumerator :: CEA_NOT_CONVERGED = 8
    end enum

    !-----------------------------------------------------------------
    ! C Interop Helpers
    !-----------------------------------------------------------------
    interface c_len
        module procedure :: c_len_cstr
        module procedure :: c_len_cptrs
    end interface

    interface c_copy
        module procedure :: c_copy_str_cstr
        module procedure :: c_copy_str_cptr
    end interface

    interface to_str
        module procedure :: to_str_cptr
    end interface

    ! Reactant input for custom species
    type, bind(c) :: cea_reactant_input
        type(c_ptr) :: name = c_null_ptr
        integer(c_int) :: num_elements = 0
        type(c_ptr) :: elements = c_null_ptr
        type(c_ptr) :: coefficients = c_null_ptr
        logical(c_bool) :: has_molecular_weight = .false.
        real(c_double) :: molecular_weight = 0.0d0
        logical(c_bool) :: has_enthalpy = .false.
        real(c_double) :: enthalpy = 0.0d0
        type(c_ptr) :: enthalpy_units = c_null_ptr
        logical(c_bool) :: has_temperature = .false.
        real(c_double) :: temperature = 0.0d0
        type(c_ptr) :: temperature_units = c_null_ptr
    end type

    ! Custom types for optional argument support
    type, bind(c) :: cea_solver_opts
        real(c_double) :: trace = -1.0d0
        logical(c_bool) :: ions = .false.
        logical(c_bool) :: transport = .false.
        type(c_ptr) :: reactants = c_null_ptr
        integer(c_int) :: ninsert = 0
        type(c_ptr) :: insert = c_null_ptr
    end type

    !-----------------------------------------------------------------
    ! Module Global Data
    !-----------------------------------------------------------------
    ! NOTE: Initialization is not thread safe. Assumes one thread
    ! (root proc) calls cea_init_* before other threads start working.
    type(ThermoDB) :: global_thermodb
    type(TransportDB) :: global_transdb
    logical :: thermo_initialized = .false.
    logical :: trans_initialized = .false.
    character(:), allocatable :: thermo_path
    character(:), allocatable :: trans_path

contains

    function cea_solver_opts_init(opts) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(cea_solver_opts), intent(out) :: opts
        ierr = CEA_SUCCESS
        opts%trace = -1.0d0
        opts%ions = .false.
        opts%transport = .false.
        opts%reactants = c_null_ptr
        opts%ninsert = 0
        opts%insert = c_null_ptr
    end function

    function cea_species_name_len(name_len) result(ierr) bind(c)
        integer(c_int) :: ierr
        integer(c_int), intent(out) :: name_len
        ierr = CEA_SUCCESS
        name_len = snl
    end function

    !-----------------------------------------------------------------
    ! Version Information
    !-----------------------------------------------------------------
    function cea_version_major(major) result(ierr) bind(c)
        integer(c_int) :: ierr
        integer(c_int), intent(out) :: major
        ierr  = CEA_SUCCESS
        major = version_major
    end function

    function cea_version_minor(minor) result(ierr) bind(c)
        integer(c_int) :: ierr
        integer(c_int), intent(out):: minor
        ierr  = CEA_SUCCESS
        minor = version_minor
    end function

    function cea_version_patch(patch) result(ierr) bind(c)
        integer(c_int) :: ierr
        integer(c_int), intent(out) :: patch
        ierr  = CEA_SUCCESS
        patch = version_patch
    end function


    !-----------------------------------------------------------------
    ! Logging Control
    !-----------------------------------------------------------------
    function cea_set_log_level(level) result(ierr) bind(c)
        integer(c_int) :: ierr
        integer(c_int), intent(in), value :: level
        ierr = CEA_SUCCESS
        call set_log_level(level)
    end function

    !-----------------------------------------------------------------
    ! Initialization (Not thread safe!)
    !-----------------------------------------------------------------
    function cea_init() result(ierr) bind(c)
        ! Load default thermo.lib and trans.lib into globals
        integer(c_int) :: ierr
        if (thermo_initialized) then
            ierr = CEA_SUCCESS
            return
        end if
        ierr = cea_init_thermo('thermo.lib'//c_null_char)
        if (ierr /= CEA_SUCCESS) return
        ierr = cea_init_trans('trans.lib'//c_null_char)
    end function

    function cea_init_thermo(cthermofile) result(ierr) bind(c)
        ! Initialize global ThermoDB using specified file

        integer(c_int) :: ierr
        character(c_char), intent(in) :: cthermofile(*)
        character(:), allocatable :: thermofile
        character(:), allocatable :: resolved
        character(:), allocatable :: search_dirs(:)

        ierr = CEA_SUCCESS

        call c_copy(cthermofile, thermofile)
        call get_data_search_dirs(search_dirs)
        resolved = locate(thermofile, search_dirs)
        if (is_empty(resolved)) then
            ierr = CEA_INVALID_FILENAME
            call log_error('Could not locate thermo database file: '//thermofile)
            return
        end if
        if (thermo_initialized) then
            if (resolved /= thermo_path) then
                ierr = CEA_INVALID_FILENAME
                call log_error('Thermo database already initialized with: '//thermo_path)
                return
            end if
            return
        end if

        global_thermodb = read_thermo(resolved)
        thermo_path = resolved
        thermo_initialized = .true.

    end function

    function cea_init_trans(ctransfile) result(ierr) bind(c)
        ! Initialize global TransDB using specified file

        integer(c_int) :: ierr
        character(c_char), intent(in) :: ctransfile(*)
        character(:), allocatable :: transfile
        character(:), allocatable :: resolved
        character(:), allocatable :: search_dirs(:)

        ierr = CEA_SUCCESS

        call c_copy(ctransfile, transfile)
        call get_data_search_dirs(search_dirs)
        resolved = locate(transfile, search_dirs)
        if (is_empty(resolved)) then
            ierr = CEA_INVALID_FILENAME
            call log_error('Could not locate transport database file: '//transfile)
            return
        end if
        if (trans_initialized) then
            if (resolved /= trans_path) then
                ierr = CEA_INVALID_FILENAME
                call log_error('Transport database already initialized with: '//trans_path)
                return
            end if
            return
        end if

        global_transdb = read_transport(resolved)
        trans_path = resolved
        trans_initialized = .true.

    end function

    function cea_is_initialized(initialized) result(ierr) bind(c)
        integer(c_int) :: ierr
        integer(c_int), intent(out) :: initialized
        ierr = CEA_SUCCESS
        if (thermo_initialized) then
            initialized = 1
        else
            initialized = 0
        end if
    end function


    !-----------------------------------------------------------------
    ! Mixture
    !-----------------------------------------------------------------
    function cea_mixture_create(mptr, nspecies, cspecies) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: mptr
        integer(c_int), value :: nspecies
        type(c_ptr), intent(in) :: cspecies(*)
        type(Mixture), pointer :: mix
        character(snl) :: species(nspecies)
        character(:), allocatable :: name
        integer :: n

        ierr = CEA_SUCCESS

        ! Convert names from C strings to fortran strings
        do n = 1,nspecies
            call c_copy(cspecies(n), name)
            if (len(name) > snl) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            species(n) = name
        end do

        ! Construct object
        allocate(mix)
        mix = Mixture(global_thermodb, species)
        mptr = c_loc(mix)
        call log_info('BINDC: Created Mixture object at '//to_str(mptr))

    end function

    function cea_mixture_create_w_ions(mptr, nspecies, cspecies) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: mptr
        integer(c_int), value :: nspecies
        type(c_ptr), intent(in) :: cspecies(*)
        type(Mixture), pointer :: mix
        character(snl) :: species(nspecies)
        character(:), allocatable :: name
        integer :: n

        ierr = CEA_SUCCESS

        ! Convert names from C strings to fortran strings
        do n = 1,nspecies
            call c_copy(cspecies(n), name)
            if (len(name) > snl) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            species(n) = name
        end do

        ! Construct object
        allocate(mix)
        mix = Mixture(global_thermodb, species, ions=.true.)
        mptr = c_loc(mix)
        call log_info('BINDC: Created Mixture object at '//to_str(mptr))

    end function

    function cea_mixture_create_from_reactants(mptr, nreac, creac, nomit, comit) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: mptr
        integer(c_int), value    :: nreac
        type(c_ptr), intent(in)  :: creac(*)
        integer(c_int), value    :: nomit
        type(c_ptr), intent(in)  :: comit(*)
        type(Mixture), pointer :: mix
        character(snl) :: reac(nreac), omit(nomit)
        character(:), allocatable :: name
        integer :: n

        ierr = CEA_SUCCESS

        ! Convert names from C strings to fortran strings
        do n = 1,nreac
            call c_copy(creac(n), name)
            if (len(name) > snl) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            reac(n) = name
        end do
        do n = 1,nomit
            call c_copy(comit(n), name)
            if (len(name) > snl) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            omit(n) = name
        end do

        ! Construct object
        allocate(mix)
        mix = Mixture(                    &
            thermo = global_thermodb,     &
            reactant_names = reac,        &
            omitted_product_names = omit  &
        )
        mptr = c_loc(mix)
        call log_info('BINDC: Create Mixture object at '//to_str(mptr))

    end function

    function cea_mixture_create_from_reactants_w_ions(mptr, nreac, creac, nomit, comit) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: mptr
        integer(c_int), value    :: nreac
        type(c_ptr), intent(in)  :: creac(*)
        integer(c_int), value    :: nomit
        type(c_ptr), intent(in)  :: comit(*)
        type(Mixture), pointer :: mix
        character(snl) :: reac(nreac), omit(nomit)
        character(:), allocatable :: name
        integer :: n

        ierr = CEA_SUCCESS

        ! Convert names from C strings to fortran strings
        do n = 1,nreac
            call c_copy(creac(n), name)
            if (len(name) > snl) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            reac(n) = name
        end do
        do n = 1,nomit
            call c_copy(comit(n), name)
            if (len(name) > snl) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            omit(n) = name
        end do

        ! Construct object
        allocate(mix)
        mix = Mixture(                    &
            thermo = global_thermodb,     &
            reactant_names = reac,        &
            omitted_product_names = omit, &
            ions= .true. &
        )
        mptr = c_loc(mix)
        call log_info('BINDC: Create Mixture object at '//to_str(mptr))

    end function

    function cea_mixture_create_from_input_reactants(mptr, nreac, creac) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: mptr
        integer(c_int), value :: nreac
        type(cea_reactant_input), intent(in) :: creac(*)
        type(Mixture), pointer :: mix
        type(ReactantInput), allocatable :: input_reactants(:)
        type(c_ptr), pointer :: c_elements(:)
        real(c_double), pointer :: c_coeffs(:)
        character(:), allocatable :: name, units, elem
        integer :: i, j, ne

        ierr = CEA_SUCCESS
        if (nreac <= 0) then
            ierr = CEA_INVALID_SIZE
            return
        end if

        allocate(input_reactants(nreac))

        do i = 1, nreac
            if (.not. c_associated(creac(i)%name)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            call c_copy(creac(i)%name, name)
            if (len_trim(name) == 0 .or. len(name) > snl) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            input_reactants(i)%name = name

            ne = creac(i)%num_elements
            if (ne < 0) then
                ierr = CEA_INVALID_SIZE
                return
            end if

            if (ne > 0) then
                if (.not. c_associated(creac(i)%elements) .or. .not. c_associated(creac(i)%coefficients)) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                allocate(input_reactants(i)%formula)
                allocate(input_reactants(i)%formula%elements(ne))
                allocate(input_reactants(i)%formula%coefficients(ne))
                call c_f_pointer(creac(i)%elements, c_elements, [ne])
                call c_f_pointer(creac(i)%coefficients, c_coeffs, [ne])
                do j = 1, ne
                    if (.not. c_associated(c_elements(j))) then
                        ierr = CEA_INVALID_SIZE
                        return
                    end if
                    call c_copy(c_elements(j), elem)
                    if (len_trim(elem) == 0 .or. len(elem) > enl) then
                        ierr = CEA_INVALID_SIZE
                        return
                    end if
                    input_reactants(i)%formula%elements(j) = elem
                    input_reactants(i)%formula%coefficients(j) = c_coeffs(j)
                end do
            else
                if (.not. thermodb_has_species(name)) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
            end if

            if (creac(i)%has_molecular_weight) then
                allocate(input_reactants(i)%molecular_weight)
                input_reactants(i)%molecular_weight = creac(i)%molecular_weight
            end if

            if (creac(i)%has_enthalpy) then
                if (.not. c_associated(creac(i)%enthalpy_units)) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                call c_copy(creac(i)%enthalpy_units, units)
                if (len_trim(units) == 0) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                allocate(input_reactants(i)%enthalpy)
                input_reactants(i)%enthalpy%name = 'h'
                input_reactants(i)%enthalpy%units = units
                allocate(input_reactants(i)%enthalpy%values(1))
                input_reactants(i)%enthalpy%values(1) = creac(i)%enthalpy
            end if

            if (creac(i)%has_temperature) then
                if (.not. c_associated(creac(i)%temperature_units)) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                call c_copy(creac(i)%temperature_units, units)
                if (len_trim(units) == 0) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                allocate(input_reactants(i)%temperature)
                input_reactants(i)%temperature%name = 't'
                input_reactants(i)%temperature%units = units
                allocate(input_reactants(i)%temperature%values(1))
                input_reactants(i)%temperature%values(1) = creac(i)%temperature
            end if
        end do

        allocate(mix)
        mix = Mixture(global_thermodb, input_reactants=input_reactants)
        mptr = c_loc(mix)
        call log_info('BINDC: Created Mixture object at '//to_str(mptr))

    end function

    function cea_mixture_create_from_input_reactants_w_ions(mptr, nreac, creac) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: mptr
        integer(c_int), value :: nreac
        type(cea_reactant_input), intent(in) :: creac(*)
        type(Mixture), pointer :: mix
        type(ReactantInput), allocatable :: input_reactants(:)
        type(c_ptr), pointer :: c_elements(:)
        real(c_double), pointer :: c_coeffs(:)
        character(:), allocatable :: name, units, elem
        integer :: i, j, ne

        ierr = CEA_SUCCESS
        if (nreac <= 0) then
            ierr = CEA_INVALID_SIZE
            return
        end if

        allocate(input_reactants(nreac))

        do i = 1, nreac
            if (.not. c_associated(creac(i)%name)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            call c_copy(creac(i)%name, name)
            if (len_trim(name) == 0 .or. len(name) > snl) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            input_reactants(i)%name = name

            ne = creac(i)%num_elements
            if (ne < 0) then
                ierr = CEA_INVALID_SIZE
                return
            end if

            if (ne > 0) then
                if (.not. c_associated(creac(i)%elements) .or. .not. c_associated(creac(i)%coefficients)) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                allocate(input_reactants(i)%formula)
                allocate(input_reactants(i)%formula%elements(ne))
                allocate(input_reactants(i)%formula%coefficients(ne))
                call c_f_pointer(creac(i)%elements, c_elements, [ne])
                call c_f_pointer(creac(i)%coefficients, c_coeffs, [ne])
                do j = 1, ne
                    if (.not. c_associated(c_elements(j))) then
                        ierr = CEA_INVALID_SIZE
                        return
                    end if
                    call c_copy(c_elements(j), elem)
                    if (len_trim(elem) == 0 .or. len(elem) > enl) then
                        ierr = CEA_INVALID_SIZE
                        return
                    end if
                    input_reactants(i)%formula%elements(j) = elem
                    input_reactants(i)%formula%coefficients(j) = c_coeffs(j)
                end do
            else
                if (.not. thermodb_has_species(name)) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
            end if

            if (creac(i)%has_molecular_weight) then
                allocate(input_reactants(i)%molecular_weight)
                input_reactants(i)%molecular_weight = creac(i)%molecular_weight
            end if

            if (creac(i)%has_enthalpy) then
                if (.not. c_associated(creac(i)%enthalpy_units)) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                call c_copy(creac(i)%enthalpy_units, units)
                if (len_trim(units) == 0) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                allocate(input_reactants(i)%enthalpy)
                input_reactants(i)%enthalpy%name = 'h'
                input_reactants(i)%enthalpy%units = units
                allocate(input_reactants(i)%enthalpy%values(1))
                input_reactants(i)%enthalpy%values(1) = creac(i)%enthalpy
            end if

            if (creac(i)%has_temperature) then
                if (.not. c_associated(creac(i)%temperature_units)) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                call c_copy(creac(i)%temperature_units, units)
                if (len_trim(units) == 0) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                allocate(input_reactants(i)%temperature)
                input_reactants(i)%temperature%name = 't'
                input_reactants(i)%temperature%units = units
                allocate(input_reactants(i)%temperature%values(1))
                input_reactants(i)%temperature%values(1) = creac(i)%temperature
            end if
        end do

        allocate(mix)
        mix = Mixture(global_thermodb, input_reactants=input_reactants, ions=.true.)
        mptr = c_loc(mix)
        call log_info('BINDC: Created Mixture object at '//to_str(mptr))

    end function

    function cea_mixture_destroy(mptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(inout) :: mptr
        type(Mixture), pointer :: mix
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        deallocate(mix)
        call log_info('BINDC: Destroyed Mixture object at '//to_str(mptr))
        mptr = c_null_ptr
    end function

    function cea_mixture_get_num_species(mptr, num_species) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: mptr
        integer(c_int), intent(out) :: num_species
        type(Mixture), pointer :: mix
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        num_species = mix%num_species
    end function

    function cea_mixture_get_species_name(mptr, i_species, cspecies) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in) :: mptr
        integer(c_int), intent(in), value :: i_species
        type(c_ptr), intent(out) :: cspecies
        type(Mixture), pointer :: mix
        character(snl+1), pointer :: species
        integer :: len
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        if (i_species < 0 .or. i_species >= mix%num_species) then
            ierr = CEA_INVALID_INDEX
            return
        end if
        allocate(species)
        ! Convert names from fortran strings to C strings
        species = mix%species_names(i_species+1)
        len = len_trim(species)
        species(len+1:len+1) = c_null_char
        cspecies = c_loc(species)
    end function

    function cea_mixture_get_species_names(mptr, nspecies, cspecies) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in) :: mptr
        integer(c_int), intent(in), value :: nspecies
        type(c_ptr), intent(out) :: cspecies(nspecies)
        type(Mixture), pointer :: mix
        character(snl+1), pointer :: species(:)
        integer :: n
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        if (nspecies /= mix%num_species) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        allocate(species(nspecies))
        ! Convert names from fortran strings to C strings
        do n = 1,nspecies
            species(n) = trim(mix%species_names(n))//c_null_char
            cspecies(n) = c_loc(species(n))
        end do
    end function

    function cea_mixture_get_species_name_buf(mptr, i_species, cspecies, buf_len) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in) :: mptr
        integer(c_int), intent(in), value :: i_species
        character(c_char), intent(out) :: cspecies(*)
        integer(c_int), intent(in), value :: buf_len
        type(Mixture), pointer :: mix
        character(snl) :: name
        integer :: n, name_len, ncopy
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)

        if (i_species < 0 .or. i_species >= mix%num_species) then
            ierr = CEA_INVALID_INDEX
            return
        end if
        if (buf_len <= 0) then
            ierr = CEA_INVALID_SIZE
            return
        end if

        name = mix%species_names(i_species+1)
        name_len = len_trim(name)
        ncopy = min(name_len, buf_len-1)
        do n = 1, ncopy
            cspecies(n) = name(n:n)
        end do
        cspecies(ncopy+1) = c_null_char
        if (name_len + 1 > buf_len) ierr = CEA_INVALID_SIZE
    end function

    function cea_mixture_get_species_names_buf(mptr, nspecies, cspecies, stride) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in) :: mptr
        integer(c_int), intent(in), value :: nspecies
        character(c_char), intent(out) :: cspecies(*)
        integer(c_int), intent(in), value :: stride
        type(Mixture), pointer :: mix
        character(snl) :: name
        integer :: n, j, name_len, ncopy, base, nmax

        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)

        if (nspecies /= mix%num_species) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        if (stride <= 0) then
            ierr = CEA_INVALID_SIZE
            return
        end if

        nmax = mix%num_species

        do n = 1, nmax
            name = mix%species_names(n)
            name_len = len_trim(name)
            ncopy = min(name_len, stride-1)
            base = (n-1) * stride
            do j = 1, ncopy
                cspecies(base+j) = name(j:j)
            end do
            cspecies(base+ncopy+1) = c_null_char
            if (name_len + 1 > stride) ierr = CEA_INVALID_SIZE
        end do
    end function

    function cea_string_free(cspecies) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: cspecies
        character(snl+1), pointer :: species
        ierr = CEA_SUCCESS
        if (.not. c_associated(cspecies)) return
        call c_f_pointer(cspecies, species)
        if (associated(species)) deallocate(species)
    end function

    function cea_string_array_free(cspecies, nspecies) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in) :: cspecies(*)
        integer(c_int), intent(in), value :: nspecies
        character(snl+1), pointer :: species
        integer :: n
        ierr = CEA_SUCCESS
        if (nspecies <= 0) return
        do n = 1, nspecies
            if (.not. c_associated(cspecies(n))) cycle
            call c_f_pointer(cspecies(n), species)
            if (associated(species)) deallocate(species)
        end do
    end function

    function cea_mixture_moles_to_weights(mptr, len, moles, weights) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: mptr
        integer(c_int), intent(in), value :: len
        real(c_double), intent(in) :: moles(*)
        real(c_double), intent(out) :: weights(*)
        type(Mixture), pointer :: mix
        integer :: ns
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        ns = mix%num_species
        if (len < ns) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        weights(:ns) = mix%weights_from_moles(moles(:ns))
    end function

    function cea_mixture_weights_to_moles(mptr, len, weights, moles) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: mptr
        integer(c_int), intent(in), value :: len
        real(c_double), intent(in) :: weights(*)
        real(c_double), intent(out) :: moles(*)
        type(Mixture), pointer :: mix
        integer :: ns
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        ns = mix%num_species
        if (len < ns) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        moles(:ns) = mix%moles_from_weights(weights(:ns))
    end function

    function cea_mixture_per_mole_to_per_weight(mptr, len, per_mole, per_weight) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: mptr
        integer(c_int), intent(in), value :: len
        real(c_double), intent(in) :: per_mole(*)
        real(c_double), intent(out) :: per_weight(*)
        type(Mixture), pointer :: mix
        integer :: ns
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        ns = mix%num_species
        if (len < ns) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        per_weight(:ns) = mix%per_weight_from_per_mole(per_mole(:ns))
    end function

    function cea_mixture_per_weight_to_per_mole(mptr, len, per_weight, per_mole) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: mptr
        integer(c_int), intent(in), value :: len
        real(c_double), intent(in) :: per_weight(*)
        real(c_double), intent(out) :: per_mole(*)
        type(Mixture), pointer :: mix
        integer :: ns
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        ns = mix%num_species
        if (len < ns) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        per_mole(:ns) = mix%per_mole_from_per_weight(per_weight(:ns))
    end function

    function cea_mixture_chem_eq_ratio_to_of_ratio(mptr, len, oxidant_weights, fuel_weights, &
        chem_eq_ratio, of_ratio) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: mptr
        integer(c_int), intent(in), value :: len
        real(c_double), intent(in) :: oxidant_weights(*)
        real(c_double), intent(in) :: fuel_weights(*)
        real(c_double), intent(in), value :: chem_eq_ratio
        real(c_double), intent(out) :: of_ratio
        type(Mixture), pointer :: mix
        integer :: ns
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        ns = mix%num_species
        if (len < ns) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        of_ratio = mix%of_from_equivalence(oxidant_weights(:ns), fuel_weights(:ns), chem_eq_ratio)
    end function

    function cea_mixture_weight_eq_ratio_to_of_ratio(mptr, len, oxidant_weights, fuel_weights, &
        weight_eq_ratio, of_ratio) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: mptr
        integer(c_int), intent(in), value :: len
        real(c_double), intent(in) :: oxidant_weights(*)
        real(c_double), intent(in) :: fuel_weights(*)
        real(c_double), intent(in), value :: weight_eq_ratio
        real(c_double), intent(out) :: of_ratio
        type(Mixture), pointer :: mix
        integer :: ns
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        ns = mix%num_species
        if (len < ns) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        of_ratio = mix%of_from_phi(oxidant_weights(:ns), fuel_weights(:ns), weight_eq_ratio)
    end function

    function cea_mixture_of_ratio_to_weights(mptr, len, oxidant_weights, fuel_weights, &
        of_ratio, weights) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: mptr
        integer(c_int), intent(in), value :: len
        real(c_double), intent(in) :: oxidant_weights(*)
        real(c_double), intent(in) :: fuel_weights(*)
        real(c_double), intent(in), value :: of_ratio
        real(c_double), intent(out) :: weights(*)
        type(Mixture), pointer :: mix
        integer :: ns
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        ns = mix%num_species
        if (len < ns) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        weights(:ns) = mix%weights_from_of( &
            oxidant_weights(:ns), &
            fuel_weights(:ns), &
            of_ratio)
    end function

    function cea_mixture_calc_property(mptr, prop_type, len_weights, weights, &
        temperature, prop_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: mptr
        integer(c_int), intent(in), value :: prop_type
        integer(c_int), intent(in), value :: len_weights
        real(c_double), intent(in) :: weights(*)
        real(c_double), intent(in), value :: temperature
        real(c_double), intent(out) :: prop_value
        type(Mixture), pointer :: mix
        integer :: ns
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        ns = mix%num_species
        if (len_weights < ns) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        select case(prop_type)
            case (CEA_ENTHALPY)
                prop_value = mix%calc_enthalpy(weights(:ns), temperature)
            case (CEA_ENERGY)
                prop_value = mix%calc_energy(weights(:ns), temperature)
            case (CEA_FROZEN_CP)
                prop_value = mix%calc_frozen_cp(weights(:ns), temperature)
            case (CEA_FROZEN_CV)
                prop_value = mix%calc_frozen_cv(weights(:ns), temperature)
            case default
                prop_value = empty_dp
                ierr = CEA_INVALID_PROPERTY_TYPE
        end select
    end function

    function cea_mixture_calc_property_multitemp(mptr, prop_type, len_weights, &
        weights, len_temperatures, temperatures, prop_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: mptr
        integer(c_int), intent(in), value :: prop_type
        integer(c_int), intent(in), value :: len_weights
        real(c_double), intent(in) :: weights(*)
        integer(c_int), intent(in), value :: len_temperatures
        real(c_double), intent(in) :: temperatures(*)
        real(c_double), intent(out) :: prop_value
        type(Mixture), pointer :: mix
        integer :: ns
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        ns = mix%num_species
        if (len_weights < ns .or. len_temperatures < ns) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        select case(prop_type)
            case (CEA_ENTHALPY)
                prop_value = mix%calc_enthalpy(weights(:ns), temperatures(:ns))
            case (CEA_ENERGY)
                prop_value = mix%calc_energy(weights(:ns), temperatures(:ns))
            case (CEA_FROZEN_CP)
                prop_value = mix%calc_frozen_cp(weights(:ns), temperatures(:ns))
            case (CEA_FROZEN_CV)
                prop_value = mix%calc_frozen_cv(weights(:ns), temperatures(:ns))
            case default
                prop_value = empty_dp
                ierr = CEA_INVALID_PROPERTY_TYPE
        end select
    end function

    function cea_mixture_calc_property_tp(mptr, prop_type, len_weights, weights, &
        temperature, pressure, prop_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: mptr
        integer(c_int), intent(in), value :: prop_type
        integer(c_int), intent(in), value :: len_weights
        real(c_double), intent(in) :: weights(*)
        real(c_double), intent(in), value :: temperature
        real(c_double), intent(in), value :: pressure
        real(c_double), intent(out) :: prop_value
        type(Mixture), pointer :: mix
        integer :: ns
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        ns = mix%num_species
        if (len_weights < ns) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        select case(prop_type)
            case (CEA_VOLUME)
                prop_value = 1.d-2*mix%calc_pressure(weights(:ns)/sum(weights(:ns)), temperature)/pressure
            case (CEA_DENSITY)
                prop_value = 1.d2*pressure/mix%calc_pressure(weights(:ns)/sum(weights(:ns)), temperature)
            case (CEA_ENTROPY)
                prop_value = mix%calc_entropy(weights(:ns), temperature, pressure)
            case (CEA_GIBBS_ENERGY)
                prop_value = mix%calc_gibbs_energy(weights(:ns), temperature, pressure)
            case default
                prop_value = empty_dp
                ierr = CEA_INVALID_PROPERTY_TYPE
        end select
    end function

    function cea_mixture_calc_property_tp_multitemp(mptr, prop_type, len_weights, weights, &
        len_temperatures, temperatures, pressure, prop_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: mptr
        integer(c_int), intent(in), value :: prop_type
        integer(c_int), intent(in), value :: len_weights
        real(c_double), intent(in) :: weights(*)
        integer(c_int), intent(in), value :: len_temperatures
        real(c_double), intent(in) :: temperatures(*)
        real(c_double), intent(in), value :: pressure
        real(c_double), intent(out) :: prop_value
        type(Mixture), pointer :: mix
        integer :: ns
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, mix)
        ns = mix%num_species
        if (len_weights < ns .or. len_temperatures < ns) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        select case(prop_type)
            case (CEA_VOLUME)
                prop_value = mix%calc_pressure(weights(:ns)/sum(weights(:ns)), temperatures(:ns))/pressure
            case (CEA_DENSITY)
                prop_value = pressure/mix%calc_pressure(weights(:ns)/sum(weights(:ns)), temperatures(:ns))
            case (CEA_ENTROPY)
                block
                    real(wp) :: pressures(ns)
                    pressures = pressure
                    prop_value = mix%calc_entropy(weights(:ns), temperatures(:ns), pressures)
                end block
            case (CEA_GIBBS_ENERGY)
                block
                    real(wp) :: pressures(ns)
                    pressures = pressure
                    prop_value = mix%calc_gibbs_energy(weights(:ns), temperatures(:ns), pressures)
                end block
            case default
                prop_value = empty_dp
                ierr = CEA_INVALID_PROPERTY_TYPE
        end select
    end function


    !-----------------------------------------------------------------
    ! Equilibrium Solver
    !-----------------------------------------------------------------
    function cea_eqsolver_create(sptr, mptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: sptr
        type(c_ptr), intent(in), value :: mptr
        type(EqSolver), pointer :: solver
        type(Mixture), pointer :: products
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, products)
        allocate(solver)
        solver = EqSolver(products)
        sptr = c_loc(solver)
        call log_info('BINDC: Created EqSolver from product mixture at '//to_str(sptr))
    end function

    function cea_eqsolver_create_with_reactants(sptr, mptr_prod, mptr_reac) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: sptr
        type(c_ptr), intent(in), value :: mptr_prod
        type(c_ptr), intent(in), value :: mptr_reac
        type(EqSolver), pointer :: solver
        type(Mixture), pointer :: products
        type(Mixture), pointer :: reactants
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr_prod, products)
        call c_f_pointer(mptr_reac, reactants)
        allocate(solver)
        solver = EqSolver(products, reactants)
        sptr = c_loc(solver)
        call log_info('BINDC: Created EqSolver from product/reactant mixtures at '//to_str(sptr))
    end function

    function cea_eqsolver_create_with_options(sptr, mptr_prod, opts) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: sptr
        type(c_ptr), intent(in), value :: mptr_prod
        type(cea_solver_opts), intent(in), value :: opts

        ! Locals
        type(EqSolver), pointer :: solver
        type(Mixture), pointer :: products
        type(Mixture), pointer :: reactants
        type(c_ptr), pointer :: cinsert(:)
        character(snl), allocatable :: insert(:)
        character(:), allocatable :: name
        integer :: n
        logical :: ions, transport, use_trace

        ierr = CEA_SUCCESS
        call c_f_pointer(mptr_prod, products)

        if (opts%ninsert < 0) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        if (opts%ninsert > 0 .and. .not. c_associated(opts%insert)) then
            ierr = CEA_INVALID_SIZE
            return
        end if

        allocate(cinsert(opts%ninsert), insert(opts%ninsert))

        ! Handle optional reactants mixture
        if (c_associated(opts%reactants)) then
            call c_f_pointer(opts%reactants, reactants)
        else
            nullify(reactants)
        end if

        ! Handle trace option
        if (opts%trace > 0.0d0) then
            use_trace = .true.
        else
            use_trace = .false.
        end if

        ! Convert ions from C bool to Fortran logical
        ions = logical(opts%ions)
        transport = logical(opts%transport)

        ! Convert insert species names from C strings to Fortran strings
        if (opts%ninsert > 0 .and. c_associated(opts%insert)) then
            call c_f_pointer(opts%insert, cinsert, [opts%ninsert])

            do n = 1, opts%ninsert
                call c_copy(cinsert(n), name)
                if (len(name) > snl) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                insert(n) = name
            end do
        end if

        allocate(solver)
        if (associated(reactants)) then
            if (transport) then
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = EqSolver(products=products, reactants=reactants, trace=opts%trace, &
                                          ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = EqSolver(products=products, reactants=reactants, trace=opts%trace, &
                                          ions=ions, all_transport=global_transdb)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = EqSolver(products=products, reactants=reactants, &
                                          ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = EqSolver(products=products, reactants=reactants, &
                                          ions=ions, all_transport=global_transdb)
                    end if
                end if
            else
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = EqSolver(products=products, reactants=reactants, trace=opts%trace, &
                                          ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = EqSolver(products=products, reactants=reactants, trace=opts%trace, &
                                          ions=ions)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = EqSolver(products=products, reactants=reactants, &
                                          ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = EqSolver(products=products, reactants=reactants, &
                                          ions=ions)
                    end if
                end if
            end if
        else
            if (transport) then
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = EqSolver(products=products, trace=opts%trace, &
                                          ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = EqSolver(products=products, trace=opts%trace, &
                                          ions=ions, all_transport=global_transdb)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = EqSolver(products=products, &
                                          ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = EqSolver(products=products, &
                                          ions=ions, all_transport=global_transdb)
                    end if
                end if
            else
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = EqSolver(products=products, trace=opts%trace, &
                                          ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = EqSolver(products=products, trace=opts%trace, &
                                          ions=ions)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = EqSolver(products=products, &
                                          ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = EqSolver(products=products, &
                                          ions=ions)
                    end if
                end if
            end if
        end if

        sptr = c_loc(solver)
        call log_info('BINDC: Created EqSolver from product mixture with options at '//to_str(sptr))
    end function

    function cea_eqsolver_destroy(sptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(inout) :: sptr
        type(EqSolver), pointer :: solver
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr, solver)
        deallocate(solver)
        call log_info('BINDC: Destroyed EqSolver object at '//to_str(sptr))
        sptr = c_null_ptr
    end function

    function cea_eqsolver_solve(sptr, eq_type, state1, state2, amounts, slptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr),    intent(in), value :: sptr
        integer(kind=kind(CEA_TP)), intent(in), value :: eq_type
        real(c_double), intent(in), value :: state1
        real(c_double), intent(in), value :: state2
        real(c_double), intent(in) :: amounts(*)
        type(c_ptr),    intent(in), value :: slptr
        type(EqSolver), pointer :: solver
        type(EqSolution), pointer :: solution
        character(2) :: type
        integer :: nr
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr,  solver)
        call c_f_pointer(slptr, solution)
        select case(eq_type)
            case (CEA_TP); type = 'tp'
            case (CEA_HP); type = 'hp'
            case (CEA_SP); type = 'sp'
            case (CEA_TV); type = 'tv'
            case (CEA_UV); type = 'uv'
            case (CEA_SV); type = 'sv'
            case default
                ierr = CEA_INVALID_EQUILIBRIUM_TYPE
                return
        end select
        nr = solver%num_reactants
        call solver%solve(solution, type, state1, state2, amounts(:nr))
        if (.not. solution%converged) ierr = CEA_NOT_CONVERGED
    end function

    function cea_eqsolver_solve_with_partials(sptr, eq_type, state1, state2, amounts, slptr, pptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr),    intent(in), value :: sptr
        integer(kind=kind(CEA_TP)), intent(in), value :: eq_type
        real(c_double), intent(in), value :: state1
        real(c_double), intent(in), value :: state2
        real(c_double), intent(in) :: amounts(*)
        type(c_ptr),    intent(in), value :: slptr
        type(c_ptr),    intent(in), value :: pptr
        type(EqSolver), pointer :: solver
        type(EqSolution), pointer :: solution
        type(EqPartials), pointer :: partials
        character(2) :: type
        integer :: nr
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr,  solver)
        call c_f_pointer(slptr, solution)
        call c_f_pointer(pptr, partials)
        select case(eq_type)
            case (CEA_TP); type = 'tp'
            case (CEA_HP); type = 'hp'
            case (CEA_SP); type = 'sp'
            case (CEA_TV); type = 'tv'
            case (CEA_UV); type = 'uv'
            case (CEA_SV); type = 'sv'
            case default
                ierr = CEA_INVALID_EQUILIBRIUM_TYPE
                return
        end select
        nr = solver%num_reactants
        call solver%solve(solution, type, state1, state2, amounts(:nr), partials)
        if (.not. solution%converged) ierr = CEA_NOT_CONVERGED
    end function

    function cea_eqsolver_get_size(sptr, eq_variable, eq_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr),    intent(in), value :: sptr
        integer(c_int), intent(in), value :: eq_variable
        integer(c_int), intent(out) :: eq_value
        type(EqSolver), pointer :: solver
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr,  solver)
        select case(eq_variable)
            case (CEA_NUM_REACTANTS); eq_value = solver%num_reactants
            case (CEA_NUM_PRODUCTS);  eq_value = solver%num_products
            case (CEA_NUM_GAS);       eq_value = solver%num_gas
            case (CEA_NUM_CONDENSED); eq_value = solver%num_condensed
            case (CEA_NUM_ELEMENTS);  eq_value = solver%num_elements
            case (CEA_MAX_EQUATIONS); eq_value = solver%max_equations
            case default
                eq_value = -1
                ierr = CEA_INVALID_EQUILIBRIUM_SIZE_TYPE
                return
        end select
    end function

    !-----------------------------------------------------------------
    ! Rocket Solver
    !-----------------------------------------------------------------
    function cea_rocket_solver_create(sptr, mptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: sptr
        type(c_ptr), intent(in), value :: mptr
        type(RocketSolver), pointer :: solver
        type(Mixture), pointer :: products
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, products)
        allocate(solver)
        solver = RocketSolver(products)
        sptr = c_loc(solver)
        call log_info('BINDC: Created RocketSolver from product mixture at '//to_str(sptr))
    end function

    function cea_rocket_solver_create_with_reactants(sptr, mptr_prod, mptr_reac) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: sptr
        type(c_ptr), intent(in), value :: mptr_prod
        type(c_ptr), intent(in), value :: mptr_reac
        type(RocketSolver), pointer :: solver
        type(Mixture), pointer :: products
        type(Mixture), pointer :: reactants
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr_prod, products)
        call c_f_pointer(mptr_reac, reactants)
        allocate(solver)
        solver = RocketSolver(products, reactants)
        sptr = c_loc(solver)
        call log_info('BINDC: Created RocketSolver from product/reactant mixtures at '//to_str(sptr))
    end function

    function cea_rocket_solver_create_with_options(sptr, mptr_prod, opts) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: sptr
        type(c_ptr), intent(in), value :: mptr_prod
        type(cea_solver_opts), intent(in), value :: opts

        ! Locals
        type(RocketSolver), pointer :: solver
        type(Mixture), pointer :: products
        type(Mixture), pointer :: reactants
        type(c_ptr), pointer :: cinsert(:)
        character(snl), allocatable :: insert(:)
        character(:), allocatable :: name
        integer :: n
        logical :: ions, transport, use_trace

        ierr = CEA_SUCCESS
        call c_f_pointer(mptr_prod, products)

        if (opts%ninsert < 0) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        if (opts%ninsert > 0 .and. .not. c_associated(opts%insert)) then
            ierr = CEA_INVALID_SIZE
            return
        end if

        allocate(cinsert(opts%ninsert), insert(opts%ninsert))

        ! Handle optional reactants mixture
        if (c_associated(opts%reactants)) then
            call c_f_pointer(opts%reactants, reactants)
        else
            nullify(reactants)
        end if

        ! Handle trace option
        if (opts%trace > 0.0d0) then
            use_trace = .true.
        else
            use_trace = .false.
        end if

        ! Convert ions from C bool to Fortran logical
        ions = logical(opts%ions)
        transport = logical(opts%transport)

        ! Convert insert species names from C strings to Fortran strings
        if (opts%ninsert > 0 .and. c_associated(opts%insert)) then
            call c_f_pointer(opts%insert, cinsert, [opts%ninsert])

            do n = 1, opts%ninsert
                call c_copy(cinsert(n), name)
                if (len(name) > snl) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                insert(n) = name
            end do
        end if

        allocate(solver)
        if (associated(reactants)) then
            if (transport) then
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = RocketSolver(products=products, reactants=reactants, trace=opts%trace, &
                                              ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = RocketSolver(products=products, reactants=reactants, trace=opts%trace, &
                                              ions=ions, all_transport=global_transdb)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = RocketSolver(products=products, reactants=reactants, &
                                              ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = RocketSolver(products=products, reactants=reactants, &
                                              ions=ions, all_transport=global_transdb)
                    end if
                end if
            else
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = RocketSolver(products=products, reactants=reactants, trace=opts%trace, &
                                              ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = RocketSolver(products=products, reactants=reactants, trace=opts%trace, &
                                              ions=ions)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = RocketSolver(products=products, reactants=reactants, &
                                              ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = RocketSolver(products=products, reactants=reactants, &
                                              ions=ions)
                    end if
                end if
            end if
        else
            if (transport) then
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = RocketSolver(products=products, trace=opts%trace, &
                                              ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = RocketSolver(products=products, trace=opts%trace, &
                                              ions=ions, all_transport=global_transdb)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = RocketSolver(products=products, &
                                              ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = RocketSolver(products=products, &
                                              ions=ions, all_transport=global_transdb)
                    end if
                end if
            else
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = RocketSolver(products=products, trace=opts%trace, &
                                              ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = RocketSolver(products=products, trace=opts%trace, &
                                              ions=ions)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = RocketSolver(products=products, &
                                              ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = RocketSolver(products=products, &
                                              ions=ions)
                    end if
                end if
            end if
        end if
        sptr = c_loc(solver)
        call log_info('BINDC: Created RocketSolver with options at '//to_str(sptr))
    end function

    function cea_rocket_solver_destroy(sptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(inout) :: sptr
        type(RocketSolver), pointer :: solver
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr, solver)
        deallocate(solver)
        call log_info('BINDC: Destroyed RocketSolver object at '//to_str(sptr))
        sptr = c_null_ptr
    end function

    function cea_rocket_solver_get_size(sptr, eq_variable, eq_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr),    intent(in), value :: sptr
        integer(c_int), intent(in), value :: eq_variable
        integer(c_int), intent(out) :: eq_value
        type(RocketSolver), pointer :: solver
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr,  solver)
        select case(eq_variable)
            case (CEA_NUM_REACTANTS); eq_value = solver%eq_solver%num_reactants
            case (CEA_NUM_PRODUCTS);  eq_value = solver%eq_solver%num_products
            case (CEA_NUM_GAS);       eq_value = solver%eq_solver%num_gas
            case (CEA_NUM_CONDENSED); eq_value = solver%eq_solver%num_condensed
            case (CEA_NUM_ELEMENTS);  eq_value = solver%eq_solver%num_elements
            case (CEA_MAX_EQUATIONS); eq_value = solver%eq_solver%max_equations
            case default
                eq_value = -1
                ierr = CEA_INVALID_EQUILIBRIUM_SIZE_TYPE
                return
        end select
    end function

    function cea_rocket_solver_solve_iac(sptr, slptr, weights, pc, pi_p, n_pi_p, subar, nsubar, &
        supar, nsupar, n_frz, hc_or_tc, use_hc, tc_est, use_tc_est) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr),     intent(in), value :: sptr
        type(c_ptr),     intent(in), value :: slptr
        real(c_double),  intent(in)        :: weights(*)
        real(c_double),  intent(in), value :: pc
        real(c_double),  intent(in)        :: pi_p(*)
        integer(c_int),  intent(in), value :: n_pi_p
        real(c_double),  intent(in)        :: subar(*)
        integer(c_int),  intent(in), value :: nsubar
        real(c_double),  intent(in)        :: supar(*)
        integer(c_int),  intent(in), value :: nsupar
        integer(c_int),  intent(in), value :: n_frz  ! Set to 0 to ignore
        real(c_double),  intent(in), value :: hc_or_tc  ! Value of hc or tc
        logical(c_bool), intent(in), value :: use_hc  ! if true, use hc; if false, use tc
        real(c_double),  intent(in), value :: tc_est
        logical(c_bool), intent(in), value :: use_tc_est  ! Initial guess for chamber temperature; only valid if hc is used

        type(RocketSolver), pointer :: solver
        type(RocketSolution), pointer :: solution
        integer :: nr

        ierr = CEA_SUCCESS

        call c_f_pointer(sptr,  solver)
        call c_f_pointer(slptr, solution)

        nr = solver%eq_solver%num_reactants

        ! Call the Rocket Solver
        if (nsubar > 0) then
            if (nsupar > 0) then
                if (use_tc_est .eqv. .true.) then
                    if (use_hc .eqv. .true.) then
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                              supar=supar(:nsupar), n_frz=n_frz, tc_est=tc_est, hc=hc_or_tc)
                    else
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                              supar=supar(:nsupar), n_frz=n_frz, tc_est=tc_est, tc=hc_or_tc)
                    end if
                else
                    if (use_hc .eqv. .true.) then
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                              supar=supar(:nsupar), n_frz=n_frz, hc=hc_or_tc)
                    else
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                              supar=supar(:nsupar), n_frz=n_frz, tc=hc_or_tc)
                    end if
                end if
            else  ! nsupar == 0
                if (use_tc_est .eqv. .true.) then
                    if (use_hc .eqv. .true.) then
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                              n_frz=n_frz, tc_est=tc_est, hc=hc_or_tc)
                    else
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                              n_frz=n_frz, tc_est=tc_est, tc=hc_or_tc)
                    end if
                else
                    if (use_hc .eqv. .true.) then
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                              n_frz=n_frz, hc=hc_or_tc)
                    else
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                              n_frz=n_frz, tc=hc_or_tc)
                    end if
                end if
            end if
        else  ! nsubar == 0
            if (nsupar > 0) then
                if (use_tc_est .eqv. .true.) then
                    if (use_hc .eqv. .true.) then
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                              supar=supar(:nsupar), n_frz=n_frz, tc_est=tc_est, hc=hc_or_tc)
                    else
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                              supar=supar(:nsupar), n_frz=n_frz, tc_est=tc_est, tc=hc_or_tc)
                    end if
                else
                    if (use_hc .eqv. .true.) then
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                              supar=supar(:nsupar), n_frz=n_frz, hc=hc_or_tc)
                    else
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                              supar=supar(:nsupar), n_frz=n_frz, tc=hc_or_tc)
                    end if
                end if
            else  ! nsupar == 0
                if (use_tc_est .eqv. .true.) then
                    if (use_hc .eqv. .true.) then
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                              n_frz=n_frz, tc_est=tc_est, hc=hc_or_tc)
                    else
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                              n_frz=n_frz, tc_est=tc_est, tc=hc_or_tc)
                    end if
                else
                    if (use_hc .eqv. .true.) then
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                              n_frz=n_frz, hc=hc_or_tc)
                    else
                        call solver%solve_iac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                              n_frz=n_frz, tc=hc_or_tc)
                    end if
                end if
            end if
        end if
        if (.not. solution%converged) ierr = CEA_NOT_CONVERGED
    end function

    function cea_rocket_solver_solve_fac(sptr, slptr, weights, pc, pi_p, n_pi_p, subar, nsubar, &
        supar, nsupar, n_frz, hc_or_tc, use_hc, mdot_or_acat, use_mdot, tc_est, use_tc_est) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr),     intent(in), value :: sptr
        type(c_ptr),     intent(in), value :: slptr
        real(c_double),  intent(in)        :: weights(*)
        real(c_double),  intent(in), value :: pc
        real(c_double),  intent(in)        :: pi_p(*)
        integer(c_int),  intent(in), value :: n_pi_p
        real(c_double),  intent(in)        :: subar(*)
        integer(c_int),  intent(in), value :: nsubar
        real(c_double),  intent(in)        :: supar(*)
        integer(c_int),  intent(in), value :: nsupar
        integer(c_int),  intent(in), value :: n_frz       ! Set to 0 to ignore
        real(c_double),  intent(in), value :: hc_or_tc    ! Value of hc or tc
        logical(c_bool), intent(in), value :: use_hc      ! if true, use hc; if false, use tc
        real(c_double),  intent(in), value :: mdot_or_acat
        logical(c_bool), intent(in), value :: use_mdot    ! if true, use mdot; if false, use ac_at
        real(c_double),  intent(in), value :: tc_est
        logical(c_bool), intent(in), value :: use_tc_est  ! Initial guess for chamber temperature; only valid if hc is used

        type(RocketSolver), pointer :: solver
        type(RocketSolution), pointer :: solution
        integer :: nr

        ierr = CEA_SUCCESS

        call c_f_pointer(sptr,  solver)
        call c_f_pointer(slptr, solution)

        nr = solver%eq_solver%num_reactants

        ! Call the Rocket Solver
        if (nsubar > 0) then
            if (nsupar > 0) then
                if (use_tc_est .eqv. .true.) then
                    if (use_hc .eqv. .true.) then
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  supar=supar(:nsupar), n_frz=n_frz, tc_est=tc_est, hc=hc_or_tc, mdot=mdot_or_acat)
                        else  ! use ac_at
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  supar=supar(:nsupar), n_frz=n_frz, tc_est=tc_est, hc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    else  ! use tc
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  supar=supar(:nsupar), n_frz=n_frz, tc_est=tc_est, tc=hc_or_tc, mdot=mdot_or_acat)
                        else  ! use ac_at
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  supar=supar(:nsupar), n_frz=n_frz, tc_est=tc_est, tc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    end if
                else  ! use_tc_est == .false.
                    if (use_hc .eqv. .true.) then
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  supar=supar(:nsupar), n_frz=n_frz, hc=hc_or_tc, mdot=mdot_or_acat)
                        else  ! use ac_at
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  supar=supar(:nsupar), n_frz=n_frz, hc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    else  ! use tc
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  supar=supar(:nsupar), n_frz=n_frz, tc=hc_or_tc, mdot=mdot_or_acat)
                        else  ! use ac_at
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  supar=supar(:nsupar), n_frz=n_frz, tc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    end if
                end if
            else  ! nsupar == 0
                if (use_tc_est .eqv. .true.) then
                    if (use_hc .eqv. .true.) then
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  n_frz=n_frz, tc_est=tc_est, hc=hc_or_tc, mdot=mdot_or_acat)
                        else  ! use ac_at
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  n_frz=n_frz, tc_est=tc_est, hc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    else  ! use tc
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  n_frz=n_frz, tc_est=tc_est, tc=hc_or_tc, mdot=mdot_or_acat)
                        else  ! use ac_at
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  n_frz=n_frz, tc_est=tc_est, tc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    end if
                else  ! use_tc_est == .false.
                    if (use_hc .eqv. .true.) then
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  n_frz=n_frz, hc=hc_or_tc, mdot=mdot_or_acat)
                        else  ! use ac_at
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  n_frz=n_frz, hc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    else  ! use tc
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  n_frz=n_frz, tc=hc_or_tc, mdot=mdot_or_acat)
                        else  ! use ac_at
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  n_frz=n_frz, tc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    end if
                end if
            end if
        else  ! nsubar == 0
            if (nsupar > 0) then
                if (use_tc_est .eqv. .true.) then
                    if (use_hc .eqv. .true.) then
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  supar=supar(:nsupar), n_frz=n_frz, tc_est=tc_est, hc=hc_or_tc, mdot=mdot_or_acat)
                        else  ! use ac_at
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  supar=supar(:nsupar), n_frz=n_frz, tc_est=tc_est, hc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    else  ! use tc
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  supar=supar(:nsupar), n_frz=n_frz, tc_est=tc_est, tc=hc_or_tc, mdot=mdot_or_acat)
                        else  ! use ac_at
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  supar=supar(:nsupar), n_frz=n_frz, tc_est=tc_est, tc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    end if
                else  ! use_tc_est == .false.
                    if (use_hc .eqv. .true.) then
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  supar=supar(:nsupar), n_frz=n_frz, hc=hc_or_tc, mdot=mdot_or_acat)
                        else  ! use ac_at
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), subar=subar(:nsubar), &
                                                  supar=supar(:nsupar), n_frz=n_frz, hc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    else  ! use tc
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  supar=supar(:nsupar), n_frz=n_frz, tc=hc_or_tc, mdot=mdot_or_acat)
                        else  ! use ac_at
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  supar=supar(:nsupar), n_frz=n_frz, tc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    end if
                end if
            else  ! nsupar == 0
                if (use_tc_est .eqv. .true.) then
                    if (use_hc .eqv. .true.) then
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  n_frz=n_frz, tc_est=tc_est, hc=hc_or_tc, mdot=mdot_or_acat)
                        else
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  n_frz=n_frz, tc_est=tc_est, hc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    else  ! use tc
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  n_frz=n_frz, tc_est=tc_est, tc=hc_or_tc, mdot=mdot_or_acat)
                        else
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  n_frz=n_frz, tc_est=tc_est, tc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    end if
                else  ! use_tc_est == .false.
                    if (use_hc .eqv. .true.) then
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  n_frz=n_frz, hc=hc_or_tc, mdot=mdot_or_acat)
                        else
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  n_frz=n_frz, hc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    else  ! use tc
                        if (use_mdot .eqv. .true.) then
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  n_frz=n_frz, tc=hc_or_tc, mdot=mdot_or_acat)
                        else
                            call solver%solve_fac(solution, weights(:nr), pc, pi_p(:n_pi_p), &
                                                  n_frz=n_frz, tc=hc_or_tc, ac_at=mdot_or_acat)
                        end if
                    end if
                end if
            end if
        end if
        if (.not. solution%converged) ierr = CEA_NOT_CONVERGED
    end function

    !-----------------------------------------------------------------
    ! Shock Solver
    !-----------------------------------------------------------------
    function cea_shock_solver_create(sptr, mptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: sptr
        type(c_ptr), intent(in), value :: mptr
        type(ShockSolver), pointer :: solver
        type(Mixture), pointer :: products
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, products)
        allocate(solver)
        solver = ShockSolver(products)
        sptr = c_loc(solver)
        call log_info('BINDC: Created ShockSolver from product mixture at '//to_str(sptr))
    end function

    function cea_shock_solver_create_with_reactants(sptr, mptr_prod, mptr_reac) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: sptr
        type(c_ptr), intent(in), value :: mptr_prod
        type(c_ptr), intent(in), value :: mptr_reac
        type(ShockSolver), pointer :: solver
        type(Mixture), pointer :: products
        type(Mixture), pointer :: reactants
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr_prod, products)
        call c_f_pointer(mptr_reac, reactants)
        allocate(solver)
        solver = ShockSolver(products, reactants)
        sptr = c_loc(solver)
        call log_info('BINDC: Created ShockSolver from product/reactant mixtures at '//to_str(sptr))
    end function

    function cea_shock_solver_create_with_options(sptr, mptr_prod, opts) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: sptr
        type(c_ptr), intent(in), value :: mptr_prod
        type(cea_solver_opts), intent(in), value :: opts

        ! Locals
        type(ShockSolver), pointer :: solver
        type(Mixture), pointer :: products
        type(Mixture), pointer :: reactants
        type(c_ptr), pointer :: cinsert(:)
        character(snl), allocatable :: insert(:)
        character(:), allocatable :: name
        integer :: n
        logical :: ions, transport, use_trace

        ierr = CEA_SUCCESS
        call c_f_pointer(mptr_prod, products)

        if (opts%ninsert < 0) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        if (opts%ninsert > 0 .and. .not. c_associated(opts%insert)) then
            ierr = CEA_INVALID_SIZE
            return
        end if

        allocate(cinsert(opts%ninsert), insert(opts%ninsert))

        ! Handle optional reactants mixture
        if (c_associated(opts%reactants)) then
            call c_f_pointer(opts%reactants, reactants)
        else
            nullify(reactants)
        end if

        ! Handle trace option
        if (opts%trace > 0.0d0) then
            use_trace = .true.
        else
            use_trace = .false.
        end if

        ! Convert ions from C bool to Fortran logical
        ions = logical(opts%ions)
        transport = logical(opts%transport)

        ! Convert insert species names from C strings to Fortran strings
        if (opts%ninsert > 0 .and. c_associated(opts%insert)) then
            call c_f_pointer(opts%insert, cinsert, [opts%ninsert])

            do n = 1, opts%ninsert
                call c_copy(cinsert(n), name)
                if (len(name) > snl) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                insert(n) = name
            end do
        end if

        allocate(solver)
        if (associated(reactants)) then
            if (transport) then
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = ShockSolver(products=products, reactants=reactants, trace=opts%trace, &
                                             ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = ShockSolver(products=products, reactants=reactants, trace=opts%trace, &
                                             ions=ions, all_transport=global_transdb)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = ShockSolver(products=products, reactants=reactants, &
                                             ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = ShockSolver(products=products, reactants=reactants, &
                                             ions=ions, all_transport=global_transdb)
                    end if
                end if
            else
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = ShockSolver(products=products, reactants=reactants, trace=opts%trace, &
                                             ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = ShockSolver(products=products, reactants=reactants, trace=opts%trace, &
                                             ions=ions)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = ShockSolver(products=products, reactants=reactants, &
                                             ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = ShockSolver(products=products, reactants=reactants, &
                                             ions=ions)
                    end if
                end if
            end if
        else
            if (transport) then
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = ShockSolver(products=products, trace=opts%trace, &
                                             ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = ShockSolver(products=products, trace=opts%trace, &
                                             ions=ions, all_transport=global_transdb)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = ShockSolver(products=products, &
                                             ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = ShockSolver(products=products, &
                                             ions=ions, all_transport=global_transdb)
                    end if
                end if
            else
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = ShockSolver(products=products, trace=opts%trace, &
                                             ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = ShockSolver(products=products, trace=opts%trace, &
                                             ions=ions)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = ShockSolver(products=products, &
                                             ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = ShockSolver(products=products, &
                                             ions=ions)
                    end if
                end if
            end if
        end if
        sptr = c_loc(solver)
        call log_info('BINDC: Created ShockSolver with options at '//to_str(sptr))
    end function

    function cea_shock_solver_destroy(sptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(inout) :: sptr
        type(ShockSolver), pointer :: solver
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr, solver)
        deallocate(solver)
        call log_info('BINDC: Destroyed ShockSolver object at '//to_str(sptr))
        sptr = c_null_ptr
    end function

    function cea_shock_solver_get_size(sptr, eq_variable, eq_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr),    intent(in), value :: sptr
        integer(c_int), intent(in), value :: eq_variable
        integer(c_int), intent(out) :: eq_value
        type(ShockSolver), pointer :: solver
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr,  solver)
        select case(eq_variable)
            case (CEA_NUM_REACTANTS); eq_value = solver%eq_solver%num_reactants
            case (CEA_NUM_PRODUCTS);  eq_value = solver%eq_solver%num_products
            case (CEA_NUM_GAS);       eq_value = solver%eq_solver%num_gas
            case (CEA_NUM_CONDENSED); eq_value = solver%eq_solver%num_condensed
            case (CEA_NUM_ELEMENTS);  eq_value = solver%eq_solver%num_elements
            case (CEA_MAX_EQUATIONS); eq_value = solver%eq_solver%max_equations
            case default
                eq_value = -1
                ierr = CEA_INVALID_EQUILIBRIUM_SIZE_TYPE
                return
        end select
    end function

    function cea_shock_solver_solve(sptr, slptr, weights, T0, p0, mach1_or_u1, use_mach, refl, incd_froz, refl_froz) &
        result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr),     intent(in), value :: sptr
        type(c_ptr),     intent(in), value :: slptr
        real(c_double),  intent(in)        :: weights(*)
        real(c_double),  intent(in), value :: T0
        real(c_double),  intent(in), value :: p0
        real(c_double),  intent(in), value :: mach1_or_u1
        logical(c_bool), intent(in), value :: use_mach
        logical(c_bool), intent(in), value :: refl
        logical(c_bool), intent(in), value :: incd_froz
        logical(c_bool), intent(in), value :: refl_froz

        type(ShockSolver), pointer :: solver
        type(ShockSolution), pointer :: solution
        logical :: reflected, incident_frozen, reflected_frozen
        integer :: nr

        ierr = CEA_SUCCESS

        call c_f_pointer(sptr,  solver)
        call c_f_pointer(slptr, solution)

        nr = solver%eq_solver%num_reactants

        reflected = .false.
        if (refl .eqv. .true.) reflected = .true.
        incident_frozen = .false.
        if (incd_froz .eqv. .true.) incident_frozen = .true.
        reflected_frozen = .false.
        if (refl_froz .eqv. .true.) reflected_frozen = .true.

        if (use_mach .eqv. .true.) then
            solution = solver%solve(weights(:nr), T0, p0, mach1=mach1_or_u1, &
                reflected=reflected, incident_frozen=incident_frozen, reflected_frozen=reflected_frozen)
        else
            solution = solver%solve(weights(:nr), T0, p0, u1=mach1_or_u1, &
                reflected=reflected, incident_frozen=incident_frozen, reflected_frozen=reflected_frozen)
        end if
        if (.not. solution%converged) ierr = CEA_NOT_CONVERGED
    end function


    !-----------------------------------------------------------------
    ! Detonation Solver
    !-----------------------------------------------------------------
    function cea_detonation_solver_create(sptr, mptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: sptr
        type(c_ptr), intent(in), value :: mptr
        type(DetonSolver), pointer :: solver
        type(Mixture), pointer :: products
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr, products)
        allocate(solver)
        solver = DetonSolver(products)
        sptr = c_loc(solver)
        call log_info('BINDC: Created DetonSolver from product mixture at '//to_str(sptr))
    end function

    function cea_detonation_solver_create_with_reactants(sptr, mptr_prod, mptr_reac) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: sptr
        type(c_ptr), intent(in), value :: mptr_prod
        type(c_ptr), intent(in), value :: mptr_reac
        type(DetonSolver), pointer :: solver
        type(Mixture), pointer :: products
        type(Mixture), pointer :: reactants
        ierr = CEA_SUCCESS
        call c_f_pointer(mptr_prod, products)
        call c_f_pointer(mptr_reac, reactants)
        allocate(solver)
        solver = DetonSolver(products, reactants)
        sptr = c_loc(solver)
        call log_info('BINDC: Created DetonSolver from product/reactant mixtures at '//to_str(sptr))
    end function

    function cea_detonation_solver_create_with_options(sptr, mptr_prod, opts) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: sptr
        type(c_ptr), intent(in), value :: mptr_prod
        type(cea_solver_opts), intent(in), value :: opts

        ! Locals
        type(DetonSolver), pointer :: solver
        type(Mixture), pointer :: products
        type(Mixture), pointer :: reactants
        type(c_ptr), pointer :: cinsert(:)
        character(snl), allocatable :: insert(:)
        character(:), allocatable :: name
        integer :: n
        logical :: ions, transport, use_trace

        ierr = CEA_SUCCESS
        call c_f_pointer(mptr_prod, products)

        if (opts%ninsert < 0) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        if (opts%ninsert > 0 .and. .not. c_associated(opts%insert)) then
            ierr = CEA_INVALID_SIZE
            return
        end if

        allocate(cinsert(opts%ninsert), insert(opts%ninsert))

        ! Handle optional reactants mixture
        if (c_associated(opts%reactants)) then
            call c_f_pointer(opts%reactants, reactants)
        else
            nullify(reactants)
        end if

        ! Handle trace option
        if (opts%trace > 0.0d0) then
            use_trace = .true.
        else
            use_trace = .false.
        end if

        ! Convert ions from C bool to Fortran logical
        ions = logical(opts%ions)
        transport = logical(opts%transport)

        ! Convert insert species names from C strings to Fortran strings
        if (opts%ninsert > 0 .and. c_associated(opts%insert)) then
            call c_f_pointer(opts%insert, cinsert, [opts%ninsert])

            do n = 1, opts%ninsert
                call c_copy(cinsert(n), name)
                if (len(name) > snl) then
                    ierr = CEA_INVALID_SIZE
                    return
                end if
                insert(n) = name
            end do
        end if

        allocate(solver)
        if (associated(reactants)) then
            if (transport) then
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = DetonSolver(products=products, reactants=reactants, trace=opts%trace, &
                                             ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = DetonSolver(products=products, reactants=reactants, trace=opts%trace, &
                                             ions=ions, all_transport=global_transdb)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = DetonSolver(products=products, reactants=reactants, &
                                            ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = DetonSolver(products=products, reactants=reactants, &
                                            ions=ions, all_transport=global_transdb)
                    end if
                end if
            else
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = DetonSolver(products=products, reactants=reactants, trace=opts%trace, &
                                            ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = DetonSolver(products=products, reactants=reactants, trace=opts%trace, &
                                             ions=ions)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = DetonSolver(products=products, reactants=reactants, &
                                             ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = DetonSolver(products=products, reactants=reactants, &
                                             ions=ions)
                    end if
                end if
            end if
        else
            if (transport) then
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = DetonSolver(products=products, trace=opts%trace, &
                                             ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = DetonSolver(products=products, trace=opts%trace, &
                                             ions=ions, all_transport=global_transdb)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = DetonSolver(products=products, &
                                             ions=ions, all_transport=global_transdb, insert=insert(1:opts%ninsert))
                    else
                        solver = DetonSolver(products=products, &
                                             ions=ions, all_transport=global_transdb)
                    end if
                end if
            else
                if (use_trace) then
                    if (opts%ninsert > 0) then
                        solver = DetonSolver(products=products, trace=opts%trace, &
                                             ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = DetonSolver(products=products, trace=opts%trace, &
                                             ions=ions)
                    end if
                else
                    if (opts%ninsert > 0) then
                        solver = DetonSolver(products=products, &
                                             ions=ions, insert=insert(1:opts%ninsert))
                    else
                        solver = DetonSolver(products=products, &
                                             ions=ions)
                    end if
                end if
            end if
        end if
        sptr = c_loc(solver)
        call log_info('BINDC: Created DetonSolver with options at '//to_str(sptr))
    end function

    function cea_detonation_solver_destroy(sptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(inout) :: sptr
        type(DetonSolver), pointer :: solver
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr, solver)
        deallocate(solver)
        call log_info('BINDC: Destroyed DetonSolver object at '//to_str(sptr))
        sptr = c_null_ptr
    end function

    function cea_detonation_solver_get_size(sptr, eq_variable, eq_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr),    intent(in), value :: sptr
        integer(c_int), intent(in), value :: eq_variable
        integer(c_int), intent(out) :: eq_value
        type(DetonSolver), pointer :: solver
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr,  solver)
        select case(eq_variable)
            case (CEA_NUM_REACTANTS); eq_value = solver%eq_solver%num_reactants
            case (CEA_NUM_PRODUCTS);  eq_value = solver%eq_solver%num_products
            case (CEA_NUM_GAS);       eq_value = solver%eq_solver%num_gas
            case (CEA_NUM_CONDENSED); eq_value = solver%eq_solver%num_condensed
            case (CEA_NUM_ELEMENTS);  eq_value = solver%eq_solver%num_elements
            case (CEA_MAX_EQUATIONS); eq_value = solver%eq_solver%max_equations
            case default
                eq_value = -1
                ierr = CEA_INVALID_EQUILIBRIUM_SIZE_TYPE
                return
        end select
    end function

    function cea_detonation_solver_solve(sptr, slptr, weights, T1, p1, frozen) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr),     intent(in), value :: sptr
        type(c_ptr),     intent(in), value :: slptr
        real(c_double),  intent(in)        :: weights(*)
        real(c_double),  intent(in), value :: T1
        real(c_double),  intent(in), value :: p1
        logical(c_bool), intent(in), value :: frozen

        type(DetonSolver), pointer :: solver
        type(DetonSolution), pointer :: solution
        logical :: frozen_
        integer :: nr

        ierr = CEA_SUCCESS

        call c_f_pointer(sptr,  solver)
        call c_f_pointer(slptr, solution)

        nr = solver%eq_solver%num_reactants

        frozen_ = .false.
        if (frozen .eqv. .true.) frozen_ = .true.

        solution = solver%solve(weights(:nr), T1, p1, frozen=frozen_)
        if (.not. solution%converged) ierr = CEA_NOT_CONVERGED
    end function

    !-----------------------------------------------------------------
    ! Equilibrium Solution
    !-----------------------------------------------------------------
    function cea_eqsolution_create(slptr, sptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: slptr
        type(c_ptr), intent(in), value  :: sptr
        type(EqSolution), pointer :: solution
        type(EqSolver), pointer :: solver
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr, solver)
        allocate(solution)
        solution = EqSolution(solver)
        slptr = c_loc(solution)
        call log_info('BINDC: Created EqSolution object at '//to_str(slptr))
    end function

    function cea_eqsolution_destroy(slptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(inout) :: slptr
        type(EqSolution), pointer :: solution
        ierr = CEA_SUCCESS
        if (.not. c_associated(slptr)) then
            slptr = c_null_ptr
            return
        end if
        call c_f_pointer(slptr, solution)
        if (associated(solution)) deallocate(solution)
        slptr = c_null_ptr
        call log_info('BINDC: Destroyed SolutionEq object at '//to_str(slptr))
    end function

    function cea_eqsolution_get_property(slptr, prop_type, prop_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: prop_type
        real(c_double), intent(out) :: prop_value
        type(EqSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        select case(prop_type)
            case (CEA_TEMPERATURE)
                prop_value = solution%T
            case (CEA_PRESSURE)
                prop_value = solution%pressure
            case (CEA_VOLUME)
                prop_value = solution%volume
            case (CEA_DENSITY)
                prop_value = solution%density
            case (CEA_M)
                prop_value = solution%M
            case (CEA_MW)
                prop_value = solution%MW
            case (CEA_ENTHALPY)
                prop_value = solution%enthalpy
            case (CEA_ENERGY)
                prop_value = solution%energy
            case (CEA_ENTROPY)
                prop_value = solution%entropy
            case (CEA_GIBBS_ENERGY)
                prop_value = solution%gibbs_energy
            case (CEA_GAMMA_S)
                prop_value = solution%gamma_s
            case (CEA_FROZEN_CP)
                prop_value = solution%cp_fr
            case (CEA_FROZEN_CV)
                prop_value = solution%cv_fr
            case (CEA_EQUILIBRIUM_CP)
                prop_value = solution%cp_eq
            case (CEA_EQUILIBRIUM_CV)
                prop_value = solution%cv_eq
            ! Transport properties
            case (CEA_VISCOSITY)
                prop_value = solution%viscosity
            case (CEA_FROZEN_CONDUCTIVITY)
                prop_value = solution%conductivity_fr
            case (CEA_EQUILIBRIUM_CONDUCTIVITY)
                prop_value = solution%conductivity_eq
            case (CEA_FROZEN_PRANDTL)
                prop_value = solution%Pr_fr
            case (CEA_EQUILIBRIUM_PRANDTL)
                prop_value = solution%Pr_eq
            case default
                prop_value = empty_dp
                ierr = CEA_INVALID_PROPERTY_TYPE
        end select
    end function

    function cea_eqsolution_get_weights(slptr, np, weights, log) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: np
        real(c_double), intent(out) :: weights(*)
        logical(c_bool), intent(in), value :: log
        type(EqSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        if (log .eqv. .true.) then
            weights(:np) = solution%ln_nj(:np)
        else
            weights(:np) = solution%nj(:np)
        end if
    end function

    function cea_eqsolution_set_T(slptr, T) result(ierr) bind(c, name="cea_eqsolution_set_T")
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        real(c_double), intent(in), value :: T
        type(EqSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        solution%T = T
    end function

    function cea_eqsolution_set_nj(slptr, sptr, np, nj) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        type(c_ptr), intent(in), value :: sptr
        integer(c_int), intent(in), value :: np
        real(c_double), intent(in) :: nj(*)
        type(EqSolution), pointer :: solution
        type(EqSolver), pointer :: solver
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        call c_f_pointer(sptr, solver)
        if (np /= solver%num_products) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        call solution%set_nj(solver, nj(:np))
    end function

    function cea_eqsolution_get_species_amounts(slptr, np, amounts, mass) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: np
        real(c_double), intent(out) :: amounts(*)
        logical(c_bool), intent(in), value :: mass
        type(EqSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        if (mass .eqv. .true.) then
            amounts(:np) = solution%mass_fractions(:np)
        else
            amounts(:np) = solution%mole_fractions(:np)
        end if
    end function

    function cea_eqsolution_get_moles(slptr, n) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        real(c_double), intent(out) :: n
        type(EqSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        n = solution%n
    end function

    function cea_eqsolution_get_converged(slptr, converged) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(out) :: converged
        type(EqSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        if (solution%converged) then
            converged = 1
        else
            converged = 0
        end if
    end function

    !-----------------------------------------------------------------
    ! Equilibrium Partials
    !-----------------------------------------------------------------
    function cea_eqpartials_create(pptr, sptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: pptr
        type(c_ptr), intent(in), value  :: sptr
        type(EqPartials), pointer :: partials
        type(EqSolver), pointer :: solver
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr, solver)
        allocate(partials)
        pptr = c_loc(partials)
        call log_info('BINDC: Created EqPartials object at '//to_str(pptr))
    end function

    function cea_eqpartials_destroy(pptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(inout) :: pptr
        type(EqPartials), pointer :: partials
        ierr = CEA_SUCCESS
        call c_f_pointer(pptr, partials)
        deallocate(partials)
        call log_info('BINDC: Destroyed EqPartials object at '//to_str(pptr))
        pptr = c_null_ptr
    end function

    !-----------------------------------------------------------------
    ! Rocket Solution
    !-----------------------------------------------------------------
    function cea_rocket_solution_create(slptr, sptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: slptr
        type(c_ptr), intent(in), value  :: sptr
        type(RocketSolution), pointer :: solution
        type(RocketSolver), pointer :: solver
        ierr = CEA_SUCCESS
        call c_f_pointer(sptr, solver)
        allocate(solution)
        solution = RocketSolution(solver)
        slptr = c_loc(solution)
        call log_info('BINDC: Created RocketSolution object at '//to_str(slptr))
    end function

    function cea_rocket_solution_destroy(slptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(inout) :: slptr
        type(RocketSolution), pointer :: solution
        ierr = CEA_SUCCESS
        if (.not. c_associated(slptr)) then
            slptr = c_null_ptr
            return
        end if
        call c_f_pointer(slptr, solution)
        if (associated(solution)) deallocate(solution)
        call log_info('BINDC: Destroyed RocketSolution object at '//to_str(slptr))
        slptr = c_null_ptr
    end function

    function cea_rocket_solution_get_size(slptr, num_pts) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr),    intent(in), value :: slptr
        integer(c_int), intent(out) :: num_pts
        type(RocketSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        num_pts = solution%num_pts
    end function

    function cea_rocket_solution_get_property(slptr, prop_type, len, prop_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: prop_type
        integer(c_int), intent(in), value :: len
        real(c_double), intent(out) :: prop_value(*)
        type(RocketSolution), pointer :: solution
        integer(c_int) :: num_pts
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        num_pts = solution%num_pts
        if (len < num_pts) then
            ierr = CEA_INVALID_SIZE
            return
        end if

        select case(prop_type)
            case (CEA_ROCKET_TEMPERATURE)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%T
            case (CEA_ROCKET_PRESSURE)
                prop_value(:num_pts) = solution%pressure(:num_pts)
            case (CEA_ROCKET_VOLUME)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%volume
            case (CEA_ROCKET_DENSITY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%density
            case (CEA_ROCKET_M)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%M
            case (CEA_ROCKET_MW)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%MW
            case (CEA_ROCKET_ENTHALPY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%enthalpy
            case (CEA_ROCKET_ENERGY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%energy
            case (CEA_ROCKET_ENTROPY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%entropy
            case (CEA_ROCKET_GIBBS_ENERGY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%gibbs_energy
            case (CEA_ROCKET_GAMMA_S)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%gamma_s
            case (CEA_ROCKET_FROZEN_CP)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%cp_fr
            case (CEA_ROCKET_FROZEN_CV)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%cv_fr
            case (CEA_ROCKET_EQUILIBRIUM_CP)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%cp_eq
            case (CEA_ROCKET_EQUILIBRIUM_CV)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%cv_eq
            case (CEA_MACH)
                prop_value(:num_pts) = solution%mach(:num_pts)
            case (CEA_SONIC_VELOCITY)
                prop_value(:num_pts) = solution%v_sonic(:num_pts)
            case (CEA_AE_AT)
                prop_value(:num_pts) = solution%ae_at(:num_pts)
            case (CEA_C_STAR)
                prop_value(:num_pts) = solution%c_star(:num_pts)
            case (CEA_COEFFICIENT_OF_THRUST)
                prop_value(:num_pts) = solution%cf(:num_pts)
            case (CEA_ISP)
                prop_value(:num_pts) = solution%i_sp(:num_pts)
            case (CEA_ISP_VACUUM)
                prop_value(:num_pts) = solution%i_vac(:num_pts)
            case (CEA_ROCKET_VISCOSITY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%viscosity
            case (CEA_ROCKET_FROZEN_CONDUCTIVITY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%conductivity_fr
            case (CEA_ROCKET_EQUILIBRIUM_CONDUCTIVITY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%conductivity_eq
            case (CEA_ROCKET_FROZEN_PRANDTL)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%Pr_fr
            case (CEA_ROCKET_EQUILIBRIUM_PRANDTL)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%Pr_eq
            case default
                prop_value(:num_pts) = empty_dp
                ierr = CEA_INVALID_PROPERTY_TYPE
        end select
    end function

    function cea_rocket_solution_get_weights(slptr, np, station, weights, log) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: np
        integer(c_int), intent(in), value :: station
        real(c_double), intent(out) :: weights(*)
        logical(c_bool), intent(in), value :: log
        type(RocketSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)

        ! Check if index is valid
        if (station < 1 .or. station > size(solution%eq_soln)) then
            ierr = CEA_INVALID_INDEX
            return
        end if

        if (log .eqv. .true.) then
            if (np > size(solution%eq_soln(station)%ln_nj)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            weights(:np) = solution%eq_soln(station)%ln_nj(:np)
        else
            if (np > size(solution%eq_soln(station)%nj)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            weights(:np) = solution%eq_soln(station)%nj(:np)
        end if
    end function

    function cea_rocket_solution_get_species_amounts(slptr, np, station, amounts, mass) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: np
        integer(c_int), intent(in), value :: station
        real(c_double), intent(out) :: amounts(*)
        logical(c_bool), intent(in), value :: mass
        type(RocketSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)

        ! Check if index is valid
        if (station < 1 .or. station > size(solution%eq_soln)) then
            ierr = CEA_INVALID_INDEX
            return
        end if

        if (mass .eqv. .true.) then
            if (np > size(solution%eq_soln(station)%mass_fractions)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            amounts(:np) = solution%eq_soln(station)%mass_fractions(:np)
        else
            if (np > size(solution%eq_soln(station)%mole_fractions)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            amounts(:np) = solution%eq_soln(station)%mole_fractions(:np)
        end if
    end function

    function cea_rocket_solution_get_moles(slptr, n) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        real(c_double), intent(out) :: n(*)
        type(RocketSolution), pointer :: solution
        integer(c_int) :: num_pts
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        num_pts = solution%num_pts
        n(:num_pts) = solution%eq_soln(:num_pts)%n
    end function

    function cea_rocket_solution_get_converged(slptr, converged) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(out) :: converged
        type(RocketSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        if (solution%converged) then
            converged = 1
        else
            converged = 0
        end if
    end function

    ! function cea_rocket_solution_get_eq_solutions(slptr, npts, eqslptrs) result(ierr) bind(c)
    !     integer(c_int) :: ierr
    !     type(c_ptr), intent(in), value :: slptr
    !     integer(c_int), intent(in), value :: npts
    !     type(c_ptr), intent(out) :: eqslptrs(npts)
    !     type(RocketSolution), pointer :: solution
    !     integer :: n
    !     ierr = CEA_SUCCESS
    !     call c_f_pointer(slptr, solution)
    !     do n = 1, npts
    !         eqslptrs(n) = c_loc(solution%eq_soln(n))
    !     end do
    ! end function

    ! function cea_rocket_solution_destroy_eq_solutions(slptr, npts, eqslptrs) result(ierr) bind(c)
    !     integer(c_int) :: ierr
    !     type(c_ptr), intent(in), value :: slptr
    !     integer(c_int), intent(in), value :: npts
    !     type(c_ptr), intent(out) :: eqslptrs(npts)
    !     type(RocketSolution), pointer :: solution
    !     integer :: n
    !     ierr = CEA_SUCCESS
    !     call c_f_pointer(slptr, solution)
    !     deallocate(solution%eq_soln)
    !     do n = 1, npts
    !         eqslptrs(n) = c_null_ptr
    !     end do
    ! end function

    !-----------------------------------------------------------------
    ! Shock Solution
    !-----------------------------------------------------------------
    function cea_shock_solution_create(slptr, num_pts) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: slptr
        integer(c_int), intent(in), value :: num_pts
        type(ShockSolution), pointer :: solution
        ierr = CEA_SUCCESS
        allocate(solution)
        solution = ShockSolution(num_pts)
        slptr = c_loc(solution)
        call log_info('BINDC: Created ShockSolution object at '//to_str(slptr))
    end function

    function cea_shock_solution_destroy(slptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(inout) :: slptr
        type(ShockSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        deallocate(solution)
        call log_info('BINDC: Destroyed ShockSolution object at '//to_str(slptr))
        slptr = c_null_ptr
    end function

    function cea_shock_solution_get_property(slptr, prop_type, len, prop_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: prop_type
        integer(c_int), intent(in), value :: len
        real(c_double), intent(out) :: prop_value(*)
        type(ShockSolution), pointer :: solution
        integer(c_int) :: num_pts
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        num_pts = solution%num_pts
        if (len < num_pts) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        select case(prop_type)
            case (CEA_SHOCK_TEMPERATURE)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%T
            case (CEA_SHOCK_PRESSURE)
                prop_value(:num_pts) = solution%pressure(:num_pts)
            case (CEA_SHOCK_VELOCITY)
                prop_value(:num_pts) = solution%u(:num_pts)
            case (CEA_SHOCK_MACH)
                prop_value(:num_pts) = solution%mach(:num_pts)
            case (CEA_SHOCK_SONIC_VELOCITY)
                prop_value(:num_pts) = solution%v_sonic(:num_pts)
            case (CEA_SHOCK_VOLUME)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%volume
            case (CEA_SHOCK_DENSITY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%density
            case (CEA_SHOCK_M)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%M
            case (CEA_SHOCK_MW)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%MW
            case (CEA_SHOCK_ENTHALPY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%enthalpy
            case (CEA_SHOCK_ENERGY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%energy
            case (CEA_SHOCK_ENTROPY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%entropy
            case (CEA_SHOCK_GIBBS_ENERGY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%gibbs_energy
            case (CEA_SHOCK_GAMMA_S)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%gamma_s
            case (CEA_SHOCK_FROZEN_CP)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%cp_fr
            case (CEA_SHOCK_FROZEN_CV)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%cv_fr
            case (CEA_SHOCK_EQUILIBRIUM_CP)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%cp_eq
            case (CEA_SHOCK_EQUILIBRIUM_CV)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%cv_eq
            case (CEA_SHOCK_VISCOSITY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%viscosity
            case (CEA_SHOCK_FROZEN_CONDUCTIVITY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%conductivity_fr
            case (CEA_SHOCK_EQUILIBRIUM_CONDUCTIVITY)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%conductivity_eq
            case (CEA_SHOCK_FROZEN_PRANDTL)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%Pr_fr
            case (CEA_SHOCK_EQUILIBRIUM_PRANDTL)
                prop_value(:num_pts) = solution%eq_soln(:num_pts)%Pr_eq
            case default
                prop_value(:num_pts) = empty_dp
                ierr = CEA_INVALID_PROPERTY_TYPE
        end select
    end function

    function cea_shock_solution_get_scalar_property(slptr, prop_type, prop_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: prop_type
        real(c_double), intent(out) :: prop_value
        type(ShockSolution), pointer :: solution
        integer(c_int) :: num_pts
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        num_pts = solution%num_pts
        select case(prop_type)
            case (CEA_SHOCK_RHO12)
                prop_value = solution%rho12
            case (CEA_SHOCK_RHO52)
                prop_value = solution%rho52
            case (CEA_SHOCK_P21)
                prop_value = solution%p21
            case (CEA_SHOCK_P52)
                prop_value = solution%p52
            case (CEA_SHOCK_T21)
                prop_value = solution%T21
            case (CEA_SHOCK_T52)
                prop_value = solution%T52
            case (CEA_SHOCK_M21)
                prop_value = solution%M21
            case (CEA_SHOCK_M52)
                prop_value = solution%M52
            case (CEA_SHOCK_V2)
                prop_value = solution%v2
            case (CEA_SHOCK_U5_P_V2)
                prop_value = solution%u5_p_v2
            case default
                prop_value = empty_dp
                ierr = CEA_INVALID_PROPERTY_TYPE
        end select
    end function

    function cea_shock_solution_get_weights(slptr, np, index, weights, log) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: np
        integer(c_int), intent(in), value :: index
        real(c_double), intent(out) :: weights(*)
        logical(c_bool), intent(in), value :: log
        type(ShockSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)

        ! Check if index is valid
        if (index < 1 .or. index > size(solution%eq_soln)) then
            ierr = CEA_INVALID_INDEX
            return
        end if

        if (log .eqv. .true.) then
            if (np > size(solution%eq_soln(index)%ln_nj)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            weights(:np) = solution%eq_soln(index)%ln_nj(:np)
        else
            if (np > size(solution%eq_soln(index)%nj)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            weights(:np) = solution%eq_soln(index)%nj(:np)
        end if
    end function

    function cea_shock_solution_get_species_amounts(slptr, np, index, amounts, mass) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: np
        integer(c_int), intent(in), value :: index
        real(c_double), intent(out) :: amounts(*)
        logical(c_bool), intent(in), value :: mass
        type(ShockSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)

        ! Check if index is valid
        if (index < 1 .or. index > size(solution%eq_soln)) then
            ierr = CEA_INVALID_INDEX
            return
        end if

        if (mass .eqv. .true.) then
            if (np > size(solution%eq_soln(index)%mass_fractions)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            amounts(:np) = solution%eq_soln(index)%mass_fractions(:np)
        else
            if (np > size(solution%eq_soln(index)%mole_fractions)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            amounts(:np) = solution%eq_soln(index)%mole_fractions(:np)
        end if
    end function

    function cea_shock_solution_get_moles(slptr, n) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        real(c_double), intent(out) :: n(*)
        type(ShockSolution), pointer :: solution
        integer(c_int) :: num_pts
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        num_pts = solution%num_pts
        n(:num_pts) = solution%eq_soln(:num_pts)%n
    end function

    function cea_shock_solution_get_converged(slptr, converged) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(out) :: converged
        type(ShockSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        if (solution%converged) then
            converged = 1
        else
            converged = 0
        end if
    end function

    !-----------------------------------------------------------------
    ! Detonation Solution
    !-----------------------------------------------------------------
    function cea_detonation_solution_create(slptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(out) :: slptr
        type(DetonSolution), pointer :: solution
        ierr = CEA_SUCCESS
        allocate(solution)
        solution = DetonSolution()
        slptr = c_loc(solution)
        call log_info('BINDC: Created DetonSolution object at '//to_str(slptr))
    end function

    function cea_detonation_solution_destroy(slptr) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(inout) :: slptr
        type(DetonSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        deallocate(solution)
        call log_info('BINDC: Destroyed DetonSolution object at '//to_str(slptr))
        slptr = c_null_ptr
    end function

    function cea_detonation_solution_get_property(slptr, prop_type, len, prop_value) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: prop_type
        integer(c_int), intent(in), value :: len
        real(c_double), intent(out) :: prop_value
        type(DetonSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        if (len < 1) then
            ierr = CEA_INVALID_SIZE
            return
        end if
        select case(prop_type)
            case (CEA_DETONATION_P1)
                prop_value = solution%P1
            case (CEA_DETONATION_T1)
                prop_value = solution%T1
            case (CEA_DETONATION_H1)
                prop_value = solution%H1
            case (CEA_DETONATION_M1)
                prop_value = solution%M1
            case (CEA_DETONATION_GAMMA1)
                prop_value = solution%gamma1
            case (CEA_DETONATION_V_SONIC1)
                prop_value = solution%v_sonic1
            case (CEA_DETONATION_PRESSURE)
                prop_value = solution%pressure
            case (CEA_DETONATION_TEMPERATURE)
                prop_value = solution%eq_soln%T
            case (CEA_DETONATION_DENSITY)
                prop_value = solution%eq_soln%density
            case (CEA_DETONATION_ENTHALPY)
                prop_value = solution%eq_soln%enthalpy
            case (CEA_DETONATION_ENERGY)
                prop_value = solution%eq_soln%energy
            case (CEA_DETONATION_GIBBS_ENERGY)
                prop_value = solution%eq_soln%gibbs_energy
            case (CEA_DETONATION_ENTROPY)
                prop_value = solution%eq_soln%entropy
            case (CEA_DETONATION_MACH)
                prop_value = solution%mach
            case (CEA_DETONATION_VELOCITY)
                prop_value = solution%velocity
            case (CEA_DETONATION_SONIC_VELOCITY)
                prop_value = solution%sonic_velocity
            case (CEA_DETONATION_GAMMA)
                prop_value = solution%gamma
            case (CEA_DETONATION_P_P1)
                prop_value = solution%P_P1
            case (CEA_DETONATION_T_T1)
                prop_value = solution%T_T1
            case (CEA_DETONATION_M_M1)
                prop_value = solution%M_M1
            case (CEA_DETONATION_RHO_RHO1)
                prop_value = solution%rho_rho1
            case (CEA_DETONATION_FROZEN_CP)
                prop_value = solution%eq_soln%cp_fr
            case (CEA_DETONATION_FROZEN_CV)
                prop_value = solution%eq_soln%cv_fr
            case (CEA_DETONATION_EQUILIBRIUM_CP)
                prop_value = solution%eq_soln%cp_eq
            case (CEA_DETONATION_EQUILIBRIUM_CV)
                prop_value = solution%eq_soln%cv_eq
            case (CEA_DETONATION_M)
                prop_value = solution%eq_soln%M
            case (CEA_DETONATION_MW)
                prop_value = solution%eq_soln%MW
            case (CEA_DETONATION_VISCOSITY)
                prop_value = solution%eq_soln%viscosity
            case (CEA_DETONATION_FROZEN_CONDUCTIVITY)
                prop_value = solution%eq_soln%conductivity_fr
            case (CEA_DETONATION_EQUILIBRIUM_CONDUCTIVITY)
                prop_value = solution%eq_soln%conductivity_eq
            case (CEA_DETONATION_FROZEN_PRANDTL)
                prop_value = solution%eq_soln%Pr_fr
            case (CEA_DETONATION_EQUILIBRIUM_PRANDTL)
                prop_value = solution%eq_soln%Pr_eq
            case default
                prop_value = empty_dp
                ierr = CEA_INVALID_PROPERTY_TYPE
        end select
    end function

    function cea_detonation_solution_get_weights(slptr, np, weights, log) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: np
        real(c_double), intent(out) :: weights(*)
        logical(c_bool), intent(in), value :: log
        type(DetonSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)

        if (log .eqv. .true.) then
            if (np > size(solution%eq_soln%ln_nj)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            weights(:np) = solution%eq_soln%ln_nj(:np)
        else
            if (np > size(solution%eq_soln%nj)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            weights(:np) = solution%eq_soln%nj(:np)
        end if
    end function

    function cea_detonation_solution_get_species_amounts(slptr, np, amounts, mass) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(in), value :: np
        real(c_double), intent(out) :: amounts(*)
        logical(c_bool), intent(in), value :: mass
        type(DetonSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)

        if (mass .eqv. .true.) then
            if (np > size(solution%eq_soln%mass_fractions)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            amounts(:np) = solution%eq_soln%mass_fractions(:np)
        else
            if (np > size(solution%eq_soln%mole_fractions)) then
                ierr = CEA_INVALID_SIZE
                return
            end if
            amounts(:np) = solution%eq_soln%mole_fractions(:np)
        end if
    end function

    function cea_detonation_solution_get_moles(slptr, n) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        real(c_double), intent(out) :: n
        type(DetonSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        n = solution%eq_soln%n
    end function

    function cea_detonation_solution_get_converged(slptr, converged) result(ierr) bind(c)
        integer(c_int) :: ierr
        type(c_ptr), intent(in), value :: slptr
        integer(c_int), intent(out) :: converged
        type(DetonSolution), pointer :: solution
        ierr = CEA_SUCCESS
        call c_f_pointer(slptr, solution)
        if (solution%converged) then
            converged = 1
        else
            converged = 0
        end if
    end function

    !-----------------------------------------------------------------
    ! Helper Functions
    !-----------------------------------------------------------------
    logical function thermodb_has_species(name) result(found)
        character(*), intent(in) :: name
        integer :: i

        found = .false.
        do i = 1, global_thermodb%num_products
            if (names_match(name, global_thermodb%product_name_list(i))) then
                found = .true.
                return
            end if
        end do
        do i = 1, global_thermodb%num_reactants
            if (names_match(name, global_thermodb%reactant_name_list(i))) then
                found = .true.
                return
            end if
        end do
    end function

    function c_len_cstr(cstr, max_len) result(n)
        ! Find the length of a null-terminated c-string
        character(c_char), intent(in) :: cstr(*)
        integer, intent(in), optional :: max_len
        integer :: n, nmax
        nmax = 2048
        if (present(max_len)) nmax = max_len
        do n = 1,nmax+1
            if (cstr(n) == c_null_char) exit
        end do
        n = n-1
        if (n == nmax) call log_warning('C-string truncated to max length of '//to_str(nmax))
    end function

    function c_len_cptrs(cptrs, max_len) result(n)
        ! Find the length of a null-terminated array of c_ptrs
        type(c_ptr), intent(in) :: cptrs(*)
        integer, intent(in), optional :: max_len
        integer :: n, nmax
        nmax = 2048
        if (present(max_len)) nmax = max_len
        do n = 1,nmax+1
            if (.not. c_associated(cptrs(n))) exit
        end do
        n = n-1
        if (n == nmax) call log_warning('Null-terminated c_ptr array truncated to max length of '//to_str(nmax))
    end function

    subroutine c_copy_str_cstr(cstr, fstr, max_len)
        ! Create a Fortran string from null-terminated c-string
        character(c_char), intent(in) :: cstr(*)
        character(:), allocatable, intent(out) :: fstr
        integer, intent(in), optional :: max_len
        integer :: n,len
        len = c_len(cstr, max_len)
        allocate(character(len) :: fstr)
        do n = 1,len
            fstr(n:n) = cstr(n)
        end do
    end subroutine

    subroutine c_copy_str_cptr(cptr, fstr, max_len)
        ! Create a Fortran string from a c_ptr to a null-terminated c-string
        type(c_ptr), intent(in) :: cptr
        character(:), allocatable, intent(out) :: fstr
        integer, intent(in), optional :: max_len
        character(c_char), pointer :: cstr(:)
        integer :: nmax
        nmax = 2048
        if (present(max_len)) nmax = max_len
        call c_f_pointer(cptr, cstr, [nmax+1])
        call c_copy(cstr, fstr, max_len)
    end subroutine

    function to_str_cptr(cptr) result(str)
        ! Return string representation of a C pointer
        type(c_ptr), intent(in) :: cptr
        character(:), allocatable :: str
        character(256) :: buffer
        integer(c_int) :: address
        address = transfer(cptr, address)
        write(buffer,'(z0)') address
        str = '0x'//trim(buffer)
    end function

    ! -----------------------------------------------------------------
    ! Function inteface for Matlab and Excel
    ! -----------------------------------------------------------------
    ! function solve_eq(problem_type, state1, state2, nreactants, creactants, amounts, nproducts, cproducts, ninsert, insert, set_trace, trace, transport, ions) result(slptr) bind(c)
    !     ! Function-based equilibrium solver interface without computing partial derivatives.
    !     ! Notes:
    !     ! - This function is designed to be called from Matlab or Excel.
    !     ! - It does not use the EqSolution or EqPartials objects.
    !     ! - It returns a pointer to the EqSolution object.
    !     ! - It requires reactant weights to be passed as amounts.
    !     !   Converting from o/f ratios (or other amount types) should be done using other functions prior to calling this function.
    !     ! - It requires passing in "state1" (e.g. temperature, enthalpy, entropy, or energy).
    !     !   Enthalpy/entropy/energy can be computed from reactant weights and temperatures using other functions prior to calling this function.
    !     ! - state1 and state2 are required to be in SI units.
    !     type(c_ptr) :: slptr
    !     integer(kind=kind(CEA_TP)), intent(in), value :: problem_type
    !     real(c_double),  intent(in), value :: state1
    !     real(c_double),  intent(in), value :: state2
    !     integer(c_int),  intent(in), value :: nreactants
    !     type(c_ptr),     intent(in)        :: creactants(*)
    !     real(c_double),  intent(in)        :: amounts(*)
    !     integer(c_int),  intent(in), value :: nproducts
    !     type(c_ptr),     intent(in)        :: cproducts(*)
    !     integer(c_int),  intent(in), value :: ninsert
    !     type(c_ptr),     intent(in)        :: insert(*)
    !     logical(c_bool), intent(in), value :: set_trace
    !     real(c_double),  intent(in), value :: trace
    !     logical(c_bool), intent(in), value :: transport
    !     logical(c_bool), intent(in), value :: ions

    !     ! Locals
    !     type(Mixture), pointer :: prod, reac
    !     type(EqSolver), pointer :: solver
    !     type(EqSolution), pointer :: solution
    !     character(snl) :: reactant_names(nreactants)
    !     character(snl) :: product_names(nproducts)
    !     character(:), allocatable :: name
    !     character(2) :: eqtype
    !     integer :: n, ierr

    !     ! Convert names from C strings to fortran strings
    !     do n = 1,nreactants
    !         call c_copy(creactants(n), name)
    !         call assert(len(name) <= snl, 'Species name is longer than CEA max length')
    !         reactant_names(n) = name
    !     end do

    !     do n = 1,nproducts
    !         call c_copy(cproducts(n), name)
    !         call assert(len(name) <= snl, 'Species name is longer than CEA max length')
    !         product_names(n) = name
    !     end do

    !     ! Initialize the Mixture objects
    !     allocate(reac)
    !     allocate(prod)
    !     reac = Mixture(thermo=global_thermodb, reactant_names=reactant_names, ions=ions)
    !     prod = Mixture(thermo=global_thermodb, reactant_names=product_names)

    !     ! Initialize the EqSolver object
    !     allocate(solver)
    !     if (set_trace) then
    !         if (transport) then
    !             solver = EqSolver(prod, reac, trace, ions=ions, &
    !                               all_transport=global_transdb, insert=insert(:ninsert))
    !         else
    !             solver = EqSolver(prod, reac, trace, ions=ions, &
    !                               insert=insert(:ninsert))
    !         end if
    !     else
    !         if (transport) then
    !             solver = EqSolver(prod, reac, ions=ions, &
    !                               all_transport=global_transdb, insert=insert(:ninsert))
    !         else
    !             solver = EqSolver(prod, reac, ions=ions, &
    !                               insert=insert(:ninsert))
    !         end if
    !     end if

    !     ! Initialize the EqSolution object
    !     allocate(solution)
    !     solution = EqSolution(solver)

    !     ! Set the problem type
    !     select case(problem_type)
    !         case (CEA_TP); eqtype = 'tp'
    !         case (CEA_HP); eqtype = 'hp'
    !         case (CEA_SP); eqtype = 'sp'
    !         case (CEA_TV); eqtype = 'tv'
    !         case (CEA_UV); eqtype = 'uv'
    !         case (CEA_SV); eqtype = 'sv'
    !         case default
    !             ierr = CEA_INVALID_EQUILIBRIUM_TYPE
    !             return
    !     end select

    !     ! Solve the equilibrium problem
    !     call solver%solve(solution, eqtype, state1, state2, amounts(:nreactants))
    !     slptr = c_loc(solution)

    ! end function

    ! function solve_rocket() result(eqsol) bind(c)

    ! end function

    ! function solve_shock() result(eqsol) bind(c)

    ! end function

    ! function solve_deton() result(eqsol) bind(c)

    ! end function

end module
