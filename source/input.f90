module cea_input
    !! Module for parsing legacy input files

    use cea_param, only: dp, empty_dp, stdout, &
                         gas_constant, &
                         sn => species_name_len, &
                         en => element_name_len
    use fb_utils, only: abort, assert, substring, startswith, endswith, to_str
    use fb_string_scanner
    use fb_logging
    implicit none

    type :: Schedule
        !! A generic schedule of values with associated name and units
        character(:), allocatable :: name
        character(:), allocatable :: units
        real(dp), allocatable :: values(:)
    end type

    type :: Formula
        !! Chemical formula with elements and coefficients
        character(en), allocatable :: elements(:)
        real(dp), allocatable :: coefficients(:)
    end type

    type :: ReactantInput
        !! Reactant input data

        ! Required items
        character(:), allocatable :: type
            !! Type of reactant: 'ox', 'fu', or 'na'
        character(:), allocatable :: name
            !! Name of the reactant
        type(Schedule), allocatable :: amount
            !! Amount of the reactant

        ! Optional items
        type(Formula),  allocatable :: formula
            !! Chemical formula of the reactant (optional)
        real(dp),  allocatable :: molecular_weight
            !! Molecular weight of the reactant (optional)
        type(Schedule), allocatable :: temperature
            !! Temperature of the reactant (optional)
        type(Schedule), allocatable :: enthalpy
            !! Reference enthalpy of the reactant (optional)
        type(Schedule), allocatable :: density
            !! Density of the reactant (optional)

    end type

    type :: ProblemDataset
        !! Dataset containing the problem specification

        ! Case info
        character(:), allocatable :: name
            !! Name of the problem
        character(:), allocatable :: type
            !! Type of problem: 'tp', 'hp', 'sp', 'tv', 'uv', 'sv', 'det', 'rkt', or 'shk'

        ! General options
        logical :: debug        = .false.
            !! Enable debugging output
        logical :: frozen       = .false.
            !! Frozen composition mode
        logical :: equilibrium  = .false.
            !! Equilibrium composition
        logical :: include_ions = .false.
            !! Include ions in equilibrium calculations

        ! Rocket options
        logical :: rkt_finite_area = .false.
            !! Flag if using finite area combustor
        integer :: rkt_nfrozen     = 1
            !! Rocket station past which composition is frozen
        real(dp), allocatable :: ac_at
            !! (rkt) Chamber/throat contraction ratio
        real(dp), allocatable :: mdot
            !! (rkt) Mass flow rate
        real(dp), allocatable :: tc_est
            !! (rkt) Estimate chamber temperature

        ! Shock options
        logical :: shk_incident  = .false.
            !! Compute incident shock values
        logical :: shk_reflected = .false.
            !! Compute reflected shock values

        ! Assigned state values (problem dependent)
        type(Schedule), allocatable :: of_schedule
            !! Oxidizer/fuel mixture
        type(Schedule), allocatable :: t_schedule
            !! Temperature
        type(Schedule), allocatable :: h_schedule
            !! Enthalpy
        type(Schedule), allocatable :: u_schedule
            !! Energy
        type(Schedule), allocatable :: s_schedule
            !! Entropy
        type(Schedule), allocatable :: p_schedule
            !! Pressure
        type(Schedule), allocatable :: v_schedule
            !! Specific volume
        type(Schedule), allocatable :: u1_schedule
            !! (shk) Initial shock velocity
        type(Schedule), allocatable :: mach1_schedule
            !! (shk) Initial shock Mach number
        type(Schedule), allocatable :: pcp_schedule
            !! (rkt) Chamber/exit pressure ratio
        type(Schedule), allocatable :: subar_schedule
            !! (rkt) Subsonic nozzle expansion ratio
        type(Schedule), allocatable :: supar_schedule
            !! (rkt) Supersonic nozzle expansion ratio

    end type

    type :: OutputDataset
        !! Dataset containing the output options
        real(dp), allocatable :: trace
            !! Only output species with concentrations greater than `trace`
        logical :: transport = .false.
            !! Optional flag to calculate mixture transport properties
        logical :: mass_fractions = .false.
            !! Output mass fractions; output mole fractions if false
        logical :: debug = .false.
            !! Print debugging output
        logical :: siunit = .true.
            !! Output using SI units
    end type

    type :: ProblemDB
        !! Problem info database

        type(ProblemDataset) :: problem
            !! Problem specification
        type(OutputDataset)  :: output
            !! Output options
        type(ReactantInput), allocatable :: reactants(:)
            !! Reactant info
        character(sn), allocatable :: omit(:)
            !! Species to omit from product mixture
        character(sn), allocatable :: only(:)
            !! Exact species to include in product mixture
        character(sn), allocatable :: insert(:)
            !! Species to insert into product mixture
    end type

contains

    function read_input(filename) result(problems)
        ! Reads all problem specifications from the CEA input file

        ! Inputs
        character(*), intent(in) :: filename

        ! Return
        type(ProblemDB), allocatable :: problems(:)

        ! Locals
        integer, parameter :: max_problems = 100
        integer :: fin, n, ierr

        call log_info('Parsing input file: '//trim(filename))
        open(newunit=fin, file=filename, &
             status="old", action="read", form="formatted")
        allocate(problems(max_problems))

        n = 1
        do
            call log_info('Parsing problem specification '//to_str(n))
            call assert(n <= max_problems, "Maximum problem count exceeded.")
            problems(n) = read_problem(fin, ierr)
            if (ierr /= 0) exit
            n = n+1
        end do

        problems = problems(1:n-1)
        call log_info('Parsed '//to_str(n-1)//' problems from '//trim(filename))

        close(fin)

        return
    end function

    function read_problem(fin, ierr) result(problem)
        ! Reads a complete ProblemDB from the input stream.
        ! Aborts program if malformed input is discovered.
        ! Returns ierr < 0 if no further problems in the stream.

        integer, intent(in) :: fin
        integer, intent(out) :: ierr
        type(ProblemDB) :: problem

        character(512) :: line
        character(:), allocatable :: buffer, dsname
        logical :: has_data, has_prob, has_reac

        has_data = .false.
        has_prob = .false.
        has_reac = .false.
        do

            ! Get next good line of input
            read(fin, '(a)', iostat=ierr) line
            if (ierr /= 0) then
                if (has_data) then
                    call abort('Problem is missing end keyword.')
                else
                    call log_info('No problem found. Reached EOF.')
                    return
                end if
            end if
            line = adjustl(line)
            if (is_empty(line)) cycle
            if (startswith(line,'end')) then
                call assert(has_prob, 'Problem is missing prob dataset.')
                call assert(has_reac, 'Problem is missing reac dataset.')
                return
            end if

            ! Extract dataset name
            call assert(is_keyword(line), 'Problem has non-empty line outside of dataset.')
            dsname = line(1:4)
            has_data = .true.

            if (dsname == 'ther' .or. dsname == 'tran') then
                ! These need special handling. Shouldn't read into a buffer.
                call abort('Thermo/transport databases not yet implemented')
            end if

            ! Read until next keyword, saving input
            buffer = trim(line)
            do
                read(fin, '(A)', iostat=ierr) line
                if (ierr /= 0) call abort('Problem has incomplete dataset: '//dsname)
                line = adjustl(line)
                if (is_empty(line)) cycle
                if (is_keyword(line)) then
                    backspace(fin)
                    exit
                else
                    buffer = buffer // ' ' // trim(line)
                end if
            end do

            ! Process buffer
            call log_debug('Parsing dataset '//dsname)
            select case(dsname)
                case('prob')
                    has_prob = .true.
                    problem%problem = parse_prob(buffer)
                case('reac')
                    has_reac = .true.
                    problem%reactants = parse_reac(buffer)
                case('only')
                    problem%only = parse_species(buffer)
                case('omit')
                    problem%omit = parse_species(buffer)
                case('inse')
                    problem%insert = parse_species(buffer)
                case('outp')
                    problem%output = parse_outp(buffer)
                case default
                    call abort('Encountered unrecognized dataset: '//dsname)
            end select

        end do
    end function

    function parse_prob(buffer) result(prob)
        ! Parse the prob dataset for a given problem specification

        character(*), intent(in) :: buffer
        type(ProblemDataset) :: prob

        type(string_scanner) :: scanner
        character(15) :: token
        logical :: match
        integer :: ierr

        scanner = string_scanner(replace_delimiters(buffer))
        token = scanner%read_word()  ! Skip "prob" keyword
        do
            token = scanner%read_word(ierr)
            if (ierr < 0) exit  ! Buffer exhausted
            call log_debug('Parsing prob literal '//token)

            ! Four letter keywords
            match = .true.
            select case(token(:4))
                case('ions');  prob%include_ions = .true.
                case('tces');  prob%tc_est = scanner%read_real(ierr)!parse_schedule(scanner, token)
                case('case');  prob%name = scanner%read_word()
                case('mach');  prob%mach1_schedule = parse_schedule(scanner, token)
                case('mdot');  prob%mdot = scanner%read_real(ierr)!parse_schedule(scanner, token)
                case default
                    match = .false.
            end select
            if (match) cycle

            ! Three letter keywords
            match = .true.
            select case(token(:3))
                case('det');       prob%type = 'det'
                case('dgb','deb'); prob%debug = .true.
                case('f/o','f/a'); prob%of_schedule = parse_schedule(scanner, token)
                case('fac');       prob%rkt_finite_area = .true.
                case('h/r');       prob%h_schedule = parse_schedule(scanner, token)
                case('inc');       prob%shk_incident = .true.
                case('nfz','nfr'); prob%rkt_nfrozen = scanner%read_int(ierr)
                case('o/f');       prob%of_schedule = parse_schedule(scanner, token)
                case('phi');       prob%of_schedule = parse_schedule(scanner, token)
                case('ref');       prob%shk_reflected = .true.
                case('rho');       prob%v_schedule = parse_schedule(scanner, token)
                case('rkt');       prob%type = 'rkt'
                case('sub');       prob%subar_schedule = parse_schedule(scanner, token)
                case('sup');       prob%supar_schedule = parse_schedule(scanner, token)
                case('s/r');       prob%s_schedule = parse_schedule(scanner, token)
                case('u/r');       prob%u_schedule = parse_schedule(scanner, token)
                case default
                    match = .false.
            end select
            if (match) cycle

            ! Two letter keywords
            match = .true.
            select case(token(:2))
                case('%f');       prob%of_schedule = parse_schedule(scanner, token)
                case('ac');       prob%ac_at = scanner%read_real(ierr)
                case('eq');       prob%equilibrium = .true.
                case('fr','fz');  prob%frozen = .true.
                case('hp','ph');  prob%type = 'hp'
                case('pc','pi');  prob%pcp_schedule = parse_schedule(scanner, token)
                case('ro');       prob%type = 'rkt'
                case('sh');       prob%type = 'shk'
                case('sp','ps');  prob%type = 'sp'
                case('sv','vs');  prob%type = 'sv'
                case('tp','pt');  prob%type = 'tp'
                case('tv','vt');  prob%type = 'tv'
                case('uv','vu');  prob%type = 'uv'
                case('u1');       prob%u1_schedule = parse_schedule(scanner, token)
                case('ma');       prob%mdot = scanner%read_real(ierr)
                case default
                    match = .false.
            end select
            if (match) cycle

            ! One-letter tokens
            match = .true.
            select case(token(:1))
                case('t');  prob%t_schedule = parse_schedule(scanner, token)
                case('p');  prob%p_schedule = parse_schedule(scanner, token)
                case('v');  prob%v_schedule = parse_schedule(scanner, token)
                case('r');  prob%of_schedule = parse_schedule(scanner, token)
                case default
                    match = .false.
            end select
            if (match) cycle

            ! Abort if no token match found
            call abort('prob dataset contains unrecognized token: '//token)

        end do

        ! Validate input data
        call assert(allocated(prob%type), 'prob dataset missing problem type.')

        if (allocated(prob%ac_at)) then
            call assert(prob%ac_at /= empty_dp, 'prob dataset missing value for ac.')
        end if
        if (allocated(prob%mdot)) then
            call assert(prob%mdot /= empty_dp, 'prob dataset missing value for mdot.')
        end if
        if (allocated(prob%tc_est)) then
            call assert(prob%tc_est /= empty_dp, 'prob dataset missing value for tces.')
        end if

        if (allocated(prob%of_schedule)) then
            call assert(size(prob%of_schedule%values) > 0, 'prob dataset missing o/f values.')
        end if
        if (allocated(prob%t_schedule)) then
            call assert(size(prob%t_schedule%values) > 0, 'prob dataset missing t values.')
        end if
        if (allocated(prob%h_schedule)) then
            call assert(size(prob%h_schedule%values) > 0, 'prob dataset missing h values.')
        end if
        if (allocated(prob%u_schedule)) then
            call assert(size(prob%u_schedule%values) > 0, 'prob dataset missing u values.')
        end if
        if (allocated(prob%s_schedule)) then
            call assert(size(prob%s_schedule%values) > 0, 'prob dataset missing s values.')
        end if
        if (allocated(prob%p_schedule)) then
            call assert(size(prob%p_schedule%values) > 0, 'prob dataset missing p values.')
        end if
        if (allocated(prob%v_schedule)) then
            call assert(size(prob%v_schedule%values) > 0, 'prob dataset missing v values.')
        end if
        if (allocated(prob%u1_schedule)) then
            call assert(size(prob%u1_schedule%values) > 0, 'prob dataset missing u1 values.')
        end if
        if (allocated(prob%mach1_schedule)) then
            call assert(size(prob%mach1_schedule%values) > 0, 'prob dataset missing mach values.')
        end if
        if (allocated(prob%pcp_schedule)) then
            call assert(size(prob%pcp_schedule%values) > 0, 'prob dataset missing pcp values.')
        end if
        if (allocated(prob%subar_schedule)) then
            call assert(size(prob%subar_schedule%values) > 0, 'prob dataset missing sub values.')
        end if
        if (allocated(prob%supar_schedule)) then
            call assert(size(prob%supar_schedule%values) > 0, 'prob dataset missing sup values.')
        end if

        return
    end function

    function parse_reac(buffer) result(reac)
        ! Pars the reac dataset for a given problem specification

        character(*), intent(in) :: buffer
        !type(ReactantDataset) :: reac
        type(ReactantInput), allocatable :: reac(:)

        type(string_scanner) :: scanner
        character(15) :: token
        integer :: ierr, n, capacity, i
        logical :: has_na, has_fuel_oxid

        scanner = string_scanner(replace_delimiters(buffer, replace_commas=.false.))
        token = scanner%read_word()  ! Skip 'reac' keyword

        n = 0
        capacity = 32
        allocate(reac(capacity))
        do
            token = scanner%read_word(ierr)
            if (ierr < 0) exit  ! Buffer exhausted

            ! Start new reactant definition
            select case(token(:2))
                case ('fu','ox','na')
                    n = n+1
                    if (n > capacity) then
                        reac = [reac, reac]
                        capacity = size(reac)
                    end if
                    reac(n)%type = token(:2)
                    reac(n)%name = scanner%read_word()
                    call log_debug('Parsing parameters for reactant '//reac(n)%name)
                    cycle
            end select

            ! Parse reactant parameters
            if (n == 0) then
                call abort('reac dataset missing reactant definition before token: '//token)
            end if
            select case(token(:1))
                case ('m','w');  reac(n)%amount = parse_schedule(scanner, token)
                case ('t');      reac(n)%temperature = parse_schedule(scanner, token)
                case ('h','u');  reac(n)%enthalpy = parse_schedule(scanner, token)
                case ('r');      reac(n)%density = parse_schedule(scanner, token)
                case ('A':'Z');  reac(n)%formula = parse_formula(scanner, token)
                case default
                    call abort('reac dataset contains unrecognized token: '//token)
            end select

        end do

        reac = reac(:n)

        ! Validate input data
        call assert(n > 0, 'reac dataset missing reactant definitions.')
        has_na = .false.
        has_fuel_oxid = .false.
        do i = 1, n
            call assert(allocated(reac(i)%type), 'reac dataset missing reactant type.')
            call assert(len_trim(reac(i)%type) > 0, 'reac dataset has empty reactant type.')
            call assert(allocated(reac(i)%name), 'reac dataset missing reactant name.')
            call assert(len_trim(reac(i)%name) > 0, 'reac dataset has empty reactant name.')

            if (reac(i)%type == 'na') has_na = .true.
            if (reac(i)%type == 'fu' .or. reac(i)%type == 'ox') has_fuel_oxid = .true.

            if (reac(i)%type == 'na') then
                call assert(allocated(reac(i)%amount), 'reac dataset missing amount for reactant #'//to_str(i))
            end if

            if (allocated(reac(i)%amount)) then
                call assert(size(reac(i)%amount%values) > 0, &
                    'reac dataset missing amount values for reactant #'//to_str(i))
            end if
            if (allocated(reac(i)%temperature)) then
                call assert(size(reac(i)%temperature%values) > 0, &
                    'reac dataset missing temperature values for reactant #'//to_str(i))
            end if
            if (allocated(reac(i)%enthalpy)) then
                call assert(size(reac(i)%enthalpy%values) > 0, &
                    'reac dataset missing enthalpy values for reactant #'//to_str(i))
            end if
            if (allocated(reac(i)%density)) then
                call assert(size(reac(i)%density%values) > 0, &
                    'reac dataset missing density values for reactant #'//to_str(i))
            end if
        end do

        if (has_na .and. has_fuel_oxid) then
            call abort("reac dataset cannot mix 'name' reactants with 'fuel'/'oxid' reactants. " // &
                       "Use only name=... for all reactants, or only fuel=/oxid= for all reactants.")
        end if

    end function

    function parse_species(buffer) result(species)
        ! Parses species names from only/omit/insert datasets

        character(*), intent(in) :: buffer
        character(sn), allocatable :: species(:)

        type(string_scanner) :: scanner
        character(:), allocatable :: dsname
        integer :: ierr, n, capacity

        capacity = 64
        allocate(species(capacity))
        scanner = string_scanner(buffer)
        dsname = scanner%read_word()  ! Skip dataset name, e.g. omit

        n = 1
        do
            species(n) = scanner%read_word(ierr)
            if (ierr /= 0) exit
            n = n+1
            if (n > capacity) then
                ! Double array size when we hit capacity
                species = [species, species]
                capacity = size(species)
            end if
        end do

        call log_debug('Parsed '//to_str(n-1)//' species from '//dsname)
        species = species(:n-1)

        return
    end function

    function parse_outp(buffer) result(output)
        character(*), intent(in) :: buffer
        type(OutputDataset) :: output

        type(string_scanner) :: scanner
        character(15) :: token
        logical :: match
        integer :: ierr

        scanner = string_scanner(replace_delimiters(buffer))
        token = scanner%read_word()  ! Skip "prob" keyword
        do
            token = scanner%read_word(ierr)
            if (ierr < 0) exit  ! Buffer exhausted
            call log_debug('Parsing prob literal '//token)

            ! Five letter keywords
            match = .true.
            select case(token(:5))
                case('massf'); output%mass_fractions = .true.
                case default
                    match = .false.
            end select
            if (match) cycle

            ! Four letter keywords
            match = .true.
            select case(token(:4))
                case('tran'); output%transport = .true.
                case('trac'); output%trace = scanner%read_real()
                case default
                    match = .false.
            end select
            if (match) cycle

            ! Three letter keywords
            ! "deb" or "dbg"
            match = .true.
            select case(token(:3))
                case('deb', 'dbg'); output%debug = .true.
                case('cal'); output%siunit = .false.
                case default
                    match = .false.
            end select
            if (match) cycle

        end do

    end function

    function parse_schedule(scanner, token) result(sched)
        ! Parse a schedule object from a string buffer

        type(string_scanner), intent(inout) :: scanner
        character(*), intent(in) :: token
        type(Schedule) :: sched

        integer, parameter :: max_values = 64
        integer :: i,ierr
        real(dp) :: val

        allocate(sched%values(max_values))

        ! Try reading first schedule value; if fails, its a unit string
        val = scanner%peek_real(ierr)
        if (ierr == 0) then
            call parse_embedded(token, sched%name, sched%units)
        else
            sched%name = trim(token)
            sched%units = scanner%read_word(ierr)
        end if

        ! Read values from the schedule
        ! Need to peek then read so we don't consume next token
        do i = 1,max_values
            sched%values(i) = scanner%peek_real(ierr)
            if (ierr == 0) then
                sched%values(i) = scanner%read_real(ierr)
            else
                sched%values = sched%values(:i-1)
                exit
            end if
        end do

        return
    end function

    function parse_formula(scanner, token) result(f)
        type(string_scanner), intent(inout) :: scanner
        character(*), intent(in) :: token

        type(Formula) :: f
        integer, parameter :: max_values = 16
        integer :: i, n, ierr
        character(:), allocatable :: word

        allocate(f%elements(max_values))
        allocate(f%coefficients(max_values))

        n = 1
        if (len_trim(token) > en) then
            call abort('parse_formula: element symbol too long: '//trim(token))
        end if
        f%elements(n) = token
        f%coefficients(n) = scanner%read_real()

        do i = 2,size(f%elements)
            word = scanner%peek_word(ierr)
            if (ierr /= 0) exit  ! Buffer empty
            select case(word(1:1))
                case('A':'Z')
                    word = scanner%read_word()
                    if (len_trim(word) > en) then
                        call abort('parse_formula: element symbol too long: '//trim(word))
                    end if
                    f%elements(i) = word
                    f%coefficients(i) = scanner%read_real()
                    n = i
                case default
                    exit  ! Start of new keyword
            end select
        end do

        f%elements = f%elements(:n)
        f%coefficients = f%coefficients(:n)

        return
    end function

    subroutine parse_embedded(token, name, units)
        ! Parse an embedded literal to extract name and units

        character(*), intent(in) :: token
        character(:), allocatable, intent(out) :: name
        character(:), allocatable, intent(out) :: units
        character(:), allocatable :: txt
        integer :: i,j,l

        l = len(token)

        ! Try to get name,units via delimiter splitting
        ! This will often give better names
        do i = 1,l
            select case(token(i:i))
                case('-',':',',',';','(',')')
                    exit
            end select
        end do
        do j = i+1,l
            select case(token(j:j))
                case('-',':',',',';','(',')')
                    exit
            end select
        end do
        name  = token(:i-1)
        units = token(i+1:j-1)

        ! Token-specific matching rules
        ! Will override units found by delimiter-splitting
        if (l > 4) then
            txt = token(5:)
            select case (token(:4))
                case ('tces')
                    units = 'k'
                    if (substring(txt,'r')) units = 'r'
                    if (substring(txt,'c')) units = 'c'
                    if (substring(txt,'f')) units = 'f'
                    return
                case ('mdot')
                    units = 'kg/s/m**2'
                    return
                case ('mach')
                    return
            end select
        end if

        if (l > 3) then
            txt = token(4:)
            select case (token(:3))
                case ('rho')
                    name = 'rho'
                    if (substring(txt, 'kg')) units = 'kg/m**3'
                    if (substring(txt,  'g')) units = 'g/cm**3'
                    return
                case ('f/o','f/a','h/r','o/f','phi','sub','sup','u/r')
                    return
            end select
        end if

        if (l > 2) then
            select case (token(:2))
                case ('%f')
                    units = '%'
                    return
                case ('u1')
                    units = 'm/s'
                    return
                case ('ma')
                    units = 'kg/s/m**2'
                    return
            end select
        end if

        if (l > 1) then
            txt = token(2:)
            select case (token(:1))
                case ('m')
                    name = 'mole_frac'
                    if (endswith(token,'%')) units = '%'
                case ('w')
                    name = 'weight_frac'
                    if (endswith(token,'%')) units = '%'
                case ('t')
                    units = 'k'
                    if (substring(txt,'r')) units = 'r'
                    if (substring(txt,'c')) units = 'c'
                    if (substring(txt,'f')) units = 'f'
                case ('p')
                    units = 'bar'
                    if (substring(txt,'atm')) units = 'atm'
                    if (substring(txt,'psi')) units = 'psi'
                    if (substring(txt,'mmh')) units = 'mmhg'
                case ('h')
                    units = 'j'
                    if (substring(txt,'kj')) units = 'kj/mole'
                    if (substring(txt,'kc')) units = 'kcal/mole'
                    if (substring(txt, 'c')) units = 'cal/mole'
                case ('u')
                    units = 'j'
                    if (substring(txt,'kj')) units = 'kj/mole'
                    if (substring(txt,'kc')) units = 'kcal/mole'
                    if (substring(txt, 'c')) units = 'cal/mole'
                case ('v')
                    if (substring(txt, 'kg')) units = 'm**3/kg'
                    if (substring(txt,  'g')) units = 'cm**3/g'
            end select
        end if

        return
    end subroutine

    function is_empty(line) result(tf)
        character(*), intent(in) :: line
        logical :: tf
        tf = len_trim(line) == 0 .or. &
             startswith(line,'#') .or. &
             startswith(line,'!')
    end function

    function is_keyword(line) result(tf)
        character(*), intent(in) :: line
        logical :: tf
        character(4) :: word
        tf = startswith(line,'end')
        if (.not. tf .and. len_trim(line) >= 4) then
            word = line(1:4)
            tf = (word == 'prob') .or. &
                 (word == 'reac') .or. &
                 (word == 'only') .or. &
                 (word == 'omit') .or. &
                 (word == 'inse') .or. &
                 (word == 'outp') .or. &
                 (word == 'ther') .or. &
                 (word == 'tran')
        end if
    end function

    function replace_delimiters(line, replace_commas) result(out)
        character(*), intent(in)  :: line
        character(:), allocatable :: out
        logical, intent(in), optional :: replace_commas
        character(1), parameter :: tab = achar(9)
        character(1) :: c
        integer :: i
        logical :: commas

        commas = .true.
        if (present(replace_commas)) commas = replace_commas

        out = trim(line)
        do i = 1,len(out)
            c = out(i:i)
            select case(c)
                case(',')
                    if (commas) out(i:i) = ' '
                case('=',tab)
                    out(i:i) = ' '
            end select
        end do

        return
    end function

end module
