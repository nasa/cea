program cea
    use cea_param, snl=>species_name_len, &
                   enl=>element_name_len, &
                   R=>gas_constant
    use cea_thermo, only: ThermoDB, read_thermo
    use cea_transport, only: TransportDB, read_transport
    use cea_input, only: ProblemDB, read_input
    use cea_equilibrium, only: EqSolver, EqSolution, EqPartials
    use cea_rocket, only: RocketSolver, RocketSolution
    use cea_shock, only: ShockSolver, ShockSolution
    use cea_detonation, only: DetonSolver, DetonSolution
    use cea_db_compile, only: compile_thermo_database, compile_transport_database
    use cea_mixture
    use cea_units
    use fb_logging
    use fb_utils
    implicit none

    ! Locals
    character(:), allocatable :: input_file, thermo_file, trans_file
    character(:), allocatable :: compile_thermo_input, compile_trans_input
    type(ThermoDB) :: all_thermo
    type(TransportDB) :: all_transport
    type(ProblemDB), allocatable :: problems(:)
    type(ProblemDB) :: prob
    type(EqSolver) :: eq_solver
    type(EqSolution), allocatable :: eq_solutions(:,:,:)
    type(EqPartials), allocatable :: eq_partials(:,:,:)
    type(RocketSolver) :: rkt_solver
    type(RocketSolution), allocatable :: rkt_solutions(:,:,:)
    type(ShockSolver) :: shk_solver
    type(ShockSolution), allocatable :: shk_solutions(:,:,:)
    type(DetonSolver) :: det_solver
    type(DetonSolution), allocatable :: det_solutions(:,:,:)
    integer :: n
    logical :: ok

    call parse_arguments(input_file, thermo_file, trans_file, compile_thermo_input, compile_trans_input)
    call log_info('CEA Version: '//version_string)

    if (allocated(compile_thermo_input) .or. allocated(compile_trans_input)) then
        if (allocated(compile_thermo_input)) then
            call compile_thermo_database(compile_thermo_input, ok)
            if (.not. ok) call abort('Thermo database compilation failed.')
        end if
        if (allocated(compile_trans_input)) then
            call compile_transport_database(compile_trans_input, ok)
            if (.not. ok) call abort('Transport database compilation failed.')
        end if
        stop
    end if

    call log_info('Input File:  '//input_file)
    call log_info('Thermo File: '//thermo_file)

    if (.not. exists(input_file//'.inp')) then
        call abort('Could not locate input file: '//input_file)
    end if

    thermo_file = locate(thermo_file, data_dirs)
    if (len(thermo_file) == 0) then
        call abort('Could not locate thermo file: '//thermo_file)
    end if

    trans_file = locate(trans_file, data_dirs)
    if (len(trans_file) == 0) then
        call log_info('Could not locate transport file: '//trans_file)
    else
        all_transport = read_transport(trans_file)
    end if

    ! Read and parse the input and data files
    problems = read_input(input_file//'.inp')
    all_thermo = read_thermo(thermo_file)

    ! Initialize the output file
    open(1, file=input_file(1:len_trim(input_file))//".out", status="replace")

    ! Loop over all problem in the input file
    do n = 1, size(problems)
        write(log_buffer,'("Executing problem ",i0,": name=",a, ", type=",a)') &
            n, problems(n)%problem%name, problems(n)%problem%type
        call log_info(trim(log_buffer))

        ! Setup
        prob = problems(n)

        select case(prob%problem%type)
            case ("tp", "hp", "sp", "tv", "uv", "sv")
                call log_info('Solving equilibrium problem:')

                call run_thermo_problem(prob, all_thermo, eq_solver, eq_solutions, eq_partials)
                call thermo_output(1, prob, eq_solver, eq_solutions, eq_partials)
                deallocate(eq_solutions, eq_partials)

            case ("rkt")
                call log_info('Solving rocket problem:')

                call run_rocket_problem(prob, all_thermo, rkt_solver, rkt_solutions)
                call rocket_output(1, prob, rkt_solver, rkt_solutions)
                deallocate(rkt_solutions)

            case ("shk")
                call log_info('Solving shock problem:')

                call run_shock_problem(prob, all_thermo, shk_solver, shk_solutions)
                call shock_output(1, prob, shk_solver, shk_solutions)
                deallocate(shk_solutions)

            case ("det")
                call log_info('Solving detonation problem:')

                call run_detonation_problem(prob, all_thermo, det_solver, det_solutions)
                call deton_output(1, prob, det_solver, det_solutions)
                deallocate(det_solutions)

            case default
                call log_error('Problem type '//prob%problem%type//' is not supported')
                call abort
        end select

    end do

    ! Close the output file
    close(1)

contains

    subroutine parse_arguments(input_file, thermo_file, trans_file, compile_thermo_input, compile_trans_input)
        character(:), allocatable, intent(out) :: input_file
        character(:), allocatable, intent(out) :: thermo_file
        character(:), allocatable, intent(out) :: trans_file
        character(:), allocatable, intent(out) :: compile_thermo_input
        character(:), allocatable, intent(out) :: compile_trans_input
        character(:), allocatable :: arg
        integer :: n,nargs

        ! Defaults
        thermo_file = 'thermo.lib'
        trans_file = 'trans.lib'

        nargs = command_argument_count()
        !if (nargs == 0) then
        !    ! TODO: Add CEA2 style interactive mode
        !end if

        n = 1
        do while (n <= nargs)
            arg = pop_argument(n)
            select case(arg)
                case ('-v')
                    call set_log_level(log_levels%info)
                case ('-d')
                    call set_log_level(log_levels%debug)
                case ('-t','--thermo_file')
                    call log_info('Reading thermo_file from command line argument')
                    thermo_file = pop_argument(n)
                case ('-r','--trans_file')
                    call log_info('Reading trans_file from command line argument')
                    trans_file = pop_argument(n)
                case ('--compile-thermo')
                    call log_info('Reading thermo input file for compilation')
                    compile_thermo_input = pop_argument(n)
                case ('--compile-trans')
                    call log_info('Reading transport input file for compilation')
                    compile_trans_input = pop_argument(n)
                case ('-h')
                    call display_help
                    stop
                case default
                    call log_info('Reading input_file from first positional argument')
                    input_file = arg
            end select
        end do

        if (allocated(compile_thermo_input) .or. allocated(compile_trans_input)) then
            if (allocated(input_file)) then
                call log_error('input_file is not allowed with --compile-thermo/--compile-trans')
                call display_help
                call abort
            end if
            return
        end if

        if (.not. allocated(input_file)) then
            call log_error('Required argument not specified: input_file')
            call display_help
            call abort
        end if

        return
    end subroutine

    subroutine display_help
        write(stdout,'(a)') &
            'usage: cea [options] input_file',&
            '       cea --compile-thermo thermo.inp',&
            '       cea --compile-trans trans.inp',&
            '',&
            'Arguments:',&
            '  input_file   CEA free-form program input file (see NASA RP-1311)',&
            '',&
            'Options:',&
            '  -h   Display command line help message',&
            '  -v   Activate verbose logging mode',&
            '  -d   Activate debug logging mode',&
            '  -t   Specify name of thermodynamic database to read (def: thermo.lib)',&
            '  -r   Specify name of transport database to read (def: trans.lib)',&
            '  --compile-thermo  Compile a thermo.inp file into thermo.lib',&
            '  --compile-trans   Compile a trans.inp file into trans.lib',&
            ''
    end subroutine

    function pop_argument(n) result(arg)
        integer, intent(inout) :: n
        character(:), allocatable :: arg
        character(256) :: buffer
        integer :: arglen, stat
        call get_command_argument(n,buffer,arglen,stat)
        if (stat /= 0) call abort('Error reading program argument #'//to_str(n))
        arg = buffer(:arglen)
        n = n+1
    end function

    subroutine run_thermo_problem(prob, thermo, solver, solutions, partials)
        ! Loop over problem state values and solve the thermodynamic equilibrium problems

        ! Arguments
        type(ProblemDB), intent(in) :: prob
        type(ThermoDB), intent(in) :: thermo
        type(EqSolver), intent(out) :: solver
        type(EqSolution), allocatable, intent(out) :: solutions(:, :, :)
        type(EqPartials), allocatable, intent(out) :: partials(:, :, :)

        ! Locals
        type(Mixture) :: reactants, products
        type(EqSolution) :: solution
        type(EqPartials) :: problem_partials
        real(dp) :: state1, state2
        real(dp), allocatable :: weights(:)
        integer :: i, j, k, num_state1, num_state2, num_of
        character(snl), allocatable :: product_names(:), reactant_names(:)

        ! Get the reactants Mixture object
        reactants = Mixture(thermo, input_reactants=prob%reactants, ions=prob%problem%include_ions)
        allocate(weights(reactants%num_species))

        ! Get the products Mixture object
        if (allocated(prob%only)) then
            allocate(product_names(size(prob%only)))
            product_names = prob%only
        else if (allocated(prob%omit)) then
            product_names = reactants%get_products(thermo, prob%omit)
        else
            product_names = reactants%get_products(thermo)
        end if
        products = Mixture(thermo, product_names)

        ! Get the loop sizes
        num_state1 = 1
        select case(prob%problem%type)
            case("tp", "tv")
                if (allocated(prob%problem%T_schedule)) then
                    num_state1 = size(prob%problem%T_schedule%values)
                end if
            case("hp")
                if (allocated(prob%problem%h_schedule)) then
                    num_state1 = size(prob%problem%h_schedule%values)
                end if
            case("sp", "sv")
                if (allocated(prob%problem%s_schedule)) then
                    num_state1 = size(prob%problem%s_schedule%values)
                end if
            case("uv")
                if (allocated(prob%problem%u_schedule)) then
                    num_state1 = size(prob%problem%u_schedule%values)
                end if
        end select

        num_state2 = 1
        select case(prob%problem%type)
            case("tp", "hp", "sp")
                num_state2 = size(prob%problem%p_schedule%values)
            case("tv", "uv", "sv")
                num_state2 = size(prob%problem%v_schedule%values)
        end select

        num_of = 1
        if (allocated(prob%problem%of_schedule)) then
            num_of = size(prob%problem%of_schedule%values)
        end if

        ! Loop over problem state values
        allocate(solutions(num_state1, num_state2, num_of))
        allocate(partials(num_state1, num_state2, num_of))

        ! Initialize the EqSolver and EqSolution objects
        if (allocated(prob%output%trace)) then
            if (prob%output%transport) then
                solver = EqSolver(products, reactants, prob%output%trace, ions=prob%problem%include_ions, &
                                  all_transport=all_transport, insert=prob%insert)
            else
                solver = EqSolver(products, reactants, prob%output%trace, ions=prob%problem%include_ions, &
                                  insert=prob%insert)
            end if
        else
            if (prob%output%transport) then
                solver = EqSolver(products, reactants, ions=prob%problem%include_ions, &
                                  all_transport=all_transport, insert=prob%insert)
            else
                solver = EqSolver(products, reactants, ions=prob%problem%include_ions, &
                                  insert=prob%insert)
            end if
        end if

        ! Initialize the solution object
        solution = EqSolution(solver)

        do i = 1, num_state1 ! State 1 values (temperature, enthalpy, energy, or entropy)
            do j = 1, num_state2 ! State 2 values (pressure or volume)
                state2 = get_state2(prob, j)

                do k = 1, num_of ! Oxidant-to-fuel (or equivalent) ratio values

                    ! Get the reactant weights
                    weights = get_problem_weights(prob, reactants, k)

                    ! Get the fixed state values
                    state1 = get_state1(prob, reactants, weights, i)

                    ! Solve the thermodynamic equilibrium problem
                    ! solution = prev_solution
                    call solver%solve(solution, prob%problem%type, state1, state2, weights, problem_partials)

                    !call output(prob, state1, state2, weights, solution)
                    solutions(i, j, k) = solution
                    partials(i, j, k) = problem_partials

                end do
            end do
        end do

    end subroutine

    subroutine run_rocket_problem(prob, thermo, solver, solutions)
        ! Loop over problem state values and solve rocket problems

        ! Arguments
        type(ProblemDB), intent(in) :: prob
        type(ThermoDB), intent(in) :: thermo
        type(RocketSolver), intent(out) :: solver
        type(RocketSolution), allocatable, intent(out) :: solutions(:, :, :)

        ! Locals
        type(Mixture) :: reactants, products
        type(RocketSolution) :: solution
        real(dp), allocatable :: weights(:)
        integer :: i, j, k, num_pc, num_of
        real(dp) :: pc, hc
        real(dp), allocatable :: tc, tc_est, mdot, ac_at
        real(dp), allocatable :: subar(:), supar(:), pi_p(:)
        logical :: fac, frz, eql, ions, need_hc
        integer :: nfrz
        character(snl), allocatable :: product_names(:), reactant_names(:)

        ! Initialize
        need_hc = .false.

        ! Get the reactants Mixture object
        reactants = Mixture(thermo, input_reactants=prob%reactants, ions=prob%problem%include_ions)
        allocate(weights(reactants%num_species))

        ! Get the products Mixture object
        if (allocated(prob%only)) then
            allocate(product_names(size(prob%only)))
            product_names = prob%only
        else if (allocated(prob%omit)) then
            product_names = reactants%get_products(thermo, prob%omit)
        else
            product_names = reactants%get_products(thermo)
        end if
        products = Mixture(thermo, product_names, ions=prob%problem%include_ions)

        ! Get the loop sizes
        num_pc = 1
        if (allocated(prob%problem%p_schedule)) then
            num_pc = size(prob%problem%p_schedule%values)
        else
            call abort("Chamber pressure not supplied for rocket problem")
        end if

        num_of = 1
        if (allocated(prob%problem%of_schedule)) then
            num_of = size(prob%problem%of_schedule%values)
        end if

        ! Get the rocket variables
        pi_p = prob%problem%pcp_schedule%values
        if (allocated(prob%problem%subar_schedule)) subar = prob%problem%subar_schedule%values
        if (allocated(prob%problem%supar_schedule)) supar = prob%problem%supar_schedule%values
        if (allocated(prob%problem%mdot)) mdot = prob%problem%mdot
        if (allocated(prob%problem%ac_at)) ac_at = prob%problem%ac_at
        if (allocated(prob%problem%tc_est)) tc_est = prob%problem%tc_est
        if (allocated(prob%problem%h_schedule)) then
            hc = prob%problem%h_schedule%values(1) ! Only allow a single value for the rocket problem
        else if (allocated(prob%problem%t_schedule)) then
            tc = prob%problem%t_schedule%values(1) ! Only allow a single value for the rocket problem
        else
            need_hc = .true.
        end if

        ! Parse the boolean arguments
        eql = prob%problem%equilibrium
        frz = prob%problem%frozen
        fac = prob%problem%rkt_finite_area
        nfrz = prob%problem%rkt_nfrozen

        if (.not. eql .and. .not. frz) then
            eql = .true.
        end if

        ! Size the solution output
        if (eql .and. frz) then
            allocate(solutions(num_pc, num_of, 2))
        else
            allocate(solutions(num_pc, num_of, 1))
        end if

        ! Initialize the RocketSolver object
        if (allocated(prob%output%trace)) then
            if (prob%output%transport) then
                solver = RocketSolver(products, reactants, prob%output%trace, ions=prob%problem%include_ions, &
                                      all_transport=all_transport, insert=prob%insert)
            else
                solver = RocketSolver(products, reactants, prob%output%trace, ions=prob%problem%include_ions, insert=prob%insert)
            end if
        else
            if (prob%output%transport) then
                solver = RocketSolver(products, reactants, ions=prob%problem%include_ions, &
                                      all_transport=all_transport, insert=prob%insert)
            else
                solver = RocketSolver(products, reactants, ions=prob%problem%include_ions, insert=prob%insert)
            end if
        end if

        ! Loop over the input parameters
        k = 1
        if (eql) then
            do i = 1, num_pc
                ! Get the chamber pressure value
                pc = get_state2(prob, i)
                do j = 1, num_of
                    weights = get_problem_weights(prob, reactants, j)

                    ! Compute hc from the reactants if needed
                    if (need_hc) hc = compute_reactant_enthalpy(prob, reactants, weights)

                    ! Call the rocket solver
                    solution = solver%solve(weights, pc, pi_p, fac=fac, subar=subar, supar=supar, &
                                            mdot=mdot, ac_at=ac_at, tc_est=tc_est, hc=hc, tc=tc)

                    ! Set the solution
                    solutions(i, j, k) = solution

                end do
            end do
            k = k + 1
        end if

        if (frz) then
            do i = 1, num_pc
                ! Get the chamber pressure value
                pc = get_state2(prob, i)
                do j = 1, num_of
                    weights = get_problem_weights(prob, reactants, j)

                    ! Compute hc from the reactants if needed
                    if (need_hc) hc = compute_reactant_enthalpy(prob, reactants, weights)

                    ! Call the rocket solver
                    solution = solver%solve(weights, pc, pi_p, fac=fac, subar=subar, supar=supar, &
                        mdot=mdot, ac_at=ac_at, n_frz=nfrz, tc_est=tc_est, hc=hc, tc=tc)

                    ! Set the solution
                    solutions(i, j, k) = solution

                end do
            end do
        end if

    end subroutine

    subroutine run_shock_problem(prob, thermo, solver, solutions)
        ! Loop over problem state values and solve shock problems

        ! Arguments
        type(ProblemDB), intent(in) :: prob
        type(ThermoDB), intent(in) :: thermo
        type(ShockSolver), intent(out) :: solver
        type(ShockSolution), allocatable, intent(out) :: solutions(:, :, :)

        ! Locals
        type(Mixture) :: reactants, products
        type(ShockSolution) :: solution
        real(dp) :: T0, P0, u1, mach1
        real(dp), allocatable :: weights(:)
        logical :: incident, input_reflected, frozen, equilibrium, incident_frozen, reflected_frozen, reflected, use_mach
        integer :: i, k, idx_P, idx_T, npts, num_u1, num_P, num_T
        character(snl), allocatable :: product_names(:), reactant_names(:)

        ! Get the reactants Mixture object
        reactants = Mixture(thermo, input_reactants=prob%reactants, ions=prob%problem%include_ions)
        allocate(weights(reactants%num_species))

        ! Get the products Mixture object
        if (allocated(prob%only)) then
            allocate(product_names(size(prob%only)))
            product_names = prob%only
        else if (allocated(prob%omit)) then
            product_names = reactants%get_products(thermo, prob%omit)
        else
            product_names = reactants%get_products(thermo)
        end if
        products = Mixture(thermo, product_names)

        ! Get the problem flags
        if (prob%problem%shk_incident) then
            incident = .true.
        else
            incident = .false.
        end if

        if (prob%problem%shk_reflected) then
            input_reflected = .true.
        else
            input_reflected = .false.
            if (.not. incident) then
                incident = .true.  ! Incident by default
            end if
        end if

        if (prob%problem%equilibrium) then
            equilibrium = .true.
        else
            equilibrium = .false.
        end if

        if (prob%problem%frozen) then
            frozen = .true.
        else
            frozen = .false.
            if (.not. equilibrium) then
                frozen = .true.  ! Frozem by default
            end if
        end if

        ! Total number of shock problem permutations and set solve flags
        npts = 1
        if (equilibrium .and. frozen) then
            if (incident .and. input_reflected) then
                npts = 4
            else
                npts = 2
            end if
        else ! npts = 1
            reflected = input_reflected
            if (equilibrium) then
                if (incident) then
                    incident_frozen = .false.
                else
                    incident_frozen = .true.
                end if
                reflected_frozen = .false.
            else if (frozen) then
                incident_frozen = .true.
                if (reflected) then
                    reflected_frozen = .true.
                else
                    reflected_frozen = .false.
                end if
            else
                ! Weird but true: frozen and equilibrium are both false, only run incident frozen, even if reflected is true
                incident_frozen = .true.
                reflected_frozen = .false.
            end if
        end if

        ! Get the loop sizes
        num_u1 = 1
        if (allocated(prob%problem%u1_schedule)) then
            num_u1 = size(prob%problem%u1_schedule%values)
            use_mach = .false.
        else if (allocated(prob%problem%mach1_schedule)) then
            num_u1 = size(prob%problem%mach1_schedule%values)
            use_mach = .true.
        else
            call abort("Initial velocity or Mach number not supplied for shock problem")
        end if

        num_P = size(prob%problem%p_schedule%values)

        num_T = 1
        if (allocated(prob%problem%t_schedule)) then
            num_T = size(prob%problem%t_schedule%values)
        end if

        if (num_u1 < num_P .or. num_u1 < num_T) then
            call abort("Number of initial velocities must be greater than or equal to number of pressures and temperatures")
        end if

        ! Loop over problem state values
        allocate(solutions(num_u1, 1, npts))

        ! Initialize the ShockSolver and ShockSolution objects
        if (allocated(prob%output%trace)) then
            if (prob%output%transport) then
                solver = ShockSolver(products, reactants, prob%output%trace, ions=prob%problem%include_ions, &
                                     all_transport=all_transport, insert=prob%insert)
            else
                solver = ShockSolver(products, reactants, prob%output%trace, ions=prob%problem%include_ions, &
                                     insert=prob%insert)
            end if
        else
            if (prob%output%transport) then
                solver = ShockSolver(products, reactants, ions=prob%problem%include_ions, &
                                     all_transport=all_transport, insert=prob%insert)
            else
                solver = ShockSolver(products, reactants, ions=prob%problem%include_ions, &
                                     insert=prob%insert)
            end if
        end if

        do k = 1, npts ! Loop over the number of permutations

            ! Set the problem flags if npts > 1 (equilibrium and frozen both true)
            if (npts > 1) then
                if (k == 1) then
                    if (incident) then
                        incident_frozen = .false.
                        reflected_frozen = input_reflected
                        reflected = input_reflected
                    else
                        incident_frozen = .true.
                        reflected_frozen = .true.
                    end if
                else if (k == 2) then
                    if (incident .and. input_reflected) then
                        incident_frozen = .false.
                        reflected = .true.
                        reflected_frozen = .false.
                    else if (incident) then
                        incident_frozen = .true.
                        reflected = .false.
                        reflected_frozen = .false.
                    else
                        incident_frozen = .true.
                        reflected = .true.
                        reflected_frozen = .false.
                    end if
                else if (k == 3) then
                    incident_frozen = .true.
                    reflected = .true.
                    reflected_frozen = .true.
                else  ! k == 4
                    incident_frozen = .true.
                    reflected = .true.
                    reflected_frozen = .false.
                end if
            end if

            do i = 1, num_u1 ! Loop over initial velocity or Mach values
                ! NOTE: Initial velocity or Mach number determines the number of output points. If
                !       there are fewer pressures or temperatures, the last pressure/temperature value is used
                !       for the remaining output points.

                ! Get the initial velocity value
                if (use_mach) then
                    mach1 = prob%problem%mach1_schedule%values(i)
                else
                    u1 = prob%problem%u1_schedule%values(i)
                end if

                ! Get the pressure value
                if (i > num_P) then
                    idx_P = num_P
                else
                    idx_P = i
                end if

                P0 = get_state2(prob, idx_P)

                ! Get the temperature value
                if (i > num_T) then
                    idx_T = num_T
                else
                    idx_T = i
                end if

                if (allocated(prob%problem%t_schedule)) then
                    T0 = prob%problem%t_schedule%values(idx_T)
                else if (allocated(prob%reactants(1)%temperature)) then
                    T0 = prob%reactants(1)%temperature%values(1)
                else
                    call abort("Temperature not supplied for shock problem")
                end if

                ! Get the reactant weights
                weights = get_problem_weights(prob, reactants, 1)

                ! Solve the shock problem
                if (use_mach) then
                    solution = solver%solve(weights, T0, P0, mach1=mach1, reflected=reflected, &
                                            reflected_frozen=reflected_frozen, incident_frozen=incident_frozen)
                else
                    solution = solver%solve(weights, T0, P0, u1=u1, reflected=reflected, &
                                            reflected_frozen=reflected_frozen, incident_frozen=incident_frozen)
                end if
                solutions(i, 1, k) = solution
            end do
        end do

    end subroutine

    subroutine run_detonation_problem(prob, thermo, solver, solutions)
        ! Loop over problem state values and solve detonation problems

        ! Arguments
        type(ProblemDB), intent(in) :: prob
        type(ThermoDB), intent(in) :: thermo
        type(DetonSolver), intent(out) :: solver
        type(DetonSolution), allocatable, intent(out) :: solutions(:, :, :)

        ! Locals
        type(Mixture) :: reactants, products
        type(DetonSolution) :: solution
        real(dp) :: T0, P0
        real(dp), allocatable :: weights(:)
        logical :: frozen, equilibrium
        integer :: i, j, k, num_T, num_P, num_of
        character(snl), allocatable :: product_names(:), reactant_names(:)

        ! Get the reactants Mixture object
        reactants = Mixture(thermo, input_reactants=prob%reactants, ions=prob%problem%include_ions)
        allocate(weights(reactants%num_species))

        ! Get the products Mixture object
        if (allocated(prob%only)) then
            allocate(product_names(size(prob%only)))
            product_names = prob%only
        else if (allocated(prob%omit)) then
            product_names = reactants%get_products(thermo, prob%omit)
        else
            product_names = reactants%get_products(thermo)
        end if
        products = Mixture(thermo, product_names)

        ! Get the problem flags
        if (prob%problem%frozen) then
            frozen = .true.
        else
            frozen = .false.
        end if

        ! Get the problem sizes
        num_T = 1
        if (allocated(prob%problem%T_schedule)) then
            num_T = size(prob%problem%T_schedule%values)
        end if

        if (allocated(prob%problem%P_schedule)) then
            num_P = size(prob%problem%P_schedule%values)
        else
            call abort("Detonation problem requires pressure values for undetonated gas.")
        end if

        num_of = 1
        if (allocated(prob%problem%of_schedule)) then
            num_of = size(prob%problem%of_schedule%values)
        end if

        ! Allocate solutions
        allocate(solutions(num_T, num_P, num_of))

        ! Initial the DetonSolver
        if (allocated(prob%output%trace)) then
            if (prob%output%transport) then
                solver = DetonSolver(products, reactants, prob%output%trace, ions=prob%problem%include_ions, &
                                  all_transport=all_transport, insert=prob%insert)
            else
                solver = DetonSolver(products, reactants, prob%output%trace, ions=prob%problem%include_ions, &
                                  insert=prob%insert)
            end if
        else
            if (prob%output%transport) then
                solver = DetonSolver(products, reactants, ions=prob%problem%include_ions, &
                                  all_transport=all_transport, insert=prob%insert)
            else
                solver = DetonSolver(products, reactants, ions=prob%problem%include_ions, &
                                  insert=prob%insert)
            end if
        end if

        do i = 1, num_T ! Temperature values
            do j = 1, num_P ! Pressure values
                P0 = get_state2(prob, j)

                do k = 1, num_of ! Oxidant-to-fuel (or equivalent) ratio values

                    ! Get the reactant weights
                    weights = get_problem_weights(prob, reactants, k)

                    ! Get the values
                    T0 = get_state1(prob, reactants, weights, i)

                    ! Solve the thermodynamic equilibrium problem
                    solution = solver%solve(weights, T0, P0, frozen)

                    !call output(prob, state1, state2, weights, solution)
                    solutions(i, j, k) = solution

                end do
            end do
        end do

    end subroutine

    subroutine shock_output(ioout, prob, solver, solutions)
        ! Write out an output file for shock problems

        ! Arguments
        integer, intent(in) :: ioout
        type(ProblemDB), intent(in) :: prob
        type(ShockSolver), intent(in) :: solver
        type(ShockSolution), intent(in) :: solutions(:, :, :)

        ! Locals
        integer :: i, j, k, m, npts, idx, last_row_cols, nrows, ncols, num_trace
        integer, parameter :: max_cols = 6
        real(dp) :: trace
        logical :: incident, reflected, equilibrium, frozen, input_reflected
        logical :: write_incd_frz, write_refl_frz, write_incd_eql, write_refl_eql
        logical :: use_mach
        character(snl), allocatable :: trace_names(:)
        logical, allocatable :: is_trace(:)
        character(80) :: eq_fmt
        character(4) :: mass_or_mole
        character(11) :: refl_type, incd_type

        ! Initialization
        m = size(solutions, 1)  ! Number of initial velocities or Mach numbers
        k = 1
        allocate(trace_names(solver%eq_solver%num_products), is_trace(solver%eq_solver%num_products))

        if (allocated(prob%problem%u1_schedule)) then
            use_mach = .false.
        else if (allocated(prob%problem%mach1_schedule)) then
            use_mach = .true.
        end if

        mass_or_mole = "MOLE"
        if (prob%output%mass_fractions) then
            mass_or_mole = "MASS"
        end if

        if (prob%problem%shk_incident) then
            incident = .true.
        else
            incident = .false.
        end if

        if (prob%problem%shk_reflected) then
            input_reflected = .true.
        else
            input_reflected = .false.
            if (.not. incident) then
                incident = .true.  ! Incident by default
            end if
        end if

        if (prob%problem%equilibrium) then
            equilibrium = .true.
        else
            equilibrium = .false.
        end if

        if (prob%problem%frozen) then
            frozen = .true.
        else
            frozen = .false.
            if (.not. equilibrium) then
                frozen = .true.  ! Frozem by default
            end if
        end if

        ! Total number of shock problem permutations and set solve flags
        npts = size(solutions, 3)

        ! Set output flags when npts = 1
        if (npts == 1) then
            if (equilibrium) then
                write_refl_eql = input_reflected
                write_refl_frz = .false.
                refl_type = "EQUILIBRIUM"
                if (incident) then
                    write_incd_frz = .false.
                    write_incd_eql = .true.
                    incd_type = "EQUILIBRIUM"
                end if
            else if (frozen) then
                write_incd_frz = .true.
                write_incd_eql = .false.
                incd_type = "FROZEN     "
                refl_type = "FROZEN     "
                write_refl_eql = .false.
                if (reflected) then
                    write_refl_frz = .true.
                else
                    write_refl_frz = .false.
                end if
            else
                write_incd_frz = .true.
                write_incd_eql = .false.
                incd_type = "FROZEN     "
                write_refl_eql = .false.
                write_refl_frz = .false.
            end if
        end if

        do k = 1, npts

            ! Set the output flags when npts > 1 (equilibrium and frozen both true)
            if (npts > 1) then
                if (k == 1) then
                    if (incident .and. reflected) then
                        write_incd_frz = .false.
                        write_incd_eql = .true.
                        incd_type = "EQUILIBRIUM"
                        write_refl_frz = .true.
                        write_refl_eql = .false.
                        refl_type = "FROZEN     "
                    else if (incident) then
                        write_incd_frz = .false.
                        write_incd_eql = .true.
                        incd_type = "EQUILIBRIUM"
                        write_refl_frz = .false.
                        write_refl_eql = .false.
                    else
                        write_incd_frz = .true.
                        write_incd_eql = .false.
                        incd_type = "FROZEN     "
                        write_refl_frz = .true.
                        write_refl_eql = .false.
                        refl_type = "FROZEN     "
                    end if
                else if (k == 2) then
                    if (incident .and. reflected) then
                        write_incd_frz = .false.
                        write_incd_eql = .false.
                        write_refl_frz = .false.
                        write_refl_eql = .true.
                        refl_type = "EQUILIBRIUM"
                    else if (incident) then
                        write_incd_frz = .true.
                        write_incd_eql = .false.
                        incd_type = "FROZEN     "
                        write_refl_frz = .false.
                        write_refl_eql = .false.
                    else
                        write_incd_frz = .false.
                        write_incd_eql = .false.
                        write_refl_frz = .false.
                        write_refl_eql = .true.
                        refl_type = "EQUILIBRIUM"
                    end if
                else if (k == 3) then
                    write_incd_frz = .true.
                    write_incd_eql = .false.
                    incd_type = "FROZEN     "
                    refl_type = "FROZEN     "
                    write_refl_eql = .false.
                    write_refl_frz = .true.
                else  ! k == 4
                    write_incd_frz = .false.
                    write_incd_eql = .false.
                    refl_type = "EQUILIBRIUM"
                    write_refl_eql = .true.
                    write_refl_frz = .false.
                end if
            end if

            ! -------------------------------------------------------------------
            ! Write results for unshocked gas (equilibrium)
            ! -------------------------------------------------------------------

            if (write_incd_eql .or. write_incd_frz) then
                ! Write the header
                write(ioout, '(A)') "SHOCK WAVE PARAMETERS ASSUMING"
                write(ioout, '(A)') "EQUILIBRIUM COMPOSITION FOR INCIDENT SHOCKED CONDITIONS"
                write(ioout, '(A)') ""
                write(ioout, '(A)') "INITIAL GAS (1)"
                write(ioout, '(A, 14F9.3)') "Mach1            ", (solutions(i, 1, k)%mach(1), i=1,m)
                write(ioout, '(A, 14F9.3)') "u1, m/s          ", (solutions(i, 1, k)%u(1), i=1,m)
                write(ioout, '(A, 14F9.3)') "P, bar           ", (solutions(i, 1, k)%pressure(1), i=1,m)
                write(ioout, '(A, 14F9.3)') "T, K             ", (solutions(i, 1, k)%eq_soln(1)%T, i=1,m)
                write(ioout, '(A, 14E9.3e1)') "rho1, kg/m^3     ", (solutions(i, 1, k)%eq_soln(1)%density, i=1,m)
                write(ioout, '(A, 14F9.2)') "H, kJ/kg         ", (solutions(i, 1, k)%eq_soln(1)%enthalpy, i=1,m)
                write(ioout, '(A, 14F9.2)') "U, kJ/kg         ", (solutions(i, 1, k)%eq_soln(1)%energy, i=1,m)
                write(ioout, '(A, 14F9.1)') "G, kJ/kg         ", (solutions(i, 1, k)%eq_soln(1)%gibbs_energy, i=1,m)
                write(ioout, '(A, 14F9.2)') "S, kJ/(kg-K)     ", (solutions(i, 1, k)%eq_soln(1)%entropy, i=1,m)
                write(ioout, '(A)') ""
                write(ioout, '(A, 14F9.3)') "M, (1/n)         ", (1.0d0/solutions(i, 1, k)%eq_soln(1)%n, i=1,m)
                write(ioout, '(A, 14F9.3)') "Cp, kJ/(kg-K)    ", (solutions(i, 1, k)%eq_soln(1)%cp_eq, i=1,m)
                write(ioout, '(A, 14F9.3)') "Gamma_s          ", (solutions(i, 1, k)%eq_soln(1)%gamma_s, i=1,m)
                write(ioout, '(A, 14F9.3)') "Son. Vel., m/s   ", (solutions(i, 1, k)%v_sonic(1), i=1,m)

                ! -------------------------------------------------------------------
                ! Write results for the incident shock (equilibrium)
                ! -------------------------------------------------------------------
                ! If not converged, write a warning message
                do i = 1, m
                    if (.not. solutions(i, 1, k)%converged) then
                        write(ioout, '(A)') ""
                        if (use_mach) then
                            write(ioout, '(A, 14F9.3)') "WARNING: No convergence for Mach1=",solutions(i, 1, k)%mach(1)
                        else
                            write(ioout, '(A, 14F9.3)') "WARNING: No convergence for u1=",solutions(i, 1, k)%u(1)
                        end if
                        write(ioout, '(A)') "ANSWERS NOT RELIABLE, SOLUTION MAY NOT EXIST"
                    end if
                end do

                ! Write the shocked gas properties
                write(ioout, '(A)') ""
                write(ioout, '(A)') "SHOCKED GAS (2)--INCIDENT--"//incd_type
                write(ioout, '(A, 14F9.3)') "u2, m/s          ", (solutions(i, 1, k)%u(2), i=1,m)
                write(ioout, '(A, 14F9.3)') "P, bar           ", (solutions(i, 1, k)%pressure(2), i=1,m)
                write(ioout, '(A, 14F9.3)') "T, K             ", (solutions(i, 1, k)%eq_soln(2)%T, i=1,m)
                write(ioout, '(A, 14E9.3e1)') "rho, kg/m^3      ", (solutions(i, 1, k)%eq_soln(2)%density, i=1,m)
                write(ioout, '(A, 14F9.2)') "H, kJ/kg         ", (solutions(i, 1, k)%eq_soln(2)%enthalpy, i=1,m)
                write(ioout, '(A, 14F9.2)') "U, kJ/kg         ", (solutions(i, 1, k)%eq_soln(2)%energy, i=1,m)
                write(ioout, '(A, 14F9.1)') "G, kJ/kg         ", (solutions(i, 1, k)%eq_soln(2)%gibbs_energy, i=1,m)
                write(ioout, '(A, 14F9.2)') "S, kJ/(kg-K)     ", (solutions(i, 1, k)%eq_soln(2)%entropy, i=1,m)
                write(ioout, '(A)') ""
                write(ioout, '(A, 14F9.3)') "M, (1/n)         ", (1.0d0/solutions(i, 1, k)%eq_soln(2)%n, i=1,m)
                if (write_incd_eql) then
                    write(ioout, '(A, 14F9.3)') "(dln(V)/dln(P))t ", (solutions(i, 1, k)%eq_partials(2)%dlnV_dlnP, i=1,m)
                    write(ioout, '(A, 14F9.3)') "(dln(V)/dln(T))p ", (solutions(i, 1, k)%eq_partials(2)%dlnV_dlnT, i=1,m)
                end if
                write(ioout, '(A, 14F9.3)') "Cp, kJ/(kg-K)    ", (solutions(i, 1, k)%eq_soln(2)%cp_eq, i=1,m)
                write(ioout, '(A, 14F9.3)') "Gamma_s          ", (solutions(i, 1, k)%eq_partials(2)%gamma_s, i=1,m)
                write(ioout, '(A, 14F9.3)') "Son. Vel., m/s   ", (solutions(i, 1, k)%v_sonic(2), i=1,m)
                write(ioout, '(A)') ""

                ! Write the ratios
                write(ioout, '(A, 14F9.5)') "P2/P1            ", (solutions(i, 1, k)%p21, i=1,m)
                write(ioout, '(A, 14F9.5)') "T2/T1            ", (solutions(i, 1, k)%t21, i=1,m)
                write(ioout, '(A, 14F9.5)') "M2/M1            ", (solutions(i, 1, k)%M21, i=1,m)
                write(ioout, '(A, 14F9.5)') "rho2/rho1        ", (1.0d0/solutions(i, 1, k)%rho12, i=1,m)
                write(ioout, '(A, 14F9.3)') "v2, m/s          ", (solutions(i, 1, k)%v2, i=1,m)
                write(ioout, '(A)') ""

                ! Set the trace output value
                trace = 5.d-9
                if (allocated(prob%output%trace)) trace = prob%output%trace

                ! Get the list of trace species
                num_trace = 0
                do idx = 1, solver%eq_solver%num_products
                    ! Check if this is a trace species
                    is_trace(idx) = .true.
                    do i = 1,m
                        if (is_trace(idx) .eqv. .false.) exit
                        if (prob%output%mass_fractions) then
                            if (solutions(i, 1, k)%eq_soln(2)%mass_fractions(idx) > trace) then
                                is_trace(idx) = .false.
                                exit
                            end if
                        else
                            if (solutions(i, 1, k)%eq_soln(2)%mole_fractions(idx) > trace) then
                                is_trace(idx) = .false.
                                exit
                            end if
                        end if
                    end do
                    if (is_trace(idx)) then
                        num_trace = num_trace + 1
                        trace_names(num_trace) = solver%eq_solver%products%species_names(idx)
                    end if
                end do

                ! Print the mole or mass fractions
                write(ioout, *) ""
                write(ioout, '(A)') mass_or_mole//" FRACTIONS"
                write(ioout, *) ""
                do idx = 1, solver%eq_solver%num_products
                    if (is_trace(idx) .eqv. .false.) then
                        eq_fmt = get_shock_species_format(solutions, idx, 1, k, m, 2, prob%output%mass_fractions, trace)
                        if (prob%output%mass_fractions) then
                            write(ioout, eq_fmt) solver%eq_solver%products%species_names(idx), &
                                (solutions(i, 1, k)%eq_soln(2)%mass_fractions(idx), i=1,m)
                        else
                            write(ioout, eq_fmt) solver%eq_solver%products%species_names(idx), &
                                (solutions(i, 1, k)%eq_soln(2)%mole_fractions(idx), i=1,m)
                        end if
                    end if
                end do
                write(ioout, *) ""

                ! Print the list of products with negligible amounts
                if (write_incd_eql) then
                    write(ioout, '(A)') "PRODUCTS WHICH WERE CONSIDERED BUT WHOSE "// mass_or_mole //" FRACTIONS"
                    write(ioout, '(A, 1PE13.6 ,A)') "WERE LESS THAN", trace, " FOR ALL ASSIGNED CONDITIONS"
                    write(ioout, '(A)') ""
                    ncols = 5
                    nrows = (num_trace - 1)/ncols + 1
                    do i = 1, nrows
                        last_row_cols = merge(mod(num_trace, ncols), ncols, mod(num_trace, ncols) /= 0 .and. i == nrows)
                        write(ioout, '(5A16)') (trace_names((i-1)*ncols + j), j=1,last_row_cols)
                    end do
                    write(ioout, '(A)') ""
                    write(ioout, '(A)') ""
                end if

            end if

            if (write_refl_eql .or. write_refl_frz) then
                write(ioout, '(A)') "SHOCKED GAS (5)--REFLECTED--"//refl_type
                write(ioout, '(A, 14F9.3)') "u5, m/s          ", (solutions(i, 1, k)%u(3), i=1,m)
                write(ioout, '(A, 14F9.3)') "P, bar           ", (solutions(i, 1, k)%pressure(3), i=1,m)
                write(ioout, '(A, 14F9.3)') "T, K             ", (solutions(i, 1, k)%eq_soln(3)%T, i=1,m)
                write(ioout, '(A, 14E9.3e1)') "rho, kg/m^3      ", (solutions(i, 1, k)%eq_soln(3)%density, i=1,m)
                write(ioout, '(A, 14F9.2)') "H, kJ/kg         ", (solutions(i, 1, k)%eq_soln(3)%enthalpy, i=1,m)
                write(ioout, '(A, 14F9.2)') "U, kJ/kg         ", (solutions(i, 1, k)%eq_soln(3)%energy, i=1,m)
                write(ioout, '(A, 14F9.1)') "G, kJ/kg         ", (solutions(i, 1, k)%eq_soln(3)%gibbs_energy, i=1,m)
                write(ioout, '(A, 14F9.2)') "S, kJ/(kg-K)     ", (solutions(i, 1, k)%eq_soln(3)%entropy, i=1,m)
                write(ioout, '(A)') ""
                write(ioout, '(A, 14F9.3)') "M, (1/n)         ", (1.0d0/solutions(i, 1, k)%eq_soln(3)%n, i=1,m)
                if (write_refl_eql) then
                    write(ioout, '(A, 14F9.3)') "(dln(V)/dln(P))t ", (solutions(i, 1, k)%eq_partials(3)%dlnV_dlnP, i=1,m)
                    write(ioout, '(A, 14F9.3)') "(dln(V)/dln(T))p ", (solutions(i, 1, k)%eq_partials(3)%dlnV_dlnT, i=1,m)
                end if
                write(ioout, '(A, 14F9.3)') "Cp, kJ/(kg-K)    ", (solutions(i, 1, k)%eq_soln(3)%cp_eq, i=1,m)
                write(ioout, '(A, 14F9.3)') "Gamma_s          ", (solutions(i, 1, k)%eq_partials(3)%gamma_s, i=1,m)
                write(ioout, '(A, 14F9.3)') "Son. Vel., m/s   ", (solutions(i, 1, k)%v_sonic(3), i=1,m)
                write(ioout, '(A)') ""

                ! Write the ratios
                write(ioout, '(A, 14F9.5)') "P5/P2            ", (solutions(i, 1, k)%p52, i=1,m)
                write(ioout, '(A, 14F9.5)') "T5/T2            ", (solutions(i, 1, k)%t52, i=1,m)
                write(ioout, '(A, 14F9.5)') "M5/M2            ", (solutions(i, 1, k)%M52, i=1,m)
                write(ioout, '(A, 14F9.5)') "rho5/rho2        ", (solutions(i, 1, k)%rho52, i=1,m)
                write(ioout, '(A, 14F9.3)') "u5+v2, m/s       ", (solutions(i, 1, k)%u5_p_v2, i=1,m)
                write(ioout, '(A)') ""

                ! Set the trace output value
                trace = 5.d-9
                if (allocated(prob%output%trace)) trace = prob%output%trace

                ! Get the list of trace species
                num_trace = 0
                do idx = 1, solver%eq_solver%num_products
                    ! Check if this is a trace species
                    is_trace(idx) = .true.
                    do i = 1,m
                        if (is_trace(idx) .eqv. .false.) exit
                        if (prob%output%mass_fractions) then
                            if (solutions(i, 1, k)%eq_soln(3)%mass_fractions(idx) > trace) then
                                is_trace(idx) = .false.
                                exit
                            end if
                        else
                            if (solutions(i, 1, k)%eq_soln(3)%mole_fractions(idx) > trace) then
                                is_trace(idx) = .false.
                                exit
                            end if
                        end if
                    end do
                    if (is_trace(idx)) then
                        num_trace = num_trace + 1
                        trace_names(num_trace) = solver%eq_solver%products%species_names(idx)
                    end if
                end do

                ! Print the mole or mass fractions
                write(ioout, *) ""
                write(ioout, '(A)') mass_or_mole//" FRACTIONS"
                write(ioout, *) ""
                do idx = 1, solver%eq_solver%num_products
                    if (is_trace(idx) .eqv. .false.) then
                        eq_fmt = get_shock_species_format(solutions, idx, 1, k, m, 3, prob%output%mass_fractions, trace)
                        if (prob%output%mass_fractions) then
                            write(ioout, eq_fmt) solver%eq_solver%products%species_names(idx), &
                                (solutions(i, 1, k)%eq_soln(3)%mass_fractions(idx), i=1,m)
                        else
                            write(ioout, eq_fmt) solver%eq_solver%products%species_names(idx), &
                                (solutions(i, 1, k)%eq_soln(3)%mole_fractions(idx), i=1,m)
                        end if
                    end if
                end do
                write(ioout, *) ""

                ! Print the list of products with negligible amounts
                if (write_refl_eql) then
                    write(ioout, '(A)') "PRODUCTS WHICH WERE CONSIDERED BUT WHOSE "// mass_or_mole //" FRACTIONS"
                    write(ioout, '(A, 1PE13.6 ,A)') "WERE LESS THAN", trace, " FOR ALL ASSIGNED CONDITIONS"
                    write(ioout, '(A)') ""
                    ncols = 5
                    nrows = (num_trace - 1)/ncols + 1
                    do i = 1, nrows
                        last_row_cols = merge(mod(num_trace, ncols), ncols, mod(num_trace, ncols) /= 0 .and. i == nrows)
                        write(ioout, '(5A16)') (trace_names((i-1)*ncols + j), j=1,last_row_cols)
                    end do
                    write(ioout, '(A)') ""
                    write(ioout, '(A)') ""
                end if

            end if

        end do

    end subroutine

    subroutine deton_output(ioout, prob, solver, solutions)
        ! Write out an output file for detonation problems

        ! Arguments
        integer, intent(in) :: ioout
        type(ProblemDB), intent(in) :: prob
        type(DetonSolver), intent(in) :: solver
        type(DetonSolution), intent(in) :: solutions(:, :, :)

        ! Locals
        integer :: i, j, k, m, n, idx, last_row_cols, nrows, ncols, num_trace
        integer, parameter :: max_cols = 6
        real(dp) :: trace
        logical :: frozen, transport
        character(snl), allocatable :: trace_names(:)
        logical, allocatable :: is_trace(:)
        character(80) :: eq_fmt
        character(4) :: mass_or_mole

        ! Initialization
        m = size(solutions, 1)
        n = size(solutions, 2)
        allocate(trace_names(solver%eq_solver%num_products), is_trace(solver%eq_solver%num_products))

        mass_or_mole = "MOLE"
        if (prob%output%mass_fractions) then
            mass_or_mole = "MASS"
        end if

        do k = 1, size(solutions, 3)  ! Loop over o/f ratio
            ! Print the thermodynamic properties
            write(ioout, *) ""
            write(ioout, '(A)') ' UNBURNED GAS'
            write(ioout, *) ""

            if (prob%output%siunit) then

                ! Print the values of each variable in the required format
                write(ioout, '(A, 14F9.4)') ' P1, bar          ', ((solutions(i, j, k)%P1,       j=1,n), i=1,m)
                write(ioout, '(A, 14F9.2)') ' T1, K            ', ((solutions(i, j, k)%T1,       j=1,n), i=1,m)
                write(ioout, '(A, 14F9.2)') ' H1, kJ/kg        ', ((solutions(i, j, k)%H1,       j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' M1, (1/n)        ', ((1.0/solutions(i, j, k)%M1,   j=1,n), i=1,m)
                write(ioout, '(A, 14F9.4)') ' Gamma1           ', ((solutions(i, j, k)%gamma1,   j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' Son. Vel.1, m/s  ', ((solutions(i, j, k)%v_sonic1, j=1,n), i=1,m)

            else

                ! Print the values of each variable in the required format
                write(ioout, '(A, 14F9.4)') ' P1, atm          ', ((solutions(i, j, k)%P1/1.01325d0, j=1,n), i=1,m)
                write(ioout, '(A, 14F9.2)') ' T1, K            ', ((solutions(i, j, k)%T1,           j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' H1, cal/g        ', ((solutions(i, j, k)%H1/4.184d0,   j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' M1, (1/n)        ', ((solutions(i, j, k)%M1,           j=1,n), i=1,m)
                write(ioout, '(A, 14F9.4)') ' Gamma1           ', ((solutions(i, j, k)%gamma1,       j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' Son. Vel.1, m/s  ', ((solutions(i, j, k)%v_sonic1,     j=1,n), i=1,m)

            end if

            write(ioout, *) ""
            write(ioout, '(A)') ' BURNED GAS'
            write(ioout, *) ""

            if (prob%output%siunit) then

                ! Print the values of each variable in the required format
                write(ioout, '(A, 14F9.4)') ' P, bar          ', ((solutions(i, j, k)%pressure,             j=1,n), i=1,m)
                write(ioout, '(A, 14F9.2)') ' T, K            ', ((solutions(i, j, k)%eq_soln%T,            j=1,n), i=1,m)
                write(ioout, '(A, 14E9.3e1)') ' Density, kg/m^3   ', ((solutions(i, j, k)%eq_soln%density,  j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' H, kJ/kg        ', ((solutions(i, j, k)%eq_soln%enthalpy,     j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' U, kJ/kg        ', ((solutions(i, j, k)%eq_soln%energy,       j=1,n), i=1,m)
                write(ioout, '(A, 14F9.1)') ' G, kJ/kg        ', ((solutions(i, j, k)%eq_soln%gibbs_energy, j=1,n), i=1,m)
                write(ioout, '(A, 14F9.4)') ' S, kJ/kg-K      ', ((solutions(i, j, k)%eq_soln%entropy,      j=1,n), i=1,m)
                write(ioout, *) ""
                write(ioout, '(A, 14F9.3)') ' M, (1/n)        ', ((1.0/solutions(i, j, k)%eq_soln%n,         j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' (dln(V)/dln(P))t', ((solutions(i, j, k)%eq_partials%dlnV_dlnP, j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' (dln(V)/dln(T))p', ((solutions(i, j, k)%eq_partials%dlnV_dlnT, j=1,n), i=1,m)
                write(ioout, '(A, 14F9.4)') ' Cp, kJ/kg-K     ', ((solutions(i, j, k)%eq_soln%cp_eq,         j=1,n), i=1,m)
                write(ioout, '(A, 14F9.4)') ' Gamma_s         ', ((solutions(i, j, k)%eq_partials%gamma_s,   j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' Son. Vel., m/s  ', &
                    ((sqrt(solutions(i, j, k)%eq_soln%n*R*solutions(i, j, k)%eq_partials%gamma_s*solutions(i, j, k)%eq_soln%T),&
                     j=1,n), i=1,m)

            else

                ! Print the values of each variable in the required format
                write(ioout, '(A, 14F9.4)') ' P, atm          ', ((solutions(i, j, k)%pressure/1.01325d0,   j=1,n), i=1,m)
                write(ioout, '(A, 14F9.2)') ' T, K            ', ((solutions(i, j, k)%eq_soln%T,                    j=1,n), i=1,m)
                write(ioout, '(A, 14E9.3e1)') ' Density, g/cc   ', ((solutions(i, j, k)%eq_soln%density/1.d3,       j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' H, cal/g        ', ((solutions(i, j, k)%eq_soln%enthalpy/4.184d0,     j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' U, cal/g        ', ((solutions(i, j, k)%eq_soln%energy/4.184d0,       j=1,n), i=1,m)
                write(ioout, '(A, 14F9.1)') ' G, cal/g        ', ((solutions(i, j, k)%eq_soln%gibbs_energy/4.184d0, j=1,n), i=1,m)
                write(ioout, '(A, 14F9.4)') ' S, cal/g-K      ', ((solutions(i, j, k)%eq_soln%entropy/4.184d0,      j=1,n), i=1,m)
                write(ioout, *) ""
                write(ioout, '(A, 14F9.3)') ' M, (1/n)        ', ((1.0/solutions(i, j, k)%eq_soln%n,         j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' (dln(V)/dln(P))t', ((solutions(i, j, k)%eq_partials%dlnV_dlnP, j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' (dln(V)/dln(T))p', ((solutions(i, j, k)%eq_partials%dlnV_dlnT, j=1,n), i=1,m)
                write(ioout, '(A, 14F9.4)') ' Cp, cal/g-K     ', ((solutions(i, j, k)%eq_soln%cp_eq/4.184d0, j=1,n), i=1,m)
                write(ioout, '(A, 14F9.4)') ' Gamma_s         ', ((solutions(i, j, k)%eq_partials%gamma_s,   j=1,n), i=1,m)
                write(ioout, '(A, 14F9.3)') ' Son. Vel., m/s  ', &
                    ((sqrt(solutions(i, j, k)%eq_soln%n*R*solutions(i, j, k)%eq_partials%gamma_s*solutions(i, j, k)%eq_soln%T),&
                     j=1,n), i=1,m)

            end if

            ! Write out transport properties
            if (prob%output%transport) then
                write(ioout, *) ""
                write(ioout, '(A)') " TRANSPORT PROPERTIES (GASES ONLY)"
                write(ioout, '(A)') "    CONDUCTIVITY IN UNITS OF MILLICALORIES/(CM)(K)(SEC)"
                write(ioout, *) ""

                ! Viscosity
                write(ioout, '(A, 14F9.4)') " Visc, Millipoise", ((solutions(i, j, k)%eq_soln%viscosity, j=1,n), i=1,m)
                write(ioout, *) ""

                ! Equilibrium properies
                write(ioout, '(A)') " WITH EQUILIBRIUM REACTIONS"
                if (prob%output%siunit) then
                    write(ioout, '(A, 14F9.4)') " Cp, kJ/kg-K     ", ((solutions(i, j, k)%eq_soln%cp_eq, j=1,n), i=1,m)
                    write(ioout, '(A, 14F9.4)') " Conductivity    ", ((solutions(i, j, k)%eq_soln%conductivity_eq, j=1,n), i=1,m)
                else
                    write(ioout, '(A, 14F9.4)') " Cp, cal/g-K     ", ((solutions(i, j, k)%eq_soln%cp_eq/4.184d0, j=1,n), i=1,m)
                    write(ioout, '(A, 14F9.4)') " Conductivity    ", &
                        ((solutions(i, j, k)%eq_soln%conductivity_eq/4.184d0, j=1,n), i=1,m)
                end if
                write(ioout, '(A, 14F9.4)') " Prandtl Number  ", ((solutions(i, j, k)%eq_soln%Pr_eq, j=1,n), i=1,m)
                write(ioout, *) ""

                ! Frozen properties
                write(ioout, '(A)') " WITH FROZEN REACTIONS"
                if (prob%output%siunit) then
                    write(ioout, '(A, 14F9.4)') " Cp, kJ/kg-K     ", ((solutions(i, j, k)%eq_soln%cp_fr, j=1,n), i=1,m)
                    write(ioout, '(A, 14F9.4)') " Conductivity    ", ((solutions(i, j, k)%eq_soln%conductivity_fr, j=1,n), i=1,m)
                else
                    write(ioout, '(A, 14F9.4)') " Cp, cal/g-K     ", ((solutions(i, j, k)%eq_soln%cp_fr/4.184d0, j=1,n), i=1,m)
                    write(ioout, '(A, 14F9.4)') " Conductivity    ", &
                        ((solutions(i, j, k)%eq_soln%conductivity_fr/4.184d0, j=1,n), i=1,m)
                end if
                write(ioout, '(A, 14F9.4)') " Prandtl Number  ", ((solutions(i, j, k)%eq_soln%Pr_fr, j=1,n), i=1,m)
                write(ioout, *) ""

            end if

            ! Write out the detonation parameters
            write(ioout,'(A)') " DETONATION PARAMETERS"
            write(ioout, '(A, 14F9.4)') " P/P1             ", ((solutions(i, j, k)%P_P1,     j=1,n), i=1,m)
            write(ioout, '(A, 14F9.4)') " T/T1             ", ((solutions(i, j, k)%T_T1,     j=1,n), i=1,m)
            write(ioout, '(A, 14F9.4)') " M/M1             ", ((solutions(i, j, k)%M_M1,     j=1,n), i=1,m)
            write(ioout, '(A, 14F9.4)') " rho/rho1         ", ((solutions(i, j, k)%rho_rho1, j=1,n), i=1,m)
            write(ioout, '(A, 14F9.4)') " Det. Mach Number ", ((solutions(i, j, k)%mach,     j=1,n), i=1,m)
            write(ioout, '(A, 14F9.3)') " Det. Vel., m/s   ", ((solutions(i, j, k)%velocity, j=1,n), i=1,m)

            ! Set the trace output value
            trace = 5.d-6
            if (allocated(prob%output%trace)) trace = prob%output%trace

            ! Get the list of trace species
            num_trace = 0
            do idx = 1, solver%eq_solver%num_products
                ! Check if this is a trace species
                is_trace(idx) = .true.
                do i = 1,m
                    if (is_trace(idx) .eqv. .false.) exit
                    do j = 1,n
                        if (solutions(i, j, k)%eq_soln%mole_fractions(idx) > trace) then
                            is_trace(idx) = .false.
                            exit
                        end if
                    end do
                end do
                if (is_trace(idx)) then
                    num_trace = num_trace + 1
                    trace_names(num_trace) = solver%eq_solver%products%species_names(idx)
                end if
            end do

            ! Print the mole or mass fractions
            write(ioout, *) ""
            write(ioout, '(A)') mass_or_mole//" FRACTIONS"
            write(ioout, *) ""
            do idx = 1, solver%eq_solver%num_products
                if (.not. is_trace(idx)) then
                    eq_fmt = get_deton_species_format(solutions, idx, k, m, n, prob%output%mass_fractions, trace)
                    if (prob%output%mass_fractions) then
                        write(ioout, eq_fmt) solver%eq_solver%products%species_names(idx), &
                            ((solutions(i, j, k)%eq_soln%mass_fractions(idx), j=1,n), i=1,m)
                    else
                        write(ioout, eq_fmt) solver%eq_solver%products%species_names(idx), &
                            ((solutions(i, j, k)%eq_soln%mole_fractions(idx), j=1,n), i=1,m)
                    end if
                end if
            end do

            ! Print the list of products with negligible amounts
            write(ioout, '(A)') ""
            write(ioout, '(A)') "PRODUCTS WHICH WERE CONSIDERED BUT WHOSE "// mass_or_mole //" FRACTIONS"
            write(ioout, '(A, 1PE13.6 ,A)') "WERE LESS THAN", trace, " FOR ALL ASSIGNED CONDITIONS"
            write(ioout, '(A)') ""
            ncols = 5
            nrows = (num_trace - 1)/ncols + 1
            do i = 1, nrows
                last_row_cols = merge(mod(num_trace, ncols), ncols, mod(num_trace, ncols) /= 0 .and. i == nrows)
                write(ioout, '(5A16)') (trace_names((i-1)*ncols + j), j=1,last_row_cols)
            end do
            write(ioout, '(A)') ""
            write(ioout, '(A)') ""

        end do

    end subroutine

    function get_state1(prob, reactants, weights, idx) result(state1)
        ! Get the fixed temperature/enthalpy/energy/entropy state from the problem data

        ! Arguments
        type(ProblemDB), intent(in) :: prob
        type(Mixture), intent(in) :: reactants
        real(dp), intent(in) :: weights(:)
        integer, intent(in), optional :: idx

        ! Results
        real(dp) :: state1

        ! Locals
        integer :: idx_
        real(dp), allocatable :: reac_temps(:)

        idx_ = 1
        if (present(idx)) idx_ = idx

        allocate(reac_temps(reactants%num_species))

        ! Get state1: temperature, enthalpy, entropy, or energy
        ! Note: The state value MUST be given in a schedule unless it is an
        ! enthalpy problem, which can be computed from the reactant weights and temperatures
        select case(prob%problem%type)
            case("tp", "tv", "det")
                state1 = convert_units_to_si(prob%problem%T_schedule%values(idx_), prob%problem%T_schedule%units)
            case("hp")
                state1 = compute_reactant_enthalpy(prob, reactants, weights, idx)
            case("sp", "sv")
                if (allocated(prob%problem%s_schedule)) then
                    state1 = prob%problem%s_schedule%values(idx_)
                else  ! Compute the entropy from the reactants
                    call abort("Entropy calculation from reactants not yet supported")
                end if
            case("uv")
                if (allocated(prob%problem%u_schedule)) then
                    state1 = prob%problem%u_schedule%values(idx_)
                else  ! Compute the internal energy from the reactants
                    call abort("Internal energy calculation from reactants not yet supported")
                end if
        end select

    end function

    function get_state2(prob, idx) result(state2)
        ! Get the fixed pressure/volume state from the problem data

        ! Arguments
        type(ProblemDB), intent(in) :: prob
        integer, intent(in), optional :: idx

        ! Results
        real(dp) :: state2

        ! Locals
        integer :: idx_
        real(dp) :: val

        idx_ = 1
        if (present(idx)) idx_ = idx

        select case(prob%problem%type)
            case("tp", "hp", "sp")
                state2 = convert_units_to_si(prob%problem%p_schedule%values(idx_), prob%problem%p_schedule%units)
            case("tv", "uv", "sv")
                select case(prob%problem%v_schedule%units)
                    case('m**3/kg', 'cm**3/g', 'cc/g')
                        state2 = convert_units_to_si(prob%problem%v_schedule%values(idx_), prob%problem%v_schedule%units)
                    case('kg/m**3', 'g/cm**3', 'g/cc')
                        val = convert_units_to_si(prob%problem%v_schedule%values(idx_), prob%problem%v_schedule%units)
                        state2 = 1.0d0/val
                    case default
                        call abort("Volume units not recognized: "//prob%problem%v_schedule%units)
                end select
            case default
                state2 = convert_units_to_si(prob%problem%p_schedule%values(idx_), prob%problem%p_schedule%units)
        end select

    end function

    function compute_reactant_enthalpy(prob, reactants, weights, idx) result(h0)

        ! Arguments
        type(ProblemDB), intent(in) :: prob
        type(Mixture), intent(in) :: reactants
        real(dp), intent(in) :: weights(:)
        integer, intent(in), optional :: idx

        ! Results
        real(dp) :: h0

        ! Locals
        integer :: i, idx_
        real(dp), allocatable :: reac_temps(:)

        idx_ = 1
        if (present(idx)) idx_ = idx

        allocate(reac_temps(reactants%num_species))

        h0 = 0.0d0
        if (allocated(prob%problem%h_schedule)) then
            h0 = prob%problem%h_schedule%values(idx_)
        else  ! Compute the enthalpy from the reactants

            ! Get the reactant temperatures
            do i = 1, reactants%num_species
                reac_temps(i) = 0.0d0
                if (allocated(prob%reactants(1)%temperature)) then
                    reac_temps(i) = convert_units_to_si(prob%reactants(i)%temperature%values(1), &
                                                        prob%reactants(i)%temperature%units)
                end if
            end do

            ! Compute the reactant enthalpy
            h0 = reactants%calc_enthalpy(weights, reac_temps)/R

        end if

    end function

    function get_problem_weights(prob, reactants, idx) result(weights)
        ! Get the reactant weights from the problem data

        ! Arguments
        type(ProblemDB), intent(in) :: prob
        type(Mixture), intent(in) :: reactants
        integer, optional, intent(in) :: idx  ! Index of fuel ratio to use
        real(dp), allocatable :: weights(:)

        ! Locals
        integer :: i, idx_
        real(dp) :: of_ratio, ratio_val
        real(dp), allocatable :: moles(:)
        real(dp), allocatable :: fuel_moles(:), oxidant_moles(:)
        real(dp), allocatable :: fuel_weights(:), oxidant_weights(:)

        allocate(moles(reactants%num_species), weights(reactants%num_species), &
                 fuel_moles(reactants%num_species), oxidant_moles(reactants%num_species), &
                 fuel_weights(reactants%num_species), oxidant_weights(reactants%num_species))

        idx_ = 1
        if (present(idx)) idx_ = idx

        fuel_weights = 0.0d0
        oxidant_weights = 0.0d0
        fuel_moles = 0.0d0
        oxidant_moles = 0.0d0
        weights = 0.0d0

        ! If fuel and oxidant are not specified separately, use the provided weights,
        ! or convert moles to weights
        if (prob%reactants(1)%type == "na") then

            ! If weights are given, return those
            if (prob%reactants(1)%amount%name == "weight_frac") then
                do i = 1, size(prob%reactants)
                    weights(i) = prob%reactants(i)%amount%values(1)
                end do
            end if

            ! If moles are given, convert to weights
            if (prob%reactants(1)%amount%name == "mole_frac") then
                do i = 1, size(prob%reactants)
                    moles(i) = prob%reactants(i)%amount%values(1)
                end do
                weights = reactants%weights_from_moles(moles)
            end if

        ! If fuel and oxidant are specified separately, compute fuel weights and oxidant weights
        else

            if (allocated(prob%reactants(1)%amount) .eqv. .false.) then
                call log_warning("Reactant amounts not specified; assuming 100% for each.")
                do i = 1, size(prob%reactants)
                    if (prob%reactants(i)%type == "fu") then
                        fuel_weights(i) = 1.0
                    else if (prob%reactants(i)%type == "ox") then
                        oxidant_weights(i) = 1.0
                    end if
                end do

            else if (prob%reactants(1)%amount%name == "weight_frac") then
                do i = 1, size(prob%reactants)
                    if (prob%reactants(i)%type == "fu") then
                        fuel_weights(i) = prob%reactants(i)%amount%values(1)
                    else if (prob%reactants(i)%type == "ox") then
                        oxidant_weights(i) = prob%reactants(i)%amount%values(1)
                    end if
                end do

            else
                do i = 1, size(prob%reactants)
                    if (prob%reactants(i)%type == "fu") then
                        fuel_moles(i) = prob%reactants(i)%amount%values(1)
                    else if (prob%reactants(i)%type == "ox") then
                        oxidant_moles(i) = prob%reactants(i)%amount%values(1)
                    end if
                end do
                fuel_weights = reactants%weights_from_moles(fuel_moles)
                oxidant_weights = reactants%weights_from_moles(oxidant_moles)

            end if

            ! If an o/f schedule is set, get the o/f ratio and compute the weights
            if (allocated(prob%problem%of_schedule)) then
                ratio_val = prob%problem%of_schedule%values(idx_)

                ! Get the o/f ratio
                select case(prob%problem%of_schedule%name)
                    case ("f/o", "f/a")
                        of_ratio = 1.0d0/ratio_val

                    case ("%f", "%fuel")
                        of_ratio = (100.0d0-ratio_val)/ratio_val

                    case ("phi")
                        of_ratio = reactants%of_from_phi(oxidant_weights, fuel_weights, ratio_val)

                    case ("r")
                        of_ratio = reactants%of_from_equivalence(oxidant_weights, fuel_weights, ratio_val)

                    case default
                        of_ratio = ratio_val

                end select

                ! Compute the weights from the o/f ratio
                weights = reactants%weights_from_of(oxidant_weights, fuel_weights, of_ratio)

            else

                ! If no o/f schedule is set, add the oxidant and fuel weights
                weights = fuel_weights + oxidant_weights

            end if

        end if

    end function

    subroutine thermo_output(ioout, prob, solver, solutions, partials)
        ! Write out an output file

        ! Arguments
        integer, intent(in) :: ioout
        type(ProblemDB), intent(in) :: prob
        type(EqSolver), intent(in) :: solver
        type(EqSolution), intent(in) :: solutions(:, :, :)
        type(EqPartials), intent(in) :: partials(:, :, :)

        ! Locals
        integer :: i, j, k, idx, m, n, num_trace, nrows, ncols, last_row_cols
        integer, parameter :: max_cols = 6
        real(dp) :: of_ratio, pct_fuel, r_eq, phi_eq
        real(dp) :: trace
        character(snl), allocatable :: trace_names(:)
        character(4) :: mass_or_mole
        logical, allocatable :: is_trace(:)
        character(70) :: eq_fmt

        ! TODO: Print the problem input

        ! Initialization
        m = size(solutions, 1)
        n = size(solutions, 2)
        allocate(trace_names(solver%num_products), is_trace(solver%num_products))

        mass_or_mole = "MOLE"
        if (prob%output%mass_fractions) then
            mass_or_mole = "MASS"
        end if

        do k = 1, size(solutions, 3)  ! Loop over o/f ratio

            ! TODO: Print the reactant amounts and temperatures
            write(ioout, *) ""
            write(ioout, "(A20, A10, A10, A10)") "REACTANT", "MOLES", "ENERGY", "TEMP"
            write(ioout, "(A20, A10, A10, A10)") "", "", "KJ/MOL", "K"
            do i = 1, solver%num_reactants
                write(ioout, "(A20, F10.3, F10.3, F10.2)") solver%reactants%species_names(i), -1.0, -1.0, -1.0
            end do

            ! TODO: Print the o/f ratio and equivalent values
            call compute_fuel_ratios(prob, solver%reactants, k, of_ratio, pct_fuel, r_eq, phi_eq)
            write(ioout, '(A, F10.5, A, F8.5, A, F8.5, A, F8.5)') 'O/F = ', &
                of_ratio, '    % Fuel = ', pct_fuel, '    r, Eq. Ratio = ', r_eq, '    phi, Eq. Ratio = ', phi_eq

            ! Print the thermodynamic properties
            write(ioout, *) ""
            write(ioout, '(A)') ' THERMODYNAMIC PROPERTIES'
            write(ioout, *) ""

            if (prob%output%siunit) then

                ! Print the values of each variable in the required format
                write(ioout, '(A, 18F14.5)') ' P, bar          ', ((solutions(i, j, k)%pressure,     i=1,m), j=1,n)
                write(ioout, '(A, 18F14.2)') ' T, K            ', ((solutions(i, j, k)%T,            i=1,m), j=1,n)
                write(ioout, '(A, 18E14.4e1)') ' Density, kg/m^3 ', ((solutions(i, j, k)%density,    i=1,m), j=1,n)
                write(ioout, '(A, 18F14.3)') ' H, kJ/kg        ', ((solutions(i, j, k)%enthalpy,     i=1,m), j=1,n)
                write(ioout, '(A, 18F14.3)') ' U, kJ/kg        ', ((solutions(i, j, k)%energy,       i=1,m), j=1,n)
                write(ioout, '(A, 18F14.3)') ' G, kJ/kg        ', ((solutions(i, j, k)%gibbs_energy, i=1,m), j=1,n)
                write(ioout, '(A, 18F14.4)') ' S, kJ/kg-K      ', ((solutions(i, j, k)%entropy,      i=1,m), j=1,n)
                write(ioout, *) ""
                write(ioout, '(A, 18F14.5)') ' M, (1/n)        ', ((1.0/solutions(i, j, k)%n,        i=1,m), j=1,n)
                write(ioout, '(A, 18F14.5)') ' (dln(V)/dln(P))t', ((partials(i, j, k)%dlnV_dlnP,     i=1,m), j=1,n)
                write(ioout, '(A, 18F14.5)') ' (dln(V)/dln(T))p', ((partials(i, j, k)%dlnV_dlnT,     i=1,m), j=1,n)
                write(ioout, '(A, 18F14.5)') ' Cp, kJ/kg-K     ', ((solutions(i, j, k)%cp_eq,        i=1,m), j=1,n)
                write(ioout, '(A, 18F14.5)') ' Gamma_s         ', ((partials(i, j, k)%gamma_s,       i=1,m), j=1,n)
                write(ioout, '(A, 18F14.3)') ' Son. Vel., m/s  ', &
                    ((sqrt(solutions(i, j, k)%n * R * partials(i, j, k)%gamma_s * solutions(i, j, k)%T), i=1,m), j=1,n)

            else

                ! Print the values of each variable in the required format
                write(ioout, '(A, 18F14.5)') ' P, atm          ', ((solutions(i, j, k)%pressure/1.01325d0,   i=1,m), j=1,n)
                write(ioout, '(A, 18F14.2)') ' T, K            ', ((solutions(i, j, k)%T,                    i=1,m), j=1,n)
                write(ioout, '(A, 18E14.4e1)') ' Density, g/cc   ', ((solutions(i, j, k)%density/1.d3,       i=1,m), j=1,n)
                write(ioout, '(A, 18F14.3)') ' H, cal/g        ', ((solutions(i, j, k)%enthalpy/4.184d0,     i=1,m), j=1,n)
                write(ioout, '(A, 18F14.3)') ' U, cal/g        ', ((solutions(i, j, k)%energy/4.184d0,       i=1,m), j=1,n)
                write(ioout, '(A, 18F14.3)') ' G, cal/g        ', ((solutions(i, j, k)%gibbs_energy/4.184d0, i=1,m), j=1,n)
                write(ioout, '(A, 18F14.4)') ' S, cal/g-K      ', ((solutions(i, j, k)%entropy/4.184d0,      i=1,m), j=1,n)
                write(ioout, *) ""
                write(ioout, '(A, 18F14.5)') ' M, (1/n)        ', ((1.0/solutions(i, j, k)%n,         i=1,m), j=1,n)
                write(ioout, '(A, 18F14.5)') ' (dln(V)/dln(P))t', ((partials(i, j, k)%dlnV_dlnP,      i=1,m), j=1,n)
                write(ioout, '(A, 18F14.5)') ' (dln(V)/dln(T))p', ((partials(i, j, k)%dlnV_dlnT,      i=1,m), j=1,n)
                write(ioout, '(A, 18F14.5)') ' Cp, cal/g-K     ', ((solutions(i, j, k)%cp_eq/4.184d0, i=1,m), j=1,n)
                write(ioout, '(A, 18F14.5)') ' Gamma_s         ', ((partials(i, j, k)%gamma_s,        i=1,m), j=1,n)
                write(ioout, '(A, 18F14.3)') ' Son. Vel., m/s  ', &
                    ((sqrt(solutions(i, j, k)%n * R * partials(i, j, k)%gamma_s * solutions(i, j, k)%T), i=1,m), j=1,n)

            end if

            ! Write out transport properties
            if (prob%output%transport) then
                write(ioout, *) ""
                write(ioout, '(A)') " TRANSPORT PROPERTIES (GASES ONLY)"
                write(ioout, '(A)') "    CONDUCTIVITY IN UNITS OF MILLICALORIES/(CM)(K)(SEC)"
                write(ioout, *) ""

                ! Viscosity
                write(ioout, '(A, 18F14.4)') " Visc, Millipoise", ((solutions(i, j, k)%viscosity, i=1,m), j=1,n)
                write(ioout, *) ""

                ! Equilibrium properies
                write(ioout, '(A)') " WITH EQUILIBRIUM REACTIONS"
                if (prob%output%siunit) then
                    write(ioout, '(A, 18F14.4)') " Cp, kJ/kg-K     ", ((solutions(i, j, k)%cp_eq, i=1,m), j=1,n)
                    write(ioout, '(A, 18F14.4)') " Conductivity    ", ((solutions(i, j, k)%conductivity_eq, i=1,m), j=1,n)
                else
                    write(ioout, '(A, 18F14.4)') " Cp, cal/g-K     ", ((solutions(i, j, k)%cp_eq/4.184d0, i=1,m), j=1,n)
                    write(ioout, '(A, 18F14.4)') " Conductivity    ", ((solutions(i, j, k)%conductivity_eq/4.184d0, i=1,m), j=1,n)
                end if
                write(ioout, '(A, 18F14.4)') " Prandtl Number  ", ((solutions(i, j, k)%Pr_eq, i=1,m), j=1,n)
                write(ioout, *) ""

                ! Frozen properties
                write(ioout, '(A)') " WITH FROZEN REACTIONS"
                if (prob%output%siunit) then
                    write(ioout, '(A, 18F14.4)') " Cp, kJ/kg-K     ", ((solutions(i, j, k)%cp_fr, i=1,m), j=1,n)
                    write(ioout, '(A, 18F14.4)') " Conductivity    ", ((solutions(i, j, k)%conductivity_fr, i=1,m), j=1,n)
                else
                    write(ioout, '(A, 18F14.4)') " Cp, cal/g-K     ", ((solutions(i, j, k)%cp_fr/4.184d0, i=1,m), j=1,n)
                    write(ioout, '(A, 18F14.4)') " Conductivity    ", ((solutions(i, j, k)%conductivity_fr/4.184d0, i=1,m), j=1,n)
                end if
                write(ioout, '(A, 18F14.4)') " Prandtl Number  ", ((solutions(i, j, k)%Pr_fr, i=1,m), j=1,n)
                write(ioout, *) ""

            end if

            ! Set the trace output value
            trace = 5.d-6
            if (allocated(prob%output%trace)) trace = prob%output%trace

            ! Get the list of trace species
            num_trace = 0
            do idx = 1, solver%num_products
                ! Check if this is a trace species
                is_trace(idx) = .true.
                do i = 1,m
                    if (is_trace(idx) .eqv. .false.) exit
                    do j = 1,n
                        if (solutions(i, j, k)%mole_fractions(idx) > trace) then
                            is_trace(idx) = .false.
                            exit
                        end if
                    end do
                end do
                if (is_trace(idx)) then
                    num_trace = num_trace + 1
                    trace_names(num_trace) = solver%products%species_names(idx)
                end if
            end do

            ! Print the mole or mass fractions
            write(ioout, *) ""
            write(ioout, '(A)') mass_or_mole//" FRACTIONS"
            write(ioout, *) ""
            do idx = 1, solver%num_products
                if (is_trace(idx) .eqv. .false.) then
                    eq_fmt = get_eq_species_format(solutions, idx, k, m, n, prob%output%mass_fractions, trace)
                    if (prob%output%mass_fractions) then
                        write(ioout, eq_fmt) solver%products%species_names(idx), &
                            ((solutions(i, j, k)%mass_fractions(idx), i=1,m), j=1,n)
                    else
                        write(ioout, eq_fmt) solver%products%species_names(idx), &
                            ((solutions(i, j, k)%mole_fractions(idx), i=1,m), j=1,n)
                    end if
                end if
            end do

            ! Print the list of products with negligible amounts
            write(ioout, '(A)') ""
            write(ioout, '(A)') "PRODUCTS WHICH WERE CONSIDERED BUT WHOSE "// mass_or_mole //" FRACTIONS"
            write(ioout, '(A, 1PE13.6 ,A)') "WERE LESS THAN", trace, " FOR ALL ASSIGNED CONDITIONS"
            write(ioout, '(A)') ""
            ncols = 5
            nrows = (num_trace - 1)/ncols + 1
            do i = 1, nrows
                last_row_cols = merge(mod(num_trace, ncols), ncols, mod(num_trace, ncols) /= 0 .and. i == nrows)
                write(ioout, '(5A16)') (trace_names((i-1)*ncols + j), j=1,last_row_cols)
            end do
            write(ioout, '(A)') ""
            write(ioout, '(A)') ""

            ! Clear the trace names for next pass
            do i = 1, size(trace_names)
                trace_names(i) = "               "
            end do

        end do

    end subroutine

    subroutine rocket_output(ioout, prob, solver, solutions)
        ! Write out an output file

        ! Arguments
        integer, intent(in) :: ioout
        type(ProblemDB), intent(in) :: prob
        type(RocketSolver), intent(in) :: solver
        type(RocketSolution), intent(inout) :: solutions(:, :, :)

        ! Locals
        integer :: i, j, k, idx, ii, jj, m, n, num_trace, nrows, ncols, np, ne, nc
        integer :: max_exit, exit_extra, x, y, last_row_cols, nfrz
        integer, parameter :: max_cols = 6
        integer, parameter :: infty_idx = 2
        real(dp) :: of_ratio, pct_fuel, r_eq, phi_eq
        real(dp) :: trace
        character(12) :: fac_print, eql_print
        character(60) :: station_fmt, thermo_fmt, aeat_fmt, cstar_fmt, cf_fmt, isp_fmt
        character(snl), allocatable :: trace_names(:)
        logical, allocatable :: is_trace(:)
        logical :: frozen
        character(4) :: mass_or_mole
        character(60) :: pinj_fmt
        character(60) :: p_fmt
        character(60) :: t_fmt
        character(60) :: h_fmt
        character(60) :: u_fmt
        character(60) :: g_fmt
        character(60) :: s_fmt
        character(60) :: rho_fmt
        character(60) :: vsonic_fmt
        character(10) :: ffmt
        character(90) :: spec_fmt

        ! TODO: Print the problem input again

        ! Initialization
        m = size(solutions, 1)
        n = size(solutions, 2)
        exit_extra = 0
        frozen = .false.
        nfrz = prob%problem%rkt_nfrozen

        allocate(trace_names(solver%eq_solver%num_products),&
                 is_trace(solver%eq_solver%num_products))

        mass_or_mole = "MOLE"
        if (prob%output%mass_fractions) then
            mass_or_mole = "MASS"
        end if

        if (prob%problem%rkt_finite_area) then
            max_exit = 5
            fac_print = "FINITE"
            station_fmt = "(/,21X,'INJECTOR   COMB END     THROAT',5(7X,A4))"
            pinj_fmt = "(' Pinj/P           ',3(F13.4), 5(F13.3))"
            t_fmt =    "(' T, K             ',3(F13.2), 5(F13.2))"
            thermo_fmt = "(1x, A16, 1x,       3(F13.4), 5(F13.4))"
            vsonic_fmt = "(1x, A16, 1x,       3(F13.2), 5(F13.2))"

            ! perf_fmt  = "(1x, A9, 19x,         F13.3,   5(F13.3))"
            aeat_fmt  = "(1x, A9, 19x,         7(F13.4))"
            cstar_fmt = "(1x, A9, 19x,         7(F13.2))"
            cf_fmt    = "(1x, A9, 19x,         7(F13.4))"
            isp_fmt   = "(1x, A9, 19x,         7(F13.2))"

            ffmt = "(5(F13.3))"
            if (prob%output%siunit) then
                p_fmt =    "(' P, bar           ',8(F13.4))"
                rho_fmt =  "(' Density, kg/m^3  ',8(E13.4e1))"
                h_fmt =    "(' H, kJ/kg         ',8(F13.3))"
                u_fmt =    "(' U, kJ/kg         ',8(F13.3))"
                g_fmt =    "(' G, kJ/kg         ',8(F13.2))"
                s_fmt =    "(' S, kJ/kg-K       ',8(F13.3))"
            else
                p_fmt =    "(' P, atm           ',8(F13.4))"
                rho_fmt =  "(' Density, g/cc    ',8(E13.4e1))"
                h_fmt =    "(' H, cal/g         ',8(F13.3))"
                u_fmt =    "(' U, cal/g         ',8(F13.3))"
                g_fmt =    "(' G, cal/g         ',8(F13.2))"
                s_fmt =    "(' S, cal/g-K       ',8(F13.3))"
            end if
        else
            max_exit = 6
            fac_print = "INFINITE"
            station_fmt = "(/,22X,'CHAMBER     THROAT',6(7X,A4))"
            pinj_fmt = "(' Pinf/P           ',2(F13.4), 6(F13.3))"
            t_fmt =    "(' T, K             ',2(F13.2), 6(F13.2))"
            thermo_fmt = "(1x, A16,  1x,      2(F13.4), 6(F13.4))"
            vsonic_fmt = "(1x, A16,  1x,      2(F13.2), 6(F13.2))"

            aeat_fmt  = "(1x, A9, 19x,         F13.4,   6(F13.4))"
            cstar_fmt = "(1x, A9, 19x,         F13.2,   6(F13.2))"
            cf_fmt    = "(1x, A9, 19x,         F13.4,   6(F13.4))"
            isp_fmt   = "(1x, A9, 19x,         F13.2,   6(F13.2))"

            ffmt = "(6(F13.4))"
            if (prob%output%siunit) then
                p_fmt =    "(' P, bar           ',8(F13.4))"
                rho_fmt =  "(' Density, kg/m^3  ',8(E13.4e1))"
                h_fmt =    "(' H, kJ/kg         ',8(F13.3))"
                u_fmt =    "(' U, kJ/kg         ',8(F13.3))"
                g_fmt =    "(' G, kJ/kg         ',8(F13.2))"
                s_fmt =    "(' S, kJ/kg-K       ',8(F13.3))"
            else
                p_fmt =    "(' P, atm           ',8(F13.4))"
                rho_fmt =  "(' Density, g/cc    ',8(E13.4e1))"
                h_fmt =    "(' H, cal/g         ',8(F13.3))"
                u_fmt =    "(' U, cal/g         ',8(F13.3))"
                g_fmt =    "(' G, cal/g         ',8(F13.2))"
                s_fmt =    "(' S, cal/g-K       ',8(F13.3))"
            end if
        end if

        do k = 1, size(solutions, 3)  ! Loop over eql/frz

            if (k == 2) then
                eql_print = "FROZEN"
                frozen = .true.
            else
                if (prob%problem%equilibrium) then
                    eql_print = "EQUILIBRIUM"
                    frozen = .false.
                else if (prob%problem%frozen) then
                    eql_print = "FROZEN"
                    frozen = .true.
                else
                    eql_print = "EQUILIBRIUM"
                    frozen = .false.
                end if
            end if

            ! Write the preamble
            write(ioout, '(A)') "THEORETICAL ROCKET PERFORMANCE ASSUMING "//trim(eql_print)
            write(ioout, '(A)') "COMPOSITION DURING EXPANSION FROM "//trim(fac_print)//" AREA COMBUSTOR"
            write(ioout, *) ""

            do i = 1, m  ! Loop over chamber/injector pressures
                do j = 1, n  ! Loop over o/f ratios

                    ! Get the number of stations ("np"), and the number of exit stations ("ne")
                    np = solutions(i,j,k)%num_pts

                    ! Remove the solution at "inifnity" if this is an FAC problem
                    if (prob%problem%rkt_finite_area) then
                        do idx = 1, np
                            if (idx <= infty_idx) then
                                continue
                            else
                                solutions(i,j,k)%station(idx-1) = solutions(i,j,k)%station(idx)
                                solutions(i,j,k)%eq_soln(idx-1) = solutions(i,j,k)%eq_soln(idx)
                                solutions(i,j,k)%eq_partials(idx-1) = solutions(i,j,k)%eq_partials(idx)
                                solutions(i,j,k)%pressure(idx-1) = solutions(i,j,k)%pressure(idx)
                                solutions(i,j,k)%mach(idx-1) = solutions(i,j,k)%mach(idx)
                                solutions(i,j,k)%gamma_s(idx-1) = solutions(i,j,k)%gamma_s(idx)
                                solutions(i,j,k)%v_sonic(idx-1) = solutions(i,j,k)%v_sonic(idx)
                                solutions(i,j,k)%ae_at(idx-1) = solutions(i,j,k)%ae_at(idx)
                                solutions(i,j,k)%c_star(idx-1) = solutions(i,j,k)%c_star(idx)
                                solutions(i,j,k)%cf(idx-1) = solutions(i,j,k)%cf(idx)
                                solutions(i,j,k)%i_sp(idx-1) = solutions(i,j,k)%i_sp(idx)
                                solutions(i,j,k)%i_vac(idx-1) = solutions(i,j,k)%i_vac(idx)
                            end if
                        end do
                        solutions(i,j,k)%num_pts = solutions(i,j,k)%num_pts - 1
                        np = np - 1
                    end if

                    ! Get the number of exit stations ("ne")
                    ne = np - 2
                    if (prob%problem%rkt_finite_area) then
                        ne = ne - 1
                    end if
                    if (ne > max_exit) then
                        exit_extra = ne - max_exit
                        ne = max_exit
                        np = np - exit_extra
                        if (exit_extra > max_exit) then
                            call log_warning("Too many exit conditions; truncating")
                            exit_extra = max_exit
                        end if
                    end if

                    ! Print the injector/chamber pressure
                    write(ioout, "(A, F6.2, A)") "Pc = ",solutions(i,j,k)%pressure(1)," bar"
                    write(ioout,*) ""

                    ! TODO: Print the area/pressure ratio values

                    ! TODO: Print the reactant amounts and temperatures
                    write(ioout, *) ""
                    write(ioout, "(A20, A10, A10, A10)") "REACTANT", "MOLES", "ENERGY", "TEMP"
                    write(ioout, "(A20, A10, A10, A10)") "", "", "KJ/MOL", "K"
                    do idx = 1, solver%eq_solver%num_reactants
                        write(ioout, "(A20, F10.3, F10.3, F10.2)") solver%eq_solver%reactants%species_names(idx), -1.0, -1.0, -1.0
                    end do
                    write(ioout,*) ""

                    ! Print the o/f ratio and equivalent values
                    call compute_fuel_ratios(prob, solver%eq_solver%reactants, j, of_ratio, pct_fuel, r_eq, phi_eq)
                    write(ioout, '(A, F10.5, A, F8.5, A, F8.5, A, F8.5)') &
                        'O/F = ', of_ratio, '    % Fuel = ', pct_fuel, '    r, Eq. Ratio = ', r_eq, '    phi, Eq. Ratio = ', phi_eq
                    write(ioout,*) ""

                    ! Print the header with the station names
                    write(ioout, station_fmt) (("EXIT"), idx=1,ne)

                    ! Print the rocket properties
                    write(ioout, pinj_fmt) (solutions(i,j,k)%pressure(1)/solutions(i,j,k)%pressure(idx), idx=1,np)
                    if (prob%output%siunit) then
                        write(ioout, p_fmt) (solutions(i,j,k)%pressure(idx), idx=1,np)
                        write(ioout, t_fmt) (solutions(i,j,k)%eq_soln(idx)%T, idx=1,np)
                        write(ioout, rho_fmt) (solutions(i,j,k)%eq_soln(idx)%density, idx=1,np)
                        write(ioout, h_fmt) (solutions(i,j,k)%eq_soln(idx)%enthalpy, idx=1,np)
                        write(ioout, u_fmt) (solutions(i,j,k)%eq_soln(idx)%energy, idx=1,np)
                        write(ioout, g_fmt) (solutions(i,j,k)%eq_soln(idx)%gibbs_energy, idx=1,np)
                        write(ioout, s_fmt) (solutions(i,j,k)%eq_soln(idx)%entropy, idx=1,np)
                        write(ioout, *) ""
                        write(ioout, thermo_fmt) 'M, (1/n)        ', (1.0/solutions(i, j, k)%eq_soln(idx)%n         , idx=1,np)
                        if (frozen .eqv. .false.) then
                            write(ioout, thermo_fmt) '(dln(V)/dln(P))t', (solutions(i, j, k)%eq_partials(idx)%dlnV_dlnP , idx=1,np)
                            write(ioout, thermo_fmt) '(dln(V)/dln(T))p', (solutions(i, j, k)%eq_partials(idx)%dlnV_dlnT , idx=1,np)
                        end if
                        write(ioout, thermo_fmt) 'Cp, kJ/(kg-K)   ', (solutions(i, j, k)%eq_soln(idx)%cp_eq         , idx=1,np)
                        write(ioout, thermo_fmt) 'Gamma_s         ', (solutions(i, j, k)%eq_partials(idx)%gamma_s, idx=1,np)
                        write(ioout, vsonic_fmt) 'Son. Vel., m/s  ', (solutions(i, j, k)%v_sonic(idx), idx=1,np)
                        write(ioout, thermo_fmt) 'Mach            ', (solutions(i, j, k)%mach(idx), idx=1,np)
                    else
                        write(ioout, p_fmt) (solutions(i,j,k)%pressure(idx)/1.01325d0, idx=1,np)
                        write(ioout, t_fmt) (solutions(i,j,k)%eq_soln(idx)%T, idx=1,np)
                        write(ioout, rho_fmt) (solutions(i,j,k)%eq_soln(idx)%density/1.d3, idx=1,np)
                        write(ioout, h_fmt) (solutions(i,j,k)%eq_soln(idx)%enthalpy/4.184d0, idx=1,np)
                        write(ioout, u_fmt) (solutions(i,j,k)%eq_soln(idx)%energy/4.184d0, idx=1,np)
                        write(ioout, g_fmt) (solutions(i,j,k)%eq_soln(idx)%gibbs_energy/4.184d0, idx=1,np)
                        write(ioout, s_fmt) (solutions(i,j,k)%eq_soln(idx)%entropy/4.184d0, idx=1,np)
                        write(ioout, *) ""
                        write(ioout, thermo_fmt) 'M, (1/n)        ', (1.0/solutions(i, j, k)%eq_soln(idx)%n         , idx=1,np)
                        if (frozen .eqv. .false.) then
                            write(ioout, thermo_fmt) '(dln(V)/dln(P))t', (solutions(i, j, k)%eq_partials(idx)%dlnV_dlnP , idx=1,np)
                            write(ioout, thermo_fmt) '(dln(V)/dln(T))p', (solutions(i, j, k)%eq_partials(idx)%dlnV_dlnT , idx=1,np)
                        end if
                        write(ioout, thermo_fmt) 'Cp, cal/(g-K)   ', (solutions(i, j, k)%eq_soln(idx)%cp_eq/4.184d0, idx=1,np)
                        write(ioout, thermo_fmt) 'Gamma_s         ', (solutions(i, j, k)%eq_partials(idx)%gamma_s, idx=1,np)
                        write(ioout, vsonic_fmt) 'Son. Vel., m/s  ', (solutions(i, j, k)%v_sonic(idx), idx=1,np)
                        write(ioout, thermo_fmt) 'Mach            ', (solutions(i, j, k)%mach(idx), idx=1,np)
                    end if

                    ! Write out transport properties
                    if (prob%output%transport) then
                        write(ioout, *) ""
                        write(ioout, '(A)') " TRANSPORT PROPERTIES (GASES ONLY)"
                        write(ioout, '(A)') "    CONDUCTIVITY IN UNITS OF MILLICALORIES/(CM)(K)(SEC)"
                        write(ioout, *) ""

                        ! Viscosity
                        write(ioout, thermo_fmt) " Visc, Millipoise", (solutions(i, j, k)%eq_soln(idx)%viscosity, idx=1,np)
                        write(ioout, *) ""

                        ! Equilibrium properies
                        write(ioout, '(A)') " WITH EQUILIBRIUM REACTIONS"
                        if (prob%output%siunit) then
                            write(ioout, thermo_fmt) " Cp, kJ/(kg-K)   ", (solutions(i, j, k)%eq_soln(idx)%cp_eq, idx=1,np)
                            write(ioout, thermo_fmt) " Conductivity    ", &
                                (solutions(i, j, k)%eq_soln(idx)%conductivity_eq, idx=1,np)
                        else
                            write(ioout, thermo_fmt) " Cp, cal/(g-K)   ", (solutions(i, j, k)%eq_soln(idx)%cp_eq/4.184d0, idx=1,np)
                            write(ioout, thermo_fmt) " Conductivity    ", &
                                (solutions(i, j, k)%eq_soln(idx)%conductivity_eq/4.184d0, idx=1,np)
                        end if
                        write(ioout, thermo_fmt) " Prandtl Number  ", (solutions(i, j, k)%eq_soln(idx)%Pr_eq, idx=1,np)
                        write(ioout, *) ""

                        ! Frozen properties
                        write(ioout, '(A)') " WITH FROZEN REACTIONS"
                        if (prob%output%siunit) then
                            write(ioout, thermo_fmt) " Cp, kJ/(kg-K)   ", (solutions(i, j, k)%eq_soln(idx)%cp_fr, idx=1,np)
                            write(ioout, thermo_fmt) " Conductivity    ", &
                                (solutions(i, j, k)%eq_soln(idx)%conductivity_fr, idx=1,np)
                        else
                            write(ioout, thermo_fmt) " Cp, cal/(g-K)   ", (solutions(i, j, k)%eq_soln(idx)%cp_fr/4.184d0, idx=1,np)
                            write(ioout, thermo_fmt) " Conductivity    ", &
                                (solutions(i, j, k)%eq_soln(idx)%conductivity_fr/4.184d0, idx=1,np)
                        end if
                        write(ioout, thermo_fmt) " Prandtl Number  ", (solutions(i, j, k)%eq_soln(idx)%Pr_fr, idx=1,np)
                        write(ioout, *) ""

                    end if

                    ! Print the performance parameters
                    write(ioout, *) ""
                    write(ioout, '(A)') " PERFORMANCE PARAMETERS"
                    write(ioout, *) ""

                    write(ioout, aeat_fmt) 'Ae/At    ', (solutions(i, j, k)%ae_at(idx), idx=2,np)
                    if (prob%output%siunit) then
                        write(ioout, cstar_fmt) 'C*, m/s  ', (solutions(i, j, k)%c_star(idx), idx=2,np)
                    else
                        write(ioout, cstar_fmt) 'C*, ft/s ', (solutions(i, j, k)%c_star(idx)*3.2808349d0, idx=2,np)
                    end if
                    write(ioout, cf_fmt) 'Cf       ', (solutions(i, j, k)%cf(idx), idx=2,np)
                    if (prob%output%siunit) then
                        write(ioout, isp_fmt) 'Ivac, m/s', (solutions(i, j, k)%i_vac(idx), idx=2,np)
                        write(ioout, isp_fmt) 'Isp, m/s ', (solutions(i, j, k)%i_sp(idx), idx=2,np)
                    else
                        write(ioout, isp_fmt) 'Ivac, lb-s/lb', (solutions(i, j, k)%i_vac(idx)/9.80665d0, idx=2,np)
                        write(ioout, isp_fmt) 'Isp, lb-s/lb ', (solutions(i, j, k)%i_sp(idx)/9.80665d0, idx=2,np)
                    end if

                    ! Set the trace output value
                    trace = 5.d-6
                    if (allocated(prob%output%trace)) trace = prob%output%trace

                    ! Get the list of trace species
                    num_trace = 0
                    do ii = 1, solver%eq_solver%num_products
                        ! Check if this is a trace species
                        is_trace(ii) = .true.
                        if (is_trace(ii) .eqv. .false.) exit
                        do idx = 1,np
                            if (solutions(i, j, k)%eq_soln(idx)%mole_fractions(ii) > trace) then
                                is_trace(ii) = .false.
                                exit
                            end if
                        end do
                        if (is_trace(ii)) then
                            num_trace = num_trace + 1
                            trace_names(num_trace) = solver%eq_solver%products%species_names(ii)
                        end if
                    end do
                    trace_names = trace_names(:num_trace)

                    ! Print the mole or mass fractions
                    write(ioout, *) ""
                    write(ioout, '(A)') mass_or_mole//" FRACTIONS"
                    write(ioout, *) ""
                    do ii = 1, solver%eq_solver%num_products
                        if (is_trace(ii) .eqv. .false.) then
                            if (frozen) then
                                spec_fmt = get_rocket_species_format(solutions, ii, i, j, k, np, &
                                    prob%output%mass_fractions, trace, .true., nfrz)
                            else
                                spec_fmt = get_rocket_species_format(solutions, ii, i, j, k, np, &
                                    prob%output%mass_fractions, trace)
                            end if
                            if (prob%output%mass_fractions) then
                                if (frozen) then
                                    if (solutions(i, j, k)%eq_soln(nfrz)%mass_fractions(ii) > trace) then
                                        write(ioout, spec_fmt) solver%eq_solver%products%species_names(ii), &
                                            solutions(i, j, k)%eq_soln(nfrz)%mass_fractions(ii)
                                    end if
                                else
                                    write(ioout, spec_fmt) solver%eq_solver%products%species_names(ii), &
                                        (solutions(i, j, k)%eq_soln(idx)%mass_fractions(ii), idx=1,np)
                                end if
                            else
                                if (frozen) then
                                    if (solutions(i, j, k)%eq_soln(nfrz)%mole_fractions(ii) > trace) then
                                        write(ioout, spec_fmt) solver%eq_solver%products%species_names(ii), &
                                             solutions(i, j, k)%eq_soln(nfrz)%mole_fractions(ii)
                                    end if
                                else
                                    write(ioout, spec_fmt) solver%eq_solver%products%species_names(ii), &
                                        (solutions(i, j, k)%eq_soln(idx)%mole_fractions(ii), idx=1,np)
                                end if
                            end if
                        end if
                    end do

                    ! Print the list of products with negligible amounts
                    write(ioout, '(A)') ""
                    write(ioout, '(A)') "PRODUCTS WHICH WERE CONSIDERED BUT WHOSE "// mass_or_mole //" FRACTIONS"
                    write(ioout, '(A, 1PE13.6 ,A)') "WERE LESS THAN", trace, " FOR ALL ASSIGNED CONDITIONS"
                    write(ioout, '(A)') ""
                    ncols = 5
                    nrows = (num_trace-1)/ncols + 1
                    do ii = 1, nrows
                        last_row_cols = merge(mod(num_trace, ncols), ncols, mod(num_trace, ncols) /= 0 .and. ii == nrows)
                        write(ioout, '(5A16)') (trace_names((ii-1)*ncols + jj), jj=1,last_row_cols)
                    end do
                    write(ioout, '(A)') ""
                    write(ioout, '(A)') "NOTE. WEIGHT FRACTION OF FUEL IN TOTAL FUELS AND OF OXIDANT IN TOTAL OXIDANTS"
                    write(ioout, '(A)') ""

                    ! Write out the additional exit conditions if any are in overflow
                    if (exit_extra > 0) then
                        ne = exit_extra
                        np = ne + 2
                        if (prob%problem%rkt_finite_area) np = np + 1
                        nc = np - ne

                        ! Print the header with the station names
                        write(ioout, station_fmt) (("EXIT"), idx=1,ne)

                        ! Print the rocket properties
                        x = nc+max_exit+1
                        y = nc+max_exit+ne

                        write(ioout, pinj_fmt, advance="no") (solutions(i,j,k)%pressure(1)/solutions(i,j,k)%pressure(idx), idx=1,nc)
                        write(ioout, ffmt) (solutions(i,j,k)%pressure(1)/solutions(i,j,k)%pressure(idx), idx=x,y)
                        if (prob%output%siunit) then
                            write(ioout, p_fmt, advance="no") (solutions(i,j,k)%pressure(idx), idx=1,nc)
                            write(ioout, ffmt) (solutions(i,j,k)%pressure(idx), idx=x,y)

                            write(ioout, t_fmt, advance="no") (solutions(i,j,k)%eq_soln(idx)%T, idx=1,nc)
                            write(ioout, ffmt) (solutions(i,j,k)%eq_soln(idx)%T, idx=x,y)

                            write(ioout, rho_fmt, advance="no") (solutions(i,j,k)%eq_soln(idx)%density, idx=1,nc)
                            write(ioout, "(6(E13.4e1))") (solutions(i,j,k)%eq_soln(idx)%density, idx=x,y)

                            write(ioout, h_fmt, advance="no") (solutions(i,j,k)%eq_soln(idx)%enthalpy, idx=1,nc)
                            write(ioout, ffmt) (solutions(i,j,k)%eq_soln(idx)%enthalpy, idx=x,y)

                            write(ioout, u_fmt, advance="no") (solutions(i,j,k)%eq_soln(idx)%energy, idx=1,nc)
                            write(ioout, ffmt) (solutions(i,j,k)%eq_soln(idx)%energy, idx=x,y)

                            write(ioout, g_fmt, advance="no") (solutions(i,j,k)%eq_soln(idx)%gibbs_energy, idx=1,nc)
                            write(ioout, "(6(F13.2))") (solutions(i,j,k)%eq_soln(idx)%gibbs_energy, idx=x,y)

                            write(ioout, s_fmt, advance="no") (solutions(i,j,k)%eq_soln(idx)%entropy, idx=1,nc)
                            write(ioout, ffmt) (solutions(i,j,k)%eq_soln(idx)%entropy, idx=x,y)

                            write(ioout, *) ""
                            write(ioout, thermo_fmt, advance="no") 'M, (1/n)        ', &
                                (1.0/solutions(i, j, k)%eq_soln(idx)%n, idx=1,nc)
                            write(ioout, ffmt) (1.0/solutions(i, j, k)%eq_soln(idx)%n, idx=x,y)

                            if (frozen .eqv. .false.) then
                                write(ioout, thermo_fmt, advance="no") '(dln(V)/dln(P))t', &
                                    (solutions(i, j, k)%eq_partials(idx)%dlnV_dlnP, idx=1,nc)
                                write(ioout, ffmt) (solutions(i, j, k)%eq_partials(idx)%dlnV_dlnP, idx=x,y)

                                write(ioout, thermo_fmt, advance="no") '(dln(V)/dln(T))p', &
                                    (solutions(i, j, k)%eq_partials(idx)%dlnV_dlnT, idx=1,nc)
                                write(ioout, ffmt) (solutions(i, j, k)%eq_partials(idx)%dlnV_dlnT, idx=x,y)
                            end if

                            write(ioout, thermo_fmt, advance="no") 'Cp, kJ/(kg-K)   ', &
                                (solutions(i, j, k)%eq_soln(idx)%cp_eq, idx=1,nc)
                            write(ioout, ffmt) (solutions(i, j, k)%eq_soln(idx)%cp_eq, idx=x,y)
                        else
                            write(ioout, p_fmt, advance="no") (solutions(i,j,k)%pressure(idx)/1.01325d0, idx=1,nc)
                            write(ioout, ffmt) (solutions(i,j,k)%pressure(idx), idx=x,y)

                            write(ioout, t_fmt, advance="no") (solutions(i,j,k)%eq_soln(idx)%T, idx=1,nc)
                            write(ioout, ffmt) (solutions(i,j,k)%eq_soln(idx)%T, idx=x,y)

                            write(ioout, rho_fmt, advance="no") (solutions(i,j,k)%eq_soln(idx)%density/1.d3, idx=1,nc)
                            write(ioout, "(6(E13.4e1))") (solutions(i,j,k)%eq_soln(idx)%density, idx=x,y)

                            write(ioout, h_fmt, advance="no") (solutions(i,j,k)%eq_soln(idx)%enthalpy/4.184d0, idx=1,nc)
                            write(ioout, ffmt) (solutions(i,j,k)%eq_soln(idx)%enthalpy, idx=x,y)

                            write(ioout, u_fmt, advance="no") (solutions(i,j,k)%eq_soln(idx)%energy/4.184d0, idx=1,nc)
                            write(ioout, ffmt) (solutions(i,j,k)%eq_soln(idx)%energy, idx=x,y)

                            write(ioout, g_fmt, advance="no") (solutions(i,j,k)%eq_soln(idx)%gibbs_energy/4.184d0, idx=1,nc)
                            write(ioout, "(6(F13.2))") (solutions(i,j,k)%eq_soln(idx)%gibbs_energy, idx=x,y)

                            write(ioout, s_fmt, advance="no") (solutions(i,j,k)%eq_soln(idx)%entropy/4.184d0, idx=1,nc)
                            write(ioout, ffmt) (solutions(i,j,k)%eq_soln(idx)%entropy, idx=x,y)

                            write(ioout, *) ""
                            write(ioout, thermo_fmt, advance="no") 'M, (1/n)        ', &
                                (1.0/solutions(i, j, k)%eq_soln(idx)%n, idx=1,nc)
                            write(ioout, ffmt) (1.0/solutions(i, j, k)%eq_soln(idx)%n, idx=x,y)

                            if (frozen .eqv. .false.) then
                                write(ioout, thermo_fmt, advance="no") '(dln(V)/dln(P))t', &
                                    (solutions(i, j, k)%eq_partials(idx)%dlnV_dlnP, idx=1,nc)
                                write(ioout, ffmt) (solutions(i, j, k)%eq_partials(idx)%dlnV_dlnP, idx=x,y)

                                write(ioout, thermo_fmt, advance="no") '(dln(V)/dln(T))p', &
                                     (solutions(i, j, k)%eq_partials(idx)%dlnV_dlnT, idx=1,nc)
                                write(ioout, ffmt) (solutions(i, j, k)%eq_partials(idx)%dlnV_dlnT, idx=x,y)
                            end if

                            write(ioout, thermo_fmt, advance="no") 'Cp, cal/(g-K)   ', &
                                (solutions(i, j, k)%eq_soln(idx)%cp_eq/4.184d0, idx=1,nc)
                            write(ioout, ffmt) (solutions(i, j, k)%eq_soln(idx)%cp_eq/4.184d0, idx=x,y)
                        end if

                        write(ioout, thermo_fmt, advance="no") 'Gamma_s         ', &
                            (solutions(i, j, k)%eq_partials(idx)%gamma_s, idx=1,nc)
                        write(ioout, ffmt) (solutions(i, j, k)%eq_partials(idx)%gamma_s, idx=x,y)

                        write(ioout, vsonic_fmt, advance="no") 'Son. Vel., m/s  ', &
                            (solutions(i, j, k)%v_sonic(idx), idx=1,nc)
                        write(ioout, ffmt) (solutions(i, j, k)%v_sonic(idx), idx=x,y)

                        write(ioout, thermo_fmt, advance="no") 'Mach            ', &
                            (solutions(i, j, k)%mach(idx), idx=1,nc)
                        write(ioout, ffmt) (solutions(i, j, k)%mach(idx), idx=x,y)

                        ! Write out transport properties
                        if (prob%output%transport) then
                            write(ioout, *) ""
                            write(ioout, '(A)') " TRANSPORT PROPERTIES (GASES ONLY)"
                            write(ioout, '(A)') "    CONDUCTIVITY IN UNITS OF MILLICALORIES/(CM)(K)(SEC)"
                            write(ioout, *) ""

                            ! Viscosity
                            write(ioout, '(A, F13.3)', advance="no") " Visc, Millipoise", &
                                (solutions(i, j, k)%eq_soln(idx)%viscosity, idx=1,nc)
                            write(ioout, '(F13.3)') (solutions(i, j, k)%eq_soln(idx)%viscosity, idx=x,y)
                            write(ioout, *) ""

                            ! Equilibrium properies
                            write(ioout, '(A)') " WITH EQUILIBRIUM REACTIONS"
                            if (prob%output%siunit) then
                                write(ioout, '(A, F13.4)', advance="no") " Cp, kJ/(kg-K)   ", &
                                    (solutions(i, j, k)%eq_soln(idx)%cp_eq, idx=1,nc)
                                write(ioout, '(F13.4)') (solutions(i, j, k)%eq_soln(idx)%cp_eq, idx=x,y)

                                write(ioout, '(A, F13.3)', advance="no") " Conductivity    ", &
                                    (solutions(i, j, k)%eq_soln(idx)%conductivity_eq, idx=1,nc)
                                write(ioout, '(F13.3)') (solutions(i, j, k)%eq_soln(idx)%conductivity_eq, idx=x,y)
                            else
                                write(ioout, '(A, F13.4)', advance="no") " Cp, cal/(g-K)   ", &
                                    (solutions(i, j, k)%eq_soln(idx)%cp_eq/4.184d0, idx=1,nc)
                                write(ioout, '(F13.4)') (solutions(i, j, k)%eq_soln(idx)%cp_eq/4.184d0, idx=x,y)

                                write(ioout, '(A, F13.3)', advance="no") " Conductivity    ", &
                                    (solutions(i, j, k)%eq_soln(idx)%conductivity_eq/4.184d0, idx=1,nc)
                                write(ioout, '(F13.3)') (solutions(i, j, k)%eq_soln(idx)%conductivity_eq/4.184d0, idx=x,y)
                            end if
                            write(ioout, '(A, F13.4)', advance="no") " Prandtl Number  ", &
                                (solutions(i, j, k)%eq_soln(idx)%Pr_eq, idx=1,nc)
                            write(ioout, '(F13.4)') (solutions(i, j, k)%eq_soln(idx)%Pr_eq, idx=x,y)
                            write(ioout, *) ""

                            ! Frozen properties
                            write(ioout, '(A)') " WITH FROZEN REACTIONS"
                            if (prob%output%siunit) then
                                write(ioout, '(A, F13.4)', advance="no") " Cp, kJ/(kg-K)   ", &
                                    (solutions(i, j, k)%eq_soln(idx)%cp_fr, idx=1,nc)
                                write(ioout, '(F13.4)') (solutions(i, j, k)%eq_soln(idx)%cp_fr, idx=x,y)

                                write(ioout, '(A, F13.3)', advance="no") " Conductivity    ", &
                                    (solutions(i, j, k)%eq_soln(idx)%conductivity_fr, idx=1,nc)
                                write(ioout, '(F13.3)') (solutions(i, j, k)%eq_soln(idx)%conductivity_fr, idx=x,y)
                            else
                                write(ioout, '(A, F13.4)', advance="no") " Cp, cal/(g-K)   ", &
                                    (solutions(i, j, k)%eq_soln(idx)%cp_fr/4.184d0, idx=1,nc)
                                write(ioout, '(F13.4)') (solutions(i, j, k)%eq_soln(idx)%cp_fr/4.184d0, idx=x,y)

                                write(ioout, '(A, F13.3)', advance="no") " Conductivity    ", &
                                    (solutions(i, j, k)%eq_soln(idx)%conductivity_fr/4.184d0, idx=1,nc)
                                write(ioout, '(F13.3)') (solutions(i, j, k)%eq_soln(idx)%conductivity_fr/4.184d0, idx=x,y)
                            end if
                            write(ioout, '(A, F13.4)', advance="no") " Prandtl Number  ", &
                                (solutions(i, j, k)%eq_soln(idx)%Pr_fr, idx=1,nc)
                            write(ioout, '(F13.4)') (solutions(i, j, k)%eq_soln(idx)%Pr_fr, idx=x,y)
                            write(ioout, *) ""

                        end if

                        ! Print the performance parameters
                        write(ioout, *) ""
                        write(ioout, '(A)') " PERFORMANCE PARAMETERS"
                        write(ioout, *) ""
                        write(ioout, aeat_fmt, advance="no") 'Ae/At    ', (solutions(i, j, k)%ae_at(idx), idx=2,nc)
                        write(ioout, "(5(F13.4))") (solutions(i, j, k)%ae_at(idx), idx=x,y)

                        write(ioout, cstar_fmt, advance="no") 'C*, m/s  ', (solutions(i, j, k)%c_star(idx), idx=2,nc)
                        write(ioout, "(5(F13.2))") (solutions(i, j, k)%c_star(idx), idx=x,y)

                        write(ioout, cf_fmt, advance="no") 'Cf       ', (solutions(i, j, k)%cf(idx), idx=2,nc)
                        write(ioout, "(5(F13.4))") (solutions(i, j, k)%cf(idx), idx=x,y)

                        write(ioout, isp_fmt, advance="no") 'Ivac, m/s', (solutions(i, j, k)%i_vac(idx), idx=2,nc)
                        write(ioout, "(5(F13.2))") (solutions(i, j, k)%i_vac(idx), idx=x,y)

                        write(ioout, isp_fmt, advance="no") 'Isp, m/s ', (solutions(i, j, k)%i_sp(idx), idx=2,nc)
                        write(ioout, "(5(F13.2))") (solutions(i, j, k)%i_sp(idx), idx=x,y)

                        ! Get the list of trace species
                        num_trace = 0
                        do ii = 1, solver%eq_solver%num_products
                            ! Check if this is a trace species
                            is_trace(ii) = .true.
                            if (is_trace(ii) .eqv. .false.) exit
                            do idx = 1,nc  ! Check the chamber and throat conditions
                                if (solutions(i, j, k)%eq_soln(idx)%mole_fractions(ii) > trace) then
                                    is_trace(ii) = .false.
                                    exit
                                end if
                            end do
                            do idx = x,y  ! Check the additional exit conditions
                                if (solutions(i, j, k)%eq_soln(idx)%mole_fractions(ii) > trace) then
                                    is_trace(ii) = .false.
                                    exit
                                end if
                            end do
                            if (is_trace(ii)) then
                                num_trace = num_trace + 1
                                trace_names(num_trace) = solver%eq_solver%products%species_names(ii)
                            end if
                        end do
                        trace_names = trace_names(:num_trace)

                        ! Print the mole or mass fractions
                        write(ioout, *) ""
                        write(ioout, '(A)') mass_or_mole//" FRACTIONS"
                        write(ioout, *) ""
                        do ii = 1, solver%eq_solver%num_products
                            if (is_trace(ii) .eqv. .false.) then
                                if (prob%output%mass_fractions) then
                                    if (frozen) then
                                        if (solutions(i, j, k)%eq_soln(nfrz)%mass_fractions(ii) > trace) then
                                            write(ioout, '(1x, A15, 2x, F13.5)') &
                                                solver%eq_solver%products%species_names(ii), &
                                                solutions(i, j, k)%eq_soln(nfrz)%mass_fractions(ii)
                                        end if
                                    else
                                        write(ioout, '(1x, A15, 2x, 3(F13.5))', advance="no") &
                                            solver%eq_solver%products%species_names(ii), &
                                            (solutions(i, j, k)%eq_soln(idx)%mass_fractions(ii), idx=1,nc)
                                        write(ioout, '(3(F13.5)))') (solutions(i, j, k)%eq_soln(idx)%mass_fractions(ii), idx=x,y)
                                    end if
                                else
                                    if (frozen) then
                                        if (solutions(i, j, k)%eq_soln(nfrz)%mole_fractions(ii) > trace) then
                                            write(ioout, '(1x, A15, 2x, F13.5)') &
                                                solver%eq_solver%products%species_names(ii), &
                                                solutions(i, j, k)%eq_soln(nfrz)%mole_fractions(ii)
                                        end if
                                    else
                                        write(ioout, '(1x, A15, 2x, 3(F13.5))', advance="no") &
                                            solver%eq_solver%products%species_names(ii), &
                                            (solutions(i, j, k)%eq_soln(idx)%mole_fractions(ii), idx=1,nc)
                                        write(ioout, '(3(F13.5)))') (solutions(i, j, k)%eq_soln(idx)%mole_fractions(ii), idx=x,y)
                                    end if
                                end if
                            end if
                        end do

                        ! Print the list of products with negligible amounts
                        write(ioout, '(A)') ""
                        write(ioout, '(A)') "PRODUCTS WHICH WERE CONSIDERED BUT WHOSE "// mass_or_mole //" FRACTIONS"
                        write(ioout, '(A, 1PE13.6 ,A)') "WERE LESS THAN", trace, " FOR ALL ASSIGNED CONDITIONS"
                        write(ioout, '(A)') ""
                        ncols = 5
                        nrows = (num_trace - 1)/ncols + 1
                        do ii = 1, nrows
                            last_row_cols = merge(mod(num_trace, ncols), ncols, mod(num_trace, ncols) /= 0 .and. ii == nrows)
                            write(ioout, '(5A16)') (trace_names((ii-1)*ncols + jj), jj=1,last_row_cols)
                        end do
                        write(ioout, '(A)') ""
                        write(ioout, '(A)') "NOTE. WEIGHT FRACTION OF FUEL IN TOTAL FUELS AND OF OXIDANT IN TOTAL OXIDANTS"
                        write(ioout, '(A)') ""

                    end if

                end do
            end do
        end do

    end subroutine

    subroutine compute_fuel_ratios(prob, reactants, idx, of_ratio, pct_fuel, r_eq, phi_eq)
        ! Compute all of the fuel ratio values (o/f ratio, % fuel, r eq. ratio, phi eq. ratio) at once

        ! Arguments
        type(ProblemDB), intent(in) :: prob
        type(Mixture), intent(in) :: reactants
        integer, intent(in) :: idx
        real(dp), intent(inout) :: of_ratio
        real(dp), intent(inout) :: pct_fuel
        real(dp), intent(inout) :: r_eq
        real(dp), intent(inout) :: phi_eq

        ! Locals
        integer :: i
        real(dp) :: ratio_val
        real(dp), allocatable :: fuel_moles(:), oxidant_moles(:)
        real(dp), allocatable :: fuel_weights(:), oxidant_weights(:)

        allocate(fuel_moles(reactants%num_species), oxidant_moles(reactants%num_species), &
                 fuel_weights(reactants%num_species), oxidant_weights(reactants%num_species))

        fuel_weights = 0.0d0
        oxidant_weights = 0.0d0
        fuel_moles = 0.0d0
        oxidant_moles = 0.0d0

        ! Set default values
        of_ratio = 0.0d0
        pct_fuel = 100.0d0
        r_eq = 0.0d0
        phi_eq = 0.0d0

        ! Ratio values are meaningless if fuel and oxidant are not specified separately
        if (prob%reactants(1)%type == "name") then
            return

        ! If fuel and oxidant are specified separately, compute fuel weights and oxidant weights
        else

            if (allocated(prob%reactants(1)%amount) .eqv. .false.) then
                do i = 1, size(prob%reactants)
                    if (prob%reactants(i)%type == "fu") then
                        fuel_weights(i) = 1.0
                    else if (prob%reactants(i)%type == "ox") then
                        oxidant_weights(i) = 1.0
                    end if
                end do

            else if (prob%reactants(1)%amount%name == "weight_frac") then
                do i = 1, size(prob%reactants)
                    if (prob%reactants(i)%type == "fu") then
                        fuel_weights(i) = prob%reactants(i)%amount%values(1)
                    else if (prob%reactants(i)%type == "ox") then
                        oxidant_weights(i) = prob%reactants(i)%amount%values(1)
                    end if
                end do

            else
                do i = 1, size(prob%reactants)
                    if (prob%reactants(i)%type == "fu") then
                        fuel_moles(i) = prob%reactants(i)%amount%values(1)
                    else if (prob%reactants(i)%type == "ox") then
                        oxidant_moles(i) = prob%reactants(i)%amount%values(1)
                    end if
                end do
                fuel_weights = reactants%weights_from_moles(fuel_moles)
                oxidant_weights = reactants%weights_from_moles(oxidant_moles)

            end if

            ! If an o/f schedule is set, compute the fuel ratios
            if (allocated(prob%problem%of_schedule)) then
                ratio_val = prob%problem%of_schedule%values(idx)

                select case(prob%problem%of_schedule%name)
                    case ("o/f")
                        of_ratio = ratio_val
                        pct_fuel = 100.0d0/(1.0d0 + ratio_val)
                        r_eq = reactants%equivalence_from_of(oxidant_weights, fuel_weights, of_ratio)
                        phi_eq = reactants%phi_from_of(oxidant_weights, fuel_weights, of_ratio)

                    case ("f/o", "f/a")
                        of_ratio = 1.0d0/ratio_val
                        pct_fuel = 100.0d0/(1.0d0 + of_ratio)
                        r_eq = reactants%equivalence_from_of(oxidant_weights, fuel_weights, of_ratio)
                        phi_eq = reactants%phi_from_of(oxidant_weights, fuel_weights, of_ratio)

                    case ("%f")
                        pct_fuel = ratio_val
                        of_ratio = (100.0d0-ratio_val)/ratio_val
                        r_eq = reactants%equivalence_from_of(oxidant_weights, fuel_weights, of_ratio)
                        phi_eq = reactants%phi_from_of(oxidant_weights, fuel_weights, of_ratio)

                    case ("phi")
                        phi_eq = ratio_val
                        of_ratio = reactants%of_from_phi(oxidant_weights, fuel_weights, ratio_val)
                        pct_fuel = 100.0d0/(1.0d0 + of_ratio)
                        r_eq = reactants%equivalence_from_of(oxidant_weights, fuel_weights, of_ratio)

                    case ("r")
                        r_eq = ratio_val
                        of_ratio = reactants%of_from_equivalence(oxidant_weights, fuel_weights, ratio_val)
                        pct_fuel = 100.0d0/(1.0d0 + of_ratio)
                        phi_eq = reactants%phi_from_of(oxidant_weights, fuel_weights, of_ratio)

                end select

            end if

        end if

    end subroutine

    function get_eq_species_format(solutions, idx, k, m, n, mass_frac, trace) result(eq_fmt)
        ! Write out the line formatting string for mole or mass fraction output

        ! Arguments
        type(EqSolution), intent(in) :: solutions(:, :, :)
        integer, intent(in) :: idx  ! Value of index 3 of solutions
        integer, intent(in) :: k  ! o/f ratio index
        integer, intent(in) :: m  ! Size of index 1 of solutions
        integer, intent(in) :: n  ! Size of index 2 of solutions
        logical, intent(in) :: mass_frac  ! If true, use mass fractions; if false, use mole fractions
        real(dp), intent(in) :: trace  ! Threshold amount

        ! Result
        character(70) :: eq_fmt

        ! Locals
        integer :: i, j, idx1, idx2
        real(dp) :: amount

        ! Set the format
        idx1 = 1
        idx2 = 9
        eq_fmt(1:) = " "
        eq_fmt(idx1:idx2) = '(A16, 1x'
        do j = 1, n
            do i = 1, m
                if (mass_frac) then
                    amount = solutions(i, j, k)%mass_fractions(idx)
                else
                    amount = solutions(i, j, k)%mole_fractions(idx)
                end if

                idx1 = idx2

                if (amount > 1.d-3 .or. amount < 1.d-20) then
                    idx2 = idx1 + 6
                    eq_fmt(idx1:idx2) = ', F9.6'
                elseif (amount < 1.d-9) then
                    idx2 = idx1 + 9
                    eq_fmt(idx1:idx2) = ', ES9.2E2'
                else
                    idx2 = idx1 + 9
                    eq_fmt(idx1:idx2) = ', ES9.3E1'
                end if
            end do
        end do
        eq_fmt(idx2:idx2+1) = ')'

    end function

    function get_rocket_species_format(solutions, idx, i, j, k, np, mass_frac, trace, frozen, nfrz) result(spec_fmt)
        ! Write out the line formatting string for mole or mass fraction output

        ! Arguments
        type(RocketSolution), intent(in) :: solutions(:, :, :)
        integer, intent(in) :: idx ! Species index
        integer, intent(in) :: i  ! index 1 of solutions
        integer, intent(in) :: j  ! index 2 of solutions
        integer, intent(in) :: k  ! o/f ratio index
        integer, intent(in) :: np  ! Number of product species to write out on this line
        logical, intent(in) :: mass_frac  ! If true, use mass fractions; if false, use mole fractions
        real(dp), intent(in) :: trace  ! Threshold amount
        logical, intent(in), optional :: frozen  ! If true, use frozen values; if false, use equilibrium values
        integer, intent(in), optional :: nfrz  ! Index of frozen solution

        ! Result
        character(90) :: spec_fmt

        ! Locals
        integer :: idx1, idx2, ii, nfrz_
        real(dp) :: amount
        logical :: frozen_

        ! Defaults
        frozen_ = .false.
        if (present(frozen)) then
            frozen_ = frozen
            if (present(nfrz)) then
                nfrz_ = nfrz
            else
                nfrz_ = 2
            end if
        end if

        ! Set the format
        idx1 = 1
        idx2 = 13
        spec_fmt(1:) = " "
        spec_fmt(idx1:idx2) = '(1x, A15, 2x'
        if (frozen_) then
            if (mass_frac) then
                amount = solutions(i, j, k)%eq_soln(nfrz_)%mass_fractions(idx)
            else
                amount = solutions(i, j, k)%eq_soln(nfrz_)%mole_fractions(idx)
            end if

            idx1 = idx2

            if (amount > 1.d-3 .or. amount < 1.d-20) then
                idx2 = idx1 + 6
                spec_fmt(idx1:idx2) = ', F9.6'
            else if (amount < 1.d-9) then
                idx2 = idx1 + 9
                spec_fmt(idx1:idx2) = ', ES9.2E2'
            else
                idx2 = idx1 + 9
                spec_fmt(idx1:idx2) = ', ES9.3E1'
            end if
        else
            do ii = 1, np
                if (mass_frac) then
                    amount = solutions(i, j, k)%eq_soln(ii)%mass_fractions(idx)
                else
                    amount = solutions(i, j, k)%eq_soln(ii)%mole_fractions(idx)
                end if

                idx1 = idx2

                if (amount > 1.d-3 .or. amount < 1.d-20) then
                    idx2 = idx1 + 6
                    spec_fmt(idx1:idx2) = ', F9.6'
                elseif (amount < 1.d-9) then
                    idx2 = idx1 + 9
                    spec_fmt(idx1:idx2) = ', ES9.2E2'
                else
                    idx2 = idx1 + 9
                    spec_fmt(idx1:idx2) = ', ES9.3E1'
                end if
            end do
        end if
        spec_fmt(idx2:idx2+1) = ')'

    end function

    function get_shock_species_format(solutions, idx, j, k, m, n, mass_frac, trace) result(spec_fmt)
        ! Write out the line formatting string for mole or mass fraction output

        ! Arguments
        type(ShockSolution), intent(in) :: solutions(:, :, :)
        integer, intent(in) :: idx ! Species index
        integer, intent(in) :: j  ! Should be fixed = 1
        integer, intent(in) :: k  ! eql/frozen index
        integer, intent(in) :: m  ! Number of initial shock velocities or Mach numbers
        integer, intent(in) :: n  ! Shock station index
        logical, intent(in) :: mass_frac  ! If true, use mass fractions; if false, use mole fractions
        real(dp), intent(in) :: trace  ! Threshold amount

        ! Result
        character(80) :: spec_fmt

        ! Locals
        integer :: i, idx1, idx2
        real(dp) :: amount

        ! Set the format
        idx1 = 1
        idx2 = 13
        spec_fmt(1:) = " "
        spec_fmt(idx1:idx2) = '(1x, A15, 2x'
        do i = 1, m
            if (mass_frac) then
                amount = solutions(i, j, k)%eq_soln(n)%mass_fractions(idx)
            else
                amount = solutions(i, j, k)%eq_soln(n)%mole_fractions(idx)
            end if

            idx1 = idx2

            if (amount > 1.d-3 .or. amount < 1.d-20) then
                idx2 = idx1 + 6
                spec_fmt(idx1:idx2) = ', F9.6'
            elseif (amount < 1.d-9) then
                idx2 = idx1 + 9
                spec_fmt(idx1:idx2) = ', ES9.2E2'
            else
                idx2 = idx1 + 9
                spec_fmt(idx1:idx2) = ', ES9.3E1'
            end if
        end do
        spec_fmt(idx2:idx2+1) = ')'

    end function

    function get_deton_species_format(solutions, idx, k, m, n, mass_frac, trace) result(eq_fmt)
        ! Write out the line formatting string for mole or mass fraction output

        ! Arguments
        type(DetonSolution), intent(in) :: solutions(:, :, :)
        integer, intent(in) :: idx  ! Value of index 3 of solutions
        integer, intent(in) :: k  ! o/f ratio index
        integer, intent(in) :: m  ! Size of index 1 of solutions
        integer, intent(in) :: n  ! Size of index 2 of solutions
        logical, intent(in) :: mass_frac  ! If true, use mass fractions; if false, use mole fractions
        real(dp), intent(in) :: trace  ! Threshold amount

        ! Result
        character(70) :: eq_fmt

        ! Locals
        integer :: i, j, idx1, idx2
        real(dp) :: amount

        ! Set the format
        idx1 = 1
        idx2 = 9
        eq_fmt(1:) = " "
        eq_fmt(idx1:idx2) = '(A16, 1x'
        do j = 1, n
            do i = 1, m
                if (mass_frac) then
                    amount = solutions(i, j, k)%eq_soln%mass_fractions(idx)
                else
                    amount = solutions(i, j, k)%eq_soln%mole_fractions(idx)
                end if

                idx1 = idx2

                if (amount > 1.d-3 .or. amount < 1.d-20) then
                    idx2 = idx1 + 6
                    eq_fmt(idx1:idx2) = ', F9.6'
                elseif (amount < 1.d-9) then
                    idx2 = idx1 + 9
                    eq_fmt(idx1:idx2) = ', ES9.2E2'
                else
                    idx2 = idx1 + 9
                    eq_fmt(idx1:idx2) = ', ES9.3E1'
                end if
            end do
        end do
        eq_fmt(idx2:idx2+1) = ')'

    end function

end program
