module cea_rocket
    !! CEA Rocket Solver

    use cea_param, only: dp, empty_dp, &
                         R=>gas_constant, &
                         snl=>species_name_len
    use cea_mixture, only: Mixture, MixtureThermo
    use cea_equilibrium, only: EqSolution, EqSolver, EqPartials
    use cea_transport, only: TransportDB
    use fb_utils
    implicit none

    type :: RocketSolver
        !! Rocket solver class

        type(EqSolver) :: eq_solver
            !! Equilibrium solver

    contains

        procedure :: solve => RocketSolver_solve
        procedure :: solve_iac => RocketSolver_solve_iac
        procedure :: solve_fac => RocketSolver_solve_fac
        procedure :: solve_throat => RocketSolver_solve_throat
        procedure :: solve_throat_frozen => RocketSolver_solve_throat_frozen
        procedure :: solve_pi_p => RocketSolver_solve_pi_p
        procedure :: solve_pi_p_frozen => RocketSolver_solve_pi_p_frozen
        procedure :: solve_subar => RocketSolver_solve_subar
        procedure :: solve_supar => RocketSolver_solve_supar
        procedure :: solve_supar_frozen => RocketSolver_solve_supar_frozen
        procedure :: frozen => RocketSolver_frozen
        procedure :: post_process => RocketSolver_post_process
        procedure :: set_init_state => RocketSolver_set_init_state

    end type
    interface RocketSolver
        module procedure :: RocketSolver_init
    end interface

    type :: RocketSolution
        !! Rocket solution class

        integer :: num_pts
            !! Number of evaluation points
        character(8), allocatable :: station(:)
            !! Name of each evaluation point: infinity, chamber, throat, exit
        integer :: throat_idx
            !! Index of the throat station

        ! Initial guess variables
        integer :: i_save
            !! Determines which previous solution index to use as the initial guess
        type(EqSolution) :: eq_soln_save
            !! Stored initial guess for future solves

        ! Thermo solution at each evaluation point
        type(EqSolution), allocatable :: eq_soln(:)
            !! Mole fractions and temperature
        type(EqPartials), allocatable :: eq_partials(:)
            !! Partial derivatives of the equilibrium solution

        ! States
        ! real(dp), allocatable :: Pinf_P_ratio(:)
        real(dp), allocatable :: pressure(:)
            !! Pressure [bar]
        real(dp), allocatable :: mach(:)
            !! Mach number
        real(dp), allocatable :: gamma_s(:)
            !! Isentropic exponent (Eq. 2.71)
        real(dp), allocatable :: v_sonic(:)
            !! Speed of sound [m/s]

        ! Performance parameters
        real(dp), allocatable :: ae_at(:)
            !! Area ratio Ae/At
        real(dp), allocatable :: c_star(:)
            !! Characteristic velocity [m/s]
        real(dp), allocatable :: cf(:)
            !! Thrust coefficient
        real(dp), allocatable :: i_vac(:)
            !! Specific impulse (at vacuum) [m/s]
        real(dp), allocatable :: i_sp(:)
            !! Specific impulse (ambient pressure) [m/s]

        ! Convergence variables
        logical :: converged = .false.
            !! Convergence flag

    end type
    interface RocketSolution
        module procedure :: RocketSolution_init
    end interface

contains

    !-----------------------------------------------------------------------
    ! RocketSolver
    !-----------------------------------------------------------------------
    function RocketSolver_init(products, reactants, trace, ions, all_transport, insert) result(self)
        type(RocketSolver) :: self
        type(Mixture), intent(in) :: products
        type(Mixture), intent(in), optional :: reactants
        real(dp), intent(in), optional :: trace
        logical, intent(in), optional :: ions
        type(TransportDB), intent(in), optional :: all_transport
        character(*), intent(in), optional :: insert(:)  ! List of condensed species to insert

        call log_debug("Initializing RocketSolver")

        self%eq_solver = EqSolver(products, reactants, trace=trace, ions=ions, &
                                  all_transport=all_transport, insert=insert)

        call log_debug("RocketSolver initialized successfully")

    end function

    subroutine RocketSolver_frozen(self, soln, idx, n_frz)

        ! Arguments
        class(RocketSolver), intent(in) :: self
        type(RocketSolution), intent(inout) :: soln
        integer, intent(in) :: idx
        integer, intent(in) :: n_frz

        ! Locals
        integer :: i, j                      ! Loop index
        integer, parameter :: max_iter_frozen = 8  ! Maximum number of iterations for frozen conditions
        integer :: ng                        ! Number of gas species
        logical :: convg                     ! Convergence flag
        logical :: finalized                 ! True if frozen properties were finalized in-loop
        logical :: in_range                  ! True if frozen point is within condensed species temperature ranges
        real(dp) :: dlpm
        real(dp), parameter :: tol = 0.5d-4  ! Tolerance for frozen convergence
        real(dp), parameter :: approx_zero_tol = 1.0d-12  ! Tolerance to check if a value is approximately zero
        real(dp), parameter :: phase_gap = 50.0d0  ! Condensed phase guard band [K]
        real(dp) :: cpsum, ssum              ! Temporary variables for mixture properties
        real(dp) :: cpj, sj                  ! Temporary variables for species properties
        real(dp) :: dlnt                     ! Update variable for log-temperature
        real(dp) :: T_low, T_high            ! Condensed species temperature bounds [K]

        call log_debug("Starting frozen calculations")

        ! Shorthand
        ng = self%eq_solver%num_gas

        ! Set the equilibrium solution values
        soln%eq_soln(idx)%nj = soln%eq_soln(n_frz)%nj
        soln%eq_soln(idx)%n = soln%eq_soln(n_frz)%n
        soln%eq_soln(idx)%mole_fractions = soln%eq_soln(n_frz)%mole_fractions
        soln%eq_soln(idx)%mass_fractions = soln%eq_soln(n_frz)%mass_fractions
        soln%eq_soln(idx)%M = soln%eq_soln(n_frz)%M
        soln%eq_soln(idx)%MW = soln%eq_soln(n_frz)%MW

        dlpm = log(soln%pressure(idx)/soln%eq_soln(n_frz)%n)
        do j = 1, ng
            if ( abs(soln%eq_soln(n_frz)%nj(j)) > approx_zero_tol) then
                soln%eq_soln(idx)%dln_nj(j) = -(log(soln%eq_soln(n_frz)%nj(j))+dlpm)
            end if
        end do

        convg = .false.
        finalized = .false.
        do i = 1, max_iter_frozen
            cpsum = 0.0d0
            ssum = 0.0d0
            do j = 1, self%eq_solver%num_gas
                cpj = self%eq_solver%products%species(j)%calc_cp(soln%eq_soln(idx)%T)
                sj = self%eq_solver%products%species(j)%calc_entropy(soln%eq_soln(idx)%T)
                cpsum = cpsum + soln%eq_soln(n_frz)%nj(j)*cpj
                ssum = ssum + soln%eq_soln(n_frz)%nj(j)*(sj+soln%eq_soln(idx)%dln_nj(j))
            end do
            do j = 1, self%eq_solver%num_condensed
                if (.not. soln%eq_soln(n_frz)%is_active(j)) cycle
                cpj = self%eq_solver%products%species(ng+j)%calc_cp(soln%eq_soln(idx)%T)
                sj = self%eq_solver%products%species(ng+j)%calc_entropy(soln%eq_soln(idx)%T)
                cpsum = cpsum + soln%eq_soln(n_frz)%nj(ng+j)*cpj
                ssum = ssum + soln%eq_soln(n_frz)%nj(ng+j)*sj
            end do

            if (convg) then
                soln%eq_partials(idx)%gamma_s = cpsum/(cpsum - soln%eq_soln(n_frz)%n)
                soln%eq_soln(idx)%gamma_s = soln%eq_partials(idx)%gamma_s
                soln%gamma_s(idx) = soln%eq_partials(idx)%gamma_s
                soln%eq_partials(idx)%dlnV_dlnP = -1.0d0
                soln%eq_partials(idx)%dlnV_dlnT = 1.0d0
                call self%eq_solver%products%calc_thermo(soln%eq_soln(idx)%thermo, soln%eq_soln(idx)%T)

                ! Stop if any frozen condensed species is outside its
                ! valid temperature range by more than 50 K.
                in_range = .true.
                do j = 1, self%eq_solver%num_condensed
                    if (abs(soln%eq_soln(n_frz)%nj(ng+j)) <= approx_zero_tol) cycle
                    T_low = minval(self%eq_solver%products%species(ng+j)%T_fit(:, 1))
                    T_high = maxval(self%eq_solver%products%species(ng+j)%T_fit(:, 2))
                    if (soln%eq_soln(idx)%T < (T_low-phase_gap) .or. soln%eq_soln(idx)%T > (T_high+phase_gap)) then
                        in_range = .false.
                        exit
                    end if
                end do

                if (.not. in_range) then
                    call log_warning("Frozen calculations stopped: temperature is more than 50 K outside "// &
                        "the range of a condensed species")
                    soln%converged = .false.
                    return
                end if

                finalized = .true.
                exit
            else
                dlnt = (1.d3*soln%eq_soln(n_frz)%entropy/R - ssum)/cpsum
                soln%eq_soln(idx)%T = exp(log(soln%eq_soln(idx)%T) + dlnt)
                if (abs(dlnt) < tol) convg = .true.
            end if

        end do

        if (.not. finalized) then
            call log_warning("Frozen calculations did not converge in 8 iterations")
            soln%converged = .false.
            return
        end if

        ! Compute and save the mixture properties
        soln%eq_soln(idx)%pressure = soln%pressure(idx)
        soln%eq_soln(idx)%volume = (1.d-5 * R * soln%eq_soln(idx)%n * soln%eq_soln(idx)%T) / soln%pressure(idx)
        soln%eq_soln(idx)%density = 1.0d0/soln%eq_soln(idx)%volume
        soln%eq_soln(idx)%cp_fr = dot_product(soln%eq_soln(idx)%thermo%cp, soln%eq_soln(idx)%nj) * R / 1.d3
        soln%eq_soln(idx)%cp_eq = soln%eq_soln(idx)%cp_fr
        soln%eq_soln(idx)%enthalpy = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%enthalpy)*R*soln%eq_soln(idx)%T/1.d3
        soln%eq_soln(idx)%energy = soln%eq_soln(idx)%enthalpy - soln%eq_soln(idx)%n*soln%eq_soln(idx)%T*R/1.d3
        soln%eq_soln(idx)%entropy = soln%eq_soln(n_frz)%entropy
        soln%eq_soln(idx)%gibbs_energy = (soln%eq_soln(idx)%enthalpy - soln%eq_soln(idx)%T*soln%eq_soln(idx)%entropy)

    end subroutine

    subroutine RocketSolver_solve_throat(self, soln, idx, pc, h_inf, s0, weights, awt)

        ! Arguments
        class(RocketSolver), intent(in) :: self
        type(RocketSolution), intent(inout) :: soln
        integer, intent(in) :: idx                     ! Station index
        real(dp), intent(in) :: pc                     ! Chamber pressure [Pa]
        real(dp), intent(in) :: h_inf                  ! Enthalpy at infinity
        real(dp), intent(in) :: s0                     ! Fixed entropy for equilibrium solve
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: awt                   ! Mass flow per area in the throat

        ! Locals
        integer :: i                         ! Loop index
        integer, parameter :: max_iter_throat = 22  ! Maximum number of iterations for throat conditions
        real(dp), parameter :: ut_tol = 0.4d-4  ! Tolerance for throat velocity convergence
        real(dp) :: p, delta_p               ! Temporary pressure variable
        real(dp) :: usq, asq                 ! velocity squared; sonic velocity squared
        real(dp) :: h                        ! Enthalpy at any other station (temporary)
        real(dp) :: gamma_s                  ! Temp variable for isentropic exponent gamma_s
        real(dp) :: T_melt                   ! Melting temperature of the condensed species
        real(dp) :: dlt                      ! Log(T_melt/T)
        real(dp) :: cp                       ! Mixture specific heat at constant pressure
        real(dp) :: dlnV_dlnT                ! Partial derivative of ln(V) with respect to ln(T)

        call log_debug("Starting equilibrium throat calculations")

        ! Initialization
        soln%station(idx) = "throat  "
        T_melt = 0.0d0

        ! Initial estimate pressure ratio (Eq. 6.15)
        gamma_s = soln%eq_partials(1)%gamma_s
        soln%pressure(idx) = pc/((0.5d0*gamma_s+0.5d0)**(gamma_s/(gamma_s-1.0d0)))
        soln%ae_at(idx) = 1.0d0

        ! Set the initial guess for the equilibrium solve
        soln%eq_soln(idx) = EqSolution(self%eq_solver)!, T_init=soln%eq_soln(1)%T, nj_init=soln%eq_soln(1)%nj)
        call self%set_init_state(soln, idx)

        do i = 1, max_iter_throat

            ! Solve the equilibrium problem using the previous solution as the initial guess
            call self%eq_solver%solve(soln%eq_soln(idx), "sp", s0, soln%pressure(idx), weights, partials=soln%eq_partials(idx))

            ! Compute throat properties
            h = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%enthalpy)*soln%eq_soln(idx)%T
            gamma_s = soln%eq_partials(idx)%gamma_s
            asq = soln%eq_soln(idx)%n*R*gamma_s*soln%eq_soln(idx)%T
            usq = 2.0d0*(h_inf-h)*R
            soln%v_sonic(idx) = sqrt(asq)
            soln%mach(idx) = sqrt(usq/asq)

            awt = soln%eq_soln(idx)%n*soln%eq_soln(idx)%T/(soln%pressure(idx)*sqrt(usq))

            ! Check throat convergence tolerance (Eq. 6.16)
            if (abs(usq - asq)/usq <= ut_tol) exit

            ! Update estimate pressure if not converged (Eq. 6.17)
            p = soln%pressure(idx)*((1.0d0 + gamma_s*(usq/asq))/(1.0d0 + gamma_s))
            if (i <= 3) then
                if (soln%eq_soln(idx)%j_sol /= 0) then
                    T_melt = soln%eq_soln(idx)%T
                    soln%pressure(idx) = p
                else if (abs(T_melt) < 1.0d-6) then
                    soln%pressure(idx) = p
                else  ! T_melt /= 0
                    dlt = log(T_melt/soln%eq_soln(idx)%T)
                    cp = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%cp)
                    dlnV_dlnT = soln%eq_partials(idx)%dlnV_dlnT
                    soln%pressure(idx) = soln%pressure(idx)*exp(dlt*cp/(soln%eq_soln(idx)%n*dlnV_dlnT))
                    exit
                end if
            else
                delta_p = abs(p-soln%pressure(idx))/20.d0
                soln%pressure(idx) = max(p, soln%pressure(idx)) - delta_p
            end if

        end do

        soln%i_save = idx
        ! Update i_save if solid and liquid are both present
        if (soln%eq_soln(idx)%j_liq /= 0 .and. soln%i_save > 0) then
            soln%i_save = 0
        end if

    end subroutine

    subroutine RocketSolver_solve_throat_frozen(self, soln, idx, n_frz, pc, h_inf, awt)

        ! Arguments
        class(RocketSolver), intent(in) :: self
        type(RocketSolution), intent(inout) :: soln
        integer, intent(in) :: idx                     ! Station index
        integer, intent(in) :: n_frz                   ! Frozen index
        real(dp), intent(in) :: pc                     ! Chamber pressure [Pa]
        real(dp), intent(in) :: h_inf                  ! Enthalpy at infinity
        real(dp), intent(out) :: awt                   ! Mass flow per area in the throat

        ! Locals
        integer :: i                         ! Loop index
        integer, parameter :: max_iter_throat = 22  ! Maximum number of iterations for throat conditions
        real(dp), parameter :: ut_tol = 0.4d-4  ! Tolerance for throat velocity convergence
        real(dp) :: p, delta_p               ! Temporary pressure variable
        real(dp) :: usq, asq                 ! velocity squared; sonic velocity squared
        real(dp) :: h                        ! Enthalpy at any other station (temporary)
        real(dp) :: gamma_s                  ! Temp variable for isentropic exponent gamma_s

        call log_debug("Starting frozen throat calculations")
        soln%station(idx) = "throat  "
        awt = 0.0d0

        ! Initial estimate pressure ratio (Eq. 6.15)
        gamma_s = soln%gamma_s(n_frz)
        soln%pressure(idx) = pc/((0.5d0*gamma_s+0.5d0)**(gamma_s/(gamma_s-1.0d0)))
        soln%ae_at(idx) = 1.0d0

        ! Set the initial guess for the equilibrium solve
        soln%eq_soln(idx) = EqSolution(self%eq_solver, T_init=soln%eq_soln(n_frz)%T)

        ! Set the equilibrium solution as the frozen point
        soln%eq_soln(idx) = soln%eq_soln(n_frz)
        soln%eq_partials(idx) = soln%eq_partials(n_frz)

        do i = 1, max_iter_throat

            ! Compute the frozen properties
            call self%frozen(soln, idx, n_frz)
            if (.not. soln%converged) return

            ! Compute throat properties
            h = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%enthalpy)*soln%eq_soln(idx)%T
            gamma_s = soln%eq_partials(idx)%gamma_s
            asq = soln%eq_soln(idx)%n*R*gamma_s*soln%eq_soln(idx)%T
            usq = 2.0d0*(h_inf-h)*R
            soln%v_sonic(idx) = sqrt(asq)
            soln%mach(idx) = sqrt(usq/asq)
            awt = soln%eq_soln(idx)%n*soln%eq_soln(idx)%T/(soln%pressure(idx)*sqrt(usq))

            ! Check throat convergence tolerance (Eq. 6.16)
            if (abs(usq - asq)/usq <= ut_tol) exit

            ! Update estimate pressure if not converged (Eq. 6.17)
            p = soln%pressure(idx)*((1.0d0 + gamma_s*(usq/asq))/(1.0d0 + gamma_s))
            if (i > 3) then
                delta_p = dabs(soln%pressure(idx)-p)/20.d0
                soln%pressure(idx) = max(p, soln%pressure(idx)) - delta_p
            else
                soln%pressure(idx) = p
            end if

        end do

    end subroutine

    subroutine RocketSolver_solve_pi_p(self, soln, idx, pc, pi_p, h_inf, s0, weights)

        ! Arguments
        class(RocketSolver), intent(in) :: self
        type(RocketSolution), intent(inout) :: soln
        integer, intent(inout) :: idx                ! Initial station index
        real(dp), intent(in) :: pc                   ! Chamber pressure
        real(dp), intent(in) :: pi_p(:)              ! Exit pressure ratio
        real(dp), intent(in) :: h_inf                ! Enthalpy at infinity
        real(dp), intent(in) :: s0                   ! Fixed entropy for equilibrium solve
        real(dp), intent(in) :: weights(:)           ! Reactant weights for equilibrium solve

        ! Locals
        integer :: i                         ! Loop index
        real(dp) :: usq, asq                 ! velocity squared; sonic velocity squared
        real(dp) :: awt                      ! Throat area per unit mass flow rate
        real(dp) :: h                        ! Enthalpy at any other station (temporary)
        real(dp) :: gamma_s                  ! Temp variable for isentropic exponent gamma_s

        call log_debug("Starting equilibrium pi/p calculations")

        ! Set the initial guess for the equilibrium solve based on throat conditions
        soln%eq_soln(idx) = EqSolution(self%eq_solver)!, T_init=soln%eq_soln(T_idx)%T, nj_init=soln%eq_soln(nj_idx)%nj)

        awt = soln%eq_soln(soln%throat_idx)%n*soln%eq_soln(soln%throat_idx)%T/ &
            (soln%pressure(soln%throat_idx)*soln%v_sonic(soln%throat_idx))
        do i = 1, size(pi_p)
            soln%station(idx) = "exit    "

            ! Get the pressure from the pressure rato
            soln%pressure(idx) = pc/pi_p(i)

            ! Solve the equilibrium problem
            call self%set_init_state(soln, idx)
            call self%eq_solver%solve(soln%eq_soln(idx), "sp", s0, soln%pressure(idx), weights, partials=soln%eq_partials(idx))

            ! Compute exit properties
            h = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%enthalpy)*soln%eq_soln(idx)%T
            gamma_s = soln%eq_partials(idx)%gamma_s
            asq = soln%eq_soln(idx)%n*R*gamma_s*soln%eq_soln(idx)%T
            usq = 2.0d0*(h_inf-h)*R
            soln%v_sonic(idx) = sqrt(asq)
            soln%mach(idx) = sqrt(usq/asq)
            soln%ae_at(idx) = soln%eq_soln(idx)%n*soln%eq_soln(idx)%T/(soln%pressure(idx)*usq**0.5*awt)

            ! Set the initial conditions for the next point
            if (i < size(pi_p)) then
                soln%eq_soln(idx+1) = EqSolution(self%eq_solver, T_init=soln%eq_soln(idx)%T, nj_init=soln%eq_soln(idx)%nj)
            end if

            ! Update i_save if solid and liquid are both present
            soln%i_save = idx
            if (soln%eq_soln(idx)%j_liq /= 0 .and. soln%i_save > 0) then
                soln%i_save = 0
            end if

            ! Iterate
            idx = idx + 1

        end do

    end subroutine

    subroutine RocketSolver_solve_pi_p_frozen(self, soln, idx, n_frz, pc, pi_p, h_inf, T_idx)

        ! Arguments
        class(RocketSolver), intent(in) :: self
        type(RocketSolution), intent(inout) :: soln
        integer, intent(inout) :: idx                ! Initial station index
        integer, intent(in) :: n_frz                 ! Frozen index
        real(dp), intent(in) :: pc                   ! Chamber pressure
        real(dp), intent(in) :: pi_p(:)              ! Exit pressure ratio
        real(dp), intent(in) :: h_inf                ! Enthalpy at infinity
        integer, intent(in) :: T_idx                 ! Solution index of initial temperature for equilibrium solve

        ! Locals
        integer :: i                         ! Loop index
        real(dp) :: usq, asq                 ! velocity squared; sonic velocity squared
        real(dp) :: awt                      ! Throat area per unit mass flow rate
        real(dp) :: h                        ! Enthalpy at any other station (temporary)
        real(dp) :: gamma_s                  ! Temp variable for isentropic exponent gamma_s
        real(dp) :: pip_nf                   ! Pressure ratio at freeze point

        call log_debug("Starting frozen pi/p calculations")

        awt = soln%eq_soln(soln%throat_idx)%n*soln%eq_soln(soln%throat_idx)%T/ &
            (soln%pressure(soln%throat_idx)*soln%v_sonic(soln%throat_idx))
        pip_nf = pc/soln%pressure(n_frz)
        do i = 1, size(pi_p)
            ! Legacy frozen scheduling: omit assigned pressure ratios lower than
            ! the value at the freeze point.
            if (pi_p(i) < pip_nf) then
                call log_info('RocketSolver: WARNING!!  FOR FROZEN PERFORMANCE, POINT OMITTED BECAUSE '// &
                    'ASSIGNED pi/p IS LESS THAN VALUE AT nfz='//to_str(n_frz))
                cycle
            end if

            soln%station(idx) = "exit    "

            ! Get the pressure from the pressure rato
            soln%pressure(idx) = pc/pi_p(i)

            ! Set the initial guess for the equilibrium solve based on throat conditions
            soln%eq_soln(idx) = EqSolution(self%eq_solver, T_init=soln%eq_soln(T_idx)%T)

            ! Compute the frozen solution
            call self%frozen(soln, idx, n_frz)
            if (.not. soln%converged) return

            ! Compute exit properties
            h = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%enthalpy)*soln%eq_soln(idx)%T
            gamma_s = soln%eq_partials(idx)%gamma_s
            asq = soln%eq_soln(idx)%n*R*gamma_s*soln%eq_soln(idx)%T
            usq = 2.0d0*(h_inf-h)*R
            soln%v_sonic(idx) = sqrt(asq)
            soln%mach(idx) = sqrt(usq/asq)
            soln%ae_at(idx) = soln%eq_soln(idx)%n*soln%eq_soln(idx)%T/(soln%pressure(idx)*usq**0.5*awt)

            idx = idx + 1

        end do

    end subroutine

    subroutine RocketSolver_solve_subar(self, soln, idx, pc, subar, h_inf, s0, weights, T_idx, nj_idx, ln_pinf_pt, awt)

        ! Arguments
        class(RocketSolver), intent(in) :: self
        type(RocketSolution), intent(inout) :: soln
        integer, intent(inout) :: idx                  ! Initial station index
        real(dp), intent(in) :: pc                     ! Chamber pressure
        real(dp), intent(in) :: subar(:)               ! Subsonic area ratio
        real(dp), intent(in) :: h_inf                  ! Enthalpy at infinity
        real(dp), intent(in) :: s0                     ! Fixed entropy for equilibrium solve
        real(dp), intent(in) :: weights(:)             ! Reactant weights for equilibrium solve
        integer, intent(in) :: T_idx                   ! Initial temperature for equilibrium solve
        integer, intent(in) :: nj_idx                  ! Initial mole fractions for equilibrium solve
        real(dp), intent(in) :: ln_pinf_pt             ! ln(Pinf/Pt)
        real(dp), intent(in) :: awt                    ! Mass flow per area in the throat

        ! Locals
        integer :: i, j                      ! Loop index
        integer, parameter :: max_iter_area = 10  ! Maximum number of iterations for exit condition using area ratio
        real(dp), parameter :: area_tol = 4.0d-5  ! Area-ratio convergence tolerance
        real(dp) :: usq, asq                 ! velocity squared; sonic velocity squared
        real(dp) :: h                        ! Enthalpy at any other station (temporary)
        real(dp) :: gamma_s                  ! Temp variable for isentropic exponent gamma_s
        real(dp) :: ln_pinf_pe               ! Temporary variable for ln(Pinf/Pe)
        real(dp) :: dln_pinf_pe_dln_aeat     ! Partial derivative ∂ln(Pinf/Pe)/∂ln(Ae/At)
        real(dp) :: dln_pinf_pe              ! Update to ln(Pinf/Pe)

        call log_debug("Starting subar calculations")

        do i = 1, size(subar)
            soln%station(idx) = "exit    "

            ! Set the initial guess for the equilibrium solve based on throat conditions
            soln%eq_soln(idx) = EqSolution(self%eq_solver)!, T_init=soln%eq_soln(T_idx)%T, nj_init=soln%eq_soln(nj_idx)%nj)
            call self%set_init_state(soln, idx)

            ! Check for valid area ratio
            if (subar(i) .le. 1.0d0) then
                call abort("Subsonic area ratio must be greater than 1.0")
            end if

            ! Compute intial estimate of pressure ratio (Eq. 6.19/6.20)
            ln_pinf_pe = ln_pinf_pt/(subar(i) + 10.587d0*(log(subar(i))**3.0d0) + 9.454d0*log(subar(i)))
            if (subar(i) < 1.09d0) then
                ln_pinf_pe = 0.9d0*ln_pinf_pe
            else if (subar(i) > 10.0d0) then
                ln_pinf_pe = ln_pinf_pe/subar(i)
            end if

            do j = 1, max_iter_area

                ! Solve the equilibrium problem
                soln%pressure(idx) = pc/exp(ln_pinf_pe)
                call self%eq_solver%solve(soln%eq_soln(idx), "sp", s0, soln%pressure(idx), weights, partials=soln%eq_partials(idx))

                ! Compute exit properties
                h = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%enthalpy)*soln%eq_soln(idx)%T
                gamma_s = soln%eq_partials(idx)%gamma_s
                asq = soln%eq_soln(idx)%n*R*gamma_s*soln%eq_soln(idx)%T
                usq = 2.0d0*(h_inf-h)*R
                soln%v_sonic(idx) = sqrt(asq)
                soln%mach(idx) = sqrt(usq/asq)
                soln%ae_at(idx) = soln%eq_soln(idx)%n*soln%eq_soln(idx)%T/(soln%pressure(idx)*sqrt(usq)*awt)

                ! Compute updated pressure ratio estimate
                dln_pinf_pe_dln_aeat = gamma_s*usq/(usq - asq)  ! (Eq. 6.23)
                dln_pinf_pe = dln_pinf_pe_dln_aeat*(log(subar(i)) - log(soln%ae_at(idx)))

                ! Convergence test for assigned area ratio:
                ! relative Ae/At error OR small pressure-ratio update.
                if (abs(soln%ae_at(idx)-subar(i))/subar(i) <= area_tol) exit
                if (abs(dln_pinf_pe) < area_tol) exit
                ln_pinf_pe = ln_pinf_pe + dln_pinf_pe  ! If not converged, update estimate

            end do

            ! Update i_save if solid and liquid are both present
            soln%i_save = idx
            if (soln%eq_soln(idx)%j_liq /= 0 .and. soln%i_save > 0) then
                soln%i_save = 0
            end if

            ! Iterate to the next exit condition
            idx = idx + 1

        end do

    end subroutine

    subroutine RocketSolver_solve_supar(self, soln, idx, pc, supar, h_inf, s0, weights, ln_pinf_pt, awt)

        ! Arguments
        class(RocketSolver), intent(in) :: self
        type(RocketSolution), intent(inout) :: soln
        integer, intent(inout) :: idx                  ! Initial station index
        real(dp), intent(in) :: pc                     ! Chamber pressure
        real(dp), intent(in) :: supar(:)               ! Supersonic area ratio
        real(dp), intent(in) :: h_inf                  ! Enthalpy at infinity
        real(dp), intent(in) :: s0                     ! Fixed entropy for equilibrium solve
        real(dp), intent(in) :: weights(:)             ! Reactant weights for equilibrium solve
        real(dp), intent(in) :: ln_pinf_pt             ! ln(Pinf/Pt)
        real(dp), intent(in) :: awt                    ! Mass flow per area in the throat

        ! Locals
        integer :: i, j                      ! Loop index
        integer, parameter :: max_iter_area = 10  ! Maximum number of iterations for exit condition using area ratio
        real(dp), parameter :: area_tol = 4.0d-5  ! Area-ratio convergence tolerance
        real(dp) :: usq, asq                 ! velocity squared; sonic velocity squared
        real(dp) :: h                        ! Enthalpy at any other station (temporary)
        real(dp) :: gamma_s                  ! Temp variable for isentropic exponent gamma_s
        real(dp) :: ln_pinf_pe               ! Temporary variable for ln(Pinf/Pe)
        real(dp) :: dln_pinf_pe_dln_aeat     ! Partial derivative ∂ln(Pinf/Pe)/∂ln(Ae/At)
        real(dp) :: dln_pinf_pe              ! Update to ln(Pinf/Pe)

        call log_debug("Starting equilibrium supar calculations")

        do i = 1, size(supar)
            soln%station(idx) = "exit    "

            ! Set the initial guess for the equilibrium solve based on throat conditions
            soln%eq_soln(idx) = EqSolution(self%eq_solver)
            call self%set_init_state(soln, idx)

            ! Check for valid area ratio
            if (supar(i) .le. 1.0d0) then
                call abort("Supersonic area ratio must be greater than 1.0")
            end if

            ! Compute intial estimate of pressure ratio (Eq. 6.21/6.22)
            if (supar(i) < 2.0d0) then
                ln_pinf_pe = ln_pinf_pt + sqrt(3.294d0*(log(supar(i))**2.0d0) + 1.535d0*log(supar(i)))
            else if (supar(i) >= 2.0d0) then
                ln_pinf_pe = soln%eq_partials(2)%gamma_s + 1.4d0*log(supar(i))
            end if

            do j = 1, max_iter_area

                ! Solve the equilibrium problem
                soln%pressure(idx) = pc/exp(ln_pinf_pe)
                call self%eq_solver%solve(soln%eq_soln(idx), "sp", s0, soln%pressure(idx), weights, partials=soln%eq_partials(idx))

                ! Compute exit properties
                h = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%enthalpy)*soln%eq_soln(idx)%T
                gamma_s = soln%eq_partials(idx)%gamma_s
                asq = soln%eq_soln(idx)%n*R*gamma_s*soln%eq_soln(idx)%T
                usq = 2.0d0*(h_inf-h)*R
                soln%v_sonic(idx) = sqrt(asq)
                soln%mach(idx) = sqrt(usq/asq)
                soln%ae_at(idx) = soln%eq_soln(idx)%n*soln%eq_soln(idx)%T/(soln%pressure(idx)*sqrt(usq)*awt)

                ! Compute updated pressure ratio estimate
                dln_pinf_pe_dln_aeat = gamma_s*usq/(usq - asq)  ! (Eq. 6.23)
                dln_pinf_pe = dln_pinf_pe_dln_aeat*(log(supar(i)) - log(soln%ae_at(idx)))

                ! Convergence test for assigned area ratio:
                ! relative Ae/At error OR small pressure-ratio update.
                if (abs(soln%ae_at(idx)-supar(i))/supar(i) <= area_tol) exit
                if (abs(dln_pinf_pe) < area_tol) exit
                ln_pinf_pe = ln_pinf_pe + dln_pinf_pe  ! If not converged, update estimate

            end do

            ! Update i_save if solid and liquid are both present
            soln%i_save = idx
            if (soln%eq_soln(idx)%j_liq /= 0 .and. soln%i_save > 0) then
                soln%i_save = 0
            end if

            idx = idx + 1

        end do

    end subroutine

    subroutine RocketSolver_solve_supar_frozen(self, soln, idx, n_frz, pc, supar, h_inf, weights, T_idx, ln_pinf_pt, awt)

        ! Arguments
        class(RocketSolver), intent(in) :: self
        type(RocketSolution), intent(inout) :: soln
        integer, intent(inout) :: idx                  ! Initial station index
        integer, intent(in) :: n_frz                   ! Frozen index
        real(dp), intent(in) :: pc                     ! Chamber pressure
        real(dp), intent(in) :: supar(:)               ! Supersonic area ratio
        real(dp), intent(in) :: h_inf                  ! Enthalpy at infinity
        real(dp), intent(in) :: weights(:)             ! Reactant weights for equilibrium solve
        integer, intent(in) :: T_idx                   ! Initial temperature for equilibrium solve
        real(dp), intent(in) :: ln_pinf_pt             ! ln(Pinf/Pt)
        real(dp), intent(in) :: awt                    ! Mass flow per area in the throat

        ! Locals
        integer :: i, j                      ! Loop index
        integer, parameter :: max_iter_area = 10  ! Maximum number of iterations for exit condition using area ratio
        real(dp), parameter :: area_tol = 4.0d-5  ! Area-ratio convergence tolerance
        real(dp) :: usq, asq                 ! velocity squared; sonic velocity squared
        real(dp) :: h                        ! Enthalpy at any other station (temporary)
        real(dp) :: gamma_s                  ! Temp variable for isentropic exponent gamma_s
        real(dp) :: ln_pinf_pe               ! Temporary variable for ln(Pinf/Pe)
        real(dp) :: dln_pinf_pe_dln_aeat     ! Partial derivative ∂ln(Pinf/Pe)/∂ln(Ae/At)
        real(dp) :: dln_pinf_pe              ! Update to ln(Pinf/Pe)

        call log_debug("Starting frozen supar calculations")

        do i = 1, size(supar)
            ! Frozen scheduling: omit assigned supersonic area ratios
            ! that are not greater than the value at the freeze point.
            if (n_frz >= 3 .and. supar(i) <= soln%ae_at(n_frz)) then
                call log_info('RocketSolver: WARNING!!  FOR FROZEN PERFORMANCE, POINT OMITTED BECAUSE '// &
                    'ASSIGNED Ae/At IS LESS THAN OR EQUAL TO VALUE AT nfz='//to_str(n_frz))
                cycle
            end if

            soln%station(idx) = "exit    "

            ! Set the initial guess for the equilibrium solve based on throat conditions
            soln%eq_soln(idx) = EqSolution(self%eq_solver, T_init=soln%eq_soln(T_idx)%T)

            ! Check for valid area ratio
            if (supar(i) .le. 1.0d0) then
                call abort("Supersonic area ratio must be greater than 1.0")
            end if

            ! Compute intial estimate of pressure ratio (Eq. 6.21/6.22)
            if (supar(i) < 2.0d0) then
                ln_pinf_pe = ln_pinf_pt + sqrt(3.294d0*(log(supar(i))**2.0d0) + 1.535d0*log(supar(i)))
            else if (supar(i) >= 2.0d0) then
                ln_pinf_pe = soln%eq_partials(2)%gamma_s + 1.4d0*log(supar(i))
            end if

            do j = 1, max_iter_area

                ! Compute the frozen solution
                soln%pressure(idx) = pc/exp(ln_pinf_pe)
                call self%frozen(soln, idx, n_frz)
                if (.not. soln%converged) return

                ! Compute exit properties
                h = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%enthalpy)*soln%eq_soln(idx)%T
                gamma_s = soln%eq_partials(idx)%gamma_s
                asq = soln%eq_soln(idx)%n*R*gamma_s*soln%eq_soln(idx)%T
                usq = 2.0d0*(h_inf-h)*R
                soln%v_sonic(idx) = sqrt(asq)
                soln%mach(idx) = sqrt(usq/asq)
                soln%ae_at(idx) = soln%eq_soln(idx)%n*soln%eq_soln(idx)%T/(soln%pressure(idx)*sqrt(usq)*awt)

                ! Compute updated pressure ratio estimate
                dln_pinf_pe_dln_aeat = gamma_s*usq/(usq - asq)  ! (Eq. 6.23)
                dln_pinf_pe = dln_pinf_pe_dln_aeat*(log(supar(i)) - log(soln%ae_at(idx)))

                ! Convergence test for assigned area ratio:
                ! relative Ae/At error OR small pressure-ratio update.
                if (abs(soln%ae_at(idx)-supar(i))/supar(i) <= area_tol) exit
                if (abs(dln_pinf_pe) < area_tol) exit
                ln_pinf_pe = ln_pinf_pe + dln_pinf_pe  ! If not converged, update estimate

            end do

            idx = idx + 1

        end do

    end subroutine

    subroutine RocketSolver_solve_iac(self, soln, reactant_weights, pc, pi_p, subar, supar, n_frz, tc_est, hc, tc)
        ! Solve the inifit-area combustor (IAC) rocket problem

        ! Arguments
        class(RocketSolver) :: self
        type(RocketSolution), intent(inout) :: soln
        real(dp), intent(in) :: reactant_weights(:)
        real(dp), intent(in) :: pc                   ! Chamber pressure [Pa]
        real(dp), intent(in) :: pi_p(:)              ! Ratio of chamber pressure to exit pressure [unitless]
        real(dp), intent(in), optional :: subar(:)   ! Subsonic area ratio [unitless]
        real(dp), intent(in), optional :: supar(:)   ! Supersonic area ratio [unitless]
        integer,  intent(in), optional :: n_frz      ! Station where the composition should be frozen
        real(dp), intent(in), optional :: tc_est     ! Initial chamber temperature estimate [K]
        real(dp), intent(in), optional :: hc         ! Assigned chamber enthalpy [unitless]
        real(dp), intent(in), optional :: tc         ! Assigned chamber temperature [K]

        ! Locals
        integer :: idx                       ! Station index
        integer :: n_frz_                    ! Temporary variable for frozen index
        integer :: num_pts                   ! Total number of evaluation points
        character(len=2) :: prob_type        ! Equilibrium problem type
        real(dp) :: state1                   ! Chamber temperature or enthalpy, or entropy at other stations
        real(dp) :: gamma_s                  ! Temp variable for isentropic exponent gamma_s
        real(dp) :: h_inf                    ! Enthalpy at infinity (chamber); enthalpy at any other station (temporary)
        real(dp) :: ln_pinf_pt               ! Temporary variable for ln(Pinf/Pt)
        real(dp) :: awt                      ! Mass flow per area in the throat
        logical :: frozen                    ! Flag to determine if frozen composition is used

        call log_debug("Starting rocket IAC solve")

        ! Set the total number of evaluation points
        num_pts = 2 + size(pi_p)
        if (present(subar)) num_pts = num_pts + size(subar)
        if (present(supar)) num_pts = num_pts + size(supar)

        ! Initialize the solution variables
        soln = RocketSolution(self, num_pts=num_pts)
        call log_debug("Initialized RocketSolution")
        soln%throat_idx = 2
        soln%converged = .true.

        ! Set the equilibrium problem type
        if (present(hc)) then
            prob_type = "hp"
            state1 = hc
            h_inf = hc
        else
            prob_type = "tp"
            state1 = tc
        end if

        ! Set the frozen composition flag
        if (present(n_frz)) then
            if (n_frz > 0 .and. n_frz <= num_pts) then
                frozen = .true.
                n_frz_ = n_frz
            else
                frozen = .false.
                n_frz_ = 100  ! Set to a large number to ensure it is not used
                if (n_frz /= 0) call log_info('RocketSolver: WARNING!!  INVALID FROZEN INDEX.  FROZEN MODE WILL NOT BE USED')
            end if
        else
            frozen = .false.
            n_frz_ = 100  ! Set to a large number to ensure it is not used
        end if

        ! -----------------------------------------------
        ! Chamber conditions (infinity)
        ! -----------------------------------------------
        if (present(tc_est)) then
            soln%eq_soln(1) = EqSolution(self%eq_solver, T_init=tc_est)
        else
            soln%eq_soln(1) = EqSolution(self%eq_solver)
        end if

        call log_debug("Starting chamber eqsolve")
        call self%eq_solver%solve(soln%eq_soln(1), prob_type, state1, pc, reactant_weights, partials=soln%eq_partials(1))
        soln%i_save = -1

        ! Set the states
        soln%station(1) = "infinity"
        soln%pressure(1) = pc
        soln%mach(1) = 0.0d0

        if (frozen .and. n_frz_ <= 1) then
            gamma_s = (soln%eq_soln(1)%cp_fr/(R*1.d-3))/((soln%eq_soln(1)%cp_fr/(R*1.d-3))-soln%eq_soln(1)%n)
            soln%gamma_s(1) = gamma_s
            soln%eq_partials(1)%gamma_s = gamma_s
            soln%eq_soln(1)%gamma_s = soln%eq_partials(1)%gamma_s
            soln%v_sonic(1) = sqrt(soln%eq_soln(1)%n*R*gamma_s*soln%eq_soln(1)%T)
        else
            gamma_s = soln%eq_partials(1)%gamma_s
            soln%gamma_s(1) = gamma_s
            soln%eq_partials(1)%gamma_s = gamma_s
            soln%eq_soln(1)%gamma_s = soln%eq_partials(1)%gamma_s
            soln%v_sonic(1) = sqrt(soln%eq_soln(1)%n*R*soln%gamma_s(1)*soln%eq_soln(1)%T)
        end if

        ! Save some chamber solution variables for later use
        state1 = soln%eq_soln(1)%calc_entropy_sum(self%eq_solver)  ! Combustor entropy
        if (prob_type == "tp") h_inf = dot_product(soln%eq_soln(1)%nj, soln%eq_soln(1)%thermo%enthalpy)*soln%eq_soln(1)%T  ! Combustor enthalpy

        ! -----------------------------------------------
        ! Throat conditions
        ! -----------------------------------------------
        idx = 2

        if (frozen .and. idx > n_frz_) then
            call self%solve_throat_frozen(soln, idx, n_frz, pc, h_inf, awt)
        else
            call self%solve_throat(soln, idx, pc, h_inf, state1, reactant_weights, awt)
        end if
        if (.not. soln%converged) return
        ln_pinf_pt = log(soln%pressure(1)/soln%pressure(2))

        ! -----------------------------------------------
        ! Exit conditions: pressure ratio
        ! -----------------------------------------------
        idx = 3

        if (frozen .and. idx > n_frz_) then
            call self%solve_pi_p_frozen(soln, idx, n_frz, pc, pi_p, h_inf, 1)
        else
            call self%solve_pi_p(soln, idx, pc, pi_p, h_inf, state1, reactant_weights)
        end if
        if (.not. soln%converged) return

        ! -----------------------------------------------
        ! Exit conditions: subsonic area ratio
        ! -----------------------------------------------

        if (present(subar)) then

            if (frozen .and. n_frz_ > 1) then
                call log_info('RocketSolver: WARNING!!  FREEZING IS NOT ALLOWED AT A SUBSONIC PRESSURE RATIO')
            else
                call self%solve_subar(soln, idx, pc, subar, h_inf, state1, reactant_weights, idx-1, 2, ln_pinf_pt, awt)
            end if

        end if

        ! -----------------------------------------------
        ! Exit conditions: supersonic area ratio
        ! -----------------------------------------------

        if (present(supar)) then

            if (frozen .and. idx > n_frz_) then
                call self%solve_supar_frozen(soln, idx, n_frz, pc, supar, h_inf, reactant_weights, idx-1, ln_pinf_pt, awt)
            else
                call self%solve_supar(soln, idx, pc, supar, h_inf, state1, reactant_weights, ln_pinf_pt, awt)
            end if
            if (.not. soln%converged) return

        end if

        ! Omitted frozen schedule points do not consume output indices.
        soln%num_pts = idx - 1

        ! Compute performance parameters
        call self%post_process(soln, .false.)

    end subroutine

    subroutine RocketSolver_solve_fac(self, soln, reactant_weights, pc, pi_p, subar, supar, ac_at, mdot, n_frz, tc_est, hc, tc)
        ! Solve the finite-area combustor (FAC) rocket problem

        ! Arguments
        class(RocketSolver) :: self
        type(RocketSolution), intent(inout) :: soln
        real(dp), intent(in) :: reactant_weights(:)
        real(dp), intent(in) :: pc                   ! Injector pressure [Pa]
        real(dp), intent(in) :: pi_p(:)              ! Ratio of chamber pressure to exit pressure [unitless]
        real(dp), intent(in), optional :: subar(:)   ! Subsonic area ratio [unitless]
        real(dp), intent(in), optional :: supar(:)   ! Supersonic area ratio [unitless]
        real(dp), intent(in), optional :: ac_at      ! Contraction ratio: ratio of finite chamber area to throat area, Ac/At (FAC only) [unitless]
        real(dp), intent(in), optional :: mdot       ! Mass flow rate per chamber area (FAC only) [(kg/s)/m^2]
        integer,  intent(in), optional :: n_frz      ! Station where the composition should be frozen
        real(dp), intent(in), optional :: tc_est     ! Initial chamber temperature estimate [K]
        real(dp), intent(in), optional :: hc         ! Assigned chamber enthalpy [unitless]
        real(dp), intent(in), optional :: tc         ! Assigned chamber temperature [K]

        ! Locals
        integer :: i, j, k                   ! Loop index
        integer :: idx                       ! Station index
        integer :: n_frz_                    ! Temporary variable for frozen index
        integer :: num_pts                   ! Total number of evaluation points
        integer, parameter :: max_iter_area = 10  ! Maximum number of iterations for exit condition using area ratio
        real(dp), parameter :: area_tol = 4.0d-5  ! Area-ratio convergence tolerance
        character(len=2) :: prob_type        ! Equilibrium problem type
        real(dp) :: state1                   ! Chamber temperature or enthalpy, or entropy at other stations
        real(dp) :: S_ref                    ! Reference entropy
        real(dp), parameter :: c_tol = 2.d-5 ! Tolerance for chamber conditions convergence
        real(dp) :: gamma_s                  ! Temp variable for isentropic exponent gamma_s
        real(dp) :: usq, asq                 ! velocity squared; sonic velocity squared
        real(dp) :: h_inj, h                 ! Enthalpy at injector; enthalpy at any other station (temporary)
        real(dp) :: ln_pinf_pt               ! Temporary variable for ln(Pinf/Pt)
        real(dp) :: ln_pinf_pc               ! Temporary variable for ln(Pinf/Pc)
        real(dp) :: dln_pinf_pc_dln_acat     ! Partial derivative ∂ln(Pinf/Pc)/∂ln(Ac/At)
        real(dp) :: dln_pinf_pc              ! Update to ln(Pinf/Pc)
        real(dp) :: ac_at_                   ! Contraction ratio: ratio of finite chamber area to throat area, Ac/At (FAC only) [unitless]
        real(dp) :: awt                      ! Mass flow per area in the throat
        real(dp) :: p_inf                    ! Pressure at infinity
        real(dp) :: p_inj_check              ! Computed value of injector pressure (used to check convergence)
        real(dp) :: volume                   ! Volume temporary variable
        logical :: use_acat                  ! Flag to use the contraction ratio for chamber solution (if false, use mdot)
        logical :: frozen                    ! Flag to determine if frozen is used
        real(dp) :: acatsv, pratsv, mat, prat, pjrat, pr, pracat

        ! Index of solution:
        ! 1: injector
        ! 2: infinity
        ! 3: combustor end
        ! 4: throat
        ! 5+: exit

        call log_debug("Starting rocket FAC solve")

        ! Set the total number of evaluation points
        num_pts = 4 + size(pi_p)
        if (present(subar)) num_pts = num_pts + size(subar)
        if (present(supar)) num_pts = num_pts + size(supar)

        ! Initialize the solution variables
        soln = RocketSolution(self, num_pts=num_pts)
        call log_debug("Initialized RocketSolution")
        soln%throat_idx = 4
        soln%converged = .true.

        ! Set a flag to determine whether to use the contraction ratio or mdot for the chamber solution
        if (present(ac_at)) then
            use_acat = .true.
            ac_at_ = ac_at
        else if (present(mdot)) then
            use_acat = .false.
        else
            call abort("Either Ac/At or mdot/At must be specified for FAC rocket problem")
        end if

        ! Set the equilibrium problem type
        if (present(hc)) then
            prob_type = "hp"
            state1 = hc
            h_inj = hc
        else
            prob_type = "tp"
            state1 = tc
        end if

        ! Set the frozen composition flag
        if (present(n_frz)) then
            if (n_frz > 0 .and. n_frz <= num_pts) then
                frozen = .true.
                n_frz_ = n_frz
            else
                frozen = .false.
                n_frz_ = 100  ! Set to a large number to ensure it is not used
                if (n_frz /= 0) call log_info('RocketSolver: WARNING!!  INVALID FROZEN INDEX.  FROZEN MODE WILL NOT BE USED')
            end if
        else
            frozen = .false.
            n_frz_ = 100  ! Set to a large number to ensure it is not used
        end if

        ! -----------------------------------------------
        ! Injector conditions
        ! -----------------------------------------------
        idx = 1
        soln%station(idx) = "injector"

        if (present(tc_est)) then
            soln%eq_soln(idx) = EqSolution(self%eq_solver, T_init=tc_est)
            soln%eq_soln(2) = EqSolution(self%eq_solver, T_init=tc_est)  ! Initialize soln at infinity too
        else
            soln%eq_soln(idx) = EqSolution(self%eq_solver)
            soln%eq_soln(2) = EqSolution(self%eq_solver)  ! Initialize soln at infinity too
        end if
        call log_debug("Starting injector eqsolve")
        call self%eq_solver%solve(soln%eq_soln(idx), prob_type, state1, pc, reactant_weights, partials=soln%eq_partials(idx))

        ! Set the states
        soln%pressure(idx) = pc
        soln%mach(idx) = 0.0d0
        gamma_s = soln%eq_partials(idx)%gamma_s  ! Save gamma_s as the initial value for the next station
        soln%gamma_s(1) = gamma_s
        soln%v_sonic(idx) = sqrt(soln%eq_soln(idx)%n*R*soln%gamma_s(idx)*soln%eq_soln(idx)%T)

        ! Save the injector enthalpy to compute conditions at infinity, which is isenthropic with the injector
        if (prob_type == "tp") h_inj = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%enthalpy)*soln%eq_soln(idx)%T

        ! Initial estimate of pressure at infinity (Eq. 6.30)
        if (use_acat .eqv. .false.) ac_at_ = 2.0d0
        do i = 1, 4
            prat = (1.0257d0 - 1.2318d0*ac_at_)/(1.0d0 - 1.26505d0*ac_at_)
            p_inf = pc*prat
            if (use_acat) exit
            ac_at_ = 1.d5*p_inf/(mdot*2350.d0)
            if (ac_at_ < 1) then
                call abort("Input value of mdot is too large, resulting in a contraction ratio less than 1.")
            end if
        end do

        ! P_c initial guess
        soln%pressure(3) = p_inf

        ! Iterate until combustor conditions converge
        do i = 1, 4 ! Max outer iterations

            soln%pressure(2) = p_inf

            ! Solve conditions at infinity one time
            call self%eq_solver%solve(soln%eq_soln(2), prob_type, state1, p_inf, reactant_weights, partials=soln%eq_partials(2))

            ! Save the reference enthalpy at infinity
            S_ref = soln%eq_soln(2)%calc_entropy_sum(self%eq_solver)  ! Entropy at infinity

            ! -----------------------------------------------
            ! Throat conditions
            ! -----------------------------------------------
            idx = 4
            ! Legacy SETEN-style initialization for FAC throat:
            ! seed the throat solve from point 2 (infinity) by saving/using it.
            soln%eq_soln(3) = soln%eq_soln(2)
            soln%i_save = -2
            call self%solve_throat(soln, idx, p_inf, h_inj, S_ref, reactant_weights, awt)

            ! -----------------------------------------------
            ! Solve combustor once and check convergence
            ! -----------------------------------------------
            idx = 3

            ! Update initial guess of Pc
            if (i == 1) then
                ! Set equilibrium initial guess
                soln%eq_soln(idx) = EqSolution(self%eq_solver, T_init=soln%eq_soln(1)%T, nj_init=soln%eq_soln(1)%nj)

                ln_pinf_pt = log(soln%pressure(2)/soln%pressure(4))
                ln_pinf_pc = ln_pinf_pt/(ac_at_ + 10.587d0*(log(ac_at_)**3.0d0) + 9.454d0*log(ac_at_))
                if (ac_at_ < 1.09d0) then
                    ln_pinf_pc = 0.9d0*ln_pinf_pc
                else if (ac_at_ > 10.0d0) then
                    ln_pinf_pc = ln_pinf_pc/ac_at_
                end if
                ! Update the pressure
                soln%pressure(idx) = p_inf/exp(ln_pinf_pc)
            end if

            ! -----------------------------------------------
            ! Iterate at the combustor end until convergence for assigned area ratio
            ! -----------------------------------------------
            do j = 1, max_iter_area

                ! Solve equilibrium at combustor end (isentropic with infinity)
                call self%eq_solver%solve(soln%eq_soln(idx), "sp", S_ref, soln%pressure(idx), reactant_weights, &
                    partials=soln%eq_partials(idx))

                ! Compute combustor properties
                h = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%enthalpy)*soln%eq_soln(idx)%T
                gamma_s = soln%eq_partials(idx)%gamma_s
                asq = soln%eq_soln(idx)%n*R*gamma_s*soln%eq_soln(idx)%T
                usq = 2.0d0*(h_inj-h)*R
                soln%v_sonic(idx) = sqrt(asq)
                soln%mach(idx) = sqrt(usq/asq)

                soln%ae_at(idx) = soln%eq_soln(idx)%n*soln%eq_soln(idx)%T/(soln%pressure(idx)*sqrt(usq)*awt)

                ! Compute updated pressure ratio estimate
                dln_pinf_pc_dln_acat = gamma_s*usq/(usq - asq)  ! (Eq. 6.23)
                dln_pinf_pc = dln_pinf_pc_dln_acat*(log(ac_at_) - log(soln%ae_at(idx)))

                ! Convergence test for assigned area ratio:
                ! relative Ac/At error OR small pressure-ratio update.
                if (abs(soln%ae_at(idx)-ac_at_)/ac_at_ <= area_tol) exit
                if (abs(dln_pinf_pc) < area_tol) exit
                ln_pinf_pc = ln_pinf_pc + dln_pinf_pc  ! If not converged, update estimate
                if (ln_pinf_pc < 0.0d0) ln_pinf_pc = 0.000001d0

                ! Update the pressure
                soln%pressure(idx) = p_inf/exp(ln_pinf_pc)

            end do

            volume = R * soln%eq_soln(idx)%n * soln%eq_soln(idx)%T / soln%pressure(idx)
            p_inj_check = soln%pressure(3) + usq/volume

            if (use_acat) then

                prat = (p_inj_check/pc)

                ! Check combustor convergence (Eq. 6.29)
                if (abs(pc - p_inj_check)/p_inj_check <= c_tol) exit

                ! Update estimate
                p_inf = p_inf*prat

            else  ! Use mdot
                acatsv = soln%ae_at(idx)
                pratsv = prat
                mat = 1.d5/(awt*R)
                ac_at_ = mat/mdot
                prat = (1.0257d0 - 1.2318d0*ac_at_)/(1.0d0 - 1.26505d0*ac_at_)

                ! Check combustor convergence (Eq. 6.29)
                if (abs(pc - p_inj_check)/p_inj_check <= c_tol) exit

                pjrat = p_inj_check/pc
                do k = 1, 2
                    pracat = pratsv/prat
                    pr = pjrat*pracat
                    p_inf = p_inf/pr
                    ac_at_ = ac_at_/pr
                    soln%ae_at(idx) = ac_at_
                    pratsv = prat
                    pjrat = 1.0d0
                    prat = (1.0257d0 - 1.2318d0*ac_at_)/(1.0d0 - 1.26505d0*ac_at_)
                end do
            end if

        end do

        ! Recompute throat under frozen assumptions after chamber convergence,
        ! so frozen expansion uses a consistent throat state.
        idx = 4
        if (frozen .and. idx > n_frz_) then
            call self%solve_throat_frozen(soln, idx, n_frz_, p_inf, h_inj, awt)
        end if
        if (.not. soln%converged) return

        ! -----------------------------------------------
        ! Exit conditions: pressure ratio
        ! -----------------------------------------------
        idx = 5

        if (frozen .and. idx > n_frz_) then
            call self%solve_pi_p_frozen(soln, idx, n_frz_, pc, pi_p, h_inj, 2)
        else
            call self%solve_pi_p(soln, idx, pc, pi_p, h_inj, S_ref, reactant_weights)
        end if
        if (.not. soln%converged) return

        ! -----------------------------------------------
        ! Exit conditions: subsonic area ratio
        ! -----------------------------------------------

        if (present(subar)) then

            if (frozen .and. n_frz_ > 1) then
                call log_info('RocketSolver: WARNING!!  FREEZING IS NOT ALLOWED AT A SUBSONIC PRESSURE RATIO')
            else
                ln_pinf_pt = log(soln%pressure(2)/soln%pressure(4))
                call self%solve_subar(soln, idx, soln%pressure(2), subar, h_inj, S_ref, reactant_weights, &
                    idx-1, 2, ln_pinf_pt, awt)
            end if

        end if

        ! -----------------------------------------------
        ! Exit conditions: supersonic area ratio
        ! -----------------------------------------------

        if (present(supar)) then
            ln_pinf_pt = log(soln%pressure(2)/soln%pressure(4))

            if (frozen .and. idx > n_frz_) then
                call self%solve_supar_frozen(soln, idx, n_frz_, soln%pressure(2), supar, h_inj, reactant_weights, &
                    idx-1, ln_pinf_pt, awt)
            else
                call self%solve_supar(soln, idx, soln%pressure(2), supar, h_inj, S_ref, reactant_weights, &
                    ln_pinf_pt, awt)
            end if
            if (.not. soln%converged) return
        end if

        ! Omitted frozen schedule points do not consume output indices.
        soln%num_pts = idx - 1

        ! Compute performance parameters
        call self%post_process(soln, .true.)

    end subroutine

    function RocketSolver_solve(self, reactant_weights, pc, pi_p, fac, subar, supar, mdot, ac_at, n_frz, tc_est, hc, tc) &
        result(soln)
        ! Solve the rocket proclem

        ! Arguments
        class(RocketSolver), intent(in) :: self
        real(dp), intent(in) :: reactant_weights(:)  !
        real(dp), intent(in) :: pc                   ! Chamber pressure [bar]
        real(dp), intent(in) :: pi_p(:)              ! Ratio of chamber pressure to exit pressure [unitless]
        logical,  intent(in), optional :: fac        ! Finite-area combustor flag
        real(dp), intent(in), optional :: subar(:)   ! Subsonic area ratio [unitless]
        real(dp), intent(in), optional :: supar(:)   ! Supersonic area ratio [unitless]
        real(dp), intent(in), optional :: mdot       ! Mass flow rate per chamber area (FAC only) [(kg/s)/m^2]
        real(dp), intent(in), optional :: ac_at      ! Contraction ratio: ratio of finite chamber area to throat area, Ac/At (FAC only) [unitless]
        integer,  intent(in), optional :: n_frz      ! Station where the composition should be frozen
        real(dp), intent(in), optional :: tc_est     ! Initial chamber temperature estimate [K]
        real(dp), intent(in), optional :: hc         ! Assigned chamber enthalpy [unitless]
        real(dp), intent(in), optional :: tc         ! Assigned chamber temperature [K]

        ! Result
        type(RocketSolution) :: soln

        ! Locals
        logical :: fac_

        ! Set defaults for optional arguments
        fac_ = .false.
        if (present(fac)) fac_ = fac

        call log_info("Starting rocket solve")

        ! Call either the FAC or IAC solver
        if (fac_) then
            call self%solve_fac(soln, reactant_weights, pc, pi_p, subar, supar, ac_at, mdot, n_frz, tc_est, hc, tc)
        else  ! IAC
            call self%solve_iac(soln, reactant_weights, pc, pi_p, subar=subar, supar=supar, n_frz=n_frz, &
                tc_est=tc_est, hc=hc, tc=tc)
        end if

    end function

    subroutine RocketSolver_post_process(self, soln, fac)
        ! Arguments
        class(RocketSolver), intent(in) :: self
        type(RocketSolution), intent(inout) :: soln
        logical, intent(in) :: fac

        ! Locals
        integer :: i, start
        real(dp) :: aw, h0, hi

        start = 2
        if (fac) start = 3

        h0 = dot_product(soln%eq_soln(1)%nj, soln%eq_soln(1)%thermo%enthalpy)*soln%eq_soln(1)%T

        do i = start, soln%num_pts

            ! Isp
            hi = dot_product(soln%eq_soln(i)%nj, soln%eq_soln(i)%thermo%enthalpy)*soln%eq_soln(i)%T
            soln%i_sp(i) = (2.0d0*R*(h0 - hi))**0.5d0 ! /agv

            ! Common term
            aw = R*soln%eq_soln(i)%T/(soln%pressure(i)*soln%i_sp(i)/soln%eq_soln(i)%n)  ! / agv**2

            ! c*
            if (i == soln%throat_idx) then
                soln%c_star(:) = aw*soln%pressure(1)
                if (fac) then
                    soln%c_star(:) = aw*soln%pressure(2)
                end if
            end if

            ! Isp_vac
            soln%i_vac(i) = soln%i_sp(i) + soln%pressure(i)*aw

        end do

        ! CF needs c* so compute it last
        do i = start, soln%num_pts
            ! CF
            soln%cf(i) = soln%i_sp(i)/soln%c_star(i) ! *gc
        end do

    end subroutine

    subroutine RocketSolver_set_init_state(self, soln, idx)
        ! Set the initial state for the equilibrium solver

        ! Arguments
        class(RocketSolver), intent(in) :: self
        type(RocketSolution), intent(inout) :: soln
        integer, intent(in) :: idx  ! New index

        ! Locals
        real(dp), allocatable :: nj_init(:)
        integer :: ng, nc, i
        logical :: const_t

        allocate(nj_init(self%eq_solver%num_products))

        ! Shorthand
        ng = self%eq_solver%num_gas
        nc = self%eq_solver%num_condensed
        const_t = soln%eq_soln(idx)%constraints%is_constant_temperature()

        ! Set the initial temperature and composition
        if (soln%i_save < 0) then
            soln%i_save = -soln%i_save
            soln%eq_soln_save%T = soln%eq_soln(soln%i_save)%T
            soln%eq_soln_save%nj = soln%eq_soln(idx-1)%nj
            soln%eq_soln_save%ln_nj = soln%eq_soln(idx-1)%ln_nj
            soln%eq_soln_save%n = soln%eq_soln(idx-1)%n

            soln%eq_soln(idx)%nj(:ng) = soln%eq_soln(idx-1)%nj(:ng)
            soln%eq_soln(idx)%ln_nj = soln%eq_soln(idx-1)%ln_nj
            soln%eq_soln(idx)%nj(ng+1:) = soln%eq_soln(soln%i_save)%nj(ng+1:)
            soln%eq_soln(idx)%is_active = soln%eq_soln(idx-1)%is_active
            do i = 1, nc
                if (.not. soln%eq_soln(soln%i_save)%is_active(i)) cycle
                if (i == soln%eq_soln(idx-1)%j_liq) then
                    soln%eq_soln(idx)%nj(ng+soln%eq_soln(idx-1)%j_sol) = &
                        soln%eq_soln(soln%i_save)%nj(ng+soln%eq_soln(idx-1)%j_liq) + &
                        soln%eq_soln(soln%i_save)%nj(ng+soln%eq_soln(idx-1)%j_sol)
                    soln%eq_soln(idx)%nj(ng+soln%eq_soln(idx-1)%j_liq) = 0.0d0
                    soln%eq_soln(idx)%j_liq = 0
                    soln%eq_soln(idx)%j_sol = 0
                    soln%eq_soln_save%T = soln%eq_soln_save%T - 5.0d0
                    soln%eq_soln_save%nj(ng+i) = 0.0d0
                    soln%eq_soln(idx)%is_active(i) = .false.
                end if
            end do
            soln%eq_soln(idx)%T = soln%eq_soln(soln%i_save)%T
        else if (soln%i_save > 0) then
            ! Use composition from i_save; everything else from previous point
            soln%eq_soln(idx)%T = soln%eq_soln(idx-1)%T
            soln%eq_soln(idx)%nj = soln%eq_soln(soln%i_save)%nj !eq_soln(idx-1)%nj
            soln%eq_soln(idx)%ln_nj = soln%eq_soln(idx-1)%ln_nj !eq_soln(idx-1)%ln_nj
            soln%eq_soln(idx)%n = soln%eq_soln(idx-1)%n !eq_soln(idx-1)%n
            soln%eq_soln(idx)%is_active = soln%eq_soln(idx-1)%is_active !eq_soln(idx-1)%is_active
            soln%eq_soln(idx)%j_liq = soln%eq_soln(idx-1)%j_liq
            soln%eq_soln(idx)%j_sol = soln%eq_soln(idx-1)%j_sol
        else  ! i_save = 0
            soln%eq_soln(idx)%j_liq = 0
            soln%eq_soln(idx)%j_sol = 0
            soln%eq_soln(idx)%n = soln%eq_soln_save%n
            soln%eq_soln(idx)%ln_nj = soln%eq_soln_save%ln_nj
            soln%eq_soln(idx)%nj(ng+1:) = soln%eq_soln_save%nj(ng+1:)
            soln%eq_soln(idx)%is_active = soln%eq_soln_save%is_active
            do i = 1, ng
                if (soln%eq_soln(idx)%ln_nj(i) - log(soln%eq_soln(idx)%n) + 18.5d0 > 0.0d0) then
                    soln%eq_soln(idx)%nj(i) = exp(soln%eq_soln(idx)%ln_nj(i))
                else
                    soln%eq_soln(idx)%nj(i) = 0.0d0
                end if
            end do
            if (.not. const_t) soln%eq_soln(idx)%T = soln%eq_soln_save%T
        end if


    end subroutine

    !-----------------------------------------------------------------------
    ! RocketSolution
    !-----------------------------------------------------------------------
    function RocketSolution_init(solver, num_pts) result(self)
        type(RocketSolution) :: self
        type(RocketSolver), intent(in) :: solver
        integer, optional :: num_pts
        integer :: i

        if (present(num_pts)) then
            self%num_pts = num_pts
        else
            self%num_pts = 10  ! Default oversized for most cases
        end if

        ! Allocate data structures
        allocate(self%station(self%num_pts))
        allocate(self%eq_soln(self%num_pts))
        allocate(self%eq_partials(self%num_pts))
        allocate(self%pressure(self%num_pts), source=0.0d0)
        allocate(self%mach(self%num_pts), source=0.0d0)
        allocate(self%gamma_s(self%num_pts), source=0.0d0)
        allocate(self%v_sonic(self%num_pts), source=0.0d0)
        allocate(self%ae_at(self%num_pts), source=0.0d0)
        allocate(self%c_star(self%num_pts), source=0.0d0)
        allocate(self%cf(self%num_pts), source=0.0d0)
        allocate(self%i_sp(self%num_pts), source=0.0d0)
        allocate(self%i_vac(self%num_pts), source=0.0d0)

        ! Initalize the station names as blank
        do i = 1,self%num_pts
            self%station(i) = "        "
        end do

        ! Initialize the saved EqSolution for future initial guesses
        self%i_save = 0
        self%eq_soln_save = EqSolution(solver%eq_solver)

    end function

end module
