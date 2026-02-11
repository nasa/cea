module cea_equilibrium
    !! Equilibrium Solver Module

    use cea_param, only: dp, empty_dp, R=>gas_constant, &
                         snl=>species_name_len, &
                         enl=>element_name_len, &
                         Avgdr=>avogadro, &
                         Boltz=>boltzmann, &
                         pi
    use cea_mixture, only: Mixture, MixtureThermo
    use cea_transport, only: TransportDB, get_mixture_transport
    use fb_findloc, only: findloc
    use fb_utils
    implicit none

    type :: EqSolver
        !! Equilibrium Solver Type

        ! Sizing variables
        integer :: num_reactants = 0
            !! Number of species in the product mixture
        integer :: num_products = 0
            !! Number of species in the product mixture
        integer :: num_gas = 0
            !! Number of gas species in the product mixture
        integer :: num_condensed = 0
            !! Number of condensed species in the product mixture
        integer :: num_elements = 0
            !! Number of elements in the product mixture
        integer :: max_equations = 0
            !! Maximum possible number of equations in the matrix system

        ! Thermodynamic database
        type(Mixture) :: reactants
            !! Reactant thermodynamic database
        type(Mixture) :: products
            !! Product thermodynamic database
        type(MixtureThermo) :: thermo
            !! Product heat capacity, enthalpy, entropy

        ! Transport database
        type(TransportDB) :: transport_db
            !! Transport properties database

        ! Options
        logical :: ions = .false.
            !! Flag if ions should be included
        logical :: active_ions = .true.
            !! Flag if ions are currently active
        integer :: reduced_elements = 0
            !! Number of temporarily reduced element equations in singular recovery
        logical :: transport = .false.
            !! Flag if transport properties should be computed
        character(snl), allocatable :: insert(:)

        ! Solver parameters
        ! TODO: Give these better names
        real(dp) :: size  = 18.420681d0
            !! log(1.d-8)
        real(dp) :: xsize = 25.328436d0
            !! log(1.d-11)
        real(dp) :: esize = 32.236191d0
            !! log(1.d-14)
        real(dp) :: tsize = 18.420681d0
            !! log(1.d-8)
        real(dp) :: trace = 0.0d0
            !! Threshold for trace species
        real(dp) :: log_min = -87.0d0
            !! After convergence, convert all log(n) to n where log(n) > LOG_MIN
        integer  :: max_iterations = 50
            !! Maximum number of iterations
        integer  :: max_converged  = 10
            !! Maximum number of times the problem can converge without establishing a set of condensed species
        real(dp) :: T_min = 160.0d0
            !! Minimum mixture temperature (K)
        real(dp) :: T_max = 6600.0d0
            !! Maximum mixture temperature (K)

    contains

        procedure :: num_active_elements => EqSolver_num_active_elements
        procedure :: compute_damped_update_factor => EqSolver_compute_damped_update_factor
        procedure :: get_solution_vars => EqSolver_get_solution_vars
        procedure :: update_solution => EqSolver_update_solution
        procedure :: check_convergence => EqSolver_check_convergence
        procedure :: check_condensed_phases =>EqSolver_check_condensed_phases
        procedure :: test_condensed => EqSolver_test_condensed
        procedure :: correct_singular => EqSolver_correct_singular
        procedure :: assemble_matrix => EqSolver_assemble_matrix
        procedure :: post_process => EqSolver_post_process
        procedure :: solve => EqSolver_solve

    end type
    interface EqSolver
        module procedure :: EqSolver_init
    end interface

    type :: EqConstraints
        !! Equilibrium Constraints Type
        character(2) :: type
            !! Constraint type, e.g. 'tp', 'hp', etc.
        real(dp) :: state1
            !! First thermodynamic state constraint (t, h, s, or u)
        real(dp) :: state2
            !! Second thermodynamic state constraint (p or v)
        real(dp), allocatable :: b0(:)
            !! Element abundance constraints
    contains
        procedure :: set => EqConstraints_set
        procedure :: is_constant_temperature => EqConstraints_is_constant_temperature
        procedure :: is_constant_enthalpy    => EqConstraints_is_constant_enthalpy
        procedure :: is_constant_energy      => EqConstraints_is_constant_energy
        procedure :: is_constant_entropy     => EqConstraints_is_constant_entropy
        procedure :: is_constant_pressure    => EqConstraints_is_constant_pressure
        procedure :: is_constant_volume      => EqConstraints_is_constant_volume
    end type
    interface EqConstraints
        module procedure :: EqConstraints_init
        module procedure :: EqConstraints_alloc
    end interface

    type :: EqSolution
        !! Equilibrium Solution Type

        ! Mixture data
        real(dp) :: T
            !! Mixture temperautre
        real(dp), allocatable :: nj(:)
            !! Product concentrations [kg-mol per kg mixture]
        real(dp) :: n
            !! Total moles of mixture
        type(MixtureThermo) :: thermo
            !! Product heat capacity, enthalpy, entropy
        real(dp), allocatable :: ln_nj(:)
            !! Log(nj)

        ! Solution update variables
        real(dp), allocatable :: pi(:)
            !! Modified Lagrance multipliers, 𝛑 = -λ/RT
        real(dp), allocatable :: pi_prev(:)
            !! pi variable from previous iteration
        real(dp) :: pi_e = 0.0d0
            !! pi variable for ionized species
        real(dp) :: dpi_e = 0.0d0
            !! 𝛥pi variable for ionized species
        real(dp), allocatable :: dln_nj(:)
            !! 𝛥ln(nj) (gas-phase only)
        real(dp), allocatable :: dnj_c(:)
            !! Change in condensed species concentrations, 𝛥nj_c
        real(dp) :: dln_n
            !! Change in total log moles, 𝛥ln(n)
        real(dp) :: dln_T
            !! Change in log temperature, 𝛥ln(T)

        ! Algorithm workspace
        real(dp), allocatable :: G(:,:)
            !! Augmented Newton iteration matrix (oversized!)
        type(EqConstraints) :: constraints
            !! State and element constraints
        logical, allocatable :: is_active(:)
            !! True if condensed species included in G
        integer, allocatable :: active_rank(:)
            !! Active condensed ordering rank (1 = newest/front), 0 = inactive
        integer :: j_liq = 0
            !! Index of liquid phase if two phases are present for same species
        integer :: j_sol = 0
            !! Index of solid phase if two phases are present for same species
        integer :: j_switch = 0
            !! Index of the condensed species that was last removed from the active set
        integer :: last_cond_idx = 0
            !! Index of the condensed species that was added most recently

        ! Convenience variables
        logical :: gas_converged         = .false.
            !! Flag if gas species concentrations have converged
        logical :: condensed_converged   = .false.
            !! Flag if condensed species concentrations have converged
        logical :: moles_converged       = .false.
            !! Flag if total moles have converged
        logical :: element_converged     = .false.
            !! Flag if element concentrations have converged
        logical :: temperature_converged = .false.
            !! Flag if temperature has converged
        logical :: entropy_converged     = .false.
            !! Flag if entropy has converged
        logical :: pi_converged          = .false.
            !! Flag if modified Lagrance multipliers have converged
        logical :: ions_converged        = .false.
            !! Flag if ionized species have converged
        logical :: converged             = .false.
            !! Flag if the solution has converged
        integer :: times_converged       = 0
            !! Number of times the solution has converged without establishing a set of condensed species

        ! Other solution variables
        real(dp), allocatable :: mole_fractions(:)
            !! Mole fractions of products
        real(dp), allocatable :: mass_fractions(:)
            !! Mass fractions of products

        ! Mixture properties
        real(dp) :: density = 0.0d0
            !! TODO: Mixture density (kg/m^3)
        real(dp) :: pressure = 0.0d0
            !! Mixture pressure (bar)
        real(dp) :: volume = 0.0d0
            !! Mixture specific volume (m^3/kg)
        real(dp) :: M = 0.0d0
            !! Molecular weight, (1/n) (Eq. 2.3a/b)
        real(dp) :: MW = 0.0d0
            !! Molecular weight (Eq. 2.4a/b)
        real(dp) :: enthalpy = 0.0d0
            !! Mixture enthalpy (kJ/kg)
        real(dp) :: energy = 0.0d0
            !! Mixture internal energy (kJ/kg)
        real(dp) :: gibbs_energy = 0.0d0
            !! Mixture Gibb's energy (kJ/kg)
        real(dp) :: entropy = 0.0d0
            !! Mixture entropy (kJ/kg-K)

        ! Transport properties
        real(dp) :: gamma_s = 0.0d0
            !! Isentropic exponent (Eq. 2.71/2.73)
        real(dp) :: viscosity = 0.0d0
            !! Viscosity (millipose)
        real(dp) :: cp_fr = 0.0d0
            !! Heat capacity at constant pressure, frozen (kJ/kg-K)
        real(dp) :: cp_eq = 0.0d0
            !! Heat capacity at constant pressure, equilibrium (kJ/kg-K)
        real(dp) :: cv_fr = 0.0d0
            !! Heat capacity at constant volume, frozen (kJ/kg-K)
        real(dp) :: cv_eq = 0.0d0
            !! Heat capacity at constant volume, equilibrium (kJ/kg-K)
        real(dp) :: conductivity_fr = 0.0d0
            !! Thermal conductivity, frozen (mW/cm-K)
        real(dp) :: conductivity_eq = 0.0d0
            !! Thermal conductivity, equilibrium (mW/cm-K)
        real(dp) :: Pr_fr = 0.0d0
            !! Prandtl number, frozen (unitless)
        real(dp) :: Pr_eq = 0.0d0
            !! Prandtl number, equilibrium (unitless)

        ! Note:
        !   The augmented iteration matrix is an n-by-n+1 matrix that appends
        !   the function residual vector as an extra column after the Jacobian.
        !   Thus, given F(x)=0 with corresponding Newton iteration:
        !
        !      x(n+1) = x(n) + inv(dF/dx)*F(x(n))
        !
        !   The augmented iteration matrix is defined as:
        !
        !      G(x) = [ dF/dx, F(x) ]
        !
        !   However, the number of equations that appears in this matrix depends
        !   on the type of equilibrium problem solved and the number of condensed
        !   species that are "active", i.e. assumed to be present in the mixture.
        !   The CEA solution algorithm requires activating/deactivating condensed
        !   species as the solution process proceeds. Thus, to have avoid having
        !   to reallocate this matrix repeatedly during soluiton, we over-size it
        !   and then only work with the "active" subset of the G array.
    contains
        procedure :: set_nj => EqSolution_set_nj
        procedure :: activate_condensed => EqSolution_activate_condensed
        procedure :: activate_condensed_front => EqSolution_activate_condensed_front
        procedure :: deactivate_condensed => EqSolution_deactivate_condensed
        procedure :: replace_active_condensed => EqSolution_replace_active_condensed
        procedure :: active_condensed_indices => EqSolution_active_condensed_indices
        procedure :: num_equations => EqSolution_num_equations
        procedure :: calc_pressure => EqSolution_calc_pressure
        procedure :: calc_volume => EqSolution_calc_volume
        procedure :: calc_entropy_sum => EqSolution_calc_entropy_sum
    end type
    interface EqSolution
        module procedure :: EqSolution_init
    end interface

    type :: EqPartials
        !! Equilibrium Partial Derivatives Type
        real(dp), allocatable :: dpi_dlnT(:)
            !! Partial derivative of 𝛑 wrt ln(T) (const P)
        real(dp), allocatable :: dnc_dlnT(:)
            !! Partial derivative of nc wrt ln(T) (const P)
        real(dp)              :: dn_dlnT
            !! Partial derivative of n wrt ln(T) (const P)
        real(dp)              :: dlnV_dlnT
            !! Partial derivative of ln(V) wrt ln(T) (const P)

        real(dp), allocatable :: dpi_dlnP(:)
            !! Partial derivative of 𝛑 wrt ln(P) (const T)
        real(dp), allocatable :: dnc_dlnP(:)
            !! Partial derivative of nc wrt ln(P) (const T)
        real(dp)              :: dn_dlnP
            !! Partial derivative of n wrt ln(P) (const T)
        real(dp)              :: dlnV_dlnP
            !! Partial derivative of ln(V) wrt ln(P) (const T)

        real(dp)              :: cp_eq
            !! Equilibrium heat capacity [J/kg-K] (Eq. 2.59)
        real(dp)              :: gamma_s
            !! Isentropic exponent (Eq. 2.71/2.73)
    contains
        procedure :: assemble_partials_matrix_const_p => EqPartials_assemble_partials_matrix_const_p
        procedure :: assemble_partials_matrix_const_t => EqPartials_assemble_partials_matrix_const_t
        procedure :: compute_partials => EqPartials_compute_partials
    end type
    interface EqPartials
        module procedure :: EqPartials_init
    end interface

contains

    !-----------------------------------------------------------------------
    ! EquilibriumSolver
    !-----------------------------------------------------------------------
    function EqSolver_init(products, reactants, trace, ions, all_transport, insert) result(self)
        type(EqSolver) :: self
        type(Mixture), intent(in) :: products
        type(Mixture), intent(in), optional :: reactants
        real(dp), intent(in), optional :: trace
        logical, intent(in), optional :: ions
        type(TransportDB), intent(in), optional :: all_transport
        character(*), intent(in), optional :: insert(:)  ! List of condensed species to insert
        integer :: i
        integer :: ngc_equiv

        ! Initialize reactant data
        self%products = products
        if (present(reactants)) then
            self%reactants = reactants
            if (reactants%num_species > products%num_species) then
                call abort("EqSolver_init: num_reactants > num_products. Arguments for products and reactants may be flipped.")
            end if
        else
            self%reactants = products
        end if
        call assert( &
            all(self%products%element_names == self%reactants%element_names), &
            'eqsolver_init: Element lists for products and reactants must match.' &
        )

        self%num_reactants = self%reactants%num_species
        self%num_products  = self%products%num_species
        self%num_gas       = self%products%num_gas
        self%num_condensed = self%products%num_condensed
        self%num_elements  = self%products%num_elements
        self%max_equations = self%num_elements + self%num_condensed + 2

        ! Optional argument handling
        if (present(trace)) self%trace = trace
        if (present(ions)) self%ions = ions
        if (present(all_transport)) self%transport = .true.

        ! Update size parameters
        if (self%trace > 0.0d0) then
            ! Match CEA2 scaling (maxitn = 50 + Ngc/2), where Ngc counts condensed
            ! species by temperature interval. Approximate Ngc from this data model.
            ngc_equiv = self%num_gas
            do i = 1, self%num_condensed
                ngc_equiv = ngc_equiv + max(1, self%products%species(self%num_gas+i)%num_intervals)
            end do
            self%max_iterations = 50 + ngc_equiv/2
            self%xsize = -log(self%trace)
            if (self%xsize < self%size) self%xsize = self%size + 0.1d0
        end if

        if (self%xsize > 80.0d0) self%xsize = 80.0d0

        ! Match legacy CEA2 behavior: always derive esize from xsize.
        self%esize = min(80.0d0, self%xsize + 6.90775528d0)

        ! Set the max number of times that the solution can converge without establishing a set of condensed species
        self%max_converged = 3*self%products%num_elements

        ! Initialize transport database
        if (self%transport) self%transport_db = get_mixture_transport(all_transport, products, ions=self%ions)

        ! Store the insert species
        if (present(insert)) then
            do i = 1, size(insert)
                if (len_trim(insert(i)) > snl) then
                    call abort('EqSolver_init: insert species name too long: '//trim(insert(i)))
                end if
            end do
            self%insert = insert
        end if

    end function

    function EqSolver_num_active_elements(self) result(ne)
        class(EqSolver), intent(in) :: self
        integer :: ne

        ne = self%num_elements - self%reduced_elements
        if (self%ions .and. .not. self%active_ions) ne = max(0, ne-1)
        ne = max(0, ne)
    end function

    subroutine EqSolver_swap_elements(self, soln, i, j)
        ! Swap two element equations/columns in the solver state.
        class(EqSolver), intent(inout), target :: self
        type(EqSolution), intent(inout), target :: soln
        integer, intent(in) :: i
        integer, intent(in) :: j
        real(dp) :: tmp_col(self%num_products)
        real(dp) :: tmp
        character(enl) :: tmp_name

        if (i == j) return
        if (i < 1 .or. i > self%num_elements .or. j < 1 .or. j > self%num_elements) then
            call abort('EqSolver_swap_elements: index out of bounds.')
        end if

        tmp_col = self%products%stoich_matrix(:, i)
        self%products%stoich_matrix(:, i) = self%products%stoich_matrix(:, j)
        self%products%stoich_matrix(:, j) = tmp_col

        tmp_name = self%products%element_names(i)
        self%products%element_names(i) = self%products%element_names(j)
        self%products%element_names(j) = tmp_name

        tmp = soln%constraints%b0(i)
        soln%constraints%b0(i) = soln%constraints%b0(j)
        soln%constraints%b0(j) = tmp

        tmp = soln%pi(i)
        soln%pi(i) = soln%pi(j)
        soln%pi(j) = tmp

        tmp = soln%pi_prev(i)
        soln%pi_prev(i) = soln%pi_prev(j)
        soln%pi_prev(j) = tmp
    end subroutine

    subroutine EqSolver_restore_reduced_elements(self, soln, num_swaps, swap_from, swap_to)
        ! Restore element ordering after solve-local component reduction.
        class(EqSolver), intent(inout), target :: self
        type(EqSolution), intent(inout), target :: soln
        integer, intent(in) :: num_swaps
        integer, intent(in) :: swap_from(:)
        integer, intent(in) :: swap_to(:)
        integer :: i

        if (num_swaps > 0) then
            do i = num_swaps, 1, -1
                call EqSolver_swap_elements(self, soln, swap_from(i), swap_to(i))
            end do
        end if
        self%reduced_elements = 0
    end subroutine

    function EqSolver_compute_damped_update_factor(self, soln) result(lambda)
        ! Compute the damped update factor, lambda, for the Newton solver

        ! Arguments
        class(EqSolver), intent(in), target :: self
        type(EqSolution), intent(inout), target :: soln

        ! Result
        real(dp) :: lambda

        ! Locals
        integer :: ng                                 ! Number of gas species
        integer :: ne                                 ! Number of elements
        real(dp), pointer :: nj_g(:)                  ! Total/gas species concentrations [kmol-per-kg]
        real(dp) :: n                                 ! Total moles of mixture
        real(dp), pointer :: ln_nj(:)                 ! Log of the product concentrations
        real(dp), pointer :: dln_nj(:)                ! 𝛥ln(nj)
        integer :: i                                  ! Indices
        logical :: const_t                            ! Flag that is true if problem is constant temperature
        type(EqConstraints), pointer :: cons          ! Abbreviation for soln%constraints
        real(dp) :: dln_n                             ! 𝛥ln(n)
        real(dp) :: dln_T                             ! 𝛥ln(T)
        real(dp) :: lambda1, lambda2                  ! Candidate damped update factors
        real(dp) :: l1_denom, l2_denom, temp_l2       ! Temporary variables
        real(dp), parameter :: FACTOR = -9.2103404d0  ! log(1.d-4)

        ! Define shorthand
        ng = self%num_gas
        ne = self%num_active_elements()
        cons => soln%constraints
        ln_nj => soln%ln_nj
        dln_nj => soln%dln_nj
        n = soln%n
        const_t = cons%is_constant_temperature()
        dln_n = soln%dln_n
        dln_T = soln%dln_T

        ! Associate subarray pointers
        nj_g => soln%nj(:ng)

        ! Compute lambda1 (Eq. 3.1) and lambda2 (Eq. 3.2)
        ! TODO: what happens for a TV problem when these are both not defined?
        l1_denom = max(5.0d0*abs(dln_T), 5.0d0*abs(dln_n))

        lambda1 = 1.0d0
        lambda2 = 1.0d0
        do i = 1, ng
            if (dln_nj(i) > 0.0d0) then !.and. (nj_g(i) > 0.0d0)) then
                if (ln_nj(i) - log(n) + self%size <= 0.0d0) then
                    l2_denom = abs(dln_nj(i) - dln_n)
                    if (l2_denom >= (self%size + FACTOR)) then
                        temp_l2 = abs(FACTOR - ln_nj(i) + log(n))/l2_denom
                        lambda2 = min(lambda2, temp_l2)
                    end if
                else if (dln_nj(i) > l1_denom) then
                    l1_denom = dln_nj(i)
                end if
            end if
        end do
        if (l1_denom > 2.0d0) lambda1 = 2.0d0/l1_denom

        ! Compute lamba (Eq. 3.3)
        lambda = min(1.0d0, lambda1, lambda2)

    end function

    subroutine EqSolver_get_solution_vars(self, soln)
        ! Get the solution variables from the solution vector X

        ! Arguments
        class(EqSolver), target :: self
        type(EqSolution), intent(inout), target :: soln

        ! Locals
        integer :: i                          ! Indices
        integer :: ng                         ! Number of gas species
        integer :: ne                         ! Number of elements
        integer :: ne_full                    ! Total number of elements (including electron)
        integer :: na                         ! Number of active condensed species
        integer :: num_eqn                    ! Number of equations in the matrix system
        real(dp), pointer :: x(:)             ! Solution vector
        logical :: const_p, const_t           ! Flags enabling/disabling matrix equations
        type(EqConstraints), pointer :: cons  ! Abbreviation for soln%constraints
        real(dp) :: P                         ! Pressure state (bar)
        real(dp) :: mu_g(self%num_gas)        ! Gas phase chemical potentials [unitless]
        real(dp), pointer :: h_g(:)           ! Gas enthalpies [unitless]
        real(dp), pointer :: s_g(:)           ! Gas entropies [unitless]
        real(dp), pointer :: A_g(:,:), A(:,:) ! Gas/total stoichiometric matrices
        real(dp) :: n                         ! Total moles of mixture
        real(dp), pointer :: ln_nj(:)         ! Log of the product concentrations

        ! Define shorthand
        ng = self%num_gas
        ne = self%num_active_elements()
        ne_full = self%num_elements
        na = count(soln%is_active)
        num_eqn = soln%num_equations(self)
        cons => soln%constraints
        x => soln%G(:,num_eqn+1)
        ln_nj => soln%ln_nj
        n = soln%n
        const_p = cons%is_constant_pressure()
        const_t = cons%is_constant_temperature()

        ! Associate subarray pointers
        A_g => self%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A   => self%products%stoich_matrix(:,:)
        h_g => soln%thermo%enthalpy(:ng)
        s_g => soln%thermo%entropy(:ng)

        ! Set the solution variables
        x => soln%G(:num_eqn,num_eqn+1)

        ! Get the mixture pressure
        P = soln%calc_pressure()

        ! Compute gas phase chemical potentials
        mu_g = h_g - s_g + ln_nj + log(P/n)

        ! Get the pi variables
        soln%pi = 0.0d0
        if (ne > 0) soln%pi(:ne) = x(:ne)

        ! Check on removing ions: if all ionized species have a concentration of 0, remove ions entirely
        if (self%ions .and. self%active_ions .and. ne_full > 0) then
            self%active_ions = .false.
            do i = 1, self%num_products
                if (A(i, ne_full) /= 0.0d0) then
                    if (soln%nj(i) > 0.0d0) self%active_ions = .true.
                    if (soln%nj(i) > 0.0d0) exit  ! * I get a compile error when this is included with the above line *
                end if
            end do
            if (self%active_ions .and. soln%converged .and. .not. soln%ions_converged) then
                soln%pi_e = x(ne)
            else
                soln%pi_e = 0.0d0
            end if
        else
            soln%pi_e = 0.0d0
        end if

        ! Get the updates to the condensed species concentrations
        do i = 1, na
            soln%dnj_c(i) = x(ne+i)
        end do

        ! Get the update to the total log moles
        if (const_p) soln%dln_n = x(ne+na+1)

        ! Get the update to the log of temperature
        if (.not. const_t) then
            if (const_p) then
                soln%dln_T = x(ne+na+2)
            else
                soln%dln_T = x(ne+na+1)
            end if
        end if

        ! Get the update to the log of gas-species concentrations, 𝛥ln(nj)
        do i = 1, ng
            ! TODO: Skip the update here for species with removed elements

            soln%dln_nj(i) = -mu_g(i) + soln%dln_n + dot_product(A_g(i, :ne), soln%pi(:ne)) + soln%dln_T*h_g(i)
            if (.not. const_p) soln%dln_nj(i) = soln%dln_nj(i) - soln%dln_T

            ! Ionized species update
            if (self%ions .and. self%active_ions .and. ne_full > 0 .and. soln%pi_e /= 0.0d0) then
                soln%dln_nj(i) = soln%dln_nj(i) + A_g(i, ne_full)*soln%pi_e
            end if
        end do

    end subroutine

    subroutine EqSolver_update_solution(self, soln)
        ! Update the variables in the solution database using the solution vector X

        ! Arguments
        class(EqSolver), intent(in), target :: self
        type(EqSolution), intent(inout), target :: soln

        ! Locals
        integer  :: ng                        ! Number of gas species
        integer  :: nc                        ! Number of condensed species
        integer  :: na                        ! Number of active condensed species
        integer  :: ne                        ! Number of elements
        integer  :: ne_full                   ! Total number of elements (including electron)
        integer  :: num_eqn                   ! Number of equations in the matrix system
        real(dp), pointer :: nj_g(:), nj_c(:) ! Gas/condensed species concentrations [kmol-per-kg]
        real(dp) :: n                         ! Total moles of mixture
        real(dp), pointer :: ln_nj(:)         ! Log of the product concentrations
        real(dp) :: ln_n                      ! Log of the total concentration
        real(dp) :: T                         ! Temperature state
        real(dp), pointer :: h_g(:)           ! Gas enthalpies [unitless]
        real(dp), pointer :: s_g(:)           ! Gas entropies [unitless]
        real(dp), pointer :: A_g(:,:)         ! Gas stoichiometric matrices
        real(dp), pointer :: dln_nj(:)        ! 𝛥ln(nj)
        real(dp), pointer :: dnj_c(:)         ! 𝛥nj_c
        real(dp) :: dln_n                     ! 𝛥ln(n)
        real(dp) :: dln_T                     ! 𝛥ln(T)
        integer :: i, idx_c                   ! Indices
        integer, allocatable :: active_idx(:) ! Active condensed indices in legacy order
        logical :: const_p, const_t           ! Flags enabling/disabling matrix equations
        type(EqConstraints), pointer :: cons  ! Abbreviation for soln%constraints
        real(dp) :: lambda                    ! Damped update factor

        ! Get the solution update variables (pi, dnj_c, dln_n, dln_T, dln_nj)
        call self%get_solution_vars(soln)

        ! Define shorthand
        ng = self%num_gas
        nc = self%num_condensed
        na = count(soln%is_active)
        ne = self%num_active_elements()
        ne_full = self%num_elements
        num_eqn = soln%num_equations(self)
        cons => soln%constraints
        dln_nj => soln%dln_nj
        dnj_c => soln%dnj_c
        dln_n = soln%dln_n
        dln_T = soln%dln_T
        ln_nj => soln%ln_nj
        n = soln%n
        ln_n = log(n)
        const_p = cons%is_constant_pressure()
        const_t = cons%is_constant_temperature()
        T = soln%T

        ! Associate subarray pointers
        A_g => self%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        nj_g => soln%nj(:ng)
        nj_c => soln%nj(ng+1:)
        h_g => soln%thermo%enthalpy(:ng)
        s_g => soln%thermo%entropy(:ng)

        ! Compute the damped update factor
        lambda = self%compute_damped_update_factor(soln)

        ! Update gas species concentration
        do i = 1, ng
            ln_nj(i) = ln_nj(i) + lambda*dln_nj(i)
            if (ln_nj(i) - ln_n + self%tsize > 0.0d0) then
                nj_g(i) = exp(ln_nj(i))
            else
                nj_g(i) = 0.0d0
            end if
        end do

        ! Use a lower threshold for ionized species before truncating the concentrations
        if (self%ions .and. self%active_ions .and. ne_full > 0) then
            do i = 1, ng
                if (A_g(i, ne_full) /= 0.0d0 .and. nj_g(i) == 0.0d0) then
                    if (ln_nj(i) - ln_n + self%esize > 0.0d0) then
                        nj_g(i) = exp(ln_nj(i))
                    end if
                end if
            end do
        end if

        ! Condensed species concentrations
        active_idx = soln%active_condensed_indices()
        do idx_c = 1, size(active_idx)
            i = active_idx(idx_c)
            nj_c(i) = nj_c(i) + lambda*dnj_c(idx_c)
        end do

        ! Total moles
        if (const_p) then
            soln%n = exp(ln_n + lambda*dln_n)
        else
            soln%n = sum(nj_g)
        end if

        ! Temperature
        if (.not. const_t) then
            soln%T = exp(log(T) + lambda*dln_T)
        end if

        ! Update thermodynamic properties
        call self%products%calc_thermo(soln%thermo, soln%T, condensed=.true.)

    end subroutine

    subroutine EqSolver_check_convergence(self, soln)
        ! Check if the problem is converged

        ! Arguments
        class(EqSolver), target :: self
        type(EqSolution), intent(inout), target :: soln

        ! Locals
        integer  :: ng                        ! Number of gas species
        integer  :: na                        ! Number of active condensed species
        integer  :: ne                        ! Number of elements
        integer  :: ne_full                   ! Total number of elements (including electron)
        real(dp) :: b_delta(self%num_elements)! Residual for element contraints
        real(dp) :: s_delta                   ! Residual for entropy state
        real(dp), pointer :: b0(:)            ! Fixed element concentrations
        real(dp), pointer :: nj(:), nj_g(:)   ! Total/gas species concentrations [kmol-per-kg]
        real(dp) :: n, ln_n                   ! Total moles of mixture, log of total moles
        real(dp), pointer :: ln_nj(:)         ! Log of the product concentrations
        real(dp) :: P                         ! Mixture pressure (bar)
        real(dp), pointer :: h_g(:)           ! Gas enthalpies [unitless]
        real(dp), pointer :: s_g(:), s_c(:)   ! Gas entropies [unitless]
        real(dp), pointer :: A_g(:,:)         ! Gas stoichiometric matrices
        integer :: i, j                       ! Indices
        logical :: const_p, const_t, const_s  ! Flags enabling/disabling matrix equations
        type(EqConstraints), pointer :: cons  ! Abbreviation for soln%constraints
        real(dp), pointer :: pi(:)            ! 𝛑_j (k-th iteration)
        real(dp), pointer :: pi_prev(:)       ! 𝛑_j (k-1 iteration)
        real(dp), pointer :: dln_nj(:)        ! 𝛥ln(nj)
        real(dp), pointer :: dnj_c(:)         ! 𝛥nj_c
        real(dp) :: dln_n                     ! 𝛥ln(n)
        real(dp) :: dln_T                     ! 𝛥ln(T)
        real(dp), parameter :: nj_tol = 0.5d-5   ! Tolerance for species concentrations
        real(dp), parameter :: b_tol = 1.0d-6    ! Tolerance for the element concentrations
        real(dp), parameter :: T_tol = 1.0d-4    ! Tolerance for the temperature
        real(dp), parameter :: s_tol = 0.5d-4    ! Tolerance for the entropy
        real(dp), parameter :: pi_tol = 1.0d-3   ! Tolerance for modified lagrance multipliers
        real(dp), parameter :: ion_tol = 1.0d-4  ! Tolerance for ionized species
        real(dp) :: sum1, sum2, aa, temp         ! Temporary variables for ionized species

        ! Define shorthand
        ng = self%num_gas
        ne = self%num_active_elements()
        ne_full = self%num_elements
        na = count(soln%is_active)
        cons => soln%constraints
        nj => soln%nj
        ln_nj => soln%ln_nj
        n = soln%n
        ln_n = log(n)
        const_p = cons%is_constant_pressure()
        const_t = cons%is_constant_temperature()
        const_s = cons%is_constant_entropy()
        b0 => cons%b0
        pi => soln%pi
        pi_prev => soln%pi_prev
        dln_nj => soln%dln_nj
        dnj_c => soln%dnj_c
        dln_n = soln%dln_n
        dln_T = soln%dln_T

        ! Associate subarray pointers
        A_g => self%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        nj_g => soln%nj(:ng)
        h_g => soln%thermo%enthalpy(:ng)
        s_g => soln%thermo%entropy(:ng)
        s_c => soln%thermo%entropy(ng+1:)

        ! Get the mixture pressure
        P = soln%calc_pressure()

        ! Evalutate constraint residuals
        b_delta = b0 - self%products%elements_from_species(nj)
        if (const_s) s_delta = cons%state1 - dot_product(nj_g, s_g-ln_nj-log(P/n)) - dot_product(nj(ng+1:), s_c)

        ! Initialize the convergence to true
        soln%gas_converged         = .true.
        soln%condensed_converged   = .true.
        soln%moles_converged       = .true.
        soln%element_converged     = .true.
        soln%temperature_converged = .true.
        soln%entropy_converged     = .true.
        soln%pi_converged          = .true.
        soln%ions_converged        = .true.
        soln%converged             = .false.

        ! Check the converge of the species concentrations: Equation (3.5)
        ! ----------------------------------------------------------------

        ! Check gas species updates
        do i = 1, ng
            if ((nj_g(i)*abs(dln_nj(i))/sum(nj)) > nj_tol) then
                soln%gas_converged = .false.
                return
            end if
        end do

        ! Check condensed species updates
        do i = 1, na
            if (abs(dnj_c(i))/sum(nj) > nj_tol) then
                soln%condensed_converged = .false.
                return
            end if
        end do

        ! Check total moles update
        if (const_p) then
            if (abs(n*dln_n/sum(nj(1:ng))) > nj_tol) then
                soln%moles_converged = .false.
                return
            end if
        end if

        ! Check the convergence of the element concentrations: Equation (3.6a)
        ! --------------------------------------------------------------------
        do i = 1, ne
            if (b0(i) > b_tol) then
                if (abs(b_delta(i)) > b_tol*maxval(b0)) then
                    soln%element_converged = .false.
                    return
                end if
            end if
        end do

        ! Check the convergence of the temperature: Equation (3.6b)
        ! ---------------------------------------------------------
        if (.not. const_t) then
            if (abs(dln_T) > T_tol) then
                soln%temperature_converged = .false.
                return
            end if
        end if

        ! Check the convergence of entropy: Equation (3.6c)
        ! -------------------------------------------------
        if (const_s) then
            if (abs(s_delta) > s_tol) then
                soln%entropy_converged = .false.
                return
            end if
        end if

        ! Check the convergence of the modified lagrange multipliers: Equation (3.6d)
        ! ---------------------------------------------------------------------------
        if (self%trace > 0.0d0) then
            do i = 1, ne
                ! Match CEA2 behavior: if denominator is zero, skip ratio check.
                if (abs(pi(i)) > tiny(1.0d0)) then
                    if (abs((pi_prev(i) - pi(i))/pi(i)) > pi_tol) then
                        soln%pi_converged = .false.
                        return
                    end if
                end if
            end do
        end if

        ! Check total convergence
        soln%converged = .true.
        soln%times_converged = soln%times_converged + 1

        ! Update tsize after initial convergence, and adjust species concentrations
        self%tsize = self%xsize
        do i = 1, ng
            if (ln_nj(i) - ln_n + self%tsize > 0.0d0) then
                nj_g(i) = exp(ln_nj(i))
            else
                nj_g(i) = 0.0d0
            end if
        end do

        ! If everything converged, check ion convergence: Equation (3.14)
        ! ---------------------------------------------------------------
        soln%ions_converged = .false.
        if (self%ions .and. self%active_ions .and. ne_full > 0) then
            ! Check on electron balance
            do i = 1, 80  ! Max iterations
                sum1 = 0.0d0
                sum2 = 0.0d0
                do j = 1, ng
                    if (A_g(j, ne_full) /= 0.0d0) then
                        soln%nj(j) = 0.0d0
                        temp = 0.0d0
                        if (soln%ln_nj(j) > -87.0d0) temp = exp(soln%ln_nj(j))
                        if (soln%ln_nj(j) - log(n) + self%tsize > 0.0d0) then
                            !soln%pi_e = 0.0d0
                            soln%nj(j) = temp
                        end if
                        aa = A_g(j, ne_full)*temp
                        sum1 = sum1 + aa
                        sum2 = sum2 + aa*A_g(j, ne_full)
                    end if
                end do
                if (sum2 /= 0.0d0) then
                    soln%dpi_e = -sum1/sum2
                    do j = 1, ng
                        if (A_g(j, ne_full) /= 0.0d0) then
                            soln%ln_nj(j) = soln%ln_nj(j) + A_g(j, ne_full)*soln%dpi_e
                        end if
                    end do
                end if

                if (abs(soln%dpi_e) > ion_tol) then
                    soln%pi_e = soln%pi_e + soln%dpi_e
                    continue
                else
                    soln%ions_converged = .true.
                    return
                end if
            end do
            if (.not. soln%ions_converged) soln%converged = .false.
        end if
        soln%pi_e = 0.0d0

    end subroutine

    subroutine EqSolver_assemble_matrix(self, soln)
        class(EqSolver), intent(in), target :: self
        class(EqSolution), intent(inout), target :: soln

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: na                          ! Number of active condensed species
        integer  :: ne                          ! Number of elements
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: tmp(self%num_gas)           ! Common sub-expression storage
        real(dp) :: mu_g(self%num_gas)          ! Gas phase chemical potentials [unitless]
        real(dp) :: b_delta(self%num_elements)  ! Residual for element contraints
        real(dp) :: n_delta                     ! Residual for total moles / pressure constraint
        real(dp) :: hsu_delta                   ! Residual for enthalpy / entropy constraint
        real(dp) :: n                           ! Total moles of mixture
        real(dp) :: P                           ! Pressure of mixture (bar)
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp), pointer :: cp(:), cv(:)       ! Species heat capacities [unitless]
        real(dp), pointer :: h_g(:), h_c(:)     ! Gas/condensed enthalpies [unitless]
        real(dp), pointer :: s_g(:), s_c(:)     ! Gas/condensed entropies [unitless]
        real(dp), pointer :: u_g(:), u_c(:)     ! Gas/condensed energies [unitless]
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices
        real(dp), pointer :: G(:,:)             ! Augmented Newton iteration matrix
        real(dp), pointer :: h_or_s_or_u(:)     ! For evaluating Eq 2.27/2.28
        integer :: r, c                         ! Iteration matrix row/column indices
        integer :: i, j                         ! Loop counters
        integer, allocatable :: active_idx(:)   ! Active condensed indices in legacy order
        logical :: const_p, const_t, const_s, const_h, const_u  ! Flags enabling/disabling matrix equations
        type(EqConstraints), pointer :: cons    ! Abbreviation for soln%constraints

        ! Define shorthand
        ng = self%num_gas
        nc = self%num_condensed
        ne = self%num_active_elements()
        na = count(soln%is_active)
        active_idx = soln%active_condensed_indices()
        num_eqn = soln%num_equations(self)
        cons => soln%constraints
        const_p = cons%is_constant_pressure()
        const_t = cons%is_constant_temperature()
        const_s = cons%is_constant_entropy()
        const_h = cons%is_constant_enthalpy()
        const_u = cons%is_constant_energy()

        ! Associate subarray pointers
        G   => soln%G(:num_eqn, :num_eqn+1)
        A_g => self%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => self%products%stoich_matrix(ng+1:,:)
        n = soln%n
        nj  => soln%nj
        nj_g => soln%nj(:ng)
        ln_nj => soln%ln_nj
        cp  => soln%thermo%cp
        cv  => soln%thermo%cv
        h_g => soln%thermo%enthalpy(:ng)
        h_c => soln%thermo%enthalpy(ng+1:)
        s_g => soln%thermo%entropy(:ng)
        s_c => soln%thermo%entropy(ng+1:)
        u_g => soln%thermo%energy(:ng)
        u_c => soln%thermo%energy(ng+1:)

        ! Get the mixture pressure
        P = soln%calc_pressure()

        ! Compute gas phase chemical potentials
        mu_g = h_g - s_g + ln_nj + log(P/n)

        ! Evalutate constraint residuals
        b_delta = cons%b0 - self%products%elements_from_species(nj)
        n_delta = n - sum(nj_g)
        if (const_s) then
            hsu_delta = (cons%state1 - soln%calc_entropy_sum(self))
        else if (const_h) then
            hsu_delta = (cons%state1/soln%T - dot_product(nj, soln%thermo%enthalpy))
        else if (const_u) then
            hsu_delta = (cons%state1/soln%T - dot_product(nj, soln%thermo%energy))
        end if

        ! Initialize the iteration matrix
        G = 0.0d0
        r = 0
        c = 0

        !-------------------------------------------------------
        ! Equation (2.24/2.45): Element constraints
        !-------------------------------------------------------
        do i = 1,ne
            tmp = nj_g*A_g(:,i)
            r = r+1
            c = 0

            ! Pi derivatives
            do j = 1,ne
                c = c+1
                G(r,c) = dot_product(tmp, A_g(:,j))
            end do

            ! Condensed derivatives
            ! Handled in (2.25) below (symmetric)
            c = c+na

            ! Delta ln(n) derivative
            ! Symmetric with (2.26) pi derivative
            if (const_p) then
                c = c+1
                G(r,c) = sum(tmp)
                G(c,r) = G(r,c)
            end if

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
                if (const_p) then
                    G(r,c) = dot_product(tmp, h_g)
                else
                    G(r,c) = dot_product(tmp, u_g)
                end if
            end if

            ! Right hand side
            G(r,c+1) = b_delta(i) + dot_product(tmp, mu_g)

        end do

        !-------------------------------------------------------
        ! Equation (2.25/2.46): Condensed phase constraints
        !-------------------------------------------------------
        do j = 1, na
            i = active_idx(j)
            r = r+1
            c = 0

            ! Pi derivatives
            ! Symmetric with (2.24) condensed derivatives
            G(r,:ne) = A_c(i,:ne)
            G(:ne,r) = A_c(i,:ne)
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            ! Handled via zero initialization
            if (const_p) c = c+1

            ! Delta ln(T) derivative
            if (.not. const_t) then
                c = c+1
                if (const_p) then
                    G(r,c) = h_c(i)
                else
                    G(r,c) = u_c(i)
                end if
            end if

            ! Right hand size
            G(r,c+1) = h_c(i) - s_c(i) ! = mu_c(i)

        end do

        !-------------------------------------------------------
        ! Equation (2.26)
        !-------------------------------------------------------
        if (const_p) then
            r = r+1
            c = 0

            ! Pi derivatives
            ! Handled in Eq (2.24) above.
            c = c+ne

            ! Condensed derivatives
            ! Handled via zero initialization of G
            c = c+na

            ! Delta ln(n) derivative
            c = c+1
            G(r,c) = -n_delta

            ! Delta ln(T) derviative
            if (.not. const_t) then
                c = c+1
                G(r,c) = dot_product(nj_g, h_g)
            end if

            ! Right-hand-side
            G(r,c+1) = n_delta + dot_product(nj_g, mu_g)

        end if

        !---------------------------------------------------------
        ! Equation (2.27)/(2.28)/(2.47)/(2.48): Energy constraints
        !---------------------------------------------------------
        if (.not. const_t) then
            r = r+1
            c = 0

            ! Select entropy/enthalpy constraint
            if (const_s) then
                tmp = nj_g*(h_g-mu_g)!(s_g-ln_nj-log(P/n))
                h_or_s_or_u => soln%thermo%entropy(ng+1:)
            else if (const_h) then
                tmp = nj_g*h_g
                h_or_s_or_u => soln%thermo%enthalpy(ng+1:)
            else if (const_u) then
                tmp = nj_g*u_g
                h_or_s_or_u => soln%thermo%energy(ng+1:)
            end if

            ! Pi derivatives
            do j = 1,ne
                c = c+1
                G(r,c) = dot_product(tmp, A_g(:,j))
                if (.not. const_p .and. const_s) then
                    G(r,c) = G(r,c) - dot_product(nj_g, A_g(:, j))
                end if
            end do

            ! Condensed derivatives
            do j = 1, na
                i = active_idx(j)
                c = c+1
                G(r,c) = h_or_s_or_u(i)
            end do

            ! Delta ln(n) derivative
            if (const_p) then
                c = c+1
                G(r,c) = sum(tmp)
            end if

            ! Delta ln(T) derivative
            c = c+1
            if (const_p) then
                G(r,c) = dot_product(nj, cp) + dot_product(tmp, h_g)
            else
                G(r,c) = dot_product(nj, cv) + dot_product(tmp, u_g)
                if (const_s) then
                    G(r,c) = G(r,c) - dot_product(nj_g, u_g)
                end if
            end if

            ! Right-hand-side
            G(r,c+1) = hsu_delta + dot_product(tmp, mu_g)
            if (const_s) then
                if (const_p) then
                    G(r,c+1) = G(r,c+1) + n_delta
                else
                    G(r,c+1) = G(r,c+1) - dot_product(nj_g, mu_g)
                end if
            end if

        end if

    end subroutine

    subroutine EqSolver_check_condensed_phases(self, soln, iter, made_change)
        ! Check if condensed species phases makes sense for the current temperature

        ! Arguments
        class(EqSolver), intent(in), target :: self
        type(EqSolution), intent(inout), target :: soln
        integer, intent(inout) :: iter  ! Global solver iteration counter; reset to -1 if condensed species are added or removed
        logical, intent(out) :: made_change  ! Flag to indicate if a species was added or removed

        ! Locals
        integer :: ng                              ! Number of gas species
        integer :: nc                              ! Number of condensed species
        integer :: na                              ! Number of active condensed species
        integer :: i, j, idx_c                     ! Index
        integer, allocatable :: active_idx(:)      ! Active condensed indices in legacy order
        real(dp) :: T_low_i, T_high_i              ! Low and high temperature limits for a species [K]
        real(dp) :: T_low_j, T_high_j              ! Low and high temperature limits for a species [K]
        real(dp) :: max_T_j                        ! Max melting temperature of the candidate phase [K]
        integer, allocatable :: idx_other_phase(:) ! Indices of other phases of the same species
        real(dp), parameter :: dT_phase = 50.0d0   ! Temperature gap (above/below transition) to force phase change [K]
        real(dp), parameter :: dT_limit = 1.2d0    ! Factor of safety on high temperature range for a species to exist [unitless]
        real(dp), parameter :: T_tol = 1.d-3       ! Tolerance for temperature range comparison [K]
        integer :: d_phase                         ! Integer difference in phase between the two species

        ! Shorthand
        ng = self%num_gas
        nc = self%num_condensed
        na = count(soln%is_active)

        made_change = .false.

        ! Legacy CEA applies condensed-phase validity checks during TP solves too.
        if (na == 0) return

        ! Update condensed thermodynamic properties
        ! call self%products%calc_thermo(soln%thermo, soln%T, condensed=.true.)

        active_idx = soln%active_condensed_indices()
        do idx_c = 1, na
            i = active_idx(idx_c)

            if (i == soln%j_sol .or. i == soln%j_liq) cycle

            T_low_i  = minval(self%products%species(ng+i)%T_fit(:, 1))
            T_high_i = maxval(self%products%species(ng+i)%T_fit(:, 2))

            ! Get the indices of the other phases of the same species
            idx_other_phase = get_species_other_phases(self%products%species_names(ng+i), &
                                                       self%products)

            ! Reset
            max_T_j = 0.0d0

            do j = 1, size(idx_other_phase)

                d_phase = self%products%species(ng+i)%i_phase - self%products%species(ng+idx_other_phase(j))%i_phase

                ! Get the temperature range of the cantidate phase
                T_low_j  = minval(self%products%species(ng+idx_other_phase(j))%T_fit(:, 1))
                T_high_j = maxval(self%products%species(ng+idx_other_phase(j))%T_fit(:, 2))

                if (T_high_j > max_T_j) max_T_j = T_high_j

                ! Check that the temperature is less than the upper bound of the cantidate phase
                if (soln%T <= T_high_j) then
                    if (d_phase /= 0) then

                        ! Check if the cantidate phase differs by more than 1
                        if (abs(d_phase) > 1) then
                            ! Switch phase
                            call log_info("Phase change: replace "//trim(self%products%species_names(ng+i))//&
                                          " with "//self%products%species_names(ng+idx_other_phase(j)))
                            call soln%replace_active_condensed(i, idx_other_phase(j))
                            soln%nj(ng+idx_other_phase(j)) = soln%nj(ng+i)
                            soln%nj(ng+i) = 0.0d0
                            soln%converged = .false.
                            soln%j_switch = i
                            soln%j_sol = 0
                            soln%j_liq = 0
                            iter = -1
                            made_change = .true.
                            return
                        end if

                        ! If the cantidate species is the one we just removed, keep both phases
                        if (idx_other_phase(j) == soln%j_switch) then
                            call log_debug("Cantidate phase ("//trim(self%products%species_names(ng+idx_other_phase(j)))//&
                                           ") is the one that was removed last ("&
                                           //trim(self%products%species_names(ng+i))//"). Keep both.")
                            soln%T = min(T_high_i, T_high_j)  ! Set T as the melting temperature
                            if (T_high_i > T_high_j) then
                                soln%j_sol = idx_other_phase(j)
                                soln%j_liq = i
                            else
                                soln%j_sol = i
                                soln%j_liq = idx_other_phase(j)
                            end if
                            call soln%activate_condensed_front(idx_other_phase(j))
                            soln%nj(ng+idx_other_phase(j)) = 0.5d0*soln%nj(ng+i)
                            soln%nj(ng+i)                  = 0.5d0*soln%nj(ng+i)
                            soln%converged = .false.
                            iter = -1
                            made_change = .true.
                            return
                        end if

                        ! The solution temperature is outside the allowable range for the existing phase
                        if (soln%T < (T_low_i-dT_phase) .or. soln%T > (T_high_i+dT_phase)) then
                            ! Switch phase
                            call log_info("Phase change: replace "//trim(self%products%species_names(ng+i))//&
                                          " with "//self%products%species_names(ng+idx_other_phase(j)))
                            call soln%replace_active_condensed(i, idx_other_phase(j))
                            soln%nj(ng+idx_other_phase(j)) = soln%nj(ng+i)
                            soln%nj(ng+i) = 0.0d0
                            soln%converged = .false.
                            soln%j_switch = i
                            soln%j_sol = 0
                            soln%j_liq = 0
                            iter = -1
                            made_change = .true.
                            return
                        end if

                        ! else
                        call log_debug("Adding "//self%products%species_names(ng+idx_other_phase(j)))
                        soln%T = min(T_high_i, T_high_j)  ! Set T as the melting temperature
                        if (T_high_i > T_high_j) then
                            soln%j_sol = idx_other_phase(j)
                            soln%j_liq = i
                        else
                            soln%j_sol = i
                            soln%j_liq = idx_other_phase(j)
                        end if
                        call soln%activate_condensed_front(idx_other_phase(j))
                        soln%nj(ng+idx_other_phase(j)) = 0.5d0*soln%nj(ng+i)
                        soln%nj(ng+i)                  = 0.5d0*soln%nj(ng+i)
                        soln%converged = .false.
                        iter = -1
                        made_change = .true.
                        return
                    else
                        ! Exit the inner loop
                        exit
                    end if
                end if
            end do

            if (soln%T > 1.2d0*max_T_j) then
                ! Remove condensed species
                call log_info("Removing condensed species: "//self%products%species_names(ng+i))
                call soln%deactivate_condensed(i)
                soln%nj(ng+i) = 0.0d0
                soln%converged = .false.
                soln%j_switch = 0
                soln%j_sol = 0
                soln%j_liq = 0
                iter = -1
                made_change = .true.
                return
            end if

        end do

    end subroutine

    subroutine EqSolver_test_condensed(self, soln, iter, singular_index)
        ! Test if adding condensed species improves the solution

        ! Arguments
        class(EqSolver), intent(in), target :: self
        type(EqSolution), intent(inout), target :: soln
        integer, intent(inout) :: iter  ! Global solver iteration counter; reset to -1 if condensed species are added or removed
        integer, intent(in), optional :: singular_index

        ! Locals
        integer :: ng                             ! Number of gas species
        integer :: ne                             ! Number of active elements
        integer :: na                             ! Number of active condensed species
        integer :: nc                             ! Number of total condensed species
        real(dp), pointer :: nj_c(:)              ! Condensed species concentrations [kmol-per-kg]
        real(dp), pointer :: cp_c(:)              ! Condensed pecies heat capacities [unitless]
        integer :: i, j                           ! Indices
        integer :: cond_idx                       ! Index of condensed species to add/remove
        integer :: singular_index_                ! Index of condensed species that caused a singular matrix
        real(dp) :: temp                          ! Temp value to select condensed species
        real(dp) :: delg                          ! ∂g (Gibb's energy) of a cantidate condensed species [unitless]
        real(dp) :: min_delg                      ! Largest negative ∂g (Gibb's energy) of a cantidate condensed species [unitless]
        real(dp) :: delg_singular                 ! ∂g of the condensed species with singular_index
        real(dp), pointer :: h_c(:)               ! Condensed enthalpies [unitless]
        real(dp), pointer :: s_c(:)               ! Condensed entropies [unitless]
        real(dp), pointer :: A_c(:,:)             ! Condensed stoichiometric matrices
        real(dp), pointer :: pi(:)                ! 𝛑_j (k-th iteration)
        integer, allocatable :: active_idx(:)     ! Active condensed indices in legacy order
        logical :: made_change                    ! Flag to indicate if a species was added or removed (used for other subroutine calls)
        real(dp), parameter :: T_min = 200.0d0    ! Minimum gas temperature defined in thermo data [K]
        real(dp), parameter :: tol = 1d-12

        ! Shorthand
        ng = self%num_gas
        ne = self%num_active_elements()
        na = count(soln%is_active)
        nc = self%num_condensed

        ! Associate subarray pointers
        nj_c => soln%nj(ng+1:)
        cp_c => soln%thermo%cp(ng+1:)
        A_c => self%products%stoich_matrix(ng+1:,:)
        h_c => soln%thermo%enthalpy(ng+1:)
        s_c => soln%thermo%entropy(ng+1:)
        pi => soln%pi_prev

        ! Set default value for optional argument
        singular_index_ = 0
        if (present(singular_index)) singular_index_ = singular_index

        ! If no condensed species to consider, return
        if (nc == 0) return

        ! Remove the condensed species with the largest negative value of nj*Cp
        if (na > 0) then
            temp = 0.0d0
            cond_idx = 0
            active_idx = soln%active_condensed_indices()

            do i = 1, size(active_idx)
                j = active_idx(i)
                if (nj_c(j)*cp_c(j) <= temp) then
                    temp = nj_c(j)*cp_c(j)
                    cond_idx = j
                end if
            end do

            if (cond_idx > 0) then
                if (cond_idx==soln%j_sol .or. cond_idx==soln%j_liq) then
                    soln%j_sol = 0
                    soln%j_liq = 0
                end if
                call soln%deactivate_condensed(cond_idx)
                soln%nj(ng+cond_idx) = 0.0d0
                soln%converged = .false.
                iter = -1
                call log_info("Removing "//self%products%species_names(ng+cond_idx))
                return
            end if

        end if

        ! Check if any existing condensed species need to change phase
        call self%check_condensed_phases(soln, iter, made_change)
        if (made_change) return

        ! Check if adding any condensed species improves the solution
        min_delg = 0.0d0
        delg_singular = 0.0d0
        cond_idx = 0
        do i = 1, nc

            if (soln%is_active(i)) cycle

            ! Check if this species can be present at the current temperature
            if (soln%T >= minval(self%products%species(ng+i)%T_fit(:, 1)) .or. &
                T_min == minval(self%products%species(ng+i)%T_fit(:, 1))) then
                if (soln%T <= maxval(self%products%species(ng+i)%T_fit(:, 2))) then

                    temp = 0.0d0
                    if (ne > 0) temp = dot_product(A_c(i,:ne), pi(:ne))
                    delg = (h_c(i) - s_c(i) - temp)/self%products%species(ng+i)%molecular_weight

                    if (delg < min_delg .and. delg < 0.0d0) then
                        if (i /= singular_index_) then
                            min_delg = delg
                            cond_idx = i
                        else
                            delg_singular = delg
                        end if
                    end if

                end if
            end if

        end do

        if (abs(min_delg) < tol .and. abs(delg_singular) < tol) return  ! Converged with no condensed changes

        ! Insert the selected condensed species
        if (abs(min_delg) > tol) then
            call soln%activate_condensed_front(cond_idx)
            soln%converged = .false.
            iter = -1
            soln%last_cond_idx = cond_idx
            call log_info("Adding "//self%products%species_names(ng+cond_idx))
            return
        else
            call abort('EqSolver_test_condensed: Re-insertion of '// &
                self%products%species_names(ng+singular_index_)//' likely to cause singular matrix.')
        end if

    end subroutine

    subroutine EqSolver_correct_singular(self, soln, iter, ierr, singular_index, reduced_from, reduced_to)
        ! Try to correct the singular Jacobian matrix

        ! Arguments
        class(EqSolver), target :: self
        type(EqSolution), intent(inout), target :: soln
        integer, intent(inout) :: iter
        integer, intent(in) :: ierr
        integer, intent(out), optional :: singular_index
        integer, intent(out), optional :: reduced_from
        integer, intent(out), optional :: reduced_to

        ! Locals
        integer :: i, j, k                   ! Iterators
        real(dp) :: temp                     ! Temporary summation value
        integer :: idx                       ! Condensed species index
        integer :: idx_active                ! Active condensed row index in the Jacobian
        integer :: ng                        ! Number of gas species
        integer :: nc                        ! Number of condensed species
        integer :: ne                        ! Number of elements
        integer :: na                        ! Number of active condensed species
        integer, allocatable :: active_idx(:)! Active condensed indices in legacy order
        real(dp), pointer :: A(:,:)          ! Stoichiometric matrix
        real(dp), pointer :: A_all(:,:)      ! Full stoichiometric matrix
        real(dp), parameter :: tol = 1.d-8   ! Tolerance to check if value ~0
        real(dp), parameter :: smalno = 1.0d-6
        real(dp), parameter :: smnol = -13.815511d0
        logical :: made_change


        ! Shorthand
        ng = self%num_gas
        nc = self%num_condensed
        ne = self%num_active_elements()
        na = count(soln%is_active)
        active_idx = soln%active_condensed_indices()
        A => self%products%stoich_matrix(ng+1:,:)
        A_all => self%products%stoich_matrix(:,:)

        self%xsize = 80.0d0
        self%tsize = 80.0d0
        if (present(singular_index)) singular_index = 0
        if (present(reduced_from)) reduced_from = 0
        if (present(reduced_to)) reduced_to = 0
        made_change = .false.

        if (ierr > ne .and. iter < 1 .and. na > 1 &
            .and. soln%last_cond_idx > 0) then

            temp = 1000.0d0
            idx = 0  ! Condensed species index selected to correct singular matrix
            do k = 1, size(active_idx)
                i = active_idx(k)

                if (i /= soln%last_cond_idx) then
                    do j = 1, ne
                        if (abs(A(soln%last_cond_idx, j)) > tol .and. abs(A(i, j)) > tol) then
                            if (soln%nj(ng+i) <= temp) then
                                temp = soln%nj(ng+i)
                                idx = i
                            end if
                        end if
                    end do
                end if
            end do
            ! Remove condensed species contributing to singular matrix
            if (idx > 0) then
                call log_info("Removing condensed species "//self%products%species_names(ng+idx)// &
                              " to correct singular matrix")
                call soln%deactivate_condensed(idx)
                soln%nj(ng+idx) = 0.0d0
                soln%converged = .false.
                soln%j_switch = idx
                if (present(singular_index)) singular_index = idx
                iter = -1
                made_change = .true.
            end if

        ! Legacy-inspired recovery for element-row singularities:
        ! remove the smallest active condensed species that participates
        ! in the singular element equation.
        else if (ierr >= 1 .and. ierr <= ne .and. na > 0) then
            temp = huge(1.0d0)
            idx = 0
            do k = 1, size(active_idx)
                i = active_idx(k)
                if (abs(A(i, ierr)) <= tol) cycle
                if (soln%nj(ng+i) <= temp) then
                    temp = soln%nj(ng+i)
                    idx = i
                end if
            end do

            if (idx > 0) then
                call log_info("Removing condensed species "//self%products%species_names(ng+idx)// &
                              " to correct element-row singularity")
                call soln%deactivate_condensed(idx)
                soln%nj(ng+idx) = 0.0d0
                soln%converged = .false.
                soln%j_switch = idx
                if (present(singular_index)) singular_index = idx
                iter = -1
                made_change = .true.
            end if

        ! Remove condensed species to correct singularity
        else if (ierr > ne .and. ierr <= na+ne) then
            ! Map Jacobian active condensed row index to condensed species index.
            idx_active = ierr - ne
            idx = 0
            if (idx_active >= 1 .and. idx_active <= size(active_idx)) idx = active_idx(idx_active)

            if (idx > 0 .and. soln%is_active(idx)) then
                call soln%deactivate_condensed(idx)
                soln%nj(self%num_gas+idx) = 0.0d0
                soln%converged = .false.
                soln%j_switch = idx
                if (present(singular_index)) singular_index = idx
                iter = -1
                made_change = .true.
            end if
        end if

        ! Legacy ion-row fallback: if the electron equation is singular, remove
        ! ionized species from the active iterate and disable ion solving.
        if (.not. made_change .and. ierr >= 1 .and. ierr <= ne .and. &
            self%ions .and. self%active_ions .and. self%num_elements > 0) then
            if (trim(self%products%element_names(self%num_elements)) == 'E' .and. ierr == ne) then
                do i = 1, ng
                    if (A_all(i, self%num_elements) /= 0.0d0) then
                        soln%nj(i) = 0.0d0
                        soln%ln_nj(i) = smnol
                        made_change = .true.
                    end if
                end do
                if (made_change) then
                    self%active_ions = .false.
                    soln%pi_e = 0.0d0
                    soln%dpi_e = 0.0d0
                    soln%converged = .false.
                    iter = -1
                end if
            end if
        end if

        ! Legacy component-reduction fallback for persistent element-row singularities.
        if (.not. made_change .and. ierr >= 1 .and. ierr <= ne .and. iter < 1 .and. &
            ne > 1 .and. .not. (self%ions .and. self%active_ions)) then
            call log_info("Reducing active element equations after singular restart on "// &
                          trim(self%products%element_names(ierr)))
            if (ierr /= ne) call EqSolver_swap_elements(self, soln, ierr, ne)
            self%reduced_elements = self%reduced_elements + 1
            soln%pi(ne) = 0.0d0
            soln%pi_prev(ne) = 0.0d0
            soln%converged = .false.
            iter = -1
            made_change = .true.
            if (present(reduced_from)) reduced_from = ierr
            if (present(reduced_to)) reduced_to = ne
        end if

        ! Legacy fallback path: seed trace gas species to break persistent singularity.
        if (.not. made_change) then
            do i = 1, ng
                if (soln%nj(i) <= 0.0d0) then
                    soln%nj(i) = smalno
                    soln%ln_nj(i) = smnol
                    made_change = .true.
                end if
            end do
            if (made_change) then
                soln%converged = .false.
                iter = -1
            end if
        end if

    end subroutine

    subroutine EqSolver_post_process(self, soln, computed_partials)
        ! Arguments
        class(EqSolver), intent(in), target :: self
        type(EqSolution), intent(inout), target :: soln
        logical, intent(in), optional :: computed_partials

        ! Locals
        integer :: i
        logical :: computed_partials_

        computed_partials_ = .false.
        if (present(computed_partials)) computed_partials_ = computed_partials

        ! Compute the mole and mass fractions
        soln%mole_fractions = soln%nj / sum(soln%nj)
        soln%mass_fractions = soln%nj * self%products%species%molecular_weight / &
            sum(soln%nj * self%products%species%molecular_weight)

        ! Add mixture properties
        soln%pressure = soln%calc_pressure()
        soln%volume   = soln%calc_volume()
        soln%density  = 1.0d0/soln%volume

        soln%enthalpy = dot_product(soln%nj, soln%thermo%enthalpy) * R * soln%T / 1.d3
        soln%energy   = soln%enthalpy - soln%n*soln%T*R/1.d3
        soln%entropy  = soln%calc_entropy_sum(self) * R / 1.d3
        soln%gibbs_energy = (soln%enthalpy - soln%T*soln%entropy)

        if (soln%cp_fr < 1.d-10) soln%cp_fr = dot_product(soln%thermo%cp, soln%nj) * R / 1.d3
        if (.not. computed_partials_ .and. soln%cp_eq < 1.d-10) then
            soln%cp_eq = dot_product(soln%thermo%cp, soln%nj) * R / 1.d3
        end if

        ! Calculate molecular weights
        soln%M = 1.0d0/soln%n
        soln%MW = 1.0d0
        do i = 1, self%num_condensed
            if (soln%is_active(i)) then
                soln%MW = soln%MW - soln%mole_fractions(self%num_gas+i)
            end if
        end do
        soln%MW = soln%M*soln%MW

        ! Calculate Cv
        soln%cv_fr = soln%cp_fr - soln%n*R/1.d3

        ! Gamma_s
        if (soln%gamma_s < 1.d-10) soln%gamma_s = soln%cp_eq/(soln%cp_eq - soln%n)

    end subroutine

    subroutine EqSolver_solve(self, soln, type, state1, state2, reactant_weights, partials)

        ! Arguments
        class(EqSolver), target :: self
        type(EqSolution), intent(inout), target :: soln
        character(2), intent(in) :: type
        real(dp), intent(in) :: state1
        real(dp), intent(in) :: state2
        real(dp), intent(in) :: reactant_weights(:)
        type(EqPartials), intent(out), optional :: partials

        ! Locals
        integer :: i, iter, ierr, num_eqn, times_singular
        integer :: cond_idx
        integer :: singular_index, singular_index_iter
        integer :: reduced_from_iter, reduced_to_iter
        integer :: num_reduced
        integer :: reduced_from(self%num_elements), reduced_to(self%num_elements)
        integer :: phase_iter, phase_pass
        real(dp) :: gas_moles, xi, xln
        real(dp), pointer :: G(:, :)
        type(EqPartials) :: partials_
        logical :: made_change, max_iter_fallback_used

        call log_debug("Starting Eq. Solve.")

        ! Set problem type, fixed-state values, and element amounts
        call soln%constraints%set( &
            type, state1, state2, &
            self%reactants%element_amounts_from_weights(reactant_weights) &
        )

        ! If fixed-temperature, set it
        if (soln%constraints%is_constant_temperature()) then
            soln%T = state1
        end if

        ! Initialize values
        self%tsize = 18.420681d0  ! Re-set in case solver is being re-used
        times_singular = 0  ! Number of times a singular matrix was encountered ("ixsing" in CEA2)
        soln%times_converged = 0  ! Number of times initial convergence was established
        soln%j_switch = 0  ! Make sure this is reset every time
        self%active_ions = self%ions
        self%reduced_elements = 0
        soln%pi_e = 0.0d0
        num_reduced = 0
        reduced_from = 0
        reduced_to = 0

        ! Pre-check active condensed phases before the first Newton matrix build.
        phase_iter = 0
        do phase_pass = 1, self%num_condensed + 1
            call self%check_condensed_phases(soln, phase_iter, made_change)
            if (.not. made_change) exit
        end do

        ! Initial call of the thermodynamic properties.
        ! Compute condensed thermo too so the first Newton build does not use
        ! stale condensed values when active condensed species are present.
        call self%products%calc_thermo(soln%thermo, soln%T, condensed=.true.)

        ierr = 0
        iter = 0
        singular_index = 0
        max_iter_fallback_used = .false.
        do while (self%max_iterations > iter)

            iter = iter + 1

            ! Assemble the matrix
            call self%assemble_matrix(soln)

            ! Get the size of the matrix that we need
            num_eqn = soln%num_equations(self)
            G => soln%G(:num_eqn, :num_eqn+1)

            call gauss(G, ierr)

            if (ierr == 0) then
                call self%update_solution(soln)
                call self%check_convergence(soln)

                ! Update pi_prev
                soln%pi_prev = soln%pi

                if (.not. soln%converged) cycle

            else
                call log_warning('Singular update matrix encountered at iteration '//to_str(iter))

                times_singular = times_singular + 1
                if (times_singular > 8) then
                    call EqSolver_restore_reduced_elements(self, soln, num_reduced, reduced_from, reduced_to)
                    call abort('EqSolver_solve: Too many singular matrices encountered.')
                end if

                ! Try to correct the singular matrix
                singular_index_iter = 0
                reduced_from_iter = 0
                reduced_to_iter = 0
                call self%correct_singular(soln, iter, ierr, singular_index_iter, reduced_from_iter, reduced_to_iter)
                if (singular_index_iter > 0) singular_index = singular_index_iter
                if (reduced_to_iter > 0) then
                    num_reduced = num_reduced + 1
                    reduced_from(num_reduced) = reduced_from_iter
                    reduced_to(num_reduced) = reduced_to_iter
                end if

                ! Start next iteration
                cycle

            end if

            if (soln%times_converged > 3*self%num_active_elements()) then
                soln%converged = .false.
                call EqSolver_restore_reduced_elements(self, soln, num_reduced, reduced_from, reduced_to)
                call abort("Convergence failed to establish set of condensed species.")
            end if

            ! Initial convergence; check on adding or removing condensed species
            call self%test_condensed(soln, iter, singular_index)

            if (soln%converged .or. (iter == self%max_iterations)) then

                ! Compute final species concentrations
                ! * NOTE: post-processing uses a lower threshold when computing nj = exp(ln(nj))
                do i = 1,self%num_gas
                    if (soln%ln_nj(i) > self%log_min) soln%nj(i) = exp(soln%ln_nj(i))
                end do

                if (.not. soln%converged) then
                    ! Legacy-style fallback for high-temperature condensed edge cases.
                    gas_moles = sum(soln%nj(:self%num_gas))
                    if (.not. max_iter_fallback_used) then
                        if (self%num_gas > 0 .and. &
                            (.not. soln%constraints%is_constant_enthalpy() .or. soln%T > 100.0d0) .and. &
                            (count(soln%is_active) == 1) .and. (gas_moles <= 1.0d-4)) then
                            max_iter_fallback_used = .true.
                            soln%n = 0.1d0
                            xi = soln%n / self%num_gas
                            xln = log(xi)
                            do i = 1, self%num_gas
                                soln%nj(i) = xi
                                soln%ln_nj(i) = xln
                            end do
                            do cond_idx = 1, self%num_condensed
                                if (soln%is_active(cond_idx)) then
                                    call soln%deactivate_condensed(cond_idx)
                                    soln%nj(self%num_gas+cond_idx) = 0.0d0
                                    soln%j_switch = cond_idx
                                    exit
                                end if
                            end do
                            soln%j_sol = 0
                            soln%j_liq = 0
                            soln%converged = .false.
                            soln%times_converged = 0
                            iter = -1
                            call self%products%calc_thermo(soln%thermo, soln%T, condensed=.false.)
                            cycle
                        end if
                    end if

                    call self%post_process(soln, .false.)
                    call EqSolver_restore_reduced_elements(self, soln, num_reduced, reduced_from, reduced_to)
                    call abort('EqSolver_solve: Maximum iterations reached without convergence')
                end if

                ! Compute the partial derivatives
                if (present(partials) .or. self%transport) then

                    if (present(partials)) then
                        partials = EqPartials(self%num_elements, count(soln%is_active))
                        call partials%compute_partials(self, soln)
                    else  ! Transport; partials are required
                        partials_ = EqPartials(self%num_elements, count(soln%is_active))
                        call partials_%compute_partials(self, soln)
                    end if

                end if

                ! Compute transport properties
                if (self%transport) call compute_transport_properties(self, soln)

                ! Compute post-processing solution values
                call self%post_process(soln, present(partials))

                ! Check for temperature outside of bounds
                if (soln%T > self%T_max .or. soln%T < self%T_min) then
                    call log_warning("Mixture temperature outside of allowable bounds.")
                    soln%converged = .false.
                end if

                call EqSolver_restore_reduced_elements(self, soln, num_reduced, reduced_from, reduced_to)
                return
            end if

        end do

        call EqSolver_restore_reduced_elements(self, soln, num_reduced, reduced_from, reduced_to)

    end subroutine


    !-----------------------------------------------------------------------
    !  EqConstraint Implementation
    !-----------------------------------------------------------------------
    function EqConstraints_init(type, state1, state2, element_moles) result(self)
        ! Construct a new EqConstraints object
        type(EqConstraints) :: self
        character(2), intent(in) :: type
        real(dp), intent(in) :: state1
        real(dp), intent(in) :: state2
        real(dp), intent(in) :: element_moles(:)
        call self%set(type, state1, state2, element_moles)
    end function

    function EqConstraints_alloc(num_elements) result(self)
        ! Pre-allocate and empty-initialize a EqConstraints object.
        ! This makes "set" method allocation-free, provided num_elements = size(element_moles)
        type(EqConstraints) :: self
        integer, intent(in) :: num_elements
        self%type   = '  '
        self%state1 = empty_dp
        self%state2 = empty_dp
        self%b0     = empty_dp * ones(num_elements)
    end function

    subroutine EqConstraints_set(self, type, state1, state2, element_moles)
        class(EqConstraints), intent(inout) :: self
        character(2), intent(in) :: type
        real(dp), intent(in) :: state1
        real(dp), intent(in) :: state2
        real(dp), intent(in) :: element_moles(:)
        character(2) :: ltype

        ltype = lower(type)

        select case (ltype)

            case ('tp','hp','sp','tv','uv','sv')
                self%type   = ltype
                self%state1 = state1
                self%state2 = state2

            case ('pt','ph','ps','vt','vu','vs')
                ! Flip state order to match "default"
                self%type(1:1) = ltype(2:2)
                self%type(2:2) = ltype(1:1)
                self%state1 = state2
                self%state2 = state1

            case default
                call abort('Invalid equilibrium problem type: '//type)

        end select

        self%b0 = element_moles

    end subroutine

    function EqConstraints_is_constant_temperature(self) result(tf)
        class(EqConstraints), intent(in) :: self
        logical :: tf
        tf = (self%type(1:1) == 't')
    end function

    function EqConstraints_is_constant_enthalpy(self) result(tf)
        class(EqConstraints), intent(in) :: self
        logical :: tf
        tf = (self%type(1:1) == 'h')
    end function

    function EqConstraints_is_constant_energy(self) result(tf)
        class(EqConstraints), intent(in) :: self
        logical :: tf
        tf = (self%type(1:1) == 'u')
    end function

    function EqConstraints_is_constant_entropy(self) result(tf)
        class(EqConstraints), intent(in) :: self
        logical :: tf
        tf = (self%type(1:1) == 's')
    end function

    function EqConstraints_is_constant_pressure(self) result(tf)
        class(EqConstraints), intent(in) :: self
        logical :: tf
        tf = (self%type(2:2) == 'p')
    end function

    function EqConstraints_is_constant_volume(self) result(tf)
        class(EqConstraints), intent(in) :: self
        logical :: tf
        tf = (self%type(2:2) == 'v')
    end function


    !-----------------------------------------------------------------------
    ! EquilibriumSolution
    !-----------------------------------------------------------------------
    function EqSolution_init(solver, T_init, nj_init) result(self)
        type(EqSolution) :: self
        type(EqSolver), intent(in) :: solver
        real(dp), intent(in), optional :: T_init
        real(dp), intent(in), optional :: nj_init(:)

        ! Locals
        integer :: i, j

        ! Allocate data structures
        allocate(self%nj(solver%num_products), source=0.0d0)
        allocate(self%ln_nj(solver%num_gas), source=0.0d0)
        allocate(self%G(solver%max_equations, solver%max_equations+1), source=empty_dp)
        allocate(self%is_active(solver%num_condensed), source=.false.)
        allocate(self%active_rank(solver%num_condensed), source=0)
        self%constraints = EqConstraints(solver%num_elements)

        ! Set initial guess
        ! From CEA2: Assume a temperature of 3800K with a total molar
        ! concentration of 0.1d0. The total mole count is split evenly between
        ! all gas species + and condensed species requested via the "insert"
        ! keyword. For simplicity, we ignore the inserted species here; that
        ! modification can be introduced at to solver level.

        ! Set the initial temperature
        if (present(T_init)) then
            self%T = T_init
        else
            self%T = 3800.0d0
        end if

        ! Set the initial mole fractions
        if (present(nj_init)) then
            call self%set_nj(solver, nj_init)
        else
            self%n = 0.1d0
            self%nj(:solver%num_gas) = self%n / solver%num_gas
            self%ln_nj = log(self%nj(:solver%num_gas))
        end if

        ! Set the inserted species as active
        if (allocated(solver%insert)) then
            do i = 1, size(solver%insert)
                j = findloc(solver%products%species_names, solver%insert(i), 1)
                if (j > 0) then
                    ! Only count this as an "insert" if it is condensed; no effect otherwise
                    if (solver%products%species(j)%i_phase > 0) then
                        call log_info("Inserting "//solver%products%species_names(j))
                        call self%activate_condensed_front(j-solver%num_gas)
                    end if
                end if
            end do
            call solver%products%calc_thermo(self%thermo, self%T, condensed=.true.)
        else
            call solver%products%calc_thermo(self%thermo, self%T, condensed=.false.)
        end if

        ! Allocate solution update variables
        allocate(self%pi(solver%num_elements), source=0.0d0)
        allocate(self%pi_prev(solver%num_elements), source=0.0d0)
        allocate(self%dln_nj(solver%num_gas), source=0.0d0)
        allocate(self%dnj_c(solver%num_condensed), source=0.0d0)
        self%dln_n = 0.0d0
        self%dln_T = 0.0d0

        ! Allocate solution variables
        allocate(self%mole_fractions(solver%num_products), source=0.0d0)
        allocate(self%mass_fractions(solver%num_products), source=0.0d0)

    end function

    subroutine EqSolution_set_nj(self, solver, nj_init)
        ! Set the initial mole fractions of the solution

        ! Arguments
        class(EqSolution), intent(inout), target :: self
        type(EqSolver), intent(in), target :: solver
        real(dp), intent(in) :: nj_init(:)
        integer :: i
        real(dp), parameter :: smalno = 1.0d-6
        real(dp), parameter :: smnol = -13.815511d0

        ! Check if nj_init is allocated and has the correct size
        if (size(nj_init) == solver%num_products) then
            self%nj = nj_init
            do i = 1, solver%num_gas
                if (self%nj(i) <= 0.0d0) then
                    self%nj(i) = smalno
                    self%ln_nj(i) = smnol
                else
                    self%ln_nj(i) = log(self%nj(i))
                end if
            end do
            self%n = sum(self%nj(:solver%num_gas))
        else
            call abort('EqSolution_set_nj: nj_init must be allocated and have size equal to num_products.')
        end if

    end subroutine

    subroutine EqSolution_activate_condensed(self, idx, rank)
        class(EqSolution), intent(inout) :: self
        integer, intent(in) :: idx
        integer, intent(in), optional :: rank
        integer :: i, na, target_rank

        if (idx < 1 .or. idx > size(self%is_active)) then
            call abort('EqSolution_activate_condensed: idx out of bounds.')
        end if

        if (self%is_active(idx)) call self%deactivate_condensed(idx)

        na = count(self%is_active)
        target_rank = 1
        if (present(rank)) target_rank = rank
        target_rank = max(1, min(target_rank, na+1))

        do i = 1, size(self%is_active)
            if (self%is_active(i) .and. self%active_rank(i) >= target_rank) then
                self%active_rank(i) = self%active_rank(i) + 1
            end if
        end do

        self%is_active(idx) = .true.
        self%active_rank(idx) = target_rank
    end subroutine

    subroutine EqSolution_activate_condensed_front(self, idx)
        class(EqSolution), intent(inout) :: self
        integer, intent(in) :: idx
        call self%activate_condensed(idx, 1)
    end subroutine

    subroutine EqSolution_deactivate_condensed(self, idx)
        class(EqSolution), intent(inout) :: self
        integer, intent(in) :: idx
        integer :: i, old_rank

        if (idx < 1 .or. idx > size(self%is_active)) then
            call abort('EqSolution_deactivate_condensed: idx out of bounds.')
        end if

        if (.not. self%is_active(idx)) return

        old_rank = max(1, self%active_rank(idx))
        self%is_active(idx) = .false.
        self%active_rank(idx) = 0

        do i = 1, size(self%is_active)
            if (self%is_active(i) .and. self%active_rank(i) > old_rank) then
                self%active_rank(i) = self%active_rank(i) - 1
            end if
        end do
    end subroutine

    subroutine EqSolution_replace_active_condensed(self, old_idx, new_idx)
        class(EqSolution), intent(inout) :: self
        integer, intent(in) :: old_idx, new_idx
        integer :: old_rank

        old_rank = 1
        if (old_idx >= 1 .and. old_idx <= size(self%is_active)) then
            if (self%is_active(old_idx)) old_rank = max(1, self%active_rank(old_idx))
        end if

        if (old_idx >= 1 .and. old_idx <= size(self%is_active)) then
            call self%deactivate_condensed(old_idx)
        end if
        call self%activate_condensed(new_idx, old_rank)
    end subroutine

    function EqSolution_active_condensed_indices(self) result(active_idx)
        class(EqSolution), intent(in) :: self
        integer, allocatable :: active_idx(:)
        integer, allocatable :: used(:)
        integer :: na, nc, i, r, next_slot

        na = count(self%is_active)
        nc = size(self%is_active)
        allocate(active_idx(na))
        if (na == 0) return
        active_idx = 0

        do i = 1, nc
            if (.not. self%is_active(i)) cycle
            r = self%active_rank(i)
            if (r >= 1 .and. r <= na) then
                if (active_idx(r) == 0) active_idx(r) = i
            end if
        end do

        if (any(active_idx == 0)) then
            allocate(used(nc), source=0)
            do i = 1, na
                if (active_idx(i) > 0) used(active_idx(i)) = 1
            end do
            next_slot = 1
            do i = 1, nc
                if (.not. self%is_active(i)) cycle
                if (used(i) /= 0) cycle
                do while (next_slot <= na .and. active_idx(next_slot) /= 0)
                    next_slot = next_slot + 1
                end do
                if (next_slot > na) exit
                active_idx(next_slot) = i
            end do
        end if
    end function

    function EqSolution_num_equations(self, solver) result(num_equations)
        ! Compute the number of equations in the current equilibrium problem

        ! Arguments
        class(EqSolution), intent(in), target :: self
        type(EqSolver), intent(in), target :: solver

        ! Result
        integer :: num_equations

        ! Locals
        integer :: na                         ! Number of active condensed species
        integer :: ne                         ! Number of elements
        logical :: const_p, const_t           ! Flags enabling/disabling matrix equations
        type(EqConstraints), pointer :: cons  ! Abbreviation for soln%constraints

        ! Shorthand
        ne = solver%num_active_elements()
        na = count(self%is_active)
        cons => self%constraints
        const_p = cons%is_constant_pressure()
        const_t = cons%is_constant_temperature()

        ! Get the number of equations
        num_equations = ne + na
        if (.not. const_t) num_equations = num_equations + 1
        if (const_p) num_equations = num_equations + 1

    end function

    function EqSolution_calc_pressure(self) result(pressure)
        ! Calculate the pressure of the mixture

        ! Arguments
        class(EqSolution), intent(in), target :: self

        ! Result
        real(dp) :: pressure  ! Mixture pressure (bar)

        if (self%constraints%is_constant_pressure()) then
            pressure = self%constraints%state2
        else
            pressure = 1.d-5 * R * self%n * self%T / self%constraints%state2
        end if

    end function

    function EqSolution_calc_volume(self) result(volume)
        ! Calculate the specific volume of the mixture

        ! Arguments
        class(EqSolution), intent(in), target :: self

        ! Result
        real(dp) :: volume  ! Mixture specific volume (m**3/kg)

        if (self%constraints%is_constant_pressure()) then
            volume = (1.d-5 * R * self%n * self%T) / self%constraints%state2
        else
            volume = self%constraints%state2
        end if

    end function

    function EqSolution_calc_entropy_sum(self, solver) result(S)
        ! Calculate the total entropy of the product mixture
        ! NOTE: This lives here because it depends on solver state (ln_nj, soln%n,
        !       and pressure during iterations) to preserve numerical behavior.

        ! Arguments
        class(EqSolution), intent(in), target :: self
        type(EqSolver), intent(in), target :: solver

        ! Result
        real(dp) :: S  ! Mixture entropy

        ! Locals
        integer  :: ng                       ! Number of gas species
        real(dp) :: n                        ! Total moles of mixture
        real(dp) :: P                        ! Pressure of mixture (bar)
        real(dp), pointer :: nj(:)           ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)        ! Log of gas species concentrations [kmol-per-kg]
        real(dp), pointer :: s_g(:), s_c(:)  ! Gas/condensed entropies [unitless]

        ! Shorthand
        ng = solver%num_gas
        n = self%n
        nj => self%nj
        ln_nj => self%ln_nj
        s_g => self%thermo%entropy(:ng)
        s_c => self%thermo%entropy(ng+1:)
        P = self%calc_pressure()

        S = dot_product(nj(:ng), s_g-ln_nj-log(P/n)) + dot_product(nj(ng+1:), s_c)

    end function


    !-----------------------------------------------------------------------
    ! EqPartials Implementation
    !-----------------------------------------------------------------------
    function EqPartials_init(num_elements, num_active) result(self)
        ! Construct a new EqPartials object
        type(EqPartials) :: self
        integer, intent(in) :: num_elements
        integer, intent(in) :: num_active

        ! NOTE: EqPartials only need to be computed after convergence, so
        !       "num_active" is used for sizing instead of "num_condensed"
        !       because it is known and fixed at this point

        ! Allocate the partials
        allocate(self%dpi_dlnT(num_elements), &
                 self%dnc_dlnT(num_active), &
                 self%dpi_dlnP(num_elements), &
                 self%dnc_dlnP(num_active))
    end function

    subroutine EqPartials_assemble_partials_matrix_const_p(self, solver, soln, J)
        ! Assemble the matrix for evaluating derivatives with respect to
        ! log(T) at constant P (RP-1311 Table 2.3)

        ! Arguments
        class(EqPartials), intent(in) :: self
        class(EqSolver), intent(in), target :: solver
        class(EqSolution), intent(in), target :: soln
        real(dp), intent(out), allocatable :: J(:,:)

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: na                          ! Number of active condensed species
        integer  :: ne                          ! Number of elements
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: tmp(solver%num_gas)         ! Common sub-expression storage
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp), pointer :: h_g(:), h_c(:)     ! Gas/condensed enthalpies [unitless]
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices
        integer :: r, c                         ! Iteration matrix row/column indices
        integer :: i, k, ic                     ! Loop counters
        integer, allocatable :: active_idx(:)   ! Active condensed indices in legacy order

        ! Define shorthand
        ng = solver%num_gas
        nc = solver%num_condensed
        ne = solver%num_active_elements()
        na = count(soln%is_active)
        active_idx = soln%active_condensed_indices()
        num_eqn = ne+na+1

        ! Associate subarray pointers
        allocate(J(num_eqn, num_eqn+1))
        A_g => solver%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => solver%products%stoich_matrix(ng+1:,:)
        nj  => soln%nj
        nj_g => soln%nj(:ng)
        ln_nj => soln%ln_nj
        h_g => soln%thermo%enthalpy(:ng)
        h_c => soln%thermo%enthalpy(ng+1:)

        ! Initialize the iteration matrix
        J = 0.0d0
        r = 0
        c = 0

        !-------------------------------------------------------
        ! Equation (2.56)
        !-------------------------------------------------------
        do ic = 1,ne
            tmp = nj_g*A_g(:,ic)
            r = r+1
            c = 0

            ! ∂𝛑_i/∂lnT
            do k = 1,ne
                c = c+1
                J(r,c) = dot_product(tmp, A_g(:,k))
            end do

            ! ∂n,c_i/∂lnT
            do k = 1,na
                i = active_idx(k)
                c = c+1
                J(r,c) = A_c(i, ic)
                J(c,r) = J(r,c)  ! Symmetric
            end do

            ! ∂ln(n)/∂lnT
            c = c+1
            J(r,c) = sum(tmp)
            J(c,r) = J(r,c)  ! Symmetric

            ! Right hand side
            c = c+1
            J(r,c) = -dot_product(tmp, h_g)
        end do

        !-------------------------------------------------------
        ! Equation (2.57)
        !-------------------------------------------------------
        do k = 1,na
            i = active_idx(k)
            r = r+1

            ! Right hand size
            J(r,c) = -h_c(i)
        end do

        !-------------------------------------------------------
        ! Equation (2.58)
        !-------------------------------------------------------
        ! Right hand side
        r = r+1
        J(r,c) = -dot_product(nj_g, h_g)

    end subroutine

    subroutine EqPartials_assemble_partials_matrix_const_t(self, solver, soln, J)
        ! Assemble the matrix for evaluating derivatives with respect to
        ! log(P) at constant T (RP-1311 Table 2.4)

        ! Arguments
        class(EqPartials), intent(in) :: self
        class(EqSolver), intent(in), target :: solver
        class(EqSolution), intent(in), target :: soln
        real(dp), intent(out), allocatable :: J(:,:)

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: na                          ! Number of active condensed species
        integer  :: ne                          ! Number of elements
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: tmp(solver%num_gas)           ! Common sub-expression storage
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices
        integer :: r, c                         ! Iteration matrix row/column indices
        integer :: i, k, ic                     ! Loop counters
        integer, allocatable :: active_idx(:)   ! Active condensed indices in legacy order

        ! Define shorthand
        ng = solver%num_gas
        nc = solver%num_condensed
        ne = solver%num_active_elements()
        na = count(soln%is_active)
        active_idx = soln%active_condensed_indices()
        num_eqn = ne+na+1

        ! Associate subarray pointers
        allocate(J(num_eqn, num_eqn+1))
        A_g => solver%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => solver%products%stoich_matrix(ng+1:,:)
        nj  => soln%nj
        nj_g => soln%nj(:ng)
        ln_nj => soln%ln_nj

        ! Initialize the iteration matrix
        J = 0.0d0
        r = 0
        c = 0

        !-------------------------------------------------------
        ! Equation (2.64)
        !-------------------------------------------------------
        do ic = 1,ne
            tmp = nj_g*A_g(:,ic)
            r = r+1
            c = 0

            ! ∂𝛑_i/∂lnP
            do k = 1,ne
                c = c+1
                J(r,c) = dot_product(tmp, A_g(:,k))
            end do

            ! ∂n,c_i/∂lnP
            do k = 1,na
                i = active_idx(k)
                c = c+1
                J(r,c) = A_c(i, ic)
                J(c,r) = J(r,c)  ! Symmetric
            end do

            ! ∂ln(n)/∂lnP
            c = c+1
            J(r,c) = sum(tmp)
            J(c,r) = J(r,c)  ! Symmetric

            ! Right hand side
            c = c+1
            J(r,c) = sum(tmp)
        end do

        !-------------------------------------------------------
        ! Equation (2.65)
        !-------------------------------------------------------
        ! *None*
        r = r+na

        !-------------------------------------------------------
        ! Equation (2.58)
        !-------------------------------------------------------
        ! Right hand side
        r = r+1
        J(r,c) = sum(nj_g)

    end subroutine

    subroutine EqPartials_compute_partials(self, solver, soln)

        ! Arguments
        class(EqPartials), intent(inout), target :: self
        class(EqSolver), intent(in), target :: solver
        class(EqSolution), intent(inout), target :: soln

        ! Locals
        real(dp), allocatable :: J(:,:)
        real(dp) :: nj_solid ! Temporary variables for condensed species
        integer :: ng, ne, nc, na, ierr, i, idx, liq_rank
        integer, allocatable :: active_idx(:)
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: cp(:)              ! Species heat capacity [unitless]
        real(dp), pointer :: h_g(:), h_c(:)     ! Gas/condensed enthalpies [unitless]
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices

        ! Shorthand
        ng = solver%num_gas
        nc = solver%num_condensed
        ne = solver%num_active_elements()
        na = count(soln%is_active)
        A_g => solver%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => solver%products%stoich_matrix(ng+1:,:)
        nj  => soln%nj
        nj_g => soln%nj(:ng)
        cp  => soln%thermo%cp
        h_g => soln%thermo%enthalpy(:ng)
        h_c => soln%thermo%enthalpy(ng+1:)

        ! Compute the equilibrium heat capacity (Eq. 2.59)
        ! cp,eq = ∑_i^ne (∑_j^ng a_ij*n_j*H_j/RT) (∂𝛑_i/∂lnT)_P +  ∑_j^nc H_j/RT (∂n_i/∂lnT)_P +
        !         (∑_j^ng n_j*H_j/RT)(∂ln(n)/∂lnT)_P + ∑_j^ns n_j*Cp_j/R + ∑_j^ng n_j*H_j^2/(RT^2)

        ! Initialize values
        self%cp_eq = 0.0d0

        ! Term 4: ∑_j^ns n_j*Cp_j/R
        self%cp_eq = self%cp_eq + dot_product(soln%thermo%cp, soln%nj)
        soln%cp_fr = dot_product(soln%thermo%cp, soln%nj)*(R/1.d3)

        ! If both solid and liquid present, temporarily remove liquid to prevent singular matrix
        if (soln%j_sol /= 0) then
            call log_info("Solid and liquid phases present; removing liquid for partial derivative calculation")
            nj_solid = soln%nj(ng+soln%j_sol)
            soln%nj(ng+soln%j_sol) = soln%nj(ng+soln%j_sol) + soln%nj(ng+soln%j_liq)
            soln%nj(ng+soln%j_liq) = 0.0d0
            liq_rank = max(1, soln%active_rank(soln%j_liq))
            call soln%deactivate_condensed(soln%j_liq)
            na = count(soln%is_active)
            self%dlnV_dlnT = 0.0d0
            self%cp_eq = 0.0d0
        else
            ! Solve the constant pressure partials
            ierr = 0
            call self%assemble_partials_matrix_const_p(solver, soln, J)
            call gauss(J, ierr)
            if (ierr == 0) then
                self%dpi_dlnT = 0.0d0
                if (ne > 0) self%dpi_dlnT(:ne) = J(:ne, ne+na+2)
                self%dnc_dlnT = J(ne+1:ne+na, ne+na+2)
                self%dn_dlnT = J(ne+na+1, ne+na+2)
                self%dlnV_dlnT = 1.0d0 + self%dn_dlnT

                ! Term 1: ∑_i^ne (∑_j^ng a_ij*n_j*H_j/RT) (∂𝛑_i/∂lnT)_P
                do i = 1,ne
                    self%cp_eq = self%cp_eq + dot_product(nj_g*A_g(:,i), h_g)*self%dpi_dlnT(i)
                end do

                ! Term 2: ∑_j^nc H_j/RT (∂n_i/∂lnT)_P
                active_idx = soln%active_condensed_indices()
                do idx = 1, size(active_idx)
                    i = active_idx(idx)
                    self%cp_eq = self%cp_eq + h_c(i)*self%dnc_dlnT(idx)
                end do

                ! Term 3: (∑_j^ng n_j*H_j/RT)(∂ln(n)/∂lnT)_P
                self%cp_eq = self%cp_eq + dot_product(nj_g, h_g)*self%dn_dlnT

                ! Term 5: ∑_j^ng n_j*H_j^2/(RT^2)
                self%cp_eq = self%cp_eq + dot_product(nj_g*h_g, h_g)

            else
                call log_info("Singular matrix encountered for constant pressure partial derivatives")
                self%dlnV_dlnP = -1.0d0
                self%dlnV_dlnT = 1.0d0
                soln%cp_eq = self%cp_eq*(R/1.d3)
                self%gamma_s = -1.0d0/(self%dlnV_dlnP+(self%dlnV_dlnT**2)*soln%n/self%cp_eq)
                call log_warning('Singular update matrix encountered solver for constant temperature partial derivatives')
            end if
        end if

        ! Solve the constant temperature partials
        ierr = 0
        call self%assemble_partials_matrix_const_t(solver, soln, J)
        call gauss(J, ierr)
        if (ierr == 0) then
            self%dpi_dlnP = 0.0d0
            if (ne > 0) self%dpi_dlnP(:ne) = J(:ne, ne+na+2)
            self%dnc_dlnP = J(ne+1:ne+na, ne+na+2)
            self%dn_dlnP = J(ne+na+1, ne+na+2)
            self%dlnV_dlnP = -1.0d0 + self%dn_dlnP

            ! γ_s := (∂ln(P)/∂ln(ρ))_s
            if (soln%j_liq == 0) then
                self%gamma_s = -1./(self%dlnV_dlnP+(self%dlnV_dlnT**2)*soln%n/self%cp_eq)
            else
                self%gamma_s = -1.0d0/self%dlnV_dlnP
                soln%nj(ng+soln%j_liq) = soln%nj(ng+soln%j_sol) - nj_solid
                soln%nj(ng+soln%j_sol) = nj_solid
                call soln%activate_condensed(soln%j_liq, liq_rank)
            end if

        else
            call log_info("Singular matrix encountered for constant temperature partial derivatives")
            self%dlnV_dlnP = -1.0d0
            self%gamma_s = -1.0d0/(self%dlnV_dlnP+(self%dlnV_dlnT**2)*soln%n/self%cp_eq)
            call log_warning('Singular update matrix encountered solver for constant temperature partial derivatives')
        end if

        ! Store gamma_s in the EqSolution too
        soln%gamma_s = self%gamma_s

        ! Final scaling of heat capacity
        soln%cp_eq = self%cp_eq*(R/1.d3)

        ! Compute Cv
        soln%cv_eq = -soln%cp_eq/(self%gamma_s * self%dlnV_dlnP)

    end subroutine


    !-----------------------------------------------------------------------
    ! Helper Functions
    !-----------------------------------------------------------------------
    subroutine gauss(G, ierr)
        ! Gaussian elimination solver
        !
        ! Solves linear system A*x = b via Gaussian elimination with
        ! partial pivoting. Input is the augmented system G = [A, b].
        ! Solution is done in-place, and on output the augmented
        ! system matrix contains the solution vector in the last column.

        ! Inputs
        real(dp), intent(inout) :: G(:, :)
        integer,  intent(out), optional :: ierr

        ! Locals
        integer :: i, j, k, n, nrow
        real(dp) :: tmp

        nrow = size(G,1)
        if (present(ierr)) ierr = 0

        ! Forward elimination
        do n = 1,nrow-1

            ! Perform partial pivoting
            i = find_pivot(G,n)
            if (i == 0) then
                if (present(ierr)) ierr = n
                return
            end if
            do j = n,nrow+1
                tmp = G(i,j)
                G(i,j) = G(n,j)
                G(n,j) = tmp
            end do

            ! Use pivot to elimate n-th unknown from remaining equations
            G(n,n+1:) = G(n,n+1:)/G(n,n)
            do i = n+1,nrow
                do j = n+1,nrow+1
                    G(i,j) = G(i,j) - G(i,n)*G(n,j)
                end do
            end do

        end do

        ! Backsolve for the variables
        G(nrow,nrow:) = G(nrow,nrow:)/G(nrow,nrow)
        do k = nrow-1,1,-1
            do i = k+1,nrow
                G(k,nrow+1) = G(k,nrow+1) - G(k,i)*G(i,nrow+1)
            end do
        end do

    end subroutine

    function find_pivot(G,n) result(ipivot)
        ! Locate best row, i, to use as pivot for the n-th unknown.

        ! This appears to use a modified form of "classic" partial pivoting.
        ! Instead of simply using the row with largest |G(i,n)| as the pivot,
        ! it uses the row where max(|G(i,n:)|)/|G(i,n)| is the smallest, i.e.,
        ! use the row where x(n) has the greatest influence relative to the
        ! other unknowns in the problem.

        real(dp), intent(in) :: G(:,:)   ! Augmented linear system
        integer,  intent(in) :: n        ! Equation index to be eliminated
        integer :: ipivot                ! Index of pivot row for elimination

        real(dp), parameter :: bigno = 1.d25
        real(dp) :: gn, row_ratio, min_row_ratio
        integer :: nrow,i,j

        nrow = size(G,1)
        ipivot = 0
        min_row_ratio = bigno

        do i = n, nrow
            gn = abs(G(i,n))
            row_ratio = bigno

            ! Find the largest influence ratio along the row
            if (gn /= 0.0d0) then
                row_ratio = 0.0d0
                do j = n+1, nrow+1
                    row_ratio = max(row_ratio, abs(G(i,j)))
                end do
                row_ratio = row_ratio/gn
            end if

            ! Use if lower than current minimum
            if (row_ratio < min_row_ratio) then
                min_row_ratio = row_ratio
                ipivot = i
            end if

        end do

    end function

    function trim_phase(name) result(trim_name)
        ! Take a species name and return the name with the phase suffix removed
        ! e.g., "H2O(L)" -> "H2O"

        character(snl), intent(in) :: name
        character(:), allocatable :: trim_name

        integer :: idx1, idx2

        idx1 = scan(name, '(', back=.true.)
        idx2 = scan(name, ')', back=.true.)
        if (idx2 == len_trim(name)) then
            trim_name = trim(name(1:idx1-1))
        else
            trim_name = trim(name)
        end if

    end function

    function get_species_other_phases(name, products) result(index_list)
        ! Take a species name and return a list of indices to all other instances
        ! of the same species with a different phase. The resulting index list is
        ! for the species_names array.

        character(*), intent(in) :: name
        type(Mixture), intent(in) :: products
        integer, allocatable :: index_list(:), tmp(:)

        ! Locals
        character(:), allocatable :: trim_name, test_name
        integer :: i, j, num_match, ng
        real(dp) :: T_melt_i, T_melt_j

        ! Shorthand
        ng = products%num_gas

        ! If "name" has a species phase suffix, remove it
        trim_name = trim_phase(name)

        allocate(tmp(size(products%species_names(ng+1:))), &
                 index_list(size(products%species_names(ng+1:))))
        num_match = 0
        do i = 1, size(products%species_names(ng+1:))

            ! Keep the current phase in this list to preserve legacy phase-pair logic.
            ! Downstream checks rely on seeing the active phase while iterating candidates.

            test_name = trim_phase(products%species_names(ng+i))
            if (trim_name == test_name) then
                num_match = num_match + 1
                tmp(num_match) = i
            end if

        end do
        tmp = tmp(:num_match)

        ! Sort the species by melting temperature, low to high
        do i = 1, num_match
            T_melt_i = maxval(products%species(ng+tmp(i))%T_fit(:, 2))
            do j = 1, num_match
                if (i == j) cycle
                T_melt_j = maxval(products%species(ng+tmp(j))%T_fit(:, 2))
                if (T_melt_i < T_melt_j) then
                    index_list(i) = tmp(j)
                    tmp(j) = tmp(i)
                    tmp(i) = index_list(i)
                end if
            end do
        end do
        index_list = tmp(:num_match)

    end function


    !-----------------------------------------------------------------------
    ! Transport property calculations
    !-----------------------------------------------------------------------

    subroutine compute_transport_properties(eq_solver, eq_soln, frozen_shock)
        ! Compute the transport properties of a mixture for a given equilibrium solution

        ! Arguments
        class(EqSolver), target :: eq_solver
        type(EqSolution), intent(inout) :: eq_soln
        logical, intent(in), optional :: frozen_shock  ! TODO: Estimate mole fractions for a frozen shock problem

        ! Locals
        real(dp), allocatable :: psi(:,:)     ! Viscosity interaction matrix
        real(dp), allocatable :: phi(:,:)     ! Conductivity interaction matrix
        real(dp), allocatable :: eta(:,:)     ! Binary interaction matrix
        real(dp), allocatable :: cond(:)      ! Conductivity array
        integer, allocatable :: idx_list(:)   ! List of indices (in solver order) to compute transport properties for
        integer, allocatable :: pure_idx(:)   ! List of indices (into transport order) to compute transport properties for
        integer, allocatable :: bin_idx(:)    ! List of indices (into transport order) to compute transport properties for
        integer :: bin_count                  ! Total number of binary pairs to consider
        integer :: ng                         ! Number of gas species
        integer :: ne                         ! Number of elements
        integer :: np, nb                     ! Number of pure, binary species
        integer :: nm                         ! Number of of gaseous species for thermal transport calculations
        integer :: nr                         ! Number of chemical reactions
        integer :: i, j, k, k1, k2, ii, m     ! Index counters
        integer :: idx(1), idx1(1), idx2(1)   ! Temporary findloc index
        integer :: max_elem_idx               ! Number of elements (minus electron, if applicaple)
        real(dp) :: cfit_val                  ! Value computed by curve-fit function
        real(dp) :: te, ekt, qc, xsel, debye, ionic, lambda  ! Variables for ionized species interactions
        real(dp), parameter :: tol = 1.d-8    ! Tolerance to test if a value is ~ zero
        real(dp) :: test_tot, test_nj
        real(dp) :: nj_el                     ! Electron concentration
        real(dp) :: nj_cutoff                 ! Minimum species concentration to include
        real(dp) :: viscns                    ! Viscosity constant
        real(dp), pointer :: A(:, :)          ! Stoichiometric matrix
        real(dp) :: wmol, wmol1, wmol2        ! Species molecular weight
        real(dp), allocatable :: cp(:)        ! Heat capacity of each species in the mixture
        real(dp) :: omega
        logical :: ion1, ion2, elc1, elc2     ! Flags for ionized interactions
        real(dp) :: wmred                     ! Reduced molecular weight
        real(dp) :: ratio                     ! Ratio of molecular weights
        real(dp) :: sumc, sumv                ! Temporary sum for conductivity/viscosity calculations
        real(dp) :: total                     ! Total moles of species with transport properties
        real(dp), allocatable :: xs(:)        ! Mole fraction of species with transport properties
        real(dp), allocatable, target :: G(:, :)  ! Matrix to compute reaction properties
        real(dp), allocatable :: gmat(:, :)   ! Temporary storage of G values
        real(dp), allocatable :: delh(:)      ! Reaction enthalpy
        real(dp), allocatable :: rtpd(:, :)   ! ??
        real(dp), allocatable :: xsij(:, :)   ! Cross-product of mole fractions
        real(dp), allocatable :: alpha(:, :)  ! Stoichiometrix matrix for the chemical reactions
        real(dp), allocatable :: stcf(:, :)   ! Stores some coefficients
        real(dp), allocatable :: stcoef(:)    ! Stores some coefficients
        logical :: change                     ! Flag to switch value in stoichiometric matrix
        real(dp) :: coeff                     ! Temporary coefficient value
        integer, allocatable :: tmp(:)        ! Temporrary indexing array
        real(dp), allocatable :: stx(:), stxij(:, :)  !
        real(dp), pointer :: x(:)             ! Solution vector
        real(dp) :: cpreac, cp_eq             ! Reaction, equilibrium heat capacity
        real(dp) :: reacon                    ! Reaction conductivity
        integer :: ierr                       ! Gauss solver error index
        real(dp) :: wtmol                     ! Total molecular weight
        integer, parameter :: max_tr = 40     ! Maximum allowable transport species

        ! Define shorthand
        np = eq_solver%transport_db%num_pure
        nb = eq_solver%transport_db%num_binary
        ng = eq_solver%num_gas
        ne = eq_solver%num_elements
        A => eq_solver%products%stoich_matrix

        ! Allocate
        allocate(psi(ng, ng), phi(ng, ng), eta(ng, ng), cond(ng), &
                 idx_list(max_tr), pure_idx(np), bin_idx(nb), &
                 cp(ng), xs(ng), G(ng, ng), rtpd(ng, ng), &
                 xsij(ng, ng), delh(ng), alpha(ng, ng), &
                 stcf(ng, ng), stcoef(ng), tmp(max_tr), gmat(ng, ng), &
                 stx(ng), stxij(ng, ng))

        ! cond can be used without being initialized, and uninitialized elements
        ! could be used later if all of the species aren't found in the transport database
        cond = 0.0d0

        ! Build the list of relevant mixture species, starting with monoatomic gasses
        nm = 0
        total = 0.0d0
        wtmol = 1.0/sum(eq_soln%nj)
        nj_cutoff = 1.d-11/wtmol
        test_tot = 0.999999999d0/wtmol
        max_elem_idx = eq_solver%num_elements
        if (eq_solver%ions) max_elem_idx = max_elem_idx - 1
        do i = 1, ng
            ! Check if this is a monoatomic gas
            if ((sum(abs(A(i, :)))-1.0d0) < tol) then
                if (eq_soln%nj(i) <= 0.0d0) then
                    if (eq_soln%ln_nj(i) - log(eq_soln%n) + eq_solver%xsize > 0.0d0) then
                        eq_soln%nj(i) = exp(eq_soln%ln_nj(i))
                    end if
                end if
                nm = nm + 1
                idx_list(nm) = i
                total = total + eq_soln%nj(i)
                if (eq_solver%products%species(i)%molecular_weight < 1.0d0) nj_el = eq_soln%nj(i)
                eq_soln%nj(i) = -eq_soln%nj(i)
            end if
        end do
        test_nj = 1.0d0/(ng*wtmol)

        ! Add the remaining species that meet the minimum size threshold
        do i = 1, ng
            if (total <= test_tot .and. nm < max_tr) then
                test_nj = test_nj / 10.0d0
                do j = 1, ng
                    if (eq_soln%nj(j) >= test_nj) then
                        if (nm >= max_tr) then
                            call log_info("Reached maximum number of allowable transport species.")
                            exit
                        else
                            total = total + eq_soln%nj(j)
                            nm = nm + 1
                            idx_list(nm) = j
                            eq_soln%nj(j) = -eq_soln%nj(j)
                        end if
                    end if
                end do
                if (test_nj < nj_cutoff) then
                    exit
                end if
            else
                exit
            end if
        end do
        idx_list = idx_list(:nm)

        ! Undo the negative species concentrations
        do i = 1, ng
            if (eq_soln%nj(i) < 0.0d0) then
                eq_soln%nj(i) = -eq_soln%nj(i)
            end if
        end do

        ! Build the list of pure species index
        j = 0
        do i = 1, nm
            idx = findloc(eq_solver%transport_db%pure_species, eq_solver%products%species_names(idx_list(i)))
            if (idx(1) > 0) then
                j = j + 1
                pure_idx(idx(1)) = i
            else
                call log_info('compute_transport_properties: Species '//eq_solver%products%species_names(idx_list(i))//&
                              ' not found in transport database.')
            end if
        end do
        pure_idx = pure_idx(:j)
        np = j

        ! Remove any binary pairs with negligible concentrations
        bin_count = 0
        do i = 1, nb
            idx1 = findloc(eq_solver%products%species_names, eq_solver%transport_db%binary_species(i,1))
            idx2 = findloc(eq_solver%products%species_names, eq_solver%transport_db%binary_species(i,2))
            if (idx1(1) > 0 .and. idx2(1) > 0) then
                if (eq_soln%nj(idx1(1)) > 0.0d0 .and. eq_soln%nj(idx2(1)) > 0.0d0) then
                    bin_count = bin_count + 1
                    bin_idx(bin_count) = i
                end if
            else
                call log_info('compute_transport_properties: Binary species'//eq_solver%transport_db%binary_species(i,1)// &
                              ' or '//eq_solver%transport_db%binary_species(i,2)//'not found in products mixture.')
            end if

        end do
        bin_idx = bin_idx(:bin_count)

        ! Compute moles of species with transport data
        xs = eq_soln%nj(idx_list)/total

        ! --------------------------------------------------------------
        ! Build the eta matrix
        ! --------------------------------------------------------------

        ! Find transport data for important interactions
        eta = eta(:nm, :nm)
        eta = 0.0d0
        do i = 1, bin_count
            idx1 = findloc(eq_solver%transport_db%pure_species, eq_solver%transport_db%binary_species(bin_idx(i), 1))
            idx2 = findloc(eq_solver%transport_db%pure_species, eq_solver%transport_db%binary_species(bin_idx(i), 2))
            if (idx1(1) > 0 .and. idx2(1) > 0) then
                cfit_val = exp(eq_solver%transport_db%binary_transport(bin_idx(i))%calc_eta(eq_soln%T))
                eta(pure_idx(idx1(1)), pure_idx(idx2(1))) = cfit_val
                eta(pure_idx(idx2(1)), pure_idx(idx1(1))) = cfit_val
            else
                call log_info('compute_transport_properties: Binary species'//eq_solver%transport_db%binary_species(bin_idx(i),1)//&
                ' or '//eq_solver%transport_db%binary_species(bin_idx(i),2)//'not found in transport database')
            end if
        end do

        ! Add the diagonal terms
        do i = 1, np
            eta(pure_idx(i),pure_idx(i)) = exp(eq_solver%transport_db%pure_transport(i)%calc_eta(eq_soln%T))
        end do

        ! Build the conductivity array
        cond = cond(:nm)
        do i = 1, np
            cond(pure_idx(i)) = exp(eq_solver%transport_db%pure_transport(i)%calc_lambda(eq_soln%T))
        end do

        ! Build the stoichiometrix matrix for the chemical reactions
        alpha = alpha(:, :nm)
        alpha = 0.0d0
        nr = nm - ne
        k = 1
        do i = (ne+1), nm
            alpha(k, i) = -1.0d0
            j = idx_list(i)
            do m = 1, ne
                alpha(k, m) = A(j, m)
            end do
            k = k + 1
        end do

        do i = 1, nm
            if (xs(i) < 1.d-10) then
                m = 1
                change = .false.
                do j = 1, nr
                    coeff = alpha(j, i)
                    if (abs(coeff) > 1.d-5) then
                        if (.not. change) then
                            change = .true.
                            do k = 1, nm
                                stcoef(k) = alpha(j, k)/coeff
                            end do
                            continue
                        else
                            do k = 1, nm
                                alpha(j, k) = (alpha(j, k)/coeff) - stcoef(k)
                            end do
                        end if
                    end if
                    do k = 1, nm
                        stcf(m, k) = alpha(j, k)
                    end do
                    m = m + 1
                end do
                do ii = 1, nm
                    do j = 1, nr
                        alpha(j, ii) = stcf(j, ii)
                    end do
                end do
                nr = m - 1
            end if
        end do
        alpha = alpha(:nr, :nm)

        ! Make estimates for interactions with missing data
        if (eq_solver%ions) then
            te = eq_soln%T/1000.d0
            ekt = 4.8032d0**2/(Boltz*te)
            qc = 100.d0*(ekt**2)
            xsel = nj_el/total
            IF ( xsel < 1.0d-12 ) xsel = 1.0d-12
            debye = ((22.5d0/Pi)*(R/Avgdr*100.d0)*(te/xsel))/ekt**3
            ionic = ((810.d0/(4.0d0*Pi))*(R/Avgdr*100d0)*(te/xsel))**(2.0/3.0)/ekt**2
            lambda = sqrt(debye+ionic)
            lambda = max(lambda, 2.71828183d0)
        end if

        ! Fill in missing diagonals
        viscns = 0.3125*sqrt(1.d5*Boltz/(pi*Avgdr))
        cp = eq_soln%thermo%cp(idx_list)
        cp = cp(:nm)
        do i = 1, nm
            k = idx_list(i)
            if (.not. (eq_solver%ions .and. abs(A(k, ne)-1.0d0) < tol) &
                .and. abs(eta(i, i)) < tol) then
                if (abs(eta(i, i)) < tol) then
                    wmol = eq_solver%products%species(k)%molecular_weight
                    omega = log(50.0d0*wmol**4.6/eq_soln%T**1.4)
                    omega = max(omega, 1.0d0)
                    eta(i, i) = viscns*sqrt(wmol*eq_soln%T)/omega
                end if
                if (abs(cond(i)) < tol) cond(i) = &
                    eta(i, i)*R*(0.00375d0 + 0.00132d0*(cp(i) - 2.5d0))/wmol
            end if
        end do

        ! Fill in missing off-diagonals
        do i = 1, nm
            k1 = idx_list(i)
            wmol1 = eq_solver%products%species(k1)%molecular_weight
            do j = 1, nm
                ion1 = .false.
                ion2 = .false.
                elc1 = .false.
                elc2 = .false.
                omega = 0.0d0
                if (eta(i, j) == 0.0d0) eta(i, j) = eta(j, i)
                if (eta(j, i) == 0.0d0) eta(j, i) = eta(i, j)
                if (abs(eta(i, j)) < tol) then
                    k2 = idx_list(j)
                    wmol2 = eq_solver%products%species(k2)%molecular_weight
                    if (eq_solver%ions) then
                        ! Estimate for ions
                        if (abs(A(k1, ne)) == 1.0d0) ion1 = .true.
                        if (abs(A(k2, ne)) == 1.0d0) ion2 = .true.
                        if (wmol1 < 1.0d0) elc1 = .true.
                        if (wmol2 < 1.0d0) elc2 = .true.
                        if (ion1 .and. ion2) omega = 1.36d0*qc*log(lambda)
                        if ((ion1 .and. elc2) .or. (ion2 .and. elc1)) &
                            omega = 1.29d0*qc*log(lambda)
                        if ((ion1 .and. .not. ion2) .or. (ion2 .and. .not. ion1)) &
                            omega =  exp(6.776-0.4*log(eq_soln%T))
                        if (abs(omega) > tol) then
                            wmred = sqrt(2.0*eq_soln%T*wmol1*wmol2/(wmol1+wmol2))
                            eta(i, j) = viscns*wmred*pi/omega
                            eta(j, i) = eta(i, j)
                            if (i == j) then
                                cond(i) = eta(i, i)*R*(0.00375d0 + 0.00132d0*(cp(i) - 2.5d0))/wmol1
                            end if
                        else
                            ratio = sqrt(wmol2/wmol1)
                            eta(i, j) = 5.656854d0*eta(i, i)*sqrt(wmol2/(wmol1 + wmol2))
                            eta(i, j) = eta(i, j)/(1.d0 + sqrt(ratio*eta(i, i)/eta(j, j)))**2
                            eta(j, i) = eta(i, j)
                        end if
                    else
                        ! Estimate for unlike interactions from rigid sphere analogy
                        ratio = sqrt(wmol2/wmol1)
                        eta(i, j) = 5.656854d0*eta(i, i)*sqrt(wmol2/(wmol1 + wmol2))
                        eta(i, j) = eta(i, j)/(1.d0 + sqrt(ratio*eta(i, i)/eta(j, j)))**2
                        eta(j, i) = eta(i, j)
                    end if
                end if
            end do
        end do

        ! --------------------------------------------------------------
        ! Calculate viscosity and frozen thermal conductivity
        ! --------------------------------------------------------------

        ! Build phi, psi matrices
        phi = phi(:nm, :nm)
        psi = psi(:nm, :nm)
        do i = 1,nm
            rtpd(i, i) = 0.0d0
            phi(i, i) = 1.0d0
            psi(i, i) = 1.0d0
        end do

        do i = 1, (nm-1)
            k1 = idx_list(i)
            wmol1 = eq_solver%products%species(k1)%molecular_weight
            do j = (i+1), nm
                k2 = idx_list(j)
                wmol2 = eq_solver%products%species(k2)%molecular_weight
                sumc = 2.d0 / (eta(i, j)*(wmol1+wmol2))
                phi(i, j) = sumc*wmol2*eta(i, i)
                phi(j, i) = sumc*wmol1*eta(j, j)
                sumc = (wmol1 + wmol2)**2
                psi(i, j) = phi(i, j) * (1.d0 + 2.41d0*(wmol1 - wmol2) * &
                            (wmol1 - 0.142d0*wmol2)/sumc)
                psi(j, i) = phi(j, i) * (1.d0 + 2.41d0*(wmol2 - wmol1) * &
                            (wmol2 - 0.142d0*wmol1)/sumc)
            end do
        end do

        ! Calculate viscsosity and frozen conductivity
        eq_soln%viscosity = 0.0d0
        eq_soln%conductivity_fr = 0.0d0
        do i = 1, nm
            sumc = 0.0d0
            sumv = 0.0d0
            do j = 1, nm
                sumc = sumc + psi(i, j)*xs(j)
                sumv = sumv + phi(i, j)*xs(j)
            end do
            eq_soln%viscosity = eq_soln%viscosity + eta(i,i)*xs(i)/sumv
            eq_soln%conductivity_fr = eq_soln%conductivity_fr + cond(i)*xs(i)/sumc
        end do
        eq_soln%viscosity = eq_soln%viscosity/1.d3  ! Set viscosity to millipoise
        eq_soln%conductivity_fr = eq_soln%conductivity_fr/1.d3  ! Set conductivity to W/m-K

        ! --------------------------------------------------------------
        ! Calculate reaction heat capacity and thermal conductivity
        ! --------------------------------------------------------------

        delh = delh(:nr)
        G = G(:nr, :nr+1)
        do i = 1, nr
            delh(i) = 0.0d0
            do j = 1, nm
                delh(i) = delh(i) + alpha(i, j)*eq_soln%thermo%enthalpy(idx_list(j))
            end do
            G(i, nr+1) = delh(i)
        end do

        do i = 1, nr
            do j = 1, nm
                if (abs(alpha(i, j)) < 1.d-6) alpha(i, j) = 0.0d0
            end do
        end do

        xsij = xsij(:nm, :nm)
        rtpd = rtpd(:nm, :nm)
        do i = 1, (nm-1)
            k1 = idx_list(i)
            wmol1 = eq_solver%products%species(k1)%molecular_weight
            do j = (i+1), nm
                k2 = idx_list(j)
                wmol2 = eq_solver%products%species(k2)%molecular_weight
                rtpd(i, j) = wmol1*wmol2/(1.1d0 * eta(i, j) * (wmol1 + wmol2))
                xsij(i, j) = xs(i)*xs(j)
                xsij(j, i) = xsij(i, j)
                rtpd(j, i) = rtpd(i, j)
            end do
        end do

        do i = 1, nr
            do j = 1, nr
                G(i, j) = 0.0d0
                gmat(i, j) = 0.0d0
            end do
        end do

        do k = 1, (nm-1)
            do m = (k+1), nm
                if (xs(k) >= 1.d-10 .and. xs(m) >= 1.d-10) then
                    do j = 1, nr
                        if ((alpha(j, k) == 0.0d0) .and. (alpha(j, m) == 0.0d0)) then
                            stx(j) = 0.0d0
                        else
                            stx(j) = xs(m)*alpha(j, k) - xs(k)*alpha(j, m)
                        end if
                    end do
                    do i = 1, nr
                        do j = 1, nr
                            stxij(i, j) = stx(i)*stx(j)/xsij(k, m)
                            G(i, j) = G(i, j) + stxij(i, j)
                            gmat(i, j) = gmat(i, j) + stxij(i, j)*rtpd(k, m)
                        end do
                    end do
                end if
            end do
        end do

        m = 1 + nr
        do i = 1, nr
            do j = 1, nr
                G(j, i) = G(i, j)
            end do
            G(i, m) = delh(i)
        end do

        call gauss(G, ierr)
        x => G(:, m)

        cpreac = 0.0d0
        do i = 1, nr
            cpreac = cpreac + (R*1.d-3)*delh(i)*x(i)
            G(i, m) = delh(i)  ! *** "x" is a pointer, so this updates "x(i)" as well ***
            do j = i, nr
                G(i, j) = gmat(i, j)
                G(j, i) = G(i, j)
            end do
        end do

        call gauss(G, ierr)
        x => G(:, m)

        reacon = 0.0d0
        do i = 1, nr
            reacon = reacon + (R*1.d-3)*delh(i)*x(i)
        end do
        reacon = 0.6d0*reacon

        wtmol = 0.0d0
        eq_soln%cp_fr = 0.0d0
        do i = 1, nm
            eq_soln%cp_fr = eq_soln%cp_fr + cp(i)*xs(i)
            wtmol = wtmol + xs(i)*eq_solver%products%species(idx_list(i))%molecular_weight
        end do
        eq_soln%cp_fr = eq_soln%cp_fr*(R*1.d-3)/wtmol

        ! Compute the remaining properties
        eq_soln%pr_fr = eq_soln%viscosity*eq_soln%cp_fr/eq_soln%conductivity_fr
        cpreac = cpreac/wtmol
        cp_eq = eq_soln%cp_fr + cpreac
        eq_soln%conductivity_eq = eq_soln%conductivity_fr + (reacon*1.d-3)
        eq_soln%pr_eq = eq_soln%viscosity*cp_eq/eq_soln%conductivity_eq

    end subroutine

end module
