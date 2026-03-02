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
        logical :: smooth_truncation = .false.
            !! Enable smooth logistic gating instead of hard truncation
        real(dp) :: truncation_width = 0.25d0
            !! Logistic gate width in log-space (w); only used when smooth_truncation is true
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
        procedure :: update_transport_basis => EqSolver_update_transport_basis
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

        ! Inputs
        real(dp), allocatable :: w0(:)
            !! Reactant weights

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
        real(dp) :: T_seed = 0.0d0
            !! Last stable warm-start temperature seed
        real(dp) :: n_seed = 0.0d0
            !! Last stable warm-start total-moles seed
        real(dp), allocatable :: nj_seed(:)
            !! Last stable warm-start species concentrations
        real(dp), allocatable :: ln_nj_seed(:)
            !! Last stable warm-start log gas concentrations
        logical, allocatable :: is_active_seed(:)
            !! Last stable warm-start active condensed flags
        integer, allocatable :: active_rank_seed(:)
            !! Last stable warm-start condensed ordering
        integer :: j_liq_seed = 0
            !! Last stable warm-start liquid index
        integer :: j_sol_seed = 0
            !! Last stable warm-start solid index
        integer :: j_switch_seed = 0
            !! Last stable warm-start removed condensed index
        integer :: last_cond_idx_seed = 0
            !! Last stable warm-start added condensed index

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

        ! Legacy-style transport component basis cached after convergence.
        integer :: transport_basis_rows = 0
            !! Number of active basis rows used by transport reaction assembly
        integer, allocatable :: transport_component_idx(:)
            !! Component species indices (global gas-species indexing) for each basis row
        real(dp), allocatable :: transport_basis_matrix(:, :)
            !! Row-reduced element-by-gas coefficient matrix used by transport reactions

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

        type :: EqDerivatives
        !! Equilibrium Total Derivatives Type
        !!
        !! Computes the total derivatives of the solution variables x with
        !! respect to the input variables u using the direct method, where
        !! x: [P0/V0, T0/H0/S0/U0, b0] \in R^n
        !! u: [𝛑, nc, ln(n), ln(T)] \in R^m,(ln(T) ommitted for const T, ln(n) ommitted for const V)

        !! Size variables
        integer :: m = 0
            !! Number of equations in the matrix system = number of equations in the active set
        integer :: n = 0
            !! Number of variables in the solution vector = number of elements + 2

        !! Solver workspace
        real(dp), allocatable :: R(:)
            !! Nonlinear residuals (m)
        real(dp), allocatable :: J(:, :)
            !! Jacobian matrix of the nonlinear residuals (m x m)
        real(dp), allocatable :: Rx(:, :)
            !! Partial derivatives of the nonlinear residuals wrt inputs (m x n)
        real(dp), allocatable :: dudx(:, :)
            !! Total derivatives of the solution variables wrt inputs (m x n)
        real(dp), allocatable :: delta_check(:, :)
            !! Delta = J*dudx + Rx = 0, used to check the correctness of the computed derivatives (m x n)

        !! Final unpacked derivatives
        real(dp) :: dT_dstate1
            !! Total derivative of T wrt state1 (T0/H0/S0/U0)
        real(dp) :: dT_dstate2
            !! Total derivative of T wrt state2 (P0/V0)
        real(dp), allocatable :: dT_dw0(:)
            !! Total derivative of T wrt input weights

        real(dp) :: dn_dstate1
            !! Total derivative of n wrt state1 (T0/H0/S0/U0)
        real(dp) :: dn_dstate2
            !! Total derivative of n wrt state2 (P0/V0)
        real(dp), allocatable :: dn_dw0(:)
            !! Total derivative of n wrt input weights

        real(dp), allocatable :: dnj_dstate1(:)
            !! Total derivative of species concentrations wrt state1 (T0/H0/S0/U0)
        real(dp), allocatable :: dnj_dstate2(:)
            !! Total derivative of species concentrations wrt state2 (P0/V0)
        real(dp), allocatable :: dnj_dw0(:,:)
            !! Total derivative of species concentrations wrt input weights

        real(dp) :: dH_dstate1
            !! Total derivative of enthalpy wrt state1 (T0/H0/S0/U0)
        real(dp) :: dH_dstate2
            !! Total derivative of enthalpy wrt state2 (P0/V0)
        real(dp), allocatable :: dH_dw0(:)
            !! Total derivative of enthalpy wrt input weights

        real(dp) :: dU_dstate1
            !! Total derivative of energy wrt state1 (T0/H0/S0/U0)
        real(dp) :: dU_dstate2
            !! Total derivative of energy wrt state2 (P0/V0)
        real(dp), allocatable :: dU_dw0(:)
            !! Total derivative of energy wrt input weights

        real(dp) :: dG_dstate1
            !! Total derivative of Gibbs free energy wrt state1 (T0/H0/S0/U0)
        real(dp) :: dG_dstate2
            !! Total derivative of Gibbs free energy wrt state2 (P0/V0)
        real(dp), allocatable :: dG_dw0(:)
            !! Total derivative of Gibbs free energy wrt input weights

        real(dp) :: dS_dstate1
            !! Total derivative of entropy wrt state1 (T0/H0/S0/U0)
        real(dp) :: dS_dstate2
            !! Total derivative of entropy wrt state2 (P0/V0)
        real(dp), allocatable :: dS_dw0(:)
            !! Total derivative of entropy wrt input weights

        real(dp) :: dCp_fr_dstate1
            !! Total derivative of frozen heat capacity wrt state1 (T0/H0/S0/U0)
        real(dp) :: dCp_fr_dstate2
            !! Total derivative of frozen heat capacity wrt state2 (P0/V0)
        real(dp), allocatable :: dCp_fr_dw0(:)
            !! Total derivative of frozen heat capacity wrt input weights

        !! Finite-difference derivatives (for verification)
        real(dp) :: dT_dstate1_fd
        real(dp) :: dT_dstate2_fd
        real(dp), allocatable :: dT_dw0_fd(:)

        real(dp) :: dn_dstate1_fd
        real(dp) :: dn_dstate2_fd
        real(dp), allocatable :: dn_dw0_fd(:)

        real(dp), allocatable :: dnj_dstate1_fd(:)
        real(dp), allocatable :: dnj_dstate2_fd(:)
        real(dp), allocatable :: dnj_dw0_fd(:,:)

        real(dp) :: dH_dstate1_fd
        real(dp) :: dH_dstate2_fd
        real(dp), allocatable :: dH_dw0_fd(:)

        real(dp) :: dU_dstate1_fd
        real(dp) :: dU_dstate2_fd
        real(dp), allocatable :: dU_dw0_fd(:)

        real(dp) :: dG_dstate1_fd
        real(dp) :: dG_dstate2_fd
        real(dp), allocatable :: dG_dw0_fd(:)

        real(dp) :: dS_dstate1_fd
        real(dp) :: dS_dstate2_fd
        real(dp), allocatable :: dS_dw0_fd(:)

        real(dp) :: dCp_fr_dstate1_fd
        real(dp) :: dCp_fr_dstate2_fd
        real(dp), allocatable :: dCp_fr_dw0_fd(:)

    contains

        procedure :: assemble_jacobian => EqDerivatives_assemble_jacobian
        procedure :: assemble_Rx => EqDerivatives_assemble_Rx
        procedure :: compute_residual => EqDerivatives_compute_residual
        procedure :: check_closure_defect => EqDerivatives_check_closure_defect
        procedure :: unpack_values => EqDerivatives_unpack_values
        procedure :: compute_derivatives => EqDerivatives_compute_derivatives
        procedure :: compute_fd => EqDerivatives_compute_fd

    end type
    interface EqDerivatives
        module procedure :: EqDerivatives_init
    end interface

contains

    pure subroutine sigmoid_stable(x, g, log_g)
        ! Numerically stable sigmoid and log-sigmoid.
        real(dp), intent(in) :: x
        real(dp), intent(out) :: g
        real(dp), intent(out), optional :: log_g
        real(dp) :: ex

        if (x >= 0.0d0) then
            ex = exp(-x)
            g = 1.0d0 / (1.0d0 + ex)
            if (present(log_g)) log_g = -log(1.0d0 + ex)
        else
            ex = exp(x)
            g = ex / (1.0d0 + ex)
            if (present(log_g)) log_g = x - log(1.0d0 + ex)
        end if
    end subroutine

    pure subroutine compute_nj_effective(ln_nj, ln_threshold, smooth_truncation, width, nj_eff, g, dg_dln, ln_nj_eff, &
                                         dln_nj_eff_dln_nj, is_hard_active)
        ! Map log-species amount to the effective physical amount used in thermo/properties.
        ! In this solver, ln_threshold = log(n) - tsize.
        real(dp), intent(in) :: ln_nj
        real(dp), intent(in) :: ln_threshold
        logical, intent(in) :: smooth_truncation
        real(dp), intent(in) :: width
        real(dp), intent(out) :: nj_eff
        real(dp), intent(out), optional :: g
        real(dp), intent(out), optional :: dg_dln
        real(dp), intent(out), optional :: ln_nj_eff
        real(dp), intent(out), optional :: dln_nj_eff_dln_nj
        logical, intent(out), optional :: is_hard_active

        real(dp) :: gate
        real(dp) :: gate_log
        real(dp) :: gate_prime
        real(dp) :: x
        logical :: hard_active

        hard_active = (ln_nj > ln_threshold)
        if (present(is_hard_active)) is_hard_active = hard_active

        if (.not. smooth_truncation) then
            if (hard_active) then
                nj_eff = exp(ln_nj)
            else
                nj_eff = 0.0d0
            end if
            if (present(g)) then
                if (hard_active) then
                    g = 1.0d0
                else
                    g = 0.0d0
                end if
            end if
            if (present(dg_dln)) dg_dln = 0.0d0
            if (present(ln_nj_eff)) ln_nj_eff = ln_nj
            if (present(dln_nj_eff_dln_nj)) dln_nj_eff_dln_nj = 1.0d0
            return
        end if

        x = (ln_nj - ln_threshold) / width
        call sigmoid_stable(x, gate, gate_log)
        gate_prime = (1.0d0 / width) * gate * (1.0d0 - gate)

        nj_eff = exp(ln_nj) * gate
        if (present(g)) g = gate
        if (present(dg_dln)) dg_dln = gate_prime
        if (present(ln_nj_eff)) ln_nj_eff = ln_nj + gate_log
        if (present(dln_nj_eff_dln_nj)) then
            if (gate > 0.0d0) then
                dln_nj_eff_dln_nj = 1.0d0 + gate_prime / gate
            else
                dln_nj_eff_dln_nj = 1.0d0
            end if
        end if
    end subroutine

    pure real(dp) function gas_amount_ln_threshold(ln_n, tsize, esize, ion_species) result(ln_threshold)
        ! Species-amount truncation threshold: ions use esize, non-ions use tsize.
        real(dp), intent(in) :: ln_n
        real(dp), intent(in) :: tsize
        real(dp), intent(in) :: esize
        logical, intent(in) :: ion_species

        if (ion_species) then
            ln_threshold = ln_n - esize
        else
            ln_threshold = ln_n - tsize
        end if
    end function

    !-----------------------------------------------------------------------
    ! EquilibriumSolver
    !-----------------------------------------------------------------------
    function EqSolver_init(products, reactants, trace, ions, all_transport, insert, &
            smooth_truncation, truncation_width) result(self)
        type(EqSolver) :: self
        type(Mixture), intent(in) :: products
        type(Mixture), intent(in), optional :: reactants
        real(dp), intent(in), optional :: trace
        logical, intent(in), optional :: ions
        type(TransportDB), intent(in), optional :: all_transport
        character(*), intent(in), optional :: insert(:)  ! List of condensed species to insert
        logical, intent(in), optional :: smooth_truncation
        real(dp), intent(in), optional :: truncation_width
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
        if (present(smooth_truncation)) self%smooth_truncation = smooth_truncation
        if (present(truncation_width)) self%truncation_width = truncation_width

        if (self%smooth_truncation .and. self%truncation_width <= 0.0d0) then
            call abort("EqSolver_init: truncation_width must be > 0 when smooth_truncation is enabled.")
        end if

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
        real(dp) :: log_n                             ! Log of total moles
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
        log_n = log(n)

        ! Associate subarray pointers
        nj_g => soln%nj(:ng)

        ! Compute lambda1 (Eq. 3.1) and lambda2 (Eq. 3.2)
        ! TODO: what happens for a TV problem when these are both not defined?
        l1_denom = max(5.0d0*abs(dln_T), 5.0d0*abs(dln_n))

        lambda1 = 1.0d0
        lambda2 = 1.0d0
        do i = 1, ng
            if (dln_nj(i) > 0.0d0) then !.and. (nj_g(i) > 0.0d0)) then
                if (ln_nj(i) - log_n + self%size <= 0.0d0) then
                    l2_denom = abs(dln_nj(i) - dln_n)
                    if (l2_denom >= (self%size + FACTOR)) then
                        temp_l2 = abs(FACTOR - ln_nj(i) + log_n)/l2_denom
                        lambda2 = min(lambda2, temp_l2)
                    end if
                else if (dln_nj(i) > l1_denom) then
                    l1_denom = dln_nj(i)
                end if
            end if
        end do
        if (l1_denom > 2.0d0) lambda1 = 2.0d0/l1_denom

        ! Compute lambda (Eq. 3.3)
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
        real(dp) :: ln_n                      ! Log total moles
        real(dp) :: ln_threshold              ! Truncation threshold in log-space
        real(dp), pointer :: ln_nj(:)         ! Log of the product concentrations
        real(dp) :: ln_nj_eff(self%num_gas)   ! Effective log-species concentrations
        real(dp) :: dln_nj_eff_dln_nj(self%num_gas)  ! d(ln(nj_eff))/d(ln_nj)
        real(dp) :: nj_eff_tmp                ! Temporary effective species amount
        logical :: ion_species                ! True if gas species is charged and ions are active

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
        ln_n = log(n)
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
        do i = 1, ng
            ion_species = self%ions .and. self%active_ions .and. ne_full > 0 .and. A_g(i, ne_full) /= 0.0d0
            ln_threshold = gas_amount_ln_threshold(ln_n, self%tsize, self%esize, ion_species)
            call compute_nj_effective(ln_nj(i), ln_threshold, self%smooth_truncation, self%truncation_width, &
                                      nj_eff=nj_eff_tmp, ln_nj_eff=ln_nj_eff(i), &
                                      dln_nj_eff_dln_nj=dln_nj_eff_dln_nj(i))
        end do
        mu_g = h_g - s_g + ln_nj_eff + log(P/n)

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
            if (self%smooth_truncation) then
                soln%dln_nj(i) = soln%dln_nj(i) / max(dln_nj_eff_dln_nj(i), tiny(1.0d0))
            end if

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
        logical :: ion_species                ! True if gas species is charged and ions are active
        logical :: const_p, const_t           ! Flags enabling/disabling matrix equations
        type(EqConstraints), pointer :: cons  ! Abbreviation for soln%constraints
        real(dp) :: lambda                    ! Damped update factor
        real(dp) :: ln_threshold              ! Truncation threshold in log-space

        allocate(active_idx(0))

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
        ln_threshold = ln_n - self%tsize
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

        ! Update gas species concentration. Use esize for ions and tsize otherwise.
        do i = 1, ng
            ln_nj(i) = ln_nj(i) + lambda*dln_nj(i)
            ion_species = self%ions .and. self%active_ions .and. ne_full > 0 .and. A_g(i, ne_full) /= 0.0d0
            ln_threshold = gas_amount_ln_threshold(ln_n, self%tsize, self%esize, ion_species)
            call compute_nj_effective(ln_nj(i), ln_threshold, self%smooth_truncation, self%truncation_width, &
                                      nj_eff=nj_g(i))
        end do

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
        real(dp) :: ln_threshold              ! Truncation threshold in log-space
        real(dp) :: ln_nj_eff(self%num_gas)   ! Effective log-species concentrations
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
        real(dp) :: sum1, sum2, aa, temp_raw, temp_eff, gate, x_gate
            ! Temporary variables for ionized species
        real(dp) :: nj_eff_tmp                    ! Temporary species amount
        real(dp) :: sum_nj                        ! Total species amount used in convergence checks
        logical :: ion_species                    ! True if gas species is charged and ions are active
        logical :: hard_active                    ! Legacy hard-threshold activity flag

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
        sum_nj = sum(nj)
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
        if (const_s) then
            do i = 1, ng
                ion_species = self%ions .and. self%active_ions .and. ne_full > 0 .and. A_g(i, ne_full) /= 0.0d0
                ln_threshold = gas_amount_ln_threshold(ln_n, self%tsize, self%esize, ion_species)
                call compute_nj_effective(ln_nj(i), ln_threshold, self%smooth_truncation, self%truncation_width, &
                                          nj_eff=nj_eff_tmp, ln_nj_eff=ln_nj_eff(i))
            end do
            s_delta = cons%state1 - dot_product(nj_g, s_g-ln_nj_eff-log(P/n)) - dot_product(nj(ng+1:), s_c)
        end if

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
            ion_species = self%ions .and. self%active_ions .and. ne_full > 0 .and. A_g(i, ne_full) /= 0.0d0
            ln_threshold = gas_amount_ln_threshold(ln_n, self%tsize, self%esize, ion_species)
            if (self%smooth_truncation) then
                hard_active = (ln_nj(i) > ln_threshold)
                if (.not. hard_active) cycle
            end if
            if ((nj_g(i)*abs(dln_nj(i))/sum_nj) > nj_tol) then
                soln%gas_converged = .false.
                return
            end if
        end do

        ! Check condensed species updates
        do i = 1, na
            if (abs(dnj_c(i))/sum_nj > nj_tol) then
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
            ion_species = self%ions .and. self%active_ions .and. ne_full > 0 .and. A_g(i, ne_full) /= 0.0d0
            ln_threshold = gas_amount_ln_threshold(ln_n, self%tsize, self%esize, ion_species)
            call compute_nj_effective(ln_nj(i), ln_threshold, self%smooth_truncation, self%truncation_width, &
                                      nj_eff=nj_g(i))
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
                        ! NOTE(smooth_truncation): Preserve legacy hard-threshold semantics
                        ! when smooth truncation is disabled. Smooth mode uses
                        ! sigmoid-gated amounts in this ion-convergence loop.
                        temp_raw = 0.0d0
                        if (soln%ln_nj(j) > -87.0d0) temp_raw = exp(soln%ln_nj(j))

                        if (self%smooth_truncation) then
                            ln_threshold = gas_amount_ln_threshold(ln_n, self%tsize, self%esize, .true.)
                            x_gate = (soln%ln_nj(j) - ln_threshold) / self%truncation_width
                            call sigmoid_stable(x_gate, gate)
                            temp_eff = temp_raw*gate
                            soln%nj(j) = temp_eff
                        else
                            temp_eff = temp_raw
                            soln%nj(j) = 0.0d0
                            ln_threshold = gas_amount_ln_threshold(ln_n, self%tsize, self%esize, .true.)
                            if (soln%ln_nj(j) > ln_threshold) then
                                soln%nj(j) = temp_raw
                            end if
                        end if

                        aa = A_g(j, ne_full)*temp_eff
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
        integer  :: ne_full                     ! Total number of elements (including electron)
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: tmp(self%num_gas)           ! Common sub-expression storage
        real(dp) :: mu_g(self%num_gas)          ! Gas phase chemical potentials [unitless]
        real(dp) :: b_delta(self%num_elements)  ! Residual for element contraints
        real(dp) :: n_delta                     ! Residual for total moles / pressure constraint
        real(dp) :: hsu_delta                   ! Residual for enthalpy / entropy constraint
        real(dp) :: n                           ! Total moles of mixture
        real(dp) :: P                           ! Pressure of mixture (bar)
        real(dp), pointer :: nj(:)             ! Total species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp) :: ln_nj_eff(self%num_gas)     ! Effective log-species concentrations
        real(dp) :: nj_eff_g(self%num_gas)      ! Smooth-mapped gas species concentrations
        real(dp) :: nj_linear(self%num_gas)     ! d(nj_eff)/d(ln_nj) weights used in Newton linearization
        real(dp) :: dln_nj_eff_dln_nj(self%num_gas)  ! d(ln(nj_eff))/d(ln_nj)
        real(dp) :: ln_threshold                ! Truncation threshold in log-space
        real(dp) :: ln_n                        ! Log of total moles
        real(dp) :: nj_eval(self%num_products)  ! Physical amounts consistent with smooth mapping
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
        logical :: ion_species                  ! True if gas species is charged and ions are active
        logical :: const_p, const_t, const_s, const_h, const_u  ! Flags enabling/disabling matrix equations
        type(EqConstraints), pointer :: cons    ! Abbreviation for soln%constraints

        allocate(active_idx(0))

        ! Define shorthand
        ng = self%num_gas
        nc = self%num_condensed
        ne = self%num_active_elements()
        ne_full = self%num_elements
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
        ln_n = log(n)
        nj  => soln%nj
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
        do i = 1, ng
            ion_species = self%ions .and. self%active_ions .and. ne_full > 0 .and. A_g(i, ne_full) /= 0.0d0
            ln_threshold = gas_amount_ln_threshold(ln_n, self%tsize, self%esize, ion_species)
            call compute_nj_effective(ln_nj(i), ln_threshold, self%smooth_truncation, self%truncation_width, &
                                      nj_eff=nj_eff_g(i), ln_nj_eff=ln_nj_eff(i), &
                                      dln_nj_eff_dln_nj=dln_nj_eff_dln_nj(i))
            nj_linear(i) = nj_eff_g(i)
        end do
        mu_g = h_g - s_g + ln_nj_eff + log(P/n)
        nj_eval = nj
        nj_eval(:ng) = nj_eff_g

        ! Evalutate constraint residuals
        b_delta = cons%b0 - self%products%elements_from_species(nj_eval)
        n_delta = n - sum(nj_eff_g)
        if (const_s) then
            hsu_delta = (cons%state1 - soln%calc_entropy_sum(self))
        else if (const_h) then
            hsu_delta = (cons%state1/soln%T - dot_product(nj_eval, soln%thermo%enthalpy))
        else if (const_u) then
            hsu_delta = (cons%state1/soln%T - dot_product(nj_eval, soln%thermo%energy))
        end if

        ! Initialize the iteration matrix
        G = 0.0d0
        r = 0
        c = 0

        !-------------------------------------------------------
        ! Equation (2.24/2.45): Element constraints
        !-------------------------------------------------------
        do i = 1,ne
            tmp = nj_linear*A_g(:,i)
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
                G(r,c) = dot_product(nj_linear, h_g)
            end if

            ! Right-hand-side
            G(r,c+1) = n_delta + dot_product(nj_linear, mu_g)

        end if

        !---------------------------------------------------------
        ! Equation (2.27)/(2.28)/(2.47)/(2.48): Energy constraints
        !---------------------------------------------------------
        if (.not. const_t) then
            r = r+1
            c = 0

            ! Select entropy/enthalpy constraint
            if (const_s) then
                tmp = nj_linear*(h_g-mu_g)
                h_or_s_or_u => soln%thermo%entropy(ng+1:)
            else if (const_h) then
                tmp = nj_linear*h_g
                h_or_s_or_u => soln%thermo%enthalpy(ng+1:)
            else if (const_u) then
                tmp = nj_linear*u_g
                h_or_s_or_u => soln%thermo%energy(ng+1:)
            end if

            ! Pi derivatives
            do j = 1,ne
                c = c+1
                G(r,c) = dot_product(tmp, A_g(:,j))
                if (.not. const_p .and. const_s) then
                    G(r,c) = G(r,c) - dot_product(nj_linear, A_g(:, j))
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
                G(r,c) = dot_product(nj_eval, cp) + dot_product(tmp, h_g)
            else
                G(r,c) = dot_product(nj_eval, cv) + dot_product(tmp, u_g)
                if (const_s) then
                    G(r,c) = G(r,c) - dot_product(nj_linear, u_g)
                end if
            end if

            ! Right-hand-side
            G(r,c+1) = hsu_delta + dot_product(tmp, mu_g)
            if (const_s) then
                if (const_p) then
                    G(r,c+1) = G(r,c+1) + n_delta
                else
                    G(r,c+1) = G(r,c+1) - dot_product(nj_linear, mu_g)
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
        allocate(active_idx(0))

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

        allocate(active_idx(0))

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

        allocate(active_idx(0))


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
                ! NOTE(smooth_truncation): Legacy singular-recovery path intentionally seeds
                ! only non-positive species to preserve historical restart behavior.
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

    subroutine EqSolver_update_transport_basis(self, soln)
        ! Cache a legacy-style component basis for transport-reaction assembly.
        class(EqSolver), intent(in), target :: self
        type(EqSolution), intent(inout), target :: soln

        integer :: ng, ne, nn
        integer :: i, j, irow, jbx, jex, jb, l, njc
        integer, allocatable :: jx(:), jcm(:), lcs(:)
        real(dp), allocatable :: a_rows(:, :), xs_work(:)
        real(dp), pointer :: A(:, :)
        logical :: newcom, same_col, accept_comp
        real(dp) :: tem
        real(dp), parameter :: tol = 1.d-8
        real(dp), parameter :: smalno = 1.d-10

        ng = self%num_gas
        ne = self%num_active_elements()
        if (ng <= 0 .or. ne <= 0) then
            soln%transport_basis_rows = 0
            if (allocated(soln%transport_component_idx)) soln%transport_component_idx = 0
            if (allocated(soln%transport_basis_matrix)) soln%transport_basis_matrix = 0.0d0
            return
        end if

        A => self%products%stoich_matrix
        nn = ne
        if (self%ions .and. self%active_ions) nn = max(1, ne-1)

        allocate(a_rows(nn, ng), xs_work(ng), jx(nn), jcm(nn), lcs(nn))
        do irow = 1, nn
            do j = 1, ng
                a_rows(irow, j) = A(j, irow)
            end do
        end do

        jx = 0
        jcm = 0
        do irow = 1, nn
            do j = 1, ng
                if (abs(abs(A(j, irow)) - 1.0d0) < tol .and. abs(sum(abs(A(j, :ne))) - 1.0d0) < tol) then
                    jx(irow) = j
                    exit
                end if
            end do
            if (jx(irow) == 0) then
                do j = 1, ng
                    if (abs(a_rows(irow, j)) > smalno) then
                        jx(irow) = j
                        exit
                    end if
                end do
            end if
            jcm(irow) = jx(irow)
        end do

        xs_work = max(soln%nj(:ng), 0.0d0)
        lcs = 0
        njc = 0
        newcom = .false.
        do
            if (njc >= nn) exit
            jbx = maxloc(xs_work(:ng), dim=1)
            if (jbx <= 0) exit
            if (xs_work(jbx) <= 0.0d0) exit

            if (self%ions .and. self%active_ions) then
                if (abs(A(jbx, ne)) > smalno) then
                    xs_work(jbx) = -1.0d0
                    cycle
                end if
            end if

            do irow = 1, nn
                if (jbx == 0) jbx = jx(irow)
                if (jbx == 0) cycle
                if (a_rows(irow, jbx) <= smalno) cycle

                accept_comp = .true.
                if (njc /= 0) then
                    do i = 1, njc
                        l = lcs(i)
                        if (l == irow) then
                            accept_comp = .false.
                            exit
                        end if
                        if (l == 0) cycle
                        jb = jcm(l)
                        if (jb == 0) cycle
                        same_col = .true.
                        do j = 1, nn
                            if (a_rows(j, jbx) /= a_rows(j, jb)) then
                                same_col = .false.
                                exit
                            end if
                        end do
                        if (same_col) then
                            accept_comp = .false.
                            exit
                        end if
                    end do
                end if
                if (.not. accept_comp) cycle

                do i = 1, nn
                    if (i == irow) cycle
                    jex = jx(i)
                    if (jex == 0) cycle
                    if (abs(a_rows(irow, jbx)*a_rows(i, jex) - a_rows(irow, jex)*a_rows(i, jbx)) <= smalno) then
                        accept_comp = .false.
                        exit
                    end if
                end do
                if (.not. accept_comp) cycle

                njc = njc + 1
                if (jbx /= jcm(irow)) newcom = .true.
                jcm(irow) = jbx
                lcs(njc) = irow
                exit
            end do

            xs_work(jbx) = -1.0d0
        end do

        do irow = 1, nn
            if (jcm(irow) == 0) jcm(irow) = jx(irow)
        end do

        if (newcom) then
            do irow = 1, nn
                jb = jcm(irow)
                if (jb == 0) cycle
                if (a_rows(irow, jb) == 0.0d0) then
                    jb = jx(irow)
                    jcm(irow) = jb
                end if
                tem = a_rows(irow, jb)
                if (tem == 0.0d0) cycle
                if (tem /= 1.0d0) a_rows(irow, :ng) = a_rows(irow, :ng)/tem
                do i = 1, nn
                    if (i == irow) cycle
                    if (a_rows(i, jb) /= 0.0d0) then
                        tem = a_rows(i, jb)
                        a_rows(i, :ng) = a_rows(i, :ng) - a_rows(irow, :ng)*tem
                        do j = 1, ng
                            if (abs(a_rows(i, j)) < 1.0d-5) a_rows(i, j) = 0.0d0
                        end do
                    end if
                end do
            end do
        end if

        soln%transport_basis_rows = nn
        soln%transport_component_idx = 0
        soln%transport_basis_matrix = 0.0d0
        soln%transport_component_idx(:nn) = jcm(:nn)
        soln%transport_basis_matrix(:nn, :ng) = a_rows(:nn, :ng)

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

    subroutine EqSolution_reset_iteration_state(soln)
        ! Reset transient Newton-update state before each solve call.
        ! This allows EqSolution instances to be safely reused across solves.
        type(EqSolution), intent(inout) :: soln

        if (allocated(soln%dln_nj)) soln%dln_nj = 0.0d0
        if (allocated(soln%dnj_c)) soln%dnj_c = 0.0d0
        soln%dln_n = 0.0d0
        soln%dln_T = 0.0d0
        soln%dpi_e = 0.0d0

        soln%gas_converged = .false.
        soln%condensed_converged = .false.
        soln%moles_converged = .false.
        soln%element_converged = .false.
        soln%temperature_converged = .false.
        soln%entropy_converged = .false.
        soln%pi_converged = .false.
        soln%ions_converged = .false.
        soln%converged = .false.
    end subroutine

    subroutine EqSolution_save_seed(soln)
        type(EqSolution), intent(inout) :: soln

        if (allocated(soln%nj_seed)) soln%nj_seed = soln%nj
        if (allocated(soln%ln_nj_seed)) soln%ln_nj_seed = soln%ln_nj
        if (allocated(soln%is_active_seed)) soln%is_active_seed = soln%is_active
        if (allocated(soln%active_rank_seed)) soln%active_rank_seed = soln%active_rank

        soln%T_seed = soln%T
        soln%n_seed = soln%n
        soln%j_liq_seed = soln%j_liq
        soln%j_sol_seed = soln%j_sol
        soln%j_switch_seed = soln%j_switch
        soln%last_cond_idx_seed = soln%last_cond_idx
    end subroutine

    subroutine EqSolution_restore_seed(soln)
        type(EqSolution), intent(inout) :: soln

        if (allocated(soln%nj_seed)) soln%nj = soln%nj_seed
        if (allocated(soln%ln_nj_seed)) soln%ln_nj = soln%ln_nj_seed
        if (allocated(soln%is_active_seed)) soln%is_active = soln%is_active_seed
        if (allocated(soln%active_rank_seed)) soln%active_rank = soln%active_rank_seed

        soln%T = soln%T_seed
        soln%n = soln%n_seed
        soln%j_liq = soln%j_liq_seed
        soln%j_sol = soln%j_sol_seed
        soln%j_switch = soln%j_switch_seed
        soln%last_cond_idx = soln%last_cond_idx_seed
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
        logical :: made_change, max_iter_fallback_used, was_converged

        call log_debug("Starting Eq. Solve.")

        ! If the prior solve did not converge, restore the last stable iterate
        ! seed before applying new constraints for this solve call.
        was_converged = soln%converged
        if (.not. was_converged) then
            call log_debug("Restoring last stable warm-start seed after non-converged solve.")
            call EqSolution_restore_seed(soln)
        end if

        ! Set problem type, fixed-state values, and element amounts
        soln%w0 = reactant_weights
        call soln%constraints%set( &
            type, state1, state2, &
            self%reactants%element_amounts_from_weights(reactant_weights) &
        )

        ! If fixed-temperature, set it
        if (soln%constraints%is_constant_temperature()) then
            soln%T = state1
        end if

        ! Reset transient update/convergence fields on every solve so a reused
        ! EqSolution starts from a clean Newton step history.
        call EqSolution_reset_iteration_state(soln)

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
                    soln%converged = .false.
                    call EqSolver_restore_reduced_elements(self, soln, num_reduced, reduced_from, reduced_to)
                    call log_warning('EqSolver_solve: Too many singular matrices encountered.')
                    return
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
                call log_warning("Convergence failed to establish set of condensed species.")
                return
            end if

            ! Initial convergence; check on adding or removing condensed species
            call self%test_condensed(soln, iter, singular_index)

            if (soln%converged .or. (iter == self%max_iterations)) then

                ! Compute final species concentrations
                ! * NOTE: post-processing uses a lower threshold when computing nj = exp(ln(nj))
                do i = 1, self%num_gas
                    if (self%smooth_truncation) then
                        call compute_nj_effective(soln%ln_nj(i), log(soln%n)-self%tsize, self%smooth_truncation, &
                                                  self%truncation_width, nj_eff=soln%nj(i))
                    else
                        if (soln%ln_nj(i) > self%log_min) soln%nj(i) = exp(soln%ln_nj(i))
                    end if
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
                    call log_warning('EqSolver_solve: Maximum iterations reached without convergence')
                    return
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
                if (self%transport) then
                    call self%update_transport_basis(soln)
                    call compute_transport_properties(self, soln)
                end if

                ! Compute post-processing solution values
                call self%post_process(soln, present(partials))

                ! Check for temperature outside of bounds
                if (soln%T > self%T_max .or. soln%T < self%T_min) then
                    call log_warning("Mixture temperature outside of allowable bounds.")
                    soln%converged = .false.
                end if

                if (soln%converged) call EqSolution_save_seed(soln)

                call EqSolver_restore_reduced_elements(self, soln, num_reduced, reduced_from, reduced_to)
                return
            end if

        end do

        call EqSolver_restore_reduced_elements(self, soln, num_reduced, reduced_from, reduced_to)
        soln%converged = .false.

    end subroutine

    !-----------------------------------------------------------------------
    ! EqDerivatives
    !-----------------------------------------------------------------------
    function EqDerivatives_init(solver, solution) result(self)

        ! Arguments
        type(EqSolver), intent(in) :: solver
        type(EqSolution), intent(in) :: solution
        type(EqDerivatives) :: self

        ! Locals
        integer :: m, n, nr, ns

        m = solution%num_equations(solver)
        n = solver%num_elements + 2
        nr = solver%num_reactants
        ns = solver%num_gas + count(solution%is_active)  ! Number of species (gas + active condensed)
        self%m = m
        self%n = n

        allocate(self%R(m), source=empty_dp)
        allocate(self%J(m, m), source=empty_dp)
        allocate(self%Rx(m, n), source=empty_dp)
        allocate(self%dudx(m, n), source=empty_dp)

        allocate(self%delta_check(m, n), source=empty_dp)

        allocate(self%dT_dw0(nr), source=empty_dp)
        allocate(self%dn_dw0(nr), source=empty_dp)
        allocate(self%dnj_dstate1(ns), source=empty_dp)
        allocate(self%dnj_dstate2(ns), source=empty_dp)
        allocate(self%dnj_dw0(ns, nr), source=empty_dp)
        allocate(self%dH_dw0(nr), source=empty_dp)
        allocate(self%dU_dw0(nr), source=empty_dp)
        allocate(self%dG_dw0(nr), source=empty_dp)
        allocate(self%dS_dw0(nr), source=empty_dp)
        allocate(self%dCp_fr_dw0(nr), source=empty_dp)

        allocate(self%dT_dw0_fd(nr), source=empty_dp)
        allocate(self%dn_dw0_fd(nr), source=empty_dp)
        allocate(self%dnj_dstate1_fd(ns), source=empty_dp)
        allocate(self%dnj_dstate2_fd(ns), source=empty_dp)
        allocate(self%dnj_dw0_fd(ns, nr), source=empty_dp)
        allocate(self%dH_dw0_fd(nr), source=empty_dp)
        allocate(self%dU_dw0_fd(nr), source=empty_dp)
        allocate(self%dG_dw0_fd(nr), source=empty_dp)
        allocate(self%dS_dw0_fd(nr), source=empty_dp)
        allocate(self%dCp_fr_dw0_fd(nr), source=empty_dp)

    end function

    subroutine EqDerivatives_assemble_jacobian(self, solver, solution)

        ! Arguments
        class(EqDerivatives), intent(inout) :: self
        type(EqSolver), intent(in) :: solver
        type(EqSolution), intent(inout) :: solution

        call solver%assemble_matrix(solution)
        self%J = solution%G(:self%m, :self%m)
    end subroutine

    subroutine EqDerivatives_assemble_Rx(self, solver, solution)

        ! Arguments
        class(EqDerivatives), intent(inout), target :: self
        type(EqSolver), intent(in), target :: solver
        type(EqSolution), intent(in), target :: solution

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: ne                          ! Number of elements
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: tmp(solver%num_gas)         ! Common sub-expression storage
        real(dp) :: dtmp_dP(solver%num_gas)     ! d/dP of common sub-expression storage
        real(dp) :: mu_g(solver%num_gas)        ! Gas phase chemical potentials [unitless]
        real(dp) :: dhsu_delta_dP               ! d/dP of residual for enthalpy / entropy constraint
        real(dp) :: n                           ! Total moles of mixture
        real(dp) :: P                           ! Pressure of mixture (bar)
        real(dp) :: T                           ! Temperature of mixture (K)
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp) :: ln_nj_eff(solver%num_gas)   ! Effective log-species concentrations
        real(dp) :: ln_threshold                ! Truncation threshold in log-space
        real(dp), pointer :: h_g(:), h_c(:)     ! Gas enthalpies [unitless]
        real(dp), pointer :: s_g(:)             ! Gas entropies [unitless]
        real(dp), pointer :: u_g(:)             ! Gas energies [unitless]
        real(dp) :: dh_g_dT(solver%num_gas)     ! d/dT of gas enthalpies [unitless]
        real(dp) :: ds_g_dT(solver%num_gas)     ! d/dT of gas entropies [unitless]
        real(dp) :: dh_c_dT(solver%num_condensed)  ! d/dT of condensed enthalpies [unitless]
        real(dp) :: ds_c_dT(solver%num_condensed)  ! d/dT of condensed entropies [unitless]
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices
        integer :: r, c                         ! Iteration matrix row/column indices
        integer :: i                            ! Loop counters
        logical :: const_p, const_t, const_s, const_h, const_u  ! Flags enabling/disabling matrix equations
        type(EqConstraints), pointer :: cons    ! Abbreviation for soln%constraints

        ! Define shorthand
        ng = solver%num_gas
        nc = solver%num_condensed
        ne = solver%num_elements
        num_eqn = solution%num_equations(solver)
        cons => solution%constraints
        const_p = cons%is_constant_pressure()
        const_t = cons%is_constant_temperature()
        const_s = cons%is_constant_entropy()
        const_h = cons%is_constant_enthalpy()
        const_u = cons%is_constant_energy()

        ! Associate subarray pointers
        A_g => solver%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => solver%products%stoich_matrix(ng+1:,:)
        n = solution%n
        nj  => solution%nj
        nj_g => solution%nj(:ng)
        ln_nj => solution%ln_nj
        h_g => solution%thermo%enthalpy(:ng)
        h_c => solution%thermo%enthalpy(ng+1:)
        s_g => solution%thermo%entropy(:ng)
        u_g => solution%thermo%energy(:ng)

        ! Get the mixture pressure and temperature
        P = solution%calc_pressure()
        T = solution%T
        ln_threshold = log(n) - solver%tsize

        ! Compute gas phase chemical potentials
        do i = 1, ng
            call compute_nj_effective(ln_nj(i), ln_threshold, solver%smooth_truncation, solver%truncation_width, &
                                      nj_eff=tmp(i), ln_nj_eff=ln_nj_eff(i))
        end do
        mu_g = h_g - s_g + ln_nj_eff + log(P/n)

        ! Evalutate constraint residuals
        dhsu_delta_dP = 0.0d0
        if (const_s .and. const_p) then
            dhsu_delta_dP = sum(nj_g)/P
        end if

        ! Compute intermediate derivatives
        do i = 1, ng
            dh_g_dT(i) = solver%products%species(i)%calc_denthalpy_dT(T)/T - h_g(i)/T
            ds_g_dT(i) = solver%products%species(i)%calc_dentropy_dT(T)
        end do
        do i = 1, nc
            dh_c_dT(i) = solver%products%species(ng+i)%calc_denthalpy_dT(T)/T - h_c(i)/T
            ds_c_dT(i) = solver%products%species(ng+i)%calc_dentropy_dT(T)
        end do

        ! Initialize the iteration matrix
        self%Rx = 0.0d0
        r = 0
        c = 0

        !-------------------------------------------------------
        ! Equation (2.24/2.45): Element constraints
        !-------------------------------------------------------
        do i = 1,ne
            tmp = nj_g*A_g(:,i)
            r = r+1
            c = 0

            ! dR/dx1 (x1: fixed pressure or volume)
            c = c+1
            self%Rx(r, c) = -sum(tmp)/P
            if (.not. const_p) then
                self%Rx(r, c) = self%Rx(r, c)*(-P/cons%state2)
            end if

            ! dR/dx2 (x2: fixed temperature/enthalpy/entropy/energy)
            c = c+1
            if (const_t) then
                self%Rx(r, c) = -dot_product(tmp, dh_g_dT-ds_g_dT)
                if (.not. const_p) then
                    self%Rx(r, c) = self%Rx(r, c) - sum(tmp)/T
                end if
            else
                self%Rx(r, c) = 0.0d0
                if (.not. const_p) then
                    self%Rx(r, c) = self%Rx(r, c) - sum(tmp)/T
                end if
            end if

            ! dR/dx3...n (x3...n: element amounts)
            c = i+2
            self%Rx(r, c) = -1.0d0 ! Other contributions = 0
        end do

        !-------------------------------------------------------
        ! Equation (2.25/2.46): Condensed phase constraints
        !-------------------------------------------------------
        do i = 1,nc
            if (.not. solution%is_active(i)) cycle
            r = r+1
            c = 0

            ! dR/dx1 (x1: fixed pressure or volume) are 0

            ! dR/dx2 (x2: fixed temperature/enthalpy/entropy/energy)
            if (const_t) then
                self%Rx(r, 2) = -dh_c_dT(i) + ds_c_dT(i)
            end if

            ! dR/dx3...n (x3...n: element amounts) are 0

        end do

        !-------------------------------------------------------
        ! Equation (2.26)
        !-------------------------------------------------------
        if (const_p) then
            r = r+1
            c = 0

            ! dR/dx1 (x1: fixed pressure or volume)
            c = c+1
            self%Rx(r, c) = -sum(nj_g)/P
            if (.not. const_p) then
                self%Rx(r, c) = self%Rx(r, c)*(-P/cons%state2)
            end if

            ! dR/dx2 (x2: fixed temperature/enthalpy/entropy/energy)
            c = c+1
            if (const_t) then
                self%Rx(r, c) = -dot_product(nj_g, dh_g_dT-ds_g_dT)
            else
                self%Rx(r, c) = 0.0d0
            end if

            ! dR/dx3...n (x3...n: element amounts) are 0

        end if

        !---------------------------------------------------------
        ! Equation (2.27)/(2.28)/(2.47)/(2.48): Energy constraints
        !---------------------------------------------------------
        if (.not. const_t) then
            r = r+1
            c = 0

            ! Select entropy/enthalpy constraint
            if (const_s) then
                tmp = nj_g*(h_g-mu_g)
            else if (const_h) then
                tmp = nj_g*h_g
            else if (const_u) then
                tmp = nj_g*u_g
            end if

            dtmp_dP = 0.0d0

            ! dR/dx1 (x1: fixed pressure or volume)
            c = c+1
            self%Rx(r, c) = -dhsu_delta_dP - dot_product(dtmp_dP, mu_g) - sum(tmp)/P
            if (.not. const_p) then
                self%Rx(r, c) = self%Rx(r, c)*(-P/cons%state2)
            end if

            ! dR/dx2 (x2: fixed temperature/enthalpy/entropy/energy)
            c = c+1
            if (const_s) then
                self%Rx(r, c) = -1.0d0
            else
                self%Rx(r, c) = -1.0d0/T
            end if

            ! dR/dx3...n (x3...n: element amounts) are 0

        end if

    end subroutine

    subroutine EqDerivatives_compute_residual(self, solver, solution)

        ! Arguments
        class(EqDerivatives), intent(inout), target :: self
        type(EqSolver), intent(in), target :: solver
        type(EqSolution), intent(in), target :: solution

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: ne                          ! Number of elements
        integer  :: r                           ! Residual vector row index
        real(dp) :: b_delta(solver%num_elements)  ! Residual for element constraints
        real(dp) :: n_delta                     ! Residual for total moles / pressure constraint
        real(dp) :: hsu_delta                   ! Residual for enthalpy / entropy constraint
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: h_c(:)             ! Condensed enthalpies [unitless]
        real(dp), pointer :: s_c(:)             ! Condensed entropies [unitless]
        real(dp), pointer :: A_c(:,:)           ! Condensed stoichiometric matrix
        real(dp), pointer :: pi(:)              ! Modified Lagrange multipliers
        integer :: i
        logical :: const_p, const_t, const_s, const_h, const_u
        type(EqConstraints), pointer :: cons

        ! Define shorthand
        ng = solver%num_gas
        nc = solver%num_condensed
        ne = solver%num_elements
        cons => solution%constraints
        const_p = cons%is_constant_pressure()
        const_t = cons%is_constant_temperature()
        const_s = cons%is_constant_entropy()
        const_h = cons%is_constant_enthalpy()
        const_u = cons%is_constant_energy()

        ! Associate subarray pointers
        A_c => solver%products%stoich_matrix(ng+1:, :)
        nj => solution%nj
        nj_g => solution%nj(:ng)
        h_c => solution%thermo%enthalpy(ng+1:)
        s_c => solution%thermo%entropy(ng+1:)
        pi => solution%pi

        ! Evaluate constraint residuals
        b_delta = cons%b0 - solver%products%elements_from_species(nj)
        n_delta = solution%n - sum(nj_g)
        if (const_s) then
            hsu_delta = cons%state1 - solution%calc_entropy_sum(solver)
        else if (const_h) then
            hsu_delta = cons%state1/solution%T - dot_product(nj, solution%thermo%enthalpy)
        else if (const_u) then
            hsu_delta = cons%state1/solution%T - dot_product(nj, solution%thermo%energy)
        else
            hsu_delta = 0.0d0
        end if

        ! Assemble the residual vector
        self%R = 0.0d0
        r = 0

        ! Element residuals
        do i = 1, ne
            r = r + 1
            self%R(r) = b_delta(i)
        end do

        ! Condensed species residuals
        do i = 1, nc
            if (.not. solution%is_active(i)) cycle
            r = r + 1
            self%R(r) = h_c(i) - s_c(i) - dot_product(A_c(i, :), pi)
        end do

        ! Total moles residual
        if (const_p) then
            r = r + 1
            self%R(r) = n_delta
        end if

        ! Energy residual
        if (.not. const_t) then
            r = r + 1
            self%R(r) = hsu_delta
        end if

    end subroutine

    subroutine EqDerivatives_check_closure_defect(self, verbose)
        ! Arguments
        class(EqDerivatives), intent(inout) :: self
        logical, intent(in), optional :: verbose

        ! Locals
        integer :: i
        logical :: lverbose

        lverbose = .false.
        if (present(verbose)) lverbose = verbose

        do i = 1, self%n
            self%delta_check(:, i) = matmul(self%J, self%dudx(:, i)) + self%Rx(:, i)
            if (lverbose) then
                write(*,*) "max|delta| row=", maxloc(abs(self%delta_check(:, i))), " val=", maxval(abs(self%delta_check(:, i)))
                write(*,*) "delta_check(:, ", i, ") = ", self%delta_check(:, i)
            end if
        end do

    end subroutine

    subroutine EqDerivatives_unpack_values(self, solver, solution)

        ! Arguments
        class(EqDerivatives), intent(inout), target :: self
        type(EqSolver), intent(in), target :: solver
        type(EqSolution), intent(in), target :: solution

        ! Locals
        integer  :: ng                          ! Number of gas species
        integer  :: nc                          ! Number of condensed species
        integer  :: na                          ! Number of active condensed species
        integer  :: ne                          ! Number of elements
        integer  :: nr                          ! Number of reactant species
        integer  :: num_eqn                     ! Active number of equations
        real(dp) :: mu_g(solver%num_gas)        ! Gas phase chemical potentials [unitless]
        real(dp) :: n                           ! Total moles of mixture
        real(dp) :: P                           ! Pressure of mixture (bar)
        real(dp) :: T                           ! Temperature of mixture (K)
        real(dp), pointer :: nj(:), nj_g(:)     ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: nj_c(:)            ! Condensed species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)           ! Log of gas species concentrations [kmol-per-kg]
        real(dp), pointer :: h_g(:)             ! Gas enthalpies [unitless]
        real(dp), pointer :: s_g(:)             ! Gas entropies [unitless]
        real(dp), pointer :: u_g(:)             ! Gas energies [unitless]
        real(dp), pointer :: h_c(:)             ! Condensed enthalpies [unitless]
        real(dp), pointer :: s_c(:)             ! Condensed entropies [unitless]
        real(dp), pointer :: cp(:)              ! Species heat capacities [unitless]
        real(dp) :: dh_g_dT(solver%num_gas)     ! d/dT of gas enthalpies [unitless]
        real(dp) :: ds_g_dT(solver%num_gas)     ! d/dT of gas entropies [unitless]
        real(dp) :: dh_c_dT(solver%num_condensed)  ! d/dT of condensed enthalpies [unitless]
        real(dp) :: ds_c_dT(solver%num_condensed)  ! d/dT of condensed entropies [unitless]
        real(dp) :: dcp_g_dT(solver%num_gas)    ! d/dT of gas heat capacities [unitless]
        real(dp) :: dcp_c_dT(solver%num_condensed)  ! d/dT of condensed heat capacities [unitless]
        real(dp) :: s_g_minus(solver%num_gas)   ! s_g - ln(nj) - log(P/n)
        real(dp) :: db0_dw0(solver%num_elements, solver%num_reactants)
        real(dp) :: inv_w_sum, w_sum
        real(dp), pointer :: A_g(:,:), A_c(:,:) ! Gas/condensed stoichiometric matrices
        integer :: i, j, idx_c                  ! Loop counter
        integer :: lnT_idx, lnn_idx             ! Indices for ln(T) and ln(n) derivatives in du/dx
        integer, allocatable :: active_cond_idx(:)
        real(dp) :: ln_n
        real(dp) :: log_p_over_n
        real(dp) :: ln_threshold
        real(dp) :: ln_threshold_nj
        real(dp) :: species_size
        real(dp) :: dlogP_over_n_state1
        real(dp) :: dlogP_over_n_state2
        real(dp), allocatable :: dlogP_over_n_db0(:)
        real(dp), allocatable :: dlogP_over_n_dw0(:)
        real(dp), allocatable :: dT_db0(:)
        real(dp), allocatable :: dn_db0(:)
        real(dp), allocatable :: dln_nj_dstate1(:)
        real(dp), allocatable :: dln_nj_dstate2(:)
        real(dp), allocatable :: dln_nj_eff_dstate1(:)
        real(dp), allocatable :: dln_nj_eff_dstate2(:)
        real(dp), allocatable :: dln_nj_db0(:,:)
        real(dp), allocatable :: dln_nj_dw0(:,:)
        real(dp), allocatable :: dln_nj_eff_dw0(:,:)
        real(dp), allocatable :: dnj_db0(:,:)
        real(dp), allocatable :: dS_sum_dw0(:)
        real(dp), allocatable :: ln_nj_eff(:)
        real(dp), allocatable :: nj_g_eff(:)
        real(dp), allocatable :: dnj_dln_nj(:)
        real(dp), allocatable :: dln_nj_eff_dln_nj(:)
        real(dp), allocatable :: dln_nj_amount_dln_nj(:)
        real(dp) :: sum_h
        real(dp) :: sum_dh_dT
        real(dp) :: sum_h_dnj
        real(dp) :: sum_cp_dnj
        real(dp) :: sum_dcp_dT_term
        real(dp) :: dS_sum_state1
        real(dp) :: dS_sum_state2
        real(dp) :: entropy_sum
        real(dp) :: entropy_dim
        real(dp) :: temp_dT
        real(dp) :: threshold_value
        real(dp) :: threshold_margin
        real(dp) :: fac
        logical :: ion_species
        logical :: const_p, const_t, const_s, const_h, const_u  ! Flags enabling/disabling matrix equations
        type(EqConstraints), pointer :: cons    ! Abbreviation for soln%constraints

        ! Define shorthand
        ng = solver%num_gas
        nc = solver%num_condensed
        ne = solver%num_elements
        nr = solver%num_reactants
        na = count(solution%is_active)
        num_eqn = solution%num_equations(solver)
        cons => solution%constraints
        const_p = cons%is_constant_pressure()
        const_t = cons%is_constant_temperature()
        const_s = cons%is_constant_entropy()
        const_h = cons%is_constant_enthalpy()
        const_u = cons%is_constant_energy()

        ! Associate subarray pointers
        A_g => solver%products%stoich_matrix(:ng,:) ! NOTE: A is transpose of a_ij in RP-1311
        A_c => solver%products%stoich_matrix(ng+1:,:)
        n = solution%n
        nj  => solution%nj
        nj_g => solution%nj(:ng)
        nj_c => solution%nj(ng+1:)
        ln_nj => solution%ln_nj
        h_g => solution%thermo%enthalpy(:ng)
        h_c => solution%thermo%enthalpy(ng+1:)
        s_g => solution%thermo%entropy(:ng)
        s_c => solution%thermo%entropy(ng+1:)
        u_g => solution%thermo%energy(:ng)
        cp => solution%thermo%cp(:)

        ! Get the mixture pressure and temperature
        P = solution%calc_pressure()
        T = solution%T

        ! Compute intermediate derivatives
        do i = 1, ng
            dh_g_dT(i) = solver%products%species(i)%calc_denthalpy_dT(T)/T - h_g(i)/T
            ds_g_dT(i) = solver%products%species(i)%calc_dentropy_dT(T)
            dcp_g_dT(i) = solver%products%species(i)%calc_dcp_dT(T)
        end do
        do i = 1, nc
            dh_c_dT(i) = solver%products%species(ng+i)%calc_denthalpy_dT(T)/T - h_c(i)/T
            ds_c_dT(i) = solver%products%species(ng+i)%calc_dentropy_dT(T)
            dcp_c_dT(i) = solver%products%species(ng+i)%calc_dcp_dT(T)
        end do

        ! Set indices for ln(T) and ln(n) derivatives in du/dx
        if (const_t) then
            lnT_idx = 0
        else
            if (const_p) then
                lnT_idx = ne+na+2
            else
                lnT_idx = ne+na+1
            end if
        end if

        if (const_p) then
            lnn_idx = ne+na+1
        else
            lnn_idx = 0
        end if

        ! Pre-compute db0/dw0 as a common term
        w_sum = sum(solution%w0)
        inv_w_sum = 1.0d0 / w_sum
        do j = 1, nr
            db0_dw0(:, j) = solver%reactants%stoich_matrix(j, :) / &
                solver%reactants%species(j)%molecular_weight
            db0_dw0(:, j) = (db0_dw0(:, j) - cons%b0) * inv_w_sum
        end do

        fac = R / 1.d3
        ln_n = log(n)
        log_p_over_n = log(P/n)
        ln_threshold = ln_n - solver%tsize

        allocate(ln_nj_eff(ng), nj_g_eff(ng), dnj_dln_nj(ng), dln_nj_eff_dln_nj(ng), dln_nj_amount_dln_nj(ng))
        do i = 1, ng
            call compute_nj_effective(ln_nj(i), ln_threshold, solver%smooth_truncation, solver%truncation_width, &
                                      nj_eff=nj_g_eff(i), ln_nj_eff=ln_nj_eff(i), dln_nj_eff_dln_nj=dln_nj_eff_dln_nj(i))
            ion_species = solver%ions .and. solver%active_ions .and. ne > 0 .and. A_g(i, ne) /= 0.0d0
            ln_threshold_nj = gas_amount_ln_threshold(ln_n, solver%tsize, solver%esize, ion_species)
            call compute_nj_effective(ln_nj(i), ln_threshold_nj, solver%smooth_truncation, solver%truncation_width, &
                                      nj_eff=nj_g_eff(i), dln_nj_eff_dln_nj=dln_nj_amount_dln_nj(i))
            if (solver%smooth_truncation) then
                dnj_dln_nj(i) = exp(ln_nj(i)) * dln_nj_amount_dln_nj(i)
            else
                dnj_dln_nj(i) = nj_g_eff(i)
            end if
        end do

        ! Compute gas phase chemical potentials
        mu_g = h_g - s_g + ln_nj_eff + log(P/n)

        if (na > 0) then
            allocate(active_cond_idx(na))
            idx_c = 0
            do i = 1, nc
                if (.not. solution%is_active(i)) cycle
                idx_c = idx_c + 1
                active_cond_idx(idx_c) = i
            end do
        else
            allocate(active_cond_idx(0))
        end if

        ! ---------------------------------------------------------
        ! dT/dx
        ! ---------------------------------------------------------

        ! dT/dstate1:
        if (const_t) then
            self%dT_dstate1 = 1.0d0
        else
            self%dT_dstate1 = T*self%dudx(lnT_idx, 2)  ! dT/dstate1 = T*d(ln(T))/dstate1
        end if

        ! dT/dstate2:
        if (const_t) then
            self%dT_dstate2 = 0.0d0
        else
            self%dT_dstate2 = T*self%dudx(lnT_idx, 1)  ! dT/dstate2 = T*d(ln(T))/dstate2
        end if

        ! dT/dw0:
        if (const_t) then
            self%dT_dw0 = 0.0d0
        else
            self%dT_dw0 = T*matmul(self%dudx(lnT_idx, 3:ne+2), db0_dw0)
        end if

        ! ---------------------------------------------------------
        ! dn/dx
        ! ---------------------------------------------------------

        ! dn/dstate1:
        if (const_p) then
            self%dn_dstate1 = n*self%dudx(lnn_idx, 2)  ! dn/dstate1 = n*d(ln(n))/dstate1
        else
            self%dn_dstate1 = 0.0d0
        end if

        ! dn/dstate2:
        if (const_p) then
            self%dn_dstate2 = n*self%dudx(lnn_idx, 1)  ! dn/dstate2 = n*d(ln(n))/dstate2
        else
            self%dn_dstate2 = 0.0d0
        end if

        ! dn/dw0:
        if (const_p) then
            self%dn_dw0 = n*matmul(self%dudx(lnn_idx, 3:ne+2), db0_dw0)  ! dn/dw0 = n*d(ln(n))/dw0
        else
            self%dn_dw0 = 0.0d0
        end if

        ! ---------------------------------------------------------
        ! dnj/dx
        ! ---------------------------------------------------------

        allocate(dT_db0(ne), dn_db0(ne))
        allocate(dlogP_over_n_db0(ne), dlogP_over_n_dw0(nr))
        allocate(dln_nj_dstate1(ng), dln_nj_dstate2(ng))
        allocate(dln_nj_eff_dstate1(ng), dln_nj_eff_dstate2(ng))
        allocate(dln_nj_db0(ng, ne), dln_nj_dw0(ng, nr))
        allocate(dln_nj_eff_dw0(ng, nr))
        allocate(dnj_db0(ng+na, ne))
        allocate(dS_sum_dw0(nr))

        if (const_t) then
            dT_db0 = 0.0d0
        else
            dT_db0 = T*self%dudx(lnT_idx, 3:ne+2)
        end if

        if (const_p) then
            dn_db0 = n*self%dudx(lnn_idx, 3:ne+2)
        else
            dn_db0 = 0.0d0
        end if

        if (const_p) then
            dlogP_over_n_state1 = -self%dn_dstate1 / n
            dlogP_over_n_state2 = 1.0d0/P - self%dn_dstate2 / n
            dlogP_over_n_db0 = -dn_db0 / n
            dlogP_over_n_dw0 = -self%dn_dw0 / n
        else
            dlogP_over_n_state1 = self%dT_dstate1 / T
            dlogP_over_n_state2 = self%dT_dstate2 / T - 1.0d0/cons%state2
            dlogP_over_n_db0 = dT_db0 / T
            dlogP_over_n_dw0 = self%dT_dw0 / T
        end if

        dln_nj_dstate1 = 0.0d0
        dln_nj_dstate2 = 0.0d0
        dln_nj_eff_dstate1 = 0.0d0
        dln_nj_eff_dstate2 = 0.0d0
        dln_nj_db0 = 0.0d0
        dln_nj_eff_dw0 = 0.0d0
        dnj_db0 = 0.0d0
        self%dnj_dstate1 = 0.0d0
        self%dnj_dstate2 = 0.0d0

        ! ln(nj) = dot(A_g, pi) - h_g + s_g - log(P/n) [+ A_g(i, ne)*pi_e if ions]
        ! pi_e is treated as constant in the total derivatives.
        do i = 1, ng
            temp_dT = ds_g_dT(i) - dh_g_dT(i)

            dln_nj_dstate1(i) = dot_product(A_g(i, :), self%dudx(1:ne, 2)) &
                + temp_dT*self%dT_dstate1 - dlogP_over_n_state1
            dln_nj_dstate2(i) = dot_product(A_g(i, :), self%dudx(1:ne, 1)) &
                + temp_dT*self%dT_dstate2 - dlogP_over_n_state2

            do j = 1, ne
                dln_nj_db0(i, j) = dot_product(A_g(i, :), self%dudx(1:ne, j+2)) &
                    + temp_dT*dT_db0(j) - dlogP_over_n_db0(j)
            end do
        end do

        dln_nj_dw0 = matmul(dln_nj_db0, db0_dw0)
        do i = 1, ng
            dln_nj_eff_dstate1(i) = dln_nj_eff_dln_nj(i) * dln_nj_dstate1(i)
            dln_nj_eff_dstate2(i) = dln_nj_eff_dln_nj(i) * dln_nj_dstate2(i)
            dln_nj_eff_dw0(i, :) = dln_nj_eff_dln_nj(i) * dln_nj_dw0(i, :)
        end do

        do i = 1, ng
            ion_species = solver%ions .and. solver%active_ions .and. ne > 0 .and. A_g(i, ne) /= 0.0d0
            if (ion_species) then
                species_size = solver%esize
            else
                species_size = solver%tsize
            end if
            threshold_value = ln_nj(i) - ln_n + species_size
            threshold_margin = 0.05d0*species_size
            if (solver%smooth_truncation) then
                self%dnj_dstate1(i) = dnj_dln_nj(i)*dln_nj_dstate1(i)
                self%dnj_dstate2(i) = dnj_dln_nj(i)*dln_nj_dstate2(i)
                dnj_db0(i, :) = dnj_dln_nj(i)*dln_nj_db0(i, :)
            else if (threshold_value > 0.0d0) then
                self%dnj_dstate1(i) = nj_g_eff(i)*dln_nj_dstate1(i)
                self%dnj_dstate2(i) = nj_g_eff(i)*dln_nj_dstate2(i)
                dnj_db0(i, :) = nj_g_eff(i)*dln_nj_db0(i, :)
            else
                if (threshold_value >= -threshold_margin) then
                    call log_warning("EqDerivatives_unpack_values: "// &
                        trim(solver%products%species_names(i))// &
                        " not in the active-set so derivatives are 0, but they are close to the threshold")
                end if
            end if
        end do

        do idx_c = 1, na
            i = active_cond_idx(idx_c)
            self%dnj_dstate1(ng+idx_c) = self%dudx(ne+idx_c, 2)
            self%dnj_dstate2(ng+idx_c) = self%dudx(ne+idx_c, 1)
            dnj_db0(ng+idx_c, :) = self%dudx(ne+idx_c, 3:ne+2)
        end do

        self%dnj_dw0 = matmul(dnj_db0, db0_dw0)

        ! For volume-constrained problems, n is the sum of gas species.
        if (.not. const_p) then
            self%dn_dstate1 = sum(self%dnj_dstate1(:ng))
            self%dn_dstate2 = sum(self%dnj_dstate2(:ng))
            do j = 1, nr
                self%dn_dw0(j) = sum(self%dnj_dw0(:ng, j))
            end do
        end if

        ! ---------------------------------------------------------
        ! dH/dx
        ! ---------------------------------------------------------

        sum_h = dot_product(nj, solution%thermo%enthalpy)
        sum_dh_dT = dot_product(nj_g_eff, dh_g_dT)
        if (nc > 0) sum_dh_dT = sum_dh_dT + dot_product(nj_c, dh_c_dT)

        do i = 1, ng
            s_g_minus(i) = s_g(i) - ln_nj_eff(i) - log_p_over_n
        end do

        entropy_sum = dot_product(nj_g_eff, s_g_minus)
        if (nc > 0) entropy_sum = entropy_sum + dot_product(nj_c, s_c)
        entropy_dim = fac*entropy_sum

        ! dH/dstate1:
        if (const_h) then
            self%dH_dstate1 = fac
        else
            sum_h_dnj = dot_product(h_g, self%dnj_dstate1(:ng))
            do idx_c = 1, na
                i = active_cond_idx(idx_c)
                sum_h_dnj = sum_h_dnj + h_c(i)*self%dnj_dstate1(ng+idx_c)
            end do
            self%dH_dstate1 = fac*(T*sum_h_dnj + (sum_h + T*sum_dh_dT)*self%dT_dstate1)
        end if

        ! dH/dstate2:
        if (const_h) then
            self%dH_dstate2 = 0.0d0
        else
            sum_h_dnj = dot_product(h_g, self%dnj_dstate2(:ng))
            do idx_c = 1, na
                i = active_cond_idx(idx_c)
                sum_h_dnj = sum_h_dnj + h_c(i)*self%dnj_dstate2(ng+idx_c)
            end do
            self%dH_dstate2 = fac*(T*sum_h_dnj + (sum_h + T*sum_dh_dT)*self%dT_dstate2)
        end if

        ! dH/dw0:
        if (const_h) then
            self%dH_dw0 = 0.0d0
        else
            do j = 1, nr
                sum_h_dnj = dot_product(h_g, self%dnj_dw0(:ng, j))
                do idx_c = 1, na
                    i = active_cond_idx(idx_c)
                    sum_h_dnj = sum_h_dnj + h_c(i)*self%dnj_dw0(ng+idx_c, j)
                end do
                self%dH_dw0(j) = fac*(T*sum_h_dnj + (sum_h + T*sum_dh_dT)*self%dT_dw0(j))
            end do
        end if

        ! ---------------------------------------------------------
        ! dU/dx
        ! ---------------------------------------------------------

        ! dU/dstate1:
        self%dU_dstate1 = self%dH_dstate1 - fac*(n*self%dT_dstate1 + T*self%dn_dstate1)

        ! dU/dstate2:
        self%dU_dstate2 = self%dH_dstate2 - fac*(n*self%dT_dstate2 + T*self%dn_dstate2)

        ! dU/dw0:
        do j = 1, nr
            self%dU_dw0(j) = self%dH_dw0(j) - fac*(n*self%dT_dw0(j) + T*self%dn_dw0(j))
        end do

        ! ---------------------------------------------------------
        ! dS/dx
        ! ---------------------------------------------------------

        ! dS/dstate1:
        if (const_s) then
            ! state1 is entropy/R; convert to dimensional entropy derivative
            self%dS_dstate1 = fac
        else
            dS_sum_state1 = 0.0d0
            do i = 1, ng
                dS_sum_state1 = dS_sum_state1 + self%dnj_dstate1(i)*s_g_minus(i)
                dS_sum_state1 = dS_sum_state1 + nj_g_eff(i)* &
                    (ds_g_dT(i)*self%dT_dstate1 - dln_nj_eff_dstate1(i))
            end do
            ! Add the log(P/n) derivative terms once (not once per species)
            dS_sum_state1 = dS_sum_state1 - n*dlogP_over_n_state1
            do idx_c = 1, na
                i = active_cond_idx(idx_c)
                dS_sum_state1 = dS_sum_state1 + self%dnj_dstate1(ng+idx_c)*s_c(i)
                dS_sum_state1 = dS_sum_state1 + nj_c(i)*ds_c_dT(i)*self%dT_dstate1
            end do
            self%dS_dstate1 = fac*dS_sum_state1
        end if

        ! dS/dstate2:
        if (const_s) then
            self%dS_dstate2 = 0.0d0
        else
            dS_sum_state2 = 0.0d0
            do i = 1, ng
                dS_sum_state2 = dS_sum_state2 + self%dnj_dstate2(i)*s_g_minus(i)
                dS_sum_state2 = dS_sum_state2 + nj_g_eff(i)* &
                    (ds_g_dT(i)*self%dT_dstate2 - dln_nj_eff_dstate2(i))
            end do
            ! Add the log(P/n) derivative terms once (not once per species)
            dS_sum_state2 = dS_sum_state2 - n*dlogP_over_n_state2
            do idx_c = 1, na
                i = active_cond_idx(idx_c)
                dS_sum_state2 = dS_sum_state2 + self%dnj_dstate2(ng+idx_c)*s_c(i)
                dS_sum_state2 = dS_sum_state2 + nj_c(i)*ds_c_dT(i)*self%dT_dstate2
            end do
            self%dS_dstate2 = fac*dS_sum_state2
        end if

        ! dS/dw0:
        if (const_s) then
            self%dS_dw0 = 0.0d0
        else
            do j = 1, nr
                dS_sum_dw0(j) = 0.0d0
                do i = 1, ng
                    dS_sum_dw0(j) = dS_sum_dw0(j) + self%dnj_dw0(i, j)*s_g_minus(i)
                    dS_sum_dw0(j) = dS_sum_dw0(j) + nj_g_eff(i)* &
                        (ds_g_dT(i)*self%dT_dw0(j) - dln_nj_eff_dw0(i, j))
                end do
                ! Add the log(P/n) derivative terms once (not once per species)
                dS_sum_dw0(j) = dS_sum_dw0(j) - n*dlogP_over_n_dw0(j)
                do idx_c = 1, na
                    i = active_cond_idx(idx_c)
                    dS_sum_dw0(j) = dS_sum_dw0(j) + self%dnj_dw0(ng+idx_c, j)*s_c(i)
                    dS_sum_dw0(j) = dS_sum_dw0(j) + nj_c(i)*ds_c_dT(i)*self%dT_dw0(j)
                end do
                self%dS_dw0(j) = fac*dS_sum_dw0(j)
            end do
        end if

        ! ---------------------------------------------------------
        ! dG/dx
        ! ---------------------------------------------------------

        ! dG/dstate1:
        self%dG_dstate1 = self%dH_dstate1 - entropy_dim*self%dT_dstate1 - T*self%dS_dstate1

        ! dG/dstate2:
        self%dG_dstate2 = self%dH_dstate2 - entropy_dim*self%dT_dstate2 - T*self%dS_dstate2

        ! dG/dw0:
        do j = 1, nr
            self%dG_dw0(j) = self%dH_dw0(j) - entropy_dim*self%dT_dw0(j) - T*self%dS_dw0(j)
        end do

        ! ---------------------------------------------------------
        ! dCp_fr/dx
        ! ---------------------------------------------------------

        ! dCp_fr/dstate1:
        ! dCp_fr/dx = sum_j(dnj/dx * Cp_j) + sum_j(nj * dCp_j/dT * dT/dx)
        sum_cp_dnj = dot_product(cp(:ng), self%dnj_dstate1(:ng))
        do idx_c = 1, na
            i = active_cond_idx(idx_c)
            sum_cp_dnj = sum_cp_dnj + cp(ng+i)*self%dnj_dstate1(ng+idx_c)
        end do
        ! Add temperature derivative terms: sum_j(nj * dCp_j/dT) * dT/dstate1
        sum_dcp_dT_term = 0.0d0
        do i = 1, ng
            sum_dcp_dT_term = sum_dcp_dT_term + nj_g_eff(i) * dcp_g_dT(i)
        end do
        do idx_c = 1, na
            i = active_cond_idx(idx_c)
            sum_dcp_dT_term = sum_dcp_dT_term + nj_c(i) * dcp_c_dT(i)
        end do
        self%dCp_fr_dstate1 = fac*(sum_cp_dnj + sum_dcp_dT_term*self%dT_dstate1)

        ! dCp_fr/dstate2:
        sum_cp_dnj = dot_product(cp(:ng), self%dnj_dstate2(:ng))
        do idx_c = 1, na
            i = active_cond_idx(idx_c)
            sum_cp_dnj = sum_cp_dnj + cp(ng+i)*self%dnj_dstate2(ng+idx_c)
        end do
        ! Add temperature derivative terms: sum_j(nj * dCp_j/dT) * dT/dstate2
        sum_dcp_dT_term = 0.0d0
        do i = 1, ng
            sum_dcp_dT_term = sum_dcp_dT_term + nj_g_eff(i) * dcp_g_dT(i)
        end do
        do idx_c = 1, na
            i = active_cond_idx(idx_c)
            sum_dcp_dT_term = sum_dcp_dT_term + nj_c(i) * dcp_c_dT(i)
        end do
        self%dCp_fr_dstate2 = fac*(sum_cp_dnj + sum_dcp_dT_term*self%dT_dstate2)

        ! dCp_fr/dw0:
        do j = 1, nr
            sum_cp_dnj = dot_product(cp(:ng), self%dnj_dw0(:ng, j))
            do idx_c = 1, na
                i = active_cond_idx(idx_c)
                sum_cp_dnj = sum_cp_dnj + cp(ng+i)*self%dnj_dw0(ng+idx_c, j)
            end do
            ! Add temperature derivative terms: sum_j(nj * dCp_j/dT) * dT/dw0
            sum_dcp_dT_term = 0.0d0
            do i = 1, ng
                sum_dcp_dT_term = sum_dcp_dT_term + nj_g_eff(i) * dcp_g_dT(i)
            end do
            do idx_c = 1, na
                i = active_cond_idx(idx_c)
                sum_dcp_dT_term = sum_dcp_dT_term + nj_c(i) * dcp_c_dT(i)
            end do
            self%dCp_fr_dw0(j) = fac*(sum_cp_dnj + sum_dcp_dT_term*self%dT_dw0(j))
        end do

    end subroutine

    subroutine EqDerivatives_compute_derivatives(self, solver, solution, check_closure_defect)

        ! Arguments
        class(EqDerivatives), intent(inout) :: self
        class(EqSolver), intent(in) :: solver
        class(EqSolution), intent(inout) :: solution
        logical, intent(in), optional :: check_closure_defect

        ! Locals
        integer :: i, ierr
        real(dp) :: G(self%m, self%m+1)

        call log_debug("Starting compute_derivatives")

        ! Check if the solution is converged; raise a warning if not
        if (.not. solution%converged) then
            call log_warning("Computing derivatives for a non-converged solution.")
        end if

        ! Compute the Jacobian and Rx matrices
        call self%assemble_jacobian(solver, solution)
        call self%assemble_Rx(solver, solution)

        ! Compute the derivatives: du/dx = -J^-1 * Rx
        do i = 1, self%n
            ierr = 0
            G(:, :self%m) = self%J
            G(:, self%m+1) = -self%Rx(:, i)
            call gauss(G, ierr)
            if (ierr /= 0) then
                call log_warning("Singular matrix in total derivatives")
                self%dudx(:, i) = 0.0d0
                cycle
            end if
            self%dudx(:, i) = G(:, self%m+1)
        end do

        ! Check output for closure defect: delta = J*du/dx + Rx = 0
        if (present(check_closure_defect)) then
            if (check_closure_defect) then
                 call log_debug("Checking closure defect: ||J*du/dx + Rx|| should be close to 0.")
                 call self%check_closure_defect()
            end if
        end if

    end subroutine

    subroutine EqDerivatives_compute_fd(self, solver, solution, h, verbose, central)

        ! Arguments
        class(EqDerivatives), intent(inout) :: self
        class(EqSolver), intent(in), target :: solver
        class(EqSolution), intent(inout), target :: solution
        real(dp), intent(in) :: h
        logical, intent(in), optional :: verbose
        logical, intent(in), optional :: central

        ! Locals
        type(EqSolution), target :: pert_soln, pert_soln_minus
        real(dp), allocatable :: base_nj(:)
        real(dp), allocatable :: pert_nj(:), pert_nj_minus(:)
        real(dp), allocatable :: w0(:)
        integer, allocatable :: active_cond_idx(:)
        integer :: ng, na, nr, ns
        integer :: i, j, idx_c
        real(dp) :: base_T, base_n
        real(dp) :: base_H, base_U, base_G, base_S, base_Cp_fr
        real(dp) :: abs_err, rel_err
        logical :: verbose_, central_
        character(2) :: ctype
        real(dp) :: state1, state2
        real(dp) :: h_state1, h_state2, h_w
        integer :: ne
        integer :: lnT_idx, lnn_idx
        integer :: idx_max_state1, idx_max_state2, idx_max_dw0, j_max_dw0
        real(dp) :: max_err_dw0
        real(dp) :: base_P, base_logP_over_n, pert_logP_over_n, pert_logP_over_n_minus
        real(dp) :: dlogP_over_n_state1_fd, dlogP_over_n_state2_fd
        real(dp) :: dlogP_over_n_dw0_fd, dlogP_over_n_dw0_fd_max
        real(dp) :: dlogP_over_n_state1_an, dlogP_over_n_state2_an
        real(dp) :: dh_g_dT(solver%num_gas), ds_g_dT(solver%num_gas)
        real(dp), allocatable :: dln_nj_state1_fd(:)
        real(dp), allocatable :: dln_nj_state2_fd(:)
        real(dp), allocatable :: dln_nj_dw0_fd_max(:)
        real(dp), pointer :: A_g(:,:)
        logical :: const_p, const_t
        type(EqConstraints), pointer :: cons

        ! NOTE: EqDerivatives_compute_derivatives and EqDerivatives_unpack_values should be called first.

        verbose_ = .true.
        if (present(verbose)) verbose_ = verbose

        central_ = .false.
        if (present(central)) central_ = central

        ng = solver%num_gas
        ne = solver%num_elements
        na = count(solution%is_active)
        nr = solver%num_reactants
        ns = ng + na
        cons => solution%constraints
        const_p = cons%is_constant_pressure()
        const_t = cons%is_constant_temperature()

        if (const_t) then
            lnT_idx = 0
        else
            if (const_p) then
                lnT_idx = ne+na+2
            else
                lnT_idx = ne+na+1
            end if
        end if

        if (const_p) then
            lnn_idx = ne+na+1
        else
            lnn_idx = 0
        end if

        allocate(base_nj(ns), pert_nj(ns))
        if (central_) allocate(pert_nj_minus(ns))
        if (na > 0) then
            allocate(active_cond_idx(na))
            idx_c = 0
            do i = 1, solver%num_condensed
                if (.not. solution%is_active(i)) cycle
                idx_c = idx_c + 1
                active_cond_idx(idx_c) = i
            end do
        else
            allocate(active_cond_idx(0))
        end if

        base_nj(:ng) = solution%nj(:ng)
        do idx_c = 1, na
            i = active_cond_idx(idx_c)
            base_nj(ng+idx_c) = solution%nj(ng+i)
        end do

        base_T = solution%T
        base_n = solution%n
        base_P = solution%calc_pressure()
        base_logP_over_n = log(base_P/base_n)
        base_H = solution%enthalpy
        base_U = solution%energy
        base_G = solution%gibbs_energy
        base_S = solution%entropy
        base_Cp_fr = solution%cp_fr
        if (verbose_) then
            allocate(dln_nj_state1_fd(ng), dln_nj_state2_fd(ng), dln_nj_dw0_fd_max(ng))
            dln_nj_state1_fd = 0.0d0
            dln_nj_state2_fd = 0.0d0
            dln_nj_dw0_fd_max = 0.0d0
            idx_max_state1 = 1
            idx_max_state2 = 1
            idx_max_dw0 = 1
            j_max_dw0 = 1
            max_err_dw0 = -1.0d0
            dlogP_over_n_state1_fd = 0.0d0
            dlogP_over_n_state2_fd = 0.0d0
            dlogP_over_n_dw0_fd_max = 0.0d0
            A_g => solver%products%stoich_matrix(:ng,:)
            do i = 1, ng
                dh_g_dT(i) = solver%products%species(i)%calc_denthalpy_dT(base_T)/base_T - &
                    solution%thermo%enthalpy(i)/base_T
                ds_g_dT(i) = solver%products%species(i)%calc_dentropy_dT(base_T)
            end do
        end if

        ctype = solution%constraints%type
        state1 = solution%constraints%state1
        state2 = solution%constraints%state2
        h_state1 = h * max(1.0d0, abs(state1))
        h_state2 = h * max(1.0d0, abs(state2))

        allocate(w0(nr))
        w0 = solution%w0

        ! state1 perturbation
        pert_soln = solution
        pert_soln%cp_fr = 0.0d0
        call solver%solve(pert_soln, ctype, state1 + h_state1, state2, w0)
        pert_nj(:ng) = pert_soln%nj(:ng)
        do idx_c = 1, na
            i = active_cond_idx(idx_c)
            pert_nj(ng+idx_c) = pert_soln%nj(ng+i)
        end do

        if (central_) then
            ! Central difference: compute backward perturbation
            pert_soln_minus = solution
            pert_soln_minus%cp_fr = 0.0d0
            call solver%solve(pert_soln_minus, ctype, state1 - h_state1, state2, w0)
            pert_nj_minus(:ng) = pert_soln_minus%nj(:ng)
            do idx_c = 1, na
                i = active_cond_idx(idx_c)
                pert_nj_minus(ng+idx_c) = pert_soln_minus%nj(ng+i)
            end do

            self%dT_dstate1_fd = (pert_soln%T - pert_soln_minus%T) / (2.0d0*h_state1)
            self%dn_dstate1_fd = (pert_soln%n - pert_soln_minus%n) / (2.0d0*h_state1)
            self%dnj_dstate1_fd = (pert_nj - pert_nj_minus) / (2.0d0*h_state1)
            self%dH_dstate1_fd = (pert_soln%enthalpy - pert_soln_minus%enthalpy) / (2.0d0*h_state1)
            self%dU_dstate1_fd = (pert_soln%energy - pert_soln_minus%energy) / (2.0d0*h_state1)
            self%dG_dstate1_fd = (pert_soln%gibbs_energy - pert_soln_minus%gibbs_energy) / (2.0d0*h_state1)
            self%dS_dstate1_fd = (pert_soln%entropy - pert_soln_minus%entropy) / (2.0d0*h_state1)
            self%dCp_fr_dstate1_fd = (pert_soln%cp_fr - pert_soln_minus%cp_fr) / (2.0d0*h_state1)
            if (verbose_) then
                pert_logP_over_n = log(pert_soln%calc_pressure()/pert_soln%n)
                pert_logP_over_n_minus = log(pert_soln_minus%calc_pressure()/pert_soln_minus%n)
                dlogP_over_n_state1_fd = (pert_logP_over_n - pert_logP_over_n_minus) / (2.0d0*h_state1)
                do i = 1, ng
                    if (pert_nj(i) > 0.0d0 .and. pert_nj_minus(i) > 0.0d0) then
                        dln_nj_state1_fd(i) = (log(pert_nj(i)) - log(pert_nj_minus(i))) / (2.0d0*h_state1)
                    else
                        dln_nj_state1_fd(i) = 0.0d0
                    end if
                end do
                idx_max_state1 = maxloc(abs(self%dnj_dstate1_fd(:ng) - self%dnj_dstate1(:ng)), dim=1)
            end if
        else
            ! Forward difference
            self%dT_dstate1_fd = (pert_soln%T - base_T) / h_state1
            self%dn_dstate1_fd = (pert_soln%n - base_n) / h_state1
            self%dnj_dstate1_fd = (pert_nj - base_nj) / h_state1
            self%dH_dstate1_fd = (pert_soln%enthalpy - base_H) / h_state1
            self%dU_dstate1_fd = (pert_soln%energy - base_U) / h_state1
            self%dG_dstate1_fd = (pert_soln%gibbs_energy - base_G) / h_state1
            self%dS_dstate1_fd = (pert_soln%entropy - base_S) / h_state1
            self%dCp_fr_dstate1_fd = (pert_soln%cp_fr - base_Cp_fr) / h_state1
            if (verbose_) then
                pert_logP_over_n = log(pert_soln%calc_pressure()/pert_soln%n)
                dlogP_over_n_state1_fd = (pert_logP_over_n - base_logP_over_n) / h_state1
                do i = 1, ng
                    if (base_nj(i) > 0.0d0 .and. pert_nj(i) > 0.0d0) then
                        dln_nj_state1_fd(i) = (log(pert_nj(i)) - log(base_nj(i))) / h_state1
                    else
                        dln_nj_state1_fd(i) = 0.0d0
                    end if
                end do
                idx_max_state1 = maxloc(abs(self%dnj_dstate1_fd(:ng) - self%dnj_dstate1(:ng)), dim=1)
            end if
        end if

        ! state2 perturbation
        pert_soln = solution
        pert_soln%cp_fr = 0.0d0
        call solver%solve(pert_soln, ctype, state1, state2 + h_state2, w0)
        pert_nj(:ng) = pert_soln%nj(:ng)
        do idx_c = 1, na
            i = active_cond_idx(idx_c)
            pert_nj(ng+idx_c) = pert_soln%nj(ng+i)
        end do

        if (central_) then
            ! Central difference: compute backward perturbation
            pert_soln_minus = solution
            pert_soln_minus%cp_fr = 0.0d0
            call solver%solve(pert_soln_minus, ctype, state1, state2 - h_state2, w0)
            pert_nj_minus(:ng) = pert_soln_minus%nj(:ng)
            do idx_c = 1, na
                i = active_cond_idx(idx_c)
                pert_nj_minus(ng+idx_c) = pert_soln_minus%nj(ng+i)
            end do

            self%dT_dstate2_fd = (pert_soln%T - pert_soln_minus%T) / (2.0d0*h_state2)
            self%dn_dstate2_fd = (pert_soln%n - pert_soln_minus%n) / (2.0d0*h_state2)
            self%dnj_dstate2_fd = (pert_nj - pert_nj_minus) / (2.0d0*h_state2)
            self%dH_dstate2_fd = (pert_soln%enthalpy - pert_soln_minus%enthalpy) / (2.0d0*h_state2)
            self%dU_dstate2_fd = (pert_soln%energy - pert_soln_minus%energy) / (2.0d0*h_state2)
            self%dG_dstate2_fd = (pert_soln%gibbs_energy - pert_soln_minus%gibbs_energy) / (2.0d0*h_state2)
            self%dS_dstate2_fd = (pert_soln%entropy - pert_soln_minus%entropy) / (2.0d0*h_state2)
            self%dCp_fr_dstate2_fd = (pert_soln%cp_fr - pert_soln_minus%cp_fr) / (2.0d0*h_state2)
            if (verbose_) then
                pert_logP_over_n = log(pert_soln%calc_pressure()/pert_soln%n)
                pert_logP_over_n_minus = log(pert_soln_minus%calc_pressure()/pert_soln_minus%n)
                dlogP_over_n_state2_fd = (pert_logP_over_n - pert_logP_over_n_minus) / (2.0d0*h_state2)
                do i = 1, ng
                    if (pert_nj(i) > 0.0d0 .and. pert_nj_minus(i) > 0.0d0) then
                        dln_nj_state2_fd(i) = (log(pert_nj(i)) - log(pert_nj_minus(i))) / (2.0d0*h_state2)
                    else
                        dln_nj_state2_fd(i) = 0.0d0
                    end if
                end do
                idx_max_state2 = maxloc(abs(self%dnj_dstate2_fd(:ng) - self%dnj_dstate2(:ng)), dim=1)
            end if
        else
            ! Forward difference
            self%dT_dstate2_fd = (pert_soln%T - base_T) / h_state2
            self%dn_dstate2_fd = (pert_soln%n - base_n) / h_state2
            self%dnj_dstate2_fd = (pert_nj - base_nj) / h_state2
            self%dH_dstate2_fd = (pert_soln%enthalpy - base_H) / h_state2
            self%dU_dstate2_fd = (pert_soln%energy - base_U) / h_state2
            self%dG_dstate2_fd = (pert_soln%gibbs_energy - base_G) / h_state2
            self%dS_dstate2_fd = (pert_soln%entropy - base_S) / h_state2
            self%dCp_fr_dstate2_fd = (pert_soln%cp_fr - base_Cp_fr) / h_state2
            if (verbose_) then
                pert_logP_over_n = log(pert_soln%calc_pressure()/pert_soln%n)
                dlogP_over_n_state2_fd = (pert_logP_over_n - base_logP_over_n) / h_state2
                do i = 1, ng
                    if (base_nj(i) > 0.0d0 .and. pert_nj(i) > 0.0d0) then
                        dln_nj_state2_fd(i) = (log(pert_nj(i)) - log(base_nj(i))) / h_state2
                    else
                        dln_nj_state2_fd(i) = 0.0d0
                    end if
                end do
                idx_max_state2 = maxloc(abs(self%dnj_dstate2_fd(:ng) - self%dnj_dstate2(:ng)), dim=1)
            end if
        end if

        ! weight perturbations
        do j = 1, nr
            h_w = h * max(1.0d0, abs(w0(j)))

            if (central_) then
                ! Central difference: compute forward perturbation
                w0(j) = w0(j) + h_w
                pert_soln = solution
                pert_soln%cp_fr = 0.0d0
                call solver%solve(pert_soln, ctype, state1, state2, w0)
                pert_nj(:ng) = pert_soln%nj(:ng)
                do idx_c = 1, na
                    i = active_cond_idx(idx_c)
                    pert_nj(ng+idx_c) = pert_soln%nj(ng+i)
                end do
                w0(j) = w0(j) - h_w

                ! Compute backward perturbation
                w0(j) = w0(j) - h_w
                pert_soln_minus = solution
                pert_soln_minus%cp_fr = 0.0d0
                call solver%solve(pert_soln_minus, ctype, state1, state2, w0)
                pert_nj_minus(:ng) = pert_soln_minus%nj(:ng)
                do idx_c = 1, na
                    i = active_cond_idx(idx_c)
                    pert_nj_minus(ng+idx_c) = pert_soln_minus%nj(ng+i)
                end do
                w0(j) = w0(j) + h_w

                self%dT_dw0_fd(j) = (pert_soln%T - pert_soln_minus%T) / (2.0d0*h_w)
                self%dn_dw0_fd(j) = (pert_soln%n - pert_soln_minus%n) / (2.0d0*h_w)
                self%dnj_dw0_fd(:, j) = (pert_nj - pert_nj_minus) / (2.0d0*h_w)
                self%dH_dw0_fd(j) = (pert_soln%enthalpy - pert_soln_minus%enthalpy) / (2.0d0*h_w)
                self%dU_dw0_fd(j) = (pert_soln%energy - pert_soln_minus%energy) / (2.0d0*h_w)
                self%dG_dw0_fd(j) = (pert_soln%gibbs_energy - pert_soln_minus%gibbs_energy) / (2.0d0*h_w)
                self%dS_dw0_fd(j) = (pert_soln%entropy - pert_soln_minus%entropy) / (2.0d0*h_w)
                self%dCp_fr_dw0_fd(j) = (pert_soln%cp_fr - pert_soln_minus%cp_fr) / (2.0d0*h_w)
                if (verbose_) then
                    pert_logP_over_n = log(pert_soln%calc_pressure()/pert_soln%n)
                    pert_logP_over_n_minus = log(pert_soln_minus%calc_pressure()/pert_soln_minus%n)
                    dlogP_over_n_dw0_fd = (pert_logP_over_n - pert_logP_over_n_minus) / (2.0d0*h_w)
                    abs_err = maxval(abs(self%dnj_dw0_fd(:ng, j) - self%dnj_dw0(:ng, j)))
                    if (abs_err > max_err_dw0) then
                        max_err_dw0 = abs_err
                        idx_max_dw0 = maxloc(abs(self%dnj_dw0_fd(:ng, j) - self%dnj_dw0(:ng, j)), dim=1)
                        j_max_dw0 = j
                        do i = 1, ng
                            if (pert_nj(i) > 0.0d0 .and. pert_nj_minus(i) > 0.0d0) then
                                dln_nj_dw0_fd_max(i) = (log(pert_nj(i)) - log(pert_nj_minus(i))) / (2.0d0*h_w)
                            else
                                dln_nj_dw0_fd_max(i) = 0.0d0
                            end if
                        end do
                        dlogP_over_n_dw0_fd_max = dlogP_over_n_dw0_fd
                    end if
                end if
            else
                ! Forward difference
                w0(j) = w0(j) + h_w
                pert_soln = solution
                pert_soln%cp_fr = 0.0d0
                call solver%solve(pert_soln, ctype, state1, state2, w0)
                pert_nj(:ng) = pert_soln%nj(:ng)
                do idx_c = 1, na
                    i = active_cond_idx(idx_c)
                    pert_nj(ng+idx_c) = pert_soln%nj(ng+i)
                end do

                self%dT_dw0_fd(j) = (pert_soln%T - base_T) / h_w
                self%dn_dw0_fd(j) = (pert_soln%n - base_n) / h_w
                self%dnj_dw0_fd(:, j) = (pert_nj - base_nj) / h_w
                self%dH_dw0_fd(j) = (pert_soln%enthalpy - base_H) / h_w
                self%dU_dw0_fd(j) = (pert_soln%energy - base_U) / h_w
                self%dG_dw0_fd(j) = (pert_soln%gibbs_energy - base_G) / h_w
                self%dS_dw0_fd(j) = (pert_soln%entropy - base_S) / h_w
                self%dCp_fr_dw0_fd(j) = (pert_soln%cp_fr - base_Cp_fr) / h_w
                if (verbose_) then
                    pert_logP_over_n = log(pert_soln%calc_pressure()/pert_soln%n)
                    dlogP_over_n_dw0_fd = (pert_logP_over_n - base_logP_over_n) / h_w
                    abs_err = maxval(abs(self%dnj_dw0_fd(:ng, j) - self%dnj_dw0(:ng, j)))
                    if (abs_err > max_err_dw0) then
                        max_err_dw0 = abs_err
                        idx_max_dw0 = maxloc(abs(self%dnj_dw0_fd(:ng, j) - self%dnj_dw0(:ng, j)), dim=1)
                        j_max_dw0 = j
                        do i = 1, ng
                            if (base_nj(i) > 0.0d0 .and. pert_nj(i) > 0.0d0) then
                                dln_nj_dw0_fd_max(i) = (log(pert_nj(i)) - log(base_nj(i))) / h_w
                            else
                                dln_nj_dw0_fd_max(i) = 0.0d0
                            end if
                        end do
                        dlogP_over_n_dw0_fd_max = dlogP_over_n_dw0_fd
                    end if
                end if

                w0(j) = w0(j) - h_w
            end if
        end do

        if (verbose_) then
            write(*,*) "EqDerivatives_compute_fd: FD vs analytic derivatives"

            abs_err = abs(self%dT_dstate1_fd - self%dT_dstate1)
            rel_err = abs_err / max(abs(self%dT_dstate1), 1.0d-30)
            write(*,*) "dT/dstate1: abs=", abs_err, " rel=", rel_err

            abs_err = abs(self%dT_dstate2_fd - self%dT_dstate2)
            rel_err = abs_err / max(abs(self%dT_dstate2), 1.0d-30)
            write(*,*) "dT/dstate2: abs=", abs_err, " rel=", rel_err

            if (nr > 0) then
                abs_err = maxval(abs(self%dT_dw0_fd - self%dT_dw0))
                rel_err = maxval(abs(self%dT_dw0_fd - self%dT_dw0) / max(abs(self%dT_dw0), 1.0d-30))
                write(*,*) "dT/dw0 (max): abs=", abs_err, " rel=", rel_err
            end if

            abs_err = abs(self%dn_dstate1_fd - self%dn_dstate1)
            rel_err = abs_err / max(abs(self%dn_dstate1), 1.0d-30)
            write(*,*) "dn/dstate1: abs=", abs_err, " rel=", rel_err

            abs_err = abs(self%dn_dstate2_fd - self%dn_dstate2)
            rel_err = abs_err / max(abs(self%dn_dstate2), 1.0d-30)
            write(*,*) "dn/dstate2: abs=", abs_err, " rel=", rel_err

            if (nr > 0) then
                abs_err = maxval(abs(self%dn_dw0_fd - self%dn_dw0))
                rel_err = maxval(abs(self%dn_dw0_fd - self%dn_dw0) / max(abs(self%dn_dw0), 1.0d-30))
                write(*,*) "dn/dw0 (max): abs=", abs_err, " rel=", rel_err
            end if

            do i = 1, ng
                if (solution%nj(i) > 1.0d-10) then
                    abs_err = abs(self%dnj_dstate1_fd(i) - self%dnj_dstate1(i))
                    rel_err = abs_err / max(abs(self%dnj_dstate1(i)), 1.0d-30)
                    write(*,*) "dnj/dstate1 (", solver%products%species_names(i), "): abs=", abs_err, " rel=", rel_err

                    abs_err = abs(self%dnj_dstate2_fd(i) - self%dnj_dstate2(i))
                    rel_err = abs_err / max(abs(self%dnj_dstate2(i)), 1.0d-30)
                    write(*,*) "dnj/dstate2 (", solver%products%species_names(i), "): abs=", abs_err, " rel=", rel_err

                    if (nr > 0) then
                        abs_err = maxval(abs(self%dnj_dw0_fd(i, :) - self%dnj_dw0(i, :)))
                        rel_err = maxval(abs(self%dnj_dw0_fd(i, :) - self%dnj_dw0(i, :)) / max(abs(self%dnj_dw0(i, :)), 1.0d-30))
                        write(*,*) "dnj/dw0 (", solver%products%species_names(i), ") (max): abs=", abs_err, " rel=", rel_err
                    end if
                end if
            end do

            abs_err = maxval(abs(self%dnj_dstate1_fd - self%dnj_dstate1))
            rel_err = maxval(abs(self%dnj_dstate1_fd - self%dnj_dstate1) / max(abs(self%dnj_dstate1), 1.0d-30))
            write(*,*) "dnj/dstate1 (max): abs=", abs_err, " rel=", rel_err

            abs_err = maxval(abs(self%dnj_dstate2_fd - self%dnj_dstate2))
            rel_err = maxval(abs(self%dnj_dstate2_fd - self%dnj_dstate2) / max(abs(self%dnj_dstate2), 1.0d-30))
            write(*,*) "dnj/dstate2 (max): abs=", abs_err, " rel=", rel_err

            if (nr > 0) then
                abs_err = maxval(abs(self%dnj_dw0_fd - self%dnj_dw0))
                rel_err = maxval(abs(self%dnj_dw0_fd - self%dnj_dw0) / max(abs(self%dnj_dw0), 1.0d-30))
                write(*,*) "dnj/dw0 (max): abs=", abs_err, " rel=", rel_err
            end if

            if (const_p) then
                dlogP_over_n_state1_an = -self%dn_dstate1 / base_n
                dlogP_over_n_state2_an = 1.0d0/base_P - self%dn_dstate2 / base_n
            else
                dlogP_over_n_state1_an = self%dT_dstate1 / base_T
                dlogP_over_n_state2_an = self%dT_dstate2 / base_T - 1.0d0/cons%state2
            end if

            abs_err = abs(self%dH_dstate1_fd - self%dH_dstate1)
            rel_err = abs_err / max(abs(self%dH_dstate1), 1.0d-30)
            write(*,*) "dH/dstate1: abs=", abs_err, " rel=", rel_err

            abs_err = abs(self%dH_dstate2_fd - self%dH_dstate2)
            rel_err = abs_err / max(abs(self%dH_dstate2), 1.0d-30)
            write(*,*) "dH/dstate2: abs=", abs_err, " rel=", rel_err

            if (nr > 0) then
                abs_err = maxval(abs(self%dH_dw0_fd - self%dH_dw0))
                rel_err = maxval(abs(self%dH_dw0_fd - self%dH_dw0) / max(abs(self%dH_dw0), 1.0d-30))
                write(*,*) "dH/dw0 (max): abs=", abs_err, " rel=", rel_err
            end if

            abs_err = abs(self%dU_dstate1_fd - self%dU_dstate1)
            rel_err = abs_err / max(abs(self%dU_dstate1), 1.0d-30)
            write(*,*) "dU/dstate1: abs=", abs_err, " rel=", rel_err

            abs_err = abs(self%dU_dstate2_fd - self%dU_dstate2)
            rel_err = abs_err / max(abs(self%dU_dstate2), 1.0d-30)
            write(*,*) "dU/dstate2: abs=", abs_err, " rel=", rel_err

            if (nr > 0) then
                abs_err = maxval(abs(self%dU_dw0_fd - self%dU_dw0))
                rel_err = maxval(abs(self%dU_dw0_fd - self%dU_dw0) / max(abs(self%dU_dw0), 1.0d-30))
                write(*,*) "dU/dw0 (max): abs=", abs_err, " rel=", rel_err
            end if

            abs_err = abs(self%dG_dstate1_fd - self%dG_dstate1)
            rel_err = abs_err / max(abs(self%dG_dstate1), 1.0d-30)
            write(*,*) "dG/dstate1: abs=", abs_err, " rel=", rel_err

            abs_err = abs(self%dG_dstate2_fd - self%dG_dstate2)
            rel_err = abs_err / max(abs(self%dG_dstate2), 1.0d-30)
            write(*,*) "dG/dstate2: abs=", abs_err, " rel=", rel_err

            if (nr > 0) then
                abs_err = maxval(abs(self%dG_dw0_fd - self%dG_dw0))
                rel_err = maxval(abs(self%dG_dw0_fd - self%dG_dw0) / max(abs(self%dG_dw0), 1.0d-30))
                write(*,*) "dG/dw0 (max): abs=", abs_err, " rel=", rel_err
            end if

            abs_err = abs(self%dS_dstate1_fd - self%dS_dstate1)
            rel_err = abs_err / max(abs(self%dS_dstate1), 1.0d-30)
            write(*,*) "dS/dstate1: abs=", abs_err, " rel=", rel_err

            abs_err = abs(self%dS_dstate2_fd - self%dS_dstate2)
            rel_err = abs_err / max(abs(self%dS_dstate2), 1.0d-30)
            write(*,*) "dS/dstate2: abs=", abs_err, " rel=", rel_err

            if (nr > 0) then
                abs_err = maxval(abs(self%dS_dw0_fd - self%dS_dw0))
                rel_err = maxval(abs(self%dS_dw0_fd - self%dS_dw0) / max(abs(self%dS_dw0), 1.0d-30))
                write(*,*) "dS/dw0 (max): abs=", abs_err, " rel=", rel_err
            end if

            abs_err = abs(self%dCp_fr_dstate1_fd - self%dCp_fr_dstate1)
            rel_err = abs_err / max(abs(self%dCp_fr_dstate1), 1.0d-30)
            write(*,*) "dCp_fr/dstate1: abs=", abs_err, " rel=", rel_err

            abs_err = abs(self%dCp_fr_dstate2_fd - self%dCp_fr_dstate2)
            rel_err = abs_err / max(abs(self%dCp_fr_dstate2), 1.0d-30)
            write(*,*) "dCp_fr/dstate2: abs=", abs_err, " rel=", rel_err

            if (nr > 0) then
                abs_err = maxval(abs(self%dCp_fr_dw0_fd - self%dCp_fr_dw0))
                rel_err = maxval(abs(self%dCp_fr_dw0_fd - self%dCp_fr_dw0) / max(abs(self%dCp_fr_dw0), 1.0d-30))
                write(*,*) "dCp_fr/dw0 (max): abs=", abs_err, " rel=", rel_err
            end if
        end if

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
        allocate(self%b0(num_elements))
        self%b0 = empty_dp
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
        allocate(self%nj_seed(solver%num_products), source=0.0d0)
        allocate(self%ln_nj_seed(solver%num_gas), source=0.0d0)
        allocate(self%G(solver%max_equations, solver%max_equations+1), source=empty_dp)
        allocate(self%is_active(solver%num_condensed), source=.false.)
        allocate(self%w0(solver%num_reactants), source=0.0d0)
        allocate(self%active_rank(solver%num_condensed), source=0)
        allocate(self%is_active_seed(solver%num_condensed), source=.false.)
        allocate(self%active_rank_seed(solver%num_condensed), source=0)
        allocate(self%transport_component_idx(solver%num_elements), source=0)
        allocate(self%transport_basis_matrix(solver%num_elements, solver%num_gas), source=0.0d0)
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

        call EqSolution_save_seed(self)

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
                ! NOTE(smooth_truncation): Initialization keeps the legacy hard floor for
                ! non-positive seeds to preserve stable and backward-compatible starts.
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
        integer  :: ne_full                  ! Total number of elements (including electron)
        real(dp) :: n                        ! Total moles of mixture
        real(dp) :: ln_n                     ! Log total moles
        real(dp) :: P                        ! Pressure of mixture (bar)
        real(dp), pointer :: nj(:)           ! Total/gas species concentrations [kmol-per-kg]
        real(dp), pointer :: ln_nj(:)        ! Log of gas species concentrations [kmol-per-kg]
        real(dp) :: ln_nj_eff(solver%num_gas)
        real(dp) :: ln_threshold
        real(dp), pointer :: s_g(:), s_c(:)  ! Gas/condensed entropies [unitless]
        integer :: i
        real(dp) :: nj_eff_tmp
        logical :: ion_species

        ! Shorthand
        ng = solver%num_gas
        ne_full = solver%num_elements
        n = self%n
        ln_n = log(n)
        nj => self%nj
        ln_nj => self%ln_nj
        s_g => self%thermo%entropy(:ng)
        s_c => self%thermo%entropy(ng+1:)
        P = self%calc_pressure()

        do i = 1, ng
            ion_species = solver%ions .and. solver%active_ions .and. ne_full > 0 .and. &
                          solver%products%stoich_matrix(i, ne_full) /= 0.0d0
            ln_threshold = gas_amount_ln_threshold(ln_n, solver%tsize, solver%esize, ion_species)
            call compute_nj_effective(ln_nj(i), ln_threshold, solver%smooth_truncation, solver%truncation_width, &
                                      nj_eff=nj_eff_tmp, ln_nj_eff=ln_nj_eff(i))
        end do

        S = dot_product(nj(:ng), s_g-ln_nj_eff-log(P/n)) + dot_product(nj(ng+1:), s_c)

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

        allocate(active_idx(0))

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

        allocate(active_idx(0))

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

        allocate(active_idx(0))

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
        integer, allocatable :: selected_transport_pure_idx(:)  ! Selected indices into transport pure species list
        integer, allocatable :: selected_local_idx(:)           ! Local transport indices (into idx_list) for selected pure species
        integer, allocatable :: transport_to_local(:)           ! Reverse map from transport pure index to local transport index
        integer, allocatable :: bin_idx(:)    ! List of indices (into transport order) to compute transport properties for
        integer :: bin_count                  ! Total number of binary pairs to consider
        integer :: ng                         ! Number of gas species
        integer :: ne                         ! Number of elements
        integer :: np, nb                     ! Number of pure, binary species
        integer :: nm                         ! Number of of gaseous species for thermal transport calculations
        integer :: nr                         ! Number of chemical reactions
        integer :: i, j, k, k1, k2, ii, m     ! Index counters
        integer :: idx(1), idx1(1), idx2(1)   ! Temporary findloc index
        integer :: local_idx1, local_idx2
        integer :: best_idx
        integer :: max_elem_idx               ! Number of elements (minus electron, if applicaple)
        integer :: ncomp
        integer :: nseed
        integer, allocatable :: comp_local_idx(:), comp_basis_row(:)
        logical, allocatable :: is_component(:)
        real(dp) :: cfit_val                  ! Value computed by curve-fit function
        real(dp) :: best_nj
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
        logical, allocatable :: selected_species(:)

        if (present(frozen_shock)) then
            continue  ! Placeholder for future frozen-shock transport handling.
        end if

        ! Define shorthand
        np = eq_solver%transport_db%num_pure
        nb = eq_solver%transport_db%num_binary
        ng = eq_solver%num_gas
        ne = eq_solver%num_elements
        A => eq_solver%products%stoich_matrix

        ! Allocate
        allocate(psi(ng, ng), phi(ng, ng), eta(ng, ng), cond(ng), &
                 idx_list(max_tr), selected_transport_pure_idx(np), selected_local_idx(np), &
                 transport_to_local(np), bin_idx(nb), selected_species(ng), &
                 cp(ng), xs(ng), G(ng, ng), rtpd(ng, ng), &
                 xsij(ng, ng), delh(ng), alpha(ng, ng), &
                 stcf(ng, ng), stcoef(ng), tmp(max_tr), gmat(ng, ng), &
                 stx(ng), stxij(ng, ng), comp_local_idx(max_tr), &
                 comp_basis_row(max_tr), is_component(max_tr))

        ! cond can be used without being initialized, and uninitialized elements
        ! could be used later if all of the species aren't found in the transport database
        cond = 0.0d0

        ! Build the list of relevant mixture species.
        ! Legacy CEA seeds transport from a basis-like set (Jcm/Lsave), then expands by abundance.
        ! We approximate that behavior by seeding one dominant carrier per active element, then
        ! adding monoatomic species and finally expanding by threshold.
        nm = 0
        total = 0.0d0
        selected_species = .false.
        wtmol = 1.0/sum(eq_soln%nj)
        nj_cutoff = 1.d-11/wtmol
        test_tot = 0.999999999d0/wtmol
        max_elem_idx = eq_solver%num_elements
        if (eq_solver%ions) max_elem_idx = max_elem_idx - 1
        nj_el = 0.0d0

        do i = 1, ng
            ! TODO(smooth_truncation): transport species screening currently uses hard-zero semantics.
            ! Revisit whether smooth mode should use a practical-zero cutoff instead.
            if (eq_soln%nj(i) <= 0.0d0) then
                if (eq_soln%ln_nj(i) - log(eq_soln%n) + eq_solver%xsize > 0.0d0) then
                    eq_soln%nj(i) = exp(eq_soln%ln_nj(i))
                end if
            end if
        end do

        ! Seed transport species with the cached equilibrium component basis.
        nseed = 0
        if (eq_soln%transport_basis_rows > 0) then
            do m = 1, min(eq_soln%transport_basis_rows, size(eq_soln%transport_component_idx))
                i = eq_soln%transport_component_idx(m)
                if (i < 1 .or. i > ng) cycle
                if (selected_species(i)) cycle
                if (nm >= max_tr) exit
                nm = nm + 1
                idx_list(nm) = i
                selected_species(i) = .true.
                total = total + eq_soln%nj(i)
                if (eq_solver%products%species(i)%molecular_weight < 1.0d0) nj_el = eq_soln%nj(i)
                eq_soln%nj(i) = -eq_soln%nj(i)
                nseed = nseed + 1
            end do
        end if

        ! Fallback when no cached basis is available.
        if (nseed == 0) then
            do m = 1, max_elem_idx
                best_idx = 0
                best_nj = -1.0d0
                do i = 1, ng
                    if (A(i, m) > tol .and. eq_soln%nj(i) > best_nj) then
                        best_idx = i
                        best_nj = eq_soln%nj(i)
                    end if
                end do
                if (best_idx > 0 .and. .not. selected_species(best_idx) .and. nm < max_tr) then
                    nm = nm + 1
                    idx_list(nm) = best_idx
                    selected_species(best_idx) = .true.
                    total = total + eq_soln%nj(best_idx)
                    if (eq_solver%products%species(best_idx)%molecular_weight < 1.0d0) nj_el = eq_soln%nj(best_idx)
                    eq_soln%nj(best_idx) = -eq_soln%nj(best_idx)
                end if
            end do

            do i = 1, ng
                if ((sum(abs(A(i, :)))-1.0d0) < tol .and. .not. selected_species(i)) then
                    if (nm >= max_tr) then
                        call log_info("Reached maximum number of allowable transport species.")
                        exit
                    end if
                    nm = nm + 1
                    idx_list(nm) = i
                    selected_species(i) = .true.
                    total = total + eq_soln%nj(i)
                    if (eq_solver%products%species(i)%molecular_weight < 1.0d0) nj_el = eq_soln%nj(i)
                    eq_soln%nj(i) = -eq_soln%nj(i)
                end if
            end do
        end if
        test_nj = 1.0d0/(ng*wtmol)

        ! Add the remaining species that meet the minimum size threshold
        do i = 1, ng
            if (total <= test_tot .and. nm < max_tr) then
                test_nj = test_nj / 10.0d0
                do j = 1, ng
                    if (eq_soln%nj(j) >= test_nj .and. .not. selected_species(j)) then
                        if (nm >= max_tr) then
                            call log_info("Reached maximum number of allowable transport species.")
                            exit
                        else
                            total = total + eq_soln%nj(j)
                            nm = nm + 1
                            idx_list(nm) = j
                            selected_species(j) = .true.
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

        if (nm <= 0 .or. total <= 0.0d0) return

        ! Align electron concentration with the finalized transport species set.
        do i = 1, nm
            if (eq_solver%products%species(idx_list(i))%molecular_weight < 1.0d0) then
                nj_el = eq_soln%nj(idx_list(i))
            end if
        end do

        ! Build aligned pure-species mappings.
        j = 0
        transport_to_local = 0
        do i = 1, nm
            idx = findloc(eq_solver%transport_db%pure_species, eq_solver%products%species_names(idx_list(i)))
            if (idx(1) > 0) then
                j = j + 1
                selected_transport_pure_idx(j) = idx(1)
                selected_local_idx(j) = i
                transport_to_local(idx(1)) = i
            else
                call log_info('compute_transport_properties: Species '//eq_solver%products%species_names(idx_list(i))//&
                              ' not found in transport database.')
            end if
        end do
        selected_transport_pure_idx = selected_transport_pure_idx(:j)
        selected_local_idx = selected_local_idx(:j)
        np = j

        ! Remove any binary pairs with negligible concentrations
        bin_count = 0
        do i = 1, nb
            idx1 = findloc(eq_solver%transport_db%pure_species, eq_solver%transport_db%binary_species(i,1))
            idx2 = findloc(eq_solver%transport_db%pure_species, eq_solver%transport_db%binary_species(i,2))
            if (idx1(1) > 0 .and. idx2(1) > 0) then
                local_idx1 = transport_to_local(idx1(1))
                local_idx2 = transport_to_local(idx2(1))
                if (local_idx1 > 0 .and. local_idx2 > 0) then
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
                local_idx1 = transport_to_local(idx1(1))
                local_idx2 = transport_to_local(idx2(1))
                if (local_idx1 > 0 .and. local_idx2 > 0) then
                    cfit_val = exp(eq_solver%transport_db%binary_transport(bin_idx(i))%calc_eta(eq_soln%T))
                    eta(local_idx1, local_idx2) = cfit_val
                    eta(local_idx2, local_idx1) = cfit_val
                end if
            else
                call log_info('compute_transport_properties: Binary species'//eq_solver%transport_db%binary_species(bin_idx(i),1)//&
                ' or '//eq_solver%transport_db%binary_species(bin_idx(i),2)//'not found in transport database')
            end if
        end do

        ! Add the diagonal terms
        do i = 1, np
            eta(selected_local_idx(i), selected_local_idx(i)) = &
                exp(eq_solver%transport_db%pure_transport(selected_transport_pure_idx(i))%calc_eta(eq_soln%T))
        end do

        ! Build the conductivity array
        cond = cond(:nm)
        do i = 1, np
            cond(selected_local_idx(i)) = &
                exp(eq_solver%transport_db%pure_transport(selected_transport_pure_idx(i))%calc_lambda(eq_soln%T))
        end do

        ! Build the stoichiometric matrix for the chemical reactions.
        ! Prefer the cached equilibrium component basis to match legacy TRANIN/TRANP behavior.
        alpha = alpha(:, :nm)
        alpha = 0.0d0
        is_component = .false.
        comp_local_idx = 0
        comp_basis_row = 0
        ncomp = 0

        if (eq_soln%transport_basis_rows > 0 .and. &
            allocated(eq_soln%transport_component_idx) .and. &
            allocated(eq_soln%transport_basis_matrix)) then

            do m = 1, min(eq_soln%transport_basis_rows, size(eq_soln%transport_component_idx))
                j = eq_soln%transport_component_idx(m)
                if (j < 1 .or. j > ng) cycle
                do i = 1, nm
                    if (idx_list(i) /= j) cycle
                    if (.not. is_component(i)) then
                        ncomp = ncomp + 1
                        if (ncomp > max_tr) exit
                        comp_local_idx(ncomp) = i
                        comp_basis_row(ncomp) = m
                        is_component(i) = .true.
                    end if
                    exit
                end do
            end do
        end if

        ! Legacy TRANIN uses Lsave components, which can include the electron row
        ! when ions are active. Append the electron species as an extra component
        ! if it is present in the selected transport set but not already included.
        if (eq_solver%ions) then
            do i = 1, nm
                if (eq_solver%products%species(idx_list(i))%molecular_weight < 1.0d0) then
                    if (.not. is_component(i) .and. ncomp < max_tr) then
                        ncomp = ncomp + 1
                        comp_local_idx(ncomp) = i
                        comp_basis_row(ncomp) = 0
                        is_component(i) = .true.
                    end if
                    exit
                end if
            end do
        end if

        if (ncomp > 0 .and. ncomp < nm) then
            nr = 0
            do i = 1, nm
                if (is_component(i)) cycle
                nr = nr + 1
                alpha(nr, i) = -1.0d0
                j = idx_list(i)
                do k = 1, ncomp
                    m = comp_local_idx(k)
                    if (comp_basis_row(k) > 0) then
                        alpha(nr, m) = eq_soln%transport_basis_matrix(comp_basis_row(k), j)
                    else
                        alpha(nr, m) = A(j, ne)
                    end if
                end do
            end do
        else
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
        end if
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
            if (.not. (eq_solver%ions .and. abs(abs(A(k, ne))-1.0d0) < tol .and. &
                abs(eta(i, i)) < tol)) then
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
            do j = i, nm
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

        if (nr > 0) then
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
        else
            cpreac = 0.0d0
            reacon = 0.0d0
        end if

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
