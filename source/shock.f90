module cea_shock
    !! CEA Shock Solver module

    use cea_param, only: dp, empty_dp, R=>gas_constant
    use cea_mixture, only: Mixture, MixtureThermo
    use cea_transport, only: TransportDB
    use cea_equilibrium, only: EqSolution, EqSolver, EqPartials, compute_transport_properties
    use fb_findloc, only: findloc
    use fb_utils
    implicit none

    integer, parameter :: SHOCK_SOLVE_STATUS_NOT_CONVERGED = 0
    integer, parameter :: SHOCK_SOLVE_STATUS_CONVERGED = 1
    integer, parameter :: SHOCK_SOLVE_STATUS_LAST_VALID = 2

    type :: ShockSolver
        !! Shock solver object

        type(EqSolver) :: eq_solver
            !! Equilibrium solver object

    contains

        procedure :: solve => ShockSolver_solve
        procedure :: solve_incident => ShockSolver_solve_incident
        procedure :: solve_incident_frozen => ShockSolver_solve_incident_frozen
        procedure :: solve_reflected => ShockSolver_solve_reflected
        procedure :: solve_reflected_frozen => ShockSolver_solve_reflected_frozen
        procedure :: update_solution => ShockSolver_update_solution

    end type
    interface ShockSolver
        module procedure :: ShockSolver_init
    end interface

    type :: ShockSolution
        !! Shock solution object

        integer :: num_pts
            !! Number of evaluation points

        ! Thermo solution at each evaluation point
        type(EqSolution), allocatable :: eq_soln(:)
            !! Equilibrium solution object at each evaluation point
        type(EqPartials), allocatable :: eq_partials(:)
            !! Partial derivatives of the equilibrium solution at each point

        ! States
        real(dp), allocatable :: pressure(:)
            !! Pressure [bar]
        real(dp), allocatable :: mach(:)
            !! u/v_sonic at stations [1, 2, 5], in the same shock frames as u.
        real(dp), allocatable :: u(:)
            !! Gas speed [m/s]: stations 1 and 2 in the incident-shock frame;
            !! station 5 in the reflected-shock frame. Station 5 is at rest
            !! in the wall frame, so u(3) also equals the reflected wave-speed
            !! magnitude. These are speeds, not signed wall-frame velocities.
        real(dp), allocatable :: v_sonic(:)
            !! Reported sound speed [m/s], sqrt(R*T*gamma_s/M).
            !! Initial/incident-frozen: frozen Cp; equilibrium: equilibrium
            !! isentropic exponent. Reflected-frozen retains the legacy
            !! incident-state frozen Cp basis at the reflected temperature.

        ! Solution variables
        real(dp) :: rho12, rho52 = 0.0d0
            !! Ratios of density across the incident and reflected shocks
        real(dp) :: p21, p52 = 0.0d0
            !! Pressure ratios across the incident and reflected shocks
        real(dp) :: t21, t52 = 0.0d0
            !! Temperature ratios across the incident and reflected shocks
        real(dp) :: M21, M52 = 0.0d0
            !! Molecular-weight ratios M2/M1 and M5/M2 (legacy CEA2 output), not Mach ratios.
        real(dp) :: v2 = 0.0d0
            !! Station-2 gas speed in the wall/unshocked-gas frame: u(1)-u(2) [m/s].
        real(dp) :: u5_p_v2 = 0.0d0
            !! Upstream gas speed in the reflected-shock frame: u(3)+v2 [m/s].
            !! Also the reflected wave speed relative to station-2 gas, not the wall.

        ! Convergence variables
        logical :: converged = .false.
            !! Convergence flag
        integer :: solve_status = SHOCK_SOLVE_STATUS_NOT_CONVERGED
            !! Internal solve outcome used by bindings to distinguish strict
            !! non-convergence from a retained last-valid incident state.

    end type
    interface ShockSolution
        module procedure :: ShockSolution_init
    end interface

contains

    !-----------------------------------------------------------------------
    ! ShockSolver
    !-----------------------------------------------------------------------
    function ShockSolver_init(products, reactants, trace, ions, all_transport, insert, &
            smooth_truncation, truncation_width) result(self)
        type(ShockSolver) :: self
        type(Mixture), intent(in) :: products
        type(Mixture), intent(in), optional :: reactants
        real(dp), intent(in), optional :: trace
        logical, intent(in), optional :: ions
        type(TransportDB), intent(in), optional :: all_transport
        character(*), intent(in), optional :: insert(:)  ! List of condensed species to insert
        logical, intent(in), optional :: smooth_truncation
        real(dp), intent(in), optional :: truncation_width

        ! Initialize the equilibrium solver
        if (present(trace)) then
            self%eq_solver = EqSolver(products, reactants, trace=trace, ions=ions, all_transport=all_transport, &
                                      insert=insert, smooth_truncation=smooth_truncation, &
                                      truncation_width=truncation_width)
        else
            ! Use a default trace of 5e-9 for shock problems.
            self%eq_solver = EqSolver(products, reactants, trace=5.d-9, ions=ions, all_transport=all_transport, &
                                      insert=insert, smooth_truncation=smooth_truncation, &
                                      truncation_width=truncation_width)
        end if
    end function

    subroutine ShockSolver_update_solution(self, soln, X1, X2, p21, t21, iter)
        ! Compute the damped update factor, apply the solution update, and check convergence

        ! Arguments
        class(ShockSolver), intent(in) :: self
        type(ShockSolution), intent(inout) :: soln
        real(dp), intent(inout) :: X1, X2           ! Solution update variables
        real(dp), intent(inout) :: p21, t21         ! Pressure/temperature ratio across shock
        integer, intent(in) :: iter                 ! Iteration number

        ! Locals
        real(dp) :: cormax, ax  ! Max correction factor; correction coefficient

        ax = max(abs(X1), abs(X2))
        if (ax >= 0.00005d0) then  ! Shock threshold for update damping
            cormax = .40546511d0
            if (iter > 4) then
                cormax = .22314355d0
            end if
            if (iter > 12) then
                cormax = .09531018d0
            end if
            if (iter > 20) then
                cormax = .04879016d0
            end if
            ax = ax/cormax
            if (ax > 1.d0) then
                X1 = X1/ax
                X2 = X2/ax
            end if

            ! Apply the update to the solution
            p21 = EXP(log(p21) + X1)
            t21 = EXP(log(t21) + X2)
        else  ! Converged
            soln%converged = .true.
            soln%solve_status = SHOCK_SOLVE_STATUS_CONVERGED
        end if

    end subroutine

    logical function ShockSolver_state_valid(soln, idx) result(valid)
        type(ShockSolution), intent(in) :: soln
        integer, intent(in) :: idx

        valid = .false.
        if (idx < 1 .or. idx > soln%num_pts) return

        ! Continue with the last nonzero equilibrium state after
        ! singular-recovery exits, and let the shock layer decide whether to
        ! finalize or reject the state.
        valid = soln%eq_soln(idx)%T > 0.0d0 .and. &
                soln%eq_soln(idx)%n > 0.0d0
    end function

    subroutine ShockSolver_fail_state(soln, idx, message)
        type(ShockSolution), intent(inout) :: soln
        integer, intent(in) :: idx
        character(*), intent(in) :: message

        if (len_trim(message) > 0) call log_warning(trim(message))
        call log_warning("Shock calculation stopped after the last valid point.")

        soln%converged = .false.
        soln%solve_status = SHOCK_SOLVE_STATUS_NOT_CONVERGED
        soln%pressure(idx) = 0.0d0
        soln%mach(idx) = 0.0d0
        soln%u(idx) = 0.0d0
        soln%v_sonic(idx) = 0.0d0
        soln%eq_soln(idx)%T = 0.0d0
        soln%eq_soln(idx)%n = 0.0d0
        soln%eq_soln(idx)%converged = .false.
        soln%eq_soln(idx)%density = 0.0d0
        soln%eq_soln(idx)%enthalpy = 0.0d0
        soln%eq_soln(idx)%energy = 0.0d0
        soln%eq_soln(idx)%gibbs_energy = 0.0d0
        soln%eq_soln(idx)%entropy = 0.0d0
        soln%eq_soln(idx)%cp_fr = 0.0d0
        soln%eq_soln(idx)%cp_eq = 0.0d0
        soln%eq_soln(idx)%cp_eq_transport = 0.0d0
        soln%eq_soln(idx)%cv_eq = 0.0d0
        soln%eq_soln(idx)%gamma_s = 0.0d0
        soln%eq_soln(idx)%viscosity = 0.0d0
        soln%eq_soln(idx)%conductivity_fr = 0.0d0
        soln%eq_soln(idx)%conductivity_eq = 0.0d0
        soln%eq_soln(idx)%Pr_fr = 0.0d0
        soln%eq_soln(idx)%Pr_eq = 0.0d0
        soln%eq_partials(idx)%dlnV_dlnP = 0.0d0
        soln%eq_partials(idx)%dlnV_dlnT = 0.0d0
        soln%eq_partials(idx)%gamma_s = 0.0d0

        select case (idx)
            case (2)
                soln%rho12 = 0.0d0
                soln%p21 = 0.0d0
                soln%t21 = 0.0d0
                soln%M21 = 0.0d0
                soln%v2 = 0.0d0
                if (soln%num_pts >= 3) then
                    soln%pressure(3) = 0.0d0
                    soln%mach(3) = 0.0d0
                    soln%u(3) = 0.0d0
                    soln%v_sonic(3) = 0.0d0
                    soln%eq_soln(3)%T = 0.0d0
                    soln%eq_soln(3)%n = 0.0d0
                    soln%eq_soln(3)%converged = .false.
                    soln%eq_soln(3)%density = 0.0d0
                    soln%eq_soln(3)%enthalpy = 0.0d0
                    soln%eq_soln(3)%energy = 0.0d0
                    soln%eq_soln(3)%gibbs_energy = 0.0d0
                    soln%eq_soln(3)%entropy = 0.0d0
                    soln%eq_soln(3)%cp_fr = 0.0d0
                    soln%eq_soln(3)%cp_eq = 0.0d0
                    soln%eq_soln(3)%cp_eq_transport = 0.0d0
                    soln%eq_soln(3)%cv_eq = 0.0d0
                    soln%eq_soln(3)%gamma_s = 0.0d0
                    soln%eq_soln(3)%viscosity = 0.0d0
                    soln%eq_soln(3)%conductivity_fr = 0.0d0
                    soln%eq_soln(3)%conductivity_eq = 0.0d0
                    soln%eq_soln(3)%Pr_fr = 0.0d0
                    soln%eq_soln(3)%Pr_eq = 0.0d0
                    soln%eq_partials(3)%dlnV_dlnP = 0.0d0
                    soln%eq_partials(3)%dlnV_dlnT = 0.0d0
                    soln%eq_partials(3)%gamma_s = 0.0d0
                    soln%rho52 = 0.0d0
                    soln%p52 = 0.0d0
                    soln%t52 = 0.0d0
                    soln%M52 = 0.0d0
                    soln%u5_p_v2 = 0.0d0
                end if
            case (3)
                soln%rho52 = 0.0d0
                soln%p52 = 0.0d0
                soln%t52 = 0.0d0
                soln%M52 = 0.0d0
                soln%u5_p_v2 = 0.0d0
        end select
    end subroutine

    subroutine ShockSolver_update_transport(self, eq_soln, frozen_shock, update_basis)
        class(ShockSolver), intent(in) :: self
        type(EqSolution), intent(inout) :: eq_soln
        logical, intent(in), optional :: frozen_shock
        logical, intent(in), optional :: update_basis
        logical :: do_update_basis
        logical :: frozen_shock_

        if (.not. self%eq_solver%transport) return

        do_update_basis = .true.
        if (present(update_basis)) do_update_basis = update_basis

        frozen_shock_ = .false.
        if (present(frozen_shock)) frozen_shock_ = frozen_shock

        if (frozen_shock_) then
            call compute_transport_properties(self%eq_solver, eq_soln, frozen_shock=.true.)
        else
            if (do_update_basis) then
                call self%eq_solver%update_transport_basis(eq_soln)
            end if
            call compute_transport_properties(self%eq_solver, eq_soln)
        end if
    end subroutine

    subroutine ShockSolver_apply_atomic_transport_basis(self, eq_soln)
        ! Force transport basis rows to elemental gas species.
        class(ShockSolver), intent(in) :: self
        type(EqSolution), intent(inout) :: eq_soln

        integer :: ne, nn, ng
        integer :: i, lc, fallback
        real(dp), parameter :: tol = 1.0d-8

        ng = self%eq_solver%num_gas
        ne = self%eq_solver%num_active_elements()
        nn = ne
        if (self%eq_solver%ions .and. self%eq_solver%active_ions) nn = max(1, ne-1)
        if (nn <= 0 .or. ng <= 0) return

        eq_soln%transport_basis_rows = nn
        eq_soln%transport_component_idx = 0
        eq_soln%transport_basis_matrix = 0.0d0

        do lc = 1, nn
            fallback = 0
            do i = 1, ng
                if (self%eq_solver%products%stoich_matrix(i, lc) > tol .and. fallback == 0) fallback = i
                if (abs(self%eq_solver%products%stoich_matrix(i, lc) - 1.0d0) < tol .and. &
                    abs(sum(abs(self%eq_solver%products%stoich_matrix(i, :ne))) - 1.0d0) < tol) then
                    eq_soln%transport_component_idx(lc) = i
                    exit
                end if
            end do
            if (eq_soln%transport_component_idx(lc) == 0) eq_soln%transport_component_idx(lc) = fallback
            eq_soln%transport_basis_matrix(lc, :ng) = self%eq_solver%products%stoich_matrix(:ng, lc)
        end do
    end subroutine

    subroutine ShockSolver_finalize_equilibrium_state(self, soln, idx, update_transport)
        class(ShockSolver), intent(in) :: self
        type(ShockSolution), intent(inout) :: soln
        integer, intent(in) :: idx
        logical, intent(in), optional :: update_transport
        logical :: update_transport_

        if (.not. ShockSolver_state_valid(soln, idx)) return

        update_transport_ = .true.
        if (present(update_transport)) update_transport_ = update_transport

        soln%eq_partials(idx) = EqPartials(self%eq_solver%num_elements, count(soln%eq_soln(idx)%is_active))
        call soln%eq_partials(idx)%compute_partials(self%eq_solver, soln%eq_soln(idx))
        if (self%eq_solver%transport .and. update_transport_) call ShockSolver_update_transport(self, soln%eq_soln(idx))
        call self%eq_solver%post_process(soln%eq_soln(idx), .true.)
    end subroutine

    subroutine ShockSolver_finalize_reflected_frozen_state(self, soln, idx)
        class(ShockSolver), intent(in) :: self
        type(ShockSolution), intent(inout) :: soln
        integer, intent(in) :: idx

        integer :: i, ng
        real(dp) :: pressure, cp_dimless, wmx, pmn
        real(dp) :: entropy_sum

        if (.not. ShockSolver_state_valid(soln, idx)) return

        ng = self%eq_solver%num_gas
        pressure = soln%pressure(idx)
        wmx = soln%eq_soln(2)%M

        soln%eq_soln(idx)%mole_fractions = soln%eq_soln(idx)%nj / sum(soln%eq_soln(idx)%nj)
        soln%eq_soln(idx)%mass_fractions = soln%eq_soln(idx)%nj * self%eq_solver%products%species%molecular_weight / &
            sum(soln%eq_soln(idx)%nj * self%eq_solver%products%species%molecular_weight)

        soln%eq_soln(idx)%pressure = pressure
        soln%eq_soln(idx)%volume = soln%eq_soln(idx)%calc_volume()
        soln%eq_soln(idx)%density = 1.0d0/soln%eq_soln(idx)%volume

        soln%eq_soln(idx)%enthalpy = dot_product(soln%eq_soln(idx)%nj(:ng), soln%eq_soln(idx)%thermo%enthalpy(:ng)) * &
                                     R * soln%eq_soln(idx)%T / 1.d3
        soln%eq_soln(idx)%energy = soln%eq_soln(idx)%enthalpy - soln%eq_soln(idx)%n*soln%eq_soln(idx)%T*R/1.d3

        entropy_sum = 0.0d0
        do i = 1, size(soln%eq_soln(idx)%nj)
            if (soln%eq_soln(idx)%nj(i) <= 0.0d0) cycle
            pmn = pressure*wmx*soln%eq_soln(idx)%nj(i)
            entropy_sum = entropy_sum + soln%eq_soln(idx)%nj(i)*(soln%eq_soln(idx)%thermo%entropy(i) - log(pmn))
        end do
        soln%eq_soln(idx)%entropy = entropy_sum * R / 1.d3
        soln%eq_soln(idx)%gibbs_energy = soln%eq_soln(idx)%enthalpy - soln%eq_soln(idx)%T*soln%eq_soln(idx)%entropy

        ! For reflected-frozen header thermo, use the carried incident-state
        ! frozen-gas Cp basis with incident molecular weight.
        cp_dimless = dot_product(soln%eq_soln(2)%nj(:ng), soln%eq_soln(2)%thermo%cp(:ng))
        soln%eq_soln(idx)%cp_eq = cp_dimless * R / 1.d3
        soln%eq_soln(idx)%cv_eq = soln%eq_soln(idx)%cp_eq - soln%eq_soln(idx)%n*R/1.d3
        soln%eq_partials(idx)%gamma_s = cp_dimless/(cp_dimless - 1.0d0/wmx)
        soln%eq_soln(idx)%gamma_s = soln%eq_partials(idx)%gamma_s
        soln%v_sonic(idx) = sqrt(R*soln%eq_soln(idx)%T*soln%eq_soln(idx)%gamma_s/wmx)
        ! This legacy sound-speed basis overwrites the iteration value.
        soln%mach(idx) = soln%u(idx)/soln%v_sonic(idx)

        soln%eq_soln(idx)%M = wmx
        soln%eq_soln(idx)%MW = wmx*(1.0d0 - sum(soln%eq_soln(idx)%mole_fractions(ng+1:)))
    end subroutine

    subroutine ShockSolver_solve_incident(self, soln, weights, T0, P0)
        ! Solve the incident shock conditions

        ! Arguments
        class(ShockSolver), intent(in) :: self
        class(ShockSolution), intent(inout) :: soln
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: T0                     ! Initial reactant temperature [K]
        real(dp), intent(in) :: P0                     ! Initial reactant pressure [bar]

        ! Locals
        integer :: idx  ! Solution index for the incident conditions
        real(dp) :: mach1, u1             ! Incident Mach and velocity
        integer :: i                      ! Loop index
        real(dp) :: gamma1                ! Ratio of specific heats at initial condition
        real(dp) :: cp                    ! Mixture heat capacity
        real(dp) :: wm, wm_k              ! Mixture molecular weight (initial, k-th iteration)
        real(dp) :: h_init, h0            ! Mixture enthalpy (initial, <all other points>)
        real(dp) :: T2                    ! Temperature after incident, reflected shocks [K]
        real(dp) :: p21, t21, ttmax       ! Pressure/temperature ratio across the incident shock
        real(dp) :: G(2, 3)               ! Solution matrix
        real(dp) :: X(3)                  ! Solution vector
        real(dp) :: dlnV_dlnP, dlnV_dlnT  ! Partial derivatives
        real(dp) :: mu12rt, rho12         ! Ratios of chemical potential and density across the incident shock
        real(dp) :: tmp                   ! Intermediate variabls
        real(dp), allocatable :: nj_g(:)  ! Total/gas species concentrations [kmol-per-kg]
        integer, parameter :: max_iter = 60
        real(dp), parameter :: T_gas_max = 20000.d0  ! Max gas temperature in the thermo database [K]
        type(EqSolution) :: last_valid_eq
        type(EqPartials) :: last_valid_partials
        real(dp) :: last_rho12, last_p21, last_t21, last_T2, last_wm_k
        logical :: have_last_valid

        ! Initialize
        idx = 2
        G = 0.d0
        soln%converged = .false.
        soln%solve_status = SHOCK_SOLVE_STATUS_NOT_CONVERGED
        have_last_valid = .false.
        u1 = soln%u(1)
        mach1 = soln%mach(1)
        ! Seed incident-equilibrium solve from the unshocked reactant composition
        soln%eq_soln(idx) = EqSolution(self%eq_solver, T_init=T0, nj_init=soln%eq_soln(1)%nj)

        ! Compute the molecular weight of the initial mixture
        wm = sum(weights)

        ! Compute properties of the initial mixture
        cp = self%eq_solver%reactants%calc_frozen_cp(weights, T0)/R
        gamma1 = cp/(cp - 1.0/wm)
        h_init = self%eq_solver%reactants%calc_enthalpy(weights, T0)/R

        ! Compute the solution for the incident shock
        p21 = (2.0d0*gamma1*mach1**2.d0 - gamma1 + 1.d0)/(gamma1 + 1.d0)
        h0 = h_init + u1**2/(2.d0*R)
        mu12rt = wm*u1**2/(R*T0)

        soln%pressure(idx) = p21*P0
        call self%eq_solver%solve(soln%eq_soln(idx), "hp", h0, soln%pressure(idx), weights, partials=soln%eq_partials(idx))
        if (.not. ShockSolver_state_valid(soln, idx)) then
            call ShockSolver_fail_state(soln, idx, &
                                        "ShockSolver_solve_incident: incident equilibrium initialization failed.")
            return
        end if

        t21 = soln%eq_soln(idx)%T/T0
        ttmax = 1.05*T_gas_max/T0
        t21 = min(t21, ttmax)

        do i = 1, max_iter
            ! Update the pressure
            soln%pressure(idx) = p21*P0
            T2 = t21*T0

            call self%eq_solver%solve(soln%eq_soln(idx), "tp", T2, soln%pressure(idx), weights, partials=soln%eq_partials(idx))
            if (.not. ShockSolver_state_valid(soln, idx)) then
                if (have_last_valid) then
                    soln%eq_soln(idx) = last_valid_eq
                    soln%eq_partials(idx) = last_valid_partials
                    rho12 = last_rho12
                    p21 = last_p21
                    t21 = last_t21
                    T2 = last_T2
                    wm_k = last_wm_k
                    exit
                else
                    call ShockSolver_fail_state(soln, idx, &
                                                "ShockSolver_solve_incident: incident equilibrium iteration failed.")
                    return
                end if
            end if
            if (.not. soln%eq_soln(idx)%converged) then
                call ShockSolver_finalize_equilibrium_state(self, soln, idx, update_transport=.false.)
            end if

            ! Update properties after the equilibrium shock
            wm_k = 1.0d0/soln%eq_soln(idx)%n
            nj_g = soln%eq_soln(idx)%nj(1:self%eq_solver%num_gas)
            cp = soln%eq_soln(idx)%cp_eq/(R*1.d-3)
            h0 = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%enthalpy)*T2
            dlnV_dlnP = soln%eq_partials(idx)%dlnV_dlnP
            dlnV_dlnT = soln%eq_partials(idx)%dlnV_dlnT

            rho12 = wm*t21/(p21*wm_k)
            last_valid_eq = soln%eq_soln(idx)
            last_valid_partials = soln%eq_partials(idx)
            last_rho12 = rho12
            last_p21 = p21
            last_t21 = t21
            last_T2 = T2
            last_wm_k = wm_k
            have_last_valid = .true.

            tmp = rho12*mu12rt

            ! Compute the solution matrix
            G(1,1) = -tmp*dlnV_dlnP - p21
            G(1,2) = -tmp*dlnV_dlnT
            G(1,3) = p21 - 1.0d0 + tmp - mu12rt

            tmp = tmp*T0/wm
            tmp = tmp*rho12

            G(2,1) = -tmp*dlnV_dlnP + T2*(dlnV_dlnT-1.d0)/wm_k
            G(2,2) = -tmp*dlnV_dlnT - T2*cp
            tmp = 1.0d0 - rho12**2
            G(2,3) = h0 - h_init - u1**2*tmp/(2.0d0*R)

            ! Solve the solution vector directly
            X(3) = G(1,1)*G(2,2) - G(1,2)*G(2,1)
            X(1) = (G(1,3)*G(2,2)-G(2,3)*G(1,2))/X(3)
            X(2) = (G(1,1)*G(2,3)-G(2,1)*G(1,3))/X(3)

            ! Compute the damped update factor, apply the solution update, and check convergence
            call self%update_solution(soln, X(1), X(2), p21, t21, i)

            if (i == 1 .and. .not. soln%converged .and. t21 >= ttmax) then
                call ShockSolver_fail_state(soln, idx, "ShockSolver_solve_incident: first-iteration update hit " // &
                                            "temperature cap; marking incident point as failed.")
                return
            end if

            ! Convergence check and exit
            if (soln%converged) then
                soln%rho12 = rho12
                soln%p21 = p21
                soln%t21 = t21
                soln%u(idx) = u1*rho12
                soln%M21 = wm_k/wm
                soln%v2 = u1 - soln%u(idx)
                soln%v_sonic(idx) = (R*T2*soln%eq_soln(idx)%gamma_s/wm_k)**0.5d0
                soln%mach(idx) = soln%u(idx)/soln%v_sonic(idx)
                exit
            end if

        end do

        ! For non-converged incident-equilibrium states, recompute transport
        ! using elemental component rows, even if shock-level variables converge.
        if (ShockSolver_state_valid(soln, idx)) then
            if (self%eq_solver%transport .and. .not. soln%eq_soln(idx)%converged) then
                call ShockSolver_apply_atomic_transport_basis(self, soln%eq_soln(idx))
                call ShockSolver_update_transport(self, soln%eq_soln(idx), update_basis=.false.)
            end if
        end if

        ! Keep the last valid incident-equilibrium state after
        ! singular-recovery exits and continue with warning semantics.
        if (.not. soln%converged) then
            if (ShockSolver_state_valid(soln, idx)) then
                ! Preserve the retained incident state for inspection, but do
                ! not mark it as converged; callers must check converged before
                ! treating this point as a successful shock solve.
                soln%solve_status = SHOCK_SOLVE_STATUS_LAST_VALID
                soln%pressure(idx) = p21*P0
                soln%rho12 = rho12
                soln%p21 = p21
                soln%t21 = t21
                soln%u(idx) = u1*rho12
                soln%M21 = wm_k/wm
                soln%v2 = u1 - soln%u(idx)
                soln%v_sonic(idx) = (R*T2*soln%eq_soln(idx)%gamma_s/wm_k)**0.5d0
                soln%mach(idx) = soln%u(idx)/soln%v_sonic(idx)
            else
                call ShockSolver_fail_state(soln, idx, "ShockSolver_solve_incident: incident equilibrium solve did not converge.")
            end if
        end if

    end subroutine

    subroutine ShockSolver_solve_incident_frozen(self, soln, weights, T0, P0)
        ! Solve the incident shock conditions

        ! Arguments
        class(ShockSolver), intent(in) :: self
        class(ShockSolution), intent(inout) :: soln
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: T0                     ! Initial reactant temperature [K]
        real(dp), intent(in) :: P0                     ! Initial reactant pressure [bar]

        ! Locals
        integer :: idx  ! Solution index for the incident conditions
        real(dp) :: mach1, u1             ! Incident Mach and velocity
        integer :: i, j                   ! Loop index
        real(dp) :: gamma1, gamma2        ! Ratio of specific heats at initial/shocked conditions
        real(dp) :: cp                    ! Mixture heat capacity
        real(dp) :: wm, wm_k              ! Mixture molecular weight (initial, k-th iteration)
        real(dp) :: h_init, h0            ! Mixture enthalpy (initial, <all other points>)
        real(dp) :: T2                    ! Temperature after incident, reflected shocks [K]
        real(dp) :: p21, t21, ttmax       ! Pressure/temperature ratio across the incident shock
        real(dp) :: G(2, 3)               ! Solution matrix
        real(dp) :: X(3)                  ! Solution vector
        real(dp) :: dlnV_dlnP, dlnV_dlnT  ! Partial derivatives
        real(dp) :: mu12rt, rho12         ! Ratios of chemical potential and density across the incident shock
        real(dp) :: tmp                   ! Intermediate variabls
        integer, parameter :: max_iter = 60
        real(dp), parameter :: T_gas_max = 20000.d0  ! Max gas temperature in the thermo database [K]
        character(len=2), parameter :: type="tp"

        ! Initialize
        idx = 2
        G = 0.d0
        soln%converged = .false.
        soln%solve_status = SHOCK_SOLVE_STATUS_NOT_CONVERGED
        u1 = soln%u(1)
        mach1 = soln%mach(1)
        soln%eq_soln(idx) = EqSolution(self%eq_solver, T_init=T0)
        call soln%eq_soln(idx)%constraints%set(type, T0, P0, &
            self%eq_solver%reactants%element_amounts_from_weights(weights))

        ! Set the reactant weights as the species amount
        soln%eq_soln(idx)%converged = .true.
        soln%eq_soln(idx)%nj(:) = 0.0d0
        soln%eq_soln(idx)%ln_nj(:) = self%eq_solver%log_min
        do i = 1, self%eq_solver%num_reactants
            j = findloc(self%eq_solver%products%species_names, self%eq_solver%reactants%species_names(i), 1)
            if (j > 0) then
                soln%eq_soln(idx)%nj(j) = (weights(i)/self%eq_solver%reactants%species(i)%molecular_weight)/sum(weights)
                soln%eq_soln(idx)%ln_nj(j) = log(soln%eq_soln(idx)%nj(j))
            else
                call log_warning("ShockSolver_solve_incident_frozen: Reactant not found in products.")
            end if
        end do

        soln%eq_partials(idx)%dlnV_dlnP = -1.0d0
        soln%eq_partials(idx)%dlnV_dlnT = 1.0d0

        ! Compute the molecular weight of the initial mixture
        wm = sum(weights)
        wm_k = wm
        soln%eq_soln(idx)%n = 1.0d0/wm

        ! Compute properties of the initial mixture
        cp = self%eq_solver%reactants%calc_frozen_cp(weights, T0)/R
        gamma1 = cp/(cp - 1.0/wm)
        h_init = self%eq_solver%reactants%calc_enthalpy(weights, T0)/R

        ! Compute the solution for the incident shock
        p21 = (2.0d0*gamma1*mach1**2.d0 - gamma1 + 1.d0)/(gamma1 + 1.d0)
        h0 = h_init + u1**2/(2.d0*R)
        mu12rt = wm*u1**2/(R*T0)

        soln%pressure(idx) = p21*P0

        t21 = p21*(2.0d0/mach1**2+gamma1 - 1.0d0)/(gamma1 + 1.0d0)
        ttmax = 1.05*T_gas_max/T0
        t21 = min(t21, ttmax)

        do i = 1, max_iter
            ! Update the pressure
            soln%pressure(idx) = p21*P0
            T2 = t21*T0
            soln%eq_soln(idx)%T = T2
            soln%eq_soln(idx)%constraints%state2 = soln%pressure(idx)

            ! Update properties after the equilibrium shock
            cp = self%eq_solver%reactants%calc_frozen_cp(weights, T2)/R
            h0 = self%eq_solver%reactants%calc_enthalpy(weights, T2)/R
            dlnV_dlnP = soln%eq_partials(idx)%dlnV_dlnP
            dlnV_dlnT = soln%eq_partials(idx)%dlnV_dlnT

            rho12 = wm*t21/(p21*wm_k)
            soln%rho12 = rho12
            tmp = rho12*mu12rt

            ! Compute the solution matrix
            G(1,1) = -tmp*dlnV_dlnP - p21
            G(1,2) = -tmp*dlnV_dlnT
            G(1,3) = p21 - 1. + tmp - mu12rt

            tmp = tmp*T0/wm
            tmp = tmp*rho12

            G(2,1) = -tmp*dlnV_dlnP + T2*(dlnV_dlnT-1.)/wm_k
            G(2,2) = -tmp*dlnV_dlnT - T2*cp
            tmp = 1. - rho12**2
            G(2,3) = h0 - h_init - u1**2*tmp/(2.*R)

            ! Solve the solution vector directly
            X(3) = G(1,1)*G(2,2) - G(1,2)*G(2,1)
            X(1) = (G(1,3)*G(2,2)-G(2,3)*G(1,2))/X(3)
            X(2) = (G(1,1)*G(2,3)-G(2,1)*G(1,3))/X(3)

            ! Compute the damped update factor, apply the solution update, and check convergence
            call self%update_solution(soln, X(1), X(2), p21, t21, i)

            if (i == 1 .and. .not. soln%converged .and. t21 >= ttmax) then
                call ShockSolver_fail_state(soln, idx, "ShockSolver_solve_incident_frozen: first-iteration update hit " // &
                                            "temperature cap; marking incident point as failed.")
                return
            end if

            ! Convergence check
            if (soln%converged) then
                call self%eq_solver%products%calc_thermo(soln%eq_soln(idx)%thermo, soln%eq_soln(idx)%T, condensed=.false.)
                call ShockSolver_update_transport(self, soln%eq_soln(idx), frozen_shock=.true.)
                gamma2 = cp/(cp - 1.0d0/wm_k)
                soln%rho12 = rho12
                soln%p21 = p21
                soln%t21 = t21
                soln%u(idx) = u1*rho12
                soln%M21 = wm_k/wm
                soln%v2 = u1 - soln%u(idx)
                soln%eq_partials(idx)%gamma_s = gamma2
                soln%eq_soln(idx)%gamma_s = gamma2
                soln%v_sonic(idx) = (R*T2*gamma2/wm_k)**0.5d0
                soln%mach(idx) = soln%u(idx)/soln%v_sonic(idx)
                exit
            end if

        end do

        ! Not converged; compute shock properties
        if (.not. soln%converged) then
            call ShockSolver_fail_state(soln, idx, "ShockSolver_solve_incident_frozen: incident frozen solve did not converge.")
        end if

    end subroutine

    subroutine ShockSolver_solve_reflected(self, soln, weights, T0, P0)
        ! Solve the reflected shock conditions

        ! Arguments
        class(ShockSolver), intent(in) :: self
        class(ShockSolution), intent(inout) :: soln
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: T0                     ! Initial reactant temperature [K]
        real(dp), intent(in) :: P0                     ! Initial reactant pressure [bar]

        ! Locals
        integer :: idx  ! Solution index for the incident conditions
        integer :: i                      ! Loop index
        real(dp) :: gamma1                ! Ratio of specific heats at incident condition
        real(dp) :: cp                    ! Mixture heat capacity
        real(dp) :: wm, wm_k              ! Mixture molecular weight (initial, k-th iteration)
        real(dp) :: h_init, h0            ! Mixture enthalpy (initial, <all other points>)
        real(dp) :: u1                    ! Incident shock velocity
        real(dp) :: a1                    ! Provisional sound speed
        real(dp) :: T2, T5                ! Temperature after incident, reflected shocks [K]
        real(dp) :: p52, t52, ttmax       ! Pressure/temperature ratio across the reflected shock
        real(dp) :: b5                    ! Intermediate variable for reflected shock initial state
        real(dp) :: G(2, 3)               ! Solution matrix
        real(dp) :: X(3)                  ! Solution vector
        real(dp) :: dlnV_dlnP, dlnV_dlnT  ! Partial derivatives
        real(dp) :: rho12                 ! Incident-shock density ratio carried into reflected solve
        real(dp) :: mu25rt, rho52         ! Ratios of chemical potential and density across the reflected shock
        real(dp) :: tmp                   ! Intermediate variabls
        real(dp), allocatable :: nj_g(:)  ! Total/gas species concentrations [kmol-per-kg]
        integer, parameter :: max_iter = 60
        real(dp), parameter :: T_gas_max = 20000.d0  ! Max gas temperature in the thermo database [K]

        ! Initialize
        if (.false.) print *, P0
        idx = 3
        G = 0.0d0  ! Reset the matrix
        soln%converged = .false.
        soln%solve_status = SHOCK_SOLVE_STATUS_NOT_CONVERGED
        soln%eq_soln(idx) = EqSolution(self%eq_solver, T_init=soln%eq_soln(2)%T, nj_init=soln%eq_soln(2)%nj)

        ! Retrieve values from the incident solution
        u1 = soln%u(1)
        T2 = soln%eq_soln(2)%T
        rho12 = soln%rho12

        ! Initialize the solution for the reflected shock
        wm = 1.0d0/soln%eq_soln(2)%n
        gamma1 = soln%eq_soln(2)%gamma_s
        h_init = dot_product(soln%eq_soln(2)%nj, soln%eq_soln(2)%thermo%enthalpy)*T2
        mu25rt = wm*(u1 - u1*rho12)**2/(R*soln%eq_soln(2)%T)
        t52 = 2.0d0
        b5 = (-1.d0 - mu25rt - t52)/2.0d0
        p52 = -b5 + sqrt(b5**2 - t52)
        ttmax = 1.05*T_gas_max/T2
        t52 = min(t52, ttmax)

        do i = 1, max_iter
            ! Update the pressure
            soln%pressure(idx) = p52*soln%pressure(2)
            T5 = t52*T2

            call self%eq_solver%solve(soln%eq_soln(idx), "tp", T5, soln%pressure(idx), weights, partials=soln%eq_partials(idx))
            if ((.not. ShockSolver_state_valid(soln, idx)) .or. (.not. soln%eq_soln(idx)%converged)) then
                call ShockSolver_fail_state(soln, idx, &
                                            "ShockSolver_solve_reflected: reflected equilibrium iteration failed.")
                return
            end if

            ! Update properties after the equilibrium shock
            wm_k = 1.0d0/soln%eq_soln(idx)%n
            nj_g = soln%eq_soln(idx)%nj(1:self%eq_solver%num_gas)
            cp = soln%eq_soln(idx)%cp_eq/(R*1.d-3)
            call self%eq_solver%products%calc_thermo(soln%eq_soln(idx)%thermo, soln%eq_soln(idx)%T, condensed=.true.)
            h0 = dot_product(soln%eq_soln(idx)%nj, soln%eq_soln(idx)%thermo%enthalpy)*T5
            dlnV_dlnP = soln%eq_partials(idx)%dlnV_dlnP
            dlnV_dlnT = soln%eq_partials(idx)%dlnV_dlnT

            ! Retain the legacy provisional values during iteration. Final
            ! reflected-frame speed, sound speed, and Mach are set on convergence.
            a1 = (R*gamma1*T0/wm)**0.5d0
            soln%u(idx) = u1*rho12
            soln%v_sonic(idx) = a1

            rho52 = 1./(wm*t52/(p52*wm_k))
            tmp = -mu25rt*rho52/(rho52 - 1.0d0)**2

            ! Compute the solution matrix
            G(1,1) = -tmp*dlnV_dlnP - p52
            G(1,2) = -tmp*dlnV_dlnT
            G(1,3) = p52 - 1.0d0 + tmp*(rho52-1.)

            tmp = tmp*T2/wm
            G(2,1) = -tmp*dlnV_dlnP + T5*(dlnV_dlnT-1.)/wm_k
            G(2,2) = -tmp*dlnV_dlnT - T5*cp
            tmp = (rho52 + 1.0d0)/(rho52 - 1.0d0)
            G(2,3) = h0 - h_init - (u1 - u1*rho12)**2*tmp/(2.*R)

            ! Solve the solution vector directly
            X(3) = G(1,1)*G(2,2) - G(1,2)*G(2,1)
            X(1) = (G(1,3)*G(2,2)-G(2,3)*G(1,2))/X(3)
            X(2) = (G(1,1)*G(2,3)-G(2,1)*G(1,3))/X(3)

            ! Compute the damped update factor
            call self%update_solution(soln, X(1), X(2), p52, t52, i)

            if (i == 1 .and. .not. soln%converged .and. t52 >= ttmax) then
                call ShockSolver_fail_state(soln, idx, "ShockSolver_solve_reflected: first-iteration update hit " // &
                                            "temperature cap; marking reflected point as failed.")
                return
            end if

            ! Convergence check
            if (soln%converged) then
                soln%rho52 = rho52
                soln%p52 = p52
                soln%t52 = t52
                soln%u(idx) = (u1 - u1*rho12)/(rho52 - 1.0d0)
                soln%M52 = wm_k/wm
                soln%u5_p_v2 = (u1 - u1*rho12)*rho52/(rho52 - 1.0d0)
                soln%v_sonic(idx) = (R*T5*soln%eq_soln(idx)%gamma_s/wm_k)**0.5d0
                soln%mach(idx) = soln%u(idx)/soln%v_sonic(idx)
                exit
            end if

        end do

        ! Not converged; compute shock properties
        if (.not. soln%converged) then
            call ShockSolver_fail_state(soln, idx, "ShockSolver_solve_reflected: reflected equilibrium solve did not converge.")
        end if

    end subroutine

    subroutine ShockSolver_solve_reflected_frozen(self, soln, weights, T0, P0)
        ! Solve the incident shock conditions

        ! Arguments
        class(ShockSolver), intent(in) :: self
        class(ShockSolution), intent(inout) :: soln
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: T0                     ! Initial reactant temperature [K]
        real(dp), intent(in) :: P0                     ! Initial reactant pressure [bar]

        ! Locals
        integer :: idx  ! Solution index for the incident conditions
        integer :: i                      ! Loop index
        integer :: ng                     ! Number of gas species
        real(dp) :: cp, gamma5            ! Mixture heat capacity, ratio of specific heats
        real(dp) :: wm, wm_k              ! Mixture molecular weight (initial, k-th iteration)
        real(dp) :: h_init, h0            ! Mixture enthalpy (initial, <all other points>)
        real(dp) :: u1                    ! Incident shock velocity
        real(dp) :: u2, u1u2              ! Reflected shock velocity, difference
        real(dp) :: T2, T5                ! Temperature after incident, reflected shocks [K]
        real(dp) :: p52, t52, ttmax       ! Pressure/temperature ratio across the reflected shock
        real(dp) :: b5                    ! Intermediate variable for reflected shock initial state
        real(dp) :: G(2, 3)               ! Solution matrix
        real(dp) :: X(3)                  ! Solution vector
        real(dp) :: dlnV_dlnP, dlnV_dlnT  ! Partial derivatives
        real(dp) :: rho12                 ! Ratios of chemical potential and density across the incident shock
        real(dp) :: mu25rt, rho52         ! Ratios of chemical potential and density across the reflected shock
        real(dp) :: rho25_inv             ! Reflected-shock inverse density ratio used in the Newton system
        real(dp) :: tmp                   ! Intermediate variabls
        integer, parameter :: max_iter = 60
        real(dp), parameter :: T_gas_max = 20000.d0  ! Max gas temperature in the thermo database [K]
        character(len=2), parameter :: type="tp"

        call log_debug("Calling ShockSolver_solve_reflected_frozen")

        ! Initialize
        idx = 3
        G = 0.0d0  ! Reset the matrix
        soln%eq_soln(idx) = EqSolution(self%eq_solver, T_init=T0)
        call soln%eq_soln(idx)%constraints%set(type, T0, P0, &
            self%eq_solver%reactants%element_amounts_from_weights(weights))
        soln%eq_partials(idx)%dlnV_dlnP = -1.0d0
        soln%eq_partials(idx)%dlnV_dlnT = 1.0d0
        soln%converged = .false.
        soln%solve_status = SHOCK_SOLVE_STATUS_NOT_CONVERGED

        ! Set the reactant weights as the species amount
        soln%eq_soln(idx)%converged = .true.
        soln%eq_soln(idx)%nj = soln%eq_soln(2)%nj
        soln%eq_soln(idx)%ln_nj = soln%eq_soln(2)%ln_nj
        soln%eq_soln(idx)%n = soln%eq_soln(2)%n
        soln%eq_soln(idx)%is_active = soln%eq_soln(2)%is_active
        soln%eq_soln(idx)%active_rank = soln%eq_soln(2)%active_rank
        soln%eq_soln(idx)%j_liq = soln%eq_soln(2)%j_liq
        soln%eq_soln(idx)%j_sol = soln%eq_soln(2)%j_sol
        soln%eq_soln(idx)%j_switch = soln%eq_soln(2)%j_switch

        ! Retrieve values from the incident solution
        u1 = soln%u(1)
        T2 = soln%eq_soln(2)%T
        rho12 = soln%rho12

        ! Compute the molecular weight of the initial mixture
        wm = 1.0d0/soln%eq_soln(2)%n
        wm_k = wm
        ng = self%eq_solver%num_gas

        ! Initialize the solution for the reflected shock
        h_init = dot_product(soln%eq_soln(2)%nj, soln%eq_soln(2)%thermo%enthalpy)*T2
        u2 = u1*rho12
        u1u2 = soln%v2
        mu25rt = wm*(u1 - u1*rho12)**2/(R*soln%eq_soln(2)%T)
        t52 = 2.0d0
        b5 = (-1.d0 - mu25rt - t52)/2.0d0
        p52 = -b5 + sqrt(b5**2 - t52)
        ttmax = 1.05*T_gas_max/T2
        t52 = min(t52, ttmax)

        do i = 1, max_iter
            ! Update the pressure
            soln%pressure(idx) = p52*soln%pressure(2)
            T5 = t52*T2
            soln%eq_soln(idx)%T = T5
            soln%eq_soln(idx)%constraints%state2 = soln%pressure(idx)
            call self%eq_solver%products%calc_thermo(soln%eq_soln(idx)%thermo, soln%eq_soln(idx)%T, condensed=.false.)

            ! Update frozen properties from the incident (state-2) frozen composition.
            cp = dot_product(soln%eq_soln(idx)%nj(:ng), soln%eq_soln(idx)%thermo%cp(:ng))
            h0 = dot_product(soln%eq_soln(idx)%nj(:ng), soln%eq_soln(idx)%thermo%enthalpy(:ng))*T5
            dlnV_dlnP = soln%eq_partials(idx)%dlnV_dlnP
            dlnV_dlnT = soln%eq_partials(idx)%dlnV_dlnT

            rho25_inv = wm*t52/(wm_k*p52)
            rho52 = 1./rho25_inv
            soln%rho52 = rho52
            tmp = -mu25rt*rho52/(rho52 - 1.0d0)**2

            ! Compute the solution matrix
            G(1,1) = -tmp*dlnV_dlnP - p52
            G(1,2) = -tmp*dlnV_dlnT
            G(1,3) = p52 - 1.0d0 + tmp*(rho52 - 1.0d0)

            tmp = tmp*T2/wm
            G(2,1) = -tmp*dlnV_dlnP + T5*(dlnV_dlnT-1.)/wm_k
            G(2,2) = -tmp*dlnV_dlnT - T5*cp
            tmp = (rho52 + 1.0d0)/(rho52 - 1.0d0)
            G(2,3) = h0 - h_init - soln%v2**2*tmp/(2.0d0*R)

            ! Solve the solution vector directly
            X(3) = G(1,1)*G(2,2) - G(1,2)*G(2,1)
            X(1) = (G(1,3)*G(2,2)-G(2,3)*G(1,2))/X(3)
            X(2) = (G(1,1)*G(2,3)-G(2,1)*G(1,3))/X(3)

            ! Compute the damped update factor
            call self%update_solution(soln, X(1), X(2), p52, t52, i)

            if (i == 1 .and. .not. soln%converged .and. t52 >= ttmax) then
                call ShockSolver_fail_state(soln, idx, "ShockSolver_solve_reflected_frozen: first-iteration update hit " // &
                                            "temperature cap; marking reflected point as failed.")
                return
            end if

            ! Convergence check
            if (soln%converged) then
                call ShockSolver_update_transport(self, soln%eq_soln(idx), frozen_shock=.true.)
                gamma5 = cp/(cp - 1.0d0/wm_k)
                soln%rho52 = rho52
                soln%p52 = p52
                soln%t52 = t52
                soln%u(idx) = soln%v2/(rho52 - 1.0d0)
                soln%M52 = wm_k/wm
                soln%u5_p_v2 = (u1 - u1*rho12)*rho52/(rho52 - 1.0d0)
                soln%eq_partials(idx)%gamma_s = gamma5
                soln%eq_soln(idx)%gamma_s = gamma5
                soln%v_sonic(idx) = (R*T5*gamma5/wm_k)**0.5d0
                soln%mach(idx) = soln%u(idx)/soln%v_sonic(idx)
                exit
            end if

        end do

        ! Not converged; compute shock properties
        if (.not. soln%converged) then
            call ShockSolver_fail_state(soln, idx, "ShockSolver_solve_reflected_frozen: reflected frozen solve did not converge.")
        end if

    end subroutine

    function ShockSolver_solve(self, reactant_weights, T0, P0, u1, mach1, reflected, incident_frozen, reflected_frozen) result(soln)
        ! Solve the moving shock problem

        ! Arguments
        class(ShockSolver), intent(in) :: self
        real(dp), intent(in) :: reactant_weights(:)
        real(dp), intent(in) :: T0                          ! Initial reactant temperature [K]
        real(dp), intent(in) :: P0                          ! Initial reactant pressure [bar]
        real(dp), intent(in), optional :: u1                ! Shock velocities [m/s]
        real(dp), intent(in), optional :: mach1             ! Shock Mach
        logical,  intent(in), optional :: reflected         ! Compute the solution for a relfected shock
        logical,  intent(in), optional :: incident_frozen   ! Use frozen analysis for the incident shock
        logical,  intent(in), optional :: reflected_frozen  ! Use frozen analysis for the reflected shock

        ! Result
        type(ShockSolution) :: soln

        ! Locals
        logical :: reflected_, incident_frozen_, reflected_frozen_    ! Problem flags
        real(dp) :: mach1_, u1_           ! Initial mach and velocity
        integer :: npts                   ! Number of problem types to solve
        integer :: i, j
        real(dp) :: gamma1                ! Ratio of specific heats at initial condition
        real(dp) :: cp                    ! Mixture heat capacity
        real(dp) :: wm                    ! Mixture molecular weight (initial, k-th iteration)
        real(dp) :: a1                    ! Initial speed of sound
        character(len=2), parameter :: type="tp"

        call log_info("Calling ShockSolver_solve")

        ! Index:
        ! 1: Unshocked gas
        ! 2: Incident shock
        ! 3: Reflected shock

        ! NOTE: Always solve the incident problem
        ! There are 6 possible solution options:
        ! 1. Incident equilibrium
        ! 2. Incident frozen
        ! 3. Incident equilibrium + reflected equilibrium
        ! 4. Incident equilibrium + reflected frozen
        ! 5. Incident frozen + reflected frozen
        ! 5. Incident frozen + reflected equilbrium

        ! NOTE: solution defaults to equilibrium analysis, which differs from CEA2 default.

        ! Input handling
        if ((present(u1)) .and. (present(mach1))) then
            call abort("ShockSolver_solve: u1 and mach1 cannot both be present")
        else if ((.not. present(u1)) .and. (.not. present(mach1))) then
            call abort("ShockSolver_solve: one of u1 or mach1 must be present")
        end if

        ! Set defaults for optional arguments
        if (present(reflected)) then
            reflected_ = reflected
        else
            reflected_ = .false.
        end if

        if (present(incident_frozen)) then
            incident_frozen_ = incident_frozen
        else
            incident_frozen_ = .true.
        end if

        if (present(reflected_frozen)) then
            reflected_frozen_ = reflected_frozen
        else
            reflected_frozen_ = .false.
        end if

        ! Get the number of points in the solution
        npts = 2
        if (reflected_) npts = npts + 1

        ! Initialize the solution
        soln = ShockSolution_init(npts)
        soln%solve_status = SHOCK_SOLVE_STATUS_NOT_CONVERGED

        ! Solve the problem
        ! --------------------------------------------------------------------
        ! First, compute the conditions of the unshocked gas

        ! Compute the molecular weight of the initial mixture
        wm = sum(reactant_weights)

        ! Compute properties of the initial mixture
        cp = self%eq_solver%reactants%calc_frozen_cp(reactant_weights, T0)/R
        gamma1 = cp/(cp - 1.0/wm)
        a1 = (R*gamma1*T0/wm)**0.5d0
        if (present(u1)) then
            u1_ = u1
            mach1_ = u1/a1
        else
            mach1_ = mach1
            u1_ = a1*mach1
        end if

        ! Store the unshocked gas solution
        soln%eq_soln(1) = EqSolution(self%eq_solver, T_init=T0)
        soln%pressure(1) = P0
        call soln%eq_soln(1)%constraints%set(type, T0, P0, &
            self%eq_solver%reactants%element_amounts_from_weights(reactant_weights))
        soln%eq_soln(1)%n = 1.0d0/wm
        soln%eq_partials(1)%gamma_s = gamma1
        soln%eq_soln(1)%gamma_s = gamma1

        ! Store the shocked velocity and Mach number
        soln%mach(1) = mach1_
        soln%u(1) = u1_
        soln%v_sonic(1) = a1

        ! Set the reactant weights as the species amount
        soln%eq_soln(1)%converged = .true.
        soln%eq_soln(1)%nj(:) = 0.0d0
        soln%eq_soln(1)%ln_nj(:) = self%eq_solver%log_min
        do i = 1, self%eq_solver%num_reactants
            j = findloc(self%eq_solver%products%species_names, self%eq_solver%reactants%species_names(i), 1)
            if (j > 0) then
                soln%eq_soln(1)%nj(j) = (reactant_weights(i)/ &
                    self%eq_solver%reactants%species(i)%molecular_weight)/sum(reactant_weights)
                soln%eq_soln(1)%ln_nj(j) = log(soln%eq_soln(1)%nj(j))
            else
                call log_warning("ShockSolver_solve_incident_frozen: Reactant not found in products.")
            end if
        end do

        ! Compute properties of the unshocked gas
        call self%eq_solver%post_process(soln%eq_soln(1))

        ! Compute the incident shock solution
        if (incident_frozen_) then
            call self%solve_incident_frozen(soln, reactant_weights, T0, P0)
        else  ! Equilibrium
            call self%solve_incident(soln, reactant_weights, T0, P0)
        end if
        if (soln%eq_soln(2)%T <= 0.0d0) then
            return
        end if
        if (soln%eq_soln(2)%converged) then
            call self%eq_solver%post_process(soln%eq_soln(2))
        end if
        if (.not. incident_frozen_ .and. reflected_ .and. reflected_frozen_ .and. self%eq_solver%transport) then
            call ShockSolver_apply_atomic_transport_basis(self, soln%eq_soln(2))
            call ShockSolver_update_transport(self, soln%eq_soln(2), update_basis=.false.)
        end if
        ! Compute the reflected shock solution
        if (reflected_) then
            if (reflected_frozen_) then
                call self%solve_reflected_frozen(soln, reactant_weights, T0, P0)
            else  ! Equilibrium
                call self%solve_reflected(soln, reactant_weights, T0, P0)
            end if
            if (soln%eq_soln(3)%T <= 0.0d0) then
                return
            end if
            if (reflected_frozen_) then
                call self%eq_solver%products%calc_thermo(soln%eq_soln(3)%thermo, soln%eq_soln(3)%T, condensed=.true.)
                call ShockSolver_finalize_reflected_frozen_state(self, soln, 3)
            else
                call self%eq_solver%post_process(soln%eq_soln(3))
            end if
        end if

    end function

    !-----------------------------------------------------------------------
    ! ShockSolution
    !-----------------------------------------------------------------------
    function ShockSolution_init(num_pts) result(self)
        type(ShockSolution) :: self
        integer, intent(in) :: num_pts

        self%num_pts = num_pts
        self%solve_status = SHOCK_SOLVE_STATUS_NOT_CONVERGED

        ! Allocate data structures
        allocate(self%eq_soln(num_pts))
        allocate(self%eq_partials(num_pts))
        allocate(self%pressure(num_pts))
        allocate(self%mach(num_pts))
        allocate(self%u(num_pts))
        allocate(self%v_sonic(num_pts))
    end function

end module
