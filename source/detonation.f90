module cea_detonation
    !! Detonation solver module

    use cea_param, only: dp, empty_dp, R=>gas_constant
    use cea_mixture, only: Mixture, MixtureThermo
    use cea_transport, only: TransportDB
    use cea_equilibrium, only: EqSolution, EqSolver, EqPartials
    use fb_utils
    implicit none

    type :: DetonSolver
        !! Detonation solver object

        type(EqSolver) :: eq_solver
            !! Equilibrium solver

    contains

        procedure :: solve_eq => DetonSolver_solve_eq
        ! procedure :: solve_frozen => DetonSolver_solve_frozen
        procedure :: solve => DetonSolver_solve

    end type
    interface DetonSolver
        module procedure :: DetonSolver_init
    end interface

    type :: DetonSolution
        !! Detonation solution object

        ! Thermo solution at each evaluation point
        type(EqSolution) :: eq_soln
            !! Equilibrium solution at each evaluation point
        type(EqPartials) :: eq_partials
            !! Partial derivatives of the equilibrium solution at each point

        ! States (unburned gas)
        real(dp) :: P1
        real(dp) :: T1
        real(dp) :: H1
        real(dp) :: M1
        real(dp) :: gamma1
        real(dp) :: v_sonic1

        ! Detonation parameters
        real(dp) :: pressure
            !! Pressure [bar]
        real(dp) :: sonic_velocity
            !! Speed of sound [m/s]
        real(dp) :: velocity
            !! Detonation velocity [m/s]
        real(dp) :: mach
            !! Detonation Mach number
        real(dp) :: gamma
            !! Ratio of specific heats
        real(dp) :: enthalpy
            !! Enthalpy of the detonation
        real(dp) :: P_P1
            !! P/P1
        real(dp) :: T_T1
            !! T/T1
        real(dp) :: M_M1
            !! M/M1
        real(dp) :: rho_rho1
            !! rho/rho1

        ! Convergence variables
        logical :: converged = .false.
            !! Convergence flag

    end type
    interface DetonSolution
        module procedure :: DetonSolution_init
    end interface

contains

    !-----------------------------------------------------------------------
    ! DetonSolver
    !-----------------------------------------------------------------------
    function DetonSolver_init(products, reactants, trace, ions, all_transport, insert) result(self)
        type(DetonSolver) :: self
        type(Mixture), intent(in) :: products
        type(Mixture), intent(in), optional :: reactants
        real(dp), intent(in), optional :: trace
        logical, intent(in), optional :: ions
        type(TransportDB), intent(in), optional :: all_transport
        character(*), intent(in), optional :: insert(:)  ! List of condensed species to insert

        ! Initialize the equilibrium solver
        self%eq_solver = EqSolver(products, reactants, trace, ions, all_transport, insert)

    end function

    subroutine DetonSolver_solve_eq(self, soln, weights, t1, p1)
        ! Solve the detonation problem with equilibrium analysis

        ! Arguments
        class(DetonSolver), intent(in) :: self
        class(DetonSolution), intent(inout) :: soln
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: t1                     ! Initial reactant temperature [K]
        real(dp), intent(in) :: p1                     ! Initial reactant pressure [bar]

        ! Locals
        integer :: i                      ! Loop index
        real(dp) :: wm, wm_k              ! Mixture molecular weight (initial, current iteration)
        real(dp) :: h0, h_ref             ! Enthalpy of initial mixture, enthalpy of detonation for hp equilibrium
        real(dp) :: pp1, tt1              ! Pressure/temperature ratios
        real(dp) :: t2                    ! Temperature of the detonation
        real(dp) :: temp, amm, alfa, rk, rr1, alam  ! Gibberish variable names copied from CEA2
        real(dp) :: a11, a12, a21, a22, b1, b2, d, x1, x2  ! Newton-Raphson update variables
        real(dp) :: cp, gamma
        real(dp) :: dlnV_dlnP, dlnV_dlnT  ! Partial derivatives

        integer, parameter :: max_iter = 8

        call log_info("DETON: Starting detonation solve with equilibrium analysis")

        ! Compute unburned gas properties
        soln%T1 = t1
        soln%eq_soln%T = t1
        soln%P1 = p1
        soln%H1 = (self%eq_solver%reactants%calc_enthalpy(weights, t1) - &
            self%eq_solver%reactants%calc_enthalpy(weights, 298.15d0))/1.d3
        cp = self%eq_solver%reactants%calc_frozen_cp(weights, t1)
        wm = sum(weights)/sum(self%eq_solver%reactants%moles_from_weights(weights))
        soln%M1 = wm
        soln%gamma1 = cp/(cp-R/wm)
        soln%v_sonic1 = (R*soln%gamma1*t1/wm)**0.5d0

        ! Compute initial estimate of enthalpy after combustion
        h0 = soln%H1/(R/1.d3)
        pp1 = 15.0
        soln%pressure = p1*pp1
        h_ref = h0 + 0.75d0*t1*pp1/wm
        call self%eq_solver%solve(soln%eq_soln, "hp", h_ref, soln%pressure, weights, partials=soln%eq_partials)

        gamma = soln%eq_partials%gamma_s
        cp = soln%eq_partials%cp_eq
        wm_k = 1.0d0/soln%eq_soln%n
        t2 = soln%eq_soln%T
        tt1 = t2/t1
        temp = tt1 - 0.75d0*pp1/(cp*wm)
        amm = wm_k/wm

        ! Improve initial estimate using recursion scheme
        do i = 1, 3
            alfa = amm/tt1
            pp1 = (1.+gamma)*(1.+(1.-4.*gamma*alfa/(1.+gamma)**2)**.5)/(2.*gamma*alfa)
            rk = pp1*alfa
            tt1 = temp + .5*pp1*gamma*(rk*rk-1.)/(wm*cp*rk)
        end do

        t2 = t1*tt1
        rr1 = pp1*amm/tt1

        ! Newton-Raphson iterations until convergence
        do i = 1, max_iter

            soln%pressure = p1*pp1

            call log_debug("DETON: Starting new iteration")

            call self%eq_solver%solve(soln%eq_soln, "tp", t2, soln%pressure, weights, partials=soln%eq_partials)

            if (soln%eq_soln%converged) then
                call log_debug("DETON: Equilibrium solver converged")
            else
                call log_error("DETON: Equilibrium solver did not converge")
            end if

            ! Solve the update the log of the pressure and pressure ratios
            gamma = soln%eq_partials%gamma_s
            cp = soln%eq_partials%cp_eq
            wm_k = 1.0d0/soln%eq_soln%n
            h_ref = dot_product(soln%eq_soln%nj, soln%eq_soln%thermo%enthalpy)*t2
            amm = wm_k/wm
            rr1 = pp1*amm/tt1
            dlnV_dlnP = soln%eq_partials%dlnV_dlnP
            dlnV_dlnT = soln%eq_partials%dlnV_dlnT

            a11 = 1./pp1 + gamma*rr1*dlnV_dlnP
            a12 = gamma*rr1*dlnV_dlnT
            a21 = .5*gamma*(rr1**2-1.-dlnV_dlnP*(1.+rr1**2)) + dlnV_dlnT - 1.
            a22 = -.5*gamma*dlnV_dlnT*(rr1**2+1.) - wm_k*cp
            b1 = (1./pp1) - 1. + (gamma*(rr1-1.) )
            b2 = wm_k*(h_ref-h0)/t2 - .5*gamma*(rr1*rr1-1.)
            d = a11*a22 - a12*a21
            x1 = (a22*b1-a12*b2)/d
            x2 = (a11*b2-a21*b1)/d

            alam = 1.
            temp = x1
            IF ( temp < 0.0d0 ) temp = -temp
            IF ( x2 > temp ) temp = x2
            IF ( -x2 > temp ) temp = -x2
            IF ( temp > 0.4054652 ) alam = .4054652/temp
            pp1 = pp1*exp(x1*alam)
            tt1 = tt1*exp(x2*alam)

            t2 = t1*tt1

            ! Compute detonation parameters
            soln%velocity = rr1*(R*gamma*t2/wm_k)**.5
            soln%sonic_velocity = (gamma*R*t2/wm_k)**.5
            soln%mach = soln%velocity/(soln%gamma1*R*soln%T1/wm)**.5

            ! Store solution values
            soln%eq_soln%T = t2
            soln%gamma = gamma
            soln%enthalpy = h_ref

            ! Check convergence
            if (temp < 0.5d-4) then
                soln%P_P1 = pp1
                soln%T_T1 = tt1
                soln%M_M1 = (1.0d0/soln%eq_soln%n)/soln%M1
                soln%rho_rho1 = rr1
                soln%converged = .true.
                exit
            end if

        end do

        ! Not converged; compute detonation parameters before exiting
        soln%P_P1 = pp1
        soln%T_T1 = tt1
        soln%M_M1 = (1.0d0/soln%eq_soln%n)/soln%M1
        soln%rho_rho1 = rr1

    end subroutine

    ! subroutine DetonSolver_solve_frozen(self, soln, weights, t1, p1, v, mach)
    !     ! Solve the detonation problem with frozen composition

    !     ! Arguments
    !     class(DetonSolver), intent(in) :: self
    !     class(DetonSolution), intent(inout) :: soln
    !     real(dp), intent(in) :: weights(:)
    !     real(dp), intent(in) :: t1                     ! Initial reactant temperature [K]
    !     real(dp), intent(in) :: p1                     ! Initial reactant pressure [bar]
    !     real(dp), intent(in), optional :: v            ! Detonation velocities [m/s]
    !     real(dp), intent(in), optional :: mach         ! Detontion Mach

    !     ! Locals
    !     integer :: i                      ! Loop index
    !     real(dp) :: wm, wm_k              ! Mixture molecular weight (initial, current iteration)
    !     real(dp) :: h0, h_ref             ! Enthalpy of initial mixture, enthalpy of detonation for hp equilibrium
    !     real(dp) :: pp1, tt1              ! Pressure/temperature ratios
    !     real(dp) :: t2                    ! Temperature of the detonation
    !     real(dp) :: temp, amm, alfa, rk, rr1, alam  ! Gibberish variable names copied from CEA2
    !     real(dp) :: a11, a12, a21, a22, b1, b2, d, x1, x2  ! Newton-Raphson update variables
    !     real(dp) :: cp, gamma
    !     real(dp) :: dlnV_dlnP, dlnV_dlnT  ! Partial derivatives
    !     logical :: frozen_

    !     integer, parameter :: max_iter = 8

    !     call log_info("Starting Detonation Frozen solve.")

    ! end subroutine

    function DetonSolver_solve(self, reactant_weights, t1, p1, frozen) result(soln)
        ! Solve the detonation problem

        ! Arguments
        class(DetonSolver), intent(in) :: self
        real(dp), intent(in) :: reactant_weights(:)
        real(dp), intent(in) :: t1                     ! Initial reactant temperature [K]
        real(dp), intent(in) :: p1                     ! Initial reactant pressure [bar]
        logical, intent(in),  optional :: frozen       ! Frozen detonation flag

        ! Result
        type(DetonSolution) :: soln

        ! Locals
        logical :: frozen_

        call log_info("Starting Detonation solve.")

        ! TODO: detonation problems are limited to gaseous reactants; add a check

        ! Initialize the solution
        soln = DetonSolution_init()
        soln%eq_soln = EqSolution(self%eq_solver)

        ! Set the frozen flag
        if (present(frozen)) then
            frozen_ = frozen
        else
            frozen_ = .false.
        end if

        ! Call the solver
        if (frozen_) then
            call log_info('DetonSolver: frozen composition not supported yet')
            ! call self%solve_frozen(soln, reactant_weights, t1, p1)
        else
            call self%solve_eq(soln, reactant_weights, t1, p1)
        end if

    end function

    !-----------------------------------------------------------------------
    ! DetonSolution
    !-----------------------------------------------------------------------
    function DetonSolution_init() result(self)
        type(DetonSolution) :: self
    end function

end module
