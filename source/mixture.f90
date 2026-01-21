module cea_mixture
    !! Mixture module to compute and store properties of the reactant or product mixtures

    use cea_param, snl=>species_name_len, &
                   enl=>element_name_len
    use cea_thermo, only: ThermoDB, SpeciesThermo, build_elem_list
    use cea_input, only: Formula, ReactantInput
    use cea_atomic_data, only: get_atom_valence, get_atom_weight
    use cea_units, only: convert_units_to_si
    use iso_c_binding
    use fb_findloc, only: findloc
    use fb_utils, only: abort, is_empty
    implicit none


    !-----------------------------------------------------------------------
    ! Mixture Type
    !-----------------------------------------------------------------------
    type :: Mixture
        !! Computes properties of a mixture of gaseous/condensed species

        ! Sizes
        integer :: num_species   = empty_int
            !! Total number of species in the mixture
        integer :: num_elements  = empty_int
            !! Number of elements in the mixture
        integer :: num_gas       = empty_int
            !! Number of gaseous species in the mixture
        integer :: num_condensed = empty_int
            !! Number of condensed species in the mixture

        ! Species Data
        type(SpeciesThermo), allocatable :: species(:)
            !! Array of species thermo in the mixture
        logical, allocatable :: is_condensed(:)
            !! True if species is condensed
        real(dp), allocatable :: stoich_matrix(:,:)
            !! Stoichiometric matrix of species in the mixture

        ! Meta Data
        ! WORKAROUND: gfortran 9.3 fails unit tests if these are character(:)
        ! TODO: Add null terminators to names so they can be used in the C API
        character(snl), allocatable :: species_names(:)
            !! Names of species in the mixture
        character(enl), allocatable :: element_names(:)
            !! Names of elements in the mixture

        ! Options
        logical :: ions = .false.
            !! True if the mixture includes ions

    contains
        procedure :: elements_from_species => mixture_elements_from_species
        procedure :: calc_thermo => mixture_calc_thermo
        procedure :: get_products => mixture_get_products
        procedure :: weights_from_of => mixture_weights_from_of
        procedure :: of_from_equivalence => mixture_chem_eq_ratio_to_of_ratio
        procedure :: equivalence_from_of => mixture_of_ratio_to_chem_eq_ratio
        procedure :: of_from_phi => mixture_weight_eq_ratio_to_of_ratio
        procedure :: phi_from_of => mixture_of_ratio_to_weight_eq_ratio
        procedure :: weights_from_moles => mixture_weights_from_moles
        procedure :: moles_from_weights => mixture_moles_from_weights
        procedure :: per_weight_from_per_mole => mixture_per_weight_from_per_mole
        procedure :: per_mole_from_per_weight => mixture_per_mole_from_per_weight
        procedure :: get_valence => mixture_get_valence
        procedure :: element_amounts_from_weights => mixture_element_amounts_from_weights

        generic   :: calc_enthalpy     => mixture_calc_enthalpy_single, &
                                          mixture_calc_enthalpy_multi
        generic   :: calc_energy       => mixture_calc_energy_single, &
                                          mixture_calc_energy_multi
        generic   :: calc_entropy      => mixture_calc_entropy_single, &
                                          mixture_calc_entropy_multi
        generic   :: calc_gibbs_energy => mixture_calc_gibbs_energy_single, &
                                          mixture_calc_gibbs_energy_multi
        generic   :: calc_frozen_cp    => mixture_calc_frozen_cp_single, &
                                          mixture_calc_frozen_cp_multi
        generic   :: calc_frozen_cv    => mixture_calc_frozen_cv_single, &
                                          mixture_calc_frozen_cv_multi
        generic   :: calc_pressure     => mixture_calc_pressure_single, &
                                          mixture_calc_pressure_multi

        procedure :: mixture_calc_enthalpy_single,      mixture_calc_enthalpy_multi
        procedure :: mixture_calc_energy_single,        mixture_calc_energy_multi
        procedure :: mixture_calc_entropy_single,       mixture_calc_entropy_multi
        procedure :: mixture_calc_gibbs_energy_single,  mixture_calc_gibbs_energy_multi
        procedure :: mixture_calc_frozen_cp_single,     mixture_calc_frozen_cp_multi
        procedure :: mixture_calc_frozen_cv_single,     mixture_calc_frozen_cv_multi
        procedure :: mixture_calc_pressure_single,      mixture_calc_pressure_multi

    end type
    interface Mixture
        module procedure :: mixture_init
    end interface

    !-----------------------------------------------------------------------
    ! MixtureThermo
    !-----------------------------------------------------------------------
    type :: MixtureThermo
        !! Storage for thermodynamic properties computed by Mixture type
        integer :: num_species = empty_int
            !! Number of species in the mixture
        real(dp), allocatable :: cp(:)
            !! Specific heat at constant pressure (frozen)
        real(dp), allocatable :: cv(:)
            !! Specific heat at constant volume
        real(dp), allocatable :: enthalpy(:)
            !! Mixture enthalpy
        real(dp), allocatable :: entropy(:)
            !! Mixture entropy
        real(dp), allocatable :: energy(:)
            !! Mixture internal energy
    end type
    interface MixtureThermo
        module procedure :: mixturethermo_init
    end interface

contains

    !-----------------------------------------------------------------------
    ! Mixture Implementation
    !-----------------------------------------------------------------------
    recursive function mixture_init(thermo, species_names, element_names, reactant_names, &
        input_reactants, omitted_product_names, sort_condensed, ions) result(self)
        ! Create a Mixture from a subset of species in a ThermoDB
        ! If element_names is specified, it must be a superset of elements in species_names
        ! This allows constructing multiple Mixtures with consistent element lists

        ! TODO: Should sorting condensed be default? I think yes.
        !       As-is, user of Fortran API likely to pass bad order to solver

        ! Arguments
        type(ThermoDB), intent(in) :: thermo
        character(*), intent(in), optional :: species_names(:)
        character(*), intent(in), optional :: element_names(:)
        character(*), intent(in), optional :: reactant_names(:)
        type(ReactantInput), intent(in), optional :: input_reactants(:)
        character(*), intent(in), optional :: omitted_product_names(:)
        logical, intent(in), optional :: sort_condensed
        logical, intent(in), optional :: ions

        ! Result
        type(Mixture) :: self

        ! Locals
        logical :: sort_condensed_
        logical, allocatable :: found_db(:)
        integer :: i, j, k, ns
        character(enl), allocatable :: enames(:,:)
        character(:), allocatable :: slist(:), elist(:)
        real(dp) :: h_val, T_val
        type(Mixture) :: reactants

        ! Optional argument handling
        sort_condensed_ = .false.
        if (present(sort_condensed)) sort_condensed_ = sort_condensed
        if (present(ions)) self%ions = ions

        ! Build the species list
        if (present(species_names)) then
            call check_name_list_len(species_names, snl, 'mixture_init species')
            slist = species_names
        else if (present(input_reactants)) then
            block
                character(snl), allocatable :: names(:)
                allocate(names(size(input_reactants)))
                do i = 1, size(input_reactants)
                    call check_name_len(input_reactants(i)%name, snl, 'mixture_init reactant')
                    names(i) = input_reactants(i)%name
                end do
                slist = names
            end block
        else if (present(reactant_names)) then
            reactants = Mixture(thermo, reactant_names, ions=ions)
            slist = reactants%get_products(thermo, omitted_product_names)
        else
            call abort("Must specify either species_names or reactant_names")
        end if

        ! Populate the species data
        ns = size(slist)
        if (ns == 0) then
            call abort('mixture_init: empty species list')
        end if
        allocate(self%species(ns))
        allocate(found_db(ns))
        do i = 1,ns
            self%species(i) = get_species(thermo, slist(i), found_db(i))
            if (.not. found_db(i)) then
                self%species(i)%name = slist(i)
                self%species(i)%i_phase = 0
                self%species(i)%num_intervals = 0
                self%species(i)%molecular_weight = empty_dp
            end if
        end do
        self%species_names = slist
        if (.not. present(input_reactants)) then
            do i = 1,ns
                if (.not. found_db(i)) then
                    call abort('mixture_init: Species not found in ThermoDB: '//trim(slist(i)))
                end if
            end do
        end if
        self%is_condensed = self%species%is_condensed()
        self%num_condensed = count(self%is_condensed)
        self%num_species = ns
        self%num_gas = ns - self%num_condensed

        ! Re-order species if requested
        if (sort_condensed_ .and. self%num_condensed > 0) then
            self%species = [ &
                pack(self%species, .not. self%is_condensed), &
                pack(self%species,       self%is_condensed)  &
            ]
            self%species_names = [ &
                pack(self%species_names, .not. self%is_condensed), &
                pack(self%species_names,       self%is_condensed)  &
            ]
            self%is_condensed = self%species%is_condensed()
        end if

        ! Add the info from the ReactantInput
        if (present(input_reactants)) then
            ! Loop over each reactant, checking if any additional info is provided
            do i = 1, size(input_reactants)
                if (.not. found_db(i) .and. .not. allocated(input_reactants(i)%formula)) then
                    call abort('mixture_init: Reactant not found in ThermoDB and no formula provided: '// &
                               trim(input_reactants(i)%name))
                end if

                ! Name
                self%species(i)%name = input_reactants(i)%name

                ! Formula
                if (allocated(input_reactants(i)%formula)) then
                    self%species(i)%formula = input_reactants(i)%formula
                end if

                ! Enthalpy
                if (allocated(input_reactants(i)%enthalpy)) then
                    h_val = convert_units_to_si(input_reactants(i)%enthalpy%values(1), input_reactants(i)%enthalpy%units)
                    self%species(i)%enthalpy_ref = h_val
                end if

                ! Temperature
                if (allocated(input_reactants(i)%temperature)) then
                    T_val = convert_units_to_si(input_reactants(i)%temperature%values(1), input_reactants(i)%temperature%units)
                    self%species(i)%T_ref = T_val
                end if

                ! TODO: is this needed or not?
                ! ! Density
                ! if (allocated(input_reactants(i)%density)) then
                !     rho_val = 0.0! convert_to_si(input_reactants(i)%enthalpy%values(1))
                !     self%species(i)%enthalpy_ref = rho_val
                ! end if

                ! Molecular weight
                if (allocated(input_reactants(i)%molecular_weight)) then
                    self%species(i)%molecular_weight = input_reactants(i)%molecular_weight
                else if (allocated(input_reactants(i)%formula)) then
                    self%species(i)%molecular_weight = molecular_weight_from_formula(input_reactants(i)%formula)
                end if
            end do
        end if

        ! Validate formulas before building the element list
        do i = 1, ns
            if (.not. allocated(self%species(i)%formula)) then
                call abort('mixture_init: Missing formula for species '//trim(self%species_names(i)))
            end if
            if (.not. allocated(self%species(i)%formula%elements)) then
                call abort('mixture_init: Missing formula elements for species '//trim(self%species_names(i)))
            end if
            if (.not. allocated(self%species(i)%formula%coefficients)) then
                call abort('mixture_init: Missing formula coefficients for species '//trim(self%species_names(i)))
            end if
            if (size(self%species(i)%formula%elements) /= size(self%species(i)%formula%coefficients)) then
                call abort('mixture_init: Formula element/coeff size mismatch for species '//trim(self%species_names(i)))
            end if
        end do

        ! Construct the element list
        if (present(element_names)) then

            ! If element_names specified, use that
            call check_name_list_len(element_names, enl, 'mixture_init element')
            self%element_names = element_names

            ! Allow flexible naming of electron
            ! TODO: In ThermoDB/Formula, translate 'E' to 'e-'
            k = findloc(self%element_names, 'e-', 1)
            if (k /= 0) self%element_names(k) = 'E'

        else

            ! Construct element_names from the species formulas
            k = 0
            do i = 1,ns
                k = max(k, size(self%species(i)%formula%elements))
            end do
            if (k == 0) then
                call abort('mixture_init: empty formula element list')
            end if
            allocate(enames(ns,k))
            enames = ' '
            do i = 1,ns
                enames(i,1:size(self%species(i)%formula%elements)) = self%species(i)%formula%elements
            end do

            call build_elem_list(enames, elist)
            self%element_names = elist
        end if

        if (size(self%element_names) == 0) then
            call abort('mixture_init: empty element list')
        end if

        ! Move "E" to the end of the element list if present
        if (self%ions) then
            self%element_names = [ &
                pack(self%element_names, .not. self%element_names == 'E'), &
                pack(self%element_names,       self%element_names == 'E')  &
            ]
        end if

        ! If this is an ionized problem, add "E" to the end of the element list if not present
        if (self%ions .and. self%element_names(size(self%element_names)) /= 'E') then
            self%element_names = [self%element_names, "E "]
        end if

        ! Set the number of elements
        self%num_elements = size(self%element_names)

        ! Construct the stoichiometric matrix
        allocate(self%stoich_matrix(ns,self%num_elements))
        self%stoich_matrix = 0.0d0
        do i = 1,ns
            associate (f => self%species(i)%formula)
                do j = 1,size(f%elements)
                    if (is_empty(f%elements(j))) exit ! No more elements in formula
                    k = findloc(self%element_names, f%elements(j), 1)
                    if (k == 0) call abort('mixture_init: Species element not in element_names')
                    self%stoich_matrix(i,k) = f%coefficients(j)
                end do
            end associate
        end do

    end function

    function mixture_elements_from_species(self, n_species) result(n_elements)
        ! Compute element concentrations from species concentrations
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: n_species(:)
        real(dp) :: n_elements(self%num_elements)
        integer :: e
        call check_array_len(size(n_species), self%num_species, 'mixture_elements_from_species n_species')
        do e = 1,self%num_elements
            n_elements(e) = dot_product(self%stoich_matrix(:,e), n_species)
        end do
    end function

    subroutine mixture_calc_thermo(self, thermo, temperature, condensed)
        ! Evaluate all mixture thermodynamic functions at given temperature
        class(Mixture), intent(in) :: self
        type(MixtureThermo), intent(inout) :: thermo
        real(dp), intent(in) :: temperature
        logical, intent(in), optional :: condensed

        ! Locals
        integer :: n, ng, nc
        logical :: condensed_

        ng = self%num_gas
        nc = self%num_condensed
        condensed_ = .true.
        if (present(condensed)) condensed_ = condensed

        if (self%num_species /= thermo%num_species) then
            thermo = MixtureThermo(self%num_species)
        end if

        do n = 1,ng
            ! TODO: Optimize implementation; lots of extraneous temp interval checking.
            ! TODO: Make enthalpy non-dimensional in ThermoFit, etc.
            thermo%cp(n)       = self%species(n)%calc_cp(temperature)
            thermo%cv(n)       = self%species(n)%calc_cv(temperature)
            thermo%enthalpy(n) = self%species(n)%calc_enthalpy(temperature)/temperature
            thermo%entropy(n)  = self%species(n)%calc_entropy(temperature)
            thermo%energy(n)   = self%species(n)%calc_energy(temperature)/temperature
        end do
        if (condensed_) then
            do n = 1,nc
                thermo%cp(ng+n)       = self%species(ng+n)%calc_cp(temperature)
                thermo%cv(ng+n)       = self%species(ng+n)%calc_cv(temperature)
                thermo%enthalpy(ng+n) = self%species(ng+n)%calc_enthalpy(temperature)/temperature
                thermo%entropy(ng+n)  = self%species(ng+n)%calc_entropy(temperature)
                thermo%energy(ng+n)   = self%species(ng+n)%calc_energy(temperature)/temperature
            end do
        end if

    end subroutine

    function mixture_get_products(self, thermo, omit) result(product_names)
        ! Get the list of possible products from the reactants

        ! Arguments
        class(Mixture), intent(in) :: self
        type(ThermoDB), intent(in) :: thermo
        character(*), intent(in), optional :: omit(:)

        ! Return
        character(snl), allocatable :: product_names(:)

        ! Locals
        integer :: n, i, j, idx(1), np
        logical :: is_omitted, is_product
        real(dp), parameter :: tol = 1.d-10

        if (present(omit)) then
            call check_name_list_len(omit, snl, 'mixture_get_products omit')
        end if

        np = thermo%num_products

        allocate(product_names(np))

        ! Get the list of possible products from the elements
        n = 0
        do i = 1, np
            if (.not. allocated(thermo%product_thermo(i)%formula)) then
                call abort('mixture_get_products: missing formula for product '//trim(thermo%product_name_list(i)))
            end if
            if (.not. allocated(thermo%product_thermo(i)%formula%elements)) then
                call abort('mixture_get_products: missing formula elements for product '//trim(thermo%product_name_list(i)))
            end if
            if (.not. allocated(thermo%product_thermo(i)%formula%coefficients)) then
                call abort('mixture_get_products: missing formula coefficients for product '//trim(thermo%product_name_list(i)))
            end if
            if (size(thermo%product_thermo(i)%formula%elements) /= &
                size(thermo%product_thermo(i)%formula%coefficients)) then
                call abort('mixture_get_products: formula element/coeff size mismatch for product '// &
                           trim(thermo%product_name_list(i)))
            end if
            ! Exclude "omit" names
            ! TODO: pop omit name every time it is found to speed up search
            if (present(omit)) then
                call check_name_list_len(omit, snl, 'mixture_get_products omit')
                is_omitted = .false.
                do j = 1, size(omit)
                    if (names_match(thermo%product_name_list(i), omit(j))) then
                        is_omitted = .true.
                        exit
                    end if
                end do
                if (is_omitted) cycle
            end if

            ! Check that all elements in this species formula are in the reactant element list
            is_product = .true.
            do j = 1, size(thermo%product_thermo(i)%formula%elements)
                if (abs(thermo%product_thermo(i)%formula%coefficients(j)) < tol) cycle
                idx = findloc(self%element_names, thermo%product_thermo(i)%formula%elements(j))
                if (idx(1) == 0) is_product = .false.
            end do

            if (is_product) then
                n = n + 1
                product_names(n) = thermo%product_name_list(i)
            end if
        end do

        product_names = product_names(:n)

    end function

    function mixture_weights_from_of(self, oxidant_weights, fuel_weights, of_ratio) result(weights)
        ! Compute mixture weights from oxidant and fuel weights and O/F ratio

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: oxidant_weights(:), fuel_weights(:)
        real(dp), intent(in) :: of_ratio

        ! Return
        real(dp) :: weights(size(oxidant_weights))

        ! Locals
        real(dp) :: ow(size(weights)), fw(size(weights))  ! Normalized weights
        real(dp) :: mw_ox, mw_fu  ! Molecular weight of total oxidant and fuel

        call check_array_len(size(oxidant_weights), self%num_species, 'mixture_weights_from_of oxidant_weights')
        call check_array_len(size(fuel_weights), self%num_species, 'mixture_weights_from_of fuel_weights')

        ! Compute molecular weights of total oxidant and fuel
        mw_ox = dot_product(oxidant_weights, self%species(:)%molecular_weight)
        mw_fu = dot_product(fuel_weights,    self%species(:)%molecular_weight)

        ! oxidant_weights and fuel_weights should sum to 1; normalize if not
        ow = oxidant_weights / sum(oxidant_weights)
        fw = fuel_weights    / sum(fuel_weights)

        weights = (mw_ox + mw_fu) * (fw + of_ratio*ow) / (1.0d0 + of_ratio)

    end function

    function mixture_chem_eq_ratio_to_of_ratio(self, oxidant_weights, fuel_weights, eq_ratio) result(of_ratio)
        ! Compute o/f ratio from oxidant and fuel weights and chemical equivalence ratio "r"

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: oxidant_weights(:), fuel_weights(:)
        real(dp), intent(in) :: eq_ratio  ! "r": Chemical equivalence ratio in terms of valences

        ! Return
        real(dp) :: of_ratio

        ! Locals
        real(dp) :: vm(self%num_elements), vp(self%num_elements)  ! Valences (m = minus, p = plus)
        real(dp) :: b0_ox(self%num_elements), b0_fu(self%num_elements)
        real(dp) :: vm_fu, vp_fu, vm_ox, vp_ox  ! Total fuel and oxidant valences (m = minus, p = plus)

        call check_array_len(size(oxidant_weights), self%num_species, 'mixture_chem_eq_ratio_to_of_ratio oxidant_weights')
        call check_array_len(size(fuel_weights), self%num_species, 'mixture_chem_eq_ratio_to_of_ratio fuel_weights')

        ! Get the species valences
        call self%get_valence(vm, vp)

        ! Get the element amounts of the total fuel and oxidant
        b0_ox = self%element_amounts_from_weights(oxidant_weights)
        b0_fu = self%element_amounts_from_weights(fuel_weights)

        ! Get the total fuel and oxidant valences
        vm_ox = dot_product(vm, b0_ox)
        vp_ox = dot_product(vp, b0_ox)
        vm_fu = dot_product(vm, b0_fu)
        vp_fu = dot_product(vp, b0_fu)

        ! Compute the o/f ratio: derived from RP-1311 Pt. 1 Equation (9.18)
        ! (o/f) =  -(V_fu(+) + rV_fu(-)) / (V_ox(+) + rV_ox(-))
        of_ratio = -(vp_fu + eq_ratio*vm_fu)/(vp_ox + eq_ratio*vm_ox)

    end function

    function mixture_of_ratio_to_chem_eq_ratio(self, oxidant_weights, fuel_weights, of_ratio) result(eq_ratio)
        ! Compute o/f ratio from oxidant and fuel weights and chemical equivalence ratio "r"

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: oxidant_weights(:), fuel_weights(:)
        real(dp), intent(in) :: of_ratio  ! oxidant-to-fuel ratio value

        ! Return
        real(dp) :: eq_ratio  ! "r": Chemical equivalence ratio in terms of valences

        ! Locals
        real(dp) :: vm(self%num_elements), vp(self%num_elements)  ! Valences (m = minus, p = plus)
        real(dp) :: b0_ox(self%num_elements), b0_fu(self%num_elements)
        real(dp) :: vm_fu, vp_fu, vm_ox, vp_ox  ! Total fuel and oxidant valences (m = minus, p = plus)

        call check_array_len(size(oxidant_weights), self%num_species, 'mixture_of_ratio_to_chem_eq_ratio oxidant_weights')
        call check_array_len(size(fuel_weights), self%num_species, 'mixture_of_ratio_to_chem_eq_ratio fuel_weights')

        ! Get the species valences
        call self%get_valence(vm, vp)

        ! Get the element amounts of the total fuel and oxidant
        b0_ox = self%element_amounts_from_weights(oxidant_weights)
        b0_fu = self%element_amounts_from_weights(fuel_weights)

        ! Get the total fuel and oxidant valences
        vm_ox = dot_product(vm, b0_ox)
        vp_ox = dot_product(vp, b0_ox)
        vm_fu = dot_product(vm, b0_fu)
        vp_fu = dot_product(vp, b0_fu)

        ! Compute the eq. ratio: (Eq. 9.18)
        eq_ratio = -(vp_fu + vp_ox*of_ratio)/(vm_fu + vm_ox*of_ratio)

    end function

    function mixture_weight_eq_ratio_to_of_ratio(self, oxidant_weights, fuel_weights, phi) result(of_ratio)
        ! Compute o/f ratio from oxidant and fuel weights and weight equivalence ratio "phi"

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: oxidant_weights(:), fuel_weights(:)
        real(dp) :: b0_ox(self%num_elements), b0_fu(self%num_elements)
        real(dp), intent(in) :: phi  ! Equivalence ratio

        ! Return
        real(dp) :: of_ratio

        ! Locals
        real(dp) :: vm(self%num_elements), vp(self%num_elements)  ! Valences (m = minus, p = plus)
        real(dp) :: vm_fu, vp_fu, vm_ox, vp_ox

        call check_array_len(size(oxidant_weights), self%num_species, 'mixture_weight_eq_ratio_to_of_ratio oxidant_weights')
        call check_array_len(size(fuel_weights), self%num_species, 'mixture_weight_eq_ratio_to_of_ratio fuel_weights')

        ! Get the species valences
        call self%get_valence(vm, vp)

        ! Get the element amounts of the total fuel and oxidant
        b0_ox = self%element_amounts_from_weights(oxidant_weights)
        b0_fu = self%element_amounts_from_weights(fuel_weights)

        ! Get the total fuel and oxidant valences
        vm_ox = dot_product(vm, b0_ox)
        vp_ox = dot_product(vp, b0_ox)
        vm_fu = dot_product(vm, b0_fu)
        vp_fu = dot_product(vp, b0_fu)

        ! Compute the o/f ratio: derived from RP-1311 Pt. 1 Equation (9.18)
        of_ratio = -(vm_fu + vp_fu)/(phi*(vm_ox + vp_ox))

    end function

    function mixture_of_ratio_to_weight_eq_ratio(self, oxidant_weights, fuel_weights, of_ratio) result(phi)
        ! Compute weight equivalence ratio "phi" from o/f ratio + oxidant and fuel weights

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: oxidant_weights(:), fuel_weights(:)
        real(dp) :: b0_ox(self%num_elements), b0_fu(self%num_elements)
        real(dp), intent(in) :: of_ratio

        ! Return
        real(dp) :: phi  ! Equivalence ratio

        ! Locals
        real(dp) :: vm(self%num_elements), vp(self%num_elements)  ! Valences (m = minus, p = plus)
        real(dp) :: vm_fu, vp_fu, vm_ox, vp_ox

        call check_array_len(size(oxidant_weights), self%num_species, 'mixture_of_ratio_to_weight_eq_ratio oxidant_weights')
        call check_array_len(size(fuel_weights), self%num_species, 'mixture_of_ratio_to_weight_eq_ratio fuel_weights')

        ! Get the species valences
        call self%get_valence(vm, vp)

        ! Get the element amounts of the total fuel and oxidant
        b0_ox = self%element_amounts_from_weights(oxidant_weights)
        b0_fu = self%element_amounts_from_weights(fuel_weights)

        ! Get the total fuel and oxidant valences
        vm_ox = dot_product(vm, b0_ox)
        vp_ox = dot_product(vp, b0_ox)
        vm_fu = dot_product(vm, b0_fu)
        vp_fu = dot_product(vp, b0_fu)

        ! Compute the o/f ratio (Eq. 9.18)
        phi = -(vm_fu + vp_fu)/(of_ratio*(vm_ox + vp_ox))

    end function

    function mixture_weights_to_moles(self, weights) result(moles)
        ! Convert weight fractions to mole fractions

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)

        ! Return
        real(dp) :: moles(size(weights))

        ! Locals
        integer :: i
        real(dp) :: mw  ! Molecular weights

        call check_array_len(size(weights), self%num_species, 'mixture_weights_to_moles weights')

        ! Get the species molecular weights
        do i = 1, self%num_species
            mw = self%species(i)%molecular_weight
            moles(i) = weights(i) / mw
        end do

    end function

    function mixture_weights_from_moles(self, moles) result(weights)
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: moles(:)
        real(dp) :: weights(self%num_species)
        call check_array_len(size(moles), self%num_species, 'mixture_weights_from_moles moles')
        weights = moles * self%species%molecular_weight
    end function

    function mixture_moles_from_weights(self, weights) result(moles)
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp) :: moles(self%num_species)
        call check_array_len(size(weights), self%num_species, 'mixture_moles_from_weights weights')
        moles = weights / self%species%molecular_weight
    end function

    function mixture_per_weight_from_per_mole(self, per_mole) result(per_weight)
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: per_mole(:)
        real(dp) :: per_weight(self%num_species)
        call check_array_len(size(per_mole), self%num_species, 'mixture_per_weight_from_per_mole per_mole')
        per_weight = per_mole / self%species%molecular_weight
    end function

    function mixture_per_mole_from_per_weight(self, per_weight) result(per_mole)
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: per_weight(:)
        real(dp) :: per_mole(self%num_species)
        call check_array_len(size(per_weight), self%num_species, 'mixture_per_mole_from_per_weight per_weight')
        per_mole = per_weight * self%species%molecular_weight
    end function

    subroutine mixture_get_valence(self, vm, vp)
        ! Get the valences of the species in the mixture

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(out) :: vm(:), vp(:)  ! Valences (m = minus, p = plus)

        ! Locals
        integer :: i, ne
        real(dp) :: v

        ! Shorthand
        ne = self%num_elements

        ! Loop over the species
        vm = 0.0d0; vp = 0.0d0
        do i = 1, ne
            v = get_atom_valence(self%element_names(i))
            if (v < 0.0d0) then
                vm(i) = vm(i) + v
            else
                vp(i) = vp(i) + v
            end if
        end do

    end subroutine

    function mixture_calc_enthalpy_single(self, weights, temperature) result(h)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: temperature

        ! Return
        real(dp) :: h

        ! Locals
        integer :: j
        real(dp) :: hj, nj

        call check_array_len(size(weights), self%num_species, 'mixture_calc_enthalpy_single weights')

        h = 0.0d0
        do j = 1, self%num_species
            nj = weights(j)/self%species(j)%molecular_weight/sum(weights)
            hj = self%species(j)%calc_enthalpy(temperature)
            h = h + nj*hj
        end do
        h = h * gas_constant

    end function

    function mixture_calc_enthalpy_multi(self, weights, temperatures) result(h)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: temperatures(:)

        ! Return
        real(dp) :: h

        ! Locals
        integer :: j
        real(dp) :: hj, nj

        call check_array_len(size(weights), self%num_species, 'mixture_calc_enthalpy_multi weights')
        call check_array_len(size(temperatures), self%num_species, 'mixture_calc_enthalpy_multi temperatures')

        h = 0.0d0
        do j = 1, self%num_species
            nj = weights(j)/self%species(j)%molecular_weight/sum(weights)
            hj = self%species(j)%calc_enthalpy(temperatures(j))
            h = h + nj*hj
        end do
        h = h * gas_constant

    end function

    function mixture_calc_energy_single(self, weights, temperature) result(e)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: temperature

        ! Return
        real(dp) :: e

        call check_array_len(size(weights), self%num_species, 'mixture_calc_energy_single weights')

        e = self%calc_enthalpy(weights, temperature) &
          - self%calc_pressure(weights, temperature)/sum(weights) ! PV since using weights, not densities.

    end function

    function mixture_calc_energy_multi(self, weights, temperatures) result(e)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: temperatures(:)

        ! Return
        real(dp) :: e

        call check_array_len(size(weights), self%num_species, 'mixture_calc_energy_multi weights')
        call check_array_len(size(temperatures), self%num_species, 'mixture_calc_energy_multi temperatures')

        e = self%calc_enthalpy(weights, temperatures) &
          - self%calc_pressure(weights, temperatures)/sum(weights) ! PV since using weights, not densities.

    end function

    function mixture_calc_entropy_single(self, weights, temperature, pressure) result(s)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: temperature
        real(dp), intent(in) :: pressure

        ! Return
        real(dp) :: s

        ! Locals
        integer :: j
        real(dp) :: sj, nj, n

        call check_array_len(size(weights), self%num_species, 'mixture_calc_entropy_single weights')

        n = 0.0d0
        s = 0.0d0
        do j = 1, self%num_species
            nj = weights(j)/self%species(j)%molecular_weight/sum(weights)
            if (nj < 1e-35) cycle
            sj = self%species(j)%calc_entropy(temperature)
            if (self%is_condensed(j)) then
                s = s + nj*sj
            else
                n = n + nj
                s = s + nj*(sj - log(nj))
            end if
        end do
        s = s - n*log(pressure/n)
        s = s * gas_constant !/ 1.d3

    end function

    function mixture_calc_entropy_multi(self, weights, temperatures, pressures) result(s)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: temperatures(:)
        real(dp), intent(in) :: pressures(:) ! NOTE: should this be a scalar?

        ! Return
        real(dp) :: s

        ! Locals
        integer :: j
        real(dp) :: sj, nj, n

        call check_array_len(size(weights), self%num_species, 'mixture_calc_entropy_multi weights')
        call check_array_len(size(temperatures), self%num_species, 'mixture_calc_entropy_multi temperatures')
        call check_array_len(size(pressures), self%num_species, 'mixture_calc_entropy_multi pressures')

        n = 0.0d0
        s = 0.0d0
        do j = 1, self%num_species
            nj = weights(j)/self%species(j)%molecular_weight/sum(weights)
            if (nj < 1e-35) cycle
            sj = self%species(j)%calc_entropy(temperatures(j))
            if (self%is_condensed(j)) then
                s = s + nj*sj
            else
                n = n + nj
                s = s + nj*(sj - log(nj*pressures(j)))
            end if
        end do
        s = s + n*log(n)
        s = s * gas_constant / 1.d3

    end function

    function mixture_calc_gibbs_energy_single(self, weights, temperature, pressure) result(g)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: temperature
        real(dp), intent(in) :: pressure

        ! Return
        real(dp) :: g

        call check_array_len(size(weights), self%num_species, 'mixture_calc_gibbs_energy_single weights')

        g = self%calc_enthalpy(weights, temperature) &
          - self%calc_entropy(weights, temperature, pressure) * temperature

    end function

    function mixture_calc_gibbs_energy_multi(self, weights, temperatures, pressures) result(g)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: temperatures(:)
        real(dp), intent(in) :: pressures(:)

        ! Return
        real(dp) :: g

        ! Locals
        integer :: j
        real(dp) :: gj, nj, n, nr, pr

        call check_array_len(size(weights), self%num_species, 'mixture_calc_gibbs_energy_multi weights')
        call check_array_len(size(temperatures), self%num_species, 'mixture_calc_gibbs_energy_multi temperatures')
        call check_array_len(size(pressures), self%num_species, 'mixture_calc_gibbs_energy_multi pressures')

        ! Must pre-compute because cannot factor out of the
        ! sum when consitituent temperatures can vary.
        n = 0.0d0
        do j = 1, self%num_species
            if (self%is_condensed(j)) cycle
            nj = weights(j)/self%species(j)%molecular_weight
            n = n + nj
        end do

        g = 0.0d0
        do j = 1, self%num_species
            nj = weights(j)/self%species(j)%molecular_weight
            gj = self%species(j)%calc_gibbs_energy(temperatures(j))
            if (self%is_condensed(j)) then
                g = g + nj*gj
            else
                nr = nj/n
                pr = pressures(j)/std_pressure
                g = g + nj*(gj + temperatures(j)*log(nr*pr))
            end if
        end do
        g = g * gas_constant

    end function

    function mixture_calc_frozen_cp_single(self, weights, temperature) result(cp)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: temperature

        ! Return
        real(dp) :: cp

        ! Locals
        integer :: j
        real(dp) :: cpj, nj

        call check_array_len(size(weights), self%num_species, 'mixture_calc_frozen_cp_single weights')

        cp = 0.0d0
        do j = 1, self%num_species
            nj  = weights(j)/self%species(j)%molecular_weight/sum(weights)
            cpj = self%species(j)%calc_cp(temperature)
            cp  = cp + nj*cpj
        end do
        cp = cp * gas_constant

    end function

    function mixture_calc_frozen_cp_multi(self, weights, temperatures) result(cp)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: temperatures(:)

        ! Return
        real(dp) :: cp

        ! Locals
        integer :: j
        real(dp) :: cpj, nj

        call check_array_len(size(weights), self%num_species, 'mixture_calc_frozen_cp_multi weights')
        call check_array_len(size(temperatures), self%num_species, 'mixture_calc_frozen_cp_multi temperatures')

        cp = 0.0d0
        do j = 1, self%num_species
            nj  = weights(j)/self%species(j)%molecular_weight/sum(weights)
            cpj = self%species(j)%calc_cp(temperatures(j))
            cp  = cp + nj*cpj
        end do
        cp = cp * gas_constant

    end function

    function mixture_calc_frozen_cv_single(self, weights, temperature) result(cv)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: temperature

        ! Return
        real(dp) :: cv

        ! Locals
        integer :: j
        real(dp) :: cvj, nj

        call check_array_len(size(weights), self%num_species, 'mixture_calc_frozen_cv_single weights')

        cv = 0.0d0
        do j = 1, self%num_species
            nj  = weights(j)/self%species(j)%molecular_weight/sum(weights)
            cvj = self%species(j)%calc_cv(temperature)
            cv  = cv + nj*cvj
        end do
        cv = cv * gas_constant

    end function

    function mixture_calc_frozen_cv_multi(self, weights, temperatures) result(cv)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)
        real(dp), intent(in) :: temperatures(:)

        ! Return
        real(dp) :: cv

        ! Locals
        integer :: j
        real(dp) :: cvj, nj

        call check_array_len(size(weights), self%num_species, 'mixture_calc_frozen_cv_multi weights')
        call check_array_len(size(temperatures), self%num_species, 'mixture_calc_frozen_cv_multi temperatures')

        cv = 0.0d0
        do j = 1, self%num_species
            nj  = weights(j)/self%species(j)%molecular_weight/sum(weights)
            cvj = self%species(j)%calc_cv(temperatures(j))
            cv  = cv + nj*cvj
        end do
        cv = cv * gas_constant

    end function

    function mixture_calc_pressure_single(self, densities, temperature) result(p)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: densities(:)
        real(dp), intent(in) :: temperature

        ! Return
        real(dp) :: p

        ! Locals
        integer :: j
        real(dp) :: n, nj

        call check_array_len(size(densities), self%num_species, 'mixture_calc_pressure_single densities')

        n = 0.0d0
        do j = 1, self%num_species
            if (self%is_condensed(j)) cycle
            nj = densities(j)/self%species(j)%molecular_weight
            n = n + nj
        end do
        p = n * gas_constant * temperature

    end function

    function mixture_calc_pressure_multi(self, densities, temperatures) result(p)

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: densities(:)
        real(dp), intent(in) :: temperatures(:)

        ! Return
        real(dp) :: p

        ! Locals
        integer :: j
        real(dp) :: nT, nj

        call check_array_len(size(densities), self%num_species, 'mixture_calc_pressure_multi densities')
        call check_array_len(size(temperatures), self%num_species, 'mixture_calc_pressure_multi temperatures')

        nT = 0.0d0
        do j = 1, self%num_species
            if (self%is_condensed(j)) cycle
            nj = densities(j)/self%species(j)%molecular_weight
            nT = nT + nj * temperatures(j)
        end do
        p = nT * gas_constant

    end function

    function mixture_element_amounts_from_weights(self, weights) result(b0)
        ! For a given set of weights, compute the amount of each element in the mixture

        ! Arguments
        class(Mixture), intent(in) :: self
        real(dp), intent(in) :: weights(:)

        ! Return
        real(dp) :: b0(self%num_elements)

        ! Locals
        integer :: i
        real(dp) :: moles(size(weights))

        call check_array_len(size(weights), self%num_species, 'mixture_element_amounts_from_weights weights')

        ! Convert weights to moles and normalize by *weights*
        ! TODO: Do we want the weight normalization baked into this function?
        moles = self%moles_from_weights(weights)
        moles = moles / sum(weights)

        do i = 1, self%num_elements
            b0(i) = dot_product(self%stoich_matrix(:, i), moles)
        end do

    end function

    function molecular_weight_from_formula(fm) result(molecular_weight)
        ! Compute the molecular weight of a species based on its chemical formula

        ! Inputs
        type(Formula), intent(in) :: fm

        ! Result
        real(dp) :: molecular_weight

        ! Locals
        integer :: i
        real(dp) :: mw

        molecular_weight = 0.0d0
        do i = 1, size(fm%coefficients)
            mw = get_atom_weight(fm%elements(i))
            molecular_weight = molecular_weight + mw*fm%coefficients(i)
        end do

    end function

    !-----------------------------------------------------------------------
    ! MixtureThermo Implementation
    !-----------------------------------------------------------------------
    function mixturethermo_init(num_species) result(self)
        ! Pre-allocate and empty-initialize a MixtureThermo object.
        ! This makes Mixture%calc_thermo method allocation-free
        integer, intent(in) :: num_species
        type(MixtureThermo) :: self
        self%num_species = num_species
        allocate(self%cp(num_species), source=empty_dp)
        allocate(self%cv(num_species), source=empty_dp)
        allocate(self%enthalpy(num_species), source=empty_dp)
        allocate(self%entropy(num_species), source=empty_dp)
        allocate(self%energy(num_species), source=empty_dp)
    end function

    !-----------------------------------------------------------------------
    ! Helper Functions
    !-----------------------------------------------------------------------
    subroutine check_name_len(name, max_len, context)
        character(*), intent(in) :: name
        integer, intent(in) :: max_len
        character(*), intent(in) :: context

        if (len_trim(name) > max_len) then
            call abort(trim(context)//' name too long: '//trim(name))
        end if
    end subroutine

    subroutine check_name_list_len(names, max_len, context)
        character(*), intent(in) :: names(:)
        integer, intent(in) :: max_len
        character(*), intent(in) :: context
        integer :: i

        do i = 1, size(names)
            call check_name_len(names(i), max_len, context)
        end do
    end subroutine

    subroutine check_array_len(n, expected, context)
        integer, intent(in) :: n
        integer, intent(in) :: expected
        character(*), intent(in) :: context

        if (n /= expected) then
            call abort(trim(context)//' size mismatch')
        end if
    end subroutine

    function get_species(thermo, name, found) result(species)
        ! Search a ThermoDB for a species of a given name
        ! Linear search for now. Lots of opportunities to be smarter
        ! TODO: Make this a method on the ThermoDB

        ! Arguments
        type(ThermoDB), intent(in), target :: thermo
        character(*), intent(in) :: name
        logical, intent(out), optional :: found

        ! Return
        type(SpeciesThermo) :: species

        ! Locals
        type(SpeciesThermo), pointer :: candidate
        integer :: i

        if (present(found)) found = .false.

        ! Search products
        do i = 1,thermo%num_products
            candidate => thermo%product_thermo(i)
            if (names_match(name, candidate%name)) then
                species = candidate
                if (present(found)) found = .true.
                return
            end if
        end do

        ! Reactant search
        do i = 1,thermo%num_reactants
            candidate => thermo%reactant_thermo(i)
            if (names_match(name, candidate%name)) then
                species = candidate
                if (present(found)) found = .true.
                return
            end if
        end do

        ! call log_info('Species '//trim(name)//' not found in ThermoDB')

    end function

    logical function names_match(name, dbname)
        ! Check if species name matches database entry
        ! TODO: Remove name decorations in ThermoDB
        character(*), intent(in) :: name, dbname
        if (dbname(1:1) == '*') then
            names_match = (name == dbname(2:))
        else
            names_match = (name == dbname)
        end if
    end function

end module
