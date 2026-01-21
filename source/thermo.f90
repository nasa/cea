module cea_thermo
    !! Module containing the thermodynamic curve fit data and functions

    ! ** Should we remove the "*" from names? It means that the high temperature range has data, but
    ! if not it's interpolated for that range anyway

    use cea_param, only: dp, empty_dp, stdout, &
                         gas_constant, &
                         sn => species_name_len, &
                         en => element_name_len
    use cea_input, only: Formula
    use cea_fits, only: ThermoFit
    use fb_algorithms, only: sort, unique
    use fb_utils, only: abort, startswith, to_str, is_empty
    use fb_logging
    implicit none

    ! Sizing paramters
    integer, parameter :: num_coefs = 9  ! Num. coeffs. for a given thermo fit
    integer, parameter :: num_fit_g = 3  ! Num. fits for each gas species
    integer, parameter :: max_elem_per_species = 5  ! Maximum possible elements in a species formula

    type :: SpeciesThermo
        !! Container for individual species thermo data

        ! Required items
        character(sn) :: name
            !! Species name
        type(Formula), allocatable :: formula
            !! Species chemical formula
        integer :: i_phase
            !! Phase of the species (0: gas, 1+: condensed)
        integer :: num_intervals
            !! Number of curve-fit temperature intervals
        real(dp) :: molecular_weight
            !! Species molecular weight

        ! Thermo coefficient fits
        real(dp), allocatable :: T_fit(:, :)
            !! Temperature range for each fit
        type(ThermoFit), allocatable :: fits(:)
            !! Curve fit coefficients

        ! Reactants only
        real(dp) :: enthalpy_ref = 0.0
            !! Assigned enthalpy
        real(dp) :: T_ref = 0.0
            !! Reference temperature for assigned enthalpy

    contains
        procedure :: is_condensed => st_is_condensed
        procedure :: calc_cv => st_calc_cv
        procedure :: calc_cp => st_calc_cp
        procedure :: calc_energy => st_calc_energy
        procedure :: calc_enthalpy => st_calc_enthalpy
        procedure :: calc_entropy => st_calc_entropy
        procedure :: calc_gibbs_energy => st_calc_gibbs_energy
        procedure :: calc_potential => st_calc_potential
    end type

    type :: ThermoDB
        !! Thermodynamic data

        integer :: num_gas
            !! Number of gas species
        integer :: num_condensed
            !! Number of *unique* condensed species
            ! Note: not one per temperature interval, as is the case in CEA2
        integer :: num_reactants
            !! Number of reactants
        integer :: num_products
            !! Total number of product species
        integer :: num_elems
            !! Total number of unique elements

        character(sn), allocatable :: product_name_list(:)
            !! List of product names
        character(sn), allocatable :: reactant_name_list(:)
            !! List of reactant names
        character(en), allocatable :: element_name_list(:)
            !! List of element names
        type(SpeciesThermo), allocatable :: product_thermo(:)
            !! Product species thermo data
        type(SpeciesThermo), allocatable :: reactant_thermo(:)
            !! Reactant species thermo data

    end type
    interface ThermoDB
        module procedure :: tdb_init
    end interface

contains

    function tdb_init(ng, nc, nr, np, ne) result(tdb)
        !! Initialize a ThermoDB object
        integer, intent(in) :: ng, nc, nr, np, ne
        type(ThermoDB) :: tdb
        tdb%num_gas       = ng
        tdb%num_condensed = nc
        tdb%num_reactants = nr
        tdb%num_products  = np
        tdb%num_elems     = ne
    end function

    elemental function st_is_condensed(self) result(tf)
        class(SpeciesThermo), intent(in) :: self
        logical :: tf
        tf = (self%i_phase > 0)
    end function

    elemental function st_calc_cv(self, T) result(cv)
        class(SpeciesThermo), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: cv
        integer :: i, idx

        if (.not. allocated(self%fits)) then
            cv = 0.0d0
            return
        end if

        ! Select temperature range
        idx = 1
        do i = 1,self%num_intervals
            if (T > self%T_fit(i, 1)) then
                idx = i
            end if
        end do

        ! Evaluate selected fit
        cv = self%fits(idx)%calc_cv(T)

    end function

    elemental function st_calc_cp(self, T) result(cp)
        class(SpeciesThermo), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: cp
        integer :: i, idx

        if (.not. allocated(self%fits)) then
            cp = 0.0d0
            return
        end if

        ! Select temperature range
        idx = 1
        do i = 1,self%num_intervals
            if (T > self%T_fit(i, 1)) then
                idx = i
            end if
        end do

        ! Evaluate selected fit
        cp = self%fits(idx)%calc_cp(T)

    end function

    elemental function st_calc_energy(self, T) result(u_R)
        class(SpeciesThermo), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: u_R
        integer :: i, idx

        if (.not. allocated(self%fits)) then
            u_R = self%calc_enthalpy(T) - T
            return
        end if

        ! Select temperature range
        idx = 1
        do i = 1,self%num_intervals
            if (T > self%T_fit(i, 1)) then
                idx = i
            end if
        end do

        ! Evaluate selected fit
        u_R = self%fits(idx)%calc_energy(T, log(T))

    end function

    elemental function st_calc_enthalpy(self, T) result(h_R)
        class(SpeciesThermo), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: h_R
        integer :: i, idx

        if (allocated(self%fits)) then

            ! Select temperature range
            idx = 1
            do i = 1,self%num_intervals
                if (T > self%T_fit(i, 1)) then
                    idx = i
                end if
            end do

            ! Evaluate selected fit
            h_R = self%fits(idx)%calc_enthalpy(T, log(T))

        else
            h_R = 1.d3*self%enthalpy_ref/gas_constant
        end if

    end function

    elemental function st_calc_entropy(self, T) result(s_R)
        class(SpeciesThermo), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: s_R
        integer :: i, idx

        if (.not. allocated(self%fits)) then
            s_R = 0.0d0
            return
        end if

        ! Select temperature range
        idx = 1
        do i = 1,self%num_intervals
            if (T > self%T_fit(i, 1)) then
                idx = i
            end if
        end do

        ! Evaluate selected fit
        s_R = self%fits(idx)%calc_entropy(T, log(T))

    end function

    elemental function st_calc_gibbs_energy(self, T) result(g)
        class(SpeciesThermo), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: g
        integer :: i, idx

        if (.not. allocated(self%fits)) then
            g = self%calc_enthalpy(T)
            return
        end if

        ! Select temperature range
        idx = 1
        do i = 1,self%num_intervals
            if (T > self%T_fit(i, 1)) then
                idx = i
            end if
        end do

        ! Evaluate selected fit
        g = self%fits(idx)%calc_gibbs_energy(T, log(T))

    end function

    elemental function st_calc_potential(self, T, log_nj, P, n) result(mu)

        class(SpeciesThermo), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp), intent(in), optional :: log_nj, P, n
        real(dp) :: mu, h, s

        real(dp) :: log_nj_, P_, n_

        ! Set the provided refrence state values
        log_nj_ = 0.0d0; P_ = 1.0d0; n_ = 1.0d0
        if (present(log_nj)) log_nj_ = log_nj
        if (present(P)) P_ = P
        if (present(n)) n_ = n

        ! TODO: Check to make sure log P/n are > 0
        ! TODO: Warn if species is gas and log(nj)/P/n are not provided

        h = self%calc_enthalpy(T)
        s = self%calc_entropy(T)
        mu = h/T - s

        ! Gas phase
        if (self%i_phase <= 0) then
            mu = mu + log_nj_ + log(P_/n_)
        end if

    end function

    subroutine build_product_list(names, ng, nc, species, name_map)
        !! Build a unique list of product names
        !! The list is sorted by phase (gas then condensed), and then alphabetically

        ! Inputs
        character(*), intent(in) :: names(:)
            !! Initial list of product names (not sorted or de-duped)
        integer, intent(in) :: ng
            !! Number of gas phase species
        integer, intent(in) :: nc
            !! Number of condensed phase species
        character(:), allocatable, intent(out) :: species(:)
            !! Unique list of sorted product names (gas then condensed, then alphabetically within each)
        integer, intent(out), allocatable, optional :: name_map(:)
            !! Index mapping from the input list to the unique list

        ! Locals
        integer :: sort_map(ng+nc), uniq_map(ng+nc)
        integer, allocatable :: gas_map(:), condensed_map(:), &
                                gas_uniq_map(:), condensed_uniq_map(:)
        integer :: num_gas, num_condensed
        character(:), allocatable :: gas_species(:), condensed_species(:)
        call log_info('Building list of unique species')

        ! Group the species by phase, and then sort each group alphabetaically
        allocate(species, source=names)  ! Only used for initial sizing
        allocate(gas_species, source=names(1:ng))
        allocate(condensed_species, source=names(ng+1:ng+nc))

        call sort(gas_species, gas_map)
        if (nc > 0) call sort(condensed_species, condensed_map)

        ! Remove duplicates from each list
        num_condensed = 0
        call unique(gas_species, gas_uniq_map, num_gas)
        if (nc > 0) call unique(condensed_species, condensed_uniq_map, num_condensed)

        ! Combine: species = [gas_species, condensed_species]
        species = species(1:num_gas+num_condensed)  ! Resize
        species(1:num_gas) = gas_species(1:num_gas)
        species(num_gas+1:num_gas+num_condensed) = condensed_species(1:num_condensed)

        call log_info('Found '//to_str(size(species))//' unique species')

        if (present(name_map)) then
            sort_map(1:ng) = gas_map
            if (nc > 0) sort_map(ng+1:ng+nc) = condensed_map + ng

            uniq_map(1:ng) = gas_uniq_map
            if (nc > 0) uniq_map(ng+1:ng+nc) = condensed_uniq_map + num_gas

            allocate(name_map(num_gas+nc))
            name_map = uniq_map
            name_map(sort_map) = name_map
        end if

    end subroutine

    subroutine build_reactant_list(names, species, name_map)
        ! Build a unique list of reactant names

        ! Inputs
        character(*), intent(in) :: names(:)
        character(:), allocatable, intent(out) :: species(:)
        integer, allocatable, intent(out), optional :: name_map(:)

        ! Locals
        integer, allocatable :: sort_map(:), uniq_map(:)
        integer :: num_unique

        ! This should work (and does in ifort) but doesn't in gfortran
        ! species = sorted(names, sort_map)
        ! species = uniqued(species, uniq_map)

        call log_info('Building list of unique species')
        allocate(species, source=names)
        call sort(species, sort_map)
        call unique(species, uniq_map, num_unique)
        species = species(1:num_unique)
        call log_info('Found '//to_str(size(species))//' unique species')

        if (present(name_map)) then
            name_map = uniq_map
            name_map(sort_map) = name_map
        end if

    end subroutine

    subroutine build_elem_list(sym, elem, sym_map)
        ! Build a unique list of elements

        ! Inputs
        character(*), intent(in) :: sym(:,:)
        character(:), intent(out), allocatable :: elem(:)
        integer, allocatable, intent(out), optional :: sym_map(:,:)

        ! Locals
        integer, allocatable :: sort_map(:), uniq_map(:)
        integer :: num_unique, skip, i

        call log_info('Building list of unique elements')
        elem = reshape(sym, [size(sym)])
        call sort(elem, sort_map)
        call unique(elem, uniq_map, num_unique)

        skip = 0
        do i = 1, num_unique
            if (.not. is_empty(elem(i))) exit
            skip = skip + 1
        end do

        elem = elem(1+skip:num_unique)
        call log_info('Found '//to_str(size(elem))//' unique elements')

        if (present(sym_map)) then
            uniq_map(sort_map) = uniq_map - skip
            sym_map = reshape(uniq_map, shape(sym))
        end if

    end subroutine

    subroutine build_formulas(species_thermo, name_map, sym, fno)

        ! Inputs
        type(SpeciesThermo), intent(inout) :: species_thermo(:)
        integer, intent(in) :: name_map(:)
        character(en), intent(in) :: sym(:, :)
        real(dp), intent(in) :: fno(:, :)

        ! Locals
        integer :: i

        do i = 1,size(name_map)

            ! *** Should I allocate this per-specieis based on non-zero size?
            if (.not. allocated(species_thermo(name_map(i))%formula)) then
                allocate(species_thermo(name_map(i))%formula, &
                         species_thermo(name_map(i))%formula%elements(max_elem_per_species), &
                         species_thermo(name_map(i))%formula%coefficients(max_elem_per_species))
            end if

            species_thermo(name_map(i))%formula%elements = sym(i, 1:max_elem_per_species)
            species_thermo(name_map(i))%formula%coefficients = fno(i, 1:max_elem_per_species)

        end do

    end subroutine

    subroutine build_thermo_coeffs(species_thermo, num_species, name_map, tl, thermo, ntl)

        ! Inputs
        type(SpeciesThermo), intent(inout) :: species_thermo(:)
        integer, intent(in) :: num_species
        integer, intent(in) :: name_map(:)
        real(dp), intent(in) :: tl(:, :)
        real(dp), intent(in) :: thermo(:, :, :)
        integer, intent(in) :: ntl(:)

        ! Locals
        integer :: i, j
        integer :: num_set(num_species)  ! Counter for how many condensed phase temperature intervals have been set
        ! TODO: convert this from "Tg" that gets read in
        real(dp), parameter :: Tg(num_fit_g, 2) = reshape([ &
                                                  200.0d0,  1000.0d0, 6000.0d0, &
                                                  1000.0d0, 6000.0d0, 20000.0d0], shape(Tg))
        integer, parameter :: max_intervals = 10

        num_set = 0
        do i = 1,size(name_map)

            ! Skip reactants with no thermo data defined
            if (ntl(i) == 0) cycle

            ! If this is a gas phase, allocate all thermo coefficients
            if (species_thermo(name_map(i))%i_phase <= 0) then

                if (.not. allocated(species_thermo(name_map(i))%fits)) then
                           allocate(species_thermo(name_map(i))%fits(num_fit_g))
                end if

                species_thermo(name_map(i))%T_fit = Tg
                do j = 1,num_fit_g
                    species_thermo(name_map(i))%fits(j)%a1 = thermo(i, 1, j)
                    species_thermo(name_map(i))%fits(j)%a2 = thermo(i, 2, j)
                    species_thermo(name_map(i))%fits(j)%a3 = thermo(i, 3, j)
                    species_thermo(name_map(i))%fits(j)%a4 = thermo(i, 4, j)
                    species_thermo(name_map(i))%fits(j)%a5 = thermo(i, 5, j)
                    species_thermo(name_map(i))%fits(j)%a6 = thermo(i, 6, j)
                    species_thermo(name_map(i))%fits(j)%a7 = thermo(i, 7, j)
                    species_thermo(name_map(i))%fits(j)%b1 = thermo(i, 8, j)
                    species_thermo(name_map(i))%fits(j)%b2 = thermo(i, 9, j)
                end do

            ! If this is a condensed phase, allocate thermo coefficients for each temperature interval
            else
                if (.not. allocated(species_thermo(name_map(i))%T_fit)) then
                           allocate(species_thermo(name_map(i))%T_fit(max_intervals, 2))
                end if

                if (.not. allocated(species_thermo(name_map(i))%fits)) then
                           allocate(species_thermo(name_map(i))%fits(max_intervals))
                end if

                num_set(name_map(i)) = num_set(name_map(i)) + 1
                species_thermo(name_map(i))%T_fit(num_set(name_map(i)), :) = tl(i, :)

                species_thermo(name_map(i))%fits(num_set(name_map(i)))%a1 = thermo(i, 1, 1)
                species_thermo(name_map(i))%fits(num_set(name_map(i)))%a2 = thermo(i, 2, 1)
                species_thermo(name_map(i))%fits(num_set(name_map(i)))%a3 = thermo(i, 3, 1)
                species_thermo(name_map(i))%fits(num_set(name_map(i)))%a4 = thermo(i, 4, 1)
                species_thermo(name_map(i))%fits(num_set(name_map(i)))%a5 = thermo(i, 5, 1)
                species_thermo(name_map(i))%fits(num_set(name_map(i)))%a6 = thermo(i, 6, 1)
                species_thermo(name_map(i))%fits(num_set(name_map(i)))%a7 = thermo(i, 7, 1)
                species_thermo(name_map(i))%fits(num_set(name_map(i)))%b1 = thermo(i, 8, 1)
                species_thermo(name_map(i))%fits(num_set(name_map(i)))%b2 = thermo(i, 9, 1)

            end if

        end do

        do i = 1,size(name_map)

            ! Re-size arrays and assign number of temperature intervals to the SpeciesThermo
            if (species_thermo(name_map(i))%i_phase > 0) then

                species_thermo(name_map(i))%num_intervals = num_set(name_map(i))
                if (num_set(name_map(i)) > 0) then
                    species_thermo(name_map(i))%fits = &
                    species_thermo(name_map(i))%fits(:num_set(name_map(i)))
                    species_thermo(name_map(i))%T_fit = &
                    species_thermo(name_map(i))%T_fit(:num_set(name_map(i)), :)
                end if

            else
                species_thermo(name_map(i))%num_intervals = num_fit_g
            end if

        end do

    end subroutine

    function read_thermo(filename) result(db)

        ! Inputs
        character(*), intent(in) :: filename

        ! Return
        type(ThermoDB) :: db

        ! Locals
        integer :: ng   ! Number of gas records
        integer :: ngc  ! Number of gas+cond records
        integer :: ngcr ! Number of gas+cond+reactant records
        character(sn), allocatable :: name(:)
        character(sn), allocatable :: reac_name(:)
        character(:), allocatable :: species(:)       ! unique version of name
        character(:), allocatable :: reac_species(:)  ! unique version of reac_name
        character(6), allocatable  :: date(:)
        character(en), allocatable :: sym(:, :)
        character(en), allocatable :: reac_sym(:, :)
        integer, allocatable :: ntl(:)
        integer, allocatable :: ifaz(:)
        integer, allocatable :: reac_ntl(:)
        integer, allocatable :: reac_ifaz(:)
        real(dp), allocatable :: fno(:, :)
        real(dp), allocatable :: tl(:, :)
        real(dp), allocatable :: mwt(:)
        real(dp), allocatable :: thermo(:, :, :)
        real(dp), allocatable :: reac_fno(:, :)
        real(dp), allocatable :: reac_tl(:, :)
        real(dp), allocatable :: reac_mwt(:)
        real(dp), allocatable :: enthalpy(:)
        real(dp), allocatable :: reac_thermo(:, :, :)
        real(dp) :: Tg(num_fit_g+1)
        character(:), allocatable :: elements(:)
        integer, allocatable :: name_map(:)
        integer, allocatable :: reac_name_map(:)
        integer, allocatable :: sym_map(:,:)
        integer :: i, j, fin
        real(dp) :: tl_temp(2), thermo_temp(num_coefs, num_fit_g)

        call log_info('Reading thermo database file: '//filename)
        open(newunit=fin, file=filename, &
             status="old", action="read", form="unformatted")

        ! Size input arrays
        read(fin) Tg, ng, ngc, ngcr

        allocate(name(ngc), date(ngcr), ntl(ngc),             &
                 sym(ngc, max_elem_per_species),              &
                 fno(ngc, max_elem_per_species),              &
                 ifaz(ngc), tl(ngc,2), mwt(ngc),              &
                 thermo(ngc, num_coefs, num_fit_g),           &
                 name_map(ngc),                               &
                 reac_name(ngcr-ngc), reac_ntl(ngcr-ngc),     &
                 reac_sym(ngcr-ngc, max_elem_per_species),    &
                 reac_fno(ngcr-ngc, max_elem_per_species),    &
                 reac_ifaz(ngcr-ngc), reac_tl(ngcr-ngc,2),    &
                 reac_mwt(ngcr-ngc),                          &
                 reac_thermo(ngcr-ngc, num_coefs, num_fit_g), &
                 reac_name_map(ngcr-ngc), enthalpy(ngcr-ngc))

        ! Some records don't set all values in thermo, so initialize with
        ! an "empty" value so we can check later if unset data leaked into
        ! the final data structure.
        thermo = empty_dp

        ! Load gaseous products
        call log_info('Reading '//to_str(ng)//' gaseous products')
        do i = 1,ng
            read(fin) name(i), ntl(i), date(i), &
                      (sym(i, j),fno(i, j), j=1,max_elem_per_species), &
                      ifaz(i), tl_temp, mwt(i), thermo_temp
            tl(i, :) = tl_temp
            thermo(i, :, :) = thermo_temp

            ! Remove the '*' ahead of species names with extended temperature
            ! data. This is not needed, and unique to the .lib format
            if (startswith(name(i),'*')) name(i) = name(i)(2:)
        end do

        ! Load condensed products
        call log_info('Reading '//to_str(ngc-ng)//' condensed products')
        do i = ng+1,ngc
            read(fin) name(i), ntl(i), date(i), &
                      (sym(i, j),fno(i, j), j=1,max_elem_per_species), &
                      ifaz(i), tl_temp, mwt(i), thermo_temp(:, 1)
            tl(i, :) = tl_temp
            thermo(i, :, 1) = thermo_temp(:,1)
        end do

        ! Load reactants
        call log_info('Reading '//to_str(ngcr-ngc)//' reactants')
        do i = 1, (ngcr-ngc)
            read(fin) reac_name(i), reac_ntl(i), date(i+ngc), &
                      (reac_sym(i, j),reac_fno(i, j), j=1,max_elem_per_species), &
                      reac_ifaz(i), tl_temp, reac_mwt(i), enthalpy(i)
            reac_tl(i, :) = tl_temp
            if (reac_ntl(i) > 0) then
                read(fin) thermo_temp
                reac_thermo(i, :, :) = thermo_temp
            end if
        end do

        ! Cleanup
        close(fin)

        ! Build the stochiometric coefficent matrix
        call build_product_list(name, ng, ngc-ng, species, name_map)
        call build_reactant_list(reac_name, reac_species, reac_name_map)
        call build_elem_list(sym, elements, sym_map)

        db = ThermoDB(ng, size(species)-ng, size(reac_species), size(species), size(elements))
        db%product_name_list = species
        db%reactant_name_list = reac_species
        db%element_name_list = elements

        ! Initialize SpeciesThermo arrays
        allocate(db%product_thermo(db%num_products))
        allocate(db%reactant_thermo(db%num_reactants))

        ! Set the species name
        do i = 1,size(species)
            db%product_thermo(i)%name = species(i)
        end do

        do i = 1,size(reac_species)
            db%reactant_thermo(i)%name = reac_species(i)
        end do

        ! Assign formulas
        call build_formulas(db%product_thermo, name_map, sym, fno)
        call build_formulas(db%reactant_thermo, reac_name_map, reac_sym, reac_fno)

        ! Assign phase values
        do i = 1,size(name_map)
            db%product_thermo(name_map(i))%i_phase = ifaz(i)
        end do

        do i = 1,size(reac_name_map)
            db%reactant_thermo(reac_name_map(i))%i_phase = reac_ifaz(i)
        end do

        ! Assign molecular weights
        do i = 1,size(name_map)
            db%product_thermo(name_map(i))%molecular_weight = mwt(i)
        end do

        do i = 1,size(reac_name_map)
            db%reactant_thermo(reac_name_map(i))%molecular_weight = reac_mwt(i)
        end do

        ! Assign thermo coefficients
        call build_thermo_coeffs(db%product_thermo, db%num_products, name_map, tl, thermo, ntl)
        call build_thermo_coeffs(db%reactant_thermo, db%num_reactants, reac_name_map, &
                                 reac_tl, reac_thermo, reac_ntl)

        ! Assign reference temperatures and enthalpy values
        do i = 1,size(reac_name_map)
            if (reac_ntl(i) == 0) then
                db%reactant_thermo(reac_name_map(i))%T_ref = reac_tl(i, 1)
                db%reactant_thermo(reac_name_map(i))%enthalpy_ref = enthalpy(i)
            end if
        end do

    end function

end module
