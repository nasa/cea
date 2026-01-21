module cea_transport
    !! Transport properties module

    use cea_param, only: dp, empty_dp, stdout, &
                         sn => species_name_len
    use cea_fits, only: TransportFit
    use cea_mixture, only: Mixture
    use fb_findloc, only: findloc
    use fb_logging
    implicit none

    type :: SpeciesTransport
        !! Transport coeffient data for individual species

        ! Required items
        character(16) :: name
            !! Species name
        integer :: num_intervals
            !! Number of temperature intervals

        ! Coefficients
        ! NOTE: This assumes temperature intervals are the same for  eta and lambda, which should be true
        real(dp), allocatable :: T_fit(:, :)
            !! Temperature intervals
        type(TransportFit), allocatable :: eta_coef(:)
            !! Viscosity curve fit coefficients
        type(TransportFit), allocatable :: lambda_coef(:)
            !! Thermal conductivity curve fit coefficients
    contains
        procedure :: calc_eta => st_calc_eta
        procedure :: calc_lambda => st_calc_lambda
    end type

    type :: BinaryTransport
        !! Transport coeffient data for binary interactions

        ! Required items
        character(16) :: name(2)
            !! Binary pair names
        integer :: num_intervals
            !! Number of temperature intervals

        ! Coefficients
        ! NOTE: This assumes temperature intervals are the same for  eta and lambda, which should be true
        real(dp), allocatable :: T_fit(:, :)
            !! Temperature intervals
        type(TransportFit), allocatable :: eta_coef(:)
            !! Viscosity curve fit coefficients
    contains
        procedure :: calc_eta => bt_calc_eta
    end type

    type :: TransportDB
        !! Transport database

        integer :: num_pure
            !! Number of pure species
        integer :: num_binary
            !! Number of binary interactions

        character(16), allocatable :: pure_species(:)
            !! Pure species names
        character(16), allocatable :: binary_species(:, :)
            !! Binary interaction pairs

        type(SpeciesTransport), allocatable :: pure_transport(:)
            !! Transport curve fit data for pure species
        type(BinaryTransport), allocatable :: binary_transport(:)
            !! Transport curve fit data for binary interactions

    end type

contains

    ! function transport_init(num_records) result(tdb)
    !     integer, intent(in) :: num_records
    !     type(TransportDB) :: tdb
    ! end function

    elemental function st_calc_eta(self, T) result(eta)
        class(SpeciesTransport), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: eta
        integer :: i, idx

        ! Select temperature range
        idx = 1
        do i = 1,self%num_intervals
            if (T > self%T_fit(i, 1)) then
                idx = i
            end if
        end do

        ! Evaluate selected fit
        eta = self%eta_coef(idx)%calc_transport_value(T)
    end function

    elemental function st_calc_lambda(self, T) result(lambda)
        class(SpeciesTransport), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: lambda
        integer :: i, idx

        ! Select temperature range
        idx = 1
        do i = 1,self%num_intervals
            if (T > self%T_fit(i, 1)) then
                idx = i
            end if
        end do

        ! Evaluate selected fit
        lambda = self%lambda_coef(idx)%calc_transport_value(T)
    end function

    elemental function bt_calc_eta(self, T) result(eta)
        class(BinaryTransport), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: eta
        integer :: i, idx

        ! Select temperature range
        idx = 1
        do i = 1,self%num_intervals
            if (T > self%T_fit(i, 1)) then
                idx = i
            end if
        end do

        ! Evaluate selected fit
        eta = self%eta_coef(idx)%calc_transport_value(T)
    end function

    function read_transport(filename) result(db)

        ! Inputs
        character(*), intent(in) :: filename

        ! Return
        type(TransportDB) :: db

        ! Locals
        integer :: num_records        ! Number of records
        integer :: fin                ! File unit
        integer :: i, j               ! Index
        character(16) :: species(2)   ! Species names
        real(dp) :: tr_data(6, 3, 2)  ! Transport data
        integer :: n_int              ! Temporary variable for number of intervals


        call log_info('Reading transport database file: '//filename)
        open(newunit=fin, file=filename, &
             status="old", action="read", form="unformatted")

        ! Size input arrays
        read(fin) num_records

        ! Oversize the arrays
        allocate(db%pure_species(num_records), db%pure_transport(num_records), &
                 db%binary_species(num_records, 2), db%binary_transport(num_records))

        db%num_pure = 0
        db%num_binary = 0

        ! Load transport data
        do i = 1, num_records
            read(fin) species(:), tr_data

            if (len_trim(species(2)) == 0) then
                ! Pure species
                db%num_pure = db%num_pure + 1
                db%pure_species(db%num_pure) = species(1)

                ! Oversize the curvefit arrays
                allocate(db%pure_transport(db%num_pure)%T_fit(3, 2), &
                         db%pure_transport(db%num_pure)%eta_coef(3), &
                         db%pure_transport(db%num_pure)%lambda_coef(3))

                db%pure_transport(db%num_pure)%name = species(1)

                ! Get the number of valid intervals for these coefficients
                n_int = 3
                do j = 1,3
                    if (abs(tr_data(1, j, 1)) < 1.0d-6) then
                        n_int = j-1
                        exit
                    end if
                end do
                db%pure_transport(db%num_pure)%num_intervals = n_int

                ! Set the coefficients
                do j = 1, n_int
                    db%pure_transport(db%num_pure)%T_fit(j, :) = tr_data(:2, j, 1)

                    db%pure_transport(db%num_pure)%eta_coef(j)%A    = tr_data(3, j, 1)
                    db%pure_transport(db%num_pure)%eta_coef(j)%B    = tr_data(4, j, 1)
                    db%pure_transport(db%num_pure)%eta_coef(j)%C    = tr_data(5, j, 1)
                    db%pure_transport(db%num_pure)%eta_coef(j)%D    = tr_data(6, j, 1)

                    db%pure_transport(db%num_pure)%lambda_coef(j)%A = tr_data(3, j, 2)
                    db%pure_transport(db%num_pure)%lambda_coef(j)%B = tr_data(4, j, 2)
                    db%pure_transport(db%num_pure)%lambda_coef(j)%C = tr_data(5, j, 2)
                    db%pure_transport(db%num_pure)%lambda_coef(j)%D = tr_data(6, j, 2)
                end do
                db%pure_transport(db%num_pure)%T_fit = db%pure_transport(db%num_pure)%T_fit(:n_int, :)
                db%pure_transport(db%num_pure)%eta_coef = db%pure_transport(db%num_pure)%eta_coef(:n_int)
                db%pure_transport(db%num_pure)%lambda_coef = db%pure_transport(db%num_pure)%lambda_coef(:n_int)

            else
                ! Binary species
                db%num_binary = db%num_binary + 1
                db%binary_species(db%num_binary, :) = species

                ! Oversize the curvefit arrays
                allocate(db%binary_transport(db%num_binary)%T_fit(3, 2), &
                         db%binary_transport(db%num_binary)%eta_coef(3))

                db%binary_transport(db%num_binary)%name = species

                ! Get the number of valid intervals for these coefficients
                n_int = 3
                do j = 1,3
                    if (abs(tr_data(1, j, 1)) < 1.0d-6) then
                        n_int = j-1
                        exit
                    end if
                end do
                db%binary_transport(db%num_binary)%num_intervals = n_int

                ! Set the coefficients
                do j = 1, n_int
                    db%binary_transport(db%num_binary)%T_fit(j, :) = tr_data(:2, j, 1)

                    db%binary_transport(db%num_binary)%eta_coef(j)%A = tr_data(3, j, 1)
                    db%binary_transport(db%num_binary)%eta_coef(j)%B = tr_data(4, j, 1)
                    db%binary_transport(db%num_binary)%eta_coef(j)%C = tr_data(5, j, 1)
                    db%binary_transport(db%num_binary)%eta_coef(j)%D = tr_data(6, j, 1)
                end do
                db%binary_transport(db%num_binary)%T_fit = db%binary_transport(db%num_binary)%T_fit(:n_int, :)
                db%binary_transport(db%num_binary)%eta_coef = db%binary_transport(db%num_binary)%eta_coef(:n_int)

            end if

        end do

        ! Resize arrays
        db%pure_species = db%pure_species(1:db%num_pure)
        db%binary_species = db%binary_species(1:db%num_binary, :)

        ! Cleanup
        close(fin)

    end function

    function get_mixture_transport(all_transport, products, ions) result(transport_db)
        ! Get the transport database for the subset of species in the mixture
        type(TransportDB), intent(in) :: all_transport
        type(Mixture), intent(in) :: products
        logical, intent(in), optional :: ions
        type(TransportDB) :: transport_db

        ! Locals
        integer :: i             ! Index
        logical :: ions_         ! Flag for ions
        integer, allocatable :: pure_indices(:), binary_indices(:)  ! Indices of species in the transport database
        integer :: check_idx(1)  ! Index of the species in the transport database
        integer :: idx           ! Iterator to build the index array

        ! Oversize the index array intially
        allocate(pure_indices(all_transport%num_pure), &
                 binary_indices(all_transport%num_binary))

        ! Handle the optional ions flag
        ions_ = .false.
        if (present(ions)) ions_ = ions

        ! Find the indices of the pure species in the transport database
        ! This should be the faster sorting direction, I think?
        idx = 0
        do i = 1, products%num_gas
            ! Check if the gas product is in the transport database
            check_idx = findloc(all_transport%pure_species, products%species_names(i))
            if (check_idx(1) > 0) then
                idx = idx + 1
                pure_indices(idx) = check_idx(1)
            end if
        end do
        if (ions_) then
            check_idx = findloc(all_transport%pure_species, 'e-')
            if (check_idx(1) > 0) then
                idx = idx + 1
                pure_indices(idx) = check_idx(1)
            end if
        end if
        pure_indices = pure_indices(:idx)

        ! Find the indices of the binary species in the transport database:
        ! We need to look through the binary species in all_transport and see if
        ! both species are in the mixture
        idx = 0
        do i = 1, all_transport%num_binary
            ! Check if the first binary species is in the transport database
            check_idx = findloc(products%species_names, all_transport%binary_species(i, 1))
            if (check_idx(1) > 0) then
                ! Check if the second binary species is in the transport database
                check_idx = findloc(products%species_names, all_transport%binary_species(i, 2))
                if (check_idx(1) > 0) then
                    idx = idx + 1
                    binary_indices(idx) = i
                end if
            else if (ions_) then
                ! Check if the first binary species is "e-"
                ! *** NOTE: This assumes "e-" is always the first species in the name pair
                check_idx = findloc(products%species_names, 'e-')
                if (check_idx(1) > 0) then
                    ! Check if the second binary species is in the transport database
                    check_idx = findloc(products%species_names, all_transport%binary_species(i, 2))
                    if (check_idx(1) > 0) then
                        idx = idx + 1
                        binary_indices(idx) = i
                    end if
                end if
            end if
        end do
        binary_indices = binary_indices(:idx)

        ! Set the transport database for the mixture
        transport_db%num_pure = size(pure_indices)
        transport_db%num_binary = size(binary_indices, 1)
        transport_db%pure_species = all_transport%pure_species(pure_indices)
        transport_db%binary_species = all_transport%binary_species(binary_indices, :)
        transport_db%pure_transport = all_transport%pure_transport(pure_indices)
        transport_db%binary_transport = all_transport%binary_transport(binary_indices)

    end function

end module
