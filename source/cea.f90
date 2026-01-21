module cea
    ! Top level import for the CEA Fortran Application Programming Interface
    use cea_param, only: real_kind => dp,  &
                         species_name_len, &
                         element_name_len, &
                         version_string,   &
                         version_major,    &
                         version_minor,    &
                         version_patch
    use cea_thermo, only: ThermoDB, read_thermo
    use cea_transport, only: TransportDB, read_transport
    use cea_mixture, only: Mixture
    use cea_equilibrium, only: EqSolver, EqSolution, EqPartials
    use cea_rocket, only: RocketSolver, RocketSolution
    use cea_shock, only: ShockSolver, ShockSolution
    use cea_detonation, only: DetonSolver, DetonSolution
    implicit none
end module
