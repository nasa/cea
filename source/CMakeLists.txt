
configure_file(param.f90.in param.f90)


#=====================================================================
# Core library
#=====================================================================
add_installed_library(cea_core
    ${CMAKE_CURRENT_BINARY_DIR}/param.f90
    atomic_data.f90
    cea.f90
    equilibrium.f90
    fits.f90
    input.f90
    mixture.f90
    database_compile.f90
    rocket.f90
    shock.f90
    detonation.f90
    thermo.f90
    transport.f90
    units.f90
)
# Call function to correctly process compiler flags to check Fortran 
# standard compliance based on compiler / OS
set(CMAKE_Fortran_STANDARD 2008)
set(CMAKE_Fortran_STANDARD_REQUIRED ON)
include(../cmake/CompilerFlags.cmake)
project_enable_fortran_std(cea_core)

target_link_libraries(cea_core PRIVATE fbasics::core)

if(CEA_BUILD_TESTING)
    add_pfunit_ctest(cea_core_test
        TEST_SOURCES
            thermo_test.pf
            transport_test.pf
            fits_test.pf
            input_test.pf
            mixture_test.pf
            equilibrium_test.pf
            rocket_test.pf
            shock_test.pf
            detonation_test.pf
        LINK_LIBRARIES cea::core fbasics::core
    )
    set_tests_properties(cea_core_test PROPERTIES
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
endif()


#=====================================================================
# Main application
#=====================================================================
add_executable(cea main.f90)
project_enable_fortran_std(cea)
target_link_libraries(cea PRIVATE cea::core fbasics::core)
install(TARGETS cea DESTINATION ${CMAKE_INSTALL_BINDIR})

if(WIN32)
    # on windows, several dll either need to be installed on the target
    # machine, or they need to be statically linked. linking seemed easier than
    # finding the path to the various DLLs and installing them.
    # If the intel compilers are used on windows, CMAKE_MSVC_RUNTIME_LIBRARY is set in the root
    # CMakeLists.txt to force the use of static runtime libraries for all targets
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        target_link_options(cea
                            PUBLIC "-libs:static")
    else()
        target_link_options(cea
                            PUBLIC "-static-libgfortran"
                            PUBLIC "-static-libgcc")
    endif()
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" OR CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    # Static link runtime libraries on linux if using Intel compilers
    target_link_options(cea
                        PUBLIC "-static-intel")
endif()

if(CEA_BUILD_TESTING)

    # Primary test: Run all example problems from NASA RP-1311
    add_test(
        NAME cea_main_test
        COMMAND cea -v samples/example1
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )

    # Secondary tests: verify CLI behavior
    add_test(NAME cea_main_test_empty  COMMAND cea   )  # Prints help,  returns 1
    add_test(NAME cea_main_test_help   COMMAND cea -h)  # Prints help,  returns 0
    add_test(NAME cea_main_test_noinp  COMMAND cea -v)  # Prints error, returns 1
    set_tests_properties(
        cea_main_test_empty
        cea_main_test_noinp
      PROPERTIES
        WILL_FAIL TRUE
    )

endif()


#=====================================================================
# Multi-Language Bindings
#=====================================================================
add_subdirectory(bind)
