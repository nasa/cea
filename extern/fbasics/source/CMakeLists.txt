set(mod_build_dir   "${CMAKE_CURRENT_BINARY_DIR}/modules")
set(mod_install_dir "${CMAKE_INSTALL_INCLUDEDIR}/fbasics")
configure_file(parameters.f90.in parameters.f90)

# Core library + tests
add_library(fbasics_core
    ${CMAKE_CURRENT_BINARY_DIR}/parameters.f90
    algorithms.f90
    c_binding.f90
    findloc.f90
    geometry.f90
    logging.f90
    math.f90
    string_scanner.f90
    timing.f90
    utils.f90
)

add_library("fbasics::core" ALIAS fbasics_core)
set_target_properties(fbasics_core PROPERTIES
    EXPORT_NAME core
    Fortran_MODULE_DIRECTORY ${mod_build_dir}
)
target_include_directories(fbasics_core PUBLIC
    $<BUILD_INTERFACE:${mod_build_dir}>
    $<INSTALL_INTERFACE:${mod_install_dir}>
)

install(TARGETS fbasics_core
    EXPORT fbasics-config
    DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY "${mod_build_dir}/"  # Trailing slash => copy contents
    DESTINATION ${mod_install_dir})

if(FB_BUILD_TESTING)
    add_pfunit_ctest(fb_core_test
        TEST_SOURCES
            algorithms_test.pf
            c_binding_test.pf
            geometry_test.pf
            math_test.pf
            string_scanner_test.pf
            timing_test.pf
            utils_test.pf
        LINK_LIBRARIES fbasics::core
    )
    set_tests_properties(fb_core_test PROPERTIES
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
endif()
