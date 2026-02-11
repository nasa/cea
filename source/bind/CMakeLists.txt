option(CEA_ENABLE_BIND_C      "Enable C language bindings"       ON)
option(CEA_ENABLE_BIND_CXX    "Enable C++ language bindings"     OFF)
option(CEA_ENABLE_BIND_PYTHON "Enable Python language bindings"  ON)
option(CEA_ENABLE_BIND_MATLAB "Enable MatLab language bindings"  ON)
option(CEA_ENABLE_BIND_EXCEL  "Enable Excel language bindings"   ON)

if (CEA_ENABLE_BIND_MATLAB)
    set(CEA_ENABLE_BIND_PYTHON ON CACHE BOOL "Enable Python language bindings" FORCE)
endif()

if (CEA_ENABLE_BIND_CXX OR
    CEA_ENABLE_BIND_PYTHON OR
    CEA_ENABLE_BIND_MATLAB OR
    CEA_ENABLE_BIND_EXCEL)
    set(CEA_ENABLE_BIND_C ON CACHE BOOL "Enable C language bindings" FORCE)
endif()

if (CEA_ENABLE_BIND_C)
    enable_language(C)
    add_subdirectory(c)
endif()

#=====================================================================
# Python Bindings
#=====================================================================
if (CEA_ENABLE_BIND_PYTHON)
    set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)

    find_package(PythonExtensions REQUIRED)
    find_package(NumPy REQUIRED)
    find_package(Cython REQUIRED)

    execute_process(
        COMMAND "${PYTHON_EXECUTABLE}"
        -c "import numpy; print(numpy.get_include())"
        OUTPUT_VARIABLE NumPy_INCLUDE_DIRS
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    include_directories(python)

    add_cython_target(libcea libcea.pyx)
    add_library(libcea MODULE ${libcea})
    python_extension_module(libcea)
    if(APPLE)
        # Ensure Python symbols are resolved by the interpreter at load time.
        target_link_options(libcea PRIVATE "-undefined" "dynamic_lookup")
    endif()

    # Ensure proper linking on Windows
    if (WIN32)
        set_target_properties(libcea PROPERTIES
            WINDOWS_EXPORT_ALL_SYMBOLS ON
            PREFIX ""
        )
        target_link_libraries(libcea ${PYTHON_LIBRARIES})
        # numpy headers include a workaround for intel compilers, at least as of 2.4.1
        # They are triggered if __INTEL_COMPILER is defined, but that isn't defined by
        # default if using ifc
        if (CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM" OR CMAKE_C_COMPILER_ID STREQUAL "Intel")
            target_compile_options(libcea PUBLIC "-D__INTEL_COMPILER")
            target_link_options(libcea PUBLIC "/static")
        endif()
    endif()

    # set rpath for dynmaic linkers on non-windows platforms
    if(APPLE)
        set_target_properties(libcea PROPERTIES INSTALL_RPATH @loader_path)
    elseif(LINUX)
        set_target_properties(libcea PROPERTIES INSTALL_RPATH $ORIGIN)
    endif()

    target_link_libraries(libcea cea::bindc)
    target_include_directories(libcea PUBLIC
        ${PYTHON_INCLUDE_DIRS}
        ${NumPy_INCLUDE_DIRS}
        c)
    if(TARGET compile_databases)
        add_dependencies(libcea compile_databases)
    endif()
    # the compiled library
    install(TARGETS libcea DESTINATION lib)
    install(TARGETS cea_bindc DESTINATION lib)
    # default data files (compiled during build)
    install(FILES ${CMAKE_BINARY_DIR}/thermo.lib DESTINATION data)
    install(FILES ${CMAKE_BINARY_DIR}/trans.lib DESTINATION data)

    # rest of the pure python package files handled by pip & scikit-build-core

endif()

#=====================================================================
# Matlab Bindings
#=====================================================================
# if (CEA_ENABLE_BIND_MATLAB)

#     find_package(PythonExtensions REQUIRED)
#     find_package(NumPy REQUIRED)
#     find_package(Cython REQUIRED)

#     execute_process(
#         COMMAND "${PYTHON_EXECUTABLE}"
#         -c "import numpy; print(numpy.get_include())"
#         OUTPUT_VARIABLE NumPy_INCLUDE_DIRS
#         OUTPUT_STRIP_TRAILING_WHITESPACE
#     )

#     include_directories(matlab)

#     add_cython_target(cea_matlab cea_matlab.pyx)
#     add_library(cea_matlab MODULE ${cea_matlab})
#     python_extension_module(cea_matlab)

#     # Ensure proper linking on Windows
#     if (WIN32)
#         set_target_properties(cea_matlab PROPERTIES
#             WINDOWS_EXPORT_ALL_SYMBOLS ON
#             PREFIX ""
#         )
#         target_link_libraries(cea_matlab ${PYTHON_LIBRARIES})
#     endif()

#     target_link_libraries(cea_matlab cea::bindc)
#     target_include_directories(cea_matlab PUBLIC
#         ${PYTHON_INCLUDE_DIRS}
#         ${NumPy_INCLUDE_DIRS}
#         c)
#     install(TARGETS cea_matlab DESTINATION cea_matlab)

# endif()
