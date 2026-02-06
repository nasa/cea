cmake_minimum_required(VERSION 3.19)

# Force static linking of runtime libraries
# for Intel & MS compilers on windows
cmake_policy(SET CMP0091 NEW)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)
project(CEA
    VERSION 3.0.0
    LANGUAGES Fortran
)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(GNUInstallDirs)
include(cea_utilities)


#-----------------------------------------------------------------------
# Configuration
#-----------------------------------------------------------------------

# Establish if we are part of a bigger project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(PROJECT_IS_TOP_LEVEL TRUE)
else()
    set(PROJECT_IS_TOP_LEVEL FALSE)
endif()

if(PROJECT_IS_TOP_LEVEL)

     # Default Install Path: Local install in current build dir
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install"
            CACHE PATH "CEA installation directory"
            FORCE)
        set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT FALSE)
        message(STATUS "Setting install path to '${CMAKE_INSTALL_PREFIX}' since none specified.")
    endif()

    # Default Build Type: Standard release build
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
            "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
        message(STATUS "Setting build type to 'Release' since none specified.")
    endif()

endif()

# Testing
option(CEA_BUILD_TESTING "Build test harness for CEA" ${PROJECT_IS_TOP_LEVEL})
if(CEA_BUILD_TESTING)
    find_package(PFUNIT REQUIRED)
    include(CTest)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# change some build systems on windows
if(WIN32)
    # find and alternative to 'ar' archive utility
    find_program(AR_EXECUTABLE NAMES ar)
    if(AR_EXECUTABLE)
        message(STATUS "Found archive utility: ${AR_EXECUTABLE}")
        set(CMAKE_AR ${AR_EXECUTABLE} CACHE FILEPATH "Path to a program." FORCE)
    else()
        find_program(AR_EXECUTABLE NAMES gcc-ar)
        if(AR_EXECUTABLE)
            message(STATUS "Found archive utility: ${AR_EXECUTABLE}")
            set(CMAKE_AR ${AR_EXECUTABLE} CACHE FILEPATH "Path to a program." FORCE)
        else()
            message(WARNING "Archive utility not found.")
        endif()
    endif()
endif()


#-----------------------------------------------------------------------
# Targets
#-----------------------------------------------------------------------
add_subdirectory(extern/fbasics)
add_subdirectory(source)

# Build the thermodynamic properties database
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/thermo.lib
    COMMAND $<TARGET_FILE:cea> --compile-thermo ${CMAKE_SOURCE_DIR}/data/thermo.inp
    DEPENDS cea ${CMAKE_SOURCE_DIR}/data/thermo.inp
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Compiling thermodynamic database from thermo.inp"
)

# Build the transport properties database
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/trans.lib
    COMMAND $<TARGET_FILE:cea> --compile-trans ${CMAKE_SOURCE_DIR}/data/trans.inp
    DEPENDS cea ${CMAKE_SOURCE_DIR}/data/trans.inp
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Compiling transport database from trans.inp"
)

# Create a custom target that depends on both database files
add_custom_target(compile_databases ALL
    DEPENDS ${CMAKE_BINARY_DIR}/thermo.lib ${CMAKE_BINARY_DIR}/trans.lib
)

if(PROJECT_IS_TOP_LEVEL)
    add_custom_command(TARGET compile_databases POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_BINARY_DIR}/thermo.lib ${CMAKE_SOURCE_DIR}/data/thermo.lib
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_BINARY_DIR}/trans.lib ${CMAKE_SOURCE_DIR}/data/trans.lib
        COMMAND ${CMAKE_COMMAND} -E make_directory
                ${CMAKE_SOURCE_DIR}/source/bind/python/cea/data
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_BINARY_DIR}/thermo.lib
                ${CMAKE_SOURCE_DIR}/source/bind/python/cea/data/thermo.lib
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                ${CMAKE_BINARY_DIR}/trans.lib
                ${CMAKE_SOURCE_DIR}/source/bind/python/cea/data/trans.lib
        COMMENT "Copying compiled databases to source data directory"
    )
endif()


#-----------------------------------------------------------------------
# Installation
#-----------------------------------------------------------------------
install(
    FILES ${CMAKE_BINARY_DIR}/thermo.lib ${CMAKE_BINARY_DIR}/trans.lib
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/cea
)
install(
    FILES ${CMAKE_BINARY_DIR}/thermo.lib ${CMAKE_BINARY_DIR}/trans.lib
    DESTINATION ${CMAKE_INSTALL_PREFIX}/data
)
install(EXPORT cea-config
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cea
    NAMESPACE cea::
)
