#!/bin/bash
set -e
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXTERN="${ROOT}/extern"

#=======================================================================
# Arguments
#=======================================================================
JOBS=4
GENERATOR="Unix Makefiles"
BUILD_GFE=true
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            echo "usage: develop.sh [-h] [-j JOBS] [-G GENERATOR] [--skip-gfe] [PRESET]"
            echo
            echo "  Configures and builds CEA for development: downloads "
            echo "  test dependencies, builds in debug mode with warning "
            echo "  flags, etc."
            echo
            echo "  Arguments:"
            echo "    PRESET     Preset from CMakePresets.json to use [def: dev]"
            echo "    JOBS       Number of parallel compilation jobs [def: 4]"
            echo "    GENERATOR  CMake build system generators to use [def: "Unix Makefiles"]"
            echo
            echo "  Options:"
            echo "    skip_gfe   Don't build GFE dependency (must be pre-installed)"
            echo
            exit 0
            ;;
        -j|--jobs)
            JOBS="$2"
            shift
            shift
            ;;
        -G|--generator)
            GENERATOR="$2"
            shift
            shift
            ;;
        --skip-gfe)
            BUILD_GFE=false
            shift
            ;;
        *)
            break
            ;;
    esac
done
PRESET=${1:-dev}


#=======================================================================
# Build configuration
#=======================================================================
BUILD="build-${PRESET}" # Must be consistent with CMakePresets.json
INSTALL="${BUILD}/install"

PRESETS="$(cmake --preset ${PRESET} -N)"
PRESET_CC="$( echo "$PRESETS" | grep C_COMP       | cut -d= -f2 | tr -d \")"
PRESET_CXX="$(echo "$PRESETS" | grep CXX_COMP     | cut -d= -f2 | tr -d \")"
PRESET_FC="$( echo "$PRESETS" | grep Fortran_COMP | cut -d= -f2 | tr -d \")"
CC=${CC:-$PRESET_CC}
CXX=${CXX:-$PRESET_CXX}
FC=${FC:-$PRESET_FC}

echo "CMake preset:     ${PRESET}"
echo "CMake generator:  ${GENERATOR}"
echo "Build directory:  ${BUILD}"
echo "Build jobs:       ${JOBS}"
echo "Install prefix:   ${INSTALL}"
echo "C Compiler:       ${CC}"
echo "C++ Compiler:     ${CXX}"
echo "Fortran Compiler: ${FC}"
echo


#=======================================================================
# Check of prior build
#=======================================================================
if [[ -d "${BUILD}" ]]; then
    echo "Previous build exists! Action? [abort,delete,continue]"
    read -e response
    case ${response} in
        abort)
            echo "Aborting"
            exit
            ;;
        delete)
            rm -rf "${BUILD}"
            ;;
        continue)
            ;;
        *)
            echo "Invalid input. Aborting."
            exit 1
    esac
fi


#=======================================================================
# Development Dependency: Goddard Fortran Ecosystem
#=======================================================================
if [ "$BUILD_GFE" = true ]; then
    GFE="${EXTERN}/gfe"
    GFE_BUILD="${BUILD}/extern/gfe"
    if [[ ! -d "${GFE}" ]]; then
        git clone https://github.com/Goddard-Fortran-Ecosystem/GFE "${GFE}"
        (cd "${GFE}" && git submodule update --init)
    fi
    GFE_CMAKE_ARGS=(
        -DCMAKE_CXX_COMPILER="${CXX}"
        -DCMAKE_Fortran_COMPILER="${FC}"
        -DCMAKE_INSTALL_PREFIX="${INSTALL}"
        -DCMAKE_BUILD_TYPE=Release
        -DSKIP_MPI=YES
        -DSKIP_OPENMP=YES
        -DSKIP_FHAMCREST=YES
        -DSKIP_ESMF=YES
        -DSKIP_ROBUST=YES
    )
    case "$(uname -s)" in
        MINGW*|MSYS*|CYGWIN*)
            # gfortran's -cpp preprocessor does not predefine _WIN32 the way a
            # C/C++ compiler does, but pFUnit's FUnit.F90 guards its
            # `use pf_RegexFilter` on `#ifndef _WIN32` while pFUnit's own
            # CMakeLists excludes RegexFilter.F90 from the build via CMake's
            # own (always-correct) WIN32 variable. Without this define,
            # FUnit.F90 fails to compile: it tries to use a module that was
            # never built.
            GFE_CMAKE_ARGS+=(-DCMAKE_Fortran_FLAGS=-D_WIN32)
            ;;
    esac
    cmake -S "${GFE}" -B "${GFE_BUILD}" -G "${GENERATOR}" "${GFE_CMAKE_ARGS[@]}"
    cmake --build ${GFE_BUILD} --target install -j ${JOBS}
fi


#=======================================================================
# CEA
#=======================================================================
cmake --preset="${PRESET}" -G "${GENERATOR}" -DCMAKE_PREFIX_PATH="${INSTALL}"
cmake --build "${BUILD}" -j ${JOBS}  #--target install