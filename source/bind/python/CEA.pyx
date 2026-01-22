# Import the necessary Cython modules
from libc.stdlib cimport malloc, free
from cpython.bytes cimport PyBytes_AsString
import cython
import ctypes
import array
import warnings
import os
import importlib.resources as importlib_resources
from typing import Optional

# Import numpy
cimport numpy as np
import numpy as np

np.import_array()

# Import the definitions
from cea_def cimport *


cdef class _CString:
    """Keep a UTF-8 encoded representation alive for the duration of a C call."""

    cdef bytes data
    cdef cea_string ptr

    def __cinit__(self, object value, str context):
        cdef bytes encoded
        if isinstance(value, bytes):
            encoded = <bytes>value
        elif isinstance(value, str):
            encoded = (<str>value).encode('utf-8')
        else:
            raise TypeError(f"{context} must be str or bytes, not {type(value).__name__}")
        self.data = encoded
        self.ptr = <cea_string>PyBytes_AsString(self.data)


cdef inline str cea_convert_chars_to_str(const char* s):
    if s == NULL:
        return None
    elif PY_MAJOR_VERSION >= 3:
        return s.decode('utf-8', 'strict')
    else:
        return str(s)


def _err_name(cea_err ierr):
    if ierr == SUCCESS:
        return "CEA_SUCCESS"
    if ierr == CEA_INVALID_FILENAME:
        return "CEA_INVALID_FILENAME"
    if ierr == CEA_INVALID_PROPERTY_TYPE:
        return "CEA_INVALID_PROPERTY_TYPE"
    if ierr == CEA_INVALID_EQUILIBRIUM_TYPE:
        return "CEA_INVALID_EQUILIBRIUM_TYPE"
    if ierr == CEA_INVALID_INDEX:
        return "CEA_INVALID_INDEX"
    if ierr == CEA_INVALID_SIZE:
        return "CEA_INVALID_SIZE"
    if ierr == CEA_NOT_CONVERGED:
        return "CEA_NOT_CONVERGED"
    return f"CEA_ERR_{int(ierr)}"


def _check_ierr(cea_err ierr, str context, bint allow_not_converged=True):
    if ierr == CEA_NOT_CONVERGED and allow_not_converged:
        warnings.warn(f"{context}: { _err_name(ierr) }", RuntimeWarning)
        return
    if ierr != SUCCESS:
        raise RuntimeError(f"{context} failed with {_err_name(ierr)}")

# Alias the error types
SUCCESS                  = CEA_SUCCESS
INVALID_FILENAME         = CEA_INVALID_FILENAME
INVALID_PROPERTY_TYPE    = CEA_INVALID_PROPERTY_TYPE
INVALID_EQUILIBRIUM_TYPE = CEA_INVALID_EQUILIBRIUM_TYPE

# Alias the log levels
LOG_CRITICAL = CEA_LOG_CRITICAL
LOG_ERROR    = CEA_LOG_ERROR
LOG_WARNING  = CEA_LOG_WARNING
LOG_INFO     = CEA_LOG_INFO
LOG_DEBUG    = CEA_LOG_DEBUG
LOG_NONE     = CEA_LOG_NONE

_py_log_level = LOG_WARNING

# Alias the equilibrium problem types
TP = CEA_TP
HP = CEA_HP
SP = CEA_SP
TV = CEA_TV
UV = CEA_UV
SV = CEA_SV

# Alias the equilibrium problem sizes
NUM_REACTANTS = CEA_NUM_REACTANTS
NUM_PRODUCTS  = CEA_NUM_PRODUCTS
NUM_GAS       = CEA_NUM_GAS
NUM_CONDENSED = CEA_NUM_CONDENSED
NUM_ELEMENTS  = CEA_NUM_ELEMENTS
MAX_EQUATIONS = CEA_MAX_EQUATIONS

# Alias the equilibrium property types
TEMPERATURE              = CEA_TEMPERATURE
PRESSURE                 = CEA_PRESSURE
VOLUME                   = CEA_VOLUME
DENSITY                  = CEA_DENSITY
ENTHALPY                 = CEA_ENTHALPY
ENERGY                   = CEA_ENERGY
ENTROPY                  = CEA_ENTROPY
GIBBS_ENERGY             = CEA_GIBBS_ENERGY
GAMMA_S                  = CEA_GAMMA_S
FROZEN_CP                = CEA_FROZEN_CP
FROZEN_CV                = CEA_FROZEN_CV
EQUILIBRIUM_CP           = CEA_EQUILIBRIUM_CP
EQUILIBRIUM_CV           = CEA_EQUILIBRIUM_CV
VISCOSITY                = CEA_VISCOSITY
FROZEN_CONDUCTIVITY      = CEA_FROZEN_CONDUCTIVITY
EQUILIBRIUM_CONDUCTIVITY = CEA_EQUILIBRIUM_CONDUCTIVITY
FROZEN_PRANDTL           = CEA_FROZEN_PRANDTL
EQUILIBRIUM_PRANDTL      = CEA_EQUILIBRIUM_PRANDTL

# Alias the rocket property types
ROCKET_TEMPERATURE              = CEA_ROCKET_TEMPERATURE
ROCKET_PRESSURE                 = CEA_ROCKET_PRESSURE
ROCKET_VOLUME                   = CEA_ROCKET_VOLUME
ROCKET_DENSITY                  = CEA_ROCKET_DENSITY
ROCKET_M                        = CEA_ROCKET_M
ROCKET_MW                       = CEA_ROCKET_MW
ROCKET_ENTHALPY                 = CEA_ROCKET_ENTHALPY
ROCKET_ENERGY                   = CEA_ROCKET_ENERGY
ROCKET_ENTROPY                  = CEA_ROCKET_ENTROPY
ROCKET_GIBBS_ENERGY             = CEA_ROCKET_GIBBS_ENERGY
ROCKET_GAMMA_S                  = CEA_ROCKET_GAMMA_S
ROCKET_FROZEN_CP                = CEA_ROCKET_FROZEN_CP
ROCKET_FROZEN_CV                = CEA_ROCKET_FROZEN_CV
ROCKET_EQUILIBRIUM_CP           = CEA_ROCKET_EQUILIBRIUM_CP
ROCKET_EQUILIBRIUM_CV           = CEA_ROCKET_EQUILIBRIUM_CV
MACH                            = CEA_MACH
SONIC_VELOCITY                  = CEA_SONIC_VELOCITY
AE_AT                           = CEA_AE_AT
C_STAR                          = CEA_C_STAR
COEFFICIENT_OF_THRUST           = CEA_COEFFICIENT_OF_THRUST
ISP                             = CEA_ISP
ISP_VACUUM                      = CEA_ISP_VACUUM
ROCKET_VISCOSITY                = CEA_ROCKET_VISCOSITY
ROCKET_FROZEN_CONDUCTIVITY      = CEA_ROCKET_FROZEN_CONDUCTIVITY
ROCKET_EQUILIBRIUM_CONDUCTIVITY = CEA_ROCKET_EQUILIBRIUM_CONDUCTIVITY
ROCKET_FROZEN_PRANDTL           = CEA_ROCKET_FROZEN_PRANDTL
ROCKET_EQUILIBRIUM_PRANDTL      = CEA_ROCKET_EQUILIBRIUM_PRANDTL

# Alias the shock property types
SHOCK_TEMPERATURE              = CEA_SHOCK_TEMPERATURE
SHOCK_PRESSURE                 = CEA_SHOCK_PRESSURE
SHOCK_VELOCITY                 = CEA_SHOCK_VELOCITY
SHOCK_MACH                     = CEA_SHOCK_MACH
SHOCK_SONIC_VELOCITY           = CEA_SHOCK_SONIC_VELOCITY
SHOCK_RHO12                    = CEA_SHOCK_RHO12
SHOCK_RHO52                    = CEA_SHOCK_RHO52
SHOCK_P21                      = CEA_SHOCK_P21
SHOCK_P52                      = CEA_SHOCK_P52
SHOCK_T21                      = CEA_SHOCK_T21
SHOCK_T52                      = CEA_SHOCK_T52
SHOCK_M21                      = CEA_SHOCK_M21
SHOCK_M52                      = CEA_SHOCK_M52
SHOCK_V2                       = CEA_SHOCK_V2
SHOCK_U5_P_V2                  = CEA_SHOCK_U5_P_V2
SHOCK_VOLUME                   = CEA_SHOCK_VOLUME
SHOCK_DENSITY                  = CEA_SHOCK_DENSITY
SHOCK_M                        = CEA_SHOCK_M
SHOCK_MW                       = CEA_SHOCK_MW
SHOCK_ENTHALPY                 = CEA_SHOCK_ENTHALPY
SHOCK_ENERGY                   = CEA_SHOCK_ENERGY
SHOCK_ENTROPY                  = CEA_SHOCK_ENTROPY
SHOCK_GIBBS_ENERGY             = CEA_SHOCK_GIBBS_ENERGY
SHOCK_GAMMA_S                  = CEA_SHOCK_GAMMA_S
SHOCK_FROZEN_CP                = CEA_SHOCK_FROZEN_CP
SHOCK_FROZEN_CV                = CEA_SHOCK_FROZEN_CV
SHOCK_EQUILIBRIUM_CP           = CEA_SHOCK_EQUILIBRIUM_CP
SHOCK_EQUILIBRIUM_CV           = CEA_SHOCK_EQUILIBRIUM_CV
SHOCK_VISCOSITY                = CEA_SHOCK_VISCOSITY
SHOCK_FROZEN_CONDUCTIVITY      = CEA_SHOCK_FROZEN_CONDUCTIVITY
SHOCK_EQUILIBRIUM_CONDUCTIVITY = CEA_SHOCK_EQUILIBRIUM_CONDUCTIVITY
SHOCK_FROZEN_PRANDTL           = CEA_SHOCK_FROZEN_PRANDTL
SHOCK_EQUILIBRIUM_PRANDTL      = CEA_SHOCK_EQUILIBRIUM_PRANDTL

# Alias the detonation property types
DETONATION_P1                       = CEA_DETONATION_P1
DETONATION_T1                       = CEA_DETONATION_T1
DETONATION_H1                       = CEA_DETONATION_H1
DETONATION_M1                       = CEA_DETONATION_M1
DETONATION_GAMMA1                   = CEA_DETONATION_GAMMA1
DETONATION_V_SONIC1                 = CEA_DETONATION_V_SONIC1
DETONATION_PRESSURE                 = CEA_DETONATION_PRESSURE
DETONATION_TEMPERATURE              = CEA_DETONATION_TEMPERATURE
DETONATION_DENSITY                  = CEA_DETONATION_DENSITY
DETONATION_ENTHALPY                 = CEA_DETONATION_ENTHALPY
DETONATION_ENERGY                   = CEA_DETONATION_ENERGY
DETONATION_GIBBS_ENERGY             = CEA_DETONATION_GIBBS_ENERGY
DETONATION_ENTROPY                  = CEA_DETONATION_ENTROPY
DETONATION_MACH                     = CEA_DETONATION_MACH
DETONATION_VELOCITY                 = CEA_DETONATION_VELOCITY
DETONATION_SONIC_VELOCITY           = CEA_DETONATION_SONIC_VELOCITY
DETONATION_GAMMA                    = CEA_DETONATION_GAMMA
DETONATION_SONIC_VELOCITY           = CEA_DETONATION_SONIC_VELOCITY
DETONATION_P_P1                     = CEA_DETONATION_P_P1
DETONATION_T_T1                     = CEA_DETONATION_T_T1
DETONATION_M_M1                     = CEA_DETONATION_M_M1
DETONATION_RHO_RHO1                 = CEA_DETONATION_RHO_RHO1
DETONATION_FROZEN_CP                = CEA_DETONATION_FROZEN_CP
DETONATION_FROZEN_CV                = CEA_DETONATION_FROZEN_CV
DETONATION_EQUILIBRIUM_CP           = CEA_DETONATION_EQUILIBRIUM_CP
DETONATION_EQUILIBRIUM_CV           = CEA_DETONATION_EQUILIBRIUM_CV
DETONATION_M                        = CEA_DETONATION_M
DETONATION_MW                       = CEA_DETONATION_MW
DETONATION_VISCOSITY                = CEA_DETONATION_VISCOSITY
DETONATION_FROZEN_CONDUCTIVITY      = CEA_DETONATION_FROZEN_CONDUCTIVITY
DETONATION_EQUILIBRIUM_CONDUCTIVITY = CEA_DETONATION_EQUILIBRIUM_CONDUCTIVITY
DETONATION_FROZEN_PRANDTL           = CEA_DETONATION_FROZEN_PRANDTL
DETONATION_EQUILIBRIUM_PRANDTL      = CEA_DETONATION_EQUILIBRIUM_PRANDTL

# Version
def _version():
    """
    Get the complete version string of the CEA library.

    Returns
    -------
    str
        Version string in format "major.minor.patch"
    """
    return "{0:d}.{1:d}.{2:d}".format(_version_major(), _version_minor(), _version_patch())

def _version_major():
    """
    Get the major version number of the CEA library.

    Returns
    -------
    int
        Major version number
    """
    cdef cea_err ierr
    cdef cea_int major
    ierr = cea_version_major(&major)
    return major

def _version_minor():
    """
    Get the minor version number of the CEA library.

    Returns
    -------
    int
        Minor version number
    """
    cdef cea_err ierr
    cdef cea_int minor
    ierr = cea_version_minor(&minor)
    return minor

def _version_patch():
    """
    Get the patch version number of the CEA library.

    Returns
    -------
    int
        Patch version number
    """
    cdef cea_err ierr
    cdef cea_int patch
    ierr = cea_version_patch(&patch)
    return patch

def set_log_level(level):
    """
    Set the logging level for CEA library output.

    Parameters
    ----------
    level : int (enum)
        Log level constant (LOG_CRITICAL, LOG_ERROR, LOG_WARNING, LOG_INFO, LOG_DEBUG, LOG_NONE)
    """
    cdef cea_err ierr
    ierr = cea_set_log_level(level)
    global _py_log_level
    _py_log_level = int(level)
    return

def _maybe_print_init_path(label, path):
    if _py_log_level == LOG_NONE:
        return
    if path:
        print(f"Loaded {label} from: {path}")

def init(path=None, trans=True):
    """
    Initialize the CEA library with thermodynamic and transport data files.

    Parameters
    ----------
    path : str, optional
        Path to directory containing thermo.lib and trans.lib files. If None, searches the
        current directory, then CEA_DATA_DIR, then the packaged cea/data directory, then
        the repo data directory if present.
    trans : bool, default True
        Whether to initialize transport properties. Transport data are optional.

    Raises
    ------
    TypeError
        If path is not a string
    ValueError
        If thermo.lib file not found at specified path
    """
    cdef cea_err ierr
    cdef _CString cthermo_path
    cdef _CString ctrans_path
    if path is None:
        candidate_dirs = [os.getcwd()]
        env_dir = os.environ.get("CEA_DATA_DIR", "")
        if env_dir:
            candidate_dirs.append(env_dir)

        dev_data_dir = None

        pkg_data = None
        try:
            pkg_data = importlib_resources.files("cea").joinpath("data")
        except Exception:
            pkg_data = None

        if pkg_data is not None:
            with importlib_resources.as_file(pkg_data) as pkg_path:
                if dev_data_dir is None:
                    try:
                        # Find repo root by locating data/thermo.inp for dev-tree runs.
                        search_root = os.path.dirname(pkg_path)
                        for _ in range(6):
                            dev_candidate = os.path.join(search_root, "data")
                            if os.path.isfile(os.path.join(dev_candidate, "thermo.inp")):
                                dev_data_dir = dev_candidate
                                break
                            search_root = os.path.dirname(search_root)
                    except Exception:
                        dev_data_dir = None
                pkg_dirs = candidate_dirs + [str(pkg_path)]
                if dev_data_dir and dev_data_dir not in pkg_dirs:
                    pkg_dirs.append(dev_data_dir)
                thermo_path = None
                trans_path = None
                for dirname in pkg_dirs:
                    candidate = os.path.join(dirname, "thermo.lib")
                    if os.path.isfile(candidate):
                        thermo_path = candidate
                        break
                for dirname in pkg_dirs:
                    candidate = os.path.join(dirname, "trans.lib")
                    if os.path.isfile(candidate):
                        trans_path = candidate
                        break

                if thermo_path is None:
                    raise ValueError("thermo.lib not found. Searched: " + ", ".join(pkg_dirs))

                cthermo_path = _CString(thermo_path, "init path")
                ierr = cea_init_thermo(cthermo_path.ptr)
                if ierr != SUCCESS:
                    raise ValueError("Invalid path; thermo.lib not found.")
                _maybe_print_init_path("thermo.lib", thermo_path)

                if trans:
                    if trans_path is None:
                        warnings.warn(
                            "trans.lib not found; continuing without transport properties.",
                            RuntimeWarning,
                        )
                    else:
                        ctrans_path = _CString(trans_path, "init path")
                        ierr = cea_init_trans(ctrans_path.ptr)
                        if ierr != SUCCESS:
                            warnings.warn(
                                "trans.lib could not be initialized; continuing without transport properties.",
                                RuntimeWarning,
                            )
                        else:
                            _maybe_print_init_path("trans.lib", trans_path)
                return

        if dev_data_dir and dev_data_dir not in candidate_dirs:
            candidate_dirs.append(dev_data_dir)

        thermo_path = None
        trans_path = None
        for dirname in candidate_dirs:
            candidate = os.path.join(dirname, "thermo.lib")
            if os.path.isfile(candidate):
                thermo_path = candidate
                break
        for dirname in candidate_dirs:
            candidate = os.path.join(dirname, "trans.lib")
            if os.path.isfile(candidate):
                trans_path = candidate
                break

        if thermo_path is None:
            raise ValueError("thermo.lib not found. Searched: " + ", ".join(candidate_dirs))

        cthermo_path = _CString(thermo_path, "init path")
        ierr = cea_init_thermo(cthermo_path.ptr)
        if ierr != SUCCESS:
            raise ValueError("Invalid path; thermo.lib not found.")
        _maybe_print_init_path("thermo.lib", thermo_path)

        if trans:
            if trans_path is None:
                warnings.warn(
                    "trans.lib not found; continuing without transport properties.",
                    RuntimeWarning,
                )
            else:
                ctrans_path = _CString(trans_path, "init path")
                ierr = cea_init_trans(ctrans_path.ptr)
                if ierr != SUCCESS:
                    warnings.warn(
                        "trans.lib could not be initialized; continuing without transport properties.",
                        RuntimeWarning,
                    )
                else:
                    _maybe_print_init_path("trans.lib", trans_path)
        return

    if type(path) is not str:
        raise TypeError("init path must be a string to the location of thermo.lib, and optionally trans.lib")
    thermo_path = os.path.join(path, "thermo.lib")
    trans_path = os.path.join(path, "trans.lib")
    cthermo_path = _CString(thermo_path, "init path")
    ierr = cea_init_thermo(cthermo_path.ptr)
    if ierr != SUCCESS:
        raise ValueError("Invalid path; thermo.lib not found.")
    _maybe_print_init_path("thermo.lib", thermo_path)
    if trans:
        ctrans_path = _CString(trans_path, "init path")
        ierr = cea_init_trans(ctrans_path.ptr)
        if ierr != SUCCESS:
            warnings.warn(
                "trans.lib not found; continuing without transport properties.",
                RuntimeWarning,
            )
        else:
            _maybe_print_init_path("trans.lib", trans_path)
    return

def init_thermo(thermofile):
    """
    Initialize CEA with thermodynamic data from specified file.

    Parameters
    ----------
    thermofile : str
        Path to thermodynamic data file (typically thermo.lib)

    Raises
    ------
    ValueError
        If thermofile path is invalid or file not found
    """
    cdef cea_err ierr
    cdef _CString cthermo
    cthermo = _CString(thermofile, "init_thermo thermofile")
    ierr = cea_init_thermo(cthermo.ptr)
    if ierr != SUCCESS:
        raise ValueError("Invalid thermofile name.")
    return

def init_trans(transfile):
    """
    Initialize CEA with transport property data from specified file.

    Parameters
    ----------
    transfile : str
        Path to transport data file (typically trans.lib)

    Raises
    ------
    ValueError
        If transfile path is invalid or file not found
    """
    cdef cea_err ierr
    cdef _CString ctrans
    ctrans = _CString(transfile, "init_trans transfile")
    ierr = cea_init_trans(ctrans.ptr)
    if ierr != SUCCESS:
        raise ValueError("Invalid transfile name.")
    return

def is_initialized():
    """
    Check whether the thermodynamic database has been initialized.

    Returns
    -------
    bool
        True if thermo database is initialized, False otherwise
    """
    cdef cea_err ierr
    cdef cea_int initialized
    ierr = cea_is_initialized(&initialized)
    if ierr != SUCCESS:
        raise RuntimeError("Failed to query initialization status.")
    return bool(initialized)


cdef class Mixture:
    """
    Mixture base class for managing chemical species compositions.

    Represents a mixture of chemical species and provides methods for
    converting between moles and weights, calculating mixture properties,
    and converting equivalence ratios.

    Parameters
    ----------
    species : list of str
        List of chemical species names
    products_from_reactants : bool, default False
        If True, treat species as reactants and generate corresponding products
    omit : list of str, default []
        List of species to omit from products (only used when products_from_reactants=True)
    ions : bool, default False
        If True, include ionized species in the mixture
    """
    cdef cea_mixture ptr
    cdef public int num_species
    cdef object _keepalive_species
    cdef object _keepalive_omit
    def __init__(self, list species=None, bint products_from_reactants=False, list omit=[], bint ions=False):
        """
        Create the Mixture object
        If products_from_reactants is True, the species are treated as reactants, and the
        returned object is the corresponding product Mixture
        """

        cdef cea_err ierr
        cdef cea_int num_products
        self.num_species = <int>len(species)
        cdef cea_string* cea_species = NULL
        cdef int nomit = <int>len(omit)
        cdef cea_string* cea_omit = NULL
        cdef list _species_keepalive = []
        cdef list _omit_keepalive = []
        cdef _CString _encoded

        cea_species = <cea_string*>malloc(sizeof(cea_string)*len(species))
        if cea_species == NULL:
            raise MemoryError("Failed to allocate species pointer buffer")
        cea_omit = <cea_string*>malloc(sizeof(cea_string)*len(omit))
        if cea_omit == NULL:
            free(cea_species)
            raise MemoryError("Failed to allocate omit pointer buffer")

        try:
            for i, val in enumerate(species):
                _encoded = _CString(val, "Mixture species entries")
                _species_keepalive.append(_encoded)
                cea_species[i] = _encoded.ptr
            for i, val in enumerate(omit):
                _encoded = _CString(val, "Mixture omit entries")
                _omit_keepalive.append(_encoded)
                cea_omit[i] = _encoded.ptr
            self._keepalive_species = _species_keepalive
            self._keepalive_omit = _omit_keepalive

            if ions:
                if products_from_reactants:
                    # Create the products mixture from the provided reactants
                    ierr = cea_mixture_create_from_reactants_w_ions(&self.ptr, <cea_int>self.num_species, cea_species, <cea_int>nomit, cea_omit)
                    ierr = cea_mixture_get_num_species(self.ptr, &num_products)
                    self.num_species = num_products
                else:
                    # Create the mixture directly from the provided species (products or reactants)
                    ierr = cea_mixture_create_w_ions(&self.ptr, <cea_int>self.num_species, cea_species)
            else:
                if products_from_reactants:
                    # Create the products mixture from the provided reactants
                    ierr = cea_mixture_create_from_reactants(&self.ptr, <cea_int>self.num_species, cea_species, <cea_int>nomit, cea_omit)
                    ierr = cea_mixture_get_num_species(self.ptr, &num_products)
                    self.num_species = num_products
                else:
                    # Create the mixture directly from the provided species (products or reactants)
                    ierr = cea_mixture_create(&self.ptr, <cea_int>self.num_species, cea_species)
        finally:
            free(cea_species)
            free(cea_omit)

        return

    def __dealloc__(self):
        if self.ptr:
            cea_mixture_destroy(&self.ptr)
        return

    property species_names:
        """
        Get list of species names in the mixture.

        Returns
        -------
        list of str
            Species names
        """
        def __get__(self):
            cdef list species = [None] * self.num_species
            for i in range(self.num_species):
                species[i] = self._get_species_name(i)

            return species

    def _get_species_name(self, int i):
        cdef cea_err ierr
        cdef cea_int name_len
        cdef char *cea_species

        ierr = cea_species_name_len(&name_len)
        if ierr != CEA_SUCCESS:
            raise RuntimeError("Failed to query species name length")

        cea_species = <char *>malloc((name_len + 1) * sizeof(char))
        if cea_species == NULL:
            raise MemoryError("Failed to allocate species name buffer")

        ierr = cea_mixture_get_species_name_buf(&self.ptr, i, cea_species, name_len + 1)
        if ierr != CEA_SUCCESS:
            free(cea_species)
            raise RuntimeError("Failed to fetch species name")

        species = cea_convert_chars_to_str(cea_species)
        free(cea_species)

        return species

    def moles_to_weights(self, np.ndarray moles):
        """
        Convert molar amounts to mass fractions.

        Parameters
        ----------
        moles : np.ndarray
            Molar amounts of each species

        Returns
        -------
        np.ndarray
            Mass fractions of each species
        """
        cdef cea_err ierr
        cdef int nspecies = <int>len(moles)
        cdef cea_real *reac_moles = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if reac_moles == NULL:
            raise MemoryError("Failed to allocate moles buffer")
        cdef cea_real *reac_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if reac_weights == NULL:
            free(reac_moles)
            raise MemoryError("Failed to allocate weights buffer")
        weights = np.zeros(nspecies, dtype=np.double)

        try:
            for i in range(nspecies):
                reac_moles[i] = moles[i]

            ierr = cea_mixture_moles_to_weights(self.ptr, nspecies, reac_moles, reac_weights)
            _check_ierr(ierr, "Mixture.moles_to_weights")

            # Assign the weights to the numpy array
            for i in range(nspecies):
                weights[i] = reac_weights[i]
        finally:
            free(reac_moles)
            free(reac_weights)

        return weights

    def weights_to_moles(self, np.ndarray weights):
        """
        Convert mass fractions to molar amounts.

        Parameters
        ----------
        weights : np.ndarray
            Mass fractions of each species

        Returns
        -------
        np.ndarray
            Molar amounts of each species
        """
        cdef cea_err ierr
        cdef int nspecies = <int>len(weights)
        cdef cea_real *reac_moles = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if reac_moles == NULL:
            raise MemoryError("Failed to allocate moles buffer")
        cdef cea_real *reac_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if reac_weights == NULL:
            free(reac_moles)
            raise MemoryError("Failed to allocate weights buffer")
        moles = np.zeros(nspecies, dtype=np.double)

        try:
            for i in range(nspecies):
                reac_weights[i] = weights[i]

            ierr = cea_mixture_weights_to_moles(self.ptr, nspecies, reac_weights, reac_moles)
            _check_ierr(ierr, "Mixture.weights_to_moles")

            # Assign the weights to the numpy array
            for i in range(nspecies):
                moles[i] = reac_moles[i]
        finally:
            free(reac_moles)
            free(reac_weights)

        return moles

    def chem_eq_ratio_to_of_ratio(self, np.ndarray oxidant_weights, np.ndarray fuel_weights, float chem_eq_ratio):
        """
        Convert chemical equivalence ratio (r,eq) to oxidizer/fuel mass ratio.

        Parameters
        ----------
        oxidant_weights : np.ndarray
            Mass fractions of oxidizer species
        fuel_weights : np.ndarray
            Mass fractions of fuel species
        chem_eq_ratio : float
            Chemical equivalence ratio (r,eq)

        Returns
        -------
        float
            Oxidizer/fuel mass ratio
        """
        # 'r,eq' to o/f ratio
        cdef cea_err ierr
        cdef int nspecies = <int>len(oxidant_weights)
        cdef cea_real of_ratio

        cdef cea_real *ox_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if ox_weights == NULL:
            raise MemoryError("Failed to allocate oxidant weights buffer")
        cdef cea_real *fu_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if fu_weights == NULL:
            free(ox_weights)
            raise MemoryError("Failed to allocate fuel weights buffer")

        try:
            for i in range(nspecies):
                ox_weights[i] = oxidant_weights[i]
                fu_weights[i] = fuel_weights[i]

            ierr = cea_mixture_chem_eq_ratio_to_of_ratio(self.ptr, nspecies, ox_weights, fu_weights, chem_eq_ratio, &of_ratio)
            _check_ierr(ierr, "Mixture.chem_eq_ratio_to_of_ratio")
        finally:
            free(ox_weights)
            free(fu_weights)

        return of_ratio

    def weight_eq_ratio_to_of_ratio(self, np.ndarray oxidant_weights, np.ndarray fuel_weights, float weight_eq_ratio):
        """
        Convert weight equivalence ratio (phi) to oxidizer/fuel mass ratio.

        Parameters
        ----------
        oxidant_weights : np.ndarray
            Mass fractions of oxidizer species
        fuel_weights : np.ndarray
            Mass fractions of fuel species
        weight_eq_ratio : float
            Weight equivalence ratio (phi)

        Returns
        -------
        float
            Oxidizer/fuel mass ratio
        """
        # 'phi' to o/f ratio
        cdef cea_err ierr
        cdef int nspecies = <int>len(oxidant_weights)
        cdef cea_real of_ratio

        cdef cea_real *ox_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if ox_weights == NULL:
            raise MemoryError("Failed to allocate oxidant weights buffer")
        cdef cea_real *fu_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if fu_weights == NULL:
            free(ox_weights)
            raise MemoryError("Failed to allocate fuel weights buffer")

        try:
            for i in range(nspecies):
                ox_weights[i] = oxidant_weights[i]
                fu_weights[i] = fuel_weights[i]

            ierr = cea_mixture_weight_eq_ratio_to_of_ratio(self.ptr, nspecies, ox_weights, fu_weights, weight_eq_ratio, &of_ratio)
            _check_ierr(ierr, "Mixture.weight_eq_ratio_to_of_ratio")
        finally:
            free(ox_weights)
            free(fu_weights)

        return of_ratio

    def of_ratio_to_weights(self, np.ndarray oxidant_weights, np.ndarray fuel_weights, float of_ratio):
        """
        Convert oxidizer/fuel ratio to reactant mass fractions.

        Parameters
        ----------
        oxidant_weights : np.ndarray
            Mass fractions of oxidizer species
        fuel_weights : np.ndarray
            Mass fractions of fuel species
        of_ratio : float
            Oxidizer/fuel mass ratio

        Returns
        -------
        np.ndarray
            Combined reactant mass fractions
        """
        cdef cea_err ierr
        cdef int nspecies = <int>len(oxidant_weights)

        cdef cea_real *ox_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if ox_weights == NULL:
            raise MemoryError("Failed to allocate oxidant weights buffer")
        cdef cea_real *fu_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if fu_weights == NULL:
            free(ox_weights)
            raise MemoryError("Failed to allocate fuel weights buffer")
        cdef cea_real *reac_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if reac_weights == NULL:
            free(ox_weights)
            free(fu_weights)
            raise MemoryError("Failed to allocate reactant weights buffer")
        weights = np.zeros(nspecies, dtype=np.double)

        try:
            for i in range(nspecies):
                ox_weights[i] = oxidant_weights[i]
                fu_weights[i] = fuel_weights[i]

            ierr = cea_mixture_of_ratio_to_weights(self.ptr, nspecies, ox_weights, fu_weights, of_ratio, reac_weights)
            _check_ierr(ierr, "Mixture.of_ratio_to_weights")

            # Assign the weights to the numpy array
            for i in range(nspecies):
                weights[i] = reac_weights[i]
        finally:
            free(ox_weights)
            free(fu_weights)
            free(reac_weights)

        return weights

    def calc_property(self, cea_property_type prop_type, np.ndarray weights,
                      temperature: float | list | np.ndarray, pressure: Optional[float] = None):
        """
        Calculate thermodynamic property for mixture at specified conditions.

        Parameters
        ----------
        prop_type : int
            Property type constant (ENTHALPY, ENERGY, FROZEN_CP, FROZEN_CV)
        weights : np.ndarray
            Mass fractions for each species
        temperature : float, list, or np.ndarray
            Temperature(s) in K. If list/array, must match species length for multi-temperature calculation
        pressure : float, optional
            Pressure in bar (required for enthalpy and Gibbs energy calculations)

        Returns
        -------
        float
            Calculated property value

        Raises
        ------
        ValueError
            If property type not supported or temperature format invalid
        """
        cdef cea_err ierr
        cdef cea_real value
        cdef int nspecies = <int>len(weights)
        cdef cea_real *reac_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if reac_weights == NULL:
            raise MemoryError("Failed to allocate weights buffer")

        for i in range(nspecies):
            reac_weights[i] = weights[i]

        if prop_type not in [VOLUME, DENSITY, ENTHALPY, ENERGY, FROZEN_CP, FROZEN_CV, ENTROPY, GIBBS_ENERGY]:
            raise ValueError("Property type not supported for mixture calculations")

        if prop_type in [ENTROPY, GIBBS_ENERGY, VOLUME, DENSITY]:
            if pressure is None:
                raise ValueError("Pressure must be provided for enthalpy and Gibbs energy calculations")

        # Handle the case where temperature is a list of numpy array
        cdef cea_real *reac_temps = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if reac_temps == NULL:
            free(reac_weights)
            raise MemoryError("Failed to allocate temperature buffer")

        try:
            if prop_type in [ENTROPY, GIBBS_ENERGY, VOLUME, DENSITY]:
                if isinstance(temperature, float):
                    ierr = cea_mixture_calc_property_tp(self.ptr, prop_type, nspecies, reac_weights, temperature, pressure, &value)
                elif type(temperature) in [list, np.ndarray]:
                    for i in range(nspecies):
                        reac_temps[i] = temperature[i]
                    ierr = cea_mixture_calc_property_tp_multitemp(self.ptr, prop_type, nspecies, reac_weights, nspecies, reac_temps, pressure, &value)
                else:
                    raise ValueError("Mixiture.calc_property: temperature must be a float, list, or np.ndarray")
            else:
                if isinstance(temperature, float):
                    ierr = cea_mixture_calc_property(self.ptr, prop_type, nspecies, reac_weights, temperature, &value)
                elif type(temperature) in [list, np.ndarray]:
                    for i in range(nspecies):
                        reac_temps[i] = temperature[i]
                    ierr = cea_mixture_calc_property_multitemp(self.ptr, prop_type, nspecies, reac_weights, nspecies, reac_temps, &value)
                else:
                    raise ValueError("Mixiture.calc_property: temperature must be a float, list, or np.ndarray")
            _check_ierr(ierr, "Mixture.calc_property")
        finally:
            free(reac_weights)
            free(reac_temps)

        return value


cdef class EqSolver:
    """
    Equilibrium solver for chemical equilibrium calculations.

    Solves for chemical equilibrium compositions and properties given
    thermodynamic constraints and reactant compositions.

    Parameters
    ----------
    products : Mixture
        Mixture object containing product species
    **kwargs
        reactants : Mixture, optional
            Mixture object containing reactant species
        transport : bool, default False
            Enable transport property calculations
        ions : bool, default False
            Include ionized species in calculations
        trace : float, default -1.0
            Trace species threshold value; values < 0.0 uses default value
        insert : list of str, default []
            Additional species to insert into calculation; used to start initial guess with condensed species
    """
    cdef cea_eqsolver ptr
    cdef Mixture products
    cdef object _keepalive_insert

    def __init__(self, Mixture products, **kwargs):

        # Extract keyword arguments with defaults
        cdef bint transport = kwargs.get('transport', False)
        cdef bint ions = kwargs.get('ions', False)
        cdef double trace_val = kwargs.get('trace', -1.0)
        cdef cea_string* cea_insert = NULL
        insert = kwargs.get('insert', [])
        cdef list _insert_keepalive = []
        cdef _CString _encoded

        cdef cea_err ierr
        cdef cea_mixture prods = products.ptr
        cdef cea_solver_opts opts
        cdef cea_mixture reacs = NULL
        cdef Mixture py_reactants
        self.products = products

        cea_solver_opts_init(&opts)

        if len(insert) > 0:
            cea_insert = <cea_string*>malloc(sizeof(cea_string)*len(insert))
            if cea_insert == NULL:
                raise MemoryError("Failed to allocate insert pointer buffer")
            for i, val in enumerate(insert):
                _encoded = _CString(val, "EqSolver insert species")
                _insert_keepalive.append(_encoded)
                cea_insert[i] = _encoded.ptr
        self._keepalive_insert = _insert_keepalive

        # Set options
        opts.transport = transport
        opts.ions = ions
        opts.trace = trace_val
        opts.ninsert = <int>len(insert)
        opts.insert = cea_insert

        if "reactants" in kwargs:
            reactants = kwargs["reactants"]
            if not isinstance(reactants, Mixture):
                raise TypeError("reactants must be a Mixture object")
            py_reactants = <Mixture>reactants
            opts.reactants = py_reactants.ptr

        try:
            ierr = cea_eqsolver_create_with_options(&self.ptr, prods, opts)
        finally:
            if cea_insert != NULL:
                free(cea_insert)

        return

    def __dealloc__(self):
        if self.ptr:
            cea_eqsolver_destroy(&self.ptr)
        return

    property num_reactants:
        """
        Number of reactant species.

        Returns
        -------
        int
            Number of reactants
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_eqsolver_get_size(self.ptr, NUM_REACTANTS, &size_val)
            return size_val

    property num_products:
        """
        Number of product species.

        Returns
        -------
        int
            Number of products
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_eqsolver_get_size(self.ptr, NUM_PRODUCTS, &size_val)
            return size_val

    property num_gas:
        """
        Number of gas phase species.

        Returns
        -------
        int
            Number of gas species
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_eqsolver_get_size(self.ptr, NUM_GAS, &size_val)
            return size_val

    property num_condensed:
        """
        Number of condensed phase species.

        Returns
        -------
        int
            Number of condensed species
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_eqsolver_get_size(self.ptr, NUM_CONDENSED, &size_val)
            return size_val

    property num_elements:
        """
        Number of chemical elements.

        Returns
        -------
        int
            Number of elements
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_eqsolver_get_size(self.ptr, NUM_ELEMENTS, &size_val)
            return size_val

    property max_equations:
        """
        Maximum number of equations in system; num_elements + num_condensed + (1 for TP/TV problems, 2 otherwise).

        Returns
        -------
        int
            Maximum equations
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_eqsolver_get_size(self.ptr, MAX_EQUATIONS, &size_val)
            return size_val

    def solve(self, EqSolution soln, cea_equilibrium_type eq_type, float state1, float state2, np.ndarray amounts):
        """
        Solve equilibrium problem with specified constraints.
        NOTE: only mass fractions are accepted; other values must be converted to mass fractions first.

        Parameters
        ----------
        soln : EqSolution
            Solution object to store results
        eq_type : int
            Equilibrium type constant (TP, HP, SP, TV, UV, SV)
        state1 : float
            First thermodynamic state variable (e.g., T, H, S, or U)
        state2 : float
            Second thermodynamic state variable (e.g., P or V)
        amounts : np.ndarray
            Initial reactant mass fractions
        """
        cdef cea_err ierr
        cdef cea_eqpartials partials
        cdef int namounts = <int>len(amounts)
        cdef cea_array amts = <cea_array>malloc(namounts * sizeof(double))
        if amts == NULL:
            raise MemoryError("Failed to allocate reactant amounts buffer")

        try:
            for i in range(namounts):
                amts[i] = amounts[i]

            ierr = cea_eqpartials_create(&partials, self.ptr)
            _check_ierr(ierr, "EqSolver.solve: create partials")
            try:
                ierr = cea_eqsolver_solve_with_partials(self.ptr, eq_type, <cea_real>state1, <cea_real>state2, amts, soln.ptr, partials)
                soln.last_error = <int>ierr
                _check_ierr(ierr, "EqSolver.solve")
            finally:
                ierr = cea_eqpartials_destroy(&partials)
                _check_ierr(ierr, "EqSolver.solve: destroy partials")
        finally:
            free(amts)

        return

cdef class EqSolution:
    """
    Solution object containing equilibrium calculation results.

    Stores thermodynamic properties, species compositions, and other
    results from equilibrium calculations.

    Parameters
    ----------
    solver : EqSolver
        Equilibrium solver instance
    T_init : float, optional
        Initial temperature guess in K
    nj_init : np.ndarray, optional
        Initial species mole fractions
    """
    cdef cea_eqsolution ptr
    cdef EqSolver solver
    cdef public int last_error
    def __cinit__(self, EqSolver solver, T_init: Optional[double] = None, nj_init: Optional[np.ndarray] = None):

        cdef cea_err ierr
        cdef int nj_len
        cdef np.ndarray[np.float64_t, ndim=1, mode="c"] nj_arr
        ierr = cea_eqsolution_create(&self.ptr, solver.ptr)
        _check_ierr(ierr, "EqSolution.__cinit__")
        self.solver = solver
        self.last_error = <int>SUCCESS

        if T_init is not None:
            ierr = cea_eqsolution_set_T(self.ptr, <cea_real>T_init)
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "EqSolution.__cinit__: set_T")
        if nj_init is not None:
            nj_arr = np.ascontiguousarray(nj_init, dtype=np.float64)
            if nj_arr.ndim != 1:
                raise ValueError("EqSolution.__cinit__: nj_init must be a 1D array")
            nj_len = <int>nj_arr.shape[0]
            if nj_len != self.solver.num_products:
                raise ValueError("EqSolution.__cinit__: nj_init length must match solver.num_products")
            ierr = cea_eqsolution_set_nj(self.ptr, self.solver.ptr, nj_len, <cea_array>nj_arr.data)
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "EqSolution.__cinit__: set_nj")

        return

    def __dealloc__(self):
        if self.ptr:
            cea_eqsolution_destroy(&self.ptr)
        return

    property T:
        """
        Temperature in K.

        Returns
        -------
        float
            Temperature
        """
        def __get__(self):
            return self._get_property(TEMPERATURE)
        # def __set__(self, double T_init):
        #     cdef cea_err ierr
        #     ierr = cea_eqsolution_set_T(self.ptr, <cea_real>T_init)
        #     return

    property P:
        """
        Pressure in bar.

        Returns
        -------
        float
            Pressure
        """
        def __get__(self):
            return self._get_property(PRESSURE)

    property volume:
        """
        Specific volume in m³/kg.

        Returns
        -------
        float
            Specific volume
        """
        def __get__(self):
            return self._get_property(VOLUME)

    property density:
        """
        Density in kg/m³.

        Returns
        -------
        float
            Density
        """
        def __get__(self):
            return self._get_property(DENSITY)

    property M:
        """
        Molecular weight, gasses only (1/n).

        Returns
        -------
        float
            Molecular weight (1/n)
        """
        def __get__(self):
            return self._get_property(CEA_M)

    property MW:
        """
        Molecular weight, including condensed species.

        Returns
        -------
        float
            Molecular weight
        """
        def __get__(self):
            return self._get_property(CEA_MW)

    property enthalpy:
        """
        Specific enthalpy in kJ/kg.

        Returns
        -------
        float
            Specific enthalpy
        """
        def __get__(self):
            return self._get_property(ENTHALPY)

    property energy:
        """
        Specific internal energy in kJ/kg.

        Returns
        -------
        float
            Specific internal energy
        """
        def __get__(self):
            return self._get_property(ENERGY)

    property entropy:
        """
        Specific entropy in kJ/(kg·K).

        Returns
        -------
        float
            Specific entropy
        """
        def __get__(self):
            return self._get_property(ENTROPY)

    property gibbs_energy:
        """
        Specific Gibbs energy in kJ/kg.

        Returns
        -------
        float
            Specific Gibbs energy
        """
        def __get__(self):
            return self._get_property(GIBBS_ENERGY)

    property gamma_s:
        """
        Isentropic specific heat ratio.

        Returns
        -------
        float
            Specific heat ratio
        """
        def __get__(self):
            return self._get_property(GAMMA_S)

    property cp_fr:
        """
        Frozen specific heat at constant pressure in kJ/(kg·K).

        Returns
        -------
        float
            Frozen cp
        """
        def __get__(self):
            return self._get_property(FROZEN_CP)

    property cp_eq:
        """
        Equilibrium specific heat at constant pressure in kJ/(kg·K).

        Returns
        -------
        float
            Equilibrium cp
        """
        def __get__(self):
            return self._get_property(EQUILIBRIUM_CP)

    property cp:
        """
        Equilibrium specific heat at constant pressure in kJ/(kg·K).

        Returns
        -------
        float
            Equilibrium cp (alias for cp_eq)
        """
        def __get__(self):
            return self._get_property(EQUILIBRIUM_CP)

    property cv_fr:
        """
        Frozen specific heat at constant volume in kJ/(kg·K).

        Returns
        -------
        float
            Frozen cv
        """
        def __get__(self):
            return self._get_property(FROZEN_CV)

    property cv_eq:
        """
        Equilibrium specific heat at constant volume in kJ/(kg·K).

        Returns
        -------
        float
            Equilibrium cv
        """
        def __get__(self):
            return self._get_property(EQUILIBRIUM_CV)

    property cv:
        """
        Equilibrium specific heat at constant volume in kJ/(kg·K).

        Returns
        -------
        float
            Equilibrium Cv (alias for cv_eq)
        """
        def __get__(self):
            return self._get_property(EQUILIBRIUM_CV)

    property viscosity:
        """
        Viscosity in micropoise (\mu P).

        Returns
        -------
        float
            Viscosity
        """
        def __get__(self):
            return self._get_property(VISCOSITY)

    property conductivity_fr:
        """
        Frozen thermal conductivity in \mu W/(cm·K).

        Returns
        -------
        float
            Frozen thermal conductivity
        """
        def __get__(self):
            return self._get_property(FROZEN_CONDUCTIVITY)

    property conductivity_eq:
        """
        Equilibrium thermal conductivity in \mu W/(cm·K).

        Returns
        -------
        float
            Equilibrium thermal conductivity
        """
        def __get__(self):
            return self._get_property(EQUILIBRIUM_CONDUCTIVITY)

    property Pr_fr:
        """
        Frozen Prandtl number.

        Returns
        -------
        float
            Frozen Prandtl number
        """
        def __get__(self):
            return self._get_property(FROZEN_PRANDTL)

    property Pr_eq:
        """
        Equilibrium Prandtl number.

        Returns
        -------
        float
            Equilibrium Prandtl number
        """
        def __get__(self):
            return self._get_property(EQUILIBRIUM_PRANDTL)

    property nj:
        """
        Species concentrations (kg-mol of species j per kg of mixture).

        Returns
        -------
        np.ndarray
            Concentrations for each species
        """
        def __get__(self):
            return self._get_weights()
        # def __set__(self, np.ndarray nj_init):
        #     cdef cea_err ierr
        #     cdef int nj_len = <int>len(nj_init)
        #     cdef cea_array c_nj_init = <cea_array>malloc(nj_len * sizeof(double))

        #     for i in range(nj_len):
        #         c_nj_init[i] = nj_init[i]
        #     ierr = cea_eqsolution_set_nj(self.ptr, self.solver.ptr, nj_len, c_nj_init)
        #     return

    property ln_nj:
        """
        Natural logarithm of gas phase species concentrations (gasses only).
        NOTE: ln_nj /= log(nj) because nj is truncated for trace species, while ln_nj preserves small amounts.

        Returns
        -------
        np.ndarray
            Log concentrations for gas species
        """
        def __get__(self):
            return self._get_weights(log=True)

    property n:
        """
        Total moles.

        Returns
        -------
        float
            Total moles
        """
        def __get__(self):
            cdef cea_err ierr
            cdef cea_real value
            ierr = cea_eqsolution_get_moles(self.ptr, &value)
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "EqSolution.n")
            return value

    property mass_fractions:
        """
        Species mass fractions.

        Returns
        -------
        dict
            Dictionary mapping species names to mass fractions
        """
        def __get__(self):
            species_names = self.solver.products.species_names
            vals = self._get_species_amounts(mass_fraction=True)
            return dict(zip(species_names, vals))

    property mole_fractions:
        """
        Species mole fractions.

        Returns
        -------
        dict
            Dictionary mapping species names to mole fractions
        """
        def __get__(self):
            species_names = self.solver.products.species_names
            vals = self._get_species_amounts(mass_fraction=False)
            return dict(zip(species_names, vals))

    property converged:
        """
        Convergence status of the solution.

        Returns
        -------
        bool
            True if solution converged, False otherwise
        """
        def __get__(self):
            cdef cea_err ierr
            cdef bint value  # Use bint for boolean values
            ierr = cea_eqsolution_get_converged(self.ptr, <bint *>(&value))  # Explicit cast to bint *
            return bool(value)

    def _get_property(self, cea_property_type prop_type):
        cdef cea_err ierr
        cdef cea_real value

        ierr = cea_eqsolution_get_property(self.ptr, prop_type, &value)
        if ierr != SUCCESS:
            self.last_error = <int>ierr
        _check_ierr(ierr, "EqSolution._get_property")

        return value

    def _get_weights(self, bint log=False):
        cdef cea_err ierr
        cdef int nspecies

        if log:
            nspecies = self.solver.num_gas
        else:
            nspecies = self.solver.num_products
        cdef cea_real *reac_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if reac_weights == NULL:
            raise MemoryError("Failed to allocate weights buffer")

        try:
            if log:
                weights = np.zeros(nspecies, dtype=np.double)
                ierr = cea_eqsolution_get_weights(self.ptr, nspecies, reac_weights, log)
            else:
                weights = np.zeros(nspecies, dtype=np.double)
                ierr = cea_eqsolution_get_weights(self.ptr, nspecies, reac_weights, log)
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "EqSolution._get_weights")

            # Assign the weights to the numpy array
            for i in range(nspecies):
                weights[i] = reac_weights[i]
        finally:
            free(reac_weights)

        return weights

    def _get_species_amounts(self, bint mass_fraction=False):

        cdef cea_err ierr
        cdef int nspecies

        nprod = self.solver.num_products
        cdef cea_real *cea_amounts = <cea_real *>malloc(nprod * sizeof(cea_real))
        if cea_amounts == NULL:
            raise MemoryError("Failed to allocate species amounts buffer")

        try:
            amounts = np.zeros(nprod, dtype=np.double)
            ierr = cea_eqsolution_get_species_amounts(self.ptr, nprod, cea_amounts, mass_fraction)
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "EqSolution._get_species_amounts")

            # Assign the amounts to the numpy array
            for i in range(nprod):
                amounts[i] = cea_amounts[i]
        finally:
            free(cea_amounts)

        return amounts


cdef class RocketSolver:
    """
    Rocket performance solver.

    Solves for rocket nozzle flow properties including thrust coefficient, specific impulse,
    and exit conditions for both infinite area combustor (IAC) and finite area
    combustor (FAC) assumptions.

    Parameters
    ----------
    products : Mixture
        Mixture object containing product species
    **kwargs
        reactants : Mixture, optional
            Mixture object containing reactant species
        transport : bool, default False
            Enable transport property calculations
        ions : bool, default False
            Include ionized species in calculations
        trace : float, default -1.0
            Trace species threshold value; values < 0.0 uses default value
        insert : list of str, default []
            Additional species to insert into calculation; used to start initial guess with condensed species
    """
    cdef cea_rocket_solver ptr
    cdef Mixture products
    cdef object _keepalive_insert

    def __init__(self, Mixture products, **kwargs):

        # Extract keyword arguments with defaults
        cdef bint transport = kwargs.get('transport', False)
        cdef bint ions = kwargs.get('ions', False)
        cdef double trace_val = kwargs.get('trace', -1.0)
        cdef cea_string* cea_insert = NULL
        insert = kwargs.get('insert', [])
        cdef list _insert_keepalive = []
        cdef _CString _encoded

        cdef cea_err ierr
        cdef cea_mixture prods = products.ptr
        cdef cea_solver_opts opts
        cdef cea_mixture reacs = NULL
        cdef Mixture py_reactants
        self.products = products

        cea_solver_opts_init(&opts)

        if len(insert) > 0:
            cea_insert = <cea_string*>malloc(sizeof(cea_string)*len(insert))
            if cea_insert == NULL:
                raise MemoryError("Failed to allocate insert pointer buffer")
            for i, val in enumerate(insert):
                _encoded = _CString(val, "RocketSolver insert species")
                _insert_keepalive.append(_encoded)
                cea_insert[i] = _encoded.ptr
        self._keepalive_insert = _insert_keepalive

        # Set options
        opts.transport = transport
        opts.ions = ions
        opts.trace = trace_val
        opts.ninsert = <int>len(insert)
        opts.insert = cea_insert

        if "reactants" in kwargs:
            reactants = kwargs["reactants"]
            if not isinstance(reactants, Mixture):
                raise TypeError("reactants must be a Mixture object")
            py_reactants = <Mixture>reactants
            opts.reactants = py_reactants.ptr

        try:
            ierr = cea_rocket_solver_create_with_options(&self.ptr, prods, opts)
        finally:
            if cea_insert != NULL:
                free(cea_insert)

        return

    def __dealloc__(self):
        if self.ptr:
            cea_rocket_solver_destroy(&self.ptr)
        return

    property num_reactants:
        """
        Number of reactant species.

        Returns
        -------
        int
            Number of reactants
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_rocket_solver_get_size(self.ptr, NUM_REACTANTS, &size_val)
            return size_val

    property num_products:
        """
        Number of product species.

        Returns
        -------
        int
            Number of products
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_rocket_solver_get_size(self.ptr, NUM_PRODUCTS, &size_val)
            return size_val

    property num_gas:
        """
        Number of gas phase species.

        Returns
        -------
        int
            Number of gas species
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_rocket_solver_get_size(self.ptr, NUM_GAS, &size_val)
            return size_val

    property num_condensed:
        """
        Number of condensed phase species.

        Returns
        -------
        int
            Number of condensed species
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_rocket_solver_get_size(self.ptr, NUM_CONDENSED, &size_val)
            return size_val

    property num_elements:
        """
        Number of chemical elements.

        Returns
        -------
        int
            Number of elements
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_rocket_solver_get_size(self.ptr, NUM_ELEMENTS, &size_val)
            return size_val

    def solve(self, RocketSolution soln, np.ndarray weights, float pc, pi_p=None, subar=None, supar=None,
              iac=True, n_frz=None, hc=None, tc=None, mdot=None, ac_at=None, tc_est=None):
        """
        Solve rocket nozzle flow problem.

        Parameters
        ----------
        soln : RocketSolution
            Solution object to store results
        weights : np.ndarray
            Reactant mass fractions
        pc : float
            Chamber pressure in bar
        pi_p : float, list, or np.ndarray, optional
            Pressure ratios (pc/pe) for exit conditions
        subar : float, list, or np.ndarray, optional
            Subsonic area ratios (Ae/At)
        supar : float, list, or np.ndarray, optional
            Supersonic area ratios (Ae/At)
        iac : bool, default True
            Use infinite area combustor (IAC) model. If False, uses finite area combustor (FAC)
        n_frz : int, optional
            Point at which flow becomes frozen (for FAC calculations); turns on the frozen calculation option
        hc : float, optional
            Chamber enthalpy in (kJ/kg)/R (mutually exclusive with tc)
        tc : float, optional
            Chamber temperature in K (mutually exclusive with hc)
        mdot : float, optional
            Ratio of mass flow to chamber area in (kg/s/m^2) (for FAC, mutually exclusive with ac_at)
        ac_at : float, optional
            Chamber to throat area ratio (for FAC, mutually exclusive with mdot)
        tc_est : float, optional
            Initial chamber temperature estimate in K. Not relevant when tc is provided.

        Raises
        ------
        ValueError
            If exit conditions not specified, or conflicting parameters provided
        """
        cdef cea_err ierr
        cdef int nweights = <int>len(weights)

        cdef cea_array wts = <cea_array>malloc(nweights * sizeof(double))
        if wts == NULL:
            raise MemoryError("Failed to allocate weights buffer")

        for i in range(nweights):
            wts[i] = weights[i]

        # Convert pi_p, subar, and supar to arrays
        cdef int npi_p = 0
        cdef int nsubar = 0
        cdef int nsupar = 0

        if pi_p is not None:
            if type(pi_p) in [list, np.ndarray]:
                npi_p = len(pi_p)
                if npi_p == 0:
                    raise ValueError("RocketSolver.solve: pi_p cannot be empty")
            elif type(pi_p) is float:
                pi_p = [pi_p]
                npi_p = 1
            else:
                raise ValueError("RocketSolver.solve: pi_p must be a list, np.ndarray, or float")

        if subar is not None:
            if type(subar) in [list, np.ndarray]:
                nsubar = len(subar)
                if nsubar == 0:
                    raise ValueError("RocketSolver.solve: subar cannot be empty")
            elif type(subar) is float:
                subar = [subar]
                nsubar = 1
            else:
                raise ValueError("RocketSolver.solve: subar must be a list, np.ndarray, or float")

        if supar is not None:
            if type(supar) in [list, np.ndarray]:
                nsupar = len(supar)
                if nsupar == 0:
                    raise ValueError("RocketSolver.solve: supar cannot be empty")
            elif type(supar) is float:
                supar = [supar]
                nsupar = 1
            else:
                raise ValueError("RocketSolver.solve: supar must be a list, np.ndarray, or float")

        if ((npi_p + nsubar + nsupar) == 0):
            raise ValueError("Exit condition not specified.")

        cdef cea_array pi_p_c = <cea_array>malloc(npi_p * sizeof(double))
        cdef cea_array subar_c = <cea_array>malloc(nsubar * sizeof(double))
        cdef cea_array supar_c = <cea_array>malloc(nsupar * sizeof(double))
        if (npi_p > 0 and pi_p_c == NULL) or (nsubar > 0 and subar_c == NULL) or (nsupar > 0 and supar_c == NULL):
            free(wts)
            if pi_p_c != NULL:
                free(pi_p_c)
            if subar_c != NULL:
                free(subar_c)
            if supar_c != NULL:
                free(supar_c)
            raise MemoryError("Failed to allocate exit condition buffers")

        for i in range(npi_p):
            pi_p_c[i] = pi_p[i]

        for i in range(nsubar):
            subar_c[i] = subar[i]

        for i in range(nsupar):
            supar_c[i] = supar[i]

        # Process the other optional inputs
        if (hc is not None) & (tc is not None):
            raise ValueError("RocketSolver.solve: hc and tc cannot both be provided")
        if (hc is None) & (tc is None):
            raise ValueError("RocketSolver.solve: either hc or tc must be provided")

        cdef cea_int n_frz_c = 0
        if n_frz is not None:
            n_frz_c = <cea_int>n_frz

        cdef bint use_hc
        cdef cea_real hc_or_tc
        if hc is not None:
            use_hc = True
            hc_or_tc = <cea_real>hc
        else:
            use_hc = False
            hc_or_tc = <cea_real>tc

        cdef bint use_tc_est = False
        cdef cea_real tc_est_c
        if tc_est is not None:
            use_tc_est = True
            tc_est_c = <cea_real>tc_est

        cdef bint use_mdot
        cdef cea_real mdot_or_acat
        if not iac:
            if (mdot is not None) & (ac_at is not None):
                raise ValueError("RocketSolver.solve: mdot and ac_at cannot both be provided")
            if (mdot is None) & (ac_at is None):
                raise ValueError("RocketSolver.solve: one of mdot or ac_at must be provided for FAC problems")

            if mdot is not None:
                use_mdot = True
                mdot_or_acat = <cea_real>mdot
            else:
                use_mdot = False
                mdot_or_acat = <cea_real>ac_at

        try:
            # Call the solver
            if iac:
                ierr = cea_rocket_solver_solve_iac(self.ptr, soln.ptr, wts, <cea_real>pc, pi_p_c, npi_p, subar_c, nsubar, supar_c, nsupar, n_frz_c, hc_or_tc, use_hc, tc_est_c, use_tc_est)
            else:
                ierr = cea_rocket_solver_solve_fac(self.ptr, soln.ptr, wts, <cea_real>pc, pi_p_c, npi_p, subar_c, nsubar, supar_c, nsupar, n_frz_c, hc_or_tc, use_hc, mdot_or_acat, use_mdot, tc_est_c, use_tc_est)
            soln.last_error = <int>ierr
            _check_ierr(ierr, "RocketSolver.solve")
        finally:
            free(wts)
            if pi_p_c != NULL:
                free(pi_p_c)
            if subar_c != NULL:
                free(subar_c)
            if supar_c != NULL:
                free(supar_c)

        return


cdef class RocketSolution:
    """
    Solution object containing rocket nozzle flow calculation results.

    Stores flow properties, performance parameters, and species compositions
    at multiple stations along the nozzle.

    Parameters
    ----------
    solver : RocketSolver
        Rocket solver instance
    """
    cdef cea_rocket_solution ptr
    cdef RocketSolver solver
    cdef public int last_error
    def __cinit__(self, RocketSolver solver):

        cdef cea_err ierr
        ierr = cea_rocket_solution_create(&self.ptr, solver.ptr)
        _check_ierr(ierr, "RocketSolution.__cinit__")
        self.solver = solver
        self.last_error = <int>SUCCESS

        return

    def __dealloc__(self):

        # Free the allocated pointers array
        if self.ptr:
            cea_rocket_solution_destroy(&self.ptr)

        return

    def _get_size(self):
        cdef cea_err ierr
        cdef cea_int num_pts
        ierr = cea_rocket_solution_get_size(self.ptr, &num_pts)
        if ierr != SUCCESS:
            self.last_error = <int>ierr
        _check_ierr(ierr, "RocketSolution._get_size")
        return num_pts

    property num_pts:
        """
        Number of calculation points (stations).

        Returns
        -------
        int
            Number of stations
        """
        def __get__(self):
            return self._get_size()

    property T:
        """
        Temperature at each station in K.

        Returns
        -------
        np.ndarray
            Temperature array
        """
        def __get__(self):
            return self._get_property(ROCKET_TEMPERATURE)

    property P:
        """
        Pressure at each station in bar.

        Returns
        -------
        np.ndarray
            Pressure array
        """
        def __get__(self):
            return self._get_property(ROCKET_PRESSURE)

    property volume:
        """
        Specific volume at each station in m³/kg.

        Returns
        -------
        np.ndarray
            Specific volume array
        """
        def __get__(self):
            return self._get_property(ROCKET_VOLUME)

    property density:
        """
        Density at each station in kg/m³.

        Returns
        -------
        np.ndarray
            Density array
        """
        def __get__(self):
            return self._get_property(ROCKET_DENSITY)

    property M:
        """
        Molecular weight (1/n) (gasses only) at each station.

        Returns
        -------
        np.ndarray
            Molecular weight array
        """
        def __get__(self):
            return self._get_property(ROCKET_M)

    property MW:
        """
        Molecular weight (including condensed species) at each station.

        Returns
        -------
        np.ndarray
            Molecular weight array
        """
        def __get__(self):
            return self._get_property(ROCKET_MW)

    property enthalpy:
        """
        Specific enthalpy at each station in kJ/kg.

        Returns
        -------
        np.ndarray
            Specific enthalpy array
        """
        def __get__(self):
            return self._get_property(ROCKET_ENTHALPY)

    property energy:
        """
        Specific internal energy at each station in kJ/kg.

        Returns
        -------
        np.ndarray
            Specific internal energy array
        """
        def __get__(self):
            return self._get_property(ROCKET_ENERGY)

    property entropy:
        """
        Specific entropy at each station in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Specific entropy array
        """
        def __get__(self):
            return self._get_property(ROCKET_ENTROPY)

    property gibbs_energy:
        """
        Specific Gibbs energy at each station in kJ/kg.

        Returns
        -------
        np.ndarray
            Specific Gibbs energy array
        """
        def __get__(self):
            return self._get_property(ROCKET_GIBBS_ENERGY)

    property gamma_s:
        """
        Isentropic exponent \gamma_s at each station.

        Returns
        -------
        np.ndarray
            Isentropic exponent \gamma_s array
        """
        def __get__(self):
            return self._get_property(ROCKET_GAMMA_S)

    property cp_fr:
        """
        Frozen specific heat at constant pressure in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Frozen cp array
        """
        def __get__(self):
            return self._get_property(ROCKET_FROZEN_CP)

    property cp_eq:
        """
        Equilibrium specific heat at constant pressure in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Equilibrium cp array
        """
        def __get__(self):
            return self._get_property(ROCKET_EQUILIBRIUM_CP)

    property cp:
        """
        Equilibrium specific heat at constant pressure in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Equilibrium cp array (alias for cp_eq)
        """
        def __get__(self):
            return self._get_property(ROCKET_EQUILIBRIUM_CP)

    property cv_fr:
        """
        Frozen specific heat at constant volume in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Frozen cv array
        """
        def __get__(self):
            return self._get_property(ROCKET_FROZEN_CV)

    property cv_eq:
        """
        Equilibrium specific heat at constant volume in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Equilibrium cv array
        """
        def __get__(self):
            return self._get_property(ROCKET_EQUILIBRIUM_CV)

    property cv:
        """
        Equilibrium specific heat at constant volume in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Equilibrium cv array (alias for cv_eq)
        """
        def __get__(self):
            return self._get_property(ROCKET_EQUILIBRIUM_CV)

    property Mach:
        """
        Mach number at each station.

        Returns
        -------
        np.ndarray
            Mach number array
        """
        def __get__(self):
            return self._get_property(MACH)

    property sonic_velocity:
        """
        Sonic velocity at each station in m/s.

        Returns
        -------
        np.ndarray
            Sonic velocity array
        """
        def __get__(self):
            return self._get_property(SONIC_VELOCITY)

    property ae_at:
        """
        Area ratio (Ae/At) at each station.

        Returns
        -------
        np.ndarray
            Area ratio array
        """
        def __get__(self):
            return self._get_property(AE_AT)

    property c_star:
        """
        Characteristic velocity in m/s.

        Returns
        -------
        np.ndarray
            Characteristic velocity array
        """
        def __get__(self):
            return self._get_property(C_STAR)

    property coefficient_of_thrust:
        """
        Thrust coefficient at each station.

        Returns
        -------
        np.ndarray
            Thrust coefficient array
        """
        def __get__(self):
            return self._get_property(COEFFICIENT_OF_THRUST)

    property Isp:
        """
        Specific impulse at each station in m/s.

        Returns
        -------
        np.ndarray
            Specific impulse array
        """
        def __get__(self):
            return self._get_property(ISP)

    property Isp_vacuum:
        """
        Vacuum specific impulse at each station in m/s.

        Returns
        -------
        np.ndarray
            Vacuum specific impulse array
        """
        def __get__(self):
            return self._get_property(ISP_VACUUM)

    property viscosity:
        """
        Viscosity in millipoise (\mu P).

        Returns
        -------
        float
            Viscosity
        """
        def __get__(self):
            return self._get_property(ROCKET_VISCOSITY)

    property conductivity_fr:
        """
        Frozen thermal conductivity in \mu W/(cm·K).

        Returns
        -------
        float
            Frozen thermal conductivity
        """
        def __get__(self):
            return self._get_property(ROCKET_FROZEN_CONDUCTIVITY)

    property conductivity_eq:
        """
        Equilibrium thermal conductivity in \mu W/(cm·K).

        Returns
        -------
        float
            Equilibrium thermal conductivity
        """
        def __get__(self):
            return self._get_property(ROCKET_EQUILIBRIUM_CONDUCTIVITY)

    property Pr_fr:
        """
        Frozen Prandtl number.

        Returns
        -------
        float
            Frozen Prandtl number
        """
        def __get__(self):
            return self._get_property(ROCKET_FROZEN_PRANDTL)

    property Pr_eq:
        """
        Equilibrium Prandtl number.

        Returns
        -------
        float
            Equilibrium Prandtl number
        """
        def __get__(self):
            return self._get_property(ROCKET_EQUILIBRIUM_PRANDTL)

    property nj:
        """
        Species concentrations (kg-mol of species j per kg of mixture) at each station.

        Returns
        -------
        np.ndarray
            2D array of concentrations (stations × species)
        """
        def __get__(self):
            return self._get_weights()

    property ln_nj:
        """
        Natural logarithm of gas phase species concentrations at each station (gasses only).
        NOTE: ln_nj /= log(nj) because nj is truncated for trace species, while ln_nj preserves small amounts.

        Returns
        -------
        np.ndarray
            2D array of log concentrations (stations × gas species)
        """
        def __get__(self):
            return self._get_weights(log=True)

    property n:
        """
        Total moles at each station.

        Returns
        -------
        np.ndarray
            Total moles array
        """
        def __get__(self):
            cdef cea_err ierr
            cdef int num_pts = self._get_size()
            cdef cea_real *value = <cea_real *>malloc(num_pts * sizeof(cea_real))
            if value == NULL:
                raise MemoryError("Failed to allocate moles buffer")

            try:
                ierr = cea_rocket_solution_get_moles(self.ptr, value)
                if ierr != SUCCESS:
                    self.last_error = <int>ierr
                _check_ierr(ierr, "RocketSolution.n")

                val = np.zeros(num_pts, dtype=np.double)
                for i in range(num_pts):
                    val[i] = value[i]
            finally:
                free(value)

            return val

    property mass_fractions:
        """
        Species mass fractions at each station.

        Returns
        -------
        dict
            Dictionary mapping species names to 2D arrays (stations × species)
        """
        def __get__(self):
            species_names = self.solver.products.species_names
            vals = np.transpose(self._get_species_amounts(mass_fraction=True))
            return dict(zip(species_names, vals))

    property mole_fractions:
        """
        Species mole fractions at each station.

        Returns
        -------
        dict
            Dictionary mapping species names to 2D arrays (stations × species)
        """
        def __get__(self):
            species_names = self.solver.products.species_names
            vals = np.transpose(self._get_species_amounts(mass_fraction=False))
            return dict(zip(species_names, vals))

    property converged:
        """
        Convergence status of the solution.

        Returns
        -------
        bool
            True if solution converged, False otherwise
        """
        def __get__(self):
            cdef cea_err ierr
            cdef bint value  # Use bint for boolean values
            ierr = cea_rocket_solution_get_converged(self.ptr, <bint *>(&value))  # Explicit cast to bint *
            return bool(value)

    def _get_property(self, cea_rocket_property_type prop_type):
        cdef cea_err ierr
        cdef int num_pts = self._get_size()
        cdef cea_real *value = <cea_real *>malloc(num_pts * sizeof(cea_real))
        if value == NULL:
            raise MemoryError("Failed to allocate property buffer")

        try:
            ierr = cea_rocket_solution_get_property(self.ptr, prop_type, num_pts, value)
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "RocketSolution._get_property")

            prop_val = np.zeros(num_pts, dtype=np.double)
            for i in range(num_pts):
                prop_val[i] = value[i]
        finally:
            free(value)

        return prop_val

    def _get_weights(self, bint log=False):
        cdef cea_err ierr
        cdef int num_pts = self._get_size()
        cdef int nspecies

        if log:
            nspecies = self.solver.num_gas
        else:
            nspecies = self.solver.num_products
        cdef cea_real *reac_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if reac_weights == NULL:
            raise MemoryError("Failed to allocate weights buffer")

        try:
            if log:
                weights = np.zeros((num_pts, nspecies), dtype=np.double)
                for i in range(1,num_pts+1):
                    ierr = cea_rocket_solution_get_weights(self.ptr, nspecies, <cea_int>i, reac_weights, log)
                    for j in range(nspecies):
                        weights[i-1, j] = reac_weights[j]
            else:
                weights = np.zeros((num_pts, nspecies), dtype=np.double)
                for i in range(1,num_pts+1):
                    ierr = cea_rocket_solution_get_weights(self.ptr, nspecies, <cea_int>i, reac_weights, log)
                    for j in range(nspecies):
                        weights[i-1, j] = reac_weights[j]
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "RocketSolution._get_weights")
        finally:
            free(reac_weights)

        return weights

    def _get_species_amounts(self, bint mass_fraction=False):

        cdef cea_err ierr
        cdef int num_pts = self._get_size()
        cdef int nspecies = self.solver.num_products

        cdef cea_real *cea_amounts = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if cea_amounts == NULL:
            raise MemoryError("Failed to allocate species amounts buffer")

        try:
            amounts = np.zeros((num_pts, nspecies), dtype=np.double)
            for i in range(1,num_pts+1):
                ierr = cea_rocket_solution_get_species_amounts(self.ptr, nspecies, <cea_int>i, cea_amounts, mass_fraction)
                for j in range(nspecies):
                    amounts[i-1, j] = cea_amounts[j]
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "RocketSolution._get_species_amounts")
        finally:
            free(cea_amounts)

        return amounts

    # def get_EqSolutions(self):
    #     cdef cea_err ierr
    #     cdef cea_int num_pts = self.num_pts
    #     self.eqslptrs = <cea_eqsolution **>malloc(num_pts * sizeof(cea_eqsolution*))

    #     ierr = cea_rocket_solution_get_eq_solutions(self.ptr, num_pts, self.eqslptrs)

    #     eqsolns = []
    #     for i in range(num_pts):
    #         eqsoln = EqSolution(self.solver.eqsolver)
    #         eqsoln.ptr = <cea_eqsolution_t *>self.eqslptrs[i]
    #         eqsolns.append(eqsoln)

    #     return eqsolns

cdef class ShockSolver:
    """
    Shock wave solver for calculating shock and detonation properties.

    Solves for normal shock wave properties including incident and
    reflected shock conditions.

    Parameters
    ----------
    products : Mixture
        Mixture object containing product species
    **kwargs
        reactants : Mixture, optional
            Mixture object containing reactant species
        transport : bool, default False
            Enable transport property calculations
        ions : bool, default False
            Include ionized species in calculations
        trace : float, default -1.0
            Trace species threshold value; values < 0.0 uses default value
        insert : list of str, default []
            Additional species to insert into calculation; used to start initial guess with condensed species
    """
    cdef cea_shock_solver ptr
    cdef Mixture products
    cdef object _keepalive_insert

    def __init__(self, Mixture products, **kwargs):

        # Extract keyword arguments with defaults
        cdef bint transport = kwargs.get('transport', False)
        cdef bint ions = kwargs.get('ions', False)
        cdef double trace_val = kwargs.get('trace', -1.0)
        cdef cea_string* cea_insert = NULL
        insert = kwargs.get('insert', [])
        cdef list _insert_keepalive = []
        cdef _CString _encoded

        cdef cea_err ierr
        cdef cea_mixture prods = products.ptr
        cdef cea_solver_opts opts
        cdef cea_mixture reacs = NULL
        cdef Mixture py_reactants
        self.products = products

        cea_solver_opts_init(&opts)

        if len(insert) > 0:
            cea_insert = <cea_string*>malloc(sizeof(cea_string)*len(insert))
            if cea_insert == NULL:
                raise MemoryError("Failed to allocate insert pointer buffer")
            for i, val in enumerate(insert):
                _encoded = _CString(val, "ShockSolver insert species")
                _insert_keepalive.append(_encoded)
                cea_insert[i] = _encoded.ptr
        self._keepalive_insert = _insert_keepalive

        # Set options
        opts.transport = transport
        opts.ions = ions
        opts.trace = trace_val
        opts.ninsert = <int>len(insert)
        opts.insert = cea_insert

        if "reactants" in kwargs:
            reactants = kwargs["reactants"]
            if not isinstance(reactants, Mixture):
                raise TypeError("reactants must be a Mixture object")
            py_reactants = <Mixture>reactants
            opts.reactants = py_reactants.ptr

        try:
            ierr = cea_shock_solver_create_with_options(&self.ptr, prods, opts)
        finally:
            if cea_insert != NULL:
                free(cea_insert)

        return

    def __dealloc__(self):
        if self.ptr:
            cea_shock_solver_destroy(&self.ptr)
        return

    property num_reactants:
        """
        Number of reactant species.

        Returns
        -------
        int
            Number of reactants
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_shock_solver_get_size(self.ptr, NUM_REACTANTS, &size_val)
            return size_val

    property num_products:
        """
        Number of product species.

        Returns
        -------
        int
            Number of products
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_shock_solver_get_size(self.ptr, NUM_PRODUCTS, &size_val)
            return size_val

    property num_gas:
        """
        Number of gas phase species.

        Returns
        -------
        int
            Number of gas species
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_shock_solver_get_size(self.ptr, NUM_GAS, &size_val)
            return size_val

    property num_condensed:
        """
        Number of condensed phase species.

        Returns
        -------
        int
            Number of condensed species
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_shock_solver_get_size(self.ptr, NUM_CONDENSED, &size_val)
            return size_val

    property num_elements:
        """
        Number of chemical elements.

        Returns
        -------
        int
            Number of elements
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_shock_solver_get_size(self.ptr, NUM_ELEMENTS, &size_val)
            return size_val

    def solve(self, ShockSolution soln, np.ndarray weights, float T0, float p0,
              u1=None, Mach1=None, bint reflected=True, bint incident_frozen=False, bint reflected_frozen=False):
        """
        Solve shock wave problem.

        Parameters
        ----------
        soln : ShockSolution
            Solution object to store results
        weights : np.ndarray
            Initial composition mass fractions
        T0 : float
            Unshocked temperature in K
        p0 : float
            Unshocked pressure in bar
        u1 : float, optional
            Shock velocity in m/s (mutually exclusive with Mach1)
        Mach1 : float, optional
            Shock Mach number (mutually exclusive with u1)
        reflected : bool, default True
            Calculate reflected shock conditions
        incident_frozen : bool, default False
            Use frozen chemistry for incident shock
        reflected_frozen : bool, default False
            Use frozen chemistry for reflected shock

        Raises
        ------
        ValueError
            If both u1 and Mach1 provided, or neither provided
        """
        cdef cea_err ierr
        cdef int nweights = <int>len(weights)

        cdef cea_array wts = <cea_array>malloc(nweights * sizeof(double))
        if wts == NULL:
            raise MemoryError("Failed to allocate weights buffer")

        for i in range(nweights):
            wts[i] = weights[i]

        # Make sure one of u1 or Mach1 is provided (and only one)
        if (u1 is not None) & (Mach1 is not None):
            raise ValueError("ShockSolver.solve: u1 and Mach1 cannot both be provided")
        if (u1 is None) & (Mach1 is None):
            raise ValueError("ShockSolver.solve: either u1 or Mach1 must be provided")

        # Convert the values for u1 or Mach1
        cdef bint use_mach
        cdef cea_real u1_or_mach1
        if u1 is not None:
            use_mach = False
            u1_or_mach1 = <cea_real>u1
        else:
            use_mach = True
            u1_or_mach1 = <cea_real>Mach1

        try:
            ierr = cea_shock_solver_solve(self.ptr, soln.ptr, wts, <cea_real>T0, <cea_real>p0, u1_or_mach1, use_mach,
                                          reflected, incident_frozen, reflected_frozen)
            soln.last_error = <int>ierr
            _check_ierr(ierr, "ShockSolver.solve")
        finally:
            free(wts)

        return

cdef class ShockSolution:
    """
    Solution object containing shock wave calculation results.

    Stores shock properties at initial, incident, and optionally
    reflected shock states.

    Parameters
    ----------
    reflected : bool, default True
        Include reflected shock calculations (3 points vs 2 points)
    """
    cdef cea_shock_solution ptr
    cdef ShockSolver solver
    cdef cea_int num_pts
    cdef public int last_error
    def __cinit__(self, ShockSolver solver, reflected=True):

        self.solver = solver
        self.num_pts = 2
        if reflected:
            self.num_pts = 3

        cdef cea_err ierr
        ierr = cea_shock_solution_create(&self.ptr, self.num_pts)
        _check_ierr(ierr, "ShockSolution.__cinit__")
        self.last_error = <int>SUCCESS

        return

    def __dealloc__(self):
        if self.ptr:
            cea_shock_solution_destroy(&self.ptr)
        return

    property T:
        """
        Temperature at each shock state in K.

        Returns
        -------
        np.ndarray
            Temperature array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_TEMPERATURE)

    property P:
        """
        Pressure at each shock state in bar.

        Returns
        -------
        np.ndarray
            Pressure array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_PRESSURE)

    property velocity:
        """
        Velocity at each shock state in m/s.

        Returns
        -------
        np.ndarray
            Velocity array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_VELOCITY)

    property Mach:
        """
        Mach number at each shock state.

        Returns
        -------
        np.ndarray
            Mach number array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_MACH)

    property sonic_velocity:
        """
        Sonic velocity at each shock state in m/s.

        Returns
        -------
        np.ndarray
            Sonic velocity array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_SONIC_VELOCITY)

    property rho12:
        """
        Density ratio across the incident shock (ρ2/ρ1).

        Returns
        -------
        float
            Density ratio across the incident shock
        """
        def __get__(self):
            return self._get_scalar_property(SHOCK_RHO12)

    property rho52:
        """
        Density ratio across the reflected shock (ρ5/ρ2).

        Returns
        -------
        float
            Density ratio across the reflected shock
        """
        def __get__(self):
            return self._get_scalar_property(SHOCK_RHO52)

    property P21:
        """
        Pressure ratio across the incident shock (P2/P1).

        Returns
        -------
        float
            Pressure ratio across the incident shock
        """
        def __get__(self):
            return self._get_scalar_property(SHOCK_P21)

    property P52:
        """
        Pressure ratio across the reflected shock (P5/P2).

        Returns
        -------
        float
            Pressure ratio across the reflected shock
        """
        def __get__(self):
            return self._get_scalar_property(SHOCK_P52)

    property T21:
        """
        Temperature ratio across the incident shock (T2/T1).

        Returns
        -------
        float
            Temperature ratio across the incident shock
        """
        def __get__(self):
            return self._get_scalar_property(SHOCK_T21)

    property T52:
        """
        Temperature ratio across the reflected shock (T5/T2).

        Returns
        -------
        float
            Temperature ratio across the reflected shock
        """
        def __get__(self):
            return self._get_scalar_property(SHOCK_T52)

    property M21:
        """
        Mach number ratio across the incident shock (M2/M1).

        Returns
        -------
        float
            Mach number ratio across the incident shock
        """
        def __get__(self):
            return self._get_scalar_property(SHOCK_M21)

    property M52:
        """
        Mach number ratio across the reflected shock (M5/M2).

        Returns
        -------
        float
            Mach number ratio across the reflected shock
        """
        def __get__(self):
            return self._get_scalar_property(SHOCK_M52)

    property v2:
        """
        Velocity at the incident shock state in m/s.

        Returns
        -------
        float
            Velocity at the incident shock state
        """
        def __get__(self):
            return self._get_scalar_property(SHOCK_V2)

    property u5_p_v2:
        """
        Velocity of the reflected shock wave in m/s.

        Returns
        -------
        float
            u5 + v2
        """
        def __get__(self):
            return self._get_scalar_property(SHOCK_U5_P_V2)

    property volume:
        """
        Specific volume at each shock state in m³/kg.

        Returns
        -------
        np.ndarray
            Specific volume array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_VOLUME)

    property density:
        """
        Density at each shock state in kg/m³.

        Returns
        -------
        np.ndarray
            Density array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_DENSITY)

    property M:
        """
        Molecular weight (1/n) at each shock state.

        Returns
        -------
        np.ndarray
            Molecular weight (1/n) array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_M)

    property MW:
        """
        Molecular weight (including condensed species) at each shock state in kg/kmol.

        Returns
        -------
        np.ndarray
            Molecular weight array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_MW)

    property enthalpy:
        """
        Specific enthalpy at each shock state in kJ/kg.

        Returns
        -------
        np.ndarray
            Specific enthalpy array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_ENTHALPY)

    property energy:
        """
        Specific internal energy at each shock state in kJ/kg.

        Returns
        -------
        np.ndarray
            Specific internal energy array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_ENERGY)

    property entropy:
        """
        Specific entropy at each shock state in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Specific entropy array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_ENTROPY)

    property gibbs_energy:
        """
        Gibbs energy at each shock state in kJ/kg.

        Returns
        -------
        np.ndarray
            Gibbs energy array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_GIBBS_ENERGY)

    property gamma_s:
        """
        Isentropic exponent \gamma_s at each shock state.

        Returns
        -------
        np.ndarray
            Isentropic exponent \gamma_s array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_GAMMA_S)

    property cp_fr:
        """
        Frozen specific heat at constant pressure at each shock state in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Frozen C_p array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_FROZEN_CP)

    property cp_eq:
        """
        Equilibrium specific heat at constant pressure at each shock state in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Equilibrium C_p array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_EQUILIBRIUM_CP)

    property cp:
        """
        Equilibrium specific heat at constant pressure at each shock state in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Equilibrium C_p array (alias for cp_eq) [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_EQUILIBRIUM_CP)

    property cv_fr:
        """
        Frozen specific heat at constant volume at each shock state in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Frozen C_v array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_FROZEN_CV)

    property cv_eq:
        """
        Equilibrium specific heat at constant volume at each shock state in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Equilibrium C_v array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_EQUILIBRIUM_CV)

    property cv:
        """
        Equilibrium specific heat at constant volume at each shock state in kJ/(kg·K).

        Returns
        -------
        np.ndarray
            Equilibrium C_v array (alias for cv_eq) [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_EQUILIBRIUM_CV)

    property viscosity:
        """
        Viscosity at each shock state in micropoise (\mu P).

        Returns
        -------
        np.ndarray
            Viscosity array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_VISCOSITY)

    property conductivity_fr:
        """
        Frozen thermal conductivity at each shock state in \mu W/(cm·K).

        Returns
        -------
        np.ndarray
            Frozen thermal conductivity array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_FROZEN_CONDUCTIVITY)

    property conductivity_eq:
        """
        Equilibrium thermal conductivity at each shock state in \mu W/(cm·K).

        Returns
        -------
        np.ndarray
            Equilibrium thermal conductivity array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_EQUILIBRIUM_CONDUCTIVITY)

    property Pr_fr:
        """
        Frozen Prandtl number at each shock state.

        Returns
        -------
        np.ndarray
            Frozen Prandtl number array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_FROZEN_PRANDTL)

    property Pr_eq:
        """
        Equilibrium Prandtl number at each shock state.

        Returns
        -------
        np.ndarray
            Equilibrium Prandtl number array [initial, incident, reflected (if applicable)]
        """
        def __get__(self):
            return self._get_property(SHOCK_EQUILIBRIUM_PRANDTL)

    property nj:
        """
        Species concentrations (kg-mol of species j per kg of mixture) at each index.

        Returns
        -------
        np.ndarray
            2D array of concentrations (index × species)
        """
        def __get__(self):
            return self._get_weights()

    property ln_nj:
        """
        Natural logarithm of gas phase species concentrations at each station (gasses only).
        NOTE: ln_nj /= log(nj) because nj is truncated for trace species, while ln_nj preserves small amounts.

        Returns
        -------
        np.ndarray
            2D array of log concentrations (index × gas species)
        """
        def __get__(self):
            return self._get_weights(log=True)

    property n:
        """
        Total moles at each station.

        Returns
        -------
        np.ndarray
            Total moles array
        """
        def __get__(self):
            cdef cea_err ierr
            cdef int num_pts = self._get_size()
            cdef cea_real *value = <cea_real *>malloc(num_pts * sizeof(cea_real))
            if value == NULL:
                raise MemoryError("Failed to allocate moles buffer")

            try:
                ierr = cea_shock_solution_get_moles(self.ptr, value)
                if ierr != SUCCESS:
                    self.last_error = <int>ierr
                _check_ierr(ierr, "ShockSolution.n")

                val = np.zeros(num_pts, dtype=np.double)
                for i in range(num_pts):
                    val[i] = value[i]
            finally:
                free(value)

            return val

    property mass_fractions:
        """
        Species mass fractions at each station.

        Returns
        -------
        dict
            Dictionary mapping species names to 2D arrays (index × species)
        """
        def __get__(self):
            species_names = self.solver.products.species_names
            vals = np.transpose(self._get_species_amounts(mass_fraction=True))
            return dict(zip(species_names, vals))

    property mole_fractions:
        """
        Species mole fractions at each station.

        Returns
        -------
        dict
            Dictionary mapping species names to 2D arrays (index × species)
        """
        def __get__(self):
            species_names = self.solver.products.species_names
            vals = np.transpose(self._get_species_amounts(mass_fraction=False))
            return dict(zip(species_names, vals))

    property converged:
        """
        Convergence status of the solution.

        Returns
        -------
        bool
            True if solution converged, False otherwise
        """
        def __get__(self):
            cdef cea_err ierr
            cdef bint value  # Use bint for boolean values
            ierr = cea_shock_solution_get_converged(self.ptr, <bint *>(&value))  # Explicit cast to bint *
            return bool(value)

    def _get_property(self, cea_shock_property_type prop_type):
        cdef cea_err ierr
        cdef int num_pts = self.num_pts
        cdef cea_real *value = <cea_real *>malloc(num_pts * sizeof(cea_real))
        if value == NULL:
            raise MemoryError("Failed to allocate property buffer")

        try:
            ierr = cea_shock_solution_get_property(self.ptr, prop_type, num_pts, value)
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "ShockSolution._get_property")

            prop_val = np.zeros(num_pts, dtype=np.double)
            for i in range(num_pts):
                prop_val[i] = value[i]
        finally:
            free(value)

        return prop_val

    def _get_scalar_property(self, cea_shock_property_type prop_type):
        cdef cea_err ierr
        cdef cea_real value

        ierr = cea_shock_solution_get_scalar_property(self.ptr, prop_type, &value)
        if ierr != SUCCESS:
            self.last_error = <int>ierr
        _check_ierr(ierr, "ShockSolution._get_scalar_property")

        return value

    def _get_weights(self, bint log=False):
        cdef cea_err ierr
        cdef int num_pts = self._get_size()
        cdef int nspecies

        if log:
            nspecies = self.solver.num_gas
        else:
            nspecies = self.solver.num_products
        cdef cea_real *reac_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if reac_weights == NULL:
            raise MemoryError("Failed to allocate weights buffer")

        try:
            if log:
                weights = np.zeros((num_pts, nspecies), dtype=np.double)
                for i in range(1,num_pts+1):
                    ierr = cea_shock_solution_get_weights(self.ptr, nspecies, <cea_int>i, reac_weights, log)
                    for j in range(nspecies):
                        weights[i-1, j] = reac_weights[j]
            else:
                weights = np.zeros((num_pts, nspecies), dtype=np.double)
                for i in range(1,num_pts+1):
                    ierr = cea_shock_solution_get_weights(self.ptr, nspecies, <cea_int>i, reac_weights, log)
                    for j in range(nspecies):
                        weights[i-1, j] = reac_weights[j]
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "ShockSolution._get_weights")
        finally:
            free(reac_weights)

        return weights

    def _get_species_amounts(self, bint mass_fraction=False):

        cdef cea_err ierr
        cdef int nspecies = self.solver.num_products

        cdef cea_real *cea_amounts = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if cea_amounts == NULL:
            raise MemoryError("Failed to allocate species amounts buffer")

        try:
            amounts = np.zeros((self.num_pts, nspecies), dtype=np.double)
            for i in range(1,self.num_pts+1):
                ierr = cea_shock_solution_get_species_amounts(self.ptr, nspecies, <cea_int>i, cea_amounts, mass_fraction)
                for j in range(nspecies):
                    amounts[i-1, j] = cea_amounts[j]
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "ShockSolution._get_species_amounts")
        finally:
            free(cea_amounts)

        return amounts


cdef class DetonationSolver:
    """
    Detonation solver for calculating Chapman-Jouguet detonation properties.

    Solves for detonation wave properties including detonation velocity,
    pressure, and temperature.

    Parameters
    ----------
    products : Mixture
        Mixture object containing product species
    **kwargs
        reactants : Mixture, optional
            Mixture object containing reactant species
        transport : bool, default False
            Enable transport property calculations
        ions : bool, default False
            Include ionized species in calculations
        trace : float, default -1.0
            Trace species threshold value; values < 0.0 uses default value
        insert : list of str, default []
            Additional species to insert into calculation; used to start initial guess with condensed species
    """
    cdef cea_detonation_solver ptr
    cdef Mixture products
    cdef object _keepalive_insert

    def __init__(self, Mixture products, **kwargs):

        # Extract keyword arguments with defaults
        cdef bint transport = kwargs.get('transport', False)
        cdef bint ions = kwargs.get('ions', False)
        cdef double trace_val = kwargs.get('trace', -1.0)
        cdef cea_string* cea_insert = NULL
        insert = kwargs.get('insert', [])
        cdef list _insert_keepalive = []
        cdef _CString _encoded

        cdef cea_err ierr
        cdef cea_mixture prods = products.ptr
        cdef cea_solver_opts opts
        cdef cea_mixture reacs = NULL
        cdef Mixture py_reactants
        self.products = products

        cea_solver_opts_init(&opts)

        if len(insert) > 0:
            cea_insert = <cea_string*>malloc(sizeof(cea_string)*len(insert))
            if cea_insert == NULL:
                raise MemoryError("Failed to allocate insert pointer buffer")
            for i, val in enumerate(insert):
                _encoded = _CString(val, "DetonationSolver insert species")
                _insert_keepalive.append(_encoded)
                cea_insert[i] = _encoded.ptr
        self._keepalive_insert = _insert_keepalive

        # Set options
        opts.transport = transport
        opts.ions = ions
        opts.trace = trace_val
        opts.ninsert = <int>len(insert)
        opts.insert = cea_insert

        if "reactants" in kwargs:
            reactants = kwargs["reactants"]
            if not isinstance(reactants, Mixture):
                raise TypeError("reactants must be a Mixture object")
            py_reactants = <Mixture>reactants
            opts.reactants = py_reactants.ptr

        try:
            ierr = cea_detonation_solver_create_with_options(&self.ptr, prods, opts)
        finally:
            if cea_insert != NULL:
                free(cea_insert)

        return

    def __dealloc__(self):
        if self.ptr:
            cea_detonation_solver_destroy(&self.ptr)
        return

    property num_reactants:
        """
        Number of reactant species.

        Returns
        -------
        int
            Number of reactants
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_detonation_solver_get_size(self.ptr, NUM_REACTANTS, &size_val)
            return size_val

    property num_products:
        """
        Number of product species.

        Returns
        -------
        int
            Number of products
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_detonation_solver_get_size(self.ptr, NUM_PRODUCTS, &size_val)
            return size_val

    property num_gas:
        """
        Number of gas phase species.

        Returns
        -------
        int
            Number of gas species
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_detonation_solver_get_size(self.ptr, NUM_GAS, &size_val)
            return size_val

    property num_condensed:
        """
        Number of condensed phase species.

        Returns
        -------
        int
            Number of condensed species
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_detonation_solver_get_size(self.ptr, NUM_CONDENSED, &size_val)
            return size_val

    property num_elements:
        """
        Number of chemical elements.

        Returns
        -------
        int
            Number of elements
        """
        def __get__(self):
            cdef cea_int size_val
            ierr = cea_detonation_solver_get_size(self.ptr, NUM_ELEMENTS, &size_val)
            return size_val

    def solve(self, DetonationSolution soln, np.ndarray weights, float T1, float p1,
              bint frozen=False):
        """
        Solve detonation problem.

        Parameters
        ----------
        soln : DetonationSolution
            Solution object to store results
        weights : np.ndarray
            Initial composition mass fractions
        T1 : float
            Initial temperature in K
        p1 : float
            Initial pressure in bar
        frozen : bool, default False
            Use frozen calculations for detonation
        """
        cdef cea_err ierr
        cdef int nweights = <int>len(weights)

        cdef cea_array wts = <cea_array>malloc(nweights * sizeof(double))
        if wts == NULL:
            raise MemoryError("Failed to allocate weights buffer")

        for i in range(nweights):
            wts[i] = weights[i]

        try:
            ierr = cea_detonation_solver_solve(self.ptr, soln.ptr, wts, <cea_real>T1, <cea_real>p1, frozen)
            soln.last_error = <int>ierr
            _check_ierr(ierr, "DetonationSolver.solve")
        finally:
            free(wts)

        return

cdef class DetonationSolution:
    """
    Solution object containing detonation calculation results.

    Stores Chapman-Jouguet detonation properties including velocity,
    pressure, temperature, and thermodynamic states.

    Parameters
    ----------
    solver : DetonationSolver
        Detonation solver instance
    """
    cdef cea_detonation_solution ptr
    cdef DetonationSolver solver
    cdef public int last_error

    def __cinit__(self, DetonationSolver solver):

        cdef cea_err ierr
        ierr = cea_detonation_solution_create(&self.ptr)
        _check_ierr(ierr, "DetonationSolution.__cinit__")
        self.solver = solver
        self.last_error = <int>SUCCESS

        return

    def __dealloc__(self):
        if self.ptr:
            cea_detonation_solution_destroy(&self.ptr)
        return

    property P1:
        """
        Initial pressure in bar.

        Returns
        -------
        float
            Pressure
        """
        def __get__(self):
            return self._get_property(DETONATION_P1)

    property T1:
        """
        Initial temperature in K.

        Returns
        -------
        float
            Temperature
        """
        def __get__(self):
            return self._get_property(DETONATION_T1)

    property H1:
        """
        Initial specific enthalpy in kJ/kg.

        Returns
        -------
        float
            Enthalpy
        """
        def __get__(self):
            return self._get_property(DETONATION_H1)

    property M1:
        """
        Initial molecular weight (1/n).

        Returns
        -------
        float
            Molecular weight
        """
        def __get__(self):
            return self._get_property(DETONATION_M1)

    property gamma1:
        """
        Initial isentropic exponent \gamma_1.

        Returns
        -------
        float
            Isentropic exponent
        """
        def __get__(self):
            return self._get_property(DETONATION_GAMMA1)

    property sonic_velocity1:
        """
        Initial sonic velocity in m/s.

        Returns
        -------
        float
            Sonic velocity
        """
        def __get__(self):
            return self._get_property(DETONATION_V_SONIC1)

    property P:
        """
        Detonation pressure in bar.

        Returns
        -------
        float
            Pressure of detonation products
        """
        def __get__(self):
            return self._get_property(DETONATION_PRESSURE)

    property T:
        """
        Detonation temperature in K.

        Returns
        -------
        float
            Temperature of detonation products
        """
        def __get__(self):
            return self._get_property(DETONATION_TEMPERATURE)

    property density:
        """
        Detonation density in kg/m³.

        Returns
        -------
        float
            Density of detonation products
        """
        def __get__(self):
            return self._get_property(DETONATION_DENSITY)

    property enthalpy:
        """
        Detonation specific enthalpy in kJ/kg.

        Returns
        -------
        float
            Specific enthalpy of detonation products
        """
        def __get__(self):
            return self._get_property(DETONATION_ENTHALPY)

    property energy:
        """
        Detonation specific internal energy in kJ/kg.

        Returns
        -------
        float
            Specific internal energy of detonation products
        """
        def __get__(self):
            return self._get_property(DETONATION_ENERGY)

    property gibbs_energy:
        """
        Detonation Gibbs energy in kJ/kg.

        Returns
        -------
        float
            Gibbs energy of detonation products
        """
        def __get__(self):
            return self._get_property(DETONATION_GIBBS_ENERGY)

    property entropy:
        """
        Detonation specific entropy in kJ/(kg·K).

        Returns
        -------
        float
            Specific entropy of detonation products
        """
        def __get__(self):
            return self._get_property(DETONATION_ENTROPY)

    property Mach:
        """
        Detonation Mach number.

        Returns
        -------
        float
            Mach number of the detonation wave
        """
        def __get__(self):
            return self._get_property(DETONATION_MACH)

    property velocity:
        """
        Detonation velocity in m/s.

        Returns
        -------
        float
            Velocity of the detonation wave
        """
        def __get__(self):
            return self._get_property(DETONATION_VELOCITY)

    property sonic_velocity:
        """
        Detonation sonic velocity in m/s.

        Returns
        -------
        float
            Sonic velocity of the detonation products
        """
        def __get__(self):
            return self._get_property(DETONATION_SONIC_VELOCITY)

    property gamma_s:
        """
        Detonation isentropic exponent \gamma_s.

        Returns
        -------
        float
            Isentropic exponent of the detonation products
        """
        def __get__(self):
            return self._get_property(DETONATION_GAMMA)

    property P_P1:
        """
        Pressure ratio P/P1.

        Returns
        -------
        float
            Pressure ratio of detonation products to initial pressure
        """
        def __get__(self):
            return self._get_property(DETONATION_P_P1)

    property T_T1:
        """
        Temperature ratio T/T1.

        Returns
        -------
        float
            Temperature ratio of detonation products to initial temperature
        """
        def __get__(self):
            return self._get_property(DETONATION_T_T1)

    property M_M1:
        """
        Molecular weight ratio M/M1.

        Returns
        -------
        float
            Molecular weight ratio of detonation products to initial molecular weight
        """
        def __get__(self):
            return self._get_property(DETONATION_M_M1)

    property rho_rho1:
        """
        Density ratio ρ/ρ1.

        Returns
        -------
        float
            Density ratio of detonation products to initial density
        """
        def __get__(self):
            return self._get_property(DETONATION_RHO_RHO1)

    property cp_fr:
        """
        Frozen specific heat at constant pressure in kJ/(kg·K).

        Returns
        -------
        float
            Frozen specific heat at constant pressure
        """
        def __get__(self):
            return self._get_property(DETONATION_FROZEN_CP)

    property cv_fr:
        """
        Frozen specific heat at constant volume in kJ/(kg·K).

        Returns
        -------
        float
            Frozen specific heat at constant volume
        """
        def __get__(self):
            return self._get_property(DETONATION_FROZEN_CV)

    property cp_eq:
        """
        Equilibrium specific heat at constant pressure in kJ/(kg·K).

        Returns
        -------
        float
            Equilibrium specific heat at constant pressure
        """
        def __get__(self):
            return self._get_property(DETONATION_EQUILIBRIUM_CP)

    property cv_eq:
        """
        Equilibrium specific heat at constant volume in kJ/(kg·K).

        Returns
        -------
        float
            Equilibrium specific heat at constant volume
        """
        def __get__(self):
            return self._get_property(DETONATION_EQUILIBRIUM_CV)

    property M:
        """
        Molecular weight (1/n) (gasses only) of detonation products in kg/kmol.

        Returns
        -------
        float
            Molecular weight of detonation products
        """
        def __get__(self):
            return self._get_property(DETONATION_M)

    property MW:
        """
        Molecular weight (including condensed species) of detonation products in kg/kmol.

        Returns
        -------
        float
            Molecular weight of detonation products including condensed species
        """
        def __get__(self):
            return self._get_property(DETONATION_MW)

    property viscosity:
        """
        Viscosity in micropoise (\mu P).

        Returns
        -------
        float
            Viscosity
        """
        def __get__(self):
            return self._get_property(DETONATION_VISCOSITY)

    property conductivity_fr:
        """
        Frozen thermal conductivity in \mu W/(cm·K).

        Returns
        -------
        float
            Frozen thermal conductivity
        """
        def __get__(self):
            return self._get_property(DETONATION_FROZEN_CONDUCTIVITY)

    property conductivity_eq:
        """
        Equilibrium thermal conductivity in \mu W/(cm·K).

        Returns
        -------
        float
            Equilibrium thermal conductivity
        """
        def __get__(self):
            return self._get_property(DETONATION_EQUILIBRIUM_CONDUCTIVITY)

    property Pr_fr:
        """
        Frozen Prandtl number.

        Returns
        -------
        float
            Frozen Prandtl number
        """
        def __get__(self):
            return self._get_property(DETONATION_FROZEN_PRANDTL)

    property Pr_eq:
        """
        Equilibrium Prandtl number.

        Returns
        -------
        float
            Equilibrium Prandtl number
        """
        def __get__(self):
            return self._get_property(DETONATION_EQUILIBRIUM_PRANDTL)

    property nj:
        """
        Species concentrations (kg-mol of species j per kg of mixture).

        Returns
        -------
        np.ndarray
            Concentrations for each species
        """
        def __get__(self):
            return self._get_weights()

    property ln_nj:
        """
        Natural logarithm of gas phase species concentrations (gasses only).
        NOTE: ln_nj /= log(nj) because nj is truncated for trace species, while ln_nj preserves small amounts.

        Returns
        -------
        np.ndarray
            Log concentrations for gas species
        """
        def __get__(self):
            return self._get_weights(log=True)

    property n:
        """
        Total moles.

        Returns
        -------
        float
            Total moles
        """
        def __get__(self):
            cdef cea_err ierr
            cdef cea_real value
            ierr = cea_detonation_solution_get_moles(self.ptr, &value)
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "DetonationSolution.n")
            return value

    property mass_fractions:
        """
        Species mass fractions.

        Returns
        -------
        dict
            Dictionary mapping species names to mass fractions
        """
        def __get__(self):
            species_names = self.solver.products.species_names
            vals = self._get_species_amounts(mass_fraction=True)
            return dict(zip(species_names, vals))

    property mole_fractions:
        """
        Species mole fractions.

        Returns
        -------
        dict
            Dictionary mapping species names to mole fractions
        """
        def __get__(self):
            species_names = self.solver.products.species_names
            vals = self._get_species_amounts(mass_fraction=False)
            return dict(zip(species_names, vals))

    property converged:
        """
        Convergence status of the solution.

        Returns
        -------
        bool
            True if solution converged, False otherwise
        """
        def __get__(self):
            cdef cea_err ierr
            cdef bint value  # Use bint for boolean values
            ierr = cea_detonation_solution_get_converged(self.ptr, <bint *>(&value))  # Explicit cast to bint *
            return bool(value)

    def _get_property(self, cea_detonation_property_type prop_type):
        cdef cea_err ierr
        cdef cea_real value

        ierr = cea_detonation_solution_get_property(self.ptr, prop_type, 1, &value)
        if ierr != SUCCESS:
            self.last_error = <int>ierr
        _check_ierr(ierr, "DetonationSolution._get_property")

        return value

    def _get_weights(self, bint log=False):
        cdef cea_err ierr
        cdef int nspecies

        if log:
            nspecies = self.solver.num_gas
        else:
            nspecies = self.solver.num_products
        cdef cea_real *reac_weights = <cea_real *>malloc(nspecies * sizeof(cea_real))
        if reac_weights == NULL:
            raise MemoryError("Failed to allocate weights buffer")

        try:
            if log:
                weights = np.zeros(nspecies, dtype=np.double)
                ierr = cea_detonation_solution_get_weights(self.ptr, nspecies, reac_weights, log)
            else:
                weights = np.zeros(nspecies, dtype=np.double)
                ierr = cea_detonation_solution_get_weights(self.ptr, nspecies, reac_weights, log)
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "DetonationSolution._get_weights")

            # Assign the weights to the numpy array
            for i in range(nspecies):
                weights[i] = reac_weights[i]
        finally:
            free(reac_weights)

        return weights

    def _get_species_amounts(self, bint mass_fraction=False):

        cdef cea_err ierr
        cdef int nspecies

        nprod = self.solver.num_products
        cdef cea_real *cea_amounts = <cea_real *>malloc(nprod * sizeof(cea_real))
        if cea_amounts == NULL:
            raise MemoryError("Failed to allocate species amounts buffer")

        try:
            amounts = np.zeros(nprod, dtype=np.double)
            ierr = cea_detonation_solution_get_species_amounts(self.ptr, nprod, cea_amounts, mass_fraction)
            if ierr != SUCCESS:
                self.last_error = <int>ierr
            _check_ierr(ierr, "DetonationSolution._get_species_amounts")

            # Assign the amounts to the numpy array
            for i in range(nprod):
                amounts[i] = cea_amounts[i]
        finally:
            free(cea_amounts)

        return amounts

# ----------------
# Matlab interface
# ----------------

def eq_solve(cea_equilibrium_type eq_type,
             list reactants,
             T: Optional[double] = None, H: Optional[double] = None, S: Optional[double] = None, U: Optional[double] = None,
             P: Optional[double] = None, V: Optional[double] = None,
             T_reac: Optional[list | np.ndarray | double] = None, fuel_amounts: Optional[list | np.ndarray] = None, oxid_amounts: Optional[list | np.ndarray] = None, bint moles=False,
             of_ratio: Optional[list | np.ndarray | double] = None, phi: Optional[list | np.ndarray | double] = None,
             r_eq: Optional[list | np.ndarray | double] = None, pct_fuel: Optional[list | np.ndarray | double] = None,
             only: Optional[list] = None, omit: Optional[list] = None, insert: Optional[list] = None,
             trace: Optional[double] = None, bint transport=False, bint ions=False):
    """
    Solve equilibrium problem using CEA methodology (Matlab-style interface).

    High-level interface function that handles mixture creation, solver setup,
    and solution calculation for equilibrium problems.

    Parameters
    ----------
    eq_type : int
        Equilibrium type constant (TP, HP, SP, TV, UV, SV)
    reactants : list of str
        List of reactant species names
    T : float, optional
        Temperature value in K (mutually exclusive with H, S, U)
    H : float, optional
        Enthalpy value in J/kg (mutually exclusive with T, S, U)
    S : float, optional
        Entropy value in J/(kg·K) (mutually exclusive with T, H, U)
    U : float, optional
        Internal energy value in J/kg (mutually exclusive with T, H, S)
    P : float, optional
        Pressure value in bar (mutually exclusive with V)
    V : float, optional
        Specific volume value in m³/kg (mutually exclusive with P)
    T_reac : float, list, or np.ndarray, optional
        Reactant temperature(s) in K. If list/array, must match reactants length
    fuel_amounts : list or np.ndarray, optional
        Fuel species amounts (mass or mole fractions)
    oxid_amounts : list or np.ndarray, optional
        Oxidizer species amounts (mass or mole fractions)
    moles : bool, default False
        Whether fuel_amounts and oxid_amounts are in moles (True) or mass (False)
    of_ratio : float, list, or np.ndarray, optional
        Oxidizer/fuel mass ratio
    phi : float, list, or np.ndarray, optional
        Weight equivalence ratio
    r_eq : float, list, or np.ndarray, optional
        Chemical equivalence ratio
    pct_fuel : float, list, or np.ndarray, optional
        Fuel percentage by mass
    only : list of str, optional
        Restrict products to only these species
    omit : list of str, optional
        Omit these species from products
    insert : list of str, optional
        Insert additional species into calculation
    trace : float, optional
        Trace species threshold value
    transport : bool, default False
        Enable transport property calculations
    ions : bool, default False
        Include ionic species in calculations

    Returns
    -------
    EqSolution
        Solution object containing equilibrium results

    Raises
    ------
    TypeError
        If required parameters missing or fuel/oxidizer amounts not defined
    ValueError
        If T_reac length doesn't match reactants, or invalid parameter combinations
    """

    cdef cea_err ierr

    omit_list = [] if omit is None else omit
    insert_list = [] if insert is None else insert
    only_list = [] if only is None else only

    # Create the reactants mixture from input reactants (following main.f90 pattern)
    reactants_mix = Mixture(reactants, ions=ions)

    # Get the products Mixture object (following main.f90 pattern)
    if len(only_list) > 0:
        # If only specified, use those species
        products_mix = Mixture(only_list, ions=ions)
    else:
        products_mix = Mixture(reactants, products_from_reactants=True, omit=omit_list, ions=ions)

    # Check thermodynamic state values
    state1 = 0.0
    needs_state1 = (T is None) and (H is None) and (S is None) and (U is None)
    if (T is not None):
        state1 = T
    elif (H is not None):
        state1 = H
    elif (S is not None):
        state1 = S
    elif (U is not None):
        state1 = U

    if (P is not None):
        state2 = P
    elif (V is not None):
        state2 = V
    else:
        raise TypeError("Pressure or volume state not defined")

    # Check that the reactant weights are defined somehow
    weights = np.zeros(len(reactants))
    if (fuel_amounts is not None) and (oxid_amounts is not None):
        # Convert fuel/oxidant amounts to total weights
        fuel_amounts = np.asarray(fuel_amounts, dtype=np.double)
        oxid_amounts = np.asarray(oxid_amounts, dtype=np.double)
        if moles:
            fuel_weights = reactants_mix.moles_to_weights(fuel_amounts)
            oxid_weights = reactants_mix.moles_to_weights(oxid_amounts)
        else:
            fuel_weights = fuel_amounts
            oxid_weights = oxid_amounts

        if of_ratio is not None:
            weights = reactants_mix.of_ratio_to_weights(oxid_weights, fuel_weights, of_ratio)
        elif phi is not None:
            of_ratio_calc = reactants_mix.weight_eq_ratio_to_of_ratio(oxid_weights, fuel_weights, phi)
            weights = reactants_mix.of_ratio_to_weights(oxid_weights, fuel_weights, of_ratio_calc)
        elif r_eq is not None:
            of_ratio_calc = reactants_mix.chem_eq_ratio_to_of_ratio(oxid_weights, fuel_weights, r_eq)
            weights = reactants_mix.of_ratio_to_weights(oxid_weights, fuel_weights, of_ratio_calc)
        elif pct_fuel is not None:
            of_ratio_calc = (100.0 - pct_fuel) / pct_fuel
            weights = reactants_mix.of_ratio_to_weights(oxid_weights, fuel_weights, of_ratio_calc)
        else:
            # If no ratio specified, add the oxidant and fuel weights
            weights = fuel_weights + oxid_weights
    else:
        raise TypeError("Fuel and oxidizer amounts not defined")

    # Compute the fixed thermodynamic state if needed
    if needs_state1:
        if T_reac is not None:
            if type(T_reac) in [list, np.ndarray]:
                if len(T_reac) == 0:
                    raise ValueError("T_reac cannot be empty")
                if len(T_reac) != len(reactants):
                    raise ValueError("T_reac must have the same length as reactants")
            elif not isinstance(T_reac, (float, int)):
                raise ValueError("T_reac must be a float, list, or np.ndarray")

            if eq_type == TP or eq_type == TV:
                state1 = T_reac[0] if type(T_reac) in [list, np.ndarray] else T_reac
                warnings.warn("Problem temperature not defined; using first reactant temperature")
            elif eq_type == HP:
                state1 = reactants_mix.calc_property(ENTHALPY, weights, T_reac)
            elif eq_type == SP or eq_type == SV:
                state1 = reactants_mix.calc_property(ENTROPY, weights, T_reac)
            elif eq_type == UV:
                state1 = reactants_mix.calc_property(ENERGY, weights, T_reac)
        else:
            raise TypeError("Reactant temperature not defined")

    # Initialize the solver (following main.f90 pattern)
    solver_kwargs = {
        "reactants": reactants_mix,
        "transport": transport,
        "ions": ions,
        "insert": insert_list,
    }
    if trace is not None:
        solver_kwargs["trace"] = trace
    solver = EqSolver(products_mix, **solver_kwargs)

    # Initialize the solution
    soln = EqSolution(solver)

    # Call the solver
    solver.solve(soln, eq_type, state1, state2, weights)

    return soln


# Public API: export non-underscore names, excluding helper modules/types.
_public_exclude = {"cython", "ctypes", "array", "warnings", "np", "Optional"}
__all__ = [name for name in globals() if not name.startswith("_") and name not in _public_exclude]
