#define CEA_LOG_LEVEL_ENUM \
    CEA_LOG_CRITICAL = 50, \
    CEA_LOG_ERROR    = 40, \
    CEA_LOG_WARNING  = 30, \
    CEA_LOG_INFO     = 20, \
    CEA_LOG_DEBUG    = 10, \
    CEA_LOG_NONE     = 0

#define CEA_EQUILIBRIUM_TYPE_ENUM \
    CEA_TP = 0, \
    CEA_HP = 1, \
    CEA_SP = 2, \
    CEA_TV = 3, \
    CEA_UV = 4, \
    CEA_SV = 5

#define CEA_DERIVATIVE_METHOD_ENUM \
    CEA_DERIV_ANALYTIC = 0, \
    CEA_DERIV_FD = 1

#define CEA_EQDERIV_SCALAR_ENUM \
    CEA_DERIV_DT_DSTATE1 = 0, \
    CEA_DERIV_DT_DSTATE2 = 1, \
    CEA_DERIV_DN_DSTATE1 = 2, \
    CEA_DERIV_DN_DSTATE2 = 3, \
    CEA_DERIV_DH_DSTATE1 = 4, \
    CEA_DERIV_DH_DSTATE2 = 5, \
    CEA_DERIV_DU_DSTATE1 = 6, \
    CEA_DERIV_DU_DSTATE2 = 7, \
    CEA_DERIV_DG_DSTATE1 = 8, \
    CEA_DERIV_DG_DSTATE2 = 9, \
    CEA_DERIV_DS_DSTATE1 = 10, \
    CEA_DERIV_DS_DSTATE2 = 11

#define CEA_EQDERIV_ARRAY_ENUM \
    CEA_DERIV_DT_DW0 = 0, \
    CEA_DERIV_DN_DW0 = 1, \
    CEA_DERIV_DNJ_DSTATE1 = 2, \
    CEA_DERIV_DNJ_DSTATE2 = 3, \
    CEA_DERIV_DH_DW0 = 4, \
    CEA_DERIV_DU_DW0 = 5, \
    CEA_DERIV_DG_DW0 = 6, \
    CEA_DERIV_DS_DW0 = 7

#define CEA_EQDERIV_MATRIX_ENUM \
    CEA_DERIV_DNJ_DW0 = 0

#define CEA_FUEL_RATIO_TYPE_ENUM \
    CEA_RATIO_NONE = 0, \
    CEA_OF_RATIO = 1, \
    CEA_FUEL_PERCENT = 2, \
    CEA_PHI = 3, \
    CEA_R_EQ = 4

#define CEA_EQUILIBRIUM_SIZE_ENUM \
    CEA_NUM_REACTANTS = 0, \
    CEA_NUM_PRODUCTS  = 1, \
    CEA_NUM_GAS       = 2, \
    CEA_NUM_CONDENSED = 3, \
    CEA_NUM_ELEMENTS  = 4, \
    CEA_MAX_EQUATIONS = 5

#define CEA_PROPERTY_TYPE_ENUM \
    CEA_TEMPERATURE              = 0, \
    CEA_PRESSURE                 = 1, \
    CEA_VOLUME                   = 2, \
    CEA_DENSITY                  = 3, \
    CEA_M                        = 4, \
    CEA_MW                       = 5, \
    CEA_ENTHALPY                 = 6, \
    CEA_ENERGY                   = 7, \
    CEA_ENTROPY                  = 8, \
    CEA_GIBBS_ENERGY             = 9, \
    CEA_GAMMA_S                  = 10, \
    CEA_FROZEN_CP                = 11, \
    CEA_FROZEN_CV                = 12, \
    CEA_EQUILIBRIUM_CP           = 13, \
    CEA_EQUILIBRIUM_CV           = 14, \
    CEA_VISCOSITY                = 15, \
    CEA_FROZEN_CONDUCTIVITY      = 16, \
    CEA_EQUILIBRIUM_CONDUCTIVITY = 17, \
    CEA_FROZEN_PRANDTL           = 18, \
    CEA_EQUILIBRIUM_PRANDTL      = 19

#define CEA_ROCKET_PROPERTY_TYPE_ENUM \
    CEA_ROCKET_TEMPERATURE              = 0, \
    CEA_ROCKET_PRESSURE                 = 1, \
    CEA_ROCKET_VOLUME                   = 2, \
    CEA_ROCKET_DENSITY                  = 3, \
    CEA_ROCKET_M                        = 4, \
    CEA_ROCKET_MW                       = 5, \
    CEA_ROCKET_ENTHALPY                 = 6, \
    CEA_ROCKET_ENERGY                   = 7, \
    CEA_ROCKET_ENTROPY                  = 8, \
    CEA_ROCKET_GIBBS_ENERGY             = 9, \
    CEA_ROCKET_GAMMA_S                  = 10, \
    CEA_ROCKET_FROZEN_CP                = 11, \
    CEA_ROCKET_FROZEN_CV                = 12, \
    CEA_ROCKET_EQUILIBRIUM_CP           = 13, \
    CEA_ROCKET_EQUILIBRIUM_CV           = 14, \
    CEA_MACH                            = 15, \
    CEA_SONIC_VELOCITY                  = 16, \
    CEA_AE_AT                           = 17, \
    CEA_C_STAR                          = 18, \
    CEA_COEFFICIENT_OF_THRUST           = 19, \
    CEA_ISP                             = 20, \
    CEA_ISP_VACUUM                      = 21, \
    CEA_ROCKET_VISCOSITY                = 22, \
    CEA_ROCKET_FROZEN_CONDUCTIVITY      = 23, \
    CEA_ROCKET_EQUILIBRIUM_CONDUCTIVITY = 24, \
    CEA_ROCKET_FROZEN_PRANDTL           = 25, \
    CEA_ROCKET_EQUILIBRIUM_PRANDTL      = 26

#define CEA_SHOCK_PROPERTY_TYPE_ENUM \
    CEA_SHOCK_TEMPERATURE              = 0, \
    CEA_SHOCK_PRESSURE                 = 1, \
    CEA_SHOCK_VELOCITY                 = 2, \
    CEA_SHOCK_MACH                     = 3, \
    CEA_SHOCK_SONIC_VELOCITY           = 4, \
    CEA_SHOCK_RHO12                    = 5, \
    CEA_SHOCK_RHO52                    = 6, \
    CEA_SHOCK_P21                      = 7, \
    CEA_SHOCK_P52                      = 8, \
    CEA_SHOCK_T21                      = 9, \
    CEA_SHOCK_T52                      = 10, \
    CEA_SHOCK_M21                      = 11, \
    CEA_SHOCK_M52                      = 12, \
    CEA_SHOCK_V2                       = 13, \
    CEA_SHOCK_U5_P_V2                  = 14, \
    CEA_SHOCK_VOLUME                   = 15, \
    CEA_SHOCK_DENSITY                  = 16, \
    CEA_SHOCK_M                        = 17, \
    CEA_SHOCK_MW                       = 18, \
    CEA_SHOCK_ENTHALPY                 = 19, \
    CEA_SHOCK_ENERGY                   = 20, \
    CEA_SHOCK_ENTROPY                  = 21, \
    CEA_SHOCK_GIBBS_ENERGY             = 22, \
    CEA_SHOCK_GAMMA_S                  = 23, \
    CEA_SHOCK_FROZEN_CP                = 24, \
    CEA_SHOCK_FROZEN_CV                = 25, \
    CEA_SHOCK_EQUILIBRIUM_CP           = 26, \
    CEA_SHOCK_EQUILIBRIUM_CV           = 27, \
    CEA_SHOCK_VISCOSITY                = 28, \
    CEA_SHOCK_FROZEN_CONDUCTIVITY      = 29, \
    CEA_SHOCK_EQUILIBRIUM_CONDUCTIVITY = 30, \
    CEA_SHOCK_FROZEN_PRANDTL           = 31, \
    CEA_SHOCK_EQUILIBRIUM_PRANDTL      = 32

#define CEA_DETONATION_PROPERTY_TYPE_ENUM \
    CEA_DETONATION_P1                       = 0, \
    CEA_DETONATION_T1                       = 1, \
    CEA_DETONATION_H1                       = 2, \
    CEA_DETONATION_M1                       = 3, \
    CEA_DETONATION_GAMMA1                   = 4, \
    CEA_DETONATION_V_SONIC1                 = 5, \
    CEA_DETONATION_PRESSURE                 = 6, \
    CEA_DETONATION_TEMPERATURE              = 7, \
    CEA_DETONATION_DENSITY                  = 8, \
    CEA_DETONATION_ENTHALPY                 = 9, \
    CEA_DETONATION_ENERGY                   = 10, \
    CEA_DETONATION_GIBBS_ENERGY             = 11, \
    CEA_DETONATION_ENTROPY                  = 12, \
    CEA_DETONATION_MACH                     = 13, \
    CEA_DETONATION_VELOCITY                 = 14, \
    CEA_DETONATION_SONIC_VELOCITY           = 15, \
    CEA_DETONATION_GAMMA                    = 16, \
    CEA_DETONATION_P_P1                     = 17, \
    CEA_DETONATION_T_T1                     = 18, \
    CEA_DETONATION_M_M1                     = 19, \
    CEA_DETONATION_RHO_RHO1                 = 20, \
    CEA_DETONATION_FROZEN_CP                = 21, \
    CEA_DETONATION_FROZEN_CV                = 22, \
    CEA_DETONATION_EQUILIBRIUM_CP           = 23, \
    CEA_DETONATION_EQUILIBRIUM_CV           = 24, \
    CEA_DETONATION_M                        = 25, \
    CEA_DETONATION_MW                       = 26, \
    CEA_DETONATION_VISCOSITY                = 27, \
    CEA_DETONATION_FROZEN_CONDUCTIVITY      = 28, \
    CEA_DETONATION_EQUILIBRIUM_CONDUCTIVITY = 29, \
    CEA_DETONATION_FROZEN_PRANDTL           = 30, \
    CEA_DETONATION_EQUILIBRIUM_PRANDTL      = 31

#define CEA_ERROR_CODE_ENUM \
    CEA_SUCCESS                  = 0, \
    CEA_INVALID_FILENAME         = 1, \
    CEA_INVALID_PROPERTY_TYPE    = 2, \
    CEA_INVALID_EQUILIBRIUM_TYPE = 3, \
    CEA_INVALID_ROCKET_TYPE      = 4, \
    CEA_INVALID_EQUILIBRIUM_SIZE_TYPE = 5, \
    CEA_INVALID_INDEX            = 6, \
    CEA_INVALID_SIZE             = 7, \
    CEA_NOT_CONVERGED            = 8
