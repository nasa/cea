module cea_fits

    use cea_param, only: dp
    implicit none

    type :: ThermoFit
        !! Container for thermodynamic curve-fit coeffients

        real(dp) :: a1
        real(dp) :: a2
        real(dp) :: a3
        real(dp) :: a4
        real(dp) :: a5
        real(dp) :: a6
        real(dp) :: a7
        real(dp) :: b1
        real(dp) :: b2
    contains
        procedure :: calc_cv
        procedure :: calc_cp
        procedure :: calc_energy
        procedure :: calc_enthalpy
        procedure :: calc_entropy
        procedure :: calc_gibbs_energy
        procedure :: calc_dcv_dT
        procedure :: calc_dcp_dT
        procedure :: calc_denergy_dT
        procedure :: calc_denthalpy_dT
        procedure :: calc_dentropy_dT
    end type

    type :: TransportFit
        !! Container for transport curve-fit coefficients

        real(dp) :: A
        real(dp) :: B
        real(dp) :: C
        real(dp) :: D
    contains
        procedure :: calc_transport_value
    end type

contains

    elemental function calc_cv(self,T) result(cv)
        class(ThermoFit), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: cv
        cv = self%calc_cp(T) - 1.0d0
    end function

    elemental function calc_cp(self,T) result(cp)
        class(ThermoFit), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: cp
        cp = self%a7
        cp = T*cp + self%a6
        cp = T*cp + self%a5
        cp = T*cp + self%a4
        cp = T*cp + self%a3
        cp = T*cp + self%a2
        cp = T*cp + self%a1
        cp = cp/(T*T)
    end function

    elemental function calc_energy(self,T,logT) result(e)
        class(ThermoFit), intent(in) :: self
        real(dp), intent(in) :: T, logT
        real(dp) :: e
        e = self%calc_enthalpy(T,logT) - T
    end function

    elemental function calc_enthalpy(self,T,logT) result(h)
        class(ThermoFit), intent(in) :: self
        real(dp), intent(in) :: T, logT
        real(dp) :: h  ! h/R
        h = self%a7/5.0d0
        h = T*h + self%a6/4.0d0
        h = T*h + self%a5/3.0d0
        h = T*h + self%a4/2.0d0
        h = T*h + self%a3
        h = T*h + self%b1
        h = T*h - self%a1
        h = h/T + self%a2*logT
    end function

    elemental function calc_entropy(self,T,logT) result(s)
        class(ThermoFit), intent(in) :: self
        real(dp), intent(in) :: T, logT
        real(dp) :: s
        s = self%a7/4.0d0
        s = T*s + self%a6/3.0d0
        s = T*s + self%a5/2.0d0
        s = T*s + self%a4
        s = T*s + self%b2
        s = T*s - self%a2
        s = T*s - self%a1/2.0d0
        s = s/(T*T) + self%a3*logT
    end function

    elemental function calc_gibbs_energy(self,T,logT) result(g)
        class(ThermoFit), intent(in) :: self
        real(dp), intent(in) :: T, logT
        real(dp) :: g
        g = -self%a7/20.0d0
        g = T*g - self%a6/12.0d0
        g = T*g - self%a5/6.0d0
        g = T*g - self%a4/2.0d0
        g = T*g + self%a3 - self%b2
        g = T*g + self%a2 + self%b1
        g = T*g - self%a1/2.0d0
        g = g/T + (self%a2 - T*self%a3)*logT
    end function

    elemental function calc_dcv_dT(self,T) result(cv)
        class(ThermoFit), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: cv
        cv = self%calc_cp(T) - 1.0d0
    end function

    elemental function calc_dcp_dT(self,T) result(dcp_dT)
        class(ThermoFit), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: dcp_dT
        dcp_dT = 4.0d0 * self%a7
        dcp_dT = T*dcp_dT + 3.0d0 * self%a6
        dcp_dT = T*dcp_dT + 2.0d0 * self%a5
        dcp_dT = T*dcp_dT + self%a4
        dcp_dT = T*dcp_dT
        dcp_dT = T*dcp_dT - self%a2
        dcp_dT = T*dcp_dT - 2.0d0 * self%a1
        dcp_dT = dcp_dT/(T*T*T)
    end function

    elemental function calc_denergy_dT(self,T) result(de_dT)
        class(ThermoFit), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: de_dT
        de_dT = self%calc_denthalpy_dT(T) - 1.0d0
    end function

    elemental function calc_denthalpy_dT(self,T) result(dh_dT)
        class(ThermoFit), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: dh_dT  ! h/R
        dh_dT = self%a7
        dh_dT = T*dh_dT + self%a6
        dh_dT = T*dh_dT + self%a5
        dh_dT = T*dh_dT + self%a4
        dh_dT = T*dh_dT + self%a3
        dh_dT = T*dh_dT + self%a2
        dh_dT = T*dh_dT + 2.0d0*self%a1
        dh_dT = dh_dT/(T*T)
    end function

    elemental function calc_dentropy_dT(self,T) result(ds_dT)
        class(ThermoFit), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: ds_dT
        ds_dT = self%a7
        ds_dT = T*ds_dT + self%a6
        ds_dT = T*ds_dT + self%a5
        ds_dT = T*ds_dT + self%a4
        ds_dT = T*ds_dT + self%a3
        ds_dT = T*ds_dT + self%a2
        ds_dT = T*ds_dT + self%a1
        ds_dT = ds_dT/(T*T*T)
    end function

    elemental function calc_transport_value(self,T) result(value)
        class(TransportFit), intent(in) :: self
        real(dp), intent(in) :: T
        real(dp) :: value
        value = self%A*log(T) + self%B/T + self%C/(T*T) + self%D
    end function

end module