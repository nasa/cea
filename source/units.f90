module cea_units
    !! Unit conversion helper functions

    use cea_param
    use fb_logging
    implicit none

contains

    ! Temperature conversions
    ! -----------------------

    function celcius_to_kelvin(celcius) result(kelvin)
        real(dp), intent(in) :: celcius
        real(dp) :: kelvin
        kelvin = celcius + 273.15d0
    end function

    function fahrenheit_to_kelvin(fahrenheit) result(celcius)
        real(dp), intent(in) :: fahrenheit
        real(dp) :: celcius
        celcius = (fahrenheit - 32.0d0)/1.8d0 + 273.15d0
    end function

    function rankine_to_kelvin(rankine) result(kelvin)
        real(dp), intent(in) :: rankine
        real(dp) :: kelvin
        kelvin = rankine/1.8d0
    end function

    ! Pressure conversions
    ! --------------------

    function atm_to_bar(atm) result(bar)
        real(dp), intent(in) :: atm
        real(dp) :: bar
        bar = atm*1.01325d0
    end function

    function psi_to_bar(psi) result(bar)
        real(dp), intent(in) :: psi
        real(dp) :: bar
        bar = psi*0.0689475729d0
    end function

    function mmhg_to_bar(mmhg) result(bar)
        real(dp), intent(in) :: mmhg
        real(dp) :: bar
        bar = mmhg*0.001333d0
    end function

    ! Energy conversions
    ! ------------------

    function kilojoules_to_joules(kilojoules) result(joules)
        real(dp), intent(in) :: kilojoules
        real(dp) :: joules
        joules = kilojoules*1.0d3
    end function

    function kilocal_to_joules(kilocal) result(joules)
        real(dp), intent(in) :: kilocal
        real(dp) :: joules
        joules = kilocal*4.184d3
    end function

    function cal_to_joules(cal) result(joules)
        real(dp), intent(in) :: cal
        real(dp) :: joules
        joules = cal*4.184d0
    end function

    function joules_to_cal(joules) result(cal)
        real(dp), intent(in) :: joules
        real(dp) :: cal
        cal = joules/4.184d0
    end function

    ! Density conversions
    ! -------------------

    function g_per_cm3_to_kg_per_m3(g_per_cm3) result(kg_per_m3)
        real(dp), intent(in) :: g_per_cm3
        real(dp) :: kg_per_m3
        kg_per_m3 = g_per_cm3*1.0d3
    end function

    ! Volume conversions
    ! ------------------

    function cm3_per_g_to_m3_per_kg(cm3_per_g) result(m3_per_kg)
        real(dp), intent(in) :: cm3_per_g
        real(dp) :: m3_per_kg
        m3_per_kg = cm3_per_g*1.0d-3
    end function

    ! Converter helper function
    ! -------------------------
    function convert_units_to_si(val, units) result(si_val)

        ! Arguments
        real(dp), intent(in) :: val
        character(*), intent(in) :: units

        ! Results
        real(dp) :: si_val

        select case(units)

            ! Temperature
            case('k')
                si_val = val
            case('c')
                si_val = celcius_to_kelvin(val)
            case('f')
                si_val = fahrenheit_to_kelvin(val)
            case('r')
                si_val = rankine_to_kelvin(val)

            ! Pressure
            case('bar')
                si_val = val
            case('atm')
                si_val = atm_to_bar(val)
            case('psi', 'psia')
                si_val = psi_to_bar(val)
            case('mmhg')
                si_val = mmhg_to_bar(val)

            ! Energy
            case('j/mole')
                si_val = val
            case('kj/mole')
                si_val = kilojoules_to_joules(val)
            case('kcal/mole')
                si_val = kilocal_to_joules(val)
            case('cal/mole')
                si_val = cal_to_joules(val)

            ! Density
            case('kg/m**3')
                si_val = val
            case('g/cm**3', 'g/cc')
                si_val = g_per_cm3_to_kg_per_m3(val)

            ! Volume
            case('m**3/kg')
                si_val = val
            case('cm**3/g', 'cc/g')
                si_val = cm3_per_g_to_m3_per_kg(val)

            ! Unknown units
            case default
                call log_error('Unknown units in unit conversion: '//trim(units))
                si_val = val

        end select

    end function

end module