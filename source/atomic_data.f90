module cea_atomic_data
    !! Module to store atomic data (atomic weights and valences)

    use cea_param, enl => element_name_len
    use fb_utils, only: abort, is_empty
    implicit none

    type Atomic
        !! Atomic data object

        character(enl) :: name
            !! Atomic symbol
        real(dp) :: atomic_weight = empty_dp
            !! Atomic weight
        real(dp) :: valence = empty_dp
            !! Atom valence

    end type

    type(Atomic), parameter :: AtomicData(103) = [ &
        !! Standard atomic weights are IUPAC/CIAAW's current values (accessed
        !! 2026-08-24): https://www.ciaaw.org/abridged-atomic-weights.htm
        !! Deliberately kept at CIAAW's own 5-significant-digit abridged
        !! precision throughout (their recommendation since 2013) rather than
        !! the more precise unabridged table -- see the PR for why.
        !!
        !! Interval-notation elements (e.g. Li, H, Cl) reflect real natural
        !! variation, not measurement uncertainty -- don't "fix" one by averaging
        !! it. The value here is IUPAC's own conventional recommendation, not
        !! the interval's midpoint.
        !!
        !! Elements with no standard atomic weight (plus D) use CIAAW's
        !! recommended isotope's mass instead
        !! (https://www.ciaaw.org/radioactive-elements.htm), at AME2020's full
        !! published precision (https://www-nds.iaea.org/amdc/ame2020/mass_1.mas20.txt)
        !! -- unlike standard atomic weights, isotope masses have no abridged
        !! form to defer to, so full precision is used rather than an invented
        !! shorter one; see the PR for why.
        Atomic("H ", 1.008d0     ,  1.), &
        Atomic("IH", 1.008d0     ,  1.), &   ! Inert hydrogen
        Atomic("D ", 2.014101777844d0 ,  1.), &   ! Deuterium atomic mass (not a standard atomic weight); ciaaw.org/hydrogen.htm
        Atomic("HE", 4.0026d0    ,  0.), &
        Atomic("LI", 6.94d0      ,  1.), &
        Atomic("BE", 9.0122d0    ,  2.), &
        Atomic("B ", 10.81d0     ,  3.), &
        Atomic("C ", 12.011d0    ,  4.), &
        Atomic("IC", 12.011d0    ,  4.), &   ! Inert carbon
        Atomic("N ", 14.007d0    ,  0.), &
        Atomic("O ", 15.999d0    , -2.), &
        Atomic("IO", 15.999d0    , -2.), &   ! Inert oxygen
        Atomic("F ", 18.998d0    , -1.), &
        Atomic("NE", 20.180d0    ,  0.), &
        Atomic("NA", 22.990d0    ,  1.), &
        Atomic("MG", 24.305d0    ,  2.), &
        Atomic("AL", 26.982d0    ,  3.), &
        Atomic("SI", 28.085d0    ,  4.), &
        Atomic("P ", 30.974d0    ,  5.), &
        Atomic("S ", 32.06d0     ,  4.), &
        Atomic("CL", 35.45d0     , -1.), &
        Atomic("AR", 39.95d0     ,  0.), &
        Atomic("K ", 39.098d0    ,  1.), &
        Atomic("CA", 40.078d0    ,  2.), &
        Atomic("SC", 44.956d0    ,  3.), &
        Atomic("TI", 47.867d0    ,  4.), &
        Atomic("V ", 50.942d0    ,  5.), &
        Atomic("CR", 51.996d0    ,  3.), &
        Atomic("MN", 54.938d0    ,  2.), &
        Atomic("FE", 55.845d0    ,  3.), &
        Atomic("CO", 58.933d0    ,  2.), &
        Atomic("NI", 58.693d0    ,  2.), &
        Atomic("CU", 63.546d0    ,  2.), &
        Atomic("ZN", 65.38d0     ,  2.), &
        Atomic("GA", 69.723d0    ,  3.), &
        Atomic("GE", 72.630d0    ,  4.), &
        Atomic("AS", 74.922d0    ,  3.), &
        Atomic("SE", 78.971d0    ,  4.), &
        Atomic("BR", 79.904d0    , -1.), &
        Atomic("KR", 83.798d0    ,  0.), &
        Atomic("RB", 85.468d0    ,  1.), &
        Atomic("SR", 87.62d0     ,  2.), &
        Atomic("Y ", 88.906d0    ,  3.), &
        Atomic("ZR", 91.222d0    ,  4.), &
        Atomic("NB", 92.906d0    ,  5.), &
        Atomic("MO", 95.95d0     ,  6.), &
        Atomic("TC", 97.907211206d0 ,  7.), &   ! No standard weight; mass of Tc-98 (CIAAW: Tc-97/Tc-98 equally stable)
        Atomic("RU", 101.07d0    ,  3.), &
        Atomic("RH", 102.91d0    ,  3.), &
        Atomic("PD", 106.42d0    ,  2.), &
        Atomic("AG", 107.87d0    ,  1.), &
        Atomic("CD", 112.41d0    ,  2.), &
        Atomic("IN", 114.82d0    ,  3.), &
        Atomic("SN", 118.71d0    ,  4.), &
        Atomic("SB", 121.76d0    ,  3.), &
        Atomic("TE", 127.6d0     ,  4.), &
        Atomic("I ", 126.90d0    , -1.), &
        Atomic("XE", 131.29d0    ,  0.), &
        Atomic("CS", 132.91d0    ,  1.), &
        Atomic("BA", 137.33d0    ,  2.), &
        Atomic("LA", 138.91d0    ,  3.), &
        Atomic("CE", 140.12d0    ,  3.), &
        Atomic("PR", 140.91d0    ,  3.), &
        Atomic("ND", 144.24d0    ,  3.), &
        Atomic("PM", 144.912755748d0 ,  3.), &   ! No standard weight; mass of Pm-145
        Atomic("SM", 150.36d0    ,  3.), &
        Atomic("EU", 151.96d0    ,  3.), &
        Atomic("GD", 157.25d0    ,  3.), &
        Atomic("TB", 158.93d0    ,  3.), &
        Atomic("DY", 162.50d0    ,  3.), &
        Atomic("HO", 164.93d0    ,  3.), &
        Atomic("ER", 167.26d0    ,  3.), &
        Atomic("TM", 168.93d0    ,  3.), &
        Atomic("YB", 173.05d0    ,  3.), &
        Atomic("LU", 174.97d0    ,  3.), &
        Atomic("HF", 178.49d0    ,  4.), &
        Atomic("TA", 180.95d0    ,  5.), &
        Atomic("W ", 183.84d0    ,  6.), &
        Atomic("RE", 186.21d0    ,  7.), &
        Atomic("OS", 190.23d0    ,  4.), &
        Atomic("IR", 192.22d0    ,  4.), &
        Atomic("PT", 195.08d0    ,  4.), &
        Atomic("AU", 196.97d0    ,  3.), &
        Atomic("HG", 200.59d0    ,  2.), &
        Atomic("TL", 204.38d0    ,  1.), &
        Atomic("PB", 207.2d0     ,  2.), &
        Atomic("BI", 208.98d0    ,  3.), &
        Atomic("PO", 208.982430361d0 ,  2.), &   ! No standard weight; mass of Po-209
        Atomic("AT", 209.987147423d0 , -1.), &   ! No standard weight; mass of At-210
        Atomic("RN", 222.017576017d0 ,  0.), &   ! No standard weight; mass of Rn-222
        Atomic("FR", 223.019734241d0 ,  1.), &   ! No standard weight; mass of Fr-223
        Atomic("RA", 226.025408186d0 ,  2.), &   ! No standard weight; mass of Ra-226
        Atomic("AC", 227.027750594d0 ,  3.), &   ! No standard weight; mass of Ac-227
        Atomic("TH", 232.04d0    ,  4.), &
        Atomic("PA", 231.04d0    ,  5.), &
        Atomic("U ", 238.03d0    ,  6.), &
        Atomic("NP", 237.048171640d0 ,  5.), &   ! No standard weight; mass of Np-237
        Atomic("PU", 244.064204401d0 ,  4.), &   ! No standard weight; mass of Pu-244
        Atomic("AM", 243.061379889d0 ,  3.), &   ! No standard weight; mass of Am-243
        Atomic("CM", 247.070352678d0 ,  3.), &   ! No standard weight; mass of Cm-247
        Atomic("BK", 247.070305889d0 ,  3.), &   ! No standard weight; mass of Bk-247
        Atomic("CF", 251.079587171d0 ,  3.), &   ! No standard weight; mass of Cf-251
        Atomic("ES", 252.082979173d0 ,  3.) &    ! No standard weight; mass of Es-252
    ]

contains

    function get_atom(symbol) result(atom)
        ! Takes an atomic symbol and returns an Atomic object

        ! Arguments
        character(enl), intent(in) :: symbol

        ! Return
        type(Atomic) :: atom

        ! Locals
        integer :: idx(1)
        integer :: i

        ! Find the index of the requested symbol
        idx = 0
        do i = 1, size(AtomicData)
            if (AtomicData(i)%name == symbol) then
                idx(1) = i
                exit
            end if
        end do

        ! Error if symbol not found
        if (idx(1) < 1) then
            call abort("Symbol "//symbol//" not found in the element list.")
        end if

        ! Select the Atomic object
        atom = AtomicData(idx(1))

    end function

    function get_atom_weight(symbol) result(weight)
        ! Takes an atomic symbol and returns its atomic weight

        ! Arguments
        character(enl), intent(in) :: symbol

        ! Return
        real(dp) :: weight

        ! Locals
        type(Atomic) :: atom

        weight = 0.0d0

        ! Skip any blank elements
        if (is_empty(symbol)) return

        ! Get the atomic object
        atom = get_atom(symbol)

        ! Return the weight
        weight = atom%atomic_weight

    end function

    function get_atom_valence(symbol) result(valence)
        ! Takes an atomic symbol and returns its valence

        ! Arguments
        character(enl), intent(in) :: symbol

        ! Return
        real(dp) :: valence

        ! Locals
        type(Atomic) :: atom

        valence = 0.0d0

        ! Skip any blank elements
        if (is_empty(symbol)) return

        ! Get the atomic object
        atom = get_atom(symbol)

        ! Return the valence
        valence = atom%valence

    end function

end module
