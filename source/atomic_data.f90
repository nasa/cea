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

    type(Atomic), parameter :: AtomicData(106) = [ &
        Atomic("H ", 1.00794d0   ,  1.), &
        Atomic("IH", 1.00794d0   ,  1.), &
        Atomic("D ", 2.014102d0  ,  1.), &
        Atomic("HE", 4.002602d0  ,  0.), &
        Atomic("LI", 6.941d0     ,  1.), &
        Atomic("BE", 9.012182d0  ,  2.), &
        Atomic("B ", 10.811d0    ,  3.), &
        Atomic("C ", 12.0107d0   ,  4.), &
        Atomic("IC", 12.0107d0   ,  4.), &
        Atomic("N ", 14.0067d0   ,  0.), &
        Atomic("O ", 15.9994d0   , -2.), &
        Atomic("IO", 15.9994d0   , -2.), &
        Atomic("F ", 18.9984032d0, -1.), &
        Atomic("NE", 20.1797d0   ,  0.), &
        Atomic("NA", 22.989770d0 ,  1.), &
        Atomic("MG", 24.305d0    ,  2.), &
        Atomic("AL", 26.981538d0 ,  3.), &
        Atomic("SI", 28.0855d0   ,  4.), &
        Atomic("P ", 30.973761d0 ,  5.), &
        Atomic("S ", 32.065d0    ,  4.), &
        Atomic("CL", 35.453d0    , -1.), &
        Atomic("AR", 39.948d0    ,  0.), &
        Atomic("K ", 39.0983d0   ,  1.), &
        Atomic("CA", 40.078d0    ,  2.), &
        Atomic("SC", 44.95591d0  ,  3.), &
        Atomic("TI", 47.867d0    ,  4.), &
        Atomic("V ", 50.9415d0   ,  5.), &
        Atomic("CR", 51.9961d0   ,  3.), &
        Atomic("MN", 54.938049d0 ,  2.), &
        Atomic("FE", 55.845d0    ,  3.), &
        Atomic("CO", 58.933200d0 ,  2.), &
        Atomic("NI", 58.6934d0   ,  2.), &
        Atomic("CU", 63.546d0    ,  2.), &
        Atomic("ZN", 65.39d0     ,  2.), &
        Atomic("GA", 69.723d0    ,  3.), &
        Atomic("GE", 72.64d0     ,  4.), &
        Atomic("AS", 74.92160d0  ,  3.), &
        Atomic("SE", 78.96d0     ,  4.), &
        Atomic("BR", 79.904d0    , -1.), &
        Atomic("KR", 83.80d0     ,  0.), &
        Atomic("RB", 85.4678d0   ,  1.), &
        Atomic("SR", 87.62d0     ,  2.), &
        Atomic("Y ", 88.90585d0  ,  3.), &
        Atomic("ZR", 91.224d0    ,  4.), &
        Atomic("NB", 92.90638d0  ,  5.), &
        Atomic("MO", 95.94d0     ,  6.), &
        Atomic("TC", 97.9072d0   ,  7.), &
        Atomic("RU", 101.07d0    ,  3.), &
        Atomic("RH", 102.9055d0  ,  3.), &
        Atomic("PD", 106.42d0    ,  2.), &
        Atomic("AG", 107.8682d0  ,  1.), &
        Atomic("CD", 112.411d0   ,  2.), &
        Atomic("IN", 114.818d0   ,  3.), &
        Atomic("SN", 118.710d0   ,  4.), &
        Atomic("SB", 121.760d0   ,  3.), &
        Atomic("TE", 127.6d0     ,  4.), &
        Atomic("I ", 126.90447d0 , -1.), &
        Atomic("XE", 131.293d0   ,  0.), &
        Atomic("CS", 132.90545d0 ,  1.), &
        Atomic("BA", 137.327d0   ,  2.), &
        Atomic("LA", 138.9055d0  ,  3.), &
        Atomic("CE", 140.116d0   ,  3.), &
        Atomic("PR", 140.90765d0 ,  3.), &
        Atomic("ND", 144.9127d0  ,  3.), &
        Atomic("PM", 145.d0      ,  3.), &
        Atomic("SM", 150.36d0    ,  3.), &
        Atomic("EU", 151.964d0   ,  3.), &
        Atomic("GD", 157.25d0    ,  3.), &
        Atomic("TB", 158.92534d0 ,  3.), &
        Atomic("DY", 162.50d0    ,  3.), &
        Atomic("HO", 164.93032d0 ,  3.), &
        Atomic("ER", 167.259d0   ,  3.), &
        Atomic("TM", 168.93421d0 ,  3.), &
        Atomic("YB", 173.04d0    ,  3.), &
        Atomic("LU", 174.967d0   ,  3.), &
        Atomic("HF", 178.49d0    ,  4.), &
        Atomic("TA", 180.9479d0  ,  5.), &
        Atomic("W ", 183.84d0    ,  6.), &
        Atomic("RE", 186.207d0   ,  7.), &
        Atomic("OS", 190.23d0    ,  4.), &
        Atomic("IR", 192.217d0   ,  4.), &
        Atomic("PT", 195.078d0   ,  4.), &
        Atomic("AU", 196.96655d0 ,  3.), &
        Atomic("HG", 200.59d0    ,  2.), &
        Atomic("TL", 204.3833d0  ,  1.), &
        Atomic("PB", 207.2d0     ,  2.), &
        Atomic("BI", 208.98038d0 ,  3.), &
        Atomic("PO", 208.9824d0  ,  2.), &
        Atomic("AT", 209.9871d0  , -1.), &
        Atomic("RN", 222.0176d0  ,  0.), &
        Atomic("FR", 223.0197d0  ,  1.), &
        Atomic("RA", 226.0254d0  ,  2.), &
        Atomic("AC", 227.0278d0  ,  3.), &
        Atomic("TH", 232.0381d0  ,  4.), &
        Atomic("PA", 231.03588d0 ,  5.), &
        Atomic("U ", 238.02891d0 ,  6.), &
        Atomic("NP", 237.0482d0  ,  5.), &
        Atomic("PU", 244.0642d0  ,  4.), &
        Atomic("AM", 243.0614d0  ,  3.), &
        Atomic("CM", 247.0703d0  ,  3.), &
        Atomic("BK", 247.0703d0  ,  3.), &
        Atomic("CF", 251.0587d0  ,  3.), &
        Atomic("ES", 252.083d0   ,  3.), &
        Atomic("IH", 1.00794d0   ,  1.), &
        Atomic("IC", 12.0107d0   ,  4.), &
        Atomic("IO", 15.9994d0   , -2.) &
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
