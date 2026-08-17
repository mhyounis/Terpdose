Module constants
    Implicit None

!  ========================================================
!    Terpdose works in the unit system:
!      Energy : MeV
!      Length : cm
!      Mass   : g
!      Time   : s
!      Charge : e
!  ========================================================

Integer, Parameter :: KREAL = KIND(0.0d0)

! Numbers
Real (KREAL), Parameter :: ZERO    = 0.0e+0_KREAL
Real (KREAL), Parameter :: ONE     = 1.0e+0_KREAL
Real (KREAL), Parameter :: TWO     = 2.0e+0_KREAL
Real (KREAL), Parameter :: THREE   = 3.0e+0_KREAL
Real (KREAL), Parameter :: FOUR    = 4.0e+0_KREAL
Real (KREAL), Parameter :: FIVE    = 5.0e+0_KREAL
Real (KREAL), Parameter :: SIX     = 6.0e+0_KREAL
Real (KREAL), Parameter :: SEVEN   = 7.0e+0_KREAL
Real (KREAL), Parameter :: EIGHT   = 8.0e+0_KREAL
Real (KREAL), Parameter :: NINE    = 9.0e+0_KREAL
Real (KREAL), Parameter :: TEN     = 1.0e+1_KREAL
Real (KREAL), Parameter :: HUNDRED = 1.0e+2_KREAL
Real (KREAL), Parameter :: HALF    = 0.5e+0_KREAL
Real (KREAL), Parameter :: FOURTH  = 2.5e-1_KREAL
Real (KREAL), Parameter :: THIRD   = ONE / THREE
Real (KREAL), Parameter :: SIXTH   = ONE / SIX
Real (KREAL), Parameter :: EIGTH   = ONE / EIGHT

! Math constants
Real (KREAL), Parameter :: PI      = 3.14159265358979323846_KREAL
Real (KREAL), Parameter :: TWOPI   = TWO * PI
Real (KREAL), Parameter :: FOURPI  = TWO * TWOPI
Real (KREAL), Parameter :: HALFPI  = HALF * PI
Real (KREAL), Parameter :: LNTWO   = LOG(TWO)
Real (KREAL), Parameter :: SQRTTWO = SQRT(TWO)

! Unit Conversions
Real (KREAL), Parameter :: EV2ELM  = ONE / 510998.95069_KREAL     ! Electron-Volts         to    Electron mass energies
Real (KREAL), Parameter :: EV2MEV  = 1.0e-6_KREAL                 ! Electron-Volts         to    Mega electron-Volts
Real (KREAL), Parameter :: MEV2ELM = EV2ELM * 1.0e+6_KREAL        ! Mega electron-Volts    to    Electron mass energies
Real (KREAL), Parameter :: J2ELM   = ONE / 8.1871057880e-14_KREAL ! Joules                 to    Electron mass energies
Real (KREAL), Parameter :: CM2ELR  = ONE / 2.817939e-13_KREAL     ! Centimeters            to    Electron radii
Real (KREAL), Parameter :: G2ELM   = ONE / 9.1093837139e-28_KREAL ! Grams                  to    Electron masses
Real (KREAL), Parameter :: B2SQCM  = 1.0e-24_KREAL                ! Barns                  to    Square centimeters

! Physical constants
Real (KREAL), Parameter :: FSA  = 0.0072973525693e+0_KREAL ! Fine structure constant
Real (KREAL), Parameter :: ELM  = 0.51099895069e+0_KREAL   ! Electron mass energy
Real (KREAL), Parameter :: HBAR = 6.582119569e-22_KREAL    ! Planck's reduced constant
Real (KREAL), Parameter :: ELR  = 2.8179403205e-13_KREAL   ! Classical electron radius
Real (KREAL), Parameter :: AVO  = 6.02214076e+23_KREAL     ! Avogadro's constant

End Module constants