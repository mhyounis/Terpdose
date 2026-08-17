Module math
    use constants
    Implicit None
    
Contains

Function CROSS_PRODUCT (v, w) Result (u)
    Implicit None
    Real (KREAL), Intent (In) :: v (3)
    Real (KREAl), Intent (In) :: w (3)
    
    Real (KREAL)              :: u (3)
    
    u(1) = v(2) * w(3) - v(3) * w(2)
    u(2) = v(3) * w(1) - v(1) * w(3)
    u(3) = v(1) * w(2) - v(2) * w(1)
    
End Function

End Module