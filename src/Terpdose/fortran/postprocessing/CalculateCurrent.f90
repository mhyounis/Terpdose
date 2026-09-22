Subroutine CalculateCurrent (d, NENK, NI, w, khat, psi, J)
    use constants
    Implicit None
    Integer,      Intent (In)  :: d
    Integer,      Intent (In)  :: NENK ! Number of spatial dofs
    Integer,      Intent (In)  :: NI   ! Number of angles
    Real (KREAL), Intent (In)  :: w    (NI)
    Real (KREAL), Intent (In)  :: khat (NI, d)
    Real (KREAL), Intent (In)  :: psi  (NENK, NI)
    
    Real (KREAL), Intent (Out) :: J    (NENK, d)
    
    Integer                    :: iDir
    Real (KREAL)               :: tmp  (NI)
    
    do iDir = 1, 3
        tmp = w * khat(1:NI, iDir)
        J(1:NENK, iDir) = MATMUL(psi, tmp)
    end do
    
End Subroutine