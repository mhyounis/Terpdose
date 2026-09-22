Subroutine CalculateFluence (NENK, NI, w, psi, fl)
    use constants
    Implicit None
    Integer,      Intent (In)  :: NENK ! Number of spatial dofs
    Integer,      Intent (In)  :: NI   ! Number of angles
    Real (KREAL), Intent (In)  :: w   (NI)
    Real (KREAL), Intent (In)  :: psi (NENK, NI)
    
    Real (KREAL), Intent (Out) :: fl  (NENK)
    
    Integer                    :: i
    Integer                    :: start
    
    ! This and TallyDeposition are actually rather simple operations.
    ! Is it possible numpy/scipy or some other package will have quicker and cleaner
    ! ways of doing them?
    
    fl = MATMUL(psi, w)
    
End Subroutine