Subroutine TallyDeposition (mNENK, NENK, mat2sd, depXS, fl, DEP)
    use constants
    Implicit None
    Integer,      Intent (In)    :: mNENK
    Integer,      Intent (In)    :: NENK
    Integer,      Intent (In)    :: mat2sd (mNENK) ! Gives spatial dofs corresponding to the desired material
    Real (KREAL), Intent (In)    :: depXS
    Real (KREAL), Intent (In)    :: fl     (NENK)  ! Fluence                                        ! (e x k)
    
    Real (KREAL), Intent (InOut) :: DEP    (NENK)  ! Deposition                                     ! (e x k)
    
    ! Should I vectorize over energy? LATER, WHEN OPTIMIZING.
    DEP(mat2sd) = DEP(mat2sd) + depXS * fl(mat2sd)
    
End Subroutine