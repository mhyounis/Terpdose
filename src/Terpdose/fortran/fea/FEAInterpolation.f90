Subroutine FEAInterpolation (n, NE, NENK, offset, refgrid, elements, mask, arr, arrm)
    use constants
    use ShapeFunctions
    Implicit None
    Integer,      Intent (In)  :: n
    Integer,      Intent (In)  :: NE
    Integer,      Intent (In)  :: NENK
    Integer,      Intent (In)  :: offset       (NE + 1) ! Mesh offset
    Real (KREAL), Intent (In)  :: refgrid      (3, n)
    Integer,      Intent (In)  :: elements     (n)
    Logical,      Intent (In)  :: mask         (n)
    Real (KREAL), Intent (In)  :: arr          (NENK)   ! Array to be interpolated
    
    Real (KREAL), Intent (Out) :: arrm         (n)      ! Array after interpolation
    
    Integer                    :: e, i, k
    Integer                    :: NKe
    Real (KREAL)               :: u
    Real (KREAL)               :: val
    
    arrm = ZERO
    
    do i = 1, n
        if (mask(i)) cycle ! arrm(i) will be culled in Python anyway
        
        e = elements(i)
        
        NKe = offset(e + 1) - offset(e)
        
        do k = 1, NKe
            u   = TetrahedralSF (k, refgrid(1:3,i))
            val = arr(offset(e) + k)
            
            arrm(i) = arrm(i) + u * val
        end do
        
    end do
    
End Subroutine