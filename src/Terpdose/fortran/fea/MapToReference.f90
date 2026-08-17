Subroutine MapToReference (n, NKg, NE, NENK, rg, connectivity, offset, grid, elements, mask, refgrid)
    use constants
    use math
    
    !  ==========================================================================================
    !    This routine takes a grid, including information about which element a grid point lies
    !    within, and determines the grid point in the frame of the reference element.
    !    
    !    The idea is that interpolation is performed using reference shape elements, and these
    !    accept the reference grid points.
    !  ==========================================================================================
    
    Implicit None
    Integer,      Intent (In)  :: n
    Integer,      Intent (In)  :: NKg
    Integer,      Intent (In)  :: NE
    Integer,      Intent (In)  :: NENK
    Real (KREAL), Intent (In)  :: rg           (3, NKg)
    Integer,      Intent (In)  :: connectivity (NENK)   ! Mesh connectivity
    Integer,      Intent (In)  :: offset       (NE + 1) ! Mesh offset
    Real (KREAL), Intent (In)  :: grid         (3, n)
    Integer,      Intent (In)  :: elements     (n)
    Logical,      Intent (In)  :: mask         (n)
    
    Real (KREAL), Intent (Out) :: refgrid      (3, n)   ! Shape functions evaluated at grid points ! MUST ZERO-PAD IN PYTHON...
    
    Integer                    :: e, i, iDir, k
    Integer                    :: NKe
    Real (KREAL)               :: invd
    Real (KREAL)               :: r1 (3)
    Real (KREAL)               :: R0 (3)
    Real (KREAL)               :: r  (3, 8)
    Real (KREAL)               :: ev (3, 3)
    Real (KREAL)               :: cv (3, 3)
    
    do i = 1, n
        e = elements(i)
        
        if (mask(i)) cycle
        
        NKe = offset(e + 1) - offset(e)
        
        do k = 1, NKe
            r(1:3, k) = rg(1:3, connectivity(offset(e) + k))
        end do
        
        ! BEGIN TETRAHEDRON-SPECIFIC CODE
        
        r1 = r(1:3, 1)
        
        do k = 1, 3
            ev(1:3, k) = r(1:3, k + 1) - r1
        end do
        
        cv(1:3, 1) = CROSS_PRODUCT (ev(1:3, 2), ev(1:3, 3))
        cv(1:3, 2) = CROSS_PRODUCT (ev(1:3, 3), ev(1:3, 1))
        cv(1:3, 3) = CROSS_PRODUCT (ev(1:3, 1), ev(1:3, 2))
        
        invd = ONE / DOT_PRODUCT(ev(1:3, 1), cv(1:3, 1))
        
        R0 = grid(1:3, i) - r1
        
        do iDir = 1, 3
            refgrid(iDir, i) = DOT_PRODUCT(R0, cv(1:3, iDir)) * invd
        end do
        ! END TETRAHEDRON-SPECIFIC CODE
        
    end do
    
End Subroutine