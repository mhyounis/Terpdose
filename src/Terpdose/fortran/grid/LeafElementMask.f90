Subroutine LeafElementMask (NKg, NES, NE, NENK, rg, connectivity, offset, elementset, xb, yb, zb, keep)
    use constants
    Implicit None
    Integer,      Intent (In)  :: NKg                   ! Number of global nodes
    Integer,      Intent (In)  :: NES                   ! Number of elements to be checked
    Integer,      Intent (In)  :: NE                    ! Number of mesh elements (total)
    Integer,      Intent (In)  :: NENK                  ! Number of spatial degrees of freedom (total)
    Real (KREAL), Intent (In)  :: rg           (3, NKg) ! Global nodes
    Integer,      Intent (In)  :: connectivity (NENK)   ! Mesh connectivity
    Integer,      Intent (In)  :: offset       (NE + 1) ! Mesh offset
    Integer,      Intent (In)  :: elementset   (NES)
    Real (KREAL), Intent (In)  :: xb           (2)      ! Leaf bounds in direction x (upper/lower)
    Real (KREAL), Intent (In)  :: yb           (2)      ! Leaf bounds in direction x (upper/lower)
    Real (KREAL), Intent (In)  :: zb           (2)      ! Leaf bounds in direction x (upper/lower)
    
    Logical,      Intent (Out) :: keep         (NES)    ! Gives true for elements in the leaf, false otherwise
    
    Integer                    :: e, ec, iDir, k
    Integer                    :: NKe                   ! Number of nodes in an element
    Real (KREAL)               :: TOL = 1.0e-12
    Real (KREAL)               :: xbe          (2)      ! Element bounds
    Real (KREAL)               :: ybe          (2)      ! Element bounds
    Real (KREAL)               :: zbe          (2)      ! Element bounds
    Real (KREAL)               :: r            (3, 8)   ! Local copy of the nodes of an element. Initialized with up to 8 nodes, but in practice only 1:NKe are used ! (dir, k)
    
    ! Approach is to embed the element in a box, and then check if that box lives within the leaf bounds.
    ! 
    ! There are MANY edge cases in which an element will NOT actually be in a box, but will be kept
    ! by this routine.
    ! 
    ! However, the goal of these leaves is just to reduce the set of candidates in which a grid point
    ! may lie, so we accept this as a minor slowdown.
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! MUST CONSIDER - is the indirect addressing of elementset going to slow me down, as opposed to
    ! sending in a truncated offset and connectivity and then re-mapping the keep array in python?
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do ec = 1, NES
        e = elementset(ec)
        
        NKe = offset(e + 1) - offset(e)
        ! Identify the nodes in the element
        do k = 1, NKe
            r(1:3, k) = rg(1:3,connectivity(offset(e) + k))
        end do
        
        ! Now compare element bounds and leaf bounds
        xbe(1) = MINVAL(r(1,1:NKe))
        xbe(2) = MAXVAL(r(1,1:NKe))
        ybe(1) = MINVAL(r(2,1:NKe))
        ybe(2) = MAXVAL(r(2,1:NKe))
        zbe(1) = MINVAL(r(3,1:NKe))
        zbe(2) = MAXVAL(r(3,1:NKe))
        
        keep(ec) = .not. ( xb(2) < xbe(1) - tol .or. xb(1) > xbe(2) + tol .or. &
                           yb(2) < ybe(1) - tol .or. yb(1) > ybe(2) + tol .or. &
                           zb(2) < zbe(1) - tol .or. zb(1) > zbe(2) + tol        )
        
    end do
    
End Subroutine