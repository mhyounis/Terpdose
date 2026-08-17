Subroutine FindElement (NKg, NES, NE, NENK, rg, connectivity, offset, elementset, r0, efound)
    use constants
    use math
    ! TO SPEED THIS UP, I CAN ACCEPT FAR MORE r0's AT ONCE... BUT - I would need to
    ! re-address indexing of leaves. Since I always supply the leaf that indeed contains
    ! the point. It would be tricky but absolutely doable.
    ! PLAN - basically provide r(3, n), leaf(n) (leaf of point n), and then elementset(leaf) somehow.
    ! Will likely rely on indirect addressing and flattening.
    ! Do this all later. Speed penalty is fine for now.
    
    Implicit None
    Integer,      Intent (In)  :: NKg                   ! Number of global nodes
    Integer,      Intent (In)  :: NES                   ! Number of elements to be checked
    Integer,      Intent (In)  :: NE                    ! Number of mesh elements (total)
    Integer,      Intent (In)  :: NENK                  ! Number of spatial degrees of freedom (total)
    Real (KREAL), Intent (In)  :: rg           (3, NKg) ! Global nodes
    Integer,      Intent (In)  :: connectivity (NENK)   ! Mesh connectivity
    Integer,      Intent (In)  :: offset       (NE + 1) ! Mesh offset
    Integer,      Intent (In)  :: elementset   (NES)
    Real (KREAL), Intent (In)  :: r0           (3)      ! Point to be localized
    
    Integer,      Intent (Out) :: efound
    
    Integer                    :: e, ec, k, f
    Logical                    :: crit
    Integer                    :: NKe
    Integer                    :: NFe = 4       ! Must be made dynamic when other finite element types are implemented
    Integer                    :: f2ks (3, 4)   ! This is bad practice. Temporary, until I implement for elements other than tetrahedra
    Real (KREAL)               :: TOL = 1.0e-12 ! For now not controllable
    Real (KREAL)               :: normfac
    Real (KREAL)               :: C    (3) ! Element center of mass
    Real (KREAL)               :: v1   (3) ! Vertex 1
    Real (KREAL)               :: v2   (3) ! Vertex 2
    Real (KREAL)               :: n    (3) ! Normal vector
    Real (KREAL)               :: r1   (3)
    Real (KREAL)               :: rf   (3, 6)
    Real (KREAL)               :: r    (3, 8)
    
    ! For now I've only implemented tetrahedra (NKe = 4)
    ! Approach for tetrahedra:
    !       Get the normal vectors and check, for every face, that the point r0
    !       is 'behind' the face (with normal vectors pointing outwards.)
    
    f2ks(:,1) = [2, 3, 4]
    f2ks(:,2) = [1, 3, 4]
    f2ks(:,3) = [1, 2, 4]
    f2ks(:,4) = [1, 2, 3]
    
    do ec = 1, NES
        e = elementset(ec)
        
        NKe = offset(e + 1) - offset(e)
        do k = 1, NKe
            r(1:3, k) = rg(1:3,connectivity(offset(e) + k))
        end do
        
        C = SUM(r(1:3,1:NKe), dim=2) / NKe
        
        ! Everything in this routine works in the frame of the center of element mass
        do k = 1, NKe
            r(1:3,k) = r(1:3,k) - C
        end do
        
        ! Now make the centers of the faces.
        ! Since we are, in this version of Terpdose/Lionbolt, only using tetrahedra,
        ! we just form the faces using all possible triplets out of k = 1, 2, 3, 4.
        ! Here is otherwise where we'd generalize to other finite element types.
        ! FINALLY NOTE, our convention here for face 'f' is that face f doesn't contain
        ! node k = f. However this is NOT the general convention in Lionbolt (via gmsh's conventions).
        ! Not a major note because no face information is saved here, but it's worth pointing out.
        do f = 1, NFe
            rf(1:3,f) = SUM(r(1:3,f2ks(1:3,f)), dim=2) / THREE
        end do
        
        ! Initialize the criterion as true
        crit = .TRUE.
        
        ! Visit each face and check the criterion
        do f = 1, NFe
            ! Get the normal vector
            v1 = r(1:3, f2ks(2,f)) - r(1:3, f2ks(1,f))
            v2 = r(1:3, f2ks(3,f)) - r(1:3, f2ks(1,f))
            
            n = CROSS_PRODUCT (v1, v2)
            
            ! We use the sign such that n points outwards. Magnitude doesn't matter
            normfac = SIGN(ONE, DOT_PRODUCT(n, rf(1:3,f)))
            
            n = n * normfac
            
            ! Now we must find that the dot product of the normal vector with
            ! r0 (in the frame of the face center) is negative, if r0 is in this element.
            r1 = (r0 - C) - rf(1:3, f)
            
            crit = crit .and. DOT_PRODUCT(n, r1) <= TOL
            if (.not. crit) exit ! If any face doesn't satisfy the criterion, the element does not contain r0.
        end do
        
        if (crit) then
            efound = e
            return
        end if
        
    end do
    
    ! If this point is reached, r0 was not placed in an element. Return 0 (which will later be taken to -1 in Python)
    efound = 0
    
    ! ! --- CODE FOR VALIDATING ---
    
    ! ! print *, '(', r0(1), ',', r0(2), ',', r0(3), ')'
    ! ! do ec = 1, NES
    ! !     e = elementset(ec)
    ! !     do k = 1, 4
    ! !         print *, '(', rg(1,connectivity(offset(e) + k)), ',', rg(2,connectivity(offset(e) + k)), ',', rg(3,connectivity(offset(e) + k)), ')'
    ! !     end do
    ! ! end do
    
    ! ! Now do diagnostics to see why an element was not found
    ! print *, '(', r0(1), ',', r0(2), ',', r0(3), ')'
    ! do ec = 1, NES
    !     e = elementset(ec)
        
    !     NKe = offset(e + 1) - offset(e)
    !     do k = 1, NKe
    !         r(1:3, k) = rg(1:3,connectivity(offset(e) + k))
    !     end do
    !     ! print *, '==================='
    !     do k = 1, NKe
    !         print *, '(', r(1,k), ',', r(2,k), ',', r(3,k), ')'
    !     end do
    !     ! print *, '==================='
    ! end do
    
    ! stop
    
    ! Much older code for validating. Would be a pain to write it all back out again
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! if (e == 14629) then
        !     print *, '-----------'
        !     print *, 'FACE ORIGIN : '
        !     do f = 1, NFe
        !         print *, '(', rf(1,f), ',', rf(2,f), ',', rf(3,f), ')'
        !     end do
        !     print *, '-----------'
        !     print *, 'NORMAL : '
        !     do f = 1, NFe
        !         ! Get the normal vector
        !         v1 = r(1:3, f2ks(2,f)) - r(1:3, f2ks(1,f))
        !         v2 = r(1:3, f2ks(3,f)) - r(1:3, f2ks(1,f))
              
        !         n(1) = v1(2) * v2(3) - v1(3) * v2(2)
        !         n(2) = v1(3) * v2(1) - v1(1) * v2(3)
        !         n(3) = v1(1) * v2(2) - v1(2) * v2(1)
              
        !         ! We use the sign such that n points outwards. Magnitude doesn't matter
        !         normfac = SIGN(ONE, DOT_PRODUCT(n, rf(1:3,f))) / NORM2(n)
              
        !         n = n * normfac
              
        !         print *, '(', n(1), ',', n(2), ',', n(3), ')'
        !     end do
        !     print *, '-----------'
        !     print *, 'DISPLACEMENT : '
        !     do f = 1, NFe
        !         r1 = (r0 - C - rf(1:3, f)) / NORM2(r0 - C - rf(1:3, f))
        !         print *, '(', r1(1), ',', r1(2), ',', r1(3), ')'
        !     end do
        !     print *, '-----------'
        !     print *, 'TETRAHEDRA : '
        !     do k = 1, NKe
        !         print *, '(', r(1,k), ',', r(2,k), ',', r(3,k), ')'
        !     end do
        !     print *, '-----------'
        !     print *, 'EVAL PT. : '
        !     print *, '(', r0(1) - C(1), ',', r0(2) - C(2), ',', r0(3) - C(3), ')'
        !     print *, '-----------'
        !     print *, 'CRIT : '
        !     do f = 1, NFe
        !         ! Get the normal vector
        !         v1 = r(1:3, f2ks(2,f)) - r(1:3, f2ks(1,f))
        !         v2 = r(1:3, f2ks(3,f)) - r(1:3, f2ks(1,f))
              
        !         n(1) = v1(2) * v2(3) - v1(3) * v2(2)
        !         n(2) = v1(3) * v2(1) - v1(1) * v2(3)
        !         n(3) = v1(1) * v2(2) - v1(2) * v2(1)
              
        !         ! We use the sign such that n points outwards. Magnitude doesn't matter
        !         normfac = SIGN(ONE, DOT_PRODUCT(n, rf(1:3,f))) / NORM2(n)
              
        !         n = n * normfac
              
        !         ! Now we must find that the dot product of the normal vector with
        !         ! r0 (in the frame of the face center) is negative, if r0 is in this element.
        !         r1 = (r0 - C - rf(1:3, f)) / NORM2(r0 - C - rf(1:3, f))
              
        !         print *, DOT_PRODUCT(n, r1)
        !     end do
        !     print *, '-----------'
        !     stop
        ! end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    
End Subroutine