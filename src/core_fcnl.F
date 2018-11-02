!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Calculation of the non-local pseudopotential type contribution to the core Hamiltonian
!>        from the frozen core states
!>         <a|FC(non-local)|b> = <a|core(i)*h(i)*<core(i)|b>
!> \par History
!>      - copied from core_ppnl [jhu, 2018-11-01]
! **************************************************************************************************
MODULE core_fcnl
   USE ai_contraction,                  ONLY: block_add,&
                                              contraction
   USE ai_overlap,                      ONLY: overlap_ab
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE dbcsr_api,                       ONLY: dbcsr_add,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type
   USE kinds,                           ONLY: dp
   USE orbital_pointers,                ONLY: init_orbital_pointers,&
                                              nco,&
                                              ncoset
   USE particle_types,                  ONLY: particle_type
   USE qs_force_types,                  ONLY: qs_force_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE sap_kind_types,                  ONLY: alist_type,&
                                              clist_type,&
                                              get_alist,&
                                              release_sap_int,&
                                              sap_int_type,&
                                              sap_sort
   USE virial_methods,                  ONLY: virial_pair_force
   USE virial_types,                    ONLY: virial_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'core_fcnl'

   PUBLIC :: build_core_fcnl

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param matrix_h ...
!> \param matrix_p ...
!> \param force ...
!> \param virial ...
!> \param calculate_forces ...
!> \param use_virial ...
!> \param nder ...
!> \param qs_kind_set ...
!> \param atomic_kind_set ...
!> \param particle_set ...
!> \param sab_orb ...
!> \param sac_fc ...
!> \param eps_ppnl ...
!> \param nimages ...
!> \param cell_to_index ...
!> \param basis_type ...
! **************************************************************************************************
   SUBROUTINE build_core_fcnl(matrix_h, matrix_p, force, virial, calculate_forces, use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_fc, eps_ppnl, &
                              nimages, cell_to_index, basis_type)

      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(virial_type), POINTER                         :: virial
      LOGICAL, INTENT(IN)                                :: calculate_forces
      LOGICAL                                            :: use_virial
      INTEGER                                            :: nder
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_fc
      REAL(KIND=dp), INTENT(IN)                          :: eps_ppnl
      INTEGER, INTENT(IN)                                :: nimages
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      CHARACTER(LEN=*), INTENT(IN)                       :: basis_type

      CHARACTER(LEN=*), PARAMETER :: routineN = 'build_core_fcnl', &
         routineP = moduleN//':'//routineN

      INTEGER :: atom_a, atom_b, atom_c, first_col, handle, i, iab, iac, iatom, ibc, icol, ikind, &
         ilist, img, inode, irow, iset, j, jatom, jkind, jneighbor, jset, kac, katom, kbc, kkind, &
         l, lc_max, lc_min, ldai, ldsab, maxcfc, maxco, maxder, maxl, maxlfc, maxlgto, maxsfc, &
         maxsgf, mepos, n1, n2, na, natom, nb, nkind, nlist, nneighbor, nnl, nnode, np, nprjc, &
         nseta, nsetb, nsgfa, nsgfb, nthread, prjc, sgfa, sgfb
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: atom_of_kind
      INTEGER, DIMENSION(3)                              :: cell_b, cell_c
      INTEGER, DIMENSION(:), POINTER                     :: la_max, la_min, lb_max, lb_min, npgfa, &
                                                            npgfb, nsgf_seta, nsgf_setb
      INTEGER, DIMENSION(:, :), POINTER                  :: first_sgfa, first_sgfb
      LOGICAL                                            :: dogth, dokp, fc_present, found
      LOGICAL, DIMENSION(0:9)                            :: is_nonlocal
      REAL(KIND=dp)                                      :: dac, f0
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: escale, owork
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: oint
      REAL(KIND=dp), DIMENSION(3)                        :: fa, fb, rab, rac, rbc
      REAL(KIND=dp), DIMENSION(:), POINTER               :: a_nl, set_radius_a, set_radius_b
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: h_block, h_nl, hshift, p_block, rpgfa, &
                                                            rpgfb, scon_a, scon_b, zeta, zetb
      REAL(KIND=dp), DIMENSION(:, :, :), POINTER         :: achint, acint, bchint, bcint, c_nl
      TYPE(alist_type), POINTER                          :: alist_ac, alist_bc
      TYPE(clist_type), POINTER                          :: clist
      TYPE(gto_basis_set_p_type), DIMENSION(:), POINTER  :: basis_set, fc_basis
      TYPE(gto_basis_set_type), POINTER                  :: fc_basis_set, orb_basis_set
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(sap_int_type), DIMENSION(:), POINTER          :: sap_int

      IF (calculate_forces) THEN
         CALL timeset(routineN//"_forces", handle)
      ELSE
         CALL timeset(routineN, handle)
      ENDIF

      fc_present = ASSOCIATED(sac_fc)

      IF (fc_present) THEN

         nkind = SIZE(atomic_kind_set)
         natom = SIZE(particle_set)

         dokp = (nimages > 1)

         ALLOCATE (atom_of_kind(natom))
         CALL get_atomic_kind_set(atomic_kind_set, atom_of_kind=atom_of_kind)

         IF (calculate_forces) THEN
            IF (SIZE(matrix_p, 1) == 2) THEN
               DO img = 1, nimages
                  CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                                 alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
                  CALL dbcsr_add(matrix_p(2, img)%matrix, matrix_p(1, img)%matrix, &
                                 alpha_scalar=-2.0_dp, beta_scalar=1.0_dp)
               END DO
            END IF
         END IF

         maxder = ncoset(nder)

         CALL get_qs_kind_set(qs_kind_set, maxco=maxco, maxlgto=maxlgto, maxsgf=maxsgf, basis_type=basis_type)
         CALL get_qs_kind_set(qs_kind_set, maxco=maxcfc, maxlgto=maxlfc, maxsgf=maxsfc, basis_type="FC")

         maxl = MAX(maxlgto, maxlfc)
         CALL init_orbital_pointers(maxl+nder+1)

         ldsab = MAX(maxco, maxcfc, maxsgf, maxsfc)
         ldai = ncoset(maxl+nder+1)

         !sap_int needs to be shared as multiple threads need to access this
         ALLOCATE (sap_int(nkind*nkind))
         DO i = 1, nkind*nkind
            NULLIFY (sap_int(i)%alist, sap_int(i)%asort, sap_int(i)%aindex)
            sap_int(i)%nalist = 0
         END DO

         !set up direct access to the basis sets
         ALLOCATE (basis_set(nkind), fc_basis(nkind))
         ALLOCATE (escale(maxsfc, nkind))
         escale = 0.0_dp
         DO ikind = 1, nkind
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set, basis_type=basis_type)
            IF (ASSOCIATED(orb_basis_set)) THEN
               basis_set(ikind)%gto_basis_set => orb_basis_set
            ELSE
               NULLIFY (basis_set(ikind)%gto_basis_set)
            END IF
            CALL get_qs_kind(qs_kind_set(ikind), basis_set=fc_basis_set, basis_type="FC")
            IF (ASSOCIATED(fc_basis_set)) THEN
               fc_basis(ikind)%gto_basis_set => fc_basis_set
               hshift => qs_kind_set(ikind)%hshift
               DO i = 1, fc_basis_set%nset
                  DO j = 1, fc_basis_set%nshell(i)
                     n1 = fc_basis_set%first_sgf(j, i)
                     n2 = fc_basis_set%last_sgf(j, i)
                     escale(n1:n2, ikind) = hshift(j, i)
                  END DO
               END DO
            ELSE
               NULLIFY (fc_basis(ikind)%gto_basis_set)
            END IF
         END DO

         nthread = 1
!$       nthread = omp_get_max_threads()

         !calculate the overlap integrals <a|p>
         CALL neighbor_list_iterator_create(nl_iterator, sac_fc, nthread=nthread)
!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, basis_set, fc_basis, maxder, ncoset, &
!$OMP          sap_int, nkind, ldsab, ldai, nder, nco, escale ) &
!$OMP PRIVATE (mepos, oint, owork, ikind, kkind, iatom, katom, nlist, ilist, nneighbor, jneighbor, &
!$OMP          cell_c, rac, iac, first_sgfa, la_max, la_min, npgfa, nseta, nsgfa, nsgf_seta, &
!$OMP          scon_a, zeta, first_sgfb, lb_max, lb_min, npgfb, nsetb, nsgfb, nsgf_setb, &
!$OMP          scon_b, zetb, n1, n2, &
!$OMP          clist, iset, jset, sgfa, sgfb, prjc, nprjc, &
!$OMP          rpgfa, first_col, vprj_ppnl, i, j, l, dogth, &
!$OMP          set_radius_a, rprjc, dac, lc_max, lc_min, &
!$OMP          nnl, is_nonlocal, a_nl, h_nl, c_nl)
         mepos = 0
!$       mepos = omp_get_thread_num()

         ALLOCATE (oint(ldsab, ldsab, maxder), owork(ldsab, ldsab))

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=kkind, iatom=iatom, &
                                   jatom=katom, nlist=nlist, ilist=ilist, nnode=nneighbor, &
                                   inode=jneighbor, cell=cell_c, r=rac)
            iac = ikind+nkind*(kkind-1)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            IF (.NOT. ASSOCIATED(fc_basis(kkind)%gto_basis_set)) CYCLE
            ! get definition of basis set
            first_sgfa => basis_set(ikind)%gto_basis_set%first_sgf
            la_max => basis_set(ikind)%gto_basis_set%lmax
            la_min => basis_set(ikind)%gto_basis_set%lmin
            npgfa => basis_set(ikind)%gto_basis_set%npgf
            nseta = basis_set(ikind)%gto_basis_set%nset
            nsgfa = basis_set(ikind)%gto_basis_set%nsgf
            nsgf_seta => basis_set(ikind)%gto_basis_set%nsgf_set
            rpgfa => basis_set(ikind)%gto_basis_set%pgf_radius
            set_radius_a => basis_set(ikind)%gto_basis_set%set_radius
            scon_a => basis_set(ikind)%gto_basis_set%scon
            zeta => basis_set(ikind)%gto_basis_set%zet
            ! get definition of FC orbitals
            first_sgfb => fc_basis(kkind)%gto_basis_set%first_sgf
            lb_max => fc_basis(kkind)%gto_basis_set%lmax
            lb_min => fc_basis(kkind)%gto_basis_set%lmin
            npgfb => fc_basis(kkind)%gto_basis_set%npgf
            nsetb = fc_basis(kkind)%gto_basis_set%nset
            nsgfb = fc_basis(kkind)%gto_basis_set%nsgf
            nsgf_setb => fc_basis(kkind)%gto_basis_set%nsgf_set
            rpgfb => fc_basis(kkind)%gto_basis_set%pgf_radius
            set_radius_b => fc_basis(kkind)%gto_basis_set%set_radius
            scon_b => fc_basis(kkind)%gto_basis_set%scon
            zetb => fc_basis(kkind)%gto_basis_set%zet
!$OMP CRITICAL(sap_int_critical)
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) THEN
               sap_int(iac)%a_kind = ikind
               sap_int(iac)%p_kind = kkind
               sap_int(iac)%nalist = nlist
               ALLOCATE (sap_int(iac)%alist(nlist))
               DO i = 1, nlist
                  NULLIFY (sap_int(iac)%alist(i)%clist)
                  sap_int(iac)%alist(i)%aatom = 0
                  sap_int(iac)%alist(i)%nclist = 0
               END DO
            END IF
            IF (.NOT. ASSOCIATED(sap_int(iac)%alist(ilist)%clist)) THEN
               sap_int(iac)%alist(ilist)%aatom = iatom
               sap_int(iac)%alist(ilist)%nclist = nneighbor
               ALLOCATE (sap_int(iac)%alist(ilist)%clist(nneighbor))
               DO i = 1, nneighbor
                  sap_int(iac)%alist(ilist)%clist(i)%catom = 0
               END DO
            END IF
!$OMP END CRITICAL(sap_int_critical)
            dac = SQRT(SUM(rac*rac))
            clist => sap_int(iac)%alist(ilist)%clist(jneighbor)
            clist%catom = katom
            clist%cell = cell_c
            clist%rac = rac
            ALLOCATE (clist%acint(nsgfa, nsgfb, maxder), &
                      clist%achint(nsgfa, nsgfb, maxder))
            clist%acint = 0._dp
            clist%achint = 0._dp
            clist%nsgf_cnt = 0
            NULLIFY (clist%sgf_list)
            DO iset = 1, nseta
               n1 = npgfa(iset)*(ncoset(la_max(iset))-ncoset(la_min(iset)-1))
               sgfa = first_sgfa(1, iset)
               DO jset = 1, nsetb
                  IF (set_radius_a(iset)+set_radius_b(jset) < dac) CYCLE
                  n2 = npgfb(jset)*(ncoset(lb_max(jset))-ncoset(lb_min(jset)-1))
                  sgfb = first_sgfb(1, jset)
                  ! calculate integrals and derivatives
                  SELECT CASE (nder)
                  CASE (0)
                     CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                     lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                                     rab, sab=oint(:, :, 1))
                  CASE (1)
                     CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                     lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                                     rab, sab=oint(:, :, 1), dab=oint(:, :, 2:4))
                  CASE (2)
                     CALL overlap_ab(la_max(iset), la_min(iset), npgfa(iset), rpgfa(:, iset), zeta(:, iset), &
                                     lb_max(jset), lb_min(jset), npgfb(jset), rpgfb(:, jset), zetb(:, jset), &
                                     rab, sab=oint(:, :, 1), dab=oint(:, :, 2:4), ddab=oint(:, :, 5:10))
                  CASE DEFAULT
                     CPABORT("")
                  END SELECT
                  !
                  DO i = 1, maxder
                     CALL contraction(oint(:, :, i), owork, ca=scon_a(:, sgfa:), na=n1, ma=nsgf_seta(iset), &
                                      cb=scon_b(:, sgfb:), nb=n2, mb=nsgf_setb(jset))
!$OMP CRITICAL(blockadd)
                     CALL block_add("IN", owork, nsgf_seta(iset), nsgf_setb(jset), clist%acint(:, :, i), sgfa, sgfb)
!$OMP END CRITICAL(blockadd)
                  END DO
               END DO
            END DO
            ! scale by orbital shift
            DO i = 1, maxder
               DO j = 1, nsgfb
                  clist%achint(:, j, i) = clist%acint(:, j, i)*escale(j, kkind)
               END DO
            END DO
            !
            clist%maxac = MAXVAL(ABS(clist%acint(:, :, 1)))
            clist%maxach = MAXVAL(ABS(clist%achint(:, :, 1)))
         END DO

         DEALLOCATE (oint, owork)
!$OMP END PARALLEL
         CALL neighbor_list_iterator_release(nl_iterator)

         ! *** Set up a sorting index
         CALL sap_sort(sap_int)
         ! *** All integrals needed have been calculated and stored in sap_int
         ! *** We now calculate the Hamiltonian matrix elements
         CALL neighbor_list_iterator_create(nl_iterator, sab_orb, nthread=nthread)

!$OMP PARALLEL &
!$OMP DEFAULT (NONE) &
!$OMP SHARED  (nl_iterator, dokp, basis_set, atom_of_kind, matrix_h, cell_to_index,&
!$OMP          matrix_p, sap_int, nkind, eps_ppnl, force, virial, use_virial, calculate_forces) &
!$OMP PRIVATE (mepos, ikind, jkind, iatom, jatom, nlist, ilist, nnode, inode, cell_b, rab, &
!$OMP          iab, atom_a, atom_b, f0, irow, icol, h_block, &
!$OMP          found,p_block, iac, ibc, alist_ac, alist_bc, acint, bcint, &
!$OMP          achint, bchint, na, np, nb, katom, atom_c, j, fa, fb, rbc, rac, &
!$OMP          kkind, kac, kbc, i, img)

         mepos = 0
!$       mepos = omp_get_thread_num()

         DO WHILE (neighbor_list_iterate(nl_iterator, mepos=mepos) == 0)
            CALL get_iterator_info(nl_iterator, mepos=mepos, ikind=ikind, jkind=jkind, iatom=iatom, &
                                   jatom=jatom, nlist=nlist, ilist=ilist, nnode=nnode, inode=inode, cell=cell_b, r=rab)
            IF (.NOT. ASSOCIATED(basis_set(ikind)%gto_basis_set)) CYCLE
            IF (.NOT. ASSOCIATED(basis_set(jkind)%gto_basis_set)) CYCLE
            iab = ikind+nkind*(jkind-1)
            atom_a = atom_of_kind(iatom)
            atom_b = atom_of_kind(jatom)

            ! *** Use the symmetry of the first derivatives ***
            IF (iatom == jatom) THEN
               f0 = 1.0_dp
            ELSE
               f0 = 2.0_dp
            END IF

            IF (dokp) THEN
               img = cell_to_index(cell_b(1), cell_b(2), cell_b(3))
            ELSE
               img = 1
            END IF

            ! *** Create matrix blocks for a new matrix block column ***
            IF (iatom <= jatom) THEN
               irow = iatom
               icol = jatom
            ELSE
               irow = jatom
               icol = iatom
            END IF
            NULLIFY (h_block)
            CALL dbcsr_get_block_p(matrix_h(1, img)%matrix, irow, icol, h_block, found)
            IF (calculate_forces) THEN
               NULLIFY (p_block)
               CALL dbcsr_get_block_p(matrix_p(1, img)%matrix, irow, icol, p_block, found)
            END IF

            ! loop over all kinds for projector atom
            IF (ASSOCIATED(h_block)) THEN
               DO kkind = 1, nkind
                  iac = ikind+nkind*(kkind-1)
                  ibc = jkind+nkind*(kkind-1)
                  IF (.NOT. ASSOCIATED(sap_int(iac)%alist)) CYCLE
                  IF (.NOT. ASSOCIATED(sap_int(ibc)%alist)) CYCLE
                  CALL get_alist(sap_int(iac), alist_ac, iatom)
                  CALL get_alist(sap_int(ibc), alist_bc, jatom)
                  IF (.NOT. ASSOCIATED(alist_ac)) CYCLE
                  IF (.NOT. ASSOCIATED(alist_bc)) CYCLE
                  DO kac = 1, alist_ac%nclist
                     DO kbc = 1, alist_bc%nclist
                        IF (alist_ac%clist(kac)%catom /= alist_bc%clist(kbc)%catom) CYCLE
                        IF (ALL(cell_b+alist_bc%clist(kbc)%cell-alist_ac%clist(kac)%cell == 0)) THEN
                           IF (alist_ac%clist(kac)%maxac*alist_bc%clist(kbc)%maxach < eps_ppnl) CYCLE
                           acint => alist_ac%clist(kac)%acint
                           bcint => alist_bc%clist(kbc)%acint
                           achint => alist_ac%clist(kac)%achint
                           bchint => alist_bc%clist(kbc)%achint
                           na = SIZE(acint, 1)
                           np = SIZE(acint, 2)
                           nb = SIZE(bcint, 1)
!$OMP CRITICAL(h_block_critical)
                           IF (iatom <= jatom) THEN
                              CALL dgemm("N", "T", na, nb, np, 1._dp, achint(1, 1, 1), na, &
                                         bcint(1, 1, 1), nb, 1.0_dp, h_block, SIZE(h_block, 1))
                           ELSE
                              CALL dgemm("N", "T", nb, na, np, 1.0_dp, bchint(1, 1, 1), nb, &
                                         acint(1, 1, 1), na, 1.0_dp, h_block, SIZE(h_block, 1))
                           END IF
!$OMP END CRITICAL(h_block_critical)
                           IF (calculate_forces) THEN
                              IF (ASSOCIATED(p_block)) THEN
                                 katom = alist_ac%clist(kac)%catom
                                 atom_c = atom_of_kind(katom)
                                 DO i = 1, 3
                                    j = i+1
                                    IF (iatom <= jatom) THEN
                                       fa(i) = SUM(p_block(1:na, 1:nb)* &
                                                   MATMUL(acint(1:na, 1:np, j), TRANSPOSE(bchint(1:nb, 1:np, 1))))
                                       fb(i) = SUM(p_block(1:na, 1:nb)* &
                                                   MATMUL(achint(1:na, 1:np, 1), TRANSPOSE(bcint(1:nb, 1:np, j))))
                                    ELSE
                                       fa(i) = SUM(p_block(1:nb, 1:na)* &
                                                   MATMUL(bchint(1:nb, 1:np, 1), TRANSPOSE(acint(1:na, 1:np, j))))
                                       fb(i) = SUM(p_block(1:nb, 1:na)* &
                                                   MATMUL(bcint(1:nb, 1:np, j), TRANSPOSE(achint(1:na, 1:np, 1))))
                                    END IF
!$OMP CRITICAL(force_critical)
                                    force(ikind)%gth_ppnl(i, atom_a) = force(ikind)%gth_ppnl(i, atom_a)+f0*fa(i)
                                    force(kkind)%gth_ppnl(i, atom_c) = force(kkind)%gth_ppnl(i, atom_c)-f0*fa(i)
                                    force(jkind)%gth_ppnl(i, atom_b) = force(jkind)%gth_ppnl(i, atom_b)+f0*fb(i)
                                    force(kkind)%gth_ppnl(i, atom_c) = force(kkind)%gth_ppnl(i, atom_c)-f0*fb(i)
!$OMP END CRITICAL(force_critical)
                                 END DO

                                 IF (use_virial) THEN
                                    rac = alist_ac%clist(kac)%rac
                                    rbc = alist_bc%clist(kbc)%rac
!$OMP CRITICAL(virial_critical)
                                    CALL virial_pair_force(virial%pv_virial, f0, fa, rac)
                                    CALL virial_pair_force(virial%pv_virial, f0, fb, rbc)
!$OMP END CRITICAL(virial_critical)
                                 END IF
                              ENDIF
                           END IF
                           EXIT ! We have found a match and there can be only one single match
                        END IF
                     END DO
                  END DO
               END DO
            ENDIF
         END DO
!$OMP END PARALLEL
         CALL neighbor_list_iterator_release(nl_iterator)

         CALL release_sap_int(sap_int)

         DEALLOCATE (atom_of_kind)
         DEALLOCATE (basis_set, fc_basis)

         IF (calculate_forces) THEN
            ! *** If LSD, then recover alpha density and beta density     ***
            ! *** from the total density (1) and the spin density (2)     ***
            IF (SIZE(matrix_p, 1) == 2) THEN
               DO img = 1, nimages
                  CALL dbcsr_add(matrix_p(1, img)%matrix, matrix_p(2, img)%matrix, &
                                 alpha_scalar=0.5_dp, beta_scalar=0.5_dp)
                  CALL dbcsr_add(matrix_p(2, img)%matrix, matrix_p(1, img)%matrix, &
                                 alpha_scalar=-1.0_dp, beta_scalar=1.0_dp)
               END DO
            END IF
         END IF

      END IF !fc_present

      CALL timestop(handle)

   END SUBROUTINE build_core_fcnl

! **************************************************************************************************

END MODULE core_fcnl
