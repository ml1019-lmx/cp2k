!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Frozen Core
! **************************************************************************************************
MODULE qs_fc_methods
   USE atom_kind_orbitals,              ONLY: calculate_atomic_energies
   USE atom_types,                      ONLY: atom_basis_type,&
                                              init_atom_basis_cgto,&
                                              release_atom_basis
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE kinds,                           ONLY: dp
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type,&
                                              set_qs_kind
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters (only in this module)
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fc_methods'

   PUBLIC :: set_fc_kind

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_kind_set ...
! **************************************************************************************************
   SUBROUTINE set_fc_kind(qs_kind_set)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CHARACTER(len=*), PARAMETER :: routineN = 'set_fc_kind', routineP = moduleN//':'//routineN

      INTEGER                                            :: i, ii, ikind, iounit, l, nkind, nn, ns, &
                                                            nset, nsmax
      LOGICAL                                            :: fcpresent, frozen_core
      REAL(KIND=dp)                                      :: ecore, efc, ekin, exchange, fc, ortho
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: pmat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hshift
      TYPE(atom_basis_type), POINTER                     :: basis
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: fc_basis
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL get_qs_kind_set(qs_kind_set, frozen_core_present=fcpresent)

      IF (fcpresent) THEN

         logger => cp_get_default_logger()
         iounit = cp_logger_get_default_io_unit(logger)

         nkind = SIZE(qs_kind_set)
         DO ikind = 1, nkind
            qs_kind => qs_kind_set(ikind)
            CALL get_qs_kind(qs_kind, frozen_core=frozen_core)
            IF (frozen_core) THEN
               CALL get_qs_kind(qs_kind, basis_set=fc_basis, basis_type="FC")
               IF (.NOT. ASSOCIATED(fc_basis)) THEN
                  CALL cp_abort(__LOCATION__, qs_kind%name//": Atoms with frozen core need a FC basis defined")
               END IF
               CALL get_qs_kind(qs_kind, grid_atom=grid_atom)
               CPASSERT(ASSOCIATED(grid_atom))
               ! initialize atomic basis
               ALLOCATE (basis)
               NULLIFY (basis%grid)
               CALL init_atom_basis_cgto(basis, fc_basis%nset, fc_basis%lmin, fc_basis%lmax, fc_basis%npgf, &
                                         fc_basis%nshell, fc_basis%l, fc_basis%zet, fc_basis%gcc, &
                                         grid_atom%quadrature, grid_atom%nr)
               ! calculate integrals
               nn = MAXVAL(basis%nbas)
               ALLOCATE (pmat(nn, nn, 0:3))
               pmat = 0.0_dp
               DO l = 0, 3
                  DO ii = 1, basis%nbas(l)
                     pmat(ii, ii, l) = REAL(2*(2*l+1), KIND=dp)
                  END DO
               END DO
               ALLOCATE (hshift(nn, 0:3))
               ! needs some programming for exchange
               exchange = 0.0_dp
               CALL calculate_atomic_energies(qs_kind, basis, pmat=pmat, ekin=ekin, ecore=ecore, &
                                              hshift=hshift, ortho=ortho)
               efc = ekin+ecore+exchange
               CALL set_qs_kind(qs_kind, efc=efc)
               !
               nset = fc_basis%nset
               nsmax = MAXVAL(fc_basis%nshell)
               ! orbital projector shifts (nshell,nset)
               ALLOCATE (qs_kind%hshift(nsmax, nset))
               qs_kind%hshift = 0.0_dp
               DO l = 0, 3
                  ns = 0
                  DO i = 1, fc_basis%nset
                     IF (l >= fc_basis%lmin(i) .AND. l <= fc_basis%lmax(i)) THEN
                        DO ii = 1, fc_basis%nshell(i)
                           IF (fc_basis%l(ii, i) == l) THEN
                              ns = ns+1
                              qs_kind%hshift(ii, i) = hshift(ns, l)
                           END IF
                        END DO
                     END IF
                  END DO
               END DO
               ! core density
               ns = grid_atom%nr
               ALLOCATE (qs_kind%rho_fc(ns), qs_kind%drho_fc(ns), qs_kind%tau_fc(ns))
               qs_kind%rho_fc = 0.0_dp
               qs_kind%drho_fc = 0.0_dp
               qs_kind%tau_fc = 0.0_dp
               DO l = 0, 3
                  fc = REAL(2*(2*l+1), KIND=dp)
                  DO i = 1, basis%nbas(l)
                     qs_kind%rho_fc(:) = qs_kind%rho_fc(:)+fc*basis%bf(:, i, l)**2
                     qs_kind%drho_fc(:) = qs_kind%drho_fc(:)+2.0_dp*fc*basis%dbf(:, i, l)*basis%bf(:, i, l)
                     qs_kind%tau_fc(:) = qs_kind%tau_fc(:)+0.5_dp*fc*(basis%dbf(:, i, l)**2+ &
                                                                      REAL(l*(l+1), dp)*basis%bf(:, i, l)**2/grid_atom%rad2(:))
                  END DO
               END DO
!deb
               WRITE (6, *) "FROZEN ENERGY          ", efc
               WRITE (6, *) "FROZEN ENERGY kin      ", ekin
               WRITE (6, *) "FROZEN ENERGY core     ", ecore
               WRITE (6, *) "FROZEN SHIFT  s        ", hshift(:, 0)
!deb
               DEALLOCATE (pmat, hshift)
               !
               CALL release_atom_basis(basis)
               DEALLOCATE (basis)
               !
               IF (ortho > 1.e-10) THEN
                  IF (iounit > 0) THEN
                     WRITE (iounit, "(A,I3,T66,F15.10)") " FROZEN CORE| Basis orthogonality for Kind ", ikind, ortho
                  END IF
               END IF
               !
            END IF

         END DO
      END IF

   END SUBROUTINE set_fc_kind

! **************************************************************************************************

END MODULE qs_fc_methods
