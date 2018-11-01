!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Frozen Core
! **************************************************************************************************
MODULE qs_fc_methods
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind,&
                                              get_atomic_kind_set
   USE atom_kind_orbitals,              ONLY: calculate_atomic_energies
   USE atom_types,                      ONLY: atom_basis_type,&
                                              init_atom_basis_cgto,&
                                              release_atom_basis
   USE atom_utils,                      ONLY: coulomb_potential_numeric,&
                                              integrate_grid
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: fourpi
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_grid_atom,                    ONLY: grid_atom_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type,&
                                              set_qs_kind
   USE qs_local_rho_types,              ONLY: get_local_rho,&
                                              local_rho_type,&
                                              rhoz_type
   USE qs_rho0_types,                   ONLY: rho0_mpole_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   ! Global parameters (only in this module)
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_fc_methods'

   PUBLIC :: set_fc_kind, set_fc_rho0

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
      LOGICAL                                            :: frozen_core
      REAL(KIND=dp)                                      :: ecore, efc, ekin, exchange, fc, ortho, &
                                                            zcore
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: pmat
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: hshift
      TYPE(atom_basis_type), POINTER                     :: basis
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(gto_basis_set_type), POINTER                  :: fc_basis
      TYPE(qs_kind_type), POINTER                        :: qs_kind

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
            zcore = SUM(pmat)
            ALLOCATE (hshift(nn, 0:3))
            ! needs some programming for exchange
            exchange = 0.0_dp
            CALL calculate_atomic_energies(qs_kind, basis, pmat=pmat, ekin=ekin, ecore=ecore, &
                                           hshift=hshift, ortho=ortho)
            efc = ekin+ecore+exchange
            CALL set_qs_kind(qs_kind, zcore=zcore, efc=efc)
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
            fc = 1.0_dp/SQRT(fourpi)
            qs_kind%rho_fc(:) = fc*qs_kind%rho_fc(:)
            qs_kind%drho_fc(:) = fc*qs_kind%drho_fc(:)
            qs_kind%tau_fc(:) = fc*qs_kind%tau_fc(:)
!deb
            WRITE (6, *) "FROZEN DENSITY         ", integrate_grid(qs_kind%rho_fc(:), grid_atom)
            WRITE (6, *) "FROZEN DENSITY ZCORE   ", zcore
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

   END SUBROUTINE set_fc_kind

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE set_fc_rho0(qs_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'set_fc_rho0', routineP = moduleN//':'//routineN

      INTEGER                                            :: ikind, natom, nkind, nn, ns, iat, nat, iatom
      INTEGER, DIMENSION(:), POINTER                     :: atom_list
      LOGICAL                                            :: frozen_core
      REAL(KIND=dp)                                      :: zcore
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: cpot
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(grid_atom_type), POINTER                      :: grid_atom
      TYPE(local_rho_type), POINTER                      :: local_rho_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(rho0_mpole_type), POINTER                     :: rho0_mpole
      TYPE(rhoz_type), DIMENSION(:), POINTER             :: rhoz_set

      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
      local_rho_set => qs_env%local_rho_set

      nkind = SIZE(qs_kind_set)
      DO ikind = 1, nkind
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind, frozen_core=frozen_core)
         IF (frozen_core) THEN
            CALL get_atomic_kind(atomic_kind_set(ikind), atom_list=atom_list, natom=nat)
            CALL get_qs_kind(qs_kind, grid_atom=grid_atom)
            ns = grid_atom%nr
            ! initialize atomic charge terms for GAPW Hartree one-center terms
            CALL get_local_rho(local_rho_set, &
                               rho0_mpole=rho0_mpole, &
                               rhoz_set=rhoz_set)
            rhoz_set(ikind)%r_coef(1:ns) = rhoz_set(ikind)%r_coef(1:ns)+qs_kind%rho_fc(1:ns)
            rhoz_set(ikind)%dr_coef(1:ns) = rhoz_set(ikind)%dr_coef(1:ns)+ &
                                            grid_atom%rad(1:ns)*qs_kind%drho_fc(1:ns)
            ALLOCATE (cpot(ns))
            CALL coulomb_potential_numeric(cpot, qs_kind%rho_fc(1:ns), grid_atom)
            rhoz_set(ikind)%vr_coef(1:ns) = rhoz_set(ikind)%vr_coef(1:ns)+cpot(1:ns)
            DEALLOCATE (cpot)
            natom = qs_kind%natom
            CALL get_qs_kind(qs_kind, zcore=zcore)
            local_rho_set%rhoz_tot = local_rho_set%rhoz_tot-zcore*REAL(natom, KIND=dp)
            ! reset Qlm_z
            DO iat = 1, nat
               iatom = atom_list(iat)
               rho0_mpole%mp_rho(iatom)%Qlm_z = rho0_mpole%mp_rho(iatom)%Qlm_z + qs_kind%zcore*SQRT(fourpi)
            END DO
!deb
            nn = grid_atom%nr
            WRITE (6, *) "RHOZ Qlm    ", rho0_mpole%mp_rho(:)%Qlm_z/SQRT(fourpi)
            WRITE (6, *) "RHOZ set    ", SUM(rhoz_set(ikind)%r_coef(1:nn)*grid_atom%wr(1:nn))
            WRITE (6, *) "RHOZ Qtot   ", local_rho_set%rhoz_tot
!deb
            !
         END IF

      END DO

   END SUBROUTINE set_fc_rho0

! **************************************************************************************************

END MODULE qs_fc_methods
