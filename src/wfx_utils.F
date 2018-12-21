!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief  Functions handling the AIM extended wavefunction format.
!> \author Sergey Chulkov, 10.2018
! **************************************************************************************************
MODULE wfx_utils
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_types,                     ONLY: cp_fm_get_submatrix,&
                                              cp_fm_type
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              default_string_length,&
                                              dp
   USE mathconstants,                   ONLY: pi
   USE orbital_pointers,                ONLY: nco,&
                                              nso
!, so
   USE orbital_transformation_matrices, ONLY: orbtramat
   USE particle_types,                  ONLY: particle_type
   USE periodic_table,                  ONLY: get_ptable_info
   USE qs_energy_types,                 ONLY: qs_energy_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE util,                            ONLY: find_boundary_binary
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE
   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'wfx_utils'
   LOGICAL, PARAMETER                   :: debug_this_module = .TRUE.

   PUBLIC :: write_mos_wfx

   TYPE uncontracted_basis_set_map
      ! number of uncontracted Cartesian basis functions
      INTEGER                                            :: ncgf
      ! number of primitive basis functions
      INTEGER                                            :: npgf
      !> (1:npgf); the index of the set which contains the given primitive GTO
      !> basis function within the contracted orbital basis set
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: iset
      !> (1:npgf); index of the primitive GTO within the set in contracted basis set
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ipgf_set
      !> Angular momentum of the primitive GTO
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: l_pgf
      !> (1:npgf); number of shells in contracted basis set which contain
      !> the primitive GTO with the given exponent and angular momentum
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nshell_pgf
      !> (1:nshell_pgf(i), 1:npgf); list of shell indices which contain
      !> the primitive GTO with the given exponent and angular momentum
      INTEGER, ALLOCATABLE, DIMENSION(:, :)               :: ishell_pgf
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)         :: orbtramat_uncontract
   END TYPE uncontracted_basis_set_map

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mos ...
!> \param qs_energy ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param dft_control ...
!> \param print_section ...
! **************************************************************************************************
   SUBROUTINE write_mos_wfx(mos, qs_energy, qs_kind_set, particle_set, dft_control, print_section)
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_energy_type), POINTER                      :: qs_energy
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(section_vals_type), POINTER                   :: print_section

      CHARACTER(len=*), PARAMETER :: print_section_name = "MOS_WFX", routineN = 'write_mos_wfx', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=2)                                   :: element_symbol
      CHARACTER(len=default_path_length)                 :: filename
      CHARACTER(len=default_string_length)               :: atomic_kind_name, fmt_str_es, fmt_str_i, &
                                                            fmt_str_max_es
      INTEGER :: atomic_number, handle, iatom, ikind, imo, imo_alpha_beta, ispin, natoms, ncgf, &
         ndigits, nkinds, npgf, nsgf, nspins, total_nuclear_charge, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: i1d
      LOGICAL                                            :: ghost
      REAL(kind=dp)                                      :: effective_nuclear_charge
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:)           :: r1d
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: mo_coeff_sgf, mo_coeff_ucgf
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(uncontracted_basis_set_map), ALLOCATABLE, &
         DIMENSION(:)                                    :: pbasis_map_set

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      IF (BTEST(cp_print_key_should_output(logger%iter_info, print_section, print_section_name), cp_p_file)) THEN
         unit_nr = cp_print_key_unit_nr(logger, print_section, print_section_name, &
                                        extension='.wfx', file_status='REPLACE', fout=filename)

         natoms = SIZE(particle_set)
         nkinds = SIZE(qs_kind_set)
         nspins = dft_control%nspins

         CALL get_qs_kind_set(qs_kind_set, ncgf=ncgf, nsgf=nsgf)

         DO ispin = 1, nspins
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
!            CALL cp_fm_get_info(mos(ispin)%mo_set%mo_coeff, nrow_global=nao, ncol_global=nmo)
         END DO

         IF (unit_nr > 0) THEN
            CALL section_vals_val_get(print_section, print_section_name//'%NDIGITS', i_val=ndigits)
            IF (ndigits < 2) ndigits = 2

            ! construct format string based on the number of significiant digits
            !      ndigits-1
            ! ' +d.dddddddddE+000'
            !  12 3         45678   extra characters
            WRITE (fmt_str_es, '(A2,I0,A1,I0,A2)') "ES", ndigits+8, ".", ndigits-1, "E3"
            fmt_str_es = ADJUSTL(fmt_str_es)
            iatom = default_string_length/(ndigits+8)
            IF (iatom < 1) iatom = 1
            WRITE (fmt_str_max_es, '(I0)') iatom
            fmt_str_max_es = TRIM(ADJUSTL(fmt_str_max_es))//fmt_str_es

            ! title
            WRITE (unit_nr, '(A,/,A,/,A)') "<Title>", TRIM(ADJUSTL(filename)), "</Title>"

            ! keywords
            WRITE (unit_nr, '(A,/,A,/,A)') "<Keywords>", "GTO", "</Keywords>"

            ! number of nuclei
            WRITE (unit_nr, '(A,/,I0,/,A)') "<Number of Nuclei>", SIZE(particle_set), "</Number of Nuclei>"

            ! atomic labels; label atoms with the name of the corresponding atomic kind + atomic number;
            ! AIMAll expects nuclear names to be unique
            WRITE (unit_nr, '(A)') "<Nuclear Names>"
            DO iatom = 1, SIZE(particle_set)
               CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                    name=atomic_kind_name)
               WRITE (unit_nr, '(A,I0)') TRIM(ADJUSTL(atomic_kind_name)), iatom
            END DO
            WRITE (unit_nr, '(A)') "<Nuclear Names>"

            ! atomic numbers
            WRITE (unit_nr, '(A)') "<Atomic Numbers>"
            total_nuclear_charge = 0
            DO iatom = 1, SIZE(particle_set)
               CALL get_atomic_kind(atomic_kind=particle_set(iatom)%atomic_kind, &
                                    kind_number=ikind, element_symbol=element_symbol)
               CALL get_ptable_info(element_symbol, number=atomic_number)
               CALL get_qs_kind(qs_kind_set(ikind), ghost=ghost)
               IF (ghost) atomic_number = 0
               total_nuclear_charge = total_nuclear_charge+atomic_number
               WRITE (unit_nr, '(I0)') atomic_number
            END DO
            WRITE (unit_nr, '(A)') "</Atomic Numbers>"

            ! effective nuclear charge = atomic number - number of electrons modelled using pseudo-potentials
            WRITE (unit_nr, '(A)') "<Nuclear Charges>"
            DO iatom = 1, SIZE(particle_set)
               CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
               CALL get_qs_kind(qs_kind_set(ikind), zeff=effective_nuclear_charge)
               WRITE (unit_nr, '('//TRIM(fmt_str_es)//')') effective_nuclear_charge
            END DO
            WRITE (unit_nr, '(A)') "</Nuclear Charges>"

            ! cartesian coordinates (in a.u.)
            WRITE (unit_nr, '(A)') "<Nuclear Cartesian Coordinates>"
            DO iatom = 1, SIZE(particle_set)
               WRITE (unit_nr, '(3'//TRIM(fmt_str_es)//')') particle_set(iatom)%r(:)
            END DO
            WRITE (unit_nr, '(A)') "</Nuclear Cartesian Coordinates>"

            ! number of "perturbations"; should be 0 for "non-magnetic wave-functions"
            WRITE (unit_nr, '(A,/,I0,/,A)') "<Number of Perturbations>", 0, "</Number of Perturbations>"

            ! number of electrons
            WRITE (unit_nr, '(A)') "<Number of Electrons>"
            IF (dft_control%lsd) THEN
               WRITE (unit_nr, '(I0)') mos(1)%mo_set%nelectron+mos(2)%mo_set%nelectron
            ELSE
               WRITE (unit_nr, '(I0)') mos(1)%mo_set%nelectron
            END IF
            WRITE (unit_nr, '(A)') "</Number of Electrons>"

            WRITE (unit_nr, '(A)') "<Number of Alpha Electrons>"
            IF (dft_control%lsd) THEN
               WRITE (unit_nr, '(I0)') mos(1)%mo_set%nelectron
            ELSE
               WRITE (unit_nr, '(I0)') mos(1)%mo_set%nelectron/2
            END IF
            WRITE (unit_nr, '(A)') "</Number of Alpha Electrons>"

            WRITE (unit_nr, '(A)') "<Number of Beta Electrons>"
            IF (dft_control%lsd) THEN
               WRITE (unit_nr, '(I0)') mos(2)%mo_set%nelectron
            ELSE
               WRITE (unit_nr, '(I0)') mos(1)%mo_set%nelectron/2
            END IF
            WRITE (unit_nr, '(A)') "</Number of Beta Electrons>"

            WRITE (unit_nr, '(A)') "<Number of Core Electrons>"
            IF (dft_control%lsd) THEN
               WRITE (unit_nr, '(I0)') total_nuclear_charge- &
                  (dft_control%charge+mos(1)%mo_set%nelectron+mos(2)%mo_set%nelectron)
            ELSE
               WRITE (unit_nr, '(I0)') total_nuclear_charge-(dft_control%charge+mos(1)%mo_set%nelectron)
            END IF
            WRITE (unit_nr, '(A)') "</Number of Core Electrons>"

            ! charge of the system
            WRITE (unit_nr, '(A,/,'//TRIM(fmt_str_es)//',/,A)') &
               "<Net Charge>", REAL(dft_control%charge, kind=dp), "</Net Charge>"

            ! multiplicity
            WRITE (unit_nr, '(A,/,I0,/,A)') "<Electronic Spin Multiplicity>", &
               dft_control%multiplicity, "</Electronic Spin Multiplicity>"

            ! number of Occupied Molecular Orbitals
            WRITE (unit_nr, '(A)') "<Number of Occupied Molecular Orbitals>"
            IF (dft_control%lsd) THEN
               WRITE (unit_nr, '(I0)') mos(1)%mo_set%homo+mos(2)%mo_set%homo
            ELSE
               WRITE (unit_nr, '(I0)') mos(1)%mo_set%homo
            END IF
            WRITE (unit_nr, '(A)') "</Number of Occupied Molecular Orbitals>"

            ! generate mapping between the contracted and primitive basis sets
            ALLOCATE (pbasis_map_set(nkinds))

            DO ikind = 1, nkinds
               NULLIFY (orb_basis_set)
               CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)

               CALL uncontract_basis_set(qs_kind_set, orb_basis_set, pbasis_map_set(ikind))
            END DO

            ! <Number of Primitives>, <Primitive Centers>, <Primitive Types>, <Primitive Exponents>
            npgf = get_uncontracted_nao(particle_set, pbasis_map_set) ! -> nao, ncgf ?
            WRITE (unit_nr, '(A,/,I0,/,A)') "<Number of Primitives>", npgf, "</Number of Primitives>"
            ALLOCATE (i1d(npgf))

            CALL get_uncontracted_cgf_centers(particle_set, pbasis_map_set, i1d)
            WRITE (fmt_str_i, '(I0)') MAXVAL(i1d)
            iatom = LEN_TRIM(ADJUSTL(fmt_str_i))
            WRITE (fmt_str_i, '(A,I0,A)') "(10I", iatom+1, ")"
            WRITE (unit_nr, '(A)') "<Primitive Centers>"
            WRITE (unit_nr, fmt_str_i) i1d
            WRITE (unit_nr, '(A)') "</Primitive Centers>"

            CALL get_uncontracted_cgf_types(qs_kind_set, particle_set, pbasis_map_set, i1d)
            WRITE (fmt_str_i, '(I0)') MAXVAL(i1d)
            iatom = LEN_TRIM(ADJUSTL(fmt_str_i))
            WRITE (fmt_str_i, '(A,I0,A)') "(10I", iatom+1, ")"
            WRITE (unit_nr, '(A)') "<Primitive Types>"
            WRITE (unit_nr, fmt_str_i) i1d
            WRITE (unit_nr, '(A)') "</Primitive Types>"

            DEALLOCATE (i1d)
            ALLOCATE (r1d(npgf))

            CALL get_uncontracted_cgf_zet(qs_kind_set, particle_set, pbasis_map_set, r1d)
            WRITE (unit_nr, '(A)') "<Primitive Exponents>"
            WRITE (unit_nr, '('//TRIM(fmt_str_max_es)//')') r1d
            WRITE (unit_nr, '(A)') "</Primitive Exponents>"

            DEALLOCATE (r1d)

            ! total energy
            WRITE (unit_nr, '(A)') "<Energy = T + Vne + Vee + Vnn>"
            WRITE (unit_nr, '('//TRIM(fmt_str_es)//')') qs_energy%total
            WRITE (unit_nr, '(A)') "</Energy = T + Vne + Vee + Vnn>"

            ! virial ratio -E/T + 1
            WRITE (unit_nr, '(A)') "<Virial Ratio (-V/T)>"
            WRITE (unit_nr, '('//TRIM(fmt_str_es)//')') 1.0_dp-qs_energy%total/qs_energy%kinetic
            WRITE (unit_nr, '(A)') "/Virial Ratio (-V/T)"

            ! <Molecular Orbital Occupation Numbers>
            WRITE (unit_nr, '(A)') "<Molecular Orbital Occupation Numbers>"
            DO ispin = 1, nspins
               WRITE (unit_nr, '('//TRIM(fmt_str_es)//')') mos(ispin)%mo_set%occupation_numbers(1:mos(ispin)%mo_set%nmo)
            END DO
            WRITE (unit_nr, '(A)') "</Molecular Orbital Occupation Numbers>"

            ! <Molecular Orbital Energies>
            WRITE (unit_nr, '(A)') "<Molecular Orbital Energies>"
            DO ispin = 1, nspins
               WRITE (unit_nr, '('//TRIM(fmt_str_es)//')') mos(ispin)%mo_set%eigenvalues(1:mos(ispin)%mo_set%nmo)
            END DO
            WRITE (unit_nr, '(A)') "</Molecular Orbital Energies>"

            ! <Molecular Orbital Spin Types>
            WRITE (unit_nr, '(A)') "<Molecular Orbital Spin Types>"
            SELECT CASE (nspins)
            CASE (1)
               DO imo = 1, mos(1)%mo_set%nmo
                  WRITE (unit_nr, '(A)') "Alpha and Beta"
               END DO
            CASE (2)
               DO imo = 1, mos(1)%mo_set%nmo
                  WRITE (unit_nr, '(A)') "Alpha"
               END DO
               DO imo = 1, mos(2)%mo_set%nmo
                  WRITE (unit_nr, '(A)') "Beta"
               END DO
            END SELECT
            WRITE (unit_nr, '(A)') "</Molecular Orbital Spin Types>"

            ! <Molecular Orbital Primitive Coefficients>
            WRITE (unit_nr, '(A)') "<Molecular Orbital Primitive Coefficients>"
         END IF

         imo_alpha_beta = 0
         DO ispin = 1, nspins
            ALLOCATE (mo_coeff_sgf(mos(ispin)%mo_set%nao, mos(ispin)%mo_set%nmo))
            CALL cp_fm_get_submatrix(mos(ispin)%mo_set%mo_coeff, mo_coeff_sgf)

            IF (unit_nr > 0) THEN
               ALLOCATE (mo_coeff_ucgf(npgf, mos(ispin)%mo_set%nmo))
               CALL uncontract_mocoeff(mo_coeff_sgf, mo_coeff_ucgf, qs_kind_set, particle_set, pbasis_map_set)

               DO imo = 1, mos(ispin)%mo_set%nmo
                  imo_alpha_beta = imo_alpha_beta+1
                  WRITE (unit_nr, '(A,/,I0,/,A)') "<MO Number>", imo_alpha_beta, "</MO Number>"
                  WRITE (unit_nr, '('//TRIM(fmt_str_max_es)//')') mo_coeff_ucgf(:, imo)
               END DO

               DEALLOCATE (mo_coeff_ucgf)
            END IF

            DEALLOCATE (mo_coeff_sgf)
         END DO

         IF (unit_nr > 0) THEN
            WRITE (unit_nr, '(A)') "</Molecular Orbital Primitive Coefficients>"

            ! release mapping
            DO ikind = 1, nkinds
               NULLIFY (orb_basis_set)
               CALL uncontract_basis_set(qs_kind_set, orb_basis_set, pbasis_map_set(ikind))
            END DO
            DEALLOCATE (pbasis_map_set)
         END IF

         CALL cp_print_key_finished_output(unit_nr, logger, print_section, print_section_name)
      END IF

      CALL timestop(handle)
   END SUBROUTINE write_mos_wfx

! **************************************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param orb_basis_set ...
!> \param pbasis_map ...
! **************************************************************************************************
   SUBROUTINE uncontract_basis_set(qs_kind_set, orb_basis_set, pbasis_map)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set
      TYPE(uncontracted_basis_set_map), INTENT(inout)    :: pbasis_map

      REAL(kind=dp), PARAMETER :: threshold = 16.0_dp*EPSILON(0.0_dp)

      INTEGER :: icgf, icgf_orb, icgf_pgf, ipgf, ipgf_set, iset, ishell, lmax_set, lmin_set, &
         lshell, ncgf, nco_pgf, npgf, npgf_list, nshell_pgf, nshell_set, orb_ncgf, orb_nset
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: cgf_offset, ipgf_list
      INTEGER, DIMENSION(:), POINTER                     :: orb_lmax, orb_lmin, orb_npgf, orb_nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: orb_first_cgf, orb_l, orb_last_cgf
      REAL(kind=dp)                                      :: expzet, gcca, prefac, zeta
      REAL(kind=dp), DIMENSION(:), POINTER               :: orb_norm_cgf
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: orb_zet
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: orb_gcc

      IF (ASSOCIATED(orb_basis_set)) THEN
         CALL get_gto_basis_set(orb_basis_set, nset=orb_nset, nshell=orb_nshell, npgf=orb_npgf, &
                                l=orb_l, lmin=orb_lmin, lmax=orb_lmax, gcc=orb_gcc)
         CALL get_gto_basis_set(orb_basis_set, first_cgf=orb_first_cgf, last_cgf=orb_last_cgf, zet=orb_zet, norm_cgf=orb_norm_cgf)

         ! compute the total number of primitive Gaussian functions for all angular momenta
         npgf = 0
         DO iset = 1, orb_nset
            npgf = npgf+orb_npgf(iset)*(orb_lmax(iset)-orb_lmin(iset)+1)
         END DO

         ! map primitive GTO-s between uncontracted and contracted basis sets
         pbasis_map%npgf = npgf
         ALLOCATE (pbasis_map%iset(npgf), pbasis_map%ipgf_set(npgf), pbasis_map%l_pgf(npgf))
         ALLOCATE (pbasis_map%nshell_pgf(npgf), pbasis_map%ishell_pgf(MAXVAL(orb_nshell(1:orb_nset)), npgf))
         pbasis_map%ishell_pgf = 0

         ipgf = 0
         DO iset = 1, orb_nset
            nshell_set = orb_nshell(iset)
            lmin_set = orb_lmin(iset)
            lmax_set = orb_lmax(iset)
            DO ipgf_set = 1, orb_npgf(iset)
               DO lshell = lmin_set, lmax_set
                  ipgf = ipgf+1
                  pbasis_map%iset(ipgf) = iset
                  pbasis_map%ipgf_set(ipgf) = ipgf_set
                  pbasis_map%l_pgf(ipgf) = lshell
                  nshell_pgf = 0

                  DO ishell = 1, nshell_set
                     IF (orb_l(ishell, iset) == lshell) THEN
                        IF (ABS(orb_gcc(ipgf_set, ishell, iset)) >= threshold) THEN
                           nshell_pgf = nshell_pgf+1
                           pbasis_map%ishell_pgf(nshell_pgf, ipgf) = ishell
                        END IF
                     END IF
                  END DO

                  pbasis_map%nshell_pgf(ipgf) = nshell_pgf
               END DO
            END DO
         END DO

         IF (debug_this_module) THEN
            CPASSERT(ipgf == npgf)
         END IF

         ncgf = 0
         DO ipgf = 1, npgf
            IF (pbasis_map%nshell_pgf(ipgf) > 0) &
               ncgf = ncgf+nco(pbasis_map%l_pgf(ipgf))
         END DO
         pbasis_map%ncgf = ncgf

         CALL get_qs_kind_set(qs_kind_set, ncgf=orb_ncgf)
         ALLOCATE (pbasis_map%orbtramat_uncontract(ncgf, orb_ncgf))
         pbasis_map%orbtramat_uncontract = 0.0_dp

         ALLOCATE (cgf_offset(npgf))
         icgf = 0
         DO ipgf = 1, npgf
            cgf_offset(ipgf) = icgf
            IF (pbasis_map%nshell_pgf(ipgf) > 0) THEN
               icgf = icgf+nco(pbasis_map%l_pgf(ipgf))
            END IF
         END DO

         IF (debug_this_module) THEN
            CPASSERT(icgf == ncgf)
         END IF

         icgf_orb = 0
         ALLOCATE (ipgf_list(MAXVAL(orb_npgf)))

         DO iset = 1, orb_nset
            nshell_set = orb_nshell(iset)

            DO ishell = 1, nshell_set
               lshell = orb_l(ishell, iset)
               expzet = 0.25_dp*REAL(2*lshell+3, dp)
               prefac = 2.0_dp**lshell*(2.0_dp/pi)**0.75_dp
               DO icgf_pgf = orb_first_cgf(ishell, iset), orb_last_cgf(ishell, iset)
                  icgf_orb = icgf_orb+1
                  CALL get_list_of_uncontracted_pgfs(pbasis_map, iset, ishell, orb_l(ishell, iset), npgf_list, ipgf_list)

                  DO ipgf_set = 1, npgf_list
                     ipgf = ipgf_list(ipgf_set)

                     gcca = orb_gcc(pbasis_map%ipgf_set(ipgf), ishell, iset)
                     zeta = orb_zet(pbasis_map%ipgf_set(ipgf), iset)

                     icgf = cgf_offset(ipgf)
                     IF (pbasis_map%nshell_pgf(ipgf) > 0) THEN
                        nco_pgf = nco(pbasis_map%l_pgf(ipgf))
                        DO lshell = 1, pbasis_map%nshell_pgf(ipgf)
                           IF (pbasis_map%ishell_pgf(lshell, ipgf) == ishell) THEN
                              icgf = icgf+icgf_pgf
                              EXIT
                           ELSE
                              icgf = icgf+nco_pgf
                           END IF
                        END DO
                     END IF

!                     pbasis_map%orbtramat_uncontract(icgf, icgf_orb) = 1.0_dp/(prefac*zeta**expzet*gcca*orb_norm_cgf(icgf_orb))
                     pbasis_map%orbtramat_uncontract(icgf, icgf_orb) = orb_norm_cgf(icgf_orb)/(prefac*zeta**expzet*gcca)
                  END DO
               END DO
            END DO
         END DO
         DEALLOCATE (cgf_offset, ipgf_list)
      ELSE
         pbasis_map%ncgf = 0
         pbasis_map%npgf = 0

         IF (ALLOCATED(pbasis_map%iset)) DEALLOCATE (pbasis_map%iset)
         IF (ALLOCATED(pbasis_map%ipgf_set)) DEALLOCATE (pbasis_map%ipgf_set)
         IF (ALLOCATED(pbasis_map%l_pgf)) DEALLOCATE (pbasis_map%l_pgf)
         IF (ALLOCATED(pbasis_map%nshell_pgf)) DEALLOCATE (pbasis_map%nshell_pgf)
         IF (ALLOCATED(pbasis_map%ishell_pgf)) DEALLOCATE (pbasis_map%ishell_pgf)
         IF (ALLOCATED(pbasis_map%orbtramat_uncontract)) DEALLOCATE (pbasis_map%orbtramat_uncontract)
      END IF
   END SUBROUTINE uncontract_basis_set

! **************************************************************************************************
!> \brief ...
!> \param pbasis_map ...
!> \param iset_orb ...
!> \param ishell_orb ...
!> \param l_orb ...
!> \param npgf ...
!> \param ipgf_list ...
! **************************************************************************************************
   SUBROUTINE get_list_of_uncontracted_pgfs(pbasis_map, iset_orb, ishell_orb, l_orb, npgf, ipgf_list)
      TYPE(uncontracted_basis_set_map), INTENT(in)       :: pbasis_map
      INTEGER, INTENT(in)                                :: iset_orb, ishell_orb, l_orb
      INTEGER, INTENT(out)                               :: npgf
      INTEGER, DIMENSION(:), INTENT(out)                 :: ipgf_list

      INTEGER                                            :: ipgf, ipgf_first, ipgf_last, ishell

      npgf = 0

      CALL find_boundary_binary(pbasis_map%iset, iset_orb, ipgf_first, ipgf_last)
      DO ipgf = ipgf_first, ipgf_last
         IF (pbasis_map%l_pgf(ipgf) == l_orb .AND. pbasis_map%nshell_pgf(ipgf) > 0) THEN

            ! the size of the ishell_pgf(:, ipgf) array is very small (typically <= 5 items),
            ! so locate() has no advantange against the naive linear search algorithm
            DO ishell = 1, pbasis_map%nshell_pgf(ipgf)
               IF (pbasis_map%ishell_pgf(ishell, ipgf) == ishell_orb) THEN
                  IF (debug_this_module) THEN
                     CPASSERT(npgf < SIZE(ipgf_list))
                  END IF

                  npgf = npgf+1
                  ipgf_list(npgf) = ipgf
                  EXIT
               END IF
            END DO
         END IF
      END DO
   END SUBROUTINE get_list_of_uncontracted_pgfs

! **************************************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param pbasis_map_set ...
!> \return ...
! **************************************************************************************************
   FUNCTION get_uncontracted_nao(particle_set, pbasis_map_set) RESULT(nao)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(uncontracted_basis_set_map), DIMENSION(:), &
         INTENT(in)                                      :: pbasis_map_set
      INTEGER                                            :: nao

      INTEGER                                            :: iatom, ikind

      ! compute the number of atomic orbitals in terms of uncontracted Cartesian basis set
      nao = 0
      DO iatom = 1, SIZE(particle_set)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         nao = nao+pbasis_map_set(ikind)%ncgf
      END DO
   END FUNCTION get_uncontracted_nao

! **************************************************************************************************
!> \brief ...
!> \param pbasis_map_set ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION maxval_ncgf(pbasis_map_set) RESULT(ncgf)
      TYPE(uncontracted_basis_set_map), DIMENSION(:), &
         INTENT(in)                                      :: pbasis_map_set
      INTEGER                                            :: ncgf

      INTEGER                                            :: ikind

      ncgf = 0
      DO ikind = 1, SIZE(pbasis_map_set)
         ncgf = MAX(pbasis_map_set(ikind)%ncgf, ncgf)
      END DO
   END FUNCTION maxval_ncgf

! **************************************************************************************************
!> \brief ...
!> \param particle_set ...
!> \param pbasis_map_set ...
!> \param cgf_centers ...
! **************************************************************************************************
   SUBROUTINE get_uncontracted_cgf_centers(particle_set, pbasis_map_set, cgf_centers)
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(uncontracted_basis_set_map), DIMENSION(:), &
         INTENT(in)                                      :: pbasis_map_set
      INTEGER, DIMENSION(:), INTENT(out)                 :: cgf_centers

      INTEGER                                            :: iatom, ikind, nao, ncgf_atom

      IF (debug_this_module) THEN
         nao = get_uncontracted_nao(particle_set, pbasis_map_set)
         CPASSERT(nao == SIZE(cgf_centers))
      END IF

      nao = 0
      DO iatom = 1, SIZE(particle_set)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         ncgf_atom = pbasis_map_set(ikind)%ncgf
         cgf_centers(nao+1:nao+ncgf_atom) = iatom
         nao = nao+ncgf_atom
      END DO

   END SUBROUTINE get_uncontracted_cgf_centers

! **************************************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param pbasis_map_set ...
!> \param cgf_types ...
! **************************************************************************************************
   SUBROUTINE get_uncontracted_cgf_types(qs_kind_set, particle_set, pbasis_map_set, cgf_types)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(uncontracted_basis_set_map), DIMENSION(:), &
         INTENT(in)                                      :: pbasis_map_set
      INTEGER, DIMENSION(:), INTENT(out)                 :: cgf_types

      INTEGER                                            :: iatom, icgf_orb, ikind, ipgf, ipgf_set, &
                                                            iset, ishell, ncgf, ncgf_atom, nkinds
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: cgf_types_kind
      INTEGER, DIMENSION(:), POINTER                     :: orb_lx, orb_ly, orb_lz
      INTEGER, DIMENSION(:, :), POINTER                  :: orb_first_cgf, orb_last_cgf
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      IF (debug_this_module) THEN
         ncgf = get_uncontracted_nao(particle_set, pbasis_map_set)
         CPASSERT(ncgf == SIZE(cgf_types))
      END IF

      nkinds = SIZE(pbasis_map_set)
      ALLOCATE (cgf_types_kind(maxval_ncgf(pbasis_map_set), nkinds))

      DO ikind = 1, nkinds
         NULLIFY (orb_basis_set)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)

         IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(orb_basis_set, first_cgf=orb_first_cgf, last_cgf=orb_last_cgf, &
                                   lx=orb_lx, ly=orb_ly, lz=orb_lz)

            ncgf = 0
            DO ipgf = 1, pbasis_map_set(ikind)%npgf
               iset = pbasis_map_set(ikind)%iset(ipgf)
               ipgf_set = pbasis_map_set(ikind)%ipgf_set(ipgf)

               IF (pbasis_map_set(ikind)%nshell_pgf(ipgf) > 0) THEN
                  ishell = pbasis_map_set(ikind)%ishell_pgf(1, ipgf)

                  DO icgf_orb = orb_first_cgf(ishell, iset), orb_last_cgf(ishell, iset)
                     ncgf = ncgf+1
                     cgf_types_kind(ncgf, ikind) = wfx_primitive_type(orb_lx(icgf_orb), orb_ly(icgf_orb), orb_lz(icgf_orb))
                  END DO
               END IF
            END DO

            IF (debug_this_module) THEN
               CPASSERT(ncgf == pbasis_map_set(ikind)%ncgf)
            END IF
         END IF
      END DO

      ncgf = 0
      DO iatom = 1, SIZE(particle_set)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         ncgf_atom = pbasis_map_set(ikind)%ncgf
         cgf_types(ncgf+1:ncgf+ncgf_atom) = cgf_types_kind(1:ncgf_atom, ikind)
         ncgf = ncgf+ncgf_atom
      END DO

      DEALLOCATE (cgf_types_kind)
   END SUBROUTINE get_uncontracted_cgf_types

! **************************************************************************************************
!> \brief ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param pbasis_map_set ...
!> \param cgf_zet ...
! **************************************************************************************************
   SUBROUTINE get_uncontracted_cgf_zet(qs_kind_set, particle_set, pbasis_map_set, cgf_zet)
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(uncontracted_basis_set_map), DIMENSION(:), &
         INTENT(in)                                      :: pbasis_map_set
      REAL(kind=dp), DIMENSION(:), INTENT(out)           :: cgf_zet

      INTEGER                                            :: iatom, ikind, ipgf, ipgf_set, iset, &
                                                            ncgf, ncgf_atom, nco_pgf, nkinds
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: cgf_zet_kind
      REAL(kind=dp), DIMENSION(:, :), POINTER            :: orb_zet
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      IF (debug_this_module) THEN
         ncgf = get_uncontracted_nao(particle_set, pbasis_map_set)
         CPASSERT(ncgf == SIZE(cgf_zet))
      END IF

      nkinds = SIZE(pbasis_map_set)
      ALLOCATE (cgf_zet_kind(maxval_ncgf(pbasis_map_set), nkinds))

      DO ikind = 1, nkinds
         NULLIFY (orb_basis_set)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)

         IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(orb_basis_set, zet=orb_zet)

            ncgf = 0
            DO ipgf = 1, pbasis_map_set(ikind)%npgf
               IF (pbasis_map_set(ikind)%nshell_pgf(ipgf) > 0) THEN
                  iset = pbasis_map_set(ikind)%iset(ipgf)
                  ipgf_set = pbasis_map_set(ikind)%ipgf_set(ipgf)
                  nco_pgf = nco(pbasis_map_set(ikind)%l_pgf(ipgf))

                  cgf_zet_kind(ncgf+1:ncgf+nco_pgf, ikind) = orb_zet(ipgf_set, iset)

                  ncgf = ncgf+nco_pgf
               END IF
            END DO

            IF (debug_this_module) THEN
               CPASSERT(ncgf == pbasis_map_set(ikind)%ncgf)
            END IF
         END IF
      END DO

      ncgf = 0
      DO iatom = 1, SIZE(particle_set)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         ncgf_atom = pbasis_map_set(ikind)%ncgf
         cgf_zet(ncgf+1:ncgf+ncgf_atom) = cgf_zet_kind(1:ncgf_atom, ikind)
         ncgf = ncgf+ncgf_atom
      END DO

      DEALLOCATE (cgf_zet_kind)
   END SUBROUTINE get_uncontracted_cgf_zet

! **************************************************************************************************
!> \brief ...
!> \param lx ...
!> \param ly ...
!> \param lz ...
!> \return ...
! **************************************************************************************************
   PURE FUNCTION wfx_primitive_type(lx, ly, lz) RESULT(pgf_type)
      INTEGER, INTENT(in)                                :: lx, ly, lz
      INTEGER                                            :: pgf_type

      !  1 S
      !  2 PX,     3 PY,     4 PZ
      !  5 DXX,    6 DYY,    7 DZZ,    8 DXY,    9 DXZ,   10 DYZ
      ! 11 FXXX,  12 FYYY,  13 FZZZ,  14 FXXY,  15 FXXZ,  16 FYYZ,  17 FXYY,  18 FXZZ, 19 FYZZ, 20 FXYZ
      ! 21 GXXXX, 22 GYYYY, 23 GZZZZ, 24 GXXXY, 25 GXXXZ, 26 GXYYY, 27 GYYYZ, 28 GXZZZ
      ! 29 GYZZZ, 30 GXXYY, 31 GXXZZ, 32 GYYZZ, 33 GXXYZ, 34 GXYYZ, 35 GXYZZ
      INTEGER, DIMENSION(0:3, 0:3, 0:3), PARAMETER :: lxlylz_to_wfxtype = RESHAPE(&
         (/ 1,  2,  5, 11,  3,  8, 14, 24,  6, 17, 30,  0, 12, 26,  0,  0,&
            4,  9, 15, 25, 10, 20, 33,  0, 16, 34,  0,  0, 27,  0,  0,  0,&
            7, 18, 31,  0, 19, 35,  0,  0, 32,  0,  0,  0,  0,  0,  0,  0,&
           13, 28,  0,  0, 29,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0/), (/4, 4, 4/))

      INTEGER                                            :: l

      l = lx+ly+lz

      IF (l >= 0) THEN
         ! conventional enumeration scheme (l <= 4)
         IF (l <= 3) THEN
            pgf_type = lxlylz_to_wfxtype(lx, ly, lz)
         ELSE IF (l == 4) THEN
            IF (lx == 4) THEN
               pgf_type = 21 ! GXXXX
            ELSE IF (ly == 4) THEN
               pgf_type = 22 ! GYYYY
            ELSE IF (lz == 4) THEN
               pgf_type = 23 ! GZZZZ
            ELSE
               pgf_type = lxlylz_to_wfxtype(lx, ly, lz)
            END IF
         ELSE
            ! alternative enumeration scheme is in use starting from H-type GTOs (l >= 5):
            ! DO LX = 0, L
            !    DO LY = 0, (L-LX)
            !       LZ = L - LX - LY
            ! GTO index for the given angular momentum l : 0 .. (l+1)(l+2)/2-1
            pgf_type = ly
            IF (lx > 0) pgf_type = pgf_type+(2*l+3-lx)*lx/2

            ! add the number of Cartesian GTOs with angular momenta <l [\sum_{c=0}^{l-1} (c+1)*(c+2)/2];
            ! also add +1, as WFX_TYPE should start from 1 instead of 0
            !
            ! + \sum_{c=5}^{l-1} (c+1)*(c+2)/2
            DO WHILE (l > 5)
               pgf_type = pgf_type+l*(l+1)/2
               l = l-1
            END DO
            pgf_type = pgf_type+36 ! 36 = 1 + \sum_{c=0}^{4} (c+1)*(c+2)/2
         END IF
      ELSE
         pgf_type = 0
      END IF

   END FUNCTION wfx_primitive_type

! **************************************************************************************************
!> \brief ...
!> \param mo_coeff_sgf ...
!> \param mo_coeff_ucgf ...
!> \param qs_kind_set ...
!> \param particle_set ...
!> \param pbasis_map_set ...
! **************************************************************************************************
   SUBROUTINE uncontract_mocoeff(mo_coeff_sgf, mo_coeff_ucgf, qs_kind_set, particle_set, pbasis_map_set)
      REAL(kind=dp), DIMENSION(:, :), INTENT(in)         :: mo_coeff_sgf
      REAL(kind=dp), DIMENSION(:, :), INTENT(out)        :: mo_coeff_ucgf
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(uncontracted_basis_set_map), DIMENSION(:), &
         INTENT(in)                                      :: pbasis_map_set

      INTEGER :: iatom, icgf, ikind, iset, isgf, ishell, lshell, ncgf_offset, ncgf_orb_offset, &
         nmo, orb_ncgf, orb_ncgf_atom, orb_nset, orb_nsgf
      INTEGER, DIMENSION(:), POINTER                     :: orb_nshell
      INTEGER, DIMENSION(:, :), POINTER                  :: orb_l
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: mo_coeff_cgf
      TYPE(gto_basis_set_type), POINTER                  :: orb_basis_set

      nmo = SIZE(mo_coeff_sgf, 2)

      CALL get_qs_kind_set(qs_kind_set, ncgf=orb_ncgf, nsgf=orb_nsgf)
      IF (debug_this_module) THEN
         CPASSERT(orb_nsgf == SIZE(mo_coeff_sgf, 1))
      END IF

      ALLOCATE (mo_coeff_cgf(orb_ncgf, nmo))

      ! Transform spherical MOs to Cartesian MOs
      icgf = 1
      isgf = 1
      ncgf_offset = 0
      ncgf_orb_offset = 0
      DO iatom = 1, SIZE(particle_set)
         NULLIFY (orb_basis_set)
         CALL get_atomic_kind(particle_set(iatom)%atomic_kind, kind_number=ikind)
         CALL get_qs_kind(qs_kind_set(ikind), basis_set=orb_basis_set)
         IF (ASSOCIATED(orb_basis_set)) THEN
            CALL get_gto_basis_set(gto_basis_set=orb_basis_set, nset=orb_nset, nshell=orb_nshell, l=orb_l, ncgf=orb_ncgf_atom)

            DO iset = 1, orb_nset
               DO ishell = 1, orb_nshell(iset)
                  lshell = orb_l(ishell, iset)

                  CALL dgemm("T", "N", nco(lshell), nmo, nso(lshell), 1.0_dp, &
                             orbtramat(lshell)%s2c, nso(lshell), &
                             mo_coeff_sgf(isgf, 1), orb_nsgf, 0.0_dp, &
                             mo_coeff_cgf(icgf, 1), orb_ncgf)

                  icgf = icgf+nco(lshell)
                  isgf = isgf+nso(lshell)
               END DO
            END DO

            ! transform contraction coefficients from spherical GTO to uncontracted GTO basis set
            ! using a single matrix-matrix multiplication (sgf -> ucgf, as opposite to sgf -> cgf -> ucgf)
            CALL dgemm("N", "N", pbasis_map_set(ikind)%ncgf, nmo, orb_ncgf_atom, 1.0_dp, &
                       pbasis_map_set(ikind)%orbtramat_uncontract(1, 1), &
                       pbasis_map_set(ikind)%ncgf, &
                       mo_coeff_cgf(ncgf_orb_offset+1, 1), orb_ncgf, &
                       0.0_dp, mo_coeff_ucgf(ncgf_offset+1, 1), SIZE(mo_coeff_ucgf, 1))
            ncgf_offset = ncgf_offset+pbasis_map_set(ikind)%ncgf
            ncgf_orb_offset = ncgf_orb_offset+orb_ncgf_atom
         END IF
      END DO

      IF (debug_this_module) THEN
         CPASSERT(ncgf_offset == SIZE(mo_coeff_ucgf, 1))
         CPASSERT(ncgf_orb_offset == orb_ncgf)
      END IF

      DEALLOCATE (mo_coeff_cgf)

   END SUBROUTINE uncontract_mocoeff
END MODULE wfx_utils
