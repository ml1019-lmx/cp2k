!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief NEGF based quantum transport calculations
! **************************************************************************************************

MODULE negf_io
   USE cp_files,                        ONLY: close_file,&
                                              open_file
   USE cp_fm_types,                     ONLY: cp_fm_get_info,&
                                              cp_fm_get_submatrix,&
                                              cp_fm_p_type,&
                                              cp_fm_set_submatrix
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_get_lval,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE message_passing,                 ONLY: mp_bcast
   USE particle_list_types,             ONLY: particle_list_type
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_copy,&
                                              pw_scale
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_p_type,&
                                              pw_pool_type
   USE pw_types,                        ONLY: REALDATA3D,&
                                              REALSPACE,&
                                              pw_p_type,&
                                              pw_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: qs_kind_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_io'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   CHARACTER(len=8), PARAMETER, PRIVATE :: restart_magic = "NEGFHMAT"
   INTEGER, PARAMETER, PRIVATE          :: max_bcast_items_dp = 134217728 ! 2^27 = 1GB

   PUBLIC :: negf_read_restart, negf_write_restart
   PUBLIC :: negf_print_density_cube, negf_print_vhartree_cube

CONTAINS

! **************************************************************************************************
!> \brief Read the initial Kohn-Sham matrix of the scattering region from a NEGF restart file.
!>        If restart file does not exist or has an incorrect format, the content of the variable
!>        'h_s' is not altered.
!> \param h_s                Kohn-Sham matrix to be read
!> \param atomlist_s         atoms belonging to the scattering region
!> \param subsys             QuickStep subsystem of the device force environment
!> \param negf_section       NEGF input section
!> \param para_env_global    global parallel environment
!> \param log_unit           output I/O unit
!> \return                   whether or not the restart file has been read successfully
! **************************************************************************************************
   FUNCTION negf_read_restart(h_s, atomlist_s, subsys, negf_section, para_env_global, log_unit) RESULT(is_restarted)
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(inout)    :: h_s
      INTEGER, DIMENSION(:), INTENT(in)                  :: atomlist_s
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(section_vals_type), POINTER                   :: negf_section
      TYPE(cp_para_env_type), POINTER                    :: para_env_global
      INTEGER, INTENT(in)                                :: log_unit
      LOGICAL                                            :: is_restarted

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_read_restart', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=8)                                   :: restart_magic_read
      CHARACTER(len=default_path_length)                 :: filename
      INTEGER :: handle, i, ispin, maxcols, n_rep_val, natoms, natoms_read, ncols, ncols_batch, &
         ncols_read, nparticles, nrows, nrows_read, nspins, nspins_read, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nsgfs, nsgfs_local, nsgfs_read
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: r2d
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      ! generate restart file name
      CALL section_vals_val_get(negf_section, "E_DENSITY_RESTART_FILE_NAME", n_rep_val=n_rep_val)
      IF (n_rep_val > 0) THEN
         CALL section_vals_val_get(negf_section, "E_DENSITY_RESTART_FILE_NAME", c_val=filename)
      ELSE
         is_restarted = .FALSE.
         CALL timestop(handle)
         RETURN
      END IF

      IF (para_env_global%ionode) THEN
         INQUIRE (FILE=filename, exist=is_restarted)
      END IF

      CALL mp_bcast(is_restarted, para_env_global%source, para_env_global%group)

      IF (.NOT. is_restarted) THEN
         CALL cp_warn(__LOCATION__, &
                      "Initial electron density of the scattering region cannot be read from the file '"//TRIM(filename)// &
                      "' which does not exist.")
         CALL timestop(handle)
         RETURN
      END IF

      unit_nr = 0
      IF (para_env_global%ionode) THEN
         CALL open_file(file_name=filename, file_action="READ", file_form="UNFORMATTED", file_status="OLD", unit_number=unit_nr)
         READ (unit_nr) restart_magic_read
         is_restarted = (restart_magic_read == restart_magic)
         IF (.NOT. is_restarted) &
            CALL close_file(unit_number=unit_nr)
      END IF

      CALL mp_bcast(is_restarted, para_env_global%source, para_env_global%group)

      IF (.NOT. is_restarted) THEN
         CALL cp_warn(__LOCATION__, &
                      "Initial electron density of the scattering region cannot be read from the file '"//TRIM(filename)// &
                      "' which is not a NEGF restart file.")
         CALL timestop(handle)
         RETURN
      END IF

      IF (para_env_global%ionode) THEN
         READ (unit_nr) nspins_read, nrows_read, natoms_read
      END IF

      CALL mp_bcast(nspins_read, para_env_global%source, para_env_global%group)
      CALL mp_bcast(nrows_read, para_env_global%source, para_env_global%group)
      CALL mp_bcast(natoms_read, para_env_global%source, para_env_global%group)

      nspins = SIZE(h_s)
      natoms = SIZE(atomlist_s)
      CALL cp_fm_get_info(h_s(1)%matrix, nrow_global=nrows, ncol_global=ncols)

      CPASSERT(nspins_read == nspins)
      CPASSERT(natoms_read == natoms)
      CPASSERT(nrows_read == nrows)
      CPASSERT(nrows == ncols) ! square matrix

      ALLOCATE (nsgfs_local(natoms), nsgfs_read(natoms))

      CALL qs_subsys_get(subsys, particle_set=particle_set, qs_kind_set=qs_kind_set)
      nparticles = SIZE(particle_set)

      ALLOCATE (nsgfs(nparticles))
      CALL get_particle_set(particle_set, qs_kind_set, nsgf=nsgfs)
      DO i = 1, natoms
         nsgfs_local(i) = nsgfs(atomlist_s(i))
      END DO
      DEALLOCATE (nsgfs)

      IF (debug_this_module) THEN
         CPASSERT(SUM(nsgfs_local(:)) == nrows)
      END IF

      IF (para_env_global%ionode) THEN
         READ (unit_nr) nsgfs_read
         nsgfs_read = nsgfs_read - nsgfs_local
         is_restarted = (COUNT(nsgfs_read /= 0) == 0)
         IF (.NOT. is_restarted) &
            CALL close_file(unit_number=unit_nr)
      END IF
      DEALLOCATE (nsgfs_local, nsgfs_read)

      CALL mp_bcast(is_restarted, para_env_global%source, para_env_global%group)

      IF (.NOT. is_restarted) THEN
         CALL cp_warn(__LOCATION__, &
                      "Initial electron density of the scattering region cannot be read from the file '"//TRIM(filename)// &
                      "' which is the NEGF restart file for a different system.")
         CALL timestop(handle)
         RETURN
      END IF

      ! read matrix elements
      maxcols = max_bcast_items_dp/nrows
      IF (maxcols < 1) maxcols = 1
      IF (maxcols > ncols) maxcols = ncols

      ALLOCATE (r2d(nrows, maxcols))
      DO ispin = 1, nspins
         ncols_read = 0

         DO WHILE (ncols_read < ncols)
            ncols_batch = ncols - ncols_read
            IF (ncols_batch > maxcols) ncols_batch = maxcols

            IF (para_env_global%ionode) THEN
               DO i = 1, ncols_batch
                  READ (unit_nr) r2d(:, i)
               END DO
            END IF

            CALL mp_bcast(r2d(:, 1:ncols_batch), para_env_global%source, para_env_global%group)

            CALL cp_fm_set_submatrix(h_s(ispin)%matrix, r2d, &
                                     start_row=1, start_col=ncols_read + 1, n_rows=nrows, n_cols=ncols_batch)

            ncols_read = ncols_read + ncols_batch
         END DO
      END DO
      DEALLOCATE (r2d)

      IF (para_env_global%ionode) &
         CALL close_file(unit_number=unit_nr)

      IF (log_unit > 0) THEN
         WRITE (log_unit, '(/,T2,A)') "NEGF| Initial electron density has been read from "//TRIM(filename)
      END IF

      CALL timestop(handle)
   END FUNCTION negf_read_restart

! **************************************************************************************************
!> \brief Write the Kohn-Sham matrix of the scattering region into a NEGF restart file.
!>        It may be a good idea to save the electron density matrix instead.
!> \param h_s                Kohn-Sham matrix to be saved
!> \param atomlist_s         atoms belonging to the scattering region
!> \param subsys             QuickStep subsystem of the device force environment
!> \param logger             I/O logger
!> \param restart_print_key  NEGF/PRINT/RESTART input section
! **************************************************************************************************
   SUBROUTINE negf_write_restart(h_s, atomlist_s, subsys, logger, restart_print_key)
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(in)       :: h_s
      INTEGER, DIMENSION(:), INTENT(in)                  :: atomlist_s
      TYPE(qs_subsys_type), POINTER                      :: subsys
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(section_vals_type), POINTER                   :: restart_print_key

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_write_restart', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ispin, maxcols, natoms, &
                                                            ncols, ncols_batch, ncols_written, &
                                                            nparticles, nrows, nspins, unit_nr
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: nsgfs, nsgfs_local
      REAL(kind=dp), ALLOCATABLE, DIMENSION(:, :)        :: r2d
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set

      CALL timeset(routineN, handle)

      IF (BTEST(cp_print_key_should_output(logger%iter_info, restart_print_key), cp_p_file)) THEN
         nspins = SIZE(h_s)
         natoms = SIZE(atomlist_s)
         ALLOCATE (nsgfs_local(natoms))

         CALL qs_subsys_get(subsys, particle_set=particle_set, qs_kind_set=qs_kind_set)
         nparticles = SIZE(particle_set)

         ALLOCATE (nsgfs(nparticles))
         CALL get_particle_set(particle_set, qs_kind_set, nsgf=nsgfs)
         DO i = 1, natoms
            nsgfs_local(i) = nsgfs(atomlist_s(i))
         END DO
         DEALLOCATE (nsgfs)

         CALL cp_fm_get_info(h_s(1)%matrix, nrow_global=nrows, ncol_global=ncols)

         IF (debug_this_module) THEN
            CPASSERT(nrows == ncols)
            CPASSERT(SUM(nsgfs_local(:)) == nrows)
         END IF

         unit_nr = cp_print_key_unit_nr(logger, restart_print_key, &
                                        extension=".negfpmat", file_status="REPLACE", file_action="WRITE", &
                                        do_backup=.TRUE., file_form="UNFORMATTED")
         IF (unit_nr > 0) THEN
            WRITE (unit_nr) restart_magic
            WRITE (unit_nr) nspins, nrows, natoms
            WRITE (unit_nr) nsgfs_local(1:natoms)
         END IF
         DEALLOCATE (nsgfs_local)

         maxcols = max_bcast_items_dp/nrows
         IF (maxcols < 1) maxcols = 1
         IF (maxcols > ncols) maxcols = ncols

         ALLOCATE (r2d(nrows, maxcols))
         DO ispin = 1, nspins
            ncols_written = 0

            DO WHILE (ncols_written < ncols)
               ncols_batch = ncols - ncols_written
               IF (ncols_batch > maxcols) ncols_batch = maxcols

               CALL cp_fm_get_submatrix(h_s(ispin)%matrix, r2d, &
                                        start_row=1, start_col=ncols_written + 1, &
                                        n_rows=nrows, n_cols=ncols_batch)

               IF (unit_nr > 0) THEN
                  DO i = 1, ncols_batch
                     WRITE (unit_nr) r2d(:, i)
                  END DO
               END IF

               ncols_written = ncols_written + ncols_batch
            END DO
         END DO

         DEALLOCATE (r2d)
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_write_restart

! **************************************************************************************************
!> \brief Print electron density cube file.
!> \param e_density_print_key    input section
!> \param qs_env                 Quickstep environment
!>                               (accessed components: pw_env, rho, subsys)
!> \param log_unit               output unit
! **************************************************************************************************
   SUBROUTINE negf_print_density_cube(e_density_print_key, qs_env, log_unit)
      TYPE(section_vals_type), POINTER                   :: e_density_print_key
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(in)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_print_density_cube', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=6)                                   :: my_pos_cube
      CHARACTER(len=default_path_length)                 :: filename, mpi_filename
      INTEGER                                            :: handle, nspins, unit_nr
      LOGICAL                                            :: append_cube, mpi_io
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rho_elec_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: rho_r
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pw_pools
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_rho_type), POINTER                         :: rho_struct
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      IF (BTEST(cp_print_key_should_output(logger%iter_info, e_density_print_key), cp_p_file)) THEN
         append_cube = section_get_lval(e_density_print_key, "APPEND")
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         ELSE
            my_pos_cube = "REWIND"
         END IF

         CALL get_qs_env(qs_env, rho=rho_struct, pw_env=pw_env, subsys=subsys)
         CALL qs_rho_get(rho_struct, rho_r=rho_r)
         CALL qs_subsys_get(subsys, particles=particles)
         nspins = SIZE(rho_r)

         IF (nspins > 1) THEN
            CALL pw_env_get(pw_env=pw_env, auxbas_pw_pool=auxbas_pw_pool, pw_pools=pw_pools)
            CALL pw_pool_create_pw(pool=auxbas_pw_pool, pw=rho_elec_rspace%pw, use_data=REALDATA3D, in_space=REALSPACE)
            CALL pw_copy(rho_r(1)%pw, rho_elec_rspace%pw)
            CALL pw_axpy(rho_r(2)%pw, rho_elec_rspace%pw)

            filename = "ELECTRON_DENSITY"
            mpi_io = .TRUE.
            unit_nr = cp_print_key_unit_nr(logger, e_density_print_key, &
                                           extension=".cube", middle_name=TRIM(filename), &
                                           file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                           fout=mpi_filename)
            IF (log_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=unit_nr, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF

               WRITE (log_unit, "(/,T2,A,/,/,T2,A)") &
                  "The sum of alpha and beta density is written in cube file format to the file:", &
                  TRIM(filename)
            END IF
            CALL cp_pw_to_cube(rho_elec_rspace%pw, unit_nr, "SUM OF ALPHA AND BETA DENSITY", &
                               particles=particles, stride=section_get_ivals(e_density_print_key, "STRIDE"), &
                               mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, e_density_print_key, mpi_io=mpi_io)
            CALL pw_copy(rho_r(1)%pw, rho_elec_rspace%pw)
            CALL pw_axpy(rho_r(2)%pw, rho_elec_rspace%pw, alpha=-1.0_dp)
            filename = "SPIN_DENSITY"
            ! mpi_io = .TRUE.
            mpi_io = .FALSE.
            unit_nr = cp_print_key_unit_nr(logger, e_density_print_key, &
                                           extension=".cube", middle_name=TRIM(filename), &
                                           file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                           fout=mpi_filename)
            IF (log_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=unit_nr, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (log_unit, "(/,T2,A,/,/,T2,A)") &
                  "The spin density is written in cube file format to the file:", &
                  TRIM(filename)
            END IF
            CALL cp_pw_to_cube(rho_elec_rspace%pw, unit_nr, "SPIN DENSITY", particles=particles, &
                               stride=section_get_ivals(e_density_print_key, "STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, e_density_print_key, mpi_io=mpi_io)
            CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_elec_rspace%pw)
         ELSE
            filename = "ELECTRON_DENSITY"
            !mpi_io = .TRUE.
            mpi_io = .FALSE.
            unit_nr = cp_print_key_unit_nr(logger, e_density_print_key, &
                                           extension=".cube", middle_name=TRIM(filename), &
                                           file_position=my_pos_cube, log_filename=.FALSE., mpi_io=mpi_io, &
                                           fout=mpi_filename)
            IF (log_unit > 0) THEN
               IF (.NOT. mpi_io) THEN
                  INQUIRE (UNIT=unit_nr, NAME=filename)
               ELSE
                  filename = mpi_filename
               END IF
               WRITE (log_unit, "(/,T2,A,/,/,T2,A)") &
                  "The electron density is written in cube file format to the file:", &
                  TRIM(filename)
            END IF
            CALL cp_pw_to_cube(rho_r(1)%pw, unit_nr, "ELECTRON DENSITY", particles=particles, &
                               stride=section_get_ivals(e_density_print_key, "STRIDE"), mpi_io=mpi_io)
            CALL cp_print_key_finished_output(unit_nr, logger, e_density_print_key, mpi_io=mpi_io)
         END IF
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_print_density_cube

! **************************************************************************************************
!> \brief Print Hartree potential cube file.
!> \param v_hartree_print_key    input section
!> \param qs_env                 Quickstep environment
!>                               (accessed components: pw_env, subsys, v_hartree_rspace)
!> \param log_unit               output unit
! **************************************************************************************************
   SUBROUTINE negf_print_vhartree_cube(v_hartree_print_key, qs_env, log_unit)
      TYPE(section_vals_type), POINTER                   :: v_hartree_print_key
      TYPE(qs_environment_type), POINTER                 :: qs_env
      INTEGER, INTENT(in)                                :: log_unit

      CHARACTER(LEN=*), PARAMETER :: routineN = 'negf_print_vhartree_cube', &
         routineP = moduleN//':'//routineN

      CHARACTER(len=6)                                   :: my_pos_cube
      CHARACTER(len=default_path_length)                 :: filename, mpi_filename
      INTEGER                                            :: handle, unit_nr
      LOGICAL                                            :: append_cube, mpi_io
      REAL(kind=dp)                                      :: udvol
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(particle_list_type), POINTER                  :: particles
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: aux_r
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(pw_type), POINTER                             :: v_hartree_rspace
      TYPE(qs_subsys_type), POINTER                      :: subsys

      CALL timeset(routineN, handle)
      logger => cp_get_default_logger()

      IF (BTEST(cp_print_key_should_output(logger%iter_info, v_hartree_print_key), cp_p_file)) THEN
         append_cube = section_get_lval(v_hartree_print_key, "APPEND")
         IF (append_cube) THEN
            my_pos_cube = "APPEND"
         ELSE
            my_pos_cube = "REWIND"
         END IF

         CALL get_qs_env(qs_env=qs_env, pw_env=pw_env, v_hartree_rspace=v_hartree_rspace, subsys=subsys)
         CALL qs_subsys_get(subsys, particles=particles)
         CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool)
         CALL pw_pool_create_pw(auxbas_pw_pool, aux_r%pw, use_data=REALDATA3D, in_space=REALSPACE)

         mpi_io = .FALSE.
         unit_nr = cp_print_key_unit_nr(logger, v_hartree_print_key, &
                                        extension=".cube", middle_name="negf_v_hartree", &
                                        file_position=my_pos_cube, log_filename=.FALSE., &
                                        mpi_io=mpi_io, fout=mpi_filename)

         IF (log_unit > 0) THEN
            IF (.NOT. mpi_io) THEN
               INQUIRE (UNIT=unit_nr, NAME=filename)
            ELSE
               filename = mpi_filename
            END IF
            WRITE (log_unit, "(/,T2,A,/,/,T2,A)") &
               "The Hartree potential is written in cube file format to the file:", &
               TRIM(filename)
         END IF

         CALL pw_copy(v_hartree_rspace, aux_r%pw)
         udvol = 1.0_dp/v_hartree_rspace%pw_grid%dvol
         CALL pw_scale(aux_r%pw, udvol)

         CALL cp_pw_to_cube(aux_r%pw, unit_nr, "HARTREE POTENTIAL", particles=particles, &
                            stride=section_get_ivals(v_hartree_print_key, "STRIDE"), mpi_io=mpi_io)
         CALL cp_print_key_finished_output(unit_nr, logger, v_hartree_print_key, mpi_io=mpi_io)

         CALL pw_pool_give_back_pw(auxbas_pw_pool, aux_r%pw)
      END IF

      CALL timestop(handle)
   END SUBROUTINE negf_print_vhartree_cube
END MODULE negf_io
