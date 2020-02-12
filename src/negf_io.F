!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief NEGF based quantum transport calculations
! **************************************************************************************************

MODULE negf_io
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_type
   USE cp_output_handling,              ONLY: cp_p_file,&
                                              cp_print_key_finished_output,&
                                              cp_print_key_should_output,&
                                              cp_print_key_unit_nr
   USE cp_realspace_grid_cube,          ONLY: cp_pw_to_cube
   USE input_section_types,             ONLY: section_get_ivals,&
                                              section_get_lval,&
                                              section_vals_type
   USE kinds,                           ONLY: default_path_length,&
                                              dp
   USE particle_list_types,             ONLY: particle_list_type
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
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_subsys_types,                 ONLY: qs_subsys_get,&
                                              qs_subsys_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'negf_io'
   LOGICAL, PARAMETER, PRIVATE          :: debug_this_module = .TRUE.

   PUBLIC :: negf_print_density_cube, negf_print_vhartree_cube

CONTAINS

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
