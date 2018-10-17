!********************************************************************************
!   Cassandra - An open source atomistic Monte Carlo software package
!   developed at the University of Notre Dame.
!   http://cassandra.nd.edu
!   Prof. Edward Maginn <ed@nd.edu>
!   Copyright (2013) University of Notre Dame du Lac
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
!*******************************************************************************

SUBROUTINE GCMC_Driver
  !****************************************************************************
  !
  ! The subroutine performs GEMC moves. 
  ! 
  ! Called by
  !
  !   main.f90
  !
  ! Revision history
  !
  !   12/10/13 : Beta Release
!*******************************************************************************

  USE Global_Variables
  USE Random_Generators
  USE File_Names
  USE Energy_Routines
  USE Read_Write_Checkpoint

  IMPLICIT NONE

!  !$ include 'omp_lib.h'

  INTEGER :: i, ibox, is 

  REAL(DP) :: rand_no
  REAL(DP) :: time_start, now_time, thermo_time, coord_time, block_avg_time

  LOGICAL :: overlap, complete, write_flag

  ! The total number of trial move array may not have been set if this
  ! is a fresh run i.e. start_type == make_config. Otherwise this array
  ! is set in read_checkpoint subroutine in the module Read_Write_Checkpoint
  IF(.NOT. ALLOCATED(ntrials)) ALLOCATE(ntrials(nspecies,nbr_boxes))
  IF(.NOT. ALLOCATED(tot_trials)) ALLOCATE(tot_trials(nbr_boxes))

  complete = .FALSE.
  openmp_flag = .FALSE.
  thermo_time = 0.0
  coord_time = 0.0
  block_avg_time = 0.0

!$ openmp_flag = .TRUE.

  WRITE(*,*) 'openmp_flag = ', openmp_flag

  IF(.NOT. openmp_flag) THEN
     CALL cpu_time(time_start)
  ELSE
!$  time_start = omp_get_wtime()
  END IF

  i_mcstep = initial_mcstep

  DO WHILE (.NOT. complete)

     i_mcstep = i_mcstep + 1

     ! Change mode from equilibration to production if specified
     IF (change_to_production .AND. i_mcstep > nsteps_until_prod) THEN
       int_run_type = run_prod
     END IF

     !*****************************************************************************
     ! select a move from Golden Sampling scheme
     !*****************************************************************************
  
     rand_no = rranf()
 
     IF (rand_no <= cut_trans) THEN

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF

        CALL Translate

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_trans) = movetime(imove_trans) + time_e - time_s

     ELSE IF ( rand_no <= cut_rot) THEN
 
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF

        CALL Rotate

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_rot) = movetime(imove_rot) + time_e - time_s

     ELSE IF (rand_no <= cut_torsion) THEN
 
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF

        CALL Rotate_Dihedral

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_dihedral) = movetime(imove_dihedral) + time_e - time_s

     ELSE IF (rand_no <= cut_angle) THEN
 
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF

        CALL Angle_Distortion

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_angle) = movetime(imove_angle) + time_e - time_s

    ELSE IF (rand_no <= cut_insertion) THEN
 
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF
      
        CALL Insertion

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_insert) = movetime(imove_insert) + time_e - time_s

     ELSE IF (rand_no <= cut_deletion) THEN
 
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF

        CALL Deletion

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_delete) = movetime(imove_delete) + time_e - time_s

     ELSE IF ( rand_no <= cut_regrowth) THEN
        
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF

        CALL Cut_N_Grow

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_regrowth) = movetime(imove_regrowth) + time_e - time_s

     ELSE IF (rand_no <= cut_atom_displacement) THEN

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF

        CALL Atom_Displacement

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_atom_displacement) = movetime(imove_atom_displacement) + time_e - time_s

     END IF

     IF(.NOT. openmp_flag) THEN
        CALL cpu_time(now_time)
     ELSE
!$      now_time = omp_get_wtime()
     END IF

     now_time = ((now_time - time_start) / 60.0_DP) 
     IF(.NOT. timed_run) THEN
        IF(i_mcstep == n_mcsteps) complete = .TRUE.
     ELSE
        IF(now_time .GT. n_mcsteps) complete = .TRUE.
     END IF

     !*****************************************************************************
     ! check if compute properties this step
     !*****************************************************************************
     IF (echeck) THEN
        IF(MOD(i_mcstep,echeck_freq) == 0) THEN
           DO ibox = 1,nbr_boxes
              CALL Check_System_Energy(ibox)
           END DO
        END IF
     END IF

     write_flag = .FALSE.
     IF(.NOT. timed_run) THEN
        IF (MOD(i_mcstep,nthermo_freq) == 0) write_flag = .TRUE.
     ELSE
        now_time = now_time - thermo_time
        IF(now_time .GT. nthermo_freq) THEN
           thermo_time = thermo_time + nthermo_freq
           write_flag = .TRUE.
        END IF
     END IF

     ! Write the information to various files at regular intervals
     IF (write_flag) THEN
        IF (.NOT. block_avg ) THEN
           ! write instantaneous properties
           DO ibox = 1, nbr_boxes
              CALL Write_Properties(ibox)
           END DO
        ELSE
           ! block averages

           ! Accumulate averages
           DO ibox = 1, nbr_boxes
              CALL Accumulate(ibox)
           END DO
              
           ! Check if write block avgs this step
           write_flag = .FALSE.
           IF(.NOT. timed_run) THEN
              IF (MOD(i_mcstep,block_avg_freq) == 0) write_flag = .TRUE.
           ELSE
              now_time = now_time - block_avg_time
              IF(now_time .GT. block_avg_freq) THEN
                 block_avg_time = block_avg_time + block_avg_freq
                 write_flag = .TRUE.
              END IF
           END IF

           ! Write block avgs
           DO ibox = 1, nbr_boxes
              IF(write_flag) THEN
                 CALL Write_Properties(ibox)
              END IF
           END DO
              
        END IF
     END IF

     !*****************************************************************************
     ! Check if write coords this step
     !*****************************************************************************
     write_flag = .FALSE.
     IF (.NOT. timed_run) THEN
        IF ( MOD(i_mcstep,ncoord_freq) == 0) write_flag = .TRUE.
     ELSE
        now_time = now_time - coord_time
        IF(now_time .GT. ncoord_freq) THEN
           coord_time = coord_time + nthermo_freq
           write_flag = .TRUE.
        END IF
     END IF

     IF (write_flag) THEN
        CALL Write_Checkpoint
        DO ibox = 1, nbr_boxes
           CALL Write_Coords(ibox)
        END DO
     END IF

  END DO

  CLOSE(50)


  END SUBROUTINE GCMC_Driver
