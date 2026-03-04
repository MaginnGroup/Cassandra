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
!********************************************************************************

!********************************************************************************
SUBROUTINE Pregen_Driver
  !******************************************************************************
  ! The subroutine calls load_next_frame every step and checks for end of
  ! simulation and checks when to output properties
  !
  ! CALLED BY
  !
  !        main
  !
  ! CALLS
  !
  !
  !*******************************************************************************

  USE Global_Variables
  USE Random_Generators
  USE File_Names
  USE Energy_Routines
  USE Read_Write_Checkpoint
  USE Simulation_Properties
  USE Trajectory_Reader_Routines
  USE XTC_Routines, ONLY : Close_XTC

  IMPLICIT NONE

!  !$ include 'omp_lib.h'

  INTEGER :: i,j,k, ibox, is, ifrag, ireac
  INTEGER :: ii,jj

  REAL(DP) :: rand_no
  REAL(DP) :: time_start, now_time, thermo_time, coord_time, block_avg_time

  LOGICAL :: overlap
  LOGICAL :: write_flag, complete

  TYPE(Energy_Class) :: energy_old

  ! The total number of trial move array may not have been set if this
  ! is a fresh run i.e. start_type == make_config. Otherwise this array
  ! is set in read_checkpoint subroutine in the module Read_Write_Checkpoint
  IF(.NOT. ALLOCATED(ntrials)) ALLOCATE(ntrials(nspecies,nbr_boxes))
  IF(.NOT. ALLOCATED(tot_trials)) ALLOCATE(tot_trials(nbr_boxes))

  thermo_time = 0.0
  coord_time = 0.0
  block_avg_time = 0.0
  openmp_flag = .FALSE.
  write_flag = .FALSE.
  complete = .FALSE.
  early_end = .FALSE.
  chpot(:,:) = 0.0_DP
  chpotid(:,:) = 0.0_DP
  nvolumes(:) = 0
  nvol_success(:) = 0
  ivol_success(:) = 0

  i_mcstep = initial_mcstep

  DO ibox = 1,nbr_boxes
     nsuccess(:,ibox)%displacement = 0
     nsuccess(:,ibox)%rotation = 0
     nsuccess(:,ibox)%displacement_e = 0
     nsuccess(:,ibox)%rotation_e = 0
     nsuccess(:,ibox)%dihedral = 0
     nsuccess(:,ibox)%angle = 0
     nsuccess(:,ibox)%insertion = 0
     nsuccess(:,ibox)%deletion = 0
     nsuccess(:,ibox)%switch = 0
     nsuccess(:,ibox)%disp_atom = 0
     ntrials(:,ibox)%displacement = 0
     ntrials(:,ibox)%rotation = 0
     ntrials(:,ibox)%dihedral = 0
     ntrials(:,ibox)%angle = 0
     ntrials(:,ibox)%insertion = 0
     ntrials(:,ibox)%deletion = 0
     ntrials(:,ibox)%switch = 0
     ntrials(:,ibox)%disp_atom = 0
     ntrials(:,ibox)%cpcalc = 0
     tot_trials(ibox) = 0
  END DO


!$ openmp_flag = .TRUE.

  WRITE(*,*) 'openmp_flag = ', openmp_flag

  IF(.NOT. openmp_flag) THEN
     CALL cpu_time(time_start)
  ELSE
!$  time_start = omp_get_wtime()
  END IF

  DO WHILE (.NOT. complete)

     !*****************************************************************************
     ! Load and advance to next frame
     !*****************************************************************************

     CALL Load_Next_Frame

     IF (early_end) THEN
             WRITE(logunit,*)
             WRITE(logunit,*) 'WARNING: Simulation ended early after ', i_mcstep, &
                     'frames due to reaching the end of the trajectory file'
             WRITE(logunit,*)
             EXIT
     END IF

     i_mcstep = i_mcstep + 1


     ! do widom insertions, if applicable to this simulation and step
     IF (widom_flag) THEN

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF

        CALL Widom_Subdriver

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_widom) = movetime(imove_widom) + time_e - time_s

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
        IF (need_energy .AND. (int_sim_type == sim_pregen)) THEN
           DO ibox = 1, nbr_boxes
              CALL Compute_System_Total_Energy(ibox,.TRUE.,overlap)
              IF (overlap) THEN
                 err_msg = ''
                 err_msg(1) = 'Atomic overlap in the configuration'
                 CALL Clean_Abort(err_msg,'Check_System_Energy')
              END IF
           END DO
        END IF


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

     write_flag = .FALSE.
     IF(.NOT. timed_run) THEN
        IF (MOD(i_mcstep,ncoord_freq) == 0) write_flag = .TRUE.
     ELSE
        now_time = now_time - thermo_time
        IF(now_time .GT. ncoord_freq) THEN
           coord_time = coord_time + ncoord_freq
           write_flag = .TRUE.
        END IF
     END IF

     IF (write_flag) THEN
        CALL Write_Checkpoint
        IF (int_coord_style .EQ. 1) THEN
          DO ibox = 1, nbr_boxes
            CALL Write_Coords_XYZ(ibox)
          END DO
        ELSE IF (int_coord_style .EQ. 2) THEN
          CALL Write_Coords_Custom
        ELSE
          err_msg = ''
          err_msg(1) = 'Invalid coordinate style'
          CALL Clean_Abort(err_msg,'Pregen_Driver')
        END IF
     END IF

  END DO
  DO ibox = 1, nbr_boxes
        IF (xtc_is_open(ibox)) CALL Close_XTC(ibox)
  END DO


END SUBROUTINE Pregen_Driver
