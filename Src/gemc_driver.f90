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

SUBROUTINE GEMC_Driver

  !***************************************************************************
  ! The subroutine performs GEMC Simulations
  ! 
  ! Called by
  !
  !   main.f90
  !
  ! Revision history
  !
  !   12/10/13 : Beta Release
  !****************************************************************************

  USE Global_Variables
  USE Random_Generators
  USE Read_Write_Checkpoint
  USE Energy_Routines, ONLY: Compute_System_Total_Energy
  USE File_Names
  USE IO_Utilities, ONLY: Int_To_String

  IMPLICIT NONE

!$ include 'omp_lib.h'

  INTEGER :: i,j, this_box, ibox, is, other_box, which_step

  REAL(DP) :: rand_no
  REAL(DP) :: time_start, now_time, thermo_time, coord_time

  LOGICAL :: volume_move, insertion_move, overlap
  LOGICAL :: write_flag, complete

  LOGICAL, DIMENSION(:), ALLOCATABLE :: next_write, next_rdf_write

  ALLOCATE(next_write(nbr_boxes))
  ALLOCATE(next_rdf_write(nbr_boxes))
  next_write(:) = .false.
  next_rdf_write(:) = .false.
  volume_move = .false.
  insertion_move = .false.
  thermo_time = 0.0
  coord_time = 0.0
  openmp_flag = .FALSE.
  write_flag = .FALSE.
  complete = .FALSE.
  i = 0

  ! The total number of trial move array may not have been set if this
  ! is a fresh run i.e. start_type == make_config. Otherwise this array
  ! is set in read_checkpoint subroutine in the module Read_Write_Checkpoint
  IF(.NOT. ALLOCATED(ntrials)) ALLOCATE(ntrials(nspecies,nbr_boxes))
  IF(.NOT. ALLOCATED(tot_trials)) ALLOCATE(tot_trials(nbr_boxes))

!$ openmp_flag = .TRUE.


  IF(.NOT. openmp_flag) THEN
     CALL cpu_time(time_start)
  ELSE
!$  time_start = omp_get_wtime()
  END IF

  DO WHILE (.NOT. complete)

     i = i + 1

     ! We will select a move from Golden Sampling scheme

     rand_no = rranf()
     which_step = i

     IF (rand_no <= cut_trans) THEN
 
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF
        
        CALL Translate(this_box)
        
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

        CALL Rotate(this_box)

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
        
        CALL Rigid_Dihedral_Change(this_box)

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_dihedral) = movetime(imove_dihedral) + time_e - time_s

     ELSE IF (rand_no <= cut_volume) THEN
 
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF

        IF(freev .GT. 1 .OR. int_sim_type == sim_gemc_npt) THEN
 
          CALL Volume_Change(this_box)

        ELSE

          CALL GEMC_NVT_Volume(this_box, other_box)
          volume_move = .TRUE.

        END IF

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_volume) = movetime(imove_volume) + time_e - time_s
        
     ELSE IF (rand_no <= cut_angle) THEN
 
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF
        
        CALL Angle_Distortion(this_box)

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_angle) = movetime(imove_angle) + time_e - time_s
                
     ELSE IF ( rand_no <= cut_swap) THEN
 
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF
        
        CALL GEMC_Particle_Transfer

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        insertion_move = .TRUE.

        movetime(imove_swap) = movetime(imove_swap) + time_e - time_s

     ELSE IF ( rand_no <= cut_regrowth) THEN
 
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF

        CALL Cut_N_Grow(this_box)

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

        CALL Atom_Displacement(this_box)

        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_e)
        ELSE
!$         time_e = omp_get_wtime()
        END IF

        movetime(imove_atom_displacement) = movetime(imove_atom_displacement) + time_e - time_s

     END IF

     ! Accumulate averages. Note that the averages must
     ! be updated for both the boxes for the moves that
     ! affect both the boxes. These moves are volume
     ! fluctuation and particle swap

     IF(.NOT. openmp_flag) THEN
        CALL cpu_time(now_time)
     ELSE
!$      now_time = omp_get_wtime()
     END IF

     now_time = ((now_time - time_start) / 60.0_DP) 
     IF(.NOT. timed_run) THEN
        IF(i == n_mcsteps) complete = .TRUE.
     ELSE
        IF(now_time .GT. n_mcsteps) complete = .TRUE.
     END IF

     ! Accumulate the average for all the boxes
     
     DO ibox = 1, nbr_boxes
        CALL Accumulate(ibox)
     END DO

     ! next_write is true for all the boxes
     next_write(:) = .TRUE.
     next_rdf_write(:) = .TRUE.
   
     IF ( volume_move .OR. insertion_move ) THEN
        ! now set all the flags to false. In practice only one
        ! flag needs to be turned to false but it is a simpler
        ! solution than figuring out which move was performed
        volume_move = .false.
        insertion_move = .false.
     END IF
     
     ! Write the information to various files at regular intervals
     
     ! We will check this for all the boxes

     IF ( .NOT. block_average ) THEN

        ! instantaneous values are to be printed

        IF(.NOT. timed_run) THEN
           IF ( MOD(i,nthermo_freq) == 0) write_flag = .TRUE.
        ELSE
           now_time = now_time - thermo_time
           IF(now_time .GT. nthermo_freq) THEN
              thermo_time = thermo_time + nthermo_freq
              write_flag = .TRUE.
           END IF
        END IF

        IF(write_flag) THEN


           DO ibox = 1, nbr_boxes
              
              CALL Write_Properties(i,ibox)
              CALL Reset(ibox)


               DO j = 1,nspecies

                 IF(cpcollect) WRITE(*,*) 'Chem Pot spec:',j,'=',chpot(j,ibox) / ntrials(j,ibox)%cpcalc
               END DO

           END DO

        END IF
 
        write_flag = .FALSE.

        IF(.NOT. timed_run) THEN
           IF ( MOD(i,ncoord_freq) == 0) write_flag = .TRUE.
        ELSE
           now_time = now_time - coord_time
           IF(now_time .GT. ncoord_freq) THEN
              coord_time = coord_time + nthermo_freq
              write_flag = .TRUE.
           END IF
        END IF

        IF ( write_flag ) THEN


           CALL Write_Checkpoint(i)
           DO ibox = 1, nbr_boxes
              
              CALL Write_Coords(ibox)

           END DO

        END IF
        
        write_flag = .FALSE.

     ELSE
     
        DO ibox = 1, nbr_boxes
           
           IF(.NOT. timed_run) THEN
              IF ( MOD(i,nthermo_freq) == 0) write_flag = .TRUE.
           ELSE
              now_time = now_time - thermo_time
              IF(now_time .GT. nthermo_freq) THEN
                 IF(ibox == 1) thermo_time = thermo_time + nthermo_freq
                 write_flag = .TRUE.
              END IF
           END IF
           
           IF(write_flag) THEN
              IF (next_write(ibox)) THEN
                 CALL Write_Properties(tot_trials(ibox),ibox)
                 CALL Reset(ibox)
                 next_write(ibox) = .false.
              END IF
              IF(ibox == 1) CALL Write_Checkpoint(i)
           END IF
           
           write_flag = .FALSE.
           
           IF(.NOT. timed_run) THEN
              IF ( MOD(i,ncoord_freq) == 0) write_flag = .TRUE.
           ELSE
              now_time = now_time - coord_time
              IF(now_time .GT. ncoord_freq) THEN
                 IF(ibox == 1) coord_time = coord_time + nthermo_freq
                 write_flag = .TRUE.
              END IF
           END IF
           
           IF (write_flag) THEN
              IF (next_rdf_write(ibox)) THEN
                 CALL Write_Coords(ibox)
                 next_rdf_write(ibox) = .false.
              END IF
           END IF
           
           write_flag = .FALSE.
           
        END DO
     
     END IF

  END DO
  

END SUBROUTINE GEMC_Driver
