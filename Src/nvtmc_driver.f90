
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

SUBROUTINE NVTMC_Driver
  !******************************************************************************
  ! The subroutine performs NVT MC moves. 
  ! 
  ! CALLED BY
  !
  !        main
  !
  ! CALLS
  !        CPU_Time
  !        translate
  !        rotate
  !        Rigid_Dihedral_Change
  !        Angle_Distortion
  !        Cut_N_Grow
  !        Atom_Displacement
  !        Chempot
  !        Accumulate
  !        Write_Checkpoint
  !        Write_Properties
  !        Reset
  !        Write_Coords
  !        Compute_Total_System_Energy
  !        Write_Trials_Success
  !
  !  08/07/13 : Created beta version

  !*******************************************************************************

  USE Run_Variables
  USE Random_Generators
  USE File_Names
  USE Energy_Routines
  USE Read_Write_Checkpoint

  IMPLICIT NONE

!  !$ include 'omp_lib.h'

  INTEGER :: i,j, this_box, ibox, ireac, which_step, is, im
  INTEGER :: alive1, alive2, im1
  INTEGER, ALLOCATABLE, DIMENSION(:) :: n_inside_old

  REAL(DP) :: rand_no
  REAL(DP) :: time_start, now_time, thermo_time, coord_time

  LOGICAL :: overlap, aok, write_flag, complete, inside_ch
  LOGICAL, DIMENSION(:), ALLOCATABLE :: next_write, next_rdf_write

  ! The total number of trial move array may not have been set if this
  ! is a fresh run i.e. start_type == make_config. Otherwise this array
  ! is set in read_checkpoint subroutine in the module Read_Write_Checkpoint
  IF(.NOT. ALLOCATED(ntrials)) ALLOCATE(ntrials(nspecies,nbr_boxes))
  IF(.NOT. ALLOCATED(tot_trials)) ALLOCATE(tot_trials(nbr_boxes))
  
  ALLOCATE(next_write(nbr_boxes))
  ALLOCATE(next_rdf_write(nbr_boxes))
  next_write(:) = .false.
  next_rdf_write(:) = .false.
  thermo_time = 0.0
  coord_time = 0.0
  write_flag = .FALSE.
  complete = .FALSE.
  openmp_flag = .FALSE.

  i = 0

!$ openmp_flag = .TRUE.

  write(*,*) 'openmp_flag = ', openmp_flag

  IF(.NOT. openmp_flag) THEN
     CALL cpu_time(time_start)
  ELSE
!$  time_start = omp_get_wtime()
  END IF

  DO WHILE (.NOT. complete)

     i = i + 1
     which_step = i


     ! We will select a move from Golden Sampling scheme
  
     rand_no = rranf()
    
     IF (rand_no <= cut_trans) THEN
 
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF

        CALL Translate(this_box,which_step)

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

     ELSE IF (rand_no <= cut_regrowth) THEN
 
        IF(.NOT. openmp_flag) THEN
           CALL cpu_time(time_s)
        ELSE
!$        time_s = omp_get_wtime()
        END IF
        
        CALL Cut_N_Grow(this_box,i)

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

     ! Accumulate averages

     IF(cpcollect) THEN
        DO is = 1, nspecies
           CALL Chempot(this_box, is)
        END DO
     END IF

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

     CALL Accumulate(this_box)
     next_write(this_box) = .true.
     next_rdf_write(this_box) = .true.

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

           CALL Write_Checkpoint(i)

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

           DO ibox = 1, nbr_boxes
              
              CALL Write_Coords(ibox)

           END DO

        END IF
        
        write_flag = .FALSE.

     ELSE
     
        DO ibox = 1, nbr_boxes
           
           IF (tot_trials(ibox) /= 0) THEN
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
              
           END IF
           
        END DO
     
     END IF

     DO is = 1,nspecies

        IF(species_list(is)%int_insert == int_igas) THEN
           IF(mod(i,n_igas_moves(is)) == 0) CALL Update_Reservoir(is)
        END IF

     END DO
        
  END DO

  ! let us check if at the end of the simulation, the energies are properly updated

  write(logunit,*) '*********** Ending simulation *****************'
  write(logunit,*)
  write(logunit,*)

  DO this_box = 1, nbr_boxes

    ! Display the components of the energy.

     WRITE(logunit,*) 
     WRITE(logunit,*) '*****************************************'
     WRITE(logunit,'(A36,2X,I2)') ' Ending energy components for box', this_box
     WRITE(logunit,*) ' Atomic units - Extensive'
     WRITE(logunit,*) '*****************************************'
     WRITE(logunit,*)

     write(logunit,'(A,T30,F20.3)') 'Total system energy is' , energy(this_box)%total
     write(logunit,'(A,T30,F20.3)') 'Intra molecular energy is', energy(this_box)%intra
     write(logunit,'(A,T30,F20.3)') 'Bond energy is', energy(this_box)%bond
     write(logunit,'(A,T30,F20.3)') 'Angle energy is', energy(this_box)%angle
     write(logunit,'(A,T30,F20.3)') 'Dihedral energy is', energy(this_box)%dihedral
     write(logunit,'(A,T30,F20.3)') 'Improper energy is', energy(this_box)%improper
     write(logunit,'(A,T30,F20.3)') 'Intra nonbond vdw is', energy(this_box)%intra_vdw
     write(logunit,'(A,T30,F20.3)') 'Intra nonbond elec is', energy(this_box)%intra_q
     write(logunit,'(A,T30,F20.3)') 'Inter molecule vdw is', energy(this_box)%inter_vdw
     write(logunit,'(A,T30,F20.3)') 'Long range correction is', energy(this_box)%lrc
     write(logunit,'(A,T30,F20.3)') 'Inter molecule q is', energy(this_box)%inter_q
     write(logunit,'(A,T30,F20.3)') 'Reciprocal ewald is', energy(this_box)%ewald_reciprocal
     write(logunit,'(A,T30,F20.3)') 'Self ewald is', energy(this_box)%ewald_self
     write(logunit,*) '**************************************************'

     write(logunit,*)

     CALL Compute_Total_System_Energy(this_box,.TRUE., overlap)
     write(logunit,*)
     WRITE(logunit,*) '*****************************************'
     write(logunit,'(A52,2X,I2)') 'Components of energy from total energy call for box', this_box
     WRITE(logunit,*) 'Atomic units-Extensive'
     WRITE(logunit,*) '*****************************************'
     WRITE(logunit,*)
     

     write(logunit,'(A,T30,F20.3)') 'Total system energy is' , energy(this_box)%total
     write(logunit,'(A,T30,F20.3)') 'Intra molecular energy is', energy(this_box)%intra
     write(logunit,'(A,T30,F20.3)') 'Bond energy is', energy(this_box)%bond
     write(logunit,'(A,T30,F20.3)') 'Angle energy is', energy(this_box)%angle
     write(logunit,'(A,T30,F20.3)') 'Dihedral energy is', energy(this_box)%dihedral
     write(logunit,'(A,T30,F20.3)') 'Improper energy is', energy(this_box)%improper
     write(logunit,'(A,T30,F20.3)') 'Intra nonbond vdw is', energy(this_box)%intra_vdw
     write(logunit,'(A,T30,F20.3)') 'Intra nonbond elec is', energy(this_box)%intra_q
     write(logunit,'(A,T30,F20.3)') 'Inter molecule vdw is', energy(this_box)%inter_vdw
     write(logunit,'(A,T30,F20.3)') 'Long range correction is', energy(this_box)%lrc
     write(logunit,'(A,T30,F20.3)') 'Inter molecule q is', energy(this_box)%inter_q
     write(logunit,'(A,T30,F20.3)') 'Reciprocal ewald is', energy(this_box)%ewald_reciprocal
     write(logunit,'(A,T30,F20.3)') 'Self ewald is', energy(this_box)%ewald_self
     write(logunit,*)

  END DO

  IF(int_run_style == run_test) THEN

     OPEN(75,FILE='compare.dat',POSITION="APPEND")
     WRITE(75,"(T20,A,A)") testname, 'test in the nvt ensemble'
     WRITE(75,"(A,F24.12)") 'The total system energy is:', energy(1)%total
     WRITE(75,*)
     CLOSE(75)
  
  END IF   

  CALL Write_Trials_Success
     
END SUBROUTINE NVTMC_Driver
