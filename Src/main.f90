!*******************************************************************************
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

!*******************************************************************************
PROGRAM Main
!*******************************************************************************
  ! Driver routine for CASSANDRA code.
  !
  ! CALLED BY
  !
  !        None
  !
  ! CALLS
  !
  !        date_and_time
  !        cpu_time
  !        clean_abort
  !        getarg
  !        Get_Runtime
  !        Name_Files
  !        HOSTNM
  !        NVTMC_Control
  !        NPTMC_Control
  !        GCMC_Control
  !        GEMC_Control
  !        Fragment_Control
  !        MCF_Control
  !        Get_Start_Type
  !        Get_Run_Type
  !        Initialize
  !        Reset
  !        Make_Config
  !        Read_Config
  !        Read_Checkpoint
  !        Ewald_Reciprocal_Lattice_Vector_Setup
  !        init_seeds
  !        Compute_Beads
  !        Get_Interal_Coords
  !        Get_COM
  !        Compute_Max_COM_Distance
  !        Fold_Molecule
  !        Compute_System_Total_Energy
  !        Get_Molecules_Species
  !        Angle_Distortion
  !        Rigid_Dihedral_Change
  !        NVTMC_Driver
  !        NPTMC_Driver
  !        GCMC_Driver
  !        GEMC_Driver
  !        Fragment_Driver
  !        Ring_Fragment
  !        Write_Subroutine_Times
  !
  !  08/07/13 : Created beta version
  !
  !
!*******************************************************************************

  USE Global_Variables
  USE File_Names
  USE IO_Utilities
  USE Input_Routines
  USE Read_Write_Checkpoint
  USE Energy_Routines
  USE Simulation_Properties
  USE Fragment_Growth

  IMPLICIT NONE

!  !$ include 'omp_lib.h'

  INTEGER(4) :: count
  INTEGER :: i, j, is, im, ia, lm, ibox
  INTEGER :: t_num

  INTEGER :: nyears, nmonths, ndays, nhours, nmin, nsec, nms
  CHARACTER(120) :: version
  CHARACTER(120) filename1
  CHARACTER(80) :: name

  LOGICAL :: overlap, check_charge

  REAL(DP) :: q_box

  REAL(DP) :: month_time, day_time, hour_time, min_time, sec_time, ms_time

  INTEGER :: IARGC

!********************************************************************************
! Code name and version. Change as updates are made.
  version = 'Cassandra Version 1.2 20160804'
! Get starting time information (intrinsic function)
  CALL DATE_AND_TIME(date,time,zone,begin_values)
  CALL cpu_time(start_time)
 
!Get the input file name as a command line parameter
  count = IARGC()
  IF (count < 1) THEN 
     err_msg = ""
     err_msg(1) = "Please specify input file &
          &name as command line argument."
     CALL Clean_Abort(err_msg,'Main')
  END IF

  CALL GETARG(1,filename1)
  inputfile = filename1

! Open up the input file
  OPEN(UNIT=inputunit,FILE=inputfile,STATUS="OLD",IOSTAT=openstatus,ACTION="READ")
  IF (openstatus .NE. 0) THEN
     err_msg = ""
     err_msg(1) = "Unable to open input file."
     CALL Clean_Abort(err_msg,'read_inputfile')
  ENDIF

! Now read input file and get the run_name
  CALL Get_Run_Name

! Create log file and write out some initial information
  CALL Name_Files(run_name,'.log',logfile)

  OPEN(unit=logunit,file=logfile,IOSTAT=openstatus)
  IF (openstatus .NE. 0) THEN
     err_msg = ""
     err_msg(1) = 'Unable to open logfile'
     err_msg(2) = logfile
     CALL Clean_Abort(err_msg,'Read_Inputfile')
  ENDIF

  WRITE(logunit,'(A80)')'********************************************************************************'
  WRITE(logunit,'(A80)')'               ______                                __                        '
  WRITE(logunit,'(A80)')'              / ____/___ _______________ _____  ____/ /________ _              '
  WRITE(logunit,'(A80)')'             / /   / __ `/ ___/ ___/ __ `/ __ \/ __  / ___/ __ `/              '
  WRITE(logunit,'(A80)')'            / /___/ /_/ (__  |__  ) /_/ / / / / /_/ / /  / /_/ /               '
  WRITE(logunit,'(A80)')'            \____/\__,_/____/____/\__,_/_/ /_/\__,_/_/   \__,_/                '
  WRITE(logunit,'(A80)')'                                                                               '
  WRITE(logunit,'(A80)')'********************************************************************************'
! Create a checkpoint file to periodically write system information

  CALL Name_Files(run_name,'.chk',checkpointfile)
  
  OPEN(unit=chkptunit,file=checkpointfile,IOSTAT=openstatus)
  IF (openstatus .NE. 0) THEN
     err_msg = ""
     err_msg(1) = 'Unable to open checkpointfile'
     err_msg(2) = logfile
     CALL Clean_Abort(err_msg,'Read_Checkpointfile')
  ENDIF

  WRITE(logunit,*)
  WRITE(logunit,'(A)') 'Run info'
  WRITE(logunit,'(A80)') '********************************************************************************'
  WRITE(logunit,'(a,a)') 'version: ', TRIM(version)
  WRITE(logunit,'(a,a)') 'inputfile: ', TRIM(inputfile)
  WRITE(logunit,'(a,a)') 'name of run: ',TRIM(run_name)

  CALL DATE_AND_TIME(date,time)

  WRITE(logunit,'(a18,x,a2,a1,a2,a3,a2)') 'date (mm/dd/yyyy):', date(5:6), '/', date(7:8), '/20', date(3:4)
  WRITE(logunit,'(a5,1x,a2,a1,a2,a1,a2)') 'time:',time(1:2),':',time(3:4),':',time(5:6)

! Get name of computer running code using intrinsic function.

  CALL HOSTNM(name)
  WRITE(logunit,'(a,a)') 'machine: ', TRIM(name)
  WRITE(logunit,'(A80)') '********************************************************************************'

  ! Standard level of output to logfile, or verbose output
  CALL Get_Verbosity_Info  

! Determine the simulation type, and then read in all the necessary information 
! from the input file for starting up that type of simulation
  CALL Get_Sim_Type

  ! Check what kind of simulation this is, and then call the appropriate routine
  ! that will load all the relevant information in from the input file
  IF (int_sim_type == sim_nvt .OR. int_sim_type == sim_nvt_min) THEN
     CALL NVTMC_Control
  ELSE IF (int_sim_type == sim_npt) THEN
     CALL NPTMC_Control
  ELSE IF (int_sim_type == sim_gcmc) THEN
     CALL GCMC_Control
  ELSE IF (int_sim_type == sim_gemc .OR. int_sim_type == sim_gemc_npt) THEN
     CALL GEMC_Control
  ELSE IF (int_sim_type == sim_frag .OR. int_sim_type == sim_ring) THEN
     CALL Fragment_Control
  ELSE IF (int_sim_type == sim_mcf) THEN
     CALL MCF_Control
  ELSE
     err_msg = ""
     err_msg(1) = 'Sim_Type unknown'
     err_msg(2) = logfile
     CALL Clean_Abort(err_msg,'Main')
  ENDIF
 
  ! Determine if it is equilibration or production or test
  CALL Get_Run_Type

  ! initialize counters
  CALL Initialize

  ! initialize random number generator
  CALL Init_Seeds(iseed1, iseed3)

  ! Initialize the atom list and molecule list arrays
  molecule_list(:,:)%live = .FALSE.
  molecule_list(:,:)%which_box = 0
  atom_list(:,:,:)%exist = .FALSE.
  molecule_list(:,:)%frac = 0.0_DP

  cbmc_flag = .FALSE.

  WRITE(*,*) 'Begin Cassandra simulation'
  WRITE(*,*) 

  WRITE(logunit,*)
  WRITE(logunit,'(A80)') '********************************************************************************'
  WRITE(logunit,'(A80)') '************************ Begin Cassandra simulation ****************************'
  WRITE(logunit,'(A80)') '********************************************************************************'

  WRITE(logunit,*)
  WRITE(logunit,'(A)') 'Initial configuration'
  WRITE(logunit,'(A80)') '********************************************************************************'

  DO ibox = 1, nbr_boxes
    IF (start_type(ibox) == 'make_config') THEN
       ! Insert molecules using CBMC
       CALL Make_Config(ibox)
       initial_mcstep = 0

    ELSEIF (start_type(ibox) == 'read_config') THEN
       ! Read in coordinates
       CALL Read_Config(ibox)
       initial_mcstep = 0

    ELSEIF (start_type(ibox) == 'add_to_config') THEN
       ! Add molecules using CBMC to configuration read from file
       CALL Read_Config(ibox)
       CALL Make_Config(ibox)
       initial_mcstep = 0
       
    ELSEIF (start_type(ibox) == 'checkpoint') THEN
       ! Restart from checkpoint. Shall we verify match between stuff in cpt 
       ! and stuff input, or override with cpt info?
       CALL Read_Checkpoint
       ! Read in info for all boxes, so exit loop
       EXIT
       
    ENDIF
  END DO
  WRITE(logunit,'(A80)') '********************************************************************************'

  ! Add LOCATE for any unplaced molecules in array under BOX==0, 
  ! in reverse order
  DO is = 1, nspecies
    DO im = max_molecules(is), SUM(nmols(is,1:nbr_boxes)) + 1, -1
      nmols(is,0) = nmols(is,0) + 1
      locate(nmols(is,0),is,0) = im
    END DO
  END DO

  ! compute the charge on each species and the box charge
  check_charge = .FALSE.
  DO ibox = 1, nbr_boxes
    IF (int_charge_style(ibox) /= charge_none) check_charge = .TRUE.
  END DO
  IF (check_charge) THEN
     WRITE(logunit,*)
     WRITE(logunit,'(A)') 'Charge neutrality check'
     WRITE(logunit,'(A80)') '********************************************************************************'

     DO is = 1, nspecies
        WRITE(logunit,'(X,A,T35,4x,f12.8)')'Species ' // TRIM(Int_To_String(is)) // ' has charge:', &
           species_list(is)%total_charge
     END DO
     WRITE(logunit,*)

     DO ibox = 1, nbr_boxes
        q_box = 0.0_DP
        DO is = 1, nspecies
           q_box = q_box + REAL(nmols(is,ibox),DP) * species_list(is)%total_charge
        END DO
        WRITE(logunit,'(X,A,T35,4X,f12.8)')'Box     ' // TRIM(Int_To_String(ibox)) // ' has charge:', q_box

        IF (ABS(q_box) > tiny_number .AND. &
           (int_charge_sum_style(ibox) == charge_ewald .OR. int_charge_sum_style(ibox) == charge_dsf)) THEN
           err_msg = ''
           err_msg(1) = 'Long-range electrostatics cannot be computed for box ' // TRIM(Int_To_String(ibox)) &
                     // ' due to its net charge'
           CALL Clean_Abort(err_msg,'main.f90')
        END IF
     END DO

     WRITE(logunit,'(A80)') '********************************************************************************'
  END IF

  ! Ewald stuff
  IF ( int_charge_sum_style(1) == charge_ewald) THEN
     
     ALLOCATE(hx(maxk,nbr_boxes),hy(maxk,nbr_boxes),hz(maxk,nbr_boxes), &
          hsq(maxk,nbr_boxes), Cn(maxk,nbr_boxes),Stat=AllocateStatus)
     ALLOCATE(nvecs(nbr_boxes))
     
     IF (AllocateStatus /=0) THEN
        err_msg = ""
        err_msg(1) = "Memory Could not be allocated for reciprocal vectors"
        CALL Clean_Abort(err_msg,'precalculate')
     END IF
     
     DO ibox = 1, nbr_boxes
        CALL  Ewald_Reciprocal_Lattice_Vector_Setup(ibox)
     END DO
     
     ! Here we can allocate the memory for cos_sum, sin_sum etc
     
     ALLOCATE(cos_sum(MAXVAL(nvecs),nbr_boxes))
     ALLOCATE(sin_sum(MAXVAL(nvecs),nbr_boxes))
     ALLOCATE(cos_sum_old(MAXVAL(nvecs),nbr_boxes))
     ALLOCATE(sin_sum_old(MAXVAL(nvecs),nbr_boxes))
     ALLOCATE(cos_sum_start(MAXVAL(nvecs),nbr_boxes))
     ALLOCATE(sin_sum_start(MAXVAL(nvecs),nbr_boxes))
     ALLOCATE(cos_mol(MAXVAL(nvecs), SUM(max_molecules)))
     ALLOCATE(sin_mol(MAXVAL(nvecs), SUM(max_molecules)))
     ! initialize these vectors
     cos_mol(:,:) = 0.0_DP
     sin_mol(:,:) = 0.0_DP
 
  END IF

  ! Compute total number of beads in each box

!  DO ibox = 1, nbr_boxes
!     CALL Compute_Beads(ibox)
!  END DO

  ! set up internal coordinates, the maximum distance from the COM and compute
  ! the total system energy to be displayed to the logfile.

  ! Internal coordinates
  CALL Get_Internal_Coords

  ! Calculate COM and distance of the atom farthest to the COM.
  DO ibox = 1, nbr_boxes
     DO is = 1, nspecies
        DO im = 1, nmols(is,ibox)
           lm = locate(im,is,ibox)
           CALL Get_COM(lm,is)
           CALL Compute_Max_Com_Distance(lm,is)
           CALL Fold_Molecule(lm,is,ibox)
        END DO
     END DO
  END DO

  ! Initialize accumulators
  CALL Init_Accumulators

  WRITE(logunit,*)
  WRITE(logunit,'(A)') 'Compute total energy'
  WRITE(logunit,'(A80)') '********************************************************************************'

  ! compute total system energy
  overlap = .FALSE.
     
  DO i = 1, nbr_boxes
       
     CALL Compute_System_Total_Energy(i,.TRUE.,overlap)

     IF (overlap) THEN
        ! overlap was detected between two atoms so abort the program
        err_msg = ''
        err_msg(1) = 'Overlap detected in the starting structure'
        err_msg(2) = 'Start type '//start_type(i)
        CALL Clean_Abort(err_msg,'Main')
     END IF
  END DO
  
  ! write out the initial energy components to the log file
  
  DO ibox = 1,nbr_boxes

     WRITE(logunit,'(X,A,X,I1)') 'Starting energy components for box', ibox
     WRITE(logunit,'(2X,A)') 'Atomic units-Extensive'
     WRITE(logunit,'(X,A59)') '-----------------------------------------------------------'
     WRITE(logunit,'(X,A,T30,F20.3)') 'Total system energy' , energy(ibox)%total
     WRITE(logunit,'(X,A,T30,F20.3)') 'Intra molecular energy', energy(ibox)%intra
     WRITE(logunit,'(3X,A,T30,F20.3)') 'Bond energy', energy(ibox)%bond
     WRITE(logunit,'(3X,A,T30,F20.3)') 'Bond angle energy', energy(ibox)%angle
     WRITE(logunit,'(3X,A,T30,F20.3)') 'Dihedral angle energy', energy(ibox)%dihedral
     WRITE(logunit,'(3X,A,T30,F20.3)') 'Improper angle energy', energy(ibox)%improper
     WRITE(logunit,'(X,A,T30,F20.3)') 'Intra molecule vdw', energy(ibox)%intra_vdw
     WRITE(logunit,'(X,A,T30,F20.3)') 'Intra molecule q', energy(ibox)%intra_q
     WRITE(logunit,'(X,A,T30,F20.3)') 'Inter molecule vdw', energy(ibox)%inter_vdw
     IF (int_vdw_sum_style(ibox) == vdw_cut_tail) &
        WRITE(logunit,'(X,A,T30,F20.3)') 'Long range correction', energy(ibox)%lrc
     WRITE(logunit,'(X,A,T30,F20.3)') 'Inter molecule q', energy(ibox)%inter_q
     IF (int_charge_sum_style(ibox) == charge_ewald) THEN
        WRITE(logunit,'(X,A,T30,I20)') 'Number of vectors', nvecs(ibox)
        WRITE(logunit,'(X,A,T30,F20.3)') 'Reciprocal ewald', energy(ibox)%ewald_reciprocal
        WRITE(logunit,'(X,A,T30,F20.3)') 'Self ewald', energy(ibox)%self
     ELSE IF (int_charge_sum_style(ibox) == charge_dsf) THEN
        WRITE(logunit,'(X,A,T30,F20.3)') 'Self DSF', energy(ibox)%self
     END IF

     WRITE(logunit,'(X,A59)') '-----------------------------------------------------------'
     WRITE(logunit,*)

  END DO
  WRITE(logunit,'(A80)') '********************************************************************************'

  ! Write initial properties, if needed
  DO ibox = 1, nbr_boxes
    IF (start_type(ibox) == 'make_config' .OR. n_mcsteps <= initial_mcstep) THEN
      IF (block_average) THEN
         ! need to write instantaneous props for this first step
         block_average = .FALSE.
         CALL Write_Properties(ibox)
         block_average = .TRUE.
         ! need to rewrite the header next call so additional props will be
         ! labeled as block_averages
         first_open(:,ibox) = .TRUE.
      ELSE
         CALL Write_Properties(ibox)
      END IF
    END IF
  END DO

  ! End program if no moves specified
  IF (n_mcsteps <= initial_mcstep) THEN
    WRITE(logunit,*)
    WRITE(logunit,'(A80)') '********************************************************************************'
    WRITE(logunit,'(A80)') '************************ Cassandra simulation complete *************************'
    WRITE(logunit,'(A80)') '********************************************************************************'
    WRITE(*,*)
    WRITE(*,*) 'Cassandra simulation complete'
    STOP
  END IF

  ! Now begin simulation by calling apropriate driver routines
  ! Add routines here...

  WRITE(logunit,*)
  WRITE(logunit,'(A)') 'Run simulation'
  WRITE(logunit,'(A80)') '********************************************************************************'
  WRITE(logunit,'(X,A9,X,A10,X,A5,X,A3,X,A3,X,A8,X,A9)') 'Step', 'Move', 'Mol', 'Spc', 'Box', 'Success', 'MaxWidth'

  IF (int_sim_type == sim_nvt .OR. int_sim_type == sim_nvt_min) THEN
     
     CALL NVTMC_Driver
     
  ELSE IF (int_sim_type == sim_npt) THEN
     
     CALL NPTMC_Driver

  ELSE IF (int_sim_type == sim_gcmc) THEN

     CALL GCMC_Driver

  ELSE IF (int_sim_type == sim_gemc .OR. &
           int_sim_type == sim_gemc_ig .OR. &
           int_sim_type == sim_gemc_npt) THEN
     
     CALL GEMC_Driver

  ELSE IF (int_sim_type == sim_frag) THEN

     CALL Fragment_Driver

  ELSE IF (int_sim_type == sim_ring) THEN

     CALL Ring_Fragment_Driver

  END IF
  WRITE(logunit,'(A80)') '********************************************************************************'

  ! write error of the mean
  IF (block_average) THEN
     DO ibox = 1, nbr_boxes
        CALL Write_Mean_Error(ibox)
     END DO
  END IF

  WRITE(logunit,*)
  WRITE(logunit,'(A)') 'Report ending energy & then re-compute from scratch'
  WRITE(logunit,'(A80)') '********************************************************************************'

  !***************************************************************************
  ! Report current energies and compute the energies from scratch

  DO ibox = 1, nbr_boxes
    WRITE(logunit,*)
    WRITE(logunit,*)

    ! Write the current components of the energy to log
    WRITE(logunit,'(X,A,X,I1)') 'Initial energy + deltas for box', ibox
    WRITE(logunit,'(2X,A)') 'Atomic units-Extensive'
    WRITE(logunit,'(X,A59)') '-----------------------------------------------------------'
    WRITE(logunit,'(X,A,T30,F20.3)') 'Total system energy' , energy(ibox)%total
    WRITE(logunit,'(X,A,T30,F20.3)') 'Intra molecular energy', energy(ibox)%intra
    WRITE(logunit,'(3X,A,T30,F20.3)') 'Bond energy', energy(ibox)%bond
    WRITE(logunit,'(3X,A,T30,F20.3)') 'Angle energy', energy(ibox)%angle
    WRITE(logunit,'(3X,A,T30,F20.3)') 'Dihedral energy', energy(ibox)%dihedral
    WRITE(logunit,'(3X,A,T30,F20.3)') 'Improper angle energy', energy(ibox)%improper
    WRITE(logunit,'(X,A,T30,F20.3)') 'Intra molecule vdw', energy(ibox)%intra_vdw
    WRITE(logunit,'(X,A,T30,F20.3)') 'Intra molecule q', energy(ibox)%intra_q
    WRITE(logunit,'(X,A,T30,F20.3)') 'Inter molecule vdw', energy(ibox)%inter_vdw
    IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN
       WRITE(logunit,'(X,A,T30,F20.3)') 'Long range correction', energy(ibox)%lrc
    END IF
    WRITE(logunit,'(X,A,T30,F20.3)') 'Inter molecule q', energy(ibox)%inter_q
    IF (int_charge_sum_style(ibox) == charge_ewald) THEN
       WRITE(logunit,'(X,A,T30,F20.3)') 'Reciprocal ewald', energy(ibox)%ewald_reciprocal
       WRITE(logunit,'(X,A,T30,F20.3)') 'Self ewald', energy(ibox)%self
    ELSE IF (int_charge_sum_style(ibox) == charge_dsf) THEN
       WRITE(logunit,'(X,A,T30,F20.3)') 'Self DSF', energy(ibox)%self
    END IF
    WRITE(logunit,'(X,A59)') '-----------------------------------------------------------'

    ! Write the current total energy to stdout
    WRITE(*,*)
    WRITE(*,'(X,A)') 'Energy of final configuration, box ' // TRIM(Int_To_String(ibox))
    WRITE(*,"(2X,A,T30,F24.12)") 'Initial energy + deltas = ', energy(ibox)%total

    ! Compute the energies from scratch
    CALL Compute_System_Total_Energy(ibox,.TRUE.,overlap)

    ! Write the recomputed total energy to stdout
    WRITE(*,"(2X,A,T30,F24.12)") 'Energy from scratch = ', energy(ibox)%total

    ! Write the recomputed energy components to log
    WRITE(logunit,*)
    WRITE(logunit,*)
    WRITE(logunit,'(X,A,X,I1)') 'Recomputed energy from scratch for box', ibox
    WRITE(logunit,'(2X,A)') 'Atomic units-Extensive'
    WRITE(logunit,'(X,A59)') '-----------------------------------------------------------'
    WRITE(logunit,'(X,A,T30,F20.3)') 'Total system energy' , energy(ibox)%total
    WRITE(logunit,'(X,A,T30,F20.3)') 'Intra molecular energy', energy(ibox)%intra
    WRITE(logunit,'(3X,A,T30,F20.3)') 'Bond energy', energy(ibox)%bond
    WRITE(logunit,'(3X,A,T30,F20.3)') 'Angle energy', energy(ibox)%angle
    WRITE(logunit,'(3X,A,T30,F20.3)') 'Dihedral energy', energy(ibox)%dihedral
    WRITE(logunit,'(3X,A,T30,F20.3)') 'Improper angle energy', energy(ibox)%improper
    WRITE(logunit,'(X,A,T30,F20.3)') 'Intra molecule vdw', energy(ibox)%intra_vdw
    WRITE(logunit,'(X,A,T30,F20.3)') 'Intra molecule q', energy(ibox)%intra_q
    WRITE(logunit,'(X,A,T30,F20.3)') 'Inter molecule vdw', energy(ibox)%inter_vdw
    IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN
      WRITE(logunit,'(X,A,T30,F20.3)') 'Long range correction', energy(ibox)%lrc
    END IF
    WRITE(logunit,'(X,A,T30,F20.3)') 'Inter molecule q', energy(ibox)%inter_q
    IF (int_charge_sum_style(ibox) == charge_ewald) THEN
       WRITE(logunit,'(X,A,T30,F20.3)') 'Reciprocal ewald', energy(ibox)%ewald_reciprocal
       WRITE(logunit,'(X,A,T30,F20.3)') 'Self ewald', energy(ibox)%self
    ELSE IF (int_charge_sum_style(ibox) == charge_dsf) THEN
       WRITE(logunit,'(X,A,T30,F20.3)') 'Self DSF', energy(ibox)%self
    END IF
    WRITE(logunit,'(X,A59)') '-----------------------------------------------------------'
  END DO
  WRITE(logunit,'(A80)') '********************************************************************************'

  ! Write success ratios of each move type
  CALL Write_Trials_Success

!$OMP PARALLEL

!  t_num = omp_get_thread_num()
!  write(logunit,*) 'Thread',t_num,"reporting for duty!"

!$OMP END PARALLEL


  IF(cpcollect) THEN
    DO is = 1,nspecies
      DO ibox = 1, nbr_boxes
        IF(ntrials(is,ibox)%cpcalc > 0) THEN
          WRITE(logunit,'(A,I2,A,I2,A,F24.12)') 'Chemical potential for species', &
            is,'in box',ibox,'is',chpot(is,ibox) / ntrials(is,i)%cpcalc
          WRITE(logunit,'(A,I2,A,I2,A,F24.12)') &
            'Ideal Chemical potential for species',is,'in box',i, 'is', &
            chpotid(is,i) / ntrials(is,i)%cpcalc
        END IF
      END DO
    END DO
  END IF


  CALL DATE_AND_TIME(date,time,zone,end_values)
  CALL cpu_time(tot_time)

  tot_time = tot_time - start_time

  CALL Write_Subroutine_Times

  WRITE(logunit,*)
  WRITE(logunit,'(A80)') '********************************************************************************'
  WRITE(logunit,'(A80)') '************************ Cassandra simulation complete *************************'
  WRITE(logunit,'(A80)') '********************************************************************************'
  WRITE(*,*)
  WRITE(*,*) 'Cassandra simulation complete'

  WRITE(logunit,*)
  WRITE(logunit,*) 'Program execution time'
  WRITE(logunit,'(A80)') '********************************************************************************'

  nyears = INT(tot_time/3600.0_DP/24.0_DP/365.0_DP)
  month_time = tot_time - REAL(nyears,DP) * 3600.0_DP * 24.0_DP * 365.0_DP
  nmonths = INT(month_time/3600_DP/24.0_DP/30.0_DP)
  day_time = month_time - REAL(nmonths,DP) * 3600_DP * 24.0_DP * 30.0_DP
  ndays = INT(day_time/3600.0_DP/24.0_DP)
  hour_time = day_time - REAL(ndays,DP) * 3600.0_DP * 24.0_DP
  nhours = INT(hour_time/3600.0_DP)
  min_time = hour_time - REAL(nhours,DP) * 3600.0_DP
  nmin = INT(min_time/60.0_DP)
  sec_time = min_time - REAL(nmin,DP) * 60.0_DP
  nsec = INT(sec_time)
  ms_time = sec_time - REAL(nsec,DP)
  nms = INT(ms_time*1000_DP)

  WRITE(logunit,'(I10,1x,a10)') nyears, 'Years'
  WRITE(logunit,'(I10,1x,a10)') nmonths, 'Months'
  WRITE(logunit,'(I10,1x,a10)') ndays,'Days'
  WRITE(logunit,'(I10,1x,a10)') nhours,'Hours'
  WRITE(logunit,'(I10,1x,a10)') nmin,'Minutes'
  WRITE(logunit,'(I10,1X,a10)') nsec,'Seconds'
  WRITE(logunit,'(I10,1x,a10)') nms,'ms'
  WRITE(logunit,'(A80)') '********************************************************************************'
  
END PROGRAM Main
