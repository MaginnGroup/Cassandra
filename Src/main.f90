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
PROGRAM Main
!********************************************************************************
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
  !        NVT_MC_Fragment_Control
  !        MCF_Control
  !        Get_Start_Type
  !        Get_Run_Type
  !        Initialize
  !        Reset
  !        Grow_Molecules
  !        Restart_From_Old
  !        Read_Checkpoint
  !        Ewald_Reciprocal_Lattice_Vector_Setup
  !        init_seeds
  !        Compute_Beads
  !        Get_Interal_Coords
  !        Get_COM
  !        Compute_Max_COM_Distance
  !        Fold_Molecule
  !        Compute_Total_System_Energy
  !        Get_Molecules_Species
  !        Update_Reservoirs
  !        Angle_Distortion
  !        Rigid_Dihedral_Change
  !        NVTMC_Driver
  !        NPTMC_Driver
  !        GCMC_Driver
  !        GEMC_Driver
  !        NVT_MC_Fragment_Driver
  !        NVT_MC_Ring_Fragment
  !        Write_Subroutine_Times
  !
  !  08/07/13 : Created beta version
  !
  !
!********************************************************************************

  USE Run_Variables
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
  INTEGER :: i, j, is, im, ia, this_im, ibox, nmol_is, this_box, int_phi
  INTEGER :: alive, t_num

  INTEGER :: nyears, nmonths, ndays, nhours, nmin, nsec, nms
  CHARACTER(120) :: version
  CHARACTER(120) filename1
  CHARACTER(80) :: name

  LOGICAL :: overlap, cbmc_overlap, superbad, inside_ch

  REAL(DP) :: attempt_prob, phi
  REAL(DP) :: E_st_vdw, E_st_qq, W_st_vdw, W_st_qq, e_lrc, w_lrc
  REAL(DP) :: q_tot_sys, q_mol

  REAL(DP) :: month_time, day_time, hour_time, min_time, sec_time, ms_time

  TYPE(Energy_Class) :: inrg, e_start

  INTEGER :: IARGC

  INTEGER, ALLOCATABLE, DIMENSION(:) :: frag_order
!********************************************************************************
! Code name and version. Change as updates are made.
  version = 'Cassandra 1.0'
  e_start%inter_vdw = 0.0_DP
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
  CALL Get_Runname

! Create log file and write out some initial information
  CALL Name_Files(run_name,'.log',logfile)

  OPEN(unit=logunit,file=logfile,IOSTAT=openstatus)
  IF (openstatus .NE. 0) THEN
     err_msg = ""
     err_msg(1) = 'Unable to open logfile'
     err_msg(2) = logfile
     CALL Clean_Abort(err_msg,'Read_Inputfile')
  ENDIF

! Create a checkpoint file to periodically write system information

  CALL Name_Files(run_name,'.chk',checkpointfile)
  
  OPEN(unit=chkptunit,file=checkpointfile,IOSTAT=openstatus)
  IF (openstatus .NE. 0) THEN
     err_msg = ""
     err_msg(1) = 'Unable to open checkpointfile'
     err_msg(2) = logfile
     CALL Clean_Abort(err_msg,'Read_Checkpointfile')
  ENDIF

  WRITE(logunit,'(a,/)') '*** Beginning Cassandra simulation ***'

  WRITE(logunit,'(a,a,/)') 'Version ',version
  WRITE(logunit,'(a,a,/)') 'inputfile ',inputfile

  i = LEN_TRIM(run_name)
  WRITE(logunit,'(a,a,/)') 'name of run: ',run_name(1:i)

  CALL DATE_AND_TIME(date,time)

  WRITE(logunit,'(a7,T10,a7,T20,a7)') 'Month','Day', 'Year'
  WRITE(logunit,'(a7,T10,a7,T20,a7)') date(5:6), date(7:8), '20'//date(3:4) 

  WRITE(logunit,*)
  WRITE(logunit,'(a7,1x,a2,a1,a2,a1,a2)') 'Time: ',time(1:2),':',time(3:4),':',time(5:6)
  WRITE(logunit,*)


! Get name of computer running code using intrinsic function.

  CALL HOSTNM(name)
  write(logunit,'(a,a,//)') 'machine: ',name(1:LEN_TRIM(name))

! Determine the simulation type, and then read in all the necessary information 
! from the input file for starting up that type of simulation
  CALL Get_Sim_Type

  WRITE(logunit,'(a,a,/)') 'Simulation type: ',sim_type

  ! Check what kind of simulation this is, and then call the appropriate routine
  ! that will load all the relevant information in from the input file
  IF (int_sim_type == sim_nvt .OR. int_sim_type == sim_nvt_min) THEN
     CALL NVTMC_Control
  ELSE IF (int_sim_type == sim_npt) THEN
     CALL NPTMC_Control
  ELSE IF (int_sim_type == sim_gcmc) THEN
     CALL GCMC_Control
  ELSE IF (int_sim_type == sim_gemc.OR. int_sim_type == sim_gemc_ig .OR. &
       int_sim_type == sim_gemc_npt) THEN
     CALL GEMC_Control
  ELSE IF (int_sim_type == sim_frag .OR. int_sim_type == sim_ring) THEN
     CALL NVT_MC_Fragment_Control
  ELSE IF (int_sim_type == sim_mcf) THEN
     CALL MCF_Control
  ELSE
     err_msg = ""
     err_msg(1) = 'Sim_Type unknown'
     err_msg(2) = logfile
     CALL Clean_Abort(err_msg,'Main')
  ENDIF

  ! Determine type of start
  CALL Get_Start_Type

  ! Determine if it is equilibration or production or test
  CALL Get_Run_Type

  WRITE(logunit,'(a,a,/)') 'starting type ',start_type

  IF( int_run_style == run_test ) THEN

             testname = molfile_name(1)
             IF(testname(:1) == 'c') testname = 'methane'
             IF(testname(:2) == 'pr') testname = 'propane'
             IF(testname(:7) == 'pentane') testname = 'pentane'
             IF(testname(:7) == 'pentano') testname = 'pentanol'

  END IF

  ! Initialize the counters for simulation
  
  DO i = 1, nbr_boxes
     CALL Initialize(i)
     CALL Reset(i)
  END DO

  ! Initialize the atom list and molecule list arrays
  molecule_list(:,:)%live = .FALSE.
  molecule_list(:,:)%which_box = 0
  
  atom_list(:,:,:)%exist = .FALSE.
  molecule_list(:,:)%cfc_lambda = 0.0_DP

  CBMC_Flag = .FALSE.

  DO is = 1, nspecies
     DO im = 1, nmolecules(is)
        locate(im,is) = im
     END DO
  END DO

! Initialize the counters for simulation

  DO i = 1, nbr_boxes
     CALL Initialize(i)
  END DO

  WRITE(*,*) 'Beginning Cassandra Simulation'
  WRITE(*,*) 

  IF (start_type == 'make_config') THEN
     ! Get required info from inputfile, then grow molecules using CBMC
     CALL Grow_Molecules

  ELSEIF (start_type == 'read_old') THEN
     ! Read in old coordinates and restart a new simulation, 
     ! Note that the counters have already been set to zero by the call to
     ! initialize and reset above.

     CALL Restart_From_Old

  ELSEIF (start_type == 'checkpoint') THEN
     ! Restart from checkpoint. Shall we verify match between stuff in cpt 
     ! and stuff input, or override with cpt info?

     CALL Read_Checkpoint
     
  ELSE
     err_msg = ""
     err_msg(1) = "make_config start type improperly specificed."
     CALL Clean_Abort(err_msg,'Main')
     
  ENDIF

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
     
     do i = 1, nbr_boxes
        CALL  Ewald_Reciprocal_Lattice_Vector_Setup(i)
     end do
     
     ! Here we can allocate the memory for cos_sum, sin_sum etc
     
     ALLOCATE(cos_sum(MAXVAL(nvecs),nbr_boxes))
     ALLOCATE(sin_sum(MAXVAL(nvecs),nbr_boxes))
     ALLOCATE(cos_sum_old(MAXVAL(nvecs),nbr_boxes))
     ALLOCATE(sin_sum_old(MAXVAL(nvecs),nbr_boxes))
     ALLOCATE(cos_sum_start(MAXVAL(nvecs),nbr_boxes))
     ALLOCATE(sin_sum_start(MAXVAL(nvecs),nbr_boxes))
     ALLOCATE(cos_mol(  MAXVAL(nvecs), SUM(nmolecules)))
     ALLOCATE(sin_mol(  MAXVAL(nvecs), SUM(nmolecules)))
     ! initialize these vectors
     cos_mol(:,:) = 0.0_DP
     sin_mol(:,:) = 0.0_DP
 
  END IF

  DO i = 1, nbr_boxes
     CALL Reset(i)
  END DO

  ! NR: I believe we have all the information to compute the 
  !     system charge and charge on each species

   Write(logunit,*) '******** Charge Neutrality Check *********#'

   q_tot_sys = 0.0_DP; q_mol=0.0_DP
   
   DO is = 1, nspecies
      q_mol = 0.0_DP 
      Do ia = 1, natoms(is)     
         q_mol = q_mol + nonbond_list(ia,is)%charge 
      END DO 
      Write(logunit,'(A,T15,2X,I4,4x,A,T45,4x,f12.8)')'Species', is, 'has charge', q_mol 
   ENDDO

   DO ibox = 1,nbr_boxes
      q_tot_sys = 0.0_DP 
      DO is = 1,nspecies
         DO im = 1,nmolecules(is)
            this_im  = locate(im,is)
            this_box = molecule_list(this_im,is)%which_box
            IF (ibox .eq. this_box) THEN
               DO ia = 1,natoms(is)
                  q_tot_sys = q_tot_sys + nonbond_list(ia,is)%charge 
               END DO 
            END IF
         END DO
      END DO
      Write(logunit,'(A,T13,4X,I4,4X,A,T45,4X,f12.8)')'Box ', ibox, 'has charge', q_tot_sys
      WRITE(logunit,*)

      IF (ABS(q_tot_sys) .gt. 0.000001) THEN
         IF ( .NOT. ((int_sim_type /=  sim_frag) .OR. (int_sim_type /= sim_ring)) ) THEN
            err_msg = ''
            err_msg(1) = 'Box has net charge'
            err_msg(2) = Int_To_String(ibox)
            CALL Clean_Abort(err_msg,'main.f90')
         END IF
      END IF
   END DO

   Write(logunit,*) '#********Charge Neutrality Check*********#'
 ! NR: At this point we can initialize random number generator

  CALL init_seeds(iseed1, iseed3)


  ! Assign a locate id to all the molecules
  ! This goes to the maximum number of molecules

  ! Compute total number of beads in each box

!  DO ibox = 1, nbr_boxes
!     CALL Compute_Beads(ibox)
!  END DO

  ! set up internal coordinates, the maximum distance from the COM and compute
  ! the total system energy to be displayed to the logfile.

  ! Internal coordinates
  CALL Get_Internal_Coords

 ! Calculate COM and distance of the atom farthest to the COM.

  DO is = 1, nspecies
     DO im = 1, nmolecules(is)
        this_im = locate(im,is)
        IF( .NOT. molecule_list(this_im,is)%live) CYCLE
        CALL Get_COM(this_im,is)
        CALL Compute_Max_Com_Distance(this_im,is)
     END DO
  END DO

 IF ((start_type == 'read_old') .or. (start_type == 'checkpoint')) THEN
  ! Fold the molecules. 
    DO is = 1, nspecies
     DO im = 1, nmolecules(is)
        this_im = locate(im,is)
        IF( .NOT. molecule_list(this_im,is)%live) CYCLE
        this_box = molecule_list(this_im,is)%which_box
        CALL Fold_Molecule(this_im,is,this_box) 
     END DO
    END DO
  END IF

  ! compute total system energy
  overlap = .FALSE.
     
  DO i = 1, nbr_boxes
       
  IF ( start_type == 'make_config' .and. int_vdw_sum_style(i) == vdw_cut_tail) &
                CALL Compute_Beads(i)
  		CALL Compute_Total_System_Energy(i,.TRUE.,overlap)
  END DO
           

  IF (overlap) THEN
     ! overlap was detected between two molecules so abort the program
     err_msg = ''
     err_msg(1) = 'Overlap detected in the starting structure'
     err_msg(2) = 'Start type '//start_type
     CALL Clean_Abort(err_msg,'Main')
  END IF
  
  ! write out the initial energy components to the log file
  
  DO i = 1,nbr_boxes
     
     WRITE(logunit,*) '*****************************************'
     WRITE(logunit,'(A36,2X,I2)') ' Starting energy components for box', i
     WRITE(logunit,*) ' Atomic units-Extensive'
     WRITE(logunit,*) '*****************************************'
     WRITE(logunit,*)
     
     WRITE(logunit,'(A,T30,F20.3)') 'Total system energy is' , energy(i)%total
     WRITE(logunit,'(A,T30,F20.3)') 'Intra molecular energy is', energy(i)%intra
     WRITE(logunit,'(A,T30,F20.3)') 'Bond energy is', energy(i)%bond
     WRITE(logunit,'(A,T30,F20.3)') 'Bond angle energy is', energy(i)%angle
     WRITE(logunit,'(A,T30,F20.3)') 'Dihedral angle energy is', energy(i)%dihedral
     WRITE(logunit,'(A,T30,F20.3)') 'Improper angle energy is', energy(i)%improper
     WRITE(logunit,'(A,T30,F20.3)') 'Intra nonbond vdw is', energy(i)%intra_vdw
     WRITE(logunit,'(A,T30,F20.3)') 'Intra nonbond elec is', energy(i)%intra_q
     WRITE(logunit,'(A,T30,F20.3)') 'Inter molecule vdw is', energy(i)%inter_vdw
     IF (int_vdw_sum_style(i) == vdw_cut_tail) THEN
        WRITE(logunit,'(A,T30,F20.3)') 'Long range correction is', energy(i)%lrc
     END IF
     WRITE(logunit,'(A,T30,F20.3)') 'Inter molecule q is', energy(i)%inter_q
     WRITE(logunit,'(A,T30,F20.3)') 'Reciprocal ewald is', energy(i)%ewald_reciprocal
     WRITE(logunit,'(A,T30,F20.3)') 'Self ewald is', energy(i)%ewald_self

     IF( int_charge_sum_style(i) == charge_ewald) WRITE(logunit,'(A,T30,I20)') 'Number of vectors is', nvecs(i)
     WRITE(logunit,*) '*******************************************'
     WRITE(logunit,*)

  END DO

  ! obtain the number of molecules of each species in a simulation box
  ! only if it isn't a new configuration.

  IF( (start_type /= 'make_config') .OR. (start_type /= 'read_old') ) THEN

    DO ibox = 1,nbr_boxes
       DO is = 1, nspecies
          CALL Get_Nmolecules_Species(ibox,is,nmol_is)
          nmols(is,ibox) = nmol_is
       END DO
    END DO

  END IF
  
  ! Populate the reservoir box if necessary

  DO is = 1, nspecies

     IF(species_list(is)%int_insert == int_igas) THEN
        first_res_update = .TRUE.
        CALL Update_Reservoir(is)
        first_res_update = .FALSE.
     END IF

  END DO


  ! Now begin simulation by calling apropriate driver routines
  ! Add routines here...

  IF (int_run_style == run_test .AND. n_mcsteps == 1) THEN

     DO
        
        CALL Angle_Distortion(ibox)
        IF(nsuccess(1,1)%angle .NE. 0) EXIT
        
     END DO
     
     IF(ndihedrals(1) .GT. 0) THEN
        
        DO
           
           CALL Rigid_Dihedral_Change(ibox)
           IF(nsuccess(1,1)%dihedral .NE. 0) EXIT
           
        END DO

     END IF

     CALL Compute_Total_System_Energy(1,.TRUE.,overlap)
     WRITE(logunit,'(A,T30,F20.3)')'Intra molecular energy is:', energy(1)%intra
     WRITE(logunit,'(A,T30,F20.3)')'Intra nonbond vdw is:', energy(1)%intra_vdw
     WRITE(logunit,'(A,T30,F20.3)')'Intra nonbond elec is:', energy(1)%intra_q
     
     OPEN(75,FILE='compare.dat')
     WRITE(75,'(T20,A,A)') 'Energy for a single', testname
     WRITE(75,'(A,T30,F20.3)')'Intra molecular energy is:', energy(1)%intra
     WRITE(75,'(A,T30,F20.3)')'Intra nonbond vdw is:', energy(1)%intra_vdw
     WRITE(75,'(A,T30,F20.3)')'Intra nonbond elec is:', energy(1)%intra_q
     CLOSE(75)
     
  ELSE IF (int_sim_type == sim_nvt .OR. int_sim_type == sim_nvt_min) THEN
     
     CALL NVTMC_Driver
     
  ELSE IF (int_sim_type == sim_npt) THEN
     
     CALL NPTMC_Driver
  

  ELSE IF (int_sim_type == sim_gcmc) THEN

     CALL GCMC_Driver

  ELSE IF (int_sim_type == sim_gemc .OR. int_sim_type == sim_gemc_ig .OR. &
     int_sim_type == sim_gemc_npt) THEN
     
     CALL GEMC_Driver

  ELSE IF (int_sim_type == sim_frag) THEN

     CALL NVT_MC_Fragment_Driver

  ELSE IF (int_sim_type == sim_ring) THEN

     CALL NVT_MC_Ring_Fragment

  END IF

  !***************************************************************************
  ! let us check if at the end of the simulation, the energies are properly updated
!$OMP PARALLEL

!  t_num = omp_get_thread_num()
!  write(logunit,*) 'Thread',t_num,"reporting for duty!"

!$OMP END PARALLEL

  write(logunit,*) '*********** Ending simulation *****************'
  write(logunit,*)
  write(logunit,*)


  IF(cpcollect) THEN
     DO is = 1,nspecies
        DO i = 1, nbr_boxes
           IF(ntrials(is,i)%cpcalc > 0) THEN
             write(logunit,'(A,I2,A,I2,A,F24.12)') 'Chemical potential for species',is,'in box',i, 'is', &
                                                 chpot(is,i) / ntrials(is,i)%cpcalc
             write(logunit,'(A,I2,A,I2,A,F24.12)') 'Ideal Chemical potential for species',is,'in box',i, 'is', &
                                                 chpotid(is,i) / ntrials(is,i)%cpcalc
           END IF
        END DO
     END DO
  END IF


  CALL DATE_AND_TIME(date,time,zone,end_values)
  CALL cpu_time(tot_time)

  tot_time = tot_time - start_time

  CALL Write_Subroutine_Times

  WRITE(logunit,*) '********************************************'
  WRITE(logunit,*)
  WRITE(logunit,*)
  
  WRITE(logunit,*) 'Program execution time'

  WRITE(*,*)
  WRITE(*,*) 'Cassandra simulation complete'
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
  
END PROGRAM Main
