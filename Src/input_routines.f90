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

MODULE Input_Routines

  !***************************************************************************
  ! This module contains a collection of subroutines used to obtain information from
  ! the input file for all the different types of simulations
  !
  ! Called by
  !
  !    gcmc_control
  !    gemc_control
  !    main
  !    mcf_control
  !    nptmc_control
  !    nvtmc_control
  !    nvt_mc_fragment_control
  !
  ! Revision history
  !
  !    12/10/13  : Beta version
  !***************************************************************************

  USE Run_Variables
  USE IO_Utilities
  USE File_Names
  USE Type_Definitions
  USE Random_Generators, ONLY : s1


  IMPLICIT NONE

CONTAINS

!********************************************************************************
SUBROUTINE Get_Runname
!********************************************************************************
! This routine opens the input file and determines the name of the run.
! Input file name format is 
! # Key_Word
! Variables associated with Key_Word

! The routine searches for the keyword "# Run_Name" and then reads the necessary 
! information underneath the key word. Spaces, blank lines and ! characters are 
! ignored.
!********************************************************************************

  INTEGER :: ierr,line_nbr,nbr_entries
  CHARACTER(120) :: line_string, line_array(20)

  
!********************************************************************************
! Determine the name of the run from input file.
! All output files will have this name.
!********************************************************************************
  REWIND(inputunit)

  ierr = 0
  line_nbr = 0

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(inputunit,line_string,ierr)

     IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading name of run."
        CALL Clean_Abort(err_msg,'Read_inputfile')
     END IF

     IF (line_string(1:10) == '# Run_Name') THEN
        line_nbr = line_nbr + 1
        
        CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

! Assign the first entry on the line to the name of the run
        run_name = line_array(1)
        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

! No name specified so abort
        err_msg = ""
        err_msg(1) = 'No run name specified in inputfile'
        CALL Clean_Abort(err_msg,'Get_Runname')

        EXIT

     ENDIF

  ENDDO

END SUBROUTINE Get_Runname

!*************************************************************************************
SUBROUTINE Get_Nspecies
!*************************************************************************************
!
! This routine reads in the number of species to be simulated form the input file.
! It then allocates all arrays that depend only on nspecies
! 
!*************************************************************************************

  INTEGER :: ierr,line_nbr,nbr_entries, i
  CHARACTER(120) :: line_string, line_array(20)
!*************************************************************************************
  REWIND(inputunit)

! determine the number of species to be simulated
  ierr = 0
  line_nbr = 0

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(inputunit,line_string,ierr)

     IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading number of species."
        CALL Clean_Abort(err_msg,'Get_Nspecies')
     END IF

     IF (line_string(1:13) == '# Nbr_Species') THEN
        line_nbr = line_nbr + 1
        
        CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

! Assign the first entry on the line to the number of species.
        nspecies = String_To_Int(line_array(1))

        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

! No name specified so abort
        err_msg = ""
        err_msg(1) = 'Number of species not specified with keyword # Nbr_Species'
        CALL Clean_Abort(err_msg,'Get_Nspecies')

        EXIT

     ENDIF

  ENDDO

  ALLOCATE( molfile_name(nspecies),Stat = AllocateStatus )
  IF (AllocateStatus /= 0 ) THEN
     write(*,*)'memory could no tbe allocated for molfile_name array'
     write(*,*)'stopping'
     STOP
  END IF
  ALLOCATE( nmolecules(nspecies), natoms(nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0 ) THEN
     write(*,*)'memory could no tbe allocated for nmolecules or natoms array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE (nring_atoms(nspecies), nexo_atoms(nspecies), Stat = AllocateStatus)
  IF (AllocateStatus /= 0 ) THEN
     write(*,*)'memory could no tbe allocated for nmolecules or natoms array'
     write(*,*)'stopping'
     STOP
  END IF


  ALLOCATE( nbonds(nspecies), nangles(nspecies),Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for nbonds or nangles array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE( ndihedrals(nspecies), nimpropers(nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for ndihedrals or nimpropers array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE(nbr_bond_params(nspecies),Stat = AllocateStatus)
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for nbr_bond_params array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE(nbr_angle_params(nspecies),Stat = AllocateStatus)

  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for molfile_name array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE(nbr_dihedral_params(nspecies),Stat = AllocateStatus)
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for nbr_dihedral_params array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE(nbr_improper_params(nspecies),Stat = AllocateStatus)
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for nbr_improper_params array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE(nbr_vdw_params(nspecies),Stat = AllocateStatus)
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for nbr_vdw_params array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE(nfragments(nspecies),Stat = AllocateStatus)
  IF (AllocateStatus /=0 ) THEN
     write(*,*) 'memory could not be allocated for nfragments array'
     write(*,*) 'stopping'
     STOP
  END IF

  ALLOCATE(fragment_bonds(nspecies), Stat = AllocateStatus)
  IF (AllocateStatus /= 0 ) THEN
     write(*,*) 'memroy could not be allocated for fragment_bonds array'
     write(*,*) 'stopping'
     STOP
  END IF

! Initialize everything
  molfile_name = ''
  nmolecules = 0
  natoms = 0
  nbonds = 0
  nangles = 0
  ndihedrals = 0
  nimpropers = 0
  nbr_bond_params = 0
  nbr_angle_params = 0
  nbr_dihedral_params = 0
  nbr_improper_params = 0
  nbr_vdw_params = 0
  nfragments = 0
  fragment_bonds = 0

END SUBROUTINE Get_Nspecies
!********************************************************************************


!********************************************************************************
SUBROUTINE Get_Sim_Type
!********************************************************************************
! This routine opens the input file and determines the type of simulation.
!
! Allowed simulation types:
! NVT_MC: canonical ensmeble Monte Carlo
! GCMC: grand canonical ensemble Monte Carlo
! HMC: hybrid Monte Carlo
! Others to be added...

! Input file name format is 
! # Key_Word
! Variables associated with Key_Word

! The routine searches for the keyword "# Sim_Type" and then reads the necessary 
! information underneath the key word. Spaces, blank lines and ! characters are 
! ignored.
!********************************************************************************
  INTEGER :: ierr,line_nbr,nbr_entries
  CHARACTER(120) :: line_string, line_array(20)

!********************************************************************************
  REWIND(inputunit)

  ierr = 0
  line_nbr = 0

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(inputunit,line_string,ierr)

     IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading simulation type."
        CALL Clean_Abort(err_msg,'Get_Sim_Type')
     END IF

     IF (line_string(1:10) == '# Sim_Type') THEN
        line_nbr = line_nbr + 1
         
        CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

! Assign the first entry on the line to simulation type
        sim_type = line_array(1)

        line_nbr = line_nbr + 1
   
        IF(sim_type == 'GCMC') THEN

           int_sim_type = sim_gcmc
     !     CALL Parse_String(inputunit,line_nbr,2,nbr_entries,line_array,ierr)
           !     tmmc_update = String_To_Int(line_array(1))
     !     tmmc_input = line_array(2)
           
        END IF

        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

! No sim type specified so abort
        err_msg = ""
        err_msg(1) = 'No simulation start type  specified in inputfile'
        CALL Clean_Abort(err_msg,'Get_Sim_Type')

        EXIT

     ENDIF

  ENDDO

  IF(sim_type == 'NVT_MC') THEN
     int_sim_type = sim_nvt
  ELSEIF(sim_type == 'NVT_MIN') THEN
     int_sim_type = sim_nvt_min
  ELSEIF(sim_type == 'NPT_MC') THEN
     int_sim_type = sim_npt
  ELSEIF(sim_type == 'GEMC') THEN
     int_sim_type = sim_gemc
  ELSEIF(sim_type == 'GEMC_NPT') THEN
     int_sim_type = sim_gemc_npt
  ELSEIF(sim_type == 'NVT_MC_Fragment') THEN
     int_sim_type = sim_frag
  ELSEIF(sim_type == 'NVT_MC_Ring_Fragment') THEN
     int_sim_type = sim_ring
  ELSEIF(sim_type == 'MCF_Gen') THEN
     int_sim_type = sim_mcf
  END IF
  
END SUBROUTINE Get_Sim_Type



!********************************************************************************
SUBROUTINE Get_Pair_Style
!********************************************************************************
! This routine opens the input file and determines the type of vdw and charge 
! interaction models to use.

! The routine searches for the keywords "# VDW_Style", "# Charge_Style", and
! "# Neighbor_Style then reads the necessary information underneath the key word. 
! Spaces, blank lines and ! characters are ignored.
!
! 09/08/10 (JS) : The pair style also reads in if '# Pair_Energy' is true or false.
!                 This will alert the code if pair interaction energy arrays
!                 need to be stored. 
!********************************************************************************
  INTEGER :: ierr,line_nbr,nbr_entries, iassign, ibox, k
  CHARACTER(120) :: line_string, line_array(20)

  REAL(DP), ALLOCATABLE :: ewald_tol(:)

!********************************************************************************
  REWIND(inputunit)

  ierr = 0
  line_nbr = 0
  iassign = 0
  roff_charmm(:) = 0.0_DP
  roff_switch(:) = 0.0_DP
  rcut_vdw(:) = 0.0_DP
  rcut_coul(:) = 0.0_DP
  rcut9(:) = 0.0_DP
  rcut3(:) = 0.0_DP

  l_pair_nrg = .FALSE.

  ALLOCATE(l_half_len_cutoff(nbr_boxes))
  l_half_len_cutoff = .FALSE.

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(inputunit,line_string,ierr)

     IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading pairstyle."
        CALL Clean_Abort(err_msg,'Get_Pairstyle')
     END IF

     IF (line_string(1:11) == '# VDW_Style') THEN
        DO ibox = 1,nbr_boxes
           line_nbr = line_nbr + 1
        
           ! Read the type of VDW model, the cutoff method, and the parameters associated 
           ! with this cutoff method. Minimum of 3 values must be listed. 

           CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

           ! Assign the first entry on the line to the name of the potential, then next to the 
           ! way it will be summed / truncated, and the remaining to parameters associated with
           ! the sum method
           vdw_style(ibox) = line_array(1)
           WRITE(logunit,'(A,2x,A,A,I3)') '   VDW style used is: ',vdw_style(ibox), 'in box:', ibox

           IF (vdw_style(ibox) /= 'NONE') THEN
              int_vdw_style(ibox) = vdw_lj
              vdw_sum_style(ibox) = line_array(2)
              WRITE(logunit,'(A,2x,A,A,I3)') '   VDW sum style is: ',vdw_sum_style(ibox), 'in box:', ibox

              IF (vdw_sum_style(ibox) == 'CHARMM') THEN
                 int_vdw_sum_style(ibox) = vdw_charmm
                 ron_charmm(ibox) = String_To_Double(line_array(3))
                 roff_charmm(ibox) = String_To_Double(line_array(4))
                 WRITE(logunit,'(A,2x,F7.3, A)') '    r_on = ',ron_charmm(ibox), '   Angstrom'
                 WRITE(logunit,'(A,2x,F7.3,A)') '    r_off = ',roff_charmm(ibox), '   Angstrom'

              ELSEIF (vdw_sum_style(ibox) == 'cut_switch') THEN
                 int_vdw_sum_style(ibox) = vdw_cut_switch
                 ron_switch(ibox) = String_To_Double(line_array(3))
                 roff_switch(ibox) = String_To_Double(line_array(4))
                 WRITE(logunit,'(A,2x,F7.3,A)') '   r_on =', ron_switch(ibox), ' Angstrom'
                 WRITE(logunit,'(A,2x,F7.3,A)') '  r_off =', roff_switch(ibox), ' Angstrom'

              ELSEIF (vdw_sum_style(ibox) == 'cut') THEN
                 int_vdw_sum_style(ibox) = vdw_cut
                 rcut_vdw(ibox) = String_To_Double(line_array(3))
                 IF ( nbr_entries == 4 ) THEN
                    ! a fourth entry exists indicating whether the cutoff is half of
                    ! the box length 
                    
                    IF (line_array(4) == 'TRUE' .OR. line_array(4) == 'true') THEN
                       
                       l_half_len_cutoff(ibox) = .TRUE.
                       
                       ! for now assume that the box is cubic
                       rcut_vdw(ibox) = 0.5_DP * box_list(ibox)%length(1,1)
                       
                       WRITE(logunit,*)
                       WRITE(logunit,'(A,2x,I5)') 'For box ', ibox
                       WRITE(logunit,*) 'Cutoffs are set to half of the box length'
                       
                    END IF

                 END IF

                 WRITE(logunit,'(A,2x,F7.3, A)') '    rcut = ',rcut_vdw(ibox), '   Angstrom'

              ELSEIF (vdw_sum_style(ibox) == 'cut_tail') THEN
                 int_vdw_sum_style(ibox) = vdw_cut_tail
                 rcut_vdw(ibox) = String_To_Double(line_array(3))
                
                 IF ( nbr_entries == 4 ) THEN
                    ! a fourth entry exists indicating whether the cutoff is half of
                    ! the box length 
                    
                    IF (line_array(4) == 'TRUE' .OR. line_array(4) == 'true') THEN
                       
                       l_half_len_cutoff(ibox) = .TRUE.
                       
                       ! for now assume that the box is cubic
                       rcut_vdw(ibox) = 0.5_DP * box_list(ibox)%length(1,1)
                       
                       WRITE(logunit,*)
                       WRITE(logunit,'(A,2x,I5)') 'For box ', ibox
                       WRITE(logunit,*) 'Cutoffs are set to half of the box length'
                       
                    END IF

                 END IF

                 WRITE(logunit,'(A,2x,F7.3, A)') '    rcut = ',rcut_vdw(ibox), '   Angstrom'

                 rcut3(ibox) = rcut_vdw(ibox) * rcut_vdw(ibox) * rcut_vdw(ibox)
                 rcut9(ibox) = rcut3(ibox) * rcut3(ibox) * rcut3(ibox)

              ELSEIF (vdw_sum_style(ibox) == 'cut_shift') THEN
                 int_vdw_sum_style(ibox) = vdw_cut_shift
                 rcut_vdw(ibox) = String_To_Double(line_array(3))
                 WRITE(logunit,'(A,2x,F7.3, A)') '    rcut = ',rcut_vdw(ibox), '   Angstrom'

              ELSEIF (vdw_sum_style(ibox) == 'minimum_image') THEN
                 int_vdw_sum_style(ibox) = vdw_minimum
                 WRITE(logunit,'(A)') 'Minimum image convention used for VDW'

              ELSEIF (vdw_sum_style(ibox) == 'mie') THEN
                 int_vdw_sum_style(ibox) = vdw_mie
		 rcut_vdw(ibox) = String_To_Double(line_array(3))
                 WRITE(logunit,'(A,2x,F7.3, A)') '    rcut = ',rcut_vdw(ibox), '   Angstrom'
                 WRITE(logunit,'(A)') 'Mie potential used for VDW'

              ELSE
                 err_msg(1) = 'Improper specification of vdw_sum_style'
                 CALL Clean_Abort(err_msg,'Get_Pairstyle')
              ENDIF
              IF (rcut_vdw(ibox) > MIN(box_list(ibox)%face_distance(1)/2.0_DP, &
                   box_list(ibox)%face_distance(2)/2.0_DP, box_list(ibox)%face_distance(3)/2.0_DP)) THEN

                     err_msg = ""
                     err_msg(1) = 'Initial cutoff greater than minimum box length'
                     err_msg(2) = 'For box'
                     err_msg(3) = Int_To_String(ibox)
                     CALL Clean_Abort(err_msg,'Get_Pairstyle')

    	      ENDIF



           ELSE
 
              int_vdw_style(ibox) = vdw_none

           ENDIF

           WRITE(logunit,*) 'VDW style properly input'
           WRITE(logunit,*)
           iassign = iassign + 1
           ! Test if both vdw and coulomb stuff read OK. If so, done.
           IF (iassign == 2*nbr_boxes) EXIT

        END DO

     ELSEIF (line_string(1:14) == '# Charge_Style') THEN
        
        DO ibox = 1,nbr_boxes

           line_nbr = line_nbr + 1
        
           ! Read the charge summation style. Three parameters must be specified
           ! First indicates type of the cutoff method and subsequent lines indicate the parameters.

           CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

           ! Assign the first entry on the line to the name of the way charged interactions are treated
           ! then next is the way the charges are summed. Following that are the parameters associated
           ! with this particular method of summing charges. 
           charge_style(ibox) = line_array(1)
           WRITE(logunit,'(A,2x,A,A,I3)') '   Charge style used is: ',charge_style(ibox), 'in box:', ibox

           IF (charge_style(ibox) /= 'NONE') THEN
                 int_charge_style(ibox) = charge_coul

                 charge_sum_style(ibox) = line_array(2)
                 WRITE(logunit,'(A,2x,A,A,I3)') '    Charge sum style is ',charge_sum_style(ibox), 'in box:', ibox

                 IF (charge_sum_style(ibox) == 'cut') THEN
                    int_charge_sum_style(ibox) = charge_cut
                    rcut_coul(ibox) = String_To_Double(line_array(3))

                    IF (l_half_len_cutoff(ibox)) THEN
                       rcut_coul(ibox) = rcut_vdw(ibox)
                    ELSE
                       rcut_coul(ibox) = String_To_Double(line_array(3))
                    END IF

                    WRITE(logunit,'(A,2x,F7.3, A)') '    rcut = ',rcut_coul(ibox), '   Angstrom'

                 ELSEIF (charge_sum_style(ibox) == 'Ewald') THEN
                    int_charge_sum_style(ibox) = charge_ewald
                    rcut_coul(ibox) = String_To_Double(line_array(3))

                    IF (l_half_len_cutoff(ibox)) THEN
                       rcut_coul(ibox) = rcut_vdw(ibox)
                    ELSE
                       rcut_coul(ibox) = String_To_Double(line_array(3))
                    END IF

                    IF (ibox == 1) THEN
                       ALLOCATE(ewald_tol(nbr_boxes), ewald_p_sqrt(nbr_boxes))
                       ALLOCATE(ewald_p(nbr_boxes))
                       ALLOCATE(alpha_ewald(nbr_boxes) , h_ewald_cut(nbr_boxes) )
                       ALLOCATE(alphal_ewald(nbr_boxes) )
                    END IF

                    ewald_tol(ibox) = String_To_Double(line_array(4))

                    DO k = 1, ibox - 1

                       IF ( ABS(ewald_tol(ibox) - ewald_tol(k)) > 1.0D-10 ) THEN
                          err_msg(1) = "Ewald accuracy is set differently for"
                          err_msg(2) = "Box "//Int_To_String(ibox)
                          err_msg(3) = "and box "//Int_To_String(k)
                          CALL Clean_Abort(err_msg,'Get_Pair_Style')
                       END IF

                    END DO

                    IF(ewald_tol(ibox) .GT. 1.0_DP) THEN
                       err_msg(1) = "Ewald tolerance too low for box"//Int_To_String(ibox)
                       CALL Clean_Abort(err_msg,'Get_Pair_Style')
                    END IF

                    ewald_p(ibox) = -DLOG(ewald_tol(ibox))
                    ewald_p_sqrt(ibox) = DSQRT(ewald_p(ibox))
                    
                    alpha_ewald(ibox) = ewald_p_sqrt(ibox) / rcut_coul(ibox)
                    
                    h_ewald_cut(ibox) = 2.0_DP * ewald_p(ibox) / rcut_coul(ibox)

                                        

                    WRITE(logunit,'(X,A,F7.3,A)') '   Ewald real space cutoff is ', &
                       rcut_coul(ibox), ' Angstroms.'
                    WRITE(logunit, '(X,A,F7.3,A)') ' Ewald real space parameter is ', &
                         alpha_ewald(ibox), ' inverse Angstroms'
                    WRITE(logunit,'(X,A,F7.4,A)') '   Ewald reciprocal cutoff is ', &
                         h_ewald_cut(ibox), ' inverse Angstroms'
 
               
                 ELSEIF (charge_sum_style(ibox) == 'minimum_image') THEN
                    int_charge_sum_style(ibox) = charge_minimum
                    IF (int_vdw_sum_style(ibox) /= vdw_minimum .AND. int_vdw_style(ibox) /= vdw_none) THEN
                       err_msg=""
                       err_msg(1) = 'Minimum image requires both vdw and q-q to be so-specified'
                       CALL Clean_Abort(err_msg,'Get_Pairstyle')
                    ELSE
                       WRITE(logunit,'(A)') 'Minimum image convention used for VDW'
                    ENDIF
                 ELSE
                    err_msg(1) = 'charge_sum_style not properly specified'
                    CALL Clean_Abort(err_msg,'Get_Pairstyle')
                 ENDIF

                 IF (rcut_coul(ibox) > MIN(box_list(ibox)%face_distance(1)/2.0_DP, &
                   box_list(ibox)%face_distance(2)/2.0_DP, box_list(ibox)%face_distance(3)/2.0_DP)) THEN

                     err_msg = ""
                     err_msg(1) = 'Initial cutoff greater than minimum box length'
                     err_msg(2) = 'For box'
                     err_msg(3) = Int_To_String(ibox)
                     CALL Clean_Abort(err_msg,'Get_Pairstyle')

                 ENDIF


           ELSE

              int_charge_style(ibox) = charge_none

           ENDIF

           WRITE(logunit,*) 'Charge style properly input'
           WRITE(logunit,*)                 
           iassign = iassign + 1
           ! Test if both vdw and coulomb stuff read OK. If so, done.

        END DO

        IF(ALLOCATED(ewald_tol)) DEALLOCATE(ewald_tol)

    ELSE IF (line_string(1:13) == '# Pair_Energy') THEN
        line_nbr = line_nbr + 1
        CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

        WRITE(logunit,*)
        IF (line_array(1) == 'TRUE' .OR. line_array(1) == 'true') THEN
           l_pair_nrg = .TRUE.
           WRITE(logunit,*) 'Pair interaction energy array storage enabled'
        END IF
        WRITE(logunit,*)

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN
        IF (iassign == 2*nbr_boxes) EXIT

! No pairstyle specified so abort
        err_msg = ""
        err_msg(1) = 'No simulation start type  specified in inputfile'
        CALL Clean_Abort(err_msg,'Get_Pairstyle')
     ENDIF
     
  ENDDO

  !Now determine the mixing rule to use
   CALL Get_Mixing_Rules

  IF (.NOT. l_pair_nrg) THEN
     WRITE(logunit,*)
     WRITE(logunit,*) 'Pair interaction energy arrays will not be stored'
     WRITE(logunit,*) 'Energy calculations will be done in a standard way'
     WRITE(logunit,*)
  END IF
  
END SUBROUTINE Get_Pair_Style



!********************************************************************************
SUBROUTINE Get_Mixing_Rules
!********************************************************************************
! The routine searches for the keyword "# Mixing_Rules" and then reads the necessary 
! information underneath the key word. Spaces, blank lines and ! characters are 
! ignored. If no mixing rule is specified, Lorentz-Berthelot is used as default.
!********************************************************************************
  INTEGER :: ierr,line_nbr,nbr_entries
  CHARACTER(120) :: line_string, line_array(20)

!********************************************************************************
  REWIND(inputunit)

  ierr = 0
  line_nbr = 0

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(inputunit,line_string,ierr)

     IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading mixinf rules."
        CALL Clean_Abort(err_msg,'Get_Mixing_Rules')
     END IF

     IF (line_string(1:13) == '# Mixing_Rule') THEN
        line_nbr = line_nbr + 1
        
        CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

! Assign the first entry on the line to the mixing rule
        mix_rule= line_array(1)
        
        IF (mix_rule == 'LB') THEN
           WRITE(logunit,'(A)') 'Lorentz-Berthelot mixing rule specified'
        ELSEIF (mix_rule == 'geometric') THEN
           WRITE(logunit,'(A)') 'Geometric mixing rule specified'
        ELSEIF (mix_rule == 'custom') THEN
           WRITE(logunit,'(A)') 'Custom mixing rule specified'
        ELSE
           err_msg(1) = 'Mixing rule not supported'
           err_msg(2) = mix_rule
           err_msg(3) = 'Available options are'
           err_msg(4) = 'LB or geometric'
           CALL Clean_Abort(err_msg,'Get_Mixing_Rules')
        ENDIF

        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

        WRITE(logunit,*) 'No mixing rule specified. Using Lorentz-Berthelot'
        mix_rule = 'LB'

        EXIT

     ENDIF

  ENDDO

END SUBROUTINE Get_Mixing_Rules

!********************************************************************************
SUBROUTINE Get_Molecule_Info
!********************************************************************************
! This routine opens the input file and reads connectivity information for
! each molecule. It determines the number of atoms, bonds, angles, dihedrals
! and impropers in each molecule. It the allocates associated arrays amd populates 
! most of the atom_class, bond_class, angle_class
! dihedral_class, improper_class and nonbond_class fields. 
!
! The routine searches for the keyword "# Molecule_Files" and then 
! for each species, it reads the name of the molecular connectivity file, opens
! that file and loads necessary information. Only the homegrown molecular 
! connectivity file format is supported.
!********************************************************************************

  INTEGER :: ierr,line_nbr,nbr_entries, i, openstatus, is, max_index
  INTEGER :: mcf_index(5), dummy
  CHARACTER(120) :: line_string, line_array(20)

!********************************************************************************
! determine the type of molecule input and connectivity
  REWIND(inputunit)

  ierr = 0
  line_nbr = 0
  
  input_file_loop:DO
     line_nbr = line_nbr + 1

     CALL Read_String(inputunit,line_string,ierr)

     IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading molecular info."
        CALL Clean_Abort(err_msg,'Get_Molecule_Info')
     END IF

     molecule_file_string:IF (line_string(1:16) == '# Molecule_Files') THEN
        line_nbr = line_nbr + 1

        ! next lines must contain the molecule file names of each species in order and
        ! the number of molecules. If this is an open system simulation (i.e. GCMC)
        ! the maximum expected number of molecules should be listed, as this size
        ! array will be allocated. We also determine the maximum number of molecules
        ! of a given species and use this to allocate arrays

        

           species_loop:DO i=1,nspecies

              
                 CALL Parse_String(inputunit,line_nbr,2,nbr_entries,line_array,ierr)


              IF (ierr .NE. 0) THEN
                 err_msg = ""
                 err_msg(1) = "Error reading molecular connectivity file."
                 err_msg(2) = "check that number of species and number of files match"
                 CALL Clean_Abort(err_msg,'Get_Molecule_File_Type')
              END IF

              ! assign the name of the molecular connectivity file, max number
              ! of molecules, and starting number of molecules for this species

              molfile_name(i) = line_array(1)
              nmolecules(i) = String_To_Int(line_array(2))
          
              WRITE(logunit,*)
              WRITE(logunit,*) 'Reading molecular connectivity information'
              WRITE(logunit,*) 'Species: ',i
              WRITE(logunit,*) 'Molecular connectivity file: ',molfile_name(i)

              ! Open the file and determine how many atoms, bonds, angles, dihedrals and impropers
              ! this molecule has
              OPEN(UNIT=molfile_unit,FILE=molfile_name(i),STATUS="OLD",IOSTAT=openstatus,ACTION="READ")

              IF (openstatus .NE. 0) THEN
                 err_msg = ""
                 err_msg(1) = "Unable to open molecular connectivity file."
                 CALL Clean_Abort(err_msg,'Get_Molecule_Info')
              ENDIF

              REWIND(molfile_unit)

              mcf_index = 0

              mcf_read_loop:DO 
                 CALL Read_String(molfile_unit,line_string,ierr)

                 IF (ierr == 0) THEN

                    IF (line_string(1:11) == '# Atom_Info') THEN
                       CALL Read_String(molfile_unit,line_string,ierr)
                       natoms(i) = String_To_Int(line_string)
                       WRITE(logunit,*) '  ', &
                            TRIM(Int_To_String(natoms(i))), ' atom(s) specified.'

                       mcf_index(1) = 1

                    ELSEIF (line_string(1:11) == '# Bond_Info') THEN
                       CALL Read_String(molfile_unit,line_string,ierr)
                       nbonds(i) = String_To_Int(line_string)
                       WRITE(logunit,*) '  ', &
                            TRIM(Int_To_String(nbonds(i))), ' bond(s) specified.'
                       mcf_index(2) = 1

                    ELSEIF (line_string(1:12) == '# Angle_Info') THEN
                       CALL Read_String(molfile_unit,line_string,ierr)
                       nangles(i) = String_To_Int(line_string)
                       WRITE(logunit,*) '  ', &
                            TRIM(Int_To_String(nangles(i))), ' angle(s) specified.'
                       mcf_index(3) = 1

                    ELSEIF (line_string(1:15) == '# Dihedral_Info') THEN
                       CALL Read_String(molfile_unit,line_string,ierr)
                       ndihedrals(i) = String_To_Int(line_string)
                       WRITE(logunit,*) '  ', &
                            TRIM(Int_To_String(ndihedrals(i))), ' dihedral(s) specified.'
                       mcf_index(4) = 1

                    ELSEIF (line_string(1:15) == '# Improper_Info') THEN
                       CALL Read_String(molfile_unit,line_string,ierr)
                       nimpropers(i) = String_To_Int(line_string)
                       WRITE(logunit,*) '  ', &
                            TRIM(Int_To_String(nimpropers(i))), ' improper(s) specified.'
                       mcf_index(5) = 1

                    ELSEIF (line_string(1:15) == '# Fragment_Info') THEN
                       CALL Read_String(molfile_unit,line_string,ierr)
                       nfragments(i) = String_To_Int(line_string)
                       WRITE(logunit,*) '  ', &
                            TRIM(Int_To_String(nfragments(i))), ' fragments specified.'

                    ELSEIF (line_string(1:23) == '# Fragment_Connectivity') THEN
                       CALL Read_String(molfile_unit,line_string,ierr)
                       fragment_bonds(i) = String_To_Int(line_string)
                       WRITE(logunit,*) '  ', &
                            TRIM(Int_to_String(fragment_bonds(i))), '  fragment bonds specified.'

                    END IF

                 ELSE
                    ! Make sure everything has been specified
                    IF (SUM(mcf_index) .NE. 5) THEN
                       err_msg = ""
                       err_msg(1) =  'Error! In mcf file '
                       err_msg(2) = molfile_name(i)

                       IF (mcf_index(1) .NE. 1) err_msg(3) = '   natoms  field not present'
                       IF (mcf_index(2) .NE. 1) err_msg(4) = '   nbonds  field not present'
                       IF (mcf_index(3) .NE. 1) err_msg(5) = '   nangles  field not present'
                       IF (mcf_index(4) .NE. 1) err_msg(6) = '   ndihedrals  field not present'
                       IF (mcf_index(5) .NE. 1) err_msg(7) = '   nimpropers  field not present'
                       CALL Clean_Abort(err_msg,'Get_Init_Params')
                    END IF
                    ! Everything properly specified
                    EXIT  !Exit when EOF reached
                 END IF

              END DO mcf_read_loop

              CLOSE(molfile_unit)

           ENDDO species_loop

           EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

! No connectivity file specified so abort
        err_msg = ""
        err_msg(1) = 'Molecular connectivity files improperly specified'
        CALL Clean_Abort(err_msg,'Get_Molecule_File_Type')

        EXIT

     ENDIF molecule_file_string


  ENDDO input_file_loop

  WRITE(logunit,*)
  WRITE(logunit,'(A7,2x,A13)') 'Species', 'Nbr molecules'
  WRITE(logunit,'(A7,2x,A13)') '-------', '-------------'
  DO i=1,nspecies
     WRITE(logunit,'(I6,2x,I13)') i,nmolecules(i) 
  ENDDO
  WRITE(logunit,*)

  ! Allocate arrays that depend on nmolecules, natoms, and nspecies
  ! N.B.: MAXVAL instrinsic function selects the largest value from an array

  ALLOCATE( atom_list(MAXVAL(natoms), MAXVAL(nmolecules), nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for atom_list array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE( nonbond_list(MAXVAL(natoms), nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for nonbond_list array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE( ring_atom_ids(MAXVAL(natoms), nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for ring_atom_ids array'
     write(*,*)'stopping'
     STOP
  END IF
  
  ALLOCATE( exo_atom_ids(MAXVAL(natoms), nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for exo_atom_ids array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE( bond_list(MAXVAL(nbonds), nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for bond_list array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE( angle_list(MAXVAL(nangles), nspecies),Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for angle_list array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE( dihedral_list(MAXVAL(ndihedrals), nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for dihedral_list array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE( improper_list(MAXVAL(nimpropers), nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for improper_list array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE( molecule_list(MAXVAL(nmolecules), nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for molecule_list array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE( locate(MAXVAL(nmolecules),nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for locate array'
     write(*,*)'stopping'
     STOP
  END IF

  IF (l_pair_nrg) THEN
     ALLOCATE( pair_nrg_vdw(SUM(nmolecules),SUM(nmolecules)), Stat = AllocateStatus)
     IF (AllocateStatus /= 0 ) THEN
        write(*,*) 'memmory could not be allocated for pair_nrg_vdw array'
        write(*,*) 'aborting'
        STOP
     END IF
     
     ALLOCATE( pair_nrg_qq(SUM(nmolecules),SUM(nmolecules)), Stat = AllocateStatus)
     IF (AllocateStatus /= 0 ) THEN
        write(*,*) 'memmory could not be allocated for pair_nrg_qq array'
        write(*,*) 'aborting'
        STOP
     END IF
  END IF

  ALLOCATE( species_list(nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for species_list array'
     write(*,*)'stopping'
     STOP
  END IF

  max_index = MAX(MAXVAL(nbonds),MAXVAL(nangles),MAXVAL(ndihedrals),MAXVAL(nimpropers))

  ALLOCATE( internal_coord_list(max_index, MAXVAL(nmolecules), nspecies), Stat = AllocateStatus )
  IF (AllocateStatus /= 0) THEN
     write(*,*)'memory could no tbe allocated for internal_coord_list array'
     write(*,*)'stopping'
     STOP
  END IF

  ALLOCATE(internal_coord_list_old(max_index), Stat = AllocateStatus)
  IF (AllocateStatus /= 0) THEN
     err_msg = ''
     err_msg(1) = 'memory could no tbe allocated for internal_coord_list_old array'
     CALL Clean_Abort(err_msg,'Get_Molecule_Info')
  END IF 
 

  ALLOCATE(frag_list(MAXVAL(nfragments),nspecies), Stat = AllocateStatus)
  IF (AllocateStatus /= 0) THEN
     err_msg = ''
     err_msg(1) = 'memory could not be allocated for frag_list array'
     CALL Clean_Abort(err_msg,'Get_Molecule_Info')
  END IF

  ALLOCATE(fragment_bond_list(MAXVAL(fragment_bonds),nspecies), Stat = AllocateStatus)
  IF (AllocateStatus /= 0 ) THEN
     err_msg = ''
     err_msg(1) = 'memory could not be allocated for fragment_bond_list'
     CALL Clean_Abort(err_msg,'Get_Molecule_Info')
  END IF

  ALLOCATE(res_file(MAXVAL(nfragments),nspecies), Stat = AllocateStatus)
  IF ( AllocateStatus /= 0 ) THEN
     err_msg = ''
     err_msg(1) = 'memorgy could not be alloaced for reservoir files'
  END IF

  ALLOCATE(zig_calc(nspecies))
  zig_calc(:) = .FALSE.

  ! Loop over all species and load information from mfc files into list arrays
  species_list(:)%fragment = .FALSE.
  
  DO is=1,nspecies

        OPEN(UNIT=molfile_unit,FILE=molfile_name(is),STATUS="OLD",IOSTAT=openstatus,ACTION="READ")
        IF (openstatus .NE. 0) THEN
           err_msg = ""
           err_msg(1) = "Unable to open molecular connectivity file."
           err_msg(2) = molfile_name(is)
           CALL Clean_Abort(err_msg,'Get_Molecule_Info')
        ENDIF

        CALL Get_L_Coul_CBMC(is)
        ! skip the insertion style for a fragment generation
        IF ( .NOT. (int_sim_type == sim_frag .OR. int_sim_type == sim_ring )) THEN
           CALL Get_Insertion_Style(is)
        END IF
        CALL Get_Atom_Info(is)
        CALL Get_Bond_Info(is)
        CALL Get_Angle_Info(is)
        CALL Get_Dihedral_Info(is)
        CALL Get_Improper_Info(is)
     
        IF(nfragments(is) /= 0  ) THEN
           species_list(is)%fragment = .TRUE.
           CALL Get_Fragment_Info(is)
           IF (int_sim_type /= sim_mcf) THEN
              CALL Get_Fragment_File_Info(is)
           END IF
           IF (nfragments(is) > 1 .AND. int_sim_type /= sim_mcf) THEN
              CALL Get_Fragment_Connectivity_Info(is)
           END IF
        END IF

        CLOSE(unit=molfile_unit)

 ENDDO

  
  ! Load coordinates of all the fragment conformations if 
  ! there is at least one fragment

  DO is = 1, nspecies
     IF (species_list(is)%fragment) THEN
        IF ( .NOT. (sim_type == 'NVT_MC_Fragment' .OR. sim_type == 'NVT_MC_Ring_Fragment')) THEN
           IF ( int_sim_type /= sim_mcf ) THEN
              CALL Get_Fragment_Coords
              EXIT
           END IF
        END IF
     END IF
  END DO

END SUBROUTINE Get_Molecule_Info
!********************************************************************************

!****************************************************************************

!********************************************************************************
SUBROUTINE Get_L_Coul_CBMC(is)
!********************************************************************************
! This routine opens the molfile and determines the whether to
! include Coulombic interactions duing the biased growth
! Input file name format is 
! # Key_Word
! Variables associated with Key_Word
! Added by NR 01-05-2010
!********************************************************************************

  INTEGER, INTENT(IN) :: is

  INTEGER :: ierr,line_nbr,nbr_entries
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(molfile_unit)

  ierr = 0
  line_nbr = 0

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(molfile_unit,line_string,ierr)

     IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading species type."
        CALL Clean_Abort(err_msg,'Get_L_Coul_CBMC')
     END IF

     IF (line_string(1:14) == '# L_Coul_CBMC') THEN
        line_nbr = line_nbr + 1

        CALL Parse_String(molfile_unit,line_nbr,1,nbr_entries,line_array,ierr)

! Assign the first entry on the line to the type of species
        species_list(is)%L_Coul_CBMC = String_To_Logical(line_array(1))
        WRITE(logunit,'(A,I4, 2x,L4 )') 'L_Coul_CBMC given for species ',is, &
             species_list(is)%L_Coul_CBMC

        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

! L_Coul_CBMC not specified. Default is true
        species_list(is)%L_Coul_CBMC = .true.
        WRITE(logunit,'(A,I4)') 'L_Coul_CBMC not given for species ',is
        WRITE(logunit,*) 'Defaulting to true'
        WRITE(logunit,*) 'default ',species_list(is)%L_Coul_CBMC

        EXIT

     ENDIF

  ENDDO

  WRITE(logunit,*)

END SUBROUTINE Get_L_Coul_CBMC

!********************************************************************************
SUBROUTINE Get_Species_Type(is)
!********************************************************************************
! This routine opens the molfile and determines the type of species.
! Input file name format is 
! # Key_Word
! Variables associated with Key_Word
!********************************************************************************

  INTEGER, INTENT(IN) :: is

  INTEGER :: ierr,line_nbr,nbr_entries, is_1
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(molfile_unit)

  ierr = 0
  line_nbr = 0

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(molfile_unit,line_string,ierr)

     IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading species type."
        CALL Clean_Abort(err_msg,'Get_Species_Type')
     END IF

     IF (line_string(1:14) == '# Species_Type') THEN
        line_nbr = line_nbr + 1
        
        CALL Parse_String(molfile_unit,line_nbr,1,nbr_entries,line_array,ierr)

! Assign the first entry on the line to the type of species
        species_list(is)%species_type = line_array(1)
        IF(species_list(is)%species_type == 'SORBATE') THEN 
           species_list(is)%int_species_type = int_sorbate
        ELSE
           species_list(is)%int_species_type = int_solvent
        END IF
        WRITE(logunit,'(A,I4, 2x, A)') 'Species type given for species ',is, &
             species_list(is)%species_type
  
        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

! No species type specified. Default to normal
        species_list(is)%species_type = 'normal'
        WRITE(logunit,'(A,I4)') 'No species type given for species ',is
        WRITE(logunit,*) 'Defaulting to normal'
        WRITE(logunit,*) 'assigned to type ',species_list(is)%species_type

        EXIT

     ENDIF

  ENDDO


END SUBROUTINE Get_Species_Type

!********************************************************************************
SUBROUTINE Get_Insertion_Style(is)
!********************************************************************************
! This routine opens the molfile and determines the type of species.
! Input file name format is 
! # Key_Word
! Variables associated with Key_Word
!********************************************************************************

  INTEGER, INTENT(IN) :: is

  INTEGER :: ierr,line_nbr,nbr_entries, is_1
  CHARACTER(120) :: line_string, line_array(30)

  REWIND(molfile_unit)

  ierr = 0
  line_nbr = 0

! These are default choices

  species_list(is)%lcom = .false.
  species_list(is)%insert_style = "NONE"

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(molfile_unit,line_string,ierr)

     IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading insertion style."
        CALL Clean_Abort(err_msg,'Get_Insertion_Style')
     END IF

     IF (line_string(1:24) == '# 1st_Fragment_Ins_Style') THEN
        line_nbr = line_nbr + 1
        
        CALL Parse_String(molfile_unit,line_nbr,1,nbr_entries,line_array,ierr)

! Assign the first entry on the line to the type of species
        species_list(is)%insert_style = line_array(1)
        IF((species_list(is)%insert_style == 'COM').or.(species_list(is)%insert_style == 'com')) THEN 
           species_list(is)%lcom = .true.
        ELSE IF((species_list(is)%insert_style == 'BEAD').or.(species_list(is)%insert_style == 'bead')) THEN
           species_list(is)%lcom = .false.
        ELSE
            WRITE(logunit,'(A,I4)') 'Insertion Style not provided for species', is
            err_msg = ""
            err_msg(1) = "The input for Insertion Style Not understood by the program"
            err_msg(2) = "The possible options are : COM or BEAD"
            CALL Clean_Abort(err_msg,'Get_Insertion_Style')
        END IF
        WRITE(logunit,'(A,I4, 2x, A)') 'Insertion style given for species ',is, &
             species_list(is)%insert_style
  
        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

        ! 1st Fragment style is not specified. For now, default it to 
        ! COM

        species_list(is)%insert_style = 'COM'
        species_list(is)%lcom = .true.

        WRITE(logunit,'(A48,2X,I3)') '1st_Fragment_Ins_Style not specified for species', is
        WRITE(logunit,'(A)') 'using default of "COM" for this species.'
        

!        WRITE(logunit,'(A,I4)') 'Insertion Style not given for species ',is
!        err_msg = ""
!        err_msg(1) = "Please include '# 1St_Fragment_Ins_Style' in the mcf file"
!        err_msg(2) = "Possible options are: COM or BEAD"
!        CALL Clean_Abort(err_msg, 'Get_Insertion_Style')

        EXIT

     ENDIF

  ENDDO


END SUBROUTINE Get_Insertion_Style 

!********************************************************************************
SUBROUTINE Get_Atom_Info(is)
!********************************************************************************
! This routine opens the molfile and loads information into the nonbond_list
! for species is.
!********************************************************************************

  INTEGER, INTENT(IN) :: is

  INTEGER :: ierr,line_nbr,nbr_entries, ia
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(molfile_unit)

  ierr = 0
  line_nbr = 0

  nring_atoms(is) = 0
  nexo_atoms(is) = 0
  
  ring_atom_ids(:,is) = 0
  exo_atom_ids(:,is) = 0

  IF (.NOT. ALLOCATED (has_charge)) THEN
     ALLOCATE(has_charge(nspecies))
     has_charge(:) = .FALSE.
  END IF

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(molfile_unit,line_string,ierr)

     IF (ierr /= 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading atom info."
        CALL Clean_Abort(err_msg,'Get_Atom_Info')
     END IF

     IF (line_string(1:11) == '# Atom_Info') THEN
        line_nbr = line_nbr + 1
        
        CALL Parse_String(molfile_unit,line_nbr,1,nbr_entries,line_array,ierr)
        ! Make sure the number of atoms is still the same. It should not have changed!
        IF (String_To_Int(line_array(1)) /= natoms(is)) THEN
           err_msg = ""
           err_msg(1) = 'natoms is inconsistent'
           CALL Clean_Abort(err_msg,'Get_Atom_Info')
        ENDIF

        species_list(is)%molecular_weight = 0.0_DP
        
        DO ia = 1,natoms(is)
           ! Now read the entries on the next lines. There must be at least 8 for 
           ! each atom.
           CALL Parse_String(molfile_unit,line_nbr,8,nbr_entries,line_array,ierr)

           ! Test for problems readin file
           IF (ierr /= 0) THEN
              err_msg = ""
              err_msg(1) = "Error reading atom info."
              CALL Clean_Abort(err_msg,'Get_Atom_Info')
           END IF

           ! Test to make sure atoms are listed 1,2,3,... in mcf file           
           IF (String_To_Int(line_array(1)) /= ia) THEN
              err_msg = ""
              err_msg(1) = 'atoms must be listed sequentially'
              CALL Clean_Abort(err_msg,'Get_Atom_Info')
           ENDIF
           
           ! Assign appropriate values to list elements
           nonbond_list(ia,is)%atom_name = line_array(2)
           nonbond_list(ia,is)%element = line_array(3)
           nonbond_list(ia,is)%mass = String_To_Double(line_array(4))
           nonbond_list(ia,is)%charge = String_To_Double(line_array(5))
           IF(nonbond_list(ia,is)%charge .NE. 0.0_DP) has_charge(is) = .TRUE. 
           nonbond_list(ia,is)%vdw_potential_type = line_array(6)
	   

           species_list(is)%molecular_weight = species_list(is)%molecular_weight + &
                nonbond_list(ia,is)%mass

           ! For now, force all atoms to have the same vdw pair style that was
           ! specified in the input file.  We can relax this restriction later.
           IF (nonbond_list(ia,is)%vdw_potential_type /= vdw_style(1)) THEN
              err_msg = ""
              err_msg(1) = 'vdw_potential type does not equal specified vdw_style'
              CALL Clean_Abort(err_msg,'Get_Atom_Info')
           ENDIF

           WRITE(logunit,'(A,T25,I3,1x,I3)') 'Species and atom number', is,ia
           WRITE(logunit,'(A,T25,A)') ' atom name:',nonbond_list(ia,is)%atom_name
           WRITE(logunit,'(A,T25,A)') ' element:',nonbond_list(ia,is)%element
           WRITE(logunit,'(A,T25,F10.4)') ' mass:',nonbond_list(ia,is)%mass
           WRITE(logunit,'(A,T25,F10.4)') ' charge:',nonbond_list(ia,is)%charge
           WRITE(logunit,'(A,T25,A)') ' vdw type:',nonbond_list(ia,is)%vdw_potential_type

           ! Load vdw parameters, specific for each individual type
           IF (nonbond_list(ia,is)%vdw_potential_type == 'LJ') THEN
              ! epsilon/kB in K read in
              nonbond_list(ia,is)%vdw_param(1) = String_To_Double(line_array(7))
              ! sigma = Angstrom
              nonbond_list(ia,is)%vdw_param(2) = String_To_Double(line_array(8))

              WRITE(logunit,'(A,T25,F10.4)') ' Epsilon / kB in K:', &
                   nonbond_list(ia,is)%vdw_param(1)
              WRITE(logunit,'(A,T25,F10.4)') ' Sigma in A:', &
                   nonbond_list(ia,is)%vdw_param(2)

              ! Convert epsilon to atomic units amu A^2/ps^2
              nonbond_list(ia,is)%vdw_param(1) = kboltz* nonbond_list(ia,is)%vdw_param(1) 
              ! Set number of vdw parameters
              nbr_vdw_params = 2

           ELSEIF (nonbond_list(ia,is)%vdw_potential_type == 'NONE') THEN
              WRITE(logunit,'(A,I6,1x,I6)') & 
                   'No VDW potential assigned to atom, species: ',ia,is

              ! Set number of vdw parameters
              nbr_vdw_params = 0

           ELSE
              err_msg = ""
              err_msg(1) = 'vdw_potential type improperly specified in mcf file'
              CALL Clean_Abort(err_msg,'Get_Atom_Info')
           ENDIF

           ! there may be a card for ring atoms
           nonbond_list(ia,is)%ring_atom = .FALSE.
           IF (nbr_entries == 9) THEN
              ! this atom is a ring atom
              IF (line_array(9) == 'ring') THEN
                 nring_atoms(is) = nring_atoms(is) + 1
                 ring_atom_ids(nring_atoms(is),is) = ia
                 nonbond_list(ia,is)%ring_atom = .TRUE.
                 WRITE(logunit,*) ia ,' is a ring atom'
              END IF
           ELSE
              ! this is an not a ring atom
              nexo_atoms(is) = nexo_atoms(is) + 1
              exo_atom_ids(nexo_atoms(is),is) = ia
           END IF
              

           WRITE(logunit,*)

        ENDDO
              
        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

        ! Problem reading Atom_Info
        err_msg = ""
        err_msg(1) = 'Trouble locating # Atom_Info keyword'
        CALL Clean_Abort(err_msg,'Get_Atom_Info')

        EXIT

     ENDIF

  ENDDO
  WRITE(logunit,'(A, T40, I4,A, T45, I4)') 'Total number of ring atoms in species', is, ' is', nring_atoms(is)
  WRITE(logunit,'(A, T40, I4,A, T45, I4)') 'Total number of exo atoms in species', is, ' is', nexo_atoms(is)

  WRITE(logunit,*) 'Atom ids for ring atoms', ring_atom_ids(1:nring_atoms(is),is)
  WRITE(logunit,*) 'Atom ids for exo atoms', exo_atom_ids(1:nexo_atoms(is),is)

  WRITE(logunit,*) '*** Completed assigning atom info ***'
  WRITE(logunit,*)

END SUBROUTINE Get_Atom_Info


!********************************************************************************
SUBROUTINE Get_Bond_Info(is)
!********************************************************************************
  ! This routine opens the molfile and loads information into the bond_list
  ! for species is.
  ! Written by: E. Maginn
  ! Date: Sept, 2007
  ! Revision history:
  ! Mon Nov 26 06:35:20 MST 2007: Added section that sets fixed bond lengths and
  !                               reads in fixed bond parameter.
!********************************************************************************

  INTEGER, INTENT(IN) :: is

  INTEGER :: ierr,line_nbr,nbr_entries, ib
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(molfile_unit)

  ierr = 0
  line_nbr = 0

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(molfile_unit,line_string,ierr)

     IF (ierr /= 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading atom info."
        CALL Clean_Abort(err_msg,'Get_Bond_Info')
     END IF

     IF (line_string(1:11) == '# Bond_Info') THEN
        line_nbr = line_nbr + 1
        
        line_array = ""
        CALL Parse_String(molfile_unit,line_nbr,1,nbr_entries,line_array,ierr)
        ! Make sure the number of bonds is still the same. It should not have changed!
        IF (String_To_Int(line_array(1)) /= nbonds(is)) THEN
           err_msg = ""
           err_msg(1) = 'nbonds is inconsistent'
           CALL Clean_Abort(err_msg,'Get_Bond_Info')
        ENDIF
        
        IF (nbonds(is) == 0) THEN
           WRITE(logunit,*) 'No bonds in species ',is
           EXIT
        ENDIF

        DO ib = 1,nbonds(is)
           ! Now read the entries on the next lines. There must be at least 4 for 
           ! each bond.
           line_array = ""
           CALL Parse_String(molfile_unit,line_nbr,4,nbr_entries,line_array,ierr)

           ! Test for problems readin file
           IF (ierr /= 0) THEN
              err_msg = ""
              err_msg(1) = "Error reading bond info."
              CALL Clean_Abort(err_msg,'Get_Bond_Info')
           END IF

           ! Test to make sure bonds are listed 1,2,3,... in mcf file           
           IF (String_To_Int(line_array(1)) /= ib) THEN
              err_msg = ""
              err_msg(1) = 'bonds must be listed sequentially'
              CALL Clean_Abort(err_msg,'Get_Bond_Info')
           ENDIF
           
           ! Assign appropriate values to list elements
           bond_list(ib,is)%atom1 = String_To_Int(line_array(2))
           bond_list(ib,is)%atom2 = String_To_Int(line_array(3))
           bond_list(ib,is)%bond_potential_type = line_array(4)

           WRITE(logunit,'(A,T25,I3,1x,I3)') 'Species and bond number', is,ib
           WRITE(logunit,'(A,T25,I3)') ' atom1:',bond_list(ib,is)%atom1
           WRITE(logunit,'(A,T25,I3)') ' atom2:',bond_list(ib,is)%atom2
           WRITE(logunit,'(A,T25,A)') ' bond type:',bond_list(ib,is)%bond_potential_type

           ! Load bond potential parameters, specific for each individual type
           IF (bond_list(ib,is)%bond_potential_type == 'fixed') THEN
              bond_list(ib,is)%int_bond_type = int_none    
              WRITE(logunit,'(A,I6,1x,I6, A, I4)') & 
                   'Bond fixed between atoms: ',bond_list(ib,is)%atom1, bond_list(ib,is)%atom2, &
                   ' in species', is
              ! Fixed bond length in A
              bond_list(ib,is)%bond_param(1) = String_To_Double(line_array(5))
              WRITE(logunit,'(A,T25,F10.4)') 'Fixed bond length, in A:',bond_list(ib,is)%bond_param(1)

              ! Set number of bond parameters
              nbr_bond_params = 1

           ELSE
              err_msg = ""
              err_msg(1) = 'bond_potential type improperly specified in mcf file'
              CALL Clean_Abort(err_msg,'Get_Bond_Info')
           ENDIF

           WRITE(logunit,*)

        ENDDO
              
        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

        ! Problem reading Bond_Info
        err_msg = ""
        err_msg(1) = 'Trouble locating # Bond_Info keyword'
        CALL Clean_Abort(err_msg,'Get_Bond_Info')

        EXIT

     ENDIF

  ENDDO

  WRITE(logunit,*) '*** Completed assigning bond info ***'
  WRITE(logunit,*)

END SUBROUTINE Get_Bond_Info

!********************************************************************************
SUBROUTINE Get_Angle_Info(is)
!********************************************************************************
  ! This routine opens the molfile and loads information into the angle_list
  ! for species is.
  
  ! Written by: E. Maginn
  ! Date: Sept, 2007
  ! Revision history:
  ! Mon Nov 26 06:35:20 MST 2007: Added section that sets fixed bond angles and
  !                               reads in fixed angle parameter.
!********************************************************************************

  INTEGER, INTENT(IN) :: is

  INTEGER :: ierr,line_nbr,nbr_entries, iang, nangles_linear
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(molfile_unit)

  ierr = 0
  line_nbr = 0

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(molfile_unit,line_string,ierr)

     IF (ierr /= 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading angle info."
        CALL Clean_Abort(err_msg,'Get_Angle_Info')
     END IF

     IF (line_string(1:12) == '# Angle_Info') THEN
        line_nbr = line_nbr + 1
        species_list(is)%linear = .FALSE. 
        line_array = ""
        CALL Parse_String(molfile_unit,line_nbr,1,nbr_entries,line_array,ierr)
        ! Make sure the number of angles is still the same. It should not have changed!
        IF (String_To_Int(line_array(1)) /= nangles(is)) THEN
           err_msg = ""
           err_msg(1) = 'nangles is inconsistent'
           CALL Clean_Abort(err_msg,'Get_Angle_Info')
        ENDIF
        
        IF (nangles(is) == 0) THEN
           WRITE(logunit,*) 'No angles in species ',is
           EXIT
        ENDIF

        DO iang = 1,nangles(is)
           ! Now read the entries on the next lines. There must be at least 5 for 
           ! each angle.
           line_array = ""
           CALL Parse_String(molfile_unit,line_nbr,5,nbr_entries,line_array,ierr)

           ! Test for problems readin file
           IF (ierr /= 0) THEN
              err_msg = ""
              err_msg(1) = "Error reading angle info."
              CALL Clean_Abort(err_msg,'Get_Angle_Info')
           END IF

           ! Test to make sure angles are listed 1,2,3,... in mcf file           
           IF (String_To_Int(line_array(1)) /= iang) THEN
              err_msg = ""
              err_msg(1) = 'angles must be listed sequentially'
              CALL Clean_Abort(err_msg,'Get_Angle_Info')
           ENDIF
           
           ! Assign appropriate values to list elements
           angle_list(iang,is)%atom1 = String_To_Int(line_array(2))
           angle_list(iang,is)%atom2 = String_To_Int(line_array(3))
           angle_list(iang,is)%atom3 = String_To_Int(line_array(4))
           angle_list(iang,is)%angle_potential_type = line_array(5)

           WRITE(logunit,'(A,T25,I3,1x,I3)') 'Species and angle number', is,iang
           WRITE(logunit,'(A,T25,I3)') ' atom1:',angle_list(iang,is)%atom1
           WRITE(logunit,'(A,T25,I3)') ' atom2:',angle_list(iang,is)%atom2
           WRITE(logunit,'(A,T25,I3)') ' atom3:',angle_list(iang,is)%atom3
           WRITE(logunit,'(A,T25,A)') ' angle type:',angle_list(iang,is)%angle_potential_type

           ! Load angle potential parameters, specific for each individual type
           IF (angle_list(iang,is)%angle_potential_type == 'harmonic') THEN

              IF(species_list(is)%int_species_type == int_sorbate) zig_calc(is) = .TRUE.
              angle_list(iang,is)%int_angle_type = int_harmonic
              ! K_bond/kB in K/A^2 read in
              angle_list(iang,is)%angle_param(1) = String_To_Double(line_array(6))
              ! theta0 in degrees
              angle_list(iang,is)%angle_param(2) = String_To_Double(line_array(7))

              WRITE(logunit,'(A,T25,F10.4)') ' K_angle in K/rad^2:', &
                   angle_list(iang,is)%angle_param(1)
              WRITE(logunit,'(A,T25,F10.4)') ' theta0 in degrees:', &
                   angle_list(iang,is)%angle_param(2)

              ! Convert force constant to atomic units amu A^2/(rad^2 ps^2) 
              ! so that EANGLE = amu A^2/ps^2
              angle_list(iang,is)%angle_param(1) = kboltz * angle_list(iang,is)%angle_param(1)

              ! Convert the nominal bond angle to radians
              angle_list(iang,is)%angle_param(2) = (PI/180.0_DP)*angle_list(iang,is)%angle_param(2)

              ! Set number of angle parameters
              nbr_angle_params = 2
              species_list(is)%linear = .FALSE.

           ELSEIF (angle_list(iang,is)%angle_potential_type == 'fixed') THEN
              angle_list(iang,is)%int_angle_type = int_none    
              
              angle_list(iang,is)%angle_param(1) = String_To_Double(line_array(6))

              WRITE(logunit,'(A,I6,1x,I6, 1x,I6,A, I4)') & 
                   'Angle fixed between atoms: ',angle_list(iang,is)%atom1, angle_list(iang,is)%atom2, &
                   angle_list(iang,is)%atom3,'in species', is
              WRITE(logunit,'(A,T25,F10.4)') ' fixed bond angle in degrees:', &
                   angle_list(iang,is)%angle_param(1)

              ! Set number of angle parameter = 1 for the fixed DOF
              nbr_angle_params = 1
              IF(angle_list(iang,is)%angle_param(1) .NE. 180.0_DP) species_list(is)%linear = .FALSE.

           ELSE
              err_msg = ""
              err_msg(1) = 'angle_potential type improperly specified in mcf file'
              CALL Clean_Abort(err_msg,'Get_Angle_Info')
           ENDIF

           WRITE(logunit,*)

        ENDDO
              
        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

        ! Problem reading Angle_Info
        err_msg = ""
        err_msg(1) = 'Trouble locating # Angle_Info keyword'
        CALL Clean_Abort(err_msg,'Get_Angle_Info')

        EXIT

     ENDIF

  ENDDO

  WRITE(logunit,*)
  
  ! Now loop over all the angles and the species is linear if
  ! 1. All angles are fixed and the nominal value for each of
  ! the angles is 180.0_DP

  nangles_linear = 0
  DO iang = 1, nangles(is)
     
     IF ( angle_list(iang,is)%int_angle_type == int_none) THEN
        IF ( angle_list(iang,is)%angle_param(1) == 180.0_DP) THEN
           nangles_linear = nangles_linear + 1
        END IF
     END IF
  END DO

  IF (nangles_linear == nangles(is)) species_list(is)%linear = .TRUE.

  IF(species_list(is)%linear) THEN
     WRITE(logunit,'(A,2x,I3,2x,A)') 'Molecule',is,'is defined to be linear'
  ELSE
     WRITE(logunit,'(A,2x,I3,2x,A)') 'Molecule',is,'is NOT defined to be linear'
  END IF
  WRITE(logunit,*)
  WRITE(logunit,*) '*** Completed assigning angle info ***'
  WRITE(logunit,*)

END SUBROUTINE Get_Angle_Info



!********************************************************************************
SUBROUTINE Get_Dihedral_Info(is)
!********************************************************************************
! This routine opens the molfile and loads information into the dihedral_list
! for species is.
! Modified by Dr. Amir Vahid on 12/8/2012 to include for multiple dihedrals
!********************************************************************************

  INTEGER, INTENT(IN) :: is

  INTEGER :: ierr,line_nbr,nbr_entries, idihed
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(molfile_unit)

  ierr = 0
  line_nbr = 0


  DO
     line_nbr = line_nbr + 1

     CALL Read_String(molfile_unit,line_string,ierr)

     IF (ierr /= 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading dihedral info."
        CALL Clean_Abort(err_msg,'Get_Dihedral_Info')
     END IF

     IF (line_string(1:15) == '# Dihedral_Info') THEN
        line_nbr = line_nbr + 1
        
        line_array = ""
        CALL Parse_String(molfile_unit,line_nbr,1,nbr_entries,line_array,ierr)
        ! Make sure the number of angles is still the same. It should not have changed!
        IF (String_To_Int(line_array(1)) /= ndihedrals(is)) THEN
           err_msg = ""
           err_msg(1) = 'ndihedrals is inconsistent'
           CALL Clean_Abort(err_msg,'Get_Dihedral_Info')
        ENDIF
        
        IF (ndihedrals(is) == 0) THEN
           WRITE(logunit,*) 'No dihedrals in species ',is
           EXIT
        ENDIF

        DO idihed = 1,ndihedrals(is)
           ! Now read the entries on the next lines. There must be at least 6 for 
           ! each dihedral.
           line_array = ""
           CALL Parse_String(molfile_unit,line_nbr,6,nbr_entries,line_array,ierr)

           ! Test for problems readin file
           IF (ierr /= 0) THEN
              err_msg = ""
              err_msg(1) = "Error reading dihedral info."
              CALL Clean_Abort(err_msg,'Get_Dihedral_Info')
           END IF

           ! Test to make sure angles are listed 1,2,3,... in mcf file           
           IF (String_To_Int(line_array(1)) /= idihed) THEN
              err_msg = ""
              err_msg(1) = 'dihedrals must be listed sequentially'
              CALL Clean_Abort(err_msg,'Get_Dihedral_Info')
           ENDIF
           
           ! Assign appropriate values to list elements
           dihedral_list(idihed,is)%atom1 = String_To_Int(line_array(2))
           dihedral_list(idihed,is)%atom2 = String_To_Int(line_array(3))
           dihedral_list(idihed,is)%atom3 = String_To_Int(line_array(4))
           dihedral_list(idihed,is)%atom4 = String_To_Int(line_array(5))

           dihedral_list(idihed,is)%dihedral_potential_type = line_array(6)

           WRITE(logunit,'(A,T25,I3,1x,I3)') 'Species and dihedral number', is,idihed
           WRITE(logunit,'(A,T25,I3)') ' atom1:',dihedral_list(idihed,is)%atom1
           WRITE(logunit,'(A,T25,I3)') ' atom2:',dihedral_list(idihed,is)%atom2
           WRITE(logunit,'(A,T25,I3)') ' atom3:',dihedral_list(idihed,is)%atom3
           WRITE(logunit,'(A,T25,I3)') ' atom4:',dihedral_list(idihed,is)%atom4
           WRITE(logunit,'(A,T25,A)') ' dihedral type:', &
                dihedral_list(idihed,is)%dihedral_potential_type

           ! Load dihedral potential parameters, specific for each individual type
           IF (dihedral_list(idihed,is)%dihedral_potential_type == 'OPLS') THEN

              IF(species_list(is)%int_species_type == int_sorbate) zig_calc = .TRUE.
              dihedral_list(idihed,is)%int_dipot_type = int_opls
              !a0, a1, a2, a3 in kJ/mol
              dihedral_list(idihed,is)%dihedral_param(1) = String_To_Double(line_array(7))
              dihedral_list(idihed,is)%dihedral_param(2) = String_To_Double(line_array(8))
              dihedral_list(idihed,is)%dihedral_param(3) = String_To_Double(line_array(9))
              dihedral_list(idihed,is)%dihedral_param(4) = String_To_Double(line_array(10))

              WRITE(logunit,'(A,T25,F10.4)') ' a0, kJ/mol:', &
                   dihedral_list(idihed,is)%dihedral_param(1)
              WRITE(logunit,'(A,T25,F10.4)') ' a1, kJ/mol:', &
                   dihedral_list(idihed,is)%dihedral_param(2)
              WRITE(logunit,'(A,T25,F10.4)') ' a2, kJ/mol:', &
                   dihedral_list(idihed,is)%dihedral_param(3)
              WRITE(logunit,'(A,T25,F10.4)') ' a3, kJ/mol:', &
                   dihedral_list(idihed,is)%dihedral_param(4)

              ! Convert to molecular units amu A^2/ps^2
              dihedral_list(idihed,is)%dihedral_param(1) = kjmol_to_atomic* dihedral_list(idihed,is)%dihedral_param(1)
              dihedral_list(idihed,is)%dihedral_param(2) = kjmol_to_atomic* dihedral_list(idihed,is)%dihedral_param(2)
              dihedral_list(idihed,is)%dihedral_param(3) = kjmol_to_atomic* dihedral_list(idihed,is)%dihedral_param(3)
              dihedral_list(idihed,is)%dihedral_param(4) = kjmol_to_atomic* dihedral_list(idihed,is)%dihedral_param(4)

              ! Set number of dihedral parameters
              nbr_dihedral_params = 4

              
           ELSE IF (dihedral_list(idihed,is)%dihedral_potential_type == 'CHARMM') THEN
              IF(species_list(is)%int_species_type == int_sorbate) zig_calc = .TRUE.
              dihedral_list(idihed,is)%int_dipot_type = int_charmm
              dihedral_list(idihed,is)%dihedral_param(1) = String_To_Double(line_array(7))
              dihedral_list(idihed,is)%dihedral_param(2) = String_To_Double(line_array(8))
              dihedral_list(idihed,is)%dihedral_param(3) = String_To_Double(line_array(9))
              !
              WRITE(logunit,'(A,T25,F10.4)') ' a0, kJ/mol:', &
                   dihedral_list(idihed,is)%dihedral_param(1)
              WRITE(logunit,'(A,T25,F10.4)') ' n ', &
                   dihedral_list(idihed,is)%dihedral_param(2)
              WRITE(logunit,'(A,T25,F10.4)') 'delta', &
                   dihedral_list(idihed,is)%dihedral_param(3)
              
              
              ! Convert to molecular units amu A^2/ps^2 and the delta
              ! parameter to radians
              dihedral_list(idihed,is)%dihedral_param(1) = kjmol_to_atomic* dihedral_list(idihed,is)%dihedral_param(1)
              dihedral_list(idihed,is)%dihedral_param(3) = (PI/180.0_DP)* dihedral_list(idihed,is)%dihedral_param(3)
              
              nbr_dihedral_params = 3
			  
!AV: AMBER style for dihedral multiplicity, cf. Zhong et al. JpcB, 115, 10027, 2011.
!Note that I assumed the maximum # of dihedral multiplicity is 3 and tried to avoide a 2-dimensional arrays.
           ELSE IF (dihedral_list(idihed,is)%dihedral_potential_type == 'AMBER') THEN
              IF(species_list(is)%int_species_type == int_sorbate) zig_calc = .TRUE.
              dihedral_list(idihed,is)%int_dipot_type = int_amber
              dihedral_list(idihed,is)%dihedral_param(1) = String_To_Double(line_array(7))
              dihedral_list(idihed,is)%dihedral_param(2) = String_To_Double(line_array(8))
              dihedral_list(idihed,is)%dihedral_param(3) = String_To_Double(line_array(9))
              dihedral_list(idihed,is)%dihedral_param(4) = String_To_Double(line_array(10))
              dihedral_list(idihed,is)%dihedral_param(5) = String_To_Double(line_array(11))
              dihedral_list(idihed,is)%dihedral_param(6) = String_To_Double(line_array(12))
              dihedral_list(idihed,is)%dihedral_param(7) = String_To_Double(line_array(13))
              dihedral_list(idihed,is)%dihedral_param(8) = String_To_Double(line_array(14))
              dihedral_list(idihed,is)%dihedral_param(9) = String_To_Double(line_array(15))
			  !AV: commented out b/c 3 terms is usually enough.
			  !dihedral_list(idihed,is)%dihedral_param(10) = String_To_Double(line_array(16))
			  !dihedral_list(idihed,is)%dihedral_param(11) = String_To_Double(line_array(17))
			  !dihedral_list(idihed,is)%dihedral_param(12) = String_To_Double(line_array(18))
			  
              !
              WRITE(logunit,'(A,T25,F10.4)') ' a01, kJ/mol:', &
                   dihedral_list(idihed,is)%dihedral_param(1)
              WRITE(logunit,'(A,T25,F10.4)') ' n1 ', &
                   dihedral_list(idihed,is)%dihedral_param(2)
              WRITE(logunit,'(A,T25,F10.4)') 'delta1', &
                   dihedral_list(idihed,is)%dihedral_param(3)
              WRITE(logunit,'(A,T25,F10.4)') ' a02, kJ/mol:', &
                   dihedral_list(idihed,is)%dihedral_param(4)
              WRITE(logunit,'(A,T25,F10.4)') ' n2 ', &
                   dihedral_list(idihed,is)%dihedral_param(5)
              WRITE(logunit,'(A,T25,F10.4)') 'delta2', &
                   dihedral_list(idihed,is)%dihedral_param(6)
              WRITE(logunit,'(A,T25,F10.4)') ' a03, kJ/mol:', &
                   dihedral_list(idihed,is)%dihedral_param(7)
              WRITE(logunit,'(A,T25,F10.4)') ' n3 ', &
                   dihedral_list(idihed,is)%dihedral_param(8)
              WRITE(logunit,'(A,T25,F10.4)') 'delta3', &
                   dihedral_list(idihed,is)%dihedral_param(9)
              
              
              ! Convert to molecular units amu A^2/ps^2 and the delta
              ! parameter to radians
              dihedral_list(idihed,is)%dihedral_param(1) = kjmol_to_atomic* dihedral_list(idihed,is)%dihedral_param(1)
			  dihedral_list(idihed,is)%dihedral_param(4) = kjmol_to_atomic* dihedral_list(idihed,is)%dihedral_param(4)
			  dihedral_list(idihed,is)%dihedral_param(7) = kjmol_to_atomic* dihedral_list(idihed,is)%dihedral_param(7)
              dihedral_list(idihed,is)%dihedral_param(3) = (PI/180.0_DP)* dihedral_list(idihed,is)%dihedral_param(3)
			  dihedral_list(idihed,is)%dihedral_param(6) = (PI/180.0_DP)* dihedral_list(idihed,is)%dihedral_param(6)
			  dihedral_list(idihed,is)%dihedral_param(9) = (PI/180.0_DP)* dihedral_list(idihed,is)%dihedral_param(9)
              
              nbr_dihedral_params = 9
			  

           ELSE IF (dihedral_list(idihed,is)%dihedral_potential_type == 'harmonic') THEN
              IF(species_list(is)%int_species_type == int_sorbate) zig_calc = .TRUE.
              dihedral_list(idihed,is)%int_dipot_type = int_harmonic
              ! d0 read in in units ofin K/radians^2 read in
              dihedral_list(idihed,is)%dihedral_param(1) = String_To_Double(line_array(7))
              ! theta0 in degrees
              dihedral_list(idihed,is)%dihedral_param(2) = String_To_Double(line_array(8))

              WRITE(logunit,'(A,T25,F10.4)') ' Do_angle in K/rad^2:', &
                   dihedral_list(idihed,is)%dihedral_param(1)
              WRITE(logunit,'(A,T25,F10.4)') ' theta0 in degrees:', &
                   dihedral_list(idihed,is)%dihedral_param(2)

              ! Convert force constant to atomic units amu A^2/(rad^2 ps^2) 
              ! so that Edihedral = amu A^2/ps^2v
              dihedral_list(idihed,is)%dihedral_param(1) = kboltz * dihedral_list(idihed,is)%dihedral_param(1)

              ! Convert the nominal bond angle to radians
              dihedral_list(idihed,is)%dihedral_param(2) = (PI/180.0_DP)*dihedral_list(idihed,is)%dihedral_param(2)

              ! Set number of angle parameters
              nbr_dihedral_params = 2


           ELSEIF (dihedral_list(idihed,is)%dihedral_potential_type == 'none') THEN
              dihedral_list(idihed,is)%int_dipot_type = int_none
              WRITE(logunit,'(A,4(I6,1x),A,I4)') & 
                   'No dihedral potential between atoms: ',&
                   dihedral_list(idihed,is)%atom1, dihedral_list(idihed,is)%atom2, &
                   dihedral_list(idihed,is)%atom3, dihedral_list(idihed,is)%atom4, &
                   'in species', is

              ! Set number of dihedral parameters
              nbr_dihedral_params = 0

           ELSE
              err_msg = ""
              err_msg(1) = 'dihedral_potential type improperly specified in mcf file'
              CALL Clean_Abort(err_msg,'Get_Dihedral_Info')
           ENDIF

           WRITE(logunit,*)

        ENDDO
              
        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

        ! Problem reading Dihedral_Info
        err_msg = ""
        err_msg(1) = 'Troube locating # Dihedral_Info keyword'
        CALL Clean_Abort(err_msg,'Get_Dihedral_Info')

        EXIT

     ENDIF

  ENDDO

  WRITE(logunit,*) '*** Completed assigning dihedral info ***'
  WRITE(logunit,*)

END SUBROUTINE Get_Dihedral_Info
!********************************************************************************

!********************************************************************************
SUBROUTINE Get_Improper_Info(is)
!********************************************************************************

INTEGER, INTENT(IN) :: is

  INTEGER :: ierr,line_nbr,nbr_entries, iimprop
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(molfile_unit)

  ierr = 0
  line_nbr = 0
!********************************************************************************

 DO
     line_nbr = line_nbr + 1

     CALL Read_String(molfile_unit,line_string,ierr)

     IF (ierr /= 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading improper info."
        CALL Clean_Abort(err_msg,'Get_Improper_Info')
     END IF

     IF (line_string(1:15) == '# Improper_Info') THEN
        line_nbr = line_nbr + 1
        
        line_array = ""
        CALL Parse_String(molfile_unit,line_nbr,1,nbr_entries,line_array,ierr)
        ! Make sure the number of dihedrals is still the same. It should not have changed!
        IF (String_To_Int(line_array(1)) /= nimpropers(is)) THEN
           err_msg = ""
           err_msg(1) = 'nimpropers is inconsistent'
           CALL Clean_Abort(err_msg,'Get_Imprper_Info')
        ENDIF
        
        IF (nimpropers(is) == 0) THEN
           WRITE(logunit,*) 'No impropers in species ',is
           EXIT
        ENDIF

        DO iimprop = 1,nimpropers(is)
           ! Now read the entries on the next lines. There must be at least 6 for 
           ! each improper.
           line_array = ""
           CALL Parse_String(molfile_unit,line_nbr,6,nbr_entries,line_array,ierr)

           ! Test for problems readin file
           IF (ierr /= 0) THEN
              err_msg = ""
              err_msg(1) = "Error reading improper info."
              CALL Clean_Abort(err_msg,'Get_Improper_Info')
           END IF

           ! Test to make sure angles are listed 1,2,3,... in mcf file           
           IF (String_To_Int(line_array(1)) /= iimprop) THEN
              err_msg = ""
              err_msg(1) = 'impropers must be listed sequentially'
              CALL Clean_Abort(err_msg,'Get_Improper_Info')
           ENDIF
           
           ! Assign appropriate values to list elements
           improper_list(iimprop,is)%atom1 = String_To_Int(line_array(2))
           improper_list(iimprop,is)%atom2 = String_To_Int(line_array(3))
           improper_list(iimprop,is)%atom3 = String_To_Int(line_array(4))
           improper_list(iimprop,is)%atom4 = String_To_Int(line_array(5))

           improper_list(iimprop,is)%improper_potential_type = line_array(6)

           WRITE(logunit,'(A,T25,I3,1x,I3)') 'Species and improper number', is,iimprop
           WRITE(logunit,'(A,T25,I3)') ' atom1:',improper_list(iimprop,is)%atom1
           WRITE(logunit,'(A,T25,I3)') ' atom2:',improper_list(iimprop,is)%atom2
           WRITE(logunit,'(A,T25,I3)') ' atom3:',improper_list(iimprop,is)%atom3
           WRITE(logunit,'(A,T25,I3)') ' atom4:',improper_list(iimprop,is)%atom4
           WRITE(logunit,'(A,T25,A)') ' dihedral type:', &
                improper_list(iimprop,is)%improper_potential_type

           ! Load improper potential parameters, specific for each individual type
           IF (improper_list(iimprop,is)%improper_potential_type == 'harmonic') THEN
              improper_list(iimprop,is)%int_improp_type = int_harmonic
              ! K_imp in K/rad^2l
              ! Function: V_imp = K_imp * (phi-phi0)^2

              ! Param 1 is k_imp and param2 is phi0
              improper_list(iimprop,is)%improper_param(1) = String_To_Double(line_array(7))
              improper_list(iimprop,is)%improper_param(2) = String_To_Double(line_array(8))

              WRITE(logunit,'(A,T25,F10.4)') ' K_improper, K/rad^2', &
                   improper_list(iimprop,is)%improper_param(1)

              ! Convert to molecular units of energy
              improper_list(iimprop,is)%improper_param(1) = kboltz * improper_list(iimprop,is)%improper_param(1)

              ! Convert phi0 to radians
              improper_list(iimprop,is)%improper_param(2) = (PI/180.0_DP) * improper_list(iimprop,is)%improper_param(2)

              ! Set number of improper parameters
              nbr_improper_params = 2
           ELSEIF (improper_list(iimprop,is)%improper_potential_type == 'cvff') THEN
              improper_list(iimprop,is)%int_improp_type = int_cvff
              ! Function: V_imp = K_imp * (1 + d * cos [n * phi])
              ! param 1 = K_imp, param 2 = d, and param 3 = n

              improper_list(iimprop,is)%improper_param(1) = String_To_Double(line_array(7))
              improper_list(iimprop,is)%improper_param(2) = String_To_Double(line_array(8))
              improper_list(iimprop,is)%improper_param(3) = String_To_Double(line_array(9))
              
              WRITE(logunit,'(A,T25,F10.4)') ' K_improper, KJ/mol', &
                   improper_list(iimprop,is)%improper_param(1)

              WRITE(logunit,'(A,T25,F10.4)') ' d_improper', &
                   improper_list(iimprop,is)%improper_param(2)

              WRITE(logunit,'(A,T25,F10.4)') ' n_improper', &
                   improper_list(iimprop,is)%improper_param(3)

              ! Convert to molecular units of energy
              improper_list(iimprop,is)%improper_param(1) = kjmol_to_atomic * improper_list(iimprop,is)%improper_param(1)

           ELSEIF (improper_list(iimprop,is)%improper_potential_type == 'none') THEN
              improper_list(iimprop,is)%int_improp_type = int_none    
              WRITE(logunit,'(A,4(I6,1x),A,I4)') & 
                   'No improper potential between atoms: ',&
                   improper_list(iimprop,is)%atom1, improper_list(iimprop,is)%atom2, &
                   improper_list(iimprop,is)%atom3, improper_list(iimprop,is)%atom4, &
                   'in species', is

              ! Set number of improper parameters
              nbr_improper_params = 0

           ELSE
              err_msg = ""
              err_msg(1) = 'improper_potential type improperly specified in mcf file'
              CALL Clean_Abort(err_msg,'Get_Improper_Info')
           ENDIF

           WRITE(logunit,*)

        ENDDO
              
        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

        ! Problem reading Improper_Info
        err_msg = ""
        err_msg(1) = 'Troube locating # Improper_Info keyword'
        CALL Clean_Abort(err_msg,'Get_Improper_Info')

        EXIT

     ENDIF

  ENDDO

  WRITE(logunit,*) '*** Completed assigning improper info ***'
  WRITE(logunit,*)



END SUBROUTINE Get_Improper_Info
!
!*******************************************************************************
SUBROUTINE Get_Fragment_Anchor_Info(is)
!*******************************************************************************
!
! This routine determines number of anchors in a fragment. Note that a
! ring fragment may have more than one anchor
!
! Written by Jindal Shah on 04/16/09
! 
!*******************************************************************************

  INTEGER, INTENT(IN) :: is

  INTEGER :: i, line_nbr, ierr, min_entries, nbr_entries, ianchor

  CHARACTER(120) :: line_String,line_array(20)


  line_nbr = 0
  ierr = 0

  REWIND(molfile_unit)

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(molfile_unit,line_string,ierr)
     IF ( ierr /= 0 ) THEN
        err_msg = ""
        err_msg(1) = 'Error encountered while reading Anchor Info'
        CALL Clean_Abort(err_msg,'Get_Fragment_Anchor_Info')
     END IF

     IF (line_string(1:13)  == '# Anchor_Info') THEN
        ! we found Anchor Info section

        DO i = 1, nfragments(is)

           line_nbr = line_nbr + 1
           
           CALL Parse_String(molfile_unit,line_nbr,2,nbr_entries,line_array,ierr)
           write(*,*) i, String_To_Int(line_array(1))
           IF ( String_To_Int(line_array(1)) /= i ) THEN
              ! fragments are not listed sequentially
              err_msg = ''
              err_msg(1) = 'Fragments must be listed sequentially'
              CALL Clean_Abort(err_msg,'Get_Fragment_Anchor_Info')
           END IF

           ! Now we read in the information on number of anchors and anchor atoms
           ! Note that the second entry contains number of anchors followed by
           ! anchor ids

           frag_list(i,is)%nanchors = String_To_Int(line_array(2))

           ALLOCATE(frag_list(i,is)%anchor(frag_list(i,is)%nanchors))

           min_entries = 2 + frag_list(i,is)%nanchors

           backspace(molfile_unit)

           CALL Parse_String(molfile_unit,line_nbr,min_entries,nbr_entries,line_array,ierr)

           DO ianchor = 1, frag_list(i,is)%nanchors

              frag_list(i,is)%anchor(ianchor) = String_To_Int(line_array(2+ianchor))

           END DO

           ! output information to the log file.

           WRITE(logunit,*)
           WRITE(logunit,*) 'Number of anchors for fragment ', i, '  is', frag_list(i,is)%nanchors
           WRITE(logunit,*) 'Anchor ids are', frag_list(i,is)%anchor(:)

        END DO

        EXIT

     ELSE IF (line_String(1:3) == 'END' .OR. line_nbr > 10000 ) THEN

        err_msg = ''
        err_msg(1) = 'Trouble locating Anchor Info section'
        CALL Clean_Abort(err_msg,'Get_Anchor_Info')

     END IF

  END DO

  WRITE(logunit,*)
  WRITE(logunit,*) '***** Finished loading Anchor Info **********'


END SUBROUTINE Get_Fragment_Anchor_Info
!*******************************************************************************

!********************************************************************************
SUBROUTINE Get_Fragment_Info(is)
!********************************************************************************
! The subroutine goes through the molecular connectivity file and determines
! total number of atoms in a fragment, their atomic ids, number of connections
! that a fragment makes and corresponding fragment ids to which the fragment is 
! connected. The end result is the frag_list arrays gets populated
!
! Written by Jindal Shah on 07/10/08
!
! Modified by Jindal Shah on 08/12/11 to automatically assign number of anchors 
! and anchor ids in a given fragment. So the following array gets populated
!
! frag_list(ifrag,is)%nanchors
! frag_list(ifrag,is)%anchor(1:nanchors)
!
!**********************************************************************************

  INTEGER :: is, line_nbr, ierr, nbr_entries, ifrag, min_entries, iatom
  INTEGER :: nanchors, iatoms, jatoms, ibonds, iatoms_bond
  INTEGER :: i_atom, j_atom, atom1, atom2
  INTEGER, ALLOCATABLE :: anchor_id(:)
  CHARACTER(120) :: line_string, line_array(20)

  line_nbr = 0
  ierr = 0
  REWIND(molfile_unit)

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(molfile_unit,line_string,ierr)
     IF (ierr /=0 ) THEN
        err_msg = ''
        err_msg(1) = 'Error reading mol file'
        CALL Clean_Abort(err_msg,'Get_Fragment_Info')
     END IF
    
     IF (line_string(1:15) == '# Fragment_Info') THEN
        ! we found a section on the fragment information
        line_nbr = line_nbr + 1
        CALL Parse_String(molfile_unit,line_nbr,1,nbr_entries,line_array,ierr)
        ! Make sure that number of fragments has not changed
        IF (String_To_Int(line_array(1)) /= nfragments(is)) THEN
           err_msg = ''
           err_msg(1) = 'Number of fragments is inconsistent'
           CALL Clean_Abort(err_msg,'Get_Fragment_Info')
        END IF

        IF (nfragments(is) == 0 ) THEN
           WRITE(logunit,*) 'No fragments in the species', is
           EXIT
        END IF

        ! Now read in the information for each of the fragments, number of atoms
        ! in the fragment, anchor and atom ids
        
        DO ifrag = 1, nfragments(is)
           
           line_nbr = line_nbr + 1
           ! We will first determine number of atoms in the current fragment and then
           ! allocate the array frag_list(i,j)%atoms and also read in the anchor
           CALL Parse_String(molfile_unit,line_nbr, 2, nbr_entries, line_array,ierr)
           
           IF ( ifrag /= String_To_Int(line_array(1))) THEN
              err_msg = ''
              err_msg = 'Fragments must be listed sequentially'
              CALL Clean_Abort(err_msg,'Get_Fragment_Info')
           END IF

           frag_list(ifrag,is)%natoms = String_To_Int(line_array(2))

           ! read in the identity of atoms
           backspace(molfile_unit)
           
           min_entries = 2 + frag_list(ifrag,is)%natoms

           CALL Parse_String(molfile_unit,line_nbr,min_entries,nbr_entries,line_array,ierr)

           IF (ierr /= 0) THEN
              err_msg = ''
              err_msg(1) = 'Number of atoms inconsistent in fragment'
              err_msg(2) = Int_To_String(ifrag)
              CALL Clean_Abort(err_msg,'Get_Fragment_Info')
           END IF

           ALLOCATE(frag_list(ifrag,is)%atoms(frag_list(ifrag,is)%natoms))

           DO iatom = 1, frag_list(ifrag,is)%natoms
              frag_list(ifrag,is)%atoms(iatom) = String_To_Int(line_array(iatom+2))
           END DO
          
           WRITE(logunit,*)
           WRITE(logunit,'(A38,1x,I4,A4,I4)') 'Total number of atoms in the fragment ', ifrag, 'is', &
                frag_list(ifrag,is)%natoms
           WRITE(logunit,'(A27)') 'Identity of these atoms are:'
           WRITE(logunit,*)  frag_list(ifrag,is)%atoms
          

        END DO

        EXIT

     ELSE IF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN
        
        ! Problem reading Improper_Info
        err_msg = ""
        err_msg(1) = 'Trouble locating # Fragment_Info keyword'
        CALL Clean_Abort(err_msg,'Get_Fragment_Info')
        
        EXIT
             
     END IF

  END DO

  WRITE(logunit,*) 
  WRITE(logunit,*) '******* Finished loading fragment info ********'
  WRITE(logunit,*)
  WRITE(logunit,*) '******* Generating anchor info ****************'

!  IF (nfragments(is) == 1) THEN

!     nanchors = 1
!     frag_list(1,is)%nanchors = nanchors
!     ALLOCATE(frag_list(1,is)%anchor(nanchors))
!     frag_list(1,is)%anchor(1) = 1

!  ELSE

     DO ifrag = 1, nfragments(is)

        ALLOCATE(anchor_id(frag_list(ifrag,is)%natoms))

        nanchors = 0
        anchor_id(:) = 0

        DO iatoms = 1, frag_list(ifrag,is)%natoms


           i_atom = frag_list(ifrag,is)%atoms(iatoms)
           iatoms_bond = 0
           
           DO jatoms = 1, frag_list(ifrag,is)%natoms

              ! loop over all the bonds in the species to
              ! figure out if 'iatoms' and 'jatoms' are bonded

              j_atom = frag_list(ifrag,is)%atoms(jatoms)

              IF (i_atom == j_atom) CYCLE

              DO ibonds = 1, nbonds(is)

                 atom1 = bond_list(ibonds,is)%atom1
                 atom2 = bond_list(ibonds,is)%atom2
                 
                 IF (i_atom == atom1 .AND. j_atom == atom2) &
                      iatoms_bond = iatoms_bond + 1
                 IF (i_atom == atom2 .AND. j_atom == atom1) &
                      iatoms_bond = iatoms_bond + 1

              END DO

           END DO

           IF (iatoms_bond >= 2 ) THEN
              ! this atom is connected to more than two bonds
              ! in the fragment and is an anchor
              nanchors = nanchors + 1
              anchor_id(nanchors) = i_atom
              
           END IF

        END DO

        ! at this point we have determined the total number of
        ! anchors and atom ids of the anchors put them in global
        ! arrays

        frag_list(ifrag,is)%nanchors = nanchors

        ALLOCATE(frag_list(ifrag,is)%anchor(nanchors))

        frag_list(ifrag,is)%anchor(1:nanchors) = anchor_id(1:nanchors)

        DEALLOCATE(anchor_id)

     END DO
     
!  END IF

  ! Output info

  frag_list(:,:)%ring = .FALSE.

  DO ifrag = 1,nfragments(is)
     WRITE(logunit,*)
     WRITE(logunit,'(A32,1x,I4,A4,I4)') 'Number of anchors for fragment ', ifrag, 'is', &
          frag_list(ifrag,is)%nanchors
     WRITE(logunit,'(A13,1x,I4)') 'Anchor id is:', frag_list(ifrag,is)%anchor(:)
     

     IF (frag_list(ifrag,is)%nanchors > 1) THEN
        frag_list(ifrag,is)%ring = .TRUE.
        WRITE(logunit,*)
        WRITE(logunit,*) Int_To_String(ifrag)// ' is a ring fragment'
     END IF
           
  END DO

  
END SUBROUTINE Get_Fragment_Info


!********************************************************************************
SUBROUTINE Get_Fragment_Connect_Info(is)
!********************************************************************************
! The subroutine goes through the molecular connectivity file and obtains
! information on fragment connections. It assigns the total number of connections
! for a fragment and fragmet ids of these connections. The end result is
!
! frag_list(ifrag,is)%nconnect 
! frag_list(ifrag,is)%(1:nconnect) get assigned
!
! Written by Jindal Shah on 07/11/08
!
!
!
!**********************************************************************************

  INTEGER, INTENT(IN):: is
  INTEGER :: ierr, line_nbr, nbr_entries, ifrag, min_entries, i
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(molfile_unit)
  ierr = 0
  line_nbr = 0

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(molfile_unit,line_string,ierr)
     IF ( ierr /=0 ) THEN
        err_msg = ''
        err_msg(1) = 'Error reading mol file'
        CALL Clean_Abort(err_msg,'Get_Fragment_Connect_Info')
     END IF
     
     IF (line_string(1:15) == '# Fragment_Bond') THEN

        DO ifrag = 1, nfragments(is)
           
           line_nbr = line_nbr + 1
           CALL Parse_String(molfile_unit,line_nbr,2,nbr_entries,line_array,ierr)
           
           IF (String_To_Int(line_array(1)) /= ifrag) THEN
              err_msg = ''
              err_msg(1) = 'Fragment must be listed sequentially'
              CALL Clean_Abort(err_msg,'Get_Fragment_Connect_Info')
           END IF
           
           ! The second entry represents number of fragments connected to ifrag
           frag_list(ifrag,is)%nconnect = String_To_Int(line_array(2))

           ALLOCATE(frag_list(ifrag,is)%frag_connect(frag_list(ifrag,is)%nconnect))

           min_entries = 2 + frag_list(ifrag,is)%nconnect
           
           backspace(molfile_unit)

           ! Assign ids of all the fragments connected to ifrag
           CALL Parse_String(molfile_unit,line_nbr,min_entries,nbr_entries,line_array,ierr)

           DO i = 1, frag_list(ifrag,is)%nconnect
              frag_list(ifrag,is)%frag_connect(i) = String_To_Int(line_array(2+i))
           END DO

           WRITE(logunit,*) 
           WRITE(logunit,'(A33,1X,I4,A4,I4)') 'Number of connections of fragment', ifrag, 'is', &
                frag_list(ifrag,is)%nconnect
           WRITE(logunit,'(A20)') 'These fragments are:'
           WRITE(logunit,*) frag_list(ifrag,is)%frag_connect


        END DO

        EXIT

     ELSE IF(line_string(1:3) == 'END' .or. line_nbr > 10000) THEN
        
        ! Problem reading Improper_Info
        err_msg = ""
        err_msg(1) = 'Trouble locating # Fragment_Bond keyword'
        CALL Clean_Abort(err_msg,'Get_Fragment_Connect_Info')
        
        EXIT
             
     END IF
     

  END DO

  WRITE(logunit,*) 
  WRITE(logunit,*) '**** Finished loading fragment connectivity *******'

END SUBROUTINE Get_Fragment_Connect_Info
!********************************************************************************

SUBROUTINE Get_Fragment_Connectivity_Info(is)
  !******************************************************************************
  ! This subroutine obtains the information of bonds between various fragments
  ! 
  ! Written by Jindal Shah on 05/25/08
  ! This routine is very similar to Get_Bond_Info
  !
  ! The routine was modified on 08/12/11 to automatically calculate
  ! 
  ! frag_list(ifrag,is)%nconnect
  ! frag_list(ifrag,is)%frag_connect(1:nconnect)
  !******************************************************************************

  INTEGER, INTENT(IN) :: is

  INTEGER :: ierr, line_nbr, ifrag, nbr_entries, j, ifrag_connect, frag1, frag2
  INTEGER, ALLOCATABLE :: temp_frag(:)

  CHARACTER(120) :: line_string,line_array(20)

  ierr = 0
  line_nbr = 0

  REWIND(molfile_unit)

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(molfile_unit, line_string, ierr)

     IF (ierr /= 0) THEN
        err_msg = ''
        err_msg(1) = 'Error while reading Fragment_Connectivity information'
        CALL Clean_Abort(err_msg,'Get_Fragment_Connectity_Info')
     END IF

     IF (line_string(1:23) == '# Fragment_Connectivity') THEN
        ! Encountered the part where fragment bond information is specified
        line_nbr = line_nbr + 1

        CALL Parse_String(molfile_unit,line_nbr,1,nbr_entries,line_array,ierr)

        ! Make sure the number of fragment connections is still the same.

        IF (String_To_Int(line_array(1)) /= fragment_bonds(is)) THEN
           err_msg = ""
           err_msg(1) = 'Fragment bonds is inconsistent for species ' // TRIM(Int_To_String(is))
           CALL Clean_Abort(err_msg,'Get_Fragment_Connectivity')
        END IF

        ! Now loop over all the fragment bonds and extract information on connections

        DO ifrag = 1, fragment_bonds(is)
           line_nbr = line_nbr + 1
           line_array = ""
           CALL Parse_String(molfile_unit,line_nbr,3,nbr_entries,line_array,ierr)
           
           IF (ierr /= 0 ) THEN
              err_msg = ""
              err_msg(1) = 'Error reading Fragment_Connectivity information'
              CALL Clean_Abort(err_msg,'Get_Fragment_Connectivity_Info')
           END IF

           ! Test to make sure that fragment bonds are listed as 1, 2, 3 .... in mcf file

           IF (String_To_Int(line_array(1)) /= ifrag) THEN
              err_msg = ""
              err_msg(1) = 'Fragment bonds must be listed sequentially'
              CALL Clean_Abort(err_msg,'Get_Fragment_Connectivity')
           END IF

           ! Assign appropriate values to the list elements

           fragment_bond_list(ifrag,is)%fragment1 = String_To_Int(line_array(2))
           fragment_bond_list(ifrag,is)%fragment2 = String_To_Int(line_array(3))

           WRITE(logunit,'(A32,1X,I3,1x,I3)') 'Species and fragment bond number', is,ifrag
           WRITE(logunit,'(A,T25,I3)') ' fragment 1:',fragment_bond_list(ifrag,is)%fragment1
           WRITE(logunit,'(A,T25,I3)') ' fragment 2:',fragment_bond_list(ifrag,is)%fragment2

        END DO


        EXIT

     ELSE IF (line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN

        ! problem reading Fragment Bond Connectivity

        err_msg = ''
        err_msg(1) = 'Section on Fragment_Connectivity is missing'
        CALL Clean_Abort(err_msg,'Get_Fragment_Connectivity_Info')
        
     END IF
               

  END DO


  ! Let us use this information to generate information about number of fragment
  ! connections for each fragment and ids of all the fragments that are connected
  ! to a given fragment.
  !
  ! after this the following array will be populated
  !
  ! frag_list(ifrag,is)%nconnect
  ! frag_list(ifrag,is)%frag_connect(1:nconnect)

  ALLOCATE(temp_frag(nfragments(is)))


  DO ifrag = 1, nfragments(is)

     ifrag_connect = 0

     temp_frag(:) = 0

     DO j = 1, fragment_bonds(is)

        frag1 = fragment_bond_list(j,is)%fragment1
        frag2 = fragment_bond_list(j,is)%fragment2

        IF ( ifrag == frag1 .OR. ifrag == frag2 ) THEN
           
           ifrag_connect = ifrag_connect + 1

           IF ( ifrag == frag1 ) temp_frag(ifrag_connect) = frag2
           IF ( ifrag == frag2 ) temp_frag(ifrag_connect) = frag1

        END IF

     END DO

     ! Now populate the frag_list array

     frag_list(ifrag,is)%nconnect = ifrag_connect

     IF (ifrag_connect /=0 ) THEN

        ALLOCATE(frag_list(ifrag,is)%frag_connect(ifrag_connect))

        DO j = 1, ifrag_connect

           frag_list(ifrag,is)%frag_connect(j) = temp_frag(j)

        END DO

     END IF

     WRITE(logunit,*) 
     WRITE(logunit,'(A34,1X,I4,A4,I4)') 'Number of connections of fragment', ifrag, 'is', &
          frag_list(ifrag,is)%nconnect
     WRITE(logunit,'(A21)') 'These fragments are:'
     WRITE(logunit,*) frag_list(ifrag,is)%frag_connect

  END DO

  DEALLOCATE(temp_frag)
  
 
END SUBROUTINE Get_Fragment_Connectivity_Info
!******************************************************************************


SUBROUTINE Get_Fragment_File_Info(is)
!*******************************************************************************
! This routine reads in the names of files where reservoir libraries are stored
!
! Written by Jindal Shah on 09/22/08
!
!********************************************************************************

  INTEGER :: ierr, line_nbr, i, j, ifrag, nbr_entries, is
  REAL(DP) :: vdw_cutoff, coul_cutoff
  CHARACTER(120) :: line_string, line_array(20)
  CHARACTER(4) :: ring_flag
  
  REWIND(inputunit)

  ierr = 0
  line_nbr = 0

! declare that all the fragments are not ring fragments
! if it is then, it will be assigned a true flag later on

  frag_list(:,:)%rcut_vdwsq = 0.0_DP
  frag_list(:,:)%rcut_coulsq = 0.0_DP
  frag_list(:,:)%alpha_ewald = 0.0_DP

  DO

     line_nbr = line_nbr + 1
     
     CALL Read_String(inputunit,line_string,ierr)

     IF (line_string(1: ) == '# Fragment_Files' ) THEN
        ! we found a section on fragment file information
        ! on each line of input we have name of the file corresponding to 
        ! various fragments

        ! skip the first few  lines for the reservoir file for other species

        DO i = 1, is - 1
           line_nbr = line_nbr + nfragments(i)
           DO j = 1, nfragments(i)
              CALL Read_String(inputunit,line_string,ierr)
           END DO
        END DO

        DO ifrag = 1, nfragments(is)

           line_nbr = line_nbr + 1

           CALL Parse_String(inputunit,line_nbr,2,nbr_entries,line_array,ierr)

           res_file(ifrag,is) = TRIM(line_array(1))
           ! assign a fragmet type 
           frag_list(ifrag,is)%type = String_To_Int(line_array(2))

           WRITE(logunit,*) 'Fragment file for fragment ', TRIM(Int_To_String(ifrag)), ' is'
           WRITE(logunit,*) res_file(ifrag,is)
           WRITE(logunit,*) 'Fragment type', frag_list(ifrag,is)%type

           IF (nbr_entries > 2 .AND. nbr_entries /= 6) THEN

              err_msg = ""
              err_msg(1) = 'More than two entries found for'
              err_msg(2) = 'fragment '//Int_To_String(ifrag)
              err_msg(3) = 'species '//Int_To_String(is)
              err_msg(4) = 'But does not appear that all the parameters for ring'
              err_msg(5) = 'are specified'

              CALL Clean_Abort(err_msg,'Get_Fragment_File_Info')

           END IF

!!$           IF ( frag_list(ifrag,is)%ring .AND. nbr_entries /= 6) THEN
!!$
!!$              err_msg= ''
!!$              err_msg(1) = 'The fragment '//Int_To_String(ifrag)
!!$              err_msg(2) = 'in the species '//Int_To_String(is)
!!$              err_msg(3) = 'appears to be a ring fragment'
!!$              err_msg(4) = 'However not all parameters are specified in the # Fragment_Files section'
!!$              err_msg(5) = 'But it appears that not all the parameters for ring are specified'
!!$
!!$              CALL Clean_Abort(err_msg, 'Get_Fragment_File_Info')
!!$
!!$           END IF

           IF (nbr_entries == 6) THEN
              ! there are three more entries corresponding to "ring", rcut_vdw, rcut_ewald
              ! and alpha_ewald

              ring_flag = TRIM(line_array(3))

              IF ( ring_flag /= "ring") THEN

                 err_msg = ""
                 err_msg(1) = "Cannot determine whether it's a ring fragment"
                 err_msg(2) = 'fragment '//Int_To_String(ifrag)
                 err_msg(3) = 'species '//Int_to_string(is)

                 CALL Clean_Abort(err_msg,'Get_Fragment_File_Info')

              ELSE

                 vdw_cutoff = String_To_Double(line_array(4))
                 coul_cutoff  = String_To_Double(line_array(5))
                 frag_list(ifrag,is)%rcut_vdwsq = vdw_cutoff * vdw_cutoff
                 frag_list(ifrag,is)%rcut_coulsq = coul_cutoff * coul_cutoff
                 frag_list(ifrag,is)%alpha_ewald = String_To_Double(line_array(6))

                 WRITE(logunit,*)
                 WRITE(logunit,*) 'Fragment ',TRIM(Int_To_String(ifrag)), ' of species ', TRIM(Int_To_String(is))
                 WRITE(logunit,*) 'is a ring fragment.'
                 WRITE(logunit,*)
                 WRITE(logunit,*) 'Parameters used for generating the fragment conformations'
                 WRITE(logunit,*)
                 WRITE(logunit,*) 'VDW cutoff', vdw_cutoff
                 WRITE(logunit,*) 'Ewald cutoff', coul_cutoff
                 WRITE(logunit,*) 'Ewald alpha parameter', frag_list(ifrag,is)%alpha_ewald
                 WRITE(logunit,*)

              END IF

           END IF

        END DO

        EXIT

     ELSE IF (line_string(1:3) == 'END' .OR. line_nbr > 10000) THEN

        err_msg = ''
        err_msg(1) = 'Fragment file info section is missing' 
        CALL Clean_Abort(err_msg,'Get_Fragment_File_Info')

     END IF
        

  END DO


END SUBROUTINE Get_Fragment_File_Info
!**********************************************************************************
SUBROUTINE Get_Fragment_Coords
  ! This subroutine loads in x,y,z coordinates of the fragments by fragment type
  ! 
  ! Called by:
  !
  ! Get_Molecule_Info
  !
  ! First written by Jindal Shah on 03/03/09
  !
  !
  !*********************************************************************************

  INTEGER :: nfrag_types, natoms_max, max_config, is, ifrag, ifrag_type, this_config
  INTEGER :: iconfig, ia, this_atom

  REAL(DP) :: x_this, y_this, z_this
  REAL(DP) :: this_temperature, this_nrg

  CHARACTER :: symbol*1

  LOGICAL, ALLOCATABLE :: config_read(:)

  ! Allocate arrays for frag_coords
  
  nfrag_types = MAXVAL(frag_list(:,:)%type)
  natoms_max = MAXVAL(frag_list(:,:)%natoms)
  
  ! Determine maximum number of configurations
  
  ALLOCATE(config_read(nfrag_types))
  config_read(:) = .FALSE.

  max_config = 0

  DO is = 1, nspecies
     IF (nfragments(is) /=0 ) THEN

        DO ifrag = 1, nfragments(is)
           
           ifrag_type = frag_list(ifrag,is)%type

           ! open the file and read # of configurations
           OPEN(UNIT=10,FILE=res_file(ifrag,is))
           READ(10,*) this_config

           frag_list(ifrag,is)%nconfig = this_config
           max_config = MAX(this_config,max_config)

           
           CLOSE(UNIT=10)

        END DO
     END IF
  END DO
  
  WRITE(logunit,*) 
  WRITE(logunit,*) 'Maximum configurations stored', max_config
  WRITE(logunit,*)

  ALLOCATE(frag_coords(natoms_max,max_config,nfrag_types),STAT=Allocatestatus)


  IF (Allocatestatus /= 0 ) THEN
     err_msg = ''
     err_msg(1) = 'Error allocating array for frag_coords'
     CALL Clean_Abort(err_msg,'Get_Frag_Coords')
  END IF

  ALLOCATE(nrg_frag(max_config,nfrag_types), STAT = AllocateStatus)

  IF (Allocatestatus /= 0 ) THEN
     err_msg = ''
     err_msg(1) = 'Error allocating array for nrg_frag'
     CALL Clean_Abort(err_msg,'Get_Frag_Coords')
  END IF
  

  ! Load coordinates
  
  frag_coords(:,:,:)%rxp = 0.0_DP
  frag_coords(:,:,:)%ryp = 0.0_DP
  frag_coords(:,:,:)%rzp = 0.0_DP

  config_read(:) = .FALSE.

  DO is = 1, nspecies
     
     IF(nfragments(is) /=0 ) THEN
        
        DO ifrag = 1, nfragments(is)
           
           ifrag_type = frag_list(ifrag,is)%type
           
           IF (config_read(ifrag_type)) CYCLE

           ! open the file and read # of configurations
           OPEN(UNIT=10,FILE=res_file(ifrag,is))
           READ(10,*) this_config
           
           DO iconfig = 1, this_config

              ! read in the energy of the fragment
              READ(10,*) this_temperature, this_nrg
              nrg_frag(iconfig,ifrag_type) = this_nrg
              ! read coordinates
              DO ia = 1, frag_list(ifrag,is)%natoms

                 READ(10,*) symbol, x_this, y_this, z_this
                 
                 frag_coords(ia,iconfig,ifrag_type)%rxp = x_this
                 frag_coords(ia,iconfig,ifrag_type)%ryp = y_this
                 frag_coords(ia,iconfig,ifrag_type)%rzp = z_this

                 this_atom = frag_list(ifrag,is)%atoms(ia)

 !                WRITE(12,*) nonbond_list(this_atom,is)%element, frag_coords(ia,iconfig,ifrag_type)%rxp, &
 !                     frag_coords(ia,iconfig,ifrag_type)%ryp, frag_coords(ia,iconfig,ifrag_type)%rzp
                 
              END DO

           END DO


           WRITE(logunit,*)
           WRITE(logunit,*) 'Finished loading fragment coordinates from'
           WRITE(logunit,*) res_file(ifrag,is)
           WRITE(logunit,*)
           
           config_read(ifrag_type) = .TRUE.
           
           CLOSE(UNIT=10)
           
        END DO

     END IF
           
           
  END DO

  DEALLOCATE(config_read)


END SUBROUTINE Get_Fragment_Coords
!********************************************************************************

SUBROUTINE Get_Intra_Scaling
!********************************************************************************
  INTEGER :: ierr,line_nbr,nbr_entries, iimprop, is
  CHARACTER(120) :: line_string, line_array(20)
  LOGICAL :: intrascaling_set

!********************************************************************************
  REWIND(inputunit)

  ierr = 0
  line_nbr = 0
  intrascaling_set = .false.
  ALLOCATE(scale_1_2_vdw(nspecies));ALLOCATE(scale_1_3_vdw(nspecies))
  ALLOCATE(scale_1_4_vdw(nspecies));ALLOCATE(scale_1_N_vdw(nspecies))
  ALLOCATE(scale_1_2_charge(nspecies));ALLOCATE(scale_1_3_charge(nspecies))
  ALLOCATE(scale_1_4_charge(nspecies));ALLOCATE(scale_1_N_charge(nspecies))

 DO
     line_nbr = line_nbr + 1

     CALL Read_String(inputunit,line_string,ierr)

     IF (ierr /= 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading ibtrascaling info."
        CALL Clean_Abort(err_msg,'Get_Intra_Scaling')
     END IF

     IF (line_string(1:15) == '# Intra_Scaling') THEN
        line_nbr = line_nbr + 1
        DO is = 1, nspecies
           IF (int_vdw_style(1) /= vdw_none) THEN
              ! Read vdw scaling which is listed first
              line_array = ""
              CALL Parse_String(inputunit,line_nbr,4,nbr_entries,line_array,ierr)
           
              ! Test for problems reading file
              IF (ierr /= 0) THEN
                 err_msg = ""
                 err_msg(1) = "Error reading Intra_Scaling info."
                 CALL Clean_Abort(err_msg,'Get_Intra_Scaling_Info')
              END IF
              
              ! Assign the vdw scaling
              scale_1_2_vdw(is) = String_To_Double(line_array(1))
              scale_1_3_vdw(is) = String_To_Double(line_array(2))
              scale_1_4_vdw(is) = String_To_Double(line_array(3))
              scale_1_N_vdw(is) = String_To_Double(line_array(4))
           ENDIF

           IF (int_charge_style(1) /= charge_none) THEN
              ! Read coul scaling which is listed second
              line_array = ""
              CALL Parse_String(inputunit,line_nbr,4,nbr_entries,line_array,ierr)
              
              scale_1_2_charge(is) = String_To_Double(line_array(1))
              scale_1_3_charge(is) = String_To_Double(line_array(2))
              scale_1_4_charge(is) = String_To_Double(line_array(3))
              scale_1_N_charge(is) = String_To_Double(line_array(4))
           ENDIF
        END DO
        intrascaling_set = .true.
        ! exit the loop

        EXIT

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

        ! No intrascaling set explicitly - use default.
        
        scale_1_2_vdw(:) = 0.0
        scale_1_3_vdw(:) = 0.0
        scale_1_4_vdw(:) = 0.5
        scale_1_N_vdw(:) = 1.0

        scale_1_2_charge(:) = 0.0
        scale_1_3_charge(:) = 0.0
        scale_1_4_charge(:) = 0.5
        scale_1_N_charge(:) = 1.0

        EXIT

     ENDIF

  ENDDO
  ! Report to logfile what scaling is used.
  IF (intrascaling_set) THEN
     WRITE(logunit,*) 'intramolecular scaling factors explicitly set'
  ELSE
     WRITE(logunit,*) 'Using default intramolecular scaling factors, if required'
  ENDIF

  DO is = 1, nspecies
     WRITE(logunit,'(A,T50,I7)') 'Intra molecule scaling factors for species', is 
     WRITE(logunit,'(A,T30,f7.3)') 'VDW 1-2 scaling factor', scale_1_2_vdw(is)
     WRITE(logunit,'(A,T30,f7.3)') 'VDW 1-3 scaling factor', scale_1_3_vdw(is)
     WRITE(logunit,'(A,T30,f7.3)') 'VDW 1-4 scaling factor', scale_1_4_vdw(is) 
     WRITE(logunit,'(A,T30,f7.3)') 'VDW 1-N scaling factor', scale_1_N_vdw(is) 

     WRITE(logunit,'(A,T30,f7.3)') 'Coulomb 1-2 scaling factor', scale_1_2_charge(is) 
     WRITE(logunit,'(A,T30,f7.3)') 'Coulomb 1-3 scaling factor', scale_1_3_charge(is) 
     WRITE(logunit,'(A,T30,f7.3)') 'Coulomb 1-4 scaling factor', scale_1_4_charge(is) 
     WRITE(logunit,'(A,T30,f7.3)') 'Coulomb 1-N scaling factor', scale_1_N_charge(is)
     WRITE(logunit,*) 
  END DO

  WRITE(logunit,*)
  WRITE(logunit,*) '*** Completed assigning intramolecular scaling factors ***'
  WRITE(logunit,*)


END SUBROUTINE Get_Intra_Scaling
!********************************************************************************


!********************************************************************************
SUBROUTINE Get_Box_Info
!********************************************************************************
  ! This routine determines the number of boxes, the dimensions of the boxes and 
  ! box types.
  !
  ! Allowed box types:
  ! CUBIC, ORTHOGONAL, CELL_MATRIX
!********************************************************************************
  INTEGER :: ierr,line_nbr,nbr_entries,ibox, is
  CHARACTER(120) :: line_string, line_array(20)

!********************************************************************************
  REWIND(inputunit)

  ierr = 0
  line_nbr = 0
 
  DO
     line_nbr = line_nbr + 1

     CALL Read_String(inputunit,line_string,ierr)

     IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error reading box info."
        CALL Clean_Abort(err_msg,'Get_Box_Info')
     END IF

     IF (line_string(1:10) == '# Box_Info') THEN
        line_nbr = line_nbr + 1
        
        CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

        ! Number of boxes
        nbr_boxes = String_To_Int(line_array(1))

        IF ( int_sim_type == sim_gemc_ig .AND. nbr_boxes < 3 ) THEN
           err_msg = ''
           err_msg(1) = 'GEMC simulation with intermediate box is specified'
           err_msg(2) = ' Number of boxes less than 3'
           CALL Clean_Abort(err_msg,'Get_Box_Info')
        END IF

        WRITE(logunit,*) '*** Simulation box data ***'
        WRITE(logunit,'(A,T30,I3)') 'number of simulation boxes ',nbr_boxes
        WRITE(logunit,*)

        ! Allocate arrays associated with the box variables
        ALLOCATE(box_list(nbr_boxes), STAT=AllocateStatus)
        IF (AllocateStatus /= 0) THEN
           write(*,*)'memory could no tbe allocated for box_list array'
           write(*,*)'stopping'
           STOP
        END IF

        ALLOCATE(l_cubic(nbr_boxes), STAT=AllocateStatus)
        IF (AllocateStatus /= 0) THEN
           write(*,*)'memory could no tbe allocated for l_cubic array'
           write(*,*)'stopping'
           STOP
        END IF

        l_cubic(:) = .FALSE.

        DO ibox = 1,nbr_boxes
           ! Get box type
           CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
           line_nbr = line_nbr + 1
           box_list(ibox)%box_shape = line_array(1)

           IF (box_list(ibox)%box_shape == 'CUBIC') THEN
              box_list(ibox)%int_box_shape = int_cubic
              l_cubic(ibox) = .TRUE.
              ! Read in the x,y,z box edge lengths in A
              CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
              line_nbr = line_nbr + 1
              box_list(ibox)%length(1,1) = String_To_Double(line_array(1))
              box_list(ibox)%length(2,2) = String_To_Double(line_array(1))
              box_list(ibox)%length(3,3) = String_To_Double(line_array(1))

              box_list(ibox)%hlength(1,1) = 0.5_DP * box_list(ibox)%length(1,1)
              box_list(ibox)%hlength(2,2) = 0.5_DP * box_list(ibox)%length(2,2)
              box_list(ibox)%hlength(3,3) = 0.5_DP * box_list(ibox)%length(3,3)

              WRITE (logunit,'(A,T15,I3,T25,A,T40,A)') 'Box number ',ibox,'Box Shape: ', ' Cubic '
              WRITE (logunit,'(A,T20,F10.4,T35,A)') 'Each Side Of :', box_list(ibox)%length(1,1), 'Angstrom'
              WRITE(logunit,*)
              
              ! Set off-diagonal components to zero
              box_list(ibox)%length(1,2) = 0.0
              box_list(ibox)%length(1,3) = 0.0
              box_list(ibox)%length(2,1) = 0.0
              box_list(ibox)%length(2,3) = 0.0
              box_list(ibox)%length(3,1) = 0.0
              box_list(ibox)%length(3,2) = 0.0

           ELSEIF (box_list(ibox)%box_shape == 'ORTHOGONAL') THEN
              box_list(ibox)%int_box_shape = int_ortho
              CALL Parse_String(inputunit,line_nbr,3,nbr_entries,line_array,ierr)
              line_nbr = line_nbr + 1
              box_list(ibox)%length(1,1) = String_To_Double(line_array(1))
              box_list(ibox)%length(2,2) = String_To_Double(line_array(2))
              box_list(ibox)%length(3,3) = String_To_Double(line_array(3))

              ! Set off-diagonal components to zero
              box_list(ibox)%length(1,2) = 0.0
              box_list(ibox)%length(1,3) = 0.0
              box_list(ibox)%length(2,1) = 0.0
              box_list(ibox)%length(2,3) = 0.0
              box_list(ibox)%length(3,1) = 0.0
              box_list(ibox)%length(3,2) = 0.0

              WRITE (logunit,'(A,T15,I3,T25,A,T40,A)') 'Box number ',ibox,'Box Shape: ', ' Orthorhombic '
              WRITE (logunit,'(A,T20,F10.4,T35,A)') 'X dimension :', box_list(ibox)%length(1,1), 'Angstrom'
              WRITE (logunit,'(A,T20,F10.4,T35,A)') 'Y dimension :', box_list(ibox)%length(2,2), 'Angstrom'
              WRITE (logunit,'(A,T20,F10.4,T35,A)') 'Z dimension :', box_list(ibox)%length(3,3), 'Angstrom'
              WRITE(logunit,*)

           ELSEIF (box_list(ibox)%box_shape == 'CELL_MATRIX') THEN
              box_list(ibox)%int_box_shape = int_cell
              CALL Parse_String(inputunit,line_nbr,3,nbr_entries,line_array,ierr)
              line_nbr = line_nbr + 1
              box_list(ibox)%length(1,1) = String_To_Double(line_array(1))
              box_list(ibox)%length(1,2) = String_To_Double(line_array(2))
              box_list(ibox)%length(1,3) = String_To_Double(line_array(3))

              CALL Parse_String(inputunit,line_nbr,3,nbr_entries,line_array,ierr)
              line_nbr = line_nbr + 1
              box_list(ibox)%length(2,1) = String_To_Double(line_array(1))
              box_list(ibox)%length(2,2) = String_To_Double(line_array(2))
              box_list(ibox)%length(2,3) = String_To_Double(line_array(3))

              CALL Parse_String(inputunit,line_nbr,3,nbr_entries,line_array,ierr)
              line_nbr = line_nbr + 1
              box_list(ibox)%length(3,1) = String_To_Double(line_array(1))
              box_list(ibox)%length(3,2) = String_To_Double(line_array(2))
              box_list(ibox)%length(3,3) = String_To_Double(line_array(3))

              WRITE (logunit,'(A,T15,I3,T25,A,T40,A)') 'Box number ',ibox,'Box Shape : ', ' CELL_MATRIX '
              WRITE (logunit,'(T5,f10.4,T20,f10.4,T30,f10.4)') box_list(ibox)%length(1,1:3)
              WRITE (logunit,'(T5,f10.4,T20,f10.4,T30,f10.4)') box_list(ibox)%length(2,1:3)
              WRITE (logunit,'(T5,f10.4,T20,f10.4,T30,f10.4)') box_list(ibox)%length(3,1:3)
              WRITE(logunit,*)

           ELSE
              err_msg = ""
              err_msg(1) = 'Box type not properly specified as'
              err_msg(2) = box_list(ibox)%box_shape
              err_msg(3) = 'Must be a CUBIC, ORTHOGONAL, CELL_MATRIX'
              CALL Clean_Abort(err_msg, 'Get_Box_Info')

           END IF

           ! Compute information on the simulation box and write to log. 
           CALL Compute_Cell_Dimensions(ibox)

           write(logunit,'(A,3(f10.4,3x))') 'Cell basis vector lengths in A,  ',&
                box_list(ibox)%basis_length
           write(logunit,'(A,3(f10.4,3x))') 'Cosine of angles alpha, beta, gamma ',&
                box_list(ibox)%cos_angle
           write(logunit,'(A,3(f10.4,3x))') 'Distance between box faces ',&
                box_list(ibox)%face_distance
           write(logunit,'(A,f18.4)') 'Box volume, A^3 ', box_list(ibox)%volume
           write(logunit,*)

           ! Skip 1 line between boxes
           line_nbr = line_nbr + 1
           CALL Read_String(inputunit,line_string,ierr)

        ENDDO ! End loop over nbr_boxes

        EXIT ! We found box info so we are done

     ELSEIF (line_string(1:3) == 'END' .or. line_nbr > 10000) THEN

        ! No box info specified
        err_msg = ""
        err_msg(1) = 'No box info specified in input file'
        CALL Clean_Abort(err_msg,'Get_Box_Info')

        EXIT

     ENDIF

  ENDDO !  End overall read do

  WRITE(logunit,*) '*** Finished loading box information ***'
  WRITE(logunit,*)

  ! Allocate memory for total number of mols of each species in a given box

  ALLOCATE(nmols(nspecies,nbr_boxes),Stat=Allocatestatus)
  IF (Allocatestatus /=0) THEN
     err_msg = ''
     err_msg(1) = 'Memory could not be allocated for nmols'
     CALL Clean_Abort(err_msg,'Get_Box_Info')
  END IF

  ALLOCATE(vdw_style(nbr_boxes) , charge_style(nbr_boxes))
  ALLOCATE(vdw_sum_style(nbr_boxes) , charge_sum_style(nbr_boxes))

  ALLOCATE(int_vdw_style(nbr_boxes) , int_vdw_sum_style(nbr_boxes))
  ALLOCATE(int_charge_style(nbr_boxes) , int_charge_sum_style(nbr_boxes))

  ALLOCATE(rcut_CBMC(nbr_boxes))
  ALLOCATE(rcut_vdw(nbr_boxes) , rcut_coul(nbr_boxes))
  ALLOCATE(ron_charmm(nbr_boxes) , roff_charmm(nbr_boxes))
  ALLOCATE(ron_switch(nbr_boxes) , roff_switch(nbr_boxes))
  ALLOCATE(rcut_max(nbr_boxes), rcut_vdwsq(nbr_boxes))
  ALLOCATE(ron_switch_sq(nbr_boxes) , roff_switch_sq(nbr_boxes))
  ALLOCATE(ron_charmmsq(nbr_boxes) , roff_charmmsq(nbr_boxes))
  ALLOCATE(switch_factor1(nbr_boxes) , switch_factor2(nbr_boxes))
  ALLOCATE(rcut_coulsq(nbr_boxes))
  ALLOCATE(rcut9(nbr_boxes) , rcut3(nbr_boxes))

  ALLOCATE(W_tensor_charge(3,3,nbr_boxes) , W_tensor_recip(3,3,nbr_boxes))
  ALLOCATE(W_tensor_vdw(3,3,nbr_boxes) , W_tensor_total(3,3,nbr_boxes))
  ALLOCATE(W_tensor_elec(3,3,nbr_boxes), Pressure_tensor(3,3,nbr_boxes))

  ALLOCATE(P_inst(nbr_boxes),P_ideal(nbr_boxes))

  W_tensor_charge(:,:,:) = 0.0_DP
  W_tensor_recip(:,:,:) = 0.0_DP
  W_tensor_vdw(:,:,:) = 0.0_DP
  W_tensor_total(:,:,:) = 0.0_DP
  W_tensor_elec(:,:,:) = 0.0_DP
  Pressure_tensor(:,:,:) = 0.0_DP

  P_inst(:) = 0.0_DP
  P_ideal(:) = 0.0_DP

  ALLOCATE(rcut_vdw3(nbr_boxes), rcut_vdw6(nbr_boxes))

END SUBROUTINE Get_Box_Info

! Obtain the participation matrix

SUBROUTINE Get_Temperature_Info
  ! This routine obtains temperature for all the boxes in the simulation. Make sure Get_Box_Info routine
  ! is first called before calling this routine as it requies the information on number of boxes

  IMPLICIT NONE

  INTEGER :: ierr, line_nbr, i, nbr_entries
  CHARACTER(120) :: line_string, line_array(20)

  ! Check to make sure that we have read in number of boxes if not then abort

  IF ( .NOT. ALLOCATED(box_list) ) THEN
     err_msg = ""
     err_msg(1) = 'Number of boxes has not been read yet'
     CALL Clean_Abort(err_msg,'Get_Temperature')
  END IF

  REWIND(inputunit)

  ALLOCATE(temperature(nbr_boxes))
  ALLOCATE(beta(nbr_boxes))

  ierr = 0
  line_nbr = 0

  write(logunit,*)
  write(logunit,*) '**** Simulation temperature information ****'

  outer: DO
     line_nbr = line_nbr + 1

     CALL Read_String(inputunit,line_string,ierr)

     IF ( ierr /= 0 ) THEN
        err_msg = ""
        err_msg(1) = 'Error reading temperature'
        CALL Clean_Abort(err_msg,'Get_Temperature')
     END IF

     IF(line_string(1:18) == '# Temperature_Info') THEN

        line_nbr = line_nbr + 1
        
        CALL Parse_String(inputunit,line_nbr,nbr_boxes,nbr_entries,line_array,ierr)
        
        IF ( ierr /= 0 ) THEN
           err_msg = ""
           err_msg(1) = 'Error while reading temperature info'
           CALL Clean_Abort(err_msg,'Get_Temperature_Info')
        END IF
        
        DO i = 1, nbr_boxes
           temperature(i) = String_To_Double(line_array(i))
           ! compute inverse temperature
              beta(i) = 1.0_DP / (kboltz * temperature(i))
           ! write to the logunit that temperature is specified for box
           
           write(logunit,'(A30,2X,i3,2x,A2,2X,F7.3,2X,A3)')'Temperature assigned to box ', i, 'is', temperature(i), ' K'
           write(logunit,*)
           
        END DO

        EXIT

     ELSE IF (line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN

        err_msg = ""
        err_msg(1) = 'No temperature info specified in the input file'
        CALL Clean_Abort(err_msg,'Get_Temperature_Info')
        
        EXIT outer
        
     END IF
     
  END DO outer

  write(logunit,*) '*** Finished loading information for temperature ****'
  write(logunit,*)
  
END SUBROUTINE Get_Temperature_Info
!------------------------------------------------------------------------------------------------------

SUBROUTINE Get_Pressure_Info
  ! This subroutine goes through the input file and obtains information on
  ! pressure of simulation boxes. At present, the subroutine is designed such
  ! that only different pressure can be specified.

  INTEGER :: ierr, line_nbr, nbr_entries, i
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(inputunit)

  ierr = 0
  line_nbr = 0

  IF ( .NOT. ALLOCATED(box_list)) THEN
     err_msg = ""
     err_msg = 'Box information has not yet been obtained'
     CALL Clean_Abort(err_msg,'Get_Pressure_Info')
  END IF

  ALLOCATE(pressure(nbr_boxes))

  WRITE(logunit,*) 
  WRITE(logunit,*) '************* Reading pressure Info ***********'
  

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)
     IF (ierr /=0 ) THEN
        err_msg = ''
        err_msg = 'Error reading input file.'
        CALL Clean_Abort(err_msg,'Get_Pressure_Info')
     END IF

     IF (line_string(1:15) == '# Pressure_Info') THEN

        ! pressure is specified for each of the boxes on the same line.
        line_nbr = line_nbr + 1
        CALL Parse_String(inputunit,line_nbr,nbr_boxes,nbr_entries,line_array,ierr)

        IF ( ierr /= 0 ) THEN
           err_msg = ''
           err_msg = 'Error while reading pressure info in the input file'
           CALL Clean_Abort(err_msg,'Get_Pressure_Info')
        END IF

        ! assign the pressures

        DO i = 1, nbr_boxes
           pressure(i) = String_To_Double(line_array(i))
           ! convert pressure into atomic units
           pressure(i) = pressure(i)/atomic_to_bar
           
           WRITE(logunit,*)
           WRITE(logunit,*) 'Pressure of box', i, ' is ', pressure(i),  ' amu / (A ps^2). '
           

        END DO

        EXIT
        

     ELSE IF (line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN
        err_msg = ''
        err_msg(1) = 'Pressure required but not specified in the input file'
        CALL Clean_Abort(err_msg,'Get_Pressure_Info')
        
     END IF

  END DO

  WRITE(logunit,*)
  WRITE(logunit,*) '********* Finished reading pressure info ********'
  
  
END SUBROUTINE Get_Pressure_Info
!---------------------------------------------------------------------------------------------------------

SUBROUTINE Get_Fugacity_Info
  ! this code goes through the input file and obtains information about fugacities of species
  ! for GCMC move.

  IMPLICIT NONE

  INTEGER :: line_nbr, nbr_entries, ierr, i, spec_counter, j
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(inputunit)

  ierr = 0
  line_nbr = 0
  spec_counter = 0

 ! initialize fugacities

  species_list(:)%fugacity = 0.0_DP
  species_list(:)%chem_potential = 0.0_DP
  species_list(:)%activity = 0.0_DP
  lchempot = .FALSE.
  lfugacity = .FALSE.

  inputLOOP: DO
     line_nbr = line_nbr  + 1
     CALL Read_String(inputunit,line_string,ierr)
     IF ( ierr /= 0 ) THEN
        err_msg = ''
        err_msg(1) = 'Error while reading Fugacity Info'
        CALL Clean_Abort(err_msg,'Get_Fugacity_Info')
     END IF
    
     
     IF (line_string(1:15) == '# Fugacity_Info' ) THEN
        ! we found a section that contains the information on fugacity of all the species
        lfugacity = .TRUE.
        line_nbr = line_nbr + 1
        WRITE(logunit,*)
        WRITE(logunit,*) '*********** Fugacity Info ***************'
        CALL Parse_String(inputunit,line_nbr,nspec_insert,nbr_entries,line_array,ierr)
        
        DO i = 1,nspecies

           IF(species_list(i)%int_species_type == int_sorbate) THEN
              spec_counter = spec_counter + 1
              species_list(i)%fugacity = String_To_Double(line_array(spec_counter))
              species_list(i)%fugacity = species_list(i)%fugacity / atomic_to_bar
           ELSE
              species_list(i)%fugacity = 0.0_DP
           END IF

           WRITE(logunit,*)
           WRITE(logunit,'(A,T25,I3,A,T35,E16.9,A)')'Fugacity of species', i, ' is', species_list(i)%fugacity, ' amu /(A ps^2)'
        END DO
        
        EXIT

     ELSE IF (line_string(1:25) == '# Chemical_Potential_Info'  ) THEN
        ! we found a section that contains the information on Chemical Potential of all the species
        line_nbr = line_nbr + 1
        lchempot = .TRUE.
        WRITE(logunit,*)
        WRITE(logunit,*) '*********** Fugactiy Info ***************'
        CALL Parse_String(inputunit,line_nbr,nspec_insert,nbr_entries,line_array,ierr)
        
        DO i = 1,nspecies

           ALLOCATE(species_list(i)%de_broglie(nbr_boxes))

           IF(species_list(i)%int_species_type == int_sorbate) THEN
              spec_counter = spec_counter + 1
              species_list(i)%chem_potential = String_To_Double(line_array(spec_counter))
           ELSE
              species_list(i)%chem_potential = 0.0_DP
           END IF

           ! convert the chemical potential into atomic units

           species_list(i)%chem_potential = species_list(i)%chem_potential / atomic_to_kJmol

           WRITE(logunit,*)
           WRITE(logunit,'(A,T25,I3,A,T35,E16.9,A)')'Chemical Potential of species', i, &
                ' is', species_list(i)%chem_potential, 'in atomic units'

           ! Now compute the de Broglie wavelength for this species in each box

           DO j = 1, nbr_boxes

              species_list(i)%de_broglie(j) = &
                   h_plank  * DSQRT( beta(j)/(twopi * species_list(i)%molecular_weight))

              WRITE(logunit,'(A,T35,I3,T40,A,T48,I5,T55,A)') 'de Broglie wavelength of species', i, 'in box ', j, ' is'
              WRITE(logunit,'(F14.10,2x,A)') species_list(i)%de_broglie(j), ' Angstrom'

           END DO
   
        END DO
        
        EXIT

     ELSE IF (line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN
        err_msg = ''
        err_msg(1) = 'Activity_Info section is missing from the input file'
        CALL Clean_Abort(err_msg,'Get_Fugacity_Info')
     END IF
     
  END DO inputLOOP
  
  WRITE(logunit,*)
  WRITE(logunit,*) '******* Finished reading fugacity info ***********'

END SUBROUTINE Get_Fugacity_Info


SUBROUTINE Get_Move_Probabilities
  ! This routine goes through the input file and obtains probabilities for each of the
  ! moves to be performed. At the end of the routine a check is made to ensure that
  ! all probabilities add upto 1.0_DP
  !

  IMPLICIT NONE

  INTEGER :: ierr, nbr_entries, line_nbr,i, j, ibox, is, vol_int
  INTEGER ::  kbox, this_box
  CHARACTER(120) :: line_string, line_array(30), line_string2
  CHARACTER(4) :: Symbol

  REAL(DP) :: total_mass, this_mass, prob_box_swap

  LOGICAL :: l_prob_box_swap


  ierr = 0
  line_nbr = 0

  ALLOCATE(prob_species_trans(nspecies))
  ALLOCATE(prob_species_rotate(nspecies))

  ! set all the probabilities to zero initially and make sure that they add up to 1 at the end.
  prob_trans = 0.0_DP
  prob_species_trans = 0.0_DP
  prob_rot = 0.0_DP
  prob_species_rotate = 0.0_DP
  prob_torsion = 0.0_DP
  prob_volume = 0.0_DP
  prob_angle = 0.0_DP
  prob_insertion = 0.0_DP
  prob_deletion = 0.0_DP
  prob_swap = 0.0_DP
  prob_regrowth = 0.0_DP
  prob_ring = 0.0_DP  ! sampling of ring atoms using flip move
  prob_atom_displacement = 0.0_DP ! sampling of atoms using atom displacement routine

  ALLOCATE(sorbate_file(nspecies))
  ALLOCATE(init_list(MAXVAL(natoms),1,nspecies))
  ALLOCATE(max_disp(nspecies,nbr_boxes))
  ALLOCATE(max_rot(nspecies,nbr_boxes))
  ALLOCATE(prob_rot_species(nspecies))

  REWIND(inputunit)

  ! start reading the input file

  inputLOOP: DO
     line_nbr = line_nbr + 1

     CALL Read_String(inputunit,line_string,ierr)

     IF (ierr /=0) THEN
        err_msg = ""
        err_msg(1) = 'Error while reading Move Probabilities'
        CALL Clean_Abort(err_msg,'Get_Move_Probabilities')
     END IF
     num_moves = 0
     IF (line_string(1:23) == '# Move_Probability_Info') THEN
        ! we entered a section of the input file where all the move probabilties
        ! are identified. We will read each line and then compare against possible
        ! move types that are defined.
        WRITE(logunit,*)
        WRITE(logunit,*)'********* Move Probability Info******'
        sectionLOOP: DO
           line_nbr = line_nbr + 1
           CALL Read_String(inputunit,line_string,ierr)
           IF(line_string(1:18) == '# Prob_Translation') THEN
              num_moves = num_moves + 1
              ! we found specification of translation routine
              ! parse this line to figure out what the move proability is. This line
              ! contains three fields #, Move_Translation and the third field is probability.
              line_nbr = line_nbr + 1            
              CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
              prob_trans = String_To_Double(line_array(1))

              WRITE(logunit,'(A30,2X,F9.6)') 'Translation probability is', prob_trans
              WRITE(logunit,*)

              IF (int_sim_type == sim_ring .OR. int_sim_type == sim_frag) THEN
                 ! the second line contains information on delta_cos_max and delta_phi_max
                 line_nbr = line_nbr + 1
                 CALL Parse_String(inputunit,line_nbr,2,nbr_entries,line_array,ierr)

                 delta_cos_max = String_To_Double(line_array(1))
                 delta_phi_max = String_To_Double(line_array(2))

                 WRITE(logunit,*)
                 WRITE(logunit,*) 'Maximum width in cosine of polar angle is', delta_cos_max
                 WRITE(logunit,*) 'Maximum width (degrees) in azimuthal angle is', delta_phi_max

                 ! convert delta_phi_max to radians

                 delta_phi_max = delta_phi_max * PI/180.0_DP

              ELSE
                 DO j = 1, nbr_boxes
                    line_nbr = line_nbr + 1
                    CALL Parse_String(inputunit,line_nbr,nspecies,nbr_entries,line_array,ierr)

                    ! assign the maximum displacement widths to each of the species
                    DO i = 1, nspecies
                       max_disp(i,j) = String_To_Double(line_array(i))
                       WRITE(logunit,'(A,T40,I3,A,T50,I3,T55,A,T60,F10.5)') 'Maximum displacement width for species', & 
                            i, ' in box', j, 'is', max_disp(i,j)
                       WRITE(logunit,*)
                    END DO
                 END DO

              END IF

           ELSE IF(line_string(1:15) == '# Prob_Rotation') THEN
              num_moves = num_moves + 1
              ! we found the specification for rotational probability
              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
              prob_rot = String_To_Double(line_array(1))

              WRITE(logunit,'(A30,2X,F9.6)') 'Rotation probability is', prob_rot
              WRITE(logunit,*)

              DO j = 1, nbr_boxes
                 ! get maximum rotational width for each of the species
                 line_nbr = line_nbr + 1
                 CALL Parse_String(inputunit,line_nbr,nspecies,nbr_entries,line_array,ierr)

                 DO i = 1, nspecies
                    max_rot(i,j) = String_To_Double(line_array(i))
                    ! Note that input is in degrees. Convert the displacement to radians
                    max_rot(i,j) = max_rot(i,j) * PI / 180.0_DP
                    WRITE(logunit,'(A,T40,I3,A,T50,I3,T55,A,T60,F10.4,T70,A)') 'The rotational width for the species', &
                         i, ' in box', j, ' is', max_rot(i,j), ' radians'
                    WRITE(logunit,*)
                 END DO
              END DO

           ELSE IF(line_string(1:15) == '# Prob_Dihedral') THEN
              num_moves = num_moves + 1
              ! we located specification for dihedral move probability
              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
              prob_torsion = String_To_Double(line_array(1))

              WRITE(logunit,'(A,T40,F10.4)') 'Dihedral move probability is', prob_torsion

              ! Get species dependent maximum displacements

              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,nspecies,nbr_entries,line_array,ierr)
              DO i = 1, nspecies
                 species_list(i)%max_torsion = String_To_Double(line_array(i))
                 species_list(i)%max_torsion = species_list(i)%max_torsion * PI / 180.0_DP
                 WRITE(logunit,'(A,T40,I3,A,T55,F10.4,T70,A)')'The dihedral move width for the species', i, ' is', &
                      species_list(i)%max_torsion, ' radians'
              END DO
              WRITE(logunit,*)

           ELSE IF (line_string(1:12) == '# Prob_Angle') THEN
              num_moves = num_moves + 1
              ! we located the section that describes probability to attempt angle perturbation
              line_nbr = line_nbr  + 1
              CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
              prob_angle = String_To_Double(line_array(1))

              WRITE(logunit,*)
              WRITE(logunit,'(A,T40,F10.4)') 'Angle move probability is', prob_angle

           ELSE IF (line_string(1:13) == '# Prob_Volume') THEN
              num_moves = num_moves + 1
              ! we found information for volume probability move
              ! set the flag for volume changes in log of the volume ratios to be false.

              f_dv = .true.
              f_vratio = .false.

              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

              prob_volume = String_To_Double(line_array(1))

              WRITE(logunit,*)
              WRITE(logunit,'(A40,2X,F9.6)') 'Probability for volume move is ', prob_volume

              ! Now read in information for each of the boxes for maximum displacement
              IF (int_sim_type == sim_gemc) THEN 
                 ! only one maximum volume width is specified
                 line_nbr = line_nbr + 1
                 CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

                 box_list(:)%dv_max = String_To_Double(line_array(1))

                 WRITE(logunit,*) 
                 WRITE(logunit,*) 'Maximum volume displacement for box in GEMC_NVT simulation is'
                 WRITE(logunit,"(F10.3)") box_list(1)%dv_max

              ELSE

                 DO ibox = 1,nbr_boxes

                    WRITE(logunit,*)
                    WRITE(logunit,*) 'Writing maximum volume displacement elements for box', ibox
                    WRITE(logunit,*)                 

                    line_nbr = line_nbr + 1
                    CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

                    box_list(ibox)%dv_max = String_To_Double(line_array(1))

                    WRITE(logunit,'((F24.3,2X,A))') box_list(ibox)%dv_max, ' A^3'


                 END DO

              END IF

              ! Check to see if additional line exists that indicate the type of volume
              ! move to be attempted
              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,0,nbr_entries,line_array,ierr)
              IF ( nbr_entries > 0 ) THEN
                 vol_int = String_To_Int(line_array(1))
                 IF (vol_int == 1 ) THEN
                    ! volume move to be performed in log ratio
                    f_dv = .false.
                    f_vratio = .true.
                 END IF
              END IF

              IF (f_dv) THEN

                 WRITE(logunit,*) 
                 WRITE(logunit,*) 'Volume moves will be performed in actual volumes'
              ELSE IF (f_vratio) THEN
                 WRITE(logunit,*)
                 WRITE(logunit,*) 'Volume moves will be performed in logarithm of ratio of volumes'
              END IF

           ELSE IF (line_string(1:16) == '# Prob_Insertion' .OR. line_string(1:11) == '# Prob_Swap') THEN
              ! we found information for insertion move probability or swap move probability
              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

              IF (line_string(1:16) == '# Prob_Insertion') THEN
                 num_moves = num_moves + 1
                 prob_insertion = String_To_Double(line_array(1))
                 WRITE(logunit,*)
                 WRITE(logunit,'(A,T40,F10.4)') 'Probability for insertion move is', prob_insertion
              ELSE IF (line_string(1:11) == '# Prob_Swap') THEN
                 num_moves = num_moves + 1
                 prob_swap = String_To_Double(line_array(1))
                 WRITE(logunit,*)
                 WRITE(logunit,'(A40,2X,F9.6)') 'Probability for particle swap is', prob_swap
              END IF

              ! Assign the second argument as type of insertion to sorbate species
              ! also determine the number of insertable species at this point

              nspec_insert = 0
              DO is = 1, nspecies

                 ! Also read in the file from which the configuration will be read
                 line_nbr = line_nbr + 1
                 CALL Read_String(inputunit,line_string2,ierr)

                 IF(line_string2(1:17) .NE. 'insertion method') THEN
                    err_msg(1) = "Expecting : insertion method"
                    err_msg(2) = "But found :"
                    err_msg(3) = line_string2
                    err_msg(4) = 'in '//TRIM(inputfile)
                    CALL Clean_Abort(err_msg,'Prob_Insertion')
                 END IF

                 line_nbr = line_nbr + 1
                 CALL Read_String(inputunit,line_string2,ierr)

                 IF(line_string2(1:4) == 'IGAS') THEN
                    species_list(is)%insertion = 'IGAS'
                    species_list(is)%int_insert = int_igas
                    nspec_insert = nspec_insert + 1
                    
                    IF(.NOT. ALLOCATED(n_igas)) ALLOCATE(n_igas(nspecies))
                    IF(.NOT. ALLOCATED(n_igas_update)) ALLOCATE(n_igas_update(nspecies))
                    IF(.NOT. ALLOCATED(n_igas_moves)) ALLOCATE(n_igas_moves(nspecies))
                    IF(.NOT. ALLOCATED(nzovero)) ALLOCATE(nzovero(nspecies))
                    
                    line_nbr = line_nbr + 1
                    CALL Read_String(inputunit,line_string2,ierr)

                    WRITE(logunit,*)
                    WRITE(logunit,*) 'Ideal gas configurations will be used for insertion for species ', TRIM(Int_To_String(is))

                    IF(line_string2(1:24) .NE. 'Number of igas particles') THEN
                       err_msg(1) = "Expection : Number of igas particles"
                       err_msg(2) = "But found :"
                       err_msg(3) = line_string2
                       CALL Clean_Abort(err_msg,'Prob_Insertion')
                    END IF

                    line_nbr = line_nbr + 1
                    CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
                    n_igas(is) = String_To_Int(line_array(1))

                    WRITE(logunit,*) 'Number of ideal gas particles in reservoir is',n_igas

                    line_nbr = line_nbr + 1
                    CALL Read_String(inputunit,line_string2,ierr)

                    IF(line_string2(1:15) .NE. 'Nupdate, Nmoves') THEN          
                       err_msg(1) = "Expection : Nupdate, Nmoves, Nzbyomega"
                       err_msg(2) = "But found :"
                       err_msg(3) = line_string2
                       CALL Clean_Abort(err_msg,'Prob_Insertion')
                    END IF

                    line_nbr = line_nbr + 1
                    CALL Parse_String(inputunit,line_nbr,2,nbr_entries,line_array,ierr)
                    n_igas_update(is) = String_To_Int(line_array(1))
                    n_igas_moves(is) = String_To_Int(line_array(2))
                    nzovero(is) = String_To_Int(line_array(3))

                    WRITE(logunit,*)
                    WRITE(logunit,*) 'Ideal gas reservoir will be updated every', n_igas_update, 'moves.'
                    WRITE(logunit,*) 'Each molecule will undergo',n_igas_moves, 'Monte Carlo moves.'
                    WRITE(logunit,*)

                    line_nbr = line_nbr + 1
                    CALL Read_String(inputunit,line_string2,ierr)

                    IF(line_string2(1:15) .NE. 'configuration f') THEN          
                       err_msg(1) = "Expection : configuration file"        
                       err_msg(2) = "But found :"
                       err_msg(3) = line_string2
                       CALL Clean_Abort(err_msg,'Prob_Insertion')
                    END IF

                    line_nbr = line_nbr + 1
                    CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
                    sorbate_file(is) = line_array(1)
                    
                 ELSE IF(line_string2(1:6) == 'RANDOM') THEN
                    species_list(is)%insertion = 'RANDOM'
                    nspec_insert = nspec_insert + 1

                    WRITE(logunit,*)
                    WRITE(logunit,*) ' Random insertion will be carried out for species ', TRIM(Int_To_String(is))
                    species_list(is)%int_insert = int_random

                    line_nbr = line_nbr + 1
                    CALL Read_String(inputunit,line_string2,ierr)

                    IF(line_string2(1:18) .NE. 'configuration file') THEN
                       err_msg(1) = "Expection : configuration file"      
                       err_msg(2) = "But found :"
                       err_msg(3) = line_string2
                       CALL Clean_Abort(err_msg,'Prob_Insertion')
                    END IF

                    line_nbr = line_nbr + 1
                    CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
                    sorbate_file(is) = line_array(1)

                    WRITE(logunit,*) 'Initial configuration is found in the file', TRIM(sorbate_file(is))
                    WRITE(logunit,*)

                    OPEN(unit=sorbate_unit,file=sorbate_file(is),status='old')
                    ! the file is in xyz format

                    ! obtain COM of the input geometry as well
                    total_mass = 0.0_DP

                    species_list(is)%xcom = 0.0_DP
                    species_list(is)%ycom = 0.0_DP
                    species_list(is)%zcom = 0.0_DP


                    DO j = 1, natoms(is)
                       READ(sorbate_unit,*) symbol, init_list(j,1,is)%rxp, init_list(j,1,is)%ryp, init_list(j,1,is)%rzp

                       this_mass = nonbond_list(j,is)%mass
                       total_mass = total_mass + this_mass

                       species_list(is)%xcom = species_list(is)%xcom + this_mass * &
                            init_list(j,1,is)%rxp
                       species_list(is)%ycom = species_list(is)%ycom + this_mass * &
                            init_list(j,1,is)%ryp
                       species_list(is)%zcom = species_list(is)%zcom + this_mass * &
                            init_list(j,1,is)%rzp

                    END DO

                    species_list(is)%xcom = species_list(is)%xcom / total_mass
                    species_list(is)%ycom = species_list(is)%ycom / total_mass
                    species_list(is)%zcom = species_list(is)%zcom / total_mass

                    CLOSE(UNIT=sorbate_unit)

                 ELSE IF(line_string2(1:9 ) == 'reservoir ') THEN
                    nspec_insert = nspec_insert + 1
                    species_list(is)%insertion = 'RANDOM'

                    WRITE(logunit,*) 
                    WRITE(logunit,'(A,2X,A)') 'Insertion will be carried out for species ', TRIM(Int_To_String(is))
                    WRITE(logunit,*) 'using the reservoir sampling'
                    species_list(is)%int_insert = int_random
                    species_list(is)%species_type = 'SORBATE'
                    species_list(is)%int_species_type = int_sorbate

                 ELSE IF(line_string2(1:4) == 'none') THEN
                    ! it's a species that does not get exchanged with
                    ! reservoir in the case of a GCMC simulation or
                    ! swapped in a GEMC simulation

                    species_list(is)%insertion = 'NONE'
                    species_list(is)%int_insert = int_noinsert
                    species_list(is)%species_type ='NON_EXCHANGE'
                    species_list(is)%int_species_type = int_solvent

                    WRITE(logunit,'(A50,2X,A3)') 'No insertion/swap will be carried out for species', &
                         TRIM(Int_To_String(is))
                 END IF

              END DO

              WRITE(logunit,*) 
              WRITE(logunit,"(A39,2X,I3)") 'Total number of exchangeable species is', nspec_insert
              WRITE(logunit,*)

              ! Check if individual swap species are defined for species
              l_mol_frac_swap = .TRUE.

              line_nbr = line_nbr + 1

              CALL Read_String(inputunit, line_string2,ierr)
              IF (line_string2(1:19) == '# Prob_Species_Swap') THEN
                 l_mol_frac_swap = .FALSE.
                 ! read the next line where probabilties are specified

                 ALLOCATE(prob_swap_species(nspecies))

                 line_nbr = line_nbr + 1
                 CALL Parse_String(inputunit,line_nbr,nspecies,nbr_entries,line_array,ierr)

                 DO is = 1, nspecies
                    prob_swap_species(is) = String_To_Double(line_array(is))

                    WRITE(logunit,'(A35,I4,2X,A2,2X,F10.6)') 'Cumulative probabilty of swap for species', is, 'is', &
                         prob_swap_species(is)                  
                 END DO

                 IF (ABS(prob_swap_species(nspecies) -1.0_DP) > 0.000001_DP) THEN
                    err_msg =''
                    err_msg(1) = 'Individual species swap probabilties do not add up to 1'
                    CALL Clean_Abort(err_msg,'Get_Move_Probabiltieis')
                 END IF

              ELSE

                 line_nbr = line_nbr -1
                 backspace(inputunit)

              END IF

              DO is = 1,nspecies

                 IF(species_list(is)%insertion == 'IGAS') THEN

                    ALLOCATE( molecule_list_igas(MAXVAL(nmolecules), nspecies), Stat = AllocateStatus )
                    IF (AllocateStatus /= 0) THEN
                       write(*,*)'memory could no tbe allocated for molecule_list_igas array'
                       write(*,*)'stopping'
                       STOP
                    END IF

                    ALLOCATE( atom_list_igas(MAXVAL(natoms), MAXVAL(nmolecules), nspecies), Stat = AllocateStatus )
                    IF (AllocateStatus /= 0) THEN
                       write(*,*)'memory could no tbe allocated for atom_list_igas array'
                       write(*,*)'stopping'
                       STOP
                    END IF

                    ALLOCATE( energy_igas(MAXVAL(nmolecules), nspecies), Stat = AllocateStatus )
                    IF (AllocateStatus /= 0) THEN
                       write(*,*)'memory could no tbe allocated for energy_igas array'
                       write(*,*)'stopping'
                       STOP
                    END IF

                    molecule_list_igas(:,:)%xcom = 0.0_DP
                    molecule_list_igas(:,:)%ycom = 0.0_DP
                    molecule_list_igas(:,:)%zcom = 0.0_DP
                    atom_list_igas(:,:,:)%rxp = 0.0_DP
                    atom_list_igas(:,:,:)%ryp = 0.0_DP
                    atom_list_igas(:,:,:)%rzp = 0.0_DP

                    DEALLOCATE(molecule_list)
                    DEALLOCATE(atom_list)

                    ALLOCATE(molecule_list(MAXVAL(nmolecules)+1, nspecies))
                    ALLOCATE(atom_list(MAXVAL(natoms),MAXVAL(nmolecules)+1, nspecies))

                 END IF

                 EXIT

              END DO

           ELSE IF (line_string(1:15) == '# Prob_Deletion') THEN
              num_moves = num_moves + 1
              ! we found information on the deletion move probability
              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

              prob_deletion = String_To_Double(line_array(1))

              WRITE(logunit,*)
              WRITE(logunit,'(A,T40,F10.4)') 'Probability for deletion move is', prob_deletion

!!$           ELSE IF (line_string(1:11) == '# Prob_Swap') THEN
!!$              ! we found information on the particle swap probability
!!$              line_nbr = line_nbr + 1
!!$              CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
!!$              
!!$              prob_swap = String_To_Double(line_array(1))
!!$
!!$              WRITE(logunit,*)
!!$              WRITE(logunit,'(A,T40,F10.4)')'Probability for particle swap is', prob_swap


           ELSE IF (line_string(1:15 ) == '# Prob_Regrowth') THEN
              ALLOCATE(prob_growth_species(nspecies))
              num_moves = num_moves + 1
              ! Probability for regroth of molecule is specified
              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

              prob_regrowth = String_To_Double(line_array(1))

              WRITE(logunit,*)
              WRITE(logunit,'(A30,2X,F9.6)')' Probability for regrowth is', prob_regrowth

              ! On the next line read the species probability
              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,nspecies,nbr_entries,line_array,ierr)

              DO is = 1,nspecies
                 prob_growth_species(is) = String_To_Double(line_array(is))
                 IF (is > 1) THEN
                    prob_growth_species(is) = prob_growth_species(is) + &
                         prob_growth_species(is-1)
                 END IF

                 WRITE(logunit,*)
                 WRITE(logunit,'(A50,2X,I3,2X,A2,2X,F9.6)')'Cumulative probability for regrowth of species ',is, 'is ',&
                      prob_growth_species(is) 
              END DO

              IF ( abs(prob_growth_species(nspecies) - 1.0_DP) > 0.000001_DP) THEN
                 err_msg = ''
                 err_msg(1) = 'Growth probabilities do not add up to 1.0'
                 err_msg(2) = 'Aborting'
                 CALL Clean_Abort(err_msg, 'Get_Move_Probabilities')

              END IF

           ELSE IF (line_string(1:11) == '# Prob_Ring') THEN
              num_moves = num_moves + 1
              ! Probability for ring atom displacement for fragment sampling
              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,2,nbr_entries,line_array,ierr)

              prob_ring = String_To_Double(line_array(1))

              omega_max = String_To_Double(line_array(2))

              ! convert omega_max into radians

              omega_max = omega_max * PI/180.0_DP

              WRITE(logunit,*)
              WRITE(logunit,*) 'Ring sampling enabled'
              WRITE(logunit,*) 'Probability of moving ring atoms', prob_ring
              WRITE(logunit,*) 'Maximum flip angle in radians', omega_max

           ELSE IF (line_string(1:24) == '# Prob_Atom_Displacement') THEN
              num_moves = num_moves + 1
              ! Probability for atom displacement for fragment sampling
              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

              prob_atom_displacement = String_To_Double(line_array(1))

              WRITE(logunit,*)
              WRITE(logunit,*) 'Atom displacement enabled'
              WRITE(logunit,*) 'Probability of this move', prob_atom_displacement

              ! on next line read in the information about delta_cos_max and
              ! delta_phi_max
              line_nbr = line_nbr + 1

              CALL Parse_String(inputunit,line_nbr,2,nbr_entries,line_array,ierr)

              delta_cos_max = String_To_Double(line_array(1))
              delta_phi_max = String_To_Double(line_array(2))

              WRITE(logunit,*)
              WRITE(logunit,*) 'Maximum width in cosine of polar angle is', delta_cos_max
              WRITE(logunit,*) 'Maximum width (degrees) in azimuthal angle is', delta_phi_max

              ! convert delta_phi_max to radians

              delta_phi_max = delta_phi_max * PI/180.0_DP

           ELSE IF (line_string(1:23) == '# Done_Probability_Info') THEN

              ! finished the section 

              EXIT inputLOOP

           END IF

        END DO sectionLOOP

     ELSE IF ( line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN

        err_msg = ""
        err_msg(1) = 'Move probabilities info not specified in the input'
        CALL Clean_Abort(err_msg,'Get_Move_Probabilities')

     END IF

     ! We finished the probability section

     ! Check to make sure that prob_volume is specified for an NPT simulation

  END DO inputLOOP

  IF (int_sim_type == sim_gemc .OR. int_sim_type == sim_gemc_ig .OR. &
       int_sim_type == sim_gemc_npt) THEN
     ! we will determine whether the individual probabilties are defined
     ! for transfer between the boxes
     l_prob_box_swap = .false.
     ALLOCATE(prob_swap_boxes(nbr_boxes,nbr_boxes))

     REWIND(inputunit)
     line_nbr = 0

     mainLoop:DO
        line_nbr = line_nbr + 1
        CALL Read_String(inputunit,line_string,ierr)

        IF (line_string(1:15) == '# Prob_Box_Swap') THEN
           l_prob_box_swap = .TRUE.
           ! we found a section specifying the probabilities of transfer
           ! between boxes

           prob_swap_boxes(:,:) = 0.0_DP

           boxLoop:DO ibox = 1, nbr_boxes
              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,nbr_boxes + 1, nbr_entries,line_array,ierr)

              ! first check that the box number matches up

              this_box = String_To_Int(line_array(1))

              IF (this_box /= ibox ) THEN
                 err_msg = ''
                 err_msg(1) = 'Box ids do not match up in # Prob_Box_Swap section'
                 CALL Clean_Abort(err_msg,'Get_Move_Probabilities')
              END IF

              DO kbox = 1, nbr_boxes
                 prob_swap_boxes(ibox,kbox) = String_To_Double(line_array(kbox+1))
                 IF (ibox == kbox ) THEN

                    IF (prob_swap_boxes(ibox,kbox) /= 0.0_DP) THEN
                       err_msg = ''
                       err_msg(1) = 'Nonzero probability specified for transfer in the same box for box'                               
                       err_msg(2) = 'Intra box swap not yet implemented'
                       err_msg(3) = 'Box id '//Int_To_String(kbox)
                       CALL Clean_Abort(err_msg,'Get_Move_Probabilities')
                    END IF
                 END IF

              END DO

              ! Now from these relative probabilities, obtain cumulative
              ! probabilities

              DO kbox = 1, nbr_boxes

                 IF (kbox > 1 ) THEN
                    prob_swap_boxes(ibox,kbox) = prob_swap_boxes(ibox,kbox) + &
                         prob_swap_boxes(ibox,kbox-1)
                 END IF

              END DO

              IF (ABS(prob_swap_boxes(ibox,nbr_boxes)-1.0_DP) > 0.000001_DP) THEN
                 err_msg = ''
                 err_msg(1) = 'Probabilities of box transfer do not add up to 1'
                 err_msg(2) = 'Box id'//Int_To_String(ibox)
                 CALL Clean_Abort(err_msg,'Get_Move_Probabilities')
              END IF


           END DO boxLoop

           EXIT

        ELSE IF ( line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN

           WRITE(logunit,*)
           WRITE(logunit,*) 'Probabilties of swap for individual boxes not specified'
           WRITE(logunit,*) 'Defaulting to uniform probabilities'

           EXIT

        END IF

     END DO mainLoop

     IF (.NOT. l_prob_box_swap) THEN

        prob_box_swap = 1.0_DP/REAL(nbr_boxes-1,DP)
        prob_swap_boxes(:,:) = prob_box_swap

        DO ibox = 1, nbr_boxes

           DO kbox = 2, nbr_boxes - 1
              prob_swap_boxes(ibox,kbox) = prob_swap_boxes(ibox,kbox-1) + &
                   prob_swap_boxes(ibox,kbox)
           END DO

           prob_swap_boxes(ibox,ibox) = 0.0_DP           
           prob_swap_boxes(ibox,nbr_boxes) = 1.0_DP

        END DO

     END IF

     ! log the probabilties of choosing a pair of box
     WRITE(logunit,*)
     WRITE(logunit,*) '******** Box pair selection probabilities *********'
     WRITE(logunit,*)

     DO ibox = 1, nbr_boxes

        DO kbox = 1, nbr_boxes

           WRITE(logunit,'(I3,2X,I3,2X,F10.6)') ibox, kbox, prob_swap_boxes(ibox,kbox)

        END DO

     END DO

  END IF

  IF( int_sim_type == sim_npt .OR. int_sim_type == sim_gemc .OR. &
       int_sim_type == sim_gemc_npt .OR. & 
       int_sim_type == sim_gemc_ig) THEN

     IF (prob_volume == 0.0_DP) THEN
        err_msg = ''
        err_msg(1) = 'Volume probability is not specified for'
        err_msg(2) = TRIM(sim_type)//' ensemble'
        CALL Clean_Abort(err_msg,'Get_Pressure_Info')
     END IF
  END IF

  IF ( int_sim_type == sim_gcmc) THEN
     IF ( (prob_insertion == 0.0_DP) .OR. (prob_deletion == 0.0_DP)) THEN
        err_msg = ''
        err_msg(1) = 'Either probability of insertion or deletion is not specified'
        CALL Clean_Abort(err_msg,'Get_Move_Probabilities')
     END IF

  ELSE IF (int_sim_type == sim_gemc .OR. int_sim_type == sim_gemc_ig .OR. &
       int_sim_type == sim_gemc_npt) THEN


     IF ( (prob_swap == 0.0_DP) .OR. (prob_volume == 0.0_DP)) THEN
        WRITE(logunit,*)
        WRITE(logunit,*) '**** WARNING *****      '
        WRITE(logunit,*)
        WRITE(logunit,*) 'prob_swap and prob_vol are zero for GEMC simulation'
        WRITE(logunit,*) 'If it is production run make sure they are non-zero' 
        WRITE(logunit,*)  
        WRITE(logunit,*) '**** WARNING *****      ' 
        WRITE(logunit,*)
        !        err_msg = ''
        !        err_msg(1) = 'GEMC ensemble specified'
        !        err_msg(2) = 'Either probability for volume move or particle swap is not specified'
        !        CALL Clean_Abort(err_msg,'Get_Move_Probabilities')
     END IF


  END IF

  WRITE(logunit,*) '***** Finished reading the move probability info *****'
  WRITE(logunit,'(A20,I4)') 'Number of moves is :', num_moves
  WRITE(logunit,*)

  movetime(:) = 0.0_DP
  ! For the cumulative probability

  cut_trans = prob_trans
  cut_rot = cut_trans + prob_rot
  cut_torsion = cut_rot + prob_torsion
  cut_volume = cut_torsion + prob_volume
  cut_angle = cut_volume + prob_angle
  cut_insertion = cut_angle + prob_insertion
  cut_deletion = cut_insertion + prob_deletion
  cut_swap = cut_deletion + prob_swap
  cut_regrowth = cut_swap + prob_regrowth
  cut_ring = cut_regrowth + prob_ring
  cut_atom_displacement = cut_ring + prob_atom_displacement


  IF (ABS(cut_atom_displacement-1.0_DP) > tiny_number ) THEN

     err_msg = ""
     err_msg(1) = 'Probabilities do not add upto 1.0'
     CALL Clean_Abort(err_msg,'Get_Move_Probabilities')

  END IF

END SUBROUTINE Get_Move_Probabilities


 
!*****************************************************************************************

SUBROUTINE Get_Start_Type
  ! This subroutine will read in the initial coordinates from an input file.
  ! There are three options to start a run
  ! 'make_config' --- will attempt to generate an initial configuration by random insertion
  !               --- molecules
  ! 'read_old' --- read from an exisiting file
  ! 'checkpoint'  --- read from a crash file

  INTEGER :: ierr, line_nbr, nbr_entries, i,j, ibox
  CHARACTER(120) :: line_string, line_array(20)
  CHARACTER(1) :: first_character
  CHARACTER(4) :: symbol

  REAL(DP) :: total_mass, this_mass
  
  ierr = 0
  line_nbr = 0

  REWIND(inputunit)

  ! Start reading input file to read the type of initial coordinates
  
  inputLOOP:DO
     line_nbr = line_nbr + 1

     CALL Read_String(inputunit,line_string,ierr)
     IF ( ierr /= 0 ) THEN
        err_msg = ""
        err_msg(1) = "Error encoutered while reading initial coordinate information"
        CALL Clean_Abort(err_msg,'Get_Initial_coordinates_Info')
     END IF

     IF(line_string(1:12) == '# Start_Type') THEN
        ! we entered the section of input file that contains information on
        ! initial coordinates
        WRITE(logunit,*) 
        WRITE(logunit,*) '******Initial Coordinate Info **********'
        Start_Type_LOOP: DO 
           line_nbr = line_nbr + 1
           CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
           IF (line_array(1) == 'make_config') THEN

              start_type = 'make_config'
              ALLOCATE(nmol_actual(nspecies,nbr_boxes),Stat = AllocateStatus)
              IF (AllocateStatus /= 0 ) THEN
                 write(*,*)'memory could no tbe allocated for nmol_actual array'
                 write(*,*)'stopping'
                 STOP
              END IF
              WRITE(logunit,*)
              WRITE(logunit,*) 'Initial configuration will be generated'
              
              ! read in the initial geometry to generate initial configuration
              ! for each of the species
              
              ! Allocate memory for each of the file names
              ALLOCATE(init_geomfile(nspecies))
            
              DO i = 1, nspecies
                 line_nbr = line_nbr + 1
                 CALL Parse_String(inputunit,line_nbr,nbr_boxes,nbr_entries,line_array,ierr)
                 ! Check the string for alphanumeric characters
!!$                 CALL Check_String(line_array(1),ierr)
!!$                 
!!$                 IF (ierr /= 0) THEN
!!$                    err_msg = ""
!!$                    err_msg(1) = 'An error in the input line ' // TRIM(Int_to_String(line_nbr)) &
!!$                         // ' of input file.'
!!$                    CALL Clean_Abort(err_msg,'Get_Initial_Coordinates_Info')
!!$                 END IF

!!$                 init_geomfile(i) = line_array(1)
!!$
!!$                 WRITE(logunit,*) 
!!$                 WRITE(logunit,'(A,T50,I3,A)') 'The initial geometry file for species ', i , ' is'
!!$                 WRITE(logunit,*) ADJUSTL(init_geomfile(i))

                 ! assign actual number of molecules of this species in each box
                 
                 DO ibox = 1, nbr_boxes
                    nmol_actual(i,ibox) = String_To_Int(line_array(ibox))
                    WRITE(logunit,'(A41,2x,I2,2X,A7,2X,I2,2X,A2,2X,I6)') 'Starting number of molecules of species', i, ' in box ', ibox, 'is',  nmol_actual(i,ibox)
                 END DO
 
                 species_list(i)%nmoltotal = SUM(nmol_actual(i,:))

                 IF (species_list(i)%nmoltotal > nmolecules(i)) THEN
                    err_msg = ''
                    err_msg(1) = 'Actual number of molecules of species' // INT_To_String(i) // 'is'
                    err_msg(2) = 'greater than maximum allowed'
                    CALL Clean_Abort(err_msg,'Get_Initial_Coordinates_Info')
                 END IF

                 WRITE(logunit,*) 
                 WRITE(logunit,'(A36,2X,I2,2X,A21,2X,I6)') 'Total number of molecules of species ', i , ' present initially is', species_list(i)%nmoltotal

!!$                 ! --- open the geometry file and read in the coordinates 
!!$                 OPEN(unit=init_geomunit,file=init_geomfile(i),status='old')
!!$                 ! the file is in xyz format
!!$
!!$                 ! obtain COM of the input geometry as well
!!$                 total_mass = 0.0_DP
!!$                 
!!$                 species_list(i)%xcom = 0.0_DP
!!$                 species_list(i)%ycom = 0.0_DP
!!$                 species_list(i)%zcom = 0.0_DP
!!$                 
!!$
!!$                 DO j = 1, natoms(i)
!!$                    READ(init_geomunit,*) symbol, init_list(j,1,i)%rxp, init_list(j,1,i)%ryp, init_list(j,1,i)%rzp
!!$                    
!!$                    this_mass = nonbond_list(j,i)%mass
!!$                    total_mass = total_mass + this_mass
!!$                    
!!$                    species_list(i)%xcom = species_list(i)%xcom + this_mass * &
!!$                         init_list(j,1,i)%rxp
!!$                    species_list(i)%ycom = species_list(i)%ycom + this_mass * &
!!$                         init_list(j,1,i)%ryp
!!$                    species_list(i)%zcom = species_list(i)%zcom + this_mass * &
!!$                         init_list(j,1,i)%rzp
!!$                    
!!$                 END DO
!!$                 
!!$                 CLOSE(UNIT=init_geomunit)
!!$
!!$                 species_list(i)%xcom = species_list(i)%xcom / total_mass
!!$                 species_list(i)%ycom = species_list(i)%ycom / total_mass
!!$                 species_list(i)%zcom = species_list(i)%zcom / total_mass
                 
              END DO
              
              EXIT inputLOOP

           ELSE IF (line_array(1) == 'read_old') THEN
              ! in this case we will read in the information of coordinates for
              ! all the boxes so that there must be nbr_boxes lines following the
              ! keyword.
              ALLOCATE(old_config_file(nbr_boxes))
              
              start_type = 'read_old'

              WRITE(logunit,*)
              WRITE(logunit,*) 'Configurations will be read from old files'
              DO i = 1,nbr_boxes
                 line_nbr = line_nbr + 1
                 CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
                 ! Make sure that the characters of the string are alphanumeric with
                 ! a possibility of a . (dot). The first character must be an alphabet
                 CALL Check_String(line_array(1),ierr)
                 IF (ierr /= 0 ) THEN
                    err_msg = ""
                    err_msg(1) = 'An error in the input line ' // TRIM(Int_to_String(line_nbr)) &
                         // ' of input file.'
                    CALL Clean_Abort(err_msg,'Get_Initial_Coordinates_Info')
                 END IF

                 WRITE(logunit,*)
                 WRITE(logunit,'(A,T40,I3,A,T50)')'Starting configuration for box ', i, ' is'
                 WRITE(logunit,*) ADJUSTL(line_array(1))
                 old_config_file(i) = TRIM(ADJUSTL(line_array(1)))
                 
              END DO
              EXIT inputLOOP
           ELSE IF (line_array(1) == 'checkpoint') THEN

              start_type = 'checkpoint'

              WRITE(logunit,*)
              WRITE(logunit,*) 'Starting configuration will be obtained from check point files'
              
              line_nbr = line_nbr + 1
              CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
              ! Make sure that the characters of the string are alphanumeric with
              ! a possibility of a . (dot). or _ (dash). The first character must be an alphabet
              CALL Check_String(line_array(1),ierr)
              write(*,*) line_array(1)
              IF (ierr /= 0 ) THEN
                 err_msg = ""
                 err_msg(1) = 'An error in the input line ' // TRIM(Int_to_String(line_nbr)) &
                      // ' of input file.'
                 CALL Clean_Abort(err_msg,'Get_Initial_Coordinates_Info')
              END IF
              
              restart_file = line_array(1)
              WRITE(logunit,*)
              WRITE(logunit,*)'Starting configuration is '
              WRITE(logunit,*) ADJUSTL(line_array(1))
              
              
              
           
              EXIT inputLOOP

           END IF
                       
        END DO Start_Type_LOOP
     END IF
     
  END DO InputLOOP
  WRITE(logunit,*) 
  WRITE(logunit,*) '*****Finished reading initial coordinate info *****'
END SUBROUTINE Get_Start_Type


!**************************************************************************************
SUBROUTINE Get_Run_Type
  ! The subroutine determines whether the run is equilibration or production
  ! During an equilibration run, widths of translation and rotational moves
  ! along with dihedral and angles will be changed to achieve a 50% acceptance.
  !*********************************************************************************

USE Energy_Routines
USE Rotation_Routines
USE Random_Generators

  IMPLICIT NONE

  INTEGER :: ierr, line_nbr, nbr_entries,i, ia 
  CHARACTER(120) :: line_string,line_array(20)
  LOGICAL :: overlap

  ierr = 0
  line_nbr = 0
  
  DO
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)
     IF (ierr /=0 ) THEN
        err_msg = ""
        err_msg(1) = 'Error while reading input file'
        CALL Clean_Abort(err_msg,'Get_Run_Type')
     END IF
     
     IF (line_string(1:10) == '# Run_Type') THEN
        
        WRITE(logunit,*)
        WRITE(logunit,*) '****** Reading Run_Type Information ******'

        
        line_nbr = line_nbr + 1
        CALL Parse_String(inputunit,line_nbr,2,nbr_entries,line_array,ierr)
        
        IF (line_array(1) == 'Equilibration' ) THEN
           run_style = line_array(1)
           int_run_style = run_equil

        ELSE IF (line_array(1) == 'Production') THEN
           run_style = line_array(1)
           int_run_style = run_prod
         
        ELSE IF (line_array(1) == 'Test') THEN
           run_style = line_array(1)
           int_run_style = run_test
           
        ELSE
  
           WRITE(logunit,*)
           err_msg = ""
           err_msg(1) = 'Run_Type not supported'
           CALL Clean_Abort(err_msg,'Get_Run_Type')

        END IF

        nupdate = String_To_Int(line_array(2))
        
        WRITE(logunit,*)
        WRITE(logunit,*) 'The input run type is ', TRIM(line_array(1))
        WRITE(logunit,*) 'Update frequency is ', nupdate

    IF (int_sim_type == sim_npt .OR. int_sim_type == sim_gemc .OR. &
         int_sim_type == sim_gemc_npt .OR. &
         int_sim_type == sim_gemc_ig) THEN
           
           IF (nbr_entries /= 3) THEN
              err_msg = ""
              err_msg(1) = 'Equilibration specified without the volume update'
              CALL Clean_Abort(err_msg,'Get_Run_Type')
           END IF
           
           nvol_update = String_To_Int(line_array(3))
           
           WRITE(logunit,*) 
           WRITE(logunit,*) 'Update frequency for adjusting maximum volume displacement is'
           WRITE(logunit,*) nvol_update
           
        END IF
        
        
        EXIT

     ELSE IF (line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN
        
        err_msg = ""
        err_msg(1) = 'Run_Type not specified in the input file.'
        CALL Clean_Abort(err_msg,'Get_Run_Type')
        
     END IF

  END DO
  
  WRITE(logunit,*)
  WRITE(logunit,*) '******* Finished reading run type information ******'
END SUBROUTINE Get_Run_Type

!*****************************************************************************************

SUBROUTINE Get_CBMC_Info
  ! The subroutine reads in the information on the starting seed for the simulation

  INTEGER :: ierr, line_nbr, nbr_entries,ibox
  CHARACTER(120) :: line_string,line_array(30)

  REWIND(inputunit)
  ierr = 0
  line_nbr = 0

  IF (MAXVAL(nfragments) == 0) RETURN

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)
     IF(ierr /= 0) THEN
        err_msg = ''
        err_msg(1) = 'Error while reading the Get_CBMC_Info'
        call clean_abort(err_msg,'Get_CBMC_Info')
     end if

     if(line_string(1:11) == '# CBMC_Info') then
        write(logunit,*)
        write(logunit,*) '***** reading CBMC info *********'
        line_nbr = line_nbr + 1
        call parse_string(inputunit,line_nbr,2,nbr_entries,line_array,ierr)
        kappa_ins = string_to_int(line_array(2))
        line_nbr = line_nbr + 1
        call parse_string(inputunit,line_nbr,2,nbr_entries,line_array,ierr)
        kappa_rot = string_to_int(line_array(2))
        line_nbr = line_nbr + 1
        call parse_string(inputunit,line_nbr,2,nbr_entries,line_array,ierr)
        kappa_dih = string_to_int(line_array(2))

        write(logunit,*)
        write(logunit,*) 'Writing out CBMC_Info'
        write(logunit,'(a,t35,2i12)') 'kappa for first bead insertion ', kappa_ins
        write(logunit,'(a,t35,2i12)') 'kappa for rotational bias', kappa_rot
        write(logunit,'(a,t35,2i12)') 'kappa for dihedral selection ', kappa_dih
        write(logunit,*)

        line_nbr = line_nbr + 1 

        call parse_string(inputunit,line_nbr,2,nbr_entries,line_array,ierr)

        DO ibox = 1, nbr_boxes
           rcut_CBMC(ibox) = String_To_Double(line_array(ibox+1))
           write(logunit,'(a,t35,i8,f12.2)') 'Smaller cutoff for CBMC for box ',ibox, rcut_CBMC(ibox)
        END DO
        write(logunit,*) '****** finished loading CBMC_Info ********'
        exit

     else if(line_string(1:3) == 'end' .or. line_nbr > 10000 ) then
        err_msg = ''
        err_msg(1) = 'CBMC_Info not specified'
        call clean_abort(err_msg,'Get_CBMC_Info')
     end if
  end do
end subroutine get_cbmc_info
!*****************************************************************************************
SUBROUTINE Get_Zig_By_Omega
!*************************************
! This subroutine is called when chemical potential is used in a GCMC simulation. It goes
! through the input file and looks for the keyword '# Zig_By_Omega_Info'. For all the
! species that exchange with the reserovir are included. Zig_By_Omega will be set to
! zero for rest of the species
!
! Written by Jindal Shah on 10/28/13
!*************************************

  INTEGER :: ierr, line_nbr, nbr_entries, i, spec_counter
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(inputunit)

  ierr = 0
  line_nbr = 0
  spec_counter = 0

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)
     
     IF (ierr /= 0 ) THEN
        err_msg = ''
        err_msg(1) = 'An error occurred while reading inputfile'
        err_msg(2) = inputfile
        CALL Clean_Abort(err_msg,'Get_Zig_By_Omega')

     END IF

     IF (line_string(1:19) == '# Zig_By_Omega_Info') THEN

        WRITE(logunit,*)
        WRITE(logunit,*) '***** Reading Zig_By_Omega Info ***'
        line_nbr = line_nbr + 1
        CALL Parse_String(inputunit,line_nbr,nspec_insert,nbr_entries,line_array,ierr)
        
        DO i = 1, nspecies
           IF (species_list(i)%int_species_type == int_sorbate) THEN
              spec_counter = spec_counter + 1
              species_list(i)%zig_by_omega = String_To_Double(line_array(spec_counter))
           ELSE
              species_list(i)%zig_by_omega = 0.0_DP
           END IF

           WRITE(logunit,*)
           WRITE(logunit,*) 'Zig / Omega for species ', i, ' is ', species_list(i)%zig_by_omega

        END DO

        EXIT

     ELSE IF (line_string(1:3) == 'END' .OR. line_nbr > 10000) THEN
        
        err_msg = ''
        err_msg(1) = 'Zig_By_Omega_Info section is missing from the input fiel'
        err_msg(2) = inputfile

        CALL Clean_Abort(err_msg,'Get_Zig_By_Omega')

     END IF

  END DO

END SUBROUTINE Get_Zig_By_Omega
!*****************************************************************************************

SUBROUTINE Get_Seed_Info
  ! The subroutine reads in the information on the starting seed for the simulation

  INTEGER :: ierr, line_nbr, nbr_entries
  CHARACTER(120) :: line_string,line_array(20)

  REWIND(inputunit)
  ierr = 0
  line_nbr = 0

  DO

     IF (start_type == 'checkpoint') THEN
        WRITE(logunit,*) 
        WRITE(logunit,*) 'Seed will be read from a checkpoint file'
        EXIT
     END IF
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)
     IF(ierr /= 0) THEN
        err_msg = ''
        err_msg(1) = 'Error while reading the seed info'
        CALL Clean_Abort(err_msg,'Get_Seed_Info')
     END IF

     IF(line_string(1:6) == '# Seed') THEN

        WRITE(logunit,*)
        WRITE(logunit,*) '***** Reading seed info *********'
        line_nbr = line_nbr + 1

        IF (int_sim_type == sim_mcf ) THEN 
     
           CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
           
           iseed = String_To_Int(line_array(1))

        ELSE
           
           CALL Parse_String(inputunit,line_nbr,2,nbr_entries,line_array,ierr)        
           iseed1 = String_To_Int(line_array(1))
           iseed3 = String_To_Int(line_array(2))

        END IF
        WRITE(logunit,*)
        WRITE(logunit,'(A25,1X,I12)') 'The starting seed s1 is:', s1
        WRITE(logunit,*)
        WRITE(logunit,*) '****** Finished loading the seed ********'

        EXIT
        
     ELSE IF(line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN

        err_msg = ''
        err_msg(1) = 'Seed not specified'
        CALL Clean_Abort(err_msg,'Get_Seed_Info')

     END IF        

  END DO
END SUBROUTINE Get_Seed_Info
!*********************************************************************************************

SUBROUTINE Get_Frequency_Info
  ! This subroutine obtains frequency for writing to various files

  INTEGER :: ierr, line_nbr, nbr_entries, ibox
  CHARACTER(120) :: line_string, line_array(20),movie_header_file, &
                     movie_xyz_file

  REWIND(inputunit)

  ierr = 0
  line_nbr = 0
  nthermo_freq = 0
  ncoord_freq = 0
  n_mcsteps = 0
  n_equilsteps = 0

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)
     
     IF ( ierr /= 0 ) THEN
        err_msg = ""
        err_msg(1) = "Error encoutered while reading initial coordinate information"
        CALL Clean_Abort(err_msg,'Get_Initial_coordinates_Info')
     END IF
     
     IF (line_string(1:16) == '# Frequency_Info') THEN
        ! We found a section that contains frequency info. We will read in the frequency information
        line_nbr = line_nbr + 1
        CALL Parse_String(inputunit,line_nbr,2,nbr_entries,line_array,ierr)
        IF (line_array(1) == 'freq_type') THEN

           IF (line_array(2) == 'Timed') THEN
              timed_run = .TRUE.
           ELSE
              timed_run = .FALSE.
           END IF

        ELSE  
           err_msg = ""
           err_msg(1) = 'A keyword is missing in the input file.'
           err_msg(2) = 'Check for freq_type.'
           CALL Clean_Abort(err_msg,'Get_Frequency_Info')
        END IF

        FreqLOOP: DO
           line_nbr = line_nbr + 1
           CALL Parse_String(inputunit,line_nbr,2,nbr_entries,line_array,ierr)
           IF(timed_run) THEN

              IF (line_array(1) == 'thermofreq') THEN
 
                 nthermo_freq = String_To_Int(line_array(2))

                 WRITE(logunit,*) 
                 WRITE(logunit,'(A,T50,I8,A)') 'Thermodynamic quantities will written at every', nthermo_freq, ' minutes.'

              ELSE IF (line_array(1) == 'coordfreq') THEN
              
                 ncoord_freq = String_To_Int(line_array(2))

                 WRITE(logunit,*)
                 WRITE(logunit,'(A,T50,I8,A)') 'Coordinates will be written at every', ncoord_freq, ' minutes.'

              ELSE IF (line_array(1) == 'Stop') THEN

                 n_mcsteps = String_To_Int(line_array(2))

                 WRITE(logunit,*) 
                 WRITE(logunit,'(A32,2X,I12,2X,A10)' ) 'The simulation will be run for ', n_mcsteps, ' minutes.'

              ELSE IF (line_array(2) == 'Done_Frequency_Info') THEN
              
                 WRITE(logunit,*)
                 WRITE(logunit,*)'*** Finished reading output info ******* '
              
                 EXIT FreqLOOP

              END IF

           ELSE

              IF (line_array(1) == 'Nthermofreq') THEN

                 nthermo_freq = String_To_Int(line_array(2))
              
                 WRITE(logunit,*) 
                 WRITE(logunit,'(A,T50,I8,A)') 'Thermodynamic quantities will written at every', nthermo_freq, ' MC steps.'

              ELSE IF (line_array(1) == 'Ncoordfreq') THEN
              
                 ncoord_freq = String_To_Int(line_array(2))

                 WRITE(logunit,*)
                 WRITE(logunit,'(A,T50,I8,A)') 'Coordinates will be written at every', ncoord_freq, ' MC steps.'
                 WRITE(logunit,*)

                 DO ibox = 1, nbr_boxes
                    movie_header_file = TRIM(run_name) // '.box' // TRIM(Int_To_String(ibox)) // '.H'
                    movie_xyz_file =    TRIM(run_name) // '.box' // TRIM(Int_To_String(ibox)) // '.xyz'
                    WRITE(logunit,'(A,T30,I1,A,T40,A)') 'movie header file for box ', ibox ,' is', TRIM(movie_header_file)
                    WRITE(logunit,'(A,T30,I1,A,T40,A)') 'movie_XYZ file for box ', ibox ,' is', TRIM(movie_xyz_file)
                    OPEN(unit=movie_header_unit+ibox,file=movie_header_file)
                    OPEN(unit=movie_xyz_unit+ibox,file=movie_xyz_file)
                 END DO

              ELSE IF (line_array(1) == 'MCsteps') THEN

                 n_mcsteps = String_To_Int(line_array(2))
                 WRITE(logunit,*) 
                 WRITE(logunit,'(A,T50,I10,A)' ) 'The simulation will be run for ', n_mcsteps, ' MC steps.'

                 ! # of equilibrium steps will be used only for the fragment generation
                 
              ELSE IF (line_array(1) == 'NequilSteps') THEN
                 
                 n_equilsteps = String_To_Int(line_array(2))
                 WRITE(logunit,*) 
                 WRITE(logunit, '(A,I10)') 'Number of equilibrium steps', n_equilsteps
                 

              ELSE IF (line_array(2) == 'Done_Frequency_Info') THEN
              
                 WRITE(logunit,*)
                 WRITE(logunit,*)'*** Finished reading output info ******* '
              
                 EXIT FreqLOOP

              END IF
  
           END IF

        END DO FreqLOOP
        
        EXIT
              
     ELSE IF (line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN
              
        err_msg = ""
        err_msg(1) = 'Frequency info for output not specified in the input file'
        CALL Clean_Abort(err_msg,'Get_Frequency_Info')

     END IF
           
  END DO

  ! Check to make sure that all the quantities are defined in the input file

  IF (n_mcsteps < 0 .OR. ncoord_freq <= 0 .OR. nthermo_freq <= 0) THEN
  
     err_msg = ""
     err_msg(1) = 'At least one of the keywords is missing in the input file.'

     IF(timed_run) THEN
        err_msg(2) = 'Check for coordfreq, thermofreq and Stop.'
     ELSE
        err_msg(2) = 'Check for Ncoordfreq, Nthermofreq, and MCsteps'
     END IF

     CALL Clean_Abort(err_msg,'Get_Frequency_Info')

  END IF
  
END SUBROUTINE Get_Frequency_Info
!*****************************************************************************************************

!*****************************************************************************************************
! The subroutine reads in information as to average properties or block properties will be computed
!*****************************************************************************************************
SUBROUTINE Average_Info

  CHARACTER(120) :: line_string, line_array(20)

  INTEGER :: ierr, line_nbr, this_average, nbr_entries

  REWIND(inputunit)
  line_nbr = 0
  ierr = 0

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)
     IF (ierr /= 0 ) THEN
        err_msg = ''
        err_msg(1) = 'Error encountered while reading input on average'
        CALL Clean_Abort(err_msg,'Average_Info')
     END IF

     IF (line_string(1:14) == '# Average_Info' ) THEN
        ! section on averages found
        line_nbr = line_nbr + 1
        CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
        this_average = String_To_Int(line_array(1))
        
        WRITE(logunit,*) 
        WRITE(logunit,*) '*********** Average Information ***************'
        WRITE(logunit,*)

        IF (this_average == 0 ) THEN
           block_average = .TRUE.
           WRITE(logunit,*) 'Block averages will be output'
          
        ELSE
           block_average = .FALSE.
           WRITE(logunit,*) 'Instantaneous values will be output'
        END IF
        
        WRITE(logunit,*) 
        WRITE(logunit,*) '*********** Ending Average Section ************'

        EXIT
     
     ELSE IF ( line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN

        err_msg = ''
        err_msg(1) = 'Average information is missing'
        CALL Clean_Abort(err_msg,'Average_Info')
        
     END IF
       

  END DO

END SUBROUTINE Average_Info

SUBROUTINE Get_Property_Info
  !***************************************************************************************************
  ! The subroutine obtains the information on what properties need to be written for output and
  ! how many files will be used for the output. The property section in the input file is identified
  ! with '# Property_Info'. The number of times this keyword is found indicates total number of
  ! of files that will be written. Following the keyword are the keywords for property output that
  ! will be written in respective files.
  !***************************************************************************************************

USE Run_Variables, ONLY: cpcollect

  INTEGER :: ierr, line_nbr, nbr_properties, max_properties, nbr_entries
  INTEGER :: i, j, this_box, ibox, is, average_id, ifrac
  CHARACTER(120) :: line_string, line_array(20)
  CHARACTER(12) :: extension
  CHARACTER(9) :: extension1
  CHARACTER(17) :: extension2

  REWIND(inputunit)

  ALLOCATE(nbr_prop_files(nbr_boxes))
  
  ierr = 0
  line_nbr = 0
  nbr_prop_files(:) = 0
  max_properties = 0
  cpcollect = .FALSE.

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)
     IF ( ierr /= 0 ) THEN
        err_msg = ""
        err_msg(1) = "Error encoutered while reading property information"
        CALL Clean_Abort(err_msg,'Get_Property_Info')
     END IF

     IF(line_string(1:15) == '# Property_Info') THEN
        backspace(inputunit)
        CALL Parse_String(inputunit,line_nbr,3,nbr_entries,line_array,ierr)
        ! the third entry indicates the box for which the current property
        ! file is to be written
        this_box = String_To_Int(line_array(3))
        nbr_prop_files(this_box) = nbr_prop_files(this_box) + 1
        !--- Now go through the lines following this keyword upto a point where
        !-- a blank is encountered or a '#' or '# Property_Info'
        nbr_properties = 0
        innerLOOP: DO
           line_nbr = line_nbr + 1
           CALL Parse_String(inputunit,line_nbr,0,nbr_entries,line_array,ierr)

           IF (nbr_entries == 0) EXIT innerLOOP
           IF (line_array(1) == '#' .AND. line_array(2) == 'Property_Info') THEN
              ! we encountered a line with another property file
              backspace(inputunit)
              EXIT innerLOOP
           ELSE
              ! We encountered a property
              IF (line_array(1) == 'Nmols' .OR. line_array(1) == 'Density') THEN
                 ! there are as many properties to be written as there are species
                 nbr_properties = nbr_properties + nspecies
              ELSE IF (line_array(1) == 'Chemical_Potential') THEN
                 nbr_properties = nbr_properties + nspecies
                 cpcollect = .TRUE. 

              ELSE
                 ! this is a property for the system
                 nbr_properties = nbr_properties + 1
              END IF
              
              max_properties = MAX(nbr_properties,max_properties)
           END IF
           
        END DO innerLOOP

     ELSE IF( (line_string(1:3) == 'END' .OR. line_nbr > 10000)) THEN
        EXIT
     END IF
    
  END DO

   ! Now we will figure out which properties are to be output. Allocate
  ! property files related arrays.

  ! Name of the property files
  
  ALLOCATE(prop_files(MAXVAL(nbr_prop_files),nbr_boxes))
  ! initially assign the logical variable as true. it will be assigned to true
  ! or false based on the information given.
  ! Holds the name of the properties for a given file
  ALLOCATE(prop_output(max_properties,MAXVAL(nbr_prop_files),nbr_boxes))
  ! Store number of properties for each of these files
  ALLOCATE(prop_per_file(MAXVAL(nbr_prop_files),nbr_boxes))
  ALLOCATE(first_open(MAXVAL(nbr_prop_files),nbr_boxes))

  prop_files(:,:) = " "
  prop_output(:,:,:) = ""
  prop_per_file(:,:) = 0
  first_open(:,:) = .TRUE.

 
  
  REWIND(inputunit)

  ierr = 0
  line_nbr = 0
  nbr_prop_files = 0

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)
     
     
     IF ( ierr /= 0 ) THEN
        err_msg = ""
        err_msg(1) = "Error encoutered while reading property information"
        CALL Clean_Abort(err_msg,'Get_Property_Info')
     END IF
     
     IF(line_string(1:15) == '# Property_Info') THEN
        backspace(inputunit)
        CALL Parse_String(inputunit,line_nbr,3,nbr_entries,line_array,ierr)
        this_box = String_To_Int(line_array(3))
        nbr_prop_files(this_box) = nbr_prop_files(this_box) + 1
        ! Now read in each of the lines until a blank line or a line with 'END' statement
        ! is encountered. Store the name of the property file in the prop_output array
        nbr_properties = 0
        innerLOOP1 : DO
          
           line_nbr = line_nbr + 1
           CALL Parse_String(inputunit,line_nbr,0,nbr_entries,line_array,ierr)
           
           IF (nbr_entries == 0) THEN
              prop_per_file(nbr_prop_files(this_box),this_box) = nbr_properties
          
              EXIT innerLOOP1
           END IF
           IF (line_array(1) == "#" .AND. line_array(2) ==  "Property_Info") THEN
              ! Information on another file is given
              backspace(inputunit)
              ! Store the number of properties to be written for this file
              prop_per_file(nbr_prop_files(this_box),this_box) = nbr_properties
              EXIT innerLOOP1
           ELSE
              ! we are reading a property for this file
              
              IF ( line_array(1) == 'Nmols' .OR. line_array(1) == 'Density' .OR. &
                   line_array(1) == 'Chemical_Potential') THEN
                 DO is = 1, nspecies
                    nbr_properties = nbr_properties + 1
                    prop_output(nbr_properties,nbr_prop_files(this_box),this_box) = line_array(1) 
                 END DO

              ELSE
                 nbr_properties = nbr_properties + 1
                 prop_output(nbr_properties,nbr_prop_files(this_box),this_box) = line_array(1) 
              END IF
              
           END IF

        END DO innerLOOP1
           
     ELSE IF (line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN
        EXIT

     END IF

  END DO
 ! Name the files for output
  DO ibox = 1, nbr_boxes
     DO i = 1, nbr_prop_files(ibox)
        extension = '.box' // TRIM(Int_To_String(ibox)) // '.prp' // TRIM(Int_To_String(i))
        CALL Name_Files(run_name,extension,prop_files(i,ibox))
     END DO
  END DO           

  DO ibox = 1, nbr_boxes
     IF ( nbr_prop_files(ibox) /= 0) THEN
        WRITE(logunit,*)
        WRITE(logunit,*) '**** Writing property output information *****'
        WRITE(logunit,*)
        WRITE(logunit,'(A48,2X,I3,2X,A8,I2)') 'Total number of property files to be written is ',  &
             nbr_prop_files(ibox), ' for box ', ibox
        WRITE(logunit,'(A42,2X,I2)') 'Maximum number of properties per file is ', max_properties

        WRITE(logunit,*) 'Writing the name of the property files and the corresponding property output'
        DO i = 1, nbr_prop_files(ibox)
           WRITE(logunit,*)
           WRITE(logunit,'(A15,2x,I2,2X,A3,2X,A)') 'Property file ', i, ' is ', TRIM(prop_files(i,ibox))
           WRITE(logunit,*) 'Properties output in these files are'
           WRITE(logunit,*)
           DO j = 1, prop_per_file(i,ibox)
              WRITE(logunit,*) TRIM(prop_output(j,i,ibox))
              
           END DO
        END DO
        
        WRITE(logunit,*) 
        WRITE(logunit,*) '***** Finished writing property information ******'
        
     ELSE
        WRITE(logunit,*)
        WRITE(logunit,*) 'No property output files will be generated for box ', ibox
     END IF
  END DO

  IF(cpcollect) THEN

     DEALLOCATE(locate)
     DEALLOCATE(molecule_list)
     DEALLOCATE(atom_list)

     ALLOCATE(locate(MAXVAL(nmolecules)+1,nspecies))
     ALLOCATE(molecule_list(MAXVAL(nmolecules)+1,nspecies))
     ALLOCATE(atom_list(MAXVAL(natoms),MAXVAL(nmolecules)+1,nspecies))

  END IF

END SUBROUTINE Get_Property_Info
!**********************************************************************************************

SUBROUTINE Copy_Inputfile

!**********************************************************************************************
!
! The subroutine copies the inputfile to the logfile. Thus one can easily rerun the simulation
! Only the information that is needed to repeat the simulations will be written. The comments
! will be stripped off.
!
! Written by Jindal Shah on 02/22/08
!
!***********************************************************************************************

  INTEGER :: ierr, line_nbr
  CHARACTER(120) :: line_string
  

  REWIND(inputunit)
  
  ierr = 0
  line_nbr = 0

  WRITE(logunit,*) 
  WRITE(logunit,*) '**** Copying Inputfile ****** '
  
  DO
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)

     IF ( ierr /= 0 ) THEN
        err_msg = ''
        err_msg(1) = 'Error while reading inputfile'
        CALL Clean_Abort(err_msg,'Copy_Inputfile')
     END IF

     ! Read the first character of the line_string, if it is a comment then
     ! skip the output

     IF (line_string(1:1) /= '!') THEN
        WRITE(logunit,*) TRIM(line_string)
     END IF
     IF (line_string(1:3) == 'END' .OR. line_nbr > 10000 ) EXIT
     

  END DO

  WRITE(logunit,*) 
  WRITE(logunit,*) '**** Finished writing input file *******'

END SUBROUTINE Copy_Inputfile
  

SUBROUTINE Get_Rcutoff_Low
  !**********************************************************************************************
  !
  ! The subroutine reads in the cutoff distance below which the probability of finding two atoms
  ! in a simulation is vanishingly low. It is used in the energy routines to quickly reject moves
  ! that bring two atoms closer than this distance.
  !
  ! Written by Jindal Shah on 02/25/08
  !
  !**********************************************************************************************

  INTEGER :: ierr, line_nbr, nbr_entries
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(inputunit)
  ierr = 0
  line_nbr = 0

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)
     
     IF (ierr /= 0) THEN
        err_msg = ''
        err_msg(1) = 'Error occurred while reading input file.'
        CALL Clean_Abort(err_msg,'Get_Rcutoff_Low')
     END IF
     
     ! Read each line and match it with the keyword '# Rcutoff_Low'
     
     IF (line_string(1:13) == '# Rcutoff_Low') THEN
        ! On the next line read the information on the cutoff distance
        line_nbr = line_nbr + 1
        CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
        rcut_low = String_To_Double(line_array(1))
        
        WRITE(logunit,*) 
        WRITE(logunit,'(A25,2X,F6.3,2X,A10)') 'MC low cutoff distance is ', rcut_low, ' Angstrom'
        WRITE(logunit,*)

        EXIT

     ELSE IF (line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN

        err_msg = ''
        err_msg(1) = 'MC low cutoff is not specified.'
        CALL Clean_Abort(err_msg,'Get_Rcutoff_Low')

     END IF
     

  END DO
END SUBROUTINE Get_Rcutoff_Low
!----------------------------------------------------------------------------------------

SUBROUTINE Get_File_Info
  !---------------------------------------------------------------------------------------
  ! This subroutine reads in the file information for generation of fragment library
  ! It will then be used to store information on the fragment configurations
  !
  ! Written by Jindal Shah on 09/10/08
  !
  !---------------------------------------------------------------------------------------

  USE File_Names

  INTEGER :: ierr, nbr_entries, line_nbr, is

  CHARACTER(120) :: line_array(20), line_string
  
  ierr = 0
  REWIND(inputunit)
  line_nbr = 0

  ALLOCATE(frag_file(nspecies))

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)

     IF ( ierr /= 0 ) THEN
        err_msg = ''
        err_msg(1) = 'Error reading input file'
        CALL Clean_Abort(err_msg,'Get_File_Info')
     END IF

     ! Read the input file upto # File_Info

     IF (line_string(1:11) == '# File_Info') THEN
        ! parse the string to read in the files for each species
        DO is = 1, nspecies
           line_nbr = line_nbr + 1
           CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
           frag_file(is) = line_array(1)

           WRITE(logunit,'(A35,2X,i3,A5)') 'Configuration file for species', is, ' is', TRIM(frag_file(is))
           
        END DO

        EXIT

     ELSE IF (line_nbr > 10000 .OR. line_string(1:3) == 'END') THEN
        err_msg = ''
        err_msg(1) = 'File name is not specified'
        CALL Clean_Abort(err_msg,'Get_File_Info')

     END IF
     
  END DO

END SUBROUTINE Get_File_Info


SUBROUTINE Get_Energy_Check_Info

  IMPLICIT NONE

  INTEGER :: ierr, line_nbr, nbr_entries, is, ibox
  CHARACTER(120) :: line_string, line_array(20)

  REWIND(inputunit)

! determine the number of species to be simulated
  ierr = 0
  line_nbr = 0
  echeck_flag = .FALSE.

  DO
     line_nbr = line_nbr + 1

     CALL Read_String(inputunit,line_string,ierr)

     IF (ierr .NE. 0) THEN
        err_msg = ""
        err_msg(1) = "Error getting energy check information."
        CALL Clean_Abort(err_msg,'Get_Energy_Check_Info')
     END IF

     IF (line_string(1:13) == '# Echeck_Info') THEN

        echeck_flag = .TRUE.
        line_nbr = line_nbr + 1

        CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)

        iecheck = String_To_Int(line_array(1))

        EXIT

     ELSE IF (line_string(1:3) == 'END') THEN

        EXIT

     END IF

  END DO

  END SUBROUTINE Get_Energy_Check_Info 

SUBROUTINE Get_Mie_Nonbond
  !---------------------------------------------------------------------------------------
  ! This subroutine reads in the file information for nonbond Mie potential exponents
  ! for each species type.
  !
  ! Written by Brian Yoo and Eliseo Rimoldi on 02/28/15
  !
  !---------------------------------------------------------------------------------------

  USE File_Names

  INTEGER :: ierr, nbr_entries, line_nbr, is, Mk, Mi, Mj
  CHARACTER(120) :: line_array(20), line_string
  ierr = 0
  REWIND(inputunit)
  line_nbr = 0
  Mk = 1

  ALLOCATE(mie_nlist(nspecies*(nspecies+1)/2))
  ALLOCATE(mie_mlist(nspecies*(nspecies+1)/2))
  ALLOCATE(mie_Matrix(nspecies, nspecies))

  DO
     line_nbr = line_nbr + 1
     CALL Read_String(inputunit,line_string,ierr)

     IF ( ierr /= 0 ) THEN
        err_msg = ''
        err_msg(1) = 'Error reading input file'
        CALL Clean_Abort(err_msg,'Get_Mie_Nonbond')
     END IF

     ! Read the input file up to # Mie_Nonbond

     IF (line_string(1:13) == '# Mie_Nonbond') THEN
        ! create symmetric matrix for index of species (e.g. for 3 species it will create
	! the following matrix [1,2,3;2,4,5;3,5,6]; This matrix is used to identify
	! the specified mie_n and mie_m exponents for a given species type.
        DO Mi = 1, nspecies
           DO Mj = Mi, nspecies
              mie_Matrix(Mi,Mj) = Mk
              mie_Matrix(Mj,Mi) = Mk
              Mk = Mk + 1
           END DO
        END DO

        ! parse the string to read in the files for each species
        DO is = 1, nspecies*nspecies
           line_nbr = line_nbr + 1
           CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
           mie_nlist(mie_Matrix(String_To_Int(line_array(1)),String_To_Int(line_array(2)))) = String_To_Double(line_array(3))
           mie_mlist(mie_Matrix(String_To_Int(line_array(1)),String_To_Int(line_array(2)))) = String_To_Double(line_array(4))

           WRITE(logunit,*) 'Mie exponent for ', line_array(1), 'and', line_array(2),  ' is', line_array(3), 'and', line_array(4)

        END DO

        EXIT

     ELSE IF (line_nbr > 10000 .OR. line_string(1:3) == 'END') THEN
        err_msg = ''
        err_msg(1) = 'Mie potentials not specified'
        CALL Clean_Abort(err_msg,'Get_Mie_Nonbond')

     END IF

  END DO

END SUBROUTINE Get_Mie_Nonbond


SUBROUTINE Get_Lattice_File_Info
    !-----------------------------------------------------------------
    !
    ! This soubroutine gets the name of the xyz file containing lattice
    ! coordinates
    !
    ! First written by Jindal Shah on 10/10/12
    !
    !-------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: line_nbr, ierr, nbr_entries
    CHARACTER*120 :: line_string, line_array(20)

    REWIND(inputunit)

    line_nbr = 0

    DO 
       line_nbr = line_nbr + 1
       CALL Read_String(inputunit,line_string,ierr)
       
       IF (ierr /= 0 ) THEN
          err_msg = ''
          err_msg(1) = 'Error reading lattice info in the inputfile'
          err_msg(2) = 'aborting'
          CALL Clean_Abort(err_msg,'Get_Lattice_File_Info')
       END IF

       IF (line_string(1:19) == '# Lattice_File_Info') THEN
          ! we found the section on Lattice File Info
          WRITE(logunit,*)
          WRITE(logunit,*) '***** Found section on lattice file info ****'
          WRITE(logunit,*) 

          line_nbr = line_nbr + 1
          CALL Parse_String(inputunit,line_nbr,1,nbr_entries,line_array,ierr)
          
          lattice_file = line_array(1)

          WRITE(logunit,*) 'Lattice coordinates are saved in the file'
          WRITE(logunit,*) lattice_file

          EXIT

       ELSE
          
          IF (line_nbr > 1000 .OR. line_string(1:3) == 'END') THEN
             
             err_msg = ''
             err_msg(1) = '# Lattice_File_Info section is missing in the input file'
             err_msg(2) = inputfile
             CALL Clean_Abort(err_msg,'Get_Lattice_File_Info')

          END IF

       END IF

    END DO
    

    
END SUBROUTINE Get_Lattice_File_Info


SUBROUTINE Get_Lattice_Coordinates

    !-- This subroutine obtains coordinates of the unit cell of the input
    !-- lattice.
    !
    !-- The coordinate file contains the following information
    !
    !--- first line reserved for comment
    !
    !-- second, third and fourth lines describe the cell matrix
    !
    !-- fifth line is for comment
    !-- sixth line onwards, the coordinates are specified in x,y and z. The first
    !-- column contains atomic id
    !
    !   First written by Jindal Shah on 10/10/12.
    !
    !-------------------------------------------------------------------

    IMPLICIT NONE

    INTEGER :: is, iatom
    CHARACTER*4 :: symbol
    ! Since this routine is called in the case of potential map generation
    ! and we do not call Box_Info during the map generation, set the 
    ! number of boxes to 1

    nbr_boxes = 1


    OPEN(UNIT=lattice_file_unit,FILE=lattice_file)

!    READ(lattice_file_unit,*)
!    READ(lattice_file_unit,*) box_list(1)%length(1,1), box_list(1)%length(1,2), &
!         box_list(1)%length(1,3)
!    READ(lattice_file_unit,*) box_list(1)%length(2,1), box_list(1)%length(2,2), &
!         box_list(1)%length(2,3)
!    READ(lattice_file_unit,*) box_list(1)%length(3,1), box_list(1)%length(3,2), &
!         box_list(1)%length(3,3)

    ! skip a line
    !
    ! At this point compute all the properties associated with the zeolite box
!    CALL Compute_Cell_Dimensions(1)

!    READ(lattice_file_unit,*)

    READ(lattice_file_unit,*)

    ! to read in the coordinates of the atoms, we need to first determine how
    ! many zeolite atoms there are
    
!    DO is = 1, nspecies
!       IF (species_list(is)%int_species_type == int_zeo) EXIT
!    END DO

!    IF ( is == (nspecies + 1)) THEN
       ! we could not locate a zeolite species

!       err_msg = ''
!       err_msg(1) = 'Simulation type is zeolite potential map generation'
!       err_msg(2) = 'However none of the species appears to be a zeolite species'
!       err_msg(3) = "Check the '# Species_Type' keyword in all the master mcf files"
!       CALL Clean_Abort(err_msg,'Get_Lattice_Coordinates')

!    END IF

    is = 1

    n_lat_atoms = natoms(is)

    ! loop over all the coordinates to obtain xyz coordinates
    ALLOCATE(x_lat(n_lat_atoms),y_lat(n_lat_atoms), z_lat(n_lat_atoms))

    DO iatom = 1, n_lat_atoms

       READ(lattice_file_unit,*) symbol, x_lat(iatom), y_lat(iatom), z_lat(iatom)

    END DO

    CLOSE(UNIT=lattice_file_unit)

    WRITE(logunit,*) 'Finished reading the lattice coordinates'
    
END SUBROUTINE Get_Lattice_Coordinates


END MODULE Input_Routines
