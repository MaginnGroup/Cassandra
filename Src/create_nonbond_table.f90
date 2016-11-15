!********************************************************************************
! CASSANDRA - Computational Atomistic Simulation Software at Notre Dame 
! for Research in Academia.
! http://molsim.wiki.zoho.com/
! Copyright (2007) University of Notre Dame.
! Authors: Ed Maginn (ed@nd.edu); Jindal Shah (jshah@nd.edu)
!********************************************************************************

!********************************************************************************
  SUBROUTINE Create_Nonbond_Table
!********************************************************************************
  ! This routine determines the number and identity of all unique atom types
  ! in the system. It then takes the associated vdw paramters of each atom
  ! and creates an interaction table using the specified mixing rule.
  ! Currently only supports LJ vdw type and LB or geometric mixing rules.

  ! The final product of this code are a bunch of matrices
  ! vdw_param1_table(itype,jtype), vdw_param2_table(itype,jtype), etc. 
  ! where for the LJ potential, vdw_param1_table(itype,jtype) is the epsilon
  ! to use between atoms of itype and jtype, and vdw_param2_table(itype,jtype)
  ! is the sigma to use between these two atom type. 

  ! Example: let's say you need to compute the LJ energy betweem atom 4 of species 2
  ! and atom 6 of species 1. 
  ! You first identify the unique atom type of these two atoms:
  ! itype = nonbond_list(4,2)%atom_type_number
  ! jtype = nonbond_list(6,1)%atom_type_number

  ! Then get the parameters:
  ! epsilon(itype, jtype) = vdw_param1_table(itype,jtype)
  ! sigma(itype,jtype) = vdw_param2_table(itype,jtype)
  ! Note that these matrices are symmetric, so 
  ! vdw_param1_table(itype,jtype) = vdw_param1_table(jtype,itype)

  ! Written: Sat Oct  6 08:58:06 MDT 2007
  ! Author: E. Maginn
  
  ! *** CALLS ***
  ! Clean_Abort
  ! 
  ! *** CALLED BY ***
  ! NVTMC_Control
  !
  ! Revision history:

  ! 04/08/09 (TR) : Removed any character evaluations.
  !
  ! 08/26/11 (JS) : Steele potential constants calculations
  !
  ! 03/27/14 Eliseo Rimoldi : custom (manually input) nonbond potentials added

!******************************************************************************
  USE Global_Variables
  USE Type_Definitions
  USE IO_Utilities
  USE File_Names
  
  IMPLICIT NONE

  LOGICAL :: repeat_type
  CHARACTER(6), DIMENSION(:), ALLOCATABLE :: temp_atomtypes ! same dimension as atom_name
  INTEGER :: ii, is, ia, itype, jtype, iset,jset,k,itype_1,itype_2
  REAL(DP) :: eps_i, eps_j, sig_i, sig_j
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: vdw_param_set

  ! Steele potential
  REAL(DP) :: sigma_ss, eps_ss, rho_s, delta_s, eps_sf
  !custom mixing rules
  INTEGER :: ierr,line_nbr,nbr_entries, is_1, is_2, ia_1, ia_2, itype_custom, jtype_custom
  INTEGER ::  i_type1, i_type2
  CHARACTER(120) :: line_string, line_array(20)


!******************************************************************************
  IF (verbose_log) THEN
     WRITE(logunit,*)
     WRITE(logunit,'(A)') 'Nonbond tables'
     WRITE(logunit,'(A80)') '********************************************************************************'
  END IF

  ALLOCATE(temp_atomtypes(1000), Stat=AllocateStatus)
  IF (AllocateStatus .NE. 0) THEN
     err_msg = ''
     err_msg(1) = ' ERROR: Not enough memory for temp_atomtypes '
     CALL Clean_Abort(err_msg,'create_nonbond_table')
  END IF

  ! Initialize the number of atom types and the temp_atomtypes
  nbr_atomtypes = 0
  temp_atomtypes = ''

  ! Compute the number of different atom types
  DO is = 1, nspecies
     
     DO ia = 1, natoms(is)
        
        repeat_type = .FALSE.

        IF (nonbond_list(ia,is)%vdw_type /= 'NONE') THEN

           !----------------------------------------------------------------
           ! Determine whether the atomtype has already been accounted for
           !----------------------------------------------------------------
           IF ((is .EQ. 1).AND.(ia .EQ. 1)) THEN
              ! If this is the first atomtype obtained, store it in temp_atomtypes
              nbr_atomtypes = nbr_atomtypes + 1
              temp_atomtypes(nbr_atomtypes) = nonbond_list(ia,is)%atom_name

              ! Store the unique number identifier
              nonbond_list(ia,is)%atom_type_number = nbr_atomtypes
             
           ELSE
              ! Loop over all current atomtypes to check if the atomtype has been accounted for
              !   If so, turn the fatomtype flag to true
              DO ii = 1, nbr_atomtypes
                 IF(nonbond_list(ia,is)%atom_name .EQ. temp_atomtypes(ii)) THEN

                    ! This atom name is already present. Do not advance counter
                    repeat_type = .TRUE.

                    ! Store the unique number identifier
                    nonbond_list(ia,is)%atom_type_number = ii

                 ENDIF
              ENDDO

              ! If the atomtype has not been accounted for, add it to the temp_atomtypes list
              IF(.NOT.repeat_type) THEN
                 nbr_atomtypes = nbr_atomtypes + 1
                 temp_atomtypes(nbr_atomtypes) = nonbond_list(ia,is)%atom_name

                 ! Store the unique number identifier
                 nonbond_list(ia,is)%atom_type_number = nbr_atomtypes
              ENDIF

           ENDIF
        ELSE
           ! atom has no atom_type
           nonbond_list(ia,is)%atom_type_number = 0
        ENDIF

     ENDDO

  ENDDO

  
! Write the number of different atom types to the logfile

  IF (verbose_log) THEN
     WRITE(logunit,'(A)') &
          ' There are '//TRIM(Int_To_String(nbr_atomtypes))//' different atom types in the system '
     DO ii = 1, nbr_atomtypes
        WRITE(logunit,'(3x,I3,2x,A6)') ii, temp_atomtypes(ii)
     ENDDO

     WRITE(logunit,*)

     DO is=1,nspecies
        WRITE(logunit,*)
        WRITE(logunit,'(A,T25,I3,3x,A)') 'species number and name:',is, molfile_name(is)
        WRITE(logunit,*) 'Name      number'
        WRITE(logunit,*) '------    ------'

        IF (natoms(is) < 100) THEN
          DO ia = 1, natoms(is)
             WRITE(logunit,'(A6,T10,I4)') nonbond_list(ia,is)%atom_name, &
                  nonbond_list(ia,is)%atom_type_number
          ENDDO
        ELSE
          WRITE(logunit,'(A)') 'Species has more than 100 atoms'
        END IF

     ENDDO

  END IF


  ! Create a character array containing the names of each unique atom type, with the index equal
  ! to the atom type number

  ALLOCATE(atom_type_list(nbr_atomtypes), Stat=AllocateStatus)
  IF (AllocateStatus .NE. 0) THEN
     err_msg = ''
     err_msg(1) = ' ERROR: Not enough memory for vdw atom_type_list'
     CALL Clean_Abort(err_msg,'ceate_nonbond_table')
  END IF
  atom_type_list = ""
  DO itype = 1, nbr_atomtypes
     atom_type_list(itype) = temp_atomtypes(itype)
  END DO
  IF (ALLOCATED(temp_atomtypes)) DEALLOCATE(temp_atomtypes)

  ! allocate arrays containing vdw parameters for all interaction pairs.
  ALLOCATE(vdw_param1_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
  ALLOCATE(vdw_param2_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
  ALLOCATE(vdw_param3_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
  ALLOCATE(vdw_param4_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
  ALLOCATE(vdw_param5_table(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)

  IF (AllocateStatus .NE. 0) THEN
     err_msg = ''
     err_msg(1) = ' ERROR: Not enough memory for vdw interaction table '
     CALL Clean_Abort(err_msg,'ceate_nonbond_table')
  END IF

  ! Allocate memory for total number bead types in each box
  ALLOCATE(nint_beads(nbr_atomtypes,nbr_boxes))
  
  ! allocated memory for vdw_param_set, to mark which interactions have parameters stored
  ALLOCATE(vdw_param_set(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
  IF (AllocateStatus .NE. 0) THEN
     err_msg = ''
     err_msg(1) = ' ERROR: Not enough memory for vdw_param_set'
     CALL Clean_Abort(err_msg,'create_nonbond_table')
  END IF
  vdw_param_set = 0

  ! Now determine the set of vdw parameters for each type of interaction and load them into vdw_param_table
  ! This is a brute force search - but it is fast.

  IF (verbose_log) THEN
     WRITE(logunit,*) 
     WRITE(logunit,*) 'Creating VDW interaction table'

  ! Write header for logfile output. Specific for the vdw style
     IF (int_vdw_style(1) == vdw_lj) THEN
        WRITE(logunit,'(X,A6,5X,A6,2X,A12,X,A12)') 'Atom 1','Atom 2', 'epsilon', 'sigma'
        WRITE(logunit,'(X,6X,5X,6X,2X,A12,X,A12)') 'amu A^2/ps^2', 'Ang'
     ELSEIF (int_vdw_style(1) == vdw_mie) THEN
        WRITE(logunit,'(X,A6,5X,A6,2X,A12,X,A12,X,A12,X,A12)') 'Atom 1','Atom 2', 'epsilon', 'sigma', &
              'rep-expt', 'disp-expt'
        WRITE(logunit,'(X,6X,5X,6X,2X,A12,X,A12)') 'amu A^2/ps^2', 'Ang'
     ENDIF
     WRITE(logunit,'(X,A79)') '-------------------------------------------------------------------------------'
  END IF

  ! Populate vdw_param?_table with like interactions
  DO is = 1, nspecies
     DO ia = 1, natoms(is)
        itype = nonbond_list(ia,is)%atom_type_number
        IF (itype > 0) THEN
           IF (vdw_param_set(itype,itype) == 0) THEN
              IF (nonbond_list(ia,is)%vdw_type == 'LJ') THEN
                 ! epsilon
                 IF (nonbond_list(ia,is)%vdw_param(1) <= tiny_number) THEN
                    vdw_param1_table(itype,itype) = 0.0_DP
                 ELSE
                    vdw_param1_table(itype,itype) = &
                                           nonbond_list(ia,is)%vdw_param(1)
                 END IF

                 ! sigma
                 IF (nonbond_list(ia,is)%vdw_param(2) <= tiny_number) THEN
                    vdw_param2_table(itype,itype) = 0.0_DP
                 ELSE
                    vdw_param2_table(itype,itype) = &
                                           nonbond_list(ia,is)%vdw_param(2)
                 END IF

                 ! Report parameters to logfile.
                 IF (verbose_log) THEN
                   WRITE(logunit,'(X,A6,5X,A6,2X,F12.4,X,F12.4)') &
                        atom_type_list(itype), atom_type_list(itype), &
                        vdw_param1_table(itype,itype), &
                        vdw_param2_table(itype,itype)
                 ENDIF

                 vdw_param_set(itype,itype) = 1
              ELSE IF (nonbond_list(ia,is)%vdw_type == 'Mie') THEN
                 ! epsilon
                 IF (nonbond_list(ia,is)%vdw_param(1) <= tiny_number) THEN
                    vdw_param1_table(itype,itype) = 0.0_DP
                 ELSE
                    vdw_param1_table(itype,itype) = &
                                     nonbond_list(ia,is)%vdw_param(1)
                 END IF

                 ! sigma
                 IF (nonbond_list(ia,is)%vdw_param(2) <= tiny_number) THEN
                    vdw_param2_table(itype,itype) = 0.0_DP
                 ELSE
                    vdw_param2_table(itype,itype) = &
                                     nonbond_list(ia,is)%vdw_param(2)
                 END IF

                 ! repulsive exponent
                 vdw_param3_table(itype,itype) = &
                                  nonbond_list(ia,is)%vdw_param(3)

                 ! dispersive exponent
                 vdw_param4_table(itype,itype) = &
                                  nonbond_list(ia,is)%vdw_param(4)

                 ! Report parameters to logfile.
                 IF (verbose_log) THEN
                   WRITE(logunit,'(X,A6,5X,A6,2X,F12.4,X,F12.4,X,F12.4,X,F12.4)') &
                        atom_type_list(itype), atom_type_list(itype), &
                        vdw_param1_table(itype,itype), vdw_param2_table(itype,itype), &
                        vdw_param3_table(itype,itype), vdw_param4_table(itype,itype)
                 ENDIF

                 vdw_param_set(itype,itype) = 1
              END IF
            END IF
        END IF
     END DO
  END DO

  ! Populate vdw_param?_table with mixed interactions
  IF (nbr_atomtypes > 1) THEN
     IF (verbose_log) THEN
        WRITE(logunit,'(X,A79)') '-------------------------------------------------------------------------------'
        WRITE(logunit,'(X,A,T25,A)') 'Mixing rule used is:', mix_rule
        WRITE(logunit,'(X,A79)') '-------------------------------------------------------------------------------'
     END IF

     IF (mix_rule == 'LB' .OR. mix_rule == 'geometric') THEN

        DO itype = 1, nbr_atomtypes
           DO jtype = itype+1, nbr_atomtypes

              IF (int_vdw_style(1) == vdw_lj) THEN
                 ! There are two vdw parameters that need mixing

                 ! epsilon
                 eps_i = vdw_param1_table(itype,itype)
                 eps_j = vdw_param1_table(jtype,jtype)
                 IF ( (eps_i <= tiny_number) .OR. (eps_j <= tiny_number) ) THEN
                    vdw_param1_table(itype,jtype) = 0.0_DP
                    ! for parameters with zero, avoid overflow and set to zero
                 ELSE
                    ! LB mixing rule: epsij = (epsi * epsj)^(1/2)
                    ! geometric mixing rule: epsij = (epsi * epsj)^(1/2)
                    vdw_param1_table(itype,jtype) = dsqrt(eps_i*eps_j)
                 ENDIF

                 ! sigma
                 sig_i = vdw_param2_table(itype,itype)
                 sig_j = vdw_param2_table(jtype,jtype)
                 IF ( (sig_i <= tiny_number) .OR. (sig_j <= tiny_number) ) THEN
                    vdw_param2_table(itype,jtype) = 0.0_DP
                    ! for parameters with zero, avoid overflow and set to zero
                 ELSE
                    ! LB mixing rule: sigmaij = 1/2 (sigmai + sigmaj)
                    IF (mix_rule == 'LB') &
                         vdw_param2_table(itype,jtype) = (sig_i + sig_j) * 0.5
                    ! geometric mixing rule: sigmaij = (sigmai * sigmaj)^(1/2)
                    IF (mix_rule == 'geometric') &
                         vdw_param2_table(itype,jtype) = dsqrt(sig_i*sig_j)
                 ENDIF

                 ! Report parameters to logfile.
                 IF (verbose_log) THEN
                   WRITE(logunit,'(X,A6,5X,A6,2X,F12.4,X,F12.4)') &
                        atom_type_list(itype), atom_type_list(jtype), &
                        vdw_param1_table(itype,jtype), vdw_param2_table(itype,jtype)
                 ENDIF

                 ! Mixed interactions are symmetric
                 vdw_param1_table(jtype,itype) = vdw_param1_table(itype,jtype)
                 vdw_param2_table(jtype,itype) = vdw_param2_table(itype,jtype)

                 ! Mark that the parameters have been set
                 vdw_param_set(itype,jtype) = 1
                 vdw_param_set(jtype,itype) = 1
              ELSEIF (int_vdw_style(1) == vdw_mie) THEN
                 ! Can only mix params 1 and 2 if params 3 and 4 are identical
                 IF (vdw_param3_table(itype,itype) == vdw_param3_table(jtype,jtype) .AND. &
                     vdw_param4_table(itype,itype) == vdw_param4_table(jtype,jtype)) THEN

                    ! epsilon
                    eps_i = vdw_param1_table(itype,itype)
                    eps_j = vdw_param1_table(jtype,jtype)
                    IF ( (eps_i <= tiny_number) .OR. (eps_j <= tiny_number) ) THEN
                       vdw_param1_table(itype,jtype) = 0.0_DP
                       ! for parameters with zero, avoid overflow and set to zero
                    ELSE
                       ! LB mixing rule: epsij = (epsi * epsj)^(1/2)
                       ! geometric mixing rule: epsij = (epsi * epsj)^(1/2)
                       vdw_param1_table(itype,jtype) = dsqrt(eps_i*eps_j)
                    ENDIF

                    ! sigma
                    sig_i = vdw_param2_table(itype,itype)
                    sig_j = vdw_param2_table(jtype,jtype)
                    IF ( (sig_i <= tiny_number) .OR. (sig_j <= tiny_number) ) THEN
                       vdw_param2_table(itype,jtype) = 0.0_DP
                       ! for parameters with zero, avoid overflow and set to zero
                    ELSE
                       ! LB mixing rule: sigmaij = 1/2 (sigmai + sigmaj)
                       IF (mix_rule == 'LB') &
                            vdw_param2_table(itype,jtype) = (sig_i + sig_j) * 0.5
                       ! geometric mixing rule: sigmaij = (sigmai * sigmaj)^(1/2)
                       IF (mix_rule == 'geometric') &
                            vdw_param2_table(itype,jtype) = dsqrt(sig_i*sig_j)
                    ENDIF

                    ! Exponents are the same
                    vdw_param3_table(itype,jtype) = vdw_param3_table(itype,itype)
                    vdw_param4_table(itype,jtype) = vdw_param4_table(itype,itype)

                    ! Report parameters to logfile.
                    IF (verbose_log) THEN
                      WRITE(logunit,'(X,A6,5X,A6,2X,F12.4,X,F12.4,X,F12.4,X,F12.4)') &
                           atom_type_list(itype), atom_type_list(jtype), &
                           vdw_param1_table(itype,jtype), vdw_param2_table(itype,jtype), &
                           vdw_param3_table(itype,jtype), vdw_param4_table(itype,jtype)
                    ENDIF

                    ! Mixed interactions are symmetric
                    vdw_param1_table(jtype,itype) = vdw_param1_table(itype,jtype)
                    vdw_param2_table(jtype,itype) = vdw_param2_table(itype,jtype)
                    vdw_param3_table(jtype,itype) = vdw_param3_table(itype,jtype)
                    vdw_param4_table(jtype,itype) = vdw_param4_table(itype,jtype)

                    ! Mark that the parameters have been set
                    vdw_param_set(itype,jtype) = 1
                    vdw_param_set(jtype,itype) = 1
                 ELSE
                    err_msg = ""
                    err_msg(1) = "Cross interactions for Mie potential must have same exponents"
                    CALL Clean_Abort(err_msg,'Create_Nonbond_Table')
                 END IF
              END IF

           ENDDO
        ENDDO

     ELSE IF (mix_rule == 'custom') THEN

       REWIND(inputunit)

       ierr = 0
       line_nbr = 0

       DO

         line_nbr = line_nbr + 1
         CALL Read_String(inputunit,line_string,ierr)

         IF (ierr .NE. 0) THEN
            err_msg = ""
            err_msg(1) = "Error reading mixing rules."
            CALL Clean_Abort(err_msg,'Create_Nonbond_Table')
         END IF

         IF (line_string(1:13) == '# Mixing_Rule') THEN
           ! 'custom' line
           line_nbr = line_nbr + 1
           READ(inputunit,*) 

           ! custom parameters
           CustomLOOP: DO
             line_nbr = line_nbr + 1
             CALL Parse_String(inputunit,line_nbr,0,nbr_entries,line_array,ierr)
             IF (nbr_entries < 2 .OR. line_array(1)(1:1) == '!') THEN
                EXIT CustomLOOP
             END IF

             ! find itype
             iset = 0
             DO itype = 1, nbr_atomtypes
                IF (atom_type_list(itype) == line_array(1)) THEN
                   iset = 1
                   EXIT
                END IF
             END DO
             IF (iset == 0) THEN
                WRITE(logunit,'(A)') 'Atom type ' // TRIM(line_array(1)) // ' on line ' // &
                   TRIM(Int_To_String(line_nbr)) // ' of input file not found. Skipping parameters.'
                CYCLE
             ELSE
                ! find jtype
                jset = 0
                DO jtype = 1, nbr_atomtypes
                   IF (atom_type_list(jtype) == line_array(2)) THEN
                      jset = 1
                      EXIT
                   END IF
                END DO
                IF (jset == 0) THEN
                   WRITE(logunit,'(A)') 'Atom type ' // TRIM(line_array(2)) // ' on line ' // &
                      TRIM(Int_To_String(line_nbr)) // ' of input file not found. Skipping parameters.'
                   CYCLE
                END IF
             END IF

             IF (vdw_param_set(itype,jtype) == 0) THEN
                ! Load custom parms
                IF (int_vdw_style(1) == vdw_lj) THEN
                  !Convert epsilon to atomic units amu A^2/ps^2
                  vdw_param1_table(itype,jtype) = kboltz * String_To_Double(line_array(3))
                  !Sigma
                  vdw_param2_table(itype,jtype) = String_To_Double(line_array(4))
                  IF (verbose_log) THEN
                    WRITE(logunit,'(X,A6,5X,A6,2(X,F12.4))') &
                         atom_type_list(itype), atom_type_list(jtype), &
                         vdw_param1_table(itype,jtype), &
                         vdw_param2_table(itype,jtype)
                  END IF
                  !Also define parms for reversed indices
                  vdw_param1_table(jtype,itype) = vdw_param1_table(itype,jtype)
                  vdw_param2_table(jtype,itype) = vdw_param2_table(itype,jtype)
 
                  ! Mark that the parameters have been set
                  vdw_param_set(itype,jtype) = 1
                  vdw_param_set(jtype,itype) = 1

                ELSEIF (int_vdw_style(1) == vdw_mie) THEN
                  !Convert epsilon to atomic units amu A^2/ps^2
                  vdw_param1_table(itype,jtype) = kboltz * String_To_Double(line_array(3))
                  !Sigma
                  vdw_param2_table(itype,jtype) = String_To_Double(line_array(4))
                  !Repulsive exponent
                  vdw_param3_table(itype,jtype) = String_To_Double(line_array(5))
                  !Dispersive exponent
                  vdw_param4_table(itype,jtype) = String_To_Double(line_array(6))
                  IF (verbose_log) THEN
                    WRITE(logunit,'(X,A6,5X,A6,4(X,F12.4))') &
                         atom_type_list(itype), atom_type_list(jtype), &
                         vdw_param1_table(itype,jtype), &
                         vdw_param2_table(itype,jtype), &
                         vdw_param3_table(itype,jtype), &
                         vdw_param4_table(itype,jtype)
                  END IF
                  !Also define parms for reversed indices
                  vdw_param1_table(jtype,itype) = vdw_param1_table(itype,jtype)
                  vdw_param2_table(jtype,itype) = vdw_param2_table(itype,jtype)
                  vdw_param3_table(jtype,itype) = vdw_param3_table(itype,jtype)
                  vdw_param4_table(jtype,itype) = vdw_param4_table(itype,jtype)

                  ! Mark that the parameters have been set
                  vdw_param_set(itype,jtype) = 1
                  vdw_param_set(jtype,itype) = 1
                END IF
             ELSE
                err_msg = ''
                err_msg(1) = 'Parameters for ' // TRIM(atom_type_list(itype)) // ' and ' // &
                             TRIM(atom_type_list(jtype)) // ' are repeated'
                CALL Clean_Abort(err_msg,'Create_Nonbond_Table')
             END IF
           END DO CustomLOOP
 
           ! Check that all expected atomtypes have been found
           DO itype = 1, nbr_atomtypes
             DO jtype = itype+1, nbr_atomtypes
        
               ! Check for errors
               IF (vdw_param_set(itype,jtype) == 0) THEN
                 err_msg = ''
                 err_msg(1) = 'Custom parameters not found for atom types: ' // &
                              TRIM(atom_type_list(itype)) // " " // &
                              TRIM(atom_type_list(jtype))
                 CALL Clean_Abort(err_msg,'Create_Nonbond_Table')
               END IF
            
             END DO
           END DO   
           EXIT

         ELSE IF (line_string(1:3) == 'END' .OR. line_nbr > 10000 ) THEN
                 
           err_msg = ""
           err_msg(1) = 'Section "# Mixing_Rule" missing from input file'
           CALL Clean_Abort(err_msg,'Create_Nonbond_Table')

         END IF
       END DO ! Loop over inputfile lines

     END IF ! mix_rule
  END IF ! nbr_atomtypes > 1

END SUBROUTINE Create_Nonbond_Table
