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
  REAL(DP), DIMENSION(max_nonbond_params) :: temp_param_i, temp_param_j
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

     ENDDO

  ENDDO

  
! ite the number of different atom types to the screen and logfile

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

  IF (ALLOCATED(temp_atomtypes)) DEALLOCATE(temp_atomtypes)

  ! Create a character array containing the names of each unique atom type, with the index equal
  ! to the atom type number

  ALLOCATE(atom_type_list(nbr_atomtypes), Stat=AllocateStatus)

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
  
  atom_type_list = ""

  ! Now determine the set of vdw parameters for each type of interaction and load them into vdw_param_table
  ! This is a brute force search - but it is fast.

  IF (verbose_log) THEN
     WRITE(logunit,*) 
     WRITE(logunit,*) 'Creating VDW interaction table'
     WRITE(logunit,'(X,A79)') '-------------------------------------------------------------------------------'
     WRITE(logunit,'(X,A,T25,A)') 'Mixing rule used is:', mix_rule

  ! Write header for logfile output. Specific for the vdw style
     IF (int_vdw_style(1) == vdw_lj) THEN
        WRITE(logunit,'(X,A6,5X,A6,2X,A12,X,A12)') 'Atom 1','Atom 2', 'epsilon', 'sigma'
        WRITE(logunit,'(X,6X,5X,6X,2X,A12,X,A12)') 'amu A^2/ps^2', 'Ang'
     ELSEIF (int_vdw_style(1) == vdw_mie) THEN
        WRITE(logunit,'(X,A6,5X,A6,2X,A12,X,A12,X,A12,X,A12)') 'Atom 1','Atom 2', 'epsilon', 'sigma', &
              'rep-expt', 'disp-expt'
        WRITE(logunit,'(X,6X,5X,6X,2X,A12,X,A12)') 'amu A^2/ps^2', 'Ang'
     ENDIF
  END IF

  DO itype = 1, nbr_atomtypes

    DO jtype = 1, nbr_atomtypes

      ! Flags that are tripped when a particular atomtype parameter set is located
      iset = 0
      jset = 0

      DO is = 1, nspecies

         DO ia = 1, natoms(is)

            IF (iset == 0) THEN
               ! Search for atomtype i because the parameters have not been found

               IF (nonbond_list(ia,is)%atom_type_number == itype) THEN
                  ! This type of atom has been found. Grab its vdw parameters

                  DO k=1, nbr_vdw_params(is)
                     temp_param_i(k) = nonbond_list(ia,is)%vdw_param(k)
                  ENDDO

                  ! This atomtype parameters are now determined. Trip flag
                  iset = 1
                
                  atom_type_list(itype) = nonbond_list(ia,is)%atom_name
  
               ENDIF

            ENDIF

            IF(jset == 0) THEN

            ! Now search for parameters of type j

               IF (nonbond_list(ia,is)%atom_type_number == jtype) THEN

                  DO k=1, nbr_vdw_params(is)
                     temp_param_j(k) = nonbond_list(ia,is)%vdw_param(k)
                  ENDDO
                
                  ! This atom type has been located
                  jset = 1

                  atom_type_list(jtype) = nonbond_list(ia,is)%atom_name

               ENDIF

            ENDIF

         ENDDO

      ENDDO

      ! Found i and j
      IF (mix_rule /= 'custom') THEN

        IF (int_vdw_style(1) == vdw_lj) THEN
           ! There are two vdw parameters that need mixing

           ! Set LJ epsilon
           IF ( (temp_param_i(1) <= tiny_number) .OR. (temp_param_j(1) <= tiny_number) ) THEN
              vdw_param1_table(itype,jtype) = 0.0_DP
           ! for parameters with zero, avoid overflow and set to zero

           ELSE

              ! Use specified mixing rule

           ! LB mixing rule: epsij = (epsi * epsj)^(1/2)
           ! geometric mixing rule: epsij = (epsi * epsj)^(1/2)
              vdw_param1_table(itype,jtype) = dsqrt(temp_param_i(1)*temp_param_j(1))
            
           ENDIF

           ! Set LJ sigma
           IF ( (temp_param_i(2) <= tiny_number) .OR. (temp_param_j(2) <= tiny_number) ) THEN
              vdw_param2_table(itype,jtype) = 0.0_DP

           ELSE

           ! LB mixing rule: sigmaij = 1/2 (sigmai + sigmaj)
              IF (mix_rule == 'LB') &
                   vdw_param2_table(itype,jtype) = (temp_param_i(2) + temp_param_j(2)) * 0.5
           ! geometric mixing rule: sigmaij = (sigmai * sigmaj)^(1/2)
              IF (mix_rule == 'geometric') &
                   vdw_param2_table(itype,jtype) = dsqrt(temp_param_i(2) * temp_param_j(2))

           ENDIF

           IF (verbose_log) THEN
              ! Report parameters to logfile.

             WRITE(logunit,'(X,A6,5X,A6,2X,F12.4,X,F12.4)') &
                  atom_type_list(itype), atom_type_list(jtype), &
                  vdw_param1_table(itype,jtype), vdw_param2_table(itype,jtype)

           ENDIF

        ELSEIF (int_vdw_style(1) == vdw_mie) THEN
           ! There are two vdw parameters that need mixing and two parameters
           ! that must be identical
          IF (temp_param_i(3) == temp_param_j(3) .AND. temp_param_i(4) == temp_param_j(4)) THEN

           ! Set epsilon
           IF ( (temp_param_i(1) <= tiny_number) .OR. (temp_param_j(1) <= tiny_number) ) THEN
              vdw_param1_table(itype,jtype) = 0.0_DP
           ELSE
              ! LB mixing rule: epsij = (epsi * epsj)^(1/2)
              ! geometric mixing rule: epsij = (epsi * epsj)^(1/2)
              vdw_param1_table(itype,jtype) = dsqrt(temp_param_i(1)*temp_param_j(1))
           ENDIF

           ! Set sigma
           IF ( (temp_param_i(2) <= tiny_number) .OR. (temp_param_j(2) <= tiny_number) ) THEN
              vdw_param2_table(itype,jtype) = 0.0_DP
           ELSE
              ! LB mixing rule: sigmaij = 1/2 (sigmai + sigmaj)
              IF (mix_rule == 'LB') &
                   vdw_param2_table(itype,jtype) = (temp_param_i(2) + temp_param_j(2)) * 0.5
              ! geometric mixing rule: sigmaij = (sigmai * sigmaj)^(1/2)
              IF (mix_rule == 'geometric') &
                   vdw_param2_table(itype,jtype) = dsqrt(temp_param_i(2) * temp_param_j(2))
           ENDIF

           ! Exponents are the same
           vdw_param3_table(itype,jtype) = temp_param_i(3)
           vdw_param4_table(itype,jtype) = temp_param_i(4)

           IF (verbose_log) THEN
              ! Report parameters to logfile.

             WRITE(logunit,'(X,A6,5X,A6,2X,F12.4,X,F12.4,X,F12.4,X,F12.4)') &
                  atom_type_list(itype), atom_type_list(jtype), &
                  vdw_param1_table(itype,jtype), vdw_param2_table(itype,jtype), &
                  vdw_param3_table(itype,jtype), vdw_param4_table(itype,jtype)


           ENDIF

          ELSE
            err_msg = ""
            err_msg(1) = "Cross interactions for Mie potential must have same exponents"
            CALL Clean_Abort(err_msg,'Create_Nonbond_Table')
          END IF
        ENDIF

      END IF

    ENDDO

  ENDDO

  ! atom types i and j have been found and vdw parameters loaded into temporary arrays.
  ! Apply mixing rules and load vdw_table

  ! Adapt this to other potential types and mixing rules. Custom mixing rules can 
  ! be created manually and then this routine should be bypassed.

  IF (mix_rule == 'custom') THEN

    REWIND(inputunit)

    ierr = 0
    line_nbr = 0

    ALLOCATE(vdw_param_set(nbr_atomtypes,nbr_atomtypes), Stat=AllocateStatus)
    IF (AllocateStatus .NE. 0) THEN
       err_msg = ''
       err_msg(1) = ' ERROR: Not enough memory for vdw_param_set'
       CALL Clean_Abort(err_msg,'create_nonbond_table')
    END IF
    vdw_param_set = 0

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
          IF (nbr_entries < 4) THEN
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
             err_msg = ''
             err_msg(1) = 'Atom type ' // TRIM(line_array(1)) // ' on line ' // &
                          TRIM(Int_To_String(line_nbr)) // ' of input file not found'
             CALL Clean_Abort(err_msg,'Create_Nonbond_Table')
          END IF

          ! find jtype
          jset = 0
          DO jtype = 1, nbr_atomtypes
             IF (atom_type_list(jtype) == line_array(2)) THEN
                jset = 1
                EXIT
             END IF
          END DO
          IF (jset == 0) THEN
             err_msg = ''
             err_msg(1) = 'Atom type ' // TRIM(line_array(2)) // ' on line ' // &
                          TRIM(Int_To_String(line_nbr)) // ' of input file not found'
             CALL Clean_Abort(err_msg,'Create_Nonbond_Table')
          END IF

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
        END DO CustomLOOP
 
        ! Check that all expected atomtypes have been found
        DO itype = 1, nbr_atomtypes
          DO jtype = itype, nbr_atomtypes
     
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
    END DO
  END IF

  IF (verbose_log) THEN
     WRITE(logunit,'(X,A79)') '-------------------------------------------------------------------------------'
  END IF

END SUBROUTINE Create_Nonbond_Table
