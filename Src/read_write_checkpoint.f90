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

MODULE Read_Write_Checkpoint
  !************************************************************************
  ! The module contains two subroutines
  ! Read a check point file when a simulation is restarted from a checkpoint
  ! file. Writes this checkpoint file periodically in a simulation.
  ! Note that any changes made in generating a checkpoint must mirror changes
  ! in the reading subroutine. This will be the case when additional information
  ! is written for various ensembles.
  !
  ! Revision History:
  ! 12/10/13  :: Beta version
  ! 08/12/26 (EJM) : Init/Maybe/Write_Log_Progress; reformat acceptance, move-
  !                  width, and subroutine-time logfile sections (logfile redesign)
  ! 08/13/26 (EJM) : Check_Restart_H_Consistency for read_config vs companion .H
  !**************************************************************************
  USE Global_Variables
  USE File_Names
  USE Simulation_Properties
  USE Random_Generators, ONLY : s1,s2,s3,s4,s5, rranf
  USE IO_Utilities
  USE Internal_Coordinate_Routines

  IMPLICIT NONE

  ! Next 10% progress milestone index (1..9); 10 means disabled / finished
  INTEGER :: next_log_progress_k = 1

CONTAINS

  SUBROUTINE Write_Checkpoint

    INTEGER :: ibox, is, ii, jj, im, this_im, ia, nmolecules_is, this_box
    INTEGER :: total_molecules_is, this_unit, position

    LOGICAL :: lopen
    REAL(DP) :: this_lambda = 1.0_DP

    INQUIRE(file=checkpointfile,opened=lopen)
    IF (lopen) INQUIRE(file=checkpointfile, number = this_unit)
    IF (lopen) CLOSE(unit=this_unit)

    OPEN(unit=chkptunit,file=checkpointfile)
    ! Let us write all the counters

    WRITE(chkptunit,*) '********* Translation, rotation, dihedral, angle distortion ******'

    DO ibox = 1, nbr_boxes
       DO is = 1, nspecies
          WRITE(chkptunit,'(5(I10,1x))') is, ntrials(is,ibox)%displacement, &
               ntrials(is,ibox)%rotation, ntrials(is,ibox)%dihedral, &
               ntrials(is,ibox)%angle
          WRITE(chkptunit,'(5(I10,1x))') is, nsuccess(is,ibox)%displacement, &
               nsuccess(is,ibox)%rotation, nsuccess(is,ibox)%dihedral, &
               nsuccess(is,ibox)%angle
          WRITE(chkptunit,'(E24.15)',ADVANCE='NO') max_disp(is,ibox)
          WRITE(chkptunit,'(E24.15)',ADVANCE='NO') max_rot(is,ibox)
          WRITE(chkptunit,'(E24.15)') species_list(is)%max_torsion

       END DO

       IF ( int_sim_type == sim_npt .OR. int_sim_type == sim_gemc .OR. &
            int_sim_type == sim_gemc_npt ) THEN
          WRITE(chkptunit,*) nvol_success(ibox), nvolumes(ibox)
       END IF
    END DO

    WRITE(chkptunit,*) '********** # of MC steps *********'
    WRITE(chkptunit,*) i_mcstep
    WRITE(chkptunit,*) '******** Box info ***********'

    DO ibox = 1, nbr_boxes
       WRITE(chkptunit,*) tot_trials(ibox)
       WRITE(chkptunit,*) box_list(ibox)%volume
       WRITE(chkptunit,*) box_list(ibox)%box_shape
       DO ii = 1, 3
          WRITE(chkptunit,'(3(F10.4,1X))') (box_list(ibox)%length(ii,jj), jj=1,3)
       END DO

       !--- inverse length
       DO ii = 1, 3
          WRITE(chkptunit,'(3(E12.5,1X))') (box_list(ibox)%length_inv(ii,jj), jj=1,3)
       END DO

       IF ( int_sim_type == sim_npt .OR. int_sim_type == sim_gemc .OR. &
            int_sim_type == sim_gemc_npt )  THEN

         WRITE(chkptunit,*) box_list(ibox)%dv_max

       END IF

    END DO

    WRITE(chkptunit,*) '**** SEEDS *******'
    WRITE(chkptunit,*) s1,s2,s3,s4,s5

    WRITE(chkptunit,*) '******* Info for total number of molecules'
    ! write number of molecules of each of the species
    DO is = 1, nspecies
       WRITE(chkptunit,*) is, SUM(nmols(is,1:nbr_boxes))
    END DO


    WRITE(chkptunit,*) '******* Writing coordinates for all the boxes'

    DO is = 1, nspecies
       DO ibox = 1, nbr_boxes
          DO im = 1, nmols(is,ibox)

             this_im = locate(im,is,ibox)
             this_box = molecule_list(this_im,is)%which_box

             DO ia = 1, natoms(is)
!                WRITE(chkptunit,'(A,T10,3(F15.10,1X),T70,I3)') nonbond_list(ia,is)%element, &
                WRITE(chkptunit,*) nonbond_list(ia,is)%element, &
                     atom_list(ia,this_im,is)%rxp, &
                     atom_list(ia,this_im,is)%ryp, &
                     atom_list(ia,this_im,is)%rzp, &
                     this_box
             END DO

          END DO
       END DO
    END DO

    CLOSE(unit=chkptunit)

  END SUBROUTINE Write_Checkpoint

!*******************************************************************************

SUBROUTINE Read_Checkpoint

    INTEGER :: ibox, is, ii, jj, im, this_im, ia, nmolecules_is, this_box
    INTEGER :: sp_nmoltotal(nspecies)
    INTEGER :: this_species, nfrac_global, i, this_rxnum, j, m, alive
    INTEGER :: this_unit

    INTEGER, DIMENSION(:), ALLOCATABLE :: total_molecules, n_int

    REAL(DP) :: this_lambda = 1.0_DP
    REAL(DP) :: E_self, xcom_old, ycom_old, zcom_old
    REAL(DP) :: xcom_new, ycom_new, zcom_new

    LOGICAL :: overlap
    LOGICAL :: lopen

    TYPE(Energy_Class) :: inrg

    ALLOCATE(total_molecules(nspecies))
    ALLOCATE(n_int(nspecies))
    IF(.NOT. ALLOCATED(ntrials)) ALLOCATE(ntrials(nspecies,nbr_boxes))
    IF(.NOT. ALLOCATED(tot_trials)) ALLOCATE(tot_trials(nbr_boxes))

    INQUIRE(file=restart_file,opened=lopen)
    IF (lopen) INQUIRE(file=restart_file, number = this_unit)
    IF (lopen) CLOSE(unit=this_unit)

    OPEN(unit=restartunit,file=restart_file)
    ! Let us read all the counters and count the number of molecules of
    ! each of the species in all the boxes
    nmols(:,:) = 0
    n_int(:) = 0

    READ(restartunit,*)

    DO ibox = 1, nbr_boxes
       WRITE(logunit,'(X,A)') 'Reading move parameters for box ' // TRIM(Int_To_String(ibox))

       DO is = 1, nspecies
          WRITE(logunit,'(X,A)') 'Reading move parameters for species ' // TRIM(Int_To_String(is))

          ! read information only if start_type == checkpoint

          READ(restartunit,'(5(I10,1x))') this_species, &
               ntrials(is,ibox)%displacement, &
               ntrials(is,ibox)%rotation, &
               ntrials(is,ibox)%dihedral, &
               ntrials(is,ibox)%angle

          READ(restartunit,'(5(I10,1x))') this_species, &
               nsuccess(is,ibox)%displacement, &
               nsuccess(is,ibox)%rotation, &
               nsuccess(is,ibox)%dihedral, &
               nsuccess(is,ibox)%angle

          READ(restartunit,'(3(E24.15))') max_disp(is,ibox), &
               max_rot(is,ibox), species_list(is)%max_torsion

          IF (prob_trans > 0.0_DP) THEN
            WRITE(logunit,'(2X,A,T24,F9.5)') 'max displacement', max_disp(is,ibox)
          END IF
          IF (prob_rot > 0.0_DP) THEN
            WRITE(logunit,'(2X,A,T24,F9.5)') 'max rotation', max_rot(is,ibox)
          END IF
          IF (prob_torsion > 0.0_DP) THEN
            WRITE(logunit,'(2X,A,T24,F9.5)') 'max dihedral change', species_list(is)%max_torsion
          END IF

       END DO

       IF ( int_sim_type == sim_npt .OR. int_sim_type == sim_gemc .OR. &
            int_sim_type == sim_gemc_npt ) THEN

          READ(restartunit,*) nvol_success(ibox), nvolumes(ibox)

       END IF
    END DO
    WRITE(logunit,*) 'Species move info read successfully'

    READ(restartunit,*)
    READ(restartunit,*) initial_mcstep
    READ(restartunit,*)
    WRITE(logunit,*) 'Initial MC step is ' // TRIM(Int_To_String(initial_mcstep))

    DO ibox = 1, nbr_boxes
       WRITE(logunit,'(X,A)') 'Reading info for box ' // TRIM(Int_To_String(ibox))

       READ(restartunit,*) tot_trials(ibox)
       READ(restartunit,*) box_list(ibox)%volume
       READ(restartunit,*) box_list(ibox)%box_shape

       DO ii = 1, 3
          READ(restartunit,*) (box_list(ibox)%length(ii,jj), jj=1,3)
       END DO

       !--- inverse length
       DO ii = 1, 3
          READ(restartunit,*) (box_list(ibox)%length_inv(ii,jj), jj=1,3)
       END DO

       CALL Compute_Cell_Dimensions(ibox)

       IF ( int_sim_type == sim_npt .OR. int_sim_type == sim_gemc .OR. &
            int_sim_type == sim_gemc_npt ) THEN

          READ(restartunit,*) box_list(ibox)%dv_max
          WRITE(logunit,'(2X,A,T24,F9.0)') 'max volume change', box_list(ibox)%dv_max

       END IF

    END DO
    WRITE(logunit,*) 'Box info read successfully'

    READ(restartunit,*)
    READ(restartunit,*) s1,s2,s3,s4,s5
    WRITE(logunit,*) 'Seed info read successfully'

    ! read total number of molecules of each of the species
    READ(restartunit,*)
    DO is = 1, nspecies
       READ(restartunit,*) this_species, sp_nmoltotal(is)
       IF (sp_nmoltotal(is) > 0) species_list(is)%l_solvent = .TRUE.
    END DO
    WRITE(logunit,*) 'Number of molecules read successfully'

    READ(restartunit,*)

    DO is = 1, nspecies

       DO im = 1, sp_nmoltotal(is)

          ! provide a linked number to this molecule
          molecule_list(im,is)%live = .TRUE.

          ! By default make all the molecules as integer molecules
          molecule_list(im,is)%molecule_type = int_normal

          DO ia = 1, natoms(is)
             READ(restartunit,*)nonbond_list(ia,is)%element, &
                  atom_list(ia,im,is)%rxp, &
                  atom_list(ia,im,is)%ryp, &
                  atom_list(ia,im,is)%rzp, &
                  this_box
             ! set exist flags for this atom
             atom_list(ia,im,is)%exist = .TRUE.
          END DO

          ! assign the molecule to this box
          molecule_list(im,is)%which_box = this_box
          nmols(is,this_box) = nmols(is,this_box) + 1
          locate(nmols(is,this_box),is,this_box) = SUM(nmols(is,1:nbr_boxes))

          molecule_list(im,is)%frac = this_lambda

       END DO
    END DO
    WRITE(logunit,*) 'Configuration read successfully'

    DO is = 1, nspecies
       IF(sp_nmoltotal(is) .LT. max_molecules(is)) THEN
          DO im = sp_nmoltotal(is)+1,max_molecules(is)
             molecule_list(im,is)%live = .FALSE.
             molecule_list(im,is)%frac = 1.0_DP
             molecule_list(im,is)%molecule_type = int_normal
             molecule_list(im,is)%which_box = 0
          END DO
       END IF
    END DO


    CALL Get_Internal_Coords

    ! Calculate COM and distance of the atom farthest to the COM.

    DO ibox = 1, nbr_boxes
       DO is = 1, nspecies
          DO im = 1, nmols(is,ibox)
             this_im = locate(im,is,ibox)
             IF( .NOT. molecule_list(this_im,is)%live) CYCLE
             ! Now let us ensure that the molecular COM is inside the central simulation box
             !
             CALL Get_COM(this_im,is)

             xcom_old = molecule_list(this_im,is)%xcom
             ycom_old = molecule_list(this_im,is)%ycom
             zcom_old = molecule_list(this_im,is)%zcom

             ! Apply PBC

             this_box = molecule_list(this_im,is)%which_box

             IF (l_cubic(this_box)) THEN

                CALL Apply_PBC_Anint(this_box,xcom_old,ycom_old,zcom_old, &
                     xcom_new, ycom_new, zcom_new)

!!$                IF (this_box == 2) THEN
!!$                   write(203,*) atom_list(1,this_im,is)%rxp, atom_list(1,this_im,is)%ryp, &
!!$                        atom_list(1,this_im,is)%rzp
!!$                END IF
!!$                write(*,*) 'cubic'

             ELSE

                CALL Minimum_Image_Separation(this_box,xcom_old,ycom_old,zcom_old, &
                     xcom_new, ycom_new, zcom_new)

!                write(*,*) 'minimum'

             END IF

             ! COM in the central simulation box

             molecule_list(this_im,is)%xcom = xcom_new
             molecule_list(this_im,is)%ycom = ycom_new
             molecule_list(this_im,is)%zcom = zcom_new

             ! displace atomic coordinates

             atom_list(1:natoms(is),this_im,is)%rxp = atom_list(1:natoms(is),this_im,is)%rxp + &
                  xcom_new - xcom_old
             atom_list(1:natoms(is),this_im,is)%ryp = atom_list(1:natoms(is),this_im,is)%ryp + &
                  ycom_new - ycom_old
             atom_list(1:natoms(is),this_im,is)%rzp = atom_list(1:natoms(is),this_im,is)%rzp + &
                  zcom_new - zcom_old

             CALL Compute_Max_Com_Distance(this_im,is)
          END DO
       END DO
    END DO


    DO ibox = 1, nbr_boxes
       IF(int_vdw_sum_style(ibox) == vdw_cut_tail) CALL Compute_Beads(ibox)
    END DO

    IF(ALLOCATED(total_molecules)) DEALLOCATE(total_molecules)

  END SUBROUTINE Read_Checkpoint

!*******************************************************************************

  SUBROUTINE Read_Config(ibox)
    !***************************************************************************
    ! The subroutine reads in a configuration to start a new simulation run.
    ! The format of the input
    ! file is identical to the checkpoint file in terms of atomic coordinates.
    !
    !
    !****************************************************************************

    IMPLICIT NONE

    ! Input
    INTEGER :: ibox

    ! Local
    INTEGER :: is, im, ia, nstart, nend, this_im, mols_this, nfrac_global, i, alive, j, this_rxnum
    INTEGER :: alive_new, counter, m, i_lambda, locate_base

    REAL(DP) :: this_lambda
    REAL(DP) :: E_recip, E_self, E_intra
    REAL(DP) :: E_old, xcom_old, ycom_old, zcom_old
    REAL(DP) :: xcom_new, ycom_new, zcom_new
    LOGICAL :: overlap

    Type(Energy_Class) :: inrg

    WRITE(logunit,*) 'Reading configuration for box', ibox

    OPEN(unit = old_config_unit,file=old_config_file(ibox))

    ! *.xyz format has two header lines
    READ(old_config_unit,*)
    READ(old_config_unit,*)

    ! Read in the coordinates of the molecules
    DO is = 1, nspecies

       WRITE(logunit,*) 'Reading ', nmols_to_read(is,ibox), &
          ' molecules of species ', is

       ! For box 2, need to continue numbering where box 1 left off
       locate_base = SUM(nmols(is,1:nbr_boxes))

       DO im = 1, nmols_to_read(is,ibox)
          ! provide a linked number to the molecule
          locate(im,is,ibox) = im + locate_base
          this_im = locate(im,is,ibox)
          molecule_list(this_im,is)%live = .TRUE.
          this_lambda = 1.0_DP
          ! By default all the molecules are normal
          molecule_list(this_im,is)%molecule_type = int_normal

          DO ia = 1, natoms(is)

             READ(old_config_unit,*)nonbond_list(ia,is)%element, &
                  atom_list(ia,this_im,is)%rxp, &
                  atom_list(ia,this_im,is)%ryp, &
                  atom_list(ia,this_im,is)%rzp
             ! set the frac and exist flags for this atom
             molecule_list(this_im,is)%frac = this_lambda
             atom_list(ia,this_im,is)%exist = .TRUE.

          END DO

          ! assign the molecule the box id

          molecule_list(this_im,is)%which_box = ibox
     ! Now let us ensure that the molecular COM is inside the central simulation box
       !
          CALL Get_COM(this_im,is)

          xcom_old = molecule_list(this_im,is)%xcom
          ycom_old = molecule_list(this_im,is)%ycom
          zcom_old = molecule_list(this_im,is)%zcom

          ! Apply PBC
          IF (l_cubic(ibox)) THEN
             CALL Apply_PBC_Anint(ibox,xcom_old,ycom_old,zcom_old, &
                  xcom_new, ycom_new, zcom_new)
          ELSE
             CALL Minimum_Image_Separation(ibox,xcom_old,ycom_old,zcom_old, &
                  xcom_new, ycom_new, zcom_new)
          END IF

          ! COM in the central simulation box
          molecule_list(this_im,is)%xcom = xcom_new
          molecule_list(this_im,is)%ycom = ycom_new
          molecule_list(this_im,is)%zcom = zcom_new

          ! COM in the central simulation box
          molecule_list(this_im,is)%xcom = xcom_new
          molecule_list(this_im,is)%ycom = ycom_new
          molecule_list(this_im,is)%zcom = zcom_new

          ! displace atomic coordinates
          atom_list(1:natoms(is),this_im,is)%rxp = &
               atom_list(1:natoms(is),this_im,is)%rxp + xcom_new - xcom_old
          atom_list(1:natoms(is),this_im,is)%ryp = &
               atom_list(1:natoms(is),this_im,is)%ryp + ycom_new - ycom_old
          atom_list(1:natoms(is),this_im,is)%rzp = &
               atom_list(1:natoms(is),this_im,is)%rzp + zcom_new - zcom_old

          nmols(is,ibox) = nmols(is,ibox) + 1

       END DO

    END DO

    CLOSE(unit = old_config_unit)

    CALL Get_Internal_Coords

    ! Calculate COM and distance of the atom farthest to the COM.

    DO is = 1, nspecies
       DO im = 1, nmols(is,ibox)
          this_im = locate(im,is,ibox)
          IF( .NOT. molecule_list(this_im,is)%live) CYCLE
          CALL Get_COM(this_im,is)
          CALL Compute_Max_Com_Distance(this_im,is)
       END DO
    END DO

    WRITE(logunit,*) 'Configurations read successfully'

    IF (int_vdw_sum_style(ibox) == vdw_cut_tail) CALL Compute_Beads(ibox)

    CALL Check_Restart_H_Consistency(ibox)

  END SUBROUTINE Read_Config
!*******************************************************************************

  SUBROUTINE Check_Restart_H_Consistency(ibox)
    !***************************************************************************
    ! If a companion .H exists for the read_config XYZ (same stem, .xyz → .H),
    ! compare volume / cell matrix to Box_Info and molecule counts to
    ! nmols_to_read. Abort on mismatch. Missing companion → log and return.
    !
    ! 08/13/26 (EJM) : Equil → production restart safety check
    !***************************************************************************

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ibox

    CHARACTER(FILENAME_LEN) :: companion_h, cfg
    CHARACTER(32) :: s_vol_inp, s_vol_file, s_len_inp, s_len_file
    INTEGER :: ios, is, is_file, n_spec_file, n_mol_file, ii, jj
    REAL(DP) :: vol_file, length_file(3,3), vol_tol, len_tol
    LOGICAL :: exists

    cfg = TRIM(ADJUSTL(old_config_file(ibox)))
    companion_h = cfg

    ! Derive companion: replace trailing .xyz / .XYZ with .H
    IF (LEN_TRIM(cfg) >= 4) THEN
       IF (cfg(LEN_TRIM(cfg)-3:LEN_TRIM(cfg)) == '.xyz' .OR. &
           cfg(LEN_TRIM(cfg)-3:LEN_TRIM(cfg)) == '.XYZ') THEN
          companion_h = cfg(1:LEN_TRIM(cfg)-4) // '.H'
       ELSE
          WRITE(logunit,'(A)') '  No .xyz suffix on config file; skipping restart.H check'
          RETURN
       END IF
    ELSE
       RETURN
    END IF

    INQUIRE(file=TRIM(companion_h), exist=exists)
    IF (.NOT. exists) THEN
       WRITE(logunit,'(A,A)') '  No companion restart H (skip box/N check): ', &
            TRIM(companion_h)
       RETURN
    END IF

    WRITE(logunit,'(A,A)') '  Checking Box_Info / molecule counts vs ', TRIM(companion_h)

    OPEN(unit=old_config_unit, file=TRIM(companion_h), status='OLD', &
         action='READ', iostat=ios)
    IF (ios /= 0) THEN
       err_msg = ''
       err_msg(1) = 'Unable to open companion H file: ' // TRIM(companion_h)
       CALL Clean_Abort(err_msg, 'Check_Restart_H_Consistency')
    END IF

    READ(old_config_unit, *, iostat=ios) vol_file
    IF (ios /= 0) THEN
       err_msg = ''
       err_msg(1) = 'Error reading volume from ' // TRIM(companion_h)
       CALL Clean_Abort(err_msg, 'Check_Restart_H_Consistency')
    END IF

    DO ii = 1, 3
       READ(old_config_unit, *, iostat=ios) (length_file(ii,jj), jj=1,3)
       IF (ios /= 0) THEN
          err_msg = ''
          err_msg(1) = 'Error reading cell matrix from ' // TRIM(companion_h)
          CALL Clean_Abort(err_msg, 'Check_Restart_H_Consistency')
       END IF
    END DO

    READ(old_config_unit, *, iostat=ios)  ! blank line
    READ(old_config_unit, *, iostat=ios) n_spec_file
    IF (ios /= 0) THEN
       err_msg = ''
       err_msg(1) = 'Error reading nspecies from ' // TRIM(companion_h)
       CALL Clean_Abort(err_msg, 'Check_Restart_H_Consistency')
    END IF

    IF (n_spec_file /= nspecies) THEN
       err_msg = ''
       err_msg(1) = 'nspecies in companion H does not match input'
       err_msg(2) = 'Companion H: ' // TRIM(Int_To_String(n_spec_file))
       err_msg(3) = 'Input:       ' // TRIM(Int_To_String(nspecies))
       err_msg(4) = 'File: ' // TRIM(companion_h)
       CALL Clean_Abort(err_msg, 'Check_Restart_H_Consistency')
    END IF

    DO is = 1, nspecies
       READ(old_config_unit, *, iostat=ios) is_file, n_mol_file
       IF (ios /= 0) THEN
          err_msg = ''
          err_msg(1) = 'Error reading molecule counts from ' // TRIM(companion_h)
          CALL Clean_Abort(err_msg, 'Check_Restart_H_Consistency')
       END IF
       IF (n_mol_file /= nmols_to_read(is,ibox)) THEN
          err_msg = ''
          err_msg(1) = 'Molecule count mismatch for species ' // TRIM(Int_To_String(is))
          err_msg(2) = 'Companion H nmols: ' // TRIM(Int_To_String(n_mol_file))
          err_msg(3) = 'read_config nmols:  ' // TRIM(Int_To_String(nmols_to_read(is,ibox)))
          err_msg(4) = 'File: ' // TRIM(companion_h)
          err_msg(5) = 'Fix # Start_Type counts or use the matching restart files.'
          CALL Clean_Abort(err_msg, 'Check_Restart_H_Consistency')
       END IF
    END DO

    CLOSE(unit=old_config_unit)

    ! Volume / cell vs Box_Info
    vol_tol = MAX(1.0e-6_DP * ABS(box_list(ibox)%volume), 1.0e-6_DP)
    IF (ABS(vol_file - box_list(ibox)%volume) > vol_tol) THEN
       WRITE(s_vol_inp,'(ES16.6)') box_list(ibox)%volume
       WRITE(s_vol_file,'(ES16.6)') vol_file
       err_msg = ''
       err_msg(1) = 'Box volume in # Box_Info does not match companion restart H'
       err_msg(2) = 'Box_Info volume:    ' // TRIM(ADJUSTL(s_vol_inp))
       err_msg(3) = 'Companion H volume: ' // TRIM(ADJUSTL(s_vol_file))
       err_msg(4) = 'File: ' // TRIM(companion_h)
       err_msg(5) = 'Copy final equilibration box into # Box_Info (NPT density).'
       CALL Clean_Abort(err_msg, 'Check_Restart_H_Consistency')
    END IF

    len_tol = 1.0e-6_DP
    DO ii = 1, 3
       DO jj = 1, 3
          IF (ABS(length_file(ii,jj) - box_list(ibox)%length(ii,jj)) > &
               MAX(len_tol, 1.0e-6_DP * ABS(box_list(ibox)%length(ii,jj)))) THEN
             WRITE(s_len_inp,'(ES16.6)') box_list(ibox)%length(ii,jj)
             WRITE(s_len_file,'(ES16.6)') length_file(ii,jj)
             err_msg = ''
             err_msg(1) = 'Box cell matrix in # Box_Info does not match companion restart H'
             err_msg(2) = 'Element (' // TRIM(Int_To_String(ii)) // ',' // &
                  TRIM(Int_To_String(jj)) // ')'
             err_msg(3) = 'Box_Info:    ' // TRIM(ADJUSTL(s_len_inp))
             err_msg(4) = 'Companion H: ' // TRIM(ADJUSTL(s_len_file))
             err_msg(5) = 'File: ' // TRIM(companion_h)
             CALL Clean_Abort(err_msg, 'Check_Restart_H_Consistency')
          END IF
       END DO
    END DO

    WRITE(logunit,'(A)') '  Companion restart H agrees with Box_Info and molecule counts'

  END SUBROUTINE Check_Restart_H_Consistency
!*******************************************************************************

SUBROUTINE Write_Trials_Success
  !*****************************************************************************
  ! This subroutine writes number of trials and acceptance of these trials at
  ! the end of a simulation

  IMPLICIT NONE

  INTEGER :: ibox, is, ifrag
  REAL(DP) :: x1, x2

  WRITE(logunit,*)
  WRITE(logunit,*)

  DO ibox = 1, nbr_boxes

     WRITE(logunit,'(A,I0)') '  box ', ibox
     WRITE(logunit,'(A)') '  -----'

     IF (nvolumes(ibox) /= 0 ) THEN
        WRITE(logunit,'(A20,2X,A10,2X,A10,2X,A10)') 'Move', 'Trials', 'Success', '% Success'
        WRITE(logunit,11) 'Volume', nvolumes(ibox), nvol_success(ibox), &
          100.0*dble(nvol_success(ibox))/dble(nvolumes(ibox))
     END IF

     DO is = 1, nspecies

        WRITE(logunit,*)
        WRITE(logunit,'(3X,A57)') '---------------------------------------------------------'
        WRITE(logunit,'(3X,A8,X,I2)') 'Species', is
        WRITE(logunit,*)
        WRITE(logunit,'(A20,2X,A10,2X,A10,2X,A10)') 'Move', 'Trials', 'Success', '% Success'

        ! translation

        IF (ntrials(is,ibox)%displacement /= 0 ) THEN

           WRITE(logunit,11) 'Translate', ntrials(is,ibox)%displacement, &
                nsuccess(is,ibox)%displacement, &
                100.0*dble(nsuccess(is,ibox)%displacement)/dble(ntrials(is,ibox)%displacement)
        END IF

        ! rotation

        IF (ntrials(is,ibox)%rotation /= 0 ) THEN

           WRITE(logunit,11) 'Rotate',  ntrials(is,ibox)%rotation, &
                nsuccess(is,ibox)%rotation, &
                100.0*dble(nsuccess(is,ibox)%rotation)/dble(ntrials(is,ibox)%rotation)

        END IF

        ! Angle

        IF (ntrials(is,ibox)%angle /=0 ) THEN

           WRITE(logunit,11) 'Angle',  ntrials(is,ibox)%angle, &
                nsuccess(is,ibox)%angle, &
                100.0*dble(nsuccess(is,ibox)%angle)/dble(ntrials(is,ibox)%angle)

        END IF


        ! Dihedral

        IF (ntrials(is,ibox)%dihedral /= 0 ) THEN

           WRITE(logunit,11) 'Dihedral', ntrials(is,ibox)%dihedral, &
                nsuccess(is,ibox)%dihedral, &
                100.0*dble(nsuccess(is,ibox)%dihedral)/dble(ntrials(is,ibox)%dihedral)

        END IF

        ! insertion

        IF (ntrials(is,ibox)%insertion /= 0 ) THEN

           WRITE(logunit,11) 'Insert',  ntrials(is,ibox)%insertion, &
                nsuccess(is,ibox)%insertion, &
                100.0*dble(nsuccess(is,ibox)%insertion)/dble(ntrials(is,ibox)%insertion)
        END IF

        ! deletion

        IF (ntrials(is,ibox)%deletion /= 0 ) THEN

           WRITE(logunit,11) 'Delete', ntrials(is,ibox)%deletion, &
                nsuccess(is,ibox)%deletion, &
                100.0*dble(nsuccess(is,ibox)%deletion)/dble(ntrials(is,ibox)%deletion)

        END IF

        ! identity switch

        IF (ntrials(is,ibox)%switch /= 0) THEN
           WRITE(logunit,11) 'Identity Switch', ntrials(is,ibox)%switch, &
              nsuccess(is,ibox)%switch, &
              100.0*dble(nsuccess(is,ibox)%switch)/dble(ntrials(is,ibox)%switch)

        END IF

        ! atom displacement

        IF (ntrials(is,ibox)%disp_atom /= 0 ) THEN

           WRITE(logunit,11) 'Atom Displacement', ntrials(is,ibox)%disp_atom, &
                nsuccess(is,ibox)%disp_atom, &
                100.0*dble(nsuccess(is,ibox)%disp_atom)/dble(ntrials(is,ibox)%disp_atom)

        END IF

        WRITE(logunit,'(3X,A57)') '---------------------------------------------------------'

        WRITE(logunit,*)
      END DO

   END DO

11 FORMAT(A20,2x,I10,2x,I10,2x,f10.2)
12 FORMAT(I20,2x,I10,2x,I10,2x,f10.2)

  IF (SUM(nfragments) .GT. 0) THEN
     IF (SUM(regrowth_trials(:,:)) .GT. 0) THEN

        WRITE(logunit,*)
        WRITE(logunit,'(A)') ' Fragment regrowth'
        WRITE(logunit,'(A)') ' -----------------'

        DO is = 1, nspecies

           IF (SUM(regrowth_trials(:,is)) .GT. 0) THEN

              WRITE(logunit,*)
              WRITE(logunit,'(3X,A57)') '---------------------------------------------------------'
              WRITE(logunit,'(3X,A8,X,I2)') 'Species', is
              WRITE(logunit,'(A20,2x,A10,2x,A10,2X,A10)') '#_Frags_Regrown', 'Trials', 'Success', '% Success'

              DO ifrag = 1, nfragments(is)

                 IF (regrowth_trials(ifrag,is) /= 0 ) THEN
                    WRITE(logunit,12) ifrag,regrowth_trials(ifrag,is), &
                         regrowth_success(ifrag,is), &
                         100.0_DP * dble(regrowth_success(ifrag,is))/dble(regrowth_trials(ifrag,is))
                 END IF
              END DO
              WRITE(logunit,'(3X,A57)') '---------------------------------------------------------'
           END IF

        END DO

     END IF
  END IF

END SUBROUTINE Write_Trials_Success

!*******************************************************************************
! Log progress snapshots (Cassandra V2 logfile redesign)
!*******************************************************************************

SUBROUTINE Init_Log_Progress
  ! Reset 10% progress milestones for the current driver run.
  ! Progress is measured over the remaining span (n_mcsteps - initial_mcstep).
  IMPLICIT NONE
  INTEGER :: span, k, thresh

  next_log_progress_k = 1
  IF (timed_run) THEN
     next_log_progress_k = 10
     RETURN
  END IF
  span = n_mcsteps - initial_mcstep
  IF (span < 10) THEN
     next_log_progress_k = 10
     RETURN
  END IF
  ! If restarting past some milestones of this span, skip those already reached
  DO k = 1, 9
     thresh = initial_mcstep + (k * span) / 10
     IF (initial_mcstep >= thresh) THEN
        next_log_progress_k = k + 1
     ELSE
        EXIT
     END IF
  END DO
END SUBROUTINE Init_Log_Progress

SUBROUTINE Maybe_Write_Log_Progress
  ! Call once per MC step from the drivers. Writes when i_mcstep first reaches
  ! each 10%, 20%, ..., 90% threshold of the planned step span.
  IMPLICIT NONE
  INTEGER :: span, thresh, pct

  IF (next_log_progress_k > 9) RETURN
  IF (timed_run) RETURN

  span = n_mcsteps - initial_mcstep
  IF (span < 10) RETURN

  DO WHILE (next_log_progress_k <= 9)
     thresh = initial_mcstep + (next_log_progress_k * span) / 10
     IF (i_mcstep < thresh) EXIT
     pct = next_log_progress_k * 10
     CALL Write_Log_Progress(pct)
     next_log_progress_k = next_log_progress_k + 1
  END DO
END SUBROUTINE Maybe_Write_Log_Progress

SUBROUTINE Write_Log_Progress(pct)
  ! Mid-run logfile snapshot: totals energy, N/V/density, acceptance, widths, times.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: pct
  INTEGER :: ibox, is, n_mol_box
  REAL(DP) :: e_ext, e_int, mass_density, inv_n

  WRITE(logunit,*)
  WRITE(logunit,'(A80)') '********************************************************************************'
  WRITE(logunit,'(A,I3,A,I12,A,I12)') ' Progress ', pct, '%   step ', i_mcstep, ' / ', n_mcsteps
  WRITE(logunit,'(A80)') '********************************************************************************'

  DO ibox = 1, nbr_boxes
     n_mol_box = 0
     DO is = 1, nspecies
        n_mol_box = n_mol_box + nmols(is, ibox)
     END DO

     mass_density = 0.0_DP
     IF (box_list(ibox)%volume > tiny_number) THEN
        DO is = 1, nspecies
           mass_density = mass_density + REAL(nmols(is,ibox),DP) * species_list(is)%molecular_weight
        END DO
        mass_density = mass_density / box_list(ibox)%volume * atomic_to_kgm3
     END IF

     WRITE(logunit,*)
     WRITE(logunit,'(A,I2,A,I8,A,ES14.6,A,A,F12.3,A)') &
          ' Box ', ibox, ': N = ', n_mol_box, &
          '   V = ', box_list(ibox)%volume, ' Ang^3', &
          '   density = ', mass_density, ' kg/m^3'
     WRITE(logunit,'(A)',ADVANCE='NO') '   N by species:'
     DO is = 1, nspecies
        WRITE(logunit,'(X,I0,A,I0)',ADVANCE='NO') is, '=', nmols(is,ibox)
     END DO
     WRITE(logunit,*)

     e_ext = energy(ibox)%total * atomic_to_kjmol
     IF (n_mol_box > 0) THEN
        inv_n = 1.0_DP / REAL(n_mol_box, DP)
        e_int = e_ext * inv_n
        WRITE(logunit,'(A,F16.3,A,F16.3,A)') &
             ' Total system energy: ', e_ext, ' kJ/mol   (intensive ', e_int, &
             ' kJ/mol/molecule)'
     ELSE
        WRITE(logunit,'(A,F16.3,A)') &
             ' Total system energy: ', e_ext, ' kJ/mol   (intensive n/a)'
     END IF
  END DO

  WRITE(logunit,*)
  WRITE(logunit,'(A)') ' Move acceptance (cumulative to this point)'
  CALL Write_Trials_Success

  WRITE(logunit,*)
  WRITE(logunit,'(A)') ' Current move widths'
  WRITE(logunit,'(A)') ' -------------------'
  DO ibox = 1, nbr_boxes
     IF (ALLOCATED(max_disp)) THEN
        DO is = 1, nspecies
           WRITE(logunit,'(A,I0,A,I0,A,F12.6,A)') &
                '  max_disp species ', is, ' box ', ibox, ': ', max_disp(is,ibox), ' Ang'
           WRITE(logunit,'(A,I0,A,I0,A,F12.6,A)') &
                '  max_rot  species ', is, ' box ', ibox, ': ', max_rot(is,ibox), ' rad'
        END DO
     END IF
     WRITE(logunit,'(A,I0,A,ES14.6,A)') &
          '  dv_max box ', ibox, ': ', box_list(ibox)%dv_max, ' Ang^3'
  END DO

  CALL Write_Subroutine_Times

END SUBROUTINE Write_Log_Progress

SUBROUTINE Write_Subroutine_Times

  IMPLICIT NONE

WRITE(logunit,*)
WRITE(logunit,'(A)') ' Subroutine times'
WRITE(logunit,'(A)') ' ----------------'


IF(movetime(imove_trans) .GT. 0.0_DP ) THEN

   IF(movetime(imove_trans) .LT. 60.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Translation time = ', movetime(imove_trans), ' secs.'
   ELSE IF(movetime(imove_trans) .LT. 3600.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Translation time = ', movetime(imove_trans) / 60.0_DP , ' mins.'
   ELSE
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Translation time = ', movetime(imove_trans) / 3600.0_DP , ' hrs.'
   END IF

END IF

IF(movetime(imove_rot) .GT. 0.0_DP ) THEN

   IF(movetime(imove_rot) .LT. 60.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Rotation time = ', movetime(imove_rot), ' secs.'
   ELSE IF(movetime(imove_rot) .LT. 3600.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Rotation time = ', movetime(imove_rot) / 60.0_DP , ' mins.'
   ELSE
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Rotation time = ', movetime(imove_rot) / 3600.0_DP , ' hrs.'
   END IF

END IF

IF(movetime(imove_dihedral) .GT. 0.0_DP ) THEN

   IF(movetime(imove_dihedral) .LT. 60.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Dihedral time = ', movetime(imove_dihedral), ' secs.'
   ELSE IF(movetime(imove_dihedral) .LT. 3600.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Dihedral time = ', movetime(imove_dihedral) / 60.0_DP , ' mins.'
   ELSE
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Dihedral time = ', movetime(imove_dihedral) / 3600.0_DP , ' hrs.'
   END IF

END IF

IF(movetime(imove_angle) .GT. 0.0_DP ) THEN

   IF(movetime(imove_angle) .LT. 60.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Angle change time = ', movetime(imove_angle), ' secs.'
   ELSE IF(movetime(imove_angle) .LT. 3600.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Angle change time = ', movetime(imove_angle) / 60.0_DP , ' mins.'
   ELSE
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Angle change time = ', movetime(imove_angle) / 3600.0_DP , ' hrs.'
   END IF

END IF

IF(movetime(imove_volume) .GT. 0.0_DP ) THEN

   IF(movetime(imove_volume) .LT. 60.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Volume change time = ', movetime(imove_volume), ' secs.'
   ELSE IF(movetime(imove_volume) .LT. 3600.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Volume change time = ', movetime(imove_volume) / 60.0_DP , ' mins.'
   ELSE
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Volume change time = ', movetime(imove_volume) / 3600.0_DP , ' hrs.'
   END IF

END IF

IF(movetime(imove_insert) .GT. 0.0_DP ) THEN

   IF(movetime(imove_insert) .LT. 60.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Insertion time = ', movetime(imove_insert), ' secs.'
   ELSE IF(movetime(imove_insert) .LT. 3600.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Insertion time = ', movetime(imove_insert) / 60.0_DP , ' mins.'
   ELSE
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Insertion time = ', movetime(imove_insert) / 3600.0_DP , ' hrs.'
   END IF

END IF

IF(movetime(imove_swap) .GT. 0.0_DP ) THEN

   IF(movetime(imove_swap) .LT. 60.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Swap time = ', movetime(imove_swap), ' secs.'
   ELSE IF(movetime(imove_swap) .LT. 3600.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Swap time = ', movetime(imove_swap) / 60.0_DP , ' mins.'
   ELSE
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Swap time = ', movetime(imove_swap) / 3600.0_DP , ' hrs.'
   END IF

END IF

IF(movetime(imove_delete) .GT. 0.0_DP ) THEN

   IF(movetime(imove_delete) .LT. 60.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Deletion time = ', movetime(imove_delete), ' secs.'
   ELSE IF(movetime(imove_delete) .LT. 3600.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Deletion time = ', movetime(imove_delete) / 60.0_DP , ' mins.'
   ELSE
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Deletion time = ', movetime(imove_delete) / 3600.0_DP , ' hrs.'
   END IF

END IF


IF(movetime(imove_regrowth) .GT. 0.0_DP ) THEN

   IF(movetime(imove_regrowth) .LT. 60.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Regrowth time = ', movetime(imove_regrowth), ' secs.'
   ELSE
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Regrowth time = ', movetime(imove_regrowth) / 3600.0_DP , ' hrs.'
   END IF

END IF

IF(movetime(imove_identity_switch) .GT. 0.0_DP ) THEN

   IF(movetime(imove_identity_switch) .LT. 60.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Regrowth time = ', movetime(imove_identity_switch), ' secs.'
   ELSE
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Regrowth time = ', movetime(imove_identity_switch) / 3600.0_DP , ' hrs.'
   END IF

END IF


IF(movetime(imove_widom) .GT. 0.0_DP ) THEN

   IF(movetime(imove_widom) .LT. 60.0_DP ) THEN
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Widom time = ', movetime(imove_widom), ' secs.'
   ELSE
        WRITE(logunit,'(1X,A,T25,F15.6,A)') &
       'Widom time = ', movetime(imove_widom) / 3600.0_DP , ' hrs.'
   END IF

END IF

END SUBROUTINE Write_Subroutine_Times

END MODULE Read_Write_Checkpoint

