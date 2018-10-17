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
  !**************************************************************************
  USE Global_Variables
  USE File_Names
  USE Simulation_Properties
  USE Random_Generators, ONLY : s1,s2,s3,s4,s5, rranf
  USE IO_Utilities
  
  IMPLICIT NONE

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

  END SUBROUTINE Read_Config
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

     WRITE(logunit,'(A27,X,I2)') 'Writing information for box', ibox
     WRITE(logunit,'(A80)') '********************************************************************************'

     IF (nvolumes(ibox) /= 0 ) THEN
        WRITE(logunit,'(A20,2X,A10,2X,A10,2X,A10)') 'Move', 'Trials', 'Success', '% Success'
        WRITE(logunit,11) 'Volume', nvolumes(ibox), nvol_success(ibox), &
          100.0*dble(nvol_success(ibox))/dble(nvolumes(ibox))
     END IF

     DO is = 1, nspecies
        
        WRITE(logunit,*) 
        WRITE(logunit,'(3X,A57)') '---------------------------------------------------------'
        WRITE(logunit,'(3X,A31,X,I2)') 'Writing information for species', is
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

        ! atom displacement

        IF (ntrials(is,ibox)%disp_atom /= 0 ) THEN
           
           WRITE(logunit,11) 'Atom Displacement', ntrials(is,ibox)%disp_atom, &
                nsuccess(is,ibox)%disp_atom, &
                100.0*dble(nsuccess(is,ibox)%disp_atom)/dble(ntrials(is,ibox)%disp_atom)

        END IF
        WRITE(logunit,'(3X,A57)') '---------------------------------------------------------'

        WRITE(logunit,*)
      END DO

      WRITE(logunit,'(A80)') '********************************************************************************'
        
   END DO

11 FORMAT(A20,2x,I10,2x,I10,2x,f10.2)
12 FORMAT(I20,2x,I10,2x,I10,2x,f10.2)

  IF (SUM(nfragments) .GT. 0) THEN
     IF (SUM(regrowth_trials(:,:)) .GT. 0) THEN

        WRITE(logunit,*)
        WRITE(logunit,'(A)') 'Writing information about fragments'
        WRITE(logunit,'(A80)') '********************************************************************************'
     
        DO is = 1, nspecies
           
           IF (SUM(regrowth_trials(:,is)) .GT. 0) THEN
 
              WRITE(logunit,*)
              WRITE(logunit,'(3X,A57)') '---------------------------------------------------------'
              WRITE(logunit,'(3X,A31,X,I2)') 'Writing information for species', is
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
        WRITE(logunit,'(A80)') '********************************************************************************'

     END IF
  END IF

END SUBROUTINE Write_Trials_Success

SUBROUTINE Write_Subroutine_Times

  IMPLICIT NONE

WRITE(logunit,*)
WRITE(logunit,*) 'Writing information about subroutine times'
WRITE(logunit,'(A80)') '********************************************************************************'


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

WRITE(logunit,'(A80)') '********************************************************************************'

END SUBROUTINE Write_Subroutine_Times

END MODULE Read_Write_Checkpoint

