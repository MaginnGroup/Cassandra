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

MODULE Read_Write_Checkpoint
  !************************************************************************
  ! The module contains two subroutines
  ! Read a check point file when a simulation is restarted from a checkpoint file
  ! Writes this checkpoint file periodically in a simulation.
  ! Note that any changes made in generating a checkpoint must mirror changes
  ! in the reading subroutine. This will be the case when additional information
  ! is written for various ensembles.
  !
  ! Revision History: 
  ! 12/10/13  :: Beta version 
  !**************************************************************************
  USE Run_Variables
  USE File_Names
  USE Simulation_Properties
  USE Random_Generators, ONLY : s1,s2,s3,s4,s5, rranf
  USE Energy_Routines, ONLY : Compute_Molecule_Energy
  
  IMPLICIT NONE

CONTAINS

  SUBROUTINE Write_Checkpoint(this_mc_step)

    INTEGER, INTENT(IN) :: this_mc_step

    INTEGER :: ibox, is, ii, jj, im, this_im, ia, nmolecules_is, this_box
    INTEGER :: total_molecules_is, this_unit
    
    LOGICAL :: lopen

    INQUIRE(file=checkpointfile,opened=lopen)
    IF (lopen) INQUIRE(file=checkpointfile, number = this_unit)
    IF (lopen) CLOSE(unit=this_unit)

    OPEN(unit=chkptunit,file=checkpointfile)
    ! Let us write all the counters

    WRITE(chkptunit,*) '********* Translation,rotation, dihedral, angle distortion ******'

    DO ibox = 1, nbr_boxes
       DO is = 1, nspecies
          WRITE(chkptunit,'(5(I10,1x))') is, ntrials(is,ibox)%displacement, &
               ntrials(is,ibox)%rotation, ntrials(is,ibox)%dihedral, &
               ntrials(is,ibox)%angle
          WRITE(chkptunit,'(5(I10,1x))') is, nsuccess(is,ibox)%displacement, &
               nsuccess(is,ibox)%rotation, nsuccess(is,ibox)%dihedral, &
               nsuccess(is,ibox)%angle
          WRITE(chkptunit,'(3(E24.15))') max_disp(is,ibox), max_rot(is,ibox), &
               species_list(is)%max_torsion
          
       END DO

       IF (int_sim_type == sim_npt .OR. int_sim_type == sim_gemc .OR. &
            int_sim_type == sim_gemc_npt) THEN
          WRITE(chkptunit,*) nvol_success(ibox), nvolumes(ibox)
       END IF
    END DO
    
    WRITE(chkptunit,*) '********** # of MC steps *********'
    WRITE(chkptunit,*) this_mc_step
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

       IF (int_sim_type == sim_npt .OR. int_sim_type == sim_gemc .OR. &
            int_sim_type == sim_gemc_npt)  THEN

             WRITE(chkptunit,*) box_list(ibox)%dv_max
          
       END IF
       
    END DO
    WRITE(chkptunit,*) '**** SEEDS *******'
    WRITE(chkptunit,*) s1,s2,s3,s4,s5
    
    WRITE(chkptunit,*) '******* Info for total number of molecules'
    ! write number of molecules of each of the species
    DO is = 1, nspecies
       total_molecules_is = 0
       DO ibox = 1, nbr_boxes
          CALL Get_Nmolecules_Species(ibox,is,nmolecules_is)
          total_molecules_is = total_molecules_is + nmolecules_is
       END DO
       WRITE(chkptunit,*) is,total_molecules_is
    END DO
    
    
    WRITE(chkptunit,*) '********Writing coordinates for all the boxes'
    
    DO is = 1, nspecies
       DO im = 1, nmolecules(is)
          
          this_im = locate(im,is)
          
          IF(molecule_list(this_im,is)%live) THEN
             this_box = molecule_list(this_im,is)%which_box

             DO ia = 1, natoms(is)
!                WRITE(chkptunit,'(A,T10,3(F15.10,1X),T70,I3)') nonbond_list(ia,is)%element, &
                WRITE(chkptunit,*) nonbond_list(ia,is)%element, &
                     atom_list(ia,this_im,is)%rxp, &
                     atom_list(ia,this_im,is)%ryp, &
                     atom_list(ia,this_im,is)%rzp, this_box
             END DO
          END IF
          
       END DO
    END DO

    CLOSE(unit=chkptunit)
    
  END SUBROUTINE Write_Checkpoint
!**************************************************************************************************

SUBROUTINE Read_Checkpoint

   INTEGER :: this_mc_step

    INTEGER :: ibox, is, ii, jj, im, this_im, ia, nmolecules_is, this_box, mols_this, sp_nmoltotal(nspecies)
    INTEGER :: this_species, nfrac_global, i, this_rxnum, j, m, alive
    INTEGER :: this_unit, i_lambda

    INTEGER, DIMENSION(:), ALLOCATABLE :: total_molecules, n_int

    REAL(DP) :: this_lambda, E_self, xcom_old, ycom_old, zcom_old
    REAL(DP) :: xcom_new, ycom_new, zcom_new

    LOGICAL :: f_checkpoint, f_read_old, overlap, cfc_defined
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

    f_checkpoint = .FALSE.
    f_read_old = .FALSE.
    IF (start_type == 'checkpoint') f_checkpoint = .TRUE.
    IF (start_type == 'read_old') f_read_old = .TRUE.

    DO ibox = 1, nbr_boxes

       DO is = 1, nspecies
          
          ! read information only if start_type == checkpoint
          
          IF (f_checkpoint) THEN

             READ(restartunit,'(5(I10,1x))') this_species, ntrials(is,ibox)%displacement, &
                  ntrials(is,ibox)%rotation, ntrials(is,ibox)%dihedral, &
                  ntrials(is,ibox)%angle

             READ(restartunit,'(5(I10,1x))') this_species, nsuccess(is,ibox)%displacement, &
                  nsuccess(is,ibox)%rotation, nsuccess(is,ibox)%dihedral, &
                  nsuccess(is,ibox)%angle
             READ(restartunit,'(3(E24.15))') max_disp(is,ibox), max_rot(is,ibox), &
                  species_list(is)%max_torsion

          ELSE IF (f_read_old) THEN

             READ(restartunit,*)
             READ(restartunit,*)
             READ(restartunit,*)

          END IF

       END DO
       
       IF (int_sim_type == sim_npt .OR. int_sim_type == sim_gemc .OR. &
            int_sim_type == sim_gemc_npt) THEN

          IF (f_checkpoint) THEN
             READ(restartunit,*) nvol_success(ibox), nvolumes(ibox)
          ELSE IF (f_read_old) THEN
             READ(restartunit,*)
          END IF

       END IF
    END DO
    
    READ(restartunit,*)
    READ(restartunit,*) this_mc_step
    READ(restartunit,*) 
    
    DO ibox = 1, nbr_boxes

       IF (f_checkpoint) THEN
          
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
                    
       ELSE IF (f_read_old) THEN

          READ(restartunit,*)
          READ(restartunit,*)
          READ(restartunit,*)

          ! skip box length and box length inverse information

          DO ii = 1, 3
             READ(restartunit,*)
          END DO

          DO ii = 1, 3
             READ(restartunit,*)
          END DO

       END IF

       
       IF (int_sim_type == sim_npt .OR. int_sim_type == sim_gemc .OR. &
            int_sim_type == sim_gemc_npt) THEN
          
          IF (f_checkpoint) THEN
                READ(restartunit,*) box_list(ibox)%dv_max          
          ELSE
             READ(restartunit,*)
          END IF

       END IF
       
    END DO
    
    READ(restartunit,*)
    IF (f_checkpoint) THEN
!       READ(restartunit,*) iseed1, iseed3
        READ(restartunit,*) s1,s2,s3,s4,s5
    ELSE
       READ(restartunit,*)
    END IF

    ! read total number of molecules of each of the species
    READ(restartunit,*)
    DO is = 1, nspecies
       READ(restartunit,*) this_species, sp_nmoltotal(is)
    END DO
    
    READ(restartunit,*)
    
    DO is = 1, nspecies

       mols_this = 0

       DO im = 1, sp_nmoltotal(is)
          
          ! provide a linked number to this molecule
          locate(im,is) = im + mols_this
          
          this_im = locate(im,is)
          
          molecule_list(this_im,is)%live = .TRUE.

          this_lambda = 1.0_DP
          
          ! By default make all the molecules as integer molecules
          
          molecule_list(this_im,is)%molecule_type = int_normal

          molecule_list(this_im,is)%cfc_lambda = this_lambda

          DO ia = 1, natoms(is)

             READ(restartunit,*)nonbond_list(ia,is)%element, &
                  atom_list(ia,this_im,is)%rxp, &
                  atom_list(ia,this_im,is)%ryp, &
                  atom_list(ia,this_im,is)%rzp 
             READ(restartunit,*) this_box
             ! set the cfc_lambda and exist flags for this atom
             atom_list(ia,this_im,is)%exist = .TRUE.
          END DO

          ! assign the box to this molecule
             
          molecule_list(this_im,is)%which_box = this_box
          nmols(is,this_box) = nmols(is,this_box) + 1
                
       END DO
    END DO

    DO is = 1, nspecies
       IF(sp_nmoltotal(is) .LT. nmolecules(is)) THEN
          DO im = sp_nmoltotal(is)+1,nmolecules(is)
             locate(im,is) = im
             molecule_list(im,is)%live = .FALSE.
             molecule_list(im,is)%cfc_lambda = 1.0_DP
             molecule_list(im,is)%molecule_type = int_normal
             molecule_list(im,is)%which_box = 0
          END DO
       END IF
    END DO

       
    CALL Get_Internal_Coords
    
    ! Calculate COM and distance of the atom farthest to the COM.
    
    DO is = 1, nspecies
       DO im = 1, nmolecules(is)
          this_im = locate(im,is)
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

!!$             IF (this_box == 2) THEN
!!$                write(203,*) atom_list(1,this_im,is)%rxp, atom_list(1,this_im,is)%ryp, &
!!$                     atom_list(1,this_im,is)%rzp
!!$             END IF
!!$             write(*,*) 'cubic'
             
          ELSE
             
             CALL Minimum_Image_Separation(this_box,xcom_old,ycom_old,zcom_old, &
                  xcom_new, ycom_new, zcom_new)

!             write(*,*) 'minimum'
             
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
          
          nmols(is,this_box) = nmols(is,this_box) + 1
          
          CALL Compute_Max_Com_Distance(this_im,is)
       END DO
    END DO
    
    
    IF ( f_read_old) THEN
       
       DO is = 1, nspecies
          species_list(is)%nmoltotal = SUM(nmols(is,:))
       END DO
       
       WRITE(logunit,*) 'Configurations read successfully'
    END IF


    DO ibox = 1, nbr_boxes
       IF(int_vdw_sum_style(ibox) == vdw_cut_tail) CALL Compute_Beads(ibox)
    END DO

    IF(ALLOCATED(total_molecules)) DEALLOCATE(total_molecules)

  END SUBROUTINE Read_Checkpoint

!*******************************************************************************************************
  
  SUBROUTINE Restart_From_Old
    !***************************************************************************************************
    ! The subroutine reads in a configuration to start a new simulation run. The format of the input
    ! file is identical to the checkpoint file in terms of atomic coordinates.
    !
    !
    !***************************************************************************************************
    
    IMPLICIT NONE
    
    INTEGER :: ibox, is, im, ia, nstart, nend, this_im, mols_this, nfrac_global, i, alive, j, this_rxnum
    INTEGER :: alive_new, counter, m, i_lambda
    
    INTEGER, DIMENSION(:), ALLOCATABLE ::  total_molecules_this
    INTEGER, DIMENSION(:), ALLOCATABLE :: n_int
    
    REAL(DP) :: this_lambda
    REAL(DP) :: E_recip, E_self, E_intra
    REAL(DP) :: E_old, xcom_old, ycom_old, zcom_old
    REAL(DP) :: xcom_new, ycom_new, zcom_new
    LOGICAL :: overlap
    
    Type(Energy_Class) :: inrg
  
    ! Loop over total number of boxes to read in the atomic coordinates

    ALLOCATE(total_molecules_this(nspecies))

    ALLOCATE(n_int(nspecies))

    nmols(:,:) = 0
    n_int(:) = 0
    
    DO ibox = 1, nbr_boxes
       
       OPEN(unit = old_config_unit,file=old_config_file(ibox))
       ! Read the number of molecules for each of the species
       
       READ(old_config_unit,*) (total_molecules_this(is), is = 1, nspecies)
       
       ! Read in the coordinates of the molecules
       
       DO is = 1, nspecies
          
          ! sum the total number of molecules of this species upto ibox - 1. This information
          ! will then be used to provide a locate number for molecules of species 'is', in 'ibox'

          IF (ibox /= 1) THEN
             mols_this = SUM(nmols(is,1:ibox-1))
          ELSE
             mols_this = 0

          END IF

          DO im = 1, total_molecules_this(is)
             ! provide a linked number to the molecule
             locate(im + mols_this,is) = im + mols_this
             this_im = locate(im + mols_this,is)
!             write(*,*) locate(this_im,is)
             molecule_list(this_im,is)%live = .TRUE.
             this_lambda = 1.0_DP
            ! By default all the molecules are normal 
             molecule_list(this_im,is)%molecule_type = int_normal
             
             DO ia = 1, natoms(is)
                
                READ(old_config_unit,*)nonbond_list(ia,is)%element, &
                     atom_list(ia,this_im,is)%rxp, &
                     atom_list(ia,this_im,is)%ryp, &
                     atom_list(ia,this_im,is)%rzp
                ! set the cfc_lambda and exist flags for this atom
                molecule_list(this_im,is)%cfc_lambda = this_lambda
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

             atom_list(1:natoms(is),this_im,is)%rxp = atom_list(1:natoms(is),this_im,is)%rxp + &
                  xcom_new - xcom_old
             atom_list(1:natoms(is),this_im,is)%ryp = atom_list(1:natoms(is),this_im,is)%ryp + &
                  ycom_new - ycom_old
             atom_list(1:natoms(is),this_im,is)%rzp = atom_list(1:natoms(is),this_im,is)%rzp + &
                  zcom_new - zcom_old

             nmols(is,ibox) = nmols(is,ibox) + 1
             
          END DO
          
       END DO

       CLOSE(unit = old_config_unit)
       
    END DO

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

    nfrac_global = 0
    DO is = 1, nspecies
       species_list(is)%nmoltotal = SUM(nmols(is,:))
    END DO

    WRITE(logunit,*) 'Configurations read successfully'

    DO ibox = 1, nbr_boxes
       IF (int_vdw_sum_style(ibox) == vdw_cut_tail) CALL Compute_Beads(ibox)
    END DO

    DEALLOCATE( total_molecules_this)
    
  END SUBROUTINE Restart_From_Old
!***************************************************************************************

SUBROUTINE Write_Trials_Success
  ! This subroutine writes number of trials and acceptance of these trials at the
  ! end of a simulation

  IMPLICIT NONE

  INTEGER :: ibox, is, ifrag, ireac

  WRITE(logunit,*)
  WRITE(logunit,*)

  DO ibox = 1, nbr_boxes

     WRITE(logunit,'(A28,2X,I2)') ' Writing information for box', ibox
     WRITE(logunit,*) ' *********************************************'
     WRITE(logunit,*)

     IF (nvolumes(ibox) /= 0 ) THEN
        WRITE(logunit,'(A20,2x,A10,2x,A10)') 'Move', 'Trials', 'Success'
        WRITE(logunit,11) 'Volume', nvolumes(ibox), nvol_success(ibox)
     END IF

     DO is = 1, nspecies
        
        WRITE(logunit,*) 
        WRITE(logunit,*) ' ******************************************'
        WRITE(logunit,*) ' Writing information for species', is
        WRITE(logunit,*)
        WRITE(logunit,'(A20,2x,A10,2x,A10,2x,A10)') 'Move', 'Trials', 'Success', '% Success'
        
        ! translation

        IF (ntrials(is,ibox)%displacement /= 0 ) THEN
        
           WRITE(logunit,11) 'Translate', ntrials(is,ibox)%displacement, &
                nsuccess(is,ibox)%displacement &
                ,100.0*dble(nsuccess(is,ibox)%displacement)/dble(ntrials(is,ibox)%displacement)
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
           
           WRITE(logunit,11) 'Insertion',  ntrials(is,ibox)%insertion, &
                nsuccess(is,ibox)%insertion, &
                100.0*dble(nsuccess(is,ibox)%insertion)/dble(ntrials(is,ibox)%insertion)
        END IF

        ! deletion

        IF (ntrials(is,ibox)%deletion /= 0 ) THEN
           
           WRITE(logunit,11) 'Deletion', ntrials(is,ibox)%deletion, &
                nsuccess(is,ibox)%deletion, &
                100.0*dble(nsuccess(is,ibox)%deletion)/dble(ntrials(is,ibox)%deletion)

        END IF

        ! atom displacement

        IF (ntrials(is,ibox)%disp_atom /= 0 ) THEN
           
           WRITE(logunit,11) 'Atom Displacement', ntrials(is,ibox)%disp_atom, &
                nsuccess(is,ibox)%disp_atom, &
                100.0*dble(nsuccess(is,ibox)%disp_atom)/dble(ntrials(is,ibox)%disp_atom)

        END IF
        WRITE(logunit,*) '**************************************'
        WRITE(logunit,*)
     END DO
        
  END DO

11 FORMAT(A20,2x,I10,2x,I10,2x,f10.2)
12 FORMAT(3(I10,2x),F10.2)

  IF( SUM(nfragments) .GT. 0 ) THEN

     WRITE(logunit,*) 'Writing information about fragments'
  
     DO is = 1, nspecies
        
        IF (species_list(is)%fragment) THEN
 
           WRITE(logunit,*)
           WRITE(logunit,*)'*************************************'
           WRITE(logunit,'(A31,2X,I2)') 'Writing information for species', is
           WRITE(logunit,'(A10,2x,A10,2x,A10,2X,A10)') 'Fragment', 'Trials', 'Success', '% Success'

           DO ifrag = 1, nfragments(is)
              
              IF (regrowth_trials(ifrag,is) /= 0 ) THEN
                 WRITE(logunit,12) ifrag,regrowth_trials(ifrag,is), &
                      regrowth_success(ifrag,is), &
                      100.0_DP * dble(regrowth_success(ifrag,is))/dble(regrowth_trials(ifrag,is))
              END IF
           END DO
           WRITE(logunit,*)'*********************************'
        END IF
     
     END DO

  END IF

END SUBROUTINE Write_Trials_Success

SUBROUTINE Write_Subroutine_Times

  IMPLICIT NONE

WRITE(logunit,*)
WRITE(logunit,*) 'Writing information about subroutine times'
WRITE(logunit,*)
WRITE(logunit,*) '******************************************'
WRITE(logunit,*)


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

END SUBROUTINE Write_Subroutine_Times

END MODULE Read_Write_Checkpoint

