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
SUBROUTINE Translate(this_box,mc_step)
!********************************************************************************

  !********************************************************************************
  ! This subroutine performs a COM translation move for a randomly chosen molecule.
  !
  ! 
  ! CALLED BY 
  !
  !        gcmc_driver
  !        gemc_driver
  !        nptmc_driver
  !        nvtmc_driver
  ! 
  ! CALLS
  !     
  !        Get_Index_Molecule
  !        Store_Molecule_Pair_Interaction_Arrays
  !        Compute_Molecule_Nonbond_Inter_Energy
  !        Save_Old_Cartesian_Coordinates
  !        Revert_Old_Cartestian_Coordinates
  !        Reset_Molecule_Pair_Interaction_Arrays
  !        Get_Position_Alive
  !        Compute_Ewald_Reciprocal_Energy_Difference
  !
  !
  !  12/10/13 : Beta release
  !********************************************************************************

  USE Type_Definitions
  USE Run_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE Pair_Nrg_Routines

  IMPLICIT NONE

  INTEGER  :: is, im, i, this_box, nmolecules_species, alive, total_mols, ibox, nmols_box
  INTEGER  :: dumcount, iatom, stupid_step
  INTEGER  :: mc_step, this_im, ia
  INTEGER  :: N_fracs, ind
  
  REAL(DP), ALLOCATABLE :: x_box(:), x_species(:)

  REAL(DP) :: dx, dy, dz, delta_e, delta_v, p_acc 
  REAL(DP) :: E_intra, E_vdw, E_qq
  REAL(DP) :: E_vdw_move, E_qq_move
  REAL(DP) :: E_reciprocal_move, rand_no
  REAL(DP) :: E_vdw_ch, E_qq_ch
  REAL(DP) :: old_value
  REAL(DP) :: W_vdw, W_qq, W_vdw_move, W_qq_move

  REAL(DP) :: success_ratio
  REAL(DP), ALLOCATABLE :: x_old(:), y_old(:), z_old(:)
  REAL(DP) :: rcut_small

  LOGICAL :: inter_overlap, overlap, update_flag
! Variables for prefferential sammpling

  INTEGER, ALLOCATABLE :: ni_before(:)
  REAL(DP) :: attempt_p 
  LOGICAL  :: inside_start

 ! Pair_Energy arrays and Ewald implementation

  INTEGER :: start, locate_im, count, this_species, position
!  REAL(DP), ALLOCATABLE :: pair_vdw_temp(:), pair_qq_temp(:)
  REAL(DP), ALLOCATABLE :: cos_mol_old(:), sin_mol_old(:)


! Done with that section

  E_vdw_move = 0.0_DP
  E_qq_move = 0.0_DP
  E_vdw = 0.0_DP
  E_qq = 0.0_DP
  E_reciprocal_move = 0.0_DP
  inter_overlap = .FALSE.
  attempt_p = 1.0_DP

  IF(nbr_boxes .GT. 1) THEN

    ! Choose a box based on its total mol fraction

    ALLOCATE(x_box(nbr_boxes)) 

    total_mols = SUM(nmols(:,:))
  
    ! If there are no molecules in the box then return to gcmc_driver

    IF (total_mols == 0) RETURN

    DO ibox = 1, nbr_boxes
       nmols_box = SUM(nmols(:,ibox))
       x_box(ibox) = REAL(nmols_box,DP)/REAL(total_mols,DP)
       IF ( ibox > 1 ) THEN
          x_box(ibox) = x_box(ibox) + x_box(ibox-1)
       END IF
    END DO
  
    ! Choose a box based on overall mol fraction

    rand_no = rranf()

    DO ibox = 1, nbr_boxes
       IF ( rand_no <= x_box(ibox)) EXIT
    END DO

    this_box = ibox
    DEALLOCATE(x_box)

  ELSE

    this_box = 1

  END IF

  ! Now choose species based on the mol fraction as well

  ALLOCATE(x_species(nspecies))

  nmolecules_species = SUM(nmols(:,this_box))
!     IF ( nmolecules_species /= 0 .AND. rranf() <=  prob_species_trans(is) ) EXIT
  IF( nmolecules_species == 0 ) RETURN

  DO is = 1, nspecies
     x_species(is) = REAL(nmols(is,this_box), DP) / REAL(nmolecules_species,DP)
     IF ( is > 1 ) THEN
        x_species(is) = x_species(is) + x_species(is-1)
     END IF
  END DO

  ! choose a species based on its mole fraction

  rand_no = rranf()
  
  DO is = 1, nspecies
     IF( rand_no <= x_species(is)) EXIT
  END DO

  IF (has_charge(is)) THEN
     rcut_small = MIN(rcut_vdw(this_box),rcut_coul(this_box))
  ELSE
     rcut_small = rcut_vdw(this_box)
  END IF
! number of molecules of species in this box

  nmolecules_species = nmols(is,this_box)

  DEALLOCATE(x_species)

  tot_trials(this_box) = tot_trials(this_box) + 1

  ! Select a molecule at random for displacement

  DO

     im = INT ( rranf() * nmolecules_species ) + 1
     ! Get the index of imth molecule of species is in this_box
     CALL Get_Index_Molecule(this_box,is,im,alive)

     EXIT
  END DO

  ! update the trial counter for this molecule
  
  ntrials(is,this_box)%displacement = ntrials(is,this_box)%displacement + 1

  ! Compute the LJ and real space intermolecular energies of this molecule

  IF (l_pair_nrg) THEN
     CALL Store_Molecule_Pair_Interaction_Arrays(alive,is,this_box,E_vdw,E_qq)
  ELSE
     CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_vdw,E_qq,inter_overlap)
  END IF

  CALL Save_Old_Cartesian_Coordinates(alive,is)

  IF (inter_overlap)  THEN
        WRITE(*,*) 'Disaster, overlap in the old configruation'
        WRITE(*,*) 'Translate'
        WRITE(*,*) alive, is, this_box
        WRITE(*,*) 'inter overlap', inter_overlap
!        STOP
  END IF

! Generate a random displacement vector. Note that the current formalism will
! work for cubic shaped boxes. However, it is easy to extend for nonorthorhombic
! boxes where displacements along the basis vectors. 

  IF ( box_list(this_box)%int_box_shape == int_cubic) THEN

     dx = ( 2.0_DP * rranf() - 1.0_DP) * max_disp(is,this_box)
     dy = ( 2.0_DP * rranf() - 1.0_DP) * max_disp(is,this_box)
     dz = ( 2.0_DP * rranf() - 1.0_DP) * max_disp(is,this_box)
          
  END IF

! Move atoms by the above vector dx,dy,dz and also update the COM

  atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp + dx
  atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp + dy
  atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp + dz

  molecule_list(alive,is)%xcom = molecule_list(alive,is)%xcom + dx
  molecule_list(alive,is)%ycom = molecule_list(alive,is)%ycom + dy
  molecule_list(alive,is)%zcom = molecule_list(alive,is)%zcom + dz


!**************************************************************************
     
  CALL Fold_Molecule(alive,is,this_box)
  
  CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_vdw_move,E_qq_move,inter_overlap)
  
   IF (inter_overlap) THEN ! Move is rejected
     
      CALL Revert_Old_Cartesian_Coordinates(alive,is)
      IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,this_box)

   ELSE

     delta_e = 0.0_DP

     IF ((int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is))) THEN

        ALLOCATE(cos_mol_old(nvecs(this_box)),sin_mol_old(nvecs(this_box)))
        CALL Get_Position_Alive(alive,is,position)
        
        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        cos_mol_old(:) = cos_mol(1:nvecs(this_box),position)
        sin_mol_old(:) = sin_mol(1:nvecs(this_box),position)
        !$OMP END PARALLEL WORKSHARE

        CALL Compute_Ewald_Reciprocal_Energy_Difference(alive,alive,is,this_box,int_translation,E_reciprocal_move)
        delta_e = E_reciprocal_move
        
     END IF
     
     ! Compute the difference in old and new energy
     
     delta_e = ( E_vdw_move - E_vdw ) + ( E_qq_move - E_qq ) + delta_e


     IF (int_sim_type == sim_nvt_min) THEN
        IF (delta_e  <= 0.0_DP) THEN
           p_acc = 1.0_DP
        ELSE
           p_acc = 0.0_DP
        END IF
     ELSE

         p_acc = MIN( 1.0_DP, attempt_p * DEXP(-beta(this_box) * delta_e) )

     END IF

     IF ( rranf() <= p_acc ) THEN

        ! accept the move and update the global energies
        energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + E_vdw_move - E_vdw
        energy(this_box)%inter_q   = energy(this_box)%inter_q   + E_qq_move - E_qq
        
        IF(int_charge_sum_style(this_box) == charge_ewald .AND. has_charge(is)) THEN
           energy(this_box)%ewald_reciprocal = energy(this_box)%ewald_reciprocal + E_reciprocal_move
        END IF
        
        energy(this_box)%total = energy(this_box)%total + delta_e

        ! update success counter
        
        nsuccess(is,this_box)%displacement = nsuccess(is,this_box)%displacement + 1
        nequil_success(is,this_box)%displacement = nequil_success(is,this_box)%displacement + 1

        IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
        IF (ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
        IF (ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)

     ELSE

        ! Revert to the old coordinates of atoms and com of the molecule
        CALL Revert_Old_Cartesian_Coordinates(alive,is) 

        IF ((int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is))) THEN
           ! Also reset the old cos_sum and sin_sum for reciprocal space vectors. Note
           ! that old vectors were set while difference in ewald reciprocal energy was computed.
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)           
           cos_sum(:,this_box) = cos_sum_old(:,this_box)
           sin_sum(:,this_box) = sin_sum_old(:,this_box)
           cos_mol(1:nvecs(this_box),position) =cos_mol_old(:)
           sin_mol(1:nvecs(this_box),position) =sin_mol_old(:)
           !$OMP END PARALLEL WORKSHARE
           DEALLOCATE(cos_mol_old,sin_mol_old)
        END IF
        IF (l_pair_nrg)  CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,this_box)
     ENDIF
     
  END IF

  IF ( MOD(ntrials(is,this_box)%displacement,nupdate) == 0 ) THEN
     IF ( int_run_style == run_equil ) THEN 
        success_ratio = REAL(nequil_success(is,this_box)%displacement,DP)/REAL(nupdate,DP)
        ELSE
            success_ratio = REAL(nsuccess(is,this_box)%displacement,DP)/REAL(ntrials(is,this_box)%displacement,DP)
     END IF

     WRITE(logunit,*)
     WRITE(logunit,'(A,2X,I3,2X,A,I2,2x,A,2X,F8.5)')'Translation ratio for species ', is , 'for box', this_box, 'is: ', success_ratio

     !nsuccess(is,this_box)%displacement = 0

     IF ( int_run_style == run_equil ) THEN

        ! check if the acceptace is close to 0.5

         nequil_success(is,this_box)%displacement = 0

         IF  ( success_ratio < 0.00005 ) THEN
             max_disp(is,this_box) = 0.1_DP*max_disp(is,this_box)
         ELSE
             max_disp(is,this_box) = MIN(rcut_small,2.0_DP*success_ratio*max_disp(is,this_box))
         END IF

         WRITE(logunit,'(A,1X,I1,1X,A,1X,I1)') 'Maximum width for translation of species', is,' updated in box', this_box
         WRITE(logunit,'(A,2X,F8.5)') 'new width is', max_disp(is,this_box)
         WRITE(logunit,*)
        
     END IF

  END IF

END SUBROUTINE Translate
     

     
