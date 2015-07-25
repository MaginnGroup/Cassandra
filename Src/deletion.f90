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

SUBROUTINE Deletion(this_box,mcstep,randno)
  
  !*****************************************************************************
  !
  ! PURPOSE: attempt to delete a molecule that was inserted via 
  !          configurational bias monte carlo
  !
  ! Called by
  !
  !  gcmc_driver
  ! 
  ! Revision History
  !
  !   12/10/13 : Beta Release
  !   Version 1.1
  !     04/21/15  Corrected acceptance criteria
  !     05/05/15  Documented this code
  !
  ! DESCRIPTION: This subroutine performs the following steps:
  !
  ! Step 1) Select a species with uniform probability
  ! Step 2) Select a molecule with uniform probability
  ! Step 3) Calculate the bias probability for the reverse insertion move
  ! Step 4) Calculate the change in potential energy if the molecule is deleted
  ! Step 5) Accept or reject the move
  !
  !*****************************************************************************
  
  USE Type_Definitions
  USE Run_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE Fragment_Growth
  USE IO_Utilities

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(INOUT) :: this_box ! attempt to delete a molecule in this_box
  INTEGER :: mcstep  ! not used
  REAL(DP) :: randno ! not used

  ! Local declarations
  INTEGER :: i, i_type               ! atom indices
  INTEGER :: ifrag                   ! fragment indices
  INTEGER :: im, alive               ! molecule indices
  INTEGER :: is, is_rand, is_counter ! species indices
  INTEGER :: kappa_tot, which_anchor
  INTEGER, ALLOCATABLE :: frag_order(:)
  INTEGER :: k, position

  REAL(DP) :: delta_e, delta_e_pacc
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq
  REAL(DP) :: E_reciprocal_move, E_self_move, E_lrc
  REAL(DP) :: nrg_ring_frag_tot
  REAL(DP) :: ln_pacc, P_seq, P_bias, this_lambda
  REAL(DP) :: E_intra_vdw_igas, E_intra_qq_igas

  LOGICAL :: inter_overlap, cbmc_overlap, intra_overlap
  LOGICAL :: accept, accept_or_reject

  ! Initialize variables
  ln_pacc = 0.0_DP
  P_seq = 1.0_DP
  P_bias = 1.0_DP
  this_lambda = 1.0_DP
  nrg_ring_frag_tot = 0.0_DP
  inter_overlap = .FALSE.
  cbmc_overlap = .FALSE.
  intra_overlap = .FALSE.
  
  !*****************************************************************************
  ! Step 1) Select a species with uniform probability
  !*****************************************************************************
  !
  ! All species may not be insertable. For example, in a simulation of dilute
  ! water (species 3) and CO2 (species 4) in an ionic liquid (species 1 and 2), 
  ! the number of ionic liquid molecules may be fixed and only the numbers of
  ! water and CO2 allowed to fluctuate. First, choose a random integer between 1
  ! and the number of insertable species, nspec_insert:

  is_rand = INT(rranf() * nspec_insert) + 1

  ! Now find the index 'is' that corresponds to is_rand. In the example, if
  ! is_rand == 2 a CO2 molecule will be deleted. CO2 corresponds to 'is' == 4.

  is_counter = 0
  DO is = 1, nspecies
     IF(species_list(is)%int_species_type == int_sorbate) THEN
        is_counter = is_counter + 1
     END IF
     IF(is_counter == is_rand) EXIT ! exit the loop when 'is' has been found
  END DO

  ! In the given example, now 'is' would equal 4.

  ! Cannot delete a molecule if there aren't any in the box
  IF (nmols(is,this_box) == 0) RETURN

  ! Now that a deletion will be attempted, we need to do some bookkeeping:
  !  * Increment the counters to compute success ratios

  ntrials(is,this_box)%deletion = ntrials(is,this_box)%deletion + 1  
  tot_trials(this_box) = tot_trials(this_box) + 1

  !*****************************************************************************
  ! Step 2) Select a molecule with uniform probability
  !*****************************************************************************
  !

  im = INT(rranf() * nmols(is,this_box)) + 1
  CALL Get_Index_Molecule(this_box,is,im,alive) ! sets the value of 'alive'

  ! Save the coordinates of 'alive' because Build_Molecule will erase them if
  ! cbmc_overlap is tripped.

  CALL Save_Old_Cartesian_Coordinates(alive,is)  

  ! Compute the energy of the molecule
  
  !*****************************************************************************
  ! Step 3) Calculate the bias probability for the reverse insertion move
  !*****************************************************************************
  !
  ! The bias probability, P_bias, of the reverse insertion move is required to
  ! calculate the probability of accepting the deletion. P_bias will be 
  ! calculated using the following procedure:
  ! 
  !   3.1) Select the order to insert fragments, with probability P_seq
  !   3.2) Select kappa_ins - 1 trial coordinates, each with uniform probability
  !   3.3) Calculate the probability of the fragment's current COM
  !   3.4) For each additional fragment:
  !          a) Select kappa_dih - 1 trial dihedrals, each with uniform 
  !             probability
  !          b) Calculate the probability of the fragment's current dihedral
  !
  ! These steps are implemented in the subroutine Build_Molecule
  
  IF(species_list(is)%fragment .AND. &
     (species_list(is)%int_insert .NE. int_igas)) THEN

     ! Build_Molecule places the first fragment, then calls Fragment_Placement 
     ! to place the additional fragments
     del_flag = .TRUE.      ! Don't change the coordinates of 'alive'
     get_fragorder = .TRUE. !
     ALLOCATE(frag_order(nfragments(is)))
     CALL Build_Molecule(alive,is,this_box,frag_order,this_lambda, &
             P_seq, P_bias, nrg_ring_frag_tot, cbmc_overlap)
     DEALLOCATE(frag_order)
     
     ! cbmc_overlap will only trip if the molecule being deleted had bad
     ! contacts
     IF (cbmc_overlap) THEN
        WRITE(*,*)
        WRITE(*,*) 'Warning....energy overlap detected in old configuration in deletion.f90'
        WRITE(*,*) 'molecule, species', alive, is
        WRITE(*,*)

        ! Revert to the COM and Eulerian angles for the molecule
        CALL Revert_Old_Cartesian_Coordinates(alive,is)
        atom_list(1:natoms(is),alive,is)%exist = .TRUE.
        molecule_list(alive,is)%cfc_lambda = this_lambda
        
     END IF

     ! So far P_bias only includes the probability of choosing the 
     ! insertion point from the collection of trial coordinates times the 
     ! probability of choosing each dihedral from the collection of trial 
     ! dihedrals. We need to include the number of trial coordinates, kappa_ins,
     ! and the number of trial dihedrals, kappa_dih, for each dihedral.
     kappa_tot = 1

     IF (nfragments(is) /=0 ) THEN

        kappa_tot = kappa_tot * kappa_ins

        IF (kappa_rot /= 0) THEN
           kappa_tot = kappa_tot * kappa_rot
        END IF
        
        IF (kappa_dih /=0 ) THEN
           DO ifrag = 2, nfragments(is)
              kappa_tot = kappa_tot * kappa_dih
           END DO
        END IF
     
     END IF
     
     P_bias = P_bias * REAL(kappa_tot , DP)

  END IF

  !*****************************************************************************
  ! Step 4) Calculate the change in potential energy if the molecule is deleted
  !*****************************************************************************
  !
  ! Whether the deletion will be accepted depends on the change in potential
  ! energy, delta_e. The potential energy will be computed in 5 stages:
  !   4.1) Nonbonded intermolecular energies
  !   4.2) Bonded intramolecular energies
  !   4.3) Nonbonded intramolecular energies
  !   4.4) Ewald energies
  !   4.5) Long-range energy correction
  ! 

  ! Recompute the COM
  CALL Get_COM(alive,is)

  ! Compute the distance of the atom farthest from COM
  CALL Compute_Max_COM_Distance(alive,is)

  ! 4.1) Nonbonded intermolecular energies

  IF (l_pair_nrg) THEN
     CALL Store_Molecule_Pair_Interaction_Arrays(alive,is,this_box, &
             E_inter_vdw,E_inter_qq)
  ELSE
     CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is, &
             E_inter_vdw,E_inter_qq,inter_overlap)
  END IF

  delta_e = - E_inter_vdw - E_inter_qq

  ! 4.2) Bonded intramolecular energies

  CALL Compute_Molecule_Bond_Energy(alive,is,E_bond)
  CALL Compute_Molecule_Angle_Energy(alive,is,E_angle)
  CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihedral)
  CALL Compute_Molecule_Improper_Energy(alive,is,E_improper)

  delta_e = delta_e - E_bond - E_angle - E_dihedral - E_improper  
  
  ! 4.3) Nonbonded intramolecular energies

  CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is, &
          E_intra_vdw,E_intra_qq,intra_overlap)

  delta_e = delta_e - E_intra_vdw - E_intra_qq

  ! 4.4) Ewald energies

  IF ( (int_charge_sum_style(this_box) == charge_ewald .or. int_charge_sum_style(this_box) == charge_gaussian) .AND. &
       (has_charge(is)) ) THEN

     CALL Compute_Ewald_Reciprocal_Energy_Difference(alive,alive,is,this_box, &
             int_deletion,E_reciprocal_move)
     CALL Compute_Ewald_Self_Energy_Difference(alive,is,this_box, &
             int_deletion,E_self_move)

     delta_e = delta_e + E_self_move &
                       + (E_reciprocal_move - energy(this_box)%ewald_reciprocal)
  END IF

  ! 4.5) Long-range energy correction

  IF (int_vdw_sum_style(this_box) == vdw_cut_tail .or. int_vdw_sum_style(this_box) == born_cut_tail) THEN

     ! subtract off beads for this species
     nbeads_out(:) = nint_beads(:,this_box)

     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,this_box) = nint_beads(i_type,this_box) - 1
     END DO

     CALL Compute_LR_correction(this_box,e_lrc)
     delta_e = delta_e + ( e_lrc - energy(this_box)%lrc )

  END IF  
  
  !*****************************************************************************
  ! Step 5) Accept or reject the move
  !*****************************************************************************
  !
  ! The following quantity is calculated
  !
  !                  p_m a_mn 
  !    ln_pacc = Log[--------]
  !                  p_n a_nm
  !
  ! and passed to accept_or_reject() which executes the metropolis criterion.
  ! The acceptance criterion to delete a molecule that was inserted via CBMC is
  !
  !                                                          V
  !    ln_pacc = b(dU_mn + U_frag) + b mu' + Log[-----------------------]
  !                                              P_seq P_bias N Lambda^3
  !
  !                                          b f' V
  !            = b(dU_mn + U_frag) + Log[--------------]
  !                                      P_seq P_bias N 
  !
  ! where the primes (') indicate that additional intensive terms have been
  ! absorbed into the chemical potential and fugacity, respectively.

  IF(species_list(is)%int_insert == int_igas) THEN
     igas_flag = .TRUE.
     CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is, &
             E_intra_vdw_igas,E_intra_qq_igas,intra_overlap)
     igas_flag = .FALSE. 
     ln_pacc = beta(this_box) * (delta_e + E_bond + E_angle &
                                         + E_dihedral + E_improper &
                                         + E_intra_vdw_igas + E_intra_qq_igas)
  ELSE

     IF(species_list(is)%fragment) THEN
        ln_pacc = beta(this_box) * (delta_e + E_angle + nrg_ring_frag_tot)
     END IF

  END IF

  ! P_seq and P_bias equal 1.0 unless changed by Build_Molecule
  ln_pacc = ln_pacc - DLOG(P_seq * P_bias) &
                    - DLOG(REAL(nmols(is,this_box),DP)) &
                    + DLOG(box_list(this_box)%volume)
 
  IF(lchempot) THEN
     ! chemical potential is input
     ln_pacc = ln_pacc + beta(this_box) * species_list(is)%chem_potential &
                       - 3.0_DP*DLOG(species_list(is)%de_broglie(this_box)) 
  ELSE
     ! fugacity is input
     ln_pacc = ln_pacc + DLOG(species_list(is)%fugacity) &
                       + DLOG(beta(this_box))
  END IF 
 
  accept = accept_or_reject(ln_pacc)

  IF (accept) THEN
     ! Update energies
     energy(this_box)%total = energy(this_box)%total + delta_e
     energy(this_box)%intra = energy(this_box)%intra - E_bond - E_angle &
                            - E_dihedral - E_improper
     energy(this_box)%bond = energy(this_box)%bond - E_bond
     energy(this_box)%angle = energy(this_box)%angle - E_angle
     energy(this_box)%dihedral = energy(this_box)%dihedral - E_dihedral
     energy(this_box)%improper = energy(this_box)%improper - E_improper
     energy(this_box)%intra_vdw = energy(this_box)%intra_vdw - E_intra_vdw
     energy(this_box)%intra_q = energy(this_box)%intra_q - E_intra_qq
     energy(this_box)%inter_vdw = energy(this_box)%inter_vdw - E_inter_vdw
     energy(this_box)%inter_q   = energy(this_box)%inter_q - E_inter_qq

     IF ( (int_charge_sum_style(this_box) == charge_ewald .or. int_charge_sum_style(this_box) == charge_gaussian) .AND. &
          has_charge(is)) THEN
        energy(this_box)%ewald_reciprocal = E_reciprocal_move
        energy(this_box)%ewald_self = energy(this_box)%ewald_self + E_self_move
     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail .or. int_vdw_sum_style(this_box) == born_cut_tail) THEN
        energy(this_box)%lrc = E_lrc
     END IF

     ! obtain the original position of the deleted molecule so that the
     ! linked list can be updated

     CALL Get_Position_Molecule(this_box,is,im,position)

     IF (position < SUM(nmols(is,:))) THEN
        DO k = position + 1, SUM(nmols(is,:))
           locate(k-1,is) = locate(k,is)
        END DO
     END IF
     
     ! move the deleted molecule to the end of alive molecules
     locate(nmols(is,this_box),is) = alive
     molecule_list(alive,is)%live = .FALSE.
     atom_list(:,alive,is)%exist = .FALSE.
     
     ! update the number of molecules
     nmols(is,this_box) = nmols(is,this_box) - 1

     ! Increment counter
     nsuccess(is,this_box)%deletion = nsuccess(is,this_box)%deletion + 1

!     CALL System_Energy_Check(1,mcstep,randno)
  ELSE

     IF ( (int_charge_sum_style(this_box) == charge_ewald .or. int_charge_sum_style(this_box) == charge_gaussian) .AND. &
           (has_charge(is)) ) THEN
        ! Restore cos_sum and sin_sum. Note that these were changed when
        ! difference in reciprocal energies was computed
        cos_sum(:,this_box) = cos_sum_old(:,this_box)
        sin_sum(:,this_box) = sin_sum_old(:,this_box)
     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail .or. int_vdw_sum_style(this_box) == born_cut_tail ) THEN
        ! Restore the total number of bead types
        nint_beads(:,this_box) = nbeads_out(:)
     END IF

  END IF

  IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)

END SUBROUTINE Deletion

     

