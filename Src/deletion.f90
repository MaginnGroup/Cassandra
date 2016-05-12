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

SUBROUTINE Deletion(this_box)
  
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
  USE Global_Variables
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

  ! Local declarations
  INTEGER :: i, i_type               ! atom indices
  INTEGER :: ifrag                   ! fragment indices
  INTEGER :: im                      ! molecule INDEX
  INTEGER :: lm                      ! molecule LOCATE
  INTEGER :: is, is_rand, is_counter ! species indices
  INTEGER :: kappa_tot, which_anchor
  INTEGER, ALLOCATABLE :: frag_order(:)
  INTEGER :: k, mcstep

  REAL(DP) :: delta_e, delta_e_pacc
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq, E_periodic_qq
  REAL(DP) :: E_reciprocal, E_self, E_lrc
  REAL(DP) :: nrg_ring_frag_tot
  REAL(DP) :: ln_pacc, P_seq, P_bias, this_lambda
  REAL(DP) :: E_intra_vdw_igas, E_intra_qq_igas

  LOGICAL :: inter_overlap, cbmc_overlap, intra_overlap
  LOGICAL :: accept_or_reject, fh_outside_bounds, accept

  INTEGER :: old_subensemble

  ! Initialize variables
  ln_pacc = 0.0_DP
  P_seq = 1.0_DP
  P_bias = 1.0_DP
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
  im = INT(rranf() * nmols(is,this_box)) + 1
  lm = locate(im,is,this_box)

  ! If fractional molecule, need to use lambda in reverse insertion move
  this_lambda = molecule_list(lm,is)%frac

  ! Save the coordinates of 'lm' because Build_Molecule will erase them if
  ! cbmc_overlap is tripped.

  CALL Save_Old_Cartesian_Coordinates(lm,is)  

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
     del_flag = .TRUE.      ! Don't change the coordinates of 'lm'
     get_fragorder = .TRUE. !
     ALLOCATE(frag_order(nfragments(is)))
     CALL Build_Molecule(lm,is,this_box,frag_order,this_lambda, &
             P_seq, P_bias, nrg_ring_frag_tot, cbmc_overlap)
     DEALLOCATE(frag_order)
     
     ! cbmc_overlap will only trip if the molecule being deleted had bad
     ! contacts
     IF (cbmc_overlap) THEN
        err_msg = ""
        err_msg(1) = "Energy overlap detected in existing configuration"
        err_msg(2) = "of molecule " // TRIM(Int_To_String(lm)) // " of species " // TRIM(Int_To_String(is))
        CALL Clean_Abort(err_msg, "Deletion")
     END IF

     ! So far P_bias only includes the probability of choosing the insertion 
     ! point from the collection of trial coordinates times the probability of 
     ! choosing each dihedral from the collection of trial dihedrals. We need 
     ! to include the number of trial coordinates, kappa_ins, and the number of
     ! of trial dihedrals, kappa_dih, for each dihedral.
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
  CALL Get_COM(lm,is)

  ! Compute the distance of the atom farthest from COM
  CALL Compute_Max_COM_Distance(lm,is)

  ! 4.1) Nonbonded intermolecular energies

  IF (l_pair_nrg) THEN
     CALL Store_Molecule_Pair_Interaction_Arrays(lm,is,this_box, &
             E_inter_vdw,E_inter_qq)
  ELSE
     CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is, &
             E_inter_vdw,E_inter_qq,inter_overlap)
  END IF

  delta_e = - E_inter_vdw - E_inter_qq

  ! 4.2) Bonded intramolecular energies

  CALL Compute_Molecule_Bond_Energy(lm,is,E_bond)
  CALL Compute_Molecule_Angle_Energy(lm,is,E_angle)
  CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihedral)
  CALL Compute_Molecule_Improper_Energy(lm,is,E_improper)

  delta_e = delta_e - E_bond - E_angle - E_dihedral - E_improper  
  
  ! 4.3) Nonbonded intramolecular energies

  CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is, &
          E_intra_vdw,E_intra_qq,E_periodic_qq,intra_overlap) 
  E_inter_qq = E_inter_qq + E_periodic_qq

  delta_e = delta_e - E_intra_vdw - E_intra_qq - E_periodic_qq

  ! 4.4) Ewald energies

  IF (int_charge_style(this_box) == charge_coul) THEN

      IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. &
           (has_charge(is)) ) THEN
    
         CALL Update_System_Ewald_Reciprocal_Energy(lm,is,this_box, &
                 int_deletion,E_reciprocal)

         delta_e = delta_e + (E_reciprocal - energy(this_box)%ewald_reciprocal)

      END IF

      CALL Compute_Molecule_Self_Energy(lm,is,this_box,E_self)

      delta_e = delta_e - E_self 
  END IF

  ! 4.5) Long-range energy correction

  IF (int_vdw_sum_style(this_box) == vdw_cut_tail .AND. &
			int_vdw_style(this_box) == vdw_lj ) THEN

     ! subtract off beads for this species
     nbeads_out(:) = nint_beads(:,this_box)

     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,this_box) = nint_beads(i_type,this_box) - 1
     END DO

     CALL Compute_LR_correction(this_box,e_lrc)
     delta_e = delta_e + ( e_lrc - energy(this_box)%lrc )

  ELSEIF (int_vdw_sum_style(this_box) == vdw_cut_tail .AND. &
			int_vdw_style(this_box) == vdw_mie ) THEN
     nbeads_out(:) = nint_beads_mie(is,:,this_box)
     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads_mie(is,i_type,this_box) = nint_beads_mie(is,i_type,this_box) - 1
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
     CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is, &
             E_intra_vdw_igas,E_intra_qq_igas,E_periodic_qq,intra_overlap)
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


     IF (int_charge_style(this_box) == charge_coul) THEN

         IF ( int_charge_sum_style(this_box) == charge_ewald .AND. &
              has_charge(is)) THEN
            energy(this_box)%ewald_reciprocal = E_reciprocal
         END IF

         energy(this_box)%self = energy(this_box)%self - E_self

     END IF


     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
        energy(this_box)%lrc = E_lrc
     END IF

     ! remove the deleted LOCATE from this_box's list
     IF (im < nmols(is,this_box)) THEN
        DO k = im + 1, nmols(is,this_box)
           locate(k-1,is,this_box) = locate(k,is,this_box)
        END DO
     END IF
     locate(nmols(is,this_box),is,this_box) = 0
     
     ! move LOCATE to list of unused LOCATES
     nmols(is,0) = nmols(is,0) + 1
     locate(nmols(is,0),is,0) = lm

     ! update the number of molecules
     molecule_list(lm,is)%live = .FALSE.
     atom_list(:,lm,is)%exist = .FALSE.
     nmols(is,this_box) = nmols(is,this_box) - 1

     ! Increment counter
     nsuccess(is,this_box)%deletion = nsuccess(is,this_box)%deletion + 1

!     CALL System_Energy_Check(1,mcstep,randno)
  ELSE

     IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. &
           (has_charge(is)) ) THEN
        ! Restore cos_sum and sin_sum. Note that these were changed when
        ! difference in reciprocal energies was computed
        cos_sum(:,this_box) = cos_sum_old(:,this_box)
        sin_sum(:,this_box) = sin_sum_old(:,this_box)
     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail ) THEN
        ! Restore the total number of bead types
        IF (int_vdw_style(this_box) == vdw_lj) THEN
	   nint_beads(:,this_box) = nbeads_out(:)
        ELSEIF (int_vdw_style(this_box) == vdw_mie) THEN
	   nint_beads_mie(is,:,this_box) = nbeads_out(:)
	ENDIF
     END IF

  END IF

  IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)


END SUBROUTINE Deletion

     

