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

SUBROUTINE Deletion
  
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

  ! Local declarations
  INTEGER :: i, i_type               ! atom indices
  INTEGER :: ifrag                   ! fragment indices
  INTEGER :: im                      ! molecule INDEX
  INTEGER :: lm                      ! molecule LOCATE
  INTEGER :: is, is_rand, is_counter ! species indices
  INTEGER :: ibox ! attempt to delete a molecule in ibox
  INTEGER :: which_anchor
  INTEGER, ALLOCATABLE :: frag_order(:)
  INTEGER :: k, mcstep

  REAL(DP) :: dE, dE_frag
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq, E_periodic_qq
  REAL(DP) :: E_reciprocal, E_self, E_lrc
  REAL(DP) :: nrg_ring_frag_tot
  REAL(DP) :: ln_pacc, ln_pseq, ln_pbias, this_lambda
  REAL(DP) :: E_intra_vdw_igas, E_intra_qq_igas

  LOGICAL :: inter_overlap, cbmc_overlap, intra_overlap
  LOGICAL :: accept_or_reject, fh_outside_bounds

  INTEGER :: old_subensemble

  ! Initialize variables
  ln_pacc = 0.0_DP
  ln_pseq = 0.0_DP
  ln_pbias = 0.0_DP
  nrg_ring_frag_tot = 0.0_DP
  inter_overlap = .FALSE.
  cbmc_overlap = .FALSE.
  intra_overlap = .FALSE.
  ibox = 1
  accept = .FALSE.
  
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
  IF (nmols(is,ibox) == 0) THEN
     WRITE(logunit,'(X,I9,X,A10,X,5X,X,I3,X,I3,X,L8,X,9X,X,A9)') &
           i_mcstep, 'delete' , is, ibox, .FALSE., 'no mols'
     RETURN
  END IF
 

  ! Now that a deletion will be attempted, we need to do some bookkeeping:
  !  * Increment the counters to compute success ratios
  ntrials(is,ibox)%deletion = ntrials(is,ibox)%deletion + 1  
  tot_trials(ibox) = tot_trials(ibox) + 1

  !*****************************************************************************
  ! Step 2) Select a molecule with uniform probability
  !*****************************************************************************
  im = INT(rranf() * nmols(is,ibox)) + 1
  lm = locate(im,is,ibox)

  ! If fractional molecule, need to use lambda in reverse insertion move
  this_lambda = molecule_list(lm,is)%frac

  ! Save the coordinates of 'lm' because Build_Molecule will erase them if
  ! cbmc_overlap is tripped.

  CALL Save_Old_Cartesian_Coordinates(lm,is)  

  !*****************************************************************************
  ! Step 3) Calculate the bias probability for the reverse insertion move
  !*****************************************************************************
  !
  ! The bias probability, ln_pbias, of the reverse insertion move is required to
  ! calculate the probability of accepting the deletion. ln_pbias will be 
  ! calculated using the following procedure:
  ! 
  !   3.1) Select the order to insert fragments, with probability ln_pseq
  !   3.2) Select kappa_ins - 1 trial coordinates, each with uniform probability
  !   3.3) Calculate the probability of the fragment's current COM
  !   3.4) For each additional fragment:
  !          a) Select kappa_dih - 1 trial dihedrals, each with uniform 
  !             probability
  !          b) Calculate the probability of the fragment's current dihedral
  !
  ! These steps are implemented in the subroutine Build_Molecule
  

  ! Build_Molecule places the first fragment, then calls Fragment_Placement 
  ! to place the additional fragments
  del_flag = .TRUE.      ! Don't change the coordinates of 'lm'
  get_fragorder = .TRUE. !
  ALLOCATE(frag_order(nfragments(is)))
  CALL Build_Molecule(lm,is,ibox,frag_order,this_lambda, &
          ln_pseq, ln_pbias, nrg_ring_frag_tot, cbmc_overlap)
  DEALLOCATE(frag_order)
  
  ! cbmc_overlap will only trip if the molecule being deleted had bad contacts
  IF (cbmc_overlap) THEN
     err_msg = ""
     err_msg(1) = "Attempted to delete molecule " // TRIM(Int_To_String(im)) // &
                  " of species " // TRIM(Int_To_String(is))
     err_msg(2) = "but the molecule energy is too high"
     IF (start_type(ibox) == "make_config" ) THEN
        err_msg(3) = "Try increasing Rcutoff_Low, increasing the box size, or "
        err_msg(4) = "decreasing the initial number of molecules"
     END IF
     CALL Clean_Abort(err_msg, "Deletion")
  END IF

  ! So far ln_pbias includes 
  !   * the probability of choosing the insertion 
  !     point from the collection of trial coordinates 
  !   * the probability of choosing each dihedral from the collection of trial dihedrals. 
  ! Now add
  !   * the probability of the fragment sequence, ln_pseq
  !   * the number of trial coordinates, kappa_ins
  !   * the number of trial dihedrals, kappa_dih, for each dihedral.

  ln_pbias = ln_pbias + ln_pseq
  ln_pbias = ln_pbias + DLOG(REAL(kappa_ins,DP))

  IF (kappa_rot /= 0 ) THEN
     ln_pbias = ln_pbias + DLOG(REAL(kappa_rot,DP))
  END IF

  IF (kappa_dih /= 0 ) THEN
     ln_pbias = ln_pbias + REAL(nfragments(is)-1,DP) * DLOG(REAL(kappa_dih,DP))
  END IF
  
  !*****************************************************************************
  ! Step 4) Calculate the change in potential energy if the molecule is deleted
  !*****************************************************************************
  !
  ! Whether the deletion will be accepted depends on the change in potential
  ! energy, dE. The potential energy will be computed in 5 stages:
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
     CALL Store_Molecule_Pair_Interaction_Arrays(lm,is,ibox, &
             E_inter_vdw,E_inter_qq)
  ELSE
     CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is, &
             E_inter_vdw,E_inter_qq,inter_overlap)
  END IF

  dE = - E_inter_vdw - E_inter_qq

  ! 4.2) Bonded intramolecular energies

  CALL Compute_Molecule_Bond_Energy(lm,is,E_bond)
  CALL Compute_Molecule_Angle_Energy(lm,is,E_angle)
  CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihedral)
  CALL Compute_Molecule_Improper_Energy(lm,is,E_improper)

  dE = dE - E_bond - E_angle - E_dihedral - E_improper  
  
  ! 4.3) Nonbonded intramolecular energies

  CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is, &
          E_intra_vdw,E_intra_qq,E_periodic_qq,intra_overlap) 
  E_inter_qq = E_inter_qq + E_periodic_qq

  dE = dE - E_intra_vdw - E_intra_qq - E_periodic_qq

  ! 4.4) Ewald energies

  IF (int_charge_style(ibox) == charge_coul) THEN

      IF ( (int_charge_sum_style(ibox) == charge_ewald) .AND. &
           (has_charge(is)) ) THEN
    
         CALL Update_System_Ewald_Reciprocal_Energy(lm,is,ibox, &
                 int_deletion,E_reciprocal)

         dE = dE + (E_reciprocal - energy(ibox)%ewald_reciprocal)

      END IF

      CALL Compute_Molecule_Self_Energy(lm,is,ibox,E_self)

      dE = dE - E_self 
  END IF

  ! 4.5) Long-range energy correction

  IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN

     ! subtract off beads for this species
     nbeads_out(:) = nint_beads(:,ibox)

     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,ibox) = nint_beads(i_type,ibox) - 1
     END DO

     CALL Compute_LR_correction(ibox,e_lrc)
     dE = dE + ( e_lrc - energy(ibox)%lrc )

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
  !                                                               V
  !    ln_pacc = b(dU_mn + U_frag) + b mu' + - ln_pbias + Log[----------]
  !                                                           N Lambda^3
  !
  ! where the primes (') indicate that additional intensive terms have been
  ! absorbed into the chemical potential.

  ! change in energy less energy used to bias fragment selection
  dE_frag = - E_angle - nrg_ring_frag_tot
  ln_pacc = beta(ibox) * (dE - dE_frag)

  ! chemical potential
  ln_pacc = ln_pacc + beta(ibox) * species_list(is)%chem_potential

  ! CBMC bias probability
  ln_pacc = ln_pacc - ln_pbias

  ! dimensionless density
  ln_pacc = ln_pacc + DLOG(box_list(ibox)%volume) &
                    - DLOG(REAL(nmols(is,ibox),DP)) &
                    - 3.0_DP*DLOG(species_list(is)%de_broglie(ibox)) 
 
  accept = accept_or_reject(ln_pacc)

  IF (accept) THEN
     ! Update energies
     energy(ibox)%total = energy(ibox)%total + dE
     energy(ibox)%intra = energy(ibox)%intra - E_bond - E_angle &
                            - E_dihedral - E_improper
     energy(ibox)%bond = energy(ibox)%bond - E_bond
     energy(ibox)%angle = energy(ibox)%angle - E_angle
     energy(ibox)%dihedral = energy(ibox)%dihedral - E_dihedral
     energy(ibox)%improper = energy(ibox)%improper - E_improper
     energy(ibox)%intra_vdw = energy(ibox)%intra_vdw - E_intra_vdw
     energy(ibox)%intra_q = energy(ibox)%intra_q - E_intra_qq
     energy(ibox)%inter_vdw = energy(ibox)%inter_vdw - E_inter_vdw
     energy(ibox)%inter_q   = energy(ibox)%inter_q - E_inter_qq


     IF (int_charge_style(ibox) == charge_coul) THEN

         IF ( int_charge_sum_style(ibox) == charge_ewald .AND. &
              has_charge(is)) THEN
            energy(ibox)%ewald_reciprocal = E_reciprocal
         END IF

         energy(ibox)%self = energy(ibox)%self - E_self

     END IF


     IF ( int_vdw_sum_style(ibox) == vdw_cut_tail) THEN
        energy(ibox)%lrc = E_lrc
     END IF

     ! remove the deleted LOCATE from ibox's list
     IF (im < nmols(is,ibox)) THEN
        DO k = im + 1, nmols(is,ibox)
           locate(k-1,is,ibox) = locate(k,is,ibox)
        END DO
     END IF
     locate(nmols(is,ibox),is,ibox) = 0
     
     ! move LOCATE to list of unused LOCATES
     nmols(is,0) = nmols(is,0) + 1
     locate(nmols(is,0),is,0) = lm

     ! update the number of molecules
     molecule_list(lm,is)%live = .FALSE.
     atom_list(:,lm,is)%exist = .FALSE.
     nmols(is,ibox) = nmols(is,ibox) - 1

     ! Increment counter
     nsuccess(is,ibox)%deletion = nsuccess(is,ibox)%deletion + 1

!     CALL Check_System_Energy(1,randno)
  ELSE

     IF ( (int_charge_sum_style(ibox) == charge_ewald) .AND. &
           (has_charge(is)) ) THEN
        ! Restore cos_sum and sin_sum. Note that these were changed when
        ! difference in reciprocal energies was computed
        cos_sum(:,ibox) = cos_sum_old(:,ibox)
        sin_sum(:,ibox) = sin_sum_old(:,ibox)
     END IF

     IF ( int_vdw_sum_style(ibox) == vdw_cut_tail ) THEN
        ! Restore the total number of bead types
        nint_beads(:,ibox) = nbeads_out(:)
     END IF

  END IF

  IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)

  IF (verbose_log) THEN
     WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,F9.3)') &
           i_mcstep, 'delete' , lm, is, ibox, accept, ln_pacc
  END IF


END SUBROUTINE Deletion

     

