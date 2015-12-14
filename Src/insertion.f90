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

SUBROUTINE Insertion(this_box,accept)

  !*****************************************************************************
  ! 
  ! PURPOSE: attempt to insert a molecule via configurational bias monte carlo
  !
  ! Called by
  !
  !    gcmc_driver
  !
  ! Revision history
  !
  !   12/10/13  : Beta version created
  !   Version 1.1
  !     04/21/15  Corrected acceptance criteria
  !     05/01/15  Documented this code
  !   
  ! DESCRIPTION: This subroutine performs the following steps:
  !
  ! Step 1) Select a species with uniform probability
  ! Step 2) Choose a position, orientation and conformation for the 
  !         to-be-inserted molecule
  ! Step 3) Calculate the change in potential energy if the molecule is inserted
  ! Step 4) Accept or reject the move
  !
  !*****************************************************************************

  USE Global_Variables
  USE Energy_Routines
  USE IO_Utilities
  USE Random_Generators
  USE Rotation_Routines
  USE Fragment_Growth

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  IMPLICIT NONE

  ! Arguments
  INTEGER :: this_box ! attempt to insert a molecule in this_box
  LOGICAL :: accept

  ! Local declarations
  INTEGER :: i, i_type               ! atom indices
  INTEGER :: ifrag                   ! fragment indices
  INTEGER :: im                      ! molecule INDEX
  INTEGER :: lm                      ! molecule LOCATE
  INTEGER :: is, is_rand, is_counter ! species indices
  INTEGER :: kappa_tot, which_anchor
  INTEGER, ALLOCATABLE :: frag_order(:)
  INTEGER :: rand_igas, tot_mols

  REAL(DP) :: dx, dy, dz, delta_e
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq, E_periodic_qq
  REAL(DP) :: E_reciprocal, E_self, E_lrc
  REAL(DP) :: nrg_ring_frag_tot
  REAL(DP) :: ln_pacc, P_seq, P_bias, this_lambda

  LOGICAL :: inter_overlap, cbmc_overlap, intra_overlap
  LOGICAL :: accept_or_reject

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
  ! Step 1) Randomly select a species
  !*****************************************************************************
  !
  ! All species may not be insertable. For example, in a simulation of dilute
  ! water (species 3) and CO2 (species 4) in an ionic liquid (species 1 and 2), 
  ! the number of ionic liquid molecules may be fixed and only the numbers of
  ! water and CO2 allowed to fluctuate. First, choose a random integer between 1
  ! and the number of insertable species, nspec_insert:
  is_rand = INT(rranf() * nspec_insert) + 1

  ! Now find the index 'is' that corresponds to is_rand. In the example, if
  ! is_rand == 2 a CO2 molecule will be inserted. CO2 corresponds to 'is' == 4.
  is_counter = 0
  DO is = 1, nspecies
     IF(species_list(is)%int_species_type == int_sorbate) THEN 
        is_counter = is_counter + 1
     END IF
     IF(is_counter == is_rand) EXIT ! exit the loop when 'is' has been found
  END DO
  ! In the given example, now 'is' would equal 4.

  ! Each species has a maximum allowable number of molecules specified in the 
  ! input file. The number of molecules currently in the system is
  tot_mols = SUM(nmols(is,1:nbr_boxes)) ! summed over the number of boxes

  ! Check that tot_mols is less than the maximum allowable, max_molecules(is)
  IF (tot_mols == max_molecules(is)) THEN
     err_msg = ""
     err_msg(1) = 'Number of molecule exceeds limit of ' // &
                  Int_To_String(tot_mols) 
     err_msg(2) = 'Increase molecule number limit in input file '
     CALL Clean_Abort(err_msg,'Insertion')
     ! exit if we are attempting an insertion above the maximum allowable
  END IF

  ! Now that an insertion will be attempted, we need to do some bookkeeping:

  !  * Increment the counters to compute success ratios
  ntrials(is,this_box)%insertion = ntrials(is,this_box)%insertion + 1
  tot_trials(this_box) = tot_trials(this_box) + 1

  !  * Assign a LOCATE for this molecule from the list of unused LOCATEs
  nmols(is,this_box) = nmols(is,this_box) + 1
  locate(nmols(is,this_box),is,this_box) = locate(nmols(is,0),is,0)
  locate(nmols(is,0),is,0) = 0
  nmols(is,0) = nmols(is,0) - 1

  !  * Set properties of the to-be-inserted molecule
  lm = locate(nmols(is,this_box),is,this_box)
  molecule_list(lm,is)%which_box = this_box
  molecule_list(lm,is)%frac = this_lambda
  molecule_list(lm,is)%molecule_type = int_normal

  ! With the bookkeeping completed, we are ready to attempt the insertion
  
  !*****************************************************************************
  ! Step 2) Choose a position, orientation and conformation for the 
  !         to-be-inserted molecule
  !*****************************************************************************
  !
  ! If the molecule:
  !   * has a fragment library and is not an ideal gas, then the conformation
  !     will be grown fragment-by-fragment using CBMC. The resulting 
  !     conformation will be molded to a high probability position.
  !   * has no fragment library, then there is only one conformation. Position 
  !     and orientation are random.
  !   * is an ideal gas, then molecular conformations (not fragment 
  !     conformations?) were sampled according to their Boltzmann weight. One 
  !     is chosen at random. Position and orientation are random.
  
  IF (species_list(is)%fragment .AND. &
     (species_list(is)%int_insert .NE. int_igas) ) THEN

     ! Build_Molecule places the first fragment, then calls Fragment_Placement
     ! to place the additional fragments 
     del_flag = .FALSE.     ! Change the coordinates of 'lm'
     get_fragorder = .TRUE.
     ALLOCATE(frag_order(nfragments(is)))
     CALL Build_Molecule(lm,is,this_box,frag_order,this_lambda, &
             P_seq,P_bias,nrg_ring_frag_tot,cbmc_overlap)
     DEALLOCATE(frag_order)

     ! Turn the molecule on
     molecule_list(lm,is)%live = .TRUE.
     atom_list(:,lm,is)%exist = .TRUE.

     ! So far P_bias only includes the probability of choosing the insertion 
     ! point from the collection of trial coordinates times the probability of 
     ! choosing each dihedral from the collection of trial dihedrals. We need 
     ! to include the number of trial coordinates, kappa_ins, and the number of
     ! of trial dihedrals, kappa_dih, for each dihedral.
     kappa_tot = 1

     IF (nfragments(is) /= 0 ) THEN

        kappa_tot = kappa_tot * kappa_ins

        IF (kappa_rot /= 0 ) THEN
           kappa_tot = kappa_tot * kappa_rot
        END IF

        IF (kappa_dih /= 0 ) THEN
           DO ifrag = 1, nfragments(is) - 1
              kappa_tot = kappa_tot * kappa_dih
           END DO
        END IF

     END IF

     P_bias = P_bias * REAL(kappa_tot, DP)

  ELSE
 
     ! The molecule does not have a fragment libary or is an ideal gas. The 
     ! position and orientation will be randomized independent of the molecular 
     ! conformation.

     ! Turn the molecule on
     atom_list(:,lm,is)%exist = .TRUE.
     molecule_list(lm,is)%live = .TRUE.

     ! Now we need to grab the molecule's conformation

     IF(species_list(is)%int_insert == int_random) THEN
     
        ! A rigid molecule has only one possible conformation

        molecule_list(lm,is)%xcom = species_list(is)%xcom
        molecule_list(lm,is)%ycom = species_list(is)%ycom
        molecule_list(lm,is)%zcom = species_list(is)%zcom
        
        atom_list(:,lm,is)%rxp = init_list(:,1,is)%rxp
        atom_list(:,lm,is)%ryp = init_list(:,1,is)%ryp
        atom_list(:,lm,is)%rzp = init_list(:,1,is)%rzp


     ELSE IF(species_list(is)%int_insert == int_igas) THEN

        ! An ideal gas molecule can adopt any conformation independent of its 
        ! local environment. The conformations are sampled according to their 
        ! Boltzmann weight. Read one in at random:

        rand_igas = (rranf() * n_igas(is)) + 1

        molecule_list(lm,is)%xcom = molecule_list_igas(rand_igas,is)%xcom
        molecule_list(lm,is)%ycom = molecule_list_igas(rand_igas,is)%ycom
        molecule_list(lm,is)%zcom = molecule_list_igas(rand_igas,is)%zcom

        atom_list(:,lm,is)%rxp = atom_list_igas(:,rand_igas,is)%rxp
        atom_list(:,lm,is)%ryp = atom_list_igas(:,rand_igas,is)%ryp
        atom_list(:,lm,is)%rzp = atom_list_igas(:,rand_igas,is)%rzp

     END IF   

     ! Randomize the molecule's orientation.
     
     CALL Rotate_Molecule_Eulerian(lm,is)
     
     ! Randomize the molecule's COM position anywhere in the box.

     IF ( box_list(this_box)%int_box_shape == int_cubic ) THEN
        
        molecule_list(lm,is)%xcom = &
                             (rranf() - 0.5_DP) * box_list(this_box)%length(1,1)
        molecule_list(lm,is)%ycom = &
                             (rranf() - 0.5_DP) * box_list(this_box)%length(2,2)
        molecule_list(lm,is)%zcom = &
                             (rranf() - 0.5_DP) * box_list(this_box)%length(3,3)
        
     END IF
     
     ! Move the molecule to the chosen COM position

     IF(species_list(is)%int_insert == int_random) THEN
     
        dx = molecule_list(lm,is)%xcom - species_list(is)%xcom
        dy = molecule_list(lm,is)%ycom - species_list(is)%ycom
        dz = molecule_list(lm,is)%zcom - species_list(is)%zcom

     ELSE

        dx = molecule_list(lm,is)%xcom &
           - molecule_list_igas(rand_igas,is)%xcom
        dy = molecule_list(lm,is)%ycom &
           - molecule_list_igas(rand_igas,is)%ycom
        dz = molecule_list(lm,is)%zcom &
           - molecule_list_igas(rand_igas,is)%zcom

     END IF        
     
     atom_list(:,lm,is)%rxp = atom_list(:,lm,is)%rxp + dx
     atom_list(:,lm,is)%ryp = atom_list(:,lm,is)%ryp + dy
     atom_list(:,lm,is)%rzp = atom_list(:,lm,is)%rzp + dz

  END IF

  !*****************************************************************************
  ! Step 3) Calculate the change in potential energy if the molecule is inserted
  !*****************************************************************************
  !
  ! Whether the insertion will be accepted depends on the change in potential
  ! energy, delta_e. The potential energy will be computed in 5 stages:
  !   3.1) Nonbonded energies
  !   3.2) Reject the move if there is any core overlap
  !   3.3) Bonded intramolecular energies
  !   3.4) Ewald energies
  !   3.5) Long-range energy correction
  ! 
  ! 3.1) Nonbonded energies
  ! If the inserted molecule overlaps cores with any other molecule in  
  ! the system, the move will be rejected. If the inserted molecule was grown 
  ! via CBMC and all attempted configurations overlapped cores, the cbmc_overlap
  ! flag equals .TRUE. If the molecule was inserted randomly, we still need to 
  ! detect core overlaps.

  IF (.NOT. cbmc_overlap) THEN

    ! Molecule COM may be outside the box boundary if grown via CBMC, so wrap
    ! the molecule coordinates back in the box (if needed)
    CALL Fold_Molecule(lm,is,this_box)

    ! Recompute the COM in case the molecule was wrapped
    CALL Get_COM(lm,is)

    ! Compute the distance of the atom farthest from COM
    CALL Compute_Max_COM_Distance(lm,is)

    ! Calculate the potential energy interaction between the inserted molecule
    ! and the rest of the system
    CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is, &
            E_inter_vdw,E_inter_qq,inter_overlap)

    ! Calculate the nonbonded energy interaction within the inserted molecule
    CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is, &
            E_intra_vdw,E_intra_qq,E_periodic_qq,intra_overlap)
    E_inter_qq = E_inter_qq + E_periodic_qq
 
  END IF

  ! 3.3) Reject the move if there is any core overlap
  IF (cbmc_overlap .OR. inter_overlap .OR. intra_overlap) THEN
     molecule_list(lm,is)%live = .FALSE.
     atom_list(:,lm,is)%exist = .FALSE.
     
     RETURN
  END IF

  ! There are no overlaps, so we can calculate the change in potential energy.
  !
  ! Already have the change in nonbonded energies

  delta_e = E_inter_vdw + E_inter_qq 
  delta_e = delta_e + E_intra_vdw + E_intra_qq

  ! 3.4) Bonded intramolecular energies
  ! If the molecule was grown via CBMC, we already have the intramolecular 
  ! bond energies? Otherwise we need to compute them.

  IF(species_list(is)%int_insert == int_random) THEN

     CALL Compute_Molecule_Bond_Energy(lm,is,E_bond)
     CALL Compute_Molecule_Angle_Energy(lm,is,E_angle)
     CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihedral)
     CALL Compute_Molecule_Improper_Energy(lm,is,E_improper)

  ELSE IF (species_list(is)%int_insert == int_igas) THEN

     E_bond = energy_igas(rand_igas,is)%bond
     E_angle = energy_igas(rand_igas,is)%angle
     E_dihedral = energy_igas(rand_igas,is)%dihedral
     E_improper = energy_igas(rand_igas,is)%improper

  END IF

  delta_e = delta_e + E_bond + E_angle + E_dihedral + E_improper

  ! 3.5) Ewald energies

  IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. &
       has_charge(is) ) THEN
 
     CALL Update_System_Ewald_Reciprocal_Energy(lm,is,this_box, &
             int_insertion,E_reciprocal)
     CALL Compute_Molecule_Ewald_Self_Energy(lm,is,this_box,E_self)

     delta_e = delta_e + E_self &
                       + (E_reciprocal - energy(this_box)%ewald_reciprocal)

  END IF

  ! 3.6) Long-range energy correction

  IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN

     ! increase number of integer beads
     nbeads_in = nint_beads(:,this_box)

     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,this_box) = nint_beads(i_type,this_box) + 1
     END DO

     CALL Compute_LR_correction(this_box,E_lrc)
     delta_e = delta_e + E_lrc - energy(this_box)%lrc

  END IF

  !*****************************************************************************
  ! Step 4) Accept or reject the move
  !*****************************************************************************
  !
  ! The following quantity is calculated
  !
  !                  p_m a_mn 
  !    ln_pacc = Log[--------]
  !                  p_n a_nm
  !
  ! and passed to accept_or_reject() which executes the metropolis criterion.
  ! The acceptance criterion to insert a molecule via CBMC is
  !
  !                                            P_seq P_bias (N + 1) Lambda^3
  !    ln_pacc = b(dU_mn-U_frag) - b mu' + Log[-----------------------------]  
  !                                                         V
  !
  !                                    P_seq P_bias (N + 1) 
  !            = b(dU_mn-U_frag) + Log[--------------------]
  !                                           b f' V
  !
  ! where the primes (') indicate that additional intensive terms have been
  ! absorbed into the chemical potential and fugacity, respectively.

  ! Compute the acceptance criterion

  IF(species_list(is)%int_insert == int_igas) THEN 
     ln_pacc = beta(this_box) * (delta_e - energy_igas(rand_igas,is)%total)
  ELSEIF (species_list(is)%fragment) THEN
     ln_pacc = beta(this_box) * (delta_e - E_angle - nrg_ring_frag_tot)
  ELSE
     ln_pacc = beta(this_box) * delta_e
  END IF

  ! P_seq and P_bias equal 1.0 unless changed by Build_Molecule.
  ln_pacc = ln_pacc + DLOG(P_seq * P_bias) &
                    + DLOG(REAL(nmols(is,this_box),DP)) &
                    - DLOG(box_list(this_box)%volume) 

  IF(lchempot) THEN
     ! chemical potential is input
     ln_pacc = ln_pacc - species_list(is)%chem_potential * beta(this_box) &
                       + 3.0_DP*DLOG(species_list(is)%de_broglie(this_box))
  ELSE
     ! fugacity is input
     ln_pacc = ln_pacc - DLOG(species_list(is)%fugacity) &
                       - DLOG(beta(this_box)) 
  END IF
  
  accept = accept_or_reject(ln_pacc)
  
  IF (accept) THEN
     ! accept the insertion

     ! number of molecules already incremented

     ! update the energies
     energy(this_box)%total = energy(this_box)%total + delta_e
     energy(this_box)%intra = energy(this_box)%intra + E_bond + E_angle &
                            + E_dihedral + E_improper
     energy(this_box)%bond = energy(this_box)%bond + E_bond
     energy(this_box)%angle = energy(this_box)%angle + E_angle
     energy(this_box)%dihedral = energy(this_box)%dihedral + E_dihedral
     energy(this_box)%improper = energy(this_box)%improper + E_improper
     energy(this_box)%intra_vdw = energy(this_box)%intra_vdw + E_intra_vdw
     energy(this_box)%intra_q = energy(this_box)%intra_q + E_intra_qq
     energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + E_inter_vdw
     energy(this_box)%inter_q = energy(this_box)%inter_q + E_inter_qq

     IF ( int_charge_sum_style(this_box) == charge_ewald .AND. &
          has_charge(is)) THEN
        energy(this_box)%ewald_reciprocal = E_reciprocal
        energy(this_box)%ewald_self = energy(this_box)%ewald_self + E_self
     END IF

     IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
        energy(this_box)%lrc = E_lrc
     END IF

     ! Increment counter
     nsuccess(is,this_box)%insertion = nsuccess(is,this_box)%insertion + 1

  ELSE
     ! reject the insertion 
     locate(nmols(is,this_box),is,this_box) = 0
     nmols(is,this_box) = nmols(is,this_box) - 1
     molecule_list(lm,is)%live = .FALSE.
     atom_list(:,lm,is)%exist = .FALSE.
     molecule_list(lm,is)%molecule_type = int_none

     ! move locate to the list of unused locates
     nmols(is,0) = nmols(is,0) + 1
     locate(nmols(is,0),is,0) = lm
     
     IF ( int_charge_sum_style(this_box) == charge_ewald .AND. &
          has_charge(is) ) THEN
        ! Restore cos_sum and sin_sum. Note that these were changed when the
        ! difference in reciprocal energies was computed.
        cos_sum(:,this_box) = cos_sum_old(:,this_box)
        sin_sum(:,this_box) = sin_sum_old(:,this_box)
     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail ) THEN
        ! Restore the total number of bead types
        nint_beads(:,this_box) = nbeads_in(:)
     END IF

  END IF

END SUBROUTINE Insertion
