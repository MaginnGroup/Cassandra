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

SUBROUTINE Insertion

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
  INTEGER :: ibox ! attempt to insert a molecule in ibox

  ! Local declarations
  INTEGER :: i, i_type               ! atom indices
  INTEGER :: ifrag                   ! fragment indices
  INTEGER :: im                      ! molecule INDEX
  INTEGER :: lm                      ! molecule LOCATE
  INTEGER :: is, is_rand, is_counter ! species indices
  INTEGER :: which_anchor
  INTEGER, ALLOCATABLE :: frag_order(:)
  INTEGER :: tot_mols, mcstep

  REAL(DP) :: dx, dy, dz, dE, dE_frag
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq, E_periodic_qq
  REAL(DP) :: E_reciprocal, E_self, E_lrc
  REAL(DP) :: nrg_ring_frag_tot
  REAL(DP) :: ln_pacc, ln_pseq, ln_pbias, this_lambda

  LOGICAL :: inter_overlap, cbmc_overlap, intra_overlap
  LOGICAL :: accept_or_reject

  ! Initialize variables
  ln_pacc = 0.0_DP
  ln_pseq = 0.0_DP
  ln_pbias = 0.0_DP
  this_lambda = 1.0_DP
  nrg_ring_frag_tot = 0.0_DP
  inter_overlap = .FALSE.
  cbmc_overlap = .FALSE.
  intra_overlap = .FALSE.
  ibox = 1
  accept = .FALSE.

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
  ntrials(is,ibox)%insertion = ntrials(is,ibox)%insertion + 1
  tot_trials(ibox) = tot_trials(ibox) + 1

  !  * Assign a LOCATE for this molecule from the list of unused LOCATEs
  nmols(is,ibox) = nmols(is,ibox) + 1
  locate(nmols(is,ibox),is,ibox) = locate(nmols(is,0),is,0)
  locate(nmols(is,0),is,0) = 0
  nmols(is,0) = nmols(is,0) - 1

  !  * Set properties of the to-be-inserted molecule
  im = nmols(is,ibox)
  lm = locate(im,is,ibox)
  molecule_list(lm,is)%which_box = ibox
  molecule_list(lm,is)%frac = this_lambda
  molecule_list(lm,is)%molecule_type = int_normal

  ! With the bookkeeping completed, we are ready to attempt the insertion
  
  !*****************************************************************************
  ! Step 2) Choose a position, orientation and conformation for the 
  !         to-be-inserted molecule
  !*****************************************************************************
  !
  ! Build_Molecule places the first fragment, then calls Fragment_Placement
  ! to place the additional fragments 
  del_flag = .FALSE.     ! Change the coordinates of 'lm'
  get_fragorder = .TRUE.
  ALLOCATE(frag_order(nfragments(is)))
  CALL Build_Molecule(lm,is,ibox,frag_order,this_lambda, &
          ln_pseq,ln_pbias,nrg_ring_frag_tot,cbmc_overlap)
  DEALLOCATE(frag_order)

  ! Turn the molecule on
  molecule_list(lm,is)%live = .TRUE.
  atom_list(:,lm,is)%exist = .TRUE.

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
  ! Step 3) Calculate the change in potential energy if the molecule is inserted
  !*****************************************************************************
  !
  ! Whether the insertion will be accepted depends on the change in potential
  ! energy, dE. The potential energy will be computed in 5 stages:
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
    CALL Fold_Molecule(lm,is,ibox)

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
     ! reject the insertion 
     locate(nmols(is,ibox),is,ibox) = 0
     nmols(is,ibox) = nmols(is,ibox) - 1
     molecule_list(lm,is)%live = .FALSE.
     atom_list(:,lm,is)%exist = .FALSE.
     molecule_list(lm,is)%molecule_type = int_none

     ! move locate to the list of unused locates
     nmols(is,0) = nmols(is,0) + 1
     locate(nmols(is,0),is,0) = lm
     
     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,A9)') &
             i_mcstep, 'insert' , lm, is, ibox, .FALSE., 'overlap'
     END IF

     RETURN
  END IF

  ! There are no overlaps, so we can calculate the change in potential energy.
  !
  ! Already have the change in nonbonded energies

  dE = E_inter_vdw + E_inter_qq 
  dE = dE + E_intra_vdw + E_intra_qq

  ! 3.4) Bonded intramolecular energies
  ! If the molecule was grown via CBMC, we already have the intramolecular 
  ! bond energies? Otherwise we need to compute them.

  CALL Compute_Molecule_Bond_Energy(lm,is,E_bond)
  CALL Compute_Molecule_Angle_Energy(lm,is,E_angle)
  CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihedral)
  CALL Compute_Molecule_Improper_Energy(lm,is,E_improper)

  dE = dE + E_bond + E_angle + E_dihedral + E_improper

  ! 3.5) Ewald energies

  IF (int_charge_style(ibox) == charge_coul) THEN

        IF ( (int_charge_sum_style(ibox) == charge_ewald) .AND. &
             has_charge(is) ) THEN
       
           CALL Update_System_Ewald_Reciprocal_Energy(lm,is,ibox, &
                   int_insertion,E_reciprocal)

            dE = dE + (E_reciprocal - energy(ibox)%ewald_reciprocal)
           
        END IF

        CALL Compute_Molecule_Self_Energy(lm,is,ibox,E_self)

        dE = dE + E_self

  END IF

  ! 3.6) Long-range energy correction

  IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN

     ! increase number of integer beads
     nbeads_in = nint_beads(:,ibox)

     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,ibox) = nint_beads(i_type,ibox) + 1
     END DO

     CALL Compute_LR_correction(ibox,E_lrc)
     dE = dE + E_lrc - energy(ibox)%lrc

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
  !                                                       (N + 1) Lambda^3
  !    ln_pacc = b(dU_mn-U_frag) - b mu' + ln_pbias + Log[-----------------]  
  !                                                              V
  !
  ! where the primes (') indicate that additional intensive terms have been
  ! absorbed into the chemical potential.

  ! Compute the acceptance criterion

  ! change in energy less energy used to bias selection of fragments
  dE_frag = E_angle + nrg_ring_frag_tot
  ln_pacc = beta(ibox) * (dE - dE_frag)

  ! chemical potential
  ln_pacc = ln_pacc - species_list(is)%chem_potential * beta(ibox)

  ! bias from CBMC
  ln_pacc = ln_pacc + ln_pbias

  ! density
  ln_pacc = ln_pacc + DLOG(REAL(nmols(is,ibox),DP)) &
                    + 3.0_DP*DLOG(species_list(is)%de_broglie(ibox)) &
                    - DLOG(box_list(ibox)%volume) 
  
  accept = accept_or_reject(ln_pacc)
  
  IF (accept) THEN
     ! accept the insertion

     ! number of molecules already incremented

     ! update the energies
     energy(ibox)%total = energy(ibox)%total + dE
     energy(ibox)%intra = energy(ibox)%intra + E_bond + E_angle &
                            + E_dihedral + E_improper
     energy(ibox)%bond = energy(ibox)%bond + E_bond
     energy(ibox)%angle = energy(ibox)%angle + E_angle
     energy(ibox)%dihedral = energy(ibox)%dihedral + E_dihedral
     energy(ibox)%improper = energy(ibox)%improper + E_improper
     energy(ibox)%intra_vdw = energy(ibox)%intra_vdw + E_intra_vdw
     energy(ibox)%intra_q = energy(ibox)%intra_q + E_intra_qq
     energy(ibox)%inter_vdw = energy(ibox)%inter_vdw + E_inter_vdw
     energy(ibox)%inter_q = energy(ibox)%inter_q + E_inter_qq

     IF (int_charge_style(ibox) == charge_coul) THEN

         IF ( int_charge_sum_style(ibox) == charge_ewald .AND. &
              has_charge(is)) THEN
            energy(ibox)%ewald_reciprocal = E_reciprocal
         END IF

         energy(ibox)%self = energy(ibox)%self + E_self

     END IF

     IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN
        energy(ibox)%lrc = E_lrc
     END IF

     ! Increment counter
     nsuccess(is,ibox)%insertion = nsuccess(is,ibox)%insertion + 1

  ELSE
     ! reject the insertion 
     locate(nmols(is,ibox),is,ibox) = 0
     nmols(is,ibox) = nmols(is,ibox) - 1
     molecule_list(lm,is)%live = .FALSE.
     atom_list(:,lm,is)%exist = .FALSE.
     molecule_list(lm,is)%molecule_type = int_none

     ! move locate to the list of unused locates
     nmols(is,0) = nmols(is,0) + 1
     locate(nmols(is,0),is,0) = lm
     
     IF ( int_charge_sum_style(ibox) == charge_ewald .AND. &
          has_charge(is) ) THEN
        ! Restore cos_sum and sin_sum. Note that these were changed when the
        ! difference in reciprocal energies was computed.
        cos_sum(:,ibox) = cos_sum_old(:,ibox)
        sin_sum(:,ibox) = sin_sum_old(:,ibox)
     END IF

     IF ( int_vdw_sum_style(ibox) == vdw_cut_tail ) THEN
        ! Restore the total number of bead types
        nint_beads(:,ibox) = nbeads_in(:)
     END IF

  END IF

  IF (verbose_log) THEN
     WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,F9.3)') &
           i_mcstep, 'insert' , lm, is, ibox, accept, ln_pacc
  END IF

END SUBROUTINE Insertion
