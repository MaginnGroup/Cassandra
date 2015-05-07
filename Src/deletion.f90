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
  
  !********************************************************************************
  ! This subroutine deletes a molecule a box. Thexs box is first chosen at random
  ! and then an adsorbate component is selected at random. Next a randomly chosen
  ! molecule of this component is deleted from the box
  !
  ! Called by
  !
  !  gcmc_driver
  ! 
  ! Revision History
  !
  !   12/10/13 : Beta Release
  !
  !
  !********************************************************************************
  
  USE Type_Definitions
  USE Run_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE Fragment_Growth

  IMPLICIT NONE


  INTEGER, INTENT(INOUT) :: this_box

  INTEGER :: is, im, alive, k, position, i_type, i, nsorbate, nmols_sorbate, which_anchor
  INTEGER :: is_1, is_counter, mcstep, ktothen, ifrag
  INTEGER, ALLOCATABLE :: sorbate_id(:), frag_order(:)

  REAL(DP) :: delta_e, delta_e_pacc, delta_v, E_bond, E_angle, E_dihedral, E_intra_vdw, E_intra_qq, E_inter_vdw
  REAL(DP) :: E_inter_qq, E_improper
  REAL(DP) :: E_reciprocal_move, E_self_move, e_lrc, alpha_ratio, factor
  REAL(DP) :: gnew, gold, dg ! weighting factors for the new mol number and old mol number
  REAL(DP), ALLOCATABLE :: sorbate_x(:)
  REAL(DP) :: pick_species, P_reverse, this_lambda, nrg_ring_frag_tot, randno
  REAL(DP) :: E_intra_vdw_igas, E_intra_qq_igas

  LOGICAL :: inter_overlap, accept, accept_or_reject, cbmc_overlap, intra_overlap


  REAL(DP) :: check_e, energy_old
  LOGICAL  :: superbad

  
  pacc = 0.0_DP
  paccbiased = 0.0_DP
  alpha_ratio = 1.0_DP
  nrg_ring_frag_tot = 0.0_DP
  

  tot_trials(this_box) = tot_trials(this_box) + 1

  is_1 = INT(rranf() * nspec_insert) + 1
  is_counter = 0

  DO is = 1, nspecies
     IF(species_list(is)%int_species_type == int_sorbate) is_counter = is_counter + 1
     IF(is_counter == is_1) EXIT
  END DO


  IF(nmols(is,this_box) == 0) RETURN

  ntrials(is,this_box)%deletion = ntrials(is,this_box)%deletion + 1  
  ! Determine the index of im

  im = INT(rranf() * nmols(is,this_box)) + 1

  CALL Get_Index_Molecule(this_box,is,im,alive)

  CALL Save_Old_Cartesian_Coordinates(alive,is)  

  ! Compute the energy of the molecule
  
  ! Intra molecule energy
  
  delta_e = 0.0_DP
  P_reverse = 1.0_DP

  this_lambda = 1.0_DP

  cbmc_overlap = .FALSE.


  IF(species_list(is)%fragment .AND. (species_list(is)%int_insert .NE. int_igas)) THEN
     del_Flag = .TRUE.
     get_fragorder = .TRUE.
     ALLOCATE(frag_order(nfragments(is)))


     CALL Build_Molecule(alive,is,this_box,frag_order,this_lambda, which_anchor, P_reverse, nrg_ring_frag_tot, cbmc_overlap)


     DEALLOCATE(frag_order)
     
     IF (cbmc_overlap) THEN

! Revert to the COM and Eulerian angles for the molecule

        CALL Revert_Old_Cartesian_Coordinates(alive,is)
        atom_list(1:natoms(is),alive,is)%exist = .TRUE.
        molecule_list(alive,is)%cfc_lambda = this_lambda
        
        WRITE(*,*)
        WRITE(*,*) 'Warning....energy overlap detected in old configuration in deletion.f90'
        WRITE(*,*) 'molecule, species', alive, is
        WRITE(*,*)
     END IF

     ktothen = 1

     IF (nfragments(is) /=0 ) THEN

        ktothen = ktothen * kappa_ins

        IF (kappa_rot /= 0) THEN
           ktothen = ktothen * kappa_rot
        END IF
        
        IF (kappa_dih /=0 ) THEN

           DO ifrag = 2, nfragments(is)
              ktothen = ktothen * kappa_dih
           END DO
        END IF
     
     END IF
     
     P_reverse = P_reverse * REAL(ktothen , DP)

  END IF

  CALL Get_COM(alive,is)
  CALL Compute_Max_COM_Distance(alive,is)

  CALL Compute_Molecule_Bond_Energy(alive,is,E_bond)
  CALL Compute_Molecule_Angle_Energy(alive,is,E_angle)
  CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihedral)
  CALL Compute_Molecule_Improper_Energy(alive,is,E_improper)

  delta_e = delta_e + E_bond + E_angle + E_dihedral + E_improper  
  
  ! Nonbonded energy  
  CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw,E_intra_qq,intra_overlap)

  IF (l_pair_nrg) THEN
     CALL Store_Molecule_Pair_Interaction_Arrays(alive,is,this_box,E_inter_vdw,E_inter_qq)
  ELSE
     CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_inter_vdw,E_inter_qq,inter_overlap)
  END IF

  delta_e = delta_e + E_intra_vdw + E_intra_qq
  delta_e = delta_e + E_inter_vdw + E_inter_qq

  IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is)) ) THEN
     CALL Compute_Ewald_Reciprocal_Energy_Difference(alive,alive,is,this_box,int_deletion,E_reciprocal_move)
     CALL Compute_Ewald_Self_Energy_Difference(alive,is,this_box,int_deletion,E_self_move)
     delta_e = delta_e - E_self_move
     delta_e = delta_e - (E_reciprocal_move - energy(this_box)%ewald_reciprocal)
  END IF


  IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN

     ! subtract off beads for this species

     nbeads_out(:) = nint_beads(:,this_box)

     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,this_box) = nint_beads(i_type,this_box) - 1
     END DO
     CALL Compute_LR_correction(this_box,e_lrc)
     delta_e = delta_e - ( e_lrc - energy(this_box)%lrc )

  END IF  
  

  ! calculate the factor in the acceptance rule and the weighting function in necessary
  ! min(1,exp(-factor)). Note that the energy difference is actually (e_total_new - e_total_old)
  ! delta_e computed above is (e_total_old - e_total_new) so in the factor below
  ! - delta_e is used to represent actual change in energy. See below when energies are
  ! update upon suceessful deletion


  delta_e_pacc = delta_e

  IF(species_list(is)%int_insert == int_igas) THEN
     igas_flag = .TRUE.
     CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw_igas,E_intra_qq_igas,intra_overlap)
     igas_flag = .FALSE. 
     delta_e_pacc = delta_e_pacc - E_bond - E_angle - E_dihedral - E_improper - &
                         E_intra_vdw_igas - E_intra_qq_igas
  END IF

  IF(species_list(is)%fragment .AND. species_list(is)%int_insert .NE. int_igas) THEN
     delta_e_pacc = delta_e_pacc - E_angle - nrg_ring_frag_tot
  END IF

  pacc = beta(this_box) * (-delta_e_pacc) - DLOG(alpha_ratio) - &
         DLOG(REAL(nmols(is,this_box),DP)) - DLOG(P_reverse) 
 
  IF(lchempot) THEN
     ! chemical potential is input
     pacc = pacc + beta(this_box) * species_list(is)%chem_potential + &
       DLOG(box_list(this_box)%volume) - 3.0_DP*DLOG(species_list(is)%de_broglie(this_box)) 
  ELSE
     pacc = pacc + DLOG(species_list(is)%fugacity * beta(this_box) * box_list(this_box)%volume) &
     - DLOG(species_list(is)%zig_by_omega)

  END IF 
 
  accept = accept_or_reject(pacc)


  IF(cbmc_overlap) accept = .FALSE.

 
  if (alive==80) then
        write(*,*) 'deletion', alive, accept
  end if


  IF (accept) THEN
     ! Update energies

     energy(this_box)%total = energy(this_box)%total - delta_e
     energy(this_box)%intra = energy(this_box)%intra - E_bond - E_angle - E_dihedral
     energy(this_box)%bond = energy(this_box)%bond - E_bond
     energy(this_box)%angle = energy(this_box)%angle - E_angle
     energy(this_box)%dihedral = energy(this_box)%dihedral - E_dihedral
     energy(this_box)%improper = energy(this_box)%improper - E_improper
     energy(this_box)%intra_vdw = energy(this_box)%intra_vdw - E_intra_vdw
     energy(this_box)%intra_q = energy(this_box)%intra_q - E_intra_qq
     energy(this_box)%inter_vdw = energy(this_box)%inter_vdw - E_inter_vdw
     energy(this_box)%inter_q   = energy(this_box)%inter_q - E_inter_qq

     IF ( int_charge_sum_style(this_box) == charge_ewald .AND. has_charge(is)) THEN
        
        energy(this_box)%ewald_reciprocal = E_reciprocal_move
        energy(this_box)%ewald_self = energy(this_box)%ewald_self + E_self_move
        
     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
        energy(this_box)%lrc = e_lrc
     END IF

     ! obtain the original position of the deleted molecule so that the linked list
     ! can be updated

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

     nmols(is,this_box) = nmols(is,this_box) - 1
     nsuccess(is,this_box)%deletion = nsuccess(is,this_box)%deletion + 1

!     CALL System_Energy_Check(1,mcstep,randno)
  ELSE

     IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is)) ) THEN
        ! restore cos_sum and sin_sum. Note that these were changed when difference in
        ! reciprocal energies was computed

        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)        
        cos_sum(:,this_box) = cos_sum_old(:,this_box)
        sin_sum(:,this_box) = sin_sum_old(:,this_box)
        !$OMP END PARALLEL WORKSHARE

     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail ) THEN
        ! restore the total number of bead types
        nint_beads(:,this_box) = nbeads_out(:)
     END IF

  END IF

  IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)

END SUBROUTINE Deletion

     

