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

SUBROUTINE Insertion(this_box,mcstep,randno)

  !*************************************************************************
  ! 
  ! The subroutine inserts a molecule in the system. Various configurational
  ! schemes will be enabled as the code is verified.
  !
  ! Called by
  !
  !    gcmc_driver
  !
  ! Revision history
  !
  !   12/10/13 : Beta version
  !
  !*******************************************************************************

  USE Run_Variables
  USE Energy_Routines
  USE IO_Utilities
  USE Random_Generators
  USE Rotation_Routines
  USE Fragment_Growth

  IMPLICIT NONE

  INTEGER :: is, is_1, nmolecules_is, i, alive, this_box, nmol_box, which_cell, i_type
  INTEGER :: nmols_sorbate, nsorbate, is_counter, rand_igas
  INTEGER, ALLOCATABLE :: sorbate_id(:), frag_order(:)
  INTEGER :: which_anchor,mcstep, ktothen, tot_mols, ifrag

  REAL(DP) :: dx, dy, dz, E_bond, E_angle, E_dihedral, E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq, E_improper
  REAL(DP) :: delta_e, E_reciprocal_move, E_self_move, E_lrc

  REAL(DP) :: factor, alpha_ratio, pick_species, P_forward, randno
  REAL(DP) :: gnew, gold, dg  ! weighting factors for the new mol number and old mol number
  REAL(DP), ALLOCATABLE :: sorbate_x(:)
  REAL(DP) :: this_lambda

  REAL(DP) :: nrg_ring_frag_tot
  LOGICAL :: inter_overlap, accept, accept_or_reject, cbmc_overlap, intra_overlap

  REAL(DP) :: checke, energy_old, energy_change
  LOGICAL  :: superbad

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  pacc = 0.0_DP
  paccbiased = 0.0_DP
  alpha_ratio = 1.0_DP
  inter_overlap = .FALSE.
  cbmc_overlap = .FALSE.
  P_forward = 1.0_DP
  nrg_ring_frag_tot = 0.0_DP

  IF (int_sim_type .NE. sim_gemc) THEN
     this_box = INT(rranf() * nbr_boxes) + 1
     ! else for GEMC simulation the box will be specified as input
  END IF

  tot_trials(this_box) = tot_trials(this_box) + 1  
  energy_change = energy(1)%inter_vdw

  ! Choose a species to insert

  is_1 = INT(rranf() * nspec_insert) + 1
  is_counter = 0

  DO is = 1, nspecies
     IF(species_list(is)%int_species_type == int_sorbate) is_counter = is_counter + 1
     IF(is_counter == is_1) EXIT
  END DO

  nmolecules_is = 0

!!$  IF (int_sim_type == sim_gcmc) THEN
!!$     IF ( is == 1 ) t_collect = .TRUE.
!!$  END IF

  tot_mols = SUM(nmols(is,:))

  IF (tot_mols == (nmolecules(is))) RETURN
  ! exit if we are attempting an insertion above maximum number

  ntrials(is,this_box)%insertion = ntrials(is,this_box)%insertion + 1
  tot_trials(this_box) = tot_trials(this_box) + 1

 ! Assign a locate number for this molecule

  IF ( locate(tot_mols+1,is) == 0 ) THEN
     locate(tot_mols+1,is) = tot_mols + 1
     ! otherwise we will use the locate number of a previously deleted molecule that has 
     ! been moved to the end of the array.
  END IF

  alive = locate(tot_mols+1,is)
  molecule_list(alive,is)%which_box = this_box
  this_lambda = 1.0_DP
  molecule_list(alive,is)%cfc_lambda = this_lambda
  molecule_list(alive,is)%molecule_type = int_normal

  ! Randomly insert the COM in the simulation box.
  
  IF (species_list(is)%fragment .AND. (species_list(is)%int_insert & 
       .NE. int_igas) ) THEN

     del_Flag = .FALSE.
     get_fragorder = .TRUE.
     ALLOCATE(frag_order(nfragments(is)))
     CALL Build_Molecule(alive,is,this_box,frag_order,this_lambda,which_anchor,P_forward, nrg_ring_frag_tot, cbmc_overlap)
     molecule_list(alive,is)%live = .TRUE.
     atom_list(:,alive,is)%exist = .TRUE.

     ktothen = 1

     IF (nfragments(is) /= 0 ) THEN

        ktothen = ktothen * kappa_ins

        IF (kappa_rot /= 0 ) THEN

           ktothen = ktothen * kappa_rot

        END IF

        IF (kappa_dih /= 0 ) THEN

           DO ifrag = 1, nfragments(is) - 1

              ktothen = ktothen * kappa_dih
              
           END DO

        END IF

     END IF

     P_forward = P_forward * REAL(ktothen, DP)

  ELSE
 
     atom_list(:,alive,is)%exist = .true.
     molecule_list(alive,is)%live = .true.
     IF(species_list(is)%int_insert == int_random) THEN
     
        ! COM of the species from the initial configuration is 

        molecule_list(alive,is)%xcom = species_list(is)%xcom
        molecule_list(alive,is)%ycom = species_list(is)%ycom
        molecule_list(alive,is)%zcom = species_list(is)%zcom

        
        atom_list(:,alive,is)%rxp = init_list(:,1,is)%rxp
        atom_list(:,alive,is)%ryp = init_list(:,1,is)%ryp
        atom_list(:,alive,is)%rzp = init_list(:,1,is)%rzp


     ELSE IF(species_list(is)%int_insert == int_igas) THEN

        rand_igas = (rranf() * n_igas(is)) + 1

        molecule_list(alive,is)%xcom = molecule_list_igas(rand_igas,is)%xcom
        molecule_list(alive,is)%ycom = molecule_list_igas(rand_igas,is)%ycom
        molecule_list(alive,is)%zcom = molecule_list_igas(rand_igas,is)%zcom

        atom_list(:,alive,is)%rxp = atom_list_igas(:,rand_igas,is)%rxp
        atom_list(:,alive,is)%ryp = atom_list_igas(:,rand_igas,is)%ryp
        atom_list(:,alive,is)%rzp = atom_list_igas(:,rand_igas,is)%rzp

     END IF   
     
     CALL Rotate_Molecule_Eulerian(alive,is)
     
     IF ( box_list(this_box)%int_box_shape == int_cubic ) THEN
        
        molecule_list(alive,is)%xcom = (rranf() - 0.5_DP) * box_list(this_box)%length(1,1)
        molecule_list(alive,is)%ycom = (rranf() - 0.5_DP) * box_list(this_box)%length(2,2)
        molecule_list(alive,is)%zcom = (rranf() - 0.5_DP) * box_list(this_box)%length(3,3)
        
     END IF
     
     ! Coordinates obtained are for the initial coordinates so translate it to the current position

     IF(species_list(is)%int_insert == int_random) THEN
     
        dx = molecule_list(alive,is)%xcom - species_list(is)%xcom
        dy = molecule_list(alive,is)%ycom - species_list(is)%ycom
        dz = molecule_list(alive,is)%zcom - species_list(is)%zcom

     ELSE

        dx = molecule_list(alive,is)%xcom - molecule_list_igas(rand_igas,is)%xcom
        dy = molecule_list(alive,is)%ycom - molecule_list_igas(rand_igas,is)%ycom
        dz = molecule_list(alive,is)%zcom - molecule_list_igas(rand_igas,is)%zcom

     END IF        
     
     atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp + dx
     atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp + dy
     atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp + dz

     ! Now compute the energy difference due to the insertion of this molecules
     ! make this molecule alive

     molecule_list(alive,is)%live = .TRUE.
     molecule_list(alive,is)%which_box = this_box
  
     atom_list(:,alive,is)%exist = .TRUE.
     molecule_list(alive,is)%cfc_lambda = 1.0_DP

  END IF

  ! compute the distance of atom farthest from COM
  IF(.NOT. cbmc_overlap) THEN

    CALL Fold_Molecule(alive,is,this_box)
    CALL Get_COM(alive,is)
    CALL Compute_Max_COM_Distance(alive,is)

    ! Intra molecule energy

     delta_e = 0.0_DP
     
     E_inter_vdw = 0.0_DP
     E_inter_qq = 0.0_DP
     
     CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_inter_vdw,E_inter_qq,inter_overlap)

  END IF

  IF (inter_overlap .OR. cbmc_overlap) THEN

     ! reject the move immediately
     molecule_list(alive,is)%live = .FALSE.
     atom_list(:,alive,is)%exist = .FALSE.
    
     RETURN

  END IF

  delta_e = E_inter_vdw + E_inter_qq 

  IF(species_list(is)%int_insert == int_random) THEN

     CALL Compute_Molecule_Bond_Energy(alive,is,E_bond)
     CALL Compute_Molecule_Angle_Energy(alive,is,E_angle)
     CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihedral)
     CALL Compute_Molecule_improper_Energy(alive,is,E_improper)

  ELSE IF (species_list(is)%int_insert == int_igas) THEN

     E_bond = energy_igas(rand_igas,is)%bond
     E_angle = energy_igas(rand_igas,is)%angle
     E_dihedral = energy_igas(rand_igas,is)%dihedral
     E_improper = energy_igas(rand_igas,is)%improper

  END IF

  delta_e = delta_e + E_bond + E_angle + E_dihedral + E_improper

  CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw,E_intra_qq,intra_overlap)
  delta_e = delta_e + E_intra_vdw + E_intra_qq

  IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is)) ) THEN

     CALL Compute_Ewald_Reciprocal_Energy_Difference(alive,alive,is,this_box,int_insertion,E_reciprocal_move)
     CALL Compute_Ewald_Self_Energy_Difference(alive,is,this_box,int_insertion,E_self_move)
     delta_e = delta_e + E_self_move
     delta_e = delta_e + E_reciprocal_move - energy(this_box)%ewald_reciprocal

  END IF

  IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
     ! increase number of integer beads

     nbeads_in = nint_beads(:,this_box)

     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,this_box) = nint_beads(i_type,this_box) + 1
     END DO

     CALL Compute_LR_correction(this_box,e_lrc)
     delta_e = delta_e + e_lrc - energy(this_box)%lrc

  END IF

  ! calculate weighting function and apply acceptance rule
  
    dg = 0.0_DP

  IF(species_list(is)%int_insert == int_igas) THEN 
     pacc = beta(this_box) * (delta_e - energy_igas(rand_igas,is)%total)
  ELSEIF (species_list(is)%fragment) THEN
     pacc = beta(this_box) * (delta_e - E_angle - nrg_ring_frag_tot)
  ELSE
     pacc = beta(this_box)*delta_e
  END IF

  pacc = pacc - DLOG(alpha_ratio) + DLOG(P_forward) - DLOG(species_list(is)%zig_by_omega) &
       +  DLOG(REAL(nmols(is,this_box)+1,DP)) 


  IF(lchempot) THEN
     ! chemical potential is input
     pacc = pacc - species_list(is)%chem_potential * beta(this_box)
  ELSE
     ! user input is fugacity
     pacc = pacc - DLOG(species_list(is)%fugacity * beta(this_box) * box_list(this_box)%volume)
  END IF

!!$  pacc = beta(this_box)* (delta_e) - DLOG(box_list(this_box)%volume) - &
!!$       DLOG(species_list(is)%activity) + DLOG(REAL(nmols(is,this_box) + 1, DP))
!!$
!!$  IF (species_list(is)%fragment) THEN
!!$     pacc = pacc - beta(this_box) * (E_angle - nrg_ring_frag_tot)
!!$  END IF
!!$
!!$  pacc = pacc +  DLOG(kappa_ins * P_forward)

  ! correct for ring biasing if any

  factor = pacc + dg  ! accept based on probability times weighting factor
 ! factor = pacc 
  
  accept = accept_or_reject(factor)
  
  factor = -factor
  IF( factor < 0.0_DP ) THEN
    factor = DEXP(factor)
  ELSE
    factor = 1.0_DP
  END IF

  paccbiased = factor

  IF (accept) THEN
     ! update the number of molecules
     nmols(is,this_box) = nmols(is,this_box) + 1
     ! update the energies
     energy(this_box)%total = energy(this_box)%total + delta_e
     energy(this_box)%intra = energy(this_box)%intra + E_bond + E_angle + E_dihedral + E_improper
     energy(this_box)%bond = energy(this_box)%bond + E_bond
     energy(this_box)%angle = energy(this_box)%angle + E_angle
     energy(this_box)%dihedral = energy(this_box)%dihedral + E_dihedral
     energy(this_box)%improper = energy(this_box)%improper + E_improper
     energy(this_box)%intra_vdw = energy(this_box)%intra_vdw + E_intra_vdw
     energy(this_box)%intra_q = energy(this_box)%intra_q + E_intra_qq
     energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + E_inter_vdw
     energy(this_box)%inter_q = energy(this_box)%inter_q + E_inter_qq

     IF ( int_charge_sum_style(this_box) == charge_ewald .AND. has_charge(is)) THEN

        energy(this_box)%ewald_reciprocal = E_reciprocal_move
        energy(this_box)%ewald_self = energy(this_box)%ewald_self + E_self_move
     END IF

     IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
        energy(this_box)%lrc = e_lrc
     END IF

     nsuccess(is,this_box)%insertion = nsuccess(is,this_box)%insertion + 1


  ELSE
  
     molecule_list(alive,is)%live = .FALSE.
     atom_list(:,alive,is)%exist = .FALSE.
     molecule_list(alive,is)%molecule_type = int_none
     
     IF ( (int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is)) ) THEN
        ! restore cos_sum and sin_sum. Note that these were changed when difference in
        ! reciprocal energies was computed

        cos_sum(:,this_box) = cos_sum_old(:,this_box)
        sin_sum(:,this_box) = sin_sum_old(:,this_box)

     END IF

     IF ( int_vdw_sum_style(this_box) == vdw_cut_tail ) THEN

        nint_beads(:,this_box) = nbeads_in(:)

     END IF
     
  END IF

  ! Accept or reject the move.
  ! also, update various arrays that depend on the number of molecules
  pacc = -pacc
  IF( pacc < 0.0_DP ) THEN
    pacc = DEXP(pacc)
  ELSE
    pacc = 1.0_DP
  END IF

END SUBROUTINE Insertion
