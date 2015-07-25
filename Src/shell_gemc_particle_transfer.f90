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

SUBROUTINE Shell_GEMC_Particle_Transfer(box_in, box_out)

  !********************************************************************************
  !
  ! This subroutine performs particle swaps in a GEMC ensemble when core-shell units are present in the box
  !
  !
  ! Called by
  !
  !   gemc_driver
  !*********************************************************************************

  USE Run_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE IO_Utilities
  USE Fragment_Growth
  USE Pair_Nrg_Routines
  USE Volume

  IMPLICIT NONE

  INTEGER :: box_in, box_out, i, this_species, i_type

  INTEGER :: nmols_sorbate, nsorbate, is, ibox, im ,alive, which_anchor
  INTEGER :: rand_igas 

  INTEGER, ALLOCATABLE :: sorbate_id(:), frag_order(:)

  REAL(DP), ALLOCATABLE :: sorbate_x(:)

  REAL(DP) :: pick_species, delta_e_out, delta_e_out_pacc, e_bond_out, e_angle_out, e_dihed_out
  REAL(DP) :: E_intra_vdw_out, E_intra_qq_out
  REAL(DP) :: E_inter_vdw_out, E_inter_qq_out
  REAL(DP) :: E_intra_vdw_igas, E_intra_qq_igas
  REAL(DP) :: E_reciprocal_move_out, E_self_move_out, e_lrc_out, e_improper_out
  REAL(DP) :: E_intra_vdw_in, E_intra_qq_in
  REAL(DP) :: E_inter_vdw_in, E_inter_qq_in
  REAL(DP) :: E_reciprocal_move_in, E_self_move_in, e_lrc_in, E_bond_in, E_improper_in
  REAL(DP) :: E_angle_in, E_dihed_in, delta_e_in, delta_e_in_pacc, factor, P_forward, P_reverse, potw, CP_energy, P_seq
  REAL(DP) :: lambda_for_build
  LOGICAL :: inter_overlap, accept, accept_or_reject, cbmc_overlap, forward
  LOGICAL :: intra_overlap

  REAL(DP) :: nrg_ring_frag_in, nrg_ring_frag_out

  TYPE(atom_class), ALLOCATABLE :: new_atom_list(:)
  TYPE(molecule_class) :: new_molecule_list
  TYPE(Energy_Class) :: energy_old_in, energy_old_out

 ! Variables added for l_pair_nrg and reciprocal k vector storage

  INTEGER :: position, box_number

  REAL(DP) :: time0, time1, pick_box
  REAL(DP) :: E_vdw_move_alive, E_qq_move_alive
  LOGICAL :: l_charge_in, l_charge_out, overlap
  LOGICAL :: conv_Drude1, conv_Drude2
  INTEGER :: tot_mol_box_in, tot_mol_box_out

  potw = 1.0_DP
  forward = .FALSE.
  inter_overlap = .false.
  cbmc_overlap = .false.

  P_seq = 1.0_DP
  P_forward = 1.0_DP
  P_reverse = 1.0_DP

  box_out = INT (rranf() * nbr_boxes) + 1


  pick_box = rranf()

  DO ibox = 1,nbr_boxes
     IF (ibox == box_out) CYCLE

     IF (pick_box <= prob_swap_boxes(box_out,ibox)) EXIT

  END DO

  box_in = ibox

  
  P_forward = P_forward * prob_swap_boxes(box_out,box_in)
  P_reverse = P_reverse * prob_swap_boxes(box_in,box_out)

  ! increment total number of trials for each of the boxes
  tot_trials(box_out) = tot_trials(box_out) + 1
  tot_trials(box_in) = tot_trials(box_in) + 1

  ! Choose a species based on the overall mole fraction.
  ! We will base this on the sorbate mole fraction
 
  IF (l_mol_frac_swap) THEN
     
     ! Choose a species based on the overall mole fraction.
     ! We will base this on the sorbate mole fraction
     ALLOCATE(sorbate_id(nspecies), sorbate_x(nspecies))
     
     nmols_sorbate = 0
     nsorbate = 0
     sorbate_x(:) = 0.0_DP
     
     DO is = 1, nspecies
           nsorbate = nsorbate + 1
           sorbate_id(nsorbate) = is
           nmols_sorbate = nmols_sorbate + SUM(nmols(is,:))
           sorbate_x(nsorbate) = REAL(SUM(nmols(is,:)),DP)
     END DO
  
     sorbate_x(:) = sorbate_x(:) / nmols_sorbate
     
     DO i = 2, nsorbate
        sorbate_x(i) = sorbate_x(i) + sorbate_x(i-1)
     END DO
     
     ! Now pick a species at random using Golden sampling
     
     pick_species = rranf()
     
     DO i = 1, nsorbate
        IF (pick_species < sorbate_x(i)) EXIT
     END DO
     
     this_species = sorbate_id(i)
     
  ELSE
     
     ! pick the species based on specified probability
     
     pick_species = rranf()
     
     DO is = 1, nspecies
        IF (pick_species <= prob_swap_species(is)) EXIT
     END DO
     
     this_species = is
     
  END IF


  IF(nmols(this_species,box_out) == 0) THEN
 !    IF(cpcollect)  CALL Chempot(box_in)
     RETURN
  END IF
  
  ntrials(this_species,box_out)%deletion = ntrials(this_species,box_out)%deletion + 1
  ntrials(this_species,box_in)%insertion = ntrials(this_species,box_in)%insertion + 1
  ntrials(this_species,box_in)%cpcalc = ntrials(this_species,box_in)%cpcalc + 1

  IF( box_list(box_out)%volume .lt. box_list(box_in)%volume ) forward = .TRUE.

  ! pick a molecule at random to delete

  im = INT(rranf() * nmols(this_species,box_out)) + 1

  CALL Get_Index_Molecule(box_out,this_species,im,alive)
  
  CALL Save_Cartesian_Coordinates_Box(box_out,tot_mol_box_out)
  CALL Save_Cartesian_Coordinates_Box(box_in,tot_mol_box_in)

  
  l_charge_in = .FALSE.
  l_charge_out = .FALSE.
  IF((int_charge_sum_style(box_in) == charge_ewald .or. int_charge_sum_style(box_in) == charge_gaussian) .AND. has_charge(this_species)) l_charge_in = .TRUE.
  IF((int_charge_sum_style(box_out) == charge_ewald .or. int_charge_sum_style(box_out) == charge_gaussian) .AND. has_charge(this_species)) l_charge_out = .TRUE.
  
  ! Let us first insert the molecule in other box and determine
  ! if there is an overlap. Computationally this is more efficient
  ! in the event that an overlap is detected, the move can be
  ! immediately rejected. Note that this move chanegs the box identity
  ! and coordinate related quantites of 'alive. When we calcualte
  ! the energy of removing molecule from box_out we need to restore
  !

  molecule_list(alive,this_species)%which_box = box_in

  ! set deletion flag to false and call the CBMC growth routine

  cbmc_overlap = .FALSE.

  delta_e_in = 0.0_DP

  IF (species_list(this_species)%fragment .AND. species_list(this_species)%int_insert .NE. int_igas ) THEN

     del_Flag = .FALSE.
     get_fragorder = .TRUE.
 
     ALLOCATE(frag_order(nfragments(this_species)))
 !   frag_order = 0
 !   print *, 'frag_order',frag_order 
     lambda_for_build = molecule_list(alive,this_species)%cfc_lambda
    !print *, species_list(this_species)%lcom
   
    IF (species_list(this_species)%lcom) THEN
        CALL Build_Molecule(alive,this_species,box_in,frag_order,lambda_for_build, P_seq,P_forward, &
          nrg_ring_frag_in,cbmc_overlap)
     ELSE
        CALL Build_Rigid_Fragment(alive,this_species,box_in,frag_order,lambda_for_build, P_seq,P_forward, &
          nrg_ring_frag_in,cbmc_overlap)
     END IF
  ELSE
     CALL New_Positions(box_in,alive,this_species,rand_igas)
  END IF
  
 ! print *, frag_order 
  CALL Get_COM(alive,this_species)
  CALL Compute_Max_COM_Distance(alive,this_species)

  CALL Fold_Molecule(alive,this_species,box_in)

  ! Note that the call to Build_Molecule returns cbmc_overlap flag to be
  ! true not only in the case of an overlap but also when the weight
  ! of a given fragment is zero.

  IF (.NOT. cbmc_overlap) THEN
     
        CALL compute_molecule_inter_overlap(alive,this_species,box_in, inter_overlap)

  END IF

  IF (inter_overlap .OR. cbmc_overlap) THEN
     ! change the box identity
     molecule_list(alive,this_species)%which_box = box_out
    ! Revert the coordinates
     
     CALL Revert_Old_Cartesian_Coordinates(alive,this_species)
    
     ! It may so happen that only a few atoms have been assigned the exist flag to be true
     ! during molecule building so make all the atom beads of alive as true. Also
     ! reset the cfc_lambda

     atom_list(:,alive,this_species)%exist = .TRUE.
     molecule_list(alive,this_species)%cfc_lambda = 1.0_DP

     RETURN

  END IF

    IF (int_vdw_sum_style(box_in) == vdw_cut_tail .or. int_vdw_sum_style(box_in) == born_cut_tail) THEN
     nbeads_in(:) = nint_beads(:,box_in)

     DO i = 1, natoms(this_species)
        i_type = nonbond_list(i,this_species)%atom_type_number
        nint_beads(i_type,box_in) = nint_beads(i_type,box_in) + 1
     END DO      
     
  END IF
   
   nmols(this_species,box_in) = nmols(this_species,box_in) + 1   
   energy_old_in = energy(box_in)
   
   conv_drude1 = .TRUE.

   
   IF (shell_mpm)   CALL shell_relax(box_in, conv_drude1)
   IF (conv_drude1) CALL Compute_Total_System_Energy(box_in, .TRUE., overlap)

   nmols(this_species,box_in) = nmols(this_species,box_in) - 1

  ! Now compute the energy of the molecule assuming that it has been
  ! deleted. Now if we are at this stage then we need to preserve
  ! coordinates of the molecule inserted in other box.

  ALLOCATE(new_atom_list(natoms(this_species)))

  new_atom_list(:) = atom_list(:,alive,this_species)
  new_molecule_list = molecule_list(alive,this_species)


  ! Now we will swap the particle in another box and compute the energy
  ! For this, we will change the box identity for this molecule
  ! Now revert to the original coordinates

  CALL Revert_Old_Cartesian_Coordinates(alive,this_species)  
  molecule_list(alive,this_species)%which_box = box_out
  

  delta_e_out = 0.0_DP

 
  IF ( species_list(this_species)%fragment .AND. species_list(this_species)%int_insert .NE. int_igas) THEN
     ! Note that we need to use the same fragment order in which we inserted the molecule
     ! above. so in this case frag_order becomes input to the routine. We obtain
     ! P_delete via this call. Note that, cbmc_overlap should be false as we are
     ! dealing with an existing molecule.
     del_Flag = .TRUE. 
     get_fragorder = .FALSE.

     IF (species_list(this_species)%lcom) THEN
        CALL Build_Molecule(alive,this_species,box_out,frag_order,lambda_for_build, P_seq,P_reverse, &
             nrg_ring_frag_out,cbmc_overlap)
     ELSE
        CALL Build_Rigid_Fragment(alive,this_species,box_out,frag_order,lambda_for_build, &
             P_seq, P_reverse, nrg_ring_frag_out,cbmc_overlap)
     END IF

!     DEALLOCATE(frag_order)
 
     IF (cbmc_overlap) THEN
        ! energy overlap was found in the old configuration

        atom_list(1:natoms(this_species),alive,this_species)%exist = .TRUE.
        CALL Revert_Old_Cartesian_Coordinates(alive,this_species)

        molecule_list(alive,this_species)%cfc_lambda = 1.0_DP
        
     END IF


  END IF

  molecule_list(alive,this_species)%which_box = box_in
 
  IF (int_vdw_sum_style(box_out) == vdw_cut_tail .or. int_vdw_sum_style(box_out) == born_cut_tail) THEN

     nbeads_out(:) = nint_beads(:,box_out)

     DO i = 1, natoms(this_species)
        i_type = nonbond_list(i,this_species)%atom_type_number
        nint_beads(i_type,box_out) = nint_beads(i_type,box_out) - 1
     END DO

  END IF

   nmols(this_species,box_out) = nmols(this_species,box_out) - 1
   energy_old_out = energy(box_out)
   conv_drude2 = .TRUE.
   
   IF(shell_mpm) CALL shell_relax(box_out, conv_drude2)
   IF (conv_drude2) CALL Compute_Total_System_Energy(box_out, .TRUE., overlap)
   
   nmols(this_species,box_out) = nmols(this_species,box_out) + 1

  ! Define a factor that will be used to accept or reject the move. Note that
  ! the change in energy of box_out is actually negative of delta_e_out calculated
  ! above

  !factor = beta(box_in)*delta_e_in_pacc - beta(box_out)*delta_e_out_pacc

  !factor = beta(box_in)*(energy(box_in)%total - energy_old_in%total - e_angle_in - nrg_ring_frag_in) &
  !         & + beta(box_out)*(energy(box_out)%total - energy_old_out%total - e_angle_out - nrg_ring_frag_out)


  factor = beta(box_in)*(energy(box_in)%total - energy_old_in%total )  + beta(box_out)*(energy(box_out)%total - energy_old_out%total)

  !IF(dabs(e_angle_in)> 1d-5 .or. dabs(nrg_ring_frag_in)>1d-5 .or. dabs(e_angle_out) > 1d-5 .or. dabs(nrg_ring_frag_out)>1d-5) THEN
  !  print *, 'nrg_ring_frag_in .NE. 0.0', nrg_ring_frag_in, e_angle_in,e_angle_out, nrg_ring_frag_out
  !END IF

  factor = factor - DLOG((box_list(box_in)%volume * REAL(nmols(this_species,box_out),DP))/ &
                          (box_list(box_out)%volume * REAL(nmols(this_species,box_in) + 1, DP)))
  
  factor = factor + DLOG(P_forward) - DLOG(P_reverse)

  accept = accept_or_reject(factor)

  IF (accept .and. conv_drude1 .and. conv_Drude2) THEN
     
     nmols(this_species,box_in) = nmols(this_species,box_in) + 1
     nmols(this_species,box_out) = nmols(this_species,box_out) - 1

     nsuccess(this_species,box_in)%insertion = nsuccess(this_species,box_in)%insertion + 1
     nsuccess(this_species,box_out)%deletion = nsuccess(this_species,box_out)%deletion + 1

     DO i = 1,natoms(this_species)
       atom_list(i,alive,this_species) = new_atom_list(i)
     END DO
     molecule_list(alive,this_species) = new_molecule_list

     CALL Fold_Molecule(alive,this_species,box_in)
     
     ! Restore the coordinates of the molecule due to successful insertion
 
     CALL Get_Internal_Coordinates(alive,this_species)
   

  ELSE

     molecule_list(alive,this_species)%which_box = box_out
     IF ( int_vdw_sum_style(box_in) == vdw_cut_tail .or. int_vdw_sum_style(box_in) == born_cut_tail ) THEN

        nint_beads(:,box_in) = nbeads_in(:)

     END IF

     IF ( int_vdw_sum_style(box_out) == vdw_cut_tail .or. int_vdw_sum_style(box_in) == born_cut_tail ) THEN

        nint_beads(:,box_out) = nbeads_out(:)

     END IF

       CALL Reset_Cartesian_Coordinates_Box(box_in)
       CALL Reset_Cartesian_Coordinates_Box(box_out)

       energy(box_out) = energy_old_out
       energy(box_in)  = energy_old_in

  END IF
  

  IF (ALLOCATED(sorbate_id)) DEALLOCATE(sorbate_id)
  IF (ALLOCATED(sorbate_x)) DEALLOCATE(sorbate_x)
  DEALLOCATE(new_atom_list)

END SUBROUTINE Shell_GEMC_Particle_Transfer


