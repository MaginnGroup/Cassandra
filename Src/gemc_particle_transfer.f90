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

SUBROUTINE GEMC_Particle_Transfer

  !*****************************************************************************
  !
  ! This subroutine performs particle swaps in a GEMC ensemble.
  !
  ! The algorithm to perform this move is the same as described in Rull et al,
  ! Mol. Phys. 85(3), 435-447. The algorithm used here is described on page 444.
  !
  ! Called by
  !
  !   gemc_driver
  !
  ! Revision history
  !
  !   12/10/13 : Beta Release
  !   10/10/15 : Updated to use locate(im,is,ibox)   
  !
  ! DESCRIPTION: This subrouting performs the following steps:
  !
  ! Step 1) Select a box 'box_out' from which to remove a particle
  ! Step 2) Select a species 'this_species' according to its overall mole 
  !         fraction
  ! Step 3) Select a molecule 'alive' in box_out with uniform probability
  ! Step 4) Choose a position, orientation and conformation for alive in box_in
  ! Step 5) Calculate the change in box_in's potential energy from inserting
  !         alive
  ! Step 6) Calculate the change in box_out's potential energy from deleting
  !         alive
  ! Step 7) Accept of reject the move
  !
  !*****************************************************************************

  USE Global_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE IO_Utilities
  USE Fragment_Growth
  USE Pair_Nrg_Routines

  IMPLICIT NONE

  INTEGER :: box_in, box_out, i, k, this_species, i_type

  INTEGER :: nmols_sorbate, nsorbate, is, ibox, im_in, im_out ,alive, which_anchor
  INTEGER :: rand_igas, locate_in

  INTEGER, ALLOCATABLE :: sorbate_id(:), frag_order(:)

  REAL(DP), ALLOCATABLE :: sorbate_x(:)
  REAL(DP), ALLOCATABLE :: dx(:), dy(:), dz(:)

  REAL(DP) :: pick_species, delta_e_out, delta_e_out_pacc, e_bond_out, e_angle_out, e_dihed_out
  REAL(DP) :: E_intra_vdw_out, E_intra_qq_out, E_periodic_qq
  REAL(DP) :: E_inter_vdw_out, E_inter_qq_out
  REAL(DP) :: E_intra_vdw_igas, E_intra_qq_igas
  REAL(DP) :: E_reciprocal_out, E_self_out, e_lrc_out, e_improper_out
  REAL(DP) :: E_intra_vdw_in, E_intra_qq_in
  REAL(DP) :: E_inter_vdw_in, E_inter_qq_in
  REAL(DP) :: E_reciprocal_in, E_self_in, e_lrc_in, E_bond_in, E_improper_in
  REAL(DP) :: E_angle_in, E_dihed_in, delta_e_in, delta_e_in_pacc, potw, CP_energy
  REAL(DP) :: P_seq, P_forward, P_reverse, ln_pacc
  REAL(DP) :: lambda_for_build
  LOGICAL :: inter_overlap, accept_or_reject, cbmc_overlap, forward
  LOGICAL :: intra_overlap
 ! ring biasing variables

  REAL(DP) :: nrg_ring_frag_in, nrg_ring_frag_out

  TYPE(atom_class), ALLOCATABLE :: new_atom_list(:)
  TYPE(molecule_class) :: new_molecule_list

 ! Variables added for l_pair_nrg and reciprocal k vector storage

  INTEGER :: position

  REAL(DP), ALLOCATABLE :: cos_mol_old(:), sin_mol_old(:), cos_mol_new(:), sin_mol_new(:)
  REAL(DP) :: time0, time1, pick_box
  
  LOGICAL :: l_charge_in, l_charge_out

  potw = 1.0_DP
  forward = .FALSE.
  inter_overlap = .false.
  cbmc_overlap = .false.
  accept = .false.

  P_seq = 1.0_DP
  P_forward = 1.0_DP
  P_reverse = 1.0_DP

  !*****************************************************************************
  ! Step 1) Select a box 'box_out' from which to remove a particle
  !*****************************************************************************
  box_out = INT (rranf() * nbr_boxes) + 1

  ! choose a box in which the particle will be placed
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

  IF( box_list(box_out)%volume .lt. box_list(box_in)%volume ) forward = .TRUE.

  !*****************************************************************************
  ! Step 2) Select a species 'this_species' according to its overall mole 
  !         fraction
  !*****************************************************************************
  IF (l_mol_frac_swap) THEN
     
     ! Choose a species based on the overall mole fraction.
     ! We will base this on the sorbate mole fraction
     ALLOCATE(sorbate_id(nspecies), sorbate_x(nspecies))
     
     nmols_sorbate = 0
     nsorbate = 0
     sorbate_x(:) = 0.0_DP
     
     DO is = 1, nspecies
        IF (species_list(is)%species_type == 'SORBATE') THEN
           nsorbate = nsorbate + 1
           sorbate_id(nsorbate) = is
           nmols_sorbate = nmols_sorbate + SUM(nmols(is,1:nbr_boxes))
           sorbate_x(nsorbate) = REAL(SUM(nmols(is,1:nbr_boxes)),DP)
           ! do not sum over box=0, the placeholder for unused locates
        END IF
     END DO
  
     IF (nmols_sorbate == 0) RETURN

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

  IF (nmols(this_species,box_out) == 0) THEN
     IF (cpcollect)  CALL Chempot(box_in)
     RETURN
  END IF
  
  ! Increment counters
  ntrials(this_species,box_out)%deletion = &
     ntrials(this_species,box_out)%deletion + 1
  ntrials(this_species,box_in)%insertion = &
     ntrials(this_species,box_in)%insertion + 1
  ntrials(this_species,box_in)%cpcalc = &
     ntrials(this_species,box_in)%cpcalc + 1

  !*****************************************************************************
  ! Step 3) Select a molecule 'alive' in box_out with uniform probability
  !*****************************************************************************
  ! pick a molecule INDEX at random to delete
  im_out = INT(rranf() * nmols(this_species,box_out)) + 1
  
  ! Obtain the LOCATE of this molecule
  alive = locate(im_out,this_species,box_out)

  !*****************************************************************************
  ! Step 4) Choose a position, orientation and conformation for alive in box_in
  !*****************************************************************************
  ! Insert alive into box_in before deleting alive from box_out.
  ! If there is an overlap in box_in, the move can be
  ! immediately rejected. Note that this move changes the box identity
  ! and coordinate related quantites of 'alive'. When we calculate
  ! the energy of removing molecule from box_out we need to restore these

  ! Save the coordinates of 'alive' in 'box_out'
  CALL Save_Old_Cartesian_Coordinates(alive,this_species)
  CALL Compute_Molecule_Dihedral_Energy(alive,this_species,e_dihed_out)

  ! Save the interaction energies
  IF (l_pair_nrg) CALL Store_Molecule_Pair_Interaction_Arrays(alive,this_species, &
       box_out, E_inter_vdw_out, E_inter_qq_out)
  
  ! Save the k-vectors
  
  IF (int_charge_sum_style(box_in)  == charge_ewald .AND.&
      has_charge(this_species)) THEN
     ALLOCATE(cos_mol_old(nvecs(box_out)), sin_mol_old(nvecs(box_out)))
     CALL Get_Position_Alive(alive,this_species,position)
     
     !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
     cos_mol_old(:) = cos_mol(1:nvecs(box_out),position)
     sin_mol_old(:) = sin_mol(1:nvecs(box_out),position)
     !$OMP END PARALLEL WORKSHARE
  END IF

  ! Switch the box identity of alive
  molecule_list(alive,this_species)%which_box = box_in
  nmols(this_species,box_in) = nmols(this_species,box_in) + 1
  im_in = nmols(this_species,box_in) ! INDEX of alive in box_in
  locate(im_in,this_species,box_in) = alive ! link INDEX to LOCATE

  ! set deletion flag to false and call the CBMC growth routine
  cbmc_overlap = .FALSE.

  delta_e_in = 0.0_DP

  IF (species_list(this_species)%fragment .AND. &
      species_list(this_species)%int_insert .NE. int_igas ) THEN

     del_flag = .FALSE.
     get_fragorder = .TRUE.
     ALLOCATE(frag_order(nfragments(this_species)))
     lambda_for_build = molecule_list(alive,this_species)%frac
    IF (species_list(this_species)%lcom) THEN
        CALL Build_Molecule(alive,this_species,box_in,frag_order, &
                lambda_for_build,P_seq,P_forward,nrg_ring_frag_in, &
                cbmc_overlap)
     ELSE
        CALL Build_Rigid_Fragment(alive,this_species,box_in,frag_order, &
                lambda_for_build,P_seq,P_forward,nrg_ring_frag_in, &
                cbmc_overlap)
     END IF
  ELSE
     CALL New_Positions(box_in,alive,this_species,rand_igas)
  END IF

  CALL Get_COM(alive,this_species)
  CALL Compute_Max_COM_Distance(alive,this_species)
  CALL Fold_Molecule(alive,this_species,box_in)

  ! Build_Molecule returns cbmc_overlap == .TRUE. both in case of a core
  ! overlap and also when the weight of all trials is zero.

  IF (.NOT. cbmc_overlap) THEN
        CALL Compute_Molecule_Nonbond_Inter_Energy(alive,this_species, &
             E_inter_vdw_in,E_inter_qq_in,inter_overlap)
  END IF

  IF (inter_overlap .OR. cbmc_overlap) THEN
     ! reject the swap

     ! restore the box identity
     molecule_list(alive,this_species)%which_box = box_out
     locate(im_in,this_species,box_in) = 0 ! reset LOCATE
     nmols(this_species,box_in) = nmols(this_species,box_in) - 1
     ! restore the box_out coordinates
     CALL Revert_Old_Cartesian_Coordinates(alive,this_species)

     ! All atoms will not exist if inter_overlap was tripped before the 
     ! last fragment was placed in Build_Molecule. 
     ! Set exist to TRUE for all atoms and reset frac to 1
     atom_list(:,alive,this_species)%exist = .TRUE.
     molecule_list(alive,this_species)%frac = 1.0_DP

     IF (l_pair_nrg) THEN
        CALL Reset_Molecule_Pair_Interaction_Arrays(alive,this_species,box_out)
     END IF

     IF(ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
     IF(ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)

     RETURN

  END IF

  delta_e_in = delta_e_in + E_inter_vdw_in + E_inter_qq_in 
  
  !*****************************************************************************
  ! Step 5) Calculate the change in box_in's potential energy from inserting
  !         alive
  !*****************************************************************************
  ! If here then no overlap was detected. Calculate the rest of the energies
  IF(species_list(this_species)%int_insert == int_random) THEN

     CALL Compute_Molecule_Bond_Energy(alive,this_species,E_bond_in)
     CALL Compute_Molecule_Angle_Energy(alive,this_species,E_angle_in)
     CALL Compute_Molecule_Dihedral_Energy(alive,this_species,E_dihed_in)
     CALL Compute_Molecule_Improper_Energy(alive,this_species,E_improper_in)
  
  ELSE IF(species_list(this_species)%int_insert == int_igas) THEN

     E_bond_in = energy_igas(rand_igas,this_species)%bond
     E_angle_in = energy_igas(rand_igas,this_species)%angle
     E_dihed_in = energy_igas(rand_igas,this_species)%dihedral
     E_improper_in = energy_igas(rand_igas,this_species)%improper
  
  END IF

  delta_e_in = delta_e_in + E_bond_in + E_angle_in + E_dihed_in

  CALL Compute_Molecule_Nonbond_Intra_Energy(alive,this_species, &
       E_intra_vdw_in,E_intra_qq_in,E_periodic_qq,intra_overlap)
  E_inter_qq_in = E_inter_qq_in + E_periodic_qq

  ! Already added E_inter to delta_e, so add E_periodic directly
  delta_e_in = delta_e_in + E_intra_vdw_in + E_intra_qq_in + E_periodic_qq

  call cpu_time(time0)

  IF (int_charge_style(box_in) == charge_coul .AND. has_charge(this_species)) THEN
     
     IF (int_charge_sum_style(box_in) == charge_ewald) THEN

          ! Note that this call will change cos_mol, sin_mol of alive and this
          ! will have to be restored below while computing the energy of box_out
          ! without molecule alive. 
          CALL Update_System_Ewald_Reciprocal_Energy(alive,this_species,box_in, &
               int_insertion,E_reciprocal_in)
     
          delta_e_in = delta_e_in + (E_reciprocal_in - energy(box_in)%ewald_reciprocal)
     END IF
     
     CALL Compute_Molecule_Self_Energy(alive,this_species,box_in, &
          E_self_in)
     delta_e_in = delta_e_in + E_self_in 

  END IF

  call cpu_time(time1)
  copy_time = copy_time + time1-time0

  IF (int_vdw_sum_style(box_in) == vdw_cut_tail .AND. &
		int_vdw_style(box_in) == vdw_lj) THEN
     nbeads_in(:) = nint_beads(:,box_in)

     DO i = 1, natoms(this_species)
        i_type = nonbond_list(i,this_species)%atom_type_number
        nint_beads(i_type,box_in) = nint_beads(i_type,box_in) + 1
     END DO
        
     CALL Compute_LR_Correction(box_in,e_lrc_in)
     delta_e_in = delta_e_in + e_lrc_in - energy(box_in)%lrc

  ELSEIF (int_vdw_sum_style(box_in) == vdw_cut_tail .AND. &
		int_vdw_style(box_in) == vdw_mie) THEN 
        nbeads_in(:) = nint_beads_mie(this_species,:,box_in)

     DO i = 1, natoms(this_species)
        i_type = nonbond_list(i,this_species)%atom_type_number
        nint_beads_mie(this_species,i_type,box_in) = nint_beads_mie(this_species,i_type,box_in) + 1
     END DO
        
     CALL Compute_LR_Correction(box_in,e_lrc_in)
     delta_e_in = delta_e_in + e_lrc_in - energy(box_in)%lrc
   
  END IF

  IF(cpcollect) THEN
  
     potw = 1.0_DP
     CP_energy = delta_e_in
     
     IF(species_list(this_species)%fragment) THEN
        potw = 1.0_DP / (P_forward * kappa_ins*kappa_rot*kappa_dih &
             ** (nfragments(this_species)-1)) 
        CP_energy = delta_e_in - E_angle_in 
     END IF
   
     chpot(this_species,box_in) = chpot(this_species,box_in) &
      + potw * (box_list(box_in)%volume &
      / (REAL(nmols(this_species,box_in)))) * DEXP(-beta(box_in) * CP_energy)

  END IF

  !*****************************************************************************
  ! Step 6) Calculate the change in box_out's potential energy from deleting
  !         alive
  !*****************************************************************************
  ! Need to preserve coordinates of 'alive' in box_in.
  ALLOCATE(new_atom_list(natoms(this_species)))
  new_atom_list(:) = atom_list(:,alive,this_species)
  new_molecule_list = molecule_list(alive,this_species)

  ! Change the box identity of alive back to box_out
  molecule_list(alive,this_species)%which_box = box_out
  ! Restore the box_out coordinates
  CALL Revert_Old_Cartesian_Coordinates(alive,this_species)

  delta_e_out = 0.0_DP

  ! Obtain the weight of the chain if it is made up of fragments
  call cpu_time(time0)
  IF ( species_list(this_species)%fragment .AND. &
       species_list(this_species)%int_insert .NE. int_igas) THEN
     ! The fragment order was decided when inserting alive into box_in
     ! Use the same fragment order to calculate trial insertions into box_out 
     ! 
     ! So frag_order and P_seq are inputs to the Build_Molecule routine
     ! We obtain P_reverse via this call. Note that, cbmc_overlap 
     ! must be false as we are dealing with an existing molecule.
     del_flag = .TRUE. 
     get_fragorder = .FALSE.

     IF (species_list(this_species)%lcom) THEN
        CALL Build_Molecule(alive,this_species,box_out,frag_order, &
                lambda_for_build,P_seq,P_reverse,nrg_ring_frag_out, &
                cbmc_overlap)
     ELSE
        CALL Build_Rigid_Fragment(alive,this_species,box_out,frag_order, &
                lambda_for_build,P_seq,P_reverse,nrg_ring_frag_out, &
                cbmc_overlap)
     END IF
        
     IF (cbmc_overlap) THEN
        ! If this flag gets tripped, there is an error in the code
        err_msg = ""
        err_msg(1) = "Error: existing configuration of " // & 
           "molecule " // TRIM(Int_To_String(alive))
        err_msg(2) = "of species " // TRIM(Int_To_String(this_species)) // &
           " in box " // TRIM(Int_To_String(box_out)) // &
           " tripped an overlap error"
        CALL Clean_Abort(err_msg,'GEMC_Particle_Transfer')

!        atom_list(1:natoms(this_species),alive,this_species)%exist = .TRUE.
!        CALL Revert_Old_Cartesian_Coordinates(alive,this_species)
!
!        molecule_list(alive,this_species)%frac = 1.0_DP
     END IF

  END IF

  CALL Get_COM(alive,this_species)
  CALL Compute_Max_COM_Distance(alive,this_species)
  ! debug to see if the COM and max_com_distance are identical

  ! bonded energies
  CALL Compute_Molecule_Bond_Energy(alive,this_species,e_bond_out)
  CALL Compute_Molecule_Angle_Energy(alive,this_species,e_angle_out)
  CALL Compute_Molecule_Dihedral_Energy(alive,this_species,e_dihed_out)
  CALL Compute_Molecule_Improper_Energy(alive,this_species,e_improper_out)

  delta_e_out = delta_e_out - e_bond_out - e_angle_out - e_dihed_out &
              - e_improper_out

  ! Nonbonded energy  
  CALL Compute_Molecule_Nonbond_Intra_Energy(alive,this_species, &
       E_intra_vdw_out,E_intra_qq_out,E_periodic_qq,intra_overlap)

  IF ( .NOT. l_pair_nrg) THEN
     CALL Compute_Molecule_Nonbond_Inter_Energy(alive,this_species, &
          E_inter_vdw_out,E_inter_qq_out,inter_overlap)
  END IF
  E_inter_qq_out = E_inter_qq_out + E_periodic_qq

  delta_e_out = delta_e_out - E_intra_vdw_out - E_intra_qq_out &
              - E_inter_vdw_out - E_inter_qq_out




  IF (int_charge_style(box_out) == charge_coul .AND. has_charge(this_species)) THEN
        IF (int_charge_sum_style(box_in) == charge_ewald .AND. &
            int_charge_sum_style(box_out) == charge_ewald) THEN
           ! Restore the cos_mol and sin_mol as they changed above
           ! but restoring will destroy the newly computed vector so now here allocate
           ! cos_mol_new
           ! sin_mol_new vectors so that if the move is accepted we can restore these
      
           call cpu_time(time0)
      
           ALLOCATE(cos_mol_new(nvecs(box_in)))
           ALLOCATE(sin_mol_new(nvecs(box_in)))
      
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
           cos_mol_new(:) = cos_mol(1:nvecs(box_in),position)
           sin_mol_new(:) = sin_mol(1:nvecs(box_in),position)


           cos_mol(1:nvecs(box_out),position) = cos_mol_old(1:nvecs(box_out))
           sin_mol(1:nvecs(box_out),position) = sin_mol_old(1:nvecs(box_out))
           !$OMP END PARALLEL WORKSHARE
      
           call cpu_time(time1)
      
      !     copy_time = copy_time + time1-time0
      
           CALL Update_System_Ewald_Reciprocal_Energy(alive,this_species, &
                box_out,int_deletion,E_reciprocal_out)
      
           delta_e_out = delta_e_out + (E_reciprocal_out - energy(box_out)%ewald_reciprocal)
      
        END IF
      
        CALL Compute_Molecule_Self_Energy(alive,this_species,box_out, &
          E_self_out)

        delta_e_out = delta_e_out - E_self_out

  END IF

  
  IF (int_vdw_sum_style(box_out) == vdw_cut_tail .AND. &
		int_vdw_style(box_out) == vdw_lj) THEN
     nbeads_out(:) = nint_beads(:,box_out)
     DO i = 1, natoms(this_species)
        i_type = nonbond_list(i,this_species)%atom_type_number
        nint_beads(i_type,box_out) = nint_beads(i_type,box_out) - 1
     END DO

     CALL Compute_LR_correction(box_out,e_lrc_out)
     delta_e_out = delta_e_out + ( e_lrc_out - energy(box_out)%lrc )

  ELSEIF (int_vdw_sum_style(box_out) == vdw_cut_tail .AND. &
		int_vdw_style(box_out) == vdw_mie) THEN
     nbeads_out(:) = nint_beads_mie(this_species,:,box_out)

     DO i = 1, natoms(this_species)
        i_type = nonbond_list(i,this_species)%atom_type_number
        nint_beads_mie(this_species,i_type,box_out) = nint_beads_mie(this_species,i_type,box_out) - 1
     END DO

     CALL Compute_LR_correction(box_out,e_lrc_out)
     delta_e_out = delta_e_out + ( e_lrc_out - energy(box_out)%lrc )
  END IF

  delta_e_in_pacc = delta_e_in
  delta_e_out_pacc = delta_e_out

  IF(species_list(this_species)%int_insert == int_igas) THEN
     igas_flag = .TRUE.
     CALL Compute_Molecule_Nonbond_Intra_Energy(alive,this_species, &
          E_intra_vdw_igas,E_intra_qq_igas,E_periodic_qq,intra_overlap)
     ! Ideal gas should not interact with it's periodic image, so E_periodic_qq
     ! is not added the delta_e
     igas_flag = .FALSE. 
     delta_e_out_pacc = delta_e_out_pacc + e_bond_out + e_angle_out &
                      + e_dihed_out + e_improper_out &
                      + E_intra_vdw_igas + E_intra_qq_igas
     delta_e_in_pacc = delta_e_in_pacc &
                     - energy_igas(rand_igas,this_species)%total
  END IF

  IF(species_list(this_species)%fragment .AND. &
     species_list(this_species)%int_insert .NE. int_igas) THEN
     delta_e_in_pacc  = delta_e_in_pacc  - e_angle_in  - nrg_ring_frag_in
     delta_e_out_pacc = delta_e_out_pacc + e_angle_out + nrg_ring_frag_out
  END IF

  !*****************************************************************************
  ! Step 7) Accept of reject the move
  !*****************************************************************************
  ! Define ln_pacc that will be used to accept or reject the move. Note that
  ! the change in energy of box_out is actually negative of delta_e_out 
  ! calculated above

  ln_pacc = beta(box_in)*delta_e_in_pacc + beta(box_out)*delta_e_out_pacc

  ln_pacc = ln_pacc - DLOG(box_list(box_in)%volume) &
                    + DLOG(box_list(box_out)%volume) &
                    - DLOG(REAL(nmols(this_species,box_out),DP)) &
                    + DLOG(REAL(nmols(this_species,box_in), DP))

  ! The same order of insertion is used in both the insertion and 
  ! reverse deletion, so P_seq does not factor into ln_pacc
  ln_pacc = ln_pacc + DLOG(P_forward / P_reverse)

  accept = accept_or_reject(ln_pacc)

  IF (accept) THEN
     ! accept the swap
     
     ! already updated the number of molecules in box_in

     ! remove the deleted molecule from box_out locate
     IF (im_out < nmols(this_species,box_out)) THEN
        DO k = im_out + 1, nmols(this_species,box_out)
           locate(k-1,this_species,box_out) = locate(k,this_species,box_out)
        END DO
     END IF
     locate(nmols(this_species,box_out),this_species,box_out) = 0

     ! Update the number of molecules in box_out
     nmols(this_species,box_out) = nmols(this_species,box_out) - 1

     ! Set the coordinates and properties of molecule alive to box_in values
     DO i = 1,natoms(this_species)
       atom_list(i,alive,this_species) = new_atom_list(i)
     END DO
     molecule_list(alive,this_species) = new_molecule_list
     CALL Fold_Molecule(alive,this_species,box_in)

     IF (int_charge_sum_style(box_in) == charge_ewald .AND. &
         has_charge(this_species)) THEN
        call cpu_time(time0)
        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        cos_mol(1:nvecs(box_in),position) = cos_mol_new(:)
        sin_mol(1:nvecs(box_in),position) = sin_mol_new(:)
        !$OMP END PARALLEL WORKSHARE
        
        DEALLOCATE(cos_mol_new,sin_mol_new)
        
        call cpu_time(time1)
!        copy_time = copy_time + time1-time0
     END IF

     IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
     IF (ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
     IF (ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)
     IF (ALLOCATED(cos_mol_new)) DEALLOCATE(cos_mol_new)
     IF (ALLOCATED(sin_mol_new)) DEALLOCATE(sin_mol_new)

     ! Restore the coordinates of the molecule due to successful insertion
     CALL Get_Internal_Coordinates(alive,this_species)

     ! Update energies for each box
     ! box_in
     energy(box_in)%total = energy(box_in)%total + delta_e_in
     energy(box_in)%intra = energy(box_in)%intra + e_bond_in + e_angle_in &
                          + e_dihed_in
     energy(box_in)%bond = energy(box_in)%bond + e_bond_in
     energy(box_in)%angle = energy(box_in)%angle + e_angle_in
     energy(box_in)%dihedral = energy(box_in)%dihedral + e_dihed_in
     energy(box_in)%intra_vdw = energy(box_in)%intra_vdw + e_intra_vdw_in
     energy(box_in)%intra_q = energy(box_in)%intra_q + e_intra_qq_in
     energy(box_in)%inter_vdw = energy(box_in)%inter_vdw + e_inter_vdw_in
     energy(box_in)%inter_q = energy(box_in)%inter_q + e_inter_qq_in

     IF (int_vdw_sum_style(box_in) == vdw_cut_tail) THEN
        energy(box_in)%lrc = e_lrc_in
     END IF

     ! for box_out
     energy(box_out)%total = energy(box_out)%total + delta_e_out
     energy(box_out)%intra = energy(box_out)%intra - e_bond_out - e_angle_out &
                           - e_dihed_out
     energy(box_out)%bond = energy(box_out)%bond - e_bond_out
     energy(box_out)%angle = energy(box_out)%angle - e_angle_out
     energy(box_out)%dihedral = energy(box_out)%dihedral - e_dihed_out
     energy(box_out)%intra_vdw = energy(box_out)%intra_vdw - e_intra_vdw_out
     energy(box_out)%intra_q   = energy(box_out)%intra_q - e_intra_qq_out
     energy(box_out)%inter_vdw = energy(box_out)%inter_vdw - e_inter_vdw_out
     energy(box_out)%inter_q   = energy(box_out)%inter_q - e_inter_qq_out

     IF (int_vdw_sum_style(box_out) == vdw_cut_tail) THEN
        energy(box_out)%lrc = e_lrc_out
     END IF

     IF (has_charge(this_species)) THEN
        IF (int_charge_sum_style(box_in) == charge_ewald) energy(box_in)%ewald_reciprocal = E_reciprocal_in        
        IF (int_charge_sum_style(box_out) == charge_ewald) energy(box_out)%ewald_reciprocal = E_reciprocal_out
        IF (int_charge_style(box_in) == charge_coul) energy(box_in)%self = energy(box_in)%self + E_self_in
        IF (int_charge_style(box_out) == charge_coul) energy(box_out)%self = energy(box_out)%self + E_self_out
     END IF


     ! Increment counter
     nsuccess(this_species,box_in)%insertion = &
        nsuccess(this_species,box_in)%insertion + 1
     nsuccess(this_species,box_out)%deletion = &
        nsuccess(this_species,box_out)%deletion + 1

  ELSE
     ! reject the swap. 

     ! Atomic coordinates have not changed as we used the original coordinate
     ! in box_out to calculate removal energies

     ! Reset the number of molecules and locate of box_in
     locate(im_in,this_species,box_in) = 0
     nmols(this_species,box_in) = nmols(this_species,box_in) - 1


     IF (has_charge(this_species)) THEN
         ! Restore the reciprocal space k vectors
         IF (int_charge_sum_style(box_in) == charge_ewald) THEN
            !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
            cos_sum(1:nvecs(box_in),box_in) = cos_sum_old(1:nvecs(box_in),box_in)
            sin_sum(1:nvecs(box_in),box_in) = sin_sum_old(1:nvecs(box_in),box_in)
            !$OMP END PARALLEL WORKSHARE
    
            DEALLOCATE(cos_mol_new,sin_mol_new)
         END IF
    
         IF (int_charge_sum_style(box_out) == charge_ewald) THEN
            !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
            cos_sum(1:nvecs(box_out),box_out) = cos_sum_old(1:nvecs(box_out),box_out)
            sin_sum(1:nvecs(box_out),box_out) = sin_sum_old(1:nvecs(box_out),box_out)
            
            cos_mol(1:nvecs(box_out),position) = cos_mol_old(:)
            sin_mol(1:nvecs(box_out),position) = sin_mol_old(:)
            !$OMP END PARALLEL WORKSHARE
    
            DEALLOCATE(cos_mol_old)
            DEALLOCATE(sin_mol_old)
         END IF
     END IF

     IF (l_pair_nrg) THEN
        CALL Reset_Molecule_Pair_Interaction_Arrays(alive,this_species,box_out)
     END IF

     IF ( int_vdw_sum_style(box_in) == vdw_cut_tail ) THEN
        IF (int_vdw_style(box_in) == vdw_lj) THEN
           nint_beads(:,box_in) = nbeads_in(:)
        ELSEIF (int_vdw_style(box_in) == vdw_mie) THEN
           nint_beads_mie(this_species,:,box_in) = nbeads_in(:)
        END IF
     END IF

  END IF
  
  IF (ALLOCATED(sorbate_id)) DEALLOCATE(sorbate_id)
  IF (ALLOCATED(sorbate_x)) DEALLOCATE(sorbate_x)
  DEALLOCATE(new_atom_list)

  IF (verbose_log) THEN
    WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8)') i_mcstep, 'swap_to' , alive, this_species, box_in, accept
  END IF

END SUBROUTINE GEMC_Particle_Transfer

SUBROUTINE New_Positions(this_box,alive,is,rand_igas)

  USE Global_Variables
  USE Rotation_Routines
  USE File_Names

  IMPLICIT NONE
  
  INTEGER :: this_box,alive,is, rand_igas, i
  REAL(DP) :: dx,dy,dz, this_x, this_y, this_z


  atom_list(:,alive,is)%exist = .TRUE.

  IF ( species_list(is)%int_insert == int_random ) THEN
     
     ! COM of the species from the initial configuration is 
     
     molecule_list(alive,is)%xcom = species_list(is)%xcom
     molecule_list(alive,is)%ycom = species_list(is)%ycom
     molecule_list(alive,is)%zcom = species_list(is)%zcom     
     
     atom_list(:,alive,is)%rxp = init_list(:,1,is)%rxp
     atom_list(:,alive,is)%ryp = init_list(:,1,is)%ryp
     atom_list(:,alive,is)%rzp = init_list(:,1,is)%rzp
     

     IF (natoms(is) > 1) CALL Rotate_Molecule_Eulerian(alive,is)
     
     IF ( box_list(this_box)%int_box_shape == int_cubic ) THEN
        
        molecule_list(alive,is)%xcom = (rranf() - 0.5_DP) * box_list(this_box)%length(1,1)
        molecule_list(alive,is)%ycom = (rranf() - 0.5_DP) * box_list(this_box)%length(2,2)
        molecule_list(alive,is)%zcom = (rranf() - 0.5_DP) * box_list(this_box)%length(3,3)
     END IF

     ! Coordinates obtained are for the initial coordinates so translate it to the current position
     
     dx = molecule_list(alive,is)%xcom - species_list(is)%xcom
     dy = molecule_list(alive,is)%ycom - species_list(is)%ycom
     dz = molecule_list(alive,is)%zcom - species_list(is)%zcom
     
     atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp + dx
     atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp + dy
     atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp + dz


  ELSE IF (species_list(is)%int_insert == int_igas) THEN
     
     ! obtain a ranom configuratio from the file

     rand_igas = (rranf() * n_igas(is)) + 1

     molecule_list(alive,is)%xcom = molecule_list_igas(rand_igas,is)%xcom
     molecule_list(alive,is)%ycom = molecule_list_igas(rand_igas,is)%ycom
     molecule_list(alive,is)%zcom = molecule_list_igas(rand_igas,is)%zcom

     atom_list(:,alive,is)%rxp = atom_list_igas(:,rand_igas,is)%rxp
     atom_list(:,alive,is)%ryp = atom_list_igas(:,rand_igas,is)%ryp
     atom_list(:,alive,is)%rzp = atom_list_igas(:,rand_igas,is)%rzp

     CALL Get_COM(alive,is)

     IF (natoms(is) > 1) CALL Rotate_Molecule_Eulerian(alive,is)

     ! Randomly place the molecule in the simulation cell

     IF (box_list(this_box)%int_box_shape == int_cubic) THEN
        
        this_x = (0.5_DP - rranf()) * box_list(this_box)%length(1,1)
        this_y = (0.5_DP - rranf()) * box_list(this_box)%length(2,2)
        this_z = (0.5_DP - rranf()) * box_list(this_box)%length(3,3)


     END IF

     ! Obtain atomic positions

     atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp - &
                                 molecule_list(alive,is)%xcom + this_x
     atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp - &
                                 molecule_list(alive,is)%ycom + this_y
     atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp - &
                                 molecule_list(alive,is)%zcom + this_z

     ! Now change the COM position correspondig to the random placement

     molecule_list(alive,is)%xcom = this_x
     molecule_list(alive,is)%ycom = this_y
     molecule_list(alive,is)%zcom = this_z

  END IF
             
END SUBROUTINE New_Positions


