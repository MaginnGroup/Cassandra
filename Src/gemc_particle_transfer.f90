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

SUBROUTINE GEMC_Particle_Transfer(box_in, box_out)

  !********************************************************************************
  !
  ! This subroutine performs particle swaps in a GEMC ensemble.
  !
  ! The algorithm to perform this move is the same as described in Rull et al,
  ! Mol. Phys. 85(3), 435-447. The algorithm used here is described on page 444
  ! of the paper and it relates to mixture.
  !
  ! Briefly
  !
  ! Box selection: A box is randomly selected for particle removal
  ! Species selection: Type of the molecule that is being removed is decided
  !                    according to its overall mole fraction
  ! Molecule selection : A molecule of the above species is randomly selected
  !                      from the chosen simulation box.
  ! Volume selection : this particle is randomly inserted in the other box
  ! Acceptance criterion is evaluated.
  !
  !
  ! Called by
  !
  !   gemc_driver
  !
  ! Revision history
  !
  !   12/10/13 : Beta Release
  !*********************************************************************************

  USE Run_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE IO_Utilities
  USE Fragment_Growth
  USE Pair_Nrg_Routines

  IMPLICIT NONE

  INTEGER :: box_in, box_out, i, this_species, i_type

  INTEGER :: nmols_sorbate, nsorbate, is, ibox, im ,alive, which_anchor
  INTEGER :: rand_igas 

  INTEGER, ALLOCATABLE :: sorbate_id(:), frag_order(:)

  REAL(DP), ALLOCATABLE :: sorbate_x(:)
  REAL(DP), ALLOCATABLE :: dx(:), dy(:), dz(:)

  REAL(DP) :: pick_species, delta_e_out, delta_e_out_pacc, e_bond_out, e_angle_out, e_dihed_out
  REAL(DP) :: E_intra_vdw_out, E_intra_qq_out
  REAL(DP) :: E_inter_vdw_out, E_inter_qq_out
  REAL(DP) :: E_intra_vdw_igas, E_intra_qq_igas
  REAL(DP) :: E_reciprocal_move_out, E_self_move_out, e_lrc_out, e_improper_out
  REAL(DP) :: E_intra_vdw_in, E_intra_qq_in
  REAL(DP) :: E_inter_vdw_in, E_inter_qq_in
  REAL(DP) :: E_reciprocal_move_in, E_self_move_in, e_lrc_in, E_bond_in, E_improper_in
  REAL(DP) :: E_angle_in, E_dihed_in, delta_e_in, delta_e_in_pacc, factor, P_forward, P_reverse, potw, CP_energy
  REAL(DP) :: lambda_for_build
  LOGICAL :: inter_overlap, accept, accept_or_reject, cbmc_overlap, forward
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

  P_forward = 1.0_DP
  P_reverse = 1.0_DP

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

!  write(*,*) p_forward, p_reverse

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
        IF (species_list(is)%species_type == 'SORBATE') THEN
           nsorbate = nsorbate + 1
           sorbate_id(nsorbate) = is
           nmols_sorbate = nmols_sorbate + SUM(nmols(is,:))
           sorbate_x(nsorbate) = REAL(SUM(nmols(is,:)),DP)
        END IF
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
     IF(cpcollect)  CALL Chempot(box_in)
     RETURN
  END IF
  
  ntrials(this_species,box_out)%deletion = ntrials(this_species,box_out)%deletion + 1
  ntrials(this_species,box_in)%insertion = ntrials(this_species,box_in)%insertion + 1
  ntrials(this_species,box_in)%cpcalc = ntrials(this_species,box_in)%cpcalc + 1

  IF( box_list(box_out)%volume .lt. box_list(box_in)%volume ) forward = .TRUE.

  ! pick a molecule at random to delete

  im = INT(rranf() * nmols(this_species,box_out)) + 1

  ! Obtain the index of this molecule

  CALL Get_Index_Molecule(box_out,this_species,im,alive)

  CALL Save_Old_Cartesian_Coordinates(alive,this_species)
  CALL Compute_Molecule_Dihedral_Energy(alive,this_species,e_dihed_out)

  
  IF (l_pair_nrg) CALL Store_Molecule_Pair_Interaction_Arrays(alive,this_species, &
       box_out, E_inter_vdw_out, E_inter_qq_out)
  
  
  l_charge_in = .FALSE.
  l_charge_out = .FALSE.
  IF(int_charge_sum_style(box_in) == charge_ewald .AND. has_charge(this_species)) l_charge_in = .TRUE.
  IF(int_charge_sum_style(box_out) == charge_ewald .AND. has_charge(this_species)) l_charge_out = .TRUE.
  
  IF (l_charge_out) THEN
     
     ALLOCATE(cos_mol_old(nvecs(box_out)), sin_mol_old(nvecs(box_out)))
     CALL Get_Position_Alive(alive,this_species,position)
     
     !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
     cos_mol_old(:) = cos_mol(1:nvecs(box_out),position)
     sin_mol_old(:) = sin_mol(1:nvecs(box_out),position)
     !$OMP END PARALLEL WORKSHARE
     
  END IF

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
     lambda_for_build = molecule_list(alive,this_species)%cfc_lambda
    IF (species_list(this_species)%lcom) THEN
        CALL Build_Molecule(alive,this_species,box_in,frag_order,lambda_for_build, which_anchor,P_forward, &
          nrg_ring_frag_in,cbmc_overlap)
     ELSE
        CALL Build_Rigid_Fragment(alive,this_species,box_in,frag_order,lambda_for_build, which_anchor,P_forward, &
          nrg_ring_frag_in,cbmc_overlap)
     END IF
  ELSE
     CALL New_Positions(box_in,alive,this_species,rand_igas)
  END IF

  CALL Get_COM(alive,this_species)
  CALL Compute_Max_COM_Distance(alive,this_species)
  CALL Fold_Molecule(alive,this_species,box_in)

  ! Note that the call to Build_Molecule returns cbmc_overlap flag to be
  ! true not only in the case of an overlap but also when the weight
  ! of a given fragment is zero.

  IF (.NOT. cbmc_overlap) THEN
     
        CALL Compute_Molecule_Nonbond_Inter_Energy(alive,this_species,E_inter_vdw_in, &
             E_inter_qq_in,inter_overlap)
  END IF

  IF (inter_overlap .OR. cbmc_overlap) THEN
     ! change the box identity
     molecule_list(alive,this_species)%which_box = box_out
!     mol_start(alive,this_species) = 0
     ! Revert the coordinates
     CALL Revert_Old_Cartesian_Coordinates(alive,this_species)

     ! It may so happen that only a few atoms have been assigned the exist flag to be true
     ! during molecule building so make all the atom beads of alive as true. Also
     ! reset the cfc_lambda

     atom_list(:,alive,this_species)%exist = .TRUE.
     molecule_list(alive,this_species)%cfc_lambda = 1.0_DP

     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(alive,this_species,box_out)

     IF(ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
     IF(ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)

     RETURN

  END IF


  delta_e_in = delta_e_in + E_inter_vdw_in + E_inter_qq_in 
  
  ! If here then no overlap was detected and hence calculate the energies further

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

  CALL Compute_Molecule_Nonbond_Intra_Energy(alive,this_species,E_intra_vdw_in, &
       E_intra_qq_in,intra_overlap)

  delta_e_in = delta_e_in + E_intra_vdw_in + E_intra_qq_in 

  call cpu_time(time0)

  IF (l_charge_in) THEN
     
     ! Note that this call will change cos_mol, sin_mol of alive and this
     ! will have to be restored below while computing the energy of the box_out
     ! without this molecule. 
     
     CALL Compute_Ewald_Reciprocal_Energy_Difference(alive,alive,this_species,box_in, &
          int_insertion,E_reciprocal_move_in)

     CALL Compute_Ewald_Self_Energy_Difference(alive,this_species,box_in,int_insertion, &
          E_self_move_in)
     delta_e_in = delta_e_in + E_self_move_in
     delta_e_in = delta_e_in + E_reciprocal_move_in - energy(box_in)%ewald_reciprocal

  END IF

  call cpu_time(time1)
  copy_time = copy_time + time1-time0

  IF (int_vdw_sum_style(box_in) == vdw_cut_tail) THEN
     nbeads_in(:) = nint_beads(:,box_in)

     DO i = 1, natoms(this_species)
        i_type = nonbond_list(i,this_species)%atom_type_number
        nint_beads(i_type,box_in) = nint_beads(i_type,box_in) + 1
     END DO
        
     CALL Compute_LR_Correction(box_in,e_lrc_in)
     delta_e_in = delta_e_in + e_lrc_in - energy(box_in)%lrc
     
  END IF

  IF(cpcollect) THEN
  
     potw = 1.0_DP
     CP_energy = delta_e_in
     
     IF(species_list(this_species)%fragment) THEN
        potw = 1.0_DP / (P_forward *   kappa_ins*kappa_rot*kappa_dih ** (nfragments(this_species)-1)) 
        CP_energy = delta_e_in - E_angle_in 
     END IF
   
     chpot(this_species,box_in) = chpot(this_species,box_in) + potw * (box_list(box_in)%volume / &
          (REAL(nmols(this_species,box_in)+1))) * DEXP(-beta(box_in) * CP_energy)

  END IF

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

  ! Obtain the weight of the chain if it is made up of fragments

  call cpu_time(time0)
  IF ( species_list(this_species)%fragment .AND. species_list(this_species)%int_insert .NE. int_igas) THEN
     ! Note that we need to use the same fragment order in which we inserted the molecule
     ! above. so in this case frag_order becomes input to the routine. We obtain
     ! P_delete via this call. Note that, cbmc_overlap should be false as we are
     ! dealing with an existing molecule.
     del_Flag = .TRUE. 
     get_fragorder = .FALSE.

     IF (species_list(this_species)%lcom) THEN
        CALL Build_Molecule(alive,this_species,box_out,frag_order,lambda_for_build, which_anchor,P_reverse, &
             nrg_ring_frag_out,cbmc_overlap)
     ELSE
        CALL Build_Rigid_Fragment(alive,this_species,box_out,frag_order,lambda_for_build, &
             which_anchor, P_reverse, nrg_ring_frag_out,cbmc_overlap)
     END IF
        
     IF (cbmc_overlap) THEN
        ! energy overlap was found in the old configuration

        atom_list(1:natoms(this_species),alive,this_species)%exist = .TRUE.
        CALL Revert_Old_Cartesian_Coordinates(alive,this_species)

        molecule_list(alive,this_species)%cfc_lambda = 1.0_DP
        
     END IF


  END IF

  CALL Get_COM(alive,this_species)
  CALL Compute_Max_COM_Distance(alive,this_species)
  ! debug to see if the COM and max_com_distance are identical

  CALL Compute_Molecule_Bond_Energy(alive,this_species,e_bond_out)
  CALL Compute_Molecule_Angle_Energy(alive,this_species,e_angle_out)
  CALL Compute_Molecule_Dihedral_Energy(alive,this_species,e_dihed_out)
  CALL Compute_Molecule_Improper_Energy(alive,this_species,e_improper_out)

  delta_e_out = delta_e_out + e_bond_out + e_angle_out + e_dihed_out + e_improper_out

! Nonbonded energy  

  CALL Compute_Molecule_Nonbond_Intra_Energy(alive,this_species,E_intra_vdw_out, &
       E_intra_qq_out,intra_overlap)

  IF ( .NOT. l_pair_nrg) THEN
  
     CALL Compute_Molecule_Nonbond_Inter_Energy(alive,this_species,E_inter_vdw_out, &
          E_inter_qq_out,inter_overlap)
  END IF

  delta_e_out = delta_e_out + E_intra_vdw_out + E_intra_qq_out + E_inter_vdw_out + E_inter_qq_out

  IF (l_charge_out) THEN

     ! now we will have to restore the cos_mol and sin_mol as they changed above
     ! but restoring will destroy the newly computed vector so now here allocate cos_mol_new
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

     CALL Compute_Ewald_Reciprocal_Energy_Difference(alive,alive,this_species &
          ,box_out,int_deletion,E_reciprocal_move_out)
     CALL Compute_Ewald_Self_Energy_Difference(alive,this_species,box_out, &
          int_deletion,E_self_move_out)

     delta_e_out = delta_e_out - E_self_move_out
     delta_e_out = delta_e_out - (E_reciprocal_move_out - energy(box_out)%ewald_reciprocal)

  END IF
  
  IF (int_vdw_sum_style(box_out) == vdw_cut_tail) THEN

     nbeads_out(:) = nint_beads(:,box_out)

     DO i = 1, natoms(this_species)
        i_type = nonbond_list(i,this_species)%atom_type_number
        nint_beads(i_type,box_out) = nint_beads(i_type,box_out) - 1
     END DO

     CALL Compute_LR_correction(box_out,e_lrc_out)
     delta_e_out = delta_e_out - ( e_lrc_out - energy(box_out)%lrc )

  END IF

  delta_e_in_pacc = delta_e_in
  delta_e_out_pacc = delta_e_out

  IF(species_list(this_species)%int_insert == int_igas) THEN
     igas_flag = .TRUE.
     CALL Compute_Molecule_Nonbond_Intra_Energy(alive,this_species,E_intra_vdw_igas,E_intra_qq_igas,intra_overlap)
     igas_flag = .FALSE. 
     delta_e_out_pacc = delta_e_out_pacc - e_bond_out - e_angle_out - e_dihed_out - e_improper_out - &
                         E_intra_vdw_igas - E_intra_qq_igas
     delta_e_in_pacc = delta_e_in_pacc - energy_igas(rand_igas,this_species)%total
  END IF

  IF(species_list(this_species)%fragment .AND. species_list(this_species)%int_insert .NE. int_igas) THEN
     delta_e_in_pacc = delta_e_in_pacc - e_angle_in - nrg_ring_frag_in
     delta_e_out_pacc = delta_e_out_pacc - e_angle_out - nrg_ring_frag_out
  END IF

  ! Define a factor that will be used to accept or reject the move. Note that
  ! the change in energy of box_out is actually negative of delta_e_out calculated
  ! above

  factor = beta(box_in)*delta_e_in_pacc - beta(box_out)*delta_e_out_pacc

  factor = factor - DLOG((box_list(box_in)%volume * REAL(nmols(this_species,box_out),DP))/ &
                          (box_list(box_out)%volume * REAL(nmols(this_species,box_in) + 1, DP)))
  factor = factor + DLOG(P_forward) - DLOG(P_reverse)

  accept = accept_or_reject(factor)

  IF (accept) THEN
     
     ! Change the number of molecules in each box

     nmols(this_species,box_in) = nmols(this_species,box_in) + 1
     nmols(this_species,box_out) = nmols(this_species,box_out) - 1

     nsuccess(this_species,box_in)%insertion = nsuccess(this_species,box_in)%insertion + 1
     nsuccess(this_species,box_out)%deletion = nsuccess(this_species,box_out)%deletion + 1

     DO i = 1,natoms(this_species)
       atom_list(i,alive,this_species) = new_atom_list(i)
     END DO
     molecule_list(alive,this_species) = new_molecule_list

     CALL Fold_Molecule(alive,this_species,box_in)

     IF (l_charge_in) THEN
        
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

     ! Update energies for box the boxes

     energy(box_in)%total = energy(box_in)%total + delta_e_in
     energy(box_in)%intra = energy(box_in)%intra + e_bond_in + e_angle_in + e_dihed_in
     energy(box_in)%bond = energy(box_in)%bond + e_bond_in
     energy(box_in)%angle = energy(box_in)%angle + e_angle_in
     energy(box_in)%dihedral = energy(box_in)%dihedral + e_dihed_in
     energy(box_in)%intra_vdw = energy(box_in)%intra_vdw + e_intra_vdw_in
     energy(box_in)%intra_q = energy(box_in)%intra_q + e_intra_qq_in
     energy(box_in)%inter_vdw = energy(box_in)%inter_vdw + e_inter_vdw_in
     energy(box_in)%inter_q = energy(box_in)%inter_q + e_inter_qq_in

     IF ( l_charge_in) THEN
        energy(box_in)%ewald_reciprocal = e_reciprocal_move_in
        energy(box_in)%ewald_self = energy(box_in)%ewald_self + e_self_move_in
        
     END IF

     IF ( l_charge_out) THEN

        energy(box_out)%ewald_reciprocal = e_reciprocal_move_out
        energy(box_out)%ewald_self = energy(box_out)%ewald_self + e_self_move_out
        
     END IF

     IF (int_vdw_sum_style(box_in) == vdw_cut_tail) THEN
        energy(box_in)%lrc = e_lrc_in
     END IF
     IF (int_vdw_sum_style(box_out) == vdw_cut_tail) THEN
        energy(box_out)%lrc = e_lrc_out
     END IF

     ! for box_out
     
     energy(box_out)%total = energy(box_out)%total - delta_e_out
     energy(box_out)%intra = energy(box_out)%intra - e_bond_out - e_angle_out - e_dihed_out
     energy(box_out)%bond = energy(box_out)%bond - e_bond_out
     energy(box_out)%angle = energy(box_out)%angle - e_angle_out
     energy(box_out)%dihedral = energy(box_out)%dihedral - e_dihed_out
     energy(box_out)%intra_vdw = energy(box_out)%intra_vdw - e_intra_vdw_out
     energy(box_out)%intra_q   = energy(box_out)%intra_q - e_intra_qq_out
     energy(box_out)%inter_vdw = energy(box_out)%inter_vdw - e_inter_vdw_out
     energy(box_out)%inter_q   = energy(box_out)%inter_q - e_inter_qq_out

  ELSE
     
     ! reject the move. Atomic coordinates have not changed as we used the original coordinate
     ! in box_out to calculate removal energies

     ! Restore the reciprocal space k vector


    ! Restore the reciprocal space k vector

     IF (l_charge_in) THEN
        
        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        cos_sum(1:nvecs(box_in),box_in) = cos_sum_old(1:nvecs(box_in),box_in)
        sin_sum(1:nvecs(box_in),box_in) = sin_sum_old(1:nvecs(box_in),box_in)
        !$OMP END PARALLEL WORKSHARE

        DEALLOCATE(cos_mol_new,sin_mol_new)

     END IF

     IF (l_charge_out) THEN
        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        cos_sum(1:nvecs(box_out),box_out) = cos_sum_old(1:nvecs(box_out),box_out)
        sin_sum(1:nvecs(box_out),box_out) = sin_sum_old(1:nvecs(box_out),box_out)
        
        cos_mol(1:nvecs(box_out),position) = cos_mol_old(:)
        sin_mol(1:nvecs(box_out),position) = sin_mol_old(:)
        !$OMP END PARALLEL WORKSHARE

        DEALLOCATE(cos_mol_old)
        DEALLOCATE(sin_mol_old)

     END IF

     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(alive,this_species,box_out)
     

     IF ( int_vdw_sum_style(box_in) == vdw_cut_tail ) THEN

        nint_beads(:,box_in) = nbeads_in(:)

     END IF

     IF ( int_vdw_sum_style(box_out) == vdw_cut_tail ) THEN

        nint_beads(:,box_out) = nbeads_out(:)

     END IF

  END IF
  

  ! Acceptance ratio
  ! Update various values
  ! 

  IF (ALLOCATED(sorbate_id)) DEALLOCATE(sorbate_id)
  IF (ALLOCATED(sorbate_x)) DEALLOCATE(sorbate_x)
  DEALLOCATE(new_atom_list)

END SUBROUTINE GEMC_Particle_Transfer

SUBROUTINE New_Positions(this_box,alive,is,rand_igas)

  USE Run_Variables
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


