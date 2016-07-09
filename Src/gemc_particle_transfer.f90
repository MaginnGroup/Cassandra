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
  ! Step 1) Select a box 'box_out' from which to remove a molecule:
  !          a) according to its mole fraction, and box_in with uniform prob (DEFAULT)
  !          b) with uniform probability, and box_in with probs from input file
  ! Step 2) Select a species 'is':
  !          a) according to its mole fraction in box_out (DEFAULT)
  !          b) using probabilities in the input file
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

  INTEGER :: box_in, box_out, i, k, i_type

  INTEGER :: is, ibox, im_in, im_out, alive, which_anchor
  INTEGER :: rand_igas, locate_in
  INTEGER :: nmols_tot, nmols_box(nbr_boxes)
  REAL(DP) :: x_box(nbr_boxes), x_is(nspecies)

  INTEGER, ALLOCATABLE :: frag_order(:)

  REAL(DP), ALLOCATABLE :: dx(:), dy(:), dz(:)

  REAL(DP) :: delta_e_out, delta_e_out_pacc, e_bond_out, e_angle_out, e_dihed_out
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
  LOGICAL :: inter_overlap, accept_or_reject, cbmc_overlap
  LOGICAL :: intra_overlap
 ! ring biasing variables

  REAL(DP) :: nrg_ring_frag_in, nrg_ring_frag_out

  TYPE(atom_class), ALLOCATABLE :: new_atom_list(:)
  TYPE(molecule_class) :: new_molecule_list

 ! Variables added for l_pair_nrg and reciprocal k vector storage

  INTEGER :: position

  REAL(DP), ALLOCATABLE :: cos_mol_old(:), sin_mol_old(:), cos_mol_new(:), sin_mol_new(:)
  REAL(DP) :: time0, time1, randno
  
  LOGICAL :: l_charge_in, l_charge_out

  potw = 1.0_DP
  inter_overlap = .false.
  cbmc_overlap = .false.
  accept = .false.

  P_seq = 1.0_DP
  P_forward = 1.0_DP
  P_reverse = 1.0_DP

  !*****************************************************************************
  ! Step 1) Select a box 'box_out' from which to remove a particle
  !*****************************************************************************
  IF (.NOT. l_prob_swap_from_box) THEN
     ! Choose box_out based on its mol fraction

     ! Sum the number of swappable molecules total & per box
     nmols_tot = 0 ! sum over species, box
     nmols_box = 0
     DO ibox = 1, nbr_boxes
       DO is = 1, nspecies
         ! Only count swappable species
         IF ( species_list(is)%int_insert /= int_noinsert ) THEN
           nmols_tot = nmols_tot + nmols(is,ibox)
           nmols_box(ibox) = nmols_box(ibox) + nmols(is,ibox)
         END IF
       END DO
     END DO

     ! If there are no molecules then return
     IF (nmols_tot == 0) RETURN

     ! Need cumulative mol fractions for Golden sampling
     DO ibox = 1, nbr_boxes
        x_box(ibox) = REAL(nmols_box(ibox),DP)/REAL(nmols_tot,DP)
        IF ( ibox > 1 ) THEN
           x_box(ibox) = x_box(ibox) + x_box(ibox-1)
        END IF
     END DO

     ! Ready to choose box_out
     randno = rranf()
     DO ibox = 1, nbr_boxes
        IF ( randno <= x_box(ibox)) EXIT
     END DO
     box_out = ibox

     ! If there are no molecules in this box then return
     IF( nmols_box(box_out) == 0 ) RETURN

     ! Choose box_in with uniform probability
     ! Since the number of boxes is constant, this step does not change P_forward or P_reverse
     ibox = INT(rranf() * (nbr_boxes-1)) + 1
     IF (ibox >= box_out) ibox = ibox + 1
     box_in = ibox

     P_forward = P_forward * REAL(nmols_box(box_out),DP)/REAL(nmols_tot,DP)
     P_reverse = P_reverse * REAL(nmols_box(box_in)+1,DP)/REAL(nmols_tot,DP)

  ELSE
     ! Use the probabilities in the input file

     ! choose a donor box
     randno = rranf()
     DO ibox = 1,nbr_boxes
        IF (randno <= cum_prob_swap_from_box(ibox)) EXIT
     END DO
     box_out = ibox

     ! choose an acceptor box
     randno = rranf()
     DO ibox = 1,nbr_boxes
        IF (ibox == box_out) CYCLE
        IF (randno <= cum_prob_swap_to_box(box_out,ibox)) EXIT
     END DO
     box_in = ibox

     P_forward = P_forward * prob_swap_from_box(box_out) &
                           * prob_swap_to_box(box_out,box_in)
     P_reverse = P_reverse * prob_swap_from_box(box_in) &
                           * prob_swap_to_box(box_in,box_out)

  END IF

  ! increment total number of trials for each of the boxes
  tot_trials(box_out) = tot_trials(box_out) + 1
  tot_trials(box_in) = tot_trials(box_in) + 1

  !*****************************************************************************
  ! Step 2) Select a species 'is' according to its overall mole 
  !         fraction
  !*****************************************************************************
  IF (.NOT. l_prob_swap_species) THEN
     
     ! Choose a species based on the overall mole fraction.
     ! Need cumulative mol fractions for Golden sampling
     x_is = 0.0_DP
     DO is = 1, nspecies
        IF (species_list(is)%int_insert /= int_noinsert) THEN
           x_is(is) = REAL(nmols(is,box_out),DP)/REAL(nmols_box(box_out),DP)
        END IF
        IF ( is > 1 ) THEN
           x_is(is) = x_is(is) + x_is(is-1)
        END IF
     END DO
 
     ! Ready to choose is
     randno = rranf()
     DO is = 1, nspecies
        IF (randno < x_is(is)) EXIT
     END DO
     
     P_forward = P_forward * REAL(nmols(is,box_out),DP) / REAL(nmols_box(box_out),DP)
     P_reverse = P_reverse * REAL(nmols(is,box_in)+1,DP) / REAL(nmols_box(box_in)+1,DP)

  ELSE
     
     ! pick the species based on specified probability
     randno = rranf()
     DO is = 1, nspecies
        IF (randno <= cum_prob_swap_species(is)) EXIT
     END DO
     
     ! Since prob_swap_species is constant, it will be the same for forward and reverse move
     ! and therefore does not change P_forward / P_reverse
  END IF

  IF (nmols(is,box_out) == 0) THEN
     IF (cpcollect)  CALL Chempot(box_in)
     RETURN
  END IF
  
  ! Increment counters
  ntrials(is,box_out)%deletion = &
     ntrials(is,box_out)%deletion + 1
  ntrials(is,box_in)%insertion = &
     ntrials(is,box_in)%insertion + 1
  ntrials(is,box_in)%cpcalc = &
     ntrials(is,box_in)%cpcalc + 1

  !*****************************************************************************
  ! Step 3) Select a molecule 'alive' in box_out with uniform probability
  !*****************************************************************************
  ! pick a molecule INDEX at random to delete
  im_out = INT(rranf() * nmols(is,box_out)) + 1
  ! the probability of picking im_out will be accounted for when computing ln_pacc
  
  ! Obtain the LOCATE of this molecule
  alive = locate(im_out,is,box_out)

  !*****************************************************************************
  ! Step 4) Choose a position, orientation and conformation for alive in box_in
  !*****************************************************************************
  ! Insert alive into box_in before deleting alive from box_out.
  ! If there is an overlap in box_in, the move can be
  ! immediately rejected. Note that this move changes the box identity
  ! and coordinate related quantites of 'alive'. When we calculate
  ! the energy of removing molecule from box_out we need to restore these

  ! Save the coordinates of 'alive' in 'box_out'
  CALL Save_Old_Cartesian_Coordinates(alive,is)
  CALL Compute_Molecule_Dihedral_Energy(alive,is,e_dihed_out)

  ! Save the interaction energies
  IF (l_pair_nrg) CALL Store_Molecule_Pair_Interaction_Arrays(alive,is, &
       box_out, E_inter_vdw_out, E_inter_qq_out)
  
  ! Save the k-vectors
  
  IF (int_charge_sum_style(box_in)  == charge_ewald .AND.&
      has_charge(is)) THEN
     ALLOCATE(cos_mol_old(nvecs(box_out)), sin_mol_old(nvecs(box_out)))
     CALL Get_Position_Alive(alive,is,position)
     
     !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
     cos_mol_old(:) = cos_mol(1:nvecs(box_out),position)
     sin_mol_old(:) = sin_mol(1:nvecs(box_out),position)
     !$OMP END PARALLEL WORKSHARE
  END IF

  ! Switch the box identity of alive
  molecule_list(alive,is)%which_box = box_in
  nmols(is,box_in) = nmols(is,box_in) + 1
  im_in = nmols(is,box_in) ! INDEX of alive in box_in
  locate(im_in,is,box_in) = alive ! link INDEX to LOCATE

  ! set deletion flag to false and call the CBMC growth routine
  cbmc_overlap = .FALSE.

  delta_e_in = 0.0_DP

  IF (species_list(is)%fragment .AND. &
      species_list(is)%int_insert .NE. int_igas ) THEN

     del_flag = .FALSE.
     get_fragorder = .TRUE.
     ALLOCATE(frag_order(nfragments(is)))
     lambda_for_build = molecule_list(alive,is)%frac
     CALL Build_Molecule(alive,is,box_in,frag_order, &
             lambda_for_build,P_seq,P_forward,nrg_ring_frag_in, &
             cbmc_overlap)
  ELSE
     CALL New_Positions(box_in,alive,is,rand_igas)
  END IF

  CALL Get_COM(alive,is)
  CALL Compute_Max_COM_Distance(alive,is)
  CALL Fold_Molecule(alive,is,box_in)

  ! Build_Molecule returns cbmc_overlap == .TRUE. both in case of a core
  ! overlap and also when the weight of all trials is zero.

  IF (.NOT. cbmc_overlap) THEN
        CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is, &
             E_inter_vdw_in,E_inter_qq_in,inter_overlap)
  END IF

  IF (inter_overlap .OR. cbmc_overlap) THEN
     ! reject the swap

     ! restore the box identity
     molecule_list(alive,is)%which_box = box_out
     locate(im_in,is,box_in) = 0 ! reset LOCATE
     nmols(is,box_in) = nmols(is,box_in) - 1
     ! restore the box_out coordinates
     CALL Revert_Old_Cartesian_Coordinates(alive,is)

     ! All atoms will not exist if inter_overlap was tripped before the 
     ! last fragment was placed in Build_Molecule. 
     ! Set exist to TRUE for all atoms and reset frac to 1
     atom_list(:,alive,is)%exist = .TRUE.
     molecule_list(alive,is)%frac = 1.0_DP

     IF (l_pair_nrg) THEN
        CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,box_out)
     END IF

     IF(ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
     IF(ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)

     accept = .FALSE.

  ELSE

    delta_e_in = delta_e_in + E_inter_vdw_in + E_inter_qq_in 
    
    !*****************************************************************************
    ! Step 5) Calculate the change in box_in's potential energy from inserting
    !         alive
    !*****************************************************************************
    ! If here then no overlap was detected. Calculate the rest of the energies
    IF(species_list(is)%int_insert == int_random) THEN

       CALL Compute_Molecule_Bond_Energy(alive,is,E_bond_in)
       CALL Compute_Molecule_Angle_Energy(alive,is,E_angle_in)
       CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihed_in)
       CALL Compute_Molecule_Improper_Energy(alive,is,E_improper_in)
    
    ELSE IF(species_list(is)%int_insert == int_igas) THEN

       E_bond_in = energy_igas(rand_igas,is)%bond
       E_angle_in = energy_igas(rand_igas,is)%angle
       E_dihed_in = energy_igas(rand_igas,is)%dihedral
       E_improper_in = energy_igas(rand_igas,is)%improper
    
    END IF

    delta_e_in = delta_e_in + E_bond_in + E_angle_in + E_dihed_in

    CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is, &
         E_intra_vdw_in,E_intra_qq_in,E_periodic_qq,intra_overlap)
    E_inter_qq_in = E_inter_qq_in + E_periodic_qq

    ! Already added E_inter to delta_e, so add E_periodic directly
    delta_e_in = delta_e_in + E_intra_vdw_in + E_intra_qq_in + E_periodic_qq

    call cpu_time(time0)

    IF (int_charge_style(box_in) == charge_coul .AND. has_charge(is)) THEN
       
       IF (int_charge_sum_style(box_in) == charge_ewald) THEN

            ! Note that this call will change cos_mol, sin_mol of alive and this
            ! will have to be restored below while computing the energy of box_out
            ! without molecule alive. 
            CALL Update_System_Ewald_Reciprocal_Energy(alive,is,box_in, &
                 int_insertion,E_reciprocal_in)
       
            delta_e_in = delta_e_in + (E_reciprocal_in - energy(box_in)%ewald_reciprocal)
       END IF
       
       CALL Compute_Molecule_Self_Energy(alive,is,box_in, &
            E_self_in)
       delta_e_in = delta_e_in + E_self_in 

    END IF

    call cpu_time(time1)
    copy_time = copy_time + time1-time0

    IF (int_vdw_sum_style(box_in) == vdw_cut_tail) THEN
       nbeads_in(:) = nint_beads(:,box_in)

       DO i = 1, natoms(is)
          i_type = nonbond_list(i,is)%atom_type_number
          nint_beads(i_type,box_in) = nint_beads(i_type,box_in) + 1
       END DO
          
       CALL Compute_LR_Correction(box_in,e_lrc_in)
       delta_e_in = delta_e_in + e_lrc_in - energy(box_in)%lrc

    END IF

    IF(cpcollect) THEN
    
       potw = 1.0_DP
       CP_energy = delta_e_in
       
       IF(species_list(is)%fragment) THEN
          potw = 1.0_DP / (P_forward * kappa_ins*kappa_rot*kappa_dih &
               ** (nfragments(is)-1)) 
          CP_energy = delta_e_in - E_angle_in 
       END IF
     
       chpot(is,box_in) = chpot(is,box_in) &
        + potw * (box_list(box_in)%volume &
        / (REAL(nmols(is,box_in)))) * DEXP(-beta(box_in) * CP_energy)

    END IF

    !*****************************************************************************
    ! Step 6) Calculate the change in box_out's potential energy from deleting
    !         alive
    !*****************************************************************************
    ! Need to preserve coordinates of 'alive' in box_in.
    ALLOCATE(new_atom_list(natoms(is)))
    new_atom_list(:) = atom_list(:,alive,is)
    new_molecule_list = molecule_list(alive,is)

    ! Change the box identity of alive back to box_out
    molecule_list(alive,is)%which_box = box_out
    ! Restore the box_out coordinates
    CALL Revert_Old_Cartesian_Coordinates(alive,is)

    delta_e_out = 0.0_DP

    ! Obtain the weight of the chain if it is made up of fragments
    call cpu_time(time0)
    IF ( species_list(is)%fragment .AND. &
         species_list(is)%int_insert .NE. int_igas) THEN
       ! The fragment order was decided when inserting alive into box_in
       ! Use the same fragment order to calculate trial insertions into box_out 
       ! 
       ! So frag_order and P_seq are inputs to the Build_Molecule routine
       ! We obtain P_reverse via this call. Note that, cbmc_overlap 
       ! must be false as we are dealing with an existing molecule.
       del_flag = .TRUE. 
       get_fragorder = .FALSE.

       CALL Build_Molecule(alive,is,box_out,frag_order, &
               lambda_for_build,P_seq,P_reverse,nrg_ring_frag_out, &
               cbmc_overlap)
          
       IF (cbmc_overlap) THEN
          ! If this flag gets tripped, there is an error in the code
          err_msg = ""
          err_msg(1) = "Error: existing configuration of " // & 
             "molecule " // TRIM(Int_To_String(alive))
          err_msg(2) = "of species " // TRIM(Int_To_String(is)) // &
             " in box " // TRIM(Int_To_String(box_out)) // &
             " tripped an overlap error"
          CALL Clean_Abort(err_msg,'GEMC_Particle_Transfer')

!          atom_list(1:natoms(is),alive,is)%exist = .TRUE.
!          CALL Revert_Old_Cartesian_Coordinates(alive,is)
!
!          molecule_list(alive,is)%frac = 1.0_DP
       END IF

    END IF

    CALL Get_COM(alive,is)
    CALL Compute_Max_COM_Distance(alive,is)
    ! debug to see if the COM and max_com_distance are identical

    ! bonded energies
    CALL Compute_Molecule_Bond_Energy(alive,is,e_bond_out)
    CALL Compute_Molecule_Angle_Energy(alive,is,e_angle_out)
    CALL Compute_Molecule_Dihedral_Energy(alive,is,e_dihed_out)
    CALL Compute_Molecule_Improper_Energy(alive,is,e_improper_out)

    delta_e_out = delta_e_out - e_bond_out - e_angle_out - e_dihed_out &
                - e_improper_out

    ! Nonbonded energy  
    CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is, &
         E_intra_vdw_out,E_intra_qq_out,E_periodic_qq,intra_overlap)

    IF ( .NOT. l_pair_nrg) THEN
       CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is, &
            E_inter_vdw_out,E_inter_qq_out,inter_overlap)
    END IF
    E_inter_qq_out = E_inter_qq_out + E_periodic_qq

    delta_e_out = delta_e_out - E_intra_vdw_out - E_intra_qq_out &
                - E_inter_vdw_out - E_inter_qq_out




    IF (int_charge_style(box_out) == charge_coul .AND. has_charge(is)) THEN
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
        
             CALL Update_System_Ewald_Reciprocal_Energy(alive,is, &
                  box_out,int_deletion,E_reciprocal_out)
        
             delta_e_out = delta_e_out + (E_reciprocal_out - energy(box_out)%ewald_reciprocal)
        
          END IF
        
          CALL Compute_Molecule_Self_Energy(alive,is,box_out, &
            E_self_out)

          delta_e_out = delta_e_out - E_self_out

    END IF

    
    IF (int_vdw_sum_style(box_out) == vdw_cut_tail) THEN
       nbeads_out(:) = nint_beads(:,box_out)
       DO i = 1, natoms(is)
          i_type = nonbond_list(i,is)%atom_type_number
          nint_beads(i_type,box_out) = nint_beads(i_type,box_out) - 1
       END DO

       CALL Compute_LR_correction(box_out,e_lrc_out)
       delta_e_out = delta_e_out + ( e_lrc_out - energy(box_out)%lrc )

    END IF

    delta_e_in_pacc = delta_e_in
    delta_e_out_pacc = delta_e_out

    IF(species_list(is)%int_insert == int_igas) THEN
       igas_flag = .TRUE.
       CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is, &
            E_intra_vdw_igas,E_intra_qq_igas,E_periodic_qq,intra_overlap)
       ! Ideal gas should not interact with it's periodic image, so E_periodic_qq
       ! is not added the delta_e
       igas_flag = .FALSE. 
       delta_e_out_pacc = delta_e_out_pacc + e_bond_out + e_angle_out &
                        + e_dihed_out + e_improper_out &
                        + E_intra_vdw_igas + E_intra_qq_igas
       delta_e_in_pacc = delta_e_in_pacc &
                       - energy_igas(rand_igas,is)%total
    END IF

    IF(species_list(is)%fragment .AND. &
       species_list(is)%int_insert .NE. int_igas) THEN
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
                      - DLOG(REAL(nmols(is,box_out),DP)) &
                      + DLOG(REAL(nmols(is,box_in), DP))

    ! The same order of insertion is used in both the insertion and 
    ! reverse deletion, so P_seq does not factor into ln_pacc
    ln_pacc = ln_pacc + DLOG(P_forward / P_reverse)

    accept = accept_or_reject(ln_pacc)

    IF (accept) THEN
       ! accept the swap
       
       ! already updated the number of molecules in box_in

       ! remove the deleted molecule from box_out locate
       IF (im_out < nmols(is,box_out)) THEN
          DO k = im_out + 1, nmols(is,box_out)
             locate(k-1,is,box_out) = locate(k,is,box_out)
          END DO
       END IF
       locate(nmols(is,box_out),is,box_out) = 0

       ! Update the number of molecules in box_out
       nmols(is,box_out) = nmols(is,box_out) - 1

       ! Set the coordinates and properties of molecule alive to box_in values
       DO i = 1,natoms(is)
         atom_list(i,alive,is) = new_atom_list(i)
       END DO
       molecule_list(alive,is) = new_molecule_list
       CALL Fold_Molecule(alive,is,box_in)

       IF (int_charge_sum_style(box_in) == charge_ewald .AND. &
           has_charge(is)) THEN
          call cpu_time(time0)
          !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
          cos_mol(1:nvecs(box_in),position) = cos_mol_new(:)
          sin_mol(1:nvecs(box_in),position) = sin_mol_new(:)
          !$OMP END PARALLEL WORKSHARE
          
          DEALLOCATE(cos_mol_new,sin_mol_new)
          
          call cpu_time(time1)
!          copy_time = copy_time + time1-time0
       END IF

       IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
       IF (ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
       IF (ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)
       IF (ALLOCATED(cos_mol_new)) DEALLOCATE(cos_mol_new)
       IF (ALLOCATED(sin_mol_new)) DEALLOCATE(sin_mol_new)

       ! Restore the coordinates of the molecule due to successful insertion
       CALL Get_Internal_Coordinates(alive,is)

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

       IF (has_charge(is)) THEN
          IF (int_charge_sum_style(box_in) == charge_ewald) energy(box_in)%ewald_reciprocal = E_reciprocal_in        
          IF (int_charge_sum_style(box_out) == charge_ewald) energy(box_out)%ewald_reciprocal = E_reciprocal_out
          IF (int_charge_style(box_in) == charge_coul) energy(box_in)%self = energy(box_in)%self + E_self_in
          IF (int_charge_style(box_out) == charge_coul) energy(box_out)%self = energy(box_out)%self + E_self_out
       END IF


       ! Increment counter
       nsuccess(is,box_in)%insertion = &
          nsuccess(is,box_in)%insertion + 1
       nsuccess(is,box_out)%deletion = &
          nsuccess(is,box_out)%deletion + 1

    ELSE
       ! reject the swap. 

       ! Atomic coordinates have not changed as we used the original coordinate
       ! in box_out to calculate removal energies

       ! Reset the number of molecules and locate of box_in
       locate(im_in,is,box_in) = 0
       nmols(is,box_in) = nmols(is,box_in) - 1


       IF (has_charge(is)) THEN
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
          CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,box_out)
       END IF

       IF ( int_vdw_sum_style(box_in) == vdw_cut_tail ) THEN
          nint_beads(:,box_in) = nbeads_in(:)
       END IF

       IF ( int_vdw_sum_style(box_out) == vdw_cut_tail ) THEN
          nint_beads(:,box_out) = nbeads_out(:)
       END IF


    END IF
  
  END IF 

  IF (ALLOCATED(new_atom_list)) DEALLOCATE(new_atom_list)

  IF (verbose_log) THEN
    WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8)') i_mcstep, 'swap_to' , alive, is, box_in, accept
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


