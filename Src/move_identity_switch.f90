
    !List of parameters I'll need to figure out
    !x_box, or the box
    !


    !local variables

    !variables to find elsewhere:
    !SPECIES_LIST()

    !No need to choose a box, the box is implicit


    !Step 1) Select a box

    !Step 2) Select a species 'is' that exists in that box:

    !Step 3) Select a molecule from that species in that box

    !Step 3) Select a second box:

    !Step 4) Select a second species that exists in that box

    !Step 7) Select a second molecule from that species in that box

    !a) according to its mole fraction in the one box

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

SUBROUTINE IDENTITY_EXCHANGE

  !*****************************************************************************
  !
  ! This subroutine performs particle swaps in a GEMC ensemble.
  !
  !*****************************************************************************

  !*Unchanged!
  USE Global_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE IO_Utilities
  USE Fragment_Growth
  USE Pair_Nrg_Routines
  !*/Unchanged!

  IMPLICIT NONE

  INTEGER :: box_in, box_out, i, k, i_type
  !added
  INTEGER :: box
  INTEGER :: i_alive, j_alive
  !end added

  INTEGER :: is, ibox, im_in, im_out, alive, which_anchor
  INTEGER :: rand_igas, locate_in
  INTEGER :: nmols_tot, nmols_box(nbr_boxes)
  REAL(DP) :: x_box(nbr_boxes), x_is(nspecies)

  !new
  INTEGER :: js
  !end new

  INTEGER, ALLOCATABLE :: frag_order(:)

  REAL(DP), ALLOCATABLE :: dx(:), dy(:), dz(:)

  REAL(DP) :: dE_out, dE_out_pacc
  REAL(DP) :: E_bond_out, E_angle_out, E_dihed_out, E_improper_out
  REAL(DP) :: E_intra_vdw_out, E_intra_qq_out, E_periodic_qq
  REAL(DP) :: E_inter_vdw_out, E_inter_qq_out
  REAL(DP) :: E_reciprocal_out, E_self_out, E_lrc_out
  REAL(DP) :: dE_in, dE_in_pacc
  REAL(DP) :: E_bond_in, E_angle_in, E_dihed_in, E_improper_in
  REAL(DP) :: E_intra_vdw_in, E_intra_qq_in
  REAL(DP) :: E_inter_vdw_in, E_inter_qq_in
  REAL(DP) :: E_reciprocal_in, E_self_in, E_lrc_in

  REAL(DP) :: potw, CP_energy
  REAL(DP) :: ln_pseq, ln_pfor, ln_prev, P_forward, P_reverse, ln_pacc
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

  ln_pseq = 0.0_DP
  ln_pfor = 0.0_DP
  ln_prev = 0.0_DP
  P_forward = 1.0_DP
  P_reverse = 1.0_DP

  tot_trials(box) = tot_trials(box) + 1
  !tot_trials(box_in) = tot_trials(box_in) + 1


  !similar variables to move_mol_swap
  ibox = x_box(1)

  !*****************************************************************************
  ! Step 1) Select a species 'is':
  !      -> Right now just done randomly!
  !*****************************************************************************

  ! Sum the number of swappable molecules total & per box
  nmols_tot = 0 ! sum over species, box
  DO is = 1, nspecies
    ! Only count swappable species
    IF ( species_list(is)%int_insert /= int_noinsert ) THEN
      nmols_tot = nmols_tot + nmols(is,ibox)
    END IF
  END DO

  ! check to make sure there are swappable molecules
  IF (nmols_tot == 0) THEN
     err_msg = ''
     err_msg(1) = 'No swappable molecules'
     CALL Clean_Abort(err_msg,'Identity_Switch_Move')
  END IF

  !check this
  is = INT(rranf() * nspecies) + 1


  !TODO: ONLY ERROR OUT IF SWAP ISN'T POSSIBLE
  ! check to make sure the selected species is swappable
  IF (species_list(is)%int_insert == int_noinsert) THEN
    err_msg = ''
    err_msg(1) = 'Species ' // TRIM(Int_To_String(is)) // ' is not swappable'
    CALL Clean_Abort(err_msg,'Identity_Switch_Move')
  END IF

  !*****************************************************************************
  ! Step 2) Select a molecule 'alive' from species 'is' with uniform probability
  !*****************************************************************************
  ! pick a molecule INDEX at random
  im = INT(rranf() * nmols(is,box)) + 1

  ! Obtain the LOCATE of this molecule
  i_alive = locate(im, is, box)

  !*****************************************************************************
  ! Step 3) Select a species 'js':
  !      -> Right now just done randomly!
  !*****************************************************************************

  js = INT(rranf() * (nspecies - 1)) + 1

  !avoid choosing same species as before
  IF (js >= is)
    js = js + 1
  END IF

  !*****************************************************************************
  ! Step 4) Select a molecule 'alive' from species 'js' with uniform probability
  !*****************************************************************************
  jm = INT(rranf() * nmols(js,box)) + 1

  j_alive = locate(jm, js, box)

  !*****************************************************************************
  ! Step 5) Calculate the change in box_in's potential energy from inserting
  !
  !*****************************************************************************
  ! Save the coordinates of both 'alive' molecules
  ! box?

  CALL Save_Old_Cartesian_Coordinates(i_alive,is)
  CALL Save_Old_Cartesian_Coordiantes(j_alive,js)

  !Needed?
  CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihed_i)
  CALL Compute_Molecule_Dihedral_Energy(alive,js,E_dihed_j)

  ! Save the coordinates of 'alive' in 'box'
  !CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihed_out)
  ! Returns
  !REAL(DP) :: energy_dihed

  ! Save the interaction energies
  IF (l_pair_nrg)
    CALL Store_Molecule_Pair_Interaction_Arrays(i_alive,is, &
       box, E_inter_vdw_out, E_inter_qq_out)
    CALL Store_Molecule_Pair_Interaction_Arrays(j_alive,js, &
       box, E_inter_vdw_out, E_inter_qq_out)
  END IF
  ! Save the k-vectors
  IF (int_charge_sum_style(box)  == charge_ewald .AND.&
      has_charge(is)) THEN
     ALLOCATE(cos_mol_old(nvecs(box)), sin_mol_old(nvecs(box)))
     CALL Get_Position_Alive(i_alive,is,position)

     !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
     cos_mol_old(:) = cos_mol(1:nvecs(box),position)
     sin_mol_old(:) = sin_mol(1:nvecs(box),position)
     !$OMP END PARALLEL WORKSHARE
  END IF

  IF (int_charge_sum_style(box)  == charge_ewald .AND.&
      has_charge(js)) THEN
     ALLOCATE(cos_mol_old(nvecs(box)), sin_mol_old(nvecs(box)))
     CALL Get_Position_Alive(j_alive,js,position)

     !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
     cos_mol_old(:) = cos_mol(1:nvecs(box),position)
     sin_mol_old(:) = sin_mol(1:nvecs(box),position)
     !$OMP END PARALLEL WORKSHARE
  END IF

  ! Switch the box identity of alive
  !molecule_list(alive,is)%which_box = box_in
  !nmols(is,box_in) = nmols(is,box_in) + 1
  !im_in = nmols(is,box_in) ! INDEX of alive in box_in
  !locate(im_in,is,box_in) = alive ! link INDEX to LOCATE

  ! set deletion flag to false and call the CBMC growth routine
  !remove cbmc for now
  !cbmc_overlap = .FALSE.

  dE_in = 0.0_DP


  !del_flag = .FALSE.
  !get_fragorder = .TRUE.
  !ALLOCATE(frag_order(nfragments(is)))
  !lambda_for_build = molecule_list(alive,is)%frac
  !CALL Build_Molecule(alive,is,box_in,frag_order, &
  !        lambda_for_build,ln_pseq,ln_pfor,nrg_ring_frag_in, &
  !        cbmc_overlap)

  !gets com of a molecule, stored in molecule_list(alive,is)%xcom, %ycom, %zcom
  CALL Get_COM(i_alive,is)

  !gets distance between is' COM and farthest away molecule
  !stored in molecule_list(alive,is)%max_dcom
  CALL Compute_Max_COM_Distance(i_alive,is)

  !If outside of the box, fold back into the box
  CALL Fold_Molecule(i_alive,is,box)

  CALL Get_COM(j_alive, js)
  CALL Compute_Max_COM_Distance(j_alive, js)
  CALL Fold_Molecule(j_alive, js, box)

  ! Build_Molecule returns cbmc_overlap == .TRUE. both in case of a core
  ! overlap and also when the weight of all trials is zero.

  !Need this?
  !IF (.NOT. cbmc_overlap) THEN
  !      CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is, &
  !           E_inter_vdw_in,E_inter_qq_in,inter_overlap)
  !END IF

  IF (inter_overlap .OR. cbmc_overlap) THEN
     ! reject the swap

     ! restore the box identity
     molecule_list(alive,is)%which_box = box
     locate(im_in,is,box) = 0 ! reset LOCATE
     nmols(is,box) = nmols(is,box) - 1
     ! restore the box coordinates
     CALL Revert_Old_Cartesian_Coordinates(alive,is)

     ! All atoms will not exist if inter_overlap was tripped before the
     ! last fragment was placed in Build_Molecule.
     ! Set exist to TRUE for all atoms and reset frac to 1
     atom_list(:,alive,is)%exist = .TRUE.
     molecule_list(alive,is)%frac = 1.0_DP

     IF (l_pair_nrg) THEN
        CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,box)
     END IF

     IF(ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
     IF(ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)

     accept = .FALSE.

     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I1,A1,I1,X,L8,X,9X,X,A9)') &
             i_mcstep, 'swap' , alive, is, box, '>', box, accept, 'overlap'
     END IF

  ELSE

    dE_in = dE_in + E_inter_vdw_in + E_inter_qq_in

    !*****************************************************************************
    ! Step 5) Calculate the change in box_in's potential energy from inserting
    !         alive
    !*****************************************************************************
    ! If here then no overlap was detected. Calculate the rest of the energies
    CALL Compute_Molecule_Bond_Energy(alive,is,E_bond_in)
    CALL Compute_Molecule_Angle_Energy(alive,is,E_angle_in)
    CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihed_in)
    CALL Compute_Molecule_Improper_Energy(alive,is,E_improper_in)

    dE_in = dE_in + E_bond_in + E_angle_in + E_dihed_in + E_improper_in

    CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is, &
         E_intra_vdw_in,E_intra_qq_in,E_periodic_qq,intra_overlap)
    E_inter_qq_in = E_inter_qq_in + E_periodic_qq

    ! Already added E_inter to dE, so add E_periodic directly
    dE_in = dE_in + E_intra_vdw_in + E_intra_qq_in + E_periodic_qq

    call cpu_time(time0)

    IF (int_charge_style(box_in) == charge_coul .AND. has_charge(is)) THEN

       IF (int_charge_sum_style(box_in) == charge_ewald) THEN

            ! Note that this call will change cos_mol, sin_mol of alive and this
            ! will have to be restored below while computing the energy of box_out
            ! without molecule alive.
            CALL Update_System_Ewald_Reciprocal_Energy(alive,is,box_in, &
                 int_insertion,E_reciprocal_in)

            dE_in = dE_in + (E_reciprocal_in - energy(box_in)%ewald_reciprocal)
       END IF

       CALL Compute_Molecule_Self_Energy(alive,is,box_in,E_self_in)
       dE_in = dE_in + E_self_in

    END IF

    call cpu_time(time1)
    copy_time = copy_time + time1-time0

    IF (int_vdw_sum_style(box_in) == vdw_cut_tail) THEN
       nbeads_in(:) = nint_beads(:,box_in)

       DO i = 1, natoms(is)
          i_type = nonbond_list(i,is)%atom_type_number
          nint_beads(i_type,box_in) = nint_beads(i_type,box_in) + 1
       END DO

       CALL Compute_LR_Correction(box_in,E_lrc_in)
       dE_in = dE_in + E_lrc_in - energy(box_in)%lrc

    END IF

    IF(cpcollect) THEN

       potw = 1.0_DP
       CP_energy = dE_in

       potw = 1.0_DP / (P_forward * kappa_ins*kappa_rot*kappa_dih &
            ** (nfragments(is)-1))
       CP_energy = dE_in - E_angle_in

       chpot(is,box_in) = chpot(is,box_in) &
        + potw * (box_list(box_in)%volume &
        / (REAL(nmols(is,box_in)))) * DEXP(-beta(box_in) * CP_energy)

    END IF









     !P_forward = P_forward * REAL(nmols(is,box_out),DP) / REAL(nmols_box(box_out),DP)
     !P_reverse = P_reverse * REAL(nmols(is,box_in)+1,DP) / REAL(nmols_box(box_in)+1,DP)

  !ELSE
    !

  !END IF

  ! Increment counters
  ntrials(is,box_out)%deletion = ntrials(is,box_out)%deletion + 1
  ntrials(is,box)%insertion = ntrials(is,box)%insertion + 1
  ntrials(is,box)%cpcalc = ntrials(is,box)%cpcalc + 1








END SUBROUTINE IDENTITY_EXCHANGE



