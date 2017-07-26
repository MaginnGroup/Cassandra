
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

  !*****************************************************************************
  ! Step 1) Select a box
  !*****************************************************************************
  ! If needed, choose a box based on its total mol fraction
  IF(nbr_boxes .GT. 1) THEN

    DO box = 1, nbr_boxes
       x_box(box) = REAL(nmols_box(box),DP)/REAL(nmols_tot,DP)
       IF (box > 1 ) THEN
          x_box(box) = x_box(box) + x_box(box-1)
       END IF
    END DO

    randno = rranf()
    DO box = 1, nbr_boxes
       IF ( randno <= x_box(box)) EXIT
    END DO
  ELSE
    box = 1
  END IF

  !*****************************************************************************
  ! Step 2) Select a species 'is':
  !      -> Right now just done randomly!
  !*****************************************************************************
  ! Sum the number of swappable molecules total & per box
  nmols_tot = 0 ! sum over species, box
  DO is = 1, nspecies
    ! Only count swappable species
    IF ( species_list(is)%int_insert /= int_noinsert ) THEN
      nmols_tot = nmols_tot + nmols(is,box)
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
  ! Step 3) Select a molecule 'alive' from species 'is' with uniform probability
  !*****************************************************************************
  ! pick a molecule INDEX at random
  im = INT(rranf() * nmols(is,box)) + 1

  ! Obtain the LOCATE of this molecule
  i_alive = locate(im, is, box)

  !*****************************************************************************
  ! Step 4) Select a species 'js':
  !      -> Right now just done randomly!
  !*****************************************************************************

  js = INT(rranf() * (nspecies - 1)) + 1

  !avoid choosing same species as before
  IF (js >= is)
    js = js + 1
  END IF

  !*****************************************************************************
  ! Step 5) Select a molecule 'alive' from species 'js' with uniform probability
  !*****************************************************************************
  jm = INT(rranf() * nmols(js,box)) + 1

  j_alive = locate(jm, js, box)

  !*****************************************************************************
  ! Step 6) Calculate the change in box_in's potential energy from inserting
  !
  !*****************************************************************************
  !update trials
  tot_trials(box) = tot_trials(box) + 1
  !TODO: UPDATE OTHER TRIALS TOO

  ! obtain the energy of the molecule before the move.  Note that due to
  ! this move, the interatomic energies such as vdw and electrostatics will
  ! change. Also the ewald_reciprocal energy will change but there will
  ! be no change in intramolecular energies.
  IF (l_pair_nrg) THEN
    CALL Store_Molecule_Pair_Interaction_Arrays(i_alive,is,box,E_vdw,E_qq)
    CALL Store_Molecule_Pair_Interaction_Arrays(j_alive,js,box,E_vdw,E_qq)
  ELSE
    CALL Compute_Molecule_Nonbond_Inter_Energy(i_alive,is,E_vdw,E_qq,inter_overlap_i)
    if (inter_overlap_i)
        EXIT
    END IF
    CALL Compute_Molecule_Nonbond_Inter_Energy(j_alive,js,E_vdw,E_qq,inter_overlap_j)
  END IF

  IF (inter_overlap_i .OR. inter_overlap_j)  THEN
     err_msg = ""

     err_msg(1) = "Attempted to move molecule " // TRIM(Int_To_String(i_alive)) // &
                  " of species " // TRIM(Int_To_String(is)) // "NOT RIGHT ABOUT SPECIS RIGHT NOW"
     IF (nbr_boxes > 1) err_msg(1) = err_msg(1) // " in box " // TRIM(Int_To_String(box))
     err_msg(2) = "but the molecule energy is too high"
     IF (start_type(ibox) == "make_config" ) THEN
        err_msg(3) = "Try increasing Rcutoff_Low, increasing the box size, or "
        err_msg(4) = "decreasing the initial number of molecules"
     END IF
     CALL Clean_Abort(err_msg, "Translate")
  END IF

  CALL Save_Old_Cartesian_Coordinates(i_alive,is)
  CALL Save_Old_Cartesian_Coordinates(j_alive,js)

  !******************************Actual*Switch*********************************!
  !Save COMs and COM to atom differences
  xcom_i = molecule_list(i_alive,is)%xcom
  ycom_i = molecule_list(i_alive,is)%ycom
  zcom_i = molecule_list(i_alive,is)%zcom

  dx_xcom_i = xcom_i - atom_list(:,i_alive,is)%rxp
  dy_ycom_i = ycom_i - atom_list(:,i_alive,is)%ryp
  dz_zcom_i = zcom_i - atom_list(:,i_alive,is)%rzp

  dx_xcom_j =  molecule_list(j_alive,js)%xcom - atom_list(:,j_alive,js)%rxp
  dy_ycom_j =  molecule_list(j_alive,js)%xcom - atom_list(:,j_alive,js)%ryp
  dz_zcom_j =  molecule_list(j_alive,js)%xcom - atom_list(:,j_alive,js)%rzp

  !Switch first molecule
  molecule_list(i_alive,is)%xcom = molecule_list(j_alive,js)%xcom
  molecule_list(i_alive,is)%ycom = molecule_list(j_alive,js)%ycom
  molecule_list(i_alive,is)%zcom = molecule_list(j_alive,js)%zcom
  atom_list(:,i_alive,is)%rxp = atom_list(:,i_alive,is)%com - dx_xcom_i
  atom_list(:,i_alive,is)%ryp = atom_list(:,i_alive,is)%com - dy_ycom_i
  atom_list(:,i_alive,is)%rzp = atom_list(:,i_alive,is)%com - dz_zcom_i

  !Switch second molecule
  molecule_list(j_alive,js)%xcom = xcom_i
  molecule_list(j_alive,js)%ycom = ycom_i
  molecule_list(j_alive,js)%zcom = zcom_i
  atom_list(:,i_alive,is)%rxp = xcom_i - dx_xcom_j
  atom_list(:,i_alive,is)%ryp = ycom_i - dy_ycom_j
  atom_list(:,i_alive,is)%rzp = zcom_i - dz_zcom_j

  !****************************************************************************!
  !*****************************************************************************
  ! Step 7) Calculate the change in box_in's potential energy from inserting
  !
  !*****************************************************************************

  CALL Fold_Molecule(i_alive, is, box)
  CALL Fold_Molecule(j_alive, js, box)

  CALL Compute_Molecule_Nonbond_Inter_Energy(i_alive,is,E_vdw,E_qq,inter_overlap_i)
  CALL Compute_Molecule_Nonbond_Inter_Energy(j_alive,js,E_vdw,E_qq,inter_overlap_j)

  IF (inter_overlap_i .OR. inter_overlap_j)
    CALL Revert_Old_Cartesian_Coordinates(i_alive, is)

    IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(i_alive,is,box)

    IF (verbose_log) THEN
      WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,A9)') &
            i_mcstep, 'identity switch' , i_alive, is, box, accept, 'overlap'
    END IF
  ELSE !no overlap
      dE = 0.0_DP

     !get help here!
     IF ((int_charge_sum_style(box) == charge_ewald) .AND. (has_charge(is))) THEN

        ALLOCATE(cos_mol_old(nvecs(box)),sin_mol_old(nvecs(box)))
        CALL Get_Position_Alive(i_alive,is,position)

        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        cos_mol_old(:) = cos_mol(1:nvecs(box),position)
        sin_mol_old(:) = sin_mol(1:nvecs(box),position)
        !$OMP END PARALLEL WORKSHARE

        !TODO: int_translation?
        CALL Update_System_Ewald_Reciprocal_Energy(i_alive,is,box,int_translation,E_reciprocal_move)
        dE = E_reciprocal_move - energy(box)%reciprocal + dE
     END IF

     IF ((int_charge_sum_style(box) == charge_ewald) .AND. (has_charge(js))) THEN
        ALLOCATE(cos_mol_old(nvecs(box)),sin_mol_old(nvecs(box)))
        CALL Get_Position_Alive(j_alive,js,position)

        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        cos_mol_old(:) = cos_mol(1:nvecs(box),position)
        sin_mol_old(:) = sin_mol(1:nvecs(box),position)
        !$OMP END PARALLEL WORKSHARE

        !TODO: int_translation?
        CALL Update_System_Ewald_Reciprocal_Energy(j_alive,js,box,int_translation,E_reciprocal_move)
        dE = E_reciprocal_move - energy(box)%reciprocal + dE
     END IF

     !compute difference in old and new energy
     dE = dE + (E_vdw_move - E_vdw) + (E_qq_move - E_qq)

     !TODO, find out what this does!
     IF (int_sim_type == sim_nvt_min) THEN
        IF (dE  <= 0.0_DP) THEN
           accept = .TRUE.
        END IF

     ELSE
         ln_pacc = beta(box) * dE
         accept = accept_or_reject(ln_pacc)
     END IF

     IF (accept) THEN
        !accept the move and update global energies
        energy(box)%total = energy(box)%total + dE
        energy(box)%inter = energy(box)%inter + dE
        energy(box)%inter_vdw = energy(box)%inter_vdw + E_vdw_move - E_vdw
        energy(box)%inter_q   = energy(box)%inter_q   + E_qq_move  - E_qq

        !TODO: HELP WITH EWALD
        IF(int_charge_sum_style(box) == charge_ewald .AND. has_charge(is)) THEN
           energy(box)%reciprocal = E_reciprocal_move
        END IF

        !TODO: update success counter
        !nsuccess(is,box)%displacement = nsuccess(is,box)%displacement + 1
        !nequil_success(is,ibox)%displacement = nequil_success(is,ibox)%displacement + 1

        IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
        IF (ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
        IF (ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)
     ELSE
        ! Revert to the old coordinates of atoms and com of the molecule
        CALL Revert_Old_Cartesian_Coordinates(lm,is)

        IF ((int_charge_sum_style(ibox) == charge_ewald) .AND. (has_charge(is))) THEN
           ! Also reset the old cos_sum and sin_sum for reciprocal space vectors. Note
           ! that old vectors were set while difference in ewald reciprocal energy was computed.
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
           cos_sum(:,ibox) = cos_sum_old(:,ibox)
           sin_sum(:,ibox) = sin_sum_old(:,ibox)
           cos_mol(1:nvecs(ibox),position) =cos_mol_old(:)
           sin_mol(1:nvecs(ibox),position) =sin_mol_old(:)
           !$OMP END PARALLEL WORKSHARE
           DEALLOCATE(cos_mol_old,sin_mol_old)
        END IF
        IF (l_pair_nrg)  CALL Reset_Molecule_Pair_Interaction_Arrays(lm,is,ibox)
     ENDIF

     IF (verbose_log) THEN
       WRITE(logunit,'(X,I15,X,A10,X,X,I5,X,I3,X,I3,X,L8,X,9X,X,F9.3)') &
             i_mcstep, 'identity switch' , i_alive, is, js, box, accept, ln_pacc
     END IF
  END IF !end of no overlap case

  !removed logging stuff, take another look!
END SUBROUTINE IDENTITY_EXCHANGE

