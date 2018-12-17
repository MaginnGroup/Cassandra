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

SUBROUTINE Identity_Switch

  !*****************************************************************************
  !
  ! This subroutine performs
  !
  !*****************************************************************************

  USE Global_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE IO_Utilities
  USE Fragment_Growth
  USE Pair_Nrg_Routines
  USE Type_Definitions
  USE Rotation_Routines

  IMPLICIT NONE

  !variables for selecting the box, and individual molecules
  INTEGER  :: box
  INTEGER  :: is, js            ! species index
  INTEGER  :: im, jm            ! molecule index
  INTEGER  :: i_alive, j_alive  ! molecule locate
  INTEGER  :: randint
  REAL(DP) :: x_box(nbr_boxes), x_box_i(nbr_boxes), x_box_j(nbr_boxes)
  REAL(DP) :: randno

  !variables for calculating change in energy from move
  REAL(DP) :: dE, E_vdw, E_qq, E_vdw_move, E_qq_move, E_reciprocal_move, E_reciprocal_move_i, E_reciprocal_move_j
  REAL(DP) :: E_qq_dum, E_vdw_dum
  REAL(DP) :: dE_i, dE_j, E_qq_i, E_qq_j, E_vdw_i, E_vdw_j
  REAL(DP) :: E_vdw_move_i, E_vdw_move_j, E_qq_move_i, E_qq_move_j
  INTEGER :: dum1, dum2, dum3, box_i, box_j
  LOGICAL :: inter_overlap
  REAL(DP), DIMENSION(:), ALLOCATABLE :: box_nrg_vdw_temp, box_nrg_qq_temp

  !variables for peforming the switch
  REAL(DP) :: xcom_i, ycom_i, zcom_i
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dx_xcom_i, dy_ycom_i, dz_zcom_i
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dx_xcom_j, dy_ycom_j, dz_zcom_j

  ! Initialize variables from insert
  INTEGER  :: nmols_tot ! number of molecules in the system
  INTEGER  :: nmols_box(nbr_boxes), trials
  INTEGER :: nmols_tot_i, nmols_tot_j

  !acceptance variables
  REAL(DP) :: ln_pacc, success_ratio_i, success_ratio_j
  LOGICAL :: accept_or_reject

 ! Variables added for l_pair_nrg and reciprocal k vector storage
  REAL(DP), ALLOCATABLE :: cos_mol_old_i(:), sin_mol_old_i(:), cos_mol_old_j(:), sin_mol_old_j(:)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, TARGET :: cos_sum_old_idsw, sin_sum_old_idsw
  INTEGER :: position_i, position_j

  !rotation variables
  REAL(DP) :: P_bias_rot_i, P_bias_rot_j
  LOGICAL :: rot_overlap_i, rot_overlap_j
  REAL(DP) :: P_bias
  INTEGER :: rotations, im_in, jm_in, i, j


  E_vdw_move = 0.0_DP
  E_qq_move = 0.0_DP
  E_vdw = 0.0_DP
  E_qq = 0.0_DP
  E_reciprocal_move = 0.0_DP
  inter_overlap = .FALSE.
  accept = .FALSE.
  box = 1
  rotations = 5

  P_bias_rot_i = 1
  P_bias_rot_j = 1
  rot_overlap_i = .FALSE.
  rot_overlap_j = .FALSE.

  !*****************************************************************************
  ! Steps 1 and 2) Select species i and j
  !*****************************************************************************
  !If the User has specified a limited set of species to switch, then we choose among
  !those species at random.
  IF (.NOT. default_switch) THEN
    randint = INT(rranf() * 2 * num_groups) + 1

    !Determine species
    is = swap_list(randint, 1)
    js = swap_list(randint, 2)

  ELSE
     nmols_tot = 0 ! sum over species, box
     DO is = 1, nspecies
       ! Only count swappable species
       DO box = 1, nbr_boxes
          IF ( species_list(is)%int_insert /= int_noinsert ) THEN
            nmols_tot = nmols_tot + nmols(is,box)
          END IF
       END DO
     END DO

     ! check to make sure there are swappable molecules
     IF (nmols_tot < 2) THEN
        err_msg = ''
        err_msg(1) = 'No swappable molecules'
        CALL Clean_Abort(err_msg,'Identity_Switch_Move')
     END IF

     !check this
     !old version:
     is = INT(rranf() * nspecies) + 1

     js = INT(rranf() * (nspecies - 1)) + 1

     !avoid choosing same species as before
     IF (js >= is) THEN
        js = js + 1
     END IF
  END IF

  WRITE(*,*) "species i is:"
  WRITE(*,*) is
  WRITE(*,*) "species j is:"
  WRITE(*,*) js

  !*****************************************************************************
  ! Step 3 and 4) Select a box for species i and a box for species j
  !*****************************************************************************

  nmols_tot_i = 0
  nmols_tot_j = 0
  DO box = 1, nbr_boxes
     nmols_tot_i = nmols_tot_i + nmols(is,box)
     nmols_tot_j = nmols_tot_j + nmols(js,box)
  END DO

  IF(nbr_boxes .GT. 1) THEN
     DO box = 1, nbr_boxes
        x_box_i(box) = REAL(nmols(is, box),DP)/REAL(nmols_tot_i,DP)
        x_box_j(box) = REAL(nmols(js, box),DP)/REAL(nmols_tot_j,DP)
        IF (box > 1 ) THEN
           x_box_i(box) = x_box_i(box) + x_box_i(box-1)
           x_box_j(box) = x_box_j(box) + x_box_j(box-1)
        END IF
    END DO

    randno = rranf()
    DO box_i = 1, nbr_boxes
          IF ( randno <= x_box_i(box_i)) EXIT
    END DO

    WRITE (*,*) x_box_i(1)
    WRITE (*,*) x_box_i(2)
    randno = rranf()
    DO box_j = 1, nbr_boxes
          IF ( randno <= x_box_j(box_j)) EXIT
    END DO
    WRITE (*,*) x_box_j(1)
    WRITE (*,*) x_box_j(2)

  ELSE
     box_i = 1
     box_j = 1
  END IF

  WRITE (*,*) "box_i is"
  WRITE (*,*) box_i
  WRITE (*,*) "box_j is"
  WRITE (*,*) box_j

  !*****************************************************************************
  ! Step 5) Select a molecule 'alive' from species 'is' with uniform probability
  !*****************************************************************************
  ! pick a molecule INDEX at random
  im = INT(rranf() * nmols(is,box_i)) + 1

  ! Obtain the LOCATE of this molecule
  i_alive = locate(im, is, box_i)


  !*****************************************************************************
  ! Step 6) Select a molecule 'alive' from species 'js' with uniform probability
  !*****************************************************************************
  jm = INT(rranf() * nmols(js,box_j)) + 1

  j_alive = locate(jm, js, box_j)

  WRITE (*,*) "FINISHED STEP 6"

  WRITE (*,*) "im is:"
  WRITE (*,*) im

  WRITE (*,*) "jm is:"
  WRITE (*,*) jm

  !*****************************************************************************
  ! Step 7) Calculate initial energies of each box
  !
  !*****************************************************************************
  !update trials
  WRITE(*,*) "Update Trials"
  tot_trials(box_i) = tot_trials(box_i) + 1
  tot_trials(box_j) = tot_trials(box_j) + 1
  ntrials(is,box_i)%switch = ntrials(is,box_i)%switch + 1
  ntrials(js,box_j)%switch = ntrials(js,box_j)%switch + 1
  WRITE(*,*) "End Update Trials"


  ! obtain the energy of the molecule before the move.  Note that due to
  ! this move, the interatomic energies such as vdw and electrostatics will
  ! change. Also the ewald_reciprocal energy will change but there will
  ! be no change in intramolecular energies.

  WRITE (*,*) "Start big If statement"
  IF (box_i .EQ. box_j) THEN
     WRITE (*,*) "a"
     IF (l_pair_nrg) THEN
       WRITE (*,*) "b"
       ALLOCATE(box_nrg_vdw_temp(2), box_nrg_qq_temp(2))
       WRITE (*,*) "b1"
       CALL Store_Molecule_Pair_Interaction_Arrays(dum1, dum2, dum3, E_vdw_dum, E_qq_dum, 2, &
            (/i_alive, j_alive/), (/is, js/), (/box_i, box_j/), box_nrg_vdw_temp, box_nrg_qq_temp)
       WRITE (*,*) "b2"
       E_vdw = box_nrg_vdw_temp(1)
       E_qq = box_nrg_vdw_temp(1)
       DEALLOCATE(box_nrg_vdw_temp, box_nrg_qq_temp)
       WRITE (*,*) "b3"
     WRITE (*,*) "c"
     ELSE
        WRITE (*,*) "d"
        WRITE (*,*) "Box for molecule i"
        WRITE (*,*) molecule_list(i_alive, is)%which_box

        WRITE (*,*) "Box for molecule j"
        WRITE (*,*) molecule_list(j_alive, js)%which_box


        CALL Compute_MoleculeCollection_Nonbond_Inter_Energy(2, (/i_alive, j_alive/), (/is, js/), &
             E_vdw, E_qq, inter_overlap)
        IF (inter_overlap) THEN
            WRITE(*,*) "OVERLAP ERROR!"
        END IF
        WRITE (*,*) "e"
     END IF
  !Not in the same box
  ELSE
     IF (l_pair_nrg) THEN
       ALLOCATE(box_nrg_vdw_temp(2), box_nrg_qq_temp(2))
       CALL Store_Molecule_Pair_Interaction_Arrays(dum1, dum2, dum3, E_vdw_dum, E_qq_dum, 2, &
            (/i_alive, j_alive/), (/is, js/), (/box_i, box_j/), box_nrg_vdw_temp, box_nrg_qq_temp)
       E_vdw_i = box_nrg_vdw_temp(1)
       E_vdw_j = box_nrg_vdw_temp(2)
       E_qq_i = box_nrg_vdw_temp(1)
       E_qq_j = box_nrg_vdw_temp(2)
       DEALLOCATE(box_nrg_vdw_temp, box_nrg_qq_temp)
     ELSE
        CALL Compute_Molecule_Nonbond_Inter_Energy(i_alive,is,E_vdw_i,E_qq_i,inter_overlap)
        IF (inter_overlap) THEN
            WRITE(*,*) "OVERLAP ERROR!"
        END IF
        IF (.NOT. inter_overlap) THEN
           CALL Compute_Molecule_Nonbond_Inter_Energy(j_alive,js,E_vdw_j,E_qq_j,inter_overlap)
            IF (inter_overlap) THEN
                WRITE(*,*) "OVERLAP ERROR!"
            END IF
        END IF
     END IF
  END IF
  WRITE (*,*) "End big if statement"

  IF (inter_overlap)  THEN
     err_msg = ""

     err_msg(1) = "Attempted to swap molecule " // TRIM(Int_To_String(im)) // &
                  " of species " // TRIM(Int_To_String(is)) // &
                  " and molecule " // TRIM(Int_To_String(jm)) // &
                  " of species " // TRIM(Int_To_String(js))
     IF (nbr_boxes > 1) err_msg(1) = err_msg(1) // " in box " // TRIM(Int_To_String(box))
     err_msg(2) = "but the molecule energy is too high"
     IF (start_type(box) == "make_config" ) THEN
        err_msg(3) = "Try increasing Rcutoff_Low, increasing the box size, or "
        err_msg(4) = "decreasing the initial number of molecules"
     END IF
     CALL Clean_Abort(err_msg, "Identity Switch")
  END IF

  CALL Save_Old_Cartesian_Coordinates(i_alive,is)
  CALL Save_Old_Cartesian_Coordinates(j_alive,js)

  !*****************************************************************************
  ! Step 8) Switch the two molecules by COM translation
  !
  !*****************************************************************************
  xcom_i = molecule_list(i_alive,is)%xcom
  ycom_i = molecule_list(i_alive,is)%ycom
  zcom_i = molecule_list(i_alive,is)%zcom

  WRITE(*,*) "Before swap, i has COM(x,y,z)"
  WRITE(*,*) xcom_i
  WRITE(*,*) ycom_i
  WRITE(*,*) zcom_i

  WRITE(*,*) "Individual x,y,z of atoms of i:"

  !For 3 atom example
  DO i = 1, 3
        WRITE(*,*) i
        WRITE(*,*) "x"
        WRITE(*,*) atom_list(i,i_alive,is)%rxp
        WRITE(*,*) "y"
        WRITE(*,*) atom_list(i,i_alive,is)%ryp
        WRITE(*,*) "z"
        WRITE(*,*) atom_list(i,i_alive,is)%rzp
  END DO

  WRITE(*,*) "Before swap, molecule j has COM(x,y,z)"
  WRITE(*,*) molecule_list(j_alive,js)%xcom
  WRITE(*,*)  molecule_list(j_alive,js)%ycom
  WRITE(*,*) molecule_list(j_alive,js)%zcom

  WRITE(*,*) "Individual x,y,z of atoms of j:"
  !For 3 atom example
  DO j = 1, 3
        WRITE(*,*) j
        WRITE(*,*) "x"
        WRITE(*,*) atom_list(j,j_alive,js)%rxp
        WRITE(*,*) "y"
        WRITE(*,*) atom_list(j,j_alive,js)%ryp
        WRITE(*,*) "z"
        WRITE(*,*) atom_list(j,j_alive,js)%rzp
  END DO

  ALLOCATE(dx_xcom_i(natoms(is)), dy_ycom_i(natoms(is)), dz_zcom_i(natoms(is)))
  ALLOCATE(dx_xcom_j(natoms(js)), dy_ycom_j(natoms(js)), dz_zcom_j(natoms(js)))
  dx_xcom_i = xcom_i - atom_list(:,i_alive,is)%rxp
  dy_ycom_i = ycom_i - atom_list(:,i_alive,is)%ryp
  dz_zcom_i = zcom_i - atom_list(:,i_alive,is)%rzp

  dx_xcom_j =  molecule_list(j_alive,js)%xcom - atom_list(:,j_alive,js)%rxp
  dy_ycom_j =  molecule_list(j_alive,js)%ycom - atom_list(:,j_alive,js)%ryp
  dz_zcom_j =  molecule_list(j_alive,js)%zcom - atom_list(:,j_alive,js)%rzp

  !Switch first molecule
  molecule_list(i_alive,is)%xcom = molecule_list(j_alive,js)%xcom
  molecule_list(i_alive,is)%ycom = molecule_list(j_alive,js)%ycom
  molecule_list(i_alive,is)%zcom = molecule_list(j_alive,js)%zcom
  atom_list(:,i_alive,is)%rxp = molecule_list(i_alive,is)%xcom - dx_xcom_i
  atom_list(:,i_alive,is)%ryp = molecule_list(i_alive,is)%ycom - dy_ycom_i
  atom_list(:,i_alive,is)%rzp = molecule_list(i_alive,is)%zcom - dz_zcom_i

  !Switch second molecule
  molecule_list(j_alive,js)%xcom = xcom_i
  molecule_list(j_alive,js)%ycom = ycom_i
  molecule_list(j_alive,js)%zcom = zcom_i
  atom_list(:,j_alive,js)%rxp = xcom_i - dx_xcom_j
  atom_list(:,j_alive,js)%ryp = ycom_i - dy_ycom_j
  atom_list(:,j_alive,js)%rzp = zcom_i - dz_zcom_j
  DEALLOCATE(dx_xcom_i, dy_ycom_i, dz_zcom_i)
  DEALLOCATE(dx_xcom_j, dy_ycom_j, dz_zcom_j)

  !change box identity
  molecule_list(i_alive,is)%which_box = box_j
  molecule_list(j_alive,js)%which_box = box_i

  nmols(is,box_i) = nmols(is,box_i) - 1
  nmols(js,box_j) = nmols(js,box_j) - 1
  nmols(is,box_j) = nmols(is,box_j) + 1
  nmols(js,box_i) = nmols(js,box_i) + 1

!  molecule_list(alive,is)%which_box = box_in
!  nmols(is,box_in) = nmols(is,box_in) + 1
!  im_in = nmols(is,box_in) ! INDEX of alive in box_in
!  locate(im_in,is,box_in) = alive ! link INDEX to LOCATE

  WRITE (*,*) "FOLDING!"
  CALL Fold_Molecule(i_alive,is,box_j)
  CALL Fold_Molecule(j_alive,js,box_i)
  WRITE (*,*) "DONE FOLDING!"

  !im_in = nmols(is, box_j)
  !locate(im_in,is,box_j) = i_alive

  !jm_in = nmols(js, box_i)
  !locate(jm_in,js,box_i) = j_alive

  !!Got rid of anything relating to locate
  !locate(im,is,box_j) = i_alive
  !locate(jm,js,box_i) = j_alive

  WRITE(*,*) "After swap, molecule i has COM(x,y,z)"
  WRITE(*,*) molecule_list(i_alive,is)%xcom
  WRITE(*,*)  molecule_list(i_alive,is)%ycom
  WRITE(*,*) molecule_list(i_alive,is)%zcom

  WRITE(*,*) "Individual x,y,z of atoms of i:"
  !For 3 atom example
  DO i = 1, 3
        WRITE(*,*) i
        WRITE(*,*) "x"
        WRITE(*,*) atom_list(i,i_alive,is)%rxp
        WRITE(*,*) "y"
        WRITE(*,*) atom_list(i,i_alive,is)%ryp
        WRITE(*,*) "z"
        WRITE(*,*) atom_list(i,i_alive,is)%rzp
  END DO

  WRITE(*,*) "After swap, molecule j has COM(x,y,z)"
  WRITE(*,*) molecule_list(j_alive,js)%xcom
  WRITE(*,*)  molecule_list(j_alive,js)%ycom
  WRITE(*,*) molecule_list(j_alive,js)%zcom

  WRITE(*,*) "Individual x,y,z of atoms of j:"
  !For 3 atom example
  DO j = 1, 3
        WRITE(*,*) j
        WRITE(*,*) "x"
        WRITE(*,*) atom_list(j,j_alive,js)%rxp
        WRITE(*,*) "y"
        WRITE(*,*) atom_list(j,j_alive,js)%ryp
        WRITE(*,*) "z"
        WRITE(*,*) atom_list(j,j_alive,js)%rzp
  END DO

  WRITE (*,*) "STEP 8"
  !*****************************************************************************
  ! Step 9) Rotate the molecules if desired
  !
  !*****************************************************************************
  !TODO re-add back in rotations
  rotations = 0
  IF (rotations .NE. 0) THEN
        CALL Bias_Rotate(i_alive,is,box_j,rotations,P_bias_rot_i,rot_overlap_i)
        P_bias = P_bias * P_bias_rot_i

        CALL Bias_Rotate(j_alive,js,box_i,rotations,P_bias_rot_i,rot_overlap_j)
        P_bias = P_bias * P_bias_rot_j
  END IF

  WRITE (*,*) "STEP 9"

  !*****************************************************************************
  ! Step 10) Calculate the energy change of each box after rotation
  !
  !*****************************************************************************
  !write(*,*) "Computing Nonbond Inter Energy WITH COLLECTION"
  IF (box_i .EQ. box_j) THEN
     WRITE(*,*) "1"
     WRITE(*,*) "1a"
     CALL Compute_MoleculeCollection_Nonbond_Inter_Energy(2, (/i_alive, j_alive/), (/is, js/), &
         E_vdw_move, E_qq_move, inter_overlap)
        IF (inter_overlap) THEN
            WRITE(*,*) "OVERLAP ERROR!"
        END IF
     WRITE(*,*) "1b"
  ELSE
     WRITE(*,*) "2"
     CALL Compute_Molecule_Nonbond_Inter_Energy(i_alive,is,E_vdw_move_i,E_qq_move_i,inter_overlap)
        IF (inter_overlap) THEN
            WRITE(*,*) "OVERLAP ERROR!"
        END IF
     IF (inter_overlap) THEN
        WRITE(*,*) "3"
        CALL Compute_Molecule_Nonbond_Inter_Energy(j_alive,js,E_vdw_move_j,E_qq_move_j,inter_overlap)
        IF (inter_overlap) THEN
            WRITE(*,*) "OVERLAP ERROR!"
        END IF
     END IF
  END IF
  WRITE(*,*) "4"

  IF (inter_overlap .OR. (rot_overlap_i .OR. rot_overlap_j)) THEN
    !write (*, *) "Overlap found!"
            WRITE(*,*) "OVERLAP ERROR!"
    CALL Revert_Old_Cartesian_Coordinates(i_alive, is)
    CALL Revert_Old_Cartesian_Coordinates(j_alive, js)

    !revert box identity
    molecule_list(i_alive,is)%which_box = box_i
    molecule_list(j_alive,js)%which_box = box_j
    !locate(im_in,is,box_j) = 0
    !locate(jm_in,js,box_i) = 0

    nmols(is,box_i) = nmols(is,box_i) + 1
    nmols(js,box_j) = nmols(js,box_j) + 1
    nmols(is,box_j) = nmols(is,box_j) - 1
    nmols(js,box_i) = nmols(js,box_i) - 1

    IF (l_pair_nrg) THEN
       CALL Reset_Molecule_Pair_Interaction_Arrays(i_alive,is,box_i)
       CALL Reset_Molecule_Pair_Interaction_Arrays(j_alive,js,box_j)
    END IF

    IF (verbose_log) THEN
      WRITE(logunit,'(X,I9,X,A16,X,I5,X,I3,X,I3,X,L8,X,9X,X,A9)') &
            i_mcstep, 'identity switch' , i_alive, is, box, accept, 'overlap'
    END IF
  ELSE !no overlap
     WRITE(*,*) "5"
     dE = 0.0_DP
     dE_i = 0.0_DP
     dE_j = 0.0_DP

     IF ((int_charge_sum_style(box) == charge_ewald) .AND. (has_charge(is) .OR. has_charge(js))) THEN
        !TODO:Eventually rewrite Update_System_Ewald_Reciprocal_Energy!
        !Update_System_Ewald_Reciprocal_Energy will do this, but it might override it!
        IF (box_i .EQ. box_j) THEN
           box = box_i

           ALLOCATE(cos_sum_old_idsw(nvecs(box),nbr_boxes), sin_sum_old_idsw(nvecs(box),nbr_boxes))
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
           cos_sum_old_idsw(1:nvecs(box),box) = cos_sum(1:nvecs(box),box)
           sin_sum_old_idsw(1:nvecs(box),box) = sin_sum(1:nvecs(box),box)
           !$OMP END PARALLEL WORKSHARE

           IF (has_charge(is)) THEN
              ALLOCATE(cos_mol_old_i(nvecs(box)), sin_mol_old_i(nvecs(box)))
              CALL Get_Position_Alive(i_alive, is, position_i)

              !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
              cos_mol_old_i(:) = cos_mol(1:nvecs(box),position_i)
              sin_mol_old_i(:) = sin_mol(1:nvecs(box),position_i)
              !$OMP END PARALLEL WORKSHARE

              CALL Update_System_Ewald_Reciprocal_Energy(i_alive,is,box,int_translation,E_reciprocal_move)
              dE = E_reciprocal_move + dE
           END IF
           IF (has_charge(js)) THEN
              ALLOCATE(cos_mol_old_j(nvecs(box)), sin_mol_old_j(nvecs(box)))
              CALL Get_Position_Alive(j_alive, js, position_j)

              !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
              cos_mol_old_j(:) = cos_mol(1:nvecs(box),position_j)
              sin_mol_old_j(:) = sin_mol(1:nvecs(box),position_j)
              !$OMP END PARALLEL WORKSHARE

              CALL Update_System_Ewald_Reciprocal_Energy(j_alive,js,box,int_translation,E_reciprocal_move)
              dE = E_reciprocal_move + dE
           END IF
           dE = dE - energy(box)%reciprocal
        !Can just use Update_System_Ewald_Reciprocal_Energy
        ELSE
           CALL Update_System_Ewald_Reciprocal_Energy(i_alive,is,box_j, &
              int_insertion,E_reciprocal_move_j)
           dE_j = dE_j + E_reciprocal_move_j - energy(box_j)%reciprocal
           CALL Update_System_Ewald_Reciprocal_Energy(j_alive,js,box_i, &
              int_insertion,E_reciprocal_move_i)
           dE_i = dE_i + E_reciprocal_move_i - energy(box_i)%reciprocal
        END IF

     END IF
     WRITE(*,*) "6"
     !write (*, *) "Done with big if statement"
     !compute difference in old and new energy
     IF (box_i .EQ. box_j) THEN
        dE = dE + (E_vdw_move - E_vdw) + (E_qq_move - E_qq)
        IF (int_sim_type == sim_nvt_min) THEN
           IF (dE  <= 0.0_DP) THEN
              accept = .TRUE.
           END IF
        ELSE
           ln_pacc = beta(box_i) * dE
           ln_pacc = ln_pacc + DLOG(P_bias)
           accept = accept_or_reject(ln_pacc)
        END IF
     ELSE
        dE_i = dE_i + (E_vdw_move_i - E_vdw_i) + (E_qq_move_i - E_qq_i)
        dE_j = dE_j + (E_vdw_move_j - E_vdw_j) + (E_qq_move_j - E_qq_j)

        IF (int_sim_type == sim_nvt_min) THEN
           IF ((dE_i + dE_j)  <= 0.0_DP) THEN
              accept = .TRUE.
           END IF
        ELSE
           ln_pacc = beta(box_i) * dE_i + beta(box_j) * dE_j
           ln_pacc = ln_pacc + DLOG(P_bias)
           accept = accept_or_reject(ln_pacc)
        END IF
     END IF
     WRITE(*,*) "7"

     !KEEP GOING!
     IF (accept) THEN
        !accept the move and update global energies
        !write (*,*) energy(box)%total
        IF (box_i .EQ. box_j) THEN
           energy(box_i)%total = energy(box_i)%total + dE
           energy(box_i)%inter = energy(box_i)%inter + dE
           energy(box_i)%inter_vdw = energy(box_i)%inter_vdw + E_vdw_move - E_vdw
           energy(box_i)%inter_q   = energy(box_i)%inter_q   + E_qq_move  - E_qq

           IF(int_charge_sum_style(box_i) == charge_ewald .AND. (has_charge(is) .OR. has_charge(js))) THEN
              energy(box_i)%reciprocal = E_reciprocal_move
           END IF

           nsuccess(is,box_i)%switch = nsuccess(is,box_i)%switch + 1
           nsuccess(js,box_i)%switch = nsuccess(js,box_i)%switch + 1
           nequil_success(is,box_i)%switch = nequil_success(is,box_i)%switch + 1
           nequil_success(js,box_i)%switch = nequil_success(js,box_i)%switch + 1

        ELSE
            energy(box_i)%total = energy(box_i)%total + dE_i
            energy(box_i)%inter = energy(box_i)%inter + dE_i
            energy(box_i)%inter_vdw = energy(box_i)%inter_vdw + E_vdw_move_i - E_vdw_i
            energy(box_i)%inter_q = energy(box_i)%inter_q + E_qq_move_i - E_qq_i


            energy(box_j)%total = energy(box_j)%total + dE_j
            energy(box_j)%inter = energy(box_j)%inter + dE_j
            energy(box_j)%inter_vdw = energy(box_j)%inter_vdw + E_vdw_move_j - E_vdw
            energy(box_j)%inter_q = energy(box_j)%inter_q + E_qq_move_j - E_qq_j

           IF(int_charge_sum_style(box_i) == charge_ewald .AND. (has_charge(is) .OR. has_charge(js))) THEN
              energy(box_i)%reciprocal = E_reciprocal_move_i
           END IF
           IF(int_charge_sum_style(box_j) == charge_ewald .AND. (has_charge(is) .OR. has_charge(js))) THEN
              energy(box_j)%reciprocal = E_reciprocal_move_j
           END IF

           nsuccess(is,box_i)%switch = nsuccess(is,box_i)%switch + 1
           nsuccess(js,box_j)%switch = nsuccess(js,box_j)%switch + 1
           nequil_success(is,box_i)%switch = nequil_success(is,box_i)%switch + 1
           nequil_success(js,box_j)%switch = nequil_success(js,box_j)%switch + 1

        END IF
        WRITE(*,*) "8"

        IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
        IF (ALLOCATED(cos_mol_old_i)) DEALLOCATE(cos_mol_old_i)
        IF (ALLOCATED(sin_mol_old_i)) DEALLOCATE(sin_mol_old_i)
        IF (ALLOCATED(sin_mol_old_j)) DEALLOCATE(sin_mol_old_j)
        IF (ALLOCATED(sin_mol_old_j)) DEALLOCATE(sin_mol_old_j)

     ELSE
        ! Revert to the old coordinates of atoms and com of the molecule
        CALL Revert_Old_Cartesian_Coordinates(i_alive,is)
        CALL Revert_Old_Cartesian_Coordinates(j_alive,js)

        !revert box identity
        molecule_list(i_alive,is)%which_box = box_i
        molecule_list(j_alive,js)%which_box = box_j

        nmols(is,box_i) = nmols(is,box_i) + 1
        nmols(js,box_j) = nmols(js,box_j) + 1
        nmols(is,box_j) = nmols(is,box_j) - 1
        nmols(js,box_i) = nmols(js,box_i) - 1

        !TODO: Rewrite this bit
        IF ((int_charge_sum_style(box) == charge_ewald) .AND. (has_charge(is) .OR. has_charge(js))) THEN
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
           cos_sum(1:nvecs(box), box) = cos_sum_old_idsw(1:nvecs(box),box)
           sin_sum(1:nvecs(box), box) = sin_sum_old_idsw(1:nvecs(box),box)
           !$OMP END PARALLEL WORKSHARE

           DEALLOCATE(cos_sum_old_idsw, sin_sum_old_idsw)
           IF (has_charge(is)) THEN
              !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
              cos_mol(1:nvecs(box),position_i) = cos_mol_old_i(:)
              sin_mol(1:nvecs(box),position_i) = sin_mol_old_i(:)
              !$OMP END PARALLEL WORKSHARE
              DEALLOCATE(cos_mol_old_i, sin_mol_old_i)
           END IF
           IF (has_charge(js)) THEN
              !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
              cos_mol(1:nvecs(box),position_j) = cos_mol_old_j(:)
              sin_mol(1:nvecs(box),position_j) = sin_mol_old_j(:)
              !$OMP END PARALLEL WORKSHARE
              DEALLOCATE(cos_mol_old_j, sin_mol_old_j)
           END IF
        END IF
        !END TODO
        IF (l_pair_nrg) THEN
            CALL Reset_Molecule_Pair_Interaction_Arrays(i_alive,is,box_i)
            CALL Reset_Molecule_Pair_Interaction_Arrays(j_alive,js,box_j)
        END IF
     END IF

     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A16,X,X,I5,X,I3,X,I3,X,L8,X,9X,X,F9.3)') &
             i_mcstep, 'identity switch' , i_alive, is, js, box, accept, ln_pacc
     END IF
  END IF

  !TODO: FIX Trials Calculation
  trials = ntrials(is, box_i)%switch + ntrials(is, box_j)%switch
  IF ( MOD(trials,nupdate) == 0 ) THEN
     IF ( int_run_type == run_equil ) THEN
        success_ratio_i = REAL(nequil_success(is,box_i)%switch + nequil_success(is,box_j)%switch,DP)/REAL(nupdate,DP)
        success_ratio_j = REAL(nequil_success(js,box_j)%switch + nequil_success(js,box_j)%switch,DP)/REAL(nupdate,DP)
     ELSE
        success_ratio_i = REAL(nsuccess(is,box_i)%switch+nsuccess(is,box_j)%switch,DP) / &
                          REAL(ntrials(is,box_i)%switch + REAL(ntrials(is,box_j)%switch) ,DP)
        success_ratio_j = REAL(nsuccess(js,box_i)%switch+nsuccess(js,box_j)%switch,DP) / &
                          REAL(ntrials(js,box_i)%switch + REAL(ntrials(js,box_j)%switch) ,DP)
     END IF

     WRITE(logunit,'(X,I9,X,A16,X,5X,X,I3,X,I3,X,F8.5)',ADVANCE='NO') &
           i_mcstep, 'identity switch', is, box_i, success_ratio_i

     WRITE(logunit,'(X,I9,X,A16,X,5X,X,I3,X,I3,X,F8.5)',ADVANCE='NO') &
           i_mcstep, 'identity switch', js, box_j, success_ratio_j
     !TODO: No need to change displacement like translation, but do we need to change rcut?
     WRITE(logunit,*)
  END IF

  WRITE (*,*) "FINISHED AN IS"

  CONTAINS
  !***************************************************************************
    SUBROUTINE Bias_Rotate(lm,is,ibox,rotations,P_bias,rot_overlap)
      INTEGER, INTENT(IN) :: lm, is, ibox, rotations
      REAL(DP), INTENT(OUT) :: P_bias
      LOGICAL, INTENT(OUT) :: rot_overlap
      LOGICAL :: overlap
      INTEGER :: ia, i_rot
      REAL(DP) :: nrg(rotations), nrg_rev(rotations*2), weight(rotations*2), weight_rev(rotations*2)
      REAL(DP) :: nrg_kBT, rand_no
      REAL(DP) :: P_bias_rev
      TYPE(Atom_Class) :: rtrial(natoms(is),rotations*2), rtrial_rev(natoms(is),rotations*2)

      REAL(DP) nrg_original, wt_original

      REAL(DP) :: E_inter_qq, E_inter_vdw

      ! Initialize
      P_bias = 1.0_DP
      P_bias_rev = 1.0_DP
      rot_overlap = .FALSE.
      !itrial = 0
      !ntrials = kappa_ins*kappa_rot
      nrg = 0.0_DP
      weight = 0.0_DP

      CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw, &
          E_inter_qq,overlap)

      nrg_original = E_inter_vdw + E_inter_qq
      WRITE (*,*) "Original Energy:"
      WRITE (*,*) nrg_original

      wt_original = DEXP(-1.0 * beta(ibox) * nrg_original)

      WRITE (*,*) wt_original

        DO i_rot = 1, rotations
          WRITE (*,*) "ROTATION:"
          WRITE (*,*) i_rot
          CALL Rotate_Molecule_Eulerian(lm,is)
          WRITE (*,*) "WRITING COORDINATES"
          CALL Write_Coords(ibox)

          !save rotation coordinates
          DO ia = 1, natoms(is)
            rtrial(ia,i_rot)%rxp = atom_list(ia,lm,is)%rxp
            rtrial(ia,i_rot)%ryp = atom_list(ia,lm,is)%ryp
            rtrial(ia,i_rot)%rzp = atom_list(ia,lm,is)%rzp
          END DO

          overlap = .FALSE.
          CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw, &
                  E_inter_qq,overlap)
          nrg(i_rot) = nrg(i_rot) + E_inter_vdw + E_inter_qq
          WRITE (*,*) "ENERGY:"
          WRITE (*,*) nrg(i_rot)

          IF (overlap) THEN
            weight(i_rot) = 0.0_DP
          ELSE
            nrg_kBT = beta(ibox) * nrg(i_rot)
            IF (nrg_kBT >= max_kBT) THEN
              weight(i_rot) = 0.0_DP
            ELSE
              weight(i_rot) = DEXP(-nrg_kBT)
            END IF
          END IF

          !WRITE (*,*) "ENERGY*beta:"
          !WRITE (*,*) nrg_kBT

          !WRITE (*,*) "exp(-ENERGY*beta)"
          !WRITE (*,*) weight(i_rot)

          IF (i_rot > 1) weight(i_rot) = weight(i_rot-1) + weight(i_rot)

          WRITE (*,*) "adjusted weight"
          WRITE (*,*) weight(i_rot)
        END DO

      ! Reject the move if the total weight is still zero
      IF (weight(rotations) == 0.0_DP) THEN
        rot_overlap = .TRUE.
        RETURN
      END IF

      ! Choose one from Golden sampling
      rand_no = rranf() * weight(rotations)

      DO i_rot = 1, rotations
        IF (rand_no < weight(i_rot)) EXIT
      END DO

      WRITE (*,*) "Random number"
      WRITE (*,*) rand_no

      WRITE (*,*) "Selected rotation:"
      WRITE (*,*) i_rot

      IF ( i_rot == rotations + 1 ) THEN
        ! None of the trials were picked. Could be due to the fact that all
        ! the trials had a very small cumulative weight
        rot_overlap = .TRUE.
        RETURN
      END IF

      ! Compute the weight of the selected trial coordinate
      ! save the weight of selected for the future reverse calculation
      IF (i_rot == 1) THEN
        P_bias = P_bias * weight(1) / weight(rotations)
        weight_rev(1) = weight(1)
      ELSE
        P_bias = P_bias * (weight(i_rot) - weight(i_rot-1)) / weight(rotations)
        weight_rev(1) = weight(i_rot) - weight(i_rot-1)
      END IF

      ! We chose the ith trial coordinate for placement. Store the ith trial
      ! coordinates in the atom_list array.
      DO ia = 1, natoms(is)
        atom_list(ia,lm,is)%rxp = rtrial(ia,i_rot)%rxp
        atom_list(ia,lm,is)%ryp = rtrial(ia,i_rot)%ryp
        atom_list(ia,lm,is)%rzp = rtrial(ia,i_rot)%rzp
      END DO

      WRITE (*,*) "Selected a rotation"
      CALL Write_Coords(ibox)


      !COM was not stored, so recalculate
      CALL Get_COM(lm,is)


      !Now we repeat all that junk in reverse to get the P_bias reverse
      !***************************************!

      !TODO: FIX THIS TO BE THE RIGHT FIRST ROTATION TRIAL
      !already stored the weight for the first 'rotation'
      weight_rev(1) = wt_original
        DO i_rot = 2, rotations
          WRITE (*,*) "Rotation, energy, weight"
          WRITE (*,*) i_rot
          CALL Rotate_Molecule_Eulerian(lm,is)
          !don't need to save rotation coordinates
          overlap = .FALSE.
          CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw, &
                  E_inter_qq,overlap)
          nrg_rev(i_rot) = nrg_rev(i_rot) + E_inter_vdw + E_inter_qq
          WRITE (*,*) nrg_rev(i_rot)

          IF (overlap) THEN
            weight_rev(i_rot) = 0.0_DP
          ELSE
            nrg_kBT = beta(ibox) * nrg_rev(i_rot)
            IF (nrg_kBT >= max_kBT) THEN
              weight_rev(i_rot) = 0.0_DP
            ELSE
              weight_rev(i_rot) = DEXP(-nrg_kBT)
            END IF
          END IF
          WRITE(*,*) weight_rev(i_rot)
          weight_rev(i_rot) = weight_rev(i_rot-1) + weight_rev(i_rot)
        END DO

      ! For the reverse move, we've already selected the first rotation
      P_bias_rev = P_bias_rev * weight(1) / weight(rotations)

      P_bias = P_bias / P_bias_rev
    END SUBROUTINE Bias_Rotate


END SUBROUTINE Identity_Switch

