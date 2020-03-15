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
!********************************************************************************
SUBROUTINE Identity_Switch

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

   !Variables for steps 1 and 2
   INTEGER :: randint
   INTEGER :: nmols_tot, nmols_tot_i, nmols_tot_j
   INTEGER :: is, js, ibox, box

   !Variables for steps 3 and 4
   REAL(DP) :: randno
   INTEGER :: box_i, box_j
   REAL(DP) :: x_box(nbr_boxes), x_box_i(nbr_boxes), x_box_j(nbr_boxes)

   !Variables for steps 5 and 6
   INTEGER :: im_i, lm_i
   INTEGER :: im_j, lm_j

   !Variables for step 7
   REAL(DP) :: E_vdw, E_qq, E_vdw_i, E_qq_i, E_vdw_j, E_qq_j
   REAL(DP), DIMENSION(:), ALLOCATABLE :: box_nrg_vdw_temp, box_nrg_qq_temp
   LOGICAL :: inter_overlap
   INTEGER :: dum1, dum2, dum3, position_i, position_j
   REAL(DP) :: E_qq_dum, E_vdw_dum

   !Variables for step 8
   INTEGER :: i,j,k, im_i_after, im_j_after, temp
   REAL(DP) :: xcom_i, ycom_i, zcom_i
   REAL(DP), DIMENSION(:), ALLOCATABLE :: dx_xcom_i, dy_ycom_i, dz_zcom_i
   REAL(DP), DIMENSION(:), ALLOCATABLE :: dx_xcom_j, dy_ycom_j, dz_zcom_j

   !Variables for step 10
   REAL(DP) :: E_vdw_move, E_qq_move, E_vdw_move_i, E_qq_move_i, E_vdw_move_j, E_qq_move_j
   REAL(DP) :: E_periodic_qq_i, E_periodic_qq_move_i, E_periodic_qq_j, E_periodic_qq_move_j
   REAL(DP) :: dE, dE_i, dE_j
   REAL(DP) :: E_bond_i, E_bond_j, E_angle_i, E_angle_j, E_dihed_i, E_dihed_j, E_improper_i, E_improper_j
   REAL(DP) :: E_intra_vdw_i, E_intra_qq_i, E_intra_vdw_j, E_intra_qq_j
   REAL(DP) :: dE_bond, dE_angle, dE_dihed, dE_improper, de_intra_vdw, dE_intra_qq
   REAL(DP) :: E_reciprocal_move, E_reciprocal_move_i, E_reciprocal_move_j
   REAL(DP) :: E_lrc, E_lrc_i, E_lrc_j
   REAL(DP) :: dE_lrc, dE_lrc_i, dE_lrc_j
   REAL(DP), ALLOCATABLE :: cos_mol_old_i(:), sin_mol_old_i(:), cos_mol_old_j(:), sin_mol_old_j(:)
   REAL(DP), DIMENSION(:,:), ALLOCATABLE, TARGET :: cos_sum_old_idsw, sin_sum_old_idsw


   !acceptance variables
   REAL(DP) :: ln_pacc, success_ratio_i, success_ratio_j
   LOGICAL :: accept_or_reject

   !rotation variables
   REAL(DP) :: P_bias_rot_i, P_bias_rot_j
   LOGICAL :: rot_overlap_i, rot_overlap_j
   REAL(DP) :: P_bias
   INTEGER :: im_in, jm_in

   E_vdw_move = 0.0_DP
   E_qq_move = 0.0_DP
   E_vdw = 0.0_DP
   E_qq = 0.0_DP
   E_reciprocal_move = 0.0_DP
   inter_overlap = .FALSE.
   accept = .FALSE.
   box = 1

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
       DO ibox = 1, nbr_boxes
          IF ( species_list(is)%int_insert /= int_noinsert ) THEN
            nmols_tot = nmols_tot + nmols(is,ibox)
          END IF
       END DO
     END DO

     ! check to make sure there are swappable molecules
     IF (nmols_tot < 2) THEN
        err_msg = ''
        err_msg(1) = 'No swappable molecules'
        CALL Clean_Abort(err_msg,'Identity_Switch_Move')
     END IF

     is = INT(rranf() * nspecies) + 1

     js = INT(rranf() * (nspecies - 1)) + 1

     !avoid choosing same species as before
     IF (js >= is) THEN
        js = js + 1
     END IF
   END IF

   !*****************************************************************************
   ! Step 3 and 4) Select a box for species i and a box for species j
   !*****************************************************************************

   nmols_tot_i = 0
   nmols_tot_j = 0
   DO ibox = 1, nbr_boxes
     nmols_tot_i = nmols_tot_i + nmols(is,ibox)
     nmols_tot_j = nmols_tot_j + nmols(js,ibox)
   END DO

   IF(nbr_boxes .GT. 1) THEN
     DO ibox = 1, nbr_boxes
        x_box_i(ibox) = REAL(nmols(is, ibox),DP)/REAL(nmols_tot_i,DP)
        x_box_j(ibox) = REAL(nmols(js, ibox),DP)/REAL(nmols_tot_j,DP)
        IF (ibox > 1 ) THEN
           x_box_i(ibox) = x_box_i(ibox) + x_box_i(ibox-1)
           x_box_j(ibox) = x_box_j(ibox) + x_box_j(ibox-1)
        END IF
    END DO

    randno = rranf()
    DO box_i = 1, nbr_boxes
          IF ( randno <= x_box_i(box_i)) EXIT
    END DO

    randno = rranf()
    DO box_j = 1, nbr_boxes
          IF ( randno <= x_box_j(box_j)) EXIT
    END DO

   ELSE
     box_i = 1
     box_j = 1
   END IF

   !*****************************************************************************
   ! Step 5) Select a molecule 'alive' from species 'is' with uniform probability
   !*****************************************************************************
   ! pick a molecule INDEX at random
   ! Index molecule of species i
   im_i = INT(rranf() * nmols(is,box_i)) + 1


   !Locate molecule of species i
   lm_i = locate(im_i, is, box_i)

   !*****************************************************************************
   ! Step 6) Select a molecule 'alive' from species 'js' with uniform probability
   !*****************************************************************************
   im_j = INT(rranf() * nmols(js,box_j)) + 1


   lm_j = locate(im_j, js, box_j)

   !WRITE (*,*) "FINISHED STEP 6"

   !WRITE (*,*) "im_i is:"
   !WRITE (*,*) im_i

   !WRITE (*,*) "im_j is:"
   !WRITE (*,*) im_j

   !*****************************************************************************
   ! Step 7) Calculate initial energies of each box
   !
   !*****************************************************************************
   tot_trials(box_i) = tot_trials(box_i) + 1
   tot_trials(box_j) = tot_trials(box_j) + 1
   ntrials(is,box_i)%switch = ntrials(is,box_i)%switch + 1
   ntrials(js,box_j)%switch = ntrials(js,box_j)%switch + 1

   !Same box
   IF (box_i .EQ. box_j) THEN
      IF (l_pair_nrg) THEN
         ALLOCATE(box_nrg_vdw_temp(2), box_nrg_qq_temp(2))
         CALL Store_Molecule_Pair_Interaction_Arrays(dum1, dum2, dum3, E_vdw_dum, E_qq_dum, 2, &
            (/lm_i, lm_j/), (/is, js/), (/box_i, box_j/), box_nrg_vdw_temp, box_nrg_qq_temp)
         E_vdw = box_nrg_vdw_temp(1)
         E_qq = box_nrg_vdw_temp(1)
         ! remove double counting from vdw and qq energies
         ! obtain the position of these molecules to reference vdw and qq pair energy arrays
         CALL Get_Position_Alive(lm_i, is, position_i)
         CALL Get_Position_Alive(lm_j, js, position_j)
         ! substract off the energy 
         E_vdw = E_vdw - pair_nrg_vdw(position_i,position_j)
         E_qq = E_qq - pair_nrg_qq(position_i,position_j)
         DEALLOCATE(box_nrg_vdw_temp, box_nrg_qq_temp)
      ELSE
         CALL Compute_MoleculeCollection_Nonbond_Inter_Energy(2, (/lm_i, lm_j/), (/is, js/), &
            E_vdw, E_qq, inter_overlap)
         IF (inter_overlap) THEN
         END IF
      END IF
   !Not in the same box
   ELSE
      IF (l_pair_nrg) THEN
         ALLOCATE(box_nrg_vdw_temp(2), box_nrg_qq_temp(2))
         CALL Store_Molecule_Pair_Interaction_Arrays(dum1, dum2, dum3, E_vdw_dum, E_qq_dum, 2, &
            (/lm_i, lm_j/), (/is, js/), (/box_i, box_j/), box_nrg_vdw_temp, box_nrg_qq_temp)
         E_vdw_i = box_nrg_vdw_temp(1)
         E_vdw_j = box_nrg_vdw_temp(2)
         E_qq_i = box_nrg_qq_temp(1)
         E_qq_j = box_nrg_qq_temp(2)
      
         DEALLOCATE(box_nrg_vdw_temp, box_nrg_qq_temp)
      ELSE
         CALL Compute_Molecule_Nonbond_Inter_Energy(lm_i,is,E_vdw_i,E_qq_i,inter_overlap)
         IF (inter_overlap) THEN
         END IF
         IF (.NOT. inter_overlap) THEN
            CALL Compute_Molecule_Nonbond_Inter_Energy(lm_j,js,E_vdw_j,E_qq_j,inter_overlap)
            IF (.NOT. inter_overlap) THEN
               CALL Compute_Molecule_Nonbond_Intra_Energy(lm_i,is, E_intra_vdw_i,E_intra_qq_i,E_periodic_qq_i,inter_overlap)
               CALL Compute_Molecule_Nonbond_Intra_Energy(lm_j,js, E_intra_vdw_j,E_intra_qq_j,E_periodic_qq_j,inter_overlap)
            END IF
         END IF
      END IF
   END IF


   IF (inter_overlap)  THEN
      err_msg = ""
      err_msg(1) = "Attempted to switch molecule " // TRIM(Int_To_String(im_i)) // &
                  " of species " // TRIM(Int_To_String(is)) // &
                  " and molecule " // TRIM(Int_To_String(im_j)) // &
                  " of species " // TRIM(Int_To_String(js))
      IF (nbr_boxes > 1) THEN
         err_msg(1) = err_msg(1) // " in box "  // TRIM(Int_To_String(box_i)) // &
                                    " and box " // TRIM(Int_To_String(box_j))
      END IF
      err_msg(2) = "but the molecule energy is too high"
      IF (start_type(box_i) == "make_config" .OR. start_type(box_j) == "make_config") THEN
         err_msg(3) = "Try increasing Rcutoff_Low, increasing the box size, or "
         err_msg(4) = "decreasing the initial number of molecules"
     END IF
     CALL Clean_Abort(err_msg, "Identity Switch")
   END IF

   CALL Save_Old_Cartesian_Coordinates(lm_i,is)
   CALL Save_Old_Cartesian_Coordinates(lm_j,js)

   !*****************************************************************************
   ! Step 8) Switch the two molecules by COM translation
   !
   !*****************************************************************************
   xcom_i = molecule_list(lm_i,is)%xcom
   ycom_i = molecule_list(lm_i,is)%ycom
   zcom_i = molecule_list(lm_i,is)%zcom

   ALLOCATE(dx_xcom_i(natoms(is)), dy_ycom_i(natoms(is)), dz_zcom_i(natoms(is)))
   ALLOCATE(dx_xcom_j(natoms(js)), dy_ycom_j(natoms(js)), dz_zcom_j(natoms(js)))
   dx_xcom_i = xcom_i - atom_list(:,lm_i,is)%rxp
   dy_ycom_i = ycom_i - atom_list(:,lm_i,is)%ryp
   dz_zcom_i = zcom_i - atom_list(:,lm_i,is)%rzp

   dx_xcom_j =  molecule_list(lm_j,js)%xcom - atom_list(:,lm_j,js)%rxp
   dy_ycom_j =  molecule_list(lm_j,js)%ycom - atom_list(:,lm_j,js)%ryp
   dz_zcom_j =  molecule_list(lm_j,js)%zcom - atom_list(:,lm_j,js)%rzp

   !Switch first molecule
   molecule_list(lm_i,is)%xcom = molecule_list(lm_j,js)%xcom
   molecule_list(lm_i,is)%ycom = molecule_list(lm_j,js)%ycom
   molecule_list(lm_i,is)%zcom = molecule_list(lm_j,js)%zcom
   atom_list(:,lm_i,is)%rxp = molecule_list(lm_i,is)%xcom - dx_xcom_i
   atom_list(:,lm_i,is)%ryp = molecule_list(lm_i,is)%ycom - dy_ycom_i
   atom_list(:,lm_i,is)%rzp = molecule_list(lm_i,is)%zcom - dz_zcom_i

   !Switch second molecule
   molecule_list(lm_j,js)%xcom = xcom_i
   molecule_list(lm_j,js)%ycom = ycom_i
   molecule_list(lm_j,js)%zcom = zcom_i
   atom_list(:,lm_j,js)%rxp = xcom_i - dx_xcom_j
   atom_list(:,lm_j,js)%ryp = ycom_i - dy_ycom_j
   atom_list(:,lm_j,js)%rzp = zcom_i - dz_zcom_j
   DEALLOCATE(dx_xcom_i, dy_ycom_i, dz_zcom_i)
   DEALLOCATE(dx_xcom_j, dy_ycom_j, dz_zcom_j)

   !change box identity
   molecule_list(lm_i,is)%which_box = box_j
   molecule_list(lm_j,js)%which_box = box_i

   !Delete i's locate in box_i
   IF (im_i < nmols(is,box_i)) THEN
      DO k = im_i + 1, nmols(is,box_i)
         locate(k-1,is,box_i) = locate(k,is,box_i)
      END DO
   END IF
   locate(nmols(is,box_i),is,box_i) = 0

   !Delete j's locate in box_j
   IF (im_j < nmols(js,box_j)) THEN
      DO k = im_j + 1, nmols(js,box_j)
         locate(k-1,js,box_j) = locate(k,js,box_j)
      END DO
   END IF
   locate(nmols(js,box_j),js,box_j) = 0

   !Update number of molecules
   nmols(is,box_i) = nmols(is,box_i) - 1
   nmols(js,box_j) = nmols(js,box_j) - 1
   nmols(is,box_j) = nmols(is,box_j) + 1
   nmols(js,box_i) = nmols(js,box_i) + 1

   !Add j's locate to box_i
   im_i = nmols(is, box_j)
   locate(im_i,is,box_j) = lm_i

   !Add i's locate to box_j
   im_j = nmols(js, box_i)
   locate(im_j,js,box_i) = lm_j

   CALL Fold_Molecule(lm_i,is,box_j)
   CALL Fold_Molecule(lm_j,js,box_i)

  !*****************************************************************************
  ! Step 9) Rotate the molecules if desired
  !
  !*****************************************************************************
  IF (rotations .GT. 0) THEN
      CALL Bias_Rotate(lm_i,is,box_j,rotations,P_bias_rot_i,rot_overlap_i)
       P_bias = P_bias * P_bias_rot_i

        CALL Bias_Rotate(lm_j,js,box_i,rotations,P_bias_rot_i,rot_overlap_j)
        P_bias = P_bias * P_bias_rot_j
  END IF



  !*****************************************************************************
  ! Step 10) Calculate the energy change of each box after rotation
  !
  !*****************************************************************************
   IF ( .NOT. rot_overlap_i .AND. .NOT. rot_overlap_j) THEN
      IF (box_i .EQ. box_j) THEN
         CALL Compute_MoleculeCollection_Nonbond_Inter_Energy(2, (/lm_i, lm_j/), (/is, js/), &
            E_vdw_move, E_qq_move, inter_overlap)
      ELSE
         CALL Compute_Molecule_Nonbond_Inter_Energy(lm_i,is,E_vdw_move_i,E_qq_move_i,inter_overlap)
         IF (.NOT. inter_overlap) THEN
            CALL Compute_Molecule_Nonbond_Inter_Energy(lm_j,js,E_vdw_move_j,E_qq_move_j,inter_overlap)
         END IF
      END IF
   END IF

   IF (inter_overlap .OR. rot_overlap_i .OR. rot_overlap_j) THEN
      CALL Revert_Switch
     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A16,X,X,I5,X,I5,X,I3,X,I3,X,L8,X,L8,X,9X,X,F9.3)') &
             i_mcstep, 'identity switch' , lm_i, lm_j,is, js, box_i, box_j, accept, ln_pacc
     END IF

   ELSE !no overlap
      dE = 0.0_DP
      dE_i = 0.0_DP
      dE_j = 0.0_DP

      !Ewald charge code:
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
               CALL Get_Position_Alive(lm_i, is, position_i)

               !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
               cos_mol_old_i(:) = cos_mol(1:nvecs(box),position_i)
               sin_mol_old_i(:) = sin_mol(1:nvecs(box),position_i)
               !$OMP END PARALLEL WORKSHARE

               CALL Update_System_Ewald_Reciprocal_Energy(lm_i,is,box,int_translation,E_reciprocal_move)
               dE = E_reciprocal_move + dE
            END IF
            IF (has_charge(js)) THEN
               ALLOCATE(cos_mol_old_j(nvecs(box)), sin_mol_old_j(nvecs(box)))
               CALL Get_Position_Alive(lm_j, js, position_j)

               !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
               cos_mol_old_j(:) = cos_mol(1:nvecs(box),position_j)
               sin_mol_old_j(:) = sin_mol(1:nvecs(box),position_j)
               !$OMP END PARALLEL WORKSHARE

               CALL Update_System_Ewald_Reciprocal_Energy(lm_j,js,box,int_translation,E_reciprocal_move)
               dE = E_reciprocal_move + dE
            END IF
            dE = dE - energy(box)%reciprocal
         !Can just use Update_System_Ewald_Reciprocal_Energy
         ELSE
            CALL Update_System_Ewald_Reciprocal_Energy(lm_i,is,box_j, &
               int_insertion,E_reciprocal_move_j)
            dE_j = dE_j + E_reciprocal_move_j - energy(box_j)%reciprocal
            CALL Update_System_Ewald_Reciprocal_Energy(lm_j,js,box_i, &
               int_insertion,E_reciprocal_move_i)
            dE_i = dE_i + E_reciprocal_move_i - energy(box_i)%reciprocal
         END IF

      END IF

      !Intra energies are not needed for the acceptance rule aside form E_periodic,
      !but they must be attributed to the correct box later

      CALL Compute_Molecule_Bond_Energy(lm_i,is,E_bond_i)
      CALL Compute_Molecule_Angle_Energy(lm_i,is,E_angle_i)
      CALL Compute_Molecule_Dihedral_Energy(lm_i,is,E_dihed_i)
      CALL Compute_Molecule_Improper_Energy(lm_i,is,E_improper_i)
      CALL Compute_Molecule_Nonbond_Intra_Energy(lm_i,is, E_intra_vdw_i,E_intra_qq_i,E_periodic_qq_move_i,inter_overlap)

      CALL Compute_Molecule_Bond_Energy(lm_j,js,E_bond_j)
      CALL Compute_Molecule_Angle_Energy(lm_j,js,E_angle_j)
      CALL Compute_Molecule_Dihedral_Energy(lm_j,js,E_dihed_j)
      CALL Compute_Molecule_Improper_Energy(lm_j,js,E_improper_j)
      CALL Compute_Molecule_Nonbond_Intra_Energy(lm_j,js, E_intra_vdw_j,E_intra_qq_j,E_periodic_qq_move_j,inter_overlap)

      !These dE variables will only be used in the two box case
      dE_bond = E_bond_j - E_bond_i
      dE_angle = E_angle_j - E_angle_i
      dE_dihed = E_dihed_j - E_dihed_i
      dE_improper = E_improper_j - E_improper_i
      dE_intra_vdw = E_intra_vdw_j - E_intra_vdw_i
      dE_intra_qq = E_intra_qq_j - E_intra_qq_i
      !End of Ewald Charge code


      !long range corrections
      IF (box_i .EQ. box_j) THEN
         IF (int_vdw_sum_style(box_i) == vdw_cut_tail) THEN
            CALL Compute_LR_correction(box_i,E_lrc)
            dE_lrc = E_lrc - energy(box_i)%lrc
         END IF
      ELSE
         IF (int_vdw_sum_style(box_i) == vdw_cut_tail) THEN
            CALL Compute_LR_correction(box_i, E_lrc_i)
            dE_lrc_i = E_lrc_i - energy(box_i)%lrc
         END IF
         IF (int_vdw_sum_style(box_j) == vdw_cut_tail) THEN
            CALL Compute_LR_correction(box_j, E_lrc_j)
            de_lrc_j = E_lrc_j - energy(box_j)%lrc
         END IF
      END IF

      !Compute difference with nonbonded energies only
      IF (box_i .EQ. box_j) THEN
         dE = dE + (E_vdw_move - E_vdw) + (E_qq_move - E_qq)
         dE = dE + (E_periodic_qq_move_j - E_periodic_qq_j) + (E_periodic_qq_move_i - E_periodic_qq_i)
         dE = dE + dE_lrc
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
         dE_i = dE_i + (E_periodic_qq_move_j - E_periodic_qq_i)
         dE_i = dE_i + dE_lrc_i

         dE_j = dE_j + (E_vdw_move_j - E_vdw_j) + (E_qq_move_j - E_qq_j)
         dE_j = dE_j + (E_periodic_qq_move_i - E_periodic_qq_j)
         dE_j = dE_j + dE_lrc_j

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

      !WRITE (*,*) "ln_Pacc:"
      !WRITE(*,*) ln_pacc
      IF (accept) THEN
         nsuccess(is,box_i)%switch = nsuccess(is,box_i)%switch + 1
         nsuccess(js,box_j)%switch = nsuccess(js,box_j)%switch + 1
         !accept the move and update global energies
         IF (box_i .EQ. box_j) THEN
            energy(box_i)%total = energy(box_i)%total + dE
            energy(box_i)%inter = energy(box_i)%inter + dE
            energy(box_i)%inter_vdw = energy(box_i)%inter_vdw + E_vdw_move - E_vdw
            energy(box_i)%inter_q   = energy(box_i)%inter_q   + E_qq_move  - E_qq

            IF(int_charge_sum_style(box_i) == charge_ewald .AND. (has_charge(is) .OR. has_charge(js))) THEN
               energy(box_i)%reciprocal = E_reciprocal_move
            END IF

      ELSE
         energy(box_i)%total = energy(box_i)%total + dE_i
         energy(box_i)%inter = energy(box_i)%inter + dE_i
         energy(box_i)%inter_vdw = energy(box_i)%inter_vdw + E_vdw_move_i - E_vdw_i
         energy(box_i)%inter_q = energy(box_i)%inter_q + E_qq_move_i - E_qq_i
         !Intra energies are not needed for the acceptance rule, but they must be attributed to the correct box now
         energy(box_i)%intra = energy(box_i)%intra + dE_bond + dE_angle + dE_dihed + dE_improper
         energy(box_i)%bond = energy(box_i)%bond + dE_bond
         energy(box_i)%angle = energy(box_i)%angle + dE_angle
         energy(box_i)%dihedral = energy(box_i)%dihedral + dE_dihed
         energy(box_i)%improper = energy(box_i)%improper + dE_improper
         energy(box_i)%intra_vdw = energy(box_i)%intra_vdw + dE_intra_vdw
         energy(box_i)%intra_q = energy(box_i)%intra_q + dE_intra_qq

         energy(box_j)%total = energy(box_j)%total + dE_j
         energy(box_j)%inter = energy(box_j)%inter + dE_j
         energy(box_j)%inter_vdw = energy(box_j)%inter_vdw + E_vdw_move_j - E_vdw
         energy(box_j)%inter_q = energy(box_j)%inter_q + E_qq_move_j - E_qq_j
         !Intra energies are not needed for the acceptance rule, but they must be attributed to the correct box now
         energy(box_j)%intra = energy(box_j)%intra - dE_bond - dE_angle - dE_dihed - dE_improper
         energy(box_j)%bond = energy(box_j)%bond - dE_bond
         energy(box_j)%angle = energy(box_j)%angle - dE_angle
         energy(box_j)%dihedral = energy(box_j)%dihedral - dE_dihed
         energy(box_j)%improper = energy(box_j)%improper - dE_improper
         energy(box_j)%intra_vdw = energy(box_j)%intra_vdw - dE_intra_vdw
         energy(box_j)%intra_q = energy(box_j)%intra_q - dE_intra_qq

         IF(int_charge_sum_style(box_i) == charge_ewald .AND. (has_charge(is) .OR. has_charge(js))) THEN
            energy(box_i)%reciprocal = E_reciprocal_move_i
         END IF
         IF(int_charge_sum_style(box_j) == charge_ewald .AND. (has_charge(is) .OR. has_charge(js))) THEN
            energy(box_j)%reciprocal = E_reciprocal_move_j
         END IF

        END IF

        IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
        IF (ALLOCATED(cos_mol_old_i)) DEALLOCATE(cos_mol_old_i)
        IF (ALLOCATED(sin_mol_old_i)) DEALLOCATE(sin_mol_old_i)
        IF (ALLOCATED(sin_mol_old_j)) DEALLOCATE(sin_mol_old_j)
        IF (ALLOCATED(sin_mol_old_j)) DEALLOCATE(sin_mol_old_j)

     ELSE
        CALL Revert_Switch
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
        IF (l_pair_nrg) THEN
            CALL Reset_Molecule_Pair_Interaction_Arrays(lm_i,is,box_i)
            CALL Reset_Molecule_Pair_Interaction_Arrays(lm_j,js,box_j)
        END IF
     END IF

     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A16,X,X,I5,X,I5,X,I3,X,I3,X,L8,X,L8,X,9X,X,F9.3)') &
             i_mcstep, 'identity switch' , lm_i, lm_j,is, js, box_i, box_j, accept, ln_pacc
     END IF

   END IF

   CONTAINS
   !***************************************************************************
   SUBROUTINE Revert_Switch

      CALL Revert_Old_Cartesian_Coordinates(lm_i, is)
      CALL Revert_Old_Cartesian_Coordinates(lm_j, js)

      !change box identity
      molecule_list(lm_i,is)%which_box = box_i
      molecule_list(lm_j,js)%which_box = box_j

      !Delete j's locate in box_i
      IF (im_i < nmols(is,box_j)) THEN
         DO k = im_i + 1, nmols(is,box_j)
            locate(k-1,is,box_j) = locate(k,is,box_j)
         END DO
      END IF
      locate(nmols(is,box_j),is,box_j) = 0

      !Delete i's locate in box_j
      IF (im_j < nmols(js,box_i)) THEN
         DO k = im_j + 1, nmols(js,box_i)
            locate(k-1,js,box_i) = locate(k,js,box_i)
         END DO
      END IF
      locate(nmols(js,box_i),js,box_i) = 0

      !Update number of molecules
      nmols(is,box_i) = nmols(is,box_i) + 1
      nmols(js,box_j) = nmols(js,box_j) + 1
      nmols(is,box_j) = nmols(is,box_j) - 1
      nmols(js,box_i) = nmols(js,box_i) - 1

      !Add i's locate to box_i
      im_i = nmols(is, box_i)
      locate(im_i,is,box_i) = lm_i

      !Add j's locate to box_j
      im_j = nmols(js, box_j)
      locate(im_j,js,box_j) = lm_j

      CALL Fold_Molecule(lm_i,is,box_i)
      CALL Fold_Molecule(lm_j,js,box_j)

   END SUBROUTINE Revert_Switch

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

      REAL(DP) :: E_inter_qq_rot, E_inter_vdw_rot

      ! Initialize
      P_bias = 1.0_DP
      P_bias_rev = 1.0_DP
      rot_overlap = .FALSE.
      nrg = 0.0_DP
      weight = 0.0_DP

      !WRITE (*,*) "Called Bias_Rotate"

      !CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw_rot, &
      !    E_inter_qq_rot,overlap)

      !nrg_original = E_inter_vdw_rot + E_inter_qq_rot
      !wt_original = DEXP(-1.0 * beta(ibox) * nrg_original)


        DO i_rot = 1, rotations
          CALL Rotate_Molecule_Eulerian(lm,is)

          !save rotation coordinates
          DO ia = 1, natoms(is)
            rtrial(ia,i_rot)%rxp = atom_list(ia,lm,is)%rxp
            rtrial(ia,i_rot)%ryp = atom_list(ia,lm,is)%ryp
            rtrial(ia,i_rot)%rzp = atom_list(ia,lm,is)%rzp
          END DO

          overlap = .FALSE.
          CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw_rot, &
                  E_inter_qq_rot,overlap)
          nrg(i_rot) = E_inter_vdw_rot + E_inter_qq_rot

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

          IF (i_rot > 1) weight(i_rot) = weight(i_rot-1) + weight(i_rot)

        END DO

      ! Reject the move if the total weight is still zero
      IF (weight(rotations) == 0.0_DP) THEN
        rot_overlap = .TRUE.
        RETURN
      END IF

      ! Choose one from Rosenbluth sampling
      rand_no = rranf() * weight(rotations)

      DO i_rot = 1, rotations
        IF (rand_no < weight(i_rot)) EXIT
      END DO

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

      !COM was not stored, so recalculate
      CALL Get_COM(lm,is)


      !Now we repeat all that in reverse to get the P_bias reverse
      !***************************************!

        DO i_rot = 2, rotations
          CALL Rotate_Molecule_Eulerian(lm,is)
          !don't need to save rotation coordinates
          overlap = .FALSE.
          CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw_rot, &
                  E_inter_qq_rot,overlap)
          nrg_rev(i_rot) = E_inter_vdw_rot + E_inter_qq_rot

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
          weight_rev(i_rot) = weight_rev(i_rot-1) + weight_rev(i_rot)
        END DO

      ! For the reverse move, we've already selected the first rotation
      P_bias_rev = P_bias_rev * weight_rev(1) / weight_rev(rotations)

      P_bias = P_bias_rev / P_bias
    END SUBROUTINE Bias_Rotate

END SUBROUTINE Identity_Switch

