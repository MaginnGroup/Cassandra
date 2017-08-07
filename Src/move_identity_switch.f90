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

  IMPLICIT NONE

  !variables for selecting the box, and individual molecules
  INTEGER  :: box
  INTEGER  :: is, js            ! species index
  INTEGER  :: im, jm            ! molecule index
  INTEGER  :: i_alive, j_alive  ! molecule locate
  REAL(DP) :: x_box(nbr_boxes)
  REAL(DP) :: randno

  !variables for calculating change in energy from move
  REAL(DP) :: dE, E_vdw, E_qq, E_vdw_move, E_qq_move, E_reciprocal_move
  REAL(DP) :: E_qq_dum, E_vdw_dum
  INTEGER :: dum1, dum2, dum3
  LOGICAL :: inter_overlap
  REAL(DP), DIMENSION(:), ALLOCATABLE :: box_nrg_vdw_temp, box_nrg_qq_temp

  !variables for peforming the switch
  REAL(DP) :: xcom_i, ycom_i, zcom_i
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dx_xcom_i, dy_ycom_i, dz_zcom_i
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dx_xcom_j, dy_ycom_j, dz_zcom_j

  ! Initialize variables from insert
  INTEGER  :: nmols_tot ! number of molecules in the system
  INTEGER  :: nmols_box(nbr_boxes)

  !acceptance variables
  REAL(DP) :: ln_pacc, success_ratio_i, success_ratio_j
  LOGICAL :: accept_or_reject

 ! Variables added for l_pair_nrg and reciprocal k vector storage
  REAL(DP), ALLOCATABLE :: cos_mol_old_i(:), sin_mol_old_i(:), cos_mol_old_j(:), sin_mol_old_j(:)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, TARGET :: cos_sum_old_idsw, sin_sum_old_idsw
  INTEGER :: position_i, position_j

  E_vdw_move = 0.0_DP
  E_qq_move = 0.0_DP
  E_vdw = 0.0_DP
  E_qq = 0.0_DP
  E_reciprocal_move = 0.0_DP
  inter_overlap = .FALSE.
  accept = .FALSE.
  box = 1


  !write (*, *) "Made it into identity switch! Exiting"

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
  END IF
  !box defaults to 1 otherwise

  !write (*, *) "Step 1:"
  !write (*, *) "Selected box"
  !write (*, *) box

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
  IF (nmols_tot < 2) THEN
     err_msg = ''
     err_msg(1) = 'No swappable molecules'
     CALL Clean_Abort(err_msg,'Identity_Switch_Move')
  END IF

  !check this
  is = INT(rranf() * nspecies) + 1

  !write (*, *) "Step 2:"
  !write (*, *) "Selected species:"
  !write (*, *) is


  !TODO: ONLY ERROR OUT IF SWAP ISN'T POSSIBLE
  ! check to make sure the selected species is swappable
  !IF (species_list(is)%int_insert == int_noinsert) THEN
  !  err_msg = ''
  !  err_msg(1) = 'Species ' // TRIM(Int_To_String(is)) // ' is not swappable'
  !  CALL Clean_Abort(err_msg,'Identity_Switch_Move')
  !END IF

  !*****************************************************************************
  ! Step 3) Select a molecule 'alive' from species 'is' with uniform probability
  !*****************************************************************************
  ! pick a molecule INDEX at random
  im = INT(rranf() * nmols(is,box)) + 1

  ! Obtain the LOCATE of this molecule
  i_alive = locate(im, is, box)

  !write (*, *) "Step 3:"
  !write (*, *) "Selected atom:"
  !write (*, *) im
  !write (*, *) i_alive

  !*****************************************************************************
  ! Step 4) Select a species 'js':
  !      -> Right now just done randomly!
  !*****************************************************************************

  js = INT(rranf() * (nspecies - 1)) + 1

  !avoid choosing same species as before
  IF (js >= is) THEN
    js = js + 1
  END IF

  !write (*, *) "Step 4:"
  !write (*, *) "Selected species:"
  !write (*, *) js

  !*****************************************************************************
  ! Step 5) Select a molecule 'alive' from species 'js' with uniform probability
  !*****************************************************************************
  jm = INT(rranf() * nmols(js,box)) + 1

  j_alive = locate(jm, js, box)

  !write (*, *) "Step 4:"
  !write (*, *) "Selected atom:"
  !write (*, *) jm
  !write (*, *) j_alive

  !*****************************************************************************
  ! Step 6) Calculate the change in box_in's potential energy from inserting
  !
  !*****************************************************************************
  !update trials
  tot_trials(box) = tot_trials(box) + 1
  ntrials(is,box)%switch = ntrials(is,box)%switch + 1
  ntrials(js,box)%switch = ntrials(js,box)%switch + 1

  !write (*, *) "Marker 1"

  ! obtain the energy of the molecule before the move.  Note that due to
  ! this move, the interatomic energies such as vdw and electrostatics will
  ! change. Also the ewald_reciprocal energy will change but there will
  ! be no change in intramolecular energies.

  IF (l_pair_nrg) THEN
    ALLOCATE(box_nrg_vdw_temp(2), box_nrg_qq_temp(2))
    CALL Store_Molecule_Pair_Interaction_Arrays(dum1, dum2, dum3, E_vdw_dum, E_qq_dum, 2, &
         (/i_alive, j_alive/), (/is, js/), (/box, box/), box_nrg_vdw_temp, box_nrg_qq_temp)
    E_vdw = box_nrg_vdw_temp(1)
    E_qq = box_nrg_vdw_temp(1)
    DEALLOCATE(box_nrg_vdw_temp, box_nrg_qq_temp)
  ELSE
     CALL Compute_MoleculeCollection_Nonbond_Inter_Energy(2, (/i_alive, j_alive/), (/is, js/), &
          E_vdw, E_qq, inter_overlap)
  END IF

  IF (inter_overlap)  THEN
     err_msg = ""

     err_msg(1) = "Attempted to move molecule " // TRIM(Int_To_String(im)) // &
                  " of species " // TRIM(Int_To_String(is)) // "NOT RIGHT ABOUT SPECIES RIGHT NOW"
     IF (nbr_boxes > 1) err_msg(1) = err_msg(1) // " in box " // TRIM(Int_To_String(box))
     err_msg(2) = "but the molecule energy is too high"
     IF (start_type(box) == "make_config" ) THEN
        err_msg(3) = "Try increasing Rcutoff_Low, increasing the box size, or "
        err_msg(4) = "decreasing the initial number of molecules"
     END IF
     CALL Clean_Abort(err_msg, "Identity Switch")
  END IF

  !write (*, *) "Marker 3"

  CALL Save_Old_Cartesian_Coordinates(i_alive,is)
  CALL Save_Old_Cartesian_Coordinates(j_alive,js)


  !******************************Actual*Switch*********************************!
  !Save COMs and COM to atom differences
  !write (*, *) "!************************!"
  !write (*, *) "Molecule i_alive COM (x, y, z):"
  !write (*, "(F30.25, F30.25, F30.25)") molecule_list(i_alive,is)%xcom, molecule_list(i_alive,is)%ycom, molecule_list(i_alive,is)%zcom
  !write (*, *) "Molecule i_alive Atom 1 coords (x,y,z):"
  !write (*, "(F30.25, F30.25, F30.25)") atom_list(1, i_alive,is)%rxp, atom_list(1, i_alive,is)%ryp, atom_list(1, i_alive,is)%rzp

  !write (*, *) "Molecule j_alive COM (x, y, z):"
  !write (*, "(F30.25, F30.25, F30.25)") molecule_list(j_alive,js)%xcom, molecule_list(j_alive,js)%ycom, molecule_list(j_alive,js)%zcom
  !write (*, *) "Molecule j_alive Atom 1 coords (x,y,z):"
  !write (*, "(F30.25, F30.25, F30.25)") atom_list(1, j_alive,js)%rxp, atom_list(1, j_alive,js)%ryp, atom_list(1, j_alive,js)%rzp

  !write (*, *) molecule_list(i_alive, is)%ycom
  !write (*, "(F30.25)") molecule_list(i_alive, is)%ycom

  xcom_i = molecule_list(i_alive,is)%xcom
  ycom_i = molecule_list(i_alive,is)%ycom
  zcom_i = molecule_list(i_alive,is)%zcom

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
  !write (*, *) "About to do atom_list for first molecule"
  !write (*, *) "First molecule_list xcoms"
  !write (*, "(F30.25, F30.25, F30.25)") molecule_list(i_alive,is)%xcom, molecule_list(i_alive,is)%ycom, molecule_list(i_alive,is)%zcom
  !write (*, *) "Now for dx_xcom_is"
  !write (*, "(F30.25, F30.25, F30.25)") dx_xcom_i(1), dy_ycom_i(1), dz_zcom_i(1)
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


  !write (*, *) "!************************!"
  !write (*, *) "Molecule i_alive COM (x, y, z):"
  !write (*, "(F30.25, F30.25, F30.25)") molecule_list(i_alive,is)%xcom, molecule_list(i_alive,is)%ycom, molecule_list(i_alive,is)%zcom
  !write (*, *) "Molecule i_alive Atom 1 coords (x,y,z):"
  !write (*, "(F30.25, F30.25, F30.25)") atom_list(1, i_alive,is)%rxp, atom_list(1, i_alive,is)%ryp, atom_list(1, i_alive,is)%rzp

  !write (*, *) "Molecule j_alive COM (x, y, z):"
  !write (*, "(F30.25, F30.25, F30.25)") molecule_list(j_alive,js)%xcom, molecule_list(j_alive,js)%ycom, molecule_list(j_alive,js)%zcom
  !write (*, *) "Molecule j_alive Atom 1 coords (x,y,z):"
  !write (*, "(F30.25, F30.25, F30.25)") atom_list(1, j_alive,js)%rxp, atom_list(1, j_alive,js)%ryp, atom_list(1, j_alive,js)%rzp

  !write (*, *) "!************************!"
  !write (*, *) "Finished step 6"
  !****************************************************************************!
  !*****************************************************************************
  ! Step 7) Calculate the change in box_in's potential energy from inserting
  !
  !*****************************************************************************

  !TODO: Not necessary until doing regrowth stuff
  !CALL Fold_Molecule(i_alive, is, box)
  !CALL Fold_Molecule(j_alive, js, box)

  write(*,*) "Computing Nonbond Inter Energy WITH COLLECTION"
  CALL Compute_MoleculeCollection_Nonbond_Inter_Energy(2, (/i_alive, j_alive/), (/is, js/), &
      E_vdw_move, E_qq_move, inter_overlap)
  write (*,*) "Done computing Nonbond Inter Energy"

  IF (inter_overlap) THEN
    !write (*, *) "Overlap found!"
    CALL Revert_Old_Cartesian_Coordinates(i_alive, is)
    CALL Revert_Old_Cartesian_Coordinates(j_alive, js)

    IF (l_pair_nrg) THEN
       CALL Reset_Molecule_Pair_Interaction_Arrays(i_alive,is,box)
       CALL Reset_Molecule_Pair_Interaction_Arrays(j_alive,js,box)
    END IF

    IF (verbose_log) THEN
      WRITE(logunit,'(X,I9,X,A16,X,I5,X,I3,X,I3,X,L8,X,9X,X,A9)') &
            i_mcstep, 'identity switch' , i_alive, is, box, accept, 'overlap'
    END IF
  ELSE !no overlap
     !write (*, *) "No overlap found!"
     dE = 0.0_DP

     IF ((int_charge_sum_style(box) == charge_ewald) .AND. (has_charge(is) .OR. has_charge(js))) THEN
         !write (*, *) "Doing ewald stuff!"
        !TODO:Eventually rewrite Update_System_Ewald_Reciprocal_Energy!

        !Update_System_Ewald_Reciprocal_Energy will do this, but it might override it!
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
     END IF
     !write (*, *) "Done with big if statement"


     !compute difference in old and new energy
     dE = dE + (E_vdw_move - E_vdw) + (E_qq_move - E_qq)

     write (*,*) "Moment of truth... dE is.... duh duh duh duh!"
     write (*,*) dE
     write (*,*) "How were the other stuff?"
     write (*,*) E_vdw_move
     write (*,*) E_vdw
     write (*,*) E_qq_move
     write (*,*) E_qq


     IF (int_sim_type == sim_nvt_min) THEN
        IF (dE  <= 0.0_DP) THEN
           !write (*, *) "Accepted!"
           accept = .TRUE.
        END IF
     ELSE
         ln_pacc = beta(box) * dE
         accept = accept_or_reject(ln_pacc)
         !write (*,*) accept
     END IF

     IF (accept) THEN
        !accept the move and update global energies
        write (*,*) energy(box)%total
        energy(box)%total = energy(box)%total + dE
        energy(box)%inter = energy(box)%inter + dE
        energy(box)%inter_vdw = energy(box)%inter_vdw + E_vdw_move - E_vdw
        energy(box)%inter_q   = energy(box)%inter_q   + E_qq_move  - E_qq

        IF(int_charge_sum_style(box) == charge_ewald .AND. (has_charge(is) .OR. has_charge(js))) THEN
           energy(box)%reciprocal = E_reciprocal_move
        END IF

        !TODO: update success counter
        nsuccess(is,box)%switch = nsuccess(is,box)%switch + 1
        nsuccess(js,box)%switch = nsuccess(js,box)%switch + 1
        nequil_success(is,box)%switch = nequil_success(is,box)%switch + 1
        nequil_success(js,box)%switch = nequil_success(js,box)%switch + 1

        IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
        IF (ALLOCATED(cos_mol_old_i)) DEALLOCATE(cos_mol_old_i)
        IF (ALLOCATED(sin_mol_old_i)) DEALLOCATE(sin_mol_old_i)
        IF (ALLOCATED(sin_mol_old_j)) DEALLOCATE(sin_mol_old_j)
        IF (ALLOCATED(sin_mol_old_j)) DEALLOCATE(sin_mol_old_j)
     ELSE
        ! Revert to the old coordinates of atoms and com of the molecule
        CALL Revert_Old_Cartesian_Coordinates(i_alive,is)
        CALL Revert_Old_Cartesian_Coordinates(j_alive,js)

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
        IF (l_pair_nrg) THEN
            CALL Reset_Molecule_Pair_Interaction_Arrays(i_alive,is,box)
            CALL Reset_Molecule_Pair_Interaction_Arrays(j_alive,js,box)
        END IF
     ENDIF

     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A16,X,X,I5,X,I3,X,I3,X,L8,X,9X,X,F9.3)') &
             i_mcstep, 'identity switch' , i_alive, is, js, box, accept, ln_pacc
     END IF
  END IF

  !removed logging stuff, take another look!
  IF ( MOD(ntrials(is,box)%switch,nupdate) == 0 ) THEN
     IF ( int_run_type == run_equil ) THEN
        success_ratio_i = REAL(nequil_success(is,box)%displacement,DP)/REAL(nupdate,DP)
        success_ratio_j = REAL(nequil_success(js,box)%displacement,DP)/REAL(nupdate,DP)
     ELSE
        success_ratio_i = REAL(nsuccess(js,box)%switch,DP)/REAL(ntrials(js,box)%switch,DP)
        success_ratio_j = REAL(nsuccess(js,box)%switch,DP)/REAL(ntrials(js,box)%switch,DP)
     END IF

     WRITE(logunit,'(X,I9,X,A16,X,5X,X,I3,X,I3,X,F8.5)',ADVANCE='NO') &
           i_mcstep, 'identity switch', is, box, success_ratio_i

     WRITE(logunit,'(X,I9,X,A16,X,5X,X,I3,X,I3,X,F8.5)',ADVANCE='NO') &
           i_mcstep, 'identity switch', js, box, success_ratio_j
     !TODO: No need to change displacement like translation, but do we need to change rcut?
     END IF
     WRITE(logunit,*)
  END IF
END SUBROUTINE Identity_Switch

