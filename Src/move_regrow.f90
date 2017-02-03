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

SUBROUTINE Cut_N_Grow

  !*****************************************************************************
  !
  ! The subroutine is used to sample conformations of molecules with cut and
  ! grow method proposed in Macedonia and Maginn, Mol. Phys., 1999.
  !
  ! Basic algorithm is described below.
  !
  ! A box is chosen at random depending based on its overall mole fraction.
  !
  ! A species from the box is chosen at random based on its mole fraction in the
  ! box.
  !
  ! A molecule of the species is chosen at random
  !
  ! One of the connections in the molecule is severed and randomly selected part
  ! of the molecule is deleted.
  !
  ! Configurational biasing is then used to regrow the deleted portion of the
  ! molecule
  !
  ! Called by
  !
  !   gcmc_driver
  !   gemc_driver
  !   nptmc_driver
  !   nvtmc_driver
  !
  ! Revision history
  !
  !   12/10/13 : Beta Release
  !`
  !*****************************************************************************



  USE Global_Variables
  USE Energy_Routines
  USE Fragment_Growth, ONLY : Cut_Regrow, Single_Fragment_Regrowth
  USE Random_Generators, ONLY : rranf
  USE Simulation_Properties
  USE IO_Utilities, ONLY: Int_To_String
  USE Pair_Nrg_Routines

  IMPLICIT NONE

  INTEGER :: ibox, is, nmols_tot
  INTEGER :: nfrag_species, im, lm, frag_start, frag_end, frag_total, i
  INTEGER :: dumcount, iatom, int_phi, mcstep
  INTEGER :: nmols_box(nbr_boxes)

  INTEGER, ALLOCATABLE, DIMENSION(:) :: frag_order
  REAL(DP) :: x_box(nbr_boxes), x_species(nspecies)
  REAL(DP), ALLOCATABLE :: x_old(:), y_old(:), z_old(:)
  REAL(DP), ALLOCATABLE :: dx(:), dy(:), dz(:)

  REAL(DP) :: rand_no, ln_pseq, ln_pfor, ln_prev, E_bond_n, E_angle_n, E_dihed_n
  REAL(DP) :: E_intra_vdw_n, E_intra_qq_n, E_inter_vdw_n, E_inter_qq_n
  REAL(DP) :: E_intra_vdw_o, E_intra_qq_o, E_inter_vdw_o, E_inter_qq_o
  REAL(DP) :: E_bond_o, E_angle_o, E_dihed_o, delta_e_n, delta_e_o
  REAL(DP) :: E_improper_n, E_improper_o, E_periodic_qq
  REAL(DP) :: E_reciprocal_move, ln_pacc, e_prev, delta_intra, phi
  REAL(DP) :: energy_olde, check_e
  REAL(DP) :: nrg_ring_frag_forward, nrg_ring_frag_reverse
  REAL(DP) :: lambda_for_cut

  LOGICAL ::  cbmc_overlap, accept_or_reject, update_flag, overlap
  LOGICAL :: del_overlap, cbmc_overlap_f, del_overlap_f, intra_overlap
  LOGICAL  :: inside_start, inside_finish

  TYPE(Atom_Class), ALLOCATABLE :: new_atom_list(:)
  TYPE(molecule_Class) :: new_molecule_list
  TYPE(Energy_Class) :: energy_old

  ! Pair_Energy arrays and Ewald implementation

  INTEGER :: start, locate_im, count, this_species, position, this_im
!  REAL(DP), ALLOCATABLE :: pair_vdw_temp(:), pair_qq_temp(:)
  REAL(DP), ALLOCATABLE :: cos_mol_old(:), sin_mol_old(:)

  LOGICAL :: l_charge

  E_reciprocal_move = 0.0_DP

  inside_start = .FALSE.
!  imp_Flag = .FALSE.
  nrg_ring_frag_forward = 0.0
  nrg_ring_frag_reverse = 0.0
  accept = .FALSE.

  ! Let us choose a box to pick the species and molecule to pick from.

  ! Sum the number of regrowable molecules in each box
  nmols_tot = 0
  DO ibox = 1, nbr_boxes
     nmols_box(ibox) = 0
     DO is = 1, nspecies
        IF (is == 1 .AND. prob_growth_species(is) > 0.0_DP) THEN
           nmols_tot = nmols_tot + nmols(is,ibox)
           nmols_box(ibox)  = nmols_box(ibox)  + nmols(is,ibox)
        ELSE IF (is > 1) THEN
           IF (prob_growth_species(is) > prob_growth_species(is-1)) THEN
              nmols_tot = nmols_tot + nmols(is,ibox)
              nmols_box(ibox)  = nmols_box(ibox)  + nmols(is,ibox)
           END IF
        END IF
     END DO
     x_box(ibox) = REAL(nmols_box(ibox),DP)
  END DO

  IF (nmols_tot == 0) THEN
     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,5X,X,3X,X,I3,X,L8,X,9X,X,A9)') &
             i_mcstep, 'regrow' , ibox, accept, 'no mols'
     END IF
     RETURN
  END IF

  DO ibox = 1, nbr_boxes
     x_box(ibox) = x_box(ibox)/REAL(nmols_tot,DP)
     IF (ibox > 1 ) THEN
        x_box(ibox) = x_box(ibox) + x_box(ibox-1)
     END IF
  END DO

  ! Choose a box
  rand_no = rranf()

  DO ibox = 1,nbr_boxes
     IF ( rand_no <= x_box(ibox)) EXIT
  END DO

  energy_old = energy(ibox)
  tot_trials(ibox) = tot_trials(ibox) + 1

  !----------------------------------

  ! Choose species based on the mol fraction, using Golden sampling
  DO is = 1, nspecies
     IF ( is == 1 ) THEN 
       IF ( prob_growth_species(is) > 0. ) THEN
         x_species(is) = REAL(nmols(is,ibox), DP)/REAL(nmols_box(ibox),DP)
       ELSE
         x_species(is) = 0.0_DP
       END IF
     ELSE IF ( is > 1 ) THEN
       IF ( prob_growth_species(is) > prob_growth_species(is-1) ) THEN
         x_species(is) = REAL(nmols(is,ibox), DP)/REAL(nmols_box(ibox),DP)
       ELSE
         x_species(is) = 0.0_DP
       END IF
       x_species(is) = x_species(is) + x_species(is-1)
     END IF
  END DO

  rand_no = rranf()
  DO is = 1, nspecies
     IF( rand_no <= x_species(is)) EXIT
  END DO

  ! Error check
  IF ( is == 1 ) THEN
    IF ( prob_growth_species(is) == 0. ) THEN
       err_msg = ''
       err_msg(1) = 'Probability of regrowing molecules of species ' // &
                    TRIM(Int_To_String(is)) // ' is zero'
       CALL Clean_Abort(err_msg, 'Cut_N_Grow')
    END IF
  ELSE IF ( is > 1 ) THEN
    IF ( prob_growth_species(is) == prob_growth_species(is-1) ) THEN
       err_msg = ''
       err_msg(1) = 'Probability of regrowing molecules of species ' // &
                    TRIM(Int_To_String(is)) // ' is zero'
       CALL Clean_Abort(err_msg, 'Cut_N_Grow')
    END IF
  END IF


  IF ( nmols(is,ibox) == 0 ) THEN
     err_msg = ''
     err_msg(1) = 'No regrowable molecules of species ' // TRIM(Int_To_String(is)) &
               // ' in box ' // TRIM(Int_To_String(ibox))
     CALL Clean_Abort(err_msg,'Cut_N_Grow')
  END IF

  ! Select a molecule at random for cutting 
  im = INT ( rranf() * nmols(is,ibox) ) + 1
  lm = locate(im,is,ibox)

  CALL Save_Old_Cartesian_Coordinates(lm,is)


  IF (l_pair_nrg) CALL Store_Molecule_Pair_Interaction_Arrays(lm,is,ibox,E_inter_vdw_o, &
                                                              E_inter_qq_o)
  
     
  ! We will first cut part of the molecule and then
  ! determine if a suitable position is found. If so, then
  ! the weight of the existing chain will be computed followed
  ! by the energy of the molecule

  ! set deletion flag to false and call the CBMC Cut_Regrow routine

  del_flag = .FALSE.
  cbmc_overlap = .FALSE.
  del_overlap = .FALSE.
  ln_pseq = 0.0_DP
  ln_pfor = 0.0_DP

  ALLOCATE(frag_order(nfragments(is)))

  IF ( nfragments(is) == 1 ) THEN

    lambda_for_cut = molecule_list(lm,is)%frac
    CALL Single_Fragment_Regrowth(lm,is)
    frag_total = 1

  ELSE

     nrg_ring_frag_forward = 0.0_DP

     lambda_for_cut = molecule_list(lm,is)%frac
     CALL Cut_Regrow(lm,is,frag_start,frag_end,frag_order,frag_total,lambda_for_cut, &
          e_prev,ln_pseq,ln_pfor, nrg_ring_frag_forward, cbmc_overlap, del_overlap)

  END IF

  ! The last trial in the CBMC move may have resulted in an overlap, so some
  ! of the atoms of that fragment might have 'exist' attribute as false. i
  ! Make them all exist

  atom_list(:,lm,is)%exist = .TRUE.


  ! Increment the total number of trials for growing frag_total 
 
  regrowth_trials(frag_total,is) = regrowth_trials(frag_total,is) + 1

  CALL Get_COM(lm,is)
  CALL Compute_Max_COM_Distance(lm,is)

  CALL Fold_Molecule(lm,is,ibox)

  IF (cbmc_overlap) THEN
     CALL Revert_Old_Cartesian_Coordinates(lm,is)
     ! It may so happen that only part of the molecule was grown, so make
     ! all the atoms lm and set the frac values
 
     DO iatom = 1, natoms(is)    
        atom_list(iatom,lm,is)%exist = .TRUE.
     END DO
     molecule_list(lm,is)%frac = lambda_for_cut

     DEALLOCATE(frag_order)

     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(lm,is,ibox)
     
     IF (ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
     IF (ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)
     
     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,A9)') &
             i_mcstep, 'regrow' , lm, is, ibox, .FALSE., 'overlap'
     END IF

     RETURN

  END IF

  ! If here then the molecule was succesfully grown.

  ! We will calculate the intra and intermolecular energy changes and
  ! also dihedral angle energy change.

  CALL Compute_Molecule_Bond_Energy(lm,is,E_bond_n)
  CALL Compute_Molecule_Angle_Energy(lm,is,E_angle_n)
  CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihed_n)
  CALL Compute_Molecule_Improper_Energy(lm,is,E_improper_n)
  CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is,E_intra_vdw_n, &
       E_intra_qq_n,E_periodic_qq,intra_overlap)

  CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw_n, &
       E_inter_qq_n,cbmc_overlap)
  E_inter_qq_n = E_inter_qq_n + E_periodic_qq


  ! Note that cbmc_overlap could be true when a molecule with single fragment
  ! is grown

  IF (cbmc_overlap) THEN
    
     CALL Revert_Old_Cartesian_Coordinates(lm,is)
     
     molecule_list(lm,is)%frac = lambda_for_cut
     
     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(lm,is,ibox)

     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,A9)') &
             i_mcstep, 'regrow' , lm, is, ibox, .FALSE., 'overlap'
     END IF

     RETURN

  END IF

  ! Note that the energy difference will not include angle energies as these
  ! terms cancel out in the acceptance rule. However, they are needed in 
  ! updating energies if the move is accepted

  delta_e_n = E_intra_vdw_n + E_intra_qq_n + E_inter_qq_n + E_inter_vdw_n &
            + E_dihed_n + E_improper_n

  IF (int_charge_sum_style(ibox)  == charge_ewald .AND.&
      has_charge(is)) THEN
     ! store cos_mol and sin_mol arrays
     
     ALLOCATE(cos_mol_old(nvecs(ibox)),sin_mol_old(nvecs(ibox)))
     CALL Get_Position_Alive(lm,is,position)
     
     !$OMP PARALLEL WORKSHARE  DEFAULT(SHARED)
     cos_mol_old(:) = cos_mol(1:nvecs(ibox),position)
     sin_mol_old(:) = sin_mol(1:nvecs(ibox),position)
     !$OMP END PARALLEL WORKSHARE

     ! Compute the change in Ewald reciprocal energy due to the move
     CALL Update_System_Ewald_Reciprocal_Energy(lm,is,ibox, &
          int_intra,E_reciprocal_move)
     delta_e_n = delta_e_n &
               + (E_reciprocal_move - energy(ibox)%ewald_reciprocal)

  END IF

  ! We need to compute the energy of the old positions. So need to preserve
  ! coordinates of the current molecule

  ALLOCATE(new_atom_list(natoms(is)))

  DO iatom = 1,natoms(is)
     new_atom_list(iatom) = atom_list(iatom,lm,is)
  END DO
  new_molecule_list = molecule_list(lm,is)

  CALL Revert_Old_Cartesian_Coordinates(lm,is)
  
  ! obtain weight of the old positions

  del_flag = .TRUE.
  del_overlap = .FALSE.
  cbmc_overlap = .FALSE.

  ln_prev = 0.0_DP

  IF (nfragments(is) /= 1) THEN

     lambda_for_cut = molecule_list(lm,is)%frac
     nrg_ring_frag_reverse = 0.0_DP
     CALL Cut_Regrow(lm,is,frag_start,frag_end,frag_order,frag_total, &
          lambda_for_cut, e_prev, ln_pseq, ln_prev, nrg_ring_frag_reverse, &
          cbmc_overlap, del_overlap)

     atom_list(1:natoms(is),lm,is)%exist = .true.

     IF (del_overlap .or. cbmc_overlap) THEN
        err_msg = ""
        err_msg(1) = "Attempted to delete part of molecule " // TRIM(Int_To_String(im)) // &
                     " of species " // TRIM(Int_To_String(is))
        IF (nbr_boxes > 1) THEN
           err_msg(1) = err_msg(1) // " from box " // TRIM(Int_To_String(ibox))
        END IF
        err_msg(2) = "but the molecule energy is too high"
        IF (start_type(ibox) == "make_config" ) THEN
           err_msg(3) = "Try increasing Rcutoff_Low, increasing the box size, or "
           err_msg(4) = "decreasing the initial number of molecules"
        END IF
        CALL Clean_Abort(err_msg, "Cut_N_Grow")
     END IF

  END IF

  
  CALL Get_COM(lm,is)
  CALL Compute_Max_COM_Distance(lm,is)

  DEALLOCATE(frag_order)

  IF ( .not. (del_overlap .or. cbmc_overlap) ) THEN

  ! Bonded interactions

     CALL Compute_Molecule_Bond_Energy(lm,is,E_bond_o)
     CALL Compute_Molecule_Angle_Energy(lm,is,E_angle_o)
     CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihed_o)
     CALL Compute_Molecule_Improper_Energy(lm,is,E_improper_o)
     CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is,E_intra_vdw_o, &
             E_intra_qq_o,E_periodic_qq,intra_overlap)

     IF (.NOT. l_pair_nrg) CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw_o,E_inter_qq_o,cbmc_overlap)
     E_inter_qq_o = E_inter_qq_o + E_periodic_qq

     delta_e_o = E_intra_vdw_o + E_intra_qq_o + E_inter_vdw_o + E_inter_qq_o &
               + E_dihed_o + E_improper_o

     ln_pacc = beta(ibox) * (delta_e_n - nrg_ring_frag_forward) &
             - beta(ibox) * (delta_e_o - nrg_ring_frag_reverse) &
             + ln_pfor - ln_prev
     
     accept = accept_or_reject(ln_pacc)

  END IF



  IF (accept) THEN

     ! Change the coordinates of the molecule to new 
     DO iatom = 1,natoms(is)
        atom_list(iatom,lm,is) = new_atom_list(iatom)
     END DO
     molecule_list(lm,is) = new_molecule_list

     ! update energies for the boxes
     
! Change the coordinates of the molecule to new 

     CALL Get_Internal_Coordinates(lm,is)

     ! update energies for the boxes
     ! note that we include the change in angle energy for the system
     energy(ibox)%total = energy(ibox)%total + delta_e_n - delta_e_o &
          + e_angle_n - e_angle_o
     energy(ibox)%bond = energy(ibox)%bond + e_bond_n - e_bond_o
     energy(ibox)%angle = energy(ibox)%angle + e_angle_n - e_angle_o
     energy(ibox)%dihedral = energy(ibox)%dihedral + e_dihed_n - e_dihed_o
     energy(ibox)%improper = energy(ibox)%improper + E_improper_n - E_improper_o
     energy(ibox)%intra = energy(ibox)%intra + e_bond_n - e_bond_o + &
          e_angle_n - e_angle_o + e_dihed_n - e_dihed_o + E_improper_n - E_improper_o
     energy(ibox)%intra_vdw = energy(ibox)%intra_vdw + e_intra_vdw_n - e_intra_vdw_o
     energy(ibox)%intra_q = energy(ibox)%intra_q + e_intra_qq_n - e_intra_qq_o
     energy(ibox)%inter_vdw = energy(ibox)%inter_vdw + e_inter_vdw_n - &
                                                               e_inter_vdw_o
     energy(ibox)%inter_q = energy(ibox)%inter_q + e_inter_qq_n - e_inter_qq_o

     IF (int_charge_sum_style(ibox)  == charge_ewald .AND.&
         has_charge(is)) THEN
        energy(ibox)%ewald_reciprocal = E_reciprocal_move
     END IF

     regrowth_success(frag_total,is) = regrowth_success(frag_total,is) + 1

     IF (int_charge_sum_style(ibox)  == charge_ewald .AND.&
      has_charge(is)) DEALLOCATE(cos_mol_old,sin_mol_old)
     IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)

      ! Fold the molecule in case the COM has moved out of cell boundary

  ELSE
     
     ! Positions of the chain is already set to the old position
     
     IF (int_charge_sum_style(ibox)  == charge_ewald .AND.&
         has_charge(is)) THEN
        
        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        cos_sum(1:nvecs(ibox),ibox) = cos_sum_old(1:nvecs(ibox),ibox)
        sin_sum(1:nvecs(ibox),ibox) = sin_sum_old(1:nvecs(ibox),ibox)
        cos_mol(1:nvecs(ibox),position) =cos_mol_old(:)
        sin_mol(1:nvecs(ibox),position) =sin_mol_old(:)
        !$OMP END PARALLEL WORKSHARE
        
     END IF

     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(lm,is,ibox)
     
  END IF

!  CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihed_o)

  DEALLOCATE(new_atom_list)

  IF (verbose_log) THEN
    WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,F9.3)') &
          i_mcstep, 'regrow' , lm, is, ibox, accept, ln_pacc
  END IF

END SUBROUTINE Cut_N_Grow
