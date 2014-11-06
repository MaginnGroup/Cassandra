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

SUBROUTINE Cut_N_Grow(this_box,mcstep)

  !********************************************************************************
  !
  ! The subroutine is used to sample conformations of molecules with cut and grow
  ! method proposed in Macedonia and Maginn, Mol. Phys., 1999.
  !
  ! Basic algorithm is described below.
  !
  ! A box is chosen at random depending based on its overall mole fraction.
  !
  ! A species from the box is chosen at random based on its mole fraction in the
  ! box.
  !
  ! A molecular of the species is chosen at random
  !
  ! One of the connections in the molecule is severed and randomly selected part
  ! of the molecule is deleted.
  !
  ! Configurational biasing is then used to regrow the deleted portion of the molecule
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
  !*********************************************************************************



  USE Run_Variables
  USE Energy_Routines
  USE Fragment_Growth, ONLY : Cut_Regrow, Single_Fragment_Regrowth
  USE Random_Generators, ONLY : rranf
  USE Simulation_Properties
  USE IO_Utilities, ONLY: Int_To_String
  USE Pair_Nrg_Routines

  IMPLICIT NONE

  INTEGER :: ibox, is, nmolecules_species, nmols_box, total_mols, this_box
  INTEGER :: nfrag_species, im, alive, frag_start, frag_end, frag_total, i
  INTEGER :: dumcount, iatom, int_phi, mcstep
  INTEGER :: alive_1, alive_2, locate_1, locate_2

  INTEGER, ALLOCATABLE, DIMENSION(:) :: species_id, frag_order
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: x_box, x_species
  REAL(DP), ALLOCATABLE :: x_old(:), y_old(:), z_old(:)
  REAL(DP), ALLOCATABLE :: dx(:), dy(:), dz(:)

  REAL(DP) :: rand_no, P_forward, P_reverse, E_bond_n, E_angle_n, E_dihed_n
  REAL(DP) :: E_intra_vdw_n, E_intra_qq_n, E_inter_vdw_n, E_inter_qq_n
  REAL(DP) :: E_intra_vdw_o, E_intra_qq_o, E_inter_vdw_o, E_inter_qq_o
  REAL(DP) :: E_bond_o, E_angle_o, E_dihed_o, delta_e_n, delta_e_o
  REAL(DP) :: E_improper_n, E_improper_o
  REAL(DP) :: E_selferf_n, E_selferf_o
  REAL(DP) :: E_reciprocal_move, factor, e_prev, delta_intra, phi
  REAL(DP) :: energy_olde, check_e
  REAL(DP) :: nrg_ring_frag_forward, nrg_ring_frag_reverse
  REAL(DP) :: lambda_for_cut
  REAL(DP) :: attempt_p

  LOGICAL ::  cbmc_overlap, accept, accept_or_reject, update_flag, superbad, overlap
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


  inside_start = .FALSE.
!  imp_Flag = .FALSE.
  attempt_p = 1.0_DP
  nrg_ring_frag_forward = 0.0
  nrg_ring_frag_reverse = 0.0

  ! Let us choose a box to pick the species and molecule to pick from.

  ALLOCATE(x_box(nbr_boxes), x_species(nspecies), species_id(nspecies))

  total_mols = SUM(nmols(:,:))

  IF(total_mols == 0) RETURN

  DO ibox = 1, nbr_boxes
     nmols_box = SUM(nmols(:,ibox))
     x_box(ibox) = REAL(nmols_box,DP)/REAL(total_mols,DP)
     IF (ibox > 1 ) THEN
        x_box(ibox) = x_box(ibox) + x_box(ibox-1)
     END IF
  END DO

  ! Choose a box

  rand_no = rranf()

  DO ibox = 1,nbr_boxes
     IF ( rand_no <= x_box(ibox)) EXIT
  END DO

  this_box = ibox
  energy_old = energy(this_box)
  tot_trials(this_box) = tot_trials(this_box) + 1
  DEALLOCATE(x_box)

  !----------------------------------

  ! Choose a species, Make sure that the sum is carried out
  ! only over the species that possess fragments
  
  ! NR: Changing based on probability based on species
  ! We might have situations where we don't want to do
  ! regrowth for some molecules types

  nmolecules_species = 0
  nfrag_species = 0
  x_species(:) = 0.0_DP

  DO is = 1, nspecies
     IF (nfragments(is) >= 1 ) THEN
        nfrag_species = nfrag_species + 1
        species_id(nfrag_species) = is
        nmolecules_species = nmolecules_species + nmols(is,this_box)
        x_species(nfrag_species) = REAL(nmolecules_species,DP)
     END IF
  END DO

  IF(nmolecules_species == 0) RETURN

!  x_species(:) = x_species(:)  / REAL(nmolecules_species,DP)

!  DO is = 2, nfrag_species
!     x_species(is) = x_species(is) + x_species(is-1)
!  END DO
  
  ! select a species

  rand_no = rranf()

  DO is = 1, nspecies
     IF (rand_no <= prob_growth_species(is)) EXIT
  END DO

  is = species_id(is)
  nmolecules_species = nmols(is,this_box)

  IF ( nmols(is,this_box) == 0 ) RETURN

  DEALLOCATE(x_species,species_id)

  ! Select a molecule at random for cutting 


     im = INT ( rranf() * nmolecules_species ) + 1
     ! Get the index of imth molecule of species is in this_box
     CALL Get_Index_Molecule(this_box,is,im,alive)



  ! Save the old coordinates of the molecule
!     energy_olde = energy(this_box)%inter_vdw
!     CALL Compute_Total_System_Energy(this_box,.FALSE.,overlap)
!     check_e = abs(energy(this_box)%inter_vdw - energy_olde)
!     if(check_e .GT. 0.0005) THEN
!        superbad = .true.
!        write(*,*) 'fubar'
!     end if

  CALL Save_Old_Cartesian_Coordinates(alive,is)


  IF (l_pair_nrg) CALL Store_Molecule_Pair_Interaction_Arrays(alive,is,this_box,E_inter_vdw_o, &
       E_inter_qq_o)
  
  l_charge = .FALSE.
  IF ((int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is))) l_charge = .TRUE.


     
  ! We will first cut part of the molecule and then
  ! determine if a suitable position is found. If so, then
  ! the weight of the existing chain will be computed followed
  ! by the energy of the molecule

  ! set deletion flag to false and call the CBMC Cut_Regrow routine

  del_FLAG = .FALSE.
  cbmc_overlap = .FALSE.
  del_overlap = .FALSE.

  P_forward = 1.0_DP

  ALLOCATE(frag_order(nfragments(is)))

  IF ( nfragments(is) == 1 ) THEN

    lambda_for_cut = molecule_list(alive,is)%cfc_lambda
    CALL Single_Fragment_Regrowth(alive,is)
    frag_total = 1

  ELSE

     nrg_ring_frag_forward = 0.0_DP

     lambda_for_cut = molecule_list(alive,is)%cfc_lambda
     CALL Cut_Regrow(alive,is,frag_start,frag_end,frag_order,frag_total,lambda_for_cut, &
          e_prev,P_forward, nrg_ring_frag_forward, cbmc_overlap, del_overlap)

  END IF

  ! The last trial in the CBMC move may have resulted into an overlap so some of 
  ! the atoms of that fragment might have 'exist' attribute as false. Make them alive

  atom_list(:,alive,is)%exist = .TRUE.


 ! Increment the total number of trials for growing frag_total 
 
  regrowth_trials(frag_total,is) = regrowth_trials(frag_total,is) + 1

  CALL Get_COM(alive,is)
  CALL Compute_Max_COM_Distance(alive,is)

  CALL Fold_Molecule(alive,is,this_box)

 ! write(*,*) 'cbmc_overlap', cbmc_overlap, this_box

  IF (cbmc_overlap) THEN
     CALL Revert_Old_Cartesian_Coordinates(alive,is)
     ! It may so happen that only part of the molecule was grown, so make
     ! all the atoms alive and set the cfc_lambda values
 
     DO iatom = 1, natoms(is)    
        atom_list(iatom,alive,is)%exist = .TRUE.
     END DO
     molecule_list(alive,is)%cfc_lambda = lambda_for_cut

     DEALLOCATE(frag_order)

     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,this_box)
     
     IF (ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
     IF (ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)
     
     RETURN

  END IF

  ! If here then the molecule was succesfully grown.

  ! We will calculate the intra and intermolecular energy changes and
  ! also dihedral angle energy change.

  E_selferf_n = 0.0_DP

  CALL Compute_Molecule_Bond_Energy(alive,is,E_bond_n)
  CALL Compute_Molecule_Angle_Energy(alive,is,E_angle_n)
  CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihed_n)
  CALL Compute_Molecule_Improper_Energy(alive,is,E_improper_n)
  CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw_n, E_intra_qq_n,intra_overlap)

     CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_inter_vdw_n,E_inter_qq_n,cbmc_overlap)


  ! Note that cbmc_overlap could be true when a molecule with single fragment
  ! is grown

  IF (cbmc_overlap) THEN
    
     CALL Revert_Old_Cartesian_Coordinates(alive,is)
     
     molecule_list(alive,is)%cfc_lambda = lambda_for_cut
     
     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,this_box)

     RETURN

  END IF

  ! Note that the energy difference will not include angle energies as these
  ! terms cancel out in the acceptance rule. However, they are needed in updating
  ! energies if the move is accepted

  delta_e_n = E_intra_vdw_n + E_intra_qq_n + E_inter_qq_n + E_inter_vdw_n + E_dihed_n + E_improper_n - &
              E_selferf_n

  IF (l_charge) THEN
     ! store cos_mol and sin_mol arrays
     
     ALLOCATE(cos_mol_old(nvecs(this_box)),sin_mol_old(nvecs(this_box)))
     CALL Get_Position_Alive(alive,is,position)
     
     !$OMP PARALLEL WORKSHARE  DEFAULT(SHARED)
     cos_mol_old(:) = cos_mol(1:nvecs(this_box),position)
     sin_mol_old(:) = sin_mol(1:nvecs(this_box),position)
     !$OMP END PARALLEL WORKSHARE

     ! Compute the change in Ewald reciprocal energy due to the move
     CALL Compute_Ewald_Reciprocal_Energy_Difference(alive,alive,is,this_box, &
          int_intra, E_reciprocal_move)
     delta_e_n = delta_e_n + E_reciprocal_move

  END IF

  ! We need to compute the energy of the old positions. So need to preserve
  ! coordinates of the current molecule

  ALLOCATE(new_atom_list(natoms(is)))

  DO iatom = 1,natoms(is)
     new_atom_list(iatom) = atom_list(iatom,alive,is)
  END DO
  new_molecule_list = molecule_list(alive,is)

  CALL Revert_Old_Cartesian_Coordinates(alive,is)
  
  ! obtain weight of the old positions

  del_FLAG = .TRUE.
  del_overlap = .FALSE.
  cbmc_overlap = .FALSE.

  P_reverse = 1.0_DP

  IF (nfragments(is) /= 1) THEN

     lambda_for_cut = molecule_list(alive,is)%cfc_lambda
     nrg_ring_frag_reverse = 0.0_DP
     CALL Cut_Regrow(alive,is,frag_start,frag_end,frag_order,frag_total, lambda_for_cut, & 
          e_prev, P_reverse, nrg_ring_frag_reverse, cbmc_overlap, del_overlap)

     atom_list(1:natoms(is),alive,is)%exist = .true.

     IF (del_overlap .or. cbmc_overlap) THEN
        
        ! positions of all the atoms may not have been reset to the old state
        ! so do this.
        
        CALL Revert_Old_Cartesian_Coordinates(alive,is)
        
        molecule_list(alive,is)%cfc_lambda = lambda_for_cut

        accept = .FALSE.

     END IF

  END IF

  
  CALL Get_COM(alive,is)
  CALL Compute_Max_COM_Distance(alive,is)

  DEALLOCATE(frag_order)

  IF ( .not. (del_overlap .or. cbmc_overlap) ) THEN

  ! Bonded interactions

     E_selferf_o = 0.0_DP

     CALL Compute_Molecule_Bond_Energy(alive,is,E_bond_o)
     CALL Compute_Molecule_Angle_Energy(alive,is,E_angle_o)
     CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihed_o)
     CALL Compute_Molecule_Improper_Energy(alive,is,E_improper_o)
     CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw_o, &
             E_intra_qq_o,intra_overlap)


     IF (.NOT. l_pair_nrg) CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_inter_vdw_o,E_inter_qq_o,cbmc_overlap)

     delta_e_o = E_intra_vdw_o + E_intra_qq_o + E_inter_vdw_o + E_inter_qq_o + E_dihed_o + E_improper_o - &
                 E_selferf_o

     factor = beta(this_box) * (delta_e_n - delta_e_o) + DLOG(P_forward) - DLOG(P_reverse)
     
     ! Modify factor to allow for ring biasing
     
     factor = factor + beta(this_box) * (nrg_ring_frag_reverse - nrg_ring_frag_forward)

     accept = accept_or_reject(factor)

  END IF



  IF (accept) THEN

     ! Change the coordinates of the molecule to new 
     DO iatom = 1,natoms(is)
        atom_list(iatom,alive,is) = new_atom_list(iatom)
     END DO
     molecule_list(alive,is) = new_molecule_list

     ! update energies for the boxes
     
! Change the coordinates of the molecule to new 

     CALL Get_Internal_Coordinates(alive,is)

     ! update energies for the boxes
     ! note that we include the change in angle energy for the system
     energy(this_box)%total = energy(this_box)%total + delta_e_n - delta_e_o &
          + e_angle_n - e_angle_o
     energy(this_box)%bond = energy(this_box)%bond + e_bond_n - e_bond_o
     energy(this_box)%angle = energy(this_box)%angle + e_angle_n - e_angle_o
     energy(this_box)%dihedral = energy(this_box)%dihedral + e_dihed_n - e_dihed_o
     energy(this_box)%improper = energy(this_box)%improper + E_improper_n - E_improper_o
     energy(this_box)%intra = energy(this_box)%intra + e_bond_n - e_bond_o + &
          e_angle_n - e_angle_o + e_dihed_n - e_dihed_o + E_improper_n - E_improper_o
     energy(this_box)%intra_vdw = energy(this_box)%intra_vdw + e_intra_vdw_n - e_intra_vdw_o
     energy(this_box)%intra_q = energy(this_box)%intra_q + e_intra_qq_n - e_intra_qq_o
     energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + e_inter_vdw_n - &
                                                               e_inter_vdw_o
     energy(this_box)%inter_q = energy(this_box)%inter_q + e_inter_qq_n - e_inter_qq_o

     IF (int_charge_sum_style(this_box) == charge_ewald) THEN
        energy(this_box)%ewald_reciprocal = energy(this_box)%ewald_reciprocal + &
             E_reciprocal_move
        energy(this_box)%ewald_self = energy(this_box)%ewald_self - (E_selferf_n - &
                                                                     E_selferf_o)
     END IF

     regrowth_success(frag_total,is) = regrowth_success(frag_total,is) + 1

     IF (l_charge) DEALLOCATE(cos_mol_old,sin_mol_old)
     IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)

      ! Fold the molecule in case the COM has moved out of cell boundary

  ELSE
     
     ! Positions of the chain is already set to the old position
     
     IF (l_charge) THEN
        
        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        cos_sum(1:nvecs(this_box),this_box) = cos_sum_old(1:nvecs(this_box),this_box)
        sin_sum(1:nvecs(this_box),this_box) = sin_sum_old(1:nvecs(this_box),this_box)
        cos_mol(1:nvecs(this_box),position) =cos_mol_old(:)
        sin_mol(1:nvecs(this_box),position) =sin_mol_old(:)
        !$OMP END PARALLEL WORKSHARE
        
     END IF

     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,this_box)
     
  END IF

!  CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihed_o)

  DEALLOCATE(new_atom_list)

END SUBROUTINE Cut_N_Grow


     
