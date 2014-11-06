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

SUBROUTINE Atom_Displacement(this_box)

  !********************************************************************************
  !
  ! This subroutine performs atom translation of terminal atoms (atoms with exactly
  ! one bond). This is useful to sample intramolecule degrees of freedom of molecules
  ! containing branch point.
  !
  !
  !
  !***********************************************************************************
  !
  ! Algorithm
  !
  ! A species that contains at least one terminal atom is selected on the basis
  ! of its mole fraction.
  !
  ! A molecule of species is chosen randomly.
  !
  ! A terminal atom is selected randomly.
  !
  ! Perturbation in coordinates of this atom is carried out in spherical coordinates
  !
  ! The move is accepted or rejected based on exp(-\beta \Delta E), where \Delta E is
  ! change in energy due to the atom displacment.
  !
  !*************************************************************************************
  !
  ! Called By
  !
  ! gcmc_driver
  ! gemc_driver
  ! nvtmc_driver
  ! nptmc_driver
  ! nvt_mc_ring_fragment
  !
  ! Revision history
  !
  ! 12/10/13 : Beta Release
  !
  !*************************************************************************************


  USE Run_Variables
  USE Random_Generators
  USE Energy_Routines
  USE Simulation_Properties, ONLY : Get_Index_Molecule

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: this_box

  INTEGER :: ibox, total_mols_ibox, nmolecule_species, ndisp_species, is, im
  INTEGER :: alive, this_atom, iatom, ref_atom, nmolecules_species

  INTEGER, DIMENSION(:), ALLOCATABLE :: total_mols, species_id

  REAL(DP) :: rand_no, E_intra_vdw_n, E_inter_vdw_n, E_intra_qq_n, E_inter_qq_n
  REAL(DP) :: E_inter_vdw_o, E_intra_vdw_o, E_inter_qq_o, E_intra_qq_o
  REAL(DP) :: e_recip_move,  delta_e, factor

  REAL(DP) :: e_ang_n, e_ang_o, e_dihed_n, e_dihed_o
  REAL(DP) :: e_improper_n, e_improper_o

  REAL(DP), DIMENSION(:), ALLOCATABLE :: x_box, x_species

  TYPE(Atom_Class), DIMENSION(:), ALLOCATABLE :: new_atom_list
  TYPE(Molecule_Class) :: new_molecule_list

  LOGICAL :: accept, accept_or_reject, overlap, theta_bound

  
  ! Figure out total number of molecules in each boxes that
  ! contain molecules with atoms that can be perturbed.

  ALLOCATE(total_mols(nbr_boxes), x_box(nbr_boxes))

  total_mols(:) = 0
  
  DO ibox = 1, nbr_boxes
     total_mols_ibox = 0
     DO is = 1, nspecies
        IF (species_list(is)%f_atom_disp) THEN
           total_mols_ibox = total_mols_ibox + nmols(is,ibox)
        END IF
     END DO
     total_mols(ibox) = total_mols_ibox
     IF (ibox > 1 ) THEN
        total_mols(ibox) = total_mols(ibox-1) + total_mols(ibox)
     END IF
  END DO
  
  IF (total_mols(nbr_boxes) == 0) RETURN

  ! obtain cumulative mole fractions for each box
  x_box(:) = REAL(total_mols(:),DP) / REAL(total_mols(nbr_boxes),DP)
  
  ! choose a box

  rand_no = rranf()
  
  DO ibox = 1, nbr_boxes
     IF ( rand_no <= x_box(ibox)) EXIT
  END DO

  this_box = ibox
 
  ! choose a species

 ! Choose a species, Make sure that the sum is carried out
  ! only over the species that possess fragments

  ALLOCATE(species_id(nspecies),x_species(nspecies))

  nmolecules_species = 0
  ndisp_species = 0
  x_species(:) = 0.0_DP

  DO is = 1, nspecies
     IF (species_list(is)%f_atom_disp ) THEN
        ndisp_species = ndisp_species + 1
        species_id(ndisp_species) = is
        nmolecules_species = nmolecules_species + nmols(is,this_box)
        x_species(ndisp_species) = REAL(nmolecules_species,DP)
     END IF
  END DO

  x_species(:) = x_species(:)  / REAL(nmolecules_species,DP)

  DO is = 2, ndisp_species
     x_species(is) = x_species(is) + x_species(is-1)
  END DO
  
  ! select a species

  rand_no = rranf()

  DO is = 1, ndisp_species
     IF (rand_no <= x_species(is)) EXIT
  END DO

  is = species_id(is)

  
  DEALLOCATE(x_box,x_species)
  DEALLOCATE(total_mols,species_id)

  ! pick a molecule of this species and get its index

  im = INT( rranf() * REAL(nmols(is,this_box),DP) ) + 1
  ntrials(is,this_box)%disp_atom = ntrials(is,this_box)%disp_atom + 1
  CALL Get_Index_Molecule(this_box,is,im,alive)

  ! Save coordinates

  CALL Save_Old_Cartesian_Coordinates(alive,is)

  ! choose an atom of this molecule at random

  this_atom = INT ( rranf() * REAL(species_list(is)%ndisp_atoms) ) + 1 

  iatom = species_list(is)%disp_atom_id(this_atom)
  ref_atom = species_list(is)%disp_atom_ref(this_atom)
!  write(*,*) iatom, ref_atom
!  write(*,*) atom_list(iatom,alive,is)
  atom_list(iatom,alive,is)%rxp = atom_list(iatom,alive,is)%rxp - atom_list(ref_atom,alive,is)%rxp
  atom_list(iatom,alive,is)%ryp = atom_list(iatom,alive,is)%ryp - atom_list(ref_atom,alive,is)%ryp
  atom_list(iatom,alive,is)%rzp = atom_list(iatom,alive,is)%rzp - atom_list(ref_atom,alive,is)%rzp

  ! Displace iatom

 

  CALL Change_Phi_Theta(iatom,alive,is,theta_bound)

  IF (theta_bound) THEN
     CALL Revert_Old_Cartesian_Coordinates(alive,is)
     RETURN
  END IF

  ! change coordinates of iatom to lab frame of reference

  atom_list(iatom,alive,is)%rxp = atom_list(iatom,alive,is)%rxp + atom_list(ref_atom,alive,is)%rxp
  atom_list(iatom,alive,is)%ryp = atom_list(iatom,alive,is)%ryp + atom_list(ref_atom,alive,is)%ryp
  atom_list(iatom,alive,is)%rzp = atom_list(iatom,alive,is)%rzp + atom_list(ref_atom,alive,is)%rzp
!   write(*,*) atom_list(iatom,alive,is)
  CALL Get_COM(alive,is)
  CALL Compute_Max_COM_Distance(alive,is)
  ! Compute the energy LJ interaction energy due to the move

  CALL Compute_Atom_Nonbond_Energy(iatom, alive, is, E_intra_vdw_n, &
       E_inter_vdw_n, E_intra_qq_n, E_inter_qq_n, overlap)

  IF (overlap) THEN
     CALL Revert_Old_Cartesian_Coordinates(alive,is)
     RETURN
  END IF
  

  ! If here then compute change in energy
!  write(*,*) 'new'
  CALL Compute_Molecule_Angle_Energy(alive,is,e_ang_n)
  CALL Compute_Molecule_Dihedral_Energy(alive,is,e_dihed_n)
  CALL Compute_Molecule_Improper_Energy(alive,is,e_improper_n)

  e_recip_move = 0.0_DP
  IF (int_charge_sum_style(this_box) == charge_ewald) THEN
     ! compute change in ewald reciprocal energy difference

     CALL Compute_Ewald_Reciprocal_Energy_Difference(alive,alive,is,this_box,&
          int_intra, e_recip_move)

  END IF
  
  ALLOCATE(new_atom_list(natoms(is)))
  new_atom_list(:) = atom_list(:,alive,is)
  new_molecule_list = molecule_list(alive,is)

  CALL Revert_Old_Cartesian_Coordinates(alive,is)

  ! Now compute the energy in old configuration

  CALL Compute_Molecule_Angle_Energy(alive,is,e_ang_o)
  CALL Compute_Molecule_Dihedral_Energy(alive,is,e_dihed_o)
  CALL Compute_Molecule_Improper_Energy(alive,is,e_improper_o)

  CALL Compute_Atom_Nonbond_Energy(iatom, alive, is,E_intra_vdw_o, &
       E_inter_vdw_o, E_intra_qq_o, E_inter_qq_o, overlap)
  
  
  delta_e = E_intra_vdw_n + E_inter_vdw_n + E_intra_qq_n + E_inter_qq_n + E_ang_n + E_dihed_n - &
            E_intra_vdw_o - E_inter_vdw_o - E_intra_qq_o - E_inter_qq_o - E_ang_o - E_dihed_o
  delta_e = delta_e + e_improper_n - e_improper_o
  delta_e = delta_e + e_recip_move

  factor = beta(this_box) * delta_e



  accept = accept_or_reject(factor)
!  write(*,*) accept
!!$  write(*,*) 'delta_e' , delta_e, E_intra_vdw_n, E_intra_vdw_o
!!$  write(*,*) 'angle', E_ang_n, E_ang_o
!!$  write(*,*) 'dihedral', E_dihed_n, E_dihed_o
  IF (accept) THEN
     ! update energies and transfer coordinates
     atom_list(:,alive,is) = new_atom_list(:)
     molecule_list(alive,is) = molecule_list(alive,is)
     
     energy(this_box)%total = energy(this_box)%total + delta_e
     energy(this_box)%angle = energy(this_box)%angle + E_ang_n - E_ang_o
     energy(this_box)%dihedral = energy(this_box)%dihedral + E_dihed_n - E_dihed_o
     energy(this_box)%improper = energy(this_box)%improper + e_improper_n - e_improper_o
     energy(this_box)%intra = energy(this_box)%intra + E_dihed_n - E_dihed_o + &
                                                       E_ang_n - E_ang_o
     energy(this_box)%intra_vdw = energy(this_box)%intra_vdw + E_intra_vdw_n - E_intra_vdw_o
     energy(this_box)%intra_q = energy(this_box)%intra_q + E_intra_qq_n - E_intra_qq_o
     energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + E_inter_vdw_n - E_inter_vdw_o
     energy(this_box)%inter_q = energy(this_box)%inter_q + E_inter_qq_n - E_inter_qq_o
!     write(*,*) energy
     IF (int_charge_sum_style(this_box) == charge_ewald) THEN
        energy(this_box)%ewald_reciprocal = energy(this_box)%ewald_reciprocal + e_recip_move
     END IF
     
     nsuccess(is,this_box)%disp_atom = nsuccess(is,this_box)%disp_atom + 1
     
  ELSE
     ! Reject the move, cooridnates and COM are already shifted to the old so nothing
     ! to here except restore the cos_sum and sin_sum arrays

     IF (int_charge_sum_style(this_box) == charge_ewald) THEN
        cos_sum(:,this_box) = cos_sum_old(:,this_box)
        sin_sum(:,this_box) = sin_sum_old(:,this_box)
     END IF

  END IF

!  read(*,*)
 
  DEALLOCATE(new_atom_list)

END SUBROUTINE Atom_Displacement
  

  
