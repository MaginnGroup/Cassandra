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

SUBROUTINE Atom_Displacement

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
  ! ring_fragment_driver
  !
  ! Revision history
  !
  ! 12/10/13 : Beta Release
  !
  !*************************************************************************************


  USE Global_Variables
  USE Random_Generators
  USE Energy_Routines

  IMPLICIT NONE


  INTEGER :: ibox, total_mols_ibox, nmolecule_species, ndisp_species, is, im
  INTEGER :: lm, this_atom, iatom, ref_atom, nmolecules_species, mcstep

  INTEGER, DIMENSION(:), ALLOCATABLE :: total_mols, species_id

  REAL(DP) :: rand_no, E_intra_vdw_n, E_inter_vdw_n, E_intra_qq_n, E_inter_qq_n
  REAL(DP) :: E_inter_vdw_o, E_intra_vdw_o, E_inter_qq_o, E_intra_qq_o
  REAL(DP) :: e_recip_move,  delta_e, factor

  REAL(DP) :: e_ang_n, e_ang_o, e_dihed_n, e_dihed_o
  REAL(DP) :: e_improper_n, e_improper_o

  REAL(DP), DIMENSION(:), ALLOCATABLE :: x_box, x_species

  TYPE(Atom_Class), DIMENSION(:), ALLOCATABLE :: new_atom_list
  TYPE(Molecule_Class) :: new_molecule_list

  LOGICAL :: accept_or_reject, overlap, theta_bound

  
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
        nmolecules_species = nmolecules_species + nmols(is,ibox)
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

  im = INT( rranf() * REAL(nmols(is,ibox),DP) ) + 1
  ntrials(is,ibox)%disp_atom = ntrials(is,ibox)%disp_atom + 1
  lm = locate(im,is,ibox)

  ! Save coordinates

  CALL Save_Old_Cartesian_Coordinates(lm,is)

  ! choose an atom of this molecule at random

  this_atom = INT ( rranf() * REAL(species_list(is)%ndisp_atoms) ) + 1 

  iatom = species_list(is)%disp_atom_id(this_atom)
  ref_atom = species_list(is)%disp_atom_ref(this_atom)
!  write(*,*) iatom, ref_atom
!  write(*,*) atom_list(iatom,lm,is)
  atom_list(iatom,lm,is)%rxp = atom_list(iatom,lm,is)%rxp - atom_list(ref_atom,lm,is)%rxp
  atom_list(iatom,lm,is)%ryp = atom_list(iatom,lm,is)%ryp - atom_list(ref_atom,lm,is)%ryp
  atom_list(iatom,lm,is)%rzp = atom_list(iatom,lm,is)%rzp - atom_list(ref_atom,lm,is)%rzp

  ! Displace iatom

 

  CALL Change_Phi_Theta(iatom,lm,is,theta_bound)

  IF (theta_bound) THEN
     CALL Revert_Old_Cartesian_Coordinates(lm,is)
     RETURN
  END IF

  ! change coordinates of iatom to lab frame of reference

  atom_list(iatom,lm,is)%rxp = atom_list(iatom,lm,is)%rxp + atom_list(ref_atom,lm,is)%rxp
  atom_list(iatom,lm,is)%ryp = atom_list(iatom,lm,is)%ryp + atom_list(ref_atom,lm,is)%ryp
  atom_list(iatom,lm,is)%rzp = atom_list(iatom,lm,is)%rzp + atom_list(ref_atom,lm,is)%rzp
!   write(*,*) atom_list(iatom,lm,is)
  CALL Get_COM(lm,is)
  CALL Compute_Max_COM_Distance(lm,is)
  ! Compute the energy LJ interaction energy due to the move

  CALL Compute_Atom_Nonbond_Energy(iatom, lm, is, E_intra_vdw_n, &
       E_inter_vdw_n, E_intra_qq_n, E_inter_qq_n, overlap)

  IF (overlap) THEN
     CALL Revert_Old_Cartesian_Coordinates(lm,is)
     RETURN
  END IF
  

  ! If here then compute change in energy
!  write(*,*) 'new'
  CALL Compute_Molecule_Angle_Energy(lm,is,e_ang_n)
  CALL Compute_Molecule_Dihedral_Energy(lm,is,e_dihed_n)
  CALL Compute_Molecule_Improper_Energy(lm,is,e_improper_n)

  e_recip_move = 0.0_DP
  IF (int_charge_sum_style(ibox) == charge_ewald) THEN
     ! compute change in ewald reciprocal energy difference

     CALL Update_System_Ewald_Reciprocal_Energy(lm,is,ibox,int_translation,e_recip_move)

  END IF
  
  ALLOCATE(new_atom_list(natoms(is)))
  new_atom_list(:) = atom_list(:,lm,is)
  new_molecule_list = molecule_list(lm,is)

  CALL Revert_Old_Cartesian_Coordinates(lm,is)

  ! Now compute the energy in old configuration

  CALL Compute_Molecule_Angle_Energy(lm,is,e_ang_o)
  CALL Compute_Molecule_Dihedral_Energy(lm,is,e_dihed_o)
  CALL Compute_Molecule_Improper_Energy(lm,is,e_improper_o)

  CALL Compute_Atom_Nonbond_Energy(iatom, lm, is,E_intra_vdw_o, &
       E_inter_vdw_o, E_intra_qq_o, E_inter_qq_o, overlap)
  
  
  delta_e = E_intra_vdw_n + E_inter_vdw_n + E_intra_qq_n + E_inter_qq_n + E_ang_n + E_dihed_n - &
            E_intra_vdw_o - E_inter_vdw_o - E_intra_qq_o - E_inter_qq_o - E_ang_o - E_dihed_o
  delta_e = delta_e + e_improper_n - e_improper_o
  delta_e = delta_e + (e_recip_move - energy(ibox)%ewald_reciprocal)

  factor = beta(ibox) * delta_e



  accept = accept_or_reject(factor)
!  write(*,*) accept
!!$  write(*,*) 'delta_e' , delta_e, E_intra_vdw_n, E_intra_vdw_o
!!$  write(*,*) 'angle', E_ang_n, E_ang_o
!!$  write(*,*) 'dihedral', E_dihed_n, E_dihed_o
  IF (accept) THEN
     ! update energies and transfer coordinates
     atom_list(:,lm,is) = new_atom_list(:)
     molecule_list(lm,is) = molecule_list(lm,is)
     
     energy(ibox)%total = energy(ibox)%total + delta_e
     energy(ibox)%angle = energy(ibox)%angle + E_ang_n - E_ang_o
     energy(ibox)%dihedral = energy(ibox)%dihedral + E_dihed_n - E_dihed_o
     energy(ibox)%improper = energy(ibox)%improper + e_improper_n - e_improper_o
     energy(ibox)%intra = energy(ibox)%intra + E_dihed_n - E_dihed_o + &
                                                       E_ang_n - E_ang_o
     energy(ibox)%intra_vdw = energy(ibox)%intra_vdw + E_intra_vdw_n - E_intra_vdw_o
     energy(ibox)%intra_q = energy(ibox)%intra_q + E_intra_qq_n - E_intra_qq_o
     energy(ibox)%inter_vdw = energy(ibox)%inter_vdw + E_inter_vdw_n - E_inter_vdw_o
     energy(ibox)%inter_q = energy(ibox)%inter_q + E_inter_qq_n - E_inter_qq_o
!     write(*,*) energy
     IF (int_charge_sum_style(ibox) == charge_ewald) THEN
        energy(ibox)%ewald_reciprocal = energy(ibox)%ewald_reciprocal + e_recip_move
     END IF
     
     nsuccess(is,ibox)%disp_atom = nsuccess(is,ibox)%disp_atom + 1
     
  ELSE
     ! Reject the move, cooridnates and COM are already shifted to the old so nothing
     ! to here except restore the cos_sum and sin_sum arrays

     IF (int_charge_sum_style(ibox) == charge_ewald) THEN
        cos_sum(:,ibox) = cos_sum_old(:,ibox)
        sin_sum(:,ibox) = sin_sum_old(:,ibox)
     END IF

  END IF

  IF (verbose_log) THEN
    WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8)') i_mcstep, 'atom_disp' , lm, is, ibox, accept
  END IF

 
  DEALLOCATE(new_atom_list)

END SUBROUTINE Atom_Displacement
  
!**********************************************************************************
SUBROUTINE Change_Phi_Theta(this_atom,im,is,theta_bound)
!**********************************************************************************
  USE Global_Variables
  USE Random_Generators, ONLY : rranf

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_atom,im,is

  LOGICAL, INTENT(INOUT) :: theta_bound

  REAL(DP) :: this_x, this_y, this_z, rho, bond_length, theta, phi, dcostheta
  REAL(DP) :: dphi

  theta_bound = .false.


  ! Get spherical coordinates
  this_x = atom_list(this_atom,im,is)%rxp
  this_y = atom_list(this_atom,im,is)%ryp
  this_z = atom_list(this_atom,im,is)%rzp

  rho = this_x * this_x + this_y * this_y 
  bond_length = this_z * this_z + rho

  rho = DSQRT(rho)
  bond_length = DSQRT(bond_length)

  ! azimuthal angle
  theta = DACOS(this_z/bond_length)
  
  ! polar angle
  phi = DASIN(this_y/rho)


  IF ( this_x < 0.0_DP ) THEN
     
     phi = PI - DASIN(this_y/rho)

  END IF

  ! Change theta and phi
  dcostheta = (2.0_DP * rranf() - 1.0_DP ) * delta_cos_max 
  dphi = (2.0_DP * rranf() - 1.0_DP ) * delta_phi_max 

  ! new polar and azimuthal anlges
  IF ( (ABS(DCOS(theta) + dcostheta) > 1.0_DP) ) THEN
     theta_bound = .true.
     RETURN
  END IF

  theta =  DACOS(DCOS(theta) + dcostheta)
  phi = phi + dphi
  
  ! new coordinates
  atom_list(this_atom,im,is)%rxp = bond_length * DSIN(theta) * DCOS(phi)
  atom_list(this_atom,im,is)%ryp = bond_length * DSIN(theta) * DSIN(phi)
  atom_list(this_atom,im,is)%rzp = bond_length * DCOS(theta)

END SUBROUTINE Change_Phi_Theta
  
