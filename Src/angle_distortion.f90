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

SUBROUTINE Angle_Distortion(this_box)

!*******************************************************************************
!
! This subroutine attempts to change an angle of a randomly selected molecule
!
! Called by
!
! NVTMC_driver
! NPTMC_driver
! main
! gcmc_driver
! gemc_driver
!
! Revision History
!
! 12/10/13  : Beta Release
!********************************************************************************
  USE Type_Definitions
  USE Run_Variables
  USE Angle_Dist_Pick
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE Pair_Nrg_Routines


  IMPLICIT NONE

!  !$ include 'omp_lib.h'

  INTEGER :: this_box, is, im, alive, nmolecules_species
  INTEGER :: angle_to_move, atom1, atom2, atom3, iatom1, iatom2, iatom3
  INTEGER :: natoms_to_place, this_atom, i, j

  INTEGER, DIMENSION(:), ALLOCATABLE :: atoms_to_place_list

  REAL(DP) :: theta_0, theta_new, delta_theta, prob_0, prob_new
  REAL(DP) :: iatom2_rxp, iatom2_ryp, iatom2_rzp, vec21(3), vec23(3)
  REAL(DP) :: perp_vec1(3), perp_vec2(3), aligner(3,3), hanger(3,3)
  REAL(DP) :: tempx, tempy, tempz, cos_dtheta, sin_dtheta
  REAL(DP) :: rand, p_acc
  
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper, E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq
  REAL(DP) :: E_bond_move, E_angle_move, E_dihedral_move, E_improper_move, E_intra_vdw_move
  REAL(DP) :: E_intra_qq_move, E_inter_vdw_move, E_inter_qq_move
  REAL(DP) :: E_reciprocal_move
  REAL(DP) :: delta_e

  LOGICAL ::  inter_overlap,intra_overlap

  ! Pair Energy variables

  REAL(DP), ALLOCATABLE :: cos_mol_old(:),sin_mol_old(:)
  INTEGER :: position

  intra_overlap = .false.

  ! Choose a box at random

  this_box = INT ( rranf() * nbr_boxes ) + 1

  tot_trials(this_box) = tot_trials(this_box) + 1
  
  ! Obtain number of molecules of each of the species

  DO WHILE (.true.)
     
     is = INT( rranf() * nspecies ) + 1

     ! Check to make sure that it contains at least one angle that can be perturbed

     IF ( nangles(is) == 0 ) CYCLE

     ! Choose an angle to move
     
     angle_to_move = INT ( rranf() * nangles(is) ) + 1

     ! check to see if the angle type is fixed

     IF ( angle_list(angle_to_move,is)%angle_potential_type == 'fixed') CYCLE

     ! Obtain number of molecules of this species

     CALL Get_Nmolecules_Species(this_box,is,nmolecules_species)

     ! Make sure that molecules of the selected species exist in the simulation box

     IF ( nmolecules_species /= 0 ) EXIT

  END DO

  ! Select a molecule at random for the move

  im = INT( rranf() * nmolecules_species ) + 1

  ! Get the index of imth molecule of species is in the box.

  CALL Get_Index_Molecule(this_box,is,im,alive)

  ntrials(is,this_box)%angle = ntrials(is,this_box)%angle + 1

  ! Compute the energy of the molecule before the move
  CALL Compute_Molecule_Bond_Energy(alive,is,E_bond)
  CALL Compute_Molecule_Angle_Energy(alive,is,E_angle)
  CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihedral)
  CALL Compute_Molecule_Improper_Energy(alive,is,E_improper)

  CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw,E_intra_qq,intra_overlap)
  IF (l_pair_nrg) THEN
     CALL Store_Molecule_Pair_Interaction_Arrays(alive,is,this_box,E_inter_vdw,E_inter_qq)
  ELSE
     CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_inter_vdw,E_inter_qq,inter_overlap)
  END IF

  IF (inter_overlap)  THEN
     WRITE(*,*) 'Disaster, overlap in the old configruation'
     WRITE(*,*) 'Translate'
     WRITE(*,*) alive, is, this_box
     WRITE(*,*) 'inter overlap', inter_overlap
  END IF
  ! Store the old positions of the atoms 

  CALL Save_Old_Cartesian_Coordinates(alive,is)
  CALL Save_Old_Internal_Coordinates(alive,is)

  
  ! determine the atoms that define the angle

  atom1 = angle_list(angle_to_move,is)%atom1
  atom2 = angle_list(angle_to_move,is)%atom2
  atom3 = angle_list(angle_to_move,is)%atom3


  ! We first determine the angle before the move and also the probability of observing this
  ! angle. 

  theta_0 = internal_coord_list(angle_to_move,im,is)%bond_angle_radians
  ! obtain the probability associated with this angle

  CALL Get_Theta_Prob(angle_to_move,is,theta_0,prob_0)

  CALL Pick_Angle(angle_to_move,is,theta_new,prob_new)
  
  ! Next step is to determine the position of atom_to_move such that bond lengths
  ! are not perturbed. We will implement the algorithm in which bond angle distortion
  ! is performed in the plane of three atoms atom1,atom2 and atom3 that define the angle
  ! of interest. This is equivalent to rotating one arm of the angle about a perpendicular
  ! axis passing through the plane of the angle. 

  ! Choose one of the ends to move at random and relabel the atoms

  ALLOCATE(atoms_to_place_list(MAXVAL(natoms)))
  atoms_to_place_list = 0

  rand = rranf()

  IF ( rand < 0.5_DP ) THEN

     ! atom1 is moving

     iatom1 = atom1
     iatom2 = atom2
     iatom3 = atom3

     natoms_to_place = angle_atoms_to_place_list(angle_to_move,is)%atom1_natoms

     DO j = 1, natoms_to_place
        atoms_to_place_list(j) = angle_atoms_to_place_list(angle_to_move,is)%atom1(j)
     END DO

  ELSE

     ! atom3 is moving

     iatom1 = atom3
     iatom2 = atom2
     iatom3 = atom1

     natoms_to_place = angle_atoms_to_place_list(angle_to_move,is)%atom3_natoms

     DO j = 1, natoms_to_place
        atoms_to_place_list(j) = angle_atoms_to_place_list(angle_to_move,is)%atom3(j)
     END DO

  END IF
  
  ! Move all the atoms with respect to iatom2
  
  iatom2_rxp = atom_list(iatom2,alive,is)%rxp
  iatom2_ryp = atom_list(iatom2,alive,is)%ryp
  iatom2_rzp = atom_list(iatom2,alive,is)%rzp
  
  atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp - iatom2_rxp
  atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp - iatom2_ryp
  atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp - iatom2_rzp
  
  ! We will generate a perpendicular frame at atom2 such that x - axis
  ! is aligned along iatom2 --- > iatom3 and y - axis is in the plane
  ! defined by iatom1 - iatom2 - iatom3. We will assume that iatom1 is moving.
  
  vec23(1) = atom_list(iatom3,alive,is)%rxp - atom_list(iatom2,alive,is)%rxp
  vec23(2) = atom_list(iatom3,alive,is)%ryp - atom_list(iatom2,alive,is)%ryp
  vec23(3) = atom_list(iatom3,alive,is)%rzp - atom_list(iatom2,alive,is)%rzp

  ! vector from iatom2 to iatom1

  vec21(1) = atom_list(iatom1,alive,is)%rxp - atom_list(iatom2,alive,is)%rxp
  vec21(2) = atom_list(iatom1,alive,is)%ryp - atom_list(iatom2,alive,is)%ryp
  vec21(3) = atom_list(iatom1,alive,is)%rzp - atom_list(iatom2,alive,is)%rzp

  ! Normalize these vectors

  vec23(:) = vec23(:)/DSQRT(DOT_PRODUCT(vec23,vec23))
  vec21(:) = vec21(:)/DSQRT(DOT_PRODUCT(vec21,vec21))

  ! Now cross these vectors to obtain a vector normal to the plane iatom1-
  ! iatom2, iatom3

  perp_vec1(1) =  vec23(2) * vec21(3) - vec23(3) * vec21(2)
  perp_vec1(2) = -vec23(1) * vec21(3) + vec23(3) * vec21(1)
  perp_vec1(3) =  vec23(1) * vec21(2) - vec23(2) * vec21(1)

  perp_vec1(:) = perp_vec1(:)/DSQRT(DOT_PRODUCT(perp_vec1,perp_vec1))
  
  ! Cross vec23 and perp_vec1 to get the third orthogonal vector

  perp_vec2(1) =  vec23(2) * perp_vec1(3) - vec23(3) * perp_vec1(2)
  perp_vec2(2) = -vec23(1) * perp_vec1(3) + vec23(3) * perp_vec1(1)
  perp_vec2(3) =  vec23(1) * perp_vec1(2) - vec23(2) * perp_vec1(1)

  perp_vec2(:) = perp_vec2(:)/DSQRT(DOT_PRODUCT(perp_vec2,perp_vec2))

  ! form a matrix composed of these unit vectors. Note that the 
  ! second row of the matrix is composed of perp_vec2 as it is the second
  ! axis that is in the plane of atom1,atom2 and atom3 along with vec23

  DO j = 1, 3
     aligner(1,j) = vec23(j)
     aligner(2,j) = perp_vec2(j)
     aligner(3,j) = perp_vec1(j)
  END DO

  ! use the aligner matrix to align the plane to the xy plane
  ! apply this transformation only to the atoms that are presently
  ! involved in the angle and that will move as a result of the
  ! change in dihedral

  DO j = 1, natoms_to_place
     
     this_atom = atoms_to_place_list(j)

     tempx = atom_list(this_atom,alive,is)%rxp
     tempy = atom_list(this_atom,alive,is)%ryp
     tempz = atom_list(this_atom,alive,is)%rzp

     atom_list(this_atom,alive,is)%rxp = tempx * aligner(1,1) + tempy * aligner(1,2) + &
          tempz * aligner(1,3)
     atom_list(this_atom,alive,is)%ryp = tempx * aligner(2,1) + tempy * aligner(2,2) + &
          tempz * aligner(2,3)
     atom_list(this_atom,alive,is)%rzp = tempx * aligner(3,1) + tempy * aligner(3,2) + &
          tempz * aligner(3,3)

  END DO

  ! Now we will rotate the bond formed by iatom2-iatom1 by theta_new - theta_0

  delta_theta = theta_new - theta_0

  ! This rotation takes place about z axis. If the new angle is greater than the old
  ! one that means, we need to perform rotation in the clockwise direction when the position
  ! of atom1 due to above transformation is in the 1st or 2nd quadrant. While the rotation
  ! is to be performed in the counter clockwise direction if the atom1 is positioned
  ! in the 3rd or 4th quadrant.

  IF ( atom_list(iatom1,alive,is)%ryp >= 0.0_DP ) THEN
     delta_theta = -delta_theta
  END IF

  
  cos_dtheta = DCOS(delta_theta)
  sin_dtheta = DSIN(delta_theta)

  DO j = 1, natoms_to_place

     this_atom = atoms_to_place_list(j)

     tempx = atom_list(this_atom,alive,is)%rxp
     tempy = atom_list(this_atom,alive,is)%ryp
     
     atom_list(this_atom,alive,is)%rxp = tempx * cos_dtheta + tempy * sin_dtheta
     atom_list(this_atom,alive,is)%ryp = -tempx * sin_dtheta + tempy * cos_dtheta

  END DO

  ! at this point the angle atom1-atom2-atom3 must be theta_new. This is something that
  ! needs to be checked.


  ! Apply the reverse transformation to obtain the positions of moving atoms
  ! with respect to the space fixed axes

  ! Form the hanger matrix

  DO i = 1, 3
     DO j = 1, 3
        hanger(i,j) = aligner(j,i)
     END DO
  END DO

  
  DO j = 1, natoms_to_place

     this_atom = atoms_to_place_list(j)

     tempx = atom_list(this_atom,alive,is)%rxp
     tempy = atom_list(this_atom,alive,is)%ryp
     tempz = atom_list(this_atom,alive,is)%rzp

     atom_list(this_atom,alive,is)%rxp = tempx * hanger(1,1) + tempy * hanger(1,2) + &
          tempz * hanger(1,3)
     atom_list(this_atom,alive,is)%ryp = tempx * hanger(2,1) + tempy * hanger(2,2) + &
          tempz * hanger(2,3)
     atom_list(this_atom,alive,is)%rzp = tempx * hanger(3,1) + tempy * hanger(3,2) + &
          tempz * hanger(3,3)

  END DO

  atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp + iatom2_rxp
  atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp + iatom2_ryp
  atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp + iatom2_rzp

  ! Calculate the energies after the move. First compute intramolecular and intermolecular
  ! nonbonded interactions so that the move can be immediately rejected if an overlap is detected.
  ! Since COM cutoff has been enabled. Compute the new COM of the molecule. Note that if the 
  ! move is rejected then automatically revert_cartesian_coordinates will reset the max_dcom

  CALL Get_COM(alive,is)
  CALL Compute_Max_COM_Distance(alive,is)

  CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw_move,E_intra_qq_move,intra_overlap)
  CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_inter_vdw_move,E_inter_qq_move,inter_overlap)


  IF (inter_overlap) THEN
     
     CALL Revert_Old_Cartesian_Coordinates(alive,is)
     CALL Revert_Old_Internal_Coordinates(alive,is)
     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,this_box)
     ! Note that there is no need to reset the sin and cos terms for reciprocal space Ewald
     ! as these energies have not been yet computed.

  ELSE

     delta_e = 0.0_DP

     CALL Compute_Molecule_Bond_Energy(alive,is,E_bond_move)
     CALL Compute_Molecule_Angle_Energy(alive,is,E_angle_move)
     CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihedral_move)
     CALL Compute_Molecule_Improper_Energy(alive,is,E_improper_move)


     
     IF ( int_charge_sum_style(this_box) == charge_ewald .and. has_charge(is)) THEN
        
        ALLOCATE(cos_mol_old(nvecs(this_box)),sin_mol_old(nvecs(this_box)))
        CALL Get_Position_Alive(alive,is,position)
        
        cos_mol_old(:) = cos_mol(1:nvecs(this_box),position)
        sin_mol_old(:) = sin_mol(1:nvecs(this_box),position)
        
        CALL Compute_Ewald_Reciprocal_Energy_Difference(alive,alive,is,this_box,int_intra,E_reciprocal_move)
        delta_e = E_reciprocal_move
            
     END IF
     
     ! energy difference

     delta_e = E_bond_move - E_bond + E_angle_move - E_angle + E_dihedral_move - E_dihedral + &
          E_improper_move - E_improper + E_intra_vdw_move - E_intra_vdw + E_intra_qq_move - E_intra_qq + &
          E_inter_vdw_move - E_inter_vdw + E_inter_qq_move - E_inter_qq + delta_e

     IF ( int_sim_type == sim_nvt_min ) THEN
        IF ( delta_e <= 0.0_DP ) THEN
           p_acc = 1.0_DP
        ELSE
           p_acc = 0.0_DP
        END IF
        
     ELSE
        p_acc = MIN(1.0_DP,(prob_0/prob_new)*DEXP(-beta(this_box)*delta_e))

     END IF
  
     IF ( rranf() < p_acc ) THEN
        !     write(*,*) 'accepted', theta_new
        ! accept the move and update the energies
        energy(this_box)%intra = energy(this_box)%intra + E_bond_move - E_bond + E_angle_move - E_angle + &
             E_dihedral_move - E_dihedral + E_improper_move - E_improper
        energy(this_box)%bond = energy(this_box)%bond + E_bond_move - E_bond
        energy(this_box)%angle = energy(this_box)%angle + E_angle_move - E_angle
        energy(this_box)%dihedral = energy(this_box)%dihedral + E_dihedral_move - E_dihedral
        energy(this_box)%intra_vdw = energy(this_box)%intra_vdw + E_intra_vdw_move - E_intra_vdw
        energy(this_box)%intra_q   = energy(this_box)%intra_q   + E_intra_qq_move - E_intra_qq
        energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + E_inter_vdw_move - E_inter_vdw
        energy(this_box)%inter_q   = energy(this_box)%inter_q   + E_inter_qq_move - E_inter_qq

        IF (int_charge_sum_style(this_box) == charge_ewald) THEN
           energy(this_box)%ewald_reciprocal = energy(this_box)%ewald_reciprocal + E_reciprocal_move

        END IF

        energy(this_box)%total = energy(this_box)%total + delta_e

        nsuccess(is,this_box)%angle = nsuccess(is,this_box)%angle + 1

        ! Update the virial as well
        ! Note that COM and max_dcom have already been updated

        CALL Get_Internal_Coordinates(alive,is)

        ! Fold the molecule

        CALL Fold_Molecule(alive,is,this_box)
        
        IF(ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
        IF(ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)
        
        IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
        
     ELSE 

        ! Reject the move and revert the old coordinates

        CALL Revert_Old_Cartesian_Coordinates(alive,is)
        CALL Revert_Old_Internal_Coordinates(alive,is)

        IF(l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,this_box)
        
        IF (int_charge_sum_style(this_box) == charge_ewald .AND. has_charge(is)) THEN
           ! Also reset the old cos_sum and sin_sum for reciprocal space vectors
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
           cos_sum(:,this_box) = cos_sum_old(:,this_box)
           sin_sum(:,this_box) = sin_sum_old(:,this_box)
           
           cos_mol(1:nvecs(this_box),position) = cos_mol_old(:)
           sin_mol(1:nvecs(this_box),position) = sin_mol_old(:)
           !$OMP END PARALLEL WORKSHARE
           
           DEALLOCATE(cos_mol_old,sin_mol_old)

        END IF
        
     END IF
  END IF
     
  IF(ALLOCATED(atoms_to_place_list)) DEALLOCATE(atoms_to_place_list)

END SUBROUTINE Angle_Distortion
