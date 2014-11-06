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
!********************************************************************************

!********************************************************************************
SUBROUTINE Rigid_Dihedral_Change(this_box)
!********************************************************************************

!********************************************************************************
! The subroutine performs a rigid body dihedral angle move. It randomly picks
! a molecule and chooses a dihedral angle to be perturbed. 
!
! The algorithm is as follows.:
!
!  1. Choose a species at random,
!  2. Determine if it contains any dihedral, if not return to 1.
!  3. Select a dihedral angle randomly.
!  4. Choose randomly one of the ends of the dihedral to move.
!  5. Perform rigid body rotation for the atoms bonded to the atom selected in 4.
!  6. Accept or reject the move based on Metropolis criterion
!
!
! Called by
!
!   NVTMC_Driver
!   NPTMC_Driver
!
! Calls
!
!   Get_Nmolecules_Species
!   Get_Index_Molecules
!   Save_Old_Cartesian_Coordinates
!   Save_Old_Internal_Coordinates
!   Compute_Molecule_Bond_Energy
!   Compute_Molecule_Angle_Energy
!   Compute_Molecule_Dihedral_Energy
!   Compute_Molecule_Improper_Energy
!   Compute_Molecule_Nonbond_Intra_Energy
!   Compute_Molecule_Nonbond_Inter_Energy
!   Compute_Ewald_Reciprocal_Energy_Difference
!   Get_COM
!   Get_Internal_Coordinates
!   Revert_Old_Cartesian_Coordinates
!   Revert_Old_Internal_Coordinates
!
! Revision History
!
!  12/10/13 : Beta release
!
!*********************************************************************************

  USE Type_Definitions
  USE Run_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE File_Names
  USE Energy_Routines

  IMPLICIT NONE

!  !$ include 'omp_lib.h'
  
  INTEGER :: this_box, is, nmolecules_species, im, dihedral_to_move, alive
  INTEGER :: i, j, this_atom
  INTEGER :: atom1, atom2, atom3, atom4, iatom1, iatom2, iatom3, iatom4
  INTEGER :: natoms_to_place
  INTEGER, DIMENSION(:), ALLOCATABLE :: atoms_to_place_list

  REAL(DP) :: iatom2_rxp, iatom2_ryp, iatom2_rzp, vec23(3), vec21(3)
  REAL(DP) :: perp_vec1(3), perp_vec2(3), aligner(3,3), hanger(3,3)
  REAL(DP) :: phi_trial, cosphi, sinphi, tempx, tempy, tempz

  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper, E_intra_vdw, E_intra_qq, E_inter_vdw, E_inter_qq
  REAL(DP) :: E_bond_move, E_angle_move, E_dihedral_move, E_improper_move, E_intra_vdw_move
  REAL(DP) :: E_intra_qq_move, W_intra_vdw, W_intra_qq, W_inter_vdw, W_inter_qq
  REAL(DP) :: W_intra_vdw_move, W_intra_qq_move, W_inter_vdw_move, W_inter_qq_move
  REAL(DP) :: E_reciprocal_move, E_inter_vdw_move, E_inter_qq_move, delta_e, p_acc
  REAL(DP) :: rand, W_reciprocal_move

  REAL(DP), DIMENSION(3,3) :: tvdm, tcdm, qw_di

  LOGICAL ::  inter_overlap, intra_overlap

 ! Variables associated with framework simulations

  REAL(DP) :: E_framework, E_framework_move, E_correction_move

  LOGICAL :: framework_overlap

  

  ! choose a box at random

  this_box = INT( rranf() * nbr_boxes ) + 1
  tot_trials(this_box) = tot_trials(this_box) + 1

  ! choose a species at random

  DO WHILE (.true.) 

     is = INT( rranf() * nspecies ) + 1

     ! check to make sure that it contains at least one dihedral for perturbation

     IF ( ndihedrals(is) == 0 ) CYCLE

     ! Obtain number of molecules of this species

     CALL Get_Nmolecules_Species(this_box,is, nmolecules_species)

     IF (nmolecules_species /= 0) EXIT

  END DO

  ! Choose one of the molecules at random

  im = INT( rranf() * nmolecules_species ) + 1

  ! Get the index of imth molecule of species is in the box.

  CALL Get_Index_Molecule(this_box,is,im,alive)

  ntrials(is,this_box)%dihedral = ntrials(is,this_box)%dihedral + 1

  ! Store old positions and internal cooridinates

  CALL Save_Old_Cartesian_Coordinates(alive,is)
  CALL Save_Old_Internal_Coordinates(alive,is)

  ! Compute the bonded interactions before the move
  CALL Compute_Molecule_Bond_Energy(alive,is,E_bond)
  CALL Compute_Molecule_Angle_Energy(alive,is,E_angle)
  CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihedral)
  CALL Compute_Molecule_Improper_Energy(alive,is,E_improper)
  
  ! Compute the nobonded interaction before the move
  CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw,E_intra_qq,intra_overlap)
  CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_inter_vdw,E_inter_qq,inter_overlap)
  ! compute the energy related to the framework

  IF (inter_overlap) THEN
     WRITE(*,*) 'Disaster, overlap in the old configruation'
     WRITE(*,*) 'rigid_dihedral_change.f90'
     WRITE(*,*) alive, is, this_box
     WRITE(*,*) 'inter overlap', inter_overlap
     WRITE(*,*) 'Framework overlap', framework_overlap
  END IF
     


 ! Select a dihedral at random 

  dihedral_to_move = INT( rranf() * ndihedrals(is) ) + 1

  ! Determine the atoms that form the dihedral to move
  
  atom1 = dihedral_list(dihedral_to_move,is)%atom1
  atom2 = dihedral_list(dihedral_to_move,is)%atom2
  atom3 = dihedral_list(dihedral_to_move,is)%atom3
  atom4 = dihedral_list(dihedral_to_move,is)%atom4

  
  IF (.NOT. ALLOCATED(atoms_to_place_list)) ALLOCATE(atoms_to_place_list(MAXVAL(natoms)),Stat=AllocateStatus)
  IF (AllocateStatus /=0) THEN
     err_msg =''
     err_msg(1) = 'Memory could not be allocated to atoms_to_place_list'
     CALL Clean_Abort(err_msg,'Rigid_Dihedral_Change')
  END IF
  atoms_to_place_list = 0
  ! Choose one of the ends to move at random and relable the atoms

  rand = rranf()

  IF ( rand < 0.5_DP ) THEN

     iatom1 = atom1
     iatom2 = atom2
     iatom3 = atom3
     iatom4 = atom4

     natoms_to_place = dihedral_atoms_to_place_list(dihedral_to_move,is)%atom4_natoms

     Do j = 1, natoms_to_place
        atoms_to_place_list(j) = dihedral_atoms_to_place_list(dihedral_to_move,is)%atom4(j)
     END Do

  ELSE

     iatom1 = atom4
     iatom2 = atom3
     iatom3 = atom2
     iatom4 = atom1

     natoms_to_place = dihedral_atoms_to_place_list(dihedral_to_move,is)%atom1_natoms
  
     DO j = 1, natoms_to_place
        atoms_to_place_list(j) = dihedral_atoms_to_place_list(dihedral_to_move,is)%atom1(j)
     END DO
     
  END IF

  ! Now we will align the normal to the angle iatom2, iatom3, iatom4
  ! such that these atoms are in the xy plane.

  ! Move all the atoms with respect to iatom2

  iatom2_rxp = atom_list(iatom2,alive,is)%rxp
  iatom2_ryp = atom_list(iatom2,alive,is)%ryp
  iatom2_rzp = atom_list(iatom2,alive,is)%rzp

  atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp - iatom2_rxp
  atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp - iatom2_ryp
  atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp - iatom2_rzp

  ! We will generate a perpendicular frame at atom2 such that x - axis
  ! is aligned along iatom2 --- > iatom3 and y - axis is in the plane
  ! defined by iatom1 - iatom2 - iatom3

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

  ! form a matrix composed of these unit vectors

  DO j = 1, 3
     aligner(1,j) = vec23(j)
     aligner(2,j) = perp_vec2(j)
     aligner(3,j) = perp_vec1(j)
  END DO

  
  ! use the aligner matrix to align the plane to the xy plane
  ! apply this transformation only to the atoms that are presently
  ! involved in the dihedral and that will move as a result of the
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

  !-- at this point iatom1, iatom2, iatom3 must be in xy plane. Now we will apply a rotation
  !-- around the axis atom2-atom3. Choose a dihedral angle randomly from a phi-max


  phi_trial = 2.0_DP * PI * rranf()
 
  !!!phi_trial = ( 2.0_DP * rranf() - 1.0_DP) * species_list(is)%max_torsion
 
  cosphi = DCOS(phi_trial)
  sinphi = DSIN(phi_trial)

  ! Rotate the vectors of the atoms that move due to dihedral angle change by phi_trial

  DO j = 1, natoms_to_place

     this_atom = atoms_to_place_list(j)

     tempy = atom_list(this_atom,alive,is)%ryp
     tempz = atom_list(this_atom,alive,is)%rzp

     ! apply the transformation
     
     atom_list(this_atom,alive,is)%ryp =  cosphi * tempy + sinphi * tempz
     atom_list(this_atom,alive,is)%rzp = -sinphi * tempy + cosphi * tempz

  END DO

  
  ! Now revert to the space fixed coordinates in reverse order

  ! To do this first align the perpendicular axes to the space fixed axes and then
  ! shift the origin from atom2 to (0,0,0)

  ! Form the hanger matrix

  DO j = 1, 3
     DO i = 1, 3
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

  delta_e = 0.0_DP
  ! Now compute the energy of this molecule in the new conformation. First compute the intramolecular
  ! nonbonded interactions so that if an overlap is detected, the move can be immediately rejected.

  ! Update the COM and distance of the atom farthest from COM

  ! Check to make sure that the molecule has not moved out of the slit pore


  inter_overlap = .FALSE.
     CALL Get_COM(alive,is)
     CALL Compute_Max_COM_Distance(alive,is)
     
     CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw_move,E_intra_qq_move,intra_overlap)

     IF (intra_overlap) inter_overlap = .TRUE.

     IF ( .NOT. inter_overlap) THEN
        CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_inter_vdw_move,E_inter_qq_move,inter_overlap)
     
     END IF
     
  IF (inter_overlap) THEN
     
     ! Reject the move, reset the old cartesian and internal coordinates

     CALL Revert_Old_Cartesian_Coordinates(alive,is)
     CALL Revert_Old_Internal_Coordinates(alive,is)

  ELSE
  
     CALL Compute_Molecule_Bond_Energy(alive,is,E_bond_move)
     CALL Compute_Molecule_Angle_Energy(alive,is,E_angle_move)
     CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihedral_move)
     CALL Compute_Molecule_Improper_Energy(alive,is,E_improper_move)


     IF (int_charge_sum_style(this_box) == charge_ewald) THEN
        CALL Compute_Ewald_Reciprocal_Energy_Difference(alive,alive,is,this_box,int_intra,E_reciprocal_move)
        delta_e = E_reciprocal_move

        
     END IF

     ! energy difference

     delta_e = E_bond_move - E_bond + E_angle_move - E_angle + E_dihedral_move - E_dihedral + &
          E_improper_move - E_improper + E_intra_vdw_move - E_intra_vdw + E_intra_qq_move - E_intra_qq + &
          E_inter_vdw_move - E_inter_vdw + E_inter_qq_move - E_inter_qq + delta_e


     IF ( int_sim_type == sim_nvt_min) THEN
        IF ( delta_e <= 0.0_DP ) THEN
           p_acc = 1.0_DP
        ELSE
           p_acc = 0.0_DP
        END IF

     ELSE
        
        p_acc = MIN(1.0_DP,DEXP(-beta(this_box)*delta_e))

     END IF

     IF ( rranf() < p_acc ) THEN
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

        ! update success counter

        nsuccess(is,this_box)%dihedral = nsuccess(is,this_box)%dihedral + 1

        ! Compute the COM positions

        CALL Get_Internal_Coordinates(alive,is)

     ELSE 

        ! Reject the move and revert the old coordinates

        CALL Revert_Old_Cartesian_Coordinates(alive,is)
        CALL Revert_Old_Internal_Coordinates(alive,is)

        IF (int_charge_sum_style(this_box) == charge_ewald) THEN
           ! Also reset the old cos_sum and sin_sum for reciprocal space vectors
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
           cos_sum(:,this_box) = cos_sum_old(:,this_box)
           sin_sum(:,this_box) = sin_sum_old(:,this_box)
           !$OMP END PARALLEL WORKSHARE

        END IF

     END IF

  END IF

!  DEALLOCATE(atoms_to_place_list)

END SUBROUTINE Rigid_Dihedral_Change

  
