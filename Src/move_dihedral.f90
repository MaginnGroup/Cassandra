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
SUBROUTINE Rigid_Dihedral_Change
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
!   Save_Old_Cartesian_Coordinates
!   Save_Old_Internal_Coordinates
!   Compute_Molecule_Bond_Energy
!   Compute_Molecule_Angle_Energy
!   Compute_Molecule_Dihedral_Energy
!   Compute_Molecule_Improper_Energy
!   Compute_Molecule_Nonbond_Intra_Energy
!   Compute_Molecule_Nonbond_Inter_Energy
!   Update_System_Ewald_Reciprocal_Energy
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
  USE Global_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE File_Names
  USE Energy_Routines
  USE IO_Utilities

  IMPLICIT NONE

!  !$ include 'omp_lib.h'
  
  INTEGER :: ibox, is, im, dihedral_to_move, lm
  INTEGER :: i, j, this_atom, mcstep
  INTEGER :: atom1, atom2, atom3, atom4, iatom1, iatom2, iatom3, iatom4
  INTEGER :: natoms_to_place
  INTEGER, DIMENSION(:), ALLOCATABLE :: atoms_to_place_list
  INTEGER :: nmols_tot, nmols_box(nbr_boxes)

  REAL(DP) :: x_box(nbr_boxes), x_species(nspecies), ln_pacc
  REAL(DP) :: iatom2_rxp, iatom2_ryp, iatom2_rzp, vec23(3), vec21(3)
  REAL(DP) :: perp_vec1(3), perp_vec2(3), aligner(3,3), hanger(3,3)
  REAL(DP) :: phi_trial, cosphi, sinphi, tempx, tempy, tempz

  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper, E_intra_vdw, E_intra_qq, E_inter_vdw, E_inter_qq, E_periodic_qq
  REAL(DP) :: E_bond_move, E_angle_move, E_dihedral_move, E_improper_move, E_intra_vdw_move
  REAL(DP) :: E_intra_qq_move, W_intra_vdw, W_intra_qq, W_inter_vdw, W_inter_qq
  REAL(DP) :: W_intra_vdw_move, W_intra_qq_move, W_inter_vdw_move, W_inter_qq_move
  REAL(DP) :: E_reciprocal_move, E_inter_vdw_move, E_inter_qq_move, delta_e, p_acc
  REAL(DP) :: rand_no, W_reciprocal_move

  REAL(DP), DIMENSION(3,3) :: tvdm, tcdm, qw_di

  LOGICAL ::  inter_overlap, intra_overlap, accept_or_reject

 ! Variables associated with framework simulations

  REAL(DP) :: E_framework, E_framework_move, E_correction_move

  LOGICAL :: framework_overlap

  intra_overlap = .FALSE.
  inter_overlap = .FALSE.
  accept = .FALSE.    

  ! Sum the total number of molecules 
  nmols_tot = 0 ! sum over species, box
  DO ibox = 1, nbr_boxes
    nmols_box(ibox) = 0
    DO is = 1, nspecies
      ! Only count mobile species
      IF ( ndihedrals(is) > 0 ) THEN
        nmols_tot = nmols_tot + nmols(is,ibox)
        nmols_box(ibox) = nmols_box(ibox) + nmols(is,ibox)
      END IF
    END DO
  END DO

  ! If there are no molecules then return
  IF (nmols_tot == 0) THEN
     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,5X,X,3X,X,I3,X,L8,X,9X,X,A9)') &
             i_mcstep, 'dihed' , ibox, accept, 'no mols'
     END IF
     RETURN
  END IF

  ! If needed, choose a box based on its total mol fraction
  IF(nbr_boxes .GT. 1) THEN

    DO ibox = 1, nbr_boxes
       x_box(ibox) = REAL(nmols_box(ibox),DP)/REAL(nmols_tot,DP)
       IF ( ibox > 1 ) THEN
          x_box(ibox) = x_box(ibox) + x_box(ibox-1)
       END IF
    END DO
  
    DO ibox = 1, nbr_boxes
       IF ( rranf() <= x_box(ibox)) EXIT
    END DO

  ELSE

    ibox = 1

  END IF

  ! error check
  IF( nmols_box(ibox) == 0 ) THEN
     err_msg = ''
     err_msg(1) = 'No movable molecules in box ' // TRIM(Int_To_String(ibox))
     CALL Clean_Abort(err_msg, 'Rigid_Dihedral_Change')
  END IF

  ! Choose species based on the mol fraction, using Golden sampling
  DO is = 1, nspecies
     IF (ndihedrals(is) > 0) THEN
        x_species(is) = REAL(nmols(is,ibox), DP)/REAL(nmols_box(ibox),DP)
     ELSE
        x_species(is) = 0.0_DP
     END IF
     IF ( is > 1 ) THEN
        x_species(is) = x_species(is) + x_species(is-1)
     END IF
  END DO


  rand_no = rranf()
  DO is = 1, nspecies
     IF( rand_no <= x_species(is)) EXIT
  END DO

  ! error check
  IF ( ndihedrals(is) == 0 ) THEN
     err_msg = ''
     err_msg(1) = 'Species ' // TRIM(Int_To_String(is)) // ' has no dihedrals'
     CALL Clean_Abort(err_msg, 'Rigid_Dihedral_Change')
  END IF

  ! Choose one of the molecules at random
  im = INT( rranf() * nmols(is,ibox) ) + 1
  ! Get the index of imth molecule of species is in the box.
  lm = locate(im,is,ibox)

  tot_trials(ibox) = tot_trials(ibox) + 1
  ntrials(is,ibox)%dihedral = ntrials(is,ibox)%dihedral + 1

  ! Store old positions and internal cooridinates

  CALL Save_Old_Cartesian_Coordinates(lm,is)
  CALL Save_Old_Internal_Coordinates(lm,is)

  
  ! Compute the bonded interactions before the move
  CALL Compute_Molecule_Bond_Energy(lm,is,E_bond)
  CALL Compute_Molecule_Angle_Energy(lm,is,E_angle)
  CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihedral)
  CALL Compute_Molecule_Improper_Energy(lm,is,E_improper)
  
  ! Compute the nobonded interaction before the move
  CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is,E_intra_vdw,E_intra_qq,E_periodic_qq,intra_overlap)
  IF (l_pair_nrg) THEN
     CALL Store_Molecule_Pair_Interaction_Arrays(lm,is,ibox,E_inter_vdw,E_inter_qq)
  ELSE 
     CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw,E_inter_qq,inter_overlap)
  END IF
  E_inter_qq = E_inter_qq + E_periodic_qq

  IF (inter_overlap)  THEN
     err_msg = ""
     err_msg(1) = "Attempted to change a dihedral of molecule " // TRIM(Int_To_String(im)) // &
                  " of species " // TRIM(Int_To_String(is))
     IF (nbr_boxes > 1) err_msg(1) = err_msg(1) // " in box " // TRIM(Int_To_String(ibox))
     err_msg(2) = "but the molecule energy is too high"
     IF (start_type(ibox) == "make_config" ) THEN
        err_msg(3) = "Try increasing Rcutoff_Low, increasing the box size, or "
        err_msg(4) = "decreasing the initial number of molecules"
     END IF
     CALL Clean_Abort(err_msg, "Rigid_Dihedral_Change")
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

  rand_no = rranf()

  IF ( rand_no < 0.5_DP ) THEN

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

  iatom2_rxp = atom_list(iatom2,lm,is)%rxp
  iatom2_ryp = atom_list(iatom2,lm,is)%ryp
  iatom2_rzp = atom_list(iatom2,lm,is)%rzp

  atom_list(:,lm,is)%rxp = atom_list(:,lm,is)%rxp - iatom2_rxp
  atom_list(:,lm,is)%ryp = atom_list(:,lm,is)%ryp - iatom2_ryp
  atom_list(:,lm,is)%rzp = atom_list(:,lm,is)%rzp - iatom2_rzp

  ! We will generate a perpendicular frame at atom2 such that x - axis
  ! is aligned along iatom2 --- > iatom3 and y - axis is in the plane
  ! defined by iatom1 - iatom2 - iatom3

  vec23(1) = atom_list(iatom3,lm,is)%rxp - atom_list(iatom2,lm,is)%rxp
  vec23(2) = atom_list(iatom3,lm,is)%ryp - atom_list(iatom2,lm,is)%ryp
  vec23(3) = atom_list(iatom3,lm,is)%rzp - atom_list(iatom2,lm,is)%rzp

  ! vector from iatom2 to iatom1

  vec21(1) = atom_list(iatom1,lm,is)%rxp - atom_list(iatom2,lm,is)%rxp
  vec21(2) = atom_list(iatom1,lm,is)%ryp - atom_list(iatom2,lm,is)%ryp
  vec21(3) = atom_list(iatom1,lm,is)%rzp - atom_list(iatom2,lm,is)%rzp

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

     tempx = atom_list(this_atom,lm,is)%rxp
     tempy = atom_list(this_atom,lm,is)%ryp
     tempz = atom_list(this_atom,lm,is)%rzp

     atom_list(this_atom,lm,is)%rxp = tempx * aligner(1,1) + tempy * aligner(1,2) + &
          tempz * aligner(1,3)
     atom_list(this_atom,lm,is)%ryp = tempx * aligner(2,1) + tempy * aligner(2,2) + &
          tempz * aligner(2,3)
     atom_list(this_atom,lm,is)%rzp = tempx * aligner(3,1) + tempy * aligner(3,2) + &
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

     tempy = atom_list(this_atom,lm,is)%ryp
     tempz = atom_list(this_atom,lm,is)%rzp

     ! apply the transformation
     
     atom_list(this_atom,lm,is)%ryp =  cosphi * tempy + sinphi * tempz
     atom_list(this_atom,lm,is)%rzp = -sinphi * tempy + cosphi * tempz

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

     tempx = atom_list(this_atom,lm,is)%rxp
     tempy = atom_list(this_atom,lm,is)%ryp
     tempz = atom_list(this_atom,lm,is)%rzp

     atom_list(this_atom,lm,is)%rxp = tempx * hanger(1,1) + tempy * hanger(1,2) + &
          tempz * hanger(1,3)
     atom_list(this_atom,lm,is)%ryp = tempx * hanger(2,1) + tempy * hanger(2,2) + &
          tempz * hanger(2,3)
     atom_list(this_atom,lm,is)%rzp = tempx * hanger(3,1) + tempy * hanger(3,2) + &
          tempz * hanger(3,3)

  END DO

  atom_list(:,lm,is)%rxp = atom_list(:,lm,is)%rxp + iatom2_rxp
  atom_list(:,lm,is)%ryp = atom_list(:,lm,is)%ryp + iatom2_ryp
  atom_list(:,lm,is)%rzp = atom_list(:,lm,is)%rzp + iatom2_rzp

  delta_e = 0.0_DP
  ! Now compute the energy of this molecule in the new conformation. First compute the intramolecular
  ! nonbonded interactions so that if an overlap is detected, the move can be immediately rejected.

  ! Update the COM and distance of the atom farthest from COM

  ! Check to make sure that the molecule has not moved out of the slit pore


  inter_overlap = .FALSE.
  CALL Get_COM(lm,is)
  CALL Compute_Max_COM_Distance(lm,is)
  
  CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is,E_intra_vdw_move,E_intra_qq_move,E_periodic_qq,intra_overlap)

  IF (intra_overlap) inter_overlap = .TRUE.

  IF ( .NOT. inter_overlap) THEN
     CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw_move,E_inter_qq_move,inter_overlap)
     E_inter_qq_move = E_inter_qq_move + E_periodic_qq
  
  END IF
     
  IF (inter_overlap) THEN
     
     ! Reject the move, reset the old cartesian and internal coordinates

     CALL Revert_Old_Cartesian_Coordinates(lm,is)
     CALL Revert_Old_Internal_Coordinates(lm,is)
     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(lm,is,ibox)

     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,A9)') &
             i_mcstep, 'dihed' , lm, is, ibox, accept, 'overlap'
     END IF
  ELSE
  
     CALL Compute_Molecule_Bond_Energy(lm,is,E_bond_move)
     CALL Compute_Molecule_Angle_Energy(lm,is,E_angle_move)
     CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihedral_move)
     CALL Compute_Molecule_Improper_Energy(lm,is,E_improper_move)


     IF (int_charge_sum_style(ibox) == charge_ewald) THEN
        CALL Update_System_Ewald_Reciprocal_Energy(lm,is,ibox, &
             int_intra,E_reciprocal_move)
        delta_e = E_reciprocal_move - energy(ibox)%ewald_reciprocal

        
     END IF

     ! energy difference

     delta_e = delta_e + (E_bond_move - E_bond) + (E_angle_move - E_angle) &
             + (E_dihedral_move - E_dihedral) + (E_improper_move - E_improper) &
             + (E_intra_vdw_move - E_intra_vdw) &
             + (E_intra_qq_move - E_intra_qq) &
             + (E_inter_vdw_move - E_inter_vdw) &
             + (E_inter_qq_move - E_inter_qq)

     IF ( int_sim_type == sim_nvt_min) THEN
        IF ( delta_e <= 0.0_DP ) THEN
           p_acc = 1.0_DP
        ELSE
           p_acc = 0.0_DP
        END IF

     ELSE
        
        ln_pacc = beta(ibox) * delta_e
        accept = accept_or_reject(ln_pacc)

     END IF

     IF ( accept ) THEN
        ! accept the move and update the energies
        energy(ibox)%intra = energy(ibox)%intra + E_bond_move - E_bond + E_angle_move - E_angle + &
             E_dihedral_move - E_dihedral + E_improper_move - E_improper
        energy(ibox)%bond = energy(ibox)%bond + E_bond_move - E_bond
        energy(ibox)%angle = energy(ibox)%angle + E_angle_move - E_angle
        energy(ibox)%dihedral = energy(ibox)%dihedral + E_dihedral_move - E_dihedral
        energy(ibox)%intra_vdw = energy(ibox)%intra_vdw + E_intra_vdw_move - E_intra_vdw
        energy(ibox)%intra_q   = energy(ibox)%intra_q   + E_intra_qq_move - E_intra_qq
        energy(ibox)%inter_vdw = energy(ibox)%inter_vdw + E_inter_vdw_move - E_inter_vdw
        energy(ibox)%inter_q   = energy(ibox)%inter_q   + E_inter_qq_move - E_inter_qq

        IF (int_charge_sum_style(ibox) == charge_ewald) THEN
           energy(ibox)%ewald_reciprocal = E_reciprocal_move
        END IF

        energy(ibox)%total = energy(ibox)%total + delta_e

        ! update success counter

        nsuccess(is,ibox)%dihedral = nsuccess(is,ibox)%dihedral + 1

        ! Compute the COM positions

        CALL Get_Internal_Coordinates(lm,is)
        IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
     ELSE 

        ! Reject the move and revert the old coordinates

        CALL Revert_Old_Cartesian_Coordinates(lm,is)
        CALL Revert_Old_Internal_Coordinates(lm,is)

        IF (int_charge_sum_style(ibox) == charge_ewald) THEN
           ! Also reset the old cos_sum and sin_sum for reciprocal space vectors
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
           cos_sum(:,ibox) = cos_sum_old(:,ibox)
           sin_sum(:,ibox) = sin_sum_old(:,ibox)
           !$OMP END PARALLEL WORKSHARE

        END IF
        IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(lm,is,ibox)
     END IF

     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,F9.3)') &
             i_mcstep, 'dihed' , lm, is, ibox, accept, ln_pacc
     END IF

  END IF

!  DEALLOCATE(atoms_to_place_list)

END SUBROUTINE Rigid_Dihedral_Change

  
