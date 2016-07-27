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

SUBROUTINE Angle_Distortion

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
  USE Global_Variables
  USE Angle_Dist_Pick
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE Pair_Nrg_Routines
  USE IO_Utilities


  IMPLICIT NONE

!  !$ include 'omp_lib.h'

  INTEGER :: ibox, is, im, lm, nmolecules_species
  INTEGER :: angle_to_move, atom1, atom2, atom3, iatom1, iatom2, iatom3
  INTEGER :: natoms_to_place, this_atom, i, j, mcstep

  INTEGER, DIMENSION(:), ALLOCATABLE :: atoms_to_place_list
  INTEGER :: nmols_tot, nmols_box(nbr_boxes)

  REAL(DP) :: x_box(nbr_boxes), x_species(nspecies), ln_pacc
  REAL(DP) :: theta_0, theta_new, delta_theta, prob_0, prob_new
  REAL(DP) :: iatom2_rxp, iatom2_ryp, iatom2_rzp, vec21(3), vec23(3)
  REAL(DP) :: perp_vec1(3), perp_vec2(3), aligner(3,3), hanger(3,3)
  REAL(DP) :: tempx, tempy, tempz, cos_dtheta, sin_dtheta
  REAL(DP) :: rand_no, p_acc
  
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper, E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq, E_periodic_qq
  REAL(DP) :: E_bond_move, E_angle_move, E_dihedral_move, E_improper_move, E_intra_vdw_move
  REAL(DP) :: E_intra_qq_move, E_inter_vdw_move, E_inter_qq_move
  REAL(DP) :: E_reciprocal_move
  REAL(DP) :: delta_e

  LOGICAL ::  inter_overlap,intra_overlap, accept_or_reject

  ! Pair Energy variables

  REAL(DP), ALLOCATABLE :: cos_mol_old(:),sin_mol_old(:)
  INTEGER :: position

  intra_overlap = .false.
  accept = .FALSE.

  ! Sum the total number of molecules 
  nmols_tot = 0
  DO ibox = 1, nbr_boxes
     nmols_box(ibox) = 0
     DO is = 1, nspecies
        IF (nangles(is) > nangles_fixed(is)) THEN
           nmols_tot = nmols_tot + nmols(is,ibox)
           nmols_box(ibox) = nmols_box(ibox) + nmols(is,ibox)
        END IF
     END Do
  END DO

  ! If there are no molecules then return
  IF (nmols_tot == 0) THEN
     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,F9.3)') &
             i_mcstep, 'angle' , lm, is, ibox, accept, 'no mols'
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
     err_msg(1) = 'No molecules with harmonic angles in box ' // TRIM(Int_To_String(ibox))
     CALL Clean_Abort(err_msg, 'Angle_Distortion')
  END IF

  ! Choose species based on the mol fraction, using Golden sampling
  DO is = 1, nspecies
     IF (nangles(is) > nangles_fixed(is) ) THEN
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
  IF( nangles(is) <= nangles_fixed(is) ) THEN
     err_msg = ''
     err_msg(1) = 'Species ' // TRIM(Int_To_String(is)) // ' does not have harmonic angles'
     CALL Clean_Abort(err_msg, 'Angle_Distortion')
  END IF

  tot_trials(ibox) = tot_trials(ibox) + 1
  
  ! Select a molecule at random for the move
  im = INT( rranf() * nmols(is,ibox) ) + 1
  ! Get the index of imth molecule of species is in the box.
  lm = locate(im,is,ibox)

  ntrials(is,ibox)%angle = ntrials(is,ibox)%angle + 1

  ! Compute the energy of the molecule before the move
  CALL Compute_Molecule_Bond_Energy(lm,is,E_bond)
  CALL Compute_Molecule_Angle_Energy(lm,is,E_angle)
  CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihedral)
  CALL Compute_Molecule_Improper_Energy(lm,is,E_improper)

  CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is,E_intra_vdw,E_intra_qq,E_periodic_qq,intra_overlap)
  IF (l_pair_nrg) THEN
     CALL Store_Molecule_Pair_Interaction_Arrays(lm,is,ibox,E_inter_vdw,E_inter_qq)
  ELSE
     CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw,E_inter_qq,inter_overlap)
  END IF
  E_inter_qq = E_inter_qq + E_periodic_qq

  IF (inter_overlap)  THEN
     err_msg = ""
     err_msg(1) = "Attempted to change an angle of molecule " // TRIM(Int_To_String(im)) // &
                  " of species " // TRIM(Int_To_String(is))
     IF (nbr_boxes > 1) err_msg(1) = err_msg(1) // " in box " // TRIM(Int_To_String(ibox))
     err_msg(2) = "but the molecule energy is too high"
     IF (start_type(ibox) == "make_config" ) THEN
        err_msg(3) = "Try increasing Rcutoff_Low, increasing the box size, or "
        err_msg(4) = "decreasing the initial number of molecules"
     END IF
     CALL Clean_Abort(err_msg, "Angle_Distortion")
  END IF

  ! Store the old positions of the atoms 
  CALL Save_Old_Cartesian_Coordinates(lm,is)
  CALL Save_Old_Internal_Coordinates(lm,is)

  
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

  rand_no = rranf()

  IF ( rand_no < 0.5_DP ) THEN

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
  
  iatom2_rxp = atom_list(iatom2,lm,is)%rxp
  iatom2_ryp = atom_list(iatom2,lm,is)%ryp
  iatom2_rzp = atom_list(iatom2,lm,is)%rzp
  
  atom_list(:,lm,is)%rxp = atom_list(:,lm,is)%rxp - iatom2_rxp
  atom_list(:,lm,is)%ryp = atom_list(:,lm,is)%ryp - iatom2_ryp
  atom_list(:,lm,is)%rzp = atom_list(:,lm,is)%rzp - iatom2_rzp
  
  ! We will generate a perpendicular frame at atom2 such that x - axis
  ! is aligned along iatom2 --- > iatom3 and y - axis is in the plane
  ! defined by iatom1 - iatom2 - iatom3. We will assume that iatom1 is moving.
  
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

  ! Now we will rotate the bond formed by iatom2-iatom1 by theta_new - theta_0

  delta_theta = theta_new - theta_0

  ! This rotation takes place about z axis. If the new angle is greater than the old
  ! one that means, we need to perform rotation in the clockwise direction when the position
  ! of atom1 due to above transformation is in the 1st or 2nd quadrant. While the rotation
  ! is to be performed in the counter clockwise direction if the atom1 is positioned
  ! in the 3rd or 4th quadrant.

  IF ( atom_list(iatom1,lm,is)%ryp >= 0.0_DP ) THEN
     delta_theta = -delta_theta
  END IF

  
  cos_dtheta = DCOS(delta_theta)
  sin_dtheta = DSIN(delta_theta)

  DO j = 1, natoms_to_place

     this_atom = atoms_to_place_list(j)

     tempx = atom_list(this_atom,lm,is)%rxp
     tempy = atom_list(this_atom,lm,is)%ryp
     
     atom_list(this_atom,lm,is)%rxp = tempx * cos_dtheta + tempy * sin_dtheta
     atom_list(this_atom,lm,is)%ryp = -tempx * sin_dtheta + tempy * cos_dtheta

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

  ! Calculate the energies after the move. First compute intramolecular and intermolecular
  ! nonbonded interactions so that the move can be immediately rejected if an overlap is detected.
  ! Since COM cutoff has been enabled. Compute the new COM of the molecule. Note that if the 
  ! move is rejected then automatically revert_cartesian_coordinates will reset the max_dcom

  CALL Get_COM(lm,is)
  CALL Compute_Max_COM_Distance(lm,is)

  CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is,E_intra_vdw_move,E_intra_qq_move,E_periodic_qq,intra_overlap)
  CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw_move,E_inter_qq_move,inter_overlap)
  E_inter_qq_move = E_inter_qq_move + E_periodic_qq


  IF (inter_overlap) THEN
     
     CALL Revert_Old_Cartesian_Coordinates(lm,is)
     CALL Revert_Old_Internal_Coordinates(lm,is)
     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(lm,is,ibox)
     ! Note that there is no need to reset the sin and cos terms for reciprocal space Ewald
     ! as these energies have not been yet computed.

     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,A9)') &
             i_mcstep, 'angle' , lm, is, ibox, accept, 'overlap'
     END IF
  ELSE

     delta_e = 0.0_DP

     CALL Compute_Molecule_Bond_Energy(lm,is,E_bond_move)
     CALL Compute_Molecule_Angle_Energy(lm,is,E_angle_move)
     CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihedral_move)
     CALL Compute_Molecule_Improper_Energy(lm,is,E_improper_move)

     IF ( int_charge_sum_style(ibox) == charge_ewald .and. has_charge(is)) THEN
        
        ALLOCATE(cos_mol_old(nvecs(ibox)),sin_mol_old(nvecs(ibox)))
        CALL Get_Position_Alive(lm,is,position)
        
        cos_mol_old(:) = cos_mol(1:nvecs(ibox),position)
        sin_mol_old(:) = sin_mol(1:nvecs(ibox),position)
        
        CALL Update_System_Ewald_Reciprocal_Energy(lm,is,ibox, &
             int_intra,E_reciprocal_move)
        delta_e = (E_reciprocal_move - energy(ibox)%ewald_reciprocal)
            
     END IF
     
     ! energy difference

     delta_e = (E_bond_move - E_bond) &
             + (E_angle_move - E_angle) &
             + (E_dihedral_move - E_dihedral) &
             + (E_improper_move - E_improper) &
             + (E_intra_vdw_move - E_intra_vdw) &
             + (E_intra_qq_move - E_intra_qq) &
             + (E_inter_vdw_move - E_inter_vdw) &
             + (E_inter_qq_move - E_inter_qq) + delta_e

     IF ( int_sim_type == sim_nvt_min ) THEN
        IF ( delta_e <= 0.0_DP ) THEN
           accept = .TRUE.
        ELSE
           accept = .FALSE.
        END IF
        
     ELSE

        ln_pacc = beta(ibox) * delta_e
        accept = accept_or_reject(ln_pacc)

     END IF
  
     IF ( accept ) THEN
        !     write(*,*) 'accepted', theta_new
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

        nsuccess(is,ibox)%angle = nsuccess(is,ibox)%angle + 1

        ! Update the virial as well
        ! Note that COM and max_dcom have already been updated

        CALL Get_Internal_Coordinates(lm,is)

        ! Fold the molecule

        CALL Fold_Molecule(lm,is,ibox)
        
        IF(ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
        IF(ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)
        
        IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
        
     ELSE 

        ! Reject the move and revert the old coordinates

        CALL Revert_Old_Cartesian_Coordinates(lm,is)
        CALL Revert_Old_Internal_Coordinates(lm,is)

        IF(l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(lm,is,ibox)
        
        IF (int_charge_sum_style(ibox) == charge_ewald .AND. has_charge(is)) THEN
           ! Also reset the old cos_sum and sin_sum for reciprocal space vectors
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
           cos_sum(:,ibox) = cos_sum_old(:,ibox)
           sin_sum(:,ibox) = sin_sum_old(:,ibox)
           
           cos_mol(1:nvecs(ibox),position) = cos_mol_old(:)
           sin_mol(1:nvecs(ibox),position) = sin_mol_old(:)
           !$OMP END PARALLEL WORKSHARE
           
           DEALLOCATE(cos_mol_old,sin_mol_old)

        END IF
        
     END IF

     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,F9.3)') &
             i_mcstep, 'angle' , lm, is, ibox, accept, ln_pacc
     END IF

  END IF
     
  IF(ALLOCATED(atoms_to_place_list)) DEALLOCATE(atoms_to_place_list)

END SUBROUTINE Angle_Distortion
