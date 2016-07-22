
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
!
! This file contains routines to sample ring fragment 
!
! It is based on the flip move outlined in Peristeras et al., Macromolecules, 2005, 38, 386-397
!
! Sampling of ring bond angles are accomplished by flipping a branched ring atom
! around the axis formed by neighboring ring atoms. Any exocyclic atoms attached
! to the atom being flipped will also undergo rotation about the axis. This 
! preserves the bond length constraint. Jacobian of the move cancels out for 
! forward and reverse moves. Please see the reference for further details.
!
! There are three routines:
!
! Ring_Fragment
!  
!     This is the driver routine that performs the two moves - 
!
! Flip_Move
!
!     Carries out moves to sample endocyclic degrees of freedom of the ring
!
! 08/07/13 : Created beta version
!*******************************************************************************

SUBROUTINE Ring_Fragment_Driver
  !*****************************************************************************
  !
  ! CALLED BY
  !
  !        main
  !
  ! CALLS
  !
  !        Flip_Move
  !        Atom_Displacement
  !
  !*****************************************************************************

  USE Global_Variables
  USE Random_Generators
  USE Energy_Routines
  
  IMPLICIT NONE

  INTEGER :: is, im, ibox, i, ia, nring_success, nexoring_success
  INTEGER :: nexoring_trials, nring_trials

  REAL(DP) :: rand_no 

  LOGICAL :: overlap

  ! Fragment routine is only called with 1 box, species, and molecule
  is = 1
  im = 1
  ibox = 1

  nring_success = 0
  nexoring_success = 0
  nring_trials = 0
  nexoring_trials = 0
  ! obtain energy of starting conformation
  
  OPEN(UNIT=frag_file_unit,file=frag_file(is))

  WRITE(frag_file_unit,*) n_mcsteps/nthermo_freq

  DO i_mcstep = 1, n_mcsteps
     rand_no = rranf()

     IF ( rand_no <= cut_ring) THEN
        
        nring_trials = nring_trials + 1
        CALL Flip_Move

        IF (accept) THEN
           nring_success = nring_success + 1
        END IF

     ELSE IF (rand_no <= cut_atom_displacement) THEN
        
        nexoring_trials = nexoring_trials + 1
        CALL Atom_Displacement

        IF (accept) THEN
           nexoring_success = nexoring_success + 1
        END IF

     END IF

     ! Store information with given frequency
     IF (MOD(i_mcstep,nthermo_freq) == 0) THEN
        WRITE(frag_file_unit,*) temperature(ibox), energy(ibox)%dihedral+ &
             energy(ibox)%intra_vdw + energy(ibox)%intra_q + &
             energy(ibox)%improper
!temperature(ibox), energy(ibox)%total

        DO ia = 1, natoms(is)
           WRITE(frag_file_unit,*) nonbond_list(ia,is)%element, &
                atom_list(ia,im,is)%rxp, atom_list(ia,im,is)%ryp, &
                atom_list(ia,im,is)%rzp
        END DO

     END IF

  END DO

  CLOSE(UNIT=frag_file_unit)

  WRITE(*,'(X,A,T40,I8)') 'Number of ring trials', nring_trials
  WRITE(*,'(X,A,T40,I8)') 'Number of successful ring trials', nring_success
  WRITE(*,'(X,A,T40,I8)') 'Number of exoring trials', nexoring_trials
  WRITE(*,'(X,A,T40,I8)') 'Number of successful exoring trials', nexoring_success

END SUBROUTINE Ring_Fragment_Driver
!*******************************************************************************
SUBROUTINE Flip_Move
  !*****************************************************************************
  !
  ! CALLED BY
  !
  !        nvt_mc_ring_fragment
  !
  ! CALLS
  !
  !        Save_Old_Cartesian_Coordinates
  !        Get_Aligner_Hanger
  !        Clean_Abort
  !        Compute_Molecule_Angle_Energy
  !        Compute_Molecule_Dihedral_Energy
  !        Compute_Molecule_Nonbond_Intra_Energy
  !        Update_System_Ewald_Reciprocal_Energy
  !        Revert_Old_Cartesian_Coordinates
  !*****************************************************************************

  USE Global_Variables
  USE Random_Generators
  USE Fragment_Growth, ONLY : Get_Aligner_Hanger
  USE IO_Utilities
  USE Energy_Routines

  IMPLICIT NONE

  INTEGER :: atom_num, this_atom, i, angle_id, atom1, atom2, is, im, ibox, j
  INTEGER :: natoms_to_move, atom_i
  INTEGER, DIMENSION(:), ALLOCATABLE :: atoms_to_move

  REAL(DP) :: vec1(3), vec2(3), aligner(3,3), hanger(3,3), tempy, tempz
  REAL(DP) :: domega, cosomega, sinomega, a(3), b(3), c(3)

  REAL(DP) :: dE
  REAL(DP) :: E_bond, E_angle, E_dihed, E_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq, E_periodic_qq, E_recip
  REAL(DP) :: ln_pacc

  LOGICAL :: accept_or_reject,intra_overlap

  im = 1
  is = 1
  ibox = 1

  ! save the orignal coordinates
  Call Save_Old_Cartesian_Coordinates(im,is)
  
  ! Choose a ring atom to flip
  atom_num = INT ( rranf() * nring_atoms(is)) + 1

  ! obtain id 
  this_atom = ring_atom_ids(atom_num,is)

  ! Figure out the bond angle that contains this_atom at the apex
  DO i = 1, angle_part_list(this_atom,is)%nangles

     IF (angle_part_list(this_atom,is)%position(i) == 2 ) THEN
        ! this_atom is the apex
        angle_id = angle_part_list(this_atom,is)%which_angle(i)
        
        ! check if other atoms are ring atoms
        atom1 = angle_list(angle_id,is)%atom1
        atom2 = angle_list(angle_id,is)%atom3

        IF (nonbond_list(atom1,is)%ring_atom .AND. nonbond_list(atom2,is)%ring_atom) EXIT

     END IF

  END DO

 

  ! Rotation will be achieved around the axis formed by atom1-atom2

  ! Align plane of atom1, atom2 and this_atom such that atom1 is at the origin, atom2 is 
  ! along the x-axis and this_atom is in xy plane

  vec1(1) = atom_list(atom2,im,is)%rxp - atom_list(atom1,im,is)%rxp
  vec1(2) = atom_list(atom2,im,is)%ryp - atom_list(atom1,im,is)%ryp
  vec1(3) = atom_list(atom2,im,is)%rzp - atom_list(atom1,im,is)%rzp
  
  vec2(1) = atom_list(this_atom,im,is)%rxp - atom_list(atom1,im,is)%rxp
  vec2(2) = atom_list(this_atom,im,is)%ryp - atom_list(atom1,im,is)%ryp
  vec2(3) = atom_list(this_atom,im,is)%rzp - atom_list(atom1,im,is)%rzp


  CALL Get_Aligner_Hanger(vec1,vec2,aligner,hanger)

  ! determine number of atoms to which the transformation will be applied. These atoms
  ! are 'this_atom' and exocyclic atoms attached to 'this_atom'
  
  ALLOCATE(atoms_to_move(natoms(is)), Stat = AllocateStatus)
  atoms_to_move(:) = 0

  IF ( AllocateStatus /= 0 ) THEN
     err_msg = ''
     err_msg(1) = ' Memory could be not assigned to atoms_to_move'
     CALL Clean_Abort(err_msg,'Flip_Move')
  END IF

  natoms_to_move = 1   ! corresponding to this_atom
  atoms_to_move(natoms_to_move) = this_atom

  ! Now figure out all the exocyclic atoms attached to 'this_atom'

  DO i = 1, bondpart_list(this_atom,is)%nbonds

     atom_i = bondpart_list(this_atom,is)%atom(i)
     
     IF (.NOT. nonbond_list(atom_i,is)%ring_atom) THEN

        natoms_to_move = natoms_to_move + 1
        atoms_to_move(natoms_to_move) = atom_i

     END IF

  END DO


  
  ! Apply three transformations
  ! first --- aligner 
  ! second -- rotation about the axis
  ! third -- hanger to transform coordinates to lab frame of reference

  ! obtain rotation angle before hand

  ! Rotate the atoms around axis of atom1-atom2 by dw

  domega = (2.0_DP * rranf() - 1.0_DP) * omega_max
  
  cosomega = DCOS(domega)
  sinomega = DSIN(domega)

  DO i = 1, natoms_to_move

     atom_i = atoms_to_move(i)


     a(1) = atom_list(atom_i,im,is)%rxp - atom_list(atom1,im,is)%rxp
     a(2) = atom_list(atom_i,im,is)%ryp - atom_list(atom1,im,is)%ryp
     a(3) = atom_list(atom_i,im,is)%rzp - atom_list(atom1,im,is)%rzp
     
     ! first transformation

     b = MATMUL(aligner,a)
    
     ! second transformation
     tempy = b(2)
     tempz = b(3)
     
     b(2) = cosomega * tempy + sinomega * tempz
     b(3) = -sinomega * tempy + cosomega * tempz
     
     ! third transformation

     c = MATMUL(hanger,b)

     atom_list(atom_i,im,is)%rxp = c(1) + atom_list(atom1,im,is)%rxp
     atom_list(atom_i,im,is)%ryp = c(2) + atom_list(atom1,im,is)%ryp
     atom_list(atom_i,im,is)%rzp = c(3) + atom_list(atom1,im,is)%rzp
        
  END DO

 
  DEALLOCATE(atoms_to_move)

  ! Compute energy changes 

  dE = 0.0_DP
  CALL Compute_Molecule_Bond_Energy(im,is,E_bond)
  CALL Compute_Molecule_Angle_Energy(im,is,E_angle)
  CALL Compute_Molecule_Dihedral_Energy(im,is,E_dihed)
  CALL Compute_Molecule_Improper_Energy(im,is,E_improper)

  ! Note that we will not compute the reciprocal space Ewald in the ring biasing
  ! case, nor will we use the E_periodic_qq
  CALL Compute_Molecule_Nonbond_Intra_Energy(im,is,E_intra_vdw,E_intra_qq, &
       E_periodic_qq,intra_overlap)

  IF (int_charge_sum_style(ibox) == charge_ewald) THEN
     ! Compute Ewald reciprocal energy difference
     CALL Update_System_Ewald_Reciprocal_Energy(im,is,ibox,int_intra,E_recip)
     dE = (E_recip - energy(ibox)%reciprocal)
  END IF

  ! Change in energy due to the move

  dE = dE + E_angle - energy(ibox)%angle
  dE = dE + E_dihed - energy(ibox)%dihedral
  dE = dE + E_improper - energy(ibox)%improper
  dE = dE + E_intra_vdw - energy(ibox)%intra_vdw
  dE = dE + E_intra_qq - energy(ibox)%intra_q

  ln_pacc = beta(ibox) * dE

  accept = accept_or_reject(ln_pacc)

  IF (accept) THEN

     ! update energies
     energy(ibox)%total = energy(ibox)%total + dE
     energy(ibox)%intra = energy(ibox)%intra + dE
     energy(ibox)%angle = E_angle
     energy(ibox)%dihedral = E_dihed
     energy(ibox)%improper = E_improper
     energy(ibox)%intra_vdw = E_intra_vdw
     energy(ibox)%intra_q = E_intra_qq
     
     IF (int_charge_sum_style(ibox) == charge_ewald .AND. &
         has_charge(is)) THEN
        energy(ibox)%reciprocal = E_recip
     END IF
  ELSE

     CALL Revert_Old_Cartesian_Coordinates(im,is)
   
  END IF

  IF (verbose_log) THEN
    WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8)') i_mcstep, 'flip_move' , im, is, ibox, accept
  END IF


END SUBROUTINE Flip_Move

!************************************************************************************
