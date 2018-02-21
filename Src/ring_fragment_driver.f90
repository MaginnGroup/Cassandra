
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
