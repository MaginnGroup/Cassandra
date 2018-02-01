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

MODULE Rotation_Routines

  !******************************************************************************
  ! 
  ! This module contains rotations routines that are used to perform rotation
  ! moves. It consists of rotation about a randomly chosen axis and eulerian angle
  ! change
  !
  !
  !
  !*******************************************************************************

  USE Global_Variables
  USE Random_Generators

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Rotate_Molecule_Eulerian(alive,is)

    
    !*********************************************************************************
    ! takes in the identity of the molecule and rotates the molecule with random
    ! eulerian angles
    !
    ! Gets called by
    !
    !    Insertion.f90
    !
    !
    ! 12/10/13  : Beta Release
    !*********************************************************************************

    
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: alive, is
   
   REAL(DP) :: theta, phi, psi, rot11, rot12, rot13, rot21, rot22, rot23
   REAL(DP) :: rot31, rot32, rot33, rxpnew, rypnew, rzpnew
   
   INTEGER :: i, ia
   
   ! Pick random eulerians
   
   theta = ACOS(1.0_DP-2.0_DP*rranf())
   phi  = (1.0_DP - 2.0_DP * rranf()) * PI
   psi   = (1.0_DP - 2.0_DP * rranf()) * PI
   
   ! shift the origin to the COM of the molecule, do this for atoms that are currently 
   ! present in the simulation box so that the routine can be used for partially grown
   ! molecules
   
   DO ia = 1,natoms(is)

      IF (atom_list(ia,alive,is)%exist) THEN
         
         atom_list(ia,alive,is)%rxp = atom_list(ia,alive,is)%rxp - molecule_list(alive,is)%xcom
         atom_list(ia,alive,is)%ryp = atom_list(ia,alive,is)%ryp - molecule_list(alive,is)%ycom
         atom_list(ia,alive,is)%rzp = atom_list(ia,alive,is)%rzp - molecule_list(alive,is)%zcom
         
      END IF

   END DO

   ! Construct the rotation matrix that needs to be applied to each of the vectors
   ! This is the A matrix in Goldstein notation
   
   rot11 = DCOS(psi) * DCOS(phi) - DCOS(theta) * DSIN(phi) * DSIN(psi)
   rot12 = DCOS(psi) * DSIN(phi) + DCOS(theta) * DCOS(phi) * DSIN(psi)
   rot13 = DSIN(psi) * DSIN(theta)
   
   rot21 = -DSIN(psi) * DCOS(phi) - DCOS(theta) * DSIN(phi) * DCOS(psi)
   rot22 = -DSIN(psi) * DSIN(phi) + DCOS(theta) * DCOS(phi) * DCOS(psi)
   rot23 = DCOS(psi) * DSIN(theta)

   rot31 = DSIN(theta) * DSIN(phi)
   rot32 = -DSIN(theta) * DCOS(phi)
   rot33 = DCOS(theta)

   ! Now rotate the relative vectors to obtain the new positions

   DO ia = 1, natoms(is)

      IF (atom_list(ia,alive,is)%exist) THEN
         
         rxpnew = rot11*atom_list(ia,alive,is)%rxp + rot12*atom_list(ia,alive,is)%ryp + &
              rot13*atom_list(ia,alive,is)%rzp
         rypnew = rot21*atom_list(ia,alive,is)%rxp + rot22*atom_list(ia,alive,is)%ryp + &
              rot23*atom_list(ia,alive,is)%rzp
         rzpnew = rot31*atom_list(ia,alive,is)%rxp + rot32*atom_list(ia,alive,is)%ryp + &
              rot33*atom_list(ia,alive,is)%rzp
         
         ! Shift the origin back to (0,0,0)
         
         atom_list(ia,alive,is)%rxp = rxpnew + molecule_list(alive,is)%xcom
         atom_list(ia,alive,is)%ryp = rypnew + molecule_list(alive,is)%ycom
         atom_list(ia,alive,is)%rzp = rzpnew + molecule_list(alive,is)%zcom
         
      END IF
      
   END DO

 END SUBROUTINE Rotate_Molecule_Eulerian


 SUBROUTINE Rotate_XYZ_Axes(alive,is,frag_start,lx,ly,lz,mtype) 
    
    !*********************************************************************************
    ! takes in the identity of the molecule and rotates the molecule with random
    ! eulerian angles
    !
    ! Gets called by
    !
    !    Insertion.f90
    !
    !
    ! 12/10/13 (EM): Beta Release
    !*********************************************************************************

    
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: alive, is,frag_start,  mtype
   INTEGER :: atom_orig
   LOGICAL :: lx, ly,lz
   REAL(DP) :: theta, phi, psi, rot11, rot12, rot13, rot21, rot22, rot23
   REAL(DP) :: rot31, rot32, rot33, rxpnew, rypnew, rzpnew
   
   REAL(DP) :: x_orig,y_orig,z_orig

   INTEGER :: i, ia, istart
   
   ! Pick random eulerians
   
   theta = ACOS(1.0_DP-2.0_DP*rranf())
   phi  = (1.0_DP - 2.0_DP * rranf()) * PI
   psi   = (1.0_DP - 2.0_DP * rranf()) * PI
   
   ! shift the origin to the COM of the molecule, do this for atoms that are currently 
   ! present in the simulation box so that the routine can be used for partially grown
   ! molecules

   atom_orig = frag_list(frag_start,is)%atoms(1)

   x_orig = atom_list(atom_orig,alive,is)%rxp
   y_orig = atom_list(atom_orig,alive,is)%ryp
   z_orig = atom_list(atom_orig,alive,is)%rzp

   istart = 2

   DO i = istart,frag_list(frag_start,is)%natoms
      ia = frag_list(frag_start,is)%atoms(i)
 
      IF (atom_list(ia,alive,is)%exist) THEN
         
         atom_list(ia,alive,is)%rxp = atom_list(ia,alive,is)%rxp - x_orig
         atom_list(ia,alive,is)%ryp = atom_list(ia,alive,is)%ryp - y_orig
         atom_list(ia,alive,is)%rzp = atom_list(ia,alive,is)%rzp - z_orig
         
      END IF

   END DO

   ! Construct the rotation matrix that needs to be applied to each of the vectors
   ! This is the A matrix in Goldstein notation
   
   rot11 = DCOS(psi) * DCOS(phi) - DCOS(theta) * DSIN(phi) * DSIN(psi)
   rot12 = DCOS(psi) * DSIN(phi) + DCOS(theta) * DCOS(phi) * DSIN(psi)
   rot13 = DSIN(psi) * DSIN(theta)
   
   rot21 = -DSIN(psi) * DCOS(phi) - DCOS(theta) * DSIN(phi) * DCOS(psi)
   rot22 = -DSIN(psi) * DSIN(phi) + DCOS(theta) * DCOS(phi) * DCOS(psi)
   rot23 = DCOS(psi) * DSIN(theta)

   rot31 = DSIN(theta) * DSIN(phi)
   rot32 = -DSIN(theta) * DCOS(phi)
   rot33 = DCOS(theta)

   ! Now rotate the relative vectors to obtain the new positions

   DO ia = 2, natoms(is)

      IF (atom_list(ia,alive,is)%exist) THEN
         
         rxpnew = rot11*atom_list(ia,alive,is)%rxp + rot12*atom_list(ia,alive,is)%ryp + &
              rot13*atom_list(ia,alive,is)%rzp
         rypnew = rot21*atom_list(ia,alive,is)%rxp + rot22*atom_list(ia,alive,is)%ryp + &
              rot23*atom_list(ia,alive,is)%rzp
         rzpnew = rot31*atom_list(ia,alive,is)%rxp + rot32*atom_list(ia,alive,is)%ryp + &
              rot33*atom_list(ia,alive,is)%rzp
         
         atom_list(ia,alive,is)%rxp = rxpnew
         atom_list(ia,alive,is)%ryp = rypnew
         atom_list(ia,alive,is)%rzp = rzpnew

      END IF

   END DO

   ! Shift the origin back to (0,0,0)
   DO i=istart,frag_list(frag_start,is)%natoms
      ia = frag_list(frag_start,is)%atoms(i)
      atom_list(ia,alive,is)%rxp = atom_list(ia,alive,is)%rxp + x_orig
      atom_list(ia,alive,is)%ryp = atom_list(ia,alive,is)%ryp + y_orig
      atom_list(ia,alive,is)%rzp = atom_list(ia,alive,is)%rzp + z_orig
   END DO
 END SUBROUTINE Rotate_XYZ_Axes

END MODULE Rotation_Routines
