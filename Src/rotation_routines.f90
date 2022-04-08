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

   TYPE(Molecule_Class), POINTER :: this_molecule
   TYPE(Atom_Class), POINTER :: these_atoms(:)

   IF (widom_active) THEN
           this_molecule => widom_molecule
           these_atoms => widom_atoms
   ELSE
           this_molecule => molecule_list(alive,is)
           these_atoms => atom_list(:,alive,is)
   END IF
   
   ! Pick random eulerians
   
   theta = ACOS(1.0_DP-2.0_DP*rranf())
   phi  = (1.0_DP - 2.0_DP * rranf()) * PI
   psi   = (1.0_DP - 2.0_DP * rranf()) * PI
   
   ! shift the origin to the COM of the molecule, do this for atoms that are currently 
   ! present in the simulation box so that the routine can be used for partially grown
   ! molecules
   
   DO ia = 1,natoms(is)

      IF (these_atoms(ia)%exist) THEN
         
         these_atoms(ia)%rxp = these_atoms(ia)%rxp - this_molecule%xcom
         these_atoms(ia)%ryp = these_atoms(ia)%ryp - this_molecule%ycom
         these_atoms(ia)%rzp = these_atoms(ia)%rzp - this_molecule%zcom
         
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

      IF (these_atoms(ia)%exist) THEN
         
         rxpnew = rot11*these_atoms(ia)%rxp + rot12*these_atoms(ia)%ryp + &
              rot13*these_atoms(ia)%rzp
         rypnew = rot21*these_atoms(ia)%rxp + rot22*these_atoms(ia)%ryp + &
              rot23*these_atoms(ia)%rzp
         rzpnew = rot31*these_atoms(ia)%rxp + rot32*these_atoms(ia)%ryp + &
              rot33*these_atoms(ia)%rzp
         
         ! Shift the origin back to (0,0,0)
         
         these_atoms(ia)%rxp = rxpnew + this_molecule%xcom
         these_atoms(ia)%ryp = rypnew + this_molecule%ycom
         these_atoms(ia)%rzp = rzpnew + this_molecule%zcom
         
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

   TYPE(Molecule_Class), POINTER :: this_molecule
   TYPE(Atom_Class), POINTER :: these_atoms(:)

   IF (widom_active) THEN
           this_molecule => widom_molecule
           these_atoms => widom_atoms
   ELSE
           this_molecule => molecule_list(alive,is)
           these_atoms => atom_list(:,alive,is)
   END IF
   
   
   ! Pick random eulerians
   
   theta = ACOS(1.0_DP-2.0_DP*rranf())
   phi  = (1.0_DP - 2.0_DP * rranf()) * PI
   psi   = (1.0_DP - 2.0_DP * rranf()) * PI
   
   ! shift the origin to the COM of the molecule, do this for atoms that are currently 
   ! present in the simulation box so that the routine can be used for partially grown
   ! molecules

   atom_orig = frag_list(frag_start,is)%atoms(1)

   x_orig = these_atoms(atom_orig)%rxp
   y_orig = these_atoms(atom_orig)%ryp
   z_orig = these_atoms(atom_orig)%rzp

   istart = 2

   DO i = istart,frag_list(frag_start,is)%natoms
      ia = frag_list(frag_start,is)%atoms(i)
 
      IF (these_atoms(ia)%exist) THEN
         
         these_atoms(ia)%rxp = these_atoms(ia)%rxp - x_orig
         these_atoms(ia)%ryp = these_atoms(ia)%ryp - y_orig
         these_atoms(ia)%rzp = these_atoms(ia)%rzp - z_orig
         
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

      IF (these_atoms(ia)%exist) THEN
         
         rxpnew = rot11*these_atoms(ia)%rxp + rot12*these_atoms(ia)%ryp + &
              rot13*these_atoms(ia)%rzp
         rypnew = rot21*these_atoms(ia)%rxp + rot22*these_atoms(ia)%ryp + &
              rot23*these_atoms(ia)%rzp
         rzpnew = rot31*these_atoms(ia)%rxp + rot32*these_atoms(ia)%ryp + &
              rot33*these_atoms(ia)%rzp
         
         these_atoms(ia)%rxp = rxpnew
         these_atoms(ia)%ryp = rypnew
         these_atoms(ia)%rzp = rzpnew

      END IF

   END DO

   ! Shift the origin back to (0,0,0)
   DO i=istart,frag_list(frag_start,is)%natoms
      ia = frag_list(frag_start,is)%atoms(i)
      these_atoms(ia)%rxp = these_atoms(ia)%rxp + x_orig
      these_atoms(ia)%ryp = these_atoms(ia)%ryp + y_orig
      these_atoms(ia)%rzp = these_atoms(ia)%rzp + z_orig
   END DO
 END SUBROUTINE Rotate_XYZ_Axes

END MODULE Rotation_Routines
