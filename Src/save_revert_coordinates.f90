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

!*******************************************************************************
! This file contains a collection of routines that perform the task of saving
! old coordinates of molecule that is being perturbed due to any of the following
! moves
!
! translation
! rotation
! rigid dihedral angle rotation
! angle distortion 
! bond stretching
! 
! It also contains routines that are used to reset the old cartesian and internal
! coordinates of the molecule if the attempted move is rejected.
!
! Note that in all the routines input 'im' is the linked number for the input
! molecule i.e. im = locate(molecule #,species #,box #)
!
! Written by Jindal Shah
! 
! Revision History
!
! 12/10/13 : Beta Release
!********************************************************************************

!********************************************************************************
SUBROUTINE Save_Old_Cartesian_Coordinates(im,is)
!********************************************************************************

!********************************************************************************
! This routine saves the old coordinates of the molecule im of species is.
!
! The routine gets called by
!
! translate.f90
! rotate.f90
! rigid_dihedral_change.f90
! angle_distortion_change.f90
! bond_stretching_move.f90
!
!********************************************************************************

  USE Type_Definitions
  USE Global_Variables

  IMPLICIT NONE

  INTEGER ::  is, im


! Save the parent coordinates

  atom_list(:,im,is)%rxp_old = atom_list(:,im,is)%rxp
  atom_list(:,im,is)%ryp_old = atom_list(:,im,is)%ryp
  atom_list(:,im,is)%rzp_old = atom_list(:,im,is)%rzp


! Save the COM and Eulerian angles
  
  molecule_list(im,is)%xcom_old = molecule_list(im,is)%xcom
  molecule_list(im,is)%ycom_old = molecule_list(im,is)%ycom
  molecule_list(im,is)%zcom_old = molecule_list(im,is)%zcom

  molecule_list(im,is)%euler1_old = molecule_list(im,is)%euler1
  molecule_list(im,is)%euler2_old = molecule_list(im,is)%euler2
  molecule_list(im,is)%euler3_old = molecule_list(im,is)%euler3

  molecule_list(im,is)%max_dcom_old = molecule_list(im,is)%max_dcom

END SUBROUTINE Save_Old_Cartesian_Coordinates

!********************************************************************************
SUBROUTINE Save_Old_Internal_Coordinates(im,is)
!********************************************************************************

!********************************************************************************
! Save all the internal coordinates. Note that this step is not
! really required for translation and rigid body rotation of the molecule.
!
! The routine gets called by 
!
! rigid_dihedral_change
! angle_distortion_change
! bond_stretching_move

! Revision history
!
! 12/10/13: Beta release
!
!********************************************************************************
  USE Type_Definitions
  USE Global_Variables

  IMPLICIT NONE

  INTEGER :: im, is, i, max_index

  internal_coord_list_old(:)%bond_length_angstrom = 0.0_DP
  internal_coord_list_old(:)%bond_angle_degrees = 0.0_DP
  internal_coord_list_old(:)%bond_angle_radians = 0.0_DP
  internal_coord_list_old(:)%dihedral_angle_degrees = 0.0_DP
  internal_coord_list_old(:)%dihedral_angle_radians = 0.0_DP
  internal_coord_list_old(:)%improper_angle_degrees = 0.0_DP
  internal_coord_list_old(:)%improper_angle_radians = 0.0_DP

  ! store the bonds

  DO i = 1, nbonds(is)

     internal_coord_list_old(i)%bond_length_angstrom = &
          internal_coord_list(i,im,is)%bond_length_angstrom

  END DO

  ! store the angles
  
  DO i = 1, nangles(is)
     
     internal_coord_list_old(i)%bond_angle_degrees = &
          internal_coord_list(i,im,is)%bond_angle_degrees
     
     internal_coord_list_old(i)%bond_angle_radians = &
          internal_coord_list(i,im,is)%bond_angle_radians
     
  END DO

  ! store the dihedrals

  DO i = 1, ndihedrals(is)
     
     internal_coord_list_old(i)%dihedral_angle_degrees = &
          internal_coord_list(i,im,is)%dihedral_angle_degrees

     internal_coord_list_old(i)%dihedral_angle_radians = &
          internal_coord_list(i,im,is)%dihedral_angle_radians

  END DO

  ! store the impropers

  DO i = 1, nimpropers(is)

     internal_coord_list_old(i)%improper_angle_degrees = &
          internal_coord_list(i,im,is)%improper_angle_degrees

     internal_coord_list_old(i)%improper_angle_radians = &
          internal_coord_list(i,im,is)%improper_angle_radians

  END DO

END SUBROUTINE Save_Old_Internal_Coordinates

!*******************************************************************************
SUBROUTINE Revert_Old_Cartesian_Coordinates(im,is)
!*******************************************************************************

!*******************************************************************************
! This routine performs the function of resetting the old coordinates of the
! molecule im of species is. The routine gets called when a MC move is rejected 
! so we need to revert to the old configuration of im.
! 
! All the moves that perturb the coordinates of molecules call this routine
! when the attempted move gets rejected.
!
! Called by:
!
! translation
! rotation
! rigid_dihedral_change
! angle_distortion_change
! bond_stretching_move
!
!*******************************************************************************

  USE Type_Definitions
  USE Global_Variables

  IMPLICIT NONE

  INTEGER ::  is, im
  
  
! Revert to the old x,y and z parent coordinates of the atoms

  atom_list(:,im,is)%rxp = atom_list(:,im,is)%rxp_old
  atom_list(:,im,is)%ryp = atom_list(:,im,is)%ryp_old
  atom_list(:,im,is)%rzp = atom_list(:,im,is)%rzp_old


! Revert to the COM and Eulerian angles for the molecule

  molecule_list(im,is)%xcom = molecule_list(im,is)%xcom_old
  molecule_list(im,is)%ycom = molecule_list(im,is)%ycom_old
  molecule_list(im,is)%zcom = molecule_list(im,is)%zcom_old

  molecule_list(im,is)%euler1 = molecule_list(im,is)%euler1_old
  molecule_list(im,is)%euler2 = molecule_list(im,is)%euler2_old
  molecule_list(im,is)%euler3 = molecule_list(im,is)%euler3_old

  molecule_list(im,is)%max_dcom = molecule_list(im,is)%max_dcom_old
        
END SUBROUTINE Revert_Old_Cartesian_Coordinates

!**********************************************************************************
SUBROUTINE Revert_Old_Internal_Coordinates(im,is)
!**********************************************************************************

!**********************************************************************************
! The subroutine resets the internal coordinates of im of species is after an
! internal coordinate change move is rejected.
! 
! Called by
!
! rigid_dihedral_change
! angle_distortion_change
! bond_stretching_move
!
!**********************************************************************************

  USE Type_Definitions
  USE Global_Variables

  IMPLICIT NONE

  INTEGER :: im, is, i

  ! bonds

  DO i = 1, nbonds(is)

     internal_coord_list(i,im,is)%bond_length_angstrom = &
          internal_coord_list_old(i)%bond_length_angstrom

  END DO

 ! angles

  DO i = 1, nangles(is)

     internal_coord_list(i,im,is)%bond_angle_degrees = &
          internal_coord_list_old(i)%bond_angle_degrees

     internal_coord_list(i,im,is)%bond_angle_radians = &
          internal_coord_list_old(i)%bond_angle_radians 
        
  END DO
        
  
 ! dihedrals

  DO i = 1, ndihedrals(is)

     internal_coord_list(i,im,is)%dihedral_angle_degrees = &
          internal_coord_list_old(i)%dihedral_angle_degrees 
    
     internal_coord_list(i,im,is)%dihedral_angle_radians = &
          internal_coord_list_old(i)%dihedral_angle_radians 


  END DO

  ! impropers

  DO i = 1, nimpropers(is)

     internal_coord_list(i,im,is)%improper_angle_degrees = &
          internal_coord_list_old(i)%improper_angle_degrees

     internal_coord_list(i,im,is)%improper_angle_radians = &
          internal_coord_list_old(i)%improper_angle_radians 

  END DO

 

END SUBROUTINE Revert_Old_Internal_Coordinates
