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

SUBROUTINE Get_Internal_Coords

  !****************************************************************************
  ! This module contains collection of routines that can be used to calculate 
  ! internal coordinates of molecules.
  !
  ! Called by
  !
  !   main
  !   read_write_checkpoint
  !
  ! Revision history
  !
  !   12/10/13 : Beta version
  !
  !********************************************************************************

  USE Type_Definitions
  USE Global_Variables
  
  IMPLICIT NONE

  INTEGER :: ispecies, imolecules, ibonds, iangles, idihedrals, iimpropers, alive, ibox
  REAL(DP) :: bond_length, theta, phi
  
  DO ibox = 1, nbr_boxes
     DO ispecies = 1, nspecies
        DO imolecules = 1,nmols(ispecies,ibox)

           ! obtain the linked index associated with this molecule

           alive = locate(imolecules,ispecies,ibox)

           IF ( .NOT. molecule_list(alive,ispecies)%live) CYCLE
           
           DO ibonds = 1, nbonds(ispecies)
              
              CALL Get_Bond_Length(ibonds,alive,ispecies, bond_length)

              ! Assign the bond length to the internal coordinate

              internal_coord_list(ibonds,alive,ispecies)%bond_length_angstrom = &
                   bond_length
              
           END DO
           
           DO iangles = 1, nangles(ispecies)

              CALL Get_Bond_Angle(iangles,alive,ispecies,theta)

              ! Assign this angle to internal_coord_list

              internal_coord_list(iangles,alive,ispecies)%bond_angle_radians = &
                   theta

              ! Convert the angle into degrees

              internal_coord_list(iangles,alive,ispecies)%bond_angle_degrees = &
                   theta * 180.0 / PI
      
           END DO

           DO idihedrals = 1, ndihedrals(ispecies)
              
              CALL Get_Dihedral_Angle(idihedrals,alive,ispecies,phi)

              ! Assign this value to the internal_coord_list

              internal_coord_list(idihedrals,alive,ispecies)%dihedral_angle_radians &
                   = phi
              
              ! Convert the angle into degrees
              
              internal_coord_list(idihedrals,alive,ispecies)%dihedral_angle_degrees &
                   = phi * 180.0 / PI
              
           END DO
           
           DO iimpropers = 1, nimpropers(ispecies)
              
              CAll Get_Improper_Angle(iimpropers,alive,ispecies,phi)

              ! Assign this value to the internal_coord_list

              internal_coord_list(iimpropers,alive,ispecies)%improper_angle_radians &
                   = phi

              ! Convert the angle into degrees

              internal_coord_list(iimpropers,alive,ispecies)%improper_angle_degrees &
                   = phi * 180.0 / PI
              
           END DO
           
        END DO ! do loop on imolecules
     END DO ! do loop on ispecies
  END DO ! do loop on ibox
  
 
END SUBROUTINE Get_Internal_Coords

!********************************************************************************
SUBROUTINE Get_Bond_Length(this_bond,this_molecule,is,r21)
!********************************************************************************

!********************************************************************************
! This routine computes the bond length for a given bond and assign it to the
! internal_coord_list
!
! CALLED BY
! 
!        energy_routines
!               Compute_Bond_Energy
!               Compute_Molecule_Bond_Energy
!        nvt_mc_fragment_driver
!        get_com  
!********************************************************************************


   USE Type_Definitions
   USE Global_Variables

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: this_bond, this_molecule, is
   INTEGER             :: atom1, atom2, this_box
   REAL(DP)            :: rx21, ry21, rz21, r21sq, r21

  
! Obtain the atoms involved in this_bond

   atom1 = bond_list(this_bond,is)%atom1
   atom2 = bond_list(this_bond,is)%atom2

! Calculate vectors, bond lenth

   rx21 = atom_list(atom1,this_molecule,is)%rxp - atom_list(atom2,this_molecule,is)%rxp
   ry21 = atom_list(atom1,this_molecule,is)%ryp - atom_list(atom2,this_molecule,is)%ryp
   rz21 = atom_list(atom1,this_molecule,is)%rzp - atom_list(atom2,this_molecule,is)%rzp

!   this_box = molecule_list(this_molecule,is)%which_box
!   IF (l_cubic(this_box) == .FALSE.) THEN
!   CALL Minimum_Image_Separation(this_box,rx21,ry21,rz21,rx21,ry21,rz21)
!   END IF

   r21sq = rx21 * rx21 + ry21 * ry21 + rz21 * rz21
   r21   = DSQRT(r21sq)
   
! store the bond length in the interal_coords_list

 END SUBROUTINE Get_Bond_Length


!********************************************************************************
 SUBROUTINE Get_Bond_Angle(this_angle,this_molecule,is,theta)
!********************************************************************************

!********************************************************************************
! This routine calculates the value of a given angle and assigns it to internal
! coord_list.
!
! CALLED BY
! 
!        energy_routines
!             Compute_Angle_Energy
!             Compute_Molecule_Angle_Energy
!
!        get_com

!********************************************************************************
   
   USE Type_Definitions
   USE Global_Variables
 
   IMPLICIT NONE

   INTEGER, INTENT(IN) :: this_angle, this_molecule, is

   INTEGER             :: atom1, atom2, atom3, this_box
   REAL(DP)            :: rx21, ry21, rz21, rx32, ry32, rz32
   REAL(DP)            :: C22, C23, C33, inv_rt_C3322, costheta, theta


! Figure out three atoms in the angle

   atom1 = angle_list(this_angle,is)%atom1
   atom2 = angle_list(this_angle,is)%atom2
   atom3 = angle_list(this_angle,is)%atom3

! Calculate vectors, dot products and bond angle
   
!Vector r21 points from atom 1 to atom 2. Below, the components of this vector are calculated.

   rx21 = atom_list(atom2,this_molecule,is)%rxp - atom_list(atom1,this_molecule,is)%rxp
   ry21 = atom_list(atom2,this_molecule,is)%ryp - atom_list(atom1,this_molecule,is)%ryp
   rz21 = atom_list(atom2,this_molecule,is)%rzp - atom_list(atom1,this_molecule,is)%rzp
                                                                                   
! Vector r32 points from atom 2 to atom 3. Below the components are calculated.

   rx32 = atom_list(atom3,this_molecule,is)%rxp - atom_list(atom2,this_molecule,is)%rxp
   ry32 = atom_list(atom3,this_molecule,is)%ryp - atom_list(atom2,this_molecule,is)%ryp
   rz32 = atom_list(atom3,this_molecule,is)%rzp - atom_list(atom2,this_molecule,is)%rzp
   
!   this_box = molecule_list(this_molecule,is)%which_box
!   IF (l_cubic(this_box) == .FALSE.) THEN
!   CALL Minimum_Image_Separation(this_box,rx21,ry21,rz21,rx21,ry21,rz21)
!   CALL Minimum_Image_Separation(this_box,rx32,ry32,rz32,rx32,ry32,rz32)
!   END IF

! The C values are just various dot products between vectors
!r21 and r32 and the square root of C33*C22


   C23 = rx21 * rx32 + ry21 * ry32 + rz21 * rz32
   C22 = rx21 * rx21 + ry21 * ry21 + rz21 * rz21
   C33 = rx32 * rx32 + ry32 * ry32 + rz32 * rz32

   inv_rt_C3322 = 1.0_DP / DSQRT(C33*C22)
                                                                                
! Calculate the angle cosine and guard against numerical overflow problems

   costheta = C23 * inv_rt_C3322
   costheta = MIN(costheta, 0.999999_DP)
   costheta = MAX(costheta, -0.999999_DP)

! Calculate the angle. Note that the actual bond angle is pi - theta

   theta = PI - DACOS(costheta)

 END SUBROUTINE Get_Bond_Angle


!********************************************************************************
 SUBROUTINE Get_Dihedral_Angle(this_dihedral,this_molecule,is,phi)
!********************************************************************************

!********************************************************************************
! This subroutine calculates the value of the dihedral angle passed to it.
!
! CALLED BY
!
!        energy_routines
!
!              Compute_Dihedral_Energy
!              Compute_Molecule_Dihedral_Energy
!
!        get_com
!********************************************************************************

   USE Type_Definitions
   USE Global_Variables

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: this_dihedral, this_molecule, is

   INTEGER             :: atom1, atom2, atom3, atom4, this_box

   REAL(DP)            :: rx12, ry12, rz12, rx32, ry32, rz32
   REAL(DP)            :: rx34, ry34, rz34, mx, my, mz, nx, ny, nz
   REAL(DP)            :: msq, nsq, mdn, abs_m, abs_n
   REAL(DP)            :: cosphi, phi, r12dn
   
  
! Get the atoms involved in the dihedral of interest

   atom1 = dihedral_list(this_dihedral,is)%atom1
   atom2 = dihedral_list(this_dihedral,is)%atom2
   atom3 = dihedral_list(this_dihedral,is)%atom3
   atom4 = dihedral_list(this_dihedral,is)%atom4

! Vector r12 points from atom 2 to atom 1.  Below, the components of this vector are calculated.

   rx12 = atom_list(atom1,this_molecule,is)%rxp - atom_list(atom2,this_molecule,is)%rxp
   ry12 = atom_list(atom1,this_molecule,is)%ryp - atom_list(atom2,this_molecule,is)%ryp
   rz12 = atom_list(atom1,this_molecule,is)%rzp - atom_list(atom2,this_molecule,is)%rzp
   
! Vector r32 points from atom 2 to atom 3.  Below, the components of this vector are calculated.

   rx32 = atom_list(atom3,this_molecule,is)%rxp - atom_list(atom2,this_molecule,is)%rxp
   ry32 = atom_list(atom3,this_molecule,is)%ryp - atom_list(atom2,this_molecule,is)%ryp
   rz32 = atom_list(atom3,this_molecule,is)%rzp - atom_list(atom2,this_molecule,is)%rzp

! Vector r34 points from atom 4 to atom 3. Below the components of this vector are calculated.

   rx34 = atom_list(atom3,this_molecule,is)%rxp - atom_list(atom4,this_molecule,is)%rxp
   ry34 = atom_list(atom3,this_molecule,is)%ryp - atom_list(atom4,this_molecule,is)%ryp
   rz34 = atom_list(atom3,this_molecule,is)%rzp - atom_list(atom4,this_molecule,is)%rzp


!   this_box = molecule_list(this_molecule,is)%which_box
!   IF (l_cubic(this_box) == .FALSE.) THEN
!   CALL Minimum_Image_Separation(this_box,rx12,ry12,rz12,rx12,ry12,rz12)
!   CALL Minimum_Image_Separation(this_box,rx32,ry32,rz32,rx32,ry32,rz32)
!   CALL Minimum_Image_Separation(this_box,rx34,ry34,rz34,rx34,ry34,rz34)
!   END IF

! Vector m is normal to the plane formed by atoms 1, 2 and 3

   mx =   ry12 * rz32 - ry32 * rz12
   my = - rx12 * rz32 + rz12 * rx32
   mz =   rx12 * ry32 - ry12 * rx32
   
! Vector n is normal to the plane formed by atoms 2, 3 and 4

   nx =   ry32 * rz34 - rz32 * ry34
   ny = - rx32 * rz34 + rz32 * rx34
   nz =   rx32 * ry34 - ry32 * rx34
   
   msq = mx * mx + my * my + mz * mz
   abs_m = DSQRT(msq)
   nsq = nx * nx + ny * ny + nz * nz
   abs_n = DSQRT(nsq)
                                                                                
! Calculate the dot product needed to calculate phi

   mdn = mx * nx + my * ny + mz * nz
   r12dn = rx12 * nx + ry12 * ny + rz12 * nz
                                                                                
! Determine what cosine of phi is and take care against numerical imprecision

   cosphi = mdn/(abs_m * abs_n)
   cosphi = MIN(cosphi,  1.0_DP)
   cosphi = MAX(cosphi, -1.0_DP)
   
! Calculate the value of dihedral angle

   phi = SIGN(DACOS(cosphi), r12dn)

 END SUBROUTINE Get_Dihedral_Angle

!********************************************************************************
 SUBROUTINE Get_Improper_Angle(this_improper,this_molecule,is,phi)
!********************************************************************************

!********************************************************************************
! This subroutine computes the value of the improper angle. Note that the routine
! is essentially identical to that used for dihedral angle evaluation.
!
! CALLED BY
! 
!        energy_routines
!              Compute_Improper_Energy
!              Compute_Molecule_Improper_Energy
!        get_com       
!********************************************************************************

   USE Type_Definitions
   USE Global_Variables

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: this_improper, is, this_molecule

   INTEGER             :: atom1, atom2, atom3, atom4, im, this_box

   REAL(DP)            :: rx12, ry12, rz12, rx32, ry32, rz32
   REAL(DP)            :: rx34, ry34, rz34, mx, my, mz, nx, ny, nz
   REAL(DP)            :: msq, nsq, mdn, abs_m, abs_n
   REAL(DP)            :: cosphi, phi, r12dn
   

  
! Get the atoms involved in the improper of interest

   atom1 = improper_list(this_improper,is)%atom1
   atom2 = improper_list(this_improper,is)%atom2
   atom3 = improper_list(this_improper,is)%atom3
   atom4 = improper_list(this_improper,is)%atom4

! Vector r12 points from atom 2 to atom 1.  Below, the components of this vector are calculated.

   rx12 = atom_list(atom1,this_molecule,is)%rxp - atom_list(atom2,this_molecule,is)%rxp
   ry12 = atom_list(atom1,this_molecule,is)%ryp - atom_list(atom2,this_molecule,is)%ryp
   rz12 = atom_list(atom1,this_molecule,is)%rzp - atom_list(atom2,this_molecule,is)%rzp
   
! Vector r32 points from atom 2 to atom 3.  Below, the components of this vector are calculated.

   rx32 = atom_list(atom3,this_molecule,is)%rxp - atom_list(atom2,this_molecule,is)%rxp
   ry32 = atom_list(atom3,this_molecule,is)%ryp - atom_list(atom2,this_molecule,is)%ryp
   rz32 = atom_list(atom3,this_molecule,is)%rzp - atom_list(atom2,this_molecule,is)%rzp

! Vector r34 points from atom 4 to atom 3. Below the components of this vector are calculated.

   rx34 = atom_list(atom3,this_molecule,is)%rxp - atom_list(atom4,this_molecule,is)%rxp
   ry34 = atom_list(atom3,this_molecule,is)%ryp - atom_list(atom4,this_molecule,is)%ryp
   rz34 = atom_list(atom3,this_molecule,is)%rzp - atom_list(atom4,this_molecule,is)%rzp

!   this_box = molecule_list(this_molecule,is)%which_box
!   IF (l_cubic(this_box) == .FALSE.) THEN
!   CALL Minimum_Image_Separation(this_box,rx12,ry12,rz12,rx12,ry12,rz12)
!   CALL Minimum_Image_Separation(this_box,rx32,ry32,rz32,rx32,ry32,rz32)
!   CALL Minimum_Image_Separation(this_box,rx34,ry34,rz34,rx34,ry34,rz34)
!   END IF

! Vector m is normal to the plane formed by atoms 1, 2 and 3

   mx =   ry12 * rz32 - ry32 * rz12
   my = - rx12 * rz32 + rz12 * rx32
   mz =   rx12 * ry32 - ry12 * rx32
   
! Vector n is normal to the plane formed by atoms 2, 3 and 4

   nx =   ry32 * rz34 - rz32 * ry34
   ny = - rx32 * rz34 + rz32 * rx34
   nz =   rx32 * ry34 - ry32 * rx34
   
   msq = mx * mx + my * my + mz * mz
   abs_m = DSQRT(msq)
   nsq = nx * nx + ny * ny + nz * nz
   abs_n = DSQRT(nsq)
                                                                                
! Calculate the dot product needed to calculate phi

   mdn = mx * nx + my * ny + mz * nz
   r12dn = rx12 * nx + ry12 * ny + rz12 * nz
                                                                                
! Determine what cosine of phi is and take care against numerical imprecision

   cosphi = mdn/(abs_m * abs_n)
   cosphi = MIN(cosphi,  1.0_DP)
   cosphi = MAX(cosphi, -1.0_DP)
   
! Calculate the value of dihedral angle

   phi = DACOS(cosphi) * SIGN(1.0_DP,r12dn)
   

 END SUBROUTINE Get_Improper_Angle
