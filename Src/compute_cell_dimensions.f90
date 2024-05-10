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

SUBROUTINE Compute_Cell_Dimensions(box_nbr)

  !********************************************************************************
  ! This subroutine analyzes the box cell matrix and determines 
  ! 1.  Lengths of the basis vectors - a,b,c (length(1->3)
  ! 2.  Cosines of the cell angles - alpha,beta,gamma (cos_angle(1->3))
  ! 3.  distance between each of the three sets of parallel cell faces face_distance(1->3)
  ! 4.  Box volume
  ! 5.  Inverse of the basis vectors (length_inv(3,3))
  !
  ! Called by
  ! 
  !    gemc_nvt_volume.f90
  !    input_routines.f90
  !    read_write_checkpoint.f90
  !    volume_change.f90
  !
  ! Revision history
  !
  !    12/10/13 : Beta Release
  !********************************************************************************
  USE Type_Definitions
  USE Global_Variables
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER :: box_nbr, i, j
  REAL(DP) :: axb1, axb2, axb3
  REAL(DP) :: bxc1, bxc2, bxc3
  REAL(DP) :: cxa1, cxa2, cxa3
  REAL(DP) :: det, inv_det
  REAL(DP) :: inv_H_T(3,3), unit_cross(3)
  LOGICAL :: change_basis, left_basis
  REAL(DP), DIMENSION(3,3) :: H, orig_H_inv
  REAL(DP), DIMENSION(3) :: orig_abc, orig_abcsq
  REAL(DP) :: ax,bx,by,cx,cy,cz

  ! Cassandra now automatically converts the cell matrix to the form
  !     required by LAMMPS since it allows better efficiency enhancement
  !     due to being an upper triangular matrix (with 3 elements guaranteed 
  !     to be zero). The cell matrix inverse is also upper triangular.
  !     This should especially help with vectorized conversions between real and
  !     fractional coordinates, since only 6 matrix elements must be stored in 6
  !     vector registers, leaving more vector registers available to store intermediates
  !     and coordinates or displacement components and possibly perform other operations
  !     without needing to load constant elements to vector registers every vector iteration 
  !     (there are 16 ymm vector registers in AVX2).
  !     The diagonal elements are also all positive in the new basis.
  ! The cell matrix is only converted if it does not already meet the new basis criteria.
  ! orig_length_inv and basis_converter in the box object are used elsewhere
  !     to transform the coordinates read from trajectory, checkpoint, or configuration files
  !     to the new basis if the cell matrix provided does not already meet the requirements.
  H = box_list(box_nbr)%length
  box_list(box_nbr)%orig_length = H
  left_basis = DOT_PRODUCT(Cross_Product(H(:,1),H(:,2)),H(:,3)) < 0.0_DP
  change_basis = &
          ANY(ABS((/ H(2,1),H(3,1),H(3,2) /))>tiny_number .OR. &
          (/ H(1,1), H(2,2), H(3,3) /) <= 0.0_DP) .OR. &
          left_basis
  IF (change_basis) THEN
          orig_H_inv(1,:) = Cross_Product(H(:,2),H(:,3))
          orig_H_inv(2,:) = Cross_Product(H(:,3),H(:,1))
          orig_H_inv(3,:) = Cross_Product(H(:,1),H(:,2))
          IF (left_basis) H(:,3) = -H(:,3)
          orig_abcsq(1) = DOT_PRODUCT(H(:,1),H(:,1))
          orig_abcsq(2) = DOT_PRODUCT(H(:,2),H(:,2))
          orig_abcsq(3) = DOT_PRODUCT(H(:,3),H(:,3))
          orig_abc = SQRT(orig_abcsq)
          ax = orig_abc(1)
          bx = DOT_PRODUCT(H(:,2),H(:,1))/orig_abc(1)
          by = SQRT(orig_abcsq(2)-bx*bx)
          cx = DOT_PRODUCT(H(:,3),H(:,1))/orig_abc(1)
          cy = (DOT_PRODUCT(H(:,2),H(:,3))-bx*cx)/by
          cz = SQRT(orig_abcsq(3) - cx*cx - cy*cy)
          H = 0.0_DP
          H(1,1) = ax
          H(1,2) = bx
          H(2,2) = by
          H(1,3) = cx
          H(2,3) = cy
          H(3,3) = cz
          box_list(box_nbr)%length = H
  END IF

  !Cell basis vector lengths
  box_list(box_nbr)%basis_length(1) = &
       SQRT( box_list(box_nbr)%length(1,1)*box_list(box_nbr)%length(1,1) + &
       box_list(box_nbr)%length(2,1)*box_list(box_nbr)%length(2,1) + &
       box_list(box_nbr)%length(3,1)*box_list(box_nbr)%length(3,1) )
 
  box_list(box_nbr)%basis_length(2) = &
       SQRT( box_list(box_nbr)%length(1,2)*box_list(box_nbr)%length(1,2) + &
       box_list(box_nbr)%length(2,2)*box_list(box_nbr)%length(2,2) + &
       box_list(box_nbr)%length(3,2)*box_list(box_nbr)%length(3,2) )
  
  box_list(box_nbr)%basis_length(3) = &
       SQRT( box_list(box_nbr)%length(1,3)*box_list(box_nbr)%length(1,3) + &
       box_list(box_nbr)%length(2,3)*box_list(box_nbr)%length(2,3) + &
       box_list(box_nbr)%length(3,3)*box_list(box_nbr)%length(3,3) )


  ! Compute cosine of cell angles
  box_list(box_nbr)%cos_angle(1) = ( box_list(box_nbr)%length(1,2)*box_list(box_nbr)%length(1,3) + &
       box_list(box_nbr)%length(2,2)*box_list(box_nbr)%length(2,3) + &
       box_list(box_nbr)%length(3,2)*box_list(box_nbr)%length(3,3) ) / &
       ( box_list(box_nbr)%basis_length(2)* box_list(box_nbr)%basis_length(3) )


  box_list(box_nbr)%cos_angle(2) = ( box_list(box_nbr)%length(1,1)*box_list(box_nbr)%length(1,3) + &
       box_list(box_nbr)%length(2,1)*box_list(box_nbr)%length(2,3) + &
       box_list(box_nbr)%length(3,1)*box_list(box_nbr)%length(3,3) ) / &
       ( box_list(box_nbr)%basis_length(1)* box_list(box_nbr)%basis_length(3) )

  box_list(box_nbr)%cos_angle(3) = ( box_list(box_nbr)%length(1,1)*box_list(box_nbr)%length(1,2) + &
       box_list(box_nbr)%length(2,1)*box_list(box_nbr)%length(2,2) + &
       box_list(box_nbr)%length(3,1)*box_list(box_nbr)%length(3,2) ) / &
       ( box_list(box_nbr)%basis_length(1)* box_list(box_nbr)%basis_length(2) )

  ! Compute cross products
  axb1 = box_list(box_nbr)%length(2,1)*box_list(box_nbr)%length(3,2) - &
       box_list(box_nbr)%length(3,1)*box_list(box_nbr)%length(2,2)
  axb2 = box_list(box_nbr)%length(3,1)*box_list(box_nbr)%length(1,2) - &
       box_list(box_nbr)%length(1,1)*box_list(box_nbr)%length(3,2)
  axb3 = box_list(box_nbr)%length(1,1)*box_list(box_nbr)%length(2,2) - &
       box_list(box_nbr)%length(2,1)*box_list(box_nbr)%length(1,2)

  bxc1 = box_list(box_nbr)%length(2,2)*box_list(box_nbr)%length(3,3) - &
       box_list(box_nbr)%length(3,2)*box_list(box_nbr)%length(2,3)
  bxc2 = box_list(box_nbr)%length(3,2)*box_list(box_nbr)%length(1,3) - &
       box_list(box_nbr)%length(1,2)*box_list(box_nbr)%length(3,3)
  bxc3 = box_list(box_nbr)%length(1,2)*box_list(box_nbr)%length(2,3) - &
       box_list(box_nbr)%length(2,2)*box_list(box_nbr)%length(1,3)

  cxa1 = box_list(box_nbr)%length(2,3)*box_list(box_nbr)%length(3,1) - &
       box_list(box_nbr)%length(2,1)*box_list(box_nbr)%length(3,3)
  cxa2 = box_list(box_nbr)%length(1,1)*box_list(box_nbr)%length(3,3) - &
       box_list(box_nbr)%length(3,1)*box_list(box_nbr)%length(1,3)
  cxa3 = box_list(box_nbr)%length(2,1)*box_list(box_nbr)%length(1,3) - &
       box_list(box_nbr)%length(1,1)*box_list(box_nbr)%length(2,3)

  ! Compute cell volume
  box_list(box_nbr)%volume = ABS(box_list(box_nbr)%length(1,1)*bxc1 + &
       box_list(box_nbr)%length(2,1)*bxc2 + box_list(box_nbr)%length(3,1)*bxc3)
  
  ! Compute distance between cell faces
  box_list(box_nbr)%face_distance(1) = box_list(box_nbr)%volume / &
       SQRT(bxc1*bxc1 + bxc2*bxc2 + bxc3*bxc3)

  box_list(box_nbr)%face_distance(2) = box_list(box_nbr)%volume / &
       SQRT(cxa1*cxa1 + cxa2*cxa2 + cxa3*cxa3)

  box_list(box_nbr)%face_distance(3) = box_list(box_nbr)%volume / &
       SQRT(axb1*axb1 + axb2*axb2 + axb3*axb3)

  ! Now compute the inverse of the cell basis matrix

  ! First compute the adjoint of the cell basis matrix
  box_list(box_nbr)%length_inv(1,1) =  &
       box_list(box_nbr)%length(2,2) * box_list(box_nbr)%length(3,3) - &
       box_list(box_nbr)%length(2,3) * box_list(box_nbr)%length(3,2)

  box_list(box_nbr)%length_inv(2,1) =  &
       box_list(box_nbr)%length(3,1) * box_list(box_nbr)%length(2,3) - &
       box_list(box_nbr)%length(2,1) * box_list(box_nbr)%length(3,3)

  box_list(box_nbr)%length_inv(3,1) =  &
       box_list(box_nbr)%length(2,1) * box_list(box_nbr)%length(3,2) - &
       box_list(box_nbr)%length(2,2) * box_list(box_nbr)%length(3,1)

  box_list(box_nbr)%length_inv(1,2) =  &
       box_list(box_nbr)%length(3,2) * box_list(box_nbr)%length(1,3) - &
       box_list(box_nbr)%length(1,2) * box_list(box_nbr)%length(3,3)

  box_list(box_nbr)%length_inv(2,2) = &
       box_list(box_nbr)%length(1,1)*box_list(box_nbr)%length(3,3) - &
       box_list(box_nbr)%length(3,1)*box_list(box_nbr)%length(1,3)

  box_list(box_nbr)%length_inv(3,2) = &
       box_list(box_nbr)%length(3,1)*box_list(box_nbr)%length(1,2) - &
       box_list(box_nbr)%length(1,1)*box_list(box_nbr)%length(3,2)

  box_list(box_nbr)%length_inv(1,3) = & 
       box_list(box_nbr)%length(1,2)*box_list(box_nbr)%length(2,3) - &
       box_list(box_nbr)%length(2,2)*box_list(box_nbr)%length(1,3)
 
  box_list(box_nbr)%length_inv(2,3) = & 
       box_list(box_nbr)%length(2,1)*box_list(box_nbr)%length(1,3) - &
       box_list(box_nbr)%length(1,1)*box_list(box_nbr)%length(2,3)
  
  box_list(box_nbr)%length_inv(3,3) = &
       box_list(box_nbr)%length(1,1)*box_list(box_nbr)%length(2,2) - &
       box_list(box_nbr)%length(2,1)*box_list(box_nbr)%length(1,2)

  ! Calculate the determinant of in_cell matrix. This determinant is also the volume
  ! so check 
  det = box_list(box_nbr)%length(1,1)*box_list(box_nbr)%length_inv(1,1) + &
       box_list(box_nbr)%length(1,2)*box_list(box_nbr)%length_inv(2,1) + &
       box_list(box_nbr)%length(1,3)*box_list(box_nbr)%length_inv(3,1)


  IF(ABS(det - box_list(box_nbr)%volume) > tiny_number) THEN
     err_msg = ""
     err_msg(1) = 'Determinant not equal box volume'
     err_msg(2) = 'For box'
     err_msg(3) = Int_To_String(box_nbr)
     CALL Clean_Abort(err_msg,'Compute_Cell_Dimensions')
  ENDIF

  inv_det = 1.0_DP / det

  ! The adjoint divided by the determinant is the inverse of the cell basis matrix
  box_list(box_nbr)%length_inv = box_list(box_nbr)%length_inv * inv_det
  IF (change_basis) THEN
          orig_H_inv = orig_H_inv * inv_det
          box_list(box_nbr)%orig_length_inv = orig_H_inv
          box_list(box_nbr)%basis_converter = MATMUL(H,orig_H_inv)
  ELSE
          box_list(box_nbr)%orig_length_inv = box_list(box_nbr)%length_inv
  END IF
  box_list(box_nbr)%basis_changed = change_basis

  ! Compute half of the box length to be used in Fold_Molecule

  box_list(box_nbr)%hlength(1,1) =0.5_DP * box_list(box_nbr)%basis_length(1)
  box_list(box_nbr)%hlength(2,2) =0.5_DP * box_list(box_nbr)%basis_length(2)
  box_list(box_nbr)%hlength(3,3) =0.5_DP * box_list(box_nbr)%basis_length(3)

  ! Compute face distance of transpose of inverse of cell matrix
  ! This is used for reciprocal lattice vector setup
  inv_H_T = TRANSPOSE(box_list(box_nbr)%length_inv)
  unit_cross = Cross_Product(inv_H_T(:,2),inv_H_T(:,3))
  unit_cross = unit_cross/SQRT(DOT_PRODUCT(unit_cross,unit_cross))
  box_list(box_nbr)%invT_face_distance(1) = ABS(DOT_PRODUCT(inv_H_T(:,1),unit_cross))
  unit_cross = Cross_Product(inv_H_T(:,3),inv_H_T(:,1))
  unit_cross = unit_cross/SQRT(DOT_PRODUCT(unit_cross,unit_cross))
  box_list(box_nbr)%invT_face_distance(2) = ABS(DOT_PRODUCT(inv_H_T(:,2),unit_cross))
  unit_cross = Cross_Product(inv_H_T(:,1),inv_H_T(:,2))
  unit_cross = unit_cross/SQRT(DOT_PRODUCT(unit_cross,unit_cross))
  box_list(box_nbr)%invT_face_distance(3) = ABS(DOT_PRODUCT(inv_H_T(:,3),unit_cross))

      
END SUBROUTINE Compute_Cell_Dimensions

