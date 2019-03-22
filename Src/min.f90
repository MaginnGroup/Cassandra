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


!****************************************************************************************
!
!  Contains the following routines
!
!  Minimum_Image_Separtion 
!  Apply_PBC_Int
!  Fold_Molecule
!
!  08/07/13  : Beta version created
!
!****************************************************************************************

SUBROUTINE Minimum_Image_Separation(ibox,rxijp,ryijp,rzijp,rxij,ryij,rzij)

  ! Passed a box number and the parent coordinate Cartesian separation
  ! distances, this routine returns the "minimum image" coordinate 
  ! separation distances (i.e. the closest image separation distance of
  ! atoms i and j). See Allen and Tildesley, page 28.  
  !
  !
  ! 
  !---------------------------------------------------------------------------------------------- 
  USE Global_Variables
  USE Type_Definitions

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ibox
  REAL(DP), INTENT(IN) :: rxijp,ryijp,rzijp
  REAL(DP), INTENT(OUT) :: rxij,ryij,rzij

  REAL(DP), DIMENSION(3) :: temp_vec
  REAL(DP), DIMENSION(3) :: hbox   
 
  !---------------------------------------------------------------------------------------------- 

  IF(l_cubic(ibox)) THEN

     rxij = rxijp
     ryij = ryijp
     rzij = rzijp

     IF (rxijp.gt.box_list(ibox)%hlength(1,1)) THEN
        rxij = rxijp-box_list(ibox)%length(1,1)
     ELSEIF (rxijp.lt.-box_list(ibox)%hlength(1,1)) THEN
        rxij = rxijp+box_list(ibox)%length(1,1)
     ENDIF

     IF (ryijp.gt.box_list(ibox)%hlength(2,2)) THEN
        ryij = ryijp-box_list(ibox)%length(2,2)
     ELSEIF (ryijp.lt.-box_list(ibox)%hlength(2,2)) THEN
        ryij = ryijp+box_list(ibox)%length(2,2)
     ENDIF

     IF (rzijp.gt.box_list(ibox)%hlength(3,3)) THEN
        rzij = rzijp-box_list(ibox)%length(3,3)
     ELSEIF (rzijp.lt.-box_list(ibox)%hlength(3,3)) THEN
        rzij = rzijp+box_list(ibox)%length(3,3)
     ENDIF
     
!     rxij = rxijp - box_list(ibox)%length(1,1)*REAL(NINT(rxijp/box_list(ibox)%length(1,1)),DP)
!     ryij = ryijp - box_list(ibox)%length(2,2)*REAL(NINT(ryijp/box_list(ibox)%length(2,2)),DP)
!     rzij = rzijp - box_list(ibox)%length(3,3)*REAL(NINT(rzijp/box_list(ibox)%length(3,3)),DP)

  ELSE

     ! Always use cell_matrix convention so this routine works for anything

     !First convert the parent coordinates from the Cartesian to fractional
     !coordinate system
     temp_vec(1) = box_list(ibox)%length_inv(1,1)*rxijp + &
       box_list(ibox)%length_inv(1,2)*ryijp +          &
       box_list(ibox)%length_inv(1,3)*rzijp

     temp_vec(2) = box_list(ibox)%length_inv(2,1)*rxijp + &
       box_list(ibox)%length_inv(2,2)*ryijp +          &
       box_list(ibox)%length_inv(2,3)*rzijp

     temp_vec(3) = box_list(ibox)%length_inv(3,1)*rxijp + &
       box_list(ibox)%length_inv(3,2)*ryijp +          &
       box_list(ibox)%length_inv(3,3)*rzijp

     !Apply periodic boundary conditions to the fractional distances.
     ! Recall NINT rounds and does not truncate.
     temp_vec(1) = temp_vec(1) - REAL(NINT(temp_vec(1)),DP)
     temp_vec(2) = temp_vec(2) - REAL(NINT(temp_vec(2)),DP)
     temp_vec(3) = temp_vec(3) - REAL(NINT(temp_vec(3)),DP)
  
     !Convert back to Cartesian coordinates and return the results as
     !the child coordinate separations
     rxij = box_list(ibox)%length(1,1)*temp_vec(1) + &
       box_list(ibox)%length(1,2)*temp_vec(2) +   &
       box_list(ibox)%length(1,3)*temp_vec(3)

     ryij = box_list(ibox)%length(2,1)*temp_vec(1) + &
       box_list(ibox)%length(2,2)*temp_vec(2) +   &
       box_list(ibox)%length(2,3)*temp_vec(3)

     rzij = box_list(ibox)%length(3,1)*temp_vec(1) + &
       box_list(ibox)%length(3,2)*temp_vec(2) +   &
       box_list(ibox)%length(3,3)*temp_vec(3)

  END IF

END SUBROUTINE Minimum_Image_Separation

SUBROUTINE Apply_PBC_Anint(ibox,rxijp,ryijp,rzijp,rxij,ryij,rzij)

  USE Global_Variables
  USE Type_Definitions

  IMPLICIT NONE
 
  INTEGER, INTENT(IN) :: ibox
  REAL(DP), INTENT(IN) :: rxijp, ryijp, rzijp
  REAL(DP), INTENT(OUT) :: rxij, ryij, rzij
  REAL(DP) :: fracx,fracy,fracz

  IF (l_cubic(ibox)) THEN


     rxij = rxijp - box_list(ibox)%length(1,1)* &
          REAL(ANINT( rxijp / box_list(ibox)%length(1,1)), DP)

     ryij = ryijp - box_list(ibox)%length(2,2)* &
          REAL(ANINT( rxijp / box_list(ibox)%length(2,2)), DP)

     rzij = rzijp - box_list(ibox)%length(3,3)* &
          REAL(ANINT( rxijp / box_list(ibox)%length(3,3)), DP)

  ELSE
     

     CALL Cartesian_To_Fractional(rxijp,ryijp,rzijp,fracx,fracy,fracz,ibox)

     !First convert the parent coordinates from the Cartesian to fractional
     !coordinate system

     !Apply periodic boundary conditions to the fractional distances.
     ! Recall NINT rounds and does not truncate.
     fracx = fracx - REAL(NINT(fracx),DP)
     fracy = fracy - REAL(NINT(fracy),DP)
     fracz = fracz - REAL(NINT(fracz),DP)


     !Convert back to Cartesian coordinates and return the results as
     !the child coordinate separations

     CALL Fractional_To_Cartesian(fracx,fracy,fracz,rxij,ryij,rzij,ibox)

  END IF

END SUBROUTINE Apply_PBC_Anint

SUBROUTINE Fold_Molecule(alive,is,this_box)

  USE Global_Variables
  USE Type_Definitions

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: alive, is, this_box

  REAL(DP) :: dx, dy, dz

  REAL(DP) :: thisx,thisy,thisz
  REAL(DP) :: frac_comx, frac_comy, frac_comz
  INTEGER :: i

  IF (l_cubic(this_box)) THEN
     
     IF(molecule_list(alive,is)%xcom .GT. box_list(this_box)%hlength(1,1)) THEN
        molecule_list(alive,is)%xcom = &
             molecule_list(alive,is)%xcom - box_list(this_box)%length(1,1)
        atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp - box_list(this_box)%length(1,1)
        
     ELSE IF(molecule_list(alive,is)%xcom .LT. -box_list(this_box)%hlength(1,1)) THEN
        molecule_list(alive,is)%xcom = &
             molecule_list(alive,is)%xcom + box_list(this_box)%length(1,1)
        atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp + box_list(this_box)%length(1,1)
     END IF

     IF(molecule_list(alive,is)%ycom .GT. box_list(this_box)%hlength(2,2)) THEN
        molecule_list(alive,is)%ycom = &
             molecule_list(alive,is)%ycom - box_list(this_box)%length(2,2)
        atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp - box_list(this_box)%length(2,2)

     ELSE IF(molecule_list(alive,is)%ycom .LT. -box_list(this_box)%hlength(2,2)) THEN

        molecule_list(alive,is)%ycom = &
             molecule_list(alive,is)%ycom + box_list(this_box)%length(2,2)
        atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp + box_list(this_box)%length(2,2)

     END IF

     IF(molecule_list(alive,is)%zcom .GT. box_list(this_box)%hlength(3,3)) THEN

        molecule_list(alive,is)%zcom = &
             molecule_list(alive,is)%zcom - box_list(this_box)%length(3,3)
        atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp - box_list(this_box)%length(3,3)
        
     ELSE IF(molecule_list(alive,is)%zcom .LT. -box_list(this_box)%hlength(3,3)) THEN
        
        molecule_list(alive,is)%zcom = &
             molecule_list(alive,is)%zcom + box_list(this_box)%length(3,3)
        atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp + box_list(this_box)%length(3,3)
     END IF


  ELSE

  thisx = molecule_list(alive,is)%xcom
  thisy = molecule_list(alive,is)%ycom
  thisz = molecule_list(alive,is)%zcom

  CALL Cartesian_To_Fractional(thisx,thisy,thisz,frac_comx, frac_comy, frac_comz, this_box)

     IF(frac_comx .GT. 0.5) THEN


        CALL Fold_Molecule_In_Fractional_Coords(alive,is,this_box)  


     ELSE IF(frac_comx .LT. -0.5) THEN

        CALL Fold_Molecule_In_Fractional_Coords(alive,is, this_box)

     END IF



     IF(frac_comy .GT. 0.5) THEN
        
        CALL Fold_Molecule_In_Fractional_Coords(alive,is,this_box)  

     ELSE IF(frac_comy .LT. -0.5) THEN
        CALL Fold_Molecule_In_Fractional_Coords(alive,is,this_box)  


     END IF

     IF(frac_comz .GT. 0.5) THEN

        CALL Fold_Molecule_In_Fractional_Coords(alive,is,this_box)  

     ELSE IF(frac_comz .LT. -0.5) THEN

        CALL Fold_Molecule_In_Fractional_Coords(alive,is,this_box)  
     END IF



!  ELSE

 !    err_msg = ''
  !   err_msg(1) = 'Folding attempted for box shape other than cubic or'
   !  err_msg(2) = 'slit_pore geometry'
    ! err_msg(3) = 'Folding of molecules is not supported for the the shape of the box'
     !CALL Clean_Abort(err_msg,'Fold_Molecule')

  END IF

END SUBROUTINE

SUBROUTINE Cartesian_To_Fractional(rx,ry,rz,sx,sy,sz,ibox)
USE Global_Variables
USE Type_Definitions
REAL(DP), INTENT(IN) :: rx,ry,rz
REAL(DP), INTENT(OUT) :: sx,sy,sz

     sx = box_list(ibox)%length_inv(1,1)*rx + &
box_list(ibox)%length_inv(1,2)*ry + &
box_list(ibox)%length_inv(1,3)*rz

     sy = box_list(ibox)%length_inv(2,1)*rx + &
box_list(ibox)%length_inv(2,2)*ry + &
box_list(ibox)%length_inv(2,3)*rz

     sz = box_list(ibox)%length_inv(3,1)*rx + &
box_list(ibox)%length_inv(3,2)*ry + &
box_list(ibox)%length_inv(3,3)*rz


END SUBROUTINE


SUBROUTINE Fractional_To_Cartesian(sx,sy,sz,rx,ry,rz,ibox)
USE Global_Variables
USE Type_Definitions
REAL(DP), INTENT(OUT) :: rx,ry,rz
REAL(DP), INTENT(IN) :: sx,sy,sz

     rx = box_list(ibox)%length(1,1)*sx + &
box_list(ibox)%length(1,2)*sy + &
box_list(ibox)%length(1,3)*sz

     ry = box_list(ibox)%length(2,1)*sx + &
box_list(ibox)%length(2,2)*sy + &
box_list(ibox)%length(2,3)*sz

     rz = box_list(ibox)%length(3,1)*sx + &
box_list(ibox)%length(3,2)*sy + &
box_list(ibox)%length(3,3)*sz

END SUBROUTINE

SUBROUTINE Fold_Molecule_In_Fractional_Coords(this_im, this_is,this_box)
USE Global_Variables
USE Type_Definitions
INTEGER, INTENT(IN) :: this_im, this_is
REAL(DP) :: thisx,thisy,thisz
REAL(DP) :: thisx2,thisy2,thisz2
INTEGER :: i


     thisx=molecule_list(this_im,this_is)%xcom
     thisy=molecule_list(this_im,this_is)%ycom
     thisz=molecule_list(this_im,this_is)%zcom

     CALL Apply_PBC_Anint(this_box,thisx,thisy,thisz,thisx2,thisy2,thisz2)
     
     molecule_list(this_im,this_is)%xcom = thisx2
     molecule_list(this_im,this_is)%ycom = thisy2
     molecule_list(this_im,this_is)%zcom = thisz2

     DO i = 1, natoms(this_is)

        thisx=atom_list(i,this_im,this_is)%rxp
        thisy=atom_list(i,this_im,this_is)%ryp
        thisz=atom_list(i,this_im,this_is)%rzp

        CALL Apply_PBC_Anint(this_box,thisx,thisy,thisz,thisx2,thisy2,thisz2)

        atom_list(i,this_im,this_is)%rxp = thisx2
        atom_list(i,this_im,this_is)%ryp = thisy2
        atom_list(i,this_im,this_is)%rzp = thisz2

     END DO

END SUBROUTINE
