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


MODULE Volume

  !***********************************************************************
  !
  ! This module contains collection of routines that can be accessed
  ! by routines that perform volume changes.
  !
  ! Save_Cartesian_Coordinates_Box
  !
  !      Storages of COMs and atomic positions of all the molecules of the
  !      the input box.
  !
  ! Reset_Cartesian_Coordinates_Box
  !
  !      Restores COM and atomic positions of all the molecules of the input
  !      box.
  !
  ! Scale_COM_Cartesian
  !
  !      Convert the Cartesian coordinates  COM of all the molecules of the 
  !      input box to the fractional coordinates
  !
  ! CALLED BY:
  !
  !        gemc_nvt_volume.f90
  !
  ! CALLS
  ! 
  !        Save_Old_Cartesian_Coordinates
  !        Revert_Old_Cartesian_Coordinates
  !
  ! 08/12/13 : Created the beta version
  !************************************************************************

  USE Global_Variables
  USE Type_Definitions
  
  IMPLICIT NONE

CONTAINS

  SUBROUTINE Save_Cartesian_Coordinates_Box(this_box)
    !************************************************************************
    ! store the old configuration of all atoms and COMs of molecules of the
    ! input box 
    !*************************************************************************

    INTEGER, INTENT(IN) :: this_box

    ! Local variables

    INTEGER :: is, im, lm

    DO is = 1, nspecies
       
       DO im = 1, nmols(is,this_box)
          
          lm = locate(im,is,this_box)
          
          ! store the unperturbed position of this molecule. Note that
          ! this call to the routine will save Cartesian coordinates
          ! as well as COMs of this molecule
          
          CALL Save_Old_Cartesian_Coordinates(lm,is)
          
       END DO
       
    END DO

  END SUBROUTINE Save_Cartesian_Coordinates_Box

!*************************************************************************

  SUBROUTINE Reset_Cartesian_Coordinates_Box(this_box)


    INTEGER, INTENT(IN) :: this_box
    
    ! Local variables

    INTEGER :: is, im, lm

    DO is = 1, nspecies
       DO im = 1, nmols(is,this_box)
          lm = locate(im,is,this_box)
          IF (molecule_list(lm,is)%live) THEN
             CALL Revert_Old_Cartesian_Coordinates(lm,is)
          END IF
       END DO
    END DO

  End SUBROUTINE Reset_Cartesian_Coordinates_Box

!****************************************************************************

  SUBROUTINE Scale_COM_Cartesian(this_box,box_list_old)

    ! Scales the cartesian coordinates and COM of all the molecules
    ! in the input box such that intramolecular DOFs do not change
    ! This is achieved by keeping the fractional coordinates of the COM
    ! the same before and after the move. The old cell basis vector is
    ! used for this purpose

    INTEGER, INTENT(IN) :: this_box
    TYPE(Box_Class), INTENT(IN) :: box_list_old

    ! Local variables

    INTEGER :: is, im, lm, i
    REAL(DP) :: s(3)



    DO is = 1, nspecies
       
       DO im = 1, nmols(is,this_box)
          
          lm = locate(im,is,this_box)
          
          IF (molecule_list(lm,is)%live) THEN
             
             ! obtain the new coordinates of the COM for this molecule
             
             ! first determine the fractional coordinate
             
             DO i = 1,3
                
                s(i) = box_list_old%length_inv(i,1) * molecule_list(lm,is)%xcom + &
                     box_list_old%length_inv(i,2) * molecule_list(lm,is)%ycom + &
                     box_list_old%length_inv(i,3) * molecule_list(lm,is)%zcom
             END DO
             
             
             ! now obtain the new positions of COMs
             
             
             molecule_list(lm,is)%xcom = box_list(this_box)%length(1,1) * s(1) &
                  + box_list(this_box)%length(1,2) * s(2) + &
                    box_list(this_box)%length(1,3) * s(3)
             
             molecule_list(lm,is)%ycom = box_list(this_box)%length(2,1) * s(1) &
                  + box_list(this_box)%length(2,2) * s(2) + &
                    box_list(this_box)%length(2,3) * s(3)
             
             molecule_list(lm,is)%zcom = box_list(this_box)%length(3,1) * s(1) &
                  + box_list(this_box)%length(3,2) * s(2) + &
                    box_list(this_box)%length(3,3) * s(3)
             
             ! Obtain the new positions of atoms in this molecule
             
             atom_list(:,lm,is)%rxp = atom_list(:,lm,is)%rxp + &
                  molecule_list(lm,is)%xcom - molecule_list(lm,is)%xcom_old
             
             atom_list(:,lm,is)%ryp = atom_list(:,lm,is)%ryp + &
                  molecule_list(lm,is)%ycom - molecule_list(lm,is)%ycom_old
             
             atom_list(:,lm,is)%rzp = atom_list(:,lm,is)%rzp + &
                  molecule_list(lm,is)%zcom - molecule_list(lm,is)%zcom_old
             
          END IF
          
       END DO
       
    END DO
    
  END SUBROUTINE Scale_COM_Cartesian
!************************************************************************

END MODULE Volume
