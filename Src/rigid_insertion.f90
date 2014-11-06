
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

SUBROUTINE Rigid_Insertion

  !******************************************************************************
  ! The subroutine performs a rigid insertion of a molecule into the simulation
  ! cell. This is useful when we would like to study the absorption of rigid
  ! gas molecules in IL
  !******************************************************************************
  ! The subroutine will randomly pick a species and if it is found to be a sorbate
  ! species then an attempt will be made to insert a molecule of this species in the
  ! simulation box. The placement of the molecule will follow the strategy as.
  ! 
  ! We will place randomly the COM of the species in the simulation box. In the
  ! first version of the code cubic box is assumed for the random insertion.
  ! We will obtain a random orientation of the molecule about COM and then translate
  ! the molecule at randomly chosen location of the box.
  ! The energy change due to the insertion will be computed and we will accept
  ! or reject the move based on the acceptance criteria for a given ensemble of interest.
  !
  !
  ! Revision History
  !
  ! 12/10/13 : Beta release
  !
  !***************************************************************************

  
  USE Type_Definitions
  USE Run_Variables

  IMPLICIT NONE

  
  ! We will choose a box to insert a molecule

  this_box = INT ( rranf() * nbr_boxes ) + 1
  
  DO WHILE ( .TRUE.)

     ! Choose a species at random

     is = INT ( rranf() * nspecies ) + 1

     IF (species_list(is)%species_type == 'SORBATE') EXIT

  END DO

  ! We determine the total number of molecules of type is. This number
  ! will be used to check on the bounds of array that depend on molecules.

  nmolecules_is = 0

  DO i = 1, nmolecules(is)

     alive = locate(i,is)

     IF ( molecule_list(alive,is)%alive ) THEN

        nmolecules_is = nmolecules_is + 1
        
     END IF

  END DO

  ! Provide an array bound check

  IF ( nmolecules_is + 1 .GE. MAXVAL(nmolecules) ) THEN

     ! We do not have sufficient memory to store the information for this molecule
     ! Abort simulation with a flag

     err_msg = ""
     err_msg(1) = 'Array bounds exceeded while attempting to insert a molecule of type'
     err_msg(2) = TRIM(Int_to_String(is))
     CALL Clean_Abort(err_msg,'Rigid_Insertion')

  END IF

  ! Assign a locate number for this molecule

  IF ( molecule_list(nmolecules_is+1,is) = 0 ) THEN
     locate(nmolecules_is+1,is) = nmolecules_is
     ! otherwise we will use the locate number of a previously deleted molecule that has 
     ! been moved to the end of the array.
  END IF

  alive = locate(nmolecules_is+1,is)
  
  ! in the end if the move is accepted we will make this molecule alive and assign it
  ! to this_box
  ! i. e. molecule_list(this_l

  ! Randomly insert the COM in the simulation box.

  IF ( box_list(this_box)%shape = 'CUBIC' ) THEN

     molecule_list(alive,is)%xcom = (1.0_DP - 0.5_DP * rranf()) * box_list(this_box)%length(1,1)
     molecule_list(alive,is)%ycom = (1.0_DP - 0.5_DP * rranf()) * box_list(this_box)%length(2,2)
     molecule_list(alive,is)%zcom = (1.0_DP - 0.5_DP * rranf()) * box_list(this_box)%length(3,3)

  END IF

  ! Apply random rotation to the molecule 
  ! Translate the molecule at the COM.
  ! Compute the new positions.
  ! Calculate the change in energy
  ! Accept or reject the move.
  ! also, update various arrays that depend on the number of molecules

END SUBROUTINE Rigid_Insertion
