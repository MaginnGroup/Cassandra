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
SUBROUTINE Fragment_Control
  !*******************************************************************************
  ! Reads in all the necessary information to carry out sampling of fragments
  ! using atom displacemnt 
  !
  !
  ! CALLED BY
  !        main
  !
  ! CALLS
  !
  !       Copy_Inputfile
  !       Get_Nspecies
  !       Get_Box_Info
  !       Get_Pair_Style
  !       Get_Molecule_info
  !       Get_Intra_Scaling
  !       Create_Nonbond_Table
  !       Create_Intra_Exclusion_Table
  !       Get_Seed_Info
  !       Get_Temperature_Info
  !       Get_Move_Probabilities
  !       Get_File_Info
  !       Get_Simulation_Length_Info
  !       Precalculate
  !       Participation
  !
  ! 08/07/13 : Created beta version
!********************************************************************************

  USE IO_Utilities
  USE Global_Variables
  USE Input_Routines

  IMPLICIT NONE

  ! Copy input file to the log file

  CALL Copy_Inputfile
  
  ! Number of species to simulate

  CALL Get_Nspecies

  CALL Get_Box_Info

  CALL Get_Pair_Style

  CALL Get_Molecule_Info

  ! Determine the number and identity of unique atom types, and create a vdw interaction table.
  CALL Create_Nonbond_Table
  
  ! Create the intramolecular nonbond scaling arrays.
  CALL Create_Intra_Exclusion_Table  ! Obtain information about the molecules

  CALL Get_Start_Type

  ! Random initial seed
  CALL Get_Seed_Info

  ! Temperature
  CALL Get_Temperature_Info

  ! Probabilities
  CALL Get_Move_Probabilities
  
  ! Frequency info
  CALL Get_Simulation_Length_Info

  ! Get the file names for each of the fragments
  CALL Get_File_Info

  CALL Precalculate

  CALL Participation


END SUBROUTINE Fragment_Control
