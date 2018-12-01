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

SUBROUTINE Potential_Map_Control

  !****************************************************************************
  !
  ! Reads in all necessary control variables to run a potential map calculation
  !
  ! Called by ***
  !   main
  !
!*******************************************************************************
  USE IO_Utilities
  USE Global_Variables
  USE Potential_Map_Variables
  USE Type_Definitions
  USE File_Names
  USE Input_Routines
  USE Atoms_To_Place
  USE Angle_Dist_Pick
  USE Energy_Routines
  USE Read_Write_Checkpoint

  IMPLICIT NONE

  INTEGER ::  i, ibox
!*******************************************************************************

  CALL Copy_Inputfile
 
  ! How many species to simulate?
  CALL Get_Nspecies

  ! Load box shape, number of boxes and box type. Compute various properties of the box
  ! including the volume
  CALL Get_Box_Info

  ! Determine the type of VDW and charge interaction model to use, along with 
  ! associated parameters and the vdw mixing rule.
  CALL Get_Pair_Style

  ! Get information regarding number of molecules

  CALL Get_Molecule_Info

  CALL Create_Nonbond_Table

  CALL Create_Intra_Exclusion_Table

  ! Read in the size of the grid spacing

  CALL Get_Start_Type
  
  CALL Get_Grid_Spacing

  CALL Generate_Grid

 
  ! Now read the configuration from the input file for the zeolite
    ibox = 1

  ! Load molecular conectivity and force field paramters. Note that Get_Nspecies 
  ! must be called before this routine.  

  ! Determine the number and identity of unique atom types, and create a vdw interaction table.

  ! Obtain the temperature of the simulation
  CALL Get_Temperature_Info

  CALL Get_Rcutoff_Low

  CALL Precalculate
  CALL Read_Config(ibox)

  CALL Replicate_Unit_Cell

  CALL Generate_Neighbor_List

  CALL Zeolite_Potential_Evaluation
  stop
  ! Obtain the information about lattice file
  CALL Get_Lattice_File_Info

  ! Lattice file coordinates
  CALL Get_Lattice_Coordinates

  ! Read in super cell 


END SUBROUTINE Potential_Map_Control
