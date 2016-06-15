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

SUBROUTINE GEMC_Control

  !****************************************************************************
  !
  ! Reads in all necessary control variables to run a GEMC simulation
  !
  ! Called by ***
  !   main
  !
  ! Revision history
  !
  !   12/10/13 : Beta Release
!*******************************************************************************
  USE IO_Utilities
  USE Global_Variables
  USE Type_Definitions
  USE File_Names
  USE Input_Routines
  USE Atoms_To_Place
  USE Angle_Dist_Pick
  USE Energy_Routines

  IMPLICIT NONE

  INTEGER ::  i
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

  ! Load molecular conectivity and force field paramters. Note that Get_Nspecies 
  ! must be called before this routine.  
  CALL Get_Molecule_Info

  ! If the simulation contains an intermediate box for transfer then
  ! read the number of ideal gas molecules in that box

!  IF (int_sim_type == sim_gemc_ig) THEN

 !    CALL Get_IGAS_Num

  !END IF

  ! Determine the number and identity of unique atom types, and create a vdw interaction table.
  CALL Create_Nonbond_Table

  ! Create the intramolecular nonbond scaling arrays.
  CALL Create_Intra_Exclusion_Table

  ! Start_Type
  CALL Get_Start_Type

  ! Seed info
  CALL Get_Seed_Info
 
  ! Obtain the temperature of the simulation
  CALL Get_Temperature_Info

  IF (sim_type == 'GEMC_NPT') THEN
     ! Obtain the pressure of the simulation box
     CALL Get_Pressure_Info
   
  END IF
  ! Read in the probabilities for all the moves
  CALL Get_Move_Probabilities

  CALL Get_CBMC_Info


  ! Determine the frequency with which information will be output 
  CALL Get_Simulation_Length_Info

  ! Properties to be output
  CALL Get_Property_Info

  ! Get information on the averages
  CALL Get_Average_Info

  CALL Get_Rcutoff_Low

  CALL Precalculate

  ! Determine the connectivity information for each atom such as how many bonds
  ! it participates in, what atoms are connected to the atom, similar information
  ! for angles and dihedrals
  CALL Participation

  ! Determine the atoms that need to be moved if one of the ends of a bond is
  ! perturbed
  CALL Get_Bonds_Atoms_To_Place
 
  ! Determine what angles a given atom participates and how many such angles exist
  CALL Get_Angles_Atoms_To_Place

  ! Obtain the probability for angle distribution if angle distortion moves are
  ! attempted

  IF (prob_angle > 0.0_DP) THEN
     DO i = 1, nbr_boxes
        CALL Angle_Distribution(i)
     END DO
  END IF
  
  ! Determine what dihedral angles a given atom participates and how many such
  ! angles exist
  CALL Get_Dihedral_Atoms_To_Place  

END SUBROUTINE GEMC_Control
