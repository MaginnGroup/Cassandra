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
SUBROUTINE NVTMC_Control
!*******************************************************************************
! Reads in all necessary control variables to run an NVT MC simulation.

  ! CALLED BY
  !
  !        main
  !
  ! CALLS
  !        Copy_Inputfile              
  !        Get_Nspecies                
  !        Get_Box_Info                
  !        Get_Pair_Style              
  !        Get_Molecule_Info           
  !        Get_Intra_Scaling           
  !        Create_Nonbond_Table        
  !        Create_Intra_Exclusion_Table
  !        Get_Seed_Info               
  !        Get_Temperature_Info        
  !        Get_Move_Probabilities     
  !        Get_CBMC_Info              
  !        Get_Simulation_Length_Info         
  !        Get_Property_Info          
  !        Average_Info               
  !        Get_Neighbor_Style         
  !        Get_Rcutoff_Low            
  !        Precalculate               
  !        Participation              
  !        Get_Bonds_Atoms_To_Place   
  !        Get_Angles_Atoms_To_Place   
  !        Angle_Distribution
  !        Get_Dihedral_Atoms_To_Place  
  !
  ! 08/07/13  : Created beta version
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

  INTEGER :: ierr,line_nbr,nbr_entries, i,j, ii,this_mol
  CHARACTER(120) :: line_string, line_array(20)
  REAL(DP) :: e_total_bond, e_total_angle, e_total_intra_nb, e_total_inter_nb
  REAL(DP) ::  E_bond, E_angle, E_dihedral, E_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq, E_inter_vdw, E_inter_qq, W_intra_vdw, W_intra_qq
  REAL(DP) :: W_inter_qq, W_inter_vdw

!*******************************************************************************
  ! Copy the input file to the logfile
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

  ! Determine the number and identity of unique atom types, and create a vdw interaction table.
  CALL Create_Nonbond_Table

  ! Create the intramolecular nonbond scaling arrays.
  CALL Create_Intra_Exclusion_Table


  CALL Get_Start_Type

  CALL Get_Seed_Info

  CALL Get_Temperature_Info

  CALL Get_Move_Probabilities

  CALL Get_CBMC_Info

  CALL Get_Simulation_Length_Info

  CALL Get_Property_Info

  CALL Get_Rcutoff_Low
  
  CALL Precalculate

  CALL Participation

  ! Obtain the list of atoms that change positions for moves that perturb intramolecular
  ! degrees of freedom

  CALL Get_Bonds_Atoms_To_Place

  ! Angle moves
  CALL Get_Angles_Atoms_To_Place
  
  ! Compute the angle probability list for each of the angles
  IF (prob_angle > 0.0_DP) THEN
     DO i = 1, nbr_boxes
        CALL Angle_Distribution(i)
     END DO
  END IF

  ! Dihedral moves
  CALL Get_Dihedral_Atoms_To_Place

END SUBROUTINE NVTMC_Control
