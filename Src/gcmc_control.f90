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

SUBROUTINE GCMC_Control

  !****************************************************************************
  ! Reads in all necessary control variables to run an NPT MC simulation.
  !
  ! Called by 
  !
  !   main
  !
  ! Revision history
  !
  !    12/10/13  : Beta release version
  !
!*******************************************************************************
  USE IO_Utilities
  USE Run_Variables
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

  CALL Get_Verbosity_Info  
  CALL Copy_Inputfile

  ! How many species to simulate?
  CALL Get_Nspecies
  WRITE(logunit,'(a30,1x,I5,/)') 'Number of species simulated: ',nspecies

  ! Load box shape, number of boxes and box type. Compute various properties of the box
  ! including the volume
  CALL Get_Box_Info

  ! Determine the type of VDW and charge interaction model to use, along with 
  ! associated parameters and the vdw mixing rule.
  CALL Get_Pair_Style

  ! Load molecular conectivity and force field paramters. Note that Get_Nspecies 
  ! must be called before this routine.  
  CALL Get_Molecule_Info

  ! Determine how intramoleclar scaling of vdw and coul interactions handled.
  CALL Get_Intra_Scaling

  ! Determine the number and identity of unique atom types, and create a vdw interaction table.
  CALL Create_Nonbond_Table

  ! Create the intramolecular nonbond scaling arrays.
  CALL Create_Intra_Exclusion_Table

  ! Seed info
  CALL Get_Seed_Info

  ! Obtain the temperature of the simulation
  CALL Get_Temperature_Info

  ! Read in the probabilities for all the moves
  CALL Get_Move_Probabilities

  ! Get the information on fugacities
  CALL Get_Fugacity_Info
    
  ! Get the product of z/omega for each of fragments
  IF (lfugacity) THEN
          CALL Get_Zig_By_Omega
  END IF 

  ! Determine the frequency with which information will be output 
  CALL Get_Simulation_Length_Info

  CALL Average_Info
  ! Properties to be output
  CALL Get_Property_Info

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

  CALL Get_CBMC_Info

  DO i=1, 1
  IF (int_vdw_sum_style(i) == vdw_mie .OR. int_vdw_sum_style(i) == vdw_mie_cut_shift) THEN
      CALL Get_Mie_Nonbond
  END IF
  END DO

END SUBROUTINE GCMC_Control
