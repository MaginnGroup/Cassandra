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
SUBROUTINE MCF_Control
  !*******************************************************************************
  ! Generates the mcf files needed for all the fragments in a species
  !
  ! CALLED BY
  !     
  !        main
  !
  ! CALLS
  !
  !        Copy_Inputfile
  !        Get_Nspecies
  !        Get_Box_Info
  !        Get_Pair_Style
  !        Get_Molecule_Info
  !        Get Participation
  !
  !  08/07/13 : Created beta version
  !
  !*******************************************************************************
  USE IO_Utilities
  USE Global_Variables
  USE Type_Definitions
  USE Input_Routines
  USE Atoms_To_Place

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

  CALL Get_Box_Info

  CALL Get_Pair_Style

  ! Load molecular conectivity and force field paramters. Note that Get_Nspecies 
  ! must be called before this routine.  
  CALL Get_Molecule_Info

  ! Call this routine to set up the mcf files
  CALL Participation


  WRITE(*,*) 'Finished generating fragments'
  WRITE(*,*) 'Exiting'
  WRITE(*,*)
  WRITE(logunit,*)
  WRITE(logunit,*) 'Finished generating fragments'
  WRITE(logunit,*) 'Exiting'
  STOP

END SUBROUTINE MCF_Control
