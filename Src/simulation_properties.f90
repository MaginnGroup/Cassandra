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

MODULE Simulation_Properties

  !********************************************************************************
  ! This module contains a collection of routines that can be called to obtain
  ! information such as
  ! * total number of molecules of a species in a given box
  ! * Actual index of imth molecule in the box
  ! * Position of a given molecule in global molecule_list array
  !
  ! Contains the following routines
  !
  ! Get_Nmolecules_Species
  ! Get_Index_Molecule
  ! Get_Position_Molecule
  !
  ! USED by
  !
  !   accumulate
  !   angle_distortion
  !   bond_stretching_move
  !   deletion
  !   read_write_checkpoint
  !   rigid_dihedral_change
  !   rotate
  !   translate
  !   write_properties
  !
  ! Revision History
  !
  !   08/01/13 : Beta test version
  !
  !********************************************************************************

  USE Type_Definitions
  USE Global_Variables
  USE Energy_Routines

  IMPLICIT NONE

CONTAINS

 !****************************************************************************
  SUBROUTINE Get_Index_Integer_Molecule(this_box,this_species,im,alive)

    
    IMPLICIT NONE
    
    INTEGER :: this_box, this_species, nmolecules_species

    INTEGER :: i, alive, im
    
    
    nmolecules_species = 0

    DO i = 1, nmols(this_species,this_box)
       
       alive = locate(i,this_species,this_box)
       
       IF ( molecule_list(alive,this_species)%live ) THEN
          IF (molecule_list(alive,this_species)%molecule_type == int_normal) THEN

             nmolecules_species = nmolecules_species + 1
             
          END IF
       END IF
       
       IF (nmolecules_species == im ) EXIT
       
    END DO

  END SUBROUTINE Get_Index_Integer_Molecule

  SUBROUTINE Get_Position_Locate(is,alive,position)

    IMPLICIT NONE

    INTEGER :: is, alive, position, i, this_locate, this_box

    this_box = molecule_list(alive,is)%which_box
    DO i = 1, nmols(is,this_box)

       this_locate = locate(i,is,this_box)

       IF (this_locate == alive ) EXIT

    END DO

    position = i

  END SUBROUTINE Get_Position_Locate

!  !****************************************************************************

  SUBROUTINE Compute_Beads(this_box)

    ! This subroutine computes total number of beads in a given box
    ! 
    ! Written by Jindal Shah on 03/24/09
    !
    !

    INTEGER, INTENT(IN) :: this_box

    INTEGER :: is, im, alive, ia, ia_type, ia_rxnum


    
    nint_beads(:,this_box) = 0
    DO is = 1, nspecies
       DO im = 1, nmols(is,this_box) 

          alive = locate(im,is,this_box)
          
          IF (.NOT. molecule_list(alive,is)%live ) CYCLE

          IF (molecule_list(alive,is)%molecule_type == int_normal) THEN
             
             DO ia = 1, natoms(is)
                
                ia_type = nonbond_list(ia,is)%atom_type_number
                IF (ia_type /= 0 ) nint_beads(ia_type,this_box) = &
                                   nint_beads(ia_type,this_box) + 1
             END DO
            
          END IF
 
       END DO
       
    END DO
 

    IF(this_box == 1) THEN

       IF ( .NOT. ALLOCATED(nbeads_in)) ALLOCATE(nbeads_in(nbr_atomtypes))
       IF ( .NOT. ALLOCATED(nbeads_out)) ALLOCATE(nbeads_out(nbr_atomtypes))
       IF ( .NOT. ALLOCATED(nbeadsfrac_in)) &
           ALLOCATE(nbeadsfrac_in(nbr_atomtypes))

    END IF
  END SUBROUTINE Compute_Beads

!  !****************************************************************************

  SUBROUTINE Compute_Pressure(this_box)
    !***************************************************************************
    !
    ! This subroutine calculates the pressure of the the box
    ! 
    ! CALLED BY: 
    !       Write_Properties
    !
    ! CALLS :
    !       Compute_System_Total_Force
    !
    ! Written by Ryan Mullen on 06/11/16
    !
    !***************************************************************************

    IMPLICIT NONE

    INTEGER :: this_box

    ! start with the ideal gas pressure
    pressure(this_box)%computed = SUM(nmols(:,this_box)) * temperature(this_box) &
                                * kboltz / box_list(this_box)%volume

    ! add the pressure from the virial
    CALL Compute_System_Total_Force(this_box)

    pressure_tensor(:,:,this_box) = W_tensor_total(:,:,this_box) &
                                  / box_list(this_box)%volume
    pressure(this_box)%computed = pressure(this_box)%computed &
                                + ((pressure_tensor(1,1,this_box) &
                                  + pressure_tensor(2,2,this_box) &
                                  + pressure_tensor(3,3,this_box)) / 3.0_DP)
    
    ! add pressure from tail corrections
    IF(int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
       pressure(this_box)%computed = pressure(this_box)%computed &
                                   + virial(this_box)%lrc &
                                   / box_list(this_box)%volume
    END IF

  END SUBROUTINE Compute_Pressure

!  !****************************************************************************

END MODULE Simulation_Properties


