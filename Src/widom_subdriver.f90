!*******************************************************************************
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

SUBROUTINE Widom_Subdriver

  !*****************************************************************************
  ! 
  ! PURPOSE: 
  !
  ! Called by
  !
  !    
  !
  ! Revision history
  !
  !   
  !   
  !   
  ! DESCRIPTION: This subroutine performs the following steps:
  !
  ! 
  !*****************************************************************************

  USE Global_Variables

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  IMPLICIT NONE

  ! Arguments

  ! Local declarations
  INTEGER :: is, ibox

  
  REAL(DP) :: widom_sum, widom_avg




  ! Loop over all species

  DO is = 1, nspecies
        ! Loop over all boxes
        DO ibox = 1, nbr_boxes
                ! move on to next box if this species is not used as a test particle in this box
                IF (.NOT. species_list(is)%test_particle(ibox)) CYCLE
                ! move on to next box if the step number isn't divisible by the widom_interval for this species.
                ! these are separate IF statements because we don't want to divide by zero
                IF (MOD(i_mcstep,species_list(is)%widom_interval(ibox)) .NE. 0) CYCLE
                CALL Widom_Insert(is,ibox,widom_sum)
                species_list(is)%widom_sum(ibox) = species_list(is)%widom_sum(ibox) + widom_sum


                widom_avg = widom_sum / species_list(is)%insertions_in_step(ibox)
                CALL Write_Widom_Properties(is,ibox,widom_avg)



        END DO




  END DO












END SUBROUTINE Widom_Subdriver
!*******************************************************************************
