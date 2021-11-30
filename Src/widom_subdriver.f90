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
  ! PURPOSE: This subroutine checks which species and box combinations should
  !          have Widom insertions in the current step, calls Widom_Insert to 
  !          execute Widom insertions for those species and box combinations,
  !          and calls Write_Widom_Properties to write the results to Widom 
  !          property files.
  !
  ! Called by
  !
  !   nvtmc_driver
  !   nptmc_driver
  !   gcmc_driver
  !   gemc_driver 
  !
  !*****************************************************************************

  USE Global_Variables
  USE Sector_Routines

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  IMPLICIT NONE


  ! Local declarations
  INTEGER :: is, ibox
  REAL(DP) :: widom_sum, widom_avg
  LOGICAL :: need_init
  need_init = .TRUE.
  IF (.NOT. l_sectors) need_init = .FALSE.
  ! Loop over all species
  DO is = 1, nspecies
        ! Loop over all boxes
        DO ibox = 1, nbr_boxes
                ! move on to next box if this species is not used as a test particle in this box
                IF (.NOT. species_list(is)%test_particle(ibox)) CYCLE
                ! move on to next box if the step number isn't divisible by the widom_interval for this species.
                ! these are separate IF statements because we don't want to divide by zero
                IF (MOD(i_mcstep,species_list(is)%widom_interval(ibox)) .NE. 0) CYCLE
                IF (need_init) THEN
                        CALL Sector_Setup
                        need_init = .FALSE.
                END IF
                CALL Widom_Insert(is,ibox,widom_sum)
                species_list(is)%widom_sum(ibox) = species_list(is)%widom_sum(ibox) + widom_sum
                widom_avg = widom_sum / species_list(is)%insertions_in_step(ibox)
                CALL Write_Widom_Properties(is,ibox,widom_avg)
        END DO
  END DO
END SUBROUTINE Widom_Subdriver
!*******************************************************************************
