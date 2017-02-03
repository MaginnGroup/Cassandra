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

SUBROUTINE Clean_Abort(fail_message,sub_where_failed)

  !***********************************************************************
  ! This subroutine is called when a fatal error occurs in the program
  ! and is meant to cleanly shut down the program.  Two things are 
  ! passed in:
  !       a 10-dimension array of strings containing error messages
  !       the name of the subroutine where the failure occured
  ! The error message are printed to the standard output and the
  ! logfile.  Any open files are closed down.
  !
  ! Original routine created by: David Eike on February 19, 2005
  !
  !
  ! Called by
  !
  !   angle_dist_pick
  !   compute_cell_dimensions
  !   create_intra_exclusion_table
  !   create_nonbondt_table
  !   embedded_atom
  !   energy_routines
  !   file_names
  !   fragment_growth
  !   gcmc_control
  !   gemc_control
  !   gemc_nvt_control
  !   input_routines
  !   io_utilites
  !   main
  !   minimum_image_separtion
  !   nvt_mc_ring_fragment
  !   participation
  !   rigid_dihedral_change
  !   rigid_insertion
  !   volume_change
  !
  ! Revision history
  !
  !   12/10/13 : Beta Release
  !***********************************************************************
  USE File_Names
  USE IO_Utilities

  IMPLICIT NONE

  CHARACTER(*) :: sub_where_failed
  CHARACTER(*), DIMENSION(10) :: fail_message
  INTEGER :: ii
  LOGICAL :: FileOpen, FileExist

!***********************************************************************

  !*****************************************************
  !                                                    *
  ! Write error message to standard output       
  !                                                    *
  !*****************************************************


  DO ii = 1, 10
     IF (fail_message(ii) .NE. "") THEN
        WRITE(*,'(A)') TRIM(fail_message(ii))
     END IF
  END DO

  WRITE(*,'(5A)') 'This error occurred in subroutine ',&
       TRIM(sub_where_failed), ' on step ', &
       TRIM(Int_To_String(i_mcstep)), '.'
  WRITE(*,'(A)') 'Fatal Error. Stopping program.'
  
  
  !*****************************************************
  !                                                    *
  ! Write error message to logfile if opened 
  !                                                    *
  !*****************************************************
  INQUIRE(logunit,OPENED=FileExist)
  INQUIRE(logunit,OPENED=FileOpen)
  IF (FileExist .and. FileOpen) THEN
     WRITE(logunit,*)
     WRITE(logunit,'(A)') '************ERROR************'
     DO ii = 1, 10
        IF (fail_message(ii) .NE. "") THEN
           WRITE(logunit,'(A)') TRIM(fail_message(ii))
        END IF
     END DO
     WRITE(logunit,'(5A)') 'This error occurred in subroutine ',&
          TRIM(sub_where_failed), ' on step ', &
          TRIM(Int_To_String(i_mcstep)), '.'
     WRITE(logunit,'(A)') 'Fatal Error. Stopping program.'
     CLOSE(logunit)
  END IF
  

STOP
  
END SUBROUTINE Clean_Abort
