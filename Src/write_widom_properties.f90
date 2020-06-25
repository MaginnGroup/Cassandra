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
SUBROUTINE Write_Widom_Properties(is,this_box,widom_avg)
  !
  ! CALLED BY
  !
  !     Widom_Subdriver
  !
  ! CALLS
  !
  !   None
  !
!*******************************************************************************

  USE Global_Variables
  USE File_Names

  IMPLICIT NONE

  INTEGER :: this_box, this_unit, is

  REAL(DP) :: widom_avg


     
  !-- check to see whether the file is open or not
  this_unit = wprop_file_unit(is,this_box)
     
  IF (first_open_wprop(is,this_box)) THEN
        
     OPEN(unit=this_unit,file=wprop_filenames(is,this_box))
        
     ! write the header information that indicates the properties contained
     ! in the file
     !CALL Write_Widom_Header(i)
     WRITE(this_unit,'(A12,10X,A25)') 'Step_#', 'Average widom_var in step'



     first_open_wprop(is,this_box) = .FALSE.

  END IF
        
  !CALL Write_Properties_Buffer(i)
  !widom_avg = widom_sum / species_list(is)%insertions_in_step(this_box)
  WRITE(this_unit,'(I12,10X,E25.17)') i_mcstep, widom_avg
     


 
END SUBROUTINE Write_Widom_Properties
