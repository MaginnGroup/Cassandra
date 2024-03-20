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
SUBROUTINE Write_Widom_Properties(is,this_box,widom_avg,t_cpu,n_overlaps)
  !
  ! CALLED BY
  !
  !     Widom_Subdriver
  !
  !
!*******************************************************************************

  USE Global_Variables
  USE File_Names

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box, is
  INTEGER(KIND=INT64), INTENT(IN) :: n_overlaps

  REAL(DP), INTENT(IN) :: t_cpu
  REAL(DP) :: widom_avg

  INTEGER :: this_unit

  LOGICAL :: is_sweeps

  is_sweeps = (sim_length_units == 'Sweeps')
     
  !-- check to see whether the file is open or not
  this_unit = wprop_file_unit(is,this_box)
     
  IF (first_open_wprop(is,this_box)) THEN
     OPEN(unit=this_unit,file=wprop_filenames(is,this_box))
     ! write the header information that indicates the properties contained
     ! in the file
     !CALL Write_Widom_Header(i)
     IF (is_sweeps) THEN
             WRITE(this_unit,'(A19,7X,A30,7X,A30,7X,A17)') 'Sweep_#', 'Average widom_var for sweep', 'Widom CPU time for sweep (s)', &
                     'Number of overlaps'
     ELSE
             WRITE(this_unit,'(A19,7X,A30,7X,A30,7X,A17)') 'Step_#', 'Average widom_var for step', 'Widom CPU time for step (s)', &
                     'Number of overlaps'
     END IF
     first_open_wprop(is,this_box) = .FALSE.
  END IF
  !CALL Write_Properties_Buffer(i)
  !widom_avg = widom_sum / species_list(is)%insertions_in_step(this_box)
  IF (widom_avg < 1.0e-99_DP) widom_avg = 0.0_DP
  IF (is_sweeps) THEN
          WRITE(this_unit,'(I19,7X,E30.22,7X,E30.22,7X,I19)') i_mcstep/steps_per_sweep, widom_avg, t_cpu, n_overlaps
  ELSE
          WRITE(this_unit,'(I19,7X,E30.22,7X,E30.22,7X,I19)') i_mcstep, widom_avg, t_cpu, n_overlaps
  END IF
END SUBROUTINE Write_Widom_Properties
