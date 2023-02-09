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

MODULE XTC_Routines
  !************************************************************************
  !
  !**************************************************************************
  USE Global_Variables
  USE File_Names
  USE Simulation_Properties
  USE IO_Utilities

  IMPLICIT NONE

CONTAINS

  SUBROUTINE Open_XTC(ibox)
          INTEGER :: ibox

          err_msg = ''
          err_msg(1) = 'Cassandra must be compiled with libgmxfort to support xtc file operations.'
          CALL Clean_Abort(err_msg, 'Open_XTC')

  END SUBROUTINE Open_XTC

  SUBROUTINE Close_XTC(ibox)
          INTEGER :: ibox

          err_msg = ''
          err_msg(1) = 'Cassandra must be compiled with libgmxfort to support xtc file operations.'
          CALL Clean_Abort(err_msg, 'Close_XTC')

  END SUBROUTINE Close_XTC

  INTEGER FUNCTION Skip_XTC_Frames(ibox,nframes)
          INTEGER :: ibox, nframes

          err_msg = ''
          err_msg(1) = 'Cassandra must be compiled with libgmxfort to support xtc file operations.'
          CALL Clean_Abort(err_msg, 'Skip_XTC_Frames')
          Skip_XTC_Frames = 0

  END FUNCTION Skip_XTC_Frames

  LOGICAL FUNCTION Read_XTC_Frame(ibox)
          INTEGER :: ibox

          err_msg = ''
          err_msg(1) = 'Cassandra must be compiled with libgmxfort to support xtc file operations.'
          CALL Clean_Abort(err_msg, 'Read_XTC_Frame')
          Read_XTC_Frame = .TRUE.

  END FUNCTION Read_XTC_Frame

  FUNCTION Get_XTC_Coords(ibox)
          INTEGER :: ibox
          REAL(DP), DIMENSION(natoms_to_read(ibox),3) :: Get_XTC_Coords

          err_msg = ''
          err_msg(1) = 'Cassandra must be compiled with libgmxfort to support xtc file operations.'
          CALL Clean_Abort(err_msg, 'Get_XTC_Coords')

          Get_XTC_Coords = 0.0_DP

  END FUNCTION Get_XTC_Coords

  FUNCTION Get_XTC_Box(ibox)
          INTEGER :: ibox
          REAL(DP), DIMENSION(3,3) :: Get_XTC_Box

          err_msg = ''
          err_msg(1) = 'Cassandra must be compiled with libgmxfort to support xtc file operations.'
          CALL Clean_Abort(err_msg, 'Get_XTC_Box')
          Get_XTC_Box = 0.0_DP

  END FUNCTION Get_XTC_Box




!*******************************************************************************

END MODULE XTC_Routines

