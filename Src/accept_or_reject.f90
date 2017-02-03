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

FUNCTION accept_or_reject(ln_pacc)

  !*******************************************************************
  !
  ! The function will test, based on ln_pacc, if a move
  ! is accepted.
  !
  ! ln_pacc is so calculated that the acceptance is in the
  ! form exp[-max(0,ln_pacc)]
  !
  ! 12/10/13  : Beta Release
  !*********************************************************************

  USE Type_Definitions, ONLY : DP
  USE Global_Variables, ONLY : max_kBT
  USE Random_Generators

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: ln_pacc

  REAL(DP) :: pacc

  LOGICAL accept_or_reject

  accept_or_reject = .FALSE.

  IF (ln_pacc <= 0.0_DP) THEN

     accept_or_reject = .TRUE.

  ELSE IF ( ln_pacc < max_kBT) THEN

     pacc = DEXP(-ln_pacc)

     IF ( rranf() <= pacc ) THEN
        
        accept_or_reject = .TRUE.

     END IF

  END IF

END FUNCTION accept_or_reject
