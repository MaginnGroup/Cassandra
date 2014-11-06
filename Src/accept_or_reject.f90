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

FUNCTION accept_or_reject(factor)

  !*******************************************************************
  !
  ! The function will test, based on factor, if a move
  ! is accepted.
  !
  ! factor is so calculated that the acceptance is in the
  ! form min(1,exp(-factor)) 
  !
  ! 12/10/13  : Beta Release
  !*********************************************************************

  USE Type_Definitions, ONLY : DP
  USE Run_Variables, ONLY : max_kBT
  USE Random_Generators

  IMPLICIT NONE

  REAL(DP), INTENT(IN) :: factor

  REAL(DP) :: p_acc

  LOGICAL accept_or_reject

  accept_or_reject = .FALSE.

  IF (factor <= 0.0_DP) THEN

     accept_or_reject = .TRUE.

  ELSE IF ( factor < max_kBT) THEN

     p_acc = MIN(1.0_DP, DEXP(-factor))

     IF ( rranf() <= p_acc ) THEN
        
        accept_or_reject = .TRUE.

     END IF

  END IF

END FUNCTION accept_or_reject
