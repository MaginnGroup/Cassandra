!*******************************************************************************/A
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
! This file contains three subroutines: 
!
!
! Fragment_Driver
!
!    It performs MC moves to sample angle distribution for a given fragment
!
! Change_Theta_Phi
!
!    Changes polar and azimuthal angles of a given atom in the fragment
!
! 08/07/13 : Created beta version
!    
!*******************************************************************************

SUBROUTINE Fragment_Driver
  !*****************************************************************************
  ! The subroutine carries out atoms displacement to sample conformations of a 
  ! branch point. 
  !
  ! CALLED BY
  !
  !        main
  !
  ! CALLS
  ! 
  !        Compute_Molecule_Angle_Energy
  !        Compute_Molecule_Improper_Energy
  !        Save_Old_Internal_Coordinates
  !        Save_Old_Cartesian_Coordinates
  !        Get_Interal_Coordinates
  !        Change_Phi_Theta
  !        Revert_Old_Cartesian_Coordinates
  !
  !***************************************************************************************

  USE Global_Variables
  USE Random_Generators
  USE Energy_Routines
  USE File_Names

  IMPLICIT NONE

  INTEGER :: is, im, ibox, ia

  ! Fragment routine is only called with 1 box, species, and molecule
  is = 1
  im = 1
  ibox = 1

  OPEN(UNIT=frag_file_unit,file=frag_file(is))
  WRITE(frag_file_unit,*) (n_mcsteps-n_equilsteps)/nthermo_freq

  DO i_mcstep = 1, n_mcsteps
     
     CALL Atom_Displacement

     ! Store information with given frequency
     IF ( i_mcstep > n_equilsteps ) THEN

        IF (MOD(i_mcstep,nthermo_freq) == 0) THEN
           
           WRITE(frag_file_unit,*) temperature(ibox), energy(ibox)%total
           DO ia = 1, natoms(is)
              WRITE(frag_file_unit,*) nonbond_list(ia,is)%element, &
                                      atom_list(ia,im,is)%rxp, &
                                      atom_list(ia,im,is)%ryp, &
                                      atom_list(ia,im,is)%rzp
           END DO
           
        END IF
        
     END IF

     ! revert coordinates, revert        
    
  END DO
  WRITE(*,'(A30,I10)') 'Number of Trial moves', ntrials(is,ibox)%disp_atom
  WRITE(*,'(A30,I10)') 'Accepted moves',  nsuccess(is,ibox)%disp_atom

END SUBROUTINE Fragment_Driver

