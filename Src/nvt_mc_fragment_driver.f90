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
! This file contains three subroutines: 
!
!
! NVT_MC_Fragment_Driver
!
!    It performs MC moves to sample angle distribution for a given fragment
!
! Change_Theta_Phi
!
!    Changes polar and azimuthal angles of a given atom in the fragment
!
! 08/07/13 : Created beta version
!    
!*********************************************************************************************

SUBROUTINE NVT_MC_Fragment_Driver
  !************************************************************************************
  ! The subroutine carries out atoms displacement to sample conformations of a branch
  ! point. 
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

  USE Run_Variables
  USE Random_Generators
  USE Energy_Routines
  USE File_Names

  IMPLICIT NONE

  INTEGER :: is, im, this_box, i, rand_atom, naccept, ia
  INTEGER :: naverage

  REAL(DP) :: e_angle_old, e_angle_new, delta_e_angle, p_acc, old_coord, new_coord
  REAL(DP) :: area_o, area_n
  REAL(DP) :: e_improper_n, e_improper_o, delta_e_improper, delta_e, e_total_o
  REAL(DP) :: ac_frag_energy, zig_over_omega

  LOGICAL :: theta_bound

  ! compute angle energy for the input conformation

  DO is = 1, nspecies
     
     naccept = 0

     im = locate(1,is)

     CALL Compute_Molecule_Angle_Energy(im,is,e_angle_old)
     e_total_o = e_angle_old

     IF(nimpropers(is) .GT. 0) THEN
        CALL Compute_Molecule_Improper_Energy(im,is,e_improper_o)
        e_total_o = e_total_o + e_improper_o
     END IF

     this_box = molecule_list(im,is)%which_box

     OPEN(UNIT=frag_file_unit,file=frag_file(is))
     WRITE(frag_file_unit,*) (n_mcsteps-n_equilsteps)/nthermo_freq

     ac_frag_energy = 0.0_DP
     zig_over_omega = 0.0_DP

     DO i = 1, n_mcsteps
        
        ! Note that the 1st atom is at the origin while the second atom is
        ! aligned along the x axis. So, we sample only the natoms(is) - 2 atoms

        
        rand_atom = INT ( (natoms(is) - 2) * rranf()) + 3

        ! save the coordinates

        CALL Save_Old_Cartesian_Coordinates(im,is)
        CALL Save_Old_Internal_Coordinates(im,is)

        CALL Change_Phi_Theta(rand_atom,im,is,theta_bound)
        
        CALL Compute_Molecule_Angle_Energy(im,is,e_angle_new)
        IF(nimpropers(is) .GT. 0) CALL Compute_Molecule_Improper_Energy(im,is,e_improper_n)

        delta_e_angle = e_angle_new - e_angle_old 
        IF(nimpropers(is) .GT. 0) THEN
           delta_e_improper = e_improper_n - e_improper_o
        ELSE
           delta_e_improper = 0.0_DP
        END IF
        delta_e = delta_e_angle + delta_e_improper

        p_acc = min(1.0_DP, DEXP(-beta(this_box)*delta_e))
        
        IF ( rranf() < p_acc) THEN

           ! update energies 

           e_angle_old = e_angle_old + delta_e_angle
           e_total_o = e_angle_old

           IF(nimpropers(is) .GT. 0) THEN
              e_improper_o = e_improper_o + delta_e_improper
              e_total_o = e_total_o + e_improper_o
           END IF

           naccept = naccept + 1

           CALL Get_Internal_Coordinates(im,is)
           
        ELSE

           CALL Revert_Old_Cartesian_Coordinates(im,is)

        END IF

        ! Store information with given frequency

        IF ( i > n_equilsteps ) THEN

           ac_frag_energy = ac_frag_energy + e_total_o
           zig_over_omega = zig_over_omega + EXP(beta(this_box) * e_total_o)
          
           IF (MOD(i,nthermo_freq) == 0) THEN
              !           WRITE(frag_file_unit,*) natoms(is)
              
              WRITE(frag_file_unit,*) temperature(this_box), e_angle_old
              DO ia = 1, natoms(is)
                 WRITE(frag_file_unit,*) nonbond_list(ia,is)%element, atom_list(ia,im,is)%rxp, atom_list(ia,im,is)%ryp, &
                      atom_list(ia,im,is)%rzp
              END DO
              
           END IF
           
        END IF

        ! revert coordinates, revert        
       
     END DO
     write(*,'(A30,I10)') 'Number of Trial moves', n_mcsteps
     write(*,'(A30,I10)') 'Accepted moves',  naccept
     write(*,*) 'Average intramolecular energy', ac_frag_energy / REAL(n_mcsteps,DP) * atomic_to_kJmol, 'kJ/mol'
     WRITE(*,'(A,2X,F30.10)') 'Z/Omega ',1.0_DP/ (zig_over_omega / REAL(n_mcsteps - n_equilsteps, DP))
     
  END DO



END SUBROUTINE NVT_MC_Fragment_Driver

!**********************************************************************************
SUBROUTINE Change_Phi_Theta(this_atom,im,is,theta_bound)
!**********************************************************************************
  USE Run_Variables
  USE Random_Generators, ONLY : rranf

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_atom,im,is

  LOGICAL, INTENT(INOUT) :: theta_bound

  REAL(DP) :: this_x, this_y, this_z, rho, bond_length, theta, phi, dcostheta
  REAL(DP) :: dphi

  theta_bound = .false.


!  delta_cos = 0.2_DP
!  delta_phi = 10.0_DP * PI/180.0_DP

  ! Get spherical coordinates

  this_x = atom_list(this_atom,im,is)%rxp
  this_y = atom_list(this_atom,im,is)%ryp
  this_z = atom_list(this_atom,im,is)%rzp

  rho = this_x * this_x + this_y * this_y 
  bond_length = this_z * this_z + rho

  rho = DSQRT(rho)
  bond_length = DSQRT(bond_length)

  theta = DACOS(this_z/bond_length)
  
  phi = DASIN(this_y/rho)


  IF ( this_x < 0.0_DP ) THEN
     
     phi = PI - DASIN(this_y/rho)

  END IF

  ! Change theta and phi

  dcostheta = (2.0_DP * rranf() - 1.0_DP ) * delta_cos_max 
  dphi = (2.0_DP * rranf() - 1.0_DP ) * delta_phi_max 


!  theta = DACOS(2.0_DP * rranf() - 1.0_DP)
!  phi = (2.0_DP * rranf()) * PI

  ! new polar and azimuthal anlges
  IF ( (ABS(DCOS(theta) + dcostheta) > 1.0_DP) ) THEN
     theta_bound = .true.
     RETURN
  END IF

  theta =  DACOS(DCOS(theta) + dcostheta)
  phi = phi + dphi
  
  ! new coordinates
  atom_list(this_atom,im,is)%rxp = bond_length * DSIN(theta) * DCOS(phi)
  atom_list(this_atom,im,is)%ryp = bond_length * DSIN(theta) * DSIN(phi)
  atom_list(this_atom,im,is)%rzp = bond_length * DCOS(theta)

END SUBROUTINE Change_Phi_Theta
 
     


  



 
