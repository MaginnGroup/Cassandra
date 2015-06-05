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

MODULE Energy_Routines
  !-----------------------------------------------------------------------------------
  ! This modules contains a collection of all the routines involved in computing
  ! energies and associated quantities. 
  ! 
  ! Compute_Bond_Energy: passed two atom indices, the molecule index and the species
  ! index
  !                      it returns the energy of the bond between those two atoms, 
  !                      consistent with the potential type. 
  !                      Currently supports fixed and harmonic.
  !
  ! Compute_Molecule_Bond_Energy: passed a molecule and species index, this returns
  ! the total bond energy associated with that molecule. 
  !                       Currently supports fixed and harmonic.
  ! 
  ! Compute_Angle_Energy: passed indices for three bonded atoms atom1,atom2,atom3 where atom3 
  !                       is the apex of a the angle, along with the molecule and species index
  !                       it returns the energy associated with that bond angle.
  !                       Currently supports fixed and harmonic.
  !
  ! Compute_Molecule_Angle_Energy: Passed a molecule and species index, it returns the total energy
  !                       of that molecule due to bond angles.
  !                       Currently supports fixed and harmonic.
  !
  ! Compute_Dihedral_Energy: Passed 4 contiguous atoms that form a dihedral and the molecule and
  !                        species indices, it returns the energy of that dihedral.  
  !                        Currently supports none and OPLS.
  !
  ! Compute_Molecule_Dihedral_Energy: Passed a molecule and species index, it returns the total 
  !                        dihedral energy of that molecule. 
  !                        Currently supports none and OPLS.
  !
  ! Compute_Improper_Energy: Passed the indices of 4 atoms making up an improper angle, along with
  !                        the species and molecule index, this returns the energy of that improper. 
  !                        This has not been tested!!! 
  !                        Supports improper of the harmonic form: E = k*(phi-phi0)^2
  !
  ! Compute_Molecule_Imroper_Energy: Passed molecule and species indices, this computes the 
  !                         total improper energy of the molecule. Not yet tested!!!
  !                        Supports improper of the harmonic form: E = k*(phi-phi0)^2
  !
  ! Compute_Atom_Nonbond_Energy: passed indices of an atom, molecule and species, this 
  !                        returns the vdw and either direct charge-charge or the real space
  !                        part of the Ewald energy of this atom with all existing atoms in the 
  !                        system. It accounts for intramolecular scaling of 1-2, 1-3 and 1-4.
  !                        Supports vdw_style = none or LJ and charge_style none or coul.
  !                        For LJ, it supports rcut, cut_tail and cut_shift, though TAIL CORRECTIONS
  !                        HAVE NOT YET BEEN ADDED. LJ is assumed to be 12-6 LJ.
  !                        For charge_style = coul, it supports rcut and Ewald is roughed in.
  !                        However, the Ewald parts of the code need some thought, especially in 
  !                        light of computing energy differences. This routine also returns the 
  !                        virial contribution. It needs more testing, but I believe it works. 
  !
  ! Pair_Energy: Computes the vdw and q-q pair energy between atoms ia and ja of molecules im and jm
  !              and species is and js, given their separation rijsq. I have passed each component of 
  !              separation nut right now this is unnecessary. 
  !              It also computes the real space part of the Ewald sum if necessary.
  !
  !              LJ potential: Eij = 4*epsilon(i,j) * [ (sigma(i,j)/rij)^12 - (sigma(i,j)/rij)^6 ]
  !              Wij = -rij/3 * d Eij / d rij. Use the virial in: P = NkBT + < W >
  !
  ! Ewald_Real: Real space part of Ewald sum. Not tested. Need to add reciprical, self
  !             and energy difference sin and cos sums. Contains erfc function.
  !
  ! 
  ! Ewald_Reciprocal_Lattice_Vector_Setup : Sets up lattice vectors for Ewald
  !         Summation for the input box.
  !
  ! Compute_System_Ewald_Reciprocal_Energy : Computes reciprocal space energy
  !          for a given box. 
  !
  ! Compute_Ewald_Reciprocal_Energy_Difference: Computes difference in 
  !          reciprocal space energy due to various moves. The routine 
  !          makes use of the fact that for a given move, the coordinates
  !          of only one moleule are perturbed. Hence cos_sum and sin_sum
  !          arrays can be computed by taking differences of q_i cos(k * r_i)
  !          terms in new and old configurations.
  !
  !  
  ! Compute_System_Ewald_Self_Energy: Calculation of self energy for the Ewald
  !          summation is obtained from this subroutine.
  !
  ! Compute_System_Ewald_Self_Energy_Difference: Determines the difference
  !          in self energy when the number of particles changes in the system.
  !
  ! Compute_Total_System_Energy: Computes the total system energy of the system.
  ! Compute_Molecule_Pair_Interaction: Comptues the energy of interactions
  !          between pair of molecules in a given box.
  !
  ! Compute_LR_Correction: Determines long range correction when the flag
  !          is set to 'cut_tail'.
  !
  ! Used by
  !
  !   angle_distortion
  !   atom_displacement
  !   chempot
  !   cutNgrow
  !   deletion
  !   fragment_growth
  !   gcmc_control
  !   gcmc_driver
  !   gemc_control
  !   gemc_driver
  !   gemc_nvt_volume
  !   gemc_particle_transfer
  !   grow_molecules
  !   input_routines
  !   insertion
  !   main
  !   nptmc_control
  !   nptmc_driver
  !   nvtmc_control
  !   nvtmc_driver
  !   nvt_mc_fragment_driver
  !   nvt_mc_ring_fragment
  !   precalculate
  !   rotate
  !   translate
  !   volume_change
  !   write_properties
  !   zig_by_omega
  !
  ! Revision history
  !
  !   12/10/13  : Beta Release
  !-----------------------------------------------------------------------------------

  USE Type_Definitions
  USE Run_Variables
  USE File_Names
  USE Pair_Nrg_Routines
 !$  USE OMP_LIB

  IMPLICIT NONE

CONTAINS
  
  !----------------------------------------------------------------------------------------------

  SUBROUTINE Compute_Bond_Energy(at1,at2,molecule,species,energy)
    !----------------------------------------------------------------------------------------------             
    ! Given two atoms, the molecule and the species type, this routine figures out the bond 
    ! parameters and returns the energy associated with the bond. It uses the participation matrices
    ! computed earlier to avoid having to loop over all bonds in the molecule.

    ! CALLED By:
    ! CALLS: Get_Bond_Length

    ! Written by: E. Maginn
    ! Date: Mon Nov 26 06:27:24 MST 2007
    ! Revision history:    

    ! Passed
    INTEGER :: at1,at2,molecule,species

    !Returned
    REAL(DP) :: energy

    ! local
    INTEGER :: i, nbr_bonds,ierr,ibond
    REAL(DP) :: k,l0,length
  !----------------------------------------------------------------------------------------------
    ierr = 1

    ! Compute the number of atoms bonded to at1
    nbr_bonds = bondpart_list(at1,species)%nbonds

    ! Loop over these local bonds and find a match and the global bond number
    DO i=1,nbr_bonds

       IF  (at2 == bondpart_list(at1,species)%atom(i)) THEN
          ! The global bond number is:
          ibond = bondpart_list(at1,species)%bond_num(i)
         
          ! Use the global bond number to determine the type of potential and compute the energy
          IF ( bond_list(ibond,species)%bond_potential_type == 'fixed') THEN
             energy = 0.0_DP
             ierr = 0
          ELSEIF (bond_list(ibond,species)%bond_potential_type == 'harmonic') THEN
             k=bond_list(ibond,species)%bond_param(1)
             l0 = bond_list(ibond,species)%bond_param(2)
                
             CALL Get_Bond_Length(ibond,molecule,species,length)
             energy = k*(length-l0)**2
             ierr = 0

             ! Add more potential functions here.
          ENDIF

          ! Exit once a match is found
          EXIT

       ENDIF

    ENDDO

    ! Test to make sure we successfuly found the right energy
    IF (ierr .NE. 0) THEN
       err_msg = ""
       err_msg(1) = "Unable to compute bond energy"
       CALL Clean_Abort(err_msg,'Compute_Bond_Energy')
    END IF


  END SUBROUTINE Compute_Bond_Energy
  !----------------------------------------------------------------------------------------------            


  SUBROUTINE Compute_Molecule_Bond_Energy(molecule,species,energy)
    !----------------------------------------------------------------------------------------------            
    ! Given a molecule number and species,this routine computes the total bond energy
    ! of the entire molecule.

    ! CALLED By:
    ! CALLS: Get_Bond_Length

    ! Written by: E. Maginn
    ! Date: Mon Nov 26 06:27:24 MST 2007
    ! Revision history:    

    ! Passed to
    INTEGER :: molecule,species

    ! Passed from
    INTEGER :: ibond
    REAL(DP) :: length

    ! Returned
    REAL(DP) :: energy

    ! Local
    REAL(DP) :: k,l0,eb
  !----------------------------------------------------------------------------------------------            
    energy = 0.0_DP
    DO ibond=1,nbonds(species)
       IF (bond_list(ibond,species)%int_bond_type == int_none) THEN
          eb = 0.0_DP
       ELSEIF (bond_list(ibond,species)%int_bond_type == int_harmonic) THEN
          k=bond_list(ibond,species)%bond_param(1)
          l0 = bond_list(ibond,species)%bond_param(2)
          CALL Get_Bond_Length(ibond,molecule,species,length)
          eb = k*(length-l0)**2

          ! Add more potential functions here.
       ENDIF
       energy = energy + eb
    ENDDO
  

  END SUBROUTINE Compute_Molecule_Bond_Energy
  !----------------------------------------------------------------------------------------------            


  SUBROUTINE Compute_Angle_Energy(at1,at2,at3,molecule,species,energy)
  !----------------------------------------------------------------------------------------------            
    ! This routine is passed three atoms with at2 being the central atom, along with the molecule and 
    ! species number, and computes the energy associated with the angle formed by the three atoms.
    ! It does this by finding the number of angles associated with the central atom (at2), finds the global
    ! bond angle that has the proper matches for at1 and at3, determines the right parameters to use,
    ! calls for the angle if necessary, then returns the energy. The use of the participation 
    ! matrix avoids having to loop over all angles in the molecule, and instead only requires
    ! that a very small number of angles be searched for a match. 

    ! Called by:
    ! Calls: Get_Bond_Angle
    
    ! Written by: E. Maginn
    ! Date: Mon Nov 26 06:27:24 MST 2007
    ! Revision history:    
  !----------------------------------------------------------------------------------------------            
    ! Passed
    INTEGER :: at1,at2,at3,molecule,species

    ! Returned
    REAL(DP) :: energy

    ! Local
    INTEGER :: i,ierr,nbr_angles,iangle
    REAL(DP) :: k,theta0,theta
  !----------------------------------------------------------------------------------------------            
    ierr = 1

    ! Number of angles associated with the central atom
    nbr_angles = angle_part_list(at2,species)%nangles

    DO i=1,nbr_angles
       IF (angle_part_list(at2,species)%position(i) == 2) THEN

          ! Now we know that global angle angle_part_list(at2,species)%which_angle(i)
          ! may be the angle we want. Test for this by going back to the original 
          ! global angle array and looking for matches of at1 and at3. 

          iangle = angle_part_list(at2,species)%which_angle(i)
          IF ( (angle_list(iangle,species)%atom1 == at1) .AND. &
               (angle_list(iangle,species)%atom3 == at3) ) THEN

             IF (angle_list(iangle,species)%angle_potential_type == 'fixed') THEN
                energy = 0.0_DP
                ierr = 0
             ELSEIF (angle_list(iangle,species)%angle_potential_type == 'harmonic') THEN
                k = angle_list(iangle,species)%angle_param(1)
                theta0 = angle_list(iangle,species)%angle_param(2)
                
                CALL Get_Bond_Angle(iangle,molecule,species,theta)
                energy = k*(theta-theta0)**2
                ierr = 0

                ! Add more potential functions here.                
             ENDIF
             EXIT

          ELSEIF ( (angle_list(iangle,species)%atom1 == at3) .AND. &
               (angle_list(iangle,species)%atom3 == at1) ) THEN

             IF (angle_list(iangle,species)%angle_potential_type == 'fixed') THEN
                energy = 0.0_DP
                ierr = 0
             ELSEIF (angle_list(iangle,species)%angle_potential_type == 'harmonic') THEN
                k = angle_list(iangle,species)%angle_param(1)
                theta0 = angle_list(iangle,species)%angle_param(2)
                
                CALL Get_Bond_Angle(iangle,molecule,species,theta)
                energy = k*(theta-theta0)**2
                ierr = 0
                write(logunit,*) 'angle, theta, energy  ',iangle,theta*180.0_DP/PI, energy*atomic_to_kjmol
                ! Add more potential functions here.
             ENDIF
             EXIT

          ENDIF
          
       ENDIF

    ENDDO
    
    ! Test to make sure we successfuly found the right energy
    IF (ierr .NE. 0) THEN
       err_msg = ""
       err_msg(1) = "Unable to compute angle energy"
       CALL Clean_Abort(err_msg,'Compute_Angle_Energy')
    END IF

  END SUBROUTINE Compute_Angle_Energy
  !----------------------------------------------------------------------------------------------              

  SUBROUTINE Compute_Molecule_Angle_Energy(molecule,species,energy)
    !----------------------------------------------------------------------------------------------              
    ! This routine is passed a molecule and species index. It then computes the total
    ! bond angle energy of this molecule.  

    ! Called by:
    ! Calls: Get_Bond_Angle
    
    ! Written by: E. Maginn
    ! Date: Mon Nov 26 06:27:24 MST 2007
    ! Revision history:    
    !----------------------------------------------------------------------------------------------              
    USE Random_Generators
    ! Passed to 
    INTEGER :: molecule,species
    
    ! Returns
    REAL(DP) :: energy
    
    ! Local
    INTEGER :: iangle
    REAL(DP) :: k,theta0,theta,ea
  !----------------------------------------------------------------------------------------------              

    energy = 0.0_DP
    DO iangle=1,nangles(species)
       IF (angle_list(iangle,species)%int_angle_type == int_none) THEN
          ea = 0.0_DP
       ELSEIF (angle_list(iangle,species)%int_angle_type == int_harmonic) THEN
          k=angle_list(iangle,species)%angle_param(1)
          theta0 = angle_list(iangle,species)%angle_param(2)
          CALL Get_Bond_Angle(iangle,molecule,species,theta)
          ea = k*(theta-theta0)**2
          ! Add more potential functions here.
       ENDIF
       energy = energy + ea
    ENDDO

  END SUBROUTINE Compute_Molecule_Angle_Energy
  !----------------------------------------------------------------------------------------------              

  SUBROUTINE Compute_Dihedral_Energy(at1,at2,at3,at4,molecule,species,energy_dihed)
  !----------------------------------------------------------------------------------------------            
    ! This routine is passed 4 bonded atoms, and the molecule and species index. It figures out
    ! the type of torsion potential used and the parameters, and returns the energy of this torsion
    ! angle.

    ! Called by:
    ! Calls: Get_Dihedral_Angle
    
    ! Written by: E. Maginn
    ! Date: Mon Nov 26 09:33:51 MST 2007
    ! Revision history: 
    !
    ! 01/19/09 (JS) : Added CHARMM functional form
	! 12/8/12 (AV): Added AMBER functional form for multiple dihedrals
  !----------------------------------------------------------------------------------------------            
  USE Run_Variables
    INTEGER :: at1,at2,at3,at4,molecule,species
    REAL(DP) :: energy_dihed

    INTEGER :: i,ierr, nbr_dihed, idihed
    REAL(DP) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,phi,twophi,threephi
  !----------------------------------------------------------------------------------------------            
  ierr = 1
  energy_dihed = 0.0_DP

    ! Number of angles associated with one of the two central atoms
    nbr_dihed = dihedral_part_list(at2,species)%ndihedrals

    ! Loop over and locate which of these has at2 at either position 2 or 3.
    DO i=1,nbr_dihed
       Pos2or3:IF (dihedral_part_list(at2,species)%position(i) == 2) THEN

          ! Check if at1=1, at2=2, at3=3 and at4=4 OR at1=4, at2=3, at3=2 and at4=1
          idihed = dihedral_part_list(at2,species)%which_dihedral(i)
          
          Dihed1234:IF ( (dihedral_list(idihed,species)%atom1 == at1) .AND. &
               (dihedral_list(idihed,species)%atom3 == at3) .AND. &
               (dihedral_list(idihed,species)%atom4 == at4) ) THEN

             ! Located the angle at1-at2-at3-at4

             ! Determine type of potential, parameters, then energy
             DihedPot1234:IF (dihedral_list(idihed,species)%int_dipot_type == int_opls ) THEN
                CALL Get_Dihedral_Angle(idihed,molecule,species,phi)

                ! look up params and compute using OPLS. Note these were converted to 
                ! molecular units in input_routines.
                a0 = dihedral_list(idihed,species)%dihedral_param(1)
                a1 = dihedral_list(idihed,species)%dihedral_param(2)
                a2 = dihedral_list(idihed,species)%dihedral_param(3)
                a3 = dihedral_list(idihed,species)%dihedral_param(4)

                !Calculate the potential due to the torsional angle 
                twophi = 2.0_DP*phi
                threephi = 3.0_DP*phi
                energy_dihed =  a0 + a1*(1.0_DP+COS(phi)) + &
                     a2*(1.0_DP-COS(twophi)) + a3*(1.0_DP+COS(threephi))

                ierr = 0


             ELSEIF (dihedral_list(idihed,species)%int_dipot_type == int_none ) THEN
                energy_dihed = 0.0_DP
                ierr = 0    
                ! Add more potential functions here.
             ELSEIF (dihedral_list(idihed,species)%int_dipot_type == int_charmm ) THEN
                
                a0 = dihedral_list(idihed,species)%dihedral_param(1)
                a1 = dihedral_list(idihed,species)%dihedral_param(2)
                a2 = dihedral_list(idihed,species)%dihedral_param(3)
                
                energy_dihed = a0 * (1.0_DP + DCOS(a1*phi - a2))
                ierr = 0
			!AV: AMBER dihedral form	
             ELSEIF (dihedral_list(idihed,species)%int_dipot_type == int_amber ) THEN
                
                a0 = dihedral_list(idihed,species)%dihedral_param(1)
                a1 = dihedral_list(idihed,species)%dihedral_param(2)
                a2 = dihedral_list(idihed,species)%dihedral_param(3)
                a3 = dihedral_list(idihed,species)%dihedral_param(4)
                a4 = dihedral_list(idihed,species)%dihedral_param(5)
                a5 = dihedral_list(idihed,species)%dihedral_param(6)
                a6 = dihedral_list(idihed,species)%dihedral_param(7)
                a7 = dihedral_list(idihed,species)%dihedral_param(8)
                a8 = dihedral_list(idihed,species)%dihedral_param(9)
                !AV: commented out b/c 3 terms is usually suffice.
                !a9 = dihedral_list(idihed,species)%dihedral_param(10)
                !a10 = dihedral_list(idihed,species)%dihedral_param(11)
                !a11= dihedral_list(idihed,species)%dihedral_param(12)
                
                energy_dihed = a0 * (1.0_DP + DCOS(a1*phi - a2)) + &
                     a3 * (1.0_DP + DCOS(a4*phi - a5)) + &
                     a6 * (1.0_DP + DCOS(a7*phi - a8)) !+ &
                !a9 * (1.0_DP + DCOS(a10*phi - a11))
                
			 
             ENDIF DihedPot1234

             EXIT
             
          ENDIF Dihed1234


       ELSEIF  (dihedral_part_list(at2,species)%position(i) == 3) THEN

     ! Check if at1=1, at2=2, at3=3 and at4=4 OR at1=4, at2=3, at3=2 and at4=1
          idihed = dihedral_part_list(at2,species)%which_dihedral(i)
          
          Dihed4321:IF ( (dihedral_list(idihed,species)%atom4 == at1) .AND. &
               (dihedral_list(idihed,species)%atom2 == at3) .AND. &
               (dihedral_list(idihed,species)%atom1 == at4) ) THEN

             ! Located the angle at4-at3-at2-at1

             ! Determine type of potential, parameters, then energy
             DihedPot4321:IF (dihedral_list(idihed,species)%int_dipot_type == int_opls ) THEN
                CALL Get_Dihedral_Angle(idihed,molecule,species,phi)

                ! look up params and compute using OPLS. Note these were converted to 
                ! molecular units in input_routines.
                a0 = dihedral_list(idihed,species)%dihedral_param(1)
                a1 = dihedral_list(idihed,species)%dihedral_param(2)
                a2 = dihedral_list(idihed,species)%dihedral_param(3)
                a3 = dihedral_list(idihed,species)%dihedral_param(4)

                !Calculate the potential due to the torsional angle 
                twophi = 2.0_DP*phi
                threephi = 3.0_DP*phi
                energy_dihed =  a0 + a1*(1.0_DP+COS(phi)) + &
                     a2*(1.0_DP-COS(twophi)) + a3*(1.0_DP+COS(threephi))

                ierr = 0

             ELSEIF (dihedral_list(idihed,species)%int_dipot_type == int_none ) THEN
                energy_dihed = 0.0_DP
                ierr = 0    

             ELSEIF (dihedral_list(idihed,species)%int_dipot_type == int_charmm ) THEN

                a0 = dihedral_list(idihed,species)%dihedral_param(1)
                a1 = dihedral_list(idihed,species)%dihedral_param(2)
                a2 = dihedral_list(idihed,species)%dihedral_param(3)
                
                energy_dihed = a0 * (1.0_DP + DCOS(a1*phi - a2))
                ierr = 0

                ! Add more potential functions here.
			 ELSEIF (dihedral_list(idihed,species)%int_dipot_type == int_amber ) THEN

                a0 = dihedral_list(idihed,species)%dihedral_param(1)
                a1 = dihedral_list(idihed,species)%dihedral_param(2)
                a2 = dihedral_list(idihed,species)%dihedral_param(3)
                a3 = dihedral_list(idihed,species)%dihedral_param(4)
                a4 = dihedral_list(idihed,species)%dihedral_param(5)
                a5 = dihedral_list(idihed,species)%dihedral_param(6)
                a6 = dihedral_list(idihed,species)%dihedral_param(7)
                a7 = dihedral_list(idihed,species)%dihedral_param(8)
                a8 = dihedral_list(idihed,species)%dihedral_param(9)
              
                
                energy_dihed = a0 * (1.0_DP + DCOS(a1*phi - a2)) + &
                     a3 * (1.0_DP + DCOS(a4*phi - a5)) + &
                     a6 * (1.0_DP + DCOS(a7*phi - a8)) !+ &

                ierr = 0

             ENDIF DihedPot4321

             EXIT
             
          ENDIF Dihed4321

       ENDIF Pos2or3

    ENDDO
    
    ! Test to make sure we successfuly found the right energy
    IF (ierr .NE. 0) THEN
       err_msg = ""
       err_msg(1) = "Unable to compute dihedral energy"
       CALL Clean_Abort(err_msg,'Compute_Dihedral_Energy')
    END IF

  END SUBROUTINE Compute_Dihedral_Energy
  !----------------------------------------------------------------------------------------------            

  SUBROUTINE Compute_Molecule_Dihedral_Energy(molecule,species,energy_dihed)
   !----------------------------------------------------------------------------------------------              
    ! This routine is passed a molecule and species index. It then computes the total
    ! dihedral angle energy of this molecule.  

    ! Called by:
    ! Calls: Get_Dihedral_Angle
    
    ! Written by: E. Maginn
    ! Date: Mon Nov 26 10:01:22 MST 2007
    ! Revision history:
    ! AV: Added AMBER dihedral style:  12/8/12
    !----------------------------------------------------------------------------------------------              
  USE Run_Variables  
    ! Passed to 
    INTEGER :: molecule,species
    
    ! Returns
    REAL(DP) :: energy_dihed
    
    ! Local
    INTEGER :: idihed, atom1, atom2, atom3, atom4
    REAL(DP) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,edihed,phi,twophi,threephi
  !----------------------------------------------------------------------------------------------              

    energy_dihed = 0.0_DP
    DO idihed=1,ndihedrals(species)
       IF (dihedral_list(idihed,species)%int_dipot_type == int_none ) THEN
          edihed = 0.0_DP
       ELSEIF (dihedral_list(idihed,species)%int_dipot_type == int_opls ) THEN

          ! Check to see if the atoms of this dihedral exists. This is required
          ! for CBMC moves in which only a part of the molecule is present in
          ! the simulation

          atom1 = dihedral_list(idihed,species)%atom1
          atom2 = dihedral_list(idihed,species)%atom2
          atom3 = dihedral_list(idihed,species)%atom3
          atom4 = dihedral_list(idihed,species)%atom4

          IF ( .NOT. atom_list(atom1,molecule,species)%exist) CYCLE
          IF ( .NOT. atom_list(atom2,molecule,species)%exist) CYCLE
          IF ( .NOT. atom_list(atom3,molecule,species)%exist) CYCLE
          IF ( .NOT. atom_list(atom4,molecule,species)%exist) CYCLE

          
          a0 = dihedral_list(idihed,species)%dihedral_param(1)
          a1 = dihedral_list(idihed,species)%dihedral_param(2)
          a2 = dihedral_list(idihed,species)%dihedral_param(3)
          a3 = dihedral_list(idihed,species)%dihedral_param(4)

          CALL Get_Dihedral_Angle(idihed,molecule,species,phi)

          twophi = 2.0_DP*phi
          threephi = 3.0_DP*phi
          edihed =  a0 + a1*(1.0_DP+COS(phi)) + &
               a2*(1.0_DP-COS(twophi)) + a3*(1.0_DP+COS(threephi))


       ELSEIF (dihedral_list(idihed,species)%int_dipot_type == int_charmm ) THEN

          atom1 = dihedral_list(idihed,species)%atom1
          atom2 = dihedral_list(idihed,species)%atom2
          atom3 = dihedral_list(idihed,species)%atom3
          atom4 = dihedral_list(idihed,species)%atom4

          IF ( .NOT. atom_list(atom1,molecule,species)%exist) CYCLE
          IF ( .NOT. atom_list(atom2,molecule,species)%exist) CYCLE
          IF ( .NOT. atom_list(atom3,molecule,species)%exist) CYCLE
          IF ( .NOT. atom_list(atom4,molecule,species)%exist) CYCLE
          
          a0 = dihedral_list(idihed,species)%dihedral_param(1)
          a1 = dihedral_list(idihed,species)%dihedral_param(2)
          a2 = dihedral_list(idihed,species)%dihedral_param(3)

          CALL Get_Dihedral_Angle(idihed,molecule,species,phi)
          
          edihed = a0 * (1.0_DP + DCOS(a1*phi - a2))

!AV: AMBER dihedral style
       ELSEIF (dihedral_list(idihed,species)%int_dipot_type == int_amber ) THEN

          atom1 = dihedral_list(idihed,species)%atom1
          atom2 = dihedral_list(idihed,species)%atom2
          atom3 = dihedral_list(idihed,species)%atom3
          atom4 = dihedral_list(idihed,species)%atom4

          IF ( .NOT. atom_list(atom1,molecule,species)%exist) CYCLE
          IF ( .NOT. atom_list(atom2,molecule,species)%exist) CYCLE
          IF ( .NOT. atom_list(atom3,molecule,species)%exist) CYCLE
          IF ( .NOT. atom_list(atom4,molecule,species)%exist) CYCLE
          
          a0 = dihedral_list(idihed,species)%dihedral_param(1)
          a1 = dihedral_list(idihed,species)%dihedral_param(2)
          a2 = dihedral_list(idihed,species)%dihedral_param(3)
          a3 = dihedral_list(idihed,species)%dihedral_param(4)
          a4 = dihedral_list(idihed,species)%dihedral_param(5)
          a5 = dihedral_list(idihed,species)%dihedral_param(6)
          a6 = dihedral_list(idihed,species)%dihedral_param(7)
          a7 = dihedral_list(idihed,species)%dihedral_param(8)
          a8 = dihedral_list(idihed,species)%dihedral_param(9)
		  ! AV: I comment this out b/c it is not usually necessary.
		  !a9 = dihedral_list(idihed,species)%dihedral_param(10)
		  !a10 = dihedral_list(idihed,species)%dihedral_param(11)
		  !a11 = dihedral_list(idihed,species)%dihedral_param(12)

          CALL Get_Dihedral_Angle(idihed,molecule,species,phi)
          
          edihed = a0 * (1.0_DP + DCOS(a1*phi - a2)) + &
			a3 * (1.0_DP + DCOS(a4*phi - a5)) + &
				a6 * (1.0_DP + DCOS(a7*phi - a8)) !+ &
				!a9 * (1.0_DP + DCOS(a10*phi - a11))
		  
		  
       ELSEIF (dihedral_list(idihed,species)%int_dipot_type == int_harmonic ) THEN

          atom1 = dihedral_list(idihed,species)%atom1
          atom2 = dihedral_list(idihed,species)%atom2
          atom3 = dihedral_list(idihed,species)%atom3
          atom4 = dihedral_list(idihed,species)%atom4

          IF ( .NOT. atom_list(atom1,molecule,species)%exist) CYCLE
          IF ( .NOT. atom_list(atom2,molecule,species)%exist) CYCLE
          IF ( .NOT. atom_list(atom3,molecule,species)%exist) CYCLE
          IF ( .NOT. atom_list(atom4,molecule,species)%exist) CYCLE
          
          a0 = dihedral_list(idihed,species)%dihedral_param(1)
          a1 = dihedral_list(idihed,species)%dihedral_param(2)

          CALL Get_Dihedral_Angle(idihed,molecule,species,phi)

          IF(a1 .GT. 0.0_DP .AND. phi .LT.0) phi = phi + twoPI
          
          edihed = a0 * (phi - a1)**2
          
          ! Add more potential functions here.
       ENDIF
       energy_dihed = energy_dihed + edihed
    ENDDO


  END SUBROUTINE Compute_Molecule_Dihedral_Energy
  !----------------------------------------------------------------------------------------------              


  SUBROUTINE Compute_Improper_Energy(at1,at2,at3,at4,molecule,species,energy)
  !----------------------------------------------------------------------------------------------              
    ! This routine is passed 4 atoms making up an imprioper angle, as well as the molecule 
    ! and species indices. It computes the energy assuming a harmonic improper style of the form
    !                                E = k*(phi-phi0)^2 
    ! where k is a force constant in amu A^2/(ps^2 rad^2) and phi is the improper angle in radians.  
    ! Phi0 is the nominal improper angle. 
    ! If the 4 atoms in an improper quadruplet are ordered I,J,K,L then phi is the angle between the 
    ! plane of I,J,K and the plane of J,K,L. Alternatively, you can think of atoms J,K,L as being
    ! in a plane, and atom I above the plane, and phi as a measure of how far out-of-plane I is
    ! with respect to the other 3 atoms. Note that defining 4 atoms to interact in this way
    ! does not mean that bonds necessarily exist between I-J, J-K, or K-L, as they would in a 
    ! linear dihedral. Normally, the bonds I-J, I-K, I-L would exist for an improper to be 
    ! defined between the 4 atoms. 
    
    ! Called by:
    ! Calls: Get_Improper_Angle
    
    ! Written by: E. Maginn
    ! Date: Mon Nov 26 10:01:22 MST 2007
    ! Revision history:    
  !----------------------------------------------------------------------------------------------              
    INTEGER :: at1,at2,at3,at4,molecule,species
    REAL(DP) :: energy

    INTEGER :: ierr, iimprop
    REAL(DP) :: k,phi0,phi
  !----------------------------------------------------------------------------------------------              
    ierr = 1
    energy = 0.0_DP

    ! Loop over all improper angles for this species and search for a match at1-at2-at3-at4
    ! or at4-at3-at2-at1. We don't use a participation matrix, since there are typically few
    ! impropers, so the time savings is minimal

    DO iimprop = 1,nimpropers(species)

       AtomOrder:IF  ( (improper_list(iimprop,species)%atom1== at1) .AND. &
            (improper_list(iimprop,species)%atom2== at2) .AND. &
            (improper_list(iimprop,species)%atom3== at3) .AND. &
            (improper_list(iimprop,species)%atom4== at4) ) THEN

          Improp1234:IF (improper_list(iimprop,species)%improper_potential_type == 'harmonic') THEN
             k = improper_list(iimprop,species)%improper_param(1)
             phi0 = improper_list(iimprop,species)%improper_param(2)
             CALL Get_Improper_Angle(iimprop,molecule,species,phi)
             energy = k*(phi-phi0)**2
             ierr = 0
             EXIT
          ELSEIF (improper_list(iimprop,species)%improper_potential_type == 'none') THEN
             energy = 0.0_DP
             ierr = 0
             EXIT
             ! Add more potential functions here.
          ENDIF Improp1234

       ELSEIF ( (improper_list(iimprop,species)%atom1== at4) .AND. &
            (improper_list(iimprop,species)%atom2== at3) .AND. &
            (improper_list(iimprop,species)%atom3== at2) .AND. &
            (improper_list(iimprop,species)%atom4== at1) ) THEN

          Improp4321:IF (improper_list(iimprop,species)%improper_potential_type == 'harmonic') THEN
             k = improper_list(iimprop,species)%improper_param(1)
             phi0 = improper_list(iimprop,species)%improper_param(2)
             CALL Get_Improper_Angle(iimprop,molecule,species,phi)
             energy = k*(phi-phi0)**2
             ierr = 0
             EXIT
          ELSEIF (improper_list(iimprop,species)%improper_potential_type == 'none') THEN
             energy = 0.0_DP
             ierr = 0
             EXIT
             ! Add more potential functions here.
          ENDIF Improp4321

       ENDIF AtomOrder

    ENDDO

    ! Test to make sure we successfuly found the right energy
    IF (ierr .NE. 0) THEN
       err_msg = ""
       err_msg(1) = "Unable to compute improper energy"
       CALL Clean_Abort(err_msg,'Compute_Improper_Energy')
    END IF


  END SUBROUTINE Compute_Improper_Energy
  !----------------------------------------------------------------------------------------------              

  SUBROUTINE Compute_Molecule_Improper_Energy(molecule,species,energy)
  !----------------------------------------------------------------------------------------------              
    ! This routie is passed the molecule and species index, and returns the total improper
    ! energy of that molecule. Only "none" and "harminic" types are supported.

    ! Called by:
    ! Calls: Get_Improper_Angle
    !
    ! Written by: E. Maginn
    ! Date: Mon Nov 26 11:37:40 MST 2007
    ! Revision history
  !----------------------------------------------------------------------------------------------              
    INTEGER :: molecule,species,iimprop
    REAL(DP) :: energy
    REAL(DP) :: eimprop,k,phi0,phi,n_imp,d_imp
  !----------------------------------------------------------------------------------------------              
    energy = 0.0_DP
    DO iimprop=1,nimpropers(species)
       IF (improper_list(iimprop,species)%int_improp_type == int_none) THEN
          eimprop = 0.0_DP
       ELSEIF (improper_list(iimprop,species)%int_improp_type == int_harmonic) THEN
          k = improper_list(iimprop,species)%improper_param(1)
          phi0 = improper_list(iimprop,species)%improper_param(2)
          CALL Get_Improper_Angle(iimprop,molecule,species,phi)
          eimprop = k*(phi-phi0)**2
       ELSEIF (improper_list(iimprop,species)%int_improp_type == int_cvff) THEN
          k = improper_list(iimprop,species)%improper_param(1)
          d_imp = improper_list(iimprop,species)%improper_param(2)
          n_imp = improper_list(iimprop,species)%improper_param(3)
          CALL Get_Improper_Angle(iimprop,molecule,species,phi)
          eimprop = k*(1.0_DP + d_imp*DCOS(n_imp*phi))
       ELSE
          err_msg = ""
          err_msg(1) = "Unable to compute improper energy"
          CALL Clean_Abort(err_msg,'Compute_Molecular_Improper_Energy')
       ENDIF
       energy = energy + eimprop
    ENDDO
    
  END SUBROUTINE Compute_Molecule_Improper_Energy
  !----------------------------------------------------------------------------------------------              

  SUBROUTINE Compute_Atom_Nonbond_Energy(this_atom,this_molecule,this_species,E_intra_vdw, &
       E_inter_vdw,E_intra_qq,E_inter_qq,overlap)

    ! Computes the energy components between one particular atom and ALL others in its box, accounting
    ! for exclusions, scalings and existence. It returns two components of energy and the virial for the passed atom.
    ! Note that if this is used to compute the energy of several atoms, care must
    ! be taken to avoid "double counting" the energy. 
    ! This would most typically be used to compute energies for addition of atoms during CBMC growth.

    ! Note that the VDW energy (without LRC) is returned as is the real space part of the q-q interactions
    ! (for Ewald) or the direct sum 1-1 part for charge_sum_style = cut.

    ! if vdw_style == 'NONE' and charge_style='NONE' then this routine should really not be called.
    ! However, after much effort a value of 0 will be returned for the energy!
    !
    ! CALLED BY:
    ! CALLS: 
    ! Minimum_Image_Separation
    ! Clean_Abort
    ! Pair_Energy

    ! Written by: E. Maginn
    ! Date: Wed Nov 28 14:23:16 MST 2007
    ! Revision history:
    !
    ! 01/22/09 (JS) : Modified to separate out intra and intermolecule energy terms. Also
    !                 identity of molecules is obtained via locate array
    !
    ! 03/11/11 (JS) : Note that the loops for energy calculations are unrolled in this routine (Tom Rosch's code)
    !                 Need to figure out how to do openMP here so that other routines can be called much the same
    !                 way as done in Compute_Molecule_Nonbond_Inter_Energy
    !----------------------------------------------------------------------------------------------              

    INTEGER, INTENT(IN) :: this_atom,this_molecule,this_species
    REAL(DP), INTENT(OUT) :: E_intra_vdw,E_inter_vdw,E_intra_qq,E_inter_qq
    REAL(DP) :: E_intra_vdw_new,E_inter_vdw_new,E_intra_qq_new,E_inter_qq_new
    LOGICAL, INTENT(OUT) :: overlap
    
    INTEGER :: this_box,is,im,js,ia, mol_is, itype, jtype, rinteraction
    REAL(DP) :: rxij,ryij,rzij,rijsq,rxijp,ryijp,rzijp
    REAL(DP) :: Eij_vdw,Eij_qq
    REAL(DP) :: eps, sig, SigOverRsq, SigOverR6, SigOverR12
    REAL(DP) :: qi, qj, rij, erf_val, erfc_val, qsc
    REAL(DP) :: T, x, xsq, TP
    REAL(DP) :: rcom,rx,ry,rz
    REAL(DP) :: rcut, rcutsq
    REAL(DP) :: this_lambda_lj
    REAL(DP) :: SigOverR, SigOverRn, SigOverRm, mie_coeff,  mie_n, mie_m

    LOGICAL :: get_vdw,get_qq, f_intra_nrg, get_interaction

    REAL(DP), PARAMETER :: A1 = 0.254829592_DP, A2 = -0.284496736_DP
    REAL(DP), PARAMETER :: A3 = 1.421413741_DP, A4 = -1.453152027_DP
    REAL(DP), PARAMETER :: A5 = 1.061405429_DP, P = 0.3275911_DP


    !----------------------------------------------------------------------------------------------              
    E_inter_vdw = 0.0_DP
    E_intra_vdw = 0.0_DP
    E_inter_qq = 0.0_DP
    E_intra_qq = 0.0_DP
    E_inter_vdw_new = 0.0_DP
    E_intra_vdw_new = 0.0_DP
    E_inter_qq_new = 0.0_DP
    E_intra_qq_new = 0.0_DP

    f_intra_nrg = .FALSE.

    ! Set the box number this particular atom is in.
    this_box = molecule_list(this_molecule,this_species)%which_box
   
    ! Initialize flags which force a call to pair_energy
    get_vdw = .FALSE.
    get_qq = .FALSE.

    ! Initialize the overlap flag to false to indicate no overlap between atoms.
    overlap = .FALSE.
   
    IF (.NOT. atom_list(this_atom,this_molecule,this_species)%exist ) THEN
       err_msg = ""
       err_msg(1) = 'Attempt to compute energy of an atom that does not exist'
       CALL Clean_Abort(err_msg,'Compute_Atom_Nonbond_Energy')      
    ENDIF
     IF (int_vdw_sum_style(this_box) /= vdw_cut_tail)  THEN
!	write(*,*) "test"
!    IF(int_vdw_style(this_box) == vdw_lj .AND. int_vdw_sum_style(this_box) == vdw_cut_tail .AND. &
!       int_charge_style(this_box) == charge_coul .AND. int_charge_sum_style(this_box) == charge_ewald) THEN
!       write(*,*) "test"       
       IF (cbmc_flag) THEN
          rcut = rcut_cbmc(this_box)
       ELSE
          rcut = rcut_max(this_box)
       END IF
       
       rcutsq = rcut * rcut
       
       DO is=1,nspecies
          
          !$OMP PARALLEL DO DEFAULT(SHARED) &
          !$OMP PRIVATE(im,mol_is,rxijp,ryijp,rzijp,rijsq,rinteraction,T,x,xsq,TP,erfc_val,ia,f_intra_nrg) &
          !$OMP PRIVATE(itype,jtype,eps,sig,SigOverRsq,SigOverR6,SigOverR12,rij,erf_val,Eij_vdw,Eij_qq,qsc) &
          !$OMP SCHEDULE(DYNAMIC) &
          !$OMP REDUCTION(+:E_inter_vdw,E_inter_qq,E_intra_vdw,E_intra_qq)
          
          DO mol_is=1,nmolecules(is)
             
             IF(overlap) CYCLE
             
             im = locate(mol_is,is)
             
             IF( .NOT. molecule_list(im,is)%live ) CYCLE 
             
             rxijp = molecule_list(im,is)%xcom - molecule_list(this_molecule,this_species)%xcom
             ryijp = molecule_list(im,is)%ycom - molecule_list(this_molecule,this_species)%ycom
             rzijp = molecule_list(im,is)%zcom - molecule_list(this_molecule,this_species)%zcom
             
             IF (rxijp.gt.box_list(this_box)%hlength(1,1)) THEN
                rxijp = rxijp-box_list(this_box)%length(1,1)
             ELSEIF (rxijp.lt.-box_list(this_box)%hlength(1,1)) THEN
                rxijp = rxijp+box_list(this_box)%length(1,1)
             ENDIF

             IF (ryijp.gt.box_list(this_box)%hlength(2,2)) THEN
                ryijp = ryijp-box_list(this_box)%length(2,2)
             ELSEIF (ryijp.lt.-box_list(this_box)%hlength(2,2)) THEN
                ryijp = ryijp+box_list(this_box)%length(2,2)
             ENDIF
             
             IF (rzijp.gt.box_list(this_box)%hlength(3,3)) THEN
                rzijp = rzijp-box_list(this_box)%length(3,3)
             ELSEIF (rzijp.lt.-box_list(this_box)%hlength(3,3)) THEN
                rzijp = rzijp+box_list(this_box)%length(3,3)
             ENDIF

             rijsq = rxijp*rxijp + ryijp*ryijp + rzijp*rzijp
             rinteraction = rcut + molecule_list(im,is)%max_dcom &
                  + molecule_list(this_molecule,this_species)%max_dcom
             rinteraction = rinteraction * rinteraction
             
             IF (rijsq .GT. rinteraction) CYCLE 
             
             DO ia=1,natoms(is)
                
                IF (overlap) CYCLE
                IF (is == this_species .AND. im == this_molecule) THEN
                   
                   IF (ia == this_atom) CYCLE
                   IF (.NOT. atom_list(ia,im,is)%exist ) CYCLE 
                   
                   rxijp = atom_list(ia,im,is)%rxp - atom_list(this_atom,this_molecule,this_species)%rxp
                   ryijp = atom_list(ia,im,is)%ryp - atom_list(this_atom,this_molecule,this_species)%ryp
                   rzijp = atom_list(ia,im,is)%rzp - atom_list(this_atom,this_molecule,this_species)%rzp
                   
                   rijsq = rxijp*rxijp + ryijp*ryijp + rzijp*rzijp
                   
                   IF (rijsq <= rcut_lowsq) THEN
                      overlap = .true.
                   END IF

                ELSE
                   
                   rxijp = atom_list(ia,im,is)%rxp - atom_list(this_atom,this_molecule,this_species)%rxp
                   ryijp = atom_list(ia,im,is)%ryp - atom_list(this_atom,this_molecule,this_species)%ryp
                   rzijp = atom_list(ia,im,is)%rzp - atom_list(this_atom,this_molecule,this_species)%rzp
                   
                   IF (rxijp.gt.box_list(this_box)%hlength(1,1)) THEN
                      rxijp = rxijp-box_list(this_box)%length(1,1)
                   ELSEIF (rxijp.lt.-box_list(this_box)%hlength(1,1)) THEN
                      rxijp = rxijp+box_list(this_box)%length(1,1)
                   ENDIF
                   
                   IF (ryijp.gt.box_list(this_box)%hlength(2,2)) THEN
                      ryijp = ryijp-box_list(this_box)%length(2,2)
                   ELSEIF (ryijp.lt.-box_list(this_box)%hlength(2,2)) THEN
                      ryijp = ryijp+box_list(this_box)%length(2,2)
                   ENDIF
                   
                   IF (rzijp.gt.box_list(this_box)%hlength(3,3)) THEN
                      rzijp = rzijp-box_list(this_box)%length(3,3)
                   ELSEIF (rzijp.lt.-box_list(this_box)%hlength(3,3)) THEN
                      rzijp = rzijp+box_list(this_box)%length(3,3)
                   ENDIF
                   
                   rijsq = rxijp*rxijp + ryijp*ryijp + rzijp*rzijp
                   
                   IF (rijsq < rcut_lowsq) THEN
                      overlap = .true.
                   END IF
                   
                ENDIF

                IF(.NOT. overlap .AND. rijsq <= rcutsq) THEN         
                   
                   qsc = 1.0_DP
                   itype = nonbond_list(ia,is)%atom_type_number
                   jtype = nonbond_list(this_atom,this_species)%atom_type_number
                   eps = vdw_param1_table(itype,jtype)
                   sig = vdw_param2_table(itype,jtype)
                   qi = nonbond_list(ia,is)%charge
                   qj = nonbond_list(this_atom,this_species)%charge
                   
                   this_lambda_lj = 1.0_DP
                   
                   IF(is == this_species .AND. im == this_molecule) THEN
                      f_intra_nrg = .TRUE.
                      eps = eps * vdw_intra_scale(ia,this_atom,is)
                      qsc = charge_intra_scale(ia,this_atom,is)
                   ELSE
                      rij = SQRT(rijsq)
		      !IF (int_vdw_sum_style(this_box) == vdw_mie)  THEN
		      !	WRITE(*,*) "test"
                	!mie_n = mie_nlist(mie_Matrix(is,js))
                	!mie_m = mie_mlist(mie_Matrix(is,js))
                	!mie_coeff = mie_n/(mie_n-mie_m) *(mie_n/mie_m)**(mie_m/(mie_n-mie_m))
                	!SigOverR = sig/rij
                	!SigOverRn = SigOverR ** mie_n
                	!SigOverRm = SigOverR ** mie_m
                	!Eij_vdw =  mie_coeff * eps * (SigOverRn - SigOverRm)
		      !ELSE
                        SigOverRsq = (sig**2) / rijsq
                        SigOverR6  = SigOverRsq * SigOverRsq * SigOverRsq
                        SigOverR12 = SigOverR6 * SigOverR6
                        Eij_vdw = 4.0_DP * eps * (SigOverR12 - SigOverR6)
                      !ENDIF
                      
		      x = alpha_ewald(this_box) * rij
                      T = 1.0_DP / (1.0_DP + P*x)
                      xsq = x*x
                      TP = T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * A5))))
                      erfc_val = TP * EXP(-xsq)
                      
                      erf_val = 1.0_DP - erfc_val
                      
                      Eij_qq = (qi*qj/rij)*(qsc - erf_val)*charge_factor
                      
                      
                      IF (f_intra_nrg) THEN
                         
                         E_intra_vdw = E_intra_vdw + Eij_vdw
                         E_intra_qq = E_intra_qq + Eij_qq
                         f_intra_nrg = .FALSE.
                         
                      ELSE
                         
                         E_inter_vdw = E_inter_vdw + Eij_vdw
                         E_inter_qq = E_inter_qq + Eij_qq

                      END IF
                      
                   ENDIF
             
                ENDIF

             END DO

          END DO

       !$OMP END PARALLEL DO
       
          IF(overlap) RETURN
       
       END DO
    
    ELSE
       SpeciesLoop:DO is=1,nspecies
          
          MoleculeLoop:DO mol_is=1,nmolecules(is)
             
             im = locate(mol_is,is)
             
             ! Make sure that the molecule exists in the simulation
             
             IF( .NOT. molecule_list(im,is)%live ) CYCLE MoleculeLoop
             
             ! Only allow interactions in the same box
             IF (this_box /= molecule_list(im,is)%which_box) CYCLE MoleculeLoop
             
             ! Check tos see if atom is to interact with the molecule based
             ! on COM cutoff.
             
             CALL Check_Interaction(im,is,this_molecule,this_species,get_interaction,rcom,rx,ry,rz)
             
             IF ( .NOT. get_interaction) CYCLE MoleculeLoop
             
                AtomLoop:DO ia=1,natoms(is)
                   ! Test for intramolecular interaction
                   IF (.NOT. atom_list(ia,im,is)%exist ) CYCLE AtomLoop  
                   IF (is == this_species .AND. im == this_molecule) THEN
                      
                      IF (ia == this_atom) THEN
                         ! Avoid computing energy with self
                         CYCLE AtomLoop
                      ELSE
                         ! Intra energy. Do not apply PBC
                         
                         f_intra_nrg = .TRUE.
                         
                         IF ( .NOT. atom_list(ia,im,is)%exist) CYCLE AtomLoop
                         
                         
                         ! Find distance between this atom and all others in the system
                         rxij = atom_list(ia,im,is)%rxp - atom_list(this_atom,this_molecule,this_species)%rxp
                         ryij = atom_list(ia,im,is)%ryp - atom_list(this_atom,this_molecule,this_species)%ryp
                         rzij = atom_list(ia,im,is)%rzp - atom_list(this_atom,this_molecule,this_species)%rzp
                         
                         rijsq = rxij*rxij + ryij*ryij + rzij*rzij
                         
                         IF (rijsq <= rcut_lowsq) THEN
                            IF (.not.(l_bonded(ia,this_atom,is))) THEN
                               overlap = .true.
!                               write(*,*) 'bonded overlap', ia, this_atom, im, rijsq
                               RETURN
                            ENDIF
                         END IF
                      ENDIF
                      
                   ELSE                      
                      
                      ! It is an intermolecular energy so apply pbc. First compute the parent separation
                      rxijp = atom_list(ia,im,is)%rxp - atom_list(this_atom,this_molecule,this_species)%rxp
                      ryijp = atom_list(ia,im,is)%ryp - atom_list(this_atom,this_molecule,this_species)%ryp
                      rzijp = atom_list(ia,im,is)%rzp - atom_list(this_atom,this_molecule,this_species)%rzp
                      
                      ! Now get the minimum image separation 
                      CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp,rxij,ryij,rzij)
                      
                      rijsq = rxij*rxij + ryij*ryij + rzij*rzij
                      
                      IF (rijsq < rcut_lowsq) THEN
                         overlap = .true.
                         RETURN
                      END IF
                      
                   ENDIF
                   
                   CALL Energy_Test(rijsq,get_vdw,get_qq,this_box)
                   
                   ! Compute vdw and q-q energy using if required
                   IF (get_vdw .OR. get_qq) THEN 
                      
                      CALL Pair_Energy(rxij,ryij,rzij,rijsq,is,im,ia,this_species,this_molecule,this_atom,&
                           get_vdw,get_qq,Eij_vdw,Eij_qq)
                      
                      IF (f_intra_nrg) THEN
                         ! put the energies in intramolecular
                         E_intra_vdw = E_intra_vdw + Eij_vdw
                         E_intra_qq = E_intra_qq + Eij_qq
                         
                         f_intra_nrg = .FALSE.
                         
                      ELSE
                         !IF (CBMC_flag) THEN
                         !  IF ( .not. del_flag) THEN
                         !  IF (this_mcstep == 30985 ) THEN
                         !    write(301,*) ia, im, rijsq
                         ! END IF
                         ! END IF
                         ! END IF
                         
                         E_inter_vdw = E_inter_vdw + Eij_vdw
                         E_inter_qq = E_inter_qq + Eij_qq
                         
                      END IF
                      
                   ENDIF
                   
                END DO AtomLoop
                
             END DO MoleculeLoop
             
          END DO SpeciesLoop
          
       END IF
       
     END SUBROUTINE Compute_Atom_Nonbond_Energy

  !----------------------------------------------------------------------------------------------              

  SUBROUTINE Compute_Molecule_Nonbond_Intra_Energy(im,is,E_intra_vdw,E_intra_qq,intra_overlap,E_self)

    ! The subroutine calculates the intramolecular LJ potential energy and electrostatic
    ! energy of an entire molecule. The routine is based off the above routine 'Compute
    ! _Atom_Nonbond_Intra_Energy' and takes care of double counting by looping only over
    ! i+1 to natoms for ith atom interaction. 
    !
    ! CALLS
    ! 
    ! Pair_Energy
    !
    ! CALLED BY
    !
    ! Rigid_Dihedral_Change
    ! Angle_Distortion
    !
    !
    ! Written by Jindal Shah on 12/05/07
    !***************************************************************************************

    IMPLICIT NONE

    INTEGER :: ia, ja, im, is, this_box

    REAL(DP) :: E_intra_vdw, E_intra_qq, rxij, ryij, rzij, rijsq
    REAL(DP) :: E_intra_vdw_old, E_intra_qq_old
    REAL(DP) :: Eij_vdw, Eij_qq
    REAL(DP), OPTIONAL :: E_self
    
    REAL(DP), PARAMETER :: A1 = 0.254829592_DP, A2 = -0.284496736_DP
    REAL(DP), PARAMETER :: A3 = 1.421413741_DP, A4 = -1.453152027_DP
    REAL(DP), PARAMETER :: A5 = 1.061405429_DP, P = 0.3275911_DP

    LOGICAL :: get_vdw, get_qq, intra_overlap

    E_intra_vdw = 0.0_DP
    E_intra_qq = 0.0_DP
    E_intra_vdw_old = 0.0_DP
    E_intra_qq_old = 0.0_DP
    IF(PRESENT(E_self)) E_self = 0.0_DP
    
    ! loop over all the atoms in a molecule

    this_box = molecule_list(im,is)%which_box

    DO ia = 1, natoms(is)
       
       ! check to see if this atom exists
       ! Note 'im' is the linked number of the molecule of interest i.e locate(molecule,is)
       ! The checking for existence of a molecule may be unneccessary. 
       
       IF ( atom_list(ia,im,is)%exist) THEN
             
          DO ja = ia+1,natoms(is)
                
             ! make sure that the atom is present
             
             IF ( .NOT. atom_list(ja,im,is)%exist) CYCLE
             
             ! Find distance between this atom and all others in the system
             rxij = atom_list(ia,im,is)%rxp - atom_list(ja,im,is)%rxp
             ryij = atom_list(ia,im,is)%ryp - atom_list(ja,im,is)%ryp
             rzij = atom_list(ia,im,is)%rzp - atom_list(ja,im,is)%rzp
             
             rijsq = rxij*rxij + ryij*ryij + rzij*rzij
             
             IF (rijsq <= rcut_lowsq) THEN
                IF (.not.(l_bonded(ia,ja,is))) THEN
                   intra_overlap = .true.
                   RETURN
                ENDIF
             END IF
             
             CALL Energy_Test(rijsq,get_vdw,get_qq,this_box)
             
             IF(cbmc_flag.and.species_list(is)%L_Coul_CBMC) THEN
                get_qq=.false.
             ENDIF
             
             ! Compute vdw and q-q energy using if required
             IF (get_vdw .OR. get_qq) THEN 
                
                CALL Pair_Energy(rxij,ryij,rzij,rijsq,is,im,ia,is,im,ja,get_vdw,get_qq,Eij_vdw,Eij_qq)
                
                E_intra_vdw = E_intra_vdw + Eij_vdw
                E_intra_qq  = E_intra_qq + Eij_qq
                
             END IF
             
          END DO
          
       END IF
       
    END DO
    
  END SUBROUTINE Compute_Molecule_Nonbond_Intra_Energy
  !----------------------------------------------------------------------------------------------------------

  SUBROUTINE Compute_Molecule_Nonbond_Inter_Energy(im,is,E_inter_vdw,E_inter_qq,overlap)

    !**************************************************************************************************
    ! This subroutine computes interatomic LJ and charge interactions as well as virials associated
    ! with these interactions. 
    !
    ! CALLS
    ! 
    ! Minimum_Image_Separation
    ! Pair_Energy
    ! Clean_Abort
    !
    ! CALLED BY
    !
    ! Translate
    ! Rotation
    ! Rigid_Dihedral_Change
    ! Angle_Distortion
    !
    ! Written by Jindal Shah on 12/07/07
    !***************************************************************************************************

    IMPLICIT NONE

!    !$ include 'omp_lib.h'

    INTEGER, INTENT(IN):: im, is
    REAL(DP), INTENT(OUT) :: E_inter_vdw, E_inter_qq
    LOGICAL :: overlap
    !---------------------------------------------------------------------------------------------------

    INTEGER  :: ispecies, imolecule, this_box, this_locate
    
    REAL(DP) :: Eij_vdw, Eij_qq
    REAL(DP) :: eps
    REAL(DP) :: rcom, rx, ry, rz

    REAL(DP), PARAMETER :: A1 = 0.254829592_DP, A2 = -0.284496736_DP
    REAL(DP), PARAMETER :: A3 = 1.421413741_DP, A4 = -1.453152027_DP
    REAL(DP), PARAMETER :: A5 = 1.061405429_DP, P = 0.3275911_DP

    LOGICAL :: get_interaction

    INTEGER :: locate_1, locate_2

    LOGICAL :: l_pair_store 
    LOGICAL :: my_overlap, shared_overlap

    E_inter_vdw = 0.0_DP
    E_inter_qq = 0.0_DP
    overlap = .FALSE.
    my_overlap = .FALSE.
    shared_overlap = .false.
    
    this_box = molecule_list(im,is)%which_box

    l_pair_store = .FALSE.

    IF (l_pair_nrg .AND. (.NOT. cbmc_flag)) l_pair_store = .TRUE.

    IF (l_pair_store) THEN
       ! find out the location correspoding to im
       
       IF (is == 1) THEN
          locate_1 = im
       ELSE
          locate_1 = SUM(nmolecules(1:is-1)) + im
       END IF

    END IF

    speciesLoop: DO ispecies = 1, nspecies
       
       !$OMP PARALLEL DO DEFAULT(SHARED) &
       !$OMP PRIVATE(imolecule,this_locate,locate_2,get_interaction) &
       !$OMP PRIVATE(rx,ry,rz,rcom,Eij_vdw,Eij_qq) &
       !$OMP SCHEDULE(DYNAMIC) &
       !$OMP REDUCTION(+:E_inter_vdw,E_inter_qq) & 
       !$OMP REDUCTION(.OR.:my_overlap)  

       
       moleculeLoop: DO imolecule = 1, nmolecules(ispecies)
          
          IF(shared_overlap) CYCLE
          
          this_locate = locate(imolecule,ispecies)
          
          IF (ispecies == is .AND. this_locate == im) CYCLE moleculeLOOP
          
          IF (molecule_list(this_locate,ispecies)%which_box /= this_box) CYCLE moleculeLOOP
          
          ! make sure that the molecule is currently part of the simulation
          
          IF(.NOT. molecule_list(this_locate,ispecies)%live) CYCLE moleculeLOOP
          
          IF (l_pair_store) THEN
             
             IF (ispecies == 1) THEN
                locate_2 = this_locate
             ELSE
                locate_2 = SUM(nmolecules(1:ispecies-1)) + this_locate
             END IF
             
             pair_nrg_vdw(locate_1,locate_2) = 0.0_DP
             pair_nrg_vdw(locate_2,locate_1) = 0.0_DP
             
             pair_nrg_qq(locate_1,locate_2) = 0.0_DP
             pair_nrg_qq(locate_2,locate_1) = 0.0_DP
             
          END IF
          
          ! Determine if any atoms of these two molecules will interact
          CALL Check_Interaction(im,is,this_locate,ispecies,get_interaction,rcom,rx,ry,rz) 

          IF (.NOT. get_interaction) CYCLE moleculeLOOP       
          
          
          CALL Compute_Molecule_Pair_Interaction(im,is,this_locate,ispecies,this_box, &
               Eij_vdw,Eij_qq,my_overlap)
          
          IF( my_overlap) shared_overlap = .true.
          
          E_inter_vdw = E_inter_vdw + Eij_vdw
          E_inter_qq  = E_inter_qq + Eij_qq
          
       END DO moleculeLoop
       !$OMP END PARALLEL DO
       
       IF(shared_overlap) THEN
          overlap = .true.
          RETURN
       ENDIF
       
    END DO speciesLoop
    
  END SUBROUTINE Compute_Molecule_Nonbond_Inter_Energy
  !------------------------------------------------------------------------------------------
  !------------------------------------------------------------------------------------------
  SUBROUTINE Pair_Energy &
       (rxij,ryij,rzij,rijsq,is,im,ia,js,jm,ja,get_vdw,get_qq,Eij_vdw,Eij_qq)

    ! LJ potential: Eij = 4*epsilon(i,j) * [ (sigma(i,j)/rij)^12 - (sigma(i,j)/rij)^6 ]

    ! Computes the vdw and q-q pair energy between atoms ia and ja of molecules im and jm
    ! and species is and js, given their separation rijsq. I have passed each component of 
    ! separation nut right now this is unnecessary. 
    ! It also computes the real space part of the Ewald sum if necessary.

    ! Called by: Compute_Atom_Nonbond_Energy
    ! Calls: CFC_LJ_Scaling
    !        Ewald_Real
  !------------------------------------------------------------------------------------------
    ! Passed to
    REAL(DP) :: rxij,ryij,rzij,rijsq
    INTEGER :: is,im,ia,js,jm,ja,ibox
    LOGICAL :: get_vdw,get_qq

    ! Returned
    REAL(DP) :: Eij_vdw,Eij_qq

    ! Local
    INTEGER :: itype,jtype, this_box
    REAL(DP) :: eps,sig,SigOverRsq,SigOverR6,SigOverR12
    REAL(DP) :: SigOverRsq_shift,SigOverR6_shift,SigOverR12_shift
    REAL(DP) :: roffsq_rijsq, roffsq_rijsq_sq, factor2, fscale
    REAL(DP) :: SigOverR, SigOverRn, SigOverRm, mie_coeff, rij,  mie_n, mie_m, rij_shift, SigOverR_shift, SigOverRn_shift, SigOverRm_shift, rcut_vdw
!    REAL(DP) :: Eij_vdw_check
    Real(DP) :: qi,qj, qsc
    REAL(DP) :: this_lambda, RsqOverSig, R6OverSig, factorLJ
    REAL(DP) :: RsqOverSig_Shift, RsqOverSig6_Shift, factorLJ_Shift

    LOGICAL :: fraction

  !------------------------------------------------------------------------------------------
    Eij_vdw = 0.0_DP
    Eij_qq = 0.0_DP
    fraction = .false.
    this_lambda = 1.0_DP

    ibox = molecule_list(im,is)%which_box

    ! If either atom is not yet present, then don't try to compute an energy
    ExistCheck: IF (atom_list(ia,im,is)%exist .AND. atom_list(ja,jm,js)%exist) THEN

       ! Determine atom type indices
       itype = nonbond_list(ia,is)%atom_type_number
       jtype = nonbond_list(ja,js)%atom_type_number
       
       VDW_calculation: IF (get_vdw) THEN

          LJ_12_6_calculation: IF (int_vdw_style(1) == vdw_lj) THEN
             ! For now, assume all interactions are the same. Use the lookup table created in Compute_Nonbond_Table
             eps = vdw_param1_table(itype,jtype)
             sig = vdw_param2_table(itype,jtype)

             ! Apply intramolecular scaling if necessary
             IF (is == js .AND. im == jm) THEN
                
                ! This controls 1-2, 1-3, and 1-4 interactions
                
                eps = eps * vdw_intra_scale(ia,ja,is)

             ENDIF
                
             SigOverRsq = (sig**2) / rijsq
             SigOverR6  = SigOverRsq * SigOverRsq * SigOverRsq
             SigOverR12 = SigOverR6 * SigOverR6
             
             IF (int_vdw_sum_style(ibox) == vdw_charmm) THEN

                ! use the form for modified LJ potential
                
                Eij_vdw = eps * (SigOverR12 - 2.0_DP * SigOverR6)

             ELSEIF (int_vdw_sum_style(ibox) == vdw_cut .OR. int_vdw_sum_style(ibox) == vdw_cut_tail) THEN
                 
                   Eij_vdw = 4.0_DP * eps * (SigOverR12 - SigOverR6)

             ELSEIF (int_vdw_sum_style(ibox) == vdw_cut_shift) THEN
                
                SigOverRsq = (sig**2)/rijsq
                SigOverR6 = SigOverRsq * SigOverRsq * SigOverRsq
                SigOverR12 = SigOverR6 * SigOverR6

                
                SigOverRsq_shift = sig**2/rcut_vdwsq(ibox)
                SigOverR6_shift = SigOverRsq_shift * SigOverRsq_shift * SigOverRsq_shift
                SigOverR12_shift = SigOverR6_shift * SigOverR6_shift
                
                Eij_vdw = 4.0_DP * eps * ( (SigOverR12 - SigOverR6) - (SigOverR12_shift - SigOverR6_shift) )

                
             ELSEIF (int_vdw_sum_style(ibox) == vdw_cut_switch) THEN
                
                Eij_vdw = 4.0_DP * eps * (SigOverR12 - SigOverR6)
                                            
                IF ( (rijsq < ron_switch_sq(ibox) )) THEN
                   
                   fscale = 1.0_DP
                   
                ELSEIF ( rijsq <= roff_switch_sq(ibox)) THEN
                   
                   roffsq_rijsq = roff_switch_sq(ibox) - rijsq
                   
                   roffsq_rijsq_sq = roffsq_rijsq * roffsq_rijsq
                   
                   factor2 = switch_factor2(ibox) + 2.0_DP * rijsq
                   
                   fscale = roffsq_rijsq_sq * factor2 * switch_factor1(ibox)
                   
                   Eij_vdw = fscale * Eij_vdw

                ELSE
                   
                   fscale = 0.0_DP
                   Eij_vdw = 0.0_DP
                ENDIF

             ELSEIF (int_vdw_sum_style(ibox) == vdw_mie) THEN

                rij = SQRT(rijsq)
		rcut_vdw = SQRT(rcut_vdwsq(ibox))
                
                mie_n = mie_nlist(mie_Matrix(is,js))
                mie_m = mie_mlist(mie_Matrix(is,js))
                mie_coeff = mie_n/(mie_n-mie_m) * (mie_n/mie_m)**(mie_m/(mie_n-mie_m))
                SigOverR = sig/rij
		SigOverR_shift = sig/rcut_vdw
		!use cut-shift potential
		SigOverRn_shift = SigOverR_shift ** mie_n
		SigOverRm_shift = SigOverR_shift ** mie_m
                SigOverRn = SigOverR ** mie_n
                SigOverRm = SigOverR ** mie_m
                Eij_vdw =  mie_coeff * eps * ((SigOverRn - SigOverRm) - (SigOverRn_shift - SigOverRm_shift))
                
		!print *, Eij_vdw
		!READ(*,*)
             ENDIF
             
             
          ENDIF LJ_12_6_calculation
          
       ENDIF VDW_calculation
  
       qq_calculation: IF (get_qq) THEN
          
          qi = nonbond_list(ia,is)%charge
          qj = nonbond_list(ja,js)%charge

          
          IF (int_charge_sum_style(ibox) == charge_cut .OR. igas_flag) THEN
             ! Apply charge scaling for intramolecular energies
             qsc = 1.0_DP
             IF ( is == js .AND. im == jm ) THEN
                qsc = charge_intra_scale(ia,ja,is)
             END IF
             Eij_qq = qsc*charge_factor*(qi*qj)/SQRT(rijsq)
          ELSEIF (int_charge_sum_style(ibox) == charge_ewald .AND. ( .NOT. igas_flag) ) THEN
             ! Real space Ewald part
             this_box = molecule_list(im,is)%which_box
             CALL Ewald_Real(ia,im,is,qi,ja,jm,js,qj,rijsq,Eij_qq,this_box)
             
             ! self and recipricoal parts need to be computed as total energy differences between original
             ! configuration and the perturbed configuration. These terms are thus added on after all atoms 
             ! have been moved. 

          ENDIF

       ENDIF qq_calculation
       
    ENDIF ExistCheck

  END SUBROUTINE Pair_Energy

  !------------------------------------------------------------------------------------------
  SUBROUTINE Ewald_Real(ia,im,is,qi,ja,jm,js,qj,rijsq,Eij,ibox)
  !------------------------------------------------------------------------------------------
    ! Real space part of the Ewald sum between atoms 1a and jq with 
    ! charges qi and qj.

    ! For intramolecular scaling, we want Efull - (1-w)*qiqj/rij
    ! where Efull is the non-scaled real-space sum and w is the intramolecular
    ! scaling between i and j. We thus have to subtract 
    ! (1-lambda)qiqj/rij from this to get the correct scaling. 

    ! Called by: Pair_Energy
    ! Calls: CFC_qq_scaling
  !------------------------------------------------------------------------------------------
    INTEGER :: ia,im,is,ja,jm,js,ibox
    REAL(DP) :: qi,qj,qsc,rijsq,rij,erf_val
    REAL(DP) :: Eij

    qsc = 1.0_DP
    ibox = molecule_list(im,is)%which_box 

    ! Apply intramolecular scaling if necessary
    IF (is == js .AND. im == jm) THEN
       
       ! Intramolecular charge scaling
       qsc = charge_intra_scale(ia,ja,is)
              
    ENDIF

    ! Real space part: This does the intrascaling correct. For cfc intra,
    ! we use full scaling. I think for CFC inter, we simply scale the actual 
    ! value of the charge, but do NOT scale it here for intra interactions.
    ! Come back to this later.

    rij = SQRT(rijsq)
    ! May need to protect against very small rijsq
    erf_val = 1.0_DP - erfc(alpha_ewald(ibox) * rij)
    Eij = (qi*qj/rij)*(qsc - erf_val)*charge_factor
!                   IF(en_flag) THEN
!                      WRITE(60,"(4I4,2F8.5,F24.12)") ia, ja, jm, js, qi,qj, Eij
!                   END IF

!------------------------------------------------------------------------------
  CONTAINS

    FUNCTION erfc(x)
      !**************************************************************************
      !                                                                         *
      ! Calculate the complementary error function for  a number
      !                                                                         *
      !**************************************************************************

      REAL(DP) :: erfc
      REAL(DP), PARAMETER :: A1 = 0.254829592_DP, A2 = -0.284496736_DP
      REAL(DP), PARAMETER :: A3 = 1.421413741_DP, A4 = -1.453152027_DP
      REAL(DP), PARAMETER :: A5 = 1.061405429_DP, P = 0.3275911_DP
      REAL(DP) :: T, x, xsq, TP

      T = 1.0_DP / (1.0_DP + P*x)
      xsq = x*x

      TP = T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * A5))))

      erfc = TP * EXP(-xsq)

    END FUNCTION erfc
!------------------------------------------------------------------------------

  END SUBROUTINE Ewald_Real

  !********************************************************************************************

  SUBROUTINE Ewald_Reciprocal_Lattice_Vector_Setup(this_box)
    !******************************************************************************************
    ! This subroutine sets up the reciprocal lattice vector constants required in the reciprocal
    ! space summation. Note that these constants need to be recomputed every time a volume
    ! change move is attempted. 
    ! Based on the APSS code, ewald_setup.f90 
    !
    ! Added by Jindal Shah on 12/05/07
    !
    !******************************************************************************************
    
    USE Type_Definitions
    USE Run_Variables
    
    IMPLICIT NONE
    
    INTEGER :: nx, ny, nz, this_box, kvecs, kx_max, ky_max, kz_max

    REAL(DP) :: const_val, hcutsq, x, y, z
    REAL(DP) :: hx_val, hy_val, hz_val, hsq_val

    ! Total number of k vectors
    kvecs = 1

    const_val = 1.0_DP/(4.0_DP * alpha_ewald(this_box) * alpha_ewald(this_box))
    hcutsq = h_ewald_cut(this_box) * h_ewald_cut(this_box)

    IF (box_list(this_box)%int_box_shape == int_cell .OR. box_list(this_box)%int_box_shape == int_ortho) THEN
       
       ! The most general definition for a wave-vector is h = 2*pi*TRANSPOSE(cell_matrix)^-1)*n
       ! where h is the wave vector and n is a vector of integers

       ! We will use symmetry about the x to calculate only half the wave vectors

       DO nz = -20, 20
          DO ny = -20, 20
             DO nx = 0, 20
                
                ! Exclude the possiblity for the central  simulation box where h = 0

                IF ( (nx == 0) .AND. (ny == 0) .AND. (nz == 0) ) CYCLE

                x = REAL(nx,DP)
                y = REAL(ny,DP)
                z = REAL(nz,DP)
                
                hx_val = twoPI * (box_list(this_box)%length_inv(1,1)*x + &
                     box_list(this_box)%length_inv(2,1)*y + box_list(this_box)%length_inv(3,1)*z)
                hy_val = twoPI * (box_list(this_box)%length_inv(1,2)*x + & 
                     box_list(this_box)%length_inv(2,2)*y + box_list(this_box)%length_inv(3,2)*z)
                hz_val = twoPI * (box_list(this_box)%length_inv(1,3)*x + &
                     box_list(this_box)%length_inv(2,3)*y + box_list(this_box)%length_inv(3,3)*z)

                hsq_val = hx_val * hx_val + hy_val * hy_val + hz_val * hz_val
                
                IF (hsq_val < hcutsq) THEN
                   
                   hx(kvecs,this_box) = hx_val
                   hy(kvecs,this_box) = hy_val
                   hz(kvecs,this_box) = hz_val
                   hsq(kvecs,this_box) = hsq_val
                   
                   ! if x /= 0, multipy the constant by 2 for symmetry
                   
                   IF ( nx == 0 ) THEN
                      
                      Cn(kvecs,this_box) = twoPI / box_list(this_box)%volume &
                           * DEXP ( -hsq(kvecs,this_box) * const_val ) / hsq(kvecs,this_box)
                      
                   ELSE
                      
                      Cn(kvecs,this_box) = 2.0_DP * twoPI / box_list(this_box)%volume &
                           * DEXP ( -hsq(kvecs,this_box) * const_val ) / hsq(kvecs,this_box)
                      
                   END IF
                   
                   kvecs = kvecs + 1
                   
                END IF

             END DO
          END DO
       END DO

    ELSE

       ! if it is an orthogonal box, then h vectors are simply hx = twoPI * nx / Lx and so on

       ! we will determine the number of reciprocal space vectors needed in each direction

       kz_max = INT ( (h_ewald_cut(this_box) * box_list(this_box)%basis_length(3))/twoPI ) + 1
       ky_max = INT ( (h_ewald_cut(this_box) * box_list(this_box)%basis_length(2))/twoPI ) + 1
       kx_max = INT ( (h_ewald_cut(this_box) * box_list(this_box)%basis_length(1))/twoPI ) + 1


       DO nz = -kz_max, kz_max
          DO ny = -ky_max, ky_max
             DO nx = 0, kx_max

                IF ( (nx == 0) .AND. (ny == 0) .AND. (nz == 0)) CYCLE

                hx_val = twoPI * REAL(nx,DP)/box_list(this_box)%basis_length(1)
                hy_val = twoPI * REAL(ny,DP)/box_list(this_box)%basis_length(2)
                hz_val = twoPI * REAL(nz,DP)/box_list(this_box)%basis_length(3)
                
                hsq_val = hx_val * hx_val + hy_val * hy_val + hz_val * hz_val
                
                IF (hsq_val < hcutsq ) THEN
                   
                   hx(kvecs,this_box) = hx_val
                   hy(kvecs,this_box) = hy_val
                   hz(kvecs,this_box) = hz_val
                   hsq(kvecs,this_box) = hsq_val
                   
                   !hsq(kvecs,this_box) = hx(kvecs,this_box)*hx(kvecs,this_box) + &
                   ! hy(kvecs,this_box)*hy(kvecs,this_box) + hz(kvecs,this_box)*hz(kvecs,this_box)
                   
                   ! if x /= 0, multipy the constant by 2 for symmetry
                   IF ( nx == 0 ) THEN
                      
                      Cn(kvecs,this_box) = twoPI / box_list(this_box)%volume &
                           * DEXP ( - hsq(kvecs,this_box) * const_val ) / hsq(kvecs,this_box)
                      
                   ELSE
                      
                      Cn(kvecs,this_box) = 2.0_DP * twoPI / box_list(this_box)%volume &
                           * DEXP ( - hsq(kvecs,this_box) * const_val ) / hsq(kvecs,this_box)
                      
                   END IF
                   
                   kvecs = kvecs + 1
                   
                END IF
                
                IF ( kvecs - 1 > maxk) THEN
                   err_msg = ""
                   err_msg(1) = 'Total number of k vectors exceeded'
                   CALL Clean_Abort(err_msg,'Ewald_Reciprocal_Lattice_Vector_Setup')
                END IF
                
             END DO
          END DO
       END DO
       
    END IF
    

    
    ! nvecs points to where the next wave_vector should be written, i. e. it is too high by 1
    nvecs(this_box) = kvecs - 1
    ! Note that at this point we do not allocate the memory for cos_sum and sin_sum arrays.
    ! It will have to be decided by the maximum number of k vectors encountered in all the boxes

    
  END SUBROUTINE Ewald_Reciprocal_Lattice_Vector_Setup
  !*******************************************************************************************

  SUBROUTINE Compute_System_Ewald_Reciprocal_Energy(this_box)
    !*****************************************************************************************
    ! This subroutine computes the sin and cos sum terms for the calculation of reciprocal
    ! energy of the input box. 
    !
    ! Based on APSS code reciprocal_ewald.f90
    !
    ! Added by Jindal Shah on 12/05/07
    ! Modified by Tom Rosch on 06/11/09 (See Wheeler, Mol. Phys. 1997 Vol. 92 pg. 55)
    !
    ! Modified on 04/28/09 to account for fractional molecules. The routine computes 
    !     energy of the integer and fractional particles with scaled charges. However, long
    !     range interactions of a fractional molecule with itself is computed with unscaled
    !     charges.
    ! 
    !
    !*************************************************************************
    !*****************************************************************************************
    
    USE Type_Definitions
    USE Run_Variables

    IMPLICIT NONE

!    !$ include 'omp_lib.h'

    INTEGER :: i, is, im, ia, this_locate, this_box

    REAL(DP) :: un, const_val
    REAL(DP) :: cn_val, hx_val, hy_val, hz_val, charge, hdotr, energy_temp

    REAL(DP) :: this_scaling, cos_sum_locate, sin_sum_locate


    ! individual k-space vector stuff

    INTEGER, ALLOCATABLE :: im_locate(:,:)

    ! loop over all the k vectors of this box

    const_val = 1.0_DP/(2.0_DP * alpha_ewald(this_box) * alpha_ewald(this_box))
    energy(this_box)%ewald_reciprocal = 0.0_DP
    
    
    energy_temp = 0.0_DP
    !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
    cos_sum(:,this_box) = 0.0_DP
    sin_sum(:,this_box) = 0.0_DP
    !$OMP END PARALLEL WORKSHARE

    ALLOCATE(im_locate(MAXVAL(nmolecules),nspecies))


    DO is = 1, nspecies

       DO im = 1, nmolecules(is)

          this_locate = locate(im,is)
          
          IF (.NOT. molecule_list(this_locate,is)%live) CYCLE

          IF (molecule_list(this_locate,is)%which_box /= this_box) CYCLE

          IF (is == 1) THEN
             im_locate(im,is) = this_locate
          ELSE
             im_locate(im,is) = SUM(nmolecules(1:is-1)) + this_locate
          END IF

       END DO

    END DO

    DO i = 1, nvecs(this_box)

       cn_val = Cn(i,this_box)
       hx_val = hx(i,this_box)
       hy_val = hy(i,this_box)
       hz_val = hz(i,this_box)

       ! Loop over all the atoms present in this box to calculate sin_sum and cos_sum
       ! for each of the vector

       DO is = 1, nspecies
          IF( .NOT. has_charge(is) ) CYCLE
          DO im = 1, nmolecules(is)
             this_locate = locate(im,is)

             IF( .NOT. molecule_list(this_locate,is)%live) CYCLE

             IF( molecule_list(this_locate,is)%which_box /= this_box ) CYCLE

    
             cos_sum_locate = 0.0_DP
             sin_sum_locate = 0.0_DP
             
             DO ia = 1, natoms(is)
                
                ! compute hdotr 
                
                hdotr = hx_val * atom_list(ia,this_locate,is)%rxp + &
                        hy_val * atom_list(ia,this_locate,is)%ryp + &
                        hz_val * atom_list(ia,this_locate,is)%rzp
                
                charge = nonbond_list(ia,is)%charge
                
                cos_sum_locate = cos_sum_locate + charge * DCOS(hdotr)
                sin_sum_locate = sin_sum_locate + charge * DSIN(hdotr)
                
             END DO

             ! Note that only the molecules that belong to 'this_box'
             ! has its cos_mol and sin_mol vectors changed
             
             
             cos_mol(i,im_locate(im,is)) = cos_sum_locate
             sin_mol(i,im_locate(im,is)) = sin_sum_locate
             
             ! Compute charge scaling. It needs to be done only for the first
             ! k vector. We will store this and then use it for subsequent
             ! k vectors.

             this_scaling = 1.0_DP

             
             cos_sum(i,this_box) = cos_sum(i,this_box) + this_scaling * cos_sum_locate
             sin_sum(i,this_box) = sin_sum(i,this_box) + this_scaling * sin_sum_locate
             
          END DO
       END DO
       
       un = cn_val * (cos_sum(i,this_box) * cos_sum(i,this_box) + sin_sum(i,this_box) * sin_sum(i,this_box))         
       energy_temp = energy_temp + un
       
    END DO
    
    energy(this_box)%ewald_reciprocal = energy(this_box)%ewald_reciprocal + energy_temp
    energy(this_box)%ewald_reciprocal = energy(this_box)%ewald_reciprocal * charge_factor
    
        
  END SUBROUTINE Compute_System_Ewald_Reciprocal_Energy
  !*************************************************************************************************

  SUBROUTINE Compute_Ewald_Reciprocal_Energy_Difference(im,im_prev,is,this_box,move_flag,V_recip_difference)
    !************************************************************************************************
    ! The subroutine computes the difference in Ewald reciprocal space energy for a given move.
    !/
    ! We will develop this routine for a number of moves.
    !
    ! Translation of COM
    ! Rotation about COM
    ! Angle Distortion
    ! Rigid Dihedral rotation
    ! Molecule insertion
    ! Molecule Deletion
    !***********************************************************************************************

    USE Type_Definitions
    USE Run_Variables
    
    IMPLICIT NONE

!    !$ include 'omp_lib.h'

    INTEGER, INTENT(IN) :: im, im_prev, is, this_box
    INTEGER, INTENT(IN) :: move_flag
    REAL(DP), INTENT(OUT) :: V_recip_difference

    ! Note that im_prev is included here in the anticipation that the routine will be used
    ! for CFC move. At present the input value of im_prev is immaterial.


    ! Local variables
    
    REAL(DP) :: const_val
    INTEGER :: i, ia

    REAL(DP) :: hdotr_new

    REAL(DP) :: cos_sum_im, cos_sum_im_o, sin_sum_im, sin_sum_im_o

    ! storage stuff

    INTEGER :: im_locate, im_prev_locate

    ! get the location of im and im_prev

    IF (is==1) THEN
       im_locate = im
       im_prev_locate = im_prev
    ELSE
       im_locate = SUM(nmolecules(1:is-1)) + im
       im_prev_locate = SUM(nmolecules(1:is-1)) + im_prev
    END IF

    V_recip_difference = 0.0_DP
    const_val = 1.0_DP/(2.0_DP * alpha_ewald(this_box) * alpha_ewald(this_box))

    !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
    cos_sum_old(1:nvecs(this_box),this_box) = cos_sum(1:nvecs(this_box),this_box)
    sin_sum_old(1:nvecs(this_box),this_box) = sin_sum(1:nvecs(this_box),this_box)
    !$OMP END PARALLEL WORKSHARE

    IF ( move_flag == int_translation .OR. move_flag == int_rotation .OR. move_flag == int_intra ) THEN


       ! only the particle coordinates change. Therefore, the contribution of cos(hdotr) and
       ! sin(hdotr) of the old coordinates will be subtracted off for each of reciprocal vectors
       ! and corresponding terms for the new coordinates are added.

       ! Note that the flag INTRA will refer to any of the moves that correspond to the 
       ! intramolecular DOF change.
       
           
       !$OMP PARALLEL DO DEFAULT(SHARED) &
       !$OMP PRIVATE(i,ia,cos_sum_im,sin_sum_im) &
       !$OMP PRIVATE(cos_sum_im_o, sin_sum_im_o) &
       !$OMP PRIVATE(hdotr_new) &
       !$OMP SCHEDULE(STATIC) &
       !$OMP REDUCTION(+:v_recip_difference)
       DO i = 1, nvecs(this_box)

          cos_sum_im = 0.0_DP
          sin_sum_im = 0.0_DP

          DO ia = 1,natoms(is)

             ! let us compute the old and new hdotr

             hdotr_new = hx(i,this_box) * atom_list(ia,im,is)%rxp + &
                         hy(i,this_box) * atom_list(ia,im,is)%ryp + &
                         hz(i,this_box) * atom_list(ia,im,is)%rzp
             
             cos_sum_im = cos_sum_im + nonbond_list(ia,is)%charge * DCOS(hdotr_new)
             sin_sum_im = sin_sum_im + nonbond_list(ia,is)%charge * DSIN(hdotr_new)


          END DO

          cos_sum_im_o = cos_mol(i,im_locate)
          sin_sum_im_o = sin_mol(i,im_locate)

          cos_sum(i,this_box) = cos_sum(i,this_box) + cos_sum_im - cos_sum_im_o
          sin_sum(i,this_box) = sin_sum(i,this_box) + sin_sum_im - sin_sum_im_o

          v_recip_difference = v_recip_difference + cn(i,this_box) * (cos_sum(i,this_box) * &
               cos_sum(i,this_box) + sin_sum(i,this_box) * sin_sum(i,this_box))

          ! set the molecules cos and sin terms to the one calculated here
          cos_mol(i,im_locate) = cos_sum_im
          sin_mol(i,im_locate) = sin_sum_im

       END DO
       !$OMP END PARALLEL DO

       v_recip_difference = v_recip_difference * charge_factor - energy(this_box)%ewald_reciprocal

       RETURN

    ELSE IF ( move_flag == int_deletion) THEN

       ! We need to subtract off the cos(hdotr) and sin(hdor) for each of the k vectors.

       !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
       cos_sum(1:nvecs(this_box),this_box) = cos_sum(1:nvecs(this_box),this_box) - &
            cos_mol(1:nvecs(this_box),im_locate)
       sin_sum(1:nvecs(this_box),this_box) = sin_sum(1:nvecs(this_box),this_box) - &
            sin_mol(1:nvecs(this_box),im_locate)
       !$OMP END PARALLEL WORKSHARE

       !$OMP PARALLEL DO DEFAULT(SHARED) &
       !$OMP PRIVATE(i) &
       !$OMP SCHEDULE(STATIC) &
       !$OMP REDUCTION(+:v_recip_difference)
       DO i = 1, nvecs(this_box)
             
          v_recip_difference = v_recip_difference + cn(i,this_box) * ( cos_sum(i,this_box) * cos_sum(i,this_box) + &
                                     sin_sum(i,this_box) * sin_sum(i,this_box) )

       END DO
       !$OMP END PARALLEL DO
       v_recip_difference = v_recip_difference * charge_factor
 
    ELSE IF ( move_flag == int_insertion ) THEN

       !$OMP PARALLEL DO DEFAULT(SHARED) &
       !$OMP PRIVATE(i, ia, hdotr_new) &
       !$OMP SCHEDULE(STATIC) &
       !$OMP REDUCTION(+:v_recip_difference) 

       DO i = 1, nvecs(this_box)

          cos_mol(i,im_locate) = 0.0_DP
          sin_mol(i,im_locate) = 0.0_DP
          
          DO ia = 1, natoms(is)

             ! Compute the new hdotr vector

             hdotr_new = hx(i,this_box) * atom_list(ia,im,is)%rxp + &
                         hy(i,this_box) * atom_list(ia,im,is)%ryp + &
                         hz(i,this_box) * atom_list(ia,im,is)%rzp

             cos_mol(i,im_locate) = cos_mol(i,im_locate) +  nonbond_list(ia,is)%charge * DCOS(hdotr_new)
             sin_mol(i,im_locate) = sin_mol(i,im_locate) +  nonbond_list(ia,is)%charge * DSIN(hdotr_new)

          END DO

          cos_sum(i,this_box) = cos_sum(i,this_box) + cos_mol(i,im_locate)
          sin_sum(i,this_box) = sin_sum(i,this_box) + sin_mol(i,im_locate)

  
          v_recip_difference = v_recip_difference + cn(i,this_box) * ( cos_sum(i,this_box) * cos_sum(i,this_box) + &
                                     sin_sum(i,this_box) * sin_sum(i,this_box) )

       END DO

       !$OMP END PARALLEL DO

       v_recip_difference = v_recip_difference * charge_factor

    END IF

  END SUBROUTINE Compute_Ewald_Reciprocal_Energy_Difference
  !********************************************************************************************

  SUBROUTINE Compute_System_Ewald_Self_Energy(this_box)
    ! ******************************************************************************************
    ! This subroutine calculates the constant term that arises from particles interacting with
    ! themselves in the reciprocal space. The subroutine needs to be called only once as it
    ! is a constant term as long as the particles and their charges remain the same.
    !
    !*******************************************************************************************

    USE Type_Definitions
    USE Run_Variables
    
    IMPLICIT NONE

    INTEGER :: is,im, this_locate, ia,  this_box
    REAL(DP) :: charge

    energy(this_box)%ewald_self = 0.0_DP
    
    DO is = 1, nspecies

      imLOOP: DO im = 1, nmolecules(is)

          this_locate = locate(im,is)

          ! sum only those molecules that are in this_box.
          
          IF (.NOT. molecule_list(this_locate,is)%live) CYCLE imLOOP

          IF ( molecule_list(this_locate,is)%which_box /= this_box ) CYCLE imLOOP

          DO ia = 1, natoms(is)
                
             ! obtain the charge
             charge = nonbond_list(ia,is)%charge 
             energy(this_box)%ewald_self = energy(this_box)%ewald_self + (charge*charge)
             
          END DO

       END DO imLOOP
          
    END DO

    energy(this_box)%ewald_self = energy(this_box)%ewald_self * alpha_ewald(this_box) / rootPI
    energy(this_box)%ewald_self = - energy(this_box)%ewald_self * charge_factor

    ! Note that the ewald self constant computed here is slightly different than what is
    ! computed in APSS. It has a negative sign and is already multiplied by a the charge_factor
    ! for proper unit conversion

  END SUBROUTINE Compute_System_Ewald_Self_Energy
  !********************************************************************************************

  SUBROUTINE Compute_Ewald_Self_Energy_Difference(im,is,this_box,move_flag,V_self_difference)

    !*******************************************************************************************
    ! This subroutine calculates the difference in self Ewald energy for the input molecule(s)
    ! The routine needs to be called for the following moves. The two input molecules are im and
    ! im_prev. Note that im_prev is included anticipating that CFC moves will be implemented.
    ! In all other cases, the identity
    ! of im_prev is of no consequence. The difference is computed
    ! for molecules in 'this_box' while V_self_difference is returned to the calling routine.
    !
    ! Insertion
    ! Deletion
    !
    !**********************************************************************************************

    USE Type_Definitions
    USE Run_Variables

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: im, is, this_box
    INTEGER, INTENT(IN) :: move_flag
    
    REAL(DP), INTENT(OUT) :: V_self_difference

    INTEGER :: ia
    
    V_self_difference = 0.0_DP

    IF ( move_flag == int_insertion ) THEN

       DO ia = 1, natoms(is)
          V_self_difference = V_self_difference + nonbond_list(ia,is)%charge * nonbond_list(ia,is)%charge
       END DO

       V_self_difference = - charge_factor * V_self_difference * alpha_ewald(this_box)/rootPI

    ELSE IF ( move_flag == int_deletion ) THEN

       DO ia = 1,natoms(is)
          V_self_difference = V_self_difference + nonbond_list(ia,is)%charge*nonbond_list(ia,is)%charge
       END DO

       V_self_difference = charge_factor * V_self_difference * alpha_ewald(this_box)/rootPI

    END IF

  END SUBROUTINE Compute_Ewald_Self_Energy_Difference
  !****************************************************************************************************
  !**************************************************************************************************
  SUBROUTINE Compute_Total_System_Energy(this_box,intra_flag,overlap)
    !**************************************************************************************************
    ! The subroutine calculates the total energy of a given box. The identity of the box is passed to
    ! the routine along with the intra_flag to indicate whether intramolecular computation is required.
    ! The flag will mostly be set to true except in the case of volume change move that is designed so
    ! that the intramolecular DOFs do not change. 
    !
    !***************************************************************************************************

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: this_box
    LOGICAL, INTENT(IN) :: intra_flag

    !----------------------------------------------------------------------------------------------------

    INTEGER :: im, is, this_im, im_1, im_2, is_1, is_2, this_im_1, this_im_2

    REAL(DP) :: v_molecule_bond, v_molecule_angle, v_molecule_dihedral, v_molecule_improper
    REAL(DP) :: vlj_molecule_intra, vqq_molecule_intra
    REAL(DP) :: vlj_pair, vqq_pair, e_lrc
    REAL(DP) :: rcom, rx, ry, rz
    REAL(DP) :: E_inter_vdw, E_inter_qq
    REAL(DP) :: v_selfrf, v_molecule_selfrf
    REAL(DP) :: rijsq
    REAL(DP) :: v_bond, v_angle, v_dihedral, v_intra, v_intra_vdw, v_intra_qq, v_improper
    
    LOGICAL :: overlap, get_interaction,intra_overlap

    INTEGER :: locate_1, locate_2
    LOGICAL :: l_pair_store, my_overlap, shared_overlap

    my_overlap = .FALSE.
    shared_overlap = .FALSE.
    overlap = .FALSE.

    ! Initialize the energies

    energy(this_box)%total = 0.0_DP
    energy(this_box)%inter_vdw = 0.0_DP
    energy(this_box)%inter_q = 0.0_DP
    energy(this_box)%ewald_reciprocal = 0.0_DP
    ! Compute the intramolecular energy of the system if the flag is set.

    IF (intra_flag) THEN

       energy(this_box)%intra = 0.0_DP
       energy(this_box)%bond  = 0.0_DP
       energy(this_box)%angle = 0.0_DP
       energy(this_box)%dihedral = 0.0_DP
       energy(this_box)%improper = 0.0_DP
       energy(this_box)%intra_vdw = 0.0_DP
       energy(this_box)%intra_q = 0.0_DP
       energy(this_box)%erf_self = 0.0_DP

       DO is = 1, nspecies
          v_intra = 0.0_DP
          v_bond= 0.0_DP
          v_angle= 0.0_DP
          v_dihedral= 0.0_DP
          v_improper = 0.0_DP
          v_intra_vdw= 0.0_DP
          v_intra_qq = 0.0_DP
          v_selfrf = 0.0_DP
          !$OMP PARALLEL DO DEFAULT(SHARED) &
          !$OMP SCHEDULE(DYNAMIC) &
          !$OMP PRIVATE(im, this_im, v_molecule_bond, v_molecule_angle, v_molecule_dihedral) &
          !$OMP PRIVATE(v_molecule_improper,vlj_molecule_intra,vqq_molecule_intra, intra_overlap) &
          !$OMP REDUCTION(+:v_intra,v_bond, v_angle, v_dihedral,v_improper, v_intra_vdw, v_intra_qq, v_selfrf)  
          imLoop:DO im = 1,nmolecules(is)
             
             this_im = locate(im,is)
             ! Check to see the molecule belongs to this_box
             IF( molecule_list(this_im,is)%which_box /= this_box ) CYCLE imLOOP
             
             ! Check to see if the molecule is alive 
             IF( .NOT. molecule_list(this_im,is)%live ) CYCLE imLOOP
             IF (SHARED_OVERLAP) CYCLE imLOOP

             CALL Compute_Molecule_Bond_Energy(this_im,is,v_molecule_bond)
             CALL Compute_Molecule_Angle_Energy(this_im,is,v_molecule_angle)
             CALL Compute_Molecule_Dihedral_Energy(this_im,is,v_molecule_dihedral)
             CALL Compute_Molecule_Improper_Energy(this_im,is,v_molecule_improper)

             intra_overlap = .FALSE.
             IF (int_charge_sum_style(this_box) == charge_ewald) THEN
                CALL Compute_Molecule_Nonbond_Intra_Energy(this_im,is,vlj_molecule_intra,vqq_molecule_intra,intra_overlap, &
                     v_molecule_selfrf)
             ELSE
                CALL Compute_Molecule_Nonbond_Intra_Energy(this_im,is,vlj_molecule_intra,vqq_molecule_intra,intra_overlap)
             END IF

             IF (intra_overlap) THEN
                SHARED_OVERLAP = .TRUE.
             END IF

             v_intra = v_intra + v_molecule_bond + v_molecule_angle + &
                                      v_molecule_dihedral + v_molecule_improper 
             v_bond = v_bond + v_molecule_bond
             v_angle = v_angle + v_molecule_angle
             v_dihedral = v_dihedral + v_molecule_dihedral
             v_improper = v_improper + v_molecule_improper
             v_intra_vdw = v_intra_vdw + vlj_molecule_intra 
             v_intra_qq   = v_intra_qq   + vqq_molecule_intra
             IF (int_charge_sum_style(this_box) == charge_ewald) THEN
                v_selfrf  = v_selfrf + v_molecule_selfrf
             END IF

          END DO imLoop
          !$OMP END PARALLEL DO
          IF (SHARED_OVERLAP) THEN
             overlap = .TRUE.
             RETURN
          END IF

          energy(this_box)%intra = energy(this_box)%intra +  v_intra
          energy(this_box)%bond = energy(this_box)%bond + v_bond
          energy(this_box)%angle = energy(this_box)%angle + v_angle
          energy(this_box)%dihedral = energy(this_box)%dihedral + v_dihedral
          energy(this_box)%improper = energy(this_box)%improper + v_improper
          energy(this_box)%intra_vdw = energy(this_box)%intra_vdw + v_intra_vdw
          energy(this_box)%intra_q   = energy(this_box)%intra_q   + v_intra_qq
          IF (int_charge_sum_style(this_box) == charge_ewald) THEN
             energy(this_box)%erf_self = energy(this_box)%erf_self + v_selfrf
          END IF
       END DO

    END IF

    ! The total system energy. Note that intra_flag is not used for this calculation, beacuse, if the flag
    ! is true, we already computed the energy, if not we are using the old global energy (the routine
    ! did not modify the energy).

    energy(this_box)%total = energy(this_box)%total + energy(this_box)%intra + energy(this_box)%intra_vdw + &
         energy(this_box)%intra_q
    
    ! Calculate the total intermolecular energy of the system. The calculation is divided into two parts
    ! The first part computes the interaction between the molecules of the same species, while the second
    ! bit obtains the interaction between molecules of different species.

    l_pair_store = .FALSE.
    IF (l_pair_nrg .AND. (.NOT. cbmc_flag)) l_pair_store = .TRUE.
    
    DO is = 1, nspecies
       imLOOP1: DO im_1 = 1, nmolecules(is)
          this_im_1 = locate(im_1,is)
          ! is it in this box?
          IF ( molecule_list(this_im_1,is)%which_box /= this_box ) CYCLE imLOOP1
          ! is alive?
          IF ( .NOT. molecule_list(this_im_1,is)%live ) CYCLE imLOOP1
          
          IF (l_pair_store) CALL Get_Position_Alive(this_im_1, is, locate_1) 
      
          E_inter_vdw = 0.0_DP
          E_inter_qq  = 0.0_DP
          
          !$OMP PARALLEL DO DEFAULT(SHARED) &
          !$OMP SCHEDULE(DYNAMIC) &
          !$OMP PRIVATE(im_2, this_im_2, locate_2, get_interaction) &
          !$OMP PRIVATE(rcom, rx, ry, rz, vlj_pair, vqq_pair) &
          !$OMP PRIVATE(my_overlap) &
          !$OMP REDUCTION(+: E_inter_vdw, E_inter_qq) 
          
          imLOOP2: DO im_2 = im_1 + 1, nmolecules(is)
             this_im_2 = locate(im_2,is)
             ! allow interactions only in the box
             IF ( molecule_list(this_im_2,is)%which_box /= this_box ) CYCLE imLOOP2
             ! is it alive?
             IF ( .NOT. molecule_list(this_im_2,is)%live ) CYCLE imLOOP2
             IF (SHARED_OVERLAP) CYCLE imLOOP2
             
             IF (l_pair_store) THEN
                CALL Get_Position_Alive(this_im_2,is,locate_2)
                
                pair_nrg_vdw(locate_1,locate_2) = 0.0_DP
                pair_nrg_vdw(locate_2,locate_1) = 0.0_DP
                
                pair_nrg_qq(locate_1,locate_2) = 0.0_DP
                pair_nrg_qq(locate_2,locate_1) = 0.0_DP

             END IF
             
             CALL Check_Interaction(this_im_1,is,this_im_2,is,get_interaction,rcom,rx,ry,rz)
             
            ! rijsq = rcom * rcom 
             
             IF (.NOT. get_interaction) CYCLE imLoop2
             ! Compute the intermolecular interactions between these two molecules
             
             CALL Compute_Molecule_Pair_Interaction(this_im_1,is,this_im_2,is,this_box,vlj_pair, &
                  vqq_pair,my_overlap)
             
             !             IF (overlap) RETURN
             IF (my_overlap) THEN
                SHARED_OVERLAP = .true.
             END IF
             
             E_inter_vdw  = E_inter_vdw + vlj_pair
             E_inter_qq   = E_inter_qq   + vqq_pair
             
          END DO imLOOP2
          !$OMP END PARALLEL DO
          IF (SHARED_OVERLAP) THEN
             overlap = .true.
             RETURN
          ENDIF

          energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + E_inter_vdw
          energy(this_box)%inter_q = energy(this_box)%inter_q + E_inter_qq
          
       END DO imLOOP1
    END DO
    
    ! Now compute the interaction with the molecules between different species

    DO is_1 = 1, nspecies
       imLOOP3: DO im_1 = 1, nmolecules(is_1)
          this_im_1 = locate(im_1,is_1)
          IF( molecule_list(this_im_1,is_1)%which_box /= this_box ) CYCLE imLOOP3
          IF( .NOT. molecule_list(this_im_1,is_1)%live ) CYCLE imLOOP3
          
          IF (l_pair_store) CALL Get_Position_Alive(this_im_1,is_1,locate_1)
          
          DO is_2 = is_1 + 1, nspecies
             E_inter_vdw = 0.0_DP
             E_inter_qq  = 0.0_DP
             
             !$OMP PARALLEL DO DEFAULT(SHARED) &
             !$OMP SCHEDULE(DYNAMIC) &
             !$OMP PRIVATE(im_2, this_im_2, locate_2, get_interaction) &
             !$OMP PRIVATE(rcom, rx, ry, rz, vlj_pair, vqq_pair) &
             !$OMP PRIVATE(my_overlap) &
             !$OMP REDUCTION(+: E_inter_vdw, E_inter_qq) 
             
             imLOOP4: DO im_2 = 1,nmolecules(is_2)
                this_im_2 = locate(im_2,is_2)
                IF ( molecule_list(this_im_2,is_2)%which_box /= this_box ) CYCLE imLOOP4
                IF ( .NOT. molecule_list(this_im_2,is_2)%live ) CYCLE imLOOP4

                IF (SHARED_OVERLAP) CYCLE imLOOP4
                
                IF (l_pair_store) THEN
                   CALL Get_Position_Alive(this_im_2,is_2,locate_2)
                   
                   pair_nrg_vdw(locate_1,locate_2) = 0.0_DP
                   pair_nrg_vdw(locate_2,locate_1) = 0.0_DP
                   
                   pair_nrg_qq(locate_1,locate_2) = 0.0_DP
                   pair_nrg_qq(locate_2,locate_1) = 0.0_DP
                   
                END IF
                
                ! Check to see if the interaction needs to be computed between the molecules
                CALL Check_Interaction(this_im_1,is_1,this_im_2,is_2,get_interaction,rcom,rx,ry,rz)
!                rijsq = rcom * rcom

                IF (.NOT. get_interaction ) CYCLE imLOOP4

                ! Note that this call will modify the pair interaction energies
                ! if l_pair_nrg variable is .TRUE.

                CALL Compute_Molecule_Pair_Interaction(this_im_1,is_1,this_im_2,is_2,this_box,vlj_pair, &
                     vqq_pair,my_overlap)
                
                IF (my_overlap) THEN
                   SHARED_OVERLAP = .true.
                END IF
                
                E_inter_vdw  = E_inter_vdw + vlj_pair
                E_inter_qq   = E_inter_qq   + vqq_pair                
             
             END DO imLOOP4
             !$OMP END PARALLEL DO 
             IF (SHARED_OVERLAP) THEN
                overlap = .true.
                RETURN
             ENDIF
             
             energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + E_inter_vdw
             energy(this_box)%inter_q = energy(this_box)%inter_q + E_inter_qq
             
          END DO
          
       END DO imLOOP3
    END DO
    
    energy(this_box)%total = energy(this_box)%total + energy(this_box)%inter_vdw
    energy(this_box)%total = energy(this_box)%total + energy(this_box)%inter_q

    ! Compute the reciprocal and self energy terms of the electrostatic energies if flag for Ewald is set.

    IF (int_charge_sum_style(this_box) == charge_ewald) THEN

       ! Ewald reciprocal energy -- Note that the call changes global V_ewald_reciprocal(this_box)


       CALL Compute_System_Ewald_Reciprocal_Energy2(this_box)

       energy(this_box)%total = energy(this_box)%total + energy(this_box)%ewald_reciprocal
       
       ! Ewald self energy -- Note that the call changes global V_ewald_self(this_box)
       

       CALL Compute_System_Ewald_Self_Energy(this_box)

       energy(this_box)%ewald_self = energy(this_box)%ewald_self - energy(this_box)%erf_self
       energy(this_box)%total = energy(this_box)%total + energy(this_box)%ewald_self 

    END IF

    ! Long range correction if it is required
    IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
       CALL Compute_LR_Correction(this_box,e_lrc)
       ! add to the correction to the total energy of the system
       energy(this_box)%lrc = e_lrc
       energy(this_box)%total = energy(this_box)%total + energy(this_box)%lrc
    END IF

  END SUBROUTINE Compute_Total_System_Energy
  !*********************************************************************************************************

  SUBROUTINE Compute_Molecule_Pair_Interaction(im_1,is_1,im_2,is_2,this_box,vlj_pair,vqq_pair,overlap)
    !**********************************************************************************************************
    ! The subroutine returns the interaction energy of the input molecule with another molecule. Thus,
    ! it computes the intermolecular vdw and electrostatic interactions. 
    !**********************************************************************************************************
  
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: im_1, is_1, im_2, is_2, this_box
    REAL(DP), INTENT(OUT) :: vlj_pair,vqq_pair
    LOGICAL, INTENT(OUT) :: overlap
    !----------------------------------------------------------------------------------------------------------

    INTEGER :: ia, ja 

    REAL(DP) :: rxijp, ryijp, rzijp, rxij, ryij, rzij, rijsq, Eij_vdw, Eij_qq

    LOGICAL :: get_vdw, get_qq

    INTEGER :: locate_im_1, locate_im_2

    vlj_pair = 0.0_DP
    vqq_pair = 0.0_DP
    overlap = .FALSE.

    DO ia = 1, natoms(is_1)

       IF (.NOT. atom_list(ia,im_1,is_1)%exist) CYCLE

       DO ja = 1, natoms(is_2)
          ! Obtain the minimum image separation
          
          IF ( .NOT. atom_list(ja,im_2,is_2)%exist) CYCLE

          rxijp = atom_list(ia,im_1,is_1)%rxp - atom_list(ja,im_2,is_2)%rxp
          ryijp = atom_list(ia,im_1,is_1)%ryp - atom_list(ja,im_2,is_2)%ryp
          rzijp = atom_list(ia,im_1,is_1)%rzp - atom_list(ja,im_2,is_2)%rzp
          
          ! Now get the minimum image separation 
          CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp,rxij,ryij,rzij)

          rijsq = rxij*rxij + ryij*ryij + rzij*rzij

          IF( rijsq < rcut_lowsq ) THEN
!             WRITE(*,*) im_1,is_1,im_2,is_2,rijsq,rcut_lowsq
             overlap = .TRUE.
             RETURN
          END IF
          
          ! Now figure out what needs to be computed, then call pair_energy

          CALL Energy_Test(rijsq,get_vdw,get_qq,this_box)          

          ! Compute vdw and q-q energy using if required
          IF(cbmc_flag.and.(.not.species_list(is_1)%L_Coul_CBMC)) THEN
             get_qq=.false. 
          ENDIF 
          IF (get_vdw .OR. get_qq) THEN 

             CALL Pair_Energy(rxij,ryij,rzij,rijsq,is_1,im_1,ia,is_2,im_2,ja,&
                  get_vdw,get_qq,Eij_vdw,Eij_qq)

             vlj_pair = vlj_pair + Eij_vdw
             vqq_pair = vqq_pair + Eij_qq

          END IF

       END DO

    END DO

   

    IF (l_pair_nrg) THEN

       IF ( .NOT. cbmc_flag ) THEN

       

          ! if here then, there was no overlap between im_1 and im_2
          ! update the interaction energy of the pair
          ! first find out the position of im_2 in the pair interaction energy

          CALL Get_Position_Alive(im_1,is_1,locate_im_1)
          CALL Get_Position_Alive(im_2,is_2,locate_im_2)
          
          pair_nrg_vdw(locate_im_1,locate_im_2) = vlj_pair 
          pair_nrg_vdw(locate_im_2,locate_im_1) = vlj_pair 
          
          pair_nrg_qq(locate_im_1,locate_im_2) = vqq_pair 
          pair_nrg_qq(locate_im_2,locate_im_1) = vqq_pair 
          

       END IF
       
    END IF
  
  END SUBROUTINE Compute_Molecule_Pair_Interaction

  !----------------------------------------------------------------------------------------------
  SUBROUTINE Compute_LR_Correction(this_box, e_lrc)
    !*************************************************************************************************************
    ! The subroutine calculates the long range correction for the given box.
    !
    !**************************************************************************************************************
    INTEGER, INTENT(IN) :: this_box
    REAL(DP), INTENT(OUT) :: e_lrc
    
    INTEGER ::  ia, ja
    REAL(DP) :: epsij, sigij, sigij2, sigij6, sigij12

    REAL(DP) :: e_lrc_ia_ja

    e_lrc = 0.0_DP
    
    DO ia = 1, nbr_atomtypes
       
       e_lrc_ia_ja = 0.0_DP
       
       DO ja = 1, nbr_atomtypes
          
          epsij = vdw_param1_table(ia,ja)
          sigij = vdw_param2_table(ia,ja)
          
          sigij2 = sigij*sigij
          
          sigij6 = sigij2*sigij2*sigij2 
          
          sigij12 = sigij6*sigij6

          e_lrc_ia_ja = e_lrc_ia_ja + nint_beads(ja,this_box) * &
               4.0_DP * epsij * (sigij12 /(9.0_DP*rcut9(this_box)) - &
               (sigij6 / (3.0_DP*rcut3(this_box))))
          
       END DO

       e_lrc = e_lrc + REAL( nint_beads(ia,this_box), DP ) * e_lrc_ia_ja
    END DO

    e_lrc = 2.0_DP * PI * e_lrc/box_list(this_box)%volume
    
  END SUBROUTINE Compute_LR_Correction
  
  !*****************************************************************************************************************

  SUBROUTINE Check_Interaction(im_1,is_1,im_2,is_2,get_interaction,rcom,rxcom,rycom,rzcom)
    
    REAL(DP) :: rxijp, ryijp, rzijp, rcom, rxcom, rycom, rzcom, rinteraction
    
    INTEGER :: this_box
    INTEGER :: im_1,is_1,im_2,is_2
    
    LOGICAL :: get_interaction
    
    ! Initially set the interaction to true.
    
    get_interaction = .TRUE.
    
    ! Figure out the box to be used later.
    
    this_box = molecule_list(im_1,is_1)%which_box
    
    IF(int_vdw_sum_style(this_box) == vdw_minimum) RETURN
    
    ! Parent separation
    
    rxijp = molecule_list(im_1,is_1)%xcom - molecule_list(im_2,is_2)%xcom
    ryijp = molecule_list(im_1,is_1)%ycom - molecule_list(im_2,is_2)%ycom
    rzijp = molecule_list(im_1,is_1)%zcom - molecule_list(im_2,is_2)%zcom
    
    ! Compute the minimum image distance
    
    CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp,rxcom,rycom,rzcom)
    
    rcom = DSQRT(rxcom*rxcom + rycom*rycom + rzcom*rzcom) 
    
    IF (CBMC_flag) THEN
       
       rinteraction = rcut_cbmc(this_box) + molecule_list(im_1,is_1)%max_dcom &
            + molecule_list(im_2,is_2)%max_dcom
       
       IF (rcom > rinteraction) get_interaction = .FALSE.
       
    ELSE 
       
       rinteraction = rcut_max(this_box) + molecule_list(im_1,is_1)%max_dcom &
            + molecule_list(im_2,is_2)%max_dcom
       
       IF (rcom > rinteraction) get_interaction = .FALSE.          
       
    END IF
    
 END SUBROUTINE Check_Interaction
!******************************************************************************************

!******************************************************************************************
 SUBROUTINE Energy_Test(rijsq,get_vdw,get_qq,this_box)
   
   INTEGER  :: this_box
   REAL(DP) :: rijsq,rcut_cbmcsq
   LOGICAL  :: get_vdw, get_qq
   
   rcut_cbmcsq = rcut_cbmc(this_box)*rcut_cbmc(this_box)
   
   get_vdw = .FALSE.
   get_qq = .FALSE.
   
   VDW_Test2: IF (int_vdw_style(this_box) == vdw_none) THEN
      get_vdw = .FALSE.
      
   ELSEIF (int_vdw_style(this_box) == vdw_lj) THEN
      
      IF (CBMC_flag) THEN
         IF (rijsq <= rcut_cbmcsq) THEN
            get_vdw = .TRUE.
         ELSE
            get_vdw = .FALSE.  
         ENDIF
      ELSEIF (int_vdw_sum_style(this_box) == vdw_cut .OR. int_vdw_sum_style(this_box) &
           == vdw_cut_shift .OR. int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
         
         IF (rijsq <= rcut_vdwsq(this_box)) THEN
            get_vdw = .TRUE.
         ELSE
            get_vdw = .FALSE.
         ENDIF
         
      ELSEIF (int_vdw_sum_style(this_box) == vdw_minimum) THEN
         get_vdw = .TRUE.
         
      ELSEIF (int_vdw_sum_style(this_box) == vdw_charmm) THEN
         get_vdw = .TRUE.

      ELSEIF (int_vdw_sum_style(this_box) == vdw_mie) THEN
           
         IF (rijsq <= rcut_vdwsq(this_box)) THEN 
            get_vdw = .TRUE.
         ELSE 
            get_vdw = .FALSE.
         ENDIF

      ELSEIF (int_vdw_sum_style(this_box) == vdw_cut_switch) THEN
         
         IF (rijsq <= roff_switch_sq(this_box)) THEN
            get_vdw = .TRUE.
         ELSE
            get_vdw = .FALSE.
         END IF
         
      ENDIF      
   
   
   
   ELSE
      err_msg = ""
      err_msg(1) = 'vdw_style must be NONE of LJ'
      CALL Clean_Abort(err_msg,'Compute_Atom_Nonbond_Energy')
      
   ENDIF VDW_Test2
   
   ! Charge sum tests
   IF (int_charge_style(this_box) == charge_none) THEN
      get_qq = .FALSE.
   ELSEIF (int_charge_style(this_box) == charge_coul) THEN
      
      IF (int_charge_sum_style(this_box) == charge_cut .OR. int_charge_sum_style(this_box) == charge_ewald) THEN
         IF(CBMC_flag) THEN
            IF (rijsq <= rcut_cbmcsq) THEN
               get_qq = .TRUE.
            ELSE
               get_qq = .FALSE.
            ENDIF
         ELSE
            IF (rijsq <= rcut_coulsq(this_box)) THEN
               get_qq = .TRUE.
            ELSE
               get_qq = .FALSE.
            ENDIF
         ENDIF
      ELSEIF (int_charge_sum_style(this_box) == charge_minimum) THEN
         get_qq = .TRUE.
         
      END IF
      
   ENDIF
   
   RETURN
   
 END SUBROUTINE Energy_Test
 
 SUBROUTINE Compute_Forces(this_box)
   
   !**************************************************************************************************
   ! The subroutine calculates the total forces of a given box. The identity of the box is passed to
   ! the routine. The forces are then used to compute the pressure tensor.
   !
   ! CALLS
   !
   ! CALLED BY
   !
   ! Volume_Change
   ! Main
   !
   !***************************************************************************************************

   IMPLICIT NONE
   
   INTEGER, INTENT(IN) :: this_box
   
   !----------------------------------------------------------------------------------------------------
   
   INTEGER ::  is, im_1, im_2, is_1, is_2, this_im_1, this_im_2
   
   REAL(DP) :: rcom, rx, ry, rz, w_lrc
   
   REAL(DP),DIMENSION(3,3) :: tv_pair, tc_pair, w_inter_vdw, w_inter_charge
   
   LOGICAL :: get_interaction
   
   W_tensor_vdw(:,:,this_box) = 0.0_DP
   W_tensor_charge(:,:,this_box) = 0.0_DP
   
   DO is = 1, nspecies
      imLOOP1: DO im_1 = 1, nmolecules(is)
         this_im_1 = locate(im_1,is)
         ! is it in this box?
         IF ( molecule_list(this_im_1,is)%which_box /= this_box ) CYCLE imLOOP1
         ! is alive?
         IF ( .NOT. molecule_list(this_im_1,is)%live ) CYCLE imLOOP1
         
         w_inter_vdw(:,:) = 0.0_DP
         w_inter_charge(:,:) = 0.0_DP
         
         !$OMP PARALLEL DO DEFAULT(SHARED) &
         !$OMP SCHEDULE(DYNAMIC) &
         !$OMP PRIVATE(im_2, this_im_2, get_interaction) &
         !$OMP PRIVATE(rcom, rx, ry, rz, tv_pair, tc_pair) &
         !$OMP REDUCTION(+:w_inter_vdw, w_inter_charge)
         imLOOP2: DO im_2 = im_1 + 1, nmolecules(is)
            this_im_2 = locate(im_2,is)
            ! allow interactions only in the box
            IF ( molecule_list(this_im_2,is)%which_box /= this_box ) CYCLE imLOOP2
            ! is it alive?
            IF ( .NOT. molecule_list(this_im_2,is)%live ) CYCLE imLOOP2
            
            CALL Check_Interaction(this_im_1,is,this_im_2,is,get_interaction,rcom,rx,ry,rz)
            
            IF (.NOT. Get_Interaction) CYCLE imLOOP2
            
            CALL Compute_Molecule_Pair_Force(this_im_1,is,this_im_2,is,this_box,tv_pair,tc_pair,rx,ry,rz)
            
            w_inter_vdw(:,:) = w_inter_vdw(:,:) + tv_pair(:,:)
            w_inter_charge(:,:) = w_inter_charge(:,:) + tc_pair(:,:)
            
         END DO imLOOP2
         !$OMP END PARALLEL DO 
         
         W_tensor_vdw(:,:,this_box) = W_tensor_vdw(:,:,this_box) + w_inter_vdw(:,:)
         W_tensor_charge(:,:,this_box) = W_tensor_charge(:,:,this_box) + w_inter_charge(:,:)
         
      END DO imLOOP1
   END DO
   
   DO is_1 = 1, nspecies
      imLOOP3: DO im_1 = 1, nmolecules(is_1)
         this_im_1 = locate(im_1,is_1)
         IF( molecule_list(this_im_1,is_1)%which_box /= this_box ) CYCLE imLOOP3
         IF( .NOT. molecule_list(this_im_1,is_1)%live ) CYCLE imLOOP3
         
         DO is_2 = is_1 + 1, nspecies
            
            w_inter_vdw(:,:) = 0.0_DP
            w_inter_charge(:,:) = 0.0_DP
            
            !$OMP PARALLEL DO DEFAULT(SHARED) &
            !$OMP SCHEDULE(DYNAMIC) &
            !$OMP PRIVATE(im_2, this_im_2, get_interaction) &
            !$OMP PRIVATE(rcom, rx, ry, rz, tv_pair, tc_pair) &
            !$OMP REDUCTION(+:w_inter_vdw,w_inter_charge)
            imLOOP4: DO im_2 = 1,nmolecules(is_2)
               this_im_2 = locate(im_2,is_2)
               IF ( molecule_list(this_im_2,is_2)%which_box /= this_box ) CYCLE imLOOP4
               IF ( .NOT. molecule_list(this_im_2,is_2)%live ) CYCLE imLOOP4
               
               ! Check to see if the interaction needs to be computed between the molecules
               CALL Check_Interaction(this_im_1,is_1,this_im_2,is_2,get_interaction,rcom,rx,ry,rz)
               
               IF (.NOT. get_interaction ) CYCLE imLOOP4
               
               CALL Compute_Molecule_Pair_Force(this_im_1,is_1,this_im_2,is_2,this_box,tv_pair,tc_pair,rx,ry,rz)
               
               !                W_tensor_vdw(:,:,this_box) = W_tensor_vdw(:,:,this_box) + tv_pair(:,:)
!                W_tensor_charge(:,:,this_box) = W_tensor_charge(:,:,this_box) + tc_pair(:,:)
               
               w_inter_vdw(:,:) = w_inter_vdw(:,:) + tv_pair(:,:)
               w_inter_charge(:,:) = w_inter_charge(:,:) + tc_pair(:,:)
               
            END DO imLOOP4
            
            W_tensor_vdw(:,:,this_box) = W_tensor_vdw(:,:,this_box) + w_inter_vdw(:,:)
            W_tensor_charge(:,:,this_box) = W_tensor_charge(:,:,this_box) + w_inter_charge(:,:)
         END DO
         
      END DO imLOOP3
    END DO
    
    IF (int_charge_sum_style(this_box) == charge_ewald) THEN
       
       CALL Compute_System_Ewald_Reciprocal_Force(this_box)
       
    END IF
    
    IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN

       CALL Compute_LR_Force(this_box,w_lrc)
       virial(this_box)%lrc = w_lrc
       
    END IF
    
    W_tensor_elec(:,:,this_box) = (W_tensor_charge(:,:,this_box) + W_tensor_recip(:,:,this_box))*charge_factor 
    W_tensor_total(:,:,this_box) = W_tensor_vdw(:,:,this_box) + W_tensor_elec(:,:,this_box) 
    
  END SUBROUTINE Compute_Forces
  
  SUBROUTINE Compute_Molecule_Pair_Force(im_1,is_1,im_2,is_2,this_box,tens_vdw,tens_charge,rabx,raby,rabz)
    !**********************************************************************************************************
    ! The subroutine returns the interaction force of the input molecule with another molecule. Thus,
    ! it computes the intermolecular vdw and electrostatic interactions. 
    !
    ! CALLED BY
    !
    ! Added by Jindal Shah on 12/10/07
    !**********************************************************************************************************
  
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: im_1, is_1, im_2, is_2, this_box
    !----------------------------------------------------------------------------------------------------------

    INTEGER :: ia, ja 

    REAL(DP) :: rxijp, ryijp, rzijp, rxij, ryij, rzij, rijsq, wij_vdw, wij_qq
    REAL(DP) :: rabx, raby, rabz
    REAL(DP),DIMENSION(3,3) :: tens_vdw, tens_charge

    REAL(DP) :: ffc, wxy, wxz, wyz

    LOGICAL :: get_vdw, get_qq

    tens_vdw(:,:) = 0.0_DP
    tens_charge(:,:) = 0.0_DP

    DO ia = 1, natoms(is_1)
       DO ja = 1, natoms(is_2)
          ! Obtain the minimum image separation

          rxijp = atom_list(ia,im_1,is_1)%rxp - atom_list(ja,im_2,is_2)%rxp
          ryijp = atom_list(ia,im_1,is_1)%ryp - atom_list(ja,im_2,is_2)%ryp
          rzijp = atom_list(ia,im_1,is_1)%rzp - atom_list(ja,im_2,is_2)%rzp
          
          ! Now get the minimum image separation 
          CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp,rxij,ryij,rzij)

          rijsq = rxij*rxij + ryij*ryij + rzij*rzij

          ! Now figure out what needs to be computed, then call pair_energy
          
          CALL Energy_Test(rijsq,get_vdw,get_qq,this_box)          

          ! Compute vdw and q-q energy using if required
          IF (get_vdw .OR. get_qq) THEN 

             CALL Pair_Force(rijsq,is_1,im_1,ia,is_2,im_2,ja,&
                  get_vdw,get_qq,Wij_vdw,Wij_qq)

             ffc = Wij_vdw/rijsq

             wxy = ffc*(0.5_DP*(rxij*raby+ryij*rabx))
             wxz = ffc*(0.5_DP*(rxij*rabz+rzij*rabx))
             wyz = ffc*(0.5_DP*(ryij*rabz+rzij*raby))

             tens_vdw(1,1) = tens_vdw(1,1) + ffc*rxij*rabx
             tens_vdw(1,2) = tens_vdw(1,2) + wxy
             tens_vdw(1,3) = tens_vdw(1,3) + wxz
             tens_vdw(2,1) = tens_vdw(2,1) + wxy
             tens_vdw(2,2) = tens_vdw(2,2) + ffc*ryij*raby
             tens_vdw(2,3) = tens_vdw(2,3) + wyz        
             tens_vdw(3,1) = tens_vdw(3,1) + wxz
             tens_vdw(3,2) = tens_vdw(3,2) + wyz
             tens_vdw(3,3) = tens_vdw(3,3) + ffc*rzij*rabz

             ffc = Wij_qq/rijsq

             wxy = ffc*(0.5_DP*(rxij*raby+ryij*rabx))
             wxz = ffc*(0.5_DP*(rxij*rabz+rzij*rabx))
             wyz = ffc*(0.5_DP*(ryij*rabz+rzij*raby))

             tens_charge(1,1) = tens_charge(1,1) + ffc*rxij*rabx
             tens_charge(1,2) = tens_charge(1,2) + wxy
             tens_charge(1,3) = tens_charge(1,3) + wxz
             tens_charge(2,1) = tens_charge(2,1) + wxy
             tens_charge(2,2) = tens_charge(2,2) + ffc*ryij*raby
             tens_charge(2,3) = tens_charge(2,3) + wyz        
             tens_charge(3,1) = tens_charge(3,1) + wxz
             tens_charge(3,2) = tens_charge(3,2) + wyz
             tens_charge(3,3) = tens_charge(3,3) + ffc*rzij*rabz

          END IF

       END DO

    END DO

  END SUBROUTINE Compute_Molecule_Pair_Force
  !------------------------------------------------------------------------------------------
  SUBROUTINE Pair_Force &
       (rijsq,is,im,ia,js,jm,ja,get_vdw,get_qq,Wij_vdw,Wij_qq)

    ! LJ potential:  Wij = -rij/3 * d Eij / d rij. Use the virial in: P = NkBT + < W >

    ! Computes the vdw and q-q pair force between atoms ia and ja of molecules im and jm
    ! and species is and js, given their separation rijsq. I have passed each component of 
    ! separation nut right now this is unnecessary. 
    ! It also computes the real space part of the Ewald sum if necessary.

    ! Called by: 
    ! Calls: CFC_LJ_Scaling
    !        Ewald_Real
  !------------------------------------------------------------------------------------------
    ! Passed to
    REAL(DP) :: rxij,ryij,rzij,rijsq
    INTEGER :: is,im,ia,js,jm,ja,ibox
    LOGICAL :: get_vdw,get_qq, fraction

    ! Returned
    REAL(DP) :: Wij_vdw,Wij_qq,Eij_vdw

    ! Local
    INTEGER :: itype,jtype, this_box
    REAL(DP) :: eps,sig,SigOverRsq,SigOverR6,SigOverR12
    REAL(DP) :: SigOverRsq_shift,SigOverR6_shift,SigOverR12_shift
    REAL(DP) :: SigOverR, SigOverRn, SigOverRm, mie_coeff, mie_n, mie_m
    REAL(DP) :: roffsq_rijsq, roffsq_rijsq_sq, factor2, fscale
    REAL(DP) :: qi,qj, qsc, erf_val, this_lambda
    REAL(DP) :: rij, ewald_constant, exp_const, Wij_self

  !------------------------------------------------------------------------------------------

    Wij_vdw = 0.0_DP
    Wij_qq = 0.0_DP
    fraction = .false.
    this_lambda = 1.0_DP

    ibox = molecule_list(im,is)%which_box

    ! If either atom is not yet present, then don't try to compute an energy
    ExistCheck: IF (atom_list(ia,im,is)%exist .AND. atom_list(ja,jm,js)%exist) THEN

       ! Determine atom type indices
       itype = nonbond_list(ia,is)%atom_type_number
       jtype = nonbond_list(ja,js)%atom_type_number
       
       VDW_calculation: IF (get_vdw) THEN

          LJ_12_6_calculation: IF (int_vdw_style(1) == vdw_lj) THEN
             ! For now, assume all interactions are the same. Use the lookup table created in Compute_Nonbond_Table
             eps = vdw_param1_table(itype,jtype)
             sig = vdw_param2_table(itype,jtype)

             ! Apply intramolecular scaling if necessary
             IF (is == js .AND. im == jm) THEN
                
                ! This controls 1-2, 1-3, and 1-4 interactions
                
                eps = eps * vdw_intra_scale(ia,ja,is)

             ENDIF

             IF (int_vdw_sum_style(ibox) == vdw_charmm) THEN
                SigOverRsq = (sig**2.0_DP)/rijsq
                SigOverR6  = SigOverRsq * SigOverRsq * SigOverRsq
                SigOverR12 = SigOverR6 * SigOverR6

                Wij_vdw = (12.0_DP * eps ) * (SigOverR12 - SigOverR6)

             ELSEIF (int_vdw_sum_style(ibox) == vdw_cut .OR. int_vdw_sum_style(ibox) == vdw_cut_tail) THEN
                SigOverRsq = (sig**2)/rijsq
                SigOverR6 = SigOverRsq * SigOverRsq * SigOverRsq
                SigOverR12 = SigOverR6 * SigOverR6
                
                Wij_vdw = (24.0_DP * eps) * (2.0_DP*SigOverR12 - SigOverR6)
             ELSEIF (int_vdw_sum_style(ibox) == vdw_cut_shift) THEN
                SigOverRsq = (sig**2)/rijsq
                SigOverR6 = SigOverRsq * SigOverRsq * SigOverRsq
                SigOverR12 = SigOverR6 * SigOverR6

                SigOverRsq_shift = sig**2/rcut_vdwsq(ibox)
                SigOverR6_shift = SigOverRsq_shift * SigOverRsq_shift * SigOverRsq_shift
                SigOverR12_shift = SigOverR6_shift * SigOverR6_shift

                Wij_vdw = (24.0_DP * eps) * (2.0_DP*SigOverR12 - SigOverR6) 

             ELSEIF (int_vdw_sum_style(ibox) == vdw_cut_switch) THEN
                
                SigOverRsq = (sig**2)/rijsq
                SigOverR6 = SigOverRsq * SigOverRsq * SigOverRsq
                SigOverR12 = SigOverR6 * SigOverR6
                
                IF ( (rijsq < ron_switch_sq(ibox) )) THEN
                   
                   Wij_vdw = (24.0_DP * eps) * (2.0_DP*SigOverR12 - SigOverR6)
                   fscale = 1.0_DP
                   
                ELSE IF ( rijsq <= roff_switch_sq(ibox)) THEN
                   
                   roffsq_rijsq = roff_switch_sq(ibox) - rijsq
                  
                   roffsq_rijsq_sq = roffsq_rijsq * roffsq_rijsq

                   factor2 = switch_factor2(ibox) + 2.0_DP * rijsq

                   fscale = roffsq_rijsq_sq * factor2 * switch_factor1(ibox)

                   Eij_vdw = fscale * Eij_vdw
                   Wij_vdw = fscale * ( 24.0_DP * eps / 3.0_DP ) * (2.0_DP*SigOverR12 - SigOverR6)
                   Wij_vdw = Wij_vdw + &
                        (8.0_DP * rijsq * rijsq * roffsq_rijsq * Eij_vdw * switch_factor1(ibox))/(3.0_DP)

                END IF
             ELSEIF (int_vdw_sum_style(ibox) == vdw_mie) THEN
                rij = SQRT(rijsq)

                mie_n = mie_nlist(mie_Matrix(is,js))
                mie_m = mie_mlist(mie_Matrix(is,js))
                mie_coeff = mie_n/(mie_n-mie_m) * (mie_n/mie_m)**(mie_m/(mie_n-mie_m))
                SigOverR = sig/rij
                SigOverRn = SigOverR ** mie_n
                SigOverRm = SigOverR ** mie_m
                Wij_vdw = (mie_coeff * eps) *(mie_n * SigOverRn - mie_m * SigOverRm)


             ELSE

                fscale = 0.0_DP
                Eij_vdw = 0.0_DP
                Wij_vdw = 0.0_DP
                
                
                               
             ENDIF
             
             ! Add other potential types here
          ENDIF LJ_12_6_calculation
          
       ENDIF VDW_calculation
       qq_calculation: IF (get_qq) THEN

          qi = nonbond_list(ia,is)%charge
          qj = nonbond_list(ja,js)%charge


          IF (int_charge_sum_style(ibox) == charge_cut) THEN
             ! Apply charge scaling for intramolecular energies
             qsc = 1.0_DP
             IF ( is == js .AND. im == jm ) THEN
                qsc = charge_intra_scale(ia,ja,is)
             END IF
             Wij_qq = qsc*charge_factor*(qi*qj)/SQRT(rijsq)
          ELSEIF (int_charge_sum_style(ibox) == charge_ewald) THEN
             ! Real space Ewald part
             this_box = molecule_list(im,is)%which_box
             qsc = 1.0_DP
             ibox = molecule_list(im,is)%which_box 

             ! Apply intramolecular scaling if necessary
             IF (is == js .AND. im == jm) THEN
       
                ! Intramolecular charge scaling
                qsc = charge_intra_scale(ia,ja,is)

             ENDIF

             ! Real space part: This does the intrascaling correct. For cfc intra,
             ! we use full scaling. I think for CFC inter, we simply scale the actual 
             ! value of the charge, but do NOT scale it here for intra interactions.
             ! Come back to this later.

             rij = SQRT(rijsq)
             ewald_constant = 2.0_DP * alpha_ewald(ibox) / rootPI
             exp_const = EXP(-alpha_ewald(ibox)*alpha_ewald(ibox)*rijsq) 
             ! May need to protect against very smamie_coeffsq
             erf_val = 1.0_DP - erfc(alpha_ewald(ibox) * rij)
             Wij_qq = qi*qj*( (qsc - erf_val)/rij + ewald_constant*exp_const )

             IF (is == js .AND. im == jm) THEN

                Wij_self = (qsc - 1.0_DP) * qi*qj * (erf_val/rij - ewald_constant * exp_const)
                Wij_qq = Wij_qq + Wij_self

             END IF

             ! self and recipricoal parts need to be computed as total energy differences between original
             ! configuration and the perturbed configuration. These terms are thus added on after all atoms 
             ! have been moved. 

          ENDIF

       ENDIF qq_calculation

    ENDIF ExistCheck
!------------------------------------------------------------------------------
  CONTAINS

    FUNCTION erfc(x)
      !**************************************************************************
      !                                                                         *
      ! Calculate the complementary error function for  a number
      !                                                                         *
      !**************************************************************************

      REAL(DP) :: erfc
      REAL(DP), PARAMETER :: A1 = 0.254829592_DP, A2 = -0.284496736_DP
      REAL(DP), PARAMETER :: A3 = 1.421413741_DP, A4 = -1.453152027_DP
      REAL(DP), PARAMETER :: A5 = 1.061405429_DP, P = 0.3275911_DP
      REAL(DP) :: T, x, xsq, TP

      T = 1.0_DP / (1.0_DP + P*x)
      xsq = x*x

      TP = T * (A1 + T * (A2 + T * (A3 + T * (A4 + T * A5))))

      erfc = TP * EXP(-xsq)

    END FUNCTION erfc
!------------------------------------------------------------------------------

  END SUBROUTINE Pair_Force

  !----------------------------------------------------------------------------------------------
  SUBROUTINE Compute_LR_Force(this_box, w_lrc)
    !*************************************************************************************************************
    ! The subroutine calculates the long range correction for the given box.
    !
    ! Called by
    ! 
    ! First written by Jindal Shah on 01/10/08
    !
    !
    !**************************************************************************************************************
    
    INTEGER, INTENT(IN) :: this_box
    REAL(DP), INTENT(OUT) :: w_lrc
    
    INTEGER ::   ia, ja

    REAL(DP) :: epsij, sigij
    REAL(DP) :: w_lrc_ia_ja

    w_lrc = 0.0_DP 

    DO ia = 1, nbr_atomtypes

       w_lrc_ia_ja = 0.0_DP
          
       DO ja = 1, nbr_atomtypes

          epsij = vdw_param1_table(ia,ja)
          sigij = vdw_param2_table(ia,ja)

          w_lrc_ia_ja = w_lrc_ia_ja + nint_beads(ja,this_box) * epsij * ((2.0_DP / 3.0_DP * &
                        sigij**12 / rcut9(this_box)) - (sigij**6 / rcut3(this_box)))
             
       END DO

       w_lrc = w_lrc + nint_beads(ia,this_box) * w_lrc_ia_ja

    END DO

    w_lrc = 16.0_DP / 3.0_DP * PI * w_lrc / box_list(this_box)%volume

  END SUBROUTINE Compute_LR_Force

  SUBROUTINE Compute_System_Ewald_Reciprocal_Force(this_box)
    !*****************************************************************************************
    ! This subroutine computes the long range forces due to electrostatics
    !
    ! Based on APSS code reciprocal_ewald.f90
    !
    ! Added by Tom Rosch on 06/11/09 (See Wheeler, Mol. Phys. 1997 Vol. 92 pg. 55)
    !
    !*****************************************************************************************
    
    USE Type_Definitions
    USE Run_Variables

    IMPLICIT NONE

!    !$ include 'omp_lib.h'

    INTEGER :: i, is, im, ia, this_locate, this_box

    REAL(DP) :: charge
    REAL(DP) :: qw(9), qwxy, qwxz, qwyz, un, const_val
    REAL(DP) :: xcmi, ycmi, zcmi, piix, piiy, piiz, arg, factor
    REAL(DP) :: recip_11, recip_21, recip_31, recip_22, recip_23, recip_33

    const_val = 1.0_DP/(2.0_DP * alpha_ewald(this_box) * alpha_ewald(this_box))
    qw(:) = 0.0_DP
    W_tensor_recip(:,:,this_box) = 0.0_DP

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP SCHEDULE(STATIC) &
    !$OMP PRIVATE(i, un, qwxy, qwxz, qwyz) &
    !$OMP REDUCTION(+:qw)
    DO i = 1, nvecs(this_box)

       un = Cn(i,this_box) * (cos_sum(i,this_box) * cos_sum(i,this_box) + sin_sum(i,this_box) * sin_sum(i,this_box))

       qwxy =  un * ( -2.0_DP*(1.0_DP/hsq(i,this_box) + 0.5_DP*const_val) &
              *hx(i,this_box)*hy(i,this_box) )
       qwxz =  un * ( -2.0_DP*(1.0_DP/hsq(i,this_box) + 0.5_DP*const_val) &
              *hx(i,this_box)*hz(i,this_box) )
       qwyz =  un * ( -2.0_DP*(1.0_DP/hsq(i,this_box) + 0.5_DP*const_val) &
              *hy(i,this_box)*hz(i,this_box) )

       qw(1) = qw(1) + &
               ( un * ( 1.0_DP - 2.0_DP*(1.0_DP/hsq(i,this_box) + 0.5_DP*const_val) &
               *hx(i,this_box)*hx(i,this_box)))
       qw(2) = qw(2) + qwxy
       qw(3) = qw(3) + qwxz
       qw(5) = qw(5) + &
               ( un * ( 1.0_DP - 2.0_DP*(1.0_DP/hsq(i,this_box) + 0.5_DP*const_val) &
               *hy(i,this_box)*hy(i,this_box)))
       qw(6) = qw(6) + qwyz
       qw(9) = qw(9) + &
               ( un * ( 1.0_DP - 2.0_DP*(1.0_DP/hsq(i,this_box) + 0.5_DP*const_val) &
               *hz(i,this_box)*hz(i,this_box)))

    END DO 
    !$OMP END PARALLEL DO

    W_tensor_recip(1,1,this_box) = qw(1)
    W_tensor_recip(2,1,this_box) = qw(2)
    W_tensor_recip(3,1,this_box) = qw(3)
    W_tensor_recip(2,2,this_box) = qw(5)
    W_tensor_recip(3,2,this_box) = qw(6)
    W_tensor_recip(3,3,this_box) = qw(9)

    DO is = 1, nspecies

       DO im = 1, nmolecules(is)

          this_locate = locate(im,is)

          IF( .NOT. molecule_list(this_locate,is)%live) CYCLE
          IF( molecule_list(this_locate,is)%which_box /= this_box ) CYCLE

          xcmi = molecule_list(this_locate,is)%xcom
          ycmi = molecule_list(this_locate,is)%ycom
          zcmi = molecule_list(this_locate,is)%zcom

          DO ia = 1, natoms(is)

             piix = atom_list(ia,this_locate,is)%rxp - xcmi
             piiy = atom_list(ia,this_locate,is)%ryp - ycmi
             piiz = atom_list(ia,this_locate,is)%rzp - zcmi
             charge = nonbond_list(ia,is)%charge

             recip_11 = 0.0_DP
             recip_21 = 0.0_DP
             recip_31 = 0.0_DP
             recip_22 = 0.0_DP
             recip_23 = 0.0_DP
             recip_33 = 0.0_DP

             !$OMP PARALLEL DO DEFAULT(SHARED) &
             !$OMP SCHEDULE(STATIC) &
             !$OMP PRIVATE(i,arg,factor) &
             !$OMP REDUCTION(+:recip_11, recip_21, recip_31) &
             !$OMP REDUCTION(+:recip_22, recip_23, recip_33) 
             DO i = 1, nvecs(this_box)

                arg = hx(i,this_box)*atom_list(ia,this_locate,is)%rxp + &
                      hy(i,this_box)*atom_list(ia,this_locate,is)%ryp + &
                      hz(i,this_box)*atom_list(ia,this_locate,is)%rzp

                factor = Cn(i,this_box)*2.0_DP*(-cos_sum(i,this_box)*DSIN(arg) + &
                         sin_sum(i,this_box)*DCOS(arg))*charge

                recip_11 = recip_11 + factor*hx(i,this_box)*piix
                recip_21 = recip_21 + factor* 0.5_DP*(hx(i,this_box)*piiy+hy(i,this_box)*piix)
                recip_31 = recip_31 + factor* 0.5_DP*(hx(i,this_box)*piiz+hz(i,this_box)*piix)
                recip_22 = recip_22 + factor*hy(i,this_box)*piiy
                recip_23 = recip_23 + factor* 0.5_DP*(hy(i,this_box)*piiz+hz(i,this_box)*piiy)
                recip_33 = recip_33 + factor*hz(i,this_box)*piiz

             END DO
             !$OMP END PARALLEL DO

             W_tensor_recip(1,1,this_box) = W_tensor_recip(1,1,this_box) + recip_11
             W_tensor_recip(2,1,this_box) = W_tensor_recip(2,1,this_box) + recip_21
             W_tensor_recip(3,1,this_box) = W_tensor_recip(3,1,this_box) + recip_31
             W_tensor_recip(2,2,this_box) = W_tensor_recip(2,2,this_box) + recip_22
             W_tensor_recip(2,3,this_box) = W_tensor_recip(2,3,this_box) + recip_23
             W_tensor_recip(3,3,this_box) = W_tensor_recip(3,3,this_box) + recip_33
             

          END DO

       END DO

    END DO   

  W_tensor_recip(1,2,this_box) = W_tensor_recip(2,1,this_box)
  W_tensor_recip(1,3,this_box) = W_tensor_recip(3,1,this_box)
  W_tensor_recip(3,2,this_box) = W_tensor_recip(2,3,this_box)

  END SUBROUTINE Compute_System_Ewald_Reciprocal_Force

  !***********************************************************
  SUBROUTINE Compute_Molecule_Energy(im,is,f_intra,this_nrg,overlap,f_inter, E_self)
    !****************************************************************
    !
    ! This routine computes the energy of a molecule to be used in
    ! various moves. 
    !
    ! First written by Jindal Shah on 01/30/09.
    !
    ! 10/11/10 (JS) : optional logical flag if intermolecular interactions
    !                 are to be computed. 
    !
    !*****************************************************************
    
    INTEGER, INTENT(IN) :: im,is
    LOGICAL, INTENT(IN) :: f_intra
    REAL(DP), INTENT(OUT), OPTIONAL :: E_self
    LOGICAL, INTENT(IN), OPTIONAL :: f_inter
    
    Type(Energy_Class), INTENT(OUT) :: this_nrg
    
    LOGICAL, INTENT(OUT) :: overlap
    LOGICAL :: intra_overlap, l_inter

    this_nrg%total = 0.0_DP
    this_nrg%bond = 0.0_DP
    this_nrg%angle = 0.0_DP
    this_nrg%dihedral = 0.0_DP
    this_nrg%intra = 0.0_DP
    this_nrg%intra_vdw = 0.0_DP
    this_nrg%intra_q = 0.0_DP
    this_nrg%inter_vdw = 0.0_DP
    this_nrg%inter_q = 0.0_DP
    this_nrg%lrc = 0.0_DP
    this_nrg%ewald_reciprocal = 0.0_DP
    this_nrg%ewald_self = 0.0_DP
    IF(PRESENT(E_self)) E_self = 0.0_DP
    
    overlap = .FALSE.

    l_inter = .TRUE.

    IF (present(f_inter)) THEN

       IF ( .NOT. f_inter) l_inter = .FALSE.

    END IF
    

    IF (l_inter) THEN
       CALL Compute_Molecule_Nonbond_Inter_Energy(im,is,this_nrg%inter_vdw,this_nrg%inter_q,overlap)
    END IF
    
    IF(overlap) RETURN
    
    IF (f_intra) THEN
       ! compute intramolecular bonded energy for this molecule
       CALL Compute_Molecule_Bond_Energy(im,is,this_nrg%bond)
       CALL Compute_Molecule_Angle_Energy(im,is,this_nrg%angle)
       CALL Compute_Molecule_Dihedral_Energy(im,is,this_nrg%dihedral)
       CALL Compute_Molecule_Nonbond_Intra_Energy(im,is,this_nrg%intra_vdw,this_nrg%intra_q,intra_overlap,E_self)
       
       this_nrg%intra = this_nrg%bond + this_nrg%angle + this_nrg%dihedral
    END IF
       
    
  END SUBROUTINE Compute_Molecule_Energy

  SUBROUTINE Compute_Ring_Fragment_Energy(this_frag,this_im,is,this_box,nrg_ring_frag)
    !******************************************************************************************
    !
    ! This subroutine calculates the energy of a ring fragment in its old conformation
    ! 
    ! CALLED BY: 
    !           fragment_growth.f90
    !
    ! CALLS :
    !       Compute_Molecule_Dihedral_Energy
    !       Compute_Molecule_Nonbond_Intra_Energy
    !
    ! Written by Jindal Shah on 10/03/09
    !
    !******************************************************************************************

    INTEGER, INTENT(IN) :: this_frag, this_im, is
    REAL(DP), INTENT(OUT) :: nrg_ring_frag
    INTEGER :: this_box

    ! local variables

    INTEGER :: i, this_atom

!    REAL(DP) :: rcut_vdwsq_box, rcut_coulsq_box, alpha_ewald_box
    REAL(DP) :: e_dihed,  e_improper, nrg_vdw, nrg_qq

    LOGICAL :: intra_overlap
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: exist_flag_old
    
    nrg_ring_frag = 0.0_DP

    !!! Note for now, keep this_box == 1. For flexible ring
    ! molecule change this

    this_box = molecule_list(this_im,is)%which_box


    ! first store the exist flag of the molecule

    ALLOCATE(exist_flag_old(natoms(is)))
    exist_flag_old = atom_list(1:natoms(is),this_im,is)%exist

    atom_list(1:natoms(is),this_im,is)%exist = .FALSE.
    ! Now turn on the exist flag of the ring fragment

    DO i = 1,frag_list(this_frag,is)%natoms

       this_atom = frag_list(this_frag,is)%atoms(i)
       atom_list(this_atom,this_im,is)%exist = .TRUE.

    END DO

    ! Now set the cutoffs for VDW and charge interactions along with
    ! the Ewald parameters. For now switch off the LJ and electrostatic
    ! interactions.

    ! store old cutoffs

!!$    rcut_vdwsq_box = rcut_vdwsq(this_box)
!!$
!!$    rcut_vdwsq(this_box) = frag_list(this_frag,is)%rcut_vdwsq
!!$
!!$    IF (int_charge_sum_style(this_box) == charge_ewald ) THEN
!!$
!!$       rcut_coulsq_box = rcut_coulsq(this_box)
!!$       alpha_ewald_box = alpha_ewald(this_box)
!!$
!!$       rcut_coulsq(this_box) = frag_list(this_frag,is)%rcut_coulsq
!!$       alpha_ewald(this_box) = frag_list(this_frag,is)%alpha_ewald
!!$
!!$    END IF

    ! Now compute the intramolecular energy of the fragment

    CALL Compute_Molecule_Dihedral_Energy(this_im,is,e_dihed)
    CALL Compute_Molecule_Improper_Energy(this_im,is,e_improper)
    CALL Compute_Molecule_Nonbond_Intra_Energy(this_im,is,nrg_vdw,nrg_qq,intra_overlap)


    nrg_ring_frag = nrg_vdw + nrg_qq + e_dihed + e_improper

    ! Now reset all the cutoffs back
!!$
!!$    IF (int_charge_sum_style(this_box) == charge_ewald) THEN
!!$
!!$       alpha_ewald(this_box) = alpha_ewald_box
!!$       rcut_coulsq(this_box) = rcut_coulsq_box
!!$
!!$    END IF
!!$
!!$    rcut_vdwsq(this_box) = rcut_vdwsq_box

    ! Turn the original exist flag for the atoms on

    atom_list(1:natoms(is),this_im,is)%exist = exist_flag_old

    DEALLOCATE(exist_flag_old)

  END SUBROUTINE Compute_Ring_Fragment_Energy

  SUBROUTINE System_Energy_Check(this_box,i_step,randno)

     USE Run_Variables

     INTEGER, INTENT(IN) :: this_box, i_step
     REAL(DP), INTENT(IN) :: randno

     LOGICAL :: inter_overlap
     LOGICAL :: aok

     TYPE(Energy_Class) :: e_check
     TYPE(Energy_Class) :: e_diff 

     aok = .TRUE.

     e_check%total = energy(this_box)%total
     e_check%bond = energy(this_box)%bond
     e_check%angle = energy(this_box)%angle
     e_check%dihedral = energy(this_box)%dihedral
     e_check%improper = energy(this_box)%improper
     e_check%intra_vdw = energy(this_box)%intra_vdw
     e_check%intra_q = energy(this_box)%intra_q
     e_check%inter_vdw = energy(this_box)%inter_vdw
     e_check%inter_q = energy(this_box)%inter_q
     e_check%lrc = energy(this_box)%lrc
     e_check%ewald_reciprocal = energy(this_box)%ewald_reciprocal
     e_check%ewald_self = energy(this_box)%ewald_self

     CALL Compute_Total_System_Energy(this_box,.TRUE.,inter_overlap)

     e_diff%total = ABS(energy(this_box)%total - e_check%total)
     e_diff%bond = ABS(energy(this_box)%bond - e_check%bond)
     e_diff%angle = ABS(energy(this_box)%angle - e_check%angle)
     e_diff%dihedral = ABS(energy(this_box)%dihedral - e_check%dihedral)
     e_diff%improper = ABS(energy(this_box)%improper - e_check%improper)
     e_diff%intra_vdw = ABS(energy(this_box)%intra_vdw - e_check%intra_vdw)
     e_diff%intra_q = ABS(energy(this_box)%intra_q - e_check%intra_q)
     e_diff%inter_vdw = ABS(energy(this_box)%inter_vdw - e_check%inter_vdw)
     e_diff%inter_q = ABS(energy(this_box)%inter_q - e_check%inter_q)
     e_diff%lrc = ABS(energy(this_box)%lrc - e_check%lrc)
     e_diff%ewald_reciprocal = ABS(energy(this_box)%ewald_reciprocal - e_check%ewald_reciprocal)
     e_diff%ewald_self = ABS(energy(this_box)%ewald_self - e_check%ewald_self)

     IF(e_diff%total .GT. 0.000001) THEN
        WRITE(logunit,*) 'Total energy does not match', 'box: ', this_box
        WRITE(logunit,*) 'Cumulative',e_check%total,'Total',energy(this_box)%total
        aok = .FALSE.
     END IF

     IF(e_diff%bond .GT. 0.000001) THEN
        WRITE(logunit,*) 'Bond energy does not match', 'box: ', this_box
        WRITE(logunit,*) 'Cumulative',e_check%bond,'Total',energy(this_box)%bond
        aok = .FALSE.
     END IF

     IF(e_diff%angle .GT. 0.000001) THEN
        WRITE(logunit,*) 'Angle energy does not match', 'box: ', this_box
        WRITE(logunit,*) 'Cumulative',e_check%angle,'Total',energy(this_box)%angle
        aok = .FALSE.
     END IF

     IF(e_diff%dihedral .GT. 0.00001) THEN
        WRITE(logunit,*) 'Dihedral energy does not match', 'box: ', this_box
        WRITE(logunit,*) 'Cumulative',e_check%dihedral,'Total',energy(this_box)%dihedral
        aok = .FALSE.
     END IF

     IF(e_diff%improper .GT. 0.01) THEN
        WRITE(logunit,*) 'Improper energy does not match', 'box: ', this_box
        WRITE(logunit,*) 'Cumulative',e_check%improper,'Total',energy(this_box)%improper
        aok = .FALSE.
     END IF

     IF(e_diff%intra_vdw .GT. 0.000001) THEN
        WRITE(logunit,*) 'Intra_vdw energy does not match', 'box: ', this_box
        WRITE(logunit,*) 'Cumulative',e_check%intra_vdw,'Total',energy(this_box)%intra_vdw
        aok = .FALSE.
     END IF

     IF(e_diff%intra_q .GT. 0.00001) THEN
        WRITE(logunit,*) 'Intra_q energy does not match', 'box: ', this_box
        WRITE(logunit,*) 'Cumulative',e_check%intra_q,'Total',energy(this_box)%intra_q
        aok = .FALSE.
     END IF

     IF(e_diff%inter_vdw .GT. 0.000001) THEN
        WRITE(logunit,*) 'Inter_vdw energy does not match', 'box: ', this_box
        WRITE(logunit,*) 'Cumulative',e_check%inter_vdw,'Total',energy(this_box)%inter_vdw
        aok = .FALSE.
     END IF

     IF(e_diff%inter_q .GT. 0.000001) THEN
        WRITE(logunit,*) 'Inter_vdw energy does not match', 'box: ', this_box
        WRITE(logunit,*) 'Cumulative',e_check%inter_q,'Total',energy(this_box)%inter_q
        aok = .FALSE.
     END IF
 
     IF(e_diff%lrc .GT. 0.0000001) THEN
        WRITE(logunit,*) 'LRC energy does not match', 'box: ', this_box
        WRITE(logunit,*) 'Cumulative',e_check%lrc,'Total',energy(this_box)%lrc
        aok = .FALSE.
     END IF

     IF(e_diff%ewald_reciprocal .GT. 0.000001) THEN
        WRITE(logunit,*) 'Reciprocal energy does not match', 'box: ', this_box
        WRITE(logunit,*) 'Cumulative',e_check%ewald_reciprocal,'Total',energy(this_box)%ewald_reciprocal
        aok = .FALSE.
     END IF

     IF(e_diff%ewald_self .GT. 0.000001) THEN
        WRITE(logunit,*) 'Self energy does not match', 'box: ', this_box
        WRITE(logunit,*) 'Cumulative',e_check%ewald_self,'Total',energy(this_box)%ewald_self
        aok = .FALSE.
     END IF

     IF(.NOT. aok) THEN
        IF(randno <= cut_trans) THEN
           WRITE(logunit,*) 'Problem after translation on step', i_step
        ELSE IF(randno <= cut_rot) THEN
           WRITE(logunit,*) 'Problem after rotation on step', i_step
        ELSE IF(randno <= cut_torsion) THEN
           WRITE(logunit,*) 'Problem after dihedral on step', i_step
        ELSE IF(randno <= cut_volume) THEN
           WRITE(logunit,*) 'Problem after volume on step', i_step
        ELSE IF(randno <= cut_angle) THEN
           WRITE(logunit,*) 'Problem after angle on step', i_step
        ELSE IF(randno <= cut_insertion) THEN
           WRITE(logunit,*) 'Problem after insertion on step', i_step
        ELSE IF(randno <= cut_deletion) THEN
           WRITE(logunit,*) 'Problem after deletion on step', i_step
        ELSE IF(randno <= cut_swap) THEN
           WRITE(logunit,*) 'Problem after swap on step', i_step
        ELSE IF(randno <= cut_regrowth) THEN
           WRITE(logunit,*) 'Problem after regrowth on step', i_step
        ELSE IF(randno <= cut_atom_displacement) THEN
           WRITE(logunit,*) 'Problem after atom displacement on step', i_step
        END IF
        STOP
      ELSE
         energy(this_box)%bond = e_check%bond
         energy(this_box)%angle = e_check%angle
         energy(this_box)%dihedral = e_check%dihedral
         energy(this_box)%improper = e_check%improper
         energy(this_box)%intra_vdw = e_check%intra_vdw
         energy(this_box)%intra_q = e_check%intra_q
         energy(this_box)%inter_vdw = e_check%inter_vdw
         energy(this_box)%inter_q = e_check%inter_q
         energy(this_box)%lrc = e_check%lrc
         energy(this_box)%ewald_reciprocal = e_check%ewald_reciprocal
         energy(this_box)%ewald_self = e_check%ewald_self
      END IF 

   END SUBROUTINE System_Energy_Check

   SUBROUTINE Compute_System_Ewald_Reciprocal_Energy2(this_box)
    !*****************************************************************************************
    ! This subroutine computes the sin and cos sum terms for the calculation of reciprocal
    ! energy of the input box. 
    !
    !*************************************************************************
    !*****************************************************************************************
    
    USE Type_Definitions
    USE Run_Variables

    IMPLICIT NONE

!    !$ include 'omp_lib.h'

    INTEGER :: i, is, im, ia, this_locate, this_box

    REAL(DP) :: un, const_val, this_scaling
    REAL(DP) :: charge, hdotr, energy_temp

    ! individual k-space vector stuff

    INTEGER ::  position
    INTEGER, ALLOCATABLE :: im_locate(:,:)

    ! openmp stuff

!    INTEGER :: omp_get_num_threads, omp_get_thread_num

    
    ! loop over all the k vectors of this box
          
    const_val = 1.0_DP/(2.0_DP * alpha_ewald(this_box) * alpha_ewald(this_box))
    energy(this_box)%ewald_reciprocal = 0.0_DP

    ALLOCATE(im_locate(MAXVAL(nmolecules),nspecies))
    
    
    DO is = 1, nspecies
       
       DO im = 1, nmolecules(is)
          
          this_locate = locate(im,is)
          
          IF (.NOT. molecule_list(this_locate,is)%live) CYCLE

          IF (molecule_list(this_locate,is)%which_box /= this_box) CYCLE

          IF (is == 1) THEN
             im_locate(im,is) = this_locate
          ELSE
             im_locate(im,is) = SUM(nmolecules(1:is-1)) + this_locate
          END IF

       END DO

    END DO   
    
    energy_temp = 0.0_DP
    !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
    cos_sum(:,this_box) = 0.0_DP
    sin_sum(:,this_box) = 0.0_DP
    !$OMP END PARALLEL WORKSHARE

    DO is = 1, nspecies

       IF ( .NOT. has_charge(is)) CYCLE
       
       DO im = 1, nmolecules(is)

          this_locate = locate(im,is)

          IF( .NOT. molecule_list(this_locate,is)%live) CYCLE
          
          IF( molecule_list(this_locate,is)%which_box /= this_box ) CYCLE          

          ! let us put the scaling for this molecule as 1.0. So no
          ! CFC particle. If its a fraction particle, then the scaling
          ! will be determined when the loop over k vectos is executed

          
          position = im_locate(im,is)
          this_scaling = 1.0_DP
          
          !$OMP PARALLEL DO DEFAULT(SHARED) &
          !$OMP PRIVATE(i,ia,hdotr,charge) &
          !$OMP SCHEDULE(STATIC)

          DO i = 1, nvecs(this_box)

             cos_mol(i,position) = 0.0_DP
             sin_mol(i,position) = 0.0_DP
             
             DO ia = 1, natoms(is)
                
                ! compute hdotr 
                
                hdotr = hx(i,this_box) * atom_list(ia,this_locate,is)%rxp + &
                        hy(i,this_box) * atom_list(ia,this_locate,is)%ryp + &
                        hz(i,this_box) * atom_list(ia,this_locate,is)%rzp
                
                charge = nonbond_list(ia,is)%charge
                
                cos_mol(i,position) = cos_mol(i,position) + charge * DCOS(hdotr)
                sin_mol(i,position) = sin_mol(i,position) + charge * DSIN(hdotr)
                
             END DO

             cos_sum(i,this_box) = cos_sum(i,this_box) + this_scaling * cos_mol(i,position)
             sin_sum(i,this_box) = sin_sum(i,this_box) + this_scaling * sin_mol(i,position)
             
          END DO
          
          !$OMP END PARALLEL DO
          
       END DO
       
    END DO
    
    ! At the end of all the loops we have computed cos_sum, sin_sum, cos_mol and sin_mol
    ! for each of the k-vectors. Now let us calcualte the reciprocal space energy 
    

    energy_temp = 0.0_DP

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(i,  un)  &
    !$OMP SCHEDULE(STATIC) &
    !$OMP REDUCTION(+:energy_temp)

    DO i = 1, nvecs(this_box)

       un =  cos_sum(i,this_box) * cos_sum(i,this_box) + & 
            sin_sum(i,this_box) * sin_sum(i,this_box)

       energy_temp = energy_temp + Cn(i,this_box) * un

    END DO

    !$OMP END PARALLEL DO
    
    energy(this_box)%ewald_reciprocal = energy_temp * charge_factor

  END SUBROUTINE Compute_System_Ewald_Reciprocal_Energy2

END MODULE Energy_Routines

