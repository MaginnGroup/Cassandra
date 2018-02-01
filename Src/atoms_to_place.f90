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


MODULE Atoms_To_Place

!********************************************************************
  ! This module contains a collection of subroutines that are used
  ! to determine number of the atoms and what atoms need to be regrown
  ! as a result of an internal coordinate move.
  !
  ! The moves are:
  !
  ! bond displacements
  ! angle fluctuations
  ! dihedral change
  !
  ! Used By
  !
  ! gcmc_control
  ! gemc_control
  ! nvtmc_control
  ! nptmc_control
  !
  ! Revision History
  !
  ! 12/10/13  : Beta Release 
  !********************************************************************

  USE Type_Definitions
  USE Global_Variables
  USE File_Names


  IMPLICIT NONE

  INTEGER, ALLOCATABLE, DIMENSION(:) :: alive, deadend, central
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atoms_to_place_list(:)
  
 
  
CONTAINS

!************************************************************************************************
  SUBROUTINE Get_Bonds_Atoms_To_Place
!************************************************************************************************
    ! This routine identifies the number of atoms and the indices of these atoms that will
    ! have to move when each atom on either side of a bond length is perturbed. It is assumed
    ! that only one of the bonded atoms will move and the displacement will be along the bond
    ! vector.
!************************************************************************************************

    !********************************************************************************************
    ! The routine loops over all the species and all the bonds of these species. First, atoms
    ! that form the bond are identified as atom1 and atom2. In the first pass, it is assumed
    ! that atom1 will be perturbed while atom2 is held fixed. A call to recursive_placement
    ! is then made to determine number of atoms that will move due to atom1 moving. The number
    ! of atoms and identity of these atoms are then stored in
    !
    ! bond_atoms_to_place_list(this_bond,this_species)%atom1_natoms
    ! bond_atoms_to_place_list(this_bond,this_species)%atom1(1:atom1_natoms)
    !
    ! The above procedure is next repeated for atom2 while atom1 assumed to be held fixed.
    ! The number of atoms that move and their identity are stroed in
    !
    ! bond_atoms_to_place_list(this_bond,this_species)%atom2_natoms
    ! bond_atoms_to_place_list(this_bond,this_species)%atom2(1:atom2_natoms)
    !********************************************************************************************


    IMPLICIT NONE
    
    INTEGER :: ibonds, atom1, atom2, ispecies
    INTEGER :: alive_atoms, atomplaced, iatoms,k, i, j
   
    ALLOCATE(atoms_to_place_list(MAXVAL(natoms)))
    ALLOCATE(alive(MAXVAL(natoms)))
    ALLOCATE(deadend(MAXVAL(natoms)))
    ALLOCATE(central(MAXVAL(natoms)))
    ALLOCATE(bond_atoms_to_place_list(MAXVAL(nbonds),nspecies),Stat=AllocateStatus)
    IF (AllocateStatus /=0 ) STOP
    
    DO ispecies = 1, nspecies
       
       DO ibonds = 1, nbonds(ispecies)
          
          ALLOCATE(bond_atoms_to_place_list(ibonds,ispecies)%atom1(MAXVAL(natoms)),Stat=AllocateStatus)
          IF (AllocateStatus /= 0 ) STOP
          ALLOCATE(bond_atoms_to_place_list(ibonds,ispecies)%atom2(MAXVAL(natoms)),Stat=AllocateStatus)
          IF (AllocateStatus /= 0 ) STOP

          ! identify atom1 and atom2 associated with global bond number ibonds in species ispecies
          atom1 = bond_list(ibonds,ispecies)%atom1
          atom2 = bond_list(ibonds,ispecies)%atom2

          ! initially set the total number of atoms to be placed equal to
          ! zero for both of these atoms
          alive_atoms = 0
          
          ! Initialize the array that contains the indices of all the atoms that 
          ! will move if atom1 and atom2 are displaced
          bond_atoms_to_place_list(ibonds,ispecies)%atom1(:) = 0
          bond_atoms_to_place_list(ibonds,ispecies)%atom2(:) = 0
          atoms_to_place_list = 0
          
          ! atom2 will already be placed when atom1 side molecule is
          ! regrown
          alive_atoms = alive_atoms + 1
          atoms_to_place_list(alive_atoms) = atom1

          ! bond_atoms_to_place_list(ibonds,ispecies)%atom1(alive_atom1) = atom1
          
          ! initialize a few arrays
          alive = 0
          deadend = 0
          central = 0

          ! make the other atom dead indicating that this part of the 
          ! molecule is not to be grown
          deadend(atom2) = 1  
          alive(atom1) = 1

          CALL recursive_placement(atom1,ispecies,alive_atoms,alive,central,deadend,atoms_to_place_list)
          
          ! transfer the number of atoms to be placed and the list of
          ! atoms to a permanent array

          bond_atoms_to_place_list(ibonds,ispecies)%atom1_natoms = alive_atoms
          
          DO iatoms = 1, alive_atoms
             
             bond_atoms_to_place_list(ibonds,ispecies)%atom1(iatoms) = &
                  atoms_to_place_list(iatoms)
             
          END DO

          ! now repeat the above procedure for the other side of the bond
          
          ! Since the identity of atom1 might have changed by recursive_placement routine. We
          ! determine the identity of atom1 once again.

          atom1 = bond_list(ibonds,ispecies)%atom1

          alive_atoms = 0
          alive = 0
          deadend = 0
          central = 0
          atoms_to_place_list  = 0
          
          alive_atoms = alive_atoms + 1
          atoms_to_place_list(alive_atoms) = atom2
          
          deadend(atom1) = 1
          alive(atom2) = 1

          CALL recursive_placement(atom2,ispecies,alive_atoms,alive,central,deadend,atoms_to_place_list)
          
          bond_atoms_to_place_list(ibonds,ispecies)%atom2_natoms = alive_atoms
          
          DO iatoms = 1, alive_atoms
             
             bond_atoms_to_place_list(ibonds,ispecies)%atom2(iatoms) = &
                  atoms_to_place_list(iatoms)
          END DO

       END DO

    END DO

    ! Free up the memory associated with alive, central and deadend

    IF(ALLOCATED(alive)) DEALLOCATE(alive)
    IF(ALLOCATED(deadend)) DEALLOCATE(deadend)
    IF(ALLOCATED(central)) DEALLOCATE(central)
    IF(ALLOCATED(atoms_to_place_list)) DEALLOCATE(atoms_to_place_list)

  END SUBROUTINE Get_Bonds_Atoms_To_Place
  
    
!**********************************************************************************
  SUBROUTINE Get_Angles_Atoms_To_Place
!**********************************************************************************
    ! This routine returns the number of atoms and identity of the atoms whose
    ! positions will have to be computed when an angle is distorted. It is assumed
    ! that angle distoration will be carried out randomly on either side of an 
    ! angle and only one of the end atoms will move.
    !******************************************************************************

    !******************************************************************************
    ! The program will loop over all the species and the angles contained in these
    ! species. For each angle, three atoms atom1,atom2,atom3 that form the angle 
    ! are identified. atom1 and atom3 are end atoms.
    !
    ! First it is assumed that the atom1 will move. A call to recursive_placement
    ! is made that returns the total number of atoms that will move when atom1 is
    ! perturbed and their identity. This information is stored in
    !
    ! angle_atoms_to_place_list(this_angle,this_species)%atom1_natoms
    ! angle_atoms_to_place_list(this_angle,this_species)%atom1(1:atom1_natoms)
    ! respectively.
    !
    ! Next atom3 is chosen as the atom that will move due to bond angle distortion. 
    ! The above procedure is repeated to obtain number of atoms that change position
    ! and their identity. The information is stored in
    !
    ! angle_atoms_to_place_list(this_angle,this_species)%atom3_natoms
    ! angle_atoms_to_place_list(this_angle,this_species)%atom3_natoms(1:atom3_natoms)
    !
    !*******************************************************************************
    
    IMPLICIT NONE
    
    INTEGER :: ispecies, iangles, alive_atoms, atom1, atom2, atom3, ia, k, i, j
    
    ALLOCATE(atoms_to_place_list(MAXVAL(natoms)))
    ALLOCATE(alive(MAXVAL(natoms)))
    ALLOCATE(central(MAXVAL(natoms)))
    ALLOCATE(deadend(MAXVAL(natoms)))
    
    ALLOCATE(angle_atoms_to_place_list(MAXVAL(nangles),nspecies),Stat=AllocateStatus)
    IF (AllocateStatus /= 0 ) STOP

    IF (verbose_log .AND. prob_angle > 0) THEN
       WRITE(logunit,*)
       write(logunit,'(A)') 'Atoms displaced by angle moves'
       WRITE(logunit,'(A)') '********************************************************************************'
    END IF
    
    species_loop:DO ispecies = 1, nspecies
       
       angles_loop: DO iangles = 1, nangles(ispecies)
          
          ALLOCATE(angle_atoms_to_place_list(iangles,ispecies)%atom1(MAXVAL(natoms)),Stat=AllocateStatus)
          IF (AllocateStatus /= 0 ) STOP
          ALLOCATE(angle_atoms_to_place_list(iangles,ispecies)%atom3(MAXVAL(natoms)),Stat=AllocateStatus)
          IF (AllocateStatus /= 0 ) STOP
          
          ! set the counters to zero indicating that nothing is to be regrown
          
          angle_atoms_to_place_list(iangles,ispecies)%atom1(:) = 0
          angle_atoms_to_place_list(iangles,ispecies)%atom3(:) = 0
          alive_atoms = 0
          atoms_to_place_list(:) = 0
          
          alive(:) = 0
          deadend(:) = 0
          central(:) = 0

          ! determine the atoms that form iangles

          atom1 = angle_list(iangles,ispecies)%atom1
          atom2 = angle_list(iangles,ispecies)%atom2
          atom3 = angle_list(iangles,ispecies)%atom3

          ! when the atoms, atom1 and all the atoms connected to it are regrown
          
          alive_atoms = alive_atoms + 1
          atoms_to_place_list(alive_atoms) = atom1

          deadend(atom2) = 1
          deadend(atom3) = 1
          alive(atom1) = 1
          
          CALL Recursive_Placement(atom1,ispecies,alive_atoms,alive,central,deadend,atoms_to_place_list)

          ! since the identity of atom1 might have changed due to the recursive routine, we reset
          ! this identity

          atom1 = angle_list(iangles,ispecies)%atom1
          
          ! assign all the atoms 
          
          angle_atoms_to_place_list(iangles,ispecies)%atom1_natoms = alive_atoms
          
          DO ia = 1, alive_atoms
             
             angle_atoms_to_place_list(iangles,ispecies)%atom1(ia) = &
                  atoms_to_place_list(ia)
             
          END DO

          ! Reassign the identity of atom1
          atom1 = angle_list(iangles,ispecies)%atom1
          
          alive_atoms = 0
          atoms_to_place_list(:) = 0
          
          alive(:) = 0
          deadend(:) = 0
          central(:) = 0
          
          ! increase the counter for alive_atoms by 1 for atom3
          
          alive_atoms = alive_atoms + 1
          atoms_to_place_list(alive_atoms) = atom3
          alive(atom3) = 1
          
          deadend(atom2) = 1
          deadend(atom1) = 1
          
          CALL Recursive_Placement(atom3,ispecies,alive_atoms,alive,central,deadend,atoms_to_place_list)
          
          angle_atoms_to_place_list(iangles,ispecies)%atom3_natoms = alive_atoms
          
          DO ia = 1, alive_atoms
             
             angle_atoms_to_place_list(iangles,ispecies)%atom3(ia) = &
                  atoms_to_place_list(ia)
             
          END DO

          ! echo the information in the log file

          IF (verbose_log .AND. prob_angle > 0) THEN
             write(logunit,*)'information about angle', iangles

             write(logunit,*) 'the atom that is moved is', atom1

             write(logunit,*)'atoms that change their position due to this move'
             write(logunit,*) (angle_atoms_to_place_list(iangles,ispecies)%atom1(i),i=1, &
               angle_atoms_to_place_list(iangles,ispecies)%atom1_natoms)

             write(logunit,*) 'the atoms is moved is', atom3
             write(logunit,*) 'atoms that change their position due to this move'
             write(logunit,*) (angle_atoms_to_place_list(iangles,ispecies)%atom3(i),i=1, &
               angle_atoms_to_place_list(iangles,ispecies)%atom3_natoms)
          END IF     
             
          
       END DO angles_loop

       
    END DO species_loop
    
    ! Free up the memory associated with alive, deadend and central

    IF (ALLOCATED(alive)) DEALLOCATE(alive)
    IF (ALLOCATED(central)) DEALLOCATE(central)
    IF (ALLOCATED(deadend)) DEALLOCATE(deadend)
    IF (ALLOCATED(atoms_to_place_list)) DEALLOCATE(atoms_to_place_list)

    IF (verbose_log .AND. prob_angle > 0) THEN
       WRITE(logunit,'(A)') '********************************************************************************'
    END IF
    
  END SUBROUTINE Get_Angles_Atoms_To_Place


!******************************************************************************************************
  SUBROUTINE Get_Dihedral_Atoms_To_Place
!******************************************************************************************************
    ! This subroutine detemines the number of atoms and their identity when a dihedral change is
    ! performed on either side of a dihedral angle. It is assumed that the dihedral change involves
    ! a rigid body rotation of these atoms due to dihedral angle change. 
    !**************************************************************************************************

    !**************************************************************************************************
    ! The routine loops over all the species and dihedral angles in these species. The atoms that
    ! form the dihedral angle are designated as atom1,atom2,atom3 and atom4. In the first pass, it
    ! is assumed that atom1 moves while the positions of the other atoms are held fixed. A call
    ! to recursive_placement is then made to determine the number of atoms moving and their identity.
    ! Note that the routine is written in such a way that all the atoms connected to atom2 (except
    ! atom3) will move so the original angles with vertex at atom2 are preserved. The number of
    ! atoms and their identity are then stored in
    !
    ! dihedral_atoms_to_place_list(this_dihedral,this_species)%atom1_natoms
    ! dihedral_atoms_to_place_list(this_dihedral,this_species)%atom1(1:atom1_natoms)
    !
    ! Next, atom4 is assumed to move while the other three atoms are held in their positions. The
    ! above procedure is repeated. Care is taken that all the atoms connected to atom3 (except atom2)
    ! move to preserve the original angles with vertex at atom3. The number of atoms and their 
    ! identity are then stored in
    !
    ! dihedral_atoms_to_place_list(this_dihedral,this_species)%atom4_natoms
    ! dihedral_atoms_to_place_list(this_dihedral,this_species)%atom4(1:atom4_natoms)
    !*************************************************************************************************
    
    IMPLICIT NONE

    INTEGER :: ispecies, idihedrals, atom1, atom2, atom3, atom4, iatoms
    INTEGER :: alive_atoms, k, i, j
    

    ALLOCATE(alive(MAXVAL(natoms)))
    ALLOCATE(deadend(MAXVAL(natoms)))
    ALLOCATE(central(MAXVAL(natoms)))
    ALLOCATE(atoms_to_place_list(MAXVAL(natoms)))

    ALLOCATE(dihedral_atoms_to_place_list(MAXVAL(ndihedrals),nspecies),Stat=AllocateStatus)
          IF (AllocateStatus /= 0 ) STOP

    species_loop: DO ispecies = 1, nspecies

       dihedral_loop: DO idihedrals = 1, ndihedrals(ispecies)

          ALLOCATE(dihedral_atoms_to_place_list(idihedrals,ispecies)%atom1(MAXVAL(natoms)),Stat=AllocateStatus)
          IF (AllocateStatus /= 0 ) STOP
          ALLOCATE(dihedral_atoms_to_place_list(idihedrals,ispecies)%atom4(MAXVAL(natoms)),Stat=AllocateStatus)
          IF (AllocateStatus /= 0 ) STOP

          ! assign 0 to all the atoms indicating that nothing is to be grown

          dihedral_atoms_to_place_list(idihedrals,ispecies)%atom1(:) = 0
          dihedral_atoms_to_place_list(idihedrals,ispecies)%atom4(:) = 0

          alive_atoms = 0
          atoms_to_place_list(:) = 0

          alive(:) = 0
          deadend(:) = 0
          central(:) = 0

          ! obtain atoms in this dihedral

          atom1 = dihedral_list(idihedrals,ispecies)%atom1
          atom2 = dihedral_list(idihedrals,ispecies)%atom2
          atom3 = dihedral_list(idihedrals,ispecies)%atom3
          atom4 = dihedral_list(idihedrals,ispecies)%atom4

          ! now when the dihedral move affects atom1 and other atoms connected to atom1

          deadend(atom3) = 1
          deadend(atom4) = 1

          ! atom1 will be the atom to be grown

          alive(atom1) = 1
         
          alive_atoms = alive_atoms + 1
          atoms_to_place_list(alive_atoms) = atom1

          ! Note we will make atom2 alive so that in the recursive routine it is not counted towards
          ! the atoms that is to be regrown for a dihedral change around bond 2-3. What it does though
          ! that it allows the atoms connected to atom2 to move with atom1 so that angle involving atom1
          ! are not distorted.

          alive(atom2) = 1

          CALL Recursive_Placement(atom1,ispecies,alive_atoms,alive,central,deadend,atoms_to_place_list)

          ! assign the atoms to the permanent list

          dihedral_atoms_to_place_list(idihedrals,ispecies)%atom1_natoms = alive_atoms

          DO iatoms = 1, alive_atoms

             dihedral_atoms_to_place_list(idihedrals,ispecies)%atom1(iatoms) = atoms_to_place_list(iatoms)

          END DO

          ! now hold atom1, atom2 and atom3 fixed so that dihedral move affects atom4
          ! reassign the indentity of atom1

          atom1 = dihedral_list(idihedrals,ispecies)%atom1

          alive_atoms = 0
          atoms_to_place_list(:) = 0

          alive(:) = 0
          deadend(:) = 0
          central(:) = 0

          deadend(atom1) = 1
          deadend(atom2) = 1

          alive(atom4) = 1

          alive_atoms = alive_atoms + 1
          atoms_to_place_list(alive_atoms) = atom4

          alive(atom3) = 1 ! see the note above

          CALL Recursive_Placement(atom4,ispecies,alive_atoms,alive,central,deadend,atoms_to_place_list)

          ! assign the atoms to be regrown to the list

          dihedral_atoms_to_place_list(idihedrals,ispecies)%atom4_natoms = alive_atoms

          DO iatoms = 1, alive_atoms

             dihedral_atoms_to_place_list(idihedrals,ispecies)%atom4(iatoms) = atoms_to_place_list(iatoms)

          END DO

       END DO dihedral_loop

    END DO species_loop
    
    IF (ALLOCATED(alive)) DEALLOCATE(alive)
    IF (ALLOCATED(central)) DEALLOCATE(central)
    IF (ALLOCATED(deadend)) DEALLOCATE(deadend)
    IF (ALLOCATED(atoms_to_place_list)) DEALLOCATE(atoms_to_place_list)    

    IF (verbose_log .AND. prob_torsion > 0) THEN
       WRITE(logunit,*)
       write(logunit,'(A)') 'Atoms displaced by dihedral moves'
       WRITE(logunit,'(A)') '********************************************************************************'
       DO i = 1, nspecies
          DO j = 1,ndihedrals(i)
             write(logunit,*) 'dihedral, atom1, atom2, atom3, atom4'
             write(logunit,*) j,dihedral_list(j,i)%atom1,dihedral_list(j,i)%atom2, dihedral_list(j,i)%atom3, &
                  dihedral_list(j,i)%atom4
             write(logunit,*)
             write(logunit,*) 'natoms to place if atom 1 moves and the moving atom numbers'
             write(logunit,*) dihedral_atoms_to_place_list(j,i)%atom1_natoms, dihedral_atoms_to_place_list(j,i)%atom1
             write(logunit,*) 'natoms to place if atom 4 moves and the moving atom numbers'
             write(logunit,*) dihedral_atoms_to_place_list(j,i)%atom4_natoms, dihedral_atoms_to_place_list(j,i)%atom4
             write(logunit,*)
          ENDDO
       ENDDO       
       WRITE(logunit,'(A)') '********************************************************************************'
    END IF

  END SUBROUTINE Get_Dihedral_Atoms_To_Place

!****************************************************************************************************
  RECURSIVE SUBROUTINE Recursive_Placement(atomplaced,ispecies,alive_atoms,alive,central,deadend, &
       atoms_to_place_list)
!****************************************************************************************************
    ! This is a recursive routine that walks from a given atom to determine if there is any path
    ! that has not been explored yet and counts the number of atoms along this path.
    ! A repeated call to this routine is made until all paths emanating from the first input atom 
    ! are fully explored. 
    !
    ! The routine is used to determine number of atoms that will move and their identity when one
    ! of the following moves is attempted:
    !
    ! bond length perturbation
    ! angle distortion
    ! dihedral angle change
    ! 
    ! Correspondingly the routine is called by
    !
    !********** Used by *******************
    ! Get_Bond_Atoms_To_Place
    ! Get_Angle_Atoms_To_Place
    ! Get_Dihedral_Atoms_To_Place
    !***********************************************************************************************

    !***********************************************************************************************
    ! The first call to this routine will have an input atom identity from which a molecule
    ! is to be grown. The basic algorithm involves placing all the atoms connected to the input
    ! atom. Then we identify additional atom connected to the atoms bonded to  the original input
    ! atom and so on. The routine will return to the calling routine once all the atoms have been
    ! accounted for.
    !
    ! To accomplish this task, three arrays have been defined alive,central and deadend. 
    ! alive(atom_id) = 1 means that atom_id needs to be counted towards atoms that need to be placed.
    ! central(atom_id) = 1 implies that all the atoms bonded to atom_id have been accounted for.
    ! deadend(atom_id) = 1 implies several things:
    !   - all the atoms connected to atom_id have been placed
    !   - one or more atoms connected to atom_id are terminal atoms ( i.e. are bonded to only one atom)
    !   - atom_id itself is a terminal atom.
    !
    ! The routine will first determine if the input atom (atomplaced) is a deadend. 
    ! If not, then it will determine if the atoms connected to atomplaced are already placed, if not
    ! then we will make the alive array for these atoms to be equal to 1 indicating that these atoms
    ! need to be placed and also increment alive_atoms counter to keep track of number of atoms to be
    ! placed. Since we placed all the atoms connected to atomplaced, central(atomplaced) = 1. 
    !
    ! Now, we will choose one of the atoms bonded to atomplaced to see if we can grow molecule in the
    ! direction of this atom. This will occur only if the next atom is not a terminal atom, all the
    ! atoms connected to the next atoms are not already accounted for and it is not a deadend. 
    !
    ! If we cannot grow from any of the atoms connected to atomplaced then we encountered a deadend
    ! at atomplaced and start retracing the path until we find a new path to grow molecule from.
    !
    ! The process is repeated (recursive call to the routine) until all the paths have been exhausted
    !
    ! In the end, the routine returns number of atoms whose new positions need to be computed
    ! and their indentity due to an internal coordinate change.
    !**************************************************************************************************

    IMPLICIT NONE
    
    INTEGER :: atomplaced, i, bond_i, atom1i, atom2i, nextatom
    INTEGER :: n_connect, j, bond_j, atom1j, atom2j, ispecies, alive_atoms

    INTEGER, DIMENSION(:) :: alive, deadend, central, atoms_to_place_list

 

    IF (deadend(atomplaced) == 1) THEN
       
       ! We encountered a terminal atom or all the atoms connected to atomplaced have
       ! already been accounted for.
       ! start walking back from atomplaced. 
       
       DO i = 1, bondpart_list(atomplaced,ispecies)%nbonds
          
          bond_i = bondpart_list(atomplaced,ispecies)%bond_num(i)
          
          atom1i = bond_list(bond_i,ispecies)%atom1
          atom2i = bond_list(bond_i,ispecies)%atom2
          
          IF ( atom1i == atomplaced ) THEN
             
             nextatom = atom2i
             
          ELSE
             
             nextatom = atom1i
             
          ENDIF
          
          IF (bondpart_list(nextatom,ispecies)%nbonds /= 1) THEN
             IF (deadend(nextatom) == 0) THEN
                ! there are atoms to along a path originating at nextatom
                CALL Recursive_Placement(nextatom,ispecies,alive_atoms,alive,central,deadend,atoms_to_place_list)
             END IF
          END IF
          
       END DO
       
    ELSE
      
       IF (central(atomplaced) /= 1) THEN

          ! this is the first time we are encountering this atom. So we will place
          ! all the nearest neighbors of this atom

          central(atomplaced) = 1
          
          ! n_connect is the numbers of bonds attached to atom atomplaced (which is determined already 
          ! from the participation routine.
          n_connect = bondpart_list(atomplaced,ispecies)%nbonds

          DO j = 1, n_connect
             
             ! figure out index of the atoms connected to 'atomplaced'
             
             ! bond_j is the global bond number associated with local bond number j
             bond_j = bondpart_list(atomplaced,ispecies)%bond_num(j)
             
             ! Now find the global atom numbers on either side of global bond number bond_j
             atom1j = bond_list(bond_j,ispecies)%atom1
             atom2j = bond_list(bond_j,ispecies)%atom2
             
             IF ( atom1j == atomplaced) THEN
                
                nextatom = atom2j
                
             ELSE
                
                nextatom = atom1j
                
             END IF
             
             ! There is a possibility that this is the other side of the
             ! atom of a bond or angle or dihedral that is not to be 
             ! counted towards number of atoms to be placed due to an
             ! internal coordinate move. We have already set these
             ! atoms to be deadend in the calling program. So exclude
             ! this possibility to begin with.

             IF (deadend(nextatom) == 0) THEN

                ! if nextatom has not already been accounted for
                
                IF (alive(nextatom) == 0) THEN
                   alive(nextatom) = 1
                   alive_atoms = alive_atoms + 1
                   atoms_to_place_list(alive_atoms) = nextatom
                END IF
                
             END IF
             
          END DO
          
       END IF

       
       ! change the identity of atomplaced
       
       DO j = 1, bondpart_list(atomplaced,ispecies)%nbonds
          
          bond_j = bondpart_list(atomplaced,ispecies)%bond_num(j)
          
          atom1j = bond_list(bond_j,ispecies)%atom1
          atom2j = bond_list(bond_j,ispecies)%atom2
          
          IF (atom1j == atomplaced) THEN
             
             nextatom = atom2j
             
          ELSE
             
             nextatom = atom1j
             
          END IF
          
          
          ! -- Check to see that it is not a terminal bond
          
          IF (bondpart_list(nextatom,ispecies)%nbonds &
               /= 1) THEN
             ! This is not a terminal bond. 
             ! Make sure we have not placed all the atoms connected to it
             
             IF ( deadend(nextatom) /= 1 ) THEN
                
                ! if atoms connected to it are not already placed then
                
                IF ( central(nextatom) /= 1 ) THEN
                   
                   ! change the identity
                   
                   CALL Recursive_Placement(nextatom,ispecies,alive_atoms,alive,central,deadend,atoms_to_place_list)
                   
                END IF
                
             END IF
             
          END IF
          
       END DO
       
       ! if I am here it means that I encountered a deadend and
       ! must retrace to make sure that all the atoms are placed
       
       ! since it is a recursive call if may have already been
       ! assigned a deadend
       
       IF ( deadend(atomplaced) /= 1 ) THEN
          
          deadend(atomplaced) = 1
          
          CALL Recursive_Placement(atomplaced,ispecies,alive_atoms,alive,central,deadend,atoms_to_place_list) 
          
       END IF
       
    END IF
    
  END SUBROUTINE Recursive_Placement
  
END MODULE Atoms_To_Place

