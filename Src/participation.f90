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

!*****************************************************************************
SUBROUTINE Participation
!*****************************************************************************

!**********************************************************************
! This subroutine returns connectivity information of atoms in all
! the species. That is, given an atom, it calculates how many bonds 
! it participates in and what other atoms are connected to the atom.
! It assigns bond number indices local to an atom (and not the same
! as the bond number index for the species that is read in from the 
! connectivity file).


  ! Note, that the following assignment below:
  ! bondpart_list(iatom,ispecies)%bond_num(i) = ibonds

  ! sets a local bond number (for a given atom) i equal to the overall bond number of the 
  ! particular species (ibonds). That is, bondpart_list(iatom,ispecies)%bond_num(1), 
  ! bondpart_list(iatom,ispecies)%bond_num(2), etc. are the 1st, 2nd, etc. bonds attached
  ! to atom iatom in species ispecies. The values of these array elements are the overall
  ! bond numbers of the species, set in the bond_list array.

  ! bondpart_list(iatom,ispecies)%atom(1) = an atom number
  ! gives the glboal atom number in the species (an atom number) connected to atom via its 1st bond while
  ! bondpart_list(iatom,ispecies)%atom(2) = an atom number
  ! gives the glboal atom number in the species connected to atom via its 2nd bond, etc.
  ! bondpart_list(iatom,ispecies)%nbonds is the total number of bonds attached to atom iatom
  ! in species ispecies.

  ! Written by: Jindal Shah
  ! Date: October, 2007

  ! *** Called by *** 
  ! NVTMC_Control

  ! *** Calls ***
  ! None

  ! Revision history:
  
  ! 12/10/13  : Beta version 
!**********************************************************************

  USE Type_Definitions
  USE Run_Variables
  USE IO_Utilities
  USE File_Names
  USE Random_Generators, ONLY : Generate_Random_Sphere

  IMPLICIT NONE

  INTEGER :: ispecies, iatom, jatom, tot_bonds, ibonds, atom1, atom2, i ,j, k
  INTEGER :: atom3, atom4, idihedral, tot_dihedral, tot_angles
  INTEGER :: iangles, ndisp_atoms, is, ia, this_atom, this_bond, frag_no, ia_nangles
  INTEGER :: this_angle, first_atom, third_atom, anchor_atom, ifrag, atom_j
  INTEGER, ALLOCATABLE :: temp_atom_id(:)

  INTEGER :: line_nbr, pdbunit, ierr

  REAL(DP) :: x_this, y_this, z_this, this_l

  CHARACTER(120) :: file_name, car_file, xyz_file
  CHARACTER(120) :: line_string

  

  ! allocate memory for bondpart_list. 
  ALLOCATE(bondpart_list(MAXVAL(natoms),nspecies))
  ALLOCATE(l_bonded(MAXVAL(natoms),MAXVAL(natoms),nspecies))

  ! set l_bonded to false. set them to true based on the
  ! molecule connectivity data
  
  DO ispecies = 1, nspecies
     DO iatom = 1,natoms(ispecies)
        DO jatom = 1,natoms(ispecies)
           l_bonded(iatom,jatom,ispecies) = .false.
        ENDDO
     ENDDO 
  ENDDO

  DO ispecies = 1,nspecies
     DO ibonds = 1,nbonds(ispecies)
        iatom = bond_list(ibonds,ispecies)%atom1
        jatom = bond_list(ibonds,ispecies)%atom2
        l_bonded(iatom,jatom,ispecies)=.true.
        l_bonded(jatom,iatom,ispecies)=.true.
     ENDDO       
  ENDDO   

!  DO ispecies = 1, nspecies
!     DO iatom = 1,natoms(ispecies)
!        DO jatom = 1,natoms(ispecies)
!           WRITE(logunit,*) ispecies,iatom,jatom,l_bonded(iatom,jatom,ispecies) 
!        ENDDO
!     ENDDO
!  ENDDO


  species_list(:)%ndisp_atoms = 0
  species_list(:)%f_atom_disp = .FALSE.


  DO ispecies = 1, nspecies

     ndisp_atoms = 0
     
     ALLOCATE(species_list(ispecies)%disp_atom_id(MAXVAL(natoms)))
     ALLOCATE(species_list(ispecies)%disp_atom_ref(MAXVAL(natoms)))

     species_list(ispecies)%disp_atom_id(:) = 0
     species_list(ispecies)%disp_atom_ref(:) = 0

     DO iatom = 1, natoms(ispecies)

        ALLOCATE(bondpart_list(iatom,ispecies)%bond_num(MAXVAL(nbonds)))
        ALLOCATE(bondpart_list(iatom,ispecies)%atom(MAXVAL(nbonds)))

        tot_bonds = 0

        DO ibonds = 1, nbonds(ispecies)
           
           atom1 = bond_list(ibonds,ispecies)%atom1
           atom2 = bond_list(ibonds,ispecies)%atom2

           IF (atom1 == iatom .OR. atom2 == iatom) THEN

              ! iatom is present in the bond

              tot_bonds = tot_bonds + 1

              IF (tot_bonds > max_connectivity) THEN
                 
                 err_msg =''
                 err_msg(1) = 'Connectivity of atom'
                 err_msg(2) = Int_To_String(iatom)
                 err_msg(3) = 'exceeds the maximum allowed'

                 CALL Clean_Abort(err_msg,'Participation')
                 
              END IF
              ! assign the bond number as well
              bondpart_list(iatom,ispecies)%bond_num(tot_bonds) = ibonds
              
              IF (atom1 == iatom) THEN

                 bondpart_list(iatom,ispecies)%atom(tot_bonds) = atom2

              ELSE

                 bondpart_list(iatom,ispecies)%atom(tot_bonds) = atom1

              END IF
              
           END IF
           
        END DO

        ! Store the total number of bonds for iatoms

        bondpart_list(iatom,ispecies)%nbonds = tot_bonds
     
        ! Note that if the atom is connected to more than two atoms
        ! it is a branch point. Also if atom is connected to only one
        ! atom, it is a terminal point
        
        ! Make use of this fact to identify atoms that may be moved using
        ! spherical coordinate changes

        IF (tot_bonds == 1) THEN
           ! iatom is a terminal point. Increment number of atoms that can
           ! be moved

           ndisp_atoms = ndisp_atoms + 1
           species_list(ispecies)%disp_atom_id(ndisp_atoms) = iatom
           species_list(ispecies)%disp_atom_ref(ndisp_atoms) = &
                bondpart_list(iatom,ispecies)%atom(1)
           
        END IF
           
     END DO

     ! assign total number of atoms to be displaced

     species_list(ispecies)%ndisp_atoms = ndisp_atoms
     IF (ndisp_atoms > 0 ) THEN
        species_list(ispecies)%f_atom_disp = .TRUE.
     END IF

  END DO

  ! Echo to logfile for checking
  WRITE(logunit,*) '*** Computed bond participation information ***'
  DO ispecies = 1, nspecies

     DO iatom = 1, natoms(ispecies)

        DO ibonds = 1, bondpart_list(iatom,ispecies)%nbonds
           WRITE(logunit,*) 'species, atom, local bond nbr, global bond nbr, atom bonded '
           WRITE(logunit,'(3X,I2,7X,I3,5X,I3,10X,I3,15X,I3)') ispecies,iatom,ibonds,bondpart_list(iatom,ispecies)%bond_num(ibonds), bondpart_list(iatom,ispecies)%atom(ibonds)
        ENDDO
        WRITE(logunit,*)
     ENDDO
     WRITE(logunit,*)
  ENDDO

  DO ispecies = 1, nspecies
     IF(species_list(ispecies)%f_atom_disp) THEN
        
        WRITE(logunit,*)
        WRITE(logunit,*) 'For species', ispecies
        WRITE(logunit,*) ' Number of atoms that can be displaced', species_list(ispecies)%ndisp_atoms
        WRITE(logunit,'(2(A20,2x))') 'Atom ID', 'Reference Atom ID'
        DO iatom = 1, species_list(ispecies)%ndisp_atoms
           WRITE(logunit,'(2(I20,2x))') species_list(ispecies)%disp_atom_id(iatom), &
                species_list(ispecies)%disp_atom_ref(iatom)
        END DO
     END IF
  END DO
        
  
  WRITE(logunit,*)

  !***************************************************************************
  ! This part of the program computes the participation matrix for angles.
  ! i.e., given an atom, we will figure out how many angles the atoms 
  ! participates in and the indices of these angles.
  !
  ! The information will be stored in
  !
  ! angle_part_list(iatom,ispecies)%nangles == total number of angles iatom
  !       of ispecies participates in
  ! angle_part_list(iatom,ispecies)%which_angle(1) gives the global angle
  !       number for iatom
  ! angle_part_list(iatom,ispecies)%which_angle(2) gives the second global
  !       angle number and so on.
  ! angle_part_list(iatom,ispecies)%position(1) will give the position of the
  !       iatom in the first angle it is encountered in. The position can be 
  !       1,2 or 3.
  !****************************************************************************

  ! allocate memory for the angle participation list

  ALLOCATE(angle_part_list(MAXVAL(natoms),nspecies))

  DO ispecies = 1, nspecies

     DO iatom = 1, natoms(ispecies)

        ! loop over all the angles

        tot_angles = 0

        ALLOCATE(angle_part_list(iatom,ispecies)%which_angle(MAXVAL(nangles)))
        ALLOCATE(angle_part_list(iatom,ispecies)%position(MAXVAL(nangles)))

        DO iangles = 1, nangles(ispecies)

           ! obtain three atoms that define the angle

           atom1 = angle_list(iangles,ispecies)%atom1
           atom2 = angle_list(iangles,ispecies)%atom2
           atom3 = angle_list(iangles,ispecies)%atom3

           IF ( iatom == atom1 ) THEN

              tot_angles = tot_angles + 1
              angle_part_list(iatom,ispecies)%which_angle(tot_angles) = iangles
              angle_part_list(iatom,ispecies)%position(tot_angles) = 1
              
              CYCLE

           END IF
              
           IF ( iatom == atom2 ) THEN
    
              tot_angles = tot_angles + 1
              angle_part_list(iatom,ispecies)%which_angle(tot_angles) = iangles
              angle_part_list(iatom,ispecies)%position(tot_angles) = 2
              
              CYCLE

           END IF

           IF ( iatom == atom3 ) THEN

              tot_angles = tot_angles + 1
              angle_part_list(iatom,ispecies)%which_angle(tot_angles) = iangles
              angle_part_list(iatom,ispecies)%position(tot_angles) = 3
              
              CYCLE

           END IF

        END DO

        ! assign the total number of angles 

        angle_part_list(iatom,ispecies)%nangles = tot_angles

        ! echo the information to the log file for checking
        
        write(logunit,*)'Number of atoms in question', iatom
        write(logunit,*)'Total number of angles', angle_part_list(iatom,ispecies)%nangles
        write(logunit,*)'Angles for this atom'
        write(logunit,*)(angle_part_list(iatom,ispecies)%which_angle(i),i=1, &
             angle_part_list(iatom,ispecies)%nangles)
        write(logunit,*)'Position of the atom in these dihedrals'
        write(logunit,*)(angle_part_list(iatom,ispecies)%position(i), i=1, &
             angle_part_list(iatom,ispecies)%nangles)

     END DO

  END DO
        


  ! Now we will compute the number of dihedrals a given atoms participates in.
  ! What these dihedral angles are and position of the atom in a given dihedral.
  ! i.e. 1,2,3 and 4 based on the dihedral angle definition supplied in the input
  ! file. Note the positions 1 and 4 are equivalent so are 2 and 3. The variable
  ! dihedral_part computed in this routine will be useful in angle distortion move
  ! where we will pick a dihedral angle at random for the atom that is being moved
  ! to determine its new position. The reader is referred to angle_distortion.f90
  ! for further details.
  !
  ! dihedral_part_list(iatom,ispecies)%ndihedrals === total number of dihedrals
  !       iatom of ispecies participates in
  ! dihedral_part_list(iatom,ispecies)%which_dihedral(1) gives the global dihedral
  !       angle number for iatom of ispecies in which it is present
  ! dihedral_part_list(iatom,ispecies)%which_dihedral(2) gives the second global
  !       dihedral angle for iatom and so on.
  !
  ! dihedral_part_list(iatom,ispecies)%position(1) gives the position of the iatom
  !       in dihedral_part_list(iatom,ispecies)%which_dihedral(1) 
  !


  ! Allocate memeory for the participation list

  ALLOCATE(dihedral_part_list(MAXVAL(natoms),nspecies))

  DO ispecies = 1, nspecies

     DO iatom = 1, natoms(ispecies)

        ! loop over all the dihedrals to determine how many dihedrals contain
        ! iatom

        tot_dihedral = 0

        ALLOCATE(dihedral_part_list(iatom,ispecies)%which_dihedral(MAXVAL(ndihedrals)))
        ALLOCATE(dihedral_part_list(iatom,ispecies)%position(MAXVAL(ndihedrals)))

        DO idihedral = 1, ndihedrals(ispecies)

           ! get the four atoms in the dihedral

           atom1 = dihedral_list(idihedral,ispecies)%atom1
           atom2 = dihedral_list(idihedral,ispecies)%atom2
           atom3 = dihedral_list(idihedral,ispecies)%atom3
           atom4 = dihedral_list(idihedral,ispecies)%atom4

           IF ( iatom == atom1 ) THEN

              tot_dihedral = tot_dihedral + 1
              dihedral_part_list(iatom,ispecies)%which_dihedral(tot_dihedral) = idihedral
              dihedral_part_list(iatom,ispecies)%position(tot_dihedral) = 1
              
              CYCLE
              
           END IF

           IF ( iatom == atom2 ) THEN

              tot_dihedral = tot_dihedral + 1
              dihedral_part_list(iatom,ispecies)%which_dihedral(tot_dihedral) = idihedral
              dihedral_part_list(iatom,ispecies)%position(tot_dihedral) = 2

              CYCLE

           END IF

           IF ( iatom == atom3 ) THEN

              tot_dihedral = tot_dihedral + 1
              dihedral_part_list(iatom,ispecies)%which_dihedral(tot_dihedral) = idihedral
              dihedral_part_list(iatom,ispecies)%position(tot_dihedral) = 3

              CYCLE

           END IF

           IF ( iatom == atom4 ) THEN

              tot_dihedral = tot_dihedral + 1
              dihedral_part_list(iatom,ispecies)%which_dihedral(tot_dihedral) = idihedral
              dihedral_part_list(iatom,ispecies)%position(tot_dihedral) = 4

              CYCLE

           END IF

        END DO

        ! assign the total number of dihedrals found for iatom here

        dihedral_part_list(iatom,ispecies)%ndihedrals = tot_dihedral

        ! echo the information to the log file for checking

        write(logunit,*)'Number of atoms in question', iatom
        write(logunit,*)'Total number of dihedrals', dihedral_part_list(iatom,ispecies)%ndihedrals
        write(logunit,*)'Dihedral angles for this atom'
        write(logunit,*)(dihedral_part_list(iatom,ispecies)%which_dihedral(i),i=1, &
             dihedral_part_list(iatom,ispecies)%ndihedrals)
        write(logunit,*)'Position of the atom in angles'
        write(logunit,*)(dihedral_part_list(iatom,ispecies)%position(i), i=1, &
             dihedral_part_list(iatom,ispecies)%ndihedrals)
        
     END DO

  END DO

  IF (int_sim_type == sim_mcf) THEN
     
     ! Write the fragments mcf file if any of the atoms is a fragment
     
     
     DO is = 1, nspecies
        
        DO ifrag = 1, nfragments(is)

!           IF ( frag_list(ifrag,is)%ring) CYCLE

           ! locate the anchor of this fragment. Since non-ring
           ! fragments have only 1 anchor we will use this to
           ! identify the id of the anchor.
         
           ia = frag_list(ifrag,is)%anchor(1)
           
           anchor_atom = 1
              
           frag_no = ifrag
           file_name = 'frag_'//Int_To_String(frag_no)
           file_name = file_name(1:LEN_TRIM(file_name))//"_"//Int_To_String(is)
           file_name = file_name(1:LEN_TRIM(file_name))//".mcf"
           
           car_file = 'frag_'//Int_To_String(frag_no)
           car_file = car_file(1:LEN_TRIM(car_file))//"_"//Int_To_String(is)
           car_file = car_file(1:LEN_TRIM(car_file))//".car"
           
           xyz_file = 'frag_'//Int_To_String(frag_no)
           xyz_file = xyz_file(1:LEN_TRIM(xyz_file))//"_"//Int_To_String(is)
           xyz_file = xyz_file(1:LEN_TRIM(xyz_file))//".xyz"           
           
           
           OPEN(UNIT=201,file=file_name)
           OPEN(UNIT=202,file=car_file)
           OPEN(UNIT=203,file=xyz_file)
           
           ! write fragment file
           
           WRITE(201,'(A)') '#Species_Type'
           WRITE(201,'(A)') 'NORMAL'
           WRITE(201,*) 
           
           WRITE(201,'(A)') '# Atom_Info'

           IF (frag_list(ifrag,is)%ring) THEN

              ! The temp_atom_id list will contain the original atom ids of the fragment. For example
              ! if there are 5 atoms in the ring fragment with original atom ids in the molecule as
              ! 3 5, 8, 9 and 10...and assuming that this also the order in which they are listed in
              ! the Fragment_Info section, we will have
              ! tempp_atom_id(1) = 3, temp_atom_id(2) = 5 and so on...

              WRITE(201,*) frag_list(ifrag,is)%natoms
              ALLOCATE(temp_atom_id(frag_list(ifrag,is)%natoms))
              

              DO i = 1, frag_list(ifrag,is)%natoms

                 this_atom = frag_list(ifrag,is)%atoms(i)
                 
                 ! assign this atom to temp_atom_id array

                 temp_atom_id(i) = this_atom

                 ! Check to see if this atom is a ring fragment, if so append, 'ring' at the end

                 IF (nonbond_list(this_atom,is)%ring_atom) THEN
                    
                    WRITE(201,'(I5,2X,2(A4,2X),2(F11.7,2X),A6,2X,2(F11.7, 2X),A4)') i, &
                         nonbond_list(this_atom,is)%atom_name, nonbond_list(this_atom,is)%element, &
                         nonbond_list(this_atom,is)%mass, nonbond_list(this_atom,is)%charge, &
                         nonbond_list(this_atom,is)%vdw_potential_type, &
                         nonbond_list(this_atom,is)%vdw_param(1)/kboltz, &
                         nonbond_list(this_atom,is)%vdw_param(2), 'ring'
                 ELSE

                    WRITE(201,'(I5,2X,2(A4,2X),2(F11.7,2X),A6,2X,2(F11.7, 2X))')i, &
                         nonbond_list(this_atom,is)%atom_name, nonbond_list(this_atom,is)%element, &
                         nonbond_list(this_atom,is)%mass, nonbond_list(this_atom,is)%charge, &
                         nonbond_list(this_atom,is)%vdw_potential_type, &
                         nonbond_list(this_atom,is)%vdw_param(1)/kboltz, &
                         nonbond_list(this_atom,is)%vdw_param(2)
                    
                 END IF
                 
              END DO

           ELSE
              
              WRITE(201,*) bondpart_list(ia,is)%nbonds + 1
              WRITE(201,100) anchor_atom, nonbond_list(ia,is)%atom_name, nonbond_list(ia,is)%element, &
                   nonbond_list(ia,is)%mass, nonbond_list(ia,is)%charge, 'NONE', &
                   (nonbond_list(ia,is)%vdw_param(1))/kboltz, nonbond_list(ia,is)%vdw_param(2)
              
              ! write the 'read_old' file for the fragment
              
              WRITE(203,*) '1'  ! indicating one molecule
              
              x_this = 0.0_DP
              y_this = 0.0_DP
              z_this = 0.0_DP
              
              WRITE(202,'(A,3F11.7)') nonbond_list(ia,is)%atom_name, x_this, y_this, z_this
              WRITE(203,'(A,3F11.7)') nonbond_list(ia,is)%atom_name, x_this, y_this, z_this
              
              ! The order of the atom listing will be the same as that given in the master mcf file
              
              DO i = 2, frag_list(ifrag,is)%natoms
                 
                 this_atom = frag_list(ifrag,is)%atoms(i)
                 
                 WRITE(201,100) i, nonbond_list(this_atom,is)%atom_name, nonbond_list(this_atom,is)%element, &
                      nonbond_list(this_atom,is)%mass, nonbond_list(this_atom,is)%charge, 'NONE', &
                      (nonbond_list(this_atom,is)%vdw_param(1))/kboltz, nonbond_list(this_atom,is)%vdw_param(2)
                 
              END DO
              
           END IF  ! IF (frag_list(ifrag,is)%ring)
           
           !---------------------------------------------------------------------------------
           
           ! now let us write the Bond_Info section
           
           IF (frag_list(ifrag,is)%ring) THEN

              CALL Write_Ring_Fragment_Car_File(ifrag,is)
              CALL Write_Ring_Fragment_MCF_Bond_Info(ifrag,is)

           ELSE
           
              WRITE(201,*)
              WRITE(201,'(A)') '# Bond_Info'
              WRITE(201,*) bondpart_list(ia,is)%nbonds 
              
              ALLOCATE(temp_atom_id(bondpart_list(ia,is)%nbonds)) ! to hold atom ids of connections
              
              ! We will place the atoms such that the first atom is at the origin, the second along
              ! the x axis and rest placed randomly on a sphere. The algorithm works as follows
              
              ! for every non-anchor atom, determine the bond number from the bond_part(ia,i)
              ! use this information to find out the bond length and corresponding harmonic
              ! force constant if bond vibrations are included
              
              DO i = 2, frag_list(ifrag,is)%natoms
                 
                 this_atom = frag_list(ifrag,is)%atoms(i)
                 
                 this_bond = 0
                 
                 DO j = 1, bondpart_list(ia,is)%nbonds
                    
                    atom_j = bondpart_list(ia,is)%atom(j)
                    
                    IF (atom_j == this_atom) THEN
                       
                       ! We found the bond bond between this_atom and anchor
                       
                       this_bond = bondpart_list(ia,is)%bond_num(j)
                       temp_atom_id(i-1) = this_atom
                       EXIT
                       
                    END IF
                    
                 END DO
                 
                 IF (this_bond == 0 ) THEN
                    ! We could not locate the bond connecting this_atom and anchor, abort with an error
                    
                    err_msg = ''
                    err_msg(1) = 'Error while generating a fragment mcf file for the fragment '//Int_To_String(ifrag)
                    err_msg(2) = 'Species '//Int_To_String(is)
                    err_msg(3) = 'Fragment ' //Int_To_String(ifrag)
                    err_msg(4) = ' No bond found between the atoms'//Int_To_String(ia)
                    err_msg(5) = ' and '//Int_To_String(this_atom)
                    err_msg(6) = 'Check the fragment info section in the master mcf file for this fragment'
                    CALL Clean_Abort(err_msg,'participation.f90')
                    
                 END IF
                 
                 IF (bond_list(this_bond,is)%int_bond_type == int_harmonic) THEN
                    
                    WRITE(201,101) i-1, anchor_atom, i, "harmonic", bond_list(this_bond,is)%bond_param(1) /kboltz, &
                         bond_list(this_bond,is)%bond_param(2)
                    
                 ELSE IF (bond_list(this_bond,is)%int_bond_type == int_none) THEN
                    ! it is a fixed bond
                    
                    WRITE(201,102) i-1, anchor_atom, i, "fixed", bond_list(this_bond,is)%bond_param(1)
                    
                    ! for a fixed bond length system, generate points of this atom on a unit sphere
                    
                    this_l = bond_list(this_bond,is)%bond_param(1)
                 
                    IF ( i == 2) THEN
                       ! this is the first bond and hence the second atom, it will be placed along 
                       ! the x-axis
                       
                       x_this = this_l
                       y_this = 0.0_DP
                       z_this = 0.0_DP
                       
                    ELSE
                       
                       CALL Generate_Random_Sphere(x_this,y_this,z_this)
                       
                       ! coordinates of the 'this_atom'
                       
                       
                       x_this = this_l * x_this
                       y_this = this_l * y_this
                       z_this = this_l * z_this
                       
                    END IF
                    
                    ! Now write the coordinates to the car file/xyz file
                    
                    WRITE(202,'(A,3F11.7)') nonbond_list(this_atom,is)%atom_name, x_this, y_this, z_this
                    WRITE(203,'(A,3F11.7)') nonbond_list(this_atom,is)%atom_name, x_this, y_this, z_this
                    
                    
                 END IF
                 
              END DO
              
              CLOSE(UNIT=202)
              CLOSE(UNIT=203)
           
           END IF ! IF(frag_list(ifrag,is)%ring)
           
           ! now let us write the Angle_Info section
           
           ! first count the number of angles in which this atom is in the central
           ! position
           
           
           ! now write the section on angles

           IF (frag_list(ifrag,is)%ring) THEN
              CALL Write_Ring_Fragment_MCF_Angle_Info(ifrag,is)
           ELSE

              ia_nangles = 0
              
              DO i = 1, angle_part_list(ia,is)%nangles
                 
                 IF ( angle_part_list(ia,is)%position(i) == 2 ) THEN
                    
                    ia_nangles = ia_nangles + 1
                    
                 END IF
                 
              END DO
              
              WRITE(201,*)
              WRITE(201,'(A)') '# Angle_Info'
              WRITE(201,*) ia_nangles

              ia_nangles = 0
              
              DO i = 1, angle_part_list(ia,is)%nangles
                 
                 IF (angle_part_list(ia,is)%position(i) == 2) THEN
                    
                    ia_nangles = ia_nangles + 1
                    
                    this_angle = angle_part_list(ia,is)%which_angle(i)
                    
                    first_atom = angle_list(this_angle,is)%atom1
                    third_atom = angle_list(this_angle,is)%atom3
                    
                    DO j = 1, bondpart_list(ia,is)%nbonds
                       
                       IF (temp_atom_id(j) == first_atom) THEN
                          first_atom = j + 1
                          EXIT
                          
                       END IF
                       
                    END DO
                    
                    DO j = 1, bondpart_list(ia,is)%nbonds
                       
                       IF (temp_atom_id(j) == third_atom) THEN
                          third_atom = j + 1
                          EXIT
                       END IF
                    END DO
                    
                    IF (angle_list(this_angle,is)%int_angle_type == int_harmonic) THEN
                       
                       WRITE(201,103) ia_nangles, first_atom, anchor_atom, third_atom, "harmonic",  &
                            angle_list(this_angle,is)%angle_param(1) / kboltz, &
                            (180.0_DP/PI) * angle_list(this_angle,is)%angle_param(2)
                    
                    ELSE
                       
                       WRITE(201,104) ia_nangles, first_atom, anchor_atom, third_atom, "fixed",  &
                            angle_list(this_angle,is)%angle_param(1)
                       
                    END IF
                    
                 END IF
                 
              END DO

           END IF ! IF (frag_list(ifrag,is)%ring)

           ! --- Remaining sections

           IF (frag_list(ifrag,is)%ring) THEN

              CALL Write_Ring_Fragment_MCF_Dihedral_Info(ifrag,is)
              CALL Write_Ring_Fragment_MCF_Improper_Info(ifrag,is)

           ELSE
              
              WRITE(201,'(A)') '# Dihedral_Info'
              WRITE(201,'(A)') '0'
              
              WRITE(201,*) 
              WRITE(201,'(A)') '# Improper_Info'
              WRITE(201,'(A)') '0'
              WRITE(201,*)
              
              WRITE(201,'(A)') 'END'

           END IF
              
           DEALLOCATE(temp_atom_id)
           CLOSE(UNIT=201)
           
        END DO
     
     END DO
     !Amir To Jindal: It the 3rd # format should be 11.7 otherwise it would generate errors. 10/12/12
	 
100  FORMAT(I5,2X,A4,2X,A4,F11.7,2X,F11.7,2X,A5,2X,F11.7, 2X, F11.7)
101  FORMAT(I5,2X,I5,2X,I5,2X,A9,2X,F10.3,2X,F8.5)
102  FORMAT(I5,2X,I5,2X,I5,2X,A9,2X,F8.5)
103  FORMAT(I5,2X,I5,2X,I5,2X,I5,2X,A9,2X,F10.3,2X,F10.5)
104  FORMAT(I5,2X,I5,2X,I5,2X,I5,2X,A9,2X,F10.5)
105  FORMAT(A,2X,I5,2X,A2,2X,6(I3,2X))
     
  END IF

CONTAINS

  SUBROUTINE Write_Ring_Fragment_MCF_Bond_Info(ifrag,is)

    ! The subroutine writes the '# Bond_Info' section in the ring fragment mcf file

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ifrag, is

    ! Local variables

    INTEGER :: i, this_anchor, j, atom_j, k, this_atom, this_bond, bond_counter
    INTEGER, ALLOCATABLE, DIMENSION(:) :: atom1, atom2, bond_id
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: bond_placed

    ALLOCATE(bond_placed(nbonds(is)))
    ALLOCATE(atom1(nbonds(is)), atom2(nbonds(is)), bond_id(nbonds(is)))
    bond_placed(:) = .FALSE.
    atom1 = 0 ; atom2 =0 ; bond_id = 0
    bond_counter = 0

    ! We will loop over all the anchor atoms of the ring fragment and place all
    ! the connections of each of the anchors ensuring that there are no repetition.
    !
    ! This will be achieved by identifying number of connections of each anchor 
    ! by bondpart_list()%nbonds. We will have to determine the total number of
    ! bonds, so a counter is required. Also, for each of the bonds found in
    ! the fragment, we will keep track of the two atoms in that bond and its
    ! global bond id. 

    DO i = 1, frag_list(ifrag,is)%nanchors

       this_anchor = frag_list(ifrag,is)%anchor(i)

       DO j = 1, bondpart_list(this_anchor,is)%nbonds

          atom_j = bondpart_list(this_anchor,is)%atom(j)
          
          ! find the local id of this atom in the fragment

          DO k = 1, frag_list(ifrag,is)%natoms

             this_atom = frag_list(ifrag,is)%atoms(k)

             IF (this_atom == atom_j) THEN
                this_bond = bondpart_list(this_anchor,is)%bond_num(j)
                IF (.NOT. bond_placed(this_bond)) THEN
                   
                   ! write this bond number to the Bond Info section and quit
                   ! we will also mark this bond number as placed so that it can be skipped
                   bond_placed(this_bond) = .TRUE.
                   
                   ! increment the bond counter 
                   bond_counter = bond_counter + 1
                   atom1(bond_counter) = i
                   atom2(bond_counter) = k
                   bond_id(bond_counter) = this_bond

                   ! EXIT

                END IF
                
             END IF

          END DO ! loop over all the fragment atoms

       END DO ! loop over all the bonds for this_anchor

    END DO ! loop over all the anchors, i

    ! Now write the bond info section to the file

    WRITE(201,*)
    WRITE(201,'(A)') '# Bond_Info'
    WRITE(201,'(I3)') bond_counter

    ! Write the bond information for each of the bonds

    DO i = 1, bond_counter
       this_bond = bond_id(i)

       IF (bond_list(this_bond,is)%int_bond_type == int_harmonic) THEN
          WRITE(201,'(I5,2X,I5,2X,I5,2X,A9,2X,F10.3,2X,F8.5)') i, &
               atom1(i), atom2(i), "harmonic", bond_list(this_bond,is)%bond_param(1)/kboltz, &
               bond_list(this_bond,is)%bond_param(2)

       ELSE IF (bond_list(this_bond,is)%int_bond_type == int_none) THEN
          WRITE(201,'(I5,2X,I5,2X,I5,2X,A9,2X,F8.5)') i, &
               atom1(i), atom2(i), "fixed", bond_list(this_bond,is)%bond_param(1)

       END IF

    END DO ! loop over the bond_counter
             
  END SUBROUTINE Write_Ring_Fragment_MCF_Bond_Info
!-----------------------------------------------------------------------------------------------

  SUBROUTINE Write_Ring_Fragment_MCF_Angle_Info(ifrag,is)
    !------------------------------------------------------------------------------------------
    ! Write the '# Angle_Info' section for the ring fragment
    !
    !------------------------------------------------------------------------------------------
    
    IMPLICIT NONE 

    INTEGER, INTENT(IN) :: ifrag, is

    ! Local variables

    INTEGER :: i, this_atom, j, frag_1_nangles, frag_2_nangles
    INTEGER :: this_angle, first_atom, second_atom, third_atom
    INTEGER :: local_id_1, local_id_2, local_id_3, k, atom_k
  

    ! First find out the total number of angles for the ring fragment
    frag_1_nangles = 0

    DO i = 1, frag_list(ifrag,is)%nanchors

       this_atom = frag_list(ifrag,is)%anchor(i)

       DO j = 1, angle_part_list(this_atom,is)%nangles
          
          IF (angle_part_list(this_atom,is)%position(j) == 2) THEN
             frag_1_nangles = frag_1_nangles + 1
          END IF

       END DO ! loop over all the angles, the anchor 'i' participates in

    END DO ! loop over all the anchors

    !----------------------------------------------------------------------
    WRITE(201,*)
    WRITE(201,'(A)') '# Angle_Info'
    WRITE(201,'(I3)') frag_1_nangles
    
    ! Now identify the atoms of each of the angles, find local ids and print the
    ! information

    frag_2_nangles = 0

    DO i = 1, frag_list(ifrag,is)%nanchors

       this_atom = frag_list(ifrag,is)%anchor(i)

       DO j = 1, angle_part_list(this_atom,is)%nangles

          IF (angle_part_list(this_atom,is)%position(j) == 2) THEN

             frag_2_nangles = frag_2_nangles + 1

             this_angle = angle_part_list(this_atom,is)%which_angle(j)

             first_atom = angle_list(this_angle,is)%atom1
             second_atom = angle_list(this_angle,is)%atom2
             third_atom = angle_list(this_angle,is)%atom3

             ! Now look for these two atoms in the fragment list to find their
             ! local id

             local_id_1 = 0; local_id_2 = 0; local_id_3 = 0

             DO k = 1, frag_list(ifrag,is)%natoms

                atom_k = frag_list(ifrag,is)%atoms(k)
    
                IF (atom_k == first_atom) local_id_1 = k
                IF (atom_k == second_atom) local_id_2 = k
                IF (atom_k == third_atom) local_id_3 = k

             END DO

             IF (local_id_1 == 0 .OR. local_id_2 == 0 .OR. local_id_3 == 0) THEN
                err_msg = ''
                err_msg(1) = 'One of the local ids is zero'
                err_msg(2) = 'This occurred while writing the angle info section'
                err_msg(3) = 'for the ring fragment mcf file'
                err_msg(4) = 'Fragment '//TRIM(Int_To_String(ifrag))
                err_msg(5) = 'Species '//TRIM(Int_To_String(is))
                CALL Clean_Abort(err_msg, 'Write_Ring_Fragment_MCF_Angle_Info')

             END IF

             ! This angle will go in the mcf file

             IF (angle_list(this_angle,is)%int_angle_type == int_harmonic) THEN
                
                WRITE(201,'(I5,2X,I5,2X,I5,2X,I5,2X,A9,2X,F10.3,2X,F10.5)') frag_2_nangles, &
                     local_id_1, local_id_2, local_id_3, "harmonic", &
                     angle_list(this_angle,is)%angle_param(1)/kboltz, &
                     angle_list(this_angle,is)%angle_param(2) * (180.0_DP/PI)

             ELSE
                
                WRITE(201,'(I5,2X,I5,2X,I5,2X,I5,2X,A9,2X,F10.3)') frag_2_nangles, &
                     local_id_1, local_id_2, local_id_3, "fixed", &
                     angle_list(this_angle,is)%angle_param(1)

             END IF


          END IF ! if the anchor is at the vertex

       END DO ! Do loop over the total number of angles, the anchor 'i' participates in

    END DO ! Do loop over all the anchors

    ! Ensure that all the angles were printed

    IF (frag_1_nangles /= frag_2_nangles) THEN
       err_msg = ''
       err_msg(1) ='An error occurred while writing the angle info section'
       err_msg(2) ='for the ring fragment '//TRIM(Int_To_String(ifrag))
       err_msg(3) ='in species '//TRIM(Int_To_String(is))
       err_msg(4) = 'Number of angles in the fragment is '//TRIM(Int_To_String(frag_1_nangles))
       err_msg(5) = 'However only '//TRIM(Int_To_String(frag_2_nangles))
       err_msg(6) = 'written.'
       CALL Clean_Abort(err_msg,'Write_Ring_Fragment_MCF_Angle_Info')
    END IF

  END SUBROUTINE Write_Ring_Fragment_MCF_Angle_Info
  !--------------------------------------------------------------------------------------------
  SUBROUTINE Write_Ring_Fragment_MCF_Dihedral_Info(ifrag,is)
    !------------------------------------------------------------------------------------------
    ! Write the '# Dihderal Info' section to the ring fragment MCF file
    !------------------------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ifrag, is

    ! Local variables

    INTEGER :: idihed, atom_1, atom_2, atom_3, atom_4
    INTEGER :: ncounted, dihed_counter, i, this_atom, this_dihedral
    INTEGER, ALLOCATABLE, DIMENSION(:) :: first_atom, second_atom, third_atom, fourth_atom
    INTEGER, ALLOCATABLE, DIMENSION(:) :: dihed_id
    

    ALLOCATE(first_atom(ndihedrals(is)), second_atom(ndihedrals(is)))
    ALLOCATE(third_atom(ndihedrals(is)), fourth_atom(ndihedrals(is)))
    ALLOCATE(dihed_id(ndihedrals(is)))

    first_atom = 0 ; second_atom = 0 ; third_atom = 0 ; fourth_atom = 0 ; dihed_id = 0
    

    ! loop over all the dihedrals. If all the four atoms of the dihedral are part of the
    ! the dihedral write it out to the file
    !
    ! first determine total number of dihedrals that will be written to the file

    dihed_counter = 1

    DO idihed = 1, ndihedrals(is)

       atom_1 = dihedral_list(idihed,is)%atom1
       atom_2 = dihedral_list(idihed,is)%atom2
       atom_3 = dihedral_list(idihed,is)%atom3
       atom_4 = dihedral_list(idihed,is)%atom4

       ! loop over all the fragment atoms and see if all the four atoms are part of
       ! the fragement

       ncounted = 0

       DO i = 1, frag_list(ifrag,is)%natoms

          this_atom = frag_list(ifrag,is)%atoms(i)

          IF (this_atom == atom_1) THEN

             first_atom(dihed_counter) = i
             ncounted = ncounted + 1

          ELSE IF (this_atom == atom_2) THEN
             
             second_atom(dihed_counter) = i
             ncounted = ncounted + 1

          ELSE IF (this_atom == atom_3) THEN

             third_atom(dihed_counter) = i
             ncounted = ncounted + 1

          ELSE IF(this_atom == atom_4) THEN

             fourth_atom(dihed_counter) = i
             ncounted = ncounted + 1

          END IF

          IF (ncounted == 4) THEN
             ! All the atoms in the 'idihed' are in the fragment and hence
             ! they are to be written to the file
             dihed_id(dihed_counter) = idihed
             dihed_counter = dihed_counter + 1
             
             EXIT

          END IF
                      
       END DO

    END DO

    ! actual number of dihedrals in the fragment is 1 less
    
    dihed_counter = dihed_counter - 1

    WRITE(201,*) 
    WRITE(201,'(A)') '# Dihedral_Info'
    WRITE(201,'(I4)') dihed_counter

    ! Write all the information to the files

    DO i = 1, dihed_counter

       this_dihedral = dihed_id(i)

       IF (dihedral_list(this_dihedral,is)%int_dipot_type == int_opls) THEN

          WRITE(201,'(I5,2X,4(I4,2X),A4,2X,4(F8.3,2X))') i, &
               first_atom(i), second_atom(i), third_atom(i), fourth_atom(i), &
               'OPLS', &
               dihedral_list(this_dihedral,is)%dihedral_param(1)/kjmol_to_atomic, &
               dihedral_list(this_dihedral,is)%dihedral_param(2)/kjmol_to_atomic, &
               dihedral_list(this_dihedral,is)%dihedral_param(3)/kjmol_to_atomic, &
               dihedral_list(this_dihedral,is)%dihedral_param(4)/kjmol_to_atomic

       ELSE IF (dihedral_list(this_dihedral,is)%int_dipot_type == int_charmm) THEN
          
          WRITE(201,'(I5,2X,4(I4,2X), A6, 2X, F8.3, 2X, F8.0, 2X, F8.3)') i, &
               first_atom(i), second_atom(i), third_atom(i), fourth_atom(i), &
               'CHARMM', &
               dihedral_list(this_dihedral,is)%dihedral_param(1)/kjmol_to_atomic, &
               dihedral_list(this_dihedral,is)%dihedral_param(2), &
               dihedral_list(this_dihedral,is)%dihedral_param(3) * ( 180.0_DP/PI)

       ELSE IF (dihedral_list(this_dihedral,is)%int_dipot_type == int_amber) THEN

          WRITE(201,'(I5,2X,4(I4,2X),A5,2X,9(F8.3,2X))') i, &
               first_atom(i), second_atom(i), third_atom(i), fourth_atom(i), &
               'AMBER', &
               dihedral_list(this_dihedral,is)%dihedral_param(1)/kjmol_to_atomic, &
               dihedral_list(this_dihedral,is)%dihedral_param(2), &
               dihedral_list(this_dihedral,is)%dihedral_param(3) * (180_DP/PI), &
               dihedral_list(this_dihedral,is)%dihedral_param(4)/kjmol_to_atomic, &
               dihedral_list(this_dihedral,is)%dihedral_param(5), &
               dihedral_list(this_dihedral,is)%dihedral_param(6) * (180_DP/PI), &
               dihedral_list(this_dihedral,is)%dihedral_param(7)/kjmol_to_atomic, &
               dihedral_list(this_dihedral,is)%dihedral_param(8), &
               dihedral_list(this_dihedral,is)%dihedral_param(9) * (180_DP/PI)
               
       ELSE IF (dihedral_list(this_dihedral,is)%int_dipot_type == int_harmonic) THEN
          
          WRITE(201,'(I5,2X, 4(I4,2X), A8,2X, 2(F8.3,2X))') i, &
               first_atom(i), second_atom(i), third_atom(i), fourth_atom(i), &
               'harmonic', &
               dihedral_list(this_dihedral,is)%dihedral_param(1)/ kboltz, &
               dihedral_list(this_dihedral,is)%dihedral_param(2) * (180.0_DP/PI) 

       ELSE IF (dihedral_list(this_dihedral,is)%int_dipot_type == int_none) THEN

          WRITE(201,'(I5,2X,4(I4,2X),A4)') i, &
               first_atom(i), second_atom(i), third_atom(i), fourth_atom(i), 'none'

          
       END IF

    END DO

  END SUBROUTINE Write_Ring_Fragment_MCF_Dihedral_Info
  !------------------------------------------------------------------------------------------
  SUBROUTINE Write_Ring_Fragment_MCF_Improper_Info(ifrag,is)
    !----------------------------------------------------------------------------------------
    !
    ! The subroutine outputs improper angle info to the ring fragment mcf file
    !
    !----------------------------------------------------------------------------------------


    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ifrag, is

    INTEGER :: improper_counter, imp_ang, atom_1, atom_2, atom_3, atom_4, ncounted, this_atom
    INTEGER :: this_improper, i
    
    INTEGER, DIMENSION(:), ALLOCATABLE :: first_atom, second_atom, third_atom, fourth_atom, improper_id

    ALLOCATE(first_atom(nimpropers(is)), second_atom(nimpropers(is)))
    ALLOCATE(third_atom(nimpropers(is)), fourth_atom(nimpropers(is)))
    ALLOCATE(improper_id(nimpropers(is)))

    ! loop over all the improper angles in the molecules. If all the four atoms of a given
    ! improper angle are part of the ring fragment then the improper angle gets written.

    ! Determine the total number of impropers and all the atoms that for these impropers

    improper_counter = 1

    DO imp_ang = 1, nimpropers(is)

       atom_1 = improper_list(imp_ang,is)%atom1
       atom_2 = improper_list(imp_ang,is)%atom2
       atom_3 = improper_list(imp_ang,is)%atom3
       atom_4 = improper_list(imp_ang,is)%atom4

       ! loop over the fragment atoms and determine if all the atoms belong to the fragment

       ncounted = 0

       DO i = 1, frag_list(ifrag,is)%natoms
          
          this_atom = frag_list(ifrag,is)%atoms(i)
          
          IF (this_atom == atom_1) THEN
             
             first_atom(improper_counter) = i
             ncounted = ncounted + 1
             
          ELSE IF (this_atom == atom_2) THEN
             
             second_atom(improper_counter) = i
             ncounted = ncounted + 1

          ELSE IF (this_atom == atom_3) THEN

             third_atom(improper_counter) = i
             ncounted = ncounted + 1

          ELSE IF (this_atom == atom_4) THEN
             
             fourth_atom(improper_counter) = i
             ncounted = ncounted + 1
             
          END IF ! IF (this_atom == atom_1)

          IF (ncounted == 4 ) THEN
             ! all the atoms in the 'imp_angle' are in the fragment and
             ! hence they are to be written to the file
             improper_id(improper_counter) = imp_ang
             improper_counter = improper_counter + 1

             EXIT
             
          END IF
    
       END DO ! DO i = 1, frag_list(ifrag,is)%natoms

    END  DO ! DO imp_ang = 1, nimpropers(is)

    ! actual number of impropers in the fragment is 1 less

    improper_counter = improper_counter - 1
    
    WRITE(201,*) 
    WRITE(201,'(A)') '# Improper_Info'
    WRITE(201,*) improper_counter

    DO i = 1, improper_counter
       
       this_improper = improper_id(i)

       IF (improper_list(this_improper,is)%improper_potential_type == 'harmonic') THEN

          WRITE(201,'(I5,2X,4(I4,2X),A8,2X,2(F8.3,2X))') i, &
               first_atom(i), second_atom(i), third_atom(i), fourth_atom(i), &
               'harmonic', &
               improper_list(this_improper,is)%improper_param(1)/kboltz, &
               improper_list(this_improper,is)%improper_param(2) * (180_DP/PI)
               
       ELSE IF (improper_list(this_improper,is)%improper_potential_type == 'none') THEN
          
          WRITE(201,'(I5,2X,4(I4,2X), A4)') i, &
               first_atom(i), second_atom(i), third_atom(i), fourth_atom(i), &
               'none'

       END IF

    END DO

    WRITE(201,'(A)') 'END'

    CLOSE(UNIT=201)
    

  END SUBROUTINE Write_Ring_Fragment_MCF_Improper_Info
  !----------------------------------------------------------------------------------------
  SUBROUTINE Write_Ring_Fragment_Car_File(ifrag,is)
    !--------------------------------------------------------------------------------------
    !
    ! The subroutine writes the car file for a ring fragment
    !
    !--------------------------------------------------------------------------------------
   
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ifrag, is

    INTEGER :: i, this_atom, j
    REAL(DP) :: x, y, z
    
    CHARACTER(120) :: car_file, pdb_file, xyz_file
    CHARACTER(30) :: this_symbol

    
    pdb_file = 'molecule.pdb'

    car_file = 'frag_'//Int_To_String(ifrag)
    car_file = car_file(1:LEN_TRIM(car_file))//"_"//Int_To_String(is)
    car_file = car_file(1:LEN_TRIM(car_file))//".car"
    
    xyz_file = car_file(1:LEN_TRIM(car_file)-4)//'.xyz'

    ! loop over the fragment atoms
    OPEN(UNIT=202,file=car_file)
    OPEN(UNIT=204,file=xyz_file)
    WRITE(204,*)'1'

    DO i = 1, frag_list(ifrag,is)%natoms
       
       this_atom = frag_list(ifrag,is)%atoms(i)

       ! open the pdb file and go to the line for 'this_atom'
       ! read in the coordinates and output to the car file

       line_nbr = 0
       pdbunit = 203
       OPEN(UNIT=pdbunit,file=pdb_file)

       ! Go through the file and find out where 'this_atom' is
       
       DO

         line_nbr = line_nbr + 1

         CALL Read_String(pdbunit,line_string,ierr)
         IF (ierr .NE. 0) THEN
            err_msg = ""
            err_msg(1) = "Error reading the pdb file"
            err_msg(2) = pdb_file
            CALL Clean_Abort(err_msg,'Read_inputfile')
         END IF

         IF (line_string(1:6) == 'HETATM' .OR. line_string(1:4) == 'ATOM') THEN
            backspace(pdbunit)
            line_nbr = line_nbr + 1
            DO j = 1, this_atom-1
               READ(pdbunit,*)
            END DO
            ! Now read the coordinates for this atom
            READ(pdbunit,'(A30,3(F8.3))') this_symbol, x, y, z
            
            WRITE(202,'(A,2X,3(F11.7,2X))') nonbond_list(this_atom,is)%atom_name, x,y, z
            WRITE(204,'(A,2X,3(F11.7,2X))') nonbond_list(this_atom,is)%atom_name, x,y, z

            EXIT

         ELSEIF (line_string(1:3) == 'END') THEN

            EXIT 
            
         ELSEIF (line_nbr > 10000 ) THEN
            ! No name specified so abort
            err_msg = ""
            err_msg(1) = 'Exceeded file size in'
            err_msg(2) = pdb_file
            err_msg(3) = 'The atom'
            err_msg(4) = Int_To_String(this_atom)
            err_msg(5) = 'could not be found'
            CALL Clean_Abort(err_msg,'Get_Runname')

        ENDIF


       ENDDO 
       CLOSE(UNIT=pdbunit)
      

    END DO
    CLOSE(UNIT=202)
    CLOSE(UNIT=204)
 
 
    

  END SUBROUTINE Write_Ring_Fragment_Car_File

END SUBROUTINE Participation
