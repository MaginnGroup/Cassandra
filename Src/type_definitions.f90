!********************************************************************************
!   Cassandra open source atomistic Monte Carlo software package
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

!********************************************************************************
MODULE Type_Definitions
!********************************************************************************

  !****************************************************************************
  ! defines a few key parameters and then declares a series of variable types
  ! that are associated with particular objects.

  ! USED BY
  !
  !      atoms_to_place
  !      compute_cell_dimensions
  !      create_nonbond_table
  !      get_internal_coords
  !      input_routines
  !      io_utilities
  !      nvtmc_control
  !      participation
  !      random_generators ONLY : DP
  !      global_variables
  !      translate
  !
  ! Revision history:
  !
  !
  ! 08/02/13 : Beta release version
  !
  ! 03/17/19 (JS) : Activity field defined for species class for GCMC simulations.
  !****************************************************************************
  IMPLICIT NONE

  !Create a double precision definition, to make numbers machine independent
  INTEGER, PARAMETER :: SP = KIND(1.0e0)
  INTEGER, PARAMETER :: DP = KIND(1.0d0)

  ! Specify the limits on the parameters for various intramolecular function classes
  INTEGER, PARAMETER :: max_bond_params = 5
  INTEGER, PARAMETER :: max_angle_params = 5
  INTEGER, PARAMETER :: max_dihedral_params = 10
  INTEGER, PARAMETER :: max_improper_params = 5
  INTEGER, PARAMETER :: max_nonbond_params = 10

  ! Include the limit on connectivity
  INTEGER, PARAMETER :: max_connectivity = 6

  ! Place a limit on the number of angle evaluations for probability calculations
  INTEGER, PARAMETER :: nregions = 1000
  
  ! Define some classes to hold variables associated with different objects
  ! in the simulation. These will be converted to lists in global_variables for speed.

  !****************************************************************************
  TYPE Species_Class

     ! Basic information of this particular species. 
     ! Species list will have dimensions (nspecies)

     REAL(DP) :: molecular_weight, total_charge

     ! simulation specific information
     CHARACTER(20) :: species_type, insertion

     ! State point dependent information
     REAL(DP) :: chem_potential, activity
     REAL(DP) :: max_lambda, max_torsion
     REAL(DP), ALLOCATABLE :: de_broglie(:)

     ! COM of the specified initial configuration
     REAL(DP) :: xcom, ycom, zcom

     ! following variables are defined so that an atom in a molecule may be subjected
     ! to atom displacement via spherical coordinates.
     !
     ! ndisp_atoms : number of atoms that can be displaced by this routine
     ! disp_atom_id : list containing atom ids of ndisp_atoms
     ! disp_atom_ref : list containing reference atoms about which disp_atom_id is moved
     ! f_atom_disp : logical flag that is true if a species contains atoms to be displaced
     !              

     INTEGER :: ndisp_atoms
     INTEGER :: nmoltotal
     INTEGER :: int_species_type, int_insert
     INTEGER, DIMENSION(:), ALLOCATABLE :: disp_atom_id, disp_atom_ref
     LOGICAL :: f_atom_disp

     LOGICAL :: fragment
     LOGICAL :: linear
     ! NR: Adding to have an option not to include 
     ! Coul interaction during biased growth
     LOGICAL :: L_Coul_CBMC 
     ! N.B. natoms, max_molecules, etc. are in a separate arrays
     ! NR: for insertion style
     LOGICAL :: lcom 

  END TYPE Species_Class
  !****************************************************************************




  !****************************************************************************
  TYPE Molecule_Class

     ! The molecule list will have dimensions (max_molecules,nspecies)

     ! What kind of molecule is this? normal, fractional, fixed, etc. 
     ! Note that the following integers will be defined for the type
     ! int_normal = 0
     ! int_fractional = 1
     ! int_fixed = 2
     INTEGER :: molecule_type 
     INTEGER :: rx_num
     
     ! for open system and multi-box simulations (GEMC, parallel tempering)
     LOGICAL :: live
     LOGICAL :: inside
     INTEGER :: which_box

     ! also include com information for each of the molecules. 
     ! com and euler angles refer to the x,y and z com coordinates
     ! and 1, 2 and 3 euler angles of the molecule. The suffix
     ! old denotes old coordinates.
     ! frac is the fractional scaling parameter for the molecule
     REAL(DP)  :: xcom, ycom, zcom, euler1, euler2, euler3
     REAL(DP)  :: xcom_old, ycom_old, zcom_old, euler1_old,euler2_old,euler3_old
     REAL(DP)  :: frac
     ! This variable records the maximum distance of any psuedo atom from its
     ! COM. This is used to speed up energy calculations.

     REAL(DP) :: max_dcom, max_dcom_old
    


  END TYPE Molecule_Class
  !****************************************************************************



  !****************************************************************************
  TYPE Internal_Coord_Class

     ! The internal coordinate list will have dimensions 
     ! (MAXVAL(nbonds), MAXVAL(max_molecules), MAXVAL(nspecies))

     REAL(DP) :: bond_length_angstrom
     REAL(DP) :: bond_angle_degrees, bond_angle_radians
     REAL(DP) :: dihedral_angle_degrees, dihedral_angle_radians
     REAL(DP) :: improper_angle_degrees, improper_angle_radians
     
  END TYPE Internal_Coord_Class
  !****************************************************************************

  TYPE Internal_Coord_Class_Old
     
     ! The internal_coord_list_old will have dimension of (MAXVAL(nbonds))
     ! The type is identical to Internal_Coord_Class and is used to store
     ! the old internal coordinates of a molecule during the move.

     REAL(DP) :: bond_length_angstrom
     REAL(DP) :: bond_angle_degrees, bond_angle_radians
     REAL(DP) :: dihedral_angle_degrees, dihedral_angle_radians
     REAL(DP) :: improper_angle_degrees, improper_angle_radians
     
  END TYPE Internal_Coord_Class_Old

  !****************************************************************************
  TYPE Atom_Class    

     ! Information for each atom in system. Defines the location of the atom.  
     ! The p variables (rxp,rxyp,rzp) denote parent coordinates     
     ! exists here is used to signify if the particular atom exists yet or not
     ! in the system. For example, if we are growing a molecule, some atoms may
     ! not yet have been placed, and so exist = 'false'.

     ! atom_list has dimensions (natoms, max_molecules, nspecies)
     
     REAL(DP) :: rxp, ryp, rzp
     REAL(DP) :: rxp_nls, ryp_nls, rzp_nls  ! The starting positions for the neighbor list
     REAL(DP) :: rxp_old, ryp_old, rzp_old
     LOGICAL :: exist

  END TYPE Atom_Class
  !****************************************************************************



  !****************************************************************************
  TYPE Nonbond_Class

     ! Information on non-bonded interactions stored here. The type of vdw interaction 
     ! associated with a given atom is stored, along with its mass, charge 
     ! and vdw parameters.  Store an element and atom name. A global type number
     ! is assigned to each atom of the same name. 
     
     ! We also include here whether this atom is a ring atom 

     ! nonbond list has dimensions (MAXVAL(natoms), nspecies)

     CHARACTER(20) :: vdw_type
     REAL(DP), DIMENSION(max_nonbond_params) :: vdw_param

     CHARACTER(2) :: element
     CHARACTER(6) :: atom_name

     REAL(DP) :: mass, charge
     INTEGER :: atom_type_number

     LOGICAL :: ring_atom

  END TYPE Nonbond_Class
  !****************************************************************************



  !****************************************************************************
  TYPE Bond_Class

     ! bond list has dimensions (MAXVAL(nbonds), nspecies)
     INTEGER :: atom1, atom2, int_bond_type
     REAL(DP), DIMENSION(max_bond_params) :: bond_param
     CHARACTER(20) :: bond_potential_type

  END TYPE Bond_Class
  !****************************************************************************



  !****************************************************************************
  TYPE Angle_Class

     ! angle list has dimensions (MAXVAL(nangles,nspecies)

     ! 1 - 2 - 3 forms an angle with 2 at the apex.
     INTEGER :: atom1, atom2, atom3

     REAL(DP), DIMENSION(max_angle_params) :: angle_param
     CHARACTER(20) :: angle_potential_type
     INTEGER :: int_angle_type

  END TYPE Angle_Class
  !****************************************************************************




  !****************************************************************************
  TYPE Dihedral_Class

     ! dihedral list has dimensions (MAXVAL(ndihedrals, nspecies)

     ! Describes a standard dihedral: the four atoms are held with atom numbers
     ! defined sequentially (1-2-3-4) along dihedral angle, and parameters for the
     ! potential are held in an array. 

     INTEGER :: atom1, atom2, atom3, atom4
     REAL(DP), DIMENSION(max_dihedral_params) :: dihedral_param
     CHARACTER(20) :: dihedral_potential_type
     INTEGER :: int_dipot_type

  END TYPE Dihedral_Class
  !****************************************************************************



  !****************************************************************************
  TYPE Improper_Class

     ! improper list has dimensions (MAXVAL(nimpropers), nspecies)

     ! Describes an improper dihedral: atom 1 is the central atom, and atoms 2 
     ! through 4 are attached to atom 1, parameters are stored in an array
     
     INTEGER :: atom1, atom2, atom3, atom4
     REAL(DP), DIMENSION(max_improper_params) :: improper_param
     CHARACTER(20) :: improper_potential_type
     INTEGER :: int_improp_type
     
  END TYPE Improper_Class
  !****************************************************************************


  !****************************************************************************
  TYPE Bond_Atoms_To_Place_Class
     
     ! bond_atoms_to_place_list has dimensions of (MAXVAL(nbonds),nspecies)
     ! it is designed to hold the number of atoms and the identity of these
     ! atoms to be regrown due to a bond length move.

     INTEGER :: atom1_natoms
     INTEGER :: atom2_natoms
     INTEGER, DIMENSION(:), ALLOCATABLE :: atom1
     INTEGER, DIMENSION(:), ALLOCATABLE :: atom2
     
  END TYPE Bond_Atoms_To_Place_Class
  
  

  !****************************************************************************
  TYPE Angle_Atoms_To_Place_Class

     ! angles atoms to place list has dimensions (MAXVAL(nangles),nspecies)
     ! atom1 and atom3 will have dimensions of MAXVAL(natoms)

     ! Describes the number of atoms that need to be regrown due an internal
     ! coordinate move involving angles. atom1 and atom3 are the end atoms of
     ! the given angle.

     ! atom1 and atom3 hold the indices of the atoms that need to be regrown.
     ! atom1_num and atom3_num contain total number of atoms are regrown.

     INTEGER, DIMENSION(:), ALLOCATABLE :: atom1(:)
     INTEGER, DIMENSION(:), ALLOCATABLE :: atom3(:)
     INTEGER :: atom1_natoms
     INTEGER :: atom3_natoms

  END TYPE Angle_Atoms_To_Place_Class
  !****************************************************************************

  !****************************************************************************
  TYPE Dihedral_Atoms_to_Place_Class

     ! dihedral_atoms_to_place_list has dimensions of (MAXVAL(dihedrals), nspecies)
     ! The elements of the list, atom1 and atom4, will have diemnsions of
     ! (MAXVAL(natoms)). Atom1 holds the indices of all the atoms for which
     ! positions need to computed when dihedral is rotated on atom1 side
     ! atom2 holds the corresponding indices for the dihedral move on the other ised.

     INTEGER, DIMENSION(:), ALLOCATABLE :: atom1(:)
     INTEGER, DIMENSION(:), ALLOCATABLE :: atom4(:)
     INTEGER :: atom1_natoms
     INTEGER :: atom4_natoms

  END TYPE Dihedral_Atoms_to_Place_Class
  !****************************************************************************

  TYPE Bond_Participation_Class

     ! bondpart_list has dimensions of (MAXVAL(atoms),nspecies)
     ! The purpose of this class is to extract information as to how many
     ! bonds are shared by a given atom, the bond numbers and bonded atoms
     ! for these bonds. Atom and bond_num are allocatable arrays that can
     ! be assigned either MAXVAL(nbonds) or nbonds(ispecies)

     INTEGER :: nbonds
     INTEGER, DIMENSION(:), ALLOCATABLE :: atom
     INTEGER, DIMENSION(:), ALLOCATABLE :: bond_num
     
  END TYPE Bond_Participation_Class
  !****************************************************************************

  TYPE Angle_Participation_Class

     !angle_part_list has dimensions of (MAXVAL(natoms),nspecies)
     ! The purpose of this class is to obtain information on how many angles
     ! a given atom participates in, the identity of these angles and the
     ! position of atom in these angles, 1,2 or 3.

     INTEGER :: nangles
     INTEGER, DIMENSION(:), ALLOCATABLE :: which_angle
     INTEGER, DIMENSION(:), ALLOCATABLE :: position

  END TYPE Angle_Participation_Class
  !****************************************************************************

  TYPE Dihedral_Participation_Class

     ! dihedral_part_list has dimension of (MAXVAL(natoms),nspecies). This
     ! class defines the number of dihedral angles a given atom participates
     ! in. This information is stored in the variable ndihedrals. It also
     ! returns what these dihedrals are and what is the position of the atom
     ! in each of these dihedrals, e.g. 1,2,3 or 4.

     INTEGER :: ndihedrals
     INTEGER, DIMENSION(:), ALLOCATABLE :: which_dihedral
     INTEGER, DIMENSION(:), ALLOCATABLE :: position

  END TYPE Dihedral_Participation_Class


  ! Define a class that hold information on a box. length hold the box
  ! cell matrix, and length_inv its inverse in A and A^(-1). 
  ! basis_length is hte computed basis vector lengths
  ! cos_angle and angle are the vectors cos(alpha), cos(beta) and cos(gamma) 
  ! and alpha, beta, gamma.
  ! Face_distance is the distance in A between the faces of the box. 
  ! box_list has dimensions of the number of boxes.

 TYPE Box_Class
    CHARACTER(20) :: box_shape
    INTEGER :: int_box_shape
    REAL(DP), DIMENSION(3,3) :: length, length_inv, max_delta, hlength
    REAL(DP), DIMENSION(3) :: basis_length, cos_angle, angle, face_distance
    REAL(DP) :: volume, dv_max
 END TYPE Box_Class

  !****************************************************************************

 TYPE Angle_Probability_Class

    ! Will have dimension of (MAXVAL(nangles),nspecies)
    ! this class is used to hold the probabilty of observing a bond angle
    ! region for a given bond angle. There are three elements of this
    ! class. 
    ! prob : holds the probability of a given region. It is a vector of length
    !        nregions
    ! theta : holds the midpoint angle of a given region. It is a vector of
    !         length nregions
    ! nregions : number of regions for which there is a nonzero probability
    
    INTEGER :: nregions
    REAL(DP), DIMENSION(nregions) :: prob
    REAL(DP), DIMENSION(nregions) :: theta

 END TYPE Angle_Probability_Class
 !**************************************************************************************

 TYPE Bond_Length_Probability_Class

    ! Will have dimentions of (MAXVAL(nbonds),nspecies)
    ! The class is used to store the probability of observing a bond length region
    ! for a given bond length. There are three elements of this class 
    !
    ! prob : holds the probability of a given region. 
    ! length : bond length corresponding to the mid point of a region
    ! nregions : number of regions for which there is a nonzero probability
    ! delta_r : resolution of the probability distribution
    ! rlower : lower limit for the probability calculation

    INTEGER :: nregions
    REAL(DP) :: delta_r
    REAL(DP), DIMENSION(:), ALLOCATABLE :: prob
    REAL(DP), DIMENSION(:), ALLOCATABLE :: length
    REAL(DP) :: rlower

 END TYPE Bond_Length_Probability_Class
 !*******************************************************************************************

 TYPE Energy_Class

    ! This class will have dimensions of nbr_boxes
    ! It is defined to hold various energies of the system.
    ! 
    ! nonbond_inter_vdw : total intermolecular vdW interactions of molecules in a box
    ! nonbond_inter_q :   intermolecular charge interactions
    ! nonbond_intra_vdw : intramolecular nonbonded vdW interactions
    ! nonbond_intra_q   : intramolecular nonbonded charge interactions
    ! intra             : bonded intramolecular interactions
    ! bond              : contribution from bonds
    ! angle             : contribution from angles
    ! dihedral          : contribution from dihedrals
    ! ewald_reciprocal  : reciprocal space component of Ewald summation
    ! ewald_self        : Ewald self energy component
    ! total             : total energy of the system

    REAL(DP) :: inter_vdw, lrc, inter_q, intra_vdw, intra_q
    REAL(DP) :: intra, ewald_reciprocal, total
    REAL(DP) :: bond, angle, dihedral, improper
    REAL(DP) :: self
 END TYPE Energy_Class
 !------------------------------------------------------------------------------------------------

 TYPE MC_Moves_Class
    ! This class holds the information about number of trial MC moves attempted for
    ! each of the species in every box. The type has dimensions of (nspecies,nbr_boxes)
    
    ! displacement : COM translation move
    ! rotation : rotation around COM
    ! angle : angle distortion
    ! bond : bond fluctuation
    ! dihedral : dihedral rotation
    ! insertion : insertion of a molecule
    ! deletion : deletion of a molecule
    ! disp_atom : atom displacement move
    ! cpcalc : calculation of chemical potential
    ! cluster : move a cluster of molecules
    
    INTEGER :: displacement, rotation, angle, bond, dihedral, insertion, deletion, cluster
    INTEGER :: disp_atom, cpcalc, displacement_e, rotation_e

 END TYPE MC_Moves_Class
 !-------------------------------------------------------------------------------------------------

 TYPE Frag_Class
    ! This class has been defined to hold information about the fragments defined for species.
    ! The type has dimensions (MAXVAL(nfragments), nspecies)

    ! natoms : total number of atoms in a given fragment
    ! atoms(i) : ith atom of the fragment
    ! nconnect : number of connections for a fragment
    ! frag_connect(i) : fragment id of ith connection
    ! anchor : atom id with respect to which atoms in the fragment are located.
    ! type : fragment type (1, 2 etc)
    ! nconfig : total number of conformations in each fragment
    ! Total number of anchors and id of the anchors are stored in 'nanchors' 
    ! and 'anchor' array
    ! prob_ins : probability of inserting this fragment first
    ! cum_prob_ins : cumulative probablility of inserting this fragment first
    INTEGER :: natoms, nconnect, nanchors, type, nconfig
    INTEGER, DIMENSION(:), ALLOCATABLE :: anchor
    INTEGER, DIMENSION(:), ALLOCATABLE :: atoms
    INTEGER, DIMENSION(:), ALLOCATABLE :: frag_connect
    LOGICAL:: ring
    REAL(DP)::rcut_vdwsq, rcut_coulsq,alpha_ewald
    REAL(DP) :: prob_ins, cum_prob_ins

 END TYPE Frag_Class
!-------------------------------------------------------------------------------------------------

 TYPE Fragment_Bond_Class
    ! This class holds the information on the fragments involved in a fragment bond
    ! The type has dimensions of (MAXVAL(fragment_bonds), nspecies)

    ! frag_id of the two fragments connected via this bond
    INTEGER :: fragment1, fragment2
    ! probability of deleting fragment1 if this bond is cut
    ! don't need to store prob_del2, since it must be 1.O_DP - prob_del1
    REAL(DP) :: prob_del1

 END TYPE Fragment_Bond_Class

!-------------------------------------------------------------------------------------------------

 TYPE Library_Coords_Class
    REAL(DP) :: rxp, ryp, rzp
 END TYPE Library_Coords_Class

!-------------------------------------------------------------------------------------------------

 TYPE Energy_Fragment_Class
    REAL(DP), DIMENSION(:), ALLOCATABLE :: this_config_energy
 END TYPE Energy_Fragment_Class

!-------------------------------------------------------------------------------------------------

 TYPE Pressure_Class
    ! This class holds the pressure for each box

    ! Setpoint pressure, provided for constant pressure simulations
    REAL(DP) :: setpoint

    ! Computed pressure, kNT/V + virial
    REAL(DP) :: computed

    ! last calculation
    INTEGER :: last_calc

 END TYPE Pressure_Class

!-------------------------------------------------------------------------------------------------

END MODULE Type_Definitions
