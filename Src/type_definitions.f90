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

  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

  !Create a double precision definition, to make numbers machine independent
  INTEGER, PARAMETER :: SP = REAL32
  INTEGER, PARAMETER :: DP = REAL64

  ! Define REAL type infinity using hexadecimal literal constants
  REAL(REAL32), PARAMETER :: infinity_sp = REAL(Z'7F800000',REAL32)
  REAL(REAL64), PARAMETER :: infinity_dp = REAL(Z'7FF0000000000000',REAL64)

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


  REAL(DP), PARAMETER :: PI=3.1415926536_DP
  REAL(DP), PARAMETER :: twoPI = 6.2831853072_DP
  REAL(DP), PARAMETER :: rootPI = 1.7724538509_DP
  REAL(DP), PARAMETER :: halfPI = 0.5_DP*PI
  REAL(DP), PARAMETER :: threehalfPI_DP = 3.0_DP*halfPI
  REAL(SP), PARAMETER :: twoPI_SP = REAL(twoPI,SP)
  REAL(SP), PARAMETER :: PI_SP = REAL(PI,SP)

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
     LOGICAL :: L_Coul_CBMC = .TRUE.
     ! N.B. natoms, max_molecules, etc. are in a separate arrays
     ! NR: for insertion style
     LOGICAL :: lcom


     !!!Widom insertions
     LOGICAL, DIMENSION(:), ALLOCATABLE :: test_particle
     INTEGER (KIND=INT64), DIMENSION(:), ALLOCATABLE :: insertions_in_step, widom_interval
     REAL(DP), DIMENSION(:), ALLOCATABLE :: widom_sum

     !! Is this a solute?  Is this a solvent species?
     LOGICAL :: l_solute, l_solvent, l_wsolute

     !!Atompair_nrg_table index bases
     INTEGER :: solute_base, solvent_base, wsolute_base

     ! Pair energy array index base
     INTEGER :: superlocate_base

     ! # of RB dihedrals, # of dihedrals with nonzero energy, # of dihedrals before stacked dihedrals were combined
     INTEGER :: ndihedrals_rb, ndihedrals_energetic, ndihedrals_uncombined

     ! CBMC biasing info
     INTEGER :: kappa_ins = 0, kappa_ins_pad8, kappa_ins_pad64
     INTEGER :: kappa_dih = 0, kappa_dih_pad8, kappa_dih_pad32
     INTEGER :: kappa_rot = 0
     INTEGER :: nfragments
     REAL(DP) :: theta_step, log_kappa_ins, log_kappa_rot, ln_pbias_dih_const
     REAL(SP) :: theta_step_sp
     LOGICAL :: need_kappa_ins = .FALSE., need_kappa_dih = .FALSE.

     REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sincos_lintheta_dp
     REAL(SP), DIMENSION(:,:), ALLOCATABLE :: sincos_lintheta_sp

     CONTAINS
             PROCEDURE setup_CBMC_kappas => Setup_Species_kappas
             PROCEDURE write_CBMC_kappas => Write_Species_kappas

  END TYPE Species_Class
  !****************************************************************************




  !****************************************************************************
  TYPE Molecule_Class

     ! The molecule list will have dimensions (max_molecules,nspecies)

     REAL(DP), DIMENSION(4) :: rcom, rcom_old

     ! What kind of molecule is this? normal, fractional, fixed, etc.
     ! Note that the following integers will be defined for the type
     ! int_normal = 0
     ! int_fractional = 1
     ! int_fixed = 2
     INTEGER :: molecule_type
     INTEGER :: rx_num

     ! for open system and multi-box simulations (GEMC, parallel tempering)
     INTEGER :: which_box
     LOGICAL :: live

     REAL(DP) :: unused_dummy1, unused_dummy2 ! added for padding derived type

     ! also include com information for each of the molecules.
     ! com and euler angles refer to the x,y and z com coordinates
     ! and 1, 2 and 3 euler angles of the molecule. The suffix
     ! old denotes old coordinates.
     ! frac is the fractional scaling parameter for the molecule
     ! xcom, ycom, and zcom are now rcom(1), rcom(2), and rcom(3), respectively
     ! likewise for with the "_old" suffix
     REAL(DP)  :: euler1, euler2, euler3
     REAL(DP)  :: euler1_old,euler2_old,euler3_old
     REAL(DP)  :: frac
     ! This variable records the maximum distance of any psuedo atom from its
     ! COM. This is used to speed up energy calculations.

     ! max_dcom is now rcom(4) and max_dcom_old is now rcom_old(4)

     REAL(DP) :: min_dcom



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

     REAL(DP), DIMENSION(3) :: rp, rp_old

     !REAL(DP) :: rxp, ryp, rzp
     REAL(DP) :: rxp_nls, ryp_nls, rzp_nls  ! The starting positions for the neighbor list
     !REAL(DP) :: rxp_old, ryp_old, rzp_old
     INTEGER :: ci(3), ci_cbmc(3), ci_full(3) ! the integer coordinates of the cell containing this atom
     LOGICAL :: exist

  END TYPE Atom_Class
  !****************************************************************************

  TYPE Atom256
          REAL(DP) :: rxp, ryp, rzp
          LOGICAL(8) :: exist
  END TYPE Atom256

  TYPE VdW256
          REAL(DP) :: p1, p2, p3, p4
  END TYPE VdW256



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

     CHARACTER(6) :: element
     CHARACTER(23) :: atom_name

     REAL(DP) :: mass, charge
     INTEGER :: atom_type_number

     LOGICAL :: ring_atom, multiring_atom

  END TYPE Nonbond_Class
  !****************************************************************************



  !****************************************************************************
  TYPE Bond_Class

     ! bond list has dimensions (MAXVAL(nbonds), nspecies)
     INTEGER :: atom(2), int_bond_type
     REAL(DP), DIMENSION(max_bond_params) :: bond_param
     CHARACTER(20) :: bond_potential_type

  END TYPE Bond_Class
  !****************************************************************************



  !****************************************************************************
  TYPE Angle_Class

     ! angle list has dimensions (MAXVAL(nangles,nspecies)

     ! 1 - 2 - 3 forms an angle with 2 at the apex.
     INTEGER :: atom(3)

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

     INTEGER :: atom(4)
     ! RB torsion series constants
     REAL(DP) :: rb_c(0:5)
     REAL(DP), DIMENSION(max_dihedral_params) :: dihedral_param
     REAL(SP) :: rb_c_sp(0:5)
     REAL(SP), DIMENSION(max_dihedral_params) :: dihedral_param_sp
     CHARACTER(20) :: dihedral_potential_type
     INTEGER :: int_dipot_type
     ! Flag to tell whether dihedral is formatted as RB torsion
     LOGICAL :: l_rb_formatted
     CONTAINS
             PROCEDURE :: SP_Convert => Convert_Dihedral_DP_to_SP
             PROCEDURE :: Init => Initialize_Dihedral_Class

  END TYPE Dihedral_Class
  !****************************************************************************



  !****************************************************************************
  TYPE Improper_Class

     ! improper list has dimensions (MAXVAL(nimpropers), nspecies)

     ! Describes an improper dihedral: atom 1 is the central atom, and atoms 2
     ! through 4 are attached to atom 1, parameters are stored in an array

     INTEGER :: atom(4)
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

  !****************************************************************************

 TYPE Box_Class

  ! Define a class that hold information on a box. length hold the box
  ! cell matrix, and length_inv its inverse in A and A^(-1).
  ! basis_length is hte computed basis vector lengths
  ! cos_angle and angle are the vectors cos(alpha), cos(beta) and cos(gamma)
  ! and alpha, beta, gamma.
  ! Face_distance is the distance in A between the faces of the box.
  ! box_list has dimensions of the number of boxes.

    CHARACTER(20) :: box_shape
    INTEGER :: int_box_shape
    REAL(DP), DIMENSION(3,3) :: length, length_inv, max_delta, hlength
    REAL(DP), DIMENSION(3) :: basis_length, cos_angle, angle, face_distance
    REAL(DP) :: volume, dv_max

    REAL(SP), DIMENSION(3) :: sp_diag_length, cell_length_recip, real_length_cells
    REAL(SP), DIMENSION(3) :: cell_H_diag
    REAL(DP), DIMENSION(3) :: cell_face_distance, cell_xyzortho_bbox_length
    REAL(DP), DIMENSION(3,3) :: cell_length_inv, cell_H_dp
    REAL(SP), DIMENSION(3,3) :: cell_H_sp
    INTEGER, DIMENSION(3) :: length_cells, sectorbound

  ! Inner shape is used to define a limited region into which molecules can be
  ! inserted
    INTEGER :: int_inner_shape
    REAL(DP) :: inner_volume, inner_radius, inner_radius2, inner_zmax, inner_zmin

    INTEGER :: border_thickness(3)
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: cbmc_cell_mask

    INTEGER, DIMENSION(3) :: length_bitcells, setbit_extent
    REAL(SP), DIMENSION(3) :: bit_cell_length_recip, real_length_bitcells
    INTEGER(4), DIMENSION(:), ALLOCATABLE :: bitcell_int32_vec
    INTEGER, DIMENSION(2:3) :: bitcell_dimfactor
    REAL(DP) :: ideal_bitcell_length, rcut_low_max
    REAL(DP), DIMENSION(3) :: bitcell_face_distance, bitcell_face_distance_recip

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
    ! total             : total energy of the system
    ! intra             : intramolecular interactions
    !   bond            : contribution from bonds
    !   angle           : contribution from angles
    !   dihedral        : contribution from dihedrals
    !   improper        : contribution from impropers
    !   intra_vdw       : intramolecular nonbonded vdW interactions
    !   intra_q         : intramolecular nonbonded charge interactions
    ! inter             : intermolecular interactions
    !   inter_vdw       : total intermolecular vdW interactions of molecules in a box
    !   lrc             : long range tail correction to vdw energy
    !   inter_q         : intermolecular charge interactions
    !   reciprocal      : reciprocal space component of Ewald summation
    !   self            : Ewald self energy component

    REAL(DP) :: total
    REAL(DP) :: intra, inter
    REAL(DP) :: bond, angle, dihedral, improper, intra_vdw, intra_q
    REAL(DP) :: inter_vdw, lrc, inter_q, reciprocal, self
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

    INTEGER (KIND=INT64) :: displacement, rotation, angle, bond, dihedral, insertion, deletion, switch, cluster
    INTEGER (KIND=INT64) :: disp_atom, cpcalc, displacement_e, rotation_e

    ! widom insertions technically aren't mc moves but it's convenient to count them in ntrials
    INTEGER (KIND=INT64) :: widom


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
    REAL(DP) :: rp(3)
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

	! Ideal Pressure
	REAL(DP) :: ideal

    ! last calculation
    INTEGER :: last_calc

 END TYPE Pressure_Class

!-------------------------------------------------------------------------------------------------

  TYPE Rotation_Class
     !This class holds info in order to calculate a rotational bias

     REAL(DP) :: angle1
     REAL(DP) :: angle2
     REAL(DP) :: angle3

     ! flag when overlap triggered, avoid further computations with energy
     LOGICAL :: overlap

     !variables to calculate weights for golden sampling
     REAL(DP) :: dE
     REAL(DP) :: exp_dE
     REAL(DP) :: exp_dE_ratio
  END TYPE Rotation_Class

  PRIVATE Convert_Dihedral_DP_to_SP, Initialize_Dihedral_Class, Setup_Species_Kappas, Write_Species_Kappas

  CONTAINS
          ELEMENTAL SUBROUTINE Convert_Dihedral_DP_to_SP(this)
                  CLASS(Dihedral_Class), INTENT(INOUT) :: this
                  this%rb_c_sp = REAL(this%rb_c,SP)
                  this%dihedral_param_sp = REAL(this%dihedral_param,SP)
          END SUBROUTINE Convert_Dihedral_DP_to_SP
          ELEMENTAL SUBROUTINE Initialize_Dihedral_Class(this)
                  CLASS(Dihedral_Class), INTENT(INOUT) :: this
                  this%atom = 0
                  this%rb_c = 0.0_DP
                  this%dihedral_param = 0.0_DP
                  this%rb_c_sp = 0.0 ! Maybe this should be initalized to NaN, so there's a clear error if it wasn't converted first
                  this%dihedral_param_sp = 0.0
                  this%dihedral_potential_type = ""
                  this%int_dipot_type = 0
                  this%l_rb_formatted = .FALSE.
          END SUBROUTINE Initialize_Dihedral_Class
          ELEMENTAL SUBROUTINE Setup_Species_Kappas(this)
                  CLASS(Species_Class), INTENT(INOUT) :: this
                  INTEGER :: i
                  REAL(DP) :: theta
                  this%kappa_ins_pad8 = IAND(this%kappa_ins+7,NOT(7))
                  this%kappa_ins_pad64 = IAND(this%kappa_ins+63,NOT(63))
                  this%kappa_dih_pad8 = IAND(this%kappa_dih+7,NOT(7))
                  this%kappa_dih_pad32 = IAND(this%kappa_dih+31,NOT(31))
                  IF (this%need_kappa_ins) this%log_kappa_ins = LOG(REAL(this%kappa_ins,DP))
                  IF (this%kappa_rot > 0) this%log_kappa_rot = LOG(REAL(this%kappa_rot,DP))
                  IF (this%need_kappa_dih) THEN
                          this%ln_pbias_dih_const = REAL(this%nfragments-1,DP)*LOG(REAL(this%kappa_dih,DP))
                          IF (ALLOCATED(this%sincos_lintheta_dp)) DEALLOCATE(this%sincos_lintheta_dp)
                          IF (ALLOCATED(this%sincos_lintheta_sp)) DEALLOCATE(this%sincos_lintheta_sp)
                          ALLOCATE(this%sincos_lintheta_dp(this%kappa_dih_pad8,2))
                          ALLOCATE(this%sincos_lintheta_sp(this%kappa_dih_pad8,2))
                          this%theta_step = twoPI/this%kappa_dih
                          this%theta_step_sp = REAL(this%theta_step,SP)
                          !DIR$ VECTOR ALIGNED
                          DO i = 0, this%kappa_dih_pad8-1
                                theta = i*this%theta_step
                                this%sincos_lintheta_dp(i+1,1) = SIN(theta)
                                this%sincos_lintheta_dp(i+1,2) = COS(theta)
                          END DO
                          this%sincos_lintheta_sp = REAL(this%sincos_lintheta_dp,SP)
                  END IF
          END SUBROUTINE Setup_Species_Kappas
          SUBROUTINE Write_Species_Kappas(this,outputunit)
                  CLASS(Species_Class), INTENT(INOUT) :: this
                  INTEGER, INTENT(IN) :: outputunit
                  IF (this%need_kappa_ins) THEN
                          WRITE(outputunit,'(X,A,T35,I12)') 'Kappa for first fragment insertion ', this%kappa_ins
                  END IF
                  IF (this%need_kappa_dih) THEN
                          WRITE(outputunit,'(X,A,T35,I12)') 'Kappa for dihedral selection ', this%kappa_dih
                  END IF
          END SUBROUTINE Write_Species_Kappas

END MODULE Type_Definitions
