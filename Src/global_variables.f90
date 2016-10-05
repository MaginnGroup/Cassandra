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

!********************************************************************************
MODULE Global_Variables
!********************************************************************************
 
  ! Written by: Ed Maginn
  ! Sept., 2007

  ! *** Used by ***
  ! atoms_to_place
  ! compute_cell_dimensions
  ! create_nonbond_table
  ! get_internal_coords
  ! make_config
  ! input_routines
  ! io_utilities
  ! main
  ! minimum_image_separation
  ! nvtmc_control
  ! participation
  ! random_generators   ONLY : iseed
  ! save_revert_coordinates
  ! translate

  ! Revision history:
  !
  !
  ! 08/02/13 : Beta release version
  ! 03/17/15 (JS) : lactivity defined for GCMC simulations
!*********************************************************************************

USE Type_Definitions

  SAVE

!*********************************************************************************
  ! This section contains global variables used by many routines during the run.

  CHARACTER(120) :: run_name
  CHARACTER(15), DIMENSION(:), ALLOCATABLE :: start_type
  CHARACTER(80) :: err_msg(10)

  ! error handling variables
  INTEGER :: AllocateStatus, OpenStatus, DeAllocateStatus

  ! Timing function
  CHARACTER(15) :: hostname,date,time,zone
  INTEGER, DIMENSION(8) :: values, begin_values,end_values

  ! Type of simulation to run:
  ! Choices: NVT_MC 
  CHARACTER(20) :: sim_type
  INTEGER :: int_sim_type
  INTEGER, PARAMETER :: sim_nvt = 0
  INTEGER, PARAMETER :: sim_nvt_min = 1
  INTEGER, PARAMETER :: sim_npt = 2
  INTEGER, PARAMETER :: sim_gemc = 3
  INTEGER, PARAMETER :: sim_gcmc = 4
  INTEGER, PARAMETER :: sim_frag = 5
  INTEGER, PARAMETER :: sim_ring = 6
  INTEGER, PARAMETER :: sim_gemc_npt = 7
  INTEGER, PARAMETER :: sim_gemc_ig = 8
  INTEGER, PARAMETER :: sim_mcf = 9
  LOGICAL :: timed_run, openmp_flag, en_flag, verbose_log
  CHARACTER(10) :: sim_length_units
  INTEGER :: steps_per_sweep

  ! The starting seed for the random generator
  ! Note iseed is used for generating points on random sphere for MCF_Gen sim type.
 INTEGER (KIND=8) :: iseed, iseed1, iseed3 

  ! Variables associated with the nonbond potential
  CHARACTER(15) :: mix_rule, run_type
  CHARACTER(15), DIMENSION(:), ALLOCATABLE :: vdw_style, charge_style, vdw_sum_style, charge_sum_style
  INTEGER :: int_mix_rule, int_run_type
  INTEGER, DIMENSION(:), ALLOCATABLE :: int_vdw_style, int_vdw_sum_style
  INTEGER, DIMENSION(:), ALLOCATABLE :: int_charge_style, int_charge_sum_style
  INTEGER, PARAMETER :: run_equil = 0
  INTEGER, PARAMETER :: run_prod = 1
  INTEGER, PARAMETER :: run_test = 2
  INTEGER, PARAMETER :: vdw_none = 0
  INTEGER, PARAMETER :: vdw_lj = 1
  INTEGER, PARAMETER :: vdw_cut = 2
  INTEGER, PARAMETER :: vdw_cut_shift = 3
  INTEGER, PARAMETER :: vdw_cut_tail = 4
  INTEGER, PARAMETER :: vdw_minimum = 5
  INTEGER, PARAMETER :: vdw_charmm = 6
  INTEGER, PARAMETER :: vdw_cut_switch = 7
  INTEGER, PARAMETER :: vdw_mie = 8

  INTEGER, PARAMETER :: charge_none = 0
  INTEGER, PARAMETER :: charge_coul = 1
  INTEGER, PARAMETER :: charge_cut = 2
  INTEGER, PARAMETER :: charge_ewald = 3
  INTEGER, PARAMETER :: charge_minimum = 4
  INTEGER, PARAMETER :: charge_dsf = 5

  REAL(DP), DIMENSION(:), ALLOCATABLE :: rcut_cbmc 
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rcut_vdw, rcut_coul, ron_charmm, roff_charmm, rcut_max
  REAL(DP), DIMENSION(:), ALLOCATABLE :: ron_switch, roff_switch, roff_switch_sq, switch_factor1
  REAL(DP), DIMENSION(:), ALLOCATABLE :: switch_factor2, ron_switch_sq
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rcut_vdwsq, rcut_coulsq, ron_charmmsq, roff_charmmsq
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rcut9, rcut3
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rcut_vdw3, rcut_vdw6
  REAL(DP) :: edens_cut, rcut_clus, rcut_low, rcut_lowsq
  LOGICAL, DIMENSION(:), ALLOCATABLE :: l_half_len_cutoff

 ! Mixing Rules variables :
  CHARACTER(40), DIMENSION(:,:), ALLOCATABLE :: vdw_interaction_table
  INTEGER, DIMENSION(:,:), ALLOCATABLE ::vdw_int_table
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: vdw_param1_table
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: vdw_param2_table, vdw_param3_table
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: vdw_param4_table, vdw_param5_table
  REAL(DP), DIMENSION(:), ALLOCATABLE :: alpha_ewald, h_ewald_cut
  REAL(DP), DIMENSION(:), ALLOCATABLE :: alphal_ewald
  REAL(DP), DIMENSION(:), ALLOCATABLE :: ewald_p_sqrt, ewald_p
  
 
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nint_beads

  ! Intramolecular exclusion variables (1-2, 1-3, 1-4 exclusions/scaling)
  ! and the scaling to use for all other intramolecular terms.
  REAL(DP), DIMENSION(:), ALLOCATABLE :: scale_1_2_vdw, scale_1_3_vdw, scale_1_4_vdw, scale_1_N_vdw
  REAL(DP), DIMENSION(:), ALLOCATABLE :: scale_1_2_charge, scale_1_3_charge, scale_1_4_charge, scale_1_N_charge

  ! Dimensions (maxatomtype,maxatomtype,nspecies)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: vdw_intra_scale, charge_intra_scale
  LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: l_bonded

  ! How many simulation boxes we have. 
  INTEGER :: nbr_boxes
  INTEGER, PARAMETER :: int_cubic = 0
  INTEGER, PARAMETER :: int_ortho = 1
  INTEGER, PARAMETER :: int_cell = 2

! Atom placement variables

  ! Univ vectors pointing to 6,18,26, and 50 uniform points.
  REAL(DP), DIMENSION(3,50) :: sphere_vec

  ! Imsl

  INTEGER, PARAMETER :: irule = 3
  
 !***************************************************************
  !Conversion factors and constants

  REAL(DP), PARAMETER :: PI=3.1415926536_DP
  REAL(DP), PARAMETER :: twoPI = 6.2831853072_DP
  REAL(DP), PARAMETER :: rootPI = 1.7724538509_DP

  !KBOLTZ is Boltzmann's constant in atomic units amu A^2 / (K ps^2)
  REAL(DP), PARAMETER :: kboltz = 0.8314472_DP

  !H_PLANK is Plank's constant in atomic units amu A^2 / ps^3
  REAL(DP), PARAMETER :: h_plank = 39.9031268605_DP

  ! The value of the fundamental electron charge squared and divided by
  ! 4*pi*epsilon0, where epsilon0 is the vacuum permittivity.  This value
  ! alllows a coulombic potential of the form (qi*qj/rij) to be used.  The
  ! full form is (qi*qj*e^2/(4*pi*epsilon0*rij)).  To simplify, the extra
  ! constant terms are subsumed into the single constant described above.
  ! The units of charge factor are:   amu A^3 / ps^2
  REAL(DP), PARAMETER :: charge_factor = 138935.4558_DP

  !Factor to convert atomic pressure (amu / (A ps^2) ) to bar
  REAL(DP), PARAMETER :: atomic_to_bar = 166.054_DP

  !Factor to convert atomic energy (amu A^2/ ps^2) to kJ/mol
  REAL(DP), PARAMETER :: atomic_to_kJmol = 0.01_DP

  !Factor to conver atomic energy to K
  REAL(DP), PARAMETER :: atomic_to_K = 1.2027221933_DP

  !Factor to convert kJ/mol to atomic energy (amu A^2/ ps^2) 
  REAL(DP), PARAMETER :: kjmol_to_atomic = 100.0_DP

  !Factor to convert energy in (eV) to atomic energy (amu A^2/ps^2)
  REAL(DP), PARAMETER :: eV_to_atomic = 9648.53082148_DP

  !Factor to convert energy in kJ/mol to kcal/mol
  REAL(DP), PARAMETER :: kJmol_to_kcalmol = 0.239005736_DP

  !Factor to convert mass_density in g/mol/A^3 to kg/m^3
  REAL(DP), PARAMETER :: atomic_to_kgm3 = 1660.54_DP

  ! small number for comparison
  REAL(DP), PARAMETER :: tiny_number  = 0.0000001_DP
  REAL(DP), PARAMETER :: small_number = 0.00001_DP

  ! Upper limit to prevent underflow in exp(-beta*energy)
  ! DEXP(-708.)= 3.307553003638408E-308
  ! DEXP(-709.)= 0.000000000000000E+000
  REAL(DP), PARAMETER :: max_kBT = 708.0_DP

  ! IMSL error bounds
  REAL(DP), PARAMETER :: errabs = 0.0_DP
  REAL(DP), PARAMETER :: errel = 1.0E-5_DP

  ! Parameter identifying number of trials

  INTEGER :: kappa_ins, kappa_rot, kappa_dih

  ! Parameters identifying move in Ewald calculations
  
  INTEGER, PARAMETER :: int_insertion = 0
  INTEGER, PARAMETER :: int_deletion = 1
  INTEGER, PARAMETER :: int_translation = 2
  INTEGER, PARAMETER :: int_rotation = 3
  INTEGER, PARAMETER :: int_intra = 4

  ! Parameter for species type 

  INTEGER, PARAMETER :: int_sorbate = 0
  INTEGER, PARAMETER :: int_solvent = 1

  ! Parameters for insertion type
  INTEGER, PARAMETER :: int_noinsert = -1
  INTEGER, PARAMETER :: int_random = 0
  INTEGER, PARAMETER :: int_igas = 1

  ! Parameters for dihedral and improper type

  INTEGER, PARAMETER :: int_none = 0
  INTEGER, PARAMETER :: int_opls = 1
  INTEGER, PARAMETER :: int_charmm = 2
  INTEGER, PARAMETER :: int_harmonic = 3
  INTEGER, PARAMETER :: int_cvff = 4
  INTEGER, PARAMETER :: int_amber = 5  

  ! Define integers for molecule type
  INTEGER, PARAMETER :: int_noexist = -1
  INTEGER, PARAMETER :: int_normal = 0
  INTEGER, PARAMETER :: int_fractional = 1
  INTEGER, PARAMETER :: int_fixed = 2

  !**********************************************************************************
  ! thermodynamic state point variables
 
  REAL(DP),DIMENSION(:),ALLOCATABLE,TARGET :: temperature, beta
  TYPE(Pressure_Class), DIMENSION(:), ALLOCATABLE, TARGET :: pressure
  LOGICAL :: need_pressure
  
  ! **********************************************************************************
  ! system size integers used in memory allocation.
  ! Number of species, molecules, atoms, bonds, angles, dihedrals and impropers should 
  ! be kept as independent arrays  

  INTEGER :: nspecies, nspec_insert
  INTEGER, DIMENSION(:), ALLOCATABLE :: n_igas, n_igas_update, n_igas_moves, nzovero ! integers for ideal gas reservoir
  LOGICAL :: first_res_update, igas_flag
  LOGICAL, DIMENSION(:), ALLOCATABLE :: zig_calc
  INTEGER, DIMENSION(:), ALLOCATABLE :: max_molecules, natoms, nmol_start, nring_atoms, nexo_atoms
  INTEGER, DIMENSION(:), ALLOCATABLE :: nbonds, nangles, nangles_fixed
  INTEGER, DIMENSION(:), ALLOCATABLE :: ndihedrals, nimpropers
  INTEGER, DIMENSION(:), ALLOCATABLE :: nfragments, fragment_bonds

  ! array to hold the total number of molecules of each species in a given box

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nmols
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nmols_to_make, nmols_to_read

  ! array to hold ring atom ids and exo atom ids for a fragment
  ! will have (MAXVAL(natoms), nspecies) dimensions
  
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ring_atom_ids, exo_atom_ids

  ! force field parameter numbers - set in Input_Routines.
  ! Keep track of the number of parameters each species has.

  ! number of unique atom types
  INTEGER :: nbr_atomtypes
  INTEGER, DIMENSION(:), ALLOCATABLE :: nbeads_in, nbeads_out
  INTEGER, DIMENSION(:), ALLOCATABLE :: nbeadsfrac_in

  ! number of ideal gas particles in the intermediate box
  INTEGER :: igas_num

  ! array containing name of each atom type with idex = atomtype number.
  ! It is set and allocated to size nbr_atomtypes in Create_Nonbond_Table
  CHARACTER(6), DIMENSION(:), ALLOCATABLE :: atom_type_list

  ! Number of parameters required for various potential functions.
  INTEGER, DIMENSION(:), ALLOCATABLE :: nbr_bond_params, nbr_angle_params 
  INTEGER, DIMENSION(:), ALLOCATABLE :: nbr_improper_params, nbr_vdw_params
  INTEGER, DIMENSION(:), ALLOCATABLE :: nbr_dihedral_params

  ! Information of the position line where starts the coordinates storage of
  ! each fragment type

   INTEGER, DIMENSION(:), ALLOCATABLE :: frag_position_library

  ! **********************************************************************************
  ! Basic data structures are in the form of arrays. Derived from 

  ! type classes defined in Type_Definitions.

  ! Array with dimension (nspecies)
  TYPE(Species_Class), DIMENSION(:), ALLOCATABLE, TARGET :: species_list
    
  ! Array with dimensions (max_molecules,nspecies)
  TYPE(Molecule_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: molecule_list
  TYPE(Molecule_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: molecule_list_igas

  ! Array with dimensions (coordinate index, max_molecules, nspecies)
  TYPE(Internal_Coord_Class), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: internal_coord_list

  ! Array with dimension (nbonds)
  TYPE(Internal_Coord_Class_Old), DIMENSION(:), ALLOCATABLE, TARGET :: internal_coord_list_old

  ! Array with dimensions (natoms, max_molecules, nspecies)
  TYPE(Atom_Class), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: atom_list
  TYPE(Atom_Class), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: atom_list_igas
  ! Array with dimension (natoms,1,nspecies) Describes positions of starting gemoetry if a configuration is to be generated
  TYPE(Atom_Class), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: init_list

  ! Array with dimensions (natoms, nspecies)
  TYPE(Nonbond_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: nonbond_list
  
  ! Array with dimensions (nbonds, nspecies)
  TYPE(Bond_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: bond_list
  
  ! Array with dimensions (nangles, nspecies)
  TYPE(Angle_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: angle_list
  
  ! Array with dimensions (ndihedrals, nspecies)
  TYPE(Dihedral_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: dihedral_list
  
  ! Array with dimensions (nimpropers, nspecies)
  TYPE(Improper_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: improper_list

  ! Array with dimension (MAXVAL(natoms), nspecies)
  TYPE(Bond_Participation_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: bondpart_list

  ! Array with dimension (MAXVAL(natoms), nspecies)
  TYPE(Angle_Participation_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: angle_part_list

  ! Array with dimension (MAXVAL(natoms),nspecies)
  TYPE(Dihedral_Participation_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: dihedral_part_list

  ! Array with dimension (MAXVAL(nbonds), nspecies)
  TYPE(Bond_Atoms_To_Place_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: bond_atoms_to_place_list

  ! Array with dimension (nangles, nspecies)
  TYPE(Angle_Atoms_To_Place_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: angle_atoms_to_place_list
  
  ! Array with dimension (ndihedrals, nspecies)

  TYPE(Dihedral_Atoms_To_Place_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: dihedral_atoms_to_place_list

  ! Array with box info, dimensions (nbr_boxes)
  TYPE(Box_Class), DIMENSION(:), ALLOCATABLE, TARGET :: box_list

  ! Array with fragment info, dimension (nfragments, nspecies)
  TYPE(Frag_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: frag_list

  ! Array with fragment bond info
  TYPE(Fragment_Bond_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: fragment_bond_list

  ! Array for storing coordinates of fragments
  !TYPE(Frag_Library_Class), DIMENSION(:), ALLOCATABLE :: frag_library
  TYPE(Library_Coords_Class), DIMENSION(:), ALLOCATABLE :: library_coords

  ! Array for storing the energy of each configuration of each fragment
  ! nrg_frag has dimension (number of fragment )
  TYPE(Energy_Fragment_Class), DIMENSION(:), ALLOCATABLE, TARGET :: nrg_frag
  



  ! **********************************************************************************

  ! Linked list for open ensemble simulations, will have dimensions of (MAXVAL(max_molecules),nspecies)

  INTEGER, DIMENSION(:,:,:), ALLOCATABLE, TARGET :: locate

  ! Array with angle probability info with dimension (MAXVAL(nangles),nspecies)
  
  TYPE(Angle_Probability_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: ang_prob

  ! Array with bond probability info with dimension (MAXVAL(nbonds),nspecies)
  TYPE(Bond_Length_Probability_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: bond_length_prob

  ! Bond probability cutoff
  REAL(DP) :: bond_probability_cutoff

  !**********************************************************************************************************
  ! Will have dimension of nbr_boxes
  TYPE(Energy_Class), DIMENSION(:), ALLOCATABLE, TARGET :: energy, virial
  TYPE(Energy_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: ac_energy
!  TYPE(Energy_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: ac_virial
  
  ! Will have dimension (MAXVAL(max_molecules))
  TYPE(Energy_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: energy_igas

  ! Accumulators for thermodynamic averages,

  ! will have dimensions of nbr_boxes
  REAL(DP), DIMENSION(:,:),ALLOCATABLE,TARGET :: ac_volume, ac_enthalpy, ac_pressure, ac_mass_density
  ! will have dimension of (nspecies,nbr_boxes)
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: ac_density, ac_nmols

  LOGICAL :: block_average

  ! The following variables are defined for Ewald calculations

  ! nvecs will have dimensions of nbr_boxes
  INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: nvecs

  INTEGER, PARAMETER  :: maxk = 100000
  
  ! Dimensions of (maxk, nbr_boxes)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, TARGET :: hx, hy, hz, hsq, Cn

  ! the following arrays will have dimensions of (MAXVAL(nvecs),nbr_boxes)

  REAL(DP), DIMENSION(:,:), ALLOCATABLE, TARGET :: cos_sum, sin_sum, cos_sum_old, sin_sum_old
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, TARGET :: cos_sum_start, sin_sum_start

  !*********************************************************************************************************
  ! Information on trial and probabilities of trial moves

  ! Will have dimensions of (nspecies,nbr_boxes)
  TYPE(MC_Moves_Class), DIMENSION(:,:), ALLOCATABLE,TARGET :: ntrials, nsuccess, nequil_success
  ! Variables associated with regrowth trials dimensions ( MAXVAL(nfragments), nbr_species)
  INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET :: regrowth_trials, regrowth_success

  INTEGER :: nupdate


  ! Will have dimension of nbr_boxes
  INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: nvolumes, nvol_success, ivol_success, tot_trials 
  INTEGER :: nvol_update

  ! individual move probability
  REAL(DP) :: prob_trans, prob_rot, prob_torsion, prob_volume, prob_angle, prob_insertion
  REAL(DP) :: prob_deletion, prob_swap, prob_regrowth, prob_ring, prob_atom_displacement
  REAL(DP), DIMENSION(:), ALLOCATABLE :: prob_rot_species
  REAL(DP), DIMENSION(:), ALLOCATABLE :: prob_swap_species, cum_prob_swap_species
  REAL(DP), DIMENSION(:), ALLOCATABLE :: prob_swap_from_box, cum_prob_swap_from_box

  LOGICAL :: l_prob_swap_species, l_prob_swap_from_box

  LOGICAL :: f_dv, f_vratio

  REAL(DP) :: omega_max, disp_max, delta_cos_max, delta_phi_max
  REAL(DP), DIMENSION(:), ALLOCATABLE, TARGET :: prob_species_trans, prob_species_rotate
  REAL(DP), DIMENSION(:), ALLOCATABLE::prob_growth_species ! dimension nspecies
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: max_disp, max_rot

  REAL(DP) :: cut_trans, cut_rot, cut_torsion, cut_volume, cut_angle, cut_insertion, cut_deletion
  REAL(DP) :: cut_swap, cut_regrowth, cut_ring, cut_atom_displacement, cut_lambda
 
  !*********************************************************************************************************

  ! Timing information
  ! Initial, current and final number of steps
  INTEGER :: i_mcstep, initial_mcstep, n_mcsteps, n_equilsteps, iblock
  ! Information on the output of data
  INTEGER :: nthermo_freq, ncoord_freq, block_avg_freq, nbr_blocks
  REAL(DP) :: data_points_per_block
 
  INTEGER,DIMENSION(:),ALLOCATABLE :: nbr_prop_files

  ! Number of properties per file, will have dimension of nbr_prop_files
  
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: prop_per_file
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: first_open

  LOGICAL :: cpcollect  !  logical determining if the chemical potential info is collected

  LOGICAL :: accept

  LOGICAL :: cbmc_flag, del_flag, phi_Flag, angle_Flag, imp_Flag

  ! Some variables for reaction Monte Carlo
  LOGICAL, DIMENSION(:), ALLOCATABLE :: has_charge
  LOGICAL :: get_fragorder, l_check

  INTEGER :: imreplace, isreplace


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! variables for the neighbor list
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  LOGICAL, ALLOCATABLE :: l_cubic(:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! timing functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  INTEGER :: num_moves, count_n, count_nn, count_cell

  INTEGER, PARAMETER :: imove_trans = 0
  INTEGER, PARAMETER :: imove_rot = 1
  INTEGER, PARAMETER :: imove_dihedral = 2
  INTEGER, PARAMETER :: imove_angle = 3
  INTEGER, PARAMETER :: imove_volume = 4
  INTEGER, PARAMETER :: imove_insert = 5
  INTEGER, PARAMETER :: imove_swap = 6
  INTEGER, PARAMETER :: imove_delete = 7
  INTEGER, PARAMETER :: imove_regrowth = 8
  INTEGER, PARAMETER :: imove_check = 9
  INTEGER, PARAMETER :: imove_atom_displacement = 10


  REAL(DP) :: time_s, time_e
  REAL(DP) :: start_time, tot_time
  REAL(DP), DIMENSION(0:14)  :: movetime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: W_tensor_charge
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: W_tensor_recip
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: W_tensor_vdw
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: W_tensor_total
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: W_tensor_elec
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: pressure_tensor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Energy check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  LOGICAL :: echeck_flag
  INTEGER :: iecheck

!!!!! Pair energy arrays. These arrays hold interaction energies between pairs of molecules !!!!!

  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: pair_nrg_vdw, pair_nrg_qq
  REAL(DP) :: copy_time, recip_time

  ! cos_mol and sin_mol arrays hold k space vectors for each molecule
  ! dimensions == (SUM(max_molecules), MAX(nvecs))
  REAL(DP), ALLOCATABLE :: cos_mol(:,:) , sin_mol(:,:)
  LOGICAL :: l_pair_nrg

  REAL(DP) pacc
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: chpot, chpotid


!!!! Zeolite variables
REAL(DP), ALLOCATABLE, DIMENSION(:) :: x_lat, y_lat, z_lat
INTEGER :: n_lat_atoms

!!! Pair_Nrg_Variables
REAL(DP), ALLOCATABLE :: pair_vdw_temp(:), pair_qq_temp(:)

!!!! DSF variables
REAL(DP), ALLOCATABLE, DIMENSION(:) :: alpha_dsf
REAL(DP), ALLOCATABLE, DIMENSION(:) :: dsf_factor1, dsf_factor2
  
END MODULE Global_Variables

