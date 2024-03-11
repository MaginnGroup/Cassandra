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

USE ISO_FORTRAN_ENV
USE Type_Definitions
!$ USE OMP_LIB

  SAVE

!*********************************************************************************
  ! This section contains global variables used by many routines during the run.

  CHARACTER(120) :: run_name
  CHARACTER(15), DIMENSION(:), ALLOCATABLE :: start_type
  CHARACTER(80) :: err_msg(10)

  ! error handling variables
  INTEGER :: AllocateStatus, OpenStatus, DeAllocateStatus
  !$OMP THREADPRIVATE(AllocateStatus, OpenStatus, DeAllocateStatus)

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
  INTEGER, PARAMETER :: sim_pregen = 10
  LOGICAL :: timed_run, openmp_flag, en_flag, verbose_log, input_is_logfile, open_mc_flag
  CHARACTER(10) :: sim_length_units
  INTEGER (KIND=INT64):: steps_per_sweep

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
  LOGICAL :: change_to_production = .false.
  INTEGER (KIND=INT64):: nsteps_until_prod
  INTEGER, PARAMETER :: vdw_none = 0
  INTEGER, PARAMETER :: vdw_lj = 1
  INTEGER, PARAMETER :: vdw_cut = 2
  INTEGER, PARAMETER :: vdw_cut_shift = 3
  INTEGER, PARAMETER :: vdw_cut_tail = 4
  INTEGER, PARAMETER :: vdw_minimum = 5
  INTEGER, PARAMETER :: vdw_charmm = 6
  INTEGER, PARAMETER :: vdw_cut_switch = 7
  INTEGER, PARAMETER :: vdw_cut_shift_force = 8
  INTEGER, PARAMETER :: vdw_mie = 9

  INTEGER, PARAMETER :: charge_none = 0
  INTEGER, PARAMETER :: charge_coul = 1
  INTEGER, PARAMETER :: charge_cut = 2
  INTEGER, PARAMETER :: charge_ewald = 3
  INTEGER, PARAMETER :: charge_minimum = 4
  INTEGER, PARAMETER :: charge_dsf = 5
  INTEGER, PARAMETER :: charge_sf = 6

  LOGICAL :: cbmc_charge_sf_flag = .TRUE.

  REAL(DP), DIMENSION(:), ALLOCATABLE :: rcut_cbmc, rcut_cbmcsq
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rcut_vdw, rcut_coul, ron_charmm, roff_charmm, rcut_max
  REAL(DP), DIMENSION(:), ALLOCATABLE :: ron_switch, roff_switch, roff_switch_sq, switch_factor1
  REAL(DP), DIMENSION(:), ALLOCATABLE :: switch_factor2, ron_switch_sq
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rcut_vdwsq, inv_rcut_vdwsq, rcut_coulsq, ron_charmmsq, roff_charmmsq
  REAL(SP), DIMENSION(:), ALLOCATABLE :: rcut_vdwsq_sp, inv_rcut_vdwsq_sp
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rcut9, rcut3
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rcut_vdw3, rcut_vdw6
  REAL(DP) :: edens_cut, rcut_clus, rcut_low, rcut_lowsq
  REAL(SP) :: sp_rcut_lowsq
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

  ! Is insertion/deletion limited to an inner region in the box?
  ! int_none = 0, defined elsewhere
  INTEGER, PARAMETER :: int_sphere = 1
  INTEGER, PARAMETER :: int_cylinder = 2
  INTEGER, PARAMETER :: int_slitpore = 3
  INTEGER, PARAMETER :: int_interface = 4

! Atom placement variables

  ! Univ vectors pointing to 6,18,26, and 50 uniform points.
  REAL(DP), DIMENSION(3,50) :: sphere_vec

  ! Imsl

  INTEGER, PARAMETER :: irule = 3

 !***************************************************************
  !Conversion factors and constants

  !REAL(DP), PARAMETER :: PI=3.1415926536_DP
  !REAL(DP), PARAMETER :: twoPI = 6.2831853072_DP
  !REAL(DP), PARAMETER :: rootPI = 1.7724538509_DP
  !REAL(DP), PARAMETER :: halfPI = 0.5_DP*PI
  !REAL(DP), PARAMETER :: threehalfPI_DP = 3.0_DP*halfPI
  !REAL(SP), PARAMETER :: twoPI_SP = REAL(twoPI,SP)
  !REAL(SP), PARAMETER :: PI_SP = REAL(PI,SP)

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
  ! Moved to Species_Class as its attributes
  !INTEGER :: kappa_ins, kappa_rot, kappa_dih, kappa_dih_pad8, kappa_dih_pad32

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
  INTEGER, PARAMETER :: int_rb_torsion = 6

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
  INTEGER, DIMENSION(:), ALLOCATABLE :: max_molecules, natoms, nmol_start, nring_atoms, nexo_atoms
  INTEGER, DIMENSION(:), ALLOCATABLE :: nbonds, nangles, nangles_fixed
  INTEGER, DIMENSION(:), ALLOCATABLE :: ndihedrals, nimpropers
  INTEGER, DIMENSION(:), ALLOCATABLE :: nfragments, fragment_bonds
  INTEGER, DIMENSION(:), ALLOCATABLE :: natoms_to_read

  ! array to hold the total number of molecules of each species in a given box

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nmols
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nmols_to_make, nmols_to_read
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: atom_ibounds

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
  CHARACTER(23), DIMENSION(:), ALLOCATABLE :: atom_type_list

  INTEGER, DIMENSION(:), ALLOCATABLE :: nbr_vdw_params
  INTEGER, DIMENSION(:), ALLOCATABLE :: n_vdw_p_list

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
  TYPE(Dihedral_Class), DIMENSION(:,:), ALLOCATABLE, TARGET :: dihedral_list, uncombined_dihedral_list

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
  !TYPE(Library_Coords_Class), DIMENSION(:), ALLOCATABLE :: library_coords
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: library_coords
  ! library_coords_dim1 is the size of library_coords in the first dimension
  !      since 1, 2, and 3 on the first axis correspond to x, y, and z, it must be at least 3,
  !      but 4 or 8 might be desirable for the sake of vectorization or alignment.
  INTEGER, PARAMETER :: library_coords_dim1 = 4

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

  LOGICAL :: block_avg

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
  REAL(DP) :: prob_deletion, prob_swap, prob_regrowth, prob_ring, prob_atom_displacement, prob_identity_switch
  REAL(DP), DIMENSION(:), ALLOCATABLE :: prob_rot_species
  REAL(DP), DIMENSION(:), ALLOCATABLE :: prob_swap_species, cum_prob_swap_species
  REAL(DP), DIMENSION(:), ALLOCATABLE :: prob_swap_from_box, cum_prob_swap_from_box

  !Switching groups for identity switch
  INTEGER :: num_groups
  INTEGER, DIMENSION(:, :), ALLOCATABLE :: swap_list
  LOGICAL :: default_switch
  INTEGER :: rotations

  LOGICAL :: l_prob_swap_species, l_prob_swap_from_box

  LOGICAL :: f_dv, f_vratio

  REAL(DP) :: omega_max, disp_max, delta_cos_max, delta_phi_max
  REAL(DP), DIMENSION(:), ALLOCATABLE, TARGET :: prob_species_trans, prob_species_rotate
  REAL(DP), DIMENSION(:), ALLOCATABLE::prob_growth_species ! dimension nspecies
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: max_disp, max_rot

  REAL(DP) :: cut_trans, cut_rot, cut_torsion, cut_volume, cut_angle, cut_insertion, cut_deletion, cut_identity_switch
  REAL(DP) :: cut_swap, cut_regrowth, cut_ring, cut_atom_displacement, cut_lambda

  !*********************************************************************************************************

  ! Timing information
  ! Initial, current and final number of steps
  INTEGER (KIND=INT64) :: i_mcstep, initial_mcstep, n_mcsteps, n_equilsteps, iblock
  ! Information on the output of data
  INTEGER (KIND=INT64) :: nthermo_freq, ncoord_freq, block_avg_freq, nbr_blocks
  REAL(DP) :: data_points_per_block
  INTEGER :: int_coord_style ! 1 = xyz, 2 = custom
  INTEGER,DIMENSION(:),ALLOCATABLE :: nbr_prop_files

  ! Number of properties per file, will have dimension of nbr_prop_files

  INTEGER, DIMENSION(:,:), ALLOCATABLE :: prop_per_file
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: first_open

  LOGICAL :: accept

  LOGICAL :: cbmc_flag, del_flag, phi_Flag, angle_Flag, imp_Flag

  ! Some variables for reaction Monte Carlo
  LOGICAL, DIMENSION(:), ALLOCATABLE :: has_charge
  LOGICAL :: get_fragorder, l_check

  INTEGER :: imreplace, isreplace

  INTEGER :: atompairdim, mol_dim


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
  INTEGER, PARAMETER :: imove_identity_switch = 10
  INTEGER, PARAMETER :: imove_atom_displacement = 11
  INTEGER, PARAMETER :: imove_widom = 12


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

  LOGICAL :: echeck
  INTEGER (KIND=INT64):: echeck_freq

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
REAL(DP), ALLOCATABLE :: pair_vdw_temp(:,:), pair_qq_temp(:,:)

!!!! DSF variables
REAL(DP), ALLOCATABLE, DIMENSION(:) :: alpha_dsf
REAL(DP), ALLOCATABLE, DIMENSION(:) :: dsf_factor1, dsf_factor2

  !!!!!!!!!
  !Widom insertions
  !!!!!!!!!
  LOGICAL :: widom_flag, widom_active
  INTEGER, DIMENSION(:), ALLOCATABLE :: tp_correction
  INTEGER :: widom_locate, widom_species
  INTEGER (KIND=INT64), DIMENSION(:,:), ALLOCATABLE :: overlap_counter

  ! Pregenerated trajectory
  !LOGICAL :: has_H ! Does the simulation use H file(s)?
  LOGICAL :: need_energy

  !!!! Sectors
  ! sector_atoms is indexed by (atom index within sector, sector index)
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE, TARGET :: sector_atoms, sector_atoms_cbmc, sector_atoms_full
  ! sector_index_map is indexed by (x index, y index, z index, box index) to get sector index
  INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: sector_index_map, sector_index_map_cbmc, sector_index_map_full
  ! sector_n_atoms is indexed by (sector index) to get number of atoms in a sector
  INTEGER(4), DIMENSION(:), ALLOCATABLE, TARGET :: sector_n_atoms, sector_n_atoms_cbmc, sector_n_atoms_full
  ! sector_has_atoms is indexed by (x index, y index, z index, box index)
  LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: sector_has_atoms

  LOGICAL :: l_sectors, cbmc_cell_list_flag, full_cell_list_flag
  ! sectorbound and length_cells are indexed by (box dimension, box index)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: sectorbound, sectorbound_cbmc, sectorbound_full
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: length_cells, length_cells_cbmc, length_cells_full
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: cell_length_inv, cell_length_inv_cbmc, cell_length_inv_full


  ! indexed like (n_adj_cell_atoms,dimension in 1:4, x index, y index, z index, box index)
  ! there is no 4th dimension of coordinates; this stores charge instead
  REAL(SP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: adj_cell_rsp, cbmc_cell_rsp
  INTEGER(INT32), DIMENSION(:,:,:,:,:), ALLOCATABLE :: adj_cell_ti, cbmc_cell_ti, cbmc_cell_atomtypes
  INTEGER(4), DIMENSION(:,:,:,:), ALLOCATABLE :: n_adj_cell_atoms, cbmc_cell_n_interact
  INTEGER, DIMENSION(3) :: adj_cellmaxbound

  INTEGER :: max_adj_cell_atoms, cbmc_max_interact
  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: cell_length_recip, real_length_cells


  INTEGER, DIMENSION(3) :: sectormaxbound, sectormaxbound_cbmc, sectormaxbound_full

  INTEGER :: n_occ_sectors, n_occ_sectors_cbmc, n_occ_sectors_full
  INTEGER :: max_sector_natoms, max_sector_natoms_cbmc, max_sector_natoms_full
  INTEGER :: max_occ_sectors, max_occ_sectors_cbmc, max_occ_sectors_full
  
  TYPE(Molecule_Class), TARGET :: widom_molecule
  TYPE(Atom_Class), ALLOCATABLE, DIMENSION(:), TARGET :: widom_atoms
  !$OMP THREADPRIVATE(widom_molecule, widom_atoms, cbmc_flag)



  ! n_widom_subgroups is indexed by (species,box)
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: n_widom_subgroups

  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: widom_cpu_time, widom_wc_time



  REAL(DP) :: Eij_max
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Eij_factor
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: w_max, Eij_w_sum 
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: Eij_freq_total
  INTEGER :: Eij_ind_ubound
  LOGICAL :: est_emax
  !$OMP THREADPRIVATE(Eij_max)




!widom_timing  INTEGER(KIND=INT64) :: n_clo, n_not_clo, n_nrg_overlap
!widom_timing  REAL(DP) :: cell_list_time, normal_overlap_time, non_overlap_time, nrg_overlap_time
!widom_timing  !$OMP THREADPRIVATE(n_clo, n_not_clo, n_nrg_overlap)
!widom_timing  !$OMP THREADPRIVATE(cell_list_time, normal_overlap_time, non_overlap_time, nrg_overlap_time)

  ! widom timing variables and parameters
  LOGICAL, PARAMETER :: widom_timing = .TRUE.
  REAL(DP) :: trial_loop_ins_time, cell_list_ins_time, cell_list_cbmc_nrg_ins_time
  REAL(DP) :: noncell_cbmc_nrg_ins_time, rng_ins_time, cbmc_setup_ins_time
  REAL(DP) :: cbmc_returnzone_ins_time, cbmc_endzone_ins_time
  REAL(DP) :: cbmc_fragment_placement_time, cbmc_dih_time, bitcell_overlap_ins_time
  REAL(DP) :: widom_ewald_recip_time
  INTEGER(INT64) :: cbmc_nonoverlap_ins_count, cbmc_dih_count, bitcell_overlap_ins_checks
  INTEGER(INT64) :: cell_list_ins_checks, cell_list_cbmc_nrg_ins_checks, bitcell_overlap_ins_overlaps
  INTEGER(INT64) :: nrg_ins_overlaps
  !$OMP THREADPRIVATE(trial_loop_ins_time, cell_list_ins_time, cell_list_cbmc_nrg_ins_time)
  !$OMP THREADPRIVATE(noncell_cbmc_nrg_ins_time,rng_ins_time,cbmc_setup_ins_time)
  !$OMP THREADPRIVATE(cbmc_returnzone_ins_time, cbmc_endzone_ins_time)
  !$OMP THREADPRIVATE(bitcell_overlap_ins_time, bitcell_overlap_ins_checks, bitcell_overlap_ins_overlaps)
  !$OMP THREADPRIVATE(cbmc_nonoverlap_ins_count, cbmc_dih_count)
  !$OMP THREADPRIVATE(cell_list_ins_checks, cell_list_cbmc_nrg_ins_checks)
  !$OMP THREADPRIVATE(total_cbmc_time, widom_ewald_recip_time, nrg_ins_overlaps)
  REAL(DP) :: trial_loop_ins_time_redux, cell_list_ins_time_redux, cell_list_cbmc_nrg_ins_time_redux
  REAL(DP) :: noncell_cbmc_nrg_ins_time_redux, rng_ins_time_redux, cbmc_setup_ins_time_redux
  REAL(DP) :: cbmc_returnzone_ins_time_redux, cbmc_endzone_ins_time_redux
  REAL(DP) :: cbmc_fragment_placement_time_redux, cbmc_dih_time_redux, bitcell_overlap_ins_time_redux
  INTEGER(INT64) :: cbmc_nonoverlap_ins_count_redux, cbmc_dih_count_redux, bitcell_overlap_ins_checks_redux
  INTEGER(INT64) :: cell_list_ins_checks_redux, cell_list_cbmc_nrg_ins_checks_redux, bitcell_overlap_ins_overlaps_redux
  REAL(DP) :: noncbmc_time_total, total_cbmc_time_redux, widom_ewald_recip_time_redux, nrg_ins_overlaps_redux

  !!! atompair energy table global variables
  INTEGER :: atompair_nrg_res
  REAL(SP) :: atompair_nrg_res_sp
  LOGICAL :: precalc_atompair_nrg
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: atompair_nrg_table
  REAL(SP), DIMENSION(:,:,:), ALLOCATABLE :: atompair_nrg_table_reduced
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: typepair_nrg_table
  REAL(DP) :: rsq_step, inv_rsq_step
  REAL(SP) :: inv_rsq_step_sp
  REAL(DP) :: rsq_shifter
  INTEGER, DIMENSION(:), ALLOCATABLE :: typepair_solute_indices, typepair_solvent_indices
  INTEGER, DIMENSION(:), ALLOCATABLE :: solute_atomtypes, solvent_atomtypes
  INTEGER :: solute_ntypes, solvent_ntypes, solute_maxind, solvent_maxind
  LOGICAL :: need_solvents
  !!!!

  !!! atompair rminsq table global variables
  ! swi stands for single Widom insertion
  ! index swi_atompair_rsqmin with (solvent_base+solvent_ia,solute_ia)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, TARGET :: swi_atompair_rsqmin
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: rsqmin_atompair_w_max
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: rsqmin_atompair_w_sum
  INTEGER(KIND=INT64), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: rsqmin_atompair_freq
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE, TARGET :: atompair_rminsq_table
  REAL(SP), DIMENSION(:,:,:), ALLOCATABLE :: sp_atompair_rminsq_table
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: solvent_max_rminsq
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: solvent_min_rminsq
  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: solvent_max_rminsq_sp
  INTEGER, DIMENSION(:), ALLOCATABLE :: typepair_wsolute_indices, wsolute_atomtypes
  REAL(DP) :: maxrminsq, rsqmin_step, rsqmin_shifter
  INTEGER :: rsqmin_res, wsolute_ntypes, wsolute_maxind
  LOGICAL :: est_atompair_rminsq, read_atompair_rminsq, l_heap
  REAL(DP), DIMENSION(:), ALLOCATABLE :: tol_list
  INTEGER :: nbr_tols, solvent_maxind_d, rsqmin_res_d
  !$OMP THREADPRIVATE(swi_atompair_rsqmin)
  !

  !
  REAL(DP), DIMENSION(0:1000)  :: type_charge_min, type_charge_max
  REAL(DP), DIMENSION(:,:), ALLOCATABLE, TARGET :: rminsq_table
  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: sp_rminsq_table
  REAL(DP) :: U_max_base, max_rmin
  LOGICAL :: calc_rmin_flag
  REAL(DP), DIMENSION(:), ALLOCATABLE :: atomtype_max_rminsq
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: atomtype_min_rminsq
  REAL(SP), DIMENSION(:), ALLOCATABLE :: atomtype_max_rminsq_sp


  LOGICAL, PARAMETER :: l_vectorized = .TRUE.
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nlive

  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: ppvdwp_table, ppvdwp_table2
  REAL(SP), DIMENSION(:,:,:,:), ALLOCATABLE :: ppvdwp_table_sp, ppvdwp_table2_sp
  TYPE(VdW256), DIMENSION(:,:,:), ALLOCATABLE :: ppvdwp_list
  LOGICAL :: l_nonuniform_exponents

  LOGICAL, PARAMETER :: l_not_all_live = .FALSE.

  INTEGER :: global_nthreads

  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: chunksize_array, nthreads_used_array
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: chunks_set_nmols

  LOGICAL :: l_debug_print

  REAL(DP), DIMENSION(:), ALLOCATABLE :: zero_field, rijsq_field
  INTEGER, DIMENSION(:), ALLOCATABLE :: vec123
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: live_xcom, live_ycom, live_zcom, live_max_dcom
  !TYPE(Atom256), DIMENSION(:,:,:,:), ALLOCATABLE :: live_atom_list
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE :: live_atom_rsp
  LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: live_atom_exist
  INTEGER :: maxnmols, maxboxnatoms
  LOGICAL :: l_not_all_exist



  LOGICAL :: bitcell_flag
  REAL(DP) :: min_ideal_bitcell_length
  INTEGER :: solvents_or_types_maxind


  ! Moved to Species_Class as its attributes
  !INTEGER :: kappa_ins_pad8, kappa_ins_pad64

  LOGICAL :: l_zerotype_present




  ! Moved to Species_Class as its attributes
  !REAL(DP), DIMENSION(:,:), ALLOCATABLE :: sincos_lintheta_dp
  !REAL(SP), DIMENSION(:,:), ALLOCATABLE :: sincos_lintheta_sp


  ! Use these parameters below to efficiently round positive integers up (not down) to the nearest multiple of 8, 16, etc.
  ! If the original integer is already a multiple of 8, 16, etc., the answer is the same as the original number.
  ! It actually works as long as the correct answer isn't negative, even if the original integer is 0 or slightly negative
  ! Example: n_pad8 = IAND(n+7,pad8mask), n_pad64 = IAND(n+63,pad64mask)
  ! In those examples, n_pad8 is a multiple of 8, n_pad64 is a multiple of 64, and n is the original number to be "padded"
  ! Note that 7 is equivalent to MASKR(3) and 63 is equivalent to MASKR(6)
  ! This technique only works for padding to positive multiples of power of 2
  INTEGER(INT32), PARAMETER :: pad8mask = NOT(MASKR(3,INT32))
  INTEGER(INT32), PARAMETER :: pad16mask = NOT(MASKR(4,INT32))
  INTEGER(INT32), PARAMETER :: pad32mask = NOT(MASKR(5,INT32))
  INTEGER(INT32), PARAMETER :: pad64mask = NOT(MASKR(6,INT32))


  INTEGER(INT64), PARAMETER :: recip_sqrt_magic_number = INT(Z'5FE6EB50C7B537A9',INT64)

  INTEGER :: nspecies_present
  INTEGER, DIMENSION(:), ALLOCATABLE :: which_species_present

  INTEGER, DIMENSION(3) :: dummy3vec


  INTEGER, DIMENSION(:), ALLOCATABLE :: which_solvent_atomtypes, which_solvent_atomtypes_inv
  INTEGER, DIMENSION(:), ALLOCATABLE :: which_wsolute_atomtypes, which_wsolute_atomtypes_inv
  INTEGER :: n_solvent_atomtypes, n_wsolute_atomtypes


  INTEGER :: n_big_atoms

  TYPE(Cavity_Data_Class), DIMENSION(:,:), ALLOCATABLE :: cavdatalist

  LOGICAL :: cavity_biasing_flag

END MODULE Global_Variables

