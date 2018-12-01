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
MODULE Potential_Map_Variables
!********************************************************************************
 
    
!*********************************************************************************

USE Type_Definitions

  SAVE

!*********************************************************************************
  ! This section contains global variables used by routines for potential map generation

  REAL(DP) :: grid_spacing, xstep, ystep, zstep
  ! zeolite stuff

  LOGICAL :: l_zeolite_pot
  INTEGER :: n_lat_atoms, nx_zeo, ny_zeo, nz_zeo, n_superlat_atoms
  INTEGER :: na_grid, nb_grid, nc_grid
  INTEGER :: xbox, ybox, zbox
  INTEGER :: n_sorb_atomtypes, total_sorb_grids
  INTEGER :: nlat_types

  INTEGER, ALLOCATABLE :: link_zeo_start(:,:,:), link_zeo_end(:,:,:), zeo_atom_int(:)
  LOGICAL, ALLOCATABLE :: link_zeo(:,:,:)
  INTEGER, ALLOCATABLE :: sorbate_atomtypes(:), n_sorb_grids(:), sorb_grid_pointer(:)
  INTEGER, ALLOCATABLE :: nlat_atom_types(:)

  REAL(DP), ALLOCATABLE, DIMENSION(:) :: x_lat, y_lat, z_lat
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: x_super_lat, y_super_lat, z_super_lat
  INTEGER, ALLOCATABLE, DIMENSION(:) :: type_super_lat

  REAL(DP) :: lxbox, lybox, lzbox

  REAL(DP), ALLOCATABLE, DIMENSION(:) :: a_grid, b_grid, c_grid
  REAL(DP) :: a_spacing, b_spacing, c_spacing
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: xgrid, ygrid, zgrid
  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: zeo_LJ12_table, zeo_LJ6_table
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: grid_pot, grid_potx, grid_poty, grid_potz
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: grid_potxy, grid_potxz, grid_potyz, grid_potxyz
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: grid_potq, grid_potqx, grid_potqy, grid_potqz
  REAL(DP), ALLOCATABLE, DIMENSION(:,:,:) :: grid_potqxy, grid_potqxz, grid_potqyz, grid_potqxyz
 
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: lattice_pot, lattice_potx, lattice_poty, lattice_potz
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: lattice_potxy, lattice_potyz, lattice_potxz, lattice_potxyz

  ! electrostatic stuff

  REAL(DP), ALLOCATABLE, DIMENSION(:) :: charge_super_lat
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: lattice_q_pot, lattice_q_potx, lattice_q_poty, lattice_q_potz
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: lattice_q_potxy, lattice_q_potyz, lattice_q_potxz
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: lattice_q_potxyz, q_grid_pointer

  ! grid stuff
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: grid_q_pointer
  INTEGER :: grid_q_count
  
  ! flag to identify if a species has charge so that the molecules of this species will
  ! use electrostatic map for calculations

  LOGICAL(DP), ALLOCATABLE :: l_elec_map(:)
  REAL(DP), ALLOCATABLE :: zeo_cos_sum(:), zeo_sin_sum(:)
  Type(Box_Class) :: zeo_unit_cell

END MODULE Potential_Map_Variables
