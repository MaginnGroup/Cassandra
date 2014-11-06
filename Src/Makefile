#!********************************************************************************
#!   Cassandra - An open source atomistic Monte Carlo software package
#!   developed at the University of Notre Dame.
#!   http://cassandra.nd.edu
#!   Prof. Edward Maginn <ed@nd.edu>
#!   Copyright (2013) University of Notre Dame du Lac
#!
#!   This program is free software: you can redistribute it and/or modify
#!   it under the terms of the GNU General Public License as published by
#!   the Free Software Foundation, either version 3 of the License, or
#!   (at your option) any later version.
#!
#!   This program is distributed in the hope that it will be useful,
#!   but WITHOUT ANY WARRANTY; without even the implied warranty of
#!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#!   GNU General Public License for more details.
#!
#!   You should have received a copy of the GNU General Public License
#!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#!*******************************************************************************
##
###################################################################
##                                                               ##
##  Makefile for building the CASSANDRA simulation package       ##
##  This is the default makefile, and will be involked if you    ##
##  type 'make'. Use this as a template if you want to build     ##
##  your own makefile, or create another makefile and run it     ##
##  using make -f name_of_your_makefile                          ##
###################################################################
##
##  Invocation Options:
##
##   1. make                  Build the Cassandra executable
##   2. make clean            Delete object and executable files
##                            Note you will get a warning upon compilaton 
##			      after invoking make clean, but disregard

###################################################################
##  Master Directory Locations; Change as Needed for Local Site  ##
###################################################################
##
##  CASS_DIR	 CASS Distribution Directory
##  BIN_DIR      Hard Copies of CASS Executables go here

CASS_DIR =  .
BIN_DIR = $(CASS_DIR)
EXEC_NAME = cassandra.exe 

####################################################################
##  Known Machine Types;                                           #
####################################################################

##
##  Machine:  
##  CPU Type: 
##  Compiler: 
##

## Name of compiler and compilation options
FC = ifort
LIBS = 
INCS = 
#
#F90FLAGS = -c  -ffree-line-length-none
F90FLAGS = -c -g -check all -traceback #-warn unused #-openmp
#F90FLAGS = -c -openmp
#F90FLAGS = -p -c -O3
#F90FLAGS = -c -O3 -openmp
#OPTFLAGS = 
#LINKFLAGS = -openmp
#LINKFLAGS = -p

####################################################################
OBJS =	main.o \
        file_names.o \
	type_definitions.o \
	run_variables.o \
	io_utilities.o \
        atoms_to_place.o \
	simulation_properties.o \
        energy_routines.o \
	angle_dist_pick.o \
	input_routines.o \
        write_properties.o \
	read_write_checkpoint.o \
	rotation_routines.o \
	fragment_growth.o \
	volume.o \
        pair_nrg_routines.o \
	nvtmc_control.o \
	nvtmc_driver.o \
	nvt_mc_fragment_control.o \
	nvt_mc_fragment_driver.o \
        nptmc_control.o \
	mcf_control.o \
        nptmc_driver.o \
	create_nonbond_table.o \
        get_internal_coords.o \
	clean_abort.o \
        translate.o \
        random_generators.o \
        save_revert_coordinates.o \
        rotate.o \
        participation.o \
	grow_molecules.o \
        compute_cell_dimensions.o \
        rigid_dihedral_change.o \
	angle_distortion.o \
	volume_change.o \
        gemc_nvt_volume.o \
        create_intra_exclusion_table.o \
        precalculate.o \
        minimum_image_separation.o \
	get_com.o \
        accumulate.o \
        initialize.o \
        gcmc_control.o \
        gcmc_driver.o \
	insertion.o \
	deletion.o \
        accept_or_reject.o \
	zig_by_omega.o \
	gemc_control.o \
        gemc_driver.o \
	gemc_particle_transfer.o \
        cutNgrow.o \
        nvt_mc_ring_fragment.o \
        atom_displacement.o \
        chempot.o \
        pair_nrg_variables.o \


####################################################################

%o: %f90
	${FC} ${F90FLAGS} ${INCS} ${OPTFLAGS} ${LINK_F90} $<

exe: $(OBJS)
	$(FC) $(OPTFLAGS) -o $(BIN_DIR)/$(EXEC_NAME) $(OBJS) $(LIBS) $(LINK_F90) $(LINKFLAGS)

clean:
	rm *.o *.mod dep.mk

dep.mk:
	./get_deps

include dep.mk
