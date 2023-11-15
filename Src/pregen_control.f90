

SUBROUTINE Pregen_Control
        ! Header comments
        !
        !

        USE IO_Utilities
        USE Global_Variables
        USE Type_Definitions
        USE File_Names
        USE Input_Routines
        USE Atoms_To_Place
        USE Angle_Dist_Pick
        USE Energy_Routines
        USE Atompair_Nrg_Table_Routines

        IMPLICIT NONE

        INTEGER :: ierr


        CALL Copy_Inputfile

        ! How many species to simulate?
        CALL Get_Nspecies

        ! Load box shape, number of boxes and box type.  Compute various properties of the box
        ! including the volume
        CALL Get_Box_Info ! Modify for pregen?

        ! Determine the type of VDW and charge interaction model to use, along with
        ! associated parameters and the vdw mixing rule.
        CALL Get_Pair_Style

        ! Determine whether widom insertions are done and get relevant details if they are
        CALL Get_Widom_Info

        CALL Get_Lookup_Info

        ! Load molecular connectivity and force field parameters.  Note that Get_Nspecies
        ! must be called before this routine
        CALL Get_Molecule_Info

        ! Obtain the temperature of the simulation
        CALL Get_Temperature_Info

        CALL Get_Rcutoff_Low

        ! Determine the number and identity of unique atom types, and create a vdw interaction table.
        CALL Create_Nonbond_Table

        ! Create the intramolecular nonbond scaling arrays.
        CALL Create_Intra_Exclusion_Table

        ! No starting configuration but calling Get_Start_Type for allocation
        CALL Get_Start_Type

        ! Seed info
        CALL Get_Seed_Info

        ! All move probabilities will be zero
        CALL Get_Move_Probabilities

        CALL Get_CBMC_Info

        ! Determine the frequency with which information will be output
        CALL Get_Simulation_Length_Info

        ! Write Widom insertion info to log file, which requires units from Get_Simulation_Length_Info
        CALL Log_Widom_Info

        ! Properties to be output
        CALL Get_Property_Info ! may need modification

        CALL Precalculate

        ! Determine the connectivity information for each atom such as how many bonds
        ! it participates in, what atoms are connected to the atom, similar information
        ! for angles and dihedrals

        CALL Participation

        ! Determine the atoms that need to be moved if one of the ends of a bond is perturbed
        CALL Get_Bonds_Atoms_To_Place

        ! Determine what angles a given atom participates in and how many such angles exist
        CALL Get_Angles_Atoms_To_Place

        ! Determine what dihedral angles a given atom participates in and how many such
        ! angles exist
        CALL Get_Dihedral_Atoms_To_Place

        ! Connect to pregenerated trajectory files and get other related settings, if any
        CALL Get_Pregen_Info

        CALL Setup_Atompair_Tables
        

END SUBROUTINE Pregen_Control
