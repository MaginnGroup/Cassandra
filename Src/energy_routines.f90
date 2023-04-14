!*******************************************************************************
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
  !-----------------------------------------------------------------------------
  ! This modules contains a collection of all the routines involved in computing
  ! energies and associated quantities.
  !
  ! Compute_Molecule_Bond_Energy: passed a molecule and species index, this
  !                       returns the total bond energy associated with that
  !                       molecule
  !                       Currently supports none and harmonic.
  !
  ! Compute_Molecule_Angle_Energy: Passed a molecule and species index, it
  !                       returns the total energy of that molecule due to bond
  !                       angles.
  !                       Currently supports none and harmonic.
  !
  ! Compute_Molecule_Dihedral_Energy: Passed a molecule and species index, it
  !                       returns the total dihedral energy of that molecule.
  !                       Currently supports none and OPLS.
  !
  ! Compute_Molecule_Improper_Energy: Passed molecule and species indices, this
  !                       computes the total improper energy of the molecule.
  !                       Not yet tested!!!
  !                       Currently supports harmonic.
  !
  ! Compute_Atom_Nonbond_Energy: passed indices of an atom, molecule and
  !                       species, this returns the vdw and either direct
  !                       charge-charge or the real space part of the Ewald
  !                       energy of this atom with all existing atoms in the
  !                       system. It accounts for intramolecular scaling of 1-2,
  !                       1-3 and 1-4.
  !
  !                       Supports vdw_style = none or LJ
  !                       For LJ, it supports rcut, cut_tail and cut_shift,
  !                       though TAIL CORRECTIONS HAVE NOT YET BEEN ADDED.
  !                       LJ is assumed to be 12-6 LJ.
  !
  !                       Supports charge_style none or coul.
  !                       For charge_style = coul, it supports rcut and Ewald is
  !                       roughed in. However, the Ewald parts of the code need
  !                       some thought, especially in light of computing energy
  !                       differences. This routine also returns the virial
  !                       contribution. It needs more testing, but I believe it
  !                       works.
  !
  ! Compute_Molecule_Nonbond_Intra_Energy: passed molecule and species indices,
  !                       returns the intramolecular LJ and electrostatic energy
  !                       of the molecule.
  !
  ! Compute_Molecule_Nonbond_Inter_Energy: passed molecule and species indices,
  !                       returns the intermolecular LJ and electrostatic energy
  !                       between this molecule and all other molecules in the
  !                       system.
  !
  ! Compute_MoleculePair_Energy:
  ! Compute_MoleculePair_Force:
  !                       Computes the intermolecular energy/force between a
  !                       pair of input molecules.
  !
  ! Compute_AtomPair_Energy:
  ! Compute_AtomPair_Force:
  !                       Computes the vdw and q-q pair energy/force between i
  !                       atoms ia and ja of molecules im and jm of species is
  !                       and js, given their separation rijsq. I have passed
  !                       each component of separation nut right now this is
  !                       unnecessary.
  !                       It also computes the real space part of the Ewald sum
  !                       if necessary.
  !
  !                       LJ potential:
  !                         Eij = 4*epsilon(i,j) *
  !                                 [ (sigma(i,j)/rij)^12 - (sigma(i,j)/rij)^6 ]
  !                         Wij = -rij/3 * d Eij / d rij.
  !                         Use the virial in: P = NkBT + < W >
  !
  ! Compute_AtomPair_Ewald_Real: Real space part of Ewald sum. Need to add
  !                       reciprocal, self and energy difference sin and cos
  !                       sums. Contains erfc function.
  !
  ! Ewald_Reciprocal_Lattice_Vector_Setup : Sets up lattice vectors for Ewald
  !                       Summation for the input box.
  !
  ! Compute_System_Ewald_Reciprocal_Energy:
  ! Compute_System_Ewald_Reciprocal_Force:
  !                       Computes reciprocal space energy/force for a given box
  !
  ! Update_System_Ewald_Reciprocal_Energy:
  !                       Updates the
  !                       reciprocal space energy due to various moves. The
  !                       routine makes use of the fact that for a given move,
  !                       the coordinates of only one molecule are perturbed.
  !                       Hence cos_sum and sin_sum arrays can be computed by
  !                       taking differences of q_i cos(k * r_i) terms in new
  !                       and old configurations.
  !
  ! Compute_System_Ewald_Self_Energy: Calculation of self energy for the Ewald
  !                       summation is obtained from this subroutine.
  !
  ! Compute_Molecule_Ewald_Self_Energy: Computes the self energy of the given
  !                       molecule.
  !
  ! Compute_System_Total_Energy:
  ! Compute_System_Total_Force:
  !                       Computes the total system energy/forces within a given
  !                       box. Forces are then used to compute the pressure
  !                       tensor.
  !
  ! Compute_LR_Correction:
  ! Compute_LR_Force:
  !                       Determines long range correction when the flag is set
  !                       to 'cut_tail'.
  !
  ! Check_MoleculePair_Cutoff:
  !
  ! Check_AtomPair_Cutoff:
  !
  ! Get_Molecule_Energy: Computes the intra- and inter-molecular energy of
  !                       a given molecule interacting with all other molecules.
  !
  ! Compute_Ring_Fragment_Energy: Computes the energy of a ring fragment in its
  !                       old conformation.
  !
  ! Check_System_Energy:
  !
  !
  !
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
  !
  !   gemc_particle_transfer
  !   make_config
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
  !
  ! Revision history
  !
  !   12/10/13  : Beta Release
  !-----------------------------------------------------------------------------

  USE Type_Definitions
  USE Global_Variables
  USE File_Names
  USE Pair_Nrg_Routines
  USE IO_Utilities
 !$  USE OMP_LIB

  IMPLICIT NONE

  REAL(DP), DIMENSION(:), ALLOCATABLE :: zero_field, rijsq_field
  INTEGER, DIMENSION(:), ALLOCATABLE :: vec123
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: live_xcom, live_ycom, live_zcom, live_max_dcom
  !TYPE(Atom256), DIMENSION(:,:,:,:), ALLOCATABLE :: live_atom_list
  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE :: live_atom_rsp
  LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: live_atom_exist
  INTEGER :: maxnmols, maxboxnatoms
  LOGICAL :: l_not_all_exist
  !dir$ attributes align:32 :: zero_field, rijsq_field, vec123, live_xcom, live_ycom, live_zcom, live_max_dcom
  !dir$ assume_aligned zero_field:32, rijsq_field:32, vec123:32, live_xcom:32, live_ycom:32, live_zcom:32, live_max_dcom:32

CONTAINS

  !-----------------------------------------------------------------------------

  SUBROUTINE Compute_Molecule_Bond_Energy(im,is,energy)

    !**************************************************************************
    ! This subroutine computes the total bond energy of a selected molecule
    ! Currently, the available potential functions are none or harmonic.
    ! If none, the code will do a check for fixed bond lengths.
    !
    ! As of now, the code can only support fixed bond length simulations.
    ! Effectively, this subroutine will act as a check for fixed bond length
    ! (useful if using a restart configuration from other packages)
    !
    ! CALLED BY
    !
    !         Compute_System_Total_Energy.
    !         Angle_Distortion
    !         Deletion
    !         Rotate_Dihedral
    !         Insertion
    !         GEMC_Particle_Transfer
    !         Cut_N_Grow
    !
    ! CALLS
    !
    !         Get_Bond_Length
    !
    ! INPUT VARIABLES
    !
    !         im[INTEGER]:     LOCATE of the molecule.
    !         is[INTEGER]:     species type of the molecule.
    !
    ! OUTPUT VARIABLES
    !
    !         energy[REALDP]: total bond energy of molecule
    !
    ! RAISES
    !         It will throw an error if bond lenghts do not match
    !         the MCF specifications within a tolerance.
    !
    !
    ! DOCUMENTATION LAST UPDATED: 08/10/2016
    !
    !**************************************************************************


    INTEGER :: im,is
    REAL(DP) :: energy
    INTEGER :: ib
    REAL(DP) :: length
    REAL(DP) :: k,l0,eb,ltol
    CHARACTER(7) :: mcf_bond_length, current_bond_length
  !-----------------------------------------------------------------------------
    energy = 0.0_DP
    DO ib=1,nbonds(is)
       IF (bond_list(ib,is)%int_bond_type == int_none) THEN
          l0 = bond_list(ib,is)%bond_param(1)
          ltol = bond_list(ib,is)%bond_param(2)
          CALL Get_Bond_Length(ib,im,is,length)
          IF ((int_sim_type /= sim_pregen) .AND. (abs(l0 - length) > ltol)) THEN
             WRITE(mcf_bond_length,'(F7.3)') l0
             WRITE(current_bond_length,'(F7.3)') length
             err_msg = ''
             err_msg(1) = 'Fixed bond is broken between atoms ' &
                        // TRIM(Int_To_String(bond_list(ib,is)%atom1)) // ' and ' &
                        // TRIM(Int_To_String(bond_list(ib,is)%atom2)) &
                        // ' of molecule ' // TRIM(Int_To_String(im)) &
                        // ' of species ' // TRIM(Int_To_String(is))
             err_msg(2) = 'Bond length in MCF:  ' // mcf_bond_length
             err_msg(3) = 'Current bond length: ' // current_bond_length
             CALL Clean_Abort(err_msg, 'Compute_Molecule_Bond_Energy')
          END IF
          eb = 0.0_DP
       ELSEIF (bond_list(ib,is)%int_bond_type == int_harmonic) THEN
          k=bond_list(ib,is)%bond_param(1)
          l0 = bond_list(ib,is)%bond_param(2)
          CALL Get_Bond_Length(ib,im,is,length)
          eb = k*(length-l0)**2

       ENDIF
       energy = energy + eb
    ENDDO


  END SUBROUTINE Compute_Molecule_Bond_Energy
  !-----------------------------------------------------------------------------


  SUBROUTINE Compute_Molecule_Angle_Energy(im,is,energy)
    !**************************************************************************
    ! This subroutine is passed a molecule and species index. It then
    ! computes the total bond angle energy of this molecule.
    !
    ! Currently, the available potential functions are none or harmonic.
    ! If none, the code will do a check for fixed angles.
    !
    ! INPUT VARIABLES
    !
    !         im[INTEGER]:     LOCATE of the molecule.
    !         is[INTEGER]:     species type of the molecule.
    !
    ! OUTPUT VARIABLES
    !
    !         energy[REALDP]:      total bond energy of molecule
    !
    ! RAISES
    !         It will throw an error if angles do not match
    !         the MCF specifications within a tolerance.
    !
    !
    ! DOCUMENTATION LAST UPDATED: 08/10/2016
    !**************************************************************************

    USE Random_Generators
    INTEGER :: im,is
    REAL(DP) :: energy
    INTEGER :: ia
    REAL(DP) :: k,theta0,theta,ea,theta_tol
    CHARACTER (7) :: mcf_angle, current_angle

    energy = 0.0_DP
    DO ia=1,nangles(is)
       IF (angle_list(ia,is)%int_angle_type == int_none) THEN
          theta0 = angle_list(ia,is)%angle_param(1) ! in degrees
          theta_tol = angle_list(ia,is)%angle_param(2) ! in degrees
          CALL Get_Bond_Angle(ia,im,is,theta)
          theta = theta * 180.0_DP / PI
          IF (abs(theta0 - theta) > theta_tol) THEN
             WRITE(mcf_angle,'(F7.3)') theta0
             WRITE(current_angle,'(F7.3)') theta
             err_msg = ''
             err_msg(1) = 'Fixed angle is broken between atoms ' &
                        // TRIM(Int_To_String(angle_list(ia,is)%atom1)) // ' and ' &
                        // TRIM(Int_To_String(angle_list(ia,is)%atom2)) // ' and ' &
                        // TRIM(Int_To_String(angle_list(ia,is)%atom3)) &
                        // ' of molecule ' // TRIM(Int_To_String(im)) &
                        // ' of species ' // TRIM(Int_To_String(is))
             err_msg(2) = 'Angle in MCF:  ' // mcf_angle
             err_msg(3) = 'Current angle: ' // current_angle
             CALL Clean_Abort(err_msg, 'Compute_Molecule_Angle_Energy')
          END IF
          ea = 0.0_DP
       ELSEIF (angle_list(ia,is)%int_angle_type == int_harmonic) THEN
          k=angle_list(ia,is)%angle_param(1)
          theta0 = angle_list(ia,is)%angle_param(2)
          CALL Get_Bond_Angle(ia,im,is,theta)
          ea = k*(theta-theta0)**2
          ! Add more potential functions here.
       ENDIF
       energy = energy + ea
    ENDDO

  END SUBROUTINE Compute_Molecule_Angle_Energy
  !-----------------------------------------------------------------------------

  SUBROUTINE Compute_Molecule_Dihedral_Energy(molecule,species,energy_dihed)
    !**************************************************************************
    ! This routine is passed a molecule and species index. It then computes
    !the total dihedral angle energy of this molecule.
    !
    ! Currently, the available potential functions are OPLS, CHARMM, harmonic,
    ! and none.
    !
    ! INPUT VARIABLES
    !
    !         im[INTEGER]:     LOCATE of the molecule.
    !         is[INTEGER]:     species type of the molecule.
    !
    ! OUTPUT VARIABLES
    !
    !         energy[REALDP]:      total dihedral energy of molecule
    !
    ! RAISES
    !
    ! DOCUMENTATION LAST UPDATED: 08/10/2016
    !**************************************************************************
  USE Global_Variables
    INTEGER :: molecule,species
    REAL(DP) :: energy_dihed
    INTEGER :: idihed, atom1, atom2, atom3, atom4
    REAL(DP) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,edihed,phi,twophi,threephi

    TYPE(Atom_Class), POINTER :: these_atoms(:)

    IF (widom_active) THEN
            these_atoms => widom_atoms
    ELSE
            these_atoms => atom_list(:,molecule,species)
    END IF

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

          IF ( .NOT. these_atoms(atom1)%exist) CYCLE
          IF ( .NOT. these_atoms(atom2)%exist) CYCLE
          IF ( .NOT. these_atoms(atom3)%exist) CYCLE
          IF ( .NOT. these_atoms(atom4)%exist) CYCLE


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

          IF ( .NOT. these_atoms(atom1)%exist) CYCLE
          IF ( .NOT. these_atoms(atom2)%exist) CYCLE
          IF ( .NOT. these_atoms(atom3)%exist) CYCLE
          IF ( .NOT. these_atoms(atom4)%exist) CYCLE

          a0 = dihedral_list(idihed,species)%dihedral_param(1)
          a1 = dihedral_list(idihed,species)%dihedral_param(2)
          a2 = dihedral_list(idihed,species)%dihedral_param(3)

          CALL Get_Dihedral_Angle(idihed,molecule,species,phi)

          edihed = a0 * (1.0_DP + DCOS(a1*phi - a2))

       ELSEIF (dihedral_list(idihed,species)%int_dipot_type == int_harmonic ) THEN

          atom1 = dihedral_list(idihed,species)%atom1
          atom2 = dihedral_list(idihed,species)%atom2
          atom3 = dihedral_list(idihed,species)%atom3
          atom4 = dihedral_list(idihed,species)%atom4

          IF ( .NOT. these_atoms(atom1)%exist) CYCLE
          IF ( .NOT. these_atoms(atom2)%exist) CYCLE
          IF ( .NOT. these_atoms(atom3)%exist) CYCLE
          IF ( .NOT. these_atoms(atom4)%exist) CYCLE

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
  !-----------------------------------------------------------------------------


  SUBROUTINE Compute_Molecule_Improper_Energy(molecule,species,energy)
    !**************************************************************************
    ! This routine is passed the molecule and species index, and returns the
    ! total improper energy of that molecule. Only "none" and "harmonic" types
    ! are supported.
    !
    ! INPUT VARIABLES
    !
    !         im[INTEGER]:     LOCATE of the molecule.
    !         is[INTEGER]:     species type of the molecule.
    !
    ! OUTPUT VARIABLES
    !
    !         energy[REALDP]:      total dihedral energy of molecule
    !
    ! RAISES
    !
    ! DOCUMENTATION LAST UPDATED: 08/10/2016
    !**************************************************************************
    INTEGER :: molecule,species,iimprop
    REAL(DP) :: energy
    REAL(DP) :: eimprop,k,phi0,phi,n_imp,d_imp
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
  !-----------------------------------------------------------------------------

  SUBROUTINE Compute_Atom_Nonbond_Energy(ia,im,is, &
       E_intra_vdw,E_inter_vdw,E_intra_qq,E_inter_qq,overlap)

    !**************************************************************************
    ! Computes the energy components between one particular atom and ALL others
    ! in its box, accounting for exclusions, scalings and existence. It returns
    ! energy components.
    !
    !
    ! Note that the VDW energy (without LRC) is returned as is the real space
    ! part of the q-q interactions (for Ewald and DSF). These two contributions
    ! are categorized into intra or intermolecular energy.
    !
    ! INPUT VARIABLES
    !
    !         ia[INTEGER]:         atom number
    !         im[INTEGER]:     LOCATE of the molecule.
    !         is[INTEGER]:      species type of the molecule.
    !
    ! OUTPUT VARIABLES
    !
    !         E_intra_vdw[REALDP]:        Intramolecular vdw energy of atom
    !         E_inter_vdw[REALDP]:        Intermolecular vdw energy of atom
    !         E_intra_qq[REALDP]:         Intramolecular qq energy of atom
    !         E_inter_qq[REALDP]:         Intermolecular qq energy of atom
    !         Overlap[LOGICAL]:           Flag that gets triggered if atom
    !                                     has a core overlap with another
    !
    ! RAISES
    !
    ! DOCUMENTATION LAST UPDATED: 08/10/2016
    !**************************************************************************


    INTEGER, INTENT(IN) :: ia,im,is
    REAL(DP), INTENT(OUT) :: E_intra_vdw,E_inter_vdw,E_intra_qq,E_inter_qq
    LOGICAL, INTENT(OUT) :: overlap
    INTEGER :: this_box,js,jm,ja, mol_js, itype, jtype, rinteraction, vdw_in
    REAL(DP) :: rxij,ryij,rzij,rijsq,rxijp,ryijp,rzijp
    REAL(DP) :: Eij_intra_vdw,Eij_inter_vdw,Eij_intra_qq,Eij_inter_qq
    REAL(DP) :: eps, sig, SigOverRsq, SigOverR6, SigOverR12
    REAL(DP) :: rij, erf_val, erfc_val, qsc
    REAL(DP) :: T, x, xsq, TP
    REAL(DP) :: rcom,rx,ry,rz
    REAL(DP) :: rcut, rcutsq
    REAL(DP) :: SigOverR, SigOverRn, SigOverRm, mie_coeff,  mie_n, mie_m

    LOGICAL :: get_vdw,get_qq, get_interaction

    REAL(DP), PARAMETER :: A1 = 0.254829592_DP, A2 = -0.284496736_DP
    REAL(DP), PARAMETER :: A3 = 1.421413741_DP, A4 = -1.453152027_DP
    REAL(DP), PARAMETER :: A5 = 1.061405429_DP, P = 0.3275911_DP

    TYPE(Molecule_Class), POINTER :: this_molecule_i, this_molecule_j
    TYPE(Atom_Class), POINTER :: these_atoms_i(:), these_atoms_j(:)

    IF (widom_active) THEN 
            this_molecule_i => widom_molecule
            these_atoms_i => widom_atoms
    ELSE
            this_molecule_i => molecule_list(im,is)
            these_atoms_i => atom_list(:,im,is)
    END IF

    !---------------------------------------------------------------------------
    E_inter_vdw = 0.0_DP
    E_intra_vdw = 0.0_DP
    E_inter_qq = 0.0_DP
    E_intra_qq = 0.0_DP
    Eij_inter_vdw = 0.0_DP
    Eij_intra_vdw = 0.0_DP
    Eij_inter_qq = 0.0_DP
    Eij_intra_qq = 0.0_DP

    ! Check that ia exists
    IF (.NOT. these_atoms_i(ia)%exist ) THEN
       err_msg = ""
       err_msg(1) = 'Attempt to compute energy of an atom that does not exist'
       CALL Clean_Abort(err_msg,'Compute_Atom_Nonbond_Energy')
    ENDIF

    ! Set the box number this particular atom is in.
    this_box = this_molecule_i%which_box

    ! Initialize flags which force a call to pair_energy
    get_vdw = .FALSE.
    get_qq = .FALSE.

    ! Initialize the overlap flag to false to indicate no overlap between atoms.
    overlap = .FALSE.

    SpeciesLoop:DO js=1,nspecies

       MoleculeLoop:DO mol_js=1,nmols(js,this_box)

          jm = locate(mol_js,js,this_box) ! molecule INDEX
          IF (jm == widom_locate .AND. js == widom_species) THEN
                  this_molecule_j => widom_molecule
                  these_atoms_j => widom_atoms
          ELSE
                  this_molecule_j => molecule_list(jm,js)
                  these_atoms_j => atom_list(:,jm,js)
          END IF
          IF (.NOT. this_molecule_j%live) CYCLE MoleculeLoop

          ! Check tos see if atom is to interact with the molecule based
          ! on COM cutoff.
          CALL Check_MoleculePair_Cutoff(jm,js,im,is, &
               get_interaction,rcom,rx,ry,rz)

          IF (.NOT. get_interaction) CYCLE MoleculeLoop

          AtomLoop:DO ja=1,natoms(js)
             ! Test for intramolecular interaction
             IF (.NOT. these_atoms_j(ja)%exist ) CYCLE AtomLoop
             IF (js == is .AND. jm == im) THEN

                IF (ja == ia) THEN
                   ! Avoid computing energy with self
                   CYCLE AtomLoop
                ELSE
                   ! Intra energy. Do not apply PBC
                   IF ( .NOT. these_atoms_j(ja)%exist) CYCLE AtomLoop

                   ! Interatomic distance
                   rxij = these_atoms_j(ja)%rxp &
                        - these_atoms_i(ia)%rxp
                   ryij = these_atoms_j(ja)%ryp &
                        - these_atoms_i(ia)%ryp
                   rzij = these_atoms_j(ja)%rzp &
                        - these_atoms_i(ia)%rzp

                   rijsq = rxij*rxij + ryij*ryij + rzij*rzij

                   IF (rijsq <= rcut_lowsq) THEN
                      IF (.not.(l_bonded(ja,ia,js))) THEN
                         overlap = .true.
                         RETURN
                      ENDIF
                   END IF
                ENDIF

             ELSE
                ! Intermolecular energy so apply pbc.

                ! First compute the parent separation
                rxijp = these_atoms_j(ja)%rxp &
                      - these_atoms_i(ia)%rxp
                ryijp = these_atoms_j(ja)%ryp &
                      - these_atoms_i(ia)%ryp
                rzijp = these_atoms_j(ja)%rzp &
                      - these_atoms_i(ia)%rzp

                ! Now get the minimum image separation
                CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp, &
                     rxij,ryij,rzij)

                rijsq = rxij*rxij + ryij*ryij + rzij*rzij

                IF (rijsq < rcut_lowsq) THEN
                   overlap = .true.
                   RETURN
                END IF

             ENDIF

             CALL Check_AtomPair_Cutoff(rijsq,get_vdw,get_qq,this_box)

             ! Compute vdw and q-q energy using if required
             IF (get_vdw .OR. get_qq) THEN

                CALL Compute_AtomPair_Energy(rxij,ryij,rzij,rijsq, &
                     js,jm,ja,is,im,ia,&
                     get_vdw,get_qq, &
                     Eij_intra_vdw,Eij_intra_qq,Eij_inter_vdw,Eij_inter_qq)

                E_intra_vdw = E_intra_vdw + Eij_intra_vdw
                E_intra_qq  = E_intra_qq  + Eij_intra_qq
                E_inter_vdw = E_inter_vdw + Eij_inter_vdw
                E_inter_qq  = E_inter_qq  + Eij_inter_qq
             ENDIF

          END DO AtomLoop

       END DO MoleculeLoop

    END DO SpeciesLoop

  END SUBROUTINE Compute_Atom_Nonbond_Energy

  !-----------------------------------------------------------------------------

  SUBROUTINE Compute_Atom_Nonbond_Intra_Energy(ia,im,is, &
                  E_intra_vdw,E_intra_qq,E_inter_qq,overlap)
          INTEGER, INTENT(IN) :: ia,im,is
          REAL(DP), INTENT(OUT) :: E_intra_vdw, E_intra_qq, E_inter_qq
          LOGICAL, INTENT(OUT) :: overlap
          TYPE(Atom_Class), DIMENSION(:), POINTER :: these_atoms
          INTEGER :: ja, this_box
          LOGICAL :: get_vdw, get_qq, nonzero_vdw, nonzero_qq
          REAL(DP) :: Eij_intra_vdw, Eij_intra_qq, Eij_inter_vdw, Eij_inter_qq
          REAL(DP) :: rxij, ryij, rzij, rijsq
          E_intra_vdw = 0.0_DP
          E_intra_qq = 0.0_DP
          E_inter_qq = 0.0_DP
          overlap = .TRUE.
          IF (widom_active) THEN
                  these_atoms => widom_atoms
          ELSE
                  these_atoms => atom_list(:,im,is)
          END IF
          this_box = molecule_list(im,is)%which_box
          DO ja = 1, natoms(is)
                IF (ja == ia .OR. .NOT. these_atoms(ja)%exist) CYCLE
                nonzero_vdw = vdw_intra_scale(ia,ja,is) > 0.0_DP 
                nonzero_qq = charge_intra_scale(ia,ja,is) > 0.0_DP 
                IF (.NOT. (nonzero_vdw .OR. nonzero_qq)) CYCLE
                rxij = these_atoms(ja)%rxp - these_atoms(ia)%rxp
                ryij = these_atoms(ja)%ryp - these_atoms(ia)%ryp
                rzij = these_atoms(ja)%rzp - these_atoms(ia)%rzp
                !CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp,rxij,ryij,rzij)
                rijsq = rxij*rxij+ryij*ryij+rzij*rzij
                IF (rijsq < rcut_lowsq .AND. nonzero_vdw) RETURN
                CALL Check_AtomPair_Cutoff(rijsq,get_vdw,get_qq,this_box)
                get_vdw = get_vdw .AND. nonzero_vdw
                get_qq = get_qq .AND. nonzero_qq
                ! Compute vdw and q-q energy if required
                IF (get_vdw .OR. get_qq) THEN
                   CALL Compute_AtomPair_Energy(rxij,ryij,rzij,rijsq, &
                        is,im,ja,is,im,ia,&
                        get_vdw,get_qq, &
                        Eij_intra_vdw,Eij_intra_qq,Eij_inter_vdw,Eij_inter_qq)
                   E_intra_vdw = E_intra_vdw + Eij_intra_vdw
                   E_intra_qq  = E_intra_qq  + Eij_intra_qq
                   E_inter_qq = E_inter_qq + Eij_inter_qq
                ENDIF
          END DO
          overlap = .FALSE.
  END SUBROUTINE Compute_Atom_Nonbond_Intra_Energy

  !-----------------------------------------------------------------------------

!  SUBROUTINE Compute_Atom_Nonbond_Inter_Energy_Cells(ia,im,is, &
!       E_inter_vdw,E_inter_qq)
!        INTEGER, INTENT(IN) :: ia,im,is
!        REAL(DP), INTENT(OUT):: E_inter_vdw, E_inter_qq
!        REAL(DP) :: Eij_intra_vdw, Eij_intra_qq, Eij_inter_vdw, Eij_inter_qq
!        INTEGER :: grid_length(3), this_box, i
!        INTEGER, DIMENSION(:), POINTER :: xi, yi, zi, thisrange_cells, sector_atom_ID, these_cells
!        INTEGER :: dummy_ind, dummy, n_cells_occupied, icell, ia_cell, cell_coords(3), secind
!        !LOGICAL, DIMENSION(:,:,:), POINTER :: filtered_mask, this_mask
!        LOGICAL :: get_vdw, get_qq
!        REAL(DP) :: rijsq, rxijp, ryijp, rzijp, rxij, ryij, rzij, cp(3)
!        TYPE(Atom_Class), POINTER :: atom_ptr
!        INTEGER, DIMENSION(:), POINTER :: this_yb
!        INTEGER, DIMENSION(:,:), POINTER :: this_zb
!        INTEGER :: ix, iy, iz, iy_c, iz_c
!        E_inter_vdw = 0.0_DP
!        E_inter_qq = 0.0_DP
!        IF (widom_active) THEN
!                cp(1) = widom_atoms(ia)%rxp
!                cp(2) = widom_atoms(ia)%ryp
!                cp(3) = widom_atoms(ia)%rzp
!                this_box = widom_molecule%which_box
!        ELSE
!                cp(1) = atom_list(ia,im,is)%rxp
!                cp(2) = atom_list(ia,im,is)%ryp
!                cp(3) = atom_list(ia,im,is)%rzp
!                this_box = molecule_list(im,is)%which_box
!        END IF
!        IF (cbmc_flag) THEN
!                thisrange_cells => cbmcrange_cells(:,this_box)
!!                this_mask => cbmc_mask(-thisrange_cells(1):thisrange_cells(1), &
!!                        -thisrange_cells(2):thisrange_cells(2), &
!!                        -thisrange_cells(3):thisrange_cells(3), &
!!                        this_box)
!                this_yb => cbmc_yb(this_box, -thisrange_cells(1):thisrange_cells(1))
!                this_zb => cbmc_zb(this_box, -thisrange_cells(1):thisrange_cells(1), &
!                        -thisrange_cells(2):thisrange_cells(2))
!        ELSE
!                thisrange_cells => cutrange_cells(:,this_box)
!!                this_mask => cut_mask(-thisrange_cells(1):thisrange_cells(1), &
!!                        -thisrange_cells(2):thisrange_cells(2), &
!!                        -thisrange_cells(3):thisrange_cells(3), &
!!                        this_box)
!                this_yb => cut_yb(this_box, -thisrange_cells(1):thisrange_cells(1))
!                this_zb => cut_zb(this_box, -thisrange_cells(1):thisrange_cells(1), &
!                        -thisrange_cells(2):thisrange_cells(2))
!        END IF
!        grid_length = thisrange_cells*2+1
!        !filtered_mask => filtered_mask_super(1:grid_length(1),1:grid_length(2),1:grid_length(3))
!        cell_coords = IDNINT(cp*cell_length_inv(:,this_box))
!        xi => ci_grid(1,1:grid_length(1))
!        yi => ci_grid(2,1:grid_length(2))
!        zi => ci_grid(3,1:grid_length(3))
!        dummy = cell_coords(1) - thisrange_cells(1)
!        DO i = 1, grid_length(1)
!                xi(i) = dummy
!                dummy = dummy + 1
!        END DO
!        dummy = cell_coords(2) - thisrange_cells(2)
!        DO i = 1, grid_length(2)
!                yi(i) = dummy
!                dummy = dummy + 1
!        END DO
!        dummy = cell_coords(3) - thisrange_cells(3)
!        DO i = 1, grid_length(3)
!                zi(i) = dummy
!                dummy = dummy + 1
!        END DO
!        IF (cell_coords(1)+thisrange_cells(1)>sectorbound(1,this_box)) THEN
!                dummy_ind = grid_length(1) + 1 - (cell_coords(1)+thisrange_cells(1)-sectorbound(1,this_box))
!                xi(dummy_ind:grid_length(1)) = xi(dummy_ind:grid_length(1)) - length_cells(1,this_box)
!        ELSE IF (cell_coords(1)-thisrange_cells(1)<-sectorbound(1,this_box)) THEN
!                dummy_ind = thisrange_cells(1)-sectorbound(1,this_box) - cell_coords(1)
!                xi(1:dummy_ind) = xi(1:dummy_ind) + length_cells(1,this_box)
!        END IF
!        IF (cell_coords(2)+thisrange_cells(2)>sectorbound(2,this_box)) THEN
!                dummy_ind = grid_length(2) + 1 - (cell_coords(2)+thisrange_cells(2)-sectorbound(2,this_box))
!                yi(dummy_ind:grid_length(2)) = yi(dummy_ind:grid_length(2)) - length_cells(2,this_box)
!        ELSE IF (cell_coords(2)-thisrange_cells(2)<-sectorbound(2,this_box)) THEN
!                dummy_ind = thisrange_cells(2)-sectorbound(2,this_box) - cell_coords(2)
!                yi(1:dummy_ind) = yi(1:dummy_ind) + length_cells(2,this_box)
!        END IF
!        IF (cell_coords(3)+thisrange_cells(3)>sectorbound(3,this_box)) THEN
!                dummy_ind = grid_length(3) + 1 - (cell_coords(3)+thisrange_cells(3)-sectorbound(3,this_box))
!                zi(dummy_ind:grid_length(3)) = zi(dummy_ind:grid_length(3)) - length_cells(3,this_box)
!        ELSE IF (cell_coords(3)-thisrange_cells(3)<-sectorbound(3,this_box)) THEN
!                dummy_ind = thisrange_cells(3)-sectorbound(3,this_box) - cell_coords(3)
!                zi(1:dummy_ind) = zi(1:dummy_ind) + length_cells(3,this_box)
!        END IF
!!        filtered_mask = sector_has_atoms(xi,yi,zi,this_box) .AND. this_mask
!!        n_cells_occupied = COUNT(filtered_mask)
!!        these_cells => cell_index_vector(1:n_cells_occupied)
!!        these_cells = PACK(sector_index_map(xi,yi,zi,this_box),filtered_mask)
!!        DO icell = 1, n_cells_occupied
!!                secind = these_cells(icell)
!        iy_c = 1+thisrange_cells(2)
!        iz_c = 1+thisrange_cells(3)
!        DO ix = 1, grid_length(1)
!                DO iy = iy_c-this_yb(ix), iy_c+this_yb(ix)
!                        DO iz = iz_c - this_zb(ix,iy), iz_c+this_zb(ix,iy)
!                                IF (.NOT. sector_has_atoms(xi(ix),yi(iy),zi(iz),this_box)) CYCLE
!                                secind = sector_index_map(xi(ix),yi(iy),zi(iz),this_box)
!                                DO ia_cell = 1, sector_n_atoms(secind)
!                                        sector_atom_ID => sector_atoms(ia_cell,secind,:)
!                                        IF (sector_atom_ID(2) == im .AND. sector_atom_ID(3) == is) CYCLE
!                                        atom_ptr => atom_list(sector_atom_ID(1),sector_atom_ID(2),sector_atom_ID(3))
!                                        rxijp = atom_ptr%rxp - cp(1)
!                                        ryijp = atom_ptr%ryp - cp(2)
!                                        rzijp = atom_ptr%rzp - cp(3)
!                                        CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp,rxij,ryij,rzij)
!                                        rijsq = rxij*rxij+ryij*ryij+rzij*rzij
!                                        CALL Check_AtomPair_Cutoff(rijsq,get_vdw,get_qq,this_box)
!                                        ! Compute vdw and q-q energy using if required
!                                        IF (get_vdw .OR. get_qq) THEN
!                                           CALL Compute_AtomPair_Energy(rxij,ryij,rzij,rijsq, &
!                                                sector_atom_ID(3),sector_atom_ID(2),sector_atom_ID(1),is,im,ia,&
!                                                get_vdw,get_qq, &
!                                                Eij_intra_vdw,Eij_intra_qq,Eij_inter_vdw,Eij_inter_qq)
!                                           E_inter_vdw = E_inter_vdw + Eij_inter_vdw
!                                           E_inter_qq  = E_inter_qq  + Eij_inter_qq
!                                        END IF
!                                END DO
!                        END DO
!                END DO
!        END DO
!!        END DO
!  END SUBROUTINE Compute_Atom_Nonbond_Inter_Energy_Cells

  !-----------------------------------------------------------------------------
  SUBROUTINE Compute_Atom_Nonbond_Inter_Energy_Cells(ia,im,is, &
                  E_inter_vdw, E_inter_qq, overlap, Eij_qq)
          !
          INTEGER, INTENT(IN) :: ia, im, is
          REAL(DP), INTENT(OUT) :: E_inter_vdw, E_inter_qq
          REAL(DP), INTENT(OUT), OPTIONAL :: Eij_qq
          REAL(DP) :: Eij_qq_temp
          INTEGER :: this_box
          INTEGER :: this_locate, i_dim, secind
          INTEGER :: xi, yi, zi, i, ia_cell
          INTEGER, DIMENSION(:), POINTER :: sector_atom_ID
          !INTEGER, DIMENSION(:,:,:), POINTER :: sector_index_map_ptr
          TYPE(Atom_Class), POINTER :: atom_ptr
          INTEGER :: cell_coords(3)
          INTEGER, DIMENSION(3) :: ci_min, ci_max
          REAL(DP) :: cp(3), dx, dy, dz, dxp, dyp, dzp, rijsq
          REAL(DP) :: Eij_intra_vdw,Eij_intra_qq,Eij_inter_vdw,Eij_inter_qq
          REAL(DP), DIMENSION(:,:), POINTER :: cell_length_inv_ptr
          LOGICAL :: get_vdw, get_qq, overlap
          !
          overlap = .TRUE.
          E_inter_vdw = 0.0_DP
          E_inter_qq = 0.0_DP
          Eij_qq_temp = 0.0_DP
          IF (widom_active) THEN
                  cp(1) = widom_atoms(ia)%rxp
                  cp(2) = widom_atoms(ia)%ryp
                  cp(3) = widom_atoms(ia)%rzp
                  this_box = widom_molecule%which_box
          ELSE
                  cp(1) = atom_list(ia,im,is)%rxp
                  cp(2) = atom_list(ia,im,is)%ryp
                  cp(3) = atom_list(ia,im,is)%rzp
                  this_box = molecule_list(im,is)%which_box
          END IF
          IF (cbmc_flag) THEN
                  cell_length_inv_ptr => cell_length_inv_cbmc(:,:,this_box)
          ELSE
                  cell_length_inv_ptr => cell_length_inv_full(:,:,this_box)
          END IF
          cell_coords = IDNINT(MATMUL(cell_length_inv_ptr,cp))
          ci_min = cell_coords - 1
          ci_max = cell_coords + 1
          IF (cbmc_flag) THEN
                  DO xi = ci_min(1), ci_max(1)
                        DO yi = ci_min(2), ci_max(2)
                                DO zi = ci_min(3), ci_max(3)
                                        secind = sector_index_map_cbmc(xi,yi,zi,this_box)
                                        DO ia_cell = 1, sector_n_atoms_cbmc(secind)
                                                sector_atom_ID => sector_atoms_cbmc(ia_cell,secind,:)
                                                atom_ptr => atom_list(sector_atom_ID(1),sector_atom_ID(2),sector_atom_ID(3))
                                                dxp = atom_ptr%rxp - cp(1)
                                                dyp = atom_ptr%ryp - cp(2)
                                                dzp = atom_ptr%rzp - cp(3)
                                                CALL Minimum_Image_Separation(this_box,dxp,dyp,dzp,dx,dy,dz)
                                                rijsq = dx*dx+dy*dy+dz*dz
                                                CALL Check_AtomPair_Cutoff(rijsq,get_vdw,get_qq,this_box)
                                                IF (get_vdw .OR. get_qq) THEN
                                                   IF (rijsq < rcut_lowsq) RETURN
                                                   CALL Compute_AtomPair_Energy(dx,dy,dz,rijsq, &
                                                        sector_atom_ID(3),sector_atom_ID(2),sector_atom_ID(1),is,im,ia,&
                                                        get_vdw,get_qq, &
                                                        Eij_intra_vdw,Eij_intra_qq,Eij_inter_vdw,Eij_inter_qq)
                                                   E_inter_vdw = E_inter_vdw + Eij_inter_vdw
                                                   E_inter_qq  = E_inter_qq  + Eij_inter_qq
                                                END IF
                                        END DO
                                END DO
                        END DO
                  END DO
          ELSE
                  DO xi = ci_min(1), ci_max(1)
                        DO yi = ci_min(2), ci_max(2)
                                DO zi = ci_min(3), ci_max(3)
                                        secind = sector_index_map_full(xi,yi,zi,this_box)
                                        DO ia_cell = 1, sector_n_atoms_full(secind)
                                                sector_atom_ID => sector_atoms_full(ia_cell,secind,:)
                                                atom_ptr => atom_list(sector_atom_ID(1),sector_atom_ID(2),sector_atom_ID(3))
                                                dxp = atom_ptr%rxp - cp(1)
                                                dyp = atom_ptr%ryp - cp(2)
                                                dzp = atom_ptr%rzp - cp(3)
                                                CALL Minimum_Image_Separation(this_box,dxp,dyp,dzp,dx,dy,dz)
                                                rijsq = dx*dx+dy*dy+dz*dz
                                                CALL Check_AtomPair_Cutoff(rijsq,get_vdw,get_qq,this_box)
                                                IF (get_vdw .OR. get_qq) THEN
                                                   IF (rijsq < rcut_lowsq) RETURN
                                                   CALL Compute_AtomPair_Energy(dx,dy,dz,rijsq, &
                                                        sector_atom_ID(3),sector_atom_ID(2),sector_atom_ID(1),is,im,ia,&
                                                        get_vdw,get_qq, &
                                                        Eij_intra_vdw,Eij_intra_qq,Eij_inter_vdw,Eij_inter_qq)
                                                   E_inter_vdw = E_inter_vdw + Eij_inter_vdw
                                                   E_inter_qq  = E_inter_qq  + Eij_inter_qq
                                                   IF (PRESENT(Eij_qq)) THEN
                                                           IF (widom_active) THEN
                                                                   Eij_max = MAX(Eij_max,Eij_inter_vdw+Eij_qq)
                                                           END IF
                                                           Eij_qq_temp = Eij_qq_temp + Eij_qq
                                                   END IF
                                                END IF
                                        END DO
                                END DO
                        END DO
                  END DO
          END IF
          Eij_qq = Eij_qq_temp
          overlap = .FALSE.
          

  END SUBROUTINE Compute_Atom_Nonbond_Inter_Energy_Cells

  SUBROUTINE Compute_Molecule_Nonbond_Intra_Energy(im,is, &
    E_intra_vdw,E_intra_qq,E_inter_qq,intra_overlap)
    !---------------------------------------------------------------------------
    ! The subroutine calculates the intramolecular LJ potential energy and
    ! electrostatic energy of an entire molecule. The routine takes care of
    ! double counting by looping only over i+1 to natoms for ith atom interaction.
    !
    ! Only the minimum image electrostatic energy is stored in E_intra_qq. The
    ! periodic image electrostatic energy is stored in E_inter_qq.
    !
    ! CALLS
    !
    ! Compute_AtomPair_Energy
    !
    ! CALLED BY
    !
    ! Rotate_Dihedral
    ! Angle_Distortion
    !
    !
    ! Written by Jindal Shah on 12/05/07
    !***************************************************************************

    IMPLICIT NONE

    INTEGER :: ia, ja, im, is, this_box

    REAL(DP) :: E_intra_vdw, E_intra_qq, E_inter_qq
    REAL(DP) :: rxij, ryij, rzij, rijsq
    REAL(DP) :: E_intra_vdw_old, E_intra_qq_old
    REAL(DP) :: Eij_intra_vdw, Eij_intra_qq, Eij_inter_vdw, Eij_inter_qq

    LOGICAL :: get_vdw, get_qq, intra_overlap

    TYPE(Atom_Class), POINTER :: these_atoms(:)

    IF (widom_active) THEN
            these_atoms => widom_atoms
    ELSE
            these_atoms => atom_list(:,im,is)
    END IF

    E_intra_vdw = 0.0_DP
    E_intra_qq = 0.0_DP
    E_inter_qq = 0.0_DP
    E_intra_vdw_old = 0.0_DP
    E_intra_qq_old = 0.0_DP

    ! loop over all the atoms in a molecule

    this_box = molecule_list(im,is)%which_box

    DO ia = 1, natoms(is)

       ! check to see if this atom exists
       ! Note 'im' is the linked number of the molecule of interest i.e locate(molecule,is)
       ! The checking for existence of a molecule may be unneccessary.

       IF ( these_atoms(ia)%exist) THEN

          DO ja = ia+1,natoms(is)

             ! make sure that the atom is present

             IF ( .NOT. these_atoms(ja)%exist) CYCLE

             ! Find distance between this atom and all others in the system
             rxij = these_atoms(ia)%rxp - these_atoms(ja)%rxp
             ryij = these_atoms(ia)%ryp - these_atoms(ja)%ryp
             rzij = these_atoms(ia)%rzp - these_atoms(ja)%rzp

             rijsq = rxij*rxij + ryij*ryij + rzij*rzij

             IF (rijsq <= rcut_lowsq) THEN
                IF (.not.(l_bonded(ia,ja,is))) THEN
                   intra_overlap = .true.
                   RETURN
                ENDIF
             END IF

             CALL Check_AtomPair_Cutoff(rijsq,get_vdw,get_qq,this_box)

             IF(cbmc_flag.and.species_list(is)%L_Coul_CBMC) THEN
                get_qq=.false.
             ENDIF

             ! Compute vdw and q-q energy using if required
             IF (get_vdw .OR. get_qq) THEN

                CALL Compute_AtomPair_Energy(rxij,ryij,rzij,rijsq, &
                   is,im,ia,is,im,ja,get_vdw,get_qq, &
                   Eij_intra_vdw,Eij_intra_qq,Eij_inter_vdw,Eij_inter_qq)

                E_intra_vdw = E_intra_vdw + Eij_intra_vdw
                E_intra_qq  = E_intra_qq + Eij_intra_qq
                E_inter_qq  = E_inter_qq + Eij_inter_qq

             END IF

          END DO

       END IF

    END DO

  END SUBROUTINE Compute_Molecule_Nonbond_Intra_Energy
  !-----------------------------------------------------------------------------

  SUBROUTINE Compute_Molecule_Nonbond_Inter_Energy(im,is, &
    E_inter_vdw,E_inter_qq,overlap)
    !***************************************************************************
    ! This subroutine computes interatomic LJ and charge interactions as well as
    ! virials associated with these interactions.
    !
    ! CALLS
    !
    ! Minimum_Image_Separation
    ! Compute_MoleculePair_Energy
    ! Clean_Abort
    !
    ! CALLED BY
    !
    ! Translate
    ! Rotation
    ! Rotate_Dihedral
    ! Angle_Distortion
    ! Insertion
    ! Deletion
    ! Reaction
    !
    ! Written by Jindal Shah on 12/07/07
    !***************************************************************************

    IMPLICIT NONE

!    !$ include 'omp_lib.h'

    INTEGER, INTENT(IN):: im, is
    REAL(DP), INTENT(OUT) :: E_inter_vdw, E_inter_qq
    LOGICAL :: overlap
    !---------------------------------------------------------------------------

    INTEGER  :: ispecies, imolecule, this_box, this_locate

    REAL(DP) :: Eij_vdw, Eij_qq
    REAL(DP) :: eps
    REAL(DP) :: rcom, rx, ry, rz

    LOGICAL :: get_interaction

    INTEGER :: locate_1, locate_2

    LOGICAL :: l_pair_store
    LOGICAL :: my_overlap, shared_overlap

    E_inter_vdw = 0.0_DP
    E_inter_qq = 0.0_DP
    overlap = .FALSE.
    my_overlap = .FALSE.
    shared_overlap = .FALSE.

    this_box = molecule_list(im,is)%which_box

    l_pair_store = .FALSE.

    IF (l_pair_nrg .AND. (.NOT. cbmc_flag)) l_pair_store = .TRUE.

    IF (l_pair_store) CALL Get_Position_Alive(im,is,locate_1)

    speciesLoop: DO ispecies = 1, nspecies

       !$OMP PARALLEL DO DEFAULT(SHARED) &
       !$OMP PRIVATE(imolecule,this_locate,locate_2,get_interaction) &
       !$OMP PRIVATE(rx,ry,rz,rcom,Eij_vdw,Eij_qq) &
       !$OMP SCHEDULE(DYNAMIC) &
       !$OMP REDUCTION(+:E_inter_vdw,E_inter_qq) &
       !$OMP REDUCTION(.OR.:my_overlap)

       moleculeLoop: DO imolecule = 1, nmols(ispecies,this_box)

          IF(shared_overlap) CYCLE

          this_locate = locate(imolecule,ispecies,this_box)
          IF (.NOT. molecule_list(this_locate,ispecies)%live) CYCLE moleculeLoop
          IF (ispecies == is .AND. this_locate == im) CYCLE moleculeLoop

          ! reset pair energy, if storing energies
          IF (l_pair_store) THEN
             CALL Get_Position_Alive(this_locate,ispecies,locate_2)

             pair_nrg_vdw(locate_1,locate_2) = 0.0_DP
             pair_nrg_vdw(locate_2,locate_1) = 0.0_DP

             pair_nrg_qq(locate_1,locate_2) = 0.0_DP
             pair_nrg_qq(locate_2,locate_1) = 0.0_DP
          END IF

          ! Determine if any atoms of these two molecules will interact
          CALL Check_MoleculePair_Cutoff(im,is,this_locate,ispecies,get_interaction, &
               rcom,rx,ry,rz)

          IF (.NOT. get_interaction) CYCLE moleculeLOOP
          IF (cbmc_flag .AND. precalc_atompair_nrg) THEN
                  CALL Estimate_MoleculePair_Energy(im,is,this_locate,ispecies, &
                          this_box,Eij_vdw,my_overlap)
                  ! Eij_vdw given by estimate is actually combined vdW and qq energy
                  ! so don't add Eij_qq
          ELSE
                  CALL Compute_MoleculePair_Energy(im,is,this_locate,ispecies, &
                       this_box,Eij_vdw,Eij_qq,my_overlap)
                  E_inter_qq  = E_inter_qq + Eij_qq
          END IF
          IF (my_overlap) shared_overlap = .TRUE.

          E_inter_vdw = E_inter_vdw + Eij_vdw

       END DO moleculeLoop
       !$OMP END PARALLEL DO

       IF(shared_overlap) THEN
          overlap = .TRUE.
          RETURN
       ENDIF

    END DO speciesLoop

  END SUBROUTINE Compute_Molecule_Nonbond_Inter_Energy
  !-----------------------------------------------------------------------------

  SUBROUTINE Compute_Molecule_Nonbond_Inter_Energy_Widom(im,is, &
    E_inter,overlap)
    !***************************************************************************
    ! This subroutine computes interatomic LJ and charge interactions as well as
    ! virials associated with these interactions.
    !
    ! CALLS
    !
    ! Minimum_Image_Separation
    ! Compute_MoleculePair_Energy
    !
    ! CALLED BY
    !
    !
    !***************************************************************************

    IMPLICIT NONE

!    !$ include 'omp_lib.h'

    INTEGER, INTENT(IN):: im, is
    REAL(DP), INTENT(OUT) :: E_inter
    LOGICAL :: overlap
    !---------------------------------------------------------------------------

    INTEGER  :: ispecies, imolecule, this_box, this_locate, ia

    REAL(DP) :: Eij, Eij_vdw, Eij_qq, Eij_qq_simple, Ei_qq_simple
    REAL(DP) :: eps
    REAL(DP) :: rcom, rx, ry, rz
!molecule_priority    REAL(DP) :: hardcore_max_r, molecule_hardcore_r
    REAL(DP) :: Ei_inter_vdw, Ei_inter_qq

    LOGICAL :: get_interaction
!molecule_priority    LOGICAL, DIMENSION(MAXVAL(nmols(:,widom_molecule%which_box)),nspecies) :: shortrange, midrange




    E_inter = 0.0_DP
    overlap = .FALSE.

    this_box = widom_molecule%which_box

    IF ((cbmc_cell_list_flag .AND. cbmc_flag) .OR. full_cell_list_flag) THEN
            DO ia = 1, natoms(is)
                IF (.NOT. widom_atoms(ia)%exist) CYCLE
                IF (est_emax) THEN
                        CALL Compute_Atom_Nonbond_Inter_Energy_Cells(ia,im,is, &
                                Ei_inter_vdw, Ei_inter_qq, overlap, Ei_qq_simple)
                ELSE
                        CALL Compute_Atom_Nonbond_Inter_Energy_Cells(ia,im,is, &
                                Ei_inter_vdw, Ei_inter_qq, overlap)
                END IF
                IF (overlap) RETURN
                E_inter = E_inter + Ei_inter_vdw + Ei_inter_qq
            END DO
            RETURN
    END IF

    ! Start vectorized version
    IF (l_vectorized .AND. .NOT. (cbmc_flag .AND. precalc_atompair_nrg)) THEN
            CALL Compute_Molecule_Nonbond_Inter_Energy_Vectorized_Widom(im,is,this_box,&
                    Ei_inter_vdw,Ei_inter_qq,overlap)
            IF (.NOT. overlap) E_inter = Ei_inter_vdw + Ei_inter_qq
            RETURN
    END IF




    ! End vectorized version


!molecule_priority    IF (widom_active .AND. l_sectors) THEN
    speciesLoop0: DO ispecies = 1, nspecies
        moleculeLoop0: DO imolecule = 1, nmols(ispecies,this_box)
                this_locate = locate(imolecule,ispecies,this_box)
                IF (ispecies == is .AND. this_locate == im) CYCLE moleculeLoop0
                IF (.NOT. molecule_list(this_locate,ispecies)%live) CYCLE moleculeLoop0
                CALL Check_MoleculePair_Cutoff(im,is,this_locate,ispecies,get_interaction, &
                        rcom,rx,ry,rz)
                IF (.NOT. get_interaction) CYCLE moleculeLoop0
                IF (cbmc_flag .AND. precalc_atompair_nrg) THEN
                        CALL Estimate_MoleculePair_Energy(im,is,this_locate,ispecies, &
                                this_box, Eij, overlap)
                        E_inter = E_inter + Eij 
                ELSE
                        IF (est_emax) THEN
                                CALL Compute_MoleculePair_Energy(im,is,this_locate,ispecies, &
                                     this_box,Eij_vdw,Eij_qq,overlap,Eij_qq_simple)
                        ELSE
                                CALL Compute_MoleculePair_Energy(im,is,this_locate,ispecies, &
                                     this_box,Eij_vdw,Eij_qq,overlap)
                        END IF
                        E_inter = E_inter + Eij_vdw + Eij_qq
                END IF

                IF (overlap) RETURN ! this should never happen when using cell list

        END DO moleculeLoop0
    END DO speciesLoop0
!molecule_priority            RETURN
!molecule_priority    END IF

!molecule_priority    hardcore_max_r = widom_molecule%max_dcom + rcut_low
!molecule_priority    molecule_hardcore_r = rcut_low - widom_molecule%min_dcom
!molecule_priority
!molecule_priority
!molecule_priority    speciesLoop: DO ispecies = 1, nspecies
!molecule_priority       moleculeLoop: DO imolecule = 1, nmols(ispecies,this_box)
!molecule_priority          this_locate = locate(imolecule,ispecies,this_box)
!molecule_priority          IF (ispecies == is .AND. this_locate == im) THEN
!molecule_priority                  shortrange(imolecule, ispecies) = .FALSE.
!molecule_priority                  midrange(imolecule, ispecies) = .FALSE.
!molecule_priority                  CYCLE moleculeLoop
!molecule_priority          ELSE IF (.NOT. molecule_list(this_locate,ispecies)%live) THEN
!molecule_priority                  shortrange(imolecule, ispecies) = .FALSE.
!molecule_priority                  midrange(imolecule, ispecies) = .FALSE.
!molecule_priority                  CYCLE moleculeLoop
!molecule_priority          END IF
!molecule_priority
!molecule_priority          ! Determine whether any atoms of these two molecules will interact
!molecule_priority          CALL Check_MoleculePair_Cutoff(im,is,this_locate,ispecies,get_interaction, &
!molecule_priority               rcom,rx,ry,rz)
!molecule_priority
!molecule_priority          IF (.NOT. get_interaction) THEN
!molecule_priority                  shortrange(imolecule, ispecies) = .FALSE.
!molecule_priority                  midrange(imolecule, ispecies) = .FALSE.
!molecule_priority          ELSE IF (rcom + molecule_list(this_locate,ispecies)%min_dcom < molecule_hardcore_r) THEN
!molecule_priority                  overlap = .TRUE.
!molecule_priority                  RETURN
!molecule_priority          ELSE IF (rcom - molecule_list(this_locate,ispecies)%max_dcom > hardcore_max_r) THEN
!molecule_priority                  shortrange(imolecule, ispecies) = .FALSE.
!molecule_priority                  midrange(imolecule, ispecies) = .TRUE.
!molecule_priority          ELSE
!molecule_priority                  shortrange(imolecule, ispecies) = .TRUE.
!molecule_priority                  midrange(imolecule, ispecies) = .FALSE.
!molecule_priority          END IF
!molecule_priority       END DO moleculeLoop
!molecule_priority    END DO speciesLoop
!molecule_priority
!molecule_priority    speciesLoop2: DO ispecies = 1, nspecies
!molecule_priority       moleculeLoop2: DO imolecule = 1, nmols(ispecies,this_box)
!molecule_priority          IF (.NOT. shortrange(imolecule,ispecies)) CYCLE moleculeLoop2
!molecule_priority          this_locate = locate(imolecule,ispecies,this_box)
!molecule_priority
!molecule_priority          CALL Compute_MoleculePair_Energy(im,is,this_locate,ispecies, &
!molecule_priority               this_box,Eij_vdw,Eij_qq,overlap)
!molecule_priority
!molecule_priority          IF (overlap) RETURN
!molecule_priority
!molecule_priority          E_inter_vdw = E_inter_vdw + Eij_vdw
!molecule_priority          E_inter_qq  = E_inter_qq + Eij_qq
!molecule_priority
!molecule_priority       END DO moleculeLoop2
!molecule_priority    END DO speciesLoop2
!molecule_priority
!molecule_priority    speciesLoop3: DO ispecies = 1, nspecies
!molecule_priority       moleculeLoop3: DO imolecule = 1, nmols(ispecies,this_box)
!molecule_priority          IF (.NOT. midrange(imolecule,ispecies)) CYCLE moleculeLoop3
!molecule_priority          this_locate = locate(imolecule,ispecies,this_box)
!molecule_priority          CALL Compute_MoleculePair_Energy(im,is,this_locate,ispecies, &
!molecule_priority               this_box,Eij_vdw,Eij_qq,overlap)
!molecule_priority
!molecule_priority          IF (overlap) RETURN ! there shouldn't be overlap for midrange molecules
!molecule_priority
!molecule_priority          E_inter_vdw = E_inter_vdw + Eij_vdw
!molecule_priority          E_inter_qq  = E_inter_qq + Eij_qq
!molecule_priority
!molecule_priority       END DO moleculeLoop3
!molecule_priority    END DO speciesLoop3
  END SUBROUTINE Compute_Molecule_Nonbond_Inter_Energy_Widom
  !-----------------------------------------------------------------------------
  SUBROUTINE Field_Allocation
          INTEGER :: ibox, i
          INTEGER :: navec(nbr_boxes), na
          DO ibox = 1, nbr_boxes
                navec(ibox) = DOT_PRODUCT(nmols(:,ibox),natoms)
                maxboxnatoms = MAX(maxboxnatoms,DOT_PRODUCT(nmols(:,ibox),natoms))
          END DO
          na = MAXVAL(navec)
          maxnmols = (MAXVAL(nmols(:,1:))/8+1)*8
          maxboxnatoms = (maxboxnatoms/8+1)*8
          IF (ALLOCATED(zero_field)) THEN
                  IF (SIZE(zero_field) .LE. na) RETURN
                  DEALLOCATE(zero_field)
                  DEALLOCATE(rijsq_field)
                  DEALLOCATE(vec123)
          END IF
          ALLOCATE(zero_field(na))
          ALLOCATE(rijsq_field(na))
          ALLOCATE(vec123(na))
          DO i = 1, na
                vec123(i) = i
          END DO
          zero_field = 0.0_DP
          rijsq_field = 9.0e98_DP
  END SUBROUTINE Field_Allocation


  SUBROUTINE Livelist_Packing(know_all_live,know_unit_stride)
          LOGICAL, INTENT(IN), OPTIONAL :: know_all_live, know_unit_stride
          LOGICAL :: l_all_live, l_unit_stride, l_ortho
          INTEGER :: js, jnlive, jnmols, ibox, istart, iend, jnatoms, maxnlive
          INTEGER :: ja, jm, maxvlen, i, vlen, cjnmols
          REAL(DP) :: lx_recip, ly_recip, lz_recip
          LOGICAL, DIMENSION(maxnmols) :: spec_live
          INTEGER, DIMENSION(maxnmols,nspecies,nbr_boxes) :: live_locates
          INTEGER, DIMENSION(maxnmols) :: spec_locates
          TYPE(Atom_Class), DIMENSION(:,:), ALLOCATABLE :: j_atom_list
          TYPE(Molecule_Class), DIMENSION(:), ALLOCATABLE :: j_molecule_list
          REAL(DP) :: h11,h21,h31,h12,h22,h32,h13,h23,h33
          REAL(DP) :: rxp,ryp,rzp,sxp,syp,szp,jrp
          REAL(DP) :: rcom,xcom,ycom,zcom
          l_not_all_exist = .FALSE.
          IF (PRESENT(know_all_live)) THEN
                  l_all_live = know_all_live
          ELSE
                  l_all_live = .FALSE.
          END IF
          IF (PRESENT(know_unit_stride)) THEN
                  l_unit_stride = know_unit_stride
          ELSE
                  l_unit_stride = .FALSE.
          END IF
          nlive = 0
          IF (l_unit_stride) THEN
                  IF (l_all_live) THEN
                          nlive = nmols(:,1:)
                  ELSE
                          DO js = 1, nspecies
                                istart = 1
                                jnatoms = natoms(js)
                                DO ibox = 1, nbr_boxes
                                        jnmols = nmols(js,ibox)
                                        IF (jnmols == 0) CYCLE
                                        iend = istart + jnmols - 1
                                        nlive(js,ibox) = COUNT(molecule_list(istart:iend,js)%live)
                                        istart = istart + jnmols
                                END DO
                          END DO
                  END IF
          ELSE
                  !$OMP PARALLEL PRIVATE(jnmols,spec_locates,spec_live)
                  !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC)
                  DO ibox = 1, nbr_boxes
                        DO js = 1, nspecies
                                IF (nmols(js,ibox) == 0) CYCLE
                                jnmols = nmols(js,ibox)
                                IF (l_all_live) THEN
                                        live_locates(1:jnmols,js,ibox) = locate(1:jnmols,js,ibox)
                                        nlive(js,ibox) = jnmols
                                        CYCLE
                                END IF
                                spec_locates(1:jnmols) = locate(1:jnmols,js,ibox)
                                spec_live(1:jnmols) = molecule_list(spec_locates(1:jnmols),js)%live
                                IF (ALL(spec_live(1:jnmols))) THEN
                                        nlive(js,ibox) = jnmols
                                        live_locates(1:jnmols,js,ibox) = spec_locates(1:jnmols)
                                ELSE
                                        jnlive = COUNT(spec_live(1:jnmols))
                                        nlive(js,ibox) = jnlive
                                        live_locates(1:jnlive,js,ibox) = PACK(spec_locates(1:jnmols),spec_live(1:jnmols))
                                END IF
                        END DO
                  END DO
                  !$OMP END DO
                  !$OMP END PARALLEL
          END IF
          IF (ALLOCATED(live_xcom)) DEALLOCATE(live_xcom)
          IF (ALLOCATED(live_ycom)) DEALLOCATE(live_ycom)
          IF (ALLOCATED(live_zcom)) DEALLOCATE(live_zcom)
          IF (ALLOCATED(live_max_dcom)) DEALLOCATE(live_max_dcom)
          !IF (ALLOCATED(live_atom_list)) DEALLOCATE(live_atom_list)
          IF (ALLOCATED(live_atom_rsp)) DEALLOCATE(live_atom_rsp)
          IF (ALLOCATED(live_atom_exist)) DEALLOCATE(live_atom_exist)
          maxnlive = MAXVAL(nlive)
          !ALLOCATE(live_atom_list(MAXVAL(natoms),maxnlive,nspecies,nbr_boxes))
          ALLOCATE(live_atom_exist(MAXVAL(natoms),maxnlive,nspecies,nbr_boxes))
          ALLOCATE(live_atom_rsp(MAXVAL(natoms),maxnlive,3,nspecies,nbr_boxes))
          IF (.NOT. (l_all_live .AND. l_unit_stride)) THEN
                  ALLOCATE(j_atom_list(MAXVAL(natoms),maxnlive))
                  ALLOCATE(j_molecule_list(maxnlive))
          END IF
          maxnlive = (maxnlive/8 + 1) * 8
          ALLOCATE(live_xcom(maxnlive,nspecies,nbr_boxes))
          ALLOCATE(live_ycom(maxnlive,nspecies,nbr_boxes))
          ALLOCATE(live_zcom(maxnlive,nspecies,nbr_boxes))
          ALLOCATE(live_max_dcom(maxnlive,nspecies,nbr_boxes))
          DO js = 1, nspecies
                istart = 1
                jnatoms = natoms(js)
                cjnmols = 0
                DO ibox = 1, nbr_boxes
                        jnlive = nlive(js,ibox)
                        jnmols = nmols(js,ibox)
                        vlen = jnlive*jnatoms
                        IF (jnlive == 0) CYCLE
                        iend = istart + jnmols - 1
                        IF (l_all_live .AND. l_unit_stride) THEN
                                !$OMP PARALLEL
                                !$OMP DO SIMD COLLAPSE(2)
                                DO jm = 1, jnlive
                                        DO ja = 1, jnatoms
                                                live_atom_rsp(ja,jm,1,js,ibox) = atom_list(ja,jm+cjnmols,js)%rxp
                                                live_atom_rsp(ja,jm,2,js,ibox) = atom_list(ja,jm+cjnmols,js)%ryp
                                                live_atom_rsp(ja,jm,3,js,ibox) = atom_list(ja,jm+cjnmols,js)%rzp
                                                live_atom_exist(ja,jm,js,ibox) = atom_list(ja,jm+cjnmols,js)%exist
                                                !live_atom_list(ja,jm,js,ibox)%rxp = atom_list(ja,jm+cjnmols,js)%rxp
                                                !live_atom_list(ja,jm,js,ibox)%ryp = atom_list(ja,jm+cjnmols,js)%ryp
                                                !live_atom_list(ja,jm,js,ibox)%rzp = atom_list(ja,jm+cjnmols,js)%rzp
                                                !live_atom_list(ja,jm,js,ibox)%exist = atom_list(ja,jm+cjnmols,js)%exist
                                        END DO
                                END DO
                                !$OMP END DO SIMD
                                !$OMP DO SIMD
                                DO jm = 1, jnlive
                                        live_xcom(jm,js,ibox) = molecule_list(jm+cjnmols,js)%xcom
                                        live_ycom(jm,js,ibox) = molecule_list(jm+cjnmols,js)%ycom
                                        live_zcom(jm,js,ibox) = molecule_list(jm+cjnmols,js)%zcom
                                        live_max_dcom(jm,js,ibox) = molecule_list(jm+cjnmols,js)%max_dcom
                                END DO
                                !$OMP END DO SIMD
                                !$OMP END PARALLEL
                        ELSE
                                IF (l_unit_stride) THEN
                                        !$OMP PARALLEL WORKSHARE
                                        spec_locates(1:jnlive) = PACK(locate(1:jnmols,js,ibox), &
                                                molecule_list(istart:iend,js)%live)
                                        j_molecule_list(1:jnlive) = molecule_list(spec_locates(1:jnlive),js)
                                        j_atom_list(1:jnatoms,1:jnlive) = atom_list(1:jnatoms,spec_locates(1:jnlive),js)
                                        !$OMP END PARALLEL WORKSHARE
                                ELSE
                                        !$OMP PARALLEL WORKSHARE
                                        j_molecule_list(1:jnlive) = molecule_list(live_locates(1:jnlive,js,ibox),js)
                                        j_atom_list(1:jnatoms,1:jnlive) = atom_list(1:jnatoms,live_locates(1:jnlive,js,ibox),js)
                                        !$OMP END PARALLEL WORKSHARE
                                END IF
                                !$OMP PARALLEL
                                !$OMP DO SIMD COLLAPSE(2)
                                DO jm = 1, jnlive
                                        DO ja = 1, jnatoms
                                                live_atom_rsp(ja,jm,1,js,ibox) = j_atom_list(ja,jm)%rxp
                                                live_atom_rsp(ja,jm,2,js,ibox) = j_atom_list(ja,jm)%ryp
                                                live_atom_rsp(ja,jm,3,js,ibox) = j_atom_list(ja,jm)%rzp
                                                live_atom_exist(ja,jm,js,ibox) = j_atom_list(ja,jm)%exist
                                                !live_atom_list(ja,jm,js,ibox)%rxp = j_atom_list(ja,jm)%rxp
                                                !live_atom_list(ja,jm,js,ibox)%ryp = j_atom_list(ja,jm)%ryp
                                                !live_atom_list(ja,jm,js,ibox)%rzp = j_atom_list(ja,jm)%rzp
                                                !live_atom_list(ja,jm,js,ibox)%exist = j_atom_list(ja,jm)%exist
                                        END DO
                                END DO
                                !$OMP END DO SIMD
                                !$OMP DO SIMD
                                DO jm = 1, jnlive
                                        live_xcom(jm,js,ibox) = j_molecule_list(jm)%xcom
                                        live_ycom(jm,js,ibox) = j_molecule_list(jm)%ycom
                                        live_zcom(jm,js,ibox) = j_molecule_list(jm)%zcom
                                        live_max_dcom(jm,js,ibox) = j_molecule_list(jm)%max_dcom
                                END DO
                                !$OMP END DO SIMD
                                !$OMP END PARALLEL
                        END IF
                        !IF (.NOT. ALL(live_atom_list(1:jnatoms,1:jnlive,js,ibox)%exist)) l_not_all_exist = .TRUE.
                        IF (.NOT. ALL(live_atom_exist(1:jnatoms,1:jnlive,js,ibox))) l_not_all_exist = .TRUE.
                        istart = istart + jnmols
                        cjnmols = cjnmols + jnmols
                END DO
          END DO
          DO ibox = 1, nbr_boxes
                l_ortho = box_list(ibox)%int_box_shape <= int_ortho
                IF (l_ortho) CYCLE
                h11 = box_list(ibox)%length_inv(1,1)
                h21 = box_list(ibox)%length_inv(2,1)
                h31 = box_list(ibox)%length_inv(3,1)
                h12 = box_list(ibox)%length_inv(1,2)
                h22 = box_list(ibox)%length_inv(2,2)
                h32 = box_list(ibox)%length_inv(3,2)
                h13 = box_list(ibox)%length_inv(1,3)
                h23 = box_list(ibox)%length_inv(2,3)
                h33 = box_list(ibox)%length_inv(3,3)
                DO js = 1, nspecies
                        jnlive = nlive(js,ibox)
                        jnatoms = natoms(js)
                        !$OMP PARALLEL
                        !$OMP DO SIMD PRIVATE(xcom,ycom,zcom,rcom) &
                        !$OMP ALIGNED(live_xcom,live_ycom,live_zcom)
                        DO jm = 1, jnlive
                                rcom = live_xcom(jm,js,ibox)
                                xcom = h11*rcom
                                ycom = h21*rcom
                                zcom = h31*rcom
                                rcom = live_ycom(jm,js,ibox)
                                xcom = xcom + h12*rcom
                                ycom = ycom + h22*rcom
                                zcom = zcom + h32*rcom
                                rcom = live_zcom(jm,js,ibox)
                                xcom = xcom + h13*rcom
                                ycom = ycom + h23*rcom
                                zcom = zcom + h33*rcom
                                live_xcom(jm,js,ibox) = xcom
                                live_ycom(jm,js,ibox) = ycom
                                live_zcom(jm,js,ibox) = zcom
                        END DO
                        !$OMP END DO SIMD
                        !$OMP DO PRIVATE(jrp,sxp,syp,szp,ja) &
                        !$OMP SCHEDULE(STATIC)
                        DO jm = 1, jnlive
                                DO ja = 1, jnatoms
                                        jrp = live_atom_rsp(ja,jm,1,js,ibox)
                                        sxp = h11*jrp
                                        syp = h21*jrp
                                        szp = h31*jrp
                                        jrp = live_atom_rsp(ja,jm,2,js,ibox)
                                        sxp = sxp + h12*jrp
                                        syp = syp + h22*jrp
                                        szp = szp + h32*jrp
                                        jrp = live_atom_rsp(ja,jm,3,js,ibox)
                                        sxp = sxp + h13*jrp
                                        syp = syp + h23*jrp
                                        szp = szp + h33*jrp
                                        live_atom_rsp(ja,jm,1,js,ibox) = sxp
                                        live_atom_rsp(ja,jm,2,js,ibox) = syp
                                        live_atom_rsp(ja,jm,3,js,ibox) = szp
                                END DO
                        END DO
                        !$OMP END DO
                        !$OMP END PARALLEL
                END DO
          END DO
          maxboxnatoms = MAXVAL(SPREAD(natoms,2,nbr_boxes)*nlive)

  END SUBROUTINE Livelist_Packing

  SUBROUTINE Compute_Molecule_Nonbond_Inter_Energy_Vectorized_Widom(im,is, &
    this_box,vdw_energy,qq_energy,overlap)
    ! Arguments
    INTEGER, INTENT(IN) :: im, is, this_box
    REAL(DP), INTENT(OUT) :: vdw_energy, qq_energy
    LOGICAL, INTENT(OUT) :: overlap
    ! End Arguments
    INTEGER, DIMENSION(nspecies) :: jspec_n_interact, jspec_istart, jspec_iend
    INTEGER, DIMENSION(COUNT(widom_atoms%exist)) :: which_i_exist, ispec_atomtypes
    REAL(DP), DIMENSION(COUNT(widom_atoms%exist)) :: icharge
    REAL(DP), DIMENSION(3,COUNT(widom_atoms%exist)) :: irp, isp
    REAL(DP), DIMENSION(3) :: irp_com, isp_com
    LOGICAL(4), DIMENSION(COUNT(widom_atoms%exist)) :: i_has_vdw
    !LOGICAL, DIMENSION(MAXVAL(nmols(:,this_box))), TARGET :: spec_live_tgt
    REAL(DP), DIMENSION(maxboxnatoms) :: rijsq_4min
    LOGICAL(4), DIMENSION(MAXVAL(nlive(:,this_box))) :: interact_vec
    !INTEGER, DIMENSION(MAXVAL(nlive(:,this_box))) :: which_interact
    LOGICAL(1), DIMENSION(maxboxnatoms) :: vdw_mask, coul_mask, j_hascharge, j_hasvdw, j_exist
    INTEGER, DIMENSION(maxboxnatoms) :: jatomtype, packed_types
    REAL(DP), DIMENSION(maxboxnatoms) :: jcharge, up_nrg_vec, rijsq_packed, rijsq, jcharge_coul, jxp, jyp, jzp
    REAL(DP), DIMENSION(maxboxnatoms,3) :: jrsp
    TYPE(VdW256), DIMENSION(maxboxnatoms) :: ij_vdw_p
    REAL(DP), DIMENSION(4,maxboxnatoms) :: mie_vdw_p_table
    REAL(DP), DIMENSION(2,maxboxnatoms) :: ij_vdw_p_table
    REAL(DP), DIMENSION(maxboxnatoms,4) :: ij_vdw_p_table_T
    LOGICAL :: this_est_emax, l_get_rij_min, l_ortho
    INTEGER :: n_vdw_p, ia_counter, i, j, n_i_exist, istart, iend, jnlive, n_interact, vlen, orig_vlen, n_j_exist
    INTEGER :: bsolvent, istart_base, natoms_js, ti_solvent, ja, js, n_coul, n_vdw, ia, live_vlen, jnmols, jnatoms, itype
    REAL(DP) :: mol_rcut, max_dcom_i_const, this_vdw_rcutsq, this_coul_rcutsq
    REAL(DP) :: sigbyr2, sigbyr6, sigbyr12, rij, roffsq_rijsq, nrg_vdw, nrg_coul, i_qq_energy, mie_m, mie_n, icharge_factor

    INTEGER :: this_int_vdw_style, ibox, jtype

    LOGICAL :: l_check_coul, l_notallvdw, l_notallcoul, l_pack_separately, l_notsamecut, i_get_vdw, i_get_coul
    LOGICAL(8) :: l_interact
    REAL(DP) :: rterm, rterm2, eps, sigsq, negsigsq, negsigbyr2, this_rijsq
    REAL(DP) :: ixp, iyp, izp, dxcom, dycom, dzcom, dscom, ixcom, iycom, izcom, dxp, dyp, dzp, dsp
    REAL(DP) :: xl, yl, zl, hxl, hyl, hzl
    REAL(DP) :: h11,h21,h31,h12,h22,h32,h13,h23,h33
    !!!dir$ attributes align:32 :: interact_vec, vdw_mask, coul_mask, j_hascharge, j_hasvdw, j_exist, jatomtype, packed_types, &
    !dir$ attributes align:32 :: jxp, jyp, jzp, jrsp, ij_vdw_p_table
    !dir$ assume_aligned jxp:32, jyp:32, jzp:32, jrsp:32, ij_vdw_p_table:32, ij_vdw_p_table_T:32
    !DIR$ ASSUME_ALIGNED rijsq:32, rijsq_packed:32, jcharge:32, jcharge_coul:32
    !DIR$ ASSUME_ALIGNED vdw_mask:32, coul_mask:32, j_hascharge:32, j_hasvdw:32, j_exist:32
    !DIR$ ASSUME_ALIGNED jatomtype:32, packed_types:32, up_nrg_vec:32
    !dir$ assume (MOD(maxboxnatoms,8) .EQ. 0)

    vdw_energy = 0.0_DP
    qq_energy = 0.0_DP
    overlap = .FALSE.
    ibox = widom_molecule%which_box


    l_get_rij_min = est_atompair_rminsq .AND. .NOT. cbmc_flag
    this_est_emax = est_emax .AND. widom_active .AND. .NOT. cbmc_flag
    jspec_n_interact = 0
    n_vdw_p = n_vdw_p_list(this_box)
    l_ortho = box_list(this_box)%int_box_shape <= int_ortho
    IF (l_ortho) THEN
            xl = box_list(this_box)%length(1,1)
            yl = box_list(this_box)%length(2,2)
            zl = box_list(this_box)%length(3,3)
            hxl = 0.5 * xl
            hyl = 0.5 * yl
            hzl = 0.5 * zl
    ELSE
            h11 = box_list(this_box)%length(1,1)
            h21 = box_list(this_box)%length(2,1)
            h31 = box_list(this_box)%length(3,1)
            h12 = box_list(this_box)%length(1,2)
            h22 = box_list(this_box)%length(2,2)
            h32 = box_list(this_box)%length(3,2)
            h13 = box_list(this_box)%length(1,3)
            h23 = box_list(this_box)%length(2,3)
            h33 = box_list(this_box)%length(3,3)
    END IF
    IF (cbmc_flag) THEN
            mol_rcut = rcut_cbmc(this_box)
            this_vdw_rcutsq = rcut_cbmcsq(this_box)
            this_coul_rcutsq = rcut_cbmcsq(this_box)
    ELSE
            mol_rcut = rcut_max(this_box)
            IF (int_vdw_sum_style(this_box) == vdw_cut_switch) THEN
                    this_vdw_rcutsq = roff_switch_sq(this_box)
            ELSE
                    this_vdw_rcutsq = rcut_vdwsq(this_box)
            END IF
            this_coul_rcutsq = rcut_coulsq(this_box)
    END IF
    max_dcom_i_const = widom_molecule%max_dcom + mol_rcut
    n_i_exist = COUNT(widom_atoms%exist)
    which_i_exist = PACK(vec123(1:natoms(is)),widom_atoms%exist)
    irp(1,:) = widom_atoms(which_i_exist)%rxp
    irp(2,:) = widom_atoms(which_i_exist)%ryp
    irp(3,:) = widom_atoms(which_i_exist)%rzp
    ispec_atomtypes = nonbond_list(which_i_exist,is)%atom_type_number
    icharge = nonbond_list(which_i_exist,is)%charge * charge_factor
    IF (l_ortho) THEN
            ixcom = widom_molecule%xcom
            iycom = widom_molecule%ycom
            izcom = widom_molecule%zcom
    ELSE
            irp_com(1) = widom_molecule%xcom
            irp_com(2) = widom_molecule%ycom
            irp_com(3) = widom_molecule%zcom
            isp_com = MATMUL(box_list(this_box)%length_inv,irp_com)
            isp = MATMUL(box_list(this_box)%length_inv,irp)
            ixcom = isp_com(1)
            iycom = isp_com(2)
            izcom = isp_com(3)
    END IF
    istart = 1
    DO js = 1, nspecies
        jnlive = nlive(js,this_box)
        IF (jnlive == 0) CYCLE
        jnatoms = natoms(js)
        n_interact = 0
        IF (l_ortho) THEN
                DO i = 1, jnlive
                        dxcom = ABS(live_xcom(i,js,ibox) - ixcom)
                        dycom = ABS(live_ycom(i,js,ibox) - iycom)
                        dzcom = ABS(live_zcom(i,js,ibox) - izcom)
                        IF (dxcom > hxl) dxcom = dxcom - xl
                        IF (dycom > hyl) dycom = dycom - yl
                        IF (dzcom > hzl) dzcom = dzcom - zl
                        ! Repurposing dxcom as rijsq accumulator to enforce FMA3 instructions if supported
                        dxcom = dxcom * dxcom
                        dxcom = dxcom + dycom * dycom
                        dxcom = dxcom + dzcom * dzcom
                        l_interact = max_dcom_i_const > &
                                SQRT(dxcom) - live_max_dcom(i,js,this_box)
                        interact_vec(i) = l_interact
                        IF (l_interact) n_interact = n_interact + 1
                END DO
                !DO i = 1,3
                !        drspw(:,i) = drspw(:,i) - &
                !                ANINT(drspw(:,i)/box_list(this_box)%length(i,i))*box_list(this_box)%length(i,i)
                !END DO
        ELSE
                DO i = 1, jnlive
                        ! COM coordinates are fractional if box is triclinic
                        dscom = live_xcom(i,js,this_box)
                        dscom = dscom - ixcom
                        IF (dscom > 0.5_DP) THEN
                                dscom = dscom - 1.0_DP
                        ELSE IF (dscom < -0.5_DP) THEN
                                dscom = dscom + 1.0_DP
                        END IF
                        dxcom = h11*dscom
                        dycom = h21*dscom
                        dzcom = h31*dscom
                        dscom = live_ycom(i,js,this_box)
                        dscom = dscom - iycom
                        IF (dscom > 0.5_DP) THEN
                                dscom = dscom - 1.0_DP
                        ELSE IF (dscom < -0.5_DP) THEN
                                dscom = dscom + 1.0_DP
                        END IF
                        dxcom = dxcom + h12*dscom
                        dycom = dycom + h22*dscom
                        dzcom = dzcom + h32*dscom
                        dscom = live_zcom(i,js,this_box)
                        dscom = dscom - izcom
                        IF (dscom > 0.5_DP) THEN
                                dscom = dscom - 1.0_DP
                        ELSE IF (dscom < -0.5_DP) THEN
                                dscom = dscom + 1.0_DP
                        END IF
                        dxcom = dxcom + h13*dscom
                        dycom = dycom + h23*dscom
                        dzcom = dzcom + h33*dscom
                        ! Repurposing dxcom as rijsq accumulator to enforce FMA3 instructions if supported
                        dxcom = dxcom * dxcom
                        dxcom = dxcom + dycom * dycom
                        dxcom = dxcom + dzcom * dzcom
                        l_interact = max_dcom_i_const > &
                                SQRT(dxcom) - live_max_dcom(i,js,this_box)
                        interact_vec(i) = l_interact
                        IF (l_interact) n_interact = n_interact + 1
                END DO
                !drspw = MATMUL(drspw,length_inv_T)
                !drspw = drspw - ANINT(drspw)
                !drspw = MATMUL(drspw,length_T)
        END IF
        !which_interact(1:n_interact) = PACK(vec123(1:jnlive),interact_vec(1:jnlive))
        vlen = n_interact*natoms(js)
        iend = istart + vlen - 1
        jrsp(istart:iend,:) = RESHAPE(live_atom_rsp(1:jnatoms,PACK(vec123(1:jnlive),interact_vec(1:jnlive)),1:3,js,this_box), &
                (/ vlen, 3 /))
        IF (l_not_all_exist) THEN
                j_exist(istart:iend) = RESHAPE(&
                        live_atom_exist(1:jnatoms,PACK(vec123(1:jnlive),interact_vec(1:jnlive)),js,this_box), (/ vlen /))
        END IF
        jatomtype(istart:iend) = RESHAPE(SPREAD(nonbond_list(1:natoms(js),js)%atom_type_number,2,n_interact), (/ vlen /))
        jcharge(istart:iend) = RESHAPE(SPREAD(nonbond_list(1:natoms(js),js)%charge,2,n_interact), (/ vlen /))
        IF (l_get_rij_min) THEN
                jspec_n_interact(js) = n_interact
                jspec_istart(js) = istart
                jspec_iend(js) = iend
        END IF
        istart = istart+vlen
    END DO
    vlen = istart - 1
    orig_vlen = vlen
    ! The following IF statement condition is unlikely to be true.
    ! That is why this process is saved until now; it's unlikely to be necessary.
    IF (l_not_all_exist) THEN
            IF (.NOT. ALL(j_exist(1:vlen))) THEN
                    n_j_exist = COUNT(j_exist(1:vlen))
                    DO i = 1,3
                        jrsp(1:n_j_exist,i) = PACK(jrsp(1:vlen,i),j_exist(1:vlen))
                    END DO
                    jcharge(1:n_j_exist) = PACK(jcharge(1:vlen),j_exist(1:vlen))
                    jatomtype(1:n_j_exist) = PACK(jatomtype(1:vlen),j_exist(1:vlen))
                    vlen = n_j_exist
            END IF
    END IF
    IF (this_est_emax) THEN
            up_nrg_vec = 0.0_DP
    END IF
    j_hascharge(1:vlen) = jcharge(1:vlen) .NE. 0.0_DP
    j_hasvdw(1:vlen) = jatomtype(1:vlen) .NE. 0
    l_notallvdw = .NOT. ALL(j_hasvdw(1:vlen))
    l_notallcoul = .NOT. ALL(j_hascharge(1:vlen))
    l_notsamecut = this_vdw_rcutsq .NE. this_coul_rcutsq
    l_pack_separately = l_notsamecut .OR. l_notallcoul .OR. l_notallvdw

    DO ia = 1, n_i_exist
        itype = ispec_atomtypes(ia)
        i_get_vdw = int_vdw_style(this_box) .NE. vdw_none .AND. itype .NE. 0
        i_get_coul = icharge(ia) .NE. 0.0_DP .AND. int_charge_style(this_box) == charge_coul .AND. &
                (species_list(is)%l_coul_cbmc .OR. .NOT. cbmc_flag)
        IF (.NOT. (i_get_vdw .OR. i_get_coul)) CYCLE
        l_check_coul = (i_get_coul .AND. l_notsamecut) .OR. .NOT. i_get_vdw
        !l_pack_coul = (i_get_coul .AND. l_pack_separately) .OR. l_check_coul
        IF (l_ortho) THEN
                ixp = irp(1,ia)
                iyp = irp(2,ia)
                izp = irp(3,ia)
                DO i = 1, vlen
                        dxp = jrsp(i,1)
                        dyp = jrsp(i,2)
                        dzp = jrsp(i,3)
                        dxp = ABS(dxp - ixp)
                        dyp = ABS(dyp - iyp)
                        dzp = ABS(dzp - izp)
                        IF (dxp > hxl) dxp = dxp - xl
                        IF (dyp > hyl) dyp = dyp - yl
                        IF (dzp > hzl) dzp = dzp - zl
                        ! Repurposing dxp as rijsq accumulator to enforce FMA3 instructions if supported
                        dxp = dxp * dxp
                        dxp = dxp + dyp * dyp
                        dxp = dxp + dzp * dzp
                        rijsq(i) = dxp
                        IF (i_get_vdw) THEN
                                l_interact = dxp < this_vdw_rcutsq
                                vdw_mask(i) = l_interact
                        END IF
                        IF (l_check_coul) THEN
                                l_interact = dxp < this_coul_rcutsq
                                coul_mask(i) = l_interact
                        END IF
                END DO
        ELSE
                ixp = isp(1,ia)
                iyp = isp(2,ia)
                izp = isp(3,ia)
                DO i = 1, vlen
                        ! atom coordinates are fractional if box is triclinic
                        dsp = jrsp(i,1)
                        dsp = dsp - ixp
                        IF (dsp > 0.5_DP) THEN
                                dsp = dsp - 1.0_DP
                        ELSE IF (dsp < -0.5_DP) THEN
                                dsp = dsp + 1.0_DP
                        END IF
                        dxp = h11*dsp
                        dyp = h21*dsp
                        dzp = h31*dsp
                        dsp = jrsp(i,2)
                        dsp = dsp - iyp
                        IF (dsp > 0.5_DP) THEN
                                dsp = dsp - 1.0_DP
                        ELSE IF (dsp < -0.5_DP) THEN
                                dsp = dsp + 1.0_DP
                        END IF
                        dxp = dxp + h12*dsp
                        dyp = dyp + h22*dsp
                        dzp = dzp + h32*dsp
                        dsp = jrsp(i,3)
                        dsp = dsp - izp
                        IF (dsp > 0.5_DP) THEN
                                dsp = dsp - 1.0_DP
                        ELSE IF (dsp < -0.5_DP) THEN
                                dsp = dsp + 1.0_DP
                        END IF
                        dxp = dxp + h13*dsp
                        dyp = dyp + h23*dsp
                        dzp = dzp + h33*dsp
                        ! Repurposing dxp as rijsq accumulator to enforce FMA3 instructions if supported
                        dxp = dxp * dxp
                        dxp = dxp + dyp * dyp
                        dxp = dxp + dzp * dzp
                        rijsq(i) = dxp
                        IF (i_get_vdw) THEN
                                l_interact = dxp < this_vdw_rcutsq
                                vdw_mask(i) = l_interact
                        END IF
                        IF (l_check_coul) THEN
                                l_interact = dxp < this_coul_rcutsq
                                coul_mask(i) = l_interact
                        END IF
                END DO
        END IF
        IF (.NOT. l_sectors) THEN
                IF (ANY(rijsq(1:vlen) < rcut_lowsq)) THEN
                        overlap = .TRUE.
                        RETURN
                END IF
        END IF
        IF (l_get_rij_min) THEN
                DO js = 1, nspecies
                        IF (jspec_n_interact(js) == 0) CYCLE
                        bsolvent = species_list(js)%solvent_base
                        iend = jspec_iend(js)
                        istart_base = jspec_istart(js) - 1
                        natoms_js = natoms(js)
                        DO ja = 1, natoms(js)
                                ti_solvent = bsolvent + ja
                                istart = istart_base + ja
                                IF (vlen .NE. orig_vlen) THEN
                                        rijsq_4min(1:orig_vlen) = UNPACK(rijsq(1:vlen),j_exist(1:orig_vlen),rijsq_field(1:orig_vlen))
                                        swi_atompair_rsqmin(ti_solvent,which_i_exist(ia)) = &
                                                MINVAL(rijsq_4min(istart:iend:natoms_js))
                                ELSE
                                        swi_atompair_rsqmin(ti_solvent,which_i_exist(ia)) = &
                                                MINVAL(rijsq(istart:iend:natoms_js))
                                END IF
                        END DO
                END DO
        END IF
        IF (i_get_vdw) THEN
                IF (i_get_coul .AND. .NOT. l_notsamecut) THEN
                        IF (l_notallcoul) THEN
                                coul_mask(1:vlen) = vdw_mask(1:vlen) .AND. j_hascharge(1:vlen)
                        ELSE IF (l_notallvdw) THEN
                                coul_mask(1:vlen) = vdw_mask(1:vlen)
                        ELSE
                                n_coul = COUNT(vdw_mask(1:vlen))
                                jcharge_coul(1:n_coul) = PACK(jcharge(1:vlen),vdw_mask(1:vlen))
                        END IF
                END IF
                IF (l_notallvdw) THEN
                        vdw_mask(1:vlen) = vdw_mask(1:vlen) .AND. j_hasvdw(1:vlen)
                END IF
                n_vdw = COUNT(vdw_mask(1:vlen))
                !packed_types(1:n_vdw) = PACK(jatomtype(1:vlen),vdw_mask(1:vlen))
                !ij_vdw_p(1:n_vdw) = ppvdwp_list(packed_types(1:n_vdw),itype,this_box)
                !ij_vdw_p(1:n_vdw) = ppvdwp_list(PACK(jatomtype(1:vlen),vdw_mask(1:vlen)),itype,this_box)
                !DO i = 1, n_vdw
                !        ij_vdw_p_table(i,1:4) = ppvdwp_table(1:4,packed_types(i),itype,this_box)
                !END DO
                !DO i = 1, n_vdw
                !        ij_vdw_p(i) = ppvdwp_list(packed_types(i),itype,this_box)
                !END DO
                !ij_vdw_p_table(1:n_vdw,1:4) = TRANSPOSE(ppvdwp_table(1:4,packed_types(1:n_vdw),itype,this_box))
                !ij_vdw_p_table(1:4,1:n_vdw) = ppvdwp_table(1:4,packed_types(1:n_vdw),itype,this_box)
                !ij_vdw_p_table(1:4,1:n_vdw) = ppvdwp_table(1:4,PACK(jatomtype(1:vlen),vdw_mask(1:vlen)),itype,this_box)
                rijsq_packed(1:n_vdw) = PACK(rijsq(1:vlen),vdw_mask(1:vlen))
                this_int_vdw_style = int_vdw_style(this_box)
                IF (this_int_vdw_style == vdw_lj) THEN
                        ij_vdw_p_table(1:2,1:n_vdw) = ppvdwp_table(1:2,PACK(jatomtype(1:vlen),vdw_mask(1:vlen)),itype,ibox)
                        IF (int_vdw_sum_style(this_box) == vdw_charmm) THEN
                                DO i = 1, n_vdw
                                        eps = ij_vdw_p_table(1,i) ! epsilon
                                        sigsq = ij_vdw_p_table(2,i) ! sigma**2
                                        sigbyr2 = sigsq/rijsq_packed(i) ! sigma was already squared
                                        sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                        sigbyr12 = sigbyr6*sigbyr6
                                        nrg_vdw = eps * (sigbyr12 - 2.0_DP*sigbyr6) ! eps was not multiplied by 4
                                        vdw_energy = vdw_energy + nrg_vdw
                                        IF (this_est_emax) up_nrg_vec(i) = nrg_vdw
                                END DO
                        ELSE IF (int_vdw_sum_style(this_box) == vdw_cut_shift) THEN
                                DO i = 1, n_vdw
                                        eps = ij_vdw_p_table(1,i) ! 4*epsilon
                                        negsigsq = ij_vdw_p_table(2,i) ! -(sigma**2)
                                        negsigbyr2 = negsigsq/rijsq_packed(i)
                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                        rterm = rterm + rterm*rterm
                                        nrg_vdw = eps * rterm
                                        negsigbyr2 = negsigsq/rcut_vdwsq(this_box)
                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                        rterm = rterm + rterm*rterm
                                        nrg_vdw = nrg_vdw - eps * rterm
                                        vdw_energy = vdw_energy + nrg_vdw
                                        IF (this_est_emax) up_nrg_vec(i) = nrg_vdw
                                END DO
                        ELSE IF (int_vdw_sum_style(this_box) == vdw_cut_switch) THEN
                                DO i = 1, n_vdw
                                        eps = ij_vdw_p_table(1,i) ! 4*epsilon
                                        negsigsq = ij_vdw_p_table(2,i) ! -(sigma**2)
                                        this_rijsq = rijsq_packed(i)
                                        negsigbyr2 = negsigsq/this_rijsq
                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                        rterm = rterm + rterm*rterm
                                        nrg_vdw = eps * rterm
                                        IF (this_rijsq .GE. ron_switch_sq(this_box)) THEN
                                                roffsq_rijsq = roff_switch_sq(this_box) - this_rijsq
                                                nrg_vdw = &
                                                        roffsq_rijsq*roffsq_rijsq * &
                                                        (switch_factor2(this_box)+2.0_DP*this_rijsq)*switch_factor1(this_box) * &
                                                        nrg_vdw
                                        END IF
                                        vdw_energy = vdw_energy + nrg_vdw
                                        IF (this_est_emax) up_nrg_vec(i) = nrg_vdw
                                END DO
                        ELSE IF (int_vdw_sum_style(this_box) == vdw_cut_shift_force) THEN
                                DO i = 1, n_vdw
                                        eps = ij_vdw_p_table(1,i) ! 4*epsilon
                                        negsigsq = ij_vdw_p_table(2,i) ! -(sigma**2)
                                        this_rijsq = rijsq_packed(i)
                                        negsigbyr2 = negsigsq/this_rijsq
                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                        rterm = rterm + rterm*rterm
                                        nrg_vdw = eps * rterm
                                        negsigbyr2 = negsigsq/rcut_vdwsq(this_box)
                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                        rterm2 = rterm * rterm
                                        rterm2 = rterm + 2.0_DP * rterm2
                                        rterm = rterm + rterm*rterm
                                        nrg_vdw = nrg_vdw - eps * rterm
                                        nrg_vdw = nrg_vdw - &
                                                (SQRT(this_rijsq) - rcut_vdw(this_box)) * &
                                                -6.0_DP * eps * rterm2 / rcut_vdw(this_box)
                                        vdw_energy = vdw_energy + nrg_vdw
                                        IF (this_est_emax) up_nrg_vec(i) = nrg_vdw
                                END DO
                        ELSE
                                DO i = 1, n_vdw
                                        !eps = ppvdwp_table(1,packed_types(i),itype,ibox) ! 4*epsilon
                                        !negsigsq = ppvdwp_table(2,packed_types(i),itype,ibox) ! -(sigma**2)
                                        eps = ij_vdw_p_table(1,i) ! 4*epsilon
                                        negsigsq = ij_vdw_p_table(2,i) ! -(sigma**2)
                                        negsigbyr2 = negsigsq/rijsq_packed(i)
                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                        rterm = rterm + rterm*rterm
                                        ! nrg_vdw = eps * rterm
                                        vdw_energy = vdw_energy + eps * rterm
                                        IF (this_est_emax) up_nrg_vec(i) = eps * rterm
                                END DO
                        END IF
                ELSE IF (int_vdw_style(this_box) == vdw_mie) THEN
                        !DIR$ ASSUME (n_vdw_p_list(this_box) .EQ. 2)
                        mie_vdw_p_table(1:4,1:n_vdw) = ppvdwp_table(1:4,PACK(jatomtype(1:vlen),vdw_mask(1:vlen)),itype,this_box)
                        DO i = 1, n_vdw
                                eps = mie_vdw_p_table(1,i) ! epsilon * mie_coeff
                                sigsq = mie_vdw_p_table(2,i) ! sigma**2
                                mie_n = mie_vdw_p_table(3,i) ! already halved
                                mie_m = mie_vdw_p_table(4,i) ! already halved
                                ! exponents were already halved so there's no need for SQRT
                                sigbyr2 = sigsq / rijsq_packed(i)
                                nrg_vdw = eps * &
                                        (sigbyr2**mie_n - sigbyr2**mie_m)
                                IF (int_vdw_sum_style(this_box) == vdw_cut_shift) THEN
                                        sigbyr2 = sigsq / rcut_vdw(this_box)
                                        nrg_vdw = nrg_vdw - eps * &
                                                (sigbyr2**mie_n - sigbyr2**mie_m)
                                END IF
                                vdw_energy = vdw_energy + nrg_vdw
                                IF (this_est_emax) up_nrg_vec(i) = nrg_vdw
                        END DO
                END IF
                !vdw_energy = SUM(nrg_vdw)
                IF (this_est_emax) up_nrg_vec(1:vlen) = UNPACK(up_nrg_vec(1:n_vdw), vdw_mask(1:vlen), zero_field(1:vlen))
        END IF
        ! Electrostatics
        IF (i_get_coul) THEN
                icharge_factor = icharge(ia) ! already multiplied by charge_factor
                IF (l_pack_separately .OR. .NOT. i_get_vdw) THEN
                        IF (l_notallcoul .AND. l_check_coul) THEN
                                coul_mask(1:vlen) = coul_mask(1:vlen) .AND. j_hascharge(1:vlen)
                        END IF
                        n_coul = COUNT(coul_mask(1:vlen))
                        jcharge_coul(1:n_coul) = PACK(jcharge(1:vlen), coul_mask(1:vlen))
                        rijsq_packed(1:n_coul) = PACK(rijsq(1:vlen), coul_mask(1:vlen))
                ELSE
                        n_coul = n_vdw
                END IF
                i_qq_energy = 0.0_DP
                DO i = 1, n_coul
                        rij = SQRT(rijsq_packed(i))
                        IF (int_charge_sum_style(this_box) == charge_ewald) THEN
                                nrg_coul = ERFC(alpha_ewald(this_box)*rij) / rij * jcharge_coul(i) ! omit icharge and charge factor for now
                        ELSE IF (int_charge_sum_style(this_box) == charge_dsf) THEN
                                nrg_coul = (dsf_factor2(this_box) * (rij - rcut_coul(this_box)) - &
                                        dsf_factor1(this_box) + &
                                        ERFC(alpha_dsf(this_box)*rij) / rij) * jcharge_coul(i)
                        !ELSE IF (int_charge_sum_style(this_box) == charge_cut) THEN
                        ELSE ! implies int_charge_sum_style(this_box) == charge_cut
                                nrg_coul = jcharge_coul(i) / rij
                        END IF
                        i_qq_energy = i_qq_energy + nrg_coul
                END DO
                qq_energy = qq_energy + i_qq_energy * icharge_factor
                IF (this_est_emax) THEN
                        up_nrg_vec = up_nrg_vec + UNPACK(jcharge_coul(1:n_coul)*icharge_factor/SQRT(rijsq_packed(1:n_coul)), coul_mask(1:vlen), zero_field(1:vlen))
                END IF
        END IF
        IF (this_est_emax) THEN
                Eij_max  = MAX(Eij_max, MAXVAL(up_nrg_vec))
        END IF
    END DO
  END SUBROUTINE Compute_Molecule_Nonbond_Inter_Energy_Vectorized_Widom


  SUBROUTINE Compute_MoleculeCollection_Nonbond_Inter_Energy(n_list,lm_list,is_list, &
    E_inter_vdw,E_inter_qq,overlap)
    !***************************************************************************
    ! This subroutine computes interatomic LJ and charge interactions as well as
    ! virials associated with these interactions.
    !
    ! CALLS
    !
    ! Minimum_Image_Separation
    ! Compute_MoleculePair_Energy
    ! Clean_Abort
    !
    ! CALLED BY
    !
    ! Reaction
    !
    ! Written by Jindal Shah on 12/07/07
    !***************************************************************************

    IMPLICIT NONE

!    !$ include 'omp_lib.h'

    INTEGER, INTENT(IN):: n_list ! number of molecules in the collection
    INTEGER, INTENT(IN):: lm_list(n_list) ! locates of each molecule
    INTEGER, INTENT(IN):: is_list(n_list) ! species of each molecule
    REAL(DP), INTENT(OUT) :: E_inter_vdw, E_inter_qq
    LOGICAL :: overlap
    !---------------------------------------------------------------------------

    INTEGER  :: ibox
    INTEGER  :: is, is2, is3 ! species index
    INTEGER  :: im, im2, im3 ! molecule index
    INTEGER  :: lm, lm2, lm3 ! molecule locate for molecule_list
    INTEGER  :: locate_1, locate_2 ! molecule locate for pair_nrg arrays

    REAL(DP) :: E12_vdw, E12_qq ! pairwise energy between molecules i,j
    REAL(DP) :: eps
    REAL(DP) :: rcom, rx, ry, rz

    LOGICAL :: get_interaction

    LOGICAL :: l_pair_store
    LOGICAL :: my_overlap, shared_overlap

    E_inter_vdw = 0.0_DP
    E_inter_qq = 0.0_DP
    overlap = .FALSE.
    my_overlap = .FALSE.
    shared_overlap = .FALSE.

    l_pair_store = .FALSE.
    IF (l_pair_nrg .AND. (.NOT. cbmc_flag)) l_pair_store = .TRUE.


    ibox = molecule_list(lm_list(1),is_list(1))%which_box
    DO im = 1, n_list
      lm = lm_list(im)
      is  = is_list(im)
      IF (molecule_list(lm,is)%which_box /= ibox) THEN
        err_msg = ""
        err_msg(1) = 'Collection of molecules must be in the same box'
        CALL Clean_Abort(err_msg,'Compute_MoleculeCollection_Nonbond_Inter_Energy')
      END IF

      IF (l_pair_store) CALL Get_Position_Alive(lm,is,locate_1)

      ! loop over other molecules in the collection
      im2loop: DO im2 = im+1, n_list
        lm2 = lm_list(im2)
        is2 = is_list(im2)

        IF (overlap) CYCLE

        IF (l_pair_store) THEN
           ! find out the location correspoding to lm2 in pair_nrg
           CALL Get_Position_Alive(lm2,is2,locate_2)

          ! reset pair energy
           pair_nrg_vdw(locate_1,locate_2) = 0.0_DP
           pair_nrg_vdw(locate_2,locate_1) = 0.0_DP

           pair_nrg_qq(locate_1,locate_2) = 0.0_DP
           pair_nrg_qq(locate_2,locate_1) = 0.0_DP
        END IF

        ! Determine if any atoms of these two molecules will interact
        CALL Check_MoleculePair_Cutoff(lm,is,lm2,is2,get_interaction,rcom,rx,ry,rz)

        IF (.NOT. get_interaction) CYCLE

        CALL Compute_MoleculePair_Energy(lm,is,lm2,is2,ibox,E12_vdw,E12_qq,overlap)

        IF (overlap) RETURN

        E_inter_vdw = E_inter_vdw + E12_vdw
        E_inter_qq  = E_inter_qq + E12_qq

     END DO im2loop 

      ! loop over molecules not in the collection
      speciesLoop: DO is2 = 1, nspecies

         !$OMP PARALLEL DO DEFAULT(SHARED) &
         !$OMP PRIVATE(im2,lm2,locate_2,get_interaction) &
         !$OMP PRIVATE(rx,ry,rz,rcom,E12_vdw,E12_qq) &
         !$OMP SCHEDULE(DYNAMIC) &
         !$OMP REDUCTION(+:E_inter_vdw,E_inter_qq) &
         !$OMP REDUCTION(.OR.:my_overlap)

         moleculeLoop: DO im2 = 1, nmols(is2,ibox)

            IF(shared_overlap) CYCLE moleculeLoop

            lm2 = locate(im2,is2,ibox)
            IF (.NOT. molecule_list(lm2,is2)%live) CYCLE moleculeLoop
            ! skip molecules that are in the collection
            DO im3 = 1, n_list
              lm3 = lm_list(im3)
              is3 = is_list(im3)
              IF (is2 == is3 .AND. lm2 == lm3) CYCLE moleculeLoop
            END DO

            ! reset pair energy, if storing energies
            IF (l_pair_store) THEN
               CALL Get_Position_Alive(lm2,is2,locate_2)

               pair_nrg_vdw(locate_1,locate_2) = 0.0_DP
               pair_nrg_vdw(locate_2,locate_1) = 0.0_DP

               pair_nrg_qq(locate_1,locate_2) = 0.0_DP
               pair_nrg_qq(locate_2,locate_1) = 0.0_DP
            END IF

            ! Determine if any atoms of these two molecules will interact
            CALL Check_MoleculePair_Cutoff(lm,is,lm2,is2,get_interaction,rcom,rx,ry,rz)

            IF (.NOT. get_interaction) CYCLE moleculeLOOP

            CALL Compute_MoleculePair_Energy(lm,is,lm2,is2,ibox,E12_vdw,E12_qq,my_overlap)

            IF (my_overlap) shared_overlap = .TRUE.

            E_inter_vdw = E_inter_vdw + E12_vdw
            E_inter_qq  = E_inter_qq + E12_qq

         END DO moleculeLoop
         !$OMP END PARALLEL DO

         IF(shared_overlap) THEN
            overlap = .TRUE.
            RETURN
         ENDIF

      END DO speciesLoop
    END DO

  END SUBROUTINE Compute_MoleculeCollection_Nonbond_Inter_Energy
  !-----------------------------------------------------------------------------

  SUBROUTINE Estimate_MoleculePair_Energy(im,is,jm,js,this_box,mpnrg,overlap)
          ! Only for use during CBMC trials
          !
          !
          ! molecule i is the "solute" molecule and molecule j
          !    is the "solvent" molecule
          IMPLICIT NONE
          INTEGER, INTENT(IN) :: im, is, jm, js, this_box
          REAL(DP), INTENT(OUT) :: mpnrg
          !----
          INTEGER :: ia, ja, bsolvent, bsolute
          TYPE(Atom_Class), POINTER :: these_atoms_i(:), these_atoms_j(:)
          REAL(DP) :: rcutsq_shift, rijsq_shift
          REAL(DP) :: ixp, iyp, izp, dx, dy, dz
          REAL(DP), DIMENSION(:,:,:), POINTER :: atompair_nrg_ptr
          LOGICAL :: overlap
          !
          overlap = .TRUE.
          mpnrg = 0.0_DP

          IF (im == widom_locate .AND. is == widom_species) THEN
                  these_atoms_i => widom_atoms
          ELSE
                  these_atoms_i => atom_list(1:natoms(is),im,is)
          END IF
          IF (jm == widom_locate .AND. js == widom_species) THEN
                  these_atoms_j => widom_atoms
          ELSE
                  these_atoms_j => atom_list(1:natoms(js),jm,js)
          END IF

          bsolvent = species_list(js)%solvent_base
          bsolute = species_list(is)%solute_base
          atompair_nrg_ptr => atompair_nrg_table(:, &
                  bsolvent+1:bsolvent+natoms(js), &
                  bsolute+1:bsolute+natoms(is), &
                  this_box)
          rcutsq_shift = rcut_cbmcsq(this_box) - rsq_shifter
          DO ia = 1, natoms(is)
                IF (.NOT. these_atoms_i(ia)%exist) CYCLE
                ixp = these_atoms_i(ia)%rxp
                iyp = these_atoms_i(ia)%ryp
                izp = these_atoms_i(ia)%rzp
                DO ja = 1, natoms(js)
                        IF (.NOT. these_atoms_j(ja)%exist) CYCLE
                        CALL Minimum_Image_Separation(this_box, &
                                ixp - these_atoms_j(ja)%rxp, &
                                iyp - these_atoms_j(ja)%ryp, &
                                izp - these_atoms_j(ja)%rzp, &
                                dx, dy, dz)
                        rijsq_shift = dx*dx + dy*dy + dz*dz - rsq_shifter
                        IF (rijsq_shift >= rcutsq_shift) CYCLE
                        IF (rijsq_shift < rsq_step) RETURN
                        mpnrg = mpnrg + atompair_nrg_ptr( &
                                INT(rijsq_shift/rsq_step), ja, ia)
                END DO
          END DO
          overlap = .FALSE.
          !
  END SUBROUTINE Estimate_MoleculePair_Energy

  SUBROUTINE Compute_MoleculePair_Energy(im,is,jm,js,this_box, &
    vlj_pair,vqq_pair,overlap,Eij_qq)
    !***************************************************************************
    ! The subroutine returns the interaction energy of the input molecule with
    ! another molecule. Thus, it computes the intermolecular vdw and
    ! electrostatic interactions.
    !
    ! CALLS
    !
    ! Minimum_Image_Separation
    ! Check_AtomPair_Cutoff
    ! Compute_AtomPair_Energy
    ! Get_Position_Alive
    !
    ! CALLED BY
    !
    ! Compute_Molecule_Nonbond_Inter_Energy
    ! Compute_System_Total_Energy
    !***************************************************************************

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: im, is, jm, js, this_box
    REAL(DP), INTENT(OUT) :: vlj_pair,vqq_pair
    LOGICAL, INTENT(OUT) :: overlap
    !---------------------------------------------------------------------------

    INTEGER :: ia, ja

    REAL(DP) :: rxijp, ryijp, rzijp, rxij, ryij, rzij, rijsq
    REAL(DP) :: Eij_intra_vdw, Eij_intra_qq, Eij_inter_vdw, Eij_inter_qq

    REAL(DP) :: Eij_qq_temp
    REAL(DP), OPTIONAL :: Eij_qq

    LOGICAL :: get_vdw, get_qq

    INTEGER :: locate_im, locate_jm

    TYPE(Atom_Class), POINTER :: these_atoms_i(:), these_atoms_j(:)
    REAL(DP), DIMENSION(:,:), POINTER :: smp_atompair_rsqmin
    INTEGER :: bsolvent
    LOGICAL ::  l_get_rij_min

    l_get_rij_min = .FALSE.

    IF (im == widom_locate .AND. is == widom_species) THEN
            these_atoms_i => widom_atoms
            IF (est_atompair_rminsq .AND. .NOT. cbmc_flag) THEN
                    l_get_rij_min = .TRUE.
                    bsolvent = species_list(js)%solvent_base
                    smp_atompair_rsqmin => swi_atompair_rsqmin( &
                            bsolvent+1:bsolvent+natoms(js),:) 
            END IF
    ELSE
            these_atoms_i => atom_list(:,im,is)
    END IF
    IF (jm == widom_locate .AND. js == widom_species) THEN
            ! not valid with est_atompair_rminsq==.TRUE.
            these_atoms_j => widom_atoms
    ELSE
            these_atoms_j => atom_list(:,jm,js)
    END IF

    vlj_pair = 0.0_DP
    vqq_pair = 0.0_DP
    Eij_qq_temp = 0.0_DP
    overlap = .FALSE.

    DO ia = 1, natoms(is)

      IF (.NOT. these_atoms_i(ia)%exist) CYCLE

      DO ja = 1, natoms(js)

        IF ( .NOT. these_atoms_j(ja)%exist) CYCLE

        ! Obtain the minimum image separation
        rxijp = these_atoms_i(ia)%rxp - these_atoms_j(ja)%rxp
        ryijp = these_atoms_i(ia)%ryp - these_atoms_j(ja)%ryp
        rzijp = these_atoms_i(ia)%rzp - these_atoms_j(ja)%rzp

        ! Now get the minimum image separation
        CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp,rxij,ryij,rzij)

        rijsq = rxij*rxij + ryij*ryij + rzij*rzij

        IF( rijsq < rcut_lowsq ) THEN
          overlap = .TRUE.
          RETURN
        END IF

        ! Now figure out what needs to be computed, then call pair_energy
        CALL Check_AtomPair_Cutoff(rijsq,get_vdw,get_qq,this_box)

        IF(cbmc_flag .AND. (.NOT. species_list(is)%L_Coul_CBMC)) THEN
          get_qq=.FALSE.
        ENDIF

        ! Compute vdw and q-q energy, if required
        IF (get_vdw .OR. get_qq) THEN

          CALL Compute_AtomPair_Energy(rxij,ryij,rzij,rijsq, &
               is,im,ia,js,jm,ja,get_vdw,get_qq, &
               Eij_intra_vdw,Eij_intra_qq,Eij_inter_vdw,Eij_inter_qq,Eij_qq)

          vlj_pair = vlj_pair + Eij_inter_vdw
          vqq_pair = vqq_pair + Eij_inter_qq

          IF (l_get_rij_min) THEN
                IF (rijsq<smp_atompair_rsqmin(ja,ia)) smp_atompair_rsqmin(ja,ia)=rijsq
          END IF

          IF (PRESENT(Eij_qq)) THEN
                  IF (widom_active .AND. .NOT. cbmc_flag) THEN
                          Eij_max = MAX(Eij_max,Eij_inter_vdw+Eij_qq)
                  END IF
                  Eij_qq_temp = Eij_qq_temp + Eij_qq
          END IF

        END IF
      END DO
    END DO
    IF (PRESENT(Eij_qq)) Eij_qq = Eij_qq_temp

    IF (l_pair_nrg) THEN
      IF ( .NOT. (cbmc_flag .OR. widom_active)) THEN
        ! if here then, there was no overlap between im and jm
        ! update the interaction energy of the pair
        ! first find out the position of each im in the pair interaction energy
        CALL Get_Position_Alive(im,is,locate_im)
        CALL Get_Position_Alive(jm,js,locate_jm)

        pair_nrg_vdw(locate_im,locate_jm) = vlj_pair
        pair_nrg_vdw(locate_jm,locate_im) = vlj_pair

        pair_nrg_qq(locate_im,locate_jm) = vqq_pair
        pair_nrg_qq(locate_jm,locate_im) = vqq_pair
      END IF
    END IF

  END SUBROUTINE Compute_MoleculePair_Energy
  !-----------------------------------------------------------------------------

  FUNCTION AtomPair_VdW_Energy_Vector(rijsq,itype,jtype,ibox)
      REAL(DP), DIMENSION(atompair_nrg_res) :: atompair_vdw_energy_vector
      REAL(DP), DIMENSION(atompair_nrg_res), INTENT(IN) :: rijsq
      INTEGER, INTENT(IN) :: itype, jtype, ibox
      REAL(DP) :: rcut_vdw
      REAL(DP) :: eps, sig, dEij_dr
      REAL(DP), DIMENSION(atompair_nrg_res) :: SigByR2, SigByR6, SigByR12
      REAL(DP) :: SigByR2_shift, SigByR6_shift, SigByR12_shift
      REAL(DP) :: roffsq_rijsq
      ! Mie potential
      REAL(DP) :: mie_coeff, mie_n, mie_m
      REAL(DP), DIMENSION(atompair_nrg_res) :: SigByR, SigByRn, SigByRm
      REAL(DP) :: SigByR_shift, SigByRn_shift, SigByRm_shift
      INTEGER :: i
      IF (int_vdw_style(ibox) == vdw_lj) THEN

           ! For now, assume all interactions are the same.
           ! Use the lookup table created in Compute_Nonbond_Table
           eps = vdw_param1_table(itype,jtype)
           sig = vdw_param2_table(itype,jtype)

           SigByR2 = (sig**2) / rijsq
           SigByR6 = SigByR2 * SigByR2 * SigByR2
           SigByR12 = SigByR6 * SigByR6

           ! use standard LJ potential
           atompair_vdw_energy_vector = 4.0_DP * eps * (SigByR12 - SigByR6)

           IF (int_vdw_sum_style(ibox) == vdw_cut_shift) THEN
                 ! shift the LJ potential
                 SigByR2_shift = sig**2/rcut_vdwsq(ibox)
                 SigByR6_shift = SigByR2_shift * SigByR2_shift * SigByR2_shift
                 SigByR12_shift = SigByR6_shift * SigByR6_shift

                 atompair_vdw_energy_vector = atompair_vdw_energy_vector &
                         - 4.0_DP * eps * (SigByR12_shift - SigByR6_shift)

           ELSE IF (int_vdw_sum_style(ibox) == vdw_cut_switch) THEN
                 DO i = 1, atompair_nrg_res
                         ! scale the LJ potential
                         IF ( rijsq(i) < ron_switch_sq(ibox) ) CYCLE
                         IF ( rijsq(i) > roff_switch_sq(ibox) ) THEN
                                 atompair_vdw_energy_vector(i) = 0.0_DP
                                 CYCLE
                         END IF
                         roffsq_rijsq = roff_switch_sq(ibox) - rijsq(i)
                         atompair_vdw_energy_vector(i) = &
                                 roffsq_rijsq*roffsq_rijsq &
                                 * (switch_factor2(ibox)+2.0_DP*rijsq(i))*switch_factor1(ibox) &
                                 * atompair_vdw_energy_vector(i)
                 END DO
           ELSE IF (int_vdw_sum_style(ibox) == vdw_charmm) THEN
                 ! use the form for modified LJ potential
                 atompair_vdw_energy_vector = eps * (SigByR12 - 2.0_DP * SigByR6)
           ELSE IF (int_vdw_sum_style(ibox) == vdw_cut_shift_force) THEN
                 ! apply the shifted-force LJ potential
                 ! u_sf(r) = u_lj(r) - u_lj(rc) - (r-rc)*du_lj/dr(rc)
                 SigByR2_shift = sig**2/rcut_vdwsq(ibox)
                 SigByR6_shift = SigByR2_shift * SigByR2_shift * SigByR2_shift
                 SigByR12_shift = SigByR6_shift * SigByR6_shift

                 rcut_vdw = SQRT(rcut_vdwsq(ibox))

                 dEij_dr = - 24.0_DP * eps * ( 2.0_DP * SigByR12_shift &
                                              - SigByR6_shift) / rcut_vdw

                 atompair_vdw_energy_vector = atompair_vdw_energy_vector &
                         - 4.0_DP * eps * (SigByR12_shift - SigByR6_shift) &
                         - (SQRT(rijsq) - rcut_vdw) * dEij_dr

           END IF

      ELSE IF (int_vdw_style(ibox) == vdw_mie) THEN
           eps = vdw_param1_table(itype,jtype)
           sig = vdw_param2_table(itype,jtype)
           mie_n = vdw_param3_table(itype,jtype) ! repulsive exponent
           mie_m = vdw_param4_table(itype,jtype) ! dispersive exponent

           rcut_vdw = SQRT(rcut_vdwsq(ibox))

           mie_coeff = mie_n/(mie_n-mie_m) * (mie_n/mie_m)**(mie_m/(mie_n-mie_m))

           SigByR = sig/SQRT(rijsq)
           SigByRn = SigByR ** mie_n
           SigByRm = SigByR ** mie_m
           atompair_vdw_energy_vector =  mie_coeff * eps * (SigByRn - SigByRm)
           !use cut-shift potential
           IF (int_vdw_sum_style(ibox) == vdw_cut_shift) THEN
                 SigByR_shift = sig/rcut_vdw
                 SigByRn_shift = SigByR_shift ** mie_n
                 SigByRm_shift = SigByR_shift ** mie_m
                 atompair_vdw_energy_vector =  atompair_vdw_energy_vector - mie_coeff * eps * (SigByRn_shift - SigByRm_shift)
           END IF

      END IF

  END FUNCTION AtomPair_VdW_Energy_Vector

  SUBROUTINE Compute_AtomPair_Energy(rxij,ryij,rzij,rijsq,is,im,ia,js,jm,ja, &
    get_vdw,get_qq,E_intra_vdw,E_intra_qq,E_inter_vdw,E_inter_qq,Eij_qq_o)

    ! Computes the vdw and q-q pair energy between atoms ia and ja of molecules
    ! im and jm and species is and js, given their separation rijsq. I have
    ! passed each component of separation but right now this is unnecessary.
    !
    ! LJ potential:
    !      Eij =  4*epsilon(i,j) * [ (sigma(i,j)/rij)^12 - (sigma(i,j)/rij)^6 ]
    !
    !

    ! It also computes the real space part of the Ewald sum if necessary.

    ! Called by:
    !   Compute_Atom_Nonbond_Energy
    !   Compute_Molecule_Nonbond_Energy
    ! Calls:
    !   Compute_AtomPair_Ewald_Real
  !----------------------------------------------------------------------------
    ! Passed to
    REAL(DP) :: rxij,ryij,rzij,rijsq
    INTEGER :: is,im,ia,js,jm,ja,ibox
    LOGICAL :: get_vdw,get_qq

    ! Returned
    REAL(DP) :: E_intra_vdw,E_intra_qq
    REAL(DP) :: E_inter_vdw,E_inter_qq

    ! Local
    ! LJ potential
    INTEGER :: itype, jtype
    REAL(DP) :: rij, rcut_vdw
    REAL(DP) :: eps, sig, Eij_vdw, dEij_dr
    REAL(DP) :: SigByR2, SigByR6, SigByR12
    REAL(DP) :: SigByR2_shift, SigByR6_shift, SigByR12_shift
    REAL(DP) :: roffsq_rijsq, roffsq_rijsq_sq, factor2, fscale
    ! Mie potential
    REAL(DP) :: mie_coeff, mie_n, mie_m
    REAL(DP) :: SigByR, SigByRn, SigByRm
    REAL(DP) :: SigByR_shift, SigByRn_shift, SigByRm_shift
    ! Coulomb potential
    REAL(DP) :: qi, qj, Eij_qq
    REAL(DP), OPTIONAL :: Eij_qq_o

    LOGICAL :: atom_i_exist, atom_j_exist

    E_intra_vdw = 0.0_DP
    E_intra_qq  = 0.0_DP
    E_inter_vdw = 0.0_DP
    E_inter_qq  = 0.0_DP
  !----------------------------------------------------------------------------
    ibox = molecule_list(im,is)%which_box

    IF (im == widom_locate .AND. is == widom_species) THEN
            atom_i_exist = widom_atoms(ia)%exist
    ELSE
            atom_i_exist = atom_list(ia,im,is)%exist
    END IF
    IF (jm == widom_locate .AND. js == widom_species) THEN
            atom_j_exist = widom_atoms(ja)%exist
    ELSE
            atom_j_exist = atom_list(ja,jm,js)%exist
    END IF

    ! If either atom is not yet present, then don't try to compute an energy
    ExistCheck: &
    IF (atom_i_exist .AND. atom_j_exist) THEN

         ! Determine atom type indices
         itype = nonbond_list(ia,is)%atom_type_number
         jtype = nonbond_list(ja,js)%atom_type_number

         VDW_calc: &
         IF (get_vdw .AND. itype /= 0 .AND. jtype /=0) THEN

              IF (int_vdw_style(ibox) == vdw_lj) THEN

                   ! For now, assume all interactions are the same.
                   ! Use the lookup table created in Compute_Nonbond_Table
                   eps = vdw_param1_table(itype,jtype)
                   sig = vdw_param2_table(itype,jtype)

                   ! Apply intramolecular scaling if necessary
                   IF (is == js .AND. im == jm) THEN
                     ! This controls 1-2, 1-3, and 1-4 interactions
                     eps = eps * vdw_intra_scale(ia,ja,is)
                   ENDIF

                   SigByR2 = (sig**2) / rijsq
                   SigByR6 = SigByR2 * SigByR2 * SigByR2
                   SigByR12 = SigByR6 * SigByR6

                   ! use standard LJ potential
                   Eij_vdw = 4.0_DP * eps * (SigByR12 - SigByR6)

                   IF (int_vdw_sum_style(ibox) == vdw_cut_shift) THEN
                         ! shift the LJ potential
                         SigByR2_shift = sig**2/rcut_vdwsq(ibox)
                         SigByR6_shift = SigByR2_shift * SigByR2_shift * SigByR2_shift
                         SigByR12_shift = SigByR6_shift * SigByR6_shift

                         Eij_vdw = Eij_vdw &
                                 - 4.0_DP * eps * (SigByR12_shift - SigByR6_shift)

                   ELSE IF (int_vdw_sum_style(ibox) == vdw_cut_switch) THEN
                         ! scale the LJ potential
                         IF ( rijsq < ron_switch_sq(ibox) ) THEN
                           fscale = 1.0_DP
                         ELSEIF ( rijsq <= roff_switch_sq(ibox) ) THEN
                           roffsq_rijsq = roff_switch_sq(ibox) - rijsq
                           roffsq_rijsq_sq = roffsq_rijsq * roffsq_rijsq
                           factor2 = switch_factor2(ibox) + 2.0_DP * rijsq
                           fscale = roffsq_rijsq_sq * factor2 * switch_factor1(ibox)
                           Eij_vdw = fscale * Eij_vdw
                         ELSE
                           fscale = 0.0_DP
                           Eij_vdw = 0.0_DP
                         END IF
                   ELSE IF (int_vdw_sum_style(ibox) == vdw_charmm) THEN
                         ! use the form for modified LJ potential
                         Eij_vdw = eps * (SigByR12 - 2.0_DP * SigByR6)
                   ELSE IF (int_vdw_sum_style(ibox) == vdw_cut_shift_force) THEN
                         ! apply the shifted-force LJ potential
                         ! u_sf(r) = u_lj(r) - u_lj(rc) - (r-rc)*du_lj/dr(rc)
                         SigByR2_shift = sig**2/rcut_vdwsq(ibox)
                         SigByR6_shift = SigByR2_shift * SigByR2_shift * SigByR2_shift
                         SigByR12_shift = SigByR6_shift * SigByR6_shift

                         Eij_vdw = Eij_vdw &
                                 - 4.0_DP * eps * (SigByR12_shift - SigByR6_shift)

                         rij = SQRT(rijsq)
                         rcut_vdw = SQRT(rcut_vdwsq(ibox))

                         dEij_dr = - 24.0_DP * eps * ( 2.0_DP * SigByR12_shift / rcut_vdw &
                                                      - SigByR6_shift / rcut_vdw )

                         Eij_vdw = Eij_vdw - (rij - rcut_vdw) * dEij_dr

                   END IF

              ELSE IF (int_vdw_style(ibox) == vdw_mie) THEN
                   eps = vdw_param1_table(itype,jtype)
                   sig = vdw_param2_table(itype,jtype)
                   mie_n = vdw_param3_table(itype,jtype) ! repulsive exponent
                   mie_m = vdw_param4_table(itype,jtype) ! dispersive exponent

                   ! Apply intramolecular scaling if necessary
                   IF (is == js .AND. im == jm) THEN
                     ! This controls 1-2, 1-3, and 1-4 interactions
                     eps = eps * vdw_intra_scale(ia,ja,is)
                   ENDIF

                   rij = SQRT(rijsq)
                   rcut_vdw = SQRT(rcut_vdwsq(ibox))

                   mie_coeff = mie_n/(mie_n-mie_m) * (mie_n/mie_m)**(mie_m/(mie_n-mie_m))

                   SigByR = sig/rij
                   SigByRn = SigByR ** mie_n
                   SigByRm = SigByR ** mie_m
                   Eij_vdw =  mie_coeff * eps * (SigByRn - SigByRm)
                   !use cut-shift potential
                   IF (int_vdw_sum_style(ibox) == vdw_cut_shift) THEN
                         SigByR_shift = sig/rcut_vdw
                         SigByRn_shift = SigByR_shift ** mie_n
                         SigByRm_shift = SigByR_shift ** mie_m
                         Eij_vdw =  Eij_vdw - mie_coeff * eps * (SigByRn_shift - SigByRm_shift)
                   END IF

              END IF

              IF (is == js .AND. im == jm) THEN
                   E_intra_vdw = Eij_vdw
                   E_inter_vdw = 0.0_DP
              ELSE
                   E_intra_vdw = 0.0_DP
                   E_inter_vdw = Eij_vdw
              ENDIF

         ENDIF VDW_calc

         qq_calc: IF (get_qq) THEN

              qi = nonbond_list(ia,is)%charge
              qj = nonbond_list(ja,js)%charge


              IF (int_charge_sum_style(ibox) == charge_ewald .AND. &
                      ( .NOT. igas_flag) ) THEN
                   ! Real space Ewald part
                   CALL Compute_AtomPair_Ewald_Real(ia,im,is,qi,ja,jm,js,qj, &
                        rijsq,E_intra_qq,E_inter_qq,ibox,Eij_qq_o)

                   ! self and reciprocal parts need to be computed as total energy
                   ! differences between original configuration and the perturbed
                   ! configuration.

              ELSEIF (int_charge_sum_style(ibox) == charge_dsf) THEN
                   CALL Compute_AtomPair_DSF_Energy(ia,im,is,qi,ja,jm,js,qj,rijsq,E_intra_qq,E_inter_qq,ibox,Eij_qq_o)

              ELSEIF (int_charge_sum_style(ibox) == charge_cut .OR. &
                  int_charge_sum_style(ibox) == charge_minimum .OR. igas_flag) THEN

                   Eij_qq = charge_factor*(qi*qj)/SQRT(rijsq)
                   ! Apply charge scaling for intramolecular energies
                   IF ( is == js .AND. im == jm ) THEN
                     Eij_qq = charge_intra_scale(ia,ja,is) * Eij_qq
                     E_intra_qq = E_intra_qq + Eij_qq
                   ELSE
                     E_inter_qq = E_inter_qq + Eij_qq
                   END IF
                   IF (PRESENT(Eij_qq_o)) Eij_qq_o = Eij_qq

             ENDIF

         ENDIF qq_calc

    ENDIF ExistCheck

  END SUBROUTINE Compute_AtomPair_Energy


SUBROUTINE Compute_AtomPair_DSF_Energy(ia,im,is,qi,ja,jm,js,qj,rijsq,E_intra_qq,E_inter_qq,ibox,Eij_qq_o)
USE Global_Variables
IMPLICIT NONE
INTEGER :: ia,im,is,ja,jm,js,ibox
REAL(DP), OPTIONAL :: Eij_qq_o
REAL(DP) :: qi,qj,rijsq,rij, Eij, qsc, E_intra_qq,E_inter_qq, cfqq, Eij_qq


      rij = SQRT(rijsq)
      cfqq = qi*qj*charge_factor
      Eij = dsf_factor2(ibox)*(rij-rcut_coul(ibox)) - dsf_factor1(ibox) + erfc(alpha_dsf(ibox)*rij)/rij
      Eij = Eij*cfqq

      IF (is==js .AND. im==jm) THEN
              qsc = charge_intra_scale(ia,ja,is)
              Eij_qq = cfqq/rij
              E_intra_qq = Eij - (1.0_DP-qsc)*Eij_qq
              E_inter_qq = 0.0_DP
              IF (PRESENT(Eij_qq_o)) Eij_qq_o = Eij_qq
      ELSE
              E_intra_qq = 0.0_DP
              E_inter_qq = Eij
              IF (PRESENT(Eij_qq_o)) Eij_qq_o = cfqq/rij
      END IF


END SUBROUTINE Compute_AtomPair_DSF_Energy



  !-----------------------------------------------------------------------------

  SUBROUTINE Compute_AtomPair_Ewald_Real(ia,im,is,qi,ja,jm,js,qj,rijsq, &
    E_intra_qq,E_inter_qq,ibox,Eij_qq_o)
  !-----------------------------------------------------------------------------
    ! Real space part of the Ewald sum between atoms ia and ja with
    ! charges qi and qj.
    !
    ! Miniumum image charges interact via Coulomb's law: qi*qj/rij
    !   * Intramolecular: scaled for 1-2, 1-3 and 1-4 interactions
    !
    ! The real space part of the periodic image charge is qi*qj/rij*erf
    !
    !
    ! CALLED BY:
    !
    ! Compute_AtomPair_Energy
  !-----------------------------------------------------------------------------
    ! Arguments
    INTEGER :: ia,im,is
    REAL(DP) :: qi
    INTEGER :: ja,jm,js
    REAL(DP) :: qj
    REAL(DP) :: rijsq
    REAL(DP) :: E_intra_qq, E_inter_qq, Eij_qq
    INTEGER :: ibox
    REAL(DP), OPTIONAL :: Eij_qq_o

    ! Local variables
    REAL(DP) :: rij,erf_val

    ibox = molecule_list(im,is)%which_box

    rij = SQRT(rijsq)
    ! May need to protect against very small rijsq
    erf_val = 1.0_DP - erfc(alpha_ewald(ibox) * rij)
    Eij_qq = qi * qj / rij * charge_factor

    ! Minimum image real space energy
    IF (is == js .AND. im == jm) THEN
       ! Intramolecular charge scaling
       E_intra_qq = charge_intra_scale(ia,ja,is) * Eij_qq
       E_inter_qq = 0.0_DP
    ELSE
       E_intra_qq = 0.0_DP
       E_inter_qq = Eij_qq
    ENDIF

    ! Periodic image real space energy
    E_inter_qq = E_inter_qq - erf_val * Eij_qq

    IF (PRESENT(Eij_qq_o)) Eij_qq_o = Eij_qq


  !-----------------------------------------------------------------------------
  CONTAINS

    FUNCTION erfc(x)
      !*************************************************************************
      !
      ! Calculate the complementary error function for a number
      !
      !*************************************************************************

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
  !-----------------------------------------------------------------------------

  END SUBROUTINE Compute_AtomPair_Ewald_Real

  !*****************************************************************************

  SUBROUTINE Ewald_Reciprocal_Lattice_Vector_Setup(this_box)
    !***************************************************************************
    ! This subroutine sets up the reciprocal lattice vector constants required in the reciprocal
    ! space summation. Note that these constants need to be recomputed every time a volume
    ! change move is attempted.
    ! Based on the APSS code, ewald_setup.f90
    !
    ! Added by Jindal Shah on 12/05/07
    !
    !***************************************************************************

    USE Type_Definitions
    USE Global_Variables

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
                IF ( kvecs > maxk) THEN
                   err_msg = ""
                   err_msg(1) = 'Total number of k vectors exceeded'
                   CALL Clean_Abort(err_msg,'Ewald_Reciprocal_Lattice_Vector_Setup')
                END IF

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

             END DO
          END DO
       END DO

    END IF



    ! nvecs points to where the next wave_vector should be written, i. e. it is too high by 1
    nvecs(this_box) = kvecs - 1
    ! Note that at this point we do not allocate the memory for cos_sum and sin_sum arrays.
    ! It will have to be decided by the maximum number of k vectors encountered in all the boxes


  END SUBROUTINE Ewald_Reciprocal_Lattice_Vector_Setup
  !*****************************************************************************

  SUBROUTINE Update_System_Ewald_Reciprocal_Energy(im,is,ibox, &
    move_flag,E_reciprocal)
    !***************************************************************************
    ! The subroutine computes the difference in Ewald reciprocal space energy
    ! for a given move.
    !
    ! We will develop this routine for a number of moves.
    !
    ! Translation of COM
    ! Rotation about COM
    ! Angle distortion
    ! Rigid dihedral rotation
    ! Molecule insertion
    ! Molecule deletion
    !***************************************************************************

    USE Type_Definitions
    USE Global_Variables

    IMPLICIT NONE

!    !$ include 'omp_lib.h'

    ! Arguments
    INTEGER, INTENT(IN) :: ibox   ! box index, 1...nbr_boxes
    INTEGER, INTENT(IN) :: is     ! species index, 1...nspecies
    INTEGER, INTENT(IN) :: im     ! molecule 'locate', index to atom_list
    INTEGER, INTENT(IN) :: move_flag

    ! Returns
    REAL(DP), INTENT(OUT) :: E_reciprocal

    ! Local variables
    REAL(DP) :: q
    INTEGER :: i, ia, jm, js

    REAL(DP) :: hdotr

    REAL(DP) :: cos_mol_im, cos_mol_im_o, sin_mol_im, sin_mol_im_o

    ! storage stuff
    INTEGER :: im_locate  ! index to cos_mol, sin_mol

    ! Initialize variables
    E_reciprocal = 0.0_DP

    !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
    cos_sum_old(1:nvecs(ibox),ibox) = cos_sum(1:nvecs(ibox),ibox)
    sin_sum_old(1:nvecs(ibox),ibox) = sin_sum(1:nvecs(ibox),ibox)
    !$OMP END PARALLEL WORKSHARE

    ! get the location of im for cos_mol, sin_mol arrays
    IF (is==1) THEN
      im_locate = im
    ELSE
      im_locate = SUM(max_molecules(1:is-1)) + im
    END IF

    IF ( move_flag == int_translation .OR. move_flag == int_rotation .OR. &
         move_flag == int_intra ) THEN

      ! only the particle coordinates change. Therefore, the contribution of
      ! cos(hdotr) and sin(hdotr) of the old coordinates will be subtracted
      ! off for each of reciprocal vectors and corresponding terms for the new
      ! coordinates are added.

      ! Note that the flag INTRA will refer to any of the moves that
      ! correspond to the intramolecular DOF change.

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP PRIVATE(i,ia,cos_mol_im,sin_mol_im) &
      !$OMP PRIVATE(cos_mol_im_o, sin_mol_im_o) &
      !$OMP PRIVATE(hdotr, q) &
      !$OMP SCHEDULE(STATIC) &
      !$OMP REDUCTION(+:E_reciprocal)
      DO i = 1, nvecs(ibox)

        cos_mol_im = 0.0_DP
        sin_mol_im = 0.0_DP

        DO ia = 1,natoms(is)
          ! compute the new hdotr
          hdotr = hx(i,ibox) * atom_list(ia,im,is)%rxp + &
                  hy(i,ibox) * atom_list(ia,im,is)%ryp + &
                  hz(i,ibox) * atom_list(ia,im,is)%rzp

          q = nonbond_list(ia,is)%charge
          cos_mol_im = cos_mol_im + q * DCOS(hdotr)
          sin_mol_im = sin_mol_im + q * DSIN(hdotr)
        END DO

        cos_mol_im_o = cos_mol(i,im_locate)
        sin_mol_im_o = sin_mol(i,im_locate)

        cos_sum(i,ibox) = cos_sum(i,ibox) + (cos_mol_im - cos_mol_im_o)
        sin_sum(i,ibox) = sin_sum(i,ibox) + (sin_mol_im - sin_mol_im_o)

        E_reciprocal = E_reciprocal + cn(i,ibox) &
                     * (cos_sum(i,ibox) * cos_sum(i,ibox) &
                     + sin_sum(i,ibox) * sin_sum(i,ibox))

        ! set the molecules cos and sin terms to the one calculated here
        cos_mol(i,im_locate) = cos_mol_im
        sin_mol(i,im_locate) = sin_mol_im

      END DO
      !$OMP END PARALLEL DO

    ELSE IF ( move_flag == int_deletion) THEN

      ! We need to subtract off the cos(hdotr) and sin(hdotr) for each of the
      ! k vectors.

      !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
      cos_sum(1:nvecs(ibox),ibox) = cos_sum(1:nvecs(ibox),ibox) &
                            - cos_mol(1:nvecs(ibox),im_locate)
      sin_sum(1:nvecs(ibox),ibox) = sin_sum(1:nvecs(ibox),ibox) &
                            - sin_mol(1:nvecs(ibox),im_locate)
      !$OMP END PARALLEL WORKSHARE

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP PRIVATE(i) &
      !$OMP SCHEDULE(STATIC) &
      !$OMP REDUCTION(+:E_reciprocal)
      DO i = 1, nvecs(ibox)

        E_reciprocal = E_reciprocal + cn(i,ibox) &
                     * ( cos_sum(i,ibox) * cos_sum(i,ibox) &
                       + sin_sum(i,ibox) * sin_sum(i,ibox) )

      END DO
      !$OMP END PARALLEL DO

    ELSE IF ( move_flag == int_insertion ) THEN

      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP PRIVATE(i, ia, hdotr, q) &
      !$OMP SCHEDULE(STATIC) &
      !$OMP REDUCTION(+:E_reciprocal)

      DO i = 1, nvecs(ibox)

        cos_mol(i,im_locate) = 0.0_DP
        sin_mol(i,im_locate) = 0.0_DP

        DO ia = 1, natoms(is)
          ! Compute the new hdotr vector
          hdotr = hx(i,ibox) * atom_list(ia,im,is)%rxp + &
                  hy(i,ibox) * atom_list(ia,im,is)%ryp + &
                  hz(i,ibox) * atom_list(ia,im,is)%rzp

          q = nonbond_list(ia,is)%charge
          cos_mol(i,im_locate) = cos_mol(i,im_locate) + q * DCOS(hdotr)
          sin_mol(i,im_locate) = sin_mol(i,im_locate) + q * DSIN(hdotr)
        END DO

        cos_sum(i,ibox) = cos_sum(i,ibox) &
                        + cos_mol(i,im_locate)
        sin_sum(i,ibox) = sin_sum(i,ibox) &
                        + sin_mol(i,im_locate)

        E_reciprocal = E_reciprocal + cn(i,ibox) &
                     * ( cos_sum(i,ibox) * cos_sum(i,ibox) &
                       + sin_sum(i,ibox) * sin_sum(i,ibox) )

      END DO
      !$OMP END PARALLEL DO
    END IF

    E_reciprocal = E_reciprocal * charge_factor

  END SUBROUTINE Update_System_Ewald_Reciprocal_Energy
  !*****************************************************************************

  SUBROUTINE Update_System_Ewald_Reciprocal_Energy_Widom(im,is,ibox, &
    E_reciprocal)
    !***************************************************************************
    ! The subroutine computes the difference in Ewald reciprocal space energy
    ! for a Widom insertion.
    !
    !***************************************************************************

    USE Type_Definitions
    USE Global_Variables

    IMPLICIT NONE

!    !$ include 'omp_lib.h'

    ! Arguments
    INTEGER, INTENT(IN) :: ibox   ! box index, 1...nbr_boxes
    INTEGER, INTENT(IN) :: is     ! species index, 1...nspecies
    INTEGER, INTENT(IN) :: im     ! molecule 'locate', index to atom_list

    ! Returns
    REAL(DP), INTENT(OUT) :: E_reciprocal

    ! Local variables
    REAL(DP), DIMENSION(natoms(is)) :: q, hdotr
    REAL(DP) :: cos_sum_i, sin_sum_i
    INTEGER :: i

    q = nonbond_list(1:natoms(is),is)%charge

    ! Initialize variables
    E_reciprocal = 0.0_DP

    DO i = 1, nvecs(ibox)

      hdotr = hx(i,ibox) * widom_atoms%rxp + &
              hy(i,ibox) * widom_atoms%ryp + &
              hz(i,ibox) * widom_atoms%rzp

      cos_sum_i = cos_sum(i,ibox) + DOT_PRODUCT(q, DCOS(hdotr))
      sin_sum_i = sin_sum(i,ibox) + DOT_PRODUCT(q, DSIN(hdotr))

      E_reciprocal = E_reciprocal + cn(i,ibox) &
                   * ( cos_sum_i * cos_sum_i &
                     + sin_sum_i * sin_sum_i )

    END DO

    E_reciprocal = E_reciprocal * charge_factor

  END SUBROUTINE Update_System_Ewald_Reciprocal_Energy_Widom
  !*****************************************************************************

  SUBROUTINE Compute_System_Self_Energy(this_box)
    !***************************************************************************
    ! This subroutine calculates the constant term that arises from particles
    ! interacting with themselves in the reciprocal space. The subroutine needs
    ! to be called only once as it is a constant term as long as the particles
    ! and their charges remain the same.
    !***************************************************************************

    USE Type_Definitions
    USE Global_Variables

    IMPLICIT NONE

    ! Arguments
    INTEGER :: this_box

    ! Returns
    ! energy(this_box)%self, global variable

    ! Local Variables
    INTEGER :: is,im, this_locate, ia
    REAL(DP) :: q, E_self

    E_self = 0.0_DP


   DO is = 1, nspecies
     imLOOP: DO im = 1, nmols(is,this_box)

       this_locate = locate(im,is,this_box)
       IF (.NOT. molecule_list(this_locate,is)%live) CYCLE imLOOP

       DO ia = 1, natoms(is)
         ! obtain the charge
         q = nonbond_list(ia,is)%charge
         E_self = E_self + q * q
       END DO
     END DO imLOOP
   END DO

    IF (int_charge_sum_style(this_box) == charge_ewald) THEN
           E_self = - E_self * charge_factor * alpha_ewald(this_box) / rootPI
           energy(this_box)%self = E_self
    ELSE IF (int_charge_sum_style(this_box) == charge_dsf) THEN
           E_self = E_self * (alpha_dsf(this_box) / rootPI + dsf_factor1(this_box)/2.0_DP)
           energy(this_box)%self = - E_self * charge_factor
    END IF

END SUBROUTINE Compute_System_Self_Energy


!*****************************************************************************

SUBROUTINE Compute_Molecule_Self_Energy(im,is,this_box,E_self)
  !***************************************************************************
  ! This subroutine calculates the self Ewald energy for the
  ! input molecule.
  !
  ! CALLED BY:
  !
  ! Chempot
  ! GEMC_Particle_Transfer
  ! Insertion
  ! Deletion
  ! Reaction
  !
  !***************************************************************************

  USE Type_Definitions
  USE Global_Variables

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN) :: im, is, this_box

  ! Returns
  REAL(DP), INTENT(OUT) :: E_self

  ! Local variables
  REAL(DP) :: q(natoms(is))

  ! Initialize variables
  E_self = 0.0_DP
  q = nonbond_list(1:natoms(is),is)%charge

  IF (int_charge_sum_style(this_box) == charge_ewald) THEN
          E_self = - DOT_PRODUCT(q,q) * charge_factor * alpha_ewald(this_box) / rootPI
  ELSE IF (int_charge_sum_style(this_box) == charge_dsf) THEN
          E_self = - DOT_PRODUCT(q,q) * (alpha_dsf(this_box) / rootPI + dsf_factor1(this_box)/2.0_DP) * charge_factor
  END IF

END SUBROUTINE Compute_Molecule_Self_Energy
  !*****************************************************************************

  SUBROUTINE Compute_System_Total_Energy(this_box,intra_flag,overlap)
    !***************************************************************************
    ! The subroutine calculates the total energy of a given box. The identity of
    ! the box is passed to the routine along with the intra_flag to indicate
    ! whether intramolecular computation is required. The flag will mostly be
    ! set to true except in the case of volume change move that is designed so
    ! that the intramolecular DOFs do not change.
    !***************************************************************************

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: this_box
    LOGICAL, INTENT(IN) :: intra_flag
    LOGICAL, INTENT(OUT) :: overlap

    !---------------------------------------------------------------------------

    INTEGER :: im, is, this_im, im_1, im_2, is_1, is_2, this_im_1, this_im_2

    REAL(DP) :: v_mol_bond, v_mol_angle, v_mol_dihedral, v_mol_improper
    REAL(DP) :: v_mol_intra_vdw, v_mol_intra_qq, v_mol_inter_qq
    REAL(DP) :: vlj_pair, vqq_pair, e_lrc
    REAL(DP) :: rcom, rx, ry, rz
    REAL(DP) :: E_inter_vdw, E_inter_qq
    REAL(DP) :: v_mol_self
    REAL(DP) :: rijsq
    REAL(DP) :: v_bond, v_angle, v_dihedral, v_intra, v_improper
    REAL(DP) :: v_intra_vdw, v_intra_qq, v_inter_qq

    LOGICAL :: get_interaction,intra_overlap

    INTEGER :: locate_1, locate_2
    LOGICAL :: l_pair_store, my_overlap, shared_overlap

    my_overlap = .FALSE.
    shared_overlap = .FALSE.
    overlap = .FALSE.

    ! Initialize the energies

    energy(this_box)%total = 0.0_DP
    energy(this_box)%inter = 0.0_DP
    energy(this_box)%inter_vdw = 0.0_DP
    energy(this_box)%lrc = 0.0_DP
    energy(this_box)%inter_q = 0.0_DP
    energy(this_box)%reciprocal = 0.0_DP
    ! Compute the intramolecular energy of the system if the flag is set.

    IF (intra_flag) THEN

       energy(this_box)%intra = 0.0_DP
       energy(this_box)%bond  = 0.0_DP
       energy(this_box)%angle = 0.0_DP
       energy(this_box)%dihedral = 0.0_DP
       energy(this_box)%improper = 0.0_DP
       energy(this_box)%intra_vdw = 0.0_DP
       energy(this_box)%intra_q = 0.0_DP
       energy(this_box)%self = 0.0_DP

       DO is = 1, nspecies
          v_intra = 0.0_DP
          v_bond= 0.0_DP
          v_angle= 0.0_DP
          v_dihedral= 0.0_DP
          v_improper = 0.0_DP
          v_intra_vdw= 0.0_DP
          v_intra_qq = 0.0_DP
          v_inter_qq = 0.0_DP
          !$OMP PARALLEL DO DEFAULT(SHARED) &
          !$OMP SCHEDULE(DYNAMIC) &
          !$OMP PRIVATE(im, this_im, v_mol_bond, v_mol_angle, v_mol_dihedral) &
          !$OMP PRIVATE(v_mol_improper,v_mol_intra_vdw,v_mol_intra_qq, v_mol_inter_qq, intra_overlap) &
          !$OMP REDUCTION(+:v_intra,v_bond, v_angle, v_dihedral,v_improper, v_intra_vdw, v_intra_qq, v_inter_qq)
          imLoop:DO im = 1, nmols(is,this_box)

             this_im = locate(im,is,this_box)
             IF (.NOT. molecule_list(this_im,is)%live) CYCLE imLoop

             IF (SHARED_OVERLAP) CYCLE imLOOP

             CALL Compute_Molecule_Bond_Energy(this_im,is,v_mol_bond)
             CALL Compute_Molecule_Angle_Energy(this_im,is,v_mol_angle)
             CALL Compute_Molecule_Dihedral_Energy(this_im,is,v_mol_dihedral)
             CALL Compute_Molecule_Improper_Energy(this_im,is,v_mol_improper)

             intra_overlap = .FALSE.
             CALL Compute_Molecule_Nonbond_Intra_Energy(this_im,is, &
                     v_mol_intra_vdw,v_mol_intra_qq,v_mol_inter_qq, &
                     intra_overlap)

             IF (intra_overlap) THEN
                SHARED_OVERLAP = .TRUE.
             END IF

             v_intra = v_intra + v_mol_bond + v_mol_angle &
                     + v_mol_dihedral + v_mol_improper &
                     + v_mol_intra_qq + v_mol_intra_vdw
             v_bond = v_bond + v_mol_bond
             v_angle = v_angle + v_mol_angle
             v_dihedral = v_dihedral + v_mol_dihedral
             v_improper = v_improper + v_mol_improper
             v_intra_vdw = v_intra_vdw + v_mol_intra_vdw
             v_intra_qq = v_intra_qq + v_mol_intra_qq
             ! electrostatic energy between this molecule and its periodic image
             v_inter_qq = v_inter_qq + v_mol_inter_qq

          END DO imLoop
          !$OMP END PARALLEL DO
          IF (SHARED_OVERLAP) THEN
             overlap = .TRUE.
             RETURN
          END IF

          energy(this_box)%intra = energy(this_box)%intra + v_intra
          energy(this_box)%bond = energy(this_box)%bond + v_bond
          energy(this_box)%angle = energy(this_box)%angle + v_angle
          energy(this_box)%dihedral = energy(this_box)%dihedral + v_dihedral
          energy(this_box)%improper = energy(this_box)%improper + v_improper
          energy(this_box)%intra_vdw = energy(this_box)%intra_vdw + v_intra_vdw
          energy(this_box)%intra_q   = energy(this_box)%intra_q   + v_intra_qq
          energy(this_box)%inter_q = energy(this_box)%inter_q + v_inter_qq
       END DO

    END IF

    ! Calculate the total intermolecular energy of the system. The calculation
    ! is divided into two parts. The first part computes the interaction between
    ! the molecules of the same species, while the second
    ! bit obtains the interaction between molecules of different species.

    l_pair_store = .FALSE.
    IF (l_pair_nrg .AND. (.NOT. cbmc_flag)) l_pair_store = .TRUE.

    DO is = 1, nspecies
       imLOOP1: DO im_1 = 1, nmols(is,this_box)
          this_im_1 = locate(im_1,is,this_box)
          IF (.NOT. molecule_list(this_im_1,is)%live) CYCLE imLOOP1

          IF (l_pair_store) CALL Get_Position_Alive(this_im_1, is, locate_1)

          E_inter_vdw = 0.0_DP
          E_inter_qq  = 0.0_DP

          !$OMP PARALLEL DO DEFAULT(SHARED) &
          !$OMP SCHEDULE(DYNAMIC) &
          !$OMP PRIVATE(im_2, this_im_2, locate_2, get_interaction) &
          !$OMP PRIVATE(rcom, rx, ry, rz, vlj_pair, vqq_pair) &
          !$OMP PRIVATE(my_overlap) &
          !$OMP REDUCTION(+: E_inter_vdw, E_inter_qq)

          imLOOP2: DO im_2 = im_1 + 1, nmols(is,this_box)
             this_im_2 = locate(im_2,is,this_box)
             IF (.NOT. molecule_list(this_im_2,is)%live) CYCLE imLOOP2

             IF (SHARED_OVERLAP) CYCLE imLOOP2

             IF (l_pair_store) THEN
                CALL Get_Position_Alive(this_im_2,is,locate_2)

                pair_nrg_vdw(locate_1,locate_2) = 0.0_DP
                pair_nrg_vdw(locate_2,locate_1) = 0.0_DP

                pair_nrg_qq(locate_1,locate_2) = 0.0_DP
                pair_nrg_qq(locate_2,locate_1) = 0.0_DP

             END IF

             CALL Check_MoleculePair_Cutoff(this_im_1,is,this_im_2,is,get_interaction, &
                  rcom,rx,ry,rz)

            ! rijsq = rcom * rcom

             IF (.NOT. get_interaction) CYCLE imLoop2
             ! Compute the intermolecular interactions between these two molecules

             CALL Compute_MoleculePair_Energy(this_im_1,is,this_im_2,is, &
                  this_box,vlj_pair,vqq_pair,my_overlap)

             !             IF (overlap) RETURN
             IF (my_overlap) THEN
                SHARED_OVERLAP = .true.
             END IF

             E_inter_vdw  = E_inter_vdw + vlj_pair
             E_inter_qq   = E_inter_qq  + vqq_pair

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
       imLOOP3: DO im_1 = 1, nmols(is_1,this_box)
          this_im_1 = locate(im_1,is_1,this_box)
          IF (.NOT. molecule_list(this_im_1,is_1)%live) CYCLE imLOOP3

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

             imLOOP4: DO im_2 = 1, nmols(is_2,this_box)
                this_im_2 = locate(im_2,is_2,this_box)
                IF (.NOT. molecule_list(this_im_2,is_2)%live) CYCLE imLOOP4

                IF (SHARED_OVERLAP) CYCLE imLOOP4

                IF (l_pair_store) THEN
                   CALL Get_Position_Alive(this_im_2,is_2,locate_2)

                   pair_nrg_vdw(locate_1,locate_2) = 0.0_DP
                   pair_nrg_vdw(locate_2,locate_1) = 0.0_DP

                   pair_nrg_qq(locate_1,locate_2) = 0.0_DP
                   pair_nrg_qq(locate_2,locate_1) = 0.0_DP

                END IF

                ! Check to see if the interaction needs to be computed between
                ! the molecules
                CALL Check_MoleculePair_Cutoff(this_im_1,is_1,this_im_2,is_2, &
                     get_interaction,rcom,rx,ry,rz)
!                rijsq = rcom * rcom

                IF (.NOT. get_interaction ) CYCLE imLOOP4

                ! Note that this call will modify the pair interaction energies
                ! if l_pair_nrg variable is .TRUE.

                CALL Compute_MoleculePair_Energy(this_im_1,is_1,this_im_2,is_2,&
                     this_box,vlj_pair,vqq_pair,my_overlap)

                IF (my_overlap) THEN
                   SHARED_OVERLAP = .true.
                END IF

                E_inter_vdw  = E_inter_vdw + vlj_pair
                E_inter_qq   = E_inter_qq  + vqq_pair

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

    energy(this_box)%inter = energy(this_box)%inter_q + energy(this_box)%inter_vdw

    ! Compute the reciprocal and self energy terms of the electrostatic energies if flag for Ewald is set.
    IF (int_charge_style(this_box) == charge_coul) THEN
       IF (int_charge_sum_style(this_box) == charge_ewald) THEN

            CALL Compute_System_Ewald_Reciprocal_Energy(this_box)

            energy(this_box)%inter = energy(this_box)%inter &
                                   + energy(this_box)%reciprocal

       END IF

       CALL Compute_System_Self_Energy(this_box)

       energy(this_box)%inter = energy(this_box)%inter &
                              + energy(this_box)%self

    END IF

    ! Long range correction if it is required
    IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
       CALL Compute_LR_Correction(this_box,e_lrc)
       ! add to the correction to the inter energy of the system
       energy(this_box)%lrc = e_lrc
       energy(this_box)%inter = energy(this_box)%inter + energy(this_box)%lrc
    END IF

    ! The total system energy. Note that intra_flag is not used for this
    ! calculation, beacuse, if the flag is true, we already computed the energy,
    ! if not we are using the old global energy (the routine
    ! did not modify the energy).
    energy(this_box)%total = energy(this_box)%intra + energy(this_box)%inter

  END SUBROUTINE Compute_System_Total_Energy
  !*****************************************************************************

  SUBROUTINE Compute_LR_Correction(this_box, e_lrc)
    !***************************************************************************
    ! The subroutine calculates the long range correction for the given box.
    !
    !***************************************************************************
    INTEGER, INTENT(IN) :: this_box
    REAL(DP), INTENT(OUT) :: e_lrc

    INTEGER ::  ia, ja, is, js
    REAL(DP) :: epsij, sigij, sigij2, sigij6, sigij12, mie_n, mie_m, mie_coeff
    REAL(DP) :: SigOverRcut, SigOverRn, SigOverRm
    REAL(DP) :: e_lrc_ia_ja

    e_lrc = 0.0_DP

    IF (int_vdw_style(this_box) == vdw_lj) THEN

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

    ELSE IF (int_vdw_style(this_box) == vdw_mie) THEN
      DO ia = 1, nbr_atomtypes
          e_lrc_ia_ja = 0.0_DP

          DO ja = 1, nbr_atomtypes

             epsij = vdw_param1_table(ia,ja)
             sigij = vdw_param2_table(ia,ja)
             mie_n = vdw_param3_table(ia,ja) ! repulsive exponent
             mie_m = vdw_param4_table(ia,ja) ! dispersive exponent
             mie_coeff = mie_n/(mie_n-mie_m)*(mie_n/mie_m)**(mie_m/(mie_n-mie_m))
             SigOverRcut = sigij/rcut_vdw(this_box)
             SigOverRn = SigOverRcut ** mie_n
             SigOverRm = SigOverRcut ** mie_m

             e_lrc_ia_ja = e_lrc_ia_ja + nint_beads(ja,this_box) * &
                   mie_coeff * epsij * rcut_vdw(this_box)**3.0_DP * ((SigOverRn/(3.0_DP-mie_n)) + &
                  (SigOverRm / (mie_m -3.0_DP)))
          END DO
          e_lrc = e_lrc + REAL( nint_beads(ia,this_box), DP ) * e_lrc_ia_ja

      END DO
      e_lrc = - 2.0_DP * PI * e_lrc/box_list(this_box)%volume
    END IF

  END SUBROUTINE Compute_LR_Correction

  !*****************************************************************************

  SUBROUTINE Check_MoleculePair_Cutoff(im_1,is_1,im_2,is_2,get_interaction, &
    rcom,rxcom,rycom,rzcom)

    REAL(DP) :: rxijp, ryijp, rzijp, rcom, rxcom, rycom, rzcom, rinteraction

    INTEGER :: this_box
    INTEGER :: im_1,is_1,im_2,is_2

    LOGICAL :: get_interaction

    TYPE(Molecule_Class), POINTER :: molecule_1, molecule_2

    IF (im_1 == widom_locate .AND. is_1 == widom_species) THEN
            molecule_1 => widom_molecule
    ELSE
            molecule_1 => molecule_list(im_1,is_1)
    END IF
    IF (im_2 == widom_locate .AND. is_2 == widom_species) THEN
            molecule_2 => widom_molecule
    ELSE
            molecule_2 => molecule_list(im_2,is_2)
    END IF

    ! Initially set the interaction to true.

    get_interaction = .TRUE.

    ! Figure out the box to be used later.

    this_box = molecule_1%which_box

    IF(int_vdw_sum_style(this_box) == vdw_minimum) RETURN

    ! Parent separation

    rxijp = molecule_1%xcom - molecule_2%xcom
    ryijp = molecule_1%ycom - molecule_2%ycom
    rzijp = molecule_1%zcom - molecule_2%zcom

    ! Compute the minimum image distance

    CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp,rxcom,rycom,rzcom)

    rcom = DSQRT(rxcom*rxcom + rycom*rycom + rzcom*rzcom)

    IF (CBMC_flag) THEN

       rinteraction = rcut_cbmc(this_box) + molecule_1%max_dcom &
            + molecule_2%max_dcom

       IF (rcom > rinteraction) get_interaction = .FALSE.

    ELSE

       rinteraction = rcut_max(this_box) + molecule_1%max_dcom &
            + molecule_2%max_dcom

       IF (rcom > rinteraction) get_interaction = .FALSE.

    END IF

 END SUBROUTINE Check_MoleculePair_Cutoff
!*******************************************************************************

!*******************************************************************************
 SUBROUTINE Check_AtomPair_Cutoff(rijsq,get_vdw,get_qq,this_box)

   INTEGER  :: this_box
   REAL(DP) :: rijsq
   LOGICAL  :: get_vdw, get_qq

   get_vdw = .FALSE.
   get_qq = .FALSE.

   VDW_Test2: IF (int_vdw_style(this_box) == vdw_none) THEN
      get_vdw = .FALSE.

   ELSEIF (int_vdw_style(this_box) == vdw_lj) THEN

      IF (CBMC_flag) THEN
         IF (rijsq <= rcut_cbmcsq(this_box)) THEN
            get_vdw = .TRUE.
         ELSE
            get_vdw = .FALSE.
         ENDIF
      ELSEIF (int_vdw_sum_style(this_box) == vdw_cut .OR. int_vdw_sum_style(this_box) &
              == vdw_cut_shift .OR. int_vdw_sum_style(this_box) == vdw_cut_tail &
              .OR. int_vdw_sum_style(this_box) == vdw_cut_shift_force ) THEN

         IF (rijsq <= rcut_vdwsq(this_box)) THEN
            get_vdw = .TRUE.
         ELSE
            get_vdw = .FALSE.
         ENDIF

      ELSEIF (int_vdw_sum_style(this_box) == vdw_minimum) THEN
         get_vdw = .TRUE.

      ELSEIF (int_vdw_sum_style(this_box) == vdw_charmm) THEN
         get_vdw = .TRUE.

      ELSEIF (int_vdw_sum_style(this_box) == vdw_cut_switch) THEN

         IF (rijsq <= roff_switch_sq(this_box)) THEN
            get_vdw = .TRUE.
         ELSE
            get_vdw = .FALSE.
         END IF

      ENDIF


   ELSEIF (int_vdw_style(this_box) == vdw_mie) THEN

      IF (CBMC_flag) THEN
         IF (rijsq <= rcut_cbmcsq(this_box)) THEN
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
      END IF

   ELSE
      err_msg = ""
      err_msg(1) = 'vdw_style must be NONE, LJ or Mie'
      CALL Clean_Abort(err_msg,'Compute_Atom_Nonbond_Energy')

   ENDIF VDW_Test2

   ! Charge sum tests
   IF (int_charge_style(this_box) == charge_none) THEN
      get_qq = .FALSE.
   ELSEIF (int_charge_style(this_box) == charge_coul) THEN

      IF (int_charge_sum_style(this_box) == charge_cut .OR. &
          int_charge_sum_style(this_box) == charge_ewald .OR. &
          int_charge_sum_style(this_box) == charge_dsf) THEN

         IF(CBMC_flag) THEN
            IF (rijsq <= rcut_cbmcsq(this_box)) THEN
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

 END SUBROUTINE Check_AtomPair_Cutoff

 SUBROUTINE Compute_System_Total_Force(this_box)

   !****************************************************************************
   ! The subroutine calculates the total forces of a given box.
   ! The identity of the box is passed to the routine.
   ! The forces are then used to compute the pressure tensor.
   !
   ! CALLS
   !
   ! CALLED BY
   !
   ! Volume_Change
   ! Main
   ! Write_Properties_Buffer
   !
   !****************************************************************************

   IMPLICIT NONE

   INTEGER, INTENT(IN) :: this_box

   !----------------------------------------------------------------------------

   INTEGER ::  is, im_1, im_2, is_1, is_2, this_im_1, this_im_2

   REAL(DP) :: rcom, rx, ry, rz, w_lrc

   REAL(DP),DIMENSION(3,3) :: tv_pair, tc_pair, w_inter_vdw, w_inter_charge

   LOGICAL :: get_interaction

   W_tensor_vdw(:,:,this_box) = 0.0_DP
   W_tensor_charge(:,:,this_box) = 0.0_DP
   W_tensor_recip(:,:,this_box) = 0.0_DP
   W_tensor_elec(:,:,this_box) =  0.0_DP

   DO is = 1, nspecies
      imLOOP1: DO im_1 = 1, nmols(is,this_box)
         this_im_1 = locate(im_1,is,this_box)
         IF (.NOT. molecule_list(this_im_1,is)%live) CYCLE imLOOP1

         w_inter_vdw(:,:) = 0.0_DP
         w_inter_charge(:,:) = 0.0_DP

         !$OMP PARALLEL DO DEFAULT(SHARED) &
         !$OMP SCHEDULE(DYNAMIC) &
         !$OMP PRIVATE(im_2, this_im_2, get_interaction) &
         !$OMP PRIVATE(rcom, rx, ry, rz, tv_pair, tc_pair) &
         !$OMP REDUCTION(+:w_inter_vdw, w_inter_charge)
         imLOOP2: DO im_2 = im_1 + 1, nmols(is,this_box)
            this_im_2 = locate(im_2,is,this_box)
            IF (.NOT. molecule_list(this_im_2,is)%live) CYCLE imLOOP2

            CALL Check_MoleculePair_Cutoff(this_im_1,is,this_im_2,is,get_interaction, &
                                   rcom,rx,ry,rz)

            IF (.NOT. Get_Interaction) CYCLE imLOOP2

            CALL Compute_MoleculePair_Force(this_im_1,is,this_im_2,is, &
                   this_box,tv_pair,tc_pair,rx,ry,rz)

            w_inter_vdw(:,:) = w_inter_vdw(:,:) + tv_pair(:,:)
            w_inter_charge(:,:) = w_inter_charge(:,:) + tc_pair(:,:)

         END DO imLOOP2
         !$OMP END PARALLEL DO

         W_tensor_vdw(:,:,this_box) = W_tensor_vdw(:,:,this_box) + w_inter_vdw(:,:)
         W_tensor_charge(:,:,this_box) = W_tensor_charge(:,:,this_box) + w_inter_charge(:,:)

      END DO imLOOP1
   END DO

   DO is_1 = 1, nspecies
      imLOOP3: DO im_1 = 1, nmols(is_1,this_box)
         this_im_1 = locate(im_1,is_1,this_box)
         IF (.NOT. molecule_list(this_im_1,is_1)%live) CYCLE imLOOP3

         DO is_2 = is_1 + 1, nspecies

            w_inter_vdw(:,:) = 0.0_DP
            w_inter_charge(:,:) = 0.0_DP

            !$OMP PARALLEL DO DEFAULT(SHARED) &
            !$OMP SCHEDULE(DYNAMIC) &
            !$OMP PRIVATE(im_2, this_im_2, get_interaction) &
            !$OMP PRIVATE(rcom, rx, ry, rz, tv_pair, tc_pair) &
            !$OMP REDUCTION(+:w_inter_vdw,w_inter_charge)
            imLOOP4: DO im_2 = 1, nmols(is_2,this_box)
               this_im_2 = locate(im_2,is_2,this_box)
               IF (.NOT. molecule_list(this_im_2,is_2)%live) CYCLE imLOOP4

               ! Check to see if the interaction needs to be computed between the molecules
               CALL Check_MoleculePair_Cutoff(this_im_1,is_1,this_im_2,is_2,get_interaction,rcom,rx,ry,rz)

               IF (.NOT. get_interaction ) CYCLE imLOOP4

               CALL Compute_MoleculePair_Force(this_im_1,is_1,this_im_2,is_2,this_box,tv_pair,tc_pair,rx,ry,rz)

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
       W_tensor_elec(:,:,this_box) =  W_tensor_recip(:,:,this_box)
    END IF

    IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN

       CALL Compute_LR_Force(this_box,w_lrc)
       virial(this_box)%lrc = w_lrc

    END IF

    W_tensor_elec(:,:,this_box) = (W_tensor_elec(:,:,this_box) + W_tensor_charge(:,:,this_box)) * charge_factor
    W_tensor_total(:,:,this_box) = W_tensor_vdw(:,:,this_box) + W_tensor_elec(:,:,this_box)

  END SUBROUTINE Compute_System_Total_Force
  !-----------------------------------------------------------------------------

  SUBROUTINE Compute_MoleculePair_Force(im,is,jm,js,this_box,tens_vdw,tens_charge,rabx,raby,rabz)
    !***************************************************************************
    ! The subroutine returns the interaction force of the input molecule with
    ! another molecule. Thus,
    ! it computes the intermolecular vdw and electrostatic interactions.
    !
    ! CALLED BY
    !
    ! Added by Jindal Shah on 12/10/07
    !***************************************************************************

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: im, is, jm, js, this_box
    !---------------------------------------------------------------------------

    INTEGER :: ia, ja

    REAL(DP) :: rxijp, ryijp, rzijp, rxij, ryij, rzij, rijsq, wij_vdw, wij_qq
    REAL(DP) :: rabx, raby, rabz
    REAL(DP),DIMENSION(3,3) :: tens_vdw, tens_charge

    REAL(DP) :: ffc, wxy, wxz, wyz

    LOGICAL :: get_vdw, get_qq

    tens_vdw(:,:) = 0.0_DP
    tens_charge(:,:) = 0.0_DP

    DO ia = 1, natoms(is)
       DO ja = 1, natoms(js)

          ! Obtain the minimum image separation
          rxijp = atom_list(ia,im,is)%rxp - atom_list(ja,jm,js)%rxp
          ryijp = atom_list(ia,im,is)%ryp - atom_list(ja,jm,js)%ryp
          rzijp = atom_list(ia,im,is)%rzp - atom_list(ja,jm,js)%rzp

          ! Now get the minimum image separation
          CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp, &
                  rxij,ryij,rzij)

          rijsq = rxij*rxij + ryij*ryij + rzij*rzij

          ! Now figure out what needs to be computed, then call pair_energy
          CALL Check_AtomPair_Cutoff(rijsq,get_vdw,get_qq,this_box)

          ! Compute vdw and q-q energy using if required
          IF (get_vdw .OR. get_qq) THEN

             CALL Compute_AtomPair_Force(rijsq,is,im,ia,js,jm,ja,&
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

  END SUBROUTINE Compute_MoleculePair_Force
  !-----------------------------------------------------------------------------
  SUBROUTINE Compute_AtomPair_Force &
       (rijsq,is,im,ia,js,jm,ja,get_vdw,get_qq,Wij_vdw,Wij_qq)

    ! LJ potential:  Wij = -rij/3 * d Eij / d rij.
    ! Use the virial in: P = NkT + < W >

    ! Computes the vdw and q-q pair force between atoms ia and ja of molecules
    ! im and jm and species is and js, given their separation rijsq. I have
    ! passed each component of separation but right now this is unnecessary.
    ! It also computes the real space part of the Ewald sum if necessary.

    ! Called by: Compute_System_Total_Force
  !-----------------------------------------------------------------------------
    ! Passed to
    REAL(DP) :: rxij,ryij,rzij,rijsq
    INTEGER :: is,im,ia,js,jm,ja
    LOGICAL :: get_vdw,get_qq

    ! Returned
    REAL(DP) :: Wij_vdw,Wij_qq

    ! Local
    INTEGER :: ibox
    ! LJ potential
    INTEGER :: itype, jtype
    REAL(DP) :: eps, sig, Eij_vdw
    REAL(DP) :: rij, rcut_vdw
    REAL(DP) :: SigByR2,SigByR6,SigByR12
    REAL(DP) :: SigByR2_shift,SigByR6_shift,SigByR12_shift
    REAL(DP) :: roffsq_rijsq, roffsq_rijsq_sq, factor2, fscale
    ! Mie potential
    REAL(DP) :: SigByR, SigByRn, SigByRm, mie_coeff, mie_n, mie_m
    ! Coulomb potential
    REAL(DP) :: qi, qj, erfc_val, prefactor
    REAL(DP) :: ewald_constant, exp_const

    Wij_vdw = 0.0_DP
    Wij_qq = 0.0_DP
  !-----------------------------------------------------------------------------
    ibox = molecule_list(im,is)%which_box

    ! If either atom is not yet present, then don't try to compute an energy
    ExistCheck: &
    IF (atom_list(ia,im,is)%exist .AND. atom_list(ja,jm,js)%exist) THEN

       ! Determine atom type indices
       itype = nonbond_list(ia,is)%atom_type_number
       jtype = nonbond_list(ja,js)%atom_type_number

       VDW_calc: &
       IF (get_vdw .AND. itype /= 0 .AND. jtype /=0) THEN

         IF (int_vdw_style(ibox) == vdw_lj) THEN
           ! For now, assume all interactions are the same.
           ! Use the lookup table created in Compute_Nonbond_Table
           eps = vdw_param1_table(itype,jtype)
           sig = vdw_param2_table(itype,jtype)

           SigByR2 = (sig**2.0_DP) / rijsq
           SigByR6  = SigByR2 * SigByR2 * SigByR2
           SigByR12 = SigByR6 * SigByR6

           ! Default potential for vdw_cut, vdw_cut_tail, vdw_cut_shift
           Wij_vdw = (24.0_DP * eps) * (2.0_DP*SigByR12 - SigByR6)

           IF (int_vdw_sum_style(ibox) == vdw_cut_switch) THEN
             IF (rijsq > ron_switch_sq(ibox) .AND. &
                 rijsq <= roff_switch_sq(ibox)) THEN
               roffsq_rijsq = roff_switch_sq(ibox) - rijsq
               roffsq_rijsq_sq = roffsq_rijsq * roffsq_rijsq
               factor2 = switch_factor2(ibox) + 2.0_DP * rijsq
               fscale = roffsq_rijsq_sq * factor2 * switch_factor1(ibox)
               Eij_vdw = 4.0_DP * eps * (SigByR12 - SigByR6)
               Eij_vdw = fscale * Eij_vdw
               Wij_vdw = fscale / 3.0_DP * Wij_vdw
               Wij_vdw = Wij_vdw + 8.0_DP * rijsq * rijsq * roffsq_rijsq &
                       * Eij_vdw * switch_factor1(ibox) / 3.0_DP
             ELSE IF (rijsq > roff_switch_sq(ibox)) THEN
               Wij_vdw = 0.0_DP
             END IF
           ELSEIF (int_vdw_sum_style(ibox) == vdw_charmm) THEN
             ! Use the CHARMM LJ potential
             Wij_vdw = (12.0_DP * eps) * (SigByR12 - SigByR6)
           ELSEIF (int_vdw_sum_style(ibox) == vdw_cut_shift_force) THEN
             ! shifted-force lj potential
             ! u_sf(r) = u_lj(r) - u_lj(rc) - (r-rc)*du_lj/dr(rc)
             SigByR2_shift = sig**2/rcut_vdwsq(ibox)
             SigByR6_shift = SigByR2_shift * SigByR2_shift * SigByR2_shift
             SigByR12_shift = SigByR6_shift * SigByR6_shift
             rij = SQRT(rijsq)
             rcut_vdw = SQRT(rcut_vdwsq(ibox))

             Wij_vdw = Wij_vdw &
                       - rij * (24.0_DP * eps) &
                       * (2.0_DP * SigByR12_shift / rcut_vdw &
                       - SigByR6_shift / rcut_vdw)
           END IF
         ELSE IF (int_vdw_style(ibox) == vdw_mie) THEN
           eps = vdw_param1_table(itype,jtype)
           sig = vdw_param2_table(itype,jtype)
           mie_n = vdw_param3_table(itype,jtype) ! repulsive exponent
           mie_m = vdw_param4_table(itype,jtype) ! dispersive exponent
           rij = SQRT(rijsq)

           mie_coeff = mie_n/(mie_n-mie_m)*(mie_n/mie_m)**(mie_m/(mie_n-mie_m))
           SigByR = sig/rij
           SigByRn = SigByR ** mie_n
           SigByRm = SigByR ** mie_m
           Wij_vdw = (mie_coeff * eps) *(mie_n * SigByRn - mie_m * SigByRm)

         ! Add other potential types here
         ENDIF

       ENDIF VDW_calc

       qq_calc: IF (get_qq) THEN

         qi = nonbond_list(ia,is)%charge
         qj = nonbond_list(ja,js)%charge

         rij = SQRT(rijsq)
         prefactor = qi * qj / rij
         IF (int_charge_sum_style(ibox) == charge_ewald) THEN
           ewald_constant = 2.0_DP * alpha_ewald(ibox) / rootPI
           exp_const = DEXP(-alpha_ewald(ibox)*alpha_ewald(ibox)*rijsq)
           ! May need to protect against very small rij
           erfc_val = erfc(alpha_ewald(ibox) * rij)
           Wij_qq = ( prefactor * erfc_val &
                  + qi * qj * ewald_constant * exp_const )

         ELSE IF (int_charge_sum_style(ibox) == charge_dsf) THEN

           Wij_qq = erfc(alpha_dsf(ibox)*rij)/(rijsq) + &
                    2.0_DP * alpha_dsf(ibox)/rootPI * &
                    DEXP(-alpha_dsf(ibox)*alpha_dsf(ibox) * rijsq) / rij - &
                    dsf_factor2(ibox)
           Wij_qq = qi*qj*Wij_qq*rij


         ELSE IF (int_charge_sum_style(ibox) == charge_cut) THEN
           Wij_qq = prefactor * charge_factor
         ENDIF

       ENDIF qq_calc

    ENDIF ExistCheck
!------------------------------------------------------------------------------
  CONTAINS

    FUNCTION erfc(x)
      !*************************************************************************
      !
      ! Calculate the complementary error function for  a number
      !
      !*************************************************************************

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

  END SUBROUTINE Compute_AtomPair_Force

  !-----------------------------------------------------------------------------
  SUBROUTINE Compute_LR_Force(this_box, w_lrc)
    !***************************************************************************
    ! The subroutine calculates the long range correction for the given box.
    !
    ! Called by
    !
    ! First written by Jindal Shah on 01/10/08
    !
    !
    !***************************************************************************

    INTEGER, INTENT(IN) :: this_box
    REAL(DP), INTENT(OUT) :: w_lrc

    INTEGER ::   is, js, ia, ja
    REAL(DP) :: mie_n, mie_m, mie_coeff
    REAL(DP) :: SigOverR, SigOverRn, SigOverRm

    REAL(DP) :: epsij, sigij
    REAL(DP) :: w_lrc_ia_ja

    w_lrc = 0.0_DP

    IF (int_vdw_style(this_box) == vdw_lj) THEN
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

    ELSEIF (int_vdw_style(this_box) == vdw_mie) THEN

      DO ia = 1, nbr_atomtypes

         w_lrc_ia_ja = 0.0_DP
         DO ja = 1, nbr_atomtypes
            epsij = vdw_param1_table(ia,ja)
            sigij = vdw_param2_table(ia,ja)
            mie_n = vdw_param3_table(ia,ja) ! repulsive exponent
            mie_m = vdw_param4_table(ia,ja) ! dispersive exponent
            mie_coeff = mie_n/(mie_n-mie_m)*(mie_n/mie_m)**(mie_m/(mie_n-mie_m))
            SigOverR = sigij/rcut_vdw(this_box)
            SigOverRn = SigOverR**mie_n
            SigOverRm = SigOverR**mie_m

            w_lrc_ia_ja = w_lrc_ia_ja + nint_beads(ja,this_box) * mie_coeff * epsij &
               *rcut3(this_box) * (mie_n/(mie_n-3.0_DP) * SigOverRn + mie_m/(3.0_DP-mie_m) * SigOverRm)

         END DO

         w_lrc = w_lrc + nint_beads(ia,this_box) * w_lrc_ia_ja
      END DO

      w_lrc =  2.0_DP / 3.0_DP * PI * w_lrc / box_list(this_box)%volume
    END IF
  END SUBROUTINE Compute_LR_Force

  SUBROUTINE Compute_System_Ewald_Reciprocal_Force(this_box)
    !***************************************************************************
    ! This subroutine computes the long range forces due to electrostatics
    !
    ! Based on APSS code reciprocal_ewald.f90
    !
    ! Added by Tom Rosch on 06/11/09
    ! (See Wheeler, Mol. Phys. 1997 Vol. 92 pg. 55)
    !
    !***************************************************************************

    USE Type_Definitions
    USE Global_Variables

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

       DO im = 1, nmols(is,this_box)

          this_locate = locate(im,is,this_box)
          IF( .NOT. molecule_list(this_locate,is)%live) CYCLE

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

  !-----------------------------------------------------------------------------

  SUBROUTINE Compute_Ring_Fragment_Energy(this_frag,this_im,is,this_box, &
    nrg_ring_frag)
    !***************************************************************************
    !
    ! This subroutine calculates the energy of a ring fragment in its old
    ! conformation
    !
    ! CALLED BY:
    !       fragment_growth.f90
    !
    ! CALLS :
    !       Compute_Molecule_Dihedral_Energy
    !       Compute_Molecule_Nonbond_Intra_Energy
    !
    ! Written by Jindal Shah on 10/03/09
    !
    !***************************************************************************

    INTEGER, INTENT(IN) :: this_frag, this_im, is
    REAL(DP), INTENT(OUT) :: nrg_ring_frag
    INTEGER :: this_box

    ! local variables

    INTEGER :: i, this_atom

!    REAL(DP) :: rcut_vdwsq_box, rcut_coulsq_box, alpha_ewald_box
    REAL(DP) :: e_dihed,  e_improper, nrg_vdw, nrg_qq, nrg_inter_qq

    LOGICAL :: intra_overlap
    LOGICAL :: exist_flag_old(natoms(is))

    nrg_ring_frag = 0.0_DP

    !!! Note for now, keep this_box == 1. For flexible ring
    ! molecule change this

    this_box = molecule_list(this_im,is)%which_box


    ! first store the exist flag of the molecule

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
    CALL Compute_Molecule_Nonbond_Intra_Energy(this_im,is,nrg_vdw,nrg_qq, &
         nrg_inter_qq,intra_overlap)


    nrg_ring_frag = nrg_vdw + nrg_qq + nrg_inter_qq + e_dihed + e_improper

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


  END SUBROUTINE Compute_Ring_Fragment_Energy

  !-----------------------------------------------------------------------------

  SUBROUTINE Check_System_Energy(ibox,check_inp)

     USE Global_Variables
     USE IO_Utilities

     INTEGER, INTENT(IN) :: ibox
     LOGICAL, OPTIONAL :: check_inp

     LOGICAL :: overlap, check

     TYPE(Energy_Class) :: e_check
     TYPE(Energy_Class) :: e_diff

     IF (present(check_inp)) THEN
        check = check_inp
     ELSE
        check = .TRUE.
     END IF

     IF (check) THEN
        e_check%total = energy(ibox)%total
        e_check%intra = energy(ibox)%intra
        e_check%inter = energy(ibox)%inter
        e_check%bond = energy(ibox)%bond
        e_check%angle = energy(ibox)%angle
        e_check%dihedral = energy(ibox)%dihedral
        e_check%improper = energy(ibox)%improper
        e_check%intra_vdw = energy(ibox)%intra_vdw
        e_check%intra_q = energy(ibox)%intra_q
        e_check%inter_vdw = energy(ibox)%inter_vdw
        e_check%inter_q = energy(ibox)%inter_q
        e_check%lrc = energy(ibox)%lrc
        e_check%reciprocal = energy(ibox)%reciprocal
        e_check%self = energy(ibox)%self
        e_diff%total = 0.0_DP
        e_diff%intra = 0.0_DP
        e_diff%inter = 0.0_DP
        e_diff%bond = 0.0_DP
        e_diff%angle = 0.0_DP
        e_diff%dihedral = 0.0_DP
        e_diff%improper = 0.0_DP
        e_diff%intra_vdw = 0.0_DP
        e_diff%intra_q = 0.0_DP
        e_diff%inter_vdw = 0.0_DP
        e_diff%inter_q = 0.0_DP
        e_diff%lrc = 0.0_DP
        e_diff%reciprocal = 0.0_DP
        e_diff%self = 0.0_DP
     END IF

     CALL Compute_System_Total_Energy(ibox,.TRUE.,overlap)

     IF (overlap) THEN
        ! overlap was detected between two atoms so abort the program
        err_msg = ''
        err_msg(1) = 'Atomic overlap in the configuration'
        CALL Clean_Abort(err_msg,'Check_System_Energy')
     END IF

     ! Compare recomputed energies to original
     IF (check) THEN
        e_diff%total = ABS(energy(ibox)%total - e_check%total)
        IF (ABS(energy(ibox)%total) > tiny_number) THEN
           e_diff%total = e_diff%total / energy(ibox)%total
        END IF
        e_diff%intra = ABS(energy(ibox)%intra - e_check%intra)
        IF (ABS(energy(ibox)%intra) > tiny_number) THEN
           e_diff%intra = e_diff%intra / energy(ibox)%intra
        END IF
        e_diff%inter = ABS(energy(ibox)%inter - e_check%inter)
        IF (ABS(energy(ibox)%inter) > tiny_number) THEN
           e_diff%inter = e_diff%inter / energy(ibox)%inter
        END IF
        e_diff%bond = ABS(energy(ibox)%bond - e_check%bond)
        IF (ABS(energy(ibox)%bond) > tiny_number) THEN
           e_diff%bond = e_diff%bond / energy(ibox)%bond
        END IF
        e_diff%angle = ABS(energy(ibox)%angle - e_check%angle)
        IF (ABS(energy(ibox)%angle) > tiny_number) THEN
           e_diff%angle = e_diff%angle / energy(ibox)%angle
        END IF
        e_diff%dihedral = ABS(energy(ibox)%dihedral - e_check%dihedral)
        IF (ABS(energy(ibox)%dihedral) > tiny_number) THEN
           e_diff%dihedral = e_diff%dihedral / energy(ibox)%dihedral
        END IF
        e_diff%improper = ABS(energy(ibox)%improper - e_check%improper)
        IF (ABS(energy(ibox)%improper) > tiny_number) THEN
           e_diff%improper = e_diff%improper / energy(ibox)%improper
        END IF
        e_diff%intra_vdw = ABS(energy(ibox)%intra_vdw - e_check%intra_vdw)
        IF (ABS(energy(ibox)%intra_vdw) > tiny_number) THEN
           e_diff%intra_vdw = e_diff%intra_vdw / energy(ibox)%intra_vdw
        END IF
        e_diff%intra_q = ABS(energy(ibox)%intra_q - e_check%intra_q)
        IF (ABS(energy(ibox)%intra_q) > tiny_number) THEN
           e_diff%intra_q = e_diff%intra_q / energy(ibox)%intra_q
        END IF
        e_diff%inter_vdw = ABS(energy(ibox)%inter_vdw - e_check%inter_vdw)
        IF (ABS(energy(ibox)%inter_vdw) > tiny_number) THEN
           e_diff%inter_vdw = e_diff%inter_vdw / energy(ibox)%inter_vdw
        END IF
        e_diff%inter_q = ABS(energy(ibox)%inter_q - e_check%inter_q)
        IF (ABS(energy(ibox)%inter_q) > tiny_number) THEN
           e_diff%inter_q = e_diff%inter_q / energy(ibox)%inter_q
        END IF
        e_diff%lrc = ABS(energy(ibox)%lrc - e_check%lrc)
        IF (ABS(energy(ibox)%lrc) > tiny_number) THEN
           e_diff%lrc = e_diff%lrc / energy(ibox)%lrc
        END IF
        e_diff%reciprocal = ABS(energy(ibox)%reciprocal - e_check%reciprocal)
        IF (ABS(energy(ibox)%reciprocal) > tiny_number) THEN
           e_diff%reciprocal = e_diff%reciprocal / energy(ibox)%reciprocal
        END IF
        e_diff%self = ABS(energy(ibox)%self - e_check%self)
        IF (ABS(energy(ibox)%self) > tiny_number) THEN
           e_diff%self = e_diff%self / energy(ibox)%self
        END IF
     END IF

     ! Write the recomputed energy components to log
     WRITE(logunit,*)
     WRITE(logunit,'(X,A,X,I1,T30,A20)',ADVANCE='NO') 'Energy components for box', ibox, 'kJ/mol-Extensive'
     IF (check) WRITE(logunit,'(X,A20)',ADVANCE='NO') 'Relative_Error'
     WRITE(logunit,*)
     WRITE(logunit,'(X,A)') '---------------------------------------------------------------------'
     WRITE(logunit,'(X,A,T30,F20.3)',ADVANCE='NO') 'Total system energy', energy(ibox)%total*atomic_to_kjmol
     IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%total
     WRITE(logunit,*)
     WRITE(logunit,'(X,A,T30,F20.3)',ADVANCE='NO') 'Intra molecular energy', energy(ibox)%intra*atomic_to_kjmol
     IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%intra
     WRITE(logunit,*)
     WRITE(logunit,'(3X,A,T30,F20.3)',ADVANCE='NO') 'Bond energy',energy(ibox)%bond*atomic_to_kjmol
     IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%bond
     WRITE(logunit,*)
     WRITE(logunit,'(3X,A,T30,F20.3)',ADVANCE='NO') 'Bond angle energy',energy(ibox)%angle*atomic_to_kjmol
     IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%angle
     WRITE(logunit,*)
     WRITE(logunit,'(3X,A,T30,F20.3)',ADVANCE='NO') 'Dihedral angle energy', energy(ibox)%dihedral*atomic_to_kjmol
     IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%dihedral
     WRITE(logunit,*)
     WRITE(logunit,'(3X,A,T30,F20.3)',ADVANCE='NO') 'Improper angle energy', energy(ibox)%improper*atomic_to_kjmol
     IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%improper
     WRITE(logunit,*)
     WRITE(logunit,'(3X,A,T30,F20.3)',ADVANCE='NO') 'Intra molecule vdw', energy(ibox)%intra_vdw*atomic_to_kjmol
     IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%intra_vdw
     WRITE(logunit,*)
     WRITE(logunit,'(3X,A,T30,F20.3)',ADVANCE='NO') 'Intra molecule q',energy(ibox)%intra_q*atomic_to_kjmol
     IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%intra_q
     WRITE(logunit,*)
     WRITE(logunit,'(X,A,T30,F20.3)',ADVANCE='NO') 'Inter molecular energy', energy(ibox)%inter*atomic_to_kjmol
     IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%inter
     WRITE(logunit,*)
     WRITE(logunit,'(3X,A,T30,F20.3)',ADVANCE='NO') 'Inter molecule vdw', energy(ibox)%inter_vdw*atomic_to_kjmol
     IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%inter_vdw
     WRITE(logunit,*)
     IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN
        WRITE(logunit,'(3X,A,T30,F20.3)',ADVANCE='NO') 'Long range correction', energy(ibox)%lrc*atomic_to_kjmol
        IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%lrc
        WRITE(logunit,*)
     END IF
     WRITE(logunit,'(3X,A,T30,F20.3)',ADVANCE='NO') 'Inter molecule q',energy(ibox)%inter_q*atomic_to_kjmol
     IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%inter_q
     WRITE(logunit,*)
     IF (int_charge_sum_style(ibox) == charge_ewald) THEN
        WRITE(logunit,'(3X,A,T30,F20.3)',ADVANCE='NO') 'Reciprocal ewald',energy(ibox)%reciprocal*atomic_to_kjmol
        IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%reciprocal
        WRITE(logunit,*)
        WRITE(logunit,'(3X,A,T30,F20.3)',ADVANCE='NO') 'Self ewald',energy(ibox)%self*atomic_to_kjmol
        IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%self
        WRITE(logunit,*)
     ELSE IF (int_charge_sum_style(ibox) == charge_dsf) THEN
        WRITE(logunit,'(X,A,T30,F20.3)',ADVANCE='NO') 'Self DSF',energy(ibox)%self*atomic_to_kjmol
        IF (check) WRITE(logunit,'(X,E20.3)',ADVANCE='NO') e_diff%self
        WRITE(logunit,*)
     END IF
     WRITE(logunit,'(X,A)') '---------------------------------------------------------------------'
     IF (int_charge_sum_style(ibox) == charge_ewald) &
        WRITE(logunit,'(3X,A,T33,I17)') 'Number of reciprocal vectors',nvecs(ibox)
     WRITE(logunit,*)

     IF (check) THEN
        energy(ibox)%total = e_check%total
        energy(ibox)%intra = e_check%intra
        energy(ibox)%inter = e_check%inter
        energy(ibox)%bond = e_check%bond
        energy(ibox)%angle = e_check%angle
        energy(ibox)%dihedral = e_check%dihedral
        energy(ibox)%improper = e_check%improper
        energy(ibox)%intra_vdw = e_check%intra_vdw
        energy(ibox)%intra_q = e_check%intra_q
        energy(ibox)%inter_vdw = e_check%inter_vdw
        energy(ibox)%inter_q = e_check%inter_q
        energy(ibox)%lrc = e_check%lrc
        energy(ibox)%reciprocal = e_check%reciprocal
        energy(ibox)%self = e_check%self
     END IF

  END SUBROUTINE Check_System_Energy

  !-----------------------------------------------------------------------------

  SUBROUTINE Compute_System_Ewald_Reciprocal_Energy(this_box)
    !***************************************************************************
    ! This subroutine computes the sin and cos sum terms for the calculation of
    ! reciprocal energy of the input box.
    !***************************************************************************

    USE Type_Definitions
    USE Global_Variables

    IMPLICIT NONE

!   !$ include 'omp_lib.h'

    ! Arguments
    INTEGER :: this_box

    ! Returns
    ! GLOBAL VARIABLE :: energy(this_box)%reciprocal

    ! Local Variables
    INTEGER :: i, is, im, ia, this_locate
    REAL(DP) :: un, const_val
    REAL(DP) :: charge, hdotr, E_reciprocal

    ! individual k-space vector stuff
    INTEGER ::  position
    INTEGER, ALLOCATABLE :: im_locate(:,:)

    ! openmp stuff
!   INTEGER :: omp_get_num_threads, omp_get_thread_num

    ! Initialize variables
    const_val = 1.0_DP/(2.0_DP * alpha_ewald(this_box) * alpha_ewald(this_box))
    E_reciprocal = 0.0_DP
    !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
    cos_sum(:,this_box) = 0.0_DP
    sin_sum(:,this_box) = 0.0_DP
    !$OMP END PARALLEL WORKSHARE

    ! Create an index, im_locate, for each live molecule in this_box
    ! im_locate will be used to access cos_mol and sin_mol
    ALLOCATE(im_locate(MAXVAL(max_molecules),nspecies))
    DO is = 1, nspecies
       DO im = 1, nmols(is,this_box)
          this_locate = locate(im,is,this_box)
          IF (.NOT. molecule_list(this_locate,is)%live) CYCLE

          ! create index
          IF (is == 1) THEN
             im_locate(im,is) = this_locate
          ELSE
             im_locate(im,is) = SUM(max_molecules(1:is-1)) + this_locate
          END IF

       END DO
    END DO

    ! Loop over each species, molecule
    DO is = 1, nspecies
       ! skip nonpolar species
       IF (.NOT. has_charge(is)) CYCLE

       DO im = 1, nmols(is,this_box)
          this_locate = locate(im,is,this_box) ! index to atom_list, molecule_list
          IF( .NOT. molecule_list(this_locate,is)%live) CYCLE

          position = im_locate(im,is) ! index to cos_mol, sin_mol

          !$OMP PARALLEL DO DEFAULT(SHARED) &
          !$OMP PRIVATE(i,ia,hdotr,charge) &
          !$OMP SCHEDULE(STATIC)

          ! loop over all the k vectors of this box
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

             cos_sum(i,this_box) = cos_sum(i,this_box) &
                                 + cos_mol(i,position)
             sin_sum(i,this_box) = sin_sum(i,this_box) &
                                 + sin_mol(i,position)
          END DO

          !$OMP END PARALLEL DO

       END DO
    END DO

    ! At the end of all the loops we have computed cos_sum, sin_sum, cos_mol and
    ! sin_mol for each of the k-vectors. Now let us calculate the reciprocal
    ! space energy

    !$OMP PARALLEL DO DEFAULT(SHARED) &
    !$OMP PRIVATE(i,  un)  &
    !$OMP SCHEDULE(STATIC) &
    !$OMP REDUCTION(+:E_reciprocal)

    DO i = 1, nvecs(this_box)
       un =  cos_sum(i,this_box) * cos_sum(i,this_box) &
          +  sin_sum(i,this_box) * sin_sum(i,this_box)
       E_reciprocal = E_reciprocal + Cn(i,this_box) * un
    END DO

    !$OMP END PARALLEL DO

    energy(this_box)%reciprocal = E_reciprocal * charge_factor

  END SUBROUTINE Compute_System_Ewald_Reciprocal_Energy

  !-----------------------------------------------------------------------------

END MODULE Energy_Routines
