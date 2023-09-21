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
  USE Internal_Coordinate_Routines
 !$  USE OMP_LIB

  IMPLICIT NONE


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
                        // TRIM(Int_To_String(bond_list(ib,is)%atom(1))) // ' and ' &
                        // TRIM(Int_To_String(bond_list(ib,is)%atom(2))) &
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
                        // TRIM(Int_To_String(angle_list(ia,is)%atom(1))) // ' and ' &
                        // TRIM(Int_To_String(angle_list(ia,is)%atom(2))) // ' and ' &
                        // TRIM(Int_To_String(angle_list(ia,is)%atom(3))) &
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
    INTEGER :: idihed, idihed_rb, atom1, atom2, atom3, atom4
    REAL(DP) :: a0,a1,a2,a3,a4,a5,a6,a7,a8,edihed,phi,cosphi,r12dn,twophi,threephi
    REAL(DP) :: cosphi_vec(0:5)

    TYPE(Atom_Class), POINTER, CONTIGUOUS :: these_atoms(:)

    IF (widom_active) THEN
            these_atoms => widom_atoms
    ELSE
            these_atoms => atom_list(:,molecule,species)
    END IF

    energy_dihed = 0.0_DP
    cosphi_vec = 0.0_DP
    cosphi_vec(0) = 1.0_DP
    DO idihed_rb = 1, species_list(species)%ndihedrals_rb
        IF (.NOT. ALL(these_atoms(dihedral_list(idihed_rb,species)%atom)%exist)) CYCLE
        CALL Get_Dihedral_Angle_COS(idihed_rb,molecule,species,cosphi_vec(1))
        cosphi_vec(2) = cosphi_vec(1)*cosphi_vec(1)
        cosphi_vec(3) = cosphi_vec(1)*cosphi_vec(2)
        cosphi_vec(4) = cosphi_vec(2)*cosphi_vec(2)
        cosphi_vec(5) = cosphi_vec(2)*cosphi_vec(3)
        energy_dihed = energy_dihed + DOT_PRODUCT(cosphi_vec,dihedral_list(idihed_rb,species)%rb_c)
    END DO
    DO idihed=idihed_rb, species_list(species)%ndihedrals_energetic
       ! Verify that the atoms of this dihedral exist. This is required
       ! for CBMC moves in which only a part of the molecule is present in
       ! the simulation
       IF (.NOT. ALL(these_atoms(dihedral_list(idihed,species)%atom)%exist)) CYCLE
       CALL Get_Dihedral_Angle_COS(idihed,molecule,species,cosphi,r12dn)
       phi = SIGN(ACOS(cosphi), r12dn)
       SELECT CASE(dihedral_list(idihed,species)%int_dipot_type)
       CASE(int_opls)
          ! Note: this case should never occur now that all OPLS dihedrals are internally converted to RB torsions
          ! However, I left it here anyway just in case that changes somehow.
          IF (.NOT. ALL(these_atoms(dihedral_list(idihed,species)%atom)%exist)) CYCLE
          a0 = dihedral_list(idihed,species)%dihedral_param(1)
          a1 = dihedral_list(idihed,species)%dihedral_param(2)
          a2 = dihedral_list(idihed,species)%dihedral_param(3)
          a3 = dihedral_list(idihed,species)%dihedral_param(4)

          twophi = 2.0_DP*phi
          threephi = 3.0_DP*phi
          edihed =  a0 + a1*(1.0_DP+COS(phi)) + &
               a2*(1.0_DP-COS(twophi)) + a3*(1.0_DP+COS(threephi))
       CASE(int_charmm)
          ! Check to see if the atoms of this dihedral exists. This is required
          ! for CBMC moves in which only a part of the molecule is present in
          ! the simulation
          IF (.NOT. ALL(these_atoms(dihedral_list(idihed,species)%atom)%exist)) CYCLE

          a0 = dihedral_list(idihed,species)%dihedral_param(1)
          a1 = dihedral_list(idihed,species)%dihedral_param(2)
          a2 = dihedral_list(idihed,species)%dihedral_param(3)

          edihed = a0 * (1.0_DP + DCOS(a1*phi - a2))

       CASE(int_harmonic)
          IF (.NOT. ALL(these_atoms(dihedral_list(idihed,species)%atom)%exist)) CYCLE

          a0 = dihedral_list(idihed,species)%dihedral_param(1)
          a1 = dihedral_list(idihed,species)%dihedral_param(2)

          IF(a1 .GT. 0.0_DP .AND. phi .LT.0) phi = phi + twoPI

          edihed = a0 * (phi - a1)**2

          ! Add more potential functions here.
       END SELECT
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
                   rxij = these_atoms_j(ja)%rp(1) &
                        - these_atoms_i(ia)%rp(1)
                   ryij = these_atoms_j(ja)%rp(2) &
                        - these_atoms_i(ia)%rp(2)
                   rzij = these_atoms_j(ja)%rp(3) &
                        - these_atoms_i(ia)%rp(3)

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
                rxijp = these_atoms_j(ja)%rp(1) &
                      - these_atoms_i(ia)%rp(1)
                ryijp = these_atoms_j(ja)%rp(2) &
                      - these_atoms_i(ia)%rp(2)
                rzijp = these_atoms_j(ja)%rp(3) &
                      - these_atoms_i(ia)%rp(3)

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
                rxij = these_atoms(ja)%rp(1) - these_atoms(ia)%rp(1)
                ryij = these_atoms(ja)%rp(2) - these_atoms(ia)%rp(2)
                rzij = these_atoms(ja)%rp(3) - these_atoms(ia)%rp(3)
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
!                cp(1) = widom_atoms(ia)%rp(1)
!                cp(2) = widom_atoms(ia)%rp(2)
!                cp(3) = widom_atoms(ia)%rp(3)
!                this_box = widom_molecule%which_box
!        ELSE
!                cp(1) = atom_list(ia,im,is)%rp(1)
!                cp(2) = atom_list(ia,im,is)%rp(2)
!                cp(3) = atom_list(ia,im,is)%rp(3)
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
!                                        rxijp = atom_ptr%rp(1) - cp(1)
!                                        ryijp = atom_ptr%rp(2) - cp(2)
!                                        rzijp = atom_ptr%rp(3) - cp(3)
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
          INTEGER :: coord_limit(4)

          REAL(SP) :: cp_sp(3)
          !
          overlap = .TRUE.
          E_inter_vdw = 0.0_DP
          E_inter_qq = 0.0_DP
          Eij_qq_temp = 0.0_DP
          IF (widom_active) THEN
                  cp(1) = widom_atoms(ia)%rp(1)
                  cp(2) = widom_atoms(ia)%rp(2)
                  cp(3) = widom_atoms(ia)%rp(3)
                  this_box = widom_molecule%which_box
          ELSE
                  cp(1) = atom_list(ia,im,is)%rp(1)
                  cp(2) = atom_list(ia,im,is)%rp(2)
                  cp(3) = atom_list(ia,im,is)%rp(3)
                  this_box = molecule_list(im,is)%which_box
          END IF
          IF (cbmc_flag) THEN
                  cp_sp = REAL(cp,SP)
                  cell_coords = IDNINT(MATMUL(REAL(box_list(this_box)%cell_length_inv,SP),cp))
                  coord_limit = UBOUND(cbmc_cell_n_interact)
                  IF (ANY(ABS(cell_coords) > coord_limit(1:3))) THEN
                          !WRITE(*,*) cell_coords, coord_limit(1:3), cp_sp, MATMUL(REAL(box_list(this_box)%cell_length_inv,SP),cp)
                          !WRITE(*,*) cp
                          DO i = 1, 3
                                IF (cell_coords(i) < LBOUND(cbmc_cell_n_interact,i)) THEN
                                        cell_coords(i) = cell_coords(i) + SIZE(cbmc_cell_n_interact,i)
                                        cp_sp(i) = cp_sp(i) + box_list(this_box)%length(i,i)*0.5_DP
                                ELSE IF (cell_coords(i) > UBOUND(cbmc_cell_n_interact,i)) THEN
                                        cell_coords(i) = cell_coords(i) - SIZE(cbmc_cell_n_interact,i)
                                        cp_sp(i) = cp_sp(i) - box_list(this_box)%length(i,i)*0.5_DP
                                ELSE
                                        CYCLE
                                END IF
                                !WRITE(*,*) cp_sp(i), box_list(this_box)%length(i,i)*0.5_DP
                          END DO
                  END IF
                  ! next line includes qq energy along with vdw
                  E_inter_vdw = REAL(Compute_Cell_List_CBMC_nrg(cp_sp,cell_coords,ia,is,this_box),DP)
                  overlap = .FALSE.
                  RETURN
                  !cell_length_inv_ptr => cell_length_inv_cbmc(:,:,this_box)
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
                                                dxp = atom_ptr%rp(1) - cp(1)
                                                dyp = atom_ptr%rp(2) - cp(2)
                                                dzp = atom_ptr%rp(3) - cp(3)
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
                                                dxp = atom_ptr%rp(1) - cp(1)
                                                dyp = atom_ptr%rp(2) - cp(2)
                                                dzp = atom_ptr%rp(3) - cp(3)
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
             rxij = these_atoms(ia)%rp(1) - these_atoms(ja)%rp(1)
             ryij = these_atoms(ia)%rp(2) - these_atoms(ja)%rp(2)
             rzij = these_atoms(ia)%rp(3) - these_atoms(ja)%rp(3)

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


    this_box = molecule_list(im,is)%which_box
    IF (l_vectorized) THEN
            atompairdim = natoms(is)*MAXVAL(natoms)
            mol_dim = MAXVAL(nmols(:,this_box))
            IF (MOD(atompairdim,8) .NE. 0) atompairdim = 8*(atompairdim/8+1)
            IF (MOD(mol_dim,8) .NE. 0) mol_dim = 8*(mol_dim/8+1)
            CALL Compute_Molecule_Nonbond_Inter_Energy_Vectorized(im,is, &
                    E_inter_vdw, E_inter_qq, overlap)
            RETURN
    END IF

    E_inter_vdw = 0.0_DP
    E_inter_qq = 0.0_DP
    overlap = .FALSE.
    my_overlap = .FALSE.
    shared_overlap = .FALSE.


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

!molecule_priority    hardcore_max_r = widom_molecule%rcom(4) + rcut_low
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
!molecule_priority          ELSE IF (rcom - molecule_list(this_locate,ispecies)%rcom(4) > hardcore_max_r) THEN
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
                                                live_atom_rsp(ja,jm,1,js,ibox) = atom_list(ja,jm+cjnmols,js)%rp(1)
                                                live_atom_rsp(ja,jm,2,js,ibox) = atom_list(ja,jm+cjnmols,js)%rp(2)
                                                live_atom_rsp(ja,jm,3,js,ibox) = atom_list(ja,jm+cjnmols,js)%rp(3)
                                                live_atom_exist(ja,jm,js,ibox) = atom_list(ja,jm+cjnmols,js)%exist
                                                !live_atom_list(ja,jm,js,ibox)%rp(1) = atom_list(ja,jm+cjnmols,js)%rp(1)
                                                !live_atom_list(ja,jm,js,ibox)%rp(2) = atom_list(ja,jm+cjnmols,js)%rp(2)
                                                !live_atom_list(ja,jm,js,ibox)%rp(3) = atom_list(ja,jm+cjnmols,js)%rp(3)
                                                !live_atom_list(ja,jm,js,ibox)%exist = atom_list(ja,jm+cjnmols,js)%exist
                                        END DO
                                END DO
                                !$OMP END DO SIMD
                                !$OMP DO SIMD
                                DO jm = 1, jnlive
                                        live_xcom(jm,js,ibox) = molecule_list(jm+cjnmols,js)%rcom(1)
                                        live_ycom(jm,js,ibox) = molecule_list(jm+cjnmols,js)%rcom(2)
                                        live_zcom(jm,js,ibox) = molecule_list(jm+cjnmols,js)%rcom(3)
                                        live_max_dcom(jm,js,ibox) = molecule_list(jm+cjnmols,js)%rcom(4)
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
                                                live_atom_rsp(ja,jm,1,js,ibox) = j_atom_list(ja,jm)%rp(1)
                                                live_atom_rsp(ja,jm,2,js,ibox) = j_atom_list(ja,jm)%rp(2)
                                                live_atom_rsp(ja,jm,3,js,ibox) = j_atom_list(ja,jm)%rp(3)
                                                live_atom_exist(ja,jm,js,ibox) = j_atom_list(ja,jm)%exist
                                                !live_atom_list(ja,jm,js,ibox)%rp(1) = j_atom_list(ja,jm)%rp(1)
                                                !live_atom_list(ja,jm,js,ibox)%rp(2) = j_atom_list(ja,jm)%rp(2)
                                                !live_atom_list(ja,jm,js,ibox)%rp(3) = j_atom_list(ja,jm)%rp(3)
                                                !live_atom_list(ja,jm,js,ibox)%exist = j_atom_list(ja,jm)%exist
                                        END DO
                                END DO
                                !$OMP END DO SIMD
                                !$OMP DO SIMD
                                DO jm = 1, jnlive
                                        live_xcom(jm,js,ibox) = j_molecule_list(jm)%rcom(1)
                                        live_ycom(jm,js,ibox) = j_molecule_list(jm)%rcom(2)
                                        live_zcom(jm,js,ibox) = j_molecule_list(jm)%rcom(3)
                                        live_max_dcom(jm,js,ibox) = j_molecule_list(jm)%rcom(4)
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
    REAL(DP) :: rterm, rterm2, eps, sigsq, negsigsq, negsigbyr2, this_rijsq, epsig_n, epsig_m
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
    max_dcom_i_const = widom_molecule%rcom(4) + mol_rcut
    n_i_exist = COUNT(widom_atoms%exist)
    which_i_exist = PACK(vec123(1:natoms(is)),widom_atoms%exist)
    irp(1,:) = widom_atoms(which_i_exist)%rp(1)
    irp(2,:) = widom_atoms(which_i_exist)%rp(2)
    irp(3,:) = widom_atoms(which_i_exist)%rp(3)
    ispec_atomtypes = nonbond_list(which_i_exist,is)%atom_type_number
    icharge = nonbond_list(which_i_exist,is)%charge * charge_factor
    IF (l_ortho) THEN
            ixcom = widom_molecule%rcom(1)
            iycom = widom_molecule%rcom(2)
            izcom = widom_molecule%rcom(3)
    ELSE
            irp_com(1) = widom_molecule%rcom(1)
            irp_com(2) = widom_molecule%rcom(2)
            irp_com(3) = widom_molecule%rcom(3)
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
                        !!$OMP CRITICAL
                        !WRITE(*,*) "CRITICAL OVERLAP!", MINVAL(rijsq(1:vlen)), MAXVAL(rijsq(1:vlen)), ia, is, cbmc_flag
                        !WRITE(*,*) "check_overlap: ", check_overlap(ia,im,is)
                        !!$OMP END CRITICAL
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
                        ij_vdw_p_table(1:2,1:n_vdw) = ppvdwp_table2(1:2,PACK(jatomtype(1:vlen),vdw_mask(1:vlen)),itype,ibox)
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
                        packed_types(1:n_vdw) = PACK(jatomtype(1:vlen),vdw_mask(1:vlen))
                        mie_vdw_p_table(1:4,1:n_vdw) = ppvdwp_table2(1:4,packed_types(1:n_vdw),itype,this_box)
                        DO i = 1, n_vdw
                                this_rijsq = LOG(rijsq_packed(i)) ! actually ln(rijsq)
                                epsig_n = mie_vdw_p_table(1,i)  ! epsilon * mie_coeff * sigma ** n
                                epsig_m = mie_vdw_p_table(2,i) ! epsilon * mie_coeff * sigma ** m
                                mie_n = mie_vdw_p_table(3,i) ! already halved
                                mie_m = mie_vdw_p_table(4,i) ! already halved
                                nrg_vdw = epsig_n * EXP(this_rijsq*mie_n)
                                nrg_vdw = nrg_vdw - epsig_m * EXP(this_rijsq*mie_m)
                                IF (int_vdw_sum_style(this_box) == vdw_cut_shift) THEN
                                        nrg_vdw = nrg_vdw - ppvdwp_table(packed_types(i),itype,5,this_box)
                                        !nrg_vdw = nrg_vdw + epsig_n * rcut_vdwsq(this_box)**mie_n
                                        !nrg_vdw = nrg_vdw - epsig_m * rcut_vdwsq(this_box)**mie_m
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

  SUBROUTINE Compute_Molecule_Nonbond_Inter_Energy_Vectorized_2(im,is, &
    this_box,vdw_energy,qq_energy,overlap)
    ! Arguments
    INTEGER, INTENT(IN) :: im, is, this_box
    REAL(DP), INTENT(OUT) :: vdw_energy, qq_energy
    LOGICAL, INTENT(OUT) :: overlap
    ! End Arguments
    INTEGER, DIMENSION(nspecies) :: jspec_n_interact, jspec_istart, jspec_iend
    INTEGER, DIMENSION(COUNT(atom_list(1:natoms(is),im,is)%exist)) :: which_i_exist, ispec_atomtypes
    REAL(DP), DIMENSION(COUNT(atom_list(1:natoms(is),im,is)%exist)) :: icharge
    REAL(DP), DIMENSION(3,COUNT(atom_list(1:natoms(is),im,is)%exist)) :: irp, isp
    REAL(DP), DIMENSION(3) :: irp_com, isp_com
    LOGICAL(4), DIMENSION(COUNT(atom_list(1:natoms(is),im,is)%exist)) :: i_has_vdw
    !LOGICAL, DIMENSION(MAXVAL(nmols(:,this_box))), TARGET :: spec_live_tgt
    REAL(DP), DIMENSION(maxboxnatoms) :: rijsq_4min
    LOGICAL(4), DIMENSION(MAXVAL(nmols(:,this_box))) :: interact_vec
    INTEGER, DIMENSION(MAXVAL(nmols(:,this_box))) :: which_interact, live_locates
    REAL(DP), DIMENSION(4,MAXVAL(nmols(:,this_box))) :: live_rcom
    LOGICAL(1), DIMENSION(maxboxnatoms) :: vdw_mask, coul_mask, j_hascharge, j_hasvdw, j_exist
    INTEGER, DIMENSION(maxboxnatoms) :: jatomtype, packed_types
    REAL(DP), DIMENSION(maxboxnatoms) :: jcharge, up_nrg_vec, rijsq_packed, rijsq, jcharge_coul, jxp, jyp, jzp
    REAL(DP), DIMENSION(maxboxnatoms,3) :: jrsp
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
    REAL(DP) :: dsxcom, dsycom, dszcom, jrp, rxp, ryp , rzp, jsxp, jsyp, jszp
    REAL(DP) :: xl, yl, zl, hxl, hyl, hzl
    REAL(DP) :: h11,h21,h31,h12,h22,h32,h13,h23,h33
    REAL(DP) :: inv_h11,inv_h21,inv_h31,inv_h12,inv_h22,inv_h32,inv_h13,inv_h23,inv_h33

    TYPE(Molecule_Class) :: molecule_i
    TYPE(Atom_Class), DIMENSION(natoms(is)) :: atoms_i
    LOGICAL :: l_all_exist, jex_redux, jex
    INTEGER :: j_limit, ibase,ibase2,jloc,k
    REAL(DP) :: max_dcom
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
    molecule_i = molecule_list(im,is)
    atoms_i = atom_list(1:natoms(is),im,is)
    ibox = molecule_i%which_box


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
            inv_h11 = box_list(this_box)%length_inv(1,1)
            inv_h21 = box_list(this_box)%length_inv(2,1)
            inv_h31 = box_list(this_box)%length_inv(3,1)
            inv_h12 = box_list(this_box)%length_inv(1,2)
            inv_h22 = box_list(this_box)%length_inv(2,2)
            inv_h32 = box_list(this_box)%length_inv(3,2)
            inv_h13 = box_list(this_box)%length_inv(1,3)
            inv_h23 = box_list(this_box)%length_inv(2,3)
            inv_h33 = box_list(this_box)%length_inv(3,3)
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
    max_dcom_i_const = molecule_i%rcom(4) + mol_rcut
    n_i_exist = COUNT(atoms_i%exist)
    which_i_exist = PACK(vec123(1:natoms(is)),atoms_i%exist)
    irp(1,:) = atoms_i(which_i_exist)%rp(1)
    irp(2,:) = atoms_i(which_i_exist)%rp(2)
    irp(3,:) = atoms_i(which_i_exist)%rp(3)
    ispec_atomtypes = nonbond_list(which_i_exist,is)%atom_type_number
    icharge = nonbond_list(which_i_exist,is)%charge * charge_factor
    ixcom = molecule_i%rcom(1)
    iycom = molecule_i%rcom(2)
    izcom = molecule_i%rcom(3)
    molecule_list(im,is)%live = .FALSE.
    l_all_exist = .TRUE.
    istart = 1
    DO js = 1, nspecies
        jnmols = nmols(js,ibox)
        IF (jnmols == 0) CYCLE
        jnatoms = natoms(js)
        n_interact = 0
        IF (open_mc_flag) THEN
                live_locates(1:jnmols) = locate(1:jnmols,js,ibox)
                IF (l_not_all_live) THEN
                        jnlive = COUNT(molecule_list(live_locates(1:jnmols),js)%live)
                ELSE IF (is .EQ. js) THEN
                        jnlive = jnmols - 1
                ELSE
                        jnlive = jnmols
                END IF
                IF (jnlive == 0) CYCLE
                IF (is .EQ. js .OR. l_not_all_live) THEN
                        live_locates(1:jnlive) = &
                                PACK(live_locates(1:jnmols),molecule_list(live_locates(1:jnmols),js)%live)
                END IF
                DO i = 1, jnlive
                        live_rcom(:,i) = molecule_list(live_locates(i),js)%rcom
                END DO
        ELSE
                jnlive = jnmols
        END IF
        IF (open_mc_flag) THEN
                j_limit = jnlive
        ELSE IF (is .EQ. js) THEN
                j_limit = im - 1
                interact_vec(im) = .FALSE.
        ELSE
                j_limit = jnlive
        END IF
        IF (l_ortho) THEN
                !$OMP PARALLEL
                !$OMP DO SIMD SCHEDULE(STATIC) &
                !$OMP PRIVATE(dxcom,dycom,dzcom,max_dcom,l_interact) &
                !$OMP REDUCTION(+: n_interact)
                DO j = 1, j_limit
                        IF (open_mc_flag) THEN
                                dxcom = live_rcom(1,j)
                                dycom = live_rcom(2,j)
                                dzcom = live_rcom(3,j)
                                max_dcom = live_rcom(4,j)
                        ELSE
                                dxcom = molecule_list(j,js)%rcom(1)
                                dycom = molecule_list(j,js)%rcom(2)
                                dzcom = molecule_list(j,js)%rcom(3)
                                max_dcom = molecule_list(j,js)%rcom(4)
                        END IF
                        dxcom = ABS(dxcom - ixcom)
                        dycom = ABS(dycom - iycom)
                        dzcom = ABS(dzcom - izcom)
                        IF (dxcom > hxl) dxcom = dxcom - xl
                        IF (dycom > hyl) dycom = dycom - yl
                        IF (dzcom > hzl) dzcom = dzcom - zl
                        ! Repurposing dxcom as rijsq accumulator to enforce FMA3 instructions if supported
                        dxcom = dxcom * dxcom
                        dxcom = dxcom + dycom * dycom
                        dxcom = dxcom + dzcom * dzcom
                        l_interact = max_dcom_i_const > &
                                SQRT(dxcom) - max_dcom
                        interact_vec(j) = l_interact
                        IF (l_interact) n_interact = n_interact + 1
                END DO
                !$OMP END DO SIMD
                !$OMP DO SIMD SCHEDULE(STATIC) &
                !$OMP PRIVATE(dxcom,dycom,dzcom,max_dcom,l_interact) &
                !$OMP REDUCTION(+: n_interact)
                DO j = j_limit+2, jnlive
                        IF (open_mc_flag) THEN
                                dxcom = live_rcom(1,j)
                                dycom = live_rcom(2,j)
                                dzcom = live_rcom(3,j)
                                max_dcom = live_rcom(4,j)
                        ELSE
                                dxcom = molecule_list(j,js)%rcom(1)
                                dycom = molecule_list(j,js)%rcom(2)
                                dzcom = molecule_list(j,js)%rcom(3)
                                max_dcom = molecule_list(j,js)%rcom(4)
                        END IF
                        dxcom = ABS(dxcom - ixcom)
                        dycom = ABS(dycom - iycom)
                        dzcom = ABS(dzcom - izcom)
                        IF (dxcom > hxl) dxcom = dxcom - xl
                        IF (dycom > hyl) dycom = dycom - yl
                        IF (dzcom > hzl) dzcom = dzcom - zl
                        ! Repurposing dxcom as rijsq accumulator to enforce FMA3 instructions if supported
                        dxcom = dxcom * dxcom
                        dxcom = dxcom + dycom * dycom
                        dxcom = dxcom + dzcom * dzcom
                        l_interact = max_dcom_i_const > &
                                SQRT(dxcom) - max_dcom
                        interact_vec(j) = l_interact
                        IF (l_interact) n_interact = n_interact + 1
                END DO
                !$OMP END DO SIMD
                !$OMP END PARALLEL
        ELSE
                !$OMP PARALLEL
                !$OMP DO SIMD SCHEDULE(STATIC) &
                !$OMP PRIVATE(dxcom,dycom,dzcom,max_dcom,l_interact) &
                !$OMP PRIVATE(dsxcom,dsycom,dszcom) &
                !$OMP REDUCTION(+: n_interact)
                DO i = 1, j_limit
                        IF (open_mc_flag) THEN
                                dxcom = live_rcom(1,j)
                                dycom = live_rcom(2,j)
                                dzcom = live_rcom(3,j)
                                max_dcom = live_rcom(4,j)
                        ELSE
                                dxcom = molecule_list(j,js)%rcom(1)
                                dycom = molecule_list(j,js)%rcom(2)
                                dzcom = molecule_list(j,js)%rcom(3)
                                max_dcom = molecule_list(j,js)%rcom(4)
                        END IF
                        dxcom = dxcom - ixcom
                        dsxcom = inv_h11*dxcom
                        dsycom = inv_h21*dxcom
                        dszcom = inv_h31*dxcom
                        dycom = dycom - iycom
                        dsxcom = dsxcom + inv_h12*dycom
                        dsycom = dsycom + inv_h22*dycom
                        dszcom = dszcom + inv_h32*dycom
                        dzcom = dzcom - izcom
                        dsxcom = dsxcom + inv_h13*dzcom
                        dsycom = dsycom + inv_h23*dzcom
                        dszcom = dszcom + inv_h33*dzcom
                        IF (dsxcom > 0.5_DP) THEN
                                dsxcom = dsxcom - 1.0_DP
                        ELSE IF (dsxcom < -0.5_DP) THEN
                                dsxcom = dsxcom + 1.0_DP
                        END IF
                        IF (dsycom > 0.5_DP) THEN
                                dsycom = dsycom - 1.0_DP
                        ELSE IF (dsycom < -0.5_DP) THEN
                                dsycom = dsycom + 1.0_DP
                        END IF
                        IF (dszcom > 0.5_DP) THEN
                                dszcom = dszcom - 1.0_DP
                        ELSE IF (dszcom < -0.5_DP) THEN
                                dszcom = dszcom + 1.0_DP
                        END IF
                        dxcom = h11*dsxcom
                        dycom = h21*dsxcom
                        dzcom = h31*dsxcom
                        dxcom = dxcom + h12*dsycom
                        dycom = dycom + h22*dsycom
                        dzcom = dzcom + h32*dsycom
                        dxcom = dxcom + h13*dszcom
                        dycom = dycom + h23*dszcom
                        dzcom = dzcom + h33*dszcom
                        ! Repurposing dxcom as rijsq accumulator to enforce FMA3 instructions if supported
                        dxcom = dxcom * dxcom
                        dxcom = dxcom + dycom * dycom
                        dxcom = dxcom + dzcom * dzcom
                        l_interact = max_dcom_i_const > &
                                SQRT(dxcom) - max_dcom
                        interact_vec(i) = l_interact
                        IF (l_interact) n_interact = n_interact + 1
                END DO
                !$OMP END DO SIMD
                !$OMP DO SIMD SCHEDULE(STATIC) &
                !$OMP PRIVATE(dxcom,dycom,dzcom,max_dcom,l_interact) &
                !$OMP PRIVATE(dsxcom,dsycom,dszcom) &
                !$OMP REDUCTION(+: n_interact)
                DO i = j_limit+2, jnlive
                        IF (open_mc_flag) THEN
                                dxcom = live_rcom(1,j)
                                dycom = live_rcom(2,j)
                                dzcom = live_rcom(3,j)
                                max_dcom = live_rcom(4,j)
                        ELSE
                                dxcom = molecule_list(j,js)%rcom(1)
                                dycom = molecule_list(j,js)%rcom(2)
                                dzcom = molecule_list(j,js)%rcom(3)
                                max_dcom = molecule_list(j,js)%rcom(4)
                        END IF
                        dxcom = dxcom - ixcom
                        dsxcom = inv_h11*dxcom
                        dsycom = inv_h21*dxcom
                        dszcom = inv_h31*dxcom
                        dycom = dycom - iycom
                        dsxcom = dsxcom + inv_h12*dycom
                        dsycom = dsycom + inv_h22*dycom
                        dszcom = dszcom + inv_h32*dycom
                        dzcom = dzcom - izcom
                        dsxcom = dsxcom + inv_h13*dzcom
                        dsycom = dsycom + inv_h23*dzcom
                        dszcom = dszcom + inv_h33*dzcom
                        IF (dsxcom > 0.5_DP) THEN
                                dsxcom = dsxcom - 1.0_DP
                        ELSE IF (dsxcom < -0.5_DP) THEN
                                dsxcom = dsxcom + 1.0_DP
                        END IF
                        IF (dsycom > 0.5_DP) THEN
                                dsycom = dsycom - 1.0_DP
                        ELSE IF (dsycom < -0.5_DP) THEN
                                dsycom = dsycom + 1.0_DP
                        END IF
                        IF (dszcom > 0.5_DP) THEN
                                dszcom = dszcom - 1.0_DP
                        ELSE IF (dszcom < -0.5_DP) THEN
                                dszcom = dszcom + 1.0_DP
                        END IF
                        dxcom = h11*dsxcom
                        dycom = h21*dsxcom
                        dzcom = h31*dsxcom
                        dxcom = dxcom + h12*dsycom
                        dycom = dycom + h22*dsycom
                        dzcom = dzcom + h32*dsycom
                        dxcom = dxcom + h13*dszcom
                        dycom = dycom + h23*dszcom
                        dzcom = dzcom + h33*dszcom
                        ! Repurposing dxcom as rijsq accumulator to enforce FMA3 instructions if supported
                        dxcom = dxcom * dxcom
                        dxcom = dxcom + dycom * dycom
                        dxcom = dxcom + dzcom * dzcom
                        l_interact = max_dcom_i_const > &
                                SQRT(dxcom) - max_dcom
                        interact_vec(i) = l_interact
                        IF (l_interact) n_interact = n_interact + 1
                END DO
                !$OMP END DO SIMD
                !$OMP END PARALLEL
        END IF
        IF (open_mc_flag) THEN
                !$OMP PARALLEL WORKSHARE
                which_interact(1:n_interact) = PACK(live_locates(1:jnlive),interact_vec(1:jnlive))
                !$OMP END PARALLEL WORKSHARE
        ELSE
                !$OMP PARALLEL WORKSHARE
                which_interact(1:n_interact) = PACK(vec123(1:jnlive),interact_vec(1:jnlive))
                !$OMP END PARALLEL WORKSHARE
        END IF
        vlen = n_interact*jnatoms
        iend = istart + vlen - 1
        ibase = istart - 1
        !$OMP PARALLEL PRIVATE(ibase2,jloc,i,rxp,ryp,rzp,jex) &
        !$OMP REDUCTION(.AND.: jex_redux)
        jex_redux = .TRUE.
        !$OMP DO SCHEDULE(STATIC)
        DO k = 1, n_interact
                ibase2 = ibase + (k-1)*jnatoms
                jloc = which_interact(k)
                DO i = 1, jnatoms
                        rxp = atom_list(i,jloc,js)%rp(1)
                        ryp = atom_list(i,jloc,js)%rp(2)
                        rzp = atom_list(i,jloc,js)%rp(3)
                        jex = atom_list(i,jloc,js)%exist
                        jrsp(ibase2+i,1) = rxp
                        jrsp(ibase2+i,2) = ryp
                        jrsp(ibase2+i,3) = rzp
                        j_exist(ibase2+i) = jex
                        jex_redux = jex_redux .AND. jex
                END DO
        END DO
        !$OMP END DO
        !$OMP WORKSHARE
        jatomtype(istart:iend) = RESHAPE(SPREAD(nonbond_list(1:natoms(js),js)%atom_type_number,2,n_interact), (/ vlen /))
        jcharge(istart:iend) = RESHAPE(SPREAD(nonbond_list(1:natoms(js),js)%charge,2,n_interact), (/ vlen /))
        !$OMP END WORKSHARE
        !$OMP END PARALLEL
        l_all_exist = l_all_exist .AND. jex_redux
        IF (l_get_rij_min) THEN
                jspec_n_interact(js) = n_interact
                jspec_istart(js) = istart
                jspec_iend(js) = iend
        END IF
        istart = istart+vlen
    END DO
    molecule_list(im,is)%live = .TRUE.
    vlen = istart - 1
    orig_vlen = vlen
    ! The following IF statement condition is unlikely to be true.
    ! That is why this process is saved until now; it's unlikely to be necessary.
    IF (.NOT. l_all_exist) THEN
            !$OMP PARALLEL WORKSHARE
            n_j_exist = COUNT(j_exist(1:vlen))
            jcharge(1:n_j_exist) = PACK(jcharge(1:vlen),j_exist(1:vlen))
            jatomtype(1:n_j_exist) = PACK(jatomtype(1:vlen),j_exist(1:vlen))
            jrsp(1:n_j_exist,1) = PACK(jrsp(1:vlen,1),j_exist(1:vlen))
            jrsp(1:n_j_exist,2) = PACK(jrsp(1:vlen,2),j_exist(1:vlen))
            jrsp(1:n_j_exist,3) = PACK(jrsp(1:vlen,3),j_exist(1:vlen))
            !$OMP END PARALLEL WORKSHARE
            vlen = n_j_exist
    END IF
    IF (this_est_emax) THEN
            up_nrg_vec = 0.0_DP
    END IF
    !$OMP PARALLEL WORKSHARE
    j_hascharge(1:vlen) = jcharge(1:vlen) .NE. 0.0_DP
    j_hasvdw(1:vlen) = jatomtype(1:vlen) .NE. 0
    l_notallvdw = .NOT. ALL(j_hasvdw(1:vlen))
    l_notallcoul = .NOT. ALL(j_hascharge(1:vlen))
    !$OMP END PARALLEL WORKSHARE
    l_notsamecut = this_vdw_rcutsq .NE. this_coul_rcutsq
    l_pack_separately = l_notsamecut .OR. l_notallcoul .OR. l_notallvdw
    IF (.NOT. l_ortho) THEN
            !$OMP PARALLEL
            !$OMP DO SIMD SCHEDULE(STATIC) &
            !$OMP PRIVATE(jrp,jsxp,jsyp,jszp)
            DO i = 1, vlen
                    jrp = jrsp(i,1)
                    jsxp = inv_h11*jrp
                    jsyp = inv_h21*jrp
                    jszp = inv_h31*jrp
                    jrp = jrsp(i,2)
                    jsxp = jsxp + inv_h12*jrp
                    jsyp = jsyp + inv_h22*jrp
                    jszp = jszp + inv_h32*jrp
                    jrp = jrsp(i,3)
                    jsxp = jsxp + inv_h13*jrp
                    jsyp = jsyp + inv_h23*jrp
                    jszp = jszp + inv_h33*jrp
                    jrsp(i,1) = jsxp
                    jrsp(i,2) = jsyp
                    jrsp(i,3) = jszp
            END DO
            !$OMP END DO SIMD
            !$OMP END PARALLEL
    END IF

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
                !$OMP PARALLEL
                !$OMP DO SIMD SCHEDULE(STATIC) &
                !$OMP PRIVATE(dxp,dyp,dzp,l_interact)
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
                !$OMP END DO SIMD
                !$OMP END PARALLEL
        ELSE
                ixp = isp(1,ia)
                iyp = isp(2,ia)
                izp = isp(3,ia)
                !$OMP PARALLEL
                !$OMP DO SIMD SCHEDULE(STATIC) &
                !$OMP PRIVATE(dsp,dxp,dyp,dzp,l_interact)
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
                !$OMP END DO SIMD
                !$OMP END PARALLEL
        END IF
        !$OMP PARALLEL WORKSHARE
        overlap = ANY(rijsq(1:vlen) < rcut_lowsq)
        !$OMP END PARALLEL WORKSHARE
        IF (overlap) RETURN
        IF (i_get_vdw) THEN
                IF (i_get_coul .AND. .NOT. l_notsamecut) THEN
                        IF (l_notallcoul) THEN
                                !$OMP PARALLEL WORKSHARE
                                coul_mask(1:vlen) = vdw_mask(1:vlen) .AND. j_hascharge(1:vlen)
                                !$OMP END PARALLEL WORKSHARE
                        ELSE IF (l_notallvdw) THEN
                                !$OMP PARALLEL WORKSHARE
                                coul_mask(1:vlen) = vdw_mask(1:vlen)
                                !$OMP END PARALLEL WORKSHARE
                        ELSE
                                !$OMP PARALLEL WORKSHARE
                                n_coul = COUNT(vdw_mask(1:vlen))
                                jcharge_coul(1:n_coul) = PACK(jcharge(1:vlen),vdw_mask(1:vlen))
                                !$OMP END PARALLEL WORKSHARE
                        END IF
                END IF
                IF (l_notallvdw) THEN
                        !$OMP PARALLEL WORKSHARE
                        vdw_mask(1:vlen) = vdw_mask(1:vlen) .AND. j_hasvdw(1:vlen)
                        !$OMP END PARALLEL WORKSHARE
                END IF
                !$OMP PARALLEL WORKSHARE
                n_vdw = COUNT(vdw_mask(1:vlen))
                rijsq_packed(1:n_vdw) = PACK(rijsq(1:vlen),vdw_mask(1:vlen))
                !$OMP END PARALLEL WORKSHARE
                this_int_vdw_style = int_vdw_style(this_box)
                IF (this_int_vdw_style == vdw_lj) THEN
                        !$OMP PARALLEL WORKSHARE
                        ij_vdw_p_table(1:2,1:n_vdw) = ppvdwp_table2(1:2,PACK(jatomtype(1:vlen),vdw_mask(1:vlen)),itype,ibox)
                        !$OMP END PARALLEL WORKSHARE
                        IF (int_vdw_sum_style(this_box) == vdw_charmm) THEN
                                !$OMP PARALLEL
                                !$OMP DO SIMD SCHEDULE(STATIC) &
                                !$OMP PRIVATE(eps,sigsq,sigbyr2,sigbyr12) &
                                !$OMP REDUCTION(+: vdw_energy)
                                DO i = 1, n_vdw
                                        eps = ij_vdw_p_table(1,i) ! epsilon
                                        sigsq = ij_vdw_p_table(2,i) ! sigma**2
                                        sigbyr2 = sigsq/rijsq_packed(i) ! sigma was already squared
                                        sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                        sigbyr12 = sigbyr6*sigbyr6
                                        !nrg_vdw = eps * (sigbyr12 - 2.0_DP*sigbyr6) ! eps was not multiplied by 4
                                        !vdw_energy = vdw_energy + nrg_vdw
                                        sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
                                        vdw_energy = vdw_energy + eps * sigbyr12
                                END DO
                                !$OMP END DO SIMD
                                !$OMP END PARALLEL
                        ELSE IF (int_vdw_sum_style(this_box) == vdw_cut_shift) THEN
                                !$OMP PARALLEL
                                !$OMP DO SIMD SCHEDULE(STATIC) &
                                !$OMP PRIVATE(eps,negsigsq,negsigbyr2,rterm,nrg_vdw) &
                                !$OMP REDUCTION(+: vdw_energy)
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
                                END DO
                                !$OMP END DO SIMD
                                !$OMP END PARALLEL
                        ELSE IF (int_vdw_sum_style(this_box) == vdw_cut_switch) THEN
                                !$OMP PARALLEL
                                !$OMP DO SIMD SCHEDULE(STATIC) &
                                !$OMP PRIVATE(eps,negsigsq,this_rijsq,negsigbyr2,rterm,nrg_vdw,roffsq_rijsq) &
                                !$OMP REDUCTION(+: vdw_energy)
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
                                END DO
                                !$OMP END DO SIMD
                                !$OMP END PARALLEL
                        ELSE IF (int_vdw_sum_style(this_box) == vdw_cut_shift_force) THEN
                                !$OMP PARALLEL
                                !$OMP DO SIMD SCHEDULE(STATIC) &
                                !$OMP PRIVATE(eps,negsigsq,this_rijsq,negsigbyr2,rterm,nrg_vdw,rterm2) &
                                !$OMP REDUCTION(+: vdw_energy)
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
                                END DO
                                !$OMP END DO SIMD
                                !$OMP END PARALLEL
                        ELSE
                                !$OMP PARALLEL
                                !$OMP DO SIMD SCHEDULE(STATIC) &
                                !$OMP PRIVATE(eps,negsigsq,negsigbyr2,rterm) &
                                !$OMP REDUCTION(+: vdw_energy)
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
                                END DO
                                !$OMP END DO SIMD
                                !$OMP END PARALLEL
                        END IF
                ELSE IF (int_vdw_style(this_box) == vdw_mie) THEN
                        !DIR$ ASSUME (n_vdw_p_list(this_box) .EQ. 2)
                        !$OMP PARALLEL
                        !$OMP WORKSHARE
                        mie_vdw_p_table(1:4,1:n_vdw) = ppvdwp_table2(1:4,PACK(jatomtype(1:vlen),vdw_mask(1:vlen)),itype,this_box)
                        !$OMP END WORKSHARE
                        !$OMP DO SIMD SCHEDULE(STATIC) &
                        !$OMP PRIVATE(eps,sigsq,mie_n,mie_m,sigbyr2,nrg_vdw) &
                        !$OMP REDUCTION(+: vdw_energy)
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
                        END DO
                        !$OMP END DO SIMD
                        !$OMP END PARALLEL
                END IF
        END IF
        ! Electrostatics
        IF (i_get_coul) THEN
                icharge_factor = icharge(ia) ! already multiplied by charge_factor
                IF (l_pack_separately .OR. .NOT. i_get_vdw) THEN
                        IF (l_notallcoul .AND. l_check_coul) THEN
                                !$OMP PARALLEL WORKSHARE
                                coul_mask(1:vlen) = coul_mask(1:vlen) .AND. j_hascharge(1:vlen)
                                !$OMP END PARALLEL WORKSHARE
                        END IF
                        !$OMP PARALLEL WORKSHARE
                        n_coul = COUNT(coul_mask(1:vlen))
                        jcharge_coul(1:n_coul) = PACK(jcharge(1:vlen), coul_mask(1:vlen))
                        rijsq_packed(1:n_coul) = PACK(rijsq(1:vlen), coul_mask(1:vlen))
                        !$OMP END PARALLEL WORKSHARE
                ELSE
                        n_coul = n_vdw
                END IF
                i_qq_energy = 0.0_DP
                !$OMP PARALLEL
                !$OMP DO SIMD SCHEDULE(STATIC) &
                !$OMP PRIVATE(rij,nrg_coul) &
                !$OMP REDUCTION(+: i_qq_energy)
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
                !$OMP END DO SIMD
                !$OMP END PARALLEL
                qq_energy = qq_energy + i_qq_energy * icharge_factor
        END IF
    END DO
  END SUBROUTINE Compute_Molecule_Nonbond_Inter_Energy_Vectorized_2

  SUBROUTINE Ortho_Minimum_Image(drxp,dryp,drzp,xl,yl,zl,hxl,hyl,hzl)
          !!!$OMP DECLARE SIMD(Ortho_Minimum_Image) UNIFORM(xl,yl,zl,hxl,hyl,hzl)
          REAL(DP) :: drxp,dryp,drzp,xl,yl,zl,hxl,hyl,hzl

          drxp = ABS(drxp)
          dryp = ABS(dryp)
          drzp = ABS(drzp)
          IF (drxp > hxl) drxp = drxp - xl
          IF (dryp > hyl) dryp = dryp - yl
          IF (drzp > hzl) drzp = drzp - zl
  END SUBROUTINE Ortho_Minimum_Image


  SUBROUTINE Triclinic_Minimum_Image(drxp,dryp,drzp,inv_h,h)
          !!!$OMP DECLARE SIMD(Triclinic_Minimum_Image) UNIFORM(inv_h,h)
          REAL(DP) :: drxp,dryp,drzp,dsxp,dsyp,dszp,inv_h(3,3),h(3,3)

          dsxp = inv_h(1,1)*drxp
          dsyp = inv_h(2,1)*drxp
          dszp = inv_h(3,1)*drxp
          dsxp = dsxp + inv_h(1,2)*dryp
          dsyp = dsyp + inv_h(2,2)*dryp
          dszp = dszp + inv_h(3,2)*dryp
          dsxp = dsxp + inv_h(1,3)*drzp
          dsyp = dsyp + inv_h(2,3)*drzp
          dszp = dszp + inv_h(3,3)*drzp
          IF (dsxp > 0.5_DP) THEN
                  dsxp = dsxp - 1.0_DP
          ELSE IF (dsxp < -0.5_DP) THEN
                  dsxp = dsxp + 1.0_DP
          END IF
          IF (dsyp > 0.5_DP) THEN
                  dsyp = dsyp - 1.0_DP
          ELSE IF (dsyp < -0.5_DP) THEN
                  dsyp = dsyp + 1.0_DP
          END IF
          IF (dszp > 0.5_DP) THEN
                  dszp = dszp - 1.0_DP
          ELSE IF (dszp < -0.5_DP) THEN
                  dszp = dszp + 1.0_DP
          END IF
          drxp = h(1,1)*dsxp
          dryp = h(2,1)*dsxp
          drzp = h(3,1)*dsxp
          drxp = drxp + h(1,2)*dsyp
          dryp = dryp + h(2,2)*dsyp
          drzp = drzp + h(3,2)*dsyp
          drxp = drxp + h(1,3)*dszp
          dryp = dryp + h(2,3)*dszp
          drzp = drzp + h(3,3)*dszp
  END SUBROUTINE Triclinic_Minimum_Image

  SUBROUTINE Set_Pair_Chunks
          INTEGER :: chunksize, nthreads_used, jnlive, ibox, js
          IF (ALLOCATED(chunksize_array)) THEN
                  IF (ALL(chunks_set_nmols .EQ. nmols)) RETURN
          ELSE
                  ALLOCATE(chunks_set_nmols(nspecies,0:nbr_boxes))
                  ! js, is
                  ALLOCATE(chunksize_array(nspecies,nspecies,nbr_boxes))
                  ALLOCATE(nthreads_used_array(nspecies,nspecies,nbr_boxes))
          END IF
          chunks_set_nmols = nmols
          DO ibox = 1, nbr_boxes
                DO js = 1, nspecies
                        jnlive = nmols(js,ibox)
                        chunksize = MAX(jnlive / global_nthreads + 1,8)
                        IF (MOD(chunksize,4) .NE. 0) chunksize = 4*(chunksize/4+1)
                        nthreads_used = jnlive/chunksize
                        IF (MOD(jnlive,chunksize) .NE. 0) nthreads_used = nthreads_used + 1
                        chunksize_array(js,:,ibox) = chunksize
                        nthreads_used_array(js,:,ibox) = nthreads_used
                        jnlive = jnlive - 1
                        chunksize = MAX(jnlive / global_nthreads + 1,8)
                        IF (MOD(chunksize,4) .NE. 0) chunksize = 4*(chunksize/4+1)
                        nthreads_used = jnlive/chunksize
                        IF (MOD(jnlive,chunksize) .NE. 0) nthreads_used = nthreads_used + 1
                        chunksize_array(js,js,ibox) = chunksize
                        nthreads_used_array(js,js,ibox) = nthreads_used
                END DO
          END DO
  END SUBROUTINE Set_Pair_Chunks

  SUBROUTINE Compute_Molecule_Nonbond_Inter_Energy_Vectorized(im,is, &
    i_vdw_energy,i_qq_energy,overlap)
    ! Arguments
    INTEGER, INTENT(IN) :: im, is
    REAL(DP), INTENT(OUT) :: i_vdw_energy, i_qq_energy
    LOGICAL, INTENT(OUT) :: overlap
    ! End Arguments
    INTEGER, DIMENSION(COUNT(atom_list(1:natoms(is),im,is)%exist)) :: which_i_exist, iatomtypes
    REAL(DP), DIMENSION(COUNT(atom_list(1:natoms(is),im,is)%exist)) :: irsxp_vec, irsyp_vec, irszp_vec
    REAL(DP), DIMENSION(3) :: irp_com, isp_com
    !LOGICAL, DIMENSION(MAXVAL(nmols(:,this_box))), TARGET :: spec_live_tgt
    !INTEGER, DIMENSION(MAXVAL(nlive(:,this_box))) :: which_interact
    LOGICAL :: this_est_emax, l_get_rij_min, l_ortho
    INTEGER :: ia_counter, i, j, js, n_i_exist, istart, iend, jnlive, orig_vlen, n_j_exist
    INTEGER :: bsolvent, istart_base, natoms_js, ti_solvent, ja, n_coul, n_vdw, ia, ia0, live_vlen, jnmols, jnatoms
    REAL(DP) :: mol_rcut, max_dcom_i_const, this_vdw_rcutsq, this_coul_rcutsq, this_coul_rcut

    INTEGER :: this_int_vdw_style, ibox, this_box

    LOGICAL :: get_vdw, get_qq, i_get_qq
    LOGICAL(1) :: l_interact
    REAL(DP) :: dxcom, dycom, dzcom, dscom, ixcom, iycom, izcom
    REAL(DP) :: dsxcom, dsycom, dszcom, jrp
    REAL(DP) :: xl, yl, zl, hxl, hyl, hzl
    REAL(DP) :: h11,h21,h31,h12,h22,h32,h13,h23,h33
    REAL(DP) :: inv_h11,inv_h21,inv_h31,inv_h12,inv_h22,inv_h32,inv_h13,inv_h23,inv_h33

    INTEGER :: j_limit, j_limit_2, ibase,ibase2,jloc,k
    REAL(DP) :: max_dcom

    LOGICAL(1), DIMENSION(MAXVAL(max_molecules)) :: interact_vec, all_interact_vec
    INTEGER, DIMENSION(MAXVAL(max_molecules)) :: live_locates, which_may_interact_p, &
            which_all_interact_p, which_interact
    INTEGER :: n_all_interact, n_may_interact, n_interact, n_all_interact_p, n_may_interact_p
    REAL(DP), DIMENSION(4,MAXVAL(max_molecules)) :: live_rcom

    INTEGER, DIMENSION(MAXVAL(natoms)) :: jatomtypes

    LOGICAL :: lj_charmm, lj_cut_shift, lj_cut_switch, lj_cut_shift_force, lj_cut, mie_cut_shift, mie_cut, l_charge_ewald, l_charge_dsf
    LOGICAL :: l_charge_cut, need_sqrt, l_pair_store, ij_overlap, i_overlap, l_order_ij

    INTEGER :: jm, jsl, isl, jsl_base, j_interact
    INTEGER(2) :: natompairs, ji, n_vdw_p
    REAL(DP) :: dxp, dyp, dzp, rijsq, rij, vdw_energy, qq_energy
    REAL(DP) :: rxp, ryp, rzp, sxp, syp, szp

    REAL(DP), DIMENSION(atompairdim) :: irxp, iryp, irzp, jrxp, jryp, jrzp, cfqq_vec, rij_vec, rijsq_packed
    REAL(DP), DIMENSION(atompairdim) :: vdw_p1,vdw_p2,vdw_p3,vdw_p4,vdw_p5
    REAL(DP), DIMENSION(MAX(atompairdim,MAXVAL(max_molecules))) :: rijsq_vec, rijsq_svec
    REAL(DP), DIMENSION(atompairdim,5) :: vdw_p, vdw_p_packed, vdw_p_nonzero
    REAL(DP) :: dsxp,dsyp,dszp,eps,sigsq,sigbyr2,sigbyr6,sigbyr12,rterm,rterm2
    REAL(DP) :: jsxp, jsyp, jszp, dsp
    REAL(DP) :: negsigsq,negsigbyr2,roffsq_rijsq
    REAL(DP) :: epsig_n,epsig_m,mie_n,mie_m
    REAL(DP) :: shift_p1, shift_p2, cfqq

    INTEGER :: chunksize, nthreads_used, im_thread, jbase, ithread, n_interact_p, cni, j_end_shared
    INTEGER(2) :: ji_end,ji_start,vlen,ji_end_2,ji_base_2
    LOGICAL(1), DIMENSION(MAX(atompairdim,mol_dim)) :: vdw_interact_vec, qq_interact_vec

    REAL(DP), DIMENSION(mol_dim,MAX(0,MAXVAL(natoms,MASK=natoms .LE. nmols(:,molecule_list(im,is)%which_box))),3) :: jrp_perm
    LOGICAL :: l_molvectorized,ij_get_vdw,ij_get_qq,l_getdrcom,l_vdw_packed,l_vdw_gather,l_vdw_strided,l_masked
    LOGICAL :: l_shortvdw, l_shortcoul, l_notsamecut
    REAL(DP) :: pair_max_rcutsq,pair_min_rcutsq,nrg,dsf_const,i_max_dcom
    INTEGER :: itype,jtype,maxvlen,minvlen,vdw_vlen,qq_vlen, n_interact_qq_pairstore, n_interact_vdw_pairstore
    INTEGER, DIMENSION(atompairdim) :: which_have_vdw
    INTEGER :: atompair_vlen, ji_base, ji_stride
    REAL(DP), DIMENSION(mol_dim) :: mol_vdw_energy_p, mol_qq_energy_p
    REAL(DP), DIMENSION(:), ALLOCATABLE :: mol_vdw_energy, mol_qq_energy

    REAL(DP) :: alpharsq, alpha_sq, ewald_const, numer, denom

    INTEGER :: i_selector
    LOGICAL :: l_read_svec, switchflag
    INTEGER(1) :: ji_start_byte, ji_byte

    !DIR$ ASSUME_ALIGNED jrxp:32, jryp:32, jrzp:32, irxp:32, iryp:32, irzp:32, cfqq_vec:32, rijsq_vec:32, rij_vec:32
    !DIR$ ASSUME_ALIGNED vdw_p:32, jrp_perm:32, vdw_p_packed:32, vdw_p_nonzero:32, interact_vec:32, all_interact_vec:32
    !DIR$ ASSUME_ALIGNED rijsq_packed:32
    !DIR$ ASSUME_ALIGNED vdw_p1:32, vdw_p2:32, vdw_p3:32, vdw_p4:32, vdw_p5:32
    !DIR$ ASSUME (MOD(atompairdim,8) .EQ. 0)
    !DIR$ ASSUME (MOD(mol_dim,8) .EQ. 0)

    i_vdw_energy = 0.0_DP
    i_qq_energy = 0.0_DP
    overlap = .FALSE.
    i_overlap = .FALSE.
    this_box = molecule_list(im,is)%which_box
    ibox = this_box
    n_vdw_p = 0
    n_i_exist = 0
    DO i = 1, natoms(is)
        IF (atom_list(i,im,is)%exist) THEN
                n_i_exist = n_i_exist + 1
                which_i_exist(n_i_exist) = i
        END IF
    END DO
    l_ortho = box_list(this_box)%int_box_shape <= int_ortho

    !CALL Set_Pair_Chunks

    ithread = 1

    nthreads_used = 1



    !$OMP PARALLEL PRIVATE(jm,jrxp,jryp,jrzp,dsp,dxp,dyp,dzp,i,j,k,ji,rijsq,rij) &
    !$OMP PRIVATE(alpharsq,numer,denom) &
    !$OMP PRIVATE(vdw_energy,qq_energy,jsl,ij_overlap) &
    !$OMP PRIVATE(jsxp,jsyp,jszp,eps,sigsq,sigbyr2,sigbyr6,sigbyr12,rterm,rterm2) &
    !$OMP PRIVATE(negsigsq,negsigbyr2,roffsq_rijsq) &
    !$OMP PRIVATE(epsig_n,epsig_m,mie_n,mie_m) &
    !$OMP PRIVATE(l_interact, shift_p1, shift_p2, cfqq) &
    !$OMP PRIVATE(rijsq_vec,rijsq_svec, rij_vec) &
    !$OMP PRIVATE(js,jnmols,get_qq,natompairs,jsl_base,jnlive,chunksize,nthreads_used) &
    !$OMP PRIVATE(j_limit,j_limit_2,im_thread,jbase,ithread) &
    !$OMP PRIVATE(dxcom,dycom,dzcom,max_dcom,dsxcom,dsycom,dszcom) &
    !$OMP PRIVATE(n_may_interact_p,which_may_interact_p,cni,need_sqrt) &
    !$OMP PRIVATE(n_all_interact_p,which_all_interact_p) &
    !$OMP PRIVATE(ji_end,ji_start,vlen,ji_end_2,ji_base_2) &
    !$OMP PRIVATE(rijsq_packed,vdw_p_packed) &
    !$OMP PRIVATE(l_molvectorized,ij_get_vdw,ij_get_qq,l_getdrcom) &
    !$OMP PRIVATE(l_vdw_packed,l_vdw_gather,l_vdw_strided,l_masked) &
    !$OMP PRIVATE(rxp,ryp,rzp,pair_max_rcutsq,pair_min_rcutsq,nrg) &
    !$OMP PRIVATE(itype,jtype,maxvlen,minvlen,vdw_vlen,qq_vlen,atompair_vlen) &
    !$OMP PRIVATE(qq_interact_vec,vdw_interact_vec) &
    !$OMP PRIVATE(mol_vdw_energy_p,mol_qq_energy_p) &
    !$OMP PRIVATE(l_read_svec, i_selector, switchflag) &
    !$OMP PRIVATE(ji_start_byte,ji_byte) &
    !$OMP PRIVATE(ji_base, ji_stride) &
    !$OMP PRIVATE(vdw_p1,vdw_p2,vdw_p3,vdw_p4,vdw_p5) &
    !$OMP PRIVATE(ia,ia0,ja,ibase) &
    !$OMP REDUCTION(+: i_vdw_energy,i_qq_energy) &
    !$OMP REDUCTION(.OR.: i_overlap)

    i_overlap = .FALSE.
    i_vdw_energy = 0.0_DP
    i_qq_energy = 0.0_DP


    !$ ithread = OMP_GET_THREAD_NUM() + 1

    !$ nthreads_used = OMP_GET_NUM_THREADS()

    !DIR$ ASSUME (MOD(atompairdim,8) .EQ. 0)
    !DIR$ ASSUME (MOD(mol_dim,8) .EQ. 0)

    !$OMP SECTIONS
    !$OMP SECTION
    l_order_ij = n_i_exist .GE. 4 .AND. l_ortho
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
            inv_h11 = box_list(this_box)%length_inv(1,1)
            inv_h21 = box_list(this_box)%length_inv(2,1)
            inv_h31 = box_list(this_box)%length_inv(3,1)
            inv_h12 = box_list(this_box)%length_inv(1,2)
            inv_h22 = box_list(this_box)%length_inv(2,2)
            inv_h32 = box_list(this_box)%length_inv(3,2)
            inv_h13 = box_list(this_box)%length_inv(1,3)
            inv_h23 = box_list(this_box)%length_inv(2,3)
            inv_h33 = box_list(this_box)%length_inv(3,3)
    END IF
    DO i = 1, n_i_exist
        IF (n_i_exist == natoms(is)) THEN
                rxp = atom_list(i,im,is)%rp(1)
                ryp = atom_list(i,im,is)%rp(2)
                rzp = atom_list(i,im,is)%rp(3)
        ELSE
                ia = which_i_exist(i)
                rxp = atom_list(ia,im,is)%rp(1)
                ryp = atom_list(ia,im,is)%rp(2)
                rzp = atom_list(ia,im,is)%rp(3)
        END IF
        IF (l_ortho) THEN
                irsxp_vec(i) = rxp
                irsyp_vec(i) = ryp
                irszp_vec(i) = rzp
        ELSE
                sxp = inv_h11*rxp
                syp = inv_h21*rxp
                szp = inv_h31*rxp
                sxp = sxp + inv_h12*ryp
                syp = syp + inv_h22*ryp
                szp = szp + inv_h32*ryp
                sxp = sxp + inv_h13*rzp
                syp = syp + inv_h23*rzp
                szp = szp + inv_h33*rzp
                irsxp_vec(i) = rxp
                irsyp_vec(i) = ryp
                irszp_vec(i) = rzp
        END IF
    END DO
    !$OMP SECTION
    lj_charmm = .FALSE.
    lj_cut_shift = .FALSE.
    lj_cut_switch = .FALSE.
    lj_cut_shift_force = .FALSE.
    lj_cut = .FALSE.
    mie_cut_shift = .FALSE.
    mie_cut = .FALSE.
    get_vdw = int_vdw_style(this_box) .NE. vdw_none
    IF (get_vdw) THEN
            this_int_vdw_style = int_vdw_style(this_box)
            IF (this_int_vdw_style == vdw_lj) THEN
                    SELECT CASE (int_vdw_sum_style(this_box))
                        CASE (vdw_charmm)
                                lj_charmm = .TRUE.
                                n_vdw_p = 2
                        CASE (vdw_cut_shift)
                                lj_cut_shift = .TRUE.
                                n_vdw_p = 3
                        CASE (vdw_cut_switch)
                                lj_cut_switch = .TRUE.
                                n_vdw_p = 2
                        CASE (vdw_cut_shift_force)
                                lj_cut_shift_force = .TRUE.
                                n_vdw_p = 4
                        CASE DEFAULT
                                lj_cut = .TRUE.
                                n_vdw_p = 2
                    END SELECT
            ELSE IF (this_int_vdw_style == vdw_mie) THEN
                    IF (int_vdw_sum_style(this_box) .EQ. vdw_cut_shift) THEN
                            mie_cut_shift = .TRUE.
                            n_vdw_p = 5
                    ELSE
                            mie_cut = .TRUE.
                            n_vdw_p = 4
                    END IF
            END IF
    END IF
    l_charge_ewald = .FALSE.
    l_charge_dsf = .FALSE.
    l_charge_cut = .FALSE.
    i_get_qq = int_charge_style(this_box) .NE. charge_none .AND. has_charge(is) .AND. &
            (species_list(is)%l_coul_cbmc .OR. .NOT. cbmc_flag)
    IF (i_get_qq) THEN
            SELECT CASE (int_charge_sum_style(this_box))
                CASE (charge_ewald)
                        l_charge_ewald = .TRUE.
                        !ewald_const = alpha_ewald(this_box)/rootPI
                        !alpha_sq = alpha_ewald(this_box)*alpha_ewald(this_box)
                CASE (charge_dsf)
                        l_charge_dsf = .TRUE.
                        dsf_const = -dsf_factor1(this_box) - rcut_coul(this_box)*dsf_factor2(this_box)
                CASE DEFAULT
                        l_charge_cut = .TRUE.
            END SELECT
    END IF
    IF (cbmc_flag) THEN
            mol_rcut = rcut_cbmc(this_box)
            this_vdw_rcutsq = rcut_cbmcsq(this_box)
            this_coul_rcutsq = rcut_cbmcsq(this_box)
            this_coul_rcut = rcut_cbmc(this_box)
            l_shortvdw = .FALSE.
            l_shortcoul = .FALSE.
            l_notsamecut = .FALSE.
    ELSE
            mol_rcut = rcut_max(this_box)
            IF (int_vdw_sum_style(this_box) == vdw_cut_switch) THEN
                    this_vdw_rcutsq = roff_switch_sq(this_box)
            ELSE IF (get_vdw) THEN
                    this_vdw_rcutsq = rcut_vdwsq(this_box)
            ELSE
                    this_vdw_rcutsq = 0.0_DP
            END IF
            IF (i_get_qq) THEN
                    this_coul_rcutsq = rcut_coulsq(this_box)
            ELSE
                    this_coul_rcutsq = 0.0_DP
            END IF
            !this_coul_rcut = rcut_coul(this_box)
            l_shortcoul = this_coul_rcutsq < this_vdw_rcutsq
            l_shortvdw = this_vdw_rcutsq < this_coul_rcutsq
            l_notsamecut = l_shortcoul .OR. l_shortvdw
    END IF
    i_max_dcom = molecule_list(im,is)%rcom(4)
    max_dcom_i_const = i_max_dcom + mol_rcut
    ixcom = molecule_list(im,is)%rcom(1)
    iycom = molecule_list(im,is)%rcom(2)
    izcom = molecule_list(im,is)%rcom(3)
    l_pair_store = l_pair_nrg .AND. .NOT. cbmc_flag
    IF (l_pair_store) THEN
            isl = species_list(is)%superlocate_base + im
    END IF
    !$OMP END SECTIONS
    DO js = 1, nspecies
            jnmols = nmols(js,ibox)
            IF (jnmols == 0) CYCLE
            get_qq = i_get_qq .AND. has_charge(js)
            jnatoms = natoms(js)
            natompairs = jnatoms*n_i_exist
            l_getdrcom = natompairs > 3 ! I don't exactly know what the threshold should be. Subject to change.
            IF (open_mc_flag .OR. ibox > 1) THEN
                    IF (l_pair_store) THEN
                            jsl_base = species_list(js)%superlocate_base
                            !$OMP WORKSHARE
                            pair_nrg_vdw(locate(1:jnmols,js,ibox)+jsl_base,isl) = 0.0_DP
                            pair_nrg_qq(locate(1:jnmols,js,ibox)+jsl_base,isl) = 0.0_DP
                            pair_nrg_vdw(isl,locate(1:jnmols,js,ibox)+jsl_base) = 0.0_DP
                            pair_nrg_qq(isl,locate(1:jnmols,js,ibox)+jsl_base) = 0.0_DP
                            !$OMP END WORKSHARE
                    END IF
                    IF (l_not_all_live) THEN
                            !$OMP SINGLE
                            molecule_list(im,is)%live = .FALSE.
                            jnlive = 0
                            DO j = 1, jnmols
                                jm = locate(j,js,ibox)
                                IF (molecule_list(jm,js)%live) THEN
                                        jnlive = jnlive + 1
                                        live_locates(jnlive) = jm
                                        IF (l_getdrcom) live_rcom(:,jnlive) = molecule_list(jm,js)%rcom
                                END IF
                            END DO
                            molecule_list(im,is)%live = .TRUE.
                            !$OMP END SINGLE
                    ELSE IF (is .EQ. js) THEN
                            jnlive = jnmols - 1
                            !$OMP SINGLE
                            j_end_shared = jnmols + 1
                            !$OMP END SINGLE
                            !$OMP DO SIMD SCHEDULE(STATIC) REDUCTION(MIN:j_end_shared)
                            DO j = 1, jnmols
                                IF (locate(j,js,ibox) .EQ. im) j_end_shared = MIN(j_end_shared,j)
                            END DO
                            !$OMP END DO SIMD
                            !$OMP WORKSHARE
                            live_locates(1:j_end_shared-1) = locate(1:j_end_shared-1,js,ibox)
                            live_locates(j_end_shared:jnlive) = locate(j_end_shared+1:jnmols,js,ibox)
                            !$OMP END WORKSHARE
                            IF (l_getdrcom) THEN
                                    !$OMP DO SCHEDULE(STATIC)
                                    DO j = 1, jnlive
                                        live_rcom(:,j) = molecule_list(live_locates(j),js)%rcom
                                    END DO
                                    !$OMP END DO
                            END IF
                    ELSE
                            jnlive = jnmols
                            !$OMP WORKSHARE
                            live_locates(1:jnlive) = locate(1:jnlive,js,ibox)
                            !$OMP END WORKSHARE NOWAIT
                            IF (l_getdrcom) THEN
                                    !$OMP DO SCHEDULE(STATIC)
                                    DO j = 1, jnlive
                                            live_rcom(:,j) = molecule_list(locate(j,js,ibox),js)%rcom
                                    END DO
                                    !$OMP END DO
                            END IF
                    END IF
            ELSE
                    jnlive = jnmols
                    IF (l_pair_store) THEN
                            jsl_base = species_list(js)%superlocate_base
                            !$OMP WORKSHARE
                            pair_nrg_vdw(jsl_base+1:jsl_base+jnmols,isl) = 0.0_DP
                            pair_nrg_qq(jsl_base+1:jsl_base+jnmols,isl) = 0.0_DP
                            pair_nrg_vdw(isl,jsl_base+1:jsl_base+jnmols) = 0.0_DP
                            pair_nrg_qq(isl,jsl_base+1:jsl_base+jnmols) = 0.0_DP
                            !$OMP END WORKSHARE
                    END IF
            END IF
            IF (l_getdrcom) THEN
                    !chunksize = chunksize_array(js,is,ibox)
                    !nthreads_used = nthreads_used_array(js,is,ibox)
                    !IF (is .EQ. js .AND. .NOT. open_mc_flag) THEN
                    !        !IF (MOD(im-1,chunksize) == 0) THEN
                    !        !        im_thread = 9999
                    !        !ELSE
                    !        !        im_thread = (im-1)/chunksize + 1
                    !        !END IF
                    !        im_thread = (im-1)/chunksize + 1
                    !        im_thread = MIN(nthreads_used,im_thread)
                    !        im_thread = MAX(1,im_thread)
                    !        IF (ithread .LE. im_thread) THEN
                    !                jbase = (ithread-1)*chunksize
                    !        ELSE
                    !                jbase = jnlive - (nthreads_used - ithread + 1)*chunksize
                    !                jbase = MAX(jbase,im)
                    !        END IF
                    !        IF (ithread == im_thread) THEN
                    !                j_limit = im - 1 - jbase
                    !                j_limit_2 = jnlive - (nthreads_used - im_thread)*chunksize - jbase
                    !                interact_vec(j_limit+1) = .FALSE.
                    !                all_interact_vec(j_limit+1) = .FALSE.
                    !        ELSE
                    !                j_limit = chunksize
                    !                j_limit_2 = chunksize
                    !        END IF
                    !ELSE
                    !        im_thread = 9999
                    !        jbase = (ithread-1)*chunksize
                    !        j_limit = chunksize
                    !        j_limit_2 = chunksize
                    !END IF
                    !!IF (l_debug_print .OR. (i_mcstep > 0 .AND. i_mcstep < 10 .AND. .NOT. cbmc_flag) ) WRITE(*,*) im_thread, ithread, im, jnlive, jbase, j_limit, j_limit_2, chunksize
                    !j_limit = MIN(j_limit,jnlive-jbase)
                    !j_limit_2 = MIN(j_limit_2,jnlive-jbase)
                    !!IF (l_debug_print .OR. (i_mcstep > 0 .AND. i_mcstep < 10 .AND. .NOT. cbmc_flag) ) WRITE(*,*) im_thread, ithread, im, jnlive, jbase, j_limit, j_limit_2, chunksize
                    IF (l_ortho) THEN
                            !$OMP DO SIMD SCHEDULE(SIMD:STATIC) PRIVATE(dxcom,dycom,dzcom,max_dcom)
                            DO j = 1, jnlive
                                    IF (open_mc_flag) THEN
                                            dxcom = live_rcom(1,j)
                                            dycom = live_rcom(2,j)
                                            dzcom = live_rcom(3,j)
                                            max_dcom = live_rcom(4,j)
                                    ELSE
                                            dxcom = molecule_list(j,js)%rcom(1)
                                            dycom = molecule_list(j,js)%rcom(2)
                                            dzcom = molecule_list(j,js)%rcom(3)
                                            max_dcom = molecule_list(j,js)%rcom(4)
                                    END IF
                                    dxcom = ABS(dxcom - ixcom)
                                    dycom = ABS(dycom - iycom)
                                    dzcom = ABS(dzcom - izcom)
                                    IF (dxcom > hxl) dxcom = dxcom - xl
                                    IF (dycom > hyl) dycom = dycom - yl
                                    IF (dzcom > hzl) dzcom = dzcom - zl
                                    ! Repurposing dxcom as rijsq accumulator to enforce FMA3 instructions if supported
                                    dxcom = dxcom * dxcom
                                    dxcom = dxcom + dycom * dycom
                                    dxcom = dxcom + dzcom * dzcom
                                    dxcom = SQRT(dxcom) ! rij
                                    max_dcom = max_dcom + i_max_dcom
                                    interact_vec(j) = mol_rcut > dxcom - max_dcom
                                    all_interact_vec(j) = mol_rcut > dxcom + max_dcom
                                    !l_interact = max_dcom_i_const > &
                                    !        SQRT(dxcom) - max_dcom
                            END DO
                            !$OMP END DO SIMD
                    ELSE
                            IF (open_mc_flag) THEN
                                    !$OMP DO SIMD SCHEDULE(SIMD:STATIC) PRIVATE(dxcom,dycom,dzcom,max_dcom,dsxcom,dsycom,dszcom)
                                    DO j = 1, jnlive
                                            dxcom = live_rcom(1,j)
                                            dycom = live_rcom(2,j)
                                            dzcom = live_rcom(3,j)
                                            max_dcom = live_rcom(4,j)
                                            dxcom = dxcom - ixcom
                                            dsxcom = inv_h11*dxcom
                                            dsycom = inv_h21*dxcom
                                            dszcom = inv_h31*dxcom
                                            dycom = dycom - iycom
                                            dsxcom = dsxcom + inv_h12*dycom
                                            dsycom = dsycom + inv_h22*dycom
                                            dszcom = dszcom + inv_h32*dycom
                                            dzcom = dzcom - izcom
                                            dsxcom = dsxcom + inv_h13*dzcom
                                            dsycom = dsycom + inv_h23*dzcom
                                            dszcom = dszcom + inv_h33*dzcom
                                            IF (dsxcom > 0.5_DP) THEN
                                                    dsxcom = dsxcom - 1.0_DP
                                            ELSE IF (dsxcom < -0.5_DP) THEN
                                                    dsxcom = dsxcom + 1.0_DP
                                            END IF
                                            IF (dsycom > 0.5_DP) THEN
                                                    dsycom = dsycom - 1.0_DP
                                            ELSE IF (dsycom < -0.5_DP) THEN
                                                    dsycom = dsycom + 1.0_DP
                                            END IF
                                            IF (dszcom > 0.5_DP) THEN
                                                    dszcom = dszcom - 1.0_DP
                                            ELSE IF (dszcom < -0.5_DP) THEN
                                                    dszcom = dszcom + 1.0_DP
                                            END IF
                                            dxcom = h11*dsxcom
                                            dycom = h21*dsxcom
                                            dzcom = h31*dsxcom
                                            dxcom = dxcom + h12*dsycom
                                            dycom = dycom + h22*dsycom
                                            dzcom = dzcom + h32*dsycom
                                            dxcom = dxcom + h13*dszcom
                                            dycom = dycom + h23*dszcom
                                            dzcom = dzcom + h33*dszcom
                                            ! Repurposing dxcom as rijsq accumulator to enforce FMA3 instructions if supported
                                            dxcom = dxcom * dxcom
                                            dxcom = dxcom + dycom * dycom
                                            dxcom = dxcom + dzcom * dzcom
                                            dxcom = SQRT(dxcom) ! rij
                                            max_dcom = max_dcom + i_max_dcom
                                            interact_vec(j) = mol_rcut > dxcom - max_dcom
                                            all_interact_vec(j) = mol_rcut > dxcom + max_dcom
                                    END DO
                                    !$OMP END DO SIMD
                            ELSE
                                    !$OMP DO SIMD SCHEDULE(SIMD:STATIC) PRIVATE(dxcom,dycom,dzcom,max_dcom,dsxcom,dsycom,dszcom)
                                    DO j = 1, jnlive
                                            dxcom = molecule_list(j,js)%rcom(1)
                                            dycom = molecule_list(j,js)%rcom(2)
                                            dzcom = molecule_list(j,js)%rcom(3)
                                            max_dcom = molecule_list(j,js)%rcom(4)
                                            dxcom = dxcom - ixcom
                                            dsxcom = inv_h11*dxcom
                                            dsycom = inv_h21*dxcom
                                            dszcom = inv_h31*dxcom
                                            dycom = dycom - iycom
                                            dsxcom = dsxcom + inv_h12*dycom
                                            dsycom = dsycom + inv_h22*dycom
                                            dszcom = dszcom + inv_h32*dycom
                                            dzcom = dzcom - izcom
                                            dsxcom = dsxcom + inv_h13*dzcom
                                            dsycom = dsycom + inv_h23*dzcom
                                            dszcom = dszcom + inv_h33*dzcom
                                            IF (dsxcom > 0.5_DP) THEN
                                                    dsxcom = dsxcom - 1.0_DP
                                            ELSE IF (dsxcom < -0.5_DP) THEN
                                                    dsxcom = dsxcom + 1.0_DP
                                            END IF
                                            IF (dsycom > 0.5_DP) THEN
                                                    dsycom = dsycom - 1.0_DP
                                            ELSE IF (dsycom < -0.5_DP) THEN
                                                    dsycom = dsycom + 1.0_DP
                                            END IF
                                            IF (dszcom > 0.5_DP) THEN
                                                    dszcom = dszcom - 1.0_DP
                                            ELSE IF (dszcom < -0.5_DP) THEN
                                                    dszcom = dszcom + 1.0_DP
                                            END IF
                                            dxcom = h11*dsxcom
                                            dycom = h21*dsxcom
                                            dzcom = h31*dsxcom
                                            dxcom = dxcom + h12*dsycom
                                            dycom = dycom + h22*dsycom
                                            dzcom = dzcom + h32*dsycom
                                            dxcom = dxcom + h13*dszcom
                                            dycom = dycom + h23*dszcom
                                            dzcom = dzcom + h33*dszcom
                                            ! Repurposing dxcom as rijsq accumulator to enforce FMA3 instructions if supported
                                            dxcom = dxcom * dxcom
                                            dxcom = dxcom + dycom * dycom
                                            dxcom = dxcom + dzcom * dzcom
                                            dxcom = SQRT(dxcom) ! rij
                                            max_dcom = max_dcom + i_max_dcom
                                            interact_vec(j) = mol_rcut > dxcom - max_dcom
                                            all_interact_vec(j) = mol_rcut > dxcom + max_dcom
                                    END DO
                                    !$OMP END DO SIMD
                            END IF
                    END IF
                    IF (is == js .AND. .NOT. open_mc_flag) THEN
                            !$OMP SINGLE
                            interact_vec(im) = .FALSE.
                            all_interact_vec(im) = .FALSE.
                            !$OMP END SINGLE
                    END IF
                    IF (nthreads_used > 1) THEN
                            !$OMP WORKSHARE
                            n_interact = COUNT(interact_vec(1:jnlive))
                            n_all_interact = COUNT(all_interact_vec(1:jnlive))
                            n_may_interact = n_interact - n_all_interact
                            !$OMP END WORKSHARE
                    END IF
                    !$OMP SINGLE
                    n_all_interact_p = 0
                    n_may_interact_p = 0
                    IF (open_mc_flag) THEN
                            DO j = 1, jnlive
                                IF (all_interact_vec(j)) THEN
                                        n_all_interact_p = n_all_interact_p + 1
                                        which_interact(n_all_interact_p) = live_locates(j)
                                ELSE IF (interact_vec(j)) THEN
                                        n_may_interact_p = n_may_interact_p + 1
                                        which_may_interact_p(n_may_interact_p) = live_locates(j)
                                END IF
                            END DO
                    ELSE
                            DO j = 1, jnlive
                                IF (all_interact_vec(j)) THEN
                                        n_all_interact_p = n_all_interact_p + 1
                                        which_interact(n_all_interact_p) = j
                                ELSE IF (interact_vec(j)) THEN
                                        n_may_interact_p = n_may_interact_p + 1
                                        which_may_interact_p(n_may_interact_p) = j
                                END IF
                            END DO
                    END IF
                    ! The first n_all_interact molecules in which_interact are guaranteed to be fully in-range.
                    ! The next n_may_interact molecules might have some atoms in range.
                    which_interact(n_all_interact_p+1:n_all_interact_p+n_may_interact_p) = which_may_interact_p(1:n_may_interact_p)
                    IF (nthreads_used < 2) THEN
                            n_all_interact = n_all_interact_p
                            n_may_interact = n_may_interact_p
                            n_interact = n_all_interact_p + n_may_interact_p
                    END IF
                    !$OMP END SINGLE NOWAIT
                    l_molvectorized = n_interact > natompairs
                    IF (l_debug_print) WRITE(*,*) "n_interact = ", n_interact, n_all_interact, n_may_interact, l_molvectorized
                    IF (l_debug_print) WRITE(*,*) n_all_interact_p, n_may_interact_p
                    IF (l_debug_print) WRITE(*,*) "nthreads_used = ", nthreads_used
            ELSE
                    l_molvectorized = .TRUE.
            END IF
            IF (l_order_ij .AND. .NOT. l_molvectorized) THEN
                    !$OMP WORKSHARE
                    iatomtypes = nonbond_list(which_i_exist,is)%atom_type_number
                    jatomtypes(1:jnatoms) = nonbond_list(1:jnatoms,js)%atom_type_number
                    vdw_p(1:natompairs,1:n_vdw_p) = RESHAPE(ppvdwp_table(iatomtypes,jatomtypes(1:jnatoms),1:n_vdw_p,ibox), (/ natompairs, n_vdw_p /))
                    irxp(1:natompairs) = RESHAPE(SPREAD(irsxp_vec(1:n_i_exist),2,jnatoms), (/ natompairs /))
                    iryp(1:natompairs) = RESHAPE(SPREAD(irsyp_vec(1:n_i_exist),2,jnatoms), (/ natompairs /))
                    irzp(1:natompairs) = RESHAPE(SPREAD(irszp_vec(1:n_i_exist),2,jnatoms), (/ natompairs /))
                    !$OMP END WORKSHARE NOWAIT
                    IF (get_qq) THEN
                            !$OMP WORKSHARE
                            cfqq_vec(1:natompairs) = RESHAPE(SPREAD(nonbond_list(1:jnatoms,js)%charge,1,n_i_exist) * &
                                    SPREAD(nonbond_list(which_i_exist,is)%charge,2,jnatoms) * charge_factor, (/ natompairs /))
                            !$OMP END WORKSHARE NOWAIT
                    END IF
            ELSE IF (.NOT. l_molvectorized) THEN
                    !$OMP WORKSHARE
                    iatomtypes = nonbond_list(which_i_exist,is)%atom_type_number
                    jatomtypes(1:jnatoms) = nonbond_list(1:jnatoms,js)%atom_type_number
                    vdw_p(1:natompairs,1:n_vdw_p) = RESHAPE(ppvdwp_table(jatomtypes(1:jnatoms),iatomtypes,1:n_vdw_p,ibox), (/ natompairs, n_vdw_p /))
                    irxp(1:natompairs) = RESHAPE(SPREAD(irsxp_vec(1:n_i_exist),1,jnatoms), (/ natompairs /))
                    iryp(1:natompairs) = RESHAPE(SPREAD(irsyp_vec(1:n_i_exist),1,jnatoms), (/ natompairs /))
                    irzp(1:natompairs) = RESHAPE(SPREAD(irszp_vec(1:n_i_exist),1,jnatoms), (/ natompairs /))
                    !$OMP END WORKSHARE NOWAIT
                    IF (get_qq) THEN
                            !$OMP WORKSHARE
                            cfqq_vec(1:natompairs) = RESHAPE(SPREAD(nonbond_list(1:jnatoms,js)%charge,2,n_i_exist) * &
                                    SPREAD(nonbond_list(which_i_exist,is)%charge,1,jnatoms) * charge_factor, (/ natompairs /))
                            !$OMP END WORKSHARE NOWAIT
                    END IF
            END IF

            need_sqrt = get_qq .OR. lj_cut_shift_force

            !$OMP BARRIER
            !IF (i_mcstep < 10 .AND. .NOT. cbmc_flag) WRITE(*,*) which_interact(1:n_interact)
            IF (l_molvectorized) THEN
                    IF (l_debug_print) WRITE(*,*) "l_molvectorized"
                    IF (l_getdrcom) THEN
                            !WRITE(*,*) SHAPE(jrp_perm)
                            !WRITE(*,*) n_interact, jnatoms
                            !$OMP WORKSHARE
                            jrp_perm(1:n_interact,1:jnatoms,1) = TRANSPOSE(atom_list(1:jnatoms,which_interact(1:n_interact),js)%rp(1))
                            jrp_perm(1:n_interact,1:jnatoms,2) = TRANSPOSE(atom_list(1:jnatoms,which_interact(1:n_interact),js)%rp(2))
                            jrp_perm(1:n_interact,1:jnatoms,3) = TRANSPOSE(atom_list(1:jnatoms,which_interact(1:n_interact),js)%rp(3))
                            !$OMP END WORKSHARE
                            IF (l_debug_print) THEN
                                    WRITE(*,*) im
                                    WRITE(*,*) irsxp_vec
                                    WRITE(*,*) irsyp_vec
                                    WRITE(*,*) irszp_vec
                                    WRITE(*,*)
                                    WRITE(*,*) n_interact, jnatoms
                                    DO i = 1, n_interact
                                        WRITE(*,*) "Molecule ", which_interact(i)
                                        DO j = 1, jnatoms
                                                WRITE(*,*) jrp_perm(i,j,1:3)
                                        END DO
                                    END DO
                            END IF
                    ELSE IF (open_mc_flag) THEN
                            !$OMP WORKSHARE
                            n_interact = jnlive
                            n_all_interact = 0
                            n_may_interact = n_interact
                            jrp_perm(1:jnlive,1:jnatoms,1) = TRANSPOSE(atom_list(1:jnatoms,live_locates(1:jnlive),js)%rp(1))
                            jrp_perm(1:jnlive,1:jnatoms,2) = TRANSPOSE(atom_list(1:jnatoms,live_locates(1:jnlive),js)%rp(2))
                            jrp_perm(1:jnlive,1:jnatoms,3) = TRANSPOSE(atom_list(1:jnatoms,live_locates(1:jnlive),js)%rp(3))
                            !$OMP END WORKSHARE
                    ELSE
                            IF (is .EQ. js) THEN
                                    !$OMP SINGLE
                                    n_interact = jnlive - 1
                                    n_all_interact = 0
                                    n_may_interact = n_interact
                                    !$OMP END SINGLE NOWAIT
                                    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
                                    DO ja = 1, jnatoms
                                        DO i = 1, 3
                                                jrp_perm(1:im-1,ja,i) = atom_list(ja,1:im-1,js)%rp(i)
                                                jrp_perm(im:jnlive-1,ja,i) = atom_list(ja,im+1:jnlive,js)%rp(i)
                                        END DO
                                    END DO
                                    !$OMP END DO
                            ELSE IF (l_ortho) THEN
                                    !$OMP SINGLE
                                    n_interact = jnlive
                                    n_all_interact = 0
                                    n_may_interact = n_interact
                                    !$OMP END SINGLE NOWAIT
                                    !$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
                                    DO ja = 1, jnatoms
                                        DO i = 1, 3
                                                jrp_perm(1:jnlive,ja,i) = atom_list(ja,1:jnlive,js)%rp(i)
                                        END DO
                                    END DO
                                    !$OMP END DO
                            ELSE
                                    !$OMP SINGLE
                                    n_interact = jnlive
                                    n_all_interact = 0
                                    n_may_interact = n_interact
                                    !$OMP END SINGLE
                            END IF
                    END IF
                    IF (.NOT. l_ortho) THEN
                            switchflag = l_getdrcom .OR. open_mc_flag .OR. is .EQ. js
                            !$OMP DO SCHEDULE(STATIC)
                            DO ja = 1, jnatoms
                                DO j = 1, n_interact
                                        IF (switchflag) THEN
                                                rxp = jrp_perm(j,ja,1)
                                                ryp = jrp_perm(j,ja,2)
                                                rzp = jrp_perm(j,ja,3)
                                        ELSE
                                                rxp = atom_list(ja,j,js)%rp(1)
                                                ryp = atom_list(ja,j,js)%rp(2)
                                                rzp = atom_list(ja,j,js)%rp(3)
                                        END IF
                                        jsxp = inv_h11*rxp
                                        jsyp = inv_h21*rxp
                                        jszp = inv_h31*rxp
                                        jsxp = jsxp + inv_h12*ryp
                                        jsyp = jsyp + inv_h22*ryp
                                        jszp = jszp + inv_h32*ryp
                                        jsxp = jsxp + inv_h13*rzp
                                        jsyp = jsyp + inv_h23*rzp
                                        jszp = jszp + inv_h33*rzp
                                        jrp_perm(j,ja,1) = jsxp
                                        jrp_perm(j,ja,2) = jsyp
                                        jrp_perm(j,ja,3) = jszp
                                END DO
                            END DO
                            !$OMP END DO NOWAIT
                    END IF
                    !$OMP SINGLE
                    IF (get_vdw .AND. l_pair_store) THEN
                            n_interact_vdw_pairstore = n_interact
                    ELSE
                            n_interact_vdw_pairstore = 1 ! reduction item cannot have zero length
                    END IF
                    IF (get_qq .AND. l_pair_store) THEN
                            n_interact_qq_pairstore = n_interact
                    ELSE
                            n_interact_qq_pairstore = 1 ! reduction item cannot have zero length
                    END IF
                    IF (.NOT. ALLOCATED(mol_vdw_energy)) THEN
                            ALLOCATE(mol_vdw_energy(n_interact_vdw_pairstore))
                    ELSE IF (SIZE(mol_vdw_energy,1) < n_interact_vdw_pairstore) THEN
                            DEALLOCATE(mol_vdw_energy)
                            ALLOCATE(mol_vdw_energy(n_interact_vdw_pairstore))
                    END IF
                    IF (.NOT. ALLOCATED(mol_qq_energy)) THEN
                            ALLOCATE(mol_qq_energy(n_interact_qq_pairstore))
                    ELSE IF (SIZE(mol_qq_energy,1) < n_interact_qq_pairstore) THEN
                            DEALLOCATE(mol_qq_energy)
                            ALLOCATE(mol_qq_energy(n_interact_qq_pairstore))
                    END IF
                    mol_vdw_energy = 0.0_DP
                    mol_qq_energy = 0.0_DP
                    !$OMP END SINGLE
                    !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC) &
                    !$OMP REDUCTION(+:mol_vdw_energy,mol_qq_energy)
                    DO ja = 1, jnatoms
                        DO ia = 1, n_i_exist
                                !DIR$ ASSUME (MOD(atompairdim,8) .EQ. 0)
                                !DIR$ ASSUME (MOD(mol_dim,8) .EQ. 0)
                                IF (overlap .OR. i_overlap) CYCLE
                                ia0 = which_i_exist(ia)
                                itype = nonbond_list(ia0,is)%atom_type_number
                                jtype = nonbond_list(ja,js)%atom_type_number
                                ij_get_vdw = itype > 0 .AND. jtype > 0 .AND. get_vdw
                                IF (i_get_qq) THEN
                                        cfqq = nonbond_list(ia0,is)%charge * nonbond_list(ja,js)%charge * charge_factor
                                        ij_get_qq = ABS(cfqq) > 0.0_DP
                                ELSE
                                        ij_get_qq = .FALSE.
                                END IF
                                IF (l_ortho) THEN
                                        DO j = 1, n_interact
                                                dxp = ABS(jrp_perm(j,ja,1) - irsxp_vec(ia))
                                                dyp = ABS(jrp_perm(j,ja,2) - irsyp_vec(ia))
                                                dzp = ABS(jrp_perm(j,ja,3) - irszp_vec(ia))
                                                IF (dxp > hxl) dxp = dxp - xl
                                                IF (dyp > hyl) dyp = dyp - yl
                                                IF (dzp > hzl) dzp = dzp - zl
                                                dxp = dxp*dxp
                                                dxp = dxp + dyp*dyp
                                                dxp = dxp + dzp*dzp
                                                rijsq_vec(j) = dxp
                                                IF (ij_get_vdw) vdw_interact_vec(j) = dxp < this_vdw_rcutsq
                                                IF (ij_get_qq) qq_interact_vec(j) = dxp < this_coul_rcutsq
                                                i_overlap = i_overlap .OR. dxp < rcut_lowsq
                                        END DO
                                ELSE
                                        DO j = 1, n_interact
                                                dsp = jrp_perm(j,ja,1) - irsxp_vec(ia)
                                                IF (dsp > 0.5_DP) THEN
                                                        dsp = dsp - 1.0_DP
                                                ELSE IF (dsp < -0.5_DP) THEN
                                                        dsp = dsp + 1.0_DP
                                                END IF
                                                dxp = h11*dsp
                                                dyp = h21*dsp
                                                dzp = h31*dsp
                                                dsp = jrp_perm(j,ja,2) - irsyp_vec(ia)
                                                IF (dsp > 0.5_DP) THEN
                                                        dsp = dsp - 1.0_DP
                                                ELSE IF (dsp < -0.5_DP) THEN
                                                        dsp = dsp + 1.0_DP
                                                END IF
                                                dxp = dxp + h12*dsp
                                                dyp = dyp + h22*dsp
                                                dzp = dzp + h32*dsp
                                                dsp = jrp_perm(j,ja,3) - irszp_vec(ia)
                                                IF (dsp > 0.5_DP) THEN
                                                        dsp = dsp - 1.0_DP
                                                ELSE IF (dsp < -0.5_DP) THEN
                                                        dsp = dsp + 1.0_DP
                                                END IF
                                                dxp = dxp + h13*dsp
                                                dyp = dyp + h23*dsp
                                                dzp = dzp + h33*dsp
                                                dxp = dxp*dxp
                                                dxp = dxp + dyp*dyp
                                                dxp = dxp + dzp*dzp
                                                rijsq_vec(j) = dxp
                                                IF (ij_get_vdw) vdw_interact_vec(j) = dxp < this_vdw_rcutsq
                                                IF (ij_get_qq) qq_interact_vec(j) = dxp < this_coul_rcutsq
                                                i_overlap = i_overlap .OR. dxp < rcut_lowsq
                                        END DO
                                END IF
                                IF (l_debug_print) THEN
                                        WRITE(*,*) "i_overlap = ", i_overlap
                                        WRITE(*,*) "rijsq_vec = "
                                        WRITE(*,*) rijsq_vec(1:n_interact)
                                END IF
                                IF (i_overlap) THEN
                                        overlap = .TRUE.
                                        CYCLE
                                END IF
                                IF (ij_get_vdw .AND. ij_get_qq) THEN
                                        pair_max_rcutsq = MAX(this_vdw_rcutsq,this_coul_rcutsq)
                                        pair_min_rcutsq = MIN(this_vdw_rcutsq,this_coul_rcutsq)
                                        n_all_interact_p = n_all_interact
                                ELSE IF (ij_get_vdw) THEN
                                        pair_max_rcutsq = this_vdw_rcutsq
                                        IF (l_shortvdw) THEN
                                                n_all_interact_p = 0
                                        ELSE
                                                n_all_interact_p = n_all_interact
                                        END IF
                                ELSE IF (ij_get_qq) THEN
                                        pair_max_rcutsq = this_coul_rcutsq
                                        IF (l_shortcoul) THEN
                                                n_all_interact_p = 0
                                        ELSE
                                                n_all_interact_p = n_all_interact
                                        END IF
                                END IF
                                maxvlen = n_all_interact_p
                                DO j = n_all_interact_p+1, n_interact
                                        rijsq = rijsq_vec(j)
                                        IF (rijsq < pair_max_rcutsq) THEN
                                                maxvlen = maxvlen + 1
                                                rijsq_vec(maxvlen) = rijsq
                                        END IF
                                END DO
                                IF (ij_get_vdw .AND. ij_get_qq .AND. l_notsamecut) THEN
                                        minvlen = 0
                                        DO j = 1, maxvlen
                                                rijsq = rijsq_vec(j)
                                                IF (rijsq < pair_min_rcutsq) THEN
                                                        minvlen = minvlen + 1
                                                        rijsq_svec(minvlen) = rijsq
                                                END IF
                                        END DO
                                        IF (l_shortcoul) THEN
                                                vdw_vlen = maxvlen
                                                qq_vlen = minvlen
                                        ELSE
                                                vdw_vlen = minvlen
                                                qq_vlen = maxvlen
                                        END IF
                                ELSE
                                        IF (ij_get_vdw) vdw_vlen = maxvlen
                                        IF (ij_get_qq) qq_vlen = maxvlen
                                END IF
                                !IF (i_mcstep == 1) THEN
                                !        WRITE(*,*) l_notsamecut, l_shortcoul, l_shortvdw
                                !        WRITE(*,*) itype, jtype, cfqq
                                !        WRITE(*,*) minvlen, maxvlen, n_all_interact_p, n_interact, vdw_vlen, qq_vlen, ij_get_vdw, ij_get_qq
                                !END IF
                                IF (ij_get_vdw) THEN
                                        l_read_svec = ij_get_qq .AND. l_shortvdw
                                        IF (l_debug_print) WRITE(*,*) "l_read_svec = ", l_read_svec
                                        IF (lj_charmm) THEN
                                                eps = ppvdwp_table2(1,itype,jtype,ibox) ! epsilon
                                                sigsq = ppvdwp_table2(2,itype,jtype,ibox) ! sigma**2
                                                vdw_energy = 0.0_DP
                                                !$OMP SIMD PRIVATE(rijsq,sigbyr2,sigbyr6,sigbyr12) REDUCTION(+:vdw_energy)
                                                DO j = 1, vdw_vlen
                                                        IF (l_read_svec) THEN
                                                                rijsq = rijsq_svec(j)
                                                        ELSE
                                                                rijsq = rijsq_vec(j)
                                                        END IF
                                                        sigbyr2 = sigsq/rijsq ! sigma was already squared
                                                        sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                                        IF (l_pair_store) THEN
                                                                sigbyr12 = sigbyr6*sigbyr6
                                                                sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
                                                                mol_vdw_energy_p(j) = eps*sigbyr12
                                                                vdw_energy = vdw_energy + sigbyr12
                                                        ELSE
                                                                vdw_energy = vdw_energy + sigbyr6*sigbyr6 !- 2.0_DP *sigbyr6
                                                                vdw_energy = vdw_energy - 2.0_DP*sigbyr6
                                                        END IF
                                                END DO
                                                !$OMP END SIMD

                                                !IF (l_pair_store) THEN
                                                !        IF (ij_get_qq .AND. l_shortvdw) THEN
                                                !                DO j = 1, vdw_vlen
                                                !                        rijsq = rijsq_svec(j)
                                                !                        sigbyr2 = sigsq/rijsq ! sigma was already squared
                                                !                        sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                                !                        sigbyr12 = sigbyr6*sigbyr6
                                                !                        sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
                                                !                        mol_vdw_energy_p(j) = eps*sigbyr12
                                                !                        vdw_energy = vdw_energy + sigbyr12
                                                !                END DO
                                                !        ELSE
                                                !                DO j = 1, vdw_vlen
                                                !                        rijsq = rijsq_vec(j)
                                                !                        sigbyr2 = sigsq/rijsq ! sigma was already squared
                                                !                        sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                                !                        sigbyr12 = sigbyr6*sigbyr6
                                                !                        sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
                                                !                        mol_vdw_energy_p(j) = eps*sigbyr12
                                                !                        vdw_energy = vdw_energy + sigbyr12
                                                !                END DO
                                                !        END IF
                                                !ELSE
                                                !        IF (ij_get_qq .AND. l_shortvdw) THEN
                                                !                DO j = 1, vdw_vlen
                                                !                        rijsq = rijsq_svec(j)
                                                !                        sigbyr2 = sigsq/rijsq ! sigma was already squared
                                                !                        sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                                !                        vdw_energy = vdw_energy + sigbyr6*sigbyr6
                                                !                        vdw_energy = vdw_energy - 2.0_DP*sigbyr6
                                                !                END DO
                                                !        ELSE
                                                !                DO j = 1, vdw_vlen
                                                !                        rijsq = rijsq_vec(j)
                                                !                        sigbyr2 = sigsq/rijsq ! sigma was already squared
                                                !                        sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                                !                        vdw_energy = vdw_energy + sigbyr6*sigbyr6
                                                !                        vdw_energy = vdw_energy - 2.0_DP*sigbyr6
                                                !                END DO
                                                !        END IF
                                                !END IF
                                                i_vdw_energy = i_vdw_energy + vdw_energy * eps
                                        ELSE IF (lj_cut_shift) THEN
                                                eps = ppvdwp_table2(1,itype,jtype,ibox) ! 4*epsilon
                                                negsigsq = ppvdwp_table2(2,itype,jtype,ibox) ! -(sigma**2)
                                                shift_p1 = ppvdwp_table2(3,itype,jtype,ibox)
                                                vdw_energy = 0.0_DP
                                                IF (l_read_svec) THEN
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        vdw_energy = vdw_energy + rterm
                                                                        mol_vdw_energy_p(j) = eps*rterm - shift_p1
                                                                END DO
                                                        ELSE
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        vdw_energy = vdw_energy + rterm
                                                                END DO
                                                        END IF
                                                ELSE
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        vdw_energy = vdw_energy + rterm
                                                                        mol_vdw_energy_p(j) = eps*rterm - shift_p1
                                                                END DO
                                                        ELSE
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        vdw_energy = vdw_energy + rterm
                                                                END DO
                                                        END IF
                                                END IF
                                                !DO j = 1, vdw_vlen
                                                !        IF (l_read_svec) THEN
                                                !                rijsq = rijsq_svec(j)
                                                !        ELSE
                                                !                rijsq = rijsq_vec(j)
                                                !        END IF
                                                !        negsigbyr2 = negsigsq/rijsq
                                                !        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                !        rterm = rterm + rterm*rterm
                                                !        vdw_energy = vdw_energy + rterm
                                                !        IF (l_pair_store) THEN
                                                !                mol_vdw_energy_p(j) = eps*rterm - shift_p1
                                                !        END IF
                                                !        !vdw_energy = vdw_energy + eps * rterm - shift_p1
                                                !END DO
                                                i_vdw_energy = i_vdw_energy + vdw_energy*eps - vdw_vlen*shift_p1
                                        ELSE IF (lj_cut_switch) THEN
                                                eps = ppvdwp_table2(1,itype,jtype,ibox) ! 4*epsilon
                                                negsigsq = ppvdwp_table2(2,itype,jtype,ibox) ! -(sigma**2)
                                                vdw_energy = 0.0_DP
                                                IF (l_read_svec) THEN
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        roffsq_rijsq = roff_switch_sq(this_box) - rijsq
                                                                        roffsq_rijsq = roffsq_rijsq * roffsq_rijsq * switch_factor1(this_box) * &
                                                                                (switch_factor2(this_box)+2.0_DP*rijsq)
                                                                        roffsq_rijsq = MERGE(roffsq_rijsq, 1.0_DP, rijsq .GE. ron_switch_sq(this_box))
                                                                        vdw_energy = vdw_energy + roffsq_rijsq*rterm
                                                                        mol_vdw_energy_p(j) = eps*roffsq_rijsq*rterm
                                                                END DO
                                                        ELSE
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        roffsq_rijsq = roff_switch_sq(this_box) - rijsq
                                                                        roffsq_rijsq = roffsq_rijsq * roffsq_rijsq * switch_factor1(this_box) * &
                                                                                (switch_factor2(this_box)+2.0_DP*rijsq)
                                                                        roffsq_rijsq = MERGE(roffsq_rijsq, 1.0_DP, rijsq .GE. ron_switch_sq(this_box))
                                                                        vdw_energy = vdw_energy + roffsq_rijsq*rterm
                                                                END DO
                                                        END IF
                                                ELSE
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        roffsq_rijsq = roff_switch_sq(this_box) - rijsq
                                                                        roffsq_rijsq = roffsq_rijsq * roffsq_rijsq * switch_factor1(this_box) * &
                                                                                (switch_factor2(this_box)+2.0_DP*rijsq)
                                                                        roffsq_rijsq = MERGE(roffsq_rijsq, 1.0_DP, rijsq .GE. ron_switch_sq(this_box))
                                                                        vdw_energy = vdw_energy + roffsq_rijsq*rterm
                                                                        mol_vdw_energy_p(j) = eps*roffsq_rijsq*rterm
                                                                END DO
                                                        ELSE
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        roffsq_rijsq = roff_switch_sq(this_box) - rijsq
                                                                        roffsq_rijsq = roffsq_rijsq * roffsq_rijsq * switch_factor1(this_box) * &
                                                                                (switch_factor2(this_box)+2.0_DP*rijsq)
                                                                        roffsq_rijsq = MERGE(roffsq_rijsq, 1.0_DP, rijsq .GE. ron_switch_sq(this_box))
                                                                        vdw_energy = vdw_energy + roffsq_rijsq*rterm
                                                                END DO
                                                        END IF
                                                END IF
                                                !DO j = 1, vdw_vlen
                                                !        IF (l_read_svec) THEN
                                                !                rijsq = rijsq_svec(j)
                                                !        ELSE
                                                !                rijsq = rijsq_vec(j)
                                                !        END IF
                                                !        negsigbyr2 = negsigsq/rijsq
                                                !        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                !        rterm = rterm + rterm*rterm
                                                !        roffsq_rijsq = roff_switch_sq(this_box) - rijsq
                                                !        roffsq_rijsq = roffsq_rijsq * roffsq_rijsq * switch_factor1(this_box) * &
                                                !                (switch_factor2(this_box)+2.0_DP*rijsq)
                                                !        roffsq_rijsq = MERGE(roffsq_rijsq, 1.0_DP, rijsq .GE. ron_switch_sq(this_box))
                                                !        vdw_energy = vdw_energy + roffsq_rijsq*rterm
                                                !        IF (l_pair_store) THEN
                                                !                mol_vdw_energy_p(j) = eps*roffsq_rijsq*rterm
                                                !        END IF
                                                !END DO
                                                i_vdw_energy = i_vdw_energy + vdw_energy*eps
                                        ELSE IF (lj_cut_shift_force) THEN
                                                eps = ppvdwp_table2(1,itype,jtype,ibox) ! 4*epsilon
                                                negsigsq = ppvdwp_table2(2,itype,jtype,ibox) ! -(sigma**2)
                                                !shift_p1 = ppvdwp_table2(3,itype,jtype,ibox)
                                                !vdw_energy = shift_p1
                                                vdw_energy = ppvdwp_table2(3,itype,jtype,ibox)
                                                shift_p2 = ppvdwp_table2(4,itype,jtype,ibox)
                                                vdw_energy = vdw_energy + rcut_vdw(this_box)*shift_p2
                                                IF (l_pair_store) shift_p1 = vdw_energy
                                                vdw_energy = -vdw_vlen*vdw_energy
                                                IF (l_read_svec) THEN
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        nrg = eps*rterm + SQRT(rijsq)*shift_p2
                                                                        vdw_energy = vdw_energy + nrg
                                                                        mol_vdw_energy_p(j) = nrg - shift_p1
                                                                END DO
                                                        ELSE
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        vdw_energy = vdw_energy + eps*rterm + SQRT(rijsq)*shift_p2
                                                                END DO
                                                        END IF
                                                ELSE
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        nrg = eps*rterm + SQRT(rijsq)*shift_p2
                                                                        vdw_energy = vdw_energy + nrg
                                                                        mol_vdw_energy_p(j) = nrg - shift_p1
                                                                END DO
                                                        ELSE
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        vdw_energy = vdw_energy + eps*rterm + SQRT(rijsq)*shift_p2
                                                                END DO
                                                        END IF
                                                END IF
                                                !DO j = 1, vdw_vlen
                                                !        IF (l_read_svec) THEN
                                                !                rijsq = rijsq_svec(j)
                                                !        ELSE
                                                !                rijsq = rijsq_vec(j)
                                                !        END IF
                                                !        negsigbyr2 = negsigsq/rijsq
                                                !        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                !        rterm = rterm + rterm*rterm
                                                !        IF (l_pair_store) THEN
                                                !                nrg = eps*rterm + SQRT(rijsq)*shift_p2
                                                !                vdw_energy = vdw_energy + nrg
                                                !                mol_vdw_energy_p(j) = nrg - shift_p1
                                                !        ELSE
                                                !                vdw_energy = vdw_energy + eps*rterm + SQRT(rijsq)*shift_p2
                                                !        END IF
                                                !        !vdw_energy = vdw_energy + eps * rterm - shift_p1 - &
                                                !        !        (rcut_vdw(this_box)-SQRT(rijsq))*shift_p2
                                                !END DO
                                                i_vdw_energy = i_vdw_energy + vdw_energy
                                        ELSE IF (lj_cut) THEN
                                                eps = ppvdwp_table2(1,itype,jtype,ibox) ! 4*epsilon
                                                negsigsq = ppvdwp_table2(2,itype,jtype,ibox) ! -(sigma**2)
                                                vdw_energy = 0.0_DP
                                                IF (l_read_svec) THEN
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        vdw_energy = vdw_energy + rterm
                                                                        mol_vdw_energy_p(j) = eps*rterm 
                                                                END DO
                                                        ELSE
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        vdw_energy = vdw_energy + rterm
                                                                END DO
                                                        END IF
                                                ELSE
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        vdw_energy = vdw_energy + rterm
                                                                        mol_vdw_energy_p(j) = eps*rterm 
                                                                END DO
                                                        ELSE
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        negsigbyr2 = negsigsq/rijsq
                                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                                        rterm = rterm + rterm*rterm
                                                                        vdw_energy = vdw_energy + rterm
                                                                END DO
                                                        END IF
                                                END IF
                                                !DO j = 1, vdw_vlen
                                                !        IF (l_read_svec) THEN
                                                !                rijsq = rijsq_svec(j)
                                                !        ELSE
                                                !                rijsq = rijsq_vec(j)
                                                !        END IF
                                                !        negsigbyr2 = negsigsq/rijsq
                                                !        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                !        rterm = rterm + rterm*rterm
                                                !        vdw_energy = vdw_energy + rterm
                                                !        IF (l_pair_store) THEN
                                                !                mol_vdw_energy_p(j) = eps*rterm 
                                                !        END IF
                                                !        !vdw_energy = vdw_energy + eps * rterm
                                                !END DO
                                                i_vdw_energy = i_vdw_energy + vdw_energy*eps
                                        ELSE IF (mie_cut .OR. mie_cut_shift) THEN
                                                epsig_n = ppvdwp_table2(1,itype,jtype,ibox) ! epsilon * mie_coeff * sigma ** n
                                                epsig_m = ppvdwp_table2(2,itype,jtype,ibox) ! epsilon * mie_coeff * sigma ** m
                                                mie_n = ppvdwp_table2(3,itype,jtype,ibox) ! already halved
                                                mie_m = ppvdwp_table2(4,itype,jtype,ibox) ! already halved
                                                IF (mie_cut_shift) THEN
                                                        shift_p1 = ppvdwp_table2(5,itype,jtype,ibox)
                                                        IF (l_pair_store) THEN
                                                                vdw_energy = 0.0_DP
                                                        ELSE
                                                                vdw_energy = -shift_p1*vdw_vlen
                                                        END IF
                                                ELSE
                                                        vdw_energy = 0.0_DP
                                                        shift_p1 = 0.0_DP
                                                END IF
                                                IF (l_read_svec) THEN
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        nrg = shift_p1
                                                                        nrg = epsig_n * rijsq**mie_n - nrg
                                                                        nrg = nrg - epsig_m * rijsq**mie_m
                                                                        mol_vdw_energy_p(j) = nrg
                                                                        vdw_energy = vdw_energy + nrg
                                                                END DO
                                                        ELSE
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        vdw_energy = vdw_energy + epsig_n * rijsq**mie_n - &
                                                                                epsig_m * rijsq**mie_m
                                                                        !vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
                                                                END DO
                                                        END IF
                                                ELSE
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        nrg = shift_p1
                                                                        nrg = epsig_n * rijsq**mie_n - nrg
                                                                        nrg = nrg - epsig_m * rijsq**mie_m
                                                                        mol_vdw_energy_p(j) = nrg
                                                                        vdw_energy = vdw_energy + nrg
                                                                END DO
                                                        ELSE
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        vdw_energy = vdw_energy + epsig_n * rijsq**mie_n - &
                                                                                epsig_m * rijsq**mie_m
                                                                        !vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
                                                                        !vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
                                                                END DO
                                                        END IF
                                                END IF
                                                !!$OMP SIMD PRIVATE(rijsq,nrg) REDUCTION(+:vdw_energy)
                                                !DO j = 1, vdw_vlen
                                                !        IF (l_read_svec) THEN
                                                !                rijsq = rijsq_svec(j)
                                                !        ELSE
                                                !                rijsq = rijsq_vec(j)
                                                !        END IF
                                                !        IF (l_pair_store) THEN
                                                !                nrg = shift_p1
                                                !                nrg = epsig_n * rijsq**mie_n - nrg
                                                !                nrg = nrg - epsig_m * rijsq**mie_m
                                                !                mol_vdw_energy_p(j) = nrg
                                                !                vdw_energy = vdw_energy + nrg
                                                !        ELSE
                                                !                vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
                                                !                vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
                                                !        END IF
                                                !        !vdw_energy = vdw_energy - shift_p1
                                                !END DO
                                                !!$OMP END SIMD
                                                i_vdw_energy = i_vdw_energy + vdw_energy
                                        END IF
                                        IF (l_pair_store) THEN
                                                IF (l_getdrcom .AND. .NOT. l_shortvdw) THEN
                                                        mol_vdw_energy(1:n_all_interact_p) = &
                                                                mol_vdw_energy(1:n_all_interact_p) + &
                                                                mol_vdw_energy_p(1:n_all_interact_p)
                                                        mol_vdw_energy(n_all_interact_p+1:n_interact) = &
                                                                mol_vdw_energy(n_all_interact_p+1:n_interact) + &
                                                                UNPACK(mol_vdw_energy_p(n_all_interact_p+1:vdw_vlen), &
                                                                vdw_interact_vec(n_all_interact_p+1:n_interact), &
                                                                SPREAD(0.0_DP,1,n_may_interact))
                                                ELSE
                                                        mol_vdw_energy(1:n_interact) = &
                                                                mol_vdw_energy(1:n_interact) + &
                                                                UNPACK(mol_vdw_energy_p(1:vdw_vlen), &
                                                                vdw_interact_vec(1:n_interact), &
                                                                SPREAD(0.0_DP,1,n_interact))
                                                END IF
                                        END IF
                                END IF
                                IF (ij_get_qq) THEN
                                        l_read_svec = ij_get_vdw .AND. l_shortcoul
                                        IF (l_charge_ewald) THEN
                                                qq_energy = 0.0_DP
                                                IF (l_read_svec) THEN
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, qq_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        rij = SQRT(rijsq)
                                                                        nrg = ERFC(alpha_ewald(this_box)*rij) / rij
                                                                        mol_qq_energy_p(j) = cfqq*nrg
                                                                        qq_energy = qq_energy + nrg
                                                                END DO
                                                        ELSE
                                                                DO j = 1, qq_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        rij = SQRT(rijsq)
                                                                        qq_energy = qq_energy + ERFC(alpha_ewald(this_box)*rij) / rij
                                                                END DO
                                                        END IF
                                                ELSE
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, qq_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        rij = SQRT(rijsq)
                                                                        nrg = ERFC(alpha_ewald(this_box)*rij) / rij
                                                                        mol_qq_energy_p(j) = cfqq*nrg
                                                                        qq_energy = qq_energy + nrg
                                                                END DO
                                                        ELSE
                                                                !! experimental algorithm.  swap out with commented if it's bad
                                                                !DO j = 1, qq_vlen
                                                                !        alpharsq = alpha_sq*rijsq_vec(j)
                                                                !        nrg = EXP(-alpharsq)*ewald_const
                                                                !        !numer = 1.0_DP
                                                                !        !denom = 3.5_DP
                                                                !        !denom = denom*2.5_DP
                                                                !        !denom = denom + numer*alpharsq
                                                                !        !numer = numer*2.0_DP + denom
                                                                !        denom = 8.75_DP + 9.0_DP*alpharsq + alpharsq*alpharsq
                                                                !        !numer = 13.0_DP + 2.0_DP*alpharsq + denom
                                                                !        numer = 21.75_DP + 11.0_DP*alpharsq + alpharsq*alpharsq
                                                                !        denom = denom*1.5_DP
                                                                !        denom = denom + numer*alpharsq
                                                                !        numer = numer + denom ! even_denom
                                                                !        denom = denom*0.5_DP
                                                                !        denom = denom + numer*alpharsq
                                                                !        qq_energy = qq_energy + (numer/denom)*nrg
                                                                !END DO
                                                                !! experimental algorithm.  swap out with commented if it's bad
                                                                !DO j = 1, qq_vlen
                                                                !        alpharsq = alpha_sq*rijsq_vec(j)
                                                                !        nrg = EXP(-alpharsq)*ewald_const
                                                                !        !numer = 1.0_DP
                                                                !        !denom = 2.5_DP
                                                                !        !numer = 4.5_DP + alpharsq
                                                                !        !denom = 3.75_DP + 1.5_DP*alpharsq !denom*1.5_DP
                                                                !        !denom = denom + numer*alpharsq
                                                                !        denom = 3.75_DP + 6.0_DP*alpharsq + alpharsq*alpharsq
                                                                !        numer = 4.5_DP + alpharsq + denom ! even_denom
                                                                !        denom = denom*0.5_DP
                                                                !        denom = denom + numer*alpharsq
                                                                !        qq_energy = qq_energy + (numer/denom)*nrg
                                                                !END DO
                                                                !! experimental algorithm.  swap out with commented if it's bad
                                                                !DO j = 1, qq_vlen
                                                                !        alpharsq = alpha_sq*rijsq_vec(j)
                                                                !        nrg = EXP(-alpharsq)*ewald_const
                                                                !        numer = 1.0_DP
                                                                !        denom = 2.5_DP
                                                                !        denom = denom + numer*alpharsq
                                                                !        numer = numer*2.0_DP + denom
                                                                !        denom = denom*1.5_DP
                                                                !        denom = denom + numer*alpharsq
                                                                !        numer = numer*1.0_DP + denom ! even_denom
                                                                !        denom = denom*0.5_DP
                                                                !        denom = denom + numer*alpharsq
                                                                !        qq_energy = qq_energy + (numer/denom)*nrg
                                                                !END DO
                                                                !! experimental algorithm.  swap out with commented if it's bad
                                                                !DO j = 1, qq_vlen
                                                                !        alpharsq = alpha_sq*rijsq_vec(j)
                                                                !        nrg = EXP(-alpharsq)*ewald_const
                                                                !        odd_denom = 1.0_DP
                                                                !        odd_numer = 2.5_DP
                                                                !        odd_numer = odd_numer + odd_denom*alpharsq
                                                                !        odd_denom = odd_denom*2.0_DP + odd_numer
                                                                !        odd_numer = odd_numer*1.5_DP
                                                                !        odd_numer = odd_numer + odd_denom*alpharsq
                                                                !        odd_denom = odd_denom*1.0_DP + odd_numer ! even_denom
                                                                !        odd_numer = odd_numer*0.5_DP
                                                                !        odd_numer = odd_numer + odd_denom*alpharsq
                                                                !        qq_energy = qq_energy + (odd_denom/odd_numer)*nrg
                                                                !END DO
                                                                !! experimental algorithm.  swap out with commented if it's bad
                                                                !DO j = 1, qq_vlen
                                                                !        alpharsq = alpha_sq*rijsq_vec(j)
                                                                !        nrg = EXP(-alpharsq)*ewald_const
                                                                !        odd_denom = 1.0_DP
                                                                !        odd_numer = 2.5_DP
                                                                !        odd_numer = odd_numer + odd_denom*alpharsq
                                                                !        even_denom = odd_numer
                                                                !        even_numer = odd_denom*2.0_DP
                                                                !        even_numer = even_numer + even_denom
                                                                !        odd_denom = even_numer
                                                                !        odd_numer = even_denom*1.5_DP
                                                                !        odd_numer = odd_numer + odd_denom*alpharsq
                                                                !        even_denom = odd_numer
                                                                !        even_numer = odd_denom*1.0_DP
                                                                !        even_numer = even_numer + even_denom
                                                                !        odd_denom = even_numer
                                                                !        odd_numer = even_denom*0.5_DP
                                                                !        odd_numer = odd_numer + odd_denom*alpharsq
                                                                !        confrac = odd_denom/odd_numer
                                                                !        nrg = nrg*confrac
                                                                !        qq_energy = qq_energy + nrg
                                                                !END DO
                                                                DO j = 1, qq_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        rij = SQRT(rijsq)
                                                                        qq_energy = qq_energy + ERFC(alpha_ewald(this_box)*rij) / rij
                                                                END DO
                                                        END IF
                                                END IF
                                                !!$OMP SIMD PRIVATE(rijsq,rij,nrg) REDUCTION(+:qq_energy)
                                                !DO j = 1, qq_vlen
                                                !        IF (l_read_svec) THEN
                                                !                rijsq = rijsq_svec(j)
                                                !        ELSE
                                                !                rijsq = rijsq_vec(j)
                                                !        END IF
                                                !        rij = SQRT(rijsq)
                                                !        IF (l_pair_store) THEN
                                                !                nrg = ERFC(alpha_ewald(this_box)*rij) / rij
                                                !                mol_qq_energy_p(j) = cfqq*nrg
                                                !                qq_energy = qq_energy + nrg
                                                !        ELSE
                                                !                qq_energy = qq_energy + ERFC(alpha_ewald(this_box)*rij) / rij
                                                !        END IF
                                                !        !qq_energy = qq_energy + ERFC(alpha_ewald(this_box)*rij) / rij * cfqq
                                                !END DO
                                                !!$OMP END SIMD
                                                i_qq_energy = i_qq_energy + cfqq*qq_energy
                                        ELSE IF (l_charge_dsf) THEN
                                                IF (l_pair_store) THEN
                                                        qq_energy = 0.0_DP
                                                ELSE
                                                        qq_energy = dsf_const*qq_vlen
                                                END IF
                                                IF (l_read_svec) THEN
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, qq_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        rij = SQRT(rijsq)
                                                                        nrg = dsf_const
                                                                        nrg = nrg + dsf_factor2(this_box)*rij + &
                                                                                ERFC(alpha_dsf(this_box)*rij) / rij
                                                                        qq_energy = qq_energy + nrg
                                                                        nrg = nrg*cfqq
                                                                        mol_qq_energy_p(j) = nrg
                                                                END DO
                                                        ELSE
                                                                DO j = 1, qq_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        rij = SQRT(rijsq)
                                                                        qq_energy = qq_energy + dsf_factor2(this_box) * rij + &
                                                                                ERFC(alpha_dsf(this_box)*rij) / rij
                                                                END DO
                                                        END IF
                                                ELSE
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, qq_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        rij = SQRT(rijsq)
                                                                        nrg = dsf_const
                                                                        nrg = nrg + dsf_factor2(this_box)*rij + &
                                                                                ERFC(alpha_dsf(this_box)*rij) / rij
                                                                        qq_energy = qq_energy + nrg
                                                                        nrg = nrg*cfqq
                                                                        mol_qq_energy_p(j) = nrg
                                                                END DO
                                                        ELSE
                                                                DO j = 1, qq_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        rij = SQRT(rijsq)
                                                                        qq_energy = qq_energy + dsf_factor2(this_box) * rij + &
                                                                                ERFC(alpha_dsf(this_box)*rij) / rij
                                                                END DO
                                                        END IF
                                                END IF
                                                !!$OMP SIMD PRIVATE(rijsq,rij,nrg) REDUCTION(+:qq_energy)
                                                !DO j = 1, qq_vlen
                                                !        IF (l_read_svec) THEN
                                                !                rijsq = rijsq_svec(j)
                                                !        ELSE
                                                !                rijsq = rijsq_vec(j)
                                                !        END IF
                                                !        rij = SQRT(rijsq)
                                                !        IF (l_pair_store) THEN
                                                !                nrg = dsf_const
                                                !                nrg = nrg + dsf_factor2(this_box)*rij + &
                                                !                        ERFC(alpha_dsf(this_box)*rij) / rij
                                                !                qq_energy = qq_energy + nrg
                                                !                nrg = nrg*cfqq
                                                !                mol_qq_energy_p(j) = nrg
                                                !        ELSE
                                                !                qq_energy = qq_energy + dsf_factor2(this_box) * rij + &
                                                !                        ERFC(alpha_dsf(this_box)*rij) / rij
                                                !        END IF
                                                !        !qq_energy = qq_energy + (dsf_factor2(this_box) * (rij - rcut_coul(this_box)) - &
                                                !        !        dsf_factor1(this_box) + &
                                                !        !        ERFC(alpha_dsf(this_box)*rij) / rij) * cfqq
                                                !END DO
                                                !!$OMP END SIMD
                                                i_qq_energy = i_qq_energy + cfqq*qq_energy
                                        ELSE IF (l_charge_cut) THEN
                                                qq_energy = 0.0_DP
                                                IF (l_read_svec) THEN
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        nrg = cfqq / SQRT(rijsq)
                                                                        mol_qq_energy_p(j) = nrg
                                                                        qq_energy = qq_energy + nrg
                                                                END DO
                                                        ELSE
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_svec(j)
                                                                        qq_energy = qq_energy + 1.0_DP / SQRT(rijsq)
                                                                END DO
                                                        END IF
                                                ELSE
                                                        IF (l_pair_store) THEN
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        nrg = cfqq / SQRT(rijsq)
                                                                        mol_qq_energy_p(j) = nrg
                                                                        qq_energy = qq_energy + nrg
                                                                END DO
                                                        ELSE
                                                                DO j = 1, vdw_vlen
                                                                        rijsq = rijsq_vec(j)
                                                                        qq_energy = qq_energy + 1.0_DP / SQRT(rijsq)
                                                                END DO
                                                        END IF
                                                END IF
                                                !!$OMP SIMD PRIVATE(rijsq,nrg) REDUCTION(+:qq_energy)
                                                !DO j = 1, qq_vlen
                                                !        IF (l_read_svec) THEN
                                                !                rijsq = rijsq_svec(j)
                                                !        ELSE
                                                !                rijsq = rijsq_vec(j)
                                                !        END IF
                                                !        IF (l_pair_store) THEN
                                                !                nrg = cfqq / SQRT(rijsq)
                                                !                mol_qq_energy_p(j) = nrg
                                                !                qq_energy = qq_energy + nrg
                                                !        ELSE
                                                !                qq_energy = qq_energy + 1.0_DP / SQRT(rijsq)
                                                !        END IF
                                                !        !qq_energy = qq_energy + cfqq / rij
                                                !END DO
                                                !!$OMP END SIMD
                                                IF (l_pair_store) THEN
                                                        i_qq_energy = i_qq_energy + qq_energy
                                                ELSE
                                                        i_qq_energy = i_qq_energy + cfqq*qq_energy
                                                END IF
                                        END IF
                                        IF (l_pair_store) THEN
                                                IF (l_getdrcom .AND. .NOT. l_shortcoul) THEN
                                                        mol_qq_energy(1:n_all_interact_p) = &
                                                                mol_qq_energy(1:n_all_interact_p) + &
                                                                mol_qq_energy_p(1:n_all_interact_p)
                                                        mol_qq_energy(n_all_interact_p+1:n_interact) = &
                                                                mol_qq_energy(n_all_interact_p+1:n_interact) + &
                                                                UNPACK(mol_qq_energy_p(n_all_interact_p+1:qq_vlen), &
                                                                qq_interact_vec(n_all_interact_p+1:n_interact), &
                                                                SPREAD(0.0_DP,1,n_may_interact))
                                                ELSE
                                                        mol_qq_energy(1:n_interact) = &
                                                                mol_qq_energy(1:n_interact) + &
                                                                UNPACK(mol_qq_energy_p(1:qq_vlen), &
                                                                qq_interact_vec(1:n_interact), &
                                                                SPREAD(0.0_DP,1,n_interact))
                                                END IF
                                        END IF
                                END IF
                        END DO
                    END DO
                    !$OMP END DO
                    IF (l_pair_store) THEN
                            IF (l_getdrcom) THEN
                                    !$OMP WORKSHARE
                                    pair_nrg_vdw(which_interact(1:n_interact_vdw_pairstore)+jsl_base,isl) = mol_vdw_energy(1:n_interact_vdw_pairstore)
                                    pair_nrg_vdw(isl,which_interact(1:n_interact_vdw_pairstore)+jsl_base) = mol_vdw_energy(1:n_interact_vdw_pairstore)
                                    pair_nrg_qq(which_interact(1:n_interact_qq_pairstore)+jsl_base,isl) = mol_qq_energy(1:n_interact_qq_pairstore)
                                    pair_nrg_qq(isl,which_interact(1:n_interact_qq_pairstore)+jsl_base) = mol_qq_energy(1:n_interact_qq_pairstore)
                                    !$OMP END WORKSHARE
                            ELSE IF (open_mc_flag) THEN
                                    !$OMP WORKSHARE
                                    pair_nrg_vdw(live_locates(1:n_interact_vdw_pairstore)+jsl_base,isl) = mol_vdw_energy(1:n_interact_vdw_pairstore)
                                    pair_nrg_vdw(isl,live_locates(1:n_interact_vdw_pairstore)+jsl_base) = mol_vdw_energy(1:n_interact_vdw_pairstore)
                                    pair_nrg_qq(live_locates(1:n_interact_qq_pairstore)+jsl_base,isl) = mol_qq_energy(1:n_interact_qq_pairstore)
                                    pair_nrg_qq(isl,live_locates(1:n_interact_qq_pairstore)+jsl_base) = mol_qq_energy(1:n_interact_qq_pairstore)
                                    !$OMP END WORKSHARE
                            ELSE
                                    !$OMP WORKSHARE
                                    pair_nrg_vdw(jsl_base+1:jsl_base+n_interact_vdw_pairstore,isl) = mol_vdw_energy(1:n_interact_vdw_pairstore)
                                    pair_nrg_qq(jsl_base+1:jsl_base+n_interact_qq_pairstore,isl) = mol_qq_energy(1:n_interact_qq_pairstore)
                                    pair_nrg_vdw(isl,jsl_base+1:jsl_base+n_interact_vdw_pairstore) = mol_vdw_energy(1:n_interact_vdw_pairstore)
                                    pair_nrg_qq(isl,jsl_base+1:jsl_base+n_interact_qq_pairstore) = mol_qq_energy(1:n_interact_qq_pairstore)
                                    !$OMP END WORKSHARE
                            END IF
                    END IF
            ELSE
                    vdw_p1(1:natompairs) = vdw_p(1:natompairs,1)
                    vdw_p2(1:natompairs) = vdw_p(1:natompairs,2)
                    IF (n_vdw_p > 2) vdw_p3(1:natompairs) = vdw_p(1:natompairs,3)
                    IF (n_vdw_p > 3) vdw_p4(1:natompairs) = vdw_p(1:natompairs,4)
                    IF (n_vdw_p > 4) vdw_p5(1:natompairs) = vdw_p(1:natompairs,5)
                    IF (l_debug_print) WRITE(*,*) "n_interact = ", n_interact
                    !$OMP DO SCHEDULE(DYNAMIC)
                    DO j_interact = 1, n_interact !n_all_interact
                        !DIR$ ASSUME (MOD(atompairdim,8) .EQ. 0)
                        !DIR$ ASSUME (MOD(mol_dim,8) .EQ. 0)
                        !DIR$ ASSUME (MOD(SIZE(vdw_p,1),8) .EQ. 0)
                        !DIR$ ASSUME (MOD(SIZE(jrp_perm,1),8) .EQ. 0)
                        IF (overlap .OR. i_overlap) CYCLE
                        ij_overlap = .FALSE.
                        vdw_energy = 0.0_DP
                        qq_energy = 0.0_DP
                        jm = which_interact(j_interact) !which_all_interact(j_interact)
                        IF (l_ortho) THEN
                                IF (l_order_ij) THEN
                                        DO j = 1, jnatoms
                                                jsxp = atom_list(j,jm,js)%rp(1)
                                                jsyp = atom_list(j,jm,js)%rp(2)
                                                jszp = atom_list(j,jm,js)%rp(3)
                                                k = (j-1) * n_i_exist
                                                !$OMP SIMD
                                                DO i = 1, n_i_exist
                                                        jrxp(i+k) = jsxp
                                                        jryp(i+k) = jsyp
                                                        jrzp(i+k) = jszp
                                                END DO
                                                !$OMP END SIMD
                                        END DO
                                ELSE
                                        !$OMP SIMD PRIVATE(jsxp,jsyp,jszp)
                                        DO j = 1, jnatoms
                                                jsxp = atom_list(j,jm,js)%rp(1)
                                                jsyp = atom_list(j,jm,js)%rp(2)
                                                jszp = atom_list(j,jm,js)%rp(3)
                                                jrxp(j) = jsxp
                                                jryp(j) = jsyp
                                                jrzp(j) = jszp
                                                DO i = 1, n_i_exist - 1
                                                        k = i*jnatoms
                                                        jrxp(j+k) = jsxp
                                                        jryp(j+k) = jsyp
                                                        jrzp(j+k) = jszp
                                                END DO
                                        END DO
                                        !$OMP END SIMD
                                END IF
                        ELSE
                                !$OMP SIMD PRIVATE(jrp,jsxp,jsyp,jszp)
                                DO j = 1, jnatoms
                                        jrp = atom_list(j,jm,js)%rp(1)
                                        jsxp = inv_h11*jrp
                                        jsyp = inv_h21*jrp
                                        jszp = inv_h31*jrp
                                        jrp = atom_list(j,jm,js)%rp(2)
                                        jsxp = jsxp + inv_h12*jrp
                                        jsyp = jsyp + inv_h22*jrp
                                        jszp = jszp + inv_h32*jrp
                                        jrp = atom_list(j,jm,js)%rp(3)
                                        jsxp = jsxp + inv_h13*jrp
                                        jsyp = jsyp + inv_h23*jrp
                                        jszp = jszp + inv_h33*jrp
                                        jrxp(j) = jsxp
                                        jryp(j) = jsyp
                                        jrzp(j) = jszp
                                        DO i = 1, n_i_exist-1
                                                k = i*jnatoms
                                                jrxp(i+k) = jsxp
                                                jryp(i+k) = jsyp
                                                jrzp(i+k) = jszp
                                        END DO
                                END DO
                                !$OMP END SIMD
                        END IF
                        !jrxp = RESHAPE(SPREAD(jrxp_vec,2,n_i_exist), (/ natompairs /))
                        !jryp = RESHAPE(SPREAD(jryp_vec,2,n_i_exist), (/ natompairs /))
                        !jrzp = RESHAPE(SPREAD(jrzp_vec,2,n_i_exist), (/ natompairs /))
                        IF (l_ortho) THEN
                                DO ji = 1, natompairs
                                        dxp = ABS(jrxp(ji) - irxp(ji))
                                        dyp = ABS(jryp(ji) - iryp(ji))
                                        dzp = ABS(jrzp(ji) - irzp(ji))
                                        IF (dxp > hxl) dxp = dxp - xl
                                        IF (dyp > hyl) dyp = dyp - yl
                                        IF (dzp > hzl) dzp = dzp - zl
                                        dxp = dxp*dxp
                                        dxp = dxp + dyp*dyp
                                        dxp = dxp + dzp*dzp
                                        rijsq_vec(ji) = dxp
                                        vdw_interact_vec(ji) = dxp < this_vdw_rcutsq
                                        IF (need_sqrt) rij_vec(ji) = SQRT(dxp)
                                        i_overlap = i_overlap .OR. dxp < rcut_lowsq
                                END DO
                        ELSE
                                DO ji = 1, natompairs
                                        dsp = jrxp(ji) - irxp(ji)
                                        IF (dsp > 0.5_DP) THEN
                                                dsp = dsp - 1.0_DP
                                        ELSE IF (dsp < -0.5_DP) THEN
                                                dsp = dsp + 1.0_DP
                                        END IF
                                        dxp = h11*dsp
                                        dyp = h21*dsp
                                        dzp = h31*dsp
                                        dsp = jryp(ji) - iryp(ji)
                                        IF (dsp > 0.5_DP) THEN
                                                dsp = dsp - 1.0_DP
                                        ELSE IF (dsp < -0.5_DP) THEN
                                                dsp = dsp + 1.0_DP
                                        END IF
                                        dxp = dxp + h12*dsp
                                        dyp = dyp + h22*dsp
                                        dzp = dzp + h32*dsp
                                        dsp = jrzp(ji) - irzp(ji)
                                        IF (dsp > 0.5_DP) THEN
                                                dsp = dsp - 1.0_DP
                                        ELSE IF (dsp < -0.5_DP) THEN
                                                dsp = dsp + 1.0_DP
                                        END IF
                                        dxp = dxp + h13*dsp
                                        dyp = dyp + h23*dsp
                                        dzp = dzp + h33*dsp
                                        dxp = dxp*dxp
                                        dxp = dxp + dyp*dyp
                                        dxp = dxp + dzp*dzp
                                        rijsq_vec(ji) = dxp
                                        !vdw_interact_vec(ji) = dxp < this_vdw_rcutsq
                                        i_overlap = i_overlap .OR. dxp < rcut_lowsq
                                END DO
                        END IF
                        IF (l_debug_print) THEN
                                WRITE(*,*) "jm = ", jm
                                WRITE(*,*) "i_overlap = ", i_overlap
                                WRITE(*,*) "rijsq_vec = "
                                WRITE(*,*) rijsq_vec(1:natompairs)
                        END IF
                        IF (i_overlap) THEN
                                overlap = .TRUE.
                                CYCLE
                        END IF
                        !ji_end = 1
                        !ji_base_2 = 0
                        !DO WHILE (ji_end .LE. natompairs)
                        !        ji_start = natompairs+1
                        !        ji_base = ji_end - 1
                        !        DO WHILE (ji_base .LE. natompairs - 32)
                        !                ji_start_byte = 33
                        !                DO ji_byte = 1, 32
                        !                        IF (vdw_interact_vec(ji_byte+ji_base)) ji_start_byte = MIN(ji_start_byte,ji_byte)
                        !                END DO
                        !                IF (ji_start_byte < 33) THEN
                        !                        ji_start = ji_base + ji_start_byte
                        !                        EXIT
                        !                END IF
                        !                ji_base = ji_base + 32
                        !        END DO
                        !        IF (ji_start > natompairs) THEN
                        !                !DIR$ ASSUME (natompairs-ji_base<32)
                        !                DO ji = ji_base+1, natompairs
                        !                        IF (vdw_interact_vec(ji)) ji_start = MIN(ji_start,ji)
                        !                END DO
                        !        END IF
                        !        !DO ji = ji_end, natompairs ! ji_end
                        !        !        IF (vdw_interact_vec(ji)) ji_start = MIN(ji_start,ji)
                        !        !END DO
                        !        IF (ji_start > natompairs) EXIT
                        !        ji_end = natompairs+1
                        !        DO ji = ji_start, natompairs
                        !                IF (.NOT. vdw_interact_vec(ji)) ji_end = MIN(ji_end,ji)
                        !                !IF (.NOT. vdw_interact_vec(ji)) EXIT
                        !        END DO
                        !        vlen = ji_end - ji_start
                        !        ji_end_2 = ji_base_2 + vlen
                        !        rijsq_packed(ji_base_2+1:ji_end_2) = rijsq_vec(ji_start:ji_end-1)
                        !        vdw_p_packed(ji_base_2+1:ji_end_2,1:n_vdw_p) = vdw_p(ji_start:ji_end-1,1:n_vdw_p)
                        !        ji_base_2 = ji_end_2
                        !END DO
                        !atompair_vlen = ji_end_2
                        !IF (l_vdw_packed) THEN
                        !        i_selector = 1
                        !ELSE IF (l_vdw_gather) THEN
                        !        i_selector = 2
                        !ELSE IF (l_vdw_strided) THEN
                        !        i_selector = 3
                        !ELSE
                        !        i_selector = 4
                        !END IF
                        !i_selector = nspecies ! remove this, it's just for testing vectorization
                        !l_masked = atompair_vlen < 16 ! remove this, it's just for testing vectorization
                        l_masked = j_interact > n_all_interact .OR. l_shortvdw
                        atompair_vlen = natompairs ! temporary, replace when appropriate
                        IF (lj_charmm) THEN
                                IF (l_masked) THEN
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                IF (rijsq < this_vdw_rcutsq) THEN
                                                        eps = vdw_p1(ji) ! epsilon
                                                        sigsq = vdw_p2(ji) ! sigma**2
                                                        sigbyr2 = sigsq/rijsq ! sigma was already squared
                                                        sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                                        sigbyr12 = sigbyr6*sigbyr6
                                                        sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
                                                        vdw_energy = vdw_energy + eps*sigbyr12
                                                END IF
                                        END DO
                                ELSE
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                eps = vdw_p1(ji) ! epsilon
                                                sigsq = vdw_p2(ji) ! sigma**2
                                                sigbyr2 = sigsq/rijsq ! sigma was already squared
                                                sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                                sigbyr12 = sigbyr6*sigbyr6
                                                sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
                                                vdw_energy = vdw_energy + eps*sigbyr12
                                        END DO
                                END IF
                                !DO ji = 1, natompairs
                                !        rijsq = rijsq_vec(ji)
                                !        IF (l_masked) THEN
                                !                l_interact = rijsq < this_vdw_rcutsq
                                !        ELSE
                                !                l_interact = .TRUE.
                                !        END IF
                                !        IF (l_interact) THEN
                                !                eps = vdw_p1(ji) ! epsilon
                                !                sigsq = vdw_p2(ji) ! sigma**2
                                !                sigbyr2 = sigsq/rijsq ! sigma was already squared
                                !                sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                !                sigbyr12 = sigbyr6*sigbyr6
                                !                sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
                                !                vdw_energy = vdw_energy + eps*sigbyr12
                                !        END IF
                                !        !IF (rijsq < this_vdw_rcutsq .OR. .NOT. l_masked) THEN
                                !        !        vdw_energy = vdw_energy + eps*sigbyr12
                                !        !END IF
                                !        !ELSE
                                !        !        l_interact = .TRUE.
                                !        !END IF
                                !END DO
                                !DO ji = 1, natompairs
                                !        rijsq = rijsq_vec(ji)
                                !        eps = vdw_p1(ji) ! epsilon
                                !        sigsq = vdw_p2(ji) ! sigma**2
                                !        sigbyr2 = sigsq/rijsq ! sigma was already squared
                                !        sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                !        sigbyr12 = sigbyr6*sigbyr6
                                !        sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
                                !        vdw_energy = vdw_energy + eps*sigbyr12
                                !        !IF (rijsq < this_vdw_rcutsq .OR. .NOT. l_masked) THEN
                                !        !        vdw_energy = vdw_energy + eps*sigbyr12
                                !        !END IF
                                !        !ELSE
                                !        !        l_interact = .TRUE.
                                !        !END IF
                                !END DO
                                !DO ji = 1, atompair_vlen !natompairs
                                !        SELECT CASE (i_selector)
                                !        CASE (1)
                                !                rijsq = rijsq_packed(ji)
                                !                eps = vdw_p_packed(ji,1) ! epsilon
                                !                sigsq = vdw_p_packed(ji,2) ! sigma**2
                                !        CASE (2)
                                !                rijsq = rijsq_vec(which_have_vdw(ji))
                                !                eps = vdw_p_nonzero(ji,1) ! epsilon
                                !                sigsq = vdw_p_nonzero(ji,2) ! sigma**2
                                !        CASE (3)
                                !                rijsq = rijsq_vec(ji_base+ji_stride*ji)
                                !                eps = vdw_p_nonzero(ji,1) ! epsilon
                                !                sigsq = vdw_p_nonzero(ji,2) ! sigma**2
                                !        CASE DEFAULT
                                !                rijsq = rijsq_vec(ji)
                                !                eps = vdw_p(ji,1) ! epsilon
                                !                sigsq = vdw_p(ji,2) ! sigma**2
                                !        END SELECT
                                !        sigbyr2 = sigsq/rijsq ! sigma was already squared
                                !        sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                !        sigbyr12 = sigbyr6*sigbyr6
                                !        sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
                                !        IF (l_masked) THEN
                                !                vdw_energy = MERGE(vdw_energy + eps*sigbyr12, &
                                !                       vdw_energy, rijsq < this_vdw_rcutsq)
                                !        ELSE
                                !                vdw_energy = vdw_energy + eps*sigbyr12
                                !        END IF
                                !END DO
                                !IF (l_masked) THEN
                                !        DO ji = 1, atompair_vlen !natompairs
                                !                SELECT CASE (i_selector)
                                !                CASE (1)
                                !                        rijsq = rijsq_packed(ji)
                                !                        eps = vdw_p_packed(ji,1) ! epsilon
                                !                        sigsq = vdw_p_packed(ji,2) ! sigma**2
                                !                CASE (2)
                                !                        rijsq = rijsq_vec(which_have_vdw(ji))
                                !                        eps = vdw_p_nonzero(ji,1) ! epsilon
                                !                        sigsq = vdw_p_nonzero(ji,2) ! sigma**2
                                !                CASE (3)
                                !                        rijsq = rijsq_vec(ji_base+ji_stride*ji)
                                !                        eps = vdw_p_nonzero(ji,1) ! epsilon
                                !                        sigsq = vdw_p_nonzero(ji,2) ! sigma**2
                                !                CASE DEFAULT
                                !                        rijsq = rijsq_vec(ji)
                                !                        eps = vdw_p(ji,1) ! epsilon
                                !                        sigsq = vdw_p(ji,2) ! sigma**2
                                !                END SELECT
                                !                sigbyr2 = sigsq/rijsq ! sigma was already squared
                                !                sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                !                sigbyr12 = sigbyr6*sigbyr6
                                !                sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
                                !                vdw_energy = MERGE(vdw_energy + eps*sigbyr12, &
                                !                       vdw_energy, rijsq < this_vdw_rcutsq)
                                !        END DO
                                !ELSE
                                !        DO ji = 1, atompair_vlen !natompairs
                                !                SELECT CASE (i_selector)
                                !                CASE (1)
                                !                        rijsq = rijsq_packed(ji)
                                !                        eps = vdw_p_packed(ji,1) ! epsilon
                                !                        sigsq = vdw_p_packed(ji,2) ! sigma**2
                                !                CASE (2)
                                !                        rijsq = rijsq_vec(which_have_vdw(ji))
                                !                        eps = vdw_p_nonzero(ji,1) ! epsilon
                                !                        sigsq = vdw_p_nonzero(ji,2) ! sigma**2
                                !                CASE (3)
                                !                        rijsq = rijsq_vec(ji_base+ji_stride*ji)
                                !                        eps = vdw_p_nonzero(ji,1) ! epsilon
                                !                        sigsq = vdw_p_nonzero(ji,2) ! sigma**2
                                !                CASE DEFAULT
                                !                        rijsq = rijsq_vec(ji)
                                !                        eps = vdw_p(ji,1) ! epsilon
                                !                        sigsq = vdw_p(ji,2) ! sigma**2
                                !                END SELECT
                                !                sigbyr2 = sigsq/rijsq ! sigma was already squared
                                !                sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                !                sigbyr12 = sigbyr6*sigbyr6
                                !                sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
                                !                vdw_energy = vdw_energy + eps * sigbyr12
                                !        END DO
                                !END IF
                                !DO ji = 1, atompair_vlen !natompairs
                                !        IF (l_vdw_packed) THEN
                                !                rijsq = rijsq_packed(ji)
                                !                eps = vdw_p_packed(ji,1) ! epsilon
                                !                sigsq = vdw_p_packed(ji,2) ! sigma**2
                                !        ELSE IF (l_vdw_gather) THEN
                                !                rijsq = rijsq_vec(which_have_vdw(ji))
                                !                eps = vdw_p_nonzero(ji,1) ! epsilon
                                !                sigsq = vdw_p_nonzero(ji,2) ! sigma**2
                                !        ELSE IF (l_vdw_strided) THEN
                                !                rijsq = rijsq_vec(ji_base+ji_stride*ji)
                                !                eps = vdw_p_nonzero(ji,1) ! epsilon
                                !                sigsq = vdw_p_nonzero(ji,2) ! sigma**2
                                !        ELSE
                                !                rijsq = rijsq_vec(ji)
                                !                eps = vdw_p(ji,1) ! epsilon
                                !                sigsq = vdw_p(ji,2) ! sigma**2
                                !        END IF
                                !        sigbyr2 = sigsq/rijsq ! sigma was already squared
                                !        sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                !        sigbyr12 = sigbyr6*sigbyr6
                                !        sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
                                !        IF (l_masked) THEN
                                !                vdw_energy = MERGE(vdw_energy + eps*sigbyr12, &
                                !                       vdw_energy, rijsq < this_vdw_rcutsq)
                                !        ELSE
                                !                vdw_energy = vdw_energy + eps * sigbyr12
                                !        END IF
                                !        !CALL Compute_LJ_Charmm(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy)
                                !END DO
                        ELSE IF (lj_cut_shift) THEN
                                IF (l_masked) THEN
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                IF (rijsq < this_vdw_rcutsq) THEN
                                                        eps = vdw_p1(ji) ! 4*epsilon
                                                        negsigsq = vdw_p2(ji) ! -(sigma**2)
                                                        shift_p1 = vdw_p3(ji)
                                                        negsigbyr2 = negsigsq/rijsq
                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                        rterm = rterm + rterm*rterm
                                                        vdw_energy = vdw_energy + eps * rterm - shift_p1
                                                END IF
                                        END DO
                                ELSE
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                eps = vdw_p1(ji) ! 4*epsilon
                                                negsigsq = vdw_p2(ji) ! -(sigma**2)
                                                shift_p1 = vdw_p3(ji)
                                                negsigbyr2 = negsigsq/rijsq
                                                rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                rterm = rterm + rterm*rterm
                                                vdw_energy = vdw_energy + eps * rterm - shift_p1
                                        END DO
                                END IF
                                !DO ji = 1, natompairs
                                !        rijsq = rijsq_packed(ji)
                                !        eps = vdw_p1(ji) ! 4*epsilon
                                !        negsigsq = vdw_p2(ji) ! -(sigma**2)
                                !        shift_p1 = vdw_p3(ji)
                                !        negsigbyr2 = negsigsq/rijsq
                                !        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                !        rterm = rterm + rterm*rterm
                                !        vdw_energy = vdw_energy + eps * rterm - shift_p1
                                !        !CALL Compute_LJ_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy)
                                !END DO
                        ELSE IF (lj_cut_switch) THEN
                                IF (l_masked) THEN
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                IF (rijsq < this_vdw_rcutsq) THEN
                                                        eps = vdw_p1(ji) ! 4*epsilon
                                                        negsigsq = vdw_p2(ji) ! -(sigma**2)
                                                        negsigbyr2 = negsigsq/rijsq
                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                        rterm = rterm + rterm*rterm
                                                        roffsq_rijsq = roff_switch_sq(this_box) - rijsq
                                                        roffsq_rijsq = roffsq_rijsq * roffsq_rijsq * switch_factor1(this_box) * &
                                                                (switch_factor2(this_box)+2.0_DP*rijsq)*eps
                                                        eps = MERGE(roffsq_rijsq, eps, rijsq .GE. ron_switch_sq(this_box))
                                                        vdw_energy = vdw_energy + eps*rterm
                                                END IF
                                        END DO
                                ELSE
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                eps = vdw_p1(ji) ! 4*epsilon
                                                negsigsq = vdw_p2(ji) ! -(sigma**2)
                                                negsigbyr2 = negsigsq/rijsq
                                                rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                rterm = rterm + rterm*rterm
                                                roffsq_rijsq = roff_switch_sq(this_box) - rijsq
                                                roffsq_rijsq = roffsq_rijsq * roffsq_rijsq * switch_factor1(this_box) * &
                                                        (switch_factor2(this_box)+2.0_DP*rijsq)*eps
                                                eps = MERGE(roffsq_rijsq, eps, rijsq .GE. ron_switch_sq(this_box))
                                                vdw_energy = vdw_energy + eps*rterm
                                        END DO
                                END IF
                                !DO ji = 1, natompairs
                                !        rijsq = rijsq_packed(ji)
                                !        eps = vdw_p1(ji) ! 4*epsilon
                                !        negsigsq = vdw_p2(ji) ! -(sigma**2)
                                !        negsigbyr2 = negsigsq/rijsq
                                !        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                !        rterm = rterm + rterm*rterm
                                !        roffsq_rijsq = roff_switch_sq(this_box) - rijsq
                                !        roffsq_rijsq = roffsq_rijsq * roffsq_rijsq * switch_factor1(this_box) * &
                                !                (switch_factor2(this_box)+2.0_DP*rijsq)*eps
                                !        eps = MERGE(roffsq_rijsq, eps, rijsq .GE. ron_switch_sq(this_box))
                                !        vdw_energy = vdw_energy + eps*rterm
                                !        !CALL Compute_LJ_Cut_Switch(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy,this_box)
                                !END DO
                        ELSE IF (lj_cut_shift_force) THEN
                                IF (l_masked) THEN
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                IF (rijsq < this_vdw_rcutsq) THEN
                                                        eps = vdw_p1(ji) ! 4*epsilon
                                                        negsigsq = vdw_p2(ji) ! -(sigma**2)
                                                        shift_p1 = vdw_p3(ji)
                                                        shift_p2 = vdw_p4(ji)
                                                        negsigbyr2 = negsigsq/rijsq
                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                        rterm = rterm + rterm*rterm
                                                        vdw_energy = vdw_energy + eps * rterm - shift_p1 - &
                                                                (rcut_vdw(this_box)-SQRT(rijsq))*shift_p2
                                                END IF
                                        END DO
                                ELSE
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                eps = vdw_p1(ji) ! 4*epsilon
                                                negsigsq = vdw_p2(ji) ! -(sigma**2)
                                                shift_p1 = vdw_p3(ji)
                                                shift_p2 = vdw_p4(ji)
                                                negsigbyr2 = negsigsq/rijsq
                                                rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                rterm = rterm + rterm*rterm
                                                vdw_energy = vdw_energy + eps * rterm - shift_p1 - &
                                                        (rcut_vdw(this_box)-SQRT(rijsq))*shift_p2
                                        END DO
                                END IF
                                !DO ji = 1, natompairs
                                !        rijsq = rijsq_packed(ji)
                                !        eps = vdw_p1(ji) ! 4*epsilon
                                !        negsigsq = vdw_p2(ji) ! -(sigma**2)
                                !        shift_p1 = vdw_p3(ji)
                                !        shift_p2 = vdw_p4(ji)
                                !        negsigbyr2 = negsigsq/rijsq
                                !        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                !        rterm = rterm + rterm*rterm
                                !        vdw_energy = vdw_energy + eps * rterm - shift_p1 - &
                                !                (rcut_vdw(this_box)-SQRT(rijsq))*shift_p2
                                !        !CALL Compute_LJ_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy)
                                !END DO
                        ELSE IF (lj_cut) THEN
                                IF (l_masked) THEN
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                IF (rijsq < this_vdw_rcutsq) THEN
                                                        eps = vdw_p1(ji) ! 4*epsilon
                                                        negsigsq = vdw_p2(ji) ! -(sigma**2)
                                                        l_interact = rijsq < this_vdw_rcutsq
                                                        negsigbyr2 = negsigsq/rijsq
                                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                        rterm = rterm + rterm*rterm
                                                        vdw_energy = vdw_energy + eps * rterm
                                                END IF
                                        END DO
                                ELSE
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                eps = vdw_p1(ji) ! 4*epsilon
                                                negsigsq = vdw_p2(ji) ! -(sigma**2)
                                                l_interact = rijsq < this_vdw_rcutsq
                                                negsigbyr2 = negsigsq/rijsq
                                                rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                                rterm = rterm + rterm*rterm
                                                vdw_energy = vdw_energy + eps * rterm
                                        END DO
                                END IF
                                !DO ji = 1, natompairs
                                !        rijsq = rijsq_packed(ji)
                                !        eps = vdw_p1(ji) ! 4*epsilon
                                !        negsigsq = vdw_p2(ji) ! -(sigma**2)
                                !        l_interact = rijsq < this_vdw_rcutsq
                                !        negsigbyr2 = negsigsq/rijsq
                                !        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                !        rterm = rterm + rterm*rterm
                                !        vdw_energy = vdw_energy + eps * rterm
                                !        !CALL Compute_LJ_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy)
                                !END DO
                        ELSE IF (mie_cut_shift) THEN
                                IF (l_masked) THEN
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                IF (rijsq < this_vdw_rcutsq) THEN
                                                        epsig_n = vdw_p1(ji) ! epsilon * mie_coeff * sigma ** n
                                                        epsig_m = vdw_p2(ji) ! epsilon * mie_coeff * sigma ** m
                                                        mie_n = vdw_p3(ji) ! already halved
                                                        mie_m = vdw_p4(ji) ! already halved
                                                        shift_p1 = vdw_p5(ji)
                                                        vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
                                                        vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
                                                        vdw_energy = vdw_energy - shift_p1
                                                END IF
                                        END DO
                                ELSE
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                epsig_n = vdw_p1(ji) ! epsilon * mie_coeff * sigma ** n
                                                epsig_m = vdw_p2(ji) ! epsilon * mie_coeff * sigma ** m
                                                mie_n = vdw_p3(ji) ! already halved
                                                mie_m = vdw_p4(ji) ! already halved
                                                shift_p1 = vdw_p5(ji)
                                                vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
                                                vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
                                                vdw_energy = vdw_energy - shift_p1
                                        END DO
                                END IF
                                !DO ji = 1, natompairs
                                !        rijsq = rijsq_packed(ji)
                                !        epsig_n = vdw_p1(ji) ! epsilon * mie_coeff * sigma ** n
                                !        epsig_m = vdw_p2(ji) ! epsilon * mie_coeff * sigma ** m
                                !        mie_n = vdw_p3(ji) ! already halved
                                !        mie_m = vdw_p4(ji) ! already halved
                                !        shift_p1 = vdw_p5(ji)
                                !        vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
                                !        vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
                                !        vdw_energy = vdw_energy - shift_p1
                                !        !CALL Compute_Mie_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2), &
                                !        !        vdw_p(ji,3),vdw_p(ji,4),vdw_energy)
                                !END DO
                        ELSE IF (mie_cut) THEN
                                IF (l_masked) THEN
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                IF (rijsq < this_vdw_rcutsq) THEN
                                                        epsig_n = vdw_p1(ji) ! epsilon * mie_coeff * sigma ** n
                                                        epsig_m = vdw_p2(ji) ! epsilon * mie_coeff * sigma ** m
                                                        mie_n = vdw_p3(ji) ! already halved
                                                        mie_m = vdw_p4(ji) ! already halved
                                                        vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
                                                        vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
                                                END IF
                                        END DO
                                ELSE
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                epsig_n = vdw_p1(ji) ! epsilon * mie_coeff * sigma ** n
                                                epsig_m = vdw_p2(ji) ! epsilon * mie_coeff * sigma ** m
                                                mie_n = vdw_p3(ji) ! already halved
                                                mie_m = vdw_p4(ji) ! already halved
                                                vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
                                                vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
                                        END DO
                                END IF
                                !DO ji = 1, natompairs
                                !        rijsq = rijsq_packed(ji)
                                !        epsig_n = vdw_p1(ji) ! epsilon * mie_coeff * sigma ** n
                                !        epsig_m = vdw_p2(ji) ! epsilon * mie_coeff * sigma ** m
                                !        mie_n = vdw_p3(ji) ! already halved
                                !        mie_m = vdw_p4(ji) ! already halved
                                !        vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
                                !        vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
                                !        !CALL Compute_Mie_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2), &
                                !        !        vdw_p(ji,3),vdw_p(ji,4),vdw_energy)
                                !END DO
                        END IF
                        l_masked = j_interact > n_all_interact .OR. l_shortcoul
                        IF (l_charge_ewald) THEN
                                IF (l_masked) THEN
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                IF (rijsq < this_coul_rcutsq) THEN
                                                        cfqq = cfqq_vec(ji)
                                                        rij = SQRT(rijsq)
                                                        qq_energy = qq_energy + ERFC(alpha_ewald(this_box)*rij) / rij * cfqq
                                                END IF
                                        END DO
                                ELSE
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                cfqq = cfqq_vec(ji)
                                                rij = SQRT(rijsq)
                                                qq_energy = qq_energy + ERFC(alpha_ewald(this_box)*rij) / rij * cfqq
                                        END DO
                                END IF
                                !DO ji = 1, natompairs
                                !        cfqq = cfqq_vec(ji)
                                !        rij = SQRT(rijsq)
                                !        qq_energy = qq_energy + ERFC(alpha_ewald(this_box)*rij) / rij * cfqq
                                !END DO
                        ELSE IF (l_charge_dsf) THEN
                                IF (l_masked) THEN
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                IF (rijsq < this_coul_rcutsq) THEN
                                                        cfqq = cfqq_vec(ji)
                                                        rij = SQRT(rijsq)
                                                        nrg = ERFC(alpha_dsf(this_box)*rij)/rij + dsf_const
                                                        nrg = nrg + rij*dsf_factor2(this_box)
                                                        qq_energy = qq_energy + nrg*cfqq
                                                END IF
                                        END DO
                                ELSE
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                cfqq = cfqq_vec(ji)
                                                rij = SQRT(rijsq)
                                                nrg = ERFC(alpha_dsf(this_box)*rij)/rij + dsf_const
                                                nrg = nrg + rij*dsf_factor2(this_box)
                                                qq_energy = qq_energy + nrg*cfqq
                                        END DO
                                END IF
                                !DO ji = 1, natompairs
                                !        cfqq = cfqq_vec(ji)
                                !        rij = SQRT(rijsq)
                                !        nrg = ERFC(alpha_dsf(this_box)*rij)/rij + dsf_const
                                !        nrg = nrg + rij*dsf_factor2(this_box)
                                !        qq_energy = qq_energy + nrg*cfqq
                                !        !qq_energy = qq_energy + (dsf_factor2(this_box) * (rij - rcut_coul(this_box)) - &
                                !        !        dsf_factor1(this_box) + &
                                !        !        ERFC(alpha_dsf(this_box)*rij) / rij) * cfqq
                                !END DO
                        ELSE IF (l_charge_cut) THEN
                                IF (l_masked) THEN
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                IF (rijsq < this_coul_rcutsq) THEN
                                                        cfqq = cfqq_vec(ji)
                                                        qq_energy = qq_energy + cfqq / SQRT(rijsq)
                                                END IF
                                        END DO
                                ELSE
                                        DO ji = 1, atompair_vlen
                                                rijsq = rijsq_vec(ji)
                                                cfqq = cfqq_vec(ji)
                                                qq_energy = qq_energy + cfqq / SQRT(rijsq)
                                        END DO
                                END IF
                                !DO ji = 1, natompairs
                                !        cfqq = cfqq_vec(ji)
                                !        qq_energy = qq_energy + cfqq / SQRT(rijsq)
                                !END DO
                        END IF
                        !i_overlap = i_overlap .OR. ij_overlap
                        IF (l_pair_store) THEN
                                jsl = jsl_base + jm
                                pair_nrg_vdw(isl,jsl) = vdw_energy
                                pair_nrg_vdw(jsl,isl) = vdw_energy
                                pair_nrg_qq(isl,jsl) = qq_energy
                                pair_nrg_qq(jsl,isl) = qq_energy
                        END IF
                        i_vdw_energy = i_vdw_energy + vdw_energy
                        i_qq_energy = i_qq_energy + qq_energy
                        !IF (i_overlap) THEN
                        !        !$OMP CRITICAL
                        !        WRITE(*,*)
                        !        WRITE(*,*)
                        !        WRITE(*,*) i_overlap, l_order_ij
                        !        WRITE(*,*) im, jm, j_interact, i_mcstep
                        !        WRITE(*,*) vdw_energy, qq_energy
                        !        WRITE(*,*) rijsq_vec(1:natompairs)
                        !        WRITE(*,*)
                        !        WRITE(*,*) irxp(1:natompairs)
                        !        WRITE(*,*) jrxp(1:natompairs)
                        !        WRITE(*,*)
                        !        WRITE(*,*) iryp(1:natompairs)
                        !        WRITE(*,*) jryp(1:natompairs)
                        !        WRITE(*,*)
                        !        WRITE(*,*) irzp(1:natompairs)
                        !        WRITE(*,*) jrzp(1:natompairs)
                        !        !$OMP END CRITICAL
                        !END IF
                    END DO
                    !$OMP END DO
            END IF

            IF (overlap) EXIT


!            !$OMP DO SCHEDULE(STATIC)
!            DO j_interact = 1, n_interact
!                IF (overlap .OR. i_overlap) CYCLE
!                ij_overlap = .FALSE.
!                vdw_energy = 0.0_DP
!                qq_energy = 0.0_DP
!                jm = which_interact(j_interact)
!                IF (l_ortho) THEN
!                        IF (l_order_ij) THEN
!                                DO j = 1, jnatoms
!                                        jsxp = atom_list(j,jm,js)%rp(1)
!                                        jsyp = atom_list(j,jm,js)%rp(2)
!                                        jszp = atom_list(j,jm,js)%rp(3)
!                                        k = (j-1) * n_i_exist
!                                        !$OMP SIMD
!                                        DO i = 1, n_i_exist
!                                                jrxp(i+k) = jsxp
!                                                jryp(i+k) = jsyp
!                                                jrzp(i+k) = jszp
!                                        END DO
!                                        !$OMP END SIMD
!                                END DO
!                        ELSE
!                                !$OMP SIMD PRIVATE(jsxp,jsyp,jszp)
!                                DO j = 1, jnatoms
!                                        jsxp = atom_list(j,jm,js)%rp(1)
!                                        jsyp = atom_list(j,jm,js)%rp(2)
!                                        jszp = atom_list(j,jm,js)%rp(3)
!                                        jrxp(j) = jsxp
!                                        jryp(j) = jsyp
!                                        jrzp(j) = jszp
!                                        DO i = 1, n_i_exist - 1
!                                                k = i*jnatoms
!                                                jrxp(j+k) = jsxp
!                                                jryp(j+k) = jsyp
!                                                jrzp(j+k) = jszp
!                                        END DO
!                                END DO
!                                !$OMP END SIMD
!                        END IF
!                ELSE
!                        !$OMP SIMD PRIVATE(jrp,jsxp,jsyp,jszp)
!                        DO j = 1, jnatoms
!                                jrp = atom_list(j,jm,js)%rp(1)
!                                jsxp = inv_h11*jrp
!                                jsyp = inv_h21*jrp
!                                jszp = inv_h31*jrp
!                                jrp = atom_list(j,jm,js)%rp(2)
!                                jsxp = jsxp + inv_h12*jrp
!                                jsyp = jsyp + inv_h22*jrp
!                                jszp = jszp + inv_h32*jrp
!                                jrp = atom_list(j,jm,js)%rp(3)
!                                jsxp = jsxp + inv_h13*jrp
!                                jsyp = jsyp + inv_h23*jrp
!                                jszp = jszp + inv_h33*jrp
!                                jrxp(j) = jsxp
!                                jryp(j) = jsyp
!                                jrzp(j) = jszp
!                                DO i = 1, n_i_exist-1
!                                        k = i*jnatoms
!                                        jrxp(i+k) = jsxp
!                                        jryp(i+k) = jsyp
!                                        jrzp(i+k) = jszp
!                                END DO
!                        END DO
!                        !$OMP END SIMD
!                END IF
!                !jrxp = RESHAPE(SPREAD(jrxp_vec,2,n_i_exist), (/ natompairs /))
!                !jryp = RESHAPE(SPREAD(jryp_vec,2,n_i_exist), (/ natompairs /))
!                !jrzp = RESHAPE(SPREAD(jrzp_vec,2,n_i_exist), (/ natompairs /))
!                IF (l_ortho) THEN
!                        DO ji = 1, natompairs
!                                dxp = ABS(jrxp(ji) - irxp(ji))
!                                dyp = ABS(jryp(ji) - iryp(ji))
!                                dzp = ABS(jrzp(ji) - irzp(ji))
!                                IF (dxp > hxl) dxp = dxp - xl
!                                IF (dyp > hyl) dyp = dyp - yl
!                                IF (dzp > hzl) dzp = dzp - zl
!                                dxp = dxp*dxp
!                                dxp = dxp + dyp*dyp
!                                dxp = dxp + dzp*dzp
!                                rijsq_vec(ji) = dxp
!                                IF (need_sqrt) rij_vec(ji) = SQRT(dxp)
!                                i_overlap = i_overlap .OR. dxp < rcut_lowsq
!                        END DO
!                ELSE
!                        DO ji = 1, natompairs
!                                dsp = jrxp(ji) - irxp(ji)
!                                IF (dsp > 0.5_DP) THEN
!                                        dsp = dsp - 1.0_DP
!                                ELSE IF (dsp < -0.5_DP) THEN
!                                        dsp = dsp + 1.0_DP
!                                END IF
!                                dxp = h11*dsp
!                                dyp = h21*dsp
!                                dzp = h31*dsp
!                                dsp = jryp(ji) - iryp(ji)
!                                IF (dsp > 0.5_DP) THEN
!                                        dsp = dsp - 1.0_DP
!                                ELSE IF (dsp < -0.5_DP) THEN
!                                        dsp = dsp + 1.0_DP
!                                END IF
!                                dxp = dxp + h12*dsp
!                                dyp = dyp + h22*dsp
!                                dzp = dzp + h32*dsp
!                                dsp = jrzp(ji) - irzp(ji)
!                                IF (dsp > 0.5_DP) THEN
!                                        dsp = dsp - 1.0_DP
!                                ELSE IF (dsp < -0.5_DP) THEN
!                                        dsp = dsp + 1.0_DP
!                                END IF
!                                dxp = dxp + h13*dsp
!                                dyp = dyp + h23*dsp
!                                dzp = dzp + h33*dsp
!                                dxp = dxp*dxp
!                                dxp = dxp + dyp*dyp
!                                dxp = dxp + dzp*dzp
!                                rijsq_vec(ji) = dxp
!                                IF (need_sqrt) rij_vec(ji) = SQRT(dxp)
!                                i_overlap = i_overlap .OR. dxp < rcut_lowsq
!                        END DO
!                END IF
!                IF (i_overlap) THEN
!                        overlap = .TRUE.
!                        CYCLE
!                END IF
!                IF (lj_charmm) THEN
!                        DO ji = 1, natompairs
!                                eps = vdw_p(ji,1) ! epsilon
!                                sigsq = vdw_p(ji,2) ! sigma**2
!                                rijsq = rijsq_vec(ji)
!                                l_interact = rijsq < this_vdw_rcutsq
!                                sigbyr2 = sigsq/rijsq ! sigma was already squared
!                                sigbyr6 = sigbyr2*sigbyr2*sigbyr2
!                                sigbyr12 = sigbyr6*sigbyr6
!                                sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
!                                IF (l_interact) vdw_energy = vdw_energy + eps * sigbyr12
!                                !CALL Compute_LJ_Charmm(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy)
!                        END DO
!                ELSE IF (lj_cut_shift) THEN
!                        DO ji = 1, natompairs
!                                rijsq = rijsq_vec(ji)
!                                eps = vdw_p(ji,1) ! 4*epsilon
!                                negsigsq = vdw_p(ji,2) ! -(sigma**2)
!                                shift_p1 = vdw_p(ji,3)
!                                l_interact = rijsq < this_vdw_rcutsq
!                                negsigbyr2 = negsigsq/rijsq
!                                rterm = negsigbyr2*negsigbyr2*negsigbyr2
!                                rterm = rterm + rterm*rterm
!                                IF (l_interact) vdw_energy = vdw_energy + eps * rterm - shift_p1
!                                !CALL Compute_LJ_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy)
!                        END DO
!                ELSE IF (lj_cut_switch) THEN
!                        DO ji = 1, natompairs
!                                rijsq = rijsq_vec(ji)
!                                eps = vdw_p(ji,1) ! 4*epsilon
!                                negsigsq = vdw_p(ji,2) ! -(sigma**2)
!                                l_interact = rijsq < this_vdw_rcutsq
!                                negsigbyr2 = negsigsq/rijsq
!                                rterm = negsigbyr2*negsigbyr2*negsigbyr2
!                                rterm = rterm + rterm*rterm
!                                roffsq_rijsq = roff_switch_sq(this_box) - rijsq
!                                IF (l_interact .AND. rijsq .GE. ron_switch_sq(this_box)) THEN
!                                        vdw_energy = vdw_energy + &
!                                                roffsq_rijsq*roffsq_rijsq * &
!                                                (switch_factor2(this_box)+2.0_DP*rijsq)*switch_factor1(this_box) * &
!                                                eps*rterm
!                                ELSE IF (l_interact) THEN
!                                        vdw_energy = vdw_energy + eps*rterm
!                                END IF
!                                !CALL Compute_LJ_Cut_Switch(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy,this_box)
!                        END DO
!                ELSE IF (lj_cut_shift_force) THEN
!                        DO ji = 1, natompairs
!                                rijsq = rijsq_vec(ji)
!                                rij = rij_vec(ji)
!                                eps = vdw_p(ji,1) ! 4*epsilon
!                                negsigsq = vdw_p(ji,2) ! -(sigma**2)
!                                shift_p1 = vdw_p(ji,3)
!                                shift_p2 = vdw_p(ji,4)
!                                l_interact = rijsq < this_vdw_rcutsq
!                                negsigbyr2 = negsigsq/rijsq
!                                rterm = negsigbyr2*negsigbyr2*negsigbyr2
!                                rterm = rterm + rterm*rterm
!                                IF (l_interact) vdw_energy = vdw_energy + eps * rterm - shift_p1 - &
!                                        (rcut_vdw(this_box)-rij)*shift_p2
!                                !CALL Compute_LJ_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy)
!                        END DO
!                ELSE IF (lj_cut) THEN
!                        DO ji = 1, natompairs
!                                rijsq = rijsq_vec(ji)
!                                eps = vdw_p(ji,1) ! 4*epsilon
!                                negsigsq = vdw_p(ji,2) ! -(sigma**2)
!                                l_interact = rijsq < this_vdw_rcutsq
!                                negsigbyr2 = negsigsq/rijsq
!                                rterm = negsigbyr2*negsigbyr2*negsigbyr2
!                                rterm = rterm + rterm*rterm
!                                IF (l_interact) vdw_energy = vdw_energy + eps * rterm
!                                !CALL Compute_LJ_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy)
!                        END DO
!                ELSE IF (mie_cut_shift) THEN
!                        DO ji = 1, natompairs
!                                rijsq = rijsq_vec(ji)
!                                epsig_n = vdw_p(ji,1) ! epsilon * mie_coeff * sigma ** n
!                                epsig_m = vdw_p(ji,2) ! epsilon * mie_coeff * sigma ** m
!                                mie_n = vdw_p(ji,3) ! already halved
!                                mie_m = vdw_p(ji,4) ! already halved
!                                shift_p1 = vdw_p(ji,5)
!                                l_interact = rijsq < this_vdw_rcutsq
!                                IF (l_interact) THEN
!                                        vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
!                                        vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
!                                        vdw_energy = vdw_energy - shift_p1
!                                END IF
!                                !CALL Compute_Mie_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2), &
!                                !        vdw_p(ji,3),vdw_p(ji,4),vdw_energy)
!                        END DO
!                ELSE IF (mie_cut) THEN
!                        DO ji = 1, natompairs
!                                rijsq = rijsq_vec(ji)
!                                epsig_n = vdw_p(ji,1) ! epsilon * mie_coeff * sigma ** n
!                                epsig_m = vdw_p(ji,2) ! epsilon * mie_coeff * sigma ** m
!                                mie_n = vdw_p(ji,3) ! already halved
!                                mie_m = vdw_p(ji,4) ! already halved
!                                l_interact = rijsq < this_vdw_rcutsq
!                                IF (l_interact) THEN
!                                        vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
!                                        vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
!                                END IF
!                                !CALL Compute_Mie_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2), &
!                                !        vdw_p(ji,3),vdw_p(ji,4),vdw_energy)
!                        END DO
!                END IF
!                IF (get_qq) THEN
!                        IF (l_charge_ewald) THEN
!                                DO ji = 1, natompairs
!                                        rij = rij_vec(ji)
!                                        cfqq = cfqq_vec(ji)
!                                        l_interact = rij < this_coul_rcut
!                                        IF (l_interact) THEN
!                                                qq_energy = qq_energy + ERFC(alpha_ewald(this_box)*rij) / rij * cfqq
!                                        END IF
!                                END DO
!                        ELSE IF (l_charge_dsf) THEN
!                                DO ji = 1, natompairs
!                                        rij = rij_vec(ji)
!                                        cfqq = cfqq_vec(ji)
!                                        l_interact = rij < this_coul_rcut
!                                        IF (l_interact) THEN
!                                                qq_energy = qq_energy + (dsf_factor2(this_box) * (rij - rcut_coul(this_box)) - &
!                                                        dsf_factor1(this_box) + &
!                                                        ERFC(alpha_dsf(this_box)*rij) / rij) * cfqq
!                                        END IF
!                                END DO
!                        ELSE IF (l_charge_cut) THEN
!                                DO ji = 1, natompairs
!                                        rij = rij_vec(ji)
!                                        cfqq = cfqq_vec(ji)
!                                        l_interact = rij < this_coul_rcut
!                                        IF (l_interact) THEN
!                                                qq_energy = qq_energy + cfqq / rij
!                                        END IF
!                                END DO
!                        END IF
!                END IF
!                IF (l_pair_store) THEN
!                        jsl = jsl_base + jm
!                        pair_nrg_vdw(isl,jsl) = vdw_energy
!                        pair_nrg_vdw(jsl,isl) = vdw_energy
!                        pair_nrg_qq(isl,jsl) = qq_energy
!                        pair_nrg_qq(jsl,isl) = qq_energy
!                END IF
!                i_vdw_energy = i_vdw_energy + vdw_energy
!                i_qq_energy = i_qq_energy + qq_energy
!                !IF (i_overlap) THEN
!                !        !$OMP CRITICAL
!                !        WRITE(*,*)
!                !        WRITE(*,*)
!                !        WRITE(*,*) i_overlap, l_order_ij
!                !        WRITE(*,*) im, jm, j_interact, i_mcstep
!                !        WRITE(*,*) vdw_energy, qq_energy
!                !        WRITE(*,*) rijsq_vec(1:natompairs)
!                !        WRITE(*,*)
!                !        WRITE(*,*) irxp(1:natompairs)
!                !        WRITE(*,*) jrxp(1:natompairs)
!                !        WRITE(*,*)
!                !        WRITE(*,*) iryp(1:natompairs)
!                !        WRITE(*,*) jryp(1:natompairs)
!                !        WRITE(*,*)
!                !        WRITE(*,*) irzp(1:natompairs)
!                !        WRITE(*,*) jrzp(1:natompairs)
!                !        !$OMP END CRITICAL
!                !END IF
!            END DO
!            !$OMP END DO
!
!            !$OMP DO SCHEDULE(STATIC)
!            DO j_interact = 1, n_interact
!                IF (overlap .OR. i_overlap) CYCLE
!                ij_overlap = .FALSE.
!                vdw_energy = 0.0_DP
!                qq_energy = 0.0_DP
!                jm = which_interact(j_interact)
!                IF (l_ortho) THEN
!                        IF (l_order_ij) THEN
!                                DO j = 1, jnatoms
!                                        jsxp = atom_list(j,jm,js)%rp(1)
!                                        jsyp = atom_list(j,jm,js)%rp(2)
!                                        jszp = atom_list(j,jm,js)%rp(3)
!                                        k = (j-1) * n_i_exist
!                                        !$OMP SIMD
!                                        DO i = 1, n_i_exist
!                                                jrxp(i+k) = jsxp
!                                                jryp(i+k) = jsyp
!                                                jrzp(i+k) = jszp
!                                        END DO
!                                        !$OMP END SIMD
!                                END DO
!                        ELSE
!                                !$OMP SIMD PRIVATE(jsxp,jsyp,jszp)
!                                DO j = 1, jnatoms
!                                        jsxp = atom_list(j,jm,js)%rp(1)
!                                        jsyp = atom_list(j,jm,js)%rp(2)
!                                        jszp = atom_list(j,jm,js)%rp(3)
!                                        jrxp(j) = jsxp
!                                        jryp(j) = jsyp
!                                        jrzp(j) = jszp
!                                        DO i = 1, n_i_exist - 1
!                                                k = i*jnatoms
!                                                jrxp(j+k) = jsxp
!                                                jryp(j+k) = jsyp
!                                                jrzp(j+k) = jszp
!                                        END DO
!                                END DO
!                                !$OMP END SIMD
!                        END IF
!                ELSE
!                        !$OMP SIMD PRIVATE(jrp,jsxp,jsyp,jszp)
!                        DO j = 1, jnatoms
!                                jrp = atom_list(j,jm,js)%rp(1)
!                                jsxp = inv_h11*jrp
!                                jsyp = inv_h21*jrp
!                                jszp = inv_h31*jrp
!                                jrp = atom_list(j,jm,js)%rp(2)
!                                jsxp = jsxp + inv_h12*jrp
!                                jsyp = jsyp + inv_h22*jrp
!                                jszp = jszp + inv_h32*jrp
!                                jrp = atom_list(j,jm,js)%rp(3)
!                                jsxp = jsxp + inv_h13*jrp
!                                jsyp = jsyp + inv_h23*jrp
!                                jszp = jszp + inv_h33*jrp
!                                jrxp(j) = jsxp
!                                jryp(j) = jsyp
!                                jrzp(j) = jszp
!                                DO i = 1, n_i_exist-1
!                                        k = i*jnatoms
!                                        jrxp(i+k) = jsxp
!                                        jryp(i+k) = jsyp
!                                        jrzp(i+k) = jszp
!                                END DO
!                        END DO
!                        !$OMP END SIMD
!                END IF
!                !jrxp = RESHAPE(SPREAD(jrxp_vec,2,n_i_exist), (/ natompairs /))
!                !jryp = RESHAPE(SPREAD(jryp_vec,2,n_i_exist), (/ natompairs /))
!                !jrzp = RESHAPE(SPREAD(jrzp_vec,2,n_i_exist), (/ natompairs /))
!                IF (l_ortho) THEN
!                        DO ji = 1, natompairs
!                                dxp = ABS(jrxp(ji) - irxp(ji))
!                                dyp = ABS(jryp(ji) - iryp(ji))
!                                dzp = ABS(jrzp(ji) - irzp(ji))
!                                IF (dxp > hxl) dxp = dxp - xl
!                                IF (dyp > hyl) dyp = dyp - yl
!                                IF (dzp > hzl) dzp = dzp - zl
!                                dxp = dxp*dxp
!                                dxp = dxp + dyp*dyp
!                                dxp = dxp + dzp*dzp
!                                rijsq_vec(ji) = dxp
!                                IF (need_sqrt) rij_vec(ji) = SQRT(dxp)
!                                i_overlap = i_overlap .OR. dxp < rcut_lowsq
!                        END DO
!                ELSE
!                        DO ji = 1, natompairs
!                                dsp = jrxp(ji) - irxp(ji)
!                                IF (dsp > 0.5_DP) THEN
!                                        dsp = dsp - 1.0_DP
!                                ELSE IF (dsp < -0.5_DP) THEN
!                                        dsp = dsp + 1.0_DP
!                                END IF
!                                dxp = h11*dsp
!                                dyp = h21*dsp
!                                dzp = h31*dsp
!                                dsp = jryp(ji) - iryp(ji)
!                                IF (dsp > 0.5_DP) THEN
!                                        dsp = dsp - 1.0_DP
!                                ELSE IF (dsp < -0.5_DP) THEN
!                                        dsp = dsp + 1.0_DP
!                                END IF
!                                dxp = dxp + h12*dsp
!                                dyp = dyp + h22*dsp
!                                dzp = dzp + h32*dsp
!                                dsp = jrzp(ji) - irzp(ji)
!                                IF (dsp > 0.5_DP) THEN
!                                        dsp = dsp - 1.0_DP
!                                ELSE IF (dsp < -0.5_DP) THEN
!                                        dsp = dsp + 1.0_DP
!                                END IF
!                                dxp = dxp + h13*dsp
!                                dyp = dyp + h23*dsp
!                                dzp = dzp + h33*dsp
!                                dxp = dxp*dxp
!                                dxp = dxp + dyp*dyp
!                                dxp = dxp + dzp*dzp
!                                rijsq_vec(ji) = dxp
!                                IF (need_sqrt) rij_vec(ji) = SQRT(dxp)
!                                i_overlap = i_overlap .OR. dxp < rcut_lowsq
!                        END DO
!                END IF
!                IF (i_overlap) THEN
!                        overlap = .TRUE.
!                        CYCLE
!                END IF
!                IF (lj_charmm) THEN
!                        DO ji = 1, natompairs
!                                eps = vdw_p(ji,1) ! epsilon
!                                sigsq = vdw_p(ji,2) ! sigma**2
!                                rijsq = rijsq_vec(ji)
!                                l_interact = rijsq < this_vdw_rcutsq
!                                sigbyr2 = sigsq/rijsq ! sigma was already squared
!                                sigbyr6 = sigbyr2*sigbyr2*sigbyr2
!                                sigbyr12 = sigbyr6*sigbyr6
!                                sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
!                                IF (l_interact) vdw_energy = vdw_energy + eps * sigbyr12
!                                !CALL Compute_LJ_Charmm(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy)
!                        END DO
!                ELSE IF (lj_cut_shift) THEN
!                        DO ji = 1, natompairs
!                                rijsq = rijsq_vec(ji)
!                                eps = vdw_p(ji,1) ! 4*epsilon
!                                negsigsq = vdw_p(ji,2) ! -(sigma**2)
!                                shift_p1 = vdw_p(ji,3)
!                                l_interact = rijsq < this_vdw_rcutsq
!                                negsigbyr2 = negsigsq/rijsq
!                                rterm = negsigbyr2*negsigbyr2*negsigbyr2
!                                rterm = rterm + rterm*rterm
!                                IF (l_interact) vdw_energy = vdw_energy + eps * rterm - shift_p1
!                                !CALL Compute_LJ_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy)
!                        END DO
!                ELSE IF (lj_cut_switch) THEN
!                        DO ji = 1, natompairs
!                                rijsq = rijsq_vec(ji)
!                                eps = vdw_p(ji,1) ! 4*epsilon
!                                negsigsq = vdw_p(ji,2) ! -(sigma**2)
!                                l_interact = rijsq < this_vdw_rcutsq
!                                negsigbyr2 = negsigsq/rijsq
!                                rterm = negsigbyr2*negsigbyr2*negsigbyr2
!                                rterm = rterm + rterm*rterm
!                                roffsq_rijsq = roff_switch_sq(this_box) - rijsq
!                                IF (l_interact .AND. rijsq .GE. ron_switch_sq(this_box)) THEN
!                                        vdw_energy = vdw_energy + &
!                                                roffsq_rijsq*roffsq_rijsq * &
!                                                (switch_factor2(this_box)+2.0_DP*rijsq)*switch_factor1(this_box) * &
!                                                eps*rterm
!                                ELSE IF (l_interact) THEN
!                                        vdw_energy = vdw_energy + eps*rterm
!                                END IF
!                                !CALL Compute_LJ_Cut_Switch(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy,this_box)
!                        END DO
!                ELSE IF (lj_cut_shift_force) THEN
!                        DO ji = 1, natompairs
!                                rijsq = rijsq_vec(ji)
!                                rij = rij_vec(ji)
!                                eps = vdw_p(ji,1) ! 4*epsilon
!                                negsigsq = vdw_p(ji,2) ! -(sigma**2)
!                                shift_p1 = vdw_p(ji,3)
!                                shift_p2 = vdw_p(ji,4)
!                                l_interact = rijsq < this_vdw_rcutsq
!                                negsigbyr2 = negsigsq/rijsq
!                                rterm = negsigbyr2*negsigbyr2*negsigbyr2
!                                rterm = rterm + rterm*rterm
!                                IF (l_interact) vdw_energy = vdw_energy + eps * rterm - shift_p1 - &
!                                        (rcut_vdw(this_box)-rij)*shift_p2
!                                !CALL Compute_LJ_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy)
!                        END DO
!                ELSE IF (lj_cut) THEN
!                        DO ji = 1, natompairs
!                                rijsq = rijsq_vec(ji)
!                                eps = vdw_p(ji,1) ! 4*epsilon
!                                negsigsq = vdw_p(ji,2) ! -(sigma**2)
!                                l_interact = rijsq < this_vdw_rcutsq
!                                negsigbyr2 = negsigsq/rijsq
!                                rterm = negsigbyr2*negsigbyr2*negsigbyr2
!                                rterm = rterm + rterm*rterm
!                                IF (l_interact) vdw_energy = vdw_energy + eps * rterm
!                                !CALL Compute_LJ_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2),vdw_energy)
!                        END DO
!                ELSE IF (mie_cut_shift) THEN
!                        DO ji = 1, natompairs
!                                rijsq = rijsq_vec(ji)
!                                epsig_n = vdw_p(ji,1) ! epsilon * mie_coeff * sigma ** n
!                                epsig_m = vdw_p(ji,2) ! epsilon * mie_coeff * sigma ** m
!                                mie_n = vdw_p(ji,3) ! already halved
!                                mie_m = vdw_p(ji,4) ! already halved
!                                shift_p1 = vdw_p(ji,5)
!                                l_interact = rijsq < this_vdw_rcutsq
!                                IF (l_interact) THEN
!                                        vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
!                                        vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
!                                        vdw_energy = vdw_energy - shift_p1
!                                END IF
!                                !CALL Compute_Mie_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2), &
!                                !        vdw_p(ji,3),vdw_p(ji,4),vdw_energy)
!                        END DO
!                ELSE IF (mie_cut) THEN
!                        DO ji = 1, natompairs
!                                rijsq = rijsq_vec(ji)
!                                epsig_n = vdw_p(ji,1) ! epsilon * mie_coeff * sigma ** n
!                                epsig_m = vdw_p(ji,2) ! epsilon * mie_coeff * sigma ** m
!                                mie_n = vdw_p(ji,3) ! already halved
!                                mie_m = vdw_p(ji,4) ! already halved
!                                l_interact = rijsq < this_vdw_rcutsq
!                                IF (l_interact) THEN
!                                        vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
!                                        vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
!                                END IF
!                                !CALL Compute_Mie_Cut(rijsq,vdw_p(ji,1),vdw_p(ji,2), &
!                                !        vdw_p(ji,3),vdw_p(ji,4),vdw_energy)
!                        END DO
!                END IF
!                IF (get_qq) THEN
!                        IF (l_charge_ewald) THEN
!                                DO ji = 1, natompairs
!                                        rij = rij_vec(ji)
!                                        cfqq = cfqq_vec(ji)
!                                        l_interact = rij < this_coul_rcut
!                                        IF (l_interact) THEN
!                                                qq_energy = qq_energy + ERFC(alpha_ewald(this_box)*rij) / rij * cfqq
!                                        END IF
!                                END DO
!                        ELSE IF (l_charge_dsf) THEN
!                                DO ji = 1, natompairs
!                                        rij = rij_vec(ji)
!                                        cfqq = cfqq_vec(ji)
!                                        l_interact = rij < this_coul_rcut
!                                        IF (l_interact) THEN
!                                                qq_energy = qq_energy + (dsf_factor2(this_box) * (rij - rcut_coul(this_box)) - &
!                                                        dsf_factor1(this_box) + &
!                                                        ERFC(alpha_dsf(this_box)*rij) / rij) * cfqq
!                                        END IF
!                                END DO
!                        ELSE IF (l_charge_cut) THEN
!                                DO ji = 1, natompairs
!                                        rij = rij_vec(ji)
!                                        cfqq = cfqq_vec(ji)
!                                        l_interact = rij < this_coul_rcut
!                                        IF (l_interact) THEN
!                                                qq_energy = qq_energy + cfqq / rij
!                                        END IF
!                                END DO
!                        END IF
!                END IF
!                IF (l_pair_store) THEN
!                        jsl = jsl_base + jm
!                        pair_nrg_vdw(isl,jsl) = vdw_energy
!                        pair_nrg_vdw(jsl,isl) = vdw_energy
!                        pair_nrg_qq(isl,jsl) = qq_energy
!                        pair_nrg_qq(jsl,isl) = qq_energy
!                END IF
!                i_vdw_energy = i_vdw_energy + vdw_energy
!                i_qq_energy = i_qq_energy + qq_energy
!                !IF (i_overlap) THEN
!                !        !$OMP CRITICAL
!                !        WRITE(*,*)
!                !        WRITE(*,*)
!                !        WRITE(*,*) i_overlap, l_order_ij
!                !        WRITE(*,*) im, jm, j_interact, i_mcstep
!                !        WRITE(*,*) vdw_energy, qq_energy
!                !        WRITE(*,*) rijsq_vec(1:natompairs)
!                !        WRITE(*,*)
!                !        WRITE(*,*) irxp(1:natompairs)
!                !        WRITE(*,*) jrxp(1:natompairs)
!                !        WRITE(*,*)
!                !        WRITE(*,*) iryp(1:natompairs)
!                !        WRITE(*,*) jryp(1:natompairs)
!                !        WRITE(*,*)
!                !        WRITE(*,*) irzp(1:natompairs)
!                !        WRITE(*,*) jrzp(1:natompairs)
!                !        !$OMP END CRITICAL
!                !END IF
!            END DO
!            !$OMP END DO
    END DO
    !IF (i_mcstep < 10 .AND. .NOT. cbmc_flag) WRITE(*,*) overlap,i_overlap, i_vdw_energy, i_qq_energy, n_interact
    !IF (i_mcstep < 10 .AND. cbmc_flag) WRITE(*,*) overlap,i_overlap, i_vdw_energy, i_qq_energy, n_interact
    !$OMP END PARALLEL
    overlap = i_overlap
    !IF (i_mcstep < 10 .AND. cbmc_flag) WRITE(*,*) overlap,i_overlap, i_vdw_energy, i_qq_energy, n_interact
    !IF (i_mcstep < 10 .AND. .NOT. cbmc_flag) WRITE(*,*) overlap,i_overlap, i_vdw_energy, i_qq_energy, n_interact
  END SUBROUTINE Compute_Molecule_Nonbond_Inter_Energy_Vectorized

  FUNCTION Pack_Parallel_INT32(mask,asize,psize,a) RESULT(p)
          INTEGER, INTENT(IN) :: asize
          INTEGER, DIMENSION(asize), OPTIONAL, INTENT(IN) :: a
          INTEGER, DIMENSION(asize) :: p
          LOGICAL(1), DIMENSION(asize), INTENT(IN) :: mask
          INTEGER, INTENT(OUT) :: psize
          LOGICAL(1) :: l_a_present
          INTEGER :: i, j
          INTEGER :: ithread, bsize, nthreads_used, tcount, ibase, ctcount
          INTEGER, DIMENSION(MAX(16,asize/global_nthreads+1)) :: b
          INTEGER, DIMENSION(global_nthreads) :: tcount_vec, ibase_vec
          bsize = SIZE(b)
          nthreads_used = asize/bsize
          IF (MOD(asize,bsize) .NE. 0) nthreads_used = nthreads_used + 1
          ithread = 1
          l_a_present = PRESENT(a)
          !$OMP PARALLEL NUM_THREADS(nthreads_used) &
          !$OMP PRIVATE(b, tcount, ithread, ibase, ctcount, j)
          !$ ithread = OMP_GET_THREAD_NUM() + 1
          tcount = 0
          IF (l_a_present) THEN
                  !$OMP DO SCHEDULE(STATIC,bsize)
                  DO i = 1, asize
                        IF (mask(i)) THEN
                                tcount = tcount + 1
                                b(tcount) = a(i)
                        END IF
                  END DO
                  !$OMP END DO
          ELSE
                  !$OMP DO SCHEDULE(STATIC,bsize)
                  DO i = 1, asize
                        IF (mask(i)) THEN
                                tcount = tcount + 1
                                b(tcount) = i
                        END IF
                  END DO
                  !$OMP END DO
          END IF
          tcount_vec(ithread) = tcount
          !$OMP BARRIER
          !$OMP SINGLE
          ctcount = 0
          DO j = 1, nthreads_used
                ibase_vec(j) = ctcount
                ctcount = ctcount + tcount_vec(j)
          END DO
          psize = ctcount
          !$OMP END SINGLE
          ibase = ibase_vec(ithread)
          p(ibase+1:ibase+tcount) = b(1:tcount)
          !$OMP END PARALLEL
  END FUNCTION Pack_Parallel_INT32

  SUBROUTINE Compute_LJ_Charmm(rijsq,eps,sigsq,vdw_energy)
          !!!!!$OMP DECLARE SIMD(Compute_LJ_Charmm)
          REAL(DP) :: eps ! epsilon
          REAL(DP) :: sigsq ! sigma**2
          REAL(DP) :: rijsq, vdw_energy
          REAL(DP) :: sigbyr2, sigbyr6, sigbyr12
          sigbyr2 = sigsq/rijsq ! sigma was already squared
          sigbyr6 = sigbyr2*sigbyr2*sigbyr2
          sigbyr12 = sigbyr6*sigbyr6
          sigbyr12 = sigbyr12 - 2.0_DP*sigbyr6
          vdw_energy = vdw_energy + eps * sigbyr12
  END SUBROUTINE Compute_LJ_Charmm

  SUBROUTINE Compute_LJ_Cut_Shift(rijsq,eps,negsigsq,vdw_energy,this_box)
          !!!!$OMP DECLARE SIMD(Compute_LJ_Cut_Shift) UNIFORM(this_box)
          REAL(DP) :: eps ! 4*epsilon
          REAL(DP) :: negsigsq ! -(sigma**2)
          INTEGER :: this_box
          REAL(DP) :: rijsq, vdw_energy
          REAL(DP) :: negsigbyr2, rterm
          negsigbyr2 = negsigsq/rijsq
          rterm = negsigbyr2*negsigbyr2*negsigbyr2
          rterm = rterm + rterm*rterm
          vdw_energy = vdw_energy + eps * rterm
          negsigbyr2 = negsigsq/rcut_vdwsq(this_box)
          rterm = negsigbyr2*negsigbyr2*negsigbyr2
          rterm = rterm + rterm*rterm
          vdw_energy = vdw_energy - eps*rterm
  END SUBROUTINE Compute_LJ_Cut_Shift

  SUBROUTINE Compute_LJ_Cut_Switch(rijsq,eps,negsigsq,vdw_energy,this_box)
          !!!!!$OMP DECLARE SIMD(Compute_LJ_Cut_Switch) UNIFORM(this_box)
          REAL(DP) :: eps ! 4*epsilon
          REAL(DP) :: negsigsq ! -(sigma**2)
          INTEGER :: this_box
          REAL(DP) :: rijsq, vdw_energy
          REAL(DP) :: negsigbyr2, rterm, roffsq_rijsq
          negsigbyr2 = negsigsq/rijsq
          rterm = negsigbyr2*negsigbyr2*negsigbyr2
          rterm = rterm + rterm*rterm
          IF (rijsq .GE. ron_switch_sq(this_box)) THEN
                  roffsq_rijsq = roff_switch_sq(this_box) - rijsq
                  vdw_energy = vdw_energy + &
                          roffsq_rijsq*roffsq_rijsq * &
                          (switch_factor2(this_box)+2.0_DP*rijsq)*switch_factor1(this_box) * &
                          eps*rterm
          ELSE
                  vdw_energy = vdw_energy + eps*rterm
          END IF
  END SUBROUTINE Compute_LJ_Cut_Switch

  SUBROUTINE Compute_LJ_Cut_Shift_Force(rijsq,rij,eps,negsigsq,vdw_energy,this_box)
          !!!!$OMP DECLARE SIMD(Compute_LJ_Cut_Shift_Force) UNIFORM(this_box)
          REAL(DP) :: eps ! 4*epsilon
          REAL(DP) :: negsigsq ! -(sigma**2)
          INTEGER :: this_box
          REAL(DP) :: rijsq, rij, vdw_energy
          REAL(DP) :: negsigbyr2, rterm, rterm2
          negsigbyr2 = negsigsq/rijsq
          rterm = negsigbyr2*negsigbyr2*negsigbyr2
          rterm = rterm + rterm*rterm
          vdw_energy = vdw_energy + eps * rterm
          negsigbyr2 = negsigsq/rcut_vdwsq(this_box)
          rterm = negsigbyr2*negsigbyr2*negsigbyr2
          rterm2 = rterm * rterm
          rterm2 = rterm + 2.0_DP * rterm2
          rterm = rterm + rterm*rterm
          vdw_energy = vdw_energy - eps * rterm
          vdw_energy = vdw_energy - &
                  (rij - rcut_vdw(this_box)) * &
                  -6.0_DP * eps * rterm2 / rcut_vdw(this_box)
  END SUBROUTINE Compute_LJ_Cut_Shift_Force

  SUBROUTINE Compute_LJ_Cut(rijsq,eps,negsigsq,vdw_energy)
          !!!!$OMP DECLARE SIMD(Compute_LJ_Cut)
          REAL(DP) :: eps ! 4*epsilon
          REAL(DP) :: negsigsq ! -(sigma**2)
          REAL(DP) :: rijsq, vdw_energy
          REAL(DP) :: negsigbyr2, rterm
          negsigbyr2 = negsigsq/rijsq
          rterm = negsigbyr2*negsigbyr2*negsigbyr2
          rterm = rterm + rterm*rterm
          vdw_energy = vdw_energy + eps * rterm
  END SUBROUTINE Compute_LJ_Cut

  SUBROUTINE Compute_Mie_Cut(rijsq,epsig_n,epsig_m,mie_n,mie_m,vdw_energy)
          !!!!$OMP DECLARE SIMD(Compute_Mie_Cut)
          REAL(DP) :: epsig_n, epsig_m ! epsilon * mie_coeff * sigma**exponent
          REAL(DP) :: mie_n, mie_m ! already halved and negative
          REAL(DP) :: rijsq, vdw_energy
          ! exponents were already halved so there's no need for SQRT
          vdw_energy = vdw_energy + epsig_n * rijsq**mie_n
          vdw_energy = vdw_energy - epsig_m * rijsq**mie_m
  END SUBROUTINE Compute_Mie_Cut


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
                ixp = these_atoms_i(ia)%rp(1)
                iyp = these_atoms_i(ia)%rp(2)
                izp = these_atoms_i(ia)%rp(3)
                DO ja = 1, natoms(js)
                        IF (.NOT. these_atoms_j(ja)%exist) CYCLE
                        CALL Minimum_Image_Separation(this_box, &
                                ixp - these_atoms_j(ja)%rp(1), &
                                iyp - these_atoms_j(ja)%rp(2), &
                                izp - these_atoms_j(ja)%rp(3), &
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

  REAL(SP) FUNCTION Compute_Cell_List_CBMC_nrg(irp,cp,ia,is,this_box)
          REAL(SP) :: irp(3)
          INTEGER, INTENT(IN) :: cp(3), ia, is, this_box
          ! The following parameters are here so you can optimize speed by setting them.
          !     The best options probably depend on your processor architecture, unfortunately.
          !     They are currently not used at all
          LOGICAL, PARAMETER :: pregather_lj = .FALSE., lj_p_readfirstdim = .FALSE.
          INTEGER :: xi, yi, zi, vlen, i
          REAL(SP) :: nrg

          INTEGER :: dmult, isolute, rsqsol
          REAL(SP) :: this_shifter, rsq_shift, drp
          REAL(SP) :: dxp, dyp, dzp

          REAL(SP) :: rsq, rcutsq
          REAL(SP), DIMENSION(cbmc_max_interact) :: rsq_vec
          !REAL(SP), DIMENSION(cbmc_max_interact,2) :: ppljp
          INTEGER :: itype, jtype, this_int_vdw_style, this_int_vdw_sum_style
          REAL(SP) :: eps, sigsq, sigbyr2, sigbyr6, sigbyr12, nrg_vdw, rterm, rterm2, negsigsq, negsigbyr2, roffsq_rijsq
          REAL(SP) :: mie_m, mie_n, lnrsq, epsig_n, epsig_m
          REAL(SP) :: icharge, jcharge, invr, rij, nrg_qq, dsf_const, alpha_ewald_sp
          REAL(SP) :: const1, const2, const3, const4
          INTEGER :: this_int_charge_sum_style
          LOGICAL :: i_get_coul

          xi = cp(1)
          yi = cp(2)
          zi = cp(3)
          vlen = cbmc_cell_n_interact(xi,yi,zi,this_box)
          !DIR$ ASSUME (MOD(vlen,8) .EQ. 0)
          nrg = 0.0

          IF (precalc_atompair_nrg) THEN
                  dmult = SIZE(atompair_nrg_table,1)
                  this_shifter = -sp_rcut_lowsq !-(rsq_shifter + rsq_step)
                  isolute = species_list(is)%solute_base + ia
                  !DIR$ VECTOR ALIGNED
                  DO i = 1, vlen
                        rsq_shift = this_shifter
                        drp = cbmc_cell_rsp(i,1,xi,yi,zi,this_box) - irp(1)
                        rsq_shift = rsq_shift + drp*drp
                        drp = cbmc_cell_rsp(i,2,xi,yi,zi,this_box) - irp(2)
                        rsq_shift = rsq_shift + drp*drp
                        drp = cbmc_cell_rsp(i,3,xi,yi,zi,this_box) - irp(3)
                        rsq_shift = rsq_shift + drp*drp
                        rsqsol = INT(MIN(rsq_shift*inv_rsq_step_sp,atompair_nrg_res_sp)) + &
                                cbmc_cell_ti(i,xi,yi,zi,this_box)*dmult
                        nrg = nrg + atompair_nrg_table_reduced(rsqsol,isolute,this_box)
                  END DO
                  Compute_Cell_List_CBMC_nrg = nrg
                  RETURN
          END IF
          rcutsq = REAL(rcut_cbmcsq(this_box),SP)
          !DIR$ ASSUME_ALIGNED rsq_vec:32
          !DIR$ VECTOR ALIGNED
          DO i = 1, vlen
                drp = cbmc_cell_rsp(i,1,xi,yi,zi,this_box) - irp(1)
                rsq = drp*drp 
                drp = cbmc_cell_rsp(i,2,xi,yi,zi,this_box) - irp(2)
                rsq = rsq + drp*drp
                drp = cbmc_cell_rsp(i,3,xi,yi,zi,this_box) - irp(3)
                rsq = rsq + drp*drp
                rsq_vec(i) = rsq
          END DO
          itype = nonbond_list(ia,is)%atom_type_number
          this_int_vdw_style = MIN(itype,1) * int_vdw_style(this_box)
          this_int_vdw_sum_style = MIN(itype,1) * int_vdw_sum_style(this_box)
          icharge = nonbond_list(ia,is)%charge
          i_get_coul = icharge .NE. 0.0_SP .AND. int_charge_style(this_box) .NE. charge_none .AND. species_list(is)%l_coul_cbmc
          this_int_charge_sum_style = MERGE(int_charge_sum_style(this_box),0,i_get_coul)
          nrg_vdw = 0.0_SP
          nrg_qq = 0.0_SP
          !ppljp(1:vlen,1:2) = TRANSPOSE(ppvdwp_table2_sp(1:2,cbmc_cell_atomtypes(1:vlen,xi,yi,zi,this_box),itype,this_box))
          !!DIR$ VECTOR MULTIPLE_GATHER_SCATTER_BY_SHUFFLES, ALIGNED
          !!$OMP SIMD
          !DO i = 1, vlen
          !      jtype = cbmc_cell_atomtypes(i,xi,yi,zi,this_box)
          !      eps = ppvdwp_table2_sp(1,jtype,itype,this_box)
          !      negsigsq = ppvdwp_table2_sp(2,jtype,itype,this_box)
          !      ppljp(i,1) = eps
          !      ppljp(i,2) = negsigsq
          !END DO
          !!$OMP END SIMD
          SELECT CASE(this_int_vdw_style)
          CASE(vdw_lj)
                  SELECT CASE (this_int_vdw_sum_style)
                  CASE(vdw_charmm)
                          !DIR$ VECTOR ALIGNED
                          !$OMP SIMD PRIVATE(rsq,jtype,eps,sigsq,sigbyr2,sigbyr6) &
                          !$OMP REDUCTION(+:nrg_vdw)
                          DO i = 1, vlen
                                rsq = rsq_vec(i)
                                jtype = cbmc_cell_atomtypes(i,xi,yi,zi,this_box)
                                eps = ppvdwp_table2_sp(1,jtype,itype,this_box) ! epsilon
                                sigsq = ppvdwp_table2_sp(2,jtype,itype,this_box) ! sigma**2
                                !eps = ppljp(i,1)
                                !sigsq = ppljp(i,2)
                                sigbyr2 = sigsq/rsq ! sigma was already squared
                                sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                nrg_vdw = nrg_vdw + eps * (sigbyr6*sigbyr6 - sigbyr6 - sigbyr6) ! eps was not multiplied by 4
                                !IF (rsq < rcutsq) THEN
                                !  jtype = cbmc_cell_atomtypes(i,xi,yi,zi,this_box)
                                !  eps = ppvdwp_table2_sp(1,jtype,itype,this_box) ! epsilon
                                !  sigsq = ppvdwp_table2_sp(2,jtype,itype,this_box) ! sigma**2
                                !  sigbyr2 = sigsq/rsq ! sigma was already squared
                                !  sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                !  nrg_vdw = nrg_vdw + eps * (sigbyr6*sigbyr6 - sigbyr6 - sigbyr6) ! eps was not multiplied by 4
                                !END IF
                          END DO
                          !$OMP END SIMD
                  CASE(vdw_cut_shift)
                          !DIR$ VECTOR MULTIPLE_GATHER_SCATTER_BY_SHUFFLES
                          DO i = 1, vlen
                                rsq = rsq_vec(i)
                                jtype = cbmc_cell_atomtypes(i,xi,yi,zi,this_box)
                                IF (rsq < rcutsq) THEN
                                  eps = ppvdwp_table2_sp(1,jtype,itype,this_box) ! 4*epsilon
                                  negsigsq = ppvdwp_table2_sp(2,jtype,itype,this_box) ! -(sigma**2)
                                  negsigbyr2 = negsigsq/rsq
                                  rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                  rterm = rterm + rterm*rterm
                                  nrg_vdw = nrg_vdw + eps*rterm
                                  negsigbyr2 = negsigsq*inv_rcut_vdwsq_sp(this_box)
                                  rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                  rterm = rterm + rterm*rterm
                                  nrg_vdw = nrg_vdw - eps*rterm
                                END IF
                          END DO
                  CASE(vdw_cut_switch)
                          const1 = switch_factor1(this_box)*switch_factor2(this_box)
                          const2 = switch_factor1(this_box)*2
                          const3 = ron_switch_sq(this_box)
                          const4 = roff_switch_sq(this_box)
                          !DIR$ VECTOR MULTIPLE_GATHER_SCATTER_BY_SHUFFLES
                          DO i = 1, vlen
                                rsq = rsq_vec(i)
                                jtype = cbmc_cell_atomtypes(i,xi,yi,zi,this_box)
                                IF (rsq < rcutsq) THEN
                                  eps = ppvdwp_table2_sp(1,jtype,itype,this_box) ! 4*epsilon
                                  negsigsq = ppvdwp_table2_sp(2,jtype,itype,this_box) ! -(sigma**2)
                                  negsigbyr2 = negsigsq/rsq
                                  rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                  rterm = rterm + rterm*rterm
                                  IF (rsq .GE. const3) THEN
                                          roffsq_rijsq = const4 - rsq
                                          eps = &
                                                  roffsq_rijsq*roffsq_rijsq * &
                                                  (const1+const2*rsq) * eps
                                  END IF
                                  nrg_vdw = nrg_vdw + eps*rterm
                                END IF
                          END DO
                  CASE(vdw_cut_shift_force)
                          const1 = 36.0_SP*inv_rcut_vdwsq_sp(this_box)
                          !DIR$ VECTOR MULTIPLE_GATHER_SCATTER_BY_SHUFFLES
                          DO i = 1, vlen
                                rsq = rsq_vec(i)
                                jtype = cbmc_cell_atomtypes(i,xi,yi,zi,this_box)
                                IF (rsq < rcutsq) THEN
                                  eps = ppvdwp_table_sp(jtype,itype,1,this_box) ! 4*epsilon
                                  negsigsq = ppvdwp_table_sp(jtype,itype,2,this_box) ! -(sigma**2)
                                  !eps = ppvdwp_table2_sp(1,jtype,itype,this_box) ! 4*epsilon
                                  !negsigsq = ppvdwp_table2_sp(2,jtype,itype,this_box) ! -(sigma**2)
                                  negsigbyr2 = negsigsq/rsq
                                  rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                  rterm = rterm + rterm*rterm
                                  nrg_vdw = nrg_vdw + eps * rterm
                                  negsigbyr2 = negsigsq*inv_rcut_vdwsq_sp(this_box)
                                  rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                  rterm2 = rterm * rterm
                                  rterm2 = rterm + 2.0_SP * rterm2
                                  rterm = rterm + rterm*rterm
                                  nrg_vdw = nrg_vdw - eps * rterm
                                  nrg_vdw = nrg_vdw + &
                                          (SQRT(rsq*const1) - 6.0_SP) * &
                                          eps * rterm2
                                  !nrg_vdw = nrg_vdw - &
                                  !        (SQRT(rsq) - rcut_vdw_sp(this_box)) * &
                                  !        -6.0_SP * eps * rterm2 * inv_rcut_vdw_sp(this_box)
                                END IF
                          END DO
                  CASE(vdw_cut,vdw_cut_tail)
                          !DIR$ VECTOR ALIGNED
                          !$OMP SIMD PRIVATE(eps,negsigsq,negsigbyr2,rterm,rsq,jtype) &
                          !$OMP REDUCTION(+:nrg_vdw)
                          DO i = 1, vlen
                                rsq = rsq_vec(i)
                                jtype = cbmc_cell_atomtypes(i,xi,yi,zi,this_box)
                                !eps = ppvdwp_table_sp(jtype,itype,1,this_box) ! 4*epsilon
                                !negsigsq = ppvdwp_table_sp(jtype,itype,2,this_box) ! -(sigma**2)
                                IF (rsq < rcutsq) THEN
                                  !eps = ppvdwp_table2_sp(1,jtype,itype,this_box) ! 4*epsilon
                                  !negsigsq = ppvdwp_table2_sp(2,jtype,itype,this_box) ! -(sigma**2)
                                  eps = ppvdwp_table_sp(jtype,itype,1,this_box) ! 4*epsilon
                                  negsigsq = ppvdwp_table_sp(jtype,itype,2,this_box) ! -(sigma**2)
                                  negsigbyr2 = negsigsq/rsq
                                  rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                  rterm = rterm + rterm*rterm
                                  ! nrg = eps * rterm
                                  nrg_vdw = nrg_vdw + eps * rterm
                                END IF
                          END DO
                          !$OMP END SIMD
                  END SELECT
          CASE(vdw_mie)
                  IF (.NOT. l_nonuniform_exponents) THEN
                          mie_n = ppvdwp_table2_sp(3,itype,itype,this_box)
                          mie_m = ppvdwp_table2_sp(4,itype,itype,this_box)
                  END IF
                  DO i = 1, vlen
                        rsq = rsq_vec(i)
                        jtype = cbmc_cell_atomtypes(i,xi,yi,zi,this_box)
                        IF (rsq < rcutsq) THEN
                          lnrsq = LOG(rsq)
                          epsig_n = ppvdwp_table2_sp(1,jtype,itype,this_box)  ! epsilon * mie_coeff * sigma ** n
                          epsig_m = ppvdwp_table2_sp(2,jtype,itype,this_box) ! epsilon * mie_coeff * sigma ** m
                          IF (l_nonuniform_exponents) THEN
                                  mie_n = ppvdwp_table2_sp(3,jtype,itype,this_box) ! already halved
                                  mie_m = ppvdwp_table2_sp(4,jtype,itype,this_box) ! already halved
                          END IF
                          nrg_vdw = nrg_vdw + epsig_n*EXP(lnrsq*mie_n) - epsig_m*EXP(lnrsq*mie_m)
                          IF (int_vdw_sum_style(this_box) == vdw_cut_shift) THEN
                                  nrg_vdw = nrg_vdw - ppvdwp_table_sp(jtype,itype,5,this_box)
                                  !nrg = nrg + epsig_n * rcut_vdwsq(this_box)**mie_n
                                  !nrg = nrg - epsig_m * rcut_vdwsq(this_box)**mie_m
                          END IF
                        END IF
                  END DO
          END SELECT
          SELECT CASE (this_int_charge_sum_style)
          CASE(charge_ewald)
                  alpha_ewald_sp = REAL(alpha_ewald(this_box),SP)
                  !DIR$ VECTOR ALIGNED
                  DO i = 1, vlen
                        rsq = rsq_vec(i)
                        jcharge = cbmc_cell_rsp(i,4,xi,yi,zi,this_box)
                        IF (rsq < rcutsq) THEN
                                invr = 1.0 / SQRT(rsq)
                                rij = rsq * invr
                                nrg_qq = nrg_qq + ERFC(alpha_ewald_sp*rij)*invr*jcharge
                        END IF
                  END DO
          CASE(charge_dsf)
                  dsf_const = -dsf_factor2(this_box)*(rcut_coul(this_box)+dsf_factor1(this_box))
                  DO i = 1, vlen
                        rsq = rsq_vec(i)
                        IF (rsq < rcutsq) THEN
                                jcharge = cbmc_cell_rsp(i,4,xi,yi,zi,this_box)
                                invr = 1.0 / SQRT(rsq)
                                rij = rsq*invr
                                nrg = ERFC(alpha_dsf(this_box)*rij)*invr + dsf_const + &
                                        rij*dsf_factor2(this_box)
                                nrg_qq = nrg_qq + jcharge * nrg
                        END IF
                  END DO
          CASE(charge_cut)
                  DO i = 1, vlen
                        rsq = rsq_vec(i)
                        IF (rsq < rcutsq) THEN
                                jcharge = cbmc_cell_rsp(i,4,xi,yi,zi,this_box)
                                invr = 1.0 / SQRT(rsq)
                                nrg_qq = nrg_qq + jcharge*invr
                        END IF
                  END DO
          END SELECT
          Compute_Cell_List_CBMC_nrg = nrg_vdw + icharge*charge_factor*nrg_qq
  END FUNCTION Compute_Cell_List_CBMC_nrg

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
        rxijp = these_atoms_i(ia)%rp(1) - these_atoms_j(ja)%rp(1)
        ryijp = these_atoms_i(ia)%rp(2) - these_atoms_j(ja)%rp(2)
        rzijp = these_atoms_i(ia)%rp(3) - these_atoms_j(ja)%rp(3)

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
         ELSEIF (PRESENT(Eij_qq_o)) THEN
             Eij_qq_o = 0.0_DP
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
          hdotr = hx(i,ibox) * atom_list(ia,im,is)%rp(1) + &
                  hy(i,ibox) * atom_list(ia,im,is)%rp(2) + &
                  hz(i,ibox) * atom_list(ia,im,is)%rp(3)

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
          hdotr = hx(i,ibox) * atom_list(ia,im,is)%rp(1) + &
                  hy(i,ibox) * atom_list(ia,im,is)%rp(2) + &
                  hz(i,ibox) * atom_list(ia,im,is)%rp(3)

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
    REAL(DP), DIMENSION(natoms(is)) :: q, rxp, ryp, rzp
    REAL(DP), DIMENSION(nvecs(ibox)) :: this_cos_sum, this_sin_sum
    REAL(DP) :: hdotr, trigfactor, ihx, ihy, ihz, qia
    REAL(DP) :: cos_sum_i, sin_sum_i, E_recip_redux
    INTEGER :: i, ia
    !DIR$ ASSUME_ALIGNED q:32, rxp:32, ryp:32, rzp:32, this_cos_sum:32, this_sin_sum:32

    q = nonbond_list(1:natoms(is),is)%charge

    ! Initialize variables
    !E_reciprocal = 0.0_DP
    E_recip_redux = 0.0_DP

    rxp = widom_atoms%rp(1)
    ryp = widom_atoms%rp(2)
    rzp = widom_atoms%rp(3)

!    DO i = 1, nvecs(ibox)
!
!      hdotr = hx(i,ibox) * rxp + &
!              hy(i,ibox) * ryp + &
!              hz(i,ibox) * rzp
!
!      cos_sum_i = cos_sum(i,ibox) + DOT_PRODUCT(q, DCOS(hdotr))
!      sin_sum_i = sin_sum(i,ibox) + DOT_PRODUCT(q, DSIN(hdotr))
!
!      E_reciprocal = E_reciprocal + cn(i,ibox) &
!                   * ( cos_sum_i * cos_sum_i &
!                     + sin_sum_i * sin_sum_i )
!
!    END DO
    DO i = 1, nvecs(ibox)
        ihx = hx(i,ibox)
        ihy = hy(i,ibox)
        ihz = hz(i,ibox)
        cos_sum_i = 0.0_DP
        sin_sum_i = 0.0_DP
        DO ia = 1, natoms(is)
                hdotr = ihx*rxp(ia)
                hdotr = hdotr + ihy*ryp(ia)
                hdotr = hdotr + ihz*rzp(ia)
                qia = q(ia)
                cos_sum_i = cos_sum_i + qia*DCOS(hdotr)
                sin_sum_i = sin_sum_i + qia*DSIN(hdotr)
        END DO
        this_cos_sum(i) = cos_sum_i
        this_sin_sum(i) = sin_sum_i
    END DO
    DO i = 1, nvecs(ibox)
        cos_sum_i = this_cos_sum(i) + cos_sum(i,ibox)
        sin_sum_i = this_sin_sum(i) + sin_sum(i,ibox)
        trigfactor = cos_sum_i*cos_sum_i
        trigfactor = trigfactor + sin_sum_i*sin_sum_i
        E_recip_redux = E_recip_redux + cn(i,ibox)*trigfactor
    END DO

    E_reciprocal = E_recip_redux * charge_factor

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
             l_debug_print = .TRUE.
             CALL Compute_Molecule_Nonbond_Inter_Energy(this_im_1,is,vlj_pair,vqq_pair,my_overlap)
             l_debug_print = .FALSE.
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
                l_debug_print = .TRUE.
                CALL Compute_Molecule_Nonbond_Inter_Energy(this_im_1,is_1,vlj_pair,vqq_pair,my_overlap)
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

    rxijp = molecule_1%rcom(1) - molecule_2%rcom(1)
    ryijp = molecule_1%rcom(2) - molecule_2%rcom(2)
    rzijp = molecule_1%rcom(3) - molecule_2%rcom(3)

    ! Compute the minimum image distance

    CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp,rxcom,rycom,rzcom)

    rcom = DSQRT(rxcom*rxcom + rycom*rycom + rzcom*rzcom)

    IF (CBMC_flag) THEN

       rinteraction = rcut_cbmc(this_box) + molecule_1%rcom(4) &
            + molecule_2%rcom(4)

       IF (rcom > rinteraction) get_interaction = .FALSE.

    ELSE

       rinteraction = rcut_max(this_box) + molecule_1%rcom(4) &
            + molecule_2%rcom(4)

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
          rxijp = atom_list(ia,im,is)%rp(1) - atom_list(ja,jm,js)%rp(1)
          ryijp = atom_list(ia,im,is)%rp(2) - atom_list(ja,jm,js)%rp(2)
          rzijp = atom_list(ia,im,is)%rp(3) - atom_list(ja,jm,js)%rp(3)

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

          xcmi = molecule_list(this_locate,is)%rcom(1)
          ycmi = molecule_list(this_locate,is)%rcom(2)
          zcmi = molecule_list(this_locate,is)%rcom(3)

          DO ia = 1, natoms(is)

             piix = atom_list(ia,this_locate,is)%rp(1) - xcmi
             piiy = atom_list(ia,this_locate,is)%rp(2) - ycmi
             piiz = atom_list(ia,this_locate,is)%rp(3) - zcmi
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

                arg = hx(i,this_box)*atom_list(ia,this_locate,is)%rp(1) + &
                      hy(i,this_box)*atom_list(ia,this_locate,is)%rp(2) + &
                      hz(i,this_box)*atom_list(ia,this_locate,is)%rp(3)

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
                hdotr = hx(i,this_box) * atom_list(ia,this_locate,is)%rp(1) + &
                        hy(i,this_box) * atom_list(ia,this_locate,is)%rp(2) + &
                        hz(i,this_box) * atom_list(ia,this_locate,is)%rp(3)

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
