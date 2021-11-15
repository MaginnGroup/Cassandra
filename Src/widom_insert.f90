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

SUBROUTINE Widom_Insert(is,ibox,widom_sum)

  !*****************************************************************************
  ! 
  ! PURPOSE: Perform all Widom insertions for species is and box ibox for the
  !          current step and return widom_sum.
  !
  ! Called by
  !
  !    Widom_Subdriver
  !
  ! 
  !*****************************************************************************

  USE Global_Variables
  USE Energy_Routines
  USE IO_Utilities
  USE Random_Generators
  USE Rotation_Routines
  USE Fragment_Growth
  !$ USE OMP_LIB

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN) :: ibox ! insert test particles in ibox
  INTEGER, INTENT(IN) :: is ! species indices

  ! Local declarations
  INTEGER :: i, i_type               ! atom indices
  INTEGER :: im                      ! molecule INDEX
  INTEGER :: frag_order(nfragments(is))

  INTEGER (KIND=INT64) :: i_widom
  INTEGER (KIND=INT64) :: insertions_in_step, n_overlaps

  REAL(DP) :: dx, dy, dz
  REAL(DP) :: dE, dE_intra, dE_inter, dE_frag
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq, E_periodic_qq
  REAL(DP) :: E_reciprocal, E_self, E_lrc
  REAL(DP) :: E_ring_frag
  REAL(DP) :: ln_pacc, ln_pseq, ln_pbias, this_lambda

  REAL(DP) :: widom_prefactor, widom_var_exp, widom_sum
  REAL(DP) :: E_recip_in, lrc_diff, E_inter_constant


  LOGICAL :: inter_overlap, cbmc_overlap, intra_overlap


  this_lambda = 1.0_DP
  widom_sum = 0.0_DP
  n_overlaps = 0_INT64

  lrc_diff = 0.0_DP
  E_self = 0.0_DP
  E_recip_in = 0.0_DP

  del_flag = .FALSE.
  get_fragorder = .TRUE.
  
  nmols(is,ibox) = nmols(is,ibox)+1
  im = nmols(is,ibox)
  locate(im,is,ibox) = locate(nmols(is,0),is,0)
  locate(nmols(is,0),is,0) = 0

  !  * Set properties of the to-be-inserted molecule
  widom_species = is
  widom_locate = locate(im,is,ibox)
  molecule_list(widom_locate,is)%which_box = ibox
  molecule_list(widom_locate,is)%frac = this_lambda
  molecule_list(widom_locate,is)%molecule_type = int_normal
  
  widom_prefactor = box_list(ibox)%volume&
                  / (REAL(nmols(is,ibox),DP)*((species_list(is)%de_broglie(ibox))**3))
  insertions_in_step = species_list(is)%insertions_in_step(ibox)
  
  ! Long-range energy correction

  IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN

     ! increase number of integer beads
     nbeads_in = nint_beads(:,ibox)

     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,ibox) = nint_beads(i_type,ibox) + 1
     END DO

     CALL Compute_LR_correction(ibox,E_lrc)
     nint_beads(:,ibox) = nbeads_in
     lrc_diff = E_lrc - energy(ibox)%lrc

  END IF
  
  IF (int_charge_style(ibox) == charge_coul) THEN
          CALL Compute_Molecule_Self_Energy(widom_locate,is,ibox,E_self)
          E_recip_in = energy(ibox)%reciprocal
  END IF
  E_inter_constant = lrc_diff + E_self + E_recip_in

  widom_active = .TRUE.

  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(ln_pseq, ln_pbias, E_ring_frag, inter_overlap, cbmc_overlap, intra_overlap) &
  !$OMP PRIVATE(widom_var_exp, E_inter_qq, E_periodic_qq, E_intra_qq, E_intra_vdw, E_inter_vdw) &
  !$OMP PRIVATE(E_bond, E_angle, E_dihedral, E_improper, dE_intra, dE_inter, E_reciprocal, frag_order) &
  !$OMP REDUCTION(+:widom_sum,n_overlaps)
  IF (ALLOCATED(widom_atoms)) DEALLOCATE(widom_atoms)
  ALLOCATE(widom_atoms(natoms(is)))
  widom_molecule = molecule_list(widom_locate,is)
  widom_atoms = atom_list(1:natoms(is),widom_locate,is)
  widom_sum = 0.0_DP
  n_overlaps = 0_INT64


  !$OMP DO SCHEDULE(DYNAMIC)
  !$
  DO i_widom = 1, insertions_in_step
          ! Initialize variables
          ln_pseq = 0.0_DP
          ln_pbias = 0.0_DP
          E_ring_frag = 0.0_DP
          inter_overlap = .FALSE.
          cbmc_overlap = .FALSE.
          intra_overlap = .FALSE.
          widom_var_exp = 0.0_DP
          ! Now that an insertion will be attempted, we need to do some bookkeeping:

          !  * Increment the counters to track number of widom insertions

         
          !*****************************************************************************
          ! Choose a position, orientation and conformation for the 
          !         to-be-inserted molecule
          !*****************************************************************************
          !
          ! Build_Molecule places the first fragment, then calls Fragment_Placement
          ! to place the additional fragments 
          CALL Build_Molecule(widom_locate,is,ibox,frag_order,this_lambda, &
                  ln_pseq,ln_pbias,E_ring_frag,cbmc_overlap)

          ! Turn the molecule on
          widom_molecule%live = .TRUE.
          widom_atoms%exist = .TRUE.

          ! So far ln_pbias includes 
          !   * the probability of choosing the insertion 
          !     point from the collection of trial coordinates 
          !   * the probability of choosing each dihedral from the collection of trial dihedrals. 
          ! Now add
          !   * the probability of the fragment sequence, ln_pseq
          !   * the number of trial coordinates, kappa_ins
          !   * the number of trial dihedrals, kappa_dih, for each dihedral.

          ln_pbias = ln_pbias + ln_pseq
          ln_pbias = ln_pbias + DLOG(REAL(kappa_ins,DP))

          IF (kappa_rot /= 0 ) THEN
             ln_pbias = ln_pbias + DLOG(REAL(kappa_rot,DP))
          END IF

          IF (kappa_dih /= 0 ) THEN
             ln_pbias = ln_pbias + REAL(nfragments(is)-1,DP) * DLOG(REAL(kappa_dih,DP))
          END IF

          IF (.NOT. cbmc_overlap) THEN

            ! Molecule COM may be outside the box boundary if grown via CBMC, so wrap
            ! the molecule coordinates back in the box (if needed)
            CALL Fold_Molecule(widom_locate,is,ibox)

            ! Recompute the COM in case the molecule was wrapped
            !CALL Get_COM(widom_locate,is)

            ! Compute the distance of the atom farthest from COM
            !CALL Compute_Max_COM_Distance(widom_locate,is)

            ! Calculate the potential energy interaction between the inserted molecule
            ! and the rest of the system
            CALL Compute_Molecule_Nonbond_Inter_Energy_Widom(widom_locate,is, &
                    E_inter_vdw,E_inter_qq,inter_overlap)

            ! Calculate the nonbonded energy interaction within the inserted molecule
            CALL Compute_Molecule_Nonbond_Intra_Energy(widom_locate,is, &
                    E_intra_vdw,E_intra_qq,E_periodic_qq,intra_overlap)
            E_inter_qq = E_inter_qq + E_periodic_qq
         
          END IF

          ! Leave widom_sum unchanged if there is any core overlap
          IF (.NOT. (cbmc_overlap .OR. inter_overlap .OR. intra_overlap)) THEN
                  ! There are no overlaps, so we can calculate the change in potential energy.
                  !
                  ! Already have the change in nonbonded energies
                  dE_inter = E_inter_vdw + E_inter_qq + E_inter_constant 
                  dE_intra = E_intra_vdw + E_intra_qq

                  ! Bonded intramolecular energies
                  ! If the molecule was grown via CBMC, we already have the intramolecular 
                  ! bond energies? Otherwise we need to compute them.
                  CALL Compute_Molecule_Bond_Energy(widom_locate,is,E_bond)
                  CALL Compute_Molecule_Angle_Energy(widom_locate,is,E_angle)
                  CALL Compute_Molecule_Dihedral_Energy(widom_locate,is,E_dihedral)
                  CALL Compute_Molecule_Improper_Energy(widom_locate,is,E_improper)

                  dE_intra = dE_intra + E_bond + E_angle + E_dihedral + E_improper

                  ! Ewald energies
                  IF (int_charge_style(ibox) == charge_coul) THEN
                        IF ( (int_charge_sum_style(ibox) == charge_ewald) .AND. &
                             has_charge(is) ) THEN
                       
                            CALL Update_System_Ewald_Reciprocal_Energy_Widom(widom_locate, &
                                   is,ibox,E_reciprocal)

                            dE_inter = dE_inter + E_reciprocal
                        END IF
                  END IF

                  ! moved to before the loop
                  !widom_prefactor = box_list(ibox)%volume&
                  !        / (REAL(nmols(is,ibox),DP)*((species_list(is)%de_broglie(ibox))**3))

                  ! change in energy, less energy used to bias fragment selection
                  dE = dE_intra + dE_inter
                  dE_frag = E_angle + E_ring_frag

                  ! mu' = -(1/beta)*ln(<widom_var>)
                  widom_var_exp = DEXP(-beta(ibox) * (dE - dE_frag) - ln_pbias)
                  ! sum of all widom_var for this step; output argument
                  widom_sum = widom_sum + widom_var_exp
          ELSE
                  n_overlaps = n_overlaps + 1_INT64
          END IF
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  widom_active = .FALSE.
  widom_sum = widom_sum * widom_prefactor
  overlap_counter(is,ibox) = overlap_counter(is,ibox) + n_overlaps


  ! remove test molecule
  nmols(is,ibox) = nmols(is,ibox)-1
  locate(im,is,ibox) = 0
  molecule_list(widom_locate,is)%live = .FALSE.
  atom_list(:,widom_locate,is)%exist = .FALSE.
  molecule_list(widom_locate,is)%molecule_type = int_none

  ! move locate to the list of unused locates
  locate(nmols(is,0),is,0) = widom_locate
  widom_locate = 0
  widom_species = 0


  ntrials(is,ibox)%widom = ntrials(is,ibox)%widom + insertions_in_step

END SUBROUTINE Widom_Insert
!*******************************************************************************
