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

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  IMPLICIT NONE

  ! Arguments
  INTEGER :: ibox ! insert test particles in ibox

  ! Local declarations
  INTEGER :: i, i_type               ! atom indices
  INTEGER :: im                      ! molecule INDEX
  INTEGER :: lm                      ! molecule LOCATE
  INTEGER :: is ! species indices
  INTEGER, ALLOCATABLE :: frag_order(:)

  INTEGER :: i_widom
  INTEGER :: insertions_in_step

  REAL(DP) :: dx, dy, dz
  REAL(DP) :: dE, dE_intra, dE_inter, dE_frag
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq, E_periodic_qq
  REAL(DP) :: E_reciprocal, E_self, E_lrc
  REAL(DP) :: E_ring_frag
  REAL(DP) :: ln_pacc, ln_pseq, ln_pbias, this_lambda

  REAL(DP) :: widom_prefactor, widom_var, widom_sum


  LOGICAL :: inter_overlap, cbmc_overlap, intra_overlap


  this_lambda = 1.0_DP
  widom_sum = 0.0_DP
  
  nmols(is,ibox) = nmols(is,ibox)+1
  im = nmols(is,ibox)
  locate(im,is,ibox) = locate(nmols(is,0),is,0)
  locate(nmols(is,0),is,0) = 0

  !  * Set properties of the to-be-inserted molecule
  lm = locate(im,is,ibox)
  molecule_list(lm,is)%which_box = ibox
  molecule_list(lm,is)%frac = this_lambda
  molecule_list(lm,is)%molecule_type = int_normal
  
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

  END IF

  DO i_widom = 1, insertions_in_step
          ! Initialize variables
          ln_pseq = 0.0_DP
          ln_pbias = 0.0_DP
          E_ring_frag = 0.0_DP
          inter_overlap = .FALSE.
          cbmc_overlap = .FALSE.
          intra_overlap = .FALSE.
          widom_var = 0.0_DP
          ! Now that an insertion will be attempted, we need to do some bookkeeping:

          !  * Increment the counters to track number of widom insertions
          ntrials(is,ibox)%widom = ntrials(is,ibox)%widom + 1

         
          !*****************************************************************************
          ! Choose a position, orientation and conformation for the 
          !         to-be-inserted molecule
          !*****************************************************************************
          !
          ! Build_Molecule places the first fragment, then calls Fragment_Placement
          ! to place the additional fragments 
          del_flag = .FALSE.     ! Change the coordinates of 'lm'
          get_fragorder = .TRUE.
          ALLOCATE(frag_order(nfragments(is)))
          CALL Build_Molecule(lm,is,ibox,frag_order,this_lambda, &
                  ln_pseq,ln_pbias,E_ring_frag,cbmc_overlap)
          DEALLOCATE(frag_order)

          ! Turn the molecule on
          molecule_list(lm,is)%live = .TRUE.
          atom_list(:,lm,is)%exist = .TRUE.

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
            CALL Fold_Molecule(lm,is,ibox)

            ! Recompute the COM in case the molecule was wrapped
            CALL Get_COM(lm,is)

            ! Compute the distance of the atom farthest from COM
            CALL Compute_Max_COM_Distance(lm,is)

            ! Calculate the potential energy interaction between the inserted molecule
            ! and the rest of the system
            CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is, &
                    E_inter_vdw,E_inter_qq,inter_overlap)

            ! Calculate the nonbonded energy interaction within the inserted molecule
            CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is, &
                    E_intra_vdw,E_intra_qq,E_periodic_qq,intra_overlap)
            E_inter_qq = E_inter_qq + E_periodic_qq
         
          END IF

          ! Leave widom_sum unchanged if there is any core overlap
          IF (.NOT. (cbmc_overlap .OR. inter_overlap .OR. intra_overlap)) THEN
                  ! There are no overlaps, so we can calculate the change in potential energy.
                  !
                  ! Already have the change in nonbonded energies
                  dE_inter = E_inter_vdw + E_inter_qq 
                  dE_intra = E_intra_vdw + E_intra_qq

                  ! Bonded intramolecular energies
                  ! If the molecule was grown via CBMC, we already have the intramolecular 
                  ! bond energies? Otherwise we need to compute them.
                  CALL Compute_Molecule_Bond_Energy(lm,is,E_bond)
                  CALL Compute_Molecule_Angle_Energy(lm,is,E_angle)
                  CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihedral)
                  CALL Compute_Molecule_Improper_Energy(lm,is,E_improper)

                  dE_intra = dE_intra + E_bond + E_angle + E_dihedral + E_improper

                  ! Ewald energies
                  IF (int_charge_style(ibox) == charge_coul) THEN
                        IF ( (int_charge_sum_style(ibox) == charge_ewald) .AND. &
                             has_charge(is) ) THEN
                       
                           CALL Update_System_Ewald_Reciprocal_Energy(lm,is,ibox, &
                                   int_insertion,E_reciprocal)

                            dE_inter = dE_inter + (E_reciprocal - energy(ibox)%reciprocal)
                        END IF

                        CALL Compute_Molecule_Self_Energy(lm,is,ibox,E_self)

                        dE_inter = dE_inter + E_self

                  END IF

                  ! Long-range energy correction

                  IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN

                     dE_inter = dE_inter + E_lrc - energy(ibox)%lrc

                  END IF

                  ! moved to before the loop
                  !widom_prefactor = box_list(ibox)%volume&
                  !        / (REAL(nmols(is,ibox),DP)*((species_list(is)%de_broglie(ibox))**3))

                  ! change in energy, less energy used to bias fragment selection
                  dE = dE_intra + dE_inter
                  dE_frag = E_angle + E_ring_frag

                  ! mu' = -(1/beta)*ln(<widom_var>)
                  widom_var = widom_prefactor*DEXP(-beta(ibox) * (dE - dE_frag) - ln_pbias)
                  ! sum of all widom_var for this step; output argument
                  widom_sum = widom_sum + widom_var

                  IF ( int_charge_sum_style(ibox) == charge_ewald .AND. &
                       has_charge(is) ) THEN
                     ! Restore cos_sum and sin_sum. Note that these were changed when the
                     ! difference in reciprocal energies was computed.
                     cos_sum(:,ibox) = cos_sum_old(:,ibox)
                     sin_sum(:,ibox) = sin_sum_old(:,ibox)
                  END IF

          ELSE
                  overlap_counter(is,ibox) = overlap_counter(is,ibox) + 1_INT64
          END IF
  END DO


  IF ( int_vdw_sum_style(ibox) == vdw_cut_tail ) THEN
     ! Restore the total number of bead types
     nint_beads(:,ibox) = nbeads_in(:)
  END IF

  ! remove test molecule
  nmols(is,ibox) = nmols(is,ibox)-1
  locate(im,is,ibox) = 0
  molecule_list(lm,is)%live = .FALSE.
  atom_list(:,lm,is)%exist = .FALSE.
  molecule_list(lm,is)%molecule_type = int_none

  ! move locate to the list of unused locates
  locate(nmols(is,0),is,0) = lm

END SUBROUTINE Widom_Insert
!*******************************************************************************
