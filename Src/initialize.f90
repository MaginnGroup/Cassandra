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
!*******************************************************************************
!
!
! Contains two routines
!
!    initialize - subroutine initializes the accumulators for a given box.
!
!    reset - resets the accumulators
!
! Revision history
!
!   12/10/13 : Beta version
!********************************************************************************
SUBROUTINE Initialize(this_box)
  !**********************************************************************
  !
  ! This subroutine initializes the accumulators for a given box.
  !
  ! CALLED BY
  !
  !        main
  ! CALLS
  !        None
  !
  ! 08/06/13 (JS) : Created beta version
  !*********************************************************************
  USE Run_Variables
 
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box

  IF ( .NOT. ALLOCATED(ac_energy)) ALLOCATE(ac_energy(nbr_boxes))
  IF ( .NOT. ALLOCATED(ac_virial)) ALLOCATE(ac_virial(nbr_boxes))
  IF ( .NOT. ALLOCATED(ac_density)) ALLOCATE(ac_density(nspecies,nbr_boxes))
  IF ( .NOT. ALLOCATED(ac_nmols)) ALLOCATE(ac_nmols(nspecies,nbr_boxes))
  IF ( .NOT. ALLOCATED(ac_volume)) ALLOCATE(ac_volume(nbr_boxes))

  IF ( .NOT. ALLOCATED(nsuccess)) ALLOCATE(nsuccess(nspecies,nbr_boxes))
  ! Allocate the memory for volume counters in NPT simulation

  IF ( .NOT. ALLOCATED(nvolumes) ) ALLOCATE(nvolumes(nbr_boxes))
  IF ( .NOT. ALLOCATED(nvol_success)) ALLOCATE(nvol_success(nbr_boxes))
  IF ( .NOT. ALLOCATED(ivol_success)) ALLOCATE(ivol_success(nbr_boxes))
  IF ( .NOT. ALLOCATED(ac_enthalpy)) ALLOCATE(ac_enthalpy(nbr_boxes))
  IF ( .NOT. ALLOCATED(chpot)) ALLOCATE(chpot(nspecies,nbr_boxes))
  IF ( .NOT. ALLOCATED(chpotid)) ALLOCATE(chpotid(nspecies,nbr_boxes))

  chpot(:,:) = 0.0_DP
  chpotid(:,:) = 0.0_DP
  nvolumes(:) = 0
  nvol_success(:) = 0
  ivol_success(:) = 0
  nsuccess(:,this_box)%displacement = 0
  nsuccess(:,this_box)%rotation = 0
  nsuccess(:,this_box)%displacement_e = 0
  nsuccess(:,this_box)%rotation_e = 0
  nsuccess(:,this_box)%dihedral = 0
  nsuccess(:,this_box)%angle = 0
  nsuccess(:,this_box)%insertion = 0
  nsuccess(:,this_box)%deletion = 0
  nsuccess(:,this_box)%disp_atom = 0

  ! define success counters for equilibration runs
  IF (.NOT. ALLOCATED(nequil_success)) ALLOCATE(nequil_success(nspecies,nbr_boxes))
  nequil_success(:,this_box)%displacement = 0
  nequil_success(:,this_box)%rotation = 0
!  nvol_equil_success(:) = 0
 
  IF ( .NOT. ALLOCATED(ntrials)) ALLOCATE(ntrials(nspecies,nbr_boxes))
  IF ( .NOT. ALLOCATED(tot_trials)) ALLOCATE(tot_trials(nbr_boxes))
 
  ntrials(:,this_box)%displacement = 0
  ntrials(:,this_box)%rotation = 0
  ntrials(:,this_box)%dihedral = 0
  ntrials(:,this_box)%angle = 0
  ntrials(:,this_box)%insertion = 0
  ntrials(:,this_box)%deletion = 0
  ntrials(:,this_box)%disp_atom = 0
  ntrials(:,this_box)%cpcalc = 0
 
  IF ( SUM(nfragments) > 0 ) THEN
     IF ( .NOT. ALLOCATED(regrowth_trials))ALLOCATE(regrowth_trials(MAXVAL(nfragments),nspecies))
     IF ( .NOT. ALLOCATED(regrowth_success)) ALLOCATE(regrowth_success(MAXVAL(nfragments),nspecies))
     regrowth_trials(:,:) = 0
     regrowth_success(:,:)= 0
  END IF


  tot_trials(this_box) = 0
 
  ac_density(:,:) = 0.0_DP
  igas_flag = .FALSE.

END SUBROUTINE Initialize

SUBROUTINE Reset(this_box)
  !*********************************************************************
  ! The subroutine resets all the accumulators to zero. 
  !
  ! CALLED BY
  !
  !        gcmc_driver
  !        gemc_driver
  !        main
  !        nptmc_driver
  !        nvtmc_driver
  ! Main
  !
  ! Written by Jindal Shah on 01/09/08
  !
  !*********************************************************************
  USE Run_Variables

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box
  
  ac_energy(this_box)%total = 0.0_DP
  ac_energy(this_box)%intra = 0.0_DP
  ac_energy(this_box)%intra_vdw = 0.0_DP
  ac_energy(this_box)%intra_q = 0.0_DP
  ac_energy(this_box)%inter_vdw = 0.0_DP
  ac_energy(this_box)%lrc = 0.0_DP
  ac_energy(this_box)%inter_q = 0.0_DP
  ac_energy(this_box)%ewald_reciprocal = 0.0_DP
  ac_energy(this_box)%ewald_self = 0.0_DP


  ac_virial(this_box)%total = 0.0_DP
  ac_virial(this_box)%intra = 0.0_DP
  ac_virial(this_box)%intra_vdw = 0.0_DP
  ac_virial(this_box)%intra_q = 0.0_DP
  ac_virial(this_box)%inter_vdw = 0.0_DP
  ac_virial(this_box)%lrc = 0.0_DP
  ac_virial(this_box)%inter_q = 0.0_DP
  ac_virial(this_box)%ewald_reciprocal = 0.0_DP
  ac_virial(this_box)%ewald_self = 0.0_DP

  ac_density(:,this_box) = 0.0_DP
  ac_nmols(:,this_box) = 0.0_DP

  ac_volume(this_box) = 0.0_DP
  ac_enthalpy(this_box) = 0.0_DP
  
END SUBROUTINE Reset
  
