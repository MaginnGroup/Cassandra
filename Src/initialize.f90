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
!    initialize - subroutine initializes the counters for the simulation
!
!    init_accumulators - subroutine initializes the accumulators for the
!    simulation
!
! Revision history
!
!   12/10/13 : Beta version
!********************************************************************************
SUBROUTINE Initialize
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
  USE Global_Variables
 
  IMPLICIT NONE

  ALLOCATE(nsuccess(nspecies,nbr_boxes))
  ! Allocate the memory for volume counters in NPT simulation

  ALLOCATE(nvolumes(nbr_boxes))
  ALLOCATE(nvol_success(nbr_boxes))
  ALLOCATE(ivol_success(nbr_boxes))
  ALLOCATE(chpot(nspecies,nbr_boxes))
  ALLOCATE(chpotid(nspecies,nbr_boxes))

  chpot(:,:) = 0.0_DP
  chpotid(:,:) = 0.0_DP
  nvolumes(:) = 0
  nvol_success(:) = 0
  ivol_success(:) = 0
  nsuccess(:,:)%displacement = 0
  nsuccess(:,:)%rotation = 0
  nsuccess(:,:)%displacement_e = 0
  nsuccess(:,:)%rotation_e = 0
  nsuccess(:,:)%dihedral = 0
  nsuccess(:,:)%angle = 0
  nsuccess(:,:)%insertion = 0
  nsuccess(:,:)%deletion = 0
  nsuccess(:,:)%disp_atom = 0

  ! define success counters for equilibration runs
  IF (.NOT. ALLOCATED(nequil_success)) ALLOCATE(nequil_success(nspecies,nbr_boxes))
  nequil_success(:,:)%displacement = 0
  nequil_success(:,:)%rotation = 0
!  nvol_equil_success(:) = 0
 
  IF ( .NOT. ALLOCATED(ntrials)) ALLOCATE(ntrials(nspecies,nbr_boxes))
  IF ( .NOT. ALLOCATED(tot_trials)) ALLOCATE(tot_trials(nbr_boxes))
 
  ntrials(:,:)%displacement = 0
  ntrials(:,:)%rotation = 0
  ntrials(:,:)%dihedral = 0
  ntrials(:,:)%angle = 0
  ntrials(:,:)%insertion = 0
  ntrials(:,:)%deletion = 0
  ntrials(:,:)%disp_atom = 0
  ntrials(:,:)%cpcalc = 0
 
  IF ( SUM(nfragments) > 0 ) THEN
     IF ( .NOT. ALLOCATED(regrowth_trials))ALLOCATE(regrowth_trials(MAXVAL(nfragments),nspecies))
     IF ( .NOT. ALLOCATED(regrowth_success)) ALLOCATE(regrowth_success(MAXVAL(nfragments),nspecies))
     regrowth_trials(:,:) = 0
     regrowth_success(:,:)= 0
  END IF


  tot_trials(:) = 0
 
  igas_flag = .FALSE.

END SUBROUTINE Initialize

SUBROUTINE Init_Accumulators
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
  ! 06/22/16 (RGM) : Created beta version
  !*********************************************************************
  USE Global_Variables
 
  IMPLICIT NONE

  ! Initialize accumulators 
  IF (block_average) THEN
     nbr_blocks = (n_mcsteps - initial_mcstep - 1) / block_avg_freq + 1

     ALLOCATE(ac_energy(nbr_boxes,nbr_blocks))
!     ALLOCATE(ac_virial(nbr_boxes,nbr_blocks))
     ALLOCATE(ac_density(nspecies,nbr_boxes,nbr_blocks))
     ALLOCATE(ac_nmols(nspecies,nbr_boxes,nbr_blocks))
     ALLOCATE(ac_volume(nbr_boxes,nbr_blocks))
     ALLOCATE(ac_pressure(nbr_boxes,nbr_blocks))
     ALLOCATE(ac_enthalpy(nbr_boxes,nbr_blocks))
     ALLOCATE(ac_mass_density(nbr_boxes,nbr_blocks))

     ac_energy(:,:)%total = 0.0_DP
     ac_energy(:,:)%intra = 0.0_DP
     ac_energy(:,:)%intra_vdw = 0.0_DP
     ac_energy(:,:)%intra_q = 0.0_DP
     ac_energy(:,:)%inter_vdw = 0.0_DP
     ac_energy(:,:)%lrc = 0.0_DP
     ac_energy(:,:)%inter_q = 0.0_DP
     ac_energy(:,:)%ewald_reciprocal = 0.0_DP
     ac_energy(:,:)%self = 0.0_DP

!     ac_virial(:,:)%total = 0.0_DP
!     ac_virial(:,:)%intra = 0.0_DP
!     ac_virial(:,:)%intra_vdw = 0.0_DP
!     ac_virial(:,:)%intra_q = 0.0_DP
!     ac_virial(:,:)%inter_vdw = 0.0_DP
!     ac_virial(:,:)%lrc = 0.0_DP
!     ac_virial(:,:)%inter_q = 0.0_DP
!     ac_virial(:,:)%ewald_reciprocal = 0.0_DP
!     ac_virial(:,:)%self = 0.0_DP

     ac_density = 0.0_DP
     ac_nmols = 0.0_DP
     ac_volume = 0.0_DP
     ac_pressure = 0.0_DP
     ac_enthalpy = 0.0_DP
     ac_mass_density = 0.0_DP
  END IF

END SUBROUTINE Init_Accumulators
