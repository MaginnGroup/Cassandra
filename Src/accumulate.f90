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

SUBROUTINE Accumulate(ibox)

!*******************************************************************************
!
! This subroutine accumulates various thermodynamics properties
!
! Called by
!
! nvtmc_driver
! nptmc_driver
!
! Revision History:
!
! 12/10/13  : Beta Release
!
!********************************************************************************
  USE Global_Variables
  USE Simulation_Properties

  IMPLICIT NONE

  INTEGER :: is, ibox, nmolecules_is
  REAL(DP) :: mass_density
  
  ! want 1 ... block_avg_freq to give iblock = 1
  iblock = (i_mcstep - initial_mcstep - 1) / block_avg_freq + 1

  !--- energy accumulators

  ac_energy(ibox,iblock)%total = ac_energy(ibox,iblock)%total + energy(ibox)%total
  ac_energy(ibox,iblock)%inter_vdw = ac_energy(ibox,iblock)%inter_vdw + energy(ibox)%inter_vdw
  ac_energy(ibox,iblock)%inter_q   = ac_energy(ibox,iblock)%inter_q   + energy(ibox)%inter_q
  ac_energy(ibox,iblock)%intra_vdw = ac_energy(ibox,iblock)%intra_vdw + energy(ibox)%intra_vdw
  ac_energy(ibox,iblock)%intra_q   = ac_energy(ibox,iblock)%intra_q   + energy(ibox)%intra_q
  ac_energy(ibox,iblock)%intra     = ac_energy(ibox,iblock)%intra + energy(ibox)%intra
  ac_energy(ibox,iblock)%ewald_reciprocal = ac_energy(ibox,iblock)%ewald_reciprocal + energy(ibox)%ewald_reciprocal
  ac_energy(ibox,iblock)%self = ac_energy(ibox,iblock)%self + energy(ibox)%self

  IF(int_vdw_sum_style(ibox) == vdw_cut_tail) THEN
     ac_energy(ibox,iblock)%lrc = ac_energy(ibox,iblock)%lrc + energy(ibox)%lrc
  END IF

  !--- virial accumulators

!  ac_virial(ibox,iblock)%total = ac_virial(ibox,iblock)%total + virial(ibox)%total
!  ac_virial(ibox,iblock)%inter_vdw = ac_virial(ibox,iblock)%inter_vdw + virial(ibox)%inter_vdw
!  ac_virial(ibox,iblock)%inter_q   = ac_virial(ibox,iblock)%inter_q + virial(ibox)%inter_q
!  ac_virial(ibox,iblock)%intra_vdw = ac_virial(ibox,iblock)%intra_vdw + virial(ibox)%intra_vdw
!  ac_virial(ibox,iblock)%intra_q   = ac_virial(ibox,iblock)%intra_q + virial(ibox)%intra_q
!  ac_virial(ibox,iblock)%ewald_reciprocal = ac_virial(ibox,iblock)%ewald_reciprocal + virial(ibox)%ewald_reciprocal

  !--- thermodynamic accumulators
  ac_volume(ibox,iblock) = ac_volume(ibox,iblock) + box_list(ibox)%volume
  IF (need_pressure) THEN
     IF (pressure(ibox)%last_calc /= i_mcstep) THEN
        pressure(ibox)%last_calc = i_mcstep
        CALL Compute_Pressure(ibox)
     END IF
     ac_pressure(ibox,iblock) = ac_pressure(ibox,iblock) + pressure(ibox)%computed
     IF (int_sim_type /= sim_npt .AND. int_sim_type /= sim_gemc_npt) THEN
        ac_enthalpy(ibox,iblock) = ac_enthalpy(ibox,iblock) &
                          + energy(ibox)%total + pressure(ibox)%computed * box_list(ibox)%volume
     END IF
  END IF
  IF (int_sim_type == sim_npt .OR. int_sim_type == sim_gemc_npt) THEN
     ac_enthalpy(ibox,iblock) = ac_enthalpy(ibox,iblock) &
                       + energy(ibox)%total + pressure(ibox)%setpoint * box_list(ibox)%volume
  END IF

  !--- particle density
  
  DO is = 1, nspecies
     ac_density(is,ibox,iblock) = ac_density(is,ibox,iblock) + REAL(nmols(is,ibox),DP) / box_list(ibox)%volume
     ac_nmols(is,ibox,iblock) = ac_nmols(is,ibox,iblock) + REAL(nmols(is,ibox),DP)
  END DO

  mass_density = 0.0_DP
  DO is = 1, nspecies
     mass_density = mass_density &
                  + REAL(nmols(is,ibox),DP) * species_list(is)%molecular_weight
  END DO
  mass_density = mass_density / box_list(ibox)%volume
  ac_mass_density(ibox,iblock) = ac_mass_density(ibox,iblock) + mass_density

END SUBROUTINE Accumulate
