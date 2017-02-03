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

SUBROUTINE Chempot(this_box,is)

  !*****************************************************************************
  !
  ! The subroutine calculates the chemical potential of the species in the system.  
  !
  ! Called by
  !
  !   nvtmc_driver.f90
  !   nptmc_driver.f90
  !   gemc_particle_transfer.f90
  !
  ! Revision history
  !
  !   12/10/13  : Beta Release
  ! 
  !*****************************************************************************

  USE Global_Variables
  USE Energy_Routines
  USE IO_Utilities
  USE Random_Generators
  USE Rotation_Routines
  USE Type_Definitions
  USE Simulation_Properties
  USE Fragment_Growth

  IMPLICIT NONE

  INTEGER :: is, alive, this_box, i_type, i, anchor_dummy

  INTEGER, ALLOCATABLE :: frag_order(:)

  REAL(DP) :: dx, dy, dz 
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq, E_improper, E_periodic_qq
  REAL(DP) :: delta_e, E_reciprocal_move, E_self_move, E_lrc
  REAL(DP) :: prefact, CP_energy, nrg_ring_frag_tot

  REAL(DP) :: ln_pseq, ln_pbias, this_lambda, nrg_ring_frag_out

  LOGICAL :: inter_overlap ,cbmc_overlap, intra_overlap

  delta_e = 0.0_DP
  prefact = 1.0_DP
  nrg_ring_frag_tot = 0.0_DP
  ln_pseq = 0.0_DP
  ln_pbias = 0.0_DP

  ntrials(is,this_box)%cpcalc = ntrials(is,this_box)%cpcalc + 1
 
  
  IF ( locate(max_molecules(is)+1,is,this_box) == 0 ) THEN
     locate(max_molecules(is)+1,is,this_box) = max_molecules(is)+1
  END IF

  alive = locate(max_molecules(is)+1,is,this_box)
  molecule_list(alive,is)%which_box = this_box
  molecule_list(alive,is)%frac = 1.0_DP
  molecule_list(alive,is)%molecule_type = int_normal
  molecule_list(alive,is)%live = .TRUE.

  cbmc_overlap = .FALSE.

  IF(species_list(is)%fragment) THEN

     del_flag = .FALSE.
     get_fragorder = .TRUE.
     this_lambda = molecule_list(alive,is)%frac
     anchor_dummy = 0
     ALLOCATE(frag_order(nfragments(is)))
     CALL Build_Molecule(alive,is,this_box,frag_order,this_lambda, &
             ln_pseq,ln_pbias,nrg_ring_frag_out,cbmc_overlap)
     DEALLOCATE(frag_order)

  ELSE

     molecule_list(alive,is)%xcom = species_list(is)%xcom
     molecule_list(alive,is)%ycom = species_list(is)%ycom
     molecule_list(alive,is)%zcom = species_list(is)%zcom

     atom_list(:,alive,is)%rxp = init_list(:,1,is)%rxp
     atom_list(:,alive,is)%ryp = init_list(:,1,is)%ryp
     atom_list(:,alive,is)%rzp = init_list(:,1,is)%rzp

     atom_list(:,alive,is)%exist = .TRUE.

     CALL Rotate_Molecule_Eulerian(alive,is)

     IF ( box_list(this_box)%int_box_shape == int_cubic ) THEN
 
        molecule_list(alive,is)%xcom = (rranf() - 0.5_DP) * box_list(this_box)%length(1,1)
        molecule_list(alive,is)%ycom = (rranf() - 0.5_DP) * box_list(this_box)%length(2,2)
        molecule_list(alive,is)%zcom = (rranf() - 0.5_DP) * box_list(this_box)%length(3,3)

     END IF

     dx = molecule_list(alive,is)%xcom - species_list(is)%xcom
     dy = molecule_list(alive,is)%ycom - species_list(is)%ycom
     dz = molecule_list(alive,is)%zcom - species_list(is)%zcom

     atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp + dx
     atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp + dy
     atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp + dz

  END IF

  IF (cbmc_overlap) THEN
        ! There is nothing to be done.
        
     molecule_list(alive,is)%live = .FALSE.
     atom_list(:,alive,is)%exist = .FALSE.
     molecule_list(alive,is)%molecule_type = int_none

     RETURN

  END IF
     
  CALL Get_COM(alive,is)

  ! compute the distance of atom farthest from COM

  CALL Compute_Max_COM_Distance(alive,is)

  ! Intra molecule energy

!  CALL Compute_Molecule_Bond_Energy(alive,is,E_bond)
!  CALL Compute_Molecule_Angle_Energy(alive,is,E_angle)
!  CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihedral)
!  CALL Compute_Molecule_Improper_Energy(alive,is,E_improper)

!  delta_e = E_bond + E_angle + E_dihedral + E_improper + delta_e


  ! Nonbonded energy
  CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw,E_intra_qq, &
       E_periodic_qq,intra_overlap)
  CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_inter_vdw,E_inter_qq, &
       inter_overlap)
  E_inter_qq = E_inter_qq + E_periodic_qq

!    delta_e = delta_e + E_inter_vdw + E_inter_qq
  delta_e = delta_e + E_intra_vdw + E_intra_qq + E_inter_vdw + E_inter_qq

  IF (int_charge_style(this_box) == charge_coul) THEN
       IF ( int_charge_sum_style(this_box) == charge_ewald .AND. has_charge(is)) THEN
          CALL Update_System_Ewald_Reciprocal_Energy(alive,is,this_box, &
               int_insertion,E_reciprocal_move)

          delta_e = delta_e + (E_reciprocal_move-energy(this_box)%ewald_reciprocal)

       END IF
       CALL Compute_Molecule_Self_Energy(alive,is,this_box,E_self_move)
       delta_e = delta_e + E_self_move 
  END IF

  IF (int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
     nbeads_in(:) = nint_beads(:,this_box)

     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,this_box) = nint_beads(i_type,this_box) + 1
     END DO
     CALL Compute_LR_correction(this_box,e_lrc)
     delta_e = delta_e + e_lrc - energy(this_box)%lrc
     nint_beads(:,this_box) = nbeads_in(:)
  END IF

  CP_energy = delta_e 

  IF(int_sim_type == sim_npt) THEN
    prefact = box_list(this_box)%volume / REAL(nmols(is,this_box) + 1, DP)
  END IF

!  IF(species_list(is)%fragment) THEN
!     prefact = prefact / (P_forward * kappa ** nfragments(is))
     ! subtract off the angle energy as this was used in biasing
     ! the branch points. Also reduce the total energy by 
     ! ring biasing energy if any
!     CP_energy = CP_energy - E_angle - nrg_ring_frag_tot
!     IF(rx_Flag) CP_energy = CP_energy - E_dihedral
!  END IF

  chpot(is,this_box) = chpot(is,this_box) + prefact *  DEXP(-beta(this_box) * CP_energy)

  molecule_list(alive,is)%live = .FALSE.
  atom_list(:,alive,is)%exist = .FALSE.
  molecule_list(alive,is)%molecule_type = int_none

  IF ( int_charge_sum_style(this_box) == charge_ewald .AND. has_charge(is)) THEN
        ! restore cos_sum and sin_sum. Note that these were changed when difference in
        ! reciprocal energies was computed
     !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
     cos_sum(:,this_box) = cos_sum_old(:,this_box)
     sin_sum(:,this_box) = sin_sum_old(:,this_box)
     !$OMP END PARALLEL WORKSHARE
  END IF


  RETURN

END SUBROUTINE Chempot    

