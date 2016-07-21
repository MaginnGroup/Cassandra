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

!*******************************************************************************
SUBROUTINE Translate
!*******************************************************************************

  !*****************************************************************************
  ! This subroutine performs a COM translation move for a randomly chosen 
  ! molecule.
  !
  ! 
  ! CALLED BY 
  !
  !        gcmc_driver
  !        gemc_driver
  !        nptmc_driver
  !        nvtmc_driver
  ! 
  ! CALLS
  !     
  !        Get_Index_Molecule
  !        Store_Molecule_Pair_Interaction_Arrays
  !        Compute_Molecule_Nonbond_Inter_Energy
  !        Save_Old_Cartesian_Coordinates
  !        Revert_Old_Cartesian_Coordinates
  !        Reset_Molecule_Pair_Interaction_Arrays
  !        Get_Position_Alive
  !        Update_Ewald_Reciprocal_Energy
  !
  !
  !  12/10/13 : Beta release
  !*****************************************************************************

  USE Type_Definitions
  USE Global_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE Pair_Nrg_Routines
  USE IO_Utilities

  IMPLICIT NONE

  ! Arguments
  INTEGER  :: ibox   ! box index
  INTEGER  :: mc_step

  ! Local declarations
  INTEGER  :: is         ! species index
  INTEGER  :: im, lm     ! molecule indices
  INTEGER  :: total_mols ! number of molecules in the system
  INTEGER  :: nmols_box(nbr_boxes)

  REAL(DP) :: x_box(nbr_boxes), x_species(nspecies)
  REAL(DP) :: randno
  REAL(DP) :: dx, dy, dz
  REAL(DP) :: ln_pacc, success_ratio
  REAL(DP) :: dE, E_vdw, E_qq, E_vdw_move, E_qq_move, E_reciprocal_move
  REAL(DP) :: rcut_small

  LOGICAL :: inter_overlap, overlap, accept_or_reject

  ! Pair_Energy arrays and Ewald implementation
  INTEGER :: position
  REAL(DP), ALLOCATABLE :: cos_mol_old(:), sin_mol_old(:)

! Done with that section

  E_vdw_move = 0.0_DP
  E_qq_move = 0.0_DP
  E_vdw = 0.0_DP
  E_qq = 0.0_DP
  E_reciprocal_move = 0.0_DP
  inter_overlap = .FALSE.
  accept = .FALSE.

  ! Sum the total number of molecules 
  total_mols = 0 ! sum over species, box
  DO ibox = 1, nbr_boxes
    nmols_box(ibox) = 0
    DO is = 1, nspecies
      ! Only count mobile species
      IF ( max_disp(is,ibox) > 0. ) THEN
        total_mols = total_mols + nmols(is,ibox)
        nmols_box(ibox) = nmols_box(ibox) + nmols(is,ibox)
      END IF
    END DO
  END DO

  ! If there are no molecules then return
  IF (total_mols == 0) RETURN

  ! If needed, choose a box based on its total mol fraction
  IF(nbr_boxes .GT. 1) THEN

    DO ibox = 1, nbr_boxes
       x_box(ibox) = REAL(nmols_box(ibox),DP)/REAL(total_mols,DP)
       IF ( ibox > 1 ) THEN
          x_box(ibox) = x_box(ibox) + x_box(ibox-1)
       END IF
    END DO
  
    randno = rranf()
    DO ibox = 1, nbr_boxes
       IF ( randno <= x_box(ibox)) EXIT
    END DO

  ELSE

    ibox = 1

  END IF

  ! If there are no molecules in this box then return
  IF( nmols_box(ibox) == 0 ) RETURN

  ! Choose species based on the mol fraction, using Golden sampling
  DO is = 1, nspecies
     IF ( max_disp(is,ibox) > 0. ) THEN
       x_species(is) = REAL(nmols(is,ibox), DP)/REAL(nmols_box(ibox),DP)
     ELSE
       x_species(is) = 0.0_DP
     END IF
     IF ( is > 1 ) THEN
        x_species(is) = x_species(is) + x_species(is-1)
     END IF
  END DO

  randno = rranf()
  DO is = 1, nspecies
     IF( randno <= x_species(is)) EXIT
  END DO

  ! If the molecule can't move then return
  IF ( max_disp(is,ibox) == 0. ) RETURN

  ! Choose a molecule at random for displacement
  im = INT ( rranf() * nmols(is,ibox) ) + 1
  ! Get the LOCATE of imth molecule of species is in ibox
  lm = locate(im,is,ibox)

  ! update the trial counters
  tot_trials(ibox) = tot_trials(ibox) + 1
  ntrials(is,ibox)%displacement = ntrials(is,ibox)%displacement + 1

  ! obtain the energy of the molecule before the move.  Note that due to
  ! this move, the interatomic energies such as vdw and electrostatics will
  ! change. Also the ewald_reciprocal energy will change but there will
  ! be no change in intramolecular energies.
  IF (l_pair_nrg) THEN
    CALL Store_Molecule_Pair_Interaction_Arrays(lm,is,ibox,E_vdw,E_qq)
  ELSE
    CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_vdw,E_qq,inter_overlap)
  END IF

  IF (inter_overlap)  THEN
     err_msg = ""
     err_msg(1) = "Energy overlap detected in existing configuration"
     err_msg(2) = "of molecule " // TRIM(Int_To_String(lm)) // " of species " // TRIM(Int_To_String(is))
     CALL Clean_Abort(err_msg, "Translate")
  END IF

  ! Store the old positions of the atoms
  CALL Save_Old_Cartesian_Coordinates(lm,is)

  ! Generate a random displacement vector. Note that the current formalism will
  ! work for cubic shaped boxes. However, it is easy to extend for nonorthorhombic
  ! boxes where displacements along the basis vectors. 
  dx = ( 2.0_DP * rranf() - 1.0_DP) * max_disp(is,ibox)
  dy = ( 2.0_DP * rranf() - 1.0_DP) * max_disp(is,ibox)
  dz = ( 2.0_DP * rranf() - 1.0_DP) * max_disp(is,ibox)

  ! Move atoms by the above vector dx,dy,dz and also update the COM
  atom_list(:,lm,is)%rxp = atom_list(:,lm,is)%rxp + dx
  atom_list(:,lm,is)%ryp = atom_list(:,lm,is)%ryp + dy
  atom_list(:,lm,is)%rzp = atom_list(:,lm,is)%rzp + dz

  molecule_list(lm,is)%xcom = molecule_list(lm,is)%xcom + dx
  molecule_list(lm,is)%ycom = molecule_list(lm,is)%ycom + dy
  molecule_list(lm,is)%zcom = molecule_list(lm,is)%zcom + dz


  !**************************************************************************
  CALL Fold_Molecule(lm,is,ibox)
  
  CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_vdw_move,E_qq_move,inter_overlap)
  ! If an overlap is detected, immediately reject the move
  IF (inter_overlap) THEN ! Move is rejected
    
     CALL Revert_Old_Cartesian_Coordinates(lm,is)
     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(lm,is,ibox)

  ELSE

     dE = 0.0_DP

     IF ((int_charge_sum_style(ibox) == charge_ewald) .AND. (has_charge(is))) THEN

        ALLOCATE(cos_mol_old(nvecs(ibox)),sin_mol_old(nvecs(ibox)))
        CALL Get_Position_Alive(lm,is,position)
        
        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        cos_mol_old(:) = cos_mol(1:nvecs(ibox),position)
        sin_mol_old(:) = sin_mol(1:nvecs(ibox),position)
        !$OMP END PARALLEL WORKSHARE

        CALL Update_System_Ewald_Reciprocal_Energy(lm,is,ibox,int_translation,E_reciprocal_move)
        dE = E_reciprocal_move - energy(ibox)%reciprocal
        
     END IF
     
     ! Compute the difference in old and new energy
     dE = dE + ( E_vdw_move - E_vdw ) + ( E_qq_move - E_qq )

     IF (int_sim_type == sim_nvt_min) THEN
        IF (dE  <= 0.0_DP) THEN
           accept = .TRUE.
        END IF
     ELSE

         ln_pacc = beta(ibox) * dE
         accept = accept_or_reject(ln_pacc)

     END IF

     IF ( accept ) THEN

        ! accept the move and update the global energies
        energy(ibox)%inter_vdw = energy(ibox)%inter_vdw + E_vdw_move - E_vdw
        energy(ibox)%inter_q   = energy(ibox)%inter_q   + E_qq_move  - E_qq
        
        IF(int_charge_sum_style(ibox) == charge_ewald .AND. has_charge(is)) THEN
           energy(ibox)%reciprocal = E_reciprocal_move
        END IF
        energy(ibox)%inter = energy(ibox)%inter + dE
        energy(ibox)%total = energy(ibox)%total + dE

        ! update success counter
        nsuccess(is,ibox)%displacement = nsuccess(is,ibox)%displacement + 1
        nequil_success(is,ibox)%displacement = nequil_success(is,ibox)%displacement + 1

        IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
        IF (ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
        IF (ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)

     ELSE

        ! Revert to the old coordinates of atoms and com of the molecule
        CALL Revert_Old_Cartesian_Coordinates(lm,is) 

        IF ((int_charge_sum_style(ibox) == charge_ewald) .AND. (has_charge(is))) THEN
           ! Also reset the old cos_sum and sin_sum for reciprocal space vectors. Note
           ! that old vectors were set while difference in ewald reciprocal energy was computed.
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)           
           cos_sum(:,ibox) = cos_sum_old(:,ibox)
           sin_sum(:,ibox) = sin_sum_old(:,ibox)
           cos_mol(1:nvecs(ibox),position) =cos_mol_old(:)
           sin_mol(1:nvecs(ibox),position) =sin_mol_old(:)
           !$OMP END PARALLEL WORKSHARE
           DEALLOCATE(cos_mol_old,sin_mol_old)
        END IF
        IF (l_pair_nrg)  CALL Reset_Molecule_Pair_Interaction_Arrays(lm,is,ibox)
     ENDIF
     
  END IF

  IF (verbose_log) THEN
    WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8)') i_mcstep, 'translate' , lm, is, ibox, accept
  END IF

  IF ( MOD(ntrials(is,ibox)%displacement,nupdate) == 0 ) THEN
     IF ( int_run_type == run_equil ) THEN 
        success_ratio = REAL(nequil_success(is,ibox)%displacement,DP)/REAL(nupdate,DP)
     ELSE
        success_ratio = REAL(nsuccess(is,ibox)%displacement,DP)/REAL(ntrials(is,ibox)%displacement,DP)
     END IF

     WRITE(logunit,'(X,I9,X,A10,X,5X,X,I3,X,I3,X,F8.5)',ADVANCE='NO') &
           i_mcstep, 'translate', is, ibox, success_ratio

     !nsuccess(is,ibox)%displacement = 0

     IF ( int_run_type == run_equil ) THEN

        ! check if the acceptace is close to 0.5

         nequil_success(is,ibox)%displacement = 0

         IF  ( success_ratio < 0.00005 ) THEN
             max_disp(is,ibox) = 0.1_DP*max_disp(is,ibox)
         ELSE
             ! minimum max_disp for this species
             IF (has_charge(is) .AND. int_charge_style(ibox) /= charge_none) THEN
                rcut_small = MIN(rcut_vdw(ibox),rcut_coul(ibox))
             ELSE
                rcut_small = rcut_vdw(ibox)
             END IF
             max_disp(is,ibox) = MIN(rcut_small,2.0_DP*success_ratio*max_disp(is,ibox))
         END IF

         WRITE(logunit,'(X,F9.5)',ADVANCE='NO') max_disp(is,ibox)
        
     END IF

     WRITE(logunit,*)

  END IF

END SUBROUTINE Translate
