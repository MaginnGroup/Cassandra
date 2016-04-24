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
SUBROUTINE Translate(this_box)
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

  IMPLICIT NONE

  ! Arguments
  INTEGER  :: this_box   ! box index
  INTEGER  :: mc_step

  ! Local declarations
  INTEGER  :: ibox       ! box index
  INTEGER  :: is         ! species index
  INTEGER  :: im, alive  ! molecule indices
  INTEGER  :: total_mols ! number of molecules in the system
  
  REAL(DP) :: nmols_box(nbr_boxes)
  REAL(DP), ALLOCATABLE :: x_box(:), x_species(:)
  REAL(DP) :: rand_no
  REAL(DP) :: dx, dy, dz
  REAL(DP) :: delta_e, ln_pacc, success_ratio
  REAL(DP) :: E_vdw, E_qq, E_vdw_move, E_qq_move, E_reciprocal_move
  REAL(DP) :: rcut_small

  LOGICAL :: inter_overlap, overlap, accept, accept_or_reject

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

    ALLOCATE(x_box(nbr_boxes)) 

    DO ibox = 1, nbr_boxes
       x_box(ibox) = REAL(nmols_box(ibox),DP)/REAL(total_mols,DP)
       IF ( ibox > 1 ) THEN
          x_box(ibox) = x_box(ibox) + x_box(ibox-1)
       END IF
    END DO
  
    DO ibox = 1, nbr_boxes
       IF ( rranf() <= x_box(ibox)) EXIT
    END DO

    this_box = ibox
    DEALLOCATE(x_box)

  ELSE

    this_box = 1

  END IF

  ! If there are no molecules in this box then return
  IF( nmols_box(this_box) == 0 ) RETURN

  ! Choose species based on the mol fraction, using Golden sampling
  ALLOCATE(x_species(nspecies))

  DO is = 1, nspecies
     IF ( max_disp(is,this_box) > 0. ) THEN
       x_species(is) = REAL(nmols(is,this_box), DP)/REAL(nmols_box(this_box),DP)
     ELSE
       x_species(is) = 0.0_DP
     END IF
     IF ( is > 1 ) THEN
        x_species(is) = x_species(is) + x_species(is-1)
     END IF
  END DO

  rand_no = rranf()
  DO is = 1, nspecies
     IF( rand_no <= x_species(is)) EXIT
  END DO

  DEALLOCATE(x_species)

  ! If the molecule can't move then return
  IF ( max_disp(is,this_box) == 0. ) RETURN

  ! Choose a molecule at random for displacement
  im = INT ( rranf() * nmols(is,this_box) ) + 1
  ! Get the index of imth molecule of species is in this_box
  CALL Get_Index_Molecule(this_box,is,im,alive)

  ! update the trial counters
  tot_trials(this_box) = tot_trials(this_box) + 1
  ntrials(is,this_box)%displacement = ntrials(is,this_box)%displacement + 1

  ! obtain the energy of the molecule before the move.  Note that due to
  ! this move, the interatomic energies such as vdw and electrostatics will
  ! change. Also the ewald_reciprocal energy will change but there will
  ! be no change in intramolecular energies.
  IF (l_pair_nrg) THEN
    CALL Store_Molecule_Pair_Interaction_Arrays(alive,is,this_box,E_vdw,E_qq)
  ELSE
    CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_vdw,E_qq,inter_overlap)
  END IF

  IF (inter_overlap)  THEN
     WRITE(*,*) 'Disaster, overlap in the old configruation'
     WRITE(*,*) 'Translate'
     WRITE(*,*) alive, is, this_box
  END IF

  ! Store the old positions of the atoms
  CALL Save_Old_Cartesian_Coordinates(alive,is)

  ! Generate a random displacement vector. Note that the current formalism will
  ! work for cubic shaped boxes. However, it is easy to extend for nonorthorhombic
  ! boxes where displacements along the basis vectors. 
  dx = ( 2.0_DP * rranf() - 1.0_DP) * max_disp(is,this_box)
  dy = ( 2.0_DP * rranf() - 1.0_DP) * max_disp(is,this_box)
  dz = ( 2.0_DP * rranf() - 1.0_DP) * max_disp(is,this_box)

  ! Move atoms by the above vector dx,dy,dz and also update the COM
  atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp + dx
  atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp + dy
  atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp + dz

  molecule_list(alive,is)%xcom = molecule_list(alive,is)%xcom + dx
  molecule_list(alive,is)%ycom = molecule_list(alive,is)%ycom + dy
  molecule_list(alive,is)%zcom = molecule_list(alive,is)%zcom + dz


  !**************************************************************************
  CALL Fold_Molecule(alive,is,this_box)
  
  CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_vdw_move,E_qq_move,inter_overlap)
  ! If an overlap is detected, immediately reject the move
  IF (inter_overlap) THEN ! Move is rejected
    
     CALL Revert_Old_Cartesian_Coordinates(alive,is)
     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,this_box)

  ELSE

     delta_e = 0.0_DP

     IF ((int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is))) THEN

        ALLOCATE(cos_mol_old(nvecs(this_box)),sin_mol_old(nvecs(this_box)))
        CALL Get_Position_Alive(alive,is,position)
        
        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        cos_mol_old(:) = cos_mol(1:nvecs(this_box),position)
        sin_mol_old(:) = sin_mol(1:nvecs(this_box),position)
        !$OMP END PARALLEL WORKSHARE

        CALL Update_System_Ewald_Reciprocal_Energy(alive,is,this_box,int_translation,E_reciprocal_move)
        delta_e = E_reciprocal_move - energy(this_box)%ewald_reciprocal
        
     END IF
     
     ! Compute the difference in old and new energy
     delta_e = ( E_vdw_move - E_vdw ) + ( E_qq_move - E_qq ) + delta_e

     IF (int_sim_type == sim_nvt_min) THEN
        IF (delta_e  <= 0.0_DP) THEN
           accept = .TRUE.
        ELSE
           accept = .FALSE.
        END IF
     ELSE

         ln_pacc = beta(this_box) * delta_e
         accept = accept_or_reject(ln_pacc)

     END IF

     IF ( accept ) THEN

        ! accept the move and update the global energies
        energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + E_vdw_move - E_vdw
        energy(this_box)%inter_q   = energy(this_box)%inter_q   + E_qq_move - E_qq
        
        IF(int_charge_sum_style(this_box) == charge_ewald .AND. has_charge(is)) THEN
           energy(this_box)%ewald_reciprocal =  E_reciprocal_move
        END IF
        energy(this_box)%total = energy(this_box)%total + delta_e


        ! update success counter
        
        nsuccess(is,this_box)%displacement = nsuccess(is,this_box)%displacement + 1
        nequil_success(is,this_box)%displacement = nequil_success(is,this_box)%displacement + 1

        IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
        IF (ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
        IF (ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)

     ELSE

        ! Revert to the old coordinates of atoms and com of the molecule
        CALL Revert_Old_Cartesian_Coordinates(alive,is) 

        IF ((int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is))) THEN
           ! Also reset the old cos_sum and sin_sum for reciprocal space vectors. Note
           ! that old vectors were set while difference in ewald reciprocal energy was computed.
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)           
           cos_sum(:,this_box) = cos_sum_old(:,this_box)
           sin_sum(:,this_box) = sin_sum_old(:,this_box)
           cos_mol(1:nvecs(this_box),position) =cos_mol_old(:)
           sin_mol(1:nvecs(this_box),position) =sin_mol_old(:)
           !$OMP END PARALLEL WORKSHARE
           DEALLOCATE(cos_mol_old,sin_mol_old)
        END IF
        IF (l_pair_nrg)  CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,this_box)
     ENDIF
     
  END IF

  IF ( MOD(ntrials(is,this_box)%displacement,nupdate) == 0 ) THEN
     IF ( int_run_style == run_equil ) THEN 
        success_ratio = REAL(nequil_success(is,this_box)%displacement,DP)/REAL(nupdate,DP)
     ELSE
        success_ratio = REAL(nsuccess(is,this_box)%displacement,DP)/REAL(ntrials(is,this_box)%displacement,DP)
     END IF

     WRITE(logunit,*)
     WRITE(logunit,'(A,I3,A,I1,A,F8.5)')'Success ratio, translation of species ', is , ' in box ', this_box, ' : ', success_ratio

     !nsuccess(is,this_box)%displacement = 0

     IF ( int_run_style == run_equil ) THEN

        ! check if the acceptace is close to 0.5

         nequil_success(is,this_box)%displacement = 0

         IF  ( success_ratio < 0.00005 ) THEN
             max_disp(is,this_box) = 0.1_DP*max_disp(is,this_box)
         ELSE
             ! minimum max_disp for this species
             IF (has_charge(is)) THEN
                rcut_small = MIN(rcut_vdw(this_box),rcut_coul(this_box))
             ELSE
                rcut_small = rcut_vdw(this_box)
             END IF
             max_disp(is,this_box) = MIN(rcut_small,2.0_DP*success_ratio*max_disp(is,this_box))
         END IF

         WRITE(logunit,'(A,I3,A,I1,A,F8.5)') 'Maximum width, translation of species ', is,' in box ', this_box, ' : ', max_disp(is,this_box)
        
     END IF

  END IF

END SUBROUTINE Translate
     

     
