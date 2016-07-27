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

SUBROUTINE GEMC_NVT_Volume

  !***********************************************************************
  !
  ! The routine performs an NVT volume move for a GEMC NVT simulation.
  ! At present the simulation is restricted to two boxes. Future possibility
  ! may include inclusion of multiple boxes. 
  !
  ! Called by
  !
  !   gemc_driver
  !
  ! CALLS 
  !
  ! Revision history
  !
  !   12/05/13 : Beta Release
  !*************************************************************************

  USE Global_Variables
  USE Random_Generators
  USE Volume
  USE IO_Utilities
  USE Energy_Routines
  USE Pair_Nrg_Routines

  IMPLICIT NONE

  ! Local variables

  INTEGER :: box_grw, box_shk
  INTEGER :: ibox, nvecs_old_1, nvecs_old_2, nvecs_max
  INTEGER :: is, im, lm

  REAL(DP) :: delta_volume, ln_pacc, delta_e_1, delta_e_2
  REAL(DP) :: success_ratio
  REAL(DP), DIMENSION(maxk) :: hx_old_1, hy_old_1, hz_old_1, Cn_old_1
  REAL(DP), DIMENSION(maxk) :: hx_old_2, hy_old_2, hz_old_2, Cn_old_2  
  REAL(DP) :: v_ratio_o, v_total, vol_factor

  LOGICAL :: overlap, accept_or_reject

  TYPE(Box_Class) :: box_list_old_1, box_list_old_2
  TYPE(Energy_Class) :: energy_old_1, energy_old_2

  INTEGER :: position

  REAL(DP), ALLOCATABLE :: pair_nrg_vdw_old(:,:), pair_nrg_qq_old(:,:)
  REAL(DP), ALLOCATABLE :: cos_mol_old(:,:), sin_mol_old(:,:)

  REAL(DP) :: rcut_vdw_old_1, rcut_coul_old_1, rcut3_old_1, rcut9_old_1, alpha_ewald_old_1
  REAL(DP) :: h_ewald_cut_old_1, rcut_vdwsq_old_1, rcut_coulsq_old_1, rcut_vdw3_old_1
  REAL(DP) :: rcut_vdw6_old_1, rcut_max_old_1

  REAL(DP) :: rcut_vdw_old_2, rcut_coul_old_2, rcut3_old_2, rcut9_old_2, alpha_ewald_old_2
  REAL(DP) :: h_ewald_cut_old_2, rcut_vdwsq_old_2, rcut_coulsq_old_2, rcut_vdw3_old_2
  REAL(DP) :: rcut_vdw6_old_2, rcut_max_old_2

  accept = .FALSE.

  ! Pick the box that will grow with uniform probability
  box_grw = INT(rranf() * nbr_boxes) + 1
  
  ! Pick the box that will shrink with uniform probability
  box_shk = INT(rranf() * (nbr_boxes - 1)) + 1
  IF (box_grw <= box_shk) box_shk = box_shk + 1

  tot_trials(box_grw) = tot_trials(box_grw) + 1
  tot_trials(box_shk) = tot_trials(box_shk) + 1

  nvolumes(box_grw) = nvolumes(box_grw) + 1
  nvolumes(box_shk) = nvolumes(box_shk) + 1

  ! store old cell matrix 
  box_list_old_1 = box_list(box_grw)
  box_list_old_2 = box_list(box_shk)

  ! Store the old configurations of all atoms and COMs
  CALL Save_Cartesian_Coordinates_Box(box_grw)
  CALL Save_Cartesian_Coordinates_Box(box_shk)

  ! store the pair interactions
  IF (l_pair_nrg) THEN
     ALLOCATE(pair_nrg_vdw_old(SUM(max_molecules),SUM(max_molecules)))
     ALLOCATE(pair_nrg_qq_old(SUM(max_molecules),SUM(max_molecules)))

     !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
     pair_nrg_vdw_old(:,:) = pair_nrg_vdw(:,:)
     pair_nrg_qq_old(:,:) = pair_nrg_qq(:,:)
     !!$OMP END PARALLEL WORKSHARE
  END IF

  ! store cos_mol and sin_mol
  IF ( int_charge_sum_style(box_grw) == charge_ewald .OR. &
       int_charge_sum_style(box_shk) == charge_ewald ) THEN
     
     ALLOCATE(cos_mol_old(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)
     
     IF (AllocateStatus /= 0 ) THEN
        err_msg = ''
        err_msg(1) = 'Memory could not be allocated for cos_mol_old'
        CALL Clean_Abort(err_msg,'gemc_nvt_volume.f90')
     END IF
     
     ALLOCATE(sin_mol_old(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)

     IF (AllocateStatus /= 0 ) THEN
        err_msg = ''
        err_msg(1) = 'Memory could not be allocated for sin_mol_old'
        CALL Clean_Abort(err_msg,'gemc_nvt_volume.f90')
     END IF
     
     !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
     cos_mol_old(:,:) = cos_mol(:,:)
     sin_mol_old(:,:) = sin_mol(:,:)
     !!$OMP END PARALLEL WORKSHARE

  END IF

  IF (f_dv) THEN
     
     delta_volume = rranf() * box_list(box_grw)%dv_max
     
     ! update each of the box volumes
     box_list(box_grw)%volume = box_list(box_grw)%volume + delta_volume
     box_list(box_shk)%volume = box_list(box_shk)%volume - delta_volume
     
  ELSE IF (f_vratio) THEN

     v_ratio_o = (box_list(box_grw)%volume/box_list(box_shk)%volume)
     v_total = box_list(box_grw)%volume + box_list(box_shk)%volume

     delta_volume = rranf() * box_list(box_grw)%dv_max
     vol_factor = v_ratio_o * EXP(delta_volume)

     ! update each of the box volumes
     box_list(box_grw)%volume = v_total * vol_factor / (1.0_DP + vol_factor)
     box_list(box_shk)%volume = v_total - box_list(box_grw)%volume
  
  END IF

  ! Assumes box_grw ix cubic
  box_list(box_grw)%length(1,1) = (box_list(box_grw)%volume) ** (1.0_DP/3.0_DP)
  box_list(box_grw)%length(2,2) = box_list(box_grw)%length(1,1)
  box_list(box_grw)%length(3,3) = box_list(box_grw)%length(2,2)
  
  ! Assumes box_shk ix cubic
  box_list(box_shk)%length(1,1) = (box_list(box_shk)%volume) ** (1.0_DP/3.0_DP)
  box_list(box_shk)%length(2,2) = box_list(box_shk)%length(1,1)
  box_list(box_shk)%length(3,3) = box_list(box_shk)%length(2,2)

  CALL Compute_Cell_Dimensions(box_grw)
  CALL Compute_Cell_Dimensions(box_shk)

  IF ( l_half_len_cutoff(box_grw)) THEN

     ! store old cutoffs and other associated quantities
     rcut_vdw_old_1 = rcut_vdw(box_grw)
     rcut_coul_old_1 = rcut_coul(box_grw)
     rcut_vdwsq_old_1 = rcut_vdwsq(box_grw)
     rcut_coulsq_old_1 = rcut_coulsq(box_grw)

     rcut3_old_1 = rcut3(box_grw)
     rcut9_old_1 = rcut9(box_grw)
     rcut_vdw3_old_1 = rcut_vdw3(box_grw)
     rcut_vdw6_old_1 = rcut_vdw6(box_grw)

     rcut_max_old_1 = rcut_max(box_grw)
     
     IF (int_charge_sum_style(box_grw) == charge_ewald) THEN
!        alpha_ewald_old_1 = alpha_ewald(box_grw)
        h_ewald_cut_old_1 = h_ewald_cut(box_grw)
     END IF

     ! change cutoffs and other associated quantities
     rcut_vdw(box_grw) = 0.5_DP * box_list(box_grw)%length(1,1)
     rcut_coul(box_grw) = rcut_vdw(box_grw)
     rcut_vdwsq(box_grw) = rcut_vdw(box_grw) * rcut_vdw(box_grw)
     rcut_coulsq(box_grw) = rcut_vdwsq(box_grw)

     rcut_vdw3(box_grw) = rcut_vdwsq(box_grw) * rcut_vdw(box_grw)
     rcut_vdw6(box_grw) = rcut_vdw3(box_grw) * rcut_vdw3(box_grw)
     rcut3(box_grw) = rcut_vdw3(box_grw)
     rcut9(box_grw) = rcut3(box_grw) * rcut_vdw6(box_grw)

     rcut_max(box_grw) = rcut_vdw(box_grw)

     IF (int_charge_sum_style(box_grw) == charge_ewald) THEN
!        alpha_ewald(box_grw) = ewald_p_sqrt(box_grw) / rcut_coul(box_grw)
        h_ewald_cut(box_grw) = 2.0_DP * ewald_p(box_grw) / rcut_coul(box_grw)
     END IF

  ELSE

     IF ( 0.5_DP * box_list(box_grw)%length(1,1) < rcut_vdw(box_grw) .OR. &
          0.5_DP * box_list(box_grw)%length(1,1) < rcut_coul(box_grw) .OR. &
          0.5_DP * box_list(box_grw)%length(1,1) < roff_charmm(box_grw) .OR. &
          0.5_DP * box_list(box_grw)%length(1,1) < roff_switch(box_grw) ) THEN
        err_msg = ''
        err_msg(1) = 'Cutoff is greater than the half box length'
        err_msg(2) = Int_to_String(box_grw)
        CALL Clean_Abort(err_msg,'GEMC_NVT_Volume')
     END IF
     
  END IF

  IF ( l_half_len_cutoff(box_shk)) THEN

     ! store old cutoffs and other associated quantities
     rcut_vdw_old_2 = rcut_vdw(box_shk)
     rcut_coul_old_2 = rcut_coul(box_shk)
     rcut_vdwsq_old_2 = rcut_vdwsq(box_shk)
     rcut_coulsq_old_2 = rcut_coulsq(box_shk)

     rcut3_old_2 = rcut3(box_shk)
     rcut9_old_2 = rcut9(box_shk)
     rcut_vdw3_old_2 = rcut_vdw3(box_shk)
     rcut_vdw6_old_2 = rcut_vdw6(box_shk)

     rcut_max_old_2 = rcut_max(box_shk)
     
     IF (int_charge_sum_style(box_shk) == charge_ewald) THEN
!        alpha_ewald_old_2 = alpha_ewald(box_shk)
        h_ewald_cut_old_2 = h_ewald_cut(box_shk)
     END IF

     ! store cutoffs and other associated quantities
     rcut_vdw(box_shk) = 0.5_DP * box_list(box_shk)%length(1,1)
     rcut_coul(box_shk) = rcut_vdw(box_shk)
     rcut_vdwsq(box_shk) = rcut_vdw(box_shk) * rcut_vdw(box_shk)
     rcut_coulsq(box_shk) = rcut_vdwsq(box_shk)

     rcut_vdw3(box_shk) = rcut_vdwsq(box_shk) * rcut_vdw(box_shk)
     rcut_vdw6(box_shk) = rcut_vdw3(box_shk) * rcut_vdw3(box_shk)
     rcut3(box_shk) = rcut_vdw3(box_shk)
     rcut9(box_shk) = rcut3(box_shk) * rcut_vdw6(box_shk)

     rcut_max(box_shk) = rcut_vdw(box_shk)
     
     IF (int_charge_sum_style(box_shk) == charge_ewald) THEN
!        alpha_ewald(box_shk) = ewald_p_sqrt(box_shk) / rcut_coul(box_shk)
        h_ewald_cut(box_shk) = 2.0_DP * ewald_p(box_shk) / rcut_coul(box_shk)
     END IF

  ELSE
     
     IF ( 0.5_DP * box_list(box_shk)%length(1,1) < rcut_vdw(box_shk) .OR. &
          0.5_DP * box_list(box_shk)%length(1,1) < rcut_coul(box_shk) .OR. &
          0.5_DP * box_list(box_shk)%length(1,1) < roff_charmm(box_shk) .OR. &
          0.5_DP * box_list(box_shk)%length(1,1) < roff_switch(box_shk) ) THEN
        err_msg = ''
        err_msg(1) = 'Cutoff is greater than the half box length'
        err_msg(2) = Int_To_String(box_shk)
        CALL Clean_Abort(err_msg,'GEMC_NVT_Volume')
     END IF
     
  END IF


  ! Rescale the COM and all atomic positions

  CALL Scale_COM_Cartesian(box_grw,box_list_old_1)
  CALL Scale_COM_Cartesian(box_shk,box_list_old_2)

  ! Now let us compute the energy change due to the combined move

  energy_old_1 = energy(box_grw)
  energy_old_2 = energy(box_shk)

  IF ( int_charge_sum_style(box_grw) == charge_ewald) THEN
     ! Then we need to determine if the number of k vectors change due to this volume move
     ! we basically follow the procedure outlined in volume_change.f90 for the two boxes.
     ! Note that the number of k vectors will increase only for the box that is expanding

     ! store old terms 

     nvecs_old_1 = nvecs(box_grw)
     nvecs_old_2 = nvecs(box_shk)
     nvecs_max = MAXVAL(nvecs)

     !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED) 
     cos_sum_old(:,:) = cos_sum(:,:)
     sin_sum_old(:,:) = sin_sum(:,:)

     hx_old_1(:) = hx(:,box_grw)
     hy_old_1(:) = hy(:,box_grw)
     hz_old_1(:) = hz(:,box_grw)
     Cn_old_1(:) = Cn(:,box_grw)

     hx_old_2(:) = hx(:,box_shk)
     hy_old_2(:) = hy(:,box_shk)
     hz_old_2(:) = hz(:,box_shk)
     Cn_old_2(:) = Cn(:,box_shk)

     !!$OMP END PARALLEL WORKSHARE

     ! Determine the new k vectors for this box. The call will change Cn, hx, hy and hz and hence will
     ! change cos_sum and sin_sum.
     
     CALL Ewald_Reciprocal_Lattice_Vector_Setup(box_grw)
     CALL Ewald_Reciprocal_Lattice_Vector_Setup(box_shk)

     ! reallocate arrays

     DEALLOCATE(cos_sum,sin_sum)

     IF (ALLOCATED(cos_mol)) DEALLOCATE(cos_mol)
     IF (ALLOCATED(sin_mol)) DEALLOCATE(sin_mol)

     ALLOCATE(cos_sum(MAXVAL(nvecs),nbr_boxes),Stat=AllocateStatus)
     
     IF (AllocateStatus /= 0) THEN
        err_msg = ''
        err_msg(1) = 'Memory could not be allocated for cos_sum'
        CALL Clean_Abort(err_msg,'GEMC_NVT_VOLUME')
     END IF
     
     
     ALLOCATE(sin_sum(MAXVAL(nvecs),nbr_boxes), Stat = AllocateStatus)
     
     IF (Allocatestatus /= 0) THEN
        err_msg = ''
        err_msg(1) = 'Memory could not be allocated for sin_sum'
        CALL Clean_Abort(err_msg,'GEMC_NVT_VOLUME')
     END IF

     ALLOCATE(cos_mol(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)
     
     IF (AllocateStatus /= 0) THEN
        err_msg = ''
        err_msg(1) = 'Memory could not be allocated for cos_mol'
        CALL Clean_Abort(err_msg,'gemc_volume_change.f90')
     END IF
     
     ALLOCATE(sin_mol(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)
     
     IF (AllocateStatus /= 0) THEN
        err_msg = ''
        err_msg(1) = 'Memory could not be allocated for sin_mol'
        CALL Clean_Abort(err_msg,'gemc_volume_change.f90')
     END IF
        
  END IF
  
  ! Compute the energies of the boxes. For computational efficiency,
  ! we will determine the energy of the box whose volume decreased
  ! as this box will likely to have more overlaps

  CALL Compute_System_Total_Energy(box_shk, .TRUE., overlap)
  
   IF (overlap) THEN 
      CALL Reset_Coords

      IF (verbose_log) THEN
         WRITE(logunit,'(X,I9,X,A10,X,5X,X,3X,X,I1,A1,I1,X,L8,X,9X,X,A9)') &
               i_mcstep, 'vol_swap', box_shk, '>', box_grw, accept, 'overlap'
      END IF

   ELSE
      
      CALL Compute_System_Total_Energy(box_grw, .TRUE.,overlap)
     
      ! actually there should be no overlap for the box whose dimensions
      ! increase but we will include this check only for safety.
      
      IF (overlap) THEN
         CALL Reset_Coords

         IF (verbose_log) THEN
            WRITE(logunit,'(X,I9,X,A10,X,5X,X,3X,X,I1,A1,I1,X,L8,X,9X,X,A9)') &
                  i_mcstep, 'vol_swap', box_shk, '>', box_grw, accept, 'overlap'
         END IF

      ELSE
         ! accept or reject the move based on the acceptance rule
         
         delta_e_1 = energy(box_grw)%total - energy_old_1%total
         delta_e_2 = energy(box_shk)%total - energy_old_2%total
         
         IF (f_dv) THEN
            ln_pacc = beta(box_grw) * delta_e_1 + beta(box_shk) * delta_e_2 &
                    - REAL(SUM(nmols(:,box_grw)),DP) * DLOG( box_list(box_grw)%volume / box_list_old_1%volume) &
                    - REAL(SUM(nmols(:,box_shk)),DP) * DLOG( box_list(box_shk)%volume / box_list_old_2%volume)
            
         ELSE IF(f_vratio) THEN
            
            ln_pacc = beta(box_grw) * delta_e_1 + beta(box_shk) * delta_e_2 &
                    - REAL(SUM(nmols(:,box_grw))+1,DP) * DLOG( box_list(box_grw)%volume / box_list_old_1%volume) &
                    - REAL(SUM(nmols(:,box_shk))+1,DP) * DLOG( box_list(box_shk)%volume / box_list_old_2%volume)
            
         END IF

         accept = accept_or_reject(ln_pacc)
         
         IF (accept) THEN
            
            nvol_success(box_grw) = nvol_success(box_grw) + 1
            nvol_success(box_shk) = nvol_success(box_shk) + 1
            ivol_success(box_grw) = ivol_success(box_grw) + 1
            ivol_success(box_shk) = ivol_success(box_shk) + 1
            ! energy,positions and box dimensions are already updated
            IF (int_charge_sum_style(box_grw) == charge_ewald) THEN
               
               ! cos_sum and sin_sum were deallocated, destroying the terms
               ! for boxes other than box_grw, box_shk. so restore these
               
               DO ibox = 1, nbr_boxes
                  
                  IF ( .NOT. ( (ibox /= box_grw) .OR. (ibox /= box_shk))) THEN
                     ! transfer cos_sum and sin_sum for other boxes
                     ! Note that direct assignment of cos_sum_old to cos_sum
                     ! will result into an error as these two arrays are of
                     ! different dimensions
                     cos_sum(1:nvecs(ibox),ibox) = cos_sum_old(1:nvecs(ibox),ibox)
                     sin_sum(1:nvecs(ibox),ibox) = sin_sum_old(1:nvecs(ibox),ibox)
                     
                  END IF
                  
               END DO
               ! Now deallocate cos_sum_old and sin_sum_old so that they have the same dimensions
               ! as sin_sum and cos_sum                  
               DEALLOCATE(cos_sum_old,sin_sum_old)
               ALLOCATE(cos_sum_old(SIZE(cos_sum,1),nbr_boxes))
               ALLOCATE(sin_sum_old(SIZE(sin_sum,1),nbr_boxes))
               
               ! Now assign cos_mol and sin_mol for the molecules present in other
               ! boxes. Note that cos_mol for box_grw and box_shk have been assigned
               ! when a call to Compute_System_Total_Energy was placed.
               
               DO ibox = 1, nbr_boxes
                  IF (ibox == box_grw .OR. ibox == box_shk) CYCLE
                        
                  DO is = 1, nspecies
                     DO im = 1, nmols(is,ibox)
                        lm = locate(im,is,ibox)
                        
                        CALL Get_Position_Alive(lm,is,position)
                        
                        !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
                        cos_mol(1:nvecs(ibox),position) = cos_mol_old(1:nvecs(ibox),position)
                        sin_mol(1:nvecs(ibox),position) = sin_mol_old(1:nvecs(ibox),position)
                        
                        cos_mol(nvecs(ibox)+1:MAXVAL(nvecs),position) = 0.0_DP
                        sin_mol(nvecs(ibox)+1:MAXVAL(nvecs),position) = 0.0_DP
                        !!$OMP END PARALLEL WORKSHARE
                        
                     END DO
                  END DO
               END DO
               DEALLOCATE(cos_mol_old,sin_mol_old)
            END IF  ! ends if of 316
            
            IF (l_pair_nrg) DEALLOCATE(pair_nrg_vdw_old,pair_nrg_qq_old)
  
         ELSE
            
            CALL Reset_Coords

         END IF
        
         IF (verbose_log) THEN
            WRITE(logunit,'(X,I9,X,A10,X,5X,X,3X,X,I1,A1,I1,X,L8,X,9X,X,F9.3)') &
                  i_mcstep, 'vol_swap', box_shk, '>', box_grw, accept, ln_pacc
         END IF

      END IF

   END IF

  ! Update the maximum volume modulus of equilibration runs
   IF (MOD(nvolumes(box_grw),nvol_update) == 0 ) THEN
      IF ( int_run_type == run_equil) THEN

         success_ratio = REAL(ivol_success(box_grw),DP)/REAL(nvol_update,DP)
      
         ivol_success(box_grw) = 0
         ivol_success(box_shk) = 0
      
         IF ( success_ratio < 0.0001 ) THEN
            ! decrease the maximum displacement to 10%
            
            box_list(:)%dv_max = 0.1_DP * box_list(:)%dv_max
            
         ELSE
            
            box_list(:)%dv_max = 2.0_DP*success_ratio * box_list(:)%dv_max

         END IF

         WRITE(logunit,'(X,I9,X,A10,X,5X,X,3X,X,3X,X,F8.5,X,F9.0)') &
               i_mcstep, 'vol_swap', success_ratio, box_list(box_grw)%dv_max

      ELSE

         success_ratio = REAL(nvol_success(box_grw),DP)/REAL(nvolumes(box_grw),DP)
         WRITE(logunit,'(X,I9,X,A10,X,5X,X,3X,X,3X,X,F8.5)') &
               i_mcstep, 'vol_swap', success_ratio
         
      END IF
         
   END IF

   
 CONTAINS
   
   SUBROUTINE Reset_Coords
     
     IMPLICIT NONE
     
     CALL Reset_Cartesian_Coordinates_Box(box_grw)
     CALL Reset_Cartesian_Coordinates_Box(box_shk)
     
     ! box list and energy 
     
     box_list(box_grw) = box_list_old_1
     box_list(box_shk) = box_list_old_2
     
     energy(box_grw) = energy_old_1
     energy(box_shk) = energy_old_2

     IF (l_pair_nrg) THEN
        
        !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        pair_nrg_vdw(:,:) = pair_nrg_vdw_old(:,:)
        pair_nrg_qq(:,:) = pair_nrg_qq_old(:,:)
        !!$OMP END PARALLEL WORKSHARE
        
        DEALLOCATE(pair_nrg_vdw_old,pair_nrg_qq_old)
     END IF

     IF (l_half_len_cutoff(box_grw)) THEN
        
        rcut_vdw(box_grw) = rcut_vdw_old_1
        rcut_coul(box_grw) = rcut_coul_old_1
        rcut_vdwsq(box_grw) = rcut_vdwsq_old_1
        rcut_coulsq(box_grw) = rcut_coulsq_old_1
        
        rcut3(box_grw) = rcut3_old_1
        rcut9(box_grw) = rcut9_old_1
        rcut_vdw3(box_grw) = rcut_vdw3_old_1
        rcut_vdw6(box_grw) = rcut_vdw6_old_1
        
        rcut_max(box_grw) = rcut_max_old_1
        
        IF( int_charge_sum_style(box_grw) == charge_ewald ) THEN
           
           !           alpha_ewald(box_grw) = alpha_ewald_old_1
           h_ewald_cut(box_grw) = h_ewald_cut_old_1
           
        END IF
        
        
     END IF

     IF (l_half_len_cutoff(box_shk)) THEN
        
        rcut_vdw(box_shk) = rcut_vdw_old_2
        rcut_coul(box_shk) = rcut_coul_old_2
        rcut_vdwsq(box_shk) = rcut_vdwsq_old_2
        rcut_coulsq(box_shk) = rcut_coulsq_old_2
        
        rcut3(box_shk) = rcut3_old_2
        rcut9(box_shk) = rcut9_old_2
        rcut_vdw3(box_shk) = rcut_vdw3_old_2
        rcut_vdw6(box_shk) = rcut_vdw6_old_2
        
        rcut_max(box_shk) = rcut_max_old_2

        IF( int_charge_sum_style(box_shk) == charge_ewald ) THEN
           
!           alpha_ewald(box_shk) = alpha_ewald_old_2
           h_ewald_cut(box_shk) = h_ewald_cut_old_2
           
        END IF
        
        
     END IF
     
     IF (int_charge_sum_style(box_grw) == charge_ewald) THEN
        
       nvecs(box_grw) = nvecs_old_1
       nvecs(box_shk) = nvecs_old_2
       
       DEALLOCATE(cos_sum,sin_sum)
       DEALLOCATE(cos_mol,sin_mol)
       
       ALLOCATE(cos_sum(MAXVAL(nvecs),nbr_boxes),stat = AllocateStatus)
       
       IF (Allocatestatus /= 0) THEN
          err_msg = ''
          err_msg(1) = 'Memory could not be allocated for cos_sum'
          err_msg(2) = 'volume move rejected'
          CALL Clean_Abort(err_msg,'Volume_Change')
       END IF
       
       ALLOCATE(sin_sum(MAXVAL(nvecs),nbr_boxes),Stat = Allocatestatus)
       IF (Allocatestatus /= 0) THEN
          err_msg = ''
          err_msg(1) = 'Memory could not be allocated in the volume rejection'
          CALL Clean_Abort(err_msg,'Volume_Change')
       END IF
       
       ALLOCATE(cos_mol(MAXVAL(nvecs),SUM(max_molecules)),Stat = Allocatestatus)
       IF (Allocatestatus /= 0) THEN
          err_msg = ''
          err_msg(1) = 'Memory could not be allocated for cos_mol in the volume rejection'
          CALL Clean_Abort(err_msg,'GEMC NVT Volume_Change')
       END IF
       
       ALLOCATE(sin_mol(MAXVAL(nvecs),SUM(max_molecules)),Stat = Allocatestatus)
       IF (Allocatestatus /= 0) THEN
          err_msg = ''
          err_msg(1) = 'Memory could not be allocated for sin_mol in the volume rejection'
          CALL Clean_Abort(err_msg,'GEMC NVT Volume_Change')
       END IF
          
       !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
       
       cos_sum(:,:) = cos_sum_old(:,:)
       sin_sum(:,:) = sin_sum_old(:,:)
       
       cos_mol(1:SIZE(cos_mol_old,1),:) = cos_mol_old(:,:)
       sin_mol(1:SIZE(sin_mol_old,1),:) = sin_mol_old(:,:)
       
       hx(:,box_grw) = hx_old_1(:)
       hy(:,box_grw) = hy_old_1(:)
       hz(:,box_grw) = hz_old_1(:)
       Cn(:,box_grw) = Cn_old_1(:)
       
       hx(:,box_shk) = hx_old_2(:)
       hy(:,box_shk) = hy_old_2(:)
       hz(:,box_shk) = hz_old_2(:)
       Cn(:,box_shk) = Cn_old_2(:)

       !!$OMP END PARALLEL WORKSHARE
       
       DEALLOCATE(cos_mol_old,sin_mol_old)
       
       ! here we make sure that cos_sum_old and sin_sum_old have the same dimensions
       ! as cos_sum and sin_sum
       
       DEALLOCATE(cos_sum_old,sin_sum_old)
       ALLOCATE(cos_sum_old(SIZE(cos_sum,1),nbr_boxes),sin_sum_old(SIZE(sin_sum,1),nbr_boxes))
       
    END IF
    
  END SUBROUTINE Reset_Coords

END SUBROUTINE GEMC_NVT_Volume


