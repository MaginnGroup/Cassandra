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
  INTEGER :: ibox, nvecs_grw, nvecs_shk, nvecs_max, nvecs_max_new_p4
  INTEGER :: is, im, lm

  REAL(DP) :: delta_volume, ln_pacc, dE_grw, dE_shk
  REAL(DP) :: success_ratio
  REAL(DP) :: v_ratio_o, v_total, vol_factor

  LOGICAL :: overlap, accept_or_reject

  REAL(DP), DIMENSION(3,3) :: length_grw, length_shk, length_inv_grw, length_inv_shk
  REAL(DP) :: volume_grw, volume_shk
  TYPE(Energy_Class) :: energy_grw, energy_shk

  INTEGER :: pos
  INTEGER :: istart, iend, im_locate_shift, nboxmols
  INTEGER, DIMENSION(MAXVAL(SUM(nmols,1))) :: posvec

  REAL(DP), ALLOCATABLE :: pair_nrg_vdw_old(:,:), pair_nrg_qq_old(:,:)
  REAL(DP), ALLOCATABLE :: cos_mol_old(:,:), sin_mol_old(:,:)
  REAL(DP), ALLOCATABLE :: kspace_vectors_grw(:,:), kspace_vectors_shk(:,:)
  INTEGER, ALLOCATABLE :: kspace_vector_ints_grw(:), kspace_vector_ints_shk(:)

  REAL(DP) :: rcut_vdw_grw, rcut_coul_grw, rcut3_grw, rcut9_grw, alpha_ewald_grw
  REAL(DP) :: h_ewald_cut_grw, rcut_vdwsq_grw, rcut_coulsq_grw, rcut_vdw3_grw
  REAL(DP) :: rcut_vdw6_grw, rcut_max_grw

  REAL(DP) :: rcut_vdw_shk, rcut_coul_shk, rcut3_shk, rcut9_shk, alpha_ewald_shk
  REAL(DP) :: h_ewald_cut_shk, rcut_vdwsq_shk, rcut_coulsq_shk, rcut_vdw3_shk
  REAL(DP) :: rcut_vdw6_shk, rcut_max_shk

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

  !box_list_grw = box_list(box_grw)
  !box_list_shk = box_list(box_shk)
  length_grw = box_list(box_grw)%length
  length_inv_grw = box_list(box_grw)%length_inv
  length_shk = box_list(box_shk)%length
  length_inv_shk = box_list(box_shk)%length_inv
  volume_shk = box_list(box_shk)%volume
  volume_grw = box_list(box_grw)%volume

  ! Store the old configurations of all atoms and COMs
  CALL Save_Cartesian_Coordinates_Box(box_grw)
  CALL Save_Cartesian_Coordinates_Box(box_shk)

  ! store the pair interactions
  IF (l_pair_nrg) THEN
     CALL MOVE_ALLOC(pair_nrg_vdw,pair_nrg_vdw_old)
     CALL MOVE_ALLOC(pair_nrg_qq,pair_nrg_qq_old)
     ALLOCATE(pair_nrg_vdw(sum_max_molecules_p4,sum_max_molecules))
     ALLOCATE(pair_nrg_qq(sum_max_molecules_p4,sum_max_molecules))
     !ALLOCATE(pair_nrg_vdw_old(sum_max_molecules_p4,sum_max_molecules))
     !ALLOCATE(pair_nrg_qq_old(sum_max_molecules_p4,sum_max_molecules))

     !!!$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
     !pair_nrg_vdw_old(:,:) = pair_nrg_vdw(:,:)
     !pair_nrg_qq_old(:,:) = pair_nrg_qq(:,:)
     !!!$OMP END PARALLEL WORKSHARE
  END IF

  ! store cos_mol and sin_mol
  IF ( int_charge_sum_style(box_grw) == charge_ewald .OR. &
       int_charge_sum_style(box_shk) == charge_ewald ) THEN
     CALL MOVE_ALLOC(cos_mol,cos_mol_old)
     CALL MOVE_ALLOC(sin_mol,sin_mol_old)
     
     !ALLOCATE(cos_mol_old(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)
     !
     !IF (AllocateStatus /= 0 ) THEN
     !   err_msg = ''
     !   err_msg(1) = 'Memory could not be allocated for cos_mol_old'
     !   CALL Clean_Abort(err_msg,'gemc_nvt_volume.f90')
     !END IF
     !
     !ALLOCATE(sin_mol_old(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)

     !IF (AllocateStatus /= 0 ) THEN
     !   err_msg = ''
     !   err_msg(1) = 'Memory could not be allocated for sin_mol_old'
     !   CALL Clean_Abort(err_msg,'gemc_nvt_volume.f90')
     !END IF
     !
     !!!$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
     !cos_mol_old(:,:) = cos_mol(:,:)
     !sin_mol_old(:,:) = sin_mol(:,:)
     !!!$OMP END PARALLEL WORKSHARE

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
     rcut_vdw_grw = rcut_vdw(box_grw)
     rcut_coul_grw = rcut_coul(box_grw)
     rcut_vdwsq_grw = rcut_vdwsq(box_grw)
     rcut_coulsq_grw = rcut_coulsq(box_grw)

     rcut3_grw = rcut3(box_grw)
     rcut9_grw = rcut9(box_grw)
     rcut_vdw3_grw = rcut_vdw3(box_grw)
     rcut_vdw6_grw = rcut_vdw6(box_grw)

     rcut_max_grw = rcut_max(box_grw)
     
     IF (int_charge_sum_style(box_grw) == charge_ewald) THEN
!        alpha_ewald_grw = alpha_ewald(box_grw)
        h_ewald_cut_grw = h_ewald_cut(box_grw)
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
     rcut_vdw_shk = rcut_vdw(box_shk)
     rcut_coul_shk = rcut_coul(box_shk)
     rcut_vdwsq_shk = rcut_vdwsq(box_shk)
     rcut_coulsq_shk = rcut_coulsq(box_shk)

     rcut3_shk = rcut3(box_shk)
     rcut9_shk = rcut9(box_shk)
     rcut_vdw3_shk = rcut_vdw3(box_shk)
     rcut_vdw6_shk = rcut_vdw6(box_shk)

     rcut_max_shk = rcut_max(box_shk)
     
     IF (int_charge_sum_style(box_shk) == charge_ewald) THEN
!        alpha_ewald_shk = alpha_ewald(box_shk)
        h_ewald_cut_shk = h_ewald_cut(box_shk)
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

  CALL Scale_COM_Cartesian(box_grw,length_inv_grw)
  CALL Scale_COM_Cartesian(box_shk,length_inv_shk)

  ! Now let us compute the energy change due to the combined move

  energy_grw = energy(box_grw)
  energy_shk = energy(box_shk)

  IF ( int_charge_sum_style(box_grw) == charge_ewald) THEN
     ! Then we need to determine if the number of k vectors change due to this volume move
     ! we basically follow the procedure outlined in volume_change.f90 for the two boxes.
     ! Note that the number of k vectors will increase only for the box that is expanding

     ! store old terms 

     nvecs_grw = nvecs(box_grw)
     nvecs_shk = nvecs(box_shk)
     nvecs_max = MAXVAL(nvecs)
     IF (ALLOCATED(box_list(box_grw)%sincos_sum_old)) DEALLOCATE(box_list(box_grw)%sincos_sum_old)
     IF (ALLOCATED(box_list(box_shk)%sincos_sum_old)) DEALLOCATE(box_list(box_shk)%sincos_sum_old)

     CALL MOVE_ALLOC(box_list(box_grw)%kspace_vectors,kspace_vectors_grw)
     CALL MOVE_ALLOC(box_list(box_shk)%kspace_vectors,kspace_vectors_shk)
     CALL MOVE_ALLOC(box_list(box_grw)%kspace_vector_ints,kspace_vector_ints_grw)
     CALL MOVE_ALLOC(box_list(box_shk)%kspace_vector_ints,kspace_vector_ints_shk)
     CALL MOVE_ALLOC(box_list(box_grw)%sincos_sum,box_list(box_grw)%sincos_sum_old)
     CALL MOVE_ALLOC(box_list(box_shk)%sincos_sum,box_list(box_shk)%sincos_sum_old)
     !!!$OMP PARALLEL WORKSHARE DEFAULT(SHARED) 
     !cos_sum_old(:,:) = cos_sum(:,:)
     !sin_sum_old(:,:) = sin_sum(:,:)

     !hx_grw(:) = hx(:,box_grw)
     !hy_grw(:) = hy(:,box_grw)
     !hz_grw(:) = hz(:,box_grw)
     !Cn_grw(:) = Cn(:,box_grw)

     !hx_shk(:) = hx(:,box_shk)
     !hy_shk(:) = hy(:,box_shk)
     !hz_shk(:) = hz(:,box_shk)
     !Cn_shk(:) = Cn(:,box_shk)

     !!!$OMP END PARALLEL WORKSHARE

     ! Determine the new k vectors for both boxes. The call will change Cn, hx, hy and hz and hence will
     ! change cos_sum and sin_sum.
     
     CALL Ewald_Reciprocal_Lattice_Vector_Setup(box_grw)
     CALL Ewald_Reciprocal_Lattice_Vector_Setup(box_shk)
     nvecs_max_new_p4 = IAND(MAXVAL(nvecs)+padconst_8byte,padmask_8byte)
     ALLOCATE(sin_mol(nvecs_max_new_p4,0:SUM(max_molecules)))
     ALLOCATE(cos_mol(nvecs_max_new_p4,0:SUM(max_molecules)))

        
  END IF
  
  ! Compute the energies of the boxes. For computational efficiency,
  ! we will determine the energy of the box whose volume decreased
  ! as this box will likely to have more overlaps

  CALL Compute_System_Total_Energy(box_shk, .TRUE., overlap)
  !CALL Compute_System_Total_Energy(box_shk, .FALSE., overlap)
  
  IF (overlap) THEN 
      CALL Reset_Coords

      IF (verbose_log) THEN
         WRITE(logunit,'(X,I19,X,A10,X,5X,X,3X,X,I1,A1,I1,X,L8,X,9X,X,A9)') &
               i_mcstep, 'vol_swap', box_shk, '>', box_grw, accept, 'overlap'
      END IF

  ELSE
      
      CALL Compute_System_Total_Energy(box_grw, .TRUE.,overlap)
      !CALL Compute_System_Total_Energy(box_grw, .FALSE.,overlap)
     
      ! actually there should be no overlap for the box whose dimensions
      ! increase but we will include this check only for safety.
      
      IF (overlap) THEN
         CALL Reset_Coords

         IF (verbose_log) THEN
            WRITE(logunit,'(X,I19,X,A10,X,5X,X,3X,X,I1,A1,I1,X,L8,X,9X,X,A9)') &
                  i_mcstep, 'vol_swap', box_shk, '>', box_grw, accept, 'overlap'
         END IF

      ELSE
         ! accept or reject the move based on the acceptance rule
         
         dE_grw = energy(box_grw)%total - energy_grw%total
         dE_shk = energy(box_shk)%total - energy_shk%total
         !dE_grw = energy(box_grw)%inter - energy_grw%inter
         !dE_shk = energy(box_shk)%inter - energy_shk%inter
         
         IF (f_dv) THEN
            ln_pacc = beta(box_grw) * dE_grw + beta(box_shk) * dE_shk &
                    - REAL(SUM(nmols(:,box_grw)),DP) * DLOG( box_list(box_grw)%volume / volume_grw) &
                    - REAL(SUM(nmols(:,box_shk)),DP) * DLOG( box_list(box_shk)%volume / volume_shk)
            
         ELSE IF(f_vratio) THEN
            
            ln_pacc = beta(box_grw) * dE_grw + beta(box_shk) * dE_shk &
                    - REAL(SUM(nmols(:,box_grw))+1,DP) * DLOG( box_list(box_grw)%volume / volume_grw) &
                    - REAL(SUM(nmols(:,box_shk))+1,DP) * DLOG( box_list(box_shk)%volume / volume_shk)
            
         END IF

         accept = accept_or_reject(ln_pacc)
         
         IF (accept) THEN
            
            nvol_success(box_grw) = nvol_success(box_grw) + 1
            nvol_success(box_shk) = nvol_success(box_shk) + 1
            ivol_success(box_grw) = ivol_success(box_grw) + 1
            ivol_success(box_shk) = ivol_success(box_shk) + 1
            ! energy,positions and box dimensions are already updated
            IF ((int_charge_sum_style(box_grw) == charge_ewald .OR. l_pair_nrg) .AND. nbr_boxes>2) THEN
               
               !! cos_sum and sin_sum were deallocated, destroying the terms
               !! for boxes other than box_grw, box_shk. so restore these
               !
               !DO ibox = 1, nbr_boxes
               !   
               !   IF ( .NOT. ( (ibox /= box_grw) .OR. (ibox /= box_shk))) THEN
               !      ! transfer cos_sum and sin_sum for other boxes
               !      ! Note that direct assignment of cos_sum_old to cos_sum
               !      ! will result into an error as these two arrays are of
               !      ! different dimensions
               !      cos_sum(1:nvecs(ibox),ibox) = cos_sum_old(1:nvecs(ibox),ibox)
               !      sin_sum(1:nvecs(ibox),ibox) = sin_sum_old(1:nvecs(ibox),ibox)
               !      
               !   END IF
               !   
               !END DO
               !! Now deallocate cos_sum_old and sin_sum_old so that they have the same dimensions
               !! as sin_sum and cos_sum                  
               !DEALLOCATE(cos_sum_old,sin_sum_old)
               !ALLOCATE(cos_sum_old(SIZE(cos_sum,1),nbr_boxes))
               !ALLOCATE(sin_sum_old(SIZE(sin_sum,1),nbr_boxes))
               
               ! Now assign cos_mol and sin_mol for the molecules present in other
               ! boxes. Note that cos_mol for box_grw and box_shk have been assigned
               ! when a call to Compute_System_Total_Energy was placed.
               
               DO ibox = 1, nbr_boxes
                  IF (ibox == box_grw .OR. ibox == box_shk) CYCLE
                  nboxmols = SUM(nmols(:,ibox))
                  IF (nboxmols < 1) CYCLE
                  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(is,pos,im_locate_shift,iend,istart)
                  !$OMP DO SCHEDULE(STATIC)
                  DO is = 1, nspecies
                     IF (is==1) THEN
                             im_locate_shift = 0
                     ELSE
                             im_locate_shift = SUM(max_molecules(1:is-1))
                     END IF
                     iend = SUM(nmols(1:is,ibox))
                     istart = iend - nmols(is,ibox) + 1
                     posvec(istart:iend) = im_locate_shift + locate(1:nmols(is,ibox),is,ibox)
                  END DO
                  !$OMP END DO
                  !$OMP DO SCHEDULE(STATIC)
                  DO im = 1, nboxmols
                        pos = posvec(im)
                        IF (int_charge_sum_style(ibox) == charge_ewald) THEN
                                !DIR$ VECTOR ALIGNED
                                cos_mol(1:nvecs(ibox),pos) = cos_mol_old(1:nvecs(ibox),pos)
                                !DIR$ VECTOR ALIGNED
                                sin_mol(1:nvecs(ibox),pos) = sin_mol_old(1:nvecs(ibox),pos)
                                cos_mol(nvecs(ibox)+1:,pos) = 0.0_DP
                                sin_mol(nvecs(ibox)+1:,pos) = 0.0_DP
                        END IF
                        IF (l_pair_nrg) THEN
                                ! Technically unnecessary but generally faster (and harmless) to copy whole column 
                                ! instead of only the elements in posvec
                                !DIR$ VECTOR ALIGNED
                                pair_nrg_vdw(:,pos) = pair_nrg_vdw_old(:,pos)
                                !DIR$ VECTOR ALIGNED
                                pair_nrg_qq(:,pos) = pair_nrg_qq_old(:,pos)
                        END IF
                  END DO
                  !$OMP END DO
                  !$OMP END PARALLEL
               END DO
            END IF
            
  
         ELSE
            
            CALL Reset_Coords

         END IF
        
         IF (verbose_log) THEN
            WRITE(logunit,'(X,I19,X,A10,X,5X,X,3X,X,I1,A1,I1,X,L8,X,9X,X,F9.3)') &
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

         WRITE(logunit,'(X,I19,X,A10,X,5X,X,3X,X,3X,X,F8.5,X,F9.0)') &
               i_mcstep, 'vol_swap', success_ratio, box_list(box_grw)%dv_max

      ELSE

         success_ratio = REAL(nvol_success(box_grw),DP)/REAL(nvolumes(box_grw),DP)
         WRITE(logunit,'(X,I19,X,A10,X,5X,X,3X,X,3X,X,F8.5)') &
               i_mcstep, 'vol_swap', success_ratio
         
      END IF

  END IF

   
 CONTAINS
   
   SUBROUTINE Reset_Coords
     
     IMPLICIT NONE
     
     CALL Reset_Cartesian_Coordinates_Box(box_grw)
     CALL Reset_Cartesian_Coordinates_Box(box_shk)
     
     ! box list and energy
     box_list(box_grw)%length = length_grw
     box_list(box_shk)%length = length_shk

     CALL Compute_Cell_Dimensions(box_grw)
     CALL Compute_Cell_Dimensions(box_shk)
     
     !box_list(box_grw) = box_list_grw
     !box_list(box_shk) = box_list_shk
     
     energy(box_grw) = energy_grw
     energy(box_shk) = energy_shk

     IF (l_pair_nrg) THEN
             !pair_nrg_vdw(:,:) = pair_nrg_vdw_old(:,:)
             !pair_nrg_qq(:,:) = pair_nrg_qq_old(:,:)
             !DEALLOCATE(pair_nrg_vdw_old,pair_nrg_qq_old)
        DEALLOCATE(pair_nrg_vdw,pair_nrg_qq)
        CALL MOVE_ALLOC(pair_nrg_vdw_old,pair_nrg_vdw)
        CALL MOVE_ALLOC(pair_nrg_qq_old,pair_nrg_qq)
     END IF

     IF (l_half_len_cutoff(box_grw)) THEN
        
        rcut_vdw(box_grw) = rcut_vdw_grw
        rcut_coul(box_grw) = rcut_coul_grw
        rcut_vdwsq(box_grw) = rcut_vdwsq_grw
        rcut_coulsq(box_grw) = rcut_coulsq_grw
        
        rcut3(box_grw) = rcut3_grw
        rcut9(box_grw) = rcut9_grw
        rcut_vdw3(box_grw) = rcut_vdw3_grw
        rcut_vdw6(box_grw) = rcut_vdw6_grw
        
        rcut_max(box_grw) = rcut_max_grw
        
        IF( int_charge_sum_style(box_grw) == charge_ewald ) THEN
           
           !           alpha_ewald(box_grw) = alpha_ewald_grw
           h_ewald_cut(box_grw) = h_ewald_cut_grw
           
        END IF
        
        
     END IF

     IF (l_half_len_cutoff(box_shk)) THEN
        
        rcut_vdw(box_shk) = rcut_vdw_shk
        rcut_coul(box_shk) = rcut_coul_shk
        rcut_vdwsq(box_shk) = rcut_vdwsq_shk
        rcut_coulsq(box_shk) = rcut_coulsq_shk
        
        rcut3(box_shk) = rcut3_shk
        rcut9(box_shk) = rcut9_shk
        rcut_vdw3(box_shk) = rcut_vdw3_shk
        rcut_vdw6(box_shk) = rcut_vdw6_shk
        
        rcut_max(box_shk) = rcut_max_shk

        IF( int_charge_sum_style(box_shk) == charge_ewald ) THEN
           
!           alpha_ewald(box_shk) = alpha_ewald_shk
           h_ewald_cut(box_shk) = h_ewald_cut_shk
           
        END IF
        
        
     END IF
     
     IF (int_charge_sum_style(box_grw) == charge_ewald) THEN
        
       nvecs(box_grw) = nvecs_grw
       nvecs(box_shk) = nvecs_shk
       DEALLOCATE(box_list(box_shk)%kspace_vectors)
       DEALLOCATE(box_list(box_grw)%kspace_vectors)
       DEALLOCATE(box_list(box_shk)%kspace_vector_ints)
       DEALLOCATE(box_list(box_grw)%kspace_vector_ints)
       DEALLOCATE(box_list(box_shk)%sincos_sum)
       DEALLOCATE(box_list(box_grw)%sincos_sum)
       DEALLOCATE(sin_mol,cos_mol)
       CALL MOVE_ALLOC(sin_mol_old,sin_mol)
       CALL MOVE_ALLOC(cos_mol_old,cos_mol)
       CALL MOVE_ALLOC(kspace_vectors_grw,box_list(box_grw)%kspace_vectors)
       CALL MOVE_ALLOC(kspace_vectors_shk,box_list(box_shk)%kspace_vectors)
       CALL MOVE_ALLOC(kspace_vector_ints_grw,box_list(box_grw)%kspace_vector_ints)
       CALL MOVE_ALLOC(kspace_vector_ints_shk,box_list(box_shk)%kspace_vector_ints)
       CALL MOVE_ALLOC(box_list(box_grw)%sincos_sum_old,box_list(box_grw)%sincos_sum)
       CALL MOVE_ALLOC(box_list(box_shk)%sincos_sum_old,box_list(box_shk)%sincos_sum)
       
       
     END IF
    
  END SUBROUTINE Reset_Coords

END SUBROUTINE GEMC_NVT_Volume


