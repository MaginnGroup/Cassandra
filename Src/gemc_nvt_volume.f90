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

SUBROUTINE GEMC_NVT_Volume(box1, box2)

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

  USE Run_Variables
  USE Random_Generators
  USE Volume
  USE IO_Utilities
  USE Energy_Routines
  USE Pair_Nrg_Routines

  IMPLICIT NONE


  INTEGER, INTENT(OUT) :: box1, box2

  ! Local variables

  INTEGER :: tot_mol_box1, tot_mol_box2, first_box, second_box
  INTEGER :: ibox, nvecs_old_1, nvecs_old_2, nvecs_max
  INTEGER :: is, im

  REAL(DP) :: delta_volume, factor, delta_e_1, delta_e_2
  REAL(DP) :: success_ratio
  REAL(DP), DIMENSION(maxk) :: hx_old_1, hy_old_1, hz_old_1, Cn_old_1
  REAL(DP), DIMENSION(maxk) :: hx_old_2, hy_old_2, hz_old_2, Cn_old_2  
  REAL(DP) :: v_ratio_o, v_total, vol_factor

  LOGICAL :: overlap, accept, accept_or_reject, allocation_cos_sin

  TYPE(Box_Class) :: box_list_old_1, box_list_old_2
  TYPE(Energy_Class) :: energy_old_1, energy_old_2

  REAL(DP) :: checke
  LOGICAL :: superbad

  INTEGER :: position, alive, my_box

  REAL(DP), ALLOCATABLE :: pair_nrg_vdw_old(:,:), pair_nrg_qq_old(:,:)
  REAL(DP), ALLOCATABLE :: cos_mol_old(:,:), sin_mol_old(:,:)

  REAL(DP) :: rcut_vdw_old_1, rcut_coul_old_1, rcut3_old_1, rcut9_old_1, alpha_ewald_old_1
  REAL(DP) :: h_ewald_cut_old_1, rcut_vdwsq_old_1, rcut_coulsq_old_1, rcut_vdw3_old_1
  REAL(DP) :: rcut_vdw6_old_1, rcut_max_old_1

  REAL(DP) :: rcut_vdw_old_2, rcut_coul_old_2, rcut3_old_2, rcut9_old_2, alpha_ewald_old_2
  REAL(DP) :: h_ewald_cut_old_2, rcut_vdwsq_old_2, rcut_coulsq_old_2, rcut_vdw3_old_2
  REAL(DP) :: rcut_vdw6_old_2, rcut_max_old_2

  IF (f_dv) THEN
     ! Pick a box at random
     box1 = INT (rranf() * nbr_boxes) + 1
     
     
     ! Obtain the identity of the other box
     
     DO WHILE(.TRUE.)
        
        box2 = INT ( rranf() * nbr_boxes) + 1
        
        IF ( box2 /= box1 ) EXIT
        
     END DO

  ELSE IF (f_vratio) THEN

     box1 = 1
     box2 = 2

  END IF

  tot_trials(box1) = tot_trials(box1) + 1
  tot_trials(box2) = tot_trials(box2) + 1

  nvolumes(box1) = nvolumes(box1) + 1
  nvolumes(box2) = nvolumes(box2) + 1

  ! store old cell matrix 

  box_list_old_1 = box_list(box1)
  box_list_old_2 = box_list(box2)

  ! Store the old configurations of all atoms and COMs

  CALL Save_Cartesian_Coordinates_Box(box1,tot_mol_box1)
  CALL Save_Cartesian_Coordinates_Box(box2,tot_mol_box2)

  ! store the pair interactions

  IF (l_pair_nrg) THEN
     ALLOCATE(pair_nrg_vdw_old(SUM(nmolecules),SUM(nmolecules)))
     ALLOCATE(pair_nrg_qq_old(SUM(nmolecules),SUM(nmolecules)))

     pair_nrg_vdw_old(:,:) = pair_nrg_vdw(:,:)
     pair_nrg_qq_old(:,:) = pair_nrg_qq(:,:)

  END IF

  ! store cos_mol and sin_mol

  IF ( int_charge_sum_style(box1) == charge_ewald .OR. &
       int_charge_sum_style(box2) == charge_ewald ) THEN
     
     ALLOCATE(cos_mol_old(MAXVAL(nvecs),SUM(nmolecules)), Stat = AllocateStatus)
     
     IF (AllocateStatus /= 0 ) THEN
        err_msg = ''
        err_msg(1) = 'Memory could not be allocated for cos_mol_old'
        CALL Clean_Abort(err_msg,'gemc_nvt_volume.f90')
     END IF
     
     ALLOCATE(sin_mol_old(MAXVAL(nvecs),SUM(nmolecules)), Stat = AllocateStatus)

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
     
     delta_volume = (1.0_DP - 2.0_DP * rranf()) * box_list(box1)%dv_max
     
     ! update each of the box volumes
     
     box_list(box1)%volume = box_list(box1)%volume + delta_volume
     box_list(box2)%volume = box_list(box2)%volume - delta_volume
     
  ELSE IF (f_vratio) THEN

!!$  ! new volumes
     v_ratio_o = (box_list(box1)%volume/box_list(box2)%volume)
     v_total = box_list(box1)%volume + box_list(box2)%volume

     delta_volume =  (1.0_DP - 2.0_DP * rranf()) * box_list(box1)%dv_max
     vol_factor = v_ratio_o * EXP(delta_volume)
     box_list(box1)%volume = v_total * vol_factor / (1.0_DP + vol_factor)
     box_list(box2)%volume = v_total - box_list(box1)%volume
  
  END IF

  IF ( box_list(box1)%int_box_shape == int_cubic ) THEN

     box_list(box1)%length(1,1) = (box_list(box1)%volume) ** (1.0_DP/3.0_DP)
     box_list(box1)%length(2,2) = box_list(box1)%length(1,1)
     box_list(box1)%length(3,3) = box_list(box1)%length(2,2)

  ELSE

     err_msg = ''
     err_msg = 'Noncubic shape encountered for the box'
     err_msg = Int_To_String(box1)
     CALL Clean_Abort(err_msg, 'GEMC_NVT_Volume')

  END IF
  
  IF ( box_list(box2)%int_box_shape == int_cubic ) THEN

     box_list(box2)%length(1,1) = (box_list(box2)%volume) ** (1.0_DP/3.0_DP)
     box_list(box2)%length(2,2) = box_list(box2)%length(1,1)
     box_list(box2)%length(3,3) = box_list(box2)%length(2,2)

  ELSE

     err_msg = ''
     err_msg = 'Noncubic shape encountered for the box'
     err_msg = Int_To_String(box2)
     CALL Clean_Abort(err_msg, 'GEMC_NVT_Volume')

  END IF

  CALL Compute_Cell_Dimensions(box1)
  CALL Compute_Cell_Dimensions(box2)

  IF ( l_half_len_cutoff(box1)) THEN
     IF ( box_list(box1)%int_box_shape == int_cubic ) THEN

        ! store old cutoffs and other associated quantities

        rcut_vdw_old_1 = rcut_vdw(box1)
        rcut_coul_old_1 = rcut_coul(box1)
        rcut_vdwsq_old_1 = rcut_vdwsq(box1)
        rcut_coulsq_old_1 = rcut_coulsq(box1)

        rcut3_old_1 = rcut3(box1)
        rcut9_old_1 = rcut9(box1)
        rcut_vdw3_old_1 = rcut_vdw3(box1)
        rcut_vdw6_old_1 = rcut_vdw6(box1)

        rcut_max_old_1 = rcut_max(box1)
        
        IF (int_charge_sum_style(box1) == charge_ewald) THEN
           
!           alpha_ewald_old_1 = alpha_ewald(box1)
           h_ewald_cut_old_1 = h_ewald_cut(box1)

        END IF

        rcut_vdw(box1) = 0.5_DP * box_list(box1)%length(1,1)
        rcut_coul(box1) = rcut_vdw(box1)
        rcut_vdwsq(box1) = rcut_vdw(box1) * rcut_vdw(box1)
        rcut_coulsq(box1) = rcut_vdwsq(box1)

        rcut_vdw3(box1) = rcut_vdwsq(box1) * rcut_vdw(box1)
        rcut_vdw6(box1) = rcut_vdw3(box1) * rcut_vdw3(box1)
        rcut3(box1) = rcut_vdw3(box1)
        rcut9(box1) = rcut3(box1) * rcut_vdw6(box1)

        rcut_max(box1) = rcut_vdw(box1)

        IF (int_charge_sum_style(box1) == charge_ewald) THEN

!           alpha_ewald(box1) = ewald_p_sqrt(box1) / rcut_coul(box1)
           h_ewald_cut(box1) = 2.0_DP * ewald_p(box1) / rcut_coul(box1)

        END IF

     END IF

  ELSE

     IF ( box_list(box1)%int_box_shape == int_cubic ) THEN
        IF ( 0.5_DP * box_list(box1)%length(1,1) < rcut_vdw(box1) .OR. &
             0.5_DP * box_list(box1)%length(1,1) < rcut_coul(box1) .OR. &
             0.5_DP * box_list(box1)%length(1,1) < roff_charmm(box1) .OR. &
             0.5_DP * box_list(box1)%length(1,1) < roff_switch(box1) ) THEN
           err_msg = ''
           err_msg(1) = 'Cutoff is greater than the half box length'
           err_msg(2) = Int_to_String(box1)
           CALL Clean_Abort(err_msg,'GEMC_NVT_Volume')
        END IF
     END IF
     
  END IF

  IF ( l_half_len_cutoff(box2)) THEN
     IF ( box_list(box2)%int_box_shape == int_cubic ) THEN

        ! store old cutoffs and other associated quantities

        rcut_vdw_old_2 = rcut_vdw(box2)
        rcut_coul_old_2 = rcut_coul(box2)
        rcut_vdwsq_old_2 = rcut_vdwsq(box2)
        rcut_coulsq_old_2 = rcut_coulsq(box2)

        rcut3_old_2 = rcut3(box2)
        rcut9_old_2 = rcut9(box2)
        rcut_vdw3_old_2 = rcut_vdw3(box2)
        rcut_vdw6_old_2 = rcut_vdw6(box2)

        rcut_max_old_2 = rcut_max(box2)
        
        IF (int_charge_sum_style(box2) == charge_ewald) THEN
           
!           alpha_ewald_old_2 = alpha_ewald(box2)
           h_ewald_cut_old_2 = h_ewald_cut(box2)

        END IF

        rcut_vdw(box2) = 0.5_DP * box_list(box2)%length(1,1)
        rcut_coul(box2) = rcut_vdw(box2)
        rcut_vdwsq(box2) = rcut_vdw(box2) * rcut_vdw(box2)
        rcut_coulsq(box2) = rcut_vdwsq(box2)

        rcut_vdw3(box2) = rcut_vdwsq(box2) * rcut_vdw(box2)
        rcut_vdw6(box2) = rcut_vdw3(box2) * rcut_vdw3(box2)
        rcut3(box2) = rcut_vdw3(box2)
        rcut9(box2) = rcut3(box2) * rcut_vdw6(box2)

        rcut_max(box2) = rcut_vdw(box2)
        
        IF (int_charge_sum_style(box2) == charge_ewald) THEN

!           alpha_ewald(box2) = ewald_p_sqrt(box2) / rcut_coul(box2)
           h_ewald_cut(box2) = 2.0_DP * ewald_p(box2) / rcut_coul(box2)

        END IF

     END IF
     

  ELSE
     
     IF ( box_list(box2)%int_box_shape == int_cubic ) THEN
        IF ( 0.5_DP * box_list(box2)%length(1,1) < rcut_vdw(box2) .OR. &
             0.5_DP * box_list(box2)%length(1,1) < rcut_coul(box2) .OR. &
             0.5_DP * box_list(box2)%length(1,1) < roff_charmm(box2) .OR. &
             0.5_DP * box_list(box2)%length(1,1) < roff_switch(box2) ) THEN
           err_msg = ''
           err_msg(1) = 'Cutoff is greater than the half box length'
           err_msg(2) = Int_To_String(box2)
           CALL Clean_Abort(err_msg,'GEMC_NVT_Volume')
        END IF
     END IF
     
  END IF


  ! Rescale the COM and all atomic positions

  CALL Scale_COM_Cartesian(box1,box_list_old_1)
  CALL Scale_COM_Cartesian(box2,box_list_old_2)

  ! Now let us compute the energy change due to the combined move

  energy_old_1 = energy(box1)
  energy_old_2 = energy(box2)

  IF ( int_charge_sum_style(box1) == charge_ewald) THEN
     ! Then we need to determine if the number of k vectors change due to this volume move
     ! we basically follow the procedure outlined in volume_change.f90 for the two boxes.
     ! Note that the number of k vectors will increase only for the box that is expanding

     ! store old terms 

     nvecs_old_1 = nvecs(box1)
     nvecs_old_2 = nvecs(box2)
     nvecs_max = MAXVAL(nvecs)

     !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED) 
     cos_sum_old(:,:) = cos_sum(:,:)
     sin_sum_old(:,:) = sin_sum(:,:)

     hx_old_1(:) = hx(:,box1)
     hy_old_1(:) = hy(:,box1)
     hz_old_1(:) = hz(:,box1)
     Cn_old_1(:) = Cn(:,box1)

     hx_old_2(:) = hx(:,box2)
     hy_old_2(:) = hy(:,box2)
     hz_old_2(:) = hz(:,box2)
     Cn_old_2(:) = Cn(:,box2)

     !!$OMP END PARALLEL WORKSHARE

     ! Determine the new k vectors for this box. The call will change Cn, hx, hy and hz and hence will
     ! change cos_sum and sin_sum.
     
     CALL Ewald_Reciprocal_Lattice_Vector_Setup(box1)
     CALL Ewald_Reciprocal_Lattice_Vector_Setup(box2)

     ! check to see if we need to reallocate arrays

     allocation_cos_sin = .FALSE.

     IF ((nvecs(box1) > SIZE(cos_sum,1)) .OR. (nvecs(box2) > SIZE(cos_sum,1))) THEN
        ! we need to expand the cos_sum and sin_sum arrays
        allocation_cos_sin = .TRUE.

        DEALLOCATE(cos_sum,sin_sum)

        IF (ALLOCATED(cos_mol)) DEALLOCATE(cos_mol)
        IF (ALLOCATED(sin_mol)) DEALLOCATE(sin_mol)

        ALLOCATE(cos_sum(MAXVAL(nvecs),nbr_boxes),Stat=AllocateStatus)
        
        IF (AllocateStatus /= 0) THEN
           err_msg = ''
           err_msg(1) = 'Memory could not be allocated for cos_sum'
           err_msg(2) = 'allocation_cos_sin'
           CALL Clean_Abort(err_msg,'GEMC_NVT_VOLUME')
        END IF
        
        
        ALLOCATE(sin_sum(MAXVAL(nvecs),nbr_boxes), Stat = AllocateStatus)
        
        IF (Allocatestatus /= 0) THEN
           err_msg = ''
           err_msg(1) = 'Memory could not be allocated for sin_sum'
           err_msg(2) = 'allocation_cos_sin'
           CALL Clean_Abort(err_msg,'GEMC_NVT_VOLUME')
        END IF

        ALLOCATE(cos_mol(MAXVAL(nvecs),SUM(nmolecules)), Stat = AllocateStatus)
        
        IF (AllocateStatus /= 0) THEN
           err_msg = ''
           err_msg(1) = 'Memory could not be allocated for cos_mol'
           err_msg(2) = 'allocation_cos_sin'
           CALL Clean_Abort(err_msg,'gemc_volume_change.f90')
        END IF
        
        ALLOCATE(sin_mol(MAXVAL(nvecs),SUM(nmolecules)), Stat = AllocateStatus)
        
        IF (AllocateStatus /= 0) THEN
           err_msg = ''
           err_msg(1) = 'Memory could not be allocated for sin_mol'
           err_msg(2) = 'allocation_cos_sin'
           CALL Clean_Abort(err_msg,'gemc_volume_change.f90')
        END IF
        
     END IF
     
  END IF
  
  ! Compute the energies of the boxes. For computational efficiency,
  ! we will determine the energy of the box whose volume decreased
  ! as this box will likely to have more overlaps

  IF (box_list(box1)%volume < box_list_old_1%volume) THEN
     ! box1 is shrinking
     first_box = box1
     second_box = box2
  ELSE
     ! else box2 is shrinking
     first_box = box2
     second_box = box1
  END IF

  CALL Compute_Total_System_Energy(first_box, .FALSE., overlap)
  
   IF (overlap) THEN 

      CALL Reset_Coords

   ELSE
      
      CALL Compute_Total_System_Energy(second_box, .FALSE.,overlap)
     
      ! actually there should be no overlap for the box whose dimensions
      ! increase but we will include this check only for safety.
      
      IF (overlap) THEN
         CALL Reset_Coords
      ELSE
         ! accept or reject the move based on the acceptance rule
         
         delta_e_1 = energy(box1)%total - energy_old_1%total
         delta_e_2 = energy(box2)%total - energy_old_2%total
         
         IF (f_dv) THEN
            factor = beta(box1) * delta_e_1 + beta(box2) * delta_e_2 - &
                 tot_mol_box1 * DLOG( box_list(box1)%volume / box_list_old_1%volume) - &
                 tot_mol_box2 * DLOG( box_list(box2)%volume / box_list_old_2%volume)
            
         ELSE IF(f_vratio) THEN
            
            factor = beta(box1) * delta_e_1 + beta(box2) * delta_e_2 - &
                 REAL((tot_mol_box1 + 1),DP) * DLOG( box_list(box1)%volume / box_list_old_1%volume) - &
                 REAL((tot_mol_box2 + 1),DP) * DLOG( box_list(box2)%volume / box_list_old_2%volume)
            
         END IF

         accept = accept_or_reject(factor)
         
         IF (accept) THEN
            
            nvol_success(box1) = nvol_success(box1) + 1
            nvol_success(box2) = nvol_success(box2) + 1
            ivol_success(box1) = ivol_success(box1) + 1
            ivol_success(box2) = ivol_success(box2) + 1
            ! energy,positions and box dimensions are already updated
            IF (int_charge_sum_style(box1) == charge_ewald) THEN
               
               IF (allocation_cos_sin) THEN
                  
                  ! cos_sum and sin_sum were deallocated, destroying the terms
                  ! for boxes other than box1, box2. so restore these
                  
                  DO ibox = 1, nbr_boxes
                     
                     IF ( .NOT. ( (ibox /= box1) .OR. (ibox /= box2))) THEN
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
                  ! boxes. Note that cos_mol for box1 and box2 have been assigned
                  ! when a call to Compute_Total_System_Energy was placed.
                  
                  DO is = 1, nspecies
                     DO im = 1, nmolecules(is)
                        alive = locate(im,is)
                        
                        IF (.NOT. molecule_list(alive,is)%live) CYCLE
                        
                        my_box = molecule_list(alive,is)%which_box
                        
                        IF (my_box == box1 .OR. my_box == box2) CYCLE
                        
                        CALL Get_Position_Alive(alive,is,position)
                        
                        !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
                        cos_mol(1:nvecs(my_box),position) = cos_mol_old(1:nvecs(my_box),position)
                        sin_mol(1:nvecs(my_box),position) = sin_mol_old(1:nvecs(my_box),position)
                        
                        cos_mol(nvecs(my_box)+1:MAXVAL(nvecs),position) = 0.0_DP
                        sin_mol(nvecs(my_box)+1:MAXVAL(nvecs),position) = 0.0_DP
                        !!$OMP END PARALLEL WORKSHARE
                        
                     END DO
                  END DO
               END IF
               DEALLOCATE(cos_mol_old,sin_mol_old)
            END IF  ! ends if of 316
            
            IF (l_pair_nrg) DEALLOCATE(pair_nrg_vdw_old,pair_nrg_qq_old)
  
         ELSE
            
            CALL Reset_Coords

         END IF
        
      END IF

   END IF

  ! Update the maximum volume modulus of equilibration runs
   
   IF (MOD(nvolumes(box1),nvol_update) == 0 ) THEN

      IF ( int_run_style == run_equil) THEN

         success_ratio = REAL(ivol_success(box1),DP)/REAL(nvol_update,DP)
         ivol_success(box1) = 0
         ivol_success(box2) = 0
      
         IF (box_list(box1)%int_box_shape == int_cubic ) THEN
            IF ( success_ratio < 0.0001 ) THEN
               ! decrease the maximum displacement by 5%
               
               box_list(:)%dv_max = 0.1_DP * box_list(:)%dv_max
               
            ELSE
               
               box_list(:)%dv_max = 2.0_DP*success_ratio * box_list(:)%dv_max

            END IF
            
         END IF

      ELSE

         success_ratio = REAL(nvol_success(box1),DP)/REAL(nvolumes(box1),DP)
         
      END IF

      WRITE(logunit,*)
      WRITE(logunit,'(A35,2X,F10.4,2X,A7,2X,I2)')'Successful volume attempt ratio is ', &
           success_ratio, 'for box', box1
      
   END IF
   
 CONTAINS
   
   SUBROUTINE Reset_Coords
     
     IMPLICIT NONE
     
     CALL Reset_Cartesian_Coordinates_Box(box1)
     CALL Reset_Cartesian_Coordinates_Box(box2)
     
     ! box list and energy 
     
     box_list(box1) = box_list_old_1
     box_list(box2) = box_list_old_2
     
     energy(box1) = energy_old_1
     energy(box2) = energy_old_2

     IF (l_pair_nrg) THEN
        
        !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        pair_nrg_vdw(:,:) = pair_nrg_vdw_old(:,:)
        pair_nrg_qq(:,:) = pair_nrg_qq_old(:,:)
        !!$OMP END PARALLEL WORKSHARE
        
        DEALLOCATE(pair_nrg_vdw_old,pair_nrg_qq_old)
     END IF

     IF (l_half_len_cutoff(box1)) THEN
        
        rcut_vdw(box1) = rcut_vdw_old_1
        rcut_coul(box1) = rcut_coul_old_1
        rcut_vdwsq(box1) = rcut_vdwsq_old_1
        rcut_coulsq(box1) = rcut_coulsq_old_1
        
        rcut3(box1) = rcut3_old_1
        rcut9(box1) = rcut9_old_1
        rcut_vdw3(box1) = rcut_vdw3_old_1
        rcut_vdw6(box1) = rcut_vdw6_old_1
        
        rcut_max(box1) = rcut_max_old_1
        
        IF( int_charge_sum_style(box1) == charge_ewald ) THEN
           
           !           alpha_ewald(box1) = alpha_ewald_old_1
           h_ewald_cut(box1) = h_ewald_cut_old_1
           
        END IF
        
        
     END IF

     IF (l_half_len_cutoff(box2)) THEN
        
        rcut_vdw(box2) = rcut_vdw_old_2
        rcut_coul(box2) = rcut_coul_old_2
        rcut_vdwsq(box2) = rcut_vdwsq_old_2
        rcut_coulsq(box2) = rcut_coulsq_old_2
        
        rcut3(box2) = rcut3_old_2
        rcut9(box2) = rcut9_old_2
        rcut_vdw3(box2) = rcut_vdw3_old_2
        rcut_vdw6(box2) = rcut_vdw6_old_2
        
        rcut_max(box2) = rcut_max_old_2

        IF( int_charge_sum_style(box2) == charge_ewald ) THEN
           
!           alpha_ewald(box2) = alpha_ewald_old_2
           h_ewald_cut(box2) = h_ewald_cut_old_2
           
        END IF
        
        
     END IF
     
     IF (int_charge_sum_style(box1) == charge_ewald) THEN
        
        nvecs(box1) = nvecs_old_1
        nvecs(box2) = nvecs_old_2
        
        IF ( allocation_cos_sin) THEN
           
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
           
           ALLOCATE(cos_mol(MAXVAL(nvecs),SUM(nmolecules)),Stat = Allocatestatus)
           IF (Allocatestatus /= 0) THEN
              err_msg = ''
              err_msg(1) = 'Memory could not be allocated for cos_mol in the volume rejection'
              CALL Clean_Abort(err_msg,'GEMC NVT Volume_Change')
           END IF
           
           ALLOCATE(sin_mol(MAXVAL(nvecs),SUM(nmolecules)),Stat = Allocatestatus)
          IF (Allocatestatus /= 0) THEN
             err_msg = ''
             err_msg(1) = 'Memory could not be allocated for sin_mol in the volume rejection'
             CALL Clean_Abort(err_msg,'GEMC NVT Volume_Change')
          END IF
          
       END IF
       
       !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
       
       cos_sum(:,:) = cos_sum_old(:,:)
       sin_sum(:,:) = sin_sum_old(:,:)
       
       cos_mol(1:SIZE(cos_mol_old,1),:) = cos_mol_old(:,:)
       sin_mol(1:SIZE(sin_mol_old,1),:) = sin_mol_old(:,:)
       
       hx(:,box1) = hx_old_1(:)
       hy(:,box1) = hy_old_1(:)
       hz(:,box1) = hz_old_1(:)
       Cn(:,box1) = Cn_old_1(:)
       
       hx(:,box2) = hx_old_2(:)
       hy(:,box2) = hy_old_2(:)
       hz(:,box2) = hz_old_2(:)
       Cn(:,box2) = Cn_old_2(:)

       !!$OMP END PARALLEL WORKSHARE
       
       DEALLOCATE(cos_mol_old,sin_mol_old)
       
       ! here we make sure that cos_sum_old and sin_sum_old have the same dimensions
       ! as cos_sum and sin_sum
       
       DEALLOCATE(cos_sum_old,sin_sum_old)
       ALLOCATE(cos_sum_old(SIZE(cos_sum,1),nbr_boxes),sin_sum_old(SIZE(sin_sum,1),nbr_boxes))
       
    END IF
    
  END SUBROUTINE Reset_Coords

END SUBROUTINE GEMC_NVT_Volume


