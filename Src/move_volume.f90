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
SUBROUTINE Volume_Change
!*******************************************************************************

  !*****************************************************************************
  ! The subroutine performs a volume perturbation move. Presently, the routine
  ! is set up to perform volume changes in cubic simulation box. Extension
  ! to nonorthorhombic cells will be undertaken later.
  !
  ! Written by Jindal Shah 
  !
  ! Called by
  !
  !   NPTMC_Driver
  !   GEMC_Driver
  !
  ! Calls
  !
  !   Save_Old_Coordinates
  !   Compute_Cell_Dimensions
  !   Ewald_Reciprocal_Lattice_Vector_Setup
  !   Clean_Abort
  !   Compute_System_Total_Energy
  !   Revert_Old_Cartesian_Coordinates
  !   
  ! Revision History
  !
  ! 12/10/13 : Beta release
  !*****************************************************************************
  ! The subroutine is based on the following algorithm
  !
  ! Change the volume of the box by an amount delta V
  !
  ! Rescale the COM coordinates of all the molecules and not the actual
  ! atomic positions. This is necessary to ensure that the intramolecular
  ! energy of molecules do not contribute to the acceptance rule.
  !
  ! Compute the change in energy due to this move and accept the move
  ! with standard metropolis criterion.
  !*****************************************************************************
  
  USE Type_Definitions
  USE Global_Variables
  USE File_Names
  USE Random_Generators
  USE Energy_Routines
  USE IO_Utilities
  
  
  IMPLICIT NONE
  
  !  !$ include 'omp_lib.h'
  
  INTEGER :: is, im, alive, this_box, i, total_molecules, nvecs_old, ibox
  INTEGER :: nvecs_max, k, iatom, mcstep
  
  INTEGER :: ia
  
  REAL(DP) :: random_displacement, s(3), delta_volume, ln_pacc, delta_e, success_ratio
  REAL(DP) :: this_volume
  REAL(DP), DIMENSION(maxk) :: hx_old, hy_old, hz_old, Cn_old
  
  REAL(DP) :: dx, dy, dz, rijsq
  REAL(DP) :: E_vol_vdw, E_vol_qq, e_lrc
  REAL(DP) :: E_tot, E_vdw, E_lrc_old
  REAL(DP) :: pres_id
  REAL(DP) :: W_vol_vdw, W_vol_qq
  
  LOGICAL :: overlap, xz_change, accept_or_reject
  
  TYPE(Box_Class) :: box_list_old
  
  TYPE(Energy_Class) :: energy_old, virial_old
  
  REAL(DP), ALLOCATABLE :: pair_nrg_vdw_old(:,:), pair_nrg_qq_old(:,:)
  
  INTEGER, ALLOCATABLE :: my_species_id(:), my_locate_id(:), my_position_id(:)
  INTEGER :: position, my_box
  
  REAL(DP), ALLOCATABLE :: cos_mol_old(:,:), sin_mol_old(:,:)
  
  REAL(DP) :: rcut_vdw_old, rcut_coul_old, rcut3_old, rcut9_old, alpha_ewald_old
  REAL(DP) :: h_ewald_cut_old, rcut_vdwsq_old, rcut_coulsq_old, rcut_vdw3_old
  REAL(DP) :: rcut_vdw6_old, rcut_max_old
  
  REAL(DP) :: time1,time0
  
  CHARACTER(7) :: box_str, cutoff_str
  ! Framework related stuff
  REAL(DP) :: pore_width_old, ratio_width, area, half_pore_width_old
  
  
  ! Done with that section
  accept = .FALSE.
  
  ! Randomly choose a box for volume perturbation
  this_box = INT ( rranf() * nbr_boxes ) + 1
  
  ! increase the total number of trials for this box
  tot_trials(this_box) = tot_trials(this_box) + 1
  
  ! increment the number of volume moves for this box
  nvolumes(this_box) = nvolumes(this_box) + 1
  
  ! Perform a random walk in dv_max
  delta_volume = (1.0_DP - 2.0_DP * rranf()) * box_list(this_box)%dv_max
  this_volume = box_list(this_box)%volume + delta_volume
  
  ! Reject the move
  IF (this_volume < 0.0_DP) THEN
     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,5X,X,3X,X,I3,X,L8,X,9X,X,A9)') &
             i_mcstep, 'volume' , this_box, accept, 'neg_volume'
     END IF
  END IF
  
  ! store the old configuration of all atoms and COMs of molecules, also
  ! calculate the total number of molecules in the box to be used in
  ! the acceptance rule
  total_molecules = 0
  
  DO is = 1, nspecies
     
     DO im = 1, nmols(is,this_box)
        
        alive = locate(im,is,this_box)
  
        IF (molecule_list(alive,is)%live) THEN
  
           IF (molecule_list(alive,is)%which_box == this_box ) THEN
              
              total_molecules = total_molecules + 1
              
              ! store the unperturbed position of this molecule. Note that
              ! this call to the routine will save Cartesian coordinates
              ! as well as COMs of this molecule
              
              CALL Save_Old_Cartesian_Coordinates(alive,is)
              
           END IF
            
        END IF
  
     END DO
  
  END DO
  
  !  call cpu_time(time0)
  
  IF (l_pair_nrg) THEN
     
     ALLOCATE(pair_nrg_vdw_old(SUM(max_molecules),SUM(max_molecules)))
     ALLOCATE(pair_nrg_qq_old(SUM(max_molecules),SUM(max_molecules)))
     
     !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
     pair_nrg_vdw_old(:,:) = pair_nrg_vdw(:,:)
     pair_nrg_qq_old(:,:) = pair_nrg_qq(:,:)
     !!$OMP END PARALLEL WORKSHARE
     
  END IF
  
  IF (int_charge_sum_style(this_box) == charge_ewald) THEN
     
     ! store the cos_mol and sin_mol arrays
  
     ALLOCATE(cos_mol_old(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)
  
     IF (AllocateStatus /=0 ) THEN
        err_msg = ''
        err_msg(1) = 'cos_mol_old cannot be allocated'
        err_msg(2) = 'in Volume_Change'
        CALL Clean_Abort(err_msg,'Volume_Change')
     END IF
  
     ALLOCATE(sin_mol_old(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)
  
     IF (AllocateStatus /=0 ) THEN
        err_msg = ''
        err_msg(1) = 'sin_mol_old cannot be allocated'
        err_msg(2) = 'in Volume_Change'
        CALL Clean_Abort(err_msg,'Volume_Change')
     END IF
  
     
     cos_mol_old(:,:) = cos_mol(:,:)
     sin_mol_old(:,:) = sin_mol(:,:)
  
  END IF
  
  ! Store the box_list matrix
  box_list_old = box_list(this_box)
  
  ! Change the box size
  ! Assume box is cubic
  box_list(this_box)%length(1,1) = this_volume ** (1.0_DP/3.0_DP)
  box_list(this_box)%length(2,2) = box_list(this_box)%length(1,1)
  box_list(this_box)%length(3,3) = box_list(this_box)%length(2,2)
  
  ! compute the new components of box_list(this_box)
  CALL Compute_Cell_Dimensions(this_box)
  
  IF ( l_half_len_cutoff(this_box)) THEN
        
     ! store old cutoffs and other associated quantities
     rcut_vdw_old = rcut_vdw(this_box)
     rcut_coul_old = rcut_coul(this_box)
     rcut_vdwsq_old = rcut_vdwsq(this_box)
     rcut_coulsq_old = rcut_coulsq(this_box)
     
     rcut3_old = rcut3(this_box)
     rcut9_old = rcut9(this_box)
     rcut_vdw3_old = rcut_vdw3(this_box)
     rcut_vdw6_old = rcut_vdw6(this_box)
     
     rcut_max_old = rcut_max(this_box)
     
     IF( int_charge_sum_style(this_box) == charge_ewald ) THEN
        ! alpha_ewald_old = alpha_ewald(this_box)
        h_ewald_cut_old = h_ewald_cut(this_box)
     END IF
  
     ! change cutoffs and other associated quantities
     rcut_vdw(this_box) = 0.5_DP * box_list(this_box)%length(1,1)
     rcut_coul(this_box) = rcut_vdw(this_box)
     rcut_vdwsq(this_box) = rcut_vdw(this_box) * rcut_vdw(this_box)
     rcut_coulsq(this_box) = rcut_vdwsq(this_box)
     
     rcut_vdw3(this_box) = rcut_vdwsq(this_box) * rcut_vdw(this_box)
     rcut_vdw6(this_box) = rcut_vdw3(this_box) * rcut_vdw3(this_box)
     rcut3(this_box) = rcut_vdw3(this_box)
     rcut9(this_box) = rcut3(this_box) * rcut_vdw6(this_box)
  
     rcut_max(this_box) = rcut_vdw(this_box)
  
     IF ( int_charge_sum_style(this_box) == charge_ewald) THEN
        ! alpha_ewald(this_box) = ewald_p_sqrt(this_box) / rcut_coul(this_box)
        h_ewald_cut(this_box) = 2.0_DP * ewald_p(this_box) / rcut_coul(this_box)
     END IF
  
  ELSE
  
     ! Assumes box is cubic
     IF ( 0.5_DP * box_list(this_box)%length(1,1) < rcut_vdw(this_box) .OR. &
          0.5_DP * box_list(this_box)%length(1,1) < rcut_coul(this_box) .OR. &
          0.5_DP * box_list(this_box)%length(1,1) < roff_charmm(this_box) .OR. &
          0.5_DP * box_list(this_box)%length(1,1) < roff_switch(this_box) ) THEN
        err_msg = ''
        err_msg(1) = 'Cutoff is greater than the half box length'
        IF (nbr_boxes > 1) err_msg(2) = 'of box ' // TRIM(Int_To_String(this_box))
        CALL Clean_Abort(err_msg,'Volume_Change')
     END IF
  END IF
     
  delta_volume = box_list(this_box)%volume - box_list_old%volume
  
  ! we scale the coordinates of the COM of each of the molecules based on
  ! the old and new cell basis vectors. The idea is to keep the fractional
  ! coordinates of the COM the same before and after the move.
  
  DO is = 1, nspecies

     DO im = 1, nmols(is,this_box)
        
        alive = locate(im,is,this_box)
        
        ! obtain the new coordinates of the COM for this molecule
        
        ! first determine the fractional coordinate
        
        DO i = 1,3
           s(i) = box_list_old%length_inv(i,1) * molecule_list(alive,is)%xcom &
                + box_list_old%length_inv(i,2) * molecule_list(alive,is)%ycom &
                + box_list_old%length_inv(i,3) * molecule_list(alive,is)%zcom
        END DO
        
        ! now obtain the new positions of COMs
        molecule_list(alive,is)%xcom = box_list(this_box)%length(1,1) * s(1) &
                                     + box_list(this_box)%length(1,2) * s(2) &
                                     + box_list(this_box)%length(1,3) * s(3)
        
        molecule_list(alive,is)%ycom = box_list(this_box)%length(2,1) * s(1) &
                                     + box_list(this_box)%length(2,2) * s(2) &
                                     + box_list(this_box)%length(2,3) * s(3)
        
        molecule_list(alive,is)%zcom = box_list(this_box)%length(3,1) * s(1) &
                                     + box_list(this_box)%length(3,2) * s(2) &
                                     + box_list(this_box)%length(3,3) * s(3)
        
        ! Obtain the new positions of atoms in this molecule
        atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp + &
             molecule_list(alive,is)%xcom - molecule_list(alive,is)%xcom_old
        
        atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp + &
             molecule_list(alive,is)%ycom - molecule_list(alive,is)%ycom_old
        
        atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp + &
             molecule_list(alive,is)%zcom - molecule_list(alive,is)%zcom_old
        
     END DO
     
  END DO
  ! Energy change section,
  ! Store the old values of energy(this_box) and recompute the new components

 
  energy_old = energy(this_box)

  IF (int_charge_sum_style(this_box) == charge_ewald) THEN
     
     ! The calculation of Ewald reciprocal energy is tricky. The k vectors for this_box will change along with
     ! total number of k vectors for the box. It may so happen that the resultant k vectors may exceed the
     ! length of cos_sum and sin_sum arrays. We will handle this issue by first copying the cos_sum and sin_sum
     ! arrays in the old arrays. We will compute the new k vectors and determine if the shape of these vectors
     ! changes.
     
     ! store old terms for the Ewald reciprocal energy calculations

     nvecs_old = nvecs(this_box)
     nvecs_max = MAXVAL(nvecs)
     
     !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED) 
     cos_sum_old(:,:) = cos_sum(:,:)
     sin_sum_old(:,:) = sin_sum(:,:)
     hx_old(:) = hx(:,this_box)
     hy_old(:) = hy(:,this_box)
     hz_old(:) = hz(:,this_box)
     Cn_old(:) = Cn(:,this_box)
     !!$OMP END PARALLEL WORKSHARE


     ! Determine the new k vectors for this box. The call will change Cn, hx, hy and hz and hence will
     ! change cos_sum and sin_sum.
     CALL Ewald_Reciprocal_Lattice_Vector_Setup(this_box)

     ! cos_sum and sin_sum need to be re-allocated since the number of vectors has changed
     ! The operation destoys the cos_sum and sin_sum for other boxes but can be easily restored
     ! from cos_sum_old and sin_sum_old. Note that these terms for this_box will be calculated
     ! via the call to total system energy routine.
     ! Similar reasoning goes for cos_mol and sin_mol
     
     IF (ALLOCATED(cos_sum)) DEALLOCATE(cos_sum)
     IF (ALLOCATED(sin_sum)) DEALLOCATE(sin_sum)
     IF (ALLOCATED(cos_mol)) DEALLOCATE(cos_mol)
     IF (ALLOCATED(sin_mol)) DEALLOCATE(sin_mol)

     ALLOCATE(cos_sum(MAXVAL(nvecs),nbr_boxes), Stat = AllocateStatus)
     
     IF (Allocatestatus /= 0) THEN
        err_msg = ''
        err_msg(1) = 'Memory could not be allocated for cos_sum'
        err_msg(2) = 'allocation_cos_sin'
        CALL Clean_Abort(err_msg,'Volume_Change')
     END IF
     

     ALLOCATE(sin_sum(MAXVAL(nvecs),nbr_boxes), Stat = AllocateStatus)

     IF (Allocatestatus /= 0) THEN
        err_msg = ''
        err_msg(1) = 'Memory could not be allocated for sin_sum'
        err_msg(2) = 'allocation_cos_sin'
        CALL Clean_Abort(err_msg,'Volume_Change')
     END IF

     ALLOCATE(cos_mol(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)

     IF (Allocatestatus /= 0) THEN
        err_msg = ''
        err_msg(1) = 'Memory could not be allocated for cos_mol'
        err_msg(2) = 'allocation_cos_sin'
        CALL Clean_Abort(err_msg,'Volume_Change')
     END IF

     ALLOCATE(sin_mol(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)

     IF (Allocatestatus /= 0) THEN
        err_msg = ''
        err_msg(1) = 'Memory could not be allocated for sin_mol'
        err_msg(2) = 'allocation_cos_sin'
        CALL Clean_Abort(err_msg,'Volume_Change')
     END IF
     
  END IF

  CALL Compute_System_Total_Energy(this_box,.TRUE.,overlap)

  IF (overlap) THEN
     ! reject move
     CALL Reset_Coords
     
     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,5X,X,3X,X,I3,X,L8,X,9X,X,A9)') &
             i_mcstep, 'volume' , this_box, accept, 'overlap'
     END IF

  ELSE 
     
     ! change in the energy of the system 
     delta_e = energy(this_box)%total - energy_old%total
     
     ! based on the energy, calculate the acceptance ratio
     ln_pacc = beta(this_box) * delta_e &
             + beta(this_box) * pressure(this_box)%setpoint * delta_volume &
             - total_molecules * DLOG(box_list(this_box)%volume/box_list_old%volume)
     accept = accept_or_reject(ln_pacc)
     
     IF ( accept ) THEN

        nvol_success(this_box) = nvol_success(this_box) + 1
        ivol_success(this_box) = ivol_success(this_box) + 1
        ! energy, positions and box dimensions are already updated

        
        IF (int_charge_sum_style(this_box) == charge_ewald) THEN
           ! put the sin_sum and cos_sum for other boxes
           ! Note that hx,hy,hz and Cn are updated when the total energy routine is called above.
           ! there is no need to update these arrays.
           DO ibox = 1, nbr_boxes
              IF (ibox /= this_box ) THEN
                 !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
                 sin_sum(1:nvecs(ibox),ibox) = sin_sum_old(1:nvecs(ibox),ibox)
                 cos_sum(1:nvecs(ibox),ibox) = cos_sum_old(1:nvecs(ibox),ibox)
                 !!$OMP END PARALLEL WORKSHARE
              END IF
           END DO
           
           ! Now deallocate cos_sum_old and sin_sum_old so that they have the same dimensions
           ! as sin_sum and cos_sum 
           
           DEALLOCATE(cos_sum_old,sin_sum_old)
           ALLOCATE(cos_sum_old(SIZE(cos_sum,1),nbr_boxes))
           ALLOCATE(sin_sum_old(SIZE(sin_sum,1),nbr_boxes))
           DEALLOCATE(cos_sum_start,sin_sum_start)
           ALLOCATE(cos_sum_start(SIZE(cos_sum,1),nbr_boxes))
           ALLOCATE(sin_sum_start(SIZE(sin_sum,1),nbr_boxes))

           ! cos_mol
           
           !              CALL cpu_time(time0)
           
           DO ibox = 1, nbr_boxes
              ! skip the molecules in 'this_box'
              IF (ibox == this_box )CYCLE 
              DO is = 1, nspecies
                 DO im = 1, nmols(is,ibox)
                    alive = locate(im,is,ibox)
                    IF (.NOT. molecule_list(alive,is)%live) CYCLE
                    CALL Get_Position_Alive(alive,is,position)
                    
                    cos_mol(1:nvecs(ibox),position) = cos_mol_old(1:nvecs(ibox),position)
                    !                    cos_mol(nvecs(ibox)+1,nvecs(this_box)) = 0.0_DP
                    
                    sin_mol(1:nvecs(ibox),position) = sin_mol_old(1:nvecs(ibox),position)
                    !                   sin_mol(nvecs(ibox)+1,nvecs(this_box)) = 0.0_DP
                    
                 END DO
              END DO
           END DO
           
           DEALLOCATE(cos_mol_old,sin_mol_old)
        END IF

     ELSE
        
        ! Reject the move
        ! Reset the coordinates
        CALL Reset_Coords

     END IF
     
     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,5X,X,3X,X,I3,X,L8,X,9X,X,F9.3)') &
             i_mcstep, 'volume' , this_box, accept, ln_pacc
     END IF

  END IF

  ! Update the maximum displaement if there are nupdate_vol attempts
  
  IF (MOD(nvolumes(this_box),nvol_update) == 0 ) THEN
     IF (int_run_type == run_equil) THEN
        success_ratio = REAL(ivol_success(this_box),DP)/REAL(nvol_update,DP)
     ELSE
        success_ratio = REAL(nvol_success(this_box),DP)/REAL(nvolumes(this_box),DP)
     END IF

     WRITE(logunit,'(X,I9,X,A10,X,5X,X,3X,X,I3,X,F8.5)',ADVANCE='NO') i_mcstep, 'volume' , this_box, success_ratio

     IF (int_run_type == run_equil) THEN
        ! dv_max will be adjusted to achieve 0.5 acceptance using the formula
        !
        !     dv_max = 2 * success_ratio * dv_max
        !
        ! If the current success_ratio is too low, dv_max will be decreased.
        ! Otherwise, dv_max will be increased. 
      
        ! Assumes box is cubic, dv_max is scalar
        IF ( success_ratio < 0.0001 ) THEN
           box_list(this_box)%dv_max = 0.1_DP * box_list(this_box)%dv_max
        ELSE
           box_list(this_box)%dv_max = 2.0_DP * success_ratio * box_list(this_box)%dv_max
        END IF
        WRITE(logunit,'(X,F9.0)',ADVANCE='NO') box_list(this_box)%dv_max
        
        ivol_success(this_box) = 0
     END IF

     WRITE(logunit,*)

  END IF

  CONTAINS

    SUBROUTINE Reset_Coords
      
      DO is = 1, nspecies
         
         DO im = 1, nmols(is,this_box)
            
            alive = locate(im,is,this_box)

            IF (molecule_list(alive,is)%live) THEN
               
               CALL Revert_Old_Cartesian_Coordinates(alive,is)
               
            END IF
            
         END DO
         
      END DO
      
      ! Reset the box dimensions
      
      box_list(this_box) = box_list_old

      
      ! Reset the energy components
      
      energy(this_box) = energy_old

      IF (l_pair_nrg) THEN

         pair_nrg_vdw(:,:) = pair_nrg_vdw_old(:,:)
         pair_nrg_qq(:,:) = pair_nrg_qq_old(:,:)
         
         DEALLOCATE(pair_nrg_vdw_old,pair_nrg_qq_old)
         
      END IF
      
      IF (l_half_len_cutoff(this_box)) THEN
         
         rcut_vdw(this_box) = rcut_vdw_old
         rcut_coul(this_box) = rcut_coul_old
         rcut_vdwsq(this_box) = rcut_vdwsq_old
         rcut_coulsq(this_box) = rcut_coulsq_old
         
         rcut3(this_box) = rcut3_old
         rcut9(this_box) = rcut9_old
         rcut_vdw3(this_box) = rcut_vdw3_old
         rcut_vdw6(this_box) = rcut_vdw6_old
         
         rcut_max_old = rcut_max(this_box)
         
         IF( int_charge_sum_style(this_box) == charge_ewald ) THEN
            
           !           alpha_ewald(this_box) = alpha_ewald_old
           h_ewald_cut(this_box) = h_ewald_cut_old
           
        END IF
         

     END IF
         
     IF (int_charge_sum_style(this_box) == charge_ewald) THEN
            
         ! reset the terms related to Ewald reciprocal energy
         ! reset the total number of kvectors
         nvecs(this_box) = nvecs_old

         
         DEALLOCATE(cos_sum,sin_sum)
         DEALLOCATE(cos_mol,sin_mol)
         
         ALLOCATE(cos_sum(MAXVAL(nvecs),nbr_boxes),Stat = Allocatestatus)
         
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
         
         ALLOCATE(cos_mol(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)
         
         IF (Allocatestatus /= 0) THEN
            err_msg = ''
            err_msg(1) = 'Memory could not be in the volume rejection'
            CALL Clean_Abort(err_msg,'Volume_Change')
         END IF
         
         ALLOCATE(sin_mol(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)
         
         IF (Allocatestatus /= 0) THEN
            err_msg = ''
            err_msg(1) = 'Memory could not be allocated in the volume rejection'
            CALL Clean_Abort(err_msg,'Volume_Change')
         END IF
         
         !!$OMP PARALLEL WORKSHARE DEFAULT(SHARED) 
         cos_mol(1:SIZE(cos_mol_old,1),:) = cos_mol_old(:,:)
         sin_mol(1:SIZE(sin_mol_old,1),:) = sin_mol_old(:,:)
         cos_sum(:,:) = cos_sum_old(:,:)
         sin_sum(:,:) = sin_sum_old(:,:)
         hx(:,this_box) = hx_old(:)
         hy(:,this_box) = hy_old(:)
         hz(:,this_box) = hz_old(:)
         Cn(:,this_box) = Cn_old(:)
         !!$OMP END PARALLEL WORKSHARE
         
         ! here we make sure that cos_sum_old and sin_sum_old have the same dimensions
         ! as cos_sum and sin_sum
         DEALLOCATE(cos_mol_old,sin_mol_old)
         DEALLOCATE(cos_sum_old,sin_sum_old)
         ALLOCATE(cos_sum_old(SIZE(cos_sum,1),nbr_boxes),sin_sum_old(SIZE(sin_sum,1),nbr_boxes))
         DEALLOCATE(cos_sum_start,sin_sum_start)
         ALLOCATE(cos_sum_start(SIZE(cos_sum,1),nbr_boxes),sin_sum_start(SIZE(sin_sum,1),nbr_boxes))
         
      END IF
    END SUBROUTINE Reset_Coords

  END SUBROUTINE Volume_Change

