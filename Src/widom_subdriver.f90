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

SUBROUTINE Widom_Subdriver

  !*****************************************************************************
  ! 
  ! PURPOSE: This subroutine checks which species and box combinations should
  !          have Widom insertions in the current step, calls Widom_Insert to 
  !          execute Widom insertions for those species and box combinations,
  !          and calls Write_Widom_Properties to write the results to Widom 
  !          property files.
  !
  ! Called by
  !
  !   nvtmc_driver
  !   nptmc_driver
  !   gcmc_driver
  !   gemc_driver 
  !
  !*****************************************************************************

  USE Global_Variables
  USE Sector_Routines
  !$ USE OMP_LIB

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  IMPLICIT NONE


  ! Local declarations
  INTEGER :: is, ibox
  INTEGER(KIND=INT64) :: n_overlaps
!widom_timing  INTEGER(KIND=INT64) :: num_cell_list_overlap, num_not_cell_list_overlap, num_nrg_overlap
  REAL(DP) :: widom_sum, widom_avg!widom_timing, setup_time_s, setup_time_e, setup_time
  REAL(DP) :: t_wc_s, t_wc_e, t_cpu
!widom_timing  REAL(DP) :: r_cell_list_time, r_normal_overlap_time, r_non_overlap_time, r_nrg_overlap_time
  LOGICAL :: need_init, omp_flag
  need_init = l_sectors
  omp_flag = .FALSE.
  !$ omp_flag = .TRUE.
  ! Loop over all species
  DO is = 1, nspecies
        ! Loop over all boxes
        DO ibox = 1, nbr_boxes
                ! move on to next box if this species is not used as a test particle in this box
                IF (.NOT. species_list(is)%test_particle(ibox)) CYCLE
                ! move on to next box if the step number isn't divisible by the widom_interval for this species.
                ! these are separate IF statements because we don't want to divide by zero
                IF (MOD(i_mcstep,species_list(is)%widom_interval(ibox)) .NE. 0) CYCLE
                IF (need_init) THEN
!widom_timing                        IF (.NOT. omp_flag) CALL cpu_time(setup_time_s)
!widom_timing                        !$ setup_time_s = omp_get_wtime()
                        CALL Sector_Setup
!widom_timing                        IF (.NOT. omp_flag) CALL cpu_time(setup_time_e)
!widom_timing                        !$ setup_time_e = omp_get_wtime()
!widom_timing                        setup_time = setup_time_e - setup_time_s
!widom_timing                        WRITE(*,*) setup_time
                        need_init = .FALSE.
                END IF
!widom_timing                !$OMP PARALLEL
!widom_timing                n_clo = 0_INT64
!widom_timing                n_not_clo = 0_INT64
!widom_timing                n_nrg_overlap = 0_INT64
!widom_timing                cell_list_time = 0.0_DP
!widom_timing                normal_overlap_time = 0.0_DP
!widom_timing                non_overlap_time = 0.0_DP
!widom_timing                nrg_overlap_time = 0.0_DP
!widom_timing                !$OMP END PARALLEL
                IF (.NOT. omp_flag) CALL CPU_TIME(t_wc_s)
                !$ t_wc_s = omp_get_wtime()
                CALL Widom_Insert(is,ibox,widom_sum,t_cpu,n_overlaps)
                widom_cpu_time(is,ibox) = widom_cpu_time(is,ibox) + t_cpu
                species_list(is)%widom_sum(ibox) = species_list(is)%widom_sum(ibox) + widom_sum
                widom_avg = widom_sum / species_list(is)%insertions_in_step(ibox)
                CALL Write_Widom_Properties(is,ibox,widom_avg,t_cpu,n_overlaps)
                IF (.NOT. omp_flag) CALL CPU_TIME(t_wc_e)
                !$ t_wc_e = omp_get_wtime()
                widom_wc_time(is,ibox) = (t_wc_e - t_wc_s) + widom_wc_time(is,ibox)
!widom_timing                num_cell_list_overlap = 0_INT64
!widom_timing                num_not_cell_list_overlap = 0_INT64
!widom_timing                num_nrg_overlap = 0_INT64
!widom_timing                r_cell_list_time = 0.0_DP
!widom_timing                r_normal_overlap_time = 0.0_DP
!widom_timing                r_non_overlap_time = 0.0_DP
!widom_timing                r_nrg_overlap_time = 0.0_DP
!widom_timing                !$OMP PARALLEL DEFAULT(SHARED) REDUCTION(+:num_cell_list_overlap,num_not_cell_list_overlap) &
!widom_timing                !$OMP REDUCTION(+: r_cell_list_time, r_normal_overlap_time, r_non_overlap_time) &
!widom_timing                !$OMP REDUCTION(+: num_nrg_overlap, r_nrg_overlap_time)
!widom_timing                num_cell_list_overlap = n_clo
!widom_timing                num_not_cell_list_overlap = n_not_clo
!widom_timing                num_nrg_overlap = n_nrg_overlap
!widom_timing                r_cell_list_time = cell_list_time
!widom_timing                r_normal_overlap_time = normal_overlap_time
!widom_timing                r_non_overlap_time = non_overlap_time
!widom_timing                r_nrg_overlap_time = nrg_overlap_time
!widom_timing                !$OMP END PARALLEL
!widom_timing                WRITE(*,*) num_cell_list_overlap
!widom_timing                WRITE(*,*) num_not_cell_list_overlap
!widom_timing                WRITE(*,*) num_nrg_overlap
!widom_timing                WRITE(*,*) r_cell_list_time
!widom_timing                WRITE(*,*) r_normal_overlap_time
!widom_timing                WRITE(*,*) r_nrg_overlap_time
!widom_timing                WRITE(*,*) r_non_overlap_time
        END DO
  END DO
END SUBROUTINE Widom_Subdriver
!*******************************************************************************
