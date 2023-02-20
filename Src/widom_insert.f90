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

SUBROUTINE Widom_Insert(is,ibox,widom_sum,t_cpu, n_overlaps)

  !*****************************************************************************
  ! 
  ! PURPOSE: Perform all Widom insertions for species is and box ibox for the
  !          current step and return widom_sum.
  !
  ! Called by
  !
  !    Widom_Subdriver
  !
  ! 
  !*****************************************************************************

  USE Global_Variables
  USE Energy_Routines
  USE IO_Utilities
  USE Random_Generators
  USE Rotation_Routines
  USE Fragment_Growth
  USE File_Names
  !$ USE OMP_LIB

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  IMPLICIT NONE

  ! Arguments
  INTEGER, INTENT(IN) :: ibox ! insert test particles in ibox
  INTEGER, INTENT(IN) :: is ! species indices
  REAL(DP), INTENT(OUT) :: t_cpu

  ! Local declarations
  INTEGER :: i, i_type               ! atom indices
  INTEGER :: im                      ! molecule INDEX
  INTEGER :: frag_order(nfragments(is))
  INTEGER :: subinterval, i_interval

  INTEGER (KIND=INT64) :: i_widom
  INTEGER (KIND=INT64) :: insertions_in_step, n_overlaps

  REAL(DP) :: dx, dy, dz
  REAL(DP) :: dE, dE_intra, dE_inter, dE_frag
  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq
  REAL(DP) :: E_inter, E_periodic_qq
  REAL(DP) :: E_reciprocal, E_self, E_lrc
  REAL(DP) :: E_ring_frag
  REAL(DP) :: ln_pacc, ln_pseq, ln_pbias, this_lambda

  REAL(DP) :: thread_Eij_factor, Eij_factor_gcopy
  REAL(DP), DIMENSION(0:Eij_ind_ubound) :: frame_w_max, w_max_gcopy
  REAL(DP), DIMENSION(0:Eij_ind_ubound) :: frame_Eij_w_sum, Eij_w_sum_gcopy
  INTEGER, DIMENSION(0:Eij_ind_ubound) :: Eij_freq, Eij_freq_gcopy


  INTEGER :: changefactor, thread_changefactor, Eij_ind

  REAL(DP), DIMENSION(rsqmin_res_d,solvent_maxind_d,natoms(is)), TARGET :: frame_rsqmin_atompair_w_sum
  REAL(DP), DIMENSION(rsqmin_res_d,solvent_maxind_d,natoms(is)), TARGET :: frame_rsqmin_atompair_w_max
  INTEGER(KIND=INT64), DIMENSION(rsqmin_res_d,solvent_maxind_d,natoms(is)), TARGET :: frame_rsqmin_atompair_freq
  REAL(DP), DIMENSION(:,:,:), POINTER :: rsqmin_atompair_w_max_ptr, rsqmin_atompair_w_sum_ptr
  INTEGER(KIND=INT64), DIMENSION(:,:,:), POINTER :: rsqmin_atompair_freq_ptr
  INTEGER :: bsolute, rsq_ind, ia, ti_solvent

  REAL(DP) :: widom_prefactor, widom_var_exp, widom_sum
  REAL(DP) :: E_recip_in, lrc_diff, E_inter_constant
  REAL(DP) :: subinterval_sums(MAX(n_widom_subgroups(is,ibox),1))
  REAL(DP) :: t_cpu_e, t_cpu_s
  INTEGER :: n_subintervals
!widom_timing  REAL(DP) :: noncbmc_time_e, noncbmc_time_s, noncbmc_time
  LOGICAL :: write_wprp2
  LOGICAL :: omp_flag


  LOGICAL :: inter_overlap, cbmc_overlap, intra_overlap

  !!! 4-d allocatable shared target arrays with the last axis for thread number are used 
  !       in conjunction with private pointers and sum or maxval is used along the last axis
  !       of the target array.  This functionally imitates the behavior of having 3-d variables
  !       in OMP REDUCTION clauses, except it allocates
  !       the arrays to the heap, whereas REDUCTION variables are kept in the stack. 
  !       This is now optional, and either the stack or the heap may be used, with the stack
  !       being the default.  If you want to use the stack but get segmentation faults 
  !       due to stack buffer overflow, set environment variable OMP_STACKSIZE to increase
  !       the size of the stack for threads other than the main thread.
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: frame_rsqmin_atompair_w_sum_tgt, frame_rsqmin_atompair_w_max_tgt
  INTEGER(KIND=INT64), DIMENSION(:,:,:,:), ALLOCATABLE, TARGET :: frame_rsqmin_atompair_freq_tgt
  REAL(DP), DIMENSION(:,:,:), POINTER :: frame_rsqmin_atompair_w_sum_ptr, frame_rsqmin_atompair_w_max_ptr
  INTEGER(KIND=INT64), DIMENSION(:,:,:), POINTER :: frame_rsqmin_atompair_freq_ptr
  !!!

  INTEGER :: nbr_threads, ithread
  nbr_threads = 1
  omp_flag = .FALSE.
  !$ omp_flag = .TRUE.
  !$ nbr_threads = omp_get_max_threads()

!widom_timing  noncbmc_time = 0.0_DP


  t_cpu = 0.0_DP

  IF (est_atompair_rminsq) THEN
          bsolute = species_list(is)%wsolute_base
          rsqmin_atompair_w_max_ptr => rsqmin_atompair_w_max(:,:, &
                  bsolute+1:bsolute+natoms(is),ibox)
          rsqmin_atompair_w_sum_ptr => rsqmin_atompair_w_sum(:,:, &
                  bsolute+1:bsolute+natoms(is),ibox)
          rsqmin_atompair_freq_ptr => rsqmin_atompair_freq(:,:, &
                  bsolute+1:bsolute+natoms(is),ibox)
          IF (l_heap) THEN
                  IF (ALLOCATED(frame_rsqmin_atompair_w_sum_tgt)) THEN
                          DEALLOCATE(frame_rsqmin_atompair_w_sum_tgt)
                  END IF
                  IF (ALLOCATED(frame_rsqmin_atompair_w_max_tgt)) THEN
                          DEALLOCATE(frame_rsqmin_atompair_w_max_tgt)
                  END IF
                  IF (ALLOCATED(frame_rsqmin_atompair_freq_tgt)) THEN
                          DEALLOCATE(frame_rsqmin_atompair_freq_tgt)
                  END IF
                  ALLOCATE(frame_rsqmin_atompair_w_sum_tgt( &
                          rsqmin_res, &
                          solvent_maxind, &
                          natoms(is), &
                          nbr_threads))
                  ALLOCATE(frame_rsqmin_atompair_w_max_tgt( &
                          rsqmin_res, &
                          solvent_maxind, &
                          natoms(is), &
                          nbr_threads))
                  ALLOCATE(frame_rsqmin_atompair_freq_tgt( &
                          rsqmin_res, &
                          solvent_maxind, &
                          natoms(is), &
                          nbr_threads))
                  !$OMP PARALLEL WORKSHARE
                  frame_rsqmin_atompair_w_sum_tgt = 0.0_DP
                  frame_rsqmin_atompair_w_max_tgt = 0.0_DP
                  frame_rsqmin_atompair_freq_tgt = 0_INT64
                  !$OMP END PARALLEL WORKSHARE
          END IF
  END IF



  IF (n_widom_subgroups(is,ibox) > 0) THEN
          n_subintervals = n_widom_subgroups(is,ibox)
          write_wprp2 = .TRUE.
  ELSE
          n_subintervals = 1
          write_wprp2 = .FALSE.
  END IF

  changefactor = 1
  thread_Eij_factor = 0.0_DP
  w_max_gcopy(:) = w_max(:,is,ibox)
  Eij_w_sum_gcopy(:) = Eij_w_sum(:,is,ibox)
  Eij_freq_gcopy(:) = Eij_freq_total(:,is,ibox)
  Eij_factor_gcopy = Eij_factor(is,ibox)

  this_lambda = 1.0_DP
  widom_sum = 0.0_DP
  n_overlaps = 0_INT64
  subinterval_sums = 0.0_DP

  lrc_diff = 0.0_DP
  E_self = 0.0_DP
  E_recip_in = 0.0_DP

  del_flag = .FALSE.
  get_fragorder = .TRUE.
  
  nmols(is,ibox) = nmols(is,ibox)+1
  im = nmols(is,ibox)
  locate(im,is,ibox) = locate(nmols(is,0),is,0)
  locate(nmols(is,0),is,0) = 0

  !  * Set properties of the to-be-inserted molecule
  widom_species = is
  widom_locate = locate(im,is,ibox)
  molecule_list(widom_locate,is)%which_box = ibox
  molecule_list(widom_locate,is)%frac = this_lambda
  molecule_list(widom_locate,is)%molecule_type = int_normal
  
  widom_prefactor = box_list(ibox)%volume&
                  / (REAL(nmols(is,ibox),DP)*((species_list(is)%de_broglie(ibox))**3))
  insertions_in_step = species_list(is)%insertions_in_step(ibox)
  
  ! Long-range energy correction

  IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN

     ! increase number of integer beads
     nbeads_in = nint_beads(:,ibox)

     DO i = 1, natoms(is)
        i_type = nonbond_list(i,is)%atom_type_number
        nint_beads(i_type,ibox) = nint_beads(i_type,ibox) + 1
     END DO

     CALL Compute_LR_correction(ibox,E_lrc)
     nint_beads(:,ibox) = nbeads_in
     lrc_diff = E_lrc - energy(ibox)%lrc

  END IF
  
  IF (int_charge_style(ibox) == charge_coul .AND. has_charge(is)) THEN
          CALL Compute_Molecule_Self_Energy(widom_locate,is,ibox,E_self)
          E_recip_in = energy(ibox)%reciprocal
  END IF
  E_inter_constant = lrc_diff + E_self - E_recip_in

  widom_active = .TRUE.


  subinterval = insertions_in_step/n_subintervals

  ! add $ in next line and include in omp parallel creation below when using full widom timing
  !OMP PRIVATE(noncbmc_time_e,noncbmc_time_s) &
  ! also include noncbmc_time in omp reduction

  frame_Eij_w_sum = 0.0_DP
  Eij_freq = 0
  frame_w_max = 0.0_DP
  frame_rsqmin_atompair_w_sum = 0.0_DP
  frame_rsqmin_atompair_w_max = 0.0_DP
  frame_rsqmin_atompair_freq = 0_INT64

  !$OMP PARALLEL DEFAULT(SHARED) &
  !$OMP PRIVATE(ln_pseq, ln_pbias, E_ring_frag, inter_overlap, cbmc_overlap, intra_overlap, i_interval) &
  !$OMP PRIVATE(widom_var_exp, E_periodic_qq, E_intra_qq, E_intra_vdw, E_inter, dE_frag, dE) &
  !$OMP PRIVATE(E_bond, E_angle, E_dihedral, E_improper, dE_intra, dE_inter, E_reciprocal, frag_order) &
  !$OMP PRIVATE(t_cpu_s, t_cpu_e, thread_changefactor, Eij_ind, rsq_ind, ia, ti_solvent) &
  !$OMP PRIVATE(frame_rsqmin_atompair_w_sum_ptr,frame_rsqmin_atompair_w_max_ptr,frame_rsqmin_atompair_freq_ptr) &
  !$OMP REDUCTION(+:widom_sum,n_overlaps,subinterval_sums,t_cpu,Eij_freq,frame_Eij_w_sum) &
  !$OMP REDUCTION(+:frame_rsqmin_atompair_w_sum,frame_rsqmin_atompair_freq) &
  !$OMP REDUCTION(MAX:frame_w_max,thread_Eij_factor,frame_rsqmin_atompair_w_max)
  frame_rsqmin_atompair_w_sum = 0.0_DP
  frame_rsqmin_atompair_w_max = 0.0_DP
  frame_rsqmin_atompair_freq = 0_INT64

  ithread = 1
  !$ ithread = omp_get_thread_num() + 1
  IF (est_atompair_rminsq) THEN
          IF (l_heap) THEN
                  frame_rsqmin_atompair_w_sum_ptr => frame_rsqmin_atompair_w_sum_tgt(:,:,:,ithread)
                  frame_rsqmin_atompair_w_max_ptr => frame_rsqmin_atompair_w_max_tgt(:,:,:,ithread)
                  frame_rsqmin_atompair_freq_ptr => frame_rsqmin_atompair_freq_tgt(:,:,:,ithread)
          ELSE
                  frame_rsqmin_atompair_w_sum_ptr => frame_rsqmin_atompair_w_sum
                  frame_rsqmin_atompair_w_max_ptr => frame_rsqmin_atompair_w_max
                  frame_rsqmin_atompair_freq_ptr => frame_rsqmin_atompair_freq
          END IF
  END IF

  IF (ALLOCATED(widom_atoms)) THEN
          DEALLOCATE(widom_atoms, STAT=DeallocateStatus)
          IF (DeallocateStatus .NE. 0) THEN
                  WRITE(*,*) "Error deallocating widom_atoms"
                  STOP
          END IF
  END IF
  IF (ALLOCATED(swi_atompair_rsqmin)) THEN
          DEALLOCATE(swi_atompair_rsqmin, STAT=DeallocateStatus)
          IF (DeallocateStatus .NE. 0) THEN
                  WRITE(*,*) "Error deallocating swi_atompair_rsqmin"
                  STOP
          END IF
  END IF
  ALLOCATE(widom_atoms(natoms(is)), STAT=AllocateStatus)
  IF (AllocateStatus .NE. 0) THEN
          WRITE(*,*) "Error allocating widom_atoms"
          STOP
  END IF
  ALLOCATE(swi_atompair_rsqmin(solvent_maxind,natoms(is)),STAT=allocatestatus)
  IF (AllocateStatus .NE. 0) THEN
          WRITE(*,*) "Error allocating swi_atompair_rsqmin"
          STOP
  END IF
  widom_molecule = molecule_list(widom_locate,is)
  widom_atoms = atom_list(1:natoms(is),widom_locate,is)
  widom_sum = 0.0_DP
  n_overlaps = 0_INT64
  subinterval_sums = 0.0_DP
  t_cpu = 0.0_DP
  frame_w_max = 0.0_DP
  frame_Eij_w_sum = 0.0_DP
  Eij_freq = 0
  thread_Eij_factor = Eij_factor_gcopy
  thread_changefactor = changefactor

  !$OMP DO SCHEDULE(DYNAMIC)
  DO i_widom = 1, insertions_in_step
          IF (.NOT. omp_flag) CALL CPU_TIME(t_cpu_s)
          !$ t_cpu_s = omp_get_wtime()
          ! Initialize variables
          swi_atompair_rsqmin = 10000.0_DP
          ln_pseq = 0.0_DP
          ln_pbias = 0.0_DP
          E_ring_frag = 0.0_DP
          inter_overlap = .FALSE.
          cbmc_overlap = .FALSE.
          intra_overlap = .FALSE.
          widom_var_exp = 0.0_DP
          ! Now that an insertion will be attempted, we need to do some bookkeeping:

          !  * Increment the counters to track number of widom insertions

         
          !*****************************************************************************
          ! Choose a position, orientation and conformation for the 
          !         to-be-inserted molecule
          !*****************************************************************************
          !
          ! Build_Molecule places the first fragment, then calls Fragment_Placement
          ! to place the additional fragments 
          CALL Build_Molecule(widom_locate,is,ibox,frag_order,this_lambda, &
                  ln_pseq,ln_pbias,E_ring_frag,cbmc_overlap)

          ! Turn the molecule on
          widom_molecule%live = .TRUE.
          widom_atoms%exist = .TRUE.

          ! So far ln_pbias includes 
          !   * the probability of choosing the insertion 
          !     point from the collection of trial coordinates 
          !   * the probability of choosing each dihedral from the collection of trial dihedrals. 
          ! Now add
          !   * the probability of the fragment sequence, ln_pseq
          !   * the number of trial coordinates, kappa_ins
          !   * the number of trial dihedrals, kappa_dih, for each dihedral.

          ln_pbias = ln_pbias + ln_pseq
          ln_pbias = ln_pbias + DLOG(REAL(kappa_ins,DP))

          IF (kappa_rot /= 0 ) THEN
             ln_pbias = ln_pbias + DLOG(REAL(kappa_rot,DP))
          END IF

          IF (kappa_dih /= 0 ) THEN
             ln_pbias = ln_pbias + REAL(nfragments(is)-1,DP) * DLOG(REAL(kappa_dih,DP))
          END IF

          IF (.NOT. cbmc_overlap) THEN
!widom_timing            IF (.NOT. omp_flag) CALL cpu_time(noncbmc_time_s)
!widom_timing            !$ noncbmc_time_s = omp_get_wtime()

            ! Molecule COM may be outside the box boundary if grown via CBMC, so wrap
            ! the molecule coordinates back in the box (if needed)
            IF (nfragments(is) > 1) CALL Fold_Molecule(widom_locate,is,ibox)

            ! Recompute the COM in case the molecule was wrapped
            !CALL Get_COM(widom_locate,is)

            ! Compute the distance of the atom farthest from COM
            !CALL Compute_Max_COM_Distance(widom_locate,is)

            ! Calculate the potential energy interaction between the inserted molecule
            ! and the rest of the system
            Eij_max = 0.0_DP
            CALL Compute_Molecule_Nonbond_Inter_Energy_Widom(widom_locate,is, &
                    E_inter,inter_overlap)

            ! Calculate the nonbonded energy interaction within the inserted molecule
            IF (.NOT. inter_overlap) THEN
                    CALL Compute_Molecule_Nonbond_Intra_Energy(widom_locate,is, &
                            E_intra_vdw,E_intra_qq,E_periodic_qq,intra_overlap)
                    E_inter = E_inter + E_periodic_qq
            END IF

!widom_timing            IF (.NOT. omp_flag) CALL cpu_time(noncbmc_time_e)
!widom_timing            !$ noncbmc_time_e = omp_get_wtime()
!widom_timing            noncbmc_time = noncbmc_time + (noncbmc_time_e - noncbmc_time_s)
         
          END IF

          ! Leave widom_sum unchanged if there is any core overlap
          IF (.NOT. (cbmc_overlap .OR. inter_overlap .OR. intra_overlap)) THEN

                  ! There are no overlaps, so we can calculate the change in potential energy.
                  !
                  ! Already have the change in nonbonded energies
                  dE_inter = E_inter + E_inter_constant 
                  dE_intra = E_intra_vdw + E_intra_qq

                  ! Bonded intramolecular energies
                  ! If the molecule was grown via CBMC, we already have the intramolecular 
                  ! bond energies? Otherwise we need to compute them.
                  CALL Compute_Molecule_Bond_Energy(widom_locate,is,E_bond)
                  CALL Compute_Molecule_Angle_Energy(widom_locate,is,E_angle)
                  CALL Compute_Molecule_Dihedral_Energy(widom_locate,is,E_dihedral)
                  CALL Compute_Molecule_Improper_Energy(widom_locate,is,E_improper)

                  dE_intra = dE_intra + E_bond + E_angle + E_dihedral + E_improper

                  ! Ewald energies
                  IF (int_charge_style(ibox) == charge_coul) THEN
                        IF ( (int_charge_sum_style(ibox) == charge_ewald) .AND. &
                             has_charge(is) ) THEN
                       
                            CALL Update_System_Ewald_Reciprocal_Energy_Widom(widom_locate, &
                                   is,ibox,E_reciprocal)

                            dE_inter = dE_inter + E_reciprocal
                        END IF
                  END IF

                  ! moved to before the loop
                  !widom_prefactor = box_list(ibox)%volume&
                  !        / (REAL(nmols(is,ibox),DP)*((species_list(is)%de_broglie(ibox))**3))

                  ! change in energy, less energy used to bias fragment selection
                  dE = dE_intra + dE_inter
                  dE_frag = E_angle + E_ring_frag

                  ! mu' = -(1/beta)*ln(<widom_var>)
                  widom_var_exp = DEXP(-beta(ibox) * (dE - dE_frag) - ln_pbias)
                  ! sum of all widom_var for this step; output argument
                  widom_sum = widom_sum + widom_var_exp


                  IF (est_atompair_rminsq) THEN
                          DO ia = 1, natoms(is)
                                DO ti_solvent = 1, solvent_maxind
                                        IF (swi_atompair_rsqmin(ti_solvent,ia) >= maxrminsq) CYCLE
                                        rsq_ind = INT((swi_atompair_rsqmin(ti_solvent,ia)-rsqmin_shifter) / rsqmin_step)
                                        IF (widom_var_exp > frame_rsqmin_atompair_w_max_ptr(rsq_ind,ti_solvent,ia)) THEN
                                                frame_rsqmin_atompair_w_max_ptr(rsq_ind,ti_solvent,ia) = widom_var_exp
                                        END IF
                                        frame_rsqmin_atompair_w_sum_ptr(rsq_ind,ti_solvent,ia) = &
                                                frame_rsqmin_atompair_w_sum_ptr(rsq_ind,ti_solvent,ia) + &
                                                widom_var_exp
                                        frame_rsqmin_atompair_freq_ptr(rsq_ind,ti_solvent,ia) = &
                                                frame_rsqmin_atompair_freq_ptr(rsq_ind,ti_solvent,ia) + 1
                                END DO
                          END DO
                  END IF

                  Eij_ind = INT(Eij_max * thread_Eij_factor)
                  DO WHILE (Eij_ind > Eij_ind_ubound)
                        CALL coarsen_w_max(frame_w_max,frame_Eij_w_sum,Eij_freq,thread_Eij_factor,2)
                        thread_changefactor = thread_changefactor * 2
                        Eij_ind = INT(Eij_max * thread_Eij_factor)
                  END DO
                  frame_w_max(Eij_ind) = MAX(frame_w_max(Eij_ind), widom_var_exp)
                  frame_Eij_w_sum(Eij_ind) = frame_Eij_w_sum(Eij_ind) + widom_var_exp
                  Eij_freq(Eij_ind) = Eij_freq(Eij_ind) + 1
          ELSE
                  n_overlaps = n_overlaps + 1_INT64
          END IF
          IF (write_wprp2) THEN
                i_interval = (i_widom - 1_INT64)/subinterval + 1_INT64
                IF (i_interval <= n_subintervals) subinterval_sums(i_interval) = subinterval_sums(i_interval) + widom_var_exp
          END IF
          IF (.NOT. omp_flag) CALL CPU_TIME(t_cpu_e)
          !$ t_cpu_e = omp_get_wtime()
          t_cpu = t_cpu + t_cpu_e - t_cpu_s
  END DO
  !$OMP END DO
  !$OMP CRITICAL
  changefactor = MAX(changefactor,thread_changefactor)
  !$OMP END CRITICAL
  !$OMP BARRIER
  IF (thread_changefactor<changefactor) CALL coarsen_w_max(frame_w_max,frame_Eij_w_sum,Eij_freq,thread_Eij_factor,changefactor/thread_changefactor)
  !$OMP END PARALLEL
  IF (changefactor>1) CALL coarsen_w_max(w_max_gcopy,Eij_w_sum_gcopy,Eij_freq_gcopy,Eij_factor_gcopy,changefactor)
  Eij_factor(is,ibox) = Eij_factor_gcopy
  !$OMP PARALLEL
  !$OMP SECTIONS
  !$OMP SECTION
  w_max(:,is,ibox) = MAX(w_max_gcopy, widom_prefactor*frame_w_max)
  !$OMP SECTION
  Eij_w_sum(:,is,ibox) = Eij_w_sum_gcopy + widom_prefactor*frame_Eij_w_sum
  !$OMP SECTION
  Eij_freq_total(:,is,ibox) = Eij_freq_gcopy + Eij_freq
  !$OMP END SECTIONS NOWAIT

  IF (est_atompair_rminsq) THEN
          IF (l_heap) THEN
                  !$OMP SECTIONS
                  !$OMP SECTION
                  rsqmin_atompair_freq_ptr = rsqmin_atompair_freq_ptr + SUM(frame_rsqmin_atompair_freq_tgt,4)
                  !$OMP SECTION
                  rsqmin_atompair_w_max_ptr = MAX(rsqmin_atompair_w_max_ptr, widom_prefactor*MAXVAL(frame_rsqmin_atompair_w_max_tgt,4))
                  !$OMP SECTION
                  rsqmin_atompair_w_sum_ptr = rsqmin_atompair_w_sum_ptr + widom_prefactor*SUM(frame_rsqmin_atompair_w_sum_tgt,4)
                  !$OMP END SECTIONS
          ELSE
                  !$OMP SECTIONS
                  !$OMP SECTION
                  rsqmin_atompair_freq_ptr = rsqmin_atompair_freq_ptr + frame_rsqmin_atompair_freq
                  !$OMP SECTION
                  rsqmin_atompair_w_sum_ptr = rsqmin_atompair_w_sum_ptr + widom_prefactor*frame_rsqmin_atompair_w_sum
                  !$OMP SECTION
                  rsqmin_atompair_w_max_ptr = MAX(rsqmin_atompair_w_max_ptr,widom_prefactor*frame_rsqmin_atompair_w_max)
                  !$OMP END SECTIONS
          END IF
  END IF
  !$OMP END PARALLEL

  widom_active = .FALSE.
  widom_sum = widom_sum * widom_prefactor
  overlap_counter(is,ibox) = overlap_counter(is,ibox) + n_overlaps
  IF (write_wprp2) THEN
          subinterval_sums = subinterval_sums * widom_prefactor / subinterval
          IF (first_open_wprop2(is,ibox)) THEN
                  OPEN(unit=wprop2_file_unit(is,ibox),file=wprop2_filenames(is,ibox))
                  first_open_wprop2(is,ibox) = .FALSE.
          END IF
          DO i = 1, n_subintervals
                IF (subinterval_sums(i) < 1.0e-98_DP) subinterval_sums(i) = 0.0_DP
                WRITE(wprop2_file_unit(is,ibox), "(E30.22)", ADVANCE="NO") subinterval_sums(i)
          END DO
          WRITE(wprop2_file_unit(is,ibox),*)
  END IF


  ! remove test molecule
  nmols(is,ibox) = nmols(is,ibox)-1
  locate(im,is,ibox) = 0
  molecule_list(widom_locate,is)%live = .FALSE.
  atom_list(:,widom_locate,is)%exist = .FALSE.
  molecule_list(widom_locate,is)%molecule_type = int_none

  ! move locate to the list of unused locates
  locate(nmols(is,0),is,0) = widom_locate
  widom_locate = 0
  widom_species = 0


  ntrials(is,ibox)%widom = ntrials(is,ibox)%widom + insertions_in_step

  CONTAINS
          SUBROUTINE coarsen_w_max(wmax,wsum,Efreq,Efactor,cfactor)
                  REAL(DP), DIMENSION(0:Eij_ind_ubound), INTENT(INOUT) :: wmax,wsum
                  INTEGER, DIMENSION(0:Eij_ind_ubound), INTENT(INOUT) :: Efreq
                  REAL(DP), INTENT(INOUT) :: Efactor
                  INTEGER, INTENT(IN) :: cfactor
                  INTEGER :: gwidth, i1, i2, ic
                  Efactor = Efactor / REAL(cfactor,DP)
                  gwidth = cfactor - 1
                  ic = 0
                  DO i2 = gwidth, Eij_ind_ubound, cfactor
                        i1 = i2-gwidth
                        wmax(ic) = MAXVAL(wmax(i1:i2))
                        Efreq(ic) = SUM(Efreq(i1:i2))
                        wsum(ic) = SUM(wsum(i1:i2))
                        ic = ic + 1
                  END DO
                  wmax(ic:) = 0.0_DP
                  wsum(ic:) = 0.0_DP
                  Efreq(ic:) = 0
          END SUBROUTINE coarsen_w_max


!widom_timing  WRITE(*,*) noncbmc_time

END SUBROUTINE Widom_Insert
!*******************************************************************************
