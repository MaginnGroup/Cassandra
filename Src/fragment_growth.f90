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

MODULE Fragment_Growth

  !*****************************************************************************
  ! The module performs the tasks related to insertion, deletion, cut and 
  ! regrowth of molecules based on fragment sampling. It also contains a routine
  ! to obtain order of fragment placement.
  !
  ! CONTAINS
  !      
  ! Build_Molecule
  ! Build_Rigid_Fragment
  ! Cut_Regrow
  ! Fragment_Order
  ! Fragment_Placement
  ! Get_Aligner_Hanger
  ! Single_Fragment_Regrowth
  ! Get_Common_Fragment_Atoms
  !
  ! Used by
  !
  !   chempot
  !   cut_n_grow
  !   deletion
  !   gemc_particle_transfer
  !   make_config
  !   insertion
  !   main
  !   ring_fragment_driver
  ! 
  ! Revision history
  !
  !   12/10/13 : Beta Release
  !   Version 1.1
  !     05/15/15 Modified Fragment_Order to choose the next fragment with 
  !              uniform probability
  !     
  !*****************************************************************************

  USE Global_Variables
  USE Energy_Routines
  USE Random_Generators
  USE Read_Write_Checkpoint
  USE Sector_Routines
  !$ USE OMP_LIB

  IMPLICIT NONE

CONTAINS

!*******************************************************************************
SUBROUTINE Build_Molecule(this_im,is,this_box,frag_order,this_lambda, &
              ln_pseq,ln_pbias,nrg_ring_frag_total,cbmc_overlap,E_interfrag)
!*******************************************************************************
!
! PURPOSE: build the molecule from scratch
!
! First written by Jindal Shah on 07/18/07
!
! 05/25/08 (JS) : First committed to the repository
!
! 01/20/0 (JS) : Boltzmann weight of the trial is computed irrespective of
!                its energy. Previously, an energy cutoff was used to set
!                the Boltzmann weight
!
! 06/14/16 (RGM) : First fragment selected based on natoms_in_frag / natoms_in_all_frags
!
! DESCRIPTION: This subroutine performs the following steps
!
! Step 1) Select which fragment will be inserted first
! Step 2) Choose a conformation for the first fragment
! Step 3) Rotate the first fragment
! Step 4) Choose kappa_ins positions for the first fragment's COM
! Step 5) Choose kappa_dih orientations for each additional fragment
!
!*******************************************************************************

  USE Rotation_Routines
  USE IO_Utilities
  USE File_Names
  USE Energy_Routines

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  !DIR$ ATTRIBUTES ALIGN : 32 :: trial_atom_rp, trial_cell_coords, bitcell_overlap, xyz_rand_dp, nrg_sp_vec, xyz_rand, rtrial
  ! Arguments
  INTEGER :: this_im  ! molecule index
  INTEGER :: is       ! species index
  INTEGER :: this_box ! box index
  INTEGER, DIMENSION(1:nfragments(is)) :: frag_order
  REAL(DP), INTENT(IN) :: this_lambda ! fractional molecule?
  REAL(DP) :: randno
  REAL(DP) :: ln_pseq   ! probability of frag_order
  REAL(DP) :: ln_pbias  ! probability of placing each fragment in the box
  REAL(DP), INTENT(OUT) :: nrg_ring_frag_total ! potential energy of the
                                               ! isolated ring fragment
  LOGICAL, INTENT(INOUT) :: cbmc_overlap ! did all trials have core overlap?
  REAL(DP), INTENT(OUT), OPTIONAL :: E_interfrag ! total intramolecular, interfragment energy of molecule

  ! Local declarations
  INTEGER :: i, j, this_atom      ! atom indices
  INTEGER :: ifrag, this_fragment ! fragment indices
  INTEGER :: im                   ! molecule indices
  INTEGER :: total_frags  ! total number of conformations for this fragment in 
                          ! the reservoir
  INTEGER :: frag_start   ! random fragment to start growing from
  INTEGER :: frag_total   ! number of non-zero entries in frag_order
  INTEGER :: nl           ! number of the line where the x,y,x coords of the atom
                          ! of  the config and fragment randomly selected
                          ! were stored in the frag_position_library
  INTEGER, DIMENSION(nfragments(is)) :: live, frag_placed

  REAL(DP) :: x_anchor, y_anchor, z_anchor ! new COM coordinates for the 
                                           ! first fragment
  REAL(DP) :: xcom_old, ycom_old, zcom_old ! old COM coordinates for the
                                           ! first fragment
  REAL(DP) :: dx, dy, dz

  REAL(DP), DIMENSION(species_list(is)%kappa_ins_pad8) :: xcom_trial, ycom_trial, zcom_trial

  LOGICAL :: overlap ! TRUE if there is core overlap between a trial atom 
                     ! position and a an atom already in the box

  CHARACTER  :: this_file*120, symbol*1

  TYPE(Molecule_Class), POINTER :: this_molecule
  TYPE(Atom_Class), CONTIGUOUS, POINTER :: these_atoms(:)


  ! Variables associated with the CBMC part
  INTEGER :: itrial, trial, frag_type, ifrag_natoms
  INTEGER :: itrial_start

  REAL(DP) :: weight(species_list(is)%kappa_ins_pad8), rand_no, E_dihed
  REAL(DP) :: E_intra_vdw, E_intra_qq, E_inter_vdw, E_inter_qq, E_total
  REAL(DP) :: nrg(species_list(is)%kappa_ins_pad8), nrg_kBT, nrg_ring_frag

  LOGICAL :: del_overlap, overlap_trial(species_list(is)%kappa_ins_pad8)

  !Type(Atom_Class) :: rtrial(MAXVAL(natoms),0:MAX(kappa_ins,kappa_rot,kappa_dih))
  ! Slit pore variables

  LOGICAL :: framework_overlap
  REAL(DP) :: E_framework

  ! Inner-volume variables
  REAL(DP) :: radius, radius2, theta, phi, rsinphi, zmax, hzmax, inner_dz, zscale, log_inner_radius

!  ! DEBUGging variables
!  INTEGER :: M_XYZ_unit
  REAL(DP) :: overlap_time_s, overlap_time_e, overlap_time
  LOGICAL :: omp_flag
  LOGICAL :: need_max_dcom
  LOGICAL :: l_store_dp_trials, l_ortho

  LOGICAL :: l_get_bitcell

  INTEGER :: n_good_trials, n_good_trials_old

  REAL(DP) :: rtrial0(3,MAXVAL(frag_list(1:nfragments(is),is)%natoms))
  REAL(DP) :: rtrial(species_list(is)%kappa_ins_pad8,3,MAXVAL(frag_list(1:nfragments(is),is)%natoms))
  REAL(DP), DIMENSION(species_list(is)%kappa_ins_pad8,3) :: xyz_rand_dp
  REAL(SP), DIMENSION(species_list(is)%kappa_ins_pad8,3) :: xyz_rand
  REAL(DP) :: xl, hxl, yl, hyl, zl, hzl, length_dp(3,3), drxcom, drycom, drzcom, isp, max_dcomsq
  REAL(DP), DIMENSION(MAXVAL(frag_list(1:nfragments(is),is)%natoms)) :: drxcom_vec, drycom_vec, drzcom_vec
  REAL(DP), DIMENSION(3,MAXVAL(frag_list(1:nfragments(is),is)%natoms)) :: drcom_mat
  REAL(SP), DIMENSION(3,MAXVAL(frag_list(1:nfragments(is),is)%natoms)) :: dscom_mat
  REAL(SP) :: length_sp(3,3), dscom(3), dsxcom, dsycom, dszcom, sxcom, sycom, szcom
  REAL(SP) :: sxp, syp, szp, rsl, hrsl, rsp
  REAL(SP) :: rxp, ryp, rzp, this_atom_rp(3)
  REAL(4), DIMENSION(species_list(is)%kappa_ins_pad8,3,MAXVAL(frag_list(1:nfragments(is),is)%natoms)) :: trial_atom_rp
  INTEGER, DIMENSION(species_list(is)%kappa_ins_pad8,3,MAXVAL(frag_list(1:nfragments(is),is)%natoms)) :: trial_cell_coords

  REAL(SP) :: nrg_sp, nrg_sp_vec(species_list(is)%kappa_ins_pad8), overlap_nrg_sp
  REAL(DP) :: nrg_dp, cweight

  INTEGER :: zbcdf, zlbc, ybcdf, ylbc, xlbc, ia_frag
  INTEGER :: bitcell_bit, bitcell_int
  INTEGER, DIMENSION(species_list(is)%kappa_ins_pad8) :: which_good_trials
  INTEGER :: rlc, gtrial
  INTEGER, DIMENSION(3) :: this_atom_ci

  LOGICAL(1) :: bitcell_overlap(species_list(is)%kappa_ins_pad64,MAXVAL(frag_list(1:nfragments(is),is)%natoms))

  LOGICAL :: l_widom_cells, this_bitcell_overlap

  REAL(DP) :: overlap_nrg
  INTEGER :: i_dim, ia

  INTEGER :: bitcell_int_ior, bitcell_bit_min, bitcell_bit_max, bitcell_int1_min, bitcell_int1_max
  INTEGER :: bitcell_int2_min, bitcell_int2_max

  REAL(DP) :: t1, t2, t12, cbmc_time_s, cbmc_time_e, cbmc_time_e2, cbmc_time_e3, cbmc_time_e4
  REAL(DP) :: rng_time_s, rng_time_e, rng_time, cell_list_cbmc_nrg_time, trial_loop_time_s, trial_loop_time
  REAL(DP) :: cell_list_cbmc_nrg_time_s, cell_list_cbmc_nrg_time_e, noncell_cbmc_nrg_time_s, noncell_cbmc_nrg_time
  REAL(DP) :: noncell_cbmc_nrg_time_e

  INTEGER :: kappa_ins, kappa_ins_pad8, kappa_ins_pad64, kappa_dih, kappa_dih_pad8, kappa_dih_pad32

  REAL(DP) :: sz_rand_dp, zpart_width, zpart_shift
  LOGICAL, PARAMETER :: l_zscan = .FALSE.
  !INTEGER :: n_zparts
  !INTEGER, PARAMETER :: zpartition_kappa = 32 ! must be a multiple of 8 (or less than 8 to turn off z-partitioning)
  IF (widom_timing) THEN
          omp_flag = .FALSE.
          !$ omp_flag = .TRUE.
          IF (.NOT. omp_flag) CALL CPU_TIME(cbmc_time_s)
          !$ cbmc_time_s = omp_get_wtime()
  END IF
  kappa_ins = species_list(is)%kappa_ins
  kappa_ins_pad8 = species_list(is)%kappa_ins_pad8
  kappa_ins_pad64 = species_list(is)%kappa_ins_pad64
  kappa_dih = species_list(is)%kappa_dih
  kappa_dih_pad8 = species_list(is)%kappa_dih_pad8
  kappa_dih_pad32 = species_list(is)%kappa_dih_pad32


  !DIR$ ASSUME_ALIGNED trial_atom_rp:32, trial_cell_coords:32, bitcell_overlap:32, xyz_rand_dp:32
  !DIR$ ASSUME_ALIGNED nrg_sp_vec:32, xyz_rand:32, rtrial:32
  !DIR$ ASSUME (MOD(kappa_ins_pad8,8) .EQ. 0)
  !DIR$ ASSUME (MOD(kappa_ins_pad64,64) .EQ. 0)






  n_good_trials = species_list(is)%kappa_ins

  omp_flag = .FALSE.
  !$ omp_flag = .TRUE.

  cbmc_overlap = .TRUE.


  ! Initialize variables
  IF (widom_active) THEN
          this_molecule => widom_molecule
          these_atoms => widom_atoms
          l_store_dp_trials = .NOT. cbmc_cell_list_flag
  ELSE
          this_molecule => molecule_list(this_im,is)
          these_atoms => atom_list(:,this_im,is)
          l_store_dp_trials = .TRUE.
  END IF
  ln_pbias = 0.0_DP
  this_box = this_molecule%which_box ! which box this_im is in
  these_atoms%exist = .FALSE. ! mark all the atoms as deleted 
  this_molecule%frac = 0.0_DP
  frag_placed(:) = 0 ! =1 if fragment been placed
  nrg_ring_frag_total = 0.0_DP
  CALL Set_CBMC_Flag(.TRUE.)
  ! The energy of each trial coordinate will be stored in the array nrg.
  nrg(:) = 0.0_DP 
  overlap_trial(:) = .FALSE.
  need_max_dcom = .FALSE.

  weight(:) = 0.0_DP

  overlap_nrg = max_kBT / beta(this_box)
  overlap_nrg_sp = REAL(overlap_nrg,SP)

  l_ortho = box_list(this_box)%int_box_shape <= int_ortho
  l_widom_cells = widom_active .AND. l_sectors
  l_get_bitcell = widom_active .AND. bitcell_flag
  length_dp = box_list(this_box)%length
  length_sp = REAL(length_dp,4)

  !*****************************************************************************
  ! Step 1) Select which fragment will be inserted first
  !*****************************************************************************
  !
  ! One fragment will be inserted first, and then additional fragments will
  ! be added one at a time
  !

  ! get_fragorder is:
  !    *  .TRUE. when making an initial configuration, when inserting or 
  !       deleting a molecule in GCMC, when inserting a trial molecule to
  !       compute the chemical potential, or when transferring a molecule into 
  !       a box (GEMC insertion)
  !    *  .FALSE. when transfering a molecule out of a box (GEMC deletion), 
  !       since the frag_order was already selected as part of transferring the 
  !       molecule into the other box


  IF (get_fragorder) THEN
 
     ! Select a fragment to insert first
     randno = rranf()
     DO ifrag = 1, nfragments(is)
        IF (randno < frag_list(ifrag,is)%cum_prob_ins) EXIT
     END DO
     frag_start = ifrag
     ln_pseq = DLOG(frag_list(frag_start,is)%prob_ins) ! weighted prob of choosing ifrag
 
     ! Select the order fragments will be added to frag_start

     live(:) = 0
     frag_order(:) = 0
     frag_order(1) = frag_start ! this array will hold the order
     live(frag_start) = 1
     frag_total = 1
     
     ! If this molecule is made up of multiple fragments, select an order in 
     ! which the fragments will be grown
     IF (nfragments(is) > 1 ) THEN

        CALL Fragment_Order(frag_start,is,frag_total,frag_order,live,ln_pseq)
        
     END IF


  ELSE

     ! frag_order is already determined, just need how many fragments will be
     ! placed (answer: all of them)
     frag_total = nfragments(is)

     ! ln_pseq was also computed as part of the GEMC insertion move

  END IF

  !
  ! At this point, we have the order in which we will grow the molecule.

  !*****************************************************************************
  ! Step 2) Choose a conformation for the first fragment
  !*****************************************************************************
  !   
  ! A fragment library was previously generated via a single fragment MC 
  ! simulation. Since the fragments were sampled according to the Boltzmann 
  ! distribution, we can now pull from the reservoir with uniform probability.
  ! The reservoir of fragment conformations is stored in frag_coords.
  !

  ! If get_fragorder = .FALSE., then frag_start is not yet defined
  frag_start = frag_order(1)
  ifrag_natoms = frag_list(frag_start,is)%natoms

  ! If inserting a molecule, choose the first fragment's conformation from the 
  ! reservoir

  IF (.NOT. del_flag) THEN

     ! Pull from the reservoir with uniform probability
     !total_fragments is the number of configurations of the frag_start
     ! this_fragment is the randomly chosen configuration
     total_frags = frag_list(frag_start,is)%nconfig
     this_fragment = INT(rranf() * total_frags) + 1

    
     frag_type = frag_list(frag_start,is)%type
    
     ! Read the coordinates for every atom
     ! this_fragment is the randomly config selected

     DO i = 1, frag_list(frag_start,is)%natoms 
        
        this_atom = frag_list(frag_start,is)%atoms(i)
        
        nl = (frag_position_library(frag_type)-1) + &
                                           frag_list(frag_start,is)%natoms*(this_fragment -1) + i
        these_atoms(this_atom)%rp(1:3) = library_coords(1:3,nl)
        !these_atoms(this_atom)%rp(1) = &
        !                            !  frag_coords(i,this_fragment,frag_type)%rp(1)
        !                               library_coords(nl)%rp(1)
        !these_atoms(this_atom)%rp(2) = &
        !                             ! frag_coords(i,this_fragment,frag_type)%rp(2)
        !                                library_coords(nl)%rp(2)
        !these_atoms(this_atom)%rp(3) = &
        !                             ! frag_coords(i,this_fragment,frag_type)%rp(3)
        !                                library_coords(nl)%rp(3)
      END DO
  END IF
  ! Turn on the molecule and its individual atoms
  this_molecule%frac = this_lambda
  
  DO i =1, frag_list(frag_start,is)%natoms
     this_atom = frag_list(frag_start,is)%atoms(i)
     these_atoms(this_atom)%exist = .TRUE.
  END DO

  !*****************************************************************************
  ! Step 3) Rotate the first fragment
  !*****************************************************************************
  !
  ! At this time, only a single orientation is selected with uniform 
  ! probability. Future releases may attempt multiple orientations.
  ! 

  CALL Get_COM(this_im,is)
  
  ! If inserting the fragment, select a random orientation 
  IF ( .NOT. del_flag) CALL Rotate_Molecule_Eulerian(this_im,is)

  !*****************************************************************************
  ! Step 4) Choose kappa_ins positions for the first fragment's COM
  !*****************************************************************************
  ! 
  ! The first fragment will now be inserted into the simulation box. Multiple
  ! trial coordinates will be randomly generated. Each trial will be weighted
  ! by the Boltzmann factor of the change in potential energy. One trial will be
  ! selected from the weighted distribution.

  ! Store the conformation and orientation as 0th trial
  DO i = 1, frag_list(frag_start,is)%natoms

     this_atom = frag_list(frag_start,is)%atoms(i)
     
     !rtrial(this_atom,0)%rp(1) = these_atoms(this_atom)%rp(1)
     !rtrial(this_atom,0)%rp(2) = these_atoms(this_atom)%rp(2)
     !rtrial(this_atom,0)%rp(3) = these_atoms(this_atom)%rp(3)
     rtrial0(:,i) = these_atoms(this_atom)%rp(1:3)

  END DO

  ! Store the COM
  xcom_old = this_molecule%rcom(1)
  ycom_old = this_molecule%rcom(2)
  zcom_old = this_molecule%rcom(3)

  ! We will place this fragment based only on its external weight

  ! When is imreplace greater than 0?
  IF(imreplace .GT. 0) THEN
     CALL Compute_Max_COM_Distance(this_im,is)
     IF(.NOT. del_flag) THEN

        dx = molecule_list(imreplace,isreplace)%rcom(1) &
           - this_molecule%rcom(1)
        dy = molecule_list(imreplace,isreplace)%rcom(2) &
           - this_molecule%rcom(2)
        dz = molecule_list(imreplace,isreplace)%rcom(3) &
           - this_molecule%rcom(3)

        this_molecule%rcom(1) = molecule_list(imreplace,isreplace)%rcom(1)
        this_molecule%rcom(2) = molecule_list(imreplace,isreplace)%rcom(2)
        this_molecule%rcom(3) = molecule_list(imreplace,isreplace)%rcom(3)

        DO i = 1, frag_list(frag_start,is)%natoms

           this_atom = frag_list(frag_start,is)%atoms(i)
           these_atoms(this_atom)%rp(1) = rtrial0(1,i) + dx
           these_atoms(this_atom)%rp(2) = rtrial0(2,i) + dy
           these_atoms(this_atom)%rp(3) = rtrial0(3,i) + dz
        END DO

     END IF

  ELSE
     need_max_dcom = .TRUE.

     IF (widom_timing) THEN
             IF (.NOT. cbmc_flag) CALL CPU_TIME(rng_time_s)
             !$ rng_time_s = omp_get_wtime()
     END IF
     !CALL vector_rranf(TRANSFER(xyz_rand_dp,xyz_rand_dp)) ! Don't do this; it causes a wasteful memcopy to be generated
     IF (l_widom_cells .AND. l_zscan .AND. kappa_ins >= 8) THEN
             CALL array_boxscan_rranf(xyz_rand_dp,kappa_ins)
             !DIR$ VECTOR ALIGNED
             xyz_rand = REAL(xyz_rand_dp,SP)
     ELSE
             DO i_dim = 1, 3
                CALL vector_rranf(xyz_rand_dp(:,i_dim))
                IF (l_widom_cells) THEN
                        !DIR$ VECTOR ALIGNED
                        xyz_rand(:,i_dim) = REAL(xyz_rand_dp(:,i_dim),SP)
                END IF
                !DO itrial = 1, kappa_ins
                !        xyz_rand_dp(itrial,i_dim) = rranf()
                !END DO
             END DO
     END IF
     !IF (l_widom_cells) xyz_rand = REAL(xyz_rand_dp,SP)
     IF (widom_timing) THEN
             IF (.NOT. cbmc_flag) CALL CPU_TIME(rng_time_e)
             !$ rng_time_e = omp_get_wtime()
             rng_time = rng_time_e - rng_time_s
             rng_ins_time = rng_ins_time + rng_time
     END IF
     xl = box_list(this_box)%length(1,1)
     hxl = 0.5_DP * xl
     yl = box_list(this_box)%length(2,2)
     hyl = 0.5_DP * yl
     zl = box_list(this_box)%length(3,3)
     hzl = 0.5_DP * zl
     drxcom_vec = 0.0_DP
     drycom_vec = 0.0_DP
     drzcom_vec = 0.0_DP
     drxcom_vec(1:frag_list(frag_start,is)%natoms) = &
             rtrial0(1,1:frag_list(frag_start,is)%natoms) - xcom_old
     drycom_vec(1:frag_list(frag_start,is)%natoms) = &
             rtrial0(2,1:frag_list(frag_start,is)%natoms) - ycom_old
     drzcom_vec(1:frag_list(frag_start,is)%natoms) = &
             rtrial0(3,1:frag_list(frag_start,is)%natoms) - zcom_old
     IF (widom_timing) THEN
          IF (.NOT. omp_flag) CALL CPU_TIME(cbmc_time_e)
          !$ cbmc_time_e = omp_get_wtime()
          t12 = cbmc_time_e - cbmc_time_s
          cbmc_setup_ins_time = cbmc_setup_ins_time + t12
     END IF
     IF (species_list(is)%insertion == 'RESTRICTED' .AND. &
             box_list(this_box)%int_inner_shape /= int_none .AND. &
             .NOT. widom_active) THEN
        n_good_trials = kappa_ins
        SELECT CASE(box_list(this_box)%int_inner_shape)
        CASE(int_sphere)
                log_inner_radius = LOG(box_list(this_box)%inner_radius)
                DO itrial = 1, kappa_ins
                        theta = 2.0_DP*PI*xyz_rand_dp(itrial,1)
                        phi = ACOS(2.0_DP*xyz_rand_dp(itrial,2)-1.0_DP)
                        radius = EXP(LOG(xyz_rand_dp(itrial,3))/3.0_DP + log_inner_radius)
                        rsinphi = SIN(phi)*radius
                        xcom_trial(itrial) = COS(theta) * rsinphi
                        ycom_trial(itrial) = SIN(theta) * rsinphi
                        zcom_trial(itrial) = radius * COS(phi)
                END DO
           !radius2 = x_anchor**2 + y_anchor**2 + z_anchor**2
           !IF (radius2 > box_list(this_box)%inner_radius2) THEN
           !   theta = 2.0_DP*PI*rranf()
           !   phi = ACOS(2.0_DP*rranf()-1.0_DP)
           !   radius = (rranf())**(1.0_DP/3.0_DP) * box_list(this_box)%inner_radius
           !   x_anchor = radius * COS(theta) * SIN(phi)
           !   y_anchor = radius * SIN(theta) * SIN(phi)
           !   z_anchor = radius * COS(phi)
           !END IF
        CASE(int_cylinder)
                ! WARNING: this does not work for triclinic boxes
                DO itrial = 1, kappa_ins
                        theta = 2.0_DP*PI*xyz_rand_dp(itrial,1)
                        radius = SQRT(xyz_rand_dp(itrial,2)) * box_list(this_box)%inner_radius
                        xcom_trial(itrial) = radius * COS(theta)
                        ycom_trial(itrial) = radius * SIN(theta)
                        zcom_trial(itrial) = zl*xyz_rand_dp(itrial,3) - hzl 
                END DO
           !radius2 = x_anchor**2 + y_anchor**2
           !IF (radius2 > box_list(this_box)%inner_radius2) THEN
           !   theta = 2.0_DP*PI*rranf()
           !   radius = SQRT(rranf()) * box_list(this_box)%inner_radius
           !   x_anchor = radius * COS(theta)
           !   y_anchor = radius * SIN(theta)
           !END IF
        CASE(int_slitpore)
                zmax = box_list(this_box)%inner_zmax
                hzmax = 0.5_DP * zmax
                DO itrial = 1, kappa_ins
                        xcom_trial(itrial) = xyz_rand_dp(itrial,1)*xl - hxl
                        ycom_trial(itrial) = xyz_rand_dp(itrial,2)*yl - hyl
                        zcom_trial(itrial) = xyz_rand_dp(itrial,3)*zmax - hzmax
                END DO
           !IF (ABS(z_anchor) > box_list(this_box)%inner_zmax) THEN
           !   z_anchor = (0.5_DP - rranf()) * box_list(this_box)%inner_zmax
           !END IF
        CASE(int_interface)
                inner_dz = box_list(this_box)%inner_zmax - box_list(this_box)%inner_zmin
                DO itrial = 1, kappa_ins
                        xcom_trial(itrial) = xyz_rand_dp(itrial,1)*xl - hxl
                        ycom_trial(itrial) = xyz_rand_dp(itrial,2)*yl - hyl
                        zscale = xyz_rand_dp(itrial,3)*2.0_DP - 1.0_DP
                        zcom_trial(itrial) = SIGN(box_list(this_box)%inner_zmin,zscale) + zscale*inner_dz
                END DO
           !IF (ABS(z_anchor) > box_list(this_box)%inner_zmax .OR. &
           !    ABS(z_anchor) < box_list(this_box)%inner_zmin) THEN
           !   z_anchor = rranf() * (box_list(this_box)%inner_zmax - box_list(this_box)%inner_zmin)
           !   z_anchor = z_anchor + box_list(this_box)%inner_zmin
           !   IF (rranf() > 0.5_DP) THEN
           !      z_anchor = - z_anchor 
           !   END IF
           !END IF
        END SELECT
     ELSE! IF (l_widom_cells) THEN
             IF (l_get_bitcell) THEN
                     IF (widom_timing) THEN
                             IF (.NOT. omp_flag) CALL cpu_time(t1)
                             !$ t1 = omp_get_wtime()
                     END IF
                     zbcdf = box_list(this_box)%bitcell_dimfactor(3)
                     zlbc = box_list(this_box)%length_bitcells(3)
                     ybcdf = box_list(this_box)%bitcell_dimfactor(2)
                     ylbc = box_list(this_box)%length_bitcells(2)
                     xlbc = box_list(this_box)%length_bitcells(1)
                     !bitcell_int_ior = 0
                     !bitcell_bit_min = 100000
                     !bitcell_bit_max = -100000
                     !bitcell_int1_max = -(HUGE(bitcell_int1_min)-1)
                     !bitcell_int1_min = HUGE(bitcell_int1_max)
                     !bitcell_int2_max = -(HUGE(bitcell_int2_min)-1)
                     !bitcell_int2_min = HUGE(bitcell_int2_max)

                     IF (l_ortho) THEN
                             dscom_mat(1,:) = drxcom_vec/xl
                             dscom_mat(2,:) = drycom_vec/yl
                             dscom_mat(3,:) = drzcom_vec/zl
                     ELSE
                             drcom_mat(1,:) = drxcom_vec
                             drcom_mat(2,:) = drycom_vec
                             drcom_mat(3,:) = drzcom_vec
                             dscom_mat = REAL(MATMUL(box_list(this_box)%length_inv, &
                                     drcom_mat),SP)
                     END IF
                     SELECT CASE(ifrag_natoms)
                     CASE(1)
                             !DIR$ LOOP COUNT = 1000
                             !DIR$ VECTOR ALIGNED
                             !$OMP SIMD PRIVATE(sxp,syp,szp,bitcell_bit,bitcell_int,this_bitcell_overlap,sxcom,sycom,szcom) 
                             DO itrial = 1, kappa_ins
                                sxcom = xyz_rand(itrial,1)
                                sycom = xyz_rand(itrial,2)
                                szcom = xyz_rand(itrial,3)
                                sxp = sxcom + dscom_mat(1,1)
                                syp = sycom + dscom_mat(2,1)
                                szp = szcom + dscom_mat(3,1)
                                ! Note: it is important that the < checks and shifts come before the >= checks and shifts
                                !       due to floating point rounding.  It's fine for sxp, syp, or szp to be zero or even
                                !       just barely less than zero due to the behavior of INT(), but if they end up equal to
                                !       1.0, that can cause bitcell_int to be out of bounds or otherwise severely incorrect.
                                !       This was found out the hard way.
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = BTEST(bitcell_int,bitcell_bit)
                                bitcell_overlap(itrial,1) = .NOT. this_bitcell_overlap
                             END DO
                             !$OMP END SIMD
                     CASE(2)
                             !DIR$ LOOP COUNT = 1000
                             !DIR$ VECTOR ALIGNED
                             !$OMP SIMD PRIVATE(sxp,syp,szp,bitcell_bit,bitcell_int,this_bitcell_overlap,sxcom,sycom,szcom) 
                             DO itrial = 1, kappa_ins
                                sxcom = xyz_rand(itrial,1)
                                sycom = xyz_rand(itrial,2)
                                szcom = xyz_rand(itrial,3)
                                sxp = sxcom + dscom_mat(1,1)
                                syp = sycom + dscom_mat(2,1)
                                szp = szcom + dscom_mat(3,1)
                                ! Note: it is important that the < checks and shifts come before the >= checks and shifts
                                !       due to floating point rounding.  It's fine for sxp, syp, or szp to be zero or even
                                !       just barely less than zero due to the behavior of INT(), but if they end up equal to
                                !       1.0, that can cause bitcell_int to be out of bounds or otherwise severely incorrect.
                                !       This was found out the hard way.
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = BTEST(bitcell_int,bitcell_bit)
                                ! Iteration 2
                                sxp = sxcom + dscom_mat(1,2)
                                syp = sycom + dscom_mat(2,2)
                                szp = szcom + dscom_mat(3,2)
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = this_bitcell_overlap .OR. BTEST(bitcell_int,bitcell_bit)
                                ! Store result
                                bitcell_overlap(itrial,1) = .NOT. this_bitcell_overlap
                             END DO
                             !$OMP END SIMD
                     CASE(3)
                             !DIR$ LOOP COUNT = 1000
                             !DIR$ VECTOR ALIGNED
                             !$OMP SIMD PRIVATE(sxp,syp,szp,bitcell_bit,bitcell_int,this_bitcell_overlap,sxcom,sycom,szcom) 
                             DO itrial = 1, kappa_ins
                                sxcom = xyz_rand(itrial,1)
                                sycom = xyz_rand(itrial,2)
                                szcom = xyz_rand(itrial,3)
                                sxp = sxcom + dscom_mat(1,1)
                                syp = sycom + dscom_mat(2,1)
                                szp = szcom + dscom_mat(3,1)
                                ! Note: it is important that the < checks and shifts come before the >= checks and shifts
                                !       due to floating point rounding.  It's fine for sxp, syp, or szp to be zero or even
                                !       just barely less than zero due to the behavior of INT(), but if they end up equal to
                                !       1.0, that can cause bitcell_int to be out of bounds or otherwise severely incorrect.
                                !       This was found out the hard way.
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = BTEST(bitcell_int,bitcell_bit)
                                ! Iteration 2
                                sxp = sxcom + dscom_mat(1,2)
                                syp = sycom + dscom_mat(2,2)
                                szp = szcom + dscom_mat(3,2)
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = this_bitcell_overlap .OR. BTEST(bitcell_int,bitcell_bit)
                                ! Iteration 3
                                sxp = sxcom + dscom_mat(1,3)
                                syp = sycom + dscom_mat(2,3)
                                szp = szcom + dscom_mat(3,3)
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = this_bitcell_overlap .OR. BTEST(bitcell_int,bitcell_bit)
                                ! Store result
                                bitcell_overlap(itrial,1) = .NOT. this_bitcell_overlap
                             END DO
                             !$OMP END SIMD
                     CASE(4)
                             !DIR$ LOOP COUNT = 1000
                             !DIR$ VECTOR ALIGNED
                             !$OMP SIMD PRIVATE(sxp,syp,szp,bitcell_bit,bitcell_int,this_bitcell_overlap,sxcom,sycom,szcom) 
                             DO itrial = 1, kappa_ins
                                sxcom = xyz_rand(itrial,1)
                                sycom = xyz_rand(itrial,2)
                                szcom = xyz_rand(itrial,3)
                                sxp = sxcom + dscom_mat(1,1)
                                syp = sycom + dscom_mat(2,1)
                                szp = szcom + dscom_mat(3,1)
                                ! Note: it is important that the < checks and shifts come before the >= checks and shifts
                                !       due to floating point rounding.  It's fine for sxp, syp, or szp to be zero or even
                                !       just barely less than zero due to the behavior of INT(), but if they end up equal to
                                !       1.0, that can cause bitcell_int to be out of bounds or otherwise severely incorrect.
                                !       This was found out the hard way.
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = BTEST(bitcell_int,bitcell_bit)
                                ! Iteration 2
                                sxp = sxcom + dscom_mat(1,2)
                                syp = sycom + dscom_mat(2,2)
                                szp = szcom + dscom_mat(3,2)
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = this_bitcell_overlap .OR. BTEST(bitcell_int,bitcell_bit)
                                ! Iteration 3
                                sxp = sxcom + dscom_mat(1,3)
                                syp = sycom + dscom_mat(2,3)
                                szp = szcom + dscom_mat(3,3)
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = this_bitcell_overlap .OR. BTEST(bitcell_int,bitcell_bit)
                                ! Iteration 4
                                sxp = sxcom + dscom_mat(1,4)
                                syp = sycom + dscom_mat(2,4)
                                szp = szcom + dscom_mat(3,4)
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = this_bitcell_overlap .OR. BTEST(bitcell_int,bitcell_bit)
                                ! Store result
                                bitcell_overlap(itrial,1) = .NOT. this_bitcell_overlap
                             END DO
                             !$OMP END SIMD
                     CASE(5)
                             !DIR$ LOOP COUNT = 1000
                             !DIR$ VECTOR ALIGNED
                             !$OMP SIMD PRIVATE(sxp,syp,szp,bitcell_bit,bitcell_int,this_bitcell_overlap,sxcom,sycom,szcom) 
                             DO itrial = 1, kappa_ins
                                sxcom = xyz_rand(itrial,1)
                                sycom = xyz_rand(itrial,2)
                                szcom = xyz_rand(itrial,3)
                                sxp = sxcom + dscom_mat(1,1)
                                syp = sycom + dscom_mat(2,1)
                                szp = szcom + dscom_mat(3,1)
                                ! Note: it is important that the < checks and shifts come before the >= checks and shifts
                                !       due to floating point rounding.  It's fine for sxp, syp, or szp to be zero or even
                                !       just barely less than zero due to the behavior of INT(), but if they end up equal to
                                !       1.0, that can cause bitcell_int to be out of bounds or otherwise severely incorrect.
                                !       This was found out the hard way.
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = BTEST(bitcell_int,bitcell_bit)
                                ! Iteration 2
                                sxp = sxcom + dscom_mat(1,2)
                                syp = sycom + dscom_mat(2,2)
                                szp = szcom + dscom_mat(3,2)
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = this_bitcell_overlap .OR. BTEST(bitcell_int,bitcell_bit)
                                ! Iteration 3
                                sxp = sxcom + dscom_mat(1,3)
                                syp = sycom + dscom_mat(2,3)
                                szp = szcom + dscom_mat(3,3)
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = this_bitcell_overlap .OR. BTEST(bitcell_int,bitcell_bit)
                                ! Iteration 4
                                sxp = sxcom + dscom_mat(1,4)
                                syp = sycom + dscom_mat(2,4)
                                szp = szcom + dscom_mat(3,4)
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = this_bitcell_overlap .OR. BTEST(bitcell_int,bitcell_bit)
                                ! Iteration 5
                                sxp = sxcom + dscom_mat(1,5)
                                syp = sycom + dscom_mat(2,5)
                                szp = szcom + dscom_mat(3,5)
                                IF (sxp < 0.0) sxp = sxp + 1.0
                                IF (syp < 0.0) syp = syp + 1.0
                                IF (szp < 0.0) szp = szp + 1.0
                                IF (sxp >= 1.0) sxp = sxp - 1.0
                                IF (syp >= 1.0) syp = syp - 1.0
                                IF (szp >= 1.0) szp = szp - 1.0
                                bitcell_bit = INT(sxp*xlbc)
                                bitcell_int = ISHFT(bitcell_bit,-5) + &
                                        INT(syp*ylbc)*ybcdf + &
                                        INT(szp*zlbc)*zbcdf
                                bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                this_bitcell_overlap = this_bitcell_overlap .OR. BTEST(bitcell_int,bitcell_bit)
                                ! Store result
                                bitcell_overlap(itrial,1) = .NOT. this_bitcell_overlap
                             END DO
                             !$OMP END SIMD
                     CASE DEFAULT
                             DO ia_frag = 1, frag_list(frag_start,is)%natoms
                                     dsxcom = dscom_mat(1,ia_frag)
                                     dsycom = dscom_mat(2,ia_frag)
                                     dszcom = dscom_mat(3,ia_frag)
                                     !DIR$ LOOP COUNT = 1000
                                     !DIR$ VECTOR ALIGNED
                                     !$OMP SIMD PRIVATE(sxp,syp,szp,bitcell_bit,bitcell_int) 
                                     DO itrial = 1, kappa_ins
                                        sxp = xyz_rand(itrial,1) + dsxcom
                                        syp = xyz_rand(itrial,2) + dsycom
                                        szp = xyz_rand(itrial,3) + dszcom
                                        ! Note: it is important that the < checks and shifts come before the >= checks and shifts
                                        !       due to floating point rounding.  It's fine for sxp, syp, or szp to be zero or even
                                        !       just barely less than zero due to the behavior of INT(), but if they end up equal to
                                        !       1.0, that can cause bitcell_int to be out of bounds or otherwise severely incorrect.
                                        !       This was found out the hard way.
                                        IF (sxp < 0.0) sxp = sxp + 1.0
                                        IF (syp < 0.0) syp = syp + 1.0
                                        IF (szp < 0.0) szp = szp + 1.0
                                        IF (sxp >= 1.0) sxp = sxp - 1.0
                                        IF (syp >= 1.0) syp = syp - 1.0
                                        IF (szp >= 1.0) szp = szp - 1.0
                                        bitcell_bit = INT(sxp*xlbc)
                                        bitcell_int = ISHFT(bitcell_bit,-5) + &
                                                INT(syp*ylbc)*ybcdf + &
                                                INT(szp*zlbc)*zbcdf
                                        bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                                        bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                                        bitcell_overlap(itrial,ia_frag) = BTEST(bitcell_int,bitcell_bit)
                                     END DO
                                     !$OMP END SIMD
                             END DO
                             bitcell_overlap(:,1) = .NOT. ANY( &
                                     bitcell_overlap(:,1:frag_list(frag_start,is)%natoms),2) ! inverted so .FALSE. is overlap
                     END SELECT
                     !DO ia_frag = 1, frag_list(frag_start,is)%natoms
                     !        this_atom = frag_list(frag_start,is)%atoms(ia_frag)
                     !        drxcom = drxcom_vec(ia_frag)
                     !        drycom = drycom_vec(ia_frag)
                     !        drzcom = drzcom_vec(ia_frag)
                     !        IF (l_ortho) THEN
                     !           dsxcom = REAL(drxcom / xl,4)
                     !           dsycom = REAL(drycom / yl,4)
                     !           dszcom = REAL(drzcom / zl,4)
                     !           dscom = (/ dsxcom, dsycom, dszcom /)
                     !        ELSE
                     !           dscom = REAL(MATMUL(box_list(this_box)%length_inv, &
                     !                   (/ drxcom, drycom, drzcom /)),4)
                     !           dsxcom = dscom(1)
                     !           dsycom = dscom(2)
                     !           dszcom = dscom(3)
                     !        END IF
                     !        !!DIR$ LOOP COUNT = 1000
                     !        !!DIR$ VECTOR ALIGNED
                     !        !!$OMP SIMD PRIVATE(sxp,syp,szp,bitcell_bit,bitcell_int) 
                     !        !DO itrial = 1, kappa_ins
                     !        !   sxp = xyz_rand(itrial,1) + dsxcom
                     !        !   syp = xyz_rand(itrial,2) + dsycom
                     !        !   szp = xyz_rand(itrial,3) + dszcom
                     !        !   ! Note: it is important that the < checks and shifts come before the >= checks and shifts
                     !        !   !       due to floating point rounding.  It's fine for sxp, syp, or szp to be zero or even
                     !        !   !       just barely less than zero due to the behavior of INT(), but if they end up equal to
                     !        !   !       1.0, that can cause bitcell_int to be out of bounds or otherwise severely incorrect.
                     !        !   !       This was found out the hard way.
                     !        !   IF (sxp < 0.0) sxp = sxp + 1.0
                     !        !   IF (syp < 0.0) syp = syp + 1.0
                     !        !   IF (szp < 0.0) szp = szp + 1.0
                     !        !   IF (sxp >= 1.0) sxp = sxp - 1.0
                     !        !   IF (syp >= 1.0) syp = syp - 1.0
                     !        !   IF (szp >= 1.0) szp = szp - 1.0
                     !        !   bitcell_bit = INT(sxp*xlbc)
                     !        !   bitcell_int = ISHFT(bitcell_bit,-5) + &
                     !        !           INT(syp*ylbc)*ybcdf + &
                     !        !           INT(szp*zlbc)*zbcdf
                     !        !   bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                     !        !   !IF (bitcell_int > UBOUND(box_list(this_box)%bitcell_int32_vec,1)) THEN
                     !        !   !        WRITE(*,*) sxp, syp, szp, bitcell_bit
                     !        !   !        WRITE(*,*) xlbc, ylbc, zlbc
                     !        !   !        WRITE(*,*) ybcdf, zbcdf
                     !        !   !        WRITE(*,*) bitcell_bit
                     !        !   !END IF
                     !        !   bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                     !        !   bitcell_overlap(itrial,ia_frag) = BTEST(bitcell_int,bitcell_bit)
                     !        !END DO
                     !        !!$OMP END SIMD
                     !END DO
                     !bitcell_overlap(:,1) = .NOT. ANY( &
                     !        bitcell_overlap(:,1:frag_list(frag_start,is)%natoms),2) ! inverted so .FALSE. is overlap
                     n_good_trials = 0
                     !DIR$ LOOP COUNT = 1000
                     DO itrial = 1, kappa_ins
                        IF (bitcell_overlap(itrial,1)) THEN
                                n_good_trials = n_good_trials+1
                                which_good_trials(n_good_trials) = itrial
                                xyz_rand(n_good_trials,:) = xyz_rand(itrial,:)
                        END IF
                     END DO
                     IF (widom_timing) THEN
                             IF (.NOT. omp_flag) CALL cpu_time(t2)
                             !$ t2 = omp_get_wtime()
                             t12 = t2 - t1
                             bitcell_overlap_ins_time = bitcell_overlap_ins_time + t12
                             bitcell_overlap_ins_checks = bitcell_overlap_ins_checks + kappa_ins
                             bitcell_overlap_ins_overlaps = bitcell_overlap_ins_overlaps + (kappa_ins-n_good_trials)
                     END IF
                     IF (n_good_trials == 0) THEN
                        CALL Set_CBMC_Flag(.FALSE.)
                        RETURN
                     END IF
             ELSE
                     n_good_trials = kappa_ins
             END IF
             IF (l_widom_cells) THEN
                     IF (widom_timing) THEN
                             IF (.NOT. omp_flag) CALL cpu_time(t1)
                             !$ t1 = omp_get_wtime()
                     END IF
                     !DIR$ LOOP COUNT = 5
                     DO ia_frag = 1, frag_list(frag_start,is)%natoms
                             this_atom = frag_list(frag_start,is)%atoms(ia_frag)
                             drxcom = drxcom_vec(ia_frag)
                             drycom = drycom_vec(ia_frag)
                             drzcom = drzcom_vec(ia_frag)
                             IF (l_ortho) THEN
                                dsxcom = REAL(drxcom / length_dp(1,1),4)
                                dsycom = REAL(drycom / length_dp(2,2),4)
                                dszcom = REAL(drzcom / length_dp(3,3),4)
                                dscom = (/ dsxcom, dsycom, dszcom /)
                             ELSE
                                dscom = REAL(MATMUL(box_list(this_box)%length_inv, &
                                        (/ drxcom, drycom, drzcom /)),4)
                             END IF
                             DO i_dim = 1, 3
                                     rsl = MERGE(length_sp(i_dim,i_dim),1.0,l_ortho)
                                     hrsl = 0.5*rsl
                                     rlc = box_list(this_box)%real_length_cells(i_dim)
                                     DO itrial = 1, n_good_trials
                                        rsp = xyz_rand(itrial,i_dim) + dscom(i_dim)
                                        IF (rsp < 0.0) rsp = rsp + 1.0
                                        IF (rsp >= 1.0) rsp = rsp - 1.0
                                        trial_cell_coords(itrial,i_dim,ia_frag) = &
                                                INT(rsp*rlc) - box_list(this_box)%sectorbound(i_dim)
                                        rsp = rsp*rsl - hrsl
                                        trial_atom_rp(itrial,i_dim,ia_frag) = rsp
                                     END DO
                             END DO
                             IF (.NOT. l_ortho) THEN
                                     DO itrial = 1, n_good_trials
                                        sxp = trial_atom_rp(itrial,1,ia_frag)
                                        rxp = length_sp(1,1)*sxp
                                        ryp = length_sp(2,1)*sxp
                                        rzp = length_sp(3,1)*sxp
                                        syp = trial_atom_rp(itrial,2,ia_frag)
                                        rxp = rxp + length_sp(1,2)*syp
                                        ryp = ryp + length_sp(2,2)*syp
                                        rzp = rzp + length_sp(3,2)*syp
                                        szp = trial_atom_rp(itrial,3,ia_frag)
                                        rxp = rxp + length_sp(1,3)*szp
                                        ryp = ryp + length_sp(2,3)*szp
                                        rzp = rzp + length_sp(3,3)*szp
                                        trial_atom_rp(itrial,1,ia_frag) = rxp
                                        trial_atom_rp(itrial,2,ia_frag) = ryp
                                        trial_atom_rp(itrial,3,ia_frag) = rzp
                                     END DO
                             END IF
                     END DO
                     n_good_trials_old = n_good_trials
                     n_good_trials = 0
                     IF (widom_timing) THEN
                             IF (.NOT. omp_flag) CALL cpu_time(trial_loop_time_s)
                             !$ trial_loop_time_s = omp_get_wtime()
                             cell_list_cbmc_nrg_time = 0.0_DP
                     END IF
                     trial_loop: DO itrial = 1, n_good_trials_old
                        DO ia_frag = 1, frag_list(frag_start,is)%natoms
                                ia = frag_list(frag_start,is)%atoms(ia_frag)
                                this_atom_rp = trial_atom_rp(itrial,:,ia_frag)
                                this_atom_ci = trial_cell_coords(itrial,:,ia_frag)
                                IF (check_overlap(this_atom_rp(1:3),this_atom_ci,ia,is,this_box)) THEN
                                        CYCLE trial_loop
                                END IF
                        END DO
                        IF (cbmc_cell_list_flag) THEN
                                IF (widom_timing) THEN
                                        IF (.NOT. omp_flag) CALL cpu_time(cell_list_cbmc_nrg_time_s)
                                        !$ cell_list_cbmc_nrg_time_s = omp_get_wtime()
                                        cell_list_cbmc_nrg_ins_checks = cell_list_cbmc_nrg_ins_checks + 1
                                END IF
                                nrg_sp = SUM(Compute_Cell_List_CBMC_nrg(&
                                        trial_atom_rp(itrial,1,1:ifrag_natoms),trial_atom_rp(itrial,2,1:ifrag_natoms),& 
                                        trial_atom_rp(itrial,3,1:ifrag_natoms), &
                                        trial_cell_coords(itrial,1,1:ifrag_natoms),trial_cell_coords(itrial,2,1:ifrag_natoms), &
                                        trial_cell_coords(itrial,3,1:ifrag_natoms), &
                                        frag_list(frag_start,is)%atoms(1:ifrag_natoms),is,this_box))
                                IF (widom_timing) THEN
                                        IF (.NOT. omp_flag) CALL cpu_time(cell_list_cbmc_nrg_time_e)
                                        !$ cell_list_cbmc_nrg_time_e = omp_get_wtime()
                                        cell_list_cbmc_nrg_time = cell_list_cbmc_nrg_time + &
                                                (cell_list_cbmc_nrg_time_e - cell_list_cbmc_nrg_time_s)
                                END IF
                                IF (nrg_sp > overlap_nrg_sp) CYCLE trial_loop
                                !nrg_sp = 0.0
                                !DO ia_frag = 1, frag_list(frag_start,is)%natoms
                                !        ia = frag_list(frag_start,is)%atoms(ia_frag)
                                !        this_atom_rp = trial_atom_rp(itrial,:,ia_frag)
                                !        this_atom_ci = trial_cell_coords(itrial,:,ia_frag)
                                !        nrg_sp = nrg_sp + Compute_Cell_List_CBMC_nrg(this_atom_rp(1:3),this_atom_ci,ia,is,this_box)
                                !END DO
                                !nrg_sp_vec(n_good_trials) = nrg_sp
                        END IF
                        n_good_trials = n_good_trials + 1
                        which_good_trials(n_good_trials) = MERGE(which_good_trials(itrial),itrial,l_get_bitcell)
                        IF (cbmc_cell_list_flag) nrg_sp_vec(n_good_trials) = nrg_sp
                     END DO trial_loop 
                     IF (widom_timing) THEN
                             IF (.NOT. omp_flag) CALL cpu_time(t2)
                             !$ t2 = omp_get_wtime()
                             t12 = t2 - t1
                             trial_loop_time = t2 - trial_loop_time_s
                             trial_loop_ins_time = trial_loop_ins_time + trial_loop_time
                             cell_list_ins_time = cell_list_ins_time + t12
                             cell_list_cbmc_nrg_ins_time = cell_list_cbmc_nrg_ins_time + cell_list_cbmc_nrg_time
                             cell_list_ins_checks = cell_list_ins_checks + n_good_trials_old
                             !cell_list_cbmc_nrg_ins_checks = cell_list_cbmc_nrg_ins_checks + n_good_trials
                     END IF
                     IF (n_good_trials == 0) THEN
                        CALL Set_CBMC_Flag(.FALSE.)
                        RETURN
                     END IF
             END IF
             IF (l_ortho .AND. l_store_dp_trials) THEN
                     IF (widom_timing) THEN
                             IF (.NOT. omp_flag) CALL cpu_time(noncell_cbmc_nrg_time_s)
                             !$ noncell_cbmc_nrg_time_s = omp_get_wtime()
                             noncell_cbmc_nrg_time = 0.0_DP
                     END IF
                     IF (l_widom_cells) THEN
                             DO itrial = 1, n_good_trials
                                gtrial = which_good_trials(itrial)
                                xcom_trial(itrial) = xyz_rand_dp(gtrial,1)*xl - hxl
                                ycom_trial(itrial) = xyz_rand_dp(gtrial,2)*yl - hyl
                                zcom_trial(itrial) = xyz_rand_dp(gtrial,3)*zl - hzl
                             END DO
                     ELSE
                             xcom_trial(1:n_good_trials) = &
                                     xyz_rand_dp(1:n_good_trials,1)*xl - hxl
                             ycom_trial(1:n_good_trials) = &
                                     xyz_rand_dp(1:n_good_trials,2)*yl - hyl
                             zcom_trial(1:n_good_trials) = &
                                     xyz_rand_dp(1:n_good_trials,3)*zl - hzl
                     END IF
                     IF (widom_timing) THEN
                             IF (.NOT. omp_flag) CALL cpu_time(noncell_cbmc_nrg_time_e)
                             !$ noncell_cbmc_nrg_time_e = omp_get_wtime()
                             noncell_cbmc_nrg_time = noncell_cbmc_nrg_time + &
                                     (noncell_cbmc_nrg_time_e-noncell_cbmc_nrg_time_s)
                     END IF
             ELSE IF (l_store_dp_trials) THEN
                     IF (widom_timing) THEN
                             IF (.NOT. omp_flag) CALL cpu_time(noncell_cbmc_nrg_time_s)
                             !$ noncell_cbmc_nrg_time_s = omp_get_wtime()
                             noncell_cbmc_nrg_time = 0.0_DP
                     END IF
                     DO itrial = 1, n_good_trials
                        ! make sure optimization takes conditional out of loop
                        IF (l_widom_cells) THEN
                                gtrial = which_good_trials(itrial)
                        ELSE
                                gtrial = itrial
                        END IF
                        isp = xyz_rand_dp(gtrial,1) - 0.5_DP
                        x_anchor = length_dp(1,1)*isp
                        y_anchor = length_dp(2,1)*isp
                        z_anchor = length_dp(3,1)*isp
                        isp = xyz_rand_dp(gtrial,2) - 0.5_DP
                        x_anchor = x_anchor + length_dp(1,2)*isp
                        y_anchor = y_anchor + length_dp(2,2)*isp
                        z_anchor = z_anchor + length_dp(3,2)*isp
                        isp = xyz_rand_dp(gtrial,3) - 0.5_DP
                        x_anchor = x_anchor + length_dp(1,3)*isp
                        y_anchor = y_anchor + length_dp(2,3)*isp
                        z_anchor = z_anchor + length_dp(3,3)*isp
                        xcom_trial(itrial) = x_anchor
                        ycom_trial(itrial) = y_anchor
                        zcom_trial(itrial) = z_anchor
                     END DO
                     IF (widom_timing) THEN
                             IF (.NOT. omp_flag) CALL cpu_time(noncell_cbmc_nrg_time_e)
                             !$ noncell_cbmc_nrg_time_e = omp_get_wtime()
                             noncell_cbmc_nrg_time = noncell_cbmc_nrg_time + &
                                     (noncell_cbmc_nrg_time_e-noncell_cbmc_nrg_time_s)
                     END IF
             END IF
     END IF
     max_dcomsq = 0.0_DP
     DO i = 1, frag_list(frag_start,is)%natoms
        drxcom = drxcom_vec(i)
        drycom = drycom_vec(i)
        drzcom = drzcom_vec(i)
        drxcom = drxcom*drxcom + drycom*drycom + drzcom*drzcom
        max_dcomsq = MAX(max_dcomsq,drxcom)
     END DO
     this_molecule%rcom(4) = SQRT(max_dcomsq)
     !DIR$ ASSUME_ALIGNED nrg:32, weight:32, overlap_trial:32
     IF (l_store_dp_trials) THEN
             IF (widom_timing) THEN
                     IF (.NOT. omp_flag) CALL cpu_time(noncell_cbmc_nrg_time_s)
                     !$ noncell_cbmc_nrg_time_s = omp_get_wtime()
             END IF
             DO i = 1, frag_list(frag_start,is)%natoms
                rtrial(1:n_good_trials,1,i) = xcom_trial(1:n_good_trials) + drxcom_vec(i)
                rtrial(1:n_good_trials,2,i) = ycom_trial(1:n_good_trials) + drycom_vec(i)
                rtrial(1:n_good_trials,3,i) = zcom_trial(1:n_good_trials) + drzcom_vec(i)
             END DO
             IF (del_flag) THEN
                     rtrial(1,1:3,1:frag_list(frag_start,is)%natoms) = &
                             rtrial0(1:3,1:frag_list(frag_start,is)%natoms)
                     xcom_trial(1) = xcom_old
                     ycom_trial(1) = ycom_old
                     zcom_trial(1) = zcom_old
             END IF
             DO itrial = 1, n_good_trials
                ! Place the fragment (and all its atoms) at the trial coordinate
                DO i = 1, frag_list(frag_start,is)%natoms
                   this_atom = frag_list(frag_start,is)%atoms(i)
                   these_atoms(this_atom)%rp = rtrial(itrial,1:3,i)
                END DO

                this_molecule%rcom(1) = xcom_trial(itrial)
                this_molecule%rcom(2) = ycom_trial(itrial)
                this_molecule%rcom(3) = zcom_trial(itrial)

                !IF (need_max_dcom) THEN
                !        CALL Compute_Max_COM_Distance(this_im,is)
                !        need_max_dcom = .FALSE.
                !END IF

                ! Note that the COM position is always chosen inside the simulation box 
                ! so there is no need to call Fold_Molecule.
                   
                ! Calculate the intermolecular energy of the fragment. Note that
                ! cbmc_flag has been set to true so that the following call will compute
                ! interaction energy of the growing molecule within a small distance
                overlap = .FALSE.
                IF (widom_active) THEN
                        CALL Compute_Molecule_Nonbond_Inter_Energy_Widom(this_im,is,&
                                E_inter_vdw,overlap)
                        ! in this case, E_inter_vdw already includes qq energy
                        overlap = overlap .OR. E_inter_vdw >= overlap_nrg
                        overlap_trial(itrial) = overlap
                        cbmc_overlap = cbmc_overlap .AND. overlap
                        IF (.NOT. overlap) nrg(itrial) = E_inter_vdw
                ELSE
                        CALL Compute_Molecule_Nonbond_Inter_Energy(this_im,is,&
                                E_inter_vdw,E_inter_qq,overlap)
                        nrg_dp = MERGE(infinity_DP, E_inter_vdw + E_inter_qq, overlap)
                        overlap = overlap .OR. nrg_dp >= overlap_nrg
                        overlap_trial(itrial) = overlap
                        cbmc_overlap = cbmc_overlap .AND. overlap
                        IF (.NOT. overlap) nrg(itrial) = nrg_dp
                END IF
             END DO
             IF (widom_timing) THEN
                     IF (.NOT. omp_flag) CALL cpu_time(noncell_cbmc_nrg_time_e)
                     !$ noncell_cbmc_nrg_time_e = omp_get_wtime()
                     noncell_cbmc_nrg_time = noncell_cbmc_nrg_time + &
                             (noncell_cbmc_nrg_time_e-noncell_cbmc_nrg_time_s)
                     noncell_cbmc_nrg_ins_time = noncell_cbmc_nrg_ins_time + noncell_cbmc_nrg_time
                     nrg_ins_overlaps = nrg_ins_overlaps + COUNT(overlap_trial(1:n_good_trials))
             END IF
             ! Reject the move if all trials tripped overlap
             IF (cbmc_overlap) THEN
                CALL Set_CBMC_Flag(.FALSE.)
                RETURN
             END IF

             DO itrial = 1, n_good_trials
                nrg_dp = nrg(itrial)
                weight(itrial) = MERGE(0.0_DP,EXP(-beta(this_box)*nrg_dp),overlap_trial(itrial))
                !IF (overlap_trial(itrial)) THEN
                !        weight(itrial) = 0.0_DP
                !ELSE
                !        weight(itrial) = EXP(-beta(this_box)*nrg_dp)
                !END IF
             END DO
     ELSE
             !$OMP SIMD PRIVATE(nrg_dp,overlap)
             DO itrial = 1, n_good_trials
                nrg_dp = REAL(nrg_sp_vec(itrial),DP)
                nrg(itrial) = nrg_dp
                !overlap = nrg_dp >= overlap_nrg
                !cbmc_overlap = cbmc_overlap .AND. overlap
                weight(itrial) = EXP(-beta(this_box)*nrg_dp)
             END DO
             !$OMP END SIMD
             !! Reject the move if all trials tripped overlap
             !IF (cbmc_overlap) THEN
             !   CALL Set_CBMC_Flag(.FALSE.)
             !   RETURN
             !END IF
     END IF
     IF (widom_timing) THEN
          IF (.NOT. omp_flag) CALL CPU_TIME(cbmc_time_e2)
          !$ cbmc_time_e2 = omp_get_wtime()
          t12 = cbmc_time_e2 - cbmc_time_e
          cbmc_returnzone_ins_time = cbmc_returnzone_ins_time + t12
          cbmc_nonoverlap_ins_count = cbmc_nonoverlap_ins_count + 1_INT64
     END IF
     ! Store the cumulative weight of each trial
     cweight = weight(1)
     DO itrial = 2, n_good_trials
        cweight = cweight + weight(itrial)
        weight(itrial) = cweight
     END DO
        



     ! Select one of the trial coordinates

     IF (del_flag) THEN
        
        ! For a deletion move, we want the weight of the current position,
        ! which is stored in trial 1.
        trial = 1

     ELSE

        ! Choose one from Golden sampling for an insertion move
        rand_no = rranf() * cweight
     
        DO i = 1, n_good_trials
           IF ( rand_no < weight(i)) EXIT
        END DO
     
        trial = i

        IF ( trial == n_good_trials + 1 ) THEN
           ! None of the trials were picked. Could be due to the fact that all 
           ! the trials had a very small cumulative weight
           cbmc_overlap = .TRUE.
           CALL Set_CBMC_Flag(.FALSE.)
           RETURN
        END IF

     END IF

     ! Compute the weight of the selected trial coordinate
     ln_pbias = ln_pbias - beta(this_box) * nrg(trial) - DLOG(cweight)
  
     ! This line is not used
     !e_total = nrg(trial)

     ! We chose the ith trial coordinate for placement. Store the ith trial
     ! coordinates in the atom_list array. Note that for the deletion move,
     ! trial=1 has the current coordinates of the fragment, so the molecule
     ! is not moved.
     IF (l_store_dp_trials) THEN
             this_molecule%rcom(1) = xcom_trial(trial)
             this_molecule%rcom(2) = ycom_trial(trial)
             this_molecule%rcom(3) = zcom_trial(trial)
             DO i = 1, frag_list(frag_start,is)%natoms
                this_atom = frag_list(frag_start,is)%atoms(i)
                these_atoms(this_atom)%rp = rtrial(trial,1:3,i)
             
                !these_atoms(this_atom)%rp(1) = rtrial(this_atom,trial)%rp(1)
                !these_atoms(this_atom)%rp(2) = rtrial(this_atom,trial)%rp(2)
                !these_atoms(this_atom)%rp(3) = rtrial(this_atom,trial)%rp(3)
             
             END DO
     ELSE
             IF (l_ortho) THEN
                     this_molecule%rcom(1:3) = (xyz_rand_dp(which_good_trials(trial),1:3) - 0.5_DP) * &
                             (/ length_dp(1,1) , length_dp(2,2), length_dp(3,3) /)
             ELSE
                     this_molecule%rcom(1:3) = MATMUL(length_dp, &
                             xyz_rand_dp(which_good_trials(trial),1:3) - 0.5_DP)
             END IF
             DO i = 1, frag_list(frag_start,is)%natoms
                this_atom = frag_list(frag_start,is)%atoms(i)
                these_atoms(this_atom)%rp(1) = this_molecule%rcom(1) + drxcom_vec(i)
                these_atoms(this_atom)%rp(2) = this_molecule%rcom(2) + drycom_vec(i)
                these_atoms(this_atom)%rp(3) = this_molecule%rcom(3) + drzcom_vec(i)
             END DO
     END IF

  END IF
  
  ! Mark this fragment is placed
  frag_placed(frag_start) = 1

  ! We have our first segment placed in the system. 

  ! If this is a ring fragment,
  ! then we need intramolecular energies of the fragment used in biasing.
  ! For insertion, it is the energy in the gas phase. For deletion, it is 
  ! the energy of the old configuration. Note that this energy does not include
  ! angle bending energy as it cancels out from the acceptance rule.

  IF (frag_list(frag_start,is)%ring) THEN

     IF (del_flag) THEN
        ! compute the old energy
        CALL Compute_Ring_Fragment_Energy(frag_start,this_im,is,this_box,&
                nrg_ring_frag)
     ELSE
        nrg_ring_frag =  nrg_frag(frag_type)%this_config_energy(this_fragment) 
                        ! nrg_frag(this_fragment,frag_type)
     END IF
     nrg_ring_frag_total = nrg_ring_frag_total + nrg_ring_frag
     
  END IF

  !*****************************************************************************
  ! Step 5) Select kappa_dih orientations for each additional fragment
  !*****************************************************************************
  ! 
  ! Now we will place rest of the segments based on the initial fragment placed

  ! Why will the dihedral energy be anything other than zero?
  !CALL Compute_Molecule_Dihedral_Energy(this_im,is,E_dihed)
  !E_total = E_dihed

  ! If we've gotten this far, cbmc_overlap = FALSE
  cbmc_overlap = .FALSE.
  del_overlap = .FALSE.

  ! The first fragment in frag_order has already been placed.
  ! We will call Fragment_Placement to place the remaining fragments in
  ! frag_order, starting with fragment frag_order(2)
  frag_start = 2
  E_total = 0.0_DP
  IF (widom_timing) THEN
       IF (.NOT. omp_flag) CALL CPU_TIME(cbmc_time_e3)
       !$ cbmc_time_e3 = omp_get_wtime()
       t12 = cbmc_time_e3 - cbmc_time_e2
       cbmc_endzone_ins_time = cbmc_endzone_ins_time + t12
  END IF
  CALL Fragment_Placement(this_box,this_im,is,frag_start,frag_total, &
                          frag_order,frag_placed,this_lambda,E_total,ln_pbias, &
                          nrg_ring_frag_total,cbmc_overlap,del_overlap)
  IF (widom_timing) THEN
       IF (.NOT. omp_flag) CALL CPU_TIME(cbmc_time_e4)
       !$ cbmc_time_e4 = omp_get_wtime()
       t12 = cbmc_time_e4 - cbmc_time_e3
       cbmc_fragment_placement_time = cbmc_fragment_placement_time + t12
       IF (species_list(is)%nfragments > 1) THEN
               cbmc_dih_time = cbmc_dih_time + t12
               cbmc_dih_count = cbmc_dih_count + 1_INT64
       END IF
  END IF

  ! Note that cbmc_overlap may be TRUE and the cbmc_flag will be properly 
  ! assigned FALSE while exiting the code.

  ! Note that del_overlap may be true while computing the weight of
  ! an old configuration. So set cbmc_overlap to TRUE so that the 
  ! calling routine appropriately sets the coordinates of the atoms that were 
  ! not grown

  IF (del_overlap) cbmc_overlap = .TRUE.

  IF (PRESENT(E_interfrag)) E_interfrag = E_total

  ! Mark cbmc_flag as FALSE so that intermolecular nonbonded interactions
  ! are properly computed for the molecule

  CALL Set_CBMC_Flag(.FALSE.)
  
END SUBROUTINE Build_Molecule

!******************************************************************************
SUBROUTINE Build_Rigid_Fragment(this_im,is,this_box,frag_order,this_lambda, &
              ln_pseq,ln_pbias,nrg_ring_frag_total,cbmc_overlap)
!******************************************************************************
!
! The subroutine is based on Build_Molecule subroutine and 
! places the first fragment of the molecule. The first atom
! of the fragement in inserted in a biased manner and then orientation of the 
! the molecule is generated in a biased manner as well. Currently it will not 
! work for reactions, I will add imreplace at a later point of time. 
!
! Written by Neeraj Rai on 06/11/2010.
!
!********************************************************************************

  USE Rotation_Routines
  USE IO_Utilities
  USE File_Names
  USE Energy_Routines

  INTEGER :: this_im, is, this_box, first_atom
  INTEGER, DIMENSION(1:nfragments(is)) :: frag_order
  INTEGER :: nl, nlo      ! number of the line where the x,y,x coords of
                          ! config and fragment randomly selected
                          !were stored in the frag_position_library
  REAL(DP), INTENT(IN) :: this_lambda
  REAL(DP) :: ln_pseq
  REAL(DP) :: ln_pbias
  REAL(DP), INTENT(OUT) :: nrg_ring_frag_total
  LOGICAL, INTENT(INOUT) :: cbmc_overlap
  !-------------------------------

  TYPE(Molecule_Class), POINTER :: this_molecule
  TYPE(Atom_Class), POINTER :: these_atoms(:)


  INTEGER :: frag_start, frag_start_old, frag_total, i, this_atom
  INTEGER :: anchor, ifrag, total_connect, frag_connect
  INTEGER :: j, im
  INTEGER :: ii,jj,kk, total_frags, this_fragment
  INTEGER, DIMENSION(nfragments(is)) :: frag_placed, live
  REAL(DP) :: x_anchor, y_anchor, z_anchor, xcom_old, ycom_old, zcom_old
  REAL(DP) :: dx, dy, dz
  REAL(DP) :: x_this, y_this, z_this,tempx, tempy, tempz, temp_var, E_ang
  LOGICAL :: overlap
  CHARACTER  :: this_file*120, symbol*1
  ! Variables associated with the CBMC part

  INTEGER :: itrial, trial, frag_type
  REAL(DP) :: weight(MAX(species_list(is)%kappa_ins,species_list(is)%kappa_rot,species_list(is)%kappa_dih)), rand_no
  REAL(DP) :: e_dihed, E_intra_vdw, E_intra_qq, E_inter_vdw
  REAL(DP) :: E_inter_qq, E_total
  REAL(DP) :: nrg(MAX(species_list(is)%kappa_ins,species_list(is)%kappa_rot,species_list(is)%kappa_dih)), nrg_kBT, time0, time1, nrg_ring_frag
  LOGICAL :: del_overlap
  TYPE(Atom_Class) :: rtrial(MAXVAL(natoms),0:MAX(species_list(is)%kappa_ins,species_list(is)%kappa_rot,species_list(is)%kappa_dih))

  INTEGER :: kappa_ins, kappa_ins_pad8, kappa_ins_pad64, kappa_dih, kappa_dih_pad8, kappa_dih_pad32
  kappa_ins = species_list(is)%kappa_ins
  kappa_ins_pad8 = species_list(is)%kappa_ins_pad8
  kappa_ins_pad64 = species_list(is)%kappa_ins_pad64
  kappa_dih = species_list(is)%kappa_dih
  kappa_dih_pad8 = species_list(is)%kappa_dih_pad8
  kappa_dih_pad32 = species_list(is)%kappa_dih_pad32



  IF (widom_active) THEN
          this_molecule => widom_molecule
          these_atoms => widom_atoms
  ELSE
          this_molecule => molecule_list(this_im,is)
          these_atoms => atom_list(1:natoms(is),this_im,is)
  END IF



  weight(:)=0.0_DP
  CALL Set_CBMC_Flag(.TRUE.)
  ln_pbias = 0.0_DP
  
  ! Assign a locate number for the molecule
  this_box = this_molecule%which_box
  ! mark all the atoms as deleted 
  these_atoms%exist = .FALSE.

  frag_placed(:) = 0
  nrg_ring_frag_total = 0.0_DP

! Get the order of the fragment growth
  IF (get_fragorder) THEN
     frag_start = 1
     ln_pseq = 0.0_DP
     ! Obtain the order of fragment addition
     live(:) = 0
     frag_order(:) = 0
     frag_order(1) = frag_start
     live(frag_start) = 1
     frag_total = 1
     IF (nfragments(is) > 1 ) THEN
        CALL Fragment_Order(frag_start,is,frag_total,frag_order,live,ln_pseq)
     END IF
  ELSE
     frag_total = nfragments(is)
  END IF

  ! At this point, we have the order in which we will grow the molecule from
  ! Add the part to read in the coordinates from its file.
  ! we will make all the atoms of frag_start as part of the simulations

  these_atoms%exist = .FALSE.
  this_molecule%frac = 0.0_DP

  frag_start = frag_order(1)


  ! Note that we need to choose from the reservoir only when insertion
  ! is attempted

  this_molecule%frac = this_lambda

  ! NR: Select the location of first bead of the fragment

  first_atom = frag_list(frag_start,is)%atoms(1)
  these_atoms(first_atom)%exist = .true.


  DO itrial = 1, kappa_ins
     IF (del_flag .and.(itrial.eq.1)) THEN
        rtrial(first_atom,itrial)%rp(1) =  these_atoms(first_atom)%rp(1)
        rtrial(first_atom,itrial)%rp(2) =  these_atoms(first_atom)%rp(2)
        rtrial(first_atom,itrial)%rp(3) =  these_atoms(first_atom)%rp(3)
     ELSE
        IF (box_list(this_box)%int_box_shape == int_cubic) THEN
           
           rtrial(first_atom,itrial)%rp(1) = (0.5_DP - rranf()) * &
                box_list(this_box)%length(1,1)
           rtrial(first_atom,itrial)%rp(2) = (0.5_DP - rranf()) * &
                box_list(this_box)%length(2,2)
           rtrial(first_atom,itrial)%rp(3) = (0.5_DP - rranf()) * &
                box_list(this_box)%length(3,3)
           
        END IF

     ENDIF

  ! Store them in correct array for the energy call
     these_atoms(first_atom)%rp(1) = rtrial(first_atom,itrial)%rp(1)
     these_atoms(first_atom)%rp(2) = rtrial(first_atom,itrial)%rp(2)
     these_atoms(first_atom)%rp(3) = rtrial(first_atom,itrial)%rp(3)


     ! Compute the energy of this atom    
     this_molecule%rcom(1) = rtrial(first_atom,itrial)%rp(1) 
     this_molecule%rcom(2) = rtrial(first_atom,itrial)%rp(2) 
     this_molecule%rcom(3) = rtrial(first_atom,itrial)%rp(3)

     CALL Compute_Max_COM_Distance(this_im,is)

     CALL Compute_Molecule_Nonbond_Inter_Energy(this_im,is,E_inter_vdw,E_inter_qq,overlap)


!    WRITE(8,*) del_flag,this_box, itrial, this_im,E_inter_vdw,  E_inter_qq,overlap

  ! compute weight for itrial 

     IF (overlap) THEN
         ! the energy is too high, set the weight to zero
         weight(itrial) = 0.0_DP
     ELSE
         nrg_kBT = beta(this_box) * (E_inter_vdw + E_inter_qq )

!         nrg_kBT = beta(this_box) * (E_inter_vdw )
!!$        ! compute the weight of this trial for the reverse move and first
!!$        ! trial. The following IF-ELSE construct ensures that the weight
!!$        ! for the reverse move is always computed.
         IF ( del_flag .AND. (itrial == 1) ) THEN
            weight(itrial) = DEXP(-nrg_kBT)
         ELSE IF ( nrg_kBT >= max_kBT) THEN
            weight(itrial) = 0.0_DP
         ELSE
            weight(itrial) = DEXP(-nrg_kBT)
         END IF
      END IF
!     WRITE(8,*) del_flag, itrial, weight(itrial),nrg_kBT,beta(this_box)
      IF (itrial > 1 ) weight(itrial) = weight(itrial-1) + weight(itrial)
      
   ENDDO ! End the generation of trial locations and associating weight

  IF (weight(kappa_ins) == 0.0_DP) THEN
      cbmc_overlap = .TRUE.
      CALL Set_CBMC_Flag(.FALSE.)
      RETURN
  END IF

  IF (del_flag) THEN
      trial = 1
  ELSE
   ! Choose one from Golden sampling for an insertion move
      rand_no = rranf() * weight(kappa_ins)
      DO i = 1, kappa_ins
         IF ( rand_no < weight(i)) EXIT
      END DO
      trial = i
  END IF
!  Write(8,*) 'selected', trial

  ln_pbias = ln_pbias - beta(this_box) * nrg(trial) - DLOG(weight(kappa_ins))

! Now we have the position of the first bead of   
! Assign this position to proper these_atoms
  these_atoms(first_atom)%rp(1) = rtrial(first_atom,trial)%rp(1)
  these_atoms(first_atom)%rp(2) = rtrial(first_atom,trial)%rp(2)
  these_atoms(first_atom)%rp(3) = rtrial(first_atom,trial)%rp(3)

! NR: End selection of the first bead of the fragment


! Now set the atoms of the first fragment to true
  DO i =1, frag_list(frag_start,is)%natoms
     this_atom = frag_list(frag_start,is)%atoms(i)
     these_atoms(this_atom)%exist = .TRUE.
  END DO

! if it is insertion move, get the fragment from the reservoir

  IF (.NOT. del_flag) THEN
     ! obtain a random configuration
     total_frags = frag_list(frag_start,is)%nconfig
     ! Choose a fragment at random
     this_fragment = INT(rranf() * total_frags) + 1
     frag_type = frag_list(frag_start,is)%type
     DO i = 1, frag_list(frag_start,is)%natoms 
        this_atom = frag_list(frag_start,is)%atoms(i)
        nlo = (frag_position_library(frag_type)-1) + &
                                            frag_list(1,is)%natoms*(this_fragment-1) + 1
        nl = (frag_position_library(frag_type)-1) + &
                                           frag_list(1,is)%natoms*(this_fragment-1)+i

         rtrial(this_atom,0)%rp(1:3) = library_coords(1:3,nl) - &
                 library_coords(1:3,nlo) + &
                 these_atoms(first_atom)%rp(1:3)
         !rtrial(this_atom,0)%rp(1) = library_coords(nl)%rp(1) -&
         !                          library_coords(nlo)%rp(1) +&
         !                          these_atoms(first_atom)%rp(1)
         !rtrial(this_atom,0)%rp(2) = library_coords(nl)%rp(2) -&
         !                          library_coords(nlo)%rp(2) +& 
         !                          these_atoms(first_atom)%rp(2)
         !rtrial(this_atom,0)%rp(3) = library_coords(nl)%rp(3) -& 
         !                          library_coords(nlo)%rp(3) +& 
         !                          these_atoms(first_atom)%rp(3)
      !  rtrial(this_atom,0)%rp(1) = frag_coords(i,this_fragment,frag_type)%rp(1)-&
      !                            frag_coords(1,this_fragment,frag_type)%rp(1)+&
      !                            these_atoms(first_atom)%rp(1) 
      !  rtrial(this_atom,0)%rp(2) = frag_coords(i,this_fragment,frag_type)%rp(2)-&
      !                            frag_coords(1,this_fragment,frag_type)%rp(2)+&
      !                            these_atoms(first_atom)%rp(2)
      !  rtrial(this_atom,0)%rp(3) = frag_coords(i,this_fragment,frag_type)%rp(3)-&
      !                            frag_coords(1,this_fragment,frag_type)%rp(3)+&
      !                            these_atoms(first_atom)%rp(3) 
     END DO
  ELSE
     DO i = 1,frag_list(frag_start,is)%natoms
        this_atom = frag_list(frag_start,is)%atoms(i)
        rtrial(this_atom,0)%rp(1) = these_atoms(this_atom)%rp(1)
        rtrial(this_atom,0)%rp(2) = these_atoms(this_atom)%rp(2)
        rtrial(this_atom,0)%rp(3) = these_atoms(this_atom)%rp(3)
     END DO
  END IF

! Now fragment's first bead at the new location

! Start the rotational-bias
! Rotate this fragment about the first bead along different axes
! select one from trials

  DO itrial = 1, species_list(is)%kappa_rot
     IF(del_flag .and. (itrial.eq.1)) THEN
       DO i=1,frag_list(frag_start,is)%natoms
          this_atom = frag_list(frag_start,is)%atoms(i)
          rtrial(this_atom,itrial)%rp(1) = these_atoms(this_atom)%rp(1)
          rtrial(this_atom,itrial)%rp(2) = these_atoms(this_atom)%rp(2)
          rtrial(this_atom,itrial)%rp(3) = these_atoms(this_atom)%rp(3)
       END DO
       CALL Get_COM(this_im,is)
    ELSE
       DO i = 1,frag_list(frag_start,is)%natoms
          this_atom = frag_list(frag_start,is)%atoms(i)  
          these_atoms(this_atom)%rp(1) = rtrial(this_atom,0)%rp(1)
          these_atoms(this_atom)%rp(2) = rtrial(this_atom,0)%rp(2)
          these_atoms(this_atom)%rp(3) = rtrial(this_atom,0)%rp(3)
       ENDDO
       CALL Rotate_XYZ_Axes(this_im,is,frag_start,.true.,.true.,.true.,1)
       CALL Get_COM(this_im,is)
       !       CALL Fold_Molecule(this_im,is,this_box) 
       DO i = 1,frag_list(frag_start,is)%natoms
          this_atom = frag_list(frag_start,is)%atoms(i)
          rtrial(this_atom,itrial)%rp(1) = these_atoms(this_atom)%rp(1)
          rtrial(this_atom,itrial)%rp(2) = these_atoms(this_atom)%rp(2)
          rtrial(this_atom,itrial)%rp(3) = these_atoms(this_atom)%rp(3)
       END DO
       
    ENDIF
    !    Get COM
    CALL Compute_Max_COM_Distance(this_im,is)
    

!     Write(8,*) frag_list(frag_start,is)%natoms
!     Write(8,*)
!     DO i = 1,frag_list(frag_start,is)%natoms
!           this_atom = frag_list(frag_start,is)%atoms(i)
!           IF (this_atom .eq. 1) then 
!               Write(8,*) 'O', these_atoms(this_atom)%rp(1), &
!                               these_atoms(this_atom)%rp(2), &
!                               these_atoms(this_atom)%rp(3)
!           ELSE
!               Write(8,*) 'C', these_atoms(this_atom)%rp(1), &
!                               these_atoms(this_atom)%rp(2), &
!                               these_atoms(this_atom)%rp(3)
!           END IF
!     END DO  

!     these_atoms(first_atom)%exist = .false.



    overlap = .FALSE.
       
       CALL Compute_Molecule_Nonbond_Inter_Energy(this_im,is,E_inter_vdw,E_inter_qq,overlap)


!     these_atoms(first_atom)%exist = .true.

!     Write(8,*) 'First fragment attempt vdw qq',itrial, E_inter_vdw, E_inter_qq,this_box,del_flag,overlap

 ! compute weight for itrial 
     IF (overlap) THEN
         ! the energy is too high, set the weight to zero
         weight(itrial) = 0.0_DP
     ELSE
         nrg_kBT = beta(this_box) * (E_inter_vdw + E_inter_qq )
!         nrg_kBT = beta(this_box) * (E_inter_vdw )
!!$        ! compute the weight of this trial for the reverse move and first
!!$        ! trial. The following IF-ELSE construct ensures that the weight
!!$        ! for the reverse move is always computed.
         IF ( del_flag .AND. (itrial == 1) ) THEN
            weight(itrial) = DEXP(-nrg_kBT)
         ELSE IF ( nrg_kBT >= max_kBT) THEN
            weight(itrial) = 0.0_DP
         ELSE
            weight(itrial) = DEXP(-nrg_kBT)
         END IF
     END IF
     IF (itrial > 1 ) weight(itrial) = weight(itrial-1) + weight(itrial)
  ENDDO ! End the generation of trial orientations and associating waight

  IF (weight(species_list(is)%kappa_rot) == 0.0_DP) THEN
      cbmc_overlap = .TRUE.
      CALL Set_CBMC_Flag(.FALSE.)
      RETURN
  END IF

  IF (del_flag) THEN
      trial = 1
  ELSE
   ! Choose one from Golden sampling for an insertion move
      rand_no = rranf() * weight(species_list(is)%kappa_rot)
      DO i = 1, species_list(is)%kappa_rot
         IF ( rand_no < weight(i)) EXIT
      END DO
      trial = i
  END IF

!  Write(8,*) 'selected', trial

  ln_pbias = ln_pbias - beta(this_box) * nrg(trial) - DLOG(weight(species_list(is)%kappa_rot))

! Assign positions to these_atoms
  DO i=1,frag_list(frag_start,is)%natoms
     this_atom =  frag_list(frag_start,is)%atoms(i)
     these_atoms(this_atom)%rp(1) = rtrial(this_atom,trial)%rp(1)
     these_atoms(this_atom)%rp(2) = rtrial(this_atom,trial)%rp(2)
     these_atoms(this_atom)%rp(3) = rtrial(this_atom,trial)%rp(3)
  END DO

! IF ( .NOT. del_flag) THEN
! write(*,*) these_atoms%rp(1)
! write(*,*) these_atoms%rp(3)
! write(*,*) these_atoms%rp(2)
! NR: Now we have place the first fragment.
! write(*,*)
!	End if
  frag_placed(frag_start) = 1


  IF (frag_list(frag_start,is)%ring) THEN

     IF (del_flag) THEN
        ! compute the old energy
        CALL Compute_Ring_Fragment_Energy(frag_start,this_im,is,this_box,nrg_ring_frag)
     ELSE
        nrg_ring_frag =   nrg_frag(frag_type)%this_config_energy(this_fragment) 
                         !nrg_frag(this_fragment,frag_type)
     END IF
     nrg_ring_frag_total = nrg_ring_frag_total + nrg_ring_frag

  END IF

  ! Now we will place rest of the segments based on the initial fragment placed

  CALL Compute_Molecule_Dihedral_Energy(this_im,is,e_dihed)
  e_total = e_dihed

  CALL Fragment_Placement(this_box,this_im,is,2,frag_total,frag_order,frag_placed,this_lambda, &
       e_total,ln_pbias,nrg_ring_frag_total, cbmc_overlap, del_overlap)

  CALL Set_CBMC_Flag(.FALSE.)

END SUBROUTINE Build_Rigid_Fragment

!*******************************************************************************
SUBROUTINE Cut_Regrow(this_im,is,frag_live,frag_dead,frag_order,frag_total, &
     this_lambda,E_prev,ln_pseq,ln_pbias,nrg_ring_frag_tot,cbmc_overlap,del_overlap)
  !*****************************************************************************
  ! The subroutine cuts a bond, deletes part of a molecule and regrows using 
  ! configurational biasing. 
  !
  ! Written by Jindal Shah on 08/26/08
  !
  ! Step 1) Select a bond to cut with uniform probability & 
  !         delete fragments on one side of the bond
  ! Step 2) Compute energy
  ! Step 3) Regrow deleted fragments
  !
  !*****************************************************************************

  USE Global_Variables
  USE Random_Generators

  IMPLICIT NONE

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  ! Arguments
  INTEGER :: this_im    ! molecule index
  INTEGER :: is         ! species index
  INTEGER :: frag_live  ! the fragment with the cut bond that will be kept
  INTEGER :: frag_dead  ! the fragment with the cut bond that will be deleted
  INTEGER :: frag_order(1:nfragments(is)) ! order in which frags are added
  INTEGER :: frag_total ! number of entries in frag_order
  REAL(DP), INTENT(IN) :: this_lambda
  REAL(DP) :: E_prev
  REAL(DP) :: ln_pseq   ! probability of sequence in frag_order
  REAL(DP) :: ln_pbias  ! probability of placing each fragment in the box
  REAL(DP) :: nrg_ring_frag_tot
  LOGICAL :: cbmc_overlap
  LOGICAL :: del_overlap

  ! Local declarations
  INTEGER :: frag_bond, frag_start
  INTEGER :: this_box, i_frag, frag_no, i, this_atom, anchor, this_frag
  INTEGER :: anchor_live, anchor_dead
  
  INTEGER, DIMENSION(nfragments(is)) :: live ! fragment was not deleted OR 
                                             ! is in frag_order
  INTEGER, DIMENSION(nfragments(is)):: frag_placed ! fragment was not deleted

  REAL(DP) :: E_total, E_intra_vdw, E_intra_qq, E_inter_vdw, E_inter_qq, E_dihed

  LOGICAL :: overlap

  TYPE(Molecule_Class), POINTER :: this_molecule
  TYPE(Atom_Class), POINTER :: these_atoms(:)

  INTEGER :: kappa_ins, kappa_ins_pad8, kappa_ins_pad64, kappa_dih, kappa_dih_pad8, kappa_dih_pad32
  kappa_ins = species_list(is)%kappa_ins
  kappa_ins_pad8 = species_list(is)%kappa_ins_pad8
  kappa_ins_pad64 = species_list(is)%kappa_ins_pad64
  kappa_dih = species_list(is)%kappa_dih
  kappa_dih_pad8 = species_list(is)%kappa_dih_pad8
  kappa_dih_pad32 = species_list(is)%kappa_dih_pad32



  IF (widom_active) THEN
          this_molecule => widom_molecule
          these_atoms => widom_atoms
  ELSE
          this_molecule => molecule_list(this_im,is)
          these_atoms => atom_list(1:natoms(is),this_im,is)
  END IF


  ! Initialize variables
  CALL Set_CBMC_Flag(.TRUE.)
  del_overlap = .FALSE.
  nrg_ring_frag_tot = 0.0_DP
  ln_pseq = 0.0_DP
  ln_pbias = 0.0_DP

  !*****************************************************************************
  ! Step 1) Select a bond to cut with uniform probability & 
  !         delete fragments on one side of the bond
  !*****************************************************************************
  ! 
  ! obtain the box in which molecule is selected

  this_box = this_molecule%which_box

  ! Choose a fragment bond to cut

  ! When is del_flag == TRUE?
  IF ( .NOT. DEL_FLAG) THEN
     
     frag_bond = INT (rranf() * fragment_bonds(is)) + 1
     
     ! Kill one of the fragments
     
     IF ( rranf() < fragment_bond_list(frag_bond,is)%prob_del1) THEN
        
        frag_dead = fragment_bond_list(frag_bond,is)%fragment1
        frag_live = fragment_bond_list(frag_bond,is)%fragment2
        
     ELSE
        
        frag_dead = fragment_bond_list(frag_bond,is)%fragment2
        frag_live = fragment_bond_list(frag_bond,is)%fragment1
        
     END IF
     
     ! Technically, we should mark all fragments connected to frag_live as live
     ! (except for frag_dead) and then all fragments connected to those
     ! fragments as live, until we come to a terminal fragment. However,  
     ! Fragment_Order will only add dead fragments connected to frag_dead to 
     ! frag_order, so we can get away here with only marking frag_live as live
     live(:) = 0
     live(frag_live) = 1

     ! frag_order starts with frag_dead
     frag_order(:) = 0
     frag_order(1) = frag_dead
     frag_total = 1
     live(frag_dead) = 1 ! frag_dead is now live b/c it is in frag_order
     
     ! Obtain random order of fragments to be regrown
     CALL Fragment_Order(frag_dead,is,frag_total,frag_order,live,ln_pseq)


  END IF
     
  ! Now we need all the fragments that have been placed. To do this, mark all 
  ! the fragments as placed and unmark the ones from frag_order.
  these_atoms%exist = .TRUE.
  frag_placed(:) = 1
  
  DO i_frag = 1, frag_total

     this_frag = frag_order(i_frag)
     
     frag_placed(this_frag) = 0
     
     ! also delete all the atoms contained in this fragment
     
     DO i = 1, frag_list(this_frag,is)%natoms
        
        this_atom = frag_list(this_frag,is)%atoms(i)

        these_atoms(this_atom)%exist = .FALSE.
        
     END DO

  END DO
 
  ! Note that the anchor of frag_dead is present in the simulation so turn
  ! the exist flag true for this atom. In the process, we have also marked
  ! the anchor of frag_live as dead, so we will make that atom alive as well

  ! Get the two anchor atoms that need to be turned alive

  CALL Get_Common_Fragment_Atoms(is,frag_live,frag_dead,anchor_live,anchor_dead)

  these_atoms(anchor_live)%exist = .TRUE.
  these_atoms(anchor_dead)%exist = .TRUE.
 
  !*****************************************************************************
  ! Step 2) Compute energy
  !*****************************************************************************

  ! At this point we will calculate the energy of the fragment(s) in the 
  ! simulation. This is done so that the energy calculated for bonded and 
  ! nonbonded interactions do not cause trouble while growing the molecule. 
  ! Since dihedral, intramolecular nonbond and intermolecular nonbond energies 
  ! are used to bias the fragment placement, we will compute these energies. 
  ! Note that this is computed only when del_flag is false. We will use the 
  ! energy computed here when del_flag is true.

  IF (.NOT. del_flag) THEN

     CALL Compute_Molecule_Dihedral_Energy(this_im,is,E_dihed)
     E_prev = E_dihed

  END IF

  ! E_total is equal to the energy computed above during the growth phase.
  ! while it is the value obtained for cut_N_grow for calculating weight of the 
  ! old chain
  E_total = E_prev 

  !*****************************************************************************
  ! Step 3) Regrow deleted fragments
  !*****************************************************************************
  !
  ! We will call Fragment_Placement to place the fragments in frag_order, 
  ! starting with fragment frag_order(1)

  frag_start = 1
  CALL Fragment_Placement(this_box,this_im,is,frag_start,frag_total, &
       frag_order,frag_placed,this_lambda,E_total,ln_pbias,nrg_ring_frag_tot, &
       cbmc_overlap,del_overlap)
  CALL Set_CBMC_Flag(.FALSE.)
   

END SUBROUTINE Cut_Regrow



!*******************************************************************************
SUBROUTINE Fragment_Order(this_frag,is,frag_total,frag_order,live,ln_pseq)
!*******************************************************************************
!
! PURPOSE: Determine the order in which remaining fragments will be placed
!
! Written by Ryan Mullen on 05/08/14
!
! DESCRIPTION: This Subroutine determines the order in which the fragments will
! be (re)grown to form a completed molecule. 
!   * If the molecule is being inserted, this_frag is the only live fragment.
!     this_frag may be connected to only 1 other fragment (if this_frag is a 
!     terminal fragment) or more than 1 fragment (if this_frag is an internal
!     fragment).
!   * If the molecule is being partially regrown, a bond between two fragments
!     has been cut and part of the molecule deleted. this_frag is the only
!     deleted fragment that is connected to a live fragment. Whether this_frag
!     is terminal or internal, there is only 1 hanging connection initially in 
!     the molecule.
! 
!   Step 1) Determine the number & identity of hanging connections
!   Step 2) Select which fragment to add next
!   Step 3) Update the number & identity of hanging connections
!   
!   Loop thru steps 2 and 3 until all fragments are live.
!
!*******************************************************************************
  USE Global_Variables
  USE Random_Generators
  
  IMPLICIT NONE

  !*****************************************************************************
  ! Declare & Initialize Variables
  !*****************************************************************************

  ! Arguments
  INTEGER :: this_frag  ! fragment index
  INTEGER :: is         ! species index
  INTEGER :: frag_total ! running total of how many fragments are live
  INTEGER, DIMENSION(1:nfragments(is)) :: frag_order ! the order in which 
                                                     ! fragments will be grown
  INTEGER, DIMENSION(1:nfragments(is)) :: live       ! fragment in frag_order?
  REAL(DP) :: ln_pseq   ! the probability of frag_order

  ! Local declarations
  INTEGER :: ifrag, jfrag, ifrag_id, jfrag_id, i
  INTEGER :: nconnect ! number of entires in hanging_connections
  INTEGER :: hanging_connections(nfragments(is)) ! to hold ids of frags ready to add
  REAL(DP) :: randno, prob(nfragments(is)), cum_prob(nfragments(is))

  INTEGER :: kappa_ins, kappa_ins_pad8, kappa_ins_pad64, kappa_dih, kappa_dih_pad8, kappa_dih_pad32
  kappa_ins = species_list(is)%kappa_ins
  kappa_ins_pad8 = species_list(is)%kappa_ins_pad8
  kappa_ins_pad64 = species_list(is)%kappa_ins_pad64
  kappa_dih = species_list(is)%kappa_dih
  kappa_dih_pad8 = species_list(is)%kappa_dih_pad8
  kappa_dih_pad32 = species_list(is)%kappa_dih_pad32


  !*****************************************************************************
  !   Step 1) Determine the number & identity of hanging connections
  !*****************************************************************************

  ! If the molecule has more than 1 fragment, each fragment is connected to
  ! others. We need to know how many fragments connected to this_frag are not 
  ! yet live, and which fragments those are:

  nconnect = 0 ! to store the number of hanging connections
  hanging_connections(:) = 0 ! to store the frag ids that are ready to be added
  DO ifrag = 1, frag_list(this_frag,is)%nconnect
     ifrag_id = frag_list(this_frag,is)%frag_connect(ifrag)
     IF ( live(ifrag_id) == 0 ) THEN
        nconnect = nconnect + 1
        hanging_connections(nconnect) = ifrag_id
     END IF
  END DO

  ! Loop through steps 2 and 3 until all fragments are live.

  DO WHILE (nconnect > 0)
     !**************************************************************************
     !   Step 2) Select which fragment to add next
     !**************************************************************************
     ! Compute probabilty of adding next fragment for each frag in hanging_connections
     ! weight each fragment by the number of unique atoms in the fragment (natoms - 2)
     ! On the first loop, store the number of atoms in prob & cum_prob
     prob = 0
     cum_prob = 0
     DO ifrag = 1, nconnect
        ifrag_id = hanging_connections(ifrag)
        prob(ifrag) = frag_list(ifrag_id,is)%natoms - 2
        IF (ifrag == 1) THEN
          cum_prob(ifrag) = frag_list(ifrag_id,is)%natoms - 2
        ELSE
          cum_prob(ifrag) = cum_prob(ifrag-1)  + (frag_list(ifrag_id,is)%natoms - 2)
        END IF
     END DO
     ! On the second loop, divide the prob & cum_prob for each frag by the total number of atoms
     DO ifrag = 1, nconnect
        ifrag_id = hanging_connections(ifrag)
        prob(ifrag) = prob(ifrag) / cum_prob(nconnect)
        cum_prob(ifrag) = cum_prob(ifrag) / cum_prob(nconnect)
     END DO

     ! Choose a random fragment
     randno = rranf()
     DO ifrag = 1, nconnect
        IF (randno < cum_prob(ifrag)) EXIT
     END DO
     ifrag_id = hanging_connections(ifrag)
     ln_pseq = ln_pseq + DLOG( prob(ifrag) )

     ! Add ifrag_id to frag_order, make ifrag_id live
     frag_total = frag_total + 1
     frag_order(frag_total) = ifrag_id
     live(ifrag_id) = 1

     !**************************************************************************
     !   Step 3) Update the number & identity of hanging connections
     !**************************************************************************
     ! Remove ifrag_id from hanging_connections
     DO jfrag = ifrag, nconnect - 1
        hanging_connections(jfrag) = hanging_connections(jfrag + 1)
     END DO
     hanging_connections(nconnect) = 0
     nconnect = nconnect - 1

     ! Update nconnect and hanging_connections
     DO jfrag = 1, frag_list(ifrag_id,is)%nconnect
        jfrag_id = frag_list(ifrag_id,is)%frag_connect(jfrag)
        IF ( live(jfrag_id) == 0 ) THEN
           nconnect = nconnect + 1
           hanging_connections(nconnect) = jfrag_id
        END IF
     END DO

  END DO

END SUBROUTINE Fragment_Order

!*******************************************************************************
SUBROUTINE Fragment_Placement(this_box, this_im, is, frag_start, frag_total, &
           frag_order, frag_placed, this_lambda, e_total, ln_pbias, &
           nrg_ring_frag_tot, cbmc_overlap, del_overlap)
!*******************************************************************************
!
! PURPOSE: place the remaining fragments of the molecule
!
! DESCRIPTION: This subroutine performs the following steps
!
! Step 1) Select a fragment conformation
! Step 2) Align the fragment to the growing molecule
! Step 3) Select kappa_dih orientations for each additional fragment
! Step 4) Compute the energy of the fragment 
! Step 5) Select a trial dihedral using the weighted probabilities
!
!*******************************************************************************

  USE IO_Utilities

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  INTEGER, INTENT(IN) :: this_im, is, frag_start, frag_total, this_box
  INTEGER, DIMENSION(1:nfragments(is)), INTENT(IN) :: frag_order

  INTEGER, DIMENSION(1:nfragments(is)), INTENT(INOUT) :: frag_placed

  REAL(DP), INTENT(INOUT) :: e_total, ln_pbias, nrg_ring_frag_tot
  REAL(DP), INTENT(IN) :: this_lambda
  
  LOGICAL, INTENT(INOUT) :: cbmc_overlap, del_overlap

  ! local variable

  INTEGER :: i, ifrag, j, this_atom, frag_connect, total_connect, anchor_ifrag, is1, im, ia
  INTEGER :: anchor_frag_connect, atom_ifrag, atom_frag_connect, ii, trial
  INTEGER :: frag_type, dumcount

  INTEGER, DIMENSION(nfragments(is)) :: counted, connection

  INTEGER :: ispecies, jmol, k

  INTEGER :: nl_base ! index of library_coords column immediately before the start of the
                     ! coordinates of the chosen fragment conformation

  INTEGER :: total_frags, this_fragment, nfrag_atoms
  REAL(DP) :: x_this,y_this,z_this, vec1(3), vec2(3), aligner_ifrag(3,3)
  REAL(DP) :: hanger_ifrag(3,3), aligner_frag_connect(3,3), hanger_frag_connect(3,3)

  REAL(DP) :: tempx, tempy, tempz, theta, e_dihed
  REAL(DP) :: weight(species_list(is)%kappa_dih_pad8), nrg_dp_vec(species_list(is)%kappa_dih_pad8)
  REAL(DP) :: E_intra_qq, E_intra_vdw, prob_pick
  REAL(DP) :: temp_var, E_ang, E_inter_vdw, E_inter_qq
  REAL(DP) :: nrg_kBT, p_acc, nrg_intra_vdw, nrg_intra_qq, nrg_inter_vdw, nrg_inter_qq
  REAL(DP) :: trial_weight
  REAL(DP) :: nrg_ring_frag

  LOGICAL :: overlap, overlap_trial(species_list(is)%kappa_dih_pad8)
  LOGICAL, DIMENSION(natoms(is)) :: grown_exist, new_exist, combined_exist

  CHARACTER :: this_file*120, element*1 

  !TYPE(Atom_Class) :: config_list(natoms(is)) ! superceded by config_list_rp and new_rp
  !TYPE(Atom_Class) :: config_temp_list(natoms(is),species_list(is)%kappa_dih) ! superceded by trial_atom_rp

  TYPE(Molecule_Class), POINTER :: this_molecule
  TYPE(Atom_Class), POINTER, CONTIGUOUS :: these_atoms(:)
  INTEGER :: ifrag_outer, itrial, itrial_orig

  INTEGER :: ifrag_natoms
  INTEGER :: clrp_anchor_ifrag, clrp_anchor_frag_connect
  INTEGER, DIMENSION(natoms(is)) :: clrp_atom_id, atom_id, inv_atom_id

  INTEGER :: i_dim, ja, ia_new

  REAL(DP), DIMENSION(library_coords_dim1,natoms(is)) :: config_list_rp
  REAL(DP), DIMENSION(3,natoms(is)) :: new_rp
  REAL(SP), DIMENSION(species_list(is)%kappa_dih_pad8,3,natoms(is)) :: trial_atom_rp, trial_atom_rsp
  INTEGER, DIMENSION(species_list(is)%kappa_dih_pad8,3,natoms(is)) :: trial_cell_coords
  REAL(DP) :: theta_step, theta0, inv_nfrag_atoms_dp, max_dcomsq, jcharge_dp

  REAL(SP), DIMENSION(3,3) :: hanger_frag_connect_sp, inv_H_sp, H_sp
  INTEGER :: zbcdf,ybcdf,zlbc,ylbc,xlbc
  REAL(SP), DIMENSION(species_list(is)%kappa_dih_pad8) :: sintheta_sp_vec, costheta_sp_vec, nrg_sp_vec, nrg_dihed_vec
  REAL(DP) :: costheta0, sintheta0, costheta, sintheta
  REAL(SP) :: costheta0_sp, sintheta0_sp, sintheta_sp, costheta_sp, overlap_nrg, nrg_sp
  REAL(SP), DIMENSION(3) :: rp00, rp0, inv_l_vec, l_sp_vec, this_atom_rp
  REAL(SP) :: x0, y0, z0, tempy_sp, tempz_sp, rxp, ryp, rzp, hl_sp, irsp
  REAL(SP) :: sxp, syp, szp, rlc, l_sp, inv_l, rcutsq, jrp(3), vdwscale, chargescale, cfqq
  REAL(SP) :: rsq, drp, rsq_vec(species_list(is)%kappa_dih_pad8), rij
  REAL(SP) :: eps, sigsq, sigbyr2, sigbyr6, negsigsq, negsigbyr2, rterm
  REAL(SP) :: epsig_n, epsig_m, mie_n, mie_m, lnrsq, invr, alpha_ewald_sp
  REAL(SP) :: const1, const2, const3, const4, shift_p1, shift_p2, dsf_const, dsf_factor2_sp, sf_const1, sf_const2
  REAL(SP) :: roffsq_rijsq
  REAL(SP) :: anchor_frag_connect_rp(3), anchor_ifrag_rp(3), atom3_rp(3), m_sp(3)
  INTEGER :: atom1, atom2, atom3, atom4, dihedral_atoms(4), idihed, idihed_rb
  REAL(DP) :: r32_base(3), r32_dp(3), vecdp(3), m_dp(3), last_phi_dp, nrg_dp, cweight
  REAL(SP) :: m_normsq, rx32, ry32, rz32, rx34, ry34, rz34, ni, cosphi, invnorm
  REAL(SP) :: rb_c(5), cosphi2, cosphi3, cosphi4, cosphi5
  REAL(SP), DIMENSION(species_list(is)%kappa_dih_pad8) :: cosphi_vec, phi_vec
  REAL(SP) :: last_phi, phi_step_sp, theta_step_sp, a0, phi
  REAL(DP) :: sinlintheta, coslintheta, sincosmat(2,2)
  INTEGER :: bitcell_bit, bitcell_int, this_atom_ci(3), jtype, itype
  INTEGER :: this_int_vdw_style, this_int_vdw_sum_style, this_int_charge_sum_style
  LOGICAL :: l_ortho, l_widom_cells, l_get_bitcell, ij_get_coul
  LOGICAL(1), DIMENSION(species_list(is)%kappa_dih_pad32,natoms(is)) :: bitcell_overlap
  INTEGER :: n_good_trials, n_good_trials_pad8, n_good_trials_old
  INTEGER :: anchor_atoms(2), anchor_atoms_reversed(2)
  INTEGER, DIMENSION(species_list(is)%kappa_dih_pad8) :: which_good_trials

  REAL(DP) :: rijsq_dp, Eij_intra_vdw, Eij_intra_qq, Eij_inter_qq, Eij_inter_vdw
  LOGICAL :: get_vdw, get_qq, l_skip_dihed_vec(species_list(is)%ndihedrals_energetic)
  INTEGER :: interfrag_nonbond_pairlist(2,natoms(is)*natoms(is)/(2*species_list(is)%nfragments))
  INTEGER :: n_interfrag_nonbond_pairs

  REAL(DP) :: nrg_vdw_dp, nrg_qq_dp

  INTEGER :: kappa_ins, kappa_ins_pad8, kappa_ins_pad64, kappa_dih, kappa_dih_pad8, kappa_dih_pad32
  kappa_ins = species_list(is)%kappa_ins
  kappa_ins_pad8 = species_list(is)%kappa_ins_pad8
  kappa_ins_pad64 = species_list(is)%kappa_ins_pad64
  kappa_dih = species_list(is)%kappa_dih
  kappa_dih_pad8 = species_list(is)%kappa_dih_pad8
  kappa_dih_pad32 = species_list(is)%kappa_dih_pad32


  IF (widom_active) THEN
          this_molecule => widom_molecule
          these_atoms => widom_atoms
  ELSE
          this_molecule => molecule_list(this_im,is)
          these_atoms => atom_list(1:natoms(is),this_im,is)
  END IF
  overlap_nrg = max_kBT/beta(this_box)
  this_int_vdw_sum_style = int_vdw_sum_style(this_box)
  n_interfrag_nonbond_pairs = 0

  SELECT CASE(this_int_vdw_sum_style)
  CASE(vdw_cut_switch)
          const1 = switch_factor1(this_box)*switch_factor2(this_box)
          const2 = switch_factor1(this_box)*2
          const3 = ron_switch_sq(this_box)
          const4 = roff_switch_sq(this_box)
  END SELECT


!  ! DEBUGging variables
!  INTEGER :: M_XYZ_unit


  l_ortho = box_list(this_box)%int_box_shape <= int_ortho
  l_widom_cells = widom_active .AND. l_sectors
  l_get_bitcell = widom_active .AND. bitcell_flag
  config_list_rp = 0.0_DP
  dumcount = 0
  weight(:) = 0.0_DP

  l_skip_dihed_vec = .FALSE. ! currently marks dihedral as already being computed and able to be skipped early
  ! l_skip_dihed_vec gets inverted later to mark which dihedrals aren't interfragment dihedrals
  !DIR$ LOOP COUNT = 0, 1, 2, 3
  DO ifrag_outer = frag_start, frag_total

     ifrag = frag_order(ifrag_outer)

     grown_exist = these_atoms%exist
     ! Note that there has to be only one connection of ifrag that is already 
     ! placed. Let us find out which fragment that is

     total_connect = 0
     connection(:) = 0

     DO j = 1, frag_list(ifrag,is)%nconnect
        
        frag_connect = frag_list(ifrag,is)%frag_connect(j)
        
        IF (frag_placed(frag_connect) == 1 ) THEN
           total_connect = total_connect + 1
           connection(total_connect) = frag_connect
        END IF
        
     END DO

     IF (total_connect > 1 ) THEN
        err_msg = ''
        err_msg(1) = 'More than one connections of' // &
                     TRIM(Int_To_String(ifrag)) // 'exist'
        CALL Clean_Abort(err_msg,'Fragment_Placement')
     END IF
     
     ! If here then only one connection found

     frag_connect = connection(1)
     ! find anchor atom ids for both the fragments

     CALL Get_Common_Fragment_Atoms(is,ifrag,frag_connect,anchor_ifrag, &
             anchor_frag_connect)
     nfrag_atoms = 0
     atom_id(:) = 0
     ifrag_natoms = frag_list(ifrag,is)%natoms
     
     DO j = 1, ifrag_natoms
       
        this_atom = frag_list(ifrag,is)%atoms(j)
        IF (this_atom .EQ. anchor_ifrag) THEN
                clrp_anchor_ifrag = j
        ELSE IF (this_atom .EQ. anchor_frag_connect) THEN
                clrp_anchor_frag_connect = j
        ELSE
                nfrag_atoms = nfrag_atoms + 1
                atom_id(nfrag_atoms) = this_atom
                clrp_atom_id(nfrag_atoms) = j
        END IF
        
        !IF ( (this_atom /= anchor_ifrag) .AND. &
        !     (this_atom /= anchor_frag_connect)) THEN
        !   nfrag_atoms = nfrag_atoms + 1
        !   atom_id(nfrag_atoms) = this_atom
        !END IF
        
     END DO
    
     ! Note that frag_connect already has anchor of ifrag placed both for fixed 
     ! and variable bond length cases so all we have to do is obtain the 
     ! coordinates of these atoms in the configuration and align it to their 
     ! coordinates in the simulation box.
     
     !**************************************************************************
     ! Step 1) Select a fragment conformation
     !**************************************************************************
     !
     IF ( del_flag) THEN

        ! For a deletion move, use the existing conformation for each additional
        ! fragment

        DO j = 1, ifrag_natoms
           this_atom = frag_list(ifrag,is)%atoms(j)

           these_atoms(this_atom)%exist = .TRUE.
           config_list_rp(1:3,j) = these_atoms(this_atom)%rp(1:3)

           !config_list(this_atom)%rp(1) = these_atoms(this_atom)%rp(1)
           !config_list(this_atom)%rp(2) = these_atoms(this_atom)%rp(2)
           !config_list(this_atom)%rp(3) = these_atoms(this_atom)%rp(3)
        END DO

        ! For a ring fragment, calculate the fragment intramolecular energy

        IF (frag_list(ifrag,is)%ring) THEN
           CALL Compute_Ring_Fragment_Energy(ifrag,this_im,is,this_box, &
                   nrg_ring_frag)
           nrg_ring_frag_tot = nrg_ring_frag_tot + nrg_ring_frag
        END IF

     ELSE

        ! Select a fragment conformation from the reservoir
        ! The reservoir was populated with a Boltzmann distribution, so now we
        ! can pull from it with a uniform probability

        total_frags = frag_list(ifrag,is)%nconfig
        this_fragment = INT(rranf() * total_frags) + 1
        
        ! Read in the coordinates
        frag_type = frag_list(ifrag,is)%type

        nl_base = (frag_position_library(frag_type)-1) + &
                     frag_list(ifrag,is)%natoms*(this_fragment-1)
        config_list_rp(:,1:ifrag_natoms) = &
                library_coords(:,nl_base+1:nl_base+ifrag_natoms)
        these_atoms(frag_list(ifrag,is)%atoms(1:ifrag_natoms))%exist = .TRUE.
        !DO j = 1, frag_list(ifrag,is)%natoms
        !   
        !   this_atom = frag_list(ifrag,is)%atoms(j)
        !   
        !   these_atoms(this_atom)%exist = .TRUE.
        !   config_list(this_atom)%rp(1) = &
        !         library_coords(nl)%rp(1)
        !      !frag_coords(j,this_fragment,frag_type)%rp(1)
        !   config_list(this_atom)%rp(2) = &
        !         library_coords(nl)%rp(2)
        !      !frag_coords(j,this_fragment,frag_type)%rp(2)
        !   config_list(this_atom)%rp(3) = &
        !         library_coords(nl)%rp(3)
        !      !frag_coords(j,this_fragment,frag_type)%rp(3)
        !END DO
     
        ! For a ring fragment, access the fragment intramolecular energy 

        IF (frag_list(ifrag,is)%ring) THEN
           nrg_ring_frag =  nrg_frag(frag_type)%this_config_energy(this_fragment)
                              ! nrg_frag(this_fragment,frag_type)
           nrg_ring_frag_tot = nrg_ring_frag_tot + nrg_ring_frag
        END IF
!why was this ever here?        CLOSE(UNIT=11)

     END IF
     combined_exist = these_atoms%exist
     new_exist = combined_exist .AND. .NOT. grown_exist
     these_atoms%exist = new_exist
     
     !**************************************************************************
     ! Step 2) Align the fragment to the growing molecule
     !**************************************************************************
     !

     
     
!     anchor_ifrag = frag_list(ifrag,is)%anchor(1)
!     anchor_frag_connect = frag_list(frag_connect,is)%anchor(1)

     ! Find one atom of ifrag and frag_connect that will be used for generating 
     ! xy plane

     !DO j = 1, frag_list(ifrag,is)%natoms

     !   this_atom = frag_list(ifrag,is)%atoms(j)

     !   IF ( (this_atom /= anchor_ifrag) .AND. &
     !        (this_atom /= anchor_frag_connect)) EXIT

     !END DO

     !atom_ifrag = this_atom

     ! Similarly for frag_connect

     DO j = 1, frag_list(frag_connect,is)%natoms

        this_atom = frag_list(frag_connect,is)%atoms(j)
     
        IF ( (this_atom /= anchor_ifrag) .AND. &
             (this_atom /= anchor_frag_connect)) EXIT

     END DO

     
     atom_frag_connect = this_atom

     ! Now use three atoms to obtain aligner and hanger matrix for the two 
     ! fragments: 

     ! atom1 == origin
     ! atom2 == id of the atom along which +ve x - axis is aligned
     ! atom3 == helps to obtain y-axis so that atom1-atom2-atom3 define xy plane

     ! vec1 == r_atom2 - r_atom1
     ! vec2 == r_atom3 - r_atom1

     ! form these vectors from the configuration read in from reservoir

     ! for ifrag

     new_rp(1:3,1:nfrag_atoms) = config_list_rp(1:3,clrp_atom_id(1:nfrag_atoms))

     vec1(1:3) = config_list_rp(1:3,clrp_anchor_ifrag) - config_list_rp(1:3,clrp_anchor_frag_connect)
     vec2(1:3) = new_rp(1:3,1) - config_list_rp(1:3,clrp_anchor_frag_connect)

     !vec1(1) = config_list(anchor_ifrag)%rp(1) &
     !        - config_list(anchor_frag_connect)%rp(1)
     !vec1(2) = config_list(anchor_ifrag)%rp(2) &
     !        - config_list(anchor_frag_connect)%rp(2)
     !vec1(3) = config_list(anchor_ifrag)%rp(3) &
     !        - config_list(anchor_frag_connect)%rp(3)


     !vec2(1) = config_list(atom_ifrag)%rp(1) &
     !        - config_list(anchor_frag_connect)%rp(1)
     !vec2(2) = config_list(atom_ifrag)%rp(2) &
     !        - config_list(anchor_frag_connect)%rp(2)
     !vec2(3) = config_list(atom_ifrag)%rp(3) &
     !        - config_list(anchor_frag_connect)%rp(3)
       

     CALL Get_Aligner_Hanger(vec1, vec2, aligner_ifrag,hanger_ifrag)
     
     ! Calculate this only for inserting a molecule
     ! The following line was moved out of the following IF statement because it's used again later
     r32_base = these_atoms(anchor_ifrag)%rp(1:3) - &
             these_atoms(anchor_frag_connect)%rp(1:3)
     IF ( .NOT. del_flag) THEN
        
        
        vec1(1:3) = r32_base
        !vec1(1) = these_atoms(anchor_ifrag)%rp(1) - &
        !     these_atoms(anchor_frag_connect)%rp(1)
        !
        !vec1(2) = these_atoms(anchor_ifrag)%rp(2) - &
        !     these_atoms(anchor_frag_connect)%rp(2)
        !
        !vec1(3) = these_atoms(anchor_ifrag)%rp(3) - &
        !     these_atoms(anchor_frag_connect)%rp(3)
        vec2(1:3) = these_atoms(atom_frag_connect)%rp(1:3) - &
                these_atoms(anchor_frag_connect)%rp(1:3)
        
        !vec2(1) = these_atoms(atom_frag_connect)%rp(1) - &
        !     these_atoms(anchor_frag_connect)%rp(1)
        !
        !vec2(2) = these_atoms(atom_frag_connect)%rp(2) - &
        !     these_atoms(anchor_frag_connect)%rp(2)
        !
        !vec2(3) = these_atoms(atom_frag_connect)%rp(3) - &
        !     these_atoms(anchor_frag_connect)%rp(3)     
     

        CALL Get_Aligner_Hanger(vec1, vec2, aligner_frag_connect,hanger_frag_connect)

     END IF
        
     ! Apply aligner of ifrag and then hanger of frag_connect to join the two 
     ! fragments.
     
     ! Apply aligner of ifrag
     

     !tempx = config_list(anchor_frag_connect)%rp(1)
     !tempy = config_list(anchor_frag_connect)%rp(2)
     !tempz = config_list(anchor_frag_connect)%rp(3)
     !
     !config_list(:)%rp(1) = config_list(:)%rp(1) - tempx
     !config_list(:)%rp(2) = config_list(:)%rp(2) - tempy
     !config_list(:)%rp(3) = config_list(:)%rp(3) - tempz
     new_rp(1:3,1:nfrag_atoms) = new_rp(1:3,1:nfrag_atoms) - SPREAD(config_list_rp(1:3,clrp_anchor_frag_connect),2,nfrag_atoms)

     new_rp(1:3,1:nfrag_atoms) = MATMUL(aligner_ifrag,new_rp(1:3,1:nfrag_atoms))
         
     !DO j = 1, frag_list(ifrag,is)%natoms
     !   
     !   this_atom = frag_list(ifrag,is)%atoms(j)
     !   
     !   tempx = config_list(this_atom)%rp(1)
     !   tempy = config_list(this_atom)%rp(2)
     !   tempz = config_list(this_atom)%rp(3)
     !   
     !   config_list(this_atom)%rp(1) = tempx * aligner_ifrag(1,1) &
     !                              + tempy * aligner_ifrag(1,2) &
     !                              + tempz * aligner_ifrag(1,3)
     !   config_list(this_atom)%rp(2) = tempx * aligner_ifrag(2,1) &
     !                              + tempy * aligner_ifrag(2,2) &
     !                              + tempz * aligner_ifrag(2,3)
     !   config_list(this_atom)%rp(3) = tempx * aligner_ifrag(3,1) &
     !                              + tempy * aligner_ifrag(3,2) &
     !                              + tempz * aligner_ifrag(3,3)
     !END DO
     
     !**************************************************************************
     ! Step 3) Select kappa_dih orientations for each additional fragment
     !**************************************************************************
     ! 
     ! At this point, we can generate kappa_dih positions of the fragment as two 
     ! anchor positions are aligned. This is, in effect, equivalent to rotating 
     ! the non-anchor atoms around the x-axis. Note that the coordinates of the 
     ! anchoring atoms do not change due to this rotation. For the deletion 
     ! move, the first trial must be the one corresponding to the actual 
     ! coordinates, hence there should be no rotation about x-axis.
     theta_step_sp = species_list(is)%theta_step_sp
     theta_step = species_list(is)%theta_step

     IF ( del_flag ) THEN

        !theta = 0.0_DP
        theta0 = 0.0_DP !-threehalfPI_DP

        ! also note that we will transform the position based on hanger_ifrag so
        ! that the original positions are recovered

        hanger_frag_connect(:,:) = hanger_ifrag(:,:)

     ELSE

        ! Select a random theta with uniform probability
        
        !theta = twopi * rranf()
        theta0 = theta_step * rranf() !- threehalfPI_DP
     END IF

     ! Now that we have a starting theta, the other dihedral positions are 
     ! uniformly spaced around the 2pi disc. We need to calculate the atomic
     ! coordinates:

     ! RS: Note that both the original code and my refactor of it actually rotate the dihdral
     !     in the negative phi direction. For the dihedrals between frag_connect and ifrag,
     !     phi = phi0 - theta.  This never mattered before, but it's important when calculating dihedral phi
     !     for a given trial based on the phi and theta of another trial.

     ! Vectorized version start
     IF (.NOT. (cbmc_cell_list_flag .AND. l_widom_cells)) THEN
             inv_nfrag_atoms_dp = 1.0_DP / nfrag_atoms
             ! ignore actual mass for fragment growth.  Use average position of new atoms instead
             this_molecule%rcom(1:3) = SUM(new_rp(1:3,1:nfrag_atoms),2)*inv_nfrag_atoms_dp
             max_dcomsq = 0.0_DP
             DO j = 1, nfrag_atoms
                vec2 = new_rp(1:3,j) - this_molecule%rcom(1:3)
                max_dcomsq = MAX(max_dcomsq,DOT_PRODUCT(vec2,vec2))
             END DO
             this_molecule%rcom(4) = SQRT(max_dcomsq) ! constant for all trials
     END IF
     hanger_frag_connect_sp = REAL(hanger_frag_connect,SP)
     IF (l_get_bitcell) THEN
             zbcdf = box_list(this_box)%bitcell_dimfactor(3)
             zlbc = box_list(this_box)%length_bitcells(3)
             ybcdf = box_list(this_box)%bitcell_dimfactor(2)
             ylbc = box_list(this_box)%length_bitcells(2)
             xlbc = box_list(this_box)%length_bitcells(1)
     END IF
     ! Use Ptolemy's identities to compute SIN and COS of randomly shifted angles from precomputed array and SIN and COS of the shift
     IF (del_flag) THEN
             sintheta_sp_vec = species_list(is)%sincos_lintheta_sp(:,1)
             costheta_sp_vec = species_list(is)%sincos_lintheta_sp(:,2)
     ELSE
             costheta0 = COS(theta0)
             sintheta0 = SIN(theta0)
             costheta0_sp = REAL(costheta0,SP)
             sintheta0_sp = REAL(sintheta0,SP)
             !DIR$ VECTOR ALIGNED
             !$OMP SIMD PRIVATE(sintheta_sp,costheta_sp)
             DO itrial = 1, kappa_dih_pad8
                sintheta_sp = species_list(is)%sincos_lintheta_sp(itrial,1)
                costheta_sp = species_list(is)%sincos_lintheta_sp(itrial,2)
                sintheta_sp_vec(itrial) = sintheta_sp*costheta0_sp + costheta_sp*sintheta0_sp
                costheta_sp_vec(itrial) = costheta_sp*costheta0_sp - sintheta_sp*sintheta0_sp
             END DO
             !$OMP END SIMD
     END IF

     !! First, generate COS(theta) and SIN(theta) values for each trial
     !! The way it's done here is faster than calling both intrinsic trig functions (proven in tests)
     !! It is also faster to call SIN than COS, for some reason, so this loop calls SIN and computes COS from the result.
     !! This method only works because of the very specific range (whole period starting at -(3/2)*PI) in which theta values
     !!  are generated; otherwise the sign of COS can be wrong.  This method is still valid because it can generate theta values
     !!  anywhere in an interval that spans a whole 2*PI period.
     !!DIR$ VECTOR ALIGNED
     !DO ii = 0, kappa_dih-1
     !   theta_sp = REAL(ii,SP)
     !   theta_sp = theta0_sp + theta_step_sp*theta_sp
     !   sintheta_sp = SIN(theta_sp)
     !   sintheta_sp_vec(ii+1) = sintheta_sp
     !   costheta_sp = sintheta_sp
     !   theta_sp = theta_sp+halfPI_SP ! angle isn't actually changing, this just repurposes theta to get the sign of costheta
     !   costheta_sp = 1.0_SP - costheta_sp*costheta_sp ! costheta on RHS is actually sintheta
     !   costheta_sp = SIGN(SQRT(costheta_sp),theta_sp)
     !   costheta_sp_vec(ii+1) = costheta_sp
     !   !sintheta = SIN(theta)
     !   !costheta = COS(theta)
     !   !sinthetavec(ii) = sintheta
     !   !costhetavec(ii) = costheta
     !END DO
     rp00 = REAL(these_atoms(anchor_frag_connect)%rp(1:3),SP)
     IF (l_ortho .AND. l_widom_cells) THEN
             DO i_dim = 1,3
                inv_l_vec(i_dim) = REAL(1.0_DP/box_list(this_box)%length(i_dim,i_dim),SP)
                l_sp_vec(i_dim) = REAL(box_list(this_box)%length(i_dim,i_dim),SP)
             END DO
     ELSE IF (l_widom_cells) THEN
             inv_H_sp = REAL(box_list(this_box)%length_inv,SP)
             H_sp = REAL(box_list(this_box)%length,SP)
     END IF
     !DIR$ ASSUME (MOD(kappa_dih_pad8,8) == 0)
     !DIR$ LOOP COUNT = 3
     DO j = 1, nfrag_atoms
             x0 = REAL(new_rp(1,j),SP)
             y0 = REAL(new_rp(2,j),SP)
             z0 = REAL(new_rp(3,j),SP)
             !DIR$ ASSUME (MOD(kappa_dih_pad8,8) == 0)
             !DIR$ LOOP COUNT = 8, 16, 32
             !DIR$ VECTOR ALIGNED
             DO ii = 1, kappa_dih_pad8
                costheta_sp = costheta_sp_vec(ii)
                sintheta_sp = sintheta_sp_vec(ii)
                trial_atom_rp(ii,2,j) = costheta_sp*y0 + sintheta_sp*z0
                trial_atom_rp(ii,3,j) = costheta_sp*z0 - sintheta_sp*y0
             END DO
             rp0 = rp00
             rp0 = rp0 + x0*hanger_frag_connect_sp(1:3,1)
             !DIR$ LOOP COUNT = 8, 16, 32
             !DIR$ VECTOR ALIGNED
             DO ii = 1, kappa_dih_pad8
                tempy_sp = trial_atom_rp(ii,2,j)
                rxp = rp0(1) + tempy_sp*hanger_frag_connect_sp(1,2)
                ryp = rp0(2) + tempy_sp*hanger_frag_connect_sp(2,2)
                rzp = rp0(3) + tempy_sp*hanger_frag_connect_sp(3,2)
                tempz_sp = trial_atom_rp(ii,3,j)
                rxp = rxp + tempz_sp*hanger_frag_connect_sp(1,3)
                ryp = ryp + tempz_sp*hanger_frag_connect_sp(2,3)
                rzp = rzp + tempz_sp*hanger_frag_connect_sp(3,3)
                trial_atom_rp(ii,1,j) = rxp
                trial_atom_rp(ii,2,j) = ryp
                trial_atom_rp(ii,3,j) = rzp
             END DO
             IF (.NOT. l_widom_cells) CYCLE ! The rest of the loop body is relevant only if doing Widom insertions with cell list
             IF (l_ortho) THEN
                     DO i_dim = 1,3
                        inv_l = inv_l_vec(i_dim)
                        l_sp = l_sp_vec(i_dim)
                        hl_sp = 0.5 * l_sp
                        !DIR$ LOOP COUNT = 8, 16, 32
                        !DIR$ VECTOR ALIGNED
                        DO ii = 1, kappa_dih_pad8
                                irsp = REAL(trial_atom_rp(ii,i_dim,j),SP)
                                irsp = irsp*inv_l + 0.5
                                ! Note: it is important that the < checks and shifts come before the >= checks and shifts
                                !       due to floating point rounding.  It's fine for irsp to be zero or even
                                !       just barely less than zero due to the behavior of INT(), but if it ends up equal to
                                !       1.0, that can cause bitcell_int to be out of bounds or otherwise severely incorrect.
                                !       This was found out the hard way.
                                IF (irsp < 0.0) irsp = irsp + 1.0
                                IF (irsp >= 1.0) irsp = irsp - 1.0
                                trial_atom_rsp(ii,i_dim,j) = irsp
                                !irsp = irsp*l_sp - hl_sp
                                !sp_config_temp_list_rp(ii,this_atom,i_dim) = irsp
                        END DO
                     END DO
             ELSE
                     !DIR$ LOOP COUNT = 8, 16, 32
                     !DIR$ VECTOR ALIGNED
                     DO ii = 1, kappa_dih_pad8
                        rxp = REAL(trial_atom_rp(ii,1,j),SP)
                        ! Adding extra 0.5 because the fractional coordinates we want here assume the box is centered at 0.5
                        ! and the real coordinates are for a box centered at 0.0
                        ! Fractional coordinates elsewhere in Cassandra typically use the same centering as the real coordinates
                        ! (centered at 0.0), but it isn't convenient for bitcell overlap detection, so it's done differently here.
                        sxp = 0.5 + inv_H_sp(1,1)*rxp
                        syp = 0.5 + inv_H_sp(2,1)*rxp
                        szp = 0.5 + inv_H_sp(3,1)*rxp
                        ryp = REAL(trial_atom_rp(ii,2,j),SP)
                        sxp = sxp + inv_H_sp(1,2)*ryp
                        syp = syp + inv_H_sp(2,2)*ryp
                        szp = szp + inv_H_sp(3,2)*ryp
                        rzp = REAL(trial_atom_rp(ii,3,j),SP)
                        sxp = sxp + inv_H_sp(1,3)*rzp
                        syp = syp + inv_H_sp(2,3)*rzp
                        szp = szp + inv_H_sp(3,3)*rzp
                        IF (sxp < 0.0) sxp = sxp + 1.0
                        IF (sxp >= 1.0) sxp = sxp - 1.0
                        IF (syp < 0.0) syp = syp + 1.0
                        IF (syp >= 1.0) syp = syp - 1.0
                        IF (szp < 0.0) szp = szp + 1.0
                        IF (szp >= 1.0) szp = szp - 1.0
                        trial_atom_rsp(ii,1,j) = sxp
                        trial_atom_rsp(ii,2,j) = syp
                        trial_atom_rsp(ii,3,j) = szp
                     END DO
                     !trial_atom_rsp(1:kappa_dih,:) = &
                     !        MATMUL(sp_config_temp_list_rp(1:kappa_dih,:),TRANSPOSE(REAL(box_list(this_box)%length_inv,SP)))
             END IF
             IF (l_get_bitcell) THEN
                     !DIR$ LOOP COUNT = 8, 16, 32
                     !DIR$ VECTOR ALIGNED
                     DO ii = 1, kappa_dih
                        sxp = trial_atom_rsp(ii,1,j)
                        syp = trial_atom_rsp(ii,2,j)
                        szp = trial_atom_rsp(ii,3,j)
                        bitcell_bit = INT(sxp*xlbc)
                        bitcell_int = ISHFT(bitcell_bit,-5) + &
                                INT(syp*ylbc)*ybcdf + &
                                INT(szp*zlbc)*zbcdf
                        bitcell_bit = IAND(bitcell_bit,MASKR(5,INT32)) ! same as modulo 32
                        bitcell_int = box_list(this_box)%bitcell_int32_vec(bitcell_int)
                        bitcell_overlap(ii,j) = BTEST(bitcell_int,bitcell_bit)
                     END DO
             END IF
     END DO
     IF (l_get_bitcell) THEN
             bitcell_overlap(:,1) = .NOT. ANY(bitcell_overlap(:,1:nfrag_atoms),2)
             n_good_trials = 0
             !DIR$ LOOP COUNT = 8, 16, 32
             DO itrial = 1, kappa_dih
                IF (bitcell_overlap(itrial,1)) THEN
                        n_good_trials = n_good_trials+1
                        which_good_trials(n_good_trials) = itrial
                        trial_atom_rsp(n_good_trials,:,1:nfrag_atoms) = trial_atom_rsp(itrial,:,1:nfrag_atoms)
                END IF
             END DO
             IF (n_good_trials == 0) THEN
                cbmc_overlap = .TRUE.
                RETURN
             END IF
     ELSE
             n_good_trials = kappa_dih
     END IF
     n_good_trials_pad8 = IAND(n_good_trials+7,NOT(7))
     !DIR$ ASSUME (MOD(n_good_trials_pad8,8) == 0)
     IF (l_widom_cells) THEN
             DO j = 1, nfrag_atoms
                     IF (l_ortho) THEN
                             DO i_dim = 1, 3
                                     rlc = box_list(this_box)%real_length_cells(i_dim)
                                     l_sp = l_sp_vec(i_dim)
                                     hl_sp = 0.5 * l_sp
                                     !DIR$ LOOP COUNT = 8, 16, 32
                                     !DIR$ VECTOR ALIGNED
                                     DO itrial = 1, n_good_trials_pad8
                                        irsp = trial_atom_rsp(itrial,i_dim,j)
                                        trial_cell_coords(itrial,i_dim,j) = &
                                                INT(irsp*rlc) - box_list(this_box)%sectorbound(i_dim)
                                        irsp = irsp*l_sp - hl_sp
                                        trial_atom_rsp(itrial,i_dim,j) = irsp
                                     END DO
                             END DO
                     ELSE
                             !DIR$ LOOP COUNT = 8, 16, 32
                             !DIR$ VECTOR ALIGNED
                             DO itrial = 1, n_good_trials_pad8
                                sxp = trial_atom_rsp(itrial,1,j)
                                trial_cell_coords(itrial,1,j) = &
                                        INT(sxp*box_list(this_box)%real_length_cells(1)) - &
                                        box_list(this_box)%sectorbound(1)
                                sxp = sxp - 0.5
                                rxp = H_sp(1,1)*sxp
                                ryp = H_sp(2,1)*sxp
                                rzp = H_sp(3,1)*sxp
                                syp = trial_atom_rsp(itrial,2,j)
                                trial_cell_coords(itrial,2,j) = &
                                        INT(syp*box_list(this_box)%real_length_cells(2)) - &
                                        box_list(this_box)%sectorbound(2)
                                syp = syp - 0.5
                                rxp = rxp + H_sp(1,2)*syp
                                ryp = ryp + H_sp(2,2)*syp
                                rzp = rzp + H_sp(3,2)*syp
                                szp = trial_atom_rsp(itrial,3,j)
                                trial_cell_coords(itrial,3,j) = &
                                        INT(szp*box_list(this_box)%real_length_cells(3)) - &
                                        box_list(this_box)%sectorbound(3)
                                szp = szp - 0.5
                                rxp = rxp + H_sp(1,3)*szp
                                ryp = ryp + H_sp(2,3)*szp
                                rzp = rzp + H_sp(3,3)*szp
                                trial_atom_rsp(itrial,1,j) = rxp
                                trial_atom_rsp(itrial,2,j) = ryp
                                trial_atom_rsp(itrial,3,j) = rzp
                             END DO
                     END IF
             END DO
     END IF
     n_good_trials_old = n_good_trials
     n_good_trials = 0
     !DIR$ LOOP COUNT = 1, 2, 3, 4, 8, 16, 32
     trial_loop: DO itrial = 1, n_good_trials_old
        IF (l_widom_cells) THEN
                !DIR$ LOOP COUNT = 3
                DO j = 1, nfrag_atoms
                        ia = atom_id(j)
                        this_atom_rp = trial_atom_rsp(itrial,:,j)
                        this_atom_ci = trial_cell_coords(itrial,:,j)
                        IF (check_overlap(this_atom_rp(1:3),this_atom_ci,ia,is,this_box)) THEN
                                CYCLE trial_loop
                        END IF
                END DO
        END IF
        itrial_orig = MERGE(which_good_trials(itrial),itrial,l_get_bitcell)
        IF (cbmc_cell_list_flag .AND. l_widom_cells) THEN
                nrg_sp = SUM(Compute_Cell_List_CBMC_nrg(&
                        trial_atom_rsp(itrial,1,1:nfrag_atoms),trial_atom_rsp(itrial,2,1:nfrag_atoms),& 
                        trial_atom_rsp(itrial,3,1:nfrag_atoms), &
                        trial_cell_coords(itrial,1,1:nfrag_atoms),trial_cell_coords(itrial,2,1:nfrag_atoms), &
                        trial_cell_coords(itrial,3,1:nfrag_atoms), &
                        atom_id(1:nfrag_atoms),is,this_box))
        ELSE
                this_molecule%rcom(1:3) = 0.0_DP
                !DIR$ LOOP COUNT = 3
                DO j = 1, nfrag_atoms
                        ia = atom_id(j)
                        these_atoms(ia)%rp(1:3) = REAL(trial_atom_rp(itrial_orig,1:3,j),DP)
                        this_molecule%rcom(1:3) = this_molecule%rcom(1:3) + inv_nfrag_atoms_dp*these_atoms(ia)%rp(1:3)
                END DO
                IF (widom_active) THEN
                        CALL Compute_Molecule_Nonbond_Inter_Energy_Widom(this_im,is,nrg_vdw_dp,overlap)
                        IF (overlap) CYCLE trial_loop
                        nrg_sp = REAL(nrg_vdw_dp,SP)
                ELSE
                        CALL Compute_Molecule_Nonbond_Inter_Energy(this_im,is,nrg_vdw_dp,nrg_qq_dp,overlap)
                        IF (overlap) CYCLE trial_loop
                        nrg_sp = REAL(nrg_vdw_dp + nrg_qq_dp,SP)
                END IF
        END IF
        !!DIR$ LOOP COUNT = 3
        !DO j = 1, nfrag_atoms
        !        ia = atom_id(j)
        !        IF (cbmc_cell_list_flag .AND. l_widom_cells) THEN
        !                this_atom_rp = trial_atom_rsp(itrial,:,j)
        !                this_atom_ci = trial_cell_coords(itrial,:,j)
        !                nrg_sp = nrg_sp + Compute_Cell_List_CBMC_nrg(this_atom_rp(1:3),this_atom_ci,ia,is,this_box)
        !        ELSE
        !                these_atoms(ia)%rp(1:3) = REAL(trial_atom_rp(itrial_orig,1:3,j),DP)
        !                IF (j > 1) THEN
        !                        this_molecule%rcom(1:3) = this_molecule%rcom(1:3) + inv_nfrag_atoms_dp*these_atoms(ia)%rp(1:3)
        !                ELSE
        !                        this_molecule%rcom(1:3) = inv_nfrag_atoms_dp*these_atoms(ia)%rp(1:3)
        !                END IF
        !        END IF
        !END DO
        !IF (.NOT. widom_active) THEN
        !        CALL Compute_Molecule_Nonbond_Inter_Energy(this_im,is,nrg_vdw_dp,nrg_qq_dp,overlap)
        !        IF (overlap) CYCLE trial_loop
        !        nrg_sp = REAL(nrg_vdw_dp + nrg_qq_dp,SP)
        !ELSE IF (.NOT. cbmc_cell_list_flag) THEN
        !        CALL Compute_Molecule_Nonbond_Inter_Energy_Widom(this_im,is,nrg_vdw_dp,overlap)
        !        IF (overlap) CYCLE trial_loop
        !        nrg_sp = REAL(nrg_vdw_dp,SP)
        !END IF
        IF (nrg_sp > overlap_nrg) CYCLE trial_loop
        n_good_trials = n_good_trials + 1
        nrg_sp_vec(n_good_trials) = nrg_sp
        which_good_trials(n_good_trials) = itrial_orig
        IF (n_good_trials .NE. itrial_orig) THEN
                !DIR$ ASSUME (itrial_orig > n_good_trials)
                trial_atom_rp(n_good_trials,:,1:nfrag_atoms) = &
                        trial_atom_rp(itrial_orig,:,1:nfrag_atoms)
        END IF
        !IF (cbmc_cell_list_flag) THEN
        !        nrg_sp = 0.0
        !        DO j = 1, nfrag_atoms
        !                ia = atom_id(j)
        !                this_atom_rp = trial_atom_rsp(itrial,:,)
        !                this_atom_ci = trial_cell_coords(itrial,:,j)
        !                nrg_sp = nrg_sp + Compute_Cell_List_CBMC_nrg(this_atom_rp(1:3),this_atom_ci,ia,is,this_box)
        !        END DO
        !        nrg_sp_vec(n_good_trials) = nrg_sp
        !END IF
     END DO trial_loop
     IF (n_good_trials == 0) THEN
             cbmc_overlap = .TRUE.
             RETURN
     END IF
     DO j = 1, nfrag_atoms
        inv_atom_id(atom_id(j)) = j
     END DO
     rcutsq = REAL(rcut_cbmcsq(this_box),SP)
     anchor_atoms = (/ anchor_frag_connect, anchor_ifrag /)
     anchor_atoms_reversed = anchor_atoms(2:1:-1)
     n_good_trials_pad8 = IAND(n_good_trials+7,NOT(7))
     !DIR$ ASSUME (MOD(n_good_trials_pad8,8) .EQ. 0)
     overlap_trial(1:n_good_trials_pad8) = .FALSE.
     DO ja = 1, natoms(is)
        IF (.NOT. grown_exist(ja) .OR. ANY(ja .EQ. anchor_atoms)) CYCLE
        jrp = these_atoms(ja)%rp(1:3)
        jtype = nonbond_list(ja,is)%atom_type_number
        jcharge_dp = nonbond_list(ja,is)%charge
        !DIR$ LOOP COUNT = 3
        DO ia_new = 1, nfrag_atoms
                ia = atom_id(ia_new)
                vdwscale = vdw_intra_scale(ia,ja,is)
                chargescale = charge_intra_scale(ia,ja,is)
                cfqq = REAL(nonbond_list(ia,is)%charge * jcharge_dp * charge_factor,SP)
                itype = nonbond_list(ia,is)%atom_type_number
                this_int_vdw_style = MIN(itype*jtype,1)*int_vdw_style(this_box)
                this_int_vdw_sum_style = int_vdw_sum_style(this_box)
                ij_get_coul = cfqq .NE. 0.0 .AND. int_charge_style(this_box) .NE. charge_none .AND. &
                        species_list(is)%l_coul_cbmc
                IF ((vdwscale == 0.0 .OR. this_int_vdw_style == 0) .AND. .NOT. ij_get_coul) CYCLE
                this_int_charge_sum_style = MERGE(MERGE(charge_sf,int_charge_sum_style(this_box),cbmc_charge_sf_flag),0,ij_get_coul)
                !DIR$ LOOP COUNT = 8, 16
                !DIR$ VECTOR ALIGNED
                DO itrial = 1, n_good_trials_pad8
                        drp = trial_atom_rp(itrial,1,ia_new) - jrp(1)
                        rsq = drp*drp
                        drp = trial_atom_rp(itrial,2,ia_new) - jrp(2)
                        rsq = rsq + drp*drp
                        drp = trial_atom_rp(itrial,3,ia_new) - jrp(3)
                        rsq = rsq + drp*drp
                        overlap_trial(itrial) = overlap_trial(itrial) .OR. rsq < sp_rcut_lowsq
                        rsq_vec(itrial) = MERGE(rsq,infinity_sp,rsq<rcutsq)
                END DO
                SELECT CASE(this_int_vdw_style)
                CASE(vdw_lj)
                        SELECT CASE(this_int_vdw_sum_style)
                        CASE(vdw_charmm)
                                eps = vdwscale*ppvdwp_table2_sp(1,itype,jtype,this_box) ! epsilon
                                sigsq = ppvdwp_table2_sp(2,itype,jtype,this_box) ! sigma**2
                                !DIR$ LOOP COUNT = 8, 16
                                !DIR$ VECTOR ALIGNED
                                DO itrial = 1, n_good_trials_pad8
                                        rsq = rsq_vec(itrial)
                                        nrg_sp = nrg_sp_vec(itrial)
                                        sigbyr2 = sigsq/rsq
                                        sigbyr6 = sigbyr2*sigbyr2*sigbyr2
                                        nrg_sp = nrg_sp + eps * (sigbyr6*sigbyr6 - sigbyr6 - sigbyr6) ! eps was not multiplied by 4
                                        nrg_sp_vec(itrial) = nrg_sp
                                END DO
                        CASE(vdw_cut_shift)
                                eps = vdwscale*ppvdwp_table2_sp(1,itype,jtype,this_box) ! 4*epsilon
                                negsigsq = ppvdwp_table2_sp(2,itype,jtype,this_box) ! -(sigma**2)
                                shift_p1 = vdwscale*ppvdwp_table2_sp(3,itype,jtype,this_box)
                                !DIR$ LOOP COUNT = 8, 16
                                !DIR$ VECTOR ALIGNED
                                DO itrial = 1, n_good_trials_pad8
                                        rsq = rsq_vec(itrial)
                                        nrg_sp = nrg_sp_vec(itrial)
                                        negsigbyr2 = negsigsq/rsq
                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                        rterm = rterm + rterm*rterm
                                        nrg_sp = nrg_sp + eps*rterm - MERGE(shift_p1,0.0,rsq<rcutsq)
                                        nrg_sp_vec(itrial) = nrg_sp
                                END DO
                        CASE(vdw_cut_switch)
                                eps = vdwscale*ppvdwp_table2_sp(1,itype,jtype,this_box) ! 4*epsilon
                                negsigsq = ppvdwp_table2_sp(2,itype,jtype,this_box) ! -(sigma**2)
                                !DIR$ LOOP COUNT = 8, 16
                                !DIR$ VECTOR ALIGNED
                                DO itrial = 1, n_good_trials_pad8
                                        rsq = rsq_vec(itrial)
                                        nrg_sp = nrg_sp_vec(itrial)
                                        negsigbyr2 = negsigsq/rsq
                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                        rterm = rterm + rterm*rterm
                                        IF (rsq .GE. const3) THEN
                                                roffsq_rijsq = const4 - rsq
                                                rterm = roffsq_rijsq*roffsq_rijsq * &
                                                        (const1+const2*rsq) * rterm
                                        END IF
                                        nrg_sp = nrg_sp + eps*rterm
                                        nrg_sp_vec(itrial) = nrg_sp
                                END DO
                        CASE(vdw_cut_shift_force)
                                eps = vdwscale*ppvdwp_table2_sp(1,itype,jtype,this_box) ! 4*epsilon
                                negsigsq = ppvdwp_table2_sp(2,itype,jtype,this_box) ! -(sigma**2)
                                shift_p1 = vdwscale*ppvdwp_table2_sp(3,itype,jtype,this_box)
                                shift_p2 = vdwscale*ppvdwp_table2_sp(4,itype,jtype,this_box)
                                shift_p1 = shift_p1 + REAL(rcut_vdw(this_box),SP)*shift_p2
                                !DIR$ LOOP COUNT = 8, 16
                                !DIR$ VECTOR ALIGNED
                                DO itrial = 1, n_good_trials_pad8
                                        rsq = rsq_vec(itrial)
                                        nrg_sp = nrg_sp_vec(itrial)
                                        negsigbyr2 = negsigsq/rsq
                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                        rterm = rterm + rterm*rterm
                                        IF (rsq < rcutsq) THEN
                                                nrg_sp = nrg_sp + eps*rterm + SQRT(rsq)*shift_p2 - shift_p1
                                        END IF
                                        nrg_sp_vec(itrial) = nrg_sp
                                END DO
                        CASE(vdw_cut,vdw_cut_tail)
                                eps = vdwscale*ppvdwp_table2_sp(1,itype,jtype,this_box) ! 4*epsilon
                                negsigsq = ppvdwp_table2_sp(2,itype,jtype,this_box) ! -(sigma**2)
                                !DIR$ LOOP COUNT = 8, 16
                                !DIR$ VECTOR ALIGNED
                                DO itrial = 1, n_good_trials_pad8
                                        rsq = rsq_vec(itrial)
                                        nrg_sp = nrg_sp_vec(itrial)
                                        negsigbyr2 = negsigsq/rsq
                                        rterm = negsigbyr2*negsigbyr2*negsigbyr2
                                        rterm = rterm + rterm*rterm
                                        nrg_sp = nrg_sp + eps*rterm
                                        nrg_sp_vec(itrial) = nrg_sp
                                END DO
                        END SELECT
                CASE(vdw_mie)
                        epsig_n = vdwscale*ppvdwp_table2_sp(1,jtype,itype,this_box)
                        epsig_m = vdwscale*ppvdwp_table2_sp(2,jtype,itype,this_box)
                        mie_n = ppvdwp_table2_sp(3,jtype,itype,this_box)
                        mie_m = ppvdwp_table2_sp(4,jtype,itype,this_box)
                        SELECT CASE(this_int_vdw_sum_style)
                        CASE(vdw_cut_shift)
                                shift_p1 = vdwscale*ppvdwp_table2_sp(5,itype,jtype,this_box)
                                !DIR$ LOOP COUNT = 8, 16
                                !DIR$ VECTOR ALIGNED
                                DO itrial = 1, n_good_trials_pad8
                                        rsq = rsq_vec(itrial)
                                        nrg_sp = nrg_sp_vec(itrial)
                                        lnrsq = LOG(rsq)
                                        nrg_sp = nrg_sp + epsig_n*EXP(lnrsq*mie_n) - epsig_m*EXP(lnrsq*mie_m) - shift_p1
                                        nrg_sp_vec(itrial) = nrg_sp
                                END DO
                        CASE(vdw_cut,vdw_cut_tail)
                                !DIR$ LOOP COUNT = 8, 16
                                !DIR$ VECTOR ALIGNED
                                DO itrial = 1, n_good_trials_pad8
                                        rsq = rsq_vec(itrial)
                                        nrg_sp = nrg_sp_vec(itrial)
                                        lnrsq = LOG(rsq)
                                        nrg_sp = nrg_sp + epsig_n*EXP(lnrsq*mie_n) - epsig_m*EXP(lnrsq*mie_m)
                                        nrg_sp_vec(itrial) = nrg_sp
                                END DO
                        END SELECT
                END SELECT
                SELECT CASE(this_int_charge_sum_style)
                CASE(charge_ewald)
                        alpha_ewald_sp = REAL(alpha_ewald(this_box),SP)
                        !DIR$ LOOP COUNT = 8, 16
                        !DIR$ VECTOR ALIGNED
                        DO itrial = 1, n_good_trials_pad8
                                rsq = rsq_vec(itrial)
                                nrg_sp = nrg_sp_vec(itrial)
                                invr = 1.0/SQRT(rsq)
                                rij = rsq*invr
                                nrg_sp = nrg_sp + (chargescale - ERF(alpha_ewald_sp*rij))*cfqq*invr
                                nrg_sp_vec(itrial) = nrg_sp
                        END DO
                CASE(charge_dsf)
                        dsf_const = -dsf_factor2(this_box)*(rcut_coul(this_box)+dsf_factor1(this_box))
                        dsf_factor2_sp = REAL(dsf_factor2(this_box),SP)
                        !DIR$ LOOP COUNT = 8, 16
                        !DIR$ VECTOR ALIGNED
                        DO itrial = 1, n_good_trials_pad8
                                rsq = rsq_vec(itrial)
                                nrg_sp = nrg_sp_vec(itrial)
                                invr = 1.0/SQRT(rsq)
                                rij = rsq*invr
                                nrg_sp = dsf_const + (chargescale - ERF(alpha_dsf(this_box)*rij))*invr + dsf_factor2_sp*rij
                                nrg_sp = cfqq*nrg_sp + nrg_sp_vec(itrial)
                                nrg_sp_vec(itrial) = nrg_sp
                        END DO
                CASE(charge_cut)
                        cfqq = cfqq*chargescale
                        !DIR$ LOOP COUNT = 8, 16
                        !DIR$ VECTOR ALIGNED
                        DO itrial = 1, n_good_trials_pad8
                                rsq = rsq_vec(itrial)
                                nrg_sp = nrg_sp_vec(itrial)
                                invr = 1.0/SQRT(rsq)
                                nrg_sp = nrg_sp + cfqq*invr
                                nrg_sp_vec(itrial) = nrg_sp
                        END DO
                CASE(charge_sf)
                        sf_const1 = REAL(2.0_DP/rcut_cbmc(this_box),SP)
                        sf_const2 = 1.0/rcutsq
                        !DIR$ LOOP COUNT = 8, 16
                        !DIR$ VECTOR ALIGNED
                        DO itrial = 1, n_good_trials_pad8
                                rsq = rsq_vec(itrial)
                                nrg_sp = nrg_sp_vec(itrial)
                                invr = 1.0/SQRT(rsq)
                                rij = rsq*invr
                                rij = rij*sf_const2 - sf_const1 ! rij no longer storing atomic separation distance
                                rij = MERGE(rij,0.0,rsq<rcutsq)
                                rij = rij + chargescale*invr
                                nrg_sp = nrg_sp + cfqq*rij
                                nrg_sp_vec(itrial) = nrg_sp
                        END DO
                END SELECT
        END DO
     END DO
     ! Having intramolecular overlap is unlikely in most cases
     IF (ANY(overlap_trial(1:n_good_trials))) THEN
             IF (ALL(overlap_trial(1:n_good_trials))) THEN
                     cbmc_overlap = .TRUE.
                     RETURN
             END IF
             n_good_trials_old = n_good_trials
             n_good_trials = 0
             DO itrial = 1, n_good_trials_old
                IF (.NOT. overlap_trial(itrial)) THEN
                        n_good_trials = n_good_trials + 1
                        which_good_trials(n_good_trials) = MERGE(itrial,which_good_trials(itrial),n_good_trials_old==kappa_dih)
                        trial_atom_rp(n_good_trials,1:3,1:nfrag_atoms) = &
                                trial_atom_rp(itrial,1:3,1:nfrag_atoms)
                        nrg_sp_vec(n_good_trials) = nrg_sp_vec(itrial)
                END IF
             END DO
             n_good_trials_pad8 = IAND(n_good_trials+7,NOT(7))
     END IF
     !DIR$ ASSUME (MOD(n_good_trials_pad8,8) == 0)
     nrg_dihed_vec = 0.0
     anchor_frag_connect_rp = REAL(these_atoms(anchor_frag_connect)%rp(1:3),SP)
     anchor_ifrag_rp = REAL(these_atoms(anchor_ifrag)%rp(1:3),SP)
     !DIR$ LOOP COUNT = 0, 9, 18
     DO idihed_rb = 1, species_list(is)%ndihedrals_rb
        IF (l_skip_dihed_vec(idihed_rb)) CYCLE ! it would CYCLE anyway due to other conditionals but this is faster
        dihedral_atoms = dihedral_list(idihed_rb,is)%atom(1:4)
        IF (ALL(dihedral_atoms(2:3) .EQ. anchor_atoms) .AND. &
                grown_exist(dihedral_atoms(1)) .AND. &
                new_exist(dihedral_atoms(4))) THEN
                atom1 = dihedral_atoms(1)
                atom2 = dihedral_atoms(2)
                atom3 = dihedral_atoms(3)
                atom4 = inv_atom_id(dihedral_atoms(4))
                r32_dp = r32_base
                atom3_rp = anchor_ifrag_rp
        ELSE IF (ALL(dihedral_atoms(2:3) .EQ. anchor_atoms_reversed) .AND. &
                new_exist(dihedral_atoms(1)) .AND. &
                grown_exist(dihedral_atoms(4))) THEN
                atom1 = dihedral_atoms(4)
                atom2 = dihedral_atoms(3)
                atom3 = dihedral_atoms(2)
                atom4 = inv_atom_id(dihedral_atoms(1))
                r32_dp = -r32_base
                atom3_rp = anchor_frag_connect_rp
        ELSE
                CYCLE
        END IF
        l_skip_dihed_vec(idihed_rb) = .TRUE.
        vecdp = these_atoms(atom1)%rp(1:3) - &
                these_atoms(atom2)%rp(1:3) ! r12
        m_dp(1) = vecdp(2)*r32_dp(3) - vecdp(3)*r32_dp(2)
        m_dp(2) = vecdp(3)*r32_dp(1) - vecdp(1)*r32_dp(3)
        m_dp(3) = vecdp(1)*r32_dp(2) - vecdp(2)*r32_dp(1)
        m_normsq = REAL(DOT_PRODUCT(m_dp,m_dp),SP)
        m_sp = REAL(m_dp,SP)
        rx32 = REAL(r32_dp(1),SP)
        ry32 = REAL(r32_dp(2),SP)
        rz32 = REAL(r32_dp(3),SP)
        !DIR$ LOOP COUNT = 8, 16
        !DIR$ VECTOR ALIGNED
        DO itrial = 1, n_good_trials_pad8
                rx34 = trial_atom_rp(itrial,1,atom4)
                ry34 = trial_atom_rp(itrial,2,atom4)
                rz34 = trial_atom_rp(itrial,3,atom4)
                rx34 = atom3_rp(1) - rx34
                ry34 = atom3_rp(2) - ry34
                rz34 = atom3_rp(3) - rz34
                ni = ry32*rz34 - rz32*ry34 ! nx
                cosphi = m_sp(1)*ni ! currently accumulator corresponding to mdn
                invnorm = ni*ni ! invnorm similarly needs to undergo further operations
                ni = rz32*rx34 - rx32*rz34 ! ny
                cosphi = cosphi + m_sp(2)*ni
                invnorm = invnorm + ni*ni
                ni = rx32*ry34 - ry32*rx34 ! nz
                invnorm = invnorm + ni*ni
                cosphi = cosphi + m_sp(3)*ni
                invnorm = invnorm * m_normsq
                invnorm = 1.0_SP / SQRT(invnorm) ! finally reciprocal of norm product
                cosphi = cosphi*invnorm ! now cosphi is actually COS(phi)
                cosphi_vec(itrial) = cosphi
        END DO
        rb_c(1:5) = dihedral_list(idihed_rb,is)%rb_c_sp(1:5)
        !DIR$ LOOP COUNT = 8, 16
        !DIR$ VECTOR ALIGNED
        DO itrial = 1, n_good_trials_pad8
                cosphi = cosphi_vec(itrial)
                cosphi2 = cosphi*cosphi
                cosphi3 = cosphi*cosphi2
                cosphi4 = cosphi2*cosphi2
                cosphi5 = cosphi2*cosphi3
                nrg_dihed_vec(itrial) = nrg_dihed_vec(itrial) + &
                        rb_c(1)*cosphi + &
                        rb_c(2)*cosphi2 + &
                        rb_c(3)*cosphi3 + &
                        rb_c(4)*cosphi4 + &
                        rb_c(5)*cosphi5
        END DO
     END DO
     !DIR$ LOOP COUNT = 0, 9, 18
     DO idihed = idihed_rb, species_list(is)%ndihedrals_energetic
        IF (l_skip_dihed_vec(idihed)) CYCLE ! it would CYCLE anyway due to other conditionals but this is faster
        dihedral_atoms = dihedral_list(idihed_rb,is)%atom(1:4)
        IF (ALL(dihedral_atoms(2:3) .EQ. anchor_atoms) .AND. &
                grown_exist(dihedral_atoms(1)) .AND. &
                new_exist(dihedral_atoms(4))) THEN
                these_atoms(dihedral_atoms(4))%rp(1:3) = &
                        trial_atom_rp(n_good_trials,1:3,inv_atom_id(dihedral_atoms(4)))
        ELSE IF (ALL(dihedral_atoms(2:3) .EQ. anchor_atoms_reversed) .AND. &
                new_exist(dihedral_atoms(1)) .AND. &
                grown_exist(dihedral_atoms(4))) THEN
                these_atoms(dihedral_atoms(1))%rp(1:3) = &
                        trial_atom_rp(n_good_trials,1:3,inv_atom_id(dihedral_atoms(1)))
        ELSE
                CYCLE
        END IF
        l_skip_dihed_vec(idihed) = .TRUE.
        CALL Get_Dihedral_Angle(idihed,this_im,is,last_phi_dp)
        SELECT CASE(dihedral_list(idihed,is)%int_dipot_type)
        CASE(int_charmm)
                last_phi = REAL(last_phi_dp - dihedral_list(idihed,is)%dihedral_param(3),SP)
                phi_step_sp = REAL(theta_step * dihedral_list(idihed,is)%dihedral_param(2),SP)
        CASE(int_harmonic)
                last_phi = REAL(last_phi_dp - dihedral_list(idihed,is)%dihedral_param(2) - twoPI,SP)
                phi_step_sp = theta_step_sp
        END SELECT
        a0 = dihedral_list(idihed,is)%dihedral_param_sp(1)
        IF (l_widom_cells) THEN
                DO itrial = 1, n_good_trials-1
                        phi_vec(itrial) = (which_good_trials(n_good_trials)-which_good_trials(itrial)) * &
                                phi_step_sp + last_phi
                END DO
        ELSE
                DO itrial = 1, n_good_trials-1
                        phi_vec(itrial) = (n_good_trials-itrial) * &
                                phi_step_sp + last_phi
                END DO
        END IF
        ! I ignore any constant terms in dihedral energy because only relative energy matters.
        SELECT CASE(dihedral_list(idihed,is)%int_dipot_type)
        CASE(int_charmm)
                nrg_dihed_vec(1:n_good_trials) = nrg_dihed_vec(1:n_good_trials) + &
                        a0*COS(phi_vec(1:n_good_trials))
        CASE(int_harmonic)
                DO itrial = 1, n_good_trials
                        phi = ABS(phi_vec(itrial))
                        IF (phi > PI_SP) phi = phi - twoPI_SP
                        nrg_dihed_vec(itrial) = nrg_dihed_vec(itrial) + phi*phi*a0
                END DO
        END SELECT
     END DO
     ! Compute trial weights
     !DIR$ LOOP COUNT = 8, 16
     !DIR$ VECTOR ALIGNED
     DO itrial = 1, n_good_trials_pad8
        nrg_sp = nrg_sp_vec(itrial) + nrg_dihed_vec(itrial)
        nrg_dp = REAL(nrg_sp,DP)
        nrg_dp_vec(itrial) = nrg_dp
        weight(itrial) = EXP(-beta(this_box)*nrg_dp)
     END DO
     cweight = weight(1)
     !DIR$ LOOP COUNT = 1, 2, 3, 4, 8
     DO itrial = 1, n_good_trials
        cweight = cweight + weight(itrial)
        weight(itrial) = cweight
     END DO
     ! If the cumulative weight is 0, then all trial dihedrals had core overlap
     ! Reject the move if all trials tripped overlap
     IF (cweight == 0.0_DP ) THEN
             cbmc_overlap = .TRUE.
             RETURN
     END IF
     IF (del_flag) THEN
             trial = 1
             itrial = 1
             IF (n_good_trials .NE. kappa_dih .AND. which_good_trials(1) .NE. 1) THEN
                     cbmc_overlap = .TRUE.
                     RETURN
             END IF
     ELSE
             prob_pick = rranf() * cweight
             DO itrial = 1, n_good_trials
                IF (prob_pick <= weight(itrial)) EXIT
             END DO
             ! The above loop should always exit "early", so itrial should never be > n_good_trials
             trial = MERGE(itrial, which_good_trials(itrial), n_good_trials==kappa_dih)
     END IF
     ! place atoms for chosen trial
     sinlintheta = species_list(is)%sincos_lintheta_dp(trial,1)
     coslintheta = species_list(is)%sincos_lintheta_dp(trial,2)
     sintheta = sinlintheta*costheta0 + coslintheta*sintheta0
     costheta = coslintheta*costheta0 - sinlintheta*sintheta0
     sincosmat(1,:) = (/ costheta, sintheta /)
     sincosmat(2,:) = (/ -sintheta, costheta /)
     new_rp(2:3,1:nfrag_atoms) = MATMUL(sincosmat,new_rp(2:3,1:nfrag_atoms))
     new_rp(1:3,1:nfrag_atoms) = MATMUL(hanger_frag_connect,new_rp(1:3,1:nfrag_atoms))
     DO j = 1, nfrag_atoms
        these_atoms(atom_id(j))%rp(1:3) = new_rp(1:3,j) + these_atoms(anchor_frag_connect)%rp(1:3)
     END DO
     these_atoms%exist = combined_exist
     IF (widom_active) THEN
             ! Note, this part could be precomputed in the simulation setup instead, since the interfragment
             !      intramolecular nonbonded pairs are the same for each run of this subroutine for a given species
             !      unless you are doing a partial deletion and regrowth, in which case some of the pairs would be
             !      unneeded, but they shouldn't make the computation wrong since the extra pairs would be the same for
             !      the partial deletion and the regrowth, and the extra energies would cancel out when computing dE.
             ! I just didn't take the time to do so.
             DO ja = 1, natoms(is)
                IF ((.NOT. grown_exist(ja)) .OR. ANY(ja == anchor_atoms)) CYCLE
                interfrag_nonbond_pairlist(1,n_interfrag_nonbond_pairs+1:n_interfrag_nonbond_pairs+nfrag_atoms) = &
                        SPREAD(ja,1,nfrag_atoms)
                interfrag_nonbond_pairlist(2,n_interfrag_nonbond_pairs+1:n_interfrag_nonbond_pairs+nfrag_atoms) = &
                        atom_id(1:nfrag_atoms)
                n_interfrag_nonbond_pairs = n_interfrag_nonbond_pairs + nfrag_atoms
             END DO
     END IF
     ! Recover the individual probability for the accepted trial
     ln_pbias = ln_pbias - beta(this_box) * nrg_dp_vec(itrial) - LOG(cweight)
     ! mark this fragment as placed
     frag_placed(ifrag) = 1
     ! Vectorized version end

     ! Loop over the trial dihedrals

     !**************************************************************************
     ! Step 4) Compute the energy of the fragment 
     !**************************************************************************
     ! 

     !**************************************************************************
     ! Step 5) Select a trial dihedral using the weighted probabilities
     !**************************************************************************
     ! 
     ! 

     
     

    

  END DO
  CALL Set_CBMC_Flag(.FALSE.)
  ! apply this part to other CBMC moves, not just Widom insertions, if we determine it's correct for them
  IF (widom_active .AND. n_interfrag_nonbond_pairs > 0) THEN
          DO i = 1, n_interfrag_nonbond_pairs
                ia = interfrag_nonbond_pairlist(1,i)
                ja = interfrag_nonbond_pairlist(2,i)
                vec2 = these_atoms(ja)%rp(1:3) - these_atoms(ia)%rp(1:3)
                rijsq_dp = DOT_PRODUCT(vec2,vec2)
                get_vdw = rijsq_dp < rcut_vdwsq(this_box) .AND. ALL(nonbond_list((/ia,ja/),is)%atom_type_number .NE. 0)
                get_qq = rijsq_dp < rcut_coulsq(this_box) .AND. ALL(nonbond_list((/ia,ja/),is)%charge .NE. 0.0_DP)
                IF (.NOT. (get_vdw .OR. get_qq)) CYCLE
                CALL Compute_AtomPair_Energy(vec2(1),vec2(2),vec2(3),rijsq_dp,is,this_im,ia,is,this_im,ja, &
                        get_vdw,get_qq, Eij_intra_vdw, Eij_intra_qq, Eij_inter_vdw, Eij_inter_qq)
                E_total = E_total + Eij_intra_vdw + Eij_intra_qq + Eij_inter_qq + Eij_inter_vdw
          END DO
  END IF
  IF (widom_active .AND. ANY(l_skip_dihed_vec)) THEN
          ! Previously indicated dihedral should be skipped early for biasing energy because it was already done
          !    and there's no way it would be needed twice (at least with the present CBMC algorithms)
          l_skip_dihed_vec = .NOT. l_skip_dihed_vec ! inverting l_skip_dihed_vec
          ! l_skip_dihed_vec now indicates which dihedrals to skip when computing interfragment dihedral energy
          !    because they aren't interfragment dihedrals
          CALL Compute_Molecule_Dihedral_Energy(this_im,is,E_dihed,l_skip_dihed_vec)
          E_total = E_total + E_dihed
  END IF
  IF (frag_total > 1 .AND. .NOT. (full_cell_list_flag .AND. widom_active)) THEN
          CALL Get_COM(this_im,is)
          CALL Compute_Max_COM_Distance(this_im,is)
  END IF

 END SUBROUTINE Fragment_Placement
                 
!***************************************************************************************************
SUBROUTINE Get_Aligner_Hanger(vec1,vec2,aligner,hanger)
!***************************************************************************************************


   USE Type_Definitions, ONLY : DP
   
   IMPLICIT NONE
   
   REAL(DP) :: vec1(3), vec2(3), perp_vec1(3), perp_vec2(3)
   REAL(DP), INTENT(OUT) :: aligner(3,3), hanger(3,3)
  
   INTEGER :: i, j
   

   aligner(:,:) = 0.0_DP
   hanger(:,:) = 0.0_DP
   
   ! Normalize the vectors
   
   vec1(:) = vec1(:) / DSQRT(DOT_PRODUCT(vec1,vec1))
   vec2(:) = vec2(:) / DSQRT(DOT_PRODUCT(vec2,vec2))
   
   ! Now cross the vector to get the vector perpendicular to these vectors
   
   perp_vec1(1) =  vec1(2) * vec2(3) - vec1(3) * vec2(2)
   perp_vec1(2) = -vec1(1) * vec2(3) + vec1(3) * vec2(1)
   perp_vec1(3) =  vec1(1) * vec2(2) - vec1(2) * vec2(1)
   
   perp_vec1(:) = perp_vec1(:) / DSQRT(DOT_PRODUCT(perp_vec1,perp_vec1))
   
   
   perp_vec2(1) =  vec1(2) * perp_vec1(3) - vec1(3) * perp_vec1(2)
   perp_vec2(2) = -vec1(1) * perp_vec1(3) + vec1(3) * perp_vec1(1)
   perp_vec2(3) =  vec1(1) * perp_vec1(2) - vec1(2) * perp_vec1(1)
   
   perp_vec2(:) = perp_vec2 / DSQRT(DOT_PRODUCT(perp_vec2,perp_vec2))
   
   DO i = 1, 3
      
      aligner(1,i) = vec1(i)
      aligner(2,i) = perp_vec2(i)
      aligner(3,i) = perp_vec1(i)
      
      hanger(i,1) = vec1(i) 
      hanger(i,2) = perp_vec2(i)
      hanger(i,3) = perp_vec1(i)
      
   END DO
   
   
 END SUBROUTINE Get_Aligner_Hanger
!***********************************************************
!
!*********************************************************************
SUBROUTINE Single_Fragment_Regrowth(this_im,is)
!*********************************************************************
   ! This routine is used when a species contains only one fragment
   ! and a change in intramolecular degrees of freedom is desired. 
   !
   ! Written by Jindal K. Shah on 01/01/09
   !
   ! Algorithm
   !
   ! 1. Pick a configuration from the reservoir library
   ! 2. Insert the fragment such that the COM coincides with the fragment
   !    being "deleted"
   ! 3. Give it a random Eulerian rotation.
   !*****************************************************************

   
   USE Rotation_Routines
   
   INTEGER, INTENT(IN) :: this_im, is

   TYPE(Molecule_Class), POINTER :: this_molecule
   TYPE(Atom_Class), POINTER :: these_atoms(:)


   INTEGER :: nl, nlo      ! number of the line where start the x,y,x coords of
                           ! config and fragment randomly selected
   INTEGER :: total_frags, i, this_atom, this_fragment, frag_type

   REAL(DP) :: temp_var, E_ang,x_this, y_this, z_this

  IF (widom_active) THEN
          this_molecule => widom_molecule
          these_atoms => widom_atoms
  ELSE
          this_molecule => molecule_list(this_im,is)
          these_atoms => atom_list(1:natoms(is),this_im,is)
  END IF


   ! As this molecule has only 1 fragment, we obtain total number
   ! of fragments and atoms using 1 as the identity

   total_frags = frag_list(1,is)%nconfig
     
   ! Choose a fragment at random
   
   this_fragment = INT(rranf() * total_frags) + 1
   
   frag_type = frag_list(1,is)%type

   DO i = 1, frag_list(1,is)%natoms 
      
      this_atom = frag_list(1,is)%atoms(i)
      nl = (frag_position_library(frag_type)-1) + &
                                           frag_list(1,is)%natoms*(this_fragment -1) + i
      these_atoms(this_atom)%rp(1:3) = library_coords(1:3,nl) 
      !these_atoms(this_atom)%rp(1) =  library_coords(nl)%rp(1)
      !                                       !frag_coords(i,this_fragment,frag_type)%rp(1)
      !these_atoms(this_atom)%rp(2) =  library_coords(nl)%rp(2)
      !                                       !frag_coords(i,this_fragment,frag_type)%rp(2)
      !these_atoms(this_atom)%rp(3) =  library_coords(nl)%rp(3)
      !                                        !frag_coords(i,this_fragment,frag_type)%rp(3)
   END DO



   ! COM and max_com_distance. Note that the following calls destroy
   ! original COM and max_com_distance but we will use xcom_old etc to
   ! transfer the COM to the original position

   CALL Get_COM(this_im,is)
   CALL Compute_Max_COM_Distance(this_im,is)

   ! Give random orientation to the fragment

   CALL Rotate_Molecule_Eulerian(this_im,is)

   ! Move the atoms such that COM gets translated to original position

   DO i = 1, frag_list(1,is)%natoms
      this_atom = frag_list(1,is)%atoms(i)

      these_atoms(this_atom)%rp(1) = these_atoms(this_atom)%rp(1) + &
           this_molecule%rcom_old(1) - this_molecule%rcom(1)

      these_atoms(this_atom)%rp(2) = these_atoms(this_atom)%rp(2) + &
           this_molecule%rcom_old(2) - this_molecule%rcom(2)

      these_atoms(this_atom)%rp(3) = these_atoms(this_atom)%rp(3) + &
           this_molecule%rcom_old(3) - this_molecule%rcom(3)

   END DO
   CALL Get_COM(this_im,is)
   CALL Compute_Max_COM_Distance(this_im,is)


 END SUBROUTINE Single_Fragment_Regrowth

 !***********************************************************************************
SUBROUTINE Get_Common_Fragment_Atoms(is,frag1,frag2,atom1,atom2)
   !*********************************************************************************
   !
   ! This routine determines atoms connecting two fragments
   !
   ! written by Jindal Shah on 04/16/09
   !
   !*********************************************************************************
   
   INTEGER, INTENT(IN) :: is, frag1, frag2
   
   INTEGER, INTENT(OUT) :: atom1, atom2
   
   INTEGER :: i, j, this_atom, nanchor1, nanchor2, atom_i
   
   nanchor1 = frag_list(frag1,is)%nanchors
   nanchor2 = frag_list(frag2,is)%nanchors
   
   IF ( nanchor1 == 1 .AND. nanchor2 == 1 ) THEN
      
      atom1 = frag_list(frag1,is)%anchor(1)
      atom2 = frag_list(frag2,is)%anchor(1)
      
      
   ELSE IF ( nanchor1 == 1) THEN
      
      atom1 = frag_list(frag1,is)%anchor(1)
      
      ! for the second fragment loop over all the anchors of frag2 with atoms of frag1
      
      anchor2_loop:DO i = 1, nanchor2
        
         atom_i = frag_list(frag2,is)%anchor(i)
         
         DO j = 1, frag_list(frag1,is)%natoms
            
            this_atom = frag_list(frag1,is)%atoms(j)
            
            IF ( atom_i == this_atom) EXIT anchor2_loop
            
         END DO
         
      END DO anchor2_loop
      
      atom2 = atom_i
      
      
   ELSE IF (nanchor2 == 1 ) THEN
      
      atom2 = frag_list(frag2,is)%anchor(1)
      
      ! for the first fragment loop over all the anchors of frag1 with atoms of frag2
      
      anchor1_loop: DO i = 1, nanchor1
         
         atom_i = frag_list(frag1,is)%anchor(i)
         
         DO j = 1, frag_list(frag2,is)%natoms
            
            this_atom = frag_list(frag2,is)%atoms(j)
            
            IF ( atom_i == this_atom ) EXIT anchor1_loop
            
         END DO
         
      END DO anchor1_loop
      
      atom1 = atom_i
      
      
   ELSE
      
      ! multiple anchor atoms in both the fragments
      
      anchor1a_loop: DO i = 1, nanchor1
         
         atom_i = frag_list(frag1,is)%anchor(i)
         
         DO j = 1, frag_list(frag2,is)%natoms
            
            this_atom = frag_list(frag2,is)%atoms(j)
            
            IF ( atom_i == this_atom ) EXIT anchor1a_loop
            
         END DO
         
      END DO anchor1a_loop
      
      atom1 = atom_i
      
      
      anchor2a_loop:DO i = 1, nanchor2
         
         atom_i = frag_list(frag2,is)%anchor(i)
         
         DO j = 1, frag_list(frag1,is)%natoms
            
            this_atom = frag_list(frag1,is)%atoms(j)
            
            IF ( atom_i == this_atom) EXIT anchor2a_loop
            
         END DO
         
      END DO anchor2a_loop
      
      atom2 = atom_i
      
      
   END IF

   
 END SUBROUTINE Get_Common_Fragment_Atoms

 SUBROUTINE Set_CBMC_Flag(flag_in)
         LOGICAL, INTENT(IN) :: flag_in
         IF (widom_active) THEN
                 CALL Set_Thread_CBMC_Flag
         ELSE
                 CALL Set_All_CBMC_Flag
         END IF
         CONTAINS
                 SUBROUTINE Set_All_CBMC_Flag
                         !$OMP PARALLEL
                         cbmc_flag = flag_in
                         !$OMP END PARALLEL
                 END SUBROUTINE Set_All_CBMC_Flag

                 SUBROUTINE Set_Thread_CBMC_Flag
                         cbmc_flag = flag_in
                 END SUBROUTINE Set_Thread_CBMC_Flag
 END SUBROUTINE Set_CBMC_Flag

END MODULE Fragment_Growth
