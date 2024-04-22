
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

MODULE Sector_Routines

  !***************************************************************************
  !
  !
  !********************************************************************************
  USE Global_Variables
  USE Type_Definitions
  USE Io_Utilities
  !$ USE OMP_LIB

  IMPLICIT NONE

  LOGICAL :: l_firstframe = .FALSE.


  INTERFACE check_overlap
          MODULE PROCEDURE check_overlap_ams, check_overlap_coordinates
  END INTERFACE
  


CONTAINS
  SUBROUTINE Sector_Setup
          INTEGER, DIMENSION(3) :: adj_cellmaxbound_old
          INTEGER :: i_sector, dx, dy, dz, xshift, yshift, zshift, nsec_old, nsec, secind
          INTEGER :: sector_ID
          INTEGER :: i, ibox, is, is_present, imol, im, ia
          INTEGER, DIMENSION(3) :: sector_atom_ID
          TYPE(Atom_Class), POINTER :: atom_ptr
          REAL(DP) :: xp, yp, zp
!          INTEGER, DIMENSION(3,nbr_boxes) :: cbmc_truth_cube_bound, cut_truth_cube_bound
          INTEGER :: xi, yi, zi, cim(3), xyzi(3), i_dim
          INTEGER, DIMENSION(2,3) :: tgt_slice, src_slice, bit_tgt_slice, bit_src_slice
          INTEGER, DIMENSION(:), ALLOCATABLE :: xi_pm, yi_pm, zi_pm
          INTEGER :: dummy

          REAL(SP), DIMENSION(maxboxnatoms,4,nbr_boxes) :: sp_live_atom_rsp
          INTEGER, DIMENSION(maxboxnatoms,3,nbr_boxes) :: live_atom_cp
          INTEGER(INT32), DIMENSION(maxboxnatoms,nbr_boxes) :: live_atom_ti, live_atom_atomtypes
          INTEGER, DIMENSION(4,maxboxnatoms,nbr_boxes) :: ci_list
          INTEGER(INT32) :: bsolvent
          INTEGER :: ci(4), nca

          INTEGER :: xi2, yi2, zi2, vlen, n_i_exist, inlive, inatoms, istart, iend
          INTEGER :: cp_ub, cp_lb, lc, cp
          INTEGER :: dxi, dyi, dzi
          INTEGER, DIMENSION(nbr_boxes) :: box_vlen

          REAL(SP) :: rp_ub, rp_lb, lbox, clr, rsp, isp, rxp, ryp, rzp
          REAL(SP) :: h11,h21,h31,h12,h22,h32,h13,h23,h33
          REAL(SP), DIMENSION(3) :: boxlen, unwrap_shifter

          LOGICAL, DIMENSION(maxboxnatoms) :: i_exist
          LOGICAL :: l_ortho

          INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE :: n_cell_atoms
          INTEGER(INT32), DIMENSION(:,:,:,:,:), ALLOCATABLE :: this_cell_ti, this_cell_atomtypes
          INTEGER(INT32), DIMENSION(:,:,:,:), ALLOCATABLE :: this_cell_atom_i
          REAL(SP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: this_cell_rsp

          LOGICAL :: need_atom_ti, need_atomtypes, need_charges

          REAL(SP) :: cell_H(3,3), length_cells_recip(3), hlcr(3), hl(3)

          INTEGER(INT64), DIMENSION(-28:28,-28:28,solvents_or_types_maxind,0:n_big_atoms,nbr_boxes) :: bitcell_int64
          LOGICAL(1), DIMENSION(solvents_or_types_maxind,0:n_big_atoms) :: l_inrange_vec, l_inrange_vec_old, l_switch
          REAL(DP) :: xyzi_dp(3), dxyzi_dp(3), xyzi_dp_spread(3,3,2), drp(3), rsq, bitcell_H(3,3), bitcell_H_diag(3)
          INTEGER(INT64) :: bitmask
          INTEGER(INT64), DIMENSION(solvents_or_types_maxind,0:n_big_atoms) :: priv_bitcell_int64
          INTEGER, DIMENSION(maxboxnatoms,4,nbr_boxes) :: live_atom_bcp
          INTEGER, DIMENSION(4,maxboxnatoms) :: live_atom_bcp_T
          INTEGER :: bcp,bcpx,bcpy,bcpz,bcps,ti,lbp32(3)
          INTEGER(1), DIMENSION(:,:,:,:), ALLOCATABLE :: bitcell_int8_array
          INTEGER(1), DIMENSION(:,:,:), ALLOCATABLE :: bitcell_int8_array_2

          INTEGER :: max_adj_cell_atoms_old, cbmc_max_interact_old, max_neighbors
          INTEGER :: bt(3), sbe(3)

          REAL(SP) :: rcutsq, cell_rsxp, cell_rsyp, cell_rszp
          REAL(SP), DIMENSION(maxboxnatoms) :: rsq_vec
          REAL(SP), DIMENSION(maxboxnatoms,4) :: cbmc_cell_rsp_priv
          INTEGER, DIMENSION(maxboxnatoms) :: ti_priv, atomtype_priv, which_interact, which_i_exist

          INTEGER :: bcd(2:3)

          LOGICAL(1), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: cell_l_inrange!, adj_cell_l_inrange

          INTEGER :: isolvent, lbc, border_range(2), n_interact, adj_iend
          REAL(SP) :: bclr, dsp, drxp, dryp, drzp

          INTEGER(INT32), PARAMETER :: pad8mask = NOT(MASKR(3,INT32))
          INTEGER(INT64) :: superpopcnt

          REAL(DP), DIMENSION(3,3) :: bitcell_H_T, bitcell_H_T_sign, cell_H_T, cell_H_T_sign
          REAL(DP), DIMENSION(3,3,2) :: bitcell_H_T_spread, bitcell_H_T_posneg_sign, cell_H_T_spread, cell_H_T_posneg_sign

          REAL(DP), DIMENSION(3) :: bitcell_xyzortho_bbox_length
          INTEGER :: yi_mult, zi_mult
          REAL(DP), DIMENSION(3,3) :: cell_H_dp
          REAL(SP), DIMENSION(3) :: xyzi_sp, cell_rp
          INTEGER :: bcp_shift, xcp_base_shift_stride, xcp_base_stride, zcp_base_shift, ycp_base_shift, xcp_base_shift
          INTEGER :: zcp_base, ycp_base, xcp_base, j, xcp, ycp, zcp, xub, ylb, zlb, xlb, yub, zub
          INTEGER(INT64) :: xfer_int64
          INTEGER :: sbe_ti(2:3), sbe_ti_mat(2:3,solvents_or_types_maxind,0:n_big_atoms)
          REAL(DP) :: bfd(3), bfdr(3)

          INTEGER, DIMENSION(3) :: int8shape, int8ub, lbp16, int16shape, int16ub
          INTEGER(INT64) :: xub_int64

          INTEGER(INT16), DIMENSION(:,:,:), ALLOCATABLE :: bitcell_int16_array
          INTEGER(INT64), DIMENSION(:), ALLOCATABLE :: zcavcount
          INTEGER(INT64) :: ncavs, icav, locbase, xi_int64

          INTEGER(INT16) :: oneshift_int16, ncavs_int16

          INTEGER :: i_big_atom

          LOGICAL(2), DIMENSION(0:MASKR(15)) :: ncavlvec

          INTEGER :: big_atom_start
          INTEGER(INT16), DIMENSION(0:15), PARAMETER :: oneshift_int16_vec = &
                  ISHFT(1_INT16,(/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15/))
          INTEGER :: yi_chunkstart, yi_chunkstride, yi_chunkend, chunksize, chunksize_p16
          INTEGER(INT16), DIMENSION(ISHFT(1,12)) :: bitcell_int16_vec
          INTEGER(INT16) :: cavbits_int16

          !WRITE(*,*) "sector_setup started"

          need_atom_ti = read_atompair_rminsq .OR. (cbmc_cell_list_flag .AND. precalc_atompair_nrg)
          need_atomtypes = calc_rmin_flag .OR. (cbmc_cell_list_flag .AND. .NOT. precalc_atompair_nrg)
          need_charges = ANY(int_charge_style .NE. charge_none) .AND. cbmc_cell_list_flag .AND. .NOT. precalc_atompair_nrg

          !sp_live_atom_rsp = REAL(live_atom_rsp)
          !$OMP PARALLEL PRIVATE(ci,nca) &
          !$OMP PRIVATE(int8shape,int8ub,lbp16,int16shape,int16ub,xub_int64,ncavs,icav,locbase,xi_int64,oneshift_int16,i_big_atom) &
          !$OMP PRIVATE(ncavlvec,ncavs_int16,yi_chunkstart, yi_chunkstride, yi_chunkend, chunksize, chunksize_p16) &
          !$OMP PRIVATE(bitcell_int16_vec, cavbits_int16) &
          !$OMP PRIVATE(ibox,istart,is,inlive,inatoms,bsolvent,vlen,iend,box_vlen) &
          !$OMP PRIVATE(xyzi,xyzi_sp,xyzi_dp,xyzi_dp_spread,dxyzi_dp,drp,cell_rp,zi_mult,yi_mult) &
          !$OMP PRIVATE(l_ortho,l_inrange_vec,l_inrange_vec_old,xi,yi,zi,i_dim,l_switch,bitmask) &
          !$OMP PRIVATE(isolvent,priv_bitcell_int64,zlb,ylb,xub,xlb,yub,zub) &
          !$OMP PRIVATE(lbox,clr,bclr,cp_ub,cp_lb,lc,rp_ub,rp_lb,bcp_shift) &
          !$OMP PRIVATE(bcp,bcpx,bcpy,bcpz,bcps,ti) &
          !$OMP PRIVATE(boxlen,lbc,tgt_slice,src_slice,bit_tgt_slice,bit_src_slice,border_range) &
          !$OMP PRIVATE(bt,rcutsq,hl,length_cells_recip,hlcr) &
          !$OMP PRIVATE(h11,h21,h31,h12,h22,h32,h13,h23,h33) &
          !$OMP PRIVATE(dzi,dyi,dxi,zi2,yi2,xi2,cell_rsxp,cell_rsyp,cell_rszp) &
          !$OMP PRIVATE(dsp,drxp,dryp,drzp,rsq_vec,cbmc_cell_rsp_priv,ti_priv,atomtype_priv) &
          !$OMP PRIVATE(rxp,ryp,rzp,isp,n_interact,adj_iend,which_interact,bcd,sbe,xcp,ycp,zcp,i,j) &
          !$OMP PRIVATE(zcp_base_shift,ycp_base_shift,xcp_base_shift,zcp_base,ycp_base,xcp_base,xfer_int64) &
          !$OMP PRIVATE(sbe_ti,sbe_ti_mat,bfd,bfdr)
          IF (bitcell_flag) THEN
                  !$OMP WORKSHARE
                  bitcell_int64 = 0_INT64
                  !$OMP END WORKSHARE
          END IF
          !!$OMP SINGLE
          !WRITE(*,*) "bitcell_int64 initialized"
          !!$OMP END SINGLE
          DO ibox = 1, nbr_boxes
                istart = 1
                DO is_present = 1, nspecies_present
                        is = which_species_present(is_present)
                        inlive = nlive(is,ibox)
                        IF (inlive < 1) CYCLE
                        inatoms = natoms(is)
                        bsolvent = INT(species_list(is)%solvent_base,2)
                        vlen = inlive*inatoms
                        iend = istart + vlen - 1
                        !$OMP WORKSHARE
                        sp_live_atom_rsp(istart:iend,1:3,ibox) = RESHAPE(REAL(live_atom_rsp(1:inatoms,1:inlive,1:3,is_present,ibox),SP), &
                                (/ vlen, 3 /))
                        !$OMP END WORKSHARE
                        IF (l_not_all_exist) THEN
                                !$OMP WORKSHARE
                                i_exist(istart:iend) = RESHAPE(live_atom_exist(1:inatoms,1:inlive,is_present,ibox), (/ vlen /))
                                !$OMP END WORKSHARE
                        END IF
                        IF (need_atom_ti) THEN
                                !$OMP WORKSHARE
                                live_atom_ti(istart:iend,ibox) = RESHAPE(SPREAD(bsolvent+INT(vec123(1:inatoms),INT32),2,inlive), (/ vlen /))
                                !$OMP END WORKSHARE
                        END IF
                        IF (need_atomtypes) THEN
                                !$OMP WORKSHARE
                                live_atom_atomtypes(istart:iend,ibox) = RESHAPE(SPREAD( &
                                        nonbond_list(1:inatoms,is)%atom_type_number,2,inlive), (/ vlen /))
                                !$OMP END WORKSHARE
                        END IF
                        IF (need_charges) THEN
                                !$OMP WORKSHARE
                                sp_live_atom_rsp(istart:iend,4,ibox) = RESHAPE(SPREAD(REAL( &
                                        nonbond_list(1:inatoms,is)%charge,SP),2,inlive), (/ vlen /))
                                !$OMP END WORKSHARE
                        END IF
                        istart = istart + vlen
                END DO
                vlen = istart - 1
                box_vlen(ibox) = vlen
                IF (l_not_all_exist) THEN
                        IF (.NOT. ALL(i_exist(1:vlen))) THEN
                                !$OMP SINGLE
                                n_i_exist = 0
                                DO i = 1, vlen
                                        IF (i_exist(i)) THEN
                                                n_i_exist = n_i_exist + 1
                                                which_i_exist(n_i_exist) = i
                                        END IF
                                END DO
                                !$OMP END SINGLE
                                !$OMP WORKSHARE
                                sp_live_atom_rsp(1:n_i_exist,:,ibox) = sp_live_atom_rsp(which_i_exist(1:n_i_exist),:,ibox)
                                live_atom_ti(1:n_i_exist,ibox) = live_atom_ti(which_i_exist(1:n_i_exist),ibox)
                                live_atom_atomtypes(1:n_i_exist,ibox) = live_atom_atomtypes(which_i_exist(1:n_i_exist),ibox)
                                !$OMP END WORKSHARE
                                box_vlen(ibox) = n_i_exist
                        END IF
                END IF
          END DO

          DO ibox = 1, nbr_boxes
                l_ortho = box_list(ibox)%int_box_shape <= int_ortho
                !$OMP SINGLE
                DO i = 1, 3
                          IF (bitcell_flag) THEN
                                  IF (i == 1) THEN
                                          box_list(ibox)%length_bitcells(i) = MIN(8*INT(box_list(ibox)%face_distance(i) / &
                                                  (8.0_DP*box_list(ibox)%ideal_bitcell_length)),SHIFTL(1,12))
                                  ELSE
                                          box_list(ibox)%length_bitcells(i) = MIN(INT(box_list(ibox)%face_distance(i) / &
                                                  box_list(ibox)%ideal_bitcell_length),SHIFTL(1,12))

                                  END IF
                                  bitcell_H(:,i) = box_list(ibox)%length(:,i) / box_list(ibox)%length_bitcells(i)
                                  !bitcell_H_diag(i) = bitcell_H(i,i)
                                  box_list(ibox)%bit_cell_length_recip(i) = REAL(REAL( &
                                          box_list(ibox)%length_bitcells(i),DP)/box_list(ibox)%length(i,i),SP)
                          END IF
                          box_list(ibox)%length_cells(i) = INT(box_list(ibox)%face_distance(i)/max_rmin)
                          IF (MOD(box_list(ibox)%length_cells(i),2) .EQ. 0) box_list(ibox)%length_cells(i) = box_list(ibox)%length_cells(i) - 1
                          box_list(ibox)%cell_H_dp(:,i) = box_list(ibox)%length(:,i) / box_list(ibox)%length_cells(i)
                          !cell_H(:,i) = box_list(ibox)%length(:,i) / box_list(ibox)%length_cells(i)
                          box_list(ibox)%cell_H_diag(i) = box_list(ibox)%cell_H_dp(i,i)
                          !box_list(ibox)%cell_H_diag(i) = cell_H(i,i)
                          box_list(ibox)%cell_length_inv(i,:) = REAL(box_list(ibox)%length_cells(i),DP) * box_list(ibox)%length_inv(i,:)
                          box_list(ibox)%cell_length_recip(i) = REAL(REAL(box_list(ibox)%length_cells(i),DP)/box_list(ibox)%length(i,i),SP) ! kind 4
                          box_list(ibox)%sectorbound(i) = box_list(ibox)%length_cells(i)/2
                          box_list(ibox)%sp_diag_length(i) = REAL(box_list(ibox)%length(i,i),SP)
                END DO
                IF (bitcell_flag) THEN
                        bitcell_xyzortho_bbox_length = SUM(ABS(bitcell_H),2)
                        box_list(ibox)%real_length_bitcells = REAL(box_list(ibox)%length_bitcells,SP)
                        box_list(ibox)%bitcell_face_distance = box_list(ibox)%face_distance / box_list(ibox)%length_bitcells
                        box_list(ibox)%bitcell_face_distance_recip = box_list(ibox)%length_bitcells / box_list(ibox)%face_distance
                        bfd = box_list(ibox)%bitcell_face_distance
                        bfdr = box_list(ibox)%bitcell_face_distance_recip
                        DO i = 1, 3
                                box_list(ibox)%setbit_extent(i) = &
                                        MIN(INT(SQRT((box_list(ibox)%rcut_low_max*bfdr(i))**2 + 1.0_DP - SUM(bfd*bfdr(i))))-1,28)
                        END DO
                        !box_list(ibox)%setbit_extent = &
                        !        MIN(INT(SQRT((box_list(ibox)%rcut_low_max*box_list(ibox)%length_bitcells/box_list(ibox)%face_distance)**2-2)-1),28)
                        !box_list(ibox)%setbit_extent = &
                        !        MIN(INT(CEILING(box_list(ibox)%rcut_low_max*box_list(ibox)%length_bitcells/box_list(ibox)%face_distance)),28)
                        ! MIN(...,28) above is unlikely to be necessary but is included just to be safe.
                END IF
                box_list(ibox)%cell_H_sp = REAL(box_list(ibox)%cell_H_dp,SP)
                cell_H_dp = box_list(ibox)%cell_H_dp
                !WRITE(*,*) box_list(ibox)%length, box_list(ibox)%cell_H_dp, box_list(ibox)%cell_H_sp
                box_list(ibox)%cell_xyzortho_bbox_length = SUM(ABS(cell_H_dp),2)
                box_list(ibox)%real_length_cells = REAL(box_list(ibox)%length_cells,SP)
                box_list(ibox)%cell_face_distance = box_list(ibox)%face_distance / box_list(ibox)%length_cells
                IF (cbmc_cell_list_flag) THEN
                        box_list(ibox)%border_thickness = &
                                INT(CEILING(rcut_cbmc(ibox)*box_list(ibox)%length_cells/box_list(ibox)%face_distance))
                        IF (ALLOCATED(box_list(ibox)%cbmc_cell_mask)) THEN
                                IF (ANY(box_list(ibox)%border_thickness .NE. UBOUND(box_list(ibox)%cbmc_cell_mask))) THEN
                                        DEALLOCATE(box_list(ibox)%cbmc_cell_mask)
                                END IF
                        END IF
                        IF (.NOT. ALLOCATED(box_list(ibox)%cbmc_cell_mask)) THEN
                                ALLOCATE(box_list(ibox)%cbmc_cell_mask( &
                                        -box_list(ibox)%border_thickness(1):box_list(ibox)%border_thickness(1), &
                                        -box_list(ibox)%border_thickness(2):box_list(ibox)%border_thickness(2), &
                                        -box_list(ibox)%border_thickness(3):box_list(ibox)%border_thickness(3)))
                        END IF
                ELSE
                        box_list(ibox)%border_thickness = 1
                END IF
                !$OMP END SINGLE
                IF (cbmc_cell_list_flag) THEN
                        zub = box_list(ibox)%border_thickness(3)
                        yub = box_list(ibox)%border_thickness(2)
                        xub = box_list(ibox)%border_thickness(1)
                        zlb = -zub
                        ylb = -yub
                        xlb = -xub
                        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
                        DO zi = zlb, zub
                        DO yi = ylb, yub
                        DO xi = xlb, xub
                                xyzi = (/ xi, yi, zi /)
                                xyzi_dp = REAL(xyzi,DP)
                                IF (l_ortho) THEN
                                        !dxyzi_dp = SIGN(MAX(ABS(xyzi_dp)-1.0_DP,0.0_DP),xyzi_DP)
                                        !drp = dxyzi_dp*box_list(ibox)%cell_xyzortho_bbox_length
                                        drp = MAX(ABS(xyzi_dp)-1.0_DP,0.0_DP)*box_list(ibox)%cell_xyzortho_bbox_length
                                ELSE
                                        !cellcenter_drp = MATMUL(cell_H_dp,xyzi_dp)
                                        !drp = SIGN(MAX(ABS(cellcenter_drp)-box_list(ibox)%cell_xyzortho_bbox_length,0.0_DP),cellcenter_drp)
                                        drp = MAX(ABS(MATMUL(cell_H_dp,xyzi_dp))-box_list(ibox)%cell_xyzortho_bbox_length,0.0_DP)
                                        !xyzi_dp_spread = SPREAD(SPREAD(xyzi_dp,2,3),3,2)
                                        !drp = MINVAL( &
                                        !        ABS(SUM((xyzi_dp_spread+MERGE(SPREAD(SPREAD(SPREAD(0.0_DP,1,3),2,3),3,2), &
                                        !        cell_H_T_posneg_sign,xyzi_dp_spread .EQ. 0.0_DP))*cell_H_T_spread,1)),2)
                                END IF
                                box_list(ibox)%cbmc_cell_mask(xi,yi,zi) = &
                                        MERGE(NOT(0),0,DOT_PRODUCT(drp,drp)<rcut_cbmcsq(ibox))
                        END DO
                        END DO
                        END DO
                        !$OMP END DO
                END IF
                IF (bitcell_flag) THEN
                        vlen = solvents_or_types_maxind
                        sbe = box_list(ibox)%setbit_extent
                        zlb = MERGE(0,-sbe(3),l_ortho)
                        ylb = MERGE(0,-sbe(2),l_ortho)
                        xub = MERGE(0,sbe(1)+1,l_ortho)
                        xlb = -sbe(1)
                        yub = sbe(2)
                        zub = sbe(3)
                        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
                        DO zi = zlb, zub
                        DO yi = ylb, yub
                        priv_bitcell_int64 = 0_INT64
                        l_inrange_vec = .FALSE.
                        DO xi = xlb, xub
                                xyzi = (/ xi, yi, zi /)
                                xyzi_dp = REAL(xyzi,DP)
                                IF (l_ortho) THEN
                                        !dxyzi_dp = SIGN(ABS(xyzi_dp)+1.0_DP,xyzi_dp)
                                        !dxyzi_dp = ABS(xyzi_dp)+1.0_DP
                                        drp = ABS(xyzi_dp)*bitcell_xyzortho_bbox_length + bitcell_xyzortho_bbox_length
                                ELSE
                                        ! This is a simplified, less accurate way of finding the maximum possible atomic separation
                                        ! squared for the two triclinic bitcells for triclinic boxes.  The inaccuracy is acceptable because it can
                                        ! never underestimate the maximum separation (overestimation is safe, just less optimized),
                                        ! and the inaccuracy should be small, especially if the box is only slightly tilted.
                                        ! It basically finds the maximum possible absolute value of the x, y, and z displacements
                                        ! separately, and uses the resulting displacement vector to compute the separation distance squared.
                                        ! Inaccuracy can arise from the maximum absolute displacement vector lying outside the
                                        ! bitcell/fractional coordinate displacement inequality constraints of the multivariate
                                        ! optimization problem we're trying to solve.
                                        !cellcenter_drp = MATMUL(bitcell_H,xyzi_dp)
                                        !drp = SIGN(ABS(cellcenter_drp)+bitcell_xyzortho_bbox_length,cellcenter_drp)
                                        ! below is absolute value of actual displacement because sign doesn't matter
                                        drp = ABS(MATMUL(bitcell_H,xyzi_dp))+bitcell_xyzortho_bbox_length
                                        ! below is old way that should work but is probably slower than what we use despite giving
                                        ! the same results
                                        !drp = MAXVAL( &
                                        !        ABS(SUM((SPREAD(SPREAD(xyzi_dp,2,3),3,2)+bitcell_H_T_posneg_sign)*bitcell_H_T_spread,1)),2)
                                        ! below is old way that doesn't work for all triclinic cells
                                        ! drp = MATMUL(bitcell_H,dxyzi_dp)
                                END IF
                                rsq = DOT_PRODUCT(drp,drp)
                                l_inrange_vec_old = l_inrange_vec
                                IF (read_atompair_rminsq) THEN
                                        l_inrange_vec(1:vlen,:) = rsq<solvent_min_rminsq(:,:,ibox)
                                ELSE IF (calc_rmin_flag) THEN
                                        l_inrange_vec(1:vlen,:) = rsq<atomtype_min_rminsq
                                ELSE
                                        l_inrange_vec(1,0) = rsq < rcut_lowsq
                                END IF
                                l_switch = l_inrange_vec .AND. .NOT. l_inrange_vec_old
                                IF (ANY(l_switch)) THEN
                                        IF (l_ortho) THEN
                                                bitmask = IAND(MASKL(36-xi,INT64),MASKR(29-xi,INT64))
                                        ELSE
                                                bitmask = MASKL(36-xi,8)
                                        END IF
                                        DO i_big_atom = 0, n_big_atoms
                                                DO isolvent = 1, vlen
                                                        IF (l_switch(isolvent,i_big_atom)) THEN
                                                                priv_bitcell_int64(isolvent,i_big_atom) = bitmask
                                                        END IF
                                                END DO
                                        END DO
                                END IF
                                IF (l_ortho) CYCLE
                                l_switch = l_inrange_vec_old .AND. .NOT. l_inrange_vec
                                IF (ANY(l_switch)) THEN
                                        bitmask = MASKR(28+xi,8)
                                        DO i_big_atom = 0, n_big_atoms
                                                DO isolvent = 1, vlen
                                                        IF (l_switch(isolvent,i_big_atom)) THEN
                                                                priv_bitcell_int64(isolvent,i_big_atom) = &
                                                                        IAND(priv_bitcell_int64(isolvent,i_big_atom),bitmask)
                                                        END IF
                                                END DO
                                        END DO
                                END IF
                        END DO
                        !DO zi_mult = MERGE(-1,1,l_ortho .AND. zi>0),1,2
                        !        DO yi_mult = MERGE(-1,1,l_ortho .AND. yi>0),1,2
                        !                bitcell_int64(yi*yi_mult,zi*zi_mult,1:vlen,ibox) = priv_bitcell_int64(1:vlen)
                        !        END DO
                        !END DO
                        !DO zi_mult = MERGE(-zi,zi,l_ortho .AND. zi>0),zi,MIN(2*zi,1)
                        !        DO yi_mult = MERGE(-yi,yi,l_ortho .AND. yi>0),yi,MIN(2*yi,1)
                        !                bitcell_int64(yi_mult,zi_mult,1:vlen,ibox) = priv_bitcell_int64(1:vlen)
                        !        END DO
                        !END DO
                        bitcell_int64(yi,zi,1:vlen,:,ibox) = priv_bitcell_int64(1:vlen,:)
                        END DO
                        END DO
                        !$OMP END DO
                        IF (l_ortho) THEN
                                !$OMP WORKSHARE
                                bitcell_int64(-1:-sbe(2):-1,-1:-sbe(3):-1,:,:,ibox) = &
                                        bitcell_int64(1:sbe(2),1:sbe(3),:,:,ibox)
                                bitcell_int64(0:sbe(2),-1:-sbe(3):-1,:,:,ibox) = &
                                        bitcell_int64(0:sbe(2),1:sbe(3),:,:,ibox)
                                bitcell_int64(-1:-sbe(2):-1,0:sbe(3),:,:,ibox) = &
                                        bitcell_int64(1:sbe(2),0:sbe(3),:,:,ibox)
                                !$OMP END WORKSHARE
                        END IF
                END IF
          END DO
          !$OMP BARRIER
          !!$OMP SINGLE
          !WRITE(*,*) "bitcell_int64 filled", SHAPE(bitcell_int64)
          !!$OMP END SINGLE
          !!$OMP SINGLE
          !WRITE(*,*) "bitcell_int64 slices:"
          !WRITE(*,*) bitcell_int64(-10,-10,:,1,0)
          !WRITE(*,*) bitcell_int64(0,0,:,1,0)
          !WRITE(*,*) bitcell_int64(15,5,:,1,0)
          !WRITE(*,*) "bitcell_int64 slices LEADZ - TRAILZ:"
          !WRITE(*,*) LEADZ(bitcell_int64(-10,-10,:,1,0)) - TRAILZ(bitcell_int64(-10,-10,:,1,0))
          !WRITE(*,*) LEADZ(bitcell_int64(0,0,:,1,0)) - TRAILZ(bitcell_int64(0,0,:,1,0))
          !WRITE(*,*) LEADZ(bitcell_int64(15,5,:,1,0)) - TRAILZ(bitcell_int64(15,5,:,1,0))
          !WRITE(*,*) "bitcell_int64 type 0 8bit transfer:"
          !WRITE(*,*) TRANSFER(bitcell_int64(-10,-10,1,1,0),bitcell_int8_array)
          !WRITE(*,*) TRANSFER(bitcell_int64(0,0,1,1,0),bitcell_int8_array)
          !WRITE(*,*) TRANSFER(bitcell_int64(15,5,1,1,0),bitcell_int8_array)
          !!$OMP END SINGLE
          !IF (bitcell_flag) THEN
          !        DO i = 1, 7
          !                !$OMP WORKSHARE
          !                bitcell_int64(:,:,:,:,i) = ISHFT(bitcell_int64(:,:,:,:,0),i)
          !                !$OMP END WORKSHARE NOWAIT
          !        END DO
          !        !$OMP BARRIER
          !END IF
          !!$OMP SINGLE
          !WRITE(*,*) "bitcell_int64 slices shift 4:"
          !WRITE(*,*) bitcell_int64(-10,-10,:,1,4)
          !WRITE(*,*) bitcell_int64(0,0,:,1,4)
          !WRITE(*,*) bitcell_int64(15,5,:,1,4)
          !!$OMP END SINGLE
          DO ibox = 1, nbr_boxes
                vlen = box_vlen(ibox)
                l_ortho = box_list(ibox)%int_box_shape <= int_ortho
                DO i_dim = 1, 3
                        IF (l_ortho) THEN
                                lbox = box_list(ibox)%sp_diag_length(i_dim)
                                clr = box_list(ibox)%cell_length_recip(i_dim)
                                bclr = box_list(ibox)%bit_cell_length_recip(i_dim)
                        ELSE
                                lbox = 1.0
                                clr = box_list(ibox)%real_length_cells(i_dim)
                                bclr = box_list(ibox)%real_length_bitcells(i_dim)
                        END IF
                        cp_ub = box_list(ibox)%sectorbound(i_dim)
                        cp_lb = -cp_ub
                        lc = box_list(ibox)%length_cells(i_dim)
                        rp_ub = 0.5*lbox
                        rp_lb = -rp_ub
                        bcp_shift = MERGE(4,32,i_dim==1)
                        !DIR$ VECTOR ALIGNED
                        !$OMP DO SIMD PRIVATE(rsp,cp) SCHEDULE(SIMD:STATIC)
                        DO i = 1, vlen
                                rsp = sp_live_atom_rsp(i,i_dim,ibox)
                                IF (rsp > rp_ub) THEN
                                        rsp = rsp - lbox
                                ELSE IF (rsp < rp_lb) THEN
                                        rsp = rsp + lbox
                                END IF
                                sp_live_atom_rsp(i,i_dim,ibox) = rsp
                                cp = NINT(rsp*clr)
                                ! Wrapping cell coordinates is rarely necessary, but it is done anyway just in case
                                ! rsp is exactly on a box boundary, which would yield out-of-bounds cell coordinates
                                ! if not corrected.  The risk is increased because single-precision floats are used.
                                IF (cp > cp_ub) THEN
                                        cp = cp - lc
                                ELSE IF (cp < cp_lb) THEN
                                        cp = cp + lc
                                END IF
                                live_atom_cp(i,i_dim,ibox) = cp
                                ! bit cell int64 array placement start
                                ! make sure this conditional is brought out of the loop by compiler.
                                ! bitcell lower bound (min bitcell index in box) is 32 (int 4, bit 0 for first axis)
                                !   in order to avoid performing integer division or modulo on negative int,
                                !   which would throw off the indexing.  The methods used require
                                !   integer division to function like float division followed by floor.
                                !IF (bitcell_flag) THEN
                                !        live_atom_bcp(i,i_dim,ibox) = INT((rsp+rp_ub)*bclr) + 4
                                !END IF
                                IF (bitcell_flag) THEN
                                        live_atom_bcp(i,i_dim,ibox) = INT((rsp+rp_ub)*bclr) + bcp_shift
                                END IF
                        END DO
                        !$OMP END DO SIMD
                END DO
          END DO
          !$OMP BARRIER
          !!$OMP SINGLE
          !WRITE(*,*) "live_atom_bcp filled"!, SHAPE(live_atom_bcp)
          !dummy3vec = SHAPE(live_atom_bcp)
          !!$OMP END SINGLE
          !$OMP SINGLE
          sectormaxbound = box_list(1)%sectorbound + box_list(1)%border_thickness
          DO ibox = 2, nbr_boxes
                sectormaxbound = MAX(sectormaxbound,box_list(ibox)%sectorbound + box_list(ibox)%border_thickness)
          END DO

          !sectormaxbound = MAXVAL(sectorbound, 2)+1
          adj_cellmaxbound_old = adj_cellmaxbound
          !adj_cellmaxbound = box_list(1)%sectorbound
          DO ibox = 1, nbr_boxes
                adj_cellmaxbound = MAX(adj_cellmaxbound,box_list(ibox)%sectorbound)
          END DO
          !map_bound = sectormaxbound*3+1
          ALLOCATE(n_cell_atoms( &
                  -sectormaxbound(1):sectormaxbound(1), &
                  -sectormaxbound(2):sectormaxbound(2), &
                  -sectormaxbound(3):sectormaxbound(3), &
                  nbr_boxes), Stat=AllocateStatus)
          IF (Allocatestatus /= 0) THEN
            err_msg = ''
            err_msg(1) = 'Memory could not be allocated for n_cell_atoms'
            CALL Clean_Abort(err_msg, 'Sector_Setup')
          END IF
          IF (ALLOCATED(n_adj_cell_atoms)) THEN
                  IF (ANY(adj_cellmaxbound > adj_cellmaxbound_old)) THEN
                          DEALLOCATE(n_adj_cell_atoms, Stat=AllocateStatus)
                          IF (Allocatestatus /= 0) THEN
                            err_msg = ''
                            err_msg(1) = 'Memory could not be deallocated from n_adj_cell_atoms'
                            CALL Clean_Abort(err_msg, 'Sector_Setup')
                          END IF
                  END IF
          END IF
          IF (.NOT. ALLOCATED(n_adj_cell_atoms)) THEN
                  ALLOCATE(n_adj_cell_atoms( &
                          -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                          -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                          -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                          nbr_boxes), Stat=AllocateStatus)
                  IF (Allocatestatus /= 0) THEN
                    err_msg = ''
                    err_msg(1) = 'Memory could not be allocated for n_adj_cell_atoms'
                    CALL Clean_Abort(err_msg, 'Sector_Setup')
                  END IF
                  WRITE(*,*) "Allocated n_adj_cell_atoms with:"
                  WRITE(*,*) "SHAPE", SHAPE(n_adj_cell_atoms)
                  WRITE(*,*) "Occupying ", SIZEOF(n_adj_cell_atoms), " bytes"
                  WRITE(*,*) "on step ", i_mcstep
          END IF
          !WRITE(*,*) "n_adj_cell_atoms allocated", SHAPE(n_adj_cell_atoms)
          max_adj_cell_atoms_old = max_adj_cell_atoms
          max_adj_cell_atoms = 0
          max_neighbors = 0
          !$OMP END SINGLE
          !$OMP WORKSHARE
          n_cell_atoms = 0
          !$OMP END WORKSHARE
          !$OMP SINGLE
          DO ibox = 1, nbr_boxes
                DO i = 1, box_vlen(ibox)
                        ci(1:3) = live_atom_cp(i,:,ibox)
                        nca = n_cell_atoms(ci(1),ci(2),ci(3),ibox)+1
                        ci(4) = nca
                        n_cell_atoms(ci(1),ci(2),ci(3),ibox) = nca
                        ci_list(:,i,ibox) = ci
                END DO
          END DO
          !$OMP END SINGLE
          !$OMP WORKSHARE
          max_sector_natoms = MAXVAL(n_cell_atoms)
          !$OMP END WORKSHARE
          !$OMP SINGLE
          ALLOCATE(this_cell_rsp(max_sector_natoms, 4, &
                  -sectormaxbound(1):sectormaxbound(1), &
                  -sectormaxbound(2):sectormaxbound(2), &
                  -sectormaxbound(3):sectormaxbound(3), &
                  nbr_boxes), Stat=AllocateStatus)
          IF (Allocatestatus /= 0) THEN
            err_msg = ''
            err_msg(1) = 'Memory could not be allocated for this_cell_rsp'
            CALL Clean_Abort(err_msg, 'Sector_Setup')
          END IF
          IF (need_atom_ti) THEN
                  ALLOCATE(this_cell_ti(max_sector_natoms, &
                          -sectormaxbound(1):sectormaxbound(1), &
                          -sectormaxbound(2):sectormaxbound(2), &
                          -sectormaxbound(3):sectormaxbound(3), &
                          nbr_boxes), Stat=AllocateStatus)
                  IF (Allocatestatus /= 0) THEN
                    err_msg = ''
                    err_msg(1) = 'Memory could not be allocated for this_cell_ti'
                    CALL Clean_Abort(err_msg, 'Sector_Setup')
                  END IF
          END IF
          IF (need_atomtypes) THEN
                  ALLOCATE(this_cell_atomtypes(max_sector_natoms, &
                          -sectormaxbound(1):sectormaxbound(1), &
                          -sectormaxbound(2):sectormaxbound(2), &
                          -sectormaxbound(3):sectormaxbound(3), &
                          nbr_boxes), Stat=AllocateStatus)
                  IF (Allocatestatus /= 0) THEN
                    err_msg = ''
                    err_msg(1) = 'Memory could not be allocated for this_cell_atomtypes'
                    CALL Clean_Abort(err_msg, 'Sector_Setup')
                  END IF
          END IF
          IF (bitcell_flag) THEN
                  ALLOCATE(this_cell_atom_i(max_sector_natoms, &
                          -sectormaxbound(1):sectormaxbound(1), &
                          -sectormaxbound(2):sectormaxbound(2), &
                          -sectormaxbound(3):sectormaxbound(3)), &
                          Stat=AllocateStatus)
                  IF (Allocatestatus /= 0) THEN
                    err_msg = ''
                    err_msg(1) = 'Memory could not be allocated for this_cell_atom_i'
                    CALL Clean_Abort(err_msg, 'Sector_Setup')
                  END IF
          END IF
          !$OMP END SINGLE
          DO ibox = 1, nbr_boxes
                vlen = box_vlen(ibox)
                !$OMP DO SCHEDULE(STATIC)
                DO i = 1, vlen
                        ci = ci_list(:,i,ibox)
                        this_cell_rsp(ci(4),1:4,ci(1),ci(2),ci(3),ibox) = &
                                sp_live_atom_rsp(i,:,ibox)
                        IF (need_atom_ti) THEN
                                this_cell_ti(ci(4),ci(1),ci(2),ci(3),ibox) = &
                                        live_atom_ti(i,ibox)
                        END IF
                        IF (need_atomtypes) THEN
                                this_cell_atomtypes(ci(4),ci(1),ci(2),ci(3),ibox) = &
                                        live_atom_atomtypes(i,ibox)
                        END IF
                        IF (bitcell_flag) THEN
                                this_cell_atom_i(ci(4),ci(1),ci(2),ci(3)) = i
                        END IF
                END DO
                !$OMP END DO
                l_ortho = box_list(ibox)%int_box_shape <= int_ortho
                IF (l_ortho) THEN
                        DO i_dim = 1, 3
                                boxlen(i_dim) = REAL(box_list(ibox)%length(i_dim,i_dim),SP)
                        END DO
                ELSE
                        boxlen = 1.0
                END IF
                IF (bitcell_flag) THEN
                        sbe = box_list(ibox)%setbit_extent
                        !$OMP SINGLE
                        lbp32 = box_list(ibox)%length_bitcells
                        ALLOCATE(bitcell_int8_array( &
                                0:7+lbp32(1)/8, &
                                0:63+lbp32(2), &
                                0:63+lbp32(3),&
                                0:n_big_atoms), Stat=AllocateStatus)
                        IF (Allocatestatus /= 0) THEN
                          err_msg = ''
                          err_msg(1) = 'Memory could not be allocated for bitcell_int8_array'
                          CALL Clean_Abort(err_msg, 'Sector_Setup')
                        END IF
                        lbp32(1) = IAND(lbp32(1)+31,NOT(31))
                        !IF (MOD(lbp32(1),32) .NE. 0) THEN
                        !        lbp32(1) = (lbp32(1)/32+1)*32
                        !END IF
                        !$OMP END SINGLE
                        !$OMP WORKSHARE
                        bitcell_int8_array = 0_INT8
                        !$OMP END WORKSHARE
                        !DIR$ VECTOR ALIGNED
                        !$OMP DO SIMD PRIVATE(bcp) SCHEDULE(SIMD:STATIC)
                        DO i = 1, vlen
                                bcp = live_atom_bcp(i,1,ibox)
                                live_atom_bcp(i,1,ibox) = ISHFT(bcp,-3) ! bcp / 8
                                live_atom_bcp(i,4,ibox) = IAND(bcp,7) ! MOD(bcp,8)
                        END DO
                        !$OMP END DO SIMD
                        !!$OMP SINGLE
                        !WRITE(*,*) "bitcell_int8_array bounds: "
                        !WRITE(*,*) LBOUND(bitcell_int8_array)
                        !WRITE(*,*) UBOUND(bitcell_int8_array)
                        !WRITE(*,*) "live_atom_bcp(1:50,:,1) = "
                        !DO i = 1, 50
                        !        WRITE(*,*) live_atom_bcp(i,:,1)
                        !END DO
                        !!$OMP END SINGLE



                        !DIR$ ASSUME_ALIGNED live_atom_bcp_T:32
                        !$OMP WORKSHARE
                        live_atom_bcp_T = TRANSPOSE(live_atom_bcp(:,:,ibox))
                        xcp_base_shift_stride = INT(CEILING( &
                                72.0_DP*REAL(box_list(ibox)%length_cells(1),DP)/REAL(box_list(ibox)%length_bitcells(1),DP)))
                        xcp_base_stride = 2*xcp_base_shift_stride
                        !$OMP END WORKSHARE
                        bfd = box_list(ibox)%bitcell_face_distance
                        bfdr = box_list(ibox)%bitcell_face_distance_recip
                        IF (read_atompair_rminsq) THEN
                                sbe_ti_mat(2,:,:) = MIN(INT(SQRT(solvent_min_rminsq(:,:,ibox)*bfdr(2)*bfdr(2) - bfd(1)*bfdr(2) - &
                                        bfd(3)*bfdr(2)))-1,28)
                                sbe_ti_mat(3,:,:) = MIN(INT(SQRT(solvent_min_rminsq(:,:,ibox)*bfdr(3)*bfdr(3) - bfd(1)*bfdr(3) - &
                                        bfd(2)*bfdr(3)))-1,28)
                        ELSE IF (calc_rmin_flag) THEN
                                sbe_ti_mat(2,:,:) = MIN(INT(SQRT(atomtype_min_rminsq*bfdr(2)*bfdr(2) - bfd(1)*bfdr(2) - &
                                        bfd(3)*bfdr(2)))-1,28)
                                sbe_ti_mat(3,:,:) = MIN(INT(SQRT(atomtype_min_rminsq*bfdr(3)*bfdr(3) - bfd(1)*bfdr(3) - &
                                        bfd(2)*bfdr(3)))-1,28)
                        ELSE
                                sbe_ti = sbe(2:3)
                                ti = 1
                        END IF
                        !!$OMP SINGLE
                        !WRITE(*,*) "xcp_base_shift_stride", xcp_base_shift_stride
                        !!$OMP END SINGLE
                        DO zcp_base_shift = 0,2,2
                        DO ycp_base_shift = 0,2,2
                        DO xcp_base_shift = 0,xcp_base_shift_stride,xcp_base_shift_stride
                        zub = box_list(ibox)%sectorbound(3)
                        zlb = zcp_base_shift-zub
                        yub = box_list(ibox)%sectorbound(2)
                        ylb = ycp_base_shift-yub
                        xub = box_list(ibox)%sectorbound(1)
                        xlb = xcp_base_shift-xub
                        !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
                        DO zcp_base = zlb, zub, 4
                        DO ycp_base = ylb, yub, 4
                        DO xcp_base = xlb, xub, xcp_base_stride
                        DO zcp = zcp_base, MIN(zcp_base+1,box_list(ibox)%sectorbound(3))
                        DO ycp = ycp_base, MIN(ycp_base+1,box_list(ibox)%sectorbound(2))
                        DO xcp = xcp_base, MIN(xcp_base+xcp_base_shift_stride-1,box_list(ibox)%sectorbound(1))
                        DO j = 1, n_cell_atoms(xcp,ycp,zcp,ibox)
                                i = this_cell_atom_i(j,xcp,ycp,zcp)
                                bcpx = live_atom_bcp_T(1,i)
                                bcpy = live_atom_bcp_T(2,i)
                                bcpz = live_atom_bcp_T(3,i)
                                bcps = live_atom_bcp_T(4,i)
                                IF (read_atompair_rminsq) THEN
                                        ti = live_atom_ti(i,ibox)
                                ELSE IF (calc_rmin_flag) THEN
                                        ti = which_solvent_atomtypes_inv(live_atom_atomtypes(i,ibox))
                                END IF
                                DO i_big_atom = 0, n_big_atoms
                                        IF (read_atompair_rminsq .OR. calc_rmin_flag) sbe_ti = sbe_ti_mat(:,ti,i_big_atom)
                                        DO zi = -sbe_ti(3), sbe_ti(3)
                                                DO yi = -sbe_ti(2), sbe_ti(2)
                                                        bitmask = bitcell_int64(yi,zi,ti,i_big_atom,ibox)
                                                        IF (bitmask .EQ. 0_INT64) CYCLE ! IOR with 0 changes nothing
                                                        bitmask = ISHFT(bitmask,bcps)
                                                        xfer_int64 = TRANSFER(bitcell_int8_array(bcpx:bcpx+7,bcpy+yi,bcpz+zi,i_big_atom),xfer_int64)
                                                        bitmask = IOR(xfer_int64,bitmask)
                                                        IF (bitmask == xfer_int64) CYCLE ! No point in writing what is already there
                                                        bitcell_int8_array(bcpx:bcpx+7,bcpy+yi,bcpz+zi,i_big_atom) = &
                                                                TRANSFER(bitmask,bitcell_int8_array)
                                                        !bitcell_int8_array(bcpx:bcpx+7,bcpy+yi,bcpz+zi) = &
                                                        !        TRANSFER(IOR(TRANSFER(&
                                                        !        bitcell_int8_array(bcpx:bcpx+7,bcpy+yi,bcpz+zi),0_INT64), &
                                                        !        ISHFT(bitmask,bcps)),bitcell_int8_array)
                                                        !bitcell_int8_array(bcpx:bcpx+7,bcpy+yi,bcpz+zi) = &
                                                        !        IOR(bitcell_int8_array(bcpx:bcpx+7,bcpy+yi,bcpz+zi), &
                                                        !        TRANSFER(bitcell_int64(yi,zi,ti,ibox,bcps),bitcell_int8_array))
                                                END DO
                                        END DO
                                END DO
                        END DO
                        END DO
                        END DO
                        END DO
                        END DO
                        END DO
                        END DO
                        !$OMP END DO
                        END DO
                        END DO
                        END DO
                        !DO i = 1, vlen
                        !        !!$OMP SINGLE
                        !        !IF (ANY(live_atom_bcp(i,1:3,ibox) < LBOUND(bitcell_int8_array) .OR. &
                        !        !        live_atom_bcp(i,1:3,ibox) + (/ 7, 56, 56 /) > &
                        !        !        UBOUND(bitcell_int8_array))) THEN
                        !        !        WRITE(*,*) "ATOM ", i, " is out of bounds!"
                        !        !        WRITE(*,*) "bad coordinates are :"
                        !        !        WRITE(*,*) live_atom_bcp(i,:,ibox)
                        !        !END IF
                        !        !!$OMP END SINGLE
                        !        bcpx = live_atom_bcp(i,1,ibox)
                        !        bcpy = live_atom_bcp(i,2,ibox)
                        !        bcpz = live_atom_bcp(i,3,ibox)
                        !        bcps = live_atom_bcp(i,4,ibox)
                        !        IF (read_atompair_rminsq) THEN
                        !                ti = live_atom_ti(i,ibox)
                        !        ELSE IF (calc_rmin_flag) THEN
                        !                ti = live_atom_atomtypes(i,ibox)+1
                        !        ELSE
                        !                ti = 1
                        !        END IF
                        !        sbe = box_list(ibox)%setbit_extent
                        !        !$OMP DO COLLAPSE(2) SCHEDULE(STATIC)
                        !        DO zi = -sbe(3), sbe(3)
                        !                DO yi = -sbe(2), sbe(2)
                        !                        bitmask = bitcell_int64(yi,zi,ti,ibox,bcps)
                        !                        IF (bitmask .NE. 0_INT64) THEN
                        !                                bitcell_int8_array(bcpx:bcpx+7,bcpy+yi,bcpz+zi) = &
                        !                                        TRANSFER(IOR(TRANSFER(&
                        !                                        bitcell_int8_array(bcpx:bcpx+7,bcpy+yi,bcpz+zi),0_INT64), &
                        !                                        bitmask),bitcell_int8_array)
                        !                        END IF
                        !                        !bitcell_int8_array(bcpx:bcpx+7,bcpy+yi,bcpz+zi) = &
                        !                        !        IOR(bitcell_int8_array(bcpx:bcpx+7,bcpy+yi,bcpz+zi), &
                        !                        !        TRANSFER(bitcell_int64(yi,zi,ti,ibox,bcps),bitcell_int8_array))
                        !                END DO
                        !        END DO
                        !        !$OMP END DO
                        !END DO
                END IF
                DO xi = -1, 1
                        DO yi = -1, 1
                                DO zi = -1, 1
                                        IF (xi == 0 .AND. yi == 0 .AND. zi == 0) CYCLE
                                        xyzi = (/ xi, yi, zi /)
                                        unwrap_shifter = xyzi*boxlen
                                        DO i_dim = 1, 3
                                                IF (bitcell_flag) THEN
                                                        lbc = box_list(ibox)%length_bitcells(i_dim)
                                                        IF (i_dim == 1) lbc = ISHFT(lbc,-3) ! same as lbc / 8 since lbc>0
                                                END IF
                                                IF (xyzi(i_dim) == 0) THEN
                                                        tgt_slice(:,i_dim) = (/ -1, 1 /) * box_list(ibox)%sectorbound(i_dim)
                                                        src_slice(:,i_dim) = tgt_slice(:,i_dim)
                                                        IF (bitcell_flag) THEN
                                                                IF (i_dim == 1) THEN
                                                                        bit_tgt_slice(:,1) = (/ 4, 3+lbc /)
                                                                ELSE
                                                                        bit_tgt_slice(:,i_dim) = (/ 32, &
                                                                                31 + lbc /)
                                                                END IF
                                                                bit_src_slice(:,i_dim) = bit_tgt_slice(:,i_dim)
                                                        END IF
                                                ELSE
                                                        IF (xyzi(i_dim) == 1) THEN
                                                                border_range = (/ 1, box_list(ibox)%border_thickness(i_dim) /)
                                                        ELSE
                                                                border_range = (/ box_list(ibox)%border_thickness(i_dim), 1 /)
                                                        END IF
                                                        tgt_slice(:,i_dim) = box_list(ibox)%sectorbound(i_dim)*xyzi(i_dim) + &
                                                                xyzi(i_dim)*border_range
                                                        src_slice(:,i_dim) = -box_list(ibox)%sectorbound(i_dim)*xyzi(i_dim) + &
                                                                xyzi(i_dim)*(border_range-1)
                                                        IF (bitcell_flag) THEN
                                                                IF (i_dim == 1) THEN
                                                                        IF (xyzi(i_dim) == 1) THEN
                                                                                bit_src_slice(:,1) = lbc + (/ 4, 7 /)
                                                                                bit_tgt_slice(:,1) = (/ 4, 7 /)
                                                                        ELSE
                                                                                bit_src_slice(:,1) = (/ 0, 3 /)
                                                                                bit_tgt_slice(:,1) = lbc + (/ 0, 3 /)
                                                                        END IF
                                                                ELSE
                                                                        IF (xyzi(i_dim) == 1) THEN
                                                                                bit_src_slice(:,i_dim) = lbc + (/ 32, 63 /)
                                                                                bit_tgt_slice(:,i_dim) = (/ 32, 63 /)
                                                                        ELSE
                                                                                bit_src_slice(:,i_dim) = (/ 0, 31 /)
                                                                                bit_tgt_slice(:,i_dim) = lbc + (/ 0, 31 /)
                                                                        END IF
                                                                END IF
                                                        END IF
                                                END IF
                                        END DO
                                        IF (bitcell_flag) THEN
                                                !$OMP WORKSHARE
                                                bitcell_int8_array( &
                                                        bit_tgt_slice(1,1):bit_tgt_slice(2,1), &
                                                        bit_tgt_slice(1,2):bit_tgt_slice(2,2), &
                                                        bit_tgt_slice(1,3):bit_tgt_slice(2,3),:) &
                                                        = IOR(bitcell_int8_array( &
                                                        bit_tgt_slice(1,1):bit_tgt_slice(2,1), &
                                                        bit_tgt_slice(1,2):bit_tgt_slice(2,2), &
                                                        bit_tgt_slice(1,3):bit_tgt_slice(2,3),:), &
                                                        bitcell_int8_array( &
                                                        bit_src_slice(1,1):bit_src_slice(2,1), &
                                                        bit_src_slice(1,2):bit_src_slice(2,2), &
                                                        bit_src_slice(1,3):bit_src_slice(2,3),:))
                                                !$OMP END WORKSHARE NOWAIT
                                        END IF
                                        DO i_dim = 1, 3
                                                !$OMP WORKSHARE
                                                this_cell_rsp(:,i_dim, &
                                                        tgt_slice(1,1):tgt_slice(2,1), &
                                                        tgt_slice(1,2):tgt_slice(2,2), &
                                                        tgt_slice(1,3):tgt_slice(2,3), &
                                                        ibox) = &
                                                        this_cell_rsp(:,i_dim, &
                                                        src_slice(1,1):src_slice(2,1), &
                                                        src_slice(1,2):src_slice(2,2), &
                                                        src_slice(1,3):src_slice(2,3), &
                                                        ibox) + unwrap_shifter(i_dim)
                                                !$OMP END WORKSHARE
                                        END DO
                                        !$OMP WORKSHARE
                                        this_cell_rsp(:,4, &
                                                tgt_slice(1,1):tgt_slice(2,1), &
                                                tgt_slice(1,2):tgt_slice(2,2), &
                                                tgt_slice(1,3):tgt_slice(2,3), &
                                                ibox) = &
                                                this_cell_rsp(:,4, &
                                                src_slice(1,1):src_slice(2,1), &
                                                src_slice(1,2):src_slice(2,2), &
                                                src_slice(1,3):src_slice(2,3), &
                                                ibox) ! charge
                                        n_cell_atoms(&
                                                tgt_slice(1,1):tgt_slice(2,1), &
                                                tgt_slice(1,2):tgt_slice(2,2), &
                                                tgt_slice(1,3):tgt_slice(2,3), &
                                                ibox) = &
                                                n_cell_atoms(&
                                                src_slice(1,1):src_slice(2,1), &
                                                src_slice(1,2):src_slice(2,2), &
                                                src_slice(1,3):src_slice(2,3), &
                                                ibox)
                                        !$OMP END WORKSHARE
                                        IF (need_atom_ti) THEN
                                                !$OMP WORKSHARE
                                                this_cell_ti(:, &
                                                        tgt_slice(1,1):tgt_slice(2,1), &
                                                        tgt_slice(1,2):tgt_slice(2,2), &
                                                        tgt_slice(1,3):tgt_slice(2,3), &
                                                        ibox) = &
                                                        this_cell_ti(:, &
                                                        src_slice(1,1):src_slice(2,1), &
                                                        src_slice(1,2):src_slice(2,2), &
                                                        src_slice(1,3):src_slice(2,3), &
                                                        ibox)
                                                !$OMP END WORKSHARE
                                        END IF
                                        IF (need_atomtypes) THEN
                                                !$OMP WORKSHARE
                                                this_cell_atomtypes(:, &
                                                        tgt_slice(1,1):tgt_slice(2,1), &
                                                        tgt_slice(1,2):tgt_slice(2,2), &
                                                        tgt_slice(1,3):tgt_slice(2,3), &
                                                        ibox) = &
                                                        this_cell_atomtypes(:, &
                                                        src_slice(1,1):src_slice(2,1), &
                                                        src_slice(1,2):src_slice(2,2), &
                                                        src_slice(1,3):src_slice(2,3), &
                                                        ibox)
                                                !$OMP END WORKSHARE
                                        END IF
                                END DO
                        END DO
                END DO
                IF (bitcell_flag) THEN
                        bcd(2) = &
                                ISHFT(lbp32(1),-5)
                        bcd(3) = &
                                bcd(2) * lbp32(2)
                        vlen = bcd(3) * lbp32(3)
                        !$OMP SINGLE
                        box_list(ibox)%bitcell_dimfactor = bcd
                        IF (ALLOCATED(box_list(ibox)%bitcell_int32_vec)) THEN
                                IF (vlen > SIZE(box_list(ibox)%bitcell_int32_vec)) THEN
                                        DEALLOCATE(box_list(ibox)%bitcell_int32_vec, Stat=Allocatestatus)
                                        IF (Allocatestatus /= 0) THEN
                                          err_msg = ''
                                          err_msg(1) = 'Memory could not be deallocated from bitcell_int32_vec'
                                          CALL Clean_Abort(err_msg, 'Sector_Setup')
                                        END IF
                                END IF
                        END IF
                        IF (.NOT. ALLOCATED(box_list(ibox)%bitcell_int32_vec)) THEN
                                ALLOCATE(box_list(ibox)%bitcell_int32_vec(0:vlen-1), Stat=Allocatestatus)
                                IF (Allocatestatus /= 0) THEN
                                  err_msg = ''
                                  err_msg(1) = 'Memory could not be allocated for bitcell_int32_vec'
                                  CALL Clean_Abort(err_msg, 'Sector_Setup')
                                END IF
                                WRITE(*,*) "Allocated bitcell_int32_vec with:"
                                WRITE(*,*) "SIZE", SIZE(box_list(ibox)%bitcell_int32_vec)
                                WRITE(*,*) "Occupying ", SIZEOF(box_list(ibox)%bitcell_int32_vec), " bytes"
                                WRITE(*,*) "on step ", i_mcstep
                        END IF
                        ! bitcell_int8_array_2 shouldn't be necessary, but I sometimes get segfaults from 
                        !    an array supposedly not being allocated if I don't use it as an 
                        !    intermediary. Feel free to remove it and transfer straight from 
                        !    bitcell_int8_array (with the same slicing) if you can do so without segfaults.
                        ALLOCATE(bitcell_int8_array_2(ISHFT(lbp32(1),-3),lbp32(2),lbp32(3)), Stat=Allocatestatus)
                        IF (Allocatestatus /= 0) THEN
                          err_msg = ''
                          err_msg(1) = 'Memory could not be allocated for bitcell_int8_array_2'
                          CALL Clean_Abort(err_msg, 'Sector_Setup')
                        END IF
                        !$OMP END SINGLE
                        !$OMP WORKSHARE
                        bitcell_int8_array_2 = bitcell_int8_array( &
                                4:3+ISHFT(lbp32(1),-3), &
                                32:31+lbp32(2), &
                                32:31+lbp32(3),0)
                        box_list(ibox)%bitcell_int32_vec(0:vlen-1) = TRANSFER(bitcell_int8_array_2, &
                                box_list(ibox)%bitcell_int32_vec)
                        !superpopcnt = SUM(POPCNT(TRANSFER(box_list(ibox)%bitcell_int32_vec(0:vlen-1),bitcell_int64)))
                        !$OMP END WORKSHARE
                        IF (cavity_biasing_flag) THEN
                                int8shape = box_list(ibox)%length_bitcells
                                xub_int64 = INT(box_list(ibox)%length_bitcells(1)-1,INT64)
                                int8shape(1) = ISHFT(int8shape(1),-3)
                                int8ub = int8shape - 1
                                lbp16 = box_list(ibox)%length_bitcells
                                lbp16(1) = IAND(lbp16(1)+15,NOT(15))
                                int16shape = lbp16
                                int16shape(1) = ISHFT(int16shape(1),-4)
                                int16ub = int16shape - 1
                                yi_chunkstride = MIN(MASKR(15)/lbp16(1),lbp16(2))
                                zub = int8ub(3)
                                !$OMP SINGLE
                                big_atom_start = MERGE(1,0,read_atompair_rminsq .OR. calc_rmin_flag)
                                IF (.NOT. ALLOCATED(cavdatalist)) ALLOCATE(cavdatalist(big_atom_start:n_big_atoms,nbr_boxes))
                                DEALLOCATE(bitcell_int8_array_2)
                                ALLOCATE(bitcell_int8_array_2(&
                                        0:ISHFT(lbp16(1),-3)-1,&
                                        0:lbp16(2)-1,&
                                        0:lbp16(3)-1), Stat=Allocatestatus)
                                IF (Allocatestatus /= 0) THEN
                                  err_msg = ''
                                  err_msg(1) = 'Memory could not be allocated for bitcell_int8_array_2'
                                  CALL Clean_Abort(err_msg, 'Sector_Setup')
                                END IF
                                ALLOCATE(bitcell_int16_array(&
                                        0:ISHFT(lbp16(1),-4)-1,&
                                        0:lbp16(2)-1,&
                                        0:lbp16(3)-1), Stat=Allocatestatus)
                                IF (Allocatestatus /= 0) THEN
                                  err_msg = ''
                                  err_msg(1) = 'Memory could not be allocated for bitcell_int16_array'
                                  CALL Clean_Abort(err_msg, 'Sector_Setup')
                                END IF
                                ALLOCATE(zcavcount(&
                                        0:lbp16(3)))
                                !$OMP END SINGLE
                                DO i_big_atom = big_atom_start, n_big_atoms
                                        !$OMP WORKSHARE
                                        bitcell_int8_array_2 = 0_INT8
                                        bitcell_int8_array_2(:int8ub(1),:,:) = &
                                                NOT(bitcell_int8_array( &
                                                4:3+int8shape(1), &
                                                32:31+int8shape(2), &
                                                32:31+int8shape(3),i_big_atom))
                                        bitcell_int16_array = &
                                                RESHAPE(&
                                                TRANSFER(bitcell_int8_array_2,bitcell_int16_array),&
                                                SHAPE(bitcell_int16_array))
                                        !zcavcount(1:) = SUM(INT(SUM(POPCNT(bitcell_int16_array),1),INT64),1)
                                        !zcavcount(0) = 0_INT64
                                        !$OMP END WORKSHARE
                                        !$OMP DO SCHEDULE(STATIC)
                                        DO zi = 0, zub
                                                ncavs = 0_INT64
                                                DO yi_chunkstart = 0, int16ub(2), yi_chunkstride
                                                        yi_chunkend = MIN(yi_chunkstart + yi_chunkstride - 1, int16ub(2))
                                                        chunksize = (yi_chunkend + 1 - yi_chunkstart) * int16shape(1)
                                                        bitcell_int16_vec(1:chunksize) = &
                                                                TRANSFER(bitcell_int16_array(:,yi_chunkstart:yi_chunkend,zi),bitcell_int16_vec)
                                                        chunksize_p16 = IAND(chunksize+15,NOT(15))
                                                        IF (chunksize_p16 > chunksize) THEN
                                                                bitcell_int16_vec(chunksize+1:chunksize_p16) = 0_INT16
                                                        END IF
                                                        ncavs_int16 = 0_INT16
                                                        !DIR$ ASSUME (MOD(chunksize_p16,16) .EQ. 0)
                                                        !DIR$ VECTOR ALIGNED
                                                        !$OMP SIMD PRIVATE(cavbits_int16) REDUCTION(+:ncavs_int16)
                                                        DO i = 1, chunksize_p16
                                                                cavbits_int16 = bitcell_int16_vec(i)
                                                                ncavs_int16 = ncavs_int16 + IAND(cavbits_int16,1_INT16)
                                                                DO j = 1, 15
                                                                        ncavs_int16 = ncavs_int16 + &
                                                                                IAND(SHIFTR(cavbits_int16,j),1_INT16)
                                                                END DO
                                                        END DO
                                                        !$OMP END SIMD
                                                        ncavs = ncavs + INT(ncavs_int16,INT64)
                                                END DO
                                                !DO yi = 0, int16ub(2)
                                                !        ncavs_int16 = 0_INT16
                                                !        !DIR$ VECTOR ALIGNED
                                                !        !$OMP SIMD PRIVATE(oneshift_int16) REDUCTION(+:ncavs_int16)
                                                !        DO i = 0, 15
                                                !                oneshift_int16 = oneshift_int16_vec(i) !ISHFT(1_INT16,i)
                                                !                DO j = 0, int16ub(1)
                                                !                        IF (IAND(oneshift_int16,&
                                                !                                bitcell_int16_array(j,yi,zi)) .NE. 0_INT16) THEN
                                                !                                ncavs_int16 = ncavs_int16 + 1
                                                !                        END IF
                                                !                END DO
                                                !        END DO
                                                !        !$OMP END SIMD
                                                !        ncavs = ncavs + ncavs_int16
                                                !END DO
                                                zcavcount(zi+1) = ncavs
                                        END DO
                                        !$OMP END DO
                                        !$OMP SINGLE
                                        zcavcount(0) = 0_INT64
                                        ncavs = zcavcount(1)
                                        DO zi = 2, lbp16(3)
                                                ncavs = ncavs + zcavcount(zi)
                                                zcavcount(zi) = ncavs
                                        END DO
                                        IF (ALLOCATED(cavdatalist(i_big_atom,ibox)%cavity_locs)) THEN
                                                IF (ncavs>SIZE(cavdatalist(i_big_atom,ibox)%cavity_locs,1,INT64)) THEN
                                                        DEALLOCATE(cavdatalist(i_big_atom,ibox)%cavity_locs)
                                                END IF
                                        END IF
                                        IF (.NOT. ALLOCATED(cavdatalist(i_big_atom,ibox)%cavity_locs)) THEN
                                                ALLOCATE(cavdatalist(i_big_atom,ibox)%cavity_locs(0:ncavs-1))
                                        END IF
                                        cavdatalist(i_big_atom,ibox)%ncavs = ncavs
                                        cavdatalist(i_big_atom,ibox)%ncavs_dp = REAL(ncavs,DP)
                                        cavdatalist(i_big_atom,ibox)%ln_cavfrac = LOG(REAL(ncavs,DP)/&
                                                PRODUCT(REAL(box_list(ibox)%length_bitcells,DP)))
                                        !$OMP END SINGLE
                                        !$OMP DO SCHEDULE(STATIC)
                                        DO zi = 0, zub
                                                icav = zcavcount(zi)
                                                locbase = ISHFT(INT(zi,INT64),42)
                                                DO yi = 0, int16ub(2)
                                                        !DIR$ VECTOR ALIGNED
                                                        !$OMP SIMD PRIVATE(oneshift_int16)
                                                        DO i = 0, 15
                                                                oneshift_int16 = oneshift_int16_vec(i) !ISHFT(1_INT16,i)
                                                                DO j = 0, int16ub(1)
                                                                        ncavlvec(j*16+i) = &
                                                                                IAND(oneshift_int16,&
                                                                                bitcell_int16_array(j,yi,zi)) .EQ. 0_INT16
                                                                END DO
                                                        END DO
                                                        !$OMP END SIMD
                                                        DO xi_int64 = 0_INT64, xub_int64
                                                                IF (ncavlvec(xi_int64)) CYCLE
                                                                cavdatalist(i_big_atom,ibox)%cavity_locs(icav) = &
                                                                        IOR(locbase,xi_int64)
                                                                icav = icav + 1_INT64
                                                        END DO
                                                        locbase = locbase + ISHFT(1_INT64,21)
                                                END DO
                                        END DO
                                        !$OMP END DO
                                END DO
                        END IF
                        !$OMP SINGLE
                        IF (l_firstframe) THEN
                                DO isolvent = 1, solvents_or_types_maxind
                                      DO i_big_atom = 0, n_big_atoms
                                              OPEN(2777,file=TRIM(TRIM(run_name) // &
                                                      '.bitcell_int64.' // &
                                                      TRIM(Int_To_String(isolvent)) // '.' // &
                                                      TRIM(Int_To_String(i_big_atom))))
                                              DO zi = -28, 28
                                                      WRITE(2777,*) bitcell_int64(-28:28,zi,isolvent,i_big_atom,1)
                                              END DO
                                              CLOSE(2777)
                                      END DO
                                END DO
                                IF (cavity_biasing_flag) THEN
                                        DO i_big_atom = big_atom_start, n_big_atoms
                                                OPEN(2777,file=TRIM(TRIM(run_name) //&
                                                        '.cavity_locs.' // TRIM(Int_To_String(i_big_atom))))
                                                DO icav = 0, cavdatalist(i_big_atom,1)%ncavs
                                                        WRITE(2777,*) cavdatalist(i_big_atom,1)%cavity_locs(icav)
                                                END DO
                                                CLOSE(2777)
                                        END DO
                                END IF
                                OPEN(2777,file=TRIM(TRIM(run_name) // '.bitcell_info'))
                                WRITE(2777,*) "length_bitcells"
                                WRITE(2777,*) box_list(1)%length_bitcells
                                WRITE(2777,*) "bitcell_face_distance"
                                WRITE(2777,*) box_list(1)%bitcell_face_distance
                                WRITE(2777,*) "setbit_extent"
                                WRITE(2777,*) box_list(1)%setbit_extent
                                WRITE(2777,*) "bitcell_xyzortho_bbox_length"
                                WRITE(2777,*) bitcell_xyzortho_bbox_length
                                WRITE(2777,*) "bitcell_H"
                                WRITE(2777,*) bitcell_H
                                WRITE(2777,*) "lbp32"
                                WRITE(2777,*) lbp32
                                WRITE(2777,*) "bcd"
                                WRITE(2777,*) bcd
                                WRITE(2777,*) "n_big_atoms"
                                WRITE(2777,*) n_big_atoms
                                WRITE(2777,*) "sbe_ti_mat"
                                WRITE(2777,*) sbe_ti_mat
                                WRITE(2777,*) "cavfrac"
                                WRITE(2777,*) EXP(cavdatalist(:,1)%ln_cavfrac)
                                WRITE(2777,*) "ncavs"
                                WRITE(2777,*) cavdatalist(:,1)%ncavs
                                CLOSE(2777)
                                OPEN(2777,file=TRIM(TRIM(run_name) // '.bitcell_int32_vec'))
                                DO i = 0, vlen-1
                                        WRITE(2777,*) box_list(ibox)%bitcell_int32_vec(i)
                                END DO
                                CLOSE(2777)
                        END IF
                        DEALLOCATE(bitcell_int8_array, Stat=Allocatestatus)
                        IF (Allocatestatus /= 0) THEN
                          err_msg = ''
                          err_msg(1) = 'Memory could not be deallocated from bitcell_int8_array'
                          CALL Clean_Abort(err_msg, 'Sector_Setup')
                        END IF
                        DEALLOCATE(bitcell_int8_array_2, Stat=Allocatestatus)
                        IF (Allocatestatus /= 0) THEN
                          err_msg = ''
                          err_msg(1) = 'Memory could not be deallocated from bitcell_int8_array_2'
                          CALL Clean_Abort(err_msg, 'Sector_Setup')
                        END IF
                        !$OMP END SINGLE
                        !WRITE(*,*) "DEALLOCATING bitcell_int8_array"
                        !DEALLOCATE(bitcell_int8_array)
                        !WRITE(*,*) "DEALLOCATED bitcell_int8_array", ALLOCATED(bitcell_int8_array)
                        !!$OMP END SINGLE
                        !!$OMP BARRIER
                        !!$OMP WORKSHARE
                        !box_list(ibox)%bitcell_int32_vec(0:vlen-1) = TRANSFER(bitcell_int8_array( &
                        !        4:3+ISHFT(lbp32(1),-3), &
                        !        32:31+lbp32(2), &
                        !        32:31+lbp32(3)), &
                        !        box_list(ibox)%bitcell_int32_vec)
                        !!$OMP END WORKSHARE
                        !!$OMP BARRIER
                        !!$OMP SINGLE
                        !WRITE(*,*) "DEALLOCATING bitcell_int8_array"
                        !DEALLOCATE(bitcell_int8_array)
                        !WRITE(*,*) "DEALLOCATED bitcell_int8_array", ALLOCATED(bitcell_int8_array)
                        !!$OMP END SINGLE
                END IF
                !!!!$OMP END DO
                !!!!$OMP END PARALLEL
          END DO
          IF (cbmc_cell_list_flag) THEN
                  DO ibox = 1, nbr_boxes
                          bt = box_list(ibox)%border_thickness
                          zub = box_list(ibox)%sectorbound(3)
                          yub = box_list(ibox)%sectorbound(2)
                          xub = box_list(ibox)%sectorbound(1)
                          zlb = -zub
                          ylb = -yub
                          xlb = -xub
                          !$OMP DO COLLAPSE(3) REDUCTION(MAX:max_neighbors)
                          DO zi = zlb, zub
                                  DO yi = ylb, yub
                                          DO xi = xlb, xub
                                                  max_neighbors = MAX(SUM(IAND(n_cell_atoms( &
                                                          xi-bt(1):xi+bt(1), &
                                                          yi-bt(2):yi+bt(2), &
                                                          zi-bt(3):zi+bt(3), &
                                                          ibox), &
                                                          box_list(ibox)%cbmc_cell_mask)), &
                                                          max_neighbors)
                                          END DO
                                  END DO
                          END DO
                          !$OMP END DO
                  END DO
                  !$OMP SINGLE
                  IF (ALLOCATED(cbmc_cell_n_interact)) THEN
                          IF (ANY(adj_cellmaxbound > adj_cellmaxbound_old)) THEN
                                  DEALLOCATE(cbmc_cell_n_interact, Stat=AllocateStatus)
                                  IF (Allocatestatus /= 0) THEN
                                    err_msg = ''
                                    err_msg(1) = 'Memory could not be deallocated from cbmc_cell_n_interact'
                                    CALL Clean_Abort(err_msg, 'Sector_Setup')
                                  END IF
                          END IF
                  END IF
                  IF (.NOT. ALLOCATED(cbmc_cell_n_interact)) THEN
                          ALLOCATE(cbmc_cell_n_interact( &
                                  -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                                  -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                                  -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                                  nbr_boxes), Stat=AllocateStatus)
                          IF (Allocatestatus /= 0) THEN
                            err_msg = ''
                            err_msg(1) = 'Memory could not be allocated for cbmc_cell_n_interact'
                            CALL Clean_Abort(err_msg, 'Sector_Setup')
                          END IF
                          WRITE(*,*) "Allocated cbmc_cell_n_interact with:"
                          WRITE(*,*) "SHAPE", SHAPE(cbmc_cell_n_interact)
                          WRITE(*,*) "Occupying ", SIZEOF(cbmc_cell_n_interact), " bytes"
                          WRITE(*,*) "on step ", i_mcstep
                  END IF
                  ALLOCATE(cell_l_inrange(max_neighbors,2, &
                          -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                          -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                          -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                          nbr_boxes), Stat=AllocateStatus)
                  cbmc_max_interact_old = cbmc_max_interact
                  cbmc_max_interact = 0
                  !$OMP END SINGLE
          ELSE
                  DO ibox = 1, nbr_boxes
                          zub = box_list(ibox)%sectorbound(3)
                          yub = box_list(ibox)%sectorbound(2)
                          xub = box_list(ibox)%sectorbound(1)
                          zlb = -zub
                          ylb = -yub
                          xlb = -xub
                          !$OMP DO COLLAPSE(3) REDUCTION(MAX:max_neighbors)
                          DO zi = zlb, zub
                                  DO yi = ylb, yub
                                          DO xi = xlb, xub
                                                  max_neighbors = MAX(SUM(n_cell_atoms( &
                                                          xi-1:xi+1, &
                                                          yi-1:yi+1, &
                                                          zi-1:zi+1, &
                                                          ibox)), &
                                                          max_neighbors)
                                          END DO
                                  END DO
                          END DO
                          !$OMP END DO
                  END DO
                  !$OMP SINGLE
                  ALLOCATE(cell_l_inrange(max_neighbors,1, &
                          -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                          -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                          -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                          nbr_boxes))
                  !$OMP END SINGLE
          END IF
          !$OMP WORKSHARE
          cell_l_inrange = .FALSE.
          !$OMP END WORKSHARE



          DO ibox = 1, nbr_boxes
                  l_ortho = box_list(ibox)%int_box_shape <= int_ortho
                  rcutsq = REAL(rcut_cbmcsq(ibox),SP)
                  hl = REAL(box_list(ibox)%cell_xyzortho_bbox_length*0.5_DP,SP)
                  bt = box_list(ibox)%border_thickness
                  IF (.NOT. l_ortho) THEN
                          h11 = REAL(box_list(ibox)%length(1,1),SP)
                          h12 = REAL(box_list(ibox)%length(1,2),SP)
                          h22 = REAL(box_list(ibox)%length(2,2),SP)
                          h13 = REAL(box_list(ibox)%length(1,3),SP)
                          h23 = REAL(box_list(ibox)%length(2,3),SP)
                          h33 = REAL(box_list(ibox)%length(3,3),SP)
                          !length_cells_recip = 1.0/box_list(ibox)%real_length_cells
                          !hlcr = length_cells_recip*0.5
                          !length_cells_recip_dp = 1.0_DP/REAL(box_list(ibox)%length_cells,DP)
                          !hlcr_dp = length_cells_recip_dp*0.5_DP
                          cell_H = box_list(ibox)%cell_H_sp
                          zub = bt(3)+box_list(ibox)%sectorbound(3)
                          yub = bt(2)+box_list(ibox)%sectorbound(2)
                          xub = bt(1)+box_list(ibox)%sectorbound(1)
                          zlb = -zub
                          ylb = -yub
                          xlb = -xub
                          !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
                          DO zi = zlb, zub
                                  DO yi = ylb, yub
                                          DO xi = xlb, xub
                                                  DO i = 1, n_cell_atoms(xi,yi,zi,ibox)
                                                          isp = this_cell_rsp(i,1,xi,yi,zi,ibox)
                                                          rxp = h11*isp
                                                          isp = this_cell_rsp(i,2,xi,yi,zi,ibox)
                                                          rxp = rxp + h12*isp
                                                          ryp = h22*isp
                                                          isp = this_cell_rsp(i,3,xi,yi,zi,ibox)
                                                          rxp = rxp + h13*isp
                                                          ryp = ryp + h23*isp
                                                          rzp = h33*isp
                                                          this_cell_rsp(i,1,xi,yi,zi,ibox) = rxp
                                                          this_cell_rsp(i,2,xi,yi,zi,ibox) = ryp
                                                          this_cell_rsp(i,3,xi,yi,zi,ibox) = rzp
                                                  END DO
                                          END DO
                                  END DO
                          END DO
                          !$OMP END DO
                  END IF
                  !!$OMP SINGLE
                  !WRITE(*,*) SHAPE(cell_l_inrange)
                  !!$OMP END SINGLE
                  zub = box_list(ibox)%sectorbound(3)
                  yub = box_list(ibox)%sectorbound(2)
                  xub = box_list(ibox)%sectorbound(1)
                  zlb = -zub
                  ylb = -yub
                  xlb = -xub

                  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC) REDUCTION(MAX:cbmc_max_interact,max_adj_cell_atoms)
                  DO zi = zlb, zub
                  DO yi = ylb, yub
                  DO xi = xlb, xub
                        istart = 1
                        iend = 0
                        DO dzi = -1, 1
                                zi2 = zi+dzi
                                DO dyi = -1, 1
                                        yi2 = yi + dyi
                                        DO dxi = -1, 1
                                                xi2 = xi + dxi
                                                vlen = n_cell_atoms(xi2,yi2,zi2,ibox)
                                                IF (vlen < 1) CYCLE
                                                iend = istart + vlen - 1
                                                cbmc_cell_rsp_priv(istart:iend,1:3) = &
                                                        this_cell_rsp(1:vlen,1:3,xi2,yi2,zi2,ibox)
                                                IF (read_atompair_rminsq) THEN
                                                        ti_priv(istart:iend) = &
                                                                this_cell_ti(1:vlen,xi2,yi2,zi2,ibox)
                                                ELSE IF (calc_rmin_flag) THEN
                                                        ti_priv(istart:iend) = &
                                                                this_cell_atomtypes(1:vlen,xi2,yi2,zi2,ibox)
                                                END IF
                                                istart = istart + vlen
                                        END DO
                                END DO
                        END DO
                        adj_iend = iend
                        DO dzi = -bt(3), bt(3)
                                IF (ABS(dzi)<2) CYCLE
                                zi2 = zi+dzi
                                DO dyi = -bt(2), bt(2)
                                        IF (ABS(dyi)<2) CYCLE
                                        yi2 = yi + dyi
                                        DO dxi = -bt(1), bt(1)
                                                IF (ABS(dxi)<2) CYCLE
                                                xi2 = xi + dxi
                                                vlen = IAND(n_cell_atoms(xi2,yi2,zi2,ibox), &
                                                        box_list(ibox)%cbmc_cell_mask(dxi,dyi,dzi))
                                                IF (vlen < 1) CYCLE
                                                iend = istart + vlen - 1
                                                cbmc_cell_rsp_priv(istart:iend,1:3) = &
                                                        this_cell_rsp(1:vlen,1:3,xi2,yi2,zi2,ibox)
                                                istart = istart + vlen
                                        END DO
                                END DO
                        END DO
                        xyzi_sp = REAL((/ xi, yi, zi /),SP)
                        IF (l_ortho) THEN
                                cell_rp = xyzi_sp*box_list(ibox)%cell_H_diag
                        ELSE
                                cell_rp = MATMUL(cell_H,xyzi_sp)
                        END IF
                                !cell_rsxp = xi*box_list(ibox)%cell_H_diag(1)
                                !cell_rsyp = yi*box_list(ibox)%cell_H_diag(2)
                                !cell_rszp = zi*box_list(ibox)%cell_H_diag(3)
                        !DIR$ VECTOR ALIGNED
                        DO i = 1, iend
                                drxp = cbmc_cell_rsp_priv(i,1)-cell_rp(1)
                                drxp = MAX(ABS(drxp)-hl(1),0.0)
                                dryp = cbmc_cell_rsp_priv(i,2)-cell_rp(2)
                                dryp = MAX(ABS(dryp)-hl(2),0.0)
                                drzp = cbmc_cell_rsp_priv(i,3)-cell_rp(3)
                                drzp = MAX(ABS(drzp)-hl(3),0.0)
                                rsq_vec(i) = drxp*drxp + dryp*dryp + drzp*drzp
                        END DO
                                !DO i = 1, iend
                                !        dsp = cbmc_cell_rsp_priv(i,1)-cell_rsxp
                                !        dsp = SIGN(MAX(ABS(dsp)-hlcr(1),0.0),dsp)
                                !        drxp = h11*dsp
                                !        dryp = h21*dsp
                                !        drzp = h31*dsp
                                !        dsp = cbmc_cell_rsp_priv(i,2)-cell_rsyp
                                !        dsp = SIGN(MAX(ABS(dsp)-hlcr(2),0.0),dsp)
                                !        drxp = drxp + h12*dsp
                                !        dryp = dryp + h22*dsp
                                !        drzp = drzp + h32*dsp
                                !        dsp = cbmc_cell_rsp_priv(i,3)-cell_rszp
                                !        dsp = SIGN(MAX(ABS(dsp)-hlcr(3),0.0),dsp)
                                !        drxp = drxp + h13*dsp
                                !        dryp = dryp + h23*dsp
                                !        drzp = drzp + h33*dsp
                                !        drxp = drxp*drxp
                                !        drxp = drxp + dryp*dryp
                                !        drxp = drxp + drzp*drzp
                                !        rsq_vec(i) = drxp
                                !END DO
                        IF (read_atompair_rminsq) THEN
                                cell_l_inrange(1:adj_iend,1,xi,yi,zi,ibox) = &
                                        rsq_vec(1:adj_iend) < solvent_max_rminsq_sp(ti_priv(1:adj_iend),ibox)
                        ELSE IF (calc_rmin_flag) THEN
                                !IF (ANY( (/ adj_iend,1,xi,yi,zi,ibox /)>UBOUND(cell_l_inrange) .OR. &
                                !        (/ adj_iend,1,xi,yi,zi,ibox /)<LBOUND(cell_l_inrange))) THEN
                                !        WRITE(*,*) adj_iend,1,xi,yi,zi,ibox
                                !ELSE IF (adj_iend > UBOUND(ti_priv,1) .OR. adj_iend < LBOUND(ti_priv,1)) THEN
                                !        WRITE(*,*) "adj_iend is out of bounds", adj_iend, UBOUND(ti_priv,1)
                                !ELSE IF (ANY(ti_priv(1:adj_iend)>UBOUND(atomtype_max_rminsq_sp,1) .OR. &
                                !        ti_priv(1:adj_iend)<LBOUND(atomtype_max_rminsq_sp,1))) THEN
                                !        WRITE(*,*) "ti_priv is out of bounds", ti_priv(1:adj_iend)
                                !ELSE IF (adj_iend > UBOUND(rsq_vec,1) ) THEN
                                !        WRITE(*,*) "adj_iend is out of bounds for rsq_vec", adj_iend, UBOUND(rsq_vec)
                                !END IF
                                cell_l_inrange(1:adj_iend,1,xi,yi,zi,ibox) = &
                                        rsq_vec(1:adj_iend) < atomtype_max_rminsq_sp(ti_priv(1:adj_iend))
                        ELSE
                                cell_l_inrange(1:adj_iend,1,xi,yi,zi,ibox) = &
                                        rsq_vec(1:adj_iend) < sp_rcut_lowsq
                        END IF
                        IF (cbmc_cell_list_flag) THEN
                                cell_l_inrange(1:iend,2,xi,yi,zi,ibox) = rsq_vec(1:iend) < rcutsq
                                cbmc_max_interact = MAX(cbmc_max_interact,COUNT(cell_l_inrange(1:iend,2,xi,yi,zi,ibox)))
                        END IF
                        max_adj_cell_atoms = MAX(max_adj_cell_atoms,COUNT(cell_l_inrange(1:adj_iend,1,xi,yi,zi,ibox)))
                  END DO
                  END DO
                  END DO
                  !$OMP END DO
          END DO
          IF (cbmc_cell_list_flag) THEN
                  !$OMP WORKSHARE
                  cell_l_inrange(:,2,:,:,:,:) = cell_l_inrange(:,2,:,:,:,:) .AND. .NOT. cell_l_inrange(:,1,:,:,:,:)
                  !$OMP END WORKSHARE
          END IF
          !$OMP SINGLE
          IF (cbmc_cell_list_flag) THEN
                  cbmc_max_interact = MAX(cbmc_max_interact,cbmc_max_interact_old)
                  ! ensure 32-byte padding of first dimension
                  IF (MOD(cbmc_max_interact,8) .NE. 0) cbmc_max_interact = (cbmc_max_interact/8+1)*8
                  IF (cbmc_max_interact > cbmc_max_interact_old .OR. ANY(adj_cellmaxbound > adj_cellmaxbound_old)) THEN
                          IF (ALLOCATED(cbmc_cell_rsp)) THEN
                                  WRITE(*,*) "cbmc_cell_rsp shape was previously: "
                                  WRITE(*,*) SHAPE(cbmc_cell_rsp)
                                  WRITE(*,*) cbmc_max_interact, ">", cbmc_max_interact_old
                                  WRITE(*,*) adj_cellmaxbound
                                  WRITE(*,*) adj_cellmaxbound_old
                                  DEALLOCATE(cbmc_cell_rsp, Stat=AllocateStatus)
                          END IF
                          IF (Allocatestatus /= 0) THEN
                            err_msg = ''
                            err_msg(1) = 'Memory could not be deallocated from cbmc_cell_rsp'
                            CALL Clean_Abort(err_msg, 'Sector_Setup')
                          END IF
                          IF (ALLOCATED(cbmc_cell_atomtypes)) DEALLOCATE(cbmc_cell_atomtypes, Stat=AllocateStatus)
                          IF (Allocatestatus /= 0) THEN
                            err_msg = ''
                            err_msg(1) = 'Memory could not be deallocated from cbmc_cell_atomtypes'
                            CALL Clean_Abort(err_msg, 'Sector_Setup')
                          END IF
                          IF (ALLOCATED(cbmc_cell_ti)) DEALLOCATE(cbmc_cell_ti, Stat=AllocateStatus)
                          IF (Allocatestatus /= 0) THEN
                            err_msg = ''
                            err_msg(1) = 'Memory could not be deallocated from cbmc_cell_ti'
                            CALL Clean_Abort(err_msg, 'Sector_Setup')
                          END IF
                  END IF
                  IF (ALLOCATED(cbmc_cell_rsp)) THEN
                          IF (cbmc_max_interact > cbmc_max_interact_old .OR. ANY(adj_cellmaxbound > adj_cellmaxbound_old)) THEN
                                  DEALLOCATE(cbmc_cell_rsp, Stat=AllocateStatus)
                                  IF (Allocatestatus /= 0) THEN
                                    err_msg = ''
                                    err_msg(1) = 'Memory could not be deallocated from cbmc_cell_rsp'
                                    CALL Clean_Abort(err_msg, 'Sector_Setup')
                                  END IF
                          END IF
                  END IF
                  IF (.NOT. ALLOCATED(cbmc_cell_rsp)) THEN
                          ALLOCATE(cbmc_cell_rsp(cbmc_max_interact,4, &
                                  -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                                  -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                                  -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                                  nbr_boxes), Stat=AllocateStatus)
                          IF (Allocatestatus /= 0) THEN
                            err_msg = ''
                            err_msg(1) = 'Memory could not be allocated for cbmc_cell_rsp'
                            CALL Clean_Abort(err_msg, 'Sector_Setup')
                          END IF
                          WRITE(*,*) "Allocated cbmc_cell_rsp with:"
                          WRITE(*,*) "SHAPE", SHAPE(cbmc_cell_rsp)
                          WRITE(*,*) "Occupying ", SIZEOF(cbmc_cell_rsp), " bytes"
                          WRITE(*,*) "on step ", i_mcstep
                  END IF
                  IF (need_atom_ti) THEN
                          IF (.NOT. ALLOCATED(cbmc_cell_ti)) THEN
                                  ALLOCATE(cbmc_cell_ti(cbmc_max_interact, &
                                          -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                                          -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                                          -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                                          nbr_boxes), Stat=AllocateStatus)
                                  IF (Allocatestatus /= 0) THEN
                                    err_msg = ''
                                    err_msg(1) = 'Memory could not be allocated for cbmc_cell_ti'
                                    CALL Clean_Abort(err_msg, 'Sector_Setup')
                                  END IF
                                  WRITE(*,*) "Allocated cbmc_cell_ti with:"
                                  WRITE(*,*) "SHAPE", SHAPE(cbmc_cell_ti)
                                  WRITE(*,*) "Occupying ", SIZEOF(cbmc_cell_ti), " bytes"
                                  WRITE(*,*) "on step ", i_mcstep
                          END IF
                  END IF
                  IF (need_atomtypes) THEN
                          IF (.NOT. ALLOCATED(cbmc_cell_atomtypes)) THEN
                                  ALLOCATE(cbmc_cell_atomtypes(cbmc_max_interact, &
                                          -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                                          -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                                          -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                                          nbr_boxes), Stat=AllocateStatus)
                                  IF (Allocatestatus /= 0) THEN
                                    err_msg = ''
                                    err_msg(1) = 'Memory could not be allocated for cbmc_cell_atomtypes'
                                    CALL Clean_Abort(err_msg, 'Sector_Setup')
                                  END IF
                                  WRITE(*,*) "Allocated cbmc_cell_atomtypes with:"
                                  WRITE(*,*) "SHAPE", SHAPE(cbmc_cell_atomtypes)
                                  WRITE(*,*) "Occupying ", SIZEOF(cbmc_cell_atomtypes), " bytes"
                                  WRITE(*,*) "on step ", i_mcstep
                          END IF
                  END IF
          ELSE
                  ! ensure 32-byte padding of first dimension
                  IF (MOD(max_adj_cell_atoms,8) .NE. 0) max_adj_cell_atoms = (max_adj_cell_atoms/8+1)*8
                  IF (ALLOCATED(cbmc_cell_rsp)) THEN
                          IF (SIZE(cbmc_cell_rsp,1) < max_adj_cell_atoms .OR. ANY(adj_cellmaxbound>adj_cellmaxbound_old)) THEN
                                  DEALLOCATE(cbmc_cell_rsp, Stat=AllocateStatus)
                                  IF (Allocatestatus /= 0) THEN
                                    err_msg = ''
                                    err_msg(1) = 'Memory could not be deallocated from cbmc_cell_rsp'
                                    CALL Clean_Abort(err_msg, 'Sector_Setup')
                                  END IF
                                  IF (need_atom_ti) THEN
                                          DEALLOCATE(cbmc_cell_ti, Stat=AllocateStatus)
                                          IF (Allocatestatus /= 0) THEN
                                            err_msg = ''
                                            err_msg(1) = 'Memory could not be deallocated from cbmc_cell_ti'
                                            CALL Clean_Abort(err_msg, 'Sector_Setup')
                                          END IF
                                  END IF
                                  IF (need_atomtypes) THEN
                                          DEALLOCATE(cbmc_cell_atomtypes, Stat=AllocateStatus)
                                          IF (Allocatestatus /= 0) THEN
                                            err_msg = ''
                                            err_msg(1) = 'Memory could not be deallocated from cbmc_cell_atomtypes'
                                            CALL Clean_Abort(err_msg, 'Sector_Setup')
                                          END IF
                                  END IF
                          END IF
                  END IF
                  IF (.NOT. ALLOCATED(cbmc_cell_rsp)) THEN
                          ALLOCATE(cbmc_cell_rsp(max_adj_cell_atoms,4, &
                                  -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                                  -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                                  -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                                  nbr_boxes), Stat=AllocateStatus)
                          IF (Allocatestatus /= 0) THEN
                            err_msg = ''
                            err_msg(1) = 'Memory could not be allocated for cbmc_cell_rsp'
                            CALL Clean_Abort(err_msg, 'Sector_Setup')
                          END IF
                          WRITE(*,*) "Allocated cbmc_cell_rsp with:"
                          WRITE(*,*) "SHAPE", SHAPE(cbmc_cell_rsp)
                          WRITE(*,*) "Occupying ", SIZEOF(cbmc_cell_rsp), " bytes"
                          WRITE(*,*) "on step ", i_mcstep
                          IF (need_atom_ti) THEN
                                  ALLOCATE(cbmc_cell_ti(max_adj_cell_atoms, &
                                          -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                                          -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                                          -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                                          nbr_boxes), Stat=AllocateStatus)
                                  IF (Allocatestatus /= 0) THEN
                                    err_msg = ''
                                    err_msg(1) = 'Memory could not be allocated for cbmc_cell_ti'
                                    CALL Clean_Abort(err_msg, 'Sector_Setup')
                                  END IF
                                  WRITE(*,*) "Allocated cbmc_cell_ti with:"
                                  WRITE(*,*) "SHAPE", SHAPE(cbmc_cell_ti)
                                  WRITE(*,*) "Occupying ", SIZEOF(cbmc_cell_ti), " bytes"
                                  WRITE(*,*) "on step ", i_mcstep
                          END IF
                          IF (need_atomtypes) THEN
                                  ALLOCATE(cbmc_cell_atomtypes(max_adj_cell_atoms, &
                                          -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                                          -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                                          -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                                          nbr_boxes), Stat=AllocateStatus)
                                  IF (Allocatestatus /= 0) THEN
                                    err_msg = ''
                                    err_msg(1) = 'Memory could not be allocated for cbmc_cell_atomtypes'
                                    CALL Clean_Abort(err_msg, 'Sector_Setup')
                                  END IF
                                  WRITE(*,*) "Allocated cbmc_cell_atomtypes with:"
                                  WRITE(*,*) "SHAPE", SHAPE(cbmc_cell_atomtypes)
                                  WRITE(*,*) "Occupying ", SIZEOF(cbmc_cell_atomtypes), " bytes"
                                  WRITE(*,*) "on step ", i_mcstep
                          END IF
                  END IF
          END IF
          l_firstframe = .FALSE.
          !$OMP END SINGLE
          !$OMP WORKSHARE
          cbmc_cell_rsp(:,1:3,:,:,:,:) = infinity_sp ! guaranteed to be out of range unless overwritten
          cbmc_cell_rsp(:,4,:,:,:,:) = 0.0
          !$OMP END WORKSHARE NOWAIT
          IF (need_atom_ti) THEN
                  !$OMP WORKSHARE
                  cbmc_cell_ti = 1
                  !$OMP END WORKSHARE NOWAIT
          END IF
          IF (need_atomtypes) THEN
                  !$OMP WORKSHARE
                  cbmc_cell_atomtypes = 0
                  !$OMP END WORKSHARE NOWAIT
          END IF
          !$OMP BARRIER
          DO ibox = 1, nbr_boxes
                  zub = box_list(ibox)%sectorbound(3)
                  yub = box_list(ibox)%sectorbound(2)
                  xub = box_list(ibox)%sectorbound(1)
                  zlb = -zub
                  ylb = -yub
                  xlb = -xub
                  !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
                  DO zi = zlb, zub
                  DO yi = ylb, yub
                  DO xi = xlb, xub
                        istart = 1
                        iend = 0
                        DO dzi = -1, 1
                                zi2 = zi+dzi
                                DO dyi = -1, 1
                                        yi2 = yi + dyi
                                        DO dxi = -1, 1
                                                xi2 = xi + dxi
                                                vlen = n_cell_atoms(xi2,yi2,zi2,ibox)
                                                IF (vlen < 1) CYCLE
                                                iend = istart + vlen - 1
                                                cbmc_cell_rsp_priv(istart:iend,1:4) = &
                                                        this_cell_rsp(1:vlen,1:4,xi2,yi2,zi2,ibox)
                                                IF (need_atom_ti) THEN
                                                        ti_priv(istart:iend) = &
                                                                this_cell_ti(1:vlen,xi2,yi2,zi2,ibox)
                                                END IF
                                                IF (need_atomtypes) THEN
                                                        atomtype_priv(istart:iend) = &
                                                                this_cell_atomtypes(1:vlen,xi2,yi2,zi2,ibox)
                                                END IF
                                                istart = istart + vlen
                                        END DO
                                END DO
                        END DO
                        n_interact = 0
                        DO i = 1, iend
                                IF (cell_l_inrange(i,1,xi,yi,zi,ibox)) THEN
                                        n_interact = n_interact + 1
                                        which_interact(n_interact) = i
                                END IF
                        END DO
                        !n_interact = IAND(n_interact+7,pad8mask) ! padding to multiple of 8 or zero
                        n_adj_cell_atoms(xi,yi,zi,ibox) = IAND(n_interact+7,NOT(7))
                        IF (cbmc_cell_list_flag) THEN
                                DO dzi = -bt(3), bt(3)
                                        IF (ABS(dzi)<2) CYCLE
                                        zi2 = zi+dzi
                                        DO dyi = -bt(2), bt(2)
                                                IF (ABS(dyi)<2) CYCLE
                                                yi2 = yi + dyi
                                                DO dxi = -bt(1), bt(1)
                                                        IF (ABS(dxi)<2) CYCLE
                                                        xi2 = xi + dxi
                                                        vlen = IAND(n_cell_atoms(xi2,yi2,zi2,ibox), &
                                                                box_list(ibox)%cbmc_cell_mask(dxi,dyi,dzi))
                                                        IF (vlen < 1) CYCLE
                                                        iend = istart + vlen - 1
                                                        cbmc_cell_rsp_priv(istart:iend,1:4) = &
                                                                this_cell_rsp(1:vlen,1:4,xi2,yi2,zi2,ibox)
                                                        IF (need_atom_ti) THEN
                                                                ti_priv(istart:iend) = &
                                                                        this_cell_ti(1:vlen,xi2,yi2,zi2,ibox)
                                                        END IF
                                                        IF (need_atomtypes) THEN
                                                                atomtype_priv(istart:iend) = &
                                                                        this_cell_atomtypes(1:vlen,xi2,yi2,zi2,ibox)
                                                        END IF
                                                        istart = istart + vlen
                                                END DO
                                        END DO
                                END DO
                                !n_interact = 0
                                DO i = 1, iend
                                        IF (cell_l_inrange(i,2,xi,yi,zi,ibox)) THEN
                                                n_interact = n_interact + 1
                                                which_interact(n_interact) = i
                                        END IF
                                END DO
                                !n_interact = IAND(n_interact+7,pad8mask) ! padding to multiple of 8 or zero
                                cbmc_cell_n_interact(xi,yi,zi,ibox) = IAND(n_interact+7,NOT(7))
                        END IF
                        !DIR$ VECTOR ALIGNED
                        cbmc_cell_rsp(1:n_interact,1:4,xi,yi,zi,ibox) = &
                                cbmc_cell_rsp_priv(which_interact(1:n_interact),:)
                        IF (need_atom_ti) THEN
                                !DIR$ VECTOR ALIGNED
                                cbmc_cell_ti(1:n_interact,xi,yi,zi,ibox) = &
                                        ti_priv(which_interact(1:n_interact))
                        END IF
                        IF (need_atomtypes) THEN
                                !DIR$ VECTOR ALIGNED
                                cbmc_cell_atomtypes(1:n_interact,xi,yi,zi,ibox) = &
                                        atomtype_priv(which_interact(1:n_interact))
                        END IF
                  END DO
                  END DO
                  END DO
                  !$OMP END DO
          END DO
          !IF (cbmc_cell_list_flag .AND. precalc_atompair_nrg) THEN
          !        !$OMP WORKSHARE
          !        cbmc_cell_ti = cbmc_cell_ti - 1
          !        !$OMP END WORKSHARE
          !END IF
          !$OMP END PARALLEL

  END SUBROUTINE Sector_Setup

  LOGICAL FUNCTION check_overlap_ams(ia, im, is)
          !
          INTEGER, INTENT(IN) :: ia, im, is
          INTEGER, DIMENSION(3) :: cp, sf
          REAL(DP) :: dprp(3)
          REAL(SP) :: irp(3), dxp, dyp, dzp, rijsq, rminsq, sprp
          INTEGER :: this_box
          INTEGER :: i
          INTEGER :: vlen
          INTEGER :: xi, yi, zi
          INTEGER :: ti_solute, icp, sb, lc
          LOGICAL :: overlap, this_overlap

          this_box = widom_molecule%which_box
          dprp(1) = widom_atoms(ia)%rp(1)
          dprp(2) = widom_atoms(ia)%rp(2)
          dprp(3) = widom_atoms(ia)%rp(3)
          IF (box_list(this_box)%int_box_shape <= int_ortho) THEN
                  DO i = 1,3
                        sprp = REAL(dprp(i),SP)
                        irp(i) = sprp
                        icp = NINT(sprp*box_list(this_box)%cell_length_recip(i))
                        cp(i) = icp
                  END DO
                  IF (ANY(ABS(cp)>box_list(this_box)%sectorbound)) THEN
                          sf = cp/(box_list(this_box)%sectorbound+1)
                          cp = cp - sf*box_list(this_box)%length_cells
                          irp = irp - sf*box_list(this_box)%sp_diag_length
                  END IF
          ELSE
                  cp = IDNINT(MATMUL(box_list(this_box)%cell_length_inv,dprp))
                  IF (ANY(ABS(cp)>box_list(this_box)%sectorbound)) THEN
                          sf = cp/(box_list(this_box)%sectorbound+1)
                          cp = cp - sf*box_list(this_box)%length_cells
                          CALL Minimum_Image_Separation(this_box, &
                                  widom_atoms(ia)%rp(1), &
                                  widom_atoms(ia)%rp(2), &
                                  widom_atoms(ia)%rp(3), &
                                  dprp(1), dprp(2), dprp(3))
                  END IF
                  irp = REAL(dprp,SP)
          END IF
          check_overlap_ams = check_overlap_coordinates(irp,cp,ia,is,this_box)
  END FUNCTION check_overlap_ams

  LOGICAL FUNCTION check_overlap_coordinates(irp,cp,ia,is,this_box)
          !
          INTEGER, INTENT(IN) :: ia, is, this_box, cp(3)
          REAL(SP), DIMENSION(3), INTENT(IN) :: irp
          !REAL(SP), DIMENSION(max_adj_cell_atoms) :: rminsq_list
          REAL(SP) :: dxp, dyp, dzp, rijsq, rminsq, sprp
          INTEGER :: i
          INTEGER :: vlen
          INTEGER :: xi, yi, zi
          INTEGER :: ti_solute, icp, sb, lc
          LOGICAL :: overlap, this_overlap

          xi = cp(1)
          yi = cp(2)
          zi = cp(3)
          vlen = n_adj_cell_atoms(xi,yi,zi,this_box)

          IF (read_atompair_rminsq) THEN
                  ti_solute = species_list(is)%wsolute_base+ia
                  !rminsq_list(1:vlen) = sp_atompair_rminsq_table(adj_cell_ti(1:vlen,xi,yi,zi,this_box), &
                  !        species_list(is)%wsolute_base+ia,this_box)
          ELSE IF (calc_rmin_flag) THEN
                  ti_solute = nonbond_list(ia,is)%atom_type_number
                  !rminsq_list(1:vlen) = sp_rminsq_table(adj_cell_ti(1:vlen,xi,yi,zi,this_box), &
                  !        nonbond_list(ia,is)%atom_type_number)
          ELSE
                  rminsq = sp_rcut_lowsq
          END IF
          overlap = .FALSE.
          !DIR$ ASSUME (MOD(vlen,8) .EQ. 0)
          !DIR$ LOOP COUNT = 8, 16, 24, 32, 40
          !DIR$ VECTOR ALIGNED
          DO i = 1, vlen
                dxp = cbmc_cell_rsp(i,1,xi,yi,zi,this_box) - irp(1)
                dyp = cbmc_cell_rsp(i,2,xi,yi,zi,this_box) - irp(2)
                dzp = cbmc_cell_rsp(i,3,xi,yi,zi,this_box) - irp(3)
                rijsq = dxp*dxp
                rijsq = rijsq + dyp*dyp
                rijsq = rijsq + dzp*dzp
                IF (read_atompair_rminsq) THEN
                        rminsq = sp_atompair_rminsq_table(cbmc_cell_ti(i,xi,yi,zi,this_box),ti_solute,this_box)
                ELSE IF (calc_rmin_flag) THEN
                        rminsq = sp_rminsq_table(cbmc_cell_atomtypes(i,xi,yi,zi,this_box),ti_solute)
                END IF
                this_overlap = rijsq < rminsq
                overlap = overlap .OR. this_overlap
          END DO
          check_overlap_coordinates = overlap
          
  END FUNCTION check_overlap_coordinates

  SUBROUTINE CBMC_Cell_List_Setup
          ! This subroutine shouldn't be used anymore
          INTEGER, DIMENSION(3) :: sectormaxbound_old !, map_bound
          INTEGER, DIMENSION(3,nbr_boxes) :: sectorbound_old
          INTEGER :: i_sector, ci(3), dx, dy, dz, xshift, yshift, zshift, nsec_old, nsec, secind
          INTEGER :: sector_ID
          INTEGER :: i, ibox, is, imol, im, ia
          INTEGER, DIMENSION(3) :: sector_atom_ID
          TYPE(Atom_Class), POINTER :: atom_ptr
          REAL(DP) :: xp, yp, zp, cp(3)
!          INTEGER, DIMENSION(3,nbr_boxes) :: cbmc_truth_cube_bound, cut_truth_cube_bound
          INTEGER :: xi, yi, zi, cim(3), xyzi(3), i_dim
          INTEGER :: max_occ_sectors_old, total_atoms
          INTEGER, DIMENSION(2,3) :: tgt_slice, src_slice
          INTEGER, DIMENSION(:), ALLOCATABLE :: xi_pm, yi_pm, zi_pm
          INTEGER :: dummy
          LOGICAL :: asflag
          REAL(DP), DIMENSION(:,:), POINTER :: cell_length_inv_ptr

          
          sectorbound_old = sectorbound_cbmc
          sectormaxbound_old = sectormaxbound_cbmc
          nsec_old = MAXVAL(PRODUCT(length_cells_cbmc,1))
          DO ibox = 1, nbr_boxes
                DO i = 1, 3
                          length_cells_cbmc(i,ibox) = INT(box_list(ibox)%face_distance(i)/rcut_cbmc(ibox))
                          IF (MOD(length_cells_cbmc(i,ibox),2) .EQ. 0) length_cells_cbmc(i,ibox) = length_cells_cbmc(i,ibox) - 1
                          cell_length_inv_cbmc(i,:,ibox) = REAL(length_cells_cbmc(i,ibox),DP) * box_list(ibox)%length_inv(i,:)
                          sectorbound_cbmc(i,ibox) = length_cells_cbmc(i,ibox)/2
                END DO
          END DO
          !cell_length_cbmc = 1.0_DP / cell_length_inv_cbmc

          sectormaxbound_cbmc = MAXVAL(sectorbound_cbmc, 2)+1
          !map_bound = sectormaxbound_cbmc*3+1
          IF (.NOT. ALL( sectormaxbound_cbmc <= sectormaxbound_old)) THEN
                  IF (ALLOCATED(sector_index_map_cbmc)) DEALLOCATE(sector_index_map_cbmc)
                  ALLOCATE(sector_index_map_cbmc(-sectormaxbound_cbmc(1):sectormaxbound_cbmc(1), &
                          -sectormaxbound_cbmc(2):sectormaxbound_cbmc(2), &
                          -sectormaxbound_cbmc(3):sectormaxbound_cbmc(3), &
                          nbr_boxes))
          END IF
          max_occ_sectors_old = max_occ_sectors_cbmc
          total_atoms = DOT_PRODUCT(SUM(nmols(:,1:),2), natoms)
          max_occ_sectors_cbmc=MAX(max_occ_sectors_old,MIN(total_atoms,SUM(PRODUCT(length_cells_cbmc,1))))
          IF (max_occ_sectors_cbmc > max_occ_sectors_old) THEN
                  IF (ALLOCATED(sector_n_atoms_cbmc)) DEALLOCATE(sector_n_atoms_cbmc)
                  IF (ALLOCATED(sector_atoms_cbmc)) DEALLOCATE(sector_atoms_cbmc)
                  ALLOCATE(sector_n_atoms_cbmc(0:max_occ_sectors_cbmc))
                  ALLOCATE(sector_atoms_cbmc(max_sector_natoms_cbmc,max_occ_sectors_cbmc,3))
          END IF
          n_occ_sectors_cbmc = 0
          !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
          sector_n_atoms_cbmc = 0
          sector_index_map_cbmc = 0
          !$OMP END PARALLEL WORKSHARE
          asflag = .TRUE.
          AllocationSizeLoop: DO
          ! place atoms in sectors
          BoxLoop:DO ibox = 1, nbr_boxes
                cell_length_inv_ptr => cell_length_inv_cbmc(:,:,ibox)
                SpeciesLoop:DO is = 1, nspecies
                        sector_atom_ID(3) = is
                        !$OMP PARALLEL DO DEFAULT(SHARED) &
                        !$OMP PRIVATE(i,im,ia,atom_ptr,cp) SCHEDULE(DYNAMIC)
                        MoleculeLoop:DO imol = 1, nmols(is,ibox)
                                im = locate(imol, is, ibox)
                                IF (.NOT. molecule_list(im,is)%live) CYCLE MoleculeLoop
                                AtomLoop:DO ia = 1, natoms(is)
                                        atom_ptr => atom_list(ia,im,is)
                                        IF (.NOT. atom_ptr%exist) CYCLE AtomLoop
                                        cp(1) = atom_ptr%rp(1)
                                        cp(2) = atom_ptr%rp(2)
                                        cp(3) = atom_ptr%rp(3)
                                        atom_ptr%ci_cbmc = IDNINT(MATMUL(cell_length_inv_ptr,cp))
                                        DO i = 1,3
                                                IF (atom_ptr%ci_cbmc(i) > sectorbound_cbmc(i,ibox)) THEN
                                                        atom_ptr%ci_cbmc(i) = atom_ptr%ci_cbmc(i) - length_cells_cbmc(i,ibox)
                                                ELSE IF (atom_ptr%ci_cbmc(i) < -sectorbound_cbmc(i,ibox)) THEN
                                                        atom_ptr%ci_cbmc(i) = atom_ptr%ci_cbmc(i) + length_cells_cbmc(i,ibox)
                                                END IF
                                        END DO
                                END DO AtomLoop
                        END DO MoleculeLoop
                        !$OMP END PARALLEL DO
                        MoleculeLoop2:DO imol = 1, nmols(is,ibox)
                                im = locate(imol, is, ibox)
                                sector_atom_ID(2) = im
                                IF (.NOT. molecule_list(im,is)%live) CYCLE MoleculeLoop2
                                AtomLoop2:DO ia = 1, natoms(is)
                                        IF (.NOT. atom_list(ia,im,is)%exist) CYCLE AtomLoop2
                                        sector_atom_ID(1) = ia
                                        ci = atom_list(ia,im,is)%ci_cbmc
                                        secind = sector_index_map_cbmc(ci(1),ci(2),ci(3),ibox)
                                        IF (secind>0) THEN
                                                sector_n_atoms_cbmc(secind) = sector_n_atoms_cbmc(secind)+1
                                        ELSE
                                                n_occ_sectors_cbmc = n_occ_sectors_cbmc+1
                                                secind = n_occ_sectors_cbmc
                                                sector_n_atoms_cbmc(secind) = 1
                                                sector_index_map_cbmc(ci(1),ci(2),ci(3),ibox) = secind
                                        END IF
                                        IF (sector_n_atoms_cbmc(secind) > max_sector_natoms_cbmc) asflag = .FALSE.
                                        IF (asflag) sector_atoms_cbmc(sector_n_atoms_cbmc(secind),secind,:) = sector_atom_ID
                                END DO AtomLoop2
                        END DO MoleculeLoop2
                END DO SpeciesLoop
                !$OMP PARALLEL PRIVATE(xyzi, i_dim, tgt_slice, src_slice)
                !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
                DO xi = -1, 1
                        DO yi = -1, 1
                                DO zi = -1, 1
                                        IF (xi == 0 .AND. yi == 0 .AND. zi == 0) CYCLE
                                        xyzi = (/ xi, yi, zi /)
                                        DO i_dim = 1,3
                                                IF (xyzi(i_dim) == 0) THEN
                                                        tgt_slice(:,i_dim) = (/ -1, 1 /) * sectorbound_cbmc(i_dim,ibox)
                                                        src_slice(:,i_dim) = tgt_slice(:,i_dim)
                                                ELSE
                                                        tgt_slice(:,i_dim) = sectorbound_cbmc(i_dim,ibox)*xyzi(i_dim) + xyzi(i_dim)
                                                        src_slice(:,i_dim) = -sectorbound_cbmc(i_dim,ibox)*xyzi(i_dim)
                                                END IF
                                        END DO
                                        sector_index_map_cbmc(&
                                                tgt_slice(1,1):tgt_slice(2,1), &
                                                tgt_slice(1,2):tgt_slice(2,2), &
                                                tgt_slice(1,3):tgt_slice(2,3), &
                                                ibox) = &
                                                sector_index_map_cbmc(&
                                                src_slice(1,1):src_slice(2,1), &
                                                src_slice(1,2):src_slice(2,2), &
                                                src_slice(1,3):src_slice(2,3), &
                                                ibox)
                                END DO
                        END DO
                END DO
                !$OMP END DO
                !$OMP END PARALLEL
          END DO BoxLoop
          IF (asflag) THEN
                EXIT AllocationSizeLoop
          ELSE
                max_sector_natoms_cbmc = MAXVAL(sector_n_atoms_cbmc)
                DEALLOCATE(sector_atoms_cbmc)
                ALLOCATE(sector_atoms_cbmc(max_sector_natoms_cbmc,max_occ_sectors_cbmc,3))
                asflag = .TRUE.
                sector_n_atoms_cbmc = 0
          END IF
          END DO AllocationSizeLoop
  END SUBROUTINE CBMC_Cell_List_Setup


  SUBROUTINE Full_Cell_List_Setup
          ! use of this subroutine should be removed
          INTEGER, DIMENSION(3) :: sectormaxbound_old !, map_bound
          INTEGER, DIMENSION(3,nbr_boxes) :: sectorbound_old
          INTEGER :: i_sector, ci(3), dx, dy, dz, xshift, yshift, zshift, nsec_old, nsec, secind
          INTEGER :: sector_ID
          INTEGER :: i, ibox, is, imol, im, ia
          INTEGER, DIMENSION(3) :: sector_atom_ID
          TYPE(Atom_Class), POINTER :: atom_ptr
          REAL(DP) :: xp, yp, zp, cp(3)
          INTEGER :: xi, yi, zi, cim(3)
          INTEGER :: max_occ_sectors_old, total_atoms
          INTEGER, DIMENSION(2) :: xslice, yslice, zslice
          INTEGER, DIMENSION(:,:,:), POINTER :: box_sector_ptr
          INTEGER, DIMENSION(:), ALLOCATABLE :: xi_pm, yi_pm, zi_pm
          INTEGER :: dummy
          REAL(DP) :: max_rcut
          LOGICAL :: asflag
          REAL(DP), DIMENSION(:,:), POINTER :: cell_length_inv_ptr

          sectorbound_old = sectorbound_full
          sectormaxbound_old = sectormaxbound_full
          nsec_old = MAXVAL(PRODUCT(length_cells_full,1))
          DO ibox = 1, nbr_boxes
                max_rcut = MAX(rcut_vdw(ibox), rcut_coul(ibox))
                DO i = 1, 3
                          length_cells_full(i,ibox) = INT(box_list(ibox)%face_distance(i)/max_rcut)
                          IF (MOD(length_cells_full(i,ibox),2) .EQ. 0) length_cells_full(i,ibox) = length_cells_full(i,ibox) - 1
                          cell_length_inv_full(i,:,ibox) = REAL(length_cells_full(i,ibox),DP) * box_list(ibox)%length_inv(i,:)
                          sectorbound_full(i,ibox) = length_cells_full(i,ibox)/2
                END DO
          END DO
          ! cell_length_full = 1.0_DP / cell_length_inv_full

          sectormaxbound_full = 3*MAXVAL(sectorbound_full, 2)+1
          !map_bound = sectormaxbound_full*3+1
          IF (.NOT. ALL( sectormaxbound_full <= sectormaxbound_old)) THEN
                  IF (ALLOCATED(sector_index_map_full)) DEALLOCATE(sector_index_map_full)
                  ALLOCATE(sector_index_map_full(-sectormaxbound_full(1):sectormaxbound_full(1), &
                          -sectormaxbound_full(2):sectormaxbound_full(2), &
                          -sectormaxbound_full(3):sectormaxbound_full(3), &
                          nbr_boxes))
          END IF
          max_occ_sectors_old = max_occ_sectors_full
          total_atoms = DOT_PRODUCT(SUM(nmols(:,1:),2), natoms)
          max_occ_sectors_full=MAX(max_occ_sectors_old,MIN(total_atoms,SUM(PRODUCT(length_cells_full,1))))
          IF (max_occ_sectors_full > max_occ_sectors_old) THEN
                  IF (ALLOCATED(sector_n_atoms_full)) DEALLOCATE(sector_n_atoms_full)
                  IF (ALLOCATED(sector_atoms_full)) DEALLOCATE(sector_atoms_full)
                  ALLOCATE(sector_n_atoms_full(0:max_occ_sectors_full))
                  ALLOCATE(sector_atoms_full(max_sector_natoms_full,max_occ_sectors_full,3))
          END IF
          n_occ_sectors_full = 0
          !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
          sector_n_atoms_full = 0
          sector_index_map_full = 0
          !$OMP END PARALLEL WORKSHARE
          asflag = .TRUE.
          AllocationSizeLoop: DO
          ! place atoms in sectors
          BoxLoop:DO ibox = 1, nbr_boxes
                cell_length_inv_ptr => cell_length_inv_full(:,:,ibox)
                SpeciesLoop:DO is = 1, nspecies
                        sector_atom_ID(3) = is
                        !$OMP PARALLEL DO DEFAULT(SHARED) &
                        !$OMP PRIVATE(i,im,ia,atom_ptr,cp) SCHEDULE(DYNAMIC)
                        MoleculeLoop:DO imol = 1, nmols(is,ibox)
                                im = locate(imol, is, ibox)
                                IF (.NOT. molecule_list(im,is)%live) CYCLE MoleculeLoop
                                AtomLoop:DO ia = 1, natoms(is)
                                        atom_ptr => atom_list(ia,im,is)
                                        IF (.NOT. atom_ptr%exist) CYCLE AtomLoop
                                        cp(1) = atom_ptr%rp(1)
                                        cp(2) = atom_ptr%rp(2)
                                        cp(3) = atom_ptr%rp(3)
                                        atom_ptr%ci_full = IDNINT(MATMUL(cell_length_inv_ptr,cp))
                                        DO i = 1,3
                                                IF (atom_ptr%ci_full(i) > sectorbound_full(i,ibox)) THEN
                                                        atom_ptr%ci_full(i) = atom_ptr%ci_full(i) - length_cells_full(i,ibox)
                                                ELSE IF (atom_ptr%ci_full(i) < -sectorbound_full(i,ibox)) THEN
                                                        atom_ptr%ci_full(i) = atom_ptr%ci_full(i) + length_cells_full(i,ibox)
                                                END IF
                                        END DO
                                END DO AtomLoop
                        END DO MoleculeLoop
                        !$OMP END PARALLEL DO
                        MoleculeLoop2:DO imol = 1, nmols(is,ibox)
                                im = locate(imol, is, ibox)
                                sector_atom_ID(2) = im
                                IF (.NOT. molecule_list(im,is)%live) CYCLE MoleculeLoop2
                                AtomLoop2:DO ia = 1, natoms(is)
                                        IF (.NOT. atom_list(ia,im,is)%exist) CYCLE AtomLoop2
                                        sector_atom_ID(1) = ia
                                        ci = atom_list(ia,im,is)%ci_full
                                        secind = sector_index_map_full(ci(1),ci(2),ci(3),ibox)
                                        IF (secind>0) THEN
                                                sector_n_atoms_full(secind) = sector_n_atoms_full(secind)+1
                                        ELSE
                                                n_occ_sectors_full = n_occ_sectors_full+1
                                                secind = n_occ_sectors_full
                                                sector_n_atoms_full(secind) = 1
                                                sector_index_map_full(ci(1),ci(2),ci(3),ibox) = secind
                                        END IF
                                        IF (sector_n_atoms_full(secind) > max_sector_natoms_full) asflag = .FALSE.
                                        IF (asflag) sector_atoms_full(sector_n_atoms_full(secind),secind,:) = sector_atom_ID
                                END DO AtomLoop2
                        END DO MoleculeLoop2
                END DO SpeciesLoop
                box_sector_ptr => sector_index_map_full( &
                        -sectorbound_full(1,ibox):sectorbound_full(1,ibox), &
                        -sectorbound_full(2,ibox):sectorbound_full(2,ibox), &
                        -sectorbound_full(3,ibox):sectorbound_full(3,ibox), &
                        ibox)
                DO xi = -1, 1
                        xslice(1) = -sectorbound_full(1,ibox)+xi*length_cells_full(1,ibox)
                        xslice(2) = sectorbound_full(1,ibox)+xi*length_cells_full(1,ibox)
                        DO yi = -1, 1
                                yslice(1) = -sectorbound_full(2,ibox)+yi*length_cells_full(2,ibox)
                                yslice(2) = sectorbound_full(2,ibox)+yi*length_cells_full(2,ibox)
                                DO zi = -1, 1
                                        IF (xi == 0 .AND. yi == 0 .AND. zi == 0) CYCLE
                                        zslice(1) = -sectorbound_full(3,ibox)+zi*length_cells_full(3,ibox)
                                        zslice(2) = sectorbound_full(3,ibox)+zi*length_cells_full(3,ibox)
                                        sector_index_map_full(&
                                                xslice(1):xslice(2), &
                                                yslice(1):yslice(2), &
                                                zslice(1):zslice(2), &
                                                ibox) = box_sector_ptr
                                END DO
                        END DO
                END DO
          END DO BoxLoop
          IF (asflag) THEN
                EXIT AllocationSizeLoop
          ELSE
                max_sector_natoms_full = MAXVAL(sector_n_atoms_full)
                DEALLOCATE(sector_atoms_full)
                ALLOCATE(sector_atoms_full(max_sector_natoms_full,max_occ_sectors_full,3))
                asflag = .TRUE.
                sector_n_atoms_full = 0
          END IF
          END DO AllocationSizeLoop
  END SUBROUTINE Full_Cell_List_Setup



  SUBROUTINE CBMC_Cell_List_Setup
          INTEGER, DIMENSION(3) :: sectormaxbound_old !, map_bound
          INTEGER, DIMENSION(3,nbr_boxes) :: sectorbound_old
          INTEGER :: i_sector, ci(3), dx, dy, dz, xshift, yshift, zshift, nsec_old, nsec, secind
          INTEGER :: sector_ID
          INTEGER :: i, ibox, is, imol, im, ia
          INTEGER, DIMENSION(3) :: sector_atom_ID
          TYPE(Atom_Class), POINTER :: atom_ptr
          REAL(DP) :: xp, yp, zp, cp(3)
!          INTEGER, DIMENSION(3,nbr_boxes) :: cbmc_truth_cube_bound, cut_truth_cube_bound
          INTEGER :: xi, yi, zi, cim(3)
          INTEGER :: max_occ_sectors_old, total_atoms
          INTEGER, DIMENSION(2) :: xslice, yslice, zslice
          INTEGER, DIMENSION(:,:,:), POINTER :: box_sector_ptr
          INTEGER, DIMENSION(:), ALLOCATABLE :: xi_pm, yi_pm, zi_pm
          INTEGER :: dummy
          LOGICAL :: asflag

          
          sectorbound_old = sectorbound_cbmc
          sectormaxbound_old = sectormaxbound_cbmc
          nsec_old = MAXVAL(PRODUCT(length_cells_cbmc,1))
          DO ibox = 1, nbr_boxes
                DO i = 1, 3
                          length_cells_cbmc(i,ibox) = INT(box_list(ibox)%length(i,i)/rcut_cbmc(ibox))
                          IF (MOD(length_cells_cbmc(i,ibox),2) .EQ. 0) length_cells_cbmc(i,ibox) = length_cells_cbmc(i,ibox) - 1
                          cell_length_inv_cbmc(i,ibox) = REAL(length_cells_cbmc(i,ibox),DP) / box_list(ibox)%length(i,i)
                          sectorbound_cbmc(i,ibox) = length_cells_cbmc(i,ibox)/2
                END DO
          END DO
          !cell_length_cbmc = 1.0_DP / cell_length_inv_cbmc

          sectormaxbound_cbmc = 3*MAXVAL(sectorbound_cbmc, 2)+1
          !map_bound = sectormaxbound_cbmc*3+1
          IF (.NOT. ALL( sectormaxbound_cbmc <= sectormaxbound_old)) THEN
                  IF (ALLOCATED(sector_index_map_cbmc)) DEALLOCATE(sector_index_map_cbmc)
                  ALLOCATE(sector_index_map_cbmc(-sectormaxbound_cbmc(1):sectormaxbound_cbmc(1), &
                          -sectormaxbound_cbmc(2):sectormaxbound_cbmc(2), &
                          -sectormaxbound_cbmc(3):sectormaxbound_cbmc(3), &
                          nbr_boxes))
          END IF
          max_occ_sectors_old = max_occ_sectors_cbmc
          total_atoms = DOT_PRODUCT(SUM(nmols(:,1:),2), natoms)
          max_occ_sectors_cbmc=MAX(max_occ_sectors_old,MIN(total_atoms,SUM(PRODUCT(length_cells_cbmc,1))))
          IF (max_occ_sectors_cbmc > max_occ_sectors_old) THEN
                  IF (ALLOCATED(sector_n_atoms_cbmc)) DEALLOCATE(sector_n_atoms_cbmc)
                  IF (ALLOCATED(sector_atoms_cbmc)) DEALLOCATE(sector_atoms_cbmc)
                  ALLOCATE(sector_n_atoms_cbmc(0:max_occ_sectors_cbmc))
                  ALLOCATE(sector_atoms_cbmc(max_sector_natoms_cbmc,max_occ_sectors_cbmc,3))
          END IF
          n_occ_sectors_cbmc = 0
          !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
          sector_n_atoms_cbmc = 0
          sector_index_map_cbmc = 0
          !$OMP END PARALLEL WORKSHARE
          asflag = .TRUE.
          AllocationSizeLoop: DO
          ! place atoms in sectors
          BoxLoop:DO ibox = 1, nbr_boxes
                SpeciesLoop:DO is = 1, nspecies
                        sector_atom_ID(3) = is
                        !$OMP PARALLEL DO DEFAULT(SHARED) &
                        !$OMP PRIVATE(i,im,ia,atom_ptr,cp) SCHEDULE(DYNAMIC)
                        MoleculeLoop:DO imol = 1, nmols(is,ibox)
                                im = locate(imol, is, ibox)
                                IF (.NOT. molecule_list(im,is)%live) CYCLE MoleculeLoop
                                AtomLoop:DO ia = 1, natoms(is)
                                        atom_ptr => atom_list(ia,im,is)
                                        IF (.NOT. atom_ptr%exist) CYCLE AtomLoop
                                        cp(1) = atom_ptr%rxp
                                        cp(2) = atom_ptr%ryp
                                        cp(3) = atom_ptr%rzp
                                        atom_ptr%ci_cbmc = IDNINT(cp*cell_length_inv_cbmc(:,ibox))
                                        DO i = 1,3
                                                IF (atom_ptr%ci_cbmc(i) > sectorbound_cbmc(i,ibox)) THEN
                                                        atom_ptr%ci_cbmc(i) = atom_ptr%ci_cbmc(i) - length_cells_cbmc(i,ibox)
                                                ELSE IF (atom_ptr%ci_cbmc(i) < -sectorbound_cbmc(i,ibox)) THEN
                                                        atom_ptr%ci_cbmc(i) = atom_ptr%ci_cbmc(i) + length_cells_cbmc(i,ibox)
                                                END IF
                                        END DO
                                END DO AtomLoop
                        END DO MoleculeLoop
                        !$OMP END PARALLEL DO
                        MoleculeLoop2:DO imol = 1, nmols(is,ibox)
                                im = locate(imol, is, ibox)
                                sector_atom_ID(2) = im
                                IF (.NOT. molecule_list(im,is)%live) CYCLE MoleculeLoop2
                                AtomLoop2:DO ia = 1, natoms(is)
                                        IF (.NOT. atom_list(ia,im,is)%exist) CYCLE AtomLoop2
                                        sector_atom_ID(1) = ia
                                        ci = atom_list(ia,im,is)%ci_cbmc
                                        secind = sector_index_map_cbmc(ci(1),ci(2),ci(3),ibox)
                                        IF (secind>0) THEN
                                                sector_n_atoms_cbmc(secind) = sector_n_atoms_cbmc(secind)+1
                                        ELSE
                                                n_occ_sectors_cbmc = n_occ_sectors_cbmc+1
                                                secind = n_occ_sectors_cbmc
                                                sector_n_atoms_cbmc(secind) = 1
                                                sector_index_map_cbmc(ci(1),ci(2),ci(3),ibox) = secind
                                        END IF
                                        IF (sector_n_atoms_cbmc(secind) > max_sector_natoms_cbmc) asflag = .FALSE.
                                        IF (asflag) sector_atoms_cbmc(sector_n_atoms_cbmc(secind),secind,:) = sector_atom_ID
                                END DO AtomLoop2
                        END DO MoleculeLoop2
                END DO SpeciesLoop
                box_sector_ptr => sector_index_map_cbmc( &
                        -sectorbound_cbmc(1,ibox):sectorbound_cbmc(1,ibox), &
                        -sectorbound_cbmc(2,ibox):sectorbound_cbmc(2,ibox), &
                        -sectorbound_cbmc(3,ibox):sectorbound_cbmc(3,ibox), &
                        ibox)
                DO xi = -1, 1
                        xslice(1) = -sectorbound_cbmc(1,ibox)+xi*length_cells_cbmc(1,ibox)
                        xslice(2) = sectorbound_cbmc(1,ibox)+xi*length_cells_cbmc(1,ibox)
                        DO yi = -1, 1
                                yslice(1) = -sectorbound_cbmc(2,ibox)+yi*length_cells_cbmc(2,ibox)
                                yslice(2) = sectorbound_cbmc(2,ibox)+yi*length_cells_cbmc(2,ibox)
                                DO zi = -1, 1
                                        IF (xi == 0 .AND. yi == 0 .AND. zi == 0) CYCLE
                                        zslice(1) = -sectorbound_cbmc(3,ibox)+zi*length_cells_cbmc(3,ibox)
                                        zslice(2) = sectorbound_cbmc(3,ibox)+zi*length_cells_cbmc(3,ibox)
                                        sector_index_map_cbmc(&
                                                xslice(1):xslice(2), &
                                                yslice(1):yslice(2), &
                                                zslice(1):zslice(2), &
                                                ibox) = box_sector_ptr
                                END DO
                        END DO
                END DO
          END DO BoxLoop
          IF (asflag) THEN
                EXIT AllocationSizeLoop
          ELSE
                max_sector_natoms_cbmc = MAXVAL(sector_n_atoms_cbmc)
                DEALLOCATE(sector_atoms_cbmc)
                ALLOCATE(sector_atoms_cbmc(max_sector_natoms_cbmc,max_occ_sectors_cbmc,3))
                asflag = .TRUE.
                sector_n_atoms_cbmc = 0
          END IF
          END DO AllocationSizeLoop
  END SUBROUTINE CBMC_Cell_List_Setup


  SUBROUTINE Full_Cell_List_Setup
          INTEGER, DIMENSION(3) :: sectormaxbound_old !, map_bound
          INTEGER, DIMENSION(3,nbr_boxes) :: sectorbound_old
          INTEGER :: i_sector, ci(3), dx, dy, dz, xshift, yshift, zshift, nsec_old, nsec, secind
          INTEGER :: sector_ID
          INTEGER :: i, ibox, is, imol, im, ia
          INTEGER, DIMENSION(3) :: sector_atom_ID
          TYPE(Atom_Class), POINTER :: atom_ptr
          REAL(DP) :: xp, yp, zp, cp(3)
          INTEGER :: xi, yi, zi, cim(3)
          INTEGER :: max_occ_sectors_old, total_atoms
          INTEGER, DIMENSION(2) :: xslice, yslice, zslice
          INTEGER, DIMENSION(:,:,:), POINTER :: box_sector_ptr
          INTEGER, DIMENSION(:), ALLOCATABLE :: xi_pm, yi_pm, zi_pm
          INTEGER :: dummy
          REAL(DP) :: max_rcut
          LOGICAL :: asflag

          sectorbound_old = sectorbound_full
          sectormaxbound_old = sectormaxbound_full
          nsec_old = MAXVAL(PRODUCT(length_cells_full,1))
          DO ibox = 1, nbr_boxes
                max_rcut = MAX(rcut_vdw(ibox), rcut_coul(ibox))
                DO i = 1, 3
                          length_cells_full(i,ibox) = INT(box_list(ibox)%length(i,i)/max_rcut)
                          IF (MOD(length_cells_full(i,ibox),2) .EQ. 0) length_cells_full(i,ibox) = length_cells_full(i,ibox) - 1
                          cell_length_inv_full(i,ibox) = REAL(length_cells_full(i,ibox),DP) / box_list(ibox)%length(i,i)
                          sectorbound_full(i,ibox) = length_cells_full(i,ibox)/2
                END DO
          END DO
          ! cell_length_full = 1.0_DP / cell_length_inv_full

          sectormaxbound_full = 3*MAXVAL(sectorbound_full, 2)+1
          !map_bound = sectormaxbound_full*3+1
          IF (.NOT. ALL( sectormaxbound_full <= sectormaxbound_old)) THEN
                  IF (ALLOCATED(sector_index_map_full)) DEALLOCATE(sector_index_map_full)
                  ALLOCATE(sector_index_map_full(-sectormaxbound_full(1):sectormaxbound_full(1), &
                          -sectormaxbound_full(2):sectormaxbound_full(2), &
                          -sectormaxbound_full(3):sectormaxbound_full(3), &
                          nbr_boxes))
          END IF
          max_occ_sectors_old = max_occ_sectors_full
          total_atoms = DOT_PRODUCT(SUM(nmols(:,1:),2), natoms)
          max_occ_sectors_full=MAX(max_occ_sectors_old,MIN(total_atoms,SUM(PRODUCT(length_cells_full,1))))
          IF (max_occ_sectors_full > max_occ_sectors_old) THEN
                  IF (ALLOCATED(sector_n_atoms_full)) DEALLOCATE(sector_n_atoms_full)
                  IF (ALLOCATED(sector_atoms_full)) DEALLOCATE(sector_atoms_full)
                  ALLOCATE(sector_n_atoms_full(0:max_occ_sectors_full))
                  ALLOCATE(sector_atoms_full(max_sector_natoms_full,max_occ_sectors_full,3))
          END IF
          n_occ_sectors_full = 0
          !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
          sector_n_atoms_full = 0
          sector_index_map_full = 0
          !$OMP END PARALLEL WORKSHARE
          asflag = .TRUE.
          AllocationSizeLoop: DO
          ! place atoms in sectors
          BoxLoop:DO ibox = 1, nbr_boxes
                SpeciesLoop:DO is = 1, nspecies
                        sector_atom_ID(3) = is
                        !$OMP PARALLEL DO DEFAULT(SHARED) &
                        !$OMP PRIVATE(i,im,ia,atom_ptr,cp) SCHEDULE(DYNAMIC)
                        MoleculeLoop:DO imol = 1, nmols(is,ibox)
                                im = locate(imol, is, ibox)
                                IF (.NOT. molecule_list(im,is)%live) CYCLE MoleculeLoop
                                AtomLoop:DO ia = 1, natoms(is)
                                        atom_ptr => atom_list(ia,im,is)
                                        IF (.NOT. atom_ptr%exist) CYCLE AtomLoop
                                        cp(1) = atom_ptr%rxp
                                        cp(2) = atom_ptr%ryp
                                        cp(3) = atom_ptr%rzp
                                        atom_ptr%ci_full = IDNINT(cp*cell_length_inv_full(:,ibox))
                                        DO i = 1,3
                                                IF (atom_ptr%ci_full(i) > sectorbound_full(i,ibox)) THEN
                                                        atom_ptr%ci_full(i) = atom_ptr%ci_full(i) - length_cells_full(i,ibox)
                                                ELSE IF (atom_ptr%ci_full(i) < -sectorbound_full(i,ibox)) THEN
                                                        atom_ptr%ci_full(i) = atom_ptr%ci_full(i) + length_cells_full(i,ibox)
                                                END IF
                                        END DO
                                END DO AtomLoop
                        END DO MoleculeLoop
                        !$OMP END PARALLEL DO
                        MoleculeLoop2:DO imol = 1, nmols(is,ibox)
                                im = locate(imol, is, ibox)
                                sector_atom_ID(2) = im
                                IF (.NOT. molecule_list(im,is)%live) CYCLE MoleculeLoop2
                                AtomLoop2:DO ia = 1, natoms(is)
                                        IF (.NOT. atom_list(ia,im,is)%exist) CYCLE AtomLoop2
                                        sector_atom_ID(1) = ia
                                        ci = atom_list(ia,im,is)%ci_full
                                        secind = sector_index_map_full(ci(1),ci(2),ci(3),ibox)
                                        IF (secind>0) THEN
                                                sector_n_atoms_full(secind) = sector_n_atoms_full(secind)+1
                                        ELSE
                                                n_occ_sectors_full = n_occ_sectors_full+1
                                                secind = n_occ_sectors_full
                                                sector_n_atoms_full(secind) = 1
                                                sector_index_map_full(ci(1),ci(2),ci(3),ibox) = secind
                                        END IF
                                        IF (sector_n_atoms_full(secind) > max_sector_natoms_full) asflag = .FALSE.
                                        IF (asflag) sector_atoms_full(sector_n_atoms_full(secind),secind,:) = sector_atom_ID
                                END DO AtomLoop2
                        END DO MoleculeLoop2
                END DO SpeciesLoop
                box_sector_ptr => sector_index_map_full( &
                        -sectorbound_full(1,ibox):sectorbound_full(1,ibox), &
                        -sectorbound_full(2,ibox):sectorbound_full(2,ibox), &
                        -sectorbound_full(3,ibox):sectorbound_full(3,ibox), &
                        ibox)
                DO xi = -1, 1
                        xslice(1) = -sectorbound_full(1,ibox)+xi*length_cells_full(1,ibox)
                        xslice(2) = sectorbound_full(1,ibox)+xi*length_cells_full(1,ibox)
                        DO yi = -1, 1
                                yslice(1) = -sectorbound_full(2,ibox)+yi*length_cells_full(2,ibox)
                                yslice(2) = sectorbound_full(2,ibox)+yi*length_cells_full(2,ibox)
                                DO zi = -1, 1
                                        IF (xi == 0 .AND. yi == 0 .AND. zi == 0) CYCLE
                                        zslice(1) = -sectorbound_full(3,ibox)+zi*length_cells_full(3,ibox)
                                        zslice(2) = sectorbound_full(3,ibox)+zi*length_cells_full(3,ibox)
                                        sector_index_map_full(&
                                                xslice(1):xslice(2), &
                                                yslice(1):yslice(2), &
                                                zslice(1):zslice(2), &
                                                ibox) = box_sector_ptr
                                END DO
                        END DO
                END DO
          END DO BoxLoop
          IF (asflag) THEN
                EXIT AllocationSizeLoop
          ELSE
                max_sector_natoms_full = MAXVAL(sector_n_atoms_full)
                DEALLOCATE(sector_atoms_full)
                ALLOCATE(sector_atoms_full(max_sector_natoms_full,max_occ_sectors_full,3))
                asflag = .TRUE.
                sector_n_atoms_full = 0
          END IF
          END DO AllocationSizeLoop
  END SUBROUTINE Full_Cell_List_Setup



END MODULE Sector_Routines
