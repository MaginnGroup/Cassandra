
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
  !$ USE OMP_LIB

  IMPLICIT NONE

  


CONTAINS
  SUBROUTINE Sector_Setup
          INTEGER, DIMENSION(3) :: adj_cellmaxbound_old
          INTEGER :: i_sector, dx, dy, dz, xshift, yshift, zshift, nsec_old, nsec, secind
          INTEGER :: sector_ID
          INTEGER :: i, ibox, is, imol, im, ia
          INTEGER, DIMENSION(3) :: sector_atom_ID
          TYPE(Atom_Class), POINTER :: atom_ptr
          REAL(DP) :: xp, yp, zp
!          INTEGER, DIMENSION(3,nbr_boxes) :: cbmc_truth_cube_bound, cut_truth_cube_bound
          INTEGER :: xi, yi, zi, cim(3), xyzi(3), i_dim
          INTEGER, DIMENSION(2,3) :: tgt_slice, src_slice
          INTEGER, DIMENSION(:), ALLOCATABLE :: xi_pm, yi_pm, zi_pm
          INTEGER :: dummy

          REAL(SP), DIMENSION(maxboxnatoms,3,nbr_boxes) :: sp_live_atom_rsp
          INTEGER, DIMENSION(maxboxnatoms,3,nbr_boxes) :: live_atom_cp
          INTEGER(2), DIMENSION(maxboxnatoms,nbr_boxes) :: live_atom_ti
          INTEGER, DIMENSION(4,maxboxnatoms,nbr_boxes) :: ci_list
          INTEGER(2) :: bsolvent
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
          INTEGER(2), DIMENSION(:,:,:,:,:), ALLOCATABLE :: this_cell_ti
          REAL(SP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: this_cell_rsp

          !sp_live_atom_rsp = REAL(live_atom_rsp)
          DO ibox = 1, nbr_boxes
                istart = 1
                DO is = 1, nspecies
                        inlive = nlive(is,ibox)
                        IF (inlive < 1) CYCLE
                        inatoms = natoms(is)
                        bsolvent = INT(species_list(is)%solvent_base,2)
                        vlen = inlive*inatoms
                        iend = istart + vlen - 1
                        !$OMP PARALLEL WORKSHARE
                        sp_live_atom_rsp(istart:iend,:,ibox) = RESHAPE(REAL(live_atom_rsp(1:inatoms,1:inlive,1:3,is,ibox),SP), &
                                (/ vlen, 3 /))
                        !$OMP END PARALLEL WORKSHARE
                        IF (l_not_all_exist) THEN
                                !$OMP PARALLEL WORKSHARE
                                i_exist(istart:iend) = RESHAPE(live_atom_exist(1:inatoms,1:inlive,is,ibox), (/ vlen /))
                                !$OMP END PARALLEL WORKSHARE
                        END IF
                        IF (read_atompair_rminsq) THEN
                                !$OMP PARALLEL WORKSHARE
                                live_atom_ti(istart:iend,ibox) = RESHAPE(SPREAD(bsolvent+INT(vec123(1:inatoms),2),2,inlive), (/ vlen /))
                                !$OMP END PARALLEL WORKSHARE
                        ELSE IF (calc_rmin_flag) THEN
                                !$OMP PARALLEL WORKSHARE
                                live_atom_ti(istart:iend,ibox) = RESHAPE(SPREAD( &
                                        INT(nonbond_list(1:inatoms,is)%atom_type_number,2),2,inlive), (/ vlen /))
                                !$OMP END PARALLEL WORKSHARE
                        END IF
                        istart = istart + vlen
                END DO
                vlen = istart - 1
                box_vlen(ibox) = vlen
                IF (l_not_all_exist) THEN
                        IF (.NOT. ALL(i_exist(1:vlen))) THEN
                                n_i_exist = COUNT(i_exist(1:vlen))
                                DO i = 1,3
                                        sp_live_atom_rsp(1:n_i_exist,i,ibox) = PACK(sp_live_atom_rsp(1:vlen,i,ibox),i_exist(1:vlen))
                                END DO
                                live_atom_ti(1:n_i_exist,ibox) = PACK(live_atom_ti(1:vlen,ibox),i_exist(1:vlen))
                                box_vlen(ibox) = n_i_exist
                        END IF
                END IF
          END DO

          DO ibox = 1, nbr_boxes
                DO i = 1, 3
                          box_list(ibox)%length_cells(i) = INT(box_list(ibox)%face_distance(i)/max_rmin)
                          IF (MOD(box_list(ibox)%length_cells(i),2) .EQ. 0) box_list(ibox)%length_cells(i) = box_list(ibox)%length_cells(i) - 1
                          box_list(ibox)%cell_length_inv(i,:) = REAL(box_list(ibox)%length_cells(i),DP) * box_list(ibox)%length_inv(i,:)
                          box_list(ibox)%cell_length_recip(i) = REAL(REAL(box_list(ibox)%length_cells(i),DP)/box_list(ibox)%length(i,i),SP) ! kind 4
                          box_list(ibox)%sectorbound(i) = box_list(ibox)%length_cells(i)/2
                          box_list(ibox)%sp_diag_length(i) = REAL(box_list(ibox)%length(i,i),SP)
                END DO
                box_list(ibox)%real_length_cells = REAL(box_list(ibox)%length_cells,SP)
          END DO
          DO ibox = 1, nbr_boxes
                vlen = box_vlen(ibox)
                l_ortho = box_list(ibox)%int_box_shape <= int_ortho
                DO i_dim = 1, 3
                        IF (l_ortho) THEN
                                lbox = box_list(ibox)%sp_diag_length(i_dim)
                                clr = box_list(ibox)%cell_length_recip(i_dim)
                        ELSE
                                lbox = 1.0
                                clr = box_list(ibox)%real_length_cells(i_dim)
                        END IF
                        cp_ub = box_list(ibox)%sectorbound(i_dim)
                        cp_lb = -cp_ub
                        lc = box_list(ibox)%length_cells(i_dim)
                        rp_ub = 0.5*lbox
                        rp_lb = -rp_ub
                        !$OMP PARALLEL
                        !$OMP DO SIMD PRIVATE(rsp,cp) SCHEDULE(STATIC)
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
                        END DO
                        !$OMP END DO SIMD
                        !$OMP END PARALLEL
                END DO
          END DO
          sectormaxbound = box_list(1)%sectorbound
          DO ibox = 2, nbr_boxes
                sectormaxbound = MAX(sectormaxbound,box_list(ibox)%sectorbound)
          END DO
          sectormaxbound = sectormaxbound + 1

          !sectormaxbound = MAXVAL(sectorbound, 2)+1
          adj_cellmaxbound_old = adj_cellmaxbound
          adj_cellmaxbound = sectormaxbound-1
          !map_bound = sectormaxbound*3+1
          ALLOCATE(n_cell_atoms( &
                  -sectormaxbound(1):sectormaxbound(1), &
                  -sectormaxbound(2):sectormaxbound(2), &
                  -sectormaxbound(3):sectormaxbound(3), &
                  nbr_boxes))
          IF (ALLOCATED(n_adj_cell_atoms)) THEN
                  IF (ANY(adj_cellmaxbound>adj_cellmaxbound_old)) THEN
                          DEALLOCATE(n_adj_cell_atoms)
                  END IF
          END IF
          IF (.NOT. ALLOCATED(n_adj_cell_atoms)) ALLOCATE(n_adj_cell_atoms( &
                  -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                  -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                  -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                  nbr_boxes))
          !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
          n_cell_atoms = 0
          !$OMP END PARALLEL WORKSHARE
          DO ibox = 1, nbr_boxes
                DO i = 1, box_vlen(ibox)
                        ci(1:3) = live_atom_cp(i,:,ibox)
                        nca = n_cell_atoms(ci(1),ci(2),ci(3),ibox)+1
                        ci(4) = nca
                        n_cell_atoms(ci(1),ci(2),ci(3),ibox) = nca
                        ci_list(:,i,ibox) = ci
                END DO
          END DO
          !$OMP PARALLEL WORKSHARE
          max_sector_natoms = MAXVAL(n_cell_atoms)
          !$OMP END PARALLEL WORKSHARE
          ALLOCATE(this_cell_rsp(max_sector_natoms, 3, &
                  -sectormaxbound(1):sectormaxbound(1), &
                  -sectormaxbound(2):sectormaxbound(2), &
                  -sectormaxbound(3):sectormaxbound(3), &
                  nbr_boxes))
          ALLOCATE(this_cell_ti(max_sector_natoms, &
                  -sectormaxbound(1):sectormaxbound(1), &
                  -sectormaxbound(2):sectormaxbound(2), &
                  -sectormaxbound(3):sectormaxbound(3), &
                  nbr_boxes))
          DO ibox = 1, nbr_boxes
                !$OMP PARALLEL PRIVATE(ci)
                !$OMP DO SCHEDULE(STATIC)
                DO i = 1, box_vlen(ibox)
                        ci = ci_list(:,i,ibox)
                        this_cell_rsp(ci(4),1:3,ci(1),ci(2),ci(3),ibox) = &
                                sp_live_atom_rsp(i,:,ibox)
                        this_cell_ti(ci(4),ci(1),ci(2),ci(3),ibox) = &
                                live_atom_ti(i,ibox)
                END DO
                !$OMP END DO
                !$OMP END PARALLEL
                l_ortho = box_list(ibox)%int_box_shape <= int_ortho
                IF (l_ortho) THEN
                        DO i_dim = 1, 3
                                boxlen(i_dim) = REAL(box_list(ibox)%length(i_dim,i_dim),SP)
                        END DO
                ELSE
                        boxlen = 1.0
                        h11 = REAL(box_list(ibox)%length(1,1),SP)
                        h21 = REAL(box_list(ibox)%length(2,1),SP)
                        h31 = REAL(box_list(ibox)%length(3,1),SP)
                        h12 = REAL(box_list(ibox)%length(1,2),SP)
                        h22 = REAL(box_list(ibox)%length(2,2),SP)
                        h32 = REAL(box_list(ibox)%length(3,2),SP)
                        h13 = REAL(box_list(ibox)%length(1,3),SP)
                        h23 = REAL(box_list(ibox)%length(2,3),SP)
                        h33 = REAL(box_list(ibox)%length(3,3),SP)
                END IF
                !!!!$OMP PARALLEL PRIVATE(xyzi, i_dim, tgt_slice, src_slice)
                !!!!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
                DO xi = -1, 1
                        DO yi = -1, 1
                                DO zi = -1, 1
                                        IF (xi == 0 .AND. yi == 0 .AND. zi == 0) CYCLE
                                        xyzi = (/ xi, yi, zi /)
                                        unwrap_shifter = xyzi*boxlen
                                        DO i_dim = 1, 3
                                                IF (xyzi(i_dim) == 0) THEN
                                                        tgt_slice(:,i_dim) = (/ -1, 1 /) * box_list(ibox)%sectorbound(i_dim)
                                                        src_slice(:,i_dim) = tgt_slice(:,i_dim)
                                                ELSE
                                                        tgt_slice(:,i_dim) = box_list(ibox)%sectorbound(i_dim)*xyzi(i_dim) + xyzi(i_dim)
                                                        src_slice(:,i_dim) = -box_list(ibox)%sectorbound(i_dim)*xyzi(i_dim)
                                                END IF
                                        END DO
                                        DO i_dim = 1, 3
                                                !$OMP PARALLEL WORKSHARE
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
                                                !$OMP END PARALLEL WORKSHARE
                                        END DO
                                        !$OMP PARALLEL WORKSHARE
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
                                        !$OMP END PARALLEL WORKSHARE
                                END DO
                        END DO
                END DO
                IF (l_ortho) CYCLE
                !$OMP PARALLEL PRIVATE(rxp,ryp,rzp,isp,i)
                !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
                DO zi = -1-box_list(ibox)%sectorbound(3), 1+box_list(ibox)%sectorbound(3)
                        DO yi = -1-box_list(ibox)%sectorbound(2), 1+box_list(ibox)%sectorbound(2)
                                DO xi = -1-box_list(ibox)%sectorbound(1), 1+box_list(ibox)%sectorbound(1)
                                        DO i = 1, n_cell_atoms(xi,yi,zi,ibox)
                                                isp = this_cell_rsp(i,1,xi,yi,zi,ibox)
                                                rxp = h11*isp
                                                ryp = h21*isp
                                                rzp = h31*isp
                                                isp = this_cell_rsp(i,2,xi,yi,zi,ibox)
                                                rxp = rxp + h12*isp
                                                ryp = ryp + h22*isp
                                                rzp = rzp + h32*isp
                                                isp = this_cell_rsp(i,3,xi,yi,zi,ibox)
                                                rxp = rxp + h13*isp
                                                ryp = ryp + h23*isp
                                                rzp = rzp + h33*isp
                                                this_cell_rsp(i,1,xi,yi,zi,ibox) = rxp
                                                this_cell_rsp(i,2,xi,yi,zi,ibox) = ryp
                                                this_cell_rsp(i,3,xi,yi,zi,ibox) = rzp
                                        END DO
                                END DO
                        END DO
                END DO
                !$OMP END DO
                !$OMP END PARALLEL

                !!!!$OMP END DO
                !!!!$OMP END PARALLEL
          END DO
          DO ibox = 1, nbr_boxes
                  !$OMP PARALLEL
                  !$OMP DO COLLAPSE(3)
                  DO zi = -box_list(ibox)%sectorbound(3), box_list(ibox)%sectorbound(3)
                          DO yi = -box_list(ibox)%sectorbound(2), box_list(ibox)%sectorbound(2)
                                  DO xi = -box_list(ibox)%sectorbound(1), box_list(ibox)%sectorbound(1)
                                          n_adj_cell_atoms(xi,yi,zi,ibox) = SUM(n_cell_atoms( &
                                                  xi-1:xi+1, &
                                                  yi-1:yi+1, &
                                                  zi-1:zi+1, &
                                                  ibox))
                                  END DO
                          END DO
                  END DO
                  !$OMP END DO
                  !$OMP END PARALLEL
          END DO
          !$OMP PARALLEL WORKSHARE
          max_adj_cell_atoms = MAXVAL(n_adj_cell_atoms) ! declare in global variables
          !$OMP END PARALLEL WORKSHARE
          ! ensure 64-byte alignment
          IF (MOD(max_adj_cell_atoms,16) .NE. 0) max_adj_cell_atoms = (max_adj_cell_atoms/16+1)*16
          IF (ALLOCATED(adj_cell_rsp)) THEN
                  IF (SIZE(adj_cell_rsp,1) < max_adj_cell_atoms .OR. ANY(adj_cellmaxbound>adj_cellmaxbound_old)) THEN
                          DEALLOCATE(adj_cell_rsp)
                          DEALLOCATE(adj_cell_ti)
                  END IF
          END IF
          IF (.NOT. ALLOCATED(adj_cell_rsp)) THEN
                  ALLOCATE(adj_cell_rsp(max_adj_cell_atoms,3, &
                          -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                          -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                          -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                          nbr_boxes))
                  ALLOCATE(adj_cell_ti(max_adj_cell_atoms, &
                          -adj_cellmaxbound(1):adj_cellmaxbound(1), &
                          -adj_cellmaxbound(2):adj_cellmaxbound(2), &
                          -adj_cellmaxbound(3):adj_cellmaxbound(3), &
                          nbr_boxes))
          END IF
          DO ibox = 1, nbr_boxes
          !$OMP PARALLEL PRIVATE(istart,iend,xi2,yi2,zi2,dzi,dyi,dxi,vlen,i)
          !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
          DO zi = -box_list(ibox)%sectorbound(3), box_list(ibox)%sectorbound(3)
          DO yi = -box_list(ibox)%sectorbound(2), box_list(ibox)%sectorbound(2)
          DO xi = -box_list(ibox)%sectorbound(1), box_list(ibox)%sectorbound(1)
                istart = 1
                DO dzi = -1, 1
                        zi2 = zi+dzi
                        DO dyi = -1, 1
                                yi2 = yi + dyi
                                DO dxi = -1, 1
                                        xi2 = xi + dxi
                                        vlen = n_cell_atoms(xi2,yi2,zi2,ibox)
                                        IF (vlen < 1) CYCLE
                                        iend = istart + vlen - 1
                                        adj_cell_rsp(istart:iend,1:3,xi,yi,zi,ibox) = &
                                                this_cell_rsp(1:vlen,1:3,xi2,yi2,zi2,ibox)
                                        adj_cell_ti(istart:iend,xi,yi,zi,ibox) = &
                                                this_cell_ti(1:vlen,xi2,yi2,zi2,ibox)
                                        istart = istart + vlen
                                END DO
                        END DO
                END DO
          END DO
          END DO
          END DO
          !$OMP END DO
          !$OMP END PARALLEL
          END DO
  END SUBROUTINE Sector_Setup

  LOGICAL FUNCTION check_overlap(ia, im, is)
          !
          INTEGER, INTENT(IN) :: ia, im, is
          REAL(SP), DIMENSION(max_adj_cell_atoms) :: rminsq_list
          INTEGER, DIMENSION(3) :: cp, sf
          REAL(DP) :: dprp(3)
          REAL(SP) :: irp(3), dxp, dyp, dzp, rijsq, rminsq, sprp
          INTEGER :: this_box
          INTEGER :: i
          INTEGER :: vlen
          INTEGER :: xi, yi, zi
          INTEGER :: ti_solute, icp, sb, lc
          LOGICAL :: overlap, this_overlap

          !DIR$ ASSUME_ALIGNED rminsq_list:32

          check_overlap = .FALSE.
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
          DO i = 1, vlen
                dxp = adj_cell_rsp(i,1,xi,yi,zi,this_box) - irp(1)
                dyp = adj_cell_rsp(i,2,xi,yi,zi,this_box) - irp(2)
                dzp = adj_cell_rsp(i,3,xi,yi,zi,this_box) - irp(3)
                rijsq = dxp*dxp
                rijsq = rijsq + dyp*dyp
                rijsq = rijsq + dzp*dzp
                IF (read_atompair_rminsq) THEN
                        rminsq = sp_atompair_rminsq_table(adj_cell_ti(i,xi,yi,zi,this_box),ti_solute,this_box)
                ELSE IF (calc_rmin_flag) THEN
                        rminsq = sp_rminsq_table(adj_cell_ti(i,xi,yi,zi,this_box),ti_solute)
                END IF
                this_overlap = rijsq < rminsq
                !IF (read_atompair_rminsq .OR. calc_rmin_flag) THEN
                !        this_overlap = rijsq < rminsq_list(i)
                !ELSE
                !        this_overlap = rijsq < sp_rcut_lowsq
                !END IF
                overlap = overlap .OR. this_overlap
          END DO
          check_overlap = overlap
          
  END FUNCTION check_overlap

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



END MODULE Sector_Routines
