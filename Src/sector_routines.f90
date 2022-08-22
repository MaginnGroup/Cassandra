
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
          INTEGER, DIMENSION(:), ALLOCATABLE :: xi_pm, yi_pm, zi_pm
          INTEGER :: dummy
          LOGICAL :: asflag

          sectorbound_old = sectorbound
          sectormaxbound_old = sectormaxbound
          nsec_old = MAXVAL(PRODUCT(length_cells,1))
          DO ibox = 1, nbr_boxes
                DO i = 1, 3
                          length_cells(i,ibox) = INT(box_list(ibox)%length(i,i)/max_rmin)
                          IF (MOD(length_cells(i,ibox),2) .EQ. 0) length_cells(i,ibox) = length_cells(i,ibox) - 1
                          cell_length_inv(i,ibox) = REAL(length_cells(i,ibox),DP) / box_list(ibox)%length(i,i)
                          sectorbound(i,ibox) = length_cells(i,ibox)/2
                END DO
          END DO
          ! cell_length = 1.0_DP / cell_length_inv

          sectormaxbound = MAXVAL(sectorbound, 2)
          !map_bound = sectormaxbound*3+1
          IF (.NOT. ALL( sectormaxbound <= sectormaxbound_old)) THEN
                  IF (ALLOCATED(sector_index_map)) DEALLOCATE(sector_index_map)
                  ALLOCATE(sector_index_map(-sectormaxbound(1):sectormaxbound(1), &
                          -sectormaxbound(2):sectormaxbound(2), &
                          -sectormaxbound(3):sectormaxbound(3), &
                          nbr_boxes))
                  IF (ALLOCATED(sector_has_atoms)) DEALLOCATE(sector_has_atoms)
                  ALLOCATE(sector_has_atoms(-sectormaxbound(1):sectormaxbound(1), &
                          -sectormaxbound(2):sectormaxbound(2), &
                          -sectormaxbound(3):sectormaxbound(3), &
                          nbr_boxes))
          END IF
          max_occ_sectors_old = max_occ_sectors
          total_atoms = DOT_PRODUCT(SUM(nmols(:,1:),2), natoms)
          max_occ_sectors=MAX(max_occ_sectors_old,MIN(total_atoms,SUM(PRODUCT(length_cells,1))))
          IF (max_occ_sectors > max_occ_sectors_old) THEN
                  IF (ALLOCATED(sector_n_atoms)) DEALLOCATE(sector_n_atoms)
                  IF (ALLOCATED(sector_atoms)) DEALLOCATE(sector_atoms)
                  ALLOCATE(sector_n_atoms(0:max_occ_sectors))
                  ALLOCATE(sector_atoms(max_sector_natoms,max_occ_sectors,3))
          END IF
          !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
          sector_n_atoms = 0
          n_occ_sectors = 0
          sector_has_atoms = .FALSE.
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
                                        atom_ptr%ci = IDNINT(cp*cell_length_inv(:,ibox))
                                        DO i = 1,3
                                                IF (atom_ptr%ci(i) > sectorbound(i,ibox)) THEN
                                                        atom_ptr%ci(i) = atom_ptr%ci(i) - length_cells(i,ibox)
                                                ELSE IF (atom_ptr%ci(i) < -sectorbound(i,ibox)) THEN
                                                        atom_ptr%ci(i) = atom_ptr%ci(i) + length_cells(i,ibox)
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
                                        ci = atom_list(ia,im,is)%ci
                                        IF (sector_has_atoms(ci(1),ci(2),ci(3),ibox)) THEN
                                                secind = sector_index_map(ci(1),ci(2),ci(3),ibox)
                                                sector_n_atoms(secind) = sector_n_atoms(secind)+1
                                        ELSE
                                                n_occ_sectors = n_occ_sectors+1
                                                secind = n_occ_sectors
                                                sector_has_atoms(ci(1),ci(2),ci(3),ibox) = .TRUE.
                                                sector_n_atoms(secind) = 1
                                                sector_index_map(ci(1),ci(2),ci(3),ibox) = secind
                                        END IF
                                        IF (sector_n_atoms(secind) > max_sector_natoms) asflag = .FALSE.
                                        IF (asflag) sector_atoms(sector_n_atoms(secind),secind,:) = sector_atom_ID
                                END DO AtomLoop2
                        END DO MoleculeLoop2
                END DO SpeciesLoop
          END DO BoxLoop
          IF (asflag) THEN
                EXIT AllocationSizeLoop
          ELSE
                max_sector_natoms = MAXVAL(sector_n_atoms)
                DEALLOCATE(sector_atoms)
                ALLOCATE(sector_atoms(max_sector_natoms,max_occ_sectors,3))
                asflag = .TRUE.
                sector_n_atoms = 0
          END IF
          END DO AllocationSizeLoop
  END SUBROUTINE Sector_Setup

  LOGICAL FUNCTION check_overlap(ia, im, is)
          !
          INTEGER, INTENT(IN) :: ia, im, is
          INTEGER :: this_species, this_molecule, this_atom, this_box
          INTEGER :: this_locate, i_dim, secind
          INTEGER :: xi, yi, zi, i, ia_cell
          INTEGER, DIMENSION(:), POINTER :: sector_atom_ID
          TYPE(Atom_Class), POINTER :: atom_ptr
          INTEGER :: cell_coords(3)
          INTEGER, DIMENSION(3,3) :: ci
          REAL(DP) :: cp(3), dx, dy, dz, dxp, dyp, dzp
          LOGICAL :: need_wrapping(3)
          INTEGER, PARAMETER, DIMENSION(3) :: delta = (/0,-1,1/)
          !
          check_overlap = .TRUE.
          IF (widom_active) THEN
                  cp(1) = widom_atoms(ia)%rxp
                  cp(2) = widom_atoms(ia)%ryp
                  cp(3) = widom_atoms(ia)%rzp
                  this_box = widom_molecule%which_box
          ELSE
                  cp(1) = atom_list(ia,im,is)%rxp
                  cp(2) = atom_list(ia,im,is)%ryp
                  cp(3) = atom_list(ia,im,is)%rzp
                  this_box = molecule_list(im,is)%which_box
          END IF
          cell_coords = IDNINT(cp*cell_length_inv(:,this_box))
          need_wrapping = ABS(cell_coords) .GE. sectorbound(:,this_box)
          DO i_dim = 1,3
                ci(i_dim,:) = cell_coords(i_dim) + delta
                IF (need_wrapping(i_dim)) THEN
                        DO i = 1,3
                                IF (ci(i_dim,i) > sectorbound(i_dim,this_box)) THEN
                                        ci(i_dim,i) = ci(i_dim,i) - length_cells(i_dim,this_box)
                                ELSE IF (ci(i_dim,i) < -sectorbound(i_dim,this_box)) THEN
                                        ci(i_dim,i) = ci(i_dim,i) + length_cells(i_dim,this_box)
                                END IF
                        END DO
                END IF
          END DO
          DO xi = 1,3
                DO yi = 1,3
                        DO zi = 1,3
                                IF (.NOT. sector_has_atoms(ci(1,xi),ci(2,yi),ci(3,zi),this_box)) CYCLE
                                secind = sector_index_map(ci(1,xi),ci(2,yi),ci(3,zi),this_box)
                                DO ia_cell = 1, sector_n_atoms(secind)
                                        sector_atom_ID => sector_atoms(ia_cell,secind,:)
                                        atom_ptr => atom_list(sector_atom_ID(1),sector_atom_ID(2),sector_atom_ID(3))
                                        dxp = atom_ptr%rxp - cp(1)
                                        dyp = atom_ptr%ryp - cp(2)
                                        dzp = atom_ptr%rzp - cp(3)
                                        CALL Minimum_Image_Separation(this_box,dxp,dyp,dzp,dx,dy,dz)
                                        IF (dx*dx+dy*dy+dz*dz < rminsq_table(nonbond_list(sector_atom_ID(1),sector_atom_ID(3))%atom_type_number,nonbond_list(ia,is)%atom_type_number)) RETURN
                                END DO
                        END DO
                END DO
          END DO
          check_overlap = .FALSE.
          
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
