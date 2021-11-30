
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
          INTEGER :: sector_atom_ID(3)
          TYPE(Atom_Class), POINTER :: atom_ptr
          REAL(DP) :: xp, yp, zp, cp(3)

          
          sectorbound_old = sectorbound
          sectormaxbound_old = sectormaxbound
          nsec_old = MAXVAL(PRODUCT(length_cells,1))
          DO ibox = 1, nbr_boxes
                DO i = 1, 3
                          length_cells(i,ibox) = INT(box_list(ibox)%length(i,i)/rcut_low)
                          IF (MOD(length_cells(i,ibox),2) .EQ. 0) length_cells(i,ibox) = length_cells(i,ibox) - 1
                          cell_length_inv(i,ibox) = REAL(length_cells(i,ibox),DP) / box_list(ibox)%length(i,i)
                          sectorbound(i,ibox) = length_cells(i,ibox)/2
                END DO
          END DO
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
          !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
          n_occ_sectors = 0
          sector_has_atoms = .FALSE.
          !$OMP END PARALLEL WORKSHARE
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
                                        sector_atoms(sector_n_atoms(secind),secind,:) = sector_atom_ID
                                END DO AtomLoop2
                        END DO MoleculeLoop2
                END DO SpeciesLoop
          END DO BoxLoop
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
          IF (widom_flag) THEN
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
                                        IF (dx*dx+dy*dy+dz*dz<rcut_lowsq) RETURN
                                END DO
                        END DO
                END DO
          END DO
          check_overlap = .FALSE.
          
  END FUNCTION check_overlap

END MODULE Sector_Routines
