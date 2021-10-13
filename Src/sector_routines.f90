
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

  IMPLICIT NONE

  INTEGER :: i, ibox, is, imol, im, ia
  


CONTAINS
  SUBROUTINE Sector_Setup
          INTEGER, DIMENSION(3) :: sectormaxbound_old !, map_bound
          INTEGER, DIMENSION(3,nbr_boxes) :: sectorbound_old
          INTEGER :: i_sector, xi, yi, zi, dx, dy, dz, xshift, yshift, zshift, nsec_old, nsec, secind
          INTEGER :: sector_ID
          TYPE(ID_Class) :: sector_atom_ID
          REAL(DP) :: xp, yp, zp, xw, yp, zp

          
          sectorbound_old = sectorbound
          sectormaxbound_old = sectormaxbound
          nsec_old = MAXVAL(PRODUCT(length_cells,1))
          DO ibox = 1, nbr_boxes
                DO i = 1, 3
                        sectorbound(i,ibox) = IDNINT(box_list(ibox)%hlength(i,i)*rcut_low_inv)
                END DO
          END DO
          length_cells = sectorbound*2+1
          sectormaxbound = MAXVAL(sectorbound, 2)
          !map_bound = sectormaxbound*3+1
          IF (.NOT. ALL(sectorbound == sectorbound_old)) THEN
                  IF (.NOT. ALL( sectormaxbound <= sectormaxbound_old)) THEN
                          IF (ALLOCATED(sector_index_map)) DEALLOCATE(sector_index_map)
                          ALLOCATE(sector_index_map(-sectormaxbound(1):sectormaxbound(1), &
                                  -sectormaxbound(2):sectormaxbound(2), &
                                  -sectormaxbound(3):sectormaxbound(3), &
                                  nbr_boxes))
                  END IF
                  DO ibox = 1, nbr_boxes
                        i_sector = 0
                        DO xi = -sectorbound(1,ibox), sectorbound(1,ibox)
                                DO yi = -sectorbound(2,ibox), sectorbound(2,ibox)
                                        DO zi = -sectorbound(3,ibox), sectorbound(3,ibox)
                                                i_sector = i_sector+1
                                                sector_index_map(xi, yi, zi, ibox) = i_sector
                                        END DO
                                END DO
                        END DO
                  END DO
          END IF
          nsec = MAXVAL(PRODUCT(length_cells,1))
          IF (nsec > nsec_old) THEN
                IF (ALLOCATED(sector_ID_list)) DEALLOCATE(sector_ID_list)
                IF (ALLOCATED(sector_has_atoms)) DEALLOCATE(sector_has_atoms)
                ALLOCATE(sector_ID_list(nsec,nbr_boxes))
                ALLOCATE(sector_has_atoms(nsec,nbr_boxes))
          END IF
          n_occ_sectors = 0
          sector_has_atoms = .FALSE.
          !sector_n_atoms = 0
          ! place atoms in sectors
          BoxLoop:DO ibox = 1, nbr_boxes
                SpeciesLoop:DO is = 1, nspecies
                        sector_atom_ID%spec = is
                        MoleculeLoop:DO imol = 1, nmols(is,ibox)
                                im = locate(imol, is, ibox)
                                sector_atom_ID%mol = im
                                IF (.NOT. molecule_list(im,is)%live) CYCLE MoleculeLoop
                                AtomLoop:DO ia = 1, natoms(is)
                                        sector_atom_ID%atom = ia
                                        xp = atom_list(ia,im,is)%rxp
                                        yp = atom_list(ia,im,is)%ryp
                                        zp = atom_list(ia,im,is)%rzp
                                        CALL Minimum_Image_Separation(ibox,xp,yp,zp,xw,yw,zw)
                                        xi = IDNINT(xw)
                                        yi = IDNINT(yw)
                                        zi = IDNINT(zw)
                                        secind = sector_index_map(xi,yi,zi,ibox)
                                        IF (sector_has_atoms(secind,ibox)) THEN
                                                sector_ID = sector_ID_list(secind,ibox)
                                                sector_n_atoms(sector_ID) = sector_n_atoms(sector_ID)+1
                                        ELSE
                                                n_occ_sectors = n_occ_sectors+1
                                                sector_ID = n_occ_sectors
                                                sector_has_atoms(secind,ibox) = .TRUE.
                                                sector_n_atoms(sector_ID) = 1
                                                sector_ID_list(secind,ibox) = sector_ID
                                        END IF
                                        sector_atoms(sector_n_atoms(sector_ID),sector_ID) = sector_atom_ID
                                END DO AtomLoop
                        END DO MoleculeLoop
                END DO SpeciesLoop
          END DO BoxLoop



  END SUBROUTINE Sector_Setup

  SUBROUTINE Check_Overlap(this_species, this_molecule, this_atom, this_box, overlap_flag)
          !
          INTEGER :: this_species, this_molecule, this_atom, this_box
          INTEGER :: this_locate, i_dim
          INTEGER :: xi, yi, zi, dxi, dyi, dzi
          INTEGER(3) :: indstart, indend, cell_coords
          REAL(DP) :: xp, yp, zp, xw, yw, zw
          LOGICAL :: overlap_flag
          !
          indstart = -1
          indend = 1

          this_locate = locate(this_molecule, this_species, this_box)
          xp = atom_list(this_atom, this_locate, this_species)%rxp
          yp = atom_list(this_atom, this_locate, this_species)%ryp
          zp = atom_list(this_atom, this_locate, this_species)%rzp
          CALL Minimum_Image_Separation(this_box,xp,yp,zp,xw,yw,zw)
          cell_coords(1) = IDNINT(xw)
          cell_coords(2) = IDNINT(yw)
          cell_coords(3) = IDNINT(zw)
          DO i_dim = 1, 3
                IF (cell_coords(i_dim) == 1 - sectorbound(i_dim, this_box)) THEN
                        indstart(i_dim) = -2
                ELSE IF (cell_coords(i_dim) == sectorbound(i_dim, this_box) - 1) THEN
                        indend(i_dim) = 2
                END IF
          END DO
          DO dxi = indstart(1), indend(1)
                xi = cell_coords(1) + dxi
                IF (xi < -sectorbound(1, this_box)) THEN
                        xi = xi + length_cells(1, this_box)
                ELSE IF (xi > sectorbound(1, this_box)) THEN
                        xi = xi - length_cells(1, this_box)
                END IF
                DO dyi = indstart(2), indend(2)
                        yi = cell_coords(2) + dyi
                        IF (yi < -sectorbound(2, this_box)) THEN
                                yi = yi + length_cells(2, this_box)
                        ELSE IF (yi > sectorbound(2, this_box)) THEN
                                yi = yi - length_cells(2, this_box)
                        END IF
                        DO dzi = indstart(3), indend(3)
                                zi = cell_coords(3) + dzi
                                IF (zi < -sectorbound(3, this_box)) THEN
                                        zi = zi + length_cells(3, this_box)
                                ELSE IF (zi > sectorbound(3, this_box)) THEN
                                        zi = zi - length_cells(3, this_box)
                                END IF

                        END DO
                END DO
          END DO
          
  END SUBROUTINE Check_Overlap

  SUBROUTINE Add_To_Sector(this_species, this_molecule, this_atom, this_box)
          RETURN
  END SUBROUTINE Add_To_Sector

END MODULE IO_Utilities
