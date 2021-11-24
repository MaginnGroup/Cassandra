
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
          INTEGER :: i_sector, ci(3), dx, dy, dz, xshift, yshift, zshift, nsec_old, nsec, secind
          INTEGER :: sector_ID
          TYPE(ID_Class) :: sector_atom_ID
          REAL(DP) :: xp, yp, zp, cw(3)

          
          sectorbound_old = sectorbound
          sectormaxbound_old = sectormaxbound
          nsec_old = MAXVAL(PRODUCT(length_cells,1))
          DO ibox = 1, nbr_boxes
                DO i = 1, 3
                          length_cells(i,ibox) = INT(box_list(ibox)%length(i,i)*rcut_low_inv)
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
                          -sectormaxbound(3):sectormaxbound(3)))
                  i_sector = 0
                  DO xi = -sectormaxbound(1), sectormaxbound(1)
                          DO yi = -sectormaxbound(2), sectormaxbound(2)
                                  DO zi = -sectormaxbound(3), sectormaxbound(3)
                                          i_sector = i_sector+1
                                          sector_index_map(xi, yi, zi) = i_sector
                                  END DO
                          END DO
                  END DO
                  nsec = i_sector
                  IF (ALLOCATED(sector_ID_list)) DEALLOCATE(sector_ID_list)
                  IF (ALLOCATED(sector_has_atoms)) DEALLOCATE(sector_has_atoms)
                  ALLOCATE(sector_ID_list(nsec,nbr_boxes))
                  ALLOCATE(sector_has_atoms(nsec,nbr_boxes))
          END IF
          sector_ID_list = -1
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
                                        CALL Minimum_Image_Separation(ibox,xp,yp,zp,cw(1),cw(2),cw(3))
                                        ci = IDNINT(cw*cell_length_inv(:,ibox))
                                        secind = sector_index_map(ci(1),ci(2),ci(3))
                                        !$OMP CRITICAL
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
                                        !$OMP END CRITICAL
                                END DO AtomLoop
                        END DO MoleculeLoop
                END DO SpeciesLoop
          END DO BoxLoop



  END SUBROUTINE Sector_Setup

  SUBROUTINE Check_Overlap(xp, yp, zp, this_box, overlap_flag)
          !
          INTEGER :: this_species, this_molecule, this_atom, this_box
          INTEGER :: this_locate, i_dim
          INTEGER :: xi, yi, zi
          INTEGER(3) :: ci, cell_coords
          REAL(DP) :: xp, yp, zp, cw(3)
          LOGICAL :: overlap_flag
          !
          overlap_flag = .TRUE.
          CALL Minimum_Image_Separation(this_box,xp,yp,zp,cw(1),cw(2),cw(3))
          cell_coords = IDNINT(cw*cell_length_inv(:,ibox))
          DO xi = cell_coords(1)-1, cell_coords(1)+1
                IF (xi < -sectorbound(1, this_box)) THEN
                        ci(1) = xi + length_cells(1, this_box)
                ELSE IF (xi > sectorbound(1, this_box)) THEN
                        ci(1) = xi - length_cells(1, this_box)
                END IF
                DO yi = cell_coords(2)-1, cell_coords(2)+1
                        IF (yi < -sectorbound(2, this_box)) THEN
                                ci(2) = yi + length_cells(2, this_box)
                        ELSE IF (yi > sectorbound(2, this_box)) THEN
                                ci(2) = yi - length_cells(2, this_box)
                        END IF
                        DO zi = cell_coords(3)-1, cell_coords(3)+1
                                IF (zi < -sectorbound(3, this_box)) THEN
                                        ci(3) = zi + length_cells(3, this_box)
                                ELSE IF (zi > sectorbound(3, this_box)) THEN
                                        ci(3) = zi - length_cells(3, this_box)
                                END IF
                                secind = sector_index_map(ci(1),ci(2),ci(3),ibox)
                                IF (.NOT. sector_has_atoms(secind,ibox)) CYCLE
                                sector_ID = sector_ID_list(secind, ibox)
                                DO ia_cell = 1, sector_n_atoms(sector_ID)
                                        this_locate = sector_atoms(ia_cell,sector_ID)%mol
                                        this_species = sector_atoms(ia_cell,sector_ID)%spec
                                        this_atom = sector_atoms(ia_cell,sector_ID)%atom
                                        dxp = atom_list(this_atom,this_locate,this_species)%rxp - xp
                                        dyp = atom_list(this_atom,this_locate,this_species)%ryp - yp
                                        dzp = atom_list(this_atom,this_locate,this_species)%rzp - zp
                                        CALL Minimum_Image_Separation(this_box,dxp,dyp,dzp,dx,dy,dz)
                                        IF (dx*dx+dy*dy+dz*dz<rcut_lowsq) RETURN
                        END DO
                END DO
          END DO
          overlap_flag = .FALSE.
          
  END SUBROUTINE Check_Overlap

  SUBROUTINE Add_To_Sector(this_species, this_molecule, this_atom, this_box)
          RETURN
  END SUBROUTINE Add_To_Sector

END MODULE IO_Utilities
