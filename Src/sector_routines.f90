
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
          REAL(DP) :: xp, yp, zp, xw, yp, zp

          
          sectorbound_old = sectorbound
          sectormaxbound_old = sectormaxbound
          nsec_old = MAXVAL(PRODUCT(sectorbound*2+1,1))
          DO ibox = 1, nbr_boxes
                DO i = 1, 3
                        sectorbound(i,ibox) = IDNINT(box_list(ibox)%hlength(i,i)*rcut_low_inv)
                END DO
          END DO
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
                        DO dx = -1,1
                                xshift = dx*(2*sectorbound(1,ibox)+1)
                                DO dy = -1,1
                                        yshift = dy*(2*sectorbound(2,ibox)+1)
                                        DO dz = -1,1
                                                IF (dx == 0 .AND. dy == 0 .AND. dz == 0) CYCLE
                                                zshift = dz*(2*sectorbound(3,ibox)+1)
                                                sector_index_map((-sectorbound(1,ibox)+xshift):(sectorbound(1,ibox)+xshift), &
                                                        (-sectorbound(2,ibox)+yshift):(sectorbound(2,ibox)+yshift), &
                                                        (-sectorbound(3,ibox)+zshift):(sectorbound(3,ibox)+zshift), &
                                                        ibox) = sector_index_map(&
                                                        -sectorbound(1,ibox):sectorbound(1,ibox), &
                                                        -sectorbound(2,ibox):sectorbound(2,ibox), &
                                                        -sectorbound(3,ibox):sectorbound(3,ibox), &
                                                        ibox)
                                        END DO
                                END DO
                        END DO
                  END DO
          END IF
          nsec = MAXVAL(PRODUCT(sectorbound*2+1,1))
          IF (nsec > nsec_old) THEN
                IF (ALLOCATED(sector_ID_list)) DEALLOCATE(sector_ID_list)
                IF (ALLOCATED(sector_has_atoms)) DEALLOCATE(sector_has_atoms)
                ALLOCATE(sector_ID_list(nsec,nbr_boxes))
                ALLOCATE(sector_has_atoms(nsec,nbr_boxes))
          END IF
          n_occ_sectors = 0
          sector_has_atoms = .FALSE.
          sector_n_atoms = 0
          ! place atoms in sectors
          BoxLoop:DO ibox = 1, nbr_boxes
                SpeciesLoop:DO is = 1, nspecies
                        MoleculeLoop:DO imol = 1, nmols(is,ibox)
                                im = locate(imol, is, ibox)
                                IF (.NOT. molecule_list(im,is)%live) CYCLE MoleculeLoop
                                DO ia = 1, natoms(is)
                                        xp = atom_list(ia,im,is)%rxp
                                        yp = atom_list(ia,im,is)%ryp
                                        zp = atom_list(ia,im,is)%rzp
                                        CALL Minimum_Image_Separation(ibox,xp,yp,zp,xw,yw,zw)
                                        xi = IDNINT(xw)
                                        yi = IDNINT(yw)
                                        zi = IDNINT(zw)
                                        secind = sector_index_map(xi,yi,zi,ibox)
                                        IF (.NOT. sector_has_atoms(secind,ibox)) THEN
                                                n_occ_sectors = n_occ_sectors+1
                                        END IF
                                END DO
                        END DO
                END DO
          END DO



  END SUBROUTINE Sector_Setup

  SUBROUTINE Add_To_Sector(this_species, this_molecule, this_atom, this_box)
          RETURN
  END SUBROUTINE Add_To_Sector

END MODULE IO_Utilities
