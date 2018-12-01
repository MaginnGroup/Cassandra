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
!********************************************************************************
MODULE Hermite_Routines
  !---------------------------------------------------------------------------
  !
  ! This module contains various functions and PBC conditions for use
  ! in Hermite Interpolation
  !
  !----------------------------------------------------------------------------

  IMPLICIT NONE

CONTAINS


  FUNCTION h0(xx,ds)

    USE Type_Definitions, ONLY : DP

    REAL(DP) :: h0, xx, ds

    h0 = xx * xx * (3.0_DP - 2.0_DP * xx)

  END FUNCTION h0

  FUNCTION hd0(xx,ds)

    USE Type_Definitions, ONLY : DP

    REAL(DP) :: hd0, xx, ds

    hd0 = 6.0_DP * xx * (1.0_DP - xx)

  END FUNCTION hd0

  FUNCTION h1(xx,ds)

    USE Type_Definitions, ONLY : DP

    REAL(DP) :: h1, xx, ds

    h1 = ds * xx * xx * (xx - 1.0_DP)

  END FUNCTION h1

  FUNCTION hd1(xx,ds)

    USE Type_Definitions, ONLY : DP

    REAL(DP) :: hd1, xx, ds

    hd1 = ds * xx * (3.0_DP * xx - 2.0_DP)

  END FUNCTION hd1


  SUBROUTINE Apply_PBC_For_Hermite(xthis,ythis,zthis,this_box,nxthis,nythis,nzthis)
    !*****************************************************************************
    !
    ! Apply PBC to atomic coordinates so that they are inside the simulation box
    !
    !*****************************************************************************

    USE Type_Definitions, ONLY: DP
    USE Potential_Map_Variables, ONLY : a_spacing, b_spacing, c_spacing
!    l_slit_pore, pore_width, half_pore_width, &
!         xstep, ystep, zstep, zeo_unit_cell

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: this_box
    INTEGER, INTENT(OUT) :: nxthis, nythis, nzthis

    REAL(DP), INTENT(INOUT) :: xthis, ythis, zthis

    REAL(DP) :: xtemp, ytemp, ztemp


    xtemp = xthis
    ytemp = ythis
    ztemp = zthis
    ! transfor the fractional coordinates so that they are inside the box.
    ! Note only the molecule COM is guaranteed to be inside the simulation box
    ! A part of the molecule may be dangling outside.


    DO WHILE (xtemp < -0.5_DP)
       xtemp = xtemp + 1.0_DP
    END DO
    
    DO WHILE (xtemp > 0.5_DP)
       xtemp = xtemp - 1.0_DP
    END DO
  

    DO WHILE (ytemp < -0.5_DP)
       ytemp = ytemp + 1.0_DP
    END DO
    
    DO WHILE (ytemp > 0.5_DP)
       ytemp = ytemp - 1.0_DP
    END DO
    
    DO WHILE (ztemp < -0.5_DP)
       ztemp = ztemp + 1.0_DP
    END DO
        
    DO WHILE (ztemp > 0.5_DP)
       ztemp = ztemp - 1.0_DP
    END DO

    nxthis = INT((xtemp+0.5_DP)/a_spacing) + 1
    nythis = INT((ytemp+0.5_DP)/b_spacing) + 1
    nzthis = INT((ztemp+0.5_DP)/c_spacing) + 1
    
    xthis = xtemp
    ythis = ytemp
    zthis = ztemp

!!$ 
!!$
!!$    ! transfrom the x,y and z coordinates so that the center of
!!$    ! the box is at (lx/2,ly/2,lz/2)
!!$
!!$    xtemp = xthis + box_list(this_box)%hlength(1,1)
!!$    ytemp = ythis + box_list(this_box)%hlength(2,2)
!!$
!!$    IF (l_slit_pore(this_box)) THEN
!!$       ztemp = zthis + half_pore_width
!!$    ELSE
!!$       ztemp = zthis + box_list(this_box)%hlength(3,3)
!!$    END IF
!!$
!!$
!!$    ! pbc in x direction. We will use the zeolite unit
!!$    ! cell to transform back to the actual zeolite unit cell
!!$
!!$    DO WHILE (xtemp < 0.0_DP) 
!!$       xtemp = xtemp + zeo_unit_cell%length(1,1)
!!$    END DO
!!$
!!$    DO WHILE (xtemp > zeo_unit_cell%length(1,1))
!!$       xtemp = xtemp - zeo_unit_cell%length(1,1)
!!$    END DO
!!$
!!$
!!$    IF (xtemp < 0.0_DP) THEN
!!$       xtemp = xtemp + box_list(this_box)%length(1,1)
!!$    ELSE IF (xtemp > box_list(this_box)%length(1,1)) THEN
!!$       xtemp = xtemp - box_list(this_box)%length(1,1)
!!$    END IF

!!$    ! pbc in y direction
!!$
!!$    DO WHILE (ytemp < 0.0_DP) 
!!$       ytemp = ytemp + zeo_unit_cell%length(2,2)
!!$    END DO
!!$
!!$    DO WHILE (ytemp > zeo_unit_cell%length(2,2))
!!$       ytemp = ytemp - zeo_unit_cell%length(2,2)
!!$    END DO
!!$
!!$ 
!!$    
!!$    IF (ytemp < 0.0_DP) THEN
!!$       ytemp = ytemp + box_list(this_box)%length(2,2)
!!$    ELSE IF (ytemp > box_list(this_box)%length(2,2)) THEN
!!$       ytemp = ytemp - box_list(this_box)%length(2,2)
!!$    END IF
!!$
!!$    ! pbc in z direction
!!$    IF (l_slit_pore(this_box)) THEN
!!$
!!$       IF (ztemp < 0.0_DP) THEN
!!$          ztemp = ztemp + pore_width
!!$       ELSE IF (ztemp > pore_width ) THEN
!!$          ztemp = ztemp - pore_width
!!$       END IF
!!$
!!$    ELSE
!!$
!!$       DO WHILE (ztemp < 0.0_DP) 
!!$          ztemp = ztemp + zeo_unit_cell%length(3,3)
!!$       END DO
!!$       
!!$       DO WHILE (ztemp > zeo_unit_cell%length(3,3))
!!$          ztemp = ztemp - zeo_unit_cell%length(3,3)
!!$       END DO
!!$
!!$    END IF



!!$       IF (ztemp < 0.0_DP) THEN
!!$          ztemp = ztemp + box_list(this_box)%length(3,3)
!!$       ELSE IF (ztemp > box_list(this_box)%length(3,3)) THEN
!!$          ztemp = ztemp - box_list(this_box)%length(3,3)
!!$       END IF
!!$
!!$    END IF

    ! return the grid positions

!!$    xtemp = xstep * 10.0_DP
!!$    ytemp = ystep * 10.0_DP
!!$    ztemp = zstep * 10.0_DP
    
!!$    nxthis = INT(xtemp/xstep) + 1
!!$    nythis = INT(ytemp/ystep) + 1
!!$    nzthis = INT(ztemp/zstep) + 1
!!$
!!$    ! also return the coordinates in the new frame of reference
!!$
!!$    xthis = xtemp
!!$    ythis = ytemp
!!$    zthis = ztemp

  END SUBROUTINE Apply_PBC_For_Hermite

END MODULE Hermite_Routines

