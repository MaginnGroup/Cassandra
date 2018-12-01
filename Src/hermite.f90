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
SUBROUTINE Hermite(this_box,xcoord_in,ycoord_in,zcoord_in,this_type,coord_potential,overlap)
  !------------------------------------------------------
  !
  ! This subroutine takes in the coordinates of an atom
  ! and returns the potential values associated with
  ! the coordinates. The code is based on that developed
  ! by R. L. June and later modified by Ed Maginn and
  ! Mike Macedonia
  !
  !-----------------------------------------------------

  USE Global_Variables
  USE Hermite_Routines
  USE Potential_MAP_Variables
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_type, this_box
  REAL(DP), INTENT(IN) :: xcoord_in, ycoord_in, zcoord_in
  REAL(DP), INTENT(OUT) :: coord_potential
  LOGICAL, INTENT(INOUT) :: overlap

  INTEGER :: itype

  INTEGER :: ii, jj, kk, ii_next, jj_next, kk_next
  INTEGER :: mi(8), i, m, common_offset

!  REAL(DP) :: h0, hd0, h1, hd1
  REAL(DP) :: h00, h01, h10, h11
  REAL(DP) :: xcoord, ycoord, zcoord
  REAL(DP) :: sum1, t, u, v
  REAL(DP) :: hi(64), hj(64), hk(64), ft(64)
  REAL(DP) :: hx0(8), hx1(8), hy0(8), hy1(8), hz0(8), hz1(8)

  REAL(DP) :: x_next, y_next, z_next
  REAL(DP) :: sx, sy, sz

  ! Initialize output from the subroutine
  coord_potential = 0.0_DP
  overlap = .FALSE.


  ! obtain the leftmost grid corner of the cube that would
  ! ultimately encapsulate the point
  
  ! The input x,y and z coordinates are between (-lx/2,lx/2), (-ly/2,ly/2) and (-lz/2,lz/2)
  ! First move the coordinates so that they are between (0,lx), (0,ly) and (0,lz). At this
  ! point it is possible that an atomic coordinate is still less than 0 or greather than
  ! the length in a given direction. This is due to the fact that only the molecule COM
  ! is guaranteed to be inside the simulation box. A part of molecule may be "hanging out"
  ! of the simulation box.

  xcoord = xcoord_in
  ycoord = ycoord_in
  zcoord = zcoord_in

  ! Transform this into fractional coordinates.
  CALL Cartesian_To_Fractional(xcoord,ycoord,zcoord,sx,sy,sz,this_box)

  ! Now obtain the grid point for this point. Note that the code will return the lower corner of the cube
  ! vertices of which will be used for compute the potential energy.
  ! The eight vertices are: (ii,jj,kk)
  !                         (ii+1,jj,kk), (ii,jj+1,kk), (ii,jj,kk+1)
  !                         (ii+1,jj+1,kk), (ii+1,jj,kk+1),(ii,jj+1,kk+1)
  !                         (ii+1,jj+1,kk+1)
  ! Ensure that the fractional coordinates are within the simulation box
  CALL Apply_PBC_For_Hermite(sx,sy,sz,this_box,ii,jj,kk)

  ! xstep, ystep, zstep are grid spacings in the x, y and z directions

!!$  ii = INT(xcoord/xstep) + 1
!!$  jj = INT(ycoord/ystep) + 1
!!$  kk = INT(zcoord/zstep) + 1
!!$
!!$  IF (ii < 1) ii = 1
!!$  IF (jj > nygrid) jj = nygrid
!!$  IF (jj < 1) jj = 1
!!$  IF (kk > nzgrid) kk = nzgrid
!!$  IF (kk < 1) kk = 1

  ! obtain the indices of upper rightmost corner of the cube
 ! write(*,*) xcoord_in,ycoord_in,zcoord_in
  
  IF (ii == na_grid) THEN
     ii_next = 1
     x_next = a_grid(ii) + a_spacing
  ELSE
     ii_next = ii + 1
     x_next = a_grid(ii_next)
  END IF

  IF (jj == nb_grid) THEN
     jj_next = 1
     y_next = b_grid(jj) + b_spacing
  ELSE
     jj_next = jj + 1
    y_next = b_grid(jj_next)
  END IF

  IF (kk == nc_grid) THEN
     kk_next = 1
     z_next = c_grid(kk) + c_spacing
  ELSE
     kk_next = kk + 1
     z_next = c_grid(kk_next)
  END IF

  ! find out the index of the atom type as if there were no zeolite
  itype = sorbate_atomtypes(this_type)

  ! now check the grid ids of all the eight points
  common_offset = na_grid * nb_grid * nc_grid * (itype - 1) 
!  write(*,*) na_grid, nb_grid,nc_grid
!  write(*,*) common_offset, itype, ii, jj, kk
!  write(*,*) ii_next, jj_next, kk_next
!!$  write(*,*) xcoord, ycoord, zcoord
!!$  write(*,*) a_spacing, ystep, zstep

  ! ii, jj, kk
  mi(1) = sorb_grid_pointer(common_offset + (ii-1) * nb_grid * nc_grid + &
       (jj-1) * nc_grid + kk )
  
  ! ii, jj, kk_next
  mi(4) = sorb_grid_pointer(common_offset + (ii-1) * nb_grid * nc_grid + &
       (jj-1) * nc_grid + kk_next)
  
  ! ii, jj_next, kk
  mi(3) = sorb_grid_pointer(common_offset + (ii-1) * nb_grid * nc_grid + &
       (jj_next-1) * nc_grid + kk)
  
  ! ii_next, jj, kk
 ! write(*,*) (ii_next-1) * nb_grid*nc_grid, (jj-1) * nc_grid
  mi(2) = sorb_grid_pointer(common_offset + (ii_next-1) * nb_grid * nc_grid + &
       (jj-1) * nc_grid + kk)
  
  ! ii, jj_next, kk_next
  mi(7) = sorb_grid_pointer(common_offset + (ii-1) * nb_grid * nc_grid + &
       (jj_next-1) * nc_grid + kk_next)
  
  ! ii_next, jj, kk_next
  mi(6) = sorb_grid_pointer(common_offset + (ii_next-1) * nb_grid * nc_grid + &
       (jj-1) * nc_grid + kk_next)
  
  ! ii_next, jj_next, kk
  mi(5) = sorb_grid_pointer(common_offset + (ii_next-1) * nb_grid * nc_grid + &
       (jj_next-1) * nc_grid + kk)
  
  ! ii_next, jj_next, kk_next
  mi(8) = sorb_grid_pointer(common_offset + (ii_next-1) * nb_grid * nc_grid + &
       (jj_next-1) * nc_grid + kk_next)

  ! Find out if any of the grid indices is 0. If so, then the point
  ! is high energy point and overlap is set to true

  DO i = 1, 8
     IF ( mi(i) == 0) THEN
        overlap = .TRUE.
!        write(*,*) overlap
        RETURN
     END IF
  END DO

!!$  IF ( .NOT. CBMC_Flag) THEN
!!$
!!$        WRITE(106,*)( lattice_pot(common_offset+mi(i)), i = 1, 8)
!!$
!!$  END IF

  ! if here then obtain the potential and its derivatives
  ! at all the eight points. The values will be stored in
  ! 64 elemements of the 'ft' array

  common_offset = SUM(n_sorb_grids(1:itype-1))

  DO i = 1, 8
     m = (i-1) * 8
     ft(m+1) = lattice_pot(common_offset + mi(i))
!    write(*,*) m+1, ft(m+1)
     ft(m+2) = lattice_potx(common_offset + mi(i))
     ft(m+3) = lattice_poty(common_offset + mi(i))
     ft(m+4) = lattice_potz(common_offset + mi(i))
     ft(m+5) = lattice_potxy(common_offset + mi(i))
     ft(m+6) = lattice_potxz(common_offset + mi(i))
     ft(m+7) = lattice_potyz(common_offset + mi(i))
     ft(m+8) = lattice_potxyz(common_offset + mi(i))
    
  END DO
!stop
  ! scaled coordinates

 ! write(*,*) 

  t = (sx - a_grid(ii))/abs((x_next - a_grid(ii)))
  u = (sy - b_grid(jj))/abs((y_next - b_grid(jj)))
  v = (sz - c_grid(kk))/abs((z_next - c_grid(kk)))

  IF (abs(sx) > 0.5_DP) THEN
     write(*,*) i_mcstep, 't', sx, a_grid(ii), x_next
  END IF
  IF (abs(sy) > 0.5_DP) THEN
     write(*,*) 'u', sy, b_grid(jj), y_next
  END IF
  IF (abs(sz) > 0.5_DP) THEN
     write(*,*) 'v', sz, c_grid(kk), z_next
  END IF
  
  !look up the basis function values at the input coordinates

  h00=h0(1.0_DP-t, a_spacing)
  h01=h0(    t, a_spacing)
  h10=h1(1.0_DP-t,-a_spacing)
  h11=h1(    t, a_spacing)

 !  if ( i_mcstep == 16) write(*,*) 'hspacing', h00, h01, h10, h11

  hx0(1)=h00
  hx0(2)=h01
  hx0(3)=h00
  hx0(4)=h00
  hx0(5)=h01
  hx0(6)=h01
  hx0(7)=h00
  hx0(8)=h01
  
  hx1(1)=h10
  hx1(2)=h11
  hx1(3)=h10
  hx1(4)=h10
  hx1(5)=h11
  hx1(6)=h11
  hx1(7)=h10
  hx1(8)=h11
 

  DO i=1,8
     m=(i-1)*8
     hi(m+1)=hx0(i)
     hi(m+2)=hx1(i)
     hi(m+3)=hx0(i)
     hi(m+4)=hx0(i)
     hi(m+5)=hx1(i)
     hi(m+6)=hx1(i)
     hi(m+7)=hx0(i)
     hi(m+8)=hx1(i)
  END DO

  ! do the y coordinate
  
  h00=h0(1.0_DP-u, b_spacing)
  h01=h0(    u, b_spacing)
  h10=h1(1.0_DP-u,-b_spacing)
  h11=h1(    u, b_spacing)

!  if ( i_mcstep == 16) write(*,*) 'hbspacing', h00, h01, h10, h11
!  write(*,*) 'u', u
!  write(*,*) h00, h01, h10, h11
  
  hy0(1)=h00
  hy0(2)=h00
  hy0(3)=h01
  hy0(4)=h00
  hy0(5)=h01
  hy0(6)=h00
  hy0(7)=h01
  hy0(8)=h01
  
  hy1(1)=h10
  hy1(2)=h10
  hy1(3)=h11
  hy1(4)=h10
  hy1(5)=h11
  hy1(6)=h10
  hy1(7)=h11
  hy1(8)=h11
  
  DO  i=1,8
     m=(i-1)*8
     hj(m+1)=hy0(i)
     hj(m+2)=hy0(i)
     hj(m+3)=hy1(i)
     hj(m+4)=hy0(i)
     hj(m+5)=hy1(i)
     hj(m+6)=hy0(i)
     hj(m+7)=hy1(i)
     hj(m+8)=hy1(i)
  END DO
  
  ! do the z coordinate basis functions
  
  h00=h0(1.0_DP-v, c_spacing)
  h01=h0(    v, c_spacing)
  h10=h1(1.0_DP-v,-c_spacing)
  h11=h1(    v, c_spacing)
!  if ( i_mcstep == 16) write(*,*) 'hcspacing', h00, h01, h10, h11
!  write(*,*) 'v'
!  write(*,*) h00, h01, h10, h11
  hz0(1)=h00
  hz0(2)=h00
  hz0(3)=h00
  hz0(4)=h01
  hz0(5)=h00
  hz0(6)=h01
  hz0(7)=h01
  hz0(8)=h01
  
  hz1(1)=h10
  hz1(2)=h10
  hz1(3)=h10
  hz1(4)=h11
  hz1(5)=h10
  hz1(6)=h11
  hz1(7)=h11
  hz1(8)=h11
 
  DO i=1,8
     m=(i-1)*8
     hk(m+1)=hz0(i)
     hk(m+2)=hz0(i)
     hk(m+3)=hz0(i)
     hk(m+4)=hz1(i)
     hk(m+5)=hz0(i)
     hk(m+6)=hz1(i)
     hk(m+7)=hz1(i)
     hk(m+8)=hz1(i)
  END DO


  sum1=0.0_DP
  
  DO i=1,64
  !   write(*,'(I5,2X,4(F10.5,2X))') i, hi(i)*hj(i)*hk(i), hi(i), hj(i), hk(i)
     sum1=sum1+ft(i)*hi(i)*hj(i)*hk(i)
  END DO
 
  coord_potential = sum1

!  if (i_mcstep == 16) WRITE(*,*) 'Hermite', sum1

!!$ IF ( .NOT. CBMC_Flag) THEN
 !!$
!!$     write(*,*) t, u, v
!!$     IF ( t < 0.0_DP) THEN
!!$        write(*,*) 't', xcoord, x_grid(ii), x_grid(ii_next)
!!$     END IF
!!$
!!$     IF ( u < 0.0_DP) THEN
!!$        write(*,*) 't', ycoord, y_grid(ii), y_grid(ii_next)
!!$     END IF
!!$     IF ( v < 0.0_DP) THEN
!!$        write(*,*) 't', zcoord, z_grid(ii), z_grid(ii_next)
!!$     END IF
!!$     write(102,*) x_grid(ii), y_grid(jj), z_grid(kk)
!!$     write(103,*) x_grid(ii_next), y_grid(jj_next), z_grid(kk_next)
!!$     write(104,*) ii, jj, kk, ii_next, jj_next, kk_next
!!$     write(105,*) coord_potential
!!$
!!$  END IF

!  IF (coord_potential <= 0.0_DP) THEN
!     write(*,*) t, u, v, coord_potential
!     write(*,*) xcoord_in, ycoord_in, zcoord_in
!     write(*,*) x_next, a_grid(ii)
!     write(*,*) y_next, b_grid(jj)
!     write(*,*) z_next, c_grid(kk)
     
!  END IF
  
END SUBROUTINE hermite

