!********************************************************************************
! CASSANDRA - Computational Atomistic Simulation Software at Notre Dame 
! for Research in Academia.
! http://molsim.wiki.zoho.com/
! Copyright (2007) University of Notre Dame.
! Authors: Ed Maginn (ed@nd.edu); Jindal Shah (jshah@nd.edu)
!********************************************************************************
SUBROUTINE Replicate_Unit_Cell
  !********************************************************************************
  !
  ! replicate_unit_cell.f90
  !
  !
  ! The subroutine replicates the central simulation box in three dimensions to construct
  ! a superlattice which will be used to locate neighboring atoms for a given cubelet.
  !
  ! The outputs from the code are:
  ! 1. total number of zeolite supercell atoms
  ! 2. Corresponding atomic types, and
  ! 3. Cartesian coordinates.
  !
  !*******************************************************************************

  USE Global_Variables
  USE Potential_Map_Variables
  USE IO_Utilities
  USE File_Names

  IMPLICIT NONE

  INTEGER :: na_total, nb_total, nc_total, i, j, k, iatom
  INTEGER :: count, zeo_id, this_box, is, imol

  REAL(DP) :: sx_this, sy_this, sz_this, sx, sy, sz
  REAL(DP) :: rx_this, ry_this, rz_this, xlattice, ylattice, zlattice

  ! Store the coordinates of the lattice atoms in a one dimensional array
  is = 1
  this_box = 1
  n_lat_atoms = natoms(is) * nmols_to_read(is,this_box)
  
  
  ALLOCATE(x_lat(n_lat_atoms), y_lat(n_lat_atoms),z_lat(n_lat_atoms))
  
  WRITE(logunit,*) 
  WRITE(logunit,'(A,2X,I5)') 'Number of lattice atoms ', n_lat_atoms
  WRITE(logunit,*)

  
  ! As the input_routines.f90 already checks that the distance between two parallel faces is at least
  ! twice the cutoff, we need to expand the central simulation cell in both the positive and negative
  ! directions by a unit.

  na_total = 3
  nb_total = 3
  nc_total = 3

!!$  IF(l_slit_pore(1)) THEN
!!$     ! Replication is to be carried out only in x and y plane
!!$     nc_total = 1
!!$  END IF

  ! Total number of atoms in the zeolite supercell

  n_superlat_atoms = n_lat_atoms * na_total * nb_total * nc_total

  ! Now determine the species id of the zeolite so that it can be
  ! used to assign atomtypes to the supercell atoms. Note that
  ! this is done only for the potential map generation. For
  ! an actual simulation, we cannot identify the id of the zeolite
  ! (and we don't need to) because the mcf file for the zeolite
  ! is not provided.
!!$
!!$  IF (int_sim_type == sim_pot_map) THEN
!!$
!!$     DO i = 1, nspecies
!!$
!!$        IF (species_list(i)%int_insert == int_noinsert) THEN
!!$
!!$           zeo_id = i
!!$
!!$           EXIT
!!$
!!$        END IF
!!$
!!$     END DO
!!$
!!$     IF (i == (nspecies+1)) THEN
!!$
!!$        ! We did not find a zeolite species
!!$        ! trigger an error
!!$
!!$        err_msg = ''
!!$        err_msg(1) = 'Species id of zeolite cound not be determined'
!!$        CALL Clean_Abort(err_msg,'replicate_unit_cell.f90')
!!$
!!$     END IF
!!$
!!$  END IF

   ALLOCATE(x_super_lat(n_superlat_atoms), Stat = AllocateStatus)

  IF (AllocateStatus /= 0 ) THEN

     err_msg = ''
     err_msg(1) = 'Memory could not be allocated for x_super_lat array'
     CALL Clean_Abort(err_msg,'replicate_unit_cell.f90')

  END IF

  ALLOCATE(y_super_lat(n_superlat_atoms), Stat = AllocateStatus)

  IF (AllocateStatus /= 0 ) THEN

     err_msg = ''
     err_msg(1) = 'Memory could not be allocated for y_super_lat array'
     CALL Clean_Abort(err_msg,'replicate_unit_cell.f90')

  END IF

  ALLOCATE(z_super_lat(n_superlat_atoms), Stat = AllocateStatus)

  IF (AllocateStatus /= 0 ) THEN

     err_msg = ''
     err_msg(1) = 'Memory could not be allocated for z_super_lat array'
     CALL Clean_Abort(err_msg,'replicate_unit_cell.f90')

  END IF

  ! Allocate type_super_lat and charge_super_lat arrays to identify
  ! atom type and charges for the superlattice atoms

  ALLOCATE(type_super_lat(n_superlat_atoms), Stat = AllocateStatus)

  IF (AllocateStatus /= 0 ) THEN

     err_msg = ''
     err_msg(1) = 'Memory could not be allocated for type_super_lat array'
     CALL Clean_Abort(err_msg,'replicate_unit_cell.f90')

  END IF

  ALLOCATE(charge_super_lat(n_superlat_atoms), Stat = AllocateStatus)

  IF (AllocateStatus /= 0 ) THEN

     err_msg = ''
     err_msg(1) = 'Memory could not be allocated for charge_super_lat array'
     CALL Clean_Abort(err_msg,'replicate_unit_cell.f90')

  END IF


  count = 0

  zeo_id = 1

  OPEN(unit=21,file='replicated.xyz')
  WRITE(21,*) n_superlat_atoms
  WRITE(21,*)
  
  DO imol = 1, nmols_to_read(zeo_id,this_box)
     DO iatom = 1, natoms(zeo_id)
     ! obtain the fractional coordinates of this atom
        xlattice = atom_list(iatom,imol,this_box)%rxp
        ylattice = atom_list(iatom,imol,this_box)%ryp
        zlattice = atom_list(iatom,imol,this_box)%rzp
     
        CALL Cartesian_To_Fractional(xlattice,ylattice,zlattice,sx,sy,sz,this_box)

!!$        IF (sx < 0.0_DP .OR. sx > 1.0_DP) THEN
!!$           write(*,*) 'sx exceeded'
!!$           stop
!!$        END IF
!!$         IF (sy < 0.0_DP .OR. sy > 1.0_DP) THEN
!!$           write(*,*) 'sy exceeded'
!!$           stop
!!$        END IF
!!$         IF (sz < 0.0_DP .OR. sz > 1.0_DP) THEN
!!$            write(*,*) 'sz exceeded'
!!$            write(*,*) sz
!!$            write(*,*) xlattice, ylattice, zlattice
!!$           stop
!!$        END IF

        ! Now generate 27 replicas of this atoms
        
        DO i = -1,1
           
           DO j = -1, 1
              
              DO k = -1, 1
                 
                 count = count + 1
                 
                 sx_this = sx + REAL(i,DP)
                 sy_this = sy + REAL(j,DP)
                 sz_this = sz + REAL(k,DP)
                 
                 ! now obtain the Cartesian coordinates of the super lattice atom
                 
                 CALL Fractional_To_Cartesian(sx_this, sy_this, sz_this, &
                      rx_this, ry_this, rz_this, this_box)
                 
                 x_super_lat(count) = rx_this
                 y_super_lat(count) = ry_this
                 z_super_lat(count) = rz_this
                 
                 ! assign the vDW type and charge to this atom
                 
                 IF (int_sim_type == sim_pot_map) THEN
                    
                    type_super_lat(count) = nonbond_list(iatom,zeo_id)%atom_type_number
                    charge_super_lat(count) = nonbond_list(iatom,zeo_id)%charge
                    
                 END IF
                 
                 WRITE(21,'(A,3x,4(F15.5,3X))') nonbond_list(iatom,zeo_id)%element, &
                      x_super_lat(count), y_super_lat(count), z_super_lat(count)
                 
              END DO ! k = -1, 1
              
           END DO ! j = -1, 1
           
        END DO ! i = -1, 1
        
     END DO ! iatom = 1, natoms(zeo_id)

  END DO ! imol = 1, nmols_to_read(zeo_id,this_box)

  CLOSE(UNIT=21)
  
  IF (n_superlat_atoms /= count) THEN

     err_msg = ''
     err_msg(1) = 'Total number of zeolite atoms in the supercell does not match'
     err_msg(2) = 'Number of atoms calculated from the replications '//Int_To_String(n_superlat_atoms)
     err_msg(3) = 'Number of atoms calculated while assinging coordinates '//Int_To_String(count)
     CALL Clean_Abort(err_msg,'Replicate_Unit_Cell')

  END IF


END SUBROUTINE Replicate_Unit_Cell
!--------------------------------------------------------------
SUBROUTINE Generate_Grid
  !----------------------------------------------------------------
  !
  !
  ! This subroutine generates the coordinates of the grids in
  ! the primary unit cell for a lattice potential map generation
  !
  ! First written by Jindal Shah on 01/09/13
  !
  !-----------------------------------------------------------------


  USE Global_Variables
  USE Potential_Map_Variables
  USE File_Names
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER :: i, nxtotal, nytotal, nztotal
  INTEGER :: ix, iy, iz, n_total_grids, count, ibox
  REAL(DP) :: x_this, y_this, z_this, sx, sy, sz

  ! Determine number of grids in each of the directions

  na_grid = NINT(box_list(1)%basis_length(1)/grid_spacing) + 1
  nb_grid = NINT(box_list(1)%basis_length(2)/grid_spacing) + 1
  nc_grid = NINT(box_list(1)%basis_length(3)/grid_spacing) + 1

  ! Recalculate the grid spacing in each of the directions in fractional coordinates

  a_spacing = 1.0_DP/REAL(na_grid,DP)
!  write(*,*) a_spacing, box_list(1)%basis_length, grid_spacing
  b_spacing = 1.0_DP/REAL(nb_grid,DP)
  c_spacing = 1.0_DP/REAL(nc_grid,DP)

  
!!$  IF (l_slit_pore(1)) THEN
!!$     nzgrid = INT(pore_width/grid_spacing) + 1
!!$     zstep = pore_width/REAL(nzgrid,DP)
!!$  ELSE
!!$     nzgrid = INT(box_list(1)%basis_length(3)/grid_spacing) + 1
!!$     zstep = box_list(1)%basis_length(3)/REAL(nzgrid,DP)
!!$  END IF


  WRITE(logunit,*)
  WRITE(logunit,'(A50,F10.8)') 'Actual grid spacing in a direction is', a_spacing * box_list(1)%basis_length(1)
  WRITE(logunit, '(A50,F10.8)') 'Actual grid spacing in b  direction is', b_spacing * box_list(1)%basis_length(2)
  WRITE(logunit,'(A50,F10.8)') 'Actual grid spacing in c direction is', c_spacing * box_list(1)%basis_length(3)

  CALL Allocate_Grid_Locations
 
  
  n_total_grids = na_grid * nb_grid * nc_grid
  ALLOCATE(xgrid(n_total_grids), ygrid(n_total_grids), zgrid(n_total_grids), Stat = AllocateStatus)
  IF (AllocateStatus /=0) THEN
     err_msg =''
     err_msg(1) = 'Memoery could not be allocated for storing grid coordinates'
     CALL Clean_Abort('Generate_Grid', err_msg)
  END IF

  ! find the grid coordinates
  count = 0
  ibox = 1
  DO ix = 1, na_grid
     DO iy = 1, nb_grid
        DO iz = 1,nc_grid

           count = count + 1
           sx = a_grid(ix)
           sy = b_grid(iy)
           sz = c_grid(iz)

           CALL Fractional_To_Cartesian(sx,sy,sz,x_this,y_this,z_this,ibox)
           xgrid(count) = x_this
           ygrid(count) = y_this
           zgrid(count) = z_this

        END DO
     END DO
  END DO
  
 ! write(*,*) 'grids'
!  write(*,*) a_grid(1), a_grid(na_grid), a_spacing
!  write(*,*) b_grid(1), b_grid(nb_grid), b_spacing
!  write(*,*) c_grid(1), c_grid(nc_grid), c_spacing
!  write(*,*) 'total number of grids', na_grid * nb_grid * nc_grid


  WRITE(logunit,*) 
  WRITE(logunit,*) 'Finished generation of the grid points'
  WRITE(logunit,*)
  WRITE(logunit,*) 'Number of grid points in a direction '//TRIM(Int_To_String(na_grid))
  WRITE(logunit,*)
  WRITE(logunit,*) 'Number of grid points in b direction '//TRIM(Int_To_String(nb_grid))
  WRITE(logunit,*)
  WRITE(logunit,*) 'Number of grid points in c direction '//TRIM(Int_To_String(nc_grid))
  WRITE(logunit,*)
  WRITE(logunit,*) 'Replication in the a, b and c directions achieved'

  ! At this point, read in the information for the supercell and expand the box
  ! First we store the unit cell information into a zeolite unit cell variable

!!$  zeo_unit_cell = box_list(1)
!!$
!!$  nxtotal = 2 * nx_zeo + 1
!!$  nytotal = 2 * ny_zeo + 1
!!$  IF (l_slit_pore(1)) THEN
!!$     nztotal = 1
!!$  ELSE
!!$     nztotal = 2 * nz_zeo + 1
!!$  END IF
!!$
!!$  box_list(1)%length(1,1) = nxtotal * box_list(1)%length(1,1)
!!$  box_list(1)%length(2,2) = nytotal * box_list(1)%length(2,2)
!!$  box_list(1)%length(3,3) = nztotal * box_list(1)%length(3,3)
!!$
!!$  CALL Compute_Cell_Dimensions(1)
!!$
!!$  WRITE(logunit,*) 'Replication in the x, y and z directions achieved'


END SUBROUTINE Generate_Grid
!--------------------------------------------------------------
SUBROUTINE Generate_Neighbor_List
  !--------------------------------------------------------------
  !
  ! This subroutine will divide the zeolite unit cell
  ! into small boxes of approximately 2.5 A cubic boxes. The center
  ! of each cubic box will be used to generate a neighbor list
  ! for a given box. The neighbor list will be useful for potential
  ! evaluation for grids that are found in a given box.
  !
  ! The output from this routine are the following variables
  !
  ! xbox, ybox, zbox = number of small boxes in the a, b, and c directions
  ! lxbox, lybox, lzbox = corresponding basis lengths
  !
  ! link_zeo(i,j,k) = .TRUE. or .FALSE.
  ! link_zeo_start(i,j,k) = n_start
  ! link_zeo_end(i,j,k) = n_end
  ! zeo_atom_int(n_start) = atom_id_start
  ! zeo_atom_int(n_end) = atom_id_end
  !
  !
  ! where n_start and n_end refer to the starting and ending index for a
  ! box identified by the (i,j,k). These indices are then used to located the
  ! atoms with which all the grid points located in (i,j,k) will interact with.
  !
  !-----------------------------------------------------------------

  USE Global_Variables
  USE Potential_Map_Variables
  IMPLICIT NONE

  INTEGER :: i, j, k, n_interaction
  INTEGER :: iatom, ibox

  REAL(DP) :: pad_distance_sq, pad_distance
  REAL(DP) :: sx, sy, sz, xc, yc, zc, cutoff_sq
  REAL(DP) :: dx, dy, dz, rsq

  LOGICAL :: l_first

  Type(Box_Class) :: small_box

  ! Degbugging variables
  INTEGER :: ii, this_atom

  ! Determine number of boxes in each direction

  xbox = NINT(box_list(1)%basis_length(1)/2.5_DP)
  ybox = NINT(box_list(1)%basis_length(2)/2.5_DP)
  zbox = NINT(box_list(1)%basis_length(3)/2.5_DP)

  ! length in each of the directions

  lxbox = box_list(1)%basis_length(1)/REAL(xbox,DP)
  lybox = box_list(1)%basis_length(2)/REAL(ybox,DP)
  lzbox = box_list(1)%basis_length(3)/REAL(zbox,DP)

  ! Allocate memory for interaction linked id

  ALLOCATE(link_zeo_start(xbox,ybox,zbox), Stat = AllocateStatus)

  IF (AllocateStatus /= 0 ) THEN

     err_msg = ''
     err_msg(1) = 'Linked start list cannont be assigned memory'
     err_msg(2) = 'The error occurred while generating neighbor list'
     err_msg(2) = 'Check potential map generation module'
     CALL Clean_Abort(err_msg,'Generate neighbor list')

  END IF

  ALLOCATE(link_zeo_end(xbox,ybox,zbox), Stat = AllocateStatus)

  IF (AllocateStatus /= 0 ) THEN

     err_msg = ''
     err_msg(1) = 'Linked end list cannont be assigned memory'
     err_msg(2) = 'The error occurred while generating neighbor list'
     err_msg(2) = 'Check potential map generation module'
     CALL Clean_Abort(err_msg,'Generate neighbor list')

  END IF

  ALLOCATE(link_zeo(xbox,ybox,zbox), Stat = AllocateStatus)

  IF (AllocateStatus /= 0 ) THEN

     err_msg = ''
     err_msg(1) = 'link_zeo array cannont be assigned memory'
     err_msg(2) = 'The error occurred while generating neighbor list'
     err_msg(2) = 'Check potential map generation module'
     CALL Clean_Abort(err_msg,'Generate neighbor list')

  END IF

  ALLOCATE(zeo_atom_int(xbox*ybox*zbox*n_superlat_atoms), Stat = AllocateStatus)

  IF (AllocateStatus /= 0 ) THEN

     err_msg = ''
     err_msg(1) = 'The array zeo_atom_int cannot be assigned a memory'
     err_msg(2) = 'The error occurred while generating neighbor list'
     err_msg(3) = 'Check potential map generation module'
     CALL Clean_Abort(err_msg, 'Generate_Neighbor_List')

  END IF

  ! We will use the center of these boxes for finding the neighbors.
  ! We need to take into account the fact that computing the neighbors
  ! just based on the center may not be enough for the atoms located 
  ! at the corners of the boxes. So we allow for an additional ``padding"
  ! to the cutoff.

  ! The padding distance is equal to the half of the body diagonal for
  ! a cell with side 2.5 A.

  small_box = box_list(1)

  ! get the cell matrix for this box
  ! the a-vector
  small_box%length(:,1) = (lxbox/box_list(1)%basis_length(1)) &
       * box_list(1)%length(:,1)
  ! the b-vector
  small_box%length(:,2) = (lybox/box_list(1)%basis_length(2)) * &
       box_list(1)%length(:,2)
  ! the c-vector
  small_box%length(:,3) = (lzbox/box_list(1)%basis_length(3)) * box_list(1)%length(:,3)

  ! compute the length of the body diagonal for small_box
  pad_distance_sq = (small_box%length(1,1) + small_box%length(1,2) + &
       small_box%length(1,3)) ** 2
  pad_distance_sq = pad_distance_sq + (small_box%length(2,1) + &
       small_box%length(2,2) + small_box%length(2,3)) ** 2
  pad_distance_sq = pad_distance_sq + (small_box%length(3,1) + &
       small_box%length(3,2) + small_box%length(3,3)) ** 2

  ! The padding is half the length of the body diagonal.
  pad_distance = 0.5_DP * SQRT(pad_distance_sq)

!  write(*,*) pad_distance, 'padding'
  ! now compute the cutoff to be used for constructing neighbor list
  cutoff_sq = (rcut_vdw(1)+pad_distance) ** 2

  ! loop over the cells now

  n_interaction = 0

  link_zeo_start(:,:,:) = 0
  link_zeo_end(:,:,:) = 0

  ! By default, energy is to be computed for all the grid points
  link_zeo(:,:,:) = .TRUE.

!!$  write(*,*) 'number of boxes'
!!$  write(*,*) 'x', xbox
!!$  write(*,*) 'y', ybox
!!$  write(*,*) 'z', zbox

  ibox = 1
  
  DO i = 1, xbox

     DO j = 1, ybox

        DO k = 1, zbox

           ! Find the center of this box in terms of fractional coordinates
          
           sx = (REAL(i,DP) - 0.5_DP) * lxbox/box_list(1)%basis_length(1) - 0.5_DP
           sy = (REAL(j,DP) - 0.5_DP) * lybox/box_list(1)%basis_length(2) - 0.5_DP
           sz = (REAL(k,DP) - 0.5_DP) * lzbox/box_list(1)%basis_length(3) - 0.5_DP

           ! obtain the Cartesian coordinates of the center

           CALL Fractional_To_Cartesian(sx,sy,sz,xc,yc,zc,ibox)
!           write(23,*)'c', xc, yc, zc
           ! Now loop over all the zeolite atoms (replicated super cell)

           l_first = .TRUE.

           DO iatom = 1, n_superlat_atoms

              dx = x_super_lat(iatom) - xc
              dy = y_super_lat(iatom) - yc
              dz = z_super_lat(iatom) - zc

              rsq = dx * dx + dy * dy + dz * dz

              IF ( rsq > cutoff_sq ) CYCLE

              IF (l_first) THEN
                ! this is the first atom for interaction 
                 link_zeo_start(i,j,k) = n_interaction + 1
                 l_first = .FALSE.
              END IF

              ! if lattice parameters are zero cycle

              ! increase the number of interaction

              n_interaction = n_interaction + 1

              zeo_atom_int(n_interaction) = iatom

           END DO

           ! set the end of the linked list. It is possible that
           ! that the box may not interact with any of the lattice
           ! atoms, for example in a slit pore geometry with large
           ! pore width, points in the center may not experience
           ! the presence of the wall atoms

           IF (link_zeo_start(i,j,k) == 0 ) THEN
              link_zeo_end(i,j,k) = 0
              link_zeo(i,j,k) = .FALSE.
           ELSE
              link_zeo_end(i,j,k) = n_interaction
           END IF

           ! check the zeolite atoms around this point

!          write(24,'(5(I8,2x))') i, j, k, link_zeo_end(i,j,k) - link_zeo_start(i,j,k)
!!$
!!$           WRITE(24,*) link_zeo_end(i,j,k) - link_zeo_start(i,j,k) + 2
!!$           WRITE(24,*)
!!$           WRITE(24,*) 'C', xc, yc, zc
!!$           DO ii = link_zeo_start(i,j,k), link_zeo_end(i,j,k)
!!$              ! identify the atom
!!$              this_atom = zeo_atom_int(ii)
!!$              WRITE(24,'(A,2X,3(F10.6,2X))') 'H',  x_super_lat(this_atom), &
!!$                   y_super_lat(this_atom), z_super_lat(this_atom)
!!$           END DO
!!$           stop

        END DO

     END DO

  END DO

  
!write(*,*) 'total number of interactions to be evaluated', n_interaction

END SUBROUTINE Generate_Neighbor_List

SUBROUTINE Allocate_Grid_Locations

  USE Global_Variables, ONLY: err_msg, AllocateStatus
  USE Potential_Map_Variables

  IMPLICIT NONE
  INTEGER :: i
  
  
 ! allocate memory for grid locations
  
  IF (.NOT. ALLOCATED(a_grid))  ALLOCATE(a_grid(na_grid), Stat = AllocateStatus)
  
  IF (AllocateStatus /= 0) THEN

     err_msg = ''
     err_msg(1) = 'Memory could not be assigned to the grid array in the a direction'
     err_msg(2) = 'Probably grid size is too small'
     err_msg(3) = 'Check grid spacing in the input file'
     CALL Clean_Abort(err_msg,'Generate_Grid')

  END IF

  IF (.NOT. ALLOCATED(b_grid))  ALLOCATE(b_grid(nb_grid), Stat = AllocateStatus)

  IF (AllocateStatus /= 0) THEN

     err_msg = ''
     err_msg(1) = 'Memory could not be assigned to the grid array in the b direction'
     err_msg(2) = 'Probably grid size is too small'
     err_msg(3) = 'Check grid spacing in the input file'
     CALL Clean_Abort(err_msg,'Generate_Grid')

  END IF


  IF (.NOT. ALLOCATED(c_grid)) ALLOCATE(c_grid(nc_grid), Stat = AllocateStatus)

  IF (AllocateStatus /= 0) THEN

     err_msg = ''
     err_msg(1) = 'Memory could not be assigned to the grid array in the c direction'
     err_msg(2) = 'Probably grid size is too small'
     err_msg(3) = 'Check grid spacing in the input file'
     CALL Clean_Abort(err_msg,'Generate_Grid')

  END IF


  ! the grid coordinates start at (-0.5,-0.5,-.5) as the original box has been moved
  ! such that the fractional coordinates in any direction are between -0.5 and 0.5

  a_grid(1) = -0.5_DP
  b_grid(1) = -0.5_DP
  c_grid(1) = -0.5_DP

  ! x_grid goes from zero to alpha - xstep, the point at alpha is the
  ! same as the one at zero. Note that the grids are in fractional coordinates

  DO i = 2, na_grid

     a_grid(i) = a_grid(i-1) + a_spacing

  END DO
! write(*,*) a_grid
  DO i = 2, nb_grid 

     b_grid(i) = b_grid(i-1) + b_spacing

  END DO

  
  DO i = 2, nc_grid 

     c_grid(i) = c_grid(i-1) + c_spacing
 !    write(*,*) i, c_grid(i)

  END DO

END SUBROUTINE Allocate_Grid_Locations
