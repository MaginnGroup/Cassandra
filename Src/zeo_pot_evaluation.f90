!********************************************************************************
! CASSANDRA - Computational Atomistic Simulation Software at Notre Dame 
! for Research in Academia.
! http://molsim.wiki.zoho.com/
! Copyright (2007) University of Notre Dame.
! Authors: Ed Maginn (ed@nd.edu); Jindal Shah (jshah@nd.edu)
!********************************************************************************
SUBROUTINE Zeolite_Elec_Potential_Evaluation
  !-------------------------------------------------------------------------------
  !
  ! This subroutine evaluates the electrostatic potential energy at each of the
  ! super lattice points in the main zeolite simulation cell. This potential map
  ! will be used in conjunction with the vDW map generated for each of the atoomtypes
  ! to arrive at the final compression for the electrostatic map and that vDW map.
  !
  ! 
  !
  !-------------------------------------------------------------------------------

  USE Global_Variables
  USE Potential_Map_Variables
  USE File_Names
  IMPLICIT NONE

  INTEGER :: ix, iy, iz, ix_box, iy_box, iz_box, n_interaction, lat_atom, iatom, itype
  INTEGER :: count
  REAL(DP) :: this_alpha, this_alpha2, this_alpha4, dx, dy, dz, rsq
  REAL(DP) :: rij, r2i, r4i, r6i, factor, factor1, factor2, factor3, q_pot

  INTEGER :: ibox
  REAL(DP) :: frac_x, frac_y, frac_z, this_x, this_y, this_z

  ! Allocate all the arrays related to energy and its derivatives

  IF (.NOT. ALLOCATED(grid_potq)) ALLOCATE(grid_potq(na_grid,nb_grid,nc_grid))
  IF (.NOT. ALLOCATED(grid_potqx)) ALLOCATE(grid_potqx(na_grid,nb_grid,nc_grid))
  IF (.NOT. ALLOCATED(grid_potqy)) ALLOCATE(grid_potqy(na_grid,nb_grid,nc_grid))
  IF (.NOT. ALLOCATED(grid_potqz)) ALLOCATE(grid_potqz(na_grid,nb_grid,nc_grid))
  IF (.NOT. ALLOCATED(grid_potqxy)) ALLOCATE(grid_potqxy(na_grid,nb_grid,nc_grid))
  IF (.NOT. ALLOCATED(grid_potqyz)) ALLOCATE(grid_potqyz(na_grid,nb_grid,nc_grid))
  IF (.NOT. ALLOCATED(grid_potqxz)) ALLOCATE(grid_potqxz(na_grid,nb_grid,nc_grid))
  IF (.NOT. ALLOCATED(grid_potqxyz)) ALLOCATE(grid_potqxyz(na_grid,nb_grid,nc_grid))

  grid_potq = 0.0_DP ; grid_potqx = 0.0_DP ; grid_potqy = 0.0_DP ; grid_potqz = 0.0_DP
  grid_potqxy = 0.0_DP ; grid_potqyz = 0.0_DP ; grid_potqxz = 0.0_DP ; grid_potqxyz = 0.0_DP

  ! loop over all the grid point
  
  this_alpha = alpha_ewald(1)
  this_alpha2 = alpha_ewald(1) * alpha_ewald(1)
  this_alpha4 = this_alpha2 * this_alpha2
  count = 0

  ibox = 1
  gridx: DO ix = 1, na_grid

     ix_box = MIN(INT(a_grid(ix)/lxbox) + 1, xbox)

     gridy: DO iy = 1, nb_grid

        iy_box = MIN(INT(b_grid(iy)/lybox) + 1, ybox)

        gridz: DO iz = 1, nc_grid

           iz_box = MIN(INT(c_grid(iz)/lzbox) + 1, zbox)

           ! check to see if this grid point is
           IF ( .NOT. link_zeo(ix_box,iy_box,iz_box)) CYCLE
           ! Note that the energy and its first, second and third
           ! derivatives are zero at this grid point by virtue
           ! of the initialization.

           ! if here compute the interactions

           interaction:DO n_interaction = link_zeo_start(ix_box,iy_box,iz_box), &
                link_zeo_end(ix_box,iy_box,iz_box)

              lat_atom = zeo_atom_int(n_interaction)

              frac_x = a_grid(ix)
              frac_y = b_grid(iy)
              frac_z = c_grid(iz)

              CALL Fractional_To_Cartesian(frac_x, frac_y, frac_z, &
                   this_x, this_y, this_z, ibox)
              
              dx = this_x - x_super_lat(lat_atom)
              dy = this_y - y_super_lat(lat_atom) 
              dz = this_z - z_super_lat(lat_atom) 
              
              rsq = dx * dx + dy * dy + dz * dz
              
              IF ( rsq > rcut_coulsq(1) .OR. rsq < rcut_lowsq) CYCLE

              ! else compute the energy due to this lattice atom

              ! compute some factors first

              rij = DSQRT(rsq)
              r2i = 1.0_DP/rsq
              r4i = r2i * r2i
              r6i = r4i * r2i

              q_pot = erfc(this_alpha*rij)/rij
              
              factor = this_alpha * DEXP(-this_alpha2 * rsq)/DSQRT(pi)
              
              factor1 = -(2.0_DP * factor + q_pot) * r2i
              
              factor2 = 4.0_DP * this_alpha2 * factor * r2i + &
                        6.0_DP * factor * r4i + 3.0_DP * q_pot * r4i
              
              factor3 = 8.0_DP * this_alpha4 * factor * r2i + 20.0_DP * this_alpha2 * factor * r4i + &
                        30.0_DP * factor * r6i + 15.0_DP * q_pot * r6i
              factor3 = -1.0_DP * factor3

              grid_potq(ix,iy,iz) = grid_potq(ix,iy,iz) + q_pot * charge_super_lat(lat_atom) 
              grid_potqx(ix,iy,iz) = grid_potqx(ix,iy,iz) + factor1 * dx * charge_super_lat(lat_atom)
              grid_potqy(ix,iy,iz) = grid_potqy(ix,iy,iz) + factor1 * dy * charge_super_lat(lat_atom)
              grid_potqz(ix,iy,iz) = grid_potqz(ix,iy,iz) + factor1 * dz * charge_super_lat(lat_atom)
              grid_potqxy(ix,iy,iz) = grid_potqxy(ix,iy,iz) + factor2 * dx * dy * charge_super_lat(lat_atom)
              grid_potqyz(ix,iy,iz) = grid_potqyz(ix,iy,iz) + factor2 * dy * dz * charge_super_lat(lat_atom)
              grid_potqxz(ix,iy,iz) = grid_potqxz(ix,iy,iz) + factor2 * dz * dx * charge_super_lat(lat_atom)
              grid_potqxyz(ix,iy,iz) = grid_potqxyz(ix,iy,iz) + factor3 * dx * dy * dz * charge_super_lat(lat_atom)

           END DO interaction

        END DO gridz

     END DO gridy
     
  END DO gridx

  ! At the end of the routine, we will have computed the electric field at each of the grid
  ! points


  ! Now multiply the potential energy by the lowest_charge on the sorbate and the conversion factor

  grid_pot = grid_potq * charge_factor
  grid_potqx = grid_potqx  * charge_factor
  grid_potqy = grid_potqy  * charge_factor
  grid_potqz = grid_potqz  * charge_factor
  grid_potqxy = grid_potqxy  * charge_factor
  grid_potqyz = grid_potqyz  * charge_factor
  grid_potqxz = grid_potqxz * charge_factor
  grid_potqxyz = grid_potqxyz * charge_factor

!!$  IF (ALLOCATED(map_output_file)) DEALLOCATE(map_output_file)
!!$  ALLOCATE(map_output_file(1))
!!$  map_output_file(1) = 'charge.map'
!!$  itype = 1
!!$  CALL Compress_Zeo_Potential_Map(itype)
!!$
!!$  ! At the end of the file add the information about the charges for the
!!$  ! n_superlat_atoms, charge for these supercell lattice atoms
!!$
!!$  OPEN(UNIT=map_output_unit+itype,file='charge.map',FORM='UNFORMATTED', POSITION='APPEND')
!!$
!!$  WRITE(map_output_unit+itype) n_superlat_atoms
!!$  WRITE(map_output_unit+itype) charge_super_lat
!!$  
!!$  CLOSE(UNIT=map_output_unit)

END SUBROUTINE Zeolite_Elec_Potential_Evaluation


SUBROUTINE Zeolite_Potential_Evaluation
  !-----------------------------------------------------------------------------------
  !
  ! zeo_pot_evaluation.f90
  !
  ! This subroutine calculates the potential, the first, second and third derivative
  ! at each of the grid points of the central lattice unit cell.
  !
  ! First written by Jindal Shah on 01/10/13
  !
  !
  ! 01/14/13 (JS) : Added two more routines: Compress_Zeo_Potential_Map and Compact_Zeo_Potential_Map
  !                 committed to repository
  !
  ! 03/04/13 (JS) : Fixed a few bugs. Added number of grid points in each direction and total number
  !                 of grid points in each direction to be output in the map file. Grid_Pointer is
  !                 explicitly written out.
  !
  ! 04/03/13 (JS) : fixed how x, y and z components (dx, dy, dz) are calculated for the vector connecting
  !                 a grid point with a given lattice point.
  !
  ! 07/02/13 (JS) : The potential is at each point is calculated based on
  !                 link_zeo_start and link_zeo_end arrays.
  !
  ! 07/30/13 (JS) : If a given grid site was determined not to interact with any of the
  !                 lattice atoms, then the potential, the first derivative and the second
  !                 derivatives are set to zero. The subroutine otherwise would crash
  !                 before the 'link_zeo' variable was defined.
  !
  ! 03/13/14 (JS) : Added the subroutine to calcuate electrostatic interactions at grid points
  !                 in the central unit cell.
  !------------------------------------------------------------------------------------

  ! Identify atomtypes for which the potential map needs to be generated.
  ! Exclude atomtypes that exclusively belong to the lattice atoms

  ! The flow of the subroutine is as follows.
  ! Step 1: Identify unique atomtypes for which a potential map is to be generated.
  ! Step 2: Obtaine the repulsive and attractive part of the LJ potential for all the
  !         the pairs
  ! Step 3: Obtain the atomtypes for the zeolite
  ! Step 4: Compute and store energy and forces for each of the grid points in the
  !         arrays pot, pot_x, pot_y, pot_z, pot_xy, pot_xz, pot_yz, and pot_xyz
  !
  

  USE Global_Variables
  USE Potential_Map_Variables
  USE File_Names
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER :: ispecies, iatom, itype, nzeo_pot_atomtypes, jtype
  INTEGER :: iatoms, this_type, la_atom_type, ilat_type, lat_atom, lat_atom_type
  INTEGER :: ix, iy, iz, ix_box, iy_box, iz_box, n_interaction
  INTEGER :: this_lat_type, count

!  INTEGER, ALLOCATABLE :: nlat_atom_types(:)

  REAL(DP) :: ptemp, px_temp, py_temp, pz_temp, pxy_temp, pyz_temp, pxz_temp
  REAL(DP) :: pxyz_temp, dx, dy, dz, x_this, y_this, z_this, rsq
  REAL(DP) :: r2i, r6i, r8i, r10i, r12i, r14i, r16i, r18i

!  REAL(DP), ALLOCATABLE, DIMENSION(:,:) :: zeo_LJ12_table, zeo_LJ6_table
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: sum12, sum6, sum14x, sum14y, sum14z
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: sum8x, sum8y, sum8z
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: sum16xy, sum16yz, sum16xz, sum10xy, sum10yz, sum10xz
  REAL(DP), ALLOCATABLE, DIMENSION(:) :: sum18xyz, sum12xyz
  
  REAL(DP) :: eps, sig, sig2, sig6
  REAL(DP) :: in_time, out_time


  LOGICAL, ALLOCATABLE :: l_compute(:), l_found_type(:)
  LOGICAL :: lexist, high_energy

  INTEGER :: ibox
  REAL(DP) :: frac_x, frac_y, frac_z

  Allocate(l_compute(nbr_atomtypes))
  ALLOCATE(map_output_file(nbr_atomtypes))
  map_output_file(:) = ''

  l_compute(:) = .FALSE.
  nzeo_pot_atomtypes = 0

  ! zeolite species will also be set the first species
  DO ispecies = 2, nspecies
  
     DO iatom = 1, natoms(ispecies)
        
        itype = nonbond_list(iatom,ispecies)%atom_type_number
        
        IF (.NOT. l_compute(itype) ) THEN
           
           nzeo_pot_atomtypes = nzeo_pot_atomtypes + 1
           l_compute(itype) = .TRUE.
           map_output_file(itype) = 'sorbate_'//TRIM(nonbond_list(iatom,ispecies)%atom_name)//'.map'

           INQUIRE(file=map_output_file(itype), exist = lexist)

           IF (lexist) THEN
              ! This file already exists so abort
              err_msg = ""
              err_msg(1) = 'The output map file'
              err_msg(2) = map_output_file(itype)
              err_msg(3) = "already exists. Cannot overwrite"
              CALL Clean_Abort(err_msg,"zeo_pot_evalulation.f90")
              
           END IF

        END IF
        
     END DO

  END DO

!  WRITE(*,*) 'total number of atomtypes for computing map ', nzeo_pot_atomtypes
!  WRITE(*,*) l_compute
!!$  stop

  ! Pull out cross interaction paramters for zeolite and all other
  ! species

  ALLOCATE(zeo_LJ12_table(nbr_atomtypes,nbr_atomtypes))
  ALLOCATE(zeo_LJ6_table(nbr_atomtypes,nbr_atomtypes))

  DO itype = 1, nbr_atomtypes

     DO jtype = 1, nbr_atomtypes

        eps = vdw_param1_table(itype,jtype)

        sig = vdw_param2_table(itype,jtype)

        sig2 = sig * sig
        sig6 = sig2 * sig2 * sig2

        zeo_LJ6_table(itype,jtype) = 4.0_DP * eps * sig6
        zeo_LJ12_table(itype,jtype) = zeo_LJ6_table(itype,jtype) * sig6

     END DO
  END DO


  ! determine total number of lattice types and corresponding atom types

  nlat_types = 0
  ALLOCATE(nlat_atom_types(nbr_atomtypes))
  nlat_atom_types(:) = 0

  ! zeolite species is always the first species
  ispecies = 1

  
  ALLOCATE(l_found_type(nbr_atomtypes))
  
  l_found_type(:) = .FALSE.
  
  DO iatoms = 1, natoms(ispecies)
     
     this_type = nonbond_list(iatoms,ispecies)%atom_type_number

     IF (l_found_type(this_type)) CYCLE
     
     nlat_types = nlat_types + 1
     
     nlat_atom_types(nlat_types) = this_type
     
     l_found_type(this_type) = .TRUE.

 !    write(*,*) this_type
  END DO


  

  ALLOCATE(sum12(nlat_types), sum6(nlat_types))
  ALLOCATE(sum14x(nlat_types), sum14y(nlat_types),sum14z(nlat_types))
  ALLOCATE(sum8x(nlat_types), sum8y(nlat_types), sum8z(nlat_types))
  ALLOCATE(sum16xy(nlat_types), sum16yz(nlat_types), sum16xz(nlat_types))
  ALLOCATE(sum10xy(nlat_types), sum10yz(nlat_types), sum10xz(nlat_types))
  ALLOCATE(sum18xyz(nlat_types), sum12xyz(nlat_types))
  
  ! Here we have a list of atomtypes for which the map needs to be computed

  ALLOCATE(grid_pot(na_grid,nb_grid,nc_grid), grid_potx(na_grid,nb_grid,nc_grid))
  ALLOCATE(grid_poty(na_grid,nb_grid,nc_grid), grid_potz(na_grid,nb_grid,nc_grid))
  ALLOCATE(grid_potxy(na_grid,nb_grid,nc_grid), grid_potxz(na_grid,nb_grid,nc_grid))
  ALLOCATE(grid_potyz(na_grid,nb_grid,nc_grid), grid_potxyz(na_grid,nb_grid,nc_grid))

  ibox = 1
!  nb_grid = 100
!  nc_grid = 100
!  na_grid = 100
  DO itype = 1, nbr_atomtypes

     ! loop over all the grid points
     IF ( .NOT. l_compute(itype)) CYCLE

     grid_pot(:,:,:) = 0.0_DP
     grid_potx(:,:,:) = 0.0_DP
     grid_poty(:,:,:) = 0.0_DP
     grid_potz(:,:,:) = 0.0_DP
     grid_potxy(:,:,:) = 0.0_DP
     grid_potyz(:,:,:) = 0.0_DP
     grid_potxz(:,:,:) = 0.0_DP
     grid_potxyz(:,:,:) = 0.0_DP

     CALL CPU_TIME(in_time)
     count = 0
     gridx:DO ix = 1, na_grid

        ! Locate the box the x-coordinate is in
        ix_box = MIN(INT((a_grid(ix)+0.5_DP)*xbox)+1,xbox)
        
        gridy:DO iy = 1, nb_grid

           iy_box = MIN(INT((b_grid(iy)+0.5_DP)*ybox)+1,ybox)

           gridz: DO iz = 1, nc_grid

              iz_box = MIN(INT((c_grid(iz)+0.5_DP)*zbox)+1,zbox)

!              write(1007,*) ix, iy, iz, ix_box, iy_box, iz_box
              ! zero the temporary storage units

              sum12(:) = 0.0_DP
              sum6(:) = 0.0_DP
              sum14x(:) = 0.0_DP
              sum14y(:) = 0.0_DP
              sum14z(:) = 0.0_DP
              sum8x(:) = 0.0_DP
              sum8y(:) = 0.0_DP
              sum8z(:) = 0.0_DP
              sum16xy(:) = 0.0_DP
              sum16yz(:) = 0.0_DP
              sum16xz(:) = 0.0_DP
              sum10xy(:) = 0.0_DP
              sum10yz(:) = 0.0_DP
              sum10xz(:) = 0.0_DP
              sum18xyz(:) = 0.0_DP
              sum12xyz(:) = 0.0_DP

              ptemp = 0.0_DP
              
              ! first derivative

              px_temp = 0.0_DP
              py_temp = 0.0_DP
              pz_temp = 0.0_DP

              ! second derivative

              pxy_temp = 0.0_DP
              pyz_temp = 0.0_DP
              pxz_temp = 0.0_DP

              ! third derivative

              pxyz_temp = 0.0_DP

              ! find the x, y, and z coordinates of this grid location
              frac_x = a_grid(ix)
              frac_y = b_grid(iy)
              frac_z = c_grid(iz)

              CALL Fractional_To_Cartesian(frac_x, frac_y, frac_z, &
                   x_this, y_this, z_this, ibox)

              IF (.NOT. link_zeo(ix_box,iy_box,iz_box)) THEN
                 
                 grid_pot(ix,iy,iz) = 0.0_DP
                 
                 grid_potx(ix,iy,iz) = 0.0_DP
                 grid_poty(ix,iy,iz) = 0.0_DP
                 grid_potz(ix,iy,iz) = 0.0_DP
                 
                 grid_potxy(ix,iy,iz) = 0.0_DP
                 grid_potyz(ix,iy,iz) = 0.0_DP
                 grid_potxz(ix,iy,iz) = 0.0_DP
                 
                 grid_potxyz(ix,iy,iz) = 0.0_DP

                 ! go to the next grid point
                 CYCLE

              END IF

              
              ! Identify total number of atoms over which interaction is to
              ! be computed
              high_energy = .false.
              DO n_interaction = link_zeo_start(ix_box,iy_box,iz_box), &
                                 link_zeo_end(ix_box,iy_box,iz_box)

                 lat_atom = zeo_atom_int(n_interaction)
                 
                 dx = x_this - x_super_lat(lat_atom)
                 dy = y_this - y_super_lat(lat_atom) 
                 dz = z_this - z_super_lat(lat_atom)

                 rsq = dx * dx + dy * dy + dz * dz

                

                 IF ( rsq > rcut_vdwsq(1)) CYCLE
                 IF ( rsq < rcut_lowsq) THEN
                    grid_pot(ix,iy,iz) = 1000.0_DP * max_kBT
                    high_energy = .true.
                    EXIT
                 END IF

                 ! if here, then the interaction needs to be evaluated
                 ! identify the lattice atom type of lat_atom

                 lat_atom_type = type_super_lat(lat_atom)


                 ! obtain inverse powers of distances

                 r2i = 1.0_DP/rsq
                 r6i = r2i * r2i * r2i
                 r8i = r6i * r2i
                 r10i = r8i * r2i
                 r12i = r10i * r2i
                 r14i = r12i * r2i
                 r16i = r14i * r2i
                 r18i = r16i * r2i
                
                 ! update the summations

                 ! potentials

                 sum12(lat_atom_type) = sum12(lat_atom_type) + r12i
                 sum6(lat_atom_type) = sum6(lat_atom_type) + r6i

                 ! x, y and z derivaties

                 sum14x(lat_atom_type) = sum14x(lat_atom_type) + dx * r14i
                 sum8x(lat_atom_type) = sum8x(lat_atom_type) + dx * r8i

                 sum14y(lat_atom_type) = sum14y(lat_atom_type) + dy * r14i
                 sum8y(lat_atom_type) = sum8y(lat_atom_type) +   dy * r8i

                 sum14z(lat_atom_type) = sum14z(lat_atom_type) + dz * r14i
                 sum8z(lat_atom_type) = sum8z(lat_atom_type) +  dz * r8i

                 ! double derivates

                 sum16xy(lat_atom_type) = sum16xy(lat_atom_type) + dx*dy*r16i
                 sum10xy(lat_atom_type) = sum10xy(lat_atom_type) + dx*dy*r10i

                 sum16yz(lat_atom_type) = sum16yz(lat_atom_type) + dy*dz*r16i
                 sum10yz(lat_atom_type) = sum10yz(lat_atom_type) + dy*dz*r10i

                 sum16xz(lat_atom_type) = sum16xz(lat_atom_type) + dx*dz*r16i
                 sum10xz(lat_atom_type) = sum10xz(lat_atom_type) + dx*dz*r10i

                 ! third derivative

                 sum18xyz(lat_atom_type) = sum18xyz(lat_atom_type) + dx*dy*dz*r18i
                 sum12xyz(lat_atom_type) = sum12xyz(lat_atom_type) + dx*dy*dz*r12i

              END DO

              ! this concludes the calculation of the grid point located
              ! at ix, iy and iz

              ! Multiply the above summations with appropriate prefactors
              ! basically cross interaction parameters

              ! skip the grid point if it's a high energy point

              IF (.NOT. high_energy) THEN
              
                 DO ilat_type = 1, nlat_types
                    
                    this_lat_type = nlat_atom_types(ilat_type)

                    ptemp = ptemp + &
                         zeo_LJ12_table(this_lat_type,itype) * sum12(ilat_type) - &
                         zeo_LJ6_table(this_lat_type,itype) * sum6(ilat_type)
                    
                    px_temp = px_temp + &
                         (-12.0_DP * zeo_LJ12_table(this_lat_type,itype) * sum14x(ilat_type)) + &
                         (  6.0_DP * zeo_LJ6_table(this_lat_type,itype) * sum8x(ilat_type))
                    
                    py_temp = py_temp + &
                         (-12.0_DP * zeo_LJ12_table(this_lat_type,itype) * sum14y(ilat_type)) + &
                         (  6.0_DP * zeo_LJ6_table(this_lat_type,itype) * sum8y(ilat_type))
                    
                    pz_temp = pz_temp + &
                         (-12.0_DP * zeo_LJ12_table(this_lat_type,itype) * sum14z(ilat_type)) + &
                         (  6.0_DP * zeo_LJ6_table(this_lat_type,itype) * sum8z(ilat_type))
                    
                    pxy_temp = pxy_temp + &
                         ( 168.0_DP * zeo_LJ12_table(this_lat_type,itype) * sum16xy(ilat_type)) + &
                         ( -48.0_DP * zeo_LJ6_table(this_lat_type,itype) * sum10xy(ilat_type))
                    
                    pyz_temp = pyz_temp + &
                         ( 168.0_DP * zeo_LJ12_table(this_lat_type,itype) * sum16yz(ilat_type)) + &
                         ( -48.0_DP * zeo_LJ6_table(this_lat_type,itype) * sum10yz(ilat_type))

                    pxz_temp = pxz_temp + &
                         ( 168.0_DP * zeo_LJ12_table(this_lat_type,itype) * sum16xz(ilat_type)) + &
                         ( -48.0_DP * zeo_LJ6_table(this_lat_type,itype) * sum10xz(ilat_type))
                    
                    pxyz_temp = pxyz_temp + &
                         (-2688.0_DP * zeo_LJ12_table(this_lat_type,itype) * sum18xyz(ilat_type)) + &
                         (  480.0_DP * zeo_LJ6_table(this_lat_type,itype) * sum12xyz(ilat_type))
                    
                 END DO
                 
               
                 grid_pot(ix,iy,iz) = ptemp
                 
                 grid_potx(ix,iy,iz) = px_temp
                 grid_poty(ix,iy,iz) = py_temp
                 grid_potz(ix,iy,iz) = pz_temp
                 
                 grid_potxy(ix,iy,iz) = pxy_temp
                 grid_potyz(ix,iy,iz) = pyz_temp
                 grid_potxz(ix,iy,iz) = pxz_temp
                 
                 grid_potxyz(ix,iy,iz) = pxyz_temp

!                 write(1001,'(3(F10.5,2X),F16.8)') x_this, y_this, z_this,  grid_pot(ix,iy,iz)
!                 write(1002,*) ix, iy,iz,  grid_pot(ix,iy,iz)

              END IF
            END DO gridz
           
        END DO gridy
        
     END DO gridx
     
     CALL CPU_TIME(out_time)

     ! At this point we have arrays holding the values of potential and derivatives
     ! at the grid points. Not all the grid points are energetically favorable.
     ! So we can reduce the size of these arrays by applying a criterion that
     !
!     ALLOCATE(grid_q_pointer(na_grid,nb_grid,nc_grid))
!     grid_q_count = 0
     CALL Compress_Zeo_Potential_Map(itype)
!
  END DO

  WRITE(logunit,*) 'Time to compute, compress and write thew grid for the atom type '//Int_To_String(itype)
  WRITE(logunit,*) out_time - in_time
  stop
END SUBROUTINE Zeolite_Potential_Evaluation
!---------------------------------------------------------------------------------


!---------------------------------------------------------------------------------
SUBROUTINE Compress_Zeo_Potential_Map(itype)
  !-------------------------------------------------------------------------------
  !
  ! This subroutine calls compresses the potential and derivative related
  ! arrays so that only grid points with energies below a threshold
  ! are saved.
  !
  ! First written by Jindal Shah on 01/14/13
  !
  !-------------------------------------------------------------------------------

  USE Global_Variables
  USE Potential_Map_Variables
  IMPLICIT NONE

  INTEGER :: ix, iy, iz, count, count2, itype
  
  INTEGER, ALLOCATABLE :: grid_pointer(:,:,:)

  LOGICAL :: first_time, last_time

  ! Call the routine to determine which grids to be saved

  ! Allocate an array that will indicate if a given grid is to be retained.
  
  ALLOCATE(grid_pointer(na_grid,nb_grid,nc_grid))


  first_time = .TRUE.
  last_time = .FALSE.
  grid_pointer = 0
!  grid_q_pointer = 0

  CALL Compact_Zeo_Potential_Map(itype,grid_pot,grid_pointer,first_time,last_time)

  first_time = .FALSE.

  ! write potential to the disc

  CALL Compact_Zeo_Potential_Map(itype,grid_pot,grid_pointer,first_time,last_time)
  
  
  ! write 1st derivative to the disc

  CALL Compact_Zeo_Potential_Map(itype,grid_potx,grid_pointer,first_time,last_time)

  CALL Compact_Zeo_Potential_Map(itype,grid_poty,grid_pointer,first_time,last_time)

  CALL Compact_Zeo_Potential_Map(itype,grid_potz,grid_pointer,first_time,last_time)
  
  ! write 2nd derivative to the disc
  
  CALL Compact_Zeo_Potential_Map(itype,grid_potxy,grid_pointer,first_time,last_time)

  CALL Compact_Zeo_Potential_Map(itype,grid_potyz,grid_pointer,first_time,last_time)

  CALL Compact_Zeo_Potential_Map(itype,grid_potxz,grid_pointer,first_time,last_time)

  ! write the 3rd derivative and the grid_pointer to the disc

  last_time = .TRUE.

  CALL Compact_Zeo_Potential_Map(itype,grid_potxyz,grid_pointer,first_time,last_time)

  DEALLOCATE(grid_pointer)


END SUBROUTINE Compress_Zeo_Potential_Map
!---------------------------------------------------------------------------------

SUBROUTINE Compact_Zeo_Potential_Map(itype,input_array,grid_pointer,first_time,last_time)
  !------------------------------------------------------------------------------
  !
  ! This subroutine will determine the energy cutoff to be used for a given
  ! potential type and then will compress the zeolite potential map
  !
  ! First written by Jindal Shah on 01/11/12
  !
  !------------------------------------------------------------------------------

  USE Global_Variables
  USE Potential_Map_Variables
  USE File_Names
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER :: count, ix, iy, iz, count2, ii, jj, kk, iii, jjj, kkk, itype
  REAL(DP), INTENT(IN) :: input_array(na_grid, nb_grid, nc_grid)
  INTEGER, INTENT(INOUT) :: grid_pointer(na_grid, nb_grid, nc_grid)
  INTEGER, ALLOCATABLE :: temp_pointer(:,:,:)
  INTEGER :: ispecies, iatoms
 
  REAL(DP) :: emax
  REAL(DP), ALLOCATABLE :: map_to_disk(:)
  REAL(DP) :: lowest_charge, highest_charge, atom_charge, this_energy
  
  LOGICAL :: first_time, last_time
  IF(first_time) THEN
     
     ! Loop over all the grid points and mark the grids that are to be preserved
     !
  
     ! we will preserve all the energies such that beta * E <= max_kBT
     emax = max_kBT / beta(1)  ! note that the zeolite/solid/slit pore is always the first box 
     
     ! for evaluating the electrostatic potential, find out the lowest negative charge
     ! and the highest positive charge for atoms belonging to itype

     lowest_charge = 10.0_DP
     highest_charge = -10.0_DP
     
     DO ispecies = 2, nspecies
        DO iatoms = 1, natoms(ispecies)
           

           IF (itype == nonbond_list(iatoms,ispecies)%atom_type_number) THEN
              ! find the charge and compare with the lowest and highest
              atom_charge = nonbond_list(iatoms,ispecies)%charge
              lowest_charge = min(lowest_charge, atom_charge)
              highest_charge = max(highest_charge,atom_charge)
           END IF

        END DO
     END DO

     ! Report to the log file
     WRITE(logunit,*)
     WRITE(logunit,'(A35,2X,I2,A5,2X,F11.7)')' The highest charge for atom type', itype, ' is', highest_charge
     WRITE(logunit,'(A35,2X,I2,A5,2X,F11.7)') ' The lowest charge for the atom type', itype, ' is', lowest_charge

     IF (lowest_charge > highest_charge) THEN
        err_msg = ''
        err_msg(1) = 'The lowest charge was fount to be greater than the highest charge'
        err_msg(2) = 'Aborting...'
        CALL Clean_Abort(err_msg,'Compact_Zeo_Potential_Map')
     END IF
    
     ALLOCATE(temp_pointer(na_grid,nb_grid,nc_grid))
     temp_pointer = 0
     count = 0
     DO ix = 1, na_grid
        
        DO iy = 1, nb_grid
           
           DO iz = 1, nc_grid

              ! Determine the sign of the grid_potq(ix,iy,iz)
!              IF ( grid_potq(ix,iy,iz) < 0.0_DP) this_energy = grid_pot(ix,iy,iz) + grid_potq(ix,iy,iz) * highest_charge
!              IF ( grid_potq(ix,iy,iz) > 0.0_DP) this_energy = grid_pot(ix,iy,iz) + grid_potq(ix,iy,iz) * lowest_charge

              this_energy = grid_pot(ix,iy,iz)
              IF (this_energy < emax ) THEN
                 
                 ! keep this grid point
                 
                 count = count + 1
                 grid_q_count = grid_q_count + 1
                 grid_pointer(ix,iy,iz) = count
!                 grid_q_pointer(ix,iy,iz) = grid_q_count
                 temp_pointer(ix,iy,iz) = count

              END IF
              
           END DO
           
        END DO
        
     END DO
     
     WRITE(logunit,*) 
     WRITE(logunit,*) 'First step in compressing the map'
     WRITE(logunit,*) 'Determining the grid points to be preserved'
     WRITE(logunit,*) 'First pass.....'
     WRITE(logunit,*) 'Number of grid points to be stored '//Int_To_String(count)
     
     ! Now go through once again and save all the neighbors of the grids found in
     ! the previous step
     
     count2 = 0
     
     DO ix = 1, na_grid
        DO iy = 1, nb_grid
           DO iz = 1, nc_grid
              
              IF (temp_pointer(ix,iy,iz) == 0) CYCLE
              
              ! check neighboring grids
              
              DO ii = -1, 1
                 
                 DO jj = -1, 1
                    
                    DO kk = -1, 1
                       
                       iii = ii + ix
                       
                       IF ( iii == 0 ) THEN
                          iii = na_grid
                       ELSE IF (iii > na_grid) THEN
                          iii = 1
                       END IF
                       
                       ! similarly do for y and z directions
                       
                       jjj = jj + iy
                       
                       IF ( jjj == 0 ) THEN
                          jjj = nb_grid
                       ELSE IF (jjj > nb_grid) THEN
                          jjj = 1
                       END IF
                       
                       kkk = kk + iz
                       
                       IF ( kkk == 0 ) THEN
                          kkk = nc_grid
                       ELSE IF (kkk > nc_grid) THEN
                          kkk = 1
                       END IF
                       
                       IF (grid_pointer(iii,jjj,kkk) /= 0) CYCLE
                       ! if here then,  this grid point is not in the list, so add
                       
                       count = count + 1
                       grid_q_count = grid_q_count + 1
                       grid_pointer(iii,jjj,kkk) = count 

                       count2 = count2 + 1

                    END DO
                    
                 END DO
                 
              END DO
              
           END DO
           
        END DO

     END DO

     
     WRITE(logunit,*)
     WRITE(logunit,*) 'Completed second pass'
     WRITE(logunit,*) 'Number of points found in the second pass '//Int_To_String(count2)
     WRITE(logunit,*) 'Total number of grids to be saved '//Int_To_String(count)
     WRITE(logunit,*) 'Fraction of the grids saved ', REAL(count,DP)/REAL(na_grid * nb_grid * nc_grid, DP)
     ! write the number of grid points in each direction and total number of grids
     ! in the map file
     OPEN(unit=map_output_unit+itype,file=map_output_file(itype),FORM='UNFORMATTED')
     WRITE(map_output_unit+itype) na_grid, nb_grid, nc_grid, count

     DEALLOCATE(temp_pointer)


  ELSE

     ! grid_pointer array has already been obtained, use this
     ! array to write the arrays to the disk
     
     
     ! if it's second time then just write out the array to a file



     ALLOCATE(map_to_disk(MAXVAL(grid_pointer)))
     map_to_disk = 0.0_DP
     
     DO ix = 1, na_grid
        
        DO iy = 1, nb_grid
           
           DO iz = 1, nc_grid
              
              IF (grid_pointer(ix,iy,iz) == 0) CYCLE
              
              map_to_disk(grid_pointer(ix,iy,iz)) = input_array(ix,iy,iz)
              WRITE(map_output_unit+itype)grid_pointer(ix,iy,iz), map_to_disk(grid_pointer(ix,iy,iz))
              
              
           END DO
           
        END DO
        
     END DO
     
     WRITE(map_output_unit+itype)
  

     IF (last_time) THEN

        DO ix = 1, na_grid

           DO iy = 1, nb_grid

              DO iz = 1, nc_grid

                 
                 WRITE(map_output_unit+itype) ix, iy, iz, grid_pointer(ix,iy,iz)

              END DO

           END DO

        END DO
        CLOSE(map_output_unit+itype)

     END IF
     DEALLOCATE(map_to_disk)

     
     
  END IF

END SUBROUTINE Compact_Zeo_Potential_Map


