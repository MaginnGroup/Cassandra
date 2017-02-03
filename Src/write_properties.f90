!*******************************************************************************
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
SUBROUTINE Write_Properties(this_box)
  ! The subroutine will write desired properties to the property files. It is
  ! called by respective drivers such as.
  !
  ! CALLED BY
  !
  !        gcmc_driver
  !        gemc_driver
  !        nptmc_driver
  !        nvtmc_driver
  ! CALLS
  !
  !   None
  !
  ! 08/12/13 : Created beta version
!*******************************************************************************

  USE Global_Variables
  USE File_Names
  USE Energy_Routines

  IMPLICIT NONE

  CHARACTER(24) :: write_str
  INTEGER :: i, this_box, this_unit, is, n_start, n_end
  INTEGER :: ifrac, my_ifrac


  DO i = 1, nbr_prop_files(this_box)
     
     !-- check to see if the file is open or not
     this_unit = (this_box-1)*MAXVAL(nbr_prop_files) + i + propunit
     
     IF (first_open(i,this_box)) THEN
        
        OPEN(unit=this_unit,file=prop_files(i,this_box))
        
        ! write the header information that indicates the properties contained
        ! in the file
        CALL Write_Header(i)

        first_open(i,this_box) = .FALSE.

     END IF
        
     CALL Write_Properties_Buffer(i)
     
  END DO


CONTAINS

  SUBROUTINE Write_Header(file_number)

    IMPLICIT NONE

    INTEGER :: file_number, ii
    CHARACTER*120 :: prop_to_write
    CHARACTER*120, ALLOCATABLE :: prop_unit(:)

    ! Header line 1
    IF (block_average) THEN
       WRITE(this_unit,'(A)') '# Block averages'
    ELSE
       WRITE(this_unit,'(A)') '# Instantaneous properties'
    END IF

    ! Header line 2, property names
    write_str = ""
    IF (sim_length_units == 'Steps' .OR. sim_length_units == 'Minutes') THEN
      write_str = "# MC_STEP"
    ELSE IF (sim_length_units == 'Sweeps') THEN
      write_str = "# MC_SWEEP"
    END IF
    
    WRITE(this_unit,'(A12,2X)',ADVANCE='NO') ADJUSTL(write_str)
    DO ii = 1, prop_per_file(file_number,this_box)
       WRITE(this_unit,'(A16,2X)',ADVANCE='NO') (TRIM(prop_output(ii,file_number,this_box)))
    END DO
    WRITE(this_unit,*)

    ! Header line 3, units
    ALLOCATE(prop_unit(prop_per_file(file_number,this_box)+1))

    prop_unit(:) = ""
    prop_unit(1) ='# '

    WRITE(this_unit, '(A12,2X)',ADVANCE='NO') ADJUSTL(prop_unit(1))

    DO ii = 1, prop_per_file(file_number,this_box)

       prop_to_write = prop_output(ii,file_number,this_box)

       IF (prop_to_write(1:6) == 'Energy') THEN

          prop_unit(ii) = '(kJ/mol)-Ext'

       ELSE IF (prop_to_write == 'Pressure') THEN
          
          prop_unit(ii) = '(bar)'

       ELSE IF (prop_to_write == 'Volume') THEN

          prop_unit(ii) = '(A^3)'

       ELSE IF (prop_to_write(1:7) == 'Density') THEN

          prop_unit(ii) = '(molec/A^3)'

       ELSE if (prop_to_write == 'Mass_Density') THEN
 
          prop_unit(ii) = '(kg/m^3)'

       ELSE IF (prop_to_write(1:18) == 'Chemical_Potential') THEN

          prop_unit(ii) = '(kJ/mol)'

       END IF

       WRITE(this_unit,'(A16,2X)',ADVANCE='NO') (TRIM(prop_unit(ii)))

    END DO

    WRITE(this_unit,*)

    DEALLOCATE(prop_unit)
    
  END SUBROUTINE Write_Header

  SUBROUTINE Write_Properties_Buffer(file_number)
   !************************************************************************
   ! The subroutine fills in a line buffer based on which properties are to
   ! be written and then write the buffer to a file.
   !
   ! Writte by Jindal Shah
   !
   ! Revision History
   !
   !*************************************************************************

   USE Simulation_Properties
   
   INTEGER :: file_number, ii, is, is_dens, is_cp, is_lambda, total_frac
   REAL(DP) :: mass_density
   REAL(DP),DIMENSION(:), ALLOCATABLE :: write_buff
   CHARACTER(FILENAME_LEN) :: prop_written

   REAL(DP) :: temp_pres
   
   ALLOCATE(write_buff(prop_per_file(file_number,this_box)+1))

   !***********************************************************************
   ! Fill the elements of write_buff with each of the properties
   ! Note that these are average properties over the frequency interval
   !***********************************************************************

!   DO ii = 1, prop_per_file(file_number,this_box)

   ii = 1
   is = 1
   is_dens = 1
   is_cp = 1
   is_lambda = 1

   DO WHILE ( ii <= prop_per_file(file_number,this_box))

      prop_written = prop_output(ii,file_number,this_box)

      IF (prop_written == 'Energy_Total') THEN

         IF ( block_average) THEN
            ac_energy(this_box,iblock)%total = ac_energy(this_box,iblock)%total / data_points_per_block
            write_buff(ii+1) = ac_energy(this_box,iblock)%total
         ELSE
            write_buff(ii+1) = energy(this_box)%total
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol


      ELSE IF (prop_written == 'Temperature') THEN

         write_buff(ii+1) = 1.0_DP/(beta(this_box)*kboltz)

      ELSE IF (prop_written == 'Pressure_Setpoint') THEN

         write_buff(ii+1) = pressure(this_box)%setpoint * atomic_to_bar

      ELSE IF (prop_written == 'Energy_LJ') THEN

         IF ( block_average) THEN
            ac_energy(this_box,iblock)%inter_vdw = ac_energy(this_box,iblock)%inter_vdw / data_points_per_block
            ac_energy(this_box,iblock)%intra_vdw = ac_energy(this_box,iblock)%intra_vdw / data_points_per_block
            write_buff(ii+1) = ac_energy(this_box,iblock)%inter_vdw + ac_energy(this_box,iblock)%intra_vdw
         ELSE
            write_buff(ii+1) = energy(this_box)%inter_vdw + energy(this_box)%intra_vdw
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_Elec') THEN

         IF ( block_average) THEN
            ac_energy(this_box,iblock)%inter_q = ac_energy(this_box,iblock)%inter_q / data_points_per_block
            ac_energy(this_box,iblock)%intra_q = ac_energy(this_box,iblock)%intra_q / data_points_per_block
            write_buff(ii+1) = ac_energy(this_box,iblock)%inter_q + ac_energy(this_box,iblock)%intra_q
         ELSE
            write_buff(ii+1) = energy(this_box)%inter_q &
                             + energy(this_box)%intra_q &
                             + energy(this_box)%ewald_reciprocal &
                             + energy(this_box)%self
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_Intra') THEN

         IF (block_average) THEN
            ac_energy(this_box,iblock)%intra = ac_energy(this_box,iblock)%intra / data_points_per_block
            write_buff(ii+1) = ac_energy(this_box,iblock)%intra
         ELSE
            write_buff(ii+1) = energy(this_box)%intra
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Pressure') THEN

         IF (block_average) THEN
            ac_pressure(this_box,iblock) = ac_pressure(this_box,iblock) / data_points_per_block
            write_buff(ii+1) = ac_pressure(this_box,iblock)
         ELSE
            IF (pressure(this_box)%last_calc /= i_mcstep) THEN
               pressure(this_box)%last_calc = i_mcstep
               CALL Compute_Pressure(this_box)
            END IF

            write_buff(ii+1) = pressure(this_box)%computed
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_bar

      ELSE IF (prop_written == 'Volume') THEN
         
         IF (block_average) THEN
            ac_volume(this_box,iblock) = ac_volume(this_box,iblock) / data_points_per_block
            write_buff(ii+1) = ac_volume(this_box,iblock)
         ELSE
            write_buff(ii+1) = box_list(this_box)%volume
         END IF

      ELSE IF (prop_written == 'Enthalpy') THEN

         IF (block_average) THEN
            ac_enthalpy(this_box,iblock) = ac_enthalpy(this_box,iblock) / data_points_per_block
            write_buff(ii+1) = ac_enthalpy(this_box,iblock)
         ELSE
            IF (int_sim_type == sim_npt .OR. int_sim_type == sim_gemc_npt) THEN
               write_buff(ii+1) = energy(this_box)%total + pressure(this_box)%setpoint * box_list(this_box)%volume
            ELSE
               IF (pressure(this_box)%last_calc /= i_mcstep) THEN
                  pressure(this_box)%last_calc = i_mcstep
                  CALL Compute_Pressure(this_box)
               END IF

               write_buff(ii+1) = energy(this_box)%total + pressure(this_box)%computed * box_list(this_box)%volume
            END IF
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written(1:5) == 'Nmols') THEN

         IF (block_average) THEN
            ac_nmols(is,this_box,iblock) = ac_nmols(is,this_box,iblock) / data_points_per_block
            write_buff(ii+1) = ac_nmols(is,this_box,iblock)
         ELSE
            write_buff(ii+1) = nmols(is,this_box)
         END IF

         ! increment the species index by 1 so that if there is
         ! another species and if nmols is to be output for that
         ! species, we will have correct index
         is = is + 1
         

      ELSE IF (prop_written(1:7) == 'Density') THEN

         IF (block_average) THEN
            ac_density(is_dens,this_box,iblock) = ac_density(is_dens,this_box,iblock) / data_points_per_block
            write_buff(ii+1) = ac_density(is_dens,this_box,iblock)
         ELSE
            write_buff(ii+1) = REAL(nmols(is_dens,this_box),DP) / box_list(this_box)%volume
         END IF
         ! increment the species index by 1 for the same reason as species
         ! in 'Nmols' was incremented
         is_dens = is_dens + 1

      ELSE IF (prop_written(1:18) == 'Chemical_Potential') THEN
         write_buff(ii+1) = chpot(is_cp,this_box) / REAL(ntrials(is_cp,this_box)%cpcalc)
         is_cp = is_cp + 1
         
      ELSE IF (prop_written == "Mass_Density") THEN
         IF (block_average) THEN
            ac_mass_density(this_box,iblock) = ac_mass_density(this_box,iblock) / data_points_per_block
            write_buff(ii+1) = ac_mass_density(this_box,iblock)
         ELSE
            mass_density = 0.0_DP
            DO is = 1, nspecies
               mass_density = mass_density + REAL(nmols(is,this_box),DP) * species_list(is)%molecular_weight
            END DO
            write_buff(ii+1) = mass_density / box_list(this_box)%volume
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kgm3
      END IF
      
      ! At the end increment property counter by 1

      ii = ii + 1

   END DO

   ! write the line buffer to the property file

   IF (sim_length_units == 'Steps' .OR. sim_length_units == 'Minutes') THEN
     WRITE(this_unit,'(I12,2X)',ADVANCE='NO') i_mcstep
   ELSE IF (sim_length_units == 'Sweeps') THEN
     WRITE(this_unit,'(I12,2X)',ADVANCE='NO') i_mcstep / steps_per_sweep
   END IF

   DO ii = 1, prop_per_file(file_number,this_box)

      WRITE(this_unit,'(E16.8,2X)',ADVANCE='NO') write_buff(ii+1)

   END DO
   WRITE(this_unit,*)
   
   DEALLOCATE(write_buff)

 END SUBROUTINE Write_Properties_Buffer
 
END SUBROUTINE Write_Properties
!*******************************************************************************
SUBROUTINE Write_Coords(this_box)
  !*****************************************************************************
  ! The subroutine writes coordinates of simulation box for later analyis of
  ! RDFs. It gets called by driver routines.
  !
  ! CALLED BY
  !
  !        gcmc_driver
  !        gemc_driver
  !        nptmc_driver
  !        nvtmc_driver
  !
  ! 08/12/13 (JS) : Created beta version
  !*****************************************************************************
  
  USE Global_Variables
  USE Simulation_Properties
  USE File_Names, ONLY : movie_header_unit,movie_xyz_unit 

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: this_box

 !******************************************************************************
  
  INTEGER :: ii, jj, is, nmolecules_is, im, this_im, ia 
  INTEGER :: M_XYZ_unit,MH_unit,Num_Atoms

  ! Write the information about volume

  MH_unit = movie_header_unit + this_box
  M_XYZ_unit = movie_xyz_unit + this_box

  Num_Atoms = 0

  WRITE(MH_unit,*) box_list(this_box)%volume
  ! The cell matrix
  DO ii = 1, 3
     WRITE(MH_unit,*)(box_list(this_box)%length(ii,jj), jj=1,3)
  END DO

  WRITE(MH_unit,*)
  WRITE(MH_unit,*) nspecies
  
  !-- Number of molecules of each of the species
  DO is = 1, nspecies
     WRITE(MH_unit,*) is,nmols(is,this_box)
     Num_Atoms = Num_Atoms + nmols(is,this_box)*natoms(is)
  END DO

  !-- Write the coordinates of molecules in this box
  WRITE(M_XYZ_unit,*) Num_Atoms
  WRITE(M_XYZ_unit,*) 'MC_STEP: ', i_mcstep
  DO is = 1,nspecies
     DO im = 1, nmols(is,this_box)
        this_im = locate(im,is,this_box)
        IF(molecule_list(this_im,is)%live) THEN
           DO ia = 1, natoms(is)
              WRITE(M_XYZ_unit,*) nonbond_list(ia,is)%element, &
                   atom_list(ia,this_im,is)%rxp, &
                   atom_list(ia,this_im,is)%ryp, &
                   atom_list(ia,this_im,is)%rzp
           END DO
        END IF
     END DO
  END DO

END SUBROUTINE Write_Coords

SUBROUTINE Write_Mean_Error(ibox)
  ! The subroutine will write desired properties to the property files. It is
  ! called by respective drivers such as.
  !
  ! CALLED BY
  !
  !        gcmc_driver
  !        gemc_driver
  !        nptmc_driver
  !        nvtmc_driver
  ! CALLS
  !
  !   None
  !
  ! 06/22/16 : Created beta version
!*******************************************************************************

   USE Global_Variables
   USE File_Names
   USE Energy_Routines
   USE Simulation_Properties

   IMPLICIT NONE

   CHARACTER(24) :: write_str
   INTEGER :: i, ibox, this_unit

   INTEGER :: ii, is, is_dens, is_cp
   REAL(DP),DIMENSION(:), ALLOCATABLE :: write_mean
   REAL(DP),DIMENSION(:), ALLOCATABLE :: write_err
   CHARACTER(FILENAME_LEN) :: prop_written

   !***********************************************************************
   ! Fill the elements of write_buff with each of the properties
   ! Note that these are average properties over the frequency interval
   !***********************************************************************

   DO i = 1, nbr_prop_files(ibox)
     
      this_unit = (ibox-1)*MAXVAL(nbr_prop_files) + i + propunit
      ALLOCATE(write_mean(prop_per_file(i,ibox)+1))
      ALLOCATE(write_err(prop_per_file(i,ibox)+1))
      write_mean = 0.0_DP
      write_err  = 0.0_DP

      ii = 1
      is = 1
      is_dens = 1
      is_cp = 1

      DO WHILE ( ii <= prop_per_file(i,ibox))

         prop_written = prop_output(ii,i,ibox)

         IF (prop_written == 'Energy_Total') THEN

            write_mean(ii+1) = SUM(ac_energy(ibox,:)%total) / nbr_blocks
            DO iblock = 1, nbr_blocks
               write_err(ii+1)  = write_err(ii+1)  + (ac_energy(ibox,iblock)%total - write_mean(ii+1))**2
            END DO
            write_mean(ii+1) = write_mean(ii+1) * atomic_to_kJmol
            write_err(ii+1)  = sqrt(write_err(ii+1) / nbr_blocks)  * atomic_to_kJmol

         ELSE IF (prop_written == 'Temperature') THEN

            write_mean(ii+1) = 1.0_DP/(beta(ibox)*kboltz)
            write_err(ii+1)  = 0.0_DP

         ELSE IF (prop_written == 'Pressure_Setpoint') THEN

            write_mean(ii+1) = pressure(ibox)%setpoint * atomic_to_bar
            write_err(ii+1)  = 0.0_DP

         ELSE IF (prop_written == 'Energy_LJ') THEN

            write_mean(ii+1) = (SUM(ac_energy(ibox,:)%inter_vdw) + SUM(ac_energy(ibox,:)%intra_vdw)) / nbr_blocks
            DO iblock = 1, nbr_blocks
                write_err(ii+1)  = write_err(ii+1)  + (ac_energy(ibox,iblock)%inter_vdw + ac_energy(ibox,iblock)%intra_vdw &
                                 - write_mean(ii+1))**2
            END DO
            write_mean(ii+1) = write_mean(ii+1) * atomic_to_kJmol
            write_err(ii+1)  = sqrt(write_err(ii+1) / nbr_blocks) * atomic_to_kJmol

         ELSE IF (prop_written == 'Energy_Elec') THEN

            write_mean(ii+1) = (SUM(ac_energy(ibox,:)%inter_q) + SUM(ac_energy(ibox,:)%intra_q)) / nbr_blocks
            DO iblock = 1, nbr_blocks
               write_err(ii+1) = write_err(ii+1) + (ac_energy(ibox,iblock)%inter_q + ac_energy(ibox,iblock)%intra_q &
                               - write_mean(ii+1))**2
            END DO
            write_mean(ii+1) = write_mean(ii+1) * atomic_to_kJmol
            write_err(ii+1)  = sqrt(write_err(ii+1) / nbr_blocks) * atomic_to_kJmol

         ELSE IF (prop_written == 'Energy_Intra') THEN

            write_mean(ii+1) = SUM(ac_energy(ibox,:)%intra) / nbr_blocks
            DO iblock = 1, nbr_blocks
               write_err(ii+1) = write_err(ii+1) + (ac_energy(ibox,iblock)%intra - write_mean(ii+1))
            END DO
            write_mean(ii+1) = write_mean(ii+1) * atomic_to_kJmol
            write_err(ii+1)  = sqrt(write_err(ii+1) / nbr_blocks) * atomic_to_kJmol

         ELSE IF (prop_written == 'Pressure') THEN

            write_mean(ii+1) = SUM(ac_pressure(ibox,:)) / nbr_blocks
            DO iblock = 1, nbr_blocks
               write_err(ii+1) = write_err(ii+1) + (ac_pressure(ibox,iblock) - write_mean(ii+1))**2
            END DO
            write_mean(ii+1) = write_mean(ii+1) * atomic_to_bar
            write_err(ii+1)  = sqrt(write_err(ii+1) / nbr_blocks) * atomic_to_bar

         ELSE IF (prop_written == 'Volume') THEN
            
            write_mean(ii+1) = SUM(ac_volume(ibox,:)) / nbr_blocks
            DO iblock = 1, nbr_blocks
               write_err(ii+1) = write_err(ii+1) + (ac_volume(ibox,iblock) - write_mean(ii+1))**2
            END DO
            write_err(ii+1)  = sqrt(write_err(ii+1) / nbr_blocks)

         ELSE IF (prop_written == 'Enthalpy') THEN

            write_mean(ii+1) = SUM(ac_enthalpy(ibox,:)) / nbr_blocks
            DO iblock = 1, nbr_blocks
               write_err(ii+1) = write_err(ii+1) + (ac_enthalpy(ibox,iblock) - write_mean(ii+1))**2
            END DO
            write_mean(ii+1) = write_mean(ii+1) * atomic_to_kJmol
            write_err(ii+1)  = sqrt(write_err(ii+1) / nbr_blocks) * atomic_to_kJmol

         ELSE IF (prop_written(1:5) == 'Nmols') THEN

            write_mean(ii+1) = SUM(ac_nmols(is,ibox,:)) / nbr_blocks
            DO iblock = 1, nbr_blocks
               write_err(ii+1) = write_err(ii+1) + (ac_nmols(is,ibox,iblock) - write_mean(ii+1))**2
            END DO
            write_err(ii+1)  = sqrt(write_err(ii+1) / nbr_blocks)

            ! increment the species index by 1 so that if there is
            ! another species and if nmols is to be output for that
            ! species, we will have correct index
            is = is + 1
            
         ELSE IF (prop_written(1:7) == 'Density') THEN

            write_mean(ii+1) = SUM(ac_density(is_dens,ibox,:)) / nbr_blocks
            DO iblock = 1, nbr_blocks
               write_err(ii+1) = write_err(ii+1) + (ac_density(is_dens,ibox,iblock) - write_mean(ii+1))**2
            END DO
            write_err(ii+1)  = sqrt(write_err(ii+1) / nbr_blocks)

            ! increment the species index by 1 for the same reason as species
            ! in 'Nmols' was incremented
            is_dens = is_dens + 1

! chem_pot needs to be tested
!         ELSE IF (prop_written(1:18) == 'Chemical_Potential') THEN
!
!            write_mean(ii+1) = chpot(is_cp,ibox) / REAL(ntrials(is_cp,ibox)%cpcalc)
!            is_cp = is_cp + 1
!            
         ELSE IF (prop_written == "Mass_Density") THEN

            write_mean(ii+1) = SUM(ac_mass_density(ibox,:)) / nbr_blocks
            DO iblock = 1, nbr_blocks
               write_err(ii+1) = write_err(ii+1) + (ac_mass_density(ibox,iblock) - write_mean(ii+1))**2
            END DO
            write_mean(ii+1) = write_mean(ii+1) * atomic_to_kgm3
            write_err(ii+1)  = sqrt(write_err(ii+1) / nbr_blocks) * atomic_to_kgm3

         END IF
         
         ! At the end increment property counter by 1

         ii = ii + 1

      END DO

      ! write the line buffer to the property file
      WRITE(this_unit,'(A)',ADVANCE='NO') '#-------------'
      DO ii = 1, prop_per_file(i,ibox)
         WRITE(this_unit,'(A)',ADVANCE='NO') '------------------'
      END DO
      WRITE(this_unit,*)

      WRITE(this_unit,'(A12,2X)',ADVANCE='NO') 'mean'
      DO ii = 1, prop_per_file(i,ibox)
         WRITE(this_unit,'(E16.8,2X)',ADVANCE='NO') write_mean(ii+1)
      END DO
      WRITE(this_unit,*)
      
      WRITE(this_unit,'(A12,2X)',ADVANCE='NO') 'stdev'
      DO ii = 1, prop_per_file(i,ibox)
         WRITE(this_unit,'(E16.8,2X)',ADVANCE='NO') write_err(ii+1)
      END DO
      WRITE(this_unit,*)
      
      DEALLOCATE(write_mean)
      DEALLOCATE(write_err)

   END DO

END SUBROUTINE Write_Mean_Error
