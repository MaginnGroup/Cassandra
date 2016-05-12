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

    IF (block_average) THEN
       WRITE(this_unit,'(A)') '# Block averages'
    ELSE
       WRITE(this_unit,'(A)') '# Instantaneous properties'
    END IF

    write_str = ""
    IF (sim_length_units == 'Steps') THEN
      write_str = "# MC_STEP"
    ELSE IF (sim_length_units == 'Sweeps') THEN
      write_str = "# MC_SWEEP"
    END IF
    
    WRITE(this_unit,'(A12,2X)',ADVANCE='NO') ADJUSTL(write_str)
    
    ! Now write strings for the rest of fields.
    
    DO ii = 1, prop_per_file(file_number,this_box)-1
       
       WRITE(this_unit,'(A16,2X)',ADVANCE='NO') (TRIM(prop_output(ii,file_number,this_box)))
       
    END DO

    WRITE(this_unit,'(A16,2X)') (TRIM(prop_output(ii,file_number,this_box)))


    ALLOCATE(prop_unit(prop_per_file(file_number,this_box)+1))

    prop_unit(:) = ""
    prop_unit(1) ='# '

    WRITE(this_unit, '(A12,2X)',ADVANCE='NO') ADJUSTL(prop_unit(1))

    DO ii = 1, prop_per_file(file_number,this_box) - 1

       prop_to_write = prop_output(ii,file_number,this_box)

       IF (prop_to_write(1:6) == 'Energy') THEN

          prop_unit(ii) = '(kJ/mol)-Ext '

       ELSE IF (prop_to_write == 'Pressure') THEN
          
          prop_unit(ii) = '(bar)'

       ELSE IF (prop_to_write == 'Volume') THEN

          prop_unit(ii) = '(A^3)'

       ELSE IF (prop_to_write == 'Density') THEN

          prop_unit(ii) = '(molec/A^3)'

       END IF

       WRITE(this_unit,'(A16,2X)',ADVANCE='NO') (TRIM(prop_unit(ii)))

    END DO

    prop_to_write = prop_output(ii,file_number,this_box)
    
    IF (prop_to_write(1:6) == 'Energy') THEN
       
       prop_unit(ii) = '(kJ/mol)-Ext '
       
    ELSE IF (prop_to_write == 'Pressure') THEN
       
       prop_unit(ii) = '(bar)'
       
    ELSE IF (prop_to_write == 'Volume') THEN
       
       prop_unit(ii) = '(A^3)'
       
    ELSE IF (prop_to_write == 'Density') THEN
       
       prop_unit(ii) = '(molec/A^3)'

    END IF
    
    WRITE(this_unit,'(A16,2X)') (TRIM(prop_unit(ii)))

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
   INTEGER :: nmols_is, nmols_box, iis
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
            write_buff(ii+1) = ac_energy(this_box)%total / REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%total
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol


      ELSE IF (prop_written == 'Temperature') THEN

         write_buff(ii+1) = 1.0_DP/(beta(this_box)*kboltz)

      ELSE IF (prop_written == 'Thermodynamic_Pressure') THEN

         write_buff(ii+1) = pressure(this_box)*atomic_to_bar

      ELSE IF (prop_written == 'Energy_LJ') THEN

         IF ( block_average) THEN
            write_buff(ii+1) = (ac_energy(this_box)%inter_vdw + ac_energy(this_box)%intra_vdw) / &
                 REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%inter_vdw + energy(this_box)%intra_vdw
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_Elec') THEN

         IF ( block_average) THEN
            write_buff(ii+1) = (ac_energy(this_box)%inter_q + ac_energy(this_box)%intra_q) / &
                 REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%inter_q &
                             + energy(this_box)%intra_q &
                             + energy(this_box)%ewald_reciprocal &
                             + energy(this_box)%ewald_self
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Energy_Intra') THEN

         IF (block_average) THEN
            write_buff(ii+1) = (ac_energy(this_box)%intra)/REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%intra
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Pressure') THEN

         CALL Compute_System_Total_Force(this_box)

         Pressure_tensor(:,:,this_box) = W_tensor_total(:,:,this_box) &
                                       / box_list(this_box)%volume
         P_inst(this_box) = ((Pressure_tensor(1,1,this_box) &
                            + Pressure_tensor(2,2,this_box) &
                            + Pressure_tensor(3,3,this_box)) / 3.0_DP) &
                          * atomic_to_bar
    
         IF(int_vdw_sum_style(this_box) == vdw_cut_tail) THEN
            P_inst(this_box) = P_inst(this_box) &
                             + virial(this_box)%lrc &
                             / box_list(this_box)%volume * atomic_to_bar
         END IF

         nmols_box = SUM(nmols(:,1:nbr_boxes))
         P_ideal(this_box) = nmols_box &
                           / box_list(this_box)%volume * temperature(this_box) &
                           * p_const

         write_buff(ii+1) = P_ideal(this_box) + P_inst(this_box)

      ELSE IF (prop_written == 'Volume') THEN
         
         IF (block_average) THEN
            write_buff(ii+1) = (ac_volume(this_box))/REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = box_list(this_box)%volume
         END IF

      ELSE IF (prop_written == 'Enthalpy') THEN
         IF (block_average) THEN
            write_buff(ii+1) = (ac_enthalpy(this_box))/REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = energy(this_box)%total + pressure(this_box) * box_list(this_box)%volume
         END IF
         write_buff(ii+1) = write_buff(ii+1) * atomic_to_kJmol

      ELSE IF (prop_written == 'Nmols') THEN

         IF (block_average) THEN
               write_buff(ii+1) = ac_nmols(is,this_box) / REAL(nthermo_freq,DP)
         ELSE
               write_buff(ii+1) = nmols(is,this_box)
         END IF

         ! increment the species index by 1 so that if there is
         ! another species and if nmols is to be output for that
         ! species, we will have correct index
         is = is + 1
         

      ELSE IF (prop_written == 'Density') THEN

         IF (block_average) THEN
            write_buff(ii+1) = ac_density(is_dens,this_box)/ REAL(nthermo_freq,DP)
         ELSE
            write_buff(ii+1) = REAL(nmols(is_dens,this_box),DP) / box_list(this_box)%volume
         END IF
         ! increment the species index by 1 for the same reason as species
         ! in 'Nmols' was incremented
         is_dens = is_dens + 1

      ELSE IF (prop_written == 'Chemical_Potential') THEN
         write_buff(ii+1) = chpot(is_cp,this_box) / REAL(ntrials(is_cp,this_box)%cpcalc)
         is_cp = is_cp + 1
         
      END IF
      
      ! At the end increment property counter by 1

      ii = ii + 1

   END DO

   ! write the line buffer to the property file

   IF (sim_length_units == 'Steps') THEN
     WRITE(this_unit,'(I12,2X)',ADVANCE='NO') i_mcstep
   ELSE IF (sim_length_units == 'Sweeps') THEN
     WRITE(this_unit,'(I12,2X)',ADVANCE='NO') i_mcstep / steps_per_sweep
   END IF

   DO ii = 1, prop_per_file(file_number,this_box)-1

      WRITE(this_unit,'(E16.8,2X)',ADVANCE='NO') write_buff(ii+1)

   END DO
   WRITE(this_unit,'(E16.8,2X)') write_buff(prop_per_file(file_number,this_box)+1)
   
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
  WRITE(M_XYZ_unit,*)
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

