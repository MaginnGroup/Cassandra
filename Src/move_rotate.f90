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

!********************************************************************************
SUBROUTINE Rotate
!********************************************************************************

!********************************************************************************
! The subroutine performs rotation of a randomly chosen molecule about one
! of the axes. The energy of the molecule in old and new configurations will
! be determined and the Metropolis criterion will be used for accepting or
! rejecting the move.
!
! Called by
!
!   NVTMC_Driver
!   NPTMC_Driver
!
! Calls
!
!   Get_Nmolecules_Species
!   Get_Index_Molecule
!   Save_Old_Cartesian_Coordinates
!   Compute_Molecule_Nonbond_Inter_Energy
!   Rotate_Molecule_Axis
!   Update_System_Ewald_Reciprocal_Energy
!   Get_COM
!   Revert_Old_Cartesian_Coordinates
!
! Revision History
!
!   12/10/13 : Beta Release 
!********************************************************************************

  USE Type_Definitions
  USE Global_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE File_Names
  USE Pair_Nrg_Routines
  USE IO_Utilities

  IMPLICIT NONE

!  !$ include 'omp_lib.h'

  ! Local delcarations
  INTEGER  :: ibox   ! box index
  INTEGER  :: is         ! species index
  INTEGER  :: im, lm     ! molecule indices
  INTEGER  :: nmols_tot,mcstep ! number of molecules in the system

  REAL(DP) :: nmols_box(nbr_boxes)
  REAL(DP), ALLOCATABLE :: x_box(:), x_species(:)
  REAL(DP) :: randno
  REAL(DP), DIMENSION(:), ALLOCATABLE :: dx, dy, dz
  REAL(DP) :: delta_e, ln_pacc, success_ratio
  REAL(DP) :: E_vdw, E_qq, E_vdw_move, E_qq_move, E_reciprocal_move

  LOGICAL :: inter_overlap, overlap, accept_or_reject

  ! Pair_Energy arrays and Ewald implementation
  INTEGER :: position
  REAL(DP), ALLOCATABLE :: cos_mol_old(:), sin_mol_old(:)

  ! Framework energy related variables
  REAL(DP) :: E_framework, E_framework_move, E_correction_move
  LOGICAL :: framework_overlap

! Done with that section

  E_vdw_move = 0.0_DP
  E_qq_move = 0.0_DP
  E_vdw = 0.0_DP
  E_qq = 0.0_DP
  E_reciprocal_move = 0.0_DP
  inter_overlap = .FALSE.
  accept = .FALSE.

  ! Sum the total number of molecules 
  nmols_tot = 0 ! sum over species, box
  DO ibox = 1, nbr_boxes
    nmols_box(ibox) = 0
    DO is = 1, nspecies
      ! Only count mobile species
      IF ( max_rot(is,ibox) > 0. ) THEN
        nmols_tot = nmols_tot + nmols(is,ibox)
        nmols_box(ibox) = nmols_box(ibox) + nmols(is,ibox)
      END IF
    END DO
  END DO

  ! If there are no molecules then return
  IF (nmols_tot == 0) THEN
     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,5X,X,3X,X,I3,X,L8,X,9X,X,A9)') &
             i_mcstep, 'rotate' , ibox, accept,'no mols'
     END IF
     RETURN
  END IF

  ! If needed, choose a box based on its total mol fraction
  IF(nbr_boxes .GT. 1) THEN

    ALLOCATE(x_box(nbr_boxes))

    DO ibox = 1, nbr_boxes
       x_box(ibox) = REAL(nmols_box(ibox),DP)/REAL(nmols_tot,DP)
       IF ( ibox > 1 ) THEN
          x_box(ibox) = x_box(ibox) + x_box(ibox-1)
       END IF
    END DO

    randno = rranf()
    DO ibox = 1, nbr_boxes
       IF ( randno <= x_box(ibox)) EXIT
    END DO

    DEALLOCATE(x_box)

  ELSE

    ibox = 1

  END IF

  ! error check
  IF( nmols_box(ibox) == 0 ) THEN
     err_msg = ''
     err_msg(1) = 'No rotatable molecules in box ' // TRIM(Int_To_String(ibox))
     CALL Clean_Abort(err_msg, "Rotate")
  END IF

  ! Choose species based on the mol fraction, using Golden sampling
  ALLOCATE(x_species(nspecies))

  DO is = 1, nspecies
     IF ( max_rot(is,ibox) > 0. ) THEN
       x_species(is) = REAL(nmols(is,ibox), DP)/REAL(nmols_box(ibox),DP)
     ELSE
       x_species(is) = 0.0_DP
     END IF
     IF ( is > 1 ) THEN
       x_species(is) = x_species(is) + x_species(is-1)
     END IF
  END DO

  randno = rranf()
  DO is = 1, nspecies
     IF ( randno <= x_species(is) ) EXIT
  END DO

  DEALLOCATE(x_species)

  ! error check
  IF ( max_rot(is,ibox) == 0. ) THEN
     err_msg = ''
     err_msg(1) = 'No rotatable molecules of species ' // TRIM(Int_To_String(is)) &
               // ' in box ' // TRIM(Int_To_String(ibox))
     CALL Clean_Abort(err_msg, "Rotate")
  END IF

  ! Choose a molecule at random for rotation
  im = INT ( rranf() * nmols(is,ibox) ) + 1
  ! Get the LOCATE of imth molecule of species is in ibox
  lm = LOCATE(im,is,ibox)

  ! update the trial counters
  tot_trials(ibox) = tot_trials(ibox) + 1
  ntrials(is,ibox)%rotation = ntrials(is,ibox)%rotation + 1

  ! obtain the energy of the molecule before the move.  Note that due to
  ! this move, the interatomic energies such as vdw and electrostatics will
  ! change. Also the ewald_reciprocal energy will change but there will
  ! be no change in intramolecular energies.
  IF (l_pair_nrg) THEN
    CALL Store_Molecule_Pair_Interaction_Arrays(lm,is,ibox,E_vdw,E_qq)
  ELSE
    CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_vdw,E_qq,inter_overlap)
  END IF

  IF (inter_overlap)  THEN
     err_msg = ""
     err_msg(1) = "Attempted to rotate molecule " // TRIM(Int_To_String(im)) // &
                  " of species " // TRIM(Int_To_String(is))
     IF (nbr_boxes > 1) err_msg(1) = err_msg(1) // " in box " // TRIM(Int_To_String(ibox))
     err_msg(2) = "but the molecule energy is too high"
     IF (start_type(ibox) == "make_config" ) THEN
        err_msg(3) = "Try increasing Rcutoff_Low, increasing the box size, or "
        err_msg(4) = "decreasing the initial number of molecules"
     END IF
     CALL Clean_Abort(err_msg, "Rotate")
  END IF

  ! Store the old positions of the atoms
  CALL Save_Old_Cartesian_Coordinates(lm,is)

  ! change the Eulerian angles of the molecule
  IF(species_list(is)%linear) THEN
     CALL Rotate_Molecule_Eulerian
  ELSE
     CALL Rotate_Molecule_Axis(dx,dy,dz)
  END IF

  ! Now compute the energy of the molecule after the rotation. 
  CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_vdw_move,E_qq_move,inter_overlap)

  ! If an overlap is detected, immediately reject the move
  IF (inter_overlap) THEN
     
     CALL Revert_Old_Cartesian_Coordinates(lm,is)
     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(lm,is,ibox)
     
     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,A9)') &
             i_mcstep, 'rotate' , lm, is, ibox, accept,'overlap'
     END IF

  ELSE

     delta_e = 0.0_DP

     IF ((int_charge_sum_style(ibox) == charge_ewald) .AND. (has_charge(is))) THEN
        
        ALLOCATE(cos_mol_old(nvecs(ibox)),sin_mol_old(nvecs(ibox)))
        CALL Get_Position_Alive(lm,is,position)
     
        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        cos_mol_old(:) = cos_mol(1:nvecs(ibox),position)
        sin_mol_old(:) = sin_mol(1:nvecs(ibox),position)
        !$OMP END PARALLEL WORKSHARE

        CALL Update_System_Ewald_Reciprocal_Energy(lm,is,ibox,int_rotation,E_reciprocal_move)
        delta_e = E_reciprocal_move - energy(ibox)%ewald_reciprocal

     END IF
     
     delta_e = E_vdw_move - E_vdw + E_qq_move - E_qq + delta_e
     ! note that the difference in framework energy will be zero if 
     ! the simulation does not have a solid support, framework, wall etc.

     IF (int_sim_type == sim_nvt_min) THEN
        ! Accept only the moves that lower energy
        IF ( delta_e <= 0.0_DP) THEN
           accept = .TRUE.
        END IF
        
     ELSE

        ln_pacc = beta(ibox) * delta_e
        accept = accept_or_reject(ln_pacc)

     END IF
     
     IF ( accept ) THEN

        energy(ibox)%total = energy(ibox)%total + delta_e
        energy(ibox)%inter_vdw = energy(ibox)%inter_vdw + E_vdw_move - E_vdw
        energy(ibox)%inter_q   = energy(ibox)%inter_q + E_qq_move - E_qq

        IF ( int_charge_sum_style(ibox) == charge_ewald .AND. has_charge(is)) THEN
           energy(ibox)%ewald_reciprocal = E_reciprocal_move
        END IF

        nsuccess(is,ibox)%rotation = nsuccess(is,ibox)%rotation + 1
        nequil_success(is,ibox)%rotation = nequil_success(is,ibox)%rotation + 1
        
        ! Update the COM after rotation

        CALL Get_COM(lm,is)

        IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
        IF (ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
        IF (ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)
       
     ELSE

        ! Revert to the old coordinates of atoms and com of the molecule
        CALL Revert_Old_Cartesian_Coordinates(lm,is) 
        ! Reset the old cos_sum and sin_sum
        
        IF ((int_charge_sum_style(ibox) == charge_ewald) .AND. (has_charge(is))) THEN
           
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED) 
           cos_sum(:,ibox) = cos_sum_old(:,ibox)
           sin_sum(:,ibox) = sin_sum_old(:,ibox)
           cos_mol(1:nvecs(ibox),position) =cos_mol_old(:)
           sin_mol(1:nvecs(ibox),position) =sin_mol_old(:)
           !$OMP END PARALLEL WORKSHARE

           DEALLOCATE(cos_mol_old,sin_mol_old)

        END IF
        
        IF (l_pair_nrg)  CALL Reset_Molecule_Pair_Interaction_Arrays(lm,is,ibox)

     ENDIF

     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I3,X,L8,X,9X,X,F9.3)') &
             i_mcstep, 'rotate' , lm, is, ibox, accept, ln_pacc
     END IF

  END IF
 
  
  IF (MOD(ntrials(is,ibox)%rotation, nupdate) == 0) THEN

     IF (int_run_type == run_equil) THEN
        success_ratio = REAL(nequil_success(is,ibox)%rotation,DP)/REAL(nupdate,DP)
     ELSE
        success_ratio = REAL(nsuccess(is,ibox)%rotation,DP)/REAL(ntrials(is,ibox)%rotation,DP)
     END IF
 
     WRITE(logunit,'(X,I9,X,A10,X,5X,X,I3,X,I3,X,F8.5)',ADVANCE='NO') &
           i_mcstep, 'rotate', is, ibox, success_ratio

     IF (int_run_type == run_equil) THEN   
    
         nequil_success(is,ibox)%rotation = 0
        ! update the width of the species for equilibration run
        
        IF (success_ratio < 0.5) THEN
           ! decrease the width
           max_rot(is,ibox) = MAX(0.01_DP,0.95_DP*max_rot(is,ibox))
        ELSE
           max_rot(is,ibox) = MIN(PI,1.05_DP*max_rot(is,ibox))

        END IF
        WRITE(logunit,'(X,F9.5)',ADVANCE='NO') max_rot(is,ibox)
        
     END IF
     
     WRITE(logunit,*)

  END IF

!***************************************************************************
CONTAINS
!***************************************************************************

  SUBROUTINE Rotate_Molecule_Axis(dxrot,dyrot,dzrot)
  
    IMPLICIT NONE
    

    INTEGER :: axis, ia
    
    REAL(DP) :: psi1, psi2, psi3, cospsi1, cospsi2, cospsi3
    REAL(DP) :: sinpsi1, sinpsi2, sinpsi3
    REAL(DP) :: rot11, rot12, rot13, rot21, rot22, rot23, rot31, rot32, rot33
    REAL(DP) :: rxpnew, rypnew, rzpnew

    REAL(DP), DIMENSION(:), ALLOCATABLE :: dxrot, dyrot, dzrot

    psi1 = 0.0_DP
    psi2 = 0.0_DP
    psi3 = 0.0_DP

    ALLOCATE( dxrot(natoms(is)) , dyrot(natoms(is)) , dzrot(natoms(is)) )

    dxrot = 0.0_DP
    dyrot = 0.0_DP
    dzrot = 0.0_DP

    ! randomly choose an axis about which a rotation needs to be performed

    axis = INT ( 3.0_DP * rranf() ) + 1

    IF ( axis == 1 ) THEN
       
       psi1 = ( 2.0_DP * rranf() - 1.0_DP) * max_rot(is,ibox)
      
    ELSE IF ( axis == 2 ) THEN
       
       psi2 = ( 2.0_DP * rranf() - 1.0_DP) * max_rot(is,ibox)

    ELSE

       psi3 = ( 2.0_DP * rranf() - 1.0_DP) * max_rot(is,ibox)

    END IF

    ! Define the rotational matrix 


    cospsi1 = DCOS(psi1)
    cospsi2 = DCOS(psi2)
    cospsi3 = DCOS(psi3)

    sinpsi1 = DSIN(psi1)
    sinpsi2 = DSIN(psi2)
    sinpsi3 = DSIN(psi3)

    rot11   = cospsi1*cospsi3 - sinpsi1*cospsi2*sinpsi3
    rot12   = sinpsi1*cospsi3 + cospsi1*cospsi2*sinpsi3
    rot13   = sinpsi2*sinpsi3
    
    rot21   = -cospsi1*sinpsi3 - sinpsi1*cospsi2*cospsi3
    rot22   = -sinpsi1*sinpsi3 + cospsi1*cospsi2*cospsi3
    rot23   = sinpsi2*cospsi3
    
    rot31   = sinpsi1*sinpsi2
    rot32   = -cospsi1*sinpsi2

    rot33   = cospsi2

    ! Move the origin to the COM of this molecule

    atom_list(:,lm,is)%rxp = atom_list(:,lm,is)%rxp - molecule_list(lm,is)%xcom
    atom_list(:,lm,is)%ryp = atom_list(:,lm,is)%ryp - molecule_list(lm,is)%ycom
    atom_list(:,lm,is)%rzp = atom_list(:,lm,is)%rzp - molecule_list(lm,is)%zcom

    ! Apply the rotation matrix to these coordinates


    DO ia = 1, natoms(is)

       rxpnew = rot11*atom_list(ia,lm,is)%rxp + rot12*atom_list(ia,lm,is)%ryp + &
            rot13*atom_list(ia,lm,is)%rzp
       rypnew = rot21*atom_list(ia,lm,is)%rxp + rot22*atom_list(ia,lm,is)%ryp + &
            rot23*atom_list(ia,lm,is)%rzp
       rzpnew = rot31*atom_list(ia,lm,is)%rxp + rot32*atom_list(ia,lm,is)%ryp + &
            rot33*atom_list(ia,lm,is)%rzp

       dxrot(ia) = rxpnew - atom_list(ia,lm,is)%rxp
       dyrot(ia) = rypnew - atom_list(ia,lm,is)%ryp
       dzrot(ia) = rzpnew - atom_list(ia,lm,is)%rzp

       atom_list(ia,lm,is)%rxp = rxpnew
       atom_list(ia,lm,is)%ryp = rypnew
       atom_list(ia,lm,is)%rzp = rzpnew

    END DO
    
    
    ! Shift the origin back to the space fixed axes.

    atom_list(:,lm,is)%rxp = atom_list(:,lm,is)%rxp + molecule_list(lm,is)%xcom
    atom_list(:,lm,is)%ryp = atom_list(:,lm,is)%ryp + molecule_list(lm,is)%ycom
    atom_list(:,lm,is)%rzp = atom_list(:,lm,is)%rzp + molecule_list(lm,is)%zcom

  END SUBROUTINE Rotate_Molecule_Axis
  !-----------------------------------------------------------------------------------------------

  SUBROUTINE Rotate_Molecule_Eulerian

    ! This subroutine will rotate the molecule based on random pickings of eulerians
    ! and is based on the code from Dr. Maginn.
    
    IMPLICIT NONE

    REAL(DP) :: theta, phi, psi, rot11, rot12, rot13, rot21, rot22, rot23
    REAL(DP) :: rot31, rot32, rot33, rxpnew, rypnew, rzpnew
    
    INTEGER :: i, ia

    ! Pick random eulerians
    
    theta = ACOS(1.0_DP - 2.0_DP * rranf())
    phi   = (1.0_DP - 2.0_DP * rranf()) * PI
    psi   = (1.0_DP - 2.0_DP * rranf()) * PI

    ! shift the origin to the COM of the molecule

    atom_list(:,lm,is)%rxp = atom_list(:,lm,is)%rxp - molecule_list(lm,is)%xcom
    atom_list(:,lm,is)%ryp = atom_list(:,lm,is)%ryp - molecule_list(lm,is)%ycom
    atom_list(:,lm,is)%rzp = atom_list(:,lm,is)%rzp - molecule_list(lm,is)%zcom

    ! Construct the rotation matrix that needs to be applied to each of the vectors
    ! This is the A matrix in Goldstein notation

    rot11 = DCOS(psi) * DCOS(phi) - DCOS(theta) * DSIN(phi) * DSIN(psi)
    rot12 = DCOS(psi) * DSIN(phi) + DCOS(theta) * DCOS(phi) * DSIN(psi)
    rot13 = DSIN(psi) * DSIN(theta)

    rot21 = -DSIN(psi) * DCOS(phi) - DCOS(theta) * DSIN(phi) * DCOS(psi)
    rot22 = -DSIN(psi) * DSIN(phi) + DCOS(theta) * DCOS(phi) * DCOS(psi)
    rot23 = DCOS(psi) * DSIN(theta)

    rot31 = DSIN(theta) * DSIN(phi)
    rot32 = -DSIN(theta) * DCOS(phi)
    rot33 = DCOS(theta)

    ! Now rotate the relative vectors to obtain the new positions

    DO ia = 1, natoms(is)

       rxpnew = rot11*atom_list(ia,lm,is)%rxp & 
              + rot12*atom_list(ia,lm,is)%ryp &
              + rot13*atom_list(ia,lm,is)%rzp
       rypnew = rot21*atom_list(ia,lm,is)%rxp &
              + rot22*atom_list(ia,lm,is)%ryp &
              + rot23*atom_list(ia,lm,is)%rzp
       rzpnew = rot31*atom_list(ia,lm,is)%rxp &
              + rot32*atom_list(ia,lm,is)%ryp &
              + rot33*atom_list(ia,lm,is)%rzp

       atom_list(ia,lm,is)%rxp = rxpnew
       atom_list(ia,lm,is)%ryp = rypnew
       atom_list(ia,lm,is)%rzp = rzpnew

    END DO    

    ! Shift the origin back to (0,0,0)

    atom_list(:,lm,is)%rxp = atom_list(:,lm,is)%rxp + molecule_list(lm,is)%xcom
    atom_list(:,lm,is)%ryp = atom_list(:,lm,is)%ryp + molecule_list(lm,is)%ycom
    atom_list(:,lm,is)%rzp = atom_list(:,lm,is)%rzp + molecule_list(lm,is)%zcom
    
  END SUBROUTINE Rotate_Molecule_Eulerian

  !**********************************************************************
     
END SUBROUTINE Rotate

