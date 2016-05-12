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
SUBROUTINE Rotate(this_box)
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

  IMPLICIT NONE

!  !$ include 'omp_lib.h'
  ! Arguments
  INTEGER  :: this_box   ! box index

  ! Local delcarations
  INTEGER  :: ibox       ! box index
  INTEGER  :: is         ! species index
  INTEGER  :: im, alive  ! molecule indices
  INTEGER  :: total_mols ! number of molecules in the system

  REAL(DP) :: nmols_box(nbr_boxes)
  REAL(DP), ALLOCATABLE :: x_box(:), x_species(:)
  REAL(DP) :: rand_no
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

  ! Sum the total number of molecules 
  total_mols = 0 ! sum over species, box
  DO ibox = 1, nbr_boxes
    nmols_box(ibox) = 0
    DO is = 1, nspecies
      ! Only count mobile species
      IF ( max_rot(is,ibox) > 0. ) THEN
        total_mols = total_mols + nmols(is,ibox)
        nmols_box(ibox) = nmols_box(ibox) + nmols(is,ibox)
      END IF
    END DO
  END DO

  ! If there are no molecules then return
  IF (total_mols == 0) RETURN

  ! If needed, choose a box based on its total mol fraction
  IF(nbr_boxes .GT. 1) THEN

    ALLOCATE(x_box(nbr_boxes))

    DO ibox = 1, nbr_boxes
       x_box(ibox) = REAL(nmols_box(ibox),DP)/REAL(total_mols,DP)
       IF ( ibox > 1 ) THEN
          x_box(ibox) = x_box(ibox) + x_box(ibox-1)
       END IF
    END DO

    rand_no = rranf()
    DO ibox = 1, nbr_boxes
       IF ( rand_no <= x_box(ibox)) EXIT
    END DO

    this_box = ibox
    DEALLOCATE(x_box)

  ELSE

    this_box = 1

  END IF

  ! If there are no molecules in this box then return
  IF( nmols_box(this_box) == 0 ) RETURN

  ! Choose species based on the mol fraction, using Golden sampling
  ALLOCATE(x_species(nspecies))

  DO is = 1, nspecies
     IF ( max_rot(is,this_box) > 0. ) THEN
       x_species(is) = REAL(nmols(is,this_box), DP)/REAL(nmols_box(this_box),DP)
     ELSE
       x_species(is) = 0.0_DP
     END IF
     IF ( is > 1 ) THEN
       x_species(is) = x_species(is) + x_species(is-1)
     END IF
  END DO

  rand_no = rranf()
  DO is = 1, nspecies
     IF ( rand_no <= x_species(is) ) EXIT
  END DO

  DEALLOCATE(x_species)

  ! If the molecule can't move then return
  IF ( max_rot(is,this_box) == 0. ) RETURN

  ! Choose a molecule at random for rotation
  im = INT ( rranf() * nmols(is,this_box) ) + 1
  ! Get the index of imth molecule of species is in this_box
  CALL Get_Index_Molecule(this_box,is,im,alive)

  ! update the trial counters
  tot_trials(this_box) = tot_trials(this_box) + 1
  ntrials(is,this_box)%rotation = ntrials(is,this_box)%rotation + 1

  ! obtain the energy of the molecule before the move.  Note that due to
  ! this move, the interatomic energies such as vdw and electrostatics will
  ! change. Also the ewald_reciprocal energy will change but there will
  ! be no change in intramolecular energies.
  IF (l_pair_nrg) THEN
    CALL Store_Molecule_Pair_Interaction_Arrays(alive,is,this_box,E_vdw,E_qq)
  ELSE
    CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_vdw,E_qq,inter_overlap)
  END IF

  IF (inter_overlap)  THEN
     WRITE(*,*) 'Disaster, overlap in the old configruation'
     WRITE(*,*) 'Rotation'
     WRITE(*,*) alive, is, this_box
  END IF       

  ! Store the old positions of the atoms
  CALL Save_Old_Cartesian_Coordinates(alive,is)

  ! change the Eulerian angles of the molecule
  IF(species_list(is)%linear) THEN
     CALL Rotate_Molecule_Eulerian
  ELSE
     CALL Rotate_Molecule_Axis(dx,dy,dz)
  END IF

  ! Now compute the energy of the molecule after the rotation. 
  CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is,E_vdw_move,E_qq_move,inter_overlap)

  ! If an overlap is detected, immediately reject the move
  IF (inter_overlap) THEN
     
     CALL Revert_Old_Cartesian_Coordinates(alive,is)
     IF (l_pair_nrg) CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,this_box)
     
  ELSE

     delta_e = 0.0_DP

     IF ((int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is))) THEN
        
        ALLOCATE(cos_mol_old(nvecs(this_box)),sin_mol_old(nvecs(this_box)))
        CALL Get_Position_Alive(alive,is,position)
     
        !$OMP PARALLEL WORKSHARE DEFAULT(SHARED)
        cos_mol_old(:) = cos_mol(1:nvecs(this_box),position)
        sin_mol_old(:) = sin_mol(1:nvecs(this_box),position)
        !$OMP END PARALLEL WORKSHARE

        CALL Update_System_Ewald_Reciprocal_Energy(alive,is,this_box,int_rotation,E_reciprocal_move)
        delta_e = E_reciprocal_move - energy(this_box)%ewald_reciprocal

     END IF
     
     delta_e = E_vdw_move - E_vdw + E_qq_move - E_qq + delta_e
     ! note that the difference in framework energy will be zero if 
     ! the simulation does not have a solid support, framework, wall etc.

     IF (int_sim_type == sim_nvt_min) THEN
        ! Accept only the moves that lower energy
        IF ( delta_e <= 0.0_DP) THEN
           accept = .TRUE.
        ELSE
           accept = .FALSE.
        END IF
        
     ELSE

        ln_pacc = beta(this_box) * delta_e
        accept = accept_or_reject(ln_pacc)

     END IF
     
     IF ( accept ) THEN

        energy(this_box)%total = energy(this_box)%total + delta_e
        energy(this_box)%inter_vdw = energy(this_box)%inter_vdw + E_vdw_move - E_vdw
        energy(this_box)%inter_q   = energy(this_box)%inter_q + E_qq_move - E_qq

        IF ( int_charge_sum_style(this_box) == charge_ewald .AND. has_charge(is)) THEN
           energy(this_box)%ewald_reciprocal = energy(this_box)%ewald_reciprocal + E_reciprocal_move
        END IF

        nsuccess(is,this_box)%rotation = nsuccess(is,this_box)%rotation + 1
        nequil_success(is,this_box)%rotation = nequil_success(is,this_box)%rotation + 1
        
        ! Update the COM after rotation

        CALL Get_COM(alive,is)

        IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)
        IF (ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
        IF (ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)
       
     ELSE

        ! Revert to the old coordinates of atoms and com of the molecule
        CALL Revert_Old_Cartesian_Coordinates(alive,is) 
        ! Reset the old cos_sum and sin_sum
        
        IF ((int_charge_sum_style(this_box) == charge_ewald) .AND. (has_charge(is))) THEN
           
           !$OMP PARALLEL WORKSHARE DEFAULT(SHARED) 
           cos_sum(:,this_box) = cos_sum_old(:,this_box)
           sin_sum(:,this_box) = sin_sum_old(:,this_box)
           cos_mol(1:nvecs(this_box),position) =cos_mol_old(:)
           sin_mol(1:nvecs(this_box),position) =sin_mol_old(:)
           !$OMP END PARALLEL WORKSHARE

           DEALLOCATE(cos_mol_old,sin_mol_old)

        END IF
        
        IF (l_pair_nrg)  CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,this_box)

     ENDIF

  END IF
 
  
  IF (MOD(ntrials(is,this_box)%rotation, nupdate) == 0) THEN

     IF (int_run_style == run_equil) THEN
        success_ratio = REAL(nequil_success(is,this_box)%rotation,DP)/REAL(nupdate,DP)
     ELSE
        success_ratio = REAL(nsuccess(is,this_box)%rotation,DP)/REAL(ntrials(is,this_box)%rotation,DP)
     END IF

     WRITE(logunit,*)
     WRITE(logunit,'(A,I3,A,I1,A,F8.5)') 'Success ratio, rotation of species ', is,' in box ',this_box, ' : ', success_ratio

     IF (int_run_style == run_equil) THEN   
    
         nequil_success(is,this_box)%rotation = 0
        ! update the width of the species for equilibration run
        
        IF (success_ratio < 0.5) THEN
           ! decrease the width
           max_rot(is,this_box) = MAX(0.01_DP,0.95_DP*max_rot(is,this_box))
        ELSE
           max_rot(is,this_box) = MIN(PI,1.05_DP*max_rot(is,this_box))

        END IF
        WRITE(logunit,'(A,I3,A,I1,A,F8.5)') 'Maximum width, rotation of species ', is,' in box ', this_box, ' : ', max_rot(is,this_box)
        
     END IF
     
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
       
       psi1 = ( 2.0_DP * rranf() - 1.0_DP) * max_rot(is,this_box)
      
    ELSE IF ( axis == 2 ) THEN
       
       psi2 = ( 2.0_DP * rranf() - 1.0_DP) * max_rot(is,this_box)

    ELSE

       psi3 = ( 2.0_DP * rranf() - 1.0_DP) * max_rot(is,this_box)

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
!    IF (alive==30) THEN
 !   WRITE(*,*) 'molecule',molecule_list(alive,is)%ycom
  !  END IF
  !WRITE(*,*) molecule_list(30,1)%ycom

    ! Move the origin to the COM of this molecule

    atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp - molecule_list(alive,is)%xcom
    atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp - molecule_list(alive,is)%ycom
    atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp - molecule_list(alive,is)%zcom

    ! Apply the rotation matrix to these coordinates


    DO ia = 1, natoms(is)

       rxpnew = rot11*atom_list(ia,alive,is)%rxp + rot12*atom_list(ia,alive,is)%ryp + &
            rot13*atom_list(ia,alive,is)%rzp
       rypnew = rot21*atom_list(ia,alive,is)%rxp + rot22*atom_list(ia,alive,is)%ryp + &
            rot23*atom_list(ia,alive,is)%rzp
       rzpnew = rot31*atom_list(ia,alive,is)%rxp + rot32*atom_list(ia,alive,is)%ryp + &
            rot33*atom_list(ia,alive,is)%rzp

       dxrot(ia) = rxpnew - atom_list(ia,alive,is)%rxp
       dyrot(ia) = rypnew - atom_list(ia,alive,is)%ryp
       dzrot(ia) = rzpnew - atom_list(ia,alive,is)%rzp

       atom_list(ia,alive,is)%rxp = rxpnew
       atom_list(ia,alive,is)%ryp = rypnew
       atom_list(ia,alive,is)%rzp = rzpnew

    END DO
    
    
    ! Shift the origin back to the space fixed axes.

    atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp + molecule_list(alive,is)%xcom
    atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp + molecule_list(alive,is)%ycom
    atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp + molecule_list(alive,is)%zcom

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

    atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp - molecule_list(alive,is)%xcom
    atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp - molecule_list(alive,is)%ycom
    atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp - molecule_list(alive,is)%zcom

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

       rxpnew = rot11*atom_list(ia,alive,is)%rxp & 
              + rot12*atom_list(ia,alive,is)%ryp &
              + rot13*atom_list(ia,alive,is)%rzp
       rypnew = rot21*atom_list(ia,alive,is)%rxp &
              + rot22*atom_list(ia,alive,is)%ryp &
              + rot23*atom_list(ia,alive,is)%rzp
       rzpnew = rot31*atom_list(ia,alive,is)%rxp &
              + rot32*atom_list(ia,alive,is)%ryp &
              + rot33*atom_list(ia,alive,is)%rzp

       atom_list(ia,alive,is)%rxp = rxpnew
       atom_list(ia,alive,is)%ryp = rypnew
       atom_list(ia,alive,is)%rzp = rzpnew

    END DO    

    ! Shift the origin back to (0,0,0)

    atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp + molecule_list(alive,is)%xcom
    atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp + molecule_list(alive,is)%ycom
    atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp + molecule_list(alive,is)%zcom
    
  END SUBROUTINE Rotate_Molecule_Eulerian

  !**********************************************************************
     
END SUBROUTINE Rotate

