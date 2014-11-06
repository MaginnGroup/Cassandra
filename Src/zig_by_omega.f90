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


SUBROUTINE Zig_By_Omega

  !******************************************************************************
  !
  ! This subroutine computes \frac{z_ig}{\Omega} to be used in acceptance rule
  ! for GCMC insertion or deletion moves
  !
  ! First written by Jindal Shah on 04/29/08
  !
  ! Gets called by 
  !
  ! Gcmc_Control.f90
  !
  ! 12/10/13 : Beta version
  !******************************************************************************

  USE Run_Variables
  USE Type_Definitions
  USE Random_Generators
  USE File_Names
  USE Energy_Routines
  USE Fragment_Growth

  IMPLICIT NONE

  INTEGER :: i, is, ia, alive, frag_start, frag_end, frag_total, this_box
  INTEGER, ALLOCATABLE, DIMENSION(:) :: frag_order

  REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper, E_intra_vdw, E_intra_qq
  REAL(DP) :: d_bond, d_angle, d_dihedral, d_improper, d_intra_vdw, d_intra_qq, delta_e
  REAL(DP) :: nrg_ring_frag_forward, P_forward, lambda_for_cut, e_prev, factor

  CHARACTER(2) :: dummy_element

  LOGICAL :: cbmc_overlap, del_overlap, intra_overlap
  LOGICAL :: accept, accept_or_reject

  ! for now set all the zig_by_omega to be 1.0_DP. It will be computed for flexible
  ! molecules later on.

  DO is = 1, nspecies

     IF (species_list(is)%species_type == 'SORBATE' .AND. zig_calc(is)) THEN

       alive = 1 
       ALLOCATE(frag_order(nfragments(is)))

       OPEN(UNIT=sorbate_unit,FILE=sorbate_file(is),status='old')

       DO ia = 1, natoms(is)
                
          READ(sorbate_unit,*)dummy_element, &
                  atom_list(ia,1,is)%rxp, &
                  atom_list(ia,1,is)%ryp, &
                  atom_list(ia,1,is)%rzp
                  atom_list(ia,1,is)%exist = .TRUE.

       END DO

       CLOSE(UNIT=sorbate_unit)

       molecule_list(alive,is)%live = .TRUE. 
       this_box = 1
       molecule_list(alive,is)%which_box = this_box
       molecule_list(alive,is)%cfc_lambda = 1.0_DP

       CALL Get_COM(alive,is)

       CALL Compute_Molecule_Bond_Energy(alive,is,E_bond)
       CALL Compute_Molecule_Angle_Energy(alive,is,E_angle)
       CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihedral)
       CALL Compute_Molecule_Improper_Energy(alive,is,E_improper)
       CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw,E_intra_qq,intra_overlap)

       energy(this_box)%bond = E_bond
       energy(this_box)%angle = E_angle
       energy(this_box)%dihedral = E_dihedral
       energy(this_box)%improper = E_improper
       energy(this_box)%intra_vdw = E_intra_vdw
       energy(this_box)%intra_q = E_intra_qq

       energy(this_box)%total = E_bond + E_angle + E_dihedral + E_improper + E_intra_vdw + E_intra_qq
     
       DO i = 1, 20000000    

          CALL Save_Old_Cartesian_Coordinates(alive,is)

          nrg_ring_frag_forward = 0.0_DP
          del_FLAG = .FALSE.
          cbmc_overlap = .FALSE.
          del_overlap = .FALSE.
          P_forward = 1.0_DP

          IF(nfragments(is) == 1) THEN
  
             CALL Single_Fragment_Regrowth(alive,is)
             frag_total = 1

          ELSE

             lambda_for_cut = molecule_list(alive,is)%cfc_lambda
             CALL Cut_Regrow(alive,is,frag_start,frag_end,frag_order,frag_total,lambda_for_cut, &
                     e_prev,P_forward, nrg_ring_frag_forward, cbmc_overlap, del_overlap)


             ! if cbmc_overlap or del_overlap is true then the molecule may have been only partially grown

             IF (cbmc_overlap .OR. del_overlap) THEN
                atom_list(1:natoms(is),alive,is)%exist = .TRUE.
                accept = .FALSE.

                ! Note that since accept is false, the coordinates will be reverted back
                ! to the original coordinates during a decision taken below.

             END IF

          END IF
          IF ( .NOT. (cbmc_overlap .or. del_overlap) ) THEN
             CALL Compute_Molecule_Bond_Energy(alive,is,E_bond)
             CALL Compute_Molecule_Angle_Energy(alive,is,E_angle)
             CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihedral)
             CALL Compute_Molecule_Improper_Energy(alive,is,E_improper)
             CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw,E_intra_qq,intra_overlap)
             
             d_bond = E_bond - energy(this_box)%bond
             d_angle = E_angle - energy(this_box)%angle
             d_dihedral = E_dihedral - energy(this_box)%dihedral
             d_improper = E_improper - energy(this_box)%improper
             d_intra_vdw = E_intra_vdw - energy(this_box)%intra_vdw
             d_intra_qq = E_intra_qq - energy(this_box)%intra_q
             
             delta_e = d_bond + d_angle + d_dihedral + d_improper + d_intra_vdw + d_intra_qq
             
             factor = beta(this_box) * delta_e
             
             accept = accept_or_reject(factor)

          END IF
             
          IF(accept) THEN
             
             energy(this_box)%bond = E_bond
             energy(this_box)%angle = E_angle
             energy(this_box)%dihedral = E_dihedral
             energy(this_box)%improper = E_improper
             energy(this_box)%intra_vdw = E_intra_vdw
             energy(this_box)%intra_q = E_intra_qq
             energy(this_box)%total = energy(this_box)%total + delta_e
             
          ELSE
             
             CALL Revert_Old_Cartesian_Coordinates(alive,is)
             
          END IF
          
          species_list(is)%zig_by_omega = species_list(is)%zig_by_omega + EXP(beta(this_box) * energy(this_box)%total)
          
       END DO
       
       species_list(is)%zig_by_omega = 1.0_DP / (species_list(is)%zig_by_omega / REAL(20000000,DP))
       
    ELSE   
      
       species_list(is)%zig_by_omega = 1.0_DP

     END IF

  write(*,*) 'Z \ Omega = ', species_list(is)%zig_by_omega

  END DO

END SUBROUTINE Zig_By_Omega
