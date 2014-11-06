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

  SUBROUTINE Grow_Molecules

    !****************************************************************************
    ! The subroutine generates initial configuration of molecules for box 1.
    ! It is based on an input geometry of molecule of each of the species.
    ! First a random point is selected in the simulation box for the placement
    ! of first atom of the molecule. The entire molecule is then constructed.
    ! Random orientation of the molecule is generated. A check is performed to
    ! ensure that the atoms of this molecules are at least
    !
    ! Called by
    !
    !   main
    !
    ! Revision history
    !
    !   12/10/13 : Beta version
    ! 
    !*************************************************************************
    USE Run_Variables
    USE Type_Definitions
    USE Random_Generators
    USE File_Names
    USE Energy_Routines
    USE Fragment_Growth
    USE Simulation_Properties, ONLY : Compute_Beads, Get_Index_Molecule
    USE IO_Utilities

    IMPLICIT NONE

    ! Local
    INTEGER :: is,im,ia, this_box, alive, this_im, n_start, n_end
    INTEGER :: is2, im2, ja, ibox, which_anchor
    INTEGER :: isdstart,isdend,isistart,isiend,ireac
    INTEGER :: i,m,box_start,box_end, ii
    INTEGER :: im_this,im_other, other_box, alive_this, alive_other, this_atom


    INTEGER, ALLOCATABLE,DIMENSION(:) :: frag_order

    REAL(DP) :: rxijp, ryijp, rzijp, rxij, ryij, rzij, rsq, attempt_prob, lambda_for_build
    REAL(DP) :: rand_lambda,lambda_ins,lambda_del
    REAL(DP) :: nrg_ring_frag_tot

    LOGICAL :: overlap, cbmc_overlap

    CHARACTER*120 :: init_config_file, init_config_file1
!********************************************************************************
    ! Place all the molecules in the box using CB growth.

    ! Needs to be written.

    ! Initialize the cfc amd existence flags. 
   
    overlap = .false. 
    molecule_list(:,:)%live = .false.

    DO is = 1, nspecies
       DO im = 1, nmolecules(is)
          locate(im,is) = im
       END DO
    END DO

    n_start = 0
    n_end = 0
    del_FLAG = .FALSE. 
    DO ibox = 1, nbr_boxes

       DO is=1,nspecies         

          IF (species_list(is)%fragment) THEN
             ALLOCATE(frag_order(nfragments(is)))
          END IF
          
          IF (ibox /=1) THEN
             n_start = SUM(nmol_actual(is,1:ibox-1))
             n_end = SUM(nmol_actual(is,1:ibox))
          ELSE
             n_start = 0
             n_end = nmol_actual(is,1)
          END IF

          DO im= n_start + 1, n_end

             alive = locate(im,is)
             
             molecule_list(alive,is)%which_box = ibox
             molecule_list(alive,is)%live = .true.
             
             molecule_list(alive,is)%molecule_type = int_normal
             
             molecule_list(alive,is)%cfc_lambda = 1.0_DP
             atom_list(:,alive,is)%exist = .TRUE.
             
             this_box = molecule_list(alive,is)%which_box
             ! Now let us insert the first atom of this molecule randomly in this box
             InsertionLOOP: DO
                
                
                ! we will grow the molecules if there are fragments in the molecules
                
                IF (species_list(is)%fragment) THEN
                   atom_list(:,alive,is)%exist = .FALSE.
                   cbmc_overlap = .FALSE.
                   get_fragorder = .TRUE.
                   lambda_for_build = molecule_list(alive,is)%cfc_lambda
                   CALL Build_Molecule(alive,is,this_box,frag_order,lambda_for_build,which_anchor, &
                        attempt_prob,nrg_ring_frag_tot,cbmc_overlap)
                   IF (cbmc_overlap) CYCLE InsertionLOOP
                   atom_list(:,alive,is)%exist = .TRUE.       

                ELSE
                   
                   IF(box_list(this_box)%box_shape == 'CUBIC') THEN
                      ! -- all the cell lengths are identical,
                      atom_list(1,alive,is)%rxp = (rranf() - 0.5_DP) * box_list(this_box)%length(1,1)
                      atom_list(1,alive,is)%ryp = (rranf() - 0.5_DP) * box_list(this_box)%length(2,2)
                      atom_list(1,alive,is)%rzp = (rranf() - 0.5_DP) * box_list(this_box)%length(3,3)                      
                   END IF
                   
                   ! insert the rest of the molecules
                   
                   DO ia = 2,natoms(is)
                      
                      atom_list(ia,alive,is)%rxp = atom_list(1,alive,is)%rxp + init_list(ia,1,is)%rxp - &
                           init_list(1,1,is)%rxp
                      atom_list(ia,alive,is)%ryp = atom_list(1,alive,is)%ryp + init_list(ia,1,is)%ryp - &
                           init_list(1,1,is)%ryp
                      atom_list(ia,alive,is)%rzp = atom_list(1,alive,is)%rzp + init_list(ia,1,is)%rzp - &
                           init_list(1,1,is)%rzp
                   END DO
                   
                   ! Obtain COM of the molecule

                   CALL Get_COM(alive,is)
                   
                   CALL Rotate_Molecule_Eulerian
                   
                END IF
                
                ! Obtain the new COM of the molecule
                
                CALL Get_COM(alive,is)
                
                ! Check for any overlaps with previously inserted molecules
                
                DO is2 = 1, nspecies
                   im2LOOP: DO im2 = 1, nmolecules(is2)
                      this_im = locate(im2,is2)
                      
                      ! skip the test for this molecule
                      
                      IF ( molecule_list(this_im,is2)%live) THEN
                         
                         IF ( molecule_list(this_im,is2)%which_box == this_box ) THEN
                            
                            IF ((is2 == is) .AND. (alive==this_im)) CYCLE im2LOOP
                            
                            ! check the overlap with all the atoms of this molecule
                            
                            DO ia = 1, natoms(is)
                               
                               DO ja = 1, natoms(is2)
                                  
                                  rxijp = atom_list(ia,alive,is)%rxp - atom_list(ja,this_im,is2)%rxp
                                  ryijp = atom_list(ia,alive,is)%ryp - atom_list(ja,this_im,is2)%ryp
                                  rzijp = atom_list(ia,alive,is)%rzp - atom_list(ja,this_im,is2)%rzp
                                  
                                  CALL Minimum_Image_Separation(this_box,rxijp,ryijp,rzijp,rxij,ryij,rzij)
                                  
                                  rsq = rxij * rxij + ryij * ryij + rzij * rzij
                                  
                                  IF (rsq < rcut_lowsq) CYCLE InsertionLOOP
                                  
                                  ! reject the insertion and go back for another trial
                                  
                               END DO
                               
                            END DO
                            
                         END IF
                         
                      END IF
                      
                   END DO im2LOOP
                   
                END DO
                
                ! If here then we did not find overlap of alive with any other molecule in the system. 
                ! exit and insert another molecule
                
                ! compute the distance of the psuedoatom farthest from the COM.
                
                CALL Compute_Max_Com_Distance(alive,is)
                ! write to a file for viewing
                
                WRITE(*,*) 'successfully inserted molecule', im
                EXIT
                
             END DO InsertionLOOP
             
          ENDDO  ! this ends the do loop from line 88
          
          nmols(is,ibox) = nmol_actual(is,ibox)
          
          IF (ALLOCATED(frag_order)) DEALLOCATE(frag_order)
          
       END DO ! ends do loop from line 74
       
       ! Here the configuration for this box has been generated
       ! Let's write this out
       init_config_file = inputfile//'.init_config'
       init_config_file = init_config_file//Int_To_String(ibox)

       CALL Name_Files(run_name,'.init_config.box',init_config_file)
       init_config_file1 = TRIM(init_config_file)//Int_To_String(ibox)

       OPEN(UNIT=configunit,file=init_config_file1)

       DO is = 1,nspecies
          DO im = 1, nmolecules(is)

             alive = locate(im,is)

             IF (molecule_list(alive,is)%live) THEN

                IF (molecule_list(alive,is)%which_box == ibox) THEN

                   DO ia = 1, natoms(is)

                      write(configunit,*) nonbond_list(ia,is)%element, atom_list(ia,im,is)%rxp, &
                           atom_list(ia,im,is)%ryp,atom_list(ia,im,is)%rzp

                   END DO

                END IF

             END IF

          END DO

       END DO

       CLOSE(UNIT=configunit)
       overlap = .FALSE.

    END DO  ! ends do loop from line 72 (nbr_boxes)
    
    DO ibox = 1, nbr_boxes

       IF(int_vdw_sum_style(ibox) == vdw_cut_tail) CALL Compute_Beads(ibox)

    END DO

  CONTAINS
    
    SUBROUTINE Rotate_Molecule_Eulerian
      
      ! This subroutine will rotate the molecule based on random pickings of eulerians
      ! and is based on the code from Dr. Maginn.
      
      IMPLICIT NONE
      
      REAL(DP) :: theta, phi, psi, rot11, rot12, rot13, rot21, rot22, rot23
      REAL(DP) :: rot31, rot32, rot33, rxpnew, rypnew, rzpnew
      
      INTEGER :: i, ia
      
      ! Pick random eulerians
      
      theta = ACOS(1.0_DP-2.0_DP*rranf())
      phi  = (1.0_DP - 2.0_DP * rranf()) * PI
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
         
         rxpnew = rot11*atom_list(ia,alive,is)%rxp + rot12*atom_list(ia,alive,is)%ryp + &
              rot13*atom_list(ia,alive,is)%rzp
         rypnew = rot21*atom_list(ia,alive,is)%rxp + rot22*atom_list(ia,alive,is)%ryp + &
              rot23*atom_list(ia,alive,is)%rzp
         rzpnew = rot31*atom_list(ia,alive,is)%rxp + rot32*atom_list(ia,alive,is)%ryp + &
              rot33*atom_list(ia,alive,is)%rzp
         
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
    
  END SUBROUTINE Grow_Molecules


  SUBROUTINE Update_Reservoir(is)
    !****************************************************************************
    ! The subroutine generates atomic positions for reservoir molecules
    !                 the event that the long range corrections are to be computed.
    !
    ! 01/06/10 (JS) : Ring biasing added while growing a molecule
    !********************************************************************************
    USE Run_Variables
    USE Type_Definitions
    USE Random_Generators
    USE File_Names
    USE Energy_Routines
    USE Fragment_Growth

    IMPLICIT NONE

    ! Local
    INTEGER :: is, i, ii, ia, alive, kappa_old
    INTEGER :: frag_start, frag_end, frag_total, this_box
    INTEGER, ALLOCATABLE, DIMENSION(:) :: frag_order
    REAL(DP) :: E_bond, E_angle, E_dihedral, E_improper, E_intra_vdw, E_intra_qq
    REAL(DP) :: d_bond, d_angle, d_dihedral, d_improper, d_intra_vdw, d_intra_qq, delta_e
    REAL(DP) :: factor
    REAL(DP) :: P_forward, e_prev, lambda_for_cut, nrg_ring_frag_forward
    CHARACTER(2) :: dummy_element
    LOGICAL :: cbmc_overlap, del_overlap, intra_overlap
    LOGICAL :: accept, accept_or_reject

    igas_flag = .TRUE.

    ALLOCATE(frag_order(nfragments(is)))
    kappa_old = kappa_dih
    kappa_dih = 1
    this_box = 1

    alive = nmolecules(is) + 1

    atom_list(:,alive,is)%exist = .TRUE.
    molecule_list(alive,is)%live = .TRUE. 
    molecule_list(alive,is)%which_box = nbr_boxes + 1

    IF(first_res_update) THEN

       OPEN(UNIT=sorbate_unit,FILE=sorbate_file(is),status='old')

       DO ia = 1, natoms(is)
                
          READ(sorbate_unit,*)dummy_element, &
                  atom_list_igas(ia,1,is)%rxp, &
                  atom_list_igas(ia,1,is)%ryp, &
                  atom_list_igas(ia,1,is)%rzp

       END DO

       CLOSE(UNIT=sorbate_unit)

       atom_list(:,alive,is)%rxp = atom_list_igas(:,1,is)%rxp
       atom_list(:,alive,is)%ryp = atom_list_igas(:,1,is)%ryp
       atom_list(:,alive,is)%rzp = atom_list_igas(:,1,is)%rzp

       CALL Get_COM(alive,is)

       CALL Compute_Molecule_Bond_Energy(alive,is,E_bond)
       CALL Compute_Molecule_Angle_Energy(alive,is,E_angle)
       CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihedral)
       CALL Compute_Molecule_Improper_Energy(alive,is,E_improper)
       CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw,E_intra_qq,intra_overlap)

       energy_igas(i,is)%bond = E_bond
       energy_igas(i,is)%angle = E_angle
       energy_igas(i,is)%dihedral = E_dihedral
       energy_igas(i,is)%improper = E_improper
       energy_igas(i,is)%intra_vdw = E_intra_vdw
       energy_igas(i,is)%intra_q = E_intra_qq
       energy_igas(i,is)%total = E_bond + E_angle + E_dihedral + E_improper + E_intra_vdw + E_intra_qq

    END IF

    DO i = 1, n_igas(is)

       IF ( .NOT. first_res_update) THEN

          atom_list(:,alive,is)%rxp = atom_list_igas(:,i,is)%rxp
          atom_list(:,alive,is)%ryp = atom_list_igas(:,i,is)%ryp
          atom_list(:,alive,is)%rzp = atom_list_igas(:,i,is)%rzp

          molecule_list(alive,is)%xcom = molecule_list_igas(i,is)%xcom
          molecule_list(alive,is)%ycom = molecule_list_igas(i,is)%ycom
          molecule_list(alive,is)%zcom = molecule_list_igas(i,is)%zcom

       END IF

       DO ii = 1, n_igas_moves(is)
 
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

             IF (cbmc_overlap .or. del_overlap) THEN
                atom_list(1:natoms(is),alive,is)%exist = .TRUE.
                accept = .FALSE.
                
                ! since accept is false, the coordinates will be reverted back to the original
                ! in the decision taken below.
                
             END IF

          END IF

          IF ( .NOT. (cbmc_overlap .OR. del_overlap)) THEN
             
             CALL Compute_Molecule_Bond_Energy(alive,is,E_bond)
             CALL Compute_Molecule_Angle_Energy(alive,is,E_angle)
             CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihedral)
             CALL Compute_Molecule_Improper_Energy(alive,is,E_improper)
             CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw,E_intra_qq,intra_overlap)
             
             d_bond = E_bond - energy_igas(i,is)%bond
             d_angle = E_angle - energy_igas(i,is)%angle
             d_dihedral = E_dihedral - energy_igas(i,is)%dihedral
             d_improper = E_improper - energy_igas(i,is)%improper
             d_intra_vdw = E_intra_vdw - energy_igas(i,is)%intra_vdw
             d_intra_qq = E_intra_qq - energy_igas(i,is)%intra_q
             
             delta_e = d_bond + d_angle + d_dihedral + d_improper + d_intra_vdw + d_intra_qq
             
             factor = beta(this_box) * delta_e
             
             accept = accept_or_reject(factor)

          END IF
             
          IF(accept) THEN
             
             energy_igas(i,is)%bond = E_bond
             energy_igas(i,is)%angle = E_angle
             energy_igas(i,is)%dihedral = E_dihedral
             energy_igas(i,is)%improper = E_improper
             energy_igas(i,is)%intra_vdw = E_intra_vdw
             energy_igas(i,is)%intra_q = E_intra_qq
             energy_igas(i,is)%total = E_bond + E_angle + E_dihedral + E_improper + E_intra_vdw + E_intra_qq
 
          ELSE

             CALL Revert_Old_Cartesian_Coordinates(alive,is)

          END IF 

       END DO

       atom_list_igas(:,i,is)%rxp = atom_list(:,alive,is)%rxp
       atom_list_igas(:,i,is)%ryp = atom_list(:,alive,is)%ryp
       atom_list_igas(:,i,is)%rzp = atom_list(:,alive,is)%rzp

       CALL Get_COM(alive,is)

       molecule_list_igas(i,is)%xcom = molecule_list(alive,is)%xcom
       molecule_list_igas(i,is)%ycom = molecule_list(alive,is)%ycom
       molecule_list_igas(i,is)%zcom = molecule_list(alive,is)%zcom

    END DO

    igas_flag = .FALSE.

  END SUBROUTINE Update_Reservoir

