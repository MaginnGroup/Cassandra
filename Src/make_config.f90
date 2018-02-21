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

  SUBROUTINE Make_Config(ibox)

    !***************************************************************************
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
    !***************************************************************************
    USE Global_Variables
    USE Type_Definitions
    USE Random_Generators
    USE File_Names
    USE Energy_Routines
    USE Fragment_Growth
    USE Simulation_Properties, ONLY : Compute_Beads
    USE IO_Utilities

    IMPLICIT NONE

    ! Input
    INTEGER, INTENT(IN) :: ibox

    ! Local
    INTEGER :: is,im,ia, this_box, alive, this_im, im_base, locate_base
    INTEGER :: is2, im2, ja, which_anchor
    INTEGER :: isdstart,isdend,isistart,isiend,ireac
    INTEGER :: i,m,box_start,box_end, ii
    INTEGER :: im_this,im_other, other_box, alive_this, alive_other, this_atom
    INTEGER :: i_frac


    INTEGER, ALLOCATABLE,DIMENSION(:) :: frag_order

    REAL(DP) :: rxijp, ryijp, rzijp, rxij, ryij, rzij, rsq, lambda_for_build
    REAL(DP) :: ln_pseq, ln_pbias
    REAL(DP) :: rand_lambda,lambda_ins,lambda_del
    REAL(DP) :: nrg_ring_frag_tot

    LOGICAL :: overlap, cbmc_overlap

    CHARACTER*120 :: init_config_file, init_config_file1
!*******************************************************************************
    ! Place molecules in the box using CB growth.

    overlap = .false. 
    del_flag = .FALSE. 

    DO is=1,nspecies         

       IF (nmols_to_make(is,ibox) > 0) &
         WRITE(logunit,*) 'Inserting ', nmols_to_make(is,ibox), &
            'molecules of species ', is
       
       IF (species_list(is)%fragment) THEN
          ALLOCATE(frag_order(nfragments(is)))
       END IF
       
       ! if there are already molecules in the box, start at next number
       im_base = nmols(is,ibox)

       ! LOCATE numbering is continuous across all boxes
       locate_base = SUM(nmols(is,1:nbr_boxes))

       DO im = im_base+1, im_base+nmols_to_make(is,ibox)
          nmols(is,ibox) = nmols(is,ibox) + 1
          locate(im,is,ibox) = im -im_base + locate_base
          alive = locate(im,is,ibox)
          molecule_list(alive,is)%which_box = ibox
          molecule_list(alive,is)%live = .true.
          
          molecule_list(alive,is)%molecule_type = int_normal
          
          molecule_list(alive,is)%frac = 1.0_DP
          atom_list(:,alive,is)%exist = .TRUE.
          
          ! Now let us insert the molecule randomly in this box
          InsertionLOOP: DO
             
             ! we will grow the molecules if there are fragments in the molecules
             IF (species_list(is)%fragment) THEN
                atom_list(:,alive,is)%exist = .FALSE.
                cbmc_overlap = .FALSE.
                get_fragorder = .TRUE.
                ln_pseq = 0.0_DP
                ln_pbias = 0.0_DP
                lambda_for_build = molecule_list(alive,is)%frac
                CALL Build_Molecule(alive,is,ibox,frag_order, &
                        lambda_for_build,ln_pseq,ln_pbias, &
                        nrg_ring_frag_tot,cbmc_overlap)
                IF (cbmc_overlap) CYCLE InsertionLOOP
                atom_list(:,alive,is)%exist = .TRUE.       

             ELSE
                
                IF(box_list(ibox)%box_shape == 'cubic') THEN
                   ! -- all the cell lengths are identical,
                   atom_list(1,alive,is)%rxp = (rranf() - 0.5_DP) * box_list(ibox)%length(1,1)
                   atom_list(1,alive,is)%ryp = (rranf() - 0.5_DP) * box_list(ibox)%length(2,2)
                   atom_list(1,alive,is)%rzp = (rranf() - 0.5_DP) * box_list(ibox)%length(3,3)                      
                END IF
                
                ! insert the rest of the molecule
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
                im2LOOP: DO im2 = 1, nmols(is2,ibox)
                   this_im = locate(im2,is2,ibox)
                   
                   ! skip the test for this molecule
                   IF ((is2 == is) .AND. (alive==this_im)) CYCLE im2LOOP
                   
                   ! check the overlap with all the atoms of this molecule
                   DO ia = 1, natoms(is)
                      
                      DO ja = 1, natoms(is2)
                         
                         rxijp = atom_list(ia,alive,is)%rxp - atom_list(ja,this_im,is2)%rxp
                         ryijp = atom_list(ia,alive,is)%ryp - atom_list(ja,this_im,is2)%ryp
                         rzijp = atom_list(ia,alive,is)%rzp - atom_list(ja,this_im,is2)%rzp
                         
                         CALL Minimum_Image_Separation(ibox,rxijp,ryijp,rzijp,rxij,ryij,rzij)
                         
                         rsq = rxij * rxij + ryij * ryij + rzij * rzij
                         
                         IF (rsq < rcut_lowsq) CYCLE InsertionLOOP
                         ! reject the insertion and go back for another trial
                         
                      END DO
                      
                   END DO
                   
                END DO im2LOOP
                
             END DO
             
             ! If here then we did not find overlap of alive with any other molecule in the system. 
             ! exit and insert another molecule
             
             ! compute the distance of the psuedoatom farthest from the COM.
             CALL Compute_Max_Com_Distance(alive,is)
             
             EXIT
             
          END DO InsertionLOOP
          
       ENDDO  ! this ends the do loop from line 88
       
       IF (ALLOCATED(frag_order)) DEALLOCATE(frag_order)
       
    END DO ! ends do loop from line 74
    
    ! Here the configuration for this box has been generated
    ! Let's write this out
    CALL Write_Coords(ibox)
    overlap = .FALSE.

    IF(int_vdw_sum_style(ibox) == vdw_cut_tail) CALL Compute_Beads(ibox)

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
    
  END SUBROUTINE Make_Config



