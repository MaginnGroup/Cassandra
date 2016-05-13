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
MODULE Pair_Nrg_Routines
  !********************************************************************
  ! This module contains several routines are used when l_pair_nrg flag 
  ! is .TRUE. i.e., when pair interaction energy arrays are stored and
  ! used for efficient calculations.
  !
  ! Get_Position_Alive
  !
  !     Location of a molecule in one dimensional array
  !
  ! Store_Molecule_Pair_Interaction_Arrays
  !
  !     Pair energy interactions of a given molecule(s) is temporarily
  !     stored during a move
  !
  ! Reset_Molecule_Pair_Interaction_Arraya
  !
  !     Pair energy interactions of a given molecule(s) are restored 
  !     if the move is rejected.
  !
  ! 08/08/13 : Created the beta version
  !
  !!********************************************************************

  USE Type_Definitions
  USE Global_Variables

  IMPLICIT NONE

CONTAINS

  !*********************************************************
  SUBROUTINE Get_Position_Alive(alive,is,position)
    ! this subroutine takes in the locate id of a molecule
    ! of species 'is' and returns its position in a one
    ! dimensional array containing all the molecules of the
    ! the system
    !
    ! CALLED BY
    ! 
    !        angle_distortion
    !        cut_N_grow
    !        energy_routines
    !        gemc_nvt_volume
    !        gemc_particle_transfer 
    !        Store_Molecule_Pair_Interaction_Arrays
    !        Reset_Molecule_Pair_Interaction_Arrays
    !        rotate
    !        translate
    !        volume_change
    !******************************************************

    INTEGER, INTENT(IN) :: alive, is
    INTEGER, INTENT(OUT) :: position

    IF ( is == 1) THEN
       position = alive
    ELSE
       position = SUM(max_molecules(1:is-1)) + alive
    END IF

  END SUBROUTINE Get_Position_Alive
  !********************************************************

  !********************************************************
  SUBROUTINE Store_Molecule_Pair_Interaction_Arrays(alive,is,this_box,E_vdw,E_qq, &
       n_cls_mol, id_cls_mol, is_cls_mol, box_cls_mol, box_nrg_vdw, box_nrg_qq)

    !******************************************************
    ! This subroutine stores energy pair interaction arrays
    ! for one or more molecules. It gets called by routines that perturb
    ! one or more molecules
    !
    ! The pair interactions are stored in global temp arrays. The routine
    ! also returns total intermolecular vdw and real space electrostatic interaction
    ! energies of the molecule being perturbed
    !
    ! CALLED BY
    !
    !        angle_distortion.f90
    !        cut_N_grow
    !        gemc_particle_transfer
    !        rotate
    !        translate
    !
    ! CALLS
    !
    !        Get_Position_Alive
    !
    ! INPUT VARIABLES
    !
    ! alive:       LOCATE of the molecule.
    !              Not used if we want to get interaction energies of two
    !              molecules (use a dummy variable in this case).
    ! is:          species type of the molecule. Not used if we want to get
    !              energies for two molecules (use dummy var)
    ! this_box:    box where alive is located.
    !              Not used if we want to get energies
    !              for two molecules (use dummy var)
    ! n_cls_mol:   optional argument. This equals the number of particles 
    !              for which we want to compute interaction energies. Typically
    !              used if we are interested in getting interaction energies of
    !              two different molecules (i.e. DFC).
    ! id_cls_mol:  optional argument. This equals the ID of particles for which we want to compute
    !              interaction energies. Typically used in DFC.
    ! is_cls_mol:  optional argument. This equals the species number of particles for which we want
    !              to compute interaction energies. Typically used in DFC.
    ! box_cls_mol: optional argument. This contains the boxes where the
    !              molecules of interest are located
    !
    ! OUTPUT VARIABLES
    !
    ! E_vdw:       Intermolecular VDW interaction energy of a single molecule 
    ! E_qq:        Intermolecular electrostatic energy of a single molecule
    ! box_nrg_vdw: optional argument. VDW energies of molecules in each box
    ! box_nrg_qq:  optional argument. Electrostatic energies of molecules in each box
    !
    !
    ! USAGE EXAMPLE:
    !
    ! For a single molecule in one box, we can store the pair arrays by 
    !
    ! Store_Molecule_Pair_Interaction_Arrays(alive,is,this_box,E_vdw,E_qq)
    !
    ! Note that we'll also get interaction energies through E_vdw, E_qq.
    !
    ! For two molecules in the same box, we can store the pair energies by
    !
    ! Store_Molecule_Pair_Interaction_Arrays(dummy1,dummy2,dummy3, E_vdw,
    ! E_qq, 2, [16 73], [1 2], [1 1], box_nrg_vdw, box_nrg_qq)
    !
    ! In this case, dummy* variables are ignored. E_vdw E_qq will contain
    ! inter vdw and electrostatic interactions of the last specified molecule
    ! (i.e. 73.). '2' Stands for the number
    ! of molecules for which we want to save arrays. [16 73] is a vector
    ! that contains the alive IDS of the molecules of interest. [1 1] is a
    ! vector that contains species ID of the molecules of interest.
    ! Box_nrg_vdw and box_nrg_qq are vectors that will contain the interaction
    ! energies for each box. In this example, the numbers will be the same.
    !******************************************************

    INTEGER :: alive, is, this_box
    INTEGER, OPTIONAL:: n_cls_mol, id_cls_mol(:), is_cls_mol(:)
    INTEGER, OPTIONAL :: box_cls_mol(:)

    INTEGER :: n_mols, stride, position, imol
    REAL(DP), INTENT(OUT) :: E_vdw, E_qq

    REAL(DP), OPTIONAL, INTENT(OUT) :: box_nrg_vdw(:), box_nrg_qq(:)

    INTEGER :: locate_1, locate_2, this_species, this_im, locate_im


    IF ( .NOT. present(n_cls_mol)) THEN
       ! only a single molecule storage is necessary
       n_mols = 1
       ALLOCATE(pair_vdw_temp(SUM(max_molecules)))
       ALLOCATE(pair_qq_temp(SUM(max_molecules)))

    ELSE

       ! storage for multiple molecules is required
       ! we will form an array that is n_cls_mol * SUM(max_molecules)
       n_mols = n_cls_mol

       ALLOCATE(pair_vdw_temp(n_mols*SUM(max_molecules)))
       ALLOCATE(pair_qq_temp(n_mols*SUM(max_molecules)))

       IF ( present(box_cls_mol)) THEN
          box_nrg_vdw(:) = 0.0_DP
          box_nrg_qq(:) = 0.0_DP
       END IF
       
    END IF

    E_vdw = 0.0_DP
    E_qq = 0.0_DP
   
    pair_vdw_temp(:) = 0.0_DP
    pair_qq_temp(:) = 0.0_DP

   

    DO imol = 1, n_mols
       
       IF ( .NOT. present(n_cls_mol)) THEN
       
          CALL Get_Position_Alive(alive,is,locate_1)

       ELSE

          alive = id_cls_mol(imol)
          is = is_cls_mol(imol)

       END IF

       IF ( present(box_cls_mol) ) THEN

          this_box = box_cls_mol(imol)
          E_vdw = 0.0_DP
          E_qq = 0.0_DP

       END IF
          
       CALL Get_Position_Alive(alive,is,locate_1)

       !Get_Position_Alive is used in conjunction with pair_vdw/pair_qq arrays

       stride = (imol-1) * SUM(max_molecules)
       
       speciesLoop: DO this_species = 1, nspecies
          
          molidLoop: DO this_im = 1, nmols(this_species, this_box)
             
             locate_im = locate(this_im,this_species,this_box)
             
             IF (molecule_list(locate_im,this_species)%live) THEN
                   
                   CALL Get_Position_Alive(locate_im,this_species,locate_2)
                   
                   position = locate_2 + stride
                   
                   pair_vdw_temp(position) = pair_nrg_vdw(locate_2,locate_1)
                   pair_qq_temp(position) = pair_nrg_qq(locate_2,locate_1)
                   
                   E_vdw = E_vdw + pair_vdw_temp(position)
                   E_qq = E_qq + pair_qq_temp(position)
                   
                
             END IF
             
          END DO molidLoop
          
       END DO speciesLoop

       IF (present(box_cls_mol)) THEN
          box_nrg_vdw(this_box) = box_nrg_vdw(this_box) + E_vdw
          box_nrg_qq(this_box) = box_nrg_qq(this_box) + E_qq
       END IF

    END DO

  END SUBROUTINE Store_Molecule_Pair_Interaction_Arrays
  !****************************************************************

  SUBROUTINE Reset_Molecule_Pair_Interaction_Arrays(alive,is,this_box, n_cls_mol, &
       id_cls_mol,is_cls_mol, box_cls_mol)
    !*****************************************************************
    ! The subroutine resets the pair interaction arrays after a move
    ! involving a molecule perturbation has been rejected 
    !
    ! CALLED BY
    !
    !        angle_distortion
    !        cut_N_grow
    !        GEMC_Particle_Transfer
    !        Rotate
    !        Translate
    !
    ! CALLS
    !        Get_Position_Alive
    !
    !******************************************************************
    
    INTEGER :: alive,is,this_box

    INTEGER, OPTIONAL :: n_cls_mol, id_cls_mol(:), is_cls_mol(:)
    INTEGER, OPTIONAL :: box_cls_mol(:)

    INTEGER :: n_mols, imol, stride

    INTEGER :: locate_1, this_species, this_im, locate_im, locate_2


    IF ( present(n_cls_mol)) THEN
       n_mols = n_cls_mol
    ELSE
       n_mols = 1 
    END IF

    DO imol = 1, n_mols

       IF ( present(n_cls_mol) ) THEN
          alive = id_cls_mol(imol)
          is = is_cls_mol(imol)
       END IF

       IF ( present(box_cls_mol) ) THEN
          this_box = box_cls_mol(imol)
       END IF
       
       
       CALL Get_Position_Alive(alive,is,locate_1)

       stride = (imol - 1) * SUM(max_molecules)
       
       DO this_species = 1, nspecies
          
          DO this_im = 1, nmols(this_species, this_box)
             
             locate_im = locate(this_im,this_species,this_box)
             
             IF (molecule_list(locate_im,this_species)%live) THEN
                
                   
                   CALL Get_Position_Alive(locate_im,this_species,locate_2)
                   
                   pair_nrg_vdw(locate_1,locate_2) = pair_vdw_temp(locate_2 + stride)
                   pair_nrg_vdw(locate_2,locate_1) = pair_vdw_temp(locate_2 + stride)
                   
                   pair_nrg_qq(locate_1,locate_2) = pair_qq_temp(locate_2 + stride)
                   pair_nrg_qq(locate_2,locate_1) = pair_qq_temp(locate_2 + stride)
                
                
             END IF
             
          END DO
          
       END DO
       
    END DO
    DEALLOCATE(pair_vdw_temp)
    DEALLOCATE(pair_qq_temp)
    
  END SUBROUTINE Reset_Molecule_Pair_Interaction_Arrays
    

END MODULE Pair_Nrg_Routines
