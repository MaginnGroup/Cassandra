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

MODULE Fragment_Growth

  !*****************************************************************************
  ! The module performs the tasks related to insertion, deletion, cut and 
  ! regrowth of molecules based on fragment sampling. It also contains a routine
  ! to obtain order of fragment placement.
  !
  ! CONTAINS
  !      
  ! Build_Molecule
  ! Build_Rigid_Fragment
  ! Cut_Regrow
  ! Fragment_Order
  ! Fragment_Placement
  ! Get_Aligner_Hanger
  ! Single_Fragment_Regrowth
  ! Get_Common_Fragment_Atoms
  !
  ! Used by
  !
  !   chempot
  !   cut_n_grow
  !   deletion
  !   gemc_particle_transfer
  !   make_config
  !   insertion
  !   main
  !   ring_fragment_driver
  ! 
  ! Revision history
  !
  !   12/10/13 : Beta Release
  !   Version 1.1
  !     05/15/15 Modified Fragment_Order to choose the next fragment with 
  !              uniform probability
  !     
  !*****************************************************************************

  USE Global_Variables
  USE Energy_Routines
  USE Random_Generators
  USE Read_Write_Checkpoint

  IMPLICIT NONE

CONTAINS

!*******************************************************************************
SUBROUTINE Build_Molecule(this_im,is,this_box,frag_order,this_lambda, &
              ln_pseq,ln_pbias,nrg_ring_frag_total,cbmc_overlap)
!*******************************************************************************
!
! PURPOSE: build the molecule from scratch
!
! First written by Jindal Shah on 07/18/07
!
! 05/25/08 (JS) : First committed to the repository
!
! 01/20/0 (JS) : Boltzmann weight of the trial is computed irrespective of
!                its energy. Previously, an energy cutoff was used to set
!                the Boltzmann weight
!
! 06/14/16 (RGM) : First fragment selected based on natoms_in_frag / natoms_in_all_frags
!
! DESCRIPTION: This subroutine performs the following steps
!
! Step 1) Select which fragment will be inserted first
! Step 2) Choose a conformation for the first fragment
! Step 3) Rotate the first fragment
! Step 4) Choose kappa_ins positions for the first fragment's COM
! Step 5) Choose kappa_dih orientations for each additional fragment
!
!*******************************************************************************

  USE Rotation_Routines
  USE IO_Utilities
  USE File_Names
  USE Energy_Routines

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  ! Arguments
  INTEGER :: this_im  ! molecule index
  INTEGER :: is       ! species index
  INTEGER :: this_box ! box index
  INTEGER, DIMENSION(1:nfragments(is)) :: frag_order
  REAL(DP), INTENT(IN) :: this_lambda ! fractional molecule?
  REAL(DP) :: randno
  REAL(DP) :: ln_pseq   ! probability of frag_order
  REAL(DP) :: ln_pbias  ! probability of placing each fragment in the box
  REAL(DP), INTENT(OUT) :: nrg_ring_frag_total ! potential energy of the
                                               ! isolated ring fragment
  LOGICAL, INTENT(INOUT) :: cbmc_overlap ! did all trials have core overlap?

  ! Local declarations
  INTEGER :: i, j, this_atom      ! atom indices
  INTEGER :: ifrag, this_fragment ! fragment indices
  INTEGER :: im                   ! molecule indices
  INTEGER :: total_frags  ! total number of conformations for this fragment in 
                          ! the reservoir
  INTEGER :: frag_start   ! random fragment to start growing from
  INTEGER :: frag_total   ! number of non-zero entries in frag_order
  INTEGER :: nl           ! number of the line where the x,y,x coords of the atom
                          ! of  the config and fragment randomly selected
                          ! were stored in the frag_position_library
  INTEGER, ALLOCATABLE, DIMENSION(:) :: live
  INTEGER, ALLOCATABLE, DIMENSION(:) :: frag_placed

  REAL(DP) :: x_anchor, y_anchor, z_anchor ! new COM coordinates for the 
                                           ! first fragment
  REAL(DP) :: xcom_old, ycom_old, zcom_old ! old COM coordinates for the
                                           ! first fragment
  REAL(DP) :: dx, dy, dz

  LOGICAL :: overlap ! TRUE if there is core overlap between a trial atom 
                     ! position and a an atom already in the box

  CHARACTER  :: this_file*120, symbol*1

  TYPE(Atom_Class), ALLOCATABLE, DIMENSION(:) :: config_list

  ! Variables associated with the CBMC part
  INTEGER :: itrial, trial, frag_type, n_frag_atoms

  REAL(DP) :: weight(kappa_ins), rand_no, E_dihed
  REAL(DP) :: E_intra_vdw, E_intra_qq, E_inter_vdw, E_inter_qq, E_total
  REAL(DP) :: nrg(kappa_ins), nrg_kBT, nrg_ring_frag

  LOGICAL :: del_overlap, overlap_trial(kappa_ins)

  Type(Atom_Class) :: rtrial(MAXVAL(natoms),0:MAX(kappa_ins,kappa_rot,kappa_dih))
  ! Slit pore variables

  LOGICAL :: framework_overlap
  REAL(DP) :: E_framework

!  ! DEBUGging variables
!  INTEGER :: M_XYZ_unit

  ! Initialize variables
  n_frag_atoms = 0
  ln_pbias = 0.0_DP
  this_box = molecule_list(this_im,is)%which_box ! which box this_im is in
  atom_list(:,this_im,is)%exist = .FALSE. ! mark all the atoms as deleted 
  molecule_list(this_im,is)%frac = 0.0_DP
  IF (ALLOCATED(frag_placed)) DEALLOCATE(frag_placed)
  ALLOCATE(frag_placed(nfragments(is)))
  frag_placed(:) = 0 ! =1 if fragment been placed
  nrg_ring_frag_total = 0.0_DP
  cbmc_flag = .TRUE.
  ! The energy of each trial coordinate will be stored in the array nrg.
  nrg(:) = 0.0_DP 
  overlap_trial(:) = .FALSE.


  !*****************************************************************************
  ! Step 1) Select which fragment will be inserted first
  !*****************************************************************************
  !
  ! One fragment will be inserted first, and then additional fragments will
  ! be added one at a time
  !

  ! get_fragorder is:
  !    *  .TRUE. when making an initial configuration, when inserting or 
  !       deleting a molecule in GCMC, when inserting a trial molecule to
  !       compute the chemical potential, or when transferring a molecule into 
  !       a box (GEMC insertion)
  !    *  .FALSE. when transfering a molecule out of a box (GEMC deletion), 
  !       since the frag_order was already selected as part of transferring the 
  !       molecule into the other box


  IF (get_fragorder) THEN
 
     ! Select a fragment to insert first
     randno = rranf()
     DO ifrag = 1, nfragments(is)
        IF (randno < frag_list(ifrag,is)%cum_prob_ins) EXIT
     END DO
     frag_start = ifrag
     ln_pseq = DLOG(frag_list(frag_start,is)%prob_ins) ! weighted prob of choosing ifrag
 
     ! Select the order fragments will be added to frag_start
     IF (ALLOCATED(live)) DEALLOCATE(live)
     ALLOCATE(live(nfragments(is)))      

     live(:) = 0
     frag_order(:) = 0
     frag_order(1) = frag_start ! this array will hold the order
     live(frag_start) = 1
     frag_total = 1
     
     ! If this molecule is made up of multiple fragments, select an order in 
     ! which the fragments will be grown
     IF (nfragments(is) > 1 ) THEN

        CALL Fragment_Order(frag_start,is,frag_total,frag_order,live,ln_pseq)
        
     END IF

     DEALLOCATE(live)

  ELSE

     ! frag_order is already determined, just need how many fragments will be
     ! placed (answer: all of them)
     frag_total = nfragments(is)

     ! ln_pseq was also computed as part of the GEMC insertion move

  END IF

  !
  ! At this point, we have the order in which we will grow the molecule.

  !*****************************************************************************
  ! Step 2) Choose a conformation for the first fragment
  !*****************************************************************************
  !   
  ! A fragment library was previously generated via a single fragment MC 
  ! simulation. Since the fragments were sampled according to the Boltzmann 
  ! distribution, we can now pull from the reservoir with uniform probability.
  ! The reservoir of fragment conformations is stored in frag_coords.
  !

  ! If get_fragorder = .FALSE., then frag_start is not yet defined
  frag_start = frag_order(1)

  ! If inserting a molecule, choose the first fragment's conformation from the 
  ! reservoir

  IF (.NOT. del_flag) THEN

     ! Pull from the reservoir with uniform probability
     !total_fragments is the number of configurations of the frag_start
     ! this_frag is the randomly configuration choosen
      total_frags = frag_list(frag_start,is)%nconfig
     this_fragment = INT(rranf() * total_frags) + 1

    
     frag_type = frag_list(frag_start,is)%type
    
     ! Read the coordinates for every atom
     ! this_fragment is the randomly config selected

     DO i = 1, frag_list(frag_start,is)%natoms 
        
        this_atom = frag_list(frag_start,is)%atoms(i)
        
        nl = (frag_position_library(frag_type)-1) + &
                                           frag_list(frag_start,is)%natoms*(this_fragment -1) + i 
        atom_list(this_atom,this_im,is)%rxp = &
                                    !  frag_coords(i,this_fragment,frag_type)%rxp
                                       library_coords(nl)%rxp
        atom_list(this_atom,this_im,is)%ryp = &
                                     ! frag_coords(i,this_fragment,frag_type)%ryp
                                        library_coords(nl)%ryp
        atom_list(this_atom,this_im,is)%rzp = &
                                     ! frag_coords(i,this_fragment,frag_type)%rzp
                                        library_coords(nl)%rzp
      END DO
  END IF
  ! Turn on the molecule and its individual atoms
  molecule_list(this_im,is)%frac = this_lambda
  
  DO i =1, frag_list(frag_start,is)%natoms
     this_atom = frag_list(frag_start,is)%atoms(i)
     atom_list(this_atom,this_im,is)%exist = .TRUE.
  END DO

  !*****************************************************************************
  ! Step 3) Rotate the first fragment
  !*****************************************************************************
  !
  ! At this time, only a single orientation is selected with uniform 
  ! probability. Future releases may attempt multiple orientations.
  ! 

  CALL Get_COM(this_im,is)
  CALL Compute_Max_COM_Distance(this_im,is)
  
  ! If inserting the fragment, select a random orientation 
  IF ( .NOT. del_flag) CALL Rotate_Molecule_Eulerian(this_im,is)

  !*****************************************************************************
  ! Step 4) Choose kappa_ins positions for the first fragment's COM
  !*****************************************************************************
  ! 
  ! The first fragment will now be inserted into the simulation box. Multiple
  ! trial coordinates will be randomly generated. Each trial will be weighted
  ! by the Boltzmann factor of the change in potential energy. One trial will be
  ! selected from the weighted distribution.

  ! Store the conformation and orientation as 0th trial
  DO i = 1, frag_list(frag_start,is)%natoms

     this_atom = frag_list(frag_start,is)%atoms(i)
     
     rtrial(this_atom,0)%rxp = atom_list(this_atom,this_im,is)%rxp
     rtrial(this_atom,0)%ryp = atom_list(this_atom,this_im,is)%ryp
     rtrial(this_atom,0)%rzp = atom_list(this_atom,this_im,is)%rzp

  END DO

  ! Store the COM
  xcom_old = molecule_list(this_im,is)%xcom
  ycom_old = molecule_list(this_im,is)%ycom
  zcom_old = molecule_list(this_im,is)%zcom

  ! We will place this fragment based only on its external weight

  ! When is imreplace greater than 0?
  IF(imreplace .GT. 0) THEN

     IF(.NOT. del_flag) THEN

        dx = molecule_list(imreplace,isreplace)%xcom &
           - molecule_list(this_im,is)%xcom
        dy = molecule_list(imreplace,isreplace)%ycom &
           - molecule_list(this_im,is)%ycom
        dz = molecule_list(imreplace,isreplace)%zcom &
           - molecule_list(this_im,is)%zcom

        molecule_list(this_im,is)%xcom = molecule_list(imreplace,isreplace)%xcom
        molecule_list(this_im,is)%ycom = molecule_list(imreplace,isreplace)%ycom
        molecule_list(this_im,is)%zcom = molecule_list(imreplace,isreplace)%zcom

        DO i = 1, frag_list(frag_start,is)%natoms

           this_atom = frag_list(frag_start,is)%atoms(i)
           atom_list(this_atom,this_im,is)%rxp = rtrial(this_atom,0)%rxp + dx
           atom_list(this_atom,this_im,is)%ryp = rtrial(this_atom,0)%ryp + dy
           atom_list(this_atom,this_im,is)%rzp = rtrial(this_atom,0)%rzp + dz
        END DO

     END IF

  ELSE

     ! Loop over the multiple trial coordinates
     DO itrial = 1, kappa_ins

        IF ( del_flag .AND. (itrial == 1 )) THEN
        
           ! Use the COM of the current position
           x_anchor = xcom_old
           y_anchor = ycom_old
           z_anchor = zcom_old
        
        ELSE

        
           ! Select a random trial coordinate
           IF (box_list(this_box)%int_box_shape == int_cubic) THEN
              x_anchor = (0.5_DP - rranf()) * box_list(this_box)%length(1,1)
              y_anchor = (0.5_DP - rranf()) * box_list(this_box)%length(2,2)
              z_anchor = (0.5_DP - rranf()) * box_list(this_box)%length(3,3)

           ELSE

              !Generate random positions in fractional coordinates

              x_anchor = 0.5_DP - rranf()
              y_anchor = 0.5_DP - rranf()
              z_anchor = 0.5_DP - rranf()

              !transform back to cartesian

              x_anchor = box_list(this_box)%length(1,1)*x_anchor + &
                       box_list(this_box)%length(1,2)*y_anchor +   &
                       box_list(this_box)%length(1,3)*z_anchor

              y_anchor = box_list(this_box)%length(2,1)*x_anchor + &
                       box_list(this_box)%length(2,2)*y_anchor +   &
                       box_list(this_box)%length(2,3)*z_anchor

              z_anchor = box_list(this_box)%length(3,1)*x_anchor + &
                       box_list(this_box)%length(3,2)*y_anchor +   &
                       box_list(this_box)%length(3,3)*z_anchor


           END IF

        END IF

        ! Place the fragment (and all its atoms) at the trial coordinate
        DO i = 1, frag_list(frag_start,is)%natoms
        
           this_atom = frag_list(frag_start,is)%atoms(i)
        
           atom_list(this_atom,this_im,is)%rxp = & 
                                   rtrial(this_atom,0)%rxp - xcom_old + x_anchor
           atom_list(this_atom,this_im,is)%ryp = & 
                                   rtrial(this_atom,0)%ryp - ycom_old + y_anchor
           atom_list(this_atom,this_im,is)%rzp = &
                                   rtrial(this_atom,0)%rzp - zcom_old + z_anchor
           
           rtrial(this_atom,itrial)%rxp = atom_list(this_atom,this_im,is)%rxp
           rtrial(this_atom,itrial)%ryp = atom_list(this_atom,this_im,is)%ryp
           rtrial(this_atom,itrial)%rzp = atom_list(this_atom,this_im,is)%rzp
        
        END DO

        molecule_list(this_im,is)%xcom = x_anchor
        molecule_list(this_im,is)%ycom = y_anchor
        molecule_list(this_im,is)%zcom = z_anchor

        ! Note that the COM position is always chosen inside the simulation box 
        ! so there is no need to call Fold_Molecule.
           
        ! For the sake of completeness we calculate the distance of the atom 
        ! furtherst from the COM
        CALL Compute_Max_COM_Distance(this_im,is)
           
        ! Calculate the intermolecular energy of the fragment. Note that
        ! cbmc_flag has been set to true so that the following call will compute
        ! interaction energy of the growing molecule within a small distance
        overlap = .FALSE.
        CALL Compute_Molecule_Nonbond_Inter_Energy(this_im,is,&
                E_inter_vdw,E_inter_qq,overlap)
        nrg(itrial) = nrg(itrial) + E_inter_vdw + E_inter_qq 

        IF (overlap) THEN
           ! atoms are too close, set the weight to zero
           weight(itrial) = 0.0_DP
           overlap_trial(itrial) = .TRUE.
        ELSE
           nrg_kBT = beta(this_box) * nrg(itrial)

           IF ( nrg_kBT >= max_kBT) THEN
              ! the energy is too high, set the weight to zero
              weight(itrial) = 0.0_DP
              overlap_trial(itrial) = .TRUE.
           ELSE
              weight(itrial) = DEXP(-nrg_kBT)
           END IF
        END IF

!        ! BEGIN DEBUGGING OUTPUT
!        ! Write out the fragment coordinates for each trial position
!        M_XYZ_unit = movie_xyz_unit + this_box
!        DO i = 1, frag_list(frag_start,is)%natoms
!           this_atom = frag_list(frag_start,is)%atoms(i)
!           WRITE(M_XYZ_unit,*) &
!              TRIM(nonbond_list(this_atom,is)%element) // &
!              TRIM(int_to_string(itrial)), & 
!              rtrial(this_atom,itrial)%rxp, &
!              rtrial(this_atom,itrial)%ryp, &
!              rtrial(this_atom,itrial)%rzp
!        END DO
!        IF (itrial==1) THEN
!           WRITE(*,'(2(A,X,I5,X))'), 'i_mcstep', i_mcstep, 'lm', this_im
!           WRITE(*,'(4(A12,X))') 'POS:trial', 'energy', 'weight', 'overlap'
!        END IF
!        WRITE(*,'(I12,X,E12.6,X,E12.6,X,L12)') itrial, beta(this_box)*nrg(itrial), weight(itrial), overlap_trial(itrial)
!        ! END DEBUGGING OUTPUT

        ! Store the cumulative weight of each trial
        IF (itrial > 1 ) weight(itrial) = weight(itrial-1) + weight(itrial)
     
     END DO


     ! Reject the move if all trials tripped overlap
     IF (ALL(overlap_trial)) THEN
        cbmc_overlap = .TRUE.
        cbmc_flag = .FALSE.
        RETURN
     END IF

     ! Select one of the trial coordinates

     IF (del_flag) THEN
        
        ! For a deletion move, we want the weight of the current position,
        ! which is stored in trial 1.
        trial = 1

     ELSE

        ! Choose one from Golden sampling for an insertion move
        rand_no = rranf() * weight(kappa_ins)
     
        DO i = 1, kappa_ins
           IF ( rand_no < weight(i)) EXIT
        END DO
     
        trial = i

        IF ( trial == kappa_ins + 1 ) THEN
           ! None of the trials were picked. Could be due to the fact that all 
           ! the trials had a very small cumulative weight

           cbmc_overlap = .TRUE.
           cbmc_flag = .FALSE.
           RETURN
        END IF

     END IF

     ! Compute the weight of the selected trial coordinate
     ln_pbias = ln_pbias - beta(this_box) * nrg(trial) - DLOG(weight(kappa_ins))
  
     ! This line is not used
     e_total = nrg(trial)

     ! We chose the ith trial coordinate for placement. Store the ith trial
     ! coordinates in the atom_list array. Note that for the deletion move,
     ! trial=1 has the current coordinates of the fragment, so the molecule
     ! is not moved.
  
     DO i = 1, frag_list(frag_start,is)%natoms
     
        this_atom = frag_list(frag_start,is)%atoms(i)
     
        atom_list(this_atom,this_im,is)%rxp = rtrial(this_atom,trial)%rxp
        atom_list(this_atom,this_im,is)%ryp = rtrial(this_atom,trial)%ryp
        atom_list(this_atom,this_im,is)%rzp = rtrial(this_atom,trial)%rzp
     
     END DO

     ! Compute COM of the fully grown molecule
     CALL Get_COM(this_im,is)

  END IF
  
  ! Mark this fragment is placed
  e_total = 0.0_DP
  frag_placed(frag_start) = 1

  ! We have our first segment placed in the system. 

  ! If this is a ring fragment,
  ! then we need intramolecular energies of the fragment used in biasing.
  ! For insertion, it is the energy in the gas phase. For deletion, it is 
  ! the energy of the old configuration. Note that this energy does not include
  ! angle bending energy as it cancels out from the acceptance rule.

  IF (frag_list(frag_start,is)%ring) THEN

     IF (del_flag) THEN
        ! compute the old energy
        CALL Compute_Ring_Fragment_Energy(frag_start,this_im,is,this_box,&
                nrg_ring_frag)
     ELSE
        nrg_ring_frag =  nrg_frag(frag_type)%this_config_energy(this_fragment) 
                        ! nrg_frag(this_fragment,frag_type)
     END IF
     nrg_ring_frag_total = nrg_ring_frag_total + nrg_ring_frag
     
  END IF

  !*****************************************************************************
  ! Step 5) Select kappa_dih orientations for each additional fragment
  !*****************************************************************************
  ! 
  ! Now we will place rest of the segments based on the initial fragment placed

  ! Why will the dihedral energy be anything other than zero?
  CALL Compute_Molecule_Dihedral_Energy(this_im,is,E_dihed)
  E_total = E_dihed

  ! If we've gotten this far, cbmc_overlap = FALSE
  del_overlap = .FALSE.

  ! The first fragment in frag_order has already been placed.
  ! We will call Fragment_Placement to place the remaining fragments in
  ! frag_order, starting with fragment frag_order(2)
  frag_start = 2
  CALL Fragment_Placement(this_box,this_im,is,frag_start,frag_total, &
                          frag_order,frag_placed,this_lambda,E_total,ln_pbias, &
                          nrg_ring_frag_total,cbmc_overlap,del_overlap)

  ! Note that cbmc_overlap may be TRUE and the cbmc_flag will be properly 
  ! assigned FALSE while exiting the code.

  ! Note that del_overlap may be true while computing the weight of
  ! an old configuration. So set cbmc_overlap to TRUE so that the 
  ! calling routine appropriately sets the coordinates of the atoms that were 
  ! not grown

  IF (del_overlap) cbmc_overlap = .TRUE.

  ! Mark cbmc_flag as FALSE so that intermolecular nonbonded interactions
  ! are properly computed for the molecule

  cbmc_flag = .FALSE.
  
END SUBROUTINE Build_Molecule

!******************************************************************************
SUBROUTINE Build_Rigid_Fragment(this_im,is,this_box,frag_order,this_lambda, &
              ln_pseq,ln_pbias,nrg_ring_frag_total,cbmc_overlap)
!******************************************************************************
!
! The subroutine is based on Build_Molecule subroutine and 
! places the first fragment of the molecule. The first atom
! of the fragement in inserted in a biased manner and then orientation of the 
! the molecule is generated in a biased manner as well. Currently it will not 
! work for reactions, I will add imreplace at a later point of time. 
!
! Written by Neeraj Rai on 06/11/2010.
!
!********************************************************************************

  USE Rotation_Routines
  USE IO_Utilities
  USE File_Names
  USE Energy_Routines

  INTEGER :: this_im, is, this_box, first_atom
  INTEGER, DIMENSION(1:nfragments(is)) :: frag_order
  INTEGER :: nl, nlo      ! number of the line where the x,y,x coords of
                          ! config and fragment randomly selected
                          !were stored in the frag_position_library
  REAL(DP), INTENT(IN) :: this_lambda
  REAL(DP) :: ln_pseq
  REAL(DP) :: ln_pbias
  REAL(DP), INTENT(OUT) :: nrg_ring_frag_total
  LOGICAL, INTENT(INOUT) :: cbmc_overlap
  !-------------------------------

  INTEGER :: frag_start, frag_start_old, frag_total, i, this_atom
  INTEGER :: anchor, ifrag, total_connect, frag_connect
  INTEGER :: j, im
  INTEGER :: ii,jj,kk, total_frags, this_fragment
  INTEGER, ALLOCATABLE, DIMENSION(:) :: live, deadend, central
  INTEGER, ALLOCATABLE, DIMENSION(:) :: frag_placed
  REAL(DP) :: x_anchor, y_anchor, z_anchor, xcom_old, ycom_old, zcom_old
  REAL(DP) :: dx, dy, dz
  REAL(DP) :: x_this, y_this, z_this,tempx, tempy, tempz, temp_var, E_ang
  LOGICAL :: overlap
  CHARACTER  :: this_file*120, symbol*1
  TYPE(Atom_Class), ALLOCATABLE, DIMENSION(:) :: config_list
  ! Variables associated with the CBMC part

  INTEGER :: itrial, trial, frag_type
  REAL(DP) :: weight(MAX(kappa_ins,kappa_rot,kappa_dih)), rand_no
  REAL(DP) :: e_dihed, E_intra_vdw, E_intra_qq, E_inter_vdw
  REAL(DP) :: E_inter_qq, E_total
  REAL(DP) :: nrg(MAX(kappa_ins,kappa_rot,kappa_dih)), nrg_kBT, time0, time1, nrg_ring_frag
  LOGICAL :: del_overlap
  TYPE(Atom_Class) :: rtrial(MAXVAL(natoms),0:MAX(kappa_ins,kappa_rot,kappa_dih))



  weight(:)=0.0_DP
  cbmc_flag = .TRUE.
  ln_pbias = 0.0_DP
  
  ! Assign a locate number for the molecule
  this_box = molecule_list(this_im,is)%which_box
  ! mark all the atoms as deleted 
  atom_list(:,this_im,is)%exist = .FALSE.

  IF (ALLOCATED(frag_placed)) DEALLOCATE(frag_placed)
  ALLOCATE(frag_placed(nfragments(is)))
  frag_placed(:) = 0
  nrg_ring_frag_total = 0.0_DP

! Get the order of the fragment growth
  IF (get_fragorder) THEN
     frag_start = 1
     ln_pseq = 0.0_DP
     ! Obtain the order of fragment addition
     IF (ALLOCATED(live)) DEALLOCATE(live)
     ALLOCATE(live(nfragments(is)))      
     live(:) = 0
     frag_order(:) = 0
     frag_order(1) = frag_start
     live(frag_start) = 1
     frag_total = 1
     IF (nfragments(is) > 1 ) THEN
        CALL Fragment_Order(frag_start,is,frag_total,frag_order,live,ln_pseq)
     END IF
     DEALLOCATE(live)
  ELSE
     frag_total = nfragments(is)
  END IF

  ! At this point, we have the order in which we will grow the molecule from
  ! Add the part to read in the coordinates from its file.
  ! we will make all the atoms of frag_start as part of the simulations

  atom_list(:,this_im,is)%exist = .FALSE.
  molecule_list(this_im,is)%frac = 0.0_DP

  frag_start = frag_order(1)


  ! Note that we need to choose from the reservoir only when insertion
  ! is attempted

  molecule_list(this_im,is)%frac = this_lambda

  ! NR: Select the location of first bead of the fragment

  first_atom = frag_list(frag_start,is)%atoms(1)
  atom_list(first_atom,this_im,is)%exist = .true.


  DO itrial = 1, kappa_ins
     IF (del_flag .and.(itrial.eq.1)) THEN
        rtrial(first_atom,itrial)%rxp =  atom_list(first_atom,this_im,is)%rxp
        rtrial(first_atom,itrial)%ryp =  atom_list(first_atom,this_im,is)%ryp
        rtrial(first_atom,itrial)%rzp =  atom_list(first_atom,this_im,is)%rzp
     ELSE
        IF (box_list(this_box)%int_box_shape == int_cubic) THEN
           
           rtrial(first_atom,itrial)%rxp = (0.5_DP - rranf()) * &
                box_list(this_box)%length(1,1)
           rtrial(first_atom,itrial)%ryp = (0.5_DP - rranf()) * &
                box_list(this_box)%length(2,2)
           rtrial(first_atom,itrial)%rzp = (0.5_DP - rranf()) * &
                box_list(this_box)%length(3,3)
           
        END IF

     ENDIF

  ! Store them in correct array for the energy call
     atom_list(first_atom,this_im,is)%rxp = rtrial(first_atom,itrial)%rxp
     atom_list(first_atom,this_im,is)%ryp = rtrial(first_atom,itrial)%ryp
     atom_list(first_atom,this_im,is)%rzp = rtrial(first_atom,itrial)%rzp


     ! Compute the energy of this atom    
     molecule_list(this_im,is)%xcom = rtrial(first_atom,itrial)%rxp 
     molecule_list(this_im,is)%ycom = rtrial(first_atom,itrial)%ryp 
     molecule_list(this_im,is)%zcom = rtrial(first_atom,itrial)%rzp

     CALL Compute_Max_COM_Distance(this_im,is)

     CALL Compute_Molecule_Nonbond_Inter_Energy(this_im,is,E_inter_vdw,E_inter_qq,overlap)


!    WRITE(8,*) del_flag,this_box, itrial, this_im,E_inter_vdw,  E_inter_qq,overlap

  ! compute weight for itrial 

     IF (overlap) THEN
         ! the energy is too high, set the weight to zero
         weight(itrial) = 0.0_DP
     ELSE
         nrg_kBT = beta(this_box) * (E_inter_vdw + E_inter_qq )

!         nrg_kBT = beta(this_box) * (E_inter_vdw )
!!$        ! compute the weight of this trial for the reverse move and first
!!$        ! trial. The following IF-ELSE construct ensures that the weight
!!$        ! for the reverse move is always computed.
         IF ( del_flag .AND. (itrial == 1) ) THEN
            weight(itrial) = DEXP(-nrg_kBT)
         ELSE IF ( nrg_kBT >= max_kBT) THEN
            weight(itrial) = 0.0_DP
         ELSE
            weight(itrial) = DEXP(-nrg_kBT)
         END IF
      END IF
!     WRITE(8,*) del_flag, itrial, weight(itrial),nrg_kBT,beta(this_box)
      IF (itrial > 1 ) weight(itrial) = weight(itrial-1) + weight(itrial)
      
   ENDDO ! End the generation of trial locations and associating weight

  IF (weight(kappa_ins) == 0.0_DP) THEN
      cbmc_overlap = .TRUE.
      cbmc_flag = .FALSE.
      RETURN
  END IF

  IF (del_flag) THEN
      trial = 1
  ELSE
   ! Choose one from Golden sampling for an insertion move
      rand_no = rranf() * weight(kappa_ins)
      DO i = 1, kappa_ins
         IF ( rand_no < weight(i)) EXIT
      END DO
      trial = i
  END IF
!  Write(8,*) 'selected', trial

  ln_pbias = ln_pbias - beta(this_box) * nrg(trial) - DLOG(weight(kappa_ins))

! Now we have the position of the first bead of   
! Assign this position to proper atom_list
  atom_list(first_atom,this_im,is)%rxp = rtrial(first_atom,trial)%rxp
  atom_list(first_atom,this_im,is)%ryp = rtrial(first_atom,trial)%ryp
  atom_list(first_atom,this_im,is)%rzp = rtrial(first_atom,trial)%rzp

! NR: End selection of the first bead of the fragment


! Now set the atoms of the first fragment to true
  DO i =1, frag_list(frag_start,is)%natoms
     this_atom = frag_list(frag_start,is)%atoms(i)
     atom_list(this_atom,this_im,is)%exist = .TRUE.
  END DO

! if it is insertion move, get the fragment from the reservoir

  IF (.NOT. del_flag) THEN
     ! obtain a random configuration
     total_frags = frag_list(frag_start,is)%nconfig
     ! Choose a fragment at random
     this_fragment = INT(rranf() * total_frags) + 1
     frag_type = frag_list(frag_start,is)%type
     DO i = 1, frag_list(frag_start,is)%natoms 
        this_atom = frag_list(frag_start,is)%atoms(i)
        nlo = (frag_position_library(frag_type)-1) + &
                                            frag_list(1,is)%natoms*(this_fragment-1) + 1
        nl = (frag_position_library(frag_type)-1) + &
                                           frag_list(1,is)%natoms*(this_fragment-1)+i

                    
         rtrial(this_atom,0)%rxp = library_coords(nl)%rxp -&
                                   library_coords(nlo)%rxp +&
                                   atom_list(first_atom,this_im,is)%rxp
         rtrial(this_atom,0)%ryp = library_coords(nl)%ryp -&
                                   library_coords(nlo)%ryp +& 
                                   atom_list(first_atom,this_im,is)%ryp
         rtrial(this_atom,0)%rzp = library_coords(nl)%rzp -& 
                                   library_coords(nlo)%rzp +& 
                                   atom_list(first_atom,this_im,is)%rzp
      !  rtrial(this_atom,0)%rxp = frag_coords(i,this_fragment,frag_type)%rxp-&
      !                            frag_coords(1,this_fragment,frag_type)%rxp+&
      !                            atom_list(first_atom,this_im,is)%rxp 
      !  rtrial(this_atom,0)%ryp = frag_coords(i,this_fragment,frag_type)%ryp-&
      !                            frag_coords(1,this_fragment,frag_type)%ryp+&
      !                            atom_list(first_atom,this_im,is)%ryp
      !  rtrial(this_atom,0)%rzp = frag_coords(i,this_fragment,frag_type)%rzp-&
      !                            frag_coords(1,this_fragment,frag_type)%rzp+&
      !                            atom_list(first_atom,this_im,is)%rzp 
     END DO
  ELSE
     DO i = 1,frag_list(frag_start,is)%natoms
        this_atom = frag_list(frag_start,is)%atoms(i)
        rtrial(this_atom,0)%rxp = atom_list(this_atom,this_im,is)%rxp
        rtrial(this_atom,0)%ryp = atom_list(this_atom,this_im,is)%ryp
        rtrial(this_atom,0)%rzp = atom_list(this_atom,this_im,is)%rzp
     END DO
  END IF

! Now fragment's first bead at the new location

! Start the rotational-bias
! Rotate this fragment about the first bead along different axes
! select one from trials

  DO itrial = 1, kappa_rot
     IF(del_flag .and. (itrial.eq.1)) THEN
       DO i=1,frag_list(frag_start,is)%natoms
          this_atom = frag_list(frag_start,is)%atoms(i)
          rtrial(this_atom,itrial)%rxp = atom_list(this_atom,this_im,is)%rxp
          rtrial(this_atom,itrial)%ryp = atom_list(this_atom,this_im,is)%ryp
          rtrial(this_atom,itrial)%rzp = atom_list(this_atom,this_im,is)%rzp
       END DO
       CALL Get_COM(this_im,is)
    ELSE
       DO i = 1,frag_list(frag_start,is)%natoms
          this_atom = frag_list(frag_start,is)%atoms(i)  
          atom_list(this_atom,this_im,is)%rxp = rtrial(this_atom,0)%rxp
          atom_list(this_atom,this_im,is)%ryp = rtrial(this_atom,0)%ryp
          atom_list(this_atom,this_im,is)%rzp = rtrial(this_atom,0)%rzp
       ENDDO
       CALL Rotate_XYZ_Axes(this_im,is,frag_start,.true.,.true.,.true.,1)
       CALL Get_COM(this_im,is)
       !       CALL Fold_Molecule(this_im,is,this_box) 
       DO i = 1,frag_list(frag_start,is)%natoms
          this_atom = frag_list(frag_start,is)%atoms(i)
          rtrial(this_atom,itrial)%rxp = atom_list(this_atom,this_im,is)%rxp
          rtrial(this_atom,itrial)%ryp = atom_list(this_atom,this_im,is)%ryp
          rtrial(this_atom,itrial)%rzp = atom_list(this_atom,this_im,is)%rzp
       END DO
       
    ENDIF
    !    Get COM
    CALL Compute_Max_COM_Distance(this_im,is)
    

!     Write(8,*) frag_list(frag_start,is)%natoms
!     Write(8,*)
!     DO i = 1,frag_list(frag_start,is)%natoms
!           this_atom = frag_list(frag_start,is)%atoms(i)
!           IF (this_atom .eq. 1) then 
!               Write(8,*) 'O', atom_list(this_atom,this_im,is)%rxp, &
!                               atom_list(this_atom,this_im,is)%ryp, &
!                               atom_list(this_atom,this_im,is)%rzp
!           ELSE
!               Write(8,*) 'C', atom_list(this_atom,this_im,is)%rxp, &
!                               atom_list(this_atom,this_im,is)%ryp, &
!                               atom_list(this_atom,this_im,is)%rzp
!           END IF
!     END DO  

!     atom_list(first_atom,this_im,is)%exist = .false.



    overlap = .FALSE.
       
       CALL Compute_Molecule_Nonbond_Inter_Energy(this_im,is,E_inter_vdw,E_inter_qq,overlap)


!     atom_list(first_atom,this_im,is)%exist = .true.

!     Write(8,*) 'First fragment attempt vdw qq',itrial, E_inter_vdw, E_inter_qq,this_box,del_flag,overlap

 ! compute weight for itrial 
     IF (overlap) THEN
         ! the energy is too high, set the weight to zero
         weight(itrial) = 0.0_DP
     ELSE
         nrg_kBT = beta(this_box) * (E_inter_vdw + E_inter_qq )
!         nrg_kBT = beta(this_box) * (E_inter_vdw )
!!$        ! compute the weight of this trial for the reverse move and first
!!$        ! trial. The following IF-ELSE construct ensures that the weight
!!$        ! for the reverse move is always computed.
         IF ( del_flag .AND. (itrial == 1) ) THEN
            weight(itrial) = DEXP(-nrg_kBT)
         ELSE IF ( nrg_kBT >= max_kBT) THEN
            weight(itrial) = 0.0_DP
         ELSE
            weight(itrial) = DEXP(-nrg_kBT)
         END IF
     END IF
     IF (itrial > 1 ) weight(itrial) = weight(itrial-1) + weight(itrial)
  ENDDO ! End the generation of trial orientations and associating waight

  IF (weight(kappa_rot) == 0.0_DP) THEN
      cbmc_overlap = .TRUE.
      cbmc_flag = .FALSE.
      RETURN
  END IF

  IF (del_flag) THEN
      trial = 1
  ELSE
   ! Choose one from Golden sampling for an insertion move
      rand_no = rranf() * weight(kappa_rot)
      DO i = 1, kappa_rot
         IF ( rand_no < weight(i)) EXIT
      END DO
      trial = i
  END IF

!  Write(8,*) 'selected', trial

  ln_pbias = ln_pbias - beta(this_box) * nrg(trial) - DLOG(weight(kappa_rot))

! Assign positions to the atom_list
  DO i=1,frag_list(frag_start,is)%natoms
     this_atom =  frag_list(frag_start,is)%atoms(i)
     atom_list(this_atom,this_im,is)%rxp = rtrial(this_atom,trial)%rxp
     atom_list(this_atom,this_im,is)%ryp = rtrial(this_atom,trial)%ryp
     atom_list(this_atom,this_im,is)%rzp = rtrial(this_atom,trial)%rzp
  END DO

! IF ( .NOT. del_flag) THEN
! write(*,*) atom_list(:,this_im,is)%rxp
! write(*,*) atom_list(:,this_im,is)%rzp
! write(*,*) atom_list(:,this_im,is)%ryp
! NR: Now we have place the first fragment.
! write(*,*)
!	End if
  frag_placed(frag_start) = 1


  IF (frag_list(frag_start,is)%ring) THEN

     IF (del_flag) THEN
        ! compute the old energy
        CALL Compute_Ring_Fragment_Energy(frag_start,this_im,is,this_box,nrg_ring_frag)
     ELSE
        nrg_ring_frag =   nrg_frag(frag_type)%this_config_energy(this_fragment) 
                         !nrg_frag(this_fragment,frag_type)
     END IF
     nrg_ring_frag_total = nrg_ring_frag_total + nrg_ring_frag

  END IF

  ! Now we will place rest of the segments based on the initial fragment placed

  CALL Compute_Molecule_Dihedral_Energy(this_im,is,e_dihed)
  e_total = e_dihed

  CALL Fragment_Placement(this_box,this_im,is,2,frag_total,frag_order,frag_placed,this_lambda, &
       e_total,ln_pbias,nrg_ring_frag_total, cbmc_overlap, del_overlap)

  cbmc_flag = .FALSE.

END SUBROUTINE Build_Rigid_Fragment

!*******************************************************************************
SUBROUTINE Cut_Regrow(this_im,is,frag_live,frag_dead,frag_order,frag_total, &
     this_lambda,E_prev,ln_pseq,ln_pbias,nrg_ring_frag_tot,cbmc_overlap,del_overlap)
  !*****************************************************************************
  ! The subroutine cuts a bond, deletes part of a molecule and regrows using 
  ! configurational biasing. 
  !
  ! Written by Jindal Shah on 08/26/08
  !
  ! Step 1) Select a bond to cut with uniform probability & 
  !         delete fragments on one side of the bond
  ! Step 2) Compute energy
  ! Step 3) Regrow deleted fragments
  !
  !*****************************************************************************

  USE Global_Variables
  USE Random_Generators

  IMPLICIT NONE

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  ! Arguments
  INTEGER :: this_im    ! molecule index
  INTEGER :: is         ! species index
  INTEGER :: frag_live  ! the fragment with the cut bond that will be kept
  INTEGER :: frag_dead  ! the fragment with the cut bond that will be deleted
  INTEGER :: frag_order(1:nfragments(is)) ! order in which frags are added
  INTEGER :: frag_total ! number of entries in frag_order
  REAL(DP), INTENT(IN) :: this_lambda
  REAL(DP) :: E_prev
  REAL(DP) :: ln_pseq   ! probability of sequence in frag_order
  REAL(DP) :: ln_pbias  ! probability of placing each fragment in the box
  REAL(DP) :: nrg_ring_frag_tot
  LOGICAL :: cbmc_overlap
  LOGICAL :: del_overlap

  ! Local declarations
  INTEGER :: frag_bond, frag_start
  INTEGER :: this_box, i_frag, frag_no, i, this_atom, anchor, this_frag
  INTEGER :: anchor_live, anchor_dead
  
  INTEGER, DIMENSION(:), ALLOCATABLE :: live ! fragment was not deleted OR 
                                             ! is in frag_order
  INTEGER, DIMENSION(:), ALLOCATABLE :: frag_placed ! fragment was not deleted

  REAL(DP) :: E_total, E_intra_vdw, E_intra_qq, E_inter_vdw, E_inter_qq, E_dihed

  LOGICAL :: overlap

  ! Initialize variables
  cbmc_flag = .TRUE.
  del_overlap = .FALSE.
  nrg_ring_frag_tot = 0.0_DP
  ln_pseq = 0.0_DP
  ln_pbias = 0.0_DP

  !*****************************************************************************
  ! Step 1) Select a bond to cut with uniform probability & 
  !         delete fragments on one side of the bond
  !*****************************************************************************
  ! 
  ! obtain the box in which molecule is selected

  this_box = molecule_list(this_im,is)%which_box

  ! Choose a fragment bond to cut

  IF(ALLOCATED(frag_placed)) DEALLOCATE(frag_placed)
  ALLOCATE(frag_placed(nfragments(is)))

  ! When is del_flag == TRUE?
  IF ( .NOT. DEL_FLAG) THEN
     
     frag_bond = INT (rranf() * fragment_bonds(is)) + 1
     
     ! Kill one of the fragments
     
     IF ( rranf() < fragment_bond_list(frag_bond,is)%prob_del1) THEN
        
        frag_dead = fragment_bond_list(frag_bond,is)%fragment1
        frag_live = fragment_bond_list(frag_bond,is)%fragment2
        
     ELSE
        
        frag_dead = fragment_bond_list(frag_bond,is)%fragment2
        frag_live = fragment_bond_list(frag_bond,is)%fragment1
        
     END IF
     
     ! Technically, we should mark all fragments connected to frag_live as live
     ! (except for frag_dead) and then all fragments connected to those
     ! fragments as live, until we come to a terminal fragment. However,  
     ! Fragment_Order will only add dead fragments connected to frag_dead to 
     ! frag_order, so we can get away here with only marking frag_live as live
     ALLOCATE(live(nfragments(is)))
     live(:) = 0
     live(frag_live) = 1

     ! frag_order starts with frag_dead
     frag_order(:) = 0
     frag_order(1) = frag_dead
     frag_total = 1
     live(frag_dead) = 1 ! frag_dead is now live b/c it is in frag_order
     
     ! Obtain random order of fragments to be regrown
     CALL Fragment_Order(frag_dead,is,frag_total,frag_order,live,ln_pseq)

     DEALLOCATE(live)

  END IF
     
  ! Now we need all the fragments that have been placed. To do this, mark all 
  ! the fragments as placed and unmark the ones from frag_order.
  atom_list(:,this_im,is)%exist = .TRUE.
  frag_placed(:) = 1
  
  DO i_frag = 1, frag_total

     this_frag = frag_order(i_frag)
     
     frag_placed(this_frag) = 0
     
     ! also delete all the atoms contained in this fragment
     
     DO i = 1, frag_list(this_frag,is)%natoms
        
        this_atom = frag_list(this_frag,is)%atoms(i)

        atom_list(this_atom,this_im,is)%exist = .FALSE.
        
     END DO

  END DO
 
  ! Note that the anchor of frag_dead is present in the simulation so turn
  ! the exist flag true for this atom. In the process, we have also marked
  ! the anchor of frag_live as dead, so we will make that atom alive as well

  ! Get the two anchor atoms that need to be turned alive

  CALL Get_Common_Fragment_Atoms(is,frag_live,frag_dead,anchor_live,anchor_dead)

  atom_list(anchor_live,this_im,is)%exist = .TRUE.
  atom_list(anchor_dead,this_im,is)%exist = .TRUE.
 
  !*****************************************************************************
  ! Step 2) Compute energy
  !*****************************************************************************

  ! At this point we will calculate the energy of the fragment(s) in the 
  ! simulation. This is done so that the energy calculated for bonded and 
  ! nonbonded interactions do not cause trouble while growing the molecule. 
  ! Since dihedral, intramolecular nonbond and intermolecular nonbond energies 
  ! are used to bias the fragment placement, we will compute these energies. 
  ! Note that this is computed only when del_flag is false. We will use the 
  ! energy computed here when del_flag is true.

  IF (.NOT. del_flag) THEN

     CALL Compute_Molecule_Dihedral_Energy(this_im,is,E_dihed)
     E_prev = E_dihed

  END IF

  ! E_total is equal to the energy computed above during the growth phase.
  ! while it is the value obtained for cut_N_grow for calculating weight of the 
  ! old chain
  E_total = E_prev 

  !*****************************************************************************
  ! Step 3) Regrow deleted fragments
  !*****************************************************************************
  !
  ! We will call Fragment_Placement to place the fragments in frag_order, 
  ! starting with fragment frag_order(1)

  frag_start = 1
  CALL Fragment_Placement(this_box,this_im,is,frag_start,frag_total, &
       frag_order,frag_placed,this_lambda,E_total,ln_pbias,nrg_ring_frag_tot, &
       cbmc_overlap,del_overlap)
  cbmc_flag = .FALSE.
   

END SUBROUTINE Cut_Regrow



!*******************************************************************************
SUBROUTINE Fragment_Order(this_frag,is,frag_total,frag_order,live,ln_pseq)
!*******************************************************************************
!
! PURPOSE: Determine the order in which remaining fragments will be placed
!
! Written by Ryan Mullen on 05/08/14
!
! DESCRIPTION: This Subroutine determines the order in which the fragments will
! be (re)grown to form a completed molecule. 
!   * If the molecule is being inserted, this_frag is the only live fragment.
!     this_frag may be connected to only 1 other fragment (if this_frag is a 
!     terminal fragment) or more than 1 fragment (if this_frag is an internal
!     fragment).
!   * If the molecule is being partially regrown, a bond between two fragments
!     has been cut and part of the molecule deleted. this_frag is the only
!     deleted fragment that is connected to a live fragment. Whether this_frag
!     is terminal or internal, there is only 1 hanging connection initially in 
!     the molecule.
! 
!   Step 1) Determine the number & identity of hanging connections
!   Step 2) Select which fragment to add next
!   Step 3) Update the number & identity of hanging connections
!   
!   Loop thru steps 2 and 3 until all fragments are live.
!
!*******************************************************************************
  USE Global_Variables
  USE Random_Generators
  
  IMPLICIT NONE

  !*****************************************************************************
  ! Declare & Initialize Variables
  !*****************************************************************************

  ! Arguments
  INTEGER :: this_frag  ! fragment index
  INTEGER :: is         ! species index
  INTEGER :: frag_total ! running total of how many fragments are live
  INTEGER, DIMENSION(1:nfragments(is)) :: frag_order ! the order in which 
                                                     ! fragments will be grown
  INTEGER, DIMENSION(1:nfragments(is)) :: live       ! fragment in frag_order?
  REAL(DP) :: ln_pseq   ! the probability of frag_order

  ! Local declarations
  INTEGER :: ifrag, jfrag, ifrag_id, jfrag_id, i
  INTEGER :: nconnect ! number of entires in hanging_connections
  INTEGER :: hanging_connections(nfragments(is)) ! to hold ids of frags ready to add
  REAL(DP) :: randno, prob(nfragments(is)), cum_prob(nfragments(is))

  !*****************************************************************************
  !   Step 1) Determine the number & identity of hanging connections
  !*****************************************************************************

  ! If the molecule has more than 1 fragment, each fragment is connected to
  ! others. We need to know how many fragments connected to this_frag are not 
  ! yet live, and which fragments those are:

  nconnect = 0 ! to store the number of hanging connections
  hanging_connections(:) = 0 ! to store the frag ids that are ready to be added
  DO ifrag = 1, frag_list(this_frag,is)%nconnect
     ifrag_id = frag_list(this_frag,is)%frag_connect(ifrag)
     IF ( live(ifrag_id) == 0 ) THEN
        nconnect = nconnect + 1
        hanging_connections(nconnect) = ifrag_id
     END IF
  END DO

  ! Loop through steps 2 and 3 until all fragments are live.

  DO WHILE (nconnect > 0)
     !**************************************************************************
     !   Step 2) Select which fragment to add next
     !**************************************************************************
     ! Compute probabilty of adding next fragment for each frag in hanging_connections
     ! weight each fragment by the number of unique atoms in the fragment (natoms - 2)
     ! On the first loop, store the number of atoms in prob & cum_prob
     prob = 0
     cum_prob = 0
     DO ifrag = 1, nconnect
        ifrag_id = hanging_connections(ifrag)
        prob(ifrag) = frag_list(ifrag_id,is)%natoms - 2
        IF (ifrag == 1) THEN
          cum_prob(ifrag) = frag_list(ifrag_id,is)%natoms - 2
        ELSE
          cum_prob(ifrag) = cum_prob(ifrag-1)  + (frag_list(ifrag_id,is)%natoms - 2)
        END IF
     END DO
     ! On the second loop, divide the prob & cum_prob for each frag by the total number of atoms
     DO ifrag = 1, nconnect
        ifrag_id = hanging_connections(ifrag)
        prob(ifrag) = prob(ifrag) / cum_prob(nconnect)
        cum_prob(ifrag) = cum_prob(ifrag) / cum_prob(nconnect)
     END DO

     ! Choose a random fragment
     randno = rranf()
     DO ifrag = 1, nconnect
        IF (randno < cum_prob(ifrag)) EXIT
     END DO
     ifrag_id = hanging_connections(ifrag)
     ln_pseq = ln_pseq + DLOG( prob(ifrag) )

     ! Add ifrag_id to frag_order, make ifrag_id live
     frag_total = frag_total + 1
     frag_order(frag_total) = ifrag_id
     live(ifrag_id) = 1

     !**************************************************************************
     !   Step 3) Update the number & identity of hanging connections
     !**************************************************************************
     ! Remove ifrag_id from hanging_connections
     DO jfrag = ifrag, nconnect - 1
        hanging_connections(jfrag) = hanging_connections(jfrag + 1)
     END DO
     hanging_connections(nconnect) = 0
     nconnect = nconnect - 1

     ! Update nconnect and hanging_connections
     DO jfrag = 1, frag_list(ifrag_id,is)%nconnect
        jfrag_id = frag_list(ifrag_id,is)%frag_connect(jfrag)
        IF ( live(jfrag_id) == 0 ) THEN
           nconnect = nconnect + 1
           hanging_connections(nconnect) = jfrag_id
        END IF
     END DO

  END DO

END SUBROUTINE Fragment_Order

!*******************************************************************************
SUBROUTINE Fragment_Placement(this_box, this_im, is, frag_start, frag_total, &
           frag_order, frag_placed, this_lambda, e_total, ln_pbias, &
           nrg_ring_frag_tot, cbmc_overlap, del_overlap)
!*******************************************************************************
!
! PURPOSE: place the remaining fragments of the molecule
!
! DESCRIPTION: This subroutine performs the following steps
!
! Step 1) Select a fragment conformation
! Step 2) Align the fragment to the growing molecule
! Step 3) Select kappa_dih orientations for each additional fragment
! Step 4) Compute the energy of the fragment 
! Step 5) Select a trial dihedral using the weighted probabilities
!
!*******************************************************************************

  USE IO_Utilities

  !*****************************************************************************
  ! Declare and Initialize Variables
  !*****************************************************************************

  INTEGER, INTENT(IN) :: this_im, is, frag_start, frag_total, this_box
  INTEGER, DIMENSION(1:nfragments(is)), INTENT(IN) :: frag_order

  INTEGER, DIMENSION(1:nfragments(is)), INTENT(INOUT) :: frag_placed

  REAL(DP), INTENT(INOUT) :: e_total, ln_pbias, nrg_ring_frag_tot
  REAL(DP), INTENT(IN) :: this_lambda
  
  LOGICAL, INTENT(INOUT) :: cbmc_overlap, del_overlap

  ! local variable

  INTEGER :: i, ifrag, j, this_atom, frag_connect, total_connect, anchor_ifrag, is1, im, ia
  INTEGER :: anchor_frag_connect, atom_ifrag, atom_frag_connect, ii, trial
  INTEGER :: frag_type, dumcount

  INTEGER, DIMENSION(:), ALLOCATABLE :: counted, connection, atom_id

  INTEGER :: ispecies, jmol, k

  INTEGER :: nl, nlo      ! number of the line where start the x,y,x coords of
                          ! config and fragment randomly selected

  INTEGER :: total_frags, this_fragment, nfrag_atoms
  REAL(DP) :: x_this,y_this,z_this, vec1(3), vec2(3), aligner_ifrag(3,3)
  REAL(DP) :: hanger_ifrag(3,3), aligner_frag_connect(3,3), hanger_frag_connect(3,3)

  REAL(DP) :: tempx, tempy, tempz, theta, e_dihed
  REAL(DP) :: weight(kappa_dih), nrg(kappa_dih)
  REAL(DP) :: E_intra_qq, E_intra_vdw, prob_pick
  REAL(DP) :: e_prev, temp_var, E_ang, E_inter_vdw, E_inter_qq
  REAL(DP) :: nrg_kBT, p_acc, nrg_intra_vdw, nrg_intra_qq, nrg_inter_vdw, nrg_inter_qq
  REAL(DP) :: trial_weight
  REAL(DP) :: nrg_ring_frag, nrg_dihed(MAX(kappa_ins,kappa_rot,kappa_dih))

  LOGICAL :: overlap, overlap_trial(kappa_dih)

  CHARACTER :: this_file*120, element*1 

  TYPE(Atom_Class), ALLOCATABLE, DIMENSION(:) :: config_list
  TYPE(Atom_Class), ALLOCATABLE, DIMENSION(:,:) :: config_temp_list

!  ! DEBUGging variables
!  INTEGER :: M_XYZ_unit


  IF (.NOT. ALLOCATED(counted)) ALLOCATE(counted(nfragments(is)))
  IF (.NOT. ALLOCATED(connection)) ALLOCATE(connection(nfragments(is)))
  ALLOCATE(config_list(natoms(is)))
! Caution: NR confirm with Jindal. This could just be kappa_dih

  ALLOCATE(config_temp_list(natoms(is),MAX(kappa_ins,kappa_rot,kappa_dih)))
  ALLOCATE(atom_id(natoms(is)))

  config_list(:)%rxp = 0.0_DP
  config_list(:)%ryp = 0.0_DP
  config_list(:)%rzp = 0.0_DP
  dumcount = 0
  e_prev = e_total

  DO i = frag_start, frag_total

     ifrag = frag_order(i)
     
     !**************************************************************************
     ! Step 1) Select a fragment conformation
     !**************************************************************************
     !
      IF ( del_flag) THEN

        ! For a deletion move, use the existing conformation for each additional
        ! fragment

        DO j = 1, frag_list(ifrag,is)%natoms
           this_atom = frag_list(ifrag,is)%atoms(j)

           atom_list(this_atom,this_im,is)%exist = .TRUE.

           config_list(this_atom)%rxp = atom_list(this_atom,this_im,is)%rxp
           config_list(this_atom)%ryp = atom_list(this_atom,this_im,is)%ryp
           config_list(this_atom)%rzp = atom_list(this_atom,this_im,is)%rzp
        END DO

        ! For a ring fragment, calculate the fragment intramolecular energy

        IF (frag_list(ifrag,is)%ring) THEN
           CALL Compute_Ring_Fragment_Energy(ifrag,this_im,is,this_box, &
                   nrg_ring_frag)
           nrg_ring_frag_tot = nrg_ring_frag_tot + nrg_ring_frag
        END IF

     ELSE

        ! Select a fragment conformation from the reservoir
        ! The reservoir was populated with a Boltzmann distribution, so now we
        ! can pull from it with a uniform probability

        total_frags = frag_list(ifrag,is)%nconfig
        this_fragment = INT(rranf() * total_frags) + 1
        
        ! Select a fragment conformation from the reservoir
        ! The reservoir was populated with a Boltzmann distribution, so now we
        ! can pull from it with a uniform probability
        
        this_fragment = INT(rranf() * total_frags) + 1
        ! Read in the coordinates
        frag_type = frag_list(ifrag,is)%type

        DO j = 1, frag_list(ifrag,is)%natoms
           
           this_atom = frag_list(ifrag,is)%atoms(j)
           
           atom_list(this_atom,this_im,is)%exist = .TRUE.
           nl = (frag_position_library(frag_type)-1) + &
                        frag_list(ifrag,is)%natoms*(this_fragment-1)+ j
           config_list(this_atom)%rxp = &
                 library_coords(nl)%rxp
              !frag_coords(j,this_fragment,frag_type)%rxp
           config_list(this_atom)%ryp = &
                 library_coords(nl)%ryp
              !frag_coords(j,this_fragment,frag_type)%ryp
           config_list(this_atom)%rzp = &
                 library_coords(nl)%rzp
              !frag_coords(j,this_fragment,frag_type)%rzp
         END DO
     
        ! For a ring fragment, access the fragment intramolecular energy 

        IF (frag_list(ifrag,is)%ring) THEN
           nrg_ring_frag =  nrg_frag(frag_type)%this_config_energy(this_fragment)
                              ! nrg_frag(this_fragment,frag_type)
           nrg_ring_frag_tot = nrg_ring_frag_tot + nrg_ring_frag
        END IF
        
        CLOSE(UNIT=11)

     END IF
     
     !**************************************************************************
     ! Step 2) Align the fragment to the growing molecule
     !**************************************************************************
     !

     ! Note that there has to be only one connection of ifrag that is already 
     ! placed. Let us find out which fragment that is

     total_connect = 0
     connection(:) = 0

     DO j = 1, frag_list(ifrag,is)%nconnect
        
        frag_connect = frag_list(ifrag,is)%frag_connect(j)
        
        IF (frag_placed(frag_connect) == 1 ) THEN
           total_connect = total_connect + 1
           connection(total_connect) = frag_connect
        END IF
        
     END DO

     IF (total_connect > 1 ) THEN
        err_msg = ''
        err_msg(1) = 'More than one connections of' // &
                     TRIM(Int_To_String(ifrag)) // 'exist'
        CALL Clean_Abort(err_msg,'Fragment_Placement')
     END IF
     
     ! If here then only one connection found

     frag_connect = connection(1)
    
     ! Note that frag_connect already has anchor of ifrag placed both for fixed 
     ! and variable bond length cases so all we have to do is obtain the 
     ! coordinates of these atoms in the configuration and align it to their 
     ! coordinates in the simulation box.
     
     ! find anchor atom ids for both the fragments

     CALL Get_Common_Fragment_Atoms(is,ifrag,frag_connect,anchor_ifrag, &
             anchor_frag_connect)
     nfrag_atoms = 0
     atom_id(:) = 0
     
     DO j = 1, frag_list(ifrag,is)%natoms
       
        this_atom = frag_list(ifrag,is)%atoms(j)
        
        IF ( (this_atom /= anchor_ifrag) .AND. &
             (this_atom /= anchor_frag_connect)) THEN
           nfrag_atoms = nfrag_atoms + 1
           atom_id(nfrag_atoms) = this_atom
        END IF
        
     END DO
     
!     anchor_ifrag = frag_list(ifrag,is)%anchor(1)
!     anchor_frag_connect = frag_list(frag_connect,is)%anchor(1)

     ! Find one atom of ifrag and frag_connect that will be used for generating 
     ! xy plane

     DO j = 1, frag_list(ifrag,is)%natoms

        this_atom = frag_list(ifrag,is)%atoms(j)

        IF ( (this_atom /= anchor_ifrag) .AND. &
             (this_atom /= anchor_frag_connect)) EXIT

     END DO

     atom_ifrag = this_atom

     ! Similarly for frag_connect

     DO j = 1, frag_list(frag_connect,is)%natoms

        this_atom = frag_list(frag_connect,is)%atoms(j)
     
        IF ( (this_atom /= anchor_ifrag) .AND. &
             (this_atom /= anchor_frag_connect)) EXIT

     END DO

     
     atom_frag_connect = this_atom

     ! Now use three atoms to obtain aligner and hanger matrix for the two 
     ! fragments: 

     ! atom1 == origin
     ! atom2 == id of the atom along which +ve x - axis is aligned
     ! atom3 == helps to obtain y-axis so that atom1-atom2-atom3 define xy plane

     ! vec1 == r_atom2 - r_atom1
     ! vec2 == r_atom3 - r_atom1

     ! form these vectors from the configuration read in from reservoir

     ! for ifrag

     vec1(1) = config_list(anchor_ifrag)%rxp &
             - config_list(anchor_frag_connect)%rxp
     vec1(2) = config_list(anchor_ifrag)%ryp &
             - config_list(anchor_frag_connect)%ryp
     vec1(3) = config_list(anchor_ifrag)%rzp &
             - config_list(anchor_frag_connect)%rzp


     vec2(1) = config_list(atom_ifrag)%rxp &
             - config_list(anchor_frag_connect)%rxp
     vec2(2) = config_list(atom_ifrag)%ryp &
             - config_list(anchor_frag_connect)%ryp
     vec2(3) = config_list(atom_ifrag)%rzp &
             - config_list(anchor_frag_connect)%rzp
       

     CALL Get_Aligner_Hanger(vec1, vec2, aligner_ifrag,hanger_ifrag)
     
     ! Calculate this only for inserting a molecule
     
     IF ( .NOT. del_flag) THEN
        
        
        vec1(1) = atom_list(anchor_ifrag,this_im,is)%rxp - &
             atom_list(anchor_frag_connect,this_im,is)%rxp
        
        vec1(2) = atom_list(anchor_ifrag,this_im,is)%ryp - &
             atom_list(anchor_frag_connect,this_im,is)%ryp
        
        vec1(3) = atom_list(anchor_ifrag,this_im,is)%rzp - &
             atom_list(anchor_frag_connect,this_im,is)%rzp
        
        vec2(1) = atom_list(atom_frag_connect,this_im,is)%rxp - &
             atom_list(anchor_frag_connect,this_im,is)%rxp
        
        vec2(2) = atom_list(atom_frag_connect,this_im,is)%ryp - &
             atom_list(anchor_frag_connect,this_im,is)%ryp
        
        vec2(3) = atom_list(atom_frag_connect,this_im,is)%rzp - &
             atom_list(anchor_frag_connect,this_im,is)%rzp     
     

        CALL Get_Aligner_Hanger(vec1, vec2, aligner_frag_connect,hanger_frag_connect)

     END IF
        
     ! Apply aligner of ifrag and then hanger of frag_connect to join the two 
     ! fragments.
     
     ! Apply aligner of ifrag
     

     tempx = config_list(anchor_frag_connect)%rxp
     tempy = config_list(anchor_frag_connect)%ryp
     tempz = config_list(anchor_frag_connect)%rzp
     
     config_list(:)%rxp = config_list(:)%rxp - tempx
     config_list(:)%ryp = config_list(:)%ryp - tempy
     config_list(:)%rzp = config_list(:)%rzp - tempz
         
     DO j = 1, frag_list(ifrag,is)%natoms
        
        this_atom = frag_list(ifrag,is)%atoms(j)
        
        tempx = config_list(this_atom)%rxp
        tempy = config_list(this_atom)%ryp
        tempz = config_list(this_atom)%rzp
        
        config_list(this_atom)%rxp = tempx * aligner_ifrag(1,1) &
                                   + tempy * aligner_ifrag(1,2) &
                                   + tempz * aligner_ifrag(1,3)
        config_list(this_atom)%ryp = tempx * aligner_ifrag(2,1) &
                                   + tempy * aligner_ifrag(2,2) &
                                   + tempz * aligner_ifrag(2,3)
        config_list(this_atom)%rzp = tempx * aligner_ifrag(3,1) &
                                   + tempy * aligner_ifrag(3,2) &
                                   + tempz * aligner_ifrag(3,3)
     END DO
     
     !**************************************************************************
     ! Step 3) Select kappa_dih orientations for each additional fragment
     !**************************************************************************
     ! 
     ! At this point, we can generate kappa_dih positions of the fragment as two 
     ! anchor positions are aligned. This is, in effect, equivalent to rotating 
     ! the non-anchor atoms around the x-axis. Note that the coordinates of the 
     ! anchoring atoms do not change due to this rotation. For the deletion 
     ! move, the first trial must be the one corresponding to the actual 
     ! coordinates, hence there should be no rotation about x-axis.

     IF ( del_flag ) THEN

        theta = 0.0_DP

        ! also note that we will transform the position based on hanger_ifrag so
        ! that the original positions are recovered

        hanger_frag_connect(:,:) = hanger_ifrag(:,:)

     ELSE

        ! Select a random theta with uniform probability
        
        theta = twopi * rranf()
     END IF

     ! Now that we have a starting theta, the other dihedral positions are 
     ! uniformly spaced around the 2pi disc. We need to calculate the atomic
     ! coordinates:

     ! Loop over the trial dihedrals
     ii = 1
     DO 
!     DO ii = 1, kappa_dih

        config_temp_list(:,ii)%rxp  = config_list(:)%rxp
        config_temp_list(:,ii)%ryp  = 0.0_DP
        config_temp_list(:,ii)%rzp  = 0.0_DP
        ! Loop over atoms
        DO j = 1, frag_list(ifrag,is)%natoms
           
           this_atom = frag_list(ifrag,is)%atoms(j)
           
           tempx = config_list(this_atom)%rxp
           tempy = config_list(this_atom)%ryp
           tempz = config_list(this_atom)%rzp
           config_temp_list(this_atom,ii)%ryp =  DCOS(theta) * tempy &
                                              +  DSIN(theta) * tempz
           config_temp_list(this_atom,ii)%rzp = -DSIN(theta) * tempy &
                                              +  DCOS(theta) * tempz

        END DO

        ! Loop over atoms (again)
        

             DO j = 1, frag_list(ifrag,is)%natoms
           
           this_atom = frag_list(ifrag,is)%atoms(j)
           
           tempx = config_temp_list(this_atom,ii)%rxp
           tempy = config_temp_list(this_atom,ii)%ryp
           tempz = config_temp_list(this_atom,ii)%rzp
           
           config_temp_list(this_atom,ii)%rxp = tempx*hanger_frag_connect(1,1) &
                                              + tempy*hanger_frag_connect(1,2) &
                                              + tempz*hanger_frag_connect(1,3)
           
           config_temp_list(this_atom,ii)%ryp = tempx*hanger_frag_connect(2,1) &
                                              + tempy*hanger_frag_connect(2,2) &
                                              + tempz*hanger_frag_connect(2,3) 
           
           config_temp_list(this_atom,ii)%rzp = tempx*hanger_frag_connect(3,1) &
                                              + tempy*hanger_frag_connect(3,2) &
                                              + tempz*hanger_frag_connect(3,3)

    
           IF ( this_atom /= anchor_ifrag) THEN
              IF ( this_atom /= anchor_frag_connect) THEN
                 
                 config_temp_list(this_atom,ii)%rxp = &
                      config_temp_list(this_atom,ii)%rxp + &
                      atom_list(anchor_frag_connect,this_im,is)%rxp
                 
                 config_temp_list(this_atom,ii)%ryp = &
                      config_temp_list(this_atom,ii)%ryp + &
                      atom_list(anchor_frag_connect,this_im,is)%ryp
                 
                 config_temp_list(this_atom,ii)%rzp = &
                      config_temp_list(this_atom,ii)%rzp + &
                      atom_list(anchor_frag_connect,this_im,is)%rzp
                 
                 atom_list(this_atom,this_im,is)%rxp = &
                      config_temp_list(this_atom,ii)%rxp 
                 
                 atom_list(this_atom,this_im,is)%ryp = &
                      config_temp_list(this_atom,ii)%ryp 
                 
                 atom_list(this_atom,this_im,is)%rzp = &
                      config_temp_list(this_atom,ii)%rzp  

              END IF
           END IF
        END DO

        ! Exit the loop if we've computed atomic coords for all trial dihedrals
        IF( ii == kappa_dih ) EXIT

        ! Increment the counter and dihedral angle
        ii = ii + 1
        theta = theta + twopi / REAL(kappa_dih,DP)

     END DO

     !**************************************************************************
     ! Step 4) Compute the energy of the fragment 
     !**************************************************************************
     ! 
     ! Initialize the energies

     nrg(:) = 0.0_DP
     nrg_dihed(:) = 0.0_DP
     overlap_trial(:) = .FALSE.

     DO ii = 1, kappa_dih

        ! Reload the coordinates for the atoms of this fragment
        DO j = 1, frag_list(ifrag,is)%natoms
           
           this_atom = frag_list(ifrag,is)%atoms(j)
           
           IF (this_atom /= anchor_ifrag) THEN
              IF (this_atom /= anchor_frag_connect) THEN
                 atom_list(this_atom,this_im,is)%rxp = &
                    config_temp_list(this_atom,ii)%rxp 
                 atom_list(this_atom,this_im,is)%ryp = &
                    config_temp_list(this_atom,ii)%ryp 
                 atom_list(this_atom,this_im,is)%rzp = &
                    config_temp_list(this_atom,ii)%rzp  
              END IF
           END IF
        END DO
 
        CALL Get_COM(this_im,is)
        CALL Compute_Max_COM_Distance(this_im,is)

        ! Turn all the atoms off
        overlap = .FALSE.
        DO j = 1, nfrag_atoms
           atom_list(atom_id(j),this_im,is)%exist = .FALSE.
        END DO
           
        ! Initialize the energies
        nrg_intra_vdw = 0.0_DP
        nrg_intra_qq = 0.0_DP
        nrg_inter_vdw = 0.0_DP
        nrg_inter_qq = 0.0_DP
           
        ! Compute the atomic energies as the fragment is slowly turned on
        DO j = 1, nfrag_atoms
              
           atom_list(atom_id(j),this_im,is)%exist = .TRUE.
              
           CALL Compute_Atom_Nonbond_Energy(atom_id(j),this_im,is, &
                E_intra_vdw,E_inter_vdw,E_intra_qq,E_inter_qq,overlap)
           IF (overlap) THEN
              ! if it is the last trial, the atom exist flag may not be
              ! properly set to true
              IF ( ii == kappa_dih)  THEN                    
                 DO k = 1, nfrag_atoms
                    atom_list(atom_id(k),this_im,is)%exist = .TRUE.
                 END DO
              END IF
              
              EXIT

           END IF
           
           nrg_intra_vdw = nrg_intra_vdw + E_intra_vdw
           nrg_intra_qq  = nrg_intra_qq  + E_intra_qq
           nrg_inter_vdw = nrg_inter_vdw + E_inter_vdw
           nrg_inter_qq  = nrg_inter_qq  + E_inter_qq
              
        END DO

        IF (overlap) THEN
           ! atoms are too close, set the weight to zero
           weight(ii) = 0.0_DP
           overlap_trial(ii) = .TRUE.
        ELSE
           CALL Compute_Molecule_Dihedral_Energy(this_im,is,e_dihed)

           nrg_dihed(ii) = e_dihed

           nrg(ii) = e_dihed + nrg_intra_vdw + nrg_intra_qq  + &          
                nrg_inter_vdw + nrg_inter_qq - e_prev

           IF (frag_list(ifrag,is)%ring) THEN
              ! subtract the biasing energy used to sample intramolecular DOFs
              nrg(ii) = nrg(ii) - nrg_ring_frag
           END IF
           
           nrg_kBT = beta(this_box) * nrg(ii)
           weight(ii) = DEXP(-nrg_kBT)
           
           IF ( nrg_kBT >= max_kBT) THEN
              ! the energy is too high, set the weight to zero
              weight(ii) = 0.0_DP
           ELSE
              weight(ii) = DEXP(-nrg_kBT)
           END IF
        END IF

!        ! BEGIN DEBUGGING OUTPUT
!        ! Write out the fragment coordinates for each trial position
!        M_XYZ_unit = movie_xyz_unit + this_box
!        DO j = 1, frag_list(ifrag,is)%natoms
!           this_atom = frag_list(ifrag,is)%atoms(j)
!           WRITE(M_XYZ_unit,*) &
!              TRIM(nonbond_list(this_atom,is)%element) // &
!              TRIM(int_to_string(ii)), & 
!              atom_list(this_atom,this_im,is)%rxp, &
!              atom_list(this_atom,this_im,is)%ryp, &
!              atom_list(this_atom,this_im,is)%rzp
!        END DO
!        IF (ii==1) THEN
!           WRITE(*,'(A,X,I5)'), 'i_mcstep', i_mcstep
!           WRITE(*,'(4(A12,X))') 'DIH' // TRIM(int_to_string(i)) // ':trial ', 'energy', 'weight', 'overlap'
!        END IF
!        WRITE(*,'(I12,X,E12.6,X,E12.6,X,L12)') ii, beta(this_box)*nrg(ii), weight(ii), overlap_trial(ii)
!        ! END DEBUGGING OUTPUT
                 
        ! Track the cumulative weights for Golden sampling
        IF ( ii > 1 ) weight(ii) = weight(ii-1) + weight(ii)

!        WRITE(*,*) 'weight',overlap, ii, nrg_kBT
!        WRITE(*,*) 'energy', e_prev, nrg(ii)
           
     END DO

     !**************************************************************************
     ! Step 5) Select a trial dihedral using the weighted probabilities
     !**************************************************************************
     ! 
     ! 

     ! If the cumulative weight is 0, then all trial dihedrals had core overlap
     ! Reject the move if all trials tripped overlap
     IF (ALL(overlap_trial)) THEN
        cbmc_overlap = .TRUE.
        RETURN
     END IF


     IF (del_flag) THEN

        ! Use trial 1, which holds the current coordinates of the fragment
        trial = 1

     ELSE

        ! Select a trial from the weighted distribution
        prob_pick = rranf() * weight(kappa_dih)

        DO ii = 1, kappa_dih

           IF ( prob_pick <= weight(ii) ) EXIT

        END DO

        trial = ii

     END IF

     ! Recover the individual probability for the accepted trial
     ln_pbias = ln_pbias - beta(this_box) * nrg(trial) - DLOG(weight(kappa_dih))
     
     ! Give the coordinates of this conformation to atom_list
     DO j = 1, frag_list(ifrag,is)%natoms
        
        this_atom = frag_list(ifrag,is)%atoms(j)
        
        IF (this_atom /= anchor_ifrag) THEN
           IF (this_atom /= anchor_frag_connect) THEN
              atom_list(this_atom,this_im,is)%rxp = &
                 config_temp_list(this_atom,trial)%rxp
              atom_list(this_atom,this_im,is)%ryp = &
                 config_temp_list(this_atom,trial)%ryp
              atom_list(this_atom,this_im,is)%rzp = &
                 config_temp_list(this_atom,trial)%rzp
           END IF
        END IF
     END DO
     
     ! also store the total energy upto this point
     e_prev =  nrg_dihed(trial)
     
     ! mark this fragment as placed
     frag_placed(ifrag) = 1
    

  END DO

  DEALLOCATE(counted)
  DEALLOCATE(config_list)
  DEALLOCATE(config_temp_list)
  DEALLOCATE(connection)
 END SUBROUTINE Fragment_Placement
                 
!***************************************************************************************************
SUBROUTINE Get_Aligner_Hanger(vec1,vec2,aligner,hanger)
!***************************************************************************************************


   USE Type_Definitions, ONLY : DP
   
   IMPLICIT NONE
   
   REAL(DP) :: vec1(3), vec2(3), perp_vec1(3), perp_vec2(3)
   REAL(DP), INTENT(OUT) :: aligner(3,3), hanger(3,3)
  
   INTEGER :: i, j
   

   aligner(:,:) = 0.0_DP
   hanger(:,:) = 0.0_DP
   
   ! Normalize the vectors
   
   vec1(:) = vec1(:) / DSQRT(DOT_PRODUCT(vec1,vec1))
   vec2(:) = vec2(:) / DSQRT(DOT_PRODUCT(vec2,vec2))
   
   ! Now cross the vector to get the vector perpendicular to these vectors
   
   perp_vec1(1) =  vec1(2) * vec2(3) - vec1(3) * vec2(2)
   perp_vec1(2) = -vec1(1) * vec2(3) + vec1(3) * vec2(1)
   perp_vec1(3) =  vec1(1) * vec2(2) - vec1(2) * vec2(1)
   
   perp_vec1(:) = perp_vec1(:) / DSQRT(DOT_PRODUCT(perp_vec1,perp_vec1))
   
   
   perp_vec2(1) =  vec1(2) * perp_vec1(3) - vec1(3) * perp_vec1(2)
   perp_vec2(2) = -vec1(1) * perp_vec1(3) + vec1(3) * perp_vec1(1)
   perp_vec2(3) =  vec1(1) * perp_vec1(2) - vec1(2) * perp_vec1(1)
   
   perp_vec2(:) = perp_vec2 / DSQRT(DOT_PRODUCT(perp_vec2,perp_vec2))
   
   DO i = 1, 3
      
      aligner(1,i) = vec1(i)
      aligner(2,i) = perp_vec2(i)
      aligner(3,i) = perp_vec1(i)
      
      hanger(i,1) = vec1(i) 
      hanger(i,2) = perp_vec2(i)
      hanger(i,3) = perp_vec1(i)
      
   END DO
   
   
 END SUBROUTINE Get_Aligner_Hanger
!***********************************************************
!
!*********************************************************************
SUBROUTINE Single_Fragment_Regrowth(alive,is)
!*********************************************************************
   ! This routine is used when a species contains only one fragment
   ! and a change in intramolecular degrees of freedom is desired. 
   !
   ! Written by Jindal K. Shah on 01/01/09
   !
   ! Algorithm
   !
   ! 1. Pick a configuration from the reservoir library
   ! 2. Insert the fragment such that the COM coincides with the fragment
   !    being "deleted"
   ! 3. Give it a random Eulerian rotation.
   !*****************************************************************

   
   USE Rotation_Routines
   
   INTEGER, INTENT(IN) :: alive, is

   INTEGER :: nl, nlo      ! number of the line where start the x,y,x coords of
                           ! config and fragment randomly selected
   INTEGER :: total_frags, i, this_atom, this_fragment, frag_type

   REAL(DP) :: temp_var, E_ang,x_this, y_this, z_this

   ! As this molecule has only 1 fragment, we obtain total number
   ! of fragments and atoms using 1 as the identity

   total_frags = frag_list(1,is)%nconfig
     
   ! Choose a fragment at random
   
   this_fragment = INT(rranf() * total_frags) + 1
   
   frag_type = frag_list(1,is)%type

   DO i = 1, frag_list(1,is)%natoms 
      
      this_atom = frag_list(1,is)%atoms(i)
      nl = (frag_position_library(frag_type)-1) + &
                                           frag_list(1,is)%natoms*(this_fragment -1) + i      
      atom_list(this_atom,alive,is)%rxp =  library_coords(nl)%rxp
                                             !frag_coords(i,this_fragment,frag_type)%rxp
      atom_list(this_atom,alive,is)%ryp =  library_coords(nl)%ryp
                                             !frag_coords(i,this_fragment,frag_type)%ryp
      atom_list(this_atom,alive,is)%rzp =  library_coords(nl)%rzp
                                              !frag_coords(i,this_fragment,frag_type)%rzp
   END DO



   ! COM and max_com_distance. Note that the following calls destroy
   ! original COM and max_com_distance but we will use xcom_old etc to
   ! transfer the COM to the original position

   CALL Get_COM(alive,is)
   CALL Compute_Max_COM_Distance(alive,is)

   ! Give random orientation to the fragment

   CALL Rotate_Molecule_Eulerian(alive,is)

   ! Move the atoms such that COM gets translated to original position

   DO i = 1, frag_list(1,is)%natoms
      this_atom = frag_list(1,is)%atoms(i)

      atom_list(this_atom,alive,is)%rxp = atom_list(this_atom,alive,is)%rxp + &
           molecule_list(alive,is)%xcom_old - molecule_list(alive,is)%xcom

      atom_list(this_atom,alive,is)%ryp = atom_list(this_atom,alive,is)%ryp + &
           molecule_list(alive,is)%ycom_old - molecule_list(alive,is)%ycom

      atom_list(this_atom,alive,is)%rzp = atom_list(this_atom,alive,is)%rzp + &
           molecule_list(alive,is)%zcom_old - molecule_list(alive,is)%zcom

   END DO
   CALL Get_COM(alive,is)
   CALL Compute_Max_COM_Distance(alive,is)


 END SUBROUTINE Single_Fragment_Regrowth

 !***********************************************************************************
SUBROUTINE Get_Common_Fragment_Atoms(is,frag1,frag2,atom1,atom2)
   !*********************************************************************************
   !
   ! This routine determines atoms connecting two fragments
   !
   ! written by Jindal Shah on 04/16/09
   !
   !*********************************************************************************
   
   INTEGER, INTENT(IN) :: is, frag1, frag2
   
   INTEGER, INTENT(OUT) :: atom1, atom2
   
   INTEGER :: i, j, this_atom, nanchor1, nanchor2, atom_i
   
   nanchor1 = frag_list(frag1,is)%nanchors
   nanchor2 = frag_list(frag2,is)%nanchors
   
   IF ( nanchor1 == 1 .AND. nanchor2 == 1 ) THEN
      
      atom1 = frag_list(frag1,is)%anchor(1)
      atom2 = frag_list(frag2,is)%anchor(1)
      
      
   ELSE IF ( nanchor1 == 1) THEN
      
      atom1 = frag_list(frag1,is)%anchor(1)
      
      ! for the second fragment loop over all the anchors of frag2 with atoms of frag1
      
      anchor2_loop:DO i = 1, nanchor2
        
         atom_i = frag_list(frag2,is)%anchor(i)
         
         DO j = 1, frag_list(frag1,is)%natoms
            
            this_atom = frag_list(frag1,is)%atoms(j)
            
            IF ( atom_i == this_atom) EXIT anchor2_loop
            
         END DO
         
      END DO anchor2_loop
      
      atom2 = atom_i
      
      
   ELSE IF (nanchor2 == 1 ) THEN
      
      atom2 = frag_list(frag2,is)%anchor(1)
      
      ! for the first fragment loop over all the anchors of frag1 with atoms of frag2
      
      anchor1_loop: DO i = 1, nanchor1
         
         atom_i = frag_list(frag1,is)%anchor(i)
         
         DO j = 1, frag_list(frag2,is)%natoms
            
            this_atom = frag_list(frag2,is)%atoms(j)
            
            IF ( atom_i == this_atom ) EXIT anchor1_loop
            
         END DO
         
      END DO anchor1_loop
      
      atom1 = atom_i
      
      
   ELSE
      
      ! multiple anchor atoms in both the fragments
      
      anchor1a_loop: DO i = 1, nanchor1
         
         atom_i = frag_list(frag1,is)%anchor(i)
         
         DO j = 1, frag_list(frag2,is)%natoms
            
            this_atom = frag_list(frag2,is)%atoms(j)
            
            IF ( atom_i == this_atom ) EXIT anchor1a_loop
            
         END DO
         
      END DO anchor1a_loop
      
      atom1 = atom_i
      
      
      anchor2a_loop:DO i = 1, nanchor2
         
         atom_i = frag_list(frag2,is)%anchor(i)
         
         DO j = 1, frag_list(frag1,is)%natoms
            
            this_atom = frag_list(frag1,is)%atoms(j)
            
            IF ( atom_i == this_atom) EXIT anchor2a_loop
            
         END DO
         
      END DO anchor2a_loop
      
      atom2 = atom_i
      
      
   END IF

   
 END SUBROUTINE Get_Common_Fragment_Atoms

END MODULE Fragment_Growth
