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

  !******************************************************************************
  ! The module performs the tasks related to insertion, deletion, cut and regrowth
  ! of molecules based on fragment sampling. It also contains routine to obtain
  ! order of fragment placement.
  !
  ! CONTAINS
  !      
  ! Fragment_Order
  ! Build_Molecule
  ! Delete_Molecule
  ! Cut_Regrow
  ! Get_Aligner_Hanger
  !
  ! Used by
  !
  !   chempot
  !   cutNgrow
  !   deletion
  !   gemc_particle_transfer
  !   grow_molecules
  !   insertion
  !   main
  !   nvt_mc_ring_fragment
  !   zig_by_omega
  ! 
  ! Revision history
  !
  !   12/10/13 : Beta Release
  !*****************************************************************************

  USE Run_Variables
  USE Energy_Routines
  USE Random_Generators
  USE Read_Write_Checkpoint

  IMPLICIT NONE

CONTAINS

!******************************************************************************
SUBROUTINE Build_Molecule(this_im,is,this_box,frag_order,this_lambda,which_anchor,&
     attempt_prob, nrg_ring_frag_total, cbmc_overlap)
!******************************************************************************
!
! The subroutine builds the molecule from scratch, it also has the capability 
! to regrow part of the input molecule
!
! First written by Jindal Shah on 07/18/07
!
! 05/25/08 (JS) : First committed to the repository
!
! 01/20/0 (JS) : Boltzmann weight of the trial is computed irrespoective of
!                 its energy. Previously, an energy cutoff was used to set
!                 the Boltzmann weight
!
! 
!
!********************************************************************************

  USE Rotation_Routines
  USE IO_Utilities
  USE File_Names
  USE Energy_Routines

  INTEGER :: this_im, is, this_box
  INTEGER, DIMENSION(1:nfragments(is)) :: frag_order

  REAL(DP) :: attempt_prob
  REAL(DP), INTENT(IN) :: this_lambda
  REAL(DP), INTENT(OUT) :: nrg_ring_frag_total
  LOGICAL, INTENT(INOUT) :: cbmc_overlap

  INTEGER :: which_anchor
  !-------------------------------

  INTEGER :: is_fragments, frag_start, frag_start_old, frag_total, i, this_atom
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

  INTEGER :: itrial, trial, frag_type, n_frag_atoms

  REAL(DP) :: weight(MAX(kappa_ins,kappa_rot,kappa_dih)), rand_no, e_dihed, E_intra_vdw, &
       E_intra_qq, E_inter_vdw
  REAL(DP) :: E_inter_qq,e_total
  REAL(DP) :: nrg(MAX(kappa_ins,kappa_rot,kappa_dih)), nrg_kBT, time0, time1, nrg_ring_frag

  LOGICAL :: del_overlap

  Type(Atom_Class) :: rtrial(MAXVAL(natoms),0:MAX(kappa_ins,kappa_rot,kappa_dih))
  ! slit pore 

  LOGICAL :: framework_overlap
  REAL(DP) :: E_framework

  n_frag_atoms = 0
  CBMC_Flag = .TRUE.
  attempt_prob = 1.0_DP
  
  ! mark all the atoms as deleted 
  
  ! Assign a locate number for the molecule

  this_box = molecule_list(this_im,is)%which_box
  
  atom_list(:,this_im,is)%exist = .FALSE.

  is_fragments = nfragments(is)
  IF (ALLOCATED(frag_placed)) DEALLOCATE(frag_placed)
  ALLOCATE(frag_placed(is_fragments))
  frag_placed(:) = 0

  nrg_ring_frag_total = 0.0_DP

  IF (get_fragorder) THEN
 
     ! First obtain a random fragment to grow from
        
     frag_start = INT ( rranf() * nfragments(is)) + 1
     frag_start_old = frag_start
 
     ! Obtain the order of fragment addition
     
     IF (ALLOCATED(live)) DEALLOCATE(live)
     IF (ALLOCATED(deadend)) DEALLOCATE(deadend)
     IF (ALLOCATED(central)) DEALLOCATE(central)
     
     ALLOCATE(live(is_fragments),deadend(is_fragments),central(is_fragments))      

     live(:) = 0
     central(:) = 0
     deadend(:) = 0
     
     frag_order(:) = 0
     frag_order(1) = frag_start
     live(frag_start) = 1
     frag_total = 1
     
     IF (is_fragments > 1 ) THEN

        CALL Fragment_Order(frag_start,is,frag_total,frag_order,live,deadend,central)
        
     END IF

     ! Note that the call to Fragment_Order might have changed the identity of frag_start
     ! so restore the value calculated above
     
     frag_start = frag_start_old

     DEALLOCATE(live,deadend,central)

  END IF

  ! Note that for a deletion move, frag_total is undefined, since its an insertion move,
  ! total number of fragments that need to be placed is is_fragments

  frag_total = is_fragments

  ! At this point, we have the order in which we will grow the molecule from
     
  ! let us place the first fragment. We will use choose a configuration from the
  ! library and assign it to the atoms in the fragment. This fragment will undergo
  ! a COM rotation before its anchor is randomly inserted in the simulation box.

  ! Add the part to read in the coordinates from its file.

  ! we will make all the atoms of frag_start as part of the simulations

  atom_list(:,this_im,is)%exist = .FALSE.
  molecule_list(this_im,is)%cfc_lambda = 0.0_DP

  frag_start = frag_order(1)

  ! Note that we need to choose from the reservoir only when insertion
  ! is attempted

  IF (.NOT. del_Flag) THEN
     
     ! obtain a random configuration

     total_frags = frag_list(frag_start,is)%nconfig
     
     ! Choose a fragment at random
     
     this_fragment = INT(rranf() * total_frags) + 1

     frag_type = frag_list(frag_start,is)%type
          
     DO i = 1, frag_list(frag_start,is)%natoms 
        
        this_atom = frag_list(frag_start,is)%atoms(i)
        
        atom_list(this_atom,this_im,is)%rxp = frag_coords(i,this_fragment,frag_type)%rxp
        atom_list(this_atom,this_im,is)%ryp = frag_coords(i,this_fragment,frag_type)%ryp
        atom_list(this_atom,this_im,is)%rzp = frag_coords(i,this_fragment,frag_type)%rzp

     END DO

  END IF
  
  molecule_list(this_im,is)%cfc_lambda = this_lambda
  
  DO i =1, frag_list(frag_start,is)%natoms
     this_atom = frag_list(frag_start,is)%atoms(i)
     atom_list(this_atom,this_im,is)%exist = .TRUE.


  END DO

 ! Now apply the rotation matrix so that the resulting orientation is random

  CALL Get_COM(this_im,is)
  CALL Compute_Max_COM_Distance(this_im,is)
  
  IF ( .NOT. del_FLAG) CALL Rotate_Molecule_Eulerian(this_im,is)

  ! Store the position as 0th trial
  
  DO i = 1, frag_list(frag_start,is)%natoms

     this_atom = frag_list(frag_start,is)%atoms(i)
     
     rtrial(this_atom,0)%rxp = atom_list(this_atom,this_im,is)%rxp
     rtrial(this_atom,0)%ryp = atom_list(this_atom,this_im,is)%ryp
     rtrial(this_atom,0)%rzp = atom_list(this_atom,this_im,is)%rzp

  END DO

  ! store the COM
  
  xcom_old = molecule_list(this_im,is)%xcom
  ycom_old = molecule_list(this_im,is)%ycom
  zcom_old = molecule_list(this_im,is)%zcom

  ! We will place this fragment based only on its external weight

  nrg(:) = 0.0_DP 

  IF(imreplace .GT. 0) THEN

     IF(.NOT. del_FLAG) THEN

        dx = molecule_list(imreplace,isreplace)%xcom - molecule_list(this_im,is)%xcom
        dy = molecule_list(imreplace,isreplace)%ycom - molecule_list(this_im,is)%ycom
        dz = molecule_list(imreplace,isreplace)%zcom - molecule_list(this_im,is)%zcom

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

     DO itrial = 1, kappa_ins

        IF ( del_Flag .AND. (itrial == 1 )) THEN
        
           ! choose old center of mass for fragment placement
           x_anchor = xcom_old
           y_anchor = ycom_old
           z_anchor = zcom_old
        
        ELSE
        
           IF (box_list(this_box)%int_box_shape == int_cubic) THEN
           
              x_anchor = (0.5_DP - rranf()) * box_list(this_box)%length(1,1)
              y_anchor = (0.5_DP - rranf()) * box_list(this_box)%length(2,2)
              z_anchor = (0.5_DP - rranf()) * box_list(this_box)%length(3,3)

           END IF

        END IF
     
        ! move all the atoms of the fragment
     
        DO i = 1, frag_list(frag_start,is)%natoms
        
           this_atom = frag_list(frag_start,is)%atoms(i)
        
              ! displace the atom
           atom_list(this_atom,this_im,is)%rxp = rtrial(this_atom,0)%rxp - &
                xcom_old + x_anchor
           
           atom_list(this_atom,this_im,is)%ryp = rtrial(this_atom,0)%ryp - &
                ycom_old + y_anchor
           
           atom_list(this_atom,this_im,is)%rzp = rtrial(this_atom,0)%rzp - &
                zcom_old + z_anchor
           
           rtrial(this_atom,itrial)%rxp = atom_list(this_atom,this_im,is)%rxp
           rtrial(this_atom,itrial)%ryp = atom_list(this_atom,this_im,is)%ryp
           rtrial(this_atom,itrial)%rzp = atom_list(this_atom,this_im,is)%rzp
        
        END DO

        molecule_list(this_im,is)%xcom = x_anchor
        molecule_list(this_im,is)%ycom = y_anchor
        molecule_list(this_im,is)%zcom = z_anchor

        ! If it's a slit pore simulation, check the atomic coordinates to
        ! ensure that the entire molecule is confined inside the pore

        ! Carry out the following only if there all the atoms of the fragment are inside the 
        ! framework
        overlap = .FALSE.

           ! Note that the COM position is always chosen inside the simulation box so there is
           ! No need for Fold_Molecule call.
           
           ! For the sake of completeness we calculate the distance of the atom furtherst from the COM
           CALL Compute_Max_COM_Distance(this_im,is)
           
           ! calculate the intermolecular energy of the fragment. Note that the
           ! CBMC flag has been set to true so that the following call will compute
           ! interaction energy of the growing molecule within a small distance
                                 
           CALL Compute_Molecule_Nonbond_Inter_Energy(this_im,is,E_inter_vdw,E_inter_qq,overlap)

          
           nrg(itrial) = nrg(itrial) + E_inter_vdw + E_inter_qq 


        IF (overlap) THEN
           ! the energy is too high, set the weight to zero

           weight(itrial) = 0.0_DP

        ELSE

           nrg_kBT = beta(this_box) * nrg(itrial)
!!$
!!$        ! compute the weight of this trial for the reverse move and first
!!$        ! trial. The following IF-ELSE construct ensures that the weight
!!$        ! for the reverse move is always computed.
!!$
           IF ( del_Flag .AND. (itrial == 1) ) THEN

              weight(itrial) = DEXP(-nrg_kBT)

           ELSE IF ( nrg_kBT >= max_kBT) THEN

              weight(itrial) = 0.0_DP

           ELSE

              weight(itrial) = DEXP(-nrg_kBT)
           
           END IF


        END IF

        IF (itrial > 1 ) weight(itrial) = weight(itrial-1) + weight(itrial)
     
     
     ! Note that if the overlap flag is true, the following update is not
     ! true. But, it does not matter as the trial will not be picked due
     ! to its zero weight.
    
     END DO


     ! Reject the move is the total weight is still zero

     IF (weight(kappa_ins) == 0.0_DP) THEN

        cbmc_overlap = .TRUE.
        CBMC_Flag = .FALSE.
        RETURN

     END IF

     attempt_prob = 1.0_DP

     IF (del_Flag) THEN

        trial = 1

     ELSE

     ! Choose one from Golden sampling for an insertion move
     
        rand_no = rranf() * weight(kappa_ins)
     
        DO i = 1, kappa_ins
           IF ( rand_no < weight(i)) EXIT
        END DO
     
        trial = i

     END IF

     IF ( trial == kappa_ins + 1 ) THEN
        ! it means that none of the trials could be picked
        ! and this could be due to the fact that all the trials
        ! generated resulted in a very small weight

        cbmc_overlap = .TRUE.
        CBMC_Flag = .FALSE.
        RETURN
        
     END IF

     IF (trial == 1) THEN
        attempt_prob = attempt_prob * weight(1)/ weight(kappa_ins)
     ELSE
        attempt_prob = attempt_prob *  (weight(trial) - weight(trial-1))/weight(kappa_ins)
     END IF
  
     e_total = nrg(trial)

     ! We chose the ith trial position for placement. So give the trial positions
     ! to the atoms of the first fragment. Note that for the deletion move,
     ! the first trial is picked so that the following assignment for positions is
     ! properly done.
  
     DO i = 1, frag_list(frag_start,is)%natoms
     
        this_atom = frag_list(frag_start,is)%atoms(i)
     
        atom_list(this_atom,this_im,is)%rxp = rtrial(this_atom,trial)%rxp
        atom_list(this_atom,this_im,is)%ryp = rtrial(this_atom,trial)%ryp
        atom_list(this_atom,this_im,is)%rzp = rtrial(this_atom,trial)%rzp
     
     
     ! we compute COM of the fully grown molecule
     END DO

     CALL Get_COM(this_im,is)

  END IF
  
  ! Mark this fragment is placed
  e_total = 0.0_DP
  frag_placed(frag_start) = 1



  ! We have our first segment placed in the system If this is a ring fragment then
  ! we need intramolecular energies of the fragment used in biasing.
  ! For insertion, it is the energy in the gas phase. For, deletion move, it is the
  ! energy of the old configuration. Note that this energy does not include angle
  ! bending energy as it cancels out from the acceptance rule

  IF (frag_list(frag_start,is)%ring) THEN

     IF (del_Flag) THEN
        ! compute the old energy
        CALL Compute_Ring_Fragment_Energy(frag_start,this_im,is,this_box,nrg_ring_frag)
     ELSE
        nrg_ring_frag = nrg_frag(this_fragment,frag_type)
     END IF
     nrg_ring_frag_total = nrg_ring_frag_total + nrg_ring_frag
     
  END IF

  ! Now we will place rest of the segments based on the initial fragment placed

  CALL Compute_Molecule_Dihedral_Energy(this_im,is,e_dihed)
  e_total = e_dihed

  del_overlap = .FALSE.


!  IF ( .NOT. del_flag) THEN
!	 write(*,*) atom_list(:,this_im,is)%rxp
!         write(*,*) atom_list(:,this_im,is)%ryp
!         write(*,*) atom_list(:,this_im,is)%rzp
!  END IF

!!$  write(*,*) 'fragment placed', this_fragment
!!$  DO i = 1, frag_list(this_fragment,is)%natoms
!!$     this_atom = frag_list(this_fragment,is)%atoms(i)
!!$     write(12,*) nonbond_list(this_atom,is)%element, atom_list(this_atom,this_im,is)%rxp, &
!!$          atom_list(this_atom,this_im,is)%ryp, atom_list(this_atom,this_im,is)%rzp
!!$  END DO
!!$  read(*,*)

 CALL Fragment_Placement(this_box,this_im,is,2,frag_total,frag_order,frag_placed,this_lambda, &
       e_total,attempt_prob,nrg_ring_frag_total, cbmc_overlap, del_overlap)


  ! Note that cbmc_overlap flag may be true and the CBMC_Flag will be properly assigned FALSE while 
  ! exiting the code.

 ! Note that the 'del_overlap' flag may be true while computing the weight of an old 
 ! configuration. So set cbmc_overlap is set to true so that the calling routine 
 ! appropriately sets the coordinates of the atoms that were not grown

 IF (del_overlap) cbmc_overlap = .TRUE.

  ! Mark the CBMC flag as false so that intermolecular nonbonded interactions are properly
  ! computed for the molecule


  CBMC_Flag = .FALSE.
  
END SUBROUTINE Build_Molecule

!******************************************************************************
SUBROUTINE Build_Rigid_Fragment(this_im,is,this_box,frag_order,this_lambda, &
     which_anchor,attempt_prob, nrg_ring_frag_total, cbmc_overlap)
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

  INTEGER :: this_im, is, this_box, first_atom,which_anchor
  INTEGER, DIMENSION(1:nfragments(is)) :: frag_order
  REAL(DP) :: attempt_prob
  REAL(DP), INTENT(IN) :: this_lambda
  REAL(DP), INTENT(OUT) :: nrg_ring_frag_total
  LOGICAL, INTENT(INOUT) :: cbmc_overlap
  !-------------------------------

  INTEGER :: is_fragments, frag_start, frag_start_old, frag_total, i, this_atom
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
  REAL(DP) :: E_inter_qq,e_total
  REAL(DP) :: nrg(MAX(kappa_ins,kappa_rot,kappa_dih)), nrg_kBT, time0, time1, nrg_ring_frag
  LOGICAL :: del_overlap
  TYPE(Atom_Class) :: rtrial(MAXVAL(natoms),0:MAX(kappa_ins,kappa_rot,kappa_dih))



  weight(:)=0.0_DP
  CBMC_Flag = .TRUE.
  attempt_prob = 1.0_DP
  
  ! Assign a locate number for the molecule
  this_box = molecule_list(this_im,is)%which_box
  ! mark all the atoms as deleted 
  atom_list(:,this_im,is)%exist = .FALSE.

  is_fragments = nfragments(is)
  IF (ALLOCATED(frag_placed)) DEALLOCATE(frag_placed)
  ALLOCATE(frag_placed(is_fragments))
  frag_placed(:) = 0
  nrg_ring_frag_total = 0.0_DP

! Get the order of the fragment growth
  IF (get_fragorder) THEN
     frag_start = 1
     ! Obtain the order of fragment addition
     IF (ALLOCATED(live)) DEALLOCATE(live)
     IF (ALLOCATED(deadend)) DEALLOCATE(deadend)
     IF (ALLOCATED(central)) DEALLOCATE(central)
     ALLOCATE(live(is_fragments),deadend(is_fragments),central(is_fragments))      
     live(:) = 0
     central(:) = 0
     deadend(:) = 0
     frag_order(:) = 0
     frag_order(1) = frag_start
     live(frag_start) = 1
     frag_total = 1
     IF (is_fragments > 1 ) THEN
        CALL Fragment_Order(frag_start,is,frag_total,frag_order,live,deadend,central)
     END IF
     DEALLOCATE(live,deadend,central)
  END IF

  ! Note that for a deletion move, frag_total is undefined, since its an insertion move,
  ! total number of fragments that need to be placed is is_fragments

  frag_total = is_fragments

  ! At this point, we have the order in which we will grow the molecule from
  ! Add the part to read in the coordinates from its file.
  ! we will make all the atoms of frag_start as part of the simulations

  atom_list(:,this_im,is)%exist = .FALSE.
  molecule_list(this_im,is)%cfc_lambda = 0.0_DP

  frag_start = frag_order(1)


  ! Note that we need to choose from the reservoir only when insertion
  ! is attempted

  molecule_list(this_im,is)%cfc_lambda = this_lambda

  ! NR: Select the location of first bead of the fragment

  first_atom = frag_list(frag_start,is)%atoms(1)
  atom_list(first_atom,this_im,is)%exist = .true.


  DO itrial = 1, kappa_ins
     IF (del_FLAG .and.(itrial.eq.1)) THEN
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


!    WRITE(8,*) del_Flag,this_box, itrial, this_im,E_inter_vdw,  E_inter_qq,overlap

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
         IF ( del_Flag .AND. (itrial == 1) ) THEN
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
      CBMC_Flag = .FALSE.
      RETURN
  END IF

  attempt_prob = 1.0_DP

  IF (del_Flag) THEN
      trial = 1
!     write(*,*) weight(1), this_box, nrg
  ELSE
   ! Choose one from Golden sampling for an insertion move
      rand_no = rranf() * weight(kappa_ins)
      DO i = 1, kappa_ins
         IF ( rand_no < weight(i)) EXIT
      END DO
      trial = i
  END IF
!  Write(8,*) 'selected', trial

  IF (trial == 1) THEN
     attempt_prob = attempt_prob * weight(1)/ weight(kappa_ins)
  ELSE
     attempt_prob = attempt_prob *  (weight(trial) - weight(trial-1))/weight(kappa_ins)
  END IF

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

  IF (.NOT. del_Flag) THEN
     ! obtain a random configuration
     total_frags = frag_list(frag_start,is)%nconfig
     ! Choose a fragment at random
     this_fragment = INT(rranf() * total_frags) + 1
     frag_type = frag_list(frag_start,is)%type
     DO i = 1, frag_list(frag_start,is)%natoms 
        this_atom = frag_list(frag_start,is)%atoms(i)
        rtrial(this_atom,0)%rxp = frag_coords(i,this_fragment,frag_type)%rxp-&
                                  frag_coords(1,this_fragment,frag_type)%rxp+&
                                  atom_list(first_atom,this_im,is)%rxp 
        rtrial(this_atom,0)%ryp = frag_coords(i,this_fragment,frag_type)%ryp-&
                                  frag_coords(1,this_fragment,frag_type)%ryp+&
                                  atom_list(first_atom,this_im,is)%ryp
        rtrial(this_atom,0)%rzp = frag_coords(i,this_fragment,frag_type)%rzp-&
                                  frag_coords(1,this_fragment,frag_type)%rzp+&
                                  atom_list(first_atom,this_im,is)%rzp 
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
     IF(del_Flag .and. (itrial.eq.1)) THEN
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

!     Write(8,*) 'First fragment attempt vdw qq',itrial, E_inter_vdw, E_inter_qq,this_box,del_Flag,overlap

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
         IF ( del_Flag .AND. (itrial == 1) ) THEN
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
      CBMC_Flag = .FALSE.
      RETURN
  END IF

  IF (del_Flag) THEN
      trial = 1
!     write(*,*) weight(1), this_box, nrg
  ELSE
   ! Choose one from Golden sampling for an insertion move
      rand_no = rranf() * weight(kappa_rot)
      DO i = 1, kappa_rot
         IF ( rand_no < weight(i)) EXIT
      END DO
      trial = i
  END IF

!  Write(8,*) 'selected', trial

  IF (trial == 1) THEN
     attempt_prob = attempt_prob * weight(1)/ weight(kappa_rot)
  ELSE
     attempt_prob = attempt_prob *  (weight(trial) - weight(trial-1))/weight(kappa_rot)
  END IF

! Assign positions to the atom_list
  DO i=1,frag_list(frag_start,is)%natoms
     this_atom =  frag_list(frag_start,is)%atoms(i)
     atom_list(this_atom,this_im,is)%rxp = rtrial(this_atom,trial)%rxp
     atom_list(this_atom,this_im,is)%ryp = rtrial(this_atom,trial)%ryp
     atom_list(this_atom,this_im,is)%rzp = rtrial(this_atom,trial)%rzp
  END DO

! IF ( .NOT. del_flag) THEN
!  write(*,*) atom_list(:,this_im,is)%rxp
! write(*,*) atom_list(:,this_im,is)%rzp
! write(*,*) atom_list(:,this_im,is)%ryp
! NR: Now we have place the first fragment.
! write(*,*)
!	End if
  frag_placed(frag_start) = 1


  IF (frag_list(frag_start,is)%ring) THEN

     IF (del_Flag) THEN
        ! compute the old energy
        CALL Compute_Ring_Fragment_Energy(frag_start,this_im,is,this_box,nrg_ring_frag)
     ELSE
        nrg_ring_frag = nrg_frag(this_fragment,frag_type)
     END IF
     nrg_ring_frag_total = nrg_ring_frag_total + nrg_ring_frag

  END IF

  ! Now we will place rest of the segments based on the initial fragment placed

  CALL Compute_Molecule_Dihedral_Energy(this_im,is,e_dihed)
  e_total = e_dihed

  CALL Fragment_Placement(this_box,this_im,is,2,frag_total,frag_order,frag_placed,this_lambda, &
       e_total,attempt_prob,nrg_ring_frag_total, cbmc_overlap, del_overlap)

  CBMC_Flag = .FALSE.

END SUBROUTINE Build_Rigid_Fragment

!*****************************************************************************************
SUBROUTINE Cut_Regrow(this_im,is,frag_start,frag_end,frag_order,frag_total,this_lambda, &
     e_prev,attempt_prob,nrg_ring_frag_tot, cbmc_overlap, del_overlap)
  !**************************************************************************************
  ! The subroutine performs cuts part of a molecule and regrows using configurational
  ! biasing. 
  !
  ! Written by Jindal Shah on 08/26/08
  !**************************************************************************************

  USE Run_Variables
  USE Random_Generators

  IMPLICIT NONE

  INTEGER :: is_fragments, this_im, is, frag_bond, frag_start, frag_start_old, frag_end
  INTEGER :: frag_total, this_box, ifrag, frag_no, i, this_atom, anchor, this_frag
  INTEGER :: frag_order(1:nfragments(is)), anchor_start, anchor_end
  
  INTEGER, DIMENSION(:), ALLOCATABLE :: live, deadend, central
  INTEGER, DIMENSION(:), ALLOCATABLE :: frag_placed

  REAL(DP), INTENT(IN) :: this_lambda

  REAL(DP) :: e_total, attempt_prob, e_prev

  REAL(DP) :: E_intra_vdw, E_intra_qq, E_inter_vdw, E_inter_qq, e_dihed

  REAL(DP) :: nrg_ring_frag_tot

  LOGICAL :: cbmc_overlap, overlap, del_overlap


  CBMC_Flag = .TRUE.
  
  del_overlap = .FALSE.

  nrg_ring_frag_tot = 0.0_DP

  ! obtain the box in which molecule is selected

  this_box = molecule_list(this_im,is)%which_box

  is_fragments = nfragments(is)

  ! Choose a fragment bond to cut

  IF(ALLOCATED(frag_placed)) DEALLOCATE(frag_placed)
  ALLOCATE(frag_placed(is_fragments))

  IF ( .NOT. DEL_FLAG) THEN
     
     frag_bond = INT (rranf() * fragment_bonds(is)) + 1
     
     ! Select one of the fragments to grow from
     
     IF ( rranf() < 0.5_DP) THEN
        
        frag_start = fragment_bond_list(frag_bond,is)%fragment1
        frag_end = fragment_bond_list(frag_bond,is)%fragment2
        
     ELSE
        
        frag_start = fragment_bond_list(frag_bond,is)%fragment2
        frag_end = fragment_bond_list(frag_bond,is)%fragment1
        
     END IF
     
     ! Protect frag_start for fragment_order call
     
     frag_start_old = frag_start
     
     ! Allocate arrays 
     
     ALLOCATE(live(is_fragments),deadend(is_fragments),central(is_fragments))
 
     live(:) = 0
     central(:) = 0
     deadend(:) = 0
     frag_placed(:) = 0
     
     frag_order(:) = 0
     frag_order(1) = frag_start
     live(frag_start) = 1
     live(frag_end) = 1
     frag_total = 1
     
     deadend(frag_end) = 1

     ! Mark all the fragments connected to frag_end dead except frag_start

    
     DO ifrag = 1, frag_list(frag_end,is)%nconnect
        
        frag_no = frag_list(frag_end,is)%frag_connect(ifrag)
        
        IF ( frag_no /= frag_start) THEN
           
           deadend(frag_no) = 1
           
        END IF

     END DO
 
     ! Obtain random order of fragments to be regrown
     
     CALL Fragment_Order(frag_start,is,frag_total,frag_order,live,deadend,central)


     ! reassign the starting fragment

     frag_start = frag_start_old
     DEALLOCATE(live,deadend,central)

  END IF
     
  ! Locate all the fragments that are placed. To do this, mark all the fragemnts
  ! as placed and unmark the ones from frag_order to be placed.  
  !---
  atom_list(:,this_im,is)%exist = .TRUE.
  frag_placed(:) = 1
  
  DO ifrag = 1, frag_total

     this_frag = frag_order(ifrag)
     
     frag_placed(this_frag) = 0
     
     ! also delete all the atoms contained in this fragment
     
     DO i = 1, frag_list(this_frag,is)%natoms
        
        this_atom = frag_list(this_frag,is)%atoms(i)

        atom_list(this_atom,this_im,is)%exist = .FALSE.
        
     END DO

  END DO
 
  ! Note that the anchor of frag_start is present in the simulation so turn
  ! the exist flag true for this atom. In the process, we have also marked
  ! the anchor of frag_end as dead, so we will make that atom alive as well

 ! Get the two anchor atoms that need to be turned alive

  CALL Get_Common_Fragment_Atoms(is,frag_start,frag_end,anchor_start,anchor_end)

  atom_list(anchor_start,this_im,is)%exist = .TRUE.
  atom_list(anchor_end,this_im,is)%exist = .TRUE.
 

  attempt_prob = 1.0_DP

  ! At this point we will calculate the energy of the fragment in the simulation.
  ! this is done so that the energy calculated for bonded and nonbonded interactions
  ! do not cause trouble while growing the molecule. Since dihedral, intramolecular
  ! nonbond and intermolecular nonbond energies are used in biasing in fragment
  ! placement, we will compute these energies. Note that this is computed only when
  ! del_FLAG is false. We will use the energy computed here when del_FLAG is true.

  IF (.NOT. del_FLAG) THEN

     CALL Compute_Molecule_Dihedral_Energy(this_im,is,e_dihed)
!     CALL Compute_Molecule_Nonbond_Intra_Energy(this_im,is,E_intra_vdw,E_intra_qq)
!     CALL Compute_Molecule_Nonbond_Inter_Energy(this_im,is,E_inter_vdw,E_inter_qq,overlap)

!     e_prev = e_dihed + E_intra_vdw + E_intra_qq + E_inter_vdw + E_inter_qq
     e_prev = e_dihed

  END IF

  e_total = e_prev ! e_total is equal to the energy computed above during the growth phase.
                   ! while it is the value obtained for cut_N_grow for calculating weight of the old chain


  CALL Fragment_Placement(this_box,this_im,is,1,frag_total,frag_order,frag_placed,this_lambda, &
       e_total,attempt_prob,nrg_ring_frag_tot,  cbmc_overlap, del_overlap)

  CBMC_Flag = .FALSE.
  

END SUBROUTINE Cut_Regrow



!*******************************************************************************
RECURSIVE SUBROUTINE Fragment_Order(this_frag,is,frag_total,frag_order,live,deadend,central)
!*******************************************************************************
!
! The subroutine obtains the order in which a molecule of species 'is' will be
! (re)grown. The code is based on the recursive routine atoms_to_place.f90 and
! is adapted for fragment sampling
!
! First written by Jindal Shah on 07/17/08
!
! 08/25/08 (JS) : First added to the repository after validating against
!                 the following molecule
!
!               C
!               |
!               C
!               |
!            C- C - C
!               |
!            C- C - C
!
!*********************************************************************************
  USE Run_Variables
  USE Random_Generators
  

  IMPLICIT NONE

  INTEGER :: this_frag, is, frag_total, n_connect, frag_no, frag_id
  

  INTEGER, DIMENSION(1:nfragments(is)) :: live, deadend, central, frag_order
  INTEGER, DIMENSION(:), ALLOCATABLE :: counted

  REAL(DP) :: rand_no

  IF ( .NOT. ALLOCATED(counted)) ALLOCATE(counted(nfragments(is)))

  ! case 1, entire molecule is regrown
  ! case 2, only part of the molecule is regrown may be for this, I can form the
  ! fragment_to_place array if there is a need

  ! First obtain total number of connections in the species and randomly choose a
  ! fragment to insert

  IF ( deadend(this_frag) == 1 ) THEN

     ! We encountered a fragment from which we cannot grow so need to trace back
     ! Randomly identify a fragment to grow from

     ! since it a deadend fragment, all the fragments connected to it would have
     ! already been placed in a previous call

     
     counted(:) = 0

     n_connect = 0

     DO WHILE (n_connect < frag_list(this_frag,is)%nconnect)
        
        rand_no = rranf()
        frag_no = INT ( rand_no * frag_list(this_frag,is)%nconnect) + 1
        frag_id = frag_list(this_frag,is)%frag_connect(frag_no)

        ! make sure that we did not check this fragment before

        IF (counted(frag_id) /= 1 ) THEN
           n_connect = n_connect + 1

           counted(frag_id) = 1
           
           IF ( frag_list(frag_id,is)%nconnect /= 1 ) THEN
              IF ( deadend(frag_id) /= 1 ) THEN
                 
                 CALL Fragment_Order(frag_id,is,frag_total,frag_order,live,deadend,central)
                 
              END IF
           END IF

        END IF
           

        
     END DO
     
        
  ELSE

     ! if the fragments connected to this fragment have not already been placed then
     ! determine the order in which the fragments are to be placed

     IF (central(this_frag) /= 1 ) THEN
        
        central(this_frag) = 1
        
        n_connect = 0
        counted(:) = 0
        
        DO WHILE (n_connect < frag_list(this_frag,is)%nconnect)
                      
           frag_no = INT( rranf() * frag_list(this_frag,is)%nconnect ) + 1

           frag_id = frag_list(this_frag,is)%frag_connect(frag_no)
!           write(*,*) this_frag, frag_id
           IF (counted(frag_id) /= 1) THEN
              counted(frag_id) = 1
              n_connect = n_connect + 1
           
              ! IF this fragment has not been placed then count then make it live

              IF (live(frag_id) == 0) THEN
              
                 frag_total = frag_total + 1
                 frag_order(frag_total) = frag_id
                 live(frag_id) = 1 
                 
              END IF

           END IF
              
        END DO
        
     END IF
     
     ! Now randomly pick a connection to grow from
     
     counted(:) = 0
     n_connect = 0
     DO WHILE (n_connect < frag_list(this_frag,is)%nconnect)

        rand_no = rranf()
        
        frag_no = INT ( rand_no * frag_list(this_frag,is)%nconnect) + 1
        frag_id = frag_list(this_frag,is)%frag_connect(frag_no)
        
        
        IF ( counted(frag_id) == 0 ) THEN

           counted(frag_id) = 1

           n_connect = n_connect + 1
!           write(*,*) 'frag', frag_id
           ! check to make sure that it is not a terminal fragment
           
           IF ( frag_list(frag_id,is)%nconnect /= 1 ) THEN
              
              ! make sure that it is not a deadend
              
              IF ( deadend(frag_id) /= 1 ) THEN
                 
                 ! if the fragments connected to this are not already placed
                 
                 IF ( central(frag_id) /= 1 ) THEN
                    
                    this_frag = frag_id
                    CALL Fragment_Order(frag_id,is,frag_total,frag_order,live,deadend,central)
                    
                 END IF
                 
              END IF
              
           END IF
           
        END IF
        
     END DO
     
     ! if here then this_frag is a deadend so must trace back
     
     IF (deadend(this_frag) /= 1 ) THEN
        
        deadend(this_frag) = 1
        CALL Fragment_Order(this_frag,is,frag_total,frag_order,live,deadend,central)        
        
     END IF
     
  END IF

END SUBROUTINE Fragment_Order

!**************************************************************************************************
SUBROUTINE Fragment_Placement(this_box,this_im,is,frag_start,frag_total,frag_order,frag_placed, &
    this_lambda, e_total,attempt_prob, nrg_ring_frag_tot, cbmc_overlap, del_overlap)
!**************************************************************************************************

  USE IO_Utilities


  INTEGER, INTENT(IN) :: this_im, is, frag_start, frag_total, this_box
  INTEGER, DIMENSION(1:nfragments(is)), INTENT(IN) :: frag_order

  INTEGER, DIMENSION(1:nfragments(is)), INTENT(INOUT) :: frag_placed

  REAL(DP), INTENT(INOUT) :: e_total, attempt_prob, nrg_ring_frag_tot
  REAL(DP), INTENT(IN) :: this_lambda
  
  LOGICAL, INTENT(INOUT) :: cbmc_overlap, del_overlap

  ! local variable

  INTEGER :: i, ifrag, j, this_atom, frag_connect, total_connect, anchor_ifrag, is1, im, ia
  INTEGER :: anchor_frag_connect, atom_ifrag, atom_frag_connect, ii, trial
  INTEGER :: frag_type, dumcount

  INTEGER, DIMENSION(:), ALLOCATABLE :: counted, connection, atom_id

  INTEGER :: ispecies, jmol, k

  INTEGER :: total_frags, this_fragment, nfrag_atoms
  REAL(DP) :: x_this,y_this,z_this, vec1(3), vec2(3), aligner_ifrag(3,3)
  REAL(DP) :: hanger_ifrag(3,3), aligner_frag_connect(3,3), hanger_frag_connect(3,3)

  REAL(DP) :: tempx, tempy, tempz, theta, e_dihed
  REAL(DP) :: weight(MAX(kappa_ins,kappa_rot,kappa_dih)), nrg(MAX(kappa_ins,kappa_rot,kappa_dih))
  REAL(DP) :: E_intra_qq, E_intra_vdw, prob_pick
  REAL(DP) :: e_prev, temp_var, E_ang, E_inter_vdw, E_inter_qq
  REAL(DP) :: nrg_kBT, p_acc, nrg_intra_vdw, nrg_intra_qq, nrg_inter_vdw, nrg_inter_qq
  REAL(DP) :: trial_weight
  REAL(DP) :: nrg_ring_frag, nrg_dihed(MAX(kappa_ins,kappa_rot,kappa_dih))

  LOGICAL :: overlap, overlap_trial(MAX(kappa_ins,kappa_rot,kappa_dih))

  CHARACTER :: this_file*120, element*1 

  TYPE(Atom_Class), ALLOCATABLE, DIMENSION(:) :: config_list
  TYPE(Atom_Class), ALLOCATABLE, DIMENSION(:,:) :: config_temp_list




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



     IF ( del_Flag) THEN

        ! use the existing coordinates for the second fragment

        DO j = 1, frag_list(ifrag,is)%natoms
           this_atom = frag_list(ifrag,is)%atoms(j)

           atom_list(this_atom,this_im,is)%exist = .TRUE.


           config_list(this_atom)%rxp = atom_list(this_atom,this_im,is)%rxp
           config_list(this_atom)%ryp = atom_list(this_atom,this_im,is)%ryp
           config_list(this_atom)%rzp = atom_list(this_atom,this_im,is)%rzp
           

        END DO

        IF (frag_list(ifrag,is)%ring) THEN
           ! obtain the fragment intramolecular + dihedral angle energy
           CALL Compute_Ring_Fragment_Energy(ifrag,this_im,is,this_box,nrg_ring_frag)
           nrg_ring_frag_tot = nrg_ring_frag_tot + nrg_ring_frag
           
        END IF

     ELSE
        
        ! let us read in the coordinates from a reservoir file
        
        ! Read in the configuration from the reservoir
      
        
        total_frags = frag_list(ifrag,is)%nconfig
        frag_type = frag_list(ifrag,is)%type
        
        ! Pick one from the reservoir
        
        this_fragment = INT(rranf() * total_frags) + 1
        
        ! Read in the coordinates
        
        DO j = 1, frag_list(ifrag,is)%natoms
           
           this_atom = frag_list(ifrag,is)%atoms(j)
           
           atom_list(this_atom,this_im,is)%exist = .TRUE.

           
           config_list(this_atom)%rxp = frag_coords(j,this_fragment,frag_type)%rxp
           config_list(this_atom)%ryp = frag_coords(j,this_fragment,frag_type)%ryp
           config_list(this_atom)%rzp = frag_coords(j,this_fragment,frag_type)%rzp
           
        END DO
     
        IF (frag_list(ifrag,is)%ring) THEN
           ! access the intramolecular + dihedral angle energy for this fragment
           nrg_ring_frag = nrg_frag(this_fragment,frag_type)
           nrg_ring_frag_tot = nrg_ring_frag_tot + nrg_ring_frag
        END IF
        
        CLOSE(UNIT=11)

     END IF
     
     ! Note that there has to be only one connection of ifrag that is already placed
     
     ! Let us find out which that fragment is
     

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
        err_msg(1) = 'More than one connections of' // TRIM(Int_To_String(ifrag)) // 'exist'
        CALL Clean_Abort(err_msg,'Fragment_Placement')
     END IF
     
     ! If here then only one connection found

     frag_connect = connection(1)

    
     ! Note that frag_connect already has anchor of ifrag placed both for fixed and variable
     ! bond length cases so all we have to do is obtain the coordinates of these atoms in the
     ! configuration and align it to their coordinates in the simulation box.

     
     ! find anchor atom ids for both the fragments

     CALL Get_Common_Fragment_Atoms(is,ifrag,frag_connect,anchor_ifrag,anchor_frag_connect)
     nfrag_atoms = 0
     atom_id(:) = 0
     
     DO j = 1, frag_list(ifrag,is)%natoms
       
        this_atom = frag_list(ifrag,is)%atoms(j)
        
        IF ( (this_atom /= anchor_ifrag) .AND. (this_atom /= anchor_frag_connect)) THEN
           nfrag_atoms = nfrag_atoms + 1
           atom_id(nfrag_atoms) = this_atom
        END IF
        
     END DO
     
!     anchor_ifrag = frag_list(ifrag,is)%anchor(1)
!     anchor_frag_connect = frag_list(frag_connect,is)%anchor(1)
!     write(*,*) anchor_ifrag, anchor_frag_connect, frag_connect
!     write(*,*) anchor_ifrag, anchor_frag_connect, frag_connect
     ! find one atom of ifrag and frag_connect that will be used for generating xy plane

     DO j = 1, frag_list(ifrag,is)%natoms

        this_atom = frag_list(ifrag,is)%atoms(j)

        IF ( (this_atom /= anchor_ifrag) .AND. (this_atom /= anchor_frag_connect)) EXIT

     END DO

     atom_ifrag = this_atom

     ! Similarly for frag_connect

     DO j = 1, frag_list(frag_connect,is)%natoms

        this_atom = frag_list(frag_connect,is)%atoms(j)
     
        IF ( (this_atom /= anchor_ifrag) .AND. (this_atom /= anchor_frag_connect)) EXIT

     END DO

     
     atom_frag_connect = this_atom

     ! Now use three atoms to obtain aligner and hanger matrix for the two fragments
     ! atom1 == origin
     ! atom2 == id of the atom along which +ve x - axis is aligned
     ! atom3 == helps to obtain y-axis so that atom1-atom2-atom3 define xy plane

     ! vec1 == r_atom2 - r_atom1
     ! vec2 == r_atom3 - r_atom1

     ! form these vectors from the configuration read in from reservoir

     ! for ifrag

     vec1(1) = config_list(anchor_ifrag)%rxp - config_list(anchor_frag_connect)%rxp
     vec1(2) = config_list(anchor_ifrag)%ryp - config_list(anchor_frag_connect)%ryp
     vec1(3) = config_list(anchor_ifrag)%rzp - config_list(anchor_frag_connect)%rzp


     vec2(1) = config_list(atom_ifrag)%rxp - config_list(anchor_frag_connect)%rxp
     vec2(2) = config_list(atom_ifrag)%ryp - config_list(anchor_frag_connect)%ryp
     vec2(3) = config_list(atom_ifrag)%rzp - config_list(anchor_frag_connect)%rzp

        
     CALL Get_Aligner_Hanger(vec1, vec2, aligner_ifrag,hanger_ifrag)
     
     ! Calculate this only for inserting a molecule
     
     IF ( .NOT. del_Flag) THEN
        
        
        ! for frag_connect
        
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
        
     ! Apply aligner of ifrag and then hanger of frag_connect to join the two fragments.
     
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
        
        config_list(this_atom)%rxp = tempx * aligner_ifrag(1,1) + tempy * aligner_ifrag(1,2) + &
             tempz * aligner_ifrag(1,3)
        config_list(this_atom)%ryp = tempx * aligner_ifrag(2,1) + tempy * aligner_ifrag(2,2) + &
             tempz * aligner_ifrag(2,3)
        config_list(this_atom)%rzp = tempx * aligner_ifrag(3,1) + tempy * aligner_ifrag(3,2) + &
             tempz * aligner_ifrag(3,3)
        
     END DO
     
     ! At this point, we can generate kappa positions of the fragment as two anchor positions
     ! are aligned. This is, in effect, equivalent to rotating the non-anchor atoms around
     ! the x-axis. Note that the coordinates of the anchoring atoms do not change due to
     ! this rotation. For the deletion move, the first trial must be the one corresponding
     ! to the actual coordinates, hence there should be no rotation about x-axis.

     IF ( del_Flag ) THEN

        theta = 0.0_DP

        ! also note that we will transform the position based on hanger_ifrag so that
        ! the original positions are recovered

        hanger_frag_connect(:,:) = hanger_ifrag(:,:)

     ELSE

        ! choose a random theta 
        
        theta = twopi * rranf()
        
     END IF

     ! initialize the energies

     nrg(:) = 0.0_DP
     nrg_dihed(:) = 0.0_DP
     ii = 1

     DO 
!     DO ii = 1, kappa

        config_temp_list(:,ii)%rxp  = config_list(:)%rxp
        config_temp_list(:,ii)%ryp  = 0.0_DP
        config_temp_list(:,ii)%rzp  = 0.0_DP
        
        DO j = 1, frag_list(ifrag,is)%natoms
           
           this_atom = frag_list(ifrag,is)%atoms(j)
           
           tempx = config_list(this_atom)%rxp
           tempy = config_list(this_atom)%ryp
           tempz = config_list(this_atom)%rzp
           
           config_temp_list(this_atom,ii)%ryp =  DCOS(theta) * tempy + DSIN(theta) * tempz
           config_temp_list(this_atom,ii)%rzp = -DSIN(theta) * tempy + DCOS(theta) * tempz
           
        END DO

        DO j = 1, frag_list(ifrag,is)%natoms
           
           this_atom = frag_list(ifrag,is)%atoms(j)
           
           tempx = config_temp_list(this_atom,ii)%rxp
           tempy = config_temp_list(this_atom,ii)%ryp
           tempz = config_temp_list(this_atom,ii)%rzp
           
           config_temp_list(this_atom,ii)%rxp = tempx * hanger_frag_connect(1,1) + &
                tempy * hanger_frag_connect(1,2) + &
                tempz * hanger_frag_connect(1,3)
           
           config_temp_list(this_atom,ii)%ryp = tempx * hanger_frag_connect(2,1) + &
                tempy * hanger_frag_connect(2,2) + &
                tempz * hanger_frag_connect(2,3) 
           
           config_temp_list(this_atom,ii)%rzp = tempx * hanger_frag_connect(3,1) + &
                tempy * hanger_frag_connect(3,2) + &
                tempz * hanger_frag_connect(3,3)


           IF ( this_atom /= anchor_ifrag) THEN
              IF ( this_atom /= anchor_frag_connect) THEN
                 
                 config_temp_list(this_atom,ii)%rxp = config_temp_list(this_atom,ii)%rxp + &
                      atom_list(anchor_frag_connect,this_im,is)%rxp
                 
                 config_temp_list(this_atom,ii)%ryp = config_temp_list(this_atom,ii)%ryp + &
                      atom_list(anchor_frag_connect,this_im,is)%ryp
                 
                 config_temp_list(this_atom,ii)%rzp = config_temp_list(this_atom,ii)%rzp + &
                      atom_list(anchor_frag_connect,this_im,is)%rzp
                 
!                 write(*,*) this_atom
                 atom_list(this_atom,this_im,is)%rxp = config_temp_list(this_atom,ii)%rxp 

                 
                 atom_list(this_atom,this_im,is)%ryp = config_temp_list(this_atom,ii)%ryp 
                 
                 
                 atom_list(this_atom,this_im,is)%rzp = config_temp_list(this_atom,ii)%rzp  
!                 write(*,*) atom_list(this_atom,this_im,is)%rxp, atom_list(this_atom,this_im,is)%ryp, &
!                      atom_list(this_atom,this_im,is)%rzp
                 
              END IF
           END IF
           
        END DO




           IF( ii == kappa_dih ) EXIT
           ii = ii + 1

           theta = theta + twopi / REAL(kappa_dih,DP)


     END DO

     DO ii = 1, kappa_dih

        DO j = 1, frag_list(ifrag,is)%natoms
           
           this_atom = frag_list(ifrag,is)%atoms(j)
           
           IF (this_atom /= anchor_ifrag) THEN
              IF (this_atom /= anchor_frag_connect) THEN
!                 write(*,*) this_atom
                 atom_list(this_atom,this_im,is)%rxp = config_temp_list(this_atom,ii)%rxp 
                 atom_list(this_atom,this_im,is)%ryp = config_temp_list(this_atom,ii)%ryp 
                 atom_list(this_atom,this_im,is)%rzp = config_temp_list(this_atom,ii)%rzp  
              END IF
           END IF

        END DO
 
        CALL Get_COM(this_im,is)
        CALL Compute_Max_COM_Distance(this_im,is)

        ! if it's a slit pore simulation first check to ensure that the atoms
        ! are inside the pore


        overlap = .FALSE.
           DO j = 1, nfrag_atoms
              atom_list(atom_id(j),this_im,is)%exist = .FALSE.
           END DO
           
           nrg_intra_vdw = 0.0_DP
           nrg_intra_qq = 0.0_DP
           nrg_inter_vdw = 0.0_DP
           nrg_inter_qq = 0.0_DP
           
           DO j = 1, nfrag_atoms
              
              atom_list(atom_id(j),this_im,is)%exist = .TRUE.
              
              
              
              CALL Compute_Atom_Nonbond_Energy(atom_id(j), this_im, is, E_intra_vdw, &
                   E_inter_vdw, E_intra_qq,E_inter_qq,overlap)
              
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

       
!        CALL Compute_Molecule_Nonbond_Intra_Energy(this_im,is,E_intra_vdw,E_intra_qq)

!        CALL Compute_Molecule_Nonbond_Inter_Energy(this_im,is,E_inter_vdw, E_inter_qq, overlap)

        ! subtract off the energy upto this point. Note that e_prev is the energy
        ! of the molecule upto this point. The substraction also helps to avoid
        ! overflow.


!        write(*,*) nrg(ii), e_dihed, E_intra_vdw, e_prev

        IF (overlap) THEN

           weight(ii) = 0.0_DP
!           overlap_trial(ii) = .TRUE.
           
        ELSE
           
           CALL Compute_Molecule_Dihedral_Energy(this_im,is,e_dihed)

           nrg_dihed(ii) = e_dihed

           nrg(ii) = e_dihed + nrg_intra_vdw + nrg_intra_qq  + &          
                nrg_inter_vdw + nrg_inter_qq - e_prev

           IF (frag_list(ifrag,is)%ring) THEN
              ! subtract off the biasing energy used to sample intramolecular DOFs
              nrg(ii) = nrg(ii) - nrg_ring_frag
           END IF
           
           nrg_kBT = beta(this_box) * nrg(ii)
           weight(ii) = DEXP(-nrg_kBT)

           
        END IF

        IF ( ii > 1 ) weight(ii) = weight(ii-1) + weight(ii)

        
!        WRITE(*,*) 'weight',overlap, ii, nrg_kBT
!        WRITE(*,*) 'energy', e_prev, nrg(ii)
           

     END DO

     IF (weight(kappa_dih) == 0.0_DP) THEN
        cbmc_overlap = .TRUE.
        IF (del_flag) THEN
           write(*,*) 'overlap detected in the deletion attempt'

        END IF

        RETURN
     END IF

!     weight(:) = weight(:)/weight(kappa)

!     write(*,*) 'weight', weight

     IF (del_Flag) THEN
        trial = 1
     ELSE

        prob_pick = rranf() * weight(kappa_dih)

        DO ii = 1, kappa_dih
!!$           
           IF ( prob_pick <= weight(ii) ) EXIT
!!$        
        END DO
!!$        
        trial = ii
!!$        
     END IF

     IF (trial == 1) THEN
        attempt_prob = attempt_prob * weight(1)/weight(kappa_dih)
     ELSE
        attempt_prob = attempt_prob * (weight(trial) - weight(trial-1))/weight(kappa_dih)
          
     END IF


     IF (attempt_prob == 0.0_DP) THEN
        IF (del_FLAG) THEN
            write(*,*) 'old configuration has zero weight'
            write(*,*) 'aborting'
           del_overlap = .TRUE.
           cbmc_overlap = .TRUE.


           RETURN
           
        ELSE

           cbmc_overlap = .TRUE.
           RETURN
           
        END IF
           
     END IF
     
     
     ! Give the coordinates of this conformation to atom_list
     
     DO j = 1, frag_list(ifrag,is)%natoms
        
        this_atom = frag_list(ifrag,is)%atoms(j)
        
        IF (this_atom /= anchor_ifrag) THEN
           IF (this_atom /= anchor_frag_connect) THEN
!              write(*,*) this_atom
              atom_list(this_atom,this_im,is)%rxp = config_temp_list(this_atom,trial)%rxp
              atom_list(this_atom,this_im,is)%ryp = config_temp_list(this_atom,trial)%ryp
              atom_list(this_atom,this_im,is)%rzp = config_temp_list(this_atom,trial)%rzp
           END IF
        END IF
     END DO
     
     ! also store the total energy upto this point
     
     e_prev =  nrg_dihed(trial)
     
     ! mark this fragment as placed
     frag_placed(ifrag) = 1
    

  END DO

!!$  DO i = 1, frag_list(1,is)%natoms
!!$     this_atom = frag_list(1,is)%atoms(i)
!!$     write(13,*) nonbond_list(this_atom,is)%element, atom_list(this_atom,this_im,is)%rxp, &
!!$          atom_list(this_atom,this_im,is)%ryp, atom_list(this_atom,this_im,is)%rzp
!!$  END DO

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
   
   INTEGER :: i
   
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
      
      atom_list(this_atom,alive,is)%rxp = frag_coords(i,this_fragment,frag_type)%rxp
      atom_list(this_atom,alive,is)%ryp = frag_coords(i,this_fragment,frag_type)%ryp
      atom_list(this_atom,alive,is)%rzp = frag_coords(i,this_fragment,frag_type)%rzp
      
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

 FUNCTION cass_EXP(exp_arg)

   USE Type_Definitions, ONLY : DP

   USE Run_Variables, ONLY : max_KBT

   IMPLICIT NONE
   
   REAL(DP) :: cass_EXP
   
   REAL(DP) :: exp_arg
   
   
   
   IF (exp_arg > max_kBT) THEN
      
      cass_EXP = 0.0_DP
      
   ELSE
      
      cass_EXP = DEXP(-exp_arg)
      
   END IF
   
 END FUNCTION cass_EXP
