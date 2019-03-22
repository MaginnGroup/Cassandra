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

SUBROUTINE Reaction

  !*****************************************************************************
  !
  ! This subroutine performs reactions in the RxMC ensemble.
  ! 
  !    A + B <--> AB
  !
  ! This algorithm is executed via the following steps
  !
  !    Step 1) Select a reaction with uniform probability
  !    Step 2) Select a direction
  !    Step 3) Choose molecules to delete
  !    Step 4) Calculate the energy of each reactant before the move
  !    Step 5) Delete the reactants
  !    Step 6) Insert products
  !    Step 7) Compute energies
  !    Step 8) Accept or reject the move
  !
  ! Called by
  !
  !   rxmc_driver
  !
  ! Revision history
  !
  !   06/30/15 : Version 1.2 Release
  !*****************************************************************************

  USE Global_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE IO_Utilities
  USE Pair_Nrg_Routines
  USE Rotation_Routines
  USE Fragment_Growth

  IMPLICIT NONE

  ! Local Variables
  INTEGER :: i, j, k           ! generic indices
  INTEGER :: ifrag
  INTEGER :: irxn, idir, irxn2 ! rxn indices
  INTEGER :: ibox              ! box index
  INTEGER :: is, is2, is_ref       ! species indices
  INTEGER :: im, im2            ! molecule indices
  INTEGER :: lm, lm2, lm_ref           ! molecule locate
  INTEGER :: ia, ia1m1, ia2m1, ia3m1, ia1m2, ia2m2, ia3m2  ! atom indices

  REAL(DP) :: dx, dy, dz

  REAL(DP) :: dE
  REAL(DP) :: E_bond, E_angle, E_dihed, E_improper
  REAL(DP) :: E_intra_vdw, E_intra_qq, E_periodic_qq
  REAL(DP) :: E_inter_vdw, E_inter_qq
  REAL(DP) :: E_pair_vdw, E_pair_qq
  REAL(DP) :: E_self, E_reciprocal, E_lrc
  REAL(DP) :: dE_bond, dE_angle, dE_dihed, dE_improper
  REAL(DP) :: dE_intra, dE_intra_vdw, dE_intra_qq
  REAL(DP) :: dE_inter_vdw, dE_inter_qq
  REAL(DP) :: dE_self, dE_reciprocal, dE_lrc
  INTEGER :: kappa_tot
  INTEGER, ALLOCATABLE :: frag_order(:)

  REAL(DP) :: ln_pacc  ! acceptance criteria

  LOGICAL :: inter_overlap, intra_overlap, cbmc_overlap
  LOGICAL :: accept_or_reject

  ! Rxn Variables
  INTEGER :: ireact, ireact2, i_start, i_end ! reactant indices
  LOGICAL :: rxn_flag, repeat_flag
  REAL(DP) :: x, xnew, xmax, xmin  ! extent of reaction
  INTEGER :: ix, ixnew ! extent of reaction index
  INTEGER :: noptions_forward, noptions_reverse 
  INTEGER :: delta_nu, nu
  CHARACTER(3) :: dir(-1:1) = (/ "rev", "000", "fwd" /)
  REAL(DP) :: P_molec, P_bias, P_temp
  INTEGER :: n_del, n_ins
  INTEGER, ALLOCATABLE :: lm_del(:), is_del(:)
  INTEGER, ALLOCATABLE :: lm_ins(:), is_ins(:)

  ! ideal gas variables
  INTEGER :: rand_igas, frag_type
  REAL(DP) :: E_intra_vdw_igas, E_intra_qq_igas

  ! Initialize variables
  inter_overlap = .FALSE.
  intra_overlap = .FALSE.
  cbmc_overlap = .FALSE.
  P_bias = 1.0_DP
  P_temp = 1.0_DP
  n_del = 0
  n_ins = 0

  !*****************************************************************************
  ! Step 1) Select a reaction
  !*****************************************************************************
  irxn = INT(rranf() * nrxns) + 1

  ibox = rxn_info(irxn)%which_box

  ! Increment counters
  tot_trials(ibox) = tot_trials(ibox) + 1

  IF (debug_log) THEN
    WRITE(logunit,'(A6,X,5X)',ADVANCE='NO') ADJUSTL(Int_To_String(irxn))
  END IF

  ALLOCATE(lm_del(MAX(rxn_info(irxn)%nreactants,rxn_info(irxn)%nproducts)))
  ALLOCATE(is_del(MAX(rxn_info(irxn)%nreactants,rxn_info(irxn)%nproducts)))
  ALLOCATE(lm_ins(MAX(rxn_info(irxn)%nreactants,rxn_info(irxn)%nproducts)))
  ALLOCATE(is_ins(MAX(rxn_info(irxn)%nreactants,rxn_info(irxn)%nproducts)))
  !*****************************************************************************
  ! Step 2) Select a direction
  !*****************************************************************************

  idir = INT(rranf() * 2) * 2 - 1
  IF (debug_log) THEN
    WRITE(logunit,'(X,A3,X,I3)',ADVANCE='NO') dir(idir), ibox
  END IF

  ! Check if there are enough reactants/products for the forward move
  rxn_flag = .TRUE.
  P_molec = 1.0_DP
  DO is = 1, nspecies
    nu = -idir* rxn_info(irxn)%stoichiometry(is)
    ! For idir=1, we need reactants. The stoichiometry of reactants is
    ! negative. In this case, -idir = -1 turns the stoichiometry of 
    ! reactants positive and the stoichiometry of products negative.
    ! 
    ! For idir=-1, we need products. -idir = 1
    IF (nu > 0) THEN
      IF (nmols(is,ibox) < nu) THEN
        ! not enough of 'is' in ibox to react
        accept = .FALSE.
        IF (verbose_log) THEN
          IF (debug_log) THEN
            WRITE(logunit,'(X,L8,X,A)') accept, 'too_few_mol'
          ELSE
            WRITE(logunit,'(X,I9,X,A10,X,5X,X,A3,X,I3,X,L8,X,A)') i_mcstep, &
               'rxn_' // TRIM(Int_To_String(irxn)), dir(idir), ibox, accept, &
               'too_few_mol'
          END IF
        END IF
        RETURN
      ELSE
        DO i = 0, nu - 1 
          P_molec = P_molec / REAL(nmols(is,ibox) - i,DP)
        END DO
      END IF
    ELSE IF (nu < 0) THEN
      IF (nmols(is,0) < - nu) THEN
        ! cannot place more 'is' in the system
        accept = .FALSE.
        IF (verbose_log) THEN
          IF (debug_log) THEN
            WRITE(logunit,'(X,L8,X,A)') accept, 'too_many_mol'
          ELSE
            WRITE(logunit,'(X,I9,X,A10,X,5X,X,A3,X,I3,X,L8,X,A)') i_mcstep, &
               'rxn_' // TRIM(Int_To_String(irxn)), dir(idir), ibox, accept, &
               'too_many_mol'
          END IF
        END IF
        RETURN
      ELSE
        DO i = 0, - nu - 1 
          P_molec = P_molec * REAL(nmols(is,ibox) - nu - i,DP)
        END DO
      END IF
    END IF
  END DO


  !Increment counter
  rxn_status(irxn)%ntrials(0,idir) = rxn_status(irxn)%ntrials(0,idir) + 1

  !*****************************************************************************
  ! Step 3) Select reactants/products to delete
  !*****************************************************************************
  ! idir ==  1, need to delete reactants
  ! idir == -1, need to delete products
  i_start = 1 + (1 - idir) / 2 * rxn_info(irxn)%nreactants
  i_end = rxn_info(irxn)%nreactants + (1 - idir) / 2 * rxn_info(irxn)%nproducts

  DO ireact = i_start, i_end
    is = rxn_info(irxn)%species(ireact)
    ! may need >1 unique molecules of 'is'
    repeat_flag = .TRUE.
    DO WHILE (repeat_flag)
      ! assume we'll get a unique molec of 'is'
      repeat_flag = .FALSE.
      ! Pick molecules with uniform probability
      im = INT(rranf() * nmols(is,ibox)) + 1 ! molecule INDEX
      ! Find the index of the im-th molecule of species 'is'
      lm = locate(im,is,ibox) ! molecule LOCATE
      ! Check that molecule 'lm' has not been previously chosen
      DO irxn2 = 1, nrxns
        DO ireact2 = 1, &
                    rxn_info(irxn2)%nreactants + rxn_info(irxn2)%nproducts
          IF (irxn == irxn2 .AND. ireact == ireact2) CYCLE
          IF (rxn_info(irxn2)%species(ireact2) /= is) CYCLE
          IF (rxn_status(irxn2)%molecule(ireact2) == lm) THEN
            ! we didn't get a unique molec
            repeat_flag = .TRUE.
          END IF
        END DO
      END DO
    END DO
    rxn_status(irxn)%molecule(ireact) = lm

    ! Build list of deleted molecules for nonbond energy calc
    n_del = n_del + 1
    lm_del(n_del) = lm
    is_del(n_del) = is
  END DO

  !*****************************************************************************
  ! Step 4) Calculate the energy of each reactant/product before it's deleted
  !*****************************************************************************
  dE_intra     = 0.0_DP
  dE_bond      = 0.0_DP
  dE_angle     = 0.0_DP
  dE_dihed     = 0.0_DP
  dE_improper  = 0.0_DP
  dE_intra_vdw = 0.0_DP
  dE_intra_qq  = 0.0_DP
  dE_self      = 0.0_DP
  DO ireact = i_start, i_end
    is = rxn_info(irxn)%species(ireact)
    lm = rxn_status(irxn)%molecule(ireact)

    ! Intramolecular energy
    CALL Compute_Molecule_Bond_Energy(lm,is,E_bond)
    dE_bond = dE_bond - E_bond

    CALL Compute_Molecule_Angle_Energy(lm,is,E_angle)
    dE_angle = dE_angle - E_angle

    CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihed)
    dE_dihed = dE_dihed - E_dihed

    CALL Compute_Molecule_Improper_Energy(lm,is,E_improper)
    dE_improper = dE_improper - E_improper

    CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is,E_intra_vdw, &
         E_intra_qq,E_periodic_qq,intra_overlap)
    dE_intra_vdw = dE_intra_vdw - E_intra_vdw
    dE_intra_qq = dE_intra_qq - E_intra_qq
    dE_inter_qq = dE_inter_qq - E_periodic_qq

    ! Ewald self energies
    IF (int_charge_sum_style(ibox) == charge_ewald .AND. has_charge(is)) THEN
      CALL Compute_Molecule_Self_Energy(lm,is,ibox,E_self)
      dE_Self = dE_Self - E_self
    END IF
  END DO

  ! Intermolecular energy, of reactants and shell molecules
  E_inter_vdw = 0.0_DP
  E_inter_qq = 0.0_DP
  IF (l_pair_nrg) THEN
    CALL Store_Molecule_Pair_Interaction_Arrays(lm,is,&
         E_inter_vdw,E_inter_qq,n_del,lm_del,is_del)
  ELSE
    CALL Compute_MoleculeCollection_Nonbond_Inter_Energy(n_del,lm_del,is_del,&
         E_inter_vdw,E_inter_qq,inter_overlap)
    ! inter_overlap /= .TRUE. since this is the existing configuration
  END IF
  dE_inter_vdw = - E_inter_vdw
  dE_inter_qq  = - E_inter_qq

  ! Ewald reciprocal energy
  dE_reciprocal = - energy(ibox)%reciprocal

  ! Long range correction
  dE_lrc = - energy(ibox)%lrc

  !*************************************************************************
  ! Step 5) Delete each reactant/product
  !*************************************************************************
  DO ireact = i_start, i_end
    is = rxn_info(irxn)%species(ireact)
    lm = rxn_status(irxn)%molecule(ireact)

    ! find the insertion_bias for the reverse move
    IF (rxn_info(irxn)%insertion_method(ireact) == int_random) THEN
      del_flag = .TRUE.     ! do not change the coordinates of 'lm'
      cbmc_overlap = .FALSE.
      CALL CBMC_Insert_Rotate(lm,is,ibox,P_temp,cbmc_overlap)
      P_bias = P_bias / P_temp

      IF (cbmc_overlap) THEN
        err_msg = ""
        err_msg(1) = "Existing configuration tripped the cbmc_overlap flag"
        err_msg(2) = "molecule " // TRIM(Int_To_String(lm)) // " of species " // TRIM(Int_To_String(is))
        CALL Clean_Abort(err_msg, "Reaction")
      END IF
    END IF

    ! delete the molecule

    ! need INDEX that corresponds to the saved LOCATE
    DO im = 1,nmols(is,ibox)
       IF ( locate(im,is,ibox) == lm ) EXIT   
    END DO
    ! remove the deleted LOCATE from the array
    IF (im < nmols(is,ibox)) THEN
      DO im2 = im + 1, nmols(is,ibox)
        locate(im2-1,is,ibox) = locate(im2,is,ibox)
      END DO
    END IF
    locate(nmols(is,ibox),is,ibox) = 0

    ! move LOCATE to list of unused LOCATES
    nmols(is,0) = nmols(is,0) + 1
    locate(nmols(is,0),is,0) = lm

    ! Turn the molecule off
    molecule_list(lm,is)%live = .FALSE.
    atom_list(:,lm,is)%exist = .FALSE.
    nmols(is,ibox) = nmols(is,ibox) - 1
  END DO

  !*************************************************************************
  ! Step 6) Insert product/reactant molecules
  !*************************************************************************
  ! idir ==  1, need to insert products
  ! idir == -1, need to insert reactants
  i_start = 1 + (1 + idir) / 2 * rxn_info(irxn)%nreactants
  i_end = rxn_info(irxn)%nreactants + (1 + idir) / 2 * rxn_info(irxn)%nproducts

  ! Get the LOCATES of all molecules being inserted
  DO ireact = i_start, i_end
    is = rxn_info(irxn)%species(ireact)
    nmols(is,ibox) = nmols(is,ibox) + 1
    im = nmols(is,ibox)
    locate(im,is,ibox) = locate(nmols(is,0),is,0)
    lm = locate(im,is,ibox)
    rxn_status(irxn)%molecule(ireact) = lm

    ! Remove locate from Box 0
    locate(nmols(is,0),is,0) = 0
    nmols(is,0) = nmols(is,0) - 1

    ! Build list of inserted molecules for nonbond energy calc
    n_ins = n_ins + 1
    lm_ins(n_ins) = lm
    is_ins(n_ins) = is
  END DO

  ! Insert the molecules
  DO ireact = i_start, i_end
    is = rxn_info(irxn)%species(ireact)
    lm = rxn_status(irxn)%molecule(ireact)
    ! Turn the molecule on
    atom_list(:,lm,is)%exist = .TRUE.
    molecule_list(lm,is)%live = .TRUE.
    molecule_list(lm,is)%frac = 1.0_DP
    molecule_list(lm,is)%which_box = ibox

    ! Now we need to grab the molecule's conformation
    rand_igas = (rranf() * n_igas(is)) + 1

    atom_list(:,lm,is)%rxp = atom_list_igas(:,rand_igas,is)%rxp
    atom_list(:,lm,is)%ryp = atom_list_igas(:,rand_igas,is)%ryp
    atom_list(:,lm,is)%rzp = atom_list_igas(:,rand_igas,is)%rzp

    CALL Get_COM(lm,is)
    CALL Compute_Max_COM_Distance(lm,is)

    IF (rxn_info(irxn)%insertion_method(ireact) == int_random) THEN
      del_flag = .FALSE.     ! change the coordinates of 'lm'
      cbmc_overlap = .FALSE.
      CALL CBMC_Insert_Rotate(lm,is,ibox,P_temp,cbmc_overlap)
      P_bias = P_bias * P_temp

      IF (cbmc_overlap) EXIT
    ELSE
      ia1m1 = rxn_info(irxn)%insertion_link(ireact,1)
      ia2m1 = rxn_info(irxn)%insertion_link(ireact,2)
      ia3m1 = rxn_info(irxn)%insertion_link(ireact,3)

      ireact2 = rxn_info(irxn)%insertion_link(ireact,0)
      is2   = rxn_info(irxn)%species(ireact2)
      lm2   = rxn_status(irxn)%molecule(ireact2)
      ia1m2 = rxn_info(irxn)%insertion_link(ireact2,1)
      ia2m2 = rxn_info(irxn)%insertion_link(ireact2,2)
      ia3m2 = rxn_info(irxn)%insertion_link(ireact2,3)
      CALL Align_Molecules(is,lm,ia1m1,ia2m1,ia3m1,is2,lm2,ia1m2,ia2m2,ia3m2)
      CALL Get_COM(lm,is)
    END IF
  END DO

  IF (cbmc_overlap) THEN
     ! Reject the move
     DO ireact = 1, rxn_info(irxn)%nreactants + rxn_info(irxn)%nproducts
       is = rxn_info(irxn)%species(ireact)
       IF (idir * rxn_info(irxn)%stoichiometry(is) < 0) THEN
         ! Restore the LOCATE of the deleted molecule
         nmols(is,ibox) = nmols(is,ibox) + 1
         locate(nmols(is,ibox),is,ibox) = locate(nmols(is,0),is,0)
         locate(nmols(is,0),is,0) = 0
         nmols(is,0) = nmols(is,0) - 1

         ! Turn the molecule back on
         lm = locate(nmols(is,ibox),is,ibox)
         molecule_list(lm,is)%live = .TRUE.
         atom_list(:,lm,is)%exist = .TRUE.
       ELSE
         ! Remove the LOCATE of inserted molecule
         nmols(is,0) = nmols(is,0) + 1
         locate(nmols(is,0),is,0) = locate(nmols(is,ibox),is,ibox)
         locate(nmols(is,ibox),is,ibox) = 0
         nmols(is,ibox) = nmols(is,ibox) - 1

         ! Turn the molecule off
         lm = locate(nmols(is,0),is,0)
         molecule_list(lm,is)%live = .FALSE.
         atom_list(:,lm,is)%exist = .FALSE.
       END IF
     END DO

     accept = .FALSE.

     DO ireact = 1, rxn_info(irxn)%nreactants + rxn_info(irxn)%nproducts
       rxn_status(irxn)%molecule(ireact) = 0
     END DO

     IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)

     IF (debug_log) THEN
       WRITE(logunit,'(X,L8)') accept
       FLUSH(logunit)
     END IF
     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,5X,X,A3,X,I3,X,L8,X,A)') i_mcstep, 'rxn_' // TRIM(Int_To_String(irxn)), &
          dir(idir), ibox, accept, 'cbmc_overlap'
     END IF
     RETURN
  END IF

  !*************************************************************************
  ! Step 7) Compute energies of each product/reactant
  !*************************************************************************
  ! Intermolecular energy, of products and solvent molecules
  CALL Compute_MoleculeCollection_Nonbond_Inter_Energy(n_ins,lm_ins,is_ins, &
       E_inter_vdw,E_inter_qq,inter_overlap)

  IF (inter_overlap) THEN
    ! Reject the move
    DO ireact = 1, rxn_info(irxn)%nreactants + rxn_info(irxn)%nproducts
      is = rxn_info(irxn)%species(ireact)
      IF (idir * rxn_info(irxn)%stoichiometry(is) < 0) THEN
        ! Restore the LOCATE of the deleted molecule
        nmols(is,ibox) = nmols(is,ibox) + 1
        locate(nmols(is,ibox),is,ibox) = locate(nmols(is,0),is,0)
        locate(nmols(is,0),is,0) = 0
        nmols(is,0) = nmols(is,0) - 1

        ! Turn the molecule back on
        lm = locate(nmols(is,ibox),is,ibox)
        molecule_list(lm,is)%live = .TRUE.
        atom_list(:,lm,is)%exist = .TRUE.
      ELSE
        ! Remove the LOCATE of inserted molecule
        nmols(is,0) = nmols(is,0) + 1
        locate(nmols(is,0),is,0) = locate(nmols(is,ibox),is,ibox)
        locate(nmols(is,ibox),is,ibox) = 0
        nmols(is,ibox) = nmols(is,ibox) - 1

        ! Turn the molecule off
        lm = locate(nmols(is,0),is,0)
        molecule_list(lm,is)%live = .FALSE.
        atom_list(:,lm,is)%exist = .FALSE.
      END IF
    END DO

    accept = .FALSE.

    DO ireact = 1, rxn_info(irxn)%nreactants + rxn_info(irxn)%nproducts
      rxn_status(irxn)%molecule(ireact) = 0
    END DO

    IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)

    IF (debug_log) THEN
      WRITE(logunit,'(X,L8)') accept
      FLUSH(logunit)
    END IF
    IF (verbose_log) THEN
      WRITE(logunit,'(X,I9,X,A10,X,5X,X,A3,X,I3,X,L8,X,A)') i_mcstep, &
         'rxn_' // TRIM(Int_To_String(irxn)), dir(idir), ibox, accept, 'inter_overlap'
    END IF
    RETURN
  END IF

  dE_inter_vdw = dE_inter_vdw + E_inter_vdw
  dE_inter_qq  = dE_inter_qq  + E_inter_qq

  DO ireact = i_start, i_end
    is = rxn_info(irxn)%species(ireact)
    lm = rxn_status(irxn)%molecule(ireact)

    ! Intramolecular energy
    CALL Compute_Molecule_Bond_Energy(lm,is,E_bond)
    dE_bond = dE_bond + E_bond

    CALL Compute_Molecule_Angle_Energy(lm,is,E_angle)
    dE_angle = dE_angle + E_angle

    CALL Compute_Molecule_Dihedral_Energy(lm,is,E_dihed)
    dE_dihed = dE_dihed + E_dihed

    CALL Compute_Molecule_Improper_Energy(lm,is,E_improper)
    dE_improper = dE_improper + E_improper

    CALL Compute_Molecule_Nonbond_Intra_Energy(lm,is,E_intra_vdw, &
         E_intra_qq,E_periodic_qq,intra_overlap)
    ! intra_overlap /= .TRUE. b/c cbmc_overlap already checked for this
    dE_intra_vdw = dE_intra_vdw + E_intra_vdw
    dE_intra_qq = dE_intra_qq + E_intra_qq
    dE_inter_qq = dE_inter_qq + E_periodic_qq

    ! Ewald self energies
    IF ( (int_charge_sum_style(ibox) == charge_ewald) .AND. &
         has_charge(is) ) THEN
      CALL Compute_Molecule_Self_Energy(lm,is,ibox,E_self)
      dE_Self = dE_Self + E_self
    END IF
    
  END DO

  ! Eward reciprocal energy
  IF (int_charge_sum_style(ibox) == charge_ewald) THEN
    ! Compute cos_mol and sin_mol terms for the new molecules
    CALL Update_System_Ewald_Reciprocal_Energy(irxn,idir,ibox,int_rxn, &
         E_reciprocal)
    dE_reciprocal = dE_reciprocal + E_reciprocal
  END IF
  
  ! Long range correction
  IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN
    CALL Compute_LR_correction(ibox,E_lrc)
    dE_lrc = dE_lrc + E_lrc
  END IF

  dE_intra = dE_bond + dE_angle + dE_dihed + dE_improper
  dE = dE_intra + dE_intra_vdw + dE_intra_qq + dE_inter_vdw + dE_inter_qq + &
       dE_self + dE_reciprocal + dE_lrc

  !*************************************************************************
  ! Step 8) Accept or reject the move
  !*************************************************************************
  ! Define ln_pacc that will be used to accept or reject the move
  ! Acceptance criterion
  ln_pacc = beta(ibox) * (dE - dE_intra - dE_intra_vdw - dE_intra_qq)

  delta_nu = SUM(rxn_info(irxn)%stoichiometry(:))
  ln_pacc = ln_pacc + idir * beta(ibox) * rxn_info(irxn)%free_energy &
          - idir * delta_nu * DLOG(beta(ibox) &
          * rxn_info(irxn)%standard_pressure * box_list(ibox)%volume)

  ln_pacc = ln_pacc + DLOG(P_molec) + DLOG(P_bias)
 
  accept = accept_or_reject(ln_pacc)

  IF (debug_log) THEN
    WRITE(logunit,'(X,L8)') accept
    FLUSH(logunit)
  END IF

  IF (accept) THEN
    ! Accept move

    ! Increment counters
    rxn_status(irxn)%nsuccess(0,idir) = rxn_status(irxn)%nsuccess(0,idir) &
                                      + 1
    ! Update intermolecular energies for ibox
    energy(ibox)%total = energy(ibox)%total + dE
    energy(ibox)%intra = energy(ibox)%intra + dE_intra
    energy(ibox)%bond = energy(ibox)%bond + dE_bond
    energy(ibox)%angle = energy(ibox)%angle + dE_angle
    energy(ibox)%dihedral = energy(ibox)%dihedral + dE_dihed
    energy(ibox)%improper = energy(ibox)%improper + dE_improper
    energy(ibox)%intra_vdw = energy(ibox)%intra_vdw + dE_intra_vdw
    energy(ibox)%intra_q   = energy(ibox)%intra_q + dE_intra_qq
    energy(ibox)%inter_vdw = energy(ibox)%inter_vdw + dE_inter_vdw
    energy(ibox)%inter_q   = energy(ibox)%inter_q + dE_inter_qq

    IF (int_charge_sum_style(ibox) == charge_ewald) THEN
      energy(ibox)%self = energy(ibox)%self + dE_self
      energy(ibox)%reciprocal = E_reciprocal
    END IF

    IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN
      energy(ibox)%lrc = E_lrc
    END IF

  ELSE
    ! Reject the move
    DO ireact = 1, rxn_info(irxn)%nreactants + rxn_info(irxn)%nproducts
      is = rxn_info(irxn)%species(ireact)
      IF (idir * rxn_info(irxn)%stoichiometry(is) < 0) THEN
        ! Restore the LOCATE of the deleted molecule
        nmols(is,ibox) = nmols(is,ibox) + 1
        locate(nmols(is,ibox),is,ibox) = locate(nmols(is,0),is,0)
        locate(nmols(is,0),is,0) = 0
        nmols(is,0) = nmols(is,0) - 1

        ! Turn the molecule back on
        lm = locate(nmols(is,ibox),is,ibox)
        molecule_list(lm,is)%live = .TRUE.
        atom_list(:,lm,is)%exist = .TRUE.
      ELSE
        ! Remove the LOCATE of inserted molecule
        nmols(is,0) = nmols(is,0) + 1
        locate(nmols(is,0),is,0) = locate(nmols(is,ibox),is,ibox)
        locate(nmols(is,ibox),is,ibox) = 0
        nmols(is,ibox) = nmols(is,ibox) - 1

        ! Turn the molecule off
        lm = locate(nmols(is,0),is,0)
        molecule_list(lm,is)%live = .FALSE.
        atom_list(:,lm,is)%exist = .FALSE.
      END IF
    END DO

    ! Restore cos_sum and sin_sum. Note that these were changed when
    ! difference in reciprocal energies was computed
    IF (int_charge_sum_style(ibox) == charge_ewald) THEN
      cos_sum(:,ibox) = cos_sum_old(:,ibox)
      sin_sum(:,ibox) = sin_sum_old(:,ibox)
    END IF
  END IF

  IF (l_pair_nrg) DEALLOCATE(pair_vdw_temp,pair_qq_temp)

  DO ireact = 1, rxn_info(irxn)%nreactants + rxn_info(irxn)%nproducts
    rxn_status(irxn)%molecule(ireact) = 0
  END DO

  IF (verbose_log) THEN
    WRITE(logunit,'(X,I9,X,A10,X,5X,X,A3,X,I3,X,L8,X,9X,X,F9.3)') i_mcstep, &
       'rxn_' // TRIM(Int_To_String(irxn)), dir(idir), ibox, accept, ln_pacc
  END IF

CONTAINS

  SUBROUTINE CBMC_Insert_Rotate(lm,is,ibox,P_bias,cbmc_overlap)
    INTEGER, INTENT(IN) :: lm, is, ibox
    REAL(DP), INTENT(OUT) :: P_bias
    LOGICAL, INTENT(OUT) :: cbmc_overlap
    LOGICAL :: overlap
    INTEGER :: itrial, i_ins, i_rot, ntrials
    REAL(DP) :: xcom, ycom, zcom, x_scaled, y_scaled, z_scaled
    REAL(DP) :: nrg(kappa_ins*kappa_rot), weight(kappa_ins*kappa_rot)
    REAL(DP) :: nrg_kBT, rand_no
    TYPE(Atom_Class) :: rtrial(natoms(is),kappa_ins*kappa_rot)

    ! Initialize
    P_bias = 1.0_DP
    cbmc_overlap = .FALSE.
    cbmc_flag = .TRUE.
    itrial = 0
    ntrials = kappa_ins*kappa_rot
    nrg = 0.0_DP
    weight = 0.0_DP

    DO i_ins = 1, kappa_ins
      ! Select a trial location
      IF (del_flag .AND. i_ins == 1) THEN
        xcom = molecule_list(lm,is)%xcom
        ycom = molecule_list(lm,is)%ycom
        zcom = molecule_list(lm,is)%zcom
      ELSE 
        x_scaled = 0.5_DP - rranf()
        y_scaled = 0.5_DP - rranf()
        z_scaled = 0.5_DP - rranf()
        xcom = x_scaled * box_list(ibox)%length(1,1) &
             + y_scaled * box_list(ibox)%length(1,2) &
             + z_scaled * box_list(ibox)%length(1,3)
        ycom = x_scaled * box_list(ibox)%length(2,1) &
             + y_scaled * box_list(ibox)%length(2,2) &
             + z_scaled * box_list(ibox)%length(2,3)
        zcom = x_scaled * box_list(ibox)%length(3,1) &
             + y_scaled * box_list(ibox)%length(3,2) &
             + z_scaled * box_list(ibox)%length(3,3)
      END IF
      
      DO i_rot = 1, kappa_rot
        itrial = itrial + 1
        ! Select a trial rotation
        IF ((.NOT. del_flag .OR. itrial /= 1) .AND. natoms(is) > 1) THEN
          CALL Rotate_Molecule_Eulerian(lm,is)
        END IF
        DO ia = 1, natoms(is)
          rtrial(ia,itrial)%rxp = atom_list(ia,lm,is)%rxp &
                                - molecule_list(lm,is)%xcom + xcom
          rtrial(ia,itrial)%ryp = atom_list(ia,lm,is)%ryp &
                                - molecule_list(lm,is)%ycom + ycom
          rtrial(ia,itrial)%rzp = atom_list(ia,lm,is)%rzp &
                                - molecule_list(lm,is)%zcom + zcom
          atom_list(ia,lm,is)%rxp = rtrial(ia,itrial)%rxp
          atom_list(ia,lm,is)%ryp = rtrial(ia,itrial)%ryp
          atom_list(ia,lm,is)%rzp = rtrial(ia,itrial)%rzp
        END DO
        molecule_list(lm,is)%xcom = xcom
        molecule_list(lm,is)%ycom = ycom
        molecule_list(lm,is)%zcom = zcom
        
        overlap = .FALSE.
        CALL Compute_Molecule_Nonbond_Inter_Energy(lm,is,E_inter_vdw, &
                E_inter_qq,overlap)
        nrg(itrial) = nrg(itrial) + E_inter_vdw + E_inter_qq

        IF (overlap) THEN
          weight(itrial) = 0.0_DP
        ELSE
          nrg_kBT = beta(ibox) * nrg(itrial)
          IF (nrg_kBT >= max_kBT) THEN
            weight(itrial) = 0.0_DP
            IF (del_flag .AND. itrial == 1) THEN
              weight(itrial) = tiny_number
            END IF
          ELSE
            weight(itrial) = DEXP(-nrg_kBT)
          END IF
        END IF

        IF (itrial > 1) weight(itrial) = weight(itrial-1) + weight(itrial)
      END DO
    END DO

    ! Reject the move if the total weight is still zero
    IF (weight(ntrials) == 0.0_DP) THEN
      cbmc_overlap = .TRUE.
      cbmc_flag = .FALSE.
      RETURN
    END IF

    ! Select one of the trial coordinates
    IF (del_flag) THEN
      ! For a deletion move, we want the weight of the current position,
      ! which is stored in trial 1.
      itrial = 1
    ELSE
      ! Choose one from Golden sampling for an insertion move
      rand_no = rranf() * weight(ntrials)
    
      DO itrial = 1, ntrials
        IF (rand_no < weight(itrial)) EXIT
      END DO
    
      IF ( itrial == ntrials + 1 ) THEN
        ! None of the trials were picked. Could be due to the fact that all 
        ! the trials had a very small cumulative weight
        cbmc_overlap = .TRUE.
        cbmc_flag = .FALSE.
        RETURN
      END IF

    END IF

    ! Compute the weight of the selected trial coordinate
    IF (itrial == 1) THEN
      P_bias = P_bias * weight(1) / weight(ntrials)
    ELSE
      P_bias = P_bias * (weight(itrial) - weight(itrial-1)) / weight(ntrials)
    END IF
    ! So far P_bias only includes the probability of choosing the insertion 
    ! point from the collection of trial coordinates.
    ! We need to include the number of trial coordinates, kappa_ins*kappa_rot
    kappa_tot = kappa_ins
    IF (kappa_rot /= 0 ) THEN
      kappa_tot = kappa_tot * kappa_rot
    END IF
    P_bias = P_bias * REAL(kappa_tot, DP)
  
    ! We chose the ith trial coordinate for placement. Store the ith trial
    ! coordinates in the atom_list array. Note that for the deletion move,
    ! itrial=1 has the current coordinates of the fragment, so the molecule
    ! is not moved.
    DO ia = 1, natoms(is)
      atom_list(ia,lm,is)%rxp = rtrial(ia,itrial)%rxp
      atom_list(ia,lm,is)%ryp = rtrial(ia,itrial)%ryp
      atom_list(ia,lm,is)%rzp = rtrial(ia,itrial)%rzp
    END DO
    CALL Get_COM(lm,is)
    
    cbmc_flag = .FALSE.
  END SUBROUTINE CBMC_Insert_Rotate

END SUBROUTINE Reaction
