
  !Might be more useful later
  !implement as a simultaneous translation
  !d_comx = molecule_list(j_alive,js)%xcom - molecule_list(i_alive,is)%xcom
  !d_comy = molecule_list(j_alive,js)%ycom - molecule_list(i_alive,is)%ycom
  !d_comz = molecule_list(j_alive,js)%zcom - molecule_list(i_alive,is)%zcom

  !Translate COMs
  !molecule_list(i_alive,is)%xcom = molecule_list(i_alive,is)%xcom + dcomx
  !molecule_list(i_alive,is)%ycom = molecule_list(i_alive,is)%ycom + dcomy
  !molecule_list(i_alive,is)%zcom = molecule_list(i_alive,is)%zcom + dcomz

  !molecule_list(j_alive,js)%xcom = molecule_list(j_alive,js)%xcom - dcomx
  !molecule_list(j_alive,js)%ycom = molecule_list(j_alive,js)%ycom - dcomy
  !molecule_list(j_alive,js)%zcom = molecule_list(j_alive,js)%zcom - dcomz

  !Translate Individual Atoms
  !atom_list(:,i_alive,is)%rxp = atom_list(:,i_alive,is)%rxp + dcomx
  !atom_list(:,i_alive,is)%ryp = atom_list(:,i_alive,is)%ryp + dcomy
  !atom_list(:,i_alive,is)%rzp = atom_list(:,i_alive,is)%rzp + dcomz

  !atom_list(:,j_alive,js)%rxp = atom_list(:,j_alive,js)%rxp - dcomx
  !atom_list(:,j_alive,js)%ryp = atom_list(:,j_alive,js)%ryp - dcomy
  !atom_list(:,j_alive,js)%rzp = atom_list(:,j_alive,js)%rzp - dcomz


!______________________-__________________________________________________________
  !del_flag = .FALSE.
  !get_fragorder = .TRUE.
  !ALLOCATE(frag_order(nfragments(is)))
  !lambda_for_build = molecule_list(alive,is)%frac
  !CALL Build_Molecule(alive,is,box_in,frag_order, &
  !        lambda_for_build,ln_pseq,ln_pfor,nrg_ring_frag_in, &
  !        cbmc_overlap)

  !gets com of a molecule, stored in molecule_list(alive,is)%xcom, %ycom, %zcom
  CALL Get_COM(i_alive,is)

  !gets distance between is' COM and farthest away molecule
  !stored in molecule_list(alive,is)%max_dcom
  CALL Compute_Max_COM_Distance(i_alive,is)

  !If outside of the box, fold back into the box
  CALL Fold_Molecule(i_alive,is,box)

 CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is, &
             E_inter_vdw_i,E_inter_qq_i,inter_overlap_i)


  CALL Get_COM(j_alive, js)
  CALL Compute_Max_COM_Distance(j_alive, js)
  CALL Fold_Molecule(j_alive, js, box)

  CALL Compute_Molecule_Nonbond_Inter_Energy(alive,is, &
             E_inter_vdw_in,E_inter_qq_in,inter_overlap)

  IF (inter_overlap) THEN
     ! reject the swap

     ! restore the box coordinates
     CALL Revert_Old_Cartesian_Coordinates(alive,is)

     ! All atoms will not exist if inter_overlap was tripped before the
     ! last fragment was placed in Build_Molecule.
     ! Set exist to TRUE for all atoms and reset frac to 1
     atom_list(:,alive,is)%exist = .TRUE.
     molecule_list(alive,is)%frac = 1.0_DP

     IF (l_pair_nrg) THEN
        CALL Reset_Molecule_Pair_Interaction_Arrays(alive,is,box)
     END IF

     IF(ALLOCATED(cos_mol_old)) DEALLOCATE(cos_mol_old)
     IF(ALLOCATED(sin_mol_old)) DEALLOCATE(sin_mol_old)

     accept = .FALSE.

     IF (verbose_log) THEN
       WRITE(logunit,'(X,I9,X,A10,X,I5,X,I3,X,I1,A1,I1,X,L8,X,9X,X,A9)') &
             i_mcstep, 'swap' , alive, is, box, '>', box, accept, 'overlap'
     END IF

  ELSE
      !may be accepted
      !Where are these set?
      dE = dE + E_inter_vdw + E_inter_qq
  END IF

    !*****************************************************************************
    ! Step 5) Calculate the change in box_in's potential energy from inserting
    !         alive
    !*****************************************************************************
    ! If here then no overlap was detected. Calculate the rest of the energies
    CALL Compute_Molecule_Bond_Energy(alive,is,E_bond)
    CALL Compute_Molecule_Angle_Energy(alive,is,E_angle)
    CALL Compute_Molecule_Dihedral_Energy(alive,is,E_dihed)
    CALL Compute_Molecule_Improper_Energy(alive,is,E_improper)

    dE_in = dE_in + E_bond_in + E_angle_in + E_dihed_in + E_improper_in

    CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is, &
         E_intra_vdw_in,E_intra_qq_in,E_periodic_qq,intra_overlap)
    E_inter_qq_in = E_inter_qq_in + E_periodic_qq

    ! Already added E_inter to dE, so add E_periodic directly
    dE_in = dE_in + E_intra_vdw_in + E_intra_qq_in + E_periodic_qq

    call cpu_time(time0)

    IF (int_charge_style(box_in) == charge_coul .AND. has_charge(is)) THEN

       IF (int_charge_sum_style(box_in) == charge_ewald) THEN

            ! Note that this call will change cos_mol, sin_mol of alive and this
            ! will have to be restored below while computing the energy of box_out
            ! without molecule alive.
            CALL Update_System_Ewald_Reciprocal_Energy(alive,is,box_in, &
                 int_insertion,E_reciprocal_in)

            dE_in = dE_in + (E_reciprocal_in - energy(box_in)%ewald_reciprocal)
       END IF

       CALL Compute_Molecule_Self_Energy(alive,is,box_in,E_self_in)
       dE_in = dE_in + E_self_in

    END IF

    call cpu_time(time1)
    copy_time = copy_time + time1-time0

    IF (int_vdw_sum_style(box_in) == vdw_cut_tail) THEN
       nbeads_in(:) = nint_beads(:,box_in)

       DO i = 1, natoms(is)
          i_type = nonbond_list(i,is)%atom_type_number
          nint_beads(i_type,box_in) = nint_beads(i_type,box_in) + 1
       END DO

       CALL Compute_LR_Correction(box_in,E_lrc_in)
       dE_in = dE_in + E_lrc_in - energy(box_in)%lrc

    END IF

    IF(cpcollect) THEN

       potw = 1.0_DP
       CP_energy = dE_in

       potw = 1.0_DP / (P_forward * kappa_ins*kappa_rot*kappa_dih &
            ** (nfragments(is)-1))
       CP_energy = dE_in - E_angle_in

       chpot(is,box_in) = chpot(is,box_in) &
        + potw * (box_list(box_in)%volume &
        / (REAL(nmols(is,box_in)))) * DEXP(-beta(box_in) * CP_energy)

    END IF









     !P_forward = P_forward * REAL(nmols(is,box_out),DP) / REAL(nmols_box(box_out),DP)
     !P_reverse = P_reverse * REAL(nmols(is,box_in)+1,DP) / REAL(nmols_box(box_in)+1,DP)

  !ELSE
    !

  !END IF

  ! Increment counters
  ntrials(is,box_out)%deletion = ntrials(is,box_out)%deletion + 1
  ntrials(is,box)%insertion = ntrials(is,box)%insertion + 1
  ntrials(is,box)%cpcalc = ntrials(is,box)%cpcalc + 1


