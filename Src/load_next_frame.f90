!*****************************************************************************************
!
!
!*****************************************************************************************
SUBROUTINE Load_Next_Frame(end_reached)
        USE Global_Variables
        USE File_Names
        USE Simulation_Properties
        USE IO_Utilities
        USE Energy_Routines
        USE Type_Definitions
        USE XTC_Routines
        USE Internal_Coordinate_Routines
        !$ USE OMP_LIB
        IMPLICIT NONE

        INTEGER :: ibox, is, im
        LOGICAL :: end_reached
        REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE :: frame_xyz
        REAL(DP), DIMENSION(3,3), SAVE :: this_length

        end_reached = .FALSE.

        nmols = 0
        locate = 0
        molecule_list(:,:)%live = .FALSE.
        atom_list(:,:,:)%exist = .FALSE.
        molecule_list(:,:)%molecule_type = int_none
        molecule_list(:,:)%which_box = 0


        DO ibox = 1, nbr_boxes
                IF (has_Hfile(ibox)) THEN
                        this_length = Read_H_Frame()
                        IF (end_reached) RETURN
                ELSEIF (.NOT. ALLOCATED(frame_xyz)) THEN
                        ALLOCATE(frame_xyz(natoms_to_read(ibox),3))
                END IF
                IF (has_xyz(ibox)) THEN
                        frame_xyz = Read_xyz_Frame()
                        IF (end_reached) RETURN
                ELSEIF (has_xtc(ibox)) THEN
                        IF (Read_xtc_Frame(ibox)) THEN
                                end_reached = .TRUE.
                                EXIT
                        END IF
                        this_length = Get_xtc_Box(ibox)
                        frame_xyz = Get_xtc_Coords(ibox)
                END IF
                CALL Set_Frame_Box
                CALL Set_Frame_Coords
        END DO

        DO is = 1, nspecies
                DO im = max_molecules(is), SUM(nmols(is,1:nbr_boxes)) + 1, -1
                        nmols(is,0) = nmols(is,0) + 1
                        locate(nmols(is,0),is,0) = im
                END DO
        END DO


   CONTAINS
        SUBROUTINE Set_Frame_Box

                !REAL(DP), DIMENSION(3,3), INTENT(IN) :: this_length
                INTEGER :: nvecsmax_old
                INTEGER :: AllocateStatus

                LOGICAL :: l_size_change


                REAL(DP) :: frame_volume

                IF (end_reached) RETURN

                l_size_change = (.NOT. ALL(box_list(ibox)%length .EQ. this_length))

                IF (l_size_change) THEN
                        box_list(ibox)%length = this_length
                        CALL Compute_Cell_Dimensions(ibox)
                END IF

                IF (l_size_change .AND. l_half_len_cutoff(ibox)) THEN
                        rcut_vdw(ibox) = 0.5 * MIN(box_list(ibox)%face_distance(1), &
                                                   box_list(ibox)%face_distance(2), &
                                                   box_list(ibox)%face_distance(3))
                        rcut_vdwsq(ibox) = rcut_vdw(ibox) * rcut_vdw(ibox)
                        IF (int_charge_sum_style(ibox) /= charge_none) THEN
                                rcut_coul(ibox) = rcut_vdw(ibox)
                                rcut_coulsq(ibox) = rcut_vdwsq(ibox)
                        END IF

                        rcut_vdw3(ibox) = rcut_vdwsq(ibox) * rcut_vdw(ibox)
                        rcut_vdw6(ibox) = rcut_vdw3(ibox) * rcut_vdw3(ibox)
                        rcut3(ibox) = rcut_vdw3(ibox)
                        rcut9(ibox) = rcut3(ibox) * rcut_vdw6(ibox)

                        rcut_max(ibox) = rcut_vdw(ibox)
                        IF ( int_charge_sum_style(ibox) == charge_ewald) THEN
                                ! alpha_ewald(ibox) = ewald_p_sqrt(ibox) / rcut_coul(ibox)
                                h_ewald_cut(ibox) = 2.0_DP * ewald_p(ibox) / rcut_coul(ibox)
                        END IF
                END IF
                IF (l_size_change .AND. int_charge_sum_style(ibox) == charge_ewald) THEN
                        nvecsmax_old = MAXVAL(nvecs)
                        CALL Ewald_Reciprocal_Lattice_Vector_Setup(ibox)
                        IF (MAXVAL(nvecs) > nvecsmax_old) THEN
                                IF (ALLOCATED(cos_sum)) DEALLOCATE(cos_sum)
                                IF (ALLOCATED(sin_sum)) DEALLOCATE(sin_sum)
                                IF (ALLOCATED(cos_sum_old)) DEALLOCATE(cos_sum_old)
                                IF (ALLOCATED(sin_sum_old)) DEALLOCATE(sin_sum_old)
                                IF (ALLOCATED(cos_sum_start)) DEALLOCATE(cos_sum_start)
                                IF (ALLOCATED(sin_sum_start)) DEALLOCATE(sin_sum_start)
                                IF (ALLOCATED(cos_mol)) DEALLOCATE(cos_mol)
                                IF (ALLOCATED(sin_mol)) DEALLOCATE(sin_mol)
                                ALLOCATE(cos_sum(MAXVAL(nvecs),nbr_boxes), Stat = AllocateStatus)
                                IF (AllocateStatus /= 0) THEN
                                        err_msg = ''
                                        err_msg(1) = 'Memory could not be allocated for cos_sum'
                                        CALL Clean_Abort(err_msg,'Read_H_frame')
                                END IF
                                ALLOCATE(sin_sum(MAXVAL(nvecs),nbr_boxes), Stat = AllocateStatus)
                                IF (AllocateStatus /= 0) THEN
                                        err_msg = ''
                                        err_msg(1) = 'Memory could not be allocated for sin_sum'
                                        CALL Clean_Abort(err_msg,'Read_H_frame')
                                END IF
                                ALLOCATE(cos_mol(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)
                                IF (AllocateStatus /= 0) THEN
                                        err_msg = ''
                                        err_msg(1) = 'Memory could not be allocated for cos_mol'
                                        CALL Clean_Abort(err_msg,'Read_H_frame')
                                END IF
                                ALLOCATE(sin_mol(MAXVAL(nvecs),SUM(max_molecules)), Stat = AllocateStatus)
                                IF (AllocateStatus /= 0) THEN
                                        err_msg = ''
                                        err_msg(1) = 'Memory could not be allocated for sin_mol'
                                        CALL Clean_Abort(err_msg,'Read_H_frame')
                                END IF
                                ALLOCATE(cos_sum_old(SIZE(cos_sum,1),nbr_boxes), Stat = AllocateStatus)
                                IF (AllocateStatus /= 0) THEN
                                        err_msg = ''
                                        err_msg(1) = 'Memory could not be allocated for cos_mol_old'
                                        CALL Clean_Abort(err_msg,'Read_H_frame')
                                END IF
                                ALLOCATE(sin_sum_old(SIZE(sin_sum,1),nbr_boxes), Stat = AllocateStatus)
                                IF (AllocateStatus /= 0) THEN
                                        err_msg = ''
                                        err_msg(1) = 'Memory could not be allocated for sin_mol_old'
                                        CALL Clean_Abort(err_msg,'Read_H_frame')
                                END IF
                                ALLOCATE(cos_sum_start(SIZE(cos_sum,1),nbr_boxes), Stat = AllocateStatus)
                                IF (AllocateStatus /= 0) THEN
                                        err_msg = ''
                                        err_msg(1) = 'Memory could not be allocated for cos_mol_start'
                                        CALL Clean_Abort(err_msg,'Read_H_frame')
                                END IF
                                ALLOCATE(sin_sum_start(SIZE(sin_sum,1),nbr_boxes), Stat = AllocateStatus)
                                IF (AllocateStatus /= 0) THEN
                                        err_msg = ''
                                        err_msg(1) = 'Memory could not be allocated for sin_mol_start'
                                        CALL Clean_Abort(err_msg,'Read_H_frame')
                                END IF

                        END IF
                END IF

        END SUBROUTINE Set_Frame_Box

        FUNCTION Read_H_frame()
                REAL(DP), DIMENSION(3,3) :: Read_H_frame
                INTEGER :: nspecies_thisframe
                INTEGER :: is_H, is
                INTEGER :: nmols_H
                INTEGER :: i, io
                INTEGER :: old_natoms_to_read

                READ(pregen_H_unit(ibox),*,IOSTAT=io)
                IF (io < 0) THEN
                        end_reached = .TRUE.
                        RETURN
                END IF
                READ(pregen_H_unit(ibox),*)Read_H_frame(1,1), &
                        Read_H_frame(1,2), &
                        Read_H_frame(1,3)
                READ(pregen_H_unit(ibox),*)Read_H_frame(2,1), &
                        Read_H_frame(2,2), &
                        Read_H_frame(2,3)
                READ(pregen_H_unit(ibox),*)Read_H_frame(3,1), &
                        Read_H_frame(3,2), &
                        Read_H_frame(3,3)

                READ(pregen_H_unit(ibox),*)
                READ(pregen_H_unit(ibox),*)nspecies_thisframe
                nmols_to_read(:,ibox) = 0
                DO i = 1,nspecies_thisframe 
                        READ(pregen_H_unit(ibox),*)is_H, nmols_H
                        nmols_to_read(is_H,ibox) = nmols_H
                END DO
                atom_ibounds(2,:,ibox) = natoms*nmols_to_read(:,ibox)
                old_natoms_to_read = natoms_to_read(ibox)
                natoms_to_read(ibox) = SUM(atom_ibounds(2,:,ibox))
                DO is = 2, nspecies
                        atom_ibounds(2,is,ibox) = SUM(atom_ibounds(2,is-1:is,ibox))
                END DO
                atom_ibounds(1,1,ibox) = 1
                IF (nspecies > 1) atom_ibounds(1,2:nspecies,ibox) = atom_ibounds(2,1:(nspecies-1),ibox)+1
                IF (natoms_to_read(ibox) .NE. old_natoms_to_read) THEN
                        IF (ALLOCATED(frame_xyz)) DEALLOCATE(frame_xyz)
                        ALLOCATE(frame_xyz(natoms_to_read(ibox),3))
                END IF

        END FUNCTION Read_H_frame

        SUBROUTINE Set_Frame_Coords

                !REAL(DP), DIMENSION(natoms_to_read(ibox),3), INTENT(IN) :: frame_xyz
                INTEGER :: is, ia, imol, this_im, locate_base


                REAL(DP) :: xcom_old, ycom_old, zcom_old
                REAL(DP) :: xcom_new, ycom_new, zcom_new
                REAL(DP) :: this_lambda, e_lrc
                LOGICAL :: overlap

                TYPE(Atom_Class), POINTER :: al_ptr(:,:)
                INTEGER :: newshape(2), sloc, eloc, aib(2)

                IF (end_reached) RETURN

                this_lambda = 1.0_DP
                ! Read in the coordinates of the molecules
                DO is = 1, nspecies
                        IF (nmols_to_read(is,ibox) < 1) CYCLE
                        locate_base = SUM(nmols(is,1:nbr_boxes))
                        DO imol = 1, nmols_to_read(is,ibox)
                                locate(imol,is,ibox) = imol+locate_base
                        END DO
                        sloc = locate_base + 1
                        eloc = locate_base +nmols_to_read(is,ibox)
                        aib = atom_ibounds(:,is,ibox)
                        al_ptr => atom_list(1:natoms(is),sloc:eloc,is)
                        newshape(1) = natoms(is)
                        newshape(2) = nmols_to_read(is,ibox)
                        al_ptr%rp(1) = &
                                RESHAPE(frame_xyz(aib(1):aib(2),1), newshape) 
                        al_ptr%rp(2) = &
                                RESHAPE(frame_xyz(aib(1):aib(2),2), newshape) 
                        al_ptr%rp(3) = &
                                RESHAPE(frame_xyz(aib(1):aib(2),3), newshape) 
                        al_ptr%exist = .TRUE.
                        !$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) &
                        !$OMP PRIVATE(xcom_old, ycom_old, zcom_old, xcom_new, ycom_new, zcom_new)
                        DO this_im = sloc, eloc
                                molecule_list(this_im,is)%live = .TRUE.
                                ! By default all the molecules are normal
                                molecule_list(this_im,is)%molecule_type = int_normal
                                molecule_list(this_im,is)%frac = this_lambda
                                ! assign the molecule the box id
                                molecule_list(this_im,is)%which_box = ibox
                                ! ensure that the molecular COM is inside the central simulation box
                                ! Calculate COM and distance from the outermost atom to the COM
                                CALL Get_COM(this_im,is)
                                CALL Compute_Max_Com_Distance(this_im,is)

                                xcom_old = molecule_list(this_im,is)%rcom(1)
                                ycom_old = molecule_list(this_im,is)%rcom(2)
                                zcom_old = molecule_list(this_im,is)%rcom(3)

                                ! Apply PBC
                                CALL Apply_PBC_Anint(ibox,xcom_old,ycom_old,zcom_old, &
                                                          xcom_new, ycom_new, zcom_new)

                                ! COM in the central simulation box
                                molecule_list(this_im,is)%rcom(1) = xcom_new
                                molecule_list(this_im,is)%rcom(2) = ycom_new
                                molecule_list(this_im,is)%rcom(3) = zcom_new

                                ! displace atomic coordinates
                                atom_list(1:natoms(is),this_im,is)%rp(1) = &
                                        atom_list(1:natoms(is),this_im,is)%rp(1) + xcom_new - xcom_old
                                atom_list(1:natoms(is),this_im,is)%rp(2) = &
                                        atom_list(1:natoms(is),this_im,is)%rp(2) + ycom_new - ycom_old
                                atom_list(1:natoms(is),this_im,is)%rp(3) = &
                                        atom_list(1:natoms(is),this_im,is)%rp(3) + zcom_new - zcom_old
                        END DO
                        !$OMP END PARALLEL DO
                        nmols(is,ibox) = nmols(is,ibox) + nmols_to_read(is,ibox)
                END DO

                CALL Get_Internal_Coords

                IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN 
                        CALL Compute_Beads(ibox)
                        CALL Compute_LR_Correction(ibox,e_lrc)
                        energy(ibox)%lrc = e_lrc
                END IF
                
                IF (int_charge_sum_style(ibox) == charge_ewald) THEN
                        CALL Compute_System_Ewald_Reciprocal_Energy(ibox)
                END IF


        END SUBROUTINE Set_Frame_Coords

        FUNCTION Read_xyz_frame()
                REAL(DP), DIMENSION(natoms_to_read(ibox),3) :: Read_xyz_frame
                INTEGER :: this_unit, io, i
                CHARACTER(6) :: this_element

                this_unit = pregen_xyz_unit(ibox)

                READ(this_unit,*,IOSTAT=io)
                IF (io < 0) THEN
                        end_reached = .TRUE.
                        RETURN
                END IF
                READ(this_unit,*)
                DO i = 1, natoms_to_read(ibox)
                        READ(this_unit,*) this_element, &
                                Read_xyz_frame(i,1), &
                                Read_xyz_frame(i,2), &
                                Read_xyz_frame(i,3)
                END DO
                
        END FUNCTION Read_xyz_frame


END SUBROUTINE Load_Next_Frame

