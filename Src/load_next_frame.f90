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
        USE Internal_Coordinate_Routines
        !$ USE OMP_LIB

        IMPLICIT NONE
        
        INTEGER, DIMENSION(nbr_boxes) :: nspecies_thisframe
        INTEGER :: ibox, is, im
        LOGICAL :: end_reached

        end_reached = .FALSE.

        nmols = 0
        locate = 0
        molecule_list(:,:)%live = .FALSE.
        atom_list(:,:,:)%exist = .FALSE.
        molecule_list(:,:)%molecule_type = int_none
        molecule_list(:,:)%which_box = 0


        DO ibox = 1, nbr_boxes
                CALL Read_H_frame
                CALL Read_xyz_frame
        END DO

        DO is = 1, nspecies
                DO im = max_molecules(is), SUM(nmols(is,1:nbr_boxes)) + 1, -1
                        nmols(is,0) = nmols(is,0) + 1
                        locate(nmols(is,0),is,0) = im
                END DO
        END DO


   CONTAINS
        SUBROUTINE Read_H_frame

                INTEGER :: is_H
                INTEGER :: nmols_H
                INTEGER :: i, io
                INTEGER :: nvecsmax_old
                INTEGER :: AllocateStatus

                LOGICAL :: l_size_change


                REAL(DP) :: frame_volume
                REAL(DP), DIMENSION(3,3) :: this_length


                READ(pregen_H_unit(ibox),*,IOSTAT=io)
                IF (io < 0) THEN
                        end_reached = .TRUE.
                        RETURN
                END IF
                READ(pregen_H_unit(ibox),*)this_length(1,1), &
                        this_length(1,2), &
                        this_length(1,3)
                READ(pregen_H_unit(ibox),*)this_length(2,1), &
                        this_length(2,2), &
                        this_length(2,3)
                READ(pregen_H_unit(ibox),*)this_length(3,1), &
                        this_length(3,2), &
                        this_length(3,3)

                l_size_change = (.NOT. ALL(box_list(ibox)%length .EQ. this_length))

                IF (l_size_change) THEN
                        box_list(ibox)%length = this_length
                        CALL Compute_Cell_Dimensions(ibox)
                END IF


                READ(pregen_H_unit(ibox),*)
                READ(pregen_H_unit(ibox),*)nspecies_thisframe(ibox)
                DO i = 1,nspecies_thisframe(ibox) 
                        READ(pregen_H_unit(ibox),*)is_H, nmols_H
                        nmols_to_read(is_H,ibox) = nmols_H
                END DO

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

        END SUBROUTINE Read_H_frame

        SUBROUTINE Read_xyz_frame


                INTEGER :: is, ia, im, this_im, locate_base, this_unit, io

                CHARACTER(6) :: this_element

                REAL(DP) :: xcom_old, ycom_old, zcom_old
                REAL(DP) :: xcom_new, ycom_new, zcom_new
                REAL(DP) :: this_lambda, e_lrc
                LOGICAL :: overlap

                this_unit = pregen_xyz_unit(ibox)

                READ(this_unit,*,IOSTAT=io)
                IF (io < 0) THEN
                        end_reached = .TRUE.
                        RETURN
                END IF
                READ(this_unit,*)

                this_lambda = 1.0_DP
                ! Read in the coordinates of the molecules
                DO is = 1, nspecies
                        locate_base = SUM(nmols(is,1:nbr_boxes))

                        DO im = 1, nmols_to_read(is,ibox)
                                this_im = im + locate_base
                                locate(im,is,ibox) = this_im
                                DO ia = 1, natoms(is)
                                        READ(this_unit,*)this_element, &
                                                atom_list(ia,this_im,is)%rxp, &
                                                atom_list(ia,this_im,is)%ryp, &
                                                atom_list(ia,this_im,is)%rzp
                                        ! set the exist flag for this atom
                                        atom_list(ia,this_im,is)%exist = .TRUE.

                                END DO
                        END DO
                        !$OMP PARALLEL DO DEFAULT(SHARED) SCHEDULE(DYNAMIC) &
                        !$OMP PRIVATE(xcom_old, ycom_old, zcom_old, xcom_new, ycom_new, zcom_new)
                        DO this_im = (locate_base+1), (locate_base+nmols_to_read(is,ibox))
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

                                xcom_old = molecule_list(this_im,is)%xcom
                                ycom_old = molecule_list(this_im,is)%ycom
                                zcom_old = molecule_list(this_im,is)%zcom

                                ! Apply PBC
                                CALL Apply_PBC_Anint(ibox,xcom_old,ycom_old,zcom_old, &
                                                          xcom_new, ycom_new, zcom_new)

                                ! COM in the central simulation box
                                molecule_list(this_im,is)%xcom = xcom_new
                                molecule_list(this_im,is)%ycom = ycom_new
                                molecule_list(this_im,is)%zcom = zcom_new

                                ! displace atomic coordinates
                                atom_list(1:natoms(is),this_im,is)%rxp = &
                                        atom_list(1:natoms(is),this_im,is)%rxp + xcom_new - xcom_old
                                atom_list(1:natoms(is),this_im,is)%ryp = &
                                        atom_list(1:natoms(is),this_im,is)%ryp + ycom_new - ycom_old
                                atom_list(1:natoms(is),this_im,is)%rzp = &
                                        atom_list(1:natoms(is),this_im,is)%rzp + zcom_new - zcom_old
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


        END SUBROUTINE Read_xyz_frame


END SUBROUTINE Load_Next_Frame

