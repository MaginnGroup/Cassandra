!*****************************************************************************************
!
!
!*****************************************************************************************

SUBROUTINE Load_Next_Frame(end_reached)



        USE Global_Variables
        USE File_Names
        USE Simulation_Properties
        USE IO_Utilities

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

        IF (end_reached) RETURN

        
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


                REAL(DP) :: frame_volume


                READ(pregen_H_unit(ibox),*,IOSTAT=io)
                IF (io < 0) THEN
                        end_reached = .TRUE.
                        RETURN
                END IF
                READ(pregen_H_unit(ibox),*)box_list(ibox)%length(1,1), &
                        box_list(ibox)%length(1,2), &
                        box_list(ibox)%length(1,3)
                READ(pregen_H_unit(ibox),*)box_list(ibox)%length(2,1), &
                        box_list(ibox)%length(2,2), &
                        box_list(ibox)%length(2,3)
                READ(pregen_H_unit(ibox),*)box_list(ibox)%length(3,1), &
                        box_list(ibox)%length(3,2), &
                        box_list(ibox)%length(3,3)
                CALL Compute_Cell_Dimensions(ibox)

                READ(pregen_H_unit(ibox),*)
                READ(pregen_H_unit(ibox),*)nspecies_thisframe(ibox)
                DO i = 1,nspecies_thisframe(ibox) 
                        READ(pregen_H_unit(ibox),*)is_H, nmols_H
                        nmols_to_read(is_H,ibox) = nmols_H
                END DO

                IF (l_half_len_cutoff(ibox)) THEN
                        rcut_vdw(ibox) = 0.5 * MIN(box_list(ibox)%face_distance(1), &
                                                   box_list(ibox)%face_distance(2), &
                                                   box_list(ibox)%face_distance(3))
                        IF (int_charge_sum_style(ibox) /= charge_none) rcut_coul(ibox) = rcut_vdw(ibox)
                END IF

        END SUBROUTINE Read_H_frame

        SUBROUTINE Read_xyz_frame


                INTEGER :: is, ia, im, this_im, locate_base, this_unit, io

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

                ! Read in the coordinates of the molecules
                DO is = 1, nspecies
                        locate_base = SUM(nmols(is,1:nbr_boxes))

                        DO im = 1, nmols_to_read(is,ibox)
                                locate(im,is,ibox) = im + locate_base
                                this_im = locate(im,is,ibox)
                                molecule_list(this_im,is)%live = .TRUE.
                                this_lambda = 1.0_DP
                                ! By default all the molecules are normal
                                molecule_list(this_im,is)%molecule_type = int_normal

                                DO ia = 1, natoms(is)
                                        READ(this_unit,*)nonbond_list(ia,is)%element, &
                                                atom_list(ia,this_im,is)%rxp, &
                                                atom_list(ia,this_im,is)%ryp, &
                                                atom_list(ia,this_im,is)%rzp
                                        ! set the frac and exist flags for this atom
                                        molecule_list(this_im,is)%frac = this_lambda
                                        atom_list(ia,this_im,is)%exist = .TRUE.

                                END DO

                                ! assign the molecule the box id
                                molecule_list(this_im,is)%which_box = ibox
                                ! ensure that the molecular COM is inside the central simulation box
                                CALL Get_COM(this_im,is)

                                xcom_old = molecule_list(this_im,is)%xcom
                                ycom_old = molecule_list(this_im,is)%ycom
                                zcom_old = molecule_list(this_im,is)%zcom

                                ! Apply PBC
                                CALL Minimum_Image_Separation(ibox,xcom_old,ycom_old,zcom_old, &
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

                                nmols(is,ibox) = nmols(is,ibox) + 1

                        END DO
                END DO

                CALL Get_Internal_Coords

                ! Calculate COM and distance from the outermost atom to the COM

                DO is = 1, nspecies
                        DO im = 1, nmols(is,ibox)
                                this_im = locate(im,is,ibox)
                                IF( .NOT. molecule_list(this_im,is)%live) CYCLE
                                CALL Get_COM(this_im,is)
                                CALL Compute_Max_Com_Distance(this_im,is)
                        END DO
                END DO

                IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN 
                        CALL Compute_Beads(ibox)
                        CALL Compute_LR_Correction(ibox,e_lrc)
                        energy(this_box)%lrc = e_lrc
                END IF


        END SUBROUTINE Read_xyz_frame


END SUBROUTINE Load_Next_Frame

