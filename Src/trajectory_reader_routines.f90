!*****************************************************************************************
!
!
!*****************************************************************************************
MODULE Trajectory_Reader_Routines
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
        REAL(DP), DIMENSION(:,:), ALLOCATABLE, SAVE, PRIVATE :: frame_xyz
        REAL(DP), DIMENSION(3,3), SAVE, PRIVATE :: this_length
        INTEGER, PRIVATE :: ibox, natr_p

        CONTAINS
        SUBROUTINE Load_Next_Frame

                INTEGER :: is, im

                early_end = .FALSE.

                nmols = 0
                locate = 0
                molecule_list(:,:)%live = .FALSE.
                atom_list(:,:,:)%exist = .FALSE.
                molecule_list(:,:)%molecule_type = int_none
                molecule_list(:,:)%which_box = 0


                DO ibox = 1, nbr_boxes
                        IF (has_Hfile(ibox)) THEN
                                this_length = Read_H_Frame()
                                IF (early_end) RETURN
                        ELSEIF (.NOT. ALLOCATED(frame_xyz)) THEN
                                natr_p = IAND(MAXVAL(natoms_to_read)+7,NOT(7))
                                ALLOCATE(frame_xyz(natr_p,3))
                                !DIR$ VECTOR ALIGNED
                                frame_xyz = 0.0_DP
                        END IF
                        IF (has_xyz(ibox)) THEN
                                frame_xyz(1:natoms_to_read(ibox),:) = Read_xyz_Frame()
                                IF (early_end) RETURN
                        ELSEIF (has_xtc(ibox)) THEN
                                IF (Read_xtc_Frame(ibox)) THEN
                                        early_end = .TRUE.
                                        EXIT
                                END IF
                                this_length = Get_xtc_Box(ibox)
                                CALL Get_xtc_Coords(ibox,frame_xyz)
                                !frame_xyz(1:natoms_to_read(ibox),:) = Get_xtc_Coords(ibox)
                        END IF
                        CALL Set_Frame_Box
                        IF (.NOT. early_end) CALL Set_Frame_Coords
                END DO

                DO is = 1, nspecies
                        DO im = max_molecules(is), SUM(nmols(is,1:nbr_boxes)) + 1, -1
                                nmols(is,0) = nmols(is,0) + 1
                                locate(nmols(is,0),is,0) = im
                        END DO
                END DO
        END SUBROUTINE Load_Next_Frame


        SUBROUTINE Set_Frame_Box

                !REAL(DP), DIMENSION(3,3), INTENT(IN) :: this_length
                INTEGER :: nvecsmax_old
                INTEGER :: AllocateStatus

                LOGICAL :: l_size_change


                REAL(DP) :: frame_volume

                IF (early_end) RETURN

                l_size_change = (.NOT. ALL(box_list(ibox)%orig_length .EQ. this_length))

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
                        CALL Ewald_Reciprocal_Lattice_Vector_Setup(ibox)
                END IF

        END SUBROUTINE Set_Frame_Box

        FUNCTION Read_H_frame()
                REAL(DP), DIMENSION(3,3) :: Read_H_frame
                INTEGER :: nspecies_thisframe
                INTEGER :: is_H, is
                INTEGER :: nmols_H
                INTEGER :: i, io

                READ(pregen_H_unit(ibox),*,IOSTAT=io)
                IF (io < 0) THEN
                        early_end = .TRUE.
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
                natoms_to_read(ibox) = SUM(atom_ibounds(2,:,ibox))
                DO is = 2, nspecies
                        atom_ibounds(2,is,ibox) = SUM(atom_ibounds(2,is-1:is,ibox))
                END DO
                atom_ibounds(1,1,ibox) = 1
                IF (nspecies > 1) atom_ibounds(1,2:nspecies,ibox) = atom_ibounds(2,1:(nspecies-1),ibox)+1
                natr_p = IAND(natoms_to_read(ibox)+7,NOT(7))
                IF (ALLOCATED(frame_xyz)) THEN
                        IF (SIZE(frame_xyz,1)<natr_p) DEALLOCATE(frame_xyz)
                END IF
                IF (.NOT. ALLOCATED(frame_xyz)) THEN
                        ALLOCATE(frame_xyz(natr_p,3))
                        frame_xyz = 0.0_DP
                END IF

        END FUNCTION Read_H_frame

        SUBROUTINE Set_Frame_Coords

                !REAL(DP), DIMENSION(natoms_to_read(ibox),3), INTENT(IN) :: frame_xyz
                INTEGER :: is, ia, ja, imol, this_im, locate_base, im


                REAL(DP) :: this_lambda, e_lrc
                LOGICAL :: overlap

                INTEGER :: sloc, eloc, aib(2)
                !REAL(DP), DIMENSION(IAND(MAXVAL(nmols_to_read(:,ibox))+3,NOT(3)),MAXVAL(natoms),3) :: &
                !        frame_xyz_rs
                !REAL(DP), DIMENSION(IAND(MAXVAL(nmols_to_read(:,ibox))+3,NOT(3)),4) :: &
                !        molwrapvec, rcom_arr
                REAL(DP) :: h11,h12,h13,h21,h22,h23,h31,h32,h33
                REAL(DP) :: inv_h11,inv_h12,inv_h13,inv_h21,inv_h22,inv_h23,inv_h31,inv_h32,inv_h33
                INTEGER :: ntr, ibond, i_dim
                !LOGICAL, DIMENSION(MAXVAL(natoms)) :: l_moved
                REAL(DP) :: rxp,ryp,rzp,sxp,syp,szp,rcom,scom,rxcom,rycom,rzcom,sxcom,sycom,szcom
                REAL(DP) :: dxcom, dycom, dzcom
                REAL(DP) :: boxlen,hl,irp,jrp, drp, absdrp, inv_total_mass, massfrac
                !REAL(DP), DIMENSION(MAXVAL(natoms)) :: massfrac_vec(MAXVAL(natoms))
                REAL(DP) :: isp, xcom, ycom, zcom, dcomsq, max_dcomsq, inv_l
                REAL(DP) :: sxmwv, symwv, szmwv, mwv
                LOGICAL :: l_ortho
                INTEGER :: chunksize, chunkstart, chunkend, ithread, nthreads, chunkshift
                INTEGER :: kx, ky, kz, i

                l_ortho = box_list(ibox)%int_box_shape <= int_ortho
                !$OMP PARALLEL DEFAULT(SHARED) &
                !$OMP PRIVATE(h11,h12,h13,h21,h22,h23,h31,h32,h33) &
                !$OMP PRIVATE(inv_h11,inv_h12,inv_h13,inv_h21,inv_h22,inv_h23,inv_h31,inv_h32,inv_h33) &
                !$OMP PRIVATE(locate_base,is, ia,ja, imol,im, this_im) &
                !$OMP PRIVATE(rxp,ryp,rzp,sxp,syp,szp,rcom,scom,rxcom,rycom,rzcom,sxcom,sycom,szcom) &
                !$OMP PRIVATE(this_lambda,ntr,ibond,boxlen,hl,irp,jrp,drp,absdrp,inv_total_mass) &
                !$OMP PRIVATE(massfrac,isp,xcom, ycom, zcom, dcomsq, max_dcomsq, inv_l) &
                !$OMP PRIVATE(sxmwv, symwv, szmwv, mwv,dxcom,dycom,dzcom,i_dim) &
                !$OMP PRIVATE(chunksize,chunkstart,chunkend,ithread,nthreads,chunkshift)

                !!$OMP PARALLEL DEFAULT(PRIVATE) &
                !!$OMP SHARED(frame_xyz_rs,ibox,frame_xyz,rcom_arr,molwrapvec) & 
                !!$OMP SHARED(box_list,atom_list,molecule_list,natoms,nspecies,nmols_to_read,nonbond_list) &
                !!$OMP SHARED(bondpart_list,nmols,locate,atom_ibounds,nbr_boxes,natr_p)



                IF (.NOT. l_ortho) THEN
                        h11 = box_list(ibox)%length(1,1)
                        h21 = box_list(ibox)%length(2,1)
                        h31 = box_list(ibox)%length(3,1)
                        h12 = box_list(ibox)%length(1,2)
                        h22 = box_list(ibox)%length(2,2)
                        h32 = box_list(ibox)%length(3,2)
                        h13 = box_list(ibox)%length(1,3)
                        h23 = box_list(ibox)%length(2,3)
                        h33 = box_list(ibox)%length(3,3)
                        IF (box_list(ibox)%basis_changed) THEN
                                inv_h11 = box_list(ibox)%orig_length_inv(1,1)
                                inv_h21 = box_list(ibox)%orig_length_inv(2,1)
                                inv_h31 = box_list(ibox)%orig_length_inv(3,1)
                                inv_h12 = box_list(ibox)%orig_length_inv(1,2)
                                inv_h22 = box_list(ibox)%orig_length_inv(2,2)
                                inv_h32 = box_list(ibox)%orig_length_inv(3,2)
                                inv_h13 = box_list(ibox)%orig_length_inv(1,3)
                                inv_h23 = box_list(ibox)%orig_length_inv(2,3)
                                inv_h33 = box_list(ibox)%orig_length_inv(3,3)
                                !DIR$ ASSUME (MOD(natr_p,4) .EQ. 0)
                                !DIR$ VECTOR ALIGNED
                                !$OMP DO SIMD SCHEDULE(SIMD:STATIC) PRIVATE(rxp,ryp,rzp,sxp,syp,szp)
                                DO ia = 1, natr_p
                                        rxp = frame_xyz(ia,1)
                                        ryp = frame_xyz(ia,2)
                                        rzp = frame_xyz(ia,3)
                                        sxp = rxp*inv_h11 + ryp*inv_h12 + rzp*inv_h13
                                        syp = rxp*inv_h21 + ryp*inv_h22 + rzp*inv_h23
                                        szp = rxp*inv_h31 + ryp*inv_h32 + rzp*inv_h33
                                        frame_xyz(ia,1) = sxp
                                        frame_xyz(ia,2) = syp
                                        frame_xyz(ia,3) = szp
                                END DO
                                !$OMP END DO SIMD
                                inv_h11 = box_list(ibox)%length_inv(1,1)
                                inv_h21 = box_list(ibox)%length_inv(2,1)
                                inv_h31 = box_list(ibox)%length_inv(3,1)
                                inv_h12 = box_list(ibox)%length_inv(1,2)
                                inv_h22 = box_list(ibox)%length_inv(2,2)
                                inv_h32 = box_list(ibox)%length_inv(3,2)
                                inv_h13 = box_list(ibox)%length_inv(1,3)
                                inv_h23 = box_list(ibox)%length_inv(2,3)
                                inv_h33 = box_list(ibox)%length_inv(3,3)
                        ELSE
                                inv_h11 = box_list(ibox)%length_inv(1,1)
                                inv_h21 = box_list(ibox)%length_inv(2,1)
                                inv_h31 = box_list(ibox)%length_inv(3,1)
                                inv_h12 = box_list(ibox)%length_inv(1,2)
                                inv_h22 = box_list(ibox)%length_inv(2,2)
                                inv_h32 = box_list(ibox)%length_inv(3,2)
                                inv_h13 = box_list(ibox)%length_inv(1,3)
                                inv_h23 = box_list(ibox)%length_inv(2,3)
                                inv_h33 = box_list(ibox)%length_inv(3,3)
                                !DIR$ ASSUME (MOD(natr_p,4) .EQ. 0)
                                !DIR$ VECTOR ALIGNED
                                !$OMP DO SIMD SCHEDULE(SIMD:STATIC) PRIVATE(rxp,ryp,rzp,sxp,syp,szp)
                                DO ia = 1, natr_p
                                        rxp = frame_xyz(ia,1)
                                        ryp = frame_xyz(ia,2)
                                        rzp = frame_xyz(ia,3)
                                        sxp = rxp*inv_h11 + ryp*inv_h12 + rzp*inv_h13
                                        syp =               ryp*inv_h22 + rzp*inv_h23
                                        szp =                             rzp*inv_h33
                                        frame_xyz(ia,1) = sxp
                                        frame_xyz(ia,2) = syp
                                        frame_xyz(ia,3) = szp
                                END DO
                                !$OMP END DO SIMD
                        END IF
                END IF

                !this_lambda = 1.0_DP
                nthreads = 1
                ithread = 0
                !$ nthreads = OMP_GET_NUM_THREADS()
                !$ ithread = OMP_GET_THREAD_NUM()
                DO is = 1, nspecies
                        IF (nmols_to_read(is,ibox) < 1) CYCLE
                        ntr = IAND(nmols_to_read(is,ibox)+3,NOT(3))
                        chunksize = ntr
                        !$ chunksize = IAND((ntr+nthreads-1)/nthreads+3,NOT(3))
                        !chunkstart = 1
                        chunkshift = 0
                        chunkend = ntr
                        !$ chunkshift = ithread*chunksize
                        !!$ chunkstart = chunkshift+1
                        !$ chunkend = MIN((ithread+1)*chunksize,ntr)
                        !$ chunksize = MAX(0,chunkend-chunkshift)
                        IF (chunksize .EQ. 0) CYCLE
                        CALL Set_Species_Frame_Coords(is,natoms(is),chunksize,chunkshift,l_ortho)
                        !locate_base = SUM(nmols(is,1:nbr_boxes))
                        !!$OMP SINGLE
                        !DO imol = 1, nmols_to_read(is,ibox)
                        !        locate(imol,is,ibox) = imol+locate_base
                        !END DO
                        !sloc = locate_base + 1
                        !eloc = locate_base +nmols_to_read(is,ibox)
                        !aib = atom_ibounds(:,is,ibox)
                        !!$OMP END SINGLE
                        !!$OMP WORKSHARE
                        !frame_xyz_rs(1:nmols_to_read(is,ibox),1:natoms(is),1:3) = &
                        !        RESHAPE(frame_xyz(aib(1):aib(2),:),&
                        !        (/ nmols_to_read(is,ibox),natoms(is),3 /),ORDER=(/2,1,3/))
                        !!$OMP END WORKSHARE
                        !l_moved = .FALSE.
                        !l_moved(1) = .TRUE.
                        !!DIR$ ASSUME (MOD(chunkend,4) .EQ. 0)
                        !!DIR$ ASSUME (MOD(chunkstart,4) .EQ. 1)
                        !DO WHILE (.NOT. ALL(l_moved(1:natoms(is))))
                        !!WRITE(*,*) l_moved(1:natoms(is))
                        !DO ia = 1, natoms(is)
                        !        IF (.NOT. l_moved(ia)) CYCLE
                        !        DO ibond = 1 , bondpart_list(ia,is)%nbonds
                        !                ja = bondpart_list(ia,is)%atom(ibond)
                        !                IF (l_moved(ja)) CYCLE
                        !                l_moved(ja) = .TRUE.
                        !                DO i_dim = 1, 3
                        !                        IF (l_ortho) THEN
                        !                                boxlen = box_list(ibox)%length(i_dim,i_dim)
                        !                        ELSE
                        !                                boxlen = 1.0_DP
                        !                        END IF
                        !                        hl = 0.5_DP * boxlen
                        !                        !DIR$ VECTOR ALIGNED
                        !                        !$OMP SIMD PRIVATE(irp,jrp,drp,absdrp)
                        !                        DO im = chunkstart, chunkend
                        !                                irp = frame_xyz_rs(im,ia,i_dim)
                        !                                jrp = frame_xyz_rs(im,ja,i_dim)
                        !                                drp = jrp - irp
                        !                                absdrp = ABS(drp)
                        !                                IF (absdrp > hl) jrp = jrp - SIGN(boxlen,drp)
                        !                                frame_xyz_rs(im,ja,i_dim) = jrp
                        !                        END DO
                        !                        !$OMP END SIMD
                        !                END DO
                        !        END DO
                        !END DO
                        !END DO
                        !IF (.NOT. l_ortho) THEN
                        !        DO ia = 1, natoms(is)
                        !                !DIR$ VECTOR ALIGNED
                        !                !$OMP SIMD PRIVATE(isp,rxp,ryp,rzp)
                        !                DO im = chunkstart, chunkend
                        !                        isp = frame_xyz_rs(im,ia,1)
                        !                        rxp = h11*isp
                        !                        isp = frame_xyz_rs(im,ia,2)
                        !                        rxp = rxp + h12*isp
                        !                        ryp = h22*isp
                        !                        isp = frame_xyz_rs(im,ia,3)
                        !                        rxp = rxp + h13*isp
                        !                        ryp = ryp + h23*isp
                        !                        rzp = h33*isp
                        !                        frame_xyz_rs(im,ia,1) = rxp
                        !                        frame_xyz_rs(im,ia,2) = ryp
                        !                        frame_xyz_rs(im,ia,3) = rzp
                        !                END DO
                        !                !$OMP END SIMD
                        !        END DO
                        !END IF
                        !inv_total_mass = 1.0_DP/SUM(nonbond_list(1:natoms(is),is)%mass)
                        !massfrac_vec(1:natoms(is)) = nonbond_list(1:natoms(is),is)%mass*inv_total_mass
                        !!DIR$ VECTOR ALIGNED
                        !!$OMP SIMD &
                        !!$OMP PRIVATE(rxp,ryp,rzp,xcom,ycom,zcom,dxcom,dycom,dzcom) &
                        !!$OMP PRIVATE(max_dcomsq,dcomsq,massfrac)
                        !DO im = chunkstart, chunkend
                        !        massfrac = massfrac_vec(1)
                        !        rxp = frame_xyz_rs(im,1,1)
                        !        ryp = frame_xyz_rs(im,1,2)
                        !        rzp = frame_xyz_rs(im,1,3)
                        !        xcom = massfrac*rxp
                        !        ycom = massfrac*ryp
                        !        zcom = massfrac*rzp
                        !        DO ia = 2, natoms(is)
                        !                massfrac = massfrac_vec(ia)
                        !                rxp = frame_xyz_rs(im,ia,1)
                        !                ryp = frame_xyz_rs(im,ia,2)
                        !                rzp = frame_xyz_rs(im,ia,3)
                        !                xcom = xcom + massfrac*rxp
                        !                ycom = ycom + massfrac*ryp
                        !                zcom = zcom + massfrac*rzp
                        !        END DO
                        !        rcom_arr(im,1) = xcom
                        !        rcom_arr(im,2) = ycom
                        !        rcom_arr(im,3) = zcom
                        !        rxp = frame_xyz_rs(im,1,1)
                        !        ryp = frame_xyz_rs(im,1,2)
                        !        rzp = frame_xyz_rs(im,1,3)
                        !        dxcom = rxp - xcom
                        !        dycom = ryp - ycom
                        !        dzcom = rzp - zcom
                        !        max_dcomsq = dxcom*dxcom + dycom*dycom + dzcom*dzcom
                        !        DO ia = 2, natoms(is)
                        !                rxp = frame_xyz_rs(im,ia,1)
                        !                ryp = frame_xyz_rs(im,ia,2)
                        !                rzp = frame_xyz_rs(im,ia,3)
                        !                dxcom = rxp - xcom
                        !                dycom = ryp - ycom
                        !                dzcom = rzp - zcom
                        !                dcomsq = dxcom*dxcom + dycom*dycom + dzcom*dzcom
                        !                max_dcomsq = MAX(max_dcomsq,dcomsq)
                        !        END DO
                        !        rcom_arr(im,4) = SQRT(max_dcomsq)
                        !END DO
                        !!$OMP END SIMD
                        !IF (l_ortho) THEN
                        !        DO i_dim = 1, 3
                        !                boxlen = box_list(ibox)%length(i_dim,i_dim)
                        !                inv_l = 1.0_DP/boxlen
                        !                !DIR$ VECTOR ALIGNED
                        !                !$OMP SIMD PRIVATE(rcom,mwv,irp)
                        !                DO im = chunkstart, chunkend
                        !                        rcom = rcom_arr(im,i_dim)
                        !                        mwv = boxlen*ANINT(rcom*inv_l)
                        !                        !molwrapvec(im,i_dim) = mwv
                        !                        rcom_arr(im,i_dim) = rcom - mwv
                        !                        DO ia = 1, natoms(is)
                        !                                irp = frame_xyz_rs(im,ia,i_dim)
                        !                                irp = irp - mwv
                        !                                frame_xyz_rs(im,ia,i_dim) = irp
                        !                        END DO
                        !                END DO
                        !                !$OMP END SIMD
                        !        END DO
                        !ELSE
                        !        !DIR$ VECTOR ALIGNED
                        !        !$OMP SIMD PRIVATE(rxcom,rycom,rzcom)
                        !        DO im = chunkstart, chunkend
                        !                rxcom = rcom_arr(im,1)
                        !                rycom = rcom_arr(im,2)
                        !                rzcom = rcom_arr(im,3)
                        !                molwrapvec(im,1) = rxcom*inv_h11 + rycom*inv_h12 + rzcom*inv_h13
                        !                molwrapvec(im,2) =                 rycom*inv_h22 + rzcom*inv_h23
                        !                molwrapvec(im,3) =                                 rzcom*inv_h33
                        !        END DO
                        !        !$OMP END SIMD
                        !        !DIR$ VECTOR ALIGNED
                        !        !$OMP SIMD PRIVATE(mwv)
                        !        DO im = chunkstart, chunkend
                        !                mwv = molwrapvec(im,1)
                        !                molwrapvec(im,1) = ANINT(mwv)
                        !                mwv = molwrapvec(im,2)
                        !                molwrapvec(im,2) = ANINT(mwv)
                        !                mwv = molwrapvec(im,3)
                        !                molwrapvec(im,3) = ANINT(mwv)
                        !        END DO
                        !        !$OMP END SIMD
                        !        !molwrapvec(1:nmols_to_read(is,ibox),:) = ANINT(molwrapvec(1:nmols_to_read(is,ibox),:))
                        !        !DIR$ VECTOR ALIGNED
                        !        !$OMP SIMD PRIVATE(sxmwv,symwv,szmwv)
                        !        DO im = chunkstart, chunkend
                        !                sxmwv = molwrapvec(im,1)
                        !                symwv = molwrapvec(im,2)
                        !                szmwv = molwrapvec(im,3)
                        !                molwrapvec(im,1) = sxmwv*h11 + symwv*h12 + szmwv*h13
                        !                molwrapvec(im,2) =             symwv*h22 + szmwv*h23
                        !                molwrapvec(im,3) =                         szmwv*h33
                        !        END DO
                        !        !$OMP END SIMD
                        !        DO i_dim = 1, 3
                        !                !DIR$ VECTOR ALIGNED
                        !                !$OMP SIMD PRIVATE(rcom,mwv,irp)
                        !                DO im = chunkstart, chunkend
                        !                        rcom = rcom_arr(im,i_dim)
                        !                        mwv = molwrapvec(im,i_dim)
                        !                        rcom = rcom - mwv
                        !                        rcom_arr(im,i_dim) = rcom
                        !                        DO ia = 1, natoms(is)
                        !                                irp = frame_xyz_rs(im,ia,i_dim)
                        !                                irp = irp - mwv
                        !                                frame_xyz_rs(im,ia,i_dim) = irp
                        !                        END DO
                        !                END DO
                        !                !$OMP END SIMD
                        !        END DO
                        !END IF
                        !!DO ia = 1, natoms(is)
                        !!        DO i_dim = 1, 3
                        !!                !DIR$ VECTOR ALIGNED
                        !!                DO im = chunkstart, chunkend
                        !!                        mwv = molwrapvec(im,i_dim)
                        !!                        irp = frame_xyz_rs(im,ia,i_dim)
                        !!                        irp = irp - mwv
                        !!                        frame_xyz_rs(im,ia,i_dim) = irp
                        !!                END DO
                        !!        END DO
                        !!END DO
                        !!$OMP BARRIER
                        !!$OMP DO SCHEDULE(STATIC)
                        !DO im = sloc, eloc
                        !        molecule_list(im,is)%live = .TRUE.
                        !        molecule_list(im,is)%frac = this_lambda
                        !        molecule_list(im,is)%which_box = ibox
                        !        molecule_list(im,is)%rcom = rcom_arr(im-locate_base,:)
                        !        !atom_list(1:natoms(is),im,is)%exist = .TRUE.
                        !END DO
                        !!$OMP END DO NOWAIT
                        !!$OMP DO SCHEDULE(STATIC) COLLAPSE(2)
                        !DO ia = 1, natoms(is)
                        !DO im = sloc, eloc
                        !        atom_list(ia,im,is)%exist = .TRUE.
                        !        atom_list(ia,im,is)%rp(1:3) = frame_xyz_rs(im-locate_base,ia,1:3)
                        !END DO
                        !END DO
                        !!$OMP END DO NOWAIT
                        !!$OMP SINGLE
                        !nmols(is,ibox) = nmols(is,ibox) + nmols_to_read(is,ibox)
                        !!$OMP END SINGLE

                END DO
                IF (ithread+1 .EQ. nthreads) THEN
                        DO is = 1, nspecies
                                locate_base = SUM(nmols(is,1:nbr_boxes))
                                DO imol = 1, nmols_to_read(is,ibox)
                                        locate(imol,is,ibox) = imol+locate_base
                                END DO
                        END DO
                END IF
                !$OMP BARRIER
                !$OMP END PARALLEL
                nmols(:,ibox) = nmols(:,ibox) + nmols_to_read(:,ibox)
                !WRITE(*,*) sloc, eloc, is, ibox, this_lambda
                !WRITE(*,*) ALL(molecule_list(sloc:eloc,is)%live), COUNT(molecule_list(sloc:eloc,is)%live), &
                !        SUM(molecule_list(sloc:eloc,is)%which_box)

                !WRITE(*,*) natoms
                !WRITE(*,*) natoms_to_read
                !WRITE(*,*) nmols_to_read
                !WRITE(*,*) nmols
                !WRITE(*,*) nspecies
                !WRITE(*,*) nbr_boxes

                !CALL Get_Internal_Coords

                IF (int_vdw_sum_style(ibox) == vdw_cut_tail) THEN 
                        CALL Compute_Beads(ibox)
                        CALL Compute_LR_Correction(ibox,e_lrc)
                        energy(ibox)%lrc = e_lrc
                ELSE
                        energy(ibox)%lrc = 0.0_DP
                END IF
                
                IF (int_charge_sum_style(ibox) == charge_ewald) THEN
                        CALL Compute_System_Ewald_Reciprocal_Energy_Pregen(ibox)
                ELSE
                        energy(ibox)%reciprocal = 0.0_DP
                END IF


        END SUBROUTINE Set_Frame_Coords
        SUBROUTINE Set_Species_Frame_Coords(is,na,chunksize,chunkshift,l_ortho)
                INTEGER, INTENT(IN) :: is,na,chunksize,chunkshift
                LOGICAL, INTENT(IN) :: l_ortho


                REAL(DP), DIMENSION(chunksize,na,3) :: frame_xyz_rs
                REAL(DP), DIMENSION(chunksize,4) :: rcom_arr
                REAL(DP), DIMENSION(MERGE(chunksize,0,l_ortho),3) :: molwrapvec
                LOGICAL, DIMENSION(na) :: l_moved
                INTEGER :: aib(2), ia, ja, im, locate_base
                REAL(DP) :: h11,h12,h13,h21,h22,h23,h31,h32,h33
                REAL(DP) :: inv_h11,inv_h12,inv_h13,inv_h21,inv_h22,inv_h23,inv_h31,inv_h32,inv_h33
                INTEGER :: ntr, ibond, i_dim
                REAL(DP) :: rxp,ryp,rzp,sxp,syp,szp,rcom,scom,rxcom,rycom,rzcom,sxcom,sycom,szcom
                REAL(DP) :: dxcom, dycom, dzcom
                REAL(DP) :: boxlen,hl,irp,jrp, drp, absdrp, inv_total_mass, massfrac
                REAL(DP), DIMENSION(na) :: massfrac_vec
                REAL(DP) :: isp, xcom, ycom, zcom, dcomsq, max_dcomsq, inv_l
                REAL(DP) :: sxmwv, symwv, szmwv, mwv
                INTEGER :: loc, chunkend
                !l_ortho = box_list(ibox)%int_box_shape <= int_ortho
                aib(1) = atom_ibounds(1,is,ibox) + chunkshift*na
                chunkend = MIN(chunksize,nmols_to_read(is,ibox)-chunkshift) ! unpadded chunksize
                aib(2) = aib(1) + chunkend*na - 1
                frame_xyz_rs(1:chunkend,:,:) = &
                        RESHAPE(frame_xyz(aib(1):aib(2),:),&
                        (/ chunkend,na,3 /),ORDER=(/2,1,3/))
                IF (chunksize>chunkend) THEN
                        frame_xyz_rs(chunkend+1:chunksize,:,:) = 0.0_DP
                END IF
                l_moved = .FALSE.
                l_moved(1) = .TRUE.
                !DIR$ ASSUME (MOD(chunksize,4) .EQ. 0)
                DO WHILE (.NOT. ALL(l_moved))
                        DO ia = 1, na
                                IF (.NOT. l_moved(ia)) CYCLE
                                DO ibond = 1 , bondpart_list(ia,is)%nbonds
                                        ja = bondpart_list(ia,is)%atom(ibond)
                                        IF (l_moved(ja)) CYCLE
                                        l_moved(ja) = .TRUE.
                                        DO i_dim = 1, 3
                                                IF (l_ortho) THEN
                                                        boxlen = box_list(ibox)%length(i_dim,i_dim)
                                                ELSE
                                                        boxlen = 1.0_DP
                                                END IF
                                                hl = 0.5_DP * boxlen
                                                !DIR$ VECTOR ALIGNED
                                                !$OMP SIMD PRIVATE(irp,jrp,drp,absdrp)
                                                DO im = 1, chunksize
                                                        irp = frame_xyz_rs(im,ia,i_dim)
                                                        jrp = frame_xyz_rs(im,ja,i_dim)
                                                        drp = jrp - irp
                                                        absdrp = ABS(drp)
                                                        IF (absdrp > hl) jrp = jrp - SIGN(boxlen,drp)
                                                        frame_xyz_rs(im,ja,i_dim) = jrp
                                                END DO
                                                !$OMP END SIMD
                                        END DO
                                END DO
                        END DO
                END DO
                IF (.NOT. l_ortho) THEN
                        h11 = box_list(ibox)%length(1,1)
                        !h21 = box_list(ibox)%length(2,1)
                        !h31 = box_list(ibox)%length(3,1)
                        h12 = box_list(ibox)%length(1,2)
                        h22 = box_list(ibox)%length(2,2)
                        !h32 = box_list(ibox)%length(3,2)
                        h13 = box_list(ibox)%length(1,3)
                        h23 = box_list(ibox)%length(2,3)
                        h33 = box_list(ibox)%length(3,3)
                        DO ia = 1, natoms(is)
                                !DIR$ VECTOR ALIGNED
                                !$OMP SIMD PRIVATE(isp,rxp,ryp,rzp)
                                DO im = 1, chunksize
                                        isp = frame_xyz_rs(im,ia,1)
                                        rxp = h11*isp
                                        isp = frame_xyz_rs(im,ia,2)
                                        rxp = rxp + h12*isp
                                        ryp = h22*isp
                                        isp = frame_xyz_rs(im,ia,3)
                                        rxp = rxp + h13*isp
                                        ryp = ryp + h23*isp
                                        rzp = h33*isp
                                        frame_xyz_rs(im,ia,1) = rxp
                                        frame_xyz_rs(im,ia,2) = ryp
                                        frame_xyz_rs(im,ia,3) = rzp
                                END DO
                                !$OMP END SIMD
                        END DO
                END IF
                inv_total_mass = 1.0_DP/SUM(nonbond_list(1:natoms(is),is)%mass)
                massfrac_vec(1:natoms(is)) = nonbond_list(1:natoms(is),is)%mass*inv_total_mass
                !DIR$ VECTOR ALIGNED
                !$OMP SIMD &
                !$OMP PRIVATE(rxp,ryp,rzp,xcom,ycom,zcom,dxcom,dycom,dzcom) &
                !$OMP PRIVATE(max_dcomsq,dcomsq,massfrac)
                DO im = 1, chunksize
                        massfrac = massfrac_vec(1)
                        rxp = frame_xyz_rs(im,1,1)
                        ryp = frame_xyz_rs(im,1,2)
                        rzp = frame_xyz_rs(im,1,3)
                        xcom = massfrac*rxp
                        ycom = massfrac*ryp
                        zcom = massfrac*rzp
                        DO ia = 2, na
                                massfrac = massfrac_vec(ia)
                                rxp = frame_xyz_rs(im,ia,1)
                                ryp = frame_xyz_rs(im,ia,2)
                                rzp = frame_xyz_rs(im,ia,3)
                                xcom = xcom + massfrac*rxp
                                ycom = ycom + massfrac*ryp
                                zcom = zcom + massfrac*rzp
                        END DO
                        rcom_arr(im,1) = xcom
                        rcom_arr(im,2) = ycom
                        rcom_arr(im,3) = zcom
                        rxp = frame_xyz_rs(im,1,1)
                        ryp = frame_xyz_rs(im,1,2)
                        rzp = frame_xyz_rs(im,1,3)
                        dxcom = rxp - xcom
                        dycom = ryp - ycom
                        dzcom = rzp - zcom
                        max_dcomsq = dxcom*dxcom + dycom*dycom + dzcom*dzcom
                        DO ia = 2, na
                                rxp = frame_xyz_rs(im,ia,1)
                                ryp = frame_xyz_rs(im,ia,2)
                                rzp = frame_xyz_rs(im,ia,3)
                                dxcom = rxp - xcom
                                dycom = ryp - ycom
                                dzcom = rzp - zcom
                                dcomsq = dxcom*dxcom + dycom*dycom + dzcom*dzcom
                                max_dcomsq = MAX(max_dcomsq,dcomsq)
                        END DO
                        rcom_arr(im,4) = SQRT(max_dcomsq)
                END DO
                !$OMP END SIMD
                IF (l_ortho) THEN
                        DO i_dim = 1, 3
                                boxlen = box_list(ibox)%length(i_dim,i_dim)
                                inv_l = 1.0_DP/boxlen
                                !DIR$ VECTOR ALIGNED
                                !$OMP SIMD PRIVATE(rcom,mwv,irp)
                                DO im = 1, chunksize
                                        rcom = rcom_arr(im,i_dim)
                                        mwv = boxlen*ANINT(rcom*inv_l)
                                        !molwrapvec(im,i_dim) = mwv
                                        rcom_arr(im,i_dim) = rcom - mwv
                                        DO ia = 1, natoms(is)
                                                irp = frame_xyz_rs(im,ia,i_dim)
                                                irp = irp - mwv
                                                frame_xyz_rs(im,ia,i_dim) = irp
                                        END DO
                                END DO
                                !$OMP END SIMD
                        END DO
                ELSE
                        !DIR$ VECTOR ALIGNED
                        !$OMP SIMD PRIVATE(rxcom,rycom,rzcom)
                        DO im = 1, chunksize
                                rxcom = rcom_arr(im,1)
                                rycom = rcom_arr(im,2)
                                rzcom = rcom_arr(im,3)
                                molwrapvec(im,1) = rxcom*inv_h11 + rycom*inv_h12 + rzcom*inv_h13
                                molwrapvec(im,2) =                 rycom*inv_h22 + rzcom*inv_h23
                                molwrapvec(im,3) =                                 rzcom*inv_h33
                        END DO
                        !$OMP END SIMD
                        !DIR$ VECTOR ALIGNED
                        !$OMP SIMD PRIVATE(mwv)
                        DO im = 1, chunksize
                                mwv = molwrapvec(im,1)
                                molwrapvec(im,1) = ANINT(mwv)
                                mwv = molwrapvec(im,2)
                                molwrapvec(im,2) = ANINT(mwv)
                                mwv = molwrapvec(im,3)
                                molwrapvec(im,3) = ANINT(mwv)
                        END DO
                        !$OMP END SIMD
                        !molwrapvec(1:nmols_to_read(is,ibox),:) = ANINT(molwrapvec(1:nmols_to_read(is,ibox),:))
                        !DIR$ VECTOR ALIGNED
                        !$OMP SIMD PRIVATE(sxmwv,symwv,szmwv)
                        DO im = 1, chunksize
                                sxmwv = molwrapvec(im,1)
                                symwv = molwrapvec(im,2)
                                szmwv = molwrapvec(im,3)
                                molwrapvec(im,1) = sxmwv*h11 + symwv*h12 + szmwv*h13
                                molwrapvec(im,2) =             symwv*h22 + szmwv*h23
                                molwrapvec(im,3) =                         szmwv*h33
                        END DO
                        !$OMP END SIMD
                        DO i_dim = 1, 3
                                !DIR$ VECTOR ALIGNED
                                !$OMP SIMD PRIVATE(rcom,mwv,irp)
                                DO im = 1, chunksize
                                        rcom = rcom_arr(im,i_dim)
                                        mwv = molwrapvec(im,i_dim)
                                        rcom = rcom - mwv
                                        rcom_arr(im,i_dim) = rcom
                                        DO ia = 1, natoms(is)
                                                irp = frame_xyz_rs(im,ia,i_dim)
                                                irp = irp - mwv
                                                frame_xyz_rs(im,ia,i_dim) = irp
                                        END DO
                                END DO
                                !$OMP END SIMD
                        END DO
                END IF
                locate_base = SUM(nmols(is,1:nbr_boxes)) + chunkshift
                DO im = 1, chunkend
                        loc = im + locate_base
                        molecule_list(loc,is)%live = .TRUE.
                        molecule_list(loc,is)%frac = 1.0_DP
                        molecule_list(loc,is)%which_box = ibox
                        molecule_list(loc,is)%rcom = rcom_arr(im,:)
                        !atom_list(1:natoms(is),im,is)%exist = .TRUE.
                        DO ia = 1, na
                                atom_list(ia,loc,is)%exist = .TRUE.
                                atom_list(ia,loc,is)%rp(1:3) = frame_xyz_rs(im,ia,1:3)
                        END DO
                END DO
        END SUBROUTINE Set_Species_Frame_Coords

        FUNCTION Read_xyz_frame()
                REAL(DP), DIMENSION(natoms_to_read(ibox),3) :: Read_xyz_frame
                INTEGER :: this_unit, io, i
                CHARACTER(6) :: this_element

                this_unit = pregen_xyz_unit(ibox)

                READ(this_unit,*,IOSTAT=io)
                IF (io < 0) THEN
                        early_end = .TRUE.
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



END MODULE Trajectory_Reader_Routines
