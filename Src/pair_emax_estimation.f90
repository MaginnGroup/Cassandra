!*****************************************************************************
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
!*****************************************************************************

MODULE Pair_Emax_Estimation

  !***************************************************************************
  !***************************************************************************

  USE Global_Variables
  USE IO_Utilities
  USE File_Names
  USE Type_Definitions
  !$ USE OMP_LIB


  IMPLICIT NONE

  INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE, TARGET ::atompair_rminsq_ind_table
  INTEGER(KIND=INT64), DIMENSION(:,:,:), ALLOCATABLE :: rminsq_nskips
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: rminsq_wmax

CONTAINS
        SUBROUTINE Estimate_Pair_Emax(is,ibox,Emax,nskips)
                INTEGER :: is, ibox
                INTEGER :: nskips
                REAL(DP) :: Emax

                INTEGER :: Eij_ind
                REAL(DP) :: Eij_mult, Eij_wfrac_csum
                REAL(DP), DIMENSION(Eij_ind_ubound+1) :: Eij_wfrac

                Eij_wfrac_csum = 0.0_DP
                Eij_mult = beta(ibox) / Eij_factor(is,ibox)
                Eij_wfrac = Eij_w_sum(:,is,ibox) / species_list(is)%widom_sum(ibox)

                DO Eij_ind = Eij_ind_ubound+1, 1, -1
                        Eij_wfrac_csum = Eij_wfrac_csum + Eij_wfrac(Eij_ind)
                        IF (Eij_wfrac_csum > 1.0e-10_DP) EXIT
                END DO

                Emax = Eij_ind * Eij_mult
                IF (Eij_ind <= Eij_ind_ubound) THEN
                        nskips = SUM(Eij_freq_total((Eij_ind+1):,is,ibox))
                ELSE
                        nskips = 0
                END IF

                OPEN(unit=emax_file_unit, file=emax_filenames(is,ibox))
                WRITE(emax_file_unit, '(A9,7X,A20,7X,A20,7X,A20,7X,A20)') 'Eij_index', 'max_Eij', 'w_frac', 'w_max', 'w_freq'
                DO Eij_ind = 1, Eij_ind_ubound+1
                        IF (Eij_freq_total(Eij_ind,is,ibox) == 0) CYCLE
                        WRITE(emax_file_unit, '(I9,7X,F20.6,7X,E20.12,7X,E20.12,7X,I20)') Eij_ind, &
                                Eij_ind*Eij_mult, Eij_wfrac(Eij_ind), w_max(Eij_ind,is,ibox), &
                                Eij_freq_total(Eij_ind,is,ibox)
                END DO
                CLOSE(emax_file_unit)


        END SUBROUTINE Estimate_Pair_Emax

        SUBROUTINE Estimate_Pair_rminsq
                INTEGER :: is, ibox, i_tol
                REAL(DP) :: this_tol

                INTEGER :: rsqmin_ind, bsolute, ti_solute, ti_solvent
                REAL(DP) :: wfrac_csum, wmax
                INTEGER(KIND=INT64) :: freq_csum
                REAL(DP), DIMENSION(:), POINTER :: rsqmin_wfrac_ptr
                REAL(DP), DIMENSION(:), POINTER :: rsqmin_wmax_ptr
                REAL(DP), DIMENSION(rsqmin_res,solvent_maxind,wsolute_maxind,nbr_boxes), TARGET :: rsqmin_atompair_wfrac
                INTEGER(KIND=INT64), DIMENSION(:), POINTER :: rsqmin_freq_ptr
                INTEGER(KIND=INT64), DIMENSION(solvent_maxind,wsolute_maxind,nbr_boxes,nbr_tols) :: atompair_nskips
                REAL(DP), DIMENSION(solvent_maxind,wsolute_maxind,nbr_boxes,nbr_tols) :: atompair_wmax
                CHARACTER(FILENAME_LEN) :: this_path
                INTEGER :: this_unit

                IF (.NOT. ALLOCATED(atompair_rminsq_ind_table)) THEN
                        ALLOCATE(atompair_rminsq_ind_table(solvent_maxind,wsolute_maxind,nbr_boxes,nbr_tols))
                END IF
                IF (ALLOCATED(rminsq_nskips)) DEALLOCATE(rminsq_nskips)
                ALLOCATE(rminsq_nskips(nspecies,nbr_boxes,nbr_tols))
                IF (ALLOCATED(rminsq_wmax)) DEALLOCATE(rminsq_wmax)
                ALLOCATE(rminsq_wmax(nspecies,nbr_boxes,nbr_tols))

                DO ibox = 1, nbr_boxes
                        DO is = 1, nspecies
                                IF (.NOT. species_list(is)%test_particle(ibox)) CYCLE
                                bsolute = species_list(is)%wsolute_base
                                rsqmin_atompair_wfrac(:,:,bsolute+1:bsolute+natoms(is),ibox) = &
                                        rsqmin_atompair_w_sum(:,:,bsolute+1:bsolute+natoms(is),ibox) / &
                                        species_list(is)%widom_sum(ibox)
                        END DO
                END DO

                !$OMP PARALLEL DEFAULT(SHARED) &
                !$OMP PRIVATE(rsqmin_wfrac_ptr, wfrac_csum, rsqmin_ind, freq_csum, rsqmin_freq_ptr, this_tol) &
                !$OMP PRIVATE(wmax, rsqmin_wmax_ptr)
                !$OMP DO COLLAPSE(4) SCHEDULE(STATIC)
                DO ibox = 1, nbr_boxes
                        DO ti_solute = 1, wsolute_maxind
                                DO ti_solvent = 1, solvent_maxind
                                        DO i_tol = 1, nbr_tols
                                                rsqmin_wfrac_ptr => rsqmin_atompair_wfrac(:,ti_solvent,ti_solute,ibox)
                                                rsqmin_freq_ptr => rsqmin_atompair_freq(:,ti_solvent,ti_solute,ibox)
                                                rsqmin_wmax_ptr => rsqmin_atompair_w_max(:,ti_solvent,ti_solute,ibox)
                                                this_tol = tol_list(i_tol)
                                                wmax = 0.0_DP
                                                wfrac_csum = 0.0_DP
                                                freq_csum = 0_INT64
                                                DO rsqmin_ind = 1, rsqmin_res
                                                        wfrac_csum = wfrac_csum + rsqmin_wfrac_ptr(rsqmin_ind)
                                                        IF (wfrac_csum > this_tol) EXIT
                                                        IF (rsqmin_wmax_ptr(rsqmin_ind)>wmax) wmax = rsqmin_wmax_ptr(rsqmin_ind)
                                                        freq_csum = freq_csum + rsqmin_freq_ptr(rsqmin_ind)
                                                END DO
                                                atompair_rminsq_ind_table(ti_solvent,ti_solute,ibox,i_tol) = rsqmin_ind
                                                atompair_nskips(ti_solvent,ti_solute,ibox,i_tol) = freq_csum
                                                atompair_wmax(ti_solvent,ti_solute,ibox,i_tol) = wmax
                                        END DO
                                END DO
                        END DO
                END DO
                !$OMP END DO
                !$OMP END PARALLEL
                this_unit = emax_file_unit
                this_path = ""
                this_path = TRIM(run_name) // ".rsqmin.wmax"
                OPEN(unit=this_unit,file=this_path,ACTION='WRITE')
                DO ibox = 1, nbr_boxes
                        WRITE(this_unit,*) "Box ", ibox
                        DO ti_solute = 1, wsolute_maxind
                                WRITE(this_unit,*)
                                WRITE(this_unit,*) "Solute atom ", ti_solute
                                DO ti_solvent = 1, solvent_maxind
                                        WRITE(this_unit,*) rsqmin_atompair_w_max(:,ti_solvent,ti_solute,ibox)
                                END DO
                        END DO
                END DO
                WRITE(this_unit,*) "END"
                CLOSE(this_unit)
                DO ibox = 1, nbr_boxes
                        DO is = 1, nspecies
                                IF (.NOT. species_list(is)%test_particle(ibox)) CYCLE
                                bsolute = species_list(is)%wsolute_base
                                rsqmin_atompair_w_max(:,:,bsolute+1:bsolute+natoms(is),ibox) = &
                                        rsqmin_atompair_w_max(:,:,bsolute+1:bsolute+natoms(is),ibox) / &
                                        species_list(is)%widom_sum(ibox)
                        END DO
                END DO
                this_path = ""
                this_path = TRIM(run_name) // ".rsqmin.wmax_frac"
                OPEN(unit=this_unit,file=this_path,ACTION='WRITE')
                DO ibox = 1, nbr_boxes
                        WRITE(this_unit,*) "Box ", ibox
                        DO ti_solute = 1, wsolute_maxind
                                WRITE(this_unit,*)
                                WRITE(this_unit,*) "Solute atom ", ti_solute
                                DO ti_solvent = 1, solvent_maxind
                                        WRITE(this_unit,*) rsqmin_atompair_w_max(:,ti_solvent,ti_solute,ibox)
                                END DO
                        END DO
                END DO
                WRITE(this_unit,*) "END"
                CLOSE(this_unit)
                DO ibox = 1, nbr_boxes
                        DO is = 1, nspecies
                                IF (.NOT. species_list(is)%test_particle(ibox)) CYCLE
                                bsolute = species_list(is)%wsolute_base
                                rsqmin_atompair_w_max(:,:,bsolute+1:bsolute+natoms(is),ibox) = &
                                        rsqmin_atompair_w_max(:,:,bsolute+1:bsolute+natoms(is),ibox) * &
                                        ntrials(is,ibox)%widom
                        END DO
                END DO
                this_path = ""
                this_path = TRIM(run_name) // ".rsqmin.wmax_norm"
                OPEN(unit=this_unit,file=this_path,ACTION='WRITE')
                DO ibox = 1, nbr_boxes
                        WRITE(this_unit,*) "Box ", ibox
                        DO ti_solute = 1, wsolute_maxind
                                WRITE(this_unit,*)
                                WRITE(this_unit,*) "Solute atom ", ti_solute
                                DO ti_solvent = 1, solvent_maxind
                                        WRITE(this_unit,*) rsqmin_atompair_w_max(:,ti_solvent,ti_solute,ibox)
                                END DO
                        END DO
                END DO
                WRITE(this_unit,*) "END"
                CLOSE(this_unit)
                this_path = ""
                this_path = TRIM(run_name) // ".rsqmin.wfrac"
                OPEN(unit=this_unit,file=this_path,ACTION='WRITE')
                DO ibox = 1, nbr_boxes
                        WRITE(this_unit,*) "Box ", ibox
                        DO ti_solute = 1, wsolute_maxind
                                WRITE(this_unit,*)
                                WRITE(this_unit,*) "Solute atom ", ti_solute
                                DO ti_solvent = 1, solvent_maxind
                                        WRITE(this_unit,*) rsqmin_atompair_wfrac(:,ti_solvent,ti_solute,ibox)
                                END DO
                        END DO
                END DO
                WRITE(this_unit,*) "END"
                CLOSE(this_unit)
                this_path = ""
                this_path = TRIM(run_name) // ".rsqmin.freq"
                OPEN(unit=this_unit,file=this_path,ACTION='WRITE')
                DO ibox = 1, nbr_boxes
                        WRITE(this_unit,*) "Box ", ibox
                        DO ti_solute = 1, wsolute_maxind
                                WRITE(this_unit,*)
                                WRITE(this_unit,*) "Solute atom ", ti_solute
                                DO ti_solvent = 1, solvent_maxind
                                        WRITE(this_unit,*) rsqmin_atompair_freq(:,ti_solvent,ti_solute,ibox)
                                END DO
                        END DO
                END DO
                WRITE(this_unit,*) "END"
                CLOSE(this_unit)


                DO ibox = 1, nbr_boxes
                        DO is = 1, nspecies
                                IF (.NOT. species_list(is)%test_particle(ibox)) CYCLE
                                bsolute = species_list(is)%wsolute_base
                                DO i_tol = 1, nbr_tols
                                        rminsq_nskips(is,ibox,i_tol) = MAXVAL(atompair_nskips(:, &
                                                bsolute+1:bsolute+natoms(is), ibox, i_tol))
                                        rminsq_wmax(is,ibox,i_tol) = MAXVAL(atompair_wmax(:, &
                                                bsolute+1:bsolute+natoms(is), ibox, i_tol))
                                END DO
                        END DO
                END DO

        END SUBROUTINE Estimate_Pair_rminsq

        FUNCTION Get_rminsq_nskips(ispecies,ibox,i_tol)
                INTEGER :: get_rminsq_nskips, ispecies, ibox, i_tol
                IF (.NOT. ALLOCATED(rminsq_nskips)) CALL Estimate_Pair_rminsq
                get_rminsq_nskips = rminsq_nskips(ispecies,ibox,i_tol)
        END FUNCTION Get_rminsq_nskips

        FUNCTION Get_rminsq_wmax(ispecies,ibox,i_tol)
                REAL(DP) :: get_rminsq_wmax
                INTEGER :: ispecies, ibox, i_tol
                IF (.NOT. ALLOCATED(rminsq_wmax)) CALL Estimate_Pair_rminsq
                get_rminsq_wmax = rminsq_wmax(ispecies,ibox,i_tol)
        END FUNCTION Get_rminsq_wmax

        SUBROUTINE Write_Pair_rminsq
                INTEGER :: this_unit, nchars, ti_solvent, ti_solute, ibox, i_tol
                CHARACTER(20) :: fstring
                CHARACTER(FILENAME_LEN) :: this_path
                nchars = INT(LOG10(REAL(rsqmin_res))) + 2
                this_unit = emax_file_unit
                this_path = ""
                fstring = ""
                fstring = "(" // TRIM(Int_To_String(wsolute_maxind)) // "I" // TRIM(Int_To_String(nchars)) // ")"
                this_path = TRIM(write_rminsq_filename)
                DO i_tol = 1, nbr_tols
                        IF (nbr_tols > 1) this_path = TRIM(write_rminsq_filename) // TRIM(Int_To_String(i_tol))
                        OPEN(unit=this_unit,file=this_path,ACTION='WRITE')
                        WRITE(this_unit,*) tol_list(i_tol)
                        WRITE(this_unit,*) rcut_lowsq, rsqmin_step, rsqmin_res
                        WRITE(this_unit,*) solvent_maxind, wsolute_maxind, nbr_boxes
                        DO ibox = 1, nbr_boxes
                                DO ti_solvent = 1, solvent_maxind
                                        WRITE(this_unit, fstring) &
                                                (atompair_rminsq_ind_table(ti_solvent,ti_solute,ibox,i_tol),ti_solute=1,wsolute_maxind)
                                        !WRITE(this_unit,"(I" // TRIM(Int_To_String(nchars)) // ")") &
                                        !        (atompair_rminsq_ind_table(ti_solvent,ti_solute,ibox),ti_solute=1,wsolute_maxind)
                                END DO
                        END DO
                        CLOSE(this_unit)
                END DO
        END SUBROUTINE Write_Pair_rminsq

        SUBROUTINE Read_Pair_rminsq
                INTEGER :: this_unit, ti_solvent, ti_solute, ibox, this_rsqmin_res
                REAL(DP) :: this_rcut_lowsq, this_rsqmin_step, this_rsqmin_shifter, this_tol, min_rmin
                INTEGER :: this_solvent_maxind, this_wsolute_maxind, this_nbr_boxes
                INTEGER, DIMENSION(:), POINTER :: atompair_rminsq_ind_table_ptr
                INTEGER, DIMENSION(wsolute_maxind) :: ti_which_big_atom, big_atom_ti_list
                INTEGER :: ifrag, ia_frag, ia_frag_ti, biggest_atom, biggest_atom_ti, i_big_atom, is, wsb
                REAL(DP) :: ia_frag_rminsq_sum, biggest_atom_rminsq_sum




                this_unit = emax_file_unit
                IF (.NOT. ALLOCATED(atompair_rminsq_ind_table)) THEN
                        ALLOCATE(atompair_rminsq_ind_table(solvent_maxind,wsolute_maxind,nbr_boxes,nbr_tols))
                END IF
                OPEN(unit=this_unit,file=read_rminsq_filename,POSITION="REWIND",ACTION='READ',STATUS='OLD')
                WRITE(logunit,'(A80)') "********************************************************************************"
                READ(this_unit,*) this_tol
                WRITE(logunit,'(A,ES10.3)') "Reading atompair rminsq table made with tolerance = ", this_tol
                READ(this_unit,*) this_rcut_lowsq, this_rsqmin_step, this_rsqmin_res
                this_rsqmin_shifter = this_rcut_lowsq - this_rsqmin_step
                READ(this_unit,*) this_solvent_maxind, this_wsolute_maxind, this_nbr_boxes
                IF (this_solvent_maxind .NE. solvent_maxind .OR. &
                        this_wsolute_maxind .NE. wsolute_maxind .OR. &
                        this_nbr_boxes .NE. nbr_boxes) THEN
                        err_msg = ""
                        err_msg(1) = "atompair rminsq array must have dimensions"
                        err_msg(2) = "(" // Int_To_String(solvent_maxind) // &
                                ", " // Int_To_String(wsolute_maxind) // &
                                ", " // Int_To_String(nbr_boxes) // ")"
                        err_msg(3) = "but dimensions provided were"
                        err_msg(4) = "(" // Int_To_String(this_solvent_maxind) // &
                                ", " // Int_To_String(this_wsolute_maxind) // &
                                ", " // Int_To_String(this_nbr_boxes) // ")"
                        CALL Clean_Abort(err_msg, "Read_Pair_rminsq")
                END IF
                DO ibox = 1, nbr_boxes
                        DO ti_solvent = 1, solvent_maxind
                                atompair_rminsq_ind_table_ptr => &
                                        atompair_rminsq_ind_table(ti_solvent,:,ibox,1)
                                READ(this_unit,*) &
                                        (atompair_rminsq_ind_table_ptr(ti_solute),ti_solute=1,wsolute_maxind)
                        END DO
                END DO
                CLOSE(this_unit)
                max_rmin = SQRT(REAL(MAXVAL(atompair_rminsq_ind_table(:,:,:,1), &
                        atompair_rminsq_ind_table(:,:,:,1)<(this_rsqmin_res+1)),DP) &
                        * this_rsqmin_step + this_rsqmin_shifter)
                atompair_rminsq_table = REAL(atompair_rminsq_ind_table(:,:,:,1),DP) * this_rsqmin_step + this_rsqmin_shifter
                min_rmin = SQRT(MINVAL(atompair_rminsq_table))
                n_big_atoms = 0
                IF (cavity_biasing_flag) THEN
                        ti_which_big_atom = -999999
                        DO is = 1, nspecies
                                IF (.NOT. species_list(is)%l_wsolute) CYCLE
                                wsb = species_list(is)%wsolute_base
                                DO ifrag = 1, nfragments(is)
                                        biggest_atom = 1
                                        biggest_atom_ti = frag_list(ifrag,is)%atoms(1) + wsb
                                        biggest_atom_rminsq_sum = SUM(atompair_rminsq_table(:,biggest_atom_ti,:))
                                        DO ia_frag = 2, frag_list(ifrag,is)%natoms
                                                ia_frag_ti = frag_list(ifrag,is)%atoms(ia_frag) + wsb
                                                ia_frag_rminsq_sum = SUM(atompair_rminsq_table(:,ia_frag_ti,:))
                                                IF (ia_frag_rminsq_sum>biggest_atom_rminsq_sum) THEN
                                                        biggest_atom = ia_frag
                                                        biggest_atom_ti = ia_frag_ti
                                                        biggest_atom_rminsq_sum = ia_frag_rminsq_sum
                                                END IF
                                        END DO
                                        i_big_atom = ti_which_big_atom(biggest_atom_ti)
                                        IF (i_big_atom < 0) THEN
                                                n_big_atoms = n_big_atoms+1
                                                i_big_atom = n_big_atoms
                                                ti_which_big_atom(biggest_atom_ti) = n_big_atoms
                                                big_atom_ti_list(n_big_atoms) = biggest_atom_ti
                                        END IF
                                        frag_list(ifrag,is)%i_big_atom = i_big_atom
                                        frag_list(ifrag,is)%ia_frag_big_atom = biggest_atom
                                END DO
                        END DO
                END IF
                ALLOCATE(solvent_max_rminsq(solvent_maxind,nbr_boxes))
                ALLOCATE(solvent_min_rminsq(solvent_maxind,0:n_big_atoms,nbr_boxes))
                ALLOCATE(solvent_max_rminsq_sp(solvent_maxind,nbr_boxes))
                solvent_max_rminsq = MAXVAL(atompair_rminsq_table,2) ! 2-D result
                solvent_min_rminsq(:,0,:) = MINVAL(atompair_rminsq_table,2) ! 2-D result
                IF (cavity_biasing_flag) solvent_min_rminsq(:,1:n_big_atoms,:) = atompair_rminsq_table(:,big_atom_ti_list(1:n_big_atoms),:)
                box_list%rcut_low_max = SQRT(MAXVAL(MAXVAL(solvent_min_rminsq,1),1))
                box_list%ideal_bitcell_length = &
                        MAX(min_ideal_bitcell_length,box_list%rcut_low_max/SQRT(902.0_DP)) ! vector with one element per box
                sp_atompair_rminsq_table = REAL(atompair_rminsq_table,SP)
                solvent_max_rminsq_sp = REAL(solvent_max_rminsq,SP)
                solvents_or_types_maxind = solvent_maxind
                WRITE(logunit,'(A,F6.3,A)') "Finished reading atompair rminsq table with atompair rmin values in the interval"
                WRITE(logunit, '(8x,A,F5.3,A,F6.3,A)') "[", min_rmin, ",", max_rmin, "] Angstroms"
                IF (bitcell_flag) THEN
                        DO ibox = 1, nbr_boxes
                                WRITE(logunit, '(A,F5.3,A)') "For box " // TRIM(Int_To_String(ibox)) // ", computed ideal bitcell length = ", &
                                        box_list(ibox)%rcut_low_max/SQRT(902.0_DP), " Angstroms"
                                WRITE(logunit, '(A,F5.3,A)') "Setting box " // TRIM(Int_To_String(ibox)) // " ideal bitcell length = ", &
                                        box_list(ibox)%ideal_bitcell_length, " Angstroms"
                        END DO
                END IF
                WRITE(logunit,'(A80)') "********************************************************************************"
        END SUBROUTINE Read_Pair_rminsq


END MODULE Pair_Emax_Estimation
