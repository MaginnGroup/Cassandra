
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
!*******************************************************************************

MODULE Atompair_nrg_table_routines

  !***************************************************************************
  !
  !
  !********************************************************************************
  USE Global_Variables
  USE Type_Definitions
  USE Energy_Routines , ONLY: AtomPair_VdW_Energy_Vector, Recip_Sqrt
  USE Input_Routines , ONLY: Get_Solvent_Info
  USE Pair_Emax_Estimation, ONLY: Read_Pair_rminsq
  !$ USE OMP_LIB

  IMPLICIT NONE



CONTAINS
        SUBROUTINE Allocate_Atompair_tables
                INTEGER :: is, ia, itype, solute_nextbase, solvent_nextbase, wsolute_nextbase
                INTEGER :: solute_ntypes_counter, solvent_ntypes_counter, wsolute_ntypes_counter
                LOGICAL, DIMENSION(0:nbr_atomtypes) :: l_type_solute, l_type_solvent, l_type_wsolute
                IF (need_solvents) CALL Get_Solvent_Info
                solute_nextbase = 0
                wsolute_nextbase = 0
                solvent_nextbase = 0
                solute_ntypes_counter = 0
                wsolute_ntypes_counter = 0
                solvent_ntypes_counter = 0
                l_type_solute = .FALSE.
                l_type_solvent = .FALSE.
                l_type_wsolute = .FALSE.
                IF (ALLOCATED(atompair_nrg_table)) DEALLOCATE(atompair_nrg_table)
                IF (ALLOCATED(typepair_nrg_table)) DEALLOCATE(typepair_nrg_table)
                IF (ALLOCATED(atompair_rminsq_table)) DEALLOCATE(atompair_rminsq_table)
                IF (ALLOCATED(typepair_solute_indices)) DEALLOCATE(typepair_solute_indices)
                IF (ALLOCATED(typepair_wsolute_indices)) DEALLOCATE(typepair_wsolute_indices)
                IF (ALLOCATED(typepair_solvent_indices)) DEALLOCATE(typepair_solute_indices)
                ALLOCATE(typepair_solute_indices(0:nbr_atomtypes))
                ALLOCATE(typepair_wsolute_indices(0:nbr_atomtypes))
                ALLOCATE(typepair_solvent_indices(0:nbr_atomtypes))
                typepair_solute_indices = 0
                typepair_wsolute_indices = 0
                typepair_solvent_indices = 0
                DO is = 1, nspecies
                        IF (species_list(is)%l_wsolute) THEN
                                species_list(is)%wsolute_base = wsolute_nextbase
                                wsolute_nextbase = wsolute_nextbase + natoms(is)
                                ! slicing nonbond_list would be improper here because of atom type duplicates
                                DO ia = 1, natoms(is)
                                        l_type_wsolute(nonbond_list(ia,is)%atom_type_number) = .TRUE.
                                END DO
                        ELSE
                                species_list(is)%wsolute_base = -99999 ! intended to cause runtime error if used in table reference
                        END IF
                        IF (species_list(is)%l_solute) THEN
                                species_list(is)%solute_base = solute_nextbase
                                solute_nextbase = solute_nextbase + natoms(is)
                                ! slicing nonbond_list would be improper here because of atom type duplicates
                                DO ia = 1, natoms(is)
                                        l_type_solute(nonbond_list(ia,is)%atom_type_number) = .TRUE.
                                END DO
                        ELSE
                                species_list(is)%solute_base = -99999 ! intended to cause runtime error if used in table reference
                        END IF
                        IF (species_list(is)%l_solvent) THEN
                                species_list(is)%solvent_base = solvent_nextbase
                                solvent_nextbase = solvent_nextbase + natoms(is)
                                ! slicing nonbond_list would be improper here because of atom type duplicates
                                DO ia = 1, natoms(is)
                                        l_type_solvent(nonbond_list(ia,is)%atom_type_number) = .TRUE.
                                END DO
                        ELSE
                                species_list(is)%solvent_base = -99999 ! intended to cause runtime error if used in table reference
                        END IF
                END DO
                solute_ntypes = COUNT(l_type_solute(1:))
                wsolute_ntypes = COUNT(l_type_wsolute(1:))
                solvent_ntypes = COUNT(l_type_solvent(1:))
                IF (ALLOCATED(solute_atomtypes)) DEALLOCATE(solute_atomtypes)
                IF (ALLOCATED(wsolute_atomtypes)) DEALLOCATE(wsolute_atomtypes)
                IF (ALLOCATED(solvent_atomtypes)) DEALLOCATE(solvent_atomtypes)
                ALLOCATE(solute_atomtypes(solute_ntypes))
                ALLOCATE(wsolute_atomtypes(wsolute_ntypes))
                ALLOCATE(solvent_atomtypes(solvent_ntypes))
                DO itype = 1, nbr_atomtypes
                        IF (l_type_wsolute(itype)) THEN
                                wsolute_ntypes_counter = wsolute_ntypes_counter + 1
                                typepair_wsolute_indices(itype) = wsolute_ntypes_counter
                                wsolute_atomtypes(wsolute_ntypes_counter) = itype
                        END IF
                        IF (l_type_solute(itype)) THEN
                                solute_ntypes_counter = solute_ntypes_counter + 1
                                typepair_solute_indices(itype) = solute_ntypes_counter
                                solute_atomtypes(solute_ntypes_counter) = itype
                        END IF
                        IF (l_type_solvent(itype)) THEN
                                solvent_ntypes_counter = solvent_ntypes_counter + 1
                                typepair_solvent_indices(itype) = solvent_ntypes_counter
                                solvent_atomtypes(solvent_ntypes_counter) = itype
                        END IF
                END DO
                solvent_maxind = solvent_nextbase
                solute_maxind = solute_nextbase
                wsolute_maxind = wsolute_nextbase
                IF (precalc_atompair_nrg) THEN
                        ALLOCATE(typepair_nrg_table(atompair_nrg_res,0:solvent_ntypes,0:solute_ntypes,nbr_boxes))
                        ALLOCATE(atompair_nrg_table(atompair_nrg_res+1,solvent_nextbase,solute_nextbase,nbr_boxes))
                        ALLOCATE(atompair_nrg_table_reduced(0:(atompair_nrg_res+1)*solvent_nextbase-1,solute_nextbase,nbr_boxes))
                        typepair_nrg_table = 0.0_DP
                END IF
                IF (est_atompair_rminsq) THEN
                        maxrminsq = rsqmin_step * rsqmin_res + rcut_lowsq
                        rsqmin_shifter = rcut_lowsq - rsqmin_step
                        ALLOCATE(rsqmin_atompair_w_max(rsqmin_res,solvent_nextbase,wsolute_nextbase,nbr_boxes))
                        ALLOCATE(rsqmin_atompair_w_sum(rsqmin_res,solvent_nextbase,wsolute_nextbase,nbr_boxes))
                        ALLOCATE(rsqmin_atompair_freq(rsqmin_res,solvent_maxind,wsolute_maxind,nbr_boxes))
                        rsqmin_atompair_w_max = 0.0_DP
                        rsqmin_atompair_w_sum = 0.0_DP
                        rsqmin_atompair_freq = 0_INT64
                END IF
                IF (read_atompair_rminsq) THEN 
                        ALLOCATE(atompair_rminsq_table(solvent_nextbase,wsolute_nextbase,nbr_boxes))
                        ALLOCATE(sp_atompair_rminsq_table(solvent_nextbase,wsolute_nextbase,nbr_boxes))
                END IF
        END SUBROUTINE Allocate_Atompair_tables


        SUBROUTINE Create_Atompair_Nrg_table
                INTEGER :: ibox, i, is
                REAL(DP), DIMENSION(:,:,:), POINTER :: atompair_nrg_ptr
                !REAL(DP), DIMENSION(:,:), POINTER :: atompair_rminsq_ptr
                INTEGER, DIMENSION(nspecies) :: solutes, solvents
                INTEGER :: nsolutes, nsolvents
                INTEGER :: ti_solute, ti_solvent, ti_type_solute, ti_type_solvent
                REAL(DP), DIMENSION(solvent_maxind, solute_maxind) :: cfqq
                REAL(DP), DIMENSION(atompair_nrg_res, nbr_boxes) :: f2
                REAL(DP) :: solvent_charges(solvent_maxind), solute_charges(solute_maxind)
                INTEGER :: solvent_typeindvec(solvent_maxind), solute_typeindvec(solute_maxind)
                REAL(DP), DIMENSION(atompair_nrg_res) :: rsq_mp_vector, rsq_lb_vector, rij, inv_rij, alpha_rij
                IF (.NOT. precalc_atompair_nrg) RETURN
                nsolutes = 0
                nsolvents = 0
                rsq_step = (MAXVAL(rcut_cbmcsq)-rcut_lowsq)/atompair_nrg_res
                inv_rsq_step = 1.0_DP/rsq_step
                inv_rsq_step_sp = REAL(inv_rsq_step,SP)
                rsq_shifter = rcut_lowsq - rsq_step
                DO i = 1, atompair_nrg_res
                        rsq_lb_vector(i) = rsq_shifter + rsq_step*i
                END DO
                rsq_mp_vector = rsq_lb_vector + 0.5*rsq_step
                DO is = 1, nspecies
                        IF (species_list(is)%l_solute) THEN
                                nsolutes = nsolutes + 1
                                solutes(nsolutes) = is
                        END IF
                        IF (species_list(is)%l_solvent) THEN
                                nsolvents = nsolvents + 1
                                solvents(nsolvents) = is
                        END IF
                END DO

                DO i = 1, nsolvents
                        is = solvents(i)
                        solvent_charges( &
                                species_list(is)%solvent_base+1:species_list(is)%solvent_base+natoms(is)) &
                                = nonbond_list(1:natoms(is),is)%charge
                        solvent_typeindvec( &
                                species_list(is)%solvent_base+1:species_list(is)%solvent_base+natoms(is)) &
                                = typepair_solvent_indices(nonbond_list(1:natoms(is),is)%atom_type_number)
                END DO
                DO i = 1, nsolutes
                        is = solutes(i)
                        solute_charges( &
                                species_list(is)%solute_base+1:species_list(is)%solute_base+natoms(is)) &
                                = nonbond_list(1:natoms(is),is)%charge
                        solute_typeindvec( &
                                species_list(is)%solute_base+1:species_list(is)%solute_base+natoms(is)) &
                                = typepair_solute_indices(nonbond_list(1:natoms(is),is)%atom_type_number)
                END DO

                inv_rij = Recip_Sqrt(rsq_mp_vector)
                rij = inv_rij*rsq_mp_vector

                DO ibox = 1, nbr_boxes
                        IF (cbmc_charge_sf_flag) THEN
                                f2(:,ibox) = inv_rij - 2.0_DP/rcut_cbmc(ibox) + rij/rcut_cbmcsq(ibox)
                        ELSEIF (int_charge_sum_style(ibox) == charge_ewald) THEN
                                alpha_rij = alpha_ewald(ibox) * rij
                                f2(:,ibox) = ERFC(alpha_rij)*inv_rij
                        ELSEIF (int_charge_sum_style(ibox) == charge_dsf) THEN
                                alpha_rij = alpha_dsf(ibox)*rij
                                f2(:,ibox) = &
                                        dsf_factor2(ibox)*(rij-rcut_coul(ibox)) - &
                                        dsf_factor1(ibox) + &
                                        ERFC(alpha_rij)*inv_rij
                        ELSEIF (int_charge_sum_style(ibox) == charge_cut) THEN
                                f2(:,ibox) = inv_rij
                        ELSE
                                f2(:,ibox) = 0.0_DP
                        END IF
                END DO



                !$OMP PARALLEL DEFAULT(SHARED)
                !$OMP DO SCHEDULE(STATIC)
                DO ti_solute = 1, solute_maxind
                        cfqq(:,ti_solute) = (solute_charges(ti_solute)*charge_factor)*solvent_charges
                END DO
                !$OMP END DO
                !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
                DO ibox = 1, nbr_boxes
                DO ti_type_solvent = 1, solvent_ntypes
                        DO ti_type_solute = 1, solute_ntypes
                                typepair_nrg_table(:,ti_type_solvent,ti_type_solute,ibox) = &
                                        AtomPair_VdW_Energy_Vector(rsq_mp_vector, &
                                                solvent_atomtypes(ti_type_solvent), &
                                                solute_atomtypes(ti_type_solute), &
                                                ibox)
                        END DO
                END DO
                END DO
                !$OMP END DO
                !$OMP END PARALLEL

                !$OMP WORKSHARE
                atompair_nrg_table(1:atompair_nrg_res,:,:,:) = typepair_nrg_table(:,solvent_typeindvec,solute_typeindvec,:)
                !$OMP END WORKSHARE

                !$OMP PARALLEL DEFAULT(SHARED)
                !$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
                DO ibox = 1, nbr_boxes
                        DO ti_solute = 1, solute_maxind
                                DO ti_solvent = 1, solvent_maxind
                                        atompair_nrg_table(1:atompair_nrg_res,ti_solvent,ti_solute,ibox) = &
                                                atompair_nrg_table(1:atompair_nrg_res,ti_solvent,ti_solute,ibox) + &
                                                f2(:,ibox)*cfqq(ti_solvent,ti_solute)
                                END DO
                        END DO
                END DO
                !$OMP END DO
                !$OMP WORKSHARE
                atompair_nrg_table(atompair_nrg_res+1,:,:,:) = 0.0
                atompair_nrg_table_reduced = REAL(RESHAPE(atompair_nrg_table, SHAPE(atompair_nrg_table_reduced)),SP)
                !$OMP END WORKSHARE
                !$OMP END PARALLEL



        END SUBROUTINE Create_Atompair_Nrg_table

        SUBROUTINE Setup_Atompair_tables
                solvent_maxind = 1
                rsqmin_res_d = 1
                solvent_maxind_d = 1
                IF (.NOT. (precalc_atompair_nrg .OR. read_atompair_rminsq .OR. est_atompair_rminsq)) RETURN
                CALL Allocate_Atompair_Tables
                IF (precalc_atompair_nrg) CALL Create_Atompair_Nrg_table
                IF (read_atompair_rminsq) CALL Read_Pair_rminsq
                IF (est_atompair_rminsq .AND. .NOT. l_heap) THEN
                        rsqmin_res_d = rsqmin_res
                        solvent_maxind_d = solvent_maxind
                END IF
        END SUBROUTINE





END MODULE Atompair_nrg_table_routines
