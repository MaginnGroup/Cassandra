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


END MODULE Pair_Emax_Estimation
