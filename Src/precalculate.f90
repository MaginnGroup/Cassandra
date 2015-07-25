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
!********************************************************************************

SUBROUTINE precalculate
!*******************************************************************************  
    ! Revision history:
    !
    ! 12/10/13  : Beta version 
!********************************************************************************    
    USE Run_variables
    USE File_Names
    USE Energy_Routines
    IMPLICIT NONE
    
    
    INTEGER :: i, ibox
    ! Determine the direct sum maximum cutoff distances squared, beyond 
    ! which we ingore pairwise interactions. 
    ! Will need to alter when neighborlist added.
    
    REAL(DP) :: roffsq_ronsq


    DO ibox = 1,nbr_boxes
       IF (int_vdw_style(ibox) /= vdw_none) THEN
          IF (int_vdw_sum_style(ibox) == vdw_charmm) THEN
             ! Compute the square for use in pair routines
             ron_charmmsq(ibox) = ron_charmm(ibox)*ron_charmm(ibox)
             roff_charmmsq(ibox) = roff_charmm(ibox)*roff_charmm(ibox)

          ELSE IF (int_vdw_sum_style(ibox) == vdw_cut_switch) THEN

             ron_switch_sq(ibox) = ron_switch(ibox) * ron_switch(ibox)

             roff_switch_sq(ibox) = roff_switch(ibox) * roff_switch(ibox)

             roffsq_ronsq = roff_switch_sq(ibox) - ron_switch_sq(ibox)

             switch_factor1(ibox) = roffsq_ronsq ** 3.0_DP
             switch_factor1(ibox) = 1.0_DP /switch_factor1(ibox)

             switch_factor2(ibox) = roff_switch_sq(ibox) - 3.0_DP * ron_switch_sq(ibox)

             rcut_vdw(ibox) = roff_switch(ibox)

          
          ELSE
             ! Compute the square for use in pair routines
             rcut_vdwsq(ibox) = rcut_vdw(ibox)*rcut_vdw(ibox)
             rcut_vdw3(ibox) = rcut_vdwsq(ibox) * rcut_vdw(ibox)
             rcut_vdw6(ibox) = rcut_vdw3(ibox) * rcut_vdw3(ibox)
          ENDIF

       ENDIF
    
       IF (int_charge_style(ibox) /= charge_none) THEN
          rcut_coulsq(ibox) = rcut_coul(ibox)*rcut_coul(ibox)
       ENDIF
    
       IF ( (int_vdw_style(ibox) /= vdw_none) .AND. (int_charge_style(ibox) /= charge_none)) THEN
          rcut_max(ibox) = MAX(rcut_vdw(ibox),rcut_coul(ibox))
       ELSE IF ( int_vdw_style(ibox) /= vdw_none) THEN
          rcut_max(ibox) = rcut_vdw(ibox)
       ELSE
          rcut_max(ibox) = rcut_coul(ibox)
       END IF

    END DO

    rcut_lowsq = rcut_low * rcut_low

    ! initialize the pair_nrg array
    
    IF (l_pair_nrg) THEN
       pair_nrg_vdw(:,:) = 0.0_DP
       pair_nrg_qq(:,:) = 0.0_DP
    END IF

    ! Add other stuff like molecule mass, LRC, and anything else needed.
    ! ALLOCATE memory for the Ewald stuff

    ALLOCATE(energy(nbr_boxes),virial(nbr_boxes))

  END SUBROUTINE precalculate




SUBROUTINE generate_reaf_table(ibox)

USE Run_variables
USE Energy_Routines

IMPLICIT NONE

integer, intent(in) :: ibox
integer             :: i, j, k, bin
real                :: cut_reaf
real                :: r_rc, p_rc, p_rc_prime, du_dr, rij, rijsq
cut_reaf = rcut_coul(ibox)

DO i = 1, int(cut_reaf/reaf_width) + 1

rij = i*reaf_width
rijsq = rij*rij
r_rc = rij/rcut_coul(ibox)

p_rc = (1.0 - r_rc)
p_rc_prime = -4.0*p_rc*p_rc*p_rc*(1.0+1.6*r_rc+0.4*r_rc*r_rc)
p_rc_prime = p_rc_prime + p_rc*p_rc*p_rc*p_rc*(1.6+0.8*r_rc)
p_rc = p_rc*p_rc*p_rc*p_rc*(1.0 + 1.6_DP*r_rc + 0.4_DP*r_rc*r_rc)
dU_dr = (-1.0)/rijsq*p_rc + 1.0/sqrt(rijsq*rcut_coulsq(ibox))*p_rc_prime
du_dr = du_dr*(-1.0)*charge_factor


reaf(i,ibox) = dU_dr

END DO

END SUBROUTINE generate_reaf_table

SUBROUTINE Initialize_MP

USE Type_Definitions
USE Run_variables
USE Energy_Routines

!LOGICAL, INTENT(OUT), optional :: mpm_logical

 shell_mpm = .FALSE.

  DO i = 1, nbr_boxes
     DO is = 1, nspecies
      DO im = 1, nmolecules(is)
        DO ia = 1, natoms(is)
            atom_list(ia,im,is)%drude_type = .FALSE.
            if(nonbond_list(ia,is)%element == 'G' ) then
            atom_list(ia,im,is)%drude_type = .TRUE.
            shell_mpm = .TRUE.
            end if
        END DO
      END DO
    END DO
  END DO

  ! mpm force and torque

  IF(shell_mpm) THEN

  ALLOCATE( Mxx(nspecies, MAXVAL(nmolecules)), stat = AllocateStatus)
  ALLOCATE( Myy(nspecies, MAXVAL(nmolecules)), stat = AllocateStatus)
  ALLOCATE( Mzz(nspecies, MAXVAL(nmolecules)), stat = AllocateStatus)

  ALLOCATE( Fxx(nspecies, MAXVAL(nmolecules), MAXVAL(natoms)) , stat = AllocateStatus)
  ALLOCATE( Fyy(nspecies, MAXVAL(nmolecules), MAXVAL(natoms)) , stat = AllocateStatus)
  ALLOCATE( Fzz(nspecies, MAXVAL(nmolecules), MAXVAL(natoms)) , stat = AllocateStatus)

  ALLOCATE(drudeFx(nspecies,maxval(nmolecules),maxval(natoms)))
  ALLOCATE(drudeFy(nspecies,maxval(nmolecules),maxval(natoms)))
  ALLOCATE(drudeFz(nspecies,maxval(nmolecules),maxval(natoms)))
  ALLOCATE(dmpmx(nspecies,SUM(nmolecules)), dmpmy(nspecies,sum(nmolecules)), dmpmz(nspecies,sum(nmolecules)) )
  ALLOCATE(thetax(nspecies,SUM(nmolecules)), thetay(nspecies,sum(nmolecules)), thetaz(nspecies,sum(nmolecules)))
 
  drudeFx = 0.0_DP
  drudeFy = 0.0_DP
  drudeFz = 0.0_DP

  Fxx = 0.0; Fyy = 0.0 ; Fzz = 0.0

  Mxx = 0.0_DP; Myy = 0.0_DP; Mzz = 0.0_DP
 
  dmpmx = 0.0_DP; dmpmy = 0.0_DP; dmpmz = 0.0_DP
  thetax = 0.0_DP; thetay = 0.0_DP; thetaz = 0.0_DP

  reaf_width = 0.01

  ALLOCATE(reaf(int(maxval(rcut_coul)/reaf_width)+1,nbr_boxes))
  
  reaf = 0.0_DP
  
  ! Generate table for reaction field force
  DO ibox = 1, nbr_boxes
      CALL generate_reaf_table(ibox)
  END DO

  mpm_logical = .TRUE.

END IF

END SUBROUTINE Initialize_MP
