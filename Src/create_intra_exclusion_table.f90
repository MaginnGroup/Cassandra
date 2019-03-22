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

SUBROUTINE Create_Intra_Exclusion_Table
  
  !-----------------------------------------------------------------------------
  ! This routine takes the input vdw and charge scaling rules and computes
  ! a arrays of dimensions (natoms,natoms,nspecies) that then get multiplied
  ! by the potential interactions for those atom pairs when performing 
  ! intramolecular potential energy calculations.
  ! Net result is creation of two arrays:
  ! vdw_intra_scale(ii,jj,is)
  ! charge_intra_scale(ii,jj,is)
  ! where ii and jj refer to atom number is species is. To compute 
  ! intramolecular NB energy, we loop over all pairs and multiply 
  ! by these scaling factors. 
  !
  ! Called by
  !
  !  gcmc_control
  !  gemc_control
  !  nptmc_control
  !  nvtmc_control
  !  fragment_control
  !
  ! Revision history
  !
  !   12/10/13 : Beta Release
  !-----------------------------------------------------------------------------  
  USE Global_Variables
  USE Type_Definitions
  USE File_Names
  USE IO_Utilities

  IMPLICIT NONE

  INTEGER :: is,ii,jj,kk
!-----------------------------------------------------------------------------
  ALLOCATE(vdw_intra_scale(MAXVAL(natoms),MAXVAL(natoms),nspecies), Stat=AllocateStatus)
  ALLOCATE(charge_intra_scale(MAXVAL(natoms),MAXVAL(natoms),nspecies), Stat=AllocateStatus)

  IF (AllocateStatus .NE. 0) THEN
     err_msg = ''
     err_msg(1) = ' ERROR: Not enough memory for scaling tables '
     CALL Clean_Abort(err_msg,'create_intra_exclusion_table')
  END IF
  
  DO is=1,nspecies

     DO ii=1,natoms(is)

        DO jj = 1,natoms(is)

           ! Self interactions are automatically excluded. 
           IF (ii == jj) THEN
              vdw_intra_scale(ii,jj,is) = 0.0_DP
              charge_intra_scale(ii,jj,is) = 0.0_DP
           ELSE
              !Set all other interactions to default "1_N" scaling (normally 1.0) 
              ! and turn selected values off using otherinput scaling values
              vdw_intra_scale(ii,jj,is) = scale_1_N_vdw(is)
              charge_intra_scale(ii,jj,is) = scale_1_N_charge(is)
           ENDIF

        ENDDO

     ENDDO

  ENDDO

  ! Now apply specific 1-2, 1-3 and 1-4 scaling rules. 

  SpeciesLoop: DO is = 1,nspecies

     ! 1-4 scaling via dihedrals. Do this first. Later, if atoms
     ! connected by a dihedral and by an angle (for example) are
     ! encountered, then the priority scaling overwrites (i.e. 1-2
     ! has priority over 1-3, which has priority over 1-4).
     ! For example, in the compound below 1-2-3-4 is a 1-4 interaction 
     ! between 1 and 4,but via 1-5-4 it is a 1-3 interaction. We thus
     ! use the 1-3 interaction value. 
     ! 
     !                       1
     !                   /       \           7   10
     !                  /         \          |    |
     !                 2           5----6----9---12
     !                  \         /          |    |
     !                   3------ 4           8   11
     !

     DO kk=1,ndihedrals(is)
        ii = dihedral_list(kk,is)%atom1
        jj = dihedral_list(kk,is)%atom4

        vdw_intra_scale(ii,jj,is) = scale_1_4_vdw(is)
        vdw_intra_scale(jj,ii,is) = scale_1_4_vdw(is)
        charge_intra_scale(ii,jj,is) = scale_1_4_charge(is)
        charge_intra_scale(jj,ii,is) = scale_1_4_charge(is)
        
     ENDDO


     ! 1-3 scaling via angles
     DO kk = 1,nangles(is)

        ii = angle_list(kk,is)%atom1
        jj = angle_list(kk,is)%atom3

        vdw_intra_scale(ii,jj,is) = scale_1_3_vdw(is)
        vdw_intra_scale(jj,ii,is) = scale_1_3_vdw(is)
        charge_intra_scale(ii,jj,is) = scale_1_3_charge(is)
        charge_intra_scale(jj,ii,is) = scale_1_3_charge(is)

     ENDDO

     ! 1-2 scaling via bonds     
     DO kk=1,nbonds(is)

        ii = bond_list(kk,is)%atom1
        jj = bond_list(kk,is)%atom2

        vdw_intra_scale(ii,jj,is) = scale_1_2_vdw(is)
        vdw_intra_scale(jj,ii,is) = scale_1_2_vdw(is)
        charge_intra_scale(ii,jj,is) = scale_1_2_charge(is)
        charge_intra_scale(jj,ii,is) = scale_1_2_charge(is)

     ENDDO
     
  ENDDO SpeciesLoop

  ! report info to log

  IF (verbose_log) THEN
    WRITE(logunit,*)
    WRITE(logunit,*) 'Creating exclusion table'
    WRITE(logunit,'(X,A79)') '-------------------------------------------------------------------------------'

    IF (int_charge_style(1) == charge_none) THEN
      WRITE(logunit,'(x,4(A10))') 'species','atom1','atom2','vdw_scale'
    ELSE
      WRITE(logunit,'(x,5(A10))') 'species','atom1','atom2','vdw_scale','qq_scale'
    END IF

    DO is=1,nspecies
       IF (natoms(is) < 100) THEN
         DO ii=1,natoms(is)
            DO jj = ii+1,natoms(is)
               IF (int_charge_style(1) == charge_none) THEN
                 WRITE(logunit,'(X,3(I10),F10.3)') is,ii,jj,vdw_intra_scale(ii,jj,is)
               ELSE
                 WRITE(logunit,'(X,3(I10),2(F10.3))') is,ii,jj,vdw_intra_scale(ii,jj,is),&
                      charge_intra_scale(ii,jj,is)
               END IF
            ENDDO
         ENDDO
       ELSE
         WRITE(logunit,'(X,A)') 'Species ' // TRIM(Int_To_String(is)) // ' has more than 10,000 interactions'
       END IF
    ENDDO

    WRITE(logunit,'(X,A)') '-------------------------------------------------------------------------------'
    WRITE(logunit,'(A)') '********************************************************************************'
  END IF


END SUBROUTINE Create_Intra_Exclusion_Table
