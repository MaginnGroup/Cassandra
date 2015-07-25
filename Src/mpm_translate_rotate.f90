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




!********************************************************************************
! This subroutine performs translate/rotate move for all the particles in this_box. 
! This is a force biased move.
SUBROUTINE Translate_MP(this_box,mc_step)
!********************************************************************************

  USE Type_Definitions
  USE Run_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE Pair_Nrg_Routines

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: this_box
  INTEGER  :: is, im, i, j, nmolecules_species, alive, total_mols, ibox, nmols_box, this_species
  INTEGER  :: mc_step, ia

  REAL(DP) :: delta_e, p_acc 
  REAL(DP) :: E_intra, E_vdw, E_qq
  REAL(DP) :: E_vdw_move, E_qq_move
  REAL(DP) :: E_intra_qq, E_intra_vdw
  REAL(DP) :: E_intra_qq_move, E_intra_vdw_move
  REAL(DP) :: E_intra_qq_alive, E_intra_vdw_alive, E_intra_qq_move_alive, E_intra_vdw_move_alive
  REAL(DP) :: E_bond, E_bond_move, E_bond_alive
  REAL(DP) :: E_reciprocal_move
  REAL(DP) :: Energy_old_reciprocal
  REAL(DP) :: E_vdw_move_alive,E_qq_move_alive, E_vdw_alive, E_qq_alive

  REAL(DP) :: success_ratio
  REAL(DP) :: rcut_small

  LOGICAL :: inter_overlap, overlap, total_overlap, intra_overlap
  LOGICAL  :: energy_cal
  LOGICAL :: conv_drude

  REAL(DP) :: attempt_p
  
  REAL(DP) :: gamma_x_a , gamma_y_a , gamma_z_a 
  REAL(DP) :: S_new_old, S_old_new, S_is
  REAL(kind=8) :: g_random
  REAL(kind=8) :: dev_b, mean_a, zero, one
  REAL(kind=8) :: dispp_max

!  REAL(dp) :: g_random
!  REAL(dp) :: dev_b, mean_a
!  REAL(dp) :: dispp_max

  E_vdw_move = 0.0_DP
  E_qq_move = 0.0_DP
  E_vdw = 0.0_DP
  E_qq = 0.0_DP
  E_reciprocal_move = 0.0_DP
  attempt_p = 1.0_DP

  mean_a = 0.D0
  one = 1.D0

  total_overlap = .FALSE.

  conv_Drude = .TRUE.  



  total_mols = SUM(nmols(:,this_box))
  
  ! If there are no molecules in the box then return to driver
  IF (total_mols == 0) RETURN


  rcut_small = MIN(rcut_vdw(this_box),rcut_coul(this_box))  
  IF (rcut_small < 1.0) rcut_small = rcut_vdw(this_box)  

  ! update the trial counter for this molecule 
  
  tot_trials(this_box) = tot_trials(this_box) + 1 
  DO is = 1, nspecies   
  ntrials(is,this_box)%displacement = ntrials(is,this_box)%displacement + 1
  END DO

  E_vdw = 0.0; E_qq = 0.0; inter_overlap = .FALSE.         

  energy_cal =.FALSE.
 
   DO is = 1, nspecies
      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP SCHEDULE(STATIC) &
      !$OMP PRIVATE(im, alive)
      DO im = 1, nmolecules(is)    !nmols(is,this_box)           
        alive = locate(im,is)
        if(molecule_list(alive,is)%which_box /= this_box) cycle
        if(.NOT. molecule_list(alive,is)%live) cycle
        CALL Save_Old_Cartesian_Coordinates(alive,is)
        call Compute_Molecule_Nonbond_Inter_Force(alive, is, Fxx(is, alive, :), Fyy(is, alive, :), Fzz(is, alive, :), &
               & E_vdw_alive, E_qq_alive, inter_overlap, energy_cal)

        IF (inter_overlap) THEN
        WRITE(*,*) 'Error in the old configuration -- > overlap', im ,is
        STOP
        END IF
     END DO 
       !$OMP END PARALLEL DO
   END DO


   E_vdw = energy(this_box)%inter_vdw
   E_qq = energy(this_box)%inter_q
   E_intra_qq = energy(this_box)%intra_q
   E_intra_vdw = energy(this_box)%intra_vdw
   E_bond = energy(this_box)%bond
   E_intra = energy(this_box)%intra

  ! Compute the displacement vector for all partciles
  ! dispp_max is an adjustable parameter for multiparticle move

   S_new_old = 0.0; S_old_new = 0.0

   LOOP1:DO is = 1, nspecies
      S_is = 0.0
      dispp_max = max_disp(is,this_box) 
      dev_b = dsqrt(dispp_max*2.0_DP)
      if(dispp_max*1000.0 < 1d-7) cycle 

      LOOP2:DO im = 1, nmolecules(is) !nmols(is,this_box)     
        alive = locate(im,is)
        if(molecule_list(alive,is)%which_box /= this_box) cycle
        if(.NOT. molecule_list(alive,is)%live) cycle      
        gamma_x_a = g_random(mean_a, dev_b)  
        gamma_y_a = g_random(mean_a, dev_b) 
        gamma_z_a = g_random(mean_a, dev_b)    
        dmpmx(is,alive) = beta(this_box)*dispp_max*sum(Fxx(is, alive, :)) + gamma_x_a
        dmpmy(is,alive) = beta(this_box)*dispp_max*sum(Fyy(is, alive, :)) + gamma_y_a
        dmpmz(is,alive) = beta(this_box)*dispp_max*sum(Fzz(is, alive, :)) + gamma_z_a
        atom_list(:,alive,is)%rxp = atom_list(:,alive,is)%rxp + dmpmx(is,alive)
        atom_list(:,alive,is)%ryp = atom_list(:,alive,is)%ryp + dmpmy(is,alive)
        atom_list(:,alive,is)%rzp = atom_list(:,alive,is)%rzp + dmpmz(is,alive)
        molecule_list(alive,is)%xcom = molecule_list(alive,is)%xcom + dmpmx(is,alive)
        molecule_list(alive,is)%ycom = molecule_list(alive,is)%ycom + dmpmy(is,alive)
        molecule_list(alive,is)%zcom = molecule_list(alive,is)%zcom + dmpmz(is,alive)
        S_is = S_is + (beta(this_box)*dispp_max*sum(Fxx(is, alive, :)) - dmpmx(is,alive))*&
                                (beta(this_box)*dispp_max*sum(Fxx(is, alive, :)) - dmpmx(is,alive))
        S_is = S_is + (beta(this_box)*dispp_max*sum(Fyy(is, alive, :)) - dmpmy(is,alive))*&
                                (beta(this_box)*dispp_max*sum(Fyy(is, alive, :)) - dmpmy(is,alive))
        S_is = S_is + (beta(this_box)*dispp_max*sum(Fzz(is, alive, :)) - dmpmz(is,alive))*&
                                (beta(this_box)*dispp_max*sum(Fzz(is, alive, : )) - dmpmz(is,alive))
       
       END DO LOOP2
      S_old_new = S_old_new + S_is/dispp_max
     END DO LOOP1
   S_old_new = S_old_new /4.D0*(-1.D0)

 
  E_vdw_move = 0.0; E_qq_move = 0.0; inter_overlap = .FALSE.; total_overlap = .FALSE.
   
   LOOP4: DO is = 1, nspecies
      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP SCHEDULE(STATIC) &
      !$OMP PRIVATE(im, alive, inter_overlap) &
      !$OMP REDUCTION(.OR.:total_overlap)
      LOOP3: DO im = 1, nmolecules(is)  !nmols(is,this_box)        
        alive = locate(im,is)
        if(molecule_list(alive,is)%which_box /= this_box) cycle
        if(.NOT. molecule_list(alive,is)%live) cycle         
        CALL Fold_Molecule(alive,is,this_box)    
        CALL compute_molecule_inter_overlap(alive,is,this_box,inter_overlap)
        IF (inter_overlap) THEN
	       total_overlap = .TRUE.
        END IF
       END DO LOOP3       
     !$OMP END PARALLEL DO
   END DO LOOP4


  IF (total_overlap ) THEN
     DO is = 1, nspecies         
      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP SCHEDULE(STATIC) &
      !$OMP PRIVATE(im, alive) 
       DO im = 1, nmolecules(is)       
         alive = locate(im,is)
         if(molecule_list(alive,is)%which_box /= this_box) cycle
         if(.NOT. molecule_list(alive,is)%live) cycle 
         CALL Revert_Old_Cartesian_Coordinates(alive,is)
       END DO
       !$OMP END PARALLEL DO
     END DO
  END IF


  IF (.NOT. total_overlap) THEN 

      conv_drude = .TRUE.
      IF(shell_mpm) CALL shell_relax(this_box, conv_drude) ! perform electrostatic energy minimization

      E_vdw_move = 0.0; E_qq_move = 0.0
      E_intra_vdw_move = 0.0; E_intra_qq_move = 0.0
      E_bond_move = 0.0      

      energy_cal = .TRUE.      
      IF (conv_drude) THEN
      DO is = 1, nspecies       
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP SCHEDULE(STATIC) &
        !$OMP PRIVATE(im, alive) &
        !$OMP PRIVATE(E_vdw_move_alive, E_qq_move_alive,inter_overlap) &
        !$OMP PRIVATE(E_intra_vdw_move_alive, E_intra_qq_move_alive,E_bond_alive,intra_overlap) &
        !$OMP REDUCTION(+:E_vdw_move,E_qq_move, E_intra_qq_move, E_intra_vdw_move, E_bond_move)

        DO im = 1, nmolecules(is) !nmols(is,this_box) 
          alive = locate(im,is)
          if(molecule_list(alive,is)%which_box /= this_box) cycle
          if(.NOT. molecule_list(alive,is)%live) cycle
          CALL Compute_Molecule_Nonbond_Inter_Force(alive, is, Fxx(is, alive, :), Fyy(is, alive, :), Fzz(is, alive, :), &
                & E_vdw_move_alive, E_qq_move_alive, inter_overlap, energy_cal) 
          CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw_move_alive,E_intra_qq_move_alive,intra_overlap)
          CALL Compute_Molecule_Bond_Energy(alive,is,E_bond_alive)
          E_vdw_move = E_vdw_move_alive + E_vdw_move
          E_qq_move = E_qq_move_alive + E_qq_move
          E_intra_qq_move = E_intra_qq_move_alive + E_intra_qq_move 
          E_intra_vdw_move = E_intra_vdw_move_alive + E_intra_vdw_move
          E_bond_move = E_bond_move + E_bond_alive
         END DO

        !$OMP END PARALLEL DO

      END DO

     LOOP22: DO is = 1, nspecies
       S_is = 0.0
       dispp_max = max_disp(is,this_box)
       if(dispp_max*1000.0 < 1d-7) cycle
       !$OMP PARALLEL DO DEFAULT(SHARED) &
       !$OMP SCHEDULE(STATIC) &
       !$OMP PRIVATE(im, alive) &
       !$OMP REDUCTION(+:S_is)
       LOOP23: DO im = 1, nmolecules(is) !nmols(is,this_box)    
         alive = locate(im,is)
         if(molecule_list(alive,is)%which_box /= this_box) cycle
         if(.NOT. molecule_list(alive,is)%live) cycle        
         S_is = S_is + (beta(this_box)*dispp_max*sum(Fxx(is, alive, :)) + dmpmx(is,alive))*&
                                (beta(this_box)*dispp_max*sum(Fxx(is, alive, :)) + dmpmx(is,alive))
         S_is = S_is + (beta(this_box)*dispp_max*sum(Fyy(is, alive, :)) + dmpmy(is,alive))*&
                                (beta(this_box)*dispp_max*sum(Fyy(is, alive, :)) + dmpmy(is,alive))
         S_is = S_is + (beta(this_box)*dispp_max*sum(Fzz(is, alive, :)) + dmpmz(is,alive))*&
                                (beta(this_box)*dispp_max*sum(Fzz(is, alive, :)) + dmpmz(is,alive))       
       END DO LOOP23
       !$OMP END PARALLEL DO
       S_new_old = S_is/dispp_max + S_new_old
     END DO LOOP22     
     S_new_old = S_new_old / 4.0 * (-1.0)
    
     END IF

  END IF

  IF ((int_charge_sum_style(this_box) == charge_ewald .or. int_charge_sum_style(this_box) == charge_gaussian)  .and. .NOT. total_overlap .and. conv_drude) THEN
  Energy_old_reciprocal = energy(this_box)%ewald_reciprocal 
  CALL Compute_System_Ewald_Reciprocal_Energy2(this_box)
  END IF

  delta_e = ( E_vdw_move/2.0 - E_vdw ) + ( E_qq_move/2.0 - E_qq ) 
  delta_e = delta_e + E_bond_move - E_bond + E_intra_vdw_move + E_intra_qq_move - E_intra_qq - E_intra_vdw
  delta_e = delta_e + energy(this_box)%ewald_reciprocal - Energy_old_reciprocal

  p_acc = MIN( 1.0_DP, ( DEXP(-beta(this_box) * delta_e + S_new_old - S_old_new ) ) )

  IF ( rranf() <= p_acc .and. .NOT. total_overlap .and. conv_drude ) THEN

   energy(this_box)%inter_vdw = E_vdw_move/2.0
   energy(this_box)%inter_q   = E_qq_move/2.0
   energy(this_box)%intra_q   = E_intra_qq_move
   energy(this_box)%intra_vdw = E_intra_vdw_move
   energy(this_box)%bond      = E_bond_move 
   energy(this_box)%intra     = E_bond_move 
   energy(this_box)%total     = energy(this_box)%total + delta_e
   nsuccess(:,this_box)%displacement = nsuccess(:,this_box)%displacement + 1
   nequil_success(:,this_box)%displacement = nequil_success(:,this_box)%displacement + 1
   
  ELSE
  
   IF (.NOT. total_overlap) THEN
    DO is = 1, nspecies
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP SCHEDULE(STATIC) &
        !$OMP PRIVATE(im, alive) 
        DO im = 1, nmolecules(is) !nmols(is,this_box)
          alive = locate(im,is)
          if(molecule_list(alive,is)%which_box /= this_box) cycle
          if(.NOT. molecule_list(alive,is)%live) cycle
          
          CALL Revert_Old_Cartesian_Coordinates(alive,is)
        END DO 
        !$OMP END PARALLEL DO
     END DO
    
     IF (int_charge_sum_style(this_box) == charge_ewald .or. int_charge_sum_style(this_box) == charge_gaussian  .and. conv_drude) THEN
     Energy(this_box)%ewald_reciprocal = Energy_old_reciprocal
     END IF

   END IF
  
 END IF
  


 DO is = 1, nspecies
     IF ( MOD(ntrials(is,this_box)%displacement,nupdate) == 0 ) THEN
     IF ( int_run_style == run_equil ) THEN 
        success_ratio = REAL(nequil_success(is,this_box)%displacement,DP)/REAL(nupdate,DP)
     ELSE
        success_ratio = REAL(nsuccess(is,this_box)%displacement,DP)/REAL(ntrials(is,this_box)%displacement,DP)
     END IF
     WRITE(logunit,*)
     WRITE(logunit,'(A,2X,I3,2X,A,I2,2x,A,2X,F8.5)')' Multiparticle Translation ratio for species ', is , 'for box', this_box, 'is: ', success_ratio
     IF ( int_run_style == run_equil ) THEN
        ! check if the acceptace is close to 0.33
         nequil_success(is,this_box)%displacement = 0
         IF  ( success_ratio < 0.001 ) THEN
             max_disp(is,this_box) = 0.1_DP*max_disp(is,this_box)
         ELSE
             max_disp(is,this_box) = MIN(rcut_small, 3.0_DP*success_ratio*max_disp(is,this_box))
         END IF
         IF(max_disp(is,this_box) .gt. 25.0) max_disp(is,this_box) = 25.0
         WRITE(logunit,'(A,1X,I1,1X,A,1X,I1)') 'Parameter A for multiparticle translation of species', is,' updated in box', this_box
         WRITE(logunit,'(A,2X,F10.5)') 'is', max_disp(is,this_box)
         WRITE(logunit,*)  
     END IF
  END IF
END DO


END SUBROUTINE Translate_MP
     

     

SUBROUTINE Rotate_MP(this_box, mc_step)


  USE Type_Definitions
  USE Run_Variables
  USE Random_Generators
  USE Simulation_Properties
  USE Energy_Routines
  USE File_Names
  USE Pair_Nrg_Routines
  IMPLICIT NONE

!  !$ include 'omp_lib.h'
  INTEGER, INTENT(IN) :: this_box
  INTEGER  :: mc_step
  INTEGER  :: is, im, ia, nmolecules_species, alive, i, nmols_box, total_mols, this_species
  INTEGER  :: ibox

  REAL(DP) :: E_intra, E_vdw, E_qq
  REAL(DP) :: E_intra_qq, E_intra_vdw
  REAL(DP) :: E_vdw_alive, E_qq_alive, E_vdw_move_alive, E_qq_move_alive, Energy_old_reciprocal
  REAL(DP) :: E_intra_qq_alive, E_intra_vdw_alive, E_intra_vdw_move_alive, E_intra_qq_move_alive
  REAL(DP) :: E_intra_qq_move, E_intra_vdw_move, E_vdw_move, E_qq_move
  REAL(DP) :: E_bond, E_bond_move, E_bond_alive, E_reciprocal_move
  REAL(DP) :: delta_e, p_acc, success_ratio
  REAL(DP) :: attempt_p
  LOGICAL  :: inter_overlap
  LOGICAL  :: overlap,total_overlap, energy_cal,intra_overlap
  LOGICAL  :: conv_drude

  REAL(DP) :: S_new_old , S_old_new, S_is
  REAL(DP) :: gamma_x_a, gamma_y_a, gamma_z_a
  REAL(kind=8) :: g_random
  REAL(kind=8) :: mean_a, dev_b
  REAL(kind=8):: rott_max
  
  
  attempt_p = 1.0_DP

  inter_overlap = .FALSE.
  E_vdw_move = 0.0_DP
  E_vdw = 0.0_DP
  E_qq_move = 0.0_DP
  E_qq = 0.0_DP
  E_reciprocal_move = 0.0_DP
  mean_a = 0.D0

  total_overlap = .FALSE.
  conv_Drude = .TRUE.  

  total_mols = SUM(nmols(:,this_box))  
 ! If there are no molecules in the box then return

  IF (total_mols == 0) RETURN  

  tot_trials(this_box) = tot_trials(this_box) + 1
  DO is = 1, nspecies
   ntrials(is,this_box)%rotation = ntrials(is,this_box)%rotation + 1
  END DO

 ! Evaluate force and torque on each molecule

   E_vdw = 0.0; E_qq = 0.0; inter_overlap = .FALSE.  
   
   energy_cal = .FALSE.
  
   LOOP31: DO is = 1, nspecies
     !$OMP PARALLEL DO DEFAULT(SHARED) &
     !$OMP SCHEDULE(STATIC) &
     !$OMP PRIVATE(im, alive) 
     LOOP32:  DO im = 1, nmolecules(is) ! nmols(is,this_box)  
      alive = locate(im,is)
      if(molecule_list(alive,is)%which_box /= this_box) cycle
      if(.NOT. molecule_list(alive,is)%live) cycle
      
       CALL Save_Old_Cartesian_Coordinates(alive,is)
       CALL Compute_Molecule_Nonbond_Inter_Force(alive, is, Fxx(is, alive, :), Fyy(is, alive, :), Fzz(is, alive, :), &
           & E_vdw_alive, E_qq_alive, inter_overlap, energy_cal)
       IF (inter_overlap) THEN
          WRITE(*,*) 'Error in the old configuration "in rotation" -- > overlap'
          STOP
       END IF    
       CALL Compute_Molecule_Torque(alive, is, this_box, Fxx(is, alive,:), Fyy(is, alive,:), Fzz(is, alive,:), Mxx(is, alive), Myy(is, alive), Mzz(is,alive) )          
     END DO LOOP32     
     !$OMP END PARALLEL DO
   END DO LOOP31


    E_vdw = energy(this_box)%inter_vdw
    E_qq = energy(this_box)%inter_q
    E_intra_qq = energy(this_box)%intra_q
    E_intra_vdw = energy(this_box)%intra_vdw
    E_bond = energy(this_box)%bond
    E_intra = energy(this_box)%intra

    S_new_old = 0.0; S_old_new = 0.0

    LOOP33: DO is = 1, nspecies
     S_is = 0.0_DP
     rott_max = max_rot(is,this_box)
     dev_b = dsqrt(rott_max*2.0)
!     if(is .gt. 2) cycle
     LOOP34: DO im = 1, nmolecules(is)  !nmols(is,this_box)
        alive = locate(im,is)
        if(molecule_list(alive,is)%which_box /= this_box) cycle
        if(.NOT. molecule_list(alive,is)%live) cycle        
        gamma_x_a = g_random(mean_a, dev_b) 
        gamma_y_a = g_random(mean_a, dev_b) 
        gamma_z_a = g_random(mean_a, dev_b)
        thetax(is, alive) =  beta(this_box)*rott_max*Mxx(is, alive) + gamma_x_a
        thetay(is, alive) =  beta(this_box)*rott_max*Myy(is, alive) + gamma_y_a
        thetaz(is, alive) =  beta(this_box)*rott_max*Mzz(is, alive) + gamma_z_a
        S_is = S_is + (beta(this_box)*rott_max*Mxx(is,alive) - thetax(is,alive))*&
                                 (beta(this_box)*rott_max*Mxx(is,alive) - thetax(is,alive)) 
        S_is = S_is + (beta(this_box)*rott_max*Myy(is,alive) - thetay(is,alive))*&
                                (beta(this_box)*rott_max*Myy(is,alive) - thetay(is,alive)) 
        S_is = S_is + (beta(this_box)*rott_max*Mzz(is,alive) - thetaz(is,alive))*&
                                (beta(this_box)*rott_max*Mzz(is,alive) - thetaz(is,alive)) 
        CALL Rotate_Molecule_Axis_MPM(alive, is,  thetax(is, alive), thetay(is, alive), thetaz(is, alive)) !Rigid rotate molecule around COM, the rotation angle is biased by the torque on the molecule
      END DO LOOP34
      S_old_new = S_is/rott_max + S_old_new
   END DO LOOP33

   S_old_new = S_old_new /4.D0*(-1.D0)

   ! Check possible overlap, calculate the new E_vdw and E_qq, Caclulate force of new configuration, Fxx, Fyy, Fzz

    inter_overlap = .FALSE.    

    LOOP35: DO is = 1, nspecies
      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP SCHEDULE(STATIC) &
      !$OMP PRIVATE(im, alive, inter_overlap) &
      !$OMP REDUCTION(.OR.:total_overlap)
      LOOP36: DO im = 1, nmolecules(is) !nmols(is,this_box)             
        alive = locate(im,is)
        if(molecule_list(alive,is)%which_box /= this_box) cycle
        if(.NOT. molecule_list(alive,is)%live) cycle        
        CALL Fold_Molecule(alive,is,this_box)
        CALL compute_molecule_inter_overlap(alive,is,this_box,inter_overlap)
        IF (inter_overlap) THEN
            total_overlap = .TRUE.
        END IF
      END DO LOOP36
     !$OMP END PARALLEL DO
   END DO LOOP35
  
  IF (.NOT. total_overlap) THEN
     
      Conv_Drude = .TRUE.
      IF(shell_mpm) CALL shell_relax(this_box, conv_drude)  ! perform electrostatic energy minimization

      E_vdw_move = 0.0; E_qq_move = 0.0      
      E_intra_vdw_move = 0.0; E_intra_qq_move = 0.0
      E_bond_move = 0.0 
      
      IF (conv_drude) THEN

      energy_cal = .TRUE.
      DO is = 1, nspecies
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP SCHEDULE(STATIC) &
        !$OMP PRIVATE(im, alive) &
        !$OMP PRIVATE(E_vdw_move_alive, E_qq_move_alive,inter_overlap) &
        !$OMP PRIVATE(E_intra_vdw_move_alive, E_intra_qq_move_alive,E_bond_alive,intra_overlap) &
        !$OMP REDUCTION(+:E_vdw_move,E_qq_move, E_intra_qq_move, E_intra_vdw_move, E_bond_move)
        DO im = 1, nmolecules(is) 
          alive = locate(im,is)
          if(molecule_list(alive,is)%which_box /= this_box) cycle
          if(.NOT. molecule_list(alive,is)%live) cycle          
          CALL Compute_Molecule_Nonbond_Inter_Force(alive, is,Fxx(is, alive, :), Fyy(is, alive, :), Fzz(is, alive, :), &
               & E_vdw_move_alive, E_qq_move_alive, inter_overlap, energy_cal)
          CALL Compute_Molecule_Torque(alive, is, this_box, Fxx(is, alive,:), Fyy(is, alive,:), Fzz(is, alive,:), Mxx(is, alive), Myy(is, alive), Mzz(is,alive) )
          CALL Compute_Molecule_Nonbond_Intra_Energy(alive,is,E_intra_vdw_move_alive,E_intra_qq_move_alive,intra_overlap)
          CALL Compute_Molecule_Bond_Energy(alive,is,E_bond_alive)
          E_vdw_move = E_vdw_move_alive + E_vdw_move
          E_qq_move = E_qq_move_alive + E_qq_move
          E_intra_qq_move = E_intra_qq_move_alive + E_intra_qq_move 
          E_intra_vdw_move = E_intra_vdw_move_alive + E_intra_vdw_move
          E_bond_move = E_bond_move + E_bond_alive
         END DO
        !$OMP END PARALLEL DO
      END DO
      

      LOOP37: DO is = 1, nspecies
         
         S_is = 0.0
         rott_max = max_rot(is,this_box)
 !        if(is .gt. 2) cycle
         !$OMP PARALLEL DO DEFAULT(SHARED) &
         !$OMP SCHEDULE(STATIC) &
         !$OMP PRIVATE(im, alive) &
         !$OMP REDUCTION(+:S_is)
         LOOP38: DO im = 1, nmolecules(is)  
               alive = locate(im,is)
               if(molecule_list(alive,is)%which_box /= this_box) cycle
               if(.NOT. molecule_list(alive,is)%live) cycle
               S_is = S_is + (beta(this_box)*rott_max*Mxx(is,alive) + thetax(is,alive))*&
                                (beta(this_box)*rott_max*Mxx(is,alive) + thetax(is,alive)) 
               S_is = S_is + (beta(this_box)*rott_max*Myy(is,alive) + thetay(is,alive))*&
                                (beta(this_box)*rott_max*Myy(is,alive) + thetay(is,alive)) 
               S_is = S_is + (beta(this_box)*rott_max*Mzz(is,alive) + thetaz(is,alive))*&
                                (beta(this_box)*rott_max*Mzz(is,alive) + thetaz(is,alive))    
         END DO LOOP38
        !$OMP END PARALLEL DO
              S_new_old = S_is/rott_max + S_new_old
       END DO LOOP37  

       S_new_old = S_new_old / 4.0 * (-1.0)

      IF (int_charge_sum_style(this_box) == charge_ewald .or. int_charge_sum_style(this_box) == charge_gaussian .and. conv_drude) THEN
         Energy_old_reciprocal = energy(this_box)%ewald_reciprocal
         CALL Compute_System_Ewald_Reciprocal_Energy2(this_box)
      END IF

      END IF

    ELSE
      
      DO is = 1, nspecies
         !$OMP PARALLEL DO DEFAULT(SHARED) &
         !$OMP SCHEDULE(STATIC) &
         !$OMP PRIVATE(im, alive) 
         DO im = 1, nmolecules(is)        
          alive = locate(im,is)
          if(molecule_list(alive,is)%which_box /= this_box) cycle
          if(.NOT. molecule_list(alive,is)%live) cycle
          
          CALL Revert_Old_Cartesian_Coordinates(alive,is)    
         
         END DO
         !$OMP END PARALLEL DO
      END DO
      
    END IF   

   delta_e = ( E_vdw_move/2.0 - E_vdw ) + ( E_qq_move/2.0 - E_qq ) 
   delta_e = delta_e + E_bond_move - E_bond + E_intra_vdw_move + E_intra_qq_move - E_intra_qq - E_intra_vdw
   delta_e = delta_e + energy(this_box)%ewald_reciprocal - Energy_old_reciprocal   
   p_acc = MIN( 1.0_DP, attempt_p * ( DEXP(-beta(this_box) * delta_e + S_new_old - S_old_new ) ) )

   IF ( rranf() <= p_acc .and. .NOT. total_overlap .and. conv_drude ) THEN
     energy(this_box)%inter_vdw = E_vdw_move/2.0
     energy(this_box)%inter_q   = E_qq_move/2.0
     energy(this_box)%intra_q   = E_intra_qq_move
     energy(this_box)%intra_vdw = E_intra_vdw_move
     energy(this_box)%bond      = E_bond_move 
     energy(this_box)%intra     = E_bond_move 
     energy(this_box)%total     = energy(this_box)%total + delta_e
     nsuccess(:,this_box)%rotation = nsuccess(:,this_box)%rotation + 1
     nequil_success(:,this_box)%rotation = nequil_success(:,this_box)%rotation + 1  
   ELSE

    IF ( .NOT. total_overlap) THEN

    DO is = 1, nspecies
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP SCHEDULE(STATIC) &
        !$OMP PRIVATE(im, alive) 
        DO im = 1, nmolecules(is)  !nmols(is,this_box)
          alive = locate(im,is)
          if(molecule_list(alive,is)%which_box /= this_box) cycle
          if(.NOT. molecule_list(alive,is)%live) cycle
          CALL Revert_Old_Cartesian_Coordinates(alive,is)
        END DO 
        !$OMP END PARALLEL DO
     END DO

     IF (int_charge_sum_style(this_box) == charge_ewald .or. int_charge_sum_style(this_box) == charge_gaussian .and. conv_drude ) THEN
     Energy(this_box)%ewald_reciprocal = Energy_old_reciprocal
     END IF
    END IF   
  END IF


IF (MOD(ntrials(1,this_box)%rotation, nupdate) == 0) THEN    
     DO is = 1,nspecies
     IF (int_run_style == run_equil) THEN
           success_ratio = REAL(nequil_success(is,this_box)%rotation,DP)/REAL(nupdate,DP)
        ELSE
           success_ratio = REAL(nsuccess(is,this_box)%rotation,DP)/REAL(ntrials(is,this_box)%rotation,DP)
     END IF
     WRITE(logunit,*)
     WRITE(logunit,'(A,2x,I2,2X,A,2X,I2,2X,A,2X,F8.5)') 'Success ratio for rotation of species', is,'for box',this_box, 'is:', success_ratio
     WRITE(logunit,*) 'adjust parameter B', max_rot(is,this_box)
     WRITE(logunit,*)
     IF (int_run_style == run_equil) THEN   
         nequil_success(is,this_box)%rotation = 0
        ! update the width of the species for equilibration run
         IF  ( success_ratio < 0.01 ) THEN
             max_rot(is,this_box) = 0.1_DP*max_rot(is,this_box)
         ELSE
             max_rot(is,this_box)  = 3.0_DP*success_ratio*max_rot(is,this_box)
        END IF    
      END IF     
     ENDDO 

END IF
  
CONTAINS

  SUBROUTINE Rotate_Molecule_Axis_MPM(alive, is, u, v, w)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: alive, is
    REAL(DP), INTENT(IN) :: u, v, w
    INTEGER :: axis, ia
    REAL(DP) :: a, b, c
    REAL(DP) :: theta, L, Lsr
    REAL(DP) :: x1, y1, z1
    REAL(DP) :: x2, y2, z2
    theta = U*u + v*v + w*w
    L     = theta
    Lsr   = sqrt(L)
    theta = sqrt(theta)    
    a = molecule_list(alive,is)%xcom
    b = molecule_list(alive,is)%ycom
    c = molecule_list(alive,is)%zcom
    DO ia = 1, natoms(is)
       x1 = atom_list(ia,alive,is)%rxp
       y1 = atom_list(ia,alive,is)%ryp
       z1 = atom_list(ia,alive,is)%rzp
       x2 = (a*(v*v + w*w) - u*(b*v+c*w-u*x1-v*y1-w*z1))*(1.0 - cos(theta)) + L*x1*cos(theta) &
        & + Lsr*(-c*v+b*w-w*y1+v*z1)*sin(theta)
       y2 = (b*(u*u + w*w) - v*(a*u+c*w-u*x1-v*y1-w*z1))*(1.0 - cos(theta)) + L*y1*cos(theta) &
                & + Lsr*(c*u-a*w+w*x1-u*z1)*sin(theta)
       z2 = (c*(u*u + v*v) - w*(a*u+b*v-u*x1-v*y1-w*z1))*(1.0 - cos(theta)) + L*z1*cos(theta) &
                & + Lsr*(-b*u+a*v-v*x1+u*y1)*sin(theta)
       x2 = x2/L
       y2 = y2/L
       z2 = z2/L       
       atom_list(ia,alive,is)%rxp = x2
       atom_list(ia,alive,is)%ryp = y2
       atom_list(ia,alive,is)%rzp = z2       
   END DO  
  END SUBROUTINE Rotate_Molecule_Axis_MPM  
END SUBROUTINE Rotate_MP



! Generate random number with Gaussian distribution
function g_random(m, s)
USE Type_Definitions
USE Random_Generators
implicit none
real(kind=8)  :: g_random 
real(kind=8)   :: m , s
real(kind=8)   :: x1, x2, w, y1
real(kind=8), save   :: y2
integer(kind=4), save :: use_last = 0

IF(use_last == 1) THEN
  y1 = y2
  use_last = 0
ELSE
  w = 2.0
  DO WHILE (w >= 1.0)
    x1 = 2.0*rranf() - 1.0
    x2 = 2.0*rranf() - 1.0
    w = x1*x1 + x2*x2
    END DO
    w = sqrt(-2.0*log(w)/w)
    y1 = x1*w
    y2 = x2*w
    use_last = 1
END IF
g_random = m + y1*s
return

END function g_random
