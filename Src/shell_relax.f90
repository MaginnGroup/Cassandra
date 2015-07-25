SUBROUTINE shell_relax (this_box, conv_drude)

USE Type_Definitions
USE Run_Variables
USE Random_Generators
USE Simulation_Properties
USE Energy_Routines
!USE Pair_Nrg_Routines

IMPLICIT NONE

INTEGER, INTENT(IN) :: this_box
LOGICAL, INTENT(OUT) :: conv_drude

INTEGER :: is, im, ia, alive
REAL(DP) :: deltaL, xdold, ydold, zdold
INTEGER :: iter_drude, this_locate
REAL(DP) :: lambda_it 


lambda_it =  1.0
iter_drude = 0
conv_drude = .FALSE.

DO WHILE ( (.NOT. conv_drude) .and. iter_drude  .LE. 100)   !Perform 100 iterations

iter_drude  = iter_drude  + 1

CALL Compute_drude_cos_sin_sum(this_box)

DO is = 1, nspecies
  DO ia = 1, natoms(is)        
    IF(atom_list(ia,1,is)%drude_type) THEN   
        !$OMP PARALLEL DO DEFAULT(SHARED) &   
        !$OMP PRIVATE(im, this_locate) &
        !$OMP SCHEDULE(STATIC) 
        DO im = 1, nmolecules(is)                   
          this_locate = locate(im,is)                    
          IF( .NOT. molecule_list(this_locate,is)%live) CYCLE                   
          IF( molecule_list(this_locate,is)%which_box /= this_box ) CYCLE  
          CALL Compute_molecule_drude_force(is,this_locate,ia,this_box, drudeFx(is,this_locate,ia), drudeFy(is,this_locate,ia), drudeFz(is,this_locate,ia) )                     
        END DO
        !$OMP END PARALLEL DO
     END IF  
  END DO
END DO

!if(iter_drude > 20) lambda_it = 0.1 !if the EM cannot converge in 20 steps then use a smaller step size. This is usually due to close intermolecular distance or soft harmonic spring.


conv_drude = .TRUE.

DO is = 1, nspecies
   DO ia = 1, natoms(is)  
      IF(atom_list(ia,1,is)%drude_type) THEN
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP PRIVATE(im, this_locate,deltaL) &
        !$OMP PRIVATE(xdold, ydold, zdold) &
        !$OMP SCHEDULE(STATIC) &
        !$OMP REDUCTION(.AND.:conv_drude)
        DO im = 1, nmolecules(is)  !nmols(is,this_box) !nmolecules(is)
           this_locate = locate(im,is)
           IF( .NOT. molecule_list(this_locate,is)%live) CYCLE          
           IF( molecule_list(this_locate,is)%which_box /= this_box ) CYCLE
           
           deltaL = 0.0    
           xdold = atom_list(ia,this_locate,is)%rxp 
           ydold = atom_list(ia,this_locate,is)%ryp  
           zdold = atom_list(ia,this_locate,is)%rzp      
           
           atom_list(ia,this_locate,is)%rxp = xdold + drudeFx(is,this_locate,ia)/nonbond_list(ia,is)%pol_alpha*lambda_it
           atom_list(ia,this_locate,is)%ryp = ydold + drudeFy(is,this_locate,ia)/nonbond_list(ia,is)%pol_alpha*lambda_it
           atom_list(ia,this_locate,is)%rzp = zdold + drudeFz(is,this_locate,ia)/nonbond_list(ia,is)%pol_alpha*lambda_it
          
           deltaL = (atom_list(ia,this_locate,is)%rxp - xdold)*(atom_list(ia,this_locate,is)%rxp - xdold) +(atom_list(ia,this_locate,is)%ryp - ydold)*(atom_list(ia,this_locate,is)%ryp - ydold) +  (atom_list(ia,this_locate,is)%rzp - zdold)*(atom_list(ia,this_locate,is)%rzp - zdold)

           if(sqrt(deltaL) > 1d-3 .or. maxval(drudeFx) > 50 .or. maxval(drudeFy) > 50.0 .or. maxval(drudeFz) > 50.0 ) then !convergence check
              conv_drude = .FALSE. 
           end if        
        END DO
        !$OMP END PARALLEL DO
     END IF
  END DO
END DO

END DO

IF(.NOT. conv_drude) THEN
check_err = check_err + 1  !accumulate unconvergence counter
print *, 'EM did not converge in 100 iterations'
END IF

END SUBROUTINE shell_relax







     
