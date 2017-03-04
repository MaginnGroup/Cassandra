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

MODULE Angle_Dist_Pick

!*********************************************************************************
! The module contains two routines.
! 1. Angle_Distribution: 
!       This routine computes the probability of observing
!       a particular value of bond angle for the angle in question.
!
!*************Used by************
! nvtmc_control.f90
! nptmc_control.f90
! gcmc_control.f90
! gemc_control.f90
!********************************
! 2. Angle_Pick: 
!        The routine pick an angle for angle distortion move.
!
! Revision history
!
! 12/10/13  : Beta Release
!*********************************************************************************

  USE Type_Definitions
  USE Global_Variables
  USE File_Names


  IMPLICIT NONE

CONTAINS

!*********************************************************************************
  SUBROUTINE Angle_Distribution(this_box)
!*********************************************************************************
    ! This routine calculates the probability of observing a given bond angle
    ! value. It is assumed that the bond angle distributions are decoupled from
    ! the other degrees of freedom. hence, the distribution corresponds to that
    ! observed in ideal gas.
    !*****************************************************************************

    !*****************************************************************************
    ! Since the angle distribution is uniform over cosine of angles, we will
    ! generate a uniform distribution between [-1,1]. The interval is discretized
    ! in Nregions. We will calculate the Boltzmann probability of observing the
    ! mid point of these regions. We will exclue the points that occur with
    ! negligible probability. This is determined if the energy*kbT > max_kBT
    ! Thus out of possible Nregions, we will have calculated probability of
    ! ngood_regions. 
    !***************************************************************************

    IMPLICIT NONE


    INTEGER :: ispecies, iangles, iregions, ngood_regions, i, this_box

    REAL(DP) :: delta_cos, costheta_0, ktheta_0, theta_0, cos_iregions, theta_iregions
    REAL(DP) :: vtheta_iregions, T   

    CHARACTER(20) :: angle_potential

    ! Allocate the memory for the ang_prob list

    IF (.NOT. ALLOCATED(ang_prob)) ALLOCATE(ang_prob(MAXVAL(nangles),nspecies))

    T = temperature(this_box)

    ! Determine the width of a region

    delta_cos = 2.0_DP/nregions

    ! starting cosine value

    costheta_0 = -1.0_DP + 0.5_DP * delta_cos

    ! Now loop over all the species and all the angles of this species

    DO ispecies = 1, nspecies

       DO iangles = 1, nangles(ispecies)

          ang_prob(iangles,ispecies)%prob(:) = 0.0_DP

          ! Determine the potential type for this angle

          angle_potential = angle_list(iangles,ispecies)%angle_potential_type

          IF (angle_potential == 'fixed') CYCLE

          IF ( angle_potential == 'harmonic' ) THEN
            
             ! There are two parameters. Note that these parameters
             ! are converted to appropriate units in input_routines.f90

             ktheta_0 = angle_list(iangles,ispecies)%angle_param(1)
             theta_0 = angle_list(iangles,ispecies)%angle_param(2)

          END IF

          ! initialize the total number of good regions

          ngood_regions = 0
          
          DO iregions = 1, nregions

             cos_iregions = costheta_0 + (iregions-1) * delta_cos

             ! obtain the angle whose cosine is cos_iregions

             theta_iregions = DACOS(cos_iregions)

             ! Determine the energy of this angle
             ! Note that the energy calcualtion will be different for
             ! different angle potential type.

             IF ( angle_potential == 'harmonic' ) THEN
                vtheta_iregions = ktheta_0 * (theta_iregions - theta_0) * (theta_iregions - theta_0)
             END IF

             ! protect against the overflow
             IF (  vtheta_iregions/(kboltz*T) < max_kBT ) THEN
                ! assign the probability 
                ngood_regions = ngood_regions + 1
                ang_prob(iangles,ispecies)%prob(ngood_regions) = EXP(-vtheta_iregions/(kboltz*T))
                ang_prob(iangles,ispecies)%theta(ngood_regions) = theta_iregions
      
             END IF
     
          


          END DO


          ! let us also store the number of good regions 
          ang_prob(iangles,ispecies)%nregions = ngood_regions

          ! sum the probabilities and normalize

          DO i = 2, ngood_regions

             ang_prob(iangles,ispecies)%prob(i) = ang_prob(iangles,ispecies)%prob(i-1) + &
                  ang_prob(iangles,ispecies)%prob(i)
          END DO
          
 
          
          ! the last term is the cumulative probability
          
          ang_prob(iangles,ispecies)%prob(:) = ang_prob(iangles,ispecies)%prob(:)/ &
               ang_prob(iangles,ispecies)%prob(ngood_regions) 
          
          
          

       END DO
       
    END DO
    
  END SUBROUTINE Angle_Distribution

!************************************************************************************************
  SUBROUTINE Pick_Angle(iangles,ispecies,theta_pick,prob_theta)
!************************************************************************************************
    ! The suroutine randomly picks a value for the bond angle for the input angle of the species
    ! It also returns the probability associated with the region in which the angle is located
    !*********************************************************************************************

    !**********************************************************************************************
    ! The program takes in the angle and species information from the calling program. To pick
    ! an angle, a random number is generated. Golden sampling method is then used to determine
    ! the region that will be picked. In this region, a value of theta is picked and returned to
    ! the calling program along with the probability of the region. Note that the routine
    ! uses the following vector that is determined in angle_distribution
    !
    ! ang_prob(iangles,ispecies)%theta
    ! ang_prob(iangles,ispecies)%prob(max_nregions)
    !**********************************************************************************************

    USE Type_Definitions
    USE Global_Variables
    USE Random_Generators
    USE File_Names

    IMPLICIT NONE

    INTEGER :: i, iangles, ispecies

    REAL(DP) :: random_no, prob_theta, theta_mid, cos_mid, delta_cos, cos_pick, theta_pick

    !--- Generate a random number that will be used to pick a region

    random_no = rranf()

    DO i = 1, ang_prob(iangles,ispecies)%nregions

       IF ( random_no < ang_prob(iangles,ispecies)%prob(i) ) EXIT

    END DO

    IF ( i == 1 ) THEN

       prob_theta = ang_prob(iangles,ispecies)%prob(1)
       
    ELSE

       prob_theta = ang_prob(iangles,ispecies)%prob(i) - ang_prob(iangles,ispecies)%prob(i-1)

    END IF

    ! We want to pick an angle in the selected region. We will do so from cosine distribution
    ! of angles

    ! determine the cosine of the midpoint

    theta_mid = ang_prob(iangles,ispecies)%theta(i)

    cos_mid = DCOS(theta_mid)

    delta_cos = 2.0_DP/nregions

    ! randomly pick a cosine around cos_mid with a width of delta_cos * 0.5

    cos_pick = cos_mid + ( rranf() - 0.5_DP ) * delta_cos

    theta_pick = DACOS(cos_pick)


  END SUBROUTINE Pick_Angle

!*************************************************************************************    
  SUBROUTINE Get_Theta_Prob(iangles,ispecies,this_theta,prob_this_theta)
!*************************************************************************************
    ! This routine takes in iangles, ispecies and this_theta from the calling program
    ! It returns theta_prob, the probability of observing the value of this_theta
    ! for iangles in ispecies
    !**********************************************************************************

    !**********************************************************************************
    ! To determine the probability, first cosine of the input angle is obtained as
    ! this will enable us to locate the region in which this_theta falls. Once the
    ! region is determined, its probability is computed by equating the midpoint angle
    ! of the region to(ang_prob(iangles,ispecies)%theta), the probablity
    ! of this region is then returned. Note that we cannot directly use the index of
    ! the region to query its probability as some of the intermediate regions might
    ! have a zero probability. It is useful to follow the code angle_distribution to
    ! see why this is the case.
    !**********************************************************************************

    USE Type_Definitions
    USE Global_Variables
    USE File_Names
    
    IMPLICIT NONE

    INTEGER :: iangles, ispecies, iregion, i
    
    REAL(DP) :: this_theta, prob_this_theta, costheta, delta_cos
    REAL(DP) :: cos_region, theta_region

    LOGICAL :: prob_assigned

    ! Initialize the logical variable to indicate that so far probability has
    ! not be assigned to this_theta. The truth of the variable to check later
    ! if we located the probability for the input angle

    prob_assigned = .FALSE.
    
    ! Determine the region in which this_theta falls
    
    costheta = DCOS(this_theta)
    
    delta_cos = 2.0_DP/nregions
    
    iregion = INT( (costheta+1.0)/ delta_cos ) + 1
    
    
    cos_region = -1.0_DP + 0.5_DP * delta_cos + (iregion - 1.0_DP) * delta_cos
    theta_region = DACOS(cos_region)
    
    ! Determine the probability of this region
    
    DO i = 1, ang_prob(iangles,ispecies)%nregions
       
       IF ( theta_region == ang_prob(iangles,ispecies)%theta(i) ) THEN
          
          prob_assigned = .TRUE.
          
          IF ( i == 1 ) THEN
             
             prob_this_theta = ang_prob(iangles,ispecies)%prob(1)
             EXIT
             
             
          ELSE
             
             prob_this_theta = ang_prob(iangles,ispecies)%prob(i) - &
                  ang_prob(iangles,ispecies)%prob(i-1)
             EXIT
             
          END IF
          
       END IF
       
    END DO
    
    ! check to make sure that we actually found the probability of this angle
    
    IF ( .NOT. prob_assigned ) THEN
       
       WRITE(logunit,*) 'Probability could not be assigned to', this_theta
       WRITE(logunit,*) 'For angle', iangles, 'in species', ispecies
       err_msg = ''
       CALL Clean_Abort(err_msg,'Get_Theta_Prob')
       
    END IF
    
  END SUBROUTINE Get_Theta_Prob
!**********************************************************************************

END MODULE Angle_Dist_Pick
    
