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

MODULE random_generators
! Generates a random number between 0 and 1 outlined in:
!
! L'Ecuyer, P. (1999) `Tables of maximally equidistributed combined LFSR
! generators', Math. of Comput., 68, 261-269.
!
! The cycle length is claimed to be ~ 2^(258) 
!
! The code is a modified version of the Fortran 90 version of the original C implementation
! of the L'Ecuyer Randomb number generator
! Modified by Andrew Paluch, 1 March 2009.
IMPLICIT NONE
! The intrinsic function "selected_real_kind" takes two arguments. The first is the number 
! of digits of precision desired, and the second is the largest magnitude of the exponent of 10.
INTEGER, PARAMETER:: da = SELECTED_REAL_KIND(14, 60)
! s1, s2, s3, s4, and s5 are the seeds to the random number generator, and are given default
! values in cast the seeds are not initialized by the user
INTEGER (KIND=8), SAVE  :: s1 = 153587801, s2 = -759022222, s3 = 1288503317, &
                           s4 = -1718083407, s5 = -123456789

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE init_seeds(i1, i3)

! Initialize the seeds to the random number generator.  If a call to this subroutine
! is forgotten, the default values will be used.
!
! For our simulations, we will specify only two of the five seeds, and leave the others
! as the default values

IMPLICIT NONE
INTEGER (KIND=8), INTENT(IN) :: i1, i3
s1 = i1
s3 = i3
IF (IAND(s1,      -2) == 0) s1 = i1 - 8388607
IF (IAND(s3,   -4096) == 0) s3 = i3 - 8388607
RETURN
END SUBROUTINE init_seeds

  !*****************************************************************************
  !*****************************************************************************
FUNCTION rranf()

! Returns a random number over the interval 0 to 1.
USE Type_Definitions, ONLY : DP
USE File_Names, ONLY : logunit

IMPLICIT NONE

REAL(DP) :: rranf

INTEGER (KIND=8) :: b

b  = ISHFT( IEOR( ISHFT(s1,1), s1), -53)
s1 = IEOR( ISHFT( IAND(s1,-2), 10), b)
b  = ISHFT( IEOR( ISHFT(s2,24), s2), -50)
s2 = IEOR( ISHFT( IAND(s2,-512), 5), b)
b  = ISHFT( IEOR( ISHFT(s3,3), s3), -23)
s3 = IEOR( ISHFT( IAND(s3,-4096), 29), b)
b  = ISHFT( IEOR( ISHFT(s4,5), s4), -24)
s4 = IEOR( ISHFT( IAND(s4,-131072), 23), b)
b  = ISHFT( IEOR( ISHFT(s5,3), s5), -33)
s5 = IEOR( ISHFT( IAND(s5,-8388608), 8), b)

! pconst is the reciprocal of (2^64 - 1)
rranf = IEOR( IEOR( IEOR( IEOR(s1,s2), s3), s4), s5) *5.4210108624275221E-20_DP + 0.5_DP

IF(rranf .GE. 1.0_DP) WRITE(logunit,*) 'rranf = 1.0'

END FUNCTION rranf


  ! Identical guts to rranf, but only return the seed.
  ! NOTE: randomization of sign bit is explicitly tossed out.
  INTEGER FUNCTION rranint()
      USE Global_Variables, ONLY : iseed

      ! LOCAL
      INTEGER :: ib, ia, ibc, ida, isum, iff, ie, ix, iy, ix2, ix1

      ib= iseed/65536
      ia= iseed - ib*65536
      ibc= ib*63253
      ida= ia*24301
      isum= ibc - 2147483647 + ida

      IF( isum .LE. 0 ) THEN
          isum= isum + 2147483647
      ELSE
          isum= isum - 1
      ENDIF

      iff = isum/32768

      ie = isum - iff*32768
      ix = ie + ia
      iy = 453816691 - 2283*ia
      ix2 = ix/32768
      ix1 = ix - 32768*ix2
      iseed = ix1*65536 - 2147483647 + iy

      IF( iseed .LE. 0 ) THEN
          iseed= iseed + 2147483647
      ELSE
          iseed = iseed - 1
      ENDIF

      rranint = iseed
  END FUNCTION rranint

  INTEGER FUNCTION random_range(low, high)
      INTEGER, INTENT(IN) :: low, high
      INTEGER :: binsize, binlimit, rangesize, ranval

      rangesize = high - low + 1

      binsize = 2147483647 / rangesize

      binlimit = rangesize*binsize

      ranval = rranint()

      ! Theoretically, this could run forever, but that's not a practical concern considering:
      ! a) the sequence is known to repeat
      ! b) the sequence is known to be uniform
      ! In combination, this dictates that the sequence will produce a valid number relatively
      ! quickly.
      DO WHILE (ranval .GT. binlimit)
          ranval = rranint()
      ENDDO

      ! Note that a remainder method (which this is NOT) will stress the low-order bits.  This
      ! method will tend to stress the high-order bits.  In a good generator, it shouldn't make
      ! much difference.  The more important impact here is that we avoid non-uniformity by having
      ! each conceptual "bin" be the same size.
      random_range = (ranval / binsize) + low
  END FUNCTION random_range



  !*****************************************************************************
   !-------------------------------------------------------------------------
  SUBROUTINE Generate_Random_Sphere(x_this,y_this,z_this)
  !------------------------------------------------------------------------
    ! This function generates points randomly on a unit sphere. It uses
    ! a single algorithm that the polar angle (the angle between the vector
    ! connecting the point and the center of the sphere and z-axis) is obtained
    ! from a cosine distribution while the azimuthal angle is sampled on the
    ! range [0,2pi]
    !
    ! Written by Jindal Shah on 07/26/10
    !
    !----------------------------------------------------------------------------

    USE Type_Definitions, ONLY : DP
    USE Global_Variables, ONLY : iseed, twoPI

    IMPLICIT NONE


    REAL(DP), INTENT(OUT) :: x_this, y_this, z_this

    REAL(DP) :: cos_theta, this_theta, phi, sin_theta

    x_this = 0.0_DP
    y_this = 0.0_DP
    z_this = 0.0_DP
    
    ! generate a polar angle, sample first from -1,1

    cos_theta = 2.0_DP * rranf() - 1.0_DP
    this_theta = DACOS(cos_theta)
    sin_theta = DSIN(this_theta)

    ! generate a random azimuthal angle, between 0,twopi

    phi = twoPI * rranf()

    x_this = DCOS(phi) * sin_theta
    y_this = DSIN(phi) * sin_theta

    z_this = cos_theta

  END SUBROUTINE Generate_Random_Sphere



end module random_generators
