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

PROGRAM Lowenstein_v1

!--------------------------------------------------------------------------------------------
! Implements the Lowenstein Rule for zeolites for a specified structure, initially containing
! Si, O, and (optionally) Al.
! Utilizes insertion and swap moves to achieve the specified ratio.
! Input: data containing xyz coordinates and element of each atom in the system
! Output: Lowenstein configuration of the system for a specified ratio
!
! Currently supports .pdb files
!
! Input arguments: file.pdb; desired ratio; use random seed (Y/N); outputfilename.xyz
! 
! Output: .xyz file of the Lowenstein configuration
!--------------------------------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER, PARAMETER :: SP = KIND(1.0e0)
  INTEGER, PARAMETER :: DP = KIND(1.0d0)

  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Coordinates_Matrix
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Actual_Neighbors
  CHARACTER(2), DIMENSION(:), ALLOCATABLE :: Element_Array, Lowenstein_Configuration
  LOGICAL, DIMENSION(:), ALLOCATABLE :: Lowenstein_Mask
  INTEGER :: maxAtomNumber, Si_Count, Al_Count, OpenStatus, ijk, ntries 
  REAL(DP) :: xlength, ylength, zlength, Lowenstein_Ratio, Desired_Ratio
  REAL(DP) :: alphaAngle, betaAngle, gammaAngle
  CHARACTER(120):: newFileName

  TYPE Lowenstein_List
    INTEGER :: atomIndex
    TYPE(Lowenstein_List), POINTER :: Next
  END TYPE Lowenstein_List
  TYPE(Lowenstein_List), POINTER :: Lowenstein_Si, Lowenstein_Al
  
  

  ! retrieve the data from the file
  CALL readInput(maxAtomNumber,xlength,ylength,zlength,alphaAngle,betaAngle,gammaAngle, &
                 Coordinates_Matrix,Element_Array, Desired_Ratio, newFileName)

  
  ! identify the neighbors
  CALL findNearestNeighbors(Actual_Neighbors)

  ! create the logical array denoting silicon sites permitted to participate in insertion/swap moves
  ! also creates the initial linked lists of AL and SI and returns Lowenstein_Si, Lowenstein_Al,
  ! and identifies the initial count of SI and AL present
  CALL createLowensteinStructures(Lowenstein_Mask, Lowenstein_Si, Lowenstein_Al, Si_Count, Al_Count,&
                                  Lowenstein_Configuration)


  ! compute an initial Lowenstein Ratio
  Lowenstein_Ratio = (1.0_DP * Si_Count / Al_Count)


  IF (Desired_Ratio .LT. 1.0) THEN
    STOP '*** Invalid Lowenstein Ratio. Valid ratios range from [1.0, inf) ***'
  END IF

  ! implements the procedure that will result in the desired Lowenstein Ratio
  CALL lowensteinProtocol(Lowenstein_Mask, Lowenstein_Ratio, Si_Count, Al_Count, &
                          Lowenstein_Si, Lowenstein_Al, Lowenstein_Configuration)

  ! If we failed to reach the specified ratio, then the program returns to this location
  ! If necessary, instruct the program to try again. Tries up to 20 times
  ntries = 1

  Try_Again: DO
    IF ((Lowenstein_Ratio .GT. 1.01_DP*Desired_Ratio) .AND. (ntries .LT. 20)) THEN
      ntries = ntries + 1
      ! reset the structures; reset Lowenstein_Ratio; call the process again
      CALL createLowensteinStructures(Lowenstein_Mask, Lowenstein_Si, Lowenstein_Al, Si_Count, Al_Count, &
                                      Lowenstein_Configuration)
      Lowenstein_Ratio = 1.0_DP * Si_Count / Al_Count
      CALL lowensteinProtocol(Lowenstein_Mask, Lowenstein_Ratio, Si_Count, Al_Count, &
                              Lowenstein_Si, Lowenstein_Al, Lowenstein_Configuration)
    ELSE IF (ntries .GE. 20) THEN
      STOP '*** Could not achieve specified Lowenstein Ratio - tried 20 times ***'
    ELSE
      EXIT Try_Again
    END IF
  END DO Try_Again

  ! Now that the specified ratio has been achieved, create an output xyz file.
  WRITE(*,801), 'Desired Lowenstein Ratio: ', Desired_Ratio
  WRITE(*,801), 'Actual Lowenstein Ratio: ', Lowenstein_Ratio
  801 FORMAT(A26,2X,F9.6)


  OPEN (UNIT = 25, FILE = newFileName, STATUS = "NEW", POSITION = "ASIS", &
         ACTION = "WRITE",  IOSTAT = OpenStatus)
  WRITE(25,800), maxAtomNumber
  800 FORMAT(1X,I6.1/)
  ijk = 0
  DO ijk = 1, maxAtomNumber
    Lowenstein_Configuration(ijk) = ADJUSTR(Lowenstein_Configuration(ijk)) ! return to original style
    WRITE(25,600, ADVANCE = "NO"), Lowenstein_Configuration(ijk), Coordinates_Matrix(ijk,1), Coordinates_Matrix(ijk,2), &
                   Coordinates_Matrix(ijk,3)
  END DO
  600 FORMAT(A2, 3(3X,F17.12)/)

  
 


  CONTAINS
  !-readInputFile--------------------------------------------------------------------------
  ! Prompts the user for a file containing the xyz coordinates and element of each atom
  ! in the system.  Returns the coordinates arranged by [(atom#,x), (atom#,y), (atom#,z)],
  ! the element for each index, the total number of atoms, and the dimensions of the unit
  ! cell.
  ! 
  ! Accepts: <filename.pdb>
  ! Output: Coordinates_Matrix, Element_Array, xlength, ylength, zlength, maxAtomNumber
  !----------------------------------------------------------------------------------------

    SUBROUTINE readInput (maxIndex,xmax,ymax,zmax,alphaAngle,betaAngle,gammaAngle, &
                          Coordinates,Element, Desired_Ratio, outFile)
                     
      REAL(DP), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: Coordinates
      REAL(DP), INTENT(OUT) :: xmax, ymax, zmax, alphaAngle, betaAngle, gammaAngle
      INTEGER, INTENT(OUT) :: maxIndex
      CHARACTER(2), INTENT(OUT), DIMENSION(:), ALLOCATABLE :: Element
      REAL(DP), INTENT(OUT) :: Desired_Ratio
      CHARACTER(120), INTENT(OUT) :: outFile
      INTEGER :: nline, OpenStatus, readStatus
      INTEGER, DIMENSION(:), ALLOCATABLE :: NumAtoms 
      CHARACTER :: FileName*120
      INTEGER, DIMENSION(:), ALLOCATABLE :: randomVarSeed
      INTEGER :: nSeed=12, fjk, hqk
      CHARACTER(1) :: initSeed

      WRITE(*, '(1X, A)', ADVANCE = "NO") "Enter the name of the initial data file: "
      READ *, FileName

      OPEN(UNIT = 15, FILE = FileName, STATUS = "OLD", IOSTAT = OpenStatus)
      IF (OpenStatus > 0) STOP "*** Cannot Open File ***"
    
      WRITE(*, '(1X, A)', ADVANCE = "NO") "Enter the desired final Lowenstein Ratio ( [1.0, inf)): "
      READ *, Desired_Ratio
      WRITE(*, '(1X, A)', ADVANCE = "NO") "Use random initial seed? (Y/N): "
      READ *, initSeed
      WRITE(*, '(1X, A)', ADVANCE = "NO") "Enter the name of the output file: "
      READ *, outFile

      IF (initSeed == 'Y') THEN
        CALL RANDOM_SEED(size = nSeed)
        ALLOCATE(randomVarSeed(nSeed))
        CALL SYSTEM_CLOCK(COUNT = fjk)
        randomVarSeed = fjk + 37 * (/ (hqk - 1, hqk = 1, nSeed) /)
        CALL RANDOM_SEED(PUT = randomVarSeed)
        DEALLOCATE(randomVarSeed)
      ELSE
        WRITE(*,*), 'Response other than "Y", using default seed'
        CALL RANDOM_SEED()
      END IF

      READ(15,350), xmax, ymax, zmax, alphaAngle, betaAngle, gammaAngle
      350 FORMAT(8X,2(F8.4,1X),F8.4,3(F7.3))

      nline = 0
      ! approximate the size of the file
      DO
        READ(15,*,IOSTAT = readStatus)
        IF (readStatus /= 0) THEN
          EXIT
        ELSE
          nline = nline + 1
        END IF
      END DO

      REWIND(15)
      READ(15,*)

      nline = nline + 1000 ! avoid error in the next loop for indices
      ALLOCATE(NumAtoms(nline))

      ! Assume it is a .pdb file since pdb's provide a,b,c, alpha, beta, gamma for the cell dimensions
      nline = 1
      ReadTheFile: DO
        READ(15,FMT = '(6X,I5.1)',END = 10,IOSTAT = readStatus), NumAtoms(nline)
        IF (readStatus /= 0) THEN
          EXIT
        END IF
    
        nline = nline + 1
      END DO ReadTheFile
      10 CONTINUE

      ! find the largest atom number encountered
      maxIndex = MAXVAL(NumAtoms)

      DEALLOCATE(NumAtoms)
      ALLOCATE (Coordinates(maxAtomNumber,3))
      ALLOCATE (Element(maxAtomNumber))

      REWIND(15)
      READ(15,*)

      DO nline = 1, maxAtomNumber
        READ(15,450,END=10, IOSTAT = readStatus), Coordinates(nline,1), Coordinates(nline,2), Coordinates(nline,3), &
                                     Element(nline)
        Element(nline) = ADJUSTL(Element(nline))
        ! ADJUSTL accounts for leading blanks which will occur for elements O, H, F, etc..
      END DO
      450 FORMAT(31X,3(F7.3,1X),21X,A2)

    END SUBROUTINE readInput
    


  !-findNearestNeighbors------------------------------------------------------------------------
  ! Finds the four nearest SI/AL for a given SI/AL identified in the original Element_Array.
  ! 
  ! Requires: xyz coordinates from Coordinates_Matrix, Element_Array, and the parameters of the
  !           unit cell.
  ! Returns: four closest SI/AL for a given SI/AL, assuming periodic boundary conditions.
  !---------------------------------------------------------------------------------------------

    SUBROUTINE findNearestNeighbors(Neighbors)

      INTEGER, INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: Neighbors
    
      LOGICAL, DIMENSION(:), ALLOCATABLE :: Mask_Previous_Closest
      REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Distance_Matrix
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: Closest_Twenty_Matrix
      INTEGER :: i,j,k,q
      REAL(DP) :: deltax, deltay, deltaz
    
      ALLOCATE (Neighbors(maxAtomNumber,4))
      ALLOCATE (Closest_Twenty_Matrix(maxAtomNumber,20)) ! Arbitrary choice of 20 nearest atoms
      ALLOCATE (Mask_Previous_Closest(maxAtomNumber))
      ALLOCATE (Distance_Matrix(maxAtomNumber,maxAtomNumber))
      Mask_Previous_Closest = .TRUE.

      ! calculate the minimum image distance and store in Distance_Matrix

      ! If the following statement is true, then we have a non-rectangular box
      IF ((alphaAngle /= 90.0) .OR. (betaAngle /= 90.0) .OR. (gammaAngle /= 90.0)) THEN
        CALL Minimum_Image_Separation(Distance_Matrix)  

      ELSE
        DO i = 1,maxAtomNumber-1
          DO j = i+1, maxAtomNumber
            deltax = MIN(ABS(ABS((Coordinates_Matrix(i,1))-Coordinates_Matrix(j,1))-xlength),&
                        ABS((Coordinates_Matrix(i,1))-Coordinates_Matrix(j,1)))
            deltay = MIN(ABS(ABS((Coordinates_Matrix(i,2))-Coordinates_Matrix(j,2))-ylength),&
                        ABS((Coordinates_Matrix(i,2))-Coordinates_Matrix(j,2)))
            deltaz = MIN(ABS(ABS((Coordinates_Matrix(i,3))-Coordinates_Matrix(j,3))-zlength),&
                        ABS((Coordinates_Matrix(i,3))-Coordinates_Matrix(j,3)))
            Distance_Matrix(i,j) = (deltax**2) + (deltay**2) + (deltaz**2)
            Distance_Matrix(j,i) = Distance_Matrix(i,j)
          END DO
        END DO
      END IF

      ! Using the above information, find the actual neighboring SI/AL for each SI/AL
      ! first, retrieve the ten lowest values identified for a given atom in Distance_Matrix
      Find_Closest:  DO i = 1, maxAtomNumber
        Mask_Previous_Closest = .TRUE.! reset for next iteration 'i'
        DO j = 1,20
          Closest_Twenty_Matrix(i,j) = MINLOC(Distance_Matrix(i,:),DIM=1,MASK=Mask_Previous_Closest(:))
          Mask_Previous_Closest(Closest_Twenty_Matrix(i,j)) = .FALSE. ! do not check this index in the future
        END DO
      END DO Find_Closest

      Find_Actual_Neighbors:  DO i = 1, maxAtomNumber
        IF ((Element_Array(i) .EQ. 'SI') .OR. (Element_Array(i) .EQ. 'AL')) THEN
          k = 0 ! initialize a counter
          q = 2 ! skip the first atom in 'Closest_Ten_Matrix'
          DO WHILE (k .LT. 4) ! need 4 closest SI/AL other than itself
            IF ((Element_Array(Closest_Twenty_Matrix(i,q)) .EQ. 'SI') .OR. &
                (Element_Array(Closest_Twenty_Matrix(i,q)) .EQ. 'AL')) THEN
              k = k + 1
              Neighbors(i,k) = Closest_Twenty_Matrix(i,q)
              q = q + 1
            ELSE
              q = q + 1
            END IF
        
            IF (q .GE. 22) THEN
              STOP '*** Could not find 4 silicons nearby ***'
            END IF
          END DO

        ELSE IF (Element_Array(i) .EQ. 'O') THEN
          CYCLE
        END IF

      END DO Find_Actual_Neighbors
  
      DEALLOCATE(Distance_Matrix)
      DEALLOCATE(Closest_Twenty_Matrix)
      DEALLOCATE(Mask_Previous_Closest)
      
    END SUBROUTINE findNearestNeighbors


  !-Minimum_Image_Separation------------------------------------------------------------------------------
  ! Finds the minimum image distance for a given pair of particles in a non-rectangular box
  !
  ! Requires: Coordinates_Matrix, a,b,c,alphaAngle,betaAngle,gammaAngle
  ! 
  ! Returns: Distance_Matrix
  !-------------------------------------------------------------------------------------------------------

    SUBROUTINE Minimum_Image_Separation(Distance_Matrix)
    
      REAL(DP),INTENT(INOUT), DIMENSION(:,:), ALLOCATABLE :: Distance_Matrix
      REAL(DP), DIMENSION(3,3) :: tM, tMAdj, tMInv

      REAL(DP), PARAMETER :: pi_cst = 3.1415926536_DP
      REAL(DP), PARAMETER :: degree_to_rad = pi_cst / 180_DP
      INTEGER :: ii,jj
      REAL(DP) :: det, deltax, deltay, deltaz, inv_det
      REAL(DP) :: temp1, temp2, temp3 

      ! tM = transformation matrix
      ! tMInv = transformation matrix inverse
      ! tMAdj = transformation matrix adjoint
      ! refer to p. 170, Macedonia appendix
      tM(1,1) = ABS(xlength)
      tM(1,2) = 0.0_DP
      tM(1,3) = 0.0_DP
      tM(2,1) = (ABS(ylength)) * COS(gammaAngle*degree_to_rad)
      tM(2,2) = (ABS(ylength)) * SIN(gammaAngle*degree_to_rad)
      tM(2,3) = 0.0_DP
      tM(3,1) = (ABS(zlength)) * COS(betaAngle*degree_to_rad)
      tM(3,2) = ((ABS(ylength)*ABS(zlength)*COS(alphaAngle*degree_to_rad)) - (tM(2,1)*tM(3,1))) / tM(2,2)
      tM(3,3) = SQRT(((ABS(zlength))**2) - (tM(3,1) ** 2) - (tM(3,2)**2))

      ! check for numerical imprecision
      ! arbitrary tolerance of 10E-05
      DO ii = 1,3
        DO jj = 1,3
          IF (ABS(tM(jj,ii)) .LT. .00001_DP) THEN
            tM(jj,ii) = 0.0_DP
          END IF
        END DO
      END DO

      ! Compute the adjoint of the transformation matrix
      ! these are the values stored in box_list(box_nbr)%length_inv
      tMAdj(1,1) = (tM(2,2) * tM(3,3)) - (tM(2,3) * tM(3,2))
      tMAdj(2,1) = (tM(3,1) * tM(2,3)) - (tM(2,1) * tM(3,3))
      tMAdj(3,1) = (tM(2,1) * tM(3,2)) - (tM(2,2) * tM(3,1))
      tMAdj(1,2) = (tM(3,2) * tM(1,3)) - (tM(1,2) * tM(3,3))
      tMAdj(2,2) = (tM(1,1) * tM(3,3)) - (tM(3,1) * tM(1,3))
      tMAdj(3,2) = (tM(3,1) * tM(1,2)) - (tM(1,1) * tM(3,2))
      tMAdj(1,3) = (tM(1,2) * tM(2,3)) - (tM(2,2) * tM(1,3))
      tMAdj(2,3) = (tM(2,1) * tM(1,3)) - (tM(1,1) * tM(2,3))
      tMAdj(3,3) = (tM(1,1) * tM(2,2)) - (tM(2,1) * tM(1,2))

      det = tM(1,1) * tMAdj(1,1) + tM(1,2) * tMAdj(2,1) + tM(1,3) * tMAdj(3,1)
 
      inv_det = 1.0_DP / det
      tMInv = tMAdj * inv_det

      ! Now, implement the actual calculations to find the minimum image. return Distance_Matrix
      DO ii = 1, maxAtomNumber - 1
        DO jj = ii + 1, maxAtomNumber
          deltax = Coordinates_Matrix(ii,1) - Coordinates_Matrix(jj,1)
          deltay = Coordinates_Matrix(ii,2) - Coordinates_Matrix(jj,2)
          deltaz = Coordinates_Matrix(ii,3) - Coordinates_Matrix(jj,3)
          ! convert to fractional coordinate system
          temp1 = tMInv(1,1)*deltax + tMInv(1,2)*deltay + tMInv(1,3)*deltaz
          temp2 = tMInv(2,1)*deltax + tMInv(2,2)*deltay + tMInv(2,3)*deltaz
          temp3 = tMInv(3,1)*deltax + tMInv(3,2)*deltay + tMInv(3,3)*deltaz
          ! apply periodic boundary conditions
          temp1 = temp1 - REAL(NINT(temp1,DP))
          temp2 = temp2 - REAL(NINT(temp2,DP))
          temp3 = temp3 - REAL(NINT(temp3,DP))
          ! convert to Cartesian Coordinates
          deltax = tM(1,1)*temp1 + tM(1,2)*temp2 + tM(1,3)*temp3
          deltay = tM(2,1)*temp1 + tM(2,2)*temp2 + tM(2,3)*temp3
          deltaz = tM(3,1)*temp1 + tM(3,2)*temp2 + tM(3,3)*temp3
          ! get the distance for this pair of atoms
          Distance_Matrix(ii,jj) = SQRT((deltax**2) + (deltay**2) + (deltaz**2))
          Distance_Matrix(jj,ii) = Distance_Matrix(ii,jj)
        END DO
      END DO 

    END SUBROUTINE Minimum_Image_Separation

  !-createLowensteinStructures----------------------------------------------------------------------------
  ! Creates the linked lists for silicon and aluminum.  The linked list for aluminum contains all indices
  ! identified as aluminum sites; the linked list for silicon contains only sites identified as valid for
  ! participation in swap and insertion moves.  These sites are marked 'TRUE' in the logical mask; all 
  ! other silicons are marked false.  
  ! 
  ! Requires: Element_Array, Actual_Neighbors, maxIndex
  !          
  ! Returns: Lowenstein_Mask, Lowenstein_Si, Lowenstein_Al, Si_Count, Al_Count
  !-------------------------------------------------------------------------------------------------------

    SUBROUTINE createLowensteinStructures(logicalMask,siList,alList,nSI,nAL,config)

      INTEGER, INTENT(OUT) :: nSI, nAL
      LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: logicalMask
      TYPE(Lowenstein_List), POINTER, INTENT(OUT) :: siList, alList
      CHARACTER(2), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: config

      INTEGER :: i, j
      TYPE(Lowenstein_List), POINTER :: P1, P2

      IF (ALLOCATED(logicalMask)) THEN
        DEALLOCATE(logicalMask)
      END IF

      ! create the initial mask
      ALLOCATE(logicalMask(maxAtomNumber))
      logicalMask = .TRUE.

      adjustMask: DO i = 1, maxAtomNumber
        IF ((Element_Array(i) .EQ. 'O') .OR. (Element_Array(i) .EQ. 'AL')) THEN
          logicalMask(i) = .FALSE.
        ELSE IF (Element_Array(i) .EQ. 'SI') THEN
          DO j = 1,4
            IF (Element_Array(Actual_Neighbors(i,j)) == 'AL') THEN
              logicalMask(i) = .FALSE.
            ELSE
              CYCLE
            END IF
          END DO
        ELSE
          STOP '**** This script only recognizes elements Al, Si, and O at this time ***'
        END IF
      END DO adjustMask

      ! count the number of silicon and aluminum present and create the linked list(s)
      NULLIFY(P1)
      NULLIFY(P2)
      NULLIFY(alList)
      NULLIFY(siList)

      nSI = 0
      nAL = 0
      linkedLists: DO i = 1, maxAtomNumber
        IF (Element_Array(i) == 'SI') THEN
          nSI = nSI + 1
          ! see if we can add this SI to linked list
          IF (logicalMask(i) .EQV. .TRUE.) THEN
            IF (.NOT. ASSOCIATED(siList)) THEN
              ALLOCATE(siList)
              siList%atomIndex = i
              NULLIFY(siList%Next)
            ELSE
              ALLOCATE(P1)
              P1%atomIndex = i
              P1%Next => siList
              siList => P1
            END IF
          END IF
        ELSE IF (Element_Array(i) == 'AL') THEN
          nAL = nAL + 1
          IF (.NOT. ASSOCIATED(alList)) THEN
            ALLOCATE(alList)
            alList%atomIndex = i
            NULLIFY(alList%Next)
          ELSE
            ALLOCATE(P2)
            P2%atomIndex = i
            P2%Next => alList
            alList => P2
          END IF
        END IF
      END DO linkedLists

      ! Deallocate the pointers. note that some input files will not have any aluminum - 
      ! therefore, check for P2 otherwise a runtime error will occur
      DEALLOCATE(P1)
      IF (ASSOCIATED(P2)) THEN
        DEALLOCATE(P2)
      END IF

      IF (ALLOCATED(config)) THEN
        DEALLOCATE(config)
      END IF

      ALLOCATE(config(maxAtomNumber))
      config = Element_Array

    END SUBROUTINE createLowensteinStructures

    
  !-lowensteinProtocol---------------------------------------------------------------------------------------
  ! Inserts aluminum at valid silicon sites while the calculated Lowenstein Ratio is greater than the
  ! desired ratio as indicated by the user.  Achieves this using insertions and swap moves.
  ! Manages the iterations
  ! Accepts: Element_Array, Actual_Neighbors, Desired_Ratio, Si_Count, Al_Count, Lowenstein_Ratio
  !          Lowenstein_Si, Lowenstein_Al
  ! Returns: Lowenstein_Si, Lowenstein_Al, Lowenstein_Ratio, Lowenstein_Configuration, Si_Count, Al_Count
  !----------------------------------------------------------------------------------------------------------

    SUBROUTINE lowensteinProtocol (logicalMask, actual, nSI, nAL, listSI, listAL, config)
      
      CHARACTER(2), DIMENSION(:), INTENT(INOUT) :: config
      LOGICAL, DIMENSION(:), INTENT(INOUT) :: logicalMask
      REAL(DP), INTENT(INOUT) :: actual
      INTEGER, INTENT(INOUT) :: nSI, nAL
      TYPE(Lowenstein_List), POINTER, INTENT(INOUT) :: listSI, listAL

      INTEGER :: i, numInList, selectAL, selectSI, maxIter
      REAL(DP) :: swapInsertParam
      REAL(DP), DIMENSION(:), ALLOCATABLE :: swapOrInsert, randSelectSI, randSelectAL
      CHARACTER :: caseID*6, elementID*2
      

      maxIter = 10*maxAtomNumber !arbitrary

      IF (ALLOCATED(swapOrInsert)) THEN
        DEALLOCATE(swapOrInsert)
      END IF

      IF (ALLOCATED(randSelectSI)) THEN
        DEALLOCATE(randSelectSI)
      END IF

      IF (ALLOCATED(randSelectAL)) THEN
        DEALLOCATE(randSelectAL)
      END IF

      ALLOCATE(swapOrInsert(maxIter))
      ALLOCATE(randSelectSI(maxIter))
      ALLOCATE(randSelectAL(maxIter))

      ! set up our random numbers
      CALL RANDOM_NUMBER(swapOrInsert)
      CALL RANDOM_NUMBER(randSelectAL)
      CALL RANDOM_NUMBER(randSelectSI)
   
      ! create the configuration
      i = 0
      swapInsertParam = 1.00

      ! while the calculated ratio is greater than the desired ratio, continue swaps and insertions
      DO WHILE (actual .GT. Desired_Ratio)

        i = i + 1
         
        numInList = COUNT(logicalMask) ! counts the number of valid silicon sites for this iteration

        ! Now, change the ratio of swaps vs insertions according to the current state of the ratio.
        ! At lower ratios, we will require more swap moves.
        IF ((Lowenstein_Ratio .LT. 2.0_DP)) THEN
          swapInsertParam = .02_DP
        ELSE IF (Lowenstein_Ratio .LT. 3.0_DP) THEN
          swapInsertParam = .05_DP
        ELSE IF (Lowenstein_Ratio .LT. 5.0_DP) THEN
          swapInsertParam = .10_DP
        ELSE IF (Lowenstein_Ratio .LT. 8.0_DP) THEN
          swapInsertParam = .15_DP
        ELSE
          swapInsertParam = 1.0_DP
        END IF

        ! taking too many moves to achieve the specified ratio
        IF (i .GE. (maxIter - 1)) THEN
          WRITE(*,*) '*** Reached Maximum Number of Iterations - lowensteinProtocol subroutine ***'
          WRITE(*,*) '*** Increase maxIter in lowensteinProtocol to continue ***'
          RETURN
        END IF
        
        ! our configuration does not have any available silicon sites for swaps or insertions
        IF (numInList .EQ. 0) THEN
          RETURN
        END IF

        ! the silicon that is permitted to be swapped or participate in insertion is randomly selected
        selectSI = NINT(numInList*(randSelectSI(i)))
        ! sometimes this might be zero. if so, just set it to 1
        IF (selectSI .LT. 1) THEN
          selectSI = 1
        END IF

        ! Now, identify whether this iteration will be a swap move or an insertion
        IF ((swapOrInsert(i) .GT. swapInsertParam) .AND. (nAL .GT. 0)) THEN
          caseID = 'Swap'
        ELSE
          caseID = 'Insert'
        END IF
       
        SELECT CASE (caseID)
        
          CASE ('Swap')
            selectAL = NINT(nAL*(randSelectAL(i)))
            ! insert an aluminum
            elementID = 'SI' 
            CALL lowensteinInsert(selectSI,listSI,listAL,config,elementID,logicalMask)
            ! insert a silicon
            elementID = 'AL'
            CALL lowensteinInsert(selectAL,listSI,listAL,config,elementID,logicalMask) 

            ! we now have our two atoms; fix the lists, and repair the mask
            CALL fixStructures(selectSI,listSI,listAL,logicalMask)
            CALL fixStructures(selectAL,listSI,listAL,logicalMask)

          CASE ('Insert')
            elementID = 'SI'
            CALL lowensteinInsert(selectSI,listSI,listAL,config,elementID,logicalMask)
 
            ! now, fix the list, and repair the mask
            CALL fixStructures(selectSI,listSI,listAL,logicalMask)
            
            nSI = nSI - 1
            nAL = nAL + 1

        END SELECT

        actual = 1.0 * Si_Count / Al_Count
      END DO
      
      DEALLOCATE(swapOrInsert)
      DEALLOCATE(randSelectSI)
      DEALLOCATE(randSelectAL)

    END SUBROUTINE lowensteinProtocol

    !-lowensteinInsert------------------------------------------------------------------------------------
    ! Implements a swap/insert move, where a randomly selected aluminum is interchanged with a silicon atom
    ! selected from the list of eligible sites.  This move affects the overall configuration, but
    ! not the actual Lowenstein Ratio of the cell.  A swap move is actually two insertions at partner
    ! sites, with no change in the overall number of silicon or aluminum in the system.
    !
    ! Accepts: indexI, (selectSI or selectAL) , listSI, listAL, config, nSI, nAL, caseID,elementID, &
    !          logicalMask
    !
    ! Returns: listSI, listAL, config, logicalMask
    !---------------------------------------------------------------------------------------------------

    SUBROUTINE lowensteinInsert (selectedAtom, listSI, listAL, config, &
                                 elementID,logicalMask)

 
      TYPE(Lowenstein_List), POINTER, INTENT(INOUT) :: listSI, listAL
      CHARACTER(2), DIMENSION(:), INTENT(INOUT) :: config
      CHARACTER, INTENT(IN) :: elementID*2
      LOGICAL, DIMENSION(:), INTENT(INOUT) :: logicalMask
      INTEGER, INTENT(INOUT) :: selectedAtom
      TYPE(Lowenstein_List), POINTER :: P1

      INTEGER :: j

      NULLIFY(P1)
      ALLOCATE(P1)
      ! check the elementID, and traverse the appropriate list to the selectedAtom
      IF (elementID == 'SI') THEN
        P1 => listSI
      ELSE IF (elementID == 'AL') THEN
        P1 => listAL
      END IF

      IF (selectedAtom /=1) THEN
          DO j = 1, selectedAtom-1
            P1 => P1%Next
          END DO
      END IF

      ! P1 will now be pointing towards our selected atom's actual index within the structure
      selectedAtom = P1%atomIndex
      NULLIFY(P1)
      ! insert the atom and alter config and the logicalMask
      IF (config(selectedAtom) == 'SI') THEN
        config(selectedAtom) = 'AL'
        logicalMask(selectedAtom) = .FALSE.
      ELSE IF (config(selectedAtom) == 'AL') THEN
        config(selectedAtom) = 'SI'
        logicalMask(selectedAtom) = .TRUE. ! note that this is inherently true
      END IF

    END SUBROUTINE lowensteinInsert

    !-fixStructures-------------------------------------------------------------------------------------
    ! Prepares the Lowenstein_SI, Lowenstein_AL, and Lowenstein_Mask for the next iteration of moves
    ! Adjusts the arrays and linked lists to account for the changes in configuration.
    !
    ! Accepts: index i, selectedAtom, listSI, listAL, logicalMask
    ! Returns: fixed listSI, listAL, and logicalMask
    !---------------------------------------------------------------------------------------------------

    SUBROUTINE fixStructures (selectedAtom,listSI,listAL,logicalMask)

      INTEGER, INTENT(IN) :: selectedAtom ! this is now an atomIndex
      TYPE(Lowenstein_List), POINTER, INTENT(INOUT) :: listSI, listAL
      LOGICAL, DIMENSION(:), INTENT(INOUT) :: logicalMask
 
      INTEGER :: j,k, neighborAtom, thisNeighbor

      ! at this point, config is updated with the new SI/AL indices.  However, the lists and mask are not
      ! necessarily up to date.

      ! This atom used to be an aluminum. Therefore, we must remove from listAL, add to listSI,
      ! and check whether the neighboring atoms are now eligible to be on listSI (consult & update mask)

      IF (Lowenstein_Configuration(selectedAtom) == 'SI') THEN
        CALL removeFromList(selectedAtom,listAL) 
        logicalMask(selectedAtom) = .TRUE.
        CALL addToList(selectedAtom,listSI)
   
        ! check the neighbors of the neighboring SI/AL
        DO j = 1,4
          neighborAtom = Actual_Neighbors(selectedAtom,j)
          logicalMask(neighborAtom) = .TRUE.
          DO k = 1,4
            IF (Lowenstein_Configuration(Actual_Neighbors(neighborAtom,k)) == 'AL') THEN
              logicalMask(neighborAtom) = .FALSE.
              EXIT
            ELSE
              CYCLE
            END IF
          END DO

          IF (logicalMask(Actual_Neighbors(selectedAtom,j)) .EQV. .TRUE.) THEN
            CALL addToList(Actual_Neighbors(selectedAtom,j),listSI)
          END IF
        END DO


      ! This atom used to be a silicon.  Therefore, remove from listSI, add to listAL, and remove the
      ! neighboring silicons from listSI as necessary (consult and update mask)
      ELSE IF (Lowenstein_Configuration(selectedAtom) == 'AL') THEN 
        CALL removeFromList(selectedAtom,listSI)
        CALL addToList(selectedAtom,listAL)

        DO j = 1,4
          thisNeighbor = Actual_Neighbors(selectedAtom,j)
          IF (logicalMask(thisNeighbor) .EQV. .TRUE.) THEN
            logicalMask(thisNeighbor) = .FALSE.
            CALL removeFromList(thisNeighbor,listSI)
          END IF
        END DO
      END IF 

    END SUBROUTINE fixStructures

    !-removeFromList--------------------------------------------------------------------------------------
    ! Removes a specified atom from a provided linked list
    !
    ! Accepts: atom index, linked list
    ! Returns: amended linked list
    !-----------------------------------------------------------------------------------------------------

    SUBROUTINE removeFromList (thisAtom, thisList)
      INTEGER, INTENT(IN) :: thisAtom
      TYPE(Lowenstein_List), POINTER, INTENT(INOUT) :: thisList

      TYPE(Lowenstein_List), POINTER :: P1, P2
      INTEGER :: Counter1


      NULLIFY(P1)
      NULLIFY(P2)
      P1 => thisList
      P2 => thisList 

      Counter1 = 0
      DO
        Counter1 = Counter1 + 1 
        IF (P1%atomIndex == thisAtom) THEN 
          IF (Counter1 == 1) THEN
            thisList => P1%Next
            NULLIFY(P1)
          ELSE IF (Counter1 .GT. 1) THEN
            P2%Next => P1%Next
            NULLIFY(P1)
          END IF
          EXIT

        ELSE
         
          IF (Counter1 .GT. 1) THEN
            P2 => P1
            P1 => P1%Next
          ELSE
            P1 => P1%Next
          END IF
        END IF
      END DO

      NULLIFY(P1)
      NULLIFY(P2)


      END SUBROUTINE removeFromList

      !-addToList-----------------------------------------------------------------------------------------
      ! Adds a specified index to a provided linked list
      !
      ! Accepts: atom index, linked list
      ! Returns: amended linked list
      !---------------------------------------------------------------------------------------------------

      SUBROUTINE addToList (thisAtom, listToAppend)
        INTEGER, INTENT(IN) :: thisAtom
        TYPE(Lowenstein_List), POINTER, INTENT(INOUT) :: listToAppend
        
        TYPE(Lowenstein_list), POINTER :: P6
 
        NULLIFY(P6)
        ! assume we don't know if the list has been created yet
        IF (.NOT. ASSOCIATED(listToAppend)) THEN
          ALLOCATE(listToAppend)
          listToAPpend%atomIndex = thisAtom
          NULLIFY(listToAppend%Next)
        ELSE
          ALLOCATE(P6)
          P6%atomIndex = thisAtom
          P6%Next => listToAppend
          listToAppend => P6
        END IF
        
        IF (ASSOCIATED(P6)) THEN
          NULLIFY(P6)
        END IF

      END SUBROUTINE addToList
 
END PROGRAM Lowenstein_v1
