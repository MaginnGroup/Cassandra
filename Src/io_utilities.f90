
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

MODULE IO_Utilities

  !***************************************************************************
  ! This module contains a colection of small routines used to read information 
  ! from file and manipulate them.
  !
  ! Called by
  !
  !   chempot
  !   compute_cell_dimensions
  !   create_nonbond_table
  !   cut_n_grow
  !   fragment_growth
  !   gcmc_control
  !   gemc_control
  !   gemc_driver
  !   gemc_nvt_volume
  !   gemc_particle_transfer
  !   input_routines
  !   insertion
  !   main
  !   mcf_control
  !   nptmc_control
  !   nvtmc_control
  !   fragment_control
  !   ring_fragment_driver
  !   participation
  !   volume_change
  !
  ! Revision history
  !
  !    12/10/13  : Beta version
  !********************************************************************************
  USE Global_Variables

  IMPLICIT NONE

CONTAINS

!********************************************************************************
  SUBROUTINE Read_String(file_number,string,ierr)
!********************************************************************************
! This routine just reads a single string from a file with unit number
! equal to file_number and returns that 120 character string. It returns
! ierr .ne. 0 if it couldn't read the file.
!********************************************************************************
    INTEGER, INTENT(IN) :: file_number
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(120), INTENT(OUT) :: string
!********************************************************************************

    READ(file_number,'(A120)',IOSTAT=ierr) string
    IF (ierr .NE. 0) RETURN

  END SUBROUTINE Read_String


!********************************************************************************
  SUBROUTINE Parse_String(file_number,line_nbr,min_entries,nbr_entries,line_array,ierr) 
!********************************************************************************
! This routine reads one line from the file file_number. It reads the total number
! of entries on the line and places the entries in the character array line_array
! in consecutive order. It skips leading blanks, and determines if an entry is 
! different by detecting a space between entries. It also tests to see if the 
! minimum number of entries specified was met or not. If not, and error is returned.
!********************************************************************************
    CHARACTER(120), INTENT(OUT) :: line_array(20)
    INTEGER, INTENT(IN) :: file_number,min_entries,line_nbr
    INTEGER, INTENT(OUT) :: nbr_entries
    INTEGER, INTENT(INOUT) :: ierr
    
    CHARACTER(120) :: string
    INTEGER :: line_position,i
    LOGICAL :: space_start
!********************************************************************************

! Zero counter for number of entries found
    nbr_entries = 0

! Counter for positon on the line
    line_position = 1

! clear entry array
    line_array = ""
      
! Read the string from the file
    CALL Read_String(file_number,string,ierr)
    
    IF (string(1:1) .NE. ' ') THEN

       ! first character is an entry, so advance counter
       nbr_entries = nbr_entries + 1
    ENDIF

    space_start = .FALSE. 

    ! Recall: LEN_TRIM() is length of string, not counting trailing blanks
    DO i=1,LEN_TRIM(string) 
        IF (string(i:i) .EQ. ' ') THEN
          IF (.NOT. space_start) space_start = .TRUE.
          ! This means a new set of spaces has been found
       ELSE
          IF (space_start) THEN
             nbr_entries = nbr_entries + 1
             line_position = 1
          ENDIF
          space_start = .FALSE.
          line_array(nbr_entries)(line_position:line_position) = &
               string(i:i)
          line_position = line_position + 1
       ENDIF
       
    ENDDO

    ! Test to see if the minimum number of entries was read in      
    IF (nbr_entries < min_entries) THEN
       err_msg = ""
       err_msg(1) = 'Error attempting to parse line ' // &
                    TRIM(Int_To_String(line_nbr)) // ' of input file: '
       err_msg(2) = TRIM(string)
       err_msg(3) = 'into at least ' // TRIM(Int_To_String(min_entries)) // ' entries'
       CALL Clean_Abort(err_msg,'Parse_String')
    END IF
      
    END SUBROUTINE Parse_String

!********************************************************************************
  SUBROUTINE Parse_String_Zeolite_Frag(file_number,line_nbr,min_entries,nbr_entries,line_array_zeo,ierr) 
!********************************************************************************
! This routine reads one line from the file file_number. It reads the total number
! of entries on the line and places the entries in the character array line_array
! in consecutive order. It skips leading blanks, and determines if an entry is 
! different by detecting a space between entries. It also tests to see if the 
! minimum number of entries specified was met or not. If not, and error is returned.
!********************************************************************************
    CHARACTER(50000), INTENT(OUT) :: line_array_zeo(10000)
    INTEGER, INTENT(IN) :: file_number,min_entries,line_nbr
    INTEGER, INTENT(OUT) :: nbr_entries
    INTEGER, INTENT(INOUT) :: ierr
    
    CHARACTER(50000) :: string
    INTEGER :: line_position,i
    LOGICAL :: space_start
!********************************************************************************

! Zero counter for number of entries found
    nbr_entries = 0

! Counter for positon on the line
    line_position = 1

! clear entry array
    line_array_zeo = ""
      
! Read the string from the file
    CALL Read_String_Zeo(file_number,string,ierr)
    IF (string(1:1) .NE. ' ') THEN

       ! first character is an entry, so advance counter
       nbr_entries = nbr_entries + 1
    ENDIF

    space_start = .FALSE. 

    ! Recall: LEN_TRIM() is length of string, not counting trailing blanks
    DO i=1,LEN_TRIM(string) 
       IF (string(i:i) .EQ. ' ') THEN
          IF (.NOT. space_start) space_start = .TRUE.
          ! This means a new set of spaces has been found
       ELSE
          IF (space_start) THEN
             nbr_entries = nbr_entries + 1
             line_position = 1
          ENDIF
          space_start = .FALSE.
          line_array_zeo(nbr_entries)(line_position:line_position) = &
               string(i:i)
          line_position = line_position + 1
       ENDIF
       
    ENDDO

    ! Test to see if the minimum number of entries was read in      
    IF (nbr_entries < min_entries) THEN
       err_msg = ""
       err_msg(1) = 'Expected at least '// TRIM(Int_To_String(min_entries))//&
            ' input(s) on line '//TRIM(Int_To_String(line_nbr))//' of input file.'
       CALL Clean_Abort(err_msg,'Parse_String_Zeolite_Frag')
    END IF
      
    END SUBROUTINE Parse_String_Zeolite_Frag

!********************************************************************************
  SUBROUTINE Read_String_Zeo(file_number,string,ierr)
!********************************************************************************
! This routine just reads a single string from a file with unit number
! equal to file_number and returns that 150 character string. It returns
! ierr .ne. 0 if it couldn't read the file.
!********************************************************************************
    INTEGER, INTENT(IN) :: file_number
    INTEGER, INTENT(OUT) :: ierr
    CHARACTER(50000), INTENT(OUT) :: string
!********************************************************************************

    READ(file_number,'(A50000)',IOSTAT=ierr) string
    IF (ierr .NE. 0) RETURN

END SUBROUTINE Read_String_Zeo


!****************************************************************************
    SUBROUTINE Name_Files(prefix,suffix,new_name)
!****************************************************************************
      ! This routine takes two strings as input and outputs the concatanation
      ! of the two strings. We use it to append a suffix onto the main run
      ! filename (the prefix).
!****************************************************************************
      IMPLICIT NONE

      INTEGER :: prefix_length
      CHARACTER(120) :: prefix,new_name
      CHARACTER(*) :: suffix

      prefix_length = LEN_TRIM(prefix)

      new_name = prefix(1:prefix_length) // suffix

    END SUBROUTINE Name_Files


!****************************************************************************
    FUNCTION Int_To_String(int_in)
!****************************************************************************
      IMPLICIT NONE
  !**************************************************************************
  !                                                                         *
  ! This function takes an integer argument
  ! and returns the character equivalent
  !                                                                         *
  !**************************************************************************
  CHARACTER(40) :: int_to_string
  LOGICAL :: is_negative
  INTEGER :: int_in, chop_int, ndigits, curr_digit
  
  int_to_string = ""
  ndigits = 0
  
  !Check to see if integer is zero
  IF (int_in == 0) THEN
     int_to_string(1:1) = "0"
     RETURN
  END IF
  
  !Determine if integer is negative
  IF (int_in < 0) THEN
     is_negative = .TRUE.
     chop_int = -int_in
  ELSE
     is_negative = .FALSE.
     chop_int = int_in
  END IF
  
  !Pop of last digit of integer and fill in string from right to
  !left
  DO
     ndigits = ndigits + 1
     curr_digit = MOD(chop_int,10)
     int_to_string(41-ndigits:41-ndigits) = ACHAR(curr_digit+48)
     chop_int = INT(chop_int / 10.0d0)
     IF (chop_int == 0) EXIT
  END DO
  
  IF (is_negative) int_to_string(40-ndigits:40-ndigits) = "-"
  
  !Left justify string
  int_to_string = ADJUSTL(int_to_string)
  
END FUNCTION Int_To_String


!****************************************************************************
FUNCTION String_To_Double(string_in)
!****************************************************************************
! This function takes a character string in and outputs the equivalent DP 
  ! number.
!****************************************************************************
  USE Type_Definitions
  IMPLICIT NONE

  LOGICAL :: dec_found, exp_found, is_negative
  CHARACTER(*) :: string_in
  CHARACTER(50) :: cff_bd, cff_ad, expv
  INTEGER :: nchars, strln, exp_start,dec_start
  INTEGER :: ii, cnt, icff_bd, iexpv
  REAL(DP) :: string_to_double, div, add_num
!****************************************************************************
  !Initialize some things
  string_to_double = 0.0_DP
  dec_found = .FALSE.
  exp_found = .FALSE.
  cff_bd = ""
  cff_ad = ""
  expv = ""
  string_in = ADJUSTL(string_in)
  nchars = LEN_TRIM(string_in)
  strln = LEN(string_in)
  is_negative = .FALSE.
  IF (string_in(1:1) == "-") is_negative = .TRUE.
  
  ! Make an initial pass through the string to find
  ! if and where the decimal and exponent marker are
  exp_start = -1
  dec_start = -1
  DO ii = 1, nchars
     IF (string_in(ii:ii) == ".") dec_start = ii
     IF (string_in(ii:ii) == "D") exp_start = ii
     IF (string_in(ii:ii) == "d") exp_start = ii
     IF (string_in(ii:ii) == "E") exp_start = ii
     IF (string_in(ii:ii) == "e") exp_start = ii
  END DO
  IF (exp_start > 0) THEN
     exp_found = .TRUE.
  ELSE
     exp_start = nchars + 1
  END IF
  IF (dec_start > 0) THEN
     dec_found = .TRUE.
  ELSE
     dec_start = exp_start
  END IF

  !Based on above, break string into components
  cnt = 1
  DO ii = 1, dec_start - 1
     cff_bd(cnt:cnt) = string_in(ii:ii)
     cnt = cnt + 1
  END DO
  
  IF (dec_found) THEN
     cnt = 1
     DO ii = dec_start + 1, exp_start - 1
        cff_ad(cnt:cnt) = string_in(ii:ii)
        cnt = cnt + 1
     END DO
  ELSE
     cff_ad(1:1) = "0"
  END IF
       
  IF (exp_found) THEN
     cnt = 1
     DO ii = exp_start + 1, nchars
        expv(cnt:cnt) = string_in(ii:ii)
        cnt = cnt + 1
     END DO
  ELSE
     expv(1:1) = "0"
  END IF
  
  !Convert exponent, predecimal components to integers
  icff_bd = string_to_int(cff_bd)
  iexpv= string_to_int(expv)

  !Combine components into real number
  string_to_double = REAL(icff_bd,DP)
  div = 10.0_DP
  DO ii = 1, LEN_TRIM(cff_ad)
     add_num = (IACHAR(cff_ad(ii:ii)) - 48) / div
     IF (is_negative) add_num = -add_num
     string_to_double = string_to_double + add_num
     div = div*10.0_DP
  END DO
  string_to_double = string_to_double*10.0**iexpv

END FUNCTION String_To_Double


!****************************************************************************
FUNCTION String_To_Int(string_in)
!****************************************************************************
  ! This function takes a character string as input out returns the 
  ! equivalent integer.
!****************************************************************************

  IMPLICIT NONE

  LOGICAL :: is_negative
  INTEGER :: string_to_int, ndigits, strln
  INTEGER :: mult, digit, pos, ii
  CHARACTER(*) :: string_in
!****************************************************************************
  !Initialize some things
  string_to_int = 0
  string_in = ADJUSTL(string_in)
  ndigits = LEN_TRIM(string_in)
  strln = LEN(string_in)

  !Find out if the number is negative
  is_negative = .FALSE.
  IF (string_in(1:1) == "-") THEN
     is_negative = .TRUE.
     ndigits = ndigits - 1
  END IF
  IF (string_in(1:1) == "+") ndigits = ndigits - 1

  !Pull of digits starting at the end, multiply by
  !the correct power of ten and add to value
  string_in = ADJUSTR(string_in)
  mult = 1
  DO ii = 1, ndigits
     pos = strln - ii + 1
     digit = IACHAR(string_in(pos:pos)) - 48
     string_to_int = string_to_int + mult*digit
     mult = mult*10
  END DO
     
  !If it's negative, make it so
  IF (is_negative) string_to_int = -string_to_int

END FUNCTION String_To_Int
!****************************************************************************
ELEMENTAL FUNCTION String_To_Logical(string_in)
!****************************************************************************
! This function takes a character string as input out returns the 
! equivalent logical.
!****************************************************************************

  IMPLICIT NONE
 
  CHARACTER(*), INTENT(in) :: string_in
  LOGICAL :: string_to_logical

  String_To_Logical = (string_in(1:1) == 't' .OR. string_in(1:1) == 'T' &
       .OR. string_in(2:2) == 't' .OR. string_in(2:2) == 'T')

END FUNCTION String_To_Logical

SUBROUTINE Check_String(string_in,ierr)
  ! The subroutine checks that the first character of the input string is
  ! an alphabet. Also, it determines if all the characters are alphanumeric
  ! or a dot. 
  
  IMPLICIT NONE

  CHARACTER(*) :: string_in

  INTEGER :: ncharacters, strln, i, ierr
  
  ierr = 0

  string_in = ADJUSTL(string_in)
  ncharacters = LEN_TRIM(string_in)
  strln = LEN(string_in)

  ! Now check for the rest of the characters

  DO i = 1, ncharacters
     IF ( .NOT. ((string_in(i:i) >=  'A' .AND. (string_in(i:i) <= 'Z')) .OR. &
          (string_in(i:i) >= 'a' .AND. (string_in(i:i) <= 'z')) .OR. &
          (string_in(i:i) >= '0' .AND. (string_in(i:i) <= '9')) .OR. &
          (string_in(i:i) == '.' .OR. string_in(i:i) == '_' .OR. string_in(i:i) == '/'))) THEN
        ! character other than letters and digits found
        ierr = 1
        RETURN
     END IF
  END DO
END SUBROUTINE Check_String

END MODULE IO_Utilities
