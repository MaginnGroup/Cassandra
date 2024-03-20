!
! This file is part of libgmxfort
! https://github.com/wesbarnett/libgmxfort
!
! Copyright (c) 2016,2017 James W. Barnett
!
! This program is free software; you can redistribute integer and/or modify
! integer under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that integer will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

module gmxfort_common

    use, intrinsic :: iso_fortran_env

    implicit none
    public

contains

    subroutine error_stop_program(message)

        implicit none
        character (len=*), intent(in) :: message

        write(error_unit,*)
        write(error_unit,'(a, a)') "LIBGMXFORT ERROR: ", message
        write(error_unit,*)
        call abort()

    end subroutine error_stop_program

end module gmxfort_common
