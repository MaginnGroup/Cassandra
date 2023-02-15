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

module gmxfort_tests

    use, intrinsic :: iso_fortran_env
    use gmxfort_trajectory

    implicit none

    type(Trajectory) :: trj
    real(8), parameter :: PI = 2.0d0*acos(0.0d0)
    character (len=8), parameter :: xtcfile = "test.xtc"
    character (len=8), parameter :: ndxfile = "test.ndx"
    real :: x(3), y(3), z(3), w(3), ans(3), box(3,3), ans_box(3,3), b, c
    integer :: passed = 0, total = 0, a, ans_val

    interface check
        module procedure check_int, check_real, check_array, check_array_2d
    end interface check

contains

    subroutine do_output(total, passed, test_result)

        implicit none
        integer, intent(inout) :: total, passed
        logical, intent(in) :: test_result
        integer :: I
        character (len=6) :: passfail(0:1) = ["FAILED", "PASSED"]
        
        I = merge(1, 0, test_result)
        total = total + 1
        passed = passed + I

        write(output_unit, '(a,i2,a,a)') "TEST ", total, ": ", passfail(I)

    end subroutine do_output

    subroutine check_int(x, y, passed, total)

        implicit none

        integer, intent(inout) :: total, passed
        integer, intent(in) :: x, y

        call do_output(total, passed, x .eq. y)

    end subroutine check_int

    subroutine check_real(x, y, passed, total)

        implicit none

        integer, intent(inout) :: total, passed
        real, intent(in) :: x, y
        real :: tol = 1e-4

        call do_output(total, passed, abs(x-y) .le. tol)

    end subroutine check_real

    subroutine check_array(x, y, passed, total)

        implicit none

        integer, intent(inout) :: total, passed
        real, intent(in) :: x(:), y(:)
        real :: tol = 1e-4

        call do_output(total, passed, all(abs(x - y) .le. tol))

    end subroutine check_array

    subroutine check_array_2d(x, y, passed, total)

        implicit none

        integer, intent(inout) :: total, passed
        real, intent(in) :: x(:,:), y(:,:)
        real :: tol = 1e-6

        call do_output(total, passed, all(abs(x - y) .le. tol))

    end subroutine check_array_2d

    subroutine finished_tests(passed, total)

        implicit none

        integer, intent(inout) :: total, passed
        
        write(output_unit,*)
        write(output_unit,'(a,i0,a,i0,a)') "Passed ", passed, " out of ", total, " tests"
        write(output_unit,*)

        if (passed .ne. total) then
            write(output_unit, '(a)') "WARNING: Some tests failed!"
            stop 1
        end if

    end subroutine finished_tests

end module gmxfort_tests

