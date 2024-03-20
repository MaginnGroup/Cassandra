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

program xtcfile_test

    use gmxfort_tests

    call trj%read(xtcfile)

    ! TEST 1
    x = trj%x(1, 1)
    ans = [1.749, 0.752, 0.220]
    call check(x, ans, passed, total)

    ! TEST 2
    x = trj%x(50, 50)
    ans = [0.359, 1.999, 1.816]
    call check(x, ans, passed, total)

    ! TEST 3
    x = trj%x(trj%NFRAMES, trj%natoms())
    ans = [4.060, 0.155, 0.262]
    call check(x, ans, passed, total)

    call finished_tests(passed, total)

end program xtcfile_test
