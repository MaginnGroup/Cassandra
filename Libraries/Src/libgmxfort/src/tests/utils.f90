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

program utils_test 

    use gmxfort_utils
    use gmxfort_tests

    ! TEST 1
    x = [5.5, 5.5, 3.5]
    box = reshape((/5.0, 0.0, 0.0, &
                    0.0, 5.0, 0.0, &
                    2.5, 2.5, 3.5/), shape(box))
    x = pbc(dble(x), dble(box))
    ans = [-2.000, -2.000, 0.000]
    call check(x, ans, passed, total) 

    ! TEST 2
    x = [5.5, 5.5, 3.5]
    b = distance(dble(x), dble(ans), dble(box))
    call check(b, 0.0, passed, total)

    ! TEST 3
    x = [5.5, 5.5, 3.5]
    y = [3.6, 4.7, 5.0]
    box = reshape((/3.5, 0.0, 0.0, &
                    0.0, 4.5, 0.0, &
                    0.0, 0.0, 4.0/), shape(box))
    b = distance(dble(x), dble(y), dble(box))
    call check(b, 2.33452, passed, total)

    ! TEST 4
    b = magnitude(dble(x))
    call check(b, 8.52936, passed, total)

    ! TEST 5
    x = [0.0, 0.0, 0.0]
    y = [0.0, 1.0, 0.0]
    z = [1.0, 1.0, 0.0]
    b = bond_angle(dble(x), dble(y), dble(z), dble(box))
    call check(b, real(PI/2.0d0), passed, total)

    ! TEST 6
    w = [1.0, 1.0, 1.0]
    b = dihedral_angle(dble(x), dble(y), dble(z), dble(w), dble(box))
    call check(b, real(PI/2.0d0), passed, total)

    ! TEST 7
    w = [1.0, 1.0, -1.0]
    b = dihedral_angle(dble(x), dble(y), dble(z), dble(w), dble(box))
    call check(b, real(-PI/2.0d0), passed, total)

    ! TEST 8
    call trj%open(xtcfile, ndxfile)
    a = trj%read_next(50, "OW")
    a = trj%read_next(200, "OW")
    x = trj%x(50, 100)
    ans = [0.115, 1.048, 3.222]
    call check(x, ans, passed, total)
    call trj%close()

    call finished_tests(passed, total)

end program utils_test
