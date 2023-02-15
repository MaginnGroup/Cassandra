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

program read_next_test 

    use gmxfort_tests

    call trj%open(xtcfile, ndxfile)

    ! TEST 1
    a = trj%read_next(100)
    x = trj%x(100, 100)
    ans = [1.455, 0.374, 0.358]
    call check(x, ans, passed, total)

    ! TEST 2
    ans_val = 100
    call check(a, ans_val, passed, total) 

    ! TEST 3
    a = trj%read_next(200)
    ans_val = 1
    call check(a, ans_val, passed, total) 

    ! TEST 4
    a = trj%read_next()
    ans_val = 0
    call check(a, ans_val, passed, total) 
    call trj%close()

    call finished_tests(passed, total)

end program read_next_test
