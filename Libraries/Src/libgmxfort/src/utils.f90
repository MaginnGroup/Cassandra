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

module gmxfort_utils

    implicit none
    public
 
contains

    function pbc(a, box)

        implicit none
        real(8), intent(in) :: a(3), box(3,3)
        real(8) :: pbc(3)
        integer :: I

        pbc = a

        do I = 3, 1, -1
            pbc(1:I) = pbc(1:I) - box(1:I,I) * nint(pbc(I) / box(I,I))
        end do

    end function pbc


    function cross(a, b)

        implicit none
        real(8) :: cross(3)
        real(8), intent(in), dimension(3) :: a, b

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)

    end function cross


    function distance2(a, b, box)

        implicit none
        real(8) :: distance2
        real(8), intent(in), dimension(3) :: a, b
        real(8) :: c(3)
        real(8), intent(in) :: box(3,3)

        c = pbc(a - b, box)
        distance2 = dot_product(c, c)

    end function distance2


    function distance(a, b, box)

        implicit none
        real(8) :: distance
        real(8), intent(in), dimension(3) :: a, b
        real(8), intent(in) :: box(3,3)
        distance = dsqrt(distance2(a, b, box))

    end function distance

    function bond_vector(a, b, box)

        implicit none
        real(8) :: bond_vector(3)
        real(8), intent(in), dimension(3) :: a, b
        real(8), intent(in) :: box(3,3)

        bond_vector = pbc(a-b, box)

    end function bond_vector


    function magnitude(a)

        implicit none
        real(8) :: magnitude
        real(8), intent(in) :: a(3)

        magnitude = dsqrt(dot_product(a, a))

    end function magnitude


    function bond_angle(a, b, c, box)

        implicit none
        real(8) :: bond_angle
        real(8), intent(in), dimension(3) :: a, b, c
        real(8), intent(in) :: box(3,3)
        real(8), dimension(3) :: bond1, bond2

        bond1 = bond_vector(b, a, box)
        bond2 = bond_vector(b, c, box)

        bond_angle = dacos(dot_product(bond1, bond2)/(magnitude(bond1)*magnitude(bond2)))

    end function bond_angle


    function dihedral_angle(i, j, k, l, box)

        implicit none
        real(8) :: dihedral_angle
        real(8), intent(in), dimension(3) :: i, j, k, l
        real(8), intent(in), dimension(3,3) :: box
        real(8) :: A_mag, B_mag, G_mag
        real(8), dimension(3) :: H, G, F, A, B, cross_BA

        H = bond_vector(k, l, box)
        G = bond_vector(k, j, box)
        F = bond_vector(j, i, box)
        A = cross(F, G)
        B = cross(H, G)
        cross_BA = cross(B, A)
        A_mag = magnitude(A)
        B_mag = magnitude(B)
        G_mag = magnitude(G)

        dihedral_angle = atan2(dot_product(cross_BA, G) / (A_mag*B_mag*G_mag), dot_product(A, B) / (A_mag*B_mag))

    end function dihedral_angle

end module gmxfort_utils
