!
! This file is part of libgmxfort
! https://github.com/wesbarnett/libgmxfort
!
! Copyright (c) 2016-2019 James W. Barnett
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

module gmxfort_index

    use gmxfort_common

    implicit none
    private

    type ndxgroups
        integer, allocatable :: LOC(:)
        integer :: NUMATOMS
        character (len=:), allocatable :: title
    end type ndxgroups

    type, public :: IndexFile
        type (ndxgroups), allocatable :: group(:)
        logical :: group_warning = .true.
    contains
        procedure :: read => indexfile_read
        procedure :: get => indexfile_get
        procedure :: get_natoms => indexfile_get
    end type IndexFile
 
contains

    subroutine indexfile_read(this, filename, N)

        implicit none
        class(IndexFile), intent(inout) :: this
        character (len=*), intent(in) :: filename
        character (len=2048) :: line, NCOLS_string, fmt_string
        integer :: INDEX_FILE_UNIT, IO_STATUS, NGRPS, I, J, NCOLS
        integer, allocatable :: INDICES_TMP(:), TITLE_LOC(:), num_array(:)
        logical :: ex
        integer, intent(in) :: N

        ! Does the file exist?
        inquire(file=trim(filename), exist=ex)
        if (ex .eqv. .false.) call error_stop_program(trim(filename)//" does not exist.")

        ! Is in index file?
        open(newunit=INDEX_FILE_UNIT, file=trim(filename), status="old")
        read(INDEX_FILE_UNIT, '(a)', iostat=IO_STATUS) line
        if (index(line, "[") .eq. 0) call error_stop_program(trim(filename)//" is not a valid index file.")

        ! How many groups are in it?
        rewind INDEX_FILE_UNIT
        IO_STATUS = 0
        NGRPS = 0
        do while (IO_STATUS .eq. 0)
            read(INDEX_FILE_UNIT, '(a)', iostat=IO_STATUS) line
            if (IO_STATUS .ne. 0) goto 100
            if (index(line, "[") .ne. 0) NGRPS = NGRPS + 1
        end do

100     continue

        if (allocated(this%group)) deallocate(this%group)
        allocate(this%group(NGRPS), TITLE_LOC(NGRPS+1)) ! Add one to include end of file

        ! Now find the title locations and save their names
        rewind INDEX_FILE_UNIT
        I = 1
        J = 1
        IO_STATUS = 0
        do while (IO_STATUS .eq. 0)

            read(INDEX_FILE_UNIT, '(a)', iostat=IO_STATUS) line
            if (IO_STATUS .ne. 0) goto 200
            if (index(line, "[") .ne. 0) then
                this%group(I)%title = trim(line(index(line, "[")+2:index(line, "]")-2))
                TITLE_LOC(I) = J
                I = I + 1
            end if
            J = J + 1

        end do

200     continue

        if (this%group(1)%title .ne. "System") call error_stop_program("Index file does not have 'System' group as first group.")

        ! Index files can have a varying number of columns. This attempts to
        ! detect the correct number by reading in the second line of the file,
        ! which should be a list of indices for the "System" group.
        NCOLS = 100
        allocate(num_array(NCOLS))
        IO_STATUS = 5000
        do while (IO_STATUS .ne. 0)
            NCOLS = NCOLS - 1
            write(NCOLS_string, '(i0)') NCOLS
            write(fmt_string, '(a)') '('//trim(NCOLS_string)//'i0)'
            rewind INDEX_FILE_UNIT
            read(INDEX_FILE_UNIT, '(a)', iostat=IO_STATUS) line
            read(INDEX_FILE_UNIT, '(a)', iostat=IO_STATUS) line
            read(line, *, iostat=IO_STATUS) num_array(1:NCOLS)
        end do

        TITLE_LOC(I) = J ! End of file location

        ! Now finally get all of the indices for each group
        ! Allocate for total number of atoms in system, since that is the maximum
        allocate(INDICES_TMP(N))
        do I = 1, NGRPS

            ! Initial guess only how many items are in the group
            ! Add 1, bc loop subtracts 1 at the beginning
            this%group(I)%NUMATOMS = (TITLE_LOC(I+1)-TITLE_LOC(I)-1)*NCOLS + 1

            if (N < this%group(I)%NUMATOMS) this%group(I)%NUMATOMS = N + 1
            IO_STATUS = 5000

            do while (IO_STATUS .ne. 0)

                ! Our guess was too large if we made it back here, go to the beginning and reduce our guess by 1, try again
                rewind INDEX_FILE_UNIT
                this%group(I)%NUMATOMS = this%group(I)%NUMATOMS - 1
                if (this%group(I)%NUMATOMS .le. 0) then 
                    this%group(I)%NUMATOMS = 0
                    goto 300
                end if

                ! Read all the way to the group
                do J = 1, TITLE_LOC(I); read(INDEX_FILE_UNIT, '(a)', iostat=IO_STATUS) line; end do

                ! Attempt to read into array
                read(INDEX_FILE_UNIT, *, iostat=IO_STATUS) INDICES_TMP(1:this%group(I)%NUMATOMS)

            end do

            ! Specifying array bounds for array to be allocated is not required for F2008 but is required for F2003
            allocate(this%group(I)%LOC(1:this%group(I)%NUMATOMS), source=INDICES_TMP(1:this%group(I)%NUMATOMS))

300         cycle

        end do
        deallocate(INDICES_TMP)

        do I = 1, NGRPS-1
            do J = I+1, NGRPS
                if (I .ne. J) then
                    if (this%group(I)%title .eq. this%group(J)%title) then
                        write(error_unit,*)
                        write(error_unit,'(a, a, a)') "LIBGMXFORT WARNING: Index group ", this%group(I)%title, &
                            " was specified more than once in index file."
                        write(error_unit,*)
                    end if
                end if
            end do
        end do

        ! Error checking to see if index file goes with the trajectory file

        ! If the number of atoms is not the same in the System group (first group) and xtc file
        if (this%group(1)%numatoms .ne. N .or. this%group(1)%loc(this%group(1)%numatoms) .ne. N) then
            call error_stop_program("Index file does not match xtc file.")
        end if

        do i = 1, NGRPS

            ! If number of atoms in index group is larger than number of atoms in xtc file
            if (this%group(i)%numatoms .gt. N) call error_stop_program("Index file does not match xtc file.")

            ! If a location number is greater than number of atoms in xtc file
            do j = 1, this%group(i)%numatoms
                if (this%group(i)%loc(j) .gt. N) call error_stop_program("Index file does not match xtc file.")
            end do

        end do


        close(INDEX_FILE_UNIT)
        
    end subroutine indexfile_read

    ! Gets the number of atoms in a group. If an atom is specified, integer returns the overall index for that atom.
    function indexfile_get(this, group_name, I)

        implicit none
        integer :: indexfile_get
        class(IndexFile), intent(inout) :: this
        character (len=*), intent(in) :: group_name
        integer, intent(in), optional :: I
        integer :: J

        do J = 1, size(this%group)

            if (trim(this%group(J)%title) .eq. trim(group_name)) then

                if (present(I)) then
                        indexfile_get = this%group(J)%LOC(I)
                else
                        indexfile_get = this%group(J)%NUMATOMS
                end if
                return

            end if

        end do

        ! If user asked to get the number of atoms in an index group, and that index group is not 
        ! in the index file, just return 0. If the user specified an atom number, then throw an error,
        ! since the overall index cannot be returned in that case
        if (.not. present(I)) then
            if (this%group_warning) then
                write(error_unit, '(a)') "LIBGMXFORT WARNING: No atoms found in index group '"//trim(group_name)//"'."
                write(error_unit, '(a)') "This warning will not appear again for any other index groups."
                this%group_warning = .false.
            end if
            indexfile_get = 0
        else
            indexfile_get = -1
            write(error_unit, '(a)') "LIBGMXFORT ERROR: "//trim(group_name)//" is not in index file. The groups available are:"
            do J = 1, size(this%group)
                write(error_unit,'(a10,a,i0,a)') this%group(J)%title, " (", this%group(J)%NUMATOMS, ")"
            end do
            stop 1
        end if

    end function indexfile_get

end module gmxfort_index
