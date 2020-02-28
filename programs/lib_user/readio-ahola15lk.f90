module readIO

PUBLIC :: readFileDimensions, readArray

contains
    subroutine readFileDimensions( fileName, rows, columns, inputDelimiter)
        implicit none
        
        character(len=1000), intent(in) :: fileName
        integer, intent(out) :: rows, columns
        character, intent(in), optional  :: inputDelimiter

        character :: delimiter = " "

        integer :: ioStatus
        integer, parameter :: read_unit = 98
        character(len=1000) :: line
        character(len=1000) :: tmpString
        integer :: ind
        
        if ( present(inputDelimiter)) delimiter = inputDelimiter
        
        open( unit=read_unit, file= trim(fileName), iostat=ioStatus )
        if ( ioStatus /= 0 ) stop "Error opening file"
        
        rows = 0
        do
            read(read_unit, '(A)', iostat=ioStatus) line
            if (ioStatus /= 0) exit
            rows = rows + 1
        end do
        
        
        ! analyze how many columns
        ! 0. remove initial blanks if any
        tmpString =trim (adjustl(line) )

        ! 1. count the number substrings separated by delimiter
        columns = count( [ (tmpString(ind:ind), ind=1, len_trim(tmpString)) ] == delimiter) + 1
        
    
        close(read_unit)
        
    end subroutine readFileDimensions

    subroutine readArray(fileName, array, inputDelimiter)

        implicit none

        character(len=1000), intent(in) :: fileName
        real, dimension(:,:), intent(out), allocatable :: array
        character, intent(in), optional  :: inputDelimiter

        character :: delimiter = " "

        integer :: ioStatus
        integer, parameter :: read_unit = 99

        integer :: rows, columns

        if ( present(inputDelimiter)) then
            delimiter = inputDelimiter
        end if

        open( unit=read_unit, file= trim(fileName), iostat=ioStatus )
        if ( ioStatus /= 0 ) stop "Error opening file"

        call readFileDimensions( fileName, rows, columns, delimiter)

        print*, "File contains ", rows, "rows"
        print*, "File contains ", columns, "columns"

        allocate(array(rows,columns))

        rewind(read_unit)

        read(read_unit,*) array

        close(read_unit)
    end subroutine readArray

end module readIO
