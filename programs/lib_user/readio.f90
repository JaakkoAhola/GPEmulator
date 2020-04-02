module readIO

PUBLIC :: readFileDimensions, readArray

contains
    subroutine readFileDimensions( fileName, rows, columns, inputDelimiter, usePrinting)
        implicit none

        character(len=1000), intent(in) :: fileName
        integer, intent(out) :: rows, columns
        character, intent(in), optional  :: inputDelimiter
        logical, intent(in), optional    :: usePrinting

        character :: delimiter = " "

        integer :: ioStatus
        integer, parameter :: read_unit = 98
        character(len=1000) :: line
        character(len=1000):: tmpString
        integer :: ind = 1
        logical :: usePrintingLocal = .FALSE.

        if ( present(inputDelimiter)) then
            delimiter = inputDelimiter
        else
            delimiter = " "
        end if

        if ( present(usePrinting)) then
            usePrintingLocal = usePrinting
        else
            usePrintingLocal = .FALSE.
        end if

        if ( usePrintingLocal ) write(*,*) "readFileDimensions fileName: ", trim(fileName)
        open( unit=read_unit, file = trim(fileName), iostat=ioStatus )
        if ( ioStatus /= 0 ) stop "readFileDimensions: Error opening file "

        rows = 0
        do
            read(read_unit, '(A)', iostat=ioStatus) line
            if (ioStatus /= 0) exit
            rows = rows + 1
        end do

        if ( usePrintingLocal ) print*, trim(line)
        ! analyze how many columns
        ! 0. remove initial blanks if any
        tmpString =trim (adjustl(line) )

        ! 1. count the number substrings separated by delimiter
        columns = count( [ (tmpString(ind:ind), ind=1, len_trim(tmpString)) ] == delimiter) + 1

        if ( usePrintingLocal ) then
            print*, "rows", rows
            print*, "columns", columns
            print*, "end of subroutine readFileDimensions"
        end if

        close(read_unit)

    end subroutine readFileDimensions

    subroutine readArray(fileName, array, inputDelimiter, usePrinting)

        implicit none
        integer, parameter ::  dp = selected_real_kind(14,100)
        character(len=1000), intent(in) :: fileName
        real(dp), dimension(:,:), intent(out), allocatable :: array
        character, intent(in), optional  :: inputDelimiter
        logical, intent(in), optional    :: usePrinting

        character :: delimiter = " "

        integer :: ioStatus
        integer, parameter :: read_unit = 99

        integer :: rowInd, columnInd, rows, columns
        logical :: usePrintingLocal = .FALSE.

        if ( present(inputDelimiter)) then
            delimiter = inputDelimiter
        end if

        if ( present(usePrinting)) then
            usePrintingLocal = usePrinting
        else
            usePrintingLocal = .FALSE.
        end if

        call readFileDimensions( fileName, rows, columns, delimiter, usePrintingLocal)

        if (usePrintingLocal) then
            print*, "readArray File contains ", rows, "rows"
            print*, "readArray File contains ", columns, "columns"
        end if

        allocate(array(rows,columns))

        open( unit=read_unit, file= trim(fileName), iostat=ioStatus )
        if ( ioStatus /= 0 ) stop "Error opening file"

        do rowInd = 1,rows
            read(read_unit,*) (array(rowInd,columnInd),columnInd=1,columns)
        end do

        if (usePrintingLocal) then
            print*, "array(1,:)", array(1,:)
            print*, "array(2,:)", array(2,:)
        end if

        close(read_unit)
    end subroutine readArray

end module readIO
