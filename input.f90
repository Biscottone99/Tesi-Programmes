  program crea_input_files
    implicit none
    integer :: i, num_files
    real*8 :: min_degree, max_degree, min_radian, max_radian
    character(len=30) :: filename

    ! Numero di file da creare
    num_files = 201
    open(1,file='num.dat')
    write(1,*) num_files
    close(1)
    ! Range dei valori in gradi
    min_degree = 0.000000
    max_degree = 360.0
    
    ! Conversione dei valori in radianti
    min_radian = deg2rad(min_degree)
    max_radian = deg2rad(max_degree)

    ! Creazione dei valori in radianti e scrittura dei file
    do i = 1, num_files
        ! Genera il nome del file
        write(filename, '(a, i0, a)') 'input_', i, '.dat'

        ! Apri il file in scrittura
        open(unit=10, file=trim(filename), status='replace', action='write')

        ! Scrivi il valore in radianti nel file
        write(10, *) min_radian + (max_radian - min_radian) * real(i-1) / real(num_files - 1)

        ! Chiudi il file
        close(unit=10)
    end do

contains

    ! Funzione per convertire gradi in radianti
    real function deg2rad(degree)
        real*8, intent(in) :: degree
        deg2rad = degree * (3.14159265358979323846 / 180.0)
    end function deg2rad

end program crea_input_files
