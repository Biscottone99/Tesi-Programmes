module module
implicit none
contains
  function cross_product(v1, v2) result(result_vector)
    implicit none
    complex*16, dimension(3), intent(in) :: v1, v2
    complex*16, dimension(3) :: result_vector

    result_vector(1) = v1(2) * v2(3) - v1(3) * v2(2)
    result_vector(2) = v1(3) * v2(1) - v1(1) * v2(3)
    result_vector(3) = v1(1) * v2(2) - v1(2) * v2(1)
  end function cross_product

  subroutine copy_upper_to_lower(matrix, n)
    implicit none
    integer, intent(in) :: n
    real*8, dimension(n, n), intent(inout) :: matrix

    integer :: i, j

    do i = 1, n
       do j = i+1, n
          matrix(j, i) = matrix(i, j)
       end do
    end do
  end subroutine copy_upper_to_lower

  subroutine copy_upper_to_lower_cplx(matrix, n, l)
    implicit none
    integer, intent(in) :: n, l
    complex*16, dimension(l,n, n), intent(inout) :: matrix

    integer :: i, j,m

    do m=1,l
       do i = 1, n
          do j = i+1, n
             matrix(m,j, i) = matrix(m,i, j)
          end do
       end do
    enddo
  end subroutine copy_upper_to_lower_cplx

  subroutine copia_comp_conj(matrice, n)
    implicit none
    integer, intent(in) :: n
    complex*16, intent(inout) :: matrice(n,n)
    integer :: i, j
    ! Copia il complesso coniugato dal triangolo superiore al triangolo inferiore
    do i = 1, n
       do j = i + 1, n
          matrice(j, i) = dconjg(matrice(i, j))
       end do
    end do
  end subroutine copia_comp_conj

  recursive function binary_search(arr, target, low, high) result(index)
    implicit none
    integer, intent(in) :: arr(:), target, low, high
    integer :: index, mid

    if (low > high) then
       index = 0  ! Target not found
    else
       mid = (low + high) / 2

       if (arr(mid) == target) then
          index = mid  ! Target found
       else if (arr(mid) < target) then
          !$omp parallel default(none) shared(arr, target, mid, high, index)
          !$omp single
          index = binary_search(arr, target, mid + 1, high)
          !$omp end single
          !$omp end parallel
       else
          !$omp parallel default(none) shared(arr, target, low, mid, index)
          !$omp single
          index = binary_search(arr, target, low, mid - 1)
          !$omp end single
          !$omp end parallel
       end if
    end if
  end function binary_search

  !generates a DOUBLE PRECISION one-electron operator in second quantization
  !given a real space basis and one-electron operator in first quantization
  !nso    number of spin-orbitals
  !dim    real-space dimension
  !op_tb  operator in the first quantization (tight binding basis)
  !op     operator in the real space basis
  !basis  1D array containing the integers that describe the real space basis in bit representaion
  !NB: the same spin-orbital ordering used to describe the real space configurations
  !    must be used for the first quantization basis
  subroutine sq_oe_op_real(nso,dim,op_tb,op,basis)
    implicit none
    integer :: iso,jso,conta,i,j,istate,jstate,step,a
    integer, intent(in) :: dim,nso,basis(dim)
    double precision, intent(in) :: op_tb(nso,nso)
    double precision, intent(out) :: op(dim,dim)
    double precision :: phase

    op = 0.d0
    !$omp parallel do default(none) private(iso,jso,conta,i,j,istate,jstate,step,a,phase) shared(dim,nso,basis,op_tb,op)
    do j = 1,dim !col index
       do iso = 0,nso-1 ! creation op index
          do jso = 0,nso-1 ! annihilation op index
             jstate = basis(j)
             if (btest(jstate,jso)) then
                istate = ibclr(jstate,jso)
                if (.not.btest(istate,iso)) then
                   istate = ibset(istate,iso)
                   i = find_value(basis, istate, dim) !row index
                   if (i/=0) then
                      !determine the phase
                      !get direction from iso to jso
                      if (jso>iso) step = -1
                      if (iso>jso) step = 1

                      if (iso==jso) then
                         phase = 1.d0
                      else
                         conta = 0
                         do a = jso+step, iso-step, step
                            if (btest(istate,a)) conta = conta + 1
                         end do

                         if (mod(conta,2)==0) then
                            phase = 1.d0
                         else
                            phase = -1.d0
                         end if
                      end if

                      op(i,j) = op(i,j) + phase * op_tb(iso+1,jso+1)
                   end if
                end if
             end if

          end do
       end do
    end do
    !$omp end parallel do
  end subroutine sq_oe_op_real

  !generates a DOUBLE COMPLEX one-electron operator in second quantization
  !given a real space basis and one-electron operator in first quantization
  !nso    number of spin-orbitals
  !dim    real-space dimension
  !op_tb  operator in the first quantization (tight binding basis)
  !op     operator in the real space basis
  !basis  1D array containing the integers that describe the real space basis in bit representaion
  !NB: the same spin-orbital ordering used to describe the real space configurations
  !    must be used for the first quantization basis
  subroutine sq_oe_op_compl(nso,dim,op_tb,op,basis)
    implicit none
    integer iso,jso,conta,i,j,istate,jstate,step,a
    integer, intent (in) :: dim,nso,basis(dim)
    double complex, intent (in) :: op_tb(nso,nso)
    double complex, intent (out) :: op(dim,dim)
    double precision phase

    op = 0.d0
    !$omp parallel do default(none) private(iso,jso,conta,i,j,istate,jstate,step,a,phase) shared(dim,nso,basis,op_tb,op)
    do j = 1,dim !col index
       do iso = 0,nso-1 ! creation op index
          do jso = 0,nso-1 ! annihilation op index
             jstate = basis(j)
             if (btest(jstate,jso)) then
                istate = ibclr(jstate,jso)
                if (.not.btest(istate,iso)) then

                   istate = ibset(istate,iso)
                   i = find_value(basis, istate,dim)  !row index
                  
                   if (i/=0) then
                      !determine the phase
                      !get direction from iso to jso
                      if (jso>iso) step = -1
                      if (iso>jso) step = 1

                      if (iso==jso) then
                         phase = 1.d0
                         goto 1000
                      end if

                      conta = 0
                      do a = jso+step, iso-step, step
                         if (btest(istate,a)) conta = conta + 1
                      end do

                      if (conta/2*2==conta) then
                         phase = 1.d0
                      else
                         phase = -1.d0
                      end if

1000                  continue

                      !write(*,*) phase,i,j,(btest(istate,a),a=0,nso-1), 0,0, (btest(jstate,a),a=0,nso-1)
                      op(i,j) = op(i,j) + phase * op_tb(iso+1,jso+1)
                   end if
                end if
             end if

          end do
       end do
       !write(*,*) ''

    end do
    !$omp end parallel do
  end subroutine sq_oe_op_compl

  subroutine rotate_real_2x2(dim2,coupling,coup,ham)
    implicit none
    integer::i,l,j,k
    integer, intent(in)::dim2
    real*8,intent(out)::coupling(dim2,dim2)
    real*8,intent(in)::coup(dim2,dim2)
    complex*16,intent(in)::ham(dim2,dim2)

    !$omp parallel do default(none) private(i, l, j, k) shared(coupling, coup, ham, dim2)
    do i=1,dim2 !a
       do l=1,dim2 !b
          do j=1,dim2 !alfa 
             do k=1,dim2 !beta
                coupling(i,l) = coupling(i,l) + dconjg(ham(j,i))*ham(k,l)*coup(j,k)
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do
  end subroutine rotate_real_2x2

  subroutine rotate_cplx_2x2(dim2, coupling, coup, ham)
    implicit none
    integer :: i, l, j, k
    integer, intent(in) :: dim2
    complex*16, intent(out) :: coupling(dim2, dim2)
    complex*16, intent(in) :: coup(dim2, dim2)
    complex*16, intent(in) :: ham(dim2, dim2)

    !$omp parallel do default(none) private(i, l, j, k) shared(coupling, coup, ham, dim2)
    do i = 1, dim2
       do l = 1, dim2
          do j = 1, dim2
             do k = 1, dim2
                coupling(i, l) = coupling(i, l) + dconjg(ham(j, i)) * ham(k, l) * coup(j, k)
             end do
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine rotate_cplx_2x2

  !ruota una matrice reale non quadrata sullo spazio degli autostati
  !out is the resulting matrix
  !in is the matrix being rotated
  !rot is the complex matrix with eigenvectors
  !col is the number of columns

  subroutine rot_diag(dim2, out, in, col, rot)
    implicit none
    integer :: i, j, k
    integer, intent(in) :: dim2, col
    complex*16, intent(in) :: rot(dim2, dim2)
    real*8, intent(in) :: in(dim2, col)
    real*8, intent(out) :: out(dim2, col)

    !$omp parallel do default(none) private(i, j, k) shared(dim2, col, rot, in, out)
    do i = 1, dim2
       do j = 1, dim2
          do k = 1, col
             out(i, k) = out(i, k) + dconjg(rot(j, i)) * rot(j, i) * in(j, k)
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine rot_diag

  subroutine copia_matrice_complessa_multi(matrice_origine, matrice_destinazione, righe, colonne, riga_fissa, colonna_fissa, dim3, dim4)
    integer, intent(in) :: righe, colonne, riga_fissa, colonna_fissa, dim3, dim4
    complex*16, intent(in) :: matrice_origine(righe, colonne, dim3, dim4)
    complex*16, intent(out) :: matrice_destinazione(righe, colonne)
    integer :: i, j

    do i = 1, righe
       do j = 1, colonne
          matrice_destinazione(i, j) = matrice_origine(i, j, riga_fissa, colonna_fissa)
       end do
    end do
  end subroutine copia_matrice_complessa_multi

  subroutine copia_matrice_complessa_cont(matrice_origine, matrice_destinazione, righe, colonne, riga_fissa, colonna_fissa, dim3,dim4)
    integer, intent(in) :: righe, colonne, riga_fissa, colonna_fissa, dim3, dim4
    complex*16, intent(in) :: matrice_origine(righe, colonne)
    complex*16, intent(out) :: matrice_destinazione(righe, colonne,dim3,dim4)
    integer :: i, j

    do i = 1, righe
       do j = 1, colonne
          matrice_destinazione(i, j, riga_fissa, colonna_fissa) = matrice_origine(i, j)
       end do
    end do
  end subroutine copia_matrice_complessa_cont

  subroutine check_hermitian(matrix, n, is_hermitian)
    implicit none
    integer, intent(in) :: n
    complex*16, intent(in) :: matrix(n, n)
    logical, intent(out) :: is_hermitian

    integer :: i, j
    complex*16, allocatable:: conj_transpose_matrix(:, :)
    allocate(conj_transpose_matrix(n, n))
    conj_transpose_matrix=0

    ! Calcola la trasposta coniugata della matrice
    do i = 1, n
       do j = 1, n
          conj_transpose_matrix(i, j) = dconjg(matrix(j, i))
       end do
    end do

    ! CHeck if the matrix is Hermitian
    is_hermitian = all(zabs(matrix - conj_transpose_matrix).le.1d-14)

  end subroutine check_hermitian

  subroutine ohno(dim2, nsiti, nuclei, vecconfig, u, nz, pot)
    implicit none
    integer, intent(in) :: dim2, nsiti
    integer, intent(in) :: vecconfig(dim2), nz(nsiti)
    real*8, intent(in) :: nuclei(nsiti, 3), u(nsiti)
    real*8, intent(out) :: pot(dim2)
    integer :: n, i, j, p, sito, occupazioni(nsiti), a, b
    real*8 :: PPP, r(nsiti, nsiti), dx, dy, dz
    logical :: bool, bool1
    pot = 0.0d0
    !$omp parallel do default(none) private(n, sito, i, j, p, occupazioni, a, b, PPP, dx, dy, dz, bool, bool1) shared(dim2, nsiti, nuclei, vecconfig, u, nz, r, pot)
    do n = 1, dim2
       sito = 0
       do i = 0, 2 * nsiti - 2, 2
          sito = (i + 2) / 2
          bool = btest(vecconfig(n), i)
          bool1 = btest(vecconfig(n), i + 1)
          a = merge(0, 1, bool)
          b = merge(0, 1, bool1)
          occupazioni(sito) = a + b
       enddo
       r = 0.0d0
       do i = 1, nsiti
          do j = i + 1, nsiti
             dx = nuclei(i, 1) - nuclei(j, 1)
             dy = nuclei(i, 2) - nuclei(j, 2)
             dz = nuclei(i, 3) - nuclei(j, 3)

             r(i, j) = dsqrt(dx**2 + dy**2 + dz**2)
             r(j, i) = r(i, j)
          enddo
       enddo

       PPP = 0.0d0
       do i = 1, nsiti
          do p = 1, nsiti
             if (i /= p) PPP = PPP + (14.397d0) / dsqrt(r(i, p)**2 + (28.794d0 / (u(i) + u(p)))**2) * (nz(i) - occupazioni(i)) * (nz(p) - occupazioni(p))
             if ((i == p) .and. (occupazioni(i) == 2)) PPP = PPP + 2.0d0 * u(i)
          enddo
       enddo
       pot(n) = 0.5d0 * PPP
    enddo
    !$omp end parallel do
  end subroutine ohno
  
  subroutine site_energy_u(nso, dim2, esite, u, vecconfig, energy)
    implicit none
    integer, intent(in) :: nso, dim2
    integer, intent(in) :: vecconfig(dim2)
    real*8, intent(in) :: esite(nso/2), u(nso/2)
    real*8, intent(out) :: energy(dim2)
    integer :: n, i, sito
    logical :: bool, bool1

    energy = 0.0d0

    !$omp parallel do default(none) private(n, i, sito, bool, bool1) shared(nso, dim2, vecconfig, esite, u, energy)
    do n = 1, dim2
       do i = 0, nso - 1
          sito = (i + 2) / 2
          if (btest(vecconfig(n), i)) then
             energy(n) = energy(n) + esite(sito)
          endif
       enddo

       do i = 0, nso - 2, 2
          sito = (i + 2) / 2
          bool = btest(vecconfig(n), i)
          bool1 = btest(vecconfig(n), i + 1)
          if (bool .and. bool1) then
             energy(n) = energy(n) + u(sito)
          endif
       enddo
    enddo
    !$omp end parallel do
  end subroutine site_energy_u
 
  subroutine site_energy(nso,dim2,esite,vecconfig,energy)
    implicit none
    integer::n,i, sito
    real*8, intent(in)::esite(nso/2)
    integer,intent(in)::dim2, vecconfig(dim2), nso
    real*8, intent(out)::energy(dim2)
    logical::bool, bool1


    energy(n)=0
    do n=1,dim2
       do i=0, nso-1
          sito=(i+2)/2
          if(btest(vecconfig(n),i))energy(n)=energy(n)+esite(sito)
       enddo
    enddo
  end subroutine site_energy

  subroutine momentum_so (nsiti, nso, hop, nuclei, pp_so)
    implicit none
    integer::k, i, j
    real*8, intent(in)::hop(nsiti,nsiti), nuclei(nsiti,3)
    complex*16, intent(out)::pp_so(3,nso,nso)
    complex*16::cplx, p(3,nsiti,nsiti)
    integer, intent(in):: nso, nsiti
    cplx=cmplx(0.d0,1.d0)
    p=0
    pp_so=0
    do k=1,3
       do i=1,nsiti
          do j=1,nsiti
             p(k,i,j)=cplx*hop(i,j)*(nuclei(i,k)-nuclei(j,k))
          enddo
       enddo
    enddo

    do K=1,3
       do i=1,nso-2
          pp_so(k,i,i+2)=p(k,(i+1)/2, (i+3)/2)
          pp_so(k,i+2,i)=p(k,(i+3)/2,(i+1)/2)
       enddo
    enddo
  end subroutine momentum_so

  subroutine write_matrix (matrix, numfile, righe, colonne, till)
    integer::i,j

    integer, intent(in):: numfile, till
    integer, intent(in)::righe, colonne
    real*8,intent(in)::matrix(righe,colonne)
    ! write(numfile,*) file_name
    do i=1,till
       write(numfile,'(<colonne>(2x,f10.5))') (matrix(i,j), j=1,colonne)
    enddo
  end subroutine write_matrix

  subroutine eigenvalues(dim2,tollerance,w,state)
    integer, intent(in) :: dim2
    real*8, intent(in) :: w(dim2)
    real*8, intent(in) :: tollerance
    character(len=1), intent(out) :: state(dim2)

    integer :: i, count_result, j

    do i = 1, dim2
       ! Inizializza il contatore
       count_result = 0

       ! Conta gli elementi che soddisfano la condizione
       do j = 1, dim2
          if (abs(w(j) - w(i)) < tollerance) then
             count_result = count_result + 1
          endif
       end do

       ! Determina lo stato in base al risultato del conteggio
       if (count_result == 3) then
          state(i) = 'T'
       elseif (count_result == 5) then
          state(i) = 'Q'
       else
          state(i) = 'S'
       endif
    end do
  end subroutine eigenvalues

  subroutine charge(carica, vecconfig, nz, dim2, nso)
    implicit none
    integer, intent(in):: dim2, nso, vecconfig(dim2), nz(nso/2)
    real*8, intent(out):: carica(dim2,nso/2)
    integer:: n, i, sito,a , b
    logical::bool, bool1

    carica=0
    do n=1,dim2
       do i=0,nso-2,2
          sito=(i+2)/2
          bool=btest(vecconfig(n), i)
          bool1=btest(vecconfig(n), i+1)
          if(bool)then
             a=1
          else
             a=0
          endif

          if(bool1)then
             b=1
          else
             b=0
          endif
          carica(n,sito)=nz(sito)-(a+b)
       enddo
    enddo
  end subroutine charge

  subroutine dipole_moment(dipole, carica, nuclei, dim2, nsiti)
    implicit none
    integer, intent(in)::dim2, nsiti
    real*8, intent(in):: carica(dim2,nsiti), nuclei(nsiti,3)
    real*8, intent(out):: dipole(dim2,3)
    integer::n,k,j

    dipole = 0.0d0

    do n = 1, dim2
       do k = 1, 3
          do j = 1, nsiti
             dipole(n, k) = dipole(n, k) + carica(n, j) * nuclei(j, k)
          enddo
       enddo
    enddo
  end subroutine dipole_moment

  subroutine square_complex_matrix(n, A)
    implicit none
    integer, intent(in) :: n
    complex*16, intent(inout) :: A(n,n)
    complex*16 :: temp(n,n)
    integer :: i, j, k

    ! Calcola il prodotto di A per se stessa
    do i = 1, n
       do j = 1, n
          temp(i,j) = (0.0d0, 0.0d0)
          do k = 1, n
             temp(i,j) = temp(i,j) + A(i,k) * A(k,j)
          end do
       end do
    end do

    ! Copia il risultato nel parametro di input 
    do i = 1, n
       do j = 1, n
          A(i,j) = temp(i,j)
       end do
    end do

  end subroutine square_complex_matrix

  subroutine vec_prod(v1, v2, risultato)
    complex*16,intent(in) :: v1(3), v2(3)
    complex*16, intent(out)::risultato(3)
    integer :: i

    ! Controllo delle dimensioni dei vettori
    if (size(v1) /= 3 .or. size(v2) /= 3) then
       print *, "Errore: i vettori devono essere di dimensione 3."
       return
    endif

    ! Calcolo del prodotto vettoriale
    risultato(1) = v1(2) * v2(3) - v1(3) * v2(2)
    risultato(2) = v1(3) * v2(1) - v1(1) * v2(3)
    risultato(3) = v1(1) * v2(2) - v1(2) * v2(1)

  end subroutine vec_prod

  subroutine calcola_distanze(coordinate, distanze, nsiti)
    implicit none
    real(8), dimension(nsiti, 3), intent(in) :: coordinate
    real(8), dimension(size(coordinate, 1), size(coordinate, 1)), intent(out) :: distanze
    integer :: i, j, n
    integer, intent(in):: nsiti
    real(8) :: dx, dy, dz

    n = size(coordinate, 1)

    ! Calcola le distanze tra tutti i punti
    do i = 1, n
       do j = i+1, n
          dx = coordinate(i, 1) - coordinate(j, 1)
          dy = coordinate(i, 2) - coordinate(j, 2)
          dz = coordinate(i, 3) - coordinate(j, 3)
          distanze(i, j) = dsqrt(dx**2 + dy**2 + dz**2)
          distanze(j, i) = distanze(i, j) ! Simmetria della matrice delle distanze
       end do
    end do
  end subroutine calcola_distanze

  subroutine rotate_real_1x2(dim2, out, input, ham)
    implicit none
    integer :: i, l, j
    integer, intent(in) :: dim2
    real*8, intent(out) :: out(dim2, dim2)
    real*8, intent(in) :: input(dim2)
    complex*16, intent(in) :: ham(dim2, dim2)

    !$omp parallel do default(none) private(i, l, j) shared(input, ham, dim2, out)
    do i = 1, dim2
       do l = 1, dim2
          do j = 1, dim2
             out(i, l) = out(i, l) + dconjg(ham(j, i)) * ham(j, l) * input(j)
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine rotate_real_1x2

  subroutine rotate_rtc_1x2(dim2,out,in,ham)
    implicit none
    integer::i,l,j,k
    integer, intent(in)::dim2
    complex*16,intent(out)::out(dim2,dim2)
    real*8,intent(in)::in(dim2)
    complex*16,intent(in)::ham(dim2,dim2)

    !$omp parallel do default(none) private(i, l, j) shared(in, ham, dim2, out)
    do i=1,dim2 !a
       do l=1,dim2 !b
          do j=1,dim2 !alfa 

             out(i,l)=out(i,l)+dconjg(ham(j,i))*ham(j,l)*in(j)

          enddo
       enddo
    enddo
    !$omp end parallel do
  end subroutine rotate_rtc_1x2

  subroutine rotate_cplx_1x2(dim2,out,in,ham)
    implicit none
    integer::i,l,j,k
    integer, intent(in)::dim2
    complex*16,intent(out)::out(dim2,dim2)
    complex*16,intent(in)::in(dim2)
    complex*16,intent(in)::ham(dim2,dim2)

    !$omp parallel do default(none) private(i, l, j) shared(in, ham, dim2, out)
    do i=1,dim2 !a
       do l=1,dim2 !b
          do j=1,dim2 !alfa 

             out(i,l)=out(i,l)+dconjg(ham(j,i))*ham(j,l)*in(j)

          enddo
       enddo
    enddo
    !$omp end parallel do
  end subroutine rotate_cplx_1x2

  !ruota una matrice reale non quadrata sullo spazio degli autostati
  !out is the resulting matrix
  !in is the matrix being rotated
  !rot is the complex matrix with eigenvectors
  !col is the number of columns

  subroutine check_q(matrix, n, is_hermitian)
    implicit none
    integer, intent(in) :: n
    complex*16, intent(in) :: matrix(n,n,n,n)
    logical, intent(out) :: is_hermitian
    integer :: i, j, k, l
    complex*16 :: conj_transpose_matrix(n, n,n,n)

    ! Calcola la trasposta coniugata della matrice
    do k=1,n
       do l=1,n
          do i = 1, n
             do j = 1, n
                conj_transpose_matrix(k,l,j,i) = dconjg(matrix(j,i,k,l))
             end do
          end do
       enddo
    enddo

    ! Check if the matrix is Hermitian 
    is_hermitian = all(zabs(matrix - conj_transpose_matrix).le.1d-14)

  end subroutine check_q

  subroutine bielectron(dim2, nso, pf, vecconfig, sootb, soo)
    implicit none
    integer, intent(in) :: dim2, nso
    integer, intent(in) :: vecconfig(dim2)
    complex*16, intent(out) :: soo(dim2, dim2)
    complex*16, intent(in) :: sootb(nso, nso, nso, nso)
    complex*16 :: phase
    complex*16, intent(in) :: pf
    integer :: n, a, b, c, d, m, temp, perm, i

    soo = 0.0d0

    !$omp parallel do default(none) private(n, a, b, c, d, m, temp, perm, phase, i) shared(dim2, nso, pf, vecconfig, sootb, soo)
    do n = 1, dim2
       do d = 0, nso - 1
          do c = 0, nso - 1
             do b = 0, nso - 1
                do a = 0, nso - 1
                   perm = 0
                   if (btest(vecconfig(n), d)) then
                      if (d /= nso - 1) then
                         do i = nso - 1, d + 1, -1
                            if (btest(vecconfig(n), i)) perm = perm + 1
                         enddo
                      endif
                      temp = ibclr(vecconfig(n), d)
                      if (btest(temp, c)) then
                         if (c /= nso - 1) then
                            do i = nso - 1, c + 1, -1
                               if (btest(temp, i)) perm = perm + 1
                            enddo
                         endif
                         temp = ibclr(temp, c)
                         if (.not. btest(temp, b)) then
                            if (b /= nso - 1) then
                               do i = nso - 1, b + 1, -1
                                  if (btest(temp, i)) perm = perm + 1
                               enddo
                            endif
                            temp = ibset(temp, b)
                            if (.not. btest(temp, a)) then
                               if (a /= nso - 1) then
                                  do i = nso - 1, a + 1, -1
                                     if (btest(temp, i)) perm = perm + 1
                                  enddo
                               endif
                               temp = ibset(temp, a)
                               m =  find_value(vecconfig, temp,dim2)
                               if (m /= 0) then
                                  if (mod(perm, 2) == 0) then
                                     phase = +1
                                  else
                                     phase = -1
                                  endif
                                  soo(n, m) = soo(n, m) + pf * phase * sootb(a + 1, b + 1, c + 1, d + 1)
                               endif
                            endif
                         endif
                      endif
                   endif
                enddo
             enddo
          enddo
       enddo
    end do
    !$omp end parallel do
  end subroutine bielectron

  subroutine rotate_cplx_2x2_multi(dim2, out, in, ham, strati)
    implicit none
    integer :: i, l, j, k, s
    integer, intent(in) :: dim2, strati
    complex*16, intent(out) :: out(strati, dim2, dim2)
    complex*16, intent(in) :: in(strati, dim2, dim2)
    complex*16, intent(in) :: ham(dim2, dim2)
    logical :: local_is_hermitian

    !$omp parallel do default(none) shared(out, in, ham, dim2, strati) private(s, i, l, j, k, local_is_hermitian)
    do s = 1, strati
       do i = 1, dim2
          do l = 1, dim2
             out(s, i, l) = (0.0d0, 0.0d0)
             do j = 1, dim2
                do k = 1, dim2
                   out(s, i, l) = out(s, i, l) + dconjg(ham(j, i)) * ham(k, l) * in(s, j, k)
                end do
             end do
          end do
       end do
       ! Che if the matrix is Hermitian in layer s
       local_is_hermitian = all(zabs(out(s, :, :) - conjg(transpose(out(s, :, :)))).le.1d-12)
       if(.not. local_is_hermitian) then
          write(*,*) 'Problem in strato s=',s
       endif
    end do
    !$omp end parallel do
  end subroutine rotate_cplx_2x2_multi

  subroutine check_symmetric(matrix, n, symmetric)
  implicit none
  integer, intent(in) :: n
  real*8, intent(in) :: matrix(n, n)
  logical, intent(out) :: symmetric
  integer :: i, j

  symmetric = .true.
  do i = 1, n
     do j = 1, n
        if (matrix(i, j) /= matrix(j, i)) then
           symmetric = .false.
           return
        end if
     end do
  end do
end subroutine check_symmetric

 function find_value(array, target, size) result(position)
    implicit none
    integer, intent(in) :: array(:)   ! Vettore in input
    integer, intent(in) :: size        ! Dimensione del vettore
    integer, intent(in) :: target      ! Elemento da cercare
    integer :: position                ! Posizione del valore cercato nel vettore

    integer :: i

    position = 0  ! Inizializza la posizione a zero

    ! Scansione del vettore per trovare il valore
    do i = 1, size
        if (array(i) == target) then
            position = i  ! Fix the position if the values is found
            return       ! Esci dalla funzione
        endif
    end do
  end function find_value

  function find_real(array, target, size) result(position)
    implicit none
    real*8, intent(in) :: array(:)  ! Vettore in input
    integer, intent(in) :: size      ! Dimensione del vettore
    real*8, intent(in) :: target    ! Elemento da cercare
    integer :: position              ! Posizione del valore cercato nel vettore
    integer :: i

    position = 0  ! Inizializza la posizione a zero

    ! Scansione del vettore per trovare il valore
    do i = 1, size
       if (array(i) == target) then
          position = i  ! Fix the position if the value is found
          return       ! Esci dalla funzione
       endif
    end do
  end function find_real

  subroutine bubble_sort(array, size, old_positions)
    implicit none
    real*8, intent(inout) :: array(:)      ! Vettore da ordinare
    integer, intent(in) :: size             ! Dimensione del vettore
    integer, intent(out) :: old_positions(:) ! Vettore colonna delle vecchie posizioni
    integer :: i, j
    real*8 :: temp
    integer :: temp_index
    
    ! Inizializza il vettore delle vecchie posizioni
    old_positions = [(i, i=1,size)]
    
    ! Algoritmo di bubble sort per ordinare il vettore
    do i = 1, size - 1
        do j = 1, size - i
            if (array(j) > array(j+1)) then
                ! Scambia gli elementi se non sono in ordine crescente
                temp = array(j)
                array(j) = array(j+1)
                array(j+1) = temp
                
                ! Scambia le posizioni corrispondenti nel vettore delle vecchie posizioni
                temp_index = old_positions(j)
                old_positions(j) = old_positions(j+1)
                old_positions(j+1) = temp_index
            endif
        end do
    end do
  end subroutine bubble_sort

  function find_character_position(array, size, target) result(position)
    implicit none
    character(len=*), intent(in) :: array(:)   ! Array di caratteri
    integer, intent(in) :: size                ! Dimensione dell'array
    character(len=2), intent(in) :: target     ! Carattere da cercare
    integer :: position                         ! Posizione del carattere
    integer :: i
    
    ! Inizializza la posizione a zero
    position = 0
    
    ! Scansione dell'array per trovare il carattere
    do i = 1, size
        if (array(i) == target) then
            position = i  ! FIx the position if the character is found 
            return       ! Esci dalla funzione
        endif
    end do
end function find_character_position
  
end module module


