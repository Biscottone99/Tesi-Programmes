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



!!$  recursive function binary_search(array, target, low, high) result(pos)
!!$    implicit none
!!$    integer, intent(in) :: array(:)
!!$    integer, intent(in) :: target
!!$    integer, intent(in) :: low, high
!!$    integer :: pos
!!$    integer :: mid
!!$    ! Se il limite inferiore supera il limite superiore, l'elemento non è presente
!!$    if (low > high) then
!!$        pos = 0 ! Valore che indica che il target non è stato trovato
!!$        return
!!$    end if
!!$    
!!$    ! Calcola il punto medio
!!$    
!!$    mid = (low + high) / 2
!!$    
!!$    ! Confronta il valore centrale con il target
!!$    if (array(mid) == target) then
!!$        pos = mid
!!$    else if (array(mid) > target) then
!!$        ! Cerca nella metà sinistra
!!$        pos = binary_search(array, target, low, mid - 1)
!!$    else
!!$        ! Cerca nella metà destra
!!$        pos = binary_search(array, target, mid + 1, high)
!!$    end if
!!$end function binary_search
!!$


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

                   i = binary_search(basis, istate, 1, dim) !row index
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

                   i = binary_search(basis, istate, 1, dim) !row index
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
                coupling(i,l)=coupling(i,l)+dconjg(ham(j,i))*ham(k,l)*coup(j,k)
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
  !out è la matrice che ottieni
  !in è la matrice che ruoti
  !rot è la matrice complessa con gli autovettori
  !col è il numero di colonne
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

    subroutine rot_diag_cplx(dim2, out, in, col, rot)
    implicit none
    integer :: i, j, k
    integer, intent(in) :: dim2, col
    complex*16, intent(in) :: rot(dim2, dim2)
    complex*16, intent(in) :: in(dim2, col)
    complex*16, intent(out) :: out(dim2, col)

    !$omp parallel do default(none) private(i, j, k) shared(dim2, col, rot, in, out)
    do i = 1, dim2
       do j = 1, dim2
          do k = 1, col
             out(i, k) = out(i, k) + dconjg(rot(j, i)) * rot(j, i) * in(j, k)
          end do
       end do
    end do
    !$omp end parallel do
  end subroutine rot_diag_cplx

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

    ! Verifica se la matrice è Hermitiana
    is_hermitian = all(zabs(matrix - conj_transpose_matrix).le.1d-10)

  end subroutine check_hermitian

  subroutine ppp_diag(dim2, nsiti, nuclei, esite, vecconfig, u, nz, pot)
    implicit none
    integer, intent(in) :: dim2, nsiti
    integer, intent(in) :: vecconfig(dim2), nz(nsiti)
    real*8, intent(in) :: nuclei(nsiti, 3), u(nsiti), esite(nsiti)
    real*8, intent(out) :: pot(dim2)
    integer :: n, i, j, p, sito, occupazioni(nsiti), a, b
    real*8 :: PPP, r(nsiti, nsiti), dx, dy, dz

    ! Initialize pot array to zero
    pot = 0.0d0
    ! Initialize r array to zero
    r = 0.0d0

    do i = 1, nsiti
       do j = i + 1, nsiti
          dx = nuclei(i, 1) - nuclei(j, 1)
          dy = nuclei(i, 2) - nuclei(j, 2)
          dz = nuclei(i, 3) - nuclei(j, 3)

          r(i, j) = dsqrt(dx**2 + dy**2 + dz**2)
          r(j, i) = r(i, j)
       end do
    end do

    do n = 1, dim2
       sito = 0
       do i = 0, 2 * nsiti - 2, 2
          sito = (i + 2) / 2
          a = 0
          b = 0
          if (btest(vecconfig(n), i)) a = 1
          if (btest(vecconfig(n), i + 1)) b = 1
          occupazioni(sito) = a + b
       end do

       do i = 1, nsiti
          pot(n) = pot(n) + esite(i) * occupazioni(i)
       end do

       PPP = 0.0d0
       do i = 1, nsiti
          do p = 1, nsiti
             if (i /= p) PPP = PPP + (14.397d0 / dsqrt(r(i, p)**2 + (28.794d0 / (u(i) + u(p)))**2)) * (nz(i) - occupazioni(i)) * (nz(p) - occupazioni(p))
             if ((i == p) .and. (occupazioni(i) == 2)) PPP = PPP + 2.d0 * u(i)
          end do
       end do
       pot(n) = pot(n) + 0.5d0 * PPP
    end do
  end subroutine ppp_diag

  subroutine site_energy_u(nso, dim2, esite, u, vecconfig, energy)
    implicit none
    integer :: n, i, sito
    real*8, intent(in) :: esite(nso/2), u(nso/2)
    integer, intent(in) :: dim2, vecconfig(dim2), nso
    real*8, intent(out) :: energy(dim2)

    ! Initialize energy array to zero
    energy = 0.0d0

    do n = 1, dim2
       ! Add single site energies
       do i = 0, nso-1
          sito = (i+2)/2
          if (btest(vecconfig(n), i)) energy(n) = energy(n) + esite(sito)
       end do

       ! Add interaction energies
       do i = 0, nso-2, 2
          sito = (i+2)/2
          if (btest(vecconfig(n), i) .and. btest(vecconfig(n), i+1)) energy(n) = energy(n) + u(sito)
       end do
    end do
  end subroutine site_energy_u

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
       elseif (count_result == 4) then
          state(i) = '4'
       elseif (count_result == 2) then
          state(i) = 'D'
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
!!$
!!$    subroutine copia_quad(matrice, dimensione)
!!$      complex(16), intent(inout) :: matrice(dimensione,dimensione)
!!$      integer, intent(in) :: dimensione
!!$      integer :: i, j
!!$
!!$      do i = 1, dimensione
!!$         do j = i+1, dimensione
!!$            matrice(j, i) = (matrice(i, j))
!!$         end do
!!$      end do
!!$
!!$    end subroutine copia_quad

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
  !out è la matrice che ottieni
  !in è la matrice che ruoti
  !rot è la matrice complessa con gli autovettori
  !col è il numero di colonne

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

    ! Verifica se la matrice è Hermitiana
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
    logical:: bool

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
                               m = binary_search(vecconfig, temp, 1, dim2)
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
    ! call check_hermitian(soo, dim2, bool)
    ! if(.not.bool) write(*,*) 'Not bool'
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
       ! Verifica se la matrice è Hermitiana per lo strato s
       local_is_hermitian = all(zabs(out(s, :, :) - conjg(transpose(out(s, :, :)))).le.1d-12)
       if(.not. local_is_hermitian) then
          write(*,*) 'Problem in strato s=',s
       endif
    end do
    !$omp end parallel do
  end subroutine rotate_cplx_2x2_multi


  subroutine compute_soc_mono(nsiti, dim2,nz, nuclei, hop, vecconfig, pf, coup, coupx, coupy, coupz)
    implicit none
    integer, intent(in) :: nsiti, dim2, nz(nsiti)
    real*8, intent(in) :: nuclei(nsiti, 3), hop(nsiti, nsiti)
    integer, intent(in) :: vecconfig(dim2)
    complex*16, intent(in) :: pf
    complex*16, intent(out) :: coup(dim2, dim2), coupx(dim2,dim2), coupy(dim2,dim2), coupz(dim2,dim2)
    integer :: i, j, k, isito, sitoi, sitoj, si, sj, nso
    complex*16 :: mom(nsiti, nsiti, 3)
    complex*16,allocatable:: soc_a(:,:,:), soc_b(:,:,:), soc_mono(:,:,:), hamsoc(:,:), tbcoupx(:,:), tbcoupy(:,:), tbcoupz(:,:)
    real*8 :: dr,  radius

    logical :: bool, is_hermitian
    complex*16 :: spin(2, 2, 3), cplx,vec1(3), vec2(3), cp(3)


    ! Define number of spin orbitals
    nso = 2 * nsiti
    allocate( soc_a(3, nso, nso), soc_b(3, nso, nso), soc_mono(3, nso, nso), hamsoc(nso, nso))
    allocate(tbcoupx(nso,nso), tbcoupy(nso,nso), tbcoupz(nso,nso))
    ! Define the complex unit
    cplx = (0.0d0, 1.0d0)

    ! Define the spin matrices
    spin = 0
    spin(1,2,1) = 1d0
    spin(2,1,1) = 1d0

    spin(1,2,2) = cplx
    spin(2,1,2) = -cplx

    spin(1,1,3) = 1d0
    spin(2,2,3) = -1d0



    ! Initialize the mom array
    mom = (0.0d0, 0.0d0)
    do i = 1, nsiti
       do j = 1, nsiti
          do k = 1, 3
             dr = nuclei(i, k) - nuclei(j, k)
             mom(i, j, k) = cplx * dr * hop(i, j)
          end do
       end do
    end do

    ! Check if mom is Hermitian
    do k = 1, 3
       call check_hermitian(mom(:, :, k), nsiti, bool)
       if (.not. bool) write(*, *) 'PROBLEMI MOM K=', k
    end do

    ! Initialize soc_a, soc_b, soc_mono arrays to zero
    soc_a = (0.0d0, 0.0d0)
    soc_b = (0.0d0, 0.0d0)
    soc_mono = (0.0d0, 0.0d0)

    ! Calculate soc_a
    do i = 1, nso !elettroni
       do j = 1, nso !elettroni
          do isito = 1, nsiti !nuclei
             if (isito /= (i + 1) / 2) then    !isito è il nucleo 
                vec1 = 0.0d0
                vec2 = 0.0d0
                cp = 0.0d0
                sitoi = (i + 1) / 2 !elettroni
                sitoj = (j + 1) / 2 !elettroni
                do k = 1, 3
                   vec1(k) = nuclei(sitoi, k) - nuclei(isito, k)
                   vec2(k) = mom(sitoi, sitoj, k)
                end do
                cp = cross_product(vec1, vec2)

                radius = 0.0d0
                do k = 1, 3
                   radius = radius + dreal(vec1(k))**2
                end do
                radius = (dsqrt(radius))**3
                cp = cp*nz(isito) / radius

                if (mod(i, 2) == 0) then !check
                   si = 2
                else
                   si = 1
                end if

                if (mod(j, 2) == 0) then
                   sj = 2
                else
                   sj = 1
                end if

                do k = 1, 3
                   soc_a(k, i, j) = soc_a(k, i, j) + cp(k) * spin(si, sj, k)
                end do
             end if
          end do
       end do
    end do

    ! Calculate soc_b
    do i = 1, nso
       do j = 1, nso
          do isito = 1, nsiti
             if (isito /= (j + 1) / 2) then
                vec1 = 0.0d0
                vec2 = 0.0d0
                cp = 0.0d0
                sitoi = (i + 1) / 2
                sitoj = (j + 1) / 2
                do k = 1, 3
                   vec1(k) = nuclei(sitoj, k) - nuclei(isito, k)
                   vec2(k) = mom(sitoi, sitoj, k)
                end do
                cp  = cross_product(vec1, vec2)

                radius = 0.0d0
                do k = 1, 3
                   radius = radius + dreal(vec1(k))**2
                end do
                radius = (dsqrt(radius))**3
                cp = cp * nz(isito) / radius

                if (mod(i, 2) == 0) then
                   si = 2
                else
                   si = 1
                end if

                if (mod(j, 2) == 0) then
                   sj = 2
                else
                   sj = 1
                end if

                do k = 1, 3
                   soc_b(k, i, j) = soc_b(k, i, j) + cp(k) * spin(si, sj, k)
                end do
             end if
          end do
       end do
    end do

    ! Calculate soc_mono
    do i = 1, nso
       do j = 1, nso
          do k = 1, 3
             soc_mono(k, i, j) = 0.5d0 * (soc_a(k, i, j) + soc_b(k, i, j))
          end do
       end do
    end do

    ! Initialize and calculate hamsoc
    hamsoc = (0.0d0, 0.0d0)
    do i = 1, nso
       do j = 1, nso
          do k = 1, 3
             hamsoc(i, j) = hamsoc(i, j) + soc_mono(k, i, j)
          end do
       end do
    end do
    coupx=0
    coupy=0
    coupz=0
    tbcoupx = soc_mono(1, :, :)
    tbcoupy = soc_mono(2, :, :)
    tbcoupz = soc_mono(3, :, :)



    ! Calculate coup
    coup = 0.0d0
    call sq_oe_op_compl(nso, dim2, hamsoc, coup, vecconfig)
    call check_hermitian(coup, dim2, is_hermitian)
    if (is_hermitian) write(*, *) 'COUP HERMITIANA'

    call sq_oe_op_compl(nso, dim2, tbcoupx, coupx, vecconfig)
    call check_hermitian(coupx, dim2, is_hermitian)
    if (.not.is_hermitian) write(*, *) 'COUPX Problem'

    call sq_oe_op_compl(nso, dim2, tbcoupy, coupy, vecconfig)
    call check_hermitian(coupy, dim2, is_hermitian)
    if (.not.is_hermitian) write(*, *) 'COUPY Problem'

    call sq_oe_op_compl(nso, dim2, tbcoupz, coupz, vecconfig)
    call check_hermitian(coupz, dim2, is_hermitian)
    if (.not.is_hermitian) write(*, *) 'COUPZ Problem'

    coup = pf * coup
    coupx = pf * coupx
    coupy = pf * coupy
    coupz = pf * coupz

  end subroutine compute_soc_mono


  subroutine compute_sso(nsiti, dim2, nuclei, hop, vecconfig, pf, sso, ssox, ssoy, ssoz)
    implicit none
    integer, intent(in) :: nsiti, dim2
    real*8, intent(in) :: nuclei(nsiti, 3), hop(nsiti, nsiti)
    integer, intent(in) :: vecconfig(dim2)
    complex*16, intent(in) :: pf
    complex*16, intent(out) :: sso(dim2, dim2), ssox(dim2, dim2), ssoy(dim2, dim2), ssoz(dim2, dim2)

    complex*16,allocatable:: hssotb(:,:,:,:,:), ssotb(:,:,:,:), ppso(:,:,:),ssotbx(:,:,:,:),ssotby(:,:,:,:),ssotbz(:,:,:,:), hssotb2(:,:,:,:,:)
    complex*16 ::  spin(2, 2, 3), cplx, vec1(3), vec2(3), cp(3), cp2(3)
    real*8 :: dist(nsiti, nsiti, 3), radius
    integer :: i, j, k, a, b, c, d, si, sj, nso, asito, bsito, csito, dsito
    logical :: bool

    nso = 2 * nsiti
    allocate(hssotb(3, nso, nso, nso, nso), ssotb(nso, nso, nso, nso), ppso(3, nsiti, nsiti),hssotb2(3, nso, nso, nso, nso))
    allocate(ssotbx(nso, nso, nso, nso),ssotby(nso, nso, nso, nso),ssotbz(nso, nso, nso, nso))

    ! Define the complex unit
    cplx = (0.0d0, 1.0d0)

    ! Define the spin matrices
    spin = 0
    spin(1, 2, 1) = 1d0
    spin(2, 1, 1) = 1d0

    spin(1, 2, 2) = cplx
    spin(2, 1, 2) = -cplx

    spin(1, 1, 3) = 1d0
    spin(2, 2, 3) = -1d0

    ! Allocate and initialize the dist array
    ppso=0d0
    do k=1,3
       do i=1,nsiti
          do j=1,nsiti
             ppso(k,i,j)=-cplx*hop(i,j)*(nuclei(i,k)-nuclei(j,k))
          enddo
       enddo
    enddo

    ! Initialize hssotb array
    hssotb = (0.0d0, 0.0d0)

    ! Calculate hssotb
    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             asito = (a+1)/2
             bsito = (b+1)/2
             csito = (c+1)/2
             
             if (asito.ne.bsito)then

                vec1 = 0.0d0
                vec2 = 0.0d0
                cp = 0.0d0
                do k = 1, 3
                   vec1(k) = nuclei(asito, k) - nuclei(bsito, k) !r_ab
                   vec2(k) = ppso(k, asito, csito) !p_ac
                end do
                cp  = cross_product(vec1, vec2)

                radius = 0.0d0
                do k = 1, 3
                   radius = radius + dreal(vec1(k))**2
                end do
                radius = (dsqrt(radius))**3
                cp = cp / radius
             else
                cp =0d0
             endif

             if(bsito.ne.csito)then
                vec1 = 0.0d0
                vec2 = 0.0d0
                cp2 = 0.0d0
                do k = 1, 3
                   vec2(k) = nuclei(csito, k) - nuclei(bsito, k) !r_cb
                   vec1(k) = ppso(k, asito, csito) !p_ac
                end do
                cp2  = cross_product(vec1, vec2)

                radius = 0.0d0
                do k = 1, 3
                   radius = radius + dreal(vec2(k))**2
                end do
                radius = (dsqrt(radius))**3
                cp2 = cp2 / radius
             else
                cp2=0d0
             endif

             if (mod(a, 2) == 0) then  !spin_a
                si = 2
             else
                si = 1
             end if

             if (mod(c, 2) == 0) then !spin_b
                sj = 2
             else
                sj = 1
             end if
             
             do k = 1, 3
                hssotb(k, a, b, c, b) = hssotb(k, a, b, c, b) + (cp(k)-cp2(k)) * spin(si, sj, k)
             end do

          end do
       end do
    end do
    ! Initialize ssotb array
    ssotb = (0.0d0, 0.0d0)

    ! Calculate ssotb
   
    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             do d = 1, nso
                do k = 1, 3
                   ssotb(a, b, c, d) = ssotb(a, b, c, d) + 0.5d0 * (hssotb(k, a, b, c, d))
                end do
             !  if(zabs(ssotb(a,b,c,d)).ge.1d-10) write(1,*) a, b, c, d, ssotb(a, b, c, d)
             end do
          end do
       end do
    end do
    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             do d = 1, nso
                ssotbx(a, b, c, d) = ssotbx(a, b, c, d) + 0.5d0 * (hssotb(1, a, b, c, d))
                ssotby(a, b, c, d) = ssotby(a, b, c, d) + 0.5d0 * (hssotb(2, a, b, c, d) )
                ssotbz(a, b, c, d) = ssotbz(a, b, c, d) + 0.5d0 * (hssotb(3, a, b, c, d) )
               ! if(zabs(ssotbz(a,b,c,d)).ge.1d-10) write(1,*) a, b, c, d, ssotbz(a, b, c, d)
             end do
          end do
       end do
    end do

    ! Initialize and calculate sso
    sso = 0.0d0
    call bielectron(dim2, nso, pf, vecconfig, ssotb, sso)
    call bielectron(dim2, nso, pf, vecconfig, ssotbx, ssox)
    call bielectron(dim2, nso, pf, vecconfig, ssotby, ssoy)
    call bielectron(dim2, nso, pf, vecconfig, ssotbz, ssoz)

    call check_hermitian(sso, dim2, bool)
    if ( bool) write(*, *) 'SSO Hermitian'

    call check_hermitian(ssoy, dim2, bool)
    if ( .not.bool) write(*, *) 'SSOY Problem'

    call check_hermitian(ssox, dim2, bool)
    if ( .not.bool) write(*, *) 'SSOX Problem'

    call check_hermitian(ssoz, dim2, bool)
    if ( .not.bool) write(*, *) 'SSOZ Problem'
  end subroutine compute_sso



  subroutine compute_soo(nsiti, dim2, nuclei, hop, vecconfig, pf, soo, soox, sooy, sooz)
    implicit none
    integer, intent(in) :: nsiti, dim2
    real*8, intent(in) :: nuclei(nsiti, 3), hop(nsiti,nsiti)
    complex*16, intent(in) :: pf
    integer, intent(in) :: vecconfig(dim2)
    complex*16, intent(out) :: soo(dim2, dim2), soox(dim2,dim2), sooy(dim2,dim2), sooz(dim2,dim2)

    complex*16, allocatable:: hsootb(:,:,:,:,:), sootb(:,:,:,:), ppso(:,:,:), check(:,:),  conj_transpose_matrix(:, :, :, :), hsootb2(:,:,:,:,:)
    complex*16,allocatable::sootbx(:,:,:,:),sootby(:,:,:,:),sootbz(:,:,:,:)
    complex*16 :: vec1(3), vec2(3), cp(3), cp2(3)
    complex*16 :: spin(2, 2, 3), cplx
    real*8 :: dist(nsiti, nsiti, 3), radius
    integer :: i, j, k, a, b, c, d, si, sj, nso, asito, bsito, csito, dsito, sa, sb
    logical :: bool


    nso=nsiti*2
    allocate(hsootb(3, nso, nso, nso, nso), sootb(nso, nso, nso, nso), ppso(3, nsiti, nsiti),check(nso,nso))
    allocate(sootbx(nso, nso, nso, nso),sootby(nso, nso, nso, nso),sootbz(nso, nso, nso, nso))
    ! Define the complex unit
    cplx = (0.0d0, 1.0d0)

    ! Define the spin matrices
    spin = 0
    spin(1, 2, 1) = 1d0
    spin(2, 1, 1) = 1d0

    spin(1, 2, 2) = cplx
    spin(2, 1, 2) = -cplx

    spin(1, 1, 3) = 1d0
    spin(2, 2, 3) = -1d0

    ! Allocate arrays
    ppso=0d0
    do i =1, nsiti
       do j = 1, nsiti
          do k = 1, 3
             ppso(k,i,j)= -cplx * hop(i,j) * (nuclei(i,k)-nuclei(j,k))
          enddo
       enddo
    enddo
    
    
    ! Initialize hsootb array
    hsootb = (0.0d0, 0.0d0)

    ! Calculate hsootb
    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             do  d = 1, nso
                asito = (a+1)/2
                bsito = (b+1)/2
                csito = (c+1)/2
                dsito = (d+1)/2
                if(mod(a,2).eq.0)then
                   sa = 2
                else
                   sa = 1
                endif

                if(mod(c,2).eq.0)then
                   sb = 2
                else
                   sb = 1
                endif

                if((sa.eq.sb).and.(bsito.eq.dsito))then
                   if(asito.ne.bsito)then
                      vec1 = 0.0d0
                      vec2 = 0.0d0
                      cp = 0.0d0
                      do k = 1, 3
                         vec1(k) = nuclei(asito, k) - nuclei(bsito, k)
                         vec2(k) = ppso(k, asito,csito)
                      end do
                      cp = cross_product(vec1, vec2)

                      radius = 0.0d0
                      do k = 1, 3
                         radius = radius + dreal(vec1(k))**2
                      end do
                      radius = (dsqrt(radius))**3
                      cp = cp / radius
                   else
                      cp=0d0
                   endif

                   if(csito.ne.bsito)then
                      vec1 = 0.0d0
                      vec2 = 0.0d0
                      cp2 = 0.0d0
                      do k = 1, 3
                         vec1(k) = ppso(k, asito, csito)
                         vec2(k) =nuclei( csito, k) - nuclei( bsito, k)
                      end do
                      cp2 = cross_product(vec1, vec2)

                      radius = 0.0d0
                      do k = 1, 3
                         radius = radius + dreal(vec2(k))**2
                      end do
                      radius = (dsqrt(radius))**3
                      cp2 = cp2 / radius
                   else
                      cp2=0d0
                   endif

                   if (mod(b, 2) == 0) then
                      si = 2
                   else
                      si = 1
                   end if

                   if (mod(d, 2) == 0) then
                      sj = 2
                   else
                      sj = 1
                   end if
                   do k = 1, 3
                      hsootb(k, a, b, c, d) = hsootb(k, a, b, c, d) + (cp(k) - cp2(k)) * spin(si, sj, k)
                   end do
                endif
             enddo
          end do
       enddo
    enddo



    ! Initialize sootb array
    sootb = (0.0d0, 0.0d0)
    sootbx = (0.0d0, 0.0d0)
    sootby = (0.0d0, 0.0d0)
    sootbz = (0.0d0, 0.0d0)

    ! Calculate sootb
    
    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             do d = 1, nso
                do k = 1, 3
                   ! sootb(a, b, c, d) = sootb(a, b, c, d) + 0.5d0 * (hsootb(k, a, b, c, d) + dconjg(hsootb(k, c, d, a, b)))
                   sootb(a, b, c, d) = sootb(a, b, c, d) + 0.5*hsootb(k, a, b, c, d)
                end do
              ! if(zabs(sootb(a,b,c,d)).ge.1d-4) write(1,'(<4>(I3, 2x), 2x, <2>(f10.5, 2x))') a, b, c, d, dreal(sootb(a,b,c,d)), dimag(sootb(a,b,c,d)) 
             end do
          end do
       end do
    end do

    
    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             do d = 1, nso
                sootbx(a, b, c, d) = sootbx(a, b, c, d) + 0.5d0 * (hsootb(1, a, b, c, d))
                sootby(a, b, c, d) = sootby(a, b, c, d) + 0.5d0 * (hsootb(2, a, b, c, d))
                sootbz(a, b, c, d) = sootbz(a, b, c, d) + 0.5d0 * (hsootb(3, a, b, c, d))
             end do
          end do
       end do
    end do

    allocate(conj_transpose_matrix(nso, nso, nso, nso))
    do a = 1, nso
       do b = 1, nso
          do c = 1, nso
             do d = 1, nso
                conj_transpose_matrix(a, b, c, d) = dconjg(sootb(c, d, a, b))
             end do
          end do
       enddo
    enddo
    bool = all(zabs(sootb - conj_transpose_matrix).le.1d-8)
    IF(bool)write(*,*) 'SOOTB hermitian'
    ! Initialize and calculate soo
    soo = 0.0d0
    soox = 0d0
    sooy = 0d0
    sooz = 0d0
    call bielectron(dim2, nso, pf, vecconfig, sootb, soo)
    call bielectron(dim2, nso, pf, vecconfig, sootbx, soox)
    call bielectron(dim2, nso, pf, vecconfig, sootby, sooy)
    call bielectron(dim2, nso, pf, vecconfig, sootbz, sooz)

    call check_hermitian(soo, dim2, bool)
    if (bool) write(*, *) 'SOO hermitian'

    call check_hermitian(sooy, dim2, bool)
    if ( .not.bool) write(*, *) 'SOOY Problem'

    call check_hermitian(soox, dim2, bool)
    if ( .not.bool) write(*, *) 'SOOX Problem'

    call check_hermitian(sooz, dim2, bool)
    if ( .not.bool) write(*, *) 'SOOZ Problem'

  end subroutine compute_soo


  subroutine op_siti_2_so(hop_so, hop_use, nso, nsiti)
    implicit none
    integer, intent(in) :: nsiti,nso
    real*8, intent(in) :: hop_use(nsiti, nsiti)
    real*8, intent(out) :: hop_so(nso, nso)
    integer :: i, j, isito, jsito

    hop_so = 0.0d0
    do i = 1, nso
       do j = 1, nso
          isito = (i + 1) / 2
          jsito = (j + 1) / 2
          if ((hop_use(isito, jsito) /= 0.0d0) .and. (mod(i, 2) == mod(j, 2))) then
             hop_so(i, j) = hop_use(isito, jsito)
          end if
       end do
    end do
  end subroutine op_siti_2_so

  subroutine spectra(nev, autovalori, mu, nomega,omegastart, omegaend, unit, intensita)
    implicit none
    integer, intent(in) :: nev
    real*8, intent(in) :: autovalori(nev), omegastart, omegaend
    complex*16, intent(in) :: mu(nev, nev, 3)
    integer, intent(in) :: nomega
    real*8, intent(out) :: intensita(nomega, 2)
    character(len=*), intent(in) :: unit

    ! Dichiarazione delle variabili
    real*8 :: sigma2, omega, domega, int, pi, max_intensity
    integer :: i, j

    ! Inizializzazione dei parametri
    sigma2 = 0.1d0
    !omegastart = 0.d0
    ! omegaend = 10d0
    domega = (omegaend - omegastart) / nomega
    omega = omegastart - domega
    pi = dacos(-1.d0)

    ! Ciclo principale
    do j = 1, nomega
       omega = omega + domega
       int = 0.d0
       do i = 2, nev
          int = int + (autovalori(i) - autovalori(1)) * &
               ((zabs(mu(i,i,1)))**2 + (zabs(mu(i,i,2)))**2 + (zabs(mu(i,i,3)))**2) * &
               ((2 * pi)**(-0.5) * sigma2) * &
               dexp(-(omega - (autovalori(i) - autovalori(1)))**2 / (2 * sigma2**2))
       end do
       if (unit .eq. 'E') then
          intensita(j,1) = omega
       elseif (unit .eq. 'N') then
          intensita(j,1) = omega * 1240.d0
       elseif (unit .eq. 'C') then
          intensita(j,1) = omega * 8065.6d0
       else
          print *, 'Unrecognized unit: ', unit
          stop
       end if
       intensita(j,2) = int
    end do

    ! Normalizzazione dell'intensità
    max_intensity = maxval(intensita(:,2))
    do i = 1, nomega
       intensita(i,2) = intensita(i,2) / max_intensity
    end do

  end subroutine spectra

  subroutine sq_op_conf(nso, dim2, basis, sq)
    implicit none
    integer, intent(in) :: nso, dim2, basis(dim2)
    complex*16, intent(out) :: sq(dim2, dim2)
    complex*16 :: spin(2, 2, 3), sx(nso, nso), sy(nso, nso), sz(nso, nso)
    complex*16 :: ssqx(dim2, dim2), ssqy(dim2, dim2), ssqz(dim2, dim2)
    integer :: i, n, m
    complex*16 :: cplx

    ! Initialize the spin matrices
    spin = 0.0d0
    spin(1, 2, 1) = 1.0d0
    spin(2, 1, 1) = 1.0d0

    cplx = (0.0d0, 1.0d0)  ! Complex unit imaginary number
    spin(1, 2, 2) = cplx
    spin(2, 1, 2) = -cplx

    spin(1, 1, 3) = 1.0d0
    spin(2, 2, 3) = -1.0d0

    sx = 0.0d0
    sy = 0.0d0
    sz = 0.0d0

    ! Populate the spin matrices sx, sy, sz
    do i = 1, nso-1, 2
       sx(i, i) = spin(1, 1, 1)
       sx(i+1, i+1) = spin(2, 2, 1)
       sx(i, i+1) = spin(1, 2, 1)
       sx(i+1, i) = spin(2, 1, 1)

       sy(i, i) = spin(1, 1, 2)
       sy(i+1, i+1) = spin(2, 2, 2)
       sy(i, i+1) = spin(1, 2, 2)
       sy(i+1, i) = spin(2, 1, 2)

       sz(i, i) = spin(1, 1, 3)
       sz(i+1, i+1) = spin(2, 2, 3)
       sz(i, i+1) = spin(1, 2, 3)
       sz(i+1, i) = spin(2, 1, 3)
    end do

    ! Compute the square of each spin matrix
    call sq_oe_op_compl(nso, dim2, sx, ssqx, basis)
    call sq_oe_op_compl(nso, dim2, sy, ssqy, basis)
    call sq_oe_op_compl(nso, dim2, sz, ssqz, basis)

    ssqx = 0.5 * ssqx
    call square_complex_matrix(dim2, ssqx)

    ssqy = 0.5 * ssqy
    call square_complex_matrix(dim2, ssqy)

    ssqz = 0.5 * ssqz
    call square_complex_matrix(dim2, ssqz)

    ! Combine the results to form sq
    sq = 0.0d0
    do n = 1, dim2
       do m = 1, dim2
          sq(n, m) = ssqx(n, m) + ssqy(n, m) + ssqz(n, m)
       end do
    end do

  end subroutine sq_op_conf

  function find_real(vec, target, n) result(position)
    implicit none
    integer :: n                 ! Dimensione del vettore
    real*8, intent(in) :: vec(n)  ! Vettore di input
    real*8, intent(in) :: target  ! Numero da cercare
    integer :: position           ! Risultato: posizione del numero trovato
    integer :: i                  ! Contatore del ciclo

    position = -1                 ! Valore predefinito, se non trovato

    ! Ciclo semplice per trovare il numero nel vettore
    do i = 1, n
       if (vec(i) == target) then
          position = i              ! Memorizza la posizione del valore trovato
          return                    ! Esce immediatamente se il valore è trovato
       end if
    end do

  end function find_real

  subroutine bubble_sort(vec, n, indices, newpositions)
    implicit none
    integer, intent(in) :: n                ! Dimensione del vettore
    real*8, intent(inout) :: vec(n)         ! Vettore di input (viene ordinato in output)
    integer, intent(out) :: indices(n)      ! Vettore delle posizioni originali
    integer, intent(out) :: newpositions(n) ! Nuove posizioni dopo l'ordinamento
    real*8, allocatable :: sorted_vec(:)    ! Vettore temporaneo per il processo di ordinamento

    integer :: i, j, min_idx
    real*8 :: temp
    integer :: temp_idx

    allocate(sorted_vec(n))

    ! Copia il vettore di input nel vettore ordinato
    sorted_vec = vec

    ! Inizializza l'array delle posizioni originali
    do i = 1, n
       indices(i) = i
    end do

    ! Algoritmo di Selection Sort con traccia degli indici originali
    do i = 1, n - 1
       min_idx = i
       do j = i + 1, n
          if (sorted_vec(j) < sorted_vec(min_idx)) then
             min_idx = j
          end if
       end do

       ! Scambia gli elementi del vettore ordinato
       if (min_idx /= i) then
          temp = sorted_vec(i)
          sorted_vec(i) = sorted_vec(min_idx)
          sorted_vec(min_idx) = temp

          ! Scambia anche gli indici corrispondenti
          temp_idx = indices(i)
          indices(i) = indices(min_idx)
          indices(min_idx) = temp_idx
       end if
    end do

    ! Mappa le nuove posizioni per ogni valore originale
    do i = 1, n
       newpositions(indices(i)) = i
    end do

    ! Restituisce il vettore ordinato in `vec`
    vec = sorted_vec

    ! Dealloca la memoria allocata dinamicamente
    deallocate(sorted_vec)

  end subroutine bubble_sort

  subroutine recursive_hamiltonian(nsiti, dim, hop_use, basis, PPPflag, hubbardflag, coord, esite, u, nz, eigenvalues, hamiltonian)
    implicit none
    integer, intent(in) :: dim, nsiti
    integer, intent(in) :: nz(nsiti), basis(dim) 
    logical, intent(in) :: PPPflag, hubbardflag
    double precision, intent(in) :: hop_use(nsiti, nsiti),esite(nsiti), u(nsiti)
    double precision, intent(in) :: coord(nsiti, 3)
    double precision, intent(out) :: eigenvalues(dim)
    complex*16, intent(out):: hamiltonian(dim, dim)
    double precision,allocatable ::  pot(:)
    real*8, allocatable :: hop_so(:,:), hop(:,:), rwork(:)
    integer :: i, j, lrwork, liwork, lwork, info, nso
    complex*16, allocatable :: work(:)
    integer, allocatable :: iwork(:)
    logical :: bool

    nso = 2 * nsiti  ! Numero di orbitali di spin

    ! Allocazione di matrici e vettori
    allocate(hop(dim, dim), pot(dim) )
    allocate(hop_so(nso, nso))
    hop_so = 0.0d0

    ! Popolamento della matrice hop_so
    call op_siti_2_so(hop_so, hop_use, nso, nsiti)

    ! Popolamento della matrice hop
    call sq_oe_op_real(nso, dim, hop_so, hop, basis)

    ! Applicazione dell'Hamiltoniano PPP se il flag PPPflag è attivo
    if (PPPflag) call ppp_diag(dim, nsiti, coord, esite, basis, u, nz, pot)

    ! Applicazione dell'energia di sito Hubbard se il flag hubbardflag è attivo
    if (hubbardflag) call site_energy_u(nso, dim, esite, u, basis, pot)

    ! Inizializzazione dell'Hamiltoniano a zero
    hamiltonian = 0.0d0

    ! Popolamento della matrice di Hamiltoniano
    do i = 1, dim
       hamiltonian(i, i) = hamiltonian(i, i) + pot(i)
       do j = 1, dim
          hamiltonian(i, j) = hamiltonian(i, j) - hop(i, j)
       end do
    end do

    ! Controllo se l'Hamiltoniano è ermitiano
    call check_hermitian(hamiltonian, dim, bool)
    if (.not. bool) then
       write(*,*) 'Problem in Hamiltonian'
       stop
    end if

    ! Dimensioni delle aree di lavoro per la diagonalizzazione con zheevd
    lrwork = (1 + 5*dim + 2*dim**2)
    liwork = (3 + 5*dim)
    lwork = (2*dim + dim**2)

    ! Allocazione delle aree di lavoro
    allocate(work(max(1, lwork)), rwork(lrwork), iwork(max(1, liwork)))

    ! Diagonalizzazione dell'Hamiltoniano con zheevd
    call zheevd('V', 'U', dim, hamiltonian, dim, eigenvalues, work, lwork, rwork, lrwork, iwork, liwork, info)


    ! Deallocazione delle matrici e dei vettori
    deallocate(hop, pot, work, rwork, iwork,hop_so)

  end subroutine recursive_hamiltonian

  function binary_search(array, target, usefull, n) result(pos)
    implicit none
    integer, intent(in) :: array(:), target, usefull, n
    integer :: pos
    integer :: i

    ! Inizializza la posizione come -1 (non trovato)
    pos = 0

    ! Ciclo attraverso il vettore per trovare il target
    do i = 1, n
       if (array(i) == target) then
          pos = i
          return
       end if
    end do
  end function binary_search


  subroutine sz_on_conf(dim, nso, basis, sz)
    integer, intent(in)::dim, nso, basis(dim)
    real*8, intent(out):: sz(dim)
    integer:: i, j
    sz=0d0
    do i = 1, dim
       do j = 0, nso-1
          if(btest(basis(i),j))then
             if(mod(j,2).eq.0)then
                sz(i)=sz(i)+0.5
             else
                sz(i)=sz(i)-0.5
             end if
          end if
       end do
    enddo
  end subroutine sz_on_conf


  SUBROUTINE ordina_autovalori(n, sz1, energia, posizioni)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n                 ! Dimensione del vettore
    real*8, DIMENSION(n), INTENT(IN) :: sz1      ! Vettore dei valori di Sz
    REAL*8, DIMENSION(n), INTENT(IN) :: energia ! Vettore delle energie
    INTEGER, DIMENSION(n), INTENT(OUT) :: posizioni  ! Vettore delle nuove posizioni

    INTEGER :: i, j, temp_pos, sz(n)
    REAL :: temp_sz, temp_energia
    do i =1, n
       sz(i)=nint(sz1(i))
    end do
    ! Inizializza il vettore delle posizioni originali
    DO i = 1, n
       posizioni(i) = i
    END DO

    ! Ordinamento "bubble sort" per semplicità
    DO i = 1, n-1
       DO j = i+1, n
          IF ( (sz(posizioni(i)) > sz(posizioni(j))) .OR. &
               ((sz(posizioni(i)) == sz(posizioni(j))) .AND. (energia(posizioni(i)) > energia(posizioni(j)))) ) THEN
             ! Scambio le posizioni
             temp_pos = posizioni(i)
             posizioni(i) = posizioni(j)
             posizioni(j) = temp_pos
          END IF
       END DO
    END DO

  END SUBROUTINE ordina_autovalori

  SUBROUTINE riordina_cplx(n, mat_in, posizioni)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n                ! Dimensione della matrice (n righe, n colonne)
    COMPLEX(16), DIMENSION(n,n), INTENT(INOUT) :: mat_in  ! Matrice complessa di input/output
    INTEGER, DIMENSION(n), INTENT(IN) :: posizioni  ! Vettore delle nuove posizioni per le colonne
    COMPLEX(16), DIMENSION(n,n) :: mat_out  ! Matrice temporanea per la riordinazione

    INTEGER :: j

    ! Riordina le colonne della matrice in base al vettore di posizioni
    DO j = 1, n
       mat_out(:,j) = mat_in(:,posizioni(j))  ! Sposta la colonna j nella nuova posizione
    END DO

    ! Copia la matrice riordinata in mat_in
    mat_in = mat_out

  END SUBROUTINE riordina_cplx
  
  subroutine s2_realspace(dim, nso, basis, szo, sq)
    integer, intent(in) :: dim, nso, basis(dim)
    complex*16, intent(out) :: sq(dim, dim), szo(dim)

    integer :: i, j
    complex*16 :: cplx, spin(2,2,3)
    complex*16, allocatable :: sx(:,:), sy(:,:), sz(:,:), sxrs(:,:), syrs(:,:), szrs(:,:)
    logical :: is_hermitian
   
    ! Definisco cplx come unità immaginaria
    cplx = cmplx(0.d0, 1.d0)

    ! Inizializzo la matrice degli operatori di spin
    spin = 0
    spin(1,2,1) = 1.d0
    spin(2,1,1) = 1.d0

    spin(1,2,2) = cplx
    spin(2,1,2) = -cplx

    spin(1,1,3) = 1.d0
    spin(2,2,3) = -1.d0

    !========================= Operatori di spin =========================
    allocate(sx(nso, nso), sy(nso, nso), sz(nso, nso))
    sx = 0
    sy = 0
    sz = 0

    ! Costruisco l'operatore sx
    do i = 1, nso - 1, 2
       sx(i,i) = spin(1,1,1)
       sx(i+1,i+1) = spin(2,2,1)
       sx(i,i+1) = spin(1,2,1)
       sx(i+1,i) = spin(2,1,1)
    end do

    ! Costruisco l'operatore sy
    do i = 1, nso - 1, 2
       sy(i,i) = spin(1,1,2)
       sy(i+1,i+1) = spin(2,2,2)
       sy(i,i+1) = spin(1,2,2)
       sy(i+1,i) = spin(2,1,2)
    end do

    ! Costruisco l'operatore sz
    do i = 1, nso - 1, 2
       sz(i,i) = spin(1,1,3)
       sz(i+1,i+1) = spin(2,2,3)
       sz(i,i+1) = spin(1,2,3)
       sz(i+1,i) = spin(2,1,3)
    end do

    !========================= Trasformazione in real space =========================
    allocate(sxrs(dim, dim), syrs(dim, dim), szrs(dim,dim))

    ! Trasformazione di sx
    call sq_oe_op_compl(nso, dim, sx, sxrs, basis)
    call check_hermitian(sxrs, dim, is_hermitian)
    if (.not. is_hermitian) write(*,*) 'Problem sxrs'

    ! Trasformazione di sy
    call sq_oe_op_compl(nso, dim, sy, syrs, basis)
    call check_hermitian(syrs, dim, is_hermitian)
    if (.not. is_hermitian) write(*,*) 'Problem syrs'

    ! Trasformazione di sz
    call sq_oe_op_compl(nso, dim, sz, szrs, basis)
    call check_hermitian(szrs, dim, is_hermitian)
    if (.not. is_hermitian) write(*,*) 'Problem szrs'

    !========================= Calcolo degli operatori al quadrato =========================
    sxrs = 0.5 * sxrs
    call square_complex_matrix(dim, sxrs)

    syrs = 0.5 * syrs
    call square_complex_matrix(dim, syrs)

    szrs = 0.5 * szrs
    do i = 1, dim
       szo(i) = szrs(i,i)
    enddo
    call square_complex_matrix(dim, szrs)

    ! Controllo se le matrici al quadrato sono hermitiane
    call check_hermitian(sxrs, dim, is_hermitian)
    if (.not. is_hermitian) write(*,*) 'Not hermitian sx^2'

    call check_hermitian(syrs, dim, is_hermitian)
    if (.not. is_hermitian) write(*,*) 'Not hermitian sy^2'

    call check_hermitian(szrs, dim, is_hermitian)
    if (.not. is_hermitian) write(*,*) 'Not hermitian sz^2'

    !========================= Calcolo di S^2 =========================
    sq = 0
    do i = 1, dim
       do j = 1, dim
          sq(i,j) = sq(i,j) + sxrs(i,j) + syrs(i,j) + szrs(i,j)
       end do
    end do

  end subroutine s2_realspace


  subroutine sort_eigenvalues(n, s2, sz, eigenvalues,  sorted_values, old_positions)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: s2(n), sz(n), eigenvalues(n)
    real(8), intent(out) :: sorted_values(n)
    integer, intent(out) :: old_positions(n)

    real(8) :: s2_local(n), sz_local(n), eigenvalues_local(n)
    integer :: i, j, min_index
    real(8) :: temp_s2, temp_sz, temp_eigenvalue
    integer :: temp_pos
    real(8), parameter :: threshold = 1.0d-4  ! Treshold per il confronto

    ! Creare copie locali dei vettori
    s2_local = s2
    sz_local = sz
    eigenvalues_local = eigenvalues

    ! Inizializzare old_positions con i loro indici originali
    do i = 1, n
       old_positions(i) = i
    end do

    ! Selection sort: ordina in base a s2, poi sz
    do i = 1, n-1
       min_index = i
       do j = i+1, n
          if (s2_local(j) < s2_local(min_index) .or. &
               (abs(s2_local(j) - s2_local(min_index)) < threshold .and. &
               sz_local(j) < sz_local(min_index)) .or. &
               (abs(s2_local(j) - s2_local(min_index)) < threshold .and. &
               abs(sz_local(j) - sz_local(min_index)) < threshold .and. &
               eigenvalues_local(j) < eigenvalues_local(min_index))) then
             min_index = j
          end if
       end do

       ! Scambio dei valori solo se necessario
       if (min_index /= i) then
          ! Scambio s2
          temp_s2 = s2_local(i)
          s2_local(i) = s2_local(min_index)
          s2_local(min_index) = temp_s2

          ! Scambio sz
          temp_sz = sz_local(i)
          sz_local(i) = sz_local(min_index)
          sz_local(min_index) = temp_sz

          ! Scambio autovalori
          temp_eigenvalue = eigenvalues_local(i)
          eigenvalues_local(i) = eigenvalues_local(min_index)
          eigenvalues_local(min_index) = temp_eigenvalue

          ! Scambio posizioni originali
          temp_pos = old_positions(i)
          old_positions(i) = old_positions(min_index)
          old_positions(min_index) = temp_pos
       end if
    end do

    ! Copia gli autovalori ordinati nell'array di output
    sorted_values = eigenvalues_local

  end subroutine sort_eigenvalues

  subroutine approx_vb(dim, hamiltonian)
    implicit none
    integer, intent(in):: dim
    complex*16, intent(inout):: hamiltonian(dim,dim)

    integer::i, j, n
    complex*16, allocatable:: v1(:), v2(:)
    complex*16:: dot_product
    real*8:: temp
    logical:: bool

    do i = 1,dim
       do j =1, dim
          if(zabs(hamiltonian(i,j)).le.1d-4)hamiltonian(i,j) = 0d0
       end do
    end do
    allocate(v1(dim), v2(dim))
    do i = 1, dim
       temp=0d0
       do j=1,dim
          temp=temp+ dconjg(hamiltonian(j,i))*hamiltonian(j,i)
       end do
       do j = 1, dim
          hamiltonian(j,i) = hamiltonian(j,i)/dsqrt(temp)
       end do
    end do

    bool = .true.

    ! Loop su ogni coppia di colonne (i, j)
    do i = 1, dim
       ! Controlla la normalizzazione del vettore i (colonna i)
       dot_product=0d0
       v1=hamiltonian(:, i)
       do n = 1, dim
          dot_product = dot_product + conjg(v1(n)) * v1(n)
       end do
       temp = zabs(dot_product) - 1d0
       if (dabs(temp) .gt. 1d-14) then
          bool = .false.
       endif

       ! Controlla ortogonalità rispetto agli altri vettori
       do j = i + 1, dim
          dot_product=0d0
          v1=hamiltonian(:, i)
          v2=hamiltonian(:, j)

          do n = 1, dim
             dot_product = dot_product + conjg(v1(n)) * v2(n)
          end do
          temp = zabs(dot_product)
          if (dabs(temp) .gt. 1d-10)then
             bool = .false.
          endif
       end do
    end do
    if(bool)write(*,*) 'Ortonormal'

  end subroutine approx_vb

    subroutine assign_labels(dim2, nsiti, charge, label)
        implicit none
        integer, intent(in) :: dim2, nsiti
        real(8), intent(in) :: charge(dim2, nsiti)
        character(len=10), intent(out) :: label(dim2)
        integer :: i
        real(8), dimension(nsiti) :: state

        do i = 1, dim2
            state = charge(i, :)
            if (all(dabs(state) < 1.0e-4)) then
                label(i) = 'N'

            else if (all(dabs(state - (/1.0d0, 0.0d0, 0.0d0, -1.0d0/)) < 1.0e-4)) then
                label(i) = 'LRCT'

            else if (all(dabs(state - (/0.0d0, -1.0d0, +1.0d0, 0.0d0/)) < 1.0e-4)) then
                label(i) = 'B1B2'

            else if (all(dabs(state - (/0.0d0, +1.0d0, -1.0d0, 0.0d0/)) < 1.0e-4)) then
                label(i) = 'B1B2'

            else if (all(dabs(state - (/1.0d0, 0.0d0, -1.0d0, -0.0d0/)) < 1.0e-4)) then
                label(i) = 'DB2'

            else if (all(dabs(state - (/1.0d0, -1.0d0, 0.0d0, 0.0d0/)) < 1.0e-4)) then
                label(i) = 'DB1'

            else if (all(dabs(state - (/0.0d0, 1.0d0, 0.0d0, -1.0d0/)) < 1.0e-4)) then
                label(i) = 'B1A'

            else if (all(dabs(state - (/0.0d0, 0.0d0, 1.d0, -1.0d0/)) < 1.0e-4)) then
                label(i) = 'B2A'

            else if (all(dabs(state - (/1.0d0, 1.0d0, -1.0d0, -1.0d0/)) < 1.0e-4)) then
                label(i) = 'DB-BA'

            else if (all(dabs(state - (/1.0d0, -1.0d0, 1.0d0, -1.0d0/)) < 1.0e-4)) then
                label(i) = 'CI'

            else if (all(dabs(state - (/0.0d0, -1.0d0, 1.0d0, 0.0d0/)) < 1.0e-4)) then
                label(i) = 'B1B2'

            else if (all(dabs(state - (/2.0d0, 0.0d0, 0.0d0, -2.0d0/)) < 1.0e-4)) then
                label(i) = 'DLRCT'
            else
                label(i) = 'ALTRO'
            end if
        end do
    end subroutine assign_labels


    subroutine charge_num(dim2, nsiti, charge, label)
        implicit none
        integer, intent(in) :: dim2, nsiti
        real(8), intent(in) :: charge(dim2, nsiti)
        real(8), intent(out) :: label(dim2)
        integer :: i
        real(8), dimension(nsiti) :: state

        do i = 1, dim2
            state = charge(i, :)
            if (all(dabs(state) < 1.0e-4)) then
                label(i) = 0d0

            else if (all(dabs(state - (/1.0d0, 0.0d0, 0.0d0, -1.0d0/)) < 1.0e-4)) then
                label(i) = 3d0

            else if (all(dabs(state - (/1.0d0, 0.0d0, -1.0d0, -0.0d0/)) < 1.0e-4)) then
                label(i) = 2d0

            else if (all(dabs(state - (/1.0d0, -1.0d0, 0.0d0, 0.0d0/)) < 1.0e-4)) then
                label(i) = 1d0

            else
                label(i) = 5d0
            end if
        end do
      end subroutine charge_num

      subroutine write_j(nsiti, dim, hamiltonian, basis, corr)
        implicit none
        integer, intent(in):: nsiti, dim, basis(dim)
        complex*16, intent(in):: hamiltonian(dim, dim)
        complex*16, intent(out):: corr(nsiti-1, dim, dim)
        complex*16::    conj_corr(nsiti-1, dim, dim)

        integer:: i, j, k, n, nso
        complex*16::number1(dim ,dim), temp1(dim, dim), temp2(dim, dim), number2(dim, dim), number3(dim, dim), j12(dim,dim), j23(dim, dim), j34(dim,dim)
        complex*16:: cplx
        logical:: bool
        cplx = cmplx(0.d0, 1.d0)
        nso = nsiti*2
        number1=0
        number2=0
        number3=0
        do i = 1, dim
           if(btest(basis(i), 0)) number1(i, i) = number1( i, i) +1
           if(btest(basis(i),1)) number1( i, i) = number1( i, i) +1
           if(btest(basis(i), 2)) number2(i, i) = number2( i, i) +1
           if(btest(basis(i),3)) number2( i, i) = number2( i, i) +1
           if(btest(basis(i), 4)) number3(i, i) = number3( i, i) +1
           if(btest(basis(i),5)) number3( i, i) = number3( i, i) +1
        enddo
        temp1=0d0
        temp2=0d0
        do i = 1, dim
           do j = 1, dim
              do k = 1, dim
                 temp1(i,j) = temp1(i,j) + number1(i,k) * hamiltonian(k,j)
                 temp2(i,j) = temp2(i,j) + hamiltonian(i,k) * number1(k,j)
              end do
           end do
        end do
        
        
        j12= cplx*(temp1-temp2)

        temp1=0d0
        temp2=0d0
        do i = 1, dim
           do j = 1, dim
              do k = 1, dim
                 temp1(i,j) = temp1(i,j) + number3(i,k) * hamiltonian(k,j)
                 temp2(i,j) = temp2(i,j) + hamiltonian(i,k) * number3(k,j)
              end do
           end do
        end do

        j34 =  cplx*(temp1-temp2)
        temp1=0d0
        temp2=0d0
        do i = 1, dim
           do j = 1, dim
              do k = 1, dim
                 temp1(i,j) = temp1(i,j) + number2(i,k) * hamiltonian(k,j)
                 temp2(i,j) = temp2(i,j) + hamiltonian(i,k) * number2(k,j)
              end do
           end do
        end do
        j23 = cplx*(temp1-temp2)

       
        j23 = j23 + j12
        j34 = j34 + j23

        call check_hermitian(j12, dim, bool)
        if(bool) write(*,*) 'J12 good'
        call check_hermitian(j23, dim, bool)
        if(bool) write(*,*) 'J23 good'
        call check_hermitian(j34, dim, bool)
        if(bool) write(*,*) 'J34 good'
        corr(1,:,:) = j12
        corr(2, :,:) = j23
        corr(3, :,:) = j34

        open(1, file='j12.dat')
        open(2, file = 'j23.dat')
        open(3, file = 'j34.dat')
        do i = 1, dim
           do j = 1, dim
              write(1, '(<2>(I3, 2x), <2>(f10.5,2x))') i, j, dreal(j12(i,j)), dimag(j12(i,j))
              write(2, '(<2>(I3, 2x), <2>(f10.5,2x))') i, j, dreal(j23(i,j)), dimag(j23(i,j))
              write(3, '(<2>(I3, 2x), <2>(f10.5,2x))') i, j, dreal(j34(i,j)), dimag(j34(i,j))

           enddo
        enddo
      end subroutine write_j


!!$      subroutine OptimizeCurrent(J_goal,H, lambda, base, t, n_sites,dim, dim_lambda,whatstate)
!!$        integer, intent(in) :: n_sites,dim, dim_lambda
!!$        integer, intent(in) :: base(dim)
!!$        doublecomplex, intent(in) :: H(dim, dim)
!!$        doublecomplex, intent(inout) :: J_value(dim_lambda)
!!$        doubleprecision, intent(in) :: J_goal
!!$        doubleprecision, intent(in) :: t(n_sites,n_sites)
!!$        doubleprecision, intent(inout) :: lambda(dim_lambda)
!!$        integer, intent(in) :: whatstate
!!$        doubleprecision xii
!!$
!!$        integer i,j, count_opt, count, info, n
!!$        doubleprecision err_j, wastedouble
!!$        doubleprecision, dimension(:), allocatable :: lambdaones, vectorwaste2
!!$        doubleprecision, allocatable :: w_0(:), rwork(:), t_temp(:,:)
!!$        doublecomplex, dimension(:,:), allocatable :: H_lambda, trash(:,:)
!!$        doublecomplex, dimension(:,:,:), allocatable :: j_current
!!$        doublecomplex,allocatable :: work(:)
!!$       
!!$
!!$        xii = 0.d0
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$        !Lambda optimization!
!!$        !Let's use some good old SCF magic!
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$        allocate(lambdaones(dim_lambda), vectorwaste2(dim_lambda), t_temp(dim_lambda+1, dim_lambda+1))
!!$        allocate(j_current(n_sites, dim, dim), H_lambda(dim, dim))
!!$        if(abs(J_goal).gt.1.d-10)then
!!$           wastedouble = 111.d0
!!$           do i=1, dim_lambda
!!$11            continue
!!$              lambda(i) = 1d0
!!$              if(abs(lambda(i)).lt.1.d-11) go to 11
!!$           end do
!!$           write(*,*) 'Starting lambdas:'
!!$           write(*,'(<dim_lambda>(F7.2,X))') (lambda(i), i=1, dim_lambda)
!!$           lambdaones = 1.d0
!!$           err_j = 1.d0
!!$           count_opt= 0
!!$           call  write_j(n_sites, dim, h, base, j_current)
!!$           do while (err_j.gt.1.d-4.and.count_opt.lt.200)
!!$
!!$              H_lambda = H
!!$              do i = 1, dim_lambda


!!$                 do j = 1, dim
!!$                    do n = 1, dim
!!$                       H_lambda(j,n) = H_lambda(j,n) - j_current(i, j,n) *lambda(i)
!!$                    enddo
!!$                 enddo
!!$              enddo
!!$              allocate(w_0(dim), rwork(3*dim-2), work(2*dim))
!!$              call zheev('V', 'U', dim, H_lambda, dim, w_0, work, 3*dim, rwork, info)
!!$              deallocate(w_0, rwork, work)
!!$              allocate(trash(dim,dim))
!!$              do i=1, dim_lambda
!!$                 trash=0d0
!!$                 trash=j_current(i, :,:)
!!$                 j_value(i) = exp_value_complex(H_lambda(:,whatstate),trash,  H_lambda(:,whatstate), dim)!*1.6d-19/hbar_eV
!!$              end do
!!$              deallocate(trash)
!!$              count_opt = count_opt + 1
!!$              write(*,'(a, I3)') 'ITER: ', count_opt
!!$              err_j = dabs((dreal(j_value(1)) - J_goal)/J_goal)
!!$              count = 1
!!$              vectorwaste2(1) = err_j
!!$              do i = 2, dim_lambda
!!$                 vectorwaste2(i) = dabs((dreal(j_value(i)) - J_goal)/J_goal)
!!$                 if(err_j.lt.vectorwaste2(i))then
!!$                    err_j = vectorwaste2(i)
!!$                    count = i
!!$                 end if
!!$              end do
!!$              write(*,'(a, <dim_lambda>(e9.1, 2x),a,e8.1)') 'ERRS', (vectorwaste2(j), j=1, dim_lambda),'ERR_MAX: ', err_j
!!$              write(*,'(a, e9.1)') 'Target:', J_goal
!!$              write(*,'(a, <dim_lambda>(e9.1, 2x))') 'Partial:', (real(J_value(j)), j=1, dim_lambda)
!!$              write(*,'(a, <dim_lambda>(e9.1, 2x))') 'lambdas:', (lambda(j), j=1, dim_lambda)
!!$              write(*,*)
!!$              if(err_j.gt.1.d-4)then
!!$                 do i=1, dim_lambda
!!$                    lambda(i) = lambda(i) * j_goal/j_value(i)
!!$                    !lambda(i) = (lambda(i) * j_goal/j_value(i) + lambda(i))/2
!!$                 end do
!!$                 !              if(mod(count_opt,25).eq.0)then
!!$                 !                 write(*,*) "Urgh... what a difficult convergence...let's displace!"
!!$                 !                 lambda = lambda+1
!!$                 !              end if
!!$              end if
!!$           end do
!!$        else
!!$           lambda = 0.d0
!!$           count_opt= 0
!!$        end if
!!$        write(*,'(a, i3, a, e9.2)') 'Optimization ended at attempt: ', count_opt,' Err: ', err_j
!!$
!!$      end subroutine OptimizeCurrent


      function exp_value_complex(state1,hamiltonian,state2, dim)result(fin_value)
        integer, intent(in) :: dim
        doublecomplex, intent(in) :: state1(dim), hamiltonian(dim, dim), state2(dim)
        doublecomplex:: bra(1, dim), temp(1,dim), fin_value

        integer:: i, j, k

        temp=0d0
        do j = 1, dim
           do k = 1, dim
              temp(1,j) = temp(1,j) +  dconjg(state1(k)) * hamiltonian(k,j)
           end do
        end do

        fin_value = 0d0
        do k = 1, dim
           fin_value = fin_value + temp(1,k) * state2(k)
        end do
             
      endfunction exp_value_complex
      !======================================================================
      subroutine OptimizeCurrent2(J_goal,H, lambda, basis, nsiti,dim)
        integer, intent(in) :: nsiti,dim
        integer, intent(in) :: basis(dim)
        doublecomplex, intent(in) :: H(dim, dim)
        doubleprecision, intent(in) :: J_goal
        doubleprecision, intent(inout) :: lambda(nsiti-1)
        doubleprecision xii

        integer i,j, count_opt, info, n, dim_lambda, delta_j
        doubleprecision err_j, wastedouble, learning_rate(nsiti-1)
        doubleprecision, dimension(:), allocatable ::  vectorwaste, prev_error
        doubleprecision, allocatable :: w_0(:), rwork(:), t_temp(:,:)
        doublecomplex, dimension(:,:), allocatable :: H_lambda, trash(:,:), trash2(:,:)
        doublecomplex, dimension(:,:,:), allocatable :: j_current
        doublecomplex,allocatable :: work(:)
        doublecomplex:: j_value(nsiti-1), state1(dim), bb

        dim_lambda = nsiti-1
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !Lambda optimization!
        !Let's use some good old SCF magic!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        allocate( vectorwaste(dim_lambda), prev_error(dim_lambda))
        allocate(j_current(nsiti, dim, dim), H_lambda(dim, dim))
        allocate( trash(dim, dim))
        allocate(w_0(dim), rwork(3*dim-2), work(2*dim))
                   
        err_j = 1d0
        count_opt = 0
        do i = 1, dim_lambda
           lambda(i) = 1.d0
        enddo
        call write_j(nsiti, dim, h, basis, j_current)
        do while ((err_j.ge.1d-4).and.(count_opt.le.200))
           h_lambda = 0d0
           do n = 1, dim_lambda
              do i = 1, dim
                 do j = 1, dim
                    h_lambda(i,j) = h(i,j) - lambda(n) * j_current(n, i, j)
                 end do
              end do
           end do
           call zheev('V', 'U', dim, H_lambda, dim, w_0, work, 3*dim, rwork, info)

           do n = 1, dim_lambda
              trash = j_current(n, :, :)
              j_value(n) = exp_value_complex(h_lambda(:,1),trash,h_lambda(:,1), dim)
              vectorwaste(n) = dreal(j_value(n) - j_goal)
           enddo
           
           err_j = 0d0
           do  n = 1, dim_lambda
              if(vectorwaste(n).ge.err_j) err_j = vectorwaste(n)
           enddo

           write(*, '(a, 2x, I4)') 'Cicle =', count_opt
           write(*,'(a, <dim_lambda>(e9.1, 2x),a,e8.1)') 'ERRS', (vectorwaste(j), j=1, dim_lambda),'ERR_MAX: ', err_j
           write(*,'(a, e9.1)') 'Target:', J_goal
           write(*,'(a, <dim_lambda>(e9.1, 2x))') 'Partial:', (real(J_value(j)), j=1, dim_lambda)
           write(*,'(a, <dim_lambda>(e9.1, 2x))') 'lambdas:', (lambda(j), j=1, dim_lambda)
           write(*,*)

           
           if(err_j.ge.1d-4)then
              count_opt = count_opt+1
              do n = 1, dim_lambda
                 lambda(n) = lambda(n)  - vectorwaste(n)
              enddo
           endif
           
        enddo
        write(*,'(a, i3, a, e9.2)') 'Optimization ended at attempt: ', count_opt,' Err: ', err_j
      end subroutine OptimizeCurrent2



      
    end module
