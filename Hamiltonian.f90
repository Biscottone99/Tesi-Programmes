program bagaglio
  use module
  implicit none

  integer :: i, n, j, nsiti, nso, dim, lwork, lrwork, liwork, info, isito,jsito, npoints, nstate, max, min,m, k
  integer, allocatable :: basis(:), nz(:), iwork(:),  basism2(:), basism(:), basisz(:), basisp(:), basisp2(:),oldpositions(:), newposition(:), old(:), new(:)
  complex*16, allocatable :: hamiltonian(:,:), coup_mono(:,:), soc(:,:), sso(:,:), soo(:,:), work(:), dummy(:,:), mu(:,:,:), sq(:,:),sqrot(:,:), mono_rot(:,:), sso_rot(:,:), soo_rot(:,:), soc_rot(:,:)
  complex*16 :: imag, pf
  real*8, allocatable :: coord(:,:), u(:), esite(:), r(:,:), hop_so(:,:), pot(:), hop(:,:), eigenvalue(:), rwork(:), hop_use(:,:), carica(:,:), dipole(:,:), temporary(:), charges(:,:), spectral(:,:)
  real*8,allocatable:: num(:,:), num_rot(:,:), eigm2(:), eigm(:), eigz(:), eigp(:), eigp2(:), ss22(:), copia(:), sorted(:), numa(:), numd(:), numar(:,:), numdr(:,:), eiggg(:)
  real*8 :: length, t, dx, dy, dz, me, gs, e, e0, pi, cl, jgoal
  logical :: PPPflag, hubbardflag, SOCflag, bool, redflag
  complex*16,allocatable::hamm2(:,:), hamm(:,:), hamz(:,:), hamp(:,:), hamp2(:,:), coupx(:,:), coupy(:,:), coupz(:,:), soox(:,:), sooy(:,:), sooz(:,:), ssox(:,:), ssoy(:,:), ssoz(:,:), corr(:,:,:), number(:,:,:)
  doublecomplex, allocatable:: coupx_r(:,:), coupy_r(:,:), coupz_r(:,:), soox_r(:,:), sooy_r(:,:), sooz_r(:,:), ssox_r(:,:), ssoy_r(:,:), ssoz_r(:,:), socx_r(:,:), socy_r(:,:), socz_r(:,:), mat_out(:,:)
  character(1) :: jobz, uplo, unit
  character(2), allocatable :: state(:)
  character*1,allocatable::molt(:)


  real*8:: norm, temp

  complex*16, allocatable:: psi0(:),  ham_red(:,:), temp_array(:,:), acc(:,:), don(:,:), rho(:,:), sd(:,:), sdr(:,:), trasha(:,:), trashb(:,:), trashc(:,:), sz(:), bubu(:), temp1(:,:), temp2(:,:)
  complex*16, allocatable:: j_value(:), szrot(:)
  real*8,allocatable:: lambda(:)

  !=========================CONSTANTS=========================
  imag = cmplx(0.0, 1.0)
  me = 9.1093837015d-31
  gs = 2.00231930436256
  e = 1.602176634d-19
  e0 = 8.8541878128d-12
  pi = dacos(-1.0d0)
  cl = 299792458
  pf = ((gs * e**2) / (8 * pi * e0 * me * cl**2)) * 10.0d10
  write(*,*) pf
  npoints=1d5
  !===========================================================

  ! Reading input data
  open(1,file='dim2.dat')
  read(1,*) dim
  close(1)
  open(1, file='input-red.inp')
  read(1,*) nsiti
  read(1,*) length ! Armstrong
  read(1,*) t ! eV
  read(1,*) PPPflag
  read(1,*) hubbardflag
  read(1,*) SOCflag
  read(1,*) redflag
  read(1,*) nstate
  read(1,*) jgoal
  nso = 2 * nsiti
  allocate(coord(nsiti, 3), basis(dim), u(nsiti), esite(nsiti), r(nsiti, nsiti), nz(nsiti))
  do i = 1, nsiti
      read(1,*) u(i), esite(i), nz(i)
  end do
  close(1)
  open(2,file='geom.dat')
  do i = 1, nsiti
      read(2,*) coord(i, 1), coord(i, 2), coord(i, 3)
  end do
  close(2)

  open(11,file='output.out')
  write(1,*) 'COORDINATES'
  do i=1,nsiti
     write(11, '(I2, 2x, <3>(f10.5))') i, coord(i, 1), coord(i, 2), coord(i, 3)
  enddo
  write(11,'(A, 2x, f10.5)')'Bond lenght (angstrom) =', length
  WRITE(11,*)
  write(11,*) 'ELECTRONIC PARAMETERS (U, site energy, Z)'
  do i =1, nsiti
     write(11,'(I2,2x,<2>(f10.5,2x), I2)') i, u(i), esite(i), nz(i)
  enddo
  write(11,'(A,2x,f10.5)') 'Hopping t=', t
  WRITE(11,*)
  write(11,*) 'OTHER PARAMETERS'
  write(11,*) 'PPP', pppflag
  write(11,*) 'Hubbard', Hubbardflag
  write(11,*) 'SOC', socflag
  write(11,'(A,2x,f10.5)') 'Current', t/jgoal
  !=========================BASIS CREATION=========================

  open(2,file='basis.dat')
  do i =1, dim
      read(2,*) basis(i)
      !write(*,*) basis(i)
  enddo

  !============================================================


  allocate(sq(dim,dim), charges(dim, nsiti), carica(dim, nsiti), sz(dim))
  call  s2_realspace(dim, nso, basis, sz, sq)
  call charge(carica, basis, nz, dim, nso)

  ! Computing the distance matrix
  do i = 1, nsiti
      do j = 1, nsiti
          dx = coord(i, 1) - coord(j, 1)
          dy = coord(i, 2) - coord(j, 2)
          dz = coord(i, 3) - coord(j, 3)
          r(i, j) = dsqrt(dx**2 + dy**2 + dz**2)
      end do
  end do
  ! computing the nn-hopping matrix
  allocate(hop_use(nsiti, nsiti))
  do i = 1, nsiti
      do j = 1, nsiti
          if (dabs(r(i, j) - length) < 1.0d-8) hop_use(i, j) = t
      end do
  end do

  allocate(hop_so(nso, nso))
  hop_so = 0.0d0
  call op_siti_2_so(hop_so, hop_use, nso, nsiti)

  ! Converting from tight binding to real space basis
  !============================================================
  allocate( eigenvalue(dim),hamiltonian(dim,dim),pot(dim),hop(dim,dim), soc(dim,dim))
  if(hubbardflag)call site_energy_u(nso, dim, esite, u, basis, pot)
  if(PPPflag) call ppp_diag(dim, nsiti, coord, esite, basis, u, nz, pot)
  call sq_oe_op_real(nso,dim,hop_so,hop,basis)
  soc=0
  if (SOCflag) then
      allocate(coup_mono(dim, dim), sso(dim, dim), soo(dim, dim))
      allocate( trasha(dim,dim), trashb(dim,dim), trashc(dim,dim))
      call compute_soc_mono(nsiti, dim, nz, coord, hop_use, basis, pf, coup_mono, trasha, trashb, trashc)
      call compute_sso(nsiti, dim, coord, hop_use, basis, pf, sso, trasha, trashb, trashc)
      call compute_soo(nsiti, dim, coord, hop_use, basis, pf, soo, trasha, trashb, trashc)
      deallocate(trasha, trashb, trashc)
!!$      sso = 0d0
!!$      soo = 0d0
      do i = 1, dim
          do j = 1, dim
              soc(i, j) = coup_mono(i, j) - sso(i, j) - soo(i, j)
          end do
      end do
  end if

  hamiltonian=0d0
  do i = 1, dim
      hamiltonian(i,i) = hamiltonian(i,i) + pot(i)
      do j =1, dim
          hamiltonian(i,j) = hamiltonian(i,j) - hop(i,j) + soc(i,j)
      end do
  end do

  !============================================================

  allocate(num(dim,nso))
  allocate(numa(dim), numd(dim), sd(dim, nsiti))
  num=0
  numa=0d0
  numd=0d0
  do i = 1, dim
      do j = 0, nso-1
          if(btest(basis(i),j))num(i,j+1)=num(i,j+1)+1
      enddo
      numa(i)=num(i,nso)+num(i,nso-1)
      numd(i)=num(i,1)+num(i,2)
      do j = 1, nso-1, 2
          isito = (j+1)/2
          sd(i, isito)=num(i, j) - num(i,j+1)
      end do
  enddo

 
  allocate(lambda(nsiti-1), j_value(nsiti-1))
  allocate(corr(nsiti-1, dim, dim))
  call write_j(nsiti, dim, hamiltonian, basis, corr)
  !call OptimizeCurrent(t/1000, J_value,hamiltonian, lambda, basis, hop_use, nsiti,dim, nsiti-1,1)
  call OptimizeCurrent2(t/jgoal,Hamiltonian, lambda, basis, nsiti,dim)
  do i = 1, dim
     do j =1, dim
        do n = 1, nsiti-1
           hamiltonian(i,j) = hamiltonian(i,j) -  lambda(n)*corr(n, i, j)
        enddo
     enddo
  enddo
  !=========================Diagonalizarion========================
  lrwork = (1 + 5*dim + 2*dim**2)
  liwork = (3 + 5*dim)
  lwork = (2*dim + dim**2)
  allocate(work(max(1, lwork)), rwork(lrwork), iwork(max(1, liwork)),molt(dim), sqrot(dim, dim), szrot(dim))
  call  check_hermitian(hamiltonian, dim, bool)
  if(.not.bool)write(*,*) 'Hamiltonian error'
  call zheevd('V', 'U', dim, hamiltonian, dim, eigenvalue, work, lwork, rwork, lrwork, iwork, liwork, info)
  temp = eigenvalue(1)
  eigenvalue = eigenvalue -temp

  szrot = 0.d0
 
  call eigenvalues(dim, 1d-8, eigenvalue, molt)
  call rotate_cplx_2x2(dim, sqrot, sq, hamiltonian)

  do i = 1, dim
     do j = 1, dim
        szrot(i) = szrot(i) + dconjg(hamiltonian(j, i)) * hamiltonian(j, i) * sz(j)
     end do
  end do

  allocate(num_rot(dim,nso), numar(dim,dim), numdr(dim,dim), sdr(dim,nsiti))
  call rot_diag(dim, charges, carica, nsiti, hamiltonian)
  call rot_diag(dim, num_rot, num, nso, hamiltonian)
  call rotate_real_1x2(dim, numar, numa, hamiltonian)
  call rotate_real_1x2(dim, numdr, numd, hamiltonian)
  call rot_diag_cplx(dim, sdr, sd, nsiti, hamiltonian)

  deallocate(num, numd, numa)
  write(11,*) 'EIGENVALUES,  CHARGES, STATE'
  do i = 1, dim
     write(11, '(I5, 2x, <nsiti+1>(f10.5,2x), A)') i, eigenvalue(i), (charges(i, j), j=1,nsiti), molt(i)
  enddo
  write(1,*)
   write(11,*) 'EIGENVALUES, S^2, SZ, SPIN DENSITY'
  do i = 1, dim
     write(11, '(I3,2x,<3>(f10.5, 2x), 2x,<2>(f10.5,2x), f15.12)') i, eigenvalue(i), dreal(sqrot(i,i)), dreal(szrot(i)), dreal(sdr(i,1)), dreal(sdr(i,nsiti)), dreal(sdr(i,1)-sdr(i,nsiti))
  enddo

  !=========================REDFIELD=========================
  if(redflag)then
     allocate(psi0(dim))
     psi0=0

     do i=1,dim
        psi0(i)=(1/dsqrt(2.d0))*(hamiltonian(54,i)-hamiltonian(35,i))
     enddo
     norm=0.d0
     do i=1,dim
        norm=norm+dconjg(psi0(i))*psi0(i)
     enddo
     psi0=psi0/dsqrt(norm)

     do i=1, dim
        write(150,*) i, zabs(psi0(i))**2
     enddo


     allocate(acc(nstate,nstate), don(nstate,nstate))
     acc=0
     don=0
     do i =1, nstate
        do j=1, nstate
           write(123, *) i,j, numar(i,j)
           acc(i,j)=numar(i,j)
           don(i,j)=numdr(i,j)
        enddo
     enddo
     allocate(temp_array(nstate, nso))
     temp_array=0d0
     do i =1, nstate
        do j =1, nso
           temp_array(i,j)= num_rot(i,j)
        end do
     end do
     !===================INPUT REDFIELD=================
     open(99,file='input-red/num.bin',form="unformatted")
     write(99) temp_array
     close(99)
     deallocate(temp_array)
     allocate(temp_array(nstate, nsiti))
     do i =1, nstate
        do j = 1, nsiti
           temp_array(i,j)=sdr(i,j)
        end do
     end do
     open(77,file='input-red/spin-density.bin',form="unformatted")
     write(77) temp_array(1:nstate, 1:nsiti)
     close(77)

     open(77,file='input-red/numa.bin',form="unformatted")
     write(77) acc(1:nstate, 1:nstate)
     close(77)
     open(77,file='input-red/numd.bin',form="unformatted")
     write(77) don(1:nstate, 1:nstate)
     close(77)
     allocate(bubu(nstate), eiggg(nstate))
     eiggg=0d0
     do i =1, nstate
        eiggg(i)=eigenvalue(i)
     end do

     open(66,file='input-red/eigen.bin',form="unformatted")
     write(66) eiggg(1:nstate)
     close(66)
     bubu=0d0
     do i =1, nstate
        bubu(i)=psi0(i)
     end do

     open(88,file='input-red/psi0.bin',form="unformatted")
     write(88) bubu(1:nstate)
     close(55)
     close(56)
     close(57)
     close(88)

     open(55,file='input-red/system_input.dat')
     write(55,*) nstate
     write(55,*) nsiti
  endif

end program bagaglio
