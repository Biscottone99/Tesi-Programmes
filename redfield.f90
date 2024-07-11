 program red2
    implicit none
    integer a,b,c,d,e,i,j,k,ii,jj,kk,iii,jjj,kkk,dim,conta,info,&
            lwork,nt,dimsq,ab,cd,nkrylov,it, nsiti, isiti, nso
    double precision pi,temp,hbar,kb,kbt,norm,gamma,dt,time, eta, cutoff, pstreshold, base
    double precision, allocatable :: w(:),work(:),spectral(:,:),bedistr(:,:), sdr(:,:), num(:,:)
    complex*16:: imag,trace,overlap,entime, muxtime, muytime, muztime, omega
    complex*16, allocatable :: psi0(:),denmat(:,:),redfield(:,:,:,:),denvec(:),&
         lmat(:,:), mux(:,:), muy(:,:), muz(:,:), sdrtime(:), numtime(:),numa(:,:), numd(:,:)
    logical:: psflag, blflag

    !!! global parameters
    hbar = 6.582119569d-16 !eV*s
    hbar = hbar * 1e15 !eV*fs
    kb = 8.61733326e-5 !eV/K
    imag = (0.d0,1.d0)
    base=exp(1.0)
    pi = dacos(-1.d0)
   
    ! read global parameters
    open(1,file='global_input.inp')
    read(1,*) temp
    read(1,*)eta! gamma
    read(1,*) cutoff
    read(1,*) dt
    read(1,*) nt
    read(1,*) nkrylov
    read(1,*) psflag
    read(1,*) pstreshold
    read(1,*) blflag
    close(1)

    kbt = kb*temp !eV

    !!! end global parameters

    !!! system parameters

    ! read input
    open(1,file='../input-red/system_input.dat')
    read(1,*) dim
    write(*,*) dim
   ! dim=dim-1
    read(1,*) nsiti
    close(1)
    nso=nsiti*2
    !!! end system parameters

    !!! reading the operators
    allocate(w(dim), mux(dim,dim), muy(dim,dim), muz(dim,dim),sdr(dim,nsiti), psi0(dim), numa(dim,dim), numd(dim,dim), num(dim,nso))
    open(1, file='../input-red/mux.bin', form="unformatted")
    read(1) mux
    close(1)
    open(1, file='../input-red/muy.bin', form="unformatted")
    read(1) muy
    close(1)
    open(1, file='../input-red/muz.bin', form="unformatted")
    read(1) muz
    close(1)
    open(1, file='../input-red/numa.bin', form="unformatted")
    read(1) numa
    close(1)
    open(1, file='../input-red/numd.bin', form="unformatted")
    read(1) numd
    close(1)
    open(1, file='../input-red/eigen.bin', form="unformatted")
    read(1) w
    close(1)
    open(1, file='../input-red/psi0.bin', form="unformatted")
    read(1) psi0
    close(1)
    open(1, file='../input-red/spin-density.bin', form="unformatted")
    read(1) sdr
    close(1)
    open(1, file='../input-red/num.bin', form="unformatted")
    read(1) num
    close(1)

    
    !!! end reading operators


    !!! preparing the initial state
    !normalization
    norm = 0.d0
    do i = 1,dim
        norm = norm + dconjg(psi0(i))*psi0(i)
    end do
    norm = dsqrt(norm)
    psi0 = psi0/norm
    !writing the density matrix at t=0
    allocate(denmat(dim,dim))
    denmat = (0.d0,0.d0)
    do i = 1,dim
        do j = 1,dim
                denmat(i,j) = dconjg(psi0(i))*psi0(j)
        end do
    end do
   open(1,file='stuff/denmat.dat')
    do i = 1,dim
        write(1,*) (denmat(i,j), j = 1,dim)
    end do
    close(1)
    !!!

    !!!Bose-Einstein and spectral density
    allocate(spectral(dim,dim),bedistr(dim,dim))
    spectral = 0.d0
    bedistr = 0.d0
    do i = 1,dim
        do j = 1,dim
            if (dabs(w(i)-w(j)).gt.1.d-9) then
                !bedistr(i,j) = (dexp((w(i)-w(j))/kbt)-1.d0)**(-1)
                if (w(i).gt.w(j)) bedistr(i,j) = (dexp((w(i)-w(j))/kbt)-1.d0)**(-1)
                if (w(j).gt.w(i)) bedistr(i,j) = (dexp((w(i)-w(j))/kbt)-1.d0)**(-1)
                omega=dabs(w(i)-w(j))/cutoff
                !if (w(i).gt.w(j)) spectral(i,j) = gamma    !costante
                ! if (w(j).gt.w(i)) spectral(i,j) = -gamma  !costante
                spectral(i,j) = (hbar**2/pi)*eta* (omega/(1+omega**2))  !Debye
               ! spectral(i,j)= 2*hbar*eta*cutoff* omega/(omega**2+eta**2) !Drude-Lorentz
               ! spectral(i,j)=(eta*hbar/cutoff)* omega*base**(-omega/cutoff) !Ohmic
             end if
        end do
    end do

    open(1,file='stuff/spectral.dat')
    open(2,file='stuff/bedistr.dat')
    do i = 1,dim
        write(1,*) (spectral(i,j), j = 1,dim)
        write(2,*) (bedistr(i,j), j = 1,dim)
    end do
    close(1)
    close(2)
    !!! end Bose-Einstein and spectral density

    !!! writing redfield superoperator
    allocate(redfield(dim,dim,dim,dim))
    redfield = (0.d0,0.d0)
    do a = 1,dim
        do b = 1,dim
            do c = 1,dim
               do d = 1,dim
                  if (b==d) then
                        do e = 1,dim
                            redfield(a,b,c,d) = redfield(a,b,c,d) - &
                                    gammap (numa,numa,spectral, bedistr,w,a,e,e,c,dim)
                            redfield(a,b,c,d) = redfield(a,b,c,d) - &
                                    gammap (numd,numd,spectral, bedistr,w,a,e,e,c,dim)
                            
                        end do
                    end if

                    if (a==c) then
                        do e = 1,dim
                            redfield(a,b,c,d) = redfield(a,b,c,d) - &
                                    gammam (numa,numa,spectral, bedistr,w,d,e,e,b,dim)
                            redfield(a,b,c,d) = redfield(a,b,c,d) - &
                                    gammam (numd,numd,spectral, bedistr,w,d,e,e,b,dim)
                            
                        end do
                    end if

                    redfield(a,b,c,d) = redfield(a,b,c,d) + &
                            gammap (numa,numa,spectral, bedistr,w,d,b,a,c,dim)
                    redfield(a,b,c,d) = redfield(a,b,c,d) + &
                            gammam (numa,numa,spectral, bedistr,w,d,b,a,c,dim)
                    redfield(a,b,c,d) = redfield(a,b,c,d) + &
                            gammap (numd,numd,spectral, bedistr,w,d,b,a,c,dim)
                    redfield(a,b,c,d) = redfield(a,b,c,d) + &
                            gammam (numd,numd,spectral, bedistr,w,d,b,a,c,dim)
                end do
            end do
        end do
    end do

!!! end writing the redfield superoperator

        !!! (pseudo)-secular approximation
    if (psflag) then
        do a=1,dim
            do b = 1,dim
                do c = 1,dim
                    do d = 1,dim
                        if (dabs(w(a)-w(b)-w(c)+w(d))>pstreshold) redfield(a,b,c,d) = (0.d0, 0.d0)
                    end do
                end do
            end do
        end do
    end if
    !!! end (pseudo)-secular approximation

    !!! Bloch approximation
    if (blflag) then
        do a = 1,dim
            do b = 1,dim
                do c =1,dim
                    do d = 1,dim
                        if ((a==b.and.c==d).or.(a==c.and.b==d)) then
                            goto 1234
                        else
                            redfield(a,b,c,d) = (0.d0, 0.d0)
                        end if
                        1234 continue
                    end do
                end do
            end do
        end do
    end if
    !!! end Bloch approximation

    !!! writing denmat and redfield tensor in Liouville space
    dimsq = dim**2
    allocate(denvec(dimsq),lmat(dimsq,dimsq))
    denvec = (0.d0, 0.d0)
    lmat = (0.d0, 0.d0)
    do a = 1,dim
        do b = 1,dim
            ab = b+dim*(a-1)
            denvec(ab) = denmat(a,b)
            lmat(ab,ab) = lmat(ab,ab) - imag*(w(a)-w(b))/hbar
            do c = 1,dim
                do d  =1,dim
                    cd = d+dim*(c-1)
                    lmat(ab,cd) = lmat(ab,cd) + redfield(a,b,c,d)
!                    if (ab==cd) then
!                        lmat(ab,cd) = lmat(ab,cd) - imag*(w(a)-w(b))/hbar
!                    end if
                end do
            end do
        end do
    end do
!    do i = 1,dimsq
!        write(7777,*) (dreal(lmat(i,j)),j=1,dimsq)
!        write(6666,*) (dimag(lmat(i,j)),j=1,dimsq)
!    end do

    deallocate(redfield,denmat)
    !!! end writing denmat and redfield tensor in Liouville space

    !!! loop on time
    ! opening some trajectory files
    open(991,file='trace.dat')
    open(992,file='pop_evolution.dat')
    open(993,file='properties.dat')
    open(994,file='number.dat')
    allocate(sdrtime(nsiti),numtime(nso))
    time = 0.d0
    do it = 1,nt
        call arnoldi(denvec,lmat,dimsq,nkrylov,dt)
        !!! compute trace
        trace = (0.d0, 0.d0)
        do a = 1,dim
            trace = trace + denvec(a+dim*(a-1))
        end do
        write(991,*) time, dreal(trace), dimag(trace)
        !!! end compute trace

        !!! save populations
        write(992,*) time, (real(denvec(a+dim*(a-1))), a = 1,dim)
        !!! end save populations

        !!! compute properties
        muxtime = (0.d0, 0.d0)
        muytime = (0.d0, 0.d0)
        muztime = (0.d0, 0.d0)
        entime = (0.d0, 0.d0)
        sdrtime = (0.d0, 0.d0)
        numtime= (0.d0, 0.d0)
        !$omp parallel do default(none), &
        !$omp private(a,b,ab, isiti), &
        !$omp shared(dim,denvec,mux, muy, muz,w,sdr, nsiti, nso, num), &
        !$omp reduction(+:muxtime, muytime, muztime, sdrtime, entime, numtime)
        do a = 1,dim
            entime = entime + denvec(a+dim*(a-1))*w(a)
            do isiti=1,nsiti
                sdrtime(isiti) = sdrtime(isiti) + denvec(a+dim*(a-1))*sdr(a,isiti)
            end do
            do isiti=1,nso
               numtime(isiti) = numtime(isiti) + denvec(a+dim*(a-1))*num(a,isiti)
            end do
            do b = 1,dim
                ab = b + dim*(a-1)
                muxtime = muxtime + denvec(ab)*mux(b,a)
                muytime = muytime + denvec(ab)*muy(b,a)
                muztime = muztime + denvec(ab)*muz(b,a)
            end do
        enddo
        !$omp end parallel do
        write(993,*) time,real(entime), real(muxtime), real(muytime), real(muztime), (real(sdrtime(isiti)), isiti=1,nsiti)
        write(994,*) time,(real(numtime(isiti)), isiti=1,nso)
        !!! end compute properties
        time = time + dt
    end do
    close(991)
    close(992)
    close(993)
    close(994)
    !!! end loop on time

   
contains

!!! GAMMA^+_{ab,cd} TERM OF THE REDFIELD TENSOR
!!! OP1 and OP2 are dimensionless system operators
!!! SPECTRAL_DENSITY is the spectral density matrix in (s^-1)
!!! BOSE_EINSTEIN is the dimensionless BE distribution matrix
!!! W is the eigenvalue vector
!!! A,B,C,D are state indexes
!!! DIM is the dimension of the system
double precision function gammap (op1,op2,spectral_density,bose_einstein,w,a,b,c,d,dim)
    implicit none
    integer a,b,c,d,dim
    double precision  spectral_density(dim,dim),&
            bose_einstein(dim,dim),w(dim)
    doublecomplex op1(dim,dim), op2(dim,dim)

    gammap = 0.d0
    if (w(d).gt.w(c)) then
        gammap = op1(a,b)*op2(c,d)*spectral_density(d,c)* (bose_einstein(d,c)+1.d0)
    elseif (w(c).gt.w(d)) then
        gammap = op2(a,b)*op1(c,d)*spectral_density(c,d)* bose_einstein(c,d)
    end if
end function gammap

!!! GAMMA^-_{ab,cd} TERM OF THE REDFIELD TENSOR
!!! see above for GAMMAP
double precision function gammam (op1,op2,spectral_density,bose_einstein,w,a,b,c,d,dim)
    implicit none
    integer a,b,c,d,dim
    double precision  spectral_density(dim,dim),&
            bose_einstein(dim,dim),w(dim)
    doublecomplex op1(dim,dim), op2(dim,dim)

    gammam = 0.d0
    if (w(a).gt.w(b)) then
        gammam = op1(a,b)*op2(c,d)*spectral_density(a,b)* (bose_einstein(a,b)+1.d0)
    elseif (w(b).gt.w(a)) then
        gammam = op2(a,b)*op1(c,d)*spectral_density(b,a)* bose_einstein(b,a)
    end if
end function gammam


!!! SHORT-ITERATIVE ARNOLDI ALGORITHM
!!! DENVEC: on input, rho(t) in vector form, on output, rho(t+dt) in vector form
!!! REDFIELD: the Redfield tensor in matrix form
subroutine arnoldi(denvec,redfield,dimsq,nkrylov,dt)
    implicit none
    integer, intent(in) :: dimsq,nkrylov
    double precision, intent(in) :: dt
    doublecomplex, intent(inout), dimension(:) :: denvec
    doublecomplex, intent(in), dimension(:,:) :: redfield
    integer a,b,c,i,j,k,ikr,jkr,kkr,info
    integer, allocatable :: ipiv(:,:)
    double precision, allocatable :: rwork(:)
    doublecomplex kspace(dimsq,nkrylov),overlapc,&
            hessenberg(nkrylov,nkrylov),rhoone(dimsq)
    doublecomplex, allocatable :: work(:),wmat(:,:),winv(:,:), lambda(:),coeffs(:),&
            temp1(:),temp2(:,:),vl(:,:)

    !!! Frobenius norm of the density matrix
    overlapc = (0.d0, 0.d0)
    do a = 1,dimsq
        overlapc = overlapc + dconjg(denvec(a)) * denvec(a)
    end do
    denvec = denvec/zsqrt(overlapc)
    !!! First column of the krylov vectors
    kspace = (0.d0, 0.d0)
    do a = 1,dimsq
        kspace(a,1) = denvec(a)
    end do

    !!! BEGIN CONSTRUCTION OF KRYLOV SPACE
    hessenberg = (0.d0, 0.d0)

    do jkr = 1,nkrylov
        !building rhoone
        rhoone = (0.d0,0.d0)

        !$omp parallel do default(none), &
        !$omp private(a,b), &
        !$omp shared(dimsq,redfield,kspace,jkr), &
        !$omp reduction(+:rhoone)
        do a = 1,dimsq
            do b = 1,dimsq
                rhoone(a) = rhoone(a) + redfield(a,b) * kspace(b,jkr)
            end do
        end do
        !$omp end parallel do

        !writing the updiagonal elements of the heseenberg matrix
        !$omp parallel do default(none), &
        !$omp private(ikr,a), &
        !$omp shared(jkr,dimsq,kspace,rhoone), &
        !$omp reduction(+:hessenberg)
        do ikr = 1,jkr
            do a = 1,dimsq
                hessenberg(ikr,jkr) = hessenberg(ikr,jkr) + dconjg(kspace(a,ikr)) * rhoone(a)
            end do
        end do
        !$omp end parallel do

        !building rhotwo
        !$omp parallel do default(none), &
        !$omp private(ikr,a), &
        !$omp shared(jkr,dimsq,kspace,hessenberg), &
        !$omp reduction(+:rhoone)
        do ikr = 1,jkr
            do a = 1,dimsq
                rhoone(a) = rhoone(a) - hessenberg(ikr,jkr) * kspace(a,ikr)
            end do
        end do
        !$omp end parallel do

        if (jkr==nkrylov) go to 999

        !writing the lower diagonakspace(a,1)l element of the hessenberg matrix
        do a = 1,dimsq
            hessenberg(jkr+1,jkr) = hessenberg(jkr+1,jkr) + dconjg(rhoone(a))*rhoone(a)
        end do
        hessenberg(jkr+1,jkr) = zsqrt(hessenberg(jkr+1,jkr))

        !writing the jkr+1 krylov space
        do a = 1,dimsq
            kspace(a,jkr+1) = rhoone(a) / hessenberg(jkr+1,jkr)
        end do

        999 continue
    end do
    !!! END CONSTRUCTION OF KRYLOV SPACE
!    do ikr = 1,nkrylov
!            write(9999,*) (dreal(hessenberg(ikr,jkr)),jkr=1,nkrylov)
!    end do
!    do ikr = 1,nkrylov
!            write(8888,*) (dimag(hessenberg(ikr,jkr)),jkr=1,nkrylov)
!    end do
    !!! diagonalization of hessenberg matrix
    allocate(lambda(nkrylov),wmat(nkrylov,nkrylov), winv(nkrylov,nkrylov),vl(nkrylov,nkrylov))
    allocate(work(2*nkrylov),rwork(2*nkrylov))
    call zgeev('N','V',nkrylov,hessenberg,nkrylov, lambda,vl,1,wmat,nkrylov,work,&
            2*nkrylov,rwork,info)
    write(*,*) 'zgeev info ', info
    winv = wmat !saving wmat. winv, after zgetrf and zgetri, it will contain the inverse of wmat
    deallocate(work,rwork,vl)
    !!! end diagonalization of hessenberg matrix

    !!! matrix inversion of W (here vr)
    allocate(ipiv(nkrylov,nkrylov))
    call zgetrf(nkrylov,nkrylov,winv,nkrylov,ipiv,info)
    write(*,*) 'zgetrf info', info
    allocate(work(nkrylov))
    call zgetri(nkrylov,winv,nkrylov,ipiv,work,nkrylov,info)
    write(*,*) 'zgetri info', info
    deallocate(ipiv,work)
    !!! end matrix inversion

    !!! getting <rho_i|rho_0>
    allocate(temp1(nkrylov))
    temp1 = (0.d0, 0.d0)
    !$omp parallel do default(none), &
    !$omp private(ikr,a), &
    !$omp shared(nkrylov,dimsq,kspace), &
    !$omp reduction(+:temp1)
    do ikr = 1,nkrylov
        do a = 1,dimsq
            temp1(ikr) = temp1(ikr) + dconjg(kspace(a,ikr))*kspace(a,1)
        end do
    end do
    !$omp end parallel do
    !!!

    !!! getting W*exp(lambda*dt)*W-1
    allocate(temp2(nkrylov,nkrylov))
    temp2 = (0.d0, 0.d0)
    !$omp parallel do default(none), &
    !$omp private(ikr,jkr,kkr), &
    !$omp shared(nkrylov,wmat,lambda,dt,winv), &
    !$omp reduction(+:temp2)
    do ikr = 1,nkrylov
        do jkr = 1,nkrylov
            do kkr = 1,nkrylov
                temp2(ikr,jkr) = temp2(ikr,jkr) + wmat(ikr,kkr)*zexp(lambda(kkr)*dt)&
                                                 *winv(kkr,jkr)
            end do
        end do
    end do
    !$omp end parallel do
    deallocate(wmat,lambda,winv)
    !!!

    !!! getting the c_i(t) coefficients
    allocate(coeffs(nkrylov))
    coeffs = (0.d0, 0.d0)
    !$omp parallel do default(none), &
    !$omp private(ikr,jkr), &
    !$omp shared(nkrylov,temp2,temp1), &
    !$omp reduction(+:coeffs)
    do ikr = 1,nkrylov
        do jkr = 1,nkrylov
            coeffs(ikr) = coeffs(ikr) + temp2(ikr,jkr)*temp1(jkr)
        end do
    end do
    !$omp end parallel do
    coeffs = zsqrt(overlapc) * coeffs
    deallocate(temp1,temp2)
    !!!

    !!! writing the density matrix at time t+dt
    denvec = (0.d0, 0.d0)
    !$omp parallel do default(none), &
    !$omp shared(dimsq,nkrylov,coeffs,kspace), &
    !$omp private(a,ikr), &
    !$omp reduction(+:denvec)
    do a = 1,dimsq
        do ikr = 1,nkrylov
            denvec(a) = denvec(a) + coeffs(ikr)*kspace(a,ikr)
        end do
    end do
    !$omp end parallel do
    deallocate(coeffs)
    !!!
end subroutine arnoldi

end program red2
