program  flash
use module

  implicit none
  character*1:: jobz,uplo
  integer:: lda, lwork,isito,irwork,liwork,info, lrwork, si, sj,nso,conta, count,perm,phase,bbbb, dimred
  integer,allocatable::iwork(:)
  integer::nsiti,dim2, i, n, a, b, c, d, p,j, k, l, asito, bsito, csito, dsito, kappa
  integer::  temp, m, contasito,sitoi, sitoj, phase2(4,4)
  integer,allocatable:: vecconfig(:), occupazioni(:), NZ(:), spar(:,:)
  real*8,allocatable:: rwork(:), w(:), dist(:,:,:), hop(:,:), nuclei(:,:), hop2(:,:), spintot(:), carica(:,:), u(:),esite(:), mu(:,:,:), charges(:,:), dipole(:,:)
  real*8,allocatable:: muax(:),muay(:), muaz(:),  mubx(:), muby(:), mubz(:), numconf(:,:), numrot(:,:)
  real*8, allocatable:: muralpha(:,:,:), murbeta(:,:,:), singlet(:), triplet(:),quintet(:),w2(:), sr(:,:), tr(:,:), qr(:,:), eig(:,:), spindensity(:,:), sdr(:,:)
  complex*16,allocatable::ham(:,:),work(:), hamsoc(:,:),  soc_a(:,:,:), soc_b(:,:,:),soc_mono(:,:,:), pp(:,:),coup(:,:), COUPLING(:,:),pp2(:,:,:,:),  pp2r(:,:,:,:), now(:,:), now2(:,:), ppso(:,:,:), sx(:,:),sy(:,:), sz(:,:), ssqx(:,:),ssqy(:,:), ssqz(:,:),srot(:,:,:), hopax(:,:), hopay(:,:), hopaz(:,:), hopbx(:,:), hopby(:,:), hopbz(:,:)
  complex*16,allocatable::ppax(:,:), ppbx(:,:),ppay(:,:), ppaz(:,:), ppby(:,:), ppbz(:,:), pprax(:,:), pprbx(:,:),ppray(:,:), pprby(:,:), ppraz(:,:), pprbz(:,:), hssotb(:,:,:,:,:),ssotb(:,:,:,:),sso(:,:),mom(:,:,:), ham2(:,:), sqrot(:,:),sqrot2(:,:), hsootb(:,:,:,:,:), sootb(:,:,:,:), soo(:,:), soc(:,:), socr(:,:), mono_coupx(:,:),mono_coupy(:,:), mono_coupz(:,:), exchange(:,:,:,:), exc(:,:), numh(:), numa(:), numl(:), numhr(:,:), numar(:,:), numlr(:,:)
  complex*16,allocatable:: bi_coupx(:,:),bi_coupy(:,:), bi_coupz(:,:), temporary(:,:,:,:,:), mcrotx(:,:), bcrotx(:,:),mcroty(:,:), mcrotz(:,:), bcroty(:,:), bcrotz(:,:), soor(:,:), ssor(:,:), variable(:,:,:,:), soc_monox(:,:), soc_monoy(:,:), soc_monoz(:,:),psi0(:),mux(:,:), muy(:,:), muz(:,:),  muarx(:,:), muary(:,:), muarz(:,:), mubrx(:,:), mubry(:,:), mubrz(:,:), mutrrot(:,:) , mutrtb(:,:), mutrconf(:,:)
  real*8,allocatable:: dsite(:,:), ssite(:), spin3(:), spin2(:), pol(:,:), polr(:,:), tt(:,:), pot(:), energy(:),  dipolex(:), dipoley(:), dipolez(:), essequadro(:)
  real*8:: Uc, t1, t2, PPP, me, gs, e, e0, pi, cl, radius, j12, homo, lumo, bridge
  logical:: bool, bool1, bool2, bool3, is_hermitian
  real*8::strenght,hbar,hbarev,echarge,emass,dpg, var1, var2, delta, gamma, gamma2, tollerance, bibi,norm, mud
  complex*16::spin(2,2,3),vec1(3),vec2(3), cp(3), cp1(3)
  complex*16::cplx,pf,uff
  character*1,allocatable::state(:)
  open(1,file='input_system.dat')
  read(1,*) nsiti
  read(1,*) Uc
  read(1,*) t1
  read(1,*) t2
  read(1,*) j12
  read(1,*) homo
  read(1,*) lumo
  read(1,*) bridge
  read(1,*) delta
  read(1,*) kappa
  read(1,*) dimred
  close(1)
!!$  open(1,file='input.dat')
!!$  read(1,*) kappa
!!$  close(1)
  nso=nsiti*2+2
  me=9.1093837015d-31
  gs=2.00231930436256
  e=1.602176634d-19
  e0=8.8541878128d-12
  pi=dacos(-1.d0)
  mud=1000000
!  write(*,*) pi
  cl=299792458 
  cplx=cmplx(0.d0, 1.d0)
  uff=cmplx(1.d0, 0.d0)  

  pf=((gs*e**2)/(8*pi*e0*me*cl**2))*10d10
  write(*,*) dimred

  open(1,file='basis.dat')
  open(2,file='geom.dat')
  open(3,file='hamiltonian.dat')
  open(4,file='results.dat')
  open(6,file='check.dat')
  open(10,file='dim.dat')
  open(11,file='work.dat')
  read(10,*) dim2
  close(10)
  write(*,*) dim2
  write(*,*) nso

  allocate(vecconfig(dim2), nuclei(nsiti,3), occupazioni(nsiti), soc_a(3,nso,nso), soc_b(3,nso, nso), hop(nsiti+1,nsiti+1),hop2(nso,nso),soc_mono(3,nso,nso),hamsoc(nso,nso), spintot(dim2), carica(dim2,nsiti), nz(nsiti), u(nsiti), esite(nsiti+1), coup(dim2,dim2), COUPLING(DIM2, dim2), spar(dim2,dim2), spin3(dim2), spin2(dim2), state(dim2), dipole(dim2,3), muax(DIM2), muay(dim2), muaz(dim2),ppso(3,nso,nso),tt(dim2,dim2), pot(dim2), energy(dim2), sqrot(dim2,dim2),sqrot2(dim2,dim2),soc(dim2,dim2),socr(dim2,dim2), dipolex(dim2), dipoley(dim2), dipolez(dim2), mubx(dim2), muby(dim2), mubz(dim2), ham2(dim2,dim2),mutrtb(nso,nso), mutrconf(dim2,dim2))
!=========================LETTURA INPUT=========================
  do i=1,dim2
     read(1,*) vecconfig(i)
  enddo
  spin=0
  spin(1,2,1)=1d0
  spin(2,1,1)=1d0

  spin(1,2,2)=cplx
  spin(2,1,2)=-cplx

  spin(1,1,3)=1d0
  spin(2,2,3)=-1d0

  u(1)=uc
  do i=2,nsiti
     u(i)=uc
  enddo
  esite=bridge
  esite(1)=homo
  esite(2)=lumo
  esite(nsiti+1)=delta

  nz=0
  nz(1)=2
  nz(2)=1
  nz(3)=1
   
  do i=1,nsiti
     read(2,*) nuclei(i,1), nuclei(i,2), nuclei(i,3)
  enddo

  do i=2,nsiti+1
     if(i.eq.2)then
        write(6,*) i-1, nz(i-1), u(i-1), esite(i-1), esite(i)
     else
        write(6,*) i-1, nz(i-1), u(i-1), esite(i)
     endif
  enddo
  !=========================EXCHANGE PART=========================
  allocate(ham(dim2,dim2), exchange(nso,nso, nso,nso))
  exchange=0
  ham=0
  !This matrix is written in this order aa, ab, ba, bb
  exchange(1,3,1,3)=1d0
  exchange(1,4,1,4)=-1d0
  exchange(2,3,2,3)=-1.d0
  exchange(1,4,2,3)=2.d0
  exchange(2,3,1,4)=2.d0
  exchange(2,4,2,4)=1d0
  exchange=(j12/4)*exchange
  write(6,*) 'EXCHANGE'

  allocate(exc(dim2,dim2))
 ! exc=0
  call bielectron(dim2, nso,uff , vecconfig, exchange, exc)
  do n=1,dim2
     do m=1,dim2
        if(exc(n,m).ne.0)write(77,*) n, m
     enddo
  enddo
  !=========================CHARGES AND MOMENTUM=========================
  carica=0
  do n=1,dim2
     count=0
     do i=0,3
        if(btest(vecconfig(n),i))count=count+1
     enddo
     carica(n,1)=nz(1)-count
     count=0
     do i=4,nso-2,2
        if(btest(vecconfig(n),i))count=count+1
        if(btest(vecconfig(n),i+1))count=count+1
        carica(n,i/2)=nz(i/2)-count
        count=0
     enddo
  enddo
  call dipole_moment(dipole, carica, nuclei, dim2, nsiti)
 
  dipolex=dipole(:,1)
  dipoley=dipole(:,2)
  dipolez=dipole(:,3)

  muax=0
  muay=0
  muaz=0
  do n=1,dim2
     do i=0,nso-2,2
        if(i.le.3)then
           isito=1
        else
           isito=i/2
        endif
        bool=btest(vecconfig(n), i)
        if(bool)muax(n)=muax(n)-nuclei(isito,1)
        if(bool)muay(n)=muay(n)-nuclei(isito,2)
        if(bool)muaz(n)=muaz(n)-nuclei(isito,3)     
     enddo
  enddo
  
  mubx=0
  muby=0
  mubz=0
  do n=1,dim2
     do i=1,nso-1,2
        if(i.le.3)then
           isito=1
        else
           isito=i/2
        endif
        bool=btest(vecconfig(n), i)
        if(bool)mubx(n)=mubx(n)-nuclei(isito,1)
        if(bool)muby(n)=muby(n)-nuclei(isito,2)
        if(bool)mubz(n)=mubz(n)-nuclei(isito,3)
     enddo
  enddo

  mutrtb=0
  mutrtb(1,3)=1d0
  mutrtb(2,4)=1d0
  mutrtb(3,1)=1d0
  mutrtb(4,2)=1d0
  mutrtb=mutrtb*mud
  mutrconf=0
  
  call  sq_oe_op_compl(nso,dim2,mutrtb,mutrconf,vecconfig)
  !=========================NUM OPERATOR=========================
  allocate(numconf(dim2,nso))
  allocate(numa(dim2), numh(dim2),numl(dim2), numar(dim2,dim2), numhr(dim2,dim2), numlr(dim2,dim2))
  numconf=0
  do n=1,dim2
     do i=0,nso-1
        if(btest(vecconfig(n),i))numconf(n,i+1)=numconf(n,i+1)+1
     enddo
  enddo
  numa=0
  numl=0
  numh=0
  numar=0
  numhr=0
  numlr=0

  do n=1,dim2
     numa(n)=(numconf(n,nso-1)+numconf(n,nso))
     numh(n)=(numconf(n,1)+numconf(n,2))
     numl(n)=(numconf(n,3)+numconf(n,4))
  enddo
 
  !=========================Spin operator=========================

  allocate(sx(nso,nso), sy(nso,nso), sz(nso,nso))
  sx=0
  sy=0
  sz=0

  do i=1,nso-1,2
     sx(i,i)=spin(1,1,1)
     sx(i+1,i+1)=spin(2,2,1)
     sx(i,i+1)=spin(1,2,1)
     sx(i+1,i)=spin(2,1,1)
  enddo

  do i=1,nso-1,2
     sy(i,i)=spin(1,1,2)
     sy(i+1,i+1)=spin(2,2,2)
     sy(i,i+1)=spin(1,2,2)
     sy(i+1,i)=spin(2,1,2)
  enddo

  do i=1,nso-1,2
     sz(i,i)=spin(1,1,3)
     sz(i+1,i+1)=spin(2,2,3)
     sz(i,i+1)=spin(1,2,3)
     sz(i+1,i)=spin(2,1,3)
  enddo

  ssqx=0
  ssqy=0
  ssqz=0
  allocate (ssqx(dim2,dim2), ssqy(dim2,dim2), ssqz(dim2,dim2))
  call sq_oe_op_compl(nso,dim2,sx,ssqx,vecconfig)
  call check_hermitian(ssqx, dim2, is_hermitian)
  if(.not.is_hermitian)write(*,*) 'Problem ssqx'

  call sq_oe_op_compl(nso,dim2,sy,ssqy,vecconfig)
  call check_hermitian(ssqy, dim2, is_hermitian)
  if(.not.is_hermitian)write(*,*) 'Problem ssqy'

  call sq_oe_op_compl(nso,dim2,sz,ssqz,vecconfig)
  call check_hermitian(ssqz, dim2, is_hermitian)
  if(.not.is_hermitian)write(*,*) 'Problem ssqz'
  
  ssqx=0.5*ssqx
  call square_complex_matrix(dim2, ssqx)
  ssqy=0.5*ssqy
  call square_complex_matrix(dim2, ssqy)
  ssqz=0.5*ssqz
 call square_complex_matrix(dim2, ssqz)

  call check_hermitian(ssqx, dim2, is_hermitian)
  if(.not.is_hermitian)write(*,*) 'Not hermitian ssqx^2'
  call check_hermitian(ssqy, dim2, is_hermitian)
  if(.not.is_hermitian)write(*,*) 'Not hermitian ssqy^2'
  call check_hermitian(ssqz, dim2, is_hermitian)
  if(.not.is_hermitian)write(*,*) 'Not hermitian ssqz^2'

  sqrot=0
  do n=1,dim2
     do m=1,dim2
        sqrot(n,m)=sqrot(n,m)+ssqx(n,m)+ssqy(n,m)+ssqz(n,m)
     enddo
  enddo
  call check_hermitian(sqrot, dim2, is_hermitian)
  if(.not.is_hermitian)write(*,*) 'Sqrot not hermitian'
  
  hop=0
  hop(1,3)=t1
  hop(2,3)=t1
  hop(3,4)=t2
  hop(4,5)=t1
  call copy_upper_to_lower(hop, nsiti+1)

  ppso=0
  do k=1,3
     do i=1,nso
        if(i.le.3)then
           sitoi=1
        else
           sitoi=(i-1)/2
        endif
        do j=1,nso
           if(j.le.3)then
              sitoj=1
           else
              sitoj=(j-1)/2
           endif
           if((mod(i,2).eq.0).and.(mod(j,2).eq.0))ppso(k,i,j)=cplx * hop((i+1)/2,(j+1)/2) * (nuclei(sitoi,k)-nuclei(sitoj,k))
           if((mod(i,2).ne.0).and.(mod(j,2).ne.0))ppso(k,i,j)=cplx * hop((i+1)/2,(j+1)/2) * (nuclei(sitoi,k)-nuclei(sitoj,k))
        enddo
     enddo
  enddo
  hop2=0
  do i=1,nso
     if(i.le.3)then
        sitoi=1
     else
        sitoi=(i-1)/2
     endif
     do j=1,nso
        if(j.le.3)then
           sitoj=1
        else
           sitoj=(j-1)/2
        endif
        if((mod(i,2).eq.0).and.(mod(j,2).eq.0)) hop2(i,j)=hop((i+1)/2,(j+1)/2)
        if((mod(i,2).ne.0).and.(mod(j,2).ne.0)) hop2(i,j)=hop((i+1)/2,(j+1)/2) 

     enddo
  enddo

  write(6,*)
  write(6,*) 'HOPPING ON NSITI'
  do i=1,nsiti+1
     write(6,'(<nsiti+1>(2x,f10.5))') (hop(i,j), j=1,nsiti+1)
  enddo

  write(6,*)
  write(6,*) 'HOPPING ON SPIN ORBITAL'
  do i=1,nso
     write(6,'(<nso>(2x,f10.5))') (hop2(i,j), j=1,nso)
  enddo
  
  do k=1,3
     write(6,*) 'PPSO K=',k
     do i=1,nso
        write(6,'(<nso>(2x,f10.5))') (dimag(ppso(k,i,j)), j=1,nso)
     enddo
     call check_hermitian(mom(:,:,k), nsiti, bool)
     if(.not.bool)write(*,*) 'PROBLEMI MOM K=',k
  enddo

  allocate(hopax(nso,nso), hopbx(nso,nso), hopay(nso,nso), hopaz(nso,nso),  hopby(nso,nso), hopbz(nso,nso))
  hopax=0
  do i=1,nso-3,2
     hopax(i,i+2)=ppso(1,i,i+2)
  enddo
  call copia_comp_conj(hopax, nso)
  
  hopay=0
  do i=1,nso-3,2
     hopay(i,i+2)=ppso(2,i,i+2)
  enddo
  call copia_comp_conj(hopay, nso)

  hopaz=0
  do i=1,nso-3,2
     hopaz(i,i+2)=ppso(3,i,i+2)
  enddo
  call copia_comp_conj(hopaz, nso)

  hopbx=0
  do i=2,nso-2,2
     hopbx(i,i+2)=ppso(1,i,i+2)
  enddo
  call copia_comp_conj(hopbx, nso)

  hopby=0
  do i=2,nso-2,2
     hopby(i,i+2)=ppso(2,i,i+2)
  enddo
  call copia_comp_conj(hopby, nso)

  hopbz=0
  do i=2,nso-2,2
     hopbz(i,i+2)=ppso(3,i,i+2)
  enddo
  call copia_comp_conj(hopbz, nso)

  allocate(ppax(dim2,dim2), ppay(dim2,dim2), ppaz(dim2,dim2), ppbx(dim2,dim2), ppby(dim2,dim2), ppbz(dim2,dim2)) 
  call sq_oe_op_compl(nso, dim2,hopax,ppax,vecconfig)
  call check_hermitian(ppax, dim2, bool)
  if(.not.bool)write(*,*) 'ERROR PPAX'

  call sq_oe_op_compl(nso, dim2,hopay,ppay,vecconfig)
  call check_hermitian(ppay, dim2, bool)
  if(.not.bool)write(*,*) 'ERROR PPAY'

  call sq_oe_op_compl(nso, dim2,hopaz,ppaz,vecconfig)
  call check_hermitian(ppaz, dim2, bool)
  if(.not.bool)write(*,*) 'ERROR PPAZ'
  
  call sq_oe_op_compl(nso, dim2,hopbx,ppbx,vecconfig)
  call check_hermitian(ppbx, dim2, bool)
  if(.not.bool)write(*,*) 'ERROR PPBX'

  call sq_oe_op_compl(nso, dim2,hopby,ppby,vecconfig)
  call check_hermitian(ppby, dim2, bool)
  if(.not.bool)write(*,*) 'ERROR PPBY'

  call sq_oe_op_compl(nso, dim2,hopbz,ppbz,vecconfig)
  call check_hermitian(ppbz, dim2, bool)
  if(.not.bool)write(*,*) 'ERROR PPBZ'

  !=========================SPIN DENSITY=========================
  allocate(spindensity(dim2,nsiti),sdr(dim2,nsiti))
  spindensity=0
  do n=1,dim2
     do i=0,3,2
        bool=btest(vecconfig(n),i)
        if(bool)spindensity(n,1)=spindensity(n,1)+1
        bool1=btest(vecconfig(n),i+1)
        if(bool1)spindensity(n,1)=spindensity(n,1)-1
     enddo
     do i=4,nso-2,2
        isito=i/2
        bool=btest(vecconfig(n),i)
        if(bool)spindensity(n,isito)=spindensity(n,isito)+1
        bool1=btest(vecconfig(n),i+1)
        if(bool1)spindensity(n,isito)=spindensity(n,isito)-1
     enddo
  enddo
  do n=1,dim2
     write(77,'(I2,4F15.8,A2)') n, ( spindensity(n,i), i=1,nsiti)
  enddo
!=========================SCRITTURA HAM=========================
  energy = 0.0d0

  !$omp parallel do default(none) private(n, i, isito, bool, bool1) shared(nso, dim2, vecconfig, esite, u, energy)
  do n = 1, dim2
     do i = 0, nso - 1
        isito=(i+2)/2
        if (btest(vecconfig(n), i)) then
           energy(n) = energy(n) + esite(isito)
        endif
     enddo

     do i = 0, nso - 2, 2
        if(i.le.3)then
           isito=1
        else
           isito=i/2
        endif
        bool = btest(vecconfig(n), i)
        bool1 = btest(vecconfig(n), i + 1)
        if (bool .and. bool1) then
           energy(n) = energy(n) + u(isito)
        endif
     enddo
  enddo
  !$omp end parallel do
 
  !=========================SOC=========================
  allocate(mom(nsiti,nsiti,3))
  mom=0
  do i = 1,nsiti+1
     if(i.le.2)then
        sitoi=1
     else
        sitoi=i-1
     endif
     do j = 1,nsiti+1
        if(j.le.2)then
           sitoj=1
        else
           sitoj=j-1
        endif
        do k = 1,3
           var1 = nuclei(sitoi,k)-nuclei(sitoj,k)
           mom(sitoi,sitoj,k) = cplx * var1 * hop(i,j)
        end do
     end do
  end do

  do k=1,3
     write(6,*) 'MOM K=',k
     do i=1,nsiti
        write(6,'(<nsiti>(2x,f10.5))') (dimag(mom(i,j,k)), j=1,nsiti)
     enddo
     call check_hermitian(mom(:,:,k), nsiti, bool)
     if(.not.bool)write(*,*) 'PROBLEMI MOM K=',k
  enddo
 
  soc_a=0
  SOC_B=0
  soc_mono=0
  do i=1,nso
     if(i.le.4)then
        sitoi=1
     else
        sitoi=(i-1)/2
     endif
     do j=1,nso
        if(j.le.4)then
           sitoj=1
        else
           sitoj=(j-1)/2
        endif
        do isito = 1,nsiti
           if (isito.ne.sitoi) then 
              ! r x p cross product
              vec1 = (0.d0, 0.d0)
              vec2 = (0.d0, 0.d0)
              cp=0
              do k = 1,3
                 vec1(k) = nuclei(sitoi,k) - nuclei(isito,k) !position
                 vec2(k) = ppso(k,(i+1)/2,(j+1)/2) !momentum
              end do
              cp  = cross_product(vec1, vec2)

              ! 1/|r_aA|^3 term
              radius = 0.d0
              do k = 1,3
                 radius = radius + dreal(vec1(k))**2 
             end do
              radius = (dsqrt(radius))**3
              cp = cp / radius
              if(i/2*2.eq.i)then
                 si = 2
              else
                 si=1
              endif

              if(j/2*2.eq.j)then
                 sj = 2
              else
                 sj=1
              endif

              do k = 1,3
                 soc_a(k,i,j) = soc_a(k,i,j) + cp(k)  * spin(si,sj,k)
                 !if(soc_a(k,i,j).ne.0)write(77,*) soc_a(k,i,j)
              end do
           end if
        end do

     end do
  end do
  
  do i=1,nso
     if(i.le.3)then
        sitoi=1
     else
        sitoi=(i-1)/2
     endif
     do j=1,nso
        if(j.le.3)then
           sitoj=1
        else
           sitoj=(j-1)/2
        endif
        do isito = 1,nsiti
           if(isito.ne.sitoj)then 
              ! r x p cross product
              vec1 = (0.d0, 0.d0)
              vec2 = (0.d0, 0.d0)
              cp=0
              do k = 1,3
                 vec1(k) = nuclei(sitoj,k) - nuclei(isito,k) !position
                 vec2(k) = ppso(k,(i+1)/2,(j+1)/2) !momentum
                 ! write(*,*) vec2             
              end do
              cp  = cross_product(vec1, vec2)
              !write(*,*) cp

              ! 1/|r_aA|^3 term
              radius = 0.d0
              do k = 1,3
                 radius = radius + dreal(vec1(k))**2
              end do
              radius = (dsqrt(radius))**3

              cp = cp / radius
              
              if(i/2*2.eq.i)then
                 si = 2
              else
                 si= 1
              endif

              if(j/2*2.eq.j)then
                 sj = 2
              else
                 sj= 1
              endif

              do k = 1,3
                 soc_b(k,i,j) = soc_b(k,i,j) + cp(k)  * spin(si,sj,k)
                 !if(soc_b(k,i,j).ne.0)write(77,*) soc_a(k,i,j)
              end do
           end if
        end do

     end do
  end do

  do i = 1,nso
     do j = 1,nso
        do k = 1,3
           soc_mono(k,i,j)= 0.5d0 * (soc_a(k,i,j) + soc_b(k,i,j))
        end do
     end do
  end do
  
  hamsoc=0d0
  do i=1,nso
     do j=1,nso
        do k=1,3
           hamsoc(i,j)=hamsoc(i,j)+soc_mono(k,i,j)
        enddo
     enddo
  enddo
  write(6,*) 'SOC MONO'
  do n=1,nso
     do m=1,nso
        if(hamsoc(n,m).ne.0)write(6,*) hamsoc(n,m)
     enddo
  enddo
        
  call  check_hermitian(soc_mono, nso, is_hermitian)
  if(is_hermitian) write(*,*) 'soc_mono HERMITIANA'
  call  check_hermitian(hamsoc, nso, is_hermitian)
  if(is_hermitian) write(*,*) 'hamsoc_HERMITIANA'
  if(.not.is_hermitian) write(*,*) 'PROBLEMI'

  coup=0
  call sq_oe_op_compl(nso,dim2,hamsoc,coup,vecconfig)
  call  check_hermitian(coup, dim2, is_hermitian)
  if(is_hermitian) write(*,*) 'COUP HERMITIANA'
  
  coup=pf*coup
  !=========================TWO TERM SS0=========================
  !SSO TERM
  allocate(dist(nsiti+1,nsiti+1,k), hssotb(3,nso,nso,nso,nso), ssotb(nso,nso,nso,nso), sso(dim2,dim2))
  dist=0
  do i=1,nsiti+1
     do j=1,nsiti+1
        do k=1,3
           if(i.eq.1)then
              isito=1
           else
              isito=i-1
           endif
           if(j.eq.1)then
              sitoj=1
           else
              sitoj=j-1
           endif
           dist(i,j,k)=nuclei(isito,k)-nuclei(sitoj,k)
        enddo
     enddo
  enddo

  write(6,*) 'DISTANCES'
  do k=1,3
     write(6,*) 'K=',k
     do i=1,5
        write(6,'(<5>(2x,f10.5))') (dist(i,j,k), j=1,5)
     enddo
  enddo
  
  hssotb=0
 
  do a=1,nso
     if(a.le.4)then
        asito=1
     else
        asito=(a-1)/2
     endif
     do b=1,nso
        if(b.le.4)then
           bsito=1
        else
           bsito=(b-1)/2
        endif
        do c=1,nso
           if(c.le.4)then
              csito=1
           else
              csito=(c-1)/2
           endif
           do d=1,nso
              if(d.le.4)then
                 dsito=1
              else
                 dsito=(d-1)/2
              endif
              if((b.eq.d).and.(asito.ne.bsito))then
                 vec1=0
                 vec2=0
                 cp=0
                 do k=1,3
                    vec1(k)=dist((a+1)/2,(b+1)/2,k)
                    vec2(k)=ppso(k,a,c)
                 enddo
                 cp=cross_product(vec1,vec2)
                  !write(*,*) cp
                 ! 1/|r_aA|^3 term
                 radius = 0.d0
                 do k = 1,3
                    radius = radius + dreal(vec1(k))**2
                 end do
                 radius = (dsqrt(radius))**3
                 ! write(*,*) radius
                 cp = cp / radius
                 if(mod(a,2).eq.0)then
                    si=2
                 else
                    si=1
                 endif

                 if(mod(c,2).eq.0)then
                    sj=2
                 else
                    sj=1
                 endif

                 do k=1,3
                    hssotb(k,a,b,c,d)=hssotb(k,a,b,c,d)+cp(k)*spin(si,sj,k)
                 enddo

              endif
           enddo
        enddo
     enddo
  enddo

  ssotb=0
  allocate(temporary(3,nso,nso,nso,nso), variable(nso,nso,nso,nso))
  temporary=0

  do k=1,3
     do a=1,nso
        do b=1,nso
           do c=1,nso
              do d=1,nso
                 ssotb(a,b,c,d)=ssotb(a,b,c,d)+0.5d0 * (hssotb(k,a,b,c,d)+ dconjg(hssotb(k,c,d,a,b)))
              enddo
           enddo
        enddo
     enddo
  enddo

  !Inizio a passare dal tb al real space
  sso=0
  call bielectron(dim2,nso,pf,vecconfig,ssotb,sso)
  call check_hermitian(sso, dim2, bool)
  if(.not.bool)write(*,*)'sso problem'

  !SOO TERM
  allocate(hsootb(3,nso,nso,nso,nso), sootb(nso,nso,nso,nso), soo(dim2,dim2))
  hsootb=0

  do a=1,nso
     if(a.le.3)then
        asito=1
     else
        asito=(a-1)/2
     endif
     do b=1,nso
        if(b.le.3)then
           bsito=1
        else
           bsito=(b-1)/2
        endif
        do c=1,nso
           if(c.le.3)then
              csito=1
           else
              csito=(c-1)/2
           endif
           do d=1,nso
              if(d.le.3)then
                 dsito=1
              else
                 dsito=(d-1)/2
              endif
              if((bsito.eq.dsito).and.(asito.ne.bsito))then
                 vec1=0
                 vec2=0
                 cp=0
                 do k=1,3
                    vec1(k)=dist((a+1)/2,(b+1)/2,k)
                    vec2(k)=ppso(k,a,c)
                 enddo
                 cp=cross_product(vec1,vec2)
                 ! write(*,*) cp
                 ! 1/|r_aA|^3 term
                 radius = 0.d0
                 do k = 1,3
                    radius = radius + dreal(vec1(k))**2
                 end do
                 radius = (dsqrt(radius))**3
                 ! write(*,*) radius
                 cp = cp / radius
                 if(mod(b,2).eq.0)then
                    si=2
                 else
                    si=1
                 endif

                 if(mod(d,2).eq.0)then
                    sj=2
                 else
                    sj=1
                 endif

                 do k=1,3
                    hsootb(k,a,b,c,d)=hsootb(k,a,b,c,d)+cp(k)*spin(si,sj,k)
                    ! write(*,*) hssotb(k,a,b,c,d)
                 enddo
              endif
           enddo
        enddo
     enddo
  enddo

  sootb=0

  do k=1,3
     do a=1,nso
        do b=1,nso
           do c=1,nso
              do d=1,nso
                 sootb(a,b,c,d)=sootb(a,b,c,d)+0.5d0*(hsootb(k,a,b,c,d)+dconjg(hsootb(k,c,d,a,b)))
                 temporary(k,a,b,c,d)=temporary(k,a,b,c,d)+0.5d0*(hsootb(k,a,b,c,d)+dconjg(hsootb(k,c,d,a,b)))
              enddo
           enddo
        enddo
     enddo
  enddo

  do a=1,nso
     do b=1,nso
        do c=1,nso
           do d=1,nso
              if(zabs(sootb(a,b,c,d)-dconjg(sootb(c,d,a,b))).ge.1d-10) write(*,*) a, b, c, d, zabs(sootb(a,b,c,d)-dconjg(sootb(c,d,a,b)))
           enddo
        enddo
     enddo
  enddo

  !Inizio a passare dal tb al real space
  soo=0

  call  bielectron(dim2,nso,pf,vecconfig,sootb,soo)
  call check_hermitian(soo,dim2,bool)
  if(.not.bool)write(*,*)'soo problem'

  allocate(bi_coupx(dim2,dim2), bi_coupy(dim2,dim2), bi_coupz(dim2,dim2))
  bi_coupx=0
  bi_coupy=0
  bi_coupz=0
  
  variable=temporary(1,:,:,:,:)
  call  bielectron(dim2,nso,pf,vecconfig,variable,bi_coupx)

  variable=temporary(2,:,:,:,:)
  call  bielectron(dim2,nso,pf,vecconfig,variable,bi_coupy)

  variable=temporary(3,:,:,:,:)
  call  bielectron(dim2,nso,pf,vecconfig,variable,bi_coupz)
   
  !============================================================
 
  !! Adding to hamiltonian
  ham2=0
  call  sq_oe_op_real(nso,dim2,hop2,tt,vecconfig)
  do n=1,dim2
     ham(n,n)=ham(n,n)+energy(n)
     ham2(n,n)=ham(n,n)+ham2(n,n)
     do m=1,dim2
        ham2(n,m)=ham2(n,m)+ham(n,m)-tt(n,m)+exc(n,m)
        ham(n,m)=ham(n,m) - tt(n,m) + exc(n,m)+ (coup(n,m)-soo(n,m)-sso(n,m))*kappa
     enddo
  enddo
  !HAM BREIT PAULI HAMILTONIAN
  
  soc=0
  do n=1,dim2
     do m=1,dim2
        soc(n,m)=soc(n,m)-sso(n,m)-soo(n,m)+coup(n,m)
     enddo
  enddo
  
  call  check_hermitian(ham, dim2, is_hermitian)
  if(is_hermitian) write(*,*) 'HAM HERMITIANA'
  jobz ='V'
  uplo='U'
  lrwork=(1+5*dim2+2*dim2**2)
  liwork=(3+5*dim2)
  lwork=(2*dim2+dim2**2)
  allocate(w(dim2),work(max(1,lwork)),rwork(lrwork), iwork(max(1,liwork)),w2(dim2), eig(dim2,3))
 
!==============================DIAGONALIZATION=========================
  call zheevd (jobz, uplo, dim2, ham, dim2, w, work, lwork,rwork,lrwork,iwork,liwork,info)
 
  call eigenvalues(dim2,1d-8,w,state)
  sqrot2=0
  
  call rotate_cplx_2x2(dim2, sqrot2, sqrot, ham)
 
  call check_hermitian(sqrot2, dim2, bool)
  if(.not.bool)write(*,*) 'AAAAH'
!!!RUOTO S^2 SULLA BASE DEGLI AUTOSTATI DELL'HAMILTONIANO DI BREIT-PAULI  
  write(4,*) 'EIGENVALUES'
  do i=1,dim2
    ! write(4, '(I0.3,2x, F15.7,2x, A,2x, F15.7)') i, w(i)-w(1), state(i), dreal(sqrot2(i,i))
     write(4,*) i, w(i)-w(1), state(i), dreal(sqrot2(i,i))
  enddo
  
  !==========================SPIN DENSITY ROTATION=========================  
  sdr=0
  call rot_diag(dim2, sdr, spindensity, nsiti, ham)
  write(*,*) 'CI SONO'
  write(4,*) 'SPIN DENSITY'
  do n=1,dim2
     write(4,'(I4,4F15.8,A2)') n, sdr(n,1), sdr(n,2), sdr(n,3), sdr(n,4), state(n)
  enddo
 
!=========================OUTPUT=========================
  allocate(charges(dim2,nsiti))
  charges=0
 
  call rot_diag(dim2,charges,carica,nsiti,ham)
  
  write(4,*) 'CARICHE'
  do n=1,dim2
     write(4,'(<5>(4x,f10.5))') w(n)-w(1), ( charges(n,i), i=1,nsiti)
  enddo
  
  allocate(mux(dim2,dim2), muy(dim2,dim2), muz(dim2,dim2), pp2r(dim2,dim2,2,3),muarx(dim2,dim2),muary(dim2,dim2), muarz(dim2,dim2), mubrx(dim2,dim2), mubry(dim2,dim2), mubrz(dim2,dim2), mutrrot(dim2,dim2))

  call rotate_rtc_1x2(dim2, muarx, muax, ham)
  call rotate_rtc_1x2(dim2, muary, muay, ham)
  call rotate_rtc_1x2(dim2, muarz, muaz, ham)
  call rotate_rtc_1x2(dim2, mubrx, mubx, ham)
  call rotate_rtc_1x2(dim2, mubry, muby, ham)
  call rotate_rtc_1x2(dim2, mubrz, mubz, ham)
  
  write(4,*) 'DIPOLE ALPHA'
  write(4,'(<3>(2x,f10.5))') muarx(1,1), muary(1,1), muarz(1,1)
  write(4,*) 'DIPOLE BETA'
  write(4,'(<3>(2x,f10.5))') mubrx(1,1), mubry(1,1), mubrz(1,1)
!!!RUOTO IL MOMENTO DI DIPOLO DEGLI SPIN ALFA E DEGLI SPIN BETA SULLA BASE DEGLI AUTOSTATI BREIT-PAULI
  
!!$  write(4,*) 'DIPOLE ALPHA+BETA'
!!$  write(4,'(<3>(2x,f10.5))') (muralpha(1,1,k)+murbeta(1,1,k), k=1,3)
  call rotate_rtc_1x2(dim2, mux, dipolex, ham)
  call rotate_rtc_1x2(dim2, muy, dipoley, ham)
  call rotate_rtc_1x2(dim2, muz, dipolez, ham)
  call rotate_cplx_2x2(dim2, mutrrot, mutrconf, ham)
  muz=0

  do n=1,dim2
     do m=1,dim2
        muz(n,m)=muz(n,m)+mutrrot(n,m)
       ! muz(n,m)=mutrrot(n,m)
     enddo
  enddo
 
  allocate(pprax(dim2,dim2),pprbx(dim2,dim2),ppray(dim2,dim2), pprby(dim2,dim2), ppraz(dim2,dim2), pprbz(dim2,dim2))
  pprax=0
  pprbx=0
  ppray=0
  pprby=0
  ppraz=0
  pprbz=0
  
  call rotate_cplx_2x2(dim2, pprax, ppax, ham)
  call rotate_cplx_2x2(dim2, ppray, ppay, ham)
  call rotate_cplx_2x2(dim2, ppraz, ppaz, ham)
  call rotate_cplx_2x2(dim2, pprbx, ppbx, ham)
  call rotate_cplx_2x2(dim2, pprby, ppby, ham)
  call rotate_cplx_2x2(dim2, pprbz, ppbz, ham)
  
!!$  write(6,*) "DIPOLE MOMENT"
!!$  write(6,*) 'DIPOLE X'
!!$  call write_matrix (dreal(mux), 6, dim2,dim2, dim2)
!!$  write(6,*) 'DIPOLE Y'
!!$  call write_matrix (dreal(muy), 6, dim2,dim2, dim2)
!!$  write(6,*) 'DIPOLE Z'
!!$  call write_matrix (dreal(muz), 6, dim2,dim2, dim2)
  
  write(4,*) 'LINEAR MOMENT ALPHA'
  write(4,*) dimag(pprax(1,1)),dimag(ppray(1,1)), dimag(ppraz(1,1))
  write(4,*) 'LINEAR MOMENT BETA'
  write(4,*) dimag(pprbx(1,1)),dimag(pprby(1,1)), dimag(pprbz(1,1))
!!! RUOTO E SCRIVO IL MOMENTO DI DIPOLO, E IL MOMENTO LINEARE SULLA BASE DEGLI AUTOSTATI BREIT-PAULI

  !=========================NUM OPERATO=========================
  allocate(numrot(dim2,nso))
  numrot=0
  call rot_diag(dim2, numrot, numconf, nso, ham)
  call rotate_cplx_1x2(dim2,numar,numa,ham)
  call rotate_cplx_1x2(dim2,numhr,numh,ham)
  call rotate_cplx_1x2(dim2,numlr,numl,ham)
  call  check_hermitian(numar, dim2, bool)
  if(.not.bool)write(*,*)'No numar'
  call  check_hermitian(numhr, dim2, bool)
  if(.not.bool)write(*,*)'No numhr'
  call  check_hermitian(numLr, dim2, bool)
  if(.not.bool)write(*,*)'No numLr'  
  write(6,*) 'NUMERO'
  do i=1,dim2
     do j=1,dim2
        if(zabs(numar(i,j)).ge.1d-5)write(6,*)i,j,numar(i,j)
        if(zabs(numhr(i,j)).le.1d-5)write(6,*)i,j,numar(i,j)
        if(zabs(numlr(i,j)).le.1d-5)write(6,*)i,j,numar(i,j)
     write(6,'(I0.3, 2x, <nso+2>(2x,f10.5))'), i,w(i), dreal(sqrot2(i,i)), (numrot(i,j), j=1,nso)
     enddo
  enddo

  !=========================CHECK=========================
!!$  w=0
!!$  call zheevd (jobz, uplo, dim2, ham2, dim2, w, work, lwork,rwork,lrwork,iwork,liwork,info)
!!$
!!$  do i=1,dim2
!!$     eig(i,1)=w(i)-w(1)
!!$  enddo
!!$  w2=w
!!$  coupling=0
!!$  call rotate_cplx_2x2(dim2,coupling,coup,ham2)
!!$  call check_hermitian(coupling, dim2,bool)
!!$  if(.not.bool)write(*,*) 'Problemi coupling'
!!$
!!$  allocate(soor(dim2,dim2), ssor(dim2,dim2))
!!$  soor=0
!!$  ssor=0
!!$  call rotate_cplx_2x2(dim2,soor,soo,ham2)
!!$  call check_hermitian(soor, dim2,bool)
!!$  if(.not.bool)write(*,*) 'Problemi soor'
!!$
!!$  call rotate_cplx_2x2(dim2,ssor,sso,ham2)
!!$  call check_hermitian(ssor, dim2,bool)
!!$  if(.not.bool)write(*,*) 'Problemi ssor'
!!$  
!!$  socr=0
!!$  call rotate_cplx_2x2(dim2,socr,soc,ham2)
!!$  call check_hermitian(socr,dim2,bool)
!!$  if(.not.bool)write(*,*) 'problemi socr'
!!$  
!!$  sqrot2=0
!!$  call  rotate_cplx_2x2(dim2,sqrot2,sqrot,ham2)
!!$  call check_hermitian(sqrot2,dim2,bool)
!!$  if(.not.bool)write(*,*) 'problemi sqrot2'
!!$  
!!$!!! RUOTO IL COUPLING MONO, SOO, SSO E TOTALE SULLA BASE DEGLI AUTOSTATI DI HUBBARD
!!$  singlet=0
!!$  triplet=0
!!$  quintet=0
!!$  allocate(singlet(dim2), triplet(dim2), quintet(dim2))
!!$  do n=1,dim2
!!$     if(dreal(sqrot2(n,n)).le.1d-15)singlet(n)=1d0
!!$     if(dabs(dreal(sqrot2(n,n))-2d0).le.1d-12)triplet(n)=1d0
!!$     if(dabs(dreal(sqrot2(n,n))-6d0).le.1d-12)quintet(n)=1d0
!!$  enddo
!!$
!!$  allocate(sr(dim2,dim2), tr(dim2,dim2), qr(dim2,dim2))
!!$  sr=0
!!$  tr=0
!!$  qr=0
!!$  call rotate_real_1x2(dim2,sr,singlet,ham)
!!$  call rotate_real_1x2(dim2,tr,triplet,ham)
!!$  call rotate_real_1x2(dim2,qr,quintet,ham)
!!$
!!$  open(9876,file='GS.dat')
!!$  open(9877,file='S1.dat')
!!$  open(9878,file='mcrot_r.dat')
!!$  open(9879,file='mcrot_i.dat')
!!$  open(9883,file='bcrot_r.dat')
!!$  open(9882,file='bcrot_i.dat')
!!$  
!!$  
!!$  open(1000,file='composition-real.dat')
!!$  open(1001,file='composition-imag.dat')
!!$  open(1002,file='composition-real-s1.dat')
!!$  open(1003,file='composition-imag-s1.dat')
!!$  
!!$  do i=1,dim2
!!$     write(1000,*) dreal(coupling(1,i)), -dreal(ssor(1,i)), -dreal(soor(1,i))
!!$     write(1001,*) dimag(coupling(1,i)), -dimag(ssor(1,i)), -dimag(soor(1,i))
!!$     write(1002,*) dreal(coupling(5,i)), -dreal(ssor(5,i)), -dreal(soor(5,i))
!!$     write(1003,*) dimag(coupling(5,i)), -dimag(ssor(5,i)), -dimag(soor(5,i))
!!$  enddo

 !=========================INITIAL STATE PREPARATION=========================
 allocate(psi0(dim2))
 psi0=0
 do i=2,dim2
    psi0(i)=muz(i,1)
 enddo
 norm=0.d0
 do i=1,dim2
    norm=norm+dconjg(psi0(i))*psi0(i)
 enddo
 psi0=psi0/dsqrt(norm)

 do i=1, dim2
    write(6,*)i, (zabs(psi0(i)))**2
 enddo
 eig=0
 do i=1,dim2
    eig(i,1)=w(i)-w(1)
   ! write(*,*) eig(i,1)
 enddo

!=========================REDFIELD-INPUTS=========================
 open(55,file='input-red/mux.bin',form="unformatted")
 write(55)mux(:dimred,:dimred)
 open(56,file='input-red/muy.bin',form="unformatted")
 write(56)muy(:dimred,:dimred)
 open(57,file='input-red/muz.bin',form="unformatted")
 write(57)muz(:dimred,:dimred)
 open(66,file='input-red/eigen.bin',form="unformatted")
 write(66) eig(:dimred,1)
 open(77,file='input-red/spin-density.bin',form="unformatted")
 write(77)sdr(:dimred,:nsiti)
 open(88,file='input-red/psi0.bin',form="unformatted")
 write(88) psi0(:dimred)
 open(99,file='input-red/numh.bin',form="unformatted")
 write(99) numhr(:dimred,:dimred)
 close(99)
 open(99,file='input-red/numl.bin',form="unformatted")
 write(99) numlr(:dimred,:dimred)
 close(99)
 open(99,file='input-red/numa.bin',form="unformatted")
 write(99) numar(:dimred,:dimred)
 close(99)
 open(99,file='input-red/num.bin',form="unformatted")
 write(99) numrot(:dimred,:nso)
 close(99)
 close(55)
 close(56)
 close(57)
 close(66)
 close(77)
 close(88)

 open(55,file='input-red/system_input.dat')
 write(55,*) dimred
 write(55,*) nsiti
end program flash



