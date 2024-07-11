program hamiltonian
use module

  implicit none
  character*1:: jobz,uplo
  integer:: lda, lwork,isito,irwork,liwork,info, lrwork, si, sj,nso,conta, count,perm,phase,bbbb, dimred,dimm2, dimm,dimz,dimp,dimp2
  integer,allocatable::iwork(:), oldpositions(:)
  integer::nsiti,dim2, i, n, a, b, c, d, p,j, k, l, gst, s1, t1,tm1
  integer::   m, contasito,sitoi, sitoj
  integer,allocatable:: vecconfig(:), occupazioni(:), NZ(:), spar(:,:), basism2(:), basism(:), basisz(:), basisp(:), basisp2(:)
  real*8,allocatable:: rwork(:), w(:), dist(:,:,:), hop(:,:), nuclei(:,:), hop2(:,:), spintot(:), carica(:,:), u(:),esite(:), mu(:,:,:), charges(:,:), dipole(:,:)
  real*8,allocatable:: muax(:),muay(:), muaz(:), dist2(:,:), mubx(:), muby(:), mubz(:)
  real*8, allocatable:: muralpha(:,:,:), murbeta(:,:,:), singlet(:), triplet(:),quintet(:),w2(:), sr(:,:), tr(:,:), qr(:,:), eig(:,:), spindensity(:,:), sdr(:,:)
  complex*16,allocatable::ham(:,:),work(:), hamsoc(:,:),  soc_a(:,:,:), soc_b(:,:,:),soc_mono(:,:,:), pp(:,:),coup(:,:), COUPLING(:,:),pp2(:,:,:,:),  pp2r(:,:,:,:), now(:,:), now2(:,:), ppso(:,:,:), sx(:,:),sy(:,:), sz(:,:), ssqx(:,:),ssqy(:,:), ssqz(:,:),srot(:,:,:), hopax(:,:), hopay(:,:), hopaz(:,:), hopbx(:,:), hopby(:,:), hopbz(:,:), hamm2(:,:), hamm(:,:), hamz(:,:), hamp(:,:), hamp2(:,:), rotation(:,:)
  complex*16,allocatable::ppax(:,:), ppbx(:,:),ppay(:,:), ppaz(:,:), ppby(:,:), ppbz(:,:), pprax(:,:), pprbx(:,:),ppray(:,:), pprby(:,:), ppraz(:,:), pprbz(:,:), hssotb(:,:,:,:,:),ssotb(:,:,:,:),sso(:,:),mom(:,:,:), ham2(:,:), sqrot(:,:),sqrot2(:,:), hsootb(:,:,:,:,:), sootb(:,:,:,:), soo(:,:), soc(:,:), socr(:,:), mono_coupx(:,:),mono_coupy(:,:), mono_coupz(:,:)
  complex*16,allocatable:: bi_coupx(:,:),bi_coupy(:,:), bi_coupz(:,:), temporary(:,:,:,:,:), mcrotx(:,:), bcrotx(:,:),mcroty(:,:), mcrotz(:,:), bcroty(:,:), bcrotz(:,:), soor(:,:), ssor(:,:), tempo(:,:,:,:), soc_monox(:,:), soc_monoy(:,:), soc_monoz(:,:),psi0(:),mux(:,:), muy(:,:), muz(:,:),  muarx(:,:), muary(:,:), muarz(:,:), mubrx(:,:), mubry(:,:), mubrz(:,:)
  real*8,allocatable:: dsite(:,:), ssite(:), spin3(:), spin2(:), pol(:,:), polr(:,:), tt(:,:), pot(:), energy(:),  dipolex(:), dipoley(:), dipolez(:), energym2(:), energym(:), energyz(:), energyp(:), energyp2(:), ttm2(:,:), ttm(:,:), ttz(:,:), ttp(:,:), ttp2(:,:), eigm2(:), eigm(:), eigz(:), eigp(:), eigp2(:), eigenvalue(:)
  real*8:: Uc, t, PPP, me, gs, e, e0, pi, cl, radius
  logical:: bool, bool1, bool2, bool3, is_hermitian
  real*8::strenght,hbar,hbarev,echarge,emass,dpg, var1, var2, delta, gamma, gamma2, tollerance, bibi,norm,temp
  complex*16::spin(2,2,3),vec1(3),vec2(3), cp(3), cp1(3)
  complex*16::cplx,pf
  character*2,allocatable::state(:)
  character*1,allocatable::molt(:)
  !input reading
  open(1,file='system-input.dat')
  read(1,*) nsiti
  read(1,*) dimred
  read(1,*) uc
  read(1,*) t
  read(1,*) delta
  close(1)
  tollerance=1d-8
  nso=nsiti*2
  me=9.1093837015d-31
  gs=2.00231930436256
  e=1.602176634d-19
  e0=8.8541878128d-12
  pi=dacos(-1.d0)
!  write(*,*) pi
  cl=299792458 
  cplx=cmplx(0.d0, 1.d0)
  pf=((gs*e**2)/(8*pi*e0*me*cl**2))*10d10
  write(*,*) pf
  open(2,file='geom.dat')
  open(3,file='hamiltonian.dat')
  open(4,file='results.dat')
  open(6,file='check.dat')
  open(10,file='dim2.dat')
  open(11,file='work.dat')
  read(10,*) dim2
  read(10,*) dimm2
  read(10,*) dimm
  read(10,*) dimz
  read(10,*) dimp
  read(10,*) dimp2
  close(10)

  allocate(vecconfig(dim2), nuclei(nsiti,3),ham(dim2,dim2), occupazioni(nsiti), soc_a(3,nso,nso), soc_b(3,nso, nso), hop(nsiti,nsiti),hop2(nso,nso), soc_mono(3,nso,nso), hamsoc(nso,nso), spintot(dim2), carica(dim2,nsiti),nz(nsiti), u(nsiti), esite(nsiti), coup(dim2,dim2), coupling(DIM2, dim2), spar(dim2,dim2), spin3(dim2), spin2(dim2), state(dim2), dipole(dim2,3), muax(DIM2), muay(dim2), muaz(dim2),ppso(3,nso,nso),tt(dim2,dim2), pot(dim2), energy(dim2), sqrot(dim2,dim2),sqrot2(dim2,dim2),soc(dim2,dim2), socr(dim2,dim2), dipolex(dim2), dipoley(dim2), dipolez(dim2), mubx(dim2), muby(dim2), mubz(dim2))
!=========================LETTURA INPUT=========================
  allocate( basism2(dimm2), basism(dimm), basisz(dimz), basisp(dimp), basisp2(dimp2))
  basism2=0
  basism=0
  basisz=0
  basisp=0
  basisp2=0
  vecconfig=0
  open(1,file='basis-2.dat')
  do i=1,dimm2
     read(1,*) basism2(i)
  enddo
  close(1)
  open(1,file='basism.dat')
  do i=1,dimm
     read(1,*) basism(i)
  enddo
  close(1)
  open(1,file='basisz.dat')
  do i=1,dimz
     read(1,*) basisz(i)
  enddo
  close(1)
  open(1,file='basisp.dat')
  do i=1,dimp
     read(1,*) basisp(i)
  enddo
  close(1)
  open(1,file='basis2.dat')
  do i=1,dimp2
     read(1,*) basisp2(i)
  enddo
  close(1)
  do i=1,dim2
     if(i.eq.dimm2) vecconfig(i)=basism2(i)
     if((i.ge.(dimm2+1).and.(i.le.dimm2+dimm))) vecconfig(i) = basism(i-dimm2)
     if((i.le.(dimm2+dimm+dimz)).and.(i.ge.(dimm2+dimm+1))) vecconfig(i) = basisz(i-dimm2-dimm)
     if((i.le.(dimm2+dimm+dimz+dimp)).and.(i.ge.(dimm2+dimm+dimz+1))) vecconfig(i) = basisp(i-dimm2-dimm-dimz) 
     if((i.le.(dimm2+dimm+dimz+dimp+dimp2)). and .(i.ge.(dimm2+dimm+dimz+dimp+1))) vecconfig(i) = basisp2(i-dimm2-dimm-dimz-dimp)
  enddo
  
  spin=0
  spin(1,2,1)=1d0
  spin(2,1,1)=1d0

  spin(1,2,2)=cplx
  spin(2,1,2)=-cplx

  spin(1,1,3)=1d0
  spin(2,2,3)=-1d0
  do i=1,nsiti
     u(i)=uc
  enddo
  esite=0
 ! esite(2)=0!2d0
  esite(1)=-delta
  esite(nsiti)=+delta

  do i=nsiti-2,nsiti-1
     nz(i)=1
     esite(i)=0
  enddo
  nz(nsiti)=0
  nz(1)=2
 
  do i=1,nsiti
     read(2,*) nuclei(i,1), nuclei(i,2), nuclei(i,3)
  enddo
  do i=1,nsiti
     write(6,*) i, nz(i), u(i), esite(i)
  enddo
  allocate(dist2(nsiti, 3))
  call calcola_distanze(nuclei, dist2, nsiti)
  write(6,*) 'DISTANCES'
  call write_matrix(dist2, 6, nsiti,nsiti, nsiti)
  !=======================CHARGES AND MOMENTUM=======================
  call charge(carica, vecconfig, nz, dim2, nso)
  call dipole_moment(dipole, carica, nuclei, dim2, nsiti)
  dipolex=dipole(:,1)
  dipoley=dipole(:,2)
  dipolez=dipole(:,3)
  muax=0
  muay=0
  muaz=0
  do n=1,dim2
     do i=0,nso-2,2
        isito=(i+2)/2
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
        isito=(i+1)/2
        bool=btest(vecconfig(n), i)
        if(bool)mubx(n)=mubx(n)-nuclei(isito,1)
        if(bool)muby(n)=muby(n)-nuclei(isito,2)
        if(bool)mubz(n)=mubz(n)-nuclei(isito,3)
     enddo
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

  write(6,*) 'REAL X'
  call write_matrix(dreal(sx), 6, nso, nso, nso)
  write(6,*) 'IMAG X'
  call write_matrix(dimag(sx), 6, nso, nso, nso)

  write(6,*) 'REAL Y'
  call write_matrix(dreal(sy), 6, nso, nso, nso)
  write(6,*) 'IMAG Y'
  call write_matrix(dimag(sy), 6, nso, nso, nso)

  write(6,*) 'REAL Z'
  call write_matrix(dreal(sz), 6, nso, nso, nso)
  write(6,*) 'IMAG Z'
  call write_matrix(dimag(sz), 6, nso, nso, nso)

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
  
  hop2=0
  do i=1,nso-2
     hop2(i,i+2)=t
  enddo
  call copy_upper_to_lower(hop2, nso)
  call  sq_oe_op_real(nso,dim2,hop2,tt,vecconfig) 
  write(6,*) 'HOP'
  hop=0
  do i=1,nsiti-1
     hop(i,i+1)=-t
  enddo
  call copy_upper_to_lower(hop, nsiti)
  call momentum_so (nsiti, nso, hop, nuclei, ppso)
  call  write_matrix (hop, 6, nsiti, nsiti, nsiti)
 
  do k=1,3
    ! write(6,*) 'PPSO K=',k
    ! call  write_matrix (dimag(ppso(k,:,:)), 6, nso, nso, nso)
     call check_hermitian(ppso(k,:,:), nso, bool)
     if(.not.bool)write(*,*) 'PROBLEMI PPSO K=',k
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
     do i=0,nso-1,2
        isito=(i+2)/2
        bool=btest(vecconfig(n),i)
        if(bool)spindensity(n,isito)=spindensity(n,isito)+1
        bool1=btest(vecconfig(n),i+1)
        if(bool1)spindensity(n,isito)=spindensity(n,isito)-1
     enddo
  enddo
!=========================SCRITTURA HAM=========================
  allocate(energym2(dimm2), energym(dimm),energyz(dimz), energyp(dimp), energyp2(dimp2)) 
  ham=0
 ! call ohno (dim2,nsiti,nuclei,vecconfig, u,nz, pot)
  !call site_energy(nso,dim2,esite,vecconfig,energy)
  call site_energy_u(nso,dim2,esite,u,vecconfig,energy)
  call site_energy_u(nso,dimm2,esite, u,basism2,energym2)
  write(*,*) energym2(1)
  call site_energy_u(nso,dimm,esite, u,basism,energym)
  call site_energy_u(nso,dimz,esite, u,basisz,energyz)
  call site_energy_u(nso,dimp,esite, u,basisp,energyp)
  call site_energy_u(nso,dimp2,esite, u,basisp2,energyp2)
  write(*,*) energyp2(1)

  !off diagonal part
  allocate(ttm2(dimm2,dimm2), ttm(dimm,dimm), ttz(dimz,dimz), ttp(dimp,dimp), ttp2(dimp2,dimp2))
  call  sq_oe_op_real(nso,dimm2,hop2,ttm2,basism2)
  call check_symmetric(ttm2, dimm2, bool)
  if(.not.bool)write(*,*) 'ttm2'
  call  sq_oe_op_real(nso,dimm,hop2,ttm,basism)
  call check_symmetric(ttm, dimm, bool)
  if(.not.bool)write(*,*) 'ttm'
  call  sq_oe_op_real(nso,dimz,hop2,ttz,basisz)
  call check_symmetric(ttz, dimz, bool)
  if(.not.bool)write(*,*) 'ttz'
  call  sq_oe_op_real(nso,dimp,hop2,ttp,basisp)
  call check_symmetric(ttp, dimp, bool)
  if(.not.bool)write(*,*) 'ttp'
  call  sq_oe_op_real(nso,dimp2,hop2,ttp2,basisp2)
  call check_symmetric(ttp2, dimp2, bool)
  if(.not.bool)write(*,*) 'ttp2'
  !=======================SOC=======================
  allocate(mom(nsiti,nsiti,3))
  mom=0
  do i = 1,nsiti
     do j = 1,nsiti
        do k = 1,3
           var1 = nuclei(i,k)-nuclei(j,k)
           mom(i,j,k) = cplx * var1 * hop(i,j)
        end do
     end do
  end do
  do k=1,3
     write(6,*) 'MOM K=',k
     call  write_matrix (dimag(mom(:,:,k)), 6, nsiti, nsiti, nsiti)
     call check_hermitian(mom(:,:,k), nsiti, bool)
     if(.not.bool)write(*,*) 'PROBLEMI MOM K=',k
  enddo
  soc_a=0
  SOC_B=0
  soc_mono=0
  do i=1,nso
     do j=1,nso
        do isito = 1,nsiti
           if (isito.ne.(i+1)/2) then 
              ! r x p cross product
              vec1 = (0.d0, 0.d0)
              vec2 = (0.d0, 0.d0)
              cp=0
              sitoi=(i+1)/2
              sitoj=(j+1)/2
              do k = 1,3
                 vec1(k) = nuclei(sitoi,k) - nuclei(isito,k) !position
                 vec2(k) = mom(sitoi,sitoj,k) !momentum
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
     do j=1,nso
        do isito = 1,nsiti
           if(isito.ne.(j+1)/2)then 
              ! r x p cross product
              vec1 = (0.d0, 0.d0)
              vec2 = (0.d0, 0.d0)
              cp=0
              sitoi=(i+1)/2
              sitoj=(j+1)/2           
              do k = 1,3
                 vec1(k) = nuclei(sitoj,k) - nuclei(isito,k) !position
                 vec2(k)= mom(sitoi,sitoj,k)  !momentum
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

  allocate(soc_monox(nso,nso), soc_monoy(nso,nso), soc_monoz(nso,nso))
  soc_monox=soc_mono(1,:,:)
  soc_monoy=soc_mono(2,:,:)
  soc_monoz=soc_mono(3,:,:)

  allocate(mono_coupx(dim2,dim2), mono_coupy(dim2,dim2), mono_coupz(dim2,dim2))
  mono_coupx=0
  mono_coupz=0
  mono_coupy=0

  call sq_oe_op_compl(nso,dim2,soc_monox,mono_coupx,vecconfig)
  call check_hermitian(mono_coupx,dim2, bool)
  if(.not.bool)write(*,*) 'P1'
  call sq_oe_op_compl(nso,dim2,soc_monoy,mono_coupy,vecconfig)
  call check_hermitian(mono_coupy, dim2, bool)
  if(.not.bool)write(*,*) 'P2'
  call sq_oe_op_compl(nso,dim2,soc_monoz,mono_coupz,vecconfig)
  call check_hermitian(mono_coupz, dim2, bool)
  if(.not.bool)write(*,*) 'P3'
  mono_coupx=pf*mono_coupx
  mono_coupy=pf*mono_coupy
  mono_coupz=pf*mono_coupz
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
  !=======================TWO TERM SS0=======================
  !SSO TERM
  allocate(dist(nsiti,nsiti,k), hssotb(3,nso,nso,nso,nso), ssotb(nso,nso,nso,nso), sso(dim2,dim2))
  dist=0
  do i=1,nsiti
     do j=1,nsiti
        do k=1,3
           dist(i,j,k)=nuclei(i,k)-nuclei(j,k)
        enddo
     enddo
  enddo
  hssotb=0
 
  do a=1,nso
     do b=1,nso
        do c=1,nso
           do d=1,nso
              if((b.eq.d).and.((a+1)/2.ne.(b+1)/2))then
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
                    hssotb(k,a,b,c,d) = hssotb(k,a,b,c,d) + cp(k)*spin(si,sj,k)
                    ! write(*,*) hssotb(k,a,b,c,d)
                 enddo                 
              endif
           enddo
        enddo
     enddo
  enddo

  ssotb=0
  allocate(temporary(3,nso,nso,nso,nso), tempo(nso,nso,nso,nso))
  temporary=0
  
  do k=1,3
     do a=1,nso
        do b=1,nso
           do c=1,nso
              do d=1,nso
                 ssotb(a,b,c,d) = ssotb(a,b,c,d) + 0.5d0*(hssotb(k,a,b,c,d) + dconjg(hssotb(k,c,d,a,b)))
                 temporary(k,a,b,c,d) = temporary(k,a,b,c,d) +0.5d0*(hssotb(k,a,b,c,d) + dconjg(hssotb(k,c,d,a,b)))
              enddo
           enddo
        enddo
     enddo
  enddo

  do a=1,nso
     do b=1,nso
        do c=1,nso
           do d=1,nso
              if(zabs(ssotb(a,b,c,d)-dconjg(ssotb(c,d,a,b))).ge.1d-10) write(*,*) a, b, c, d, zabs(ssotb(a,b,c,d)-dconjg(ssotb(c,d,a,b)))
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
     do b=1,nso
        do c=1,nso
           do d=1,nso
              if(((b+1)/2.eq.(d+1)/2).and.((a+1)/2.ne.(b+1)/2))then
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
                    hsootb(k,a,b,c,d) = hsootb(k,a,b,c,d)+ cp(k)*spin(si,sj,k)
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
                 sootb(a,b,c,d) = sootb(a,b,c,d) + 0.5d0*(hsootb(k,a,b,c,d) + dconjg(hsootb(k,c,d,a,b)))
                 temporary(k,a,b,c,d) = temporary(k,a,b,c,d)+ 0.5d0*(hsootb(k,a,b,c,d)+ dconjg(hsootb(k,c,d,a,b)))
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
              if(zabs(sootb(a,b,c,d)).ge.1d-10) write(77,*) a,b,c,d
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
  
  tempo=temporary(1,:,:,:,:)
  call  bielectron(dim2,nso,pf,vecconfig,tempo,bi_coupx)
  tempo=temporary(2,:,:,:,:)
  call  bielectron(dim2,nso,pf,vecconfig,tempo,bi_coupy)
  tempo=temporary(3,:,:,:,:)
  call  bielectron(dim2,nso,pf,vecconfig,tempo,bi_coupz)
   
  !============================================================
  allocate(ham2(dim2,dim2),hamm2(dimm2,dimm2), hamm(dimm,dimm),hamz(dimz,dimz), hamp(dimp,dimp), hamp2(dimp2,dimp2))
  hamm2=0
  hamm=0
  hamz=0
  hamp=0
  hamp2=0
  !! Adding to hamiltonian
  do n=1,dimm2
     hamm2(n,n)=energym2(n)
     do m=1,dimm2
        hamm2(n,m)=hamm2(n,m)-ttm2(n,m)
     enddo
  enddo

  do n=1,dimm
     hamm(n,n)=energym(n)
     do m=1,dimm
        hamm(n,m)=hamm(n,m)-ttm(n,m)
     enddo
  enddo

  do n=1,dimz
     hamz(n,n)=energyz(n)
     do m=1,dimz
        hamz(n,m)=hamz(n,m)-ttz(n,m)
     enddo
  enddo

  do n=1,dimp
     hamp(n,n)=energyp(n)
     do m=1,dimp
        hamp(n,m)=hamp(n,m)-ttp(n,m)
     enddo
  enddo
  
  do n=1,dimp2
     hamp2(n,n)=energyp2(n)
     do m=1,dimp2
        hamp2(n,m)=hamp2(n,m)-ttp2(n,m)
     enddo
  enddo

  call check_hermitian(hamm2,dimm2,bool)
  if(.not.bool)write(*,*)'hamm2 problem'

  call check_hermitian(hamm,dimm,bool)
  if(.not.bool)write(*,*)'hamm problem'

  call check_hermitian(hamz,dimz,bool)
  if(.not.bool)write(*,*)'hamz problem'

  call check_hermitian(hamp,dimp,bool)
  if(.not.bool)write(*,*)'hamp problem'

  call check_hermitian(hamp2,dimp2,bool)
  if(.not.bool)write(*,*)'hamp2 problem'
  
  !HAM BREIT PAULI HAMILTONIAN
  
  soc=0
  do n=1,dim2
     do m=1,dim2
        soc(n,m)=soc(n,m)-sso(n,m)-soo(n,m)+coup(n,m)
     enddo
  enddo
  call check_hermitian(soc,dim2,bool)
  if(.not.bool)write(*,*)'SOC problem'
  
  allocate(eigm2(dimm2), eigm(dimm), eigz(dimz), eigp(dimp), eigp2(dimp2), eigenvalue(dim2), oldpositions(dim2), molt(dim2))
  oldpositions=0

  jobz ='V'
  uplo='U'
  !============================DIAGONALIZATION=======================
  lrwork=(1+5*dimm2+2*dimm2**2)
  liwork=(3+5*dimm2)
  lwork=(2*dimm2+dimm2**2)
  allocate(w(dimm2),work(max(1,lwork)),rwork(lrwork), iwork(max(1,liwork)))
  call zheevd (jobz, uplo, dimm2, hamm2, dimm2, w, work, lwork,rwork,lrwork,iwork,liwork,info)
  eigm2=w
  deallocate(w, work,rwork,iwork)
  
  lrwork=(1+5*dimm+2*dimm**2)
  liwork=(3+5*dimm)
  lwork=(2*dimm+dimm**2)
  allocate(w(dimm),work(max(1,lwork)),rwork(lrwork), iwork(max(1,liwork)))
  call zheevd (jobz, uplo, dimm, hamm, dimm, w, work, lwork,rwork,lrwork,iwork,liwork,info)
  eigm=w
  deallocate(w, work,rwork,iwork)
  
  lrwork=(1+5*dimz+2*dimz**2)
  liwork=(3+5*dimz)
  lwork=(2*dimz+dimz**2)
  allocate(w(dimz),work(max(1,lwork)),rwork(lrwork), iwork(max(1,liwork)))
  call zheevd (jobz, uplo, dimz, hamz, dimz, w, work, lwork,rwork,lrwork,iwork,liwork,info)
  eigz=w
  deallocate(w, work,rwork,iwork)

  lrwork=(1+5*dimp+2*dimp**2)
  liwork=(3+5*dimp)
  lwork=(2*dimp+dimp**2)
  allocate(w(dimp),work(max(1,lwork)),rwork(lrwork), iwork(max(1,liwork)))
  call zheevd (jobz, uplo, dimp, hamp, dimp, w, work, lwork,rwork,lrwork,iwork,liwork,info)
  eigp=w
  deallocate(w, work,rwork,iwork)

  lrwork=(1+5*dimp2+2*dimp2**2)
  liwork=(3+5*dimp2)
  lwork=(2*dimp2+dimp2**2)
  allocate(w(dimp2),work(max(1,lwork)),rwork(lrwork), iwork(max(1,liwork)))
  call zheevd (jobz, uplo, dimp2, hamp2, dimp2, w, work, lwork,rwork,lrwork,iwork,liwork,info)
  eigp2=w
  deallocate(w, work,rwork,iwork) 
  
  ham=0
  do n=1,dimm2
     state(n)='-2'
     eigenvalue(n)=eigm2(1)
     do m=1,dimm2
        ham(n,m)=ham(n,m)+hamm2(n,m)
     enddo
  enddo

  do n=dimm2+1,dimm2+dimm
     state(n)='-1'
     eigenvalue(n)=eigm(n-dimm2)
     do m=dimm2+1,dimm2+dimm
        ham(n,m)=ham(n,m)+hamm(n-dimm2,m-dimm2)
     enddo
  enddo

  do n=dimm2+dimm+1,dimm2+dimm+dimz
     state(n)='0'
     eigenvalue(n)=eigz(n-dimm2-dimm)
     do m=dimm2+dimm+1,dimm2+dimm+dimz
        ham(n,m)= ham(n,m)+ hamz(n-dimm2-dimm,m-dimm2-dimm)
     enddo
  enddo

  do n=dimm2+dimm+dimz+1,dimm2+dimm+dimz+dimp
     state(n)='+1'
     eigenvalue(n)=eigp(n-dimm2-dimm-dimz)
     do m=dimm2+dimm+dimz+1,dimm2+dimm+dimz+dimp
        ham(n,m)= ham(n,m)+ hamp(n-dimm2-dimm-dimz,m-dimm2-dimm-dimz)
     enddo
  enddo

  do n=dimm2+dimm+dimz+dimp+1, dimm2+dimm+dimz+dimp+dimp2
     eigenvalue(n)=eigp2(n-dimm2-dimm-dimz-dimp)
     state(n)='+2'
     do m=dimm2+dimm+dimz+dimp+1, dimm2+dimm+dimz+dimp+dimp2
        ham(n,m)=ham(n,m)+ hamp2(n-dimm2-dimm-dimz-dimp,m-dimm2-dimm-dimz-dimp)
     enddo
  enddo

  do i=1,dim2
     eigenvalue(i)=eigenvalue(i)-minval(eigz)
  enddo
  allocate(w(dim2))
  w=eigenvalue
  call bubble_sort(eigenvalue, dim2, oldpositions)
  call eigenvalues(dim2, tollerance, eigenvalue, molt)
  gst= oldpositions(1)
  !  s1=oldpositions(k)
  l=find_character_position(state, dim2, '+1')
  m=find_character_position(state, dim2, '-1')
  t1=oldpositions(4)
  tm1=oldpositions(2)

  write(4,*) 'EIGENVALUES'
  do i=1,12
     write(4,'((f10.8,1x),3x,2(a2))') eigenvalue(i), molt(i), state(oldpositions(i)) 
  enddo
  do i=1,dim2
     write(111,'(<dim2>(2x,f10.5))') (dreal(ham(i,j)), j=1,dim2)
  enddo
  
  call rotate_cplx_2x2(dim2,sqrot2,sqrot,ham)
  call check_hermitian(sqrot2, dim2, bool)
  if(.not.bool)write(*,*) 'AAAAH'

!===========================OUTPUT===========================
  allocate(charges(dim2,nsiti))
  charges=0
  call rot_diag(dim2,charges,carica,nsiti,ham)
  write(4,*) 'CARICHE'
  do i=1,10
     write(4,'(4(f15.8,1x),3x,2(a2))') charges(oldpositions(i),1), charges(oldpositions(i),2), charges(oldpositions(i),3), charges(oldpositions(i),4), molt(i), state(oldpositions(i))
  enddo
  !======================SPIN DENSITY ROTATION=====================
  sdr=0
  do n=1,dim2
     do m=1,dim2
        do l=1,nsiti
           sdr(n,l)=sdr(n,l)+spindensity(m,l)*zabs(ham(m,n))**2
        enddo
     enddo
  enddo
  write(4,*) 'SPIN DENSITY'
  do i=1,10
    write(4,'(4(f15.8,1x),3x,2(a2))') sdr(oldpositions(i),1), sdr(oldpositions(i),2), sdr(oldpositions(i),3), sdr(oldpositions(i),4), molt(i), state(oldpositions(i))
 enddo
!=====================DIPOLE MOMENT ROTATION=====================
  allocate(mux(dim2,dim2), muy(dim2,dim2), muz(dim2,dim2), pp2r(dim2,dim2,2,3),muarx(dim2,dim2),muary(dim2,dim2), muarz(dim2,dim2), mubrx(dim2,dim2), mubry(dim2,dim2), mubrz(dim2,dim2))

  call rotate_rtc_1x2(dim2, muarx, muax, ham)
  call rotate_rtc_1x2(dim2, muary, muay, ham)
  call rotate_rtc_1x2(dim2, muarz, muaz, ham)
  call rotate_rtc_1x2(dim2, mubrx, mubx, ham)
  call rotate_rtc_1x2(dim2, mubry, muby, ham)
  call rotate_rtc_1x2(dim2, mubrz, mubz, ham)
  
  write(4,*) 'DIPOLE ALPHA'
  write(4,'(<3>(2x,f10.5))') muarx(gst,gst), muary(gst,gst), muarz(gst,gst)
  write(4,*) 'DIPOLE BETA'
  write(4,'(<3>(2x,f10.5))') mubrx(gst,gst), mubry(gst,gst), mubrz(gst,gst)
!!!RUOTO IL MOMENTO DI DIPOLO DEGLI SPIN ALFA E DEGLI SPIN BETA SULLA BASE DEGLI AUTOSTATI BREIT-PAULI

  call rotate_rtc_1x2(dim2, mux, dipolex, ham)
  call rotate_rtc_1x2(dim2, muy, dipoley, ham)
  call rotate_rtc_1x2(dim2, muz, dipolez, ham)
!=====================LINEAR MOMENT ROTATION=====================
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
  write(4,*) dimag(pprax(1,1)),dimag(ppray(gst,gst)), dimag(ppraz(gst,gst))
  write(4,*) 'LINEAR MOMENT BETA'
  write(4,*) dimag(pprbx(gst,gst)),dimag(pprby(gst,gst)), dimag(pprbz(gst,gst))
  
!!! RUOTO E SCRIVO IL MOMENTO DI DIPOLO, E IL MOMENTO LINEARE SULLA BASE DEGLI AUTOSTATI BREIT-PAULI
!=========================COUPLING ANALYSIS=========================
  coupling=0
  call rotate_cplx_2x2(dim2,coupling,coup,ham)
  call check_hermitian(coupling, dim2,bool)
  if(.not.bool)write(*,*) 'Problemi coupling'

  allocate(soor(dim2,dim2), ssor(dim2,dim2))
  soor=0
  ssor=0
  call rotate_cplx_2x2(dim2,soor,soo,ham)
  call check_hermitian(soor, dim2,bool)
  if(.not.bool)write(*,*) 'Problemi soor'

  call rotate_cplx_2x2(dim2,ssor,sso,ham)
  call check_hermitian(ssor, dim2,bool)
  if(.not.bool)write(*,*) 'Problemi ssor'
  
  socr=0
  call rotate_cplx_2x2(dim2,socr,soc,ham)
  call check_hermitian(socr,dim2,bool)
  if(.not.bool)write(*,*) 'problemi socr'
  
  sqrot2=0
  call  rotate_cplx_2x2(dim2,sqrot2,sqrot,ham)
  call check_hermitian(sqrot2,dim2,bool)
  if(.not.bool)write(*,*) 'problemi sqrot2'

  allocate(mcrotx(dim2,dim2),mcroty(dim2,dim2), mcrotz(dim2,dim2),bcrotx(dim2,dim2),bcroty(dim2,dim2), bcrotz(dim2,dim2))
  mcrotx=0
  mcroty=0
  mcrotz=0
  
  call  rotate_cplx_2x2(dim2,mcrotx,mono_coupx,ham)
  call  rotate_cplx_2x2(dim2,mcrotz,mono_coupz,ham)
  call  rotate_cplx_2x2(dim2,mcroty,mono_coupy,ham)
  call  rotate_cplx_2x2(dim2,bcrotx,bi_coupx,ham)
  call  rotate_cplx_2x2(dim2,bcroty,bi_coupy,ham)
  call  rotate_cplx_2x2(dim2,bcrotz,bi_coupz,ham)
  
!!! RUOTO IL COUPLING MONO, SOO, SSO E TOTALE SULLA BASE DEGLI AUTOSTATI DI HUBBARD
  singlet=0
  triplet=0
  quintet=0
  allocate(singlet(dim2), triplet(dim2), quintet(dim2))
  do m=1,dim2
     if(dreal(sqrot2(n,n)).le.1d-15)singlet(n)=1d0
     if(dabs(dreal(sqrot2(n,n))-2d0).le.1d-12)triplet(n)=1d0
     if(dabs(dreal(sqrot2(n,n))-6d0).le.1d-12)quintet(n)=1d0
  enddo

  allocate(sr(dim2,dim2), tr(dim2,dim2), qr(dim2,dim2))
  sr=0
  tr=0
  qr=0
  call rotate_real_1x2(dim2,sr,singlet,ham)
  call rotate_real_1x2(dim2,tr,triplet,ham)
  call rotate_real_1x2(dim2,qr,quintet,ham)
 
  open(9876,file='GS.dat')
  open(9877,file='T1.dat')
  open(9875,file='T-1.dat')
  open(9878,file='mcrot_r.dat')
  open(9879,file='mcrot_i.dat')
  open(9883,file='bcrot_r.dat')
  open(9882,file='bcrot_i.dat')
   
  do i=1,dim2
     m=find_real(eigenvalue, w(i), dim2)
     temp=(minval(w))
!!$     write(9876,'(2(f24.18,2x),3x,2(a2),3x,f15.8)' )dreal(socr(gst,i)), dimag(socr(gst,i)), molt(m), state(i),w(i)-temp
!!$     write(9875,'(2(f24.18,2x),3x,2(a2),3x,f15.8)')dreal(socr(tm1,i)), dimag(socr(tm1,i)),  molt(m), state(i), w(i)-temp
!!$     write(9877,'(2(f24.18,2x),3x,2(a2),3x,f15.8)') dreal(socr(t1,i)), dimag(socr(t1,i)),  molt(m), state(i), w(i)-temp
     write(9876,'(2(f24.18,2x),3x,2(a2),3x,f15.8)' )dreal(coupling(gst,i)), dimag(coupling(gst,i)), molt(m), state(i),w(i)-temp
     write(9875,'(2(f24.18,2x),3x,2(a2),3x,f15.8)')dreal(socr(tm1,i)), dimag(socr(tm1,i)),  molt(m), state(i), w(i)-temp
     write(9877,'(2(f24.18,2x),3x,2(a2),3x,f15.8)') dreal(socr(t1,i)), dimag(socr(t1,i)),  molt(m), state(i), w(i)-temp  
  enddo
!!$
!!$  do i=1,dim2
!!$     write(9878,*)  dreal(mcrotx(gst,i)), dreal(mcroty(gst,i)), dreal(mcrotx(gst,i))
!!$     write(9879,*)   dimag(mcrotx(gst,i)), dimag(mcroty(gst,i)), dimag(mcrotz(gst,i))
!!$     write(9883,*)- dreal(bcrotx(gst,i)), - dreal(bcroty(gst,i)), - dreal(bcrotz(gst,i))
!!$     write(9882,*)- dimag(bcrotx(gst,i)), - dimag(bcroty(gst,i)), - dimag(bcrotz(gst,i))
!!$  enddo
!!$
!!$  open(1111,file='mcrot_r_s1.dat')
!!$  open(1112,file='mcrot_i_s1.dat')
!!$  open(1113,file='bcrot_r_s1.dat')
!!$  open(1114,file='bcrot_i_s1.dat')
!!$
!!$  do i=1,dim2
!!$     write(1111,*)  dreal(mcrotx(s1,i)), dreal(mcroty(s1,i)), dreal(mcrotx(s1,i))
!!$     write(1112,*)   dimag(mcrotx(s1,i)), dimag(mcroty(s1,i)), dimag(mcrotz(s1,i))
!!$     write(1113,*)- dreal(bcrotx(s1,i)), - dreal(bcroty(s1,i)), - dreal(bcrotz(s1,i))
!!$     write(1114,*)- dimag(bcrotx(s1,i)), - dimag(bcroty(s1,i)), - dimag(bcrotz(s1,i))
!!$  enddo
!!$
!!$  open(1000,file='composition-real.dat')
!!$  open(1001,file='composition-imag.dat')
!!$  open(1002,file='composition-real-s1.dat')
!!$  open(1003,file='composition-imag-s1.dat')
!!$  
!!$  do i=1,dim2
!!$     write(1000,*) dreal(coupling(gst,i)), -dreal(ssor(gst,i)), -dreal(soor(gst,i))
!!$     write(1001,*) dimag(coupling(gst,i)), -dimag(ssor(gst,i)), -dimag(soor(gst,i))
!!$     write(1002,*) dreal(coupling(s1,i)), -dreal(ssor(s1,i)), -dreal(soor(s1,i))
!!$     write(1003,*) dimag(coupling(s1,i)), -dimag(ssor(s1,i)), -dimag(soor(s1,i))
!!$  enddo
end program hamiltonian
