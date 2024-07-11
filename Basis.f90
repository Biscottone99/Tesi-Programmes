program basis
  implicit none
  integer:: nsiti,i,max,min, check1, check2, nso, count, config, nf, a,b, n, spin, ne, dimm, dimp, dimz, dimm2, dim2
 ! real*8::spin
  character:: array(8)
  logical::bool
  nsiti=4 !siti di disposizione elettroni
  nso=nsiti*2
  ne=4
  open(1,file='basis.dat')
  open(2,file='configurations.dat')
  open(3,file='dim2.dat')
  open(4,file='basism.dat')
  open(5,file='basisp.dat')
  open(6,file='basisz.dat')
  open(7,file='basis-2.dat')
  open(8,file='basis2.dat')
  max=0
  do i=nso-1, nso-ne,-1
     max=max+2**i
  enddo
  min=0
  do i=0,ne-1
     min=min+2**i
  enddo
!!!! Starting with writing of configurations
  nf=0
  dimm=0
  dimp=0
  dimz=0
  dimm2=0
  dim2=0
  do n=min,max
     count=0
     config=0
     a=0
     b=0
     spin=0
     do i=0,nso-1
        bool=btest(n,i)
        if(bool)then
           array(i+1)='1'
           count=count+1
           if(i/2*2.eq.i)then
              a=a+1
           else
              b=b+1
           endif
        else
           array(i+1)='0'          
        endif
     enddo
     spin=(a-b)*0.5d0
     if((count.eq.ne))then
        config=0
        do i=0,nso-1
           if(array(i+1).eq.'1')then
              config=config+2**i
           endif
        enddo
        write(1,*) config
        write(2,*) (array(i),i=nso,1,-1), spin
        nf=nf+1
        if(spin.eq.-1)then
           dimm=dimm+1
           write(4,*) config
        endif
        if(spin.eq.0)then
           dimz=dimz+1
           WRITE(6,*) config
        endif
        if(spin.eq.1)then
           dimp=dimp+1
           write(5,*) config
        endif
        if(spin.eq.-2)then
           dimm2=dimm2+1
           write(7,*) config
        endif
        if(spin.eq.2)then
           dim2=dim2+1
           write(8,*) config
        endif
     endif
  enddo
  write(*,*)'numero di funzioni di base uguale', nf, dimm, dimp, dimz, dim2, dimm2
  write(3,*) nf
  write(3,*) dimm2
  write(3,*) dimm
  write(3,*) dimz
  write(3,*) dimp
  write(3,*) dim2
endprogram basis
