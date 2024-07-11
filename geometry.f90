program geometry
  implicit none
  real*8:: l, alpha, d, r, theta, c, K
  integer:: i, nsiti
  real*8::pi
  real*8,allocatable::cord(:,:), dist(:,:)
  l=4.d0 !bond lenght
  c=2.d0 !helix pitch
  nsiti=4 !number of sites
  pi=dacos(-1.d0)
  alpha=pi/5 !torsional angle
  r=dsqrt((l**2-c**2)/(2*(1-dcos(alpha))))
  write(*,*)alpha
  open(1,file='geometria.dat')
  open(2,file='geom.dat')
  allocate(cord(nsiti,3), dist(nsiti,nsiti))
  do i=1,nsiti
     if(i.eq.5) write(*,*) 'yy'
     cord(i,1)=r*dcos(pi/2-(i-1)*alpha)
     cord(i,2)=r*dcos((i-1)*alpha)
     cord(i,3)=c*(i-1)
     write(1,*) i, cord(i,1), cord(i,2), cord(i,3)
     write(2,*) cord(i,1), cord(i,2), cord(i,3)
  enddo
  close(1)
end program geometry
