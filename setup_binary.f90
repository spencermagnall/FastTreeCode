module setup_binary
  use sphere_dist
 implicit none
 contains 
 subroutine init(x,v,m,np)
  integer, intent(in) :: np
  real, intent(out) :: x(3,np)
  real, intent(out) :: v(3,np)
  real, intent(out) :: m(np)
  real :: e, a 
  real :: mtotal, q, m1,m2
  real :: center1(3), center2(3)
  integer :: start, end, i 

  e = 0.0
  a = 1000.0
  x = 0.0
  v = 0.0
  m = 0.0
  !m(1) = 1.0
  !m(2) = 1.0
  mtotal = m1 + m2
  q = m2/m1
  center1 = 0.0
  center2 = 0.0

  mtotal = np
  q = 1.0 
  m(:) = 1.0

  start = 1
  end = np/2
  center1(1) = (1-e)/(1+q)*a
  call add_particles(x,np,center1,start,end)

  center2(1) = -q * center1(1)
  start = end + 1
  end = np
  call add_particles(x,np,center2,start,end)

  start = 1
  end = np/2
  do i=start, end
     v(2,i) = 1/(1+q) * sqrt(mtotal/a * ((1.0+e)/(1.-e)))
     write(*,*) "Velocity: ", v(2,i)
  enddo

  start = end + 1 
  end = np/2  
  v(2,start:np) = -q*v(2,1:end)
     

 end subroutine init

 subroutine add_particles(x,np,center,start,end)
  integer, intent(in) :: np,start,end 
  real, intent(out) :: x(3,np)
  real, intent(in) :: center(3)
  real :: halfpart
  integer i

  ! radius of the sphere 
  real :: radius

  radius = 1
  

  ! GENERATES PARTICLES IN SPHERICAL DIST AROUND CENTER 

  call setup_particles(x,np,radius,center,start,end)





 end subroutine add_particles


end module setup_binary
