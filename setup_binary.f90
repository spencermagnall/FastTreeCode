module setup_binary
  use sphere_dist
 implicit none
 contains
 subroutine init_pointmass(x,v,m,np)
  integer, intent(in) :: np
  real, intent(out) :: x(3,np)
  real, intent(out) :: v(3,np)
  real, intent(out) :: m(np)
  real :: e, a 
  real :: mtotal, q, m1,m2
  real :: center1(3), center2(3)
  integer :: start, end, i 


  e = 0.7
  a = 1000.0
  x = 0.0
  v = 0.0
  m = 0.0
  !m(1) = 1.0
  !m(2) = 1.0
 

  mtotal = np
  q = 1.0 
  m(:) = 1.0
  
  x(1,1) = (1-e)/(1+q)*a
  x(1,2) = -q*x(1,1)



  v(2,1) = 1/(1+q) * sqrt(mtotal/a * ((1.0+e)/(1.-e)))

  v(2,2) = -q*v(2,1)

end subroutine init_pointmass

 subroutine init(x,v,m,np)
  integer, intent(in) :: np
  real, intent(out) :: x(3,np)
  real, intent(out) :: v(3,np)
  real, intent(out) :: m(np)
  real :: e, a 
  real :: mtotal, q, m1,m2
  real :: center1(3), center2(3), r(3),r2,rmag,rhat(3),vhat(3)
  integer :: start, end, i 

  e = 0.5
  a = 10000.0
  x = 0.0
  v = 0.0
  m = 0.0
  !m(1) = 1.0
  !m(2) = 1.0
  mtotal = m1 + m2
  q = m2/m1
  center1 = 0.0
  center2 = 0.0
  r = 0
  r2 = 0
  rhat = 0
  rmag =0
  vhat =0

  mtotal = np
  q = 1.0 
  m(:) = 1.0

  start = 1
  end = np/2
  write(*,*) "END IS = ", end 
  center1(1) = (1-e)/(1+q)*a
  !center1(2) = center1(1)
  !center1(3) = center1(1)
  call add_particles(x,np,center1,start,end)

  ! Get unit vector in direction of r
  r = center1
  r2 = dot_product(r,r)
  rmag = sqrt(r2)
  rhat = r/rmag
  ! get  vhat
  !vhat(1) = rhat(1)
  !vhat(2) = -rhat(2)
  !vhat(3) = rhat(3)

  write(*,*) "Dot product: ", dot_product(vhat,rhat)


  center2(1) = -q * center1(1)
  start = end + 1
  write(*,*) "START IS = ", start
  end = np
  call add_particles(x,np,center2,start,end)

  start = 1
  end = np/2
  do i=start, end
     !v(:,i) = vhat
     v(2,i) = 1/(1+q) * sqrt(mtotal/a * ((1.0+e)/(1.-e)))
     !write(*,*) "Velocity: ", v(2,i)
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

  radius = 100.0
  

  ! GENERATES PARTICLES IN SPHERICAL DIST AROUND CENTER 
  write(*,*) "Add particles reached!"
  call setup_particles(x,np,radius,center,start,end)





 end subroutine add_particles


end module setup_binary
