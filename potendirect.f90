module potendirect
 use octreetype
 implicit none 
 contains 

  subroutine get_accel (x,a,m,np,particlesindex1,particlesindex2)
  ! the index's of particles (BODIES)
  integer, intent(in) :: particlesindex1(10), particlesindex2(:), np
  real, intent(in) :: x(3,np) 
  real, intent(inout) :: a(3, np)
  real, intent(in) :: m(np)
  real :: r,r2,dx(3), h
  integer :: i, j
  integer :: indexi,indexj,part2size
  real , allocatable :: particles(:) 

  

  ! THIS IS BAD 
  if (np > 128) then
    part2size = size(particlesindex2)
  else 
    part2size = np
  endif 

  allocate(particles(part2size+10))
  particles = [particlesindex1,particlesindex2]

  print*, "Working"
  !print*, x(:,:)

  ! Put all the particles into one list 

  h = 0.
  ! in this case we don't want to reset a
  ! but is should be done at the start of each step
  ! THIS IS WRONG!!
  do i=1, size(particles)
   indexi = particles(i)
   do j=1, size(particles)
    indexj = particles(j)
    if (indexj/=indexi .AND. indexi /= 0 .AND. indexj /= 0) then
     print*, indexi, indexj
     dx = x(:,indexi) - x(:,indexj)
     !print*,"dx: ", dx
     !print*, "r2: ", r2
     !print*, "m: ", m(indexj) 
     !print*, "abefore: ", a(:,indexi)
     r2 = dot_product(dx,dx)
     r  = sqrt(r2)
     a(:,indexi) = a(:,indexi) - m(indexj)*dx*(1/((r2 + h**2)**1.5))
     !print*,"Accel direct: ", a(:,indexi)
     if (isnan(a(2,indexi))) stop 
    endif
  enddo
 enddo
end subroutine get_accel

subroutine get_poten (x,poten,m,np,particlesindex1,particlesindex2)
  ! the index's of particles (BODIES)
  integer, intent(in) :: particlesindex1(10), particlesindex2(10), np
  real, intent(in) :: x(3,np)
  ! technially this is potential 
  real, intent(inout) :: poten(np)
  real, intent(in) :: m(np)
  real :: r,r2,dx(3), h
  integer :: i, j
  integer :: indexi,indexj

  print*, "Working"
  !print*, x(:,:)

  h = 50.0
  ! in this case we don't want to reset a
  do i=1, 10
   indexi = particlesindex1(i)
   do j=1, 10
    indexj = particlesindex2(j)
    if (j/=i .AND. i /= 0 .AND. j /= 0) then
     print*, indexi, indexj
     dx = x(:,indexi) - x(:,indexj)
     r2 = dot_product(dx,dx)
     r  = sqrt(r2)
     poten(indexi) = poten(indexi) - m(indexj)*(1/((r2 + h**2)))
    endif
  enddo
 enddo 
end subroutine get_poten



subroutine get_poten_node(x,poten,m,np,particlesindex,node)
 integer, intent(in) :: particlesindex(10),np
 real, intent(in) :: x(3,np)
 real, intent(inout) :: poten(np)   
 real, intent(in) :: m(np)
 type(octreenode) :: node
 real :: nodecm(3), nodemass
 real :: dx,r2,r
 integer :: i,indexi 
 real :: h 
 h = 50.0

 nodecm = node % centerofmass
 nodemass = node % totalmass

 do i=1,10
  indexi = particlesindex(i)


enddo 

end subroutine get_poten_node

end module potendirect
