module potendirect
 use octreetype
 implicit none 
 contains 

  subroutine get_accel (x,a,m,np,particlesindex1)
  ! the index's of particles (BODIES)
  integer, intent(in) :: particlesindex1(10), np
  real, intent(in) :: x(3,np) 
  real, intent(inout) :: a(3, np)
  real, intent(in) :: m(np)
  real :: r,r2,dx(3), h
  integer :: i, j
  integer :: indexi,indexj,part2size
  integer , allocatable :: particles(:) 

   ! TODO REFACTOR TO JUST TAKE IN ARRAY CONTIAINING PARTICLE INDEXS
  
   !h = 1.
  ! THIS IS BAD 
  !if (np > 128) then
    !part2size = size(particlesindex2)
  !else 
    !part2size = np
  !endif 
  
  !print*,"particles: ", particles 

  !print*, "Direct Sum"
  !print*, x(:,:)

  !open(unit=56,file="interactions.txt")
  ! Put all the particles into one list 

  h = 0.1
  ! in this case we don't want to reset a
  ! but is should be done at the start of each step
  ! THIS IS WRONG!!
  do i=1, 10
   indexi = particlesindex1(i)
   do j=1, 10
    indexj = particlesindex1(j)
    if (indexj/=indexi .AND. indexi /= 0 .AND. indexj /= 0) then
     !print*, "indexi, indexj"
     !print*, indexi, indexj
     dx = x(:,indexi) - x(:,indexj)
     !print*,"dx: ", dx
     !print*, "r2: ", r2
     !print*, "m: ", m(indexj) 
     !print*, "abefore: ", a(:,indexi)
     r2 = dot_product(dx,dx)
     r  = sqrt(r2)
     !$omp critical (accel)
     a(:,indexi) = a(:,indexi) - m(indexj)*(1/((r2 + h**2)**1.5))*dx
     !$omp end critical (accel)
     !print*,"Accel direct: ", a(:,indexi)
     !print*, "A mag: ", sqrt(dot_product(m(indexj)*(1/((r2*r)))*dx, m(indexj)*(1/((r2*r)))*dx))
     !if (isnan(a(2,indexi))) stop 
     !write(56,*) "Particle1 ",indexi, " Particle2 ",indexj, "Accel: ", m(indexj)*(1/((r2 + h**2)**1.5))*dx
    endif
  enddo
 enddo


 


end subroutine get_accel

subroutine get_accel_leafnode(x,a,m,np,particlesindex1,particlesindex2)
  integer, intent(in) :: np, particlesindex1(:), particlesindex2(:)
  real, intent(in) :: x(3,np), m(:)
  real, intent(inout) :: a(3,np)
  real :: r,r2,dx(3), h
  integer :: i, j
  integer :: indexi,indexj

   !h = 1.
  h = 0.1
  ! in this case we don't want to reset a
  ! but is should be done at the start of each step
  ! THIS IS WRONG!!

  !print*, size(particlesindex1)
  !print*, size(particlesindex2)
  !print*, "Size of x: ", size(x)
  !print*, "Size of a: ", size(a)
  do i=1, size(particlesindex1)
   !print*, i
   indexi = particlesindex1(i)
   do j=1, size(particlesindex2)
    !print*, j 
    indexj = particlesindex2(j)
    if (indexj/=indexi .AND. indexi /= 0 .AND. indexj /= 0) then
     !print*, "indexi, indexj"
     !print*, indexi, indexj
     !print*,"i, j: ", i, j
     !print*, "Nopart: ",np
     dx = x(:,indexi) - x(:,indexj)
     !print*, "Segfault here "
     !print*,"dx: ", dx
     !print*, "r2: ", r2
     !print*, "m: ", m(indexj) 
     !print*, "abefore: ", a(:,indexi)
     r2 = dot_product(dx,dx)
     r  = sqrt(r2)
     !$omp critical (accel)
     a(:,indexi) = a(:,indexi) - m(indexj)*(1/((r2 + h**2)**1.5))*dx
     !$omp end critical (accel)
     !print*, "Second segfault here "
     !print*,"Accel direct: ", a(:,indexi)
     !print*, "A mag: ", sqrt(dot_product(m(indexj)*(1/((r2*r)))*dx, m(indexj)*(1/((r2*r)))*dx))
     !if (isnan(a(2,indexi))) stop 
     !write(56,*) "Particle1 ",indexi, " Particle2 ",indexj, "Accel: ", m(indexj)*(1/((r2 + h**2)**1.5))*dx
    endif
  enddo
 enddo

 !print*, "Finshed direct sum"


end subroutine get_accel_leafnode

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

subroutine get_accel_test(x,a,m,np)
  integer, intent(in) :: np
  real, intent(in) :: x(3,np),m(np)
  real, intent(inout) :: a(3,np)
  real :: h,dx(3),r2,r
  integer :: i,j

  !h = 0.
  h = 0.1

  ! SUBOUTINE FOR TESTING THAT ACCEL IS COMPUTED CORRECTLY 
  a = 0.
  do i=1, np
   do j=1, np
    if (i/=j) then
     !print*, indexi, indexj
     dx = x(:,i) - x(:,j)
     !print*,"dx: ", dx
     !print*, "r2: ", r2
     !print*, "m: ", m(indexj) 
     !print*, "abefore: ", a(:,indexi)
     r2 = dot_product(dx,dx)
     r  = sqrt(r2)
     a(:,i) = a(:,i) - m(j)*dx*(1/((r2 + h**2)**1.5))
     !print*,"Accel direct: ", a(:,indexi)
    endif
  enddo
 enddo

end subroutine get_accel_test

end module potendirect
