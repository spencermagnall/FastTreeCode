module sphere_dist
 implicit none

 contains

  subroutine setup_particles(r,np,radius,center)
   ! Subroutine follows steps of the following algorithim:

   ! 10 Define an implict surface 
   ! 20 Use a (pseudo) random number generator to generate particle positions
   ! 30 If particle is within the surface add it, otherwise start again.
   ! 40 GOTO 10

   real, intent(out) :: r(:,:)
   integer, intent(in) :: np
   ! Radius of sphere
   real, intent(in) :: radius
   ! Center of the Sphere
   real, optional, intent(in) :: center(:)

   real :: x,y,z
   real :: implicteq
   integer :: i
   real :: c(3) 

   if(present(center)) then
    c = center
   else 
    c = 0.0
   endif

   i = 0 
   do while (i < np)
   	
   	write(*,*) i
   	x = RAND()
   	call negative_rand(x) 
    x =  x *radius + c(1)

    y =  RAND()
    call negative_rand(y)
    y = y * radius + c(2)

    z = RAND()
    call negative_rand(z)
    z = z * radius + c(3)
    


    implicteq = (x-c(1))**2 + (y-c(2))**2 + (z-c(3))**2

    if (implicteq <= radius**2) then
     i = i+1
     r(1,i) = x
     r(2,i) = y
     r(3,i) = z
     
    endif 
   enddo
   print*, "Setup: ", i ," particles" 
  end subroutine setup_particles

  subroutine negative_rand(input)
   real, intent(out) :: input
   input = 2*input - 1
   write(*,*) input  
  end subroutine negative_rand

end module sphere_dist  

