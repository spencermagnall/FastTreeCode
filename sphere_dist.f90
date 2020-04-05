module sphere_dist
 implicit none

 contains

  subroutine setup_particles(r,np,radius,center,start,end)
   ! Subroutine follows steps of the following algorithim:

   ! 10 Define an implict surface 
   ! 20 Use a (pseudo) random number generator to generate particle positions
   ! 30 If particle is within the surface add it, otherwise start again.
   ! 40 GOTO 10

   real, intent(out) :: r(3,np)
   integer, intent(in) :: np
   integer, optional, intent(in) :: start, end 
   ! Radius of sphere
   real, intent(in) :: radius
   ! Center of the Sphere
   real, optional, intent(in) :: center(:)

   real :: x,y,z
   real :: implicteq
   integer :: i
   real :: c(3) 

   x = 0.0
   y = 0.0
   z = 0.0

   if(present(center)) then
    c = center
   else 
    c = 0.0
   endif

   i = start - 1 

   write(*,*) i

   open(unit=67, file="dist", position='append') 
   do while (i < end)
   
    write(*,*) i
    x = RAND()
    write(*,*) x
    call negative_rand(x) 
    x =  x *radius + c(1)

    y =  RAND()
    write(*,*) y
    call negative_rand(y)
    y = y * radius + c(2)

    z = RAND()
    write(*,*) z
    call negative_rand(z)
    z = z * radius + c(3)
   
    write(*,*) "X is = ", x
    write(*,*) " Y is = ", y
    write(*,*) "Z is = ", y


    implicteq = (x-c(1))**2 + (y-c(2))**2 + (z-c(3))**2

    if (implicteq <= radius**2) then
     i = i+1
     r(1,i) = x
     r(2,i) = y
     r(3,i) = z
     write(67,*) r(:,i)
     write(*,*) r(:,i) 
     
    endif 
   enddo
   print*, "Setup: ", i ," particles" 
   close(67)
  end subroutine setup_particles

  subroutine negative_rand(input)
   real, intent(inout) :: input
   input = 2*input - 1
   !write(*,*) input  
  end subroutine negative_rand

end module sphere_dist  

