module cold_collapse
 implicit none 
 contains 

subroutine init(x,v,m,np)
 integer, intent(in) :: np
 real, intent(out) :: x(3,np), v(3,np), m(np)
 integer :: i


 ! Mass of unity 
 m = 1.0/np

 ! particles at rest 
 v = 0.
 


 do i=1, np
    x(1,i) = RAND()
    x(2,i) = RAND()
    x(3,i) = RAND()
 enddo 



 
end subroutine init  
end module cold_collapse