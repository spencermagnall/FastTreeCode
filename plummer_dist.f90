module plummer_dist
 implicit none

 real, parameter :: pi = 4*ATAN(1.d0) 
 
 contains 
  subroutine init(x,v,m,np)
   integer, intent(in) :: np
   real, intent(out) :: x(3,np), v(3,np),m(np)
   real :: x1, x2, x3, x4, x5,x6,x7
   real :: r,rm,rm23,r2, az, sphr
   real :: xi,yi,zi
   real :: q, g,vr
   real :: totmass,mpart,vesc
   integer :: i
   logical :: accepted

   x1 = 0.0
   x2 = 0.0
   x3 = 0.0
   r = 0.0 
   az = 0.0
   sphr = 0.0
   xi = 0.0
   yi = 0.0
   zi = 0.0
   totmass = 1.0
   mpart = totmass/float(np)

   m(:) = mpart
   
   
   do i=1,np
    ! Generate Random numbers
    x1 = RAND()
    x2 = RAND()
    x3 = RAND()
    accepted = .FALSE.

    rm = x1*0.99
    rm23 = rm **(2./3.)
    r = sqrt(rm23/(1.-rm23))
    az = pi*(2.0*x2 - 1)
    sphr = ACOS(2.0*x3 - 1)

    ! convert to cartesian 
    call sphr_to_cart(r,az,sphr,xi,yi,zi)

    
    ! Set the position of particle 
    x(1,i) = xi
    x(2,i) = yi
    x(3,i) = zi


    do while (.NOT. accepted)
     x4 = RAND()
     x5 = RAND()
     q = x4
     g = g_q(q)

     if (0.1*x5 < g) then
      accepted = .TRUE.
     endif 
    enddo 

    r2 = dot_product(x(:,i),x(:,i))
    vesc  = sqrt(2.)*(1.+r2)**(-0.25)
    vr = q * vesc

    x6 = RAND()
    x7 = RAND()
    v(1,i) = (1.-2.*x6)*vr
    v(2,i) = sqrt(vr**2 - v(1,i)**2)*cos(2*pi*x7)
    v(3,i) = sqrt(vr**2 - v(1,i)**2)*sin(2*pi*x7)


   enddo 





  end subroutine init

  subroutine sphr_to_cart(r,az,sphr,x,y,z)
   real, intent(in) :: r, az, sphr
   real, intent(out) :: x, y, z

   x = r * sin(sphr)*cos(az)
   y = r * sin(sphr)*sin(az)
   z = r * cos(sphr)

  end subroutine sphr_to_cart

  real function g_q(q)
   real, intent(in) ::  q

   g_q = q**2*(1.0-q**2)

  end function g_q

end module plummer_dist 