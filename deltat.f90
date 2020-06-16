module delta_t

implicit none 

contains 
 real function get_delta_t(eps,partaccel) result(dt)
  real, intent(in) :: eps, partaccel(3)
  real :: eta, amag

  eta = 0.1
  amag = norm2(partaccel)

  dt = eta * sqrt(eps/amag)

 end function get_delta_t




subroutine get_dtnew(eps,accel,dt,np)
integer, intent(in) :: np
real, intent(in) :: eps, accel(3,np)
real, intent(inout) :: dt
real :: dtnew
integer :: i

 do i=1,np
         dtnew = get_delta_t(eps,accel(:,i))
         if (dtnew < dt) then
             dt = dtnew
         endif
 enddo
 print*, "dtnew: ", dt



end subroutine get_dtnew

end module delta_t