module utils
 contains 

 function inner_product(center1,center2)
  real, intent(in) :: center1(3),center2(3)

  inner_product = center1(1)*center2(1) + center1(2)*center2(2) + center1(3) + center2(3)

 end function inner_product



