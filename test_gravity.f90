module testgravity
 use taylor_expansions
 use evaluate
implicit none 
contains 
 subroutine test_gravity()
  real :: c0,c1(3),c2(3,3),c3(3,3,3)
  real  :: xposi(3),xposj(3),x0(3)
  real :: dx(3),dr,quads(6),totmass
  real :: poten, phiexact
  real :: f0(3),fnode(20),fexact(3),accel(3)
  real :: startdan, stopdan,startwal,stopwal

  dx = 0.
  dr = 0.
  c0 = 0.
  c1 = 0.
  c2 = 0.0
  c3 = 0.0

  poten = 0.
  accel = 0.

  quads = 0.

  totmass = 5.
  f0 = 0.

  ! posiion of first particle
  xposi = (/0.05,-0.04,-0.05/)
  ! position of distant node 
  xposj = (/1.,1.,1./)
  ! pos of nearest node center
  x0 = 0.

  ! get seperations, vector and scalar
  call get_dx_dr(xposi,xposj,dx,dr)

  phiexact = -totmass*dr
  fexact = -(totmass*dr**3)*dx

  call get_dx_dr(x0,xposj,dx,dr)
  ! Compute tensor fields
  call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)
  call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)
  !print*, "C2 me"
  !print*, c2
  


  !call get_dx_dr(x0,xposj,dx,dr)
  ! Calculate the potential at xposi 

  call cpu_time(startwal)
  call poten_at_bodypos(xposi,x0,c0,c1,c2,c3,poten)
  call cpu_time(stopwal)

  call accel_at_bodypos(xposi,x0,c0,c1,c2,c3,accel)
  print*,'Phi exact ', phiexact
  print*, 'Force exact ', fexact
  print*,'Phi at origin (Walter) = ',-c0
  print*,'Phi with taylor series (Walter) = ', poten
  !print*, c1
  print*, 'Force (walter) = ', accel
  print*, stopwal - startwal, " Seconds"

  poten = 0.
  call get_dx_dr(x0,xposj,dx,dr)
  fnode = 0.
  quads = 0.
  call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)
  dx = xposi - x0   ! perform expansion about x0
  call cpu_time(startdan)
  call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),f0(1),f0(2),f0(3),poten)
  call cpu_time(stopdan)
  print*, 'Phi at origin (Daniel) =', fnode(20)
  print*, 'Phi with taylor series (Daniel) =',poten
  print*, "Force daniel: ", f0
  print*, stopdan - startdan, " Seconds"

 end subroutine test_gravity

 subroutine test_coeff_trans()

  real :: c0,c1(3),c2(3,3),c3(3,3,3)
  real  :: xposi(3),xposj(3),x0(3)
  real :: dx(3),dr,quads(6),totmass
  real :: poten, phiexact
  real :: f0(3),fnode(20),fexact(3),accel(3)
  real :: startdan, stopdan,startwal,stopwal
  real :: c0ex, c1ex(3),c2ex(3,3),c3ex(3,3,3)

  dx = 0.
  dr = 0.
  c0 = 0.
  c1 = 0.
  c2 = 0.0
  c3 = 0.0
  c0ex = 0.
  c1ex = 0.
  c2ex = 0.
  c3ex = 0.

  poten = 0.
  accel = 0.

  quads = 0.

  totmass = 5.
  f0 = 0.

  ! posiion of first particle
  xposi = (/0.05,-0.04,-0.05/)
  ! position of distant node 
  xposj = (/1.,1.,1./)
  ! pos of nearest node center
  x0 = 0.

  ! compute the exact value of the coeffs 
  ! get separations, vector and scalar
  call get_dx_dr(xposi,xposj,dx,dr)
  call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0ex,c1ex,c2ex,c3ex)



  call get_dx_dr(x0,xposj,dx,dr)
  ! Compute tensor fields
  call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)
  call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)
  dx = xposi - x0
  print*,"Before transform: "
  !print*,c0
  print*,c1
  print*,c2
  !print*,c3
  call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),f0(1),f0(2),f0(3),poten)
  call translate_expansion_center(x0,xposi,c0,c1,c2,c3)

  print*,"Exact values: "
  print*,c0ex
  print*,c1ex
  print*,c2ex
  print*,c3ex

  print*, "Translated values: "
  print*, c0
  print*,c1
  print*,c2
  print*,c3

  print*, "Dan accel (c1) :"
  print*,f0

  print*, "delta values: "
  print*, c0-c0ex
  print*, c1-c1ex
  print*,c2-c2ex
  print*,c3-c3ex

  print*, "Reverse transform: "
  call translate_expansion_center(xposi,x0,c0,c1,c2,c3)
  !print*,c0
  print*,c1
  print*,c2
  !print*,c3


  ! test symmetric translation 

   dx = 0.
  dr = 0.
  c0 = 0.
  c1 = 0.
  c2 = 0.0
  c3 = 0.0
  c0ex = 0.
  c1ex = 0.
  c2ex = 0.
  c3ex = 0.
  ! compute the exact value of the coeffs 
  ! get separations, vector and scalar
  call get_dx_dr(xposj,xposi,dx,dr)
  call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0ex,c1ex,c2ex,c3ex)

  call get_dx_dr(xposi,xposj,dx,dr)
  call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)

  !print*, "Coefficent symm: "
  print*,abs(c0ex) - abs(c0)
  print*,abs(c1ex) - abs(c1)
  !print*,c2ex
  !print*,c3ex
  !print*,c0
  !print*,c1
  !print*,c2
  !print*,c3






 end subroutine test_coeff_trans

 subroutine test_trans_error()

  real :: c0,c1(3),c2(3,3),c3(3,3,3)
  real  :: xposi(3),xposj(3),x0(3),x1(3),xposii(3),xposjj(3),xposiii(3)
  real :: dx(3),dr,quads(6),totmass,totmass1
  real :: poten, phiexact
  real :: f0(3),fnode(20),fexact(3),accel(3)
  real :: startdan, stopdan,startwal,stopwal
  real :: c0ex, c1ex(3),c2ex(3,3),c3ex(3,3,3),asum(3), fnode1(3),fnode2(3)
  real :: xarray(3,10)
  integer :: i

  dx = 0.
  dr = 0.
  c0 = 0.
  c1 = 0.
  c2 = 0.0
  c3 = 0.0
  c0ex = 0.
  c1ex = 0.
  c2ex = 0.
  c3ex = 0.

  poten = 0.
  accel = 0.

  quads = 0.
  asum = 0.

  fnode1 = 0.
  fnode2 = 0.
  xarray = 0.

  totmass = 10.
  totmass1 = 15.
  f0 = 0.
  ! pos of nearest node center
  x0 = 0.
  ! pos of second node center 
  x1 = 0. 

  ! posiion of first particle
  xposi = (/0.05,-0.04,-0.05/)
  xposii = (/0.1, -0.22, -0.17/)
  xposiii = (/0.09,0.1,-0.15/)
  x0 = (xposi*5. + xposii*5. + xposiii*5.)/totmass1
  print*, "Node 1 com: ",x0
  ! position of distant particle 
  xposj = (/1.05,0.96,0.95/)
  xposjj = (/1.04, 0.85,0.97/)
  x1 = (xposj*5. + xposjj*5.)/totmass
  print*, "Node 2 com: ",x1
  
  ! GET QUADS FOR NODE 2
  xarray(:,1) = xposj
  xarray(:,2) = xposjj

  do i=1, 2
    print*, i
    dx = xarray(:,i) - x1
       
    ! UNROLLED QUADS, DONT NEED ALL OF THESE BECAUSE OF SYM
    quads(1) = quads(1) + 5.*(3.*dx(1)*dx(1) - dot_product(dx,dx))
    quads(2) = quads(2) + 5.*(3.*dx(1)*dx(2))
    quads(3) = quads(3) + 5.*(3.*dx(1)*dx(3))
       !quads(2,1) = quads(2,1) + m(j)*(3.*dx(2)*dx(1))
    quads(4) = quads(4) + 5.*(3.*dx(2)*dx(2) - dot_product(dx,dx))
    quads(5) = quads(5) + 5.*(3.*dx(2)*dx(3))
    quads(6) = quads(6) + 5.*(3.*dx(3)*dx(3) - dot_product(dx,dx)) 

  enddo

  print*, "Quads: ", quads
  !stop  



  ! calculate exact accel on particle 
  call get_dx_dr(xposi,x1,dx,dr)
  fexact = -(totmass*dr**3)*dx
   ! calculate field tensors between nodes 
  call get_dx_dr(x0,x1,dx,dr)
  ! Compute tensor fields
  fnode = 0.
  call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)
  call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)

  print*, "c1 node: ",c1
  print*, "Force node 1: ",c1*totmass1
  dx = xposi - x0
  ! translate to particle 
  call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),f0(1),f0(2),f0(3),poten)
  call translate_expansion_center(x0,xposi,c0,c1,c2,c3)
  !call accel_at_bodypos(x0,xposi,c0,c1,c2,c3,c1)
  print*, "Particle 1"
  print*, "Accel: ", c1
  print*, "Exact accel: ", fexact
  print*, "Fgrav: ", f0
  print*, "Delta: ", c1 - fexact
  asum = asum + c1 
  fnode1 = fnode1 + c1*5. 


  c0 = 0.
  c1 = 0.
  c2 = 0.0
  c3 = 0.0
    ! calculate field tensors between nodes 
  call get_dx_dr(xposii,x1,dx,dr)
  fexact = -(totmass*dr**3)*dx
    ! calculate field tensors between nodes 
  call get_dx_dr(x0,x1,dx,dr)
  ! Compute tensor fields
  fnode = 0.
  !dx = xposii - x0
  call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)
  call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)

  print*, "c1 node: ",c1 

  ! translate to second particle 
  f0 = 0.
  dx = xposii - x0
  call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),f0(1),f0(2),f0(3),poten)
  call translate_expansion_center(x0,xposii,c0,c1,c2,c3)
  !call accel_at_bodypos(x0,xposii,c0,c1,c2,c3,c1)
  print*, "Particle 2"
  print*, "Accel: ", c1
  print*, "Exact accel: ", fexact
  print*, "Fgrav: ", f0
  print*, "Delta: ", c1 - fexact
  asum = asum + c1
  fnode1 = fnode1 + c1*5.


  c0 = 0.
  c1 = 0.
  c2 = 0.0
  c3 = 0.0
  fnode = 0.
  ! calculate field tensors between nodes 
  call get_dx_dr(xposiii,x1,dx,dr)
  fexact = -(totmass*dr**3)*dx
    ! calculate field tensors between nodes 
  call get_dx_dr(x0,x1,dx,dr)
  ! Compute tensor fields
  fnode = 0.
  !dx = xposii - x0
  call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)
  call compute_fnode(dx(1),dx(2),dx(3),dr,totmass,quads,fnode)

  print*, "c1 node: ",c1 

  ! translate to second particle 
  f0 = 0.
  dx = xposiii - x0 
  call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),f0(1),f0(2),f0(3),poten)
  call translate_expansion_center(x0,xposiii,c0,c1,c2,c3)
  !call accel_at_bodypos(x0,xposiii,c0,c1,c2,c3,c1)
  print*, "Particle 5"
  print*, "Accel: ", c1
  print*, "Exact accel: ", fexact
  print*, "Fgrav: ", f0
  print*, "Delta: ", c1 - fexact
  asum = asum + c1
  fnode1 = fnode1 + c1*5. 
  fnode1 = fnode1!*totmass

  !fnode1 = fnode1 * totmass
 
  c0 = 0.
  c1ex = c1
  c1 = 0.
  c2 = 0.0
  c3 = 0.0
  fnode = 0.
  quads = 0.
  xarray = 0.
  xarray(:,1) = xposi
  xarray(:,2) = xposii
  xarray(:,3) = xposiii

  do i=1, 3
    print*, i
    dx = xarray(:,i) - x0
       
    ! UNROLLED QUADS, DONT NEED ALL OF THESE BECAUSE OF SYM
    quads(1) = quads(1) + 5.*(3.*dx(1)*dx(1) - dot_product(dx,dx))
    quads(2) = quads(2) + 5.*(3.*dx(1)*dx(2))
    quads(3) = quads(3) + 5.*(3.*dx(1)*dx(3))
       !quads(2,1) = quads(2,1) + m(j)*(3.*dx(2)*dx(1))
    quads(4) = quads(4) + 5.*(3.*dx(2)*dx(2) - dot_product(dx,dx))
    quads(5) = quads(5) + 5.*(3.*dx(2)*dx(3))
    quads(6) = quads(6) + 5.*(3.*dx(3)*dx(3) - dot_product(dx,dx)) 

  enddo

  print*, "Quads: ", quads


  call get_dx_dr(xposj,x0,dx,dr)
   ! calculate field tensors between nodes 
  fexact = -(totmass*dr**3)*dx
  call get_dx_dr(x1,x0,dx,dr)
 
  ! Compute tensor fields
  call compute_coeff(dx(1),dx(2),dx(3),dr,totmass1,quads,c0,c1,c2,c3)
  call compute_fnode(dx(1),dx(2),dx(3),dr,totmass1,quads,fnode)


  ! translate to particle 
  f0 = 0.
  dx = xposj - x1
  call translate_expansion_center(x1,xposj,c0,c1,c2,c3)
  print*, "Particle 3"
  print*, "Accel: ", c1
  print*, "Exact accel: ", fexact
  call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),f0(1),f0(2),f0(3),poten)
  print*, "Fgrav: ", f0
  print*, "Delta: ", c1 - fexact
  asum = asum + c1
  fnode2 = fnode2 + c1*5.


  c0 = 0.
  c1ex = c1
  c1 = 0.
  c2 = 0.0
  c3 = 0.0
  fnode = 0.

  call get_dx_dr(xposjj, x0,dx,dr)
   ! calculate field tensors between nodes 
   fexact = -(totmass*dr**3)*dx
  call get_dx_dr(x1,x0,dx,dr)
  
  ! Compute tensor fields
  call compute_coeff(dx(1),dx(2),dx(3),dr,totmass1,quads,c0,c1,c2,c3)
  call compute_fnode(dx(1),dx(2),dx(3),dr,totmass1,quads,fnode)

  print*, "Force node 2: ", c1*totmass

  ! translate to particle 
  f0 = 0.
  dx = xposjj - x1
  call translate_expansion_center(x1,xposjj,c0,c1,c2,c3)
  print*, "Particle 4"
  print*, "Accel: ", c1
  print*, "Exact accel: ", fexact
  call expand_fgrav_in_taylor_series(fnode,dx(1),dx(2),dx(3),f0(1),f0(2),f0(3),poten)
  print*, "Fgrav: ", f0
  print*, "Delta: ", c1 - fexact
  asum = asum + c1
  fnode2 = fnode2 + c1*5.

  fnode2 = fnode2 !*totmass1


 

  print*, "Net sum: ", asum 

  print*, "Force on node 1: ", fnode1
  print*, "Force on node 2: ", fnode2 
  print*, "Delta between nodes: ", abs(fnode1) - abs(fnode2)



  print*, "START OF ACTUAL TEST"


  dx = 0.
  dr = 0.
  c0 = 0.
  c1 = 0.
  c2 = 0.0
  c3 = 0.0
  c0ex = 0.
  c1ex = 0.
  c2ex = 0.
  c3ex = 0.

  poten = 0.
  accel = 0.

  quads = 0.
  asum = 0.

  fnode1 = 0.
  fnode2 = 0.

  !totmass = 10.
  !totmass1 = 15.
  f0 = 0.
  ! pos of nearest node center
  x0 = 0.
  ! pos of second node center 
  x1 = 0. 

  c1 = (/-5.8985004080591208E-002, -7.4183265630189205E-002, -4.4014837786063528E-002/)
  c2 = (reshape((/-9.8531906134372016E-002, 5.8153958731058664E-002, 3.2276054832296362E-002, &
   5.8153958731058664E-002, -7.0359705378363752E-002, 3.7978866563416674E-002, 3.2276054832296362E-002, &
   5.8153958731058664E-002, -0.10600313922977826/), shape(c2)))
  x0 = (/0.71381594737370813, 0.77321779727935791, 0.80617229143778490/)
  xposi = (/0.71380925178527821, 0.99330163002014160, 0.92149341106414795/)
  xposii = (/ 0.68923449516296376, 0.63332358996073390, 0.93740367889404286/)
  xposiii = (/0.72435276848929275, 0.77029136248997299, 0.71698137692042763/)

  print*, "C1 node 1: "
  call translate_expansion_center(x0,xposi,c0,c1,c2,c3)
  print*, c1 
  print*, "Force node 1: ", c1*0.025
  asum = asum + c1*0.025

  c1 = (/-5.8985004080591208E-002, -7.4183265630189205E-002, -4.4014837786063528E-002/)
  c2 = (reshape((/-9.8531906134372016E-002, 5.8153958731058664E-002, 3.2276054832296362E-002, &
   5.8153958731058664E-002, -7.0359705378363752E-002, 3.7978866563416674E-002, 3.2276054832296362E-002, &
  5.8153958731058664E-002, -0.10600313922977826/), shape(c2)))

  print*, "C1 node 2: "
  call translate_expansion_center(x0,xposii,c0,c1,c2,c3)
  print*, c1 
  print*, "Force node 2: ", c1*0.0375
  asum = asum + c1*0.0375

  print*, "C1 node 3: "

  c1 = (/-5.8985004080591208E-002, -7.4183265630189205E-002, -4.4014837786063528E-002/)
  c2 = (reshape((/-9.8531906134372016E-002, 5.8153958731058664E-002, 3.2276054832296362E-002, &
   5.8153958731058664E-002, -7.0359705378363752E-002, 3.7978866563416674E-002, 3.2276054832296362E-002, &
  5.8153958731058664E-002, -0.10600313922977826/), shape(c2)))

  call translate_expansion_center(x0,xposiii,c0,c1,c2,c3)
  print*, c1
  print*, "Force node 3: ", c1*0.0875
  asum = asum + c1*0.0875


  print*, "Force sum is: ", asum
  c1 = (/-5.8985004080591208E-002, -7.4183265630189205E-002, -4.4014837786063528E-002/)
  print*, "Force node 1 is: ", c1*0.15

  print*, "Delta: ", asum - c1*0.15

 end subroutine test_trans_error


 subroutine get_dx_dr(x1,x2,dx,dr)
  real, intent(in) :: x1(3),x2(3)
  real, intent(out) :: dx(3),dr

  dx = x1 - x2
  dr = 1./sqrt(dot_product(dx,dx))

 end subroutine get_dx_dr

end module testgravity