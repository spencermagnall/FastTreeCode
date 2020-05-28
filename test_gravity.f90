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
  real :: f0(3),fnode(20)
  real :: startdan, stopdan,startwal,stopwal

  c0 = 0.
  c1 = 0.
  c2 = 0.0
  c3 = 0.0

  poten = 0.

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

  call get_dx_dr(x0,xposj,dx,dr)
  ! Compute tensor fields
  call compute_coeff(dx(1),dx(2),dx(3),dr,totmass,quads,c0,c1,c2,c3)
  !print*, c0
  


  !call get_dx_dr(x0,xposj,dx,dr)
  ! Calculate the potential at xposi 

  call cpu_time(startwal)
  call poten_at_bodypos(xposi,x0,c0,c1,c2,c3,poten)
  call cpu_time(stopwal)

  print*,'Phi exact ', phiexact
  print*,'Phi at origin (Walter) = ',-c0
  print*,'Phi with taylor series (Walter) = ', poten
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
  print*, stopdan - startdan, " Seconds"

 end subroutine test_gravity


 subroutine get_dx_dr(x1,x2,dx,dr)
  real, intent(in) :: x1(3),x2(3)
  real, intent(out) :: dx(3),dr

  dx = x1 - x2
  dr = 1./sqrt(dot_product(dx,dx))

 end subroutine get_dx_dr

end module testgravity