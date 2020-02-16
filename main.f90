program nbody
      use octree
      use sphere_dist
      implicit none 
      type(octreenode), allocatable :: nodes(:)
      integer, parameter :: nopart = 2000000
      real :: x(3, nopart)
      real :: v(3, nopart)
      real :: a(3, nopart)
      integer :: i
      real :: rand
      !x(:,1) = (/2.0,2.0,2.0/)
      !x(:,2) = (/1,1,1/)
      !x(:,3) = (/3.0,3.0,3.0/)
      !x(:,4) = (/4.0,4.0,4.0/)
      rand = 0.0
      call negative_rand(rand)
      rand = 0.5
      call negative_rand(rand)
      rand = 1.0
      call negative_rand(rand)
      call setup_particles(x,nopart,10.0)
     
      !do i=1, nopart
      !  print*, x(:,i)
      !enddo 
      
      call maketree(nodes,x,v,a,nopart)

      !write(*,*) x(:,nodes(1) % data)
      !write(*,*) x(:,nodes(2) % data)
      !call print_tree(nodes,x,0,1)
end program nbody
