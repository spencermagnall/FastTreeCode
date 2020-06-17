program nbody
      use octree
      use delta_t
      use contrivedtree
      !use setup_binary
      use plummer_dist
      use step_leapfrog
      use output
      !use poten
      use momentum
      use computemass
      use opening_criterion
      use interaction
      use evaluate
      use testgravity
      use testdirect
      use potendirect
      implicit none 
      type(octreenode), allocatable :: nodes(:)
      integer, parameter :: nopart = 100
      real :: x(3, nopart)
      real :: v(3, nopart)
      real :: a(3, nopart)
      real :: atest(3,nopart)
      real :: poten(nopart)
      real :: m(nopart)
      integer :: iter,output_freq
      real :: t,dt,tmax,dtnew
      real :: center(3)
      real :: p(3),pold(3)
      real:: pmag, pmagold, deltap 
      real :: angmom(3),rmax
      integer :: rootnode,i,node1,node2
      real :: sumMass, cm(3)
      real :: c0,c1(3),c2(3,3),c3(3,3,3)
      integer :: interactionlist(nopart,nopart)


      
      c0 = 0.
      c1 = 0.
      c2 = 0.
      c3 = 0.

      call test_gravity()
      call test_coeff_trans()
      !STOP
      !x(:,1) = (/2.0,2.0,2.0/)
      !x(:,2) = (/1,1,1/)
      !x(:,3) = (/3.0,3.0,3.0/)
      !x(:,4) = (/4.0,4.0,4.0/)

      center = 0.0
      !call setup_particles(x,nopart,10.0,center,1,nopart)
     
      !do i=1, nopart
      !  print*, x(:,i)
      !enddo 
      
      

      !write(*,*) x(:,nodes(1) % data)
      !write(*,*) x(:,nodes(2) % data)
      

      !STOP
      t = 0
      dt = 0.01
      iter = 0 
      tmax = 10
      output_freq = 100 
      rootNode = 1
      sumMass = 0.0
      cm = 0.0
      a = 0.
      atest = 0.
      interactionlist = 0.
      
      ! PUT PARTICLE SETUP HERE
      call init(x,v,m,nopart)
      call get_accel_test(x,atest,m,nopart)
      print*, "Finished setup!"
      !call test_direct(x,m,nopart)
      !STOP
      !call write_output(x,v,a,m,nopart,t)
      !STOP
      !call maketreecontrived(nodes,x,v,a,nopart)
      call maketree(nodes,x,v,a,nopart)
      print*, "Finished tree build!"
      call print_tree(nodes,x,0,1)
      !STOP 
      call get_com(x,v,m,nopart,nodes,rootNode,sumMass,cm)
      ! Find rmax for each node
      rootNode = 1
      rmax = 0.0
      cm = 0.
      call find_rmax(x,nodes,rootNode,rmax)
      call print_tree(nodes,x,0,1)
      !STOP
      a = 0.
      node1 = 1
      node2 = 1
      call interact(node1,node2,nodes,x,m,a, nopart)
      !print*,nodes(1) % c3
      !print*, "Accel:"
      !print*, a
      !print*, nodes(1) % isLeaf
      !STOP
      c0 = 0.
      c1 = 0.
      c2 = 0.
      c3 = 0.
      cm = 0.
      call evaluate_gravity(nodes(1),nodes,cm,c0,c1,c2,c3,x,a)
      !call get_accel_body(nodes(1),nodes,x,a)
      print*,"Accel: "
      do i=1, nopart
        print*, a(:,i)
      enddo
      print*, "Accel test: "
      do i=1, nopart
        print*, atest(:,i)
      enddo
      print*,"Delta accel"
      do i=1,nopart
        print*,atest(:,i) - a(:,i)
      enddo 
      !STOP
      !call get_accel(x,a,m,nopart)
      !STOP
      ! get the timestep
      !call get_dtnew(0.1,a,dt,nopart)
      !STOP
      call write_output(x,v,a,m,nopart,t)
      !STOP
       open(unit=66,file="Momentum",status="replace")
       call get_momentum(x,v,m,p,angmom,nopart)
       pmag = dot_product(p,p)
       pmag = sqrt(pmag)
        write(66,*) "t "," pmag ", " angm x ", " angm y ", " angm z ", " deltap"
        write(66,*) t,  pmag, angmom, 0.0
      ! STOP
      do while(t < tmax)
            iter  = iter + 1

            call step(x,v,a,m,dt,nopart,nodes)

            t = t + dt
            write(*,*) t
            if (mod(iter,output_freq) .EQ. 0) then
                  if (t == 10000.OR. t ==11000) then
                      !call print_tree(nodes,x,0,1)
                  !    stop
                  endif
                  call write_output(x,v,a,m,nopart,t)
            end if 
            pold = p
            pmagold = pmag
            call get_momentum(x,v,m,p,angmom,nopart)
            pmag = dot_product(p,p)
            pmag = sqrt(pmag)
            deltap = sqrt(dot_product(p-pold,p-pold)) 
            write(66,*) t,pmag, angmom, deltap
            print*, deltap
            !call print_tree(nodes,x,0,1)
            !if (deltap > abs(1.e-14)) stop
            
            !STOP            

      enddo
      close(66)

end program nbody
