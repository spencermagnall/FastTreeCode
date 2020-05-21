program nbody
      use octree
      use contrivedtree
      use setup_binary
      !use plummer_dist
      use step_leapfrog
      use output
      !use poten
      use momentum
      use computemass
      use opening_criterion
      use interaction
      implicit none 
      type(octreenode), allocatable :: nodes(:)
      integer, parameter :: nopart = 20
      real :: x(3, nopart)
      real :: v(3, nopart)
      real :: a(3, nopart)
      real :: poten(nopart)
      real :: m(nopart)
      integer :: iter,output_freq
      real :: t, dt,tmax
      real :: center(3)
      real :: p(3)
      real:: pmag, pmagold, deltap 
      real :: angmom(3),rmax
      integer :: rootnode
      real:: sumMass, cm(3)
      
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
      dt = 0.00001
      iter = 0 
      tmax = 5
      output_freq = 100 
      rootNode = 1
      sumMass = 0.0
      cm = 0.0
      ! PUT PARTICLE SETUP HERE
      call init(x,v,m,nopart)
      print*, "Finished setup!"
      !call write_output(x,v,a,m,nopart,t)
      !STOP
      !call maketreecontrived(nodes,x,v,a,nopart)
      call maketree(nodes,x,v,a,nopart)
      print*, "Finished tree build!"
      !call print_tree(nodes,x,0,1)
      !STOP 
      call get_com(x,v,m,nopart,nodes,rootNode,sumMass,cm)
      ! Find rmax for each node
      rootNode = 1
      rmax = 0.0
      call find_rmax(x,nodes,rootNode,rmax)
      call print_tree(nodes,x,0,1)
      call interact(nodes(1),nodes(1),nodes,x,m,poten, nopart)
      STOP
      call get_accel(x,a,m,nopart)
      !STOP
      call write_output(x,v,a,m,nopart,t)
      !STOP
       open(unit=66,file="Momentum",status="replace")
       call get_momentum(x,v,m,p,angmom,nopart)
       pmag = dot_product(p,p)
       pmag = sqrt(pmag)
        write(66,*) "t "," pmag ", " angm x ", " angm y ", " angm z ", " deltap"
        write(66,*) t,  pmag, angmom, 0.0
       !STOP
      do while(t < tmax)
            iter  = iter + 1

            call step(x,v,a,m,dt,nopart)

            t = t + dt
            write(*,*) t

            if (mod(iter,output_freq) .EQ. 0) then
                  if (t == 10000.OR. t ==11000) then
                      !call print_tree(nodes,x,0,1)
                  !    stop
                  endif
                  call write_output(x,v,a,m,nopart,t)
            end if 
            pmagold = pmag
            call get_momentum(x,v,m,p,angmom,nopart)
            pmag = dot_product(p,p)
            pmag = sqrt(pmag)
            deltap = pmag - pmagold
            write(66,*) t,pmag, angmom, deltap

            !STOP            

      enddo
      close(66)

end program nbody
