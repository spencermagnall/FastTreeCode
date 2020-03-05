program nbody
      use octree
      use setup_binary
      use step_leapfrog
      use output
      use poten
      use momentum
      use computemass
      implicit none 
      type(octreenode), allocatable :: nodes(:)
      integer, parameter :: nopart = 20
      real :: x(3, nopart)
      real :: v(3, nopart)
      real :: a(3, nopart)
      real :: m(nopart)
      integer :: i,iter,output_freq
      real :: rand
      real :: t, dt,tmax
      real :: center(3)
      real :: p(3)
      real:: pmag, pmagold, deltap 
      real :: angmom(3)
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
      dt = 10
      iter = 0 
      tmax = 10000*2*3.14159
      output_freq = 100 
      rootNode = 1
      sumMass = 0.0
      cm = 0.0
      ! PUT PARTICLE SETUP HERE
      call init(x,v,m,nopart)
      call maketree(nodes,x,v,a,nopart)
      call print_tree(nodes,x,0,1)
      !STOP
      call get_com(x,v,m,nopart,nodes,rootNode,sumMass,cm) 
      call print_tree(nodes,x,0,1)
      STOP
      !call get_accel(x,a,m,nopart,nodes)
      !STOP
      call write_output(x,v,m,nopart,t)
     ! STOP
       open(unit=66,file="Momentum",position="append")
       call get_momentum(x,v,m,p,angmom,nopart)
       pmag = dot_product(p,p)
       pmag = sqrt(pmag)
        write(66,*) "t "," pmag ", " angm x ", " angm y ", " angm z ", " deltap"
        write(66,*) t,  pmag, angmom, 0.0
       STOP
      do while(t < tmax)
            iter  = iter + 1

            call step(x,v,a,m,dt,nopart,nodes)

            t = t + dt
            write(*,*) t

            if (mod(iter,output_freq) .EQ. 0) then
                  call write_output(x,v,m,nopart,t)
            end if 
            pmagold = pmag
            call get_momentum(x,v,m,p,angmom,nopart)
            pmag = dot_product(p,p)
            pmag = sqrt(pmag)
            deltap = pmag - pmagold
            write(66,*) t,pmag, angmom, deltap

            !STOP            

      enddo 

end program nbody
