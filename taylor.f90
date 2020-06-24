module taylor_expansions
 contains 
 !-----------------------------------------------------------
!+
!  Compute the gravitational force between the node centres
!  along with the derivatives and second derivatives
!  required for the Taylor series expansions.
!+
!-----------------------------------------------------------
subroutine compute_fnode(dx,dy,dz,dr,totmass,quads,fnode)
 real, intent(in)    :: dx,dy,dz,dr,totmass
 real, intent(in)    :: quads(6)
 real, intent(inout) :: fnode(20)
 real :: dr3,dr4,dr5,dr6,dr3m,rx,ry,rz,qxx,qxy,qxz,qyy,qyz,qzz
 real :: dr4m3,rijQij,riQix,riQiy,riQiz,fqx,fqy,fqz
 real :: dfxdxq,dfxdyq,dfxdzq,dfydyq,dfydzq,dfzdzq
 real :: d2fxxxq,d2fxxyq,d2fxxzq,d2fxyyq,d2fxyzq
 real :: d2fxzzq,d2fyyyq,d2fyyzq,d2fyzzq,d2fzzzq

 ! note: dr == 1/sqrt(r2)
 print*, "dr: ", dr
 dr3  = dr*dr*dr
 print*, "dr3: "
 print*, dr3
 dr4  = dr*dr3
 dr5  = dr*dr4
 dr6  = dr*dr5
 dr3m  = totmass*dr3
 print*, "dr3m: ",dr3m
 dr4m3 = 3.*totmass*dr4
 rx  = dx*dr
 ry  = dy*dr
 rz  = dz*dr

 ! No quadrapole for the moment so == 1
 qxx = quads(1)
 qxy = quads(2)
 qxz = quads(3)
 qyy = quads(4)
 qyz = quads(5)
 qzz = quads(6)
 rijQij = (rx*rx*qxx + ry*ry*qyy + rz*rz*qzz + 2.*(rx*ry*qxy + rx*rz*qxz + ry*rz*qyz))
 riQix = (rx*qxx + ry*qxy + rz*qxz)
 riQiy = (rx*qxy + ry*qyy + rz*qyz)
 riQiz = (rx*qxz + ry*qyz + rz*qzz)
 fqx = dr4*(riQix - 2.5*rx*rijQij)
 fqy = dr4*(riQiy - 2.5*ry*rijQij)
 fqz = dr4*(riQiz - 2.5*rz*rijQij)
 dfxdxq = dr5*(qxx - 10.*rx*riQix - 2.5*rijQij   + 17.5*rx*rx*rijQij)
 dfxdyq = dr5*(qxy -  5.*ry*riQix - 5.0*rx*riQiy + 17.5*rx*ry*rijQij)
 dfxdzq = dr5*(qxz -  5.*rx*riQiz - 5.0*rz*riQix + 17.5*rx*rz*rijQij)
 dfydyq = dr5*(qyy - 10.*ry*riQiy - 2.5*rijQij   + 17.5*ry*ry*rijQij)
 dfydzq = dr5*(qyz -  5.*ry*riQiz - 5.0*rz*riQiy + 17.5*ry*rz*rijQij)
 dfzdzq = dr5*(qzz - 10.*rz*riQiz - 2.5*rijQij   + 17.5*rz*rz*rijQij)
 d2fxxxq = dr6*(-15.*qxx*rx + 105.*rx*rx*riQix - 15.*riQix - 157.5*rx*rx*rx*rijQij + 52.5*rx*rijQij)
 d2fxxyq = dr6*(35.*rx*rx*riQiy -  5.*qxx*ry - 5.*riQiy + 17.5*ry*rijQij - 157.5*rx*rx*ry*rijQij &
              + 70.*rx*ry*riQix - 10.*qxy*rx)
 d2fxxzq = dr6*(35.*rx*rx*riQiz -  5.*qxx*rz - 5.*riQiz + 17.5*rz*rijQij - 157.5*rx*rx*rz*rijQij &
              + 70.*rx*rz*riQix - 10.*qxz*rx)
 d2fxyyq = dr6*(70.*rx*ry*riQiy - 10.*qxy*ry - 5.*riQix + 17.5*rx*rijQij - 157.5*rx*ry*ry*rijQij &
              + 35.*ry*ry*riQix -  5.*qyy*rx)
 d2fxyzq = dr6*(35.*rx*ry*riQiz -  5.*qyz*rx  &
              + 35.*ry*rz*riQix -  5.*qxz*ry  &
              + 35.*rx*rz*riQiy -  5.*qxy*rz                             - 157.5*rx*ry*rz*rijQij)
 d2fxzzq = dr6*(70.*rx*rz*riQiz - 10.*qxz*rz - 5.*riQix + 17.5*rx*rijQij - 157.5*rx*rz*rz*rijQij &
              + 35.*rz*rz*riQix -  5.*qzz*rx)
 d2fyyyq = dr6*(-15.*qyy*ry + 105.*ry*ry*riQiy - 15.*riQiy - 157.5*ry*ry*ry*rijQij + 52.5*ry*rijQij)
 d2fyyzq = dr6*(35.*ry*ry*riQiz -  5.*qyy*rz - 5.*riQiz + 17.5*rz*rijQij - 157.5*ry*ry*rz*rijQij &
              + 70.*ry*rz*riQiy - 10.*qyz*ry)
 d2fyzzq = dr6*(70.*ry*rz*riQiz - 10.*qyz*rz - 5.*riQiy + 17.5*ry*rijQij - 157.5*ry*rz*rz*rijQij &
              + 35.*rz*rz*riQiy -  5.*qzz*ry)
 d2fzzzq = dr6*(-15.*qzz*rz + 105.*rz*rz*riQiz - 15.*riQiz - 157.5*rz*rz*rz*rijQij + 52.5*rz*rijQij)

 fnode( 1) = fnode( 1) - dx*dr3m + fqx ! fx
 fnode( 2) = fnode( 2) - dy*dr3m + fqy ! fy
 fnode( 3) = fnode( 3) - dz*dr3m + fqz ! fz
 fnode( 4) = fnode( 4) + dr3m*(3.*rx*rx - 1.) + dfxdxq ! dfx/dx
 fnode( 5) = fnode( 5) + dr3m*(3.*rx*ry)      + dfxdyq ! dfx/dy = dfy/dx
 fnode( 6) = fnode( 6) + dr3m*(3.*rx*rz)      + dfxdzq ! dfx/dz = dfz/dx
 fnode( 7) = fnode( 7) + dr3m*(3.*ry*ry - 1.) + dfydyq ! dfy/dy
 fnode( 8) = fnode( 8) + dr3m*(3.*ry*rz)      + dfydzq ! dfy/dz = dfz/dy
 fnode( 9) = fnode( 9) + dr3m*(3.*rz*rz - 1.) + dfzdzq ! dfz/dz
 fnode(10) = fnode(10) - dr4m3*(5.*rx*rx*rx - 3.*rx) + d2fxxxq ! d2fxdxdx
 fnode(11) = fnode(11) - dr4m3*(5.*rx*rx*ry - ry)    + d2fxxyq ! d2fxdxdy
 fnode(12) = fnode(12) - dr4m3*(5.*rx*rx*rz - rz)    + d2fxxzq ! d2fxdxdz
 fnode(13) = fnode(13) - dr4m3*(5.*rx*ry*ry - rx)    + d2fxyyq ! d2fxdydy
 fnode(14) = fnode(14) - dr4m3*(5.*rx*ry*rz)         + d2fxyzq ! d2fxdydz
 fnode(15) = fnode(15) - dr4m3*(5.*rx*rz*rz - rx)    + d2fxzzq ! d2fxdzdz
 fnode(16) = fnode(16) - dr4m3*(5.*ry*ry*ry - 3.*ry) + d2fyyyq ! d2fydydy
 fnode(17) = fnode(17) - dr4m3*(5.*ry*ry*rz - rz)    + d2fyyzq ! d2fydydz
 fnode(18) = fnode(18) - dr4m3*(5.*ry*rz*rz - ry)    + d2fyzzq ! d2fydzdz
 fnode(19) = fnode(19) - dr4m3*(5.*rz*rz*rz - 3.*rz) + d2fzzzq ! d2fzdzdz
 fnode(20) = fnode(20) - totmass*dr - 0.5*rijQij*dr3   ! potential


 !print*, "C1 correct"
 !print*, fnode(1:3)
 !print*, "C2 correct"
 !print*, fnode(4:9)
 

end subroutine compute_fnode

subroutine compute_coeff(dx,dy,dz,dr,totmass,quads,c0,c1,c2,c3)
 real, intent(in) :: dx,dy,dz,dr,totmass
 real, intent(in) :: quads(6)
 !real, intent(inout) :: coeff(4)
 real, intent(inout) :: c0,c1(3),c2(3,3),c3(3,3,3)
 real :: d0,d1
 real :: d2,d3
 real :: r
 !real :: d1arry(3)
 real :: rarry(3)
 real :: dr2,dr3,dr4,dr5,dr6,dr4m3
 integer :: i, j, k


 r = sqrt(dx**2 + dy**2 + dz**2)

 rarry = (/dx,dy,dz/)

 rx  = dx*dr
 ry  = dy*dr
 rz  = dz*dr

 ! Note dr = 1/r
 print*, "dr"
 d0 = dr
 dr2 = dr*dr
 dr3 = dr2*dr
 dr4 = dr3*dr
 dr5 = dr4*dr
 print*, d0
 d1 = -dr3
 print*, 'd1'
 print*, d1*totmass
 ! Why is this 3 not 2?????
 d2 = 3.*dr3*dr2
 print*, d2
 d3 = -5.*d2*dr2

 dr4m3 = 3.*totmass*dr4
 ! C0 = totmass * Greens function
 ! scalar 
 c0 = c0 + totmass * dr
 
 ! C1 = MB*Ri*D1
 ! Should be a vector
 !do i=1, 3
 !   c1(i) = c1(i) + totmass*rarry(i)*d1
 !enddo

 ! UNROLL LOOPS
 c1(1) = c1(1) + totmass*rarry(1)*d1
 c1(2) = c1(2) + totmass*rarry(2)*d1 
 c1(3) = c1(3) + totmass*rarry(3)*d1

 ! C2 = MB kronecker ij D1 + MB Ri Rj D2
 ! rank 2 tensor 
 !do j=1,3
 !  do i=1,3
 !     c2(i,j) = c2(i,j) + totmass*delta(i,j)*d1 +  totmass*rarry(i)*rarry(j)*d2
 !  enddo 
 !enddo 

 ! UNROLLLLLL

 !c2(1,1) = c2(1,1) + totmass*d1  - totmass*3.*rx*rx*d1
 !c2(1,2) = c2(1,2) - totmass*3.*rx*ry*d1
 !c2(1,3) = c2(1,3) - totmass*3.*rx*rz*d1
 !c2(2,1) = c2(2,1) - totmass*3.*ry*rx*d1
 !c2(2,2) = c2(2,2) + totmass*d1  - totmass*3.*d1*ry*ry
 !c2(2,3) = c2(2,3) - totmass*3.*ry*rx*d1
 !c2(3,1) = c2(3,1) - totmass*3.*rz*rx*d1
 !c2(3,2) = c2(3,2) - totmass*3.*rz*ry*d1
 !c2(3,3) = c2(3,3) + totmass*d1 - totmass*3.*rz*rz*d1
 
 ! C3
 ! rank3 tensor
 do k=1,3
    do j=1,3
       do i=1,3
         c3(i,j,k) =  c3(i,j,k) +  totmass*delta(i,j)*rarry(k)*d2 + totmass*delta(j,k)*rarry(i)*d2 &
          + totmass * delta(k,i)*rarry(j)*d2 + totmass*rarry(i)*rarry(j)*rarry(k)*d3
      enddo 
    enddo
 enddo 

 ! UNROLLLLLLLLLLLLLLLLLLL

 !c3(1,1,1) = c3(1,1,1) + totmass*(rarry(1) + rarry(1) + rarry(1))*d2 + totmass*rarry(1)*rarry(1)*rarry(1)*d3
 !c3(1,1,2) = c3(1,1,2) + totmass*(rarry(2))*d2 + totmass*rarry(1)*rarry(1)*rarry(2)*d3
 !c3(1,1,3) = c3(1,1,3) + totmass*(rarry(3))*d2 + totmass*rarry(1)*rarry(1)*rarry(3)*d3
 !c3(1,2,1) = c3(1,2,1) + totmass*(rarry(2))*d2 + totmass*rarry(1)*rarry(2)*rarry(1)*d3
 !c3(1,2,2) = c3(1,2,2) + totmass*(rarry(1))*d2 + totmass*rarry(1)*rarry(2)*rarry(2)*d3
 !c3(1,)

 !print*, "Coeff 0:"
 !print*, c0

 !print*, "Coeff 1:"
 !print*, c1

 !print*, "Coeff 2:"
 !print*, c2

 !print*, "Coeff 3:"
 !print*, c3 

 !print*, "dr is: "
 !print*,"Components are: "
 !print*, totmass*delta(2,1)*rarry(1)*d2
 !print*, totmass*delta(1,1)*rarry(2)*d2
 !print*, totmass*delta(1,2)*rarry(1)*d2
 !print*, dr


end subroutine compute_coeff

real function delta(i,j)
 integer, intent(in) :: i, j
 
 if (i == j) then
    delta = 1.
 else 
    delta = 0.
 endif 

end function delta

real function greens()

end function greens

end module taylor_expansions