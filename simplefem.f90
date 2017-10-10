!==============================================!
! C. thieulot ; May 2011                       !
! Version by J. Mos ; October 2017             !
!==============================================!
                                               !
program fcubed                                 !
                                               !
#include <petsc/finclude/petscksp.h>           
use petscksp                                   !
                                               !
implicit none                                  !
                                               !
integer, parameter :: m=4                      ! number of nodes which constitute an element
integer, parameter :: ndof=2                   ! number of dofs per node
integer nnx                                    ! number of grid points in the x direction
integer nny                                    ! number of grid points in the y direction
integer np                                     ! number of grid points
integer nelx                                   ! number of elements in the x direction
integer nely                                   ! number of elements in the y direction
integer nel                                    ! number of elements
integer Nfem                                   ! size of the FEM matrix 
integer NZ                                     ! total nb of nonzeroes
integer, dimension(:,:), allocatable :: icon   ! connectivity array
integer, dimension(:), allocatable :: ipvt     ! work array needed by the solver 
integer, dimension(:), allocatable :: csr_ia   ! rowpointer for csr storage
integer, dimension(:), allocatable :: csr_ja   ! csr storage array
integer, dimension(:), allocatable :: csr_ia_p   ! rowpointer for csr storage
integer, dimension(:), allocatable :: csr_ja_p   ! csr storage array
integer, dimension(:), allocatable :: nb_nz    ! nb of nonzeroes per matrix row
                                               !
integer i1,i2,j1,j2,i,j,k,iel,counter,iq,jq,ii,jj,ip,jp    !
integer ik,jk,ikk,jkk,m1,m2,k1,k2,job,jfake(1),l,nsees !
integer ij,inode
                                               !  
real(8) Lx,Ly                                  ! size of the numerical domain
real(8) viscosity                              ! dynamic viscosity $\mu$ of the material
real(8) density                                ! mass density $\rho$ of the material
real(8) gx,gy                                  ! gravity acceleration
real(8) penalty                                ! penalty parameter lambda
real(8), dimension(:),   allocatable :: x,y    ! node coordinates arrays
real(8), dimension(:),   allocatable :: u,v    ! node velocity arrays
real(8), dimension(:),   allocatable :: press  ! pressure 
real(8), dimension(:),   allocatable :: B      ! right hand side
real(8), dimension(:,:), allocatable :: A      ! FEM matrix
real(8), dimension(:),   allocatable :: work   ! work array needed by the solver
real(8), dimension(:),   allocatable :: bc_val ! array containing bc values
real(8), dimension(:),   allocatable :: csr_mat! csr storage of matrix nonzeroes
                                               !
real(8), external :: b1,b2,uth,vth,pth         ! body force and analytical solution
real(8) rq,sq,wq                               ! local coordinate and weight of qpoint
real(8) xq,yq                                  ! global coordinate of qpoint
real(8) uq,vq                                  ! velocity at qpoint
real(8) exxq,eyyq,exyq                         ! strain-rate components at qpoint  
real(8) Ael(m*ndof,m*ndof)                     ! elemental FEM matrix
real(8) Bel(m*ndof)                            ! elemental right hand side
real(8) N(m),dNdx(m),dNdy(m),dNdr(m),dNds(m)   ! shape fcts and derivatives
real(8) jcob                                   ! determinant of jacobian matrix
real(8) jcb(2,2)                               ! jacobian matrix
real(8) jcbi(2,2)                              ! inverse of jacobian matrix
real(8) Bmat(3,ndof*m)                         ! B matrix
real(8), dimension(3,3) :: Kmat                ! K matrix 
real(8), dimension(3,3) :: Cmat                ! C matrix
real(8) Aref                                   !
real(8) eps                                    !
real(8) rcond,fixt                                  !
                                               !
logical, dimension(:), allocatable :: bc_fix   ! prescribed b.c. array
                                               !
!==============================================!

logical :: use_petsc                           ! use petsc?
integer iError,iRow,iCol                       ! 
type(tMAT) :: PETSCMAAT                        ! PETSc matrix
type(tVEC) :: PETScRHS                         ! PETSc right hand side vector
type(tVEC) :: PETScSOL                         ! PETSc solution vector
type(tKSP) :: PETScCTXT                        ! PETSc KSP context
type(tPC)  :: PETScPREC                        ! PETSc PC context
PetscScalar, pointer :: LA_vec(:)

! No point putting this inside if(use_petsc) since the rest of the code does not support compilation without petsc,
! e.g. there are no guards placed around the petsc includes, petsc objects or petsc function calls
call PetscInitialize('petsc_default_options',iError)


!==============================================!
!=====[setup]==================================!
!==============================================!

Lx=1.d0
Ly=1.d0

nnx=29
nny=29
np=nnx*nny

nelx=nnx-1
nely=nny-1
nel=nelx*nely

penalty=1.d4

viscosity=1.d0
density=1.d0

Nfem=np*ndof

eps=1.d-10

Kmat(1,1)=1.d0 ; Kmat(1,2)=1.d0 ; Kmat(1,3)=0.d0  
Kmat(2,1)=1.d0 ; Kmat(2,2)=1.d0 ; Kmat(2,3)=0.d0  
Kmat(3,1)=0.d0 ; Kmat(3,2)=0.d0 ; Kmat(3,3)=0.d0  

Cmat(1,1)=2.d0 ; Cmat(1,2)=0.d0 ; Cmat(1,3)=0.d0  
Cmat(2,1)=0.d0 ; Cmat(2,2)=2.d0 ; Cmat(2,3)=0.d0  
Cmat(3,1)=0.d0 ; Cmat(3,2)=0.d0 ; Cmat(3,3)=1.d0  

NZ=((4*4+(2*(nnx-2)+2*(nny-2))*6+(nnx-2)*(nny-2)*9)) * (ndof**2) !; print *,nz

!==============================================!
!===[petsc initialization]=====================!
!==============================================!

use_petsc=.true.

if (use_petsc) then

   call VecCreate(PETSC_COMM_WORLD,PETScSOL,iError)
   call VecCreate(PETSC_COMM_WORLD,PETScRHS,iError)
   call VecSetSizes(PETScSOL,PETSC_DECIDE,Nfem,iError)
   call VecSetSizes(PETScRHS,PETSC_DECIDE,Nfem,iError)
   call VecSetFromOptions(PETScSOL,iError)
   call VecSetFromOptions(PETScRHS,iError)

end if

!==============================================!
!===[allocate memory]==========================!
!==============================================!

allocate(x(np))
allocate(y(np))
allocate(u(np))
allocate(v(np))
allocate(icon(m,nel))
allocate(A(Nfem,Nfem))
allocate(B(Nfem))
allocate(bc_fix(Nfem))
allocate(bc_val(Nfem))
allocate(press(nel))
allocate(csr_ia(Nfem+1))
allocate(csr_ja(NZ))
allocate(csr_ia_p(Nfem+1))
allocate(csr_ja_p(NZ))
allocate(nb_nz(Nfem))
allocate(csr_mat(NZ))

!==============================================!
!===[grid points setup]========================!
!==============================================!

counter=0
do j=0,nely
   do i=0,nelx
      counter=counter+1
      x(counter)=dble(i)*Lx/dble(nelx)
      y(counter)=dble(j)*Ly/dble(nely)
   end do
end do

!open(unit=123,file='OUT/gridnodes.dat',status='replace')
!write(123,'(a)') '#     xpos      ypos    node '
!do i=1,np
!   write(123,'(2f10.5,i8)') x(i),y(i),i
!end do
!close(123)

!==============================================!
!===[connectivity]=============================!
!==============================================!

counter=0
do j=1,nely
   do i=1,nelx
      counter=counter+1
      icon(1,counter)=i+(j-1)*(nelx+1)
      icon(2,counter)=i+1+(j-1)*(nelx+1)
      icon(3,counter)=i+1+j*(nelx+1)
      icon(4,counter)=i+j*(nelx+1)
   end do
end do

!open(unit=123,file='OUT/icon.dat',status='replace')
!do iel=1,nel
!   write(123,'(a)') '----------------------------'
!   write(123,'(a,i4,a)') '---element #',iel,' -----------'
!   write(123,'(a)') '----------------------------'
!   write(123,'(a,i8,a,2f20.10)') '  node 1 ', icon(1,iel),' at pos. ',x(icon(1,iel)),y(icon(1,iel))
!   write(123,'(a,i8,a,2f20.10)') '  node 2 ', icon(2,iel),' at pos. ',x(icon(2,iel)),y(icon(2,iel))
!   write(123,'(a,i8,a,2f20.10)') '  node 3 ', icon(3,iel),' at pos. ',x(icon(3,iel)),y(icon(3,iel))
!   write(123,'(a,i8,a,2f20.10)') '  node 4 ', icon(4,iel),' at pos. ',x(icon(4,iel)),y(icon(4,iel))
!end do
!close(123)

!==============================================!
!=====[fill ia,ja]=============================!
!==============================================!

NZ=0
csr_ia(1)=1
do j1=1,nny
do i1=1,nnx
   ip=(j1-1)*nnx+i1 ! node number
   do k=1,ndof
      ii=2*(ip-1) + k ! address in the matrix
      nsees=0
      do j2=-1,1 ! exploring neighbouring nodes
      do i2=-1,1
         i=i1+i2
         j=j1+j2
         if (i>=1 .and. i<= nnx .and. j>=1 .and. j<=nny) then ! if node exists
            jp=(j-1)*nnx+i  ! node number of neighbour 
            do l=1,ndof
               jj=2*(jp-1)+l  ! address in the matrix
               nz=nz+1
               csr_ja(nz)=jj
               nsees=nsees+1
            end do
         end if
      end do
      end do
      csr_ia(ii+1)=csr_ia(ii)+nsees
   end do ! loop over ndofs
end do
end do

!print *,nz
!print *,minval(csr_ia), maxval(csr_ia)
!print *,minval(csr_ja), maxval(csr_ja)

!==============================================!
!=====[fill nb_nz]=============================!
!==============================================!

nb_nz=0
do i=1,Nfem
   nb_nz(i)=(csr_ia(i+1)-csr_ia(i))
end do
if (sum(nb_nz)/=NZ) stop 'pb with nb_nz'
!print *,minval(nb_nz), maxval(nb_nz)


if (use_petsc) then
  csr_mat=0.d0
  csr_ia_p = csr_ia
  csr_ja_p = csr_ja

  csr_ia_p = csr_ia_p - 1
  csr_ja_p = csr_ja_p - 1

  ! Create petsc matrix using csr and empty values
  ! Creating it here allows you to re-use it
  ! When using this function, fortran arrays are shared with the petsc mat.
  ! Since petsc requires indices starting with 0, but your assembly uses indices starting with 1, it
  ! it simpler to simply copy the ia, ja arrays which petsc will use
  call MatCreateSeqAIJWithArrays(PETSC_COMM_WORLD,Nfem,Nfem,csr_ia_p,csr_ja_p,csr_mat,PETScMAAT,iError)

end if

!==============================================!
!=====[define bc]==============================!
!==============================================!

bc_fix=.false.

do i=1,np
   if (x(i).lt.eps) then
      bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
      bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
   endif
   if (x(i).gt.(Lx-eps)) then
      bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
      bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
   endif
   if (y(i).lt.eps) then
      bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
      bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
   endif
   if (y(i).gt.(Ly-eps) ) then
      bc_fix((i-1)*ndof+1)=.true. ; bc_val((i-1)*ndof+1)=0.d0
      bc_fix((i-1)*ndof+2)=.true. ; bc_val((i-1)*ndof+2)=0.d0
   endif
end do

!open(unit=123,file='OUT/bc_u.dat',status='replace')
!open(unit=234,file='OUT/bc_v.dat',status='replace')
!do i=1,np
!   if (bc_fix((i-1)*ndof+1)) write(123,'(3f20.10)') x(i),y(i),u(i)
!   if (bc_fix((i-1)*ndof+2)) write(234,'(3f20.10)') x(i),y(i),v(i)
!end do
!close(123)
!close(234)

!==============================================!
!=====[build FE matrix]========================!
!==============================================!

A=0.d0
B=0.d0
csr_mat=0.d0

do iel=1,nel

   Ael=0.d0
   Bel=0.d0

   do iq=-1,1,2
   do jq=-1,1,2

      rq=iq/sqrt(3.d0)
      sq=jq/sqrt(3.d0)
      wq=1.d0*1.d0

      N(1)=0.25d0*(1.d0-rq)*(1.d0-sq)
      N(2)=0.25d0*(1.d0+rq)*(1.d0-sq)
      N(3)=0.25d0*(1.d0+rq)*(1.d0+sq)
      N(4)=0.25d0*(1.d0-rq)*(1.d0+sq)

      dNdr(1)= - 0.25d0*(1.d0-sq)   ;   dNds(1)= - 0.25d0*(1.d0-rq)
      dNdr(2)= + 0.25d0*(1.d0-sq)   ;   dNds(2)= - 0.25d0*(1.d0+rq)
      dNdr(3)= + 0.25d0*(1.d0+sq)   ;   dNds(3)= + 0.25d0*(1.d0+rq)
      dNdr(4)= - 0.25d0*(1.d0+sq)   ;   dNds(4)= + 0.25d0*(1.d0-rq)

      jcb=0.d0
      do k=1,m
         jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
         jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
         jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
         jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
      enddo

      jcob=jcb(1,1)*jcb(2,2)-jcb(2,1)*jcb(1,2)

      jcbi(1,1)=    jcb(2,2) /jcob
      jcbi(1,2)=  - jcb(1,2) /jcob
      jcbi(2,1)=  - jcb(2,1) /jcob
      jcbi(2,2)=    jcb(1,1) /jcob

      xq=0.d0
      yq=0.d0
      uq=0.d0
      vq=0.d0
      exxq=0.d0
      eyyq=0.d0
      exyq=0.d0
      do k=1,m
         xq=xq+N(k)*x(icon(k,iel))
         yq=yq+N(k)*y(icon(k,iel))
         uq=uq+N(k)*u(icon(k,iel))
         vq=vq+N(k)*v(icon(k,iel))
         dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)
         exxq=exxq+ dNdx(k)*u(icon(k,iel))
         eyyq=eyyq+ dNdy(k)*v(icon(k,iel))
         exyq=exyq+ dNdx(k)*v(icon(k,iel)) *0.5d0 &
                  + dNdy(k)*u(icon(k,iel)) *0.5d0
      end do

      do i=1,m
         i1=2*i-1
         i2=2*i
         Bmat(1,i1)=dNdx(i) ; Bmat(1,i2)=0.d0
         Bmat(2,i1)=0.d0    ; Bmat(2,i2)=dNdy(i)
         Bmat(3,i1)=dNdy(i) ; Bmat(3,i2)=dNdx(i)
      end do

      Ael=Ael + matmul(transpose(Bmat),matmul(viscosity*Cmat,Bmat))*wq*jcob

      do i=1,m
      i1=2*i-1
      i2=2*i
      !Bel(i1)=Bel(i1)+N(i)*jcob*wq*density*gx
      !Bel(i2)=Bel(i2)+N(i)*jcob*wq*density*gy
      Bel(i1)=Bel(i1)+N(i)*jcob*wq*b1(xq,yq)
      Bel(i2)=Bel(i2)+N(i)*jcob*wq*b2(xq,yq)
      end do

   end do
   end do

   ! 1 point integration

      rq=0.d0
      sq=0.d0
      wq=2.d0*2.d0

      N(1)=0.25d0*(1.d0-rq)*(1.d0-sq)
      N(2)=0.25d0*(1.d0+rq)*(1.d0-sq)
      N(3)=0.25d0*(1.d0+rq)*(1.d0+sq)
      N(4)=0.25d0*(1.d0-rq)*(1.d0+sq)

      dNdr(1)= - 0.25d0*(1.d0-sq)   ;   dNds(1)= - 0.25d0*(1.d0-rq)
      dNdr(2)= + 0.25d0*(1.d0-sq)   ;   dNds(2)= - 0.25d0*(1.d0+rq)
      dNdr(3)= + 0.25d0*(1.d0+sq)   ;   dNds(3)= + 0.25d0*(1.d0+rq)
      dNdr(4)= - 0.25d0*(1.d0+sq)   ;   dNds(4)= + 0.25d0*(1.d0-rq)

      jcb=0.d0
      do k=1,m
         jcb(1,1)=jcb(1,1)+dNdr(k)*x(icon(k,iel))
         jcb(1,2)=jcb(1,2)+dNdr(k)*y(icon(k,iel))
         jcb(2,1)=jcb(2,1)+dNds(k)*x(icon(k,iel))
         jcb(2,2)=jcb(2,2)+dNds(k)*y(icon(k,iel))
      enddo

      jcob=jcb(1,1)*jcb(2,2)-jcb(2,1)*jcb(1,2)

      jcbi(1,1)=    jcb(2,2) /jcob
      jcbi(1,2)=  - jcb(1,2) /jcob
      jcbi(2,1)=  - jcb(2,1) /jcob
      jcbi(2,2)=    jcb(1,1) /jcob

      do k=1,m
         dNdx(k)=jcbi(1,1)*dNdr(k)+jcbi(1,2)*dNds(k)
         dNdy(k)=jcbi(2,1)*dNdr(k)+jcbi(2,2)*dNds(k)
      end do

      do i=1,m
         i1=2*i-1
         i2=2*i
         Bmat(1,i1)=dNdx(i) ; Bmat(1,i2)=0.d0
         Bmat(2,i1)=0.d0    ; Bmat(2,i2)=dNdy(i)
         Bmat(3,i1)=dNdy(i) ; Bmat(3,i2)=dNdx(i)
      end do

      Ael=Ael + matmul(transpose(Bmat),matmul(penalty*Kmat,Bmat))*wq*jcob

      !======================================
      !=====[impose boundary conditions]=====
      !======================================
    
      do k1=1,m  
         ik=icon(k1,iel)  
         do i1=1,ndof    
            m1=(ik-1)*ndof+i1    
            if (bc_fix(m1)) then  
            fixt=bc_val(m1) 
            i=(k1-1)*ndof+i1    
            Aref=Ael(i,i)    
            do j=1,m*ndof    
               Bel(j)=Bel(j)-Ael(j,i)*fixt 
               Ael(i,j)=0.d0    
               Ael(j,i)=0.d0    
            enddo    
            Ael(i,i)=Aref    
            Bel(i)=Aref*fixt    
            endif    
         enddo    
      enddo    

      !=====================
      !=====[assemble]======
      !=====================

      do k1=1,m
         ik=icon(k1,iel)
         do i1=1,ndof
            ikk=ndof*(k1-1)+i1
            m1=ndof*(ik-1)+i1
            do k2=1,m
               jk=icon(k2,iel)
               do i2=1,ndof
                  jkk=ndof*(k2-1)+i2
                  m2=ndof*(jk-1)+i2
                  do k=csr_ia(m1),csr_ia(m1+1)-1    
                     if (csr_ja(k)==m2) then  
                        csr_mat(k)=csr_mat(k)+Ael(ikk,jkk)  
                     end if    
                  end do    
                  A(m1,m2)=A(m1,m2)+Ael(ikk,jkk)
               end do
            end do
            B(m1)=B(m1)+Bel(ikk)
         end do
      end do

end do

!==============================================!

if (use_petsc) then

  call VecGetArrayF90(PETScRHS,LA_vec,iError) ; call check_error(iError)
  do iRow=1,Nfem
    LA_vec(iRow) = B(iRow)
  end do
  call VecRestoreArrayF90(PETScRHS,LA_vec,iError) ; call check_error(iError)
  ! Call VecAssemblyBegin(),VecAssemblyEnd() to facilitate sending data to other rank
  call VecAssemblyBegin(PETScRHS,iError) ; call check_error(iError)
  call VecAssemblyEnd(PETScRHS,iError) ; call check_error(iError)

   !call MatCreateSeqAIJ(PETSC_COMM_WORLD,Nfem,Nfem,PETSC_DECIDE,PETSC_NULL_INTEGER,PETScMAAT,iError)
   !do iRow=1,Nfem
   !do iCol=1,Nfem
   !if (C(iRow,iCol)) then
   !   call MatSetValues(PETScMAAT,1,iRow-1,1,iCol-1,A(iRow,iCol),INSERT_VALUES,iError)
   !endif
   !enddo
   !enddo

   call MatAssemblyBegin(PETScMAAT,MAT_FINAL_ASSEMBLY,iError)
   call MatAssemblyEnd(PETScMAAT,MAT_FINAL_ASSEMBLY,iError)

end if

!==============================================!
!=====[solve system]===========================!
!==============================================!

if (use_petsc) then

   call KSPCreate(PETSC_COMM_WORLD,PETScCTXT,iError)
   call KSPSetOperators(PETScCTXT,PETScMAAT,PETScMAAT,iError) ; call check_error(iError)
   call KSPSetFromOptions(PETScCTXT,iError)
   call KSPSolve(PETScCTXT,PETScRHS,PETScSOL,iError) ; call check_error(iError)

  call VecGetArrayReadF90(PETScSOL,LA_vec,iError) ; call check_error(iError)
  do i=1,np
    j = (i-1)*ndof + 1
    u(i) = LA_vec(j)
    j = (i-1)*ndof + 2
    v(i) = LA_vec(j)
  end do
  call VecRestoreArrayReadF90(PETScSOL,LA_vec,iError) ; call check_error(iError)


   print *,'v_x (m/M):',minval(u),maxval(u)
   print *,'v_y (m/M):',minval(v),maxval(v)

else ! use linpack subroutines

   job=0
   allocate(work(Nfem))
   allocate(ipvt(Nfem))
   call DGECO (A, Nfem, Nfem, ipvt, rcond, work)
   call DGESL (A, Nfem, Nfem, ipvt, B, job)
   deallocate(ipvt)
   deallocate(work)

   do i=1,np
   u(i)=B((i-1)*ndof+1)
   v(i)=B((i-1)*ndof+2)
   end do

end if

open(unit=123,file='OUT/solution_u.dat',status='replace')
open(unit=234,file='OUT/solution_v.dat',status='replace')
do i=1,np
   write(123,'(5f20.10)') x(i),y(i),u(i),uth(x(i),y(i)),u(i)-uth(x(i),y(i))
   write(234,'(5f20.10)') x(i),y(i),v(i),vth(x(i),y(i)),v(i)-vth(x(i),y(i))
end do
close(123)
close(234)

!==============================================!
!=====[destroying petsc]=======================!
!==============================================!

if (use_petsc) then
  call MatDestroy(PETScMAAT,iError)
  call VecDestroy(PETScRHS, iError)
  call VecDestroy(PETScSOL, iError)
  call KSPDestroy(PETScCTXT, iError)
end if

call PetscFinalize(iError)

end program

!==============================================================================

function uth (x,y)
real(8) uth,x,y
uth = x**2 * (1.d0-x)**2 * (2.d0*y - 6.d0*y**2 + 4*y**3)
end function

function vth (x,y)
real(8) vth,x,y
vth = -y**2 * (1.d0-y)**2 * (2.d0*x - 6.d0*x**2 + 4*x**3)
end function

function pth (x,y)
real(8) pth,x,y
pth = x*(1.d0-x)
end function

function b1 (x,y)
real(8) b1,x,y
b1 = ( (12.d0-24.d0*y)*x**4 + (-24.d0+48.d0*y)*x**3 + (-48.d0*y+72.d0*y**2-48.d0*y**3+12.d0)*x**2 &
   + (-2.d0+24.d0*y-72.d0*y**2+48.d0*y**3)*x + 1.d0-4.d0*y+12.d0*y**2-8.d0*y**3 )
end function

function b2 (x,y)
real(8) b2,x,y
b2= ( (8.d0-48.d0*y+48.d0*y**2)*x**3 + (-12.d0+72.d0*y-72*y**2)*x**2 + &
    (4.d0-24.d0*y+48.d0*y**2-48.d0*y**3+24.d0*y**4)*x - 12.d0*y**2 + 24.d0*y**3 -12.d0*y**4)
end function

!==============================================================================

subroutine check_error (ierror)
implicit none
integer ierror
if (ierror/=0) stop 'iError/=0'
end subroutine 




