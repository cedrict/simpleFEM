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
integer, dimension(:,:), allocatable :: icon   ! connectivity array
integer, dimension(:), allocatable :: ipvt     ! work array needed by the solver 
                                               !
integer i1,i2,i,j,k,iel,counter,iq,jq          !
integer ik,jk,ikk,jkk,m1,m2,k1,k2,job,jfake(1) !
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
real(8) rcond                                  !
                                               !
logical, dimension(:), allocatable :: bc_fix   ! prescribed b.c. array
logical, dimension(:,:), allocatable :: C      ! non-zero terms in FEM matrix
                                               !
!==============================================!

logical :: use_petsc                           ! use petsc?
integer iError,iRow,iCol                       ! 
type(tMAT) :: PETSCMAAT                        ! PETSc matrix
type(tVEC) :: PETScRHS                         ! PETSc right hand side vector
type(tVEC) :: PETScSOL                         ! PETSc solution vector
type(tKSP) :: PETScCTXT                        ! PETSc KSP context
type(tPC)  :: PETScPREC                        ! PETSc PC context

!==============================================!
!=====[setup]==================================!
!==============================================!

Lx=1.d0
Ly=1.d0

nnx=21
nny=21
np=nnx*nny

nelx=nnx-1
nely=nny-1
nel=nelx*nely

penalty=1.d7

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

!==============================================!
!===[petsc initialization]=====================!
!==============================================!

use_petsc=.true.

if (use_petsc) then

   call PetscInitialize('petsc_default_options',iError)
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
allocate(C(Nfem,Nfem))
allocate(bc_fix(Nfem))
allocate(bc_val(Nfem))
allocate(press(nel))

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
C=.false.

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
                  A(m1,m2)=A(m1,m2)+Ael(ikk,jkk)
                  C(m1,m2)=.true.
               end do
            end do
            B(m1)=B(m1)+Bel(ikk)
         end do
      end do

end do

!==============================================!
!=====[impose b.c.]============================!
!==============================================!

do i=1,Nfem
    if (bc_fix(i)) then
      Aref=A(i,i)
      do j=1,Nfem
         B(j)=B(j)-A(i,j)*bc_val(i)
         A(i,j)=0.d0
         A(j,i)=0.d0
      enddo
      A(i,i)=Aref
      B(i)=Aref*bc_val(i)
   endif
enddo

!open(unit=123,file='OUT/matrix.dat',status='replace')
!open(unit=234,file='OUT/rhs.dat',status='replace')
!do i=1,Nfem
!   do j=1,Nfem
!      if (C(i,j)) write(123,'(2i6,f20.10)') i,Nfem-j,A(i,j)
!   end do
!   write(234,'(i6,f20.10)') i,B(i)
!end do
!close(123)
!close(234)

!==============================================!

if (use_petsc) then

   do iRow=0,Nfem-1
   call VecSetValue(PETScRHS,iRow,B(iRow+1),INSERT_VALUES,iError) 
   end do
   call VecAssemblyBegin(PETScRHS,iError)
   call VecAssemblyEnd(PETScRHS,iError)

   !******************************************************************
   ! this will be replaced by MatCreateSeqAIJWithArrays as soon as 
   ! we put CSR storage in this code :)
   !******************************************************************

   call MatCreateSeqAIJ(PETSC_COMM_WORLD,Nfem,Nfem,PETSC_DECIDE,PETSC_NULL_INTEGER,PETScMAAT,iError)
   do iRow=1,Nfem
   do iCol=1,Nfem
   if (C(iRow,iCol)) then
      call MatSetValues(PETScMAAT,1,iRow-1,1,iCol-1,A(iRow,iCol),INSERT_VALUES,iError)
   endif
   enddo
   enddo
   call MatAssemblyBegin(PETScMAAT,MAT_FINAL_ASSEMBLY,iError)
   call MatAssemblyEnd(PETScMAAT,MAT_FINAL_ASSEMBLY,iError)

end if

!==============================================!
!=====[solve system]===========================!
!==============================================!

if (use_petsc) then

   call KSPCreate(PETSC_COMM_WORLD,PETScCTXT,iError)
   call KSPSetOperators(PETScCTXT,PETScMAAT,PETScMAAT,iError) ; call check_error(iError)
   call KSPGetPC(PETScCTXT, PETScPREC, iError)
   call PCSetFromOptions(PETScPREC, iError)
   call KSPSetFromOptions(PETScCTXT,iError)
   call KSPSolve(PETScCTXT,PETScRHS,PETScSOL,iError) ; call check_error(iError)

   do i=1,np
      j=(i-1)*ndof+1
      jfake=j-1
      call VecGetValues(PETScSOL,1,jfake,u(i),iError)
      j=(i-1)*ndof+2
      jfake=j-1
      call VecGetValues(PETScSOL,1,jfake,v(i),iError)
   end do

else ! use linpack subroutines

   job=0
   allocate(work(Nfem))
   allocate(ipvt(Nfem))
   !call DGECO (A, Nfem, Nfem, ipvt, rcond, work)
   !call DGESL (A, Nfem, Nfem, ipvt, B, job)
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
call PetscFinalize(iError)

end if
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




