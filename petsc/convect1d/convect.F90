program main
  implicit none

#include <petsc/finclude/petscsys.h>
#include <petsc/finclude/petscvec.h>
#include <petsc/finclude/petscdmda.h>
#include <petsc/finclude/petscvec.h90>
#include <petsc/finclude/petscdmda.h90>
#include <petsc/finclude/petscdm.h>

  PetscErrorCode    ierr
  DM  da
  Vec  ug,ul
  PetscInt  i,ibeg,nloc,nx
  PetscInt sw,ndof,il,nl
  PetscMPIInt rank,size
  PetscReal dx,x
  PetscScalar v
  PetscScalar, pointer :: u(:),unew(:)
  MPI_Comm comm
  real :: cfl,xmin,xmax
  real,external :: initcond,weno5
  real,allocatable,dimension(:)::res,uold
  real::dt,lam,tfinal,t,ark(3),uleft,flux
  integer::rk,j,jm1,jm2,jm3,jp1
  
  
  xmin = -1.0
  xmax = +1.0
  nx = 10
  ark = [0.0d0, 3.0d0/4.0d0, 1.0d0/3.0d0]
  sw = 3.0
  ndof = 1.0

  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  if (ierr .ne. 0) then
    print*,'Unable to initialize PETSc'
    stop
  endif
  comm = PETSC_COMM_WORLD
  call MPI_Comm_rank(comm,rank,ierr)
  call MPI_Comm_size(comm,size,ierr)
  if (rank == 0) then
    write(*,*) 'We are solving 1D convection equation using ',size,' processes.'
  endif
  call DMDACreate1d(comm, DM_BOUNDARY_PERIODIC, nx, ndof, sw,&
  &   PETSC_NULL_INTEGER, da, ierr)

  call DMSetFromOptions(da,ierr)
 
  call DMSetup(da,ierr) 

  call DMCreateGlobalVector(da,ug,ierr)  

  call DMDAGetCorners(da,ibeg,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,&
  & nloc,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)

  dx = (xmin - xmax) / (nx)

  do i=ibeg,ibeg+nloc-1
    x = xmin + i*dx
    v = initcond(x)
    call VecSetValues(ug,ndof,i,v,INSERT_VALUES,ierr)
  end do

  call VecAssemblyBegin(ug,ierr)
  call VecAssemblyEnd(ug,ierr)

  call DMGetLocalVector(da,ul,ierr)

!  call DMDAGetGhostCorners(da,il,PETSC_NULL_INTEGER,&
!  &PETSC_NULL_INTEGER,nl,PETSC_NULL_INTEGER,PETSC_NULL_INTEGER,ierr)

!  allocate (res(nloc))
!  allocate (uold(nloc))  

!  dt = cfl * dx
!  lam = dt/dx
!  tfinal = 2.0
!  t = 0.0
!  do while(t < tfinal)

!    if(t+dt > tfinal) then
!      dt = tfinal - t
!      lam = dt/dx
!    end if
!      do rk=1,3
!        call DMGlobalToLocalBegin(da,ug,INSERT_VALUES,ul,ierr)
!        call DMGlobalToLocalEnd(da,ug,INSERT_VALUES,ul,ierr)

!        call VecGetArrayRead(da,ul,u,ierr)
!        call VecGetArray(da,ug,unew,ierr)

!        if (rk == 0) then
!          do i=ibeg,ibeg+nloc-1
!            uold(i-ibeg) = u(i)
!          end do
!        end if

!        do i=0,nloc-1
!          res(i) = 0.0
!        end do

!        do i=0,nloc
!          j = il+sw+i
!          jm1 = j-1
!          jm2 = j-2
!          jm3 = j-3
!          jp1 = j+1
   
!          uleft = weno5(u(jm3),u(jm2),u(jm1),u(j),u(jp1))
!          flux = uleft

!          if (i==0) then
!            res(i) = res(i)-flux
!          elseif(i==nloc) then
!            res(i-1) = res(i-1)+flux
!          else
!            res(i) = res(i)-flux
!            res(i-1) = res(i-1)+flux
!          end if  

!        end do 

!        do i=ibeg,ibeg+nloc-1
!          unew(i) = ark(rk)*uold(i-ibeg) + (1-ark(rk))*(u(i)-&
!          & lam * res(i-ibeg)) 
!        end do

  
!        call VecRestoreArrayRead(da,ul,u,ierr)
!        call VecRestoreArray(da,ug,unew,ierr)            


!      end do

!      t = t+dt
!  end do

  call DMDestroy(da,ierr)
  call VecDestroy(ug,ierr)
  call PetscFinalize(ierr)

end program main


real function initcond(x)
  implicit none
  
  real :: x

  if( x >= -0.8 .and. x <= -0.6) then 
    initcond = exp(-log(2.0)*((x+0.7)**2/0.0009))
  elseif(x >= -0.4 .and. x <= -0.2) then
    initcond = 1.0
  elseif(x >= 0.0 .and. x <= 0.2) then
    initcond = 1.0 - abs(10*(x-0.1))
  elseif( x>= 0.4 .and. x <= 0.6) then
    initcond = sqrt(1-100*(x-0.5)**2)
  else
    initcond = 0.0
  end if
end function initcond


real function weno5(um2,um1,u0,up1,up2)
  implicit none
 
  integer :: um2,um1,up1,up2,u0
  real :: u1, u2,u3
  real :: w1,w2,w3
  real :: eps,gamma1,gamma2,gamma3
  real :: beta1,beta2,beta3

  eps = 1.0e-6
  gamma1 = 1.0/10.0
  gamma2 = 3.0/5.0
  gamma3 = 3.0/10.0

  beta1 = (13.0/12.0)*(um2 - 2.0*um1 + u0)**2 + &
         &  (1.0/4.0)*(um2 - 4.0*um1 + 3.0*u0)**2
  beta2 = (13.0/12.0)*(um1 - 2.0*u0 + up1)**2 + &
         &  (1.0/4.0)*(um1 - up1)**2
  beta3 = (13.0/12.0)*(u0 - 2.0*up1 + up2)**2 + &
         &  (1.0/4.0)*(3.0*u0 - 4.0*up1 + up2)**2

  w1 = gamma1 / (eps+beta1)**2
  w2 = gamma2 / (eps+beta2)**2
  w3 = gamma3 / (eps+beta3)**2

  u1 = (1.0/3.0)*um2 - (7.0/6.0)*um1 + (11.0/6.0)*u0
  u2 = -(1.0/6.0)*um1 + (5.0/6.0)*u0 + (1.0/3.0)*up1
  u3 = (1.0/3.0)*u0 + (5.0/6.0)*up1 - (1.0/6.0)*up2

  
  weno5= (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3)

end function weno5


