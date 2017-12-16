program main
      implicit none
      integer::n,ncount,i,j
      real::dx,error,tol,lent,pertub
!      real,allocatable,intent(in),dimension(:,:)::k
      real::func,funcd
      real,allocatable,dimension(:,:)::t,x,y,oldt,temp,k,kd,td

      !grid details
      n=100
      allocate(t(n,n))
      allocate(oldt(n,n))
      allocate(temp(n,n))
      allocate(k(n,n))
      allocate(kd(n,n))
      allocate(td(n,n))

      !solution details
      error=0.0
      tol=1e-8
      ncount=0
      pertub=0

      !initial guess and bc

     do i=1,n
        do j=1,n
                if (i == 0) then
                   t(i,j)=10
                   k(i,j)=10+pertub
                   kd(i,j)=1
                elseif (j==0)then
                   t(i,j)=20
                   kd(i,j)=1
                   k(i,j)=5+pertub
                elseif (i==n) then
                   t(i,j)=20
                   kd(i,j)=1
                   k(i,j)=1+pertub
                elseif (j==n) then
                   t(i,j)=40
                   kd(i,j)=1
                   k(i,j)=50+pertub
                else
                   t(i,j)=0.0
                   kd(i,j)=1
                   k(i,j)=60+pertub
                end if
        end do
      end do
      

     !solver

     do
     oldt=t
!     do i=2,n-1
!        do j=2,n-1
!               T(i,j)=(T(i+1,j)*k(i+1,j)+T(i-1,j)*k(i-1,j)+T(i,j+1)*&
!                       k(i,j+1)+T(i,j-1)*k(i,j-1))&
!                       /(k(i+1,j)+k(i-1,j)+k(i,j+1)+k(i,j-1))
!        end do
!      end do

     !objective func
!     func=sum(T)

!     call solver(t,k,func,n)
     call solver_d(t,td,k,kd,func,funcd,n)
     do i=2,n-1
         do j=2,n-1
               temp(i,j)=abs(t(i,j)-oldt(i,j))
         end do
     end do
     
     error=sum(temp)/sum(abs(t))

     

     print *, error
      
     if (error  <tol  ) then
             print*,"func:",func,"funcd:",funcd
             print *, "solution converged with residual:",error,"after",ncount,"iterations"

             exit

     end if
     ncount=ncount+1
     end do


end program

