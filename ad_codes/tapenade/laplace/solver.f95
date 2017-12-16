subroutine solver(T,k,func,n)
      implicit none
      integer::i,j,n
      real,intent(out)::func
      real,dimension(n,n)::t
      real,intent(in),dimension(n,n)::k
           do i=2,n-1
             do j=2,n-1
                 T(i,j)=(T(i+1,j)*k(i+1,j)+T(i-1,j)*k(i-1,j)+T(i,j+1)*k(i,j+1)+T(i,j-1)*k(i,j-1)) &
                 &      /(k(i+1,j)+k(i-1,j)+k(i,j+1)+k(i,j-1))
             end do
           end do 
           func=sum(T)
end subroutine


