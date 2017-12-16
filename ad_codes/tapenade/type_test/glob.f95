module glob
      type soln
              real::x
              integer::y,z
      end type soln
      real::func
contains

        subroutine solver(func,s)
                implicit none
                type(soln),intent(in)::s
                real,intent(out)::func
                func=s%x+s%y+s%z
        end subroutine

end module glob
