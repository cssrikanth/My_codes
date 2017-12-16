program main
  
use glob
implicit none

type(soln)::s
type(soln_d)::sd

real::funcd

read(*,*)s%x,s%y,s%z

call SOLVER_D(func, funcd, s, sd)

call solver (func, s)

        print *,func,funcd

end program main








