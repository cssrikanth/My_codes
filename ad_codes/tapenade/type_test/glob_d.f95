!        Generated by TAPENADE     (INRIA, Tropics team)
!  Tapenade 3.9 (r5096) - 24 Feb 2014 16:53
!
MODULE GLOB_D
  IMPLICIT NONE
  TYPE SOLN
      REAL :: x
      INTEGER :: y, z
  END TYPE SOLN
  TYPE SOLN_D
      REAL :: x
  END TYPE SOLN_D
  REAL :: func

CONTAINS
!  Differentiation of solver in forward (tangent) mode (with options i4 dr8 r4):
!   variations   of useful results: func
!   with respect to varying inputs: s.x
!   RW status of diff variables: s.x:in func:out
  SUBROUTINE SOLVER_D(func, funcd, s, sd)
    IMPLICIT NONE
    TYPE(SOLN), INTENT(IN) :: s
    TYPE(SOLN_D), INTENT(IN) :: sd
    REAL, INTENT(OUT) :: func
    REAL, INTENT(OUT) :: funcd
    funcd = sd%x
    func = s%x + s%y + s%z
  END SUBROUTINE SOLVER_D
  SUBROUTINE SOLVER(func, s)
    IMPLICIT NONE
    TYPE(SOLN), INTENT(IN) :: s
    REAL, INTENT(OUT) :: func
    func = s%x + s%y + s%z
  END SUBROUTINE SOLVER
END MODULE GLOB_D
