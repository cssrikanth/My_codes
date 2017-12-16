        subroutine prime(x1,x2,f)
                double precision,intent(in)::x1,x2
                double precision,intent(out)::f
                f=x1**3+x2**4+2*x1*x2
        end subroutine prime
