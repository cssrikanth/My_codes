        program main
                implicit none
                double precision::x1,x2,f,fb,x1b,x2b,touch,x1d,x2d,fd
		touch =2
		if (touch ==1) then
		print *,"adjoint used"
                read(*,*)x1,x2,fb
		call prime_B(x1,x1b,x2,x2b,f,fb)
		print *,x1b,x2b
		else
		print *,"tangent used"
		read(*,*)x1,x2,x1d,x2d
		call PRIME_D(x1, x1d, x2, x2d, f, fd)
		print*,fd
		end if
                
        end program main



  
