program convergence

use FMZM

implicit none

type(FM) :: x
x = TO_FM('0.9878')
call FM_SET(50)

call calc_zero(x)

contains

	subroutine calc_zero(x)
		type(FM) :: x, x12, cosx12, afg
		integer :: i

		WRITE (*,"(i5,A)") i,TRIM(FM_FORMAT('E35.17',x))
		do i = 1,10
			x12 = x**12
			cosx12 = cos(x12)
			afg = -12*sin(x12)*x**11
			x = x-cosx12/afg
			WRITE (*,"(i5,A)") i,TRIM(FM_FORMAT('E60.50',x))
		enddo
	end subroutine

end program
