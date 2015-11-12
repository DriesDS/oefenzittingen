program cholesky

implicit none

call ginibre()

contains

	subroutine ginibre()
		real, dimension(7,7) :: A, B
		complex, dimension(7,7) :: C, C1, C2

		call random_seed(2568)
		call random_number(A)
		call random_number(B)
		call gaussian(A,B,C1)

		write(*,'(7(e12.3,/)') A
		write(*,'(7(e12.3,/)') B
		write(*,'(7(e12.3,/)') C1

		C2 = transpose(C1)
		C2 = conjg(C2)

		C = matmul(C1,C2)

		write(*,'(7(e12.3,/)') C

	end subroutine

	ELEMENTAL subroutine gaussian(a, b, c)
		intent(in) :: a, b
		intent(out) :: c

		c = sqrt(-2*log(a))*exp(8*atan(1.0)*sqrt(-1)*b)

	end subroutine

end program
