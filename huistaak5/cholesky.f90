program cholesky

implicit none

call ginibre()

contains

	subroutine ginibre()
		integer, parameter :: M = 7
		real, dimension(M,M) :: A, B
		complex, dimension(M,M) :: C, Cher, L, Lprod
		complex :: alpha, beta, test
		integer :: i, info

		call random_seed()
		call random_number(A)
		call random_number(B)
		call gaussian(A,B,C)

		alpha = cmplx(1.0, 0.0)
		beta = cmplx(0.0, 0.0)
		call cgemm('N', 'C', M, M, M, alpha, C, M, C, M, beta, Cher, M)
		C = Cher

		call cpotrf('L',M,Cher,M,info)

		L = 0
		do i = 1,M
			L(i:,i) = Cher(i:,i)
		enddo

		call cgemm('N', 'C', M, M, M, alpha, L, M, L, M, beta, Lprod, M)

		write(*,'(a)') '! equality within a range of 1.e-5: equality within a range of 1.e-6: '
		do i = 1,M
			write(*,'(a,t15,7(L2),t48,7(L2))') '! ', abs(C(i,:)-Lprod(i,:))<1.e-5, abs(C(i,:)-Lprod(i,:))<1.e-6
		enddo
	end subroutine

	ELEMENTAL subroutine gaussian(a, b, c)
		real, intent(in) :: a, b
		complex, intent(out) :: c
		complex :: minus1

		minus1 = cmplx(-1.0, 0.0)
		c = sqrt(-2*log(a))*exp(8*atan(1.0)*sqrt(minus1)*b)

	end subroutine

end program
