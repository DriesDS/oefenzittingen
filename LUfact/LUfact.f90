program lufact

implicit none

integer,parameter :: seed = 450
call lu()

contains

	subroutine lu()
		real, dimension(3,3) :: A,B,L,U,C
		integer :: INFO, LDA, M, N, i
		integer :: IPIV(3), WORK(100)

		M = 3
		N = 3
		LDA = 3
		call random_seed()
		call RANDOM_NUMBER(A)
		B = A

		call sgetrf(M,N,A,LDA,IPIV,INFO)

		L = 0
		U = 0
		do i = 1,3
			L(i,:i) = A(i,:i)
			U(i,i:) = A(i,i:)
			L(i,i) = 1
		enddo

		write(*,'(3(3e12.3,/))') B

		C = matmul(L,U)

		call slaswp(N,C,N,1,N,IPIV,1)

		write(*,'(3(3e12.3,/))') C

		call sgetri(N,A,LDA,IPIV,WORK,100,INFO)

		write(*,'(3(3e12.3,/))') matmul(A,B)		

	end subroutine

end program
