! DMR
!   Tests of several square matrix-matrix products
! 
! HISTORY
!   20151125 PO - Shortened
!   20110331 KP - Added documentation, better structure.
!   2008---- BV - Initial version
! 
! AUTHORS
!   Peter Opsomer, KU Leuven CS, Belgium
!   Koen Poppe, Nikon Metrology, Brussels, Belgium
!   Bart Vandewoestyne, Sizing Servers Lab, Bruges, Belgium
!

program timingtest
    use matrixop
    implicit none
    !--------------------------------------------------------------------------
    ! Abstract interfaces
    !
    ! NOTE: this simplifies the timings.
    !--------------------------------------------------------------------------
    abstract interface
        subroutine a_maal_b_interface(a, b, c)
            import dp
            real(kind=dp), dimension(:,:), intent(in)  :: a, b
            real(kind=dp), dimension(:,:), intent(out) :: c
        end subroutine a_maal_b_interface
        subroutine a_maal_b_blocks_interface(a,b,c,blocksize)
            import dp
            real(kind=dp), dimension(:,:), intent(in)  :: a, b
            real(kind=dp), dimension(:,:), intent(out) :: c
            integer, intent(in) :: blocksize
        end subroutine a_maal_b_blocks_interface
    end interface
    
    !--------------------------------------------------------------------------
    ! Main timing program
    !--------------------------------------------------------------------------
    integer :: k, N, idx_i, idx_j
    integer, parameter :: blocksize=100
    real :: flops
    real :: dummy_i, dummy_j
    integer, dimension(:), allocatable :: seed
    real(kind=dp), dimension(:,:), allocatable :: a, b, c
    real(kind=dp), dimension(:,:), allocatable :: c_matmul


    write(*,'(A10,5A15))') "", "lus", "dotprod", "Blas", "block", "matmul"
    do N = 100, 2000, 100
	    ! Make sure we use the same pseudo-random numbers each time by initializing
	    ! the seed to a certain value.
	    call random_seed(size=k)
	    allocate(seed(k))
	    seed = N
	    call random_seed(put=seed)

	    ! Calculate some random indices for the element c_ij that we are going to use
	    ! to check if the matrix computation went ok.
	    call random_number(dummy_i)
	    call random_number(dummy_j)
	    idx_i = floor(N*dummy_i) + 1
	    idx_j = floor(N*dummy_j) + 1

	    ! Allocate the matrices and one reference matrix
	    allocate(a(N,N), b(N,N), c(N,N), c_matmul(N,N))
	    call random_number(a)
	    call random_number(b)
	    call a_maal_b_matmul(a,b,c_matmul) ! Reference value
	    
	    write(*,'(I10)', advance="no") N
	
	    ! 1. Three nested loops
	    call do_timing( "JKI", a_maal_b_jki )
	    
	    ! 2. Two nested loops with dot_product and explicit transpose of matrix A
	    call do_timing( "IJ TP DOT_PRODUCT", a_maal_b_transp_ij_dot_product )
	    
	    ! 3. Using BLAS
	    call do_timing( "BLAS DGEMM", a_maal_b_blas )
	    
	    ! 4. In blocks
	    call do_timing( "IN BLOCKS", method_blocks=a_maal_b_blocks )
	    
	    ! 5. Intrinsic matmul function
	    call do_timing( "MATMUL", a_maal_b_matmul )

	    write(*,*)
	    
	    ! Clean up
	    deallocate(a, b, c, c_matmul,seed)
	enddo

contains

    subroutine do_timing( name, method, method_blocks )
        character(len=*), intent(in) :: name
        procedure(a_maal_b_interface), optional :: method
        procedure(a_maal_b_blocks_interface), optional :: method_blocks
        real(kind=dp) :: mynorm
        real :: t1, t2
        ! Do the timing
        if( present(method) ) then
            call cpu_time(t1)
            call method( a, b, c )
            call cpu_time(t2)
        else
            call cpu_time(t1)
            call method_blocks( a, b, c, blocksize)
            call cpu_time(t2)
        end if
        ! Compute the Frobenius norm of the error and divide by the Frobenius norm
        ! of the exact matrixproduct to get some kind of relative error.
        ! print on stderr if the norm is to big
        mynorm = sqrt(sum((c_matmul-c)**2))/sqrt(sum(c_matmul**2))
        if (mynorm >= 1d-10) write(0,'(A,A)') "There is an error in calculating the product with ", name

        write(*,'(F15.8)',advance="no") t2-t1

    end subroutine do_timing

end program timingtest 
