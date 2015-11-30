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
    integer :: i, k, N, idx_i, idx_j
    integer, parameter :: blocksize=100
    real :: flops
    real :: dummy_i, dummy_j
    integer, dimension(6) :: Ns
    integer, dimension(:), allocatable :: seed
    integer, parameter :: nb=5
    real(kind=dp), dimension(:,:), allocatable :: a, b, c
    real(kind=dp), dimension(:,:), allocatable :: c_matmul
    real(kind=dp), dimension(6,6,3) :: timings

    Ns = (/ 100, 200, 400, 500, 1000, 2000 /)
    write(*,'(A10,6(A12))') "", "lus", "dotprod", "Blas", "block", "matmul", "eigen"
    do i = 1,6
        N = Ns(i)

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
	    call do_timing( "JKI", timings(i,1,:), a_maal_b_jki )
	    
	    ! 2. Two nested loops with dot_product and explicit transpose of matrix A
	    call do_timing( "IJ TP DOT_PRODUCT", timings(i,2,:), a_maal_b_transp_ij_dot_product )
	    
	    ! 3. Using BLAS
	    call do_timing( "BLAS DGEMM", timings(i,3,:), a_maal_b_blas )
	    
	    ! 4. In blocks
	    call do_timing( "IN BLOCKS", timings(i,4,:), method_blocks=a_maal_b_blocks )
	    
	    ! 5. Intrinsic matmul function
	    call do_timing( "MATMUL", timings(i,5,:), a_maal_b_matmul )
	
	    ! 6. eigen ontworpen function
	    call do_timing( "EIGEN", timings(i,6,:), method_blocks=a_maal_b_eigen )
	    
	    ! Clean up
	    deallocate(a, b, c, c_matmul,seed)
	
	    write(*,'(6(f12.8))') timings(i,:,2)
	
	    write(0,*) N
	
	enddo

	write(*,*) "maximum tijden:"
	do i = 1,6
		write(*,'(I10,6(f12.8))') Ns(i), timings(i,:,3)
	enddo
	write(*,*) "minimum tijden:"
	do i = 1,6
		write(*,'(I10,6(f12.8))') Ns(i), timings(i,:,1)
	enddo

contains

    subroutine do_timing( name, times, method, method_blocks )
        character(len=*), intent(in) :: name
        procedure(a_maal_b_interface), optional :: method
        procedure(a_maal_b_blocks_interface), optional :: method_blocks
        real(kind=dp) :: mynorm
        real(kind=dp), dimension(3) :: times
        real :: t1, t2, totaltime, maxtime, mintime
        integer :: i
        ! Do the timing 3 times
        totaltime = 0
        maxtime = 0
        mintime = 1d10
        do i = 1,nb
	        if( present(method) ) then
	            call cpu_time(t1)
	            call method( a, b, c )
	            call cpu_time(t2)
	        else
	            call cpu_time(t1)
	            call method_blocks( a, b, c, blocksize)
	            call cpu_time(t2)
	        end if
	        totaltime = totaltime + t2-t1
	        if (t2-t1 < mintime) mintime = t2-t1
	        if (t2-t1 > maxtime) maxtime = t2-t1
	    enddo
        ! Compute the Frobenius norm of the error and divide by the Frobenius norm
        ! of the exact matrixproduct to get some kind of relative error.
        ! print on stderr if the norm is to big
        mynorm = sqrt(sum((c_matmul-c)**2))/sqrt(sum(c_matmul**2))
        if (mynorm >= 1d-10) write(0,'(A,A,A,f15.5)') "There is an error in calculating the product with " &
                , name, "rel err: ", mynorm
        

        times(1) = mintime
        times(2) = totaltime/nb
        times(3) = maxtime

    end subroutine do_timing

end program timingtest 
