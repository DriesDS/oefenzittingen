! MATRIXOP
!   Implementation of several square matrix-matrix products
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
module matrixop
    implicit none
    integer, parameter :: dp = selected_real_kind(15,307)
contains
    !--------------------------------------------------------------------------
    ! 1. Three nested loops
    !
    ! NOTE: use the following convention for the indices
    !       i = row index of A
    !       j = column index of B
    !       k = column index of A
    !--------------------------------------------------------------------------

    subroutine a_maal_b_ijk(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N, i, j, k
        N = size(a,1)
        c=0d0
        do i = 1,N
            do j = 1,N
                do k = 1,N
                    c(i,j) = c(i,j) + a(i,k)*b(k,j)
                enddo
            enddo
        enddo
    end subroutine a_maal_b_ijk

    subroutine a_maal_b_ikj(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N, i, j, k
        N = size(a,1)
        c=0d0
        do i = 1,N
            do k = 1,N
                do j = 1,N
                    c(i,j) = c(i,j) + a(i,k)*b(k,j)
                enddo
            enddo
        enddo
    end subroutine a_maal_b_ikj

    subroutine a_maal_b_jik(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N, i, j, k
        N = size(a,1)
        c=0d0
        do j = 1,N
            do i = 1,N
                do k = 1,N
                    c(i,j) = c(i,j) + a(i,k)*b(k,j)
                enddo
            enddo
        enddo
    end subroutine a_maal_b_jik

    subroutine a_maal_b_jki(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N, i, j, k
        N = size(a,1)
        c=0d0
        do j = 1,N
            do k = 1,N
                do i = 1,N
                    c(i,j) = c(i,j) + a(i,k)*b(k,j)
                enddo
            enddo
        enddo
    end subroutine a_maal_b_jki

    subroutine a_maal_b_kij(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N, i, j, k
        N = size(a,1)
        c=0d0
        do k = 1,N
            do i = 1,N
                do j = 1,N
                    c(i,j) = c(i,j) + a(i,k)*b(k,j)
                enddo
            enddo
        enddo
    end subroutine a_maal_b_kij
    
    subroutine a_maal_b_kji(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N, i, j, k
        N = size(a,1)
        c=0d0
        do k = 1,N
            do j = 1,N
                do i = 1,N
                    c(i,j) = c(i,j) + a(i,k)*b(k,j)
                enddo
            enddo
        enddo
    end subroutine a_maal_b_kji
    !--------------------------------------------------------------------------
    ! 2. Two nested loops with vector operations
    !--------------------------------------------------------------------------

    subroutine a_maal_b_ikj_vect(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N, i, k
        N = size(a,1)
        c=0d0
        do i = 1,N
            do k = 1,N
                c(i,:) = c(i,:) + a(i,k)*b(k,:)
            enddo
        enddo
     end subroutine a_maal_b_ikj_vect

    subroutine a_maal_b_jki_vect(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N, j, k
        N = size(a,1)
        c=0d0
        do j = 1,N
            do k = 1,N
                c(:,j) = c(:,j) + a(:,k)*b(k,j)
            enddo
        enddo
    end subroutine a_maal_b_jki_vect

    subroutine a_maal_b_kij_vect(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N, i, k
        N = size(a,1)
        c=0d0
        do k = 1,N
            do i = 1,N
                c(i,:) = c(i,:) + a(i,k)*b(k,:)
            enddo
        enddo
    end subroutine a_maal_b_kij_vect

    subroutine a_maal_b_kji_vect(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N, j, k
        N = size(a,1)
        c=0d0
        do k = 1,N
            do j = 1,N
                c(:,j) = c(:,j) + a(:,k)*b(k,j)
            enddo
        enddo
    end subroutine a_maal_b_kji_vect
    !--------------------------------------------------------------------------
    ! 3. Two nested loops with dot_product
    !--------------------------------------------------------------------------

    subroutine a_maal_b_ij_dot_product(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N, j, i
        N = size(a,1)
        c=0d0
        do i = 1,N
            do j = 1,N
                c(i,j) = dot_product(a(i,:),b(:,j))
            enddo
        enddo
    end subroutine a_maal_b_ij_dot_product

    subroutine a_maal_b_ji_dot_product(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N, j, i
        N = size(a,1)
        c=0d0
        do j = 1,N
            do i = 1,N
                c(i,j) = dot_product(a(i,:),b(:,j))
            enddo
        enddo
    end subroutine a_maal_b_ji_dot_product
    !--------------------------------------------------------------------------
    ! 4. Two nested loops with dot_product and explicit transpose of matrix A
    !--------------------------------------------------------------------------

    subroutine a_maal_b_transp_ij_dot_product(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        real(kind=dp), dimension(size(a,1),size(a,1)) :: at
        integer :: N, j, i
        N = size(a,1)
        at = transpose(a)
        c=0d0
        do j = 1,N
            do i = 1,N
                c(i,j) = dot_product(at(:,i),b(:,j))
            enddo
        enddo
    end subroutine a_maal_b_transp_ij_dot_product

    subroutine a_maal_b_transp_ji_dot_product(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        real(kind=dp), dimension(size(a,1),size(a,1)) :: at
        integer :: N, j, i
        N = size(a,1)
        at = transpose(a)
        c=0d0
        do j = 1,N
            do i = 1,N
                c(i,j) = dot_product(at(:,i),b(:,j))
            enddo
        enddo
    end subroutine a_maal_b_transp_ji_dot_product
    !--------------------------------------------------------------------------
    ! 5. Using BLAS : Add library in linking phase
    !--------------------------------------------------------------------------

    subroutine a_maal_b_blas(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer :: N
        c = 0.0_dp
        N=size(a,1)
        call dgemm('N','N',N,N,N,1d0,a,N,b,N,0d0,c,N)
    end subroutine a_maal_b_blas
   
    !--------------------------------------------------------------------------
    ! 6. In blocks
    !--------------------------------------------------------------------------

    subroutine a_maal_b_blocks(a, b, c, blocksize)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer, intent(in) :: blocksize
        real(kind=dp), dimension(blocksize,blocksize) :: cloc, aloc, bloc, clocadd
        integer :: N, i, j, k
        N=size(a,1)
        c = 0.0_dp
        do i=1,N/blocksize
            do j=1,N/blocksize
                cloc=0
                do k=1,N/blocksize
                    aloc = a( (i-1)*blocksize+1:i*blocksize,(k-1)*blocksize+1:k*blocksize )
                    bloc = b( (k-1)*blocksize+1:k*blocksize,(j-1)*blocksize+1:j*blocksize )
                    clocadd=0
                    call a_maal_b_kji( aloc, bloc, clocadd )
                    cloc = cloc + clocadd
                enddo
                c((i-1)*blocksize+1:i*blocksize,(j-1)*blocksize+1:j*blocksize) = cloc
            enddo
        enddo
    end subroutine a_maal_b_blocks
    !--------------------------------------------------------------------------
    ! 7. Intrinsic matmul function
    !--------------------------------------------------------------------------

    subroutine a_maal_b_matmul(a, b, c)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        c = matmul( a, b ) ! Already completed
    end subroutine a_maal_b_matmul
    !--------------------------------------------------------------------------
    ! 8. Eigen matrixprod function
    !--------------------------------------------------------------------------

    subroutine a_maal_b_eigen(a, b, c, blocksize)
        real(kind=dp), dimension(:,:), intent(in)  :: a, b
        real(kind=dp), dimension(:,:), intent(out) :: c
        integer, intent(in) :: blocksize
        real(kind=dp), dimension(size(a,1),size(a,1)) :: at
        real(kind=dp), dimension(size(a,1), blocksize**2/size(a,1)) :: aloc, bloc
        real(kind=dp), dimension(blocksize**2/size(a,1), blocksize**2/size(a,1)) :: cloc
        integer :: N, i, j, k, l, rows
        N = size(a,1)
        rows = blocksize**2/N
        at = transpose(a)
        c = 0.0_dp
        do i=1,N/rows
            bloc = b(:,(i-1)*rows+1:i*rows)
            do j=1,N/rows
                aloc = at(:,(i-1)*rows+1:i*rows)
                cloc = 0d0
                do k = 1,rows
                    do l = 1,rows
                        cloc(l,k) = dot_product(aloc(:,l),bloc(:,k))
                    enddo
                enddo
                c((j-1)*rows+1:j*rows,(i-1)*rows+1:i*rows) = cloc
            enddo
        enddo
    end subroutine a_maal_b_eigen

end module matrixop
