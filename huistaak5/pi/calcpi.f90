! We hebben 2 rechthoekige driehoeken.
! Driehoek 1 met zijden van lengte 1/2, s_k/2 en a
! Driehoek 2 met zijden van lengte s_2k, s_k/2 en b
! met respectievelijk zijden met lengte 1/2 en s_2k de rechthoekszijden.
! we hebben 3 onbekenden: s_2k, a en b.
! we hebben 3 vergelijkingen: - a+b = 1/2                    (1)
!                             - (1/2)^2 = a^2+(s_k/2)^2      (2)
!                             - (s_2k)^2 = b^2+(s_k/2)^2     (3)
! waarbij we de twee laatste vergelijkingen bekomen door pythagoras toe
! te passen op de twee rechthoekige driehoeken.
! We kunnen dit stelsel vergelijkingen oplossen door twee maal substitutie
! toe te passen.
! 1: (2) in functie van a schrijven en te substitueren in (1)
! 2: (1) nu in functie van b schrijven en substitueren in (3)
! 
! 1: a = sqrt( 1-s_k^2 ) / 2  =>  b = (1 - sqrt( 1-s_k^2 )) / 2
! 2: s_2k^2 = ((1 - sqrt( 1-s_k^2 )) / 2)^2 + s_k^2 / 4
! 
! en door de kwadraten allemaal uit te schrijven bekomen we ten slotte:
! s_2k = sqrt( (1 - sqrt(1-s_k^2)) / 2)

program calcpi

use FMZM
use mpmodule

implicit none

integer, parameter :: N=200
double precision :: k
type(FM) :: p, pi, fmk
type(mp_real) :: mpp, mppi, mpfmk
p = TO_FM('0.5')
mpp = 
call FM_SET(1000)
pi = 4*atan(TO_FM('1.0'))

k = 6
fmk = TO_FM('6')
call calc_pi(p,k,fmk,N,pi)

contains

	subroutine calc_pi(p, k, fmk, N, pi)
		type(FM) :: p, kp, pi, fmk
		integer :: i, N
		double precision :: k

		i = 0
		write (*,'("iteratie: ", i0, "  pi= ",A)') i,TRIM(FM_FORMAT('E50.40',p))
		do i = 1,N
			p = sqrt( (1 - sqrt(1-p*p)) / 2)
			k = 2*k
			fmk = 2*fmk
			!write (*,'("iteratie: ", i0, "  pi= ",A)') i,TRIM(FM_FORMAT('E100.90',k*p))
			write (*,'("difference with pi: ", A)') TRIM(FM_FORMAT('E15.5',fmk*p-pi))
			write (*,'(e15.6,a)') k, trim(FM_FORMAT('E30.20',p))
		enddo
	end subroutine

end program
