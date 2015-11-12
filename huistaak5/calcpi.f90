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

implicit none

integer, parameter :: N=200
double precision :: k
type(FM) :: p
p = TO_FM('0.5')
call FM_SET(1000)

k = 6
call calc_zero(p,k,iters)

contains

	subroutine calc_pi(p, k, N)
		type(FM) :: p 
		integer :: i

		i = 0
		WRITE (*,'("iteratie: ", i0, "  pi= ",A)') i,TRIM(FM_FORMAT('E50.40',p))
		do i = 1,N
			p = sqrt( (1 - sqrt(1-p^2)) / 2)
			k = 2*k
			WRITE (*,'("iteratie: ", i0, "  pi= ",A)') i,TRIM(FM_FORMAT('E50.40',p*k))
		enddo
	end subroutine

end program
