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
type(FM) :: fm_p, fm_pi
type(mp_real) :: mp_p, mp_pi
call FM_SET(1000)
p = TO_FM('0.5')

call FM_PI(fm_pi)
mp_pi = mppic
k = 6

call calc_pi(fm_p, fm_pi, mp_p, mp_pi, k, N)

call calc_pi_fm(fm_p, fm_pi, k, N)
call calc_pi_mp(mp_p, mp_pi, k, N)

contains

	subroutine calc_pi(fm_p, fm_pi, mp_p, mp_pi, k, N)
		type(FM) :: fm_p, fm_pi
		type(mp_real) :: mp_p, mp_pi
		integer :: i, N
		double precision :: k

		i = 0
		do i = 1,N
			fm_p = sqrt( (1 - sqrt(1-fm_p*fm_p)) / 2)
			mp_p = sqrt( (1 - sqrt(1-mp_p*mp_p)) / 2)
			k = 2*k

			fm_cijfers = int(floor(-(log(abs(k*fm_p-fm_pi))/log(TO_FM('10')))))
			mp_cijfers = int(floor(-(log(abs(k*mp_p-mp_pi))/log('10')))
			write (*, '(2(i10)') fm_cijfers, mp_cijfers
		enddo
	end subroutine

	subroutine calc_pi_fm(fm_p, fm_pi, k, N)
		type(FM) :: fm_p, fm_pi
		integer :: i, N
		double precision :: k

		i = 0
		do i = 1,N
			fm_p = sqrt( (1 - sqrt(1-fm_p*fm_p)) / 2)
			k = 2*k

			fm_cijfers = int(floor(-(log(abs(k*fm_p-fm_pi))/log(TO_FM('10')))))
			write (*, '(i10)') fm_cijfers
		enddo
	end subroutine

	subroutine calc_pi_mp(mp_p, mp_pi, k, N)
		type(mp) :: mp_p, mp_pi
		integer :: i, N
		double precision :: k

		i = 0
		do i = 1,N
			mp_p = sqrt( (1 - sqrt(1-mp_p*mp_p)) / 2)
			k = 2*k

			mp_cijfers = int(floor(-(log(abs(k*mp_p-mp_pi))/log('10'))))
			write (*, '(i10)') mp_cijfers
		enddo
	end subroutine


end program
