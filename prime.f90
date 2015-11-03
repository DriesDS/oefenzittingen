! This program searches for the n'th prime, but not in the most efficient way...
program prime

implicit none
integer :: n, i, j
integer, parameter :: controlen=10000, controlewaarde=104729
integer, dimension(:), allocatable :: primeList
logical :: isPrime

!write(unit=*, fmt="(A)", advance="no") "Enter the value of n: "
!read *, n
n = 200000

allocate(primeList(n))

primeList(1:2) = (/2,3/)
if (n<100) write(unit=*, fmt="(A, I0, A, I0)") "Prime ", 1, ": ", primeList(1)

do i = 3,n
	primeList(i) = primeList(i-1)
	do
		! enkel over de oneven getallen itereren
		primeList(i) = primeList(i) + 2
		isPrime = .true.
		! niet meer controleren op 2, want we werken enkel met oneven getallen
		do j = 2,i-1
			if (modulo(primeList(i), primeList(j)) == 0) then
				isPrime = .false.
			endif
			! als een priemgetal groter dan dit een deler is van het getal dat we controleren
			! is het deeltal kleiner dan dit priemgetal en ook een deler 
			! (p(i)%p(j) = 0 ==> p(i)/p(j) < p(j) en p(i)/p(j) is ook een deler met een priemfactor < p(j)
			if (primeList(j)**2 > primelist(i)) exit
		end do
		if (isPrime) then
			if (n<100) write(unit=*,fmt="(A, I0, A, I0)") "Prime ", i, ": ", primeList(i)
			exit
		end if
	end do
end do

write(unit=*, fmt="(A, i0, A, I0, A, L, A)") "! Prime ", n, " is ", primeList(N), ", check = ", primeList(controlen) == controlewaarde

deallocate(primeList)

end program prime
