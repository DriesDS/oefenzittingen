! 	Dries De Samblanx	1e master wiskundige ingenieurstechnieken
! 
! 	getest met nagfor, ifort, gfortran, g95, mpif90.
! 
!	commando's: $ ifort prime.f90	$ time ./a.out n  -- (met n het priemgetal Pn dat je wilt berekenen) 
!
! 	geschatte tijd (huistaak + voorbereiden oefenzitting): 5u
! 
! ANTWOORD:
! 
! Wanneer we het programma opkuisen (enkele onnuttige berekeningen weglaten, eendimensionale
! array voor primeList, ... kunnen we berekenen tot ongeveer P8000.
! Prime 8000 is 81773, check = F
! De check voor priemgetal 10000 is hier False omdat we niet zover geraakt zijn.
!
! (1) De volgende verbetering kunnen we maken door in te zien dat vanaf 3 er enkel oneven priemgetallen
! zijn. Wanneer we dus de lijst met priemgetallen beginnen met 2 en 3 en vanaf dan altijd +2 doen
! ipv +1 om het volgende getal te testen, krijgen we al een iets efficienter programma.
! Hiermee kunnen we binnen de seconde al tot priemgetal 11000 gaan, een verbetering van 37.5%!
! Prime 11000 is 116447, check = T
! Zoals je ziet is hier de check True, wat wil zeggen dat priemgetal 10000 juist berekend was. 
!
! (2) Nog een kleine aanpassing die we kunnen doen bouwt voort op de vorige. Aangezien we enkel nog
! oneven getallen testen, moeten we ook niet meer testen op deelbaarheid door 2. Dat wil zeggen dat
! we 2 altijd mogen overslaan in de modulo-berekening. Dit zorgt niet voor een impressionante vooruitgang
! maar we kunnen toch al berekenen tot 11250 ipv 11000.
! Prime 11250 is 119429, check = T
! Ook hier is de check nog altijd True.
!
! (3) Een aanpassing die wel voor een grote vooruitgang zorgt bekomen we wanneer we inzien dat we
! een getal n maar moeten testen op delers tot sqrt(n). Wanneer n immers een deler zou hebben
! groter dan dit getal wilt dit ook zeggen dat het een deler heeft kleiner dan dit getal.
! (p(i)%p(j) = 0 ==> p(i)/p(j) < p(j) en p(i)/p(j) is ook een deler met een priemfactor < p(j)
! Met deze aanpassing kunnen we berekenen tot p225000, het 20-voud van wat we konden voor deze
! aanpassing.
! Prime 225000 is 3122321, check = T (ook hier is het 10000e priemgetal goed berekend.)
!
! (4) Een laatste aanpassing die we doen heeft enkel met de programma-structuur te maken en niets met
! de wiskunde achter het vraagstuk. Vroeger deden we bij elke modulo-berekening een test of het 
! priemgetal groter was dan sqrt(n) of niet. Wanneer we nu bijhouden tot hoever we mogen testen, en elke
! buitenste lus 1 keer testen of we nog altijd ver genoeg testen, doen we veel minder sqrt() berekeningen.
! Dit heeft zeker geen verwaarloosbaar effect, en we kunnen nu berekenen tot 240000.
! Zeker als we er rekening mee houden dat de eerste priemgetallen veel sneller berekend worden dan de
! latere is dit zeker een goede vooruitgang.
! Prime 240000 is 3346601, check = T

program prime

implicit none
integer :: n, i, j, d
integer, parameter :: controlen=10000, controlewaarde=104729
integer, dimension(:), allocatable :: primeList
logical :: isPrime
character(len=8) :: primen

call getarg(1,primen)
read(primen,*) n 

allocate(primeList(n))

primeList(1:2) = (/2,3/)

d = 1
do i = 3,n
	primeList(i) = primeList(i-1)
	do
		! (1) enkel over de oneven getallen itereren
		primeList(i) = primeList(i) + 2
		isPrime = .true.
		! (4) enkel deze test doen in de buitenste loop, en de binnenste loop maar tot d laten lopen
		if (primeList(i) > primeList(d)**2) d = d+1
		! (2) niet meer controleren op 2, want we werken enkel met oneven getallen
		do j = 2,d
			if (modulo(primeList(i), primeList(j)) == 0) then
				isPrime = .false.
			endif
			! (3) als een priemgetal groter dan dit een deler is van het getal dat we controleren
			! is het deeltal kleiner dan dit priemgetal en ook een deler 
			! if (primeList(j)**2 > primelist(i)) exit
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
