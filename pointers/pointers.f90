!	Dries De Samblanx	1e master Wiskundige ingenieurstechnieken
!
! 	getest met ifort	($ ifort pointers.f90	$ ./a.out)
!
! ANTWOORD:
! 
! Aan de hand van deze code is het verschil duidelijk te zien tussen een deelmatrix c
! en een pointer p. Dezelfde code wordt uitgevoerd met enerzijds behulp van een pointer
! en anderzijds van een deelmatrix.
! deelmatrix:
! 	Wanneer we een deelmatrix gebruiken om te verwijzen naar matrix (8:9,5:7) dan
!	wordt er een kopie gemaakt van al deze elementen. Dit wil zeggen dat we in het
!	geheugen zowel een matrix 'matrix' van 70 elementen hebben als een matrix d van
!	6 elementen. Wanneer we dan een aanpassing doen aan matrix d, doen we dit in het
!	geheugen ook enkel op deze matrix en zien we geen aanpassing in de matrix 'matrix'
!	We hebben dus te maken met een duplicatie.
! pointer:
!	Wanneer we een pointer gebruiken om naar deze deelmatrix te verwijzen, dan wordt
!	er geen kopie gemaakt. In plaats daarvan wijst de pointer-variabele naar dezelfde
!	plaats in het geheugen als de oorsprongkelijke matrix. Wanneer er nu in 1 van deze
!	twee matrices een aanpassing wordt gedaan, zien we dit ook in de andere matrix.
!	Er is hier dus minder geheugen gebruikt, maar we moeten hier goed oppassen dat we
!	geen aanpassing doen in de ene variabele, als we willen dat de andere ongewijzigd
!	blijft.

program pointers

implicit none

real, dimension(10,7), target :: matrix
real, dimension(2,3)  :: d
real,  pointer  :: p(:,:)
integer :: i,j
real, parameter :: e=10**(-5)
logical, parameter :: point=.true.

matrix = reshape((/ ( ( ((i-1)*7.0+j)/70.0, i=1,10), j=1,7) /), (/ 10,7 /))

if (point) then
	d = matrix(8:9,5:7)
else
	p => matrix(8:9,5:7)
endif

if (point) then
	write(*,*) "!", (minval(d) - (7-3+1+7*(8-1))/70.0)<e
	d(1,1) = 2
	write(*,*) "!", matrix(8,5)
else
	write(*,*) "!", (minval(p) - (7-3+1+7*(8-1))/70.0)<e
	p(1,1) = 2
	write(*,*) "!", matrix(8,5)
endif

end program
