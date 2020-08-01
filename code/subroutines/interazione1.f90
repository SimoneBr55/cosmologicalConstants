SUBROUTINE interazione1(nObs, JD, nPts, acc, debug)
  IMPLICIT NONE
  !!------------------VARIABILI IN ENTRATA--------------------!!
  INTEGER, intent(in) :: nObs
  LOGICAL, intent(in) :: debug
  REAL*8, DIMENSION(:), intent(in) :: JD(nObs)
  !!----------------------------------------------------------!!

  !!------------------VARIABILI DI LAVORO---------------------!!
  !!----------------------------------------------------------!!

  !!------------------VARIABILI IN USCITA---------------------!!
  INTEGER, intent(out) :: nPts
  REAL*8, intent(out) :: acc
  !!----------------------------------------------------------!!

  IF(debug .EQV. .TRUE.) THEN
158  CONTINUE
     PRINT*, "In quanti punti si vuole interpolare la distribuzione?"
     PRINT*, "Si consideri che ", nObs, "dati coprono l'intervallo", JD(1), "-", JD(nObs)
     PRINT*, "Si consideri, inoltre, che il programma funziona fino a ", 416139997," punti (416 milioni)." 
     READ(*,*) nPts
     IF(nPts .LE. nObs) THEN
        PRINT*, "Non mi prenda in giro... Scelga qualche punto in più..."
        GOTO 158
     ELSEIF(nPts > 416139997) THEN
        PRINT*, "Troppi punti... c'è il rischio che, a causa della precisione finita del REAL, succedano cose strane..."
        GOTO 158
     ENDIF
     acc = (JD(nObs) - JD(1)) / (nPts - 1)
     PRINT*, "Si avrà una risoluzione temporale stimata di circa", acc, "JD."
     READ(*,*)
  ELSE
     acc = 1./1000.
     nPts = INT((JD(nObs) - JD(1))/acc + 1)
     PRINT*, "Scelgo, automaticamente, ", nPts, " punti, per avere una risoluzione temporale di circa ", acc, "JD."
  ENDIF
END SUBROUTINE interazione1
