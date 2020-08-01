SUBROUTINE sorting(nObs, JDR, magR, JD, mag)
  IMPLICIT NONE
  !!-----------------VARIABILI DI INGRESSO---------------!!
  INTEGER, intent(in) :: nObs
  REAL*8, DIMENSION(:), intent(in) :: JDR(nObs), magR(nObs)
  !!-----------------------------------------------------!!
  !!-----------------VARIABILI DI LAVORO-----------------!!
  REAL*8 :: park, park2, minimo, corresp
  INTEGER :: pos_min
  INTEGER :: i, j
  !!-----------------------------------------------------!!
  !!-----------------VARIABILI DI USCITA-----------------!!
  REAL*8, DIMENSION(:), intent(out) :: JD(nObs), mag(nObs)
  !!-----------------------------------------------------!!

  !!---------------------ALLOCAZIONE---------------------!!
  
  !!-----------------------------------------------------!!
  !!---------------------AZZERAMENTO---------------------!!
  park = 0.
  park2 = 0.
  minimo = 0.
  corresp = 0.
  pos_min = 0

  JD = JDR
  mag = magR
  !!-----------------------------------------------------!!

  DO i=1,nObs-1
     minimo = JD(i)
     corresp = JD(i)
     pos_min = i
     DO j=i+1, nObs
        IF(JD(j) < minimo) THEN
           minimo = JD(j)
           corresp = mag(j)
           pos_min = j
        ENDIF
     ENDDO
     park = JD(i)
     park2 = mag(i)
     JD(i) = minimo
     mag(i) = corresp
     JD(pos_min) = park
     mag(pos_min) = park2
  ENDDO

END SUBROUTINE sorting
