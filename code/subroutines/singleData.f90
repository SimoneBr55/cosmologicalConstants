SUBROUTINE singleData(nFile, nObs, JDR, magR)
  IMPLICIT NONE
  !!-----------------VARIABILI DI INGRESSO---------------!!
  INTEGER, intent(in) :: nFile, nObs
  !!-----------------------------------------------------!!
  !!-----------------VARIABILI DI LAVORO-----------------!!
  INTEGER :: i
  !!-----------------------------------------------------!!
  !!-----------------VARIABILI DI USCITA-----------------!!
  REAL*8, DIMENSION(:) :: JDR(nObs), magR(nObs)
  !!-----------------------------------------------------!!

  !!---------------------ALLOCAZIONE---------------------!!
  
  !!-----------------------------------------------------!!
  !!---------------------AZZERAMENTO---------------------!!
  JDR = 0.
  magR = 0.
  !!-----------------------------------------------------!!

  DO i=1,nObs
     READ(nFile,*) JDR(i), magR(i)
  ENDDO
  CLOSE(nFile)

END SUBROUTINE singleData
