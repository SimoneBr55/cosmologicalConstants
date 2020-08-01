SUBROUTINE countSingleData(nFile, nObs)
  IMPLICIT NONE
  !!-----------------VARIABILI DI INGRESSO---------------!!
  INTEGER, intent(in) :: nFile
  !!-----------------------------------------------------!!
  !!-----------------VARIABILI DI LAVORO-----------------!!
  
  !!-----------------------------------------------------!!
  !!-----------------VARIABILI DI USCITA-----------------!!
  INTEGER, intent(out) :: nObs
  !!-----------------------------------------------------!!

  !!---------------------ALLOCAZIONE---------------------!!
  
  !!-----------------------------------------------------!!
  !!---------------------AZZERAMENTO---------------------!!
  nObs = 0
  !!-----------------------------------------------------!!

  DO
     READ(nFile,*, END=998)
     nObs = nObs + 1
  ENDDO

998 CONTINUE
  REWIND nFile

END SUBROUTINE countSingleData
