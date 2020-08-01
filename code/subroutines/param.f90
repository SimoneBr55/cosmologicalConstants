SUBROUTINE param(nCeph, magnCat, perCat, c1, c2, debug)
  IMPLICIT NONE
  !!-----------------VARIABILI DI INGRESSO---------------!!
  LOGICAL, intent(in) :: debug
  INTEGER, intent(in) :: nCeph
  !!-----------------------------------------------------!!
  !!-----------------VARIABILI DI LAVORO-----------------!!
  INTEGER :: i
  REAL*8, ALLOCATABLE :: logP(:)
  REAL*8 :: sx2, sy, sx, sxy, delta
  !!-----------------------------------------------------!!
  !!-----------------VARIABILI DI USCITA-----------------!!
  REAL*8, DIMENSION(:), intent(out) :: magnCat(nCeph), perCat(nCeph)
  REAL*8, intent(out) :: c1, c2
  !!-----------------------------------------------------!!

  !!---------------------ALLOCAZIONE---------------------!!
  ALLOCATE(logP(nCeph))
  !!-----------------------------------------------------!!
  !!---------------------AZZERAMENTO---------------------!!
  sx2 = 0.
  sy = 0.
  sx = 0.
  sxy = 0.
  delta = 0.
  c1 = 0. !!
  c2 = 0. !!

  perCat = 0.
  magnCat = 0.
  logP = 0.
  !!-----------------------------------------------------!!

  ! leggo i valori da file e carico in matrice

  REWIND 120

  DO i=1,nCeph
     READ(120,*) perCat(i), magnCat(i)
  ENDDO

  DO i=1,nCeph
     logP(i) = LOG10(perCat(i))
  ENDDO

  ! Fittiamo equazione M = c1*logP + c2

  DO i=1,nCeph
     sx2 = sx2 + (logP(i))**2
     sy = sy + magnCat(i)
     sx = sx + logP(i)
     sxy = sxy + logP(i)*magnCat(i)
  ENDDO
  delta = nCeph*(sx2) - (sx)**2
  IF(debug .EQV. .TRUE.) THEN
     PRINT*, ""
     PRINT*, "SOMMA DEI QUADRATI DI logP ", sx2
     PRINT*, "SOMMA DELLE M ", sy
     PRINT*, "SOMMA DEI logP ", sx
     PRINT*, "SOMMA DEI * ", sxy
     PRINT*, "DENOMINATORE ", delta
  ENDIF  

  c1 = (nCeph*sxy - sx*sy)/(delta)
  c2 = (sx2*sy - sx*sxy)/(delta)

  DEALLOCATE(logP)

END SUBROUTINE param

