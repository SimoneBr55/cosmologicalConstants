SUBROUTINE soluzSist(inc, tn, nObs, der, derC, debug)
  IMPLICIT NONE
  !!--------------------------VARIABILI IN ENTRATA--------------------------!!
  LOGICAL, intent(in) :: debug
  INTEGER, intent(in) :: nObs
  REAL*8, DIMENSION(:,:), intent(in) :: inc(nObs-2, nObs-2)
  REAL*8, DIMENSION(:), intent(in) :: tn(nObs-2)
  !!------------------------------------------------------------------------!!

  !!----------------------------VARIABILI DI LAVORO-------------------------!!
  INTEGER :: i, j, k
  REAL*8, ALLOCATABLE :: incCopy(:,:), tnCopy(:)
  REAL*8 :: sum !BE CAREFUL
  REAL*8 :: aux, fakt
  !!------------------------------------------------------------------------!!

  !!-----------------------------VARIABILI DI USCITA------------------------!!
  REAL*8, DIMENSION(:), intent(out) :: der(nObs-2), derC(nObs)
  !!------------------------------------------------------------------------!!

  !!-----------------------------AZZERAMENTO--------------------------------!!
  !!------------------------------------------------------------------------!!

  !!-----------------------------ALLOCAMENTO--------------------------------!!
  ALLOCATE(incCopy(nObs-2, nObs-2), tnCopy(nObs-2))
  !!------------------------------------------------------------------------!!

  incCopy = inc
  tnCopy = tn

  DO i=1,nObs-2
     aux = incCopy(i,i)
     DO k=1,nObs-2
        incCopy(i,k) = incCopy(i,k) / aux
     ENDDO
     tnCopy(i) = tnCopy(i)/aux

     DO j=1,nObs-2
        fakt = incCopy(j,i)
        IF(i/=j)THEN
           DO k=1,nObs-2
              incCopy(j,k) = incCopy(j,k)-fakt*incCopy(i,k)
           ENDDO
           tnCopy(j) = tnCopy(j) - fakt*tnCopy(i)
        ENDIF
     ENDDO
  ENDDO
  der = tnCopy

  derC(1) = 0.
  derC(nObs) = 0.
  DO i=2,nObs-1
     derC(i) = der(i-1)
  ENDDO
  


  !!--------------------------DEALLOCATE---------------------------!!
  DEALLOCATE(incCopy, tnCopy)
  !!---------------------------------------------------------------!!

  

END SUBROUTINE soluzSist

