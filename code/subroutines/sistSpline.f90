SUBROUTINE sistSpline(nObs, JD, mag, inc, tn, debug)
  IMPLICIT NONE
  !!-----------------VARIABILI DI INGRESSO---------------!!
  LOGICAL, intent(in) :: debug
  REAL*8, DIMENSION(:), intent(in) :: JD(nObs), mag(nObs)
  INTEGER, intent(in) :: nObs
  
  !!-----------------------------------------------------!!
  !!-----------------VARIABILI DI LAVORO-----------------!!
  INTEGER :: i, j, k
  REAL*8, ALLOCATABLE :: incC(:,:), tnC(:)
  REAL*8, ALLOCATABLE :: h(:), b(:)
  !!-----------------------------------------------------!!
  !!-----------------VARIABILI DI USCITA-----------------!!
  REAL*8, DIMENSION(:,:), intent(out) :: inc(nObs-2,nObs-2)
  REAL*8, DIMENSION(:), intent(out) :: tn(nObs-2)
  !!-----------------------------------------------------!!

  !!---------------------ALLOCAZIONE---------------------!!
  ALLOCATE(incC(nObs-2, nObs), tnC(nObs-2))
  ALLOCATE(h(nObs-1), b(nObs-1))
  !!-----------------------------------------------------!!
  !!---------------------AZZERAMENTO---------------------!!
  incC = 0.
  tnC = 0.
  h = 0.
  b = 0.
  !!-----------------------------------------------------!!

  DO i= 1,nObs-1
     h(i) = JD(i+1) - JD(i)
     b(i) = (mag(i+1) - mag(i))/h(i)
  ENDDO

  DO i=1,nObs-2
     incC(i,i) = h(i)
     incC(i,i+1) = 2*(h(i) + h(i+1))
     incC(i,i+2) = h(i+1)
     tnC(i) = 6*(b(i+1) - b(i))
  ENDDO

  tn = tnC

  correzioneMatrice: DO i=1,nObs-2
     DO j=1,nObs-2
        inc(i,j) = incC(i,j+1)
     ENDDO
  ENDDO correzioneMatrice
     
  !!--------------------DEALLOCAMENTO--------------------!!
  DEALLOCATE(incC, tnC)
  DEALLOCATE(h, b)
  !!-----------------------------------------------------!!
  
  !!----------------------AZZERAMENTO--------------------!!
  
  !!-----------------------------------------------------!!
     

END SUBROUTINE SistSpline
