SUBROUTINE interpolazione(nObs, nPts, acc, JD, mag, derC, pts, interp, debug)
  IMPLICIT NONE
  !!------------------------VARIABILI DI ENTRATA-----------------------!!
  LOGICAL, intent(in) :: debug
  INTEGER, intent(in) :: nObs, nPts
  REAL*8, DIMENSION(:), intent(in) :: JD(nObs), mag(nObs), derC(nObs)
  REAL*8, intent(in) :: acc
  !!-------------------------------------------------------------------!!

  !!------------------------VARIABILI DI LAVORO------------------------!!
  INTEGER :: i, j
  REAL*8, ALLOCATABLE :: delta(:), posDelta(:), negDelta(:)
  INTEGER, ALLOCATABLE :: pos_inf(:), pos_sup(:)
  REAL*8, ALLOCATABLE :: interp1(:), interp2(:), interp3(:), interp4(:)
  !!-------------------------------------------------------------------!!

  !!------------------------VARIABILI DI USCITA------------------------!!
  REAL*8, DIMENSION(:), intent(out) :: interp(nPts), pts(nPts)
  !!-------------------------------------------------------------------!!

  !!--------------------------AZZERAMENTO------------------------------!!
  interp = 0.
  !!-------------------------------------------------------------------!!

  pts(1) = JD(1)
  pts(nPts) = JD(nObs)
  DO i=2,nPts-1
     pts(i) = pts(1) + (DBLE(i-1)*acc)
  ENDDO

  ALLOCATE(delta(nPts), posDelta(nPts), negDelta(nPts))
  ALLOCATE(pos_inf(nPts), pos_sup(nPts))
  ALLOCATE(interp1(nPts), interp2(nPts), interp3(nPts), interp4(nPts))
  
  DO i=1,nPts
     DO j=2,nObs
        IF(pts(i) <= JD(j)) THEN
           pos_inf(i) = j-1
           pos_sup(i) = j
           EXIT
        ENDIF
     ENDDO
  ENDDO

  DO i=1,nPts
     delta(i) = JD(pos_sup(i)) - JD(pos_inf(i))
     posDelta(i) = JD(pos_sup(i)) - pts(i)
     negDelta(i) = pts(i) - JD(pos_inf(i))
  ENDDO

  DO i=1,nPts
     interp1(i) = (derC(pos_inf(i))/(6*delta(i)))*(posDelta(i))**3
     interp2(i) = (derC(pos_sup(i))/(6*delta(i)))*(negDelta(i))**3
     interp3(i) = ((mag(pos_inf(i)))/(delta(i)) - (derC(pos_inf(i))*(delta(i)))/6)*(posDelta(i))
     interp4(i) = ((mag(pos_sup(i)))/(delta(i)) - (derC(pos_sup(i))*(delta(i)))/6)*(negDelta(i))
  ENDDO

  DO i=1,nPts
     interp(i) = interp1(i) + interp2(i) + interp3(i) + interp4(i)
  ENDDO
  

  DEALLOCATE(delta, posDelta, negDelta, pos_inf, pos_sup, interp1, interp2, interp3, interp4)
  
END SUBROUTINE interpolazione
