!!! DA OTTIMIZZARE E RISCRIVERE E RENDERE PIU' ELEGANTEEEEEEE

SUBROUTINE getPeriod(nObs, nPts, pts, mag, JD, interp, period, e_period, mMag, e_mMag, debug)
  IMPLICIT NONE
  !!INPUT!!
  LOGICAL, intent(in) :: debug
  INTEGER, intent(in) :: nObs, nPts
  REAL*8, DIMENSION(:), intent(in) :: pts(nPts), interp(nPts), mag(nObs), JD(nObs)
  !!LAVORO!!
  REAL*8 :: mean, sum
  LOGICAL :: step
  INTEGER :: salita, discesa, i, j, k
  REAL*8, ALLOCATABLE :: sal(:,:), dis(:,:), JD_sal(:), JD_dis(:)
  INTEGER, ALLOCATABLE :: pos_sal(:), pos_dis(:)
  INTEGER, ALLOCATABLE :: salR(:,:), disR(:,:)
  REAL*8 :: temp, mean_sal, mean_dis, minimo_interp, err_per
  INTEGER :: minimo
  REAL*8, ALLOCATABLE :: per_sal(:), per_dis(:)
  

  !!OUTPUT!!
  REAL*8, intent(out) :: period, e_period
  REAL*8, intent(out) :: mMag, e_mMag
  
  sum = 0.

  DO i=1,nObs
     sum = sum + mag(i)
  ENDDO
  mean = sum/nObs
  PRINT*, mean

  salita = 0
  discesa = 0

  DO i=2,nObs
     IF(step .EQV. .FALSE.) THEN
        IF(mag(i-1) < mean .and. mag(i) > mean) THEN
           salita = salita + 1
           step = .TRUE.
        ELSEIF(mag(i-1) > mean .and. mag(i) < mean) THEN
           discesa = discesa + 1
           step = .TRUE.
        ENDIF
     ELSE
        step = .FALSE.
     ENDIF
           
  ENDDO

  ALLOCATE(sal(salita, salita), dis(discesa, discesa))
  ALLOCATE(salR(salita, salita), disR(discesa, discesa))
  sal = 0.
  dis = 0.
  salR = 0.
  disR = 0.

  j=1
  k=1
  DO i=2,nObs-1
     IF(step .EQV. .FALSE.) THEN
        IF(mag(i-1) < mean .and. mag(i+1) > mean) THEN
           sal(j,1) = JD(i-1)
           sal(j,2) = JD(i+1)
           j = j+1
           step = .TRUE.
        ELSEIF(mag(i-1) > mean .and. mag(i+1) < mean) THEN
           dis(k,1) = JD(i-1)
           dis(k,2) = JD(i+1)
           k = k+1
           step = .TRUE. 
        ENDIF
     ELSE
        step = .FALSE.
     ENDIF
  ENDDO

  !!FINORA ho fatto i calcoli, senza l'interpolazione...
  !!ADESSO devo prendere quelle che sono STIME di passaggi, e cercare il vero passaggio.
  !Le stime dei passaggi sono da fare attorno agli sal(i)

  ! Analizzo i passaggi in SALITA
  i=1
  j=1
  k=1
568 CONTINUE
  
  DO WHILE(pts(i) < sal(j,1) .AND. i < nPts)
     i = i + 1
  ENDDO
  
  salR(j,1) = i
  k = i
 
  DO WHILE(pts(k) < sal(j,2) .AND. k < nPts)
     k = k + 1
  ENDDO

  salR(j,2) = k
  IF(k .GE. nPts) THEN
     GOTO 570
  ELSE
     IF(debug .EQV. .TRUE.) PRINT*, "CONTINUO" !
     i = k+1
     IF(j >= salita) THEN
        GOTO 570
     ELSE
        j = j + 1
        GOTO 568
     ENDIF
  ENDIF

570 CONTINUE

  ! Analizzo i passaggi in DISCESA

  i=1
  j=1
  k=1
569 CONTINUE
  DO WHILE(pts(i) < dis(j,1) .AND. i < nPts)
     i = i + 1
  ENDDO
  disR(j,1) = i
  k = i
  DO WHILE(pts(k) < dis(j,2) .AND. k < nPts)
     k = k + 1
  ENDDO
  disR(j,2) = k
  IF(k .GE. nPts) THEN
     GOTO 571
  ELSE
     i = k + 1
     IF(j >= discesa) THEN
        GOTO 571
     ELSE
        j = j + 1
        GOTO 569
     ENDIF
  ENDIF

571 CONTINUE
  
  IF(debug .EQV. .TRUE.) THEN
     DO i=1,salita
        PRINT*, sal(i,1), sal(i,2), "/", salR(i,1), salR(i,2)
     ENDDO
     DO i=1, discesa
        PRINT*, dis(i,1), dis(i,2), "/", disR(i,1), disR(i,2)
     ENDDO
  ENDIF

  sum = 0.
  mean = 0.

  DO i=1,nPts
     sum = sum + interp(i)
  ENDDO
  mean = sum/nPts
  
  ALLOCATE(JD_sal(salita), JD_dis(discesa))
  ALLOCATE(pos_sal(salita), pos_dis(discesa))
  ALLOCATE(per_sal(salita-1), per_dis(discesa-1))

 
  DO j=1,salita
     minimo_interp = interp(salR(j,1))
     minimo = pts(salR(j,1))
     DO i=salR(j,1), salR(j,2)
        IF(ABS(interp(i)-mean) < minimo_interp) THEN
           minimo_interp = ABS(interp(i)-mean)
           minimo = i
        ENDIF
     ENDDO
     pos_sal(j) = minimo
     IF(debug .EQV. .TRUE.) PRINT*, "CONFRONTO MINIMO, sal1, sal2" !!!
  ENDDO
  IF(debug .EQV. .TRUE.) PRINT*, "NUMERO DISCESE = ", discesa !!! 
  DO j=1,discesa
     minimo_interp = interp(disR(j,1))
     minimo = pts(disR(j,1))
     DO i=disR(j,1),disR(j,2)
        IF(ABS(interp(i)-mean) < minimo_interp) THEN
           minimo_interp = ABS(interp(i) - mean)
           minimo = i
        ENDIF
     ENDDO
     pos_dis(j) = minimo
     IF(debug .EQV. .TRUE.) PRINT*, "CONFRONTO MINIMO, dis1, dis2" !!!
     IF(debug .EQV. .TRUE.) PRINT*, minimo, disR(j,1), disR(j,2) !!!
  ENDDO

  IF(debug .EQV. .TRUE.) PRINT*, "SALITE" !!!
  DO i=2, salita
     per_sal(i-1) = pts(pos_sal(i)) - pts(pos_sal(i-1))
     IF(debug .EQV. .TRUE.) PRINT*, per_sal(i-1)  !!!
  ENDDO
  IF(debug .EQV. .TRUE.) print*, "DISCESE" !!!
  DO i=2, discesa
     per_dis(i-1) = pts(pos_dis(i)) - pts(pos_dis(i-1))
     IF(debug .EQV. .TRUE.) PRINT*, per_dis(i-1) !!!
  ENDDO

  ! CALCOLO MEDIA

  sum = 0.
  mean = 0.

  DO i=1,salita-1
     sum = sum + per_sal(i)
  ENDDO
  DO i=1,discesa-1
     sum = sum + per_dis(i)
  ENDDO
  mean = sum/(salita + discesa - 2)

  sum = 0.
  DO i=1,salita-1
     sum = sum + (per_sal(i) - mean)**2
  ENDDO
  DO i=1,discesa-1
     sum = sum + (per_dis(i) - mean)**2
  ENDDO
  err_per = SQRT(sum/(salita + discesa -2))


  period = mean
  e_period = err_per

   sum=0.
   DO i=pos_sal(1), pos_sal(salita)
      sum = sum + interp(i)
   ENDDO
   mMag = sum/(pos_sal(salita) - pos_sal(1))
   !PRINT*, mMag

   sum=0.
   DO i=pos_dis(1), pos_dis(discesa)
      sum = sum + interp(i)
   ENDDO
   sum = sum/(pos_dis(discesa) - pos_dis(1))
   mMag = (mMag + sum)/2
   
   sum = 0.
   DO i=pos_dis(1), pos_dis(discesa)
      sum = sum + (interp(i) - mMag)**2
   ENDDO
   e_mMag = sum/(pos_dis(discesa) - pos_dis(1))
   !PRINT*, e_mMag
   sum = 0.
   DO i=pos_sal(1), pos_sal(salita)
      sum = sum + (interp(i) - mMag)**2
   ENDDO
   sum = sum/(pos_sal(salita) - pos_dis(1))
   e_mMag = (e_mMag + sum)/2
   

  PRINT*, "MEDIA_MAGNITUDINE = ", mMag
  PRINT*, "DS = ", e_mMag, (e_mMag/mMag)*100

  PRINT*, "PERIODO = ", period
  PRINT*, "ERR = ", e_period, (e_period/period)*100
  
  
  DEALLOCATE(sal, dis)
  DEALLOCATE(salR, disR)
  DEALLOCATE(pos_sal, pos_dis)
  DEALLOCATE(JD_sal, JD_dis)
  DEALLOCATE(per_sal, per_dis)

    

END SUBROUTINE getPeriod
