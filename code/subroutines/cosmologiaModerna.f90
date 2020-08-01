SUBROUTINE cosmologia(hC, e_hC, omegaM, omegaL, debug)
  IMPLICIT NONE
  !INPUT
  REAL*8, intent(in) :: hC, e_hC
  LOGICAL, intent(in) :: debug

  !WORKING
  REAL*8, PARAMETER :: age = 13.82, e_age=0.1382
  REAL*8 :: M, L, MS, LS, Mi, Li
  REAL*8 :: integrale, fakt = (3.086E19/3.171E16)
  REAL*8, DIMENSION(:,:) :: density(1000000,2)
  REAL*8, ALLOCATABLE :: cosmos(:,:), aperto(:,:), chiuso(:,:)
  INTEGER :: nOmega, i, j, k, nDim, ch, ap
  REAL*8, EXTERNAL :: f

  REAL*8 :: minimo1, minimo2, minimo3, minimo1ch, minimo2ch, minimo3ch
  INTEGER :: pos_min1, pos_min2, pos_min3, pos_min1ch, pos_min2ch, pos_min3ch

  !OUTPUT
  REAL*8, intent(out) :: omegaM, omegaL
  
     
        
     

  !PRIMA PARTE: GEOMETRIA PIATTA
  !! Qui si studia se omegaM + omegaL = 1
  nOmega = 1000
  ALLOCATE(cosmos(1:nOmega,5))
  DO i=1,nOmega
     M = DBLE(i)/DBLE(nOmega)
     L = DBLE(1) - M
     cosmos(i,1) = M
     cosmos(i,2) = L
     CALL calcolo(M,L,integrale) !Operazione per ciascuna coppia di M/L
     cosmos(i,3) = (fakt/hC)*integrale
     cosmos(i,4) = (fakt/(hC+e_hC))*integrale !AGE INF
     cosmos(i,5) = (fakt/(hC-e_hC))*integrale !AGE_SUP
  ENDDO

  !Controllo età universo
  CALL ageUniverse(cosmos,nOmega, M,L, MS, LS, Mi, Li)
  PRINT*, "VERIFICO"
  CALL calcolo(M,L, integrale)
  PRINT*, (fakt/hC)*integrale

  DEALLOCATE(cosmos)


  k=0
  DO i=1,1000
     DO j=1,1000
        k= k +1
        density(k,1) = DBLE(i)/DBLE(1000)
        density(k,2) = DBLE(j)/DBLE(1000)
     ENDDO
  ENDDO
  ap=0
  ch=0

  DO k=1,1000000
     IF((density(k,1) + density(k,2)) == DBLE(1)) THEN
        CYCLE
     ELSEIF((density(k,1) + density(k,2)) .GE. DBLE(1)) THEN
        !PRINT*, "CHIUSO"
        ch = ch + 1
     ELSEIF((density(k,1) + density(k,2)) .LE. DBLE(1)) THEN
        !PRINT*, "APERTO"
        ap = ap + 1
     ENDIF    
  ENDDO

  PRINT*, "ORA ALLOCO"
  ALLOCATE(aperto(ap, 5), chiuso(ch,5))
  minimo1 = 1000.
  minimo2 = 1000.
  minimo3 = 1000.
  minimo1ch = 1000.
  minimo2ch = 1000.
  minimo3ch = 1000.
  j=1
  i=1
  DO k=1,1000000
     IF((density(k,1) + density(k,2)) .LE. DBLE(1)) THEN
        aperto(j,1) = density(k,1)
        aperto(j,2) = density(k,2)
        CALL calcolo(aperto(j,1), aperto(j,2), integrale) !Operazione per ciacuna coppia di M-L
        aperto(j,3) = (fakt/hC)*integrale
        aperto(j,4) = (fakt/(hC+e_hC))*integrale !AGE_INF
        aperto(j,5) = (fakt/(hC-e_hC))*integrale !AGE_SUP
        IF(ABS(aperto(j,3) - age) .LE. e_age) THEN
           IF(ABS(aperto(j,3) - age) < minimo1) THEN
              minimo1 = ABS(aperto(j,3) - age)
              PRINT*, "AP", minimo1
              pos_min1 = j
           ENDIF
        ENDIF
        IF(ABS(aperto(j,4) - age) .LE. e_age) THEN
           IF(ABS(aperto(j,4) - age) < minimo2) THEN
              minimo2 = ABS(aperto(j,4) - age)
              pos_min2 = j
           ENDIF
        ENDIF
        IF(ABS(aperto(j,5) - age) .LE. e_age) THEN
           IF(ABS(aperto(j,5) - age) < minimo3) THEN
              minimo3 = ABS(aperto(j,5) - age)
              pos_min3 = j
           ENDIF
        ENDIF
        j= j + 1
     ELSEIF((density(k,1) + density(k,2)) .GE. DBLE(1)) THEN
        chiuso(i,1) = density(k,1)
        chiuso(i,2) = density(k,2)
        CALL calcolo(chiuso(i,1), chiuso(i,2), integrale) !Operazione per ciascuna coppia di M-L
        chiuso(i,3) = (fakt/hC)*integrale
        chiuso(i,4) = (fakt/(hC+e_hC))*integrale
        chiuso(i,5) = (fakt/(hC-e_hC))*integrale
        IF(ABS(chiuso(i,3) - age) .LE. e_age) THEN
           IF(ABS(chiuso(i,3)- age) < minimo1ch) THEN
              minimo1ch = ABS(chiuso(i,3) - age)
              PRINT*, "CH", minimo1ch
              pos_min1ch = i
           ENDIF
        ENDIF
        IF(ABS(chiuso(i,4) - age) .LE. e_age) THEN
           IF(ABS(chiuso(i,4) - age) < minimo2ch) THEN
              minimo2ch = ABS(chiuso(i,4) - age)
              pos_min2ch = i
           ENDIF
        ENDIF
        IF(ABS(chiuso(i,5) - age) .LE. e_age) THEN
           IF(ABS(chiuso(i,5) - age) < minimo3ch) THEN
              minimo3ch = ABS(chiuso(i,5) - age)
              pos_min3ch = i
           ENDIF
        ENDIF
        i = i + 1
     ENDIF
  ENDDO

  PRINT*, "UNIVERSO APERTO"
  PRINT*, "Medi", aperto(pos_min1,1), aperto(pos_min1,2)
  PRINT*, "Minimi", aperto(pos_min2,1), aperto(pos_min2,2)
  PRINT*, "Massimi", aperto(pos_min3,1), aperto(pos_min3,2)
  PRINT*, "VERIFICO"
  CALL calcolo(aperto(pos_min1,1),aperto(pos_min1,2),integrale)
  PRINT*, "Medio", (fakt/hC)*integrale
  CALL calcolo(aperto(pos_min2,1),aperto(pos_min2,2),integrale)
  PRINT*, "Minimo", (fakt/(hC+e_hC))*integrale
  CALL calcolo(aperto(pos_min3,1),aperto(pos_min3,2), integrale)
  PRINT*, "Massimo", (fakt/(hC-e_hC))*integrale

  PRINT*, "UNIVERSO CHIUSO"
  PRINT*, "Medi", chiuso(pos_min1ch,1), chiuso(pos_min1ch,2)
  PRINT*, "Minimi", chiuso(pos_min2ch,1), chiuso(pos_min2ch,2)
  PRINT*, "Massimi", chiuso(pos_min3ch,1), chiuso(pos_min3ch,2)
  PRINT*, "VERIFICO"
  CALL calcolo(chiuso(pos_min1ch,1),chiuso(pos_min1ch,2),integrale)
  PRINT*, "Medio", (fakt/hC)*integrale
  CALL calcolo(chiuso(pos_min2ch,1),chiuso(pos_min2ch,2),integrale)
  PRINT*, "Minimo", (fakt/(hC+e_hC))*integrale
  CALL calcolo(chiuso(pos_min3ch,1),chiuso(pos_min3ch,2),integrale)
  PRINT*, "Massimo", (fakt/(hC-e_hC))*integrale
  PRINT*, "FINITO"

  DEALLOCATE(aperto, chiuso)
END SUBROUTINE cosmologia

SUBROUTINE ageUniverse(cosmos, nOmega, M, L, MS, LS, Mi, Li)
  IMPLICIT NONE
  !INPUT
  REAL*8, DIMENSION(:,:), intent(in) :: cosmos(nOmega,5)
  !WORKING
  REAL*8, PARAMETER :: age = 13.82, e_age = 0.1382
  REAL*8, ALLOCATABLE :: err(:,:)
  REAL*8 :: minimo1, minimo2, minimo3
  INTEGER :: pos_min1, pos_min2, pos_min3
  INTEGER :: nOmega, i
  !OUTPUT
  REAL*8, intent(out) :: M, L, MS, LS, Mi, Li
  ALLOCATE(err(nOmega,3))
  err = HUGE(minimo1)
  DO i=1,nOmega
     IF(ABS(cosmos(i,3) - age) .LE. e_age) THEN
        err(i,1) = ABS(cosmos(i,3) - age)
        !Questo è un possibile valore MEDIO
     ENDIF
     IF(ABS(cosmos(i,4) - age) .LE. e_age) THEN
        err(i,2) = ABS(cosmos(i,4) - age)
        !Questo è un possibile valore MINIMO
     ENDIF
     IF(ABS(cosmos(i,5) - age) .LE. e_age) THEN
        err(i,3) = ABS(cosmos(i,5) - age)
        !Questo è un possibile valore MASSIMO
     ENDIF
  ENDDO

  !Ho costruito nell'array err() tutte le differenze tra il valore 'vero' dell'età dell'universo e quello calcolato da me.
  !Ora scelgo, per ciascuno, il valore minimo e riporto i lambda relativi

  minimo1 = 1000.
  minimo2 = 1000.
  minimo3 = 1000.
  
  DO i=1,nOmega
     IF(err(i,1) .LE. minimo1) THEN
        minimo1 = err(i,1)
        pos_min1 = i
     ENDIF
     IF(err(i,2) .LE. minimo2) THEN
        minimo2 = err(i,2)
        pos_min2 = i
     ENDIF
     IF(err(i,3) .LE. minimo3) THEN
        minimo3 = err(i,3)
        pos_min3 = i
     ENDIF
  ENDDO

  !A questo punto, ho in pos_min* gli indici di cosmos (in cui trovo i valori minimi)

  PRINT*, "Valori per costante di Hubble media"
  PRINT*, pos_min1, pos_min2, pos_min3
  PRINT*, cosmos(pos_min1, 1), cosmos(pos_min1, 2), pos_min1
  PRINT*, ""
  PRINT*, "Valori per costante di Hubble massima"
  PRINT*, cosmos(pos_min2,1), cosmos(pos_min2, 2), pos_min2
  PRINT*, ""
  PRINT*, "Valori per costante di Hubble minima"
  PRINT*, cosmos(pos_min3, 1), cosmos(pos_min3, 2), pos_min3
  M = cosmos(pos_min1, 1)
  L = cosmos(pos_min1, 2)
  
DEALLOCATE(err)
END SUBROUTINE ageUniverse

SUBROUTINE calcolo(M, L, integrale)
  IMPLICIT NONE
  REAL*8, intent(in) :: M, L
  REAL*8 :: risultato
  REAL*8, intent(out) :: integrale
  INTEGER :: i
  REAL*8, EXTERNAL :: f

  integrale = DBLE(0)
  
  DO i=1,1100
     CALL quadgauss(f, DBLE(i-1), DBLE(i), 4, M, L, risultato)
     integrale = integrale + risultato
  ENDDO

END SUBROUTINE calcolo

SUBROUTINE quadgauss(f, a, b, npunti, M, L, risultato)
  IMPLICIT NONE
  REAL*8, EXTERNAL :: f
  REAL*8, intent(in) :: a, b, M, L
  INTEGER, intent(in) :: npunti
  REAL*8, intent(out) :: risultato

  REAL*8 :: summa
  REAL*8, ALLOCATABLE :: c(:), x(:), xd(:)
  INTEGER :: i

  ALLOCATE(x(0:npunti-1), c(0:npunti-1), xd(0:npunti-1))

  SELECT CASE(npunti)
  CASE(2)
     c(0) = 1.
     c(1) = 1.
     xd(0) = -1./SQRT(3.)
     xd(1) = 1./SQRT(3.)
  CASE(3)
     c(0) = 5./9.
     c(1) = 8./9.
     c(2) = 5./9.
     xd(0) = -0.774596669
     xd(1) = 0.
     xd(2) = -xd(0)
  CASE(4)
     xd(0) = SQRT(3./7. - 2./7.*SQRT(6./5.))
     c(0) = (18.+SQRT(30.))/36.
     c(1) = c(0)
     xd(1) = -xd(0)
     c(2) = (18.-SQRT(30.))/36.
     c(3) = c(2)
     xd(2) = SQRT(3./7. + 2./7.*SQRT(6./5.))
     xd(3) = -xd(2)
  CASE DEFAULT
     PRINT*, "Caso non ancora implementato"
     STOP
  END SELECT

  DO i=0,npunti-1
     x(i) = 0.5*((b+a)+(b-a)*xd(i))
  ENDDO

  summa = 0.
  DO i=0,npunti-1
     summa = summa + c(i)*f(x(i), M, L)
  ENDDO
  risultato = summa*(b-a)*0.5

  DEALLOCATE(xd,x,c)
END SUBROUTINE quadgauss

REAL*8 FUNCTION f(z, M, L)
  IMPLICIT NONE
  REAL*8 :: z, M, L
  f = 1./((1+z)*SQRT(M*(1+z)**3 + (1-M-L)*(1+z)**2 + L))
END FUNCTION f
  
