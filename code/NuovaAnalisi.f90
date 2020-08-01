PROGRAM analisiCefeidi
  IMPLICIT NONE
  !!-----------------VARIABILI DI CICLO------------------!!
  INTEGER :: i, j
  INTEGER :: database, nFile
  
  !!-----------------------------------------------------!!
  
  !!----------VARIABILI DI SUBROUTINE e LAVORO-----------!!
  CHARACTER*7, DIMENSION(:) :: ceph_galaxy(19)
  INTEGER :: nCeph
  REAL*8, ALLOCATABLE :: magnCat(:), perCat(:)
  INTEGER :: nObs
  REAL*8, ALLOCATABLE :: JDR(:), magR(:)
  REAL*8, ALLOCATABLE :: inc(:,:), tn(:)
  REAL*8, ALLOCATABLE :: try(:)
  REAL*8, DIMENSION(:,:) :: hubConst(19,2)
  INTEGER :: nPts
  REAL*8, ALLOCATABLE :: pts(:)
  REAL*8 :: acc
  INTEGER :: salgo, scendo, bin
  CHARACTER*7, DIMENSION(:) :: names(19)
  REAL*8, DIMENSION(:) :: vel(19), e_vel(19)
  
  !!-----------------------------------------------------!!
  
  !!-----------------VARIABILI DI TEST-------------------!!
  LOGICAL :: debug
  LOGICAL :: test
  REAL*8 :: sum, pes
  !!-----------------------------------------------------!!

  !!------------------VARIABILI DI USCITA----------------!!
  REAL*8 :: c1, c2
  REAL*8, ALLOCATABLE :: JD(:), mag(:)
  REAL*8, ALLOCATABLE :: der(:), derC(:)
  REAL*8, ALLOCATABLE :: interp(:)
  REAL*8 :: period, e_period
  REAL*8 :: mMag, e_mMag
  REAL*8 :: dist, e_dist
  REAL*8 :: h0, e_h0
  REAL*8 :: hC, e_hC
  REAL*8 :: omegaM, omegaL

  !!---------------------ALLOCAZIONE---------------------!!
  
  !!-----------------------------------------------------!!
  

  debug = .FALSE.


  !!-----------------------GET_DATA----------------------!!
  CALL get_data(ceph_galaxy, names, vel, e_vel, nCeph, debug)
  !!-----------------------------------------------------!!

  !!-----------------------ALLOCATE----------------------!!
  ALLOCATE(magnCat(nCeph), perCat(nCeph))
  !!-----------------------------------------------------!!
  
  !!-----------------------PARAM-------------------------!!
  CALL param(nCeph, magnCat, perCat, c1, c2, debug)
  PRINT*, "PARAMETRI DI FIT M-P CALCOLATI"
  IF(debug .EQV. .TRUE.) THEN
     PRINT*, c1, c2
  ENDIF
  PRINT*, c1, c2
  !!-----------------------------------------------------!!

  !!-------------------------DEALLOCATE------------------!!
  DEALLOCATE(magnCat, perCat)
  !!-----------------------------------------------------!!

  !!-------------INIZIO CICLO PRINCIPALE-----------------!!

  
 
  analisi: DO database=1,19!3!19
     nFile = 49+database
     PRINT*, ""
     PRINT*, ""
     PRINT*, "ANALIZZO ", database, ceph_galaxy(database)

     !!----------------------COUNTSINGLEDATA-----------------!!
     CALL countSingleData(nFile, nObs)
     IF(debug .EQV. .TRUE.) THEN
        PRINT*, nObs, " osservazioni."
     ENDIF
     !!------------------------------------------------------!!

     !!------------------------ALLOCATE----------------------!!
     ALLOCATE(JDR(nObs), magR(nObs))
     !!------------------------------------------------------!!
     
     !!-----------------------SINGLEDATA---------------------!!
     CALL singleData(nFile, nObs, JDR, magR)
     IF(debug .EQV. .TRUE.) THEN
        PRINT*, "MATRICE NON CORRETTA"
        DO i=1,nObs
           PRINT*, "JDR=", JDR(i), "magR=", magR(i)
        ENDDO
     ENDIF
     !!------------------------------------------------------!!

     !!-------------------------ALLOCATE---------------------!!
     ALLOCATE(JD(nObs), mag(nObs))
     !!------------------------------------------------------!!

     !!--------------------------SORTING---------------------!!
     CALL sorting(nObs, JDR, magR, JD, mag)
     IF(debug .EQV. .TRUE.) THEN
        PRINT*, "MATRICE CORRETTA:"
        DO i=1,nObs
           PRINT*, "JD=", JD(i), "mag=", mag(i)
        ENDDO
     ELSE
        PRINT*, "OSSERVAZIONI IMPORTATE"
     ENDIF
     !!------------------------------------------------------!!

     !!-------------------------DEALLOCATE-------------------!!
     DEALLOCATE(JDR, magR)
     !!------------------------------------------------------!!

     !!-------------------------ALLOCATE---------------------!!
     ALLOCATE(inc(nObs-2, nObs-2), tn(nObs-2))
     !!------------------------------------------------------!!
     !!-----------------------AZZERAMENTO--------------------!!
     inc = 0.
     tn = 0.
     !!------------------------------------------------------!!
     
     !!-------------------------SPLINE-----------------------!!
     CALL sistSpline(nObs, JD, mag, inc, tn, debug)
     IF(debug .EQV. .TRUE.) THEN
        PRINT*, "MATRICE INCOGNITE:"
        DO i=1,nObs-2
           PRINT*, (inc(i,j), j=1,nObs-2)
        ENDDO
        PRINT*, ""
        PRINT*, "MATRICE TERMINI NOTI:"
        DO i=1,nObs-2
           PRINT*, tn(i)
        ENDDO
     ELSE
        PRINT*, "SISTEMA SPLINE COSTRUITO"
     ENDIF
     !!------------------------------------------------------!!

     !!------------------------ALLOCATE----------------------!!
     ALLOCATE(der(nObs-2), derC(nObs))
     !!------------------------------------------------------!!

     !!----------------------AZZERAMENTO---------------------!!
     der = 0.
     !!------------------------------------------------------!!

     !!--------------------SOLUZIONE SISTEMA-----------------!!
     CALL soluzSist(inc, tn, nObs, der, derC, debug)
     IF(debug .EQV. .TRUE.) THEN
        PRINT*, "SOLUZIONI DEL SISTEMA LINEARE"
        DO i=1,nObs-2
           PRINT*, der(i)
        ENDDO
        ALLOCATE(try(nObs-2)) !!!--OCCHIO ALL'ALLOCAZIONE FUORI SEZIONE--!!!
        PRINT*, "VERIFICO RISULTATI SISTEMA LINEARE"
        DO i=1,nObs-2
           try(i) = 0.
           DO j=1,nObs-2
              try(i) = try(i) + inc(i,j)*der(j)
           ENDDO
        ENDDO
        DO i=1,nObs-2
           IF(try(i) == tn(i)) THEN
              PRINT*, "OK", try(i), tn(i)
           ELSE
              PRINT*, "CHECK", try(i), tn(i)
           ENDIF
        ENDDO
     ELSE
        PRINT*,"SISTEMA SPLINE RISOLTO"
     ENDIF
     !!---------------------------------------------------------------!!
     
     !!--------------------------DEALLOCATE---------------------------!!
     DEALLOCATE(inc, tn)
     IF(debug .EQV. .TRUE.) THEN
        DEALLOCATE(try)
     ENDIF
     !!---------------------------------------------------------------!!
     

     !!--------------------------INTERAZIONE--------------------------!!
     CALL interazione1(nObs, JD, nPts, acc, debug)
     !!---------------------------------------------------------------!!

     !!---------------------------ALLOCATE----------------------------!!
     ALLOCATE(pts(nPts), interp(nPts))
     !!---------------------------------------------------------------!!

     !!------------------------INTERPOLAZIONE-------------------------!!
     CALL interpolazione(nObs, nPts, acc, JD, mag, derC, pts, interp, debug)
     DO i=1,nPts
        WRITE(database+90,*) pts(i), interp(i)
     ENDDO
     
     IF(debug .EQV. .TRUE.) THEN
        PRINT*, "INTERPOLAZIONE COMPLETATA"
        PRINT*, "x", "y"
        DO i=1,nPts
           PRINT*, pts(i), interp(i)
        ENDDO
     ELSE
        PRINT*, "INTERPOLAZIONE COMPLETATA"
     ENDIF
     !!--------------------------------------------------------------!!
     IF(debug .eqv. .true.) THEN !TEMP
     CALL sleep(2)
     salgo = 0
     scendo = 0
     DO i=2,nObs
        IF(mag(i-1) > mag(i)) THEN
           IF(test .eqv. .false.) THEN
              PRINT*, salgo
              salgo = 0
              CALL sleep(1)
           ENDIF
           test = .TRUE.
           PRINT*, "SCENDO", mag(i)
           scendo = scendo + 1
        ELSEIF(mag(i-1) < mag(i)) THEN
           IF(test .eqv. .true.) THEN
              PRINT*, scendo
              scendo = 0
              CALL sleep(1)
           ENDIF
           test = .false.
           PRINT*, "SALGO", mag(i)
           salgo = salgo + 1
        ELSEIF(mag(i-1) == mag(i)) THEN
           bin = bin + 1
        ENDIF
     ENDDO
     PRINT*, salgo, scendo, bin
  ENDIF
  
     PRINT*, nObs, nPts, REAL(nPts)/REAL(nObs)

     CALL getPeriod(nObs, nPts, pts, mag, JD, interp, period, e_period, mMag, e_mMag, debug)

     
     CALL distance(period, e_period, c1, c2, mMag, e_mMag, dist, e_dist, debug)

     
     CALL hubble(vel, e_vel, dist, e_dist, h0, e_h0, database, debug)
     
     hubConst(database, 1) = h0
     hubConst(database, 2) = e_h0

     DEALLOCATE(JD, mag, pts, interp, der, derC)

     PRINT*, debug
     
     
  ENDDO analisi
  IF(debug .EQV. .TRUE.) THEN
     DO i=1,19
        PRINT*, hubConst(i,1), hubConst(i,2)
     ENDDO
  ENDIF
  
  sum = 0.
  pes = 0.
  DO i=1,19
     sum = sum + hubConst(i,1)/(hubConst(i,2)**2)
     pes = pes + 1/(hubConst(i,2)**2)
  ENDDO

  hC = sum/pes
  e_hC = SQRT(1/pes)
  PRINT*, hC, e_hC

  CALL cosmologia(hC, e_hC, omegaM, omegaL, debug)

  PRINT*, "ORA VERIFICO CON LA COSTANTE DI HUBBLE DERIVATA DALLA MISSIONE PLANCK (CMB)"
  hC = 67.4
  e_hC = 0.5

  CALL cosmologia(hC, e_hC, omegaM, omegaL, debug)

  PRINT*, "ORA VERIFICO CON LA COSTANTE DI HUBBLE DERIVATA DAL TEAM SH0ES (SN)"
  hC = 73.5
  e_hC = 1.4
  CALL cosmologia(hC, e_hC, omegaM, omegaL, debug)
  
  

END PROGRAM analisiCefeidi

INCLUDE 'subroutines/get_data.f90'

INCLUDE 'subroutines/param.f90'

INCLUDE 'subroutines/countSingleData.f90'

INCLUDE 'subroutines/singleData.f90'

INCLUDE 'subroutines/sorting.f90'

INCLUDE 'subroutines/sistSpline.f90'

INCLUDE 'subroutines/soluzSist.f90'

INCLUDE 'subroutines/interazione1.f90'

INCLUDE 'subroutines/interpolazione.f90'

INCLUDE 'subroutines/getPeriod.f90'

INCLUDE 'subroutines/distance.f90'

INCLUDE 'subroutines/hubble.f90'

INCLUDE 'subroutines/cosmologiaModerna.f90'
