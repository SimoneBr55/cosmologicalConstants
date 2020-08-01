SUBROUTINE distance(period, e_period, c1, c2, mMag, e_mMag, dist, e_dist, debug)
  IMPLICIT NONE
  !!-------------------VARIABILI IN ENTRATA------------------!!
  LOGICAL, intent(in) :: debug
  REAL*8, intent(in) :: period, e_period
  REAL*8, intent(in) :: c1, c2
  REAL*8, intent(in) :: mMag, e_mMag
  !!---------------------------------------------------------!!

  !!------------------VARIABILI DI LAVORO--------------------!!
  REAL*8 :: logP, e_logP
  REAL*8 :: mA, e_mA
  REAL*8 :: distMod, e_distMod
  !!---------------------------------------------------------!!

  !!-------------------VARIABILI DI USCITA-------------------!!
  REAL*8, intent(out) :: dist, e_dist
  !!---------------------------------------------------------!!

  logP = LOG10(period)
  mA = c1*logP + c2

  e_logP = LOG10(e_period)
  e_mA = (c1/period)*e_logP

  PRINT*, mA, e_mA

  distMod = mMag - mA
  e_distMod = SQRT(e_mMag**2 + e_mA**2)

  PRINT*, "Modulo di Distanza = ", distMod, e_distMod, (e_distMod/distMod)*100

  dist = 10**(0.2*distMod + 1)
  dist = dist/(10**6)
  e_dist = 4.61*EXP(0.461*distMod)*e_distMod
  e_dist = e_dist/(10**6)

  PRINT*, "Distanza = ", dist, dist*10**6
  PRINT*, "ERR = ", e_dist, e_dist*10**6, ((e_dist/dist)*100)
  
END SUBROUTINE distance
