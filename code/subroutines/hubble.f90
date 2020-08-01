SUBROUTINE hubble(vel, e_vel, dist, e_dist, h0, e_h0, database, debug)
  IMPLICIT NONE
  !!--------------------VARIABILI IN ENTRATA-----------------!!
  LOGICAL, intent(in) :: debug
  INTEGER, intent(in) :: database
  REAL*8, intent(in) :: dist, e_dist
  REAL*8, DIMENSION(:), intent(in) :: vel(19), e_vel(19)
  !!---------------------------------------------------------!!

  !!--------------------VARIABILI DI LAVORO------------------!!
  INTEGER :: i
  REAL*8 :: velVal, e_velVal
  REAL*8 :: e_vel_rel, e_dist_rel
  !!---------------------------------------------------------!!

  !!--------------------VARIABILI DI USCITA------------------!!
  REAL*8, intent(out) :: h0, e_h0
  !!---------------------------------------------------------!!

  velVal = vel(database)
  e_velVal = e_vel(database)
  
  h0 = velVal/dist
  PRINT*, "Costante di Hubble =", h0

  ! Calcolo gli errori
  e_vel_rel = (e_velVal)/(ABS(velVal))
  e_dist_rel = (e_dist)/(ABS(velVal))

 
  
  e_h0 = SQRT((e_velVal/dist)**2 + ((velVal/dist**2)*e_dist)**2)

  PRINT*, "ERRORE =", e_h0
  
  
END SUBROUTINE hubble
