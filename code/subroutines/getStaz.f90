SUBROUTINE getStaz(nPts, pts, interp, period, e_period, debug)
  IMPLICIT NONE
  !!------------------------VARIABILI DI ENTRATA-----------------------!!
  LOGICAL, intent(in) :: debug
  INTEGER, intent(in) :: nPts
  REAL*8, DIMENSION(:), intent(in) :: pts(nPts), interp(nPts)
  !!-------------------------------------------------------------------!!

  !!------------------------VARIABILI DI LAVORO------------------------!!
  !!-------------------------------------------------------------------!!

  !!-------------------------VARIABILI DI USCITA-----------------------!!
  REAL*8, intent(out) :: period, e_period
  !!-------------------------------------------------------------------!!

  !!-----------------------------ALLOCATE------------------------------!!
  !!-------------------------------------------------------------------!!

  !!-----------------------------AZZERAMENTO---------------------------!!
  !!-------------------------------------------------------------------!!

END SUBROUTINE getStaz
