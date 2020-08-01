SUBROUTINE get_data(ceph_galaxy, names, vel, e_vel, nCeph, debug)
  IMPLICIT NONE
  !!----------------------INPUT-------------------------!!
  LOGICAL, intent(in) :: debug
  !!----------------------------------------------------!!

  !!-----------------------WORK-------------------------!!
  INTEGER :: i, nGal
  !!----------------------------------------------------!!

  !!-----------------------OUTPUT-----------------------!!
  CHARACTER*7, DIMENSION(:), intent(out) :: ceph_galaxy(19) !OK!
  INTEGER, intent(out) :: nCeph !OK!
  CHARACTER*7, DIMENSION(:), intent(out) :: names(19)
  REAL*8, DIMENSION(:), intent(out) :: vel(19), e_vel(19)
  !!----------------------------------------------------!!
  

  !AZZERO VARIABILI
  ceph_galaxy = '' !AZZERATO
  
  OPEN(50,file='data/ceph_NGC0925.txt')
  ceph_galaxy(1) = 'NGC0925'
  OPEN(91,file='data/01_cdlNGC0925.txt')
  
  OPEN(51,file='data/ceph_NGC1326.txt')
  ceph_galaxy(2) = 'NGC1326'
  OPEN(92,file='data/02_cdlNGC1326.txt')
  
  OPEN(52,file='data/ceph_NGC1365.txt')
  ceph_galaxy(3) = 'NGC1365'
  OPEN(93,file='data/03_cdlNGC1365.txt')
  
  OPEN(53,file='data/ceph_NGC1425.txt')
  ceph_galaxy(4) = 'NGC1425'
  OPEN(94,file='data/04_cdlNGC1425.txt')
  
  OPEN(56,file='data/ceph_NGC2090.txt')
  ceph_galaxy(7) = 'NGC2090'
  OPEN(97,file='data/07_cdlNGC2090.txt')
  
  OPEN(54,file='data/ceph_NGC2403.txt')
  ceph_galaxy(5) = 'NGC2403'
  OPEN(95,file='data/05_cdlNGC2403.txt')
  
  OPEN(55,file='data/ceph_NGC2541.txt')
  ceph_galaxy(6) = 'NGC2541'
  OPEN(96,file='data/06_cdlNGC2541.txt')
  
  OPEN(57,file='data/ceph_NGC3198.txt')
  ceph_galaxy(8) = 'NGC3198'
  OPEN(98,file='data/08_cdlNGC3198.txt')
  
  OPEN(58,file='data/ceph_NGC3351.txt')
  ceph_galaxy(9) = 'NGC3351'
  OPEN(99,file='data/09_cdlNGC3351.txt')
  
  OPEN(59,file='data/ceph_NGC3368.txt')
  ceph_galaxy(10) = 'NGC3368'
  OPEN(100,file='data/10_cdlNGC3368.txt')
  
  OPEN(60,file='data/ceph_NGC3621.txt')
  ceph_galaxy(11) = 'NGC3521'
  OPEN(101,file='data/11_cdlNGC3621.txt')
  
  OPEN(61,file='data/ceph_NGC4321.txt')
  ceph_galaxy(12) = 'NGC4321'
  OPEN(102,file='data/12_cdlNGC4321.txt')
  
  OPEN(62,file='data/ceph_NGC4496.txt')
  ceph_galaxy(13) = 'NGC4496'
  OPEN(103,file='data/13_cdlNGC4496.txt')
  
  OPEN(64,file='data/ceph_NGC4535.txt')
  ceph_galaxy(15) = 'NGC4535'
  OPEN(105,file='data/15_cdlNGC4535.txt')
  
  OPEN(65,file='data/ceph_NGC4536.txt')
  ceph_galaxy(16) = 'NGC4536'
  OPEN(106,file='data/16_cdlNGC4536.txt')
  
  OPEN(63,file='data/ceph_NGC4548.txt')
  ceph_galaxy(14) = 'NGC4548'
  OPEN(104,file='data/14_cdlNGC4548.txt')
  
  OPEN(66,file='data/ceph_NGC4639.txt')
  ceph_galaxy(17) = 'NGC4639'
  OPEN(107,file='data/17_cdlNGC4639.txt')
  
  OPEN(67,file='data/ceph_NGC4725.txt')
  ceph_galaxy(18) = 'NGC4725'
  OPEN(108,file='data/18_cdlNGC4725.txt')
  
  OPEN(68,file='data/ceph_NGC7331.txt')
  ceph_galaxy(19) = 'NGC7331'
  OPEN(109,file='data/19_cdlNGC7331.txt')

  OPEN(120, file="data/ceph_catalog.txt")
  REWIND 120
  nCeph = 0 !AZZERATO
  DO
     READ(120,*,END=998)
     nCeph = nCeph + 1
  ENDDO

  998 CONTINUE

  !! GET VELOCITIES !!

  OPEN(111, file='data/gal_vel.txt')

  REWIND 111
  nGal = 0
  DO
     READ(111,*,END=555)
     nGal = nGal + 1
  ENDDO

555 CONTINUE
  IF(debug .EQV. .TRUE.) THEN
     IF(nGal /= 19) PRINT*, "C'è qualche problema..."
  ENDIF

  REWIND 111
  DO i=1,19
     READ(111,*) names(i), vel(i), e_vel(i)
  ENDDO
  CLOSE(111)

  IF(debug .EQV. .TRUE.) THEN
     PRINT*, "Queste sono i dati sulla velocità di recessione... serviranno più avanti nello script."
     DO i=1,19
        PRINT*, names(i), vel(i), e_vel(i)
     ENDDO
  ENDIF
  

END SUBROUTINE get_data

