MODULE pippo
  IMPLICIT NONE
  REAL*8:: r,theta,phi
  SAVE
END MODULE pippo
  
PROGRAM multi_integral
  USE pippo
  IMPLICIT NONE
  REAL*8::r0, r1, risultato, exact
  REAL*8, EXTERNAL:: integranda
  REAL*8, PARAMETER:: pi=ACOS(-1.)
  
  r0=0.
  r1=1.
  CALL quad(integranda,r0,r1,risultato)
  exact=4.*pi/3.
  PRINT *, risultato, exact, 100.*ABS((risultato-exact)/exact)
  
END PROGRAM multi_integral

REAL*8 FUNCTION integranda(rr)
  USE pippo
  IMPLICIT NONE
  REAL*8::rr, theta0,theta1, risultato
  REAL*8, EXTERNAL:: integranda2
  r=rr
  theta0=-0.5*ACOS(-1.)
  theta1=-theta0
  CALL quad2(integranda2,theta0,theta1,risultato)
  integranda=risultato
END FUNCTION integranda

REAL*8 FUNCTION integranda2(ttheta)
  USE pippo
  IMPLICIT NONE
  REAL*8::ttheta, phi0,phi1, risultato
  REAL*8, EXTERNAL:: integranda3
  theta=ttheta
  phi0=0.
  phi1=2.*ACOS(-1.)
  CALL quad3(integranda3,phi0,phi1,risultato)
  integranda2=risultato
END FUNCTION integranda2

REAL*8 FUNCTION integranda3(pphi)
  USE pippo
  IMPLICIT NONE
  REAL*8:: pphi
  phi=pphi
  integranda3=(r**2)*COS(theta)
END FUNCTION integranda3

SUBROUTINE quad(f,a,b,integrale)
  IMPLICIT NONE
  REAL*8, external:: f
  REAL*8:: a,b, x0, x1, x2, x3, x4, x5, c0,c1,c2, c3, c4, &
       c5, xx0, xx1, xx2, xx3, xx4, xx5,&
       ris_esatto,err, integrale, dx

  x0=-0.932469514
  x1=-0.661209386
  x2=-0.238619186
  x3=-x2
  x4=-x1
  x5=-x0
  c0=0.171324492
  c1=0.360761573
  c2=0.467913935
  c3=c2
  c4=c1
  c5=c0

  xx0=0.5*(a+b)+0.5*(b-a)*(x0)
  xx1=0.5*(a+b)+0.5*(b-a)*(x1)
  xx2=0.5*(a+b)+0.5*(b-a)*(x2)
  xx3=0.5*(a+b)+0.5*(b-a)*(x3)
  xx4=0.5*(a+b)+0.5*(b-a)*(x4)
  xx5=0.5*(a+b)+0.5*(b-a)*(x5)
  
  dx=0.5*(b-a)
  integrale=dx*(c0*f(xx0)+c1*f(xx1)+c2*f(xx2)+&
       c3*f(xx3)+c4*f(xx4)+c5*f(xx5))
END SUBROUTINE quad

SUBROUTINE quad3(f,a,b,integrale)
  IMPLICIT NONE
  REAL*8, external:: f
  REAL*8:: a,b, x0, x1, x2, x3, x4, x5, c0,c1,c2, c3, c4, &
       c5, xx0, xx1, xx2, xx3, xx4, xx5,&
       ris_esatto,err, integrale, dx

  x0=-0.932469514
  x1=-0.661209386
  x2=-0.238619186
  x3=-x2
  x4=-x1
  x5=-x0
  c0=0.171324492
  c1=0.360761573
  c2=0.467913935
  c3=c2
  c4=c1
  c5=c0

  xx0=0.5*(a+b)+0.5*(b-a)*(x0)
  xx1=0.5*(a+b)+0.5*(b-a)*(x1)
  xx2=0.5*(a+b)+0.5*(b-a)*(x2)
  xx3=0.5*(a+b)+0.5*(b-a)*(x3)
  xx4=0.5*(a+b)+0.5*(b-a)*(x4)
  xx5=0.5*(a+b)+0.5*(b-a)*(x5)
  
  dx=0.5*(b-a)
  integrale=dx*(c0*f(xx0)+c1*f(xx1)+c2*f(xx2)+&
       c3*f(xx3)+c4*f(xx4)+c5*f(xx5))
END SUBROUTINE quad3

SUBROUTINE quad2(f,a,b,integrale)
  IMPLICIT NONE
  REAL*8, external:: f
  REAL*8:: a,b, x0, x1, x2, x3, x4, x5, c0,c1,c2, c3, c4, &
       c5, xx0, xx1, xx2, xx3, xx4, xx5,&
       ris_esatto,err, integrale, dx

  x0=-0.932469514
  x1=-0.661209386
  x2=-0.238619186
  x3=-x2
  x4=-x1
  x5=-x0
  c0=0.171324492
  c1=0.360761573
  c2=0.467913935
  c3=c2
  c4=c1
  c5=c0

  xx0=0.5*(a+b)+0.5*(b-a)*(x0)
  xx1=0.5*(a+b)+0.5*(b-a)*(x1)
  xx2=0.5*(a+b)+0.5*(b-a)*(x2)
  xx3=0.5*(a+b)+0.5*(b-a)*(x3)
  xx4=0.5*(a+b)+0.5*(b-a)*(x4)
  xx5=0.5*(a+b)+0.5*(b-a)*(x5)
  
  dx=0.5*(b-a)
  integrale=dx*(c0*f(xx0)+c1*f(xx1)+c2*f(xx2)+&
       c3*f(xx3)+c4*f(xx4)+c5*f(xx5))
END SUBROUTINE quad2
