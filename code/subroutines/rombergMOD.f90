SUBROUTINE romberg(geom, M, L, a, b, es, integrale)
  IMPLICIT NONE
  INTEGER, intent(in) :: geom
  REAL*8 :: a, b, M, L
  REAL*8, intent(out) :: integrale
  
  REAL*8, EXTERNAL:: integranda1, integranda2, integranda3
  REAL*8, intent(in) :: es
  REAL*8 :: e_a, risultato
  INTEGER::i,n,ordine, maxit, j,k
  REAL*8, ALLOCATABLE:: mat(:,:)
  
  ordine=1
  maxit=100
  
  ALLOCATE(mat(maxit,maxit))

  mat=0.
  
  e_a=1.1*es
  i=0
  DO WHILE(e_a>es.and.i<maxit)
     i=i+1
     n=2**(i-1)
     IF(geom == 0) THEN
        CALL newton(integranda1,M,L,a,b,n,3,risultato)
     ELSE
        CALL newton(integranda2,M,L, a, b, n,3, risultato)
     ENDIF
     mat(i,1)=risultato
     DO k=2,i
        j=1+i-k
        mat(j,k)=(4.**(k-1)*mat(j+1,k-1)-mat(j,k-1))/&
             (4.**(k-1)-1)
     END DO
     IF(i/=1) e_a=100.*ABS((mat(1,i)-mat(1,i-1))/mat(1,i))
     PRINT *, n, mat(1,i)
     integrale = mat(1,i)
  END DO
  
END SUBROUTINE romberg

SUBROUTINE newton(f,M,L,a,b,n,ordine,risultato)
  IMPLICIT NONE
  REAL*8:: a,b,risultato, h, area, x1, x2,xm,xm1,xm2,delta
  INTEGER::n, i, ordine
  REAL*8, EXTERNAL:: f
  REAL*8, intent(in) :: M, L
  risultato=0.d0
  h=(b-a)/n
  DO i=1,n
     x1=a+(i-1)*h
     x2=a+i*h
     SELECT CASE(ordine)
        CASE(1)
           area=h*(f(x1)+f(x2))/2.d0
        CASE(2)
           xm=(x1+x2)/2.d0
           area=h*(f(x1)+4.d0*f(xm)+f(x2))/6.d0
        CASE(3)
           delta=(x2-x1)/3.d0
           xm1=x1+delta
           xm2=x1+2.d0*delta
           area=h*(f(x1)+3.d0*f(xm1)+3.d0*f(xm2)+f(x2))/8.d0
        CASE DEFAULT
           PRINT *, "caso non previsto"
           STOP
      END SELECT      
     risultato=risultato+area
  END DO    
END SUBROUTINE newton

REAL*8 FUNCTION integranda1(x, omegaM, omegaL)
  IMPLICIT NONE
  REAL*8::x
  REAL*8 :: omegaM, omegaL
  integranda1=1/((x+1)*SQRT(omegaM*(1+x)**3 + omegaL))
END FUNCTION integranda1

REAL*8 FUNCTION integranda2(x,omegaM, omegaL)
  IMPLICIT NONE
  REAL*8 :: x
  REAL*8 :: omegaM, omegaL
  integranda2=1/((x+1)*SQRT(omegaM*(1+x)**3 + (1-omegaM-omegaL)*(1+x)**2 + omegaL))
END FUNCTION INTEGRANDA2
  
