      IMPLICIT REAL(kind=8)(A-H,O-Z)
C   PROGRAM TO SOLVE THE INHOMOGENEOUS PERTURBATION EQUATION BY THE METHOD
C   OF FROBENIUS FOR AN ARBITRARY HYDROGENIC STATE.
C
      DIMENSION A(20,4),B(20,4),D(3,3),E2(3),EX2(3),DD(3),T(202)
     1   ,FAC(202),RR(202),TL(202),C(100),CL(100),X(3),VB(20,4),A3(20)
      GAM = 100000. ! Integrals over log(r) need this term. Final result should end up being independant of this. (euler mascheroni constant)
      IC = 0
      LRGL = 1 ! Total angular momentum
      IZ = 2   ! Nuclear Charge
C   CALCULATE AND STORE FACTORIALS FAC(I) = (I-1)!, AND SUMS OF RECIPROCALS RR(I).
C   LRGL = ANGULAR MOMENTUM
C   NEIG = PRINCIPAL QUANRUM NUMBER
      FAC(1) = 1
      RR(1) = 0. ! Sums of reciprocals
      XX = 1.
      YY = 1.
      DO I=2,50
      YY = YY + 1.D0/(I-1)
      RR(I) = YY
      XX = XX*(I-1)
      FAC(I) = XX
      end do
C
  100 READ(*,*) NEIG
      IC = IC + 1
C N is the principle quantum number, NEIG is the degree of excitation above the lowest state of that angular momenta. NEIG is the number of oscillations of the wavefunction.
C In other words, NEIG is the number of zeros plus 1 in the wavefunction.
      N = NEIG + LRGL
      E0 = -1.D0/(2*N**2)
      NP1 = N + 1
      NP2 = N + 2
      NN = 2*N + 4
      XLNA = DLOG(2.D0/N)
      AC = N/2.D0
      T(1) = AC
      TL(1) = AC*(-GAM - XLNA + RR(1))
C      TL(1) = 0.
      XX = AC
C   CALCULFE AND STORE THE NEEDED INTEGRALS T(I)
C   T(I) = (I-1)!*AC**I
      DO I=2,NN
      XX = XX*AC*(I-1)
      T(I) = XX
      TL(I) = XX*(-GAM - XLNA + RR(I))
C      TL(I) = 0.
      end do

      DO 9 J=1,20
      A3(J) = 0.
      C(J) = 0.
      CL(J) = 0.
      DO 9 I=1,4
      VB(J,I) = 0.
      B(J,I) = 0.
    9 A(J,I) = 0.
C   CALCULATION OF HYDROGENIC WAVEFUNCTION PSI_0 (SAME AS IN dbigdq.f90)
C   A(M,1) IS THE COEF. OF R**(M-1) IN PSI(0).
      XX= 1
      DO 611 M=1,NEIG
      A(M+LRGL,1) = 0.5*XX*2**(LRGL+2)
     1   *DSQRT(FAC(N+LRGL+1)/(FAC(N-LRGL)*FAC(2)))
     2   /(FAC(2*LRGL+2)*N*(N)**(LRGL+1))
  611 XX = -XX*2*(NEIG-M)/((2*LRGL+M+1)*M*N)
C
C   CHECK:  CALCULATE THE NORM
      SUM = 0.
      DO 8 I=1,N
      DO 8 J=1,N
    8 SUM = SUM + A(I,1)*A(J,1)*T(I+J+1)
      WRITE(*,'('' NORM ='',D15.7)') SUM
C
C   CALCULATE -(1/2)DEL**2 PSI_0
      A(1,2) = -0.5*((2-LRGL*(LRGL+1))*A(2,1) - 2.*A(1,1)/N)
      DO 11 I=2,NP1
   11 A(I,2) = -0.5*((I*(I+1)-LRGL*(LRGL+1))*A(I+1,1) - 2.*I*A(I,1)/N
     1   + A(I-1,1)/(N*N))
C
C   CALCULATE [-(1/2)DEL**2}**2 PSI_0
      DO 12 I=1,NP1
      A(I,3) = -0.5*((I*(I+1)-LRGL*(LRGL+1))*A(I+2,2) - 2.*I*A(I+1,2)/N
     1   + A(I,2)/(N*N))
   12 CONTINUE
      an = 2*dsqrt(6.d0)*2
      write(*,*) np1,an*a(1,2),an*a(2,2),an*a(3,2)
      write(*,*)np1,an*a(1,3),an*a(2,3),an*a(3,3)
      read(*,*)
C
C   CALCULATE [-(1/2)DEL**2}**3 PSI_0
      DO 13 I=1,NP1
      A(I,4) = -0.5*((I*(I+1)-LRGL*(LRGL+1))*A(I+2,3) - 2.*I*A(I+1,3)/N
     1   + A(I,3)/(N*N))
   13 CONTINUE
C
C   CALCULATE E1.  NOTE THAT A(.,3) = V*Psi_0 AND A(.,1) = Psi_0
      E1 = 0
      DO 14 I=1,NP1
      DO 14 J=1,NP1
      E1 = E1 + A(I,1)*A(J,3)*T(I+J)
   14 CONTINUE
      E1 = -0.5*E1
      E1X = -0.25*(2.D0/N - 1.5D0/N**2)/N**2
      WRITE(*,'('' E1 ='',2D15.7)') E1,E1X
C
C   B(I,1) IS THE COEF. OF R**(I-1) IN PSI(1) FOR V = -P**4/4.
C   APPLY THE RECURSION RELATION, INCLUDING INHOMOGENEOUS TERMS A(.,3) = V*Psi_0
C   CHANGE THIS PART AS APPROPRIATE FOR OTHER PERTURBATION POTENTIALS.
      B(LRGL+1,1) = 0.
      B(LRGL+2,1) = -A(LRGL+1,3)/(0*LRGL+1)
      DO 15 I=LRGL+1,N+2
      B(I+2,1) = -2.*((1-(I+1.D0)/N)*B(I+1,1) + A(I+1,3) +2.*E1*A(I,1))
     1   /((I+1)*(I+2)-LRGL*(LRGL+1))
      write(*,*) i,B(I+2,1)
    1 FORMAT(I4,D15.7)
   15 CONTINUE
C   CHECK THAT THE SERIES TERMINATES
!      IF(DABS(B(N+2,1)/B(N+1,1)).LT.1.D-12) GO TO 20
      WRITE(*,'('' FIRST ORDER P**4 EQUATION NOT SOLVED'',2D12.5)
     1   ') B(N+2,1),B(N+1,1)
      STOP
C   CALCULATE p**4 OVERLAP INTEGRAL.
      OV = 0.
      DO 17 I=1,NP1
      DO 17 J=1,NP1
   17 OV = OV + A(I,1)*B(J,1)*T(I+J+1)
      OV = 0.5*OV
C
C   CALCULATE delta(r) OVERLAP INTEGRAL.
      OVD = 0.
      DO 37 I=1,NP1
      DO 37 J=1,NP2
   37 OVD = OVD + A(I,1)*C(J)*T(I+J) + A(I,1)*CL(J)*TL(I+J)
C
C   ORTHOGONALIZE (NORMALIZE).
      DO 24 J=1,NP1
   24 B(J,1) = B(J,1) - 2.*OV*A(J,1)
C
C   RE-CALCULATE OVERLAP INTEGRAL.
      OV = 0.
      DO 25 I=1,NP1
      DO 25 J=1,NP1
   25 OV = OV + A(I,1)*B(J,1)*T(I+J+1)
      OV = 0.5*OV
C
C   PERFORM OTHER CALCULATIONS WITH THE PERTURBED WAVE FUNCTIONS
C================================================================
      DO 39 I=1,NP2
   39 VB(I,1) = -B(I,1)
C   CALCULATE (p**2/2)*PSI(1).  B(I,2) IS THE COEF. OF R**(I-1) IN
C    (1/R)*(p**2/2)*PSI(1).
C   VB(I,1) IS THE COEF. OF R**(I-1) IN V*PSI(1) = -PSI(1)/R.
      B(1,2) = -0.5*(2.*B(2,1) - 2.*B(1,1)/N)
      VB(1,2) = -0.5*(2.*VB(3,1) - 2.*VB(2,1)/N + VB(1,1)/(N*N))
      DO 22 I=2,NP2
      VB(I,2) = -0.5*(I*(I+1)*VB(I+2,1) - 2.*I*VB(I+1,1)/N
     1   + VB(I,1)/(N*N))
   22 B(I,2) = -0.5*(I*(I+1)*B(I+1,1) - 2.*I*B(I,1)/N + B(I-1,1)/(N*N))
C
C   CALCULATE (p**4/4)*PSI(1).  B(I,3) IS THE COEF. OF R**(I-1) IN
C    (1/R)*(p**4/4)*PSI(1).
      DO 23 I=1,NP2
      B(I,3) = -0.5*(I*(I+1)*B(I+2,2) - 2.*I*B(I+1,2)/N + B(I,2)/(N*N))
   23 CONTINUE
C
C   CALCULATE (-1/4)<PSI(0)|p**4/4|PSI(1)> FOR V = (-1/2)p**4/4.
      WL = 0.
      WR = 0
      WC = 0
      WLL = 0.
      DO 16 I=1,NP1
      DO 16 J=1,NP2
      WR = WR + A(I,1)*B(J,3)*T(I+J)
      WL = WL + A(I,3)*B(J,1)*T(I+J)
      WC = WC + A(I,2)*B(J,2)*T(I+J-1)
      WLL = WLL + A(I,1)*B(J,2)*T(I+J)*E0 - A(I,1)*VB(J,2)*T(I+J)
   16 CONTINUE
      WR = -0.25*WR
      WL = -0.25*WL
      WC = -0.25*WC
      WLL = -0.25*WLL
      WRITE(*,'('' WR,WL,WC,WLL ='',4D15.7)') WR,WL,WC,WLL
C
C   CALCULATE -(1/2)<PSI(0)|p**4/4|PSI(1)> FOR V = delta(r).
      WDL = 0.
      DO 38 I=1,NP2
      DO 38 J=1,NP2
   38 WDL = WDL + A(I,3)*C(J)*T(I+J-1) + A(I,3)*CL(J)*TL(I+J-1)
      WDL = -0.5*WDL
C
      ET1 = 0.
      ET2 = 0.
      DO 18 I=1,NP1
      DO 18 J=1,NP1
      ET1 = ET1 + A(I,4)*A(J,1)*T(I+J)
   18 ET2 = ET2 + A(I,3)*A(J,2)*T(I+J-1)
C      ET = 0.5*(13.*ET1/8. - 5.*ET2/8.)
C      W = WL
C      ET = 0.5*(17.*ET1/8. - 9.*ET2/8.)
C      W = WR
      ET = 0.5*ET1
      W = WC
C
C   CALCULATE (1/R**2)d**3/dr**3 PSI(0)
      DO 50 I=1,N
      A3(I) = I*(I+1)*(I+2)*A(I+3,1) - 3*I*(I+1)*A(I+2,1)/N
     1   + 3*I*A(I+1,1)/N**2 - A(I,1)/N**3
   50 CONTINUE
C
      AA3 = 0.
      DO 51 I=1,N
      DO 51 J=1,N
   51 AA3 = AA3 + A(I,1)*A3(J)*T(I+J-1)
      WRITE(*,*) AA3
C
      ELAM = -0.1250
C      E2(IC) = WLL + ET - OV*E1
C      E2(IC) =  -WR - E1/(2*N**2) + 2.*WLL - OV*E1
      E2(IC) =  ET/2. - E1/(4*N**2) + WLL - OV*E1
     1   + ELAM*(WDL - OVD*E1 - E1D*OV + ELAM*(-4.*(GAM+XLNA)*A(1,1)**2
     2    - E1D*OVD))
C      E2(IC) = W + ET - OV*E1 + ELAM*(WDL - E1*OVD - E1D*OV
C     1   - ELAM*E1D*OVD)
      D(IC,1) = 1.D0/N**3
      D(IC,2) = (2.D0*N**2+1)/(3.D0*N**5)
      D(IC,3) = E0**2
C      D(IC,3) = AA3
      EX2(IC) = -(2.D0/N+6.D0/N**2-12.D0/N**3+5.D0/N**4)/(16*N*N)
      DD(IC) = EX2(IC) - E2(IC)
C      DD(IC) = AA3
      WRITE(*,3) W,OV,OVD,ET,E2(IC),EX2(IC)
    3 FORMAT(' W,OV,OVD,ET =',4D15.7/' E2,EX2 =',2D15.7)
      IF(IC.LE.2) GO TO 100
      CALL SOLVE(D,DD,X,CHK,3)
      WRITE(*,2) X(1),X(2),X(3),CHK
    2 FORMAT(' DELTA FCN. COEFS. =',4D15.7)
      IC = 0
      GO TO 100
      STOP
      END
      SUBROUTINE SOLVE(A,B,X,CHK,N)
C   SOLVES THE EQUATION AX = B BY GAUSS ELIMINATION.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(3,3),B(16),X(16),Y(16),IROW(16)
      DO 16 I=1,N
   16 IROW(I) = I
      DO 10 NR1=2,N
      NR = NR1 - 1
      NRX  = IROW(NR)
      AMAX = DABS(A(NRX,NRX))
      AST = AMAX
      IMAX = NR
      DO 40 I=NR1,N
      IX = IROW(I)
      IF(AMAX.GE.DABS(A(IX,IX))) GO TO 40
      AMAX = DABS(A(IX,IX))
      IMAX = I
   40 CONTINUE
      IF(IMAX.EQ.NR) GO TO 43
      NRX = IROW(IMAX)
      IROW(IMAX) = IROW(NR)
      IROW(NR) = NRX
   41 JMAX = NR
   43 DO 14 J=NR1,N
      JY = IROW(J)
   14 Y(J) = A(NRX,JY)
      DO 10 I=NR1,N
      IX = IROW(I)
      XM = A(IX,NRX)/A(NRX,NRX)
      A(IX,NRX) = XM
      DO 10 J=NR1,N
      JY = IROW(J)
   10 A(IX,JY) = A(IX,JY) - XM*Y(J)
C      WRITE(*,70) (IROW(I),I=1,N)
C      WRITE(6,70) (IROW(I),I=1,N)
   70 FORMAT(' ROW ORDERING'/(10I4))
C
   30 DO 20 NR1=2,N
      NR = NR1 - 1
      NRY = IROW(NR)
      DO 20 I=NR1,N
      IX = IROW(I)
   20 B(IX) = B(IX) - A(IX,NRY)*B(NRY)
C
      NX = IROW(N)
      I1 = IROW(1)
      X(NX) = B(NX)/A(NX,NX)
      CHK = X(NX)*A(I1,NX)
      DO 11 I=2,N
      NR = N - I + 1
      NR1 = NR + 1
      SUM = 0.
      NRX = IROW(NR)
      DO 12 J=NR1,N
      JY = IROW(J)
   12 SUM = SUM + X(JY)*A(NRX,JY)
      X(NRX) = (B(NRX) - SUM)/A(NRX,NRX)
   11 CHK = CHK + X(NRX)*A(I1,NRX)
      IF(B(I1).NE.0.) CHK = CHK/B(I1) - 1.
   13 RETURN
      END
