!===============================================================================
!     Formation and growth of ice fog particles at Fairbanks, Alaska
!     P.J.Huffman, T, Ohtake
!
!     Comments: H.Tang
!     Date: Apr 27 2020
!     Email: htang2@alaska.edu
!===============================================================================
!===============================================================================
!     * This program did not use IMPLICIT NONE, so many variables are not
!     declared explicitly. I will explain thoses variables when they are used in
!     the program. The unit system used is the C.G.S system (centimeter, gram,
!     second).
!===============================================================================
!     R,    Radius of droplets (Equation 17);
!     II,   (I, Fortran is a case-insensitive programming language, the author
!           used double letter to represent the capital version of it.)
!           The homogeneous rate of embryo formation (Equation 18);
!     ES,   Saturation vapor pressure;
!     L,    Heat of vaporization and sublimation;
!     SIGMA,Surface tension of water;
!     U,    u+v, because these two quantities always appear together in the
!           differential equation, the author used one letter to represent them;
!     SUM,  A variable used to store the integral term in Equation 17 (It is
!           NOT a good practice to use DOUBLE PRECISION here as other single
!           precision variables will cause small differences that might cause
!           numerical errors;
!     A,    A=R/M in Equation 8;
!     B,    B=M/R in Equation 8;
!     C,    C=1/R^2 in Equation 8;
!     D,    D=16*pi*eta*M^2/(3*R^3);
!     N,    Lines of data to be read.
!===============================================================================
      DIMENSION R(500),II(500),ES(500),L(500),SIGMA(500),U(500)
      REAL L, L1, L2, II
      DOUBLE PRECISION SUM 
      A=5.8E7
      B=2.16E-7
      C=2.28E20
      D=5.67E3
      N=121
!===============================================================================
!     1101: This is a loop that read all external data. Apparently "K" is a
!     constant that amplifies original data. And apparently in the "era of
!     punchcard" people had to simplify everything as much as they can so I
!     don't blame them.
!===============================================================================
      DO 1101 I=1,N
      READ(5,2)ES(I),K,L(I),SIGMA(I),U(I),KK
2     FORMAT(F4.2,I2,2X,F4.2,5X,F4.1,2X,F3.1,I2)
      ES(I)=ES(I)*10.**K
      L(I)=L(I)*10.**10
      U(I)=U(I)*10.**KK
1101  CONTINUE
!===============================================================================
!     999:  Initialization, set up all the parameters required
!     SOLD,   "old" S or S at step i, saturation ratio starts from 1;
!     II,     mentioned before;
!     DELT,   "delta t", change in time;
!     KNT,    "count", loop counter;
!===============================================================================
999   CONTINUE
      SOLD=1.
      II(1)=0.
      DELT=1.0E-6
      KNT=0
      READ(5,1)AB,TI,TF
1     FORMAT(3F10.0)
      IF(AB.EQ.O.)STOP
      WRITE(6,1001)AB,TI,TF
1001  FORMAT(10X,F8.3,2F8.0)
!===============================================================================
!     1100: Initialize the array R with zeros. It is a bad practice of FORTRAN77
!           It is better to write
!                     DO 1100 N=1,500
!                     R(N)=0.
!               1100  CONTINUE
!===============================================================================
      DO 1100 N=1,500
1100  R(N)=0.
!===============================================================================
!     100 & 101:
!     The two main loops of the program.
!     DELT,   Timestep, T(i+1) - T(i)
!     T1,     Time at step i
!     T2,     Time at step i+1
!     DELTA,  Timedifference between T2 and T1
!     KOUNT,  Loop counter for loop 100
!     KNT,    Loop counter for loop 101
!     KNT1,   Only appeared four times in this program and three of them are
!             used for output purposes. This counter is created so in loop 103
!             N >= 1. We don't need it...
!     DTTDT,  dT/dt, derived by taking derivative of Equation 14
!     TT1,    Temperature at step i
!     TT2,    Temperature at step i+1
!     A1,     Temperature difference between T1 and 229K (-44 Celsius), since
!             229K is set to be the ambient temperature of the simulation.
!     J1,     Integer part of A1
!     DEL,    Difference between the A1 and INT(A1), therefore a small value.
!===============================================================================
100   DELT=10.*DELT
      T1=10.*DELT
      IF(T1.GE.9.9)GO TO 999
      KOUNT=0
!===============================================================================
!     The first block:
!
!     TT1=(TI-TF)/((TI-TF)*AB*T1+1.)+TF
!     A1=TT1-229.
!     J1=A1
!     DEL=A1-FLOAT(J1)
!     ES1=ES(J1)+DEL*(ES(J1+1)-ES(J1))
!     SIGMA1=SIGMA(J1)+DEL*(SIGMA(J1+1)-SIGMA(J1))
!     L1=L(J1)+DEL*(L(J1+1)-L(J)
!     U1=U(J1)+DEL*(U(J1+1)-U(J1))
!
!     This block is used as the "base" of step i, if KOUNT != 0 that means
!     the "base" step is already calculated and therefore we move on and
!     compute all the "2" variables.
!
!     The reason that these lines:
!
!     DEL=A1-FLOAT(J1)
!     ES1=ES(J1)+DEL*(ES(J1+1)-ES(J1))
!     SIGMA1=SIGMA(J1)+DEL*(SIGMA(J1+1)-SIGMA(J1))
!     L1=L(J1)+DEL*(L(J1+1)-L(J)
!     U1=U(J1)+DEL*(U(J1+1)-U(J1))
!
!     exist is because the dimension of the input data is 121 and es, sigma, L
!     and u+v all depends on temperature (the index), but the dimension of the
!     outupt data is 500, so the author had to interpolate the 121 elements
!     array on to a 500 elements array.
!===============================================================================
101   T2=T1+DELT
      KOUNT=KOUNT+1
      IF(KOUNT.EQ.90)T2=100.*DELT
      DELTA=T2-T1
      KNT=KNT+1
      KNT1=KNT+1
      DTTDT=-AB*(TI-TF)**2/(TI-TF)*AB*T1+1.)**2
      TT2=(TI-TF)/((TI-TF)*AB*T2+1.)+TF
      IF(KOUNT.NE.1)GO TO 104
      TT1=(TI-TF)/((TI-TF)*AB*T1+1.)+TF
      A1=TT1-229.
      J1=A1
      DEL=A1-FLOAT(J1)
      ES1=ES(J1)+DEL*(ES(J1+1)-ES(J1))
      SIGMA1=SIGMA(J1)+DEL*(SIGMA(J1+1)-SIGMA(J1))
      L1=L(J1)+DEL*(L(J1+1)-L(J1))
      U1=U(J1)+DEL*(U(J1+1)-U(J1))
      DEBUG T1,TT1,J1,ES1,SIGMA1,L1,U1
104   CONTINUE 
      A2=TT2-229.
      J2=A2
      DEL=A2-FLOAT(J2)
      ES2=ES(J2)+DEL*(ES(J2+1)-ES(J2))
      SIGMA2=SIGMA(J2)+DEL*(SIGMA(J2+1)-SIGMA(J2))
      L2=L(J2)+DEL*(L(J2+1)-L(J2))
      U2=U(J2)+DEL*(U(J2+1)-U(J2))
102   SUM=SUM+R(N)*II(N)
      DSDT=A*TT1/ES1*(1.-SOLD)/U1*SUM-B*L1*SOLD*DTTDT/(TT1**2)
      SNEW=SOLD+DELTA*DSDT
!===============================================================================
!     Calculation of radius:
!
!     ROOT=R(N)**2+2.*(SNEW-1.)*DELTA/U2
!     IF(ROOT.LE.0.)GO TO 1103
!     R(N)=SQRT(ROOT)
!
!     These equations are derived from Equation 11 of [Rooth, 1957], the
!     original equation is:
!
!     dr/dt = K/(rho * r) * (e_a - e_0)
!
!     Notice that the Equation 10 in [Rooth, 1957] and the Equation 16 in
!     [Huffman, 1971] are comparable, they are
!
!     dm/dt = 4 * pi * K * r * (e_a - e_0) * r  /(s + r) [Rooth, 1957]
!     dm/dt = 4 * pi / (u+v) * (S   -   1) * r^2/(f + r) [Huffman, 1971]
!
!     Since S is defined as S = e_a/e_0, we can derive that the "K" in
!     [Rooth, 1957] is
!
!     K = 1 / ((u + v) * e_0)
!
!     and plug this into the ODE of r, we have
!
!     dr/dt = 1/r * 1/(rho * (u + v)) * (S - 1)
!
!     solve this ODE for r, we have
!
!     r^2 / 2 = (S - 1)t / (rho * (u + v))
!
!     where rho is the density of water, in C.G.S system rho = 1, thus
!
!     r^2 = 2 * (S - 1)t / (u + v)
!     => r2^2 - r1^2 = 2 * (S - 1) * Delta_t / (u + v)
!     => r2^2 = r1^2 + 2 * (S - 1) * Delta_t / (u + v)
!
!     This is where "ROOT=R(N)**2+2.*(SNEW-1.)*DELTA/U2" comes from.
!===============================================================================
      DO 103 N=1,KNT1
      IF(SNEW.LT.1.)SNEW=1.
      ROOT=R(N)**2+2.*(SNEW-1.)*DELTA/U2
      IF(ROOT.LE.0.)GO TO 1103
      R(N)=SQRT(ROOT)
      GO TO 103
1103  R(N)=0.
103   CONTINUE 
      II(KNT1)=C*SQRT(SIGMA2)*ES2**2*SNEW**2/TT2*EXP(-D*SIGMA2**3
      /(ALOG(SNEW)**2*TT2**3))
      IF(SNEW.LE.1.)II(KNT1)=0.
      WRITE(6,50)T1,T2,TT1,TT2,SOLD,SNEW,II(KNT),II(KNT1)
50    FORMAT(*)
49    FORMAT(*)
      T1=T2
      TT1=TT2
      J1=J2
      ES1=ES2
      L1=L2
      U1=U2
      SIGMA1=SIGMA2
      SOLD=SNEW
      IF(KOUNT.EQ.90)GO TO 100
      GO TO 101
      END
!===============================================================================
!     I will make a table here to explain how this program works:
!
!     T1        T2          DELT        KNOUT       DELTA        101    100
!     ---------------------------------------------------------------------
!     1.E-5     2.E-5       1.E-5       1           1.E-5       GOTO
!     2.E-5     3.E-5       1.E-5       2           1.E-5       GOTO
!     ...       ...         ...         ...         ...         ...
!     89.E-5    90.E-5      1.E-5       89          1.E-5       GOTO
!     90.E-5    1.E-3       1.E-5       90          1.E-5              GOTO
!
!     1.E-3     1.1E-3      1.E-4       1           1.E-4       GOTO
!     1.1E-3    1.2E-3      1.E-4       2           1.E-4       GOTO
!     ...       ...         ...         ...         ...         ...
!     1.89E-3   1.90E-3     1.E-4       89          1.E-4       GOTO
!     1.90E-3   1.E-2       1.E-4       90          1.E-4              GOTO
!
!     1.E-2     1.1E-2      1.E-3       1           1.E-3       GOTO
!     1.1E-1    1.2E-2      1.E-3       2           1.E-3       GOTO
!     ...       ...         ...         ...         ...         ...
!     1.89E-1   1.90E-3     1.E-3       89          1.E-3       GOTO
!     1.90E-1   1.E-1       1.E-3       90          1.E-3              GOTO
!
!     1.E-1     1.1E-1      1.E-2       1           1.E-2       GOTO
!     1.1E-1    1.2E-1      1.E-2       2           1.E-2       GOTO
!     ...       ...         ...         ...         ...         ...
!     1.89E-1   1.90E-1     1.E-2       89          1.E-2       GOTO
!     1.90E-1   1.          1.E-2       90          1.E-2              GOTO
!
!     1.        1.1         1.E-1       1           1.E-1       GOTO
!     1.1       1.2         1.E-1       2           1.E-1       GOTO
!     ...       ...         ...         ...         ...         ...
!     9.9       10.0        1.E-1       89          1.E-1       STOP   STOP
!===============================================================================
