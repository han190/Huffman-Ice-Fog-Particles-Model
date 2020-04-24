      DIMENSION R(500),II(500),ES(500),L(500),SIGMA(500),U(500)
      REAL L, L1, L2, II
      DOUBLE PRECISION SUM 
      A=5.8E7
      B=2.16E-7
      C=2.28E20
      D=5.67E3
      N=121
c
c     loop 1101: This is a loop that read all external data.
c     Apparently "K" is a constant that amplifies original data. And apparently
c     in the "era of punchcard" people had to simplify everything as much as 
c     they can so I don't blame them.
c
      DO 1101 I=1,N
      READ(5,2)ES(I),K,L(I),SIGMA(I),U(I),KK
2     FORMAT(F4.2,I2,2X,F4.2,5X,F4.1,2X,F3.1,I2)
      ES(I)=ES(I)*10.**K
      L(I)=L(I)*10.**10
      U(I)=U(I)*10.**KK
1101  CONTINUE
c
c     end of loop 1101
c

c
c     Initialization, set up all the parameters required
c     SOLD is the orginal S, or S at step i, correspondingly we have SNEW later
c     meaninig the new S, or S at step i+1. II means captial I and is 
c     initialized to zero. DELT "delta_t" is the initial time step. KNT "count"
c     is an integer that represents step. AB is the product of the coefficient
c     a and b mentioned in the paper.
c
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
      DO 1100 N=1,500
1100  R(N)=0.
c
c     End of initialization
c

c
c     Beginning of the main loop
c     The amazing number 229. when computing the variable A1 and A2 is set
c     because the ambient temperature of the program is -44 degree Celsius.
c     (It took me quite a while to figure out... when you write a code, please
c     write a documentation, or at least comment on something others obviously
c     don't understand.)
c     
100   DELT=10.*DELT
      T1=10.*DELT
      IF(T1.GE.9.9)GO TO 999
      KOUNT=0
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
c     DEBUG T1,TT1,J1,ES1,SIGMA1,L1,U1
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

