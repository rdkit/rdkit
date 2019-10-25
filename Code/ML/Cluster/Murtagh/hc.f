C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                            C
C  HIERARCHICAL CLUSTERING using (user-specified) criterion. C
C                                                            C
C  Parameters:                                               C
C                                                            C
C  DISS(LEN)         dissimilarities in lower half diagonal  C
C                    storage; LEN = N.N-1/2,                 C
C  IOPT              clustering criterion to be used,        C
C  IA, IB, CRIT      history of agglomerations; dimensions   C
C                    N, first N-1 locations only used,       C
C  MEMBR, NN, DISNN  vectors of length N, used to store      C 
C                    cluster cardinalities, current nearest  C
C                    neighbour, and the dissimilarity assoc. C
C                    with the latter.                        C
C  FLAG              boolean indicator of agglomerable obj./ C
C                    clusters.                               C
C                                                            C
C  F. Murtagh, ESA/ESO/STECF, Garching, February 1986.       C
C                                                            C
C------------------------------------------------------------C
      SUBROUTINE HC(N,LEN,IOPT,IA,IB,CRIT,MEMBR,NN,DISNN,
     X                FLAG,DISS)
      implicit none
      integer n,m,len,iopt
      REAL*8 MEMBR(N),DISS(LEN)
      INTEGER IA(N),IB(N)
      REAL*8 CRIT(N)
      integer nn
      real*8 disnn
      DIMENSION NN(N),DISNN(N)
      LOGICAL FLAG(N)
      REAL*8 INF
      
      integer i,j,ncl,ind,ioffset,jm,im,i2,j2,k
      integer ind1,ind2,ind3,jj
      real*8 dmin,x,xx

      DATA INF/1.E+20/
C
C  Initializations
C
      DO I=1,N
         MEMBR(I)=1.
         FLAG(I)=.TRUE.
      ENDDO
      NCL=N

      if (IOPT.EQ.1) then
         do IND=1,n*(n-1)/2
               DISS(IND)=DISS(IND)/2.
         enddo
      endif
C     (Above is done for the case of the min. var. method
C     where merging criteria are defined in terms of variances
C     rather than distances.)

C
C  Carry out an agglomeration - first create list of NNs
C
      DO I=1,N-1
         DMIN=INF
         DO J=I+1,N
            IND=IOFFSET(N,I,J)
            IF (DISS(IND).GE.DMIN) GOTO 500
               DMIN=DISS(IND)
               JM=J
  500    CONTINUE
         ENDDO
         NN(I)=JM
         DISNN(I)=DMIN
      ENDDO
C
  400 CONTINUE
C     Next, determine least diss. using list of NNs
      DMIN=INF
      DO I=1,N-1
         IF (.NOT.FLAG(I)) GOTO 600
         IF (DISNN(I).GE.DMIN) GOTO 600
            DMIN=DISNN(I)
            IM=I
            JM=NN(I)
  600    CONTINUE
      ENDDO
      NCL=NCL-1
C
C  This allows an agglomeration to be carried out.
C
      I2=MIN0(IM,JM)
      J2=MAX0(IM,JM)
      IA(N-NCL)=I2
      IB(N-NCL)=J2
      CRIT(N-NCL)=DMIN

      ind1 = ioffset(n,i2,j2)
c      write(6,*) "agglom: ",i2,j2,dmin,diss(ind1)
C     
C  Update dissimilarities from new cluster.
C
      FLAG(J2)=.FALSE.
      DMIN=INF
      DO K=1,N
         IF (.NOT.FLAG(K)) GOTO 800
         IF (K.EQ.I2) GOTO 800
         X=MEMBR(I2)+MEMBR(J2)+MEMBR(K)
         IF (I2.LT.K) THEN
                           IND1=IOFFSET(N,I2,K)
                      ELSE
                           IND1=IOFFSET(N,K,I2)
         ENDIF
         IF (J2.LT.K) THEN
                           IND2=IOFFSET(N,J2,K)
                      ELSE
                           IND2=IOFFSET(N,K,J2)
         ENDIF
         IND3=IOFFSET(N,I2,J2)
         XX=DISS(IND3)
C
C  WARD'S MINIMUM VARIANCE METHOD - IOPT=1.
C
         IF (IOPT.EQ.1) THEN
            DISS(IND1)=(MEMBR(I2)+MEMBR(K))*DISS(IND1)+
     X                 (MEMBR(J2)+MEMBR(K))*DISS(IND2)-
     X                 MEMBR(K)*XX
            DISS(IND1)=DISS(IND1)/X
         ENDIF
C
C  SINGLE LINK METHOD - IOPT=2.
C
         IF (IOPT.EQ.2) THEN
            DISS(IND1)=MIN(DISS(IND1),DISS(IND2))
         ENDIF
C
C  COMPLETE LINK METHOD - IOPT=3.
C
         IF (IOPT.EQ.3) THEN
            DISS(IND1)=MAX(DISS(IND1),DISS(IND2))
         ENDIF
C
C  AVERAGE LINK (OR GROUP AVERAGE) METHOD - IOPT=4.
C
         IF (IOPT.EQ.4) THEN
            DISS(IND1)=(MEMBR(I2)*DISS(IND1)+MEMBR(J2)*DISS(IND2))/
     X                 (MEMBR(I2)+MEMBR(J2))
         ENDIF
C
C  MCQUITTY'S METHOD - IOPT=5.
C
         IF (IOPT.EQ.5) THEN
            DISS(IND1)=0.5*DISS(IND1)+0.5*DISS(IND2)
         ENDIF
C
C  MEDIAN (GOWER'S) METHOD - IOPT=6.
C
         IF (IOPT.EQ.6) THEN
            DISS(IND1)=0.5*DISS(IND1)+0.5*DISS(IND2)-0.25*XX
         ENDIF
C
C  CENTROID METHOD - IOPT=7.
C
         IF (IOPT.EQ.7) THEN
            DISS(IND1)=(MEMBR(I2)*DISS(IND1)+MEMBR(J2)*DISS(IND2)-
     X          MEMBR(I2)*MEMBR(J2)*XX/(MEMBR(I2)+MEMBR(J2)))/
     X          (MEMBR(I2)+MEMBR(J2))
            ENDIF
C
         IF (I2.GT.K) GOTO 800
         IF (DISS(IND1).GE.DMIN) GOTO 800
            DMIN=DISS(IND1)
            JJ=K
  800    CONTINUE
      ENDDO
      MEMBR(I2)=MEMBR(I2)+MEMBR(J2)
      DISNN(I2)=DMIN
      NN(I2)=JJ
C
C  Update list of NNs insofar as this is required.
C
      DO I=1,N-1
         IF (.NOT.FLAG(I)) GOTO 900
         IF (NN(I).EQ.I2) GOTO 850
         IF (NN(I).EQ.J2) GOTO 850
         GOTO 900
  850    CONTINUE
C        (Redetermine NN of I:)
         DMIN=INF
         DO J=I+1,N
            IND=IOFFSET(N,I,J)
            IF (.NOT.FLAG(J)) GOTO 870
            IF (I.EQ.J) GOTO 870
            IF (DISS(IND).GE.DMIN) GOTO 870
               DMIN=DISS(IND)
               JJ=J
  870       CONTINUE
         ENDDO
         NN(I)=JJ
         DISNN(I)=DMIN
  900    CONTINUE
      ENDDO
C
C  Repeat previous steps until N-1 agglomerations carried out.
C
      IF (NCL.GT.1) GOTO 400
C
C
      RETURN
      END


C
C
      FUNCTION IOFFSET(N,I,J)
C  Map row I and column J of upper half diagonal symmetric matrix 
C  onto vector.
c      IOFFSET=J+(I-1)*N-(I*(I+1))/2
      if ( j > i) then
         ioffset= (j-1)*(j-2)/2+i
      else
         ioffset= (i-1)*(i-2)/2+j
      endif
      RETURN
      END

