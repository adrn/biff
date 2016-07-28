C     Compile with:
C     gfortran -ffixed-line-length-0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow computecoeff.f
      PROGRAM test

          INCLUDE 'computecoeff.h'

          CHARACTER*10 filename
          INTEGER i,n,l,m
          REAL*8 tmp1,tmp2
          CHARACTER*32 pos_filename
          CHARACTER*32 out_coeff_filename

C         Read coeff file name and output filename from command line
          CALL getarg(1, pos_filename)
          CALL getarg(2, out_coeff_filename)

C         Initialize crap
          G = 1.
          one=1.0d0
          two=2.0d0
          pi=4.0d0*ATAN(one)
          twoopi=2.d0/pi
          onesixth=1.d0/6.d0
          tiny=1.d-30
          zero=0.0d0

          OPEN(ubodsin,FILE=pos_filename,STATUS='OLD')
          READ(ubodsin,*) nbodies

          DO 10 i=1,nbodies
              READ(ubodsin,*) x(i),y(i),z(i)
              mass(i) = 1.0d0 / nbodies
 10       CONTINUE
          CLOSE(ubodsin)

C         Compute coefficients from positions

          CALL computecoeff

          OPEN(UNIT=ucoeffout,FILE=out_coeff_filename,STATUS='NEW')
          DO 977 n=0,nmax
            DO 978 l=0,lmax
              DO 979 m=0,l
                WRITE(ucoeffout,*) n, l, m, cossum(n,l,m), sinsum(n,l,m)
 979      CONTINUE
 978      CONTINUE
 977      CONTINUE

112       FORMAT(1x,(I2,I2,I2,1pe14.6,1pe14.6))

          CLOSE(ucoeffout)


       STOP
       END

C***********************************************************************
C
C
        FUNCTION FACTRL(N)
C
C
C***********************************************************************
C
C
C     A function to compute factorials.  (From numerical recipes.)
C
C
C=======================================================================

        INTEGER n,ntop,j
        REAL*8 factrl,a,gammln,arggam

        DIMENSION A(33)

        DATA NTOP,A(1)/0,1./

        IF (N.LT.0) THEN
          PAUSE 'negative factorial'
        ELSE IF (N.LE.NTOP) THEN
          FACTRL=A(N+1)
        ELSE IF (N.LE.32) THEN
          DO 11 J=NTOP+1,N
            A(J+1)=J*A(J)
11        CONTINUE
          NTOP=N
          FACTRL=A(N+1)
        ELSE
          arggam=n+1.
          FACTRL=EXP(GAMMLN(arggam))
        ENDIF

        RETURN
        END

C***********************************************************************
C
C
        FUNCTION GAMMLN(XX)
C
C
C***********************************************************************
C
C
C     A routine to compute the natural logarithm of the gamma
C     function.  (Taken from numerical recipes.)
C
C
C=======================================================================

        INTEGER j

        REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER,gammln,xx

        DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     &      -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
        DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/

        X=XX-ONE
        TMP=X+FPF
        TMP=(X+HALF)*LOG(TMP)-TMP
        SER=ONE

        DO 11 J=1,6
          X=X+ONE
          SER=SER+COF(J)/X
11      CONTINUE

        GAMMLN=TMP+LOG(STP*SER)

        RETURN
        END

C***********************************************************************
C
C
        SUBROUTINE computecoeff
C
C
C***********************************************************************
C
C
C     Subroutine to compute expansion coefficients
C
C
C=======================================================================
C

        INCLUDE 'computecoeff.h'

        INTEGER k,l,m,n
        LOGICAL firstc
        REAL*8 anltilde,knl,sinth,sinmphi,cosmphi,phinltil,deltam0,
     &         gammln,arggam,coeflm,factrl,
     &         dblfact,ttemp5,ar,ath,aphi,temp3,temp4,
     &         temp5,temp6,plm,dplm,ultrasp,ultrasp1,ultraspt,clm,
     &         dlm,elm,flm,xi,costh,phi,r,twoalpha,c1,c2,c3,un,unm1,
     &         plm1m,plm2m,cosp,sinp

        DIMENSION ultrasp(0:nmax,0:lmax),
     &            ultraspt(0:nmax,0:lmax),ultrasp1(0:nmax,0:lmax),
     &            anltilde(0:nmax,0:lmax),dblfact(lmax+1),
     &            coeflm(0:lmax,0:lmax),
     &            twoalpha(0:lmax),c1(1:nmax,0:lmax),c2(1:nmax,0:lmax),
     &            c3(1:nmax),cosmphi(0:lmax),sinmphi(0:lmax),
     &            plm(0:lmax,0:lmax),dplm(0:lmax,0:lmax)

        DATA firstc/.TRUE./

        SAVE firstc,dblfact,anltilde,coeflm,
     &       twoalpha,c1,c2,c3

C=======================================================================

        IF(firstc) THEN

           firstc=.FALSE.

           dblfact(1)=1.

           DO 5 l=2,lmax
              dblfact(l)=dblfact(l-1)*(2.*l-1.)
 5         CONTINUE

           DO 20 n=0,nmax
              DO 10 l=0,lmax
                 knl=0.5*n*(n+4.*l+3.)+(l+1.)*(2.*l+1.)
                 anltilde(n,l)=-2.**(8.*l+6.)*FACTRL(n)*(n+2.*l+1.5)
                 arggam=2.*l+1.5
                 anltilde(n,l)=anltilde(n,l)*(EXP(GAMMLN(arggam)))**2
                 anltilde(n,l)=anltilde(n,l)/(4.*pi*knl*FACTRL(n+4*l+2))
 10           CONTINUE
 20        CONTINUE

           DO 25 l=0,lmax

              twoalpha(l)=2.0*(2.*l+1.5)

              DO 23 m=0,l
                 deltam0=2.
                 IF(m.EQ.0) deltam0=1.
                 coeflm(l,m)=(2.*l+1.)*deltam0*FACTRL(l-m)/FACTRL(l+m)
 23           CONTINUE
 25        CONTINUE

           DO 30 n=1,nmax
              c3(n)=1.0/(n+1.0)

              DO 27 l=0,lmax
                 c1(n,l)=2.0*n+twoalpha(l)
                 c2(n,l)=n-1.0+twoalpha(l)
 27           CONTINUE

 30        CONTINUE

        ENDIF

        DO 60 l=0,lmax
           DO 50 m=0,l
              DO 40 n=0,nmax
                 sinsum(n,l,m)=0.0
                 cossum(n,l,m)=0.0
 40           CONTINUE
 50        CONTINUE
 60     CONTINUE

        DO 120 k=1,nbodies

           IF(1.GT.0)THEN
              r=SQRT(x(k)**2+y(k)**2+z(k)**2)
              costh=z(k)/r
              phi=ATAN2(y(k),x(k))
              xi=(r-1.)/(r+1.)

              DO 105 m=0,lmax
                 cosmphi(m)=COS(m*phi)
                 sinmphi(m)=SIN(m*phi)
 105          CONTINUE

              DO 113 l=0,lmax

                 ultrasp(0,l)=1.0
                 ultrasp(1,l)=twoalpha(l)*xi

                 un=ultrasp(1,l)
                 unm1=1.0

                 DO 111 n=1,nmax-1
                    ultrasp(n+1,l)=(c1(n,l)*xi*un-c2(n,l)*unm1)*c3(n)
                    unm1=un
                    un=ultrasp(n+1,l)
 111             CONTINUE

                 DO 112 n=0,nmax
                    ultraspt(n,l)=ultrasp(n,l)*anltilde(n,l)
 112             CONTINUE

 113          CONTINUE

              DO 1132 m=0,lmax

                 plm(m,m)=1.0
                 IF(m.GT.0) plm(m,m)=(-1.)**m*dblfact(m)*SQRT(1.-
     &                costh*costh)**m
                 plm1m=plm(m,m)
                 plm2m=0.0

                 DO 1131 l=m+1,lmax
                    plm(l,m)=(costh*(2.*l-1.)*plm1m-
     &                        (l+m-1.)*plm2m)/(l-m)
                    plm2m=plm1m
                    plm1m=plm(l,m)
 1131            CONTINUE

 1132         CONTINUE

              DO 118 l=lmin,lmax,lskip

                 temp5=r**l/((1.+r)**(2*l+1))*mass(k)

                 DO 116 m=0,l

                    ttemp5=temp5*plm(l,m)*coeflm(l,m)
                    temp3=ttemp5*sinmphi(m)
                    temp4=ttemp5*cosmphi(m)

                    DO 114 n=0,nmax
                       sinsum(n,l,m)=sinsum(n,l,m)+temp3*ultraspt(n,l)
                       cossum(n,l,m)=cossum(n,l,m)+temp4*ultraspt(n,l)
 114                CONTINUE

 116             CONTINUE
 118          CONTINUE
           END IF
 120    CONTINUE


         RETURN
         END
