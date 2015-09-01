C     Compile with:
C     gfortran -ffixed-line-length-0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow bfe.f
      PROGRAM test

          INCLUDE 'bfe.h'

          CHARACTER*10 filename
          INTEGER n,l,m
          CHARACTER*32 coeff_filename
          CHARACTER*32 output_filename

C         Read coeff file name and output filename from command line
          CALL getarg(1, coeff_filename)
          CALL getarg(2, output_filename)

C         Initialize crap
          G = 1.
          one=1.0d0
          two=2.0d0
          pi=4.0d0*ATAN(one)
          twoopi=2.d0/pi
          onesixth=1.d0/6.d0
          tiny=1.d-30
          zero=0.0d0

          OPEN(uincoef,FILE=coeff_filename,STATUS='OLD')

C         First read the time off the first line
          READ(uincoef,*) tt

C         Read in coefficients from file
          DO 130 n=0,nmax
             DO 120 l=0,lmax
                DO 110 m=0,l
                   READ(uincoef,*) sinsum(n,l,m),cossum(n,l,m)
 110            CONTINUE
 120         CONTINUE
 130      CONTINUE

          OPEN(ubodsin,FILE='positions.dat',STATUS='OLD')
          READ(ubodsin,*) nbodies

          DO 10 i=1,nbodies
              READ(ubodsin,*) x(i),y(i),z(i)
 10       CONTINUE
          CLOSE(ubodsin)

          CALL acc_bfe

          OPEN(UNIT=ubodsout,FILE=output_filename,STATUS='NEW')
          DO 30 i=1,nbodies
              WRITE(ubodsout,111) ax(i),ay(i),az(i),pot(i)
30        CONTINUE

111       FORMAT(1x,10(1pe14.6))
          CLOSE(ubodsout)

       STOP
       END

C***********************************************************************
C
C
       SUBROUTINE acc_bfe
C
C
C***********************************************************************
C
C
C     Subroutine to compute accelerations, potential, and density.
C
C
C=======================================================================

       INCLUDE 'bfe.h'

       INTEGER k,l,m,n
       LOGICAL firstc
       REAL*8 anltilde,knl,sinth,sinmphi,cosmphi,phinltil,deltam0,
     &        gammln,arggam,coeflm,factrl,
     &        dblfact,ttemp5,ar,ath,aphi,temp3,temp4,
     &        temp5,temp6,plm,dplm,ultrasp,ultrasp1,ultraspt,clm,
     &        dlm,elm,flm,xi,costh,phi,r,twoalpha,c1,c2,c3,un,unm1,
     &        plm1m,plm2m,cosp,sinp

       DIMENSION ultrasp(0:nmax,0:lmax),
     &           ultraspt(0:nmax,0:lmax),ultrasp1(0:nmax,0:lmax),
     &           anltilde(0:nmax,0:lmax),dblfact(lmax+1),
     &           coeflm(0:lmax,0:lmax),
     &           twoalpha(0:lmax),c1(1:nmax,0:lmax),c2(1:nmax,0:lmax),
     &           c3(1:nmax),cosmphi(0:lmax),sinmphi(0:lmax),
     &           plm(0:lmax,0:lmax),dplm(0:lmax,0:lmax)

       DATA firstc/.TRUE./

       SAVE firstc,dblfact,anltilde,coeflm,twoalpha,c1,c2,c3

C=======================================================================
        IF(firstc) THEN

           firstc=.FALSE.

           dblfact(1) = 1.
           DO 5 l=2,lmax
              dblfact(l) = dblfact(l-1)*(2.*l-1.)
 5         CONTINUE

           DO 20 n=0,nmax
              DO 10 l=0,lmax
                 knl = 0.5*n*(n+4.*l+3.)+(l+1.)*(2.*l+1.)
                 anltilde(n,l) = -2.**(8.*l+6.)*FACTRL(n)*(n+2.*l+1.5)
                 arggam = 2.*l+1.5
                 anltilde(n,l) = anltilde(n,l)*(EXP(GAMMLN(arggam)))**2
                 anltilde(n,l) = anltilde(n,l)/(4.*pi*knl*FACTRL(n+4*l+2))
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

        DO 200 k=1,nbodies

           r=SQRT(x(k)**2+y(k)**2+z(k)**2)
           costh=z(k)/r
           phi=ATAN2(y(k),x(k))
           xi=(r-1.)/(r+1.)

           DO 130 m=0,lmax
              cosmphi(m)=COS(m*phi)
              sinmphi(m)=SIN(m*phi)
 130       CONTINUE

           pot(k)=0.0
           ar=0.0
           ath=0.0
           aphi=0.0

           DO 148 l=0,lmax

              ultrasp(0,l)=1.0
              ultrasp(1,l)=twoalpha(l)*xi
              ultrasp1(0,l)=0.0
              ultrasp1(1,l)=1.0

              un=ultrasp(1,l)
              unm1=1.0

              DO 144 n=1,nmax-1
                 ultrasp(n+1,l)=(c1(n,l)*xi*un-c2(n,l)*unm1)*c3(n)
                 unm1=un
                 un=ultrasp(n+1,l)
                 ultrasp1(n+1,l)=((twoalpha(l)+(n+1)-1.)*unm1-(n+1)*xi*
     &                    ultrasp(n+1,l))/(twoalpha(l)*(1.-xi*xi))
 144          CONTINUE

 148       CONTINUE

           DO 1482 m=0,lmax

              plm(m,m)=1.0
              IF(m.GT.0) plm(m,m)=(-1.)**m*dblfact(m)*SQRT(1.-
     &                            costh*costh)**m
              plm1m=plm(m,m)
              plm2m=0.0

              DO 1481 l=m+1,lmax
                 plm(l,m)=(costh*(2.*l-1.)*plm1m-(l+m-1.)*plm2m)/(l-m)
                 plm2m=plm1m
                 plm1m=plm(l,m)
 1481         CONTINUE

 1482      CONTINUE

           dplm(0,0)=0.0

           DO 1486 l=1,lmax

              DO 1484 m=0,l

                 IF(l.EQ.m) THEN
                    dplm(l,m)=l*costh*plm(l,m)/(costh*costh-1.0)
                 ELSE
                    dplm(l,m)=(l*costh*plm(l,m)-(l+m)*plm(l-1,m))/
     &                        (costh*costh-1.0)
                 ENDIF

 1484         CONTINUE
 1486      CONTINUE

           DO 190 l=lmin,lmax,lskip

              temp3=0.0
              temp4=0.0
              temp5=0.0
              temp6=0.0

              DO 180 m=0,l

                 clm=0.0
                 dlm=0.0
                 elm=0.0
                 flm=0.0

                 DO 150 n=0,nmax
                    clm=clm+ultrasp(n,l)*cossum(n,l,m)
                    dlm=dlm+ultrasp(n,l)*sinsum(n,l,m)
                    elm=elm+ultrasp1(n,l)*cossum(n,l,m)
                    flm=flm+ultrasp1(n,l)*sinsum(n,l,m)
 150             CONTINUE

                 temp3=temp3+plm(l,m)*(clm*cosmphi(m)+dlm*sinmphi(m))
                 temp4=temp4-plm(l,m)*(elm*cosmphi(m)+flm*sinmphi(m))
                 temp5=temp5-dplm(l,m)*(clm*cosmphi(m)+dlm*sinmphi(m))
                 temp6=temp6-m*plm(l,m)*(dlm*cosmphi(m)-clm*sinmphi(m))
 180          CONTINUE

              phinltil=r**l/((1.+r)**(2*l+1))
              pot(k)=pot(k)+temp3*phinltil
C       potext(k)=0.0d0
              ar=ar+phinltil*(-temp3*(l/r-(2.*l+1.)/
     &              (1.+r))+temp4*4.*(2.*l+1.5)/(1.+r)**2)
              ath=ath+temp5*phinltil
              aphi=aphi+temp6*phinltil

 190        CONTINUE

       cosp=COS(phi)
       sinp=SIN(phi)

           sinth=SQRT(1.-costh**2)
           ath= -sinth*ath/r
           aphi=aphi/(r*sinth)
           ax(k)=G*(sinth*cosp*ar+costh*cosp*ath-
     &           sinp*aphi)
           ay(k)=G*(sinth*sinp*ar+costh*sinp*ath+
     &           cosp*aphi)
           az(k)=G*(costh*ar-sinth*ath)
           pot(k)=pot(k)*G

 200    CONTINUE

        RETURN
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
C        SUBROUTINE computecoeff
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
C        INCLUDE 'bfe.h'
C
C        INTEGER k,l,m,n
C        LOGICAL firstc
C        REAL*8 anltilde,knl,sinth,sinmphi,cosmphi,phinltil,deltam0,
C     &         gammln,arggam,coeflm,factrl,
C     &         dblfact,ttemp5,ar,ath,aphi,temp3,temp4,
C     &         temp5,temp6,plm,dplm,ultrasp,ultrasp1,ultraspt,clm,
C     &         dlm,elm,flm,xi,costh,phi,r,twoalpha,c1,c2,c3,un,unm1,
C     &         plm1m,plm2m,cosp,sinp
C
C        DIMENSION ultrasp(0:nmax,0:lmax),
C     &            ultraspt(0:nmax,0:lmax),ultrasp1(0:nmax,0:lmax),
C     &            anltilde(0:nmax,0:lmax),dblfact(lmax+1),
C     &            coeflm(0:lmax,0:lmax),
C     &            twoalpha(0:lmax),c1(1:nmax,0:lmax),c2(1:nmax,0:lmax),
C     &            c3(1:nmax),cosmphi(0:lmax),sinmphi(0:lmax),
C     &            plm(0:lmax,0:lmax),dplm(0:lmax,0:lmax)
C
C        DATA firstc/.TRUE./
C
C        SAVE firstc,dblfact,anltilde,coeflm,
C     &       twoalpha,c1,c2,c3
C
CC=======================================================================
C
C        IF(firstc) THEN
C
C           firstc=.FALSE.
C
C           dblfact(1)=1.
C
C           DO 5 l=2,lmax
C              dblfact(l)=dblfact(l-1)*(2.*l-1.)
C 5         CONTINUE
C
C           DO 20 n=0,nmax
C              DO 10 l=0,lmax
C                 knl=0.5*n*(n+4.*l+3.)+(l+1.)*(2.*l+1.)
C                 anltilde(n,l)=-2.**(8.*l+6.)*FACTRL(n)*(n+2.*l+1.5)
C                 arggam=2.*l+1.5
C                 anltilde(n,l)=anltilde(n,l)*(EXP(GAMMLN(arggam)))**2
C                 anltilde(n,l)=anltilde(n,l)/(4.*pi*knl*FACTRL(n+4*l+2))
C 10           CONTINUE
C 20        CONTINUE
C
C           DO 25 l=0,lmax
C
C              twoalpha(l)=2.0*(2.*l+1.5)
C
C              DO 23 m=0,l
C                 deltam0=2.
C                 IF(m.EQ.0) deltam0=1.
C                 coeflm(l,m)=(2.*l+1.)*deltam0*FACTRL(l-m)/FACTRL(l+m)
C 23           CONTINUE
C 25        CONTINUE
C
C           DO 30 n=1,nmax
C              c3(n)=1.0/(n+1.0)
C
C              DO 27 l=0,lmax
C                 c1(n,l)=2.0*n+twoalpha(l)
C                 c2(n,l)=n-1.0+twoalpha(l)
C 27           CONTINUE
C
C 30        CONTINUE
C
C           lskip=1
C           IF(zeroodd.OR.zeroeven) lskip=2
C
C           lmin=0
C           IF(zeroeven) lmin=1
C
C        ENDIF
C
C        DO 60 l=0,lmax
C           DO 50 m=0,l
C              DO 40 n=0,nmax
C                 sinsum(n,l,m)=0.0
C                 cossum(n,l,m)=0.0
C 40           CONTINUE
C 50        CONTINUE
C 60     CONTINUE
C
C        DO 120 k=1,nbodies
C
C           IF(ibound(k).GT.0)THEN
C              r=SQRT(x(k)**2+y(k)**2+z(k)**2)
C              costh=z(k)/r
C              phi=ATAN2(y(k),x(k))
C              xi=(r-1.)/(r+1.)
C
C              DO 105 m=0,lmax
C                 cosmphi(m)=COS(m*phi)
C                 sinmphi(m)=SIN(m*phi)
C 105          CONTINUE
C
C              DO 113 l=0,lmax
C
C                 ultrasp(0,l)=1.0
C                 ultrasp(1,l)=twoalpha(l)*xi
C
C                 un=ultrasp(1,l)
C                 unm1=1.0
C
C                 DO 111 n=1,nmax-1
C                    ultrasp(n+1,l)=(c1(n,l)*xi*un-c2(n,l)*unm1)*c3(n)
C                    unm1=un
C                    un=ultrasp(n+1,l)
C 111             CONTINUE
C
C                 DO 112 n=0,nmax
C                    ultraspt(n,l)=ultrasp(n,l)*anltilde(n,l)
C 112             CONTINUE
C
C 113          CONTINUE
C
C              DO 1132 m=0,lmax
C
C                 plm(m,m)=1.0
C                 IF(m.GT.0) plm(m,m)=(-1.)**m*dblfact(m)*SQRT(1.-
C     &                costh*costh)**m
C                 plm1m=plm(m,m)
C                 plm2m=0.0
C
C                 DO 1131 l=m+1,lmax
C                    plm(l,m)=(costh*(2.*l-1.)*plm1m-
C     &                        (l+m-1.)*plm2m)/(l-m)
C                    plm2m=plm1m
C                    plm1m=plm(l,m)
C 1131            CONTINUE
C
C 1132         CONTINUE
C
C              DO 118 l=lmin,lmax,lskip
C
C                 temp5=r**l/((1.+r)**(2*l+1))*mass(k)
C
C                 DO 116 m=0,l
C
C                    ttemp5=temp5*plm(l,m)*coeflm(l,m)
C                    temp3=ttemp5*sinmphi(m)
C                    temp4=ttemp5*cosmphi(m)
C
C                    DO 114 n=0,nmax
C                       sinsum(n,l,m)=sinsum(n,l,m)+temp3*ultraspt(n,l)
C                       cossum(n,l,m)=cossum(n,l,m)+temp4*ultraspt(n,l)
C 114                CONTINUE
C
C 116             CONTINUE
C 118          CONTINUE
C           END IF
C 120    CONTINUE
C
