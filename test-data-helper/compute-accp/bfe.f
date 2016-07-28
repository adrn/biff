C     Compile with:
C     gfortran -ffixed-line-length-0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow bfe.f
      PROGRAM test

          INCLUDE 'bfe.h'

          CHARACTER*10 filename
          INTEGER i,n,l,m
          REAL*8 tmp1,tmp2
          CHARACTER*32 coeff_filename
          CHARACTER*32 output_filename
          CHARACTER*32 pos_filename

C         Read coeff file name and output filename from command line
          CALL getarg(1, coeff_filename)
          CALL getarg(2, pos_filename)
          CALL getarg(3, output_filename)

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

C         First read the number of coeffs off the first line
          READ(uincoef,*) ncoeff

C         Zero out coeff arrays
          DO 867 n=0,nmax
            DO 868 l=0,lmax
              DO 869 m=0,l
                cossum(n,l,m) = zero
                sinsum(n,l,m) = zero
 869      CONTINUE
 868      CONTINUE
 867      CONTINUE

C         Read in coefficients from file
          DO 130 i=1,ncoeff
            READ(uincoef,*) n,l,m,tmp1,tmp2
            cossum(n,l,m) = tmp1
            sinsum(n,l,m) = tmp2
 130      CONTINUE

          OPEN(ubodsin,FILE=pos_filename,STATUS='OLD')
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

C           IF (k.eq.1) THEN
C             WRITE(*,*) ar,ath,aphi
C           ENDIF

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
