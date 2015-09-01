       INTEGER nbods,nmax,lmax,lskip,ibound,uincoef,nbodies,ubodsin,ubodsout
C       PARAMETER(nbods=1024,nmax=6,lmin=0,lmax=4,lskip=1)
       PARAMETER(nbods=65536,nmax=3,lmin=0,lmax=6,lskip=1)

       REAL*8 x,y,z,vx,vy,vz,mass,pot,
     &        sinsum,cossum,ax,ay,az
       REAL*8 G,one,two,pi,twoopi,onesixth,tiny,zero
       LOGICAL zeroodd,zeroeven

       COMMON /bodscom/ pot(nbods),mass(nbods),
     &        x(nbods),y(nbods),z(nbods),
     &        ax(nbods),ay(nbods),az(nbods)

       COMMON /coefcom/ sinsum(0:nmax,0:lmax,0:lmax),
     &                  cossum(0:nmax,0:lmax,0:lmax)

       COMMON /parcom/ zeroodd,zeroeven,ibound(nbods)
       COMMON /parcomi/ nbodies
       COMMON /parcomr/G,one,two,pi,twoopi,onesixth,tiny,zero

       PARAMETER(ubodsin=12,ubodsout=13)
