       INTEGER maxnbods,nmax,lmax,lskip,ibound,uincoef,nbodies,ubodsin,
     &         ucoeffout
       PARAMETER(maxnbods=65536,nmax=6,lmin=0,lmax=10,lskip=1)

       REAL*8 x,y,z,vx,vy,vz,mass,pot,
     &        sinsum,cossum,ax,ay,az
       REAL*8 G,one,two,pi,twoopi,onesixth,tiny,zero
       LOGICAL zeroodd,zeroeven

       COMMON /bodscom/ pot(maxnbods),mass(maxnbods),
     &        x(maxnbods),y(maxnbods),z(maxnbods)

       COMMON /coefcom/ sinsum(0:nmax,0:lmax,0:lmax),
     &                  cossum(0:nmax,0:lmax,0:lmax)

       COMMON /parcom/ zeroodd,zeroeven,ibound(maxnbods)
       COMMON /parcomi/ nbodies
       COMMON /parcomr/G,one,two,pi,twoopi,onesixth,tiny,zero

       PARAMETER(ubodsin=12,ucoeffout=13)
