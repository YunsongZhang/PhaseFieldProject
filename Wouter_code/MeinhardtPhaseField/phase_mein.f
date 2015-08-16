      implicit real*8 (a-h,k,o-z)
      
      parameter (nx=100)
      parameter (ny=20)

c u,v,w: used for finite difference
c a,b,c: used for phase field
      real*8 si(0:nx)
      real*8 a(-1:nx+1,-ny-1:ny+1)
      real*8 an(0:nx,-ny:ny)
      real*8 u(-1:nx+1)
      real*8 un(0:nx)
      real*8 v(-1:nx+1)
      real*8 vn(0:nx)
      real*8 w(-1:nx+1)
      real*8 wn(0:nx)
      real*8 b(-1:nx+1,-ny-1:ny+1)
      real*8 bn(0:nx,-ny:ny)
      real*8 c(-1:nx+1,-ny-1:ny+1)
      real*8 cn(0:nx,-ny:ny)
      real*8 yi(-1:nx+1,-ny-1:ny+1),yj(-1:nx+1,-ny-1:ny+1)
      real*8 ypi(-1:nx+1,-ny-1:ny+1),ypj(-1:nx+1,-ny-1:ny+1)
      real*8 y(-1:nx+1,-ny-1:ny+1),yp(-1:nx+1,-ny-1:ny+1)
      character*80 filename

c read in parameters
        filename='input_mein'
        open(unit=77,file=filename,status='unknown')
        read(77,*)sa,sc
        read(77,*)ba,bc
        read(77,*)ra,rb,rc
        read(77,*)dr,dy,rad
        read(77,*)
        read(77,*)a0
        read(77,*)b0
        read(77,*)c0
        read(77,*)
        read(77,*)dx,dt
        read(77,*)diffa,diffb,diffc
        read(77,*)tmax,twrite,tdis
        read(77,*)
        read(77,*)gamma,thresh_yp

      niter=int(tmax/dt)
      iowrite=int(twrite/dt)
      iodis=int(tdis/dt)
      pi=4.d0*datan(1.d0)
      twopi=8.d0*datan(1.d0)

c set threshold for plotting
      thresh=15.

      xl=twopi*rad
      dx=xl/(nx)
      dx2=dx*dx
      dx2=dx*dx
      afac=diffa*dt/(dx2)
      cfac=diffc*dt/(dx2)

       x0=xl/2.d0
            do i=0,nx
              x=i*dx
              si(i)=ra*(1.d0+dy*dcos(twopi*(x-x0)/xl))
              write(1,*)x,si(i),ra,dy,twopi,x0,xl
              u(i)=0.d0
              v(i)=b0
              w(i)=0.d0
            enddo


c phase field 
      do j=-ny-1,ny+1
         do i=-1,nx+1
            if (j.lt.0.) then
              xtmp=-dx*dsqrt(dfloat(j*j))
            else
              xtmp=dx*dsqrt(dfloat(j*j))
            endif
            phitmp=dtanh(gamma*(-xtmp))
            y(i,j)=0.5d0*phitmp +0.5d0
            yp(i,j)=0.25d0*gamma*gamma*(1.d0-phitmp*phitmp)**2
            if (i.eq.0) write(7,*)j,y(i,j),yp(i,j),xtmp,gamma
         enddo
      enddo

      do j=-ny-1,ny
         do i=-1,nx
            yi(i,j)=0.5d0*(y(i+1,j)+y(i,j))
            yj(i,j)=0.5d0*(y(i,j+1)+y(i,j))
            ypi(i,j)=0.5d0*(yp(i+1,j)+yp(i,j))
            ypj(i,j)=0.5d0*(yp(i,j+1)+yp(i,j))
            if (yp(i,j+1).lt.thresh_yp)ypj(i,j)=0.d0
            if (yp(i,j).lt.thresh_yp)ypj(i,j)=0.d0
         enddo
      enddo

c this is for the integral term
      inorm=0
      do j=-ny,ny
         do i=0,nx
            if (yp(i,j).gt.thresh_yp) inorm=inorm+1
         enddo
      enddo

c initial conditions
      do j=-ny-1,ny+1
         do i=0,nx
            a(i,j)=0.d0
            b(i,j)=0.d0
            c(i,j)=0.d0
         enddo
      enddo
            do i=-1,nx+1
            do j=-ny-1,ny+1
               if (yp(i,j).gt.thresh_yp) then
                 b(i,j)=b0
               endif
            enddo
            enddo

c bc
      do j=-ny,ny
         a(-1,j)=a(nx,j)
         b(-1,j)=b(nx,j)
         c(-1,j)=c(nx,j)
         a(nx+1,j)=a(0,j)
         b(nx+1,j)=b(0,j)
         c(nx+1,j)=c(0,j)
      enddo
      do i=0,nx
         a(i,-ny-1)=a(i,ny)
         b(i,-ny-1)=b(i,ny)
         c(i,-ny-1)=c(i,ny)
         a(i,ny+1)=a(i,-ny)
         b(i,ny+1)=b(i,-ny)
         c(i,ny+1)=c(i,-ny)
      enddo
         u(-1)=u(nx)
         v(-1)=v(nx)
         w(-1)=w(nx)
         u(nx+1)=u(0)
         v(nx+1)=v(0)
         w(nx+1)=w(0)

      write(6,*)inorm

c dynamics
      do iter=1,niter
         t=iter*dt

c integral term
         suma=0.
         do j=-ny,ny
            do i=0,nx
               if (yp(i,j).gt.thresh_yp) then
                 suma=suma+a(i,j)
               endif
            enddo
         enddo
         suma=suma/inorm

            sumu=0.d0
            do i=0,nx
              sumu=sumu+u(i)
            enddo
            sumu=sumu/(nx+1)

c evolution of the variables
         do j=-ny,ny
            do i=0,nx
               if (j.eq.0) then
                 xlapu=afac*(u(i+1)-2.d0*u(i)+u(i-1))
                 xlapw=cfac*(w(i+1)-2.d0*w(i)+w(i-1))
                 funcu=si(i)*(u(i)**2/v(i)+ba)/
     1                  ((sc+w(i))*(1.d0+sa*u(i)**2))-ra*u(i)
                 funcv=rb*sumu-rb*v(i)
                 funcw=bc*u(i)-rc*w(i)
                 un(i)=u(i)+ dt*funcu+xlapu
                 vn(i)=v(i)+ dt*funcv
                 wn(i)=w(i)+ dt*funcw+xlapw
               endif
               if (yp(i,j).gt.thresh_yp) then
                 xlapa=afac*(ypi(i,j)*(a(i+1,j)-a(i,j))
     &                 - ypi(i-1,j)*(a(i,j)-a(i-1,j))
     1                 + ypj(i,j)*(a(i,j+1)-a(i,j))
     &                 - ypj(i,j-1)*(a(i,j)-a(i,j-1)))/yp(i,j)
                 xlapc=cfac*(ypi(i,j)*(c(i+1,j)-c(i,j))
     &                 - ypi(i-1,j)*(c(i,j)-c(i-1,j))
     1                 + ypj(i,j)*(c(i,j+1)-c(i,j))
     &                 - ypj(i,j-1)*(c(i,j)-c(i,j-1)))/yp(i,j)
                 
                 funca=si(i)*(a(i,j)**2/b(i,j)+ba)/
     1                  ((sc+c(i,j))*(1.d0+sa*a(i,j)**2))-ra*a(i,j)
                 funcb=rb*suma-rb*b(i,j)
                 funcc=bc*a(i,j)-rc*c(i,j)
                 an(i,j)=a(i,j)+ dt*funca+xlapa
                 bn(i,j)=b(i,j)+ dt*funcb
                 cn(i,j)=c(i,j)+ dt*funcc+xlapc
c            if (i.eq.nx/2.and.j.eq.0)write(11,*)an(i,j),bn(i,j),cn(i,j)
c            if (i.eq.nx/2.and.j.eq.0)write(12,*)funca,funcb,funcc
c    1          ,si(i),ba,sc,sa,ra
c            if (i.eq.nx/2.and.j.eq.0)write(13,*)xlapa,xlapc,suma
               endif
            enddo
         enddo

c update
         do j=-ny,ny
            do i=0,nx
               if (j.eq.0) then
                    u(i)=un(i)
                    v(i)=vn(i)
                    w(i)=wn(i)
               endif
               if (yp(i,j).gt.thresh_yp) then
                    a(i,j)=an(i,j)
                    b(i,j)=bn(i,j)
                    c(i,j)=cn(i,j)
               endif
            enddo
         enddo
c bc
         u(-1)=u(nx)
         v(-1)=v(nx)
         w(-1)=w(nx)
         u(nx+1)=u(0)
         v(nx+1)=v(0)
         w(nx+1)=w(0)
      do j=-ny,ny
         a(-1,j)=a(nx,j)
         b(-1,j)=b(nx,j)
         c(-1,j)=c(nx,j)
         a(nx+1,j)=a(0,j)
         b(nx+1,j)=b(0,j)
         c(nx+1,j)=c(0,j)
      enddo
      do i=0,nx
         a(i,-ny-1)=a(i,ny)
         b(i,-ny-1)=b(i,ny)
         c(i,-ny-1)=c(i,ny)
         a(i,ny+1)=a(i,-ny)
         b(i,ny+1)=b(i,-ny)
         c(i,ny+1)=c(i,-ny)
      enddo

c various plots are produced
         if (mod(iter,iowrite).eq.0) then
            do i=0,nx
               x=i*dx
               write(55,*)x,a(i,0)
            enddo
            do j=-ny,ny
               x=j*dx
               write(56,*)x,a(0,j)
            enddo
            write(55,*)'&'
            write(56,*)'&'

            t=iter*dt
            write(65,*)t,a(0,0)

         endif
         if (mod(iter,iodis).eq.0) then
            write(6,*)iter*dt,a(0,0)
            do jj=0,nx
                if (a(jj,0).gt.thresh)write(44,*)t,jj,a(jj,0)
                if (u(jj).gt.thresh)write(34,*)t,jj,u(jj)
            enddo

         endif

      enddo

            do j=-ny,ny
               x=j*dx
               write(43,*)x,a(0,j)
            enddo

            do i=0,nx
               x=i*dx
               write(30,*)x,u(i)
               write(31,*)x,v(i)
               write(32,*)x,w(i)
               write(40,*)x,a(i,0)
               write(41,*)x,b(i,0)
               write(42,*)x,c(i,0)
               do j=-ny,ny
                  write(48,*)i,j,a(i,j)
               enddo
            enddo
               do j=-ny,ny
                  write(49,*)j,a(0,j)
               enddo

      end
