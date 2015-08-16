c 1d diffusion, no degradation or reaction
c ic: c=c0 between -dx/2<x<dx/2

      implicit real*8 (a-h,k,o-z)
      
      parameter (nx=200)
      parameter (ny=20)

c u: used for finite difference
c v: used for phase field
      real*8 u(-nx-1:nx+1)
      real*8 un(-nx-1:nx+1)
      real*8 v(-nx-1:nx+1,-ny-1:ny+1)
      real*8 vn(-nx:nx,-ny:ny)
      real*8 yi(-nx-1:nx+1,-ny:ny),yj(-nx-1:nx+1,-ny:ny)
      real*8 ypi(-nx-1:nx+1,-ny:ny),ypj(-nx-1:nx+1,-ny:ny)
      real*8 y(-nx-1:nx+1,-ny:ny),yp(-nx-1:nx+1,-ny:ny)
      character*80 filename

        filename='input_mem'
        open(unit=77,file=filename,status='unknown')
        read(77,*)x0
        read(77,*)
        read(77,*)dx,straal,dt
        read(77,*)diffu
        read(77,*)iodis,iowrite
        read(77,*)tmax
        read(77,*)gamma0,thresh_yp

      niter=int(tmax/dt)
      gamma=gamma0/straal
      pi=4.d0*datan(1.d0)

      dx2=dx*dx
      vfac=diffu*dt/(dx2)
      ufac=diffu*dt/(dx2)

c phase field 
      do j=-ny,ny
         do i=-nx,nx
            xtmp=dx*dsqrt(dfloat(j*j))
            phitmp=dtanh(gamma*(-xtmp))
            y(i,j)=0.5d0*phitmp +0.5d0
            yp(i,j)=0.25d0*gamma*gamma*(1.d0-phitmp*phitmp)**2
            if (i.eq.0) write(7,*)j,y(i,j),yp(i,j),xtmp
         enddo
      enddo

      do j=-ny,ny-1
         do i=-nx,nx-1
            yi(i,j)=0.5d0*(y(i+1,j)+y(i,j))
            yj(i,j)=0.5d0*(y(i,j+1)+y(i,j))
            ypi(i,j)=0.5d0*(yp(i+1,j)+yp(i,j))
            ypj(i,j)=0.5d0*(yp(i,j+1)+yp(i,j))
         enddo
      enddo

      summ=0.d0
      do j=-ny,ny
         do i=-nx,nx
            summ=summ+yp(i,j)
         enddo
      enddo
      summ=summ*dx*dx
      sum=summ/(2.*nx*dx)

      write(6,*)summ,sum

c initial conditions
      do j=-ny,ny
         do i=-nx,nx
            v(i,j)=0.
         enddo
      enddo
            do j=-ny,ny
               if (yp(0,j).gt.thresh_yp) then
                 v(0,j)=x0
               endif
            enddo

      do i=-nx-1,nx+1
                 u(i)=0.
      enddo
                 u(0)=x0/1.0

c dynamics
      do iter=1,niter
         t=iter*dt

c evolution of the interface variables 
         do j=-ny,ny
            do i=-nx,nx
               if (yp(i,j).gt.thresh_yp) then
                 vn(i,j)=v(i,j)
     1           +vfac*(ypi(i,j)*(v(i+1,j)-v(i,j))
     &                 - ypi(i-1,j)*(v(i,j)-v(i-1,j))
     1                 + ypj(i,j)*(v(i,j+1)-v(i,j))
     &                 - ypj(i,j-1)*(v(i,j)-v(i,j-1)))/yp(i,j)
               endif
            enddo
         enddo

         do j=-ny,ny
            do i=-nx,nx
               if (yp(i,j).gt.thresh_yp) then
                    v(i,j)=vn(i,j)
               endif
            enddo
         enddo

c now the fd way
            do i=-nx,nx
              un(i)=u(i)+ufac*(u(i+1)-2.*u(i)+u(i-1))
            enddo
            do i=-nx,nx
              u(i)=un(i)
            enddo


c various plots are produced
         if (mod(iter,iowrite).eq.0) then
            do i=-nx,nx
               x=i*dx
               write(55,*)x,v(i,0)
            enddo
            do j=-ny,ny
               x=j*dx
               write(56,*)x,v(0,j)
            enddo
            write(58,*)0.,u(0)
            write(55,*)'&'
            write(56,*)'&'

            do i=-nx,nx
               write(57,*)i*dx,u(i)
            enddo
            write(57,*)'&'
            write(58,*)'&'
            t=iter*dt
            write(65,*)t,v(0,0)
c this is the analytical solution at x=0
            arg=0.5*dx/sqrt(4.*diffu*t)
            anal=x0*erf(arg)
            write(68,*)t,u(0),anal

         endif
         if (mod(iter,iodis).eq.0) then
            write(6,*)iter*dt,v(0,0)
         endif

      enddo

            do j=-ny,ny
               x=j*dx
               write(43,*)x,v(0,j)
               write(44,*)x,v(5,j)
            enddo

            do i=-nx,nx
               x=i*dx
               write(46,*)x,v(i,0)
               write(47,*)x,u(i)
               do j=-ny,ny
                  write(48,*)i,j,v(i,j)
               enddo
            enddo
               do j=-ny,ny
                  write(49,*)j,v(0,j)
               enddo

            write(9)u,v

      end
