c 2d phase field model : 
c state variables evolve only at the interface
c cGMP diffuses only inside the cell, is producted at the interface
c cAMP does NOT diffuse inside the cell

      implicit real*8 (a-h,k,o-z)
      
      parameter (n=200)
      parameter (nphi=564)

      real*8 u(0:nphi+1)
      real*8 un(0:nphi+1)
      real*8 v(-n-1:n+1,-n-1:n+1)
      real*8 vn(-n:n,-n:n)
      real*8 yi(-n-1:n+1,-n:n),yj(-n-1:n+1,-n:n)
      real*8 ypi(-n-1:n+1,-n:n),ypj(-n-1:n+1,-n:n)
      real*8 y(-n-1:n+1,-n:n),yp(-n-1:n+1,-n:n)
      real*8 yc(-n-1:n+1,-n:n)
      character*80 filename
      real*8 ix(1000),iy(1000),ang(1000)

        filename='input_mem'
        open(unit=77,file=filename,status='unknown')
        read(77,*)x0
        read(77,*)
        read(77,*)dx,straal,dt
        read(77,*)diffv
        read(77,*)iodis,iowrite
        read(77,*)tmax
        read(77,*)gamma0,thresh_yp

c read in file to parameterize circle in rectangular grid
        filename='perimeter_200'
        open(unit=78,file=filename,status='unknown')
          read(78,*)i_perim
        do i=1,i_perim
          read(78,*)ix(i),iy(i),ang(i)
        enddo

      ibound=n/2
      dx=straal/ibound
      r0=dfloat(ibound)*dx
      niter=int(tmax/dt)
      gamma=gamma0/straal
      pi=4.d0*datan(1.d0)

      dx2=dx*dx
      vfac=diffv*dt/(dx2)
      dphi=2.*pi/nphi
      dphi2=dphi*dphi*straal*straal
      diff_renor=1./(straal*straal)

c phase field 
      do j=-n,n
         do i=-n,n
            xtmp=dx*dsqrt(dfloat(i*i+j*j))
            phitmp=dtanh(gamma*(r0-xtmp))
            y(i,j)=0.5d0*phitmp +0.5d0
            yp(i,j)=0.25d0*gamma*gamma*(1.d0-phitmp*phitmp)**2
            yc(i,j)=-0.5d0*phitmp +0.5d0
            if (i.eq.0) write(7,*)j*dx,y(i,j),yp(i,j),xtmp
         enddo
      enddo

      do j=-n,n-1
         do i=-n,n-1
            yi(i,j)=0.5d0*(y(i+1,j)+y(i,j))
            yj(i,j)=0.5d0*(y(i,j+1)+y(i,j))
            ypi(i,j)=0.5d0*(yp(i+1,j)+yp(i,j))
            ypj(i,j)=0.5d0*(yp(i,j+1)+yp(i,j))
         enddo
      enddo

      summ=0.d0
      do j=-n,n
         do i=-n,n
            summ=summ+yp(i,j)
         enddo
      enddo
      summ=summ*dx*dx
      sum=summ/(2.d0*pi*straal)

      nx=2*n+1

c initial conditions for phase field
      do j=-n,n
         do i=-n,n
          if (yp(i,j).gt.thresh_yp) then
            v(i,j)=0.
          else
            v(i,j)=0.
          endif
         enddo
      enddo

            do i=0,n
               if (yp(i,0).gt.thresh_yp) then
                 v(i,0)=x0
               endif
            enddo

c initial conditions for finite difference
      do i=0,nphi+1
                 u(i)=0.
      enddo
                 x0_fd=x0*nphi/(n*pi)
                 u(nphi/2)=x0_fd

c check that total amount is conserved
         call amount(n,y,yp,v,totv)
         write(6,*)t,totv*dx2/sum
         write(61,*)t,totv*dx2/sum

c dynamics
      do iter=1,niter
         t=iter*dt

c evolution of the interface variables 
         do j=-n,n
            do i=-n,n
               if (yp(i,j).gt.thresh_yp) then
                 vn(i,j)=v(i,j)
     1           +vfac*(ypi(i,j)*(v(i+1,j)-v(i,j))
     &                 - ypi(i-1,j)*(v(i,j)-v(i-1,j))
     1                 + ypj(i,j)*(v(i,j+1)-v(i,j))
     &                 - ypj(i,j-1)*(v(i,j)-v(i,j-1)))/yp(i,j)
               endif
            enddo
         enddo

         do j=-n,n
            do i=-n,n
               if (yp(i,j).gt.thresh_yp) then
                    v(i,j)=vn(i,j)
               endif
            enddo
         enddo

c now the finite difference way
            do i=1,nphi
              un(i)=u(i)+dt*diffv*(u(i+1)-2.*u(i)+u(i-1))/dphi2
            enddo
            do i=1,nphi
              u(i)=un(i)
            enddo
            u(0)=u(nphi)
            u(nphi+1)=u(1)

         if (mod(iter,iowrite).eq.0) then
            do i=-n,n
               x=i*dx
               write(55,*)x,v(i,0)
            enddo
            write(55,*)'&'

            do i=1,nphi
               write(57,*)i*dphi-pi,u(i)
            enddo
            do i=1,i_perim
               write(58,*)ang(i),v(ix(i),iy(i))
            enddo
            write(57,*)'&'
            write(58,*)'&'
            write(67,*)iter*dt,v(ibound,0)
            t=iter*dt

c analytical solution at phi=nphi/2
            arg=0.5*dphi/sqrt(4.*diff_renor*t)
            anal=x0_fd*erf(arg)
            write(68,*)t,u(nphi/2),anal

c check that total amount is conserved
         call amount(n,y,yp,v,totv)
         write(6,*)iter,totv*dx2/sum
         write(61,*)iter,totv*dx2/sum
         endif
         if (mod(iter,iodis).eq.0) then
            write(6,*)iter*dt,v(ibound,0)
         endif

      enddo

c cut across the domain
            do i=-n,n
               x=i*dx
               write(43,*)x,v(i,0),v(i,10)
            enddo

c profile along circle for phase field
            do i=1,i_perim
               write(46,*)ang(i),v(ix(i),iy(i))
            enddo
c profile along finite difference circle
            do i=1,nphi
               write(47,*)i*dphi-pi,u(i)
               ampl=100.
               psi=i*dphi-pi
               write(49,*)ampl*cos(psi),ampl*sin(psi),u(i)
            enddo
            do i=-n,n
            do j=-n,n
              if (v(i,j).gt.1.d-4)  write(48,*)i,j,v(i,j)
            enddo
            enddo

      end

      subroutine amount(n,y,yp,v,totv)
      implicit real*8 (a-h,k,o-z)
      real*8 v(-n-1:n+1,-n-1:n+1)
      real*8 y(-n-1:n+1,-n:n)
      real*8 yp(-n-1:n+1,-n:n)

      totv=0.
      do j=-n,n
         do i=-n,n
           totv=totv+yp(i,j)*v(i,j)
         enddo
      enddo

      return
      end
