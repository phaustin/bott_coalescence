      parameter (n=400)
      dimension c(n),x(n),z(n),z1(n),z2(n)     
      data z /n*1./
      data z1 /n*1./
      data z2 /n*1./
      open (16,file='plot1d.dat',status='unknown',form='unformatted')
      print *,'enter i0'
      read *,i0
      print *,'c=',360./i0
      do 1000 i=1,n
      c(i)=360./float(i0)
 1000 x(i)=float(i)
      z(372)=5.
      z(373)=5.
      z(374)=5.
      write (16) x,z
      z1(12)=5.
      z1(13)=5.
      z1(14)=5.
      z2(12)=5.
      z2(13)=5.
      z2(14)=5.
      do 1020 i=1,i0
      call adv4p (z1,c)
      call adv4m (z2,c)
 1020 continue
      write (16) z1,z2
      stop
      end

      subroutine adv0 (y,c)
c upstream procedure
      parameter (n=400)
      dimension fm(n),fp(n),y(n),c(n)
      do 1000 i=1,n-1
      fm(i)=-amin1(0.,c(i))*y(i+1)
      fp(i)=amax1(0.,c(i))*y(i)
 1000 continue
      do 1010 i=2,n-1
 1010 y(i)=y(i)-fm(i-1)+fp(i-1)+fm(i)-fp(i)
      return
      end

      subroutine adv2p (y,c)
c area preserving flux form; Bott (1989): Monthly Weather Review.
c second order positive definite version.
c y(i) is transport quantity, input and output.
c boundary conditions are y(1)=const, y(n)=const.
c c(i) is Courant number satisfying the CFL criterion, input.
c fm(i), fp(i) are fluxes for u(i)<0 and u(i)>0, respectively.
c a0, a1, a2, a3, a4 are coefficients of polynomials in gridbox i.
c At i=1 and i=n first order polynomial,
c at 2<=i<=n-1 second order polynomial.
c w(i) are weighting factors.
c the numerical grid is equidistant.
c the procedure is one dimensional.
c for multidimensional applications time splitting has to be used.
c the quantities c(i), fm(i), fp(i)  are given at the right
c boundary of grid cell i.
c Thus, fm(i) is flux from gridbox i+1 into gridbox i for c(i)<0,
c fp(i) is flux from gridbox i into gridbox i+1 for c(i)>0.
      parameter (n=400)
      dimension fm(n),fp(n),w(n),y(n),c(n)
      cr=amax1(0.,c(1))
      fp(1)=amin1(y(1),cr*(y(1)+(1.-cr)*(y(2)-y(1))*0.5))
      w(1)=1.
      do 1000 i=2,n-1
      a0=(26.*y(i)-y(i+1)-y(i-1))/24.
      a1=(y(i+1)-y(i-1))/16.
      a2=(y(i+1)+y(i-1)-2.*y(i))/48.
      cl=-amin1(0.,c(i-1))
      x1=1.-2.*cl
      x2=x1*x1
      fm(i-1)=amax1(0.,a0*cl-a1*(1.-x2)+a2*(1.-x1*x2))
      cr=amax1(0.,c(i))
      x1=1.-2.*cr
      x2=x1*x1
      fp(i)=amax1(0.,a0*cr+a1*(1.-x2)+a2*(1.-x1*x2))
      w(i)=y(i)/amax1(fm(i-1)+fp(i)+1.e-15,a0+2.*a2)
 1000 continue
      cl=-amin1(0.,c(n-1))
      fm(n-1)=amin1(y(n),cl*(y(n)-(1.-cl)*(y(n)-y(n-1))*0.5))
      w(n)=1.
      do 1010 i=2,n-1
 1010 y(i)=y(i)-(fm(i-1)+fp(i))*w(i)+fm(i)*w(i+1)+fp(i-1)*w(i-1)
      return
      end

      subroutine adv4p (y,c)
c area preserving flux form; Bott (1989): Monthly Weather Review.
c fourth order positive definite version.
c y(i) is transport quantity, input and output.
c boundary conditions are y(1)=const, y(n)=const.
c c(i) is Courant number satisfying the CFL criterion, input.
c fm(i), fp(i) are fluxes for u(i)<0 and u(i)>0, respectively.
c a0, a1, a2, a3, a4 are coefficients of polynomials in gridbox i.
c At i=1 and i=n first order polynomial,
c at i=2 and i=n-1 second order polynomial,
c at 3<=i<=n-2 fourth order polynomial.
c w(i) are weighting factors.
c the numerical grid is equidistant.
c the procedure is one dimensional.
c for multidimensional applications time splitting has to be used.
c the quantities c(i), fm(i), fp(i)  are given at the right
c boundary of grid cell i.
c Thus, fm(i) is flux from gridbox i+1 into gridbox i for c(i)<0,
c fp(i) is flux from gridbox i into gridbox i+1 for c(i)>0.
      parameter (n=400)
      dimension y(n),c(n),fm(n),fp(n),w(n)
      cr=amax1(0.,c(1))
      fp(1)=amin1(y(1),cr*(y(1)+(1.-cr)*(y(2)-y(1))*0.5))
      w(1)=1.
      a0=(26.*y(2)-y(3)-y(1))/24.
      a1=(y(3)-y(1))/16.
      a2=(y(3)+y(1)-2.*y(2))/48.
      cl=-amin1(0.,c(1))
      x1=1.-2.*cl
      x2=x1*x1
      fm(1)=amax1(0.,a0*cl-a1*(1.-x2)+a2*(1.-x1*x2))
      cr=amax1(0.,c(2))
      x1=1.-2.*cr
      x2=x1*x1
      fp(2)=amax1(0.,a0*cr+a1*(1.-x2)+a2*(1.-x1*x2))
      w(2)=y(2)/amax1(fm(1)+fp(2)+1.e-15,a0+2.*a2)
      do 1000 i=3,n-2
      a0=(9.*(y(i+2)+y(i-2))-116.*(y(i+1)+y(i-1))
     &   +2134.*y(i))/1920.
      a1=(-5.*(y(i+2)-y(i-2))+34.*(y(i+1)-y(i-1)))/384.
      a2=(-y(i+2)+12.*(y(i+1)+y(i-1))-22.*y(i)-y(i-2))/384.
      a3=(y(i+2)-2.*(y(i+1)-y(i-1))-y(i-2))/768.
      a4=(y(i+2)-4.*(y(i+1)+y(i-1))+6.*y(i)+y(i-2))/3840.
      cl=-amin1(0.,c(i-1))
      x1=1.-2.*cl
      x2=x1*x1
      x3=x1*x2
      fm(i-1)=amax1(0.,a0*cl-a1*(1.-x2)+a2*(1.-x3)
     &        -a3*(1.-x1*x3)+a4*(1.-x2*x3))
      cr=amax1(0.,c(i))
      x1=1.-2.*cr
      x2=x1*x1
      x3=x1*x2
      fp(i)=amax1(0.,a0*cr+a1*(1.-x2)+a2*(1.-x3)
     &      +a3*(1.-x1*x3)+a4*(1.-x2*x3))
      w(i)=y(i)/amax1(fm(i-1)+fp(i)+1.e-15,y(i))
 1000 continue
      a0=(26.*y(n-1)-y(n)-y(n-2))/24.
      a1=(y(n)-y(n-2))/16.
      a2=(y(n)+y(n-2)-2.*y(n-1))/48.
      cl=-amin1(0.,c(n-2))
      x1=1.-2.*cl
      x2=x1*x1
      fm(n-2)=amax1(0.,a0*cl-a1*(1.-x2)+a2*(1.-x1*x2))
      cr=amax1(0.,c(n-1))
      x1=1.-2.*cr
      x2=x1*x1
      fp(n-1)=amax1(0.,a0*cr+a1*(1.-x2)+a2*(1.-x1*x2))
      w(n-1)=y(n-1)/amax1(fm(n-2)+fp(n-1)+1.e-15,a0+2.*a2)
      cl=-amin1(0.,c(n-1))
      fm(n-1)=amin1(y(n),cl*(y(n)-(1.-cl)*(y(n)-y(n-1))*0.5))
      w(n)=1.
      do 1010 i=2,n-1
 1010 y(i)=y(i)-(fm(i-1)+fp(i))*w(i)+fm(i)*w(i+1)+fp(i-1)*w(i-1)
      return
      end

      subroutine adv4m (y,c)
c area preserving flux form; Bott (1989): Monthly Weather Review.
c fourth order monotone version.
c y(i) is transport quantity, input and output.
c boundary conditions are y(1)=const, y(n)=const.
c c(i) is Courant number satisfying the CFL criterion, input.
c fm(i), fp(i) are fluxes for u(i)<0 and u(i)>0, respectively.
c a0, a1, a2, a3, a4 are coefficients of polynomials in gridbox i.
c At i=1 and i=n first order polynomial,
c at i=2 and i=n-1 second order polynomial,
c at 3<=i<=n-2 fourth order polynomial.
c w(i) are weighting factors.
c the numerical grid is equidistant.
c the procedure is one dimensional.
c for multidimensional applications time splitting has to be used.
c the quantities c(i), fm(i), fp(i)  are given at the right
c boundary of grid cell i.
c Thus, fm(i) is flux from gridbox i+1 into gridbox i for c(i)<0,
c fp(i) is flux from gridbox i into gridbox i+1 for c(i)>0.
      parameter (n=400)
      dimension a0(n),a1(n),a2(n),a3(n),a4(n),
     &          y(n),c(n),fm(n),fp(n),w(n)
      a0(2)=(26.*y(2)-y(3)-y(1))/24.
      a1(2)=(y(3)-y(1))/16.
      a2(2)=(y(3)+y(1)-2.*y(2))/48.
      a3(2)=0.
      a4(2)=0.
      do 1000 i=3,n-2
      a0(i)=(9.*(y(i+2)+y(i-2))-116.*(y(i+1)+y(i-1))
     &      +2134.*y(i))/1920.
      a1(i)=(-5.*(y(i+2)-y(i-2))+34.*(y(i+1)-y(i-1)))/384.
      a2(i)=(-y(i+2)+12.*(y(i+1)+y(i-1))-22.*y(i)-y(i-2))/384.
      a3(i)=(y(i+2)-2.*(y(i+1)-y(i-1))-y(i-2))/768.
 1000 a4(i)=(y(i+2)-4.*(y(i+1)+y(i-1))+6.*y(i)+y(i-2))/3840.
      a0(n-1)=(26.*y(n-1)-y(n)-y(n-2))/24.
      a1(n-1)=(y(n)-y(n-2))/16.
      a2(n-1)=(y(n)+y(n-2)-2.*y(n-1))/48.
      a3(n-1)=0.
      a4(n-1)=0.
      w(1)=1.
      w(n)=1.
      fm(n)=0.
      cl=-amin1(0.,c(n-1))
      fm(n-1)=amin1(y(n),cl*(y(n)-(1.-cl)*(y(n)-y(n-1))*0.5))
      clm=cl
      do 1020 i=n-1,2,-1
      cl=clm
      clm=-amin1(0.,c(i-1))
      x1=1.-2.*cl
      x2=x1*x1
      x3=x1*x2
      ymin=amin1(y(i),y(i+1))
      ymax=amax1(y(i),y(i+1))
      fmim=amax1(0.,a0(i)*cl-a1(i)*(1.-x2)+a2(i)*(1.-x3)
     &     -a3(i)*(1.-x1*x3)+a4(i)*(1.-x2*x3))
      fmim=amin1(fmim,y(i)-ymin+fm(i))
      fmim=amax1(fmim,y(i)-ymax+fm(i))
      fm(i-1)=amax1(0.,fmim-(cl-clm)*y(i))
 1020 continue
      cr=amax1(0.,c(1))
      fp(1)=amin1(y(1),cr*(y(1)+(1.-cr)*(y(2)-y(1))*0.5))
      crp=cr
      do 1030 i=2,n-1
      cr=crp
      crp=amax1(0.,c(i))
      x1=1.-2.*cr
      x2=x1*x1
      x3=x1*x2
      ymin=amin1(y(i-1),y(i))
      ymax=amax1(y(i-1),y(i))
      fpi=amax1(0.,a0(i)*cr+a1(i)*(1.-x2)+a2(i)*(1.-x3)
     &    +a3(i)*(1.-x1*x3)+a4(i)*(1.-x2*x3))
      fpi=amin1(fpi,y(i)-ymin+fp(i-1))
      fpi=amax1(fpi,y(i)-ymax+fp(i-1))
      fp(i)=amax1(0.,fpi-(cr-crp)*y(i))
      w(i)=y(i)/amax1(fm(i-1)+fp(i)+1.e-15,y(i))
 1030 continue
      do 1040 i=2,n-1
 1040 y(i)=y(i)-(fm(i-1)+fp(i))*w(i)+fm(i)*w(i+1)+fp(i-1)*w(i-1)
      return
      end
