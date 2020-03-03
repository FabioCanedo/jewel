C     program readfile
C     implicit none
C     write(*,*) 'Hello World!'
C     end program

C++   This is a library built to read external medium profiles ++
C++   for JEWEL                                                ++

C++   This is the main subroutine that performs the reading    ++
C++                                                            ++
      subroutine reader(filename,np,nt,tprofile)
      implicit none
      integer i,np,nt
      integer signal
      character*100 filename
      double precision tprofile(101,101,60)
      
      open(unit=1,file=filename,iostat=signal)
      call readheader
      i=0
      do while(signal.eq.0)
      call readtimestep(signal,i,np,tprofile)
      end do
      end subroutine

C++   This routine reads the header lines and ignores them     ++
C++                                                            ++
      subroutine readheader
      implicit none
      integer i
      character *10 dummyline
      do i=1,8
      read(1,*) dummyline
      end do
      end subroutine

C++   This is the main subroutine that performs the reading    ++
C++   a single timestep                                        ++
      subroutine readtimestep(signal,k,np,tprofile)
      implicit none
      integer i,j,k,np
      integer signal
      double precision tline(np)
      double precision tgrid(np,np)
      character line(20*np)
      double precision tprofile(101,101,60)
      do i=1,np
      read(1,*) tline
      tgrid(i,:)=tline
      enddo
      call inserttimestep(k,np,tgrid,tprofile)
      end subroutine

C++   This is the main subroutine that inserts the timestep    ++
C++   into the tprofile                                        ++
      subroutine inserttimestep(k,np,tgrid,tprofile)
      implicit none
      integer i,j,k,np
      double precision tgrid(np),tprofile(101,101,60)
      double precision x,y,interpol

      do i=1,101
            do j=1,101
                  x=-15.d0+(i-1)*0.3d0
                  y=-15.d0+(j-1)*0.3d0
                  tprofile(i,j,k)=interpol(x,y,tgrid)
            enddo
      enddo
      
      endsubroutine

      double precision function interpol(x,y,np,tgrid)
      implicit none
      integer i,j,ii,jj,iii,jjj,np
      double precision tgrid(np,np),igrid(4,4),xgrid(4,4),ygrid(4,4)
      double precision xmax,xmin,dx
      double precision x,y,xa,ya
      double precision interpol,bicubic
      
      i=1+floor((x-xmin)/dx)
      j=1+floor((y-xmin)/dx)
      xa=mod(x-xmin,dx)
      ya=mod(y-xmin,dx)

      do ii=1,4
            do jj=1,4
                  iii=i+ii-2
                  jjj=j+jj-2
                  if (iii.gt.np.or.iii.lt.1) then
                        igrid(ii,jj)=0.d0
                  else if (jjj.gt.np.or.jjj.lt.1) then
                        igrid(ii,jj)=0.d0
                  else
                        igrid(ii,jj)=tgrid(iii,jjj)
                  end if
                  xgrid(ii,jj)=xmin+(ii-1)*dx
                  ygrid(ii,jj)=xmin+(jj-1)*dx
            enddo
      enddo

      interpol=bicubic(xa,ya,dx,xmin,xmax,igrid,xgrid,ygrid)
      end function
      
      double precision function bicubic(xa,ya,dx,xmin,xmax,igrid,xc,yc)
      implicit none
      integer i,j
      double precision y(2,2),y1(2,2),y2(2,2),y12(2,2)
      double precision igrid(4,4)
      double precision xc(4,4),yc(4,4),ansy
      double precision xa,ya,dx,xmin,xmax
      double precision dertwospline
      
      do i=1,2
            do j=1,2
                  y(i,j)=igrid(i+1,j+1)
                  y1(i,j)=(igrid(i+2,j+1)-igrid(i,j+1))/(xc(i+2,j+1)
     &-xc(i,j+1))
                  y2(i,j)=(igrid(i+1,j+2)-igrid(i+1,j))/(yc(i+1,j+2)
     &-yc(i+1,j))
                  y12(i,j)=(igrid(i+2,j+2)-igrid(i+2,j)-igrid(i,j+2)
     &+igrid(i,j))/(yc(i+1,j+2)-yc(i+1,j))*(xc(i+2,j+1)-xc(i,j+1))
            end do
      end do
      
      bicubic=dertwospline(xc(2:3,2:3),yc(2:3,2:3),y,y1,y2,y12,xa,ya)
      end function

      double precision function dertwospline(x1,x2,y,y1,y2,y12,xa,xb)
      implicit none
      integer i
      double precision x1(2,2)
      double precision x2(2,2)
      double precision y(2,2)
      double precision y1(2,2)
      double precision y2(2,2)
      double precision y12(2,2)
      double precision w(2),w1(2)
      double precision xa,xb
      double precision derspline

      w(1)=derspline(x1(:,1),y(:,1),y1(:,1),xa)
      w(2)=derspline(x1(:,2),y(:,2),y1(:,1),xa)
      w1(1)=derspline(x1(:,1),y2(:,1),y12(:,1),xa)
      w1(2)=derspline(x1(:,2),y2(:,2),y12(:,1),xa)
      
      dertwospline=derspline(x2(1,:),w(:),w1(:),xb)

      end function

      double precision function derspline(x,y,yprime,xval)
      implicit none
      integer i
      double precision x(2),y(2),yprime(2),c(4),yvec(4)
      double precision dx,xval,t
      double precision a(4,4)
      data a/1.d0,0.d0,-3.d0,2.d0,0.d0,0.d0,3.d0,-2.d0,0.d0,1.d0,
     &-2.d0,1.d0,0.d0,0.d0,-1.d0,1.d0/

      dx=x(2)-x(1)
      do i=1,2
            yvec(i)=y(i)
            yvec(i+2)=dx*yprime(i)
      end do
      c=matmul(a,yvec) 

      t=(xval-x(1))/dx

      derspline=0.d0
      do i=1,4
      derspline=derspline+c(i)*t**(i-1)
      end do

      end function
