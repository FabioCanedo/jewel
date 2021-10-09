C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++ Copyright (C) 2017 Korinna C. Zapp [Korinna.Zapp@cern.ch]       ++
C++                                                                 ++
C++ This file is part of JEWEL 2.2.0                                ++
C++                                                                 ++
C++ The JEWEL homepage is jewel.hepforge.org                        ++
C++                                                                 ++
C++ The medium model was partly implemented by Jochen Klein.        ++
C++ Raghav Kunnawalkam Elayavalli helped with the implementation    ++
C++ of the V+jet processes.                                         ++
C++                                                                 ++
C++ Please follow the MCnet GUIDELINES and cite Eur.Phys.J. C74     ++
C++ (2014) no.2, 2762 [arXiv:1311.0048] for the code and            ++
C++ JHEP 1303 (2013) 080 [arXiv:1212.1599] and                      ++
C++ optionally EPJC 60 (2009) 617 [arXiv:0804.3568] for the         ++
C++ physics. The reference for V+jet processes is EPJC 76 (2016)    ++
C++ no.12 695 [arXiv:1608.03099] and for recoil effects it is       ++
C++ arXiv:1707.01539.
C++                                                                 ++
C++ JEWEL relies heavily on PYTHIA 6 for the event generation. The  ++
C++ modified version of PYTHIA 6.4.25 that is distributed with      ++
C++ JEWEL is, however, not an official PYTHIA release and must not  ++
C++ be used for anything else. Please refer to results as           ++
C++ "JEWEL+PYTHIA".                                                 ++
C++                                                                 ++
C++ JEWEL also uses code provided by S. Zhang and J. M. Jing        ++
C++ (Computation of Special Functions, John Wiley & Sons, New York, ++
C++ 1996 and http://jin.ece.illinois.edu) for computing the         ++
C++ exponential integral Ei(x).                                     ++
C++                                                                 ++
C++                                                                 ++
C++ JEWEL  is free software; you can redistribute it and/or         ++
C++ modify it under the terms of the GNU General Public License     ++
C++ as published by the Free Software Foundation; either version 2  ++
C++ of the License, or (at your option) any later version.          ++
C++                                                                 ++
C++ JEWEL is distributed in the hope that it will be useful,        ++
C++ but WITHOUT ANY WARRANTY; without even the implied warranty of  ++
C++ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the    ++
C++ GNU General Public License for more details.                    ++
C++                                                                 ++
C++ You should have received a copy of the GNU General Public       ++  
C++ License along with this program; if not, write to the Free      ++
C++ Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, ++
C++ MA 02110-1301 USA                                               ++
C++                                                                 ++
C++ Linking JEWEL statically or dynamically with other modules is   ++
C++ making a combined work based on JEWEL. Thus, the terms and      ++
C++ conditions of the GNU General Public License cover the whole    ++
C++ combination.                                                    ++
C++                                                                 ++
C++ In addition, as a special exception, I give you permission to   ++
C++ combine JEWEL with the code for the computation of special      ++
C++ functions provided by S. Zhang and J. M. Jing. You may copy and ++
C++ distribute such a system following the terms of the GNU GPL for ++
C++ JEWEL and the licenses of the other code concerned, provided    ++
C++ that you include the source code of that other code when and as ++
C++ the GNU GPL requires distribution of source code.               ++
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE MEDINIT(FILE,id,etam,mass)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDFILEC/MEDFILE,NLIST,endoff
      CHARACTER*200 MEDFILE
      INTEGER NLIST
      logical endoff
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--longitudinal boost of momentum distribution
	common/boostmed/boost
	logical boost
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      common/temperature/tempfac
      double precision tempfac
C--nuclear thickness function
      COMMON /THICKFNC/ RMAX,TA(100,2)
      DOUBLE PRECISION RMAX,TA
C--geometrical cross section
      COMMON /CROSSSEC/ IMPMAX,CROSS(200,3)
      DOUBLE PRECISION IMPMAX,CROSS
C--identifier of log file
	common/logfile/logfid
	integer logfid
C--grid parameters
      common/gridpar/ gdt,gdx,gxmax,gxmin,gnx,gny,gnt
      double precision gdt,gdx,gxmax,gxmin
      integer gnx,gny,gnt

      DATA RAU/10./
      DATA D3/0.9d0/
      DATA ZETA3/1.2d0/
C--local variables
      INTEGER I,LUN,POS,IOS,id,mass
	double precision etam
      CHARACTER*100 BUFFER,LABEL,tempbuf
	CHARACTER*100 FILE
	character firstchar
	logical fileexist

	etamax2 = etam
	logfid = id

      IOS=0
      LUN=77

      NLIST=1
      endoff=.false.

C--default settings
      TAUI=0.6d0
      TI=0.36d0
      TC=0.17d0
      WOODSSAXON=.TRUE.
      MODMED=.TRUE.
      MEDFILELIST=.FALSE.
      CENTRMIN=0.d0
      CENTRMAX=10.d0
      NF=3
      A=mass
      N0=0.17d0
      D=0.54d0
      SIGMANN=6.2
	MDFACTOR=0.45d0
	MDSCALEFAC=0.9d0
      tempfac=1.0d0
	boost = .true.
      breal=0.d0

C--read settings from file
	write(logfid,*)
	inquire(file=FILE,exist=fileexist)
	if(fileexist)then
        write(logfid,*)'Reading medium parameters from ',FILE
        OPEN(unit=LUN,file=FILE,status='old',err=10)
	  do 20 i=1,1000
          READ(LUN, '(A)', iostat=ios) BUFFER
	     if (ios.ne.0) goto 30
	     firstchar = buffer(1:1)
	     if (firstchar.eq.'#') goto 20
          POS=SCAN(BUFFER,' ')
          LABEL=BUFFER(1:POS)
          BUFFER=BUFFER(POS+1:)
          IF (LABEL=="TAUI")THEN
            READ(BUFFER,*,IOSTAT=IOS) TAUI
          ELSE IF (LABEL=="TI") THEN
            READ(BUFFER,*,IOSTAT=IOS) TI
          ELSE IF (LABEL=="TC") THEN
            READ(BUFFER,*,IOSTAT=IOS) TC
          ELSE IF (LABEL=="WOODSSAXON") THEN
            READ(BUFFER,*,IOSTAT=IOS) WOODSSAXON
          ELSE IF (LABEL=="MODMED") THEN
            READ(BUFFER,*,IOSTAT=IOS) MODMED
          ELSE IF (LABEL=="MEDFILE") THEN
            READ(BUFFER,'(50A)',IOSTAT=IOS) MEDFILE
          ELSE IF (LABEL=="CENTRMIN") THEN
            READ(BUFFER,*,IOSTAT=IOS) CENTRMIN
          ELSE IF (LABEL=="CENTRMAX") THEN
            READ(BUFFER,*,IOSTAT=IOS) CENTRMAX
          ELSE IF (LABEL=="NF") THEN
            READ(BUFFER,*,IOSTAT=IOS) NF
          ELSE IF (LABEL=="N0") THEN
            READ(BUFFER,*,IOSTAT=IOS) N0
          ELSE IF (LABEL=="D") THEN
            READ(BUFFER,*,IOSTAT=IOS) D
          ELSE IF (LABEL=="SIGMANN") THEN
            READ(BUFFER,*,IOSTAT=IOS) SIGMANN
          ELSE IF (LABEL=="MDFACTOR") THEN
            READ(BUFFER,*,IOSTAT=IOS) MDFACTOR
          ELSE IF (LABEL=="MDSCALEFAC") THEN
            READ(BUFFER,*,IOSTAT=IOS) MDSCALEFAC
          ELSE IF (LABEL=="BREAL") THEN
            READ(BUFFER,*,IOSTAT=IOS) breal
	     else
       write(logfid,*)'unknown label ',label
	     endif
 20	  continue

 30	  close(LUN,status='keep')
	  write(logfid,*)'...done'
	  goto 40

 10     write(logfid,*)'Could not open medium parameter file, '//
     &	'will run with default settings.'

	else
	  write(logfid,*)'No medium parameter file found, '//
     &	'will run with default settings.'
	endif

 40   write(logfid,*)'using parameters:'
      write(logfid,*)'TAUI        =',TAUI
      write(logfid,*)'TI          =',TI
      write(logfid,*)'TC          =',TC
      write(logfid,*)'WOODSSAXON  =',WOODSSAXON
      write(logfid,*)'MODMED      =',MODMED
      write(logfid,*)'MEDFILELIST =',MEDFILELIST
      write(logfid,*)'CENTRMIN    =',CENTRMIN
      write(logfid,*)'CENTRMAX    =',CENTRMAX
      write(logfid,*)'NF          =',NF
      write(logfid,*)'A           =',A
      write(logfid,*)'N0          =',N0
      write(logfid,*)'D           =',D
      write(logfid,*)'SIGMANN     =',SIGMANN
      write(logfid,*)'MDFACTOR    =',MDFACTOR
      write(logfid,*)'MDSCALEFAC  =',MDSCALEFAC
      write(logfid,*)'BREAL       =',breal
      write(logfid,*)
      write(logfid,*)
      write(logfid,*)

C--calculate T_A(x,y)
      CALL CALCTA
C--calculate geometrical cross section
      CALL CALCXSECTION

      IF(MODMED) THEN
      CALL MYMED()
      END IF
     
      END

      SUBROUTINE MEDNEXTEVT
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED
C--geometrical cross section
      COMMON /CROSSSEC/ IMPMAX,CROSS(200,3)
      DOUBLE PRECISION IMPMAX,CROSS
C--local variables
      integer i,j
      DOUBLE PRECISION PYR,R,b1,b2,gettemp

C--pick an impact parameter
      if(.not. modmed) then
      r=(pyr(0)*(centrmax-centrmin)+centrmin)/100.
      i=0
      do 130 j=1,200
       if ((r-cross(j,3)/cross(200,3)).ge.0.) then
        i=i+1
       else 
        goto 132
       endif
 130  continue
 132  continue
      b1 = (i-1)*0.1d0
      b2 = i*0.1d0
      breal = (b2*(cross(i,3)/cross(200,3)-r)
     &      +b1*(r-cross(i+1,3)/cross(200,3)))/
     &	(cross(i,3)/cross(200,3)-cross(i+1,3)/cross(200,3))
      !write(*,*) "Parâmetro de impacto = ",breal
      centr = r;
      
      end if

      !IF(MODMED) THEN
      !write(*,*) "Reading medium"
      !CALL MYMED()
      !write(*,*) "Finished reading medium"
      !ENDIF

      END

      double precision function getcentrality()
      implicit none
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      getcentrality=centr
      end



      SUBROUTINE PICKVTX(X,Y)
      IMPLICIT NONE
      DOUBLE PRECISION X,Y
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
C--internal medium parametersa
      COMMON/TEMPMAX/TEMPMAXIMUM
      DOUBLE PRECISION TEMPMAXIMUM
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED
C--local variables
      DOUBLE PRECISION X1,X2,Y1,Y2,Z,XVAL,YVAL,ZVAL,NTHICK,PYR,MEDPART  &
     & ,rval
      integer i,j,k
      double precision s,gettemp
      logical myprob

      myprob=.true.
      !tempmaximum=0.d0
        
      IF(myprob) THEN
100   zval=pyr(0)

      !write(*,*) "zval = ",zval
      !write(*,*) "Pick vtx called"
      !write(*,*) prob(834**2)


      k=0
      do k=1,834**2
            if(prob(k).gt.zval) then
              i=k/834
              j=mod(k,834)
              xval=-25.d0+(i-1)*(50.d0/833.d0)
              yval=-25.d0+(j-1)*(50.d0/833.d0)
              rval=(xval**2+yval**2)**(0.5d0)
              if(rval.gt.1.12d0*(208.0d0**(0.333))) then
                    write(*,*) 'Too far from the center'
              go to 100
              elseif(gettemp(xval,yval,0.0d0,taui).lt.1.0d0*tc) then
                    write(*,*) 'Too cold'
              go to 100
              endif
              !write(*,*) "xvtx=",xval," yvtx=",yval
              x=xval
              y=yval
              return
            endif
      enddo



      ELSE
      X1=BREAL/2.d0-RAU
      X2=RAU-BREAL/2.d0
      Y1=-SQRT(4*RAU**2-BREAL**2)/2.d0
      Y2=SQRT(4*RAU**2-BREAL**2)/2.d0
 131  XVAL=PYR(0)*(X2-X1)+X1
      YVAL=PYR(0)*(Y2-Y1)+Y1
      IF((NTHICK(XVAL-BREAL/2.,YVAL).EQ.0.d0).OR.
     &     NTHICK(XVAL+BREAL/2.,YVAL).EQ.0.d0) GOTO 131
      ZVAL=PYR(0)*NTHICK(-BREAL/2.d0,0d0)*NTHICK(BREAL/2.d0,0d0)
      Z=NTHICK(XVAL-BREAL/2.d0,YVAL)*NTHICK(XVAL+BREAL/2.d0,YVAL)
      IF(ZVAL.GT.Z) GOTO 131
      END IF
      X=XVAL
      Y=YVAL
      END

	SUBROUTINE SETB(BVAL)
	IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
	DOUBLE PRECISION BVAL
	BREAL=BVAL
	END



      SUBROUTINE GETSCATTERER(X,Y,Z,T,TYPE,PX,PY,PZ,E,MS)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
C--internal medium parameters
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED
C--longitudinal boost of momentum distribution
	common/boostmed/boost
	logical boost
C--function calls
      DOUBLE PRECISION GETTEMP,GETMD,GETMOM,GETMS
C--identifier of log file
	common/logfile/logfid
	integer logfid
C--local variables
      DOUBLE PRECISION X,Y,Z,T,MS,PX,PY,PZ,E,MD,TEMP
      double precision u,ux,uy,px2,py2,px3,py3,e3
      double precision getux,getuy,eta
      INTEGER TYPE
      DOUBLE PRECISION R,PYR,pmax,wt,tau,theta,phi,pi,p,ys,pz2,e2
      double precision e4
      DATA PI/3.141592653589793d0/

      R=PYR(0)
      IF(R.LT.(2.*12.*NF*D3/3.)/(2.*12.*NF*D3/3.+3.*16.*ZETA3/2.))THEN
         TYPE=2
      ELSE
         TYPE=21
      ENDIF
      MS=GETMS(X,Y,Z,T)
      MD=GETMD(X,Y,Z,T)
      TEMP=GETTEMP(X,Y,Z,T)
	tau=sqrt(t**2-z**2)
	if (boost) then
  	  ys = 0.5d0*log((t+z)/(t-z))
	else
	  ys = 0.d0
	endif
	pmax = 10.*temp

      IF(TEMP.LT.1.D-2)THEN
       write(logfid,*)'asking for a scattering centre without medium:'
       write(logfid,*)'at (x,y,z,t)=',X,Y,Z,T
       write(logfid,*)'making one up to continue but '//
     &	'something is wrong!'
       TYPE=21
       PX=0.d0
       PY=0.d0
       PZ=0.d0
       MS=GETMS(0.d0,0.d0,0.d0,0.d0)
       MD=GETMD(0.d0,0.d0,0.d0,0.d0)
       E=SQRT(PX**2+PY**2+PZ**2+MS**2)
       RETURN
      ENDIF

 10	p = pyr(0)**0.3333333*pmax
	E2 = sqrt(p**2+ms**2)
	if (type.eq.2) then
	  wt = (exp(ms/temp)-1.)/(exp(E2/temp)-1.)
	else
	  wt = (exp(ms/temp)+1.)/(exp(E2/temp)+1.)
	endif
	if (wt.gt.1.) write(logfid,*)'Error in getscatterer: weight = ',wt
	if (wt.lt.0.) write(logfid,*)'Error in getscatterer: weight = ',wt
	if (pyr(0).gt.wt) goto 10
	phi = pyr(0)*2.*pi
	theta = -acos(2.*pyr(0)-1.)+pi
	px  = p*sin(theta)*cos(phi)
	py  = p*sin(theta)*sin(phi)
	pz2 = p*cos(theta)

      !Getting local velocity to perform boost
      ux=getux(x,y,z,t)
      uy=getuy(x,y,z,t)
      u=sqrt(ux**2+uy**2)
      if(ux.ne.0.d0) then
            theta=atan(uy/ux)
      else
            if(uy.gt.0.d0) then
                  theta=pi/2.d0
            else
                  theta=-pi/2.d0
            end if
      end if
      eta=0.5d0*log((1.d0+u)/(1.d0-u))

      px2=px*cos(-theta)-py*sin(-theta)
      py2=px*sin(-theta)+py*cos(-theta)
      E2=E

      px3=px2*cosh(eta)-E2*sinh(eta)
      py3=py2
      E3=E2*cosh(eta)-px2*sinh(eta)
      
      px=px3*cos(theta)-py3*sin(theta)
      py=px3*sin(theta)+py3*cos(theta)
      E4=E3

	E   = cosh(ys)*E4 + sinh(ys)*pz2
	pz  = sinh(ys)*E4 + cosh(ys)*pz2

      END




      SUBROUTINE AVSCATCEN(X,Y,Z,T,PX,PY,PZ,E,m)
      IMPLICIT NONE
C--longitudinal boost of momentum distribution
	common/boostmed/boost
	logical boost
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--local variables
	double precision x,y,z,t,px,py,pz,e,getms,m,ys
	if (boost) then
  	  ys = 0.5*log((t+z)/(t-z))
	  if ((z.eq.0.d0).and.(t.eq.0.d0)) ys =0.d0
	  if (ys.gt.etamax2) ys=etamax2
	  if (ys.lt.-etamax2) ys=-etamax2
	else
	  ys = 0.d0
	endif
	m  = getms(x,y,z,t)
	e  = m*cosh(ys)
	px = 0.d0
	py = 0.d0
	pz = m*sinh(ys)
	end


      SUBROUTINE maxscatcen(PX,PY,PZ,E,m)
      IMPLICIT NONE
C--longitudinal boost of momentum distribution
	common/boostmed/boost
	logical boost
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--local variables
	double precision px,py,pz,e,getmsmax,m,ys
	if (boost) then
  	  ys = etamax2
	else
	  ys = 0.d0
	endif
	m  = getmsmax()
	e  = m*cosh(ys)
	px = 0.d0
	py = 0.d0
	pz = m*sinh(ys)
	end
	


      DOUBLE PRECISION FUNCTION GETMD(X1,Y1,Z1,T1)
      IMPLICIT NONE
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION X1,Y1,Z1,T1,GETTEMP
      double precision getmdmin
      GETMD=MDSCALEFAC*3.*GETTEMP(X1,Y1,Z1,T1)
      !GETMD=MAX(GETMD,MDFACTOR)
      GETMD=MAX(GETMD,GETMDMIN())
      END



      DOUBLE PRECISION FUNCTION GETMS(X2,Y2,Z2,T2)
      IMPLICIT NONE
      DOUBLE PRECISION X2,Y2,Z2,T2,GETMD
      GETMS=GETMD(X2,Y2,Z2,T2)/SQRT(2.)
      END



      DOUBLE PRECISION FUNCTION GETNEFF(X3,Y3,Z3,T3)
      IMPLICIT NONE
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED
      common/temperature/tempfac
      double precision tempfac
C--   local variables
      DOUBLE PRECISION X3,Y3,Z3,T3,PI,GETTEMP,tau,cosheta
      double precision getux,getuy
      double precision vx,vy,v,gam
      DATA PI/3.141592653589793d0/

      vx=getux(x3,y3,z3,t3)
      vy=getuy(x3,y3,z3,t3)
      v=sqrt(vx**2+vy**2)
      gam=1.d0/sqrt(1-v**2)
	tau = sqrt(t3**2-z3**2)
	cosheta = t3/tau
      GETNEFF=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *GETTEMP(X3,Y3,Z3,T3)**3/PI**2
      getneff = getneff*gam !Local velocity contraction
	getneff = getneff/cosheta
      END
      
      

      DOUBLE PRECISION FUNCTION GETTEMP(X4,Y4,Z4,T4)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
      common/temperature/tempfac
      double precision tempfac
C--local variables
      DOUBLE PRECISION X4,Y4,Z4,T4,TAU,NPART,EPS0,EPSIN,TEMPIN,PI,
     &NTHICK,ys,MEDPART,interpolate
      double precision gettempmax
      DATA PI/3.141592653589793d0/

      GETTEMP=0.D0

      IF(ABS(Z4).GT.T4)RETURN

      TAU=SQRT(T4**2-Z4**2)
      if(tau.eq.0.d0) then
      return
      end if

      GETTEMP=tempfac*interpolate(X4,Y4,tau,1)
      if(gettemp.lt.tc) gettemp=0.0d0
      if(gettemp.ge.gettempmax()) gettemp=gettempmax()

      !write(*,*) "Get temp called:"
      !write(*,*) "(x,y,z,t)=",x4,y4,z4,t4,tau
      !write(*,*) "Temp=",gettemp
      !write(*,*) gettempmax()

      RETURN

      END

      double precision function getux(x,y,z,t)
            implicit none
            double precision x,y,z,t
            double precision tau,interpolate
            tau=sqrt(t**2-z**2)
            getux=interpolate(x,y,tau,2)
            return
      end function

      double precision function getuy(x,y,z,t)
            implicit none
            double precision x,y,z,t
            double precision tau,interpolate
            tau=sqrt(t**2-z**2)
            getuy=interpolate(x,y,tau,3)
            return
      end function


      DOUBLE PRECISION FUNCTION GETTEMPMAX()
      IMPLICIT NONE
C--medium parameters
      COMMON/TEMPMAX/TEMPMAXIMUM
      DOUBLE PRECISION TEMPMAXIMUM
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--function call
      DOUBLE PRECISION GETTEMP
      
      !GETTEMPMAX=GETTEMP(0.D0,0.D0,0.D0,TAUI)
      !write(*,*) "Max temp:", tempmaximum
      GETTEMPMAX=TEMPMAXIMUM
      END



      DOUBLE PRECISION FUNCTION GETMDMAX()
      IMPLICIT NONE
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION GETTEMPMAX
      GETMDMAX=MDSCALEFAC*3.*GETTEMPMAX()
      GETMDMAX=MAX(GETMDMAX,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMDMIN()
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC
      DOUBLE PRECISION GETTEMPMAX
	GETMDMIN=MDSCALEFAC*3.*TC
      GETMDMIN=MAX(GETMDMIN,MDFACTOR)
      END



      DOUBLE PRECISION FUNCTION GETMSMAX()
      IMPLICIT NONE
      DOUBLE PRECISION GETMDMAX,SQRT
      GETMSMAX=GETMDMAX()/SQRT(2.D0)
      END



	DOUBLE PRECISION FUNCTION GETNATMDMIN()
	IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--factor to vary Debye mass
	COMMON/MDFAC/MDFACTOR,MDSCALEFAC
	DOUBLE PRECISION MDFACTOR,MDSCALEFAC,PI
      DATA PI/3.141592653589793d0/
C--local variables
	DOUBLE PRECISION T,GETMDMIN
	T=GETMDMIN()/(MDSCALEFAC*3.)
      GETNATMDMIN=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *T**3/PI**2
	END



	DOUBLE PRECISION FUNCTION GETLTIMEMAX()
	IMPLICIT NONE
C--medium parameters
      COMMON/LTIME/MODLTIME
      DOUBLE PRECISION MODLTIME
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--function call
      DOUBLE PRECISION GETTEMPMAX
      !if(medfilelist.eqv..false.) then
      if(.false.) then
	  GETLTIMEMAX=TAUI*(GETTEMPMAX()/TC)**3*cosh(etamax2)
      else
      !Fabio: Putting my LTIME
      !write(*,*) "Lifetime whithout boost:",modltime
      GETLTIMEMAX=MODLTIME*cosh(etamax2)
      endif
	 END



      DOUBLE PRECISION FUNCTION GETNEFFMAX()
      IMPLICIT NONE
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
      common/temperature/tempfac
      double precision tempfac
C--max rapidity
	common/rapmax2/etamax2
	double precision etamax2
C--   local variables
      DOUBLE PRECISION PI,GETTEMPMAX
      DATA PI/3.141592653589793d0/
      GETNEFFMAX=(2.*6.*NF*D3*2./3. + 16.*ZETA3*3./2.)
     &     *GETTEMPMAX()**3/PI**2
      END
      
      

      DOUBLE PRECISION FUNCTION NPART(XX1,YY1,XX2,YY2)
      IMPLICIT NONE
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--local variables
      DOUBLE PRECISION XX1,YY1,XX2,YY2,NTHICK
      NPART = NTHICK(XX1,YY1)*(1.-EXP(-SIGMANN*NTHICK(XX2,YY2))) +
     &        NTHICK(XX2,YY2)*(1.-EXP(-SIGMANN*NTHICK(XX1,YY1)))

      END
      


      DOUBLE PRECISION FUNCTION NTHICK(X1,Y1)
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--identifier of log file
	common/logfile/logfid
	integer logfid
C--nuclear thickness function
      COMMON /THICKFNC/ RMAX,TA(100,2)
      DOUBLE PRECISION RMAX,TA
      INTEGER LINE,LMIN,LMAX,I
      DOUBLE PRECISION X1,Y1,XA(4),YA(4),Y,DY,R,C,B,DELTA
  
      R=SQRT(X1**2+Y1**2)
      IF(R.GT.TA(100,1))THEN
	 NTHICK=0.
      ELSE
	 LINE=INT(R*99.d0/TA(100,1)+1)
	 LMIN=MAX(LINE,1)
	 LMIN=MIN(LMIN,99)
	 IF((R.LT.TA(LMIN,1)).OR.(R.GT.TA(LMIN+1,1)))
     &	write(logfid,*)LINE,LMIN,R,TA(LMIN,1),TA(LMIN+1,1)
	 XA(1)=TA(LMIN,1)
	 XA(2)=TA(LMIN+1,1)
	 YA(1)=TA(LMIN,2)
	 YA(2)=TA(LMIN+1,2)
	 C=(YA(2)-YA(1))/(XA(2)-XA(1))
	 B=YA(1)-C*XA(1)
	 NTHICK=C*R+B
      ENDIF
      END



      SUBROUTINE CALCTA()
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--   nuclear thickness function
      COMMON /THICKFNC/ RMAX,TA(100,2)
      DOUBLE PRECISION RMAX,TA
C--variables for integration
      COMMON/INTEG/B,R
      DOUBLE PRECISION B,R
C--local variables
      INTEGER NSTEPS,I
      DOUBLE PRECISION EPS,HFIRST,Y

      NSTEPS=100
      EPS=1.E-4
      HFIRST=0.1D0

      R=1.12*A**(0.33333)-0.86*A**(-0.33333)
      RMAX=2.*R

      DO 10 I=1,NSTEPS
C--set transverse position
       B=(I-1)*2.D0*R/NSTEPS
       Y=0.D0
C--integrate along longitudinal line
       CALL ODEINT(Y,-2*R,2*R,EPS,HFIRST,0.d0,101)
       TA(I,1)=B
       TA(I,2)=Y
 10   CONTINUE
      END



      SUBROUTINE CALCXSECTION()
      IMPLICIT NONE
C--medium parameters
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--   geometrical cross section
      COMMON /CROSSSEC/ IMPMAX,CROSS(200,3)
      DOUBLE PRECISION IMPMAX,CROSS
C--local variables
      INTEGER IX,IY,IB
      DOUBLE PRECISION B,P,PROD,X,Y,NTHICK,NPART,pprev
      pprev=0.d0
      DO 30 IB=1,200
       B=0.1d0*IB
       PROD=1.d0
       DO 10 IX=1,100
        DO 20 IY=1,100
         X=-20.d0+IX*0.4d0
         Y=-20.d0+IY*0.4d0
         PROD=PROD*
     &EXP(-NTHICK(X+B/2.D0,Y)*SIGMANN)**(0.16d0*NTHICK(X-B/2.D0,Y))
 20     CONTINUE
 10    CONTINUE
       P=(1.D0-PROD)*8.8D0/14.D0*B
       CROSS(IB,1)=B
       CROSS(IB,2)=P
       if (ib.eq.1) then
        cross(ib,3)=0.d0
       else
        cross(ib,3)=cross(ib-1,3)+(p+pprev)/2.*0.1
       endif
       pprev=p
 30   CONTINUE
      IMPMAX=19.95d0
      END

      SUBROUTINE MYMED()
      IMPLICIT NONE
      COMMON/LTIME/MODLTIME
      DOUBLE PRECISION MODLTIME
      COMMON/TEMPMAX/TEMPMAXIMUM
      DOUBLE PRECISION TEMPMAXIMUM
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      COMMON/logfile/logfid
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
	common/rapmax2/etamax2
	double precision etamax2
      COMMON/MEDFILEC/MEDFILE,NLIST,endoff
      CHARACTER*200 MEDFILE
      INTEGER NLIST
      LOGICAL ENDOFF
      INTEGER logfid
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      double precision gridx(834,834),gridy(834,834)
      CHARACTER DUMMIE,CONTENT*100
      DOUBLE PRECISION tempsum,entropy
      INTEGER I,J,K,POS,II,kk,kkk,length
      logical ltime
      double precision hightemp
C--grid parameters
      common/gridpar/ gdt,gdx,gxmax,gxmin,gnx,gny,gnt
      double precision gdt,gdx,gxmax,gxmin,s
      integer gnx,gny,gnt,probcounter

      NX=834
      dx=50.d0/834.d0

      hightemp=0.d0


      do i=1,834
      do j=1,834
      do k=1,60
      tprofile(i,j,k)=0.d0
      end do
      gridx(i,j)=-25.d0+dx*(i-1)
      gridy(i,j)=-25.d0+dx*(j-1)
      end do
      end do
      
      call reader(medfile,834,60,timesteps,tprofile,ux,uy)
C--Loop that finds the medium lifetime and also
C--its evolution in entropy and temperature
C--as well as its highest temperature

      s=0.d0
      prob(1)=0.d0
      probcounter=2
      do k=1,60
      
      ltime=.true.
      if(.true.) then
      
      entropy=0.d0
      tempsum=0.d0
      do kk=1,NX
            do kkk=1,NX
            entropy=entropy+tprofile(kk,kkk,k)**(3.d0/4.d0)*dx**2
            tempsum=tempsum+tprofile(kk,kkk,k)**(4.d0)*dx**2
            if(k.eq.2.and.probcounter.le.834**2) then
            !if(mod(probcounter,834).eq.2) then
            !write(*,*) "Line read:",probcounter,prob(probcounter-1)
            !end if
            if(tprofile(kk,kkk,k).gt.1.0d0*tc) then
            s=s+tprofile(kk,kkk,k)**(3.d0)
            !prob(probcounter)=prob(probcounter-1)+tprofile(kk,kkk,k)
            end if
            prob(probcounter)=s
            probcounter=probcounter+1
            end if
            if(tprofile(kk,kkk,k).ge.tc) then
            ltime=.false.
            end if
            if (tprofile(kk,kkk,k) .gt. hightemp) then
                  hightemp=tprofile(kk,kkk,k)
            end if
            enddo
      enddo
      !write(*,*) "k-1=",K-1
      !write(*,*) "Avg temp=",tempsum
      !write(*,*) "Entropy=",entropy
      tempmaximum=hightemp
      if(k.eq.2) then
      do kk=1,834**2
      prob(kk)=prob(kk)/s
      end do
      end if
      
      endif

      write(*,*) "Ltime:",ltime
      if(.not.ltime.and.timesteps(k+1).gt.timesteps(k)) then
      write(*,*) "Lifetime not reached:",timesteps(k)
      nt=k+1
      modltime=timesteps(k+1)
      else
      write(*,*) "Lifetime reached:",timesteps(k)
      endif

      end do
      write(*,*) "LTIME= ",modltime
      write(*,*) "Highest temp=",hightemp

      WRITE(*,*) "Temperature profile read succesfully :)"

      END

      DOUBLE PRECISION FUNCTION MEDPART(X4,Y4,Z4,T4)
      IMPLICIT NONE
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
      COMMON/MEDPARAM/CENTRMIN,CENTRMAX,BREAL,CENTR,RAU,
     & NX,NY,NT,NF,DX,DT,XMAX,XMIN,TMAX
      INTEGER NX,NY,NT,NF
      DOUBLE PRECISION DX,DT,XMAX,XMIN,TMAX
      DOUBLE PRECISION CENTRMIN,CENTRMAX,BREAL,CENTR,RAU
      common/grid/timesteps(60),tprofile(834,834,60),prob(834**2)
      double precision timesteps,tprofile,prob
      common/gridvel/ux(834,834,60),uy(834,834,60)
      double precision ux,uy
      DOUBLE PRECISION X4,Y4,Z4,T4
      DOUBLE PRECISION STEP
      DOUBLE PRECISION TAU,interpolate
      STEP=(XMAX-XMIN)/(NX-1)
      TAU=SQRT(T4**2-Z4)
      TAU=0.0d0
      MEDPART=interpolate(X4,Y4,tau,1)
      END

      DOUBLE PRECISION FUNCTION MEDDERIV(XVAL,W)
      IMPLICIT NONE
      DOUBLE PRECISION XVAL
      INTEGER W
C--medium parameters
      COMMON/MEDPARAMINT/TAUI,TI,TC,D3,ZETA3,D,
     &N0,SIGMANN,A,WOODSSAXON,MODMED,MEDFILELIST
      DOUBLE PRECISION TAUI,TI,TC,ALPHA,BETA,GAMMA,D3,ZETA3,D,N0,
     &SIGMANN
      INTEGER A
      LOGICAL WOODSSAXON,MODMED,MEDFILELIST
C--variables for integration
      COMMON/INTEG/B,R
      DOUBLE PRECISION B,R

      IF (W.EQ.1) THEN
C--XVAL corresponds to z-coordinate
       MEDDERIV=N0/(1+EXP((SQRT(B**2+XVAL**2)-R)/D))
      ELSE 
       MEDDERIV=0.D0
      ENDIF
      END

      DOUBLE PRECISION FUNCTION INTEGRATE(TPROF)
      IMPLICIT NONE
      DOUBLE PRECISION TPROF,XMAX,XMIN,YMAX,YMIN,STEP
      DOUBLE PRECISION TERM
      DOUBLE PRECISION INTERMED(100)
      INTEGER N,I,J
      XMAX=10.d0
      XMIN=-10.d0
      YMAX=10.d0
      YMIN=-10.d0
      N=FLOOR(MAX((XMAX-XMIN+2*STEP)/(2*STEP),0.d0))
      INTEGRATE=0.d0
      DO 10 I=1,N
          TERM=0.d0
          INTERMED(I)=0.d0
          DO 20 J=1,N
              TERM=(STEP/3)*(TPROF(I,J)+TPROF(I,J+2)
     &        +4*TPROF(I,J+1))
              INTERMED(I)=INTERMED(I)+TERM
20        CONTINUE
10    CONTINUE
      TERM=0.d0
      DO 30 I=1,N
          TERM=(STEP/3)*(INTERMED(I)+INTERMED(I+2)
     &    +4*INTERMED(I+1))
          INTEGRATE=INTEGRATE+TERM
30    CONTINUE
      END
