!*********************************************************************************************************
      SUBROUTINE HEATR(lwind,sum,latit,par,cloud,qsn,qac,qw,xev,xco,tsup,LAKE_EVA,LAKE_PREC,jyear)     
!***************************************************************************************************
! 
!       performs the thermal transfers at the surface
!       calculates the time step IF running on a daily time step
!      light has been modified to sinusoidal distribution
!      the day commences at midday - instead of a constant light
!      field lasting for 12 hours and commencing at 6am. the new
!      version also uses a variable daylength, determined from the
!      time of the year. Written by Joaquim P. Losada and Geoff Schladow
!      (some of the turbulent calculations come from the McCord's work)
!     					UC-Davis October 2000
! 
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

		REAL(8), PARAMETER :: zero = 0.0d0,one = 1.0d0,two = 2.0d0
		REAL(8), PARAMETER :: four = 4.0d0,six = 6.0d0,ten = 10.0d0
		REAL(8), PARAMETER :: pt5  = 0.5d0,pt6 = 0.6d0,pt9 = 0.9d0
		REAL(8), PARAMETER :: pt17 = 0.17d0
		REAL(8), PARAMETER :: arfac  = 1.0d+6,thsnd = 1000.0d0
		REAL(8), PARAMETER :: qhrmrn = 48.0d0
		REAL(8), PARAMETER :: secqhr = 900.0d0
		REAL(8), PARAMETER :: clwout = 5.67d-8
		REAL(8), PARAMETER :: clwin = 5.313d-13 
		REAL(8), PARAMETER :: dctodk = 273.15d0
		REAL(8), PARAMETER :: secday = 86400.0d0
		REAL(8), PARAMETER :: htevap = 2.453d+9,spheat = 4181.60000
		REAL(8), PARAMETER :: exchk  = 20.0d0
		REAL(8)	::	exchk1
		REAL(8), PARAMETER :: pi=3.1415926535897932384626433832795d0
      REAL(8)  evaporation, precipitation
      REAL(8)  LAKE_EVA,LAKE_PREC, AP_P0
!     ***** new variables added for water quality *****

!		albed0=% photosynthetically active radiation scattered previous day
!		albed1=% photosynthetically active radiation scattered present day
!		et1 = extinction coefficient
!		etca = specific attenuation (m**2 mg**-1) by chla
!		etwat = attenuation coefficient of water (m**-1), inorganic particles etc.
!		expon = extinction coefficient * depth
!		latit = latitude of reservoir
!     MCCORD = controls whether bulk of turbulent Evaporative fluxes is enabled
!		par = photosynthetically active radiation in watts for water surface
!		qsn = solar radiation at water surface
!		qsw = array of solar radiation over water layer depths
!		solar1 = solar angle today
!		solary = sun's arc from sunrise to sunset for today
!		ustar = wind speed
!		mdays = days since start of year for today

		REAL(8) ::  add
		REAL(8) ::  cloud
		REAL(8) ::  c_drag
		REAL(8) ::  dadd
		REAL(8) ::  del(maxns)
!IF		REAL(8) ::  densty
		REAL(8) ::  dt
		REAL(8) ::  expon(maxns)
		REAL(8) ::  heat(maxns)
		REAL(8) ::  latit
		REAL(8) ::  par
		REAL(8) ::  qac       ! Long Wave flux from sky {W/m2}
		REAL(8) ::  qsn       ! Short Wave flux {W/m2}
		REAL(8) ::  qsw(maxns)
		REAL(8) ::  qw		 ! Long Wave flux emitted by water surface {W/m2}
		REAL(8) ::  sum		 ! Sumation of fluxes
!IF		REAL(8) ::  satvap
		REAL(8) ::  solar1
		REAL(8) ::  solary
		REAL(8) ::  svp0
		REAL(8) ::  sw_net   ! Short wave after albedo
		REAL(8) ::  t0		! beginning of time step {s, radians}
		REAL(8) ::  t1,t3	! END of time step {s, radians}	
		REAL(8) ::  ustar
		REAL(8) ::  vold
		REAL(8) ::  wlost	! water lost to evaporation {m/s}
		REAL(8) ::  wmass
		REAL(8) ::  xco		! conductive heat gain {W/m2}
		REAL(8) ::  xev		! latent (evaporative) heat gain {W/m2}
		REAL(8) ::  xlw      ! Net atmospheric long wave {W/m2}
		REAL(8) ::  RHtoSVPD ! Saturated Vapour Pressure
      REAL(8) ::  Eta_chloro, sumvol
!     Control variables for the selection algorithm
		INTEGER*4 i, j
		INTEGER*4 icode
		INTEGER*4 layno
		INTEGER*4 lwind
		INTEGER*4 mdays,jyear
		INTEGER*4 ndlta
		INTEGER*4 ucheck
		INTEGER*4 lh	

!	logical*4 pheatr
   
      REAL(8) ::  emissivity,albedol,albedo,tsup, daily_sun_angle, albedo_sun_angle
		REAL(8) ::  A1,B11,B2,B3,B4
!		REAL*8, PARAMETER  :: rain_f =(1.0-0.36)     !equivalent rain all over the lake 
!      REAL*8, PARAMETER  :: rain_f=(1.0) 
!      REAL*8, PARAMETER  :: rain_f=(1.0+0.08) 
      REAL*8, PARAMETER  :: rain_f=(1.0) !SCT reset to zero for external adjustment	 7/24/2019  
!      Hardwired de emissivity and Albedo long wave of water
! 
		emissivity = 0.97d0
		Albedol    = 0.03d0
		exchk1= -exchk
!		IF(JYEAR.EQ.1995.OR.JYEAR.EQ.1996.OR.JYEAR.EQ.1998.OR.JYEAR.EQ.2005)THEN
!		   rain_f =(1.0-0.45)
!		ELSE
!		   rain_f =(1.0-0.36)
!		ENDIF
! ***************************************************************
!      Set the attenuation coefficient
! 
                
    		mdays  = jday-(jday/1000)*1000
!        sumvol=0.0
!        Do i=1,ns
!            sumvol=sumvol+vol(i)
!            print*, depth(i), area(i),vol1(i)
!        ENDDO
!        print*,mdays, iclock, sumvol
       
!        pause
!2015/9		DO i = 1,ns
!2015/9			et1(i) = 0.0d0
!2015/9			DO j = 1,nchl
!2015/9				et1(i) = et1(i)+etca(j)*wqual(i,j)
!2015/9			ENDDO ! j
!2015/9			DO j = 1,7
!2015/9				et1(i) = et1(i)+etpart(j)*cf(i,j)
!2015/9			ENDDO ! j
!2015/9			et1(i) = et1(i)+etwat
!2015/9        ENDDO   ! i
    DO i =1,ns
!         if(mdays.ge.90.and.mdays.lt.135) then
!			   et1(i)=1.70d0   !eta=1.7/secchi depth (1 to 1.5 SD)
!		  elseif(mdays.ge.135.and.mdays.lt.155) then
!			   et1(i)=1.8d0   !eta=1.7/secchi depth (1 to 1.5 SD)
!	     elseif(mdays.ge.155.and.mdays.lt.210) then    ! 160 to 210
!			   et1(i)=1.9d0   !eta=1.7/secchi depth (1 to 1.5 SD)
!		  elseif(mdays.ge.210.and.mdays.lt.300) then    ! 210 to 250
!			   et1(i)=1.9d0   !eta=1.7/secchi depth (1 to 1.5 SD)     
!	     elseif(mdays.ge.300.and.mdays.lt.330) then
!			   et1(i)=1.8d0   !eta=1.7/secchi depth (1 to 1.5 SD)
!        elseif(mdays.ge.330.and.mdays.lt.366) then
!	         et1(i)=1.7d0   !eta=1.7/secchi depth (1 to 1.5 SD)
!        endif
       et1(i)=1.70d0/msecchi_dep
        ! if(depth(i).ge. (depth(ns)-1.00)) then
            ! et1(i)=1.70d0/msecchi_dep
        ! elseif (depth(i).gt. (depth(ns)-1.00).and.depth(i).le.(depth(ns)-11.00)) then
            ! et1(i)=1.70d0/msecchi_dep !+1.7*((depth(ns)-depth(i))/5.0)/msecchi_dep
         ! elseif (depth(i).gt. (depth(ns)-11.00).and.depth(i).le.(depth(ns)-31.00)) then
            ! et1(i)=1.70d0/msecchi_dep !+1.7*((depth(ns)-depth(i))/5.0)/msecchi_dep    
        ! elseif (depth(i).gt. (depth(ns)-31.00).and.depth(i).lt.(depth(ns)-70.00)) then
            ! et1(i)=1.70d0/msecchi_dep !+1.7*((depth(ns)-depth(i))/5.0)/msecchi_dep
        ! elseif (depth(i).ge.(depth(ns)-70.00)) then
            ! et1(i)=1.70d0/msecchi_dep !+1.70*((depth(ns)-depth(i))/5.0)/msecchi_dep
        ! endif
 !       print*, i,ns, depth(i), temp(i) 
 !       pause        
    ENDDO
	
 !     PRINT*,MDAYS,nosecs,et1(ns),sw
 !     pause
! ***********************************************************
!      CALL Light_Tahoe(latit)
! ***********************************************************	
!		print*,temp(ns),temp(1)
!		pause
 
		DO i = 1,ns
!			IF (i.eq.ns) THEN
!				et1(i)=2.0*et1(i)
!			ENDIF
      
			IF(i.eq.1)THEN
				del(i) = depth(i)
			ELSE
				del(i) = depth(i)-depth(i-1)
			ENDIF
			expon(i) = -et1(i)*del(i)		
!  if(i.eq.ns)   print*,'gbs1', expon(ns),del(i)
! 
!      Determine the reflection at the water surface and subtract
!      from the incoming short radiation (W/m2)
! 

	  albedo = 0.03d0
      !albedo = 0.08d0 + 0.02d0*sin(2.0d0*pi*mdays/365+(pi/2.0d0)) !gbs summer high, winter low
	  
      sw_net = sw*(one - albedo)
	  


!		ustar = sqrt(1.612e-6*u6x*u6x)
! 
!      New Linear definition of Cd (Smith and Banke formulae, from Henderson-Sellers
! 
		c_drag = 1.3d-3 
!     c_drag = (0.61 + 0.063 * u6x)/ 1000.
!		ucheck = 1000.0 * u6x	
!		SELECT CASE (ucheck)
!			CASE (:300)
!				c_drag = .0015			
!			CASE (301:2200)
!				c_drag = (1.08 * u6x**(-0.15))/1000.0		
!			CASE (2201:5000)
!				c_drag = (0.771 + 0.0858 * u6x)/1000.0		
!			CASE (5001:8000)
!				c_drag= (0.867 + 0.0667 * u6x)/1000.0		
!			CASE (8001:25000)
!				c_drag= (1.2 + 0.025 * u6x)/1000.0	
!			CASE (25001:50000)
!				c_drag= (0.073 * u6x)/1000.	
!			CASE default
!				c_drag= 0.0037	
!		END SELECT

		ustar = sqrt(c_drag*(1.24d0/1000.0d0)*u6x*u6x)
		

! 	DO all the bulk transfer calcs outside of the time step iteration
! 	as these fluxes DO not change during the time step, they are functions
! 	only of properties at the END of the previous time step.

! 	evaporative heat flux affects top layer only
! 	this formula is based on the marciano-harbeck formulae for l. hefner
! 	see tva report p4.11. also p 4.17. also the spigel document.
! 	units are joules/m**2/sec

! 	conductive heat gain only affects top layer
! 	this comes from the bowen ratio; see the spigel document
! 	units are joules/m**2/s

! 	note: htevap = water density*latent heat
		ENDDO
		


		IF(bulk.eq. 1) THEN
			IF(humidity .eq. 1) THEN
!      Liu, Katsaros and Businger's algorithm. For corrected neutral conditions 
				CALL bulktrans(ustar,svp0,xev,xco,wlost)
				GOTO 999 	
			ELSE
				WRITE(*,*) 'IF bulk =1, THEN humidity =1'
				STOP
			ENDIF
		ELSE
! 
!      Bulk transfer eqs
! 
!			IF(humidity.eq.1) THEN	  
!				svpd =  VAPORTVA(T4,rh)
!			ENDIF

!	CALL bulktrans(ustar,svp0,xev,xco,wlost)
       CALL heat_evaporation(svpd,xev,xco,wlost)
!       heat_evaporation(SVPW,HEATE,HEATCT,Ewat)
!			svp0 =  satvap(temp(ns))
		IF(HUMIDITY.EQ.1) THEN
			SVPD = (rh)*(10.0**(9.286-(2322.38/(T4+273.15))))
		ENDIF
			svp0 =10.0**(9.286-(2322.38/(temp(ns)+273.15)))

      ENDIF
 999  CONTINUE
 
			xev  = 3.9d0*u6x*(svpd-svp0)	
			xco  = 1.0*2.536d0*u6x*(t4-temp(ns))		  
!       WRITE(*,fmt='(i5,5f10.2)')jday,xev,xco, T4,temp(ns), U6x
 		IF (xev.gt.0.0) xev = 0.0d0										!SCT 6/5/2019 Handeled in Heat_evaporation
		wlost = xev/htevap/1.35	!SCT 6/5/2019 Handeled in Heat_evaporation
!     print*,wlost                 
!   
!     	long wave emission (ie. longwave out) affects only top layer
! 
		qw = -clwout*emissivity*((dctodk+temp(ns))**four)
! 
!     	long wave absorption (i.e. longwave in) affects only top layer
!     	directly measured longwave absorption cloud = 1, incident longwave = 2, net lw = 3
! 

    	IF(lwind.eq.1)THEN
			cloud = srat
			qac = (1.0d0-albedol)*clwin*(one+pt17*cloud*cloud)*((dctodk+t4)**six)
			xlw = qw+qac
		ELSEIF(lwind.eq.2)THEN	 
			qac = srat*(1.0d0-albedol)
			xlw = qw+qac
		ELSEIF(lwind.eq.3)THEN
			xlw = srat
		ENDIF
! 
!      Set initial time step 
! 
! 		IF(itimes.eq.1440) THEN   !variable time step
! 			nosecs=900
!     ELSE
!2015/9			nosecs=itimes*60         !fixed time step
!     ENDIF
5		CONTINUE

      IF(MET_TIMESTEP.LT.86400) THEN          
!          daily_sun_angle=15.0*(12.0-float(iclock)/float(nosecs))
!          if(daily_sun_angle.ge.90.0.or.daily_sun_angle.lt.-90.0)daily_sun_angle=90.00
!          albedo_sun_angle=1.00-cos(daily_sun_angle/180.0)
!          qsn=(1.00-albedo_sun_angle)*sw_net*1000.0d0/ real(MET_TIMESTEP)	! convert from kJ hr^-1 m^-2 to W m^-2
!         qsn=(1.0-0.060)*sw_net*1000.0d0/ real(MET_TIMESTEP)	! convert from kJ hr^-1 m^-2 to W m^-2   gbs:2016/01/19 st:2018/08/02 albedo already applied
          qsn=sw_net*1000.0d0/ real(MET_TIMESTEP)  !KJ/hr-m2 to W/m2
		  !xev=xev*1000.0d0/ real(MET_TIMESTEP)		!KJ/hr-m2 to W/m2 - Moved from below Heat Evap SCT 6/5/2019 
          !xco=xco*1000.0d0/ real(MET_TIMESTEP) 	    !KJ/hr-m2 to W/m2 - Moved from below Heat Evap SCT 6/5/2019 
		 
         qsw(ns)  = qsn
		   SELECT CASE(ns)
			CASE(1)
				WRITE(*,*)'Warning NS = 1, THEN heat balance not well defined'
				heat(ns)  = qsw(ns)*area(ns)*arfac
			CASE default	
				qsw(ns-1) = qsn*exp(expon(ns))	! Beer's Law
				heat(ns)  = (qsw(ns)*area(ns) - qsw(ns-1)*area(ns-1))*arfac
		   END SELECT

		   heat(ns) = heat(ns) + (xev+xco+xlw)*area(ns)*arfac
!		   dt = heat(ns)/(spheat*(thsnd+den(ns))*vol(ns)*thsnd)
		   GOTO 25
      ELSEIF(MET_TIMESTEP.EQ.86400)  THEN 
		    t0 = iclock
		    t1 = iclock + nosecs
! 
!         Forces the time step to END at the output time
! 
		   IF(t0.lt.itmpr.and.t1.ge.itmpr) THEN	 
			   nosecs = itmpr-t0
			   t1 = iclock + nosecs	 	
		   ENDIF

		   IF(t1.gt.secday) THEN
			   nosecs=secday-t0
			   t1 = iclock + nosecs
		   ENDIF      
! 
!      here IF no daylight	  
! 
		   IF (t0 .ge. sunset .or. t1 .le. sunrise) THEN
		      qsn = 0.0d0
		   ELSE			
! 
!      here IF any daylight at all. The 1000.0 is     
 !      to convert Kj/m^2 to j/m^2
 ! 
 			   t0 = max(t0,sunrise) - sunrise
 			   t1 = min(t1,sunset)  - sunrise
			   sw_net = sw_net * 1000.0d0 / (t1-t0)
 			   t0 = t0 / daylength * pi	! convert to radians
 			   t1 = t1 / daylength * pi
 			   qsn = sw_net * (cos(t0) - cos(t1)) / 2.0d0
		   ENDIF 
!		   qsn		= qsn / MET_TIMESTEP	! convert from J m^-2 to W m^-2
		   qsw(ns) = qsn		   		    

! 		IF(mod(iclock,14400).eq.0) THEN
! 			WRITE(121,FMT='(2I8,7F10.1)')jday,iclock,qsn,qac,qw,xev,xco,t4,temp(ns) 
! 		ENDIF
																	
!gbs          WRITE(122,12) jday, nosecs,qsn,sw_net,t0,t1
!gbs       12 format(i3,4x,i8,4f10.2)
! 	into layer ns goes qsw(ns)-qsw(ns-1) over the area common to layers
! 	ns and ns-1 and qsw(ns) over the rest of area(ns)
! 	units of heat(ns) are joules/sec; units of area(ns) are 10**6 m**2
!      IF ns = 1, qsw is not defined.

         SELECT CASE(ns)
			   CASE(1)
				   WRITE(*,*)'Warning NS = 1, THEN heat balance not well defined'
				   heat(ns)  = qsw(ns)*area(ns)*arfac
			   CASE default	
					   qsw(ns-1) = qsn*exp(expon(ns))	! Beer's Law
				   heat(ns)  = (qsw(ns)*area(ns) - qsw(ns-1)*area(ns-1))*arfac
		   END SELECT
   	
		   heat(ns) = heat(ns) + (xev+xco+xlw)*area(ns)*arfac           
               
   !     	IF time step ends at output times exit loop
   ! 
		   IF(iclock + nosecs.eq.itmpr) GOTO 20
   ! 
   !     IF set fixed time step, exit loop
   ! 
		   IF(itimes.ne.1440) GOTO 20	

   !     a 0.4 deg c maximum surface temperature change or wind strength
   !     relative to the depth of the mixed layer (htsave) limits the increments.
   !     This doesn't apply for a set time step. 

   !      convert heat inputs to temp increments ie. sp heat = 4.186*1000 j/deg/kg
   !      dt is the temperature increase per second/vol
   ! 
   ! 		modified 9 Jun 97 by bss - changed units of spheat to j kg^-1 C^-1
   ! 		i.e. spheat = 4186.
   ! 		replaced arfac with thsnd so that the calculation is done 'properly'
   ! 		rather than relying on offsetting errors
   ! 
		   dt = heat(ns)*nosecs/(spheat*(thsnd+den(ns))*vol(ns)*thsnd)

   !     avoid a division by zero  error IF ustar is zero
         IF(ustar.eq.0.0d0)ustar = 0.05d0

         IF(htsave.eq.0.0d0) THEN
			   ndlta = 86400
		   ELSE
			   ndlta=(htsave/ten)/(ustar*ustar)
			   IF(ndlta.lt.900) ndlta = 900
		   ENDIF

   !     Temperature	condition: exit IF surface layer heates too much	
         IF(int(0.4d0/abs(dt)).le.nosecs.or.iclock+nosecs.eq.86400) GOTO 20

   !     Wind condition: exit IF surface layer accelerates too much
		   IF(ndlta .le. nosecs) GOTO 20 

   ! 	Add 15 minutes since both dt and ndlta are less than the thresholds
   ! 	This is the only place nosecs is defined, so an integral multiple of
   !     15 minutes is forced on the internal time step 
         nosecs = nosecs + 900
         GOTO 5
      ENDIF 
      
15    CONTINUE
20    CONTINUE
25    CONTINUE
! 	   
!      photosynthetically active radiation (PAR) is approx. 45% of shortwave
! 
		par = 0.45d0*qsn
		DO 30 i = ns-1,1,-1
!   defeat floating underflows for deep lakes (setting qsw(i+1).lt.1.0e-305 and exchk1) 
			IF(expon(i).lt.exchk1)THEN 
				qsw(i) = zero
			ELSEIF(qsw(i+1).lt.1.0e-20) THEN  !-290
				qsw(i) = zero
			ELSE
				qsw(i) = qsw(i+1)*exp(expon(i))
			ENDIF
30		CONTINUE

!   into layer i goes qsw(i)-qsw(i-1) over the area common to layers
!   i and i-1 and qsw(i) over the rest of area(i)
!   units of heat(i) are joules/sec; units of area(i) are 10**6 m**2

		DO i = ns-1,2,-1        
			heat(i) = (area(i)*qsw(i)-area(i-1)*qsw(i-1))*arfac
		ENDDO
		heat(1) = qsw(1)*area(1)*arfac
      
!   compute the temperature increase in non-surface layers over nosecs
!mb	ammended to avoid underflow error in below when heat(i) is very
! 	small 13/5/93
		  ! write(*,fmt='(i7,5f10.2,1f10.8),3f10.2')jday,xev,xco,xlw,qsw(ns)-qsw(ns-1),xev+xco+xlw+qsw(ns)-qsw(ns-1),heat(i)*float(nosecs)/(spheat*(thsnd+den(i))*vol(i)*thsnd),float(nosecs),T4,temp(ns) 	 
		DO i = 1,ns
			IF (abs(heat(i)).ge.1d-20) THEN
!                dt = heat(i)*float(MET_TIMESTEP)/(spheat*(thsnd+den(i))*vol(i)*thsnd) 
                dt = heat(i)*float(nosecs)/(spheat*(thsnd+den(i))*vol(i)*thsnd) 
			
			ELSE
				dt = 0.0d0
			ENDIF
			temp(i) = temp(i)+dt
            !if (temp(i).lt.4.00)temp(i)=4.0 ! temperature threshold reduced from 6 to 3 sct
			!if (i.eq.ns)write(*,fmt='(i7,6f10.2)')jday,xev,xco,xlw,dt,T4,temp(ns)
			!if (i.eq.(ns-1))write(*,fmt='(i7,3f10.2)')jday,dt,T4,temp(ns)
		ENDDO
! 
!       adjust surface level for evaporation and rainfall (mm)
!       Convert mm_day to m_Deltatime
! ********************************************************************* 
		LAKE_EVA =(-wlost)*area(ns)*1.0d6*float(nosecs)
		LAKE_PREC=(rain_f*rain)*area(ns)*1.0d6*float(nosecs)/(1000.0d0*float(MET_TIMESTEP))
		IF(iclock.eq.0.0) THEN
			evaporation  =0.0
			precipitation=0.0
			evaporation=evaporation+LAKE_EVA
			precipitation=precipitation+LAKE_PREC
		ELSE
			evaporation=evaporation+LAKE_EVA
			precipitation=precipitation+LAKE_PREC
		ENDIF
! 	IF(iclock.eq.82800) THEN	
! 		WRITE(20, FMT='(I8,4f20.3,i10)')jday,evaporation,
!      &                      precipitation,area(ns)*1.0d6, rain,nosecs
!       ENDIF
!       WRITE(20, FMT='(2I8,3f20.8)')jday,iclock, wlost*float(nosecs),
!     			(rain_f*rain/1000.0d0)*(float(nosecs)/86400.0d0),rain_f*rain      
! **********************************************************************
! 
		add = wlost*float(nosecs)+(rain_f*rain/1000.0d0)*(float(nosecs)/float(MET_TIMESTEP))	
		depth(ns) = depth(ns)+add

!      print*,teminf(1)
!        if(temp(ns).lt.max(teminf(2), teminf(2)))  temp(ns)=max(teminf(2), teminf(2))        
!        temp(ns)=max(teminf(2), teminf(1))

		vold  = vol(ns)
		icode = 1     !ARRAYS OF VOLUME AND AREAS FROM DEPTHS (icode=1)!ARRAYS OF DEPTHS AND AREAS FROM VOLUME (icode=2)				   
		layno = ns

		CALL resint(icode,layno)
! 
!      Adjust Salt, WQ variables, Particles, of top layer NS
! 
		dadd    = densty(temp(ns),zero)
		wmass   = vold*(den(ns)-dadd) + vol(ns)*dadd + vol(ns)*thsnd
		sal(ns) = vold*thsnd*sal(ns)/wmass+vold*den(ns)*sal(ns)/wmass  

      DO j = 1,28
			wqual(ns,j) = wqual(ns,j)*vold/vol(ns)
		ENDDO

      DO j = 1,7
			cf(ns,j) = cf(ns,j)*vold/vol(ns)
		ENDDO
! 
!       recalculate densities
! 
		DO i = 1,ns
			den(i) = densty(temp(i),sal(i))
		ENDDO
! 
!       Print data for total period (W m**-2)
! 
      SUM = SUM + (XLW+XCO+XEV+QSN)
!      IF(iclock.ge.43200.and.iclock.lt.(43200+nosecs)) THEN


	  IF(HUMIDITY.EQ.1) THEN
        IF(PHEATR)WRITE(19,9000)jday,ICLOCK,XLW,XEV,XCO, QW,QAC,QSN,(XLW+QSN),SUM,t4,temp(ns), temp(1), &
                    (rh)*(10.0**(9.286-(2322.38/(T4+273.15)))),SVP0,-wlost*float(nosecs)*1000.0
	  ENDIF
	  IF(HUMIDITY.EQ.0) THEN
        IF(PHEATR)WRITE(19,9000)jday,ICLOCK,XLW,XEV,XCO, QW,QAC,QSN,(XLW+QSN),SUM,t4,temp(ns), temp(1), &
                    SVPD,SVP0,-wlost*float(nosecs)*1000.0
	  ENDIF
	  IF(HUMIDITY.EQ.1) THEN
		!SVPD=TVAvapor_temp_RH(T4,rh)	 Replaced by power equation below SCT 4/22/2019
			SVPD = (rh)*(10.0**(9.286-(2322.38/(T4+273.15))))
	  ENDIF
!      ENDIF
9000  FORMAT(I8,1X,I5,8f11.1,3f8.2, 2F8.2, f10.4)
! 
!      Simulated Surface temperature
! 
		tsup = temp(ns)
!gbs          close(122)
		RETURN
		END SUBROUTINE HEATR