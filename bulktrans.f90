!***********************************************************************************
		SUBROUTINE BULKTRANS(ustar,svp0,TurbLatFlux,Tflux,EVFLUX)
!**********************************************************************************
!	This has the corrected stflux and Liuflux algorithms
!	previous versions of stflux incorrectly used the air temperature rather
!	than the water surface temperature to define the specific humidity and
!	density of air immediately above the water
!
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

		REAL*8 Cpair	! specific heat capacity of air
		REAL*8 EVFXstill	! evaporative mass flux	into still air {kg/(m^2	s)}
		REAL*8 EVFLUX	! evaporative mass flux	{kg/(m^2 s)} => {m/s}
		REAL*8 fluxsens	! sensible heat flux into still air {W/m^2}
		REAL*8 htevap	! latent heat of vaporization
		REAL*8 TurbLatFlux	! turbulent latent heat flux [W/m^2]

     	REAL*8 Q		! observed specific humidity at height Qheight
			! LKB turbulent vapor flux
!	REAL*8 Qnet		! net heat flux into surface layer [W/m^2]
		REAL*8 Qsat		! saturation specific humidity at Tsurf 
		REAL*8 Qstar	! specific humidity fluctuation scale
!gbs17Jan06	REAL*8 zq		! height of humidity sensor above water surface
	
		REAL*8 Pair		! atmospheric pressure [mb] assumed constant at 1013 mb
		REAL*8 rhoair	! density of air at temp T4, rel hum RH
!IF		REAL*8 SATVAP	! function to calculate	saturated vapor pres.
		REAL*8 svp0		! saturated vapor pressure at Tsurf {mb}

		REAL*8 Tflux	! LKB turbulent sensible heat flux [W/m^2]
		REAL*8 Tstar	! temperature fluctuation scale
		REAL*8 Tsurf	! water surface temperature [C]
!gbs17Jan06	REAL*8 ztair	! height of temperature sensor above water surface

		REAL*8 ustar	! shear velocity in water
		REAL*8 UstarAir	! shear velocity in air
!gbs17Jan06	REAL*8 zu		! height of anemometer above water surface
		REAL*8 Qflux
!	PARAMETER(Pair=1013., htevap=2.445E+6, cpair=1010.0)
		PARAMETER(htevap=2.453d+9, cpair=1010.0,Pair=1013.0) 
! relhum = measured relative humidity [%]
! T4 = observed air temperature [C] at height Theight
! u6 = measured wind speed
! SVPD = vapor pressure of atmosphere

!
!	Calculate meteorological fluxes
!	get fluxes into still air and density of atmosphere
!
		Tsurf = temp(ns)

!	CALL stflux(Tair,Tsurf,RH,Pair,rhoair,SVPD,fluxsens,EVFXstill,q,qsat)
		CALL STFLUX(T4,Tsurf,Pair,rhoair,fluxsens,EVFXstill,q,qsat)

		EVFXstill = EVFXstill * (2500.0-2.39*temp(ns))*1000*(den(ns)+1000.) 
		SVP0 = SATVAP(Tsurf)
		EVFLUX = -1.712677e-3 * rhoair * u6 * 1.015161 * 0.662 * (SVP0-SVPD) / Pair * (2500.0-2.39*temp(ns))*1000*(den(ns)+1000.)

!gbs17Jan06	zu = 10.0	! anemometer always 2 m above water surface Bill say 2.0m
!gbs17Jan06	zq = 10.0	! air temp, rh measurement height  Bill say 1.7m
!gbs17Jan06	ztair = zq	! Vaisala combines air temp and humidity
		zu = 10.0
      zq = 10.0
		ztair = zq
!	CALL liuflux(UstarAir,qstar,tstar,Tair,Tsurf,Q,Qsat,wind,ztair,zq,zu,cd,ce,ch)
		CALL LIUFLUX(UstarAir,qstar,tstar,T4,Tsurf,Q,Qsat,u6,ztair,zq,zu)

!
!	Compute the fluxes, note that the - sign is removed from the calculation of the
!	sensible heat flux since a +ve flux to air (from LKB) is a heat loss from the 
!	water's point of view.
!
		Tflux = UstarAir * tstar * rhoair * Cpair
		Qflux = qstar * UstarAir * rhoair
		ustar = SQRT(UstarAir*UstarAir*rhoair/(den(ns)+1000.))
		TurbLatFlux = Qflux*(2500.0-2.39*temp(ns))*1000*(den(ns)+1000.)

		EVFLUX = EVFLUX / (den(ns)+1000.)	! convert to m/s from kg/m2-s

		END SUBROUTINE BULKTRANS	  

!	***********************************************************************
		SUBROUTINE STFLUX (Tair,Tsurf,Pair,RHOair,fluxsens,EVFLUX,q,qw)
!***************************************************************************
!	SUBROUTINE stflux calculates the latent	and sensible heat fluxes
!	and the	evaporative mass flux into still air.  Fluxes are positive
!	downwards, i.e.	into the water.	 A negative flux indicates heat
!	transfer from the water	to the air.
!
!	See Atmosphere-Ocean Dynamics by Adrian Gill p. 41
!
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

! constants
!
		REAL*8 cp		! heat capacity of air
		REAL*8 epsilon	! ratio of gas constants (dry air)/(water vapor)
		REAL*8 K		! thermal conductivity of air {W/(m C)}
		REAL*8 kappa	! thermal diffusivity of air {m^2/s}
		REAL*8 nu		! kinematic viscosity of air {m^2/s}
		REAL*8 Ra		! gas constant for dry air {J/kg-K}
		REAL*8 Rv		! gas constant for water vapor {J/kg-K}

		REAL*8 zero,one,two,three,hundrd,g,dCtodK,c1

		PARAMETER (zero=0.0,one=1.0,two=2.0,three=3.0,					&
     		c1=0.137,hundrd=100.0,g=9.81,dCtodK=273.15,					&
     		epsilon=0.62197,Ra=287.04,Rv=461.50,							&
     		K=2.54e-2,nu=1.50e-5,kappa=2.08e-5,cp=1012.0)

		REAL*8 alphasens	! sensible heat transfer coefficient {W/(m^2 C)}
		REAL*8 alphalat	! evap mass transfer coefficient {kg/(m^2 s)}
		REAL*8 deltaRHO	! density change from surface to atmosphere
		REAL*8 EVFLUX	! evaporative mass flux	{kg/(m^2 s)}
		REAL*8 fluxsens	! sensible heat flux into still air {W/m^2}
		REAL*8 Pair		! atmospheric pressure {mb}
		REAL*8 q		! observed specific humidity
		REAL*8 qw		! specific humidity for	saturated air
		REAL*8 r		! observed mixing ratio
		REAL*8 rw		! saturated mixing ratio for moist air
!IF		REAL*8 RH		! relative humidity (%)
!IF		REAL*8 RelHum	! relative humidity {decimal}
		REAL*8 RHOair	! ambient air density {kg/m^3}
		REAL*8 RHOsurf	! density of saturated air @ surface {kg/m^3}
		REAL*8 sign		! sign of cube root
!IF		REAL*8 SATVAP	! function to calculate	saturated vapor pres.
		REAL*8 SVPair	! saturated vapor pressure at Tair {mb}
!IF		REAL*8 SVPD		! ambient vapor pressure of atmosphere {mb}
		REAL*8 Tair		! air temperature (C)
		REAL*8 Tsurf	! temperature of water surface (C)
		
!
!	added 98.03.24 because of discrepancy between Liuflux and DYRESM bulk formula
!	DYRESM assumes dQ = q - qw(Tsurf)
!
!	stflux defines qw(Tair) which is passed as qsat to Liuflux
!	The comment in Liuflux implies it expects qw(Tsurf) - check it out
!
!	26 mar 98 - checking with Liu, Katsaros and Businger confirmed that Liuflux 
!	requires qsat = qw(Tsurf). All earlier versions of Liuflux that used stflux
!	to pass qw(Tair) rather than qw(Tsurf) are WRONG.
!
!	calculate specific humidities and densities from observed data
!
!	First compute the density of the air rho(Tair,RH)
!
		RelHum=RH
!	RelHum=RH/hundrd
		SVPair=SATVAP(Tair)
		rw=SVPair/(Pair-SVPair)*epsilon
		r=RelHum*rw
		qw=rw/(one+rw)
		q=r/(one+r)
		RHOair=Pair*hundrd/(Ra*(Tair+dCtodK)*(one+(-one+one/epsilon)*q))
!
!	Now compute conditions just above the water surface, i.e. water saturated air
!	at temperature Tsurf rho(Tsurf,100%)
!
		SVPair=SATVAP(Tsurf)
		rw=SVPair/(Pair-SVPair)*epsilon
		qw=rw/(one+rw)
		RHOsurf = Pair*hundrd/(Ra*(Tsurf+dCtodK)*(one+(-one+one/epsilon)*qw))
!
!	calculate atmospheric vapor pressure for evaporative fluxes
!	into moving air.
!
		SVPD=RHOair*q*Rv*(Tair+dCtodK)/hundrd
!
!	calculate sensible heat	flux into still air
!	deltaRHO > 0 means heat	is lost to the atmosphere
!
		deltaRHO=RHOair-RHOsurf
!
!	handle the possiblity of a negative cube root
!
		IF(deltaRHO .le. zero)THEN
			sign=-one
			deltaRHO=-deltaRHO
		ELSE
			sign=one
		END IF

		alphasens=sign*c1*K*(g*deltaRHO/(RHOair*nu*kappa))**(one/three)
!
!	adjust from flat plate value to	lake value
!	adjust sign of heat flux so that fluxes	from the water
!	to the air are negative
!
		alphasens=alphasens/two
		alphalat=alphasens/cp

		fluxsens=alphasens*(Tair-Tsurf)
		EVFLUX=alphalat*(q-qw)

		RETURN
		END SUBROUTINE STFLUX


!	**************************************************************************************
		SUBROUTINE LIUFLUX(ustar,qstar,tstar,Tair,Tsurf,Q,Qsat,uz,Theight,Qheight,z)
!****************************************************************************************
!	Calculates the roughness height z0 using an iterative approach
!
!	to evaluate surface fluxes, surface roughness, and stability of
!	the atmospheric surface layer from bulk parameters according to
!	liu el al. (79) jas 36 1722-1735
!	written by tim liu on 5/8/79
!
!	input:
!	uz wind speed in m/s
!	q humidity kg/m**3
!	z height of wind sensor
!	Qheight height of humidity sensor
!	po surface pressure in mb (DEFAULT to 1013.25)
!	id see SUBROUTINE drag for detail definition,id=1 (kondo),
!	 id=2 (smith), id=3 (large) (DEFAULT to 1)
!	usr,qsr scaling quantities for u,t,q
!	z0,zl roughness and stability parameters
!	rr,rt,rq roughness reynold numbers for u,t,q
	
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE		
		REAL*8 alphak	! vonKarman constant modified for T,Q
!IF		REAL*8 visc		! kinematic viscosity of air
		REAL*8 k		! vonKarman constant
		REAL*8 g		! gravitational constant
		REAL*8 pion2	! pi/2

		INTEGER*4 kondo,smith,large
		INTEGER*4 i		! loop index
		INTEGER*4 n		! iterations to converge on stability
	
		PARAMETER (kondo=1,smith=2,large=3,k=0.4,g=9.81)
 !IF     PARAMETER (visc=0.000015)
		PARAMETER (alphak=1.0/2.2)
		PARAMETER (pion2=1.57079633)

		REAL*8 cd		! drag coefficient adjusted for stability
		REAL*8 ce		! stability adjusted bulk evaporation coeff
		REAL*8 ch		! stability adjusted bulk temperature coeff
		REAL*8 dQ		! Qair - Qsat(Tsurf)
		REAL*8 dT		! = Tair - Tsurf
		REAL*8 L		! Obukhov length scale
		REAL*8 Lold		! Obukhov length from previous iteration
		
		REAL*8 psiT		! stability correction for temperature
		REAL*8 psiU		! stability correction for velocity
		REAL*8 psiQ		! stability correction for humdity

		REAL*8 Q		! observed specific humidity at height Qheight
		REAL*8 Qheight	! measurement height for humidity
		REAL*8 Qsat		! saturation specific humidity at Tsurf 
		REAL*8 Qstar	! specific humidity fluctuation scale

		REAL*8 Rr		! roughness Reynolds number
		REAL*8 Rt		! temperature scale Reynolds number
		REAL*8 Rq		! humidity scale Reynolds number

		REAL*8 Tair		! observed air temperature [C] at height Theight
		REAL*8 Tabs		! absolute observed air temperature [K]
		REAL*8 Theight	! measurement height for temperature
		REAL*8 Tstar	! temperature fluctuation scale
		REAL*8 Tsurf	! water surface temperature [C]
		REAL*8 Tvstar	! virtual temperature fluctuation scale

		REAL*8 u10		! wind speed at 10m height [m/s]
		REAL*8 ustar	! shear velocity
		REAL*8 uz		! wind speed at height z [m/s]
		REAL*8 X		! shape factor for velocity profile
		REAL*8 Y		! shape factor for temperature profile
		REAL*8 Yprime	! shape factor for humidity profile
		REAL*8 z		! height at which uz is measured
		REAL*8 z0		! initial value of roughness height
		REAL*8 zt		! temperature scale height
!ZQ		REAL*8 zq		! humidity scale height
		REAL*8 z0n		! new value of roughness height
		
		INTEGER*4 mode		
		CHARACTER*1 tab	
		INTEGER*4 ucheck		! int(10*u)
	

		tab=CHAR(9)
	
		z0 = .0005		! .0005 initial guess
		dT = Tair - Tsurf
		dQ = Q - Qsat
		Tabs = Tair + 273.15
!
!	Assume neutral stability for first time through
!
		psiT = 0.0
		psiQ = 0.0
		psiU = 0.0
		Lold = 0.0
		
		n = 0
	
5		CONTINUE

		i = 0
		n = n + 1
		IF (n .ge. 100) THEN
			WRITE(*,*)'stability did not converge w/in 100 iterations'
			STOP
		ENDIF
	
10		CONTINUE
		i = i + 1
		IF (i .ge. 100) THEN
			WRITE(*,*)'drag failed to converge within 100 iterations'
			STOP
		ENDIF
		
		ustar = k * uz / (LOG(z/z0) - psiU)

!	compute u10 for the current value of ustar
		u10 = 2.5* ustar * LOG(10.0 / z0)
		
!		compute drag coefficient for u10
!x		CALL drag(kondo, u10, cd)
!*****************
		IF (kondo-2) 3000,1000,2000
1000	cd = (0.61d0 + 0.063d0 * u10)/ 1000.
		GOTO 4000

2000	CONTINUE
		IF (u10 .lt. 11.0d0) THEN
			cd = 0.0012d0
		ELSE
			cd = (0.49d0 + 0.065d0 * u10) / 1000.0d0
		ENDIF
		GOTO 4000

3000	CONTINUE
		ucheck = 1000.0d0 * u10	
		SELECT CASE (ucheck)
			CASE (:300)
				cd = 0.0015d0			
			CASE (301:2200)
				cd = (1.08d0 * u10**(-0.15d0))/1000.0d0		
			CASE (2201:5000)
				cd = (0.771d0 + 0.0858d0 * u10)/1000.0d0		
			CASE (5001:8000)
				cd = (0.867d0 + 0.0667d0 * u10)/1000.0d0		
			CASE (8001:25000)
				cd = (1.2d0 + 0.025d0 * u10)/1000.0d0	
			CASE (25001:50000)
				cd = (0.073d0 * u10)/1000.0d0	
			CASE DEFAULT
				cd = 0.0037d0	
		END SELECT

4000	CONTINUE		

!*****************
!	compute z0 for new drag coefficient
		z0n = 10.0D0 / exp(k / SQRT(cd))
		
		IF (ABS(z0n-z0)/z0 .gt. 0.0001d0) THEN 
			z0 = z0n
			go to 10
		END IF
		
		z0 = z0n
!
!	Correct for atmospheric stability
!
		ustar = k * uz / (LOG(z/z0) - psiU)

!	Now get zq and zt using Liu et al's relation

		Rr = ustar * z0 / visc	
		mode = 1		! temperature
		CALL lkb (Rr,Rt,mode)
		zt = Rt * visc / ustar
		Tstar = alphak * dT / (LOG(Theight/zt) - psiT)
		
		mode = 2		! humidity
		CALL lkb (Rr,Rq,mode)
		zq = Rq * visc / ustar
		Qstar = alphak * dQ / (LOG(Qheight/zq) - psiQ)
		
		Tvstar = g * k * (Tstar * (1.0 + 0.61 * Q) + 0.61 * Tabs * Qstar)
		IF (Tvstar .ne. 0.0) THEN
			L = (Tabs * (1.0 + 0.61 * Q) * ustar * ustar)/Tvstar
		ELSE
	    		L = 0.0
		END IF

		IF (L .lt. 0.0) THEN
			X = (1.0 - 16.0 * z / L)**0.25
			psiU = 2.0 * LOG((X+1.0)/2.0)+ LOG((X * X + 1)/2.0) - 2.0 * atan(X) + pion2

			Y = SQRT(1.0 - 16.0 * Theight / L)
			psiT = 2.0 * LOG((1.0 + Y)/2.0)

			Yprime = SQRT(1.0 - 16.0 * Qheight / L)
			psiQ = 2.0 * LOG((1.0 + Yprime)/2.0)

		ELSEIF (L .eq. 0.0) THEN
	   	psiU = 0.0
			psiT = 0.0
			psiQ = 0.0

		ELSE
	   	psiU = -6.0 * LOG(1.0 + z/L)
			psiT = -6.0 * LOG(1.0 + Theight/L)
			psiQ = -6.0 * LOG(1.0 + Qheight/L)
		ENDIF
		
		IF (ABS((Lold - L)/(Lold + 1.0e-8)) .lt. 0.01) THEN
	
			IF(Tstar .eq. 0.0) THEN
				ch = 0.0
			ELSE
				ch = Tstar * ustar / (uz * dT)
			ENDIF
		
			IF(Qstar .eq. 0.0) THEN
				ce = 0.0
			ELSE
				ce = Qstar * ustar / (uz * dQ)
			ENDIF
				
			cd = ustar * ustar / (uz * uz)
				
			RETURN
		ENDIF
	
		Lold = L
		GOTO 5	
		END SUBROUTINE LIUFLUX


!	******************************************************************
	SUBROUTINE DRAG (id, u, cd)
!********************************************************************
!	To determine neutral drag coefficient cd from wind speed
!	at 10 m u in m/s
!	id=1 kondo(1975) blm 9 91-112
!	id=2 smith(1980) jpo 10 709-726
!	id=3 large & pond (1981) jpo 11 324-336
!	range of u specified are: kondo(.3,50),smith(6,22),l&p(4,25)
!	written by tim liu for vax on 2/10/82
!
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	
	
		INTEGER*4 id		! specify the relation
		INTEGER*4 ucheck		! int(10*u)
		
		REAL*8 cd			! neutral drag coefficient
		REAL*8 u			! wind speed at 10m

		IF (id-2) 3000,1000,2000
1000	cd = (0.61 + 0.063 * u)/ 1000.
		GOTO 4000

2000	CONTINUE
		IF (u .lt. 11.0) THEN
			cd = 0.0012
		ELSE
			cd = (.49 + .065 * u) / 1000.0
		END IF
		GOTO 4000

3000	CONTINUE
		ucheck = 1000.0 * u	
		SELECT CASE (ucheck)
			CASE (:300)
				cd = .0015			
			CASE (301:2200)
				cd = (1.08 * u**(-0.15))/1000.0		
			CASE (2201:5000)
				cd = (0.771 + 0.0858 * u)/1000.0		
			CASE (5001:8000)
				cd = (0.867 + 0.0667 * u)/1000.0		
			CASE (8001:25000)
				cd = (1.2 + 0.025 * u)/1000.0	
			CASE (25001:50000)
				cd = (0.073 * u)/1000.	
			CASE DEFAULT
				cd = 0.0037	
		END SELECT

4000	CONTINUE		
		RETURN
		END SUBROUTINE DRAG
		  
!	******************************************************************
	SUBROUTINE LKB (rr, rt, iflag)
!************************************************************************
!	to determine the lower boundary value rt of the logarithmic
!	 profiles of temperature (iflag=1) or humidity (iflag=2)
!	 in the atmosphere from roughness reynold number rr between
!	 0 and 1000. out of range rr indicated by rt=-999.
!	 based on liu et al. (1979) jas 36 1722-1723
!	 written by tim liu on 3/22/78, revised for vax on 2/10/82
!	 dyresm version by brad sherman 27/1/92
!
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	
		REAL*8 AA(8,2)
		REAL*8 BB(8,2)
		REAL*8 rank(8)
		REAL*8 rr		! input roughness Reynolds number z0*ustar/visc
		REAL*8 rt		! output temperature or humidity reynolds number

		INTEGER*4 i
		INTEGER*4 iflag		! 1 IF temperature, 2 IF humidity

		DATA AA /0.177,1.376,1.026,1.625,4.661,34.904,1667.19,5.88e5,		&
		    0.292,1.808,1.393,1.956,4.994,30.709,1448.68,2.98e5/

		DATA BB /0.,0.929,-0.599,-1.018,-1.475,-2.067,-2.907,-3.935,			&
		   0.,0.826,-0.528,-0.870,-1.297,-1.845,-2.682,-3.616/

		DATA rank /0.11,0.825,3.0,10.0,30.0,100.,300.,1000./

		i = 1
		IF (rr .le. 0.0 .OR. rr .gt. 1000.0) THEN
			rt = -999.0
			WRITE(*,*)'sub lkb Rr out of range ',rr
			pause
			RETURN
		END IF

10		CONTINUE
		IF (rr .le. rank(i)) GOTO 20
		i = i + 1
		GOTO 10

20		CONTINUE
		rt = AA(i, iflag) * rr ** BB(i, iflag)

		RETURN
		END SUBROUTINE LKB
