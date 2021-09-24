!***************************************************************
		REAL*8 FUNCTION RHtoSVPD(Tair, Pair)
!****************************************************************
!	This function converts relative humidity and air temperature to
!	vapor pressure using the relations in Gill's Atmosphere-Ocean dynamics
!
!	use this to compute SVPD   ! ambient vapor pressure of atmosphere {mb}
!
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE
!
!	constants
!
		REAL*8 epsilon	! ratio of gas constants (dry air)/(water vapor)
		REAL*8 Ra		! gas constant for dry air {J/kg-K}
		REAL*8 Rv		! gas constant for water vapor {J/kg-K}

		PARAMETER (epsilon=0.62197,Ra=287.04,Rv=461.50)

		REAL*8 one,hundrd,dCtodK

		PARAMETER (one=1.0d0,hundrd=100.0d0,dCtodK=273.15d0)

		REAL*8 Pair		! atmospheric pressure {mb}
		REAL*8 q		! observed specific humidity
		REAL*8 qw		! specific humidity for saturated air
		REAL*8 r		! observed mixing ratio
		REAL*8 rw		! saturated mixing ratio for moist air
!IF		REAL*8 RH		! relative humidity (%)
!IF		REAL*8 RelHum		! relative humidity {decimal}
		REAL*8 RHOair		! ambient air density {kg/m^3}
!		REAL*8 RHOsurf	! density of saturated air @ surface {kg/m^3}
!		REAL*8 sign		! sign of cube root
!IF		REAL*8 SATVAP		! function to calculate saturated vapor pres.
		REAL*8 SVPair		! saturated vapor pressure at Tair {mb}
!		REAL*8 SVPD		! ambient vapor pressure of atmosphere {mb}
		REAL*8 Tair		! air temperature (C)

!
!	calculate specific humidities and densities from observed data
!
		RelHum=RH
		SVPair=SATVAP(Tair)
		rw=SVPair/(Pair-SVPair)*epsilon
		r=RelHum*rw
		qw=rw/(one+rw)
		q=r/(one+r)

		RHOair=Pair*hundrd/(Ra*(Tair+dCtodK)*(one+(-one+one/epsilon)*q))
!
!	calculate atmospheric vapor pressure for evaporative fluxes
!	into moving air, i.e. SVPD is used in bulk formula for vapor flux.
!	
		RHtoSVPD = RHOair*q*Rv*(Tair+dCtodK)/hundrd

		RETURN
		END FUNCTION RHtoSVPD 



