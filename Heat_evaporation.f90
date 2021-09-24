!**************************************************************************
		SUBROUTINE heat_evaporation(SVPW,HEATE,HEATCT,Ewat)
!***************************************************************************
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	
      
		INTEGER*4 i
!IF      REAL*8 densty       !density FUNCTION
		REAL*8 Ri_no		!Richardson number
		REAL*8 F_Ri 		!FUNCTION of the Richardson number
		REAL*8 WVCW			!water vapor concentration at the water surface kg vapor/kg air
		REAL*8 WVCA			!water vapor concentration of the surround air kg vapor/kg air
      REAL*8 SVPW			!saturated vapor pressure mb at the water surface temperature
      ! REAL*8 SVPD                   !SVPD saturated vap. pres mb at the dew point temperuate (oc)
      REAL*8 TVAvapor_temp !FUNCTION
		REAL*8 TVAVAPOR_TEMP_RH !FUNCTION
		REAL*8 AP			!atmospheic pressure at height h meter above the mean sea level
		REAL*8 AP_P0		!AP/P0 ratio
      REAL*8 MIXRD		!mix ratio at dew point temperature
		REAL*8 MIXRW		!mix ratio at water surface temperature
		REAL*8 RHOAIR       !density of air at dew point temperature
		REAL*8 RHOAIRW		!desnity of air at water surface temperautre 
		REAL*8 RHOWAT(MAXNS) !density of water                   
		REAL*8 Ewat			!Evaporated water depth m/s
		REAL*8 HEATE		!Heat due to evaporation
      REAL*8 alpha_E      !value for no-wind condition variable, evaporative heat flux
      REAL*8 HEATS        !Sensible heat flux
		REAL*8 alpha_S      !value for no-wind condition variable, Sensible heat flux
      REAL*8 HEATC        !Condensation heat
      REAL*8 HEAT_AE      !Heat advection by evaporation
      REAL*8 HEAT_AC      !Heat advection by condensation
		REAL*8 HEAT_AR      !Heat advection by rain
      REAL*8 HEATCT		!Total heat by conduction and sensible heat.
		REAL*8 LATHEATV		!Latent heat of vaporization
!IF      REAL*8 SATVAP       !DLM
		REAL(8), PARAMETER:: P0=1013.0d0	!Standard atmospheric pressure at mean sea level, mb
		REAL(8), PARAMETER:: EC_N=2.26d-6	!1.89, 2.26, 2.48 Evaporation coefficient for lake data
		REAL(8), PARAMETER:: HT_AN=10.0		!wind speed anemometer height from the lake surface
		REAL(8), PARAMETER:: gamma=1.0d0	!the constant term for the FUNCTION of Richardson number
		REAL(8), PARAMETER:: CCL=0.5d0	    !correction coeficient for lake surface  ! you can change
		REAL(8), PARAMETER:: MHCA=100.0/3600.0 !Molecular heat conductivity of air
		REAL(8), PARAMETER:: SHCA=1010.0       !Specific heat capacity of air
		REAL(8), PARAMETER:: SHCW=4180.0       !Specific heat capacity of water 
		REAL(8), PARAMETER:: MDA=0.077/3600.0  !Molecular diffusivity of air m2/s
		REAL(8), PARAMETER:: KVA=0.0548/3600.0 !Kinematic viscosity of air m2/s
        REAL(8), PARAMETER:: CB=1.0d0         !Coefficient see Martin & McCutchen
!--------------------------------------------------------------------------------
        REAL*8 F_Ri_E 		!FUNCTION of the Richardson number for Evaporation
		REAL*8 F_Ri_S 		!FUNCTION of the Richardson number for sensible heat
        REAL*8 F_Ri_CD 		!FUNCTION of the Richardson number for Wind Drag
		REAL*8 CDA,CDB      !Coefficients of a and b for wind drag, and F_Ri_E and Fi_Ri_S
		REAL*8 CSN, CEN		!Neutral stability coefficient for sensible heat and evaporation
		REAL*8 CDN          !Neutral stability coefficient for drag
		REAL*8 WINDC		!Corrected wind Strub and Powell (1987)

       PARAMETER (CDN=1.3d-3,CEN=1.4d-3, CSN=1.4d-3)  !Modern equation coefficients
      ! PARAMETER (CDN=1.3d-3,CEN=2.26d-6, CSN=2.26d-6)   !TVA equation lower coefficients
       !PARAMETER (CDN=1.3d-3,CEN=2.48d-6, CSN=2.48d-6)   !TVA equation upper coefficients
	   !PARAMETER (CDN=1.3d-3,CEN=1.0d-3, CSN=1.0d-3)
      ! PARAMETER (CDN=1.3d-3,CEN=12.00d-6, CSN=0.90*36.0d-6) 
      ! (CDN=1.3d-3, CEN=12.5d-6, CSN=31.00d-6) 
      !(CDN=1.3d-3, CEN=13.00d-6, CSN=60.000d-6)
      !(CDN=1.3d-3, CEN=12.000d-6, CSN=49.000d-6) 
 !     PARAMETER (CDN=1.3d-3, CEN=1.82d-6, CSN=3.00d-6)  !WB Paper for 4900.00 GW
               !(CDN=1.3d-3, CEN=1.82d-6, CSN=3.00d-6)                
               !(CDN=1.3d-3, CEN=2.05d-6, CSN=3.50d-6)
               !(CDN=1.3d-3, CEN=1.85d-6, CSN=1.85d-6)
!--------------------------------------------------------------------------------
		DO i = 1,ns
			den(i) = densty(temp(i),sal(i))
			rhowat(i)=1000.0d0+den(i)				
		ENDDO
!-----Initialization----------------------------------------------------------------
      EWat=0.0d0
		IF(u6x.le.0.1) u6x=0.1d0
		WINDC=u6x-0.05*u6x
!--------(1) Evaporation and evaporative heat loss (TVA 1972)----------------
		!AP_P0=((288.0-0.0065*(BASE+CRL))/288.0)**5.256
		!AP=AP_P0*(P0)
		AP = (P0)*(1- (0.0065*(BASE+CRL))/((0.0065*(BASE+CRL))+T4+273.15))**5.256
			
		!AP=AP_P0*(P0)/(20.0+273.15)*(T4+273.15)
		IF(HUMIDITY.EQ.1) THEN
			!SVPD=TVAvapor_temp_RH(T4,rh)	 Replaced by power equation below SCT 4/22/2019
			SVPD = (rh)*(10.0**(9.286-(2322.38/(T4+273.15))))
		ENDIF
		!SVPW=TVAvapor_temp(temp(ns))! Saturation value correct as such sct
		SVPW=10.0**(9.286-(2322.38/(temp(ns)+273.15)))
!	SVPW=satvap(temp(ns))
		WVCW=0.622*SVPW/AP
      WVCA=0.622*SVPD/AP
      MIXRD=0.622*SVPD/(AP-SVPD)	
      MIXRW=0.622*SVPW/(AP-SVPW)
	  !RHOAIR=0.348*(AP/(T4+273.15))*((1.0+MIXRD)/(1.0+1.61*MIXRD)) SCT replaced TVA Equation that was giving wierd results, consider alternatives
	  RHOAIR = (AP-SVPD)*100.0/287.058/(T4+273.15) + SVPD*100.0/461.495/(T4+273.15)
	  RHOAIRW=0.348*(AP/(temp(ns)+273.15))*((1.0+MIXRW)/(1.0+1.61*MIXRW))
      LATHEATV=(2500.0-2.39*temp(ns))*1000.0
!		LATHEATV = (- 0.0000614342*temp(ns)**3 + 0.00158927*temp(ns)**2   &    
!						- 2.36418*temp(ns) + 2500.79)*1000.0 



     ! write(*,fmt='(i7,7f10.2)')jday,SVPD,SVPW,WVCA,WVCW,RHOAIR,RHOAIRW,AP
!		LATHEATV = 2.453d+9/RHOWAT(NS)   
      IF(u6x.GT.0.0) THEN
!			IF(EC_N.eq.1.89d-6) THEN
!				Ri_no=0.0d0
!			ELSE
			Ri_no=-9.80*(RHOAIR-RHOAIRW)*HT_AN/(RHOAIR*(WINDC**2.0))	
!			ENDIF
!--------set the maximum and minimum limit of Richardson number see TVA, 1972 pp5.15
			IF (Ri_no.le.-1.0d0) THEN
				Ri_no=-1.0d0
			ELSEIF(Ri_no.ge.2.0) THEN
				Ri_no=2.0d0
			ENDIF
			IF(Ri_no.GT.0.0) THEN
				F_Ri_E=1.0
				F_Ri_S=1.0
			ELSEIF(Ri_no.LE.0.0) THEN
				CDA=0.83*(CDN**(-0.62))
				CDB=0.25*(CDN**(-0.80))
				F_Ri_E=1.0
				F_Ri_S=1.0
			ENDIF
			!EWat=-CEN*AP_P0*u6x*(WVCW-WVCA)      !m/s)	  !TVA Eq ST
			EWat=-CEN*(RHOAIR/RHOWAT(NS))*u6x*(WVCW-WVCA)   ! Modern Eq ST
!	write(*,fmt='(4f10.2)') F_Ri_E,Ri_no, T4, temp(ns)
		! ELSEIF(u6x.EQ.0.0) THEN	! sct clearing unlikely conditions
			! IF(RHOAIR.le.RHOAIRW) THEN
				! alpha_E=0.0d0
			! ELSE	  
				! alpha_E=0.137*CCL*(MHCA/(RHOWAT(NS)*SHCA))*					&
						! (9.81*(RHOAIR-RHOAIRW)/(RHOAIR*KVA*MDA))**(1.0/3.0) 
			! ENDIF        
			! EWat=-alpha_E*(WVCW-WVCA)					!m/s
		ENDIF
!	PRINT*,Ri_no,F_Ri
      HEATE=-RHOWAT(NS)*LATHEATV*(-EWat) !to make both negative   
		IF(HEATE.GT.0.0) HEATE = 0.0	
		IF(EWat.GT.0.0) EWat = 0.0	
!--------(2) Sensible heat loss TVA, 1972----------------------------------------
			  !HEATS=-CB*RHOWAT(NS)*SHCA*CSN*AP_P0*u6x*(TEMP(NS)-T4)      !TVA Eq ST
			HEATS=-CB*RHOWAT(NS)*SHCA*CSN*(RHOAIR/RHOWAT(NS))*u6x*(TEMP(NS)-T4)    ! Modern Eq ST
			!write(*,fmt='(i5,4f10.2)')jday,HEATE,HEATS, T4,temp(ns) 	  
		! ELSEIF(u6x.EQ.0.0) THEN	! sct clearing unlikely conditions
			! IF((RHOAIR-RHOAIRW).gt.0.0) THEN
				! alpha_S=0.137*CCL*MHCA*(9.81*(RHOAIR-RHOAIRW)/(RHOAIR*KVA*MDA))**(1.0/3.0)         
				! HEATS=-alpha_S*(TEMP(NS)-T4)					!
			! ELSEIF((RHOAIR-RHOAIRW).le.0.0) THEN
				! HEATS=-0.0d0
			! ENDIF 		
! !--------(3) Condensation TVA, 1972----------------------------------------	! sct clearing unlikely conditions
! !     Condensation flux under still air is assumed to be zero.
		! HEATC=RHOWAT(NS)*LATHEATV*EC_N*AP_P0*u6x*(WVCA-WVCW)*F_Ri_E
! !--------(4) Advection heat fluxes, TVA 1972---------------------------------
! !--------(4(a)) Advection heat by evaporation, TVA 1972---------------------------------
		! HEAT_AE=-RHOWAT(NS)*SHCW*(-EWat)*TEMP(NS)
! !--------(4(b)) Advection heat by condensation, TVA 1972---------------------------------
		! HEAT_AC=-RHOWAT(NS)*SHCW*EC_N*AP_P0*u6x*(WVCW-WVCA)*T4
! !--------(4(c)) Advection heat by rainfall, TVA 1972---------------------------------
      ! HEAT_AR=RHOWAT(NS)*SHCW*(RAIN/1000.0d0)*T4
!		either HEATE or HEATC can occur at the same time.      
!		HEATCT=HEATS +HEAT_AE+HEAT_AC
!		HEATCT=HEATS +HEAT_AC
      HEATCT=HEATS
!		PRINT*, HEAT_AC, T4
!		PRINT*,HEATE, HEATS,HEAT_AE+HEAT_AC
         ! print*, jday,HEATE,HEATS,temp(ns), T4, RH, SVPD, SVPW, WVCA, WVCA, AP, rhowat(ns),RHOAIR,u6x
       ! write(*,fmt='(i7,4f10.2)')jday,HEATE,HEATS, T4,temp(ns) 	  
	  
 		RETURN
		END SUBROUTINE heat_evaporation
!************service functions and subroutines---------------------------
!------------using only temperature dew point or water or air temp ------
		REAL*8 FUNCTION TVAvapor_temp(temp)
!************************************************************************
		IMPLICIT NONE

		REAL*8 temp
		REAL(8), PARAMETER:: a0=9.5, b0=265.5, c0=0.7858
		REAL(8), PARAMETER:: a1=7.5, b1=237.3, c1=0.7858
		IF(temp.le.0.0) THEN
			TVAvapor_temp = dexp(2.3026*(c0+a0*temp/(b0+temp)))
		ELSEIF(temp.gt.0.0) THEN
			 TVAvapor_temp = dexp(2.3026*(c1+a1*temp/(b1+temp)))
		ENDIF
		RETURN
		END FUNCTION TVAvapor_temp
!****************using only air temperature and relative humidity*********
		REAL*8 FUNCTION TVAvapor_temp_RH(airtemp,rh)
!************************************************************************
		IMPLICIT NONE

		REAL*8 airtemp,RH
		REAL(8), PARAMETER:: a0=9.5, b0=265.5, c0=0.7858
		REAL(8), PARAMETER:: a1=7.5, b1=237.3, c1=0.7858
		TVAvapor_temp_RH = (rh)*(10.0**(9.286-(2322.38/(airtemp+273.15))))
		!IF(airtemp.le.0.0) THEN
		!	TVAvapor_temp_RH = (rh)*dexp(2.3026*(c0+a0*airtemp/(b0+airtemp)))
		!ELSEIF(airtemp.gt.0.0) THEN
		!	TVAvapor_temp_RH = (rh)*dexp(2.3026*(c1+a1*airtemp/(b1+airtemp)))
		!ENDIF
		RETURN
		END FUNCTION TVAvapor_temp_RH
!******************************************************