!********************************************************************
	SUBROUTINE LAKENUM(ULKN,lk)
!*******************************************************************
!
!        lakenumber
!       ---------------
!
!MC     THIS ROUTINE CALCULATES THE lakenumber FOR THE CURRENT DAY
!MC     AND SENDS THE OUTPUT TO A FILE:'Ldddhhmm.NUM'
!KLG  i.e. has the same name as the simulation file
!MC
!MC     ALGORITHM FROM PROGRAM "LAKE NUMBER PROJECT" (11/8/1991) 
!MC     
!MC     WRITTEN BY:     STEPHEN BARRY AND JOHN SNELL
!MC
!MC     FOR:            CENTRE FOR WATER RESEARCH, UNIVERSITY OF W.A.
!MC
!MC     TRANSLATED FROM C LANGUAGE BY CAMERON CALZONI 11/2/1992
!KLG  Put into Vers 6.75.2 DW by KLG July 1992
!MC
!MC     NOTES:  nlayers = NUMBER OF LAYERS OF INTERP DATA(STEP = 0.1)
!MC             ALL ARRAYS BEGIN WITH SUBSCRIPT OF ZERO
!MC             all arrays of data interp to 0.1m have subscript of zero
!MC             referring to depth of 0.1m (Top of first layer)
!MC
!MC                     CHECKED AND O.K.
!     Modified for Lake Tahoe
!     Includes variables for thermal sensitivity analysis
!     Written by Joaquim P. Losada and Geoff Schladow
!     UC-Davis, January 2001
!*******************************************************************	implicit none
		USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE			 

		INTEGER(4)     :: ulkn,uwor,nlayers
		INTEGER(4)     :: i,j,ii,jj,lk,njday,ik
		INTEGER(4)     :: uper_i,lower_i,t_min(1), t_max(1)
		LOGICAL*4 NOPROBLEMS

		REAL(8)	:: wind_speed,print_hsave 
		REAL(8)	:: lake_top,lake_bottom
		REAL(8)	:: uper,lower,XPE
		REAL(8)	:: usquared,thermdepth,lakenumber
		REAL(8)	:: xp,zc,xmass, IsoTemp, VolSum, IsoSal, Mass
!		REAL(8) :: IsoMass, IsoDen, IsoVol, VolExtra, ijk
!
!     Definitions for Statistics Analysis
!
      REAL(8)	:: ave_low  ,adev_low  ,sdev_low,var_low,skew_low,curt_low
		REAL(8)	:: ave_up   ,adev_up   ,sdev_up ,var_up ,skew_up ,curt_up
		REAL(8)	:: ave_low_D,adev_low_D,sdev_low_D
		REAL(8)	:: var_low_D,skew_low_D,curt_low_D
		REAL(8)	:: ave_up_D ,adev_up_D ,sdev_up_D 
		REAL(8)	:: var_up_D, skew_up_D, curt_up_D

!IF		REAL(8)	:: calc_XMOM
!		REAL(8)	:: Mix_depth_F,Mix_depth_G
		REAL(8), PARAMETER	:: grav = 9.81d0
!		REAL(8) :: AreaTahoe

!		REAL(8), PARAMETER	:: Amax		= 465969.4*1000.0d0
!		REAL(8), PARAMETER	:: DepthMax = 501.45d0
!		REAL(8), PARAMETER	:: Shape	= 0.51d0
!
!     Interpolate density and temperature each 0.1m
!		
		CALL INTERPOLATE_LAYER_DATA(nlayers,idepth,density,tempy,salty)

		IF (nlayers.lt.1) THEN
			NOPROBLEMS = .FALSE.
			GOTO 4406
		ELSE
			NOPROBLEMS = .true.
		ENDIF

		wind_speed	= u6x
		IF (wind_speed .lt. 0.05) wind_speed = 0.05

      usquared = 1.612e-6*wind_speed*wind_speed

		xmom = CALC_XMOM(nlayers)

!
!     Store for use in Optima
!
		IF(xmom.lt.0.0d0) THEN
			thermdepth = -99.0
			vector_thermdepth(lk) = -99.0
			vector_thermtemp (lk) = -99.0
			lower_i    = -99.0
			uper_i     = -99.0
			n_low      = -99.0
			adev_low   = -99.0
			adev_up    = -99.0
			adev_low_D = -99.0
			adev_up_D  = -99.0
			lakenumber = -99.0
			n_up       = -99.0
			n_up_D     = -99.0
			n_low_D    = -99.0
			GOTO 4406
		ELSE
			thermdepth = idepth(nlayers-1)-xmom
		ENDIF

		vector_thermdepth(lk) = xmom
		vector_day(lk)        = jday
!
!     Temperature at thermdepth.Used for sensitivity analysis
!
      lower_i    = 0
		uper_i     = 0
      n_low      = 0.0d0
		adev_low   = 0.0d0
		adev_up    = 0.0d0
		adev_low_D = 0.0d0
		adev_up_D  = 0.0d0
      n_up       = 0.0d0
      n_up_D     = 0.0d0
      n_low_D    = 0.0d0

!
!     Case when Lake Number is zero
!
		DO i = 0,nlayers-1
			IF((vector_thermdepth(lk)-idepth(i)).lt.0.1.and.		&
		  		vector_thermdepth(lk).gt.idepth(i))THEN
				lower = tempy(i)
				lower_i = i
			ENDIF
		
			IF((idepth(i)-vector_thermdepth(lk)).lt.0.1.and.		&
		  		vector_thermdepth(lk).lt.idepth(i))THEN
				uper = tempy(i)
				uper_i = i
			ENDIF
		ENDDO
	
		IF(abs(uper-lower).lt.0.1) THEN
			vector_thermtemp(lk) = lower
		ELSE 
			IF ((idepth(uper_i)-idepth(lower_i)).le.0.0d0) THEN
				vector_thermtemp(lk) = lower
			ELSE    
				vector_thermtemp(lk) = lower +(uper-lower)/(idepth(uper_i)-idepth(lower_i))*		&
     					(vector_thermdepth(lk)-idepth(lower_i))
			ENDIF
		ENDIF

		DO i = 0,lower_i 
!		DO i = 1,lower_i + 1
			adev_low   = adev_low   + (vector_thermtemp (lk) -tempy(i))**2
			adev_low_D = adev_low_D + (vector_thermdepth(lk)-idepth(i))**2
!			n_low(i)   = tempy(i-1)
!			n_low_D(i) = idepth(i-1)
		ENDDO

		DO i = uper_i,nlayers-1
!		DO i = uper_i,nlayers-1
			adev_up   = adev_up   + (vector_thermtemp (lk) -tempy(i))**2
			adev_up_D = adev_up_D + (vector_thermdepth(lk)-idepth(i))**2
!			n_up(i)   = tempy(i-1)
!			n_up_D(i) = idepth(i-1)
		ENDDO 
!
!     Calculate the deviation. # of cases in denominator already
!     account #-1, by adding 1 extra point corresponding to the T mean     
!
      adev_low   = sqrt(adev_low   /lower_i) 
		adev_low_D = sqrt(adev_low_D /lower_i) 

      adev_up   = sqrt(adev_up  /(nlayers-uper_i)) 
		adev_up_D = sqrt(adev_up_D/(nlayers-uper_i)) 

!
!     Temperature lower deviation around TermDepth
!
!      call moment(n_low,lower_i+1,ave_low,adev_low,
!     &	        sdev_low,var_low,skew_low,curt_low)
!
!     Temperature upper deviation around TermDepth
!
!    	call moment(n_up,uper_i+1,ave_up,adev_up,
!     &	        sdev_up,var_up,skew_up,curt_up)		
!
!     Thermocline Depth lower deviation around TermDepth
!
!      call moment(n_low_D,lower_i+1,ave_low_D,adev_low_D,
!     &	        sdev_low_D,var_low_D,skew_low_d,curt_low_D)
!
!     Thermocline Depth upper deviation around TermDepth
!
!    	call moment(n_up_U,uper_i+1,ave_up_D,adev_up_D,
!     &	        sdev_up_D,var_up_D,skew_up_D,curt_up_D)		

8182	CONTINUE

		CALL LNPE3(nlayers,xp,zc,xmass,xpe)
	 	
		lake_top   = -xp*(xmass*1e6)*(1.0-((idepth(nlayers-1)-thermdepth)/idepth(nlayers-1)))
		lake_bottom= 1000*usquared*((A((nlayers-1)+1)*1e6)**1.5)*(1.0-(zc/idepth(nlayers-1)))
	
		IF (lake_top.eq.0.0d0.and.lake_bottom.eq.0.0d0) THEN
			noproblems = .false.
			GOTO 4406
		ENDIF

		IF (lake_bottom.eq.0.0d0) THEN
			noproblems = .false.
			GOTO 4406
		ENDIF

		lakenumber=lake_top/lake_bottom
	
		IF(lakenumber.lt.0.0d0) lakenumber	= -99.0

4406  CONTINUE

		IF (NOPROBLEMS) THEN		 	
			WRITE(ulkn,4407) jday,lakenumber,depmx,vector_thermdepth(lk),					&
				vector_thermdepth(lk) - adev_low_D,vector_thermdepth(lk) + adev_up_D,   &   
				vector_thermtemp (lk), vector_thermtemp (lk) - adev_low,						&
				vector_thermtemp (lk) + adev_up


4407		FORMAT(I8,f8.1,x,7f10.2)

		ELSE
			WRITE(ULKN,4408) jday
4408		FORMAT(I10,'  NaN')
		ENDIF

		RETURN
		END SUBROUTINE LAKENUM

