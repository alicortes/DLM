!***********************************************************************************************
      SUBROUTINE ATMOS_DEPOSITION(jyear,redn_Part,redn_TP,redn_TN,N_Atm,P_Atm,Part_Atm)
!***********************************************************************************************
!           -------- Atmospheric Deposition --------
!    Nitrogen forms:   NH4, NO3, DON
!    Phosphorus forms: DOP, DRP, SRP, POP, THP (=BAP)
!    Constant rates are based on LTADS (2005) report and TERC DATA
!************************************************************************
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

		INTEGER*4 i,j,k,jyear
! 
!		Nutrient Budget Control
! 
		REAL*8 N_Atm,P_Atm, Part_Atm(7)
! 
!		Atmospheric Neutrients deposition rates for Lake Tahoe (mg/m2/day)
!		Units should be in mg/m2/day so that it will match with mg/m3 or ug/L
  		REAL*8 surfacelayer_depth
!---------Nitrogen Neutrients-----------------------
		REAL*8 rate_NO3, rate_NH4, rate_PON, rate_DON
		REAL*8 rate0_NO3, rate0_NH4, rate0_PON, rate0_DON
		REAL*8 rate1_NO3, rate1_NH4, rate1_PON, rate1_DON
!---------Phophorous Neutrients--------------------
      REAL*8	rate_SRP,rate_THP, rate_POP, rate_DOP
		REAL*8 	rate0_SRP,rate0_POP, rate0_DOP
		REAL*8 	rate1_SRP,rate1_POP, rate1_DOP
!---------Details in season-wise--------------------
!************START DRY Nitrogen**********************************************
!------------Winter Nitrogen Dry---------------------------------------------
      REAL*8 winter_NO3_DRY,winter_NH4_DRY,winter_DON_DRY,winter_PN_DRY     	 

      PARAMETER (winter_NO3_DRY = 0.17884, winter_NH4_DRY = 0.53652,		&
     				  winter_DON_DRY = 0.40060, winter_PN_DRY  = 0.04650)    
!-------------Spring Nitrogen Dry -------------------------------------------- 
		REAL*8 spring_NO3_DRY,spring_NH4_DRY,spring_DON_DRY,spring_PN_DRY

      PARAMETER (spring_NO3_DRY = 0.13757, spring_NH4_DRY = 0.41270,		&
     	           spring_DON_DRY = 0.22891, spring_PN_DRY  = 0.03247)

!-------------Summer Nitrogen Dry -------------------------------------------- 
		REAL*8 summer_NO3_DRY, summer_NH4_DRY,summer_DON_DRY,summer_PN_DRY      	

      PARAMETER (summer_NO3_DRY = 0.22723, summer_NH4_DRY = 0.53020,		&
     	           summer_DON_DRY = 0.15593, summer_PN_DRY  = 0.03595)							   

!-------------Fall Nitrogen Dry -------------------------------------------- 
		REAL*8 fall_NO3_DRY,fall_NH4_DRY,fall_DON_DRY,fall_PN_DRY
     

      PARAMETER (fall_NO3_DRY = 0.23201, fall_NH4_DRY = 0.66702,			&
     	           fall_DON_DRY = 0.12041, fall_PN_DRY  = 0.04489)
!********************END DRY Nitrogen***************************************
!********************START DRY Phosphrous***********************************
!--------------Winter (December to June) Phophrous Dry --------------------
		REAL*8	winter_SRP_DRY,winter_POP_DRY,winter_DOP_DRY
	
		PARAMETER (winter_SRP_DRY = 0.00593,winter_POP_DRY = 0.01423,winter_DOP_DRY = 0.00474)

!-------------Summer (July to November) Phosphorous Dry-------------------- 
		REAL*8 summer_SRP_DRY,summer_POP_DRY,summer_DOP_DRY
	
		PARAMETER (summer_SRP_DRY = 0.01341,summer_POP_DRY = 0.02850,summer_DOP_DRY = 0.01174) 
!*******************END DRY Phophrous*************************************
!*******************START of WET Nitrogen and Phosphrous*****************	
!--------------Wet deposition same for all season---------------------
		REAL*8 NO3_WET,NH4_WET,DON_WET,PN_WET,SRP_WET,POP_WET,DOP_WET

		PARAMETER (NO3_WET = 0.68980618, NH4_WET = 0.65148361,				&
     				  DON_WET = 0.82930032, PN_WET = 0.09044125)

		PARAMETER (SRP_WET = 0.03832257,POP_WET = 0.03832257,DOP_WET = 0.03065805)
!****************END of WET Nitrogen and Phosphrous**********************
!-----------------------------------------------------------------------
!
!    Atmospheric Paricles deposition rates for Lake Tahoe (# particles/day/m2) 
		REAL(8):: rate_Particles(7)
		REAL(8):: winter_dry_Particles_rate(7), winter_wet_Particles_rate(7),	&
					 spring_dry_Particles_rate(7), spring_wet_Particles_rate(7),	&
					 summer_dry_Particles_rate(7), summer_wet_Particles_rate(7),	&
					   fall_dry_Particles_rate(7),   fall_wet_Particles_rate(7)

		DATA winter_dry_Particles_rate /9.072598E+07,2.996198E+07,					&
			 3.932766E+06,3.894573E+06,6.112662E+05,1.180736E+05,0.00E+00/

		DATA winter_wet_Particles_rate /8.116250E+08,2.680367E+08,					&
			 8.154126E+06,8.072247E+06,2.389807E+05,4.776825E+04,0.00E+00/

		DATA spring_dry_Particles_rate /5.418927E+07,1.789584E+07,					&
			 3.467331E+06,3.429486E+06,5.013876E+05,9.704511E+04,0.00E+00/

		DATA spring_wet_Particles_rate /8.113205E+08,2.679361E+08,					&
			 1.936294E+07,1.918008E+07,2.583086E+05,5.241607E+04,0.00E+00/


		DATA summer_dry_Particles_rate /6.780904E+07,2.239373E+07,					&
			 3.037390E+06,2.990081E+06,6.304551E+05,1.246261E+05,0.00E+00/

		DATA summer_wet_Particles_rate /3.868652E+08,1.277610E+08,					&
			 5.181073E+06,4.890468E+06,2.202288E+05,5.227175E+04,0.00E+00/


		DATA fall_dry_Particles_rate /8.827393E+07,2.915219E+07,						&
			 3.913478E+06,3.702312E+06,4.622116E+05,8.972125E+04,0.00E+00/

		DATA fall_wet_Particles_rate /3.167317E+07,1.045997E+07,						&
			 5.675661E+05,5.514941E+05,2.555857E+04,5.609794E+03,0.00E+00/
!------------------------------------------------------------------------- 
		REAL*8 redn, redn_part, redn_TP, redn_TN

!-----------------------------------------------------------------
!-------Initialize to zero to avoid cummulating the values in every time step-----
		rate_NO3 = 0.0				
		rate_NH4 = 0.0
		rate_DON = 0.0
		rate_PON = 0.0 
		rate_SRP = 0.0
		rate_POP = 0.0
		rate_DOP = 0.0
		rate_THP = 0.0

		rate0_NO3 = 0.0				
		rate0_NH4 = 0.0
		rate0_DON = 0.0
		rate0_PON = 0.0 
		rate0_SRP = 0.0
		rate0_POP = 0.0
		rate0_DOP = 0.0

		rate1_NO3 = 0.0				
		rate1_NH4 = 0.0
		rate1_DON = 0.0
		rate1_PON = 0.0 
		rate1_SRP = 0.0
		rate1_POP = 0.0
		rate1_DOP = 0.0	   

		DO i = 1,7
			rate_Particles(i)=0.0
		ENDDO
!		WRITE(*,10)jday,jyear, redn_part,redn_TP,redn_TN
!10	FORMAT(2i8, 3f10.2, 4x, 'Atmospher')
!		WRITE(1000,*) jyear,jday,redn,'Atmos'
!-------End of initialization----------------------------------------------------
!		redn_part1=1.0d0-0.00d0/100.0d0
!		redn_TP1  =1.0d0-0.00d0/100.0d0
!		redn_TN1  =1.0d0-0.00d0/100.0d0
		WRITE(*,10)jday,jyear, redn_Part,redn_TP,redn_TN
10		FORMAT(2i8, 3f10.2, 4x, 'Atmospher')

!    Depth (m) of the surface layer
! 
      surfacelayer_depth = depth(ns)-depth(ns-1)
!***************************************************************************	
!***********Inert Atmospheric Particles load into Tahoe Lake****************
!***************************************************************************
!           Winter: (Jan, Feb, March)  90 days   1-90
!           Spring: (April, May, June) 91 days  91-181
!           Summer: (July, Aug, Sep)   92 days 182-273
!           Fall  : (Oct, Nov, Dec)    92 days 274-365
!----------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
		DO i =1,7
			IF(jday.ge.1.and.jday.le.90) THEN  !Winter
				IF(Rain.ge. 2.54)THEN	 !Wet days for rainfall>=0.1 inch or 2.54 mm
					rate_Particles(i)=1.95*winter_wet_Particles_rate(i)
				ELSE	       
					rate_Particles(i)=1.0*winter_dry_Particles_rate(i)  
				ENDIF
			ELSEIF(jday.ge.91.and.jday.le.181) THEN	  !Spring
	  			IF(Rain.ge. 2.54)THEN
	    			rate_Particles(i)=1.95*spring_wet_Particles_rate(i)   
				ELSE
					rate_Particles(i)=1.0*spring_dry_Particles_rate(i)  
				ENDIF
			ELSEIF(jday.ge.182.and.jday.le.273) THEN   !Summer
	  			IF(Rain.ge. 2.54)THEN
	       		rate_Particles(i)=1.95*summer_wet_Particles_rate(i) 
				ELSE       
					rate_Particles(i)=1.0*summer_dry_Particles_rate(i)
				ENDIF  
			ELSE	 !Fall	 
				IF(Rain.ge. 2.54)THEN
					rate_Particles(i)=1.95*fall_wet_Particles_rate(i) 
				ELSE
					rate_Particles(i)=1.0*fall_dry_Particles_rate(i)	
				ENDIF
			ENDIF
		ENDDO

		DO i = 1,7
			cf(ns,i) = cf(ns,i) + redn_Part*rate_Particles(i)/(surfacelayer_depth)		   !redn_Part 
			Part_Atm(i) = Part_Atm(i) + redn_Part*													&
							 (rate_Particles(i)/(surfacelayer_depth))*vol(ns)*1000.0d0			!redn_Part
	ENDDO
	
!********************************************************************************
!***********Water Soluble Atmospheric Nutrients load into Tahoe Lake*************
!********************************************************************************
		IF(jday.ge.1.and.jday.le.90) THEN  !Winter  (January to March)
			IF(Rain.ge. 2.54)THEN	 !Wet days for rainfall>=0.1 inch or 2.54 mm
				rate0_NO3 =	NO3_WET			
				rate0_NH4 = NH4_WET
				rate0_DON =	DON_WET
				rate0_PON =	PN_WET
				rate0_SRP = SRP_WET
				rate0_POP = POP_WET
				rate0_DOP =	DOP_WET
			ELSE
				rate0_NO3 =	winter_NO3_DRY			
				rate0_NH4 = winter_NH4_DRY
				rate0_DON =	winter_DON_DRY
				rate0_PON =	winter_PN_DRY
				rate0_SRP = winter_SRP_DRY
				rate0_POP = winter_POP_DRY
				rate0_DOP =	winter_DOP_DRY
			ENDIF
		ELSEIF(jday.ge.91.and.jday.le.181) THEN	  !Spring  (April to June)
			IF(Rain.ge. 2.54)THEN	 !Wet days for rainfall>=0.1 inch or 2.54 mm
				rate0_NO3 =	NO3_WET			
				rate0_NH4 = NH4_WET
				rate0_DON =	DON_WET
				rate0_PON =	PN_WET
				rate0_SRP = SRP_WET
				rate0_POP = POP_WET
				rate0_DOP =	DOP_WET
			ELSE
				rate0_NO3 =	spring_NO3_DRY
				rate0_NH4 = spring_NH4_DRY
				rate0_DON =	spring_DON_DRY
				rate0_PON =	spring_PN_DRY			
				rate0_SRP = winter_SRP_DRY
				rate0_POP = winter_POP_DRY
				rate0_DOP =	winter_DOP_DRY
			ENDIF 
		ELSEIF(jday.ge.182.and.jday.le.273) THEN   !Summer (July to September)
			IF(Rain.ge. 2.54)THEN	 !Wet days for rainfall>=0.1 inch or 2.54 mm
				rate0_NO3 =	NO3_WET			
				rate0_NH4 = NH4_WET
				rate0_DON =	DON_WET
				rate0_PON =	PN_WET
				rate0_SRP = SRP_WET
				rate0_POP = POP_WET
				rate0_DOP =	DOP_WET	 
			ELSE
				rate0_NO3 =	summer_NO3_DRY
				rate0_NH4 = summer_NH4_DRY
				rate0_DON =	summer_DON_DRY
				rate0_PON =	summer_PN_DRY			
				rate0_SRP = summer_SRP_DRY
				rate0_POP = summer_POP_DRY
				rate0_DOP =	summer_DOP_DRY
			ENDIF
		ELSEIF(jday.ge.274.and.jday.le.310) THEN  !Fall
   			IF(Rain.ge. 2.54)THEN	 !Wet days for rainfall>=0.1 inch or 2.54 mm
				rate0_NO3 =	NO3_WET			
				rate0_NH4 = NH4_WET
				rate0_DON =	DON_WET
				rate0_PON =	PN_WET
				rate0_SRP = SRP_WET
				rate0_POP = POP_WET
				rate0_DOP =	DOP_WET
			ELSE
				rate0_NO3 =	fall_NO3_DRY
				rate0_NH4 = fall_NH4_DRY
				rate0_DON =	fall_DON_DRY
				rate0_PON =	fall_PN_DRY			
				rate0_SRP = summer_SRP_DRY
				rate0_POP = summer_POP_DRY
				rate0_DOP =	summer_DOP_DRY
			ENDIF

		ELSEIF(jday.ge.311.and.jday.le.334) THEN 	 !Fall	EXCLDUING DECEMBER FOR DRY P (October to November)
			IF(Rain.ge. 2.54)THEN	 !Wet days for rainfall>=0.1 inch or 2.54 mm
				rate0_NO3 =	NO3_WET			
				rate0_NH4 = NH4_WET
				rate0_DON =	DON_WET
				rate0_PON =	PN_WET
				rate0_SRP = SRP_WET
				rate0_POP = POP_WET
				rate0_DOP =	DOP_WET
			ELSE
				rate0_NO3 =	fall_NO3_DRY
				rate0_NH4 = fall_NH4_DRY
				rate0_DON =	fall_DON_DRY
				rate0_PON =	fall_PN_DRY			
				rate0_SRP = summer_SRP_DRY
				rate0_POP = summer_POP_DRY
				rate0_DOP =	summer_DOP_DRY
			ENDIF
		ELSE	!Only December
				IF(Rain.ge. 2.54)THEN	 !Wet days for rainfall>=0.1 inch or 2.54 mm
				rate0_NO3 =	NO3_WET			
				rate0_NH4 = NH4_WET
				rate0_DON =	DON_WET
				rate0_PON =	PN_WET
				rate0_SRP = SRP_WET
				rate0_POP = POP_WET
				rate0_DOP =	DOP_WET
			ELSE
				rate0_NO3 =	fall_NO3_DRY
				rate0_NH4 = fall_NH4_DRY
				rate0_DON =	fall_DON_DRY
				rate0_PON =	fall_PN_DRY			
				rate0_SRP = winter_SRP_DRY
				rate0_POP = winter_POP_DRY
				rate0_DOP =	winter_DOP_DRY
			ENDIF
		ENDIF
!********************Reduction **********************************************		
		rate1_SRP = redn_TP*rate0_SRP	!redn_TP	
		rate1_POP = redn_TP*rate0_POP	!redn_TP
		rate1_DOP =	redn_TP*rate0_DOP	!redn_TP

		rate1_NO3 =	redn_TN*rate0_NO3	!redn_TN
		rate1_NH4 = redn_TN*rate0_NH4	!redn_TN
		rate1_DON =	redn_TN*rate0_DON	!redn_TN
		rate1_PON =	redn_TN*rate0_PON	!redn_TN
!************CONVERSION FOR BAP AND BAN********************************
!   -------Bioavalable nitrogen and phosphorous----------------------
!    Refractive phosphrous is not bioavailable Quim's Dissertation				   
!   -------------BAP-------------------------
		rate_THP = 0.95*rate1_SRP + 0.30*rate1_POP + 0.10*rate1_DOP

		rate_SRP = rate1_SRP - 0.95*rate1_SRP				 
		rate_POP = rate1_POP - 0.30*rate1_POP	 !30%	 
		rate_DOP = rate1_DOP - 0.10*rate1_DOP	 !10%

!    so 10% of DON and PN is biavailable
		rate_NO3 = rate1_NO3 + 0.55*rate1_PON + 0.55*rate1_DON
		rate_NH4 = rate1_NH4
		rate_PON = rate1_PON - 0.55*rate1_PON 
		rate_DON = rate1_DON - 0.55*rate1_DON	
!****************END CONVERSION****************************************

		wqual(ns,16) = wqual(ns,16) + rate_THP/(surfacelayer_depth)		!    THP
		wqual(ns,15) = wqual(ns,15) + rate_SRP/(surfacelayer_depth)		!    SRP (wqual (i,12)
      wqual(ns,19) = wqual(ns,19) + rate_POP/(surfacelayer_depth)		!    POP
      wqual(ns,20) = wqual(ns,20) + rate_DOP/(surfacelayer_depth)		!    DOP
      wqual(ns,21) = wqual(ns,21) + rate_NO3/(surfacelayer_depth)		!    NO3
	   wqual(ns,22) = wqual(ns,22) + rate_NH4/(surfacelayer_depth)		!    NH4
		wqual(ns,25) = wqual(ns,25) + rate_PON/(surfacelayer_depth)		!    PON
		wqual(ns,26) = wqual(ns,26) + rate_DON/(surfacelayer_depth)		!    DON
! 
!    Nutrient Budget Control
! 
      N_Atm=(wqual(ns,21)+wqual(ns,22)+wqual(ns,25)+wqual(ns,26))*vol(ns)
      P_Atm=(wqual(ns,15)+wqual(ns,16)+wqual(ns,19)+wqual(ns,20))*vol(ns)   
       
    	RETURN
		END SUBROUTINE ATMOS_DEPOSITION
