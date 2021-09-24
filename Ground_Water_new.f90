	SUBROUTINE GROUND_WATER(jyear,redn_TP,redn_TN,N_Ground,P_Ground,vol_gw)
!***********************************************************************!
! 
!	This Routine calculates the Ground water inputs to Lake Tahoe
! It assumes ground water coming from the sides and the bottom of 
! the Lake. The temperature and salinity are assumed to be the same
! that the layer in which ground water inserts. 
! IF gw_flux is positive means input.
! The facts are taken from Table 9.3 and 4 of Lake Tahoe Basin Framework
!  ground water evaluation Lake Tahoe basin, California and Nevada (Final )
!  (October 2003) US Army Corps of Engineers Sacramento district.
! Called at the END of each simulated day,along with INFLOW and OUTFLOW
! Better models are expected when more will be known about ground water
!
!     written by Joa P. Losada and Geoff Schladow. UC-DAVIS 2000
!		Modified to linear Vol distributed, daily value 
! 
!***********************************************************************!	

! gw_flux: estimated water flow over day on each layer {Km3_day}
! Only dissolved substances are assumed to contribute to the nutrient balance
! The concentration of the WQ state variables are affected by the change of volume
! in each layer
! gw_THP {microg/L_day or mg/m3}
! gw_DOP  {microg/L_day or mg/m3} RP is now DOP
! gw_NO3 {microg/L_day or mg/m3}  
! gw_NH4 {microg/L_day or mg/m3} 
! gw_DON {microg/L_day or mg/m3}
! gw_Nitrogen: mean concentration of total Nitrogen {microg/L or mg/m3}
! gw_Posphorus: mean concentration of total Phosphorus {microg/L or mg/m3}
! gw_zero: set to zero the unaffected WQ variables
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE		
	
		INTEGER(4) :: i,j,k,jyear
		INTEGER*4 icode
		INTEGER*4 layno
      REAL(8):: gw_flux,vol_gw
!		REAL(8):: gw_THP,gw_POP, gw_NO3,gw_DON 
!     REAL(8):: gw_NH4
		REAL(8):: N_Ground,P_Ground	  
		REAL(8) :: vol_Total,vol_check
!		REAL(8), PARAMETER ::total_gw_flux = 64099.00 !(1000 m3/year)
     REAL(8), PARAMETER ::total_gw_flux = 4900.00  !(1000 m3/year) ! Value used in LCM SCT 27-08-2020
!      REAL(8), PARAMETER ::total_gw_flux = 10400.00 !(1000 m3/year)10400
      REAL(8), PARAMETER ::GMF=64099.00/total_gw_flux
		REAL(8), PARAMETER :: gw_Nitrogen  = 0.782025d3*GMF,gw_Posphorus = 0.106632d3*GMF  !1 mg/L = 1d3 microg/L or mg/m3 
     

!  gbs Acording to the report (2003) there is no depth and seaosonal variation
!	gbs	of phosphorous. so the values can be used as it is
!    
		REAL(8), PARAMETER :: gw_SRP  = 0.088535d3*GMF,gw_DOP  = 0.018097d3*GMF 		!1 mg/L = 1d3 microg/L or mg/m3

!gbs   Acording to the report there, nitrogen concentration of shallow 
!gbs	 ground water is 2-5 times higher than the deep aquifer. So take 
!gbs   shallow water nitrogen concentration 3 time hihger than deep. However,
!gbs   there is negligible seasonal nitrogen variation. Moreover, there 
!gbs   is not enough data to support the seasonal nitrogen variation.
 
 		REAL(8), PARAMETER ::gw_NO3_shallow = 0.378376d3*GMF,	&		!1 mg/L = 1d3 microg/L or mg/m3         		
     								gw_NH4_shallow = 0.069381d3*GMF,	&		!1 mg/L = 1d3 microg/L or mg/m3
									gw_DON_shallow = 0.138762d3*GMF		   !1 mg/L = 1d3 microg/L or mg/m3
    	REAL(8), PARAMETER ::gw_NO3_deep = 0.126125d3*GMF,		&		!1 mg/L = 1d3 microg/L or mg/m3
     								gw_NH4_deep = 0.023127d3*GMF,		&		!1 mg/L = 1d3 microg/L or mg/m3
									gw_DON_deep = 0.046254d3*GMF,		&		!1 mg/L = 1d3 microg/L or mg/m3
									gw_zero=0.0d0

		REAL(8) gw_NO3_shallow1, gw_NH4_shallow1, gw_DON_shallow1
		REAL(8) gw_NO3_deep1, gw_NH4_deep1, gw_DON_deep1
		REAL(8)	gw_DOP1,gw_SRP1

	
		REAL(8) gw_NO3_shallow_new, gw_NH4_shallow_new, gw_DON_shallow_new
		REAL(8) gw_NO3_deep_new, gw_NH4_deep_new, gw_DON_deep_new
		REAL(8)	gw_THP_new,gw_DOP_new,gw_SRP_new 		  
!IF		REAL(8):: combin,combinv,densty
		REAL(8), PARAMETER:: perecent_urban_N = 0.67,perecent_urban_P = 0.55 
		REAL*8 redn, redn_TP, redn_TN

!gbs  According to the report (2003)top GW discharge of top 40 ft is significantly
!gbs  higher than deep lake. Since no data is supported, take according to volume.
!gbs  Upper lake area and volume is higher than deeper. So, take accordingly. However,
!gbs  for depth greater than 15 m, take uniform discharge. Moreover, data supports only
!gbs for upper 110 m. So insert the GW discharge in the upper 110 m.

		WRITE (*, 10) jday,jyear,redn_TP, redn_TN
10		FORMAT(2i8,10x, 2f10.2, 4x,'Groundwater')
      vol_Total=0.0d0
		vol_check=0.0d0
!		WRITE(1000,*) jyear,jday,redn,'GW'
  
		DO 1 k = 1,ns
			IF (depth(k).GE.(depth(ns)-110)) THEN   !depth(1) is bottom and depth(ns) is surface
				vol_Total = vol_Total + vol(k)
			ELSE 
				GOTO 1
			ENDIF
!	print*,k, depth(k), vol(k),vol_Total
!	pause	
 1    CONTINUE
!    pause
! 
!gbs  Estimate the groundwater flux per time step (i.e., day)
! 
		vol_gw= total_gw_flux/(365.0)       
      N_Ground = vol_gw*(gw_Nitrogen)  
		P_Ground = vol_gw*(gw_Posphorus)
!************CONVERSION FOR BAP AND BAN********************************
!gbs ------------  Initialization  -------------------------------------
		gw_THP_new = 0.0
		gw_SRP_new = 0.0
		gw_DOP_new = 0.0
		gw_NO3_shallow_new = 0.0
		gw_DON_shallow_new = 0.0

		gw_NO3_deep_new = 0.0
		gw_DON_deep_new = 0.0
!gb ------------- END of Initialization ------------------------------
!gbs--------------Reduction rates-------------------------------------
!gbs---------------Phosphorous-----------------------------------------
		gw_SRP1 = redn_TP*gw_SRP
		gw_DOP1 = redn_TP*gw_DOP
!gbs---------------Nitrogen-------------------------------------------
		gw_NO3_shallow1 = redn_TN*gw_NO3_shallow    !ug/L or mg/m3
		gw_NH4_shallow1 = redn_TN*gw_NH4_shallow	!ug/L or mg/m3
		gw_DON_shallow1 = redn_TN*gw_DON_shallow    !ug/L or mg/m3 
		gw_NO3_deep1    = redn_TN*gw_NO3_deep 	    !ug/L or mg/m3
		gw_NH4_deep1    = redn_TN*gw_NH4_deep		!ug/L or mg/m3
		gw_DON_deep1    = redn_TN*gw_DON_deep 	    !ug/L or mg/m3

!gbs-------Bioavalable nitrogen and phosphorous----------------------	
!gbs-----------Phosphorous---------------------------------------------
		gw_THP_new=(0.95*gw_SRP1 +0.15*gw_DOP1)	    !ug/L or mg/m3	!95%SRP and 5-15% DOP
		gw_SRP_new=(gw_SRP1-0.95*gw_SRP1)	        !ug/L or mg/m3
		gw_DOP_new=(gw_DOP1-0.15*gw_DOP1)	        !ug/L or mg/m3 

!**************************************Urban and Non-urban Phosphorous **************************
!!	gw_THP_new=redn*perecent_urban_P*(0.95*gw_SRP +0.15*gw_DOP) +
!!     &            (1.0-perecent_urban_P)*(0.95*gw_SRP +0.15*gw_DOP)        	!ug/L or mg/m3	!95%SRP and 5-15% DOP

!!	gw_SRP_new=redn*perecent_urban_P*(gw_SRP-0.95*gw_SRP) +
!!     &            (1.0-perecent_urban_P)*(gw_SRP-0.95*gw_SRP)	    !ug/L or mg/m3

!!	gw_DOP_new=redn*perecent_urban_P*(gw_DOP-0.15*gw_DOP) +
!!     &            (1.0-perecent_urban_P)*(gw_DOP-0.15*gw_DOP)	    !ug/L or mg/m3 	   

!gbs-----------Nitrogen------------------------------------------------
!gbs so 10% of DON and PON is biavailable
		gw_NO3_shallow_new = gw_NO3_shallow1 + 0.55*gw_DON_shallow1 !ug/L or mg/m3
		gw_NH4_shallow_new = gw_NH4_shallow1						!ug/L or mg/m3
		gw_DON_shallow_new = gw_DON_shallow1 - 0.55*gw_DON_shallow1 !ug/L or mg/m3 
		gw_NO3_deep_new    = gw_NO3_deep1    + 0.55*gw_DON_deep1	!ug/L or mg/m3
		gw_NH4_deep_new    = gw_NH4_deep1							!ug/L or mg/m3
		gw_DON_deep_new    = gw_DON_deep1    - 0.55*gw_DON_deep1	!ug/L or mg/m3

!**************************************Urban and Non-urban Nitrogen **************************
!!	gw_NO3_shallow_new = redn*perecent_urban_N*
!!     &                         (gw_NO3_shallow + 0.80*gw_DON_shallow)+ 
!!     &  (1.0-perecent_urban_N)*(gw_NO3_shallow + 0.80*gw_DON_shallow) !ug/L or mg/m3
!!	gw_DON_shallow_new = redn*perecent_urban_N*(gw_DON_shallow - 
!!     & 	                 0.80*gw_DON_shallow)+ (1.0-perecent_urban_N)*
!!     &                     (gw_DON_shallow - 0.80*gw_DON_shallow)   !ug/L or mg/m3 
!!
!!	gw_NO3_deep_new    = redn*perecent_urban_N*(gw_NO3_deep + 
!!     &	                 0.80*gw_DON_deep)+ (1.0-perecent_urban_N)*
!!     &                     (gw_NO3_deep + 0.80*gw_DON_deep)	  !ug/L or mg/m3
!!	gw_DON_deep_new    = redn*perecent_urban_N*(gw_DON_deep - 
!!     &                	 0.80*gw_DON_deep)+ (1.0-perecent_urban_N)*
!!     &                     (gw_DON_deep - 0.80*gw_DON_deep)	  !ug/L or mg/m3
!gbs----------------------------------------------------------------
!****************END CONVERSION****************************************
! Add the volum and nutrients
! 
      DO k=1,ns
! 
! Calculate the volume of water to be inserted in upper layers 
! proportional to the vol of each layer
! 
			IF (depth(k).GE.(depth(ns)-110)) THEN	 !depth(1) is bottom and depth(ns) is surface
				gw_flux = vol(k)*vol_gw/vol_Total		  
			ELSE  
	 			gw_flux =	gw_zero 
			ENDIF
			vol_check=vol_check+gw_flux
!	gw_flux= vol_gw/dfloat(ns)
! 
! Temperature and salinity (unaffected)
! 
			temp(K)=COMBIN(temp(K),vol(K),den(K),temp(K),gw_flux,den(K))
			sal(K) =COMBIN(sal(K) ,vol(K),den(K),sal(K), gw_flux,den(K))
! 
! Particles
! 
			cf(K,1)=COMBINV(cf(K,1),vol(K),gw_zero,gw_flux)
			cf(K,2)=COMBINV(cf(K,2),vol(K),gw_zero,gw_flux)
			cf(K,3)=COMBINV(cf(K,3),vol(K),gw_zero,gw_flux)
			cf(K,4)=COMBINV(cf(K,4),vol(K),gw_zero,gw_flux)
			cf(K,5)=COMBINV(cf(K,5),vol(K),gw_zero,gw_flux)
			cf(K,6)=COMBINV(cf(K,6),vol(K),gw_zero,gw_flux)
			cf(K,7)=COMBINV(cf(K,7),vol(K),gw_zero,gw_flux)
! 
! Water Quality
! 
			wqual(K,1)= COMBINV(wqual(K,1),vol(K),gw_zero,gw_flux)
			wqual(K,2)= COMBINV(wqual(K,2),vol(K),gw_zero,gw_flux)
			wqual(K,3)= COMBINV(wqual(K,3),vol(K),gw_zero,gw_flux)
			wqual(K,4)= COMBINV(wqual(K,4),vol(K),gw_zero,gw_flux)
			wqual(K,5)= COMBINV(wqual(K,5),vol(K),gw_zero,gw_flux)
			wqual(K,6)= COMBINV(wqual(K,6),vol(K),gw_zero,gw_flux)
			wqual(K,7)= COMBINV(wqual(K,7),vol(K),gw_zero,gw_flux)
			
			
			wqual(K,23)=COMBINV(wqual(K,23),vol(K),gw_zero,gw_flux)
			wqual(K,24)=COMBINV(wqual(K,24),vol(K),gw_zero,gw_flux)
			wqual(K,25)=COMBINV(wqual(K,25),vol(K),gw_zero,gw_flux)
			wqual(K,26)=COMBINV(wqual(K,26),vol(K),gw_zero,gw_flux)
			wqual(K,27)=COMBINV(wqual(K,27),vol(K),gw_zero,gw_flux)
			wqual(K,28)=COMBINV(wqual(K,28),vol(K),gw_zero,gw_flux)
			IF (depth(k).GE.(depth(ns)-15)) THEN 
!  Shallow aquifer nitrogen concentration is 3 times higher than deep 
!  depth(1) is bottom and depth(ns) is surface
	  
				wqual(K,14)=COMBINV(wqual(K,14),vol(K),gw_THP_new,gw_flux) !THP
				wqual(K,13)=COMBINV(wqual(K,13),vol(K),gw_SRP_new,gw_flux) !SRP
				wqual(K,16)=COMBINV(wqual(K,16),vol(K),gw_DOP_new,gw_flux) !DOP
				wqual(K,18)=COMBINV(wqual(K,18),vol(K),gw_NO3_shallow_new,gw_flux) !NO3
				wqual(K,19)=COMBINV(wqual(K,19),vol(K),gw_NH4_shallow_new,gw_flux) !NH4
				wqual(K,21)=COMBINV(wqual(K,21),vol(K),gw_DON_shallow_new,gw_flux) !DON
			ELSE	!deep aquifer	  
				wqual(K,14)=COMBINV(wqual(K,14),vol(K),gw_THP_new,gw_flux) !THP
				wqual(K,13)=COMBINV(wqual(K,13),vol(K),gw_SRP_new,gw_flux) !SRP
				wqual(K,16)=COMBINV(wqual(K,16),vol(K),gw_DOP_new,gw_flux) !DOP
				wqual(K,18)=COMBINV(wqual(K,18),vol(K),gw_NO3_deep_new,gw_flux) !NO3
				wqual(K,19)=COMBINV(wqual(K,19),vol(K),gw_NH4_deep_new,gw_flux) !NH4
				wqual(K,21)=COMBINV(wqual(K,21),vol(K),gw_DON_deep_new,gw_flux) !DON
			ENDIF
! 
! Density and Volume layers
! 
			den(K)=DENSTY(temp(K),sal(K))
			vol(K)=vol(K)+gw_flux
	! 
! Adjust volumes
! 
			IF (K .eq. 1)THEN
				vol1(K)=vol(K)
			ELSE
				vol1(K)=vol1(K-1) + vol(K)
			ENDIF
      ENDDO ! k
!	WRITE(21,fmt='(i5,6f15.5)')jday,vol_gw*1000.0,vol_gw*gw_SRP1/1000.0,vol_gw*gw_DOP1/1000.0,vol_gw*(gw_NO3_shallow1+gw_NO3_deep1)/1000.0, &
!                  vol_gw*(gw_NH4_shallow1+gw_NH4_deep1)/1000.0,vol_gw*(gw_DON_shallow1+gw_DON_deep1)/1000.0

!
!  Make adjustments to correct layer volumes.
! 
		icode = 2     !ARRAYS OF VOLUME AND AREAS FROM DEPTHS (icode=1)!ARRAYS OF DEPTHS AND AREAS FROM VOLUME (icode=2)				   
		layno = 1     !Staring layer number

      CALL RESINT(icode,layno) 
! 
! Check layers for vmax,vmin
! 
      CALL THICK

		RETURN
		END SUBROUTINE GROUND_WATER
