!*********************************************************************************
		SUBROUTINE NUTRIT_TAHOE(deltnho,Al_THP,Al_NH4,Al_NO3,N_anoxic)
!*********************************************************************************
!     are(maxns) = area of sediment for a layer
!     arfac = multiplicative factor for area adjustment
!     cn1 = rate of PON ==> NH4
!     cn2 = rate of PON ==> DON
!     cn3 = rate of DON ===> NH4
!     cn4 = rate of NO3 ==> 2NH3+3O2
!     cn5 = rate of NO3 ==>  N2(g)
!     cp1 = rate of POP ==> THP  
!     cp2 = rate of POP ==> RP 
!     cp3 = rate of RP ===> THP
!     deltnho = release of oxygen from nitrate utilisation (mg l**-1)
!     dep = change in depth between layers
!     death = chlorophyll lost by death and respiration
!     faco = factor (0 or 1) for reducing biological processes to zero when
!            oxygen concentrations are low (faco=0), set in subr. algae
!     growth = chlorophyll grow
!     grazing = chlorophyll lost by zoo grazing
!     hscnit = half saturation constant for limitation of nitrification by DO
!     hscden = half saturation constant for limitation of denitrification by DO
!     hscnhno = ammonia preference factor
!     jday = julian day
!     nchl = number of chla groups modelled
!     sedpo4 = sediment release rate of po4
!     sedp = adjusted sediment release rate of p that accounts for log
!     sednh3 = sediment release rate of NH3
!     sedno3 = sediment release rate of NO3
!     sedtem = sediment temperature multiplier for nutrient release
!     thetabo = temperature multiplier for decay and desorption
!     thetanh = temperature multiplier for NH4 ==> NO3 and NO3 ==> N2
!     thetat  = temperature multiplier for phytoplankton growth
!     ypchla = ratio of P to chlorophyll a
!     ynchla = ratio of N to chlorophyll a
!     ysichla = ratio of Si to chlorophyll a  (could include IF silica is important)
!     Total Phosphorus = THP + PLP + PDP + RP
!    	wqual(i,12) = Total Hidrolizable Phosphorus, [ug/L]=[mg/m3]
!    	wqual(i,11),wqual(i,12) = Particulate Living Phosphorus-1, [ug/L]=[mg/m3]
!    	wqual(i,13) = Particulate Detritus Phosphorus, [ug/L]=[mg/m3]
!    	wqual(i,14) = Refractive Phosphorus, [ug/L]=[mg/m3]
!     Total Nitrogen = NO3 + NH3 + PLN + PDN + DON
!    	wqual(i,15) = Nitrate NO3, [ug/L]=[mg/m3]
!    	wqual(i,16) = Ammonia NH3, [ug/L]=[mg/m3]
!    	wqual(i,17),wqual(i,18) = Particulate Living Nitrogen, [ug/L]=[mg/m3]
!    	wqual(i,19) = Particulate Detritus Nitrogen, [ug/L]=[mg/m3]
!    	wqual(i,20) = Dissolved Organic Nitrogen, [ug/L]=[mg/m3]
!  
!   **********************************************************************
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	
		REAL*8 are(maxns)
		REAL*8 arfac, volfac
		REAL*8 Control_Excre
		REAL*8 dep(maxns)
		REAL*8 deltnho(maxns)
		REAL*8 factor_Caca
!gbs---------------------------------------------------------------
		REAL*8 POP_to_THP, THP_from_POP,PON_to_NH4,NH4_from_PON
		REAL*8 POP_to_RP,RP_from_POP, PON_to_DON, DON_from_PON
		REAL*8 DOP_to_THP,THP_from_DOP,DON_to_NH4,NH4_from_DON
		REAL*8 NH4_to_NO3, NO3_from_NH4, NO3_to_N2
!gbs---------------------------------------------------------------
		REAL*8 C_uptake, P_Uptake,N_Uptake,Si_uptake, C_RandM, P_RandM, N_RandM,Si_RandM,C_Grazing, P_Grazing,N_Grazing, Si_Grazing
!gbs---------------------------------------------------------------
		REAL*8 N_anoxic
		REAL*8 hscden
		REAL*8 hscnhno
		REAL*8 NH4sedload	! NH4 load from sediment
		REAL*8 THPsedload	! THP load from sediment
		REAL*8 NO3sedload   ! NO3 load from sediment
		REAL*8 TtVfactor	! lumps in the temperature, time interval
				        ! and layer volume factors
      REAL*8 VolTotal,VolBottom
		REAL*8 Al_THP,Al_NH4, Al_NO3 !Nutrient Budget
		REAL*8 kdop
      INTEGER*4 control_time
		INTEGER*4 i
		INTEGER*4 j
!IF		INTEGER*4 jday
!IF		INTEGER*4 nchl

!   set half saturation constants for oxygen limitation
!   of nitrification and denitrification... mg/l
		PARAMETER (arfac = 1.0d+6, volfac=1.0d+3)
		PARAMETER (hscden = 100.0d0)
		PARAMETER (factor_Caca = 1.0d0)
!gbs---------------------------------------------------------------
      REAL*8 SOD_area (maxns), SOD_radius (maxns), SUM_SOD_area
		REAL*8 SOD_radius_L (maxns),SOD_radius_U (maxns)
		REAL*8 slant_height(maxns)
		REAL(8), PARAMETER :: pi=3.1415926535897932384626433832795d0
		REAL(8), PARAMETER :: KHnnt=5.0d0     !1.00 Cerco
		REAL(8), PARAMETER :: KHndn=0.1d0     !cerco
!gbs----------------------------------------------------------------
		Al_THP   = 0.0
		Al_NH4   = 0.0
		Al_NO3   = 0.0
		N_anoxic = 0.0
!  
!     Check for Nutrient < 0
!  	
		DO i = 1,ns
			IF(wqual(i,17).lt. 0.0) THEN
				WRITE(*,*)'Warning PhytoP < 0 Start of Nutrit',i,wqual(i,15)		!POP
!				PAUSE 
			ENDIF
			IF(wqual(i,23).lt. 0.0) THEN
				WRITE(*,*)'Warning PhytoN < 0 Start of Nutrit',i,wqual(i,20)		!PON
!				PAUSE 
			ENDIF
		ENDDO   ! i
!  
!     Calculate total Volume
!  
      VolTotal  = 0.0
      VolBottom = 0.0
       	DO i = 1, ns
			IF(depth(i).lt.depth(ns)-150) VolBottom = VolBottom + vol(i)*volfac
		ENDDO   ! i
!  
!     Initialize Zooplankton variables
!  

      IF(iclock.ge.43200.and.iclock+nosecs.le.86400) THEN  ! Day conditions
			Shit_P = nosecs*(Caca_P/(VolBottom*43200))
			Shit_N = nosecs*(Caca_N/(VolBottom*43200))
			Control_Excre = Control_Excre + Shit_P*VolBottom
      ELSE
			Shit_P = 0.0
			Shit_N = 0.0
		ENDIF
9194  CONTINUE

!  
!     First, work out the BIOTIC components
!  
        	
		DO i = 1,ns
!	    deltnho(i) = wqual(i,15)	! deltnho is used in oxygen  
!     Preferential ammonium uptake factor (pg 270 Bowie et al.1985)  
!		hscnhno = wqual(i,16)/(wqual(i,16)+wqual(i,15))
			DO j = 1,nchl
!
!gbs ------------Nutrient uptake kinetics. Used up so should be deducted --
                C_Uptake = ycchla(j)*growth(i,j)
				P_Uptake = ypchla(j)*growth(i,j)
				N_Uptake = ynchla(j)*growth(i,j)
                Si_Uptake = ySIchla(j)*growth(i,j)
                if(j.eq.4) N_Uptake=0.0

                !
!gbs-------------Respiration and mortality. Practically released --
!  
                C_RandM = ycchla(j)*death(i,j)
                P_RandM = ypchla(j)*death(i,j)
				N_RandM = ynchla(j)*death(i,j)
                Si_RandM = ysichla(j)*death(i,j)
!  
!gbs-----Grazing phytoplsnkton by zooplankton, mysis etc. So loss. should be duducted--
!  
				P_Grazing = ypchla(j)*grazing(i,j)
				N_Grazing = ynchla(j)*grazing(i,j)
!  
!     Preferential ammonium uptake factor..(pg 270 Bowie et al.1985).
!
				IF((wqual(i,21) + wqual(i,22)).le. 0.0) THEN
					hscnhno = 0.0 
				ELSE
					hscnhno = wqual(i,22)*wqual(i,21)/((hscn(j)+wqual(i,22))*(hscn(j)+wqual(i,21)))	&		
     					+wqual(i,22)*hscn(j) / ((wqual(i,21) +wqual(i,22))*(hscn(j)+wqual(i,21)))
	!			   hscnhno=(NH4*NO3)/((hs+NH4)*(hs+NO3))+NH4*Hs/((NH4+NO3)*(hs+NO3))
				ENDIF
!     Carbon
               	IF(wqual(i,10).lt.C_Uptake) THEN
					WRITE(*,*)'Warning. C < 0. caused by Biotic fators of Nutrit'
!					PAUSE
				ENDIF

                wqual(i,10)   = wqual(i,10)   - C_Uptake + 0.90*C_RandM
!   Silica
                IF(wqual(i,27).lt.Si_Uptake) THEN
					WRITE(*,*)'Warning. Si < 0. caused by Biotic fators of Nutrit'
!					PAUSE     
				ENDIF

                wqual(i,27)   = wqual(i,27)   - Si_Uptake + 0.90*Si_RandM

!     Phosphorus
				IF(wqual(i,16).lt.P_Uptake) THEN
!					WRITE(*,*)'Warning. THP < 0. caused by Biotic fators of Nutrit'
!					PAUSE
				ENDIF

				IF((wqual(i,17)+P_Uptake).lt.P_RandM) THEN
!					WRITE(*,*)'Warning. PhytoP < 0. caused by Biotic fators of Nutrit'
					WRITE(*,*) wqual(i,17)+ P_Uptake - P_RandM
!					PAUSE
                ENDIF
 				wqual(i,16)   = wqual(i,16)   - P_Uptake + 0.15*P_RandM	 !15%  since POP 					                                  
				wqual(i,17)   = wqual(i,17)   - P_Grazing +0.10*P_RandM
				wqual(i,18)   = wqual(i,18)   + 0.10* P_RandM 
!  
!     Nitrogen
				IF(wqual(i,21).lt.N_Uptake*(1-hscnhno)) THEN
!					WRITE(*,*)'Warning. NO3 < 0. caused by Biotic factors of Nutrit',		&
!							 hscnhno,wqual(i,21),N_Uptake*(1-hscnhno)
!					PAUSE
				ENDIF
				IF(wqual(i,22).lt.N_Uptake*hscnhno) THEN
!					WRITE(*,*)'Warning. NH4 < 0. caused by Biotic factors of Nutrit'
!					PAUSE
				ENDIF

				IF((wqual(i,23)+N_Uptake).lt. N_RandM) THEN
!					WRITE(*,*)'Warning. PhytoN < 0.caused by Biotic factors of Nutrit'
!					PAUSE
				ENDIF

				wqual(i,21)   = wqual(i,21)- N_Uptake*(1-hscnhno)+0.05*N_RandM 	 !PON	 65%
				wqual(i,22)   = wqual(i,22) - N_Uptake*hscnhno +0.95*N_RandM
!				wqual(i,16+j) = wqual(i,16+j)+N_Uptake-N_RandM-N_Grazing	
				wqual(i,23)   = wqual(i,23)   + N_RandM - N_Grazing -	1.00*N_RandM  !PON		65%
!  
!     Zooplankton effect
!  
				IF(itimes.eq.1440) THEN
!					wqual(i,19)   = wqual(i,19)  + N_Grazing
!					wqual(i,13)   = wqual(i,13)  + P_Grazing
					GOTO 5033
				ENDIF   	
				IF(iclock.ge.43200.and.iclock+nosecs.le.86400)THEN  ! Day conditions

					IF(depth(i).lt.(depth(ns)-150)) THEN               ! Bottom layers
!						wqual(i,19) = wqual(i,19) + N_Grazing + Shit_N         
!						wqual(i,13) = wqual(i,13) + P_Grazing + Shit_P
					ELSE	  
!						wqual(i,19) = wqual(i,19) + N_Grazing              ! Top layers
!						wqual(i,13) = wqual(i,13) + P_Grazing           
					ENDIF
				ELSE                                               ! Night conditions
					IF(depth(i).lt.(depth(ns)-150)) THEN               ! Bottom layers
!						wqual(i,19) = wqual(i,19) + N_Grazing              ! Top layers
!						wqual(i,13) = wqual(i,13) + P_Grazing        
					ELSE
						Caca_P = Caca_P + factor_Caca*P_Grazing*vol(i)*volfac
						Caca_N = Caca_N + factor_Caca*N_Grazing*vol(i)*volfac	 	
!						wqual(i,19) = wqual(i,19) + (1-factor_Caca)*N_Grazing  ! Top layers
!						wqual(i,13) = wqual(i,13) + (1-factor_Caca)*P_Grazing
					ENDIF
				ENDIF
5033			CONTINUE
			ENDDO	! j
        ENDDO	   ! i
        
!  
!     Check consistency of the results
!  
		DO i = 1,ns
			IF(wqual(i,17).lt. 0.0) THEN
!				WRITE(*,*)'Warning PhytoP < 0 Mid of Nutrit',i,wqual(i,17)
!				PAUSE 
			ENDIF
			IF(wqual(i,23).lt. 0.0) THEN
!				WRITE(*,*)'Warning PhytoN < 0 Mid of Nutrit',i,wqual(i,23)
!				PAUSE 
			ENDIF
		ENDDO
!gbs*************************************************************
!gbs***************CONVERSION OF SPECIES FROM ONE TO OTHER*******
!gbs*************************************************************

!     Second, work out the ABIOTIC factors 
!gbs-------------cp1 and cn1-----------------------------
!gbs cp1 and cn1 portion of POP and PON converts to THP and NH4,respectively.  
!gbs POP ==> THP and PON ==> NH4
!  
		DO i = 1,ns
			POP_to_THP = cp1*thetabo**(temp(i)-20.0)*wqual(i,17)*(nosecs/86400.0)
			PON_to_NH4 = cn1*thetabo**(temp(i)-20.0)*wqual(i,23)*(nosecs/86400.0)
			IF(wqual(i,17).lt. POP_to_THP) THEN
!				WRITE(*,*)'Warning. POP < 0. caused by Abiotic CP1 in Nutrit',i,j
				POP_to_THP = 0.0
!				PAUSE
			ENDIF

			IF(wqual(i,23).lt.PON_to_NH4) THEN
!				WRITE(*,*)'Warning. PON < 0. caused by Abiotic CN1 in Nutrit',i,j
				PON_to_NH4 = 0.0
!				PAUSE
			ENDIF

			THP_from_POP = POP_to_THP
			NH4_from_PON = PON_to_NH4

			wqual(i,17) = wqual(i,17) - POP_to_THP
			wqual(i,16) = wqual(i,16) + THP_from_POP

			wqual(i,23) = wqual(i,23) - PON_to_NH4
			wqual(i,22) = wqual(i,22) + NH4_from_PON
		ENDDO	   ! i
!gbs-------------cp2 and cn2-------------------------------------------
!gbs cp2 and cn2 portion of POP and PON converts to RP and DON, respectively 
!gbs POP ==> RP and PON ==> DON
!  
		DO i = 1,ns
			POP_to_RP  = cp2*thetabo**(temp(i)-20.0)*wqual(i,17)*(nosecs/86400.0)
			PON_to_DON = cn2*thetabo**(temp(i)-20.0)*wqual(i,23)*(nosecs/86400.0)

			IF(wqual(i,17).lt.POP_to_RP) THEN
!				WRITE(*,*)'Warning. POP < 0. caused by Abiotic CP2 in Nutrit',i,j
				POP_to_RP = 0.0
!				PAUSE
			ENDIF

			IF(wqual(i,23).lt.PON_to_DON) THEN
!				WRITE(*,*)'Warning. PON < 0. caused by Abiotic CN2 in Nutrit',i,j
				PON_to_DON = 0.0
!				PAUSE
			ENDIF

			RP_from_POP = POP_to_RP
      	    DON_from_PON= PON_to_DON

			wqual(i,17) = wqual(i,17) - POP_to_RP
			wqual(i,19) = wqual(i,19) + RP_from_POP

			wqual(i,23) = wqual(i,23) - PON_to_DON
			wqual(i,24) = wqual(i,24) + DON_from_PON
		ENDDO	   ! i
!gbs------------------cp3 and cn3-------------------------------------------
!gbs cp3 and cn3 portion of RP and DON converts to THP and NH4, respectively
!gbs RP ===> THP and DON ===> NH4 
!  
		DO i = 1, ns
			DO j = 1, nchl
!				kdop = 0.12*(nosecs/86400.0)+(halfp(j)/(halfp(j)+wqual(i,10)))*	&
!                   (0.01*(86400.0/nosecs))*wqual(i,j)   !0.12,	0.01
!				DOP_to_THP = kdop*exp(0.032*(temp(i)-20.0))*wqual(i,10+j)

!				DO i = 1,ns
					DOP_to_THP = cp3*thetabo**(temp(i)-20.0)*wqual(i,18)*(nosecs/86400.0)
					DON_to_NH4 = cn3*thetabo**(temp(i)-20.0)*wqual(i,24)*(nosecs/86400.0)
     	
					IF(wqual(i,14).lt.DOP_to_THP) THEN
!						WRITE(*,*)'Warning. RP < 0. caused by Abiotic CP3 in Nutrit',i,j
						DOP_to_THP = 0.0
!						stop
					ENDIF

					IF(wqual(i,20).lt.DON_to_NH4) THEN
!						WRITE(*,*)'Warning. DON < 0. caused by Abiotic CN3 in Nutrit',i,j
						DON_to_NH4 = 0.0
!						PAUSE
					ENDIF
		
					THP_from_DOP= DOP_to_THP
					NH4_from_DON=DON_to_NH4

					wqual(i,18) = wqual(i,18) - DOP_to_THP
					wqual(i,16) = wqual(i,16) + THP_from_DOP
		
					wqual(i,24) = wqual(i,24) - DON_to_NH4
				    wqual(i,22) = wqual(i,22) + NH4_from_DON
			ENDDO     !j		
		ENDDO	   ! i
!gbs----------------cn4 and cn5--------------------------------------
!gbs cn4 portion of NH4 converts to NO3
!gbs cn4 portion of NO3 converts to N2 gas
!gbs NH3 ==> 2N03 (Oxic) and NO3 ==>  N2(g) (Anoxic).  
!gbs determine the release of oxygen from nitrate uptake for use in the
!gbs oxygen SUBROUTINE: 2NO3 -> 3N2+3O2. 
!  
		DO i = 1,ns
!			deltnho(i) = (deltnho(i)-wqual(i,15))*48.0/14.0/1000.0
!			IF (faco(i).ge.0.5) THEN	! oxic
!gbs---------Nitrification process when oxygen is abundant------------------- 
         IF (wqual(i,8).ge.0.5) THEN   ! oxic
				NH4_to_NO3 = cn4*thetanh**(temp(i)-20.0)							&
     			* (wqual(i,8)/(wqual(i,8)+hscnit))									&
     			* (wqual(i,21)/(wqual(i,22)+KHnnt))*DFLOAT(nosecs)/86400.0

				IF(wqual(i,22).lt.NH4_to_NO3) THEN
!					WRITE(*,*)'Warning. Oxic Conditions. NH4 < 0 in Nutrit'
!					PAUSE
					NH4_to_NO3=0.9*wqual(i,22)
				ENDIF
!				NO3_from_NH4 = NH4_to_NO3
!                NH4_to_NO3=0.0
				wqual(i,22) = wqual(i,22) - NH4_to_NO3		   
				wqual(i,21) = wqual(i,21) +  NH4_to_NO3 		    
!gbs---Denitrification in aerobic condition specifically for lake Tahoe --------------
				IF(wqual(i,21).gt.3.0d0.and.temp(i).gt.9.0) THEN
	!				NO3_to_N2 = cn5*thetanh**(temp(i)-20.0)		&
	!     		*(wqual(i,8)/(wqual(i,8)+hscden))*wqual(i,15)*nosecs/86400.0
				ELSE
	!	         NO3_to_N2=0.0  
				ENDIF
         ELSEIF(wqual(i,8).lt.0.5) THEN						! Anoxic
	       
!1060  format(I10,x,I3,x,f10.3)          
!gbs---------De-nitrification process when oxygen is depleted------------------- 
			NO3_to_N2 = cn5*thetanh**(temp(i)-20.0)								&
     			*(wqual(i,8)/(wqual(i,8)+hscden))*									&
            (wqual(i,21)/(wqual(i,21)+KHndn))*DFLOAT(nosecs)/86400.0    		
            NO3_to_N2=0.0
			IF(wqual(i,21).lt.NO3_to_N2) THEN
				NO3_to_N2=0.0d0
			ENDIF
!				wqual(i,21) = wqual(i,21) - NO3_to_N2
			ENDIF
		ENDDO	   ! i

!gbs*******************************************************************
!gbs***********************END CONVERSION******************************
!---------------------------------------------------------------------
!gbs2010*******SEDIMENT RELEASE OF P AND NH4 IN CASE OF ANIXIA ************
!-------------------------------------------------------------------
!	TtVfactor is common to all layers 
!   includes T, time, and V effects
!   
!gbs changes the SOD area
		DO i = ns,2,-1		
			dep(i) = depth(i)-depth(i-1)
			SOD_radius (i)= DSQRT((0.5*(area(i)+area(i-1))*arfac)/PI)
			SOD_radius_L(i)=DSQRT(area(i-1)*arfac/PI)
			SOD_radius_U(i)=DSQRT(area(i)*arfac/PI)
			slant_height(i)=SQRT(dep(i)*dep(i)+(SOD_radius_U(i)-SOD_radius_L(i))**2)
			SOD_area   (i)= 2.0D0*PI*SOD_radius (i)*slant_height(i)			
		ENDDO

		are(1) = area(1)*arfac
		dep(1) = depth(1)
		SOD_area(1)=are(1)
!gbs---------------------------------         
		sum_SOD_area=0.0
		DO i = 1,ns
!-------Initialization---------
			THPsedload=0.0d0
			NH4sedload=0.0d0 
			NO3sedload=0.0d0
!			TtVfactor   = sedtem**(temp(i)-20.0)*(float(nosecs)/86400.0)*SOD_area(i)
			TtVfactor   = SOD_area(i)*(float(nosecs)/86400.0)    
!----------------------------- 
			IF (wqual(i,8).le.0.01) THEN  ! anoxic
!			IF(depth(i).le.200.0) THEN
         
!gbs     changed the old code because of SOD area is different
!gbs			TtVfactor   = sedtem**(temp(i)-20.0)*(float(nosecs)/86400.0)*		&
!gbs     					 (are(i)/(area(i)*arfac*dep(i)))				
				THPsedload  = sedthp * TtVfactor ! microg because THP from sediments/microg/m2/day
				NH4sedload  = sednh3 * TtVfactor ! microg because NH4 from sediments/microg/m2/day
				NO3sedload  = 0.0d0
		
				wqual(i,16) = wqual(i,16) + THPsedload/(vol(i)*volfac)	! THP  microg to microg/m3
				wqual(i,21) = wqual(i,21) + 0.0d0
				wqual(i,22) = wqual(i,22) + NH4sedload/(vol(i)*volfac)	! NH4  microg to microg/m3

			ELSEIF(wqual(i,8).gt.0.01) THEN     
				THPsedload =  0.0d0
				NH4sedload =  0.0d0
				NO3sedload =  sedno3 * TtVfactor ! microg because NO3 from sediments/microg/m2/day	 	
				wqual(i,16) = wqual(i,16) + 0.0d0
				wqual(i,21) = wqual(i,21) + NO3sedload/(vol(i)*volfac)
				wqual(i,22) = wqual(i,22) + 0.0d0
			ENDIF	   ! i
			Al_THP = Al_THP + THPsedload     ! SRP budget
			Al_NH4 = Al_NH4 + NH4sedload     ! NH4 budget	  
			Al_NO3 = Al_NO3 + NO3sedload     ! NO3 budget
!		print*,i,depth(i), SOD_area(i)
!		print*,area(i)*arfac,sum_SOD_area 
!		PAUSE
        ENDDO
        
!------------------for verification-------------------------------
!		WRITE(105, fmt='(i6, 3f20.2)') jday,Al_THP, Al_NH4,Al_NO3
!  ----------------------------------------------------------------
!       note that the effect of settling on the particulate P and N is
!        handled by particles
!
!gbs*********************************************************************
		CALL NEUTRIENT_SETTLING 
!gbs*********************************************************************
!      Check consistency of the results
!  
		DO i = 1,ns
			IF(wqual(i,17).lt. 0.0) THEN
				WRITE(*,*)'Warning PhytoP < 0 End of Nutrit',i,wqual(i,17)
				PAUSE 
			ENDIF
			IF(wqual(i,23).lt. 0.0) THEN
				WRITE(*,*)'Warning PhytoN < 0 End of Nutrit',i,wqual(i,23)
				PAUSE 
			ENDIF
		ENDDO

		DO i = 1,ns
			IF (wqual(i,16).lt.0.0) THEN 
				WRITE(*,*)'Warning: THP less than zero at Jday',jday,'layer',i
				wqual(i,16) = 0.0
				PAUSE
			ENDIF	      
			IF (wqual(i,21).lt.0.0) THEN
				WRITE(*,*)'Warning: NO3 less than zero at Jday',jday,'layer',i
				wqual(i,21) = 0.0
				PAUSE
			ENDIF	

			IF (wqual(i,22).lt.0.0) THEN
				WRITE(*,*)'Warning: NH4 less than zero at Jday',jday,'layer',i
				wqual(i,22) = 0.0
				PAUSE
			ENDIF
!  
!     keep totals above inorg. fractions...
!  
!			IF (wqual(i,14).le.wqual(i,12)) wqual(i,14)=wqual(i,12)+0.001
!			IF (wqual(i,20).le.(wqual(i,15)+wqual(i,16)))  &
!			    wqual(i,20) = wqual(i,15)+wqual(i,16)+0.001
		ENDDO

      IF(iclock+nosecs.eq.86400) THEN
			caca_P = 0.0
			caca_N = 0.0
			Shit_P = 0.0
			Shit_N = 0.0
			Control_Excre = 0.0
		ENDIF
	
		RETURN
		END SUBROUTINE NUTRIT_TAHOE	  
!       
!***************Neutrints settling**********************************
		SUBROUTINE NEUTRIENT_SETTLING 
!*******************************************************************
!gbs written by Goloka Sahoo, 2006
!gbs  Phosphorous Neutrients
!     wqual(i,12) = THP  [Total Hydrolizable phosphorous (mg/m3 or ug/L)]
!     wqual(i,11) = PLP1 or INP1 [Particulate Living phosphorous or internal phytoplnakton P 1 (mg/m3 or ug/L)]
!     wqual(i,12) = PLP2 or INP2 [Particulate Living phosphorous or internal phytoplnakton P 2 (mg/m3 or ug/L)]
!     wqual(i,13) = POP or PDP  [Particulate Organic or Detrius Phosphorous	(mg/m3 or ug/L)]
!     wqual(i,14) = RP   [Refractive Phosphorous (mg/m3 or ug/L)]
!     Total Phosphorous (TP) = THP + PLP1+PLP2+POP/PDP+RP
!gbs  Nitrogen Neutrients
!     wqual(i,15) = NO3 [Nitrate (mg/m3 or ug/L)]
!     wqual(i,16) = NH4 [Ammonia (mg/m3 or ug/L)]
!     wqual(i,17) = PLN1 or INN1 [Particulate Living nitrogen or internal phytoplnakton N 1 (mg/m3 or ug/L)]
!     wqual(i,18) = PLN2 or INN2 [Particulate Living nitrogen or internal phytoplnakton N 2 (mg/m3 or ug/L)]
!     wqual(i,19) = PON or PDN   [Particulate Organic or Detrius Nitrogen	(mg/m3 or ug/L)]
!     wqual(i,20) = DON          [Dissolved organic nitrogen (mg/m3 or ug/L)]
!     Total Nitrogen (TN) = NO3+NH4+PLN1+PLN2+PON+DON	

      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	
	
		INTEGER*4, PARAMETER:: wqpara=100
		INTEGER*4 i,j,k
		REAL*8 P_neutrient_settl_rate
		REAL*8 N_neutrient_settl_rate
		REAL*8 setldep(maxns, wqpara)
		REAL*8 wqualdv(maxns,wqpara), fluxin(maxns,wqpara)
		REAL*8 dep(maxns)
		REAL*8 depbotom(maxns,wqpara),deptop(maxns,wqpara)
		PARAMETER (P_neutrient_settl_rate = 0.0449315) !16.4 m/yr or 0.0449315 m/d
		PARAMETER (N_neutrient_settl_rate = 0.0328767) !12.0 m/yr or 0.0328767 m/d
		REAL(8), PARAMETER :: BOD_settl_rate = 0.08d0 !	m/d	same like Phytoplankton
		REAL(8), PARAMETER :: pi=3.1415926535897932384626433832795d0
		REAL(8), PARAMETER ::volfac = 1.0d+3,arfac=1.0D+6

!gbs depth of each layer
		DO i = ns,2,-1
			dep(i) = depth(i)-depth(i-1)
		ENDDO
		dep(1) = depth(1)  

		DO i =1, NS        !NS IS SURFACE AND 1 IS BOTTOM
! For Phosphorous neutrients
			setldep(i,15)	= 0.0 !nosecs*P_neutrient_settl_rate/86400.0d0	!THP
            setldep(i,16)	= 0.0	!PLP1
			setldep(i,17)	= nosecs*P_neutrient_settl_rate/86400.0d0	!PP or POP
			setldep(i,18)	= 0.0  !nosecs*P_neutrient_settl_rate/86400.0d0	!PP or POP 
			setldep(i,19)	= nosecs*P_neutrient_settl_rate/86400.0d0	!PP or POP
            setldep(i,20)	= 0.0
! For Nitrogen neutrients
			setldep(i,21)	= 0.0 !nosecs*N_neutrient_settl_rate/86400.0d0	!NO3
			setldep(i,22)	= 0.0 !nosecs*N_neutrient_settl_rate/86400.0d0	!NH4
			setldep(i,23)	= nosecs*N_neutrient_settl_rate/86400.0d0  
			setldep(i,24)	= 0.0 !Only particulate not dissolved are settled down
			setldep(i,25)	= nosecs*N_neutrient_settl_rate/86400.0d0	!PLN1
			setldep(i,26)	= 0.0	!PLN2
!			setldep(i,22)	= nosecs*N_neutrient_settl_rate/86400.0d0	!PN 
!			setldep(i,20)	= nosecs*N_neutrient_settl_rate/86400.0d0 	!DON
			
! FOR BOD Settling
			setldep(i,9)	= nosecs*BOD_settl_rate/86400.0d0	!BOD
            setldep(i,10)	= 0.0 !nosecs*BOD_settl_rate/86400.0d0	!BOD
            setldep(i,11)	= nosecs*BOD_settl_rate/86400.0d0
            setldep(i,12)	= 0.0
            setldep(i,13)	= nosecs*BOD_settl_rate/86400.0d0
            setldep(i,14)	= 0.0
		ENDDO
	    
		DO k = 9, 26
			DO i = ns,1,-1
				IF(i.eq.1) THEN
					wqualdv(i,k)=setldep(i,k)*wqual(i,k)*arfac*area(i)
					deptop(i,k)=depth(i)-setldep(i,k)
					IF(deptop(i,k).le.0.0) deptop(i,k)=0.0 
					depbotom(i,k)=0.0
				ELSE
					deptop(i,k)=depth(i)-setldep(i,k)
					depbotom(i,k)=depth(i-1)-setldep(i,k)		
					IF(setldep(i,k).lt.dep(i))THEN
						wqualdv(i,k) = 1.0*wqual(i,k)*setldep(i,k)*arfac*0.5*(area(i)+area(i-1)) !mg	
				
					ELSEIF(setldep(i,k).ge.dep(i))THEN
						wqualdv(i,k) = 1.0*wqual(i,k)*(dep(i)*arfac*0.5*(area(i)+area(i-1)))	!mg
					ENDIF
				ENDIF
		  		fluxin(i,k) = 0.0d0
			ENDDO	! i
		ENDDO	! k
!gbs Convert to concentration
		DO k = 9,26
			DO i = ns,1,-1
				DO j = ns,1,-1
!gbs case where all settling particles fall within top and bottom of layer 1:

					IF(j.eq.1)THEN

!gbs case where settling particles fall within layer 1 and extend further than the bottom:
						IF(deptop(i,k).le.depth(j).and.deptop(i,k).gt.0.0.and.depbotom(i,k).le.0.0)	&
     						fluxin(j,k) = fluxin(j,k)+ wqualdv(i,k)

!gbs case where settling particles extend above the top of layer 1 & extend further than the bottom:
						IF(deptop(i,k).gt.depth(j).and.depbotom(i,k).le.0.0)		&
     	    				fluxin(j,k) = fluxin(j,k)+ wqualdv(i,k)

!gbs case where all settling particles fall within layer 1:
						IF(deptop(i,k).le.depth(j).and.depbotom(i,k).gt.0.0)		&
     						fluxin(j,k) = fluxin(j,k)+wqualdv(i,k)

!gbs case where settling particles extend above layer 1 and fall within
!gbs layer 1
						IF(deptop(i,k).gt.depth(j).and.depbotom(i,k).lt.depth(j).and.depbotom(i,k).gt.0.0)	&
     	   				fluxin(j,k)=fluxin(j,k)+wqualdv(i,k)
					ELSE

!gbs case where all settling particles fall within the top and bottom of a layer
						IF(deptop(i,k).le.depth(j).and.depbotom(i,k).ge.depth(j-1)) THEN
      					fluxin(j,k) = fluxin(j,k)+wqualdv(i,k)

!gbs case where settling particles extend above the top of a layer and
!gbs below the bottom of the layer
						ELSEIF(deptop(i,k).ge.depth(j).and.depbotom(i,k).lt. depth(j).and.		&
                        depbotom(i,k).ge. depth(j-1)) THEN
      						fluxin(j,k) = fluxin(j,k)+wqualdv(i,k)
!gbs   case where the settling dep above the layer, but bottom is below the layer and
!gbs    bottom is above the below bottom layer
						ELSEIF (deptop(i,k).ge.depth(j).and.depbotom(i,k).lt. depth(j).and.		&
                                depbotom(i,k).le. depth(j-1)) THEN
								fluxin(j,k) = fluxin(j,k)+wqualdv(i,k)
						ENDIF
					ENDIF
				ENDDO ! j
			ENDDO	! i
		ENDDO ! k

      DO k = 9,22
			DO i = ns,1,-1
				IF(i.eq.1) THEN
					wqual(i,k) = wqual(i,k)+(fluxin(i,k)-wqualdv(i,k))/(1.0d0*dep(i)*area(i)*arfac)  !mg/m^3
! POP and PON of the botttom layer settle down permanently
					IF(k.eq.9.or.k.eq.11.or.k.eq.12.or.k.eq.13.or.k.eq.17.or.k.eq.18.or.k.eq.19) THEN
						wqual(i,k)=0.0d0
					ENDIF    		
				ELSE
					wqual(i,k) = wqual(i,k)+(fluxin(i,k)-wqualdv(i,k))/					&
                     (1.0d0*dep(i)*arfac*0.5*(area(i)+area(i-1)))  !mg/m^3
				ENDIF
				IF(wqual(i,k).le.0.0) wqual(i,k)=0.0001
			ENDDO
		ENDDO ! k
!  
!gbs boundary check...
!  
!		IF(wqual(i,9).lt.0.5)THEN
!			WRITE(*,*) 'Warning, BOD less than zero after Particles',jday,depth(i),wqual(i,9)
!gbs		wqual(i,9) = 0.5
!			PAUSE
!     ENDIF  
	
!		IF(wqual(i,13).lt.0.0001)THEN
!			WRITE(*,*)' Warning, after particles les than zero POP'			&
!     	      ,jday,depth(i),wqual(i,13)
!			wqual(i,13) = 0.0001
!			PAUSE
!		ENDIF
		
!		IF(wqual(i,19).lt.0.0001) THEN
!			WRITE(*,*)'Warning, after particles less than zero PON',jday,depth(i),wqual(i,19)
!			wqual(i,19) = 0.0001
!			PAUSE
!		ENDIF
		RETURN
		END SUBROUTINE NEUTRIENT_SETTLING
!****************************************************************************************
