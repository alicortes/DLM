!********************************************************************************
		SUBROUTINE oxygen(deltnho)
!*********************************************************************************
! this calculates oxygen concentration in the water column due to
! the sources of aeration and production, and sinks of nitrification,
! respiration, sediment oxygen demand and biochemical oxygen demand

! alterdo = correction to ensure that total oxygen is always > 0 (mg l**-1)
! are = change in area of successive layers (m**2)
! arfac = factor for multiplying areas and volumes to get m**2 or m**3
! K_SOD = coefficient for biochemical sediment oxygen uptake (g m**-2 d**-1)
! constnh = coefficient for conversion of nh3 to no3
! constr = coefficient for respiration (d**-1)
! deltnho = release of oxygen from nitrate utilisation (mg l**-1) [see subr. nutrit]
! deltabo = temporary variable for breakdown of individual algal groups
! dep = change in depth of successive layers (m**2)
! gromax = maximum potential growth rate of algae (d**-1)
! hscnit = half saturation constant for nitrification limitation
! by DO (mg l**-1)
! hscdo = half saturation constant for biological limitation by DO in
! the water (mg l**-1)
! KH_OSOD = half saturation constant for biological limitation by DO in
! the sediment (mg l**-1)
! nchl = number of chlorophyll groups being modelled
! obod = change in oxygen bue to biochemical oxygen demand (mg l**-1)
! DOsat = oxygen concentration at saturation
! delta = change in oxygen concentration over the time step
! DO_NITRI = change in oxygen due to nitrification (mg l**-1)
! ophoto = change in oxygen due to photosynthesis and respiration (mg l**-1)
! ustar = wind shear at the surface
! KL = oxygen transfer velocity at the surface
! segmin = minimum of limitation by light/phosphorus/nitrogen/silica
! sod = sediment oxygen demand (g m**2 d**-1), *depth
! thetat = temperature multiplier for phytoplankton growth
! thetabs = temperature multiplier for biochemical sediment oxygen demand
!  thetanh = temperature multiplier for conversion of nh3 to no3
!  ycabod = ratio of mg DO to mg chla utilised in decay of phytoplankton
!  C_CHLA = ratio of mg carbon to microg chla (see subr. algae)
!  y1chldo = ratio of mg DO to mg chla released in photosynthesis
!  y2chldo = ratio of mg DO to mg chla utilized in respiration
!  ynhdo = ratio of mg nh3 (1) consumed to mg DO (4) utilized in nitrification 
!          NH3+2O2=>HNO3+H2O 
!  zphoto = depth below the euphotic zone (m from bottom)
!sam wqual(i,4) = inflow conc., see var. 'wqins'
!sam wqual(i,5) = surface reaeration 2015/08 GBS: SOLVED IN OXYGEN MODULE
!sam wqual(i,6) = starting profile DO, model adjusts by mixing, not by rxns. or fluxes 2015/08 GBS: SOLVED IN OXYGEN MODULE AND DIFFUSE MODULE
!    DIFFUSE: !2015/08  C(i,10)=C(i,6) +C(i,7)+C(i,8)+C(i,9)
!     print*,'C(ns,6)=',C(ns,6)    !C(i,6)=wqual(i,4)---
!	   print*,'C(ns,7)=',C(ns,7)	 !C(i,7)=wqual(i,5)--- surface
!	   print*,'C(ns,8)=',C(ns,8)	 !C(i,8)=wqual(i,6)---
!	   print*,'C(ns,9)=',C(ns,9)	 !C(i,9)=wqual(i,7)---source and sink terms
!	   print*,'C(ns,10)=',C(ns,10)	 !C(i,10)=wqual(i,8)DO
!     C(i,11)=wqual(i,9)BOD
!sam wqual(i,7) = find below, for algae, SOD, BOD, nitrif/denitrif. 2015/08 GBS: SOLVED IN OXYGEN MODULE
!sam wqual(i,8) = DO concentration
!sam wqual(i,9) = BOD concentration
!sam wqual(i,16) = NH3 concentration
!sam wqual(i,28) = pH
!sam do_set = minimum DO setting for hypolimnetic oxygenation (mg L-1)
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE
		
		REAL*8 alterdo
		REAL*8 are(maxns)
		REAL*8 deltabo(maxns)
		REAL*8 deltnho(maxns)
		REAL*8 dep(maxns)
		REAL*8 obod(maxns)
		REAL*8 DOsat(maxns)
		REAL*8 sur_del_DO
		REAL*8 DO_NITRI(maxns)
		REAL*8 ophoto(maxns)
		REAL*8 KL
		REAL*8 q1,q2,q3,q4,q5,q6,q7,q8,q9,q10
		REAL*8 sod(maxns)
		REAL*8 sunny(maxns)
		REAL*8 ychlbod
		REAL*8 y1chldo
		REAL*8 y2chldo
		REAL*8 ynhdo
		REAL*8 zphoto
!sam   Sc is the Schmidt or Prandtl number for O2, rho_air is dry air density...
		REAL*8 Sc, rho_air, do_set

		INTEGER*4 i,j, time_day
		INTEGER*4 arfac,volfac
!IF		INTEGER*4 jday
!IF		INTEGER*4 nchl

		PARAMETER(q1 = 3.35e-2,q2 = 9.6e-4,q3 = 4.1e-5,q4 = 1.1e-1,q5 = 6.4022e-5)
		PARAMETER(q6 = 1.4101e+5,q7 = 5.24838e+3,q8 = 7.7117e0,q9 = 1.31403e0)
		PARAMETER(q10 = 4.593e+1)
		PARAMETER(arfac = 1.0e+6,volfac=1000.0d0)
		PARAMETER(Sc = 500, rho_air = 1.225, do_set = 1.0)
!gbs----------------------------------------------------------------------
		REAL*8 CD10,rhoa_rho0,z_0, tao, ustar, ustar_A,PN,KHnnt
		REAL(8), PARAMETER :: ustar_cr=0.11, tao_0=7.0d0
		REAL*8 NTM, KH_ONT, ACOR,FCD, AONT, Carbon_CHLA
		PARAMETER (rhoa_rho0 = 0.0012E0,  NTM= 0.075d0)
		PARAMETER (KHnnt=1.0d0,KH_ONT = 3.0d0, FCD=0.9d0)
      PARAMETER  (AONT=1.00d0)  !4.33
		REAL(8), PARAMETER :: pi=3.1415926535897932384626433832795d0
		REAL*8 denh2o(maxns),dyn_visco(maxns),kin_visco(maxns)
		REAL*8 first_term, second_term,CO2_CHLA 
!gbs-------------SOD------------------------------------------
      REAL*8 SOD_area (maxns), SOD_radius (maxns)
		REAL*8 SOD_radius_L (maxns),SOD_radius_U (maxns)
		REAL*8 slant_height(maxns)
! ----------------------------------------------------------------------
		REAL*8 atm_press(maxns), part_press_wat_vap(maxns), DO_theta(maxns), DOsat_eq(maxns),SOD_Area_450
! --------------------------------------------------------------------------
!  set the yields of phytoplankton carbon:DO, nh:DO and phytoplankton
!  breakdown (cxhyoz->co2+h2o)

		y1chldo = 32.0/12.0
		y2chldo = 32.0/12.0
		ynhdo   = 64.0/14.0
		ychlbod = 32.0/12.0
		CO2_CHLA  =40.0/100.0
        Carbon_CHLA  = 0.40d0 !Algal biomass (gC/m3)
		ACOR=2.67d0
!************* SOURCE 1:SURFACE REAERATION**************************
!  source 1: calculate the effect of surface aeration on oxygen.
!  use the method of patterson et al. (1985). the incorrect
!  coeffficients from the 1985 paper have been adjusted.
!------density------------------------------------------
      DO i=1,ns
			den(i)=densty(temp(i),sal(i))
      ENDDO

		DO i = 1,ns
			denh2o(i) = den(i)+1000.0d0	
!     Absolute viscosity as a function of tmperature
			dyn_visco(i)=(-3.86704098283854D-10*(temp(i)**5.0)+		    &
              1.28275099347771D-7*(temp(i)**4)-							&
              1.73356497821936D-5*(temp(i)**3)+							&
              1.27453807127464D-3*(temp(i)**2)-							&
              0.0587433984793861*temp(i)+								&
              1.785149374874680)*0.001d0
!-----------kinematic viscosity---------------------
			kin_visco(i)=dyn_visco(i)/denh2o(i)
		ENDDO

!  calculate the wind-air shear at the surface (from Ivey and
!   Patterson 1984)
		IF (U6X.le.4.0) THEN
			!CD10=1.124E-3
			CD10=1.13E-3
			ustar = U6X*dsqrt(CD10)			!ustar = q1*U6X		
		ELSEIF(U6X.gt.4.0) THEN							
			CD10=(9.6+0.41*U6X)*10.0E-4
			ustar = U6X*dsqrt(CD10)		    !ustar = U6X*(q2+q3*U6X)**0.5									
		ENDIF
	
		ustar_A = ustar*dsqrt(denh2o(ns)/rho_air) 
! 
!GBS calculate the oxygen transfer velocity, according to Ivey and patterson(1985) and O'Connor (1983)
! 
	
		IF(ustar.lt.ustar_cr) THEN
			KL = ((2.0E-9/kin_visco(ns))**0.66667)*dsqrt(rhoa_rho0)*		&
					(0.40**0.333333)*ustar/tao_0                             !Equation 23a
		ELSE
			tao = tao_0*(ustar/ustar_cr)*DEXP(-((ustar/ustar_cr)-1.0))     !Equation 27
			z_0=(1.0/0.0025+(8.0*ustar/kin_visco(ns))*						&
              dexp(-ustar/0.10))**(-1.0)											!Equation 15a
			first_term= ((2.0E-9/kin_visco(ns))**0.66667)*					&	!Equation 26a
               dsqrt(rhoa_rho0)*(0.40**0.333333)*ustar/tao 	   
			second_term=dsqrt((2.0E-9*ustar/(0.40*z_0))*rhoa_rho0*		&	!Equation 26a
                    (1.0E-7/kin_visco(ns)))
         KL=(1.0/first_term+1.0/second_term)**(-1.0)       !Equation 26a                                                             
		ENDIF
	
!  calculate saturated oxygen concentration for surface temperature,
!  based on the formula of Ivey and patterson(1985) and Mortimer (1981)
!	DOsat = dexp(q8-q9*dlog(temp(ns)+q10))        
		DOsat(ns) = 14.5532-0.38217*temp(ns)+0.0054258*temp(ns)*temp(ns)-	&
             (sal(ns)/1.80655)*(1.665d-4-5.866d-6*temp(ns)+9.796d-8*temp(ns)*temp(ns))
!------effect of altitude----------
      atm_press(ns)=dexp(5.25*dlog(1.0-(CRL+BASE)/(44.3*1000.0)))
		part_press_wat_vap(ns)=dexp(11.8571-3840.70/(temp(ns)+273.15)-		&
           216961.0/(temp(ns)+273.15)**2.0)
      DO_theta(ns)=0.000975-1.426d-5*temp(ns)+6.436d-8*temp(ns)**2.0

		DOsat_eq(ns)=DOsat(ns)*atm_press(ns)*(1.0-part_press_wat_vap(ns)/atm_press(ns))*	&
             (1.0-DO_theta(ns)*atm_press(ns))/((1.0-part_press_wat_vap(ns))*(1.0-DO_theta(ns)))
!     print*, DOsat,DOsat_eq,atm_press,(CRL+BASE),temp(ns)
!		pause
      
!  calculate change in oxygen for the surface layer over the time step

		IF(ns.eq.1)THEN
			sur_del_DO = KL*dfloat(nosecs)*(DOsat_eq(ns)-wqual(ns,8))/depth(ns)
		ELSE
			sur_del_DO = KL*dfloat(nosecs)*(DOsat_eq(ns)-wqual(ns,8))/(depth(ns)-depth(ns-1))
		ENDIF

!	print*,sur_del_DO,wqual(ns,8)
!  calculate the oxygen going into the surface layer through aeration	  
		wqual(ns,8) = wqual(ns,8)+sur_del_DO 
!2015/08	wqual(ns,5) = wqual(ns,5)+sur_del_DO
!************* SOURCE 2: ALGAE PHOTOSYNTHESIS Cerco report(2004) **************************
!  determine change in oxygen due to photosynthesis/respiration (ophoto),
!  algal death (deltabo) biochemical oxygen demand (obod), nitrification
!  (DO_NITRI) and sediment oxygen demand (sod).

		DO i = 1,ns
			ophoto(i)  = 0.0
			deltabo(i) = 0.0
			DO_NITRI(i)= 0.0
			obod(i)    = 0.0
			sod(i)     = 0.0
		ENDDO
!
		DO i = 1,ns
     	   DO j = 1,nchl
				IF((wqual(i,21) + wqual(i,22)).le. 0.0) THEN
               PN = 0.0 
				ELSE
					PN = wqual(i,21)*wqual(i,22)/							&	!=(NH4*NO3)/((hs+NH4)*(hs+NO3))
     				((hscn(j)+wqual(i,22))*(hscn(j)+wqual(i,21)))+	&	!+NH4*Hs/((NH4+NO3)*(hs+NO3))
     				wqual(i,22)*hscn(j) / ((wqual(i,21) +				&
     				wqual(i,22))*(hscn(j)+wqual(i,21)))
				ENDIF
				IF (PN.GE.1.00) THEN
					PN=1.00
				ELSEIF (PN.LE.0.00) THEN
					PN=0.00
                ENDIF
               
! respiration reduces DO...
!				ophoto(i)= ophoto(i)+(Carbon_CHLA*wqual(i,j))*ACOR*(1.3-0.3*PN)*growth(i,j)		&
!                  -(1.0-FCD)*constr(j)*(dfloat(nosecs)/86400.0)*wqual(i,j)*dexp(0.08*(temp(i)-30.0))

	         if(wqual (i,j).ge.0.001) then
				ophoto(i)= ophoto(i)+(Carbon_CHLA*60.0)*ACOR*(1.3-0.3*PN)*growth(i,j)		&
                   -(1.0-FCD)*constr(j)*(dfloat(nosecs)/86400.0)*wqual(i,j)*dexp(0.08*(temp(i)-30.0))
             elseif (wqual (i,j).lt.0.001) then      
               ophoto(i)=0.0
           endif 
           
  
!		write(*,fmt='(2f10.3)') growth(i,j),death (i,j)
! mortality adds BOD...
!		deltabo(i) = deltabo(i)+(constm(j)*dexp(0.04*(temp(i)-20.0)))	&
!     	*wqual(i,j)*CO2_CHLA*ychlbod*(nosecs/86400.0) 
				deltabo(i) = deltabo(i)+(constm(j)*dexp(0.05*(temp(i)-20.0)))		&
     						*wqual(i,j)*ychlbod*(dfloat(nosecs)/86400.0)   	
			ENDDO	! j
!sam  but IF there's nothing to respire...
			IF(ophoto(i).lt.0.0.and. abs(ophoto(i)).gt.wqual(i,8)) ophoto(i)=0.0            
!------------ 3. SINK: NITRIFICATION-----------------------------------
!  note the division by 1000 in the following equation since we are
!  dealing with ug N/l -> mg d.o./l

			DO_NITRI(i)=AONT*thetanh**(temp(i)-20.0)*						&
     			(wqual(i,8)/(KH_ONT+wqual(i,8)))*							&
     			(wqual(i,22)/(KHnnt+wqual(i,22)))*(dfloat(nosecs)/86400.0)*(NTM*dfloat(nosecs)/86400.0)

!------------ 4. SOURCE: BOD (BIO-LOGICAL OXYGEN DEMAND)------------------
			obod(i) = constbo*thetabo**(temp(i)-20.0)*wqual(i,11)*		&
     	        (wqual(i,8)/(hscdo+wqual(i,8)))*(dfloat(nosecs)/86400.0)

!------------ 5. SINK: DOC(DISSOLVED ORGANIC CARBON)------------------
!			DO_DOC (i)=wqual(i,9)*ACOR*K_DOC*								&
!               (wqual(i,8)/(KH_ODOC+wqual(i,8)))*(nosecs/86400.0)

		ENDDO		! i
!------------ 6. SINK: COD(CHEMICAL OXYGEN DEMAND)------------------
!			ocod(i) = constbo*thetabo**(temp(i)-20.0)*wqual(i,**)*		&
!     	        (wqual(i,8)/(KHocod+wqual(i,8)))*(dfloat(nosecs)/86400.0)
!------------ 7. SINK: SOD(SEDIMENT OXYGEN DEMAND)------------------
!  determine change in oxygen due to sediment oxygen demand (SOD). include
!  biological and chemical sediment oxygen demand (walker and snodgrass
!  1986)
		DO i = ns,2,-1		
			dep(i) = depth(i)-depth(i-1)
			SOD_radius (i)= DSQRT((0.5*(area(i)+area(i-1))*arfac)/PI)
			SOD_radius_L(i)=DSQRT(area(i-1)*arfac/PI)
			SOD_radius_U(i)=DSQRT(area(i)*arfac/PI)
			slant_height(i)=DSQRT(dep(i)*dep(i)+(SOD_radius_U(i)-SOD_radius_L(i))**2)
			SOD_area   (i)= 2.0D0*PI*SOD_radius (i)*slant_height(i)			
		ENDDO
		are(1) = area(1)*arfac
		dep(1) = depth(1)
		SOD_area(1)=are(1)
!      SOD_Area_450=0.0
!      DO i=1,ns
!         if ((depth(ns)-depth(i)).ge.450.0) then
!            SOD_Area_450=SOD_Area_450+SOD_area   (i)
!         endif
!      enddo
!   print*,SOD_Area_450
!   pause
!  sediment oxygen demand. any cod and bod in the sediments is likely to
!  be negated by photosynthesis on the sediments of the layers within
!  the euphotic depth. -> only during daylight! (SAM)
!  below this we are assuming a constant SOD caused by cod and bod.
!  we assume the chemical sediment oxygen demand is constant between
!  lakes - K_SOD, a user-defined input, will have the major effect on sod.

!sam no compensating photosynthesis in sediments at night...

		IF (iclock.lt.sunrise .or. iclock.gt.sunset) THEN
			zphoto = depth(ns)+.01
			GOTO 235
		ENDIF

		sunny(ns) = 1.0

		IF(ns.eq.1)THEN
			zphoto = depth(ns)-4.605/et1(i)
			GOTO 235
		ENDIF

		DO i = ns-1,1,-1
			sunny(i) = sunny(i+1)*DEXP((depth(i+1)-depth(i))*(-et1(i+1)))

			IF(sunny(i).le.0.01) THEN
				zphoto = depth(i+1)-DLOG(0.01/sunny(i+1))/(-et1(i+1))
				GOTO 235
			ENDIF

			IF(i.eq.1.and.sunny(i).gt.0.01) THEN
				zphoto = depth(i)-4.605/et1(i)
			ENDIF
		ENDDO

235	CONTINUE
    	DO i = 1,ns
			sod(i)=0.0
			IF(i.eq.1)THEN
				IF(zphoto.lt.0.0)THEN
					sod(i) = 0.0
				ELSE IF(zphoto.gt.0.0.and.zphoto.lt.depth(i))THEN
					sod(i) = K_SOD*(wqual(i,8)/(KH_OSOD+wqual(i,8)))*			&
		  			(thetabs**(temp(i)-25.0))*(DFLOAT(nosecs)/86400.0)*		&
		  			SOD_area(i)*(zphoto/depth(i))                       ! (g/m2-d)*m2*d=g
				ELSE IF(zphoto.gt.depth(i))THEN
				sod(i) = K_SOD*(wqual(i,8)/(KH_OSOD+wqual(i,8)))*				&
		  			(thetabs**(temp(i)-25.0))*(DFLOAT(nosecs)/86400.0)*SOD_area(i)  ! (g/m2-d)*m2*d=g
				ENDIF
				sod(i)=sod(i)/(vol(i)*volfac)    !Covert into g/m3
			ELSE
				IF(depth(i-1).gt.zphoto)THEN
					sod(i) = 0.0
				ELSEIF(depth(i-1).lt.zphoto.and.depth(i).gt.zphoto)THEN
					sod(i) = K_SOD*(wqual(i,8)/(KH_OSOD+wqual(i,8)))*			&
		  			(thetabs**(temp(i)-25.0))*(DFLOAT(nosecs)/86400.0)*		&
		  			SOD_area(i)*((zphoto-depth(i-1))/(depth(i)-depth(i-1)))      ! (g/m2-d)*m2*d=g
				ELSE IF(depth(i).lt.zphoto)THEN
					sod(i) = K_SOD*(wqual(i,8)/(KH_OSOD+wqual(i,8)))*			&
		  			(thetabs**(temp(i)-25.0))*(DFLOAT(nosecs)/86400.0)*SOD_area(i)  !(g/m2-d)*m2*d=g
				ENDIF
					sod(i)=sod(i)/(vol(i)*volfac)   ! Covert into g/m3
			ENDIF
		ENDDO
    
!  sum the components for oxygen and bod calculations. include nitrate
!  reduction by phytoplankton as a source

		DO i = 1,ns-1
!2015/08	wqual(i,7) = wqual(i,7)+ophoto(i)+deltnho(i) - DO_NITRI(i)-obod(i)-sod(i)
			wqual(i,8) = wqual(i,8)+ophoto(i)+deltnho(i) - DO_NITRI(i)-obod(i)-sod(i)-deltabo(i)  
!			wqual(i,9) = wqual(i,9)-obod(i)+deltabo(i)

!  make sure that the total oxygen is not below zero - correct biological
!  processes since these will not CONTINUE as fast at low DO
!			IF(wqual(i,7).lt.0.0) wqual(i,7) = 0.0	  
			IF(wqual(i,8).lt.0.0) wqual(i,8) = 0.0
!			IF(wqual(i,8).gt.14.0) wqual(i,8) = 14.0	  
        ENDDO
       
        DO i=1,ns   
! 2018/05 ST - is eraseing salinity balance		
        !	sal(i)=0.0
		    DOsat(i) = 14.5532-0.38217*temp(i)+0.0054258*temp(i)*temp(i)-	&
                    (sal(i)/1.80655)*(1.665d-4-5.866d-6*temp(i)+9.796d-8*temp(i)*temp(i))
!------effect of altitude----------
            atm_press(i)=dexp(5.25*dlog(1.0-(CRL+BASE)/(44.3*1000.0)))
		    part_press_wat_vap(i)=dexp(11.8571-3840.70/(temp(i)+273.15)-216961.0/(temp(i)+273.15)**2.0)
            DO_theta(i)=0.000975-1.426d-5*temp(i)+6.436d-8*temp(i)**2.0

		    DOsat_eq(i)=DOsat(i)*atm_press(i)*(1.0-part_press_wat_vap(i)/atm_press(i))*	&
             (1.0-DO_theta(i)*atm_press(i))/((1.0-part_press_wat_vap(i))*(1.0-DO_theta(i)))            
			IF(wqual(i,8).gt.(DOsat_eq(i)+5.0)) THEN !Super saturation 3-4 mg/l
				wqual(i,8) =(DOsat_eq(i)+5.0)
	        ELSE
                wqual(i,8) =wqual(i,8)
            ENDIF            
			IF(wqual(i,8).lt.0.0d0) wqual(i,8) = 0.001d0
        ENDDO
!        IF (wqual(ns,8).lt.DOsat(ns)) wqual(ns,8)=DOsat_eq(ns)+ophoto(ns)
      DO i = 1,ns
!			IF(iclock.eq.43200) THEN
				IF(i.eq.ns) THEN
					write(105,1050) jday,iclock, depth(i),wqual(i,8),DOsat(i),temp(i),sur_del_DO,		&
                                    ophoto(i),deltnho(i),-DO_NITRI(i),-obod(i),-sod(i)
1050				FORMAT(2i10,3x, 4f10.2,6f10.4)     
				ENDIF
!			ENDIF
		ENDDO
888	CONTINUE 
		RETURN 
		END SUBROUTINE oxygen
