!************************************************************************
      SUBROUTINE ALGAE_SPECIES(par,limname,SDname,latit)
!***********************************************************************  
!  Subroutine to calculate the change in chlorophyll concentration
!--------------------------------------------------------------------------
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE		

!    **** defining variables ****
!     alg_min = minimum algae conc. [ug/L]
!     constr = coefficient for respiration
!     constm = coefficient for temperature dependent mortality
!     conquan = conversion factor for light (mol photons m**-2 watts**-1)
!     conquan = conversion factor for light (mol_quanta day-1 watts-1)
!     e = base of natural log
!     et1 = light extinction coefficient
!     etca = specific attenuation (m**2 mg**-1) by Chl a
!     etwat = attenuation coefficient of water (m**-1), inorganic particles etc.
!     faco = factor (0 or 1) for reducing phytoplankton growth to zero when
!     jday = julian date for output file
!     light = light for saturation (W m-2 PAR)
!     lt_threshold = minimum light for light response of algae, now in units [W/m2 PAR]
!     oxygen concentrations are low
!     gromax = maximum phytoplankton growth rate
!     hscsi = half saturation constant for silica uptake
!     nchl = number of algal groups in model
!     paro = PAR light level in a layer [W/m2]
!     par = PAR light level at the surface [W/m2]
!     pinmin = minimum internal phosphorus to chlorophyll mass ratio
!     ninmin = minimum internal nitrogen to chlorophyll mass ratio
!     pinmax = maximum internal phosphorus to chlorophyll mass ratio (hardwired)
!     ninmax = maximum internal nitrogen to chlorophyll mass ratio (hardwired)
!     satlit = light intensity at saturation, f(T), [W/m2 PAR]
!     seglit = fractional reduction in growth due to sub-optimal light
!     segp = fractional reduction in growth due to sub-optimal phosphorus
!     segn = fractional reduction in growth due to sub-optimal nitrogen
!     segsi = fractional reduction in growth due to sub-optimal silica (diatoms only)
!     segmin = fractional total reduction in maximum phytoplankton growth
!     thetat = temperature multiplier for phytoplankton growth
!     ychlc = ratio of mg carbon to ug Chl a
!     yquanta = quantum yield (mg carbon per mol light quanta absorbed)

      CHARACTER*20 limname(1)
	   CHARACTER*20 SDname
	   INTEGER*4 i
	   INTEGER*4 j
      INTEGER*4 k,ku,kua
!IF	   INTEGER*4 nchl
!IF	   INTEGER*4 jday
	   REAL*8 conquan
	   REAL*8 e
	   REAL*8 hscnhno
      REAL*8 P_Uptake,N_Upatke,P_RandM,N_RandM,P_Grazing,N_Grazing

!	   REAL*8 etca(3)
	   REAL*8 limiting
	   REAL*8 Minimum_Nut,Minimum_Growth 	
	   REAL*8 par
	   REAL*8 satlit(maxns,maxchl),halfsatn(maxns,maxchl),halfsatp(maxns,maxchl)
	   REAL*8 seglit(maxns,maxchl)
	   REAL*8 segp(maxns,maxchl),segn(maxns,maxchl),segsi(maxns,maxchl)
	   REAL*8 yquanta  !, alg_min
	   REAL*8 lt_threshold	!!sam for light penetration output 
!	1 W/m2 SW = 2.07 uE/m2/s
      REAL*8 zoo
 
	   REAL*8 Chl_Top
	   REAL*8 Growth_Top
	   REAL*8 Graz_Top
	   REAL*8 Death_Top

	   REAL*8 Growth_Top_2
	   REAL*8 Graz_Top_2
	   REAL*8 Death_Top_2, ref_temp(maxns), ref_temp1,threshold_temp
	   REAL*8 dep_ave_temp, sum_depth, sum_depth_temp 

       REAL*8 Vol_Mix_T,Chla_Mix,Chla_Mix_Media
	   REAL*8 latit,kmax,kmax1,kmax2, temp_opt,growth_rate(maxns,maxchl)
	   REAL*8 kk1, kk2,alphaa,betaa,maxgr(maxchl),Ir(maxns,maxchl)

!sam	PARAMETER(conquan = 720.0,e = 2.71828,yquanta = 0.96418)
	   PARAMETER(conquan = 0.4, e = 2.71828, yquanta = 0.1)  ! per Brad
	   PARAMETER(lt_threshold = 15.0/2.07*.45)	!!sam
!	PARAMETER(alg_min = 0.05)
              
 	   zoo  = 0.001d0*(dfloat(nosecs)/86400.0d0) !(#/L)
      alg_min = 0.001  !   0.15  0.1

	   IF(iclock.eq.0.0) THEN
	      ref_temp1=temp(ns)
	   ENDIF
	   ref_temp(ns)=ref_temp1
	   threshold_temp=5.6
	
	   dep_ave_temp = 0.0
	   sum_depth =0.0
	   sum_depth_temp =0.0
	   DO i =1, ns
	      IF(depth(i).ge.485.0)THEN
	         sum_depth = sum_depth + depth(i)-depth(i-1)
		      sum_depth_temp = sum_depth_temp +temp(i)*(depth(i)-depth(i-1))
	      ENDIF
	   ENDDO
	   dep_ave_temp =sum_depth_temp/sum_depth
!	   print*,dep_ave_temp 
!	   print*,iclock,ref_temp,temp(ns)
!  ----------------------------------------------------
!     part 1: determining the light limitation on growth.
!     change the value of photosynthetically active radiation for
!     the surface to the surface value for an array of par

	   paro(ns) = par	! set in heatr
!     determine the attenuation coefficient for par from chla. note that
!     etwat could be replaced by a coefficient related to particle
!     concentrations at some stage
	   DO i = 1,ns
	      et1(i) = 0.0
	      DO k = 29, 29+nchl_species-1
		      et1(i) = et1(i)+etca(j)*wqual(i,j)
	      ENDDO ! j
	      DO j = 1,7
		      et1(i) = et1(i)+etpart(j)*cf(i,j)
	      ENDDO ! j
	      et1(i) = et1(i)+etwat
	   ENDDO   ! i
!  
!quim   Ted Swift clarity model for Tahoe
!  ******************LIGHT_TAHOE**********************************************
      CALL Light_Tahoe(latit)
!  ****************************************************************************
!     calculate PAR at the surface of each layer
!  
!	   et1(i)=0.50
	   IF(ns.eq.1)GOTO 15
	   DO i = ns-1,1,-1
	      IF (paro(i+1).le.1d-20) THEN
		      paro(i) = 0.0
	      ELSE
		      paro(i) = paro(i+1)*(dexp(-et1(i)*(depth(i+1)-depth(i))))
	      ENDIF
	   ENDDO
15    CONTINUE

!   calculate mean par within each layer. note that the depth of
!   the surface layer equates to the depth of the mixed layer here.
!   also determine the saturating light intensity as a fraction
!   of the mean daily light intensity. set upper and lower bounds
!   for the saturating light intensity. calculate the ratio of mg
!   carbon to ug chla (wasp 1988) and set upper and lower limits.

	   DO i = ns,1,-1
	      IF(i.eq.1)THEN
		      paro(i) =(paro(i)-paro(i)*(dexp(-et1(i)*depth(i))))/(et1(i)*depth(i))
	      ELSE
		      paro(i) =(paro(i)-paro(i-1))/(et1(i)*(depth(i)-depth(i-1)))
	      ENDIF

!sam for outputting light penetration depth, compare as W/m2 PAR...
		   IF ( paro(i).gt. lt_threshold.and. depth(i).lt.z_light ) z_light=depth(i)	    
!sam used in zoo & oxygen. Since constant, just set values for now...
	      Do j = 29, 29+nchl_species-1 	
!	         light(j)=40.0		!40 35,12    
!		      thetat(j)=1.14		!1.13 
	   	   satlit(i,j) =light(j)*thetat(j)**(temp(i)-20.0)
!	         satlit(i,j) =light(j)*thetat(j)**(temp(i)-20.0)
	      ENDDO	! j
	   ENDDO	! i
!Find the growth reduction dueto insufficient light, i.e. the P-I curve, 
!using the Steele equation modified Jassby and Platt...
 
 	   DO i = ns,1,-1
	      DO j = 29, 29+nchl_species-1
		      seglit(i,j)=1*((paro(i)/satlit(i,j))**1.0)*dexp(1.0-(paro(i)/satlit(i,j))**1.0)	
	      ENDDO	! j
	   ENDDO	! i
!  ********************************************************************************
!   part 2: calculate the fractional reduction in the maximal growth
!   due to nutrient limitation (phosphorus, and nitrogen). THEN
!   determine  the minimum of limitation by light/p/n/si. THEN calculate
!   chlorophll concentration considering only growth

!  
! Consider growth term. No growth IF there is low oxygen.
!  
	   DO i = 1,ns
		   IF (wqual(i,8) .lt. 0.5) THEN
			   faco(i) = 0.0
		   ELSE
			   faco(i) = 1.0
		   ENDIF
	   ENDDO
!  
!quim   Set nutrient limitation factor. Michaels-Menten kinetics
!quim   Liebig's law of the minimun
!quim   Avoid math. errors
!  
	 		 
	   DO i = 1,ns	   	  
	      DO j = 29, 29+nchl_species-1	
		      halfsatn(i,j) =halfn(j)*1.13**(temp(i)-20.0)
		      halfsatn(i,j) =halfn(j)*1.13**(temp(i)-20.0)
	         segp(i,j) =wqual(i,10)/(halfsatp(i,j)+wqual(i,10)) 	
	         segn(i,j) =(wqual(i,15)+wqual(i,16))/(halfsatn(i,j)+(wqual(i,15)+wqual(i,16)))
!	         print*,segn(i,j), wqual(i,15),wqual(i,16)
	         segmin(i,j) = dmin1(seglit(i,j),segp(i,j),segn(i,j))
!		      segmin(i,j) = dmin1(segp(i,j),segn(i,j))
	         IF(segmin(i,j).gt.1.0.or.segmin(i,j).lt.0.0) THEN
	            WRITE(*,*)'Warning: ALGAE, minimum value out of range)'
	            WRITE(*,*)'segmin,i,j',segmin,i,j
	            stop
	         ENDIF
! Minimum threshold value for orthophsophate 0.5 microgram/L
	         IF (wqual(i,10).LT.0.9) THEN  !0.5
	            segmin(i,j) = 0.0
		         GOTO 4098
	         ENDIF		
		      IF (wqual(i,15).LT.1.5) THEN  
	            segmin(i,j) = 0.0
		         GOTO 4098
	         ENDIF
		      IF(iclock+nosecs.eq.itmpr.or.(itimes.eq.1440.and.itmpr.eq.86400)) THEN
	            open(202,file=limname(1),status='unknown',access='append')
		         IF (segmin(i,j).eq.  segp(i,j)) THEN
	               limiting = 0.0
		            GOTO 4097
		         ENDIF
		         IF (segmin(i,j).eq.  segn(i,j)) THEN
			         limiting = 1.0
			         GOTO 4097
		         ENDIF
		         IF (segmin(i,j).eq.seglit(i,j)) THEN
	               limiting = 2.0
			         GOTO 4097
		         ENDIF
 4097          CONTINUE
!		         IF(depth(i).gt.300) THEN
!			         WRITE(202,1034) jday, depth(i), limiting
!		         ENDIF
!1034		      FORMAT(i9,x,f10.1,x,f10.1) !,x,3(f10.5,x))
1035		      FORMAT(i10,x,f10.1,x, 3(f10.5,x))
!					CLOSE(202)
	         ENDIF
 4098	      CONTINUE
	         IF (segmin(i,j) .le. 0.0) THEN
               segmin(i,j) = 0.0
	         ENDIF
!  
!gbs	grow of the algae  
!   
!gbs******MINLAKE Algorithm for temperature multiplier for Algae growth 
!            temp_opt=4.0d0   !	   5.6
            temp_opt=20.0d0   !	   5.6 
!           IF(temp(i).le.temp_opt) THEN	   !!  temp(i)
!	            kmax =0.1*dexp(-2.3*((temp(i)-temp_opt)/(temp_opt-0.0))**2.0)	 !0.33 0.4 for two min	temp(i)
!           ELSEIF(temp(i).gt.temp_opt) THEN
!	            kmax =0.1*dexp(-2.3*((temp(i)-temp_opt)/(30.0-temp_opt))**2.0	 !0.33   0.4 for two min	temp(i)
!          	ENDIF
!************************ 
            IF(temp(i).le.temp_opt) THEN	   !!  temp(i)
	            kmax =1.0*dexp(-0.01d0*((temp(i)-temp_opt)**2.0))    !0.005  0.002
	         ELSEIF(temp(i).gt.temp_opt) THEN
	            kmax =1.0*dexp(-0.01d0*((temp_opt-temp(i))**2.0))	 !0.005  0.002		     
       	   ENDIF
!-------------Peters and Eilers, 1978 Hydrobiological Bulletin 12: 127-134----	
!	         maxgr=0.23d0	!1.5
	         maxgr(j) = gromax (j)
	         betaa = et1(i)  !
	         Ir(i,j)=paro(i)/satlit(i,j)
!	         IF(Ir(i,j).gt.0.0d0.and.paro(i).le.satlit(i,j)) THEN
!	            growth_rate(i,j)=maxgr
!	         ELSE
	         growth_rate(i,j)=maxgr(j)*2.0d0*(1.0d0+betaa)*Ir(i,j)/      &
                             (Ir(i,j)*Ir(i,j)+2.0d0*betaa*Ir(i,j)+1.0d0)
!	         ENDIF
!-------------Steel (1962), Limnology and Oceanography----
!	         maxgr=1.00
!	         Ir(i,j)=paro(i)/satlit(i,j)
! 	         growth_rate(i,j)=maxgr*(1.0d0/satlit(i,j))*dexp(1-Ir(i,j))	
!**************************************************************
!	         growth(i,j) = (wqual(i,j)*growth_rate(i,j)*kmax*seglit(i,j)
	         growth(i,j) = (wqual(i,j)*growth_rate(i,j)*kmax *segmin(i,j)*faco(i)*  &
	                        (dfloat(nosecs)/86400.0))	
!	         growth(i,j) = (wqual(i,j)*gromax(j)*thetat(j)**                         &
!      	                 (temp(i)-20.0)*segmin(i,j)*faco(i)*(nosecs/86400.0))
!  
!    To avoid less than zero Nutrients!  
!    Preferential ammonium uptake factor..(pg 270 Bowie et al.1985).
!  
		      IF((wqual(i,15) + wqual(i,16)).le. 0.0) THEN
               hscnhno = 0.0 
	         ELSE
		         hscnhno = wqual(i,16)*wqual(i,15)/((hscn(j)+wqual(i,16))*(hscn(j)+wqual(i,15)))+  &
     	                   wqual(i,16)*hscn(j)/((wqual(i,15)+wqual(i,16))*(hscn(j)+wqual(i,15)))
            ENDIF
	         IF(hscnhno.gt.1.or.hscnhno.lt.0) THEN
		         WRITE(*,*)'WARNING. In Routine ALGAE'
		         WRITE(*,*)'ammonium uptake factor (hscnhno) out of range'
		         WRITE(*,*)'hscnhno,i,j',hscnhno,i,j
		         stop
	         ENDIF
!  
! Estimate the Minimum alloable value
!  
            Minimum_Nut = dmin1(wqual(i,10)-ypchla(j)*growth(i,j),               &
     	                          wqual(i,15)-ynchla(j)*growth(i,j)*(1-hscnhno),   &
                                wqual(i,16)-ynchla(j)*growth(i,j)*hscnhno)

            IF(Minimum_Nut .lt. 0.0)THEN

               growth(i,j) = 0.9*(dmin1(wqual(i,10)/ypchla(j),                   &
      	                               wqual(i,15)/(ynchla(j)*(1-hscnhno)),     &
                                        wqual(i,16)/(ynchla(j)*hscnhno)))
	         ENDIF
	      
	         Minimum_Nut = dmin1(wqual(i,10)-ypchla(j)*growth(i,j),               &
     	                          wqual(i,15)-ynchla(j)*growth(i,j)*(1-hscnhno),   &
                                wqual(i,16)-ynchla(j)*growth(i,j)*hscnhno)
	
	         IF(Minimum_Nut .lt. 0.0)THEN
	            WRITE(*,*) 'Error in Algae.Growth too high will give <0 Nut values'	
               WRITE(*,*) Minimum_Nut
               WRITE(*,*) wqual(i,10)-ypchla(j)*growth(i,j)              !THP
               WRITE(*,*) wqual(i,15)-ynchla(j)*growth(i,j)*(1-hscnhno)  !NO3
               WRITE(*,*) wqual(i,16)-ynchla(j)*growth(i,j)*hscnhno      !NH4
	            STOP
	         ENDIF

	         wqual(i,j) = wqual(i,j) + growth(i,j)
!           IF (depth(i).gt.450.0.and.depth(i).lt.495.0) THEN
!gbs		      IF (jday.gt.1.and.jday.lt.65.0) THEN
!gbs	            WRITE(12, 111) jday, depth(i),segmin(i,1), paro(i),seglit(i,1),   &
!gbs              segp(i,1),segn(i,1), kmax, growth(i,j),growth_rate(i,j,wqual(i,j),death(i,j)             
!gbs           ENDIF	
	      ENDDO	! j chl
	   ENDDO	! i layers
  111 FORMAT(i5, 9f16.10)
!	WRITE(11111,3)jday,paro(ns-1),temp(ns-1),seglit(ns-1,1),
!     &      segp(ns-1,1),segn(ns-1,1),wqual(ns-1,1)
    
  3   FORMAT(i4, 8f10.5)

!  *********************************************************************	
!     part 3: work out the factors which limit phytoplankton growth
!     respiration, temperature dependent mortality, 
!     and predation (settling is handled by subr. particles).!  
!     account for temperature-dependent respiration and mortality
!     when considering the loss of algae
!  

	   DO i = 1,ns
	      DO j = 29, 29+nchl_species-1
!	         temp_opt=5.6d0   !	   5.6  5.2 	    
!		      constr(j) = 0.01*wqual(i,j)*dexp(2.3*((temp(i)-temp_opt)/	 &
!                      (0.0-temp_opt)))*(dfloat(nosecs)/86400.0)!                  
!     		constm(j) = 0.003*(1.12**(temp(i)-20.0))*(dfloat(nosecs)/86400.0)*wqual(i,j)
!     		death(i,j)= constr(j)+constm(j)	
!	         constr(j) =0.007	 !0.007
!	         constm(j) =0.003	 !0.003 

     	      death(i,j) = (constr(j)+constm(j))*dexp(0.08d0*(temp(i)-20.0))*		&   !0.069
                         (dfloat(nosecs)/86400.0)*wqual(i,j)
!	         death(i,j) = (constr(j)+constm(j))*                               &
!     	                1.08**(temp(i)-20.0)*(dfloat(nosecs)/86400.0)*wqual(i,j)	!0.07
!	         constr(j) =0.25	
!	         constm(j) =0.0088	  

!		      death(i,j) = constr(j)*growth(i,j)+ constm(j)*dexp(0.07*(temp(i)-20.0))*wqual(i,j)*  &
!                       (dfloat(nosecs)/86400.0) 
		
!		      IF(death(i,j).lt.0.0) death(i,j) = 0.0
            IF((wqual(i,j)-death(i,j)).lt.alg_min) THEN
               death(i,j) = 0.9*(wqual(i,j)-alg_min)
		         IF(death(i,j).lt.0.0) death(i,j) = 0.0
            ENDIF
!  
!           Estimate the Minimum alloable value
!  
            Minimum_Nut = dmin1(wqual(i,j)-death(i,j),                                   &
     	                          wqual(i,11)+ ypchla(j)*growth(i,j)-ypchla(j)*death(i,j), &
     	                          wqual(i,17)+ ynchla(j)*growth(i,j)-ynchla(j)*death(i,j))

            IF(Minimum_Nut .lt. 0.0)THEN  
               death(i,j) = 0.9*(dmin1(wqual(i,j),                                        &
     	                                (wqual(i,11)+ypchla(j)*growth(i,j))/ypchla(j),      &
      	                             (wqual(i,17)+ynchla(j)*growth(i,j))/ynchla(j)))	 
	            IF(death(i,j).lt.0.0) death(i,j) = 0.0	 
	         ENDIF

            Minimum_Nut = dmin1(wqual(i,j)-death(i,j),                                    &
     	                          wqual(i,11)+ ypchla(j)*growth(i,j)-ypchla(j)*death(i,j),  &
     	                          wqual(i,17)+ ynchla(j)*growth(i,j)-ynchla(j)*death(i,j))

	         IF(Minimum_Nut .lt. 0.0)THEN   
	            WRITE(*,*) 'Error in Algae.Death too high will give <0 PhytoN '	
               WRITE(*,*) 'Minimum_Nut',Minimum_Nut
               WRITE(*,*) 'chla-death',wqual(i,j)-death(i,j)     !Chla
               WRITE(*,*) 'PhytoP-Death',wqual(i,11)-(ypchla(j)*growth(i,j)-ypchla(j)*death(i,j))  !PhytoP
               WRITE(*,*) 'PhtyoN-Death',wqual(i,17)-(ynchla(j)*growth(i,j)-ynchla(j)*death(i,j))  !PhYtoN
	            PAUSE
	         ENDIF

!	         IF((wqual(i,j)-death(i,j)).lt.alg_min) THEN
!              WRITE(*,*)'Error ALGAE, after DEATH ',i,j,wqual(i,j),death(i,j)
!              pause
!	         ELSE
	            wqual(i,j) = wqual(i,j) - death(i,j)
!	         ENDIF
	 
	      ENDDO ! j chl
	   ENDDO    ! i layers
!  
!gbs account for grazing
!  
	   DO i = 1,ns
	      DO j = 29, 29+nchl_species-1
!  
!gbs Zooplankton effect
            IF(grazing(i,j).eq.0.0) GOTO 3997
            IF((wqual(i,j)-grazing(i,j)).le.alg_min) THEN
               grazing(i,j) = 0.9*(wqual(i,j)-alg_min)
			      IF(grazing(i,j).lt.0.0) grazing(i,j) = 0.0
            ENDIF
3997        CONTINUE
	         wqual(i,j) = wqual(i,j) - grazing(i,j)!  
!   Estimate the Minimum alloable value  
            Minimum_Nut =dmin1(wqual(i,j)-grazing(i,j),                                      &
     	      wqual(i,11)+ ypchla(j)*growth(i,j)-ypchla(j)*death(i,j)-ypchla(j)*grazing(i,j),  &
     	      wqual(i,17)+ ynchla(j)*growth(i,j)-ynchla(j)*death(i,j)-ynchla(j)*grazing(i,j))
            IF(Minimum_Nut .lt. 0.0)THEN
               grazing(i,j) = 0.9*dmin1(wqual(i,j),                                          &
     	         (wqual(i,11)+ ypchla(j)*growth(i,j)-ypchla(j)*death(i,j))/ypchla(j),          &
     	         (wqual(i,17)+ ynchla(j)*growth(i,j)-ynchla(j)*death(i,j))/ynchla(j))
 	            IF(grazing(i,j).lt.0.0) grazing(i,j) = 0.0
	         ENDIF      
            Minimum_Nut = dmin1(wqual(i,j)-grazing(i,j),                                     &
     	      wqual(i,11)+ ypchla(j)*growth(i,j)-ypchla(j)*death(i,j)-ypchla(j)*grazing(i,j),  &
     	      wqual(i,17)+ ynchla(j)*growth(i,j)-ynchla(j)*death(i,j)-ynchla(j)*grazing(i,j))
	         IF(Minimum_Nut .lt. 0.0)THEN
	            WRITE(*,*) 'Error in Algae.Graz too high will give <0 PhytoN '	
               WRITE(*,*) 'Minimum_Nut',Minimum_Nut
	            PAUSE
            ENDIF
!**************Final Chlorophyl in any layer******************************
!gbs 	  wqual(i,j) = wqual(i,j) + growth(i,j)- death(i,j) - grazing(i,j) 
!**************************************************************************		
	      ENDDO ! j chl
	   ENDDO	 ! i layers
!	   WRITE(11111,3)jday,temp(ns-1),wqual(ns-1,1),growth(ns-1,1),death(ns-1,1),depth(ns)
      GOTO 5000
!sam flux algae from sediments (for someone ELSE to DO)...
!	call alg_sed(nchl)

!   Check for algae minimas...
!
	   DO i = 1,ns
	      DO j = 29, 29+nchl_species-1
	         IF(wqual(i,j) .lt. alg_min) THEN
!	            WRITE(*,*)'ERROR in Algae. Chl<0.0 Layer ',i,'Jday',jday
		         wqual(i,j) = alg_min
!	            pause
		      ENDIF
		
	         IF(wqual(i,11).lt.ypchla(j)*alg_min) THEN 
!    		      WRITE(*,*)'ERROR in Algae. P < Pmin',i,'Jday',jday,wqual(i,11)
		         wqual(i,11) = alg_min*ypchla(j)
!	            pause		
		      ENDIF	    
		      IF(wqual(i,17).lt.ynchla(j)*alg_min) THEN
!		         WRITE(*,*)'ERROR in Algae. Nmin < min',i,'Jday',jday,wqual(i,17)
		         wqual(i,17) = alg_min*ynchla(j)
!	            pause       	
		      ENDIF
	      ENDDO ! j chl
	   ENDDO	 ! i layers
5000  CONTINUE

	   Growth_Top = Growth_Top + growth(ns,1)
	   Graz_Top   = Graz_Top   + grazing(ns,1)
	   Death_Top  = Death_Top  + death(ns,1)

	   Growth_Top_2 = Growth_Top_2 + growth(ns,1)
	   Graz_Top_2   = Graz_Top_2   + grazing(ns,1)
	   Death_Top_2  = Death_Top_2  + death(ns,1)

	   IF(iclock+nosecs.eq.itmpr) THEN   
	      Chl_Top    = 0.0
	      Growth_Top = 0.0
	      Graz_Top   = 0.0
	      Death_Top  = 0.0
	      Growth_Top_2 = 0.0
	      Graz_Top_2   = 0.0
	      Death_Top_2  = 0.0
      ENDIF
!	close(11111)
	   RETURN
	   END SUBROUTINE ALGAE_SPECIES

!************************************************************************
	   SUBROUTINE ALGAE_SPECIES_SETTLE(lost)
!***********************************************************************
!     POP and PON settling are in Nutrient Settling subroutine
!     arfac = multiplicative factor to get area to m**2
!     chf(maxns,nchl) = algae cells per layer volume [cells/vol]
!     chfdv(maxns,nchl) = number of cells in a layer [cells/area]
!     dep_alg = depth of oft expanded algae layer
!     depleft(maxns,nchl) = distance of top of layer from bottom after settling
!     deplefb(maxns,nchl) = distance of bottom of layer from bottom after settling
!     depth(maxns) = depth array for layers relative to bottom
!     depthm(maxns) = mean depth of a layer (m) i.e. between i and i-1
!     dep(maxns) = array for depth between layers (m)
!     fluxin(maxns,nchl) = flux of particles into layer
!     ns = number of dyresm layers
!     nchl = number of non-zero algae modeled
!     paro = light level at top of each layer [W/m2 PAR]
!     setldep(maxchl) = vert. settling depth in a time step [m]
!     setl_vel(maxchl) = settling velocities [m/day]
      
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	
      
      	
	   INTEGER*4 i,j,k,n,ii,l
	   REAL*8 arfac
	   REAL*8 chf(maxns,3*3)
	   REAL*8 chfdv(maxns,3*3)
	   REAL*8 dep(maxns)
	   REAL*8 depleft(maxns,3)
	   REAL*8 deplefb(maxns,3)
	   REAL*8 fluxin(maxns,3*3)
	   REAL*8 setldep(maxchl)
	   REAL*8 lt_threshold, dep_alg(maxns,maxchl), delta_a
	   REAL*8 Before(3*3),After(3*3),Lost(3*3)
	   REAL*8 Delta_area
	   PARAMETER(arfac = 1.0d+6)
	   PARAMETER(lt_threshold = 15.0/2.07*.45) 	!!sam also in algae
         
!			***** initializations ******
!   get incremental depth of each layer
!
	   DO i = ns,2,-1
	      dep(i) = depth(i)-depth(i-1)
	   ENDDO
	   dep(1) = depth(1)
!  
!     Initialize the number of particles/m3 at every box. convert the
!     concentration of phytoplankton to numbers
!     wqual 11,12,17,18 are the internal and detrital P & N concentrations,
!     respectively, in [ug/L].
!     Assign a multiplicative factor to the Detrital nutrient concentrations.
!     Note tha wqual 13 and 19 stands for POP and PON respectively
!  
	   DO k = 29, 29+nchl_species-1
!     settling - setldep is the distance settled meters/dt = sec/dt*m/day*day/sec
	      setldep(k)=dfloat(nosecs)*setl_vel(k)/86400.
!	      IF (abs(setldep(k)) .lt. 1e-4) cycle
	      DO i = 1,ns
		      chf(i,k) =   wqual(i,k)    * part_fac(k)  	! Chlorophyll conc.
!gbs		   chf(i,k+3) = wqual(i,10+k) * part_fac(k)	! Phyto P 
!gbs		   chf(i,k+6) = wqual(i,16+k) * part_fac(k)	! Phyto N 

!gbs	      DO ii=1,3	! conc., int. P, int. N
!gbs	      chfdv(i,k+(ii-1)*3) = chf(i,k+(ii-1)*3)*dep(i) ! parts/m2
            chfdv(i,k) = chf(i,k)*dep(i) ! parts/m2
            IF(i.eq.1)THEN
		         depleft(i,k) = depth(1) - setldep(k)
		         deplefb(i,k) = 0.0 - setldep(k)
	         ELSE
		         depleft(i,k) = depth(i) - setldep(k)
		         deplefb(i,k) = depth(i-1) - setldep(k)
	         ENDIF
	         dep_alg(i,k)=depleft(i,k)-deplefb(i,k)
	         fluxin(i,k)=0.0			! init.
!gbs	   fluxin(i,k+(ii-1)*3)=0.0			! init.
!gbs	  ENDDO    
	      ENDDO	! i layers
	   ENDDO	! k algae
!  
!   Calculate the initial mass of WQ variables
!  
      DO k = 29, 29+nchl_species-1  
	      Before(k)   = 0.0
	      After(k)    = 0.0
	      Lost(k)     = 0.0
      ENDDO    ! k

 	   DO k = 29, 29+nchl_species-1
         DO i = 1,ns
	         Before(k)   = Before(k)   + wqual(i,k)   *vol(i)    
!gbs	      Before(k+3) = Before(k+3) + wqual(i,10+k)*vol(i) 
!gbs	      Before(k+6) = Before(k+6) + wqual(i,16+k)*vol(i) 
        ENDDO  ! i
      ENDDO    ! k

!		******** start of settling loop **********
	   DO k = 29, 29+nchl_species-1
	      IF (abs(setldep(k)) .lt. 1e-4) cycle
	      DO i = ns,1,-1		! algae layer index 
		      DO j = ns,1,-1	    ! model depth index, top to bottom 
!   Area correction factor  
               IF(i.ne.1.and.j.ne.1) THEN
                  IF((area(i) + area(i-1)).lt.(area(j) + area(j-1))) THEN
	                  Delta_area = (area(i) + area(i-1))/(area(j) + area(j-1))
	               ELSE
	                  Delta_area =  1.0
	               ENDIF
	            ELSEIF(i.eq.1.and.j.eq.1) THEN
                  Delta_area =  1.0
	            ELSEIF(i.eq.1.and.j.ne.1) THEN
	               IF(area(i).lt.(area(j) + area(j-1))) THEN
	                  Delta_area = area(i)/(area(j) + area(j-1))
	               ELSE
	                  Delta_area=  1.0
	               ENDIF
	            ELSEIF(i.ne.1.and.j.eq.1)THEN
	               IF(area(i)+area(i-1).lt.area(j)) THEN
	                  Delta_area = (area(i) + area(i-1))/area(j)
	               ELSE
	                  Delta_area =  1.0
	               ENDIF
	            ENDIF
		         IF(j.eq.1)THEN	! look at lake layer 1 (bottom)
!   case where settling particles fall within layer 1 and extend further than the bottom:
		            IF(depleft(i,k).le.depth(j).and.depleft(i,k).gt.0.0.and.deplefb(i,k).le.0.0)THEN

!gbs			         DO ii=1,3
!gbs			            fluxin(j,k+(ii-1)*3)=fluxin(j,k+(ii-1)*3)+                       &
!gbs     			         (depleft(i,k)/dep_alg(i,k))*chfdv(i,k+(ii-1)*3)*Delta_area
!gbs			         ENDDO
			            fluxin(j,k)=fluxin(j,k)+ (depleft(i,k)/dep_alg(i,k))*chfdv(i,k)*Delta_area
!   case where settling particles extend above the top of layer 1 and extend further than the bottom:
		            ELSEIF(depleft(i,k).gt.depth(j).and.deplefb(i,k).le.0.0)THEN

!gbs			         DO ii=1,3
!gbs			            fluxin(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)+           &
!gbs     			      (dep(j)/dep_alg(i,k))*chfdv(i,k+(ii-1)*3)*Delta_area
!gbs			         ENDDO
			            fluxin(j,k) = fluxin(j,k)+(dep(j)/dep_alg(i,k))*chfdv(i,k)*Delta_area
!   case where all settling particles fall within layer 1:
		            ELSEIF(depleft(i,k).le.depth(j).and.deplefb(i,k).gt.0.0)THEN
!gbs			         DO ii=1,3
!gbs			            fluxin(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)+          &
!gbs    			            chfdv(i,k+(ii-1)*3)*Delta_area
!gbs			         ENDDO
			            fluxin(j,k) = fluxin(j,k)+chfdv(i,k)*Delta_area
!   case where settling particles extend above layer 1 and fall within layer 1
		            ELSEIF(depleft(i,k).gt.depth(j).and.deplefb(i,k).lt.depth(j).and.deplefb(i,k).gt.0.0)THEN
!gbs			         DO ii=1,3
!gbs		   	         fluxin(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)+             &
!gbs     	   	         ((depth(j)-deplefb(i,k))/dep_alg(i,k))*chfdv(i,k+(ii-1)*3)*Delta_area
!gbs			         ENDDO
 		   	         fluxin(j,k) = fluxin(j,k)+((depth(j)-deplefb(i,k))/dep_alg(i,k))*chfdv(i,k)*Delta_area  
		            ENDIF

! ****** rising algae ********************************* j > 1
		         ELSEIF (depleft(i,k).gt.depth(i)) THEN
! average area difference...
		            IF (i.eq.1) THEN
			            delta_a = area(i)/(area(j)+area(j-1))
		            ELSE
			            delta_a = (area(i)+area(i-1))/(area(j)+area(j-1))
		            ENDIF
		   
		            IF( (depleft(i,k).le.depth(j) .or. j.eq.ns) .and.deplefb(i,k).ge.depth(j-1))THEN
!gbs			         DO ii=1,3
!gbs			            fluxin(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)+          &
!gbs     		               chfdv(i,k+(ii-1)*3) * delta_a
!gbs			         ENDDO
			            fluxin(j,k) = fluxin(j,k)+ chfdv(i,k) * delta_a	     		   
		            ELSEIF(depleft(i,k).ge.depth(j).and.deplefb(i,k).le.depth(j-1))THEN
		               IF (j.eq.ns) THEN
!gbs			            Do ii=1,3
!gbs			               fluxin(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)        &
!gbs     			            +(depleft(i,k)-depth(j-1))/dep_alg(i,k)*chfdv(i,k+(ii-1)*3) * delta_a
!gbs			            ENDDO
			               fluxin(j,k) = fluxin(j,k)+(depleft(i,k)-depth(j-1))   &
     			                        /dep_alg(i,k)*chfdv(i,k) * delta_a	 
		               ELSE
!gbs			            DO ii=1,3
!gbs			               fluxin(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)+       &
!gbs     			            (dep(j)/dep_alg(i,k))*chfdv(i,k+(ii-1)*3)* delta_a
!gbs			            ENDDO
			               fluxin(j,k) = fluxin(j,k)+(dep(j)/dep_alg(i,k))*chfdv(i,k)* delta_a 
		               ENDIF
		            ELSEIF(depleft(i,k).gt.depth(j).and.deplefb(i,k).lt.        &
		            		   depth(j).and.deplefb(i,k).ge.depth(j-1))THEN
!gbs			         DO ii=1,3
!gbs			            fluxin(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)+          &
!gbs     		            ((depth(j)-deplefb(i,k))/dep_alg(i,k))*chfdv(i,k+(ii-1)*3) * delta_a
!gbs			         ENDDO
			            fluxin(j,k) = fluxin(j,k)+ ((depth(j)-deplefb(i,k))      &
     		                        /dep_alg(i,k))*chfdv(i,k) * delta_a

		            ELSEIF(depleft(i,k).le.depth(j).and.depleft(i,k).gt.        &
     	                  depth(j-1).and.deplefb(i,k).le.depth(j-1))THEN
!gbs			         DO ii=1,3
!gbs			            fluxin(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)+          &
!gbs     		               ((depleft(i,k)-depth(j-1))/dep_alg(i,k))*chfdv(i,k+(ii-1)*3) * delta_a
!gbs		            ENDDO 			
			            fluxin(j,k) = fluxin(j,k)+ ((depleft(i,k)-depth(j-1))    &
      		                        /dep_alg(i,k)) *chfdv(i,k) * delta_a 
		            ENDIF

! ****** sinking algae ********
		         ELSE
!   case where all settling particles fall within the top and bottom of a layer
		            IF( depleft(i,k).le.depth(j) .and.deplefb(i,k).ge.depth(j-1))THEN
!gbs			         DO ii=1,3
!gbs			            fluxin(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)+       &
!gbs     			         chfdv(i,k+(ii-1)*3)*Delta_area
!gbs			         ENDDO
		               fluxin(j,k) = fluxin(j,k)+ chfdv(i,k)*Delta_area	    			

!   case where settling particles extend above the top of a layer and
!   below the bottom of the layer
		            ELSEIF(depleft(i,k).ge.depth(j).and.deplefb(i,k).le.depth(j-1))THEN
!gbs			         DO ii=1,3
!gbs			            fluxin(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)+       &
!gbs     			         (dep(j)/dep_alg(i,k))*chfdv(i,k+(ii-1)*3)*Delta_area
!gbs			         ENDDO
 		               fluxin(j,k) = fluxin(j,k)+	(dep(j)/dep_alg(i,k))*chfdv(i,k)*Delta_area
!   case where settling particles extend above the top of a layer but
!   bottom is within the layer

		            ELSEIF(depleft(i,k).gt.depth(j).and.deplefb(i,k).lt.depth(j)     &
    		           .and.deplefb(i,k).ge.depth(j-1))THEN
!gbs			         DO ii=1,3
!gbs			            fluxin(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)+                &
!gbs     	               ((depth(j)-deplefb(i,k))/dep_alg(i,k))*chfdv(i,k+(ii-1)*3)*Delta_area
!gbs			         ENDDO
			            fluxin(j,k) = fluxin(j,k)+((depth(j)-deplefb(i,k))             &	    
			                          /dep_alg(i,k))*chfdv(i,k)*Delta_area
!   case where top of settling particles fall within a layer but bottom
!   extends below the bottom of the layer

		            ELSEIF(depleft(i,k).le.depth(j).and.depleft(i,k).gt.depth(j-1)    &
    	              .and.deplefb(i,k).le.depth(j-1))THEN	 
!gbs		            DO ii=1,3
!gbs		               fluxin(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)+                &
!gbs     	            ((depleft(i,k)-depth(j-1))/dep_alg(i,k))*chfdv(i,k+(ii-1)*3)*Delta_area
!gbs		            ENDDO
		               fluxin(j,k) = fluxin(j,k)+ ((depleft(i,k)-depth(j-1))/         &
     	                          dep_alg(i,k))*chfdv(i,k)*Delta_area 
		            ENDIF
	            ENDIF
		      ENDDO	! j layers
	      ENDDO	! i layers
!   convert chf back to a concentration
	      DO j = 1,ns
!gbs		   Do ii=1,3
!gbs		      chf(j,k+(ii-1)*3) = fluxin(j,k+(ii-1)*3)/dep(j)
!gbs		   ENDDO
		      chf(j,k) = fluxin(j,k)/dep(j)
	      ENDDO ! j layers	
	   ENDDO	! k algae
!   convert from numbers of particles conc. to a mass conc.
!   correct bod, tp and tn to total values.
!sam set to non-zero defaults...
	   DO k = 29, 29+nchl_species-1
!   unaffected IF the settling velocity is below some threshold...
	      IF (abs(setldep(k)) .lt. 1e-4) cycle
	      DO i = 1,ns
		      wqual(i,k) = chf(i,k) / part_fac(k)
	      ENDDO	! i
	   ENDDO	! k
!  
!   Estimate the lost of mass by setling
!  
	   DO k = 29, 29+nchl_species-1
         DO i = 1,ns
	         After(k)   = After(k)   + wqual(i,k)   *vol(i) 
         ENDDO  ! i
      ENDDO    ! k
	   DO k = 1,nchl
	      Lost(k)   = Before(k)   - After(k)  
      ENDDO    ! k
!  
!   Check for algae minimas...
!  
	   DO i = 1,ns
	      DO k = 29, 29+nchl_species-1
	         IF(wqual(i,j) .lt. alg_min) THEN
!	            WRITE(*,*)'ERROR in Setling algae. Algae conc. less than min',i,jday
		         wqual(i,j) = alg_min
	            ! pause
	         ENDIF		
		   ENDDO ! j chl
	   ENDDO	 ! i layers
	   RETURN
	   END SUBROUTINE ALGAE_SPECIES_SETTLE
!**********************************************************************
	   SUBROUTINE ALGAE_SPECIES_SEDIMENT_RELEASE
!***********************************************************************
! sediment release of algae within the mixed layer
!-----------------------------------------------------------------------
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

	   INTEGER*4 i,j
	   REAL*8 arfac, are(maxns), dep(maxns)
	   REAL*8 Afactor 
	   PARAMETER (arfac = 1.0e+6)

! initializations...
	   DO i = ns,2,-1
	      are(i) = (area(i)-area(i-1))*arfac
	      dep(i) = depth(i)-depth(i-1)
	   ENDDO
	   are(1) = area(1)*arfac
	   dep(1) = depth(1)

!     'Afactor' is for the proportionality between the layer sediment area
!     and its volume
!     'sedalgae' is a flux rate in [mg/m2/day], need to set in *.wat file
!  '  depmx' set in subr. mixer
	   DO i = 1,ns
	      IF (depth(i) .ge. depmx) THEN	! only within mixed layer
		      Afactor = (dfloat(nosecs)/86400.0)*(are(i)/(area(i)*arfac*dep(i)))
		      Do j = 29, 29+nchl_species-1
!		         IF(sedalgae(j).le.alg_min)sedalgae(j) = alg_min
!		         wqual(i,j) = wqual(i,j) + sedalgae(j) * Afactor
		      ENDDO
! set BOD flux rate at 0.5 mg/m2/day...	
!		      wqual(i,9) = wqual(i,9) + 0.5 * Afactor
 	      ENDIF
	   ENDDO
	   RETURN
	   END SUBROUTINE ALGAE_SPECIES_SEDIMENT_RELEASE
!************************************************************************************
