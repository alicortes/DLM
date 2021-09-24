!************************************************************************
      SUBROUTINE ALGAE_TAHOE(par,limname,SDname,latit)
!***********************************************************************  
!  Subroutine to calculate the change in chlorophyll concentration
!--------------------------------------------------------------------------
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE		

!    **** defining variables ****
!     alg_min(j) = minimum algae conc. [ug/L]
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
      REAL*8 C_uptake, P_Uptake,N_Upatke,Si_uptake, C_RandM, P_RandM, N_RandM,Si_RandM,C_Grazing, P_Grazing,N_Grazing, Si_Grazing
!	   REAL*8 etca(3)
	   REAL*8 limiting
	   REAL*8 Minimum_Nut,Minimum_Growth 	
	   REAL*8 par
	   REAL*8 satlit(maxns,maxchl),halfsatc(maxns,maxchl), halfsatn(maxns,maxchl),halfsatp(maxns,maxchl), halfsatsi(maxns,maxchl)
	   REAL*8 seglit(maxns,maxchl)
	   REAL*8 segc(maxns,maxchl),segn(maxns,maxchl),segp(maxns,maxchl),segsi(maxns,maxchl)
       REAL*8 yquanta  !, alg_min(j)
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
       INTEGER*4 mdays,jyear
       
       REAL*8 Vol_Mix_T,Chla_Mix,Chla_Mix_Media
	   REAL*8 latit,kmax(maxchl),temp_min(maxchl),temp_max(maxchl),temp_opt(maxchl),growth_rate(maxns,maxchl),tempgrowth(maxchl), tempdeath(maxchl)
	   REAL*8 kk1, kk2,alphaa,betaa,maxgr(maxchl),Ir(maxns,maxchl)
       REAL*8,  PARAMETER ::Qmcyst_max=0.155, Qmcyst_min=0.050
!       REAL*8,  PARAMETER ::Qmcyst_max=7.6, Qmcyst_min=0.71
!sam	PARAMETER(conquan = 720.0,e = 2.71828,yquanta = 0.96418)
	   PARAMETER(conquan = 0.4, e = 2.71828, yquanta = 0.1)  ! per Brad
	   PARAMETER(lt_threshold = 15.0/2.07*.45)	!!sam
!	PARAMETER(alg_min(j) = 0.05)
 	   zoo  = 0.00000001d0*(dfloat(nosecs)/86400.0d0) !(#/L)
      
      DO j=1,MAXCHL
        IF(j.EQ.4) THEN
            alg_min(j) = 0.10
        ELSEIF(j.NE.4) THEN
            alg_min(j) = 0.10
        ENDIF
      ENDDO
      
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
!          IF(JDAY.EQ.2004170)WQUAL(NS,4)=10.0
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
       mdays  = jday-(jday/1000)*1000
	   DO i = 1,ns
	   
	      DO j = 1,nchl
		      et1(i) = et1(i)+etca(j)*wqual(i,j)
	      ENDDO ! j
	      DO j = 1,7
		      et1(i) = et1(i)+etpart(j)*cf(i,j)
	      ENDDO ! j
	      et1(i) = et1(i)+etwat    
         if(mdays.ge.90.and.mdays.lt.135) then
			   et1(i)=1.70d0   !eta=1.7/secchi depth (1 to 1.5 SD)
		  elseif(mdays.ge.135.and.mdays.lt.155) then
			   et1(i)=1.8d0   !eta=1.7/secchi depth (1 to 1.5 SD)
	     elseif(mdays.ge.155.and.mdays.lt.210) then    ! 160 to 210
			   et1(i)=1.9d0   !eta=1.7/secchi depth (1 to 1.5 SD)
		  elseif(mdays.ge.210.and.mdays.lt.300) then    ! 210 to 250
			   et1(i)=1.9d0   !eta=1.7/secchi depth (1 to 1.5 SD)     
	     elseif(mdays.ge.300.and.mdays.lt.330) then
			   et1(i)=1.8d0   !eta=1.7/secchi depth (1 to 1.5 SD)
        elseif(mdays.ge.330.and.mdays.lt.366) then
	         et1(i)=1.7d0   !eta=1.7/secchi depth (1 to 1.5 SD)
        endif

	   ENDDO   ! i
!  
!quim   Ted Swift clarity model for Tahoe
!  ******************LIGHT_TAHOE**********************************************
!      CALL Light_Tahoe(latit)
!  ****************************************************************************
!     calculate PAR at the surface of each layer
!  
!	   et1(i)=5.2
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
	      Do j = 1,nchl 	
             satlit(i,j) =light(j)*thetat(j)**(temp(i)-20.0)!	         
	      ENDDO	! j
	   ENDDO	! i
!Find the growth reduction dueto insufficient light, i.e. the P-I curve, 
!using the Steele equation modified Jassby and Platt...
 
 	   DO i = ns,1,-1
           den(i)=densty(temp(i),sal(i))
	      DO j = 1,nchl
		      seglit(i,j)=((paro(i)/satlit(i,j))**1.0)*dexp(1.0-(paro(i)/satlit(i,j))**1.0)             
          ENDDO	! j          
        ENDDO ! i
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
!	 	PRINT*,WQUAL(NS,23)	 
	   DO i = 1,ns	   	  
	      DO j = 1,nchl	
		      halfsatc(i,j) =halfc(j)*thetat(j)**(temp(i)-20.0)
              halfsatn(i,j) =halfn(j)*thetat(j)**(temp(i)-20.0)
		      halfsatp(i,j) =halfp(j)*thetat(j)**(temp(i)-20.0)
              halfsatsi(i,j) =halfsi(j)*thetat(j)**(temp(i)-20.0)
              
             segc(i,j) =wqual(i,10)/(halfsatc(i,j)+wqual(i,10)) 
	         segp(i,j) =wqual(i,16)/(halfsatp(i,j)+wqual(i,16)) 	
	         segn(i,j) =(wqual(i,21)+wqual(i,22))/(halfsatn(i,j)+(wqual(i,21)+wqual(i,22)))
             segsi(i,j) =wqual(i,27)/(halfsatsi(i,j)+wqual(i,27))
!             print*,segc(i,j),segp(i,j),segn(i,j),seglit(i,j)
!2015/08
             if(j.eq.1) then                 
                 segmin(i,j) = dmin1(seglit(i,j),segc(i,j), segp(i,j),segn(i,j),segsi(i,j))
!                  segmin(i,j) = dmin1(seglit(i,j),segp(i,j),segn(i,j))
             
             elseif(j.eq.4) then
!                 if(paro(i).gt.light(j)) seglit(i,j)=0.0
                 segmin(i,j) = dmin1(seglit(i,j),segc(i,j),segp(i,j))
             elseif(j.gt.1.and.j.ne.4) then                 
	             segmin(i,j) = dmin1(seglit(i,j),segc(i,j),segp(i,j),segn(i,j))
             endif
!		      segmin(i,j) = dmin1(segp(i,j),segn(i,j))
!             if(i.eq.ns.and.j.eq.4) then
!                 write(99,fmt='(2i8,8f15.8)')jday,iclock,segmin(i,j),segc(i,j),segp(i,j),segn(i,j),seglit(i,j),wqual(i,16),wqual(i,21),wqual(i,22)
!              endif
             
	         IF(segmin(i,j).gt.1.0.or.segmin(i,j).lt.0.0) THEN
	            WRITE(*,*)'Warning: ALGAE, minimum value out of range)'
	            WRITE(*,*)'segmin,i,j',segmin,i,j
!	            stop
	         ENDIF
! Minimum threshold value for orthophsophate 0.5 microgram/L
!2015/09	         IF (wqual(i,14).LT.0.9) THEN  !0.5
!2015/09	            segmin(i,j) = 0.0
!2015/09		         GOTO 4098
!2015/09	         ENDIF		
!2015/09		      IF (wqual(i,18).LT.1.5) THEN  
!2015/09	            segmin(i,j) = 0.0
!2015/09		         GOTO 4098
!2015/09	         ENDIF
!2015/09		      IF(iclock+nosecs.eq.itmpr.or.(itimes.eq.1440.and.itmpr.eq.86400)) THEN
!2015/09	            open(202,file=limname(1),status='unknown',access='append')
!2015/09		         IF (segmin(i,j).eq.  segp(i,j)) THEN
!2015/09	               limiting = 0.0
!2015/09		            GOTO 4097
!2015/09		         ENDIF
!2015/09		         IF (segmin(i,j).eq.  segn(i,j)) THEN
!2015/09			         limiting = 1.0
!2015/09			         GOTO 4097
!2015/09		         ENDIF
!2015/09		         IF (segmin(i,j).eq.seglit(i,j)) THEN
!2015/09	               limiting = 2.0
!2015/09			         GOTO 4097
!2015/09		         ENDIF
!2015/09 4097          CONTINUE
!		         IF(depth(i).gt.300) THEN
!			         WRITE(202,1034) jday, depth(i), limiting
!		         ENDIF
!1034		      FORMAT(i9,x,f10.1,x,f10.1) !,x,3(f10.5,x))
!2015/091035		      FORMAT(i10,x,f10.1,x, 3(f10.5,x))
!					CLOSE(202)
!2015/09	         ENDIF
!2015/09 4098	      CONTINUE
	         IF (segmin(i,j) .le. 0.0) THEN
               segmin(i,j) = 0.0
	         ENDIF
!  
!gbs	grow of the algae  
!   
!gbs******MINLAKE Algorithm for temperature multiplier for Algae growth 
!***********************2004 -----------------------------
            temp_opt(1)=21.0d0; temp_min(1)=20.0d0; temp_max(1)=29.0d0      !Diatoms  23 10  24
            temp_opt(2)=21.0d0; temp_min(2)=20.0d0; temp_max(2)=29.0d0      !Greens
            temp_opt(3)=21.0d0; temp_min(3)=20.0d0; temp_max(3)=29.0d0      !Cryptophytes 25.5 15 25.55
            temp_opt(4)=25.0d0; temp_min(4)=20.0d0; temp_max(4)=30.0d0      !Aphanozomenon 25 23 28
            temp_opt(5)=26.0d0; temp_min(5)=22.0d0; temp_max(5)=27.0d0      !Microcystis
            temp_opt(6)=25.0d0; temp_min(6)=22.0d0; temp_max(6)=30.0d0      !Microcystin
            temp_opt(7)=25.0d0; temp_min(7)=22.0d0; temp_max(7)=30.0d0      !Uknown          
!*********************************2010 to 2012 *******************************************************************            
 !           temp_opt(1)=08.0d0; temp_min(1)=03.0d0; temp_max(1)=15.0d0      !Diatoms  23 10  24
 !           temp_opt(2)=25.0d0; temp_min(2)=10.0d0; temp_max(2)=27.0d0      !Greens
 !           temp_opt(3)=20.0d0; temp_min(3)=08.0d0; temp_max(3)=27.0d0      !Cryptophytes 25.5 15 25.55
 !           temp_opt(4)=25.0d0; temp_min(4)=22.0d0; temp_max(4)=27.0d0      !Aphanozomenon 25 23 28
 !           temp_opt(5)=25.0d0; temp_min(5)=22.0d0; temp_max(5)=30.0d0      !Microcystis
 !           temp_opt(6)=25.0d0; temp_min(6)=25.0d0; temp_max(6)=30.0d0      !Microcystin
 !           temp_opt(7)=25.0d0; temp_min(7)=20.0d0; temp_max(7)=30.0d0      !Uknown                     
!******************************************************************************************************  
              if (temp(i).ge.3.0.and.temp(i).le.15.0) then                 
                   temp_opt(2)=05.0d0; temp_min(2)=4.0d0; temp_max(2)=7.0d0
                                      
               elseif(temp(i).ge.15.0.and.temp(i).le.30.0) then

                    temp_opt(2)=21.0d0; temp_min(2)=20.0d0; temp_max(2)=29.0d0               
               endif 
               
               if (temp(i).ge.20.0.and.temp(i).le.30.0) then                 
                  temp_opt(3)=21.0d0; temp_min(3)=20.0d0; temp_max(3)=29.0d0
                                      
               elseif(temp(i).ge.5.0.and.temp(i).le.25.0) then

                    temp_opt(3)=20.0d0; temp_min(3)=08.0d0; temp_max(3)=27.0d0      !Cryptophytes 25.5 15 25.55               
               endif 
               
           IF(temp(i).le.temp_opt(j)) THEN	   
	            tempgrowth(j) =dexp(-2.3*((temp(i)-temp_opt(j))/(temp_opt(j)-temp_min(j)))**2.0)	 
           ELSEIF(temp(i).gt.temp_opt(j)) THEN
	            tempgrowth(j) =dexp(-2.3*((temp(i)-temp_opt(j))/(temp_max(j)-temp_opt(j)))**2.0)	 
          	ENDIF
!************************ 
!2015/08            IF(temp(i).le.temp_opt) THEN	   !!  temp(i)
!2015/08	            tempgrowth(j)=1.0*dexp(-0.01d0*((temp(i)-temp_opt)**2.0))    !0.005  0.002
!2015/08	         ELSEIF(temp(i).gt.temp_opt) THEN
!2015/08	            tempgrowth(j)=1.0*dexp(-0.01d0*((temp_opt-temp(i))**2.0))	 !0.005  0.002		     
!2015/08       	   ENDIF
!-------------Peters and Eilers, 1978 Hydrobiological Bulletin 12: 127-134----	
!	         maxgr=0.23d0	!100.5
	         maxgr(j) = gromax (j)
	         betaa = et1(i)  !
!	         Ir(i,j)=paro(i)/satlit(i,j)
             Ir(i,j)=paro(i)/light(j)
!	         IF(Ir(i,j).gt.0.0d0.and.paro(i).le.satlit(i,j)) THEN
!	            growth_rate(i,j)=maxgr
!	         ELSE
	         growth_rate(i,j)=maxgr(j)*2.0d0*(1.0d0+betaa)*Ir(i,j)/      &
                             (Ir(i,j)*Ir(i,j)+2.0d0*betaa*Ir(i,j)+1.0d0)
             
!             if(j.eq.4.or.j.eq.5.and.paro(i).lt.light(j))  growth_rate(i,j)=0.0
             
!	         ENDIF
!            if(i.eq.ns.and.j.eq.4.and.jday.ge.2004170.and.jday.lt.2004300) then
!            write(999999,fmt='(3i8, 4f12.5)') jday,iclock,j, depth(i), paro(i),Ir(i,j),growth_rate(i,j)
!            endif
!-------------Steel (1962), Limnology and Oceanography----
!	         maxgr=1.00
!	         Ir(i,j)=paro(i)/satlit(i,j)
! 	         growth_rate(i,j)=maxgr*(1.0d0/satlit(i,j))*dexp(1-Ir(i,j))	
!**************************************************************
!	         growth(i,j) = (wqual(i,j)*growth_rate(i,j)*tempgrowth*seglit(i,j)
!2015/09/02 *************Introduction of microcystin************************
!2015/09/02  J=6 ASSIGNED MICROSYSTIN       
!            IF(J.EQ.6) THEN
!                growth(i,j) = growth_rate(i,j)*((growth_rate(i,j)*Qmcyst_max - growth_rate(i,j)*Qmcyst_min)/maxgr(j)+  &
!	                        Qmcyst_min)*tempgrowth(j)*segmin(i,j)*faco(i)*(dfloat(nosecs)/86400.0)                
!                growth(i,j) = min((wqual(i,5)*growth_rate(i,5)*(Qmcyst_max - Qmcyst_min)*tempgrowth(5)*segmin(i,j)*faco(i)*(dfloat(nosecs)/86400.0)), &
!                                    wqual(i,5)*growth_rate(i,5)*(Qmcyst_max - Qmcyst_min)*tempgrowth(5)*segmin(i,j)*faco(i)*(dfloat(nosecs)/86400.0))
!               wqual(i,j)=max(wqual(i,5)*Qmcyst_max,wqual(i,5)*Qmcyst_min)
!            ELSEIF (J.NE.6) THEN   
!               if(j.eq.4.and.jday.ge.2004220)tempgrowth(j)=tempgrowth(j)/6.0
               growth(i,j) = wqual(i,j)*growth_rate(i,j)*tempgrowth(j)*segmin(i,j)*faco(i)*(dfloat(nosecs)/86400.0)  	                        
!            ENDIF
            
!**********************Aphanozomenon's NITROGEN FIXATION FROM ATMOSPHERE**********************            
!**********************Aphanozomenon's NITROGEN FIXATION FROM ATMOSPHERE**********************            
            IF(wqual(i,4).gt.0.2) THEN
              wqual(ns,21)=wqual(ns,21)+1.26*(dfloat(nosecs)/3600.0)
              wqual(i,22)=wqual(i,22)+1.50*growth(i,4)*(dfloat(nosecs)/3600.0)
            ENDIF
   
            
!	         growth(i,j) = (wqual(i,j)*gromax(j)*thetat(j)**                         &
!      	                 (temp(i)-20.0)*segmin(i,j)*faco(i)*(nosecs/86400.0))
!  
!    To avoid less than zero Nutrients!  
!    Preferential ammonium uptake factor..(pg 270 Bowie et al.1985).
!
            IF((wqual(i,21) + wqual(i,22)).le. 0.0) THEN
               hscnhno = 0.0 
	         ELSE
		         hscnhno = wqual(i,21)*wqual(i,22)/((hscn(j)+wqual(i,22))*(hscn(j)+wqual(i,21)))+  &
     	                   wqual(i,22)*hscn(j)/((wqual(i,21)+wqual(i,22))*(hscn(j)+wqual(i,21)))
            ENDIF
	         IF(hscnhno.gt.1.or.hscnhno.lt.0) THEN
		         WRITE(*,*)'WARNING. In Routine ALGAE'
		         WRITE(*,*)'ammonium uptake factor (hscnhno) out of range'
		         WRITE(*,*)'hscnhno,i,j',hscnhno,i,j
!		         stop
	         ENDIF
!  
! Estimate the Minimum alloable value
!  
            Minimum_Nut = dmin1(wqual(i,16)-ypchla(j)*growth(i,j),               &
     	                          wqual(i,21)-ynchla(j)*growth(i,j)*(1-hscnhno),   &
                                wqual(i,22)-ynchla(j)*growth(i,j)*hscnhno)

            IF(Minimum_Nut .lt. 0.0)THEN

               growth(i,j) = 0.9*(dmin1(wqual(i,16)/ypchla(j),                   &
      	                               wqual(i,21)/(ynchla(j)*(1-hscnhno)),     &
                                        wqual(i,22)/(ynchla(j)*hscnhno)))
	         ENDIF
	      
	         Minimum_Nut = dmin1(wqual(i,16)-ypchla(j)*growth(i,j),               &
     	                          wqual(i,21)-ynchla(j)*growth(i,j)*(1-hscnhno),   &
                                wqual(i,22)-ynchla(j)*growth(i,j)*hscnhno)
	
	         IF(Minimum_Nut .lt. 0.0)THEN
	            WRITE(*,*) 'Error in Algae.Growth too high will give <0 Nut values'	
               WRITE(*,*) Minimum_Nut
               WRITE(*,*) wqual(i,16)-ypchla(j)*growth(i,j)              !THP
               WRITE(*,*) wqual(i,21)-ynchla(j)*growth(i,j)*(1-hscnhno)  !NO3
               WRITE(*,*) wqual(i,22)-ynchla(j)*growth(i,j)*hscnhno      !NH4
!	            STOP
	         ENDIF
  	       
!             if(j.eq.6) wqual(i,28)=wqual(i,28)+0.50*growth(i,j)
            IF(J.EQ.6) THEN
!                growth(i,j) = growth_rate(i,j)*((growth_rate(i,j)*Qmcyst_max - growth_rate(i,j)*Qmcyst_min)/maxgr(j)+  &
!	                        Qmcyst_min)*tempgrowth(j)*segmin(i,j)*faco(i)*(dfloat(nosecs)/86400.0)                
!                growth(i,j) = min((wqual(i,5)*growth_rate(i,5)*(Qmcyst_max - Qmcyst_min)*tempgrowth(5)*segmin(i,j)*faco(i)*(dfloat(nosecs)/86400.0)), &
!                                    wqual(i,5)*growth_rate(i,5)*(Qmcyst_max - Qmcyst_min)*tempgrowth(5)*segmin(i,j)*faco(i)*(dfloat(nosecs)/86400.0))
!               wqual(i,j)=max(wqual(i,5)*Qmcyst_max,wqual(i,5)*Qmcyst_min)
            ELSEIF (J.NE.6) THEN                
                wqual(i,j) = wqual(i,j) + growth(i,j)                        
            ENDIF
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
	      DO j = 1,nchl
            tempdeath(1)=dexp(0.08d0*(temp(i)-temp_max(j)))
            tempdeath(2)=dexp(0.08d0*(temp(i)-temp_max(j)))
            tempdeath(3)=dexp(0.08d0*(temp(i)-temp_max(j)))
            tempdeath(4)=dexp(0.10d0*(temp(i)-temp_max(j)))
            tempdeath(5)=dexp(0.08d0*(temp(i)-temp_max(j)))
            tempdeath(6)=dexp(0.08d0*(temp(i)-temp_max(j)))
            tempdeath(7)=dexp(0.08d0*(temp(i)-temp_max(j)))
 !           temp_opt(4)=25.0d0; temp_min(4)=10.0d0; temp_max(4)=27.0d0      !Aphanozomenon 25 23 28
 !           IF(j.eq.4.and.temp(i).le.temp_opt(j)) THEN	   
 !	            tempdeath(j) =1.0-dexp(-2.3*((temp(i)-temp_opt(j))/(temp_opt(j)-temp_min(j)))**2.0)	 
 !          ELSEIF(j.eq.4.and.temp(i).gt.temp_opt(j)) THEN
 !	            tempdeath(j) =1.0-dexp(-2.3*((temp(i)-temp_opt(j))/(temp_max(j)-temp_opt(j)))**2.0)	 
 !         	ENDIF
          	 
            
     	      death(i,j) = (constr(j)+constm(j))*tempdeath(j)*		&   !0.069
                         (dfloat(nosecs)/86400.0)*wqual(i,j)		
            IF((wqual(i,j)-death(i,j)).lt.alg_min(j)) THEN
               death(i,j) = 0.9*(wqual(i,j)-alg_min(j))
		         IF(death(i,j).lt.0.0) death(i,j) = 0.0
            ENDIF
!  
!           Estimate the Minimum alloable value
!  
            Minimum_Nut = dmin1(wqual(i,j)-death(i,j),                                   &
     	                          wqual(i,17)+ ypchla(j)*growth(i,j)-ypchla(j)*death(i,j), &
     	                          wqual(i,23)+ ynchla(j)*growth(i,j)-ynchla(j)*death(i,j))

            IF(Minimum_Nut .lt. 0.0)THEN  
               death(i,j) = 0.9*(dmin1(wqual(i,j),                                        &
     	                                (wqual(i,17)+ypchla(j)*growth(i,j))/ypchla(j),      &
      	                             (wqual(i,23)+ynchla(j)*growth(i,j))/ynchla(j)))	 
	            IF(death(i,j).lt.0.0) death(i,j) = 0.0	 
	         ENDIF
            Minimum_Nut = dmin1(wqual(i,j)-death(i,j),                                    &
     	                          wqual(i,17)+ ypchla(j)*growth(i,j)-ypchla(j)*death(i,j),  &
     	                          wqual(i,23)+ ynchla(j)*growth(i,j)-ynchla(j)*death(i,j))

	         IF(Minimum_Nut .lt. 0.0)THEN   
	            WRITE(*,*) 'Error in Algae.Death too high will give <0 PhytoN '	
               WRITE(*,*) 'Minimum_Nut',Minimum_Nut
               WRITE(*,*) 'chla-death',wqual(i,j)-death(i,j)     !Chla
               WRITE(*,*) 'PhytoP-Death',wqual(i,17)-(ypchla(j)*growth(i,j)-ypchla(j)*death(i,j))  !PhytoP
               WRITE(*,*) 'PhtyoN-Death',wqual(i,23)-(ynchla(j)*growth(i,j)-ynchla(j)*death(i,j))  !PhYtoN
	            PAUSE
             ENDIF
!2015/09 toxins from microcystin ********************
!2015/09 
                !2015/09/02     J=5 ASSIGNED CYANOPHYTES J=6 ASSIGNED MICROSYSTIN      
 !           IF(J.EQ.6) THEN
 !               wqual(i,28) = wqual(i,28)+0.50*(max(death(i,5)*Qmcyst_max,death(i,5)*Qmcyst_min))
 !               wqual(i,j) = wqual(i,j)- 0.50*(max(death(i,5)*Qmcyst_max,death(i,5)*Qmcyst_min))
 !               if(wqual(i,j).le.0.1) wqual(i,j)=0.1
 !           ELSEIF (J.NE.6) THEN                
 !               wqual(i,j) = wqual(i,j) - death(i,j)  
 !           ENDIF          
!	      ENDDO ! j chl
!       ENDDO    ! i layers
!  2015/10 ADSORPTION AND DEGRADATION******
!2015/09 toxins from microcystin ********************
!2015/09 
                !2015/09/02     J=5 ASSIGNED CYANOPHYTES J=6 ASSIGNED MICROSYSTIN   
           IF(WQUAL(I,28).LE.0.00D0) then
               WQUAL(I,28)=0.00
           ENDIF    
            IF(J.EQ.6) THEN
                wqual(i,28) = wqual(i,28)+0.5*death(i,5)*Qmcyst_max +0.5*growth(i,5)*Qmcyst_max
                wqual(i,j)=max(wqual(i,5)*Qmcyst_max,wqual(i,5)*Qmcyst_min)                 
            ELSEIF (J.NE.6) THEN                
                wqual(i,j) = wqual(i,j) - death(i,j)  
            ENDIF
	      ENDDO ! j chl
       ENDDO    ! i layers
       DO i = 1,ns
           WQUAL(I,28)=WQUAL(I,28)-0.05*((wqual(i,28))**0.673)*(FLOAT(nosecs)/86400.0d0)  !ADSORPTION TO SOIL 0.0008
           WQUAL(I,28)=WQUAL(I,28)-WQUAL(I,28)*DEXP(-0.3*temp(i))*(FLOAT(nosecs)/86400.0d0)    !DEGRADATION
           IF(WQUAL(I,28).LE.0.0D0) then
               WQUAL(I,28)=0.00
           ENDIF
       ENDDO
       
   
!gbs account for grazing
!  
	   DO i = 1,ns
	      DO j = 1,nchl
!  *
!gbs Zooplankton effect
            IF(grazing(i,j).eq.0.0) GOTO 3997
            IF((wqual(i,j)-grazing(i,j)).le.alg_min(j)) THEN
               grazing(i,j) = 0.9*(wqual(i,j)-alg_min(j))
			      IF(grazing(i,j).lt.0.0) grazing(i,j) = 0.0
            ENDIF
3997        CONTINUE
	         wqual(i,j) = wqual(i,j) - grazing(i,j)!  
!   Estimate the Minimum alloable value  
            Minimum_Nut =dmin1(wqual(i,j)-grazing(i,j),                                      &
     	      wqual(i,17)+ ypchla(j)*growth(i,j)-ypchla(j)*death(i,j)-ypchla(j)*grazing(i,j),  &
     	      wqual(i,23)+ ynchla(j)*growth(i,j)-ynchla(j)*death(i,j)-ynchla(j)*grazing(i,j))
            IF(Minimum_Nut .lt. 0.0)THEN
               grazing(i,j) = 0.9*dmin1(wqual(i,j),                                          &
     	         (wqual(i,17)+ ypchla(j)*growth(i,j)-ypchla(j)*death(i,j))/ypchla(j),          &
     	         (wqual(i,23)+ ynchla(j)*growth(i,j)-ynchla(j)*death(i,j))/ynchla(j))
 	            IF(grazing(i,j).lt.0.0) grazing(i,j) = 0.0
	         ENDIF      
            Minimum_Nut = dmin1(wqual(i,j)-grazing(i,j),                                     &
     	      wqual(i,17)+ ypchla(j)*growth(i,j)-ypchla(j)*death(i,j)-ypchla(j)*grazing(i,j),  &
     	      wqual(i,23)+ ynchla(j)*growth(i,j)-ynchla(j)*death(i,j)-ynchla(j)*grazing(i,j))
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
	      DO j = 1,nchl
	         IF(wqual(i,j) .lt. alg_min(j)) THEN
!	            WRITE(*,*)'ERROR in Algae. Chl<0.0 Layer ',i,'Jday',jday
		         wqual(i,j) = alg_min(j)
!	            pause
		      ENDIF
		
	         IF(wqual(i,16).lt.ypchla(j)*alg_min(j)) THEN 
!    		      WRITE(*,*)'ERROR in Algae. P < Pmin',i,'Jday',jday,wqual(i,11)
		         wqual(i,16) = alg_min(j)*ypchla(j)
!	            pause		
		      ENDIF	    
		      IF(wqual(i,21).lt.ynchla(j)*alg_min(j)) THEN
!		         WRITE(*,*)'ERROR in Algae. Nmin < min',i,'Jday',jday,wqual(i,17)
		         wqual(i,21) = alg_min(j)*ynchla(j)
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
!      pause
	   RETURN
	   END SUBROUTINE ALGAE_TAHOE

!************************************************************************
	   SUBROUTINE alg_setl(lost, par)
!***********************************************************************
!     POP and PON settling are in Nutrient Settling subroutine
!     arfac = multiplicative factor to get area to m**2
!     chf(maxns,nchl) = algae cells per layer volume [cells/vol]
!     chfdv(maxns,nchl) = number of cells in a layer [cells/area]
!     dep_alg = depth of oft expanded algae layer
!     setldeptop(maxns,nchl) = distance of top of layer from bottom after settling
!     setldepbot(maxns,nchl) = distance of bottom of layer from bottom after settling
!     depth(maxns) = depth array for layers relative to bottom
!     depthm(maxns) = mean depth of a layer (m) i.e. between i and i-1
!     dep(maxns) = array for depth between layers (m)
!     fluxin(maxns,nchl) = flux of particles into layer
!     ns = number of dyresm layers
!     nchl = number of non-zero algae modeled
!     paro = light level at top of each layer [W/m2 PAR]
!     setldep(maxchl) = vert. settling depth in a time step [m]
!     phyto_setl_vel(maxchl) = settling velocities [m/day]
      
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	
      
      	
	   INTEGER*4 i,j,k,n,ii,l, mdays,time_step,max_layer(maxchl)
	   REAL*8 arfac
	   REAL*8 chf(maxns,maxchl)
	   REAL*8 chfdv(maxns,maxchl)
	   REAL*8 dep(maxns), ref_den
       REAL*8 par, max_setl_dep(200, maxchl),last_setl_dep(maxchl),last_paro(maxns)
	   REAL*8 setldeptop(maxns,maxchl)
	   REAL*8 setldepbot(maxns,maxchl)
	   REAL*8 phyto_fluxin(maxns,maxchl)
	   REAL*8 setldep(maxns,maxchl), restingdepth(maxns,maxchl)
       REAL(8) ::  denh2o(maxns),visco(maxns), denphyto(maxchl),diaphyto(maxchl)
	   REAL*8 lt_threshold, dep_alg(maxns,maxchl), delta_a,dummywqual(maxns,maxchl),phyto_setl_flux(maxns,maxchl),phyto_rise_flux(maxns,maxchl)
	   REAL*8 Before(maxchl),After(maxchl),Lost(maxchl),  max_dummywqual(maxchl), MAX_WQUAL(MAXCHL),sumdummywqual(maxchl)
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

! Initials

	   DO k = 1,nchl
	      DO i = ns,1,-1            
            setldeptop(i,k) = 0.0
		      setldepbot(i,k) = 0.0
		      setldep(i,k) = 0.0
	         chfdv(i,k) =0.0
	         phyto_fluxin(i,k)=0.0			! initials
	      ENDDO	! i layers
	   ENDDO	! k algae
!  ----------------------------------------------------
!     part 1: determining the light limitation on growth.
!     change the value of photosynthetically active radiation for
!     the surface to the surface value for an array of par

	   paro(ns) = par	! set in heatr
!     determine the attenuation coefficient for par from chla. note that
!     etwat could be replaced by a coefficient related to particle
!     concentrations at some stage
       mdays  = jday-(jday/1000)*1000
	   DO i = 1,ns	   
	      DO j = 1,nchl
		      et1(i) = et1(i)+etca(j)*wqual(i,j)
	      ENDDO ! j
	      DO j = 1,7
		      et1(i) = et1(i)+etpart(j)*cf(i,j)
	      ENDDO ! j
	      et1(i) = et1(i)+etwat    
         if(mdays.ge.90.and.mdays.lt.135) then
			   et1(i)=1.70d0   !eta=1.7/secchi depth (1 to 1.5 SD)
		  elseif(mdays.ge.135.and.mdays.lt.155) then
			   et1(i)=1.8d0   !eta=1.7/secchi depth (1 to 1.5 SD)
	     elseif(mdays.ge.155.and.mdays.lt.210) then    ! 160 to 210
			   et1(i)=1.9d0   !eta=1.7/secchi depth (1 to 1.5 SD)
		  elseif(mdays.ge.210.and.mdays.lt.300) then    ! 210 to 250
			   et1(i)=1.9d0   !eta=1.7/secchi depth (1 to 1.5 SD)     
	     elseif(mdays.ge.300.and.mdays.lt.330) then
			   et1(i)=1.8d0   !eta=1.7/secchi depth (1 to 1.5 SD)
        elseif(mdays.ge.330.and.mdays.lt.366) then
	         et1(i)=1.7d0   !eta=1.7/secchi depth (1 to 1.5 SD)
        endif

	   ENDDO   ! i
!  
!quim   Ted Swift clarity model for Tahoe
!  ******************LIGHT_TAHOE**********************************************
!      CALL Light_Tahoe(latit)
!  ****************************************************************************
!     calculate PAR at the surface of each layer
!  
!	   et1(i)=5.2
	   IF(ns.eq.1)GOTO 2
	   DO i = ns-1,1,-1
	      IF (paro(i+1).le.1d-20) THEN
		      paro(i) = 0.0
	      ELSE
		      paro(i) = paro(i+1)*(dexp(-et1(i)*(depth(i+1)-depth(i))))
          ENDIF       
	   ENDDO
2    CONTINUE

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
       ENDDO	! i
      
      If(ICLOCK.EQ.0) THEN
         time_step=1
      ELSE
         time_step=time_step+1
      ENDIF
      
      GOTO 5     
!***************************************************************************************************************** 
! CASE 1: APHANEZOMWNON AND MICROCYSTIS WILL STAY IN THE OPTIMAL LIGHT INTENSITY ZONE AND COME UP AS SUN GO DOWN
!***************************************************************************************************************** 
      DO J=1,NCHL
         DO I=NS,1,-1
            DUMMYWQUAL(I,J)=WQUAL(I,J)
         ENDDO
      ENDDO        
  ! CASE 1 Step 1:SETTLE DOWN TO MAXIMUM LIGHT INTENSITY ZONE    
      DO 3 j=4,4
         DO i = ns,1,-1
            IF(PARO(NS) .GE. LIGHT(J)) THEN              
               if(i.eq.ns.and.light(j).lt.paro(i)) then
                  max_setl_dep(time_step,j)=depth(ns)-depth(i)                  
                  wqual(i,j)=0.0
                  dummywqual(i+1,j)=0.0 !previous layer should be zero
               elseif(i.ne.ns.and.light(j).gt.paro(i).and.light(j).le.paro(i+1)) then
                  max_setl_dep(time_step,j)=depth(ns)-depth(i+1)
                  wqual(i,j)=wqual(i+1,j)+wqual(i,j)
                  GOTO 3
               elseif(i.ne.ns.and.light(j).lt.paro(i)) then
                  max_setl_dep(time_step,j)=depth(ns)-depth(i)
                  wqual(i,j)=wqual(i+1,j)+wqual(i,j)
                  wqual(i,j)=0.0                                                                                                                                                                                            
               endif 
             ENDIF  
         ENDDO !I           
3    CONTINUE	! J 

!CASE 1 Step 2: RISE AND EQUALLY DISTRIBUTE TO ALL LAYERS
! findout the layer of maximum accumulation

      IF(ICLOCK.GE.21600) THEN      
         DO j=4,4
            max_wqual(j) = 0.0            
            DO i = ns,1,-1
               IF(paro(ns).gt.10.0.and.paro(ns).lt.light(j)) THEN              
                  if(max_wqual(j).le.wqual(i,j)) then
                     max_wqual(j)=wqual(i,j) 
                     max_layer(j)=i                 
                  else
                     max_wqual(j)=max_wqual(j)
                  endif
               ENDIF 
              if (i.eq.ns.and.jday.eq.2004204.and.iclock.ge.21600.and.iclock.le.75600) write(996,4)
              if(             jday.eq.2004204.and.iclock.ge.21600.and.iclock.le.75600) write(996, fmt='(2i8, F15.6, 3I8,5f15.6)') &
              JDAY, ICLOCK,FLOAT(ICLOCK)/3600.0, i,ns,max_layer(j),depth(ns)-depth(i),max_wqual(j),dep(i),wqual(i,j),paro(i)               
            ENDDO !I           
         ENDDO	! J
      ENDIF  
         
 ! CASE 2 Step 3: Distribute the wqual equally in all layers above the maximum
       IF(ICLOCK.GE.21600) THEN      
         DO j=4,4                       
            DO i = ns,max_layer(j),-1
               wqual(i,j)=wqual(i,j)+ max_wqual(j)/FLOAT(max_layer(j))   
              if (i.eq.ns.and.jday.eq.2004204.and.iclock.ge.21600.and.iclock.le.75600) write(996,4)
              if             (jday.eq.2004204.and.iclock.ge.21600.and.iclock.le.75600) write(996, fmt='(2i8, F15.6, 3I8,5f15.6)') &
              JDAY, ICLOCK,FLOAT(ICLOCK)/3600.0, i,ns,max_layer(j),depth(ns)-depth(i),max_wqual(j),dep(i),wqual(i,j),paro(i)               
            ENDDO !I           
         ENDDO	! J
      ENDIF  
      
4   FORMAT(4X, 'JDAY',2X,'ICLOCK', 11X,'HOUR',7X,'I', 6X,'NS', X,'MAXLAYER', 9X, 'DEPTH',6X,'MAX_WQUAL',11X,'DEEP',10X,'WQUAL',12x, 'PAR')
    !      GOTO 19
5   CONTINUE
!*************************************************************************
!***************************************************************************************************************** 
! CASE 2: APHANEZOMWNON AND MICROCYSTIS WILL SETTLE DOWN AND RISE BASED ON DENSITY DIFFERENCE AND STOKE'S LAW
!***************************************************************************************************************** 
        DO j=1,nchl
            DO i = ns,1,-1              
             if(light(j).gt.paro(i).and.light(j).le.paro(i+1)) then
                max_setl_dep(time_step,j)=depth(ns)-depth(i+1)
                GOTO 6
              elseif(light(j).gt.paro(i)) then
                max_setl_dep(time_step,j)=depth(ns)-depth(ns)
              endif 
                dummywqual(i,j)=wqual(i,j)                
            ENDDO !I
         ENDDO
    
    6    CONTINUE	! J

		DO i = 1,ns
			denh2o(i) = den(i)+1000.0d0	
!     Absolute viscosity as a function of tmperature
			visco(i)=(-3.86704098283854D-10*(temp(i)**5.0)+			&
              1.28275099347771D-7*(temp(i)**4)-						&
              1.73356497821936D-5*(temp(i)**3)+						&
              1.27453807127464D-3*(temp(i)**2)-						&
              0.0587433984793861*temp(i)+							   &
              1.785149374874680)*0.001d0 
        ENDDO
        
      DO k=1,nchl
          if (iclock.eq.21600) ref_den=denh2o(ns)  !density at 6 AM
          if (k.ge.4.and.k.le.5) then              
              denphyto(k)= ref_den
          else
              denphyto(k)=1001.0d0
          endif
      ENDDO
      
      DO i =ns,1,-1
          do k=1,nchl
              denphyto(k)=0.0
              setldep(i,k)=0.0  !initials
              diaphyto(k)=20.00d-6
            IF (k.ge.4.and.k.le.5) THEN                  
                  if(paro(i).gt.10.0) then 
 ! the case phytoplankton settle down 
                      if(paro(ns).gt.light(k)) then
                          if(paro(i).gt.light(k)) then 
                              denphyto(k)=denh2o(ns)+50.0 
                          elseif(paro(i).gt.10.0.and.paro(i).lt.light(k)) then
                              denphyto(k)=denh2o(ns)-00.0
                          endif
 ! the case phytoplankton rise                         
                      elseif(paro(ns).le.light(k)) then
                          if(paro(i).gt.10.0) then 
                              denphyto(k)=denh2o(ns)-50.0
                          else
                              denphyto(k)=denh2o(ns)
                          endif                          
                      endif
 ! During night/early morning time
                  elseif(paro(i).le.10.0) then
                      denphyto(k)=denh2o(ns)                                            
                  endif  
                  
!  settling/rise velocity estimation using Stoke's law                                  
                  if (paro(i).gt.10.0) then
                     setldep(i,k)=dfloat(nosecs)*2.0d0*9.80665d0*((20.00d-6)**2.00d0)*(denphyto(k)-denh2o(ns))/(9.0d0*visco(i))
                  else
                     setldep(i,k)=0.0
                  endif
              ELSE
                  denphyto(k)=1005.0
                  setldep(i,k)=dfloat(nosecs)*2.0d0*9.80665d0*((20.00d-6)**2.00d0)*(denphyto(k)-denh2o(ns))/(9.0d0*visco(i))
              ENDIF              
          enddo          
      ENDDO

!************************************************************************  
	   DO k = 1,nchl  !1,nchl
	      DO i = ns,1,-1 
              phyto_setl_flux(i,k)=0.0         !Initial
              phyto_rise_flux(i,k)=0.0         !Initial
              restingdepth   (i,k)=0.0         !Initial
              chfdv          (i,k)=0.0         !Initial
              
              IF(i.eq.1) THEN
		         setldeptop(i,k) = depth(1)        !- setldep(i,k)
		         setldepbot(i,k) = depth(1) - setldep(i,k)
		         if(setldepbot(i,k).le.0.0) setldepbot(i,k)=0.0
	         ELSE
		         setldeptop(i,k) = depth(i-1)      !- setldep(i,k)
		         setldepbot(i,k) = depth(i-1) - setldep(i,k)
		         if(setldepbot(i,k).le.depth(i-1)-dep(i-1)) setldepbot(i,k)=depth(i-1)-dep(i-1)
	         ENDIF
!	         PRINT*,I,DEPTH(I),    setldeptop(i,k),    setldeptop(i,k)
!	         PAUSE
!             setldep(i,k)=setldeptop(i,k)-setldepbot(i,k)
 	         IF(ABS(setldep(i,k)).le.dep(i)) THEN
	            dep_alg(i,k)=setldep(i,k)
	         ELSEIF(ABS(setldep(i,k)).gt.dep(i)) THEN
	            dep_alg(i,k)=dep(i)*(setldep(i,k)/ABS(setldep(i,k)))  !change the sign positive to negative for the case of rise
             ENDIF
             
             IF(i.eq.ns.and.dep_alg(i,k).le.0.0) dep_alg(i,k)=0.0
! the case phytoplankon settle down                          
             IF(dep_alg(i,k).gt.0.0 .and. wqual(i,k).gt.alg_min(k)) THEN
                chfdv(i,k) = wqual(i,k)*dep_alg(i,k)/dep(i)  ! mg/m3
                if(chfdv(i,k).le.alg_min(k)) chfdv(i,k)=0.0
!                chfdv(i,k)=0.0                
                restingdepth(i,k)=setldepbot(i,k)
! the case pphytoplankton rise
             ELSEIF(dep_alg(i,k).lt.0.0.and.wqual(i,k).gt.alg_min(k))THEN
                 chfdv(i,k)=wqual(i,k)*dep_alg(i,k)/dep(i)                 
                if(abs(chfdv(i,k)).le.alg_min(k)) chfdv(i,k)=0.0
                if(i.eq.ns) then
                     chfdv(i,k)=0.0
                     restingdepth(i,k)=depth(ns)
                elseif(i.ne.ns) then
                    chfdv(i,k)=0.0
!                    chfdv(i,k)=wqual(i,k)*dep_alg(i,k)/dep(i)
                    restingdepth(i,k)=setldepbot(i,k)
                endif
             ENDIF 
 	         phyto_fluxin(i,k)=0.0			! initials
	      ENDDO	! i layers
       ENDDO	! k algae
       
! Initials
 
      DO k =1,nchl                    
          do i=ns,1-1              
               phyto_setl_flux(i,k)=0.0
               phyto_rise_flux(i,k)=0.0 
          enddo          
      ENDDO

       DO k=1,nchl
           DO i=ns,1,-1
! the case phytoplankton settle down
              IF (i.eq.ns.and.dep_alg(i,k).gt.0.0) THEN
                  dummywqual(i,k)=dummywqual(i,k) - abs(chfdv(i,k))
                  DO j=i-1,1,-1
                     if(depth(j).gt.restingdepth(i,k).and.depth(j-1).le.restingdepth(i,k)) then !smaller than depth but greater than bottom means next layer depth
                        phyto_setl_flux(j,k)=abs(chfdv(i,k))
                        GOTO 11 
                     endif
                  ENDDO
 11               CONTINUE            
              ELSEIF(i.ne.ns.and.restingdepth(i,k).gt.depth(1).and.dep_alg(i,k).gt.0.0) THEN
                  dummywqual(i,k)=dummywqual(i,k) - abs(chfdv(i,k))
                  DO j=i-1,1,-1
                     if(depth(j).gt.restingdepth(i,k).and.depth(j-1).le.restingdepth(i,k)) then
                        phyto_setl_flux(j,k)=abs(chfdv(i,k))
                        GOTO 12 
                     endif
                  ENDDO
12                CONTINUE
!         Permanent settle down to the sediment
              ELSEIF(i.ne.ns.and.restingdepth(i,k).le.depth(1).and.dep_alg(i,k).gt.0.0) THEN
!                   dummywqual(i,k)=dummywqual(i,k) - abs(chfdv(i,k))
!                  phyto_setl_flux(1,k)=abs(chfdv(i,k))
              ENDIF
!the case phytoplankton rise
! Because the light intensity is low it determines only few upper layers
!So better consider deeper layers and equally distribute
!****************************
!              IF    (i.eq.ns.and.dep_alg(i,k).lt.0.0) THEN                  
!                  dummywqual(i,k)=dummywqual(i,k)-0.0
!              ELSEIF(i.ne.ns.and.dep_alg(i,k).lt.0.0) THEN
!                  dummywqual(i,k)=dummywqual(i,k)-abs(chfdv(i,k))
!                  DO j=i,1,-1
!                     if(depth(j).lt.restingdepth(i,k).and.depth(j+1).ge.restingdepth(i,k)) then
!                        phyto_rise_flux(j,k)=abs(chfdv(i,k))
!                        GOTO 13 
!                     endif
!                  ENDDO
! 13               CONTINUE       
!              ENDIF
!*********************** 
	      ENDDO	! i layers
	   ENDDO	! k algae
	    
 
       DO k=1,nchl
           DO i=ns,1,-1
              dummywqual(i,k)=wqual(i,k)-chfdv(i,k) + phyto_setl_flux(i,k) + phyto_rise_flux(i,k) 
              IF(dummywqual(i,k).LE.0.0) dummywqual(i,k)=0.0
	      ENDDO	! i layers
	   ENDDO	! k algae
!***********************************************
!   Findout the maximum phytoplankon concentration
!  
      DO k=1,nchl
          max_layer(k)=ns       !Initials
         IF(dep_alg(ns,k).lt.0.0) THEN   
            max_dummywqual(k) = 0.0            
            DO  j= ns,1,-1               
                if(dummywqual(j,k).ge.max_dummywqual(k)) then
                     max_dummywqual(k)=dummywqual(j,k)    
                     max_layer(k)=j-2              
                endif
            ENDDO !j
            
            IF(max_layer(k).LE.1)    max_layer(k)=1
            IF(max_layer(k).GE.ns)   max_layer(k)=ns
! CASE 2 Step 3: add all and distribute the wqual equally in all layers above the maximum          
               sumdummywqual(k)=0.0        
            DO j = ns,max_layer(k),-1
                   sumdummywqual(k)=sumdummywqual(k)+dummywqual(j,k)
            ENDDO
 !           PRINT*,max_layer(k),sumdummywqual(k)
 !           PAUSE
            DO i=ns,max_layer(k),-1     
                 dummywqual(i,k)=sumdummywqual(k)/DFLOAT(ns-max_layer(k)+1)
            ENDDO
        ENDIF                        
      ENDDO  !k
!*************************************************

 	   DO k=1,nchl
	      DO i=ns,1,-1	   
              if (i.eq.ns.and.jday.eq.2004230.and.iclock.ge.0.and.iclock.le.86400) write(997,26)
              if(             jday.eq.2004230.and.iclock.ge.0.and.iclock.le.86400) write(997, fmt='(i8, F15.6, I8,11f15.3, i8)') &
              JDAY, FLOAT(ICLOCK)/3600.0, I,depth(i),restingdepth(i,k),setldepbot(i,k),setldeptop(i,k)-setldepbot(i,k),dep(i),dep_alg(i,k), chfdv(i,k),wqual(i,k),dummywqual(i,k), paro(i),temp(i),max_layer(k)	   
	      ENDDO
	   ENDDO
26  FORMAT(4X, 'JDAY',11X,'HOUR',7X,'I', 10X, 'DEPTH',3X,'RESTINGDEPTH',5X, 'SETLDEPBOT',8X, 'TOP-BOT',11X,'DEEP',8X,'DEP_ALG',10X,'CHFDV',10X,'WQUAL',5x,'dummywqual',12x, 'PAR',11x,'temp', 6x,'Max_layer',6X,'XX')
  
    
	   DO k = 1,nchl     
	      DO i = 1,ns
              IF(dummywqual(i,k).LE.0.0) dummywqual(i,k)=0.0
              wqual(i,k)=dummywqual(i,k)
	      ENDDO	! i
	   ENDDO	! k    
    
    DO k = 1, nchl  
	      Before(k)   = 0.0
	      After(k)    = 0.0
	      Lost(k)     = 0.0
      ENDDO    ! k

 	   DO k = 1,nchl
         DO i = 1,ns
	         Before(k)   = Before(k)   + wqual(i,k)   *vol(i)    
        ENDDO  ! i
      ENDDO    ! k

19 CONTINUE

!  
!   Estimate the lost of mass by setling
!  
	   DO k = 1,nchl
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

! 19 CONTINUE
  
	   DO i = ns,1,-1
	      DO j = 1,nchl
	         IF(wqual(i,j) .lt. alg_min(j)) THEN
!	            WRITE(*,*)'ERROR in Setling algae. Algae conc. less than min',i,jday
		         wqual(i,j) = alg_min(j)
	            ! pause
             ENDIF           
             if (j.eq.4.and.(depth(ns)-depth(i)).le.5.0.and.mod(iclock,1800).eq.0) then
                if (i.eq.ns.and.jday.eq.2004204.and.iclock.ge.21600.and.iclock.le.75600) write(998,25) !fmt='(///)')
                if(             jday.eq.2004204.and.iclock.ge.21600.and.iclock.le.75600) write(9999, fmt='(2i8, 10f15.1)') &
                jday,iclock,float(iclock)/3600.0,depth(ns)-depth(i),paro(i),temp(i),den(i)+1000.0,denphyto(4), &
                setldep(i,4),wqual(i,j), phyto_fluxin(i,j), chfdv(i,j) 
                 
            endif           
		   ENDDO ! j chl
       ENDDO	 ! i layers
25     FORMAT (4X,'JDAY',4X,'TIME',11X,'HOUR',10X,'DEPTH',12X,'PAR',11X,'TEMP',4X,'H2O_DENSITY', 8X,'PHYTO_DEN',1X, 'SETLING_DEPTH',2X,'BLUEGREEN_ALG',11X, 'FLUX')       
       DO j = 1,nchl
           last_setl_dep(j)=max_setl_dep(time_step,j)
       ENDDO
            last_paro(ns)=paro(ns)  
	   RETURN
	   END SUBROUTINE alg_setl
!**********************************************************************
	   SUBROUTINE alg_sed
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
		      Do j=1,nchl
!		         IF(sedalgae(j).le.alg_min(j))sedalgae(j) = alg_min(j)
!		         wqual(i,j) = wqual(i,j) + sedalgae(j) * Afactor
		      ENDDO
! set BOD flux rate at 0.5 mg/m2/day...	
!		      wqual(i,9) = wqual(i,9) + 0.5 * Afactor
 	      ENDIF
	   ENDDO
	   RETURN
	   END SUBROUTINE alg_sed
!************************************************************************************
