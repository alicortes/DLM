!*************************************************************************     
      SUBROUTINE Light_Tahoe(latit)
!***********************************************************************!
!
!     This routine calculates the light extinction coefficient for the 
!     Lake Tahoe, development of the Tahoe light model of Swift (2000) 
! 
!    scatter = scattering contribution of sediments (clay)
!    a_cdom  = absortion of colored dissolved organic matter
!    a_star  = absortivity coefficient of chlorophyll-material
!    b_star  = scattering coefficient of chlorophyll-material
!    gamma   = constant on equation SD
!    SD      = Secchi Depth
!    mu_zero = cosinus dependence
!    Kd      = diffuse irradiance attenuation coeffient
!    PAR     = Photosynthetically Active Radiation
!***********************************************************************!	
      USE DLMWQ_VARIABLES                    
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE	

		INTEGER(4), PARAMETER:: size = 7
		INTEGER(4):: i,j,j_counter,k,l,month
!IF		INTEGER(4):: jday
		INTEGER(4):: exits

		REAL(8):: Kd
		REAL(8):: A_1(maxns),B_1(maxns),C_1(maxns),BSed_1(maxns)
		REAL(8):: A_mean, B_mean, C_mean, Depth_mean,Bsed_mean,Depth_S
		REAL(8):: Chloro_mean,Parts_mean(size),Volumetric
		REAL(8):: Organic_mean,Total_mean
      REAL(8):: Kd_mean
		REAL(8):: b_sed(maxns)
		REAL(8):: Scatter(size)
		REAL(8):: SD_prima
		REAL(8):: SD_measured(size)
		REAL(8):: POP_mean,PON_mean
		REAL(8):: Vol_mean
      REAL(8):: Nitrate_mean
      REAL(8):: Ammonia_mean
      REAL(8):: THP_mean
		REAL(8):: Total, Organic
		REAL(8):: Dec, radian, tau, elevation, mu_zero, latit
		REAL(8), PARAMETER :: pi=3.1415926535897932384626433832795d0
		REAL(8), PARAMETER :: arfac  = 1.0d+6,volfac = 1000.0d0
		REAL(8):: Ch_abs,Ch_sca,Wa_abs,Wa_sca,CDOMab,Pa_sca,Pa_sc1,Pa_sc2, Pa_sc3,Pa_sc4,Pa_sc5,Pa_sc6 
!
!		New Set of Parameters (06/09/2002)
!
		REAL(8),PARAMETER:: 		gamma    = 8.7000D0   			! # 1	! constant, Gamma
		REAL(8),PARAMETER:: 		a_water  = 0.0120D0				! # 5	! /m, Absroption Pure Water
		REAL(8),PARAMETER:: 		a_cdom   = 0.0380D0				! # 2	! /m, Absorption CDOM
		REAL(8),PARAMETER:: 		a_star   = 0.0250D0				! # 4	! m2/mg, Absroption Chla !low value = up the line, high = low the line
		REAL(8),PARAMETER:: 		b_water  = 0.0027D0				! # 5	! /m, Scattering Pure Water
		REAL(8),PARAMETER:: 		b_star   = 0.1050D0				! # 6 ! m2/mg^0.62, Sacttering Chla

		Scatter(1)  = 4.287D-12			! # 7-13 ! m2/particle, Particle Scattering Coefficients
		Scatter(2)  = 3.015D-11       ! m2/particle
		Scatter(3)  = 9.939D-11       ! m2/particle,
		Scatter(4)  = 3.757D-10       ! m2/particle,
		Scatter(5)  = 1.459D-09       ! m2/particle,
		Scatter(6)  = 5.831D-09       ! m2/particle,
		Scatter(7)  = 0.0d0           ! m2/particle,
4004  CONTINUE
!
!    Factor (Ino/Total) particles
!
		fraction  = 0.0D0
		Organic   = 0.0D0
		Total     = 0.0D0

		DO i = 1,size
			Total   = Total   + cf(ns,i)
		ENDDO ! i
		DO i = 1,nchl
			 Organic = Organic + DBLE(phyto_part_fac(i)*wqual(ns,i)) 
		ENDDO ! i

!
		fraction = Total /(Total+ Organic)
		IF (fraction.ge.1.0D0) fraction = 1.0D0
		IF (fraction.lt.0.0D0) fraction = 0.0D0
!     
!
!    Inicialize  variables
!
!
! Solar angle
!
		radian	 = 2.0D0 * pi/360.0D0									! Degrees to radians factor
		jday      = jday - int(jday/1000)*1000							! Jday (XXX) this omitted the year like 1999002 
		Dec       = 23.45D0 * dsin((284+jday)*2.0D0*pi/365.0D0)	! Solar Declination
		tau       = DBLE(360.0D0*12.0D0/24.0D0)						! Hour. Fix to midday
		elevation = dsin(dsin(radian*latit)*dsin(radian*dec)	&	! Solar Elevation
		  		  + dcos(radian*latit)*dcos(radian*dec)*dcos(radian*tau))
!gbs	mu_zero   = dcos(dasin(dcos(elevation)/1.33D0))				! Angle mu
		mu_zero   = dcos(dasin(dsin(elevation)/1.33D0))				! Angle mu
		exits=0
!xxxxxxxxxxxxxx LOOP Starts from here xxxxxxxxxxxxxxxxxxxxxxxx  
5		CONTINUE
  	
		Total_mean    = 0.0d0
		Organic_mean  = 0.0d0
		fraction_mean = 0.0d0
		DO i = 1,ns
			A_1(i)    = 0.0D0
			B_1(i)    = 0.0D0
			b_sed(i)  = 0.0D0
			BSed_1(i) = 0.0D0
			Kd        = 0.0D0
			DO j = 1,nchl
				A_1(i)    = DBLE(A_1(i)   + DBLE(a_star*wqual(i,j)))
				B_1(i)    = DBLE(B_1(i)   + DBLE(b_star*wqual(i,j)**0.62D0))			
			ENDDO ! j
      
			DO j = 1,7
				b_sed(i)  = DBLE(b_sed(i) +DBLE(DBLE(fraction*Scatter(j)*cf(i,j)))) 
			ENDDO  !j  
			BSed_1(i) = b_sed(i)	
			A_1(i)    = DBLE(A_1(i) + a_water + a_cdom)	  	
			B_1(i)    = DBLE(B_1(i) + b_water + b_sed(i))	
			C_1(i)    = DBLE(A_1(i) + B_1(i))	
			Kd = (A_1(i)/mu_zero) * DSQRT(1.0D0+(0.425D0*mu_zero-0.19D0)* DBLE(B_1(i)/A_1(i)))   	
			et1(i) = 1.00*Kd 
		ENDDO ! i
!
!   Estimate Secchi depth. Use Iterative Process
!
      j_counter = 1
		DO 5000 i = ns,1,-1
			A_mean       = 0.0D0
			B_mean       = 0.0D0
			Bsed_mean    = 0.0D0
			C_mean       = 0.0D0
			Chloro_mean  = 0.0D0
			Depth_mean   = 0.0D0
			Kd_mean      = 0.0D0
			POP_mean     = 0.0D0
			PON_mean     = 0.0D0
			SD_prima     = 0.0D0
			Depth_S      = 0.0D0
			Nitrate_mean = 0.0D0
			Ammonia_mean = 0.0D0
			THP_mean     = 0.0D0
			Vol_mean     = 0.0D0
			Total_mean	  = 0.0D0
			Organic_mean = 0.0d0
			fraction_mean = 0.0d0	

			DO l = 1,7
				Parts_mean(l) = 0.0D0
			ENDDO
         
			DO j = ns,ns - j_counter,-1
				IF(j-1.eq.0.or.ns-j_counter.eq.0) THEN
					Volumetric   = vol(1)*volfac
					Depth_S      = depth(1)
					POP_mean     = POP_mean     + wqual(1,17)*Volumetric 
					PON_mean     = PON_mean     + wqual(1,23)*Volumetric	   
					Nitrate_mean = Nitrate_mean + wqual(1,21)*Volumetric
					Ammonia_mean = Ammonia_mean + wqual(1,22)*Volumetric
					THP_mean     = THP_mean     + wqual(1,16)*Volumetric
					DO k = 1,nchl
						Chloro_mean = Chloro_mean  + wqual(1,k)*Volumetric 
					ENDDO    ! k
					DO k = 1,7
						Parts_mean(k) = Parts_mean(k) + cf(1,k)*Volumetric
					ENDDO    ! k	   
					A_mean =     A_mean     + A_1(1)*depth_S   ! 
					B_mean =     B_mean     + B_1(1)*depth_S   ! 
					Bsed_mean  = Bsed_mean  + BSed_1(1)*depth_S 	   
					Depth_mean = Depth_mean + Depth_S 
					Vol_mean   = Vol_mean   + Volumetric
					GOTO 3456
				ENDIF
				Volumetric   = vol(j)*volfac
				Depth_S      = depth(j)-depth(j-1)	     
				POP_mean     = POP_mean     + wqual(j,17)*Volumetric 
				PON_mean     = PON_mean     + wqual(j,23)*Volumetric	   
				Nitrate_mean = Nitrate_mean + wqual(j,21)*Volumetric
				Ammonia_mean = Ammonia_mean + wqual(j,22)*Volumetric
				THP_mean     = THP_mean     + wqual(j,16)*Volumetric
				DO k = 1,nchl
					Chloro_mean = Chloro_mean  + wqual(j,k)*Volumetric
				ENDDO    ! k
				DO k = 1,7
					Parts_mean(k) = Parts_mean(k) + cf(j,k)*Volumetric 
				ENDDO    ! k
				A_mean     = A_mean     + A_1(j)   *depth_S !Unit of Length 
				B_mean     = B_mean     + B_1(j)   *depth_S !Unit of Length
				Bsed_mean  = Bsed_mean  + BSed_1(j)*depth_S !Unit of Length	   
				Depth_mean = Depth_mean + Depth_S 
				Vol_mean   = Vol_mean   + Volumetric
			ENDDO  ! j

 3456		CONTINUE     	  


			A_mean      = DBLE(A_mean/Depth_mean)
			B_mean      = DBLE(B_mean/Depth_mean)
			Bsed_mean   = DBLE(Bsed_mean/Depth_mean)	 	 
			C_mean      = A_mean + B_mean
			Kd_mean     = A_mean * (DSQRT(1+(0.425*mu_zero-0.19D0)*		&
							  DBLE(DBLE(B_mean/A_mean))/mu_zero))

			SD_prima    = DBLE(gamma/DBLE(C_mean + Kd_mean))
			Chloro_mean = DBLE(Chloro_mean/Vol_mean)

			DO k = 1,7 
				Parts_mean(k)  = DBLE(Parts_mean(k)/Vol_mean)
			ENDDO       ! k
			Nitrate_mean = DBLE(Nitrate_mean/Vol_mean)
			Ammonia_mean = DBLE(Ammonia_mean/Vol_mean)
			THP_mean     = DBLE(THP_mean/Vol_mean)
			POP_mean     = DBLE(POP_mean/Vol_mean)
			PON_mean     = DBLE(PON_mean/Vol_mean)
	  
			IF(SD_prima.gt.(depth(ns)-depth(i))) THEN
				SD = Depth_mean
				j_counter = j_counter + 1
!				WRITE(*,fmt='(3f10.3)')SD_prima,(depth(ns)-depth(i))
			ELSE
				SD           = SD_prima
				A_Tahoe      = A_mean
				B_Tahoe      = B_mean
				Kd_Tahoe     = Kd_mean
				Bsed_Tahoe   = Bsed_mean
				Chloro_Tahoe = Chloro_mean	  
				DO k = 1,7
					Parts_Tahoe(k)  = Parts_mean(k)
				ENDDO
				POP_Tahoe     = POP_mean
				PON_Tahoe     = PON_mean
				Nitrate_Tahoe = Nitrate_mean
				Ammonia_Tahoe = Ammonia_mean
				THP_Tahoe     = THP_mean
				DO k = 1,7
					Total_mean   = Total_mean   + Parts_mean(k)
				ENDDO ! k 	  
				Organic_mean = Organic_mean + DBLE(phyto_part_fac(1)*Chloro_mean)	
				fraction_mean = Total_mean /(Total_mean+ Organic_mean)
		      IF (fraction_mean.ge.1.0D0) fraction_mean = 1.0D0        !gbs
      		IF (fraction_mean.lt.0.0D0) fraction_mean = 0.0D0	     !gbs
!gbs   Becuase the repeation takes place for some values and it never converge
!gbs   the value has changed. ELSE GOTO 5 repeats and never converge.		  
!gbs	    IF(DABS(fraction-fraction_mean).gt.0.01) THEN
!gbs       WRITE(*,fmt='(2f10.3)')fraction-fraction_mean,SD
  				IF(DABS(fraction-fraction_mean).gt.0.01) THEN
					fraction = fraction_mean	  
					exits=exits+1
					IF(exits.gt.5) GOTO 4900
					GOTO 5
				ELSEIF (DABS(fraction-fraction_mean).LE.0.01) THEN 
				   GOTO 4900 			
				ENDIF	   
				GOTO 4900
			ENDIF 
5000	ENDDO  ! i	      
4900  CONTINUE

!
! Contribution from each source
!
      Ch_abs = 0.0d0
      Ch_sca = 0.0d0
      Wa_abs = 0.0d0
      Wa_sca = 0.0d0
      CDOMab = 0.0d0
      Pa_sca = 0.0d0
      Pa_sc1 = 0.0d0
      Pa_sc2 = 0.0d0
      Pa_sc3 = 0.0d0
      Pa_sc4 = 0.0d0
      Pa_sc5 = 0.0d0
      Pa_sc6 = 0.0d0
		DO k = ns-1, i,-1
			DO j = 1,nchl
				Ch_abs    = Ch_abs   + DBLE(a_star*wqual(k,j))*(depth(k+1)-depth(k))
				Ch_sca    = Ch_sca   + DBLE(b_star*wqual(k,j)**0.62D0)*(depth(k+1)-depth(k))			
			ENDDO ! j
      
			DO j = 1,7
				   Pa_sca = Pa_sca +DBLE(Scatter(j)*cf(k,j))*(depth(k+1)-depth(k))
				 IF (J.EQ.1) THEN
				   Pa_sc1 = Pa_sc1 +DBLE(Scatter(j)*cf(k,j))*(depth(k+1)-depth(k))
				 ELSEIF (J.EQ.2) THEN
				   Pa_sc2 = Pa_sc2 +DBLE(Scatter(j)*cf(k,j))*(depth(k+1)-depth(k))				 
				 ELSEIF (J.EQ.3) THEN
				   Pa_sc3 = Pa_sc3 +DBLE(Scatter(j)*cf(k,j))*(depth(k+1)-depth(k))
				 ELSEIF (J.EQ.4) THEN
				   Pa_sc4 = Pa_sc4 +DBLE(Scatter(j)*cf(k,j))*(depth(k+1)-depth(k))
				 ELSEIF (J.EQ.5) THEN
				   Pa_sc5 = Pa_sc5 +DBLE(Scatter(j)*cf(k,j))*(depth(k+1)-depth(k))
				 ELSEIF (J.EQ.6) THEN
				   Pa_sc6 = Pa_sc6 +DBLE(Scatter(j)*cf(k,j))*(depth(k+1)-depth(k))				 
				 ENDIF
			ENDDO  !j				
			Wa_abs    =  Wa_abs + a_water*(depth(k+1)-depth(k))
			Wa_sca    =  Wa_sca + b_water*(depth(k+1)-depth(k))
			CDOMab    =  CDOMab + a_cdom*(depth(k+1)-depth(k))					
    ENDDO ! k

      IF(ISNAN(SD).or.SD<=10.0d0)SD=10.0d0
      WRITE (99,5005)JDAY,ICLOCK,SD, depth(NS)-DEPTH(I),Ch_abs,Ch_sca,Wa_abs,  &
      Wa_sca,CDOMab,Pa_sca,Pa_sc1,Pa_sc2, Pa_sc3,Pa_sc4,Pa_sc5,Pa_sc6
5005  FORMAT(2I8, 2F10.3, 12F10.6)
		RETURN
		END SUBROUTINE Light_Tahoe