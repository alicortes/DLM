!************************************************************************
		SUBROUTINE PARTICLES(Lost)
!**************************************************************************
!  arfac = multiplicative factor to get area to m**2
!  cf(maxns,7) = particle array
!  cff(maxns,sizall) = particle and water quality array
!  cffdv(maxns,sizall) = number of particles in a layer
!  coag = variable for determining particle coagulation
!  denh2o(maxns) = density of water (kg m**-3)
!  denscf(7) = density of 7 particle groups (kg m**-3)
!  densy(1) = density of total phosphorus and total nitrogen (kg m**-3)
!  densy(2) = density of bod (kg m**-3)
!  densy(3) = density of fe (kg m**-3)
!  densy(4) = density of mn (kg m**-3)
!  deptop(maxns,sizall) = distance of top of layer from bottom after settling
!  depbotom(maxns,sizall) = distance of bottom of layer from bottom after settling
!  d(sizall) = array for diameter of particles and phytoplankton
!  dd(size) = array for diameter of particles
!  depth(maxns) = depth array for layers relative to bottom
!  depthm(maxns) = mean depth of a layer (m) i.e. between i and i-1
!  dep(maxns) = array for depth between layers (m)
!  fluxin(maxns,sizall) = flux of particles into layer
!  hh(maxns),s(maxns),ge(maxns),g(maxns)
!  h1,hsig,diss,di
!  mu,ma
!  ns = number of dyresm layers
!  na(size),nu(size),pa(size)
!  ar(size),ngo(size),tot(size)
!  size = 7 particle groups
!  sizall = 10 settling arrays
!  setldep(maxns,sizall) = settling velocities
!  v(sizall) = array for volume of particles
!  visco(maxns) = absolute viscosity of water
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE
		REAL(8),	PARAMETER :: kboltz =1.3806503d-23			! Boltzmann constant = 1.3806503 × 10^(-23) m2 kg s-2 K-1
		INTEGER*4,	PARAMETER :: sizall = 10, size = 7
		INTEGER*4 :: i,j,k,l,n,f,ii,tt,kk,ttest
		INTEGER*4 :: index1,delta_time,TEP(maxns)
!IF	INTEGER*4 :: jday

		REAL(8) ::  ar(size),ngo(size),fu
		REAL(8) ::  Before(size),After(size),Lost(size)
		REAL(8) ::  cff(maxns,sizall)
		REAL(8) ::  cffdv(maxns,sizall)
		REAL(8) ::  coagu
		REAL(8) ::  d(size),dd(size),dn(size),v(size), test(size)
		REAL(8) ::  Delta_area
		REAL(8) ::  denh2o(maxns)
		REAL(8) ::  dep(maxns)
		REAL(8) ::  deptop(maxns,size)
		REAL(8) ::  depbotom(maxns,size)
		REAL(8) ::  di
		REAL(8) ::  diamin, diamax
		REAL(8) ::  fluxin(maxns,sizall)
		REAL(8) ::  hh(maxns),s(maxns,size),ge(maxns),g(maxns)
		REAL(8) ::  mu(size),ma, mu_prev
		REAL(8) ::  na(size),nu(size),pa(size)
		REAL(8), PARAMETER :: pi=3.1415926535897932384626433832795d0
		REAL(8), PARAMETER ::volfac = 1.0d+3,arfac=1.0D+6
		REAL(8) ::  setldep(maxns,size)
		REAL(8) ::  tot(size),const
		REAL(8) ::  visco(maxns),d_new(size), cha_death(maxns,size)
!gbs---------------------------------------------------------------
      REAL*8 SOD_area (maxns), SOD_radius (maxns), SUM_SOD_area
		REAL*8 SOD_radius_L (maxns),SOD_radius_U (maxns)
		REAL*8 slant_height(maxns)
!gbs-----------------------------------------------------------------
		IF(partsw) THEN
			index1=1
		ELSE
			GOTO 336
		ENDIF
	

!  depth of each layer
		DO i = ns,2,-1
			dep(i) = depth(i)-depth(i-1)
		ENDDO
		dep(1) = depth(1)
!	open(unit=9999, file='test.txt', status='unknown')
! 
!  Units are number of particles/m3 at every box
! 
		DO i = 1,ns
		   DO J =1,7
		      if(cf(i,j).le.1.0)    cf(i,j)=1.0d0
		 !     if(cf(i,1).ge.3.5d10) cf(i,1)=3.5d10
		 !     if(cf(i,2).ge.3.5d9) cf(i,2)=3.5d9
		 !     if(cf(i,3).ge.6.5d8) cf(i,3)=6.5d8
		 !     if(cf(i,4).ge.2.5d8) cf(i,4)=2.5d8
		 !     if(cf(i,5).ge.6.0d7) cf(i,5)=6.0d7
			   cff(i,j) = cf(i,j) 
			ENDDO			
		ENDDO          
!  
		DO i = 1,size
			d(1)= 0.50d-6	   !less than 1. micrometer !00.60d-6 
			d(2)= 1.00d-6	   !1 to 2					    !01.500d-6
			d(3)= 2.00d-6	   !2 to 4					    !1.8000d-6
			d(4)= 4.00D-6	   !4 to 8					    !5.0000D-6
			d(5)= 8.00D-6	   !8 to 16					    !8.800D-6
			d(6)= 16.00D-6	   !16 to 32				    !20.000D-6	   !16.0
			d(7)= 32.00D-6	   !32 to 63				    !63.000D-6	   !32.0
			dd(i) = d(i)
			dn(1)= 0.80d-6	   !less than 1. micrometer !00.80d-6 
			dn(2)= 1.6d-6	   !1 to 2					    !01.600d-6
			dn(3)= 3.2d-6	   !2 to 4					    !3.2000d-6
			dn(4)= 6.4D-6	   !4 to 8					    !6.4000D-6
			dn(5)= 12.8D-6	   !8 to 16					    !12.800D-6
			dn(6)= 25.6D-6	   !16 to 32				    !25.600D-6	   !16.0
			dn(7)= 48.00D-6	!32 to 63				    !48.000D-6	   !32.0
		ENDDO
!gbs new ****************************************************************
		DO i = ns,2,-1		
			SOD_radius (i)= DSQRT((0.5*(area(i)+area(i-1))*arfac)/PI)
			SOD_radius_L(i)=DSQRT(area(i-1)*arfac/PI)
			SOD_radius_U(i)=DSQRT(area(i)*arfac/PI)
			slant_height(i)=SQRT(dep(i)*dep(i)+(SOD_radius_U(i)-SOD_radius_L(i))**2)
			SOD_area   (i)= 2.0D0*PI*SOD_radius (i)*slant_height(i)			
		ENDDO
!gbs new*****************************************************************

!  assign particle diameters to the remaining arrays
!  determine volumes (micrometres**3) assuming a spherical shape
!		DO i=size+1,sizall
!			d(i) = 25.0
!			v(i) = pi*d(i)**3.0/6.0
!		ENDDO
! 
!    Water properties...Kg/m-s
! 			
		DO i = 1,ns
			denh2o(i) = den(i)+1000.0d0	
!			visco(i) = 10.0d0**(1301.0d0/(998.333d0 +					&
!     			8.1855d0*(temp(i)-20.0d0) +							&
!     			0.00585d0*(temp(i)-20.0d0)**2.0d0)-3.30233d0)
!			visco(i)=(-3.86704098283854D-10*(temp(i)**5.0)+1.2827509935D-7*		&
!     		(temp(i)**4)-1.73356497822D-5*(temp(i)**3)+1.274538071275D-3*		&
!     		(temp(i)**2)-0.0587433984793861*temp(i)+									&
!     		1.785149374874680)*0.001d0

!     Absolute viscosity as a function of tmperature
			visco(i)=(-3.86704098283854D-10*(temp(i)**5.0)+			&
              1.28275099347771D-7*(temp(i)**4)-						&
              1.73356497821936D-5*(temp(i)**3)+						&
              1.27453807127464D-3*(temp(i)**2)-						&
              0.0587433984793861*temp(i)+								&
              1.785149374874680)*0.001d0
              
 !			visco(i)=(-3.86704098283854D-10*(temp(ns)**5.0)+			&
 !             1.28275099347771D-7*(temp(ns)**4)-						&
 !             1.73356497821936D-5*(temp(ns)**3)+						&
 !             1.27453807127464D-3*(temp(ns)**2)-						&
 !             0.0587433984793861*temp(ns)+								&
 !             1.785149374874680)*0.001d0
		ENDDO
		
!------------------------------------------------------
!gbs******COAGULATION STUFF FOR 7 SIZES PARTICLES******
!------------------------------------------------------
!		GOTO 334
		delta_time=nosecs
		DO 1000 tt=1, (int(nosecs/delta_time)) 	
			DO i = 1,ns
				DO k = 1,size
					cf(i,k) = cff(i,k)
				ENDDO
			ENDDO

			DO n = 1,size-1
				na(n) = (1.0d0/900000.0d0)*((10.0d0)**(n-1))
				pa(n) = 1.0d0-na(n)
				nu(size-n) = ((10.0d0**n)+1.0d0)/(10.0d0**n)
			ENDDO
			DO n = 1,ns
				DO j = 1,size
					hh(n)= (2.0d0/3.0d0)*kboltz*(temp(n)+273.0d0)*delta_time/visco(n)
					s(n,j)= pi*(denscf(j) - denh2o(n))*9.8d0*delta_time/(72.0*visco(n))

					IF((depth(ns)-depth(n)).lt.h1)THEN
						di = diss
					ELSEIF((-((depth(ns)-h1-depth(n))/hsig)**2).le.-50.0) THEN !gbs: to avoid the numerical instability 'underflow'
						di = diss*exp(-50.00)
					ELSE	              
						di = diss*exp(-((depth(ns)-h1-depth(n))/hsig)**2)
					ENDIF
					ge(n) = (di/(visco(n)/denh2o(n)))**(1.0d0/2.0d0)
					g(n)  = ge(n)*delta_time		
!					IF(jday.eq.240.and.j.eq.1) THEN
!						WRITE(105, fmt='(i8,2f10.2,3f20.1)')		&
!     					jday,depth(n),temp(n),hh(n),s(n,j),g(n)
!						ENDIF
				ENDDO	   	! j
			ENDDO		! n
!print*, 'ns',depth(ns), temp(ns)
!print*, '01',depth(1), temp(1)

!  compute the particles leaving any size, ngo(n)
			DO f = 1,ns
				DO n = 1,size
					ar(n)  = 0.0d0
					ngo(n) = 0.0d0
					IF(cff(f,n).le.0.0) cff(f,n)=1.0     
				ENDDO

				DO ii = 1,size            ! size-1
					DO n = 1,size
						coag=0.0d0
						ar(n)=0.0d0
						mu(n)=0.0d0
						IF(ii.eq.n)THEN
							fu = 0.5d0          ! counted two times  1 1,2 2,3 3,4 4,5 5, 6 6,7 7
						ELSE
							 fu=0.5d0           ! counted one times   12 21, 13 31, 14 41 etc       
						ENDIF                  
!gbs2010	           WRITE(*,fmt='(2I6,f10.5)')ii,n,coag
!gbs2010	           PASUE
                  IF(wqual(f,1).le.1.0d0) THEN
                      TEP(f)=1.0d0 !*wqual(f,1)
                  ELSE
                     TEP(f)=1.0*wqual(f,1)
                  ENDIF 
                  IF (depth(f).le.depth(ns)-50.00) THEN
                     if(wqual(f,1).le.1.0d0) then
                        TEP(f)=1.0d0 !*wqual(f,1)
                     else
                        TEP(f)=1.0*wqual(f,1)
                     endif  						   
 						ENDIF
 						IF(ii.eq.1.or.n.eq.1) THEN       !5 and 5
							coag=2.0d-7*TEP(f)*DEXP(-(2.0d6/((1.0d0/dd(ii))+(1.0d0/dd(n)))))*	&	!1.-7d)
          						(DLOG10(cff(f,ii))**3.00)*(DLOG10(cff(f,n))**3.00)        
 !							coag=2.0d-6*TEP(f)*DEXP(-(2.0d6/((1.0d0/dd(ii))+(1.0d0/dd(n)))))*	&	!1.-7d)
 !         						(DLOG10(cff(f,ii))**2.00)*(DLOG10(cff(f,n))**2.00)
 !         						print*,ii,n,coag
 !         						pause
							IF(n.eq.1) THEN
								IF(cff(f, n+1).gt.0.90d0*cff(f, n)) coag=coag/10.0d0			   
							ELSE
								IF(cff(f, n).gt.0.90d0*cff(f, n-1)) coag=coag/10.0d0
							ENDIF
							IF(ii.eq.1) THEN 
								IF(cff(f,ii+1).gt.0.90d0*cff(f,ii)) coag=coag/10.0d0	   
							ELSE 
								IF(cff(f,ii).gt.0.90d0*cff(f,ii-1)) coag=coag/10.0d0
							ENDIF
              
!							IF (SD.lt.17.0d0) coag=coag/2.0d0

						ELSEIF(ii.eq.2.or.n.eq.2) THEN			          
       					coag=7.0d-6*TEP(f)*DEXP(-(2.0d6/((1.0d0/dd(ii))+(1.0d0/dd(n)))))*	&	!2.0d-7)
          						(DLOG10(cff(f,ii))**3.00)*(DLOG10(cff(f,n))**3.00) 
							IF(n.eq.2) THEN
								IF(cff(f, n+1).gt.0.90d0*cff(f, n)) coag=coag/10.0d0			   
							ELSE
								IF(cff(f, n).gt.0.90d0*cff(f, n-1)) coag=coag/10.0d0
							ENDIF
							IF(ii.eq.2) THEN 
								IF(cff(f,ii+1).gt.0.90d0*cff(f,ii)) coag=coag/10.0d0	   
							ELSE 
								IF(cff(f,ii).gt.0.90d0*cff(f,ii-1)) coag=coag/10.0d0
							ENDIF
              
!							IF (SD.lt.17.0d0) coag=coag/2.0d0
								
						ELSEIF(ii.eq.3.or.n.eq.3) THEN			          
       					coag=5.0d-5*TEP(f)*DEXP(-(2.0d6/((1.0d0/dd(ii))+(1.0d0/dd(n)))))*	&	!1.5d-7)
          						(DLOG10(cff(f,ii))**3.00)*(DLOG10(cff(f,n))**3.00) 
							IF(n.eq.3) THEN
								IF(cff(f, n+1).gt.0.90d0*cff(f, n)) coag=coag/10.0d0			   
							ELSE
								IF(cff(f, n).gt.0.90d0*cff(f, n-1)) coag=coag/10.0d0
							ENDIF
							IF(ii.eq.3) THEN 
								IF(cff(f,ii+1).gt.0.90d0*cff(f,ii)) coag=coag/10.0d0	   
							ELSE 
								IF(cff(f,ii).gt.0.90d0*cff(f,ii-1)) coag=coag/10.0d0
							ENDIF
              
!							IF (SD.lt.17.0d0) coag=coag/2.0d0
						ELSEIF(ii.eq.4.or.n.eq.4) THEN			          
       					coag=10.0d-4*TEP(f)*DEXP(-(2.0d6/((1.0d0/dd(ii))+(1.0d0/dd(n)))))*	&	!1.0d-8)
          						(DLOG10(cff(f,ii))**3.00)*(DLOG10(cff(f,n))**3.00)                         
							IF(n.eq.4) THEN
								IF(cff(f, n+1).gt.0.90d0*cff(f, n)) coag=coag/10.0d0			   
							ELSE
								IF(cff(f, n).gt.0.90d0*cff(f, n-1)) coag=coag/10.0d0
							ENDIF
							IF(ii.eq.4) THEN 
								IF(cff(f,ii+1).gt.0.90d0*cff(f,ii)) coag=coag/10.0d0	   
							ELSE 
								IF(cff(f,ii).gt.0.90d0*cff(f,ii-1)) coag=coag/10.0d0
							ENDIF
              
!							IF (SD.lt.17.0d0) coag=coag/2.0d0
						ELSEIF(ii.eq.5.or.n.eq.5) THEN			          
       					coag=3.0d0*TEP(f)*DEXP(-(2.0d6/((1.0d0/dd(ii))+(1.0d0/dd(n)))))*	&	!1.0d-10)
          						(DLOG10(cff(f,ii))**3.00)*(DLOG10(cff(f,n))**3.00)        
 !         						print*,ii,n,coag
 !         						pause
							IF(n.eq.5) THEN
								IF(cff(f, n+1).gt.0.90d0*cff(f, n)) coag=coag/10.0d0			   
							ELSE
								IF(cff(f, n).gt.0.90d0*cff(f, n-1)) coag=coag/10.0d0
							ENDIF
							IF(ii.eq.5) THEN 
								IF(cff(f,ii+1).gt.0.90d0*cff(f,ii)) coag=coag/10.0d0	   
							ELSE 
								IF(cff(f,ii).gt.0.90d0*cff(f,ii-1)) coag=coag/10.0d0
							ENDIF              
!							IF (SD.lt.17.0d0) coag=coag/2.0d0
!                    IF (SD.lt.15.0d0) coag=coag*100.0d0	
						ELSEIF(ii.eq.6.and.n.eq.6) THEN			          
       					coag=3.0d+6*TEP(f)*DEXP(-(2.0d6/((1.0d0/dd(ii))+(1.0d0/dd(n)))))*	&	!1.0d-11)
          						(DLOG10(cff(f,ii))**3.00)*(DLOG10(cff(f,n))**3.00)        
!         			          
							IF (SD.lt.12.0d0) coag=coag*150.0d0
						ELSEIF(ii.eq.7.or.n.eq.7) THEN			          
       					coag=0.0d-12*TEP(f)*DEXP(-(2.0d6/((1.0d0/dd(ii))+(1.0d0/dd(n)))))*	&	!1.0d-11)
          						(DLOG10(cff(f,ii))**3.00)*(DLOG10(cff(f,n))**3.00)   
          						
						ENDIF
!						IF(ii.le.5.and.n.le.5) THEN
!						   coag=0.010
!						   IF (SD.lt.15.0d0) coag=coag/10.0d0
!						ELSE
!						   coag=0.0
!						ENDIF


!						print*,coag
!						WRITE(*, fmt='(3i6,2x,f8.5)')jday,ii,n,coag
!						PAUsE
						mu(n)= -fu*coag*coagu(dd,f,hh,g,s,ii,n)
						d_new(n)=(dn(ii)**3.0d0+dn(n)**3.0d0)**(1.0d0/3.0d0)
!						print*,ii,n,coag,coagu(dd,f,hh,g,s,ii,n), mu(n)
!						pause
!						WRITE(*,fmt='(2i6,2x,3f12.9)')ii,n,dd(ii),dd(n),d_new(n)

!						IF(d_new(n).gt.dd(7))THEN         !6
!							mu(n)=0.0d0      !restricts for unnecessary big partilces which doesn't occur or coagulation is zero		
!						ELSE				 
!							mu(n)=mu(n)
!						ENDIF
!						ngo(ii) = ngo(ii)+mu(n)
!--------------------particles leaving from nth size--------------
						IF(ii.eq.n)THEN
							cff(f,n) = cff(f,n)+mu(n)  
							ngo(ii) = ngo(ii)+mu(n)

						ELSE
							cff(f,n) = cff(f,n)+mu(n)
							ngo(ii) = ngo(ii)+mu(n)	
						ENDIF						
!						cff(f,n) = cff(f,n)+mu(n)
!       Estimation of coagulated particles in Euclidean method
!                	print*, ii,n,d_new(n), mu(n)
!						PASUE
!					ENDDO !n
!----------------search where the new particles are fitted in -------------				
!--------------------particles leaving from nth size and entering into kk size --------------
!	            DO k = 1, size-1
!	              IF(d_new(n).gt.dd(k).and.d_new(n).le.dd(k+1))THEN                          
!	          		 ar(k)=(-mu(n))
!                      cff(f,k) = cff(f,k)+ ar(k)
!	               ENDIF
!	            ENDDO !k
	            
						IF(d_new(n).gt.0.50d-6.and.d_new(n).le.1.0d-6)THEN                          
	          			ar(1)=(-mu(n))
							cff(f,1) = cff(f,1)+ ar(1)
						ELSEIF(d_new(n).gt.1.0d-6.and.d_new(n).le.2.0d-6)THEN
	          			ar(2)=(-mu(n))
							cff(f,2) = cff(f,2)+ ar(2)
						ELSEIF(d_new(n).gt.2.0d-6.and.d_new(n).le.4.0d-6)THEN
	          			ar(3)=(-mu(n))
							cff(f,3) = cff(f,3)+ ar(3)
						ELSEIF(d_new(n).gt.4.0d-6.and.d_new(n).le.8.0d-6)THEN
	          			ar(4)=(-mu(n))
							cff(f,4) = cff(f,4)+ ar(4)
						ELSEIF(d_new(n).gt.8.0d-6.and.d_new(n).le.16.0d-6)THEN
	          			ar(5)=(-mu(n))
							cff(f,5) = cff(f,5)+ ar(5)
						ELSEIF(d_new(n).gt.16.0d-6.and.d_new(n).le.32.0d-6)THEN
	          			ar(6)=(-mu(n))
							cff(f,6) = cff(f,6)+ ar(6)
						ELSEIF(d_new(n).ge.32.0d-6)THEN
	          			ar(7)=(-mu(n))
							cff(f,7) = cff(f,7)+ ar(7)	           
						ENDIF 			  					           
					ENDDO !n			
				ENDDO	! ii 
	
				DO n = 1,size
!					tot(n) = ar(n)+mu(n)     
!					cff(f,n) = cf(f,n)+tot(n) 
               IF(ISNAN(cff(f,n))) cff(f,n)=1.0d0   			
					IF(cff(f,n).le.0.0) cff(f,n)=1.0d0
!xxxxxxxxxxxx Since partcle of size 7 is ineffective, take the value =100.0 xxxxxxxxxxxxxxx
					IF(cff(f,7).ge.1.0d5) cff(f,7)=1.0d5
				ENDDO
			ENDDO	! f
334   CONTINUE    
!--------------------------------------------------------------
!gbs********** END OF COAGULATION STUFF of 7 sizes particles***
!gbs-----------Settling of particles---------------------------
!  settling - setldep is the distance settled (settling velocity*time),
!  a negative value indicates a buoyant particle

			DO i = 1,ns 	  
				DO k = 1,size 
					IF (k.eq.1) THEN
						test(k)=dd(k)/2.0d0
						setldep(i,k)=1.0d0*delta_time*2.0d0*((test(k))**2.00d0)*		&     
     						9.80665d0*(1.0d0*denscf(k)-denh2o(i))/(9.0d0*visco(i))	            
					ELSEIF(k.eq.2) THEN
						test(k)=dd(k)/2.0d0
						setldep(i,k)=1.0d0*delta_time*2.0d0*((test(k))**2.0d0)*		&   
     					9.80665d0*(1.0d0*denscf(k)-denh2o(i))/(9.0d0*visco(i)) 		 	
					ELSEIF(k.eq.3) THEN
						test(k)=dd(k)/2.0d0
						setldep(i,k)=1.0d0*delta_time*2.0d0*((test(k))**2.00d0)*		&  
     					9.80665d0*(1.0d0*denscf(k)-denh2o(i))/(9.0d0*visco(i)) 
					ELSEIF(k.eq.4) THEN
						test(k)=dd(k)/2.0d0
						setldep(i,k)=1.0d0*delta_time*2.0d0*((test(k))**2.0d0)*		&     
     					9.80665d0*(1.0d0*denscf(k)+0.0d0*densy(1)-denh2o(i))/(9.0d0*visco(i)) 
					ELSEIF(k.eq.5) THEN
						test(k)=dd(k)/2.0d0
						setldep(i,k)=1.0d0*delta_time*2.0d0*((test(k))**2.00d0)*		&  
     					9.80665d0*(1.0d0*denscf(k)+0.0d0*densy(1)-denh2o(i))/(9.0d0*visco(i)) 
					ELSEIF(k.eq.6) THEN
						test(k)=dd(k)/2.0d0
						setldep(i,k)=1.0d0*delta_time*2.0d0*((test(k))**2.0d0)*		&	  
     					9.80665d0*(1.0d0*denscf(k)+0.0d0*densy(1)-denh2o(i))/(9.0d0*visco(i))
  					ELSEIF(k.eq.7) THEN
						test(k)=dd(k)/2.0d0
						setldep(i,k)=1.0d0*delta_time*2.0d0*((test(k))**2.0d0)*		&	  
     					9.80665d0*(1.0d0*denscf(k)+0.0d0*densy(1)-denh2o(i))/(9.0d0*visco(i))
!						print*,setldep(i,k)
					ENDIF
				ENDDO	! k
!				IF(jday.eq.240) THEN
!					WRITE(105, fmt='(i8,2f10.2,f20.10)')jday,depth(i),temp(i),setldep(i,1)
!				ENDIF
			ENDDO	! i	
!  movement of particles within a layer or out of a layer
! 
! Calculate the initial mass of WQ variables
! 
			DO k = 1,size
				Before(k)   = 0.0d0
				After(k)    = 0.0d0
				Lost(k)     = 0.0d0
			ENDDO    ! k

 			DO k = 1,size
				DO i = 1,ns
					Before(k)  = Before(k) + cff(i,k)*vol(i)   
				ENDDO   ! i
			ENDDO   ! k 	  
  
			DO k = 1,size
				DO i = ns,1,-1
					IF(i.eq.1) THEN
						cffdv(i,k)=setldep(i,k)*cff(i,k)*arfac*area(i) 
						deptop(i,k)=depth(i)-setldep(i,k)
!					IF(k.eq.7.and.i.eq.1) print*,deptop(i,k), setldep(i,k),depth(i)
!					PASUE
						IF(deptop(i,k).le.0.0) THEN
							deptop(i,k)=0.0 		!bottom layer means permanent settling
							cffdv(i,k) = 0.0		!bottom layer means permanent settling
						ENDIF
						depbotom(i,k)=0.0	      
					ELSE
						deptop(i,k)=depth(i)-setldep(i,k)
						depbotom(i,k)=depth(i-1)-setldep(i,k)	
						IF(deptop(i,k).le.depth(1))THEN
							cffdv(i,k) =0.0       !bottom layer means permanent settling
							deptop(i,k)=depth(1)	!bottom layer means permanent settling
							depbotom(i,k)=0.0
						ELSE
							IF(setldep(i,k).lt.dep(i)) THEN
								cffdv(i,k) = 1.0*cff(i,k)*(setldep(i,k)*arfac*		&
												0.5*(area(i)+area(i-1)))	!# particles               		
							ELSEIF(setldep(i,k).ge.dep(i)) THEN
								cffdv(i,k) = 1.0*cff(i,k)*(dep(i)*arfac*				&
     										0.5*(area(i)+area(i-1)))	!# particles 
							ENDIF
						ENDIF
					ENDIF
		  			fluxin(i,k) = 0.0d0
				ENDDO	! i
			ENDDO	! k
!gbs convert cff back to a particle concentration
			DO 10 k = 1,size
				DO 20 i = ns,1,-1
					ttest=0
					DO 30 j = i,1,-1     !ith layer particles must go to the layer below.So start searching after i
!gbs case where all settling particles fall within top and bottom of layer 1:

						IF(j.eq.1)THEN  
!gbs case where settling particles fall within layer 1 and extend further than the bottom:
							IF(deptop(i,k).le.depth(j).and.deptop(i,k).gt.0.0.and.depbotom(i,k).le.0.0)	&
     							fluxin(j,k) = fluxin(j,k)+ cffdv(i,k)

!gbs case where settling particles extend above the top of layer 1 and extend further than the bottom:
							IF(deptop(i,k).gt.depth(j).and.depbotom(i,k).le.0.0)			&
     	    					fluxin(j,k) = fluxin(j,k)+ cffdv(i,k)

!gbs case where all settling particles fall within layer 1:
							IF(deptop(i,k).le.depth(j).and.depbotom(i,k).gt.0.0)			&
     							fluxin(j,k) = fluxin(j,k)+cffdv(i,k)

!gbs case where settling particles extend above layer 1 and fall within
!gbs layer 1
							IF(deptop(i,k).gt.depth(j).and.depbotom(i,k).lt.			&
     							depth(j).and.depbotom(i,k).gt.0.0)							&
     	   					fluxin(j,k)=fluxin(j,k)+cffdv(i,k)
!IF(k.eq.7.and.i.eq.1) print*, fluxin(j,k) 
						ELSE
!gbs case where all settling particles fall within the top and bottom of a layer
							IF(deptop(i,k).le.depth(j).and.depbotom(i,k).ge.depth(j-1)) THEN
      						fluxin(j,k) = fluxin(j,k)+cffdv(i,k) 
								ttest=ttest+1  
								GOTO 20
!gbs case where settling particles extend above the top of a layer and
!gbs below the bottom of the layer
							ELSEIF(deptop(i,k).ge.depth(j).and.depbotom(i,k).lt. depth(j).and.	&
                                depbotom(i,k).ge. depth(j-1)) THEN
      						fluxin(j,k) = fluxin(j,k)+cffdv(i,k) 
								ttest=ttest+1  
								GOTO 20
!gbs   case where the settling dep above the layer, but bottom is below the layer and
!gbs    bottom is above the below bottom layer

							ELSEIF (deptop(i,k).ge.depth(j).and.depbotom(i,k).lt. depth(j).and.	&
                                depbotom(i,k).le. depth(j-1)) THEN
									fluxin(j,k) = fluxin(j,k)+cffdv(i,k)
									ttest=ttest+1
									GOTO 20
!case where it reaches the bottom layer
							ELSEIF (deptop(i,k).ge.depth(j).and.depbotom(i,k).lt. depth(j).and.	&
                                depbotom(i,k).le. depth(1)) THEN
									fluxin(1,k) = fluxin(1,k)+cffdv(i,k)
									ttest=ttest+1
									GOTO 20
							ENDIF
						ENDIF
  30				CONTINUE ! j
  20			CONTINUE	! i
  10 		CONTINUE ! k  
			DO k = 1,size
				DO i = ns,1,-1
					IF(i.eq.1) THEN
						cff(i,k) = cff(i,k)+(fluxin(i,k)-cffdv(i,k))/					&
     							(1.0d0*dep(i)*area(i)*arfac)  !#particles/m^3 
					ELSE
						cff(i,k) = cff(i,k)+(fluxin(i,k)-cffdv(i,k))/					&
                     (1.0d0*dep(i)*arfac*0.5*(area(i)+area(i-1)))  !#particles/m^3
					ENDIF
					IF(cff(i,k).le.0.0) cff(i,k)=1.0					               	
				ENDDO
			ENDDO ! k
			DO i = 1,ns
!gbs boundary condition: when the particles hits bottoms layers,
!gbs they settled down permanently. ELSE there is upwelling of the particles during winter mixing.
            DO k=1,7
	 			   IF(i.LE.2) THEN	   
					   cf(i,k) = 1.0  
					   cff(i,k) = 1.0     
				   ELSEIF (i.gt.1.and.depth(i).lt.(depth(ns)-20.0d0))THEN
				      cff(i,k) = cff(i,k)-SOD_area(i)*0.005d0*cff(i,k)/(0.5*dep(i)*(area(i)+area(i-1))*arfac) 				      
		            IF(cff(i,k).le.0.0) cff(i,k)=1.0
		            cf(i,k) = cff(i,k)
				   ELSE
 					   cf(i,k) = cff(i,k)				
	   		   ENDIF
	   		   IF(ISNAN(cf(i,k))) cf(i,k)=1.0D0 
	   		ENDDO	 !k
			ENDDO	! i
1000  ENDDO !tt (sub_time step loop
!	print*, cff(1,7) !*vol(200), fluxin(200,1)
! 
! Nutrients Budget
! 
 		DO k = 1,size
			DO i = 1,ns
				After(k)  = After(k) + cff(i,k)*vol(i)    
			ENDDO   ! i
		ENDDO   ! k 

        
		DO k = 1,size
			Lost(k) = Before(k) - After(k)
		ENDDO	 
 
336	CONTINUE
		RETURN
		END SUBROUTINE PARTICLES
!***********************************************************
		REAL*8 FUNCTION coagu(dd,f,hh,g,s,ii,n)
!***********************************************************
      USE DLMWQ_VARIABLES      
		IMPLICIT NONE
      SAVE
		INTEGER*4 ii,n,f
		REAL*8 dd(7),hh(5000),g(5000),s(5000) 
			coagu=(((hh(f)*((dd(ii)+dd(n))**2))/(dd(ii)*dd(n)))+(((dd(ii)+		&
			dd(n))**3)*g(f)/6.0)+(s(f)*((dd(ii)+dd(n))**3)*abs(dd(ii)-dd(n))))*cf(f,ii)*cf(f,n)
	!	print*,((hh(f)*((dd(ii)+dd(n))**2))/(dd(ii)*dd(n))),(((dd(ii)+		&
	!		dd(n))**3)*g(f)/6.0),(s(f)*((dd(ii)+dd(n))**3)*abs(dd(ii)-dd(n))),coagu
	!		pause
		RETURN
		END FUNCTION coagu
