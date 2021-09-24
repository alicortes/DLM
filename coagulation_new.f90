		SUBROUTINE COAGULATION (Lost) 
		USE DFLIB		
		USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE
!**********************IMPORTANT******************************
! ALL UNITS ARE IN METERS, GRAMS, HOURS,KELVIN  (hopefully)  *
!*************************************************************

! eps = energy dissipation rate - Jackson 2001, default value = 0.01
! kinvis = kinematic viscosity - Jackson 2001, default value = 0.01, temo. dependant, pressure dependant - needs to be taken from DLM
! pi = ration of circumfrence to diameter
! gr = gravity
! dyvis = dynamic viscosity, temp. dependant, pressure dependant - needs to be taken from DLM
! kboltz = Boltzman's constant = 1.380658 * 10 ^ -23
! temp = temperature - depends on depth and season.  Needs to be taken from DLM
! beta = overall probability that 2 particles will collide
! betash = probability that 2 particles will collide due to shear in the water
! betabr = probability that 2 particles will collide due to Brownian motion
! betads = probability that 2 particles will collide due to differential sedimentation
! coag = probability that if 2 particles collide they will stick
! change = the number of particles that form from the collision of 2 particles
! fraclen = fractal length = r ^ D = one of the conserved dimensions ( D = fractal dimension, r = radius ) - Jackson 2001
! dd = radius of particle
! ns=number of layers
! mass = mass of particle
! a = matrix that keeps track of the number of particles in each 2-D bin
! b = matrix that keeps track of excess fractal length - numbers are excess per particle
! c = matrix that keeps treack of excess mass - numbers are excess per particle
! totalinorg = total number of inorganic particles
! totalorg = total number of organic particles
! orgden = organic particles density
! inorgden = inorganic particles density
! aa = total particle 

!		REAL(8), PARAMETER :: eps = 46656.0D0, kinvis = 0.0054684D0, pi = 3.14159265358979D0, gr = 127008000.0D0, dyvis =  5648.4D0, kboltz = 1.789332768D-13, temp =  277.d0
	
!		REAL(8) :: eps , kinvis , gr, dyvis, kboltz, temp, Delta_T
		REAL(8) :: eps(maxns), kinvis(maxns), gr, dyvis(maxns), kboltz,temp1,Delta_T
		REAL(8), PARAMETER :: pi = 3.14159265358979D0
		INTEGER*4, PARAMETER :: size = 8, sizall=11
		INTEGER*4, PARAMETER :: OUT = 1000
!		INTEGER, PARAMETER :: ns = 200
!		REAL(8) :: beta, betash, betabr, betads, change, coag, rand1, 
		REAL(8) :: beta, betash, betabr, betads, change, rand1,		&
       rand2, masssum, fraclensum, masssumorg, fraclensumorg,		&
       totalinorg, totalorg,  time_begin, time_end,			      &
        newlength, newparticles, count3, count4, count5, sink,		&
       masssumorg1, fraclensumorg1
		INTEGER*4:: k, L, i,ii,j, m1, n1,count,count1,count2, i1,jj1, i2,		&
        j2, time, m, i3, j3, i4, j4, ln, L5, i7, j7, time1, bbb, plus,		&
        plus1, k8, L8, ccc, time7, i9,jj,index1, kkk, LLL
!  
		real(8) ::  denh2o(maxns), orgden, inorgden(size)
		real(8) ::  dvisco(maxns)
		real(8) ::  kvisco(maxns),TEP(maxns)
		real(8) ::  Before(sizall),After(sizall),Lost(sizall)
		real(8) ::  ediss(maxns), dep(maxns),fluxin(maxns,sizall,sizall),			&
               aadv(maxns,sizall,sizall),depleft(maxns,sizall,sizall),			&
               deplefb(maxns,sizall,sizall), b_prime(maxns,sizall,sizall),		&
               c_prime(maxns,sizall,sizall), denspart(sizall)
		real(8) ::  cff(maxns,sizall),setldep(maxns,sizall,sizall)

		real(8) ::coef_a(maxns,sizall,sizall),coef_b(maxns,sizall,sizall)
		real(8), parameter ::volfac = 1.0d+3,arfac=1.0D+6
		REAL(8), DIMENSION(size) :: fraclen, newfraclen,dd,ddn, mass, newmass,fracdim, add
		REAL(8), DIMENSION(maxns,size,size) ::  b, c, d1, bb1, c1, d2, aa  !(layer, mass, length)
!gbs---------------------------------------------------------------
      REAL*8 SOD_area (maxns), SOD_radius (maxns), SUM_SOD_area
		REAL*8 SOD_radius_L (maxns),SOD_radius_U (maxns)
		REAL*8 slant_height(maxns)
		REAL*8 ri1, ri2, rj1, rj2,	mi, mj, wi, wj, pp, Di, Dj 
!gbs-----------------------------------------------------------------
!CALL CPU_TIME(time_begin)
!fraclen = (/0.125D-18, 1.0D-15, 5.65685D-15, 1.6D-11, 6.4D-11, 2.33032218D-9, 8.114653145D-9, 999999999999.0D0/)!1.25d-13, 2.50047d-13, 5.00566d-13, 99999999.0D0/)! 1.d-12, 1.3997d-9, 3.1072d-9, 6.893d-9, 1.53d-8, 2.2463d-6, 3.01576d-6, 4.0485d-6, 5.43567D-6, 1.477d-4, 2.32d-4, 3.1d-4, 999999999.D0/) 
		dd = (/0.8d-6, 1.60d-6, 3.20d-6, 6.40D-6,			&
            12.80D-6, 25.60D-6, 48.0D-6, 99.0D0/) !32.0D-6, 999999999999.0D0/) ! m, 1.d-4, 1.414d-4, 2.d-4, 2.828d-4, 4.d-4, 4.757d-4, 5.657d-4, 6.727d-4, 8.d-4, 11.31d-4, 16.d-4, 20.d-4, 9999999999.D0/)
		ddn= (/0.5d-6, 1.00d-6, 2.0d-6, 4.0D-6,			&
            8.0D-6, 16.00D-6, 32.0D-6, 63.0D0/) !32.0D-6, 999999999999.0D0/) ! m, 1.d-4, 1.414d-4      
!		fracdim=(/3.0D0,3.0D0,3.0D0,2.7D0,2.4D0,2.2D0,2.0D0,100.0D0/)
		fracdim=(/3.0D0,3.0D0,3.0D0,3.0D0,3.0D0,3.0D0,3.0D0,100.0D0/)
!      fracdim=(/3.0D0,3.0D0,3.0D0,3.0D0,3.0D0,3.0D0,3.0D0,100.0D0/)

!	mass = (/9D-9, 1.2D-8, 1.5D-8,D0/)!4.0D-6, 5.0D-6, 6.0D-6, 7.0D-6, 8.0D-6, 9.0D-6, 10.0D-6, 11.0D-6, 12.0D-6, 13.0D-6, 14.0D-6, 15.0D-6, 9999999.0D0/) !how do we get the mass?
!	add = (/4.13D03, 4.13D03,	4.13D03, 2.72D02,	2.72D02, 6.18D01,	6.18D01, 0D0/)
! 	inorgden = 1.0D0 * 2.65D6 !density of inorganic matter
		orgden =1025.0d0 !1.0D0 * 1.025D6 !density of organic matter
!		DO i=1, size
!			IF(i.eq.size)THEN
!				inorgden(i) = 2650.0d0
!			ELSE
!				inorgden(i) = 2650.0d0
!			ENDIF
!		ENDDO
inorgden(1) = 2650.0d0
inorgden(2) = 2650.0d0
inorgden(3) = 2650.0d0
inorgden(4) = 2650.0d0
inorgden(5) = 2650.0d0
inorgden(6) = 2650.0d0
inorgden(7) = 2650.0d0
inorgden(8) = 2650.0d0

!totalinorg = 0.1D11
!totalorg = 0.9D11
!
! Calculate the initial mass of WQ variables
!
!  Units are number of particles/m3 at every box
 		DO i = ns,2,-1		
			SOD_radius (i)= DSQRT((0.5*(area(i)+area(i-1))*arfac)/PI)
			SOD_radius_L(i)=DSQRT(area(i-1)*arfac/PI)
			SOD_radius_U(i)=DSQRT(area(i)*arfac/PI)
			slant_height(i)=SQRT(dep(i)*dep(i)+(SOD_radius_U(i)-SOD_radius_L(i))**2)
			SOD_area   (i)= 2.0D0*PI*SOD_radius (i)*slant_height(i)			
		ENDDO 
		DO i = 1,ns
		   DO J =1,7
		      if(cf(i,j).le.1.0)cf(i,j)=1.0d0
		!      if(cf(i,2).ge.3.5d9) cf(i,2)=3.5d9
		!      if(cf(i,3).ge.6.0d8) cf(i,3)=6.0d8
		!      if(cf(i,4).ge.2.5d8) cf(i,4)=2.5d8
		!      if(cf(i,5).ge.6.0d7) cf(i,5)=6.0d7
			   cff(i,j) = cf(i,j)			    
			ENDDO			
		ENDDO          
	   DO k = 1,sizall
			Before(k)   = 0.0d0
			After(k)    = 0.0d0
			Lost(k)     = 0.0d0
      ENDDO    ! k

 		DO k = 1,7
			DO i = 1,ns
				Before(k)  = Before(k) + cff(i,k)*vol(i)   
			ENDDO   ! i
		ENDDO   ! k 

 
!   properties of water 
		DO i = 1,ns
			denh2o(i) = den(i)+1000.0d0						   
!     Absolute viscosity as a function of tmperature
			dvisco(i)=(-3.86704098283854D-10*(temp(i)**5.0)+		&
              1.28275099347771D-7*(temp(i)**4)-						&
              1.73356497821936D-5*(temp(i)**3)+						&
              1.27453807127464D-3*(temp(i)**2)-						&
              0.0587433984793861*temp(i)+								&
              1.785149374874680)*1.0D-3
!          print*,'gbs',temp(i),dvisco(i)    
!          dvisco(i) = 10.0d0**(1301.0d0/(998.333d0 +	&	
!     			8.1855d0*(temp(i)-20.0d0) +0.00585d0*(temp(i)-20.0d0)**2.0d0)-3.30233d0)	!Kg/m-s
!     			print*,'gbs',temp(i),dvisco(i)
!     			pause
			kvisco(i)=dvisco(i)/denh2o(i)								!m2/s
			IF((depth(ns)-depth(i)).lt.h1)THEN
				ediss(i) = diss
			ELSEIF(hsig.le.0.01) THEN       !gbs: to avoid the numerical instability 'underflow'
	        ediss(i)= diss
			ELSEIF(((-(depth(ns)-h1-depth(i))**2)/hsig).le.-50.0) THEN !gbs: to avoid the numerical instability 'underflow'
             ediss(i) = diss*exp(-50.00)
			ELSE	              
				ediss(i)= diss*exp((-(depth(ns)-h1-depth(i))**2)/hsig)
			ENDIF		
		ENDDO
! depth of each layer
		DO i = ns,2,-1
			dep(i) = depth(i)-depth(i-1)
		ENDDO
		dep(1) = depth(1)

!	Delta_T = 3600 !time step, in seconds
		Delta_T = nosecs
		
! Definition of Variables
		gr = 9.8D0 * ( Delta_T ) ** 2
		kboltz = 1.380648813D-23 * ( Delta_T ) ** 2 !kboltz=1.346d-23 J/K, 1 J= 1 kg m2/s2

      DO i = 1,ns
			eps(i)=ediss(i)*( Delta_T ) ** 3  
			kinvis(i)= kvisco(i) * ( Delta_T )
			dyvis(i) = dvisco(i) * ( Delta_T )
      ENDDO

		DO i = 1, size 
		   newfraclen(i)=(0.5d0*ddn(i)) ** fracdim(i)
		   fraclen(i)=(0.5d0*dd(i)) ** fracdim(i)
			IF ( i.le.3) THEN					
	 			mass(i) = (1.00D0 * inorgden(i)+0.0D0 * orgden)			&
							 *(4.0D0/3.0D0)*pi*(dd(i)/2.0D0)**3.0D0 !assigning mass and length to bins
	 			newmass(i) = (1.00D0 * inorgden(i)+0.0D0 * orgden)			&
							 *(4.0D0/3.0D0)*pi*(ddn(i)/2.0D0)**3.0D0 !assigning mass and length to bins		
			ELSEIF ( i.gt.3.and.i.le.4) THEN 		
				mass(i) = (0.60D0 * inorgden(i) + 0.40D0 * orgden)	&
            *(4.0D0/3.0D0)*pi*(dd(i)/2.0D0)**3.0D0 !assigning mass and length to bins
	 			newmass(i) = (0.30D0 * inorgden(i) + 0.70D0 * orgden)			&
							 *(4.0D0/3.0D0)*pi*(ddn(i)/2.0D0)**3.0D0 !assigning mass and length to bins
			ELSEIF (i.gt.4.AND.i.le.size-1) THEN					
				mass(i) = (0.30D0 * inorgden(i) + 0.70D0 * orgden)		&
            *(4.0D0/3.0D0)*pi*(dd(i)/2.0D0)**3.0D0 !assigning mass and length to bins
	 			newmass(i) = (0.30D0 * inorgden(i) + 0.70D0 * orgden)			&
							 *(4.0D0/3.0D0)*pi*(ddn(i)/2.0D0)**3.0D0 !assigning mass and length to bins	
			ELSEIF ( i.eq.size) THEN		
			  fraclen(i) = 99999.0
			  mass(i) = 99999.0
			  newmass(i)=9999.0
			ENDIF
		ENDDO
		count1 = 0
		count2 = 0
		count3 = 0	
		IF (partsw) THEN
		   DO i = 1,ns
			   DO j=1,size
				   DO k=1,size					
						IF(k.eq.j) THEN
							IF(j.eq.1.and.k.eq.1) THEN
								aa(i,j,k) = cff(i,1)
							ELSEIF(j.eq.2.and.k.eq.2) THEN
								aa(i,j,k) = cff(i,2)
							ELSEIF(j.eq.3.and.k.eq.3) THEN
								aa(i,j,k) = cff(i,3)
							ELSEIF(j.eq.4.and.k.eq.4) THEN
								aa(i,j,k) = cff(i,4)
							ELSEIF(j.eq.5.and.k.eq.5) THEN
								aa(i,j,k) = cff(i,5)
							ELSEIF(j.eq.6.and.k.eq.6) THEN
								aa(i,j,k) = cff(i,6)
							ELSEIF(j.eq.7.and.k.eq.7) THEN
								aa(i,j,k) = cff(i,7)
							ELSEIF(j.eq.8.and.k.eq.8) THEN
							  aa(i,j,k) = 9999.00
							ENDIF			   	
						ELSE
							aa(i,j,k)=0.0				
						ENDIF					  		   	       
				   ENDDO  !k	
			   ENDDO  !j 	   
		   ENDDO !i	
		ENDIF 
!   2	format(8f15.2)
!
!      Initialization
	
		IF(jday.eq.1.and.iclock.eq.0.0) THEN
			DO j = 1, size !mass, making sure all arrays are zero
				DO k = 1, size !length
					DO i = 1, ns !this is the number of layers
						b(i,j,k) = 0
						c(i,j,k) = 0
						d1(i,j,k) = 0
						d2(i,j,k) = 0		
						bb1(i,j,k) = 0
						c1(i,j,k) = 0		   
					ENDDO
				ENDDO
			ENDDO
		ENDIF	 

		DO i = 1, ns !calculations for mass and length conservation
			DO j = 1, size - 1  !mass
				DO k = 1, size - 1	!length
					masssumorg = masssumorg  +aa(i,j,k) * mass(j) + c(i,j,k)*aa(i,j,k)
					fraclensumorg = fraclensumorg +aa(i,j,k) * fraclen(k) + b(i,j,k) *aa(i,j,k)
				ENDDO
			ENDDO
		ENDDO
!	print*,b(100,1,1),c(100,1,1)
!*****************************
!*********TIME STEP BEGINS******
!**********************************
!	DO time = 1, int(nosecs/Delta_T)
		count3 = 0	  
		masssum = 0
		fraclensum = 0
		DO ln = 1,ns
!	    count3=0
!	 	DO kkk = 1, size-1
!	       DO LLL = 1, size-1
!	         IF(aa(ln, kkk, LLL).gt.1) THEN
!	            count3=count3+1
!	         ENDIF
!	 	   ENDDO
!	    ENDDO
	 		count3 = COUNT(aa > 1)/ns
 !        print*, count3
 !        pause
         IF(wqual(ln,1).le.1.0d0) THEN
             TEP(ln)=1.0d0 !*wqual(f,1)
         ELSE
             TEP(ln)=1.0*wqual(ln,1)
         ENDIF 
			DO k = 1, (size - 1) !mass
				DO L = 1, (size - 1)  !length
					DO ii = 1, (size - 1) ! mass equate  (for two similar	k==ii)
						DO j = 1, (size - 1)  ! length equate (for two similar	L==j)
							IF (aa(ln,ii,j) < 1 .OR. aa(ln,k,L) < 1) CYCLE !Dav 7/9/04	
								masssum = 0
								fraclensum = 0		
								betads = 0
								betash = 0
								betabr = 0
								beta = 0
								change = 0
								coag = 0.000
!---------------------------------------------------------------------------------------------
						IF(L==1.or.J==1) THEN       
						   IF(L==J.and.k==ii) THEN
								coag=2.0d-7*TEP(ln)*DEXP(-(2.0d6/((1.0d0/dd(L))+(1.0d0/dd(j)))))*	&	!3.5d-5
          						(DLOG10(aa(ln, k,L))**3.0)*(DLOG10(aa(ln, ii,j))**3.0)          				
							ELSE
								COAG=0.0
							ENDIF		
						ELSEIF(L==2.or.J==2) THEN       !5 and 5
							IF(L==J.and.k==ii) THEN
								coag=7.0d-5*DEXP(-(2.0d6/((1.0d0/dd(L))+(1.0d0/dd(j)))))*	&	!3.5d-5
          						(DLOG10(aa(ln, k,L))**3.0)*(DLOG10(aa(ln, ii,j))**3.0)         				
							ELSE
								COAG=0.0
							ENDIF		
 						ELSEIF(L==3.or.J==3) THEN       !5 and 5	
							IF(L.EQ.J) THEN
								coag=5.0d-3*TEP(ln)*DEXP(-(2.0d6/((1.0d0/dd(L))+(1.0d0/dd(j)))))*	&	!3.5d-5
          						(DLOG10(aa(ln, k,L))**3.0)*(DLOG10(aa(ln, ii,j))**3.0)          				
							ELSE
								COAG=0.0
							ENDIF					           
					  ELSEIF(L==4.or.J==4) THEN       !5 and 5							
							IF(L==J.and.k==ii) THEN
								coag=5.0d-2*TEP(ln)*DEXP(-(2.0d6/((1.0d0/dd(L))+(1.0d0/dd(j)))))*	&	!3.5d-5
          						(DLOG10(aa(ln, k,L))**3.0)*(DLOG10(aa(ln, ii,j))**3.0)        				
							ELSE
								COAG=0.0
							ENDIF	
					ELSEIF(L==5.or.J==5) THEN       !5 and 5					
							IF(L==J.and.k==ii) THEN
								coag=3.0d+1*TEP(ln)*DEXP(-(2.0d6/((1.0d0/dd(L))+(1.0d0/dd(j)))))*	&	!3.5d-5
          						(DLOG10(aa(ln, k,L))**3.0)*(DLOG10(aa(ln, ii,j))**3.0)          				
							ELSE
								COAG=0.0
							ENDIF	
				   ELSEIF(L==6.or.J==6) THEN										          
							IF(L==J.and.k==ii) THEN
								coag=1.0d+5*TEP(ln)*DEXP(-(2.0d6/((1.0d0/dd(L))+(1.0d0/dd(j)))))*	&  !3.5d-5
									(DLOG10(aa(ln, k,L))**3.0)*(DLOG10(aa(ln, ii,j))**3.0)
							ELSE
								COAG=0.0
							ENDIF		
               ELSEIF(L>=7.or.J>=7) THEN
							IF(L==J.and.k==ii) THEN			          
								coag=1.0d-6*TEP(ln)*DEXP(-(2.0d6/((1.0d0/dd(L))+(1.0d0/dd(j)))))*	&  !3.5d-5
									(DLOG10(aa(ln, k,L))**3.0)*(DLOG10(aa(ln, ii,j))**3.0)
								COAG=0.0
							ELSE
								COAG=0.0
							ENDIF			
				  ENDIF
!-------------------------------------------------------------------------------------		 
!						IF(L.le.6.and.J.le.6) THEN
!						   coag=0.010
!						   IF (SD.lt.15.0d0) coag=coag/10.0d0
!						ELSE
!						   coag=0.0
!						ENDIF
						ri1= (dd(L)/2 +  b(ln,k,L))   **(1.0D0/fracdim(L))
						ri2= (dd(L)/2 + (b(ln,k,L)/2) **(1.0D0/fracdim(L)))
						rj1= (dd(j)/2 +  b(ln,ii,j))  **(1.0D0/fracdim(j))
						rj2= (dd(j)/2 + (b(ln,ii,j)/2)**(1.0D0/fracdim(j)))
						mi =  mass(k)  + c(ln,k,L)
						mj =  mass(ii) + c(ln,ii,j)
						wi =  gr*mi/(6.0D0*pi*dyvis(ln)*ri2)
						wj =  gr*mj/(6.0D0*pi*dyvis(ln)*rj2)
						pp =  MIN (ri1,rj1)/MAX(ri1,rj1)
						Di =  kboltz * (temp(ln)+273.15)/(6.0D0 * pi *dyvis(ln) * ri2) 
						Dj =  kboltz * (temp(ln)+273.15)/(6.0D0 * pi *dyvis(ln) * rj2)
						
						Betash = (pp**0.88)*1.3D0*((eps(ln)/kinvis(ln))**0.5D0)*(ri2 +rj2) ** 3
							IF ((wi-wj)*(ri2-rj2).ge. 0.0 ) THEN
								Betads = (pp** 0.984)*pi *((ri2+rj2)**2.0D0)* ABS(wj-wi)
							ELSEIF ((wi-wj)*(ri2-rj2).lt. 0.0 ) THEN 
								Betads = pi *((ri2+rj2)**2.0D0)* ABS(wj-wi)
							END IF
							Betabr = 4.0D0 * pi * ( Di+Dj)* (ri2+rj2)   
							beta = betash + betads + betabr 	
!     check**************************************	 		  
							IF ( k == ii .AND. L == j) THEN
								change=0.5d0*coag*beta*aa(ln,k,L)*aa(ln,ii,j) !*(1/count3)**2							  
							ELSE
								change=1.0d0*coag*beta*aa(ln,k,L)*aa(ln,ii,j)!*(1/count3)**2 ! compared to 1.0  							                     
							ENDIF		
!				  if (ln==1) then
!					  print*,k,ii,L,J, beta*aa(ln,k,L)*aa(ln,ii,j),coag, change,aa(ln,k,L),aa(ln,ii,j)
!					  pause
!					  endif		                  
  
							IF(k == ii .AND. L == j .AND. change >= aa(ln,ii,j)*0.9)THEN !Dav 18/8/04 accounting for too much change
								change = aa(ln,ii,j) *0.9	          
								bbb = bbb + 1
							ELSE IF ( change >= MIN( aa(ln,k,L), aa(ln,ii,j))) THEN
								change = MIN( aa(ln,k,L), aa(ln,ii,j))
								ccc = ccc + 1
							ENDIF									
! remember to add condition about what happens if the 2 particles exceed the matrix	
							IF (((mass(k) + mass(ii) + c(ln,k,L) + c(ln,ii,j))>=					&
     							MAX( newmass(k+1), newmass(ii+1) ) ) .OR. ( ( fraclen(L) +				&
     							fraclen(j) + b(ln,k,L)+b(ln,ii,j))>= MAX(newfraclen(L+1),newfraclen(j+1)))) THEN 		
			 

								DO m1 = 1 , ( size - 1 )
									IF ((mass(k) + mass(ii) + c(ln,k,L) + c(ln,ii,j)) >= newmass(m1)) THEN
										count1 = count1 + 1.
									ENDIF 				
								ENDDO !m1

								DO n1 = 1 , ( size - 1 )
									IF ((fraclen(L) + fraclen(j) + b(ln,k,L) + b(ln,ii,j)) .ge. newfraclen(n1)) THEN
										count2 = count2 + 1
									ENDIF 
								ENDDO !n1
								IF(count1.le.0) THEN
									count1=1
								ENDIF

								IF(count2.le.0) THEN
									count2=1
								ENDIF
								IF ((aa(ln,count1,count2) + change).gt. 0.0 ) THEN
									bb1(ln,count1,count2) =  bb1(ln,count1,count2) +		&
     									(fraclen(L) + fraclen(j) + b(ln,k,L) + b(ln,ii,j)- fraclen(count2) )	* change
									c1(ln,count1,count2) = c1(ln,count1,count2) +			&
     									(mass(k) + mass(ii) + c(ln,k,L) + c(ln,ii,j)- mass(count1)) * change		
								ENDIF
			
								IF ( k==ii .AND. L==j ) THEN
									d2(ln,k,L) = d2(ln,k,L) - change* 2.0
								ELSE
									d2(ln,k,L) = d2(ln,k,L) - change
									d2(ln,ii,j) = d2(ln,ii,j) - change
								ENDIF
								d1(ln,count1,count2) = d1(ln,count1,count2) + change		
							ELSEIF ((j.gt.L).AND.(ii.gt.k).AND.(mass(k)+mass(ii)+				&
									c(ln,k,L) + c(ln,ii,j)) < MAX( mass(k+1), mass(ii+1))		&
									.AND. fraclen(L) + fraclen(j) + b(ln,k,L) + b(ln,ii,j)	&
									< MAX( fraclen(L+1), fraclen(j+1) ) ) THEN 
								IF (aa(ln,ii,j) > 0 ) THEN
								  bb1(ln,ii,j)= bb1(ln,ii,j)+(fraclen(L)+b(ln,k,L))*change
								  c1(ln,ii,j)=c1(ln,ii,j)+(mass(k)+c(ln,k,L))*change					
								END IF				
								d2(ln,k,L) = d2(ln,k,L) - change			
							ELSEIF ((j < L) .AND. (ii < k) .AND. mass(k) + mass(ii) +		&
								c(ln,k,L) + c(ln,ii,j) < MAX( mass(k+1), mass(ii+1) )			&
								.AND. fraclen(L) + fraclen(j) + b(ln,k,L) + b(ln,ii,j)		&
								< MAX( fraclen(L+1), fraclen(j+1) ) ) THEN
								IF (aa(ln,k,L) > 0 ) THEN
									bb1(ln,k,L)=bb1(ln,k,L)+(fraclen(j)+b(ln,ii,j))*change
									c1(ln,k,L)=c1(ln,k,L)+(mass(ii)+c(ln,ii,j))*change					
								END IF				
								d2(ln,ii,j) = d2(ln,ii,j) - change
							ELSE IF ( (j > L) .AND. (ii < k) .AND. mass(k) + mass(ii) +		&
     							c(ln,k,L) + c(ln,ii,j) < MAX(mass(k+1), mass(ii+1))			&
     							.AND. fraclen(L) + fraclen(j) + b(ln,k,L) + b(ln,ii,j)		&
     							< MAX( fraclen(L+1), fraclen(j+1))) THEN
								rand1 =  RAN( jj1 )
								IF  ( rand1 >= 0.5D0 ) THEN !IF rand >= 0.5then we will make larger mass dominate - in this case k,L
									IF (aa(ln,k,L) > 0 ) THEN
										bb1(ln,k,L)=bb1(ln,k,L)+(fraclen(j)+b(ln,ii,j))*change
										c1(ln,k,L)=c1(ln,k,L)+(mass(ii)+c(ln,ii,j))*change						
									ENDIF				 				
									d2(ln,ii,j) = d2(ln,ii,j) - change
									
								ELSE !in this case ii,j dominate
									IF (aa(ln,ii,j) > 0 ) THEN
										bb1(ln,ii,j)=bb1(ln,ii,j)+(fraclen(L)+b(ln,k,L))*change
										c1(ln,ii,j)=c1(ln,ii,j)+(mass(k)+c(ln,k,L))*change						
									ENDIF					
									d2(ln,k,L) = d2(ln,k,L) - change
					
								ENDIF
							ELSEIF ( (j<L) .AND. (ii>k) .AND. mass(k) + mass(ii) +			&
     							c(ln,k,L) + c(ln,ii,j) < MAX( mass(k+1), mass(ii+1))			& 
     			            .AND. fraclen(L)+fraclen(j)+b(ln,k,L)+b(ln,ii,j)				&
     							< MAX( fraclen(L+1), fraclen(j+1))) THEN
								rand2 =  RAN( jj1 )
								IF ( rand2 >= 0.5D0 ) THEN ! IF rand >49 THEN k,L dominate
									IF (aa(ln,k,L) > 0 ) THEN
										bb1(ln,k,L)=bb1(ln,k,L)+(fraclen(j)+b(ln,ii,j))*change
										c1(ln,k,L)=c1(ln,k,L)+(mass(ii)+c(ln,ii,j))*change						
									ENDIF					
										d2(ln,ii,j) = d2(ln,ii,j) - change	
								ELSE ! ii,j dominate
									IF (aa(ln,ii,j) > 0 ) THEN
										bb1(ln,ii,j)=bb1(ln,ii,j)+(fraclen(L)+b(ln,k,L))*change
										c1(ln,ii,j)=c1(ln,ii,j)+(mass(k)+c(ln,k,L))*change										
									ENDIF					
									d2(ln,k,L) = d2(ln,k,L) - change				
								ENDIF
							ELSEIF ( k == ii .AND. L == j .AND. mass(k) + mass(ii) +		&
     							c(ln,k,L) + c(ln,ii,j) < MAX( mass(k+1), mass(ii+1))		&
     							.AND.   fraclen(L) + fraclen(j) + b(ln,k,L) +				&
     							b(ln,ii,j)<MAX(fraclen(L+1), fraclen(j+1))) THEN ! Dav 7/7/04 - what happens when 2 particles of the same bin meet but do not go into next bin
								IF (aa(ln,ii,j) - change /= 0 ) THEN
									bb1(ln,ii,j)=bb1(ln,ii,j)+(fraclen(L)+b(ln,k,L))*change
									c1(ln,ii,j)=c1(ln,ii,j)+(mass(k)+c(ln,k,L))*change
				
								ENDIF
								d2(ln,ii,j) = d2(ln,ii,j) - change 
				
							ELSEIF ( k == ii .AND. L > j .AND. mass(k) + mass(ii) +		&
     							c(ln,k,L) + c(ln,ii,j) < MAX( mass(k+1), mass(ii+1))		&
     							.AND. fraclen(L)+fraclen(j) + b(ln,k,L)+b(ln,ii,j)			&
     							< MAX( fraclen(L+1), fraclen(j+1))) THEN !Dav 8/7/04 - taking into account all possible interactions
								IF (aa(ln,k,L) + change > 0 ) THEN
									bb1(ln,k,L)=bb1(ln,k,L)+(fraclen(j)+b(ln,ii,j))*change
									c1(ln,k,L)=c1(ln,k,L)+(mass(ii)+c(ln,ii,j))*change 									
								ENDIF
								
								d2(ln,ii,j) = d2(ln,ii,j) - change

							ELSEIF ( k == ii .AND. L < j .AND. mass(k) + mass(ii) +		&
     							c(ln,k,L) + c(ln,ii,j) < MAX( mass(k+1), mass(ii+1))		&
     							.AND.fraclen(L) + fraclen(j)+ b(ln,k,L) + b(ln,ii,j)		&
     							< MAX( fraclen(L+1), fraclen(j+1))) THEN
								IF (aa(ln,ii,j) + change > 0 ) THEN
									bb1(ln,ii,j)=bb1(ln,ii,j)+(fraclen(L)+b(ln,k,L))*change
									c1(ln,ii,j)=c1(ln,ii,j)+(mass(k)+c(ln,k,L))*change					
								ENDIF
				
								d2(ln,k,L) = d2(ln,k,L) - change

							ELSEIF ( k > ii .AND. L == j .AND. mass(k) + mass(ii) +		&
     							c(ln,k,L) + c(ln,ii,j) < MAX( mass(k+1), mass(ii+1))		&
     							.AND. fraclen(L) + fraclen(j) + b(ln,k,L)+b(ln,ii,j)		&
     							< MAX( fraclen(L+1), fraclen(j+1) ) ) THEN
								IF (aa(ln,k,L) + change > 0 ) THEN
									bb1(ln,k,L)=bb1(ln,k,L)+(fraclen(j)+b(ln,ii,j))*change
									c1(ln,k,L)=c1(ln,k,L)+(mass(ii)+c(ln,ii,j))*change										
								ENDIF
				
								d2(ln,ii,j) = d2(ln,ii,j) - change

							ELSEIF ( k < ii .AND. L == j .AND. mass(k) + mass(ii) +		&
     							c(ln,k,L) + c(ln,ii,j) < MAX( mass(k+1), mass(ii+1))		&
     							.AND. fraclen(L) + fraclen(j) + b(ln,k,L)+b(ln,ii,j)		&
     							< MAX( fraclen(L+1), fraclen(j+1) ) ) THEN
								IF (aa(ln,ii,j) + change > 0 ) THEN
									bb1(ln,ii,j)=bb1(ln,ii,j)+(fraclen(L)+b(ln,k,L))*change
									c1(ln,ii,j)=c1(ln,ii,j)+(mass(k)+c(ln,k,L))*change					
								ENDIF				
								d2(ln,k,L) = d2(ln,k,L) - change 
							ENDIF	
							count1 = 0
							count2 = 0
							DO i3 = 1, size - 1
								DO j3 = 1, size - 1
									masssum = masssum + aa(ln,i3,j3) * mass(i3) + c(ln,i3,j3) * aa(ln,i3,j3)
									fraclensum = fraclensum + aa(ln,i3,j3) * fraclen(j3) + b(ln,i3,j3) * aa(ln,i3,j3)
	 							ENDDO !j3
							ENDDO	!i3			 
						ENDDO ! j
					ENDDO !ii
				ENDDO !L
			ENDDO !k


			DO L = 1, size-1 !Dav 8/7/04 trying to count for excess mass and length after all interactions
				DO L5 = 1, size-1
					IF ((aa(ln,L,L5) + d1(ln,L,L5) + d2(ln,L,L5)) /= 0 ) THEN	
						b(ln,L,L5) = ((aa(ln,L,L5) + d2(ln,L,L5))*b(ln,L,L5) +				&
     								bb1(ln,L,L5)) / (aa(ln,L,L5) + d2(ln,L,L5) + d1(ln,L,L5))
						c(ln,L,L5) = ((aa(ln,L,L5) + d2(ln,L,L5)) * c(ln,L,L5) +				&
     								c1(ln,L,L5)) / (aa(ln,L,L5) + d2(ln,L,L5) + d1(ln,L,L5))  

						IF ( c(ln,L,L5) >= mass(L).OR. b(ln,L,L5)>= fraclen(L5)) THEN  
							IF(c(ln,L,L5) >= mass(L)) THEN
								aa(ln,L,L5+1)=	aa(ln,L,L5+1)+(c(ln,L,L5)/mass(L))*aa(ln,L,L5+1)
								c(ln,L,L5)=0.0d0
								b(ln,L,L5)=0.0d0
							ELSEIF(b(ln,L,L5) >= fraclen(L5)) THEN	        
								aa(ln,L,L5+1)=	aa(ln,L,L5+1)+(c(ln,L,L5)/mass(L))*aa(ln,L,L5+1)
								c(ln,L,L5)=0.0d0
								b(ln,L,L5)=0.0d0  
							ENDIF 
						ENDIF 
						IF(b(ln,L,L5).le.0.0d0) b(ln,L,L5)=0.0d0
						IF(c(ln,L,L5).le.0.0d0) c(ln,L,L5)=0.0d0
					ELSE 
						b(ln,L,L5) = 0.0d0
						c(ln,L,L5) = 0.0d0
					ENDIF
      
					aa(ln,L,L5) = aa(ln,L,L5)+ d1(ln,L,L5) + d2(ln,L,L5)
!               aa(ln,L,L5) = aa(ln,L,L5)+ d2(ln,L,L5)
!               aa(ln,L+1,L5+1) = aa(ln,L+1,L5+1)+ d1(ln,L,L5) 
					IF (aa(ln,L,L5) <= 0 ) THEN
						b(ln,L,L5) = 0.0d0
						c(ln,L,L5) = 0.0d0
					ENDIF				   
					d1(ln,L,L5) = 0.0d0
					d2(ln,L,L5) = 0.0d0
					bb1(ln,L,L5) = 0.0d0
					c1(ln,L,L5) = 0.0d0								
!***************Sinking depth****************************************************	
					IF (L.eq.1 .AND.L5.eq.1) THEN
						setldep(ln,L,L5) = 2.0d0*gr *(3.0D0/fracdim(L5))*			      &
     					(ddn(L5)/2.0D0 + (( b(ln,L,L5)) **										&
     					(1.0D0 / fracdim(L5))) / 2.0D0) ** 2.0D0 * 						   &   
     					(( 1.00D0 * inorgden(L)+ 0.00D0 * orgden) - denh2o(ln))/(9.0D0 * dyvis(ln))
     !				if(ln.eq.ns) then
     !					print*,'1',setldep(ln,L,L5),inorgden(L),denh2o(ln),dyvis(ln)/3600.0
     !					pause
     !					endif     				
					ELSEIF ( L.eq.2.and.L5.eq.2) THEN
						setldep(ln,L,L5) = 2.0d0*gr *(3.0D0/fracdim(L5))*			      &
     						(ddn(L5)/2.0D0 + ((b(ln,L,L5)) **										&
     						(1.0D0/fracdim(L5))) / 2.0D0) ** 2.0D0 *							&   
     						((1.00D0 * inorgden(L)+ 0.00D0 * orgden) - denh2o(ln))/(9.0D0 * dyvis(ln))
					ELSEIF ( L.eq.3.and.L5.eq.3) THEN
						setldep(ln,L,L5) = 2.0d0*gr *(3.0D0/fracdim(L5))*			      &
     						(ddn(L5)/2.0D0 + ((b(ln,L,L5)) **										&
     						(1.0D0 / fracdim(L5))) / 2.0D0) ** 2.0D0 *						&    
     						((1.0D0 * inorgden(L)+ 0.0D0 * orgden) - denh2o(ln))/(9.0D0 * dyvis(ln))      			
					ELSEIF ( L.eq.4.and.L5.eq.4) THEN
						setldep(ln,L,L5) = 2.0d0*gr *(3.0D0/fracdim(L5))*			      &
     						(ddn(L5)/2.0D0 + ((b(ln,L,L5)) **										&
     						(1.0D0 / fracdim(L5))) / 2.0D0) ** 2.0D0 *						&		
     						((0.60D0 * inorgden(L)+ 0.40D0 * orgden) - denh2o(ln))/(9.0D0 * dyvis(ln)) 
					ELSEIF ( L.eq.5.and.L5.eq.5) THEN
						setldep(ln,L,L5) = 2.0d0*gr *(3.0D0/fracdim(L5))*			      &
     						(ddn(L5)/2.0D0 + ((b(ln,L,L5)) **										&
     						(1.0D0 / fracdim(L5))) / 2.0D0) ** 2.0D0 *						&		
     						((0.30D0 * inorgden(L)+ 0.70D0 * orgden) - denh2o(ln))/(9.0D0 * dyvis(ln)) 

					ELSEIF ( L.eq.6.and.L5.eq.6) THEN
						setldep(ln,L,L5) = 2.0d0*gr*(3.0D0/fracdim(L5))*			      &
     						(ddn(L5)/2.0D0 + (( b(ln,L,L5)) **									&
     						(1.0D0 / fracdim(L5))) / 2.0D0) ** 2.00D0 *						&		
     						((0.30D0 * inorgden(L)+ 0.70D0 * orgden) - denh2o(ln))/(9.0D0 * dyvis(ln))
					ELSE
						setldep(ln,L,L5) = 2.0d0*gr*(3.0D0/fracdim(L5))*			      &
     						(ddn(L5)/2.0D0+((b(ln,L,L5))**(1.0D0/								&
     						fracdim(L5)))/2.0D0)**2.00D0* 										&		
     						((0.30D0*inorgden(L)+0.70D0*orgden) - denh2o(ln))/(9.0D0 * dyvis(ln))		
					ENDIF
 		
					IF (aa(ln,L,L5) < 1.0D0 ) THEN 
						aa(ln,L,L5) = 1.0D0
					ENDIF
				ENDDO ! L
			ENDDO !L5
			

!	 WRITE (*,*)  time
			count5 = COUNT( aa < 0 )
			IF ( count5 >= 1 ) THEN
				EXIT 
			ENDIF
			count4 = 0
		ENDDO ! ln, this is DO loop for the ln, the one that goes through the different layers
!	END DO !time, this DO loop is for the number of the iterations 
!************  PARTICLE SETTLING *********************
      IF(partsw) THEN
	      index1=1
		ELSE
	      index1=8
		ENDIF
!
      DO k = 1,sizall
			Before(k)   = 0.0d0
			After(k)    = 0.0d0
			Lost(k)     = 0.0d0
      ENDDO    ! k

 		DO k = 1,7
			DO i = 1,ns
				Before(k)  = Before(k) + cff(i,k)*vol(i)   
			ENDDO   ! i
		ENDDO   ! k 

      DO i = 1,ns
			Before(8)  = Before(8)   + wqual(i,9) *vol(i)    
			Before(9)  = Before(9)   + wqual(i,13)*vol(i) 
			Before(10) = Before(10)  + wqual(i,19)*vol(i) 
      ENDDO  ! i
      
		DO k = index1,size-1
			DO m = index1, size-1
				DO i = ns,1,-1
				if(setldep(i,k,m).ge.depth(i)) setldep(i,k,m)=depth(i) 
					IF(i.eq.1) THEN
						aadv(i,k,m)=0.0
						depleft(i,k,m)=depth(i)-setldep(i,k,m)
						IF(depleft(i,k,m).le.0.0) depleft(i,k,m)=0.0 
						deplefb(i,k,m)=0.0
					ELSE
						depleft(i,k,m)=depth(i)-setldep(i,k,m)
						deplefb(i,k,m)=depth(i-1)-setldep(i,k,m)		
						IF(setldep(i,k,m).lt.dep(i))THEN
							aadv(i,k,m) = 1.0*aa(i,k,m)*setldep(i,k,m)*arfac*0.5*(area(i)+area(i-1)) !# particles					
						ELSEIF(setldep(i,k,m).ge.dep(i))THEN
							aadv(i,k,m) = 1.0*aa(i,k,m)*(dep(i)*arfac*0.5*(area(i)+area(i-1)))	!# particlesc			  
						ENDIF
					ENDIF
!	     print*,i, aadv(i,k,m)
		  			fluxin(i,k,m) = 0.0d0
					b_prime(i,k,m)=0.0d0
					c_prime(i,k,m)=0.0d0
!					aadv(i,k,m)=0.0d0
				ENDDO ! i
			ENDDO ! m
		ENDDO ! k
! convert cff back to a particle concentration
		DO k = index1,size-1
			DO m = index1, size-1
				DO i = ns,1,-1
					DO j = i,1,-1
! case where all settling particles fall within top and bottom of layer 1:
						IF(j.eq.1)THEN
! case where settling particles fall within layer 1 and extend furtherthan the bottom:
							IF(depleft(i,k,m).le.depth(j).and.depleft(i,k,m).gt.0.0.and.deplefb(i,k,m).le.0.0)	&
     							fluxin(j,k,m) = fluxin(j,k,m)+ aadv(i,k,m)
!								b_prime(j,k,m)=b_prime(j,k,m)+aadv(i,k,m)*b(i,k,m)
!								c_prime(j,k,m)=c_prime(j,k,m)+aadv(i,k,m)*c(i,k,m)
! case where settling particles extend above the top of layer 1 and extend further than the bottom:
							IF(depleft(i,k,m).gt.depth(j).and.deplefb(i,k,m).le.0.0)			&
     	    					fluxin(j,k,m) = fluxin(j,k,m)+ aadv(i,k,m)
!								b_prime(j,k,m)=b_prime(j,k,m)+aadv(i,k,m)*b(i,k,m)
!								c_prime(j,k,m)=c_prime(j,k,m)+aadv(i,k,m)*c(i,k,m)

! case where all settling particles fall within layer 1:
							IF(depleft(i,k,m).le.depth(j).and.deplefb(i,k,m).gt.0.0)			&
     							fluxin(j,k,m) = fluxin(j,k,m)+ aadv(i,k,m)
!								b_prime(j,k,m)=b_prime(j,k,m)+aadv(i,k,m)*b(i,k,m)
!								c_prime(j,k,m)=c_prime(j,k,m)+aadv(i,k,m)*c(i,k,m)

! case where settling particles extend above layer 1 and fall within layer 1:
							IF(depleft(i,k,m).gt.depth(j).and.deplefb(i,k,m).lt.depth(j).and.deplefb(i,k,m).gt.0.0)	&
     	   					fluxin(j,k,m) = fluxin(j,k,m)+ aadv(i,k,m)
!								b_prime(j,k,m)=b_prime(j,k,m)+aadv(i,k,m)*b(i,k,m)
!								c_prime(j,k,m)=c_prime(j,k,m)+aadv(i,k,m)*c(i,k,m)
						ELSE

! case where all settling particles fall within the top and bottom of a layer
							IF(depleft(i,k,m).le.depth(j).and.deplefb(i,k,m).ge.depth(j-1)) THEN
      						fluxin(j,k,m) = fluxin(j,k,m)+ aadv(i,k,m)
!								b_prime(j,k,m)=b_prime(j,k,m)+aadv(i,k,m)*b(i,k,m)
!								c_prime(j,k,m)=c_prime(j,k,m)+aadv(i,k,m)*c(i,k,m)

! case where settling particles extend above the top of a layer and below the bottom of the layer
							ELSEIF(depleft(i,k,m).ge.depth(j).and.deplefb(i,k,m).le. depth(j).and.	&
     							deplefb(i,k,m).ge. depth(j-1)) THEN
      						fluxin(j,k,m) = fluxin(j,k,m)+ aadv(i,k,m)
!								b_prime(j,k,m)=b_prime(j,k,m)+aadv(i,k,m)*b(i,k,m)
!								c_prime(j,k,m)=c_prime(j,k,m)+aadv(i,k,m)*c(i,k,m)
							ENDIF
						ENDIF
					ENDDO ! j
				ENDDO	! i
			ENDDO !m
		ENDDO ! k
!*** EQUAL MASS AND VOLUME DISTRIBUTION IN A LAYER************
      DO k = index1,size
			DO m = index1, size 
				DO i = ns,1,-1  	    
					IF(i.eq.1) THEN
	       			aa(i,k,m) = aa(i,k,m)+(fluxin(i,k,m)-aadv(i,k,m))/(dep(i)*area(i)*arfac)  !#particles/m^3 
						IF(aa(i,k,m).le.0)THEN
	       				aa(i,k,m)=0.0
							b(i,k,m)	=0.0
							c(i,k,m)= 0.0
						ELSE
							b(i,k,m)=((aa(i,k,m)*dep(i)*area(i)*arfac-aadv(i,k,m))*b(i,k,m)+			&
     							b_prime(i,k,m))/(aa(i,k,m)*dep(i)*area(i)*arfac+fluxin(i,k,m)-aadv(i,k,m))
			
							c(i,k,m)=((aa(i,k,m)*dep(i)*area(i)*arfac-aadv(i,k,m))*c(i,k,m)+			&
     								c_prime(i,k,m))/(aa(i,k,m)*dep(i)*area(i)*arfac+fluxin(i,k,m)-aadv(i,k,m))
						ENDIF
	       
					ELSE
						aa(i,k,m) = aa(i,k,m)+(fluxin(i,k,m)-aadv(i,k,m))/(dep(i)*arfac*0.5*(area(i)+area(i-1)))      !#particles/m^3

						IF(aa(i,k,m).le.0)THEN
							aa(i,k,m)=0.0
							b(i,k,m)	=0.0
							c(i,k,m)= 0.0
						ELSE

							b(i,k,m)=((aa(i,k,m)*dep(i)*0.5*(area(i)+area(i-1))*arfac-						&
     								aadv(i,k,m))*b(i,k,m)+b_prime(i,k,m))/(aa(i,k,m)*0.5*(area(i)+			&
     								area(i-1))*arfac+fluxin(i,k,m)-aadv(i,k,m))
			
							c(i,k,m)=((aa(i,k,m)*dep(i)*0.5*(area(i)+area(i-1))*arfac-						&			
     								aadv(i,k,m))*c(i,k,m)+c_prime(i,k,m))/(aa(i,k,m)*0.5*(area(i)+			&
									area(i-1))*arfac+fluxin(i,k,m)-aadv(i,k,m))
						ENDIF
					ENDIF	  
					IF(aa(i,k,m).le.0.0) aa(i,k,m)=0.0	
				ENDDO  ! i
			ENDDO !m
		ENDDO ! k
    
		masssumorg1 = 0
		fraclensumorg1 = 0
      DO ln=1,ns  ! Layer
			DO L = 1, size-1 !mass
				DO L5 = 1, size-1 !length
					IF (aa(ln,L,L5) == 0 ) CYCLE
					IF (L == size-1 .OR. L5 == size-1 ) THEN
						IF (b(ln,L,L5) >= fraclen(L5).AND.c(ln,L,L5)>= mass(L))THEN
							b(ln,L,L5) =  ( b(ln,L,L5) - fraclen(L5) ) / 2 !Dav 7/13/04 - accounting for doubling of particles
							c(ln,L,L5) =  ( c(ln,L,L5) - mass(L) ) / 2
							aa(ln,L,L5) = aa(ln,L,L5) * 2
							bbb = bbb + 1
						END IF
					ELSE IF (b(ln,L,L5) >= fraclen(L5 + 1) .AND. c(ln,L,L5) >= mass(L + 1) ) THEN
						b(ln,L,L5) = b(ln,L,L5) - fraclen(L5 + 1)	
						c(ln,L,L5) = c(ln,L,L5) - mass(L + 1)
						aa(ln,L + 1,L5 + 1) = aa(ln,L + 1,L5 + 1) + aa(ln,L,L5)
						b(ln,L + 1,L5 + 1) = (b(ln,L+1,L5+1)/aa(ln,L+1,L5+1))
						c(ln,L + 1,L5 + 1) = (c(ln,L+1,L5+1)/aa(ln,L+1,L5+1)) 
						bbb = bbb + 1
					ELSE IF (b(ln,L,L5) >= fraclen(L5 + 1) .AND. c(ln,L,L5)>= mass(L).AND.c(ln,L,L5) < mass(L + 1) ) THEN
						b(ln,L,L5) = b(ln,L,L5) - fraclen(L5 + 1)
						c(ln,L,L5) = c(ln,L,L5) - mass(L)
						aa(ln,L,L5+1) = aa(ln,L,L5+1) + aa(ln,L,L5)
						b(ln,L,L5 + 1) = (b(ln,L,L5 + 1) / aa(ln,L,L5 + 1))
						c(ln,L,L5 + 1) = (c(ln,L,L5 + 1) / aa(ln,L,L5 +1))
						bbb = bbb + 1
					ELSE IF ( c(ln,L,L5) >= mass(L + 1) .AND. b(ln,L,L5)>= fraclen(L5) .AND. b(ln,L,L5)<fraclen(L5 + 1)) THEN
						b(ln,L,L5) = b(ln,L,L5) - fraclen(L5)
						c(ln,L,L5) = c(ln,L,L5) - mass(L + 1)
						aa(ln,L + 1,L5) = aa(ln,L + 1,L5) + aa(ln,L,L5)
						b(ln,L + 1,L5) = ( b(ln,L + 1,L5) / aa(ln,L + 1,L5))
						c(ln,L + 1,L5) = ( c(ln,L + 1,L5) / aa(ln,L + 1,L5))
						bbb = bbb + 1
					ENDIF
					IF ( c(ln,L,L5) >= mass(L).OR. b(ln,L,L5)>= fraclen(L5)) THEN
						IF(c(ln,L,L5) >= mass(L)) THEN
							aa(ln,L,L5+1)=	aa(ln,L,L5+1)+aa(ln,L,L5+1)* (c(ln,L,L5)/mass(L))
							c(ln,L,L5)=0.0
							b(ln,L,L5)=0.0
						ELSEIF(b(ln,L,L5) >= fraclen(L5)) THEN	        
							aa(ln,L,L5+1)=	aa(ln,L,L5+1)+aa(ln,L,L5+1)*(c(ln,L,L5)/mass(L))

							c(ln,L,L5)=0.0
							b(ln,L,L5)=0.0		
						ENDIF 
					END IF 
					IF(b(ln,L,L5).le.0.0d0) b(ln,L,L5)=0.0d0
					IF(c(ln,L,L5).le.0.0d0) c(ln,L,L5)=0.0d0
				ENDDO !L  mass
			ENDDO !L5 length
		ENDDO !ln layer
!	 print*,b(60,1,1),fraclen(1),aa(60,1,1)
!	pause
!	ENDDO !time loop
!  
!*******Initialization***************
		DO i = 1,ns
			DO j=1,size
		      cff(i,j) = 0.0
			ENDDO		
		ENDDO
!**************************************
		DO i = 1,ns
			DO j=1,size
			  DO k=1,size
				 IF(j.eq.1) THEN
				   cff(i,j) = cff(i,j)+ aa(i,j,k)+b(i,j,k)/fraclen(j)+c(i,j,k)/mass(k)
				 ELSEIF (j.eq.2) THEN		    	            
				   cff(i,j) = cff(i,j)+ aa(i,j,k)+b(i,j,k)/fraclen(j)+c(i,j,k)/mass(k)
				 ELSEIF (j.eq.3) THEN		    	            
				   cff(i,j) = cff(i,j)+ aa(i,j,k)+b(i,j,k)/fraclen(j)+c(i,j,k)/mass(k)
				 ELSEIF (j.eq.4) THEN		    	            
				   cff(i,j) = cff(i,j)+ aa(i,j,k)+b(i,j,k)/fraclen(j)+c(i,j,k)/mass(k)
				 ELSEIF (j.eq.5) THEN		    	            
				   cff(i,j) = cff(i,j)+ aa(i,j,k)+b(i,j,k)/fraclen(j)+c(i,j,k)/mass(k)	
				 ELSEIF (j.eq.6) THEN		    	            
				   cff(i,j) = cff(i,j)+ aa(i,j,k)+b(i,j,k)/fraclen(j)+c(i,j,k)/mass(k)
				 ELSEIF (j.eq.7) THEN		    	            
				   cff(i,j) = cff(i,j)+ aa(i,j,k)+b(i,j,k)/fraclen(j)+c(i,j,k)/mass(k)
				ENDIF	    
				ENDDO
			ENDDO		
		ENDDO
!      print*, cff(200,1)*vol(200), fluxin(200,1,1)
 		DO k = 1,7
			DO i = 1,ns
				After(k)  = After(k) + cff(i,k)*vol(i)    
			ENDDO   ! i
		ENDDO   ! k 

       
		DO k = 1,sizall-1
			Lost(k) = Before(k) - After(k)
		ENDDO

			DO i = 1,ns
!gbs boundary condition: when the particles hits bottoms layers,
!gbs they settled down permanently. ELSE there is upwelling of the particles during winter mixing.
            DO k=1,7
	 			   IF(i.LE.1) THEN						     
					   cff(i,k) = 1.0
					   IF(ISNAN(cff(i,k)).or.cff(i,k).le.0.0d0)cff(i,k)=10.0d0
					   cf(i,k) = cff(i,k)    
				   ELSEIF (i.gt.1.and.depth(i).lt.(depth(ns)-20.0d0))THEN
				      cff(i,k) = cff(i,k)-SOD_area(i)*0.005d0*cff(i,k)/(0.5*dep(i)*(area(i)+area(i-1))*arfac) 				      
		            IF(ISNAN(cff(i,k)).or.cff(i,k).le.0.0d0)cff(i,k)=10.0d0
		            cf(i,k) = cff(i,k)
				   ELSE 					   
 					   IF(ISNAN(cff(i,k)).or.cff(i,k).le.0.0d0)cff(i,k)=10.0d0
 					   cf(i,k) = cff(i,k)				
	   		   ENDIF 
	   		ENDDO	 !k
			ENDDO	! i
      RETURN
		END SUBROUTINE COAGULATION


 


 
 
