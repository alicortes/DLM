!*******************************************************************
		SUBROUTINE SCHMIDT (uwor,uwor2,lk,Z_min,yeari)
!*******************************************************************
!
!        Schmidt Work
!       ---------------
!
!		This routine calculates the Schmidt Work for the current day
!		and sends the output to a file:'Ldddhhmm.NUM'
!		i.e. has the same name as the simulation file
!
!     NOTES:  nlayers = number of layers of INTERP DATA(STEP = 0.1)
!             all arrays begin with subscrip zero
!             all arrays of data interp to 0.1m have subscript of zero
!             referring to depth of 0.1m (Top of first layer)
!
! Modified for Lake Tahoe
! Written by Joaquim P. Losada and Geoff Schladow
! UC-Davis, January 2001
!*******************************************************************
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE 
		INTEGER(4)     :: uwor,uwor2,nlayers
		INTEGER(4)     :: i,j,ii,ik,jj,ijk,lk,njday
		LOGICAL*4 NOPROBLEMS

		REAL(8)	:: thermdepth
		REAL(8)	:: xp,zc,xmass,VolSum,Mass,den4
		REAL(8)  :: IsoDen,DepthSum 
		REAL(8)	:: Heat, Contingut
		REAL(8)	:: WS,WS_min,WT,WT_min
!IF		REAL(8)  :: densty
		REAL(8)	:: Mix_depth_F,Mix_depth_C
		REAL(8)  :: AreaT, z_star, T_average,T_ave_vol,sum_vol
		REAL(8)  :: T_ave_vol_LT100, sum_vol_LT100, T_ave_vol_GT100,sum_vol_GT100
		REAL(8)  :: heat_con(maxns),depth_con(maxns),depthcon
		INTEGER(4) :: year, yeari 
		REAL(8) :: Z_min
		REAL(8), PARAMETER	:: Amax		= 465969.4*1000.0d0
		REAL(8), PARAMETER	:: Cp		= 4186.0d0
		REAL(8), PARAMETER	:: grav		= 9.81d0
!
!quim Interpolate density and temperature each 0.1m
!
			
		CALL INTERPOLATE_LAYER_DATA(nlayers,idepth,density,tempy,salty)

		IF (nlayers.lt.1) THEN
			noproblems = .false.
			GOTO 4406
		ENDIF
!
!quim  Mixing Depth detection
!
		xmom = CALC_XMOM(nlayers)
		
		IF (xmom .lt. 0.0d0) THEN
			xmom       = 0.10
			thermdepth = 0.10
			GOTO 8182
		ELSE
			thermdepth = xmom 
		ENDIF

8182	CONTINUE

		ii = nlayers-1
		DO WHILE((tempy(nlayers-1)-tempy(ii-1)).le. 0.0000001)
			ii = ii - 1
			IF(ii.eq.0) GOTO 4421
		ENDDO

4421	CONTINUE
		Mix_depth_F = idepth(ii)
	
		ii = nlayers-1
		DO WHILE((tempy(nlayers-1)-tempy(ii-1)).le. 0.01)
			ii = ii - 1
			IF(ii.eq.0) GOTO 4422
		ENDDO

4422	CONTINUE	 
		Mix_depth_C = idepth(ii)
!
!quim Schmidt Work
!
		VolSum  = 0.0d0
		Mass    = 0.0d0
		WS		= 0.0d0	
		
		AreaT   = area(1)*1D6
		VolSum  = VolSum + vol(1)*1D3
		DepthSum	= DepthSum + vol(1)*(depth(ns) - depth(1))*1D3
		Mass    = Mass   + vol(1)*(den(1)+1000.0d0)*1D3

		DO ii =2,ns	
			AreaT   = area(ii)*1D6
			VolSum  = VolSum   + vol(ii)*1D3
			DepthSum= DepthSum + vol(ii)*(depth(ns) - depth(ii))*1D3
			Mass    = Mass     + vol(ii)*(den(ii)+1000.0d0)*1D3
		ENDDO

		DepthSum = dble(DepthSum/VolSum)
		IsoDen   = dble(Mass/VolSum)
	
		ii = 0
		DO WHILE (density(ii)-(IsoDen-1000.0d0).gt.0.0d0)
			IF(ii.eq.nlayers-1) THEN
				WS = 0.0
				GOTO 4406
			ENDIF  
			ii = ii + 1
		ENDDO
		z_star = idepth(nlayers-1) - idepth(ii)
	
		ii = nlayers-1
		DO WHILE (density(ii)-(IsoDen-1000.0d0).lt.0.0d0)
			IF(ii.eq.0) THEN
				WS = 0.0
				GOTO 4406
			ENDIF  
			ii = ii - 1
		ENDDO

		z_star = (z_star + (idepth(nlayers-1) - idepth(ii)))/2.0

		WS  =  vol(1)*1d3*(den(1)+1000.0d0 - IsoDen)*(depth(ns)-depth(1)-z_star) 

		DO ii = 2,ns
			AreaT = area(ii)*1D6
     		WS   = WS + vol(ii)*1D3*(den(ii)+1000.0d0 - IsoDen)*(depth(ns)-depth(ii)-z_star) 
			IF (WS.lt.0.0d0) WS = 0.0d0
		ENDDO

! 04 June 2004 Error area(ii)!!
		WS		= grav*WS/(area(ns)*1D6)

4406  CONTINUE
!
!	Total Work
!
		den4	= 1000.0d0 + densty(4.0d0,100.0d0)
		WT		= (grav*DepthSum*(den4-IsoDen)*(VolSum/(area(ns)*1d6)))

!	WT		= (grav*DepthSum*(den4-IsoDen)*VolSum*
!     &		  (Z_red/DepthMax)**2.0d0)

4207  CONTINUE
!
!	Heat Contens
!
		heat		= 0.0d0
		contingut	= 0.0d0	

		DO ii = 1,ns
			AreaT     = area(ii)*1D3
			contingut = contingut + Cp*(den(ii)+1000.0d0)*vol(ii)*1d3*temp(ii)	!(Joules)
		ENDDO

		heat		= contingut
																	  !(Joules)
		contingut   = Cp*(den(1)+1000.0d0)*vol(1)*1d3*temp(1)	 !(Joules)
		depthcon    = depth(1)
		depth_con   = 0.0d0

		DO ii = 2,ns
			AreaT         = area(ii)*1D3
			contingut     = contingut + Cp*(den(ii)+1000.0d0)*vol(ii)*1d3*temp(ii)	!(Joules)
			heat_con(ii)  = contingut/heat
			depthcon      = depthcon  + depth(ii) - depth(ii-1)
			depth_con(ii) = depthcon /depth(ns)
		ENDDO
!
!gbs  Temperautre volume averaged
!d   
		T_ave_vol = 0.0d0
		sum_vol   = 0.0d0
		T_ave_vol_LT100 = 0.0d0
		sum_vol_LT100   = 0.0d0
		T_ave_vol_GT100 = 0.0d0
		sum_vol_GT100   = 0.0d0

		DO ii = 2,ns
			T_ave_vol = T_ave_vol+temp(ii)*vol(ii)
			sum_vol = sum_vol+vol(ii)
!	PRINT*,'GOLOKA WORK', DEPTH(NS)
!	PAUSE

			IF(DEPTH(II).GE.(DEPTH(NS)-100.0D0)) THEN	!EPILIMNION
				T_ave_vol_LT100 = T_ave_vol_LT100+temp(ii)*vol(ii)
				sum_vol_LT100 = sum_vol_LT100+vol(ii)
			ELSE
				T_ave_vol_GT100 = T_ave_vol_GT100+temp(ii)*vol(ii)	!HYPOLIMNION
				sum_vol_GT100 = sum_vol_GT100+vol(ii)
			ENDIF
		ENDDO
		T_ave_vol=T_ave_vol/sum_vol
		T_ave_vol_LT100=T_ave_vol_LT100/sum_vol_LT100
		T_ave_vol_GT100=T_ave_vol_GT100/sum_vol_GT100

!
!    Temperature Depth Averaged
!
		T_average	= 0.0d0					! Temp (oC)

		DO ii = 0,nlayers-2	   !Gbs changed "nlayers-1" to "nlayers-2"
			T_average	= T_average	+ tempy(ii)*(idepth(ii+1)-idepth(ii))
	!		print*,(idepth(ii+1)-idepth(ii))	
		ENDDO
		T_average	= T_average/idepth(nlayers-1)	!Gbs changed "ns" to "nlayers-1"
      
		WRITE(uwor,4407) jday,depth(ns)-depmx,T_average,T_ave_vol,WS,WT,	&   
								T_ave_vol_LT100,	T_ave_vol_GT100 	

!
!		Maximum Mixing Depth
!	
		year = int(jday/1000.0)
		IF(year.eq.yeari) THEN
			Z_min = min(Z_min,depmx)
			IF(Z_min.eq.depmx) THEN
				WS_min = WS
				WT_min = WT
			ENDIF
		ELSE
			WRITE(uwor2,4409) yeari,DepthMax-Z_min,WS_min,WT_min
			Z_min = 10000.0
			yeari = year
		ENDIF

4407	FORMAT(I8,f10.2,2(3x,f7.2),2(x,f12.1,x),2(8x,f7.2))
4408	FORMAT(I10,f8.1,f8.1)
4409	FORMAT(I10,f10.2,2f15.1)

		RETURN
		END SUBROUTINE SCHMIDT

