!*********************************************************************	
      SUBROUTINE plng_mxng(iriver,BSL)
!*********************************************************************
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE		
	   REAL*8::  Slope				!declare functions
	   REAL*8:: mix_factor	! used to scale gamma returned from Safaie and HandD, set to 1 for AkSt
	   REAL*8:: gammaf		! entrainment (Qtot/Qinit)
	   LOGICAL:: check=.FALSE.
	   LOGICAL:: firstcheck=.FALSE.
!	   LOGICAL:: gpnegcheck=.FALSE.
	   LOGICAL:: lakecheck=.FALSE.
	   LOGICAL:: bypass=.FALSE.
	   INTEGER*4:: max_seg_count,iriver,laycnt,i
	   REAL*8::Lngth,upDpth,dwnDpth,Q,ratio
	   REAL*8:: rivarea(MAXNS),delQ,layvol
      REAL*8:: WQ1,WQ2,WQ3,WQ4,WQ5,WQ6,WQ7,WQ8,WQ9,WQ10,WQ11
      REAL*8:: WQ12,WQ13,WQ14,WQ15,WQ16,WQ17,WQ18,WQ19,WQ20
      REAL*8:: WQ21,WQ22,WQ23,WQ24,WQ25,WQ26,WQ27,WQ28
      REAL*8:: CCFF1,CCFF2,CCFF3,CCFF4,CCFF5,CCFF6,CCFF7
	   REAL*8:: S,T,DI,BSL
!IF   REAL*8::	DENSTY,COMBIN,COMBINV,
	   REAL*8:: Q_initial,T_initial
	   REAL*8:: firstDpth, firstLngth
	   INTEGER*4:: normdepth_segment	
	
     	ratio=0.75
	
	   plngdpth(iriver) = DEPTH(NS)		! initialise to upper bound
	   max_seg_count=numseg(iriver)		! find length of valid # of segments

!     delete this section when testing complete this bypasses the plunge calculations
!	   and returns bsl and plngdpth(iriver) as DEPTH(NS) to force the execution of the original
!		inflow algorithm
	   bypass = .FALSE.
	   IF (bypass) THEN
		   upDpth = (depth(ns)+base)-BgnEle(iriver,max_seg_count-1)
		   dwnDpth = (depth(ns)+base)-BgnEle(iriver,max_seg_count)
		   Lngth = SegLngth(iriver,max_seg_count)
		   BSL = Slope(upDpth,dwnDpth,Lngth)
			Bedslope(iriver)=Slope(upDpth,dwnDpth,Lngth)
		   RETURN
	   ENDIF
!!***** delete this section when testing complete this bypasses the plunge calculations
	   Q = Qdown(iriver,icnt(iriver))/(REAL(nosecs)/1000.0d0)	!change into flow rate
	   Q_initial = Qdown(iriver,icnt(iriver))
	   T_initial = Tdown(iriver,icnt(iriver))
!     for the entrainment calcs	
!	   write(88,*) 'plngmxng iriver,lakecheck',iriver,numseg(iriver),
!        &	lakecheck,nosecs,Q,Tdown(iriver,icnt(iriver)),temp(ns)
!	   try a generic approach for allow river segmentation schemes


	   normdepth_segment = 1	
	   CALL find_normal_segment(iriver,Q,normdepth_segment,max_seg_count,firstDpth,firstLngth)	
	   CALL find_plunge_segment(iriver,Q,normdepth_segment,max_seg_count,firstDpth,firstLngth, &
	                             gammaf,mix_factor)

!     calculate entrainment into current inflow parcel
!	   Note that riverplunge and lakeplunge both cause the variable plngdpth to be defined

	   delQ=Qdown(iriver,icnt(iriver))*(gammaf-1)*mix_factor
	   CALL findplungelayer(iriver,laycnt)			!find the layer to which the plunge occurs


!find the incremental volume from each layer
!wef	t9=0.0
	   DO i=1,laycnt
		   rivarea(i)=(depth(ns+1-i)-depth(ns-i))/plngdpth(iriver)
	   ENDDO
	   rivarea(laycnt+1)=(plngdpth(iriver)-depth(ns)+depth(ns-laycnt))/plngdpth(iriver)

! establish temporary values for the mixing inflow

      Q = QDOWN(IRIVER,icnt(iriver))
      T = TDOWN(IRIVER,icnt(iriver))
      S = SDOWN(IRIVER,icnt(iriver))
      
!	write(88,*)'begin mixing q,t',iriver,icnt(iriver),q,t,delq,gammaf,
!     &			mix_factor,plngdpth(iriver)

!	WRITE(63,974)SIMDAY, iclock, iriver, icnt(iriver), Q,
!     &			 gammaf, plngdpth(iriver), T, temp(ns)
!974	FORMAT(4I10,5F16.4)
!      IF(iriver.eq.36) THEN
!	WRITE(63,974)SIMDAY, iriver, Q,plngdpth(iriver),temp(ns)
!	ENDIF
      974 format(2I10,3F16.4)   
	   DI = DENSTY(T,S)
      CCFF1 = CFDOWN(IRIVER,1,icnt(iriver))
      CCFF2 = CFDOWN(IRIVER,2,icnt(iriver))
      CCFF3 = CFDOWN(IRIVER,3,icnt(iriver))
      CCFF4 = CFDOWN(IRIVER,4,icnt(iriver))
      CCFF5 = CFDOWN(IRIVER,5,icnt(iriver))
      CCFF6 = CFDOWN(IRIVER,6,icnt(iriver))
      CCFF7 = CFDOWN(IRIVER,7,icnt(iriver))
      WQ1 = WQDOWN(IRIVER,1,icnt(iriver))
      WQ2 = WQDOWN(IRIVER,2,icnt(iriver))
      WQ3 = WQDOWN(IRIVER,3,icnt(iriver))
      WQ4 = WQDOWN(IRIVER,4,icnt(iriver))
      WQ5 = WQDOWN(IRIVER,5,icnt(iriver))
      WQ6 = WQDOWN(IRIVER,6,icnt(iriver))
      WQ7 = WQDOWN(IRIVER,7,icnt(iriver))
      WQ8 = WQDOWN(IRIVER,8,icnt(iriver))
      WQ9 = WQDOWN(IRIVER,9,icnt(iriver))
      WQ10 = WQDOWN(IRIVER,10,icnt(iriver))
      WQ11 = WQDOWN(IRIVER,11,icnt(iriver))
      WQ12 = WQDOWN(IRIVER,12,icnt(iriver))
      WQ13 = WQDOWN(IRIVER,13,icnt(iriver))
      WQ14 = WQDOWN(IRIVER,14,icnt(iriver))
      WQ15 = WQDOWN(IRIVER,15,icnt(iriver))
      WQ16 = WQDOWN(IRIVER,16,icnt(iriver))
      WQ17 = WQDOWN(IRIVER,17,icnt(iriver))
      WQ18 = WQDOWN(IRIVER,18,icnt(iriver))
      WQ19 = WQDOWN(IRIVER,19,icnt(iriver))
      WQ20 = WQDOWN(IRIVER,20,icnt(iriver))
      WQ21 = WQDOWN(IRIVER,21,icnt(iriver))
      WQ22 = WQDOWN(IRIVER,22,icnt(iriver))
      WQ23 = WQDOWN(IRIVER,23,icnt(iriver))
      WQ24 = WQDOWN(IRIVER,24,icnt(iriver))
      WQ25 = WQDOWN(IRIVER,25,icnt(iriver))
      WQ26 = WQDOWN(IRIVER,26,icnt(iriver))
      WQ27 = WQDOWN(IRIVER,27,icnt(iriver))
      WQ28 = WQDOWN(IRIVER,28,icnt(iriver))
! for each layer determine the vol entrained and mix into infow

	   DO i=0,laycnt
			
!	Trap condition where plng_mxng returns impossibly high gamma, delQ
!	Use same logic as INFLOW and limit entrainment to 90% of a layer's volume
		   layvol=min(rivarea(i+1)*delQ,0.9d0*vol(ns-i))
		   vol(ns-i)=vol(ns-i)-layvol
	      S = COMBIN(S,Q,DI,SAL(ns-i),layvol,DEN(ns-i))
		   T = COMBIN(T,Q,DI,TEMP(ns-i),layvol,DEN(ns-i))		
!		write(88,*)'plng_mxng entrain i,laycnt,layvol,delq,q,t',i,laycnt,ns,vol(ns-i),temp(ns-i),layvol,delq,q,t
!		write(*,*)	
		   DI = DENSTY(T,S)
		   CCFF1 = COMBINV(CCFF1,Q,CF(ns-i,1),layvol)
		   CCFF2 = COMBINV(CCFF2,Q,CF(ns-i,2),layvol)
		   CCFF3 = COMBINV(CCFF3,Q,CF(ns-i,3),layvol)
		   CCFF4 = COMBINV(CCFF4,Q,CF(ns-i,4),layvol)
		   CCFF5 = COMBINV(CCFF5,Q,CF(ns-i,5),layvol)
		   CCFF6 = COMBINV(CCFF6,Q,CF(ns-i,6),layvol)
		   CCFF7 = COMBINV(CCFF7,Q,CF(ns-i,7),layvol)
		   WQ1 = COMBINV(WQ1,Q,WQUAL(ns-i,1),layvol)
		   WQ2 = COMBINV(WQ2,Q,WQUAL(ns-i,2),layvol)
		   WQ3 = COMBINV(WQ3,Q,WQUAL(ns-i,3),layvol)
		   WQ4 = COMBINV(WQ4,Q,WQUAL(ns-i,4),layvol)
		   WQ5 = COMBINV(WQ5,Q,WQUAL(ns-i,5),layvol)
		   WQ6 = COMBINV(WQ6,Q,WQUAL(ns-i,6),layvol)
		   WQ7 = COMBINV(WQ7,Q,WQUAL(ns-i,7),layvol)
		   WQ8 = COMBINV(WQ8,Q,WQUAL(ns-i,8),layvol)
		   WQ9 = COMBINV(WQ9,Q,WQUAL(ns-i,9),layvol)
		   WQ10 = COMBINV(WQ10,Q,WQUAL(ns-i,10),layvol)
		   WQ11 = COMBINV(WQ11,Q,WQUAL(ns-i,11),layvol)
		   WQ12 = COMBINV(WQ12,Q,WQUAL(ns-i,12),layvol)
		   WQ13 = COMBINV(WQ13,Q,WQUAL(ns-i,13),layvol)
		   WQ14 = COMBINV(WQ14,Q,WQUAL(ns-i,14),layvol)
		   WQ15 = COMBINV(WQ15,Q,WQUAL(ns-i,15),layvol)
		   WQ16 = COMBINV(WQ16,Q,WQUAL(ns-i,16),layvol)
		   WQ17 = COMBINV(WQ17,Q,WQUAL(ns-i,17),layvol)
		   WQ18 = COMBINV(WQ18,Q,WQUAL(ns-i,18),layvol)
		   WQ19 = COMBINV(WQ19,Q,WQUAL(ns-i,19),layvol)
		   WQ20 = COMBINV(WQ20,Q,WQUAL(ns-i,20),layvol)
		   WQ21 = COMBINV(WQ21,Q,WQUAL(ns-i,21),layvol)
		   WQ22 = COMBINV(WQ22,Q,WQUAL(ns-i,22),layvol)
		   WQ23 = COMBINV(WQ23,Q,WQUAL(ns-i,23),layvol)
		   WQ24 = COMBINV(WQ24,Q,WQUAL(ns-i,24),layvol)
		   WQ25 = COMBINV(WQ25,Q,WQUAL(ns-i,25),layvol)
		   WQ26 = COMBINV(WQ26,Q,WQUAL(ns-i,26),layvol)
		   WQ27 = COMBINV(WQ27,Q,WQUAL(ns-i,27),layvol)
		   WQ28 = COMBINV(WQ28,Q,WQUAL(ns-i,28),layvol)
		   Q = Q+layvol
	   ENDDO
!*********************************************************************
!wef  output the temp profile for the dat of the injection
!*********************************************************************
!	WRITE(111,444) ns, nosecs, iriver, lakecheck, Q_initial,T_initial,
!     &		Q,T, gammaf, mix_factor
!wef	IF(WQ3>1.0)WRITE(111,555) Q, gammaf
!	IF(WQ3>1.0)WRITE(111,555) (temp(i),depth(i), i=1,ns)
444	FORMAT(3I10,1x,L1,6f10.3)
555	FORMAT(2f10.3)
!*********************************************************************
! update the downflow stacks values
      QDOWN(IRIVER,icnt(iriver)) = Q
      TDOWN(IRIVER,icnt(iriver)) = T
      SDOWN(IRIVER,icnt(iriver)) = S
      CFDOWN(IRIVER,1,icnt(iriver)) = CCFF1
      CFDOWN(IRIVER,2,icnt(iriver)) = CCFF2
      CFDOWN(IRIVER,3,icnt(iriver)) = CCFF3
      CFDOWN(IRIVER,4,icnt(iriver)) = CCFF4
      CFDOWN(IRIVER,5,icnt(iriver)) = CCFF5
      CFDOWN(IRIVER,6,icnt(iriver)) = CCFF6
      CFDOWN(IRIVER,7,icnt(iriver)) = CCFF7
      WQDOWN(IRIVER,1,icnt(iriver)) = WQ1
      WQDOWN(IRIVER,2,icnt(iriver)) = WQ2
      WQDOWN(IRIVER,3,icnt(iriver)) = WQ3
      WQDOWN(IRIVER,4,icnt(iriver)) = WQ4
      WQDOWN(IRIVER,5,icnt(iriver)) = WQ5
      WQDOWN(IRIVER,6,icnt(iriver)) = WQ6
      WQDOWN(IRIVER,7,icnt(iriver)) = WQ7
      WQDOWN(IRIVER,8,icnt(iriver)) = WQ8
      WQDOWN(IRIVER,9,icnt(iriver)) = WQ9
      WQDOWN(IRIVER,10,icnt(iriver)) = WQ10
      WQDOWN(IRIVER,11,icnt(iriver)) = WQ11
      WQDOWN(IRIVER,12,icnt(iriver)) = WQ12
      WQDOWN(IRIVER,13,icnt(iriver)) = WQ13
      WQDOWN(IRIVER,14,icnt(iriver)) = WQ14
      WQDOWN(IRIVER,15,icnt(iriver)) = WQ15
      WQDOWN(IRIVER,16,icnt(iriver)) = WQ16
      WQDOWN(IRIVER,17,icnt(iriver)) = WQ17
      WQDOWN(IRIVER,18,icnt(iriver)) = WQ18
      WQDOWN(IRIVER,19,icnt(iriver)) = WQ19
      WQDOWN(IRIVER,20,icnt(iriver)) = WQ20
      WQDOWN(IRIVER,21,icnt(iriver)) = WQ21
      WQDOWN(IRIVER,22,icnt(iriver)) = WQ22
      WQDOWN(IRIVER,23,icnt(iriver)) = WQ23
      WQDOWN(IRIVER,24,icnt(iriver)) = WQ24
      WQDOWN(IRIVER,25,icnt(iriver)) = WQ25
      WQDOWN(IRIVER,26,icnt(iriver)) = WQ26
      WQDOWN(IRIVER,27,icnt(iriver)) = WQ27
      WQDOWN(IRIVER,28,icnt(iriver)) = WQ28

! pass the bed slope for the last specified segment, typically last river channel xsection to dam wall

	   upDpth = depth(ns)+base-BgnEle(iriver,max_seg_count-1)
	   dwnDpth = depth(ns)+base-BgnEle(iriver,max_seg_count)
	   Lngth = SegLngth(iriver,max_seg_count)
	   BSL = Slope(upDpth,dwnDpth,Lngth)
      Bedslope(iriver)=Slope(upDpth,dwnDpth,Lngth)
!      print*,iriver,Bedslope(iriver),ALPHA(IRIVER, numseg(iriver))
	   RETURN
	   END SUBROUTINE plng_mxng
!*********************************************************************	
!*********************************************************************	
	   REAL*8 FUNCTION DivAngle(upWdth,dwnWdth,Lngth)
!*********************************************************************
	   IMPLICIT NONE
	   REAL*8, parameter:: PI = 3.1415926535897932384626433832795d0
	   REAL*8:: upWdth,dwnWdth,Lngth

	   DivAngle=ATAN((dwnWdth-upWdth)/2/Lngth)*180.0d0/PI

	   DivAngle = MAX(DivAngle,0.0001d0)
	   RETURN
	   END FUNCTION DivAngle
!*********************************************************************	
!*********************************************************************	
	   REAL*8 FUNCTION Froude(Q,g_prime,upWdth,upArea)
	   IMPLICIT NONE
	   REAL*8,INTENT(IN):: Q,g_prime
	   REAL*8:: hyddpth,vel,upWdth,upArea

	   hyddpth=upArea/upWdth
	   vel=Q/upArea

	   Froude = vel/SQRT(g_prime*hyddpth)
	   RETURN
	   END FUNCTION Froude
!*********************************************************************	
!*********************************************************************	
	   REAL*8 FUNCTION Frd_A(Q,g_prime,upArea)
	   IMPLICIT NONE
	   REAL*8,INTENT(IN):: Q,g_prime
	   REAL*8:: dpth,vel,upArea

	   dpth=SQRT(upArea)
	   vel=Q/upArea

	   Frd_A = vel/SQRT(g_prime*dpth)
	   RETURN
	   END FUNCTION Frd_A
!*********************************************************************	
!*********************************************************************	
	   REAL*8 FUNCTION Slope(upDpth,dwnDpth,Lngth)
	   IMPLICIT NONE
	   REAL*8::upDpth,dwnDpth,Lngth

	   Slope=(dwnDpth-upDpth)/Lngth

	   IF (Slope<0.0) Slope=0.001d0
	   RETURN
	   END FUNCTION Slope
!*********************************************************************	
!*********************************************************************	
	   SUBROUTINE AkSt(iriver,rghnss,upDpth,Sl,delta,Fr,incdist,gammaf)
!
!	returns gammaf, dist as arguments and plngdpth(iriver) through a common block
!	computes the water column depth at the plunge point, d and the horizontal distance
!	to the plunge point, dist
!
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE		
	   REAL*8,INTENT(IN):: Sl,delta,Fr,upDpth,rghnss
	   INTEGER*4, INTENT(IN):: iriver
	   REAL*8, INTENT(INOUT):: incdist,gammaf
	   REAL*8:: Sone, Stwo, K, d
	   Sone=0.25d0 
	   Stwo=0.75d0

	   IF (Fr > 10.0d0) THEN
		   IF (delta >= 25.0d0) THEN		! half-angle delta > 25 -> free jet
			   gammaf = 0.5d0*Fr**4.0d0
		   ELSE						! half-angle delta < 25 -> wall jet
			   gammaf = 0.25d0*Fr**4.0d0
		   ENDIF
	   ELSE 							! Fr <= 10
		   IF (delta <= 7.0d0) THEN
			   gammaf = 0.456d0*Fr + 0.02d0*delta - 0.012d0*Sl + 1.0d0
		   ELSE
			   gammaf = 0.223d0*Fr + 0.008d0*delta + 1.0d0
		   END IF
	   ENDIF


! bss trap unrealistic gammaf, limit to max value of 2.1 based on Fleenor fig 5-7
!	gammaf = min(gammaf,5.4)


	   IF (Sl > 0.0067d0) THEN
		   K = 1.0d0/(2.0d0*gammaf)*(((1.0d0+gammaf)/2.0d0+Sone)	+ 	      &		! Fleenor Eq 2-28, Akiyama and Stefan 1984
     	       SQRT(((1.0d0+gammaf)/2.0d0+Sone)**2-4.0d0*Sone/gammaf))
		   d = upDpth*K*(1.0d0/Sone)**(1.0d0/3.0d0)*Fr**(2.0d0/3.0d0)*gammaf

	   ELSE
		   K=1./(2.0d0*gammaf)*(((1.0d0+gammaf)/2.0d0+Stwo*Sl/rghnss)+    &
     		   SQRT(((1.0d0+gammaf)/2.0d0+Stwo*Sl/rghnss)**2.0d0 -4.0d0*Stwo*Sl/(rghnss*gammaf)))

		   d = upDpth*K*(rghnss/(Stwo*Sl))**(1.0d0/3.0d0)*Fr**(2.0d0/3.0d0)*gammaf

	   ENDIF

	   plngdpth(iriver) = d
!	   plngdpth(iriver) = Depth(ns)  !gbs 2015/10/27
	   incdist=(d - upDpth)/Sl
!	plngdpth(iriver)=depth(ns)-10.0
	   IF (plngdpth(iriver) < 1.) THEN
!		write(*,*)'plngdpth trouble',iriver,upDpth,incdist,Sl,Fr,plngdpth(iriver)
	   ENDIF

	   RETURN
	   END SUBROUTINE AkSt
!*********************************************************************	
!*********************************************************************	
	   SUBROUTINE HandD(iriver,firstDpth,Sl,upArea,Q,g_prime,Fr,As,incdist,gammaf)
!***********************************************************************************
!	returns the plunge entrainment coefficient, gammaf and the distance to 
!	plunge from the surface (IF negatively buoyant) from the assumed
!	river/reservoir interface (where normdepth = water column depth)
!-------------------------------------------------------------------------------------
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
		SAVE
	   REAL*8, INTENT(IN):: Sl,Fr,As,Q,g_prime,upArea,firstDpth
	   INTEGER*4, INTENT(IN):: iriver
	   REAL*8, INTENT(OUT):: incdist,gammaf
	   REAL*8:: alphah,c

	   alphah = 0.1d0
	   c = -0.02d0

! plunge underflow
	
	   plngdpth(iriver)=(0.774d0*Fr+1.16d0)*(AS/Sl)**(-0.25d0)*firstdpth		! depth at plunge point as per Hauenstein & Dracos 1984
	   plngdpth(iriver)=MIN(plngdpth(iriver),depth(ns))
	
	   gammaf=(2.0d0*alphah*plngdpth(iriver))/SQRT(2.0d0*upArea*alphah*Sl)     
!	   write(88,*)'gppos',Q,alphah,firstDpth,Sl,plngdpth(iriver),depth(ns),gammaf
	
	   gammaf = max(gammaf,1.0d0)	! defeat negative entrainment, old fix for Hume
	   incdist = (plngdpth(iriver)-firstDpth)/Sl
	
	   RETURN
	   END SUBROUTINE HandD
!*********************************************************************	
!*********************************************************************	
	   SUBROUTINE findplungelayer(iriver,laycnt)
!**********************************************************************
!	returns the plunge entrainment coefficient, gammaf and the distance to either
!	separation from the bottom (IF positively buoyant) or
!	plunge from the surface (IF negatively buoyant) from the assumed
!	river/reservoir interface (where normdepth = water column depth)
!--------------------------------------------------------------------------
	   USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
		SAVE
	   INTEGER*4, INTENT(IN):: iriver
	   INTEGER*4, INTENT(OUT):: laycnt
	   INTEGER*4 :: i

	   laycnt=0
	   DO i=1,ns-1	  
		   IF ((depth(ns)-depth(ns-i)) < plngdpth(iriver)) THEN
			   laycnt=laycnt+1
		   ELSE
			   exit	! when leaving (depth(ns) - depth(ns - (laycnt + 1)) contains plngdpth
		   ENDIF
		
	   ENDDO
	  
	   RETURN
	   END SUBROUTINE findplungelayer
!*********************************************************************	
!*********************************************************************	
	   REAL*8 FUNCTION normal_flow_depth(rghnss,Q,half_angle,SL)
!************************************************************************
		USE DLMWQ_VARIABLES
	   IMPLICIT NONE
!
!	returns the normal flowing depth for a triangular cross-section channel
!	Q          = flow [m^3/s]
!	half_angle = channel cross-section half-angle [degrees]
!	bedslope   = longitudinal channel bed slop [dimensionless, m/m]
!	rghnss    = specified channel roughness coefficient
!
	   REAL*8,INTENT(IN):: rghnss,Q,half_angle,SL

	   REAL*8, parameter:: PI = 3.1415926535897932384626433832795d0
	   REAL:: deg2rad

	   deg2rad = PI/180.0d0

	   normal_flow_depth = (rghnss*Q/(TAN(deg2rad*half_angle)	&
     				*((SIN(deg2rad*half_angle))**(2.0/3.0))*sqrt(SL)))**(3.0/8.0)

	   RETURN
	   END FUNCTION normal_flow_depth
!*********************************************************************	
!*****************************************************************************************************	
	   SUBROUTINE find_normal_segment(iriver,Q,normdepth_segment,max_seg_count,firstDpth,firstLngth)
!*****************************************************************************************************
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
!
!	returns the river segment in river iriver that contains the normal flowing depth
!
	   INTEGER*4, INTENT(in):: iriver, max_seg_count
	   REAL*8, INTENT(in):: Q
	
	   INTEGER*4, INTENT(inout):: normdepth_segment
	   REAL*8, INTENT(out):: firstDpth, firstLngth
	
	   REAL*8:: normdpth
	   REAL*8:: upDpth, dwnDpth, Lngth, Sl	! segment characteristics
	   REAL*8:: rghnss
	
	   INTEGER*4:: segment

	   REAL*8:: Slope, normal_flow_depth		!declare functions
	
	   segment = normdepth_segment			! assume calling routine specifies where to start looking
										            ! normally initialised to 1
	   DO WHILE (segment < max_seg_count)
	
		   Lngth = SegLngth(iriver,segment+1)			! modified for consistency with lakeplunge, must refer to downstream input
		   upDpth = depth(ns)+base-BgnEle(iriver,segment)
		   dwnDpth = depth(ns)+base-BgnEle(iriver,segment+1)
		   Sl = Slope(upDpth,dwnDpth,Lngth)
		   rghnss = cdrag(iriver,segment)     
		   
		   normDpth = normal_flow_depth(rghnss,Q,alpha(iriver,segment),Sl)
				

		   IF (normDpth .le. dwnDpth) THEN			
			   IF (normDpth .gt. upDpth) THEN
				   firstDpth=normDpth
				   firstLngth = (dwnDpth-firstDpth)/Sl
			   ELSE
				   firstDpth = upDpth
				   firstLngth = Lngth
			   ENDIF			
			   normdepth_segment = segment + 1		! normal depth occurs in segment + 1
			   RETURN			
		   ENDIF	!(normDpth .le. dwnDpth)		
		   segment = segment+1
	
	   ENDDO
!
!	If we get here, the normal flowing depth is greater than the downstream depth of the last segment
!	So play it safe and assume the last segment ends at the deepest point of the storage, assumed
!	to be the dam

	   firstDpth = Depth(ns)
	   firstLngth = 0.0d0
	   normdepth_segment = max_seg_count
	
	   RETURN
	
	   END SUBROUTINE find_normal_segment
!*********************************************************************	
!*********************************************************************************	
	   SUBROUTINE find_plunge_segment(iriver,Q,normdepth_segment,        &
     		   max_seg_count,firstDpth,firstLngth,gammaf,mix_factor)
!********************************************************************************
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
		SAVE
!
!	returns the river segment in river iriver that contains the plunge depth
!
	   INTEGER*4, INTENT(in):: iriver, max_seg_count, normdepth_segment
	   REAL*8, INTENT(in):: firstDpth, firstLngth, Q
   	
	   REAL*8, INTENT(inout):: gammaf, mix_factor
   	
	   REAL*8:: rghnss, delta
	   REAL*8:: Slope, Froude, DivAngle		! declare functions
!IF	REAL*8:: DENSTY,Gprime,
	   REAL*8:: dist											! cumulative distance from normal flow depth to plungedepth
	   REAL*8:: upWdth, dwnWdth, Lngth, upDpth, dwnDpth		! channel segment properties
	   REAL*8:: Sl, upArea							! channel slope and boundary divergence angle
	   REAL*8:: As											! channel aspect ratio
	   REAL*8:: grav, deninf, Fr
	   REAL*8:: incdist										! distance to plunge point in current segment
	   INTEGER*4:: segment
	   LOGICAL:: gpnegcheck
! wef 29jan05  modified the parameter grammer	
! wef 29jan05	REAL, parameter:: delta_threshold = 35., max_gammaf = 10.
	   REAL*8:: delta_threshold, max_gammaf
	   delta_threshold = 35.0d0
	   max_gammaf = 10.0d0
! wef 29jan05  modified the parameter grammer	

	
	   dist = 0.0
	!delta_threshold = 45.				! this should be defined by comparing riverplunge and lake plunge intersection
	   segment = normdepth_segment			! assume calling routine specifies where to start looking
										! normally initialised to 1

	   DenInf = DENSTY(TEMINF(IRIVER),SALINF(IRIVER))	! check for overflow
	   grav=Gprime(den(ns),deninf)
	   IF (grav .LT. 0.0) THEN
		   gpnegcheck = .TRUE.		! overflow
		   grav = grav * (-1.0d0)
	   ELSE
	  	   gpnegcheck = .FALSE.
	   ENDIF
	

	   DO WHILE (segment <= max_seg_count)	
		   rghnss = cdrag(iriver,segment)
		   upWdth = BgnWdth(iriver,segment - 1)
		   dwnWdth = BgnWdth(iriver,segment)
		   Lngth = SegLngth(iriver,segment)
		   upDpth = depth(ns)+base-BgnEle(iriver,segment - 1)
		   dwnDpth = depth(ns)+base-BgnEle(iriver,segment)		
		   Sl = Slope(upDpth,dwnDpth,Lngth)
		   delta = DivAngle(upWdth,dwnWdth,Lngth)

		   IF (segment .eq. normdepth_segment) THEN		! the first segment checked contains the normal depth and can be
			   upDpth = firstDpth						! shorter than a full segment
			   Lngth = firstLngth
		   ENDIF

		   upArea = upWdth*upDpth		! can we get here IF firstDpth = Depth(ns), i.e. normdpth not found?
		
		   Fr = Froude(Q,grav,upWdth,upArea)
		   As = firstDpth/upWdth

!
!	trap a buoyant overflow and process following Safaie
!
		   IF (gpnegcheck .EQ. .TRUE.) THEN	! overflow uses Safaie for all channel configurations	
			   CALL Safaie(iriver,Q,grav,upArea,As,Fr,Sl,firstDpth,incdist,gammaf)
			   mix_factor = delta/90.0d0			
		   ELSE
!	Here for underflows!
			   IF (delta < delta_threshold) THEN		! delta_threshold defines boundary between river and lake plunge algorithms				
				   CALL AkSt(iriver,rghnss,upDpth,Sl,delta,Fr,incdist,gammaf)		! it's a river
				   mix_factor = 1.0d0			
			   ELSE				
				   CALL HandD(iriver,upDpth,Sl,upArea,Q,grav,Fr,As,incdist,gammaf)						! it's a lake
				   mix_factor = delta/90.				
			   ENDIF			
		   ENDIF !gpnegcheck		  	
		   IF(Fr > 10.0d0  .or. gammaf > 10.0d0) THEN
!			   write(88,*)
!			   write(88,*)'inflow out of bounds',iriver, segment, Fr,dist, gammaf, delta
!			   write(88,*)
		   ENDIF
!wef 14dec03  set maximum gamma to a limit since excessive values cause problems in resint
		   gammaf = min(max_gammaf, gammaf)
						
		   IF (incdist <= Lngth) THEN
			   incdist = MAX(incdist,0.0d0)
			   dist = dist + Lngth									! note that dist is never referenced and can be deleted
			   plngdpth(iriver) = min(dwnDpth,plngdpth(iriver))	! AkSt may set this incorrectly deep
!		print*,dwnDpth, plngdpth(iriver)
			   RETURN
		   ENDIF		
		   segment = segment+1
		   dist = dist + Lngth						! compute aggregate distance						
	   ENDDO		! segment < max_seg_count
	
!
!	If we get here, the plunge depth is greater than the downstream depth of the last segment
!	So play it safe and assume the last segment ends at the deepest point of the storage, assumed
!	to be the dam

	   plngdpth(iriver) = Depth(ns)	
	   RETURN	
	   END SUBROUTINE find_plunge_segment
!*********************************************************************	
!*********************************************************************	
	   SUBROUTINE Safaie(iriver,Q,g_prime,upArea,As,Fr,Sl,firstDpth,incdist,gammaf)
      USE DLMWQ_VARIABLES
      USE DLMWQ_SERVICES
		IMPLICIT NONE
		SAVE
	
	   REAL*8, INTENT(IN):: Sl,Fr,As,Q,g_prime,upArea,firstDpth
	   REAL*8, INTENT(out):: gammaf, incdist
	
	   REAL*8::	alph,FrA,Frd_A	
	   REAL*8::	hratio		! ratio of separation height to upstream normal flowing depth
	   INTEGER*4, INTENT(IN):: iriver
	
		FrA = Frd_A(Q,g_prime,upArea)
		plngdpth(iriver) = 0.914d0*SQRT(Fr)*firstDpth		! depth at separation from bottom as per Safaie 1979
		
		IF (plngdpth(iriver) > depth(ns)) THEN
!			write(88,*)'plngdpth trouble for overflow'	! delete this after debugging
		ENDIF
		
		plngdpth(iriver)=MIN(plngdpth(iriver),depth(ns))

		hratio = plngdpth(iriver)/firstDpth
! 
!	test a more straightforward implementation of Safaie's formula
!	Fleenor's Eq 2-16

		alph=(FrA**2-6.25d0)/(5.22d0*Fr-6.25d0)
		
		alph = max(alph,0.0d0)			! defeat negative entrainment this happens sometimes when
									! either FrA < 2.5 or Fr < 1.2
		
		gammaf = sqrt(1.0d0 + alph*((hratio*hratio) - 1.0d0))

!		alph=(Sl/As)*((FrA**2-6.25)/(5.22*Fr-6.25))
!		gammaf=SQRT(Q**2+(2*alph*Q**2*firstDpth**2/upArea)/
!     &			(c*COS(ATAN(Sl)))*(EXP(c*COS(ATAN(Sl))*(plngdpth(iriver)+firstDpth/Sl)/(Sl*firstDpth))-1.))/Q
		
!		write(88,*)'Safaie',Q,alph,firstDpth,Sl,plngdpth(iriver),
!     &							depth(ns),gammaf
		gammaf = max(gammaf,1.0d0)	! defeat negative entrainment, old fix for Hume
		incdist=(plngdpth(iriver)-firstDpth)/Sl

	RETURN
	END SUBROUTINE Safaie
!*********************************************************************	
!*********************************************************************	
