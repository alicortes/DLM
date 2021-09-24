      MODULE DLMWQ_SERVICES
!*****************************************************
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  
      CONTAINS
!***********************************************************************
      REAL*8 FUNCTION densty(TEMPR,SALINT)
!**********************************************************************
!*
!* CALCULATES THE sigma-T OF WATER AT TEMPR DEGS C AND SALINITY
!* SALINT PPM NACL. after Millero and Poison (1981)
!*
!* note: this function returns (density - 1000)
!*       C1 accounts for this correction (=999.842594 
!*       IF we wanted density instead of sigma-T)
!*
!**********************************************************************
!
! New definitions double precision
!
!		USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE 
 
	   REAL(8) ::  dpure						   
      REAL(8) ::  TEMPR
      REAL(8) ::  SALINT    
      REAL(8) ::  den_param_T1,den_param_T2,den_param_T3,den_param_T4,den_param_T5
      REAL(8) ::  den_param_S1,den_param_s32,den_param_s2
      REAL(8) ::  csal1
      REAL(8) ::  csal32
      REAL(8) ::  csal2

      REAL(8), PARAMETER :: C1=0.157406d0,C2=6.793952d-2,C3=9.095290d-3
	   REAL(8), PARAMETER :: C4=1.001685d-4,C5=1.120083d-6,C6=6.536332d-9				  
      REAL(8), PARAMETER :: D1=8.24493d-1,D2=4.0899d-3,D3=7.6438d-5
	   REAL(8), PARAMETER :: D4=8.2467d-7, D5=5.3875d-9,D6=5.72466d-3
	   REAL(8), PARAMETER :: D7=1.0227d-4, D8=1.6546d-6,D9=4.8314d-4

      REAL(8), PARAMETER ::  onept5=1.5d0,thsnd=1000.0d0, tenthou=1.0d5
      REAL(8) ::  term(15)
											  
	   INTEGER*4 tm

!   TAKEN FROM MOD2 DENSITY FUNCTION
      den_param_t1	=	dble(tempr)
      den_param_s1	=	dble(salint)/thsnd
      tm	=	nint(den_param_T1*tenthou)
      den_param_t1	=	dble(tm)/tenthou
      den_param_t2	=	den_param_t1*den_param_t1
      den_param_t3	=	den_param_t2*den_param_t1
      den_param_t4	=	den_param_t3*den_param_t1
      den_param_t5	=	den_param_t4*den_param_t1
      den_param_s2	=	den_param_s1*den_param_s1
      den_param_s32	=	den_param_s1**onept5

      term(1)  =	-c1
      term(2)  =	c2*den_param_t1
      term(3)  =	-c3*den_param_t2
      term(4)  =	c4*den_param_t3
      term(5)  =	-c5*den_param_t4
      term(6)  =	c6*den_param_t5
      term(7)  =	d1
      term(8)  =	-d2*den_param_t1
      term(9)  =	d3*den_param_t2
      term(10) =	-d4*den_param_t3
      term(11) =	d5*den_param_t4
      term(12) =	-d6
      term(13) =	d7*den_param_t1
      term(14) =	-d8*den_param_t2
      term(15) =	d9
      
	   dpure	= term(6)+term(5)+term(4)+term(2)+term(3)+term(1)
      csal1	= (term(11)+term(10)+term(9)+term(8)+term(7))*den_param_s1
      csal32	= (term(14)+term(13)+term(12))*den_param_s32
      csal2	= term(15)*den_param_s2
      
	   densty = REAL(dpure+csal1+csal32+csal2)
	
      RETURN
      END FUNCTION DENSTY
!*******************************************************
      REAL*8 FUNCTION vaporTVA(Tem,rh)
!************************************************
! Vapor pressure at the Air temperature TVA (1972)
!--------------------------------------------------	
!		USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE 
 
	   REAL*8 tem,rh,a_tva, b, c
	   PARAMETER (a_tva=7.5, b=237.3, c=0.7858)
      vaporTVA = rh*exp(2.303*(c+(a_tva*tem/(b+tem))))
	   RETURN
	   END FUNCTION vaporTVA
!*********************************************************************
      REAL*8 FUNCTION SATVAP(T)
!***********************************************************************
! CALCULATES THE SATURATED VAPOUR PRESSURE (SATVAP) CORRESPONDING TO
! TEMPERATURE T. T = temperature (degrees C)
!                TH = temperature (degrees k)
!
!-----------------------------------------------------------------------
		USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE 

      REAL*8 AL
      REAL*8 T
      REAL*8 TH

      REAL*8 dCtodK,C1,C2,thsnd,ten

      PARAMETER (dCtodK=273.15d0,C1=9.28603523d0,C2=2.32237885d0,   &
                thsnd=1000.0d0,ten=10.0d0)

      TH=T+dCtodK
      AL=C1-C2*thsnd/TH
      SATVAP=ten**AL

      RETURN
      END FUNCTION SATVAP
 
!*******************************************************************
      REAL*8 FUNCTION SQR(x)
!*******************************************************************
!
!  This function computes the sqare of the argument
!*******************************************************************!
		USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE 

      REAL*8 x

      SQR=x*x
      RETURN
      END FUNCTION SQR
!************************************************************************
      REAL*8 FUNCTION COMBIN(C1,V1,D1,C2,V2,D2)
!*************************************************************************
! This function combines two layers and return the mean concentration
!------------------------------------------------------------------------------
		USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE 

      REAL*8 thsnd
      PARAMETER(thsnd=1000.0d0)

      REAL*8 C1
      REAL*8 C2
      REAL*8 D1
      REAL*8 D2
      REAL*8 V1
      REAL*8 V2

      REAL*8 M1big
      REAL*8 M1sml
      REAL*8 M2big
      REAL*8 M2sml
      REAL*8 MTOTAL

      M1big = V1*thsnd
      M1sml = V1*D1
      M2big = V2*thsnd
      M2sml = V2*D2

      MTOTAL=(M1sml+M2sml)+M1big+M2big

      COMBIN=(C1*M1sml+C2*M2sml)/MTOTAL+thsnd*(C1*V1+C2*V2)/MTOTAL
	
      RETURN
      END FUNCTION COMBIN
!***************************************************************************
     REAL*8 FUNCTION COMBINV(c1,v1,c2,v2)
!***************************************************************************
!  Function to combine two layers and return the mean
!  concentration taking into account volumes
		USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE 

      REAL*8 C1,C2,V1,V2
      COMBINV=(C1*V1+C2*V2)/(V1+V2)

      RETURN
      END FUNCTION COMBINV
!*****************************************************************************
      REAL*8 FUNCTION GPRIME (D1,D2)
!*****************************************************************************
! This function calculate the reduced gravity given two sigma-T values
!-----------------------------------------------------------------------------
		USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE 

      REAL*8 thsnd,two,g
      PARAMETER(thsnd=1000.0d0,two=2.0d0,g=9.81d0)

      REAL*8 D1
      REAL*8 D2

      GPRIME=(D2-D1)*g*two/((D1+D2)+two*thsnd)

      RETURN
      END FUNCTION GPRIME

!******************************************************************************
      SUBROUTINE DIFUSE 
!******************************************************************************
!      modified by b.s. sherman 18 april 1986
!      version 6.4 single precision b.s. sherman 20 december1986
!      this version uses proper concentration units
!
!     version 6.5 G. Hocking - uses new explicit scheme for the
!     Diffusion calculations - see notes.
!--------------------------------------------------------------------------------
!    
!      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  

      REAL(8), PARAMETER :: zero=0.0d0,one=1.0d0,two=2.0d0
	   REAL(8), PARAMETER :: pt001=0.001d0,pt6=0.6d0
	   REAL(8), PARAMETER :: exchka=80.0d0
	   REAL(8), PARAMETER :: epsmin=1d-5
      REAL(8), PARAMETER :: secqhr=900.0d0,hrday=24.0d0

      REAL(8) ::  c(MAXNS,MAXDIF)
	   REAL(8) ::  caca_part_I(7),caca_part_F(7)
!IF      REAL(8) ::  densty
      REAL(8) ::  eps(MAXNS)
      REAL(8) ::  expon
	   REAL(8) ::  exchk2
!IF      REAL(8) ::  gprime
	   REAL(8) ::  Nitrogen_T_1,Nitrogen_T,Part_T_1(7)
!IF      REAL(8) ::  SQR
      REAL(8) ::  TM3
      REAL(8) ::  XNSQ

      INTEGER*4 i,j,i_p,j_p
      INTEGER*4 nz
	   INTEGER*4 jday

      LOGICAL*4 flag


      exchk2=-exchka
!
! Set the time step
!
	   TM3 = NOSECS
!
! Set up concentration array
!
      CALL DUBLUP(C)
!
! Get out of here IF fully mixed, remember ns is changed in STABLE
! Note: the above statement is no long true as STABLE has been modified
!       to retain the layer structure.
!
10    CONTINUE

      nz=ns-1
      IF (nz.eq.0) RETURN
	
!
! Calculate eddy diffusivities
!
      FLAG = .FALSE.
!
! Look for two layer structure
!
	   IF (nz .eq. 1) THEN
	      FLAG = .TRUE.
	      GOTO 29
      ENDIF

	   DO 28 i=2, nz
	      IF (den(i)-den(i-1) .gt. epsmin) GOTO 29
28    CONTINUE

	   FLAG = .TRUE.
! MB 29 - 30 TAKEN FROM MOD2
! Smoothing NSQR distribution

29    DO 30 i=3, nz-1
	      XNSQ=GPRIME(den(i+2),den(i-2))/(DEPTHM(i+2)-DEPTHM(i-2))
	      IF (XNSQ .le. 1.d-6) THEN
	         XNSQ=zero 
	      ENDIF

	      IF (XNSQ.eq.0.0.OR.VEL.eq.0.0.OR.WNSQ.lt.0.0)THEN
	         ep(i)=ZERO
	         GOTO 30
	      ENDIF
	  
	      IF (XNSQ .lt. 1.D-6 .AND. VEL .lt. 1.D-6 .AND.WNSQ .lt. 1.D-6) THEN
	         ep(i)=ZERO
	         GOTO 30
	      ELSE
            ep(i)=(DISS/two)/(XNSQ+pt6*VEL*VEL*WNSQ)
	      ENDIF

	      IF (FLAG .AND. i .eq. nz) GOTO 30
	      IF (depth(i) .gt. XMOM) GOTO 30
	      IF(HSIG .le. zero) THEN
	         ep(i)=zero
	         GOTO 30
	      ENDIF

	      EXPON=(-one*SQR(depth(ns)-H1-depth(i)))/HSIG

	      IF(EXPON .lt. exchk2) THEN
	         ep(i)=zero
	      ELSE
	         ep(i)=ep(i)*(EXP(EXPON)+1.D-7)
	      ENDIF
30    CONTINUE

      IF(nz.eq.1)THEN
	      ep(nz)= 0.0d0
      ELSE
	      ep(nz)= ep(nz-1)
      ENDIF
      ep(ns) = ep(nz)
      ep(2)  = ep(3)
      ep(1)  = ep(2)
!
! Process each diffusable species
!
      DO 200 j=1,NUMDIF
	      DO 100 i=1,ns
		      eps(i)=(ep(i)+diff(j))
100      CONTINUE
	      CALL DIFCAL(TM3,eps,C(1,j),j)
200   CONTINUE

! Because the effect of DIFCAL is slightly different for components
! of oxygen relative to its total, adjust:

      DO 250 i=1,ns
!2015/08  C(i,10)=C(i,6) +C(i,7)+C(i,8)+C(i,9)
          C(i,10)=C(i,10) 
250   CONTINUE
!     print*,'C(ns,6)=',C(ns,6)    !C(i,6)=wqual(i,4)---
!	   print*,'C(ns,7)=',C(ns,7)	 !C(i,7)=wqual(i,5)--- surface
!	   print*,'C(ns,8)=',C(ns,8)	 !C(i,8)=wqual(i,6)---
!	   print*,'C(ns,9)=',C(ns,9)	 !C(i,9)=wqual(i,7)---source and sink terms
!	   print*,'C(ns,10)=',C(ns,10)	 !C(i,10)=wqual(i,8)DO
!     C(i,11)=wqual(i,9)BOD

! Get temps and salinities
! Calculate the new density structure and check for instabilities
!
      CALL CHOP(C)

      DO 300 i=1,ns
	      den(i)=densty(temp(i),sal(i))
300   CONTINUE
!
! Call STABLE to perform the mixing
!
      DO 400 i=2,ns
	      IF (den(i) .ge. den(i-1)) THEN
		      CALL STABLE
		      GOTO 500
	      ENDIF
400   CONTINUE
!
! Add heat diffusivity to ep for reference by inflow and outflo
!
500   DO 600 i=1,ns
	      ep(i)=ep(i)+DIFF(1)
600   CONTINUE

      RETURN
      END SUBROUTINE DIFUSE

!************************************************************************
      SUBROUTINE DUBLUP(C)
!**********************************************************************
! This subroutine creates an array of the concentrations
! of diffusable species
! test version - concentration array has units kg nacl/m**3 for salt
!                and kg C/m**3 for temp
! INCLUDE THE WATER QUALITY VARIABLES
!**********************************************************************
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  
      
      REAL(8)		:: thsnd=1000.0d0
      REAL(8)		:: C(MAXNS,MAXDIF)

      INTEGER(4)	:: i, j

      DO 100 j = 1,numdif
	      DO 100 i = 1,ns
		      IF(j .eq. 1) C(i,j)=temp(i)               !Temperature
		      IF(j .eq. 2) C(i,j)=sal(i)*(den(i)+thsnd) !Salinity
	         IF(j .eq. 3) C(i,j)=wqual(i,1)            !Chla a
		      IF(j .eq. 4) C(i,j)=wqual(i,2)            !Chla a
		      IF(j .eq. 5) C(i,j)=wqual(i,3)            !Chla a
		      IF(j .eq. 6) C(i,j)=wqual(i,4)            !Chla a
	         IF(j .eq. 7) C(i,j)=wqual(i,5)            !Chla a
		      IF(j .eq. 8) C(i,j)=wqual(i,6)            !Chla a
	         IF(j .eq. 9) C(i,j)=wqual(i,7)            !Chla a
	         IF(j .eq. 10) C(i,j)=wqual(i,8)           !DO
		      IF(j .eq. 11) C(i,j)=wqual(i,9)           !BOD
	         IF(j .eq. 12) C(i,j)=wqual(i,10)
		      IF(j .eq. 13) C(i,j)=wqual(i,11)
		      IF(j .eq. 14) C(i,j)=wqual(i,12)
		      IF(j .eq. 15) C(i,j)=wqual(i,13)
	         IF(j .eq. 16) C(i,j)=wqual(i,14)
		      IF(j .eq. 17) C(i,j)=wqual(i,15)
	         IF(j .eq. 18) C(i,j)=wqual(i,16)
		      IF(j .eq. 19) C(i,j)=wqual(i,17)
	         IF(j .eq. 20) C(i,j)=wqual(i,18)
	         IF(j .eq. 21) C(i,j)=wqual(i,19)
	         IF(j .eq. 22) C(i,j)=wqual(i,20)
	         IF(j .eq. 23) C(i,j)=wqual(i,21)
	         IF(j .eq. 24) C(i,j)=wqual(i,22)
	         IF(j .eq. 25) C(i,j)=wqual(i,23)
	         IF(j .eq. 26) C(i,j)=wqual(i,24)
	         IF(j .eq. 27) C(i,j)=wqual(i,25)
	         IF(j .eq. 28) C(i,j)=wqual(i,26)
	         IF(j .eq. 29) C(i,j)=wqual(i,27)
	         IF(j .eq. 30) C(i,j)=wqual(i,28)
		      IF(j .eq. 31) C(i,j)=cf(i,1)
		      IF(j .eq. 32) C(i,j)=cf(i,2)
		      IF(j .eq. 33) C(i,j)=cf(i,3)
		      IF(j .eq. 34) C(i,j)=cf(i,4)
		      IF(j .eq. 35) C(i,j)=cf(i,5)
		      IF(j .eq. 36) C(i,j)=cf(i,6)
		      IF(j .eq. 37) C(i,j)=cf(i,7)
100   CONTINUE
      RETURN
      END SUBROUTINE DUBLUP
!***********************************************************************
      SUBROUTINE CHOP(C)
!**********************************************************************
! This subroutine creates an array of the concentrations
! of diffusable species
! INCLUDE THE WATER QUALITY VARIABLES
!-------------------------------------------------------------------------C
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  
      
      INTEGER(4)	:: i,j	
      REAL(8),PARAMETER	:: thsnd = 1000.0d0
      REAL(8)				:: C(MAXNS,MAXDIF)

      DO 100 j = 1,numdif
	      DO 100 i = 1,ns
		      IF(j .eq. 1)  temp(i)=C(i,j)
	         IF(j .eq. 2)  sal(i)=C(i,j)/(den(i)+thsnd)
	         IF(j .eq. 3)  wqual(i,1) = C(i,j)
	         IF(j .eq. 4)  wqual(i,2) = C(i,j)
	         IF(j .eq. 5)  wqual(i,3) = C(i,j)
	         IF(j .eq. 6)  wqual(i,4) = C(i,j)
	         IF(j .eq. 7)  wqual(i,5) = C(i,j)
	         IF(j .eq. 8)  wqual(i,6) = C(i,j)
	         IF(j .eq. 9)  wqual(i,7) = C(i,j)
	         IF(j .eq. 10) wqual(i,8) = C(i,j)
	         IF(j .eq. 11) wqual(i,9) = C(i,j)
		      IF(j .eq. 12) wqual(i,10)= C(i,j)
		      IF(j .eq. 13) wqual(i,11)= C(i,j)
		      IF(j .eq. 14) wqual(i,12)= C(i,j)
		      IF(j .eq. 15) wqual(i,13)= C(i,j)
		      IF(j .eq. 16) wqual(i,14)= C(i,j)
		      IF(j .eq. 17) wqual(i,15)= C(i,j)
		      IF(j .eq. 18) wqual(i,16)= C(i,j)
		      IF(j .eq. 19) wqual(i,17)= C(i,j)
		      IF(j .eq. 20) wqual(i,18)= C(i,j)
		      IF(j .eq. 21) wqual(i,19)= C(i,j)
		      IF(j .eq. 22) wqual(i,20)= C(i,j)
		      IF(j .eq. 23) wqual(i,21)= C(i,j)
		      IF(j .eq. 24) wqual(i,22)= C(i,j)
		      IF(j .eq. 25) wqual(i,23)= C(i,j)
		      IF(j .eq. 26) wqual(i,24)= C(i,j)
		      IF(j .eq. 27) wqual(i,25)= C(i,j)
		      IF(j .eq. 28) wqual(i,26)= C(i,j)
		      IF(j .eq. 29) wqual(i,27)= C(i,j)
		      IF(j .eq. 30) wqual(i,28)= C(i,j)
		      IF(j .eq. 31) cf(i,1)=C(i,j)
		      IF(j .eq. 32) cf(i,2)=C(i,j)
		      IF(j .eq. 33) cf(i,3)=C(i,j)
		      IF(j .eq. 34) cf(i,4)=C(i,j)
		      IF(j .eq. 35) cf(i,5)=C(i,j)
		      IF(j .eq. 36) cf(i,6)=C(i,j)
		      IF(j .eq. 37) cf(i,7)=C(i,j)
100   CONTINUE
      RETURN
      END SUBROUTINE CHOP
!***********************************************************************************
      SUBROUTINE DIFCAL(TM3,EPS,C,N_Dif)
!***********************************************************************
! CALCULATES THE DIFFUSION IN A SINGLE TIME STEP
! in this version the diffusion is handled explicitly between each set
! of two layers. a decay term is calculated and only one sweep of the
! concentration array is made.
!
! This version uses proper concentration units
!
!***********************************************************************
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  

      REAL(8), PARAMETER :: zero  = 0.0d0,two=2.0d0,twenty=20.0d0
	   REAL(8), PARAMETER :: sixty = 60.0d0,thsnd=1000.0d0,tenM8=1.d-8
	   REAL(8), PARAMETER :: tenM12=1.0d-12
	   REAL(8), PARAMETER :: pi=3.1415926535897932384626433832795d0

      REAL(8) ::  EPS(MAXNS)
      REAL(8) ::  C(MAXNS)
      REAL(8) ::  TOL(MAXDIF)

      REAL(8) ::  CBAR
      REAL(8) ::  DELTAC
      REAL(8) ::  DEP
      REAL(8) ::  D1
      REAL(8) ::  D2
      REAL(8) ::  EXFAC
      REAL(8) ::  FACTOR
      REAL(8) ::  TM3, TM4
      REAL(8) ::  TOTC

      INTEGER*4 i
      INTEGER*4 IBOT
      INTEGER*4 IA
      INTEGER*4 IB
      INTEGER*4 IC
      INTEGER*4 k
      INTEGER*4 nz
	   INTEGER*4 N_dif

      REAL(8) ::    caca_part_I,caca_part_F
      INTEGER*4 j_p

      TOL(1)=1.0E-8
      TOL(2)=1.0
      TM4 = TM3/two
      IBOT = 1
      nz=ns-1
!
! Perform sweeps in opposite directions.
!
      DO 30 k = 1,2
	      IA = IBOT
	      IB = nz
	      IC = 1
	      IF (k.eq.2) THEN
	         IA = nz
	         IB = IBOT
	         IC = -1
	      ENDIF
 
 	      DO 30 i=IA,IB,IC
	         IF (ABS(C(i)).lt.1d-20) C(i) = 0.0d0            
	         IF (ABS(C(i+1)).lt.1d-20) C(i+1) = 0.0d0                
	
	         IF(ABS(C(i)-C(i+1)).LE.TOL(1))GOTO 30
	         TOTC = C(i)*vol(i)+C(i+1)*vol(i+1)
	         DEP = zero

	         IF (i .ne. 1) DEP = depth(i-1)
	         D1 = depth(i)-DEP
	         D2 = depth(i+1)-depth(i)
	         CBAR = (C(i)*D1+C(i+1)*D2)/(D2+D1)
	         DELTAC = C(i)-C(i+1)	
	         FACTOR = (EPS(i)+EPS(i+1))/two*(two/(D2+D1))**2*TM4

	         IF (FACTOR .gt. twenty) THEN
		         EXFAC = zero
	         ELSE
	            EXFAC = EXP(-FACTOR)
	         ENDIF

	         C(i)   = CBAR+EXFAC*D2*DELTAC/(D2+D1)
	         C(i+1) = CBAR-EXFAC*D1*DELTAC/(D2+D1)

	         IF ((vol(i)/D1-vol(i+1)/D2) .gt. zero) THEN
		         C(i)   = (TOTC-C(i+1)*vol(i+1))/vol(i)
	         ELSE
		         C(i+1) = (TOTC-C(i)*vol(i))/vol(i+1)
	         ENDIF
30    CONTINUE 
      RETURN
      END SUBROUTINE DIFCAL
!******************************************************************
      SUBROUTINE ENER
!******************************************************************
!  Calculates dissipation due to wind and inflow energy inputs
!------------------------------------------------------------------
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  
      
      REAL(8), PARAMETER :: zero=0.0d0,two=2.0d0,pt5=0.5d0
	   REAL(8), PARAMETER :: arfac=1.0D+6,g=9.81d0,secday=86400.0d0
	   REAL(8), PARAMETER :: thsnd=1000.0d0,cwnsq1=12.4d0

!      PARAMETER(zero=0.0d0,two=2.0d0,pt5=0.5d0,arfac=1.0D+6,
!     +g=9.81d0,secday=86400.0d0,thsnd=1000.0d0,cwnsq1=12.4d0)

      REAL(8) :: BFSQ(MAXNS)
      REAL(8) ::  CMsml
      REAL(8) ::  CMMIX
      REAL(8) ::  DBAR
      REAL(8) ::  DFLOC
      REAL(8) ::  DH
      REAL(8) ::  DSIG
      REAL(8) ::  EINFW
      REAL(8) ::  EW
      REAL(8) ::  EWW
!IF      REAL(8) ::  GPRIME
      REAL(8) ::  H0
      REAL(8) ::  PERES
      REAL(8) ::  Sbig
      REAL(8) ::  Ssml
      REAL(8) ::  SMLOC
      REAL(8) ::  SIGBT
      REAL(8) ::  SIGTP
      REAL(8) ::  Tbig
      REAL(8) ::  Tsml
      REAL(8) ::  TMLOC
      REAL(8) ::  USTAR
      REAL(8) ::  VMLOC
      REAL(8) ::  VMbig
      REAL(8) ::  VMsml
      REAL(8) ::  VOLMX
      REAL(8) ::  VTILDA
      REAL(8) ::  WP 
      REAL(8) ::  XMOM0
!     REAL(8) ::  XWIND
      REAL(8) ::  XZI
	   REAL(8) ::  c_drag

      INTEGER*4 i
      INTEGER*4 KL
	   INTEGER*4 ucheck
!
! Put in wind factor
!
!	U6X = U6*XWIND
	   WP=1.9344*0.24*U6X*U6X*U6X*1.E-6
!	USTAR = SQRT(1.612E-6*U6X**2)
! New Linear definition of Cd (Smith and Banke formulae, from Henderson-Sellers)
!
      c_drag = 1.3E-3 
!      c_drag = (0.61 + 0.063 * u6x)/ 1000.
!	ucheck = 1000.0 * u6x	
!	SELECT CASE (ucheck)
!		CASE (:300)
!			c_drag = .0015			
!		CASE (301:2200)
!			c_drag = (1.08 * u6x**(-0.15))/1000.0		
!		CASE (2201:5000)
!			c_drag = (0.771 + 0.0858 * u6x)/1000.0		
!		CASE (5001:8000)
!			c_drag= (0.867 + 0.0667 * u6x)/1000.0		
!		CASE (8001:25000)
!			c_drag= (1.2 + 0.025 * u6x)/1000.0	
!		CASE (25001:50000)
!			c_drag= (0.073 * u6x)/1000.	
!		CASE default
!			c_drag= 0.0037	
!	END SELECT
	   ustar = sqrt(c_drag*(1.24/1000.)*u6x*u6x)

!
! Find center of buoyancy; XMOM is 1st moment, XMOM0 is 0th moment
!
      DO 50 i=1,MAXNS
	      BFSQ(i) = zero
50    CONTINUE

      XMOM=zero
      XMOM0=zero
      DEPTHM(1)=depth(1)/two
!       WRITE(*,*)'ns =',ns
      DO 100 i=2,ns
	      DEPTHM(i)=(depth(i)+depth(i-1))/two
	      BFSQ(i)=GPRIME(den(i),den(i-1))/(DEPTHM(i)-DEPTHM(i-1))
	      XMOM=XMOM+depth(i-1)*BFSQ(i)*(vol(i)+vol(i-1))/two
	      XMOM0=XMOM0+BFSQ(i)*(vol(i)+vol(i-1))/two
!       WRITE(*,*)XMOM,XMOM0,BFSQ(i),vol(i),vol(i-1),den(i),den(i-1)
!       WRITE(*,*)DEPTHM(i),DEPTHM(i-1)
!       PAUSE

100   CONTINUE
!
! Define length scales, XMOM is the center of buoyancy ( m above bottom)
!                       H1 is the thickness of the upper mixed layer
!                       H0 is the thickness of the surface layer
!
	   IF (XMOM.eq.0.) GOTO 105
         XMOM=XMOM/XMOM0
105   CONTINUE
      H1=depth(ns)-XMOM
      H0=depth(ns)-depth(ns-1)
!
! Calculate mean reservoir properties:
!                        SM  ... mean salinity
!                        TM  ... mean temperature
!                        CMMIX  ... center of mass of fully mixed reservoir
!
      CMsml=zero
      CMMIX=zero
      Tbig=zero
      Tsml=zero
      Sbig=zero
      Ssml=zero
      VMbig=zero
      VMsml=zero
      CMsml=zero
!
! Calculate first moments
!
      DO 200 i=1,ns
	      CALL ADDLAY(VMbig,VMsml,Tbig,Tsml,Sbig,Ssml,VMLOC,TMLOC,SMLOC,DFLOC,i)
	      CMsml=CMsml+DEPTHM(i)*den(i)*vol(i)
	      CMMIX=CMMIX+DEPTHM(i)*vol(i)
200   CONTINUE 
!
! define mean properties and PERES, the potential energy of the reservoir
! due to the stratification; this calculation has no effect on the model
! DH is actually the quantity DH*VM
!
      DH = CMMIX*DF-CMsml
      PERES = g*DH*thsnd
!
!     find layer containing the center of buoyancy   (layer i)
!
      DO 300 i=1,ns
	      IF(depth(i).gt.XMOM) GOTO 310
300   CONTINUE
      i=ns
!
! Calculate the variance of the buoyancy distribution about XMOM
!               SIGBT ... 0th moment of buoyancy about XMOM
!               SIGTP ... 2nd moment of buoyancy about XMOM
!               SIGTP/SIGBT ... variance of buoyancy distribution
!               DSIG ... standard deviation of buoyancy distribution
!
310   CONTINUE
      SIGBT=zero
      SIGTP=zero
      HSIG=XMOM**two
      IF(i.eq.1) GOTO 360
      DO 350 KL=i,2,-1
	      IF(BFSQ(KL).LE.zero) BFSQ(KL)=zero
	      XZI=DEPTHM(i)-DEPTHM(KL-1)
!
! The '2' in the equation for SIGTP implies a symmetrical distribution
!
	      SIGTP=SIGTP+two*XZI**two*BFSQ(KL)*(DEPTHM(KL)-DEPTHM(KL-1))
	      SIGBT=SIGBT+BFSQ(KL)*(DEPTHM(KL)-DEPTHM(KL-1))
350   CONTINUE
      IF(SIGBT.gt.zero) HSIG=SIGTP/SIGBT
      IF(HSIG .gt. XMOM**two)HSIG=XMOM**two
360   CONTINUE
      DSIG=zero
      IF(HSIG.gt.zero) DSIG=HSIG**pt5
!
! Find layer above which 85% of BVSQ lies and the volume containing the 85%, VTILDA
!
      DO 370 i=1,ns-1
	      IF(depth(i) .gt. (XMOM-DSIG)) GOTO 380
370   CONTINUE

      i = ns-1
380   VTILDA = vol1(ns)
      IF (i .gt. 1) VTILDA = (vol1(ns)-vol1(i-1))
!
! Calculate rate of working by inflow, EINFF is the rate of energy
! released by plunging i.e. the change in potential energy (see INFLOW)
!
      VOLMX=VTILDA-vol(ns)
      DBAR=thsnd+(den(1)+den(ns))/two
      EINFW=EINFF/(VOLMX*thsnd*DBAR)
!
! Include rate of working of the wind 
!
      EW=WP*area(ns)*arfac
      EWW=EW/(VTILDA*thsnd*DBAR)
      IF ((EW + EINFF) .le. zero) GOTO 500
!
! Calculate dissipation, velocity scale and wave number squared
!
      DISS = EWW+EINFW
!
! When PARTICLES activated, math error IF diss less than zero
!
	   IF(diss.lt.zero) THEN
	      WRITE(*,*)'Warning, diss less than zero. ERROR in ENER'
	      diss = zero
!	      pause
	   ENDIF
      IF(EWW.LE.EINFW) GOTO 500
      WNSQ = cwnsq1*area(ns)*arfac/(VTILDA*thsnd*H0)
      VEL = USTAR
500   CONTINUE

      RETURN
      END SUBROUTINE ENER
!********************************************************************	
      SUBROUTINE mixer(aeco,aews,aketil,ake)
!	subroutine mixer(taeco,taews,taketil,take)
!********************************************************************
!    Performs the surface mixing due to wind forcing
!---------------------------------------------------------------------
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  
	   REAL(8)	:: taeco,taews,taketil,take      

	   REAL(8), PARAMETER :: zero   = 0.0d0,one = 1.0d0,two = 2.0d0
	   REAL(8), PARAMETER :: three  = 3.0d0,six = 6.0d0,ten = 10.0d0
	   REAL(8), PARAMETER :: twelve = 12.0d0,twfour = 24.0d0,pt25= 0.25d0
	   REAL(8), PARAMETER :: pt3 = 0.3d0,pt587 = 0.587d0,pt6 = 0.6d0
	   REAL(8), PARAMETER :: g = 9.81d0,secshr = 3600.0d0
	   REAL(8), PARAMETER :: tenm10 = 1.0d-10,tenm5 = 1.0d-5
	   REAL(8), PARAMETER :: tifac = 1.587d0, tdfac = 8.33d-4
	   REAL(8), PARAMETER :: volfac = 1.0d+3,arfac = 1.0d+6
	   REAL(8), PARAMETER :: thsnd = 1000.0d0

	   REAL(8)	:: sech,x
	   REAL(8) :: aeco,aews,ake,aketil,ase
	   REAL(8) :: bot
	   REAL(8) :: ci,cmsml,cmm
	   REAL(8) ::    caca_part_I(7),caca_part_F(7)
	   REAL(8) :: deltcm,deltax,deltsq,delu,dh,dhypl,dz
	   REAL(8) :: el,ew
	   REAL(8) :: fn
	   REAL(8) :: gdash,gpeffc
!IF   REAL*8  ::gprime
	   REAL(8) :: hb,ht,htb,htilda
	   REAL(8) :: qcub,qsq
	   REAL(8) :: sbig,ssml,slope
	   REAL(8) :: tbig,tsml,td,tieff,top
	   REAL(8) :: uasave,uavsq,ueff,uisave,ufocub,ustar,ustsq
	   REAL(8) :: vb,vmbig,vmsml
	   REAL(8) :: wbot,wth,wtop
	   REAL(8) :: xmsml
	   REAL(8) :: c_drag

	   INTEGER*4 i,j,ix,zz
	   INTEGER*4 ij	! for mixer bug fix
	   INTEGER*4 i_p,j_p
	   INTEGER*4 ucheck

! New definitions
	   REAL(8)		:: spe_al,dh_al,h_al,db_al,gdash_al

! for mixer debug file
	   sech(x) = two/(exp(x)+exp(-x))
!
!  put in wind factor, common definition of u6x
!
!	ustar = sqrt(1.612e-6*u6x*u6x)
!
! New Linear definition of Cd (Smith and Banke formulae, from Henderson-Sellers
!
      c_drag = 1.3d-3 
!      c_drag = (0.61 + 0.063 * u6x)/ 1000.
!	ucheck = 1000.0 * u6x	
!	SELECT CASE (ucheck)
!		CASE (:300)
!			c_drag = .0015			
!		CASE (301:2200)
!			c_drag = (1.08 * u6x**(-0.15))/1000.0		
!		CASE (2201:5000)
!			c_drag = (0.771 + 0.0858 * u6x)/1000.0		
!		CASE (5001:8000)
!			c_drag= (0.867 + 0.0667 * u6x)/1000.0		
!		CASE (8001:25000)
!			c_drag= (1.2 + 0.025 * u6x)/1000.0	
!		CASE (25001:50000)
!			c_drag= (0.073 * u6x)/1000.	
!		CASE default
!			c_drag= 0.0037	
!	END SELECT
	   ustar = sqrt(c_drag*(1.24d0/1000.0d0)*u6x*u6x)


	   mstep = mstep+1

!  perform mixing due to surface heat transfers, calculate
!  p.e. released by bouyancy flux and surface wind stress

	   tbig  = zero
	   tsml  = zero
	   sbig  = zero
	   ssml  = zero
	   vmbig = zero
	   vmsml = zero
	   xmsml = zero
	   cmsml = zero

	   DO 10 i = 1,ns
	      k1 = ns+1-i

   !  add layers to the upper mixed layers, df is returned from addlay
	      CALL addlay(vmbig,vmsml,tbig,tsml,sbig,ssml,vm,tm,sm,df,k1)

   !     calculate 0th and 1st moments of density about the bottom
	      IF(k1.ne.1) THEN
	         dz = depth(k1)-depth(k1-1)
	         xmsml = xmsml+den(k1)*dz
	         cmsml = cmsml+den(k1)*dz*depthm(k1)
	         IF(df.lt.den(k1-1)) GOTO 20
	      ELSE
	            xmsml = xmsml+den(1)*depth(1)
	            cmsml = cmsml+den(1)*depth(1)*depthm(1)
	      ENDIF
 10   CONTINUE

! 29/1/92 mixer bug fix for times when the wind speed is very low

15	   DO ij = k1,ns
	      temp(ij) = tm
	      sal (ij) = sm
	      den (ij) = df
	   ENDDO

20    db = zero
	   IF(k1.ne.1) db = depth(k1-1)
		cmm=(depth(ns)+db)/two
		h = depth(ns)-db

!  check for bottom, update time index IF necessary
		j1 = k1-1
		IF (j1.eq.0)THEN
			thr = thr+nosecs/3600.0d0
		   GOTO 800
	   ENDIF
	   deltcm = (cmsml-cmm*xmsml)

!  kraus turner, ufocub measures energy released by cooling
	   ufocub = g*deltcm/((df+thsnd)*nosecs)
	   IF (ufocub.lt.zero) ufocub = zero
	   ustsq = ustar*ustar

!  calculate total energy available for stirring and add to
!  amount stored from last time step.
	   aews  = ck*eta*eta*eta*ustar*ustar*ustar*nosecs/two
	   aeco  = ck*ufocub*nosecs/two
	   ase   = aeco+aews
	   taeco = taeco+aeco
	   taews = taews+aews
	   qcub  = ase/(ck*nosecs)*two
	   IF(qcub .le. zero) qcub = tenm10
	   qsq = qcub**(two/three)
	   aem = aem+ase

!  loop for stirring. check for bottom. compute energy required to
!  mix next layer and compare with available energy.
100   h  = depth(ns) - depth(j1)
	   db = zero

	   IF(j1.gt.1) db = depth(j1-1)
	   dh    = depth(j1)-db
	   gdash = gprime(df,den(j1))
	   spe   = (gdash*h+ct*qsq)*dh/two
	   IF(aem.lt.spe) GOTO 200

!  entrain layer j1
	   call addlay(vmbig,vmsml,tbig,tsml,sbig,ssml,vm,tm,sm,df,j1)
	   call aver
	   IF(j1 .gt. 0) GOTO 100

!  here IF kraus-turner mixes to bottom
	   thr = thr+nosecs/3600.0d0
	   GOTO 800
200   CONTINUE

!  prt - shear production
!  cutoff shear production IF thermocline occurs within the hole
	   IF(depth(j1) .le. hle) GOTO 900

!  calculate kraus-turner depth
	   htilda = depmx-depth(j1)
	   IF(htilda.le.zero) htilda = zero

!  calculate parameters needed for momentum computation
!  ci is the wave speed along the thermocline for a two layer fluid
!  htsave is the mixed layer thickness from previous time step
!  dhypl is the mean hypolimnion density
!  df is the upper mixed layer density
!  ht,hb are the thicknesses of the uml and hypolimnion respectively
!  calculated as volume/mean area

	   wth = zero
	   DO i = 1,j1
	      wth = vol(i)*den(i) + wth
	   ENDDO
	   dhypl = wth/vol1(j1)
	   gpeff = gprime (df,dhypl)
	   vf = vol1(ns) - vol1(j1)
	   vb = vol1(j1)
	   ht = (vf/(area(j1)+area(ns)))*(volfac/arfac)*two
	   hb = (vb/area(j1))*(volfac/arfac)*two
	   IF (gpeff .le. 0.1d-6) gpeff = 0.1d-6
	   ci = sqrt((abs(gpeff)*ht*hb)/(ht+hb))

!     adjust momentum for fluid entrained by kraus turner deepening
!     IF ht > htsave

	   uasave = uav
	   IF (uav.le.zero .or. htsave .ge. ht) GOTO 300
	   uf = uf*htsave/ht
	   ui = uf
300   CONTINUE
	   uisave = ui

!  compute effective length at thermocline level

	   db = depth(j1)
	   DO 325 i = numout+1,1,-1
	      top = crl
	      wtop = wc
	      IF (i .ne. numout+1) top	 = olev(i)
	      IF (i .ne. numout+1) wtop = owid(i)
	      bot = zero
	      wbot = zero
	      DO 330 j = 1,numout
		      IF (olev(j) .gt. bot .and. olev(j) .lt. top) THEN
		         bot  = olev(j)
		         wbot = owid(j)
		      ENDIF
330      CONTINUE
	      IF (db .gt. bot .and. db .le.top) THEN
		      ew = (depth(j1)-bot)/(top-bot)*(wtop-wbot)+wbot
	      ENDIF
325   CONTINUE
	   el = area(j1)*(arfac)/ew

!     check momentum time counters (ti > 0 indicates a current shear event)
!     tieff is the effective forcing time
!     ti is one half the seiche period
!     timei is the time (hrs since start of sim.) at the start of shear forcing
!     timef is the time (hrs since start of sim.) at the END of shear forcing
!     td is the damping time

	   IF (ti .gt. zero) GOTO 400
	   ti = el/(two*ci*secshr)
	   tieff = ti
	   IF(ustar.le.zero)THEN
		   td = zero
	   ELSE
		   htb = ht+hb
	      td = two*(htb/tdfac)*(htb/ustsq)*(hb/ht)*(gpeff*ht*hb/htb)**pt25/  &
               (sqrt(two*el)*secshr)*sqrt(visc)
	      tieff = tifac*ti
	      IF(td/ti .lt. ten)tieff=(one+pt587*(one-sech(td/ti-one)))*ti
	   ENDIF
	   timei  = thr
	   timefi = timei + tieff
400   CONTINUE

!     calculate momentum forcing parameters for current time step
!     fn is the acceleration (m/s**2) of uml by wind stress

	   fn    = ustsq/ht
	   fsum  = fsum + fn
	   slope = (fn-fo) + oldsl
	   IF (fn .eq. zero) THEN
		   slope = zero
	   ELSE
	      IF(fn.le.zero .or. abs(slope/fn).le.tenm5)slope = zero
	   ENDIF
	   IF (slope .lt. zero) slope = zero

!  check for momentum cutoff within current time step. calculate time
!  step for forcing, ftime, and reset parameters for next time step

	   thr = thr + nosecs/3600.0d0
	   IF (thr .ge. timefi) THEN
!  here IF cutoff within current time step

	      ftime = timefi - thr + nosecs/3600.0d0
	      oldsl = fn - (fsum / float(mstep))
	      IF (oldsl .lt. zero) oldsl = zero
	   ELSE
	      ftime = nosecs/3600.0d0
	      oldsl = slope
	   ENDIF
	   fo = fn

!  compute momentum increment

	   IF (ui .lt. 1e-7)    ui    = zero
	   IF (slope .lt. 1e-7) slope = zero
	   uf     = ui + slope * ftime * secshr
	   uavsq  = (uf * uf + uf * ui + ui * ui) / three
	   uav    = sqrt(uavsq)
	   ui     = uf
	   delu   = uav - uasave
	   deltsq = pt6 * uav * delu / gpeff
	   deltax = pt3 * uav * uav  / gpeff

	   IF (deltax .lt. 1.0e-10) deltax = 0.0d0
	   IF (deltsq .lt. 1.0e-10) deltsq = 0.0d0

	   aketil = cs * (uav * uav * (htilda + deltsq / six) +                    &
              uav * deltax * delu / three) /two + gdash * deltax * (deltax *  &
              htilda /(twfour * (depth(ns) - depth(j1))) - deltsq / twelve)

	   taketil = taketil + aketil
	   aem     = aem     + aketil
	   gpeffc  = gpeff * ht

!  deepening loop for shear production begins here

	   delu = zero
500   CONTINUE
!  save current values of ht and uav in CASE of mixing
	   htsave = ht
	   ueff   = uav
!  compute energy available for mixing next layer

	   db = zero
	   IF (j1 .gt. 1) db = depth(j1 - 1)
	   dh = depth(j1) - db
	   h  = depth(ns) - depth(j1)
	   ake = cs  * (ueff  * ueff * (dh + deltsq / six) + ueff  * deltax * delu / three) / two +    &
              gdash * deltax * (deltax * dh /  (twfour * h) - deltsq / twelve)

	   aem = aem + ake
!klg for mixer debug file
      take = take + ake
!  compute energy required to entrain next layer
	   spe = (gdash * h + ct * qsq) * dh / two

!  compare energy available with energy required
	   IF (aem .lt. spe) GOTO 700

!     entrain layer j1
	   call addlay(vmbig,vmsml,tbig,tsml,sbig,ssml,vm,tm,sm,df,j1)
	   call aver

	   IF (j1 .eq. 0) GOTO 800

!  adjust uf, uav for entrained mass
	   wth = zero
	   DO 600 i = 1,j1
	      wth = vol(i)*den(i)+wth
600   CONTINUE
	   dhypl = wth/vol1(j1)
	   gdash = gprime(df,den(j1))

	   vb     = vol1(j1)
	   vf     = vf + vol(j1+1)
	   ht     = (vf/(area(j1)+area(ns)))*(volfac/arfac)*two
	   uf     = uf*htsave/ht
	   uav    = sqrt((uisave*uisave + uisave*uf + uf*uf)/three)
	   delu   = uav-ueff
	   deltax = pt3*uav*uav/(gdash)
	   deltsq = pt6*uav*delu/(gdash)
	   ui     = uf
	   GOTO 500
700   CONTINUE

!     here if insufficient energy to entrain next layer
!     check momentum time counters for cutoff
	   db    = depth(j1)
	   gpeff = gpeffc/ht

!     average water quality

	   DO zz = 1,28
	      wqualm(zz) = 0.0d0
	      DO ix = j1+1,ns
		      wqualm(zz) =(wqual(ix,zz)*vol(ix))+wqualm(zz)
	      ENDDO
	      wqualm(zz) = wqualm(zz)/(vol1(ns)-vol1(j1))
	   ENDDO

	   DO zz = 1,7
	      cfm(zz) = 0.0d0
	      DO ix = j1+1,ns
		      cfm(zz) =(cf(ix,zz)*vol(ix))+cfm(zz)
	      ENDDO
	      cfm(zz) = cfm(zz)/(vol1(ns)-vol1(j1))
	   ENDDO
!
! Kelvin-H. billowing
!
	   CALL kh

	   depmx = depth(j1)

	   IF(thr .lt. timefi)  GOTO 1000
	   GOTO 900
800   CONTINUE

!  here if deepened to bottom
	   depmx = zero   ! New addition quim (Feb 2004)
	   oldsl = zero
	   aem   = zero
	   fo    = zero

	   DO zz = 1,28
	      wqualm(zz)  = 0.0d0
	      DO ix = 1,ns
			   wqualm(zz) = (wqual(ix,zz)*vol(ix))+wqualm(zz)
	      ENDDO
	      wqualm(zz)  = wqualm(zz)/vol1(ns)
	   ENDDO

	   DO zz = 1,7
		   cfm(zz) = 0.0d0
			   DO ix = 1,ns
				   cfm(zz) =(cf(ix,zz)*vol(ix))+cfm(zz)
			   ENDDO
		   cfm(zz) = cfm(zz)/vol1(ns)
	   ENDDO
900   CONTINUE

!  here if momentum cutoff

	   mstep  = 0
	   fsum   = zero
	   ti     = zero
	   uf     = zero
	   ui     = zero
	   uav    = zero
	   timei  = thr
	   timefi = thr

!  mark 2 ends here. at this stage layers j1+1,j1+2,---k1-1 are also mixed
!  renumber mixed layers

1000  depth(j1+1) = depth(ns)
	   temp(j1+1)  = tm
	   sal(j1+1)   = sm

! include water quality and particles in insertion

	   DO zz = 1,28
	      wqual(j1+1,zz) = wqualm(zz)
	   ENDDO

	   DO zz = 1,7
	      cf(j1+1,zz) = cfm(zz)
	   ENDDO

	   vol1(j1+1) = vol1(ns)
	   vb         = zero

	   IF(j1.ne.0) vb = vol1(j1)
		vol(j1+1)  = vol1(ns)-vb
		den(j1+1)  = df
		area(j1+1) = area(ns)
		ns         = j1+1

2000  CONTINUE

	   RETURN
      END SUBROUTINE Mixer

!***********************************************************************************
      SUBROUTINE ADDLAY(VMbig,VMsml,Tbig,Tsml,Sbig,Ssml,VMLOC,TMLOC, SMLOC,DFLOC,L)
!*************************************************************************************
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  

      REAL(8), PARAMETER :: thsnd=1000.0d0

      REAL(8) ::	DFLOC
!IF      REAL(8) ::  densty
      REAL(8) ::  Sbig
      REAL(8) ::  Ssml
      REAL(8) ::  SMLOC
      REAL(8) ::  Tbig
      REAL(8) ::  Tsml
      REAL(8) ::  TMLOC
      REAL(8) ::  VMLOC
      REAL(8) ::  VMbig
      REAL(8) ::  VMsml
      REAL(8) ::  WTbig
      REAL(8) ::  WTsml

      INTEGER*4 L

      WTbig = thsnd * vol(L)
      WTsml = den(L)* vol(L)
      VMbig = VMbig + WTbig
      VMsml = VMsml + WTsml
      Tbig  = Tbig  + temp(L)*WTbig
      Tsml  = Tsml  + temp(L)*WTsml
      Sbig  = Sbig  + sal(L)*WTbig
      Ssml  = Ssml  + sal(L)*WTsml

      VMLOC = VMbig + VMsml
      TMLOC = Tsml/VMLOC + Tbig/VMLOC
      SMLOC = Ssml/VMLOC + Sbig/VMLOC
      DFLOC = densty(TMLOC,SMLOC)

      RETURN
      END SUBROUTINE ADDLAY
!*********************************************************************
      SUBROUTINE AVER
!*************************************************************************
! This routine assigns the mean mixed layer properties to all of the
! layers in the epilimnion and increments the layer pointers j1,k1
!-----------------------------------------------------------------------------
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  
      INTEGER*4 i

      AEM = AEM-SPE

      DO i = j1,ns
	      temp(i) = TM
	      sal(i)	 = SM
	      den(i)  = DF
	   ENDDO
 
      j1 = j1-1
      k1 = k1-1

      RETURN
      END SUBROUTINE AVER
! ******************************************************************
      SUBROUTINE KH
!*******************************************************************
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  

      REAL(8), PARAMETER :: zero=0.0d0,two=2.0d0,three=3.0d0
	   REAL(8), PARAMETER :: six=6.0d0,ten=10.0d0,pt02=0.02d0
	   REAL(8), PARAMETER :: g=9.81d0,secshr=3600.0d0

      REAL(8) :: COSH
!IF      REAL(8) :: COMBIN
!IF      REAL(8) :: COMBINV    
      REAL(8) :: DELTA
!IF     REAL(8) :: densty
      REAL(8) :: DEPTHB
      REAL(8) :: DNL
      REAL(8) :: DSAVE
      REAL(8) :: D3
      REAL(8) :: E
      REAL(8) :: EE
      REAL(8) :: EBOT
      REAL(8) :: ETOP
      REAL(8) :: HMIN
      REAL(8) :: T
      REAL(8) :: TBIL
      REAL(8) :: TBILF
      REAL(8) :: THB
      REAL(8) :: THT
      REAL(8) :: x

      INTEGER*4 i,zz
      INTEGER*4 icode
      INTEGER*4 IR
      INTEGER*4 KL1
      INTEGER*4 KM1
      INTEGER*4 L
      INTEGER*4 LBI
      INTEGER*4 LLB
      INTEGER*4 lnu
      INTEGER*4 N
      INTEGER*4 NL

      LOGICAL*4 TOP
      LOGICAL*4 UP	
      COSH(x)=(EXP(x)+EXP(-x))/two

!  SET TOLERANCES

      E	=	pt02
      EE	=	six*E

!  COMPUTE DELTA

      TBIL	= UAV/(GPEFF*secshr)
      TBILF	= TBIL/FTIME
      DELTA	= zero
      IF (ns .le. 1 .OR. TBILF .gt. ten) GOTO 299
      DELTA = (AKH*UAV*UAV)/(GPEFF*two*COSH(TBILF))
      HMIN	= MIN(H,DB)                                                                 ! Problems with DMIN1 definition
      IF (DELTA .gt. HMIN) DELTA=HMIN
      ETOP	= DB+DELTA
      EBOT	= DB - DELTA
      TOP		= .FALSE.

      IF(ETOP.gt.depth(ns-1)) ETOP = depth(ns-1)
      IF((ETOP-DB).lt.EE/two) GOTO 299
      IF(EBOT.lt.depth(1))	EBOT = depth(1)
      IF((DB-EBOT).lt.EE/two) GOTO 299
      
	   DSAVE = depth(ns)

!  FIND LAYER INTERSECTING EBOT

      KM1 = k1-1
      DO 20 L = 1,KM1
	      IF(depth(L) .gt. EBOT) GOTO 30
20    CONTINUE
30    CONTINUE
      DEPTHB = zero
      IF (L.gt.1) DEPTHB = depth(L-1)

!  CHECK TO SEE IF EBOT COINCIDES WITH EXISTING depth VALUE

      T = ABS(EBOT-DEPTHB)
      IF (T.gt.E) GOTO 40
      EBOT = DEPTHB
      LLB = L-1
      GOTO 60
40    CONTINUE
      T = ABS(depth(L) - EBOT)
      IF (T .gt. E) GOTO 50
      EBOT = depth(L)
      LLB = L
      IF (L .eq. k1-1) GOTO 299
      GOTO 60

!  HERE IF NEW LAYER MUST BE ADDED BELOW MIXED REGION

50    CONTINUE
      LLB = L
      k1	= k1+1
      KL1 = k1-L-1

      DO 55 i = 1,KL1
	      j1 = k1-i
	      depth(j1)	= depth(j1-1)
	      den(j1)	= den(j1-1)
	      temp(j1)	= temp(j1-1)
	      sal(j1)	= sal(j1-1)
	      DO 51 zz=1,28
	         wqual(j1,zz) = wqual(j1-1,zz)
51       CONTINUE
	      DO 52 zz=1,7
	         cf(j1,zz)  = cf(j1-1,zz)
52       CONTINUE
55    CONTINUE
      depth(LLB) = EBOT
60    CONTINUE

!  HERE AFTER POSITION (EBOT) OF BOTTOM OF SHEAR ZONE HAS BEEN DETER-
!  MINED AND EXTRA LAYER ADDED IF NECESSARY

      T = ABS(ETOP-EBOT)
      IF (T.lt.EE) GOTO 299

!  CHECK NUMBER OF LAYERS IN BOTTOM HALF OF SHEAR LAYER - THERE MUST
!  BE AT LEAST THREE

      NL = k1-LLB-1
      IF (NL.GE.3) GOTO 100
      IF (NL.eq.2) GOTO 80

!  HERE IF NL=1

      NL = 3
      D3 = (depth(k1-1) - EBOT)/three
      depth(LLB+3) = depth(k1-1)
      depth(LLB+2) = EBOT + D3 + D3
      depth(LLB+1) = EBOT + D3
      DO 70 i = 2,3
	      LBI = LLB + i
	      den(LBI)	= den(k1-1)
	      temp(LBI)	= temp(k1-1)
	      sal(LBI)	= sal(k1-1)
	      DO 62 zz=1,28
	         wqual(LBI,zz)=wqual(k1-1,zz)
62       CONTINUE
	      DO 63 zz=1,7
	         cf(LBI,zz)=cf(k1-1,zz)
63       CONTINUE
70    CONTINUE
      k1 = LLB+4
      GOTO 100
80    CONTINUE

!  HERE IF NL=two

      NL=3

!  FIND THICKER LAYER

      THT = depth(k1-1) - depth(k1-2)
      THB = depth(k1-2) - EBOT
      IF (THT .gt. THB)  GOTO 90

!  HERE IF THT .le. THB - DIVIDE LAYER k1-2

      depth(LLB+1) = EBOT + THB/two
      depth(LLB+2) = EBOT + THB
      depth(LLB+3) = EBOT + THT + THB
      k1 = k1+1
      den(k1-1)	= den(k1-2)
      sal(k1-1)	= sal(k1-2)
      temp(k1-1)	= temp(k1-2)
      DO 81 zz=1,28
	      wqual(k1-1,zz)=wqual(k1-2,zz)
81    CONTINUE      
	   DO 82 zz=1,7
	      cf(k1-1,zz)=cf(k1-2,zz)
82    CONTINUE
      den(k1-2) = den(k1-3)
      temp(k1-2) = temp(k1-3)
      sal(k1-2) = sal(k1-3)
      DO 85 zz=1,28
	      wqual(k1-2,zz) = wqual(k1-3,zz)
85    CONTINUE
      DO 86 zz=1,7
	      cf(k1-2,zz) = cf(k1-3,zz)
86    CONTINUE

      GOTO 100     
90    CONTINUE

!  HERE IF THT > THB - SPLIT LAYER k1-1
	
      k1=k1+1
      depth(k1-1) = depth(k1-2)
      sal(k1-1)	= sal(k1-2)
      temp(k1-1)	= temp(k1-2)
      DO 93 zz=1,28
	      wqual(k1-1,zz) = wqual(k1-2,zz)
93    CONTINUE     
      DO 94 zz=1,7
	      cf(k1-1,zz) = cf(k1-2,zz)
94    CONTINUE
      den(k1-1)	= den(k1-2)
      depth(k1-2) = depth(k1-3) + THT/two
100   CONTINUE

!  HERE AFTER BOTTOM HALF OF SHEAR ZONE HAS AT LEAST THREE LAYERS
!  DIVIDE TOP HALF OF SHEAR ZONE INTO NL LAYERS

      DNL = (ETOP - depth(k1-1))/float(NL)
      DO 110 i = 1,NL
	      j1 = i+k1-1
	      depth(j1)	= depth(k1-1) + float(i)*DNL
	      den(j1)	= DF
	      temp(j1)	= TM
	      sal(j1)	= SM
	      DO 102 zz=1,28
	         wqual(j1,zz) = WQUALM(zz)
102      CONTINUE
	      DO 103 zz=1,7
	         cf(j1,zz) = CFM(zz)
103      CONTINUE
110   CONTINUE

!  THE NUMBER OF THE LAYER JUST BELOW THE MIXED REGION IS NOW j1=k1+NL-1
!  UNLESS TOP = .T.  IF TOP=.T.,THEN LAYER j1 IS THE MIXED REGION

      IF (TOP) GOTO 120
      ns=k1+NL
      temp(ns) = TM
      den(ns) = DF
      sal(ns) = SM
      DO 112 zz=1,28
	      wqual(ns,zz) = WQUALM(zz)
112   CONTINUE
	   DO 113 zz=1,7
	      cf(ns,zz)=CFM(zz)
113   CONTINUE
      GOTO 130
120   CONTINUE		     
      j1 = j1-1
      ns = k1+NL-1
130   CONTINUE
      depth(ns) = DSAVE
      H = depth(ns) - depth(ns-1)

!     CALCULATE VOLUMES

      icode = 1
      lnu   = 1

      CALL RESINT(icode,lnu)

!  RELAX DENSITY STRUCTURE WITHIN SHEAR ZONE

      UP = .TRUE.
      IR = NL-1
140   CONTINUE

!  MIX MIDDLE TWO LAYERS k1, k1-1

      temp(k1) = COMBIN(temp(k1),vol(k1),den(k1),temp(k1-1),vol(k1-1),den(k1-1))
      sal(k1) = COMBIN(sal(k1),vol(k1),den(k1),sal(k1-1),vol(k1-1),den(k1-1))
      den(k1) = densty(temp(k1),sal(k1))
      DO 141 zz=1,28
	      wqual(k1,zz) = COMBINV(wqual(k1,zz),vol(k1),wqual(k1-1,zz),vol(k1-1))
141   CONTINUE
      DO 142 zz=1,7
	      cf(k1,zz) = COMBINV(cf(k1,zz),vol(k1),cf(k1-1,zz),vol(k1-1))
142   CONTINUE 

      temp(k1-1) = temp(k1)
      sal(k1-1) = sal(k1)
      den(k1-1) = den(k1)
      DO 146 zz=1,28
	      wqual(k1-1,zz)=wqual(k1,zz)
146   CONTINUE
      DO 147 zz=1,7
	      cf(k1-1,zz)=cf(k1,zz)
147   CONTINUE

      IF (IR .eq. 1) GOTO 180
150   CONTINUE

!  DO LOOP MIXES UP (UP = .T.) OR DOWN (UP = .F.)

      DO 160 i = 1,IR
	      IF (UP) N=k1+i-1
	      IF (.NOT. UP) N=k1-i-1
	      temp(N) = COMBIN(temp(N),vol(N),den(N),temp(N+1),vol(N+1),den(N+1))
	      sal(N) = COMBIN(sal(N),vol(N),den(N),sal(N+1),vol(N+1),den(N+1))
	      DO 152 zz=1,28
	         wqual(N,zz) = COMBINV (wqual(N,zz),vol(N),wqual(N+1,zz),vol(N+1))
152      CONTINUE
	      DO 153 zz=1,7
	         cf(N,zz)=COMBINV(cf(N,zz),vol(N),cf(N+1,zz),vol(N+1))
153      CONTINUE
	 
	      DO 154 zz=1,28
	         wqual(N+1,zz)=wqual(N,zz)
154      CONTINUE

	      DO 155 zz=1,7
	         cf(N+1,zz)=cf(N,zz)
155      CONTINUE
160   CONTINUE

      IF (.NOT. UP) GOTO 170
	   UP = .FALSE.
	   GOTO 150
170   CONTINUE
      UP = .TRUE.
      IR = IR-1
      GOTO 140
180   CONTINUE

!  HERE AFTER RELAXATION COMPLETE. RESET MIXED REGION VARIABLES.
			
      DF = den(ns)
      TM = temp(ns)
      SM = sal(ns)
      DO 182 zz=1,28
	      WQUALM(zz)=wqual(ns,zz)
182   CONTINUE
	   DO 183 zz=1,7
	      CFM(zz)=cf(ns,zz)
183   CONTINUE

      VF = vol(ns)
      VM = VF*DF
      GOTO 300
299   CONTINUE

!     HERE IF SHEAR LAYER IS TOO THIN

      NL=0
300   CONTINUE

      RETURN
      END SUBROUTINE KH

!********************************************************************************************
      SUBROUTINE INFLOW(xinf,N_Ins,P_Ins,N_Ent,P_Ent,N_Inf,P_Inf,                &
         Part_Ent,Part_Inf,Part_Ins,Balance_N_Stack,Balance_P_Stack,             &
         Balance_Part_Stack,Part_Stack_F,Part_Stack_I,t1,t2,t3,ts_flows)
!********************************************************************************************
!  Inserts the inflow at the level of neutral buoyancy, after calculating the entrainment
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  

      REAL(8), PARAMETER :: zero=0.0d0,one=1.0d0,two=2.0d0
	   REAL(8), PARAMETER :: three=3.0d0,five=5.0d0,six=6.0d0
	   REAL(8), PARAMETER :: ten=10.0d0,pt1=0.1d0,pt15=0.15d0
	   REAL(8), PARAMETER :: pt2=0.2d0,pt21=0.21d0,pt4=0.4d0,pt9=0.9d0
      REAL(8), PARAMETER :: onept2=1.2d0,onept5=1.5d0,onept6=1.6d0
	   REAL(8), PARAMETER :: cwnsq2=39.48d0,g=9.81d0
	   REAL(8), PARAMETER :: pideg=180.0d0,secday=86400.0d0
	   REAL(8), PARAMETER :: thsnd=1000.0d0,tsecda=86.4d0
       REAL(8) tsectimestep
       
	   REAL(8), PARAMETER :: pi=3.1415926535897932384626433832795d0

      REAL*8 BAN	! Bed half angle
      REAL*8 BSL	! Bed Slope
      REAL*8 CD
      REAL*8 DELQ,DELT, norm_DELT
      REAL*8 DEP
      REAL*8 DEPT0,DEPTH0
      REAL*8 DI
!gbs 17Jan06      REAL*8 DSTART
	   REAL*8 DSTART(maxinf) !, Bedslope(maxinf)    !gbs 17Jan06
      REAL*8 DX
      REAL*8 E
      REAL*8 EINF
      REAL*8 GDASH
!IF      REAL*8 GPRIME,densty,COMBIN,COMBINV
      REAL*8 GRAV
      REAL*8 HF0
      REAL*8 HH
      REAL*8 Q
      REAL*8 RI
      REAL*8 S
!IF      REAL*8 SQR
      REAL*8 T
      REAL*8 TIMEIN, TIM
      REAL*8 VEL1
      REAL*8 WIDTH      
      REAL*8 XINF(MAXINF)
      REAL*8 CCFF1,CCFF2,CCFF3,CCFF4,CCFF5,CCFF6,CCFF7
      REAL*8 CF1(maxinf,maxpar)
      REAL*8 CF2(maxinf,maxpar)
      REAL*8 CF3(maxinf,maxpar)
      REAL*8 CF4(maxinf,maxpar)
      REAL*8 CF5(maxinf,maxpar)
      REAL*8 CF6(maxinf,maxpar)
      REAL*8 CF7(maxinf,maxpar)
      REAL*8 WQ1,WQ2,WQ3,WQ4,WQ5,WQ6,WQ7,WQ8,WQ9,WQ10,WQ11
      REAL*8 WQ12,WQ13,WQ14,WQ15,WQ16,WQ17,WQ18,WQ19,WQ20
      REAL*8 WQ21,WQ22,WQ23,WQ24,WQ25,WQ26,WQ27,WQ28
      REAL*8 WQUAL1(maxinf,maxpar)
      REAL*8 WQUAL2(maxinf,maxpar)
      REAL*8 WQUAL3(maxinf,maxpar)
      REAL*8 WQUAL4(maxinf,maxpar)
      REAL*8 WQUAL5(maxinf,maxpar)
      REAL*8 WQUAL6(maxinf,maxpar)
      REAL*8 WQUAL7(maxinf,maxpar)
      REAL*8 WQUAL8(maxinf,maxpar)
      REAL*8 WQUAL9(maxinf,maxpar)
      REAL*8 WQUAL10(maxinf,maxpar)
      REAL*8 WQUAL11(maxinf,maxpar)
      REAL*8 WQUAL12(maxinf,maxpar)
      REAL*8 WQUAL13(maxinf,maxpar)
      REAL*8 WQUAL14(maxinf,maxpar)
      REAL*8 WQUAL15(maxinf,maxpar)
      REAL*8 WQUAL16(maxinf,maxpar)
      REAL*8 WQUAL17(maxinf,maxpar)
      REAL*8 WQUAL18(maxinf,maxpar)
      REAL*8 WQUAL19(maxinf,maxpar)
      REAL*8 WQUAL20(maxinf,maxpar)
      REAL*8 WQUAL21(maxinf,maxpar)
      REAL*8 WQUAL22(maxinf,maxpar)
      REAL*8 WQUAL23(maxinf,maxpar)
      REAL*8 WQUAL24(maxinf,maxpar)
      REAL*8 WQUAL25(maxinf,maxpar)
      REAL*8 WQUAL26(maxinf,maxpar)
      REAL*8 WQUAL27(maxinf,maxpar)
      REAL*8 WQUAL28(maxinf,maxpar)
!
! Nutrient Budget
!
	   REAL*8 N_Ins,P_Ins,N_Ent,P_Ent,N_Inf,P_Inf
	   REAL*8 Part_Ent(7),Part_Inf(7),Part_Ins(7)
	   REAL*8 N_Stack_F,N_Stack_I
	   REAL*8 P_Stack_F,P_Stack_I
	   REAL*8 Part_Stack_F(7),Part_Stack_I(7)
	   REAL*8 Balance_N_Stack,Balance_P_Stack, Balance_Part_Stack(7)
	   REAL*8 Caca_Part_F(7),caca_Part_I(7)

      INTEGER*4 j,k,L,M,LL,zz,I, I_P,J_P
      INTEGER*4 IAR,IRIV,ntims,icode
      INTEGER*4 JK
      INTEGER*4 kk
      INTEGER*4 LAYER
      INTEGER*4 ln
      INTEGER*4 IRIVER, P
!IF	   INTEGER*4 jday
!********River temperature variation*********************************
      INTEGER*4 totaltimestep, lateflow,delayedtimestep,t3
      REAL*8 increment,decrease,t1,t2
	   REAL*8 mintemp(100),maxtemp(100),Avetemp(100),temperature(100,100)
!********************************************************************
      REAL*8 ts_flows     
!----------------------------------------------------------------------      
!  FLERR is a flag, which IF true indicates that the inflow doesn't 
!  fit into the assumed geometry.  This is corrected by routing the 
!  oldest inflow parcel directly to the level of neutral buoyancy.
!  Entrainment calculations are included.

      LOGICAL*4 FLERR      
	   FLERR = .FALSE.
!      open(unit=122,file='testing.txt',status='unknown',access='append') 
! Inicialize Nutrient parameters
!
      N_Ins        = 0.0D0
	   P_Ins        = 0.0D0
      N_Ent        = 0.0D0
	   P_Ent        = 0.0D0
	   N_Inf        = 0.0D0
	   P_Inf        = 0.0D0
	   N_Stack_I    = 0.0D0	
	   P_Stack_I    = 0.0D0
	   N_Stack_F    = 0.0D0	
	   P_Stack_F    = 0.0D0
      ts_flows     = 0.0D0
       tsectimestep=REAL(RIVFLO_TIMESTEP)/1000.d0   !REPLACED THE DAILY TIME SCALE NUMBER 86.4 I.E., 86400/1000
	   DO i_p=1,7
	      Part_Stack_F(i_p)       = 0.0D0
	      Part_Stack_I(i_p)       = 0.0D0
	      Balance_Part_Stack(i_p) = 0.0D0
	      Part_Ins(i_p)           = 0.0D0
	      Part_Inf(i_p)           = 0.0D0
	      Part_Ent(i_p)           = 0.0D0
	   ENDDO

	   Balance_N_Stack  = 0.0D0
	   Balance_P_Stack  = 0.0D0
!
! Nutrient Budget
!
      DO i = 1,numinf
	      DO k = 1,icnt(i)
	         N_Stack_I = N_Stack_I + (wqdown(i,21,k)+ wqdown(i,22,k) + wqdown(i,23,k)+           &
     		                           wqdown(i,24,k) + wqdown(i,25,k) + wqdown(i,26,k))*qdown(i,k)
	         P_Stack_I = P_Stack_I + (wqdown(i,15,k) + wqdown(i,16,k) + wqdown(i,17,k) +           &
     		                            wqdown(i,18,k)+ wqdown(i,19,k) + wqdown(i,20,k))*qdown(i,k)
             DO i_p = 1,7
                Part_Stack_I(i_p) = Part_Stack_I(i_p)+ cfdown(i,i_p,k)*qdown(i,k)*1000.0D0
	         ENDDO  !  i_p  
         ENDDO  !  k
	   ENDDO  !  i

      caca_part_I = 0.0
	   caca_part_F = 0.0
      
	   DO i_p = 1,7
	      DO j_p = 1,ns
            caca_part_I(i_p) = caca_part_I(i_p) + cf(j_p,i_p)*vol(j_p)  
	      ENDDO
	   ENDDO

!  calculate WNSQ and VEL for use by sub DIFUSE; they're used to calculate the dispersion coefficient.
!  Count down from smallest river, overwriting each time, in CASE river 1 has zero inflow.
      WNSQ = zero
      VEL = zero
      DO 10 IRIVER=NUMINF,1,-1 !2015/9	      ts_flows=ts_flows+FLOINF(IRIVER)*XINF(IRIVER)*(REAL(nosecs))/(REAL(86400.0D0))
 
          ts_flows=ts_flows+FLOINF(IRIVER)*XINF(IRIVER)*(REAL(nosecs))/(REAL(RIVFLO_TIMESTEP))
	      IF (FLOINF(IRIVER)*XINF(IRIVER) .eq. zero) GOTO 10
	      BAN = ALPHA(IRIVER, numseg(iriver))*pi/pideg
			BSL=Bedslope(IRIVER)
!phi	 BSL = PHI(IRIVER)*pi/pideg
	      DI = densty(TEMINF(IRIVER),SALINF(IRIVER))
	      GRAV=GPRIME(den(ns),DI)	
!wef 11dec03	 DEPT0=DEPTH(NS)-DEPTH(NS-1)
	      IF (ns .gt. 1) THEN
		      DEPT0=DEPTH(NS)-DEPTH(NS-1) 
	      ELSE
		      depth0 = DEPTH(NS)
	      ENDIF
!2015/9	      IF(GRAV.gt.zero)DEPT0=((FLOINF(IRIVER)*XINF(IRIVER)/tsecda)**pt4)/(GRAV**pt2)
          IF(GRAV.gt.zero)DEPT0=((FLOINF(IRIVER)*XINF(IRIVER)/tsectimestep)**pt4)/(GRAV**pt2)          
	      IF(DEPT0.gt.depth(ns)) DEPT0=depth(ns)
	      WNSQ=cwnsq2/SQR(dept0)
!2015/9	      VEL=pt1*FLOINF(IRIVER)*XINF(IRIVER)/(tsecda*SQR(dept0)*SIN(BAN)/COS(BAN))
          	  VEL=pt1*FLOINF(IRIVER)*XINF(IRIVER)/(tsectimestep*SQR(dept0)*SIN(BAN)/COS(BAN))
10    CONTINUE

     
      DO 20 IRIVER = 1,NUMINF
	      IF (FLOINF(IRIVER)*XINF(IRIVER).LE.zero) GOTO 20
	      ICNT(IRIVER)=ICNT(IRIVER)+1
	      IAR = ICNT(IRIVER)
	      IF (IAR .gt. MAXPAR) THEN
	         WRITE(*,12)iriver
12          FORMAT(' Downflow stack limit too small '/            &
                   ' will give incorrect results.  In river ',i3)
	         STOP
	      ENDIF
!gbs	   QDOWN(IRIVER,ICNT(IRIVER)) = FLOINF(IRIVER)*XINF(IRIVER)
         PARCNT = PARCNT + 1
!2015/9         QDOWN(IRIVER,ICNT(IRIVER)) = FLOINF(IRIVER)*XINF(IRIVER)/(REAL(86400.0D0)/REAL(nosecs))
         QDOWN(IRIVER,ICNT(IRIVER)) = FLOINF(IRIVER)*XINF(IRIVER)/(REAL(RIVFLO_TIMESTEP)/REAL(nosecs))		       
!gbs*******************Beginning of sub-daily river temperature variation******************
         GOTO 598
         Avetemp(IRIVER)=teminf(IRIVER)
	      maxtemp(IRIVER)=teminf(IRIVER)+3.0		!       +0.3*teminf(IRIVER)
         mintemp(IRIVER)=teminf(IRIVER)-3.0     !	    -0.3*teminf(IRIVER)
	      IF(mintemp(IRIVER).le.0.0d0) mintemp(IRIVER)=0.0d0
!gbs       GOTO 598
         IF(nosecs.ge.10800) GOTO 598
         totaltimestep=86400/nosecs       ! one day=86400 seconds
	      lateflow=10800				     ! 3 hours=10800 seconds
	      delayedtimestep=10800/nosecs	
      
	      IF(iclock.eq.0)THEN	
            temperature(IRIVER,t3)=mintemp(IRIVER)
	         t2=t1	  
	         GOTO 598	
	      ENDIF
	      IF(IRIVER.eq.30) THEN    !Upper Truckee Bigger river for delayed flow
	         increment=pi/(0.5*totaltimestep+delayedtimestep)
	         decrease=pi/(0.5*totaltimestep-delayedtimestep)	
	         IF (t3.le.(1+0.5*totaltimestep+delayedtimestep)) THEN	     
	            t2=t1+increment	     
	         elseif(t3.gt.(1+0.5*totaltimestep+delayedtimestep).and.t3.le.(1+totaltimestep)) THEN	    
	            t2=t1-decrease	     
	         ENDIF  
	      ELSE                    !Smaller river
            increment=pi/(0.5*totaltimestep)                   
	         IF (t3.le.(1+0.5*totaltimestep)) THEN	    
	            t2=t1+increment	     
	         elseif(t3.gt.(1+0.5*totaltimestep).and.t3.le.(1+totaltimestep))THEN	     
	            t2=t1-increment	     
	         ENDIF        
         ENDIF
	      temperature(IRIVER,t3)=temperature(IRIVER,t3-1)+(maxtemp(IRIVER)-Avetemp(IRIVER))*(cos(t2)-cos(t1))	  
!        teminf(IRIVER)=temperature(IRIVER,t3)
!        print*,Iriver, t3,teminf(IRIVER)     
 598     CONTINUE   
!gbs*******************End of sub-daily river temperature variation******************
!GBS         TDOWN(IRIVER,ICNT(IRIVER)) = temperature(IRIVER,t3)
   	      TDOWN(IRIVER,ICNT(IRIVER)) = TEMINF(IRIVER) 
	      SDOWN(IRIVER,ICNT(IRIVER)) = SALINF(IRIVER)
	      DDOWN(IRIVER,ICNT(IRIVER)) = thsnd 
	      PDOWN(IRIVER,ICNT(IRIVER)) = PARCNT                            !initialize count of this parcel      
	      DO 14 zz=1,28
	         WQDOWN(IRIVER,zz,ICNT(IRIVER)) = WQINF(IRIVER,zz)	             
14      CONTINUE
!			print*,WQINF(1,21),WQINF(1,22),WQINF(1,16)				
	      DO 15 zz=1,7
	         CFDOWN(IRIVER,zz,ICNT(IRIVER)) = CFINF(IRIVER,zz)
15       CONTINUE
         TIMDOWN(IRIVER,ICNT(IRIVER)) = float (simday)+float(iclock)/86400.0d0 !Bill 2012/06
         
!
! Nutrient Budget
!
	      N_Inf = N_Inf +(wqinf(iriver,21) + wqinf(iriver,22) +	wqinf(iriver,23)  +         &
     	                   wqinf(iriver,24) + wqinf(iriver,25) + wqinf(iriver,26))*floinf(iriver)*XINF(iriver) 

	      P_Inf = P_Inf + (wqinf(iriver,15) + wqinf(iriver,16) + wqinf(iriver,17) + wqinf(iriver,18) + &
		  wqinf(iriver,19) + wqinf(iriver,20))*floinf(iriver)*XINF(iriver)

!gbs 17Jan06 	 
	      CALL plng_mxng(IRIVER,Bedslope(IRIVER))	
	      IF(plngdpth(IRIVER).lt.0d0) THEN
	         plngdpth(IRIVER)=0.05d0 
!	         print*,'Plung depth negative', iriver,plngdpth(iriver)
!	         pause
	      ENDIF 
!	      plngdpth(IRIVER)=depth(ns) 

	      DSTART(IRIVER)=plngdpth(IRIVER)
	      write(99999,fmt='(3i10,f10.5)')jday,iclock,iriver,plngdpth(iriver)		
!	      print*,'Plung depth', jday, iclock,numinf, iriver,plngdpth(iriver)		
!gbs 17Jan06
20    CONTINUE

!      print*,temperature(2,t3),WQINF(2,1)
!gbs       WRITE(122,*)t3, temperature(1,t3)
!  Work through each element in the downflow stacks and calculate the
!  travel distance and entrainment for the present day, and whether or not
!  it reaches its level of neutral buoyancy and hence can be inserted.

      EINFF = zero
      IRIVER = 1
22    CONTINUE
			IF (IRIVER .gt. NUMINF .OR. FLERR) GOTO 62
24			CONTINUE
			BAN = ALPHA(IRIVER, numseg(iriver))*pi/pideg
		!phi     BSL = PHI(IRIVER)*pi/pideg
!gbs 17Jan06

			BSL=Bedslope(IRIVER)
			CD = CDRAG(IRIVER, numseg(iriver))
			RI = CD*(one+pt21*SQRT(CD)*SIN(BAN))/(SIN(BAN)*SIN(BSL)/COS(BSL))
			E = onept6*(CD**onept5)/RI
			NOINS(IRIVER)=0
			k = 0

!     Loop for calculation of each separate inflow element.
!	   IF(iriver.eq.30) THEN
!	      print*,WQINF(IRIVER,7),WQINF(IRIVER,10), 'step 2'
!     ENDIF
30			k = k+1
			EINF = zero
			IF (k .gt. 1 .AND. FLERR) GOTO 62
			IF (k .gt. ICNT(IRIVER)) GOTO 60
			IF (DDOWN(IRIVER,k) .eq. thsnd) THEN
				TOTIN(IRIVER) = TOTIN(IRIVER)+QDOWN(IRIVER,k)
    !            WRITE(20,FMT='(2I8, F15.5)')JDAY,ICLOCK,TOTIN(IRIVER)
	!         DDOWN(IRIVER,K) = DEPTH(NS)-DSTART(IRIVER) !DDOWN MEASURED FROM BOTTOM Bill 2012/06
	!gbs 17Jan06
				DDOWN(IRIVER,K) = DSTART(IRIVER)
			ENDIF
			TIMEIN = zero
			DOLD(IRIVER,k) = DDOWN(IRIVER,k)
			DI = DENSTY(TDOWN(IRIVER,k),SDOWN(IRIVER,k))

!  Calculate the layer in which this element starts the day.

			DO 35 ln = 1,ns      !1,ns  gbsahoo  2015/9
				IF (depth(ln) .ge. DDOWN(IRIVER,k)) GOTO 36
35			CONTINUE
			ln = ns
36			LAYER = ln

!  Loop for progression of an element through the next layer down.
! Include water quality and particles
			P = PDOWN(IRIVER,K)        
			Q = QDOWN(IRIVER,k)
			T = TDOWN(IRIVER,k)
			S = SDOWN(IRIVER,k)
			CCFF1 = CFDOWN(IRIVER,1,k)
			CCFF2 = CFDOWN(IRIVER,2,k)
			CCFF3 = CFDOWN(IRIVER,3,k)
			CCFF4 = CFDOWN(IRIVER,4,k)
			CCFF5 = CFDOWN(IRIVER,5,k)
			CCFF6 = CFDOWN(IRIVER,6,k)
			CCFF7 = CFDOWN(IRIVER,7,k)
			WQ1 = WQDOWN(IRIVER,1,k)	
			WQ2 = WQDOWN(IRIVER,2,k)
			WQ3 = WQDOWN(IRIVER,3,k)
			WQ4 = WQDOWN(IRIVER,4,k)
			WQ5 = WQDOWN(IRIVER,5,k)
			WQ6 = WQDOWN(IRIVER,6,k)
			WQ7 = WQDOWN(IRIVER,7,k)
			WQ8 = WQDOWN(IRIVER,8,k)
			WQ9 = WQDOWN(IRIVER,9,k)
			WQ10 = WQDOWN(IRIVER,10,k)
			WQ11 = WQDOWN(IRIVER,11,k)
			WQ12 = WQDOWN(IRIVER,12,k)
			WQ13 = WQDOWN(IRIVER,13,k)
			WQ14 = WQDOWN(IRIVER,14,k)
			WQ15 = WQDOWN(IRIVER,15,k)
			WQ16 = WQDOWN(IRIVER,16,k)
			WQ17 = WQDOWN(IRIVER,17,k)
			WQ18 = WQDOWN(IRIVER,18,k)
			WQ19 = WQDOWN(IRIVER,19,k)
			WQ20 = WQDOWN(IRIVER,20,k)
			WQ21 = WQDOWN(IRIVER,21,k)
			WQ22 = WQDOWN(IRIVER,22,k)
			WQ23 = WQDOWN(IRIVER,23,k)
			WQ24 = WQDOWN(IRIVER,24,k)
			WQ25 = WQDOWN(IRIVER,25,k)
			WQ26 = WQDOWN(IRIVER,26,k)
			WQ27 = WQDOWN(IRIVER,27,k)
			WQ28 = WQDOWN(IRIVER,28,k)
			DEP =  DDOWN(IRIVER,k)  !Depth of plunge
			TIM = TIMDOWN(IRIVER,K)            
!  Check IF this element lies below level of neutral buoyancy.
!  If it is call INSERT to insert it and renumber the stack elements.
!  Include water quality and particles			
			IF (DI .le. den(LAYER))THEN
				NOINS(IRIVER)=NOINS(IRIVER)+1  ! shows only 1
				INPAR(IRIVER,NOINS(IRIVER))=k
				TINS(IRIVER,NOINS(IRIVER))=T
				SINS(IRIVER,NOINS(IRIVER))=S
				DIINS(IRIVER,NOINS(IRIVER))=DI
				QINS(IRIVER,NOINS(IRIVER))=Q        
				CFINS(IRIVER,1,NOINS(IRIVER))=CCFF1
				CFINS(IRIVER,2,NOINS(IRIVER))=CCFF2
				CFINS(IRIVER,3,NOINS(IRIVER))=CCFF3
				CFINS(IRIVER,4,NOINS(IRIVER))=CCFF4
				CFINS(IRIVER,5,NOINS(IRIVER))=CCFF5
				CFINS(IRIVER,6,NOINS(IRIVER))=CCFF6
				CFINS(IRIVER,7,NOINS(IRIVER))=CCFF7
				WQINS(IRIVER,1,NOINS(IRIVER))=WQ1
				WQINS(IRIVER,2,NOINS(IRIVER))=WQ2
				WQINS(IRIVER,3,NOINS(IRIVER))=WQ3
				WQINS(IRIVER,4,NOINS(IRIVER))=WQ4
				WQINS(IRIVER,5,NOINS(IRIVER))=WQ5
				WQINS(IRIVER,6,NOINS(IRIVER))=WQ6
				WQINS(IRIVER,7,NOINS(IRIVER))=WQ7
				WQINS(IRIVER,8,NOINS(IRIVER))=WQ8
				WQINS(IRIVER,9,NOINS(IRIVER))=WQ9
				WQINS(IRIVER,10,NOINS(IRIVER))=WQ10
				WQINS(IRIVER,11,NOINS(IRIVER))=WQ11
				WQINS(IRIVER,12,NOINS(IRIVER))=WQ12
				WQINS(IRIVER,13,NOINS(IRIVER))=WQ13
				WQINS(IRIVER,14,NOINS(IRIVER))=WQ14
				WQINS(IRIVER,15,NOINS(IRIVER))=WQ15
				WQINS(IRIVER,16,NOINS(IRIVER))=WQ16
				WQINS(IRIVER,17,NOINS(IRIVER))=WQ17
				WQINS(IRIVER,18,NOINS(IRIVER))=WQ18
				WQINS(IRIVER,19,NOINS(IRIVER))=WQ19
				WQINS(IRIVER,20,NOINS(IRIVER))=WQ20
				WQINS(IRIVER,21,NOINS(IRIVER))=WQ21
				WQINS(IRIVER,22,NOINS(IRIVER))=WQ22
				WQINS(IRIVER,23,NOINS(IRIVER))=WQ23
				WQINS(IRIVER,24,NOINS(IRIVER))=WQ24
				WQINS(IRIVER,25,NOINS(IRIVER))=WQ25
				WQINS(IRIVER,26,NOINS(IRIVER))=WQ26
				WQINS(IRIVER,27,NOINS(IRIVER))=WQ27
				WQINS(IRIVER,28,NOINS(IRIVER))=WQ28
				TIMINS(IRIVER,NOINS(IRIVER)) = TIM
				PINS(IRIVER, NOINS(IRIVER))  = P
				GOTO 30
			ENDIF
!      print*,k
!  Calculate the velocity of the inflow and hence the entrainment.
40			CONTINUE
				GDASH = g*(DI-den(LAYER))/(thsnd+den(LAYER))
!2015/9				HF0   = (two*RI*(Q*COS(BAN)/(SIN(BAN)*tsecda))**2/GDASH)**pt2  !tsectimestep
                HF0   = (two*RI*(Q*COS(BAN)/(SIN(BAN)*tsectimestep))**2/GDASH)**pt2  !tsectimestep
				IF (LAYER .eq. 1) DX = DEP/SIN(BSL)	 
				IF (LAYER .ne. 1) DX = (DEP-depth(LAYER-1))/SIN(BSL)  
			
				HH   = onept2 * E * DX  +  HF0
				VEL1 = Q*thsnd*COS(BAN)/(HH**2*SIN(BAN))
				DELT = DX/VEL1
            !norm_DELT=DELT/FLOAT (nosecs) !Produces Unreasonable results SCT OCT 9 2020
				IF ((TIMEIN + DELT) .gt. one .AND. (.NOT.FLERR)) THEN
				!IF ((TIMEIN + norm_DELT) .gt. one .AND. (.NOT.FLERR)) THEN
					DX   = DX*(one-TIMEIN)/DELT
				!	DX   = DX*(one-TIMEIN)/norm_DELT
					DELT = one-TIMEIN
				!	norm_DELT = one-TIMEIN
					HH   = onept2 * E * DX  +  HF0
				ENDIF
				TIMEIN = TIMEIN+DELT
				!TIMEIN = TIMEIN+norm_DELT
				!TIM=TIM+TIMEIN*(FLOAT(nosecs)/86400.0d0)      
				DELQ = 0.2d0*Q*((HH/HF0)**(five/three)-one)
	
! Check for negative inflow layer volume

				IF(vol(LAYER).lt.0.0)THEN
					WRITE(6,*) 'vol(layer) is negative - layer no. ',layer, ' ns = ',ns
					WRITE(6,*) 'vol(layer) = ',vol(layer), 'layer no.', ns
					PAUSE
				ENDIF

				IF (DELQ .gt. pt9*vol(LAYER)) DELQ=pt9*vol(LAYER)
				S  = COMBIN(S,Q,DI,sal(LAYER),DELQ,den(LAYER ))
				T  = COMBIN(T,Q,DI,temp(LAYER),DELQ,den(LAYER))
				DI = densty(T,S)
!
! Nutrient Budget
!
				N_Ent = N_Ent + (wqual(layer,21) + wqual(layer,22) + wqual(layer,23) +       &
     								  wqual(layer,24) + wqual(layer,25) + wqual(layer,26))*DELQ 
				P_Ent = P_Ent + (wqual(layer,15) + wqual(layer,16) + wqual(layer,17) +      &
									  wqual(layer,18) +wqual(layer,19) + wqual(layer,20))*DELQ
				DO i_p = 1,7
					Part_Ent(i_p) = Part_Ent(i_p) +cf(layer,i_p)*DELQ*1000.0D0
				ENDDO
				CCFF1 = COMBINV(CCFF1,Q,cf(LAYER,1),DELQ)
				CCFF2 = COMBINV(CCFF2,Q,cf(LAYER,2),DELQ)
				CCFF3 = COMBINV(CCFF3,Q,cf(LAYER,3),DELQ)
				CCFF4 = COMBINV(CCFF4,Q,cf(LAYER,4),DELQ)
				CCFF5 = COMBINV(CCFF5,Q,cf(LAYER,5),DELQ)
				CCFF6 = COMBINV(CCFF6,Q,cf(LAYER,6),DELQ)
				CCFF7 = COMBINV(CCFF7,Q,cf(LAYER,7),DELQ)
				WQ1   = COMBINV(WQ1,Q,wqual(LAYER,1),DELQ)
				WQ2   = COMBINV(WQ2,Q,wqual(LAYER,2),DELQ)
				WQ3   = COMBINV(WQ3,Q,wqual(LAYER,3),DELQ)
				WQ4   = COMBINV(WQ4,Q,wqual(LAYER,4),DELQ)
				WQ5   = COMBINV(WQ5,Q,wqual(LAYER,5),DELQ)
				WQ6   = COMBINV(WQ6,Q,wqual(LAYER,6),DELQ)
				WQ7   = COMBINV(WQ7,Q,wqual(LAYER,7),DELQ)
				WQ8   = COMBINV(WQ8,Q,wqual(LAYER,8),DELQ)
				WQ9   = COMBINV(WQ9,Q,wqual(LAYER,9),DELQ)
				WQ10  = COMBINV(WQ10,Q,wqual(LAYER,10),DELQ)
				WQ11  = COMBINV(WQ11,Q,wqual(LAYER,11),DELQ)
				WQ12  = COMBINV(WQ12,Q,wqual(LAYER,12),DELQ)
				WQ13  = COMBINV(WQ13,Q,wqual(LAYER,13),DELQ)
				WQ14  = COMBINV(WQ14,Q,wqual(LAYER,14),DELQ)
				WQ15  = COMBINV(WQ15,Q,wqual(LAYER,15),DELQ)
				WQ16  = COMBINV(WQ16,Q,wqual(LAYER,16),DELQ)
				WQ17  = COMBINV(WQ17,Q,wqual(LAYER,17),DELQ)
				WQ18  = COMBINV(WQ18,Q,wqual(LAYER,18),DELQ)
				WQ19  = COMBINV(WQ19,Q,wqual(LAYER,19),DELQ)
				WQ20  = COMBINV(WQ20,Q,wqual(LAYER,20),DELQ)
				WQ21  = COMBINV(WQ21,Q,wqual(LAYER,21),DELQ)
				WQ22  = COMBINV(WQ22,Q,wqual(LAYER,22),DELQ)
				WQ23  = COMBINV(WQ23,Q,wqual(LAYER,23),DELQ)
				WQ24  = COMBINV(WQ24,Q,wqual(LAYER,24),DELQ)
				WQ25  = COMBINV(WQ25,Q,wqual(LAYER,25),DELQ)
				WQ26  = COMBINV(WQ26,Q,wqual(LAYER,26),DELQ)
				WQ27  = COMBINV(WQ27,Q,wqual(LAYER,27),DELQ)
				WQ28  = COMBINV(WQ28,Q,wqual(LAYER,28),DELQ)
!	PRINT*,WQ7,WQ10
!
! Perform the addition after (previously was before)! 
!
				Q = Q+DELQ
				VOL(LAYER) = VOL(LAYER)-DELQ

!  Calculate energy of inflowing streams.

!2015/9				IF (LAYER .eq. 1) THEN
!2015/9					EINF = EINF+(DI-den(LAYER))*DEP*g*Q/tsecda
!2015/9				ELSE
!2015/9					EINF = EINF+(DI-den(LAYER))*(DEP-depth(LAYER-1))*g*Q/tsecda    !tsectimestep
!2015/9				ENDIF
!  Calculate energy of inflowing streams.

				IF (LAYER .eq. 1) THEN
					EINF = EINF+(DI-den(LAYER))*DEP*g*Q/tsectimestep
				ELSE
					EINF = EINF+(DI-den(LAYER))*(DEP-depth(LAYER-1))*g*Q/tsectimestep    !tsectimestep
				ENDIF
!  Reset the downflow stacks.
!  Include water quality and particles
            TIMDOWN(IRIVER,K)=TIM
            PDOWN (IRIVER,K)= P
				QDOWN(IRIVER,k) = Q
				TDOWN(IRIVER,k) = T
				SDOWN(IRIVER,k) = S
				CFDOWN(IRIVER,1,k)=CCFF1
				CFDOWN(IRIVER,2,k)=CCFF2
				CFDOWN(IRIVER,3,k)=CCFF3
				CFDOWN(IRIVER,4,k)=CCFF4
				CFDOWN(IRIVER,5,k)=CCFF5
				CFDOWN(IRIVER,6,k)=CCFF6
				CFDOWN(IRIVER,7,k)=CCFF7
				WQDOWN(IRIVER,1,k)=WQ1
				WQDOWN(IRIVER,2,k)=WQ2
				WQDOWN(IRIVER,3,k)=WQ3
				WQDOWN(IRIVER,4,k)=WQ4
				WQDOWN(IRIVER,5,k)=WQ5
				WQDOWN(IRIVER,6,k)=WQ6
				WQDOWN(IRIVER,7,k)=WQ7
				WQDOWN(IRIVER,8,k)=WQ8
				WQDOWN(IRIVER,9,k)=WQ9
				WQDOWN(IRIVER,10,k)=WQ10
				WQDOWN(IRIVER,11,k)=WQ11
				WQDOWN(IRIVER,12,k)=WQ12
				WQDOWN(IRIVER,13,k)=WQ13
				WQDOWN(IRIVER,14,k)=WQ14
				WQDOWN(IRIVER,15,k)=WQ15
				WQDOWN(IRIVER,16,k)=WQ16
				WQDOWN(IRIVER,17,k)=WQ17
				WQDOWN(IRIVER,18,k)=WQ18
				WQDOWN(IRIVER,19,k)=WQ19
				WQDOWN(IRIVER,20,k)=WQ20
				WQDOWN(IRIVER,21,k)=WQ21
				WQDOWN(IRIVER,22,k)=WQ22
				WQDOWN(IRIVER,23,k)=WQ23
				WQDOWN(IRIVER,24,k)=WQ24
				WQDOWN(IRIVER,25,k)=WQ25
				WQDOWN(IRIVER,26,k)=WQ26
				WQDOWN(IRIVER,27,k)=WQ27
				WQDOWN(IRIVER,28,k)=WQ28
				TOTIN(IRIVER) = TOTIN(IRIVER) + DELQ
				DDOWN(IRIVER,k) = DEP-DX*SIN(BSL)
				DEP = DDOWN(IRIVER,k)

				IF (DDOWN(IRIVER,k) .lt. DLWST(IRIVER)) THEN
					DLWST(IRIVER) =  DDOWN(IRIVER,k)
				ENDIF
!  If the inflow parcel has reached the level of neutral buoyancy
!  or has reached the bottom, put it in the insertion queue.
!  Include water quality and particles.

				IF (LAYER .eq. 1) GOTO 43
				IF (DI .LE. DEN(LAYER-1)) GOTO 43
				GOTO 45
43				NOINS(IRIVER)=NOINS(IRIVER)+1
				INPAR(IRIVER,NOINS(IRIVER))=k
				TIMINS(IRIVER,NOINS(IRIVER))=TIM
				TINS(IRIVER,NOINS(IRIVER))=T
				SINS(IRIVER,NOINS(IRIVER))=S
				PINS(IRIVER,NOINS(IRIVER))=P			
				CFINS(IRIVER,1,NOINS(IRIVER))=CCFF1
				CFINS(IRIVER,2,NOINS(IRIVER))=CCFF2
				CFINS(IRIVER,3,NOINS(IRIVER))=CCFF3
				CFINS(IRIVER,4,NOINS(IRIVER))=CCFF4
				CFINS(IRIVER,5,NOINS(IRIVER))=CCFF5
				CFINS(IRIVER,6,NOINS(IRIVER))=CCFF6
				CFINS(IRIVER,7,NOINS(IRIVER))=CCFF7
				WQINS(IRIVER,1,NOINS(IRIVER))=WQ1
				WQINS(IRIVER,2,NOINS(IRIVER))=WQ2
				WQINS(IRIVER,3,NOINS(IRIVER))=WQ3
				WQINS(IRIVER,4,NOINS(IRIVER))=WQ4
				WQINS(IRIVER,5,NOINS(IRIVER))=WQ5
				WQINS(IRIVER,6,NOINS(IRIVER))=WQ6
				WQINS(IRIVER,7,NOINS(IRIVER))=WQ7
				WQINS(IRIVER,8,NOINS(IRIVER))=WQ8
				WQINS(IRIVER,9,NOINS(IRIVER))=WQ9
				WQINS(IRIVER,10,NOINS(IRIVER))=WQ10
				WQINS(IRIVER,11,NOINS(IRIVER))=WQ11
				WQINS(IRIVER,12,NOINS(IRIVER))=WQ12
				WQINS(IRIVER,13,NOINS(IRIVER))=WQ13
				WQINS(IRIVER,14,NOINS(IRIVER))=WQ14
				WQINS(IRIVER,15,NOINS(IRIVER))=WQ15
				WQINS(IRIVER,16,NOINS(IRIVER))=WQ16
				WQINS(IRIVER,17,NOINS(IRIVER))=WQ17
				WQINS(IRIVER,18,NOINS(IRIVER))=WQ18
				WQINS(IRIVER,19,NOINS(IRIVER))=WQ19
				WQINS(IRIVER,20,NOINS(IRIVER))=WQ20
				WQINS(IRIVER,21,NOINS(IRIVER))=WQ21
				WQINS(IRIVER,22,NOINS(IRIVER))=WQ22
				WQINS(IRIVER,23,NOINS(IRIVER))=WQ23
				WQINS(IRIVER,24,NOINS(IRIVER))=WQ24
				WQINS(IRIVER,25,NOINS(IRIVER))=WQ25
				WQINS(IRIVER,26,NOINS(IRIVER))=WQ26
				WQINS(IRIVER,27,NOINS(IRIVER))=WQ27
				WQINS(IRIVER,28,NOINS(IRIVER))=WQ28
				DIINS(IRIVER,NOINS(IRIVER))=DI
				QINS(IRIVER,NOINS(IRIVER))=Q
45				CONTINUE

!  If the inflow parcel has ended its days travel, reached its level of
!  Neutral buoyancy or the bottom of the reservoir go to the next parcel.

				IF (LAYER .eq. 1) GOTO 46

!gbs 17Jan06      IF (TIMEIN .ge. one .AND. (.NOT. FLERR)) GOTO 46
				IF (TIMEIN .ge. 0.999999999 .AND. (.NOT. FLERR)) GOTO 46
				IF (DI .le. den(LAYER-1)) GOTO 46
				GOTO 52
46				IF(.NOT. FLERR)EINFF = EINFF+EINF/TIMEIN  

				IF(LAYER .eq. 1)vol1(1)=vol(1)
				IF(LAYER .ne. 1)vol1(LAYER)=vol1(LAYER-1)+vol(LAYER)
      
				DO 50 kk=LAYER+1,ns
					vol1(kk) = vol1(kk-1)+vol(kk)
50				CONTINUE
				GOTO 30
52				CONTINUE
				LAYER = LAYER-1
			GOTO 40
60			CONTINUE
      IRIVER = IRIVER+1
      GOTO 22
62    CONTINUE
!     Insert all of the parcels which reached their level of NB on this day.
!     Adjust the stacking to note the removal.
!     Note that NTIMS should be NOSECS IF the timestep is variable
!     Added to take into account variable time step
      IF(itimes .eq.1440) THEN 
		   NTIMS = 86400
	   ELSE
		   NTIMS = NOSECS
	   ENDIF
      WIDTH = 0.0			
      DO 66 IRIVER=1,NUMINF
			BSL = Bedslope(IRIVER)
		   IRIV = IRIVER
         DO 66 j=1,NOINS(IRIVER)
	         CF1(IRIV,j)=CFINS(IRIV,1,j)
	         CF2(IRIV,j)=CFINS(IRIV,2,j)
	         CF3(IRIV,j)=CFINS(IRIV,3,j)
	         CF4(IRIV,j)=CFINS(IRIV,4,j)
	         CF5(IRIV,j)=CFINS(IRIV,5,j)
	         CF6(IRIV,j)=CFINS(IRIV,6,j)
	         CF7(IRIV,j)=CFINS(IRIV,7,j)
	         WQUAL1(IRIV,j)=WQINS(IRIV,1,j)
	         WQUAL2(IRIV,j)=WQINS(IRIV,2,j)
	         WQUAL3(IRIV,j)=WQINS(IRIV,3,j)
	         WQUAL4(IRIV,j)=WQINS(IRIV,4,j)
	         WQUAL5(IRIV,j)=WQINS(IRIV,5,j)
	         WQUAL6(IRIV,j)=WQINS(IRIV,6,j)
	         WQUAL7(IRIV,j)=WQINS(IRIV,7,j)
	         WQUAL8(IRIV,j)=WQINS(IRIV,8,j)
	         WQUAL9(IRIV,j)=WQINS(IRIV,9,j)
	         WQUAL10(IRIV,j)=WQINS(IRIV,10,j)
	         WQUAL11(IRIV,j)=WQINS(IRIV,11,j)
	         WQUAL12(IRIV,j)=WQINS(IRIV,12,j)
	         WQUAL13(IRIV,j)=WQINS(IRIV,13,j)
	         WQUAL14(IRIV,j)=WQINS(IRIV,14,j)
	         WQUAL15(IRIV,j)=WQINS(IRIV,15,j)
	         WQUAL16(IRIV,j)=WQINS(IRIV,16,j)
	         WQUAL17(IRIV,j)=WQINS(IRIV,17,j)
	         WQUAL18(IRIV,j)=WQINS(IRIV,18,j)
	         WQUAL19(IRIV,j)=WQINS(IRIV,19,j)
	         WQUAL20(IRIV,j)=WQINS(IRIV,20,j)
	         WQUAL21(IRIV,j)=WQINS(IRIV,21,j)
	         WQUAL22(IRIV,j)=WQINS(IRIV,22,j)
	         WQUAL23(IRIV,j)=WQINS(IRIV,23,j)
	         WQUAL24(IRIV,j)=WQINS(IRIV,24,j)
	         WQUAL25(IRIV,j)=WQINS(IRIV,25,j)
	         WQUAL26(IRIV,j)=WQINS(IRIV,26,j)
	         WQUAL27(IRIV,j)=WQINS(IRIV,27,j)
	         WQUAL28(IRIV,j)=WQINS(IRIV,28,j)
!
! Nutrient Budget
!
	         N_Ins =N_Ins + (wqual21(iriv,j) + wqual22(iriv,j)+ wqual23(iriv,j) +         &
     	                      wqual24(iriv,j) + wqual25(iriv,j)+ wqual26(iriv,j))*qins(iriv,j)
	         P_Ins =P_Ins + (wqual5(iriv,j)  + wqual16(iriv,j)+ wqual17(iriv,j)+          &
                            wqual18(iriv,j) + wqual19(iriv,j)+ wqual20(iriv,j))*qins(iriv,j)	 
!
! Particles Budget
!      
	         Part_Ins(1) = Part_Ins(1) + CF1(IRIV,j)*qins(iriv,j)*1000.0D0
	         Part_Ins(2) = Part_Ins(2) + CF2(IRIV,j)*qins(iriv,j)*1000.0D0
	         Part_Ins(3) = Part_Ins(3) + CF3(IRIV,j)*qins(iriv,j)*1000.0D0
	         Part_Ins(4) = Part_Ins(4) + CF4(IRIV,j)*qins(iriv,j)*1000.0D0
	         Part_Ins(5) = Part_Ins(5) + CF5(IRIV,j)*qins(iriv,j)*1000.0D0
	         Part_Ins(6) = Part_Ins(6) + CF6(IRIV,j)*qins(iriv,j)*1000.0D0
	         Part_Ins(7) = Part_Ins(7) + CF7(IRIV,j)*qins(iriv,j)*1000.0D0
				
	         CALL INSERT(E,QINS(IRIV,j),DIINS(IRIV,j),BSL,IRIV,TINS(IRIV,j),SINS(IRIV,j),     &
	         WQUAL1(IRIV,j),WQUAL2(IRIV,j),WQUAL3(IRIV,j),WQUAL4(IRIV,j),WQUAL5(IRIV,j),      &
	         WQUAL6(IRIV,j),WQUAL7(IRIV,j),WQUAL8(IRIV,j),WQUAL9(IRIV,j),WQUAL10(IRIV,j),     &
            WQUAL11(IRIV,j),WQUAL12(IRIV,j),WQUAL13(IRIV,j),WQUAL14(IRIV,j),WQUAL15(IRIV,j), &
            WQUAL16(IRIV,j),WQUAL17(IRIV,j),WQUAL18(IRIV,j),WQUAL19(IRIV,j),WQUAL20(IRIV,j), &
            WQUAL21(IRIV,j),WQUAL22(IRIV,j),WQUAL23(IRIV,j),WQUAL24(IRIV,j),WQUAL25(IRIV,j), &
            WQUAL26(IRIV,j),WQUAL27(IRIV,j),WQUAL28(IRIV,j),CF1(IRIV,j),CF2(IRIV,j),         &
            CF3(IRIV,j),CF4(IRIV,j),CF5(IRIV,j),CF6(IRIV,j),CF7(IRIV,j),NTIMS,WIDTH,LL)

	         DO 65 JK=INPAR(IRIV,j),ICNT(IRIV)-1
	            QDOWN(IRIV,JK) = QDOWN(IRIV,JK+1)
	            TDOWN(IRIV,JK) = TDOWN(IRIV,JK+1)
	            SDOWN(IRIV,JK) = SDOWN(IRIV,JK+1)
	            PDOWN(IRIV,JK) = PDOWN(IRIV,JK+1)
	            TIMDOWN(IRIV,JK)=TIMDOWN(IRIV,JK+1)
	            DO 63 zz=1,28
		            WQDOWN(IRIV,zz,JK) = WQDOWN(IRIV,zz,JK+1)
63             CONTINUE
	            DO 64 zz=1,7
		            CFDOWN(IRIV,zz,JK) = CFDOWN(IRIV,zz,JK+1)
64             CONTINUE
	            DDOWN(IRIV,JK) = DDOWN(IRIV,JK+1)
	            DOLD(IRIV,JK) = DOLD(IRIV,JK+1)
65          CONTINUE 
	         TOTIN(IRIV) = TOTIN(IRIV)-QINS(IRIV,j)
	         ICNT(IRIV) = ICNT(IRIV)-1
	         IF (ICNT(IRIV).eq.0) TOTIN(IRIV)=zero
	         IF (ICNT(IRIV).eq.0) DLWST(IRIV)=thsnd
	         DO 66 k = j,NOINS(IRIVER)
	            INPAR(IRIVER,k)=INPAR(IRIVER,k)-1
66    CONTINUE 

!  Reset the number of insertions per river to be zero.

      IRIV = 0
      DO 67 k=1,MAXINF
	      NOINS(k) = 0.0
67    CONTINUE

!  Calculate the front of the downflow for each river.

      DO 70 IRIVER=1,NUMINF
	      DLWST(IRIVER) = thsnd
	      DO 70 j=1,ICNT(IRIVER)
	         IF (DDOWN(IRIVER,j) .lt. DLWST(IRIVER)) THEN
	            DLWST(IRIVER) = DDOWN(IRIVER,j)
	         ENDIF
70    CONTINUE

      FLERR = .FALSE.			
		CALL NEWSTO(IRIV, FLERR)
!IF	   CALL NEWSTO(IRIV, FLERR, Bedslope)
!  If a flow fit error has occured, go back and reroute first inflow for
!  the river of concern - IRIV.
      IF (FLERR) IRIVER = IRIV
      IF (FLERR) GOTO 24
!  Make adjustments to correct layer volumes.
      CALL RESINT(2,1)
      CALL THICK
!     IF (PINSRT)THEN
!	      DO 115 L=1,NUMINF
!	         WRITE(18,*)' DOWN STACK FOR RIVER',L
!	         WRITE(18,111)
!111        FORMAT(1X,10X,'DDOWN',10X,'QDOWN',10X,'TDOWN',10X,'SDOWN')
!	         DO 115 M=1,ICNT(L)
!	            WRITE(18,112)DDOWN(L,M),QDOWN(L,M),TDOWN(L,M),SDOWN(L,M)
!112           FORMAT(1X,4F15.5)
!115     CONTINUE
!      ENDIF

!
! Nutrient Budget
!
     

      DO i = 1,numinf
	      DO k = 1,icnt(i)
	         N_Stack_F = N_Stack_F + (wqdown(i,21,k) + wqdown(i,22,k) + wqdown(i,23,k) +  &
     		                            wqdown(i,24,k) + wqdown(i,25,k)	+ wqdown(i,26,k))*qdown(i,k)
	         P_Stack_F = P_Stack_F + (wqdown(i,15,k) + wqdown(i,16,k) + wqdown(i,17,k) + &
			 wqdown(i,18,k) + wqdown(i,19,k) + wqdown(i,20,k))*qdown(i,k)
	         DO i_p =1,7
	            Part_Stack_F(i_p) = Part_Stack_F(i_p)+ CFDOWN(i,i_p,k)*qdown(i,k)*1000.0D0 
	         ENDDO ! i_p 
         ENDDO  !  k
	   ENDDO  !  i

	   Balance_N_Stack = (N_Stack_F - N_Stack_I) - N_Inf + N_Ins - N_Ent  
	   Balance_P_Stack = (P_Stack_F - P_Stack_I) - P_Inf + P_Ins - P_Ent
	
	   DO i_p = 1,7 
 	      Balance_Part_Stack(i_p) = (Part_Stack_F(i_p) - Part_Stack_I(i_p))- &
 	                                 Part_Inf(i_p) +  Part_Ins(i_p) - Part_Ent(i_p) 
       ENDDO
!gbs     close(122)
      RETURN
      END SUBROUTINE INFLOW

!**********************************************************************************************	
      SUBROUTINE INSERT(E,q,di,bsl,iriver,t,s,wq1,wq2,wq3,wq4,wq5,wq6,wq7,wq8,wq9,     &
                        wq10,wq11,wq12,wq13,wq14,wq15,wq16,wq17,wq18,wq19,wq20,wq21,   &
                        wq22,wq23,wq24,wq25,wq26,wq27,wq28,ccff1,ccff2,ccff3,ccff4,    &
                        ccff5,ccff6,ccff7,ntims,width,ll)
!**********************************************************************************************
!  This subroutine finds the level of Neutral Bouyancy for a given inflow
!  and returns the layer number (i), the half-thickness (B0), basin length
!  at the intrusion midpoint (AL), basin width (WIDTH) at the intrusion
!  midpoint, and the mean intrusion velocity (UINF) in m/s.
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  

      REAL*8 zero,one,two,three,five,six,pt44,twenty,arfac
      REAL*8 g,secda,thsnd,tsecda,sechr,hrday,pi,pt15

      PARAMETER(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,          &
               pi=3.1415926535897932384626433832795d0,five=5.0d0,    &
               six=6.0d0,pt44=0.44,twenty=20.0,hrday=24.0d0,         &
               arfac=1.0D+6,g=9.81d0,secda=86400.0d0,thsnd=1000.0d0, &
               tsecda=86.4d0,sechr=3600.0d0,pt15=0.15d0)
      REAL*8 AHLE
      REAL*8 AL
      REAL*8 ALSQ
      REAL*8 BL
      REAL*8 B0
      REAL*8 BSL
      REAL*8 CCFF1
      REAL*8 CCFF2
      REAL*8 CCFF3
      REAL*8 CCFF4
      REAL*8 CCFF5
      REAL*8 CCFF6
      REAL*8 CCFF7
!IF      REAL*8 COMBIN,COMBINV,densty
      REAL*8 DBB,DELT,DELB
      REAL*8 DHLE
      REAL*8 DI
      REAL*8 DT
      REAL*8 DVR(MAXNS)
      REAL*8 DZ
      REAL*8 F
      REAL*8 GD
      REAL*8 GR
!IF     REAL*8 GPRIME
	   REAL*8 E
      REAL*8 Q
      REAL*8 R,S
      REAL*8 T,TDASH
      REAL*8 UINF
      REAL*8 VISCOS
      REAL*8 WIDTH
      REAL*8 XN
      REAL*8 XNSQ
      REAL*8 ZP,ZT,ZB
      REAL*8 WQ1
      REAL*8 WQ2
      REAL*8 WQ3
      REAL*8 WQ4
      REAL*8 WQ5
      REAL*8 WQ6
      REAL*8 WQ7
      REAL*8 WQ8
      REAL*8 WQ9
      REAL*8 WQ10
      REAL*8 WQ11
      REAL*8 WQ12
      REAL*8 WQ13
      REAL*8 WQ14
      REAL*8 WQ15
      REAL*8 WQ16
      REAL*8 WQ17
      REAL*8 WQ18
      REAL*8 WQ19
      REAL*8 WQ20
      REAL*8 WQ21
      REAL*8 WQ22
      REAL*8 WQ23
      REAL*8 WQ24
      REAL*8 WQ25
      REAL*8 WQ26
      REAL*8 WQ27
      REAL*8 WQ28
      INTEGER*4 I,j,k,LL
      INTEGER*4 IRIVER,NTIMS,IZ
      INTEGER*4 JB,JT
      INTEGER*4 KX,KY
      INTEGER*4 NB1
      INTEGER*4 NT1
!IF	   INTEGER*4 jday

      DO 100 IZ=1,ns
	      IF(depth(IZ).gt.HLE) GOTO 200
100   CONTINUE
      IZ=ns
200   CONTINUE
      
      DHLE=depth(IZ)
      AHLE=area(IZ)
      IF (DI.LE.den(ns)) THEN
!     Here for surface overflow
	      I=NS
	      LL=I
	      B0=(depth(ns)-depth(ns-1))/two
!2012/06 AL=LC
	      AL=LC-plngdpth(iriver)/TAN(BSL)  !bill 2012/06
!2012/06	IF(WIDTH .le. 1E-7)WIDTH=area(ns)*arfac/AL
         WIDTH=area(ns)*arfac/AL
	      UINF=Q*thsnd/(two*B0*WIDTH)
	      GOTO 2000
      ENDIF

!  Find level of neutral buoyancy

      DO 300 i=ns,2,-1
	      IF(DI.LE.den(i-1)) GOTO 400
300   CONTINUE

!  Here for underflow

      i=1
      LL=i
!2013/04      AL=DHLE/SIN(BSL)
      AL=DHLE/TAN(BSL)
!2012/06 IF(WIDTH .le. 1E-7) WIDTH=AHLE*arfac/AL
      WIDTH=AHLE*arfac/AL
      B0=(DHLE+two)/two
      UINF=Q*thsnd/(two*B0*WIDTH)
      GOTO 2000

!  Here for intrusion

400   CONTINUE
      JT=i
      LL=i
      JB=i-1
!2012/06  AL=depth(i)/SIN(BSL)
      AL=depth(i)/TAN(BSL)
      ALSQ=AL**2
!2012/06      IF(WIDTH .le. 1E-7)WIDTH=area(i)*arfac/AL
      WIDTH=area(i)*arfac/AL
500   CONTINUE
      DT=DEPTHM(JT)
      DBB=DEPTHM(JB)
      DZ=DT-DBB
      XNSQ=g*(den(JB)-den(JT))/((DI+1000.0)*DZ)

!  Here for unstable stratification

      IF (XNSQ.LE. zero)THEN
	      BL=AL
	      GOTO 600
      ENDIF

!  here for stable stratification

      XN		= SQRT(XNSQ)
!2012/06      F		= Q*1000.0/(WIDTH*ntims*XN*ALSQ)
      F	= Q*1000.0/(WIDTH*nosecs*XN*ALSQ)
      VISCOS	= ep(i)*twenty
      IF (VISCOS.LE.0.) VISCOS=VISC
      GR		= XNSQ*ALSQ**2/(VISCOS**2)
      R		= F*GR**(one/three)
!2012/06      TDASH	= ntims*XN/(GR**(one/six))
      TDASH	= nosecs*XN/(GR**(one/six))
!SGS R/TDASH IS A BUG FIX. Some twit tried to normalize the equations
!sgs	about 10yrs ago, and got it wrong!!  It now agrees with 
!sgs	Imberger and Patterson (1981) - except that it is now normalized.
!sgs	Coefficient for viscous CASE changed to 0.57
!sgs  Factor of 2 removed from calculation of tdash in line above

!2013/04     R=R/TDASH
!2013/04      IF (R.gt.1.D0) THEN
!2013/04	      BL=pt44*TDASH*AL*SQRT(R)
!2013/04      ELSE
!2013/04	      BL=0.57*AL*R**(three/two)*(TDASH/R)**(five/six)
!2013/04      ENDIF

      IF (R > TDASH) THEN
	    BL = 0.44*AL*SQRT(R)*TDASH
      ELSE
	    BL = 0.57*AL*(R**(2/3))*(TDASH**(5/6))
      END IF
      
      BL=MIN(BL,AL)
      IF (BL .lt. 1.0) BL = 1.0
600   CONTINUE
!  B0 is 1/2 the intrusion thickness
!2013/04     B0=Q*thsnd/(WIDTH*BL)
      B0 = Q/(WIDTH*BL)*(1-BL/AL)*2   
      IF(B0.LE.DZ) GOTO 800
      IF((JT.eq.ns).AND.(JB.eq.1)) GOTO 800
      IF(JT.eq.ns) GOTO 700
      JT=JT+1
700   IF(JB.eq.1) GOTO 500
      JB=JB-1
      GOTO 500
800   CONTINUE

      IF (depth(i) .lt. (DHLE + one))THEN
	      AL=DHLE/SIN(BSL)	      
	      IF (WIDTH .le. 1E-7)WIDTH=AHLE*arfac/AL
	      B0=(DHLE+two)/two
      ENDIF
!2013/04      UINF=Q*thsnd/(two*B0*WIDTH)
      UINF=(Q/nosecs)/(two*B0*WIDTH)

!  Mix the inflow with the appropriate layers.

2000  NT1=i
      ZP=DEPTHM(i)
      IF (I.GT.1.AND.I.LT.NS) GOTO 2100

!  Here for underflow, overflow, or fully mixed

      NB1		= i
      DVR(i)	= Q
      GOTO 3000

!  Here for intrusion

2100  NT1=NT1+1
      IF (NT1.eq.ns) GOTO 2200
      GD=GPRIME(den(NT1),den(i))
      IF(GD.gt.zero) DELT=pt15*(UINF/float(ntims))**2/GD
      IF(DELT.gt.(DEPTHM(NT1)-ZP) .OR. GD .le. zero) GOTO 2100
      IF(depth(NT1-1) .gt. (ZP+DELT))NT1=NT1-1
2200  CONTINUE
      ZT=depth(NT1)
      NB1=i
2300  CONTINUE
      NB1=NB1-1
      IF (NB1.eq.1) GOTO 2400
      GD=GPRIME(den(i),den(NB1))
      IF(GD.gt.zero) DELB=pt15*((UINF/float(ntims))**2)/GD
      IF(DELB.gt.(ZP-DEPTHM(NB1)) .OR. GD .le. zero) GOTO 2300
      IF(depth(NB1) .lt. (ZP-DELB))NB1=NB1+1
2400  CONTINUE
      ZB=zero
      IF (NB1.gt.1) ZB=depth(NB1-1)
      IF (NB1 .eq. NT1) THEN
!  Here IF intrusion is entirely within layer i
	      DVR(NB1) = Q
	      GOTO 3000
      ENDIF
!  Aportion inflow amongst layers NB1,NB1+1,---,NT1
      DELT=ZT-ZP
      DELB=ZP-ZB
      IF (NB1.eq.i) GOTO 2600
	   DVR(NB1)=Q*(depth(NB1)-ZB+DELB*SIN(pi*(depth(NB1)-ZP)/DELB)/pi)/(ZT-ZB)
	   IF (NB1.eq.(i-1)) GOTO 2600
	   KX=NB1+1
	   KY=i-1
	   DO 2500 k=KX,KY
	      DVR(k)=Q*(depth(k)-depth(k-1)+DELB*(SIN(pi*(ZP-depth(k-1))/DELB)-   &
                SIN(pi*(ZP-depth(k))/DELB))/pi)/(ZT-ZB)
2500  CONTINUE
2600  CONTINUE
      DVR(i)=Q*(ZP-depth(i-1)+DELB*SIN(pi*(ZP-depth(i-1))/DELB)/pi)/(ZT-ZB)
      DVR(i)=DVR(i)+Q*(depth(i)-ZP+DELT*SIN(pi*(depth(i)-ZP)/DELT)/pi)/(ZT-ZB)
      IF (NT1.eq.i) GOTO 3000
      IF (NT1.eq.(i+1)) GOTO 2800
      KX=i+1
      KY=NT1-1
      DO 2700 k=KX,KY
	      DVR(k)=Q*(depth(k)-depth(k-1)+DELT*(SIN(pi*(depth(k)-ZP)/DELT)-        &
                SIN(pi*(depth(k-1)-ZP)/DELT))/pi)/(ZT-ZB)
2700  CONTINUE
2800  CONTINUE
      DVR(NT1)=Q*(ZT-depth(NT1-1)+DELT*SIN(pi*(ZP-depth(NT1-1))/DELT)/pi)/(ZT-ZB)

!  Insert inflow into reservoir and adjust layer properties
!  Include water quality and particles

3000 CONTINUE
!     print*,nb1,nt1, depth(nb1),depth(nt1), ns, depth(ns), depth(1), temp(1)
      DO 3200 k=NB1,NT1
	      temp(k)=COMBIN(temp(k),vol(k),den(k),T,DVR(k),DI)
	      sal(k)=COMBIN(sal(k),vol(k),den(k),S,DVR(k),DI)
	      cf(k,1)=COMBINV(cf(k,1),vol(k),CCFF1,DVR(k))
	      cf(k,2)=COMBINV(cf(k,2),vol(k),CCFF2,DVR(k))
	      cf(k,3)=COMBINV(cf(k,3),vol(k),CCFF3,DVR(k))
	      cf(k,4)=COMBINV(cf(k,4),vol(k),CCFF4,DVR(k))
	      cf(k,5)=COMBINV(cf(k,5),vol(k),CCFF5,DVR(k))
	      cf(k,6)=COMBINV(cf(k,6),vol(k),CCFF6,DVR(k))
	      cf(k,7)=COMBINV(cf(k,7),vol(k),CCFF7,DVR(k))
	      wqual(k,1)=COMBINV(wqual(k,1),vol(k),WQ1,DVR(k))
	      wqual(k,2)=COMBINV(wqual(k,2),vol(k),WQ2,DVR(k))
	      wqual(k,3)=COMBINV(wqual(k,3),vol(k),WQ3,DVR(k))
	      wqual(k,4)=COMBINV(wqual(k,4),vol(k),WQ4,DVR(k))
	      wqual(k,5)=COMBINV(wqual(k,5),vol(k),WQ5,DVR(k))
	      wqual(k,6)=COMBINV(wqual(k,6),vol(k),WQ6,DVR(k))
	      wqual(k,7)=COMBINV(wqual(k,7),vol(k),WQ7,DVR(k))
	      wqual(k,8)=COMBINV(wqual(k,8),vol(k),WQ8,DVR(k))
	      wqual(k,9)=COMBINV(wqual(k,9),vol(k),WQ9,DVR(k))
	      wqual(k,10)=COMBINV(wqual(k,10),vol(k),WQ10,DVR(k))
	      wqual(k,11)=COMBINV(wqual(k,11),vol(k),WQ11,DVR(k))
	      wqual(k,12)=COMBINV(wqual(k,12),vol(k),WQ12,DVR(k))
	      wqual(k,13)=COMBINV(wqual(k,13),vol(k),WQ13,DVR(k))
	      wqual(k,14)=COMBINV(wqual(k,14),vol(k),WQ14,DVR(k))
	      wqual(k,15)=COMBINV(wqual(k,15),vol(k),WQ15,DVR(k))
	      wqual(k,16)=COMBINV(wqual(k,16),vol(k),WQ16,DVR(k))
	      wqual(k,17)=COMBINV(wqual(k,17),vol(k),WQ17,DVR(k))
	      wqual(k,18)=COMBINV(wqual(k,18),vol(k),WQ18,DVR(k))
	      wqual(k,19)=COMBINV(wqual(k,19),vol(k),WQ19,DVR(k))
	      wqual(k,20)=COMBINV(wqual(k,20),vol(k),WQ20,DVR(k))
	      wqual(k,21)=COMBINV(wqual(k,21),vol(k),WQ21,DVR(k))
	      wqual(k,22)=COMBINV(wqual(k,22),vol(k),WQ22,DVR(k))
	      wqual(k,23)=COMBINV(wqual(k,23),vol(k),WQ23,DVR(k))
	      wqual(k,24)=COMBINV(wqual(k,24),vol(k),WQ24,DVR(k))
	      wqual(k,25)=COMBINV(wqual(k,25),vol(k),WQ25,DVR(k))
	      wqual(k,26)=COMBINV(wqual(k,26),vol(k),WQ26,DVR(k))
	      wqual(k,27)=COMBINV(wqual(k,27),vol(k),WQ27,DVR(k))
	      wqual(k,28)=COMBINV(wqual(k,28),vol(k),WQ28,DVR(k))
3150     CONTINUE 
	   den(k)=densty(temp(k),sal(k))
	   vol(k)=vol(k)+DVR(k)

!quim CSGS Bubflag is introduced to prevent the inflow debug file writing out
!SGS the insertion data when the bubbler is on
!    WRITE out intrusion depths to the bubbler debug file
!    and WRITE intrusion parameters, too
!      IF(PINSRT)THEN
!	 WRITE(18,780) jday 
! 780   format('Insertion for day',I10)
!       WRITE(18,790)
! 790	 format(3x,'Layer',3X,'Depth',6X,'Volume',10X,'River')
!      ENDIF

!	IF(PINSRT.and.iriver .eq.9)
!     &	WRITE(18,9000)jday,iriver,k,depth(k),DVR(k),E
!9000  format(I10,4x,I4,I6,4X,F8.2,x,f12.1,xf10.3)

!-------------Activate when needed -----------------------------------
!	   IF(iclock.eq.itmpr.or.(itimes.eq.1440.and.itmpr.eq.86400)) THEN
	      IF(PINSRT)WRITE(18,9000)SIMDAY,nosecs+iclock,IRIVER,K,DEPTH(K),DVR(K), NB1,NT1, NS
9000     FORMAT(i4,5x,i6,5x,I6,5x,I6,3X,F10.2,3X,f12.1, 3I8)
!      ENDIF
3200  CONTINUE
      
      vol1(1)=vol(1)
      IF (ns.eq.1) GOTO 3400
      DO 3300 j=2,ns
	    vol1(j)=vol1(j-1)+vol(j)
3300  CONTINUE
      CALL RESINT(2,1)
      CALL THICK
3400  CONTINUE

      RETURN
      END SUBROUTINE INSERT
!*************************************************************************
      SUBROUTINE FLDP (HF, DL, vol, depth, ALPHABED, FLERR, bsl)
!phi  SUBROUTINE FLDP (HF, DL, vol, depth, PHI, ALPHABED, FLERR)
!***********************************************************************
! Subroutine to solve the cubic for flowing depth in a triangular river valley.
!---------------------------------------------------------------------------------
!      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  
      
      REAL*8 ARG,HF,depth,DL,vol,ALPHABED,THETA,TAG,BSL,BAN,CON
!phi	REAL*8 ARG,HF,depth,DL,vol,PHI,ALPHABED,THETA,H,TAG,BSL,BAN,CON
      REAL*8 zero,pt25,one,two,three,four,six,nine,pi,pideg,thsnd

      PARAMETER (zero=0.0,one=1.0,two=2.0,three=3.0,four=4.0,nine=9.0)
      PARAMETER (six=6.0,pi=3.141593,pt25=0.25)
      PARAMETER (pideg=180.0,thsnd=1000.0)

      LOGICAL*4 FLERR
!
! Set the coefficients of the cubic.
!
      H = zero
      IF (vol .le. zero) GOTO 30
!phi	BSL = PHI*pi/pideg
      BAN = ALPHABED*pi/pideg
      TAG = TAN(BAN)/TAN(BSL)
10    CONTINUE
	    CON = depth-DL
	    ARG = ABS(one-six*vol*thsnd/TAG/CON**3)
	    IF (ARG.gt.one)GOTO 20
	    THETA = ACOS(one-six*vol*thsnd/TAG/CON**3)
	    H = (two*COS(THETA/three)+one)*CON/two
	    IF (H .gt. zero .AND. H .lt. CON) GOTO 30
	    H = (two*COS(THETA/three+two*pi/three)+one)*CON/two
	    IF (H .gt. zero .AND. H .lt. CON) GOTO 30
		 H = (two*COS(THETA/three+four*pi/three)+one)*CON/two
	    IF (H .gt. zero .AND. H .lt. CON) GOTO 30
20     CONTINUE
	    DL = DL-pt25
	    IF (DL .le. zero) THEN
		   FLERR = .TRUE.
	    ENDIF
	    IF (.NOT. FLERR) GOTO 10
!
! Then the flowing depth for this river is H.
!
30    CONTINUE
      HF = H

      RETURN
      END SUBROUTINE FLDP
!******************************************************************************************
      REAL*8 FUNCTION EXVOL (D1,D2,NRIV,SURDEP)	                
!phi  REAL*8 FUNCTION EXVOL (D1,D2,DLOW,NRIV,SURDEP,PHI,HFLOW,ALPHA)
!******************************************************************************************
!  Function to compute the volume in the inflow stacks which lies between depths D1 and D2.
!---------------------------------------------------------------------------------------------
		USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE 


      INTEGER*4 NRIV

	   REAL(8), PARAMETER :: pi=3.1415926535897932384626433832795d0
	   REAL(8), PARAMETER :: zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0
	   REAL(8), PARAMETER :: pt5=0.5d0
	   REAL(8), PARAMETER :: pideg=180.0d0,thsnd=1000.0d0

!IF	   REAL(8) ::  ALPHABED(*)
	   REAL(8) ::  BAN
	   REAL(8) ::  BSL
	   REAL(8) ::  D1
	   REAL(8) ::  D2
	   REAL(8) ::  DTOP
	   REAL(8) ::  DBOT
!IF	REAL(8) ::  DLWST(*)
!IF	REAL(8) ::  HFLOW(*)
!gbs	REAL(8) ::  PHI(*)
	   REAL(8) ::  SUM
	   REAL(8) ::  SURDEP
	   REAL(8) ::  TAG
	   REAL(8) ::  ZB

	   INTEGER*4 k

      SUM = zero
      DO 10 k=1,NRIV
	      ZB = DLWST(k)
	      IF ((D2 .le. ZB) .OR. (SURDEP .le. DLWST(k))) GOTO 10
			BSL=Bedslope(k)
!phi	    BSL = PHI(k)*pi/pideg
	      BAN = ALPHA(K,numseg(K))*pi/pideg
!gbs		BAN = ALPHA(k, maxseg)*pi/pideg			
	      TAG = TAN(BAN)/TAN(BSL)			
	      DBOT = D1
	      IF(D1 .lt. ZB) DBOT = ZB
	      DTOP = D2
	      IF(D2.gt.SURDEP)DTOP=SURDEP
	      IF (DTOP .lt. DLWST(k)+HFLOW(k)) THEN
	         SUM = SUM+TAG*((DTOP-ZB)**3-(DBOT-ZB)**3)/three
	         GOTO 10
	      END IF
	      IF (DBOT .gt. DLWST(k)+HFLOW(k))THEN
	         SUM = SUM+TAG*HFLOW(k)**2*(DTOP-DBOT)
	         GOTO 10
	      END IF
	      SUM = SUM+TAG*(HFLOW(k)**3-(DBOT-ZB)**3)/three+TAG*HFLOW(k)**2*(DTOP-DLWST(k)-HFLOW(k))
10    CONTINUE
      EXVOL = SUM/thsnd
      RETURN
      END FUNCTION EXVOL
!**************************************************************************
!IF 		SUBROUTINE NEWSTO(IRIV, FLERR, Bedslope)
		SUBROUTINE NEWSTO(IRIV, FLERR)
!***************************************************************************
!  Calculate the new temporary storage table for use in RESINT
!
!----------------------------------------------------------------------------------
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  

	   INTEGER*4 ITOP,J,K,IRIV,jsto
	   REAL*8:: zero,one,ten,pt1,thsnd
	   PARAMETER(zero=0.0,one=1.0,ten=10.0,pt1=0.1,thsnd=1000.0)
!IF	   REAL*8:: Bedslope(maxinf)	! bedslope for each stream
	   REAL*8:: bsl
	   REAL*8:: D
	   REAL*8:: DLOW
	   REAL*8:: DOL
	   REAL*8:: EXTRA
!IF	   REAL*8:: EXVOL
	   REAL*8:: LAYPRO
	   REAL*8:: LAYVOL
	   REAL*8:: VLSUM
	   REAL*8:: VOLSUM
	   REAL*8:: VVTEMP
	   LOGICAL*4 FLERR	  
	   REAL*8:: Alphabed
! Must begin by calculating the new depth of the entire reservoir.
!
	   VLSUM = zero   	
	   DO K=1,NUMINF
		   VLSUM = VLSUM+TOTIN(K)
	   END DO	
	   VOLSUM = VLSUM+VOL1(NS)	! total volume of reservoir = vol1(ns) plus downflow stack   		
	   DO j = 1, numsto
		   IF(VV(j) .gt. VOLSUM) THEN
			   jsto = j - 1		! jsto points to storage table entry just below current volume
			   exit
		   ENDIF
	   ENDDO
	   j = MIN(jsto, numsto)	! point to last entry IF loop is exited without finding vv > volsum
!GBS   j = MAX(jsto, numsto)	! point to last entry IF loop is exited without finding vv > volsum
!       PRINT*, jday, JSTO, NUMSTO
	   DEPTH(NS) = (FLOAT(J)+(VOLSUM - VV(J))/DVV(J))*pt1
	   
!
! Calculate the flowing depths of the inflowing streams.
!
	   DLOW=thsnd
	   DO K = 1,NUMINF
		   alphabed = alpha(K,numseg(K))
		   bsl = bedslope(k)
		   CALL FLDP (HFLOW(K),DLWST(K),TOTIN(K),DEPTH(NS),Alphabed, FLERR, BSL)
! wef 15dec05     &	ALPHA(K,maxseg), FLERR, bsl)!
! If flow error THEN go back to INFLOW
!
		   IF (FLERR) THEN
			   IRIV = K
			   RETURN
		   ENDIF		
		   IF (DLOW.GT.DLWST(K))DLOW=DLWST(K)
   		
		   IF (HFLOW(K) .GT. DEPTH(NS)) THEN
			   STOP ' error in height of flow'
		   ENDIF
	   ENDDO
!
!  Account for the extra volumes of the downflows.
!
	   ITOP = NINT(DEPTH(NS)*ten)
	   EXTRA = zero
	   DOL = zero
	   DO J = 1,ITOP
		   D = FLOAT(J)/ten
		   EXTRA = EXVOL(DOL,D,NUMINF,DEPTH(NS))			
		   IF (J.EQ.1) THEN
			   LAYVOL=VV(J)
		     ELSE
			   LAYVOL=VV(J)-VV(J-1)
		   ENDIF

		   LAYPRO = (FLOAT(J)/ten-DLOW)/(FLOAT(J)/ten)
		   IF(LAYPRO.LT.zero) LAYPRO = zero
		   VVTEMP = zero
		   IF(J.NE.1) VVTEMP = VVDASH(J-1)
		   VVDASH(J) = VV(J) - EXTRA
   		
		   IF (VVDASH(J).LT.(VVTEMP+(one-LAYPRO)*LAYVOL))THEN
			   VVDASH(J) = VVTEMP+(one-LAYPRO)*LAYVOL
		   ENDIF
   		
	   ENDDO ! j

	   DO J=ITOP+1,NUMSTO
		   VVDASH(J) = VV(J)-VLSUM         
	   ENDDO

	   DO J=1,NUMSTO-1
		   DVVDA(J) = VVDASH(J+1)-VVDASH(J)
	   ENDDO
   	
	   DVVDA(NUMSTO) = DVVDA(NUMSTO-1)

	   RETURN
      END SUBROUTINE NEWSTO

! ***************************************************************************************
      SUBROUTINE OUTFLO(hh,flow,extemp,exsalt,exwq,exift,exifs,exifwq,wdl,excf,exifcf)
!****************************************************************************************
!  Removes the outflow at level N
!  VARIABLE LIST
!  AL           RESERVOIR LENGTH AT HEIGHT OF WITHDRAWAL
!  DASNK        LEVEL ABOVE W.L. FROM WHICH FLUID IS DRAWN
!  DEL          INITIAL HALF WITHDRAWAL LAYER THICKNESS
!  DELT         UPPER HALF WITHDRAWAL LAYER THICKNESS
!  DELB         LOWER HALF WITHDRAWAL LAYER THICKNESS
!  DS           DENSITY DIFFERENCE
!  DV           ARRAY OF WITHDRAWAL VOLUMES FROM EACH LAYER
!  DZ           WIDTH OVER WHICH DENSITY GRADIENT IS CALCULATED
!  F            FROUDE NUMBER
!  FLOW         WITHDRAWAL VOLUME
!  FLRT         FLOW RATE
!  GR           GRASHOF NUMBER
!  HH           HEIGHT OF OUTLET
!  IB           LAYER INDEX FOR WITHDRAWAL LAYER BOTTOM
!  ISINK        INDEX
!  IT           LAYER INDEX FOR WITHDRAWAL LAYER TOP
!  N            LAYER INDEX FOR OUTLET LEVEL
!  R            OUTFLOW PARAMETER
!  TSECDA       TIMESTEP SIZE IN SECS./1000
!  UMAX         WITHDRAWAL FLUX
!  W            RESERVOIR WIDTH AT HEIGHT OF WITHDRAWAL
!  XN           BRUNT-VAISALA FREQUENCY
!  XNSQ         B-V FREQUENCY SQUARED
!  Z1           HEIGHT OF LAYER ABOVE OUTLET
!  Z2           HEIGHT OF LAYER BELOW OUTLET
!  ZT           HEIGHT OF TOP OF WITHDRAWAL LAYER
!  ZB           HEIGHT OF BOTTOM OF WITHDRAWAL LAYER

      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  

      REAL(8), PARAMETER :: zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0
	   REAL(8), PARAMETER :: six=6.0d0,twenty=20.0d0,arfac=1.0d6
      REAL(8), PARAMETER :: twopt9=2.9d0,tsecda=86.4d0
	   REAL(8), PARAMETER :: thsnd=1000.0d0,qcheck=0.0001d0,g=9.81d0
      REAL(8) ::  AL,AVDEL
!IF      REAL(8) ::  COMBIN,COMBINV,densty,DV1
      REAL(8) ::  DASNK,DEL,DELB,DELT,DS,DT,DV(MAXNS),DZ
	   REAL(8) ::  EXCF(7),EXIFCF(7)
      REAL(8) ::  EXTEMP,EXIFT,EXSALT,EXIFS,EXVL,EXDEN
      REAL(8) ::  EXWQ(28),EXIFWQ(28)
      REAL(8) ::  FLOW,FLRT,F2,F3
      REAL(8) ::  GR
      REAL(8) ::  HH
      REAL(8) ::  QSTAR
      REAL(8) ::  R
!IF      REAL(8) ::  SQR
      REAL(8) ::  TEL
      REAL(8) ::  VISCOS
      REAL(8) ::  W,WDL(2)
      REAL(8) ::  XN,XNSQ
      REAL(8) ::  ZB,ZT

      INTEGER*4 i,IB,ILP,ISINK,IT,N,zz
      INTEGER*4 LAYCO
      INTEGER*4 NNNCO,NX
	   REAL(8) ::    caca_part_I(7),caca_part_F(7)
           
!  Initialize variables
      DEL=0.0
      ILP=0
      DZ=0.0
      EXTEMP=0.0
      EXSALT=0.0
      EXIFT=0.0
      EXIFS=0.0
      DO 7 zz=1,28
	      EXWQ(zz)=0.0
	      EXIFWQ(zz)=0.0
7     CONTINUE
	   DO 8 zz=1,7
	      EXCF(zz)=0.0
	      EXIFCF(zz)=0.0
8     CONTINUE
      EXVL=0.0     
      EXDEN=0.0
      WDL(1)=0.0
      WDL(2)=0.0
      DO 10 i=1,ns
	      DV(i)=0.0
10    CONTINUE

!  FIND NUMBER OF LAYER (N) OPPOSITE OFFTAKE

      DO 1 i=1,ns
	      IF(depth(i) .ge. HH) GOTO 2
1     CONTINUE
!  return OF RESERVOIR SURFACE IS BELOW OUTLET LEVEL

      RETURN
2     CONTINUE
      N = i

!  CALCULATE THE PROPERTIES AT THE LEVEL OF THE OUTLET
!  INCLUDE WATER QUALITY AND PARTICLES

      EXIFT=temp(N)
      EXIFS=sal(N)
      DO 14 zz=1,28
	      EXIFWQ(zz)=wqual(N,zz)
14    CONTINUE
	   DO 15 zz=1,7
	      EXIFCF(zz)=cf(N,zz)
15    CONTINUE

!  DETERMINE OFFTAKE LEVEL AND FLOW
      DO 20 i=1,NUMOUT
	      IF (HH .eq. OLEV(i)) THEN
	         AL = OLEN(i)
	         W = OWID(i)
	      ENDIF
20    CONTINUE
      
      IF (HH .eq. CRL) THEN
	      AL = LC
	      W = WC
      ENDIF				!TECDA = NOMBRE DE SEGONS PER DIA
!2015/9      FLRT = FLOW/TSECDA	!TSECDA = 86.4 PERQUE FLOW EN 1E6 LITRES/DAY
!2015/9      FLRT = FLOW/(REAL(OUTFLO_TIMESTEP)/1000.0D0) !FLOWRATE SO FLOW*1000/NOSECS
      FLRT = FLOW/(REAL(NOSECS)/1000.0D0) !FLOWRATE SO FLOW*1000/NOSECS
      IF (FLOW .le. QCHECK) RETURN

!  CALCULATE WITHDRAWAL LAYER THICKNESS

      IF (ns .ne. 1) GOTO 36
      WDL(1) = depth(ns)
      WDL(2) = 0.0
      DV(1)=FLOW
      GOTO 190
36    CONTINUE
      ISINK=1
      NX=N+1
40    IF (NX .gt. ns) NX=ns
      IF (NX .lt. 1) NX=1
      DS = den(N)-den(NX)
      IF (DS*float(ISINK) .gt. zero) GOTO 44

!  CHECKS FOR ZERO GRADIENT

46    IF (ISINK.eq.1) THEN
	      DEL=depth(ns)-HH
	      GOTO 54
      ELSE
	      DEL=HH
	      GOTO 55
      ENDIF
44    CONTINUE
      DZ = DEPTHM(NX)-HH
      IF ((DZ*float(ISINK)).LE.zero) GOTO 46
      XNSQ = DS*g/((thsnd+den(N))*DZ)
      IF (XNSQ .le. zero) GOTO 46
      XN = SQRT(XNSQ)
      VISCOS = ep(N)*twenty
      IF(VISCOS .le. zero) VISCOS=VISC
      GR=XNSQ*SQR(AL)*SQR(AL)/SQR(VISCOS)
      IF (ILP.eq.1) GOTO 48

!  POINT SINK CASE

      F3=FLRT/(XN*AL*SQR(AL))
      DEL=AL*F3**(one/three)
      IF ((two*DEL .gt. W).AND.(ISINK .eq. 1)) GOTO 48
      R=F3*SQRT(GR)
      IF(R .le. one)DEL=twopt9*((VISCOS*VISC)**(one/six)*((AL/XN)**(one/three)))
      GOTO 50

!  LINE SINK CASE

48    CONTINUE
      F2=FLRT/W/XN/SQR(AL)
      ILP=1
      DEL=two*AL*SQRT(F2)
      R=F2*GR**(one/three)
      IF(R .le. one)DEL=two*AL/GR**(one/six)

!  CHECK THAT THE WITHDRAWAL LAYER DOES NOT INTERSECT THE TOP OR BOTTOM
!  OF THE RESERVOIR, AND THAT DEL IS LESS THAN THE RANGE OVER WHICH
!  THE DENSITY GRADIENT IS CALCULATED.

50    CONTINUE
      IF (ISINK.NE.1) GOTO 55
54    IF (DEL.LE.ABS(DZ).OR.NX.eq.ns) THEN
	      DELT=DEL
	      IF(DELT.gt.(depth(ns)-HH)) DELT=depth(ns)-HH
	      ISINK=-1
	      NX=N-1
	      GOTO 40
      ELSE
	      ILP=0
	      NX=NX+1
	      GOTO 40
      ENDIF
55    CONTINUE
      IF (DEL.LE.ABS(DZ).OR.NX.eq.1) THEN
	      DELB = DEL
	      IF (DELB.gt.HH-HLE) DELB=HH-HLE 
	      GOTO 60
      ELSE
	      NX=NX-1
	      GOTO 40
      ENDIF

!  CALCULATE TOP AND BOTTOM OF W.L., AND DISTANCE ABOVE
!  W.L. FROM WHICH FLUID IS DRAWN.

60    CONTINUE
!gbs      ZT = HH+DELT-(depth(ns)-depth(ns-1))  !2015/10
       ZT = HH+DELT      
!gbs      ZB = HH-DELB     !2015/9 
!gbs      ZT=DEPTH(NS)
      ZB=HH
!gbs      TEL = DELT+DELB 
      TEL = DELT
      DASNK = FLOW/W/AL*thsnd
      AVDEL = TEL/two
      ZT = ZT+DASNK
      IF (ZT.gt.depth(ns)) ZT = depth(ns)
      WDL(1) = ZT
      WDL(2) = ZB

!  FIND THE INDICES OF THE LAYERS CONTAINING
!  ZT AND ZB. LOCATE IB

      DO 80 i=1,N
	      IF(depth(i) .ge. ZB) GOTO 81
80    CONTINUE

!  IF ZB HIGHER THAN THE BOTTOM OF THE SINK, NOTE AN ERROR.

      WRITE(*,*)' Error OUTFLO - bottom of WL above outlet bottom'
81    CONTINUE
      IB = i

!  LOCATE IT.

      DO 90 i=N,ns
	      IF (depth(i) .ge. ZT) GOTO 91
90    CONTINUE
91    CONTINUE
      IT = i

!  IF ALL DRAWN FROM ONE LAYER.

      IF (IB.eq.IT) THEN
	      WDL(1)=depth(N)
	      IF (N .ne. 1) WDL(2)=depth(N-1)
	      DV(N)=FLOW
	      GOTO 190
      ENDIF

!  CALCULATE THE DV(i), THE PORTION OF FLUID DRAWN FROM THE ITH LAYER.

      DO 100 i=IB,IT
	      DB=zero
	      IF (i .ne. 1) DB=depth(i-1)
	      IF (i .eq. IB)DB=ZB
	      DT=depth(i)
	      IF (i .eq. IT)DT=ZT
	      DV(i) = DV1(DB,DT,DASNK,AVDEL,HH,DELT,DELB)
100   CONTINUE

!  PROPORTION DRAWN FROM EACH LAYER IS KNOWN. MATCH THE
!  TOTAL VOLUME TO THE ACTUAL VOLUME DRAWN OUT.

      QSTAR=zero
      DO 130 i=1,ns   !1,ns  gbsahoo 2015/9
	      QSTAR=QSTAR+DV(i)
130   CONTINUE
      DO 140 i=1,ns !1,ns   gbsahoo 2015/9
	      DV(i)=FLOW/QSTAR*DV(i)
140   CONTINUE

!  CORRECTION IF ANY LAYER SHOULD BE EMPTIED.

      DO 150 i=1,ns-1               !1,ns-1    gbsahoo 2015/9
	      IF (DV(i)-vol(i) .ge. zero) THEN
	         DV(i+1) = DV(i+1)+DV(i)-vol(i)+one
	         DV(i) = vol(i)-one
	      ENDIF
150   CONTINUE

!  NOW HAVE DV(i) FOR ALL LAYERS AND CAN REMOVE IT.

190   CONTINUE
      EXVL =0.0
      DO 200 i=1,ns !1,ns      gbsahoo 2015/9
	      IF (DV(i) .lt. 1E-5) GOTO 200
	      vol(i) = vol(i)-DV(i)

!  Calculate the properties of the withdrawn water.
!  Include water quality and particles

	      IF (EXVL .lt. 1E-5) THEN
	         EXTEMP = temp(i)
	         EXSALT = sal(i)
	         DO 193 zz=1,28
		         EXWQ(zz) = wqual(i,zz)
193         CONTINUE
		      DO 194 zz=1,7
	            EXCF(zz) = cf(i,zz)
194         CONTINUE
	         EXDEN = den(i)
	      ELSE
	         EXTEMP = COMBIN(EXTEMP,EXVL,EXDEN,temp(i),DV(i),den(i))
	         EXSALT = COMBIN(EXSALT,EXVL,EXDEN,sal(i),DV(i),den(i))
	         DO 195 zz=1,28
	            EXWQ(zz)=COMBINV(EXWQ(zz),EXVL,wqual(i,zz),DV(i))
195         CONTINUE
	         DO 196 zz=1,7
	            EXCF(zz)=COMBINV(EXCF(zz),EXVL,cf(i,zz),DV(i))
196         CONTINUE
	         EXDEN = densty(EXTEMP,EXSALT)
	      ENDIF
	      EXVL = EXVL+DV(i)
200   CONTINUE
      vol1(1) = vol(1)
      IF (ns.eq.1) GOTO 285
      DO 280 i=2,ns
	      vol1(i)=vol1(i-1)+vol(i)
280   CONTINUE
285   CONTINUE
!      if((flow-exvl).gt.1.0) print*, flow-exvl
      NNNCO=2
      LAYCO=1
      CALL RESINT(NNNCO,LAYCO)

      RETURN
      END SUBROUTINE OUTFLO
!***********************************************************
      REAL*8 FUNCTION DV1(z1,z2,da,avdel,hh,delt,delb)
!******************************************************************
!  FUNCTION TO CALCLULATE THE PROPORTION OF FLUID WITHDRAWN
!  FROM ANY LAYER, GIVEN THE depth OF ITS TOP AND BOTTOM,
!  USING A CURVE WHICH FITS THE REGION OF FLUID DRAWN IN A 
!  GIVEN TIME TO DECIDE WHICH SET OF WITHDRAWAL CURVES TO 
!  USE. IF LARGE WITHDRAWAL USE FIRST SET, OTHERWISE THE 2ND.
!------------------------------------------------------------------------
!		USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE 

      REAL(8), PARAMETER :: zero=0.0d0,one=1.0d0,two=2.0d0
	   REAL(8), PARAMETER :: four=4.0d0,eight=8.0d0
	   REAL(8), PARAMETER :: pt25=0.25d0,pt5=0.5d0,pt75=0.75d0,pt9=0.9d0
	   REAL(8), PARAMETER :: pi=3.1415926535897932384626433832795d0
	   REAL(8) ::  A
      REAL(8) ::  AVDEL
      REAL(8) ::  DA
      REAL(8) ::  DA4
      REAL(8) ::  DA7
      REAL(8) ::  DELB
      REAL(8) ::  DELT
      REAL(8) ::  HH
!IF      REAL(8) ::  SQR
      REAL(8) ::  S1
      REAL(8) ::  S2
      REAL(8) ::  S3
      REAL(8) ::  TEMP1
      REAL(8) ::  TEMP2
      REAL(8) ::  TMP1
      REAL(8) ::  TMP2
      REAL(8) ::  TMP3
      REAL(8) ::  TMP4
      REAL(8) ::  TMP5
      REAL(8) ::  TMP6
      REAL(8) ::  TMP7
      REAL(8) ::  Z1
      REAL(8) ::  Z2
      REAL(8) ::  Z3

      IF (DA .lt. pt9*AVDEL) GOTO 100

!  CURVES FOR LARGE WITHDRAWAL.

      S1=(Z1-HH)/AVDEL
      S2=(Z2-HH)/AVDEL
      A=DA/AVDEL
      DV1=zero

!  IF TOP AND BOTTOM OF LAYER FALL WITHIN LOWER CURVE.

      IF (Z1 .le. pt25*A*AVDEL+HH) THEN
	      TMP1=S1+two/A/(one+pt5*A*(DELB/AVDEL+S1))
	      TMP2=S2+two/A/(one+pt5*A*(DELB/AVDEL+S2))
	      DV1=TMP1-TMP2
	      RETURN
      ENDIF

!  IF LAYER FALLS WITHIN UPPER CURVE.

      IF (Z2 .ge. HH+pt25*AVDEL*A) THEN
	      TMP1 = EXP(-1*SQRT(two)*(DELT/AVDEL+pt75*A))
	      TMP2 = one+pt5*A*DELB/AVDEL+SQR(A)/eight
	      TMP3 = (1-one/SQR(TMP2))*SQR((1+TMP1)/(1-TMP1))
	      TMP4 = SQRT(two)*(S1-DELT/AVDEL-A)
	      TMP5 = SQRT(two)*(S2-DELT/AVDEL-A)
	      TMP6 = four/(1+EXP(TMP4))+TMP4
	      TMP7 = four/(1+EXP(TMP5))+TMP5
	      DV1 = (TMP6-TMP7)/SQRT(two)*TMP3
	      RETURN
      ENDIF

!  IF JOIN OF CURVES LIES WITHIN THE LAYER.

      S3 = pt25*A
      TMP2 = S2+two/A/(one+pt5*A*(DELB/AVDEL+S2))
      TMP1 = S3+two/A/(one+pt5*A*(DELB/AVDEL+S3))
      DV1 = TMP1-TMP2
      TMP1 = EXP(-1*SQRT(two)*(DELT/AVDEL+pt75*A))
      TMP2 = one+pt5*A*DELB/AVDEL+SQR(A)/eight
      TMP3 = (1-one/SQR(TMP2))*SQR((1+TMP1)/(1-TMP1))
      TMP4 = SQRT(two)*(S1-DELT/AVDEL-A)
      TMP5 = SQRT(two)*(S3-DELT/AVDEL-A)
      TMP6 = four/(1+EXP(TMP4))+TMP4
      TMP7 = four/(1+EXP(TMP5))+TMP5
      DV1 = DV1+(TMP6-TMP7)/SQRT(two)*TMP3
      RETURN

!  CURVES FOR SMALL WITHDRAWAL

100   CONTINUE
      DV1=zero

!  IF TOP AND BOTTOM FALL WITHIN LOWER CURVE.

      IF (Z1 .le. DA/four+HH) THEN
	      TMP1 = Z1+(DELB+DA/four)/PI*SIN(PI*(Z1-HH-DA/four)/(DELB+DA/four))
	      TMP2 = Z2+(DELB+DA/four)/PI*SIN(PI*(Z2-HH-DA/four)/(DELB+DA/four))
	      DV1 = TMP1-TMP2
	      RETURN
      ENDIF

!  IF LAYER FALLS WITHIN UPPER CURVE.

      IF (Z2 .ge. HH+DA/four)THEN
	      TEMP1 = Z1+(DELT+DA*pt75)/PI*SIN(PI*(Z1-HH-DA/four)/(DELT+DA*pt75))
	      TEMP2 = Z2+(DELT+DA*pt75)/PI*SIN(PI*(Z2-HH-DA/four)/(DELT+DA*pt75))
	      DV1 = TEMP1-TEMP2
	      RETURN
      ENDIF

!  IF JOIN OF CURVES FALLS WITHIN THE LAYER.

      Z3 = pt25*DA+HH
      DA4 = DA/four
      DA7 = DA*pt75
      TMP1 = Z3+(DELB+DA4)/PI*SIN(PI*(Z3-HH-DA4)/(DELB+DA4))
      TMP2 = Z2+(DELB+DA4)/PI*SIN(PI*(Z2-HH-DA4)/(DELB+DA4))
      DV1 = TMP1-TMP2
      TMP3 = Z3+(DELT+DA7)/PI*SIN(PI*(Z3-HH-DA4)/(DELT+DA7))
      TMP4 = Z1+(DELT+DA7)/PI*SIN(PI*(Z1-HH-DA4)/(DELT+DA7))
      DV1 = TMP4-TMP3+DV1

      RETURN
      END FUNCTION DV1
! ***************************************************************      
	SUBROUTINE RESINT(icode,lnu)
!********************************************************************
!  FROM THE GIVEN PHYSICAL DATA EVALUATES ARRAYS OF DEPTHS AND AREAS
!  CORRESPONDING TO AN ARRAY OF VOLUMES (icode=2) OR ARRAYS OF VOLUME
!  AND AREAS FROM DEPTHS (icode=1) starting at layer lnu
!  contains bug fix for large inflows
!         SGS March, 1992
!------------------------------------------------------------------------
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  
      
	   REAL(8), PARAMETER :: zero=0.0d0,one=1.0d0,two=2.0d0,ten=10.0d0
	   REAL(8), PARAMETER :: pt1=0.1d0

      REAL(8) ::  volsum
	   REAL(8) ::  x,y

!IF	   REAL(8) ::  combin, combinv, densty

      INTEGER*4 i
      INTEGER*4 icode
      INTEGER*4 ij
      INTEGER*4 j
      INTEGER*4 k
      INTEGER*4 l
      INTEGER*4 ln
      INTEGER*4 lnu
      INTEGER*4 zz

      CHARACTER*80 MESS
10    CONTINUE

      IF(icode.eq.2)GOTO 200

!     Find volumes and areas given depths; ij points to storage table entry
!     closest to but less than the given depth, y is the relative location
!     of depth between storage table entries ij and ij+1.
!     Begin by calculating the total volume in the layer structure.

      x = depth(ns)*ten
      y = MOD(x,one)              ! Problems with the Definition

      ij = IDINT(x-y) 

      IF (ij .gt. numsto) THEN
	      y = y+float(ij-numsto)
	      ij = numsto
      ENDIF

      vol1(ns) = VV(ij)+y*DVV(ij)
      area(ns) = A(ij)+y*DADZ(ij)
      volsum=0.0d0

      DO 50 k=1,NUMINF
	      volsum = volsum+TOTIN(k)
50    CONTINUE

      vol1(ns) = vol1(ns)-volsum
      DO 100 i=lnu,ns-1
	      x = depth(i)*ten
	      y = MOD(x,one)                 ! Problems with the Definition

         ij = IDINT(x-y)

	      IF (ij .gt. numsto) THEN
	         y = y+float(ij-numsto)
	         ij = numsto
          ENDIF	
	      vol1(i) = VVDASH(ij)+y*DVVDA(ij)
	      area(i) = A(ij)+y*DADZ(ij)
100   CONTINUE
!
!SGS  These lines are part of the patch for the inflow bug Include WQ and Particles
!
	   IF (ns .gt. 1) THEN
		   IF (vol1(ns) .le. vol1(ns-1)) THEN
			   vol1(ns-1) = vol1(ns)+volsum
			   area(ns-1) = area(ns)
			   depth(ns-1)= depth(ns)

	         sal(ns-1) = COMBIN(sal(ns-1),vol(ns-1),den(ns-1),sal(ns),vol(ns),den(ns))
	         temp(ns-1) = COMBIN(temp(ns-1),vol(ns-1),den(ns-1),temp(ns),vol(ns),den(ns))

	         den(ns-1) = densty(temp(ns-1),sal(ns-1))
			   ns=ns-1
            DO 141 zz=1,28
	            wqual(ns-1,zz) = COMBINV(wqual(ns-1,zz),vol(ns-1),wqual(ns,zz),vol(ns))
141         CONTINUE
            DO 142 zz=1,7
	            cf(ns-1,zz) = COMBINV (cf(ns-1,zz),vol(ns-1),cf(ns,zz),vol(ns))
142         CONTINUE 
			   IF (ns .lt. lnu) THEN
				   lnu = lnu -1
				   IF (lnu .lt. 1) THEN
				      MESS = 'lnu LESS THAN ONE CONTACT EDL at UC-DAVIS.'
				      CALL ERRORS(MESS)
				   ENDIF
			   ENDIF
			   GOTO 10
		   ENDIF
	   ENDIF
      GOTO 500

!  calculate depths given volumes; j points to storage table volume
!  closest to but not greater than the given layer volume

200   volsum = vol1(ns)
      DO 202 k=1,NUMINF
	      volsum = volsum+TOTIN(k)
202   CONTINUE
      j = 1
205   CONTINUE
      IF (j .gt. numsto) THEN
	      j = numsto
      ELSE
	      IF (volsum .gt. VV(j)) THEN
	         j = j+1
	         GOTO 205
	      ENDIF
         j = j-1
      ENDIF
      depth(ns) = (float(j)+(volsum - VV(j))/DVV(j))*pt1
      IF (lnu .gt. 1) GOTO 250
      L=1
      j=1
!  find lowest layer (L) whose volume exceeds the first table entry

210   IF (vol1(L).gt.VVDASH(1))GOTO 300
      depth(L)=(vol1(L)/VVDASH(1))/ten
      L=L+1
      GOTO 210
250   CONTINUE
      x = depth(lnu-1)*ten
      y = MOD(x,one)							!PROBLEMS WITH DEFINITION DAVIS
      j = (x-y)
      L = lnu

300   CONTINUE
      DO 350 k=L,ns-1
310      CONTINUE
	      IF (j .gt. numsto)THEN
	         j=numsto
	      ELSE
	         IF (vol1(k) .gt. VVDASH(j))THEN
	            j=j+1
	            GOTO 310
	         ENDIF
	         j=j-1
	      ENDIF
	      depth(k)=(float(j)+(vol1(k)-VVDASH(j))/DVVDA(j))/ten
350   CONTINUE
!     determine areas


      DO 450 i=lnu,ns
	      x  = depth(i)*ten
	      y  = MOD(x,one)					 															    !PROBLEMS WITH DEFINITION DAVIS
	      ij = (x-y)

	      IF (ij.gt.numsto)ij=numsto
	      IF (ij .ne. 0) area(i) = A(ij)+ y * DADZ(ij)
	      IF (ij .eq. 0) area(i) = A(1) * depth(i) * 10.0d0
450   CONTINUE

!  calculate layer volumes and mean depths

500   CONTINUE

      ln=lnu
      IF (lnu .eq. 1)THEN
	      vol(1) = vol1(1)
	      DEPTHM(1)=depth(1)/two
	      ln=lnu+1
      ENDIF
      DO 550 i=ln,ns
	      vol(i) = vol1(i)-vol1(i-1)
	      IF((depth(i)-depth(i-1)).lt. 0.)THEN
	      WRITE(6,*)' Chaos in layer elevations ',i,i-1,depth(i),depth(i-1),'ns',ns,depth(ns)
!	      pause
	      ENDIF

	      IF(vol(i) .lt. 0) THEN
	         WRITE(6,*)'Warning:in routine RESINT, Vol layer < 0'
	         WRITE(6,555) icode, lnu, i, ns, vol(i),vol(ns)
555         format(I4,I4,I4,I5,2F10.2)
	         pause
	      ENDIF

	      DEPTHM(i)=(depth(i)+depth(i-1))/two
550   CONTINUE
      RETURN
      END SUBROUTINE RESINT
!***************************************************************************
	SUBROUTINE THICK
!***********************************************************************
!
! This subroutine checks the reservoir layer structure for compliance
! with the specified volume and depth limits.  Adjustments are made
! as required.
!
!      1)   layer structure is checked for minimum limits first,
!           combining the  checked layer with the smallest adjacent
!           layer IF necessary.  This may result in the formation
!           of a layer that exceeds the maximum limits.
!
!      2)   layer structure is checked for maximum limits, splitting
!           the checked layer IF necessary.  After splitting, the 
!           resulting layers will all be greater than their minimum 
!           limits provided VMAX >= 2 VMIN.  The default value is 
!           VMAX = 2 VMIN.
!
!      3)   the do loop counter now begins at kb (the number of the 
!           layer	immediately above the most recently adjusted layer)
!           instead of at 1.  This eliminates redundant checking of 
!           layers.
!
! VARIABLE LIST
!     DMAX                maximum allowable layer thickness
!     DMIN                minimum allowable layer thickness
!     kb                  FIRST LAYER ABOVE SPLIT OR MIXED LAYERS
!     kt                  NEW TOP LAYER NUMBER
!     M                   NUMBER OF LAYERS AFTER SPLITTING
!     V                   SPLIT VOLUME
!     VDOWN               VOLUME OF LAYER BELOW AMALGAMATION
!     VMAX                MAX ALLOWED VOLUME
!     VMIN                MIN ALLOWED VOLUME
!     VUP                 VOLUME OF LAYER ABOVE AMALGAMATION
!
!-------------------------------------------------------------------------------------
!
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  

      REAL(8), PARAMETER :: zero=0.0d0,one=1.0d0,two=2.0d0
      REAL(8) ::  D
      REAL(8) ::  D1
      REAL(8) ::  DELDP
!IF      REAL(8) ::  COMBIN,COMBINV
!IF      REAL(8) ::  densty
      REAL(8) ::  V
      REAL(8) ::  VDOWN
      REAL(8) ::  VUP
      INTEGER*4 i
      INTEGER*4 IADD
      INTEGER*4 icode
      INTEGER*4 j
      INTEGER*4 zz
      INTEGER*4 JOLD
      INTEGER*4 k
      INTEGER*4 kb
      INTEGER*4 KLAST
      INTEGER*4 kt
      INTEGER*4 M
      INTEGER*4 VSUMCHK

!  CHECK AGAINST VMIN

      KLAST=1
100   CONTINUE
      DO 200 i=KLAST,ns
	      IF (i.eq.1) THEN
	         DELDP=depth(i)
	      ELSE
	         DELDP=depth(i)-depth(i-1)
	      ENDIF   
	      IF(DELDP.lt.DMIN) GOTO 300 
!	  IF((vol(i).lt.VMIN).AND.(DELDP.lt.DMIN)) GOTO 300 
200   CONTINUE								

      GOTO 600

!  LAYER i IS AMALGAMATED WITH ITS SMALLEST NEIGHBOUR

300   CONTINUE
      IF (i.eq.1) THEN
	      VUP=zero
	      VDOWN=one
	      GOTO 350
      ENDIF
      IF (i .eq. ns) THEN
	      VUP=one
	      VDOWN=zero
	      GOTO 350
      ENDIF
      VUP = vol(i+1)
      VDOWN = vol(i-1)
350   CONTINUE

      j = i
      IF (VUP.gt.VDOWN) j = i-1
      sal(j)  = COMBIN (sal(j),vol(j),den(j),sal(j+1),vol(j+1),den(j+1))
      temp(j) = COMBIN (temp(j),vol(j),den(j),temp(j+1),vol(j+1),den(j+1))
    
      DO 360 zz=1,28
	      wqual(j,zz) = COMBINV(wqual(j,zz),vol(j),wqual(j+1,zz),vol(j+1))
360   CONTINUE

      DO 361 zz=1,7
	      cf(j,zz) = COMBINV(cf(j,zz),vol(j),cf(j+1,zz),vol(j+1))
361   CONTINUE

      den(j) = densty(temp(j),sal(j))
      vol(j) = vol(j)+vol(j+1)
      depth(j) = depth(j+1)
      vol1(j) = vol1(j+1)
      area(j) = area(j+1)
      ep(j)=ep(j+1)
      KLAST=j

!  RENUMBER LAYERS j+2,j+3,---,ns

      IF(j.eq.(ns-1)) GOTO 500
      kt=ns-1
      kb=j+1
      DO 400 k=kb,kt
	      depth(k) = depth(k+1)
	      den(k) = den(k+1)
	      temp(k) = temp(k+1)
	      sal(k) = sal(k+1)
	      DO 375 zz=1,28
	         wqual(k,zz) = wqual(k+1,zz)
375      CONTINUE
	      DO 376 zz=1,7
	         cf(k,zz) = cf(k+1,zz)
376      CONTINUE
	      vol(k) = vol(k+1)
	      vol1(k) = vol1(k+1)
	      area(k) = area(k+1)
	      ep(k) = ep(k+1)
400   CONTINUE
500   CONTINUE
      ns=ns-1

      GOTO 100

!  here when all layers have been checked for VMIN, DMIN

600   CONTINUE
      IF (ns.eq.1) GOTO 750
      DO 700 i=2,ns
	      DEPTHM(i)=(depth(i)+depth(i-1))/two
700   CONTINUE
750   CONTINUE
      DEPTHM(1)=depth(1)/two

!  check layers for VMAX
!sgs Flag to prevent top layer splitting more than once
!
      VSUMCHK = 0
      KLAST=1						   

1100  CONTINUE

      IF (VSUMCHK .eq. 1) RETURN

      DO 1200 i=KLAST,ns
	      IF (i.eq.1) THEN
	         DELDP=depth(i)
	      ELSE
	         DELDP=depth(i)-depth(i-1)
	      ENDIF
!
!SGS 
!
	      IF (i .eq. ns) VSUMCHK = 1	
	      IF(DELDP.gt.DMAX) GOTO 1300	
!	IF(vol(i).gt.VMAX.OR.DELDP.gt.DMAX) GOTO 1300	
1200  CONTINUE										

!  return to calling program when all layers have been checked
      RETURN
!  LAYER i IS SPLIT INTO M LAYERS
1300  CONTINUE
      M = 2
1350  CONTINUE
      V = vol(i)/M
      D = DELDP/M
      IF (V.LE.VMAX.AND.D.LE.DMAX) GOTO 1400
      M = M+1

!  IF M+ns IS GREATER THAN THE MAX NO. OF LAYERS, A MISTAKE WILL OCCUR
!   - AN ARRAY BOUND ERROR

      IF (M+ns .gt. MAXNS) THEN
	      WRITE(*,*)' ARRAY BOUND ERROR IN THICK - TOO MANY LAYERS '
	      WRITE(*,*)' ns = ',ns, ' M = ',M
	      STOP
      ENDIF
      GOTO 1350

!  renumber layers above split. IADD is the number of added layers
!  j is the new layer number, JOLD is the old (pre-split) layer number
! Include water quality and particles

1400  CONTINUE
      IADD = M-1
      KLAST = i
      IF (i.eq.ns) GOTO 1500
      kt = ns+IADD
      kb = i+M
      DO 1450 j=kt,kb,-1
	      JOLD = j-IADD
	      vol1(j) = vol1(JOLD)
	      den(j) = den(JOLD)
	      sal(j) = sal(JOLD)
	      temp(j) = temp(JOLD)
         DO 1430 zz=1,28
	         wqual(j,zz) = wqual(JOLD,zz)
1430     CONTINUE
         DO 1431 zz=1,7
	         cf(j,zz)=cf(JOLD,zz)
1431     CONTINUE
         vol(j) = vol(JOLD)
         ep(j) = ep(JOLD)
1450  CONTINUE

!  process the added layers
!  include water quality

1500  CONTINUE

!      d1=zero
!      IF (i .ne. 1) d1 = depth(i-1)
!      hgt = depth(i) - d1

      DO 1550 k=i,i+IADD
	      vol(k) = V
	      IF (k .eq. 1)THEN
	         vol1(k)=vol(k)
	      ELSE
	         vol1(k)=vol1(k-1)+vol(k)
	      ENDIF
	      den(k) = den(i)
	      sal(k) = sal(i)
	      temp(k) = temp(i)
         DO 1520 zz=1,28
	         wqual(k,zz)=wqual(i,zz)
1520     CONTINUE
         DO 1521 zz=1,7
	         cf(k,zz)=cf(i,zz)
1521     CONTINUE
         ep(k) = ep(i)
1550  CONTINUE
      ns = ns+IADD

!  get new depths for layers i thru ns

	   icode = 2
      CALL RESINT(icode,i)	  

2     CONTINUE
      GOTO 1100      
      END SUBROUTINE THICK
!******************************************************
      SUBROUTINE STABLE
!******************************************************
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  

!IF      REAL(8) ::  COMBIN,COMBINV
!IF      REAL(8) ::  densty
      REAL(8), PARAMETER	:: two = 2.0d0
      INTEGER*4 i,k,ii,zz


100   IF(ns.eq.1)RETURN

!     find an instability

      DO 200 k=ns,2,-1
	      IF(den(k).gt.den(k-1))THEN
	         GOTO 300
	      ENDIF
200   CONTINUE
      RETURN

300   CONTINUE

!     mix the unstable layer with the layer below
!     include water quality and particles

      temp(k-1)=COMBIN(temp(k),vol(k),den(k),temp(k-1),vol(k-1),den(k-1))
      sal(k-1)=COMBIN(sal(k),vol(k),den(k),sal(k-1),vol(k-1),den(k-1))
      DO 340 zz=1,28
	      wqual(k-1,zz)=COMBINV(wqual(k,zz),vol(k),wqual(k-1,zz),vol(k-1))
340   CONTINUE
      DO 350 ii=1,7
	      cf(k-1,ii)=COMBINV(cf(k,ii),vol(k),cf(k-1,ii),vol(k-1))
350   CONTINUE
      depth(k-1)=depth(k)
      vol(k-1)=vol(k-1)+vol(k)
      vol1(k-1)=vol1(k)
      area(k-1)=area(k)
      IF((k-1).eq.1) DEPTHM(k-1)=depth(k-1)/two !BUG FIX BY QUIM PEREZ 02/99
      IF((k-1).NE.1) DEPTHM(k-1)=(depth(k-1)+depth(k-2))/two
      den(k-1)=densty(temp(k-1),sal(k-1))
      ep(k-1)=ep(k)

!  adjust layer numbering of layers above mixing
!  include water quality and particles

      ns=ns-1
      DO 400 i=k,ns
	      ep(i)=ep(i+1)
	      depth(i)=depth(i+1)
	      DEPTHM(i)=DEPTHM(i+1)
	      sal(i)=sal(i+1)
	      temp(i)=temp(i+1)
	      DO 360 zz=1,28
	         wqual(i,zz)=wqual(i+1,zz)
360      CONTINUE
	      DO 370 ii=1,7
	         cf(i,ii)=cf(i+1,ii)
370      CONTINUE
	      vol(i)=vol(i+1)
	      vol1(i)=vol1(i+1)
	      area(i)=area(i+1)
	      den(i)=den(i+1)
400   CONTINUE

      GOTO 100
      
      END SUBROUTINE STABLE

!****************************************************************************
    !  SUBROUTINE POSTPRO(JDAY,NCHL,SWITCH,DATELEN,JDATE,BITNAME)
      SUBROUTINE POSTPRO(SWITCH,thestamp)
!*******************************************************************************
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  
!       processing for the main water quality program - gets binary output

      INTEGER*4 i,j,NL,IT1,IS1,IS2,IT2,JJ,zz

      REAL*8 WQUAL1(28),WQUAL2(28)
      REAL*8 CF1(7),CF2(7)
!
!     CHANGED TO INCREASE THE NUMBER OF LAYERS
!
      REAL*8 TTEM(maxns),SSAL(maxns),WWQUAL(maxns,28),CCF(maxns,7)
      REAL*8 DDEP(maxns)
!IF	   INTEGER*4 nchl
      CHARACTER*20 FSIM
     	CHARACTER*10 thestamp

      LOGICAL*1 FLAG
      LOGICAL*4 SWITCH

!     Convert internal nutrients back to a concentration relative to chla
!     Desactivated, because new definition of IN and IP nutrients for
!     Lake Tahoe
!      
!      DO 130 i=1,ns
!	      DO 130 j=1,nchl
!	         IF(wqual(i,j).eq.0.0)wqual(i,j)=0.001
!	         wqual(i,j+10)=wqual(i,j+10)/wqual(i,j)
!	         wqual(i,j+16)=wqual(i,j+16)/wqual(i,j)
!130   CONTINUE

      FLAG = .TRUE.
      NL=1
      TTEM(1)=temp(1)
      DDEP(1)=depth(1)
      SSAL(1)=sal(1)

      DO 150 JJ=1,28
	      WWQUAL(1,JJ)=wqual(1,JJ)
150   CONTINUE
!   bug
	   DO 151 JJ=1,7
	      CCF(1,JJ)=cf(1,JJ)
151   CONTINUE       

      IT1=INT(temp(1)*20000.)
      IS1=INT(sal(1)*20000.)

      DO 220 zz=1,28
	      WQUAL1(zz)=(wqual(1,zz))
220   CONTINUE

      DO 221 zz=1,7
	      CF1(zz)=(cf(1,zz))
221   CONTINUE

      DO 1000 i=2,ns
         IT2=INT(temp(i)*20000.)
         IS2=INT(sal(i)*20000.)

         DO zz=1,28
	         WQUAL2(zz)=(wqual(i,zz))
         ENDDO

         DO zz=1,7
	         CF2(zz)=(cf(i,zz))
         ENDDO

         IF(IT2.NE.IT1.OR.IS2.NE.IS1                                    &
            .OR.WQUAL2(1).NE.WQUAL1(1).OR.WQUAL2(2).NE.WQUAL1(2)        &
            .OR.WQUAL2(3).NE.WQUAL1(3).OR.WQUAL2(4).NE.WQUAL1(4)        &
            .OR.WQUAL2(5).NE.WQUAL1(5).OR.WQUAL2(6).NE.WQUAL1(6)        &
            .OR.WQUAL2(7).NE.WQUAL1(7).OR.WQUAL2(8).NE.WQUAL1(8)        &
            .OR.WQUAL2(9).NE.WQUAL1(9).OR.WQUAL2(10).NE.WQUAL1(10)      &
            .OR.WQUAL2(11).NE.WQUAL1(11).OR.WQUAL2(12).NE.WQUAL1(12)    &
            .OR.WQUAL2(13).NE.WQUAL1(13).OR.WQUAL2(14).NE.WQUAL1(14)    &
            .OR.WQUAL2(15).NE.WQUAL1(15).OR.WQUAL2(16).NE.WQUAL1(16)    &
            .OR.WQUAL2(17).NE.WQUAL1(17).OR.WQUAL2(18).NE.WQUAL1(18)    &
            .OR.WQUAL2(19).NE.WQUAL1(19).OR.WQUAL2(20).NE.WQUAL1(20)    &
            .OR.WQUAL2(21).NE.WQUAL1(21).OR.WQUAL2(22).NE.WQUAL1(22)    &
            .OR.WQUAL2(23).NE.WQUAL1(23).OR.WQUAL2(24).NE.WQUAL1(24)    &
            .OR.WQUAL2(25).NE.WQUAL1(25).OR.WQUAL2(26).NE.WQUAL1(26)    &
            .OR.WQUAL2(27).NE.WQUAL1(27).OR.WQUAL2(28).NE.WQUAL1(28)    &
            .OR.CF2(1).NE.CF1(1).OR.CF2(2).NE.CF1(2).OR.CF2(3).NE.CF1(3)   &
            .OR.CF2(4).NE.CF1(4).OR.CF2(5).NE.CF1(5)                    &
            .OR.CF2(6).NE.CF1(6).OR.CF2(7).NE.CF1(7))THEN

            IF (FLAG) GOTO 500
!  Save the layer bounding the top of the mixed region

	         NL=NL+1
	         TTEM(NL)=temp(i-1)
	         DDEP(NL)=depth(i-1)
	         SSAL(NL)=sal(i-1)
	         DO zz=1,28
	            WWQUAL(NL,zz)=wqual(i-1,zz)
	         ENDDO
	         DO zz=1,7
	            CCF(NL,zz)=cf(i-1,zz)
	         ENDDO

	         FLAG=.TRUE.

!  save the layer that is different (i.e. IT2,IS2)

500         NL=NL+1
	         TTEM(NL)=temp(i)
	         DDEP(NL)=depth(i)
	         SSAL(NL)=sal(i)
	         DO zz=1,28
	            WWQUAL(NL,zz)=wqual(i,zz)
	         ENDDO
	         DO zz=1,7
	            CCF(NL,zz)=cf(i,zz)
	         ENDDO 

	         IT1=IT2
	         IS1=IS2

	         DO 560 zz=1,28
	            WQUAL1(zz)=WQUAL2(zz)
560         CONTINUE
	         DO 561 zz=1,7
	            CF1(zz)=CF2(zz)
561         CONTINUE
         ELSE
	         FLAG= .FALSE.
         ENDIF
1000  CONTINUE

      IF (FLAG) GOTO 1070

!  save the surface layer IF more than 1 layer in upper mixed layer

      NL=NL+1
      TTEM(NL)=temp(ns)
      DDEP(NL)=depth(ns)
      SSAL(NL)=sal(ns)
      DO 1020 zz=1,28
	      WWQUAL(NL,zz)=wqual(ns,zz)
1020  CONTINUE
      DO 1021 zz=1,7
	      CCF(NL,zz)=cf(ns,zz)
1021  CONTINUE
1070  CONTINUE

!2015/9	   IF(NMETSW.NE.1)THEN
!2015/9	      DO 2345 i=1,ns
!2015/9	         DO 2345 j=24,28
!2015/9		         WWQUAL(i,j)=0.0
!2015/92345    CONTINUE
!2015/9	   ENDIF
!2015/9      IF(ZOOSW.NE.1)THEN
!2015/9	      DO 1022 i=1,ns
!2015/9	         wqual(i,22)=0.0
!2015/9	         wqual(i,23)=0.0
!2015/91022     CONTINUE
!2015/9      ENDIF
!2015/9      IF(PARTSW.NE.1)THEN
!2015/9	      DO 1023 i=1,ns
!2015/9	         DO 1023 j=1,7
!2015/9	            cf(i,j)=0.0
!2015/91023     CONTINUE
!2015/9      ENDIF
!
!     Create a unique name for the simulation by the date and time method
!     New routine which avoid Y2K is used instead old procedure
!
      IF (SWITCH) THEN
	      fsim = 'W' // trim(adjustl(thestamp)) // '.sim'
         open(10,file=FSIM,FORM='UNFORMATTED',STATUS='UNKNOWN')
	   ENDIF

	   WRITE(*,*)'SIMULATION WRITTEN TO : ',FSIM
	   WRITE(*,*) JDAY
!            		print*,'CHECK',JDAY,WWQUAL(NL,28)
	   WRITE(10)JDAY,NL,NCHL,ZOOSW,PARTSW,NMETSW,(DDEP(j),TTEM(j),SSAL(j),     &
            WWQUAL(j,1),WWQUAL(j,2),WWQUAL(j,3),WWQUAL(j,4),WWQUAL(j,5),      &
            WWQUAL(j,6),WWQUAL(j,7),WWQUAL(j,8),WWQUAL(j,9),WWQUAL(j,10),     &
            WWQUAL(j,11),WWQUAL(j,12),WWQUAL(j,13),WWQUAL(j,14),              &
            WWQUAL(j,15),WWQUAL(j,16),WWQUAL(j,17),WWQUAL(j,18),              &
            WWQUAL(j,19),WWQUAL(j,20),WWQUAL(j,21),WWQUAL(j,22),              &
            WWQUAL(j,23),WWQUAL(j,24),WWQUAL(j,25),WWQUAL(j,26),              &
            WWQUAL(j,27),WWQUAL(j,28),CCF(j,1),CCF(j,2),CCF(j,3),             &
            CCF(j,4),CCF(j,5),CCF(j,6),CCF(j,7),j=1,NL)

!  Multiply internal nutrients by chla
!     Desactivated according to the new definition for Tahoe
!
!      DO 1200 i=1,ns
!	DO 1200 j=1,nchl
!	    IF(wqual(i,j).eq.0.001)wqual(i,j)=0.0
!	    wqual(i,j+10)=wqual(i,j+10)*wqual(i,j)
!	    wqual(i,j+16)=wqual(i,j+16)*wqual(i,j)
!1200   CONTINUE
 

5000  RETURN
      END SUBROUTINE POSTPRO
!*******************************************************************
	   SUBROUTINE INTERPOLATE_LAYER_DATA (nlayers,idepth,density,tempy,salty)
!*******************************************************************
!
!     INTERPOLATE LAYER_DATA  -- lakenumber
!
!     
!     THIS SUBROUTINE INTERPOLATES THE depth
!     
!     and density data onto a series of points of increment every
!     
!     0.1 m. The interpolation assumes that all 3 parameters are
!
!     of constant value within each layer.(ZERO'th ORDER INTERPOL.)
!
!*******************************************************************
!     Enhanced to prevent temp inversions
!      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  

	   INTEGER(4)  :: nlayers 
	   INTEGER(4)  :: para, pare
	   INTEGER(4)  :: i,ii,j
	   REAL(8)		:: densty
	   REAL(8)		:: this_depth
	   REAL(8)     :: IDEPTH(MAXNS)
	   REAL(8)     :: DENSITY(MAXNS),TEMPY(MAXNS),SALTY(MAXNS)
	   i = 1
	   this_depth = REAL(i)*0.1	
!
!     Ensure Monotonic Temp Profile
!
!	pare = 0
!	DO while(pare .eq. 0)
!        para = 0
!		DO ii = 1, ns-1
!40		CONTINUE
!			IF (temp(ii+1) .le. temp(ii)) THEN
!				temp(ii+1) = temp(ii)   + 0.0000001D0	! Remember that the Profiles may be very sensitive!!!!!
!				para = 1                            ! Remember that Sal is not corrected at this point!!!!!
!			ENDIF
!		ENDDO
!		IF(para.eq.0) THEN
!			pare = 1
!		ENDIF
!	ENDDO
!
!     Recalculates density
!
	   DO j = 1,ns
	      DO 4410 while (this_depth .le. depth(j))
	         idepth(i)  = this_depth
	         density(i) = den(j)   !densty(temp(j),sal(j))
	         tempy(i)   = temp(j)
	         salty(i)   = sal(j)
	         i = i+1
	         this_depth = REAL(i)*0.1
4410     CONTINUE      
	   ENDDO
	
	   nlayers = i-1
	
	   RETURN
      END SUBROUTINE INTERPOLATE_LAYER_DATA
!********************************************************************
	   REAL*8 FUNCTION CALC_XMOM(nlayers)
!*******************************************************************
!
!     CALCULATE XMOM  -- lakenumber
!
!     NOTE:   DVV & DADZ DONOT START AT ELEMENT ZERO (hence '+1' in reference)
!             (DVV IS THE SAME AS DVDZ)
!                     CHECKED AND O.K.
!*******************************************************************
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  
	
	   INTEGER(4):: i,nlayers
	   REAL(8)   :: bfsq(0:maxnlayers),mid_depth(0:maxnlayers)
	   REAL(8)   :: func_XMOM,func_XMOMO
!IF	   REAL(8)   :: gprime
	   
		func_XMOM  = 0.0d0
		func_XMOMO = 0.0d0

		mid_depth(0)=idepth(0)/2.0d0

	   DO i = 1,(nlayers-1)
		   mid_depth(i)=(idepth(i)+idepth(i-1))/2.0

		   bfsq(i) = gprime(density(i),density(i-1))/(mid_depth(i)-mid_depth(i-1))	
		   func_XMOM = func_XMOM + idepth(i-1)*bfsq(i)*(dvv((i)+1)+dvv((i-1)+1))/2.0d0
		   func_XMOMO = func_XMOMO + bfsq(i)*(DVV((i)+1)+dvv((i-1)+1))/2.0d0
	   ENDDO

	   IF (func_XMOMO .le. 0.0d0) THEN
		   func_XMOM = -1.0
	   ELSE 
		   func_XMOM = func_XMOM/func_XMOMO
	   ENDIF
		
	   CALC_XMOM = func_XMOM

	   RETURN
      END FUNCTION CALC_XMOM

!******************************************************************************
	SUBROUTINE LNPE3(nlayers,xpp,zcp,xmasspe3,xpe)
!****************************************************************************** 
!                                                                     ######
!     lnpe3  -- lakenumber                                            IS THIS
!                                                                     SAME AS
!     Calculates the Potential Energy stored.                         PE3 
!     All totals are relative to the bottom.                          (1dmain)
!                                                                     ######
!             CHECKED AND O.K.
!*******************************************************************    ??????
      USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  
		
	   INTEGER(4)::    nlayers 
	   REAL(8)   ::    zcp,xmasspe3,xpp
	   REAL(8)   ::    zcvp,ab,da,ht,dens,xvol,xvolp
	   REAL(8)   ::    xmassp,xpe,zcv,xmass
	   INTEGER(4)::    i_pe3,ij,il,l,ibdep

!     GET THE LAYER CONTAINING IBDEP
	   ibdep=0

!     GET THE CENTRE OF VOLUME AND POTENTIAL ENERGY, TOTAL MASS ETC.
!     LOCATE THE SURFACE IN RELATION TO THE STORAGE TABLE.  IJ IS
!     THE LAST STORAGE ENTRY BEFORE THE SURFACE, Y IS THE DISTANCE 
!     TO THE SURFACE FROM THERE.

	   ij = int(idepth(nlayers-1)*10.0d0)
!     INITIALIZE THE LOOP. IL IS THE CURRENT LAYER NUMBER, DA THE AREA 
!     GRADIENT STARTING AT IBDEP (HERE, IBDEP ALWAYS EQUALS ZERO)

	   xvol  = 0.0d0          ! Center of Volume
	   xmass = 0.0d0          ! Center of Mass
	   xpe   = 0.0d0          ! Potential Energy
	   zcv   = 0.0d0          !

	   da = dble(A(1))
	   ab = 0.0d0
	   il = ibdep
	   ht = 0.0d0
	
	   DO i_pe3 = 0,(ij-1)
		   il = il + 1	
		   dens = dble(density(il)+1000.0d0)
		   ht   = ht + 0.1
		   IF (i_pe3 .ne. 0) da = dble(dadz((i_pe3-1)+1)*10.0)
		   zcvp   = ab*((0.1*ht)-0.005)+(da/6)*((0.03*ht)-0.001)  !!!!! DUBS ABOUT THIS EQUATION SHAPE FORM!!!!!
		   zcv    = zcv + zcvp
		   xpe    = xpe + dens*zcvp
		   xvolp  = 0.1*(ab+da*0.05)
		   xvol   = xvol + xvolp
		   xmassp = xvolp*dens
		   xmass  = xmass + xmassp
		   ab     = ab + 0.1*da
	   ENDDO

!     THIS TAKES US TO THE LAST STORAGE TABLE ENTRY.  NOW GO TO THE SURFACE.

	   IF(il .eq. (nlayers-1)) THEN
		   ht     = (idepth(nlayers-1))
		   zcvp   = ab*((0.1*ht)-0.005)+(da/6)*((0.03*ht)-0.001)
		   zcv    = zcv + zcvp
		   xpe    = xpe + dens*zcvp
		   xvolp  = 0.1*(ab+da*0.05)
		   xvol   = xvol + xvolp
		   xmassp = xvolp*dens
		   xmass  = xmass + xmassp
	   ELSE
	     DO 4405 l = il,(nlayers-1)
		      dens   = density(l)
		      ht     = idepth(l)
		      zcvp   = ab*((0.1*ht)-0.005)+(da/6)*((0.03*ht)-0.001)
		      zcv    = zcv + zcvp
		      xpe    = xpe + dens*zcvp
		      xvolp  = 0.1*(ab+da*0.05)
		      xvol   = xvol + xvolp
		      xmassp = xvolp*dens
		      xmass  = xmass + xmassp
		      ab     = dble(A(l+1))
4405     CONTINUE
	   ENDIF

!     NOTE THAT UNITS OF XP ARE MILLIONS OF KG M**2/SEC**2,
!     UNITS OF XMASS ARE MILLIONS OF KG

	   zcp      = REAL(zcv/xvol)
	   xpp      = REAL((xpe/xmass-zcv/xvol)*9.810d0)
	   xmasspe3 = REAL(xmass)

	   RETURN
      END SUBROUTINE LNPE3
!********************************************************C
      SUBROUTINE moment(ndata,n,ave,adev,sdev,var,skew,curt)
!********************************************************!
!
!     Routine copy from Numerical Recipies
!     Modified for double precision
!---------------------------------------------------------------
 
 !     USE DLMWQ_VARIABLES
	   IMPLICIT NONE   
      SAVE  
      
      INTEGER*4 n
      REAL*8 adev,ave,curt,sdev,skew,var,ndata(n)
      INTEGER(4)     :: j
      REAL*8 p,s,ep
      IF(n.le.1)pause 'n must be at least 2 in moment'
      s=0.
      DO 11 j=1,n
        s=s+ndata(j)
11    CONTINUE
      ave=s/n
      adev=0.
      var=0.
      skew=0.
      curt=0.
      ep=0.
      DO 12 j=1,n
         s=ndata(j)-ave
         ep=ep+s
         adev=adev+abs(s)
         p=s*s
         var=var+p
         p=p*s
         skew=skew+p
         p=p*s
         curt=curt+p
12    CONTINUE
      adev=adev/n
      var=(var-ep**2/n)/(n-1)
      sdev=sqrt(var)
      IF(var.ne.0.)THEN
         skew=skew/(n*sdev**3)
         curt=curt/(n*var**2)-3.
      ELSE
         WRITE(*,*) 'no skew or kurtosis when zero variance in moment'
      ENDIF
      RETURN
      END SUBROUTINE MOMENT
!**************************************************************************
      END MODULE DLMWQ_SERVICES
