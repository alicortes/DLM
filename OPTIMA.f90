!**********************************************************************************************
		SUBROUTINE OPTIMA(NLAY,DDEP,VARI,VARM,ERROR,ERROR_E,ERROR_H,ERROR_C,SUMA,SUMA_E,SUMA_H)      
!**********************************************************************************************
!*****************************************************************
!*	Written by Joaquim Perez (08/01/99)
!*	Generates the FITTING function to be optimized by POWELL
!*    or Genetics Algorithm. Function from Jorgensen 1990, and Patterson, 1984
!*	Givin a layer distribution,the position of the termocline is calculated
!*	THEN, the Root Mean Squared or the normalized RMS error of the variable 
!*	can be estimated by either all of the layers, or separated
!*	into epi and hypo contributions.
!*
!*	Variable List
!*    VARM: vector with measured data
!*    VARI: vector with simulated data
!*	W_E, W_H weight functions 
!*	MEAN_FIELD volumetric mean of the measured variable.
!*    MEAN_SIMUL volumetric mean of the simulated variable
!*    ERROR = measure of the error within all layers 
!*    ERROR_VOL = volume-weighted error function
!*	ERROR_C = error made with Epilimnetic and Hypolimnetic contributions with weight paramters
!*	ERROR_E = error made with Epilimnetic contribution
!*	ERROR_H = error made with Hypolimnetic contribution
!*    
!*****************************************************************
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE

		INTEGER*4 NLAY, NLAY_E,NLAY_H
		INTEGER*4 NLAYERS
		INTEGER*4 i,j,k
		REAL*8 MEAN_FIELD,MEAN_SIMUL
		REAL*8 VOL_TOTAL,VOL_TOTAL_H,VOL_TOTAL_E,VAR_MEDIA
     	REAL*8 ERROR,ERROR_VOL,ERROR_E,ERROR_H,ERROR_C
		REAL*8 VARI(maxns),VARM(maxns),DDEP(maxns) 
		REAL*8 suma,suma_h,suma_e,suma_vol
		REAL*8 var_media_h,var_media_t
		REAL*8 mean_field_e, mean_field_h
		REAL*8 mean_simul_h, mean_simul_t
		REAL*8 lake_top,lake_bottom
		REAL*8 thermdepth	
		REAL*8 xp,zc,xmass
		REAL*8 W_E,W_H
! 
!     Position of the Termocline
! 		
		DO i =1,ndays
			IF(vector_day(i).eq.jday) THEN
				thermdepth = vector_thermdepth(i)
				GOTO 4406
			ENDIF
		ENDDO
4406  CONTINUE
! 
!     Initialize variables
! 
		NLAY_H = 0
		NLAY_E = 0
		W_E = 1.0
		W_H = 1.0
		SUMA = 0.0
		SUMA_H = 0.0
		SUMA_E = 0.0 
		SUMA_VOL = 0.0
		ERROR = 0.0
		ERROR_E = 0.0
		ERROR_H = 0.0
		ERROR_C = 0.0
		ERROR_VOL = 0.0
		MEAN_FIELD = 0.0
      MEAN_SIMUL = 0.0
		VOL_TOTAL = 0.0
		VOL_TOTAL_H = 0.0
		VOL_TOTAL_E = 0.0
      VAR_MEDIA = 0.0
		VAR_MEDIA_H = 0.0
		VAR_MEDIA_T = 0.0
		MEAN_FIELD_E = 0.0
		MEAN_FIELD_H = 0.0
		MEAN_SIMUL_H = 0.0
		MEAN_SIMUL_T = 0.0  
! 	     
!     RMS definition
! 
      IF (CALIBRATION.eq.'RMS'.OR.CALIBRATION.EQ.'rms') THEN 
			IF(thermdepth.eq.-99.0) THEN    !Comment for Variables
				DO i = 1,NLAY
					MEAN_FIELD = MEAN_FIELD + VARM(I)*VOL(I)
					MEAN_SIMUL = MEAN_SIMUL + VARI(I)*VOL(I)
					VOL_TOTAL = VOL_TOTAL+VOL(I) 
					SUMA = SUMA + ((VARI(I)-VARM(I)))**2
				ENDDO
			ELSE
				DO i = 1,NLAY
					SUMA = SUMA + ((VARI(I)-VARM(I)))**2
					IF (DDEP(I).LE. XMOM) THEN
						NLAY_H = NLAY_H + 1
						MEAN_FIELD_H = MEAN_FIELD_H + VARM(I)*VOL(I)
						MEAN_SIMUL_H = MEAN_SIMUL_H + VARI(I)*VOL(I)
						VOL_TOTAL_H = VOL_TOTAL_H+VOL(I)
						SUMA_H = SUMA_H + ((VARI(I)-VARM(I)))**2
					ELSE
						NLAY_E = NLAY_E + 1
						MEAN_FIELD_E = MEAN_FIELD_E + VARM(I)*VOL(I)
						MEAN_SIMUL_T = MEAN_SIMUL_T + VARI(I)*VOL(I)
						VOL_TOTAL_E = VOL_TOTAL_E+VOL(I)
						SUMA_E = SUMA_E + ((VARI(I)-VARM(I)))**2
					ENDIF
				ENDDO
			ENDIF

			IF(thermdepth.eq.-99.0) THEN 
				MEAN_FIELD = MEAN_FIELD/VOL_TOTAL
 				MEAN_SIMUL = MEAN_SIMUL/VOL_TOTAL
				ERROR     = DSQRT(SUMA/NLAY)
				ERROR_E   = 0.0
				ERROR_H   = 0.0
				ERROR_C   = 0.0
				ERROR_VOL = DSQRT(SUMA_VOL/VOL_TOTAL)	
			ELSE
				MEAN_FIELD =   (MEAN_FIELD_E/VOL_TOTAL_E)+ (MEAN_FIELD_H/VOL_TOTAL_H)
 				MEAN_SIMUL =   (MEAN_SIMUL_T/VOL_TOTAL_E)+ (MEAN_SIMUL_H/VOL_TOTAL_H)
				ERROR   = DSQRT(SUMA/NLAY)
				ERROR_E =  W_E*DSQRT(SUMA_E/NLAY_E) 
		 
				IF(NLAY_H.NE.0) THEN
					ERROR_H = W_H*DSQRT(SUMA_H/NLAY_H)
					ERROR_C   = W_E*DSQRT(SUMA_E/NLAY_E) + W_H*DSQRT(SUMA_H/NLAY_H)	
				ELSE
					WRITE(*,*)' THERMOCLINE HAS REACHED THE BOTTOM'
					ERROR_C = W_E*DSQRT(SUMA_E/NLAY_E) 
				ENDIF
			ENDIF
      ENDIF
! 
!     PSI definition (RMS normalized by the mean)
! 

		IF(CALIBRATION.eq.'PSI'.OR.CALIBRATION.EQ.'psi') THEN
			IF(thermdepth.eq.0) THEN    !Comment for Variables
				DO i = 1,NLAY
					MEAN_FIELD = MEAN_FIELD + VARM(I)	 
				ENDDO
				MEAN_FIELD = MEAN_FIELD/NLAY
				DO i = 1,NLAY
					SUMA = SUMA + ((VARI(I)-VARM(I))/MEAN_FIELD)**2
				ENDDO	
				ERROR = DSQRT(SUMA/NLAY)
				ERROR_E = 0.0
				ERROR_H = 0.0
				ERROR_C = 0.0 
			ELSE
				DO i = 1,NLAY
					MEAN_FIELD = MEAN_FIELD + VARM(I)		
					IF (DDEP(I).LE. XMOM) THEN
						NLAY_H = NLAY_H + 1
						MEAN_FIELD_H = MEAN_FIELD_H + VARM(I)	 
					ELSE
						NLAY_E = NLAY_E + 1
						MEAN_FIELD_E = MEAN_FIELD_E + VARM(I)
					ENDIF
				ENDDO
				IF(NLAY_H.NE.0) THEN
					MEAN_FIELD_H = MEAN_FIELD_H / NLAY_H
				ELSE
					MEAN_FIELD_H = 0.0
				ENDIF
				MEAN_FIELD_E = MEAN_FIELD_E / NLAY_E

				DO i = 1,NLAY
					SUMA = SUMA + ((VARI(I)-VARM(I))/MEAN_FIELD)**2
					IF (DDEP(I).LE. XMOM) THEN
						SUMA_H = SUMA_H + ((VARI(I)-VARM(I))/MEAN_FIELD_H)**2
					ELSE		 
						SUMA_E = SUMA_E + ((VARI(I)-VARM(I))/MEAN_FIELD_E)**2
					ENDIF
				ENDDO

				ERROR = DSQRT(SUMA/NLAY)
				ERROR_E =  W_E*DSQRT(SUMA_E/NLAY_E)

				IF(NLAY_H.NE.0) THEN
					ERROR_H = W_H*DSQRT(SUMA_H/NLAY_H)
					ERROR_C   = W_E*DSQRT(SUMA_E/NLAY_E) + W_H*DSQRT(SUMA_H/NLAY_H)	
				ELSE
					WRITE(*,*)' THERMOCLINE HAS REACHED THE BOTTOM'
					ERROR_C = W_E*DSQRT(SUMA_E/NLAY_E) 
				ENDIF
			ENDIF
		ENDIF

		RETURN				   
		END SUBROUTINE OPTIMA