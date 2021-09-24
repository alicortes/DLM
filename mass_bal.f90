!*********************************************************************	
      SUBROUTINE mass_bal(var,outflow,conc)
!*********************************************************************	
	! jday - julian day of the mass ballance calculation
	! lakemass - mass of variable in the lake
	! inmass - mass of variable in the inflow stack
	! outmass - mass of variable in the outflow
!*********************************************************************	
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE

		INTEGER*4:: i,var,k,ntot
		REAL*8:: lakemass, inmass, outmass, massbal, outflow, conc

	! calculate the mass in the lake
		lakemass=0.0d0
		DO i=1,ns
			lakemass=lakemass+wqual(i,var)*vol(i)		! integrate lake mass
		ENDDO

	!Calculate the mass in the inflow stack
		DO i=1,numinf
!		inmass = inmass+floinf(i)*wqual(i,var)		! add in the inflow masses
			DO k=1,icnt(i)
				inmass = qdown(i,k)*wqdown(i,var,k)		! add in inlfow masses including stacks
			ENDDO
		ENDDO

!Calculate the mass in the outflow
		outmass = outflow*conc							! add up outflow masses

		massbal = lakemass-outmass+inmass				! total the mass balance for the step


!		WRITE(83,94)jday,nosecs,lakemass,outmass,massbal
94		FORMAT(2I10,3F16.4)

		END SUBROUTINE mass_bal
	
!*********************************************************************	
      SUBROUTINE h2o_bal(outflow)
!*********************************************************************	
	! jday - julian day of the mass ballance calculation
	! lakemass - mass of variable in the lake
	! inmass - mass of variable in the inflow stack
	! outmass - mass of variable in the outflow
!*********************************************************************	
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE

		INTEGER*4:: i,k,ntot
		REAL*8:: lakemass, inmass, outmass, massbal, outflow

!		calculate the mass in the lake
		lakemass=0.0d0
		DO i=1,ns
			lakemass=lakemass+vol(i)		! integrate lake mass
		ENDDO

!		Calculate the mass in the inflow stack
		DO i=1,numinf
!		inmass = inmass+floinf(i)*wqual(i,var)		! add in the inflow masses
			inmass=0.0d0
			DO k=1,icnt(i)
				inmass = inmass+qdown(i,k)		! add in inlfow masses including stacks
			ENDDO
		ENDDO

!		Calculate the mass in the outflow

		outmass=outflow
		massbal = lakemass-outmass+inmass				! total the mass balance for the step


!		WRITE(83,94)jday,lakemass,inmass,outmass,massbal
94		FORMAT(I10,4F16.4)
		END SUBROUTINE h2o_bal
!*********************************************************************	
!*********************************************************************	
		SUBROUTINE injection(var,conc,top,bottom)
!*********************************************************************	
	! var - variable to inject
	! conc - concentration of variable to inject
	! top - depth from surface of top of injection layer
	! bottom - depth from surface of bottom of injection layer
!*********************************************************************	
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE

		INTEGER*4:: i,var,ntop,nbottom
		REAL*8:: conc,top,bottom

!find layers where top and bottom of injection occur
		nbottom=0							! start with finding the bottom layer
		IF (bottom == depth(ns)) THEN
			nbottom = 1						! nbottom is the bottom layer
		ELSE
			DO WHILE ((depth(ns)-depth(ns-nbottom-1))<bottom) 
				nbottom=nbottom+1
			ENDDO
			ns=ns+1							! add the additional layer
			DO i=0,nbottom
				depth(ns-i)=depth(ns-1-i)	! assign depths to layers above split
			ENDDO
			depth(ns-nbottom-1)=depth(ns)-bottom		! assign depth to the split
			nbottom=ns-nbottom				! assign actual layer number of the bottom
			CALL RESINT(1,nbottom-1)		! reset volumes, T, S, and wqual
		ENDIF
!then work on the top layer
		ntop=0								! now find the top layer
		IF (top == 0) THEN
			ntop = ns						! ntop is the top layer
		ELSE
			nbottom=nbottom-1				! bottom will be moved 1 layer down		
			DO WHILE ((depth(ns)-depth(ns-ntop-1))<top) 
				ntop=ntop+1
			ENDDO
			ns=ns+1							! add the additional layer
			DO i=0,ntop
				depth(ns-i)=depth(ns-1-i)	! assign depths to layers above split
			ENDDO
			depth(ns-ntop-1)=depth(ns)-top			! assign depth to the split
			ntop=ns-ntop-1					! assign actual layer number of the bottom
			CALL RESINT(1,ntop)				! reset volumes, T, S, and wqual
		ENDIF
!add tracer to the desired layers
		DO i=nbottom,ntop
			wqual(i,var)=conc
		ENDDO

		END SUBROUTINE injection

!*********************************************************************
