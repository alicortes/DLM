!**************************************************************
		SUBROUTINE SHORE_EROSION(jyear,redn,Part_shore)
!*****************************************************************
!gbs  Shoreline erosion particles rates are estimated by Goloka Sahoo in 2006 from 
!gbs  total sediment erosion values reported by Kenneth D. Adams (2002)
!
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE
	
		INTEGER*4 i,jyear
		
		REAL*8 Part_shore(7)
		REAL*8 prof_toplayer
    	REAL(8):: rate_Particles(7)
		REAL*8 redn
! 
!gbs shoreline erosion paricles rates for Lake Tahoe (# particles/day/m2)  
! 
		DATA rate_Particles /4.652724E+06,1.356254E+06,						&
					2.382425E+05,3.574813E+04,3.024782E+04, 6.669400E+03,0.00E+00/
! 
!gbs Depth (m) of surface layer
! 
      prof_toplayer = depth(ns)-depth(ns-1)

!	print*,jday,jyear,redn,'shoreline'
		WRITE (*, 10) jday,jyear,redn
10		FORMAT(2i8,f10.2, 24x,'Shoreline erosion')

!	WRITE(1000,*) jyear,jday,redn,'SER'
! 
!gbs All shoreline eroded particles are inserted into the surface layer
! 
		DO i = 1,7
			cf(ns,i) = cf(ns,i) + redn*rate_Particles(i)/(prof_toplayer)		   ! redn                                      
			Part_shore(i) = Part_shore(i) +(redn*rate_Particles(i)/(prof_toplayer))*vol(ns)*1000.0d0	!redn
		ENDDO
		RETURN
		END SUBROUTINE SHORE_EROSION







