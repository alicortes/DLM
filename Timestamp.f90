!*********************************************************
		SUBROUTINE timestamp(thestamp)
!************************************************************
!Return an unique date/time stamp format
!-----------------------------------------------------------
      USE DLMWQ_VARIABLES      
      USE DLMWQ_SERVICES
		IMPLICIT NONE
      SAVE
		CHARACTER*10 thestamp
		CHARACTER*10 dada_temp(3)

      CALL DATE_AND_TIME(dada_temp(1),dada_temp(2),dada_temp(3))

		thestamp(1:8)= dada_temp(1)
		thestamp(9:10) = dada_temp(2)(1:6)
		RETURN
		END SUBROUTINE timestamp
